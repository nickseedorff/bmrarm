#' PX-DA MCMC routine to sample from HBMRVAR model
#'
#' @param formula an object of class "formula"; a symbolic description of the model to be fitted
#' @param data a dataframe containing outcome variables, covariates, and a patient or subject identifier
#' @param ordinal_outcomes a character string containing the names of the ordinal outcomes
#' @param patient_var name of the patient or subject identifier
#' @param sig_prior prior variance on the regression coefficients
#' @param all_draws logical with a default of FALSE which discards burn-in
#' @param nsim positive integer, number of iterations with default of 1000
#' @param burn_in positive integer, number of iterations to remove with default of 100. Must be >= 100.
#' @param thin positive integer, specifiers the period of saving samples. Default of 20 due to the high autocorrelation of the cutpoints
#' @param seed positive integer, seed for random number generation
#' @param verbose logical, print iteration number to keep track of progress
#' @param max_iter_rej maximum number of rejection algorithm attempts for multivariate truncated normal
#' @return mcmc
#' @importFrom zoo na.approx
#' @importFrom fastDummies dummy_cols
#' @importFrom magic adiag
#' @importFrom lme4 lmer VarCorr
#' @export

baseline_bmr2 <- function(formula, data, ordinal_outcome = c("y_ord"),
                          time_var = "time", patient_var = "patient_idx",
                          random_slope = F, ar_cov = TRUE, nsim = 1000,
                          burn_in = 100, thin = 10, seed = 14, verbose = TRUE,
                          sig_prior = 1000000000, sd_vec = c(0.15, 0.30)) {

  ## Create storage
  set.seed(seed)
  bmrarm_start(env = environment())
  cont_out_var <- setdiff(out_vars, ordinal_outcome)

  ## Starting values for ordinal outcome
  y[1:N_obs, 1] <- 0.5 + res_cuts[z, 1]
  y[is.infinite(y[, 1]), 1] <- -0.5
  y[is.na(y[, 1]), 1] <- 0

  ## Starting values for continuous outcomes
  df <- data.frame(patient = samp_info$pat_idx_long, y = as.numeric(y),
                   outcome = rep(1:N_outcomes, each = N_obs)) %>%
    group_by(patient, outcome) %>%
    mutate(y_interp = na.approx(y, na.rm = FALSE),
           y_interp = ifelse(!is.na(y_interp), y_interp,
                             ifelse(row_number() == n(),
                                    lag(y_interp), lead(y_interp))))
  y <- matrix(df$y_interp, ncol = N_outcomes)
  y[is.na(y)] <- 0
  i <- 2
  samp_info$burn_in <- burn_in
  samp_info$max_iter <- 10000000

  ## Emprical bayes priors for the random effects covariance matrix
  if(random_slope) {
    form_string <- paste0("~ . + (", time_var, "|", patient_var, ")")
  } else {
    form_string <- paste0("~ . + (1|", patient_var, ")")
  }
  form_use <- update(formula, as.formula(form_string))

  ## Fit models
  ord_form <- reformulate(deparse(form_use[[3]]), response = ordinal_outcome)
  cont_form <- reformulate(deparse(form_use[[3]]), response = cont_out_var)
  ord_mod <- lmer(ord_form, data = data)
  cont_mod <- lmer(cont_form, data = data)

  ## Prior for
  priors <- c(as.data.frame(VarCorr(ord_mod))[1:N_pat_effects, "vcov"],
              as.data.frame(VarCorr(cont_mod))[1:N_pat_effects, "vcov"])
  prior_mat <- diag(priors) * N_pat_effects

  for(i in 2:nsim) {
    samp_info$num_iter <- i

    ## Regression coefficients
    vals <- bmrarm_fc_sig_beta(y, X, Z_kron, cur_draws, samp_info)
    #print("sb")
    res_beta[, i] <- cur_draws$beta <- vals$beta
    res_sig[, i] <- cur_draws$sigma <- vals$sig

    # Autoregressive parameter
    if(samp_info$ar_cov) {
      vals <- bmrarm_mh_ar(y, X, Z_kron, cur_draws, samp_info)
      #print("ar")
      res_ar[i] <- cur_draws$ar <- vals$ar
      res_accept[i, 2] <- vals$accept
    }

    ## Subject specific effects
    vals <- bmrarm_fc_patient2(y, z, X, cur_draws, samp_info, prior_mat)
    res_pat_sig[, i] <- cur_draws$pat_sig <- vals$pat_sig
    res_pat_eff[,, i] <- cur_draws$pat_effects <- vals$pat_effects
    #res_pat_sig_sd[, i] <- cur_draws$pat_sig_sd <- vals$pat_sig_sd
    #print(head(cur_draws$pat_effects))

    ## Latent values, missing values, cut points
    y_cuts <- bmrarm_fc_y_cuts(y, z, X, Z_kron, cur_draws, samp_info)
    #print("y")
    y <-  res_y[,, i] <- y_cuts$y
    res_cuts[, i] <- cur_draws$cuts <- y_cuts$cuts
    res_accept[i, 1] <- y_cuts$accept

    ## Update missing values
    y <- res_y[,, i]<- bmrarm_fc_missing(y, z, X, Z_kron, cur_draws, samp_info)
    #print("miss")

    ## Cut points
    #if(i %% 50 == 0) print(i)
    #if(i %% 100 == 50) plot(res_cuts[4, ], type = "l")
    if(i %% 100 == 50) plot(res_ar, type = "l")
    if(i %% 100 == 0) plot(res_pat_sig[1, ], type = "l")
  }

  sim_use <- seq(burn_in + 1, nsim, by = thin)
  draws <- list(
    res_cuts = res_cuts[, sim_use],
    res_beta = res_beta[, sim_use],
    res_pat_sig = res_pat_sig[, sim_use],
    res_ar = res_ar[sim_use],
    res_sigma = res_sig[, sim_use],
    res_y = res_y[,, sim_use],
    res_pat_eff = res_pat_eff[,, sim_use],
    res_accept = res_accept[sim_use, ],
    samp_info = samp_info,
    X = X,
    Z_kron = Z_kron, z = z)
  draws
}


#' PX-DA MCMC routine to sample from HBMRVAR model
#'
#' @param formula an object of class "formula"; a symbolic description of the model to be fitted
#' @param data a dataframe containing outcome variables, covariates, and a patient or subject identifier
#' @param ordinal_outcomes a character string containing the names of the ordinal outcomes
#' @param patient_var name of the patient or subject identifier
#' @param sig_prior prior variance on the regression coefficients
#' @param all_draws logical with a default of FALSE which discards burn-in
#' @param nsim positive integer, number of iterations with default of 1000
#' @param burn_in positive integer, number of iterations to remove with default of 100. Must be >= 100.
#' @param thin positive integer, specifiers the period of saving samples. Default of 20 due to the high autocorrelation of the cutpoints
#' @param seed positive integer, seed for random number generation
#' @param verbose logical, print iteration number to keep track of progress
#' @param max_iter_rej maximum number of rejection algorithm attempts for multivariate truncated normal
#' @return mcmc
#' @importFrom zoo na.approx
#' @importFrom fastDummies dummy_cols
#' @importFrom magic adiag
#' @export

baseline_bmr3 <- function(formula, data, ordinal_outcome = c("y_ord"),
                          time_var = "time", patient_var = "patient_idx",
                          random_slope = F, ar_cov = TRUE, nsim = 1000,
                          burn_in = 100, thin = 10, seed = 14, verbose = TRUE,
                          sig_prior = 1000000000, sd_vec = c(0.15, 0.30)) {

  ## Create storage
  set.seed(seed)
  bmrarm_start(env = environment())

  ## Starting values for ordinal outcome
  y[1:N_obs, 1] <- 0.5 + res_cuts[z, 1]
  y[is.infinite(y[, 1]), 1] <- -0.5
  y[is.na(y[, 1]), 1] <- 0

  ## Starting values for continuous outcomes
  df <- data.frame(patient = samp_info$pat_idx_long, y = as.numeric(y),
                   outcome = rep(1:N_outcomes, each = N_obs)) %>%
    group_by(patient, outcome) %>%
    mutate(y_interp = na.approx(y, na.rm = FALSE),
           y_interp = ifelse(!is.na(y_interp), y_interp,
                             ifelse(row_number() == n(),
                                    lag(y_interp), lead(y_interp))))
  y <- matrix(df$y_interp, ncol = N_outcomes)
  y[is.na(y)] <- 0
  i <- 2
  samp_info$burn_in <- burn_in
  samp_info$max_iter <- 10000000
  for(i in 2:nsim) {
    samp_info$num_iter <- i

    ## Regression coefficients
    vals <- bmrarm_fc_sig_beta(y, X, Z_kron, cur_draws, samp_info)
    #print("sb")
    res_beta[, i] <- cur_draws$beta <- vals$beta
    res_sig[, i] <- cur_draws$sigma <- vals$sig

    # Autoregressive parameter
    if(samp_info$ar_cov) {
      vals <- bmrarm_mh_ar(y, X, Z_kron, cur_draws, samp_info)
      #print("ar")
      res_ar[i] <- cur_draws$ar <- vals$ar
      res_accept[i, 2] <- vals$accept
    }

    ## Subject specific effects
    vals <- bmrarm_fc_patient4(y, z, X, cur_draws, samp_info)
    #res_pat_sig[, i] <- cur_draws$pat_sig <- sim_data$sig_alpha
    #res_pat_sig[, i] <- cur_draws$pat_sig <- sim_data$sig_alpha
    res_pat_eff[,, i] <- cur_draws$pat_effects <- vals$pat_effects
    res_pat_sig[, i] <- cur_draws$pat_sig <- vals$pat_sig
    res_pat_sig_sd[, i] <- cur_draws$pat_sig_sd <- vals$pat_sig_sd
    #print(head(cur_draws$pat_effects))

    ## Latent values, missing values, cut points
    y_cuts <- bmrarm_fc_y_cuts(y, z, X, Z_kron, cur_draws, samp_info)
    #print("y")
    y <-  res_y[,, i] <- y_cuts$y
    res_cuts[, i] <- cur_draws$cuts <- y_cuts$cuts
    res_accept[i, 1] <- y_cuts$accept

    ## Update missing values
    y <- res_y[,, i]<- bmrarm_fc_missing(y, z, X, Z_kron, cur_draws, samp_info)
    #print("miss")

    ## Cut points
    #if(i %% 50 == 0) print(i)
    #if(i %% 100 == 50) plot(res_cuts[4, ], type = "l")
    if(i %% 100 == 50) plot(res_ar, type = "l")
    if(i %% 100 == 0) plot(res_pat_sig[1, ], type = "l")
  }

  sim_use <- seq(burn_in + 1, nsim, by = thin)
  draws <- list(
    res_cuts = res_cuts[, sim_use],
    res_beta = res_beta[, sim_use],
    res_pat_sig = res_pat_sig[, sim_use],
    res_ar = res_ar[sim_use],
    res_sigma = res_sig[, sim_use],
    res_y = res_y[,, sim_use],
    res_pat_sig_sd = res_pat_sig_sd[, sim_use],
    res_pat_eff = res_pat_eff[,, sim_use],
    res_accept = res_accept[sim_use, ],
    samp_info = samp_info,
    X = X,
    Z_kron = Z_kron, z = z)
  draws
}

#' PX-DA MCMC routine to sample from HBMRVAR model
#'
#' @param formula an object of class "formula"; a symbolic description of the model to be fitted
#' @param data a dataframe containing outcome variables, covariates, and a patient or subject identifier
#' @param ordinal_outcomes a character string containing the names of the ordinal outcomes
#' @param patient_var name of the patient or subject identifier
#' @param sig_prior prior variance on the regression coefficients
#' @param all_draws logical with a default of FALSE which discards burn-in
#' @param nsim positive integer, number of iterations with default of 1000
#' @param burn_in positive integer, number of iterations to remove with default of 100. Must be >= 100.
#' @param thin positive integer, specifiers the period of saving samples. Default of 20 due to the high autocorrelation of the cutpoints
#' @param seed positive integer, seed for random number generation
#' @param verbose logical, print iteration number to keep track of progress
#' @param max_iter_rej maximum number of rejection algorithm attempts for multivariate truncated normal
#' @return mcmc
#' @importFrom zoo na.approx
#' @importFrom fastDummies dummy_cols
#' @importFrom magic adiag
#' @export

baseline_bmr_test2 <- function(formula, data, ordinal_outcome = c("y_ord"),
                               time_var = "time", patient_var = "patient_idx",
                               random_slope = F, ar_cov = TRUE, nsim = 1000,
                               burn_in = 100, thin = 10, seed = 14, verbose = TRUE,
                               sig_prior = 1000000000, sd_vec = c(0.15, 0.05),
                               which_fc = c(T, T, T, T, T), base0 = F) {

  ## Create storage
  set.seed(seed)
  bmrarm_start(env = environment())

  ## Starting values for ordinal outcome
  y[1:N_obs, 1] <- 0.5 + res_cuts[z, 1]
  y[is.infinite(y[, 1]), 1] <- -0.5
  y[is.na(y[, 1]), 1] <- 0

  ## Starting values for continuous outcomes
  df <- data.frame(patient = samp_info$pat_idx_long, y = as.numeric(y),
                   outcome = rep(1:N_outcomes, each = N_obs)) %>%
    group_by(patient, outcome) %>%
    mutate(y_interp = na.approx(y, na.rm = FALSE),
           y_interp = ifelse(!is.na(y_interp), y_interp,
                             ifelse(row_number() == n(),
                                    lag(y_interp), lead(y_interp))))
  y <- matrix(df$y_interp, ncol = N_outcomes)
  y[is.na(y)] <- 0
  i <- 2
  cur_draws$beta <- matrix(as.numeric(sim_data$beta), 2)
  cur_draws$sigma <- sim_data$sigma
  cur_draws$ar <- sim_data$ar
  if(base0) cur_draws$ar <- 0
  cur_draws$pat_effects <- sim_data$alpha
  cur_draws$pat_sig <- sim_data$sig_alpha
  cur_draws$cuts <- sim_data$cuts
  samp_info$burn_in <- burn_in
  samp_info$max_iter <- 10000000
  y <- cbind(sim_data$truth$y1, sim_data$truth$y2)
  for(i in 2:nsim) {
    samp_info$num_iter <- i

    ## Regression coefficients
    if(which_fc[1]) {
      vals <- bmrarm_fc_sig_beta(y, X, Z_kron, cur_draws, samp_info)
      #print("sb")
      res_beta[, i] <- cur_draws$beta <- vals$beta
      res_sig[, i] <- cur_draws$sigma <- vals$sig
      #res_sig_sd[, i] <- cur_draws$sig_sd <- vals$sig_sd
    } else {
      res_beta[, i] <- cur_draws$beta
      res_sig[, i] <- cur_draws$sigma
    }

    # Autoregressive parameter
    if(samp_info$ar_cov & which_fc[2]) {
      vals <- bmrarm_mh_ar(y, X, Z_kron, cur_draws, samp_info)
      #print("ar")
      res_ar[i] <- cur_draws$ar <- vals$ar
      res_accept[i, 2] <- vals$accept
      #cur_draws$ar <- 0
    } else {
      res_ar[i] <- cur_draws$ar
    }

    ## Subject specific effects
    if(which_fc[3]) {
      vals <- bmrarm_fc_patient4(y, z, X, cur_draws, samp_info)
      #print("pat")
      res_pat_sig[, i] <- cur_draws$pat_sig <- vals$pat_sig
      res_pat_eff[,, i] <- cur_draws$pat_effects <- vals$pat_effects
      res_pat_sig_sd[, i] <- cur_draws$pat_sig_sd <- vals$pat_sig_sd
      #res_pat_sig_q[, i] <- vals$pat_sig_q
      #res_accept[i, 3:6] <- vals$accept_vec
      #res_pat_sig[, i] <- cur_draws$pat_sig
    } else {
      res_pat_sig[, i] <- cur_draws$pat_sig
      res_pat_eff[,, i] <- cur_draws$pat_effects
    }

    ## Latent values, missing values, cut points
    if(which_fc[4]) {
      y_cuts <- bmrarm_fc_y_cuts(y, z, X, Z_kron, cur_draws, samp_info)
      #print("y")
      y <-  res_y[,, i] <- y_cuts$y
      res_cuts[, i] <- cur_draws$cuts <- y_cuts$cuts
      res_accept[i, 1] <- y_cuts$accept
    } else {
      res_y[,, i] <- y
      res_cuts[, i] <- cur_draws$cuts
    }

    ## Update missing values
    if(which_fc[5]) {
      y <- res_y[,, i]<- bmrarm_fc_missing(y, z, X, Z_kron, cur_draws, samp_info)
    }
    #print("miss")

    ## Cut points
    if(i %% 50 == 0) print(i)
    #if(i %% 100 == 50) plot(res_cuts[4, ], type = "l")
    if(i %% 100 == 50) plot(res_ar, type = "l")
    if(i %% 100 == 0) plot(res_pat_sig[1, ], type = "l")
  }

  sim_use <- seq(burn_in + 1, nsim, by = thin)
  draws <- list(
    res_cuts = res_cuts[, sim_use],
    res_beta = res_beta[, sim_use],
    res_pat_sig = res_pat_sig,
    res_ar = res_ar[sim_use],
    res_sigma = res_sig[, sim_use],
    res_pat_sig_q = res_pat_sig_q[, sim_use],
    res_pat_sig_sd = res_pat_sig_sd[, sim_use],
    res_y = res_y[,, sim_use],
    res_pat_eff = res_pat_eff[,, sim_use],
    samp_info = samp_info,
    res_accept = res_accept[sim_use, ],
    X = X,
    Z_kron = Z_kron, z = z)
  draws
}

#' Full conditional draws of the latent continuous values
#'
#' @param y matrix of multivariate observations
#' @param z matrix of ordinal voutcomes
#' @param sig covariance matrix for the VAR process
#' @param sig0 for the initial values
#' @param M transition matrix for the VAR(1) component
#' @param cuts current threshold values
#' @param miss_mat locations of missing values
#' @param samp_info information for which locations to sample
#' @param num_iter current iteration number
#' @import tmvtnorm
#' @return matrix
#' @export

bmrarm_fc_patient2 <- function(y, z, X, cur_draws, samp_info, prior_mat) {

  ## Generate full sigma matrix
  N_pat <- samp_info$N_pat
  sig_alpha_inv <- chol2inv(chol(cur_draws$pat_sig))
  sig_inv <- chol2inv(chol(cur_draws$sigma))
  resid_mat <- y - X %*% cur_draws$beta
  N_pat_eff <- ncol(samp_info$pat_z_kron[[1]])

  ## Sigma inverses only needed is autoregressive covariance matrix
  if(samp_info$ar_cov) {
    sig_list <- get_sig_list(cur_draws, samp_info)
  }

  ## Patient effects
  res <- matrix(NA, nrow = N_pat, ncol = N_pat_eff)

  for(i in 1:N_pat) {
    ## Get locations and time matrix for patient
    locs <- samp_info$pat_locs[[i]]
    time_ind <- samp_info$pat_time_ind[i]
    resid_vec <- as.numeric(resid_mat[locs, ])
    pat_Z <- samp_info$pat_z_kron[[i]]

    ## Patient specific covariance matrix
    if(samp_info$ar_cov) {
      pat_sig_inv <- sig_list$sig_inv_list[[time_ind]]
    } else {
      pat_sig_inv <- kronecker(sig_inv, diag(rep(1, samp_info$pat_N_obs[[i]])))
    }

    ## Cross products, covariance, alpha hat
    Z_sig_prod <- crossprod(pat_Z, pat_sig_inv)
    post_cov <- chol2inv(chol(Z_sig_prod %*% pat_Z + sig_alpha_inv))
    post_mean <- post_cov %*% Z_sig_prod %*% resid_vec
    L <- t(chol(post_cov))
    res[i, ] <- L %*% rnorm(length(post_mean)) + post_mean
  }

  ## New Covariance matrix for patient specific effects
  pat_sig <- rinvwishart(N_pat + N_pat_eff,
                         crossprod(res) + prior_mat)
  list(pat_effects = res, pat_sig = pat_sig)
}

#' Full conditional draws of the latent continuous values
#'
#' @param y matrix of multivariate observations
#' @param z matrix of ordinal voutcomes
#' @param sig covariance matrix for the VAR process
#' @param sig0 for the initial values
#' @param M transition matrix for the VAR(1) component
#' @param cuts current threshold values
#' @param miss_mat locations of missing values
#' @param samp_info information for which locations to sample
#' @param num_iter current iteration number
#' @import tmvtnorm
#' @return matrix
#' @export

bmrarm_fc_patient3 <- function(y, z, X, cur_draws, samp_info) {

  ## Generate full sigma matrix
  N_pat <- samp_info$N_pat
  sig_alpha_inv <- chol2inv(chol(cur_draws$pat_sig))
  sig_inv <- chol2inv(chol(cur_draws$sigma))
  resid_mat <- y - X %*% cur_draws$beta
  N_pat_eff <- ncol(samp_info$pat_z_kron[[1]])

  ## Sigma inverses only needed is autoregressive covariance matrix
  if(samp_info$ar_cov) {
    sig_list <- get_sig_list(cur_draws, samp_info)
  }

  ## Patient effects
  res <- matrix(NA, nrow = N_pat, ncol = N_pat_eff)

  for(i in 1:N_pat) {
    ## Get locations and time matrix for patient
    locs <- samp_info$pat_locs[[i]]
    time_ind <- samp_info$pat_time_ind[i]
    resid_vec <- as.numeric(resid_mat[locs, ])
    pat_Z <- samp_info$pat_z_kron[[i]]

    ## Patient specific covariance matrix
    if(samp_info$ar_cov) {
      pat_sig_inv <- sig_list$sig_inv_list[[time_ind]]
    } else {
      pat_sig_inv <- kronecker(sig_inv, diag(rep(1, samp_info$pat_N_obs[[i]])))
    }

    ## Cross products, covariance, alpha hat
    Z_sig_prod <- crossprod(pat_Z, pat_sig_inv)
    post_cov <- chol2inv(chol(Z_sig_prod %*% pat_Z + sig_alpha_inv))
    post_mean <- post_cov %*% Z_sig_prod %*% resid_vec
    L <- t(chol(post_cov))
    res[i, ] <- L %*% rnorm(length(post_mean)) + post_mean
  }

  ## Correlation matrix
  sd_inv <- diag(1 / cur_draws$pat_sig_sd)
  pat_sig_q <- rinvwishart(N_pat + N_pat_eff + 1,
                           sd_inv %*% crossprod(res) %*% sd_inv +
                             diag(rep(1, 4)))

  ## SD parameters
  accept_vec <- rep(0, N_pat_eff)
  for(i in 1:N_pat_eff) {
    ## Propose new values
    cur_draws2 <- cur_draws
    cur_draws2$pat_sig_sd[i] <- exp(rnorm(1, log(cur_draws2$pat_sig_sd[i]),
                                          sd = samp_info$sd_pat_sd))
    cur_draws2$pat_sig <- diag(cur_draws2$pat_sig_sd) %*% pat_sig_q %*% diag(cur_draws2$pat_sig_sd)

    ## Calculate comparison values
    pat_inv_old <- chol2inv(chol(cur_draws$pat_sig))
    pat_inv_new <- chol2inv(chol(cur_draws2$pat_sig))

    comp_old <- N_pat / 2 * determinant(pat_inv_old, logarithm = T)[[1]][1] -
      0.5 * sum(diag(res %*% pat_inv_old %*% t(res))) +
      dnorm(log(cur_draws$pat_sig_sd[i]), 0, sd = 10000, log = T)

    comp_new <- N_pat / 2 * determinant(pat_inv_new, logarithm = T)[[1]][1] -
      0.5 * sum(diag(res %*% pat_inv_new %*% t(res))) +
      dnorm(log(cur_draws2$pat_sig_sd[i]), 0, sd = 10000, log = T)
    compar_val <- comp_new - comp_old

    if(compar_val >= log(runif(1))) {
      cur_draws$pat_sig_sd[i] <- cur_draws2$pat_sig_sd[i]
      accept_vec[i] <- 1
    }
  }

  # ## SD parameters
  # accept_vec <- rep(0, N_pat_eff)
  # ## Propose new values
  # cur_draws2 <- cur_draws
  # cur_draws2$pat_sig_sd <- exp(rnorm(N_pat_eff, log(cur_draws2$pat_sig_sd),
  #                                       sd = samp_info$sd_pat_sd))
  # cur_draws2$pat_sig <- diag(cur_draws2$pat_sig_sd) %*% pat_sig_q %*% diag(cur_draws2$pat_sig_sd)
  #
  # ## Calculate comparison values
  # pat_inv_old <- chol2inv(chol(cur_draws$pat_sig))
  # pat_inv_new <- chol2inv(chol(cur_draws2$pat_sig))
  #
  # comp_old <- N_pat / 2 * determinant(pat_inv_old, logarithm = T)[[1]][1] -
  #   0.5 * sum(diag(res %*% pat_inv_old %*% t(res)))
  #
  # comp_new <- N_pat / 2 * determinant(pat_inv_new, logarithm = T)[[1]][1] -
  #   0.5 * sum(diag(res %*% pat_inv_new %*% t(res)))
  # compar_val <- comp_new - comp_old
  #
  # if(compar_val >= log(runif(1))) {
  #   cur_draws$pat_sig_sd <- cur_draws2$pat_sig_sd
  #   accept_vec <- 1
  # }

  cur_draws$pat_sig_sd <- c(cur_draws$pat_sig_sd[1] , 1, 1, 1)
  list(pat_effects = res,
       pat_sig = diag(cur_draws$pat_sig_sd) %*% pat_sig_q %*% diag(cur_draws$pat_sig_sd),
       pat_sig_sd = cur_draws$pat_sig_sd, accept_vec = accept_vec, pat_sig_q = pat_sig_q)
}

#' Full conditional draws of the latent continuous values
#'
#' @param y matrix of multivariate observations
#' @param z matrix of ordinal voutcomes
#' @param sig covariance matrix for the VAR process
#' @param sig0 for the initial values
#' @param M transition matrix for the VAR(1) component
#' @param cuts current threshold values
#' @param miss_mat locations of missing values
#' @param samp_info information for which locations to sample
#' @param num_iter current iteration number
#' @import tmvtnorm
#' @return matrix
#' @export

bmrarm_fc_patient4 <- function(y, z, X, cur_draws, samp_info) {

  ## Generate full sigma matrix
  N_pat <- samp_info$N_pat
  sig_alpha_inv <- chol2inv(chol(cur_draws$pat_sig))
  sig_inv <- chol2inv(chol(cur_draws$sigma))
  resid_mat <- y - X %*% cur_draws$beta
  N_pat_eff <- ncol(samp_info$pat_z_kron[[1]])

  ## Sigma inverses only needed is autoregressive covariance matrix
  if(samp_info$ar_cov) {
    sig_list <- get_sig_list(cur_draws, samp_info)
  }

  ## Patient effects
  res <- matrix(NA, nrow = N_pat, ncol = N_pat_eff)

  for(i in 1:N_pat) {
    ## Get locations and time matrix for patient
    locs <- samp_info$pat_locs[[i]]
    time_ind <- samp_info$pat_time_ind[i]
    resid_vec <- as.numeric(resid_mat[locs, ])
    pat_Z <- samp_info$pat_z_kron[[i]]

    ## Patient specific covariance matrix
    if(samp_info$ar_cov) {
      pat_sig_inv <- sig_list$sig_inv_list[[time_ind]]
    } else {
      pat_sig_inv <- kronecker(sig_inv, diag(rep(1, samp_info$pat_N_obs[[i]])))
    }

    ## Cross products, covariance, alpha hat
    Z_sig_prod <- crossprod(pat_Z, pat_sig_inv)
    post_cov <- chol2inv(chol(Z_sig_prod %*% pat_Z + sig_alpha_inv))
    post_mean <- post_cov %*% Z_sig_prod %*% resid_vec
    L <- t(chol(post_cov))
    res[i, ] <- L %*% rnorm(length(post_mean)) + post_mean
  }

  ## Update Correlation matrix
  sd_mat <- diag(cur_draws$pat_sig_sd)
  pat_sig <- rinvwishart(N_pat + N_pat_eff + 1, crossprod(res) + 4 * sd_mat)

  ## Update SDs
  pat_sig_inv <- chol2inv(chol(pat_sig))
  pat_sig_sd <- 1 / LaplacesDemon::rinvgamma(
    N_pat_eff, shape = (2 + 4) / 2,
    scale = 2 * diag(pat_sig_inv) + 1 / 100000000
  )

  list(pat_effects = res,
       pat_sig = pat_sig,
       pat_sig_sd = pat_sig_sd)
}

#' Full conditional draws of the regression coefficients
#'
#' @param subject_effects matrix of patient specific intercepts
#' @param sigma residual covariance matrix
#' @param prior_alpha prior term for shape
#' @param prior_alpha prior term for scale
#' @return scalar
#' @importFrom matrixcalc is.positive.definite
#' @importFrom LaplacesDemon rinvwishart dinvwishart
#' @importFrom MBSP matrix.normal
#' @export

bmrarm_fc_sig_beta2 <- function(y, X, Z_kron, cur_draws, samp_info) {

  ## Constant
  N_outcomes <- samp_info$N_outcomes
  N_covars <- samp_info$N_covars

  ## Mean of patient effects
  mean_vec <- rowSums(Z_kron * cur_draws$pat_effects[samp_info$pat_idx_long, ])
  resid_mat <- y - matrix(mean_vec, ncol = N_outcomes)

  ## Calculate posterior covariance matrix
  cov_vals <- matrix(0, N_covars, N_covars)
  mean_vals <- matrix(0, nrow = N_covars, ncol = N_outcomes)
  resid_vals <- matrix(0, N_outcomes, N_outcomes)

  ## Sigma inverses only needed is autoregressive covariance matrix
  if(samp_info$ar_cov) {
    sig_list <- get_sig_list(cur_draws, samp_info)
  }

  for(i in 1:samp_info$N_pat) {
    ## Get locations, partial residuals after subtracting person effects
    locs <- samp_info$pat_locs[[i]]
    X_pat <- samp_info$pat_X[[i]]
    time_ind <- samp_info$pat_time_ind[i]
    if(samp_info$ar_cov) {
      ## Summation of residuals, mean values, and covariance values
      time_inv <- sig_list$time_inv[[time_ind]]
      resid_vals <- resid_vals +
        crossprod(resid_mat[locs, ], time_inv) %*% resid_mat[locs, ]
      mean_vals <- mean_vals + crossprod(X_pat, time_inv) %*% resid_mat[locs, ]
      cov_vals <- cov_vals + crossprod(X_pat, time_inv) %*% X_pat
    } else {
      ## Summation of residuals, mean values, and covariance values
      resid_vals <- resid_vals + crossprod(resid_mat[locs, ], resid_mat[locs, ])
      mean_vals <- mean_vals + crossprod(X_pat, resid_mat[locs, ])
      cov_vals <- cov_vals + crossprod(X_pat, X_pat)
    }
  }

  ## Get posterior draw
  sd_inv <- diag(1 / cur_draws$sig_sd)
  prior <- diag(rep(0.00001, N_covars))
  x_inv <- chol2inv(chol(prior + cov_vals))
  beta_hat <- x_inv %*% mean_vals
  val <- resid_vals + 4 * sd_inv - t(mean_vals) %*% beta_hat

  ## Draw updates
  sig <- rinvwishart(samp_info$N_obs + N_outcomes + 1, val)
  beta <- matrix.normal(M = beta_hat, V = sig, U = x_inv)

  ## Update SDs
  sig_inv <- chol2inv(chol(sig))
  sig_sd <- LaplacesDemon::rinvgamma(
    N_outcomes, shape = (2 + N_outcomes) / 2,
    scale = 2 * diag(sig_inv) + 1 / 10000
  )

  list(beta = beta, sig = sig, sig_sd = sig_sd)
}

#' PX-DA MCMC routine to sample from HBMRVAR model
#'
#' @param formula an object of class "formula"; a symbolic description of the model to be fitted
#' @param data a dataframe containing outcome variables, covariates, and a patient or subject identifier
#' @param ordinal_outcomes a character string containing the names of the ordinal outcomes
#' @param patient_var name of the patient or subject identifier
#' @param sig_prior prior variance on the regression coefficients
#' @param all_draws logical with a default of FALSE which discards burn-in
#' @param nsim positive integer, number of iterations with default of 1000
#' @param burn_in positive integer, number of iterations to remove with default of 100. Must be >= 100.
#' @param thin positive integer, specifiers the period of saving samples. Default of 20 due to the high autocorrelation of the cutpoints
#' @param seed positive integer, seed for random number generation
#' @param verbose logical, print iteration number to keep track of progress
#' @param max_iter_rej maximum number of rejection algorithm attempts for multivariate truncated normal
#' @return mcmc
#' @importFrom zoo na.approx
#' @importFrom fastDummies dummy_cols
#' @importFrom magic adiag
#' @export

bmr_cv <- function(formula, data, ordinal_outcome = c("y_ord"),
                   time_var = "time", patient_var = "patient_idx",
                   random_slope = F, ar_cov = TRUE, nsim = 1000,
                   burn_in = 100, thin = 10, seed = 14, verbose = TRUE,
                   sig_prior = 1000000000, sd_vec = c(0.15, 0.30)) {

  ## Long version of dataset
  X <- model.matrix.lm(as.formula(formula), data = data, na.action = "na.pass")
  out_vars <- setdiff(all.vars(formula), colnames(X))
  cont_out_var <- setdiff(out_vars, ordinal_outcome)
  data_long <- rbind(
    select(data, patient_var, time_var, val = ordinal_outcome) %>%
      mutate(outcome = 1),
    select(data, patient_var, time_var, val = cont_out_var) %>%
      mutate(outcome = 2)) %>%
    mutate(folds = NA)

  ## Max number of observations
  max_obs <- data_long %>%
    group_by(pat_idx) %>%
    summarise(N = n()) %>%
    ungroup() %>%
    select(N) %>%
    unlist() %>%
    max()

  ## Folds
  obs_responses <- which(!is.na(data_long$val))
  set.seed(seed)
  folds <- loo::kfold_split_stratified(max(max_obs, 10),
                                       data_long$pat_idx[obs_responses])

  data_long$folds[obs_responses] <- folds
  data$folds1 <- data_long$folds[1:nrow(data)]
  data$folds2 <- data_long$folds[(1:nrow(data)) + nrow(data)]
  cv_val <- 0

  for(i in 1:max(folds)) {
    data_tmp <- data
    data_tmp[data$folds1 == i & !is.na(data$folds1), ordinal_outcome] <- NA
    data_tmp[data$folds2 == i & !is.na(data$folds2), cont_out_var] <- NA
    cv_locs <- list(all_locs = which(data_long$folds == i),
                    cont_locs = which(data$folds2 == i))

    samps <- baseline_bmr(formula = formula, data = data_tmp,
                          ordinal_outcome = ordinal_outcome,
                          patient_var = patient_var,
                          random_slope = random_slope,
                          time_var = time_var, ar_cov = ar_cov,
                          burn_in = burn_in, nsim = nsim, thin = thin,
                          seed = seed,
                          sd_vec = sd_vec)

    if(!ar_cov) samps$res_ar[] <- 0
    cv_val <- c(cv_val, get_loocv_ar(samps, cv_locs,
                                     z_true = data[, ordinal_outcome],
                                     y_true = data[, cont_out_var]))
    print(i)
    print(sum(cv_val))
  }
  list(cv_vals = cv_val[-1], sum_cv_vals = sum(cv_val))
}

#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @importFrom mvtnorm dmvnorm
#' @importFrom loo kfold_split_stratified
#' @import tidyr
#' @return scalar

strat_model_wrapper <- function(data, pat_idx, time, ord_out, cont_out) {
  ## Long version of dataset
  data_long <- rbind(
    select(data, pat_idx, time, val = ord_out) %>% mutate(outcome = 1),
    select(data, pat_idx, time, val = cont_out) %>% mutate(outcome = 2)) %>%
    mutate(folds = NA)

  ## Max number of observations
  max_obs <- data_long %>%
    group_by(pat_idx) %>%
    summarise(N = n()) %>%
    ungroup() %>%
    select(N) %>%
    unlist() %>%
    max()

  ## Folds
  obs_responses <- which(!is.na(data_long$val))

  ## Find unique folds
  correct_folds <- F
  fold_iter <- 1
  #while(fold_iter <= 100 & !correct_folds)
  folds <- loo::kfold_split_stratified(max(max_obs, 10),
                                       data_long$pat_idx[obs_responses])
  check <- data_long %>%
    group_by(pat_idx) %>%
    summarise(unique_fold = length(unique(na.omit(folds))),
              obs_responses = sum(!is.na(val)),
              correct_folds = sum(unique_fold != obs_responses)) %>%
    ungroup() %>%
    summarise(sum(correct_folds))


  data_long$folds[obs_responses] <- folds
  data$folds1 <- data_long$folds[1:nrow(data)]
  data$folds2 <- data_long$folds[(1:nrow(data)) + nrow(data)]
  cv_val <- 0

  for(i in 1:nrow(cv_vals)) {
    data_tmp <- data
    data_tmp$y_ord[data$folds1 == i] <- NA
    data_tmp$y2[data$folds2 == i] <- NA
    cv_locs <- list(all_locs = which(data_long$folds == i),
                    cont_locs = which(data$folds2 == i))

    samps <- baseline_bmr(formula = cbind(y_ord, y2) ~ time, data = data_tmp,
                          ordinal_outcome = "y_ord", patient_var = "pat_idx",
                          random_slope = T, time_var = "time", ar_cov = F,
                          burn_in = 100, nsim = 400, thin = 5, seed = 3, sd_vec = c(0.15, 0.05))
    samps$res_ar[] <- 0

    cv_val <- sum(get_loocv_ar(samps, cv_locs, z_true = data$y_ord,
                               y_true = data$y2)) + cv_val
    print(i)
    print(cv_val)
  }
  sum(cv_vals)
}
