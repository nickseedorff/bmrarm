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
#' @importFrom nlme lme VarCorr lmeControl corAR1
#' @export

baseline_bmr <- function(formula, data, ordinal_outcome = c("y_ord"),
                             time_var = "time", patient_var = "patient_idx",
                             random_slope = F, ar_cov = TRUE, nsim = 1000,
                             burn_in = 100, thin = 10, seed = 14, verbose = TRUE,
                             sig_prior = 1000000000, sd_vec = c(0.15, 0.30),
                             N_burn_trunc = 5, prior_siw_uni = c(0.1, 10)) {

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
    rand_form <- as.formula(paste0("~ ", time_var, "|", patient_var))
  } else {
    rand_form <- as.formula(paste0("~ 1|", patient_var))
  }

  ord_form <- reformulate(deparse(formula[[3]]), response = ordinal_outcome)
  cont_form <- reformulate(deparse(formula[[3]]), response = cont_out_var)

  ## Datasets and structures for models
  data_ord <- data[!is.na(data[, ordinal_outcome]), ]
  data_cont <- data[!is.na(data[, cont_out_var]), ]
  cor_struct <- if(ar_cov) {
    cor_struct <- corAR1(form = as.formula(~ 1 | pat_idx))
    #cor_struct <- NULL
  } else {
    cor_struct <- NULL
  }

  ## Fit ordinal model
  ord_mod <- tryCatch(expr = {
    lme(ord_form, data = data_ord, random = rand_form, correlation = cor_struct)
    },
    error = function(e) {
      tryCatch(expr = {
        lme(ord_form, data = data_ord, random = rand_form,
            correlation = cor_struct, control = lmeControl(opt='optim'))
        },
        error = function(e) {
          ## Resort to no AR structure
          lme(ord_form, data = data_ord, random = rand_form,
              control = lmeControl(opt='optim'))
      })
  })

  ## Fit continuous model
  cont_mod <- tryCatch(expr = {
    lme(cont_form, data = data_cont, random = rand_form,
        correlation = cor_struct)
  },
  error = function(e) {
    tryCatch(expr = {
      lme(cont_form, data = data_cont, random = rand_form,
          correlation = cor_struct, control = lmeControl(opt='optim'))
    },
    ## Resort to no AR structure
    error = function(w) {
      ## Resort to no AR structure
      lme(cont_form, data = data_cont, random = rand_form,
          control = lmeControl(opt='optim'))
    })
  })

  ## Pass prior matrices
  priors <- as.numeric(c(VarCorr(ord_mod)[, 1][1:N_pat_effects],
                         VarCorr(cont_mod)[, 1][1:N_pat_effects]))
  prior_mat <- diag(priors) * length(priors)
  res_accept <- matrix(NA, nsim, 6)

  for(i in 2:nsim) {
    samp_info$num_iter <- i

    ## Regression coefficients
    vals <- bmrarm_fc_sig_beta(y, X, Z_kron, cur_draws, samp_info)
    res_beta[, i] <- cur_draws$beta <- vals$beta
    res_sig[, i] <- cur_draws$sigma <- vals$sig

    # Autoregressive parameter
    if(samp_info$ar_cov) {
      vals <- bmrarm_mh_ar(y, X, Z_kron, cur_draws, samp_info)
      cur_draws$ar <- vals$ar
      res_accept[i, 2] <- vals$accept
    }
    res_ar[i] <- cur_draws$ar

    ## Subject specific effects
    vals <- bmrarm_fc_patient_siw(y, z, X, cur_draws, samp_info, 1, Z_kron, prior_siw_uni)
    res_pat_sig[, i] <- cur_draws$pat_sig <- vals$pat_sig
    res_pat_eff[,, i] <- cur_draws$pat_effects <- vals$pat_effects
    res_pat_sig_q[,i] <- cur_draws$pat_sig_q <- vals$pat_sig_q
    res_pat_sig_sd[,i] <- cur_draws$pat_sig_sd <- vals$pat_sig_sd
    res_accept[i, 3:6] <- vals$accept_vec

    ## Latent values, missing values, cut points
    y_cuts <- bmrarm_fc_y_cuts(y, z, X, Z_kron, cur_draws, samp_info)
    y <-  res_y[,, i] <- y_cuts$y
    res_cuts[, i] <- cur_draws$cuts <- y_cuts$cuts
    res_accept[i, 1] <- y_cuts$accept

    ## Update missing values
    y <- res_y[,, i]<- bmrarm_fc_missing(y, z, X, Z_kron, cur_draws, samp_info)

    ## Cut points
    #if(i %% 150 == 100) plot(res_cuts[4, ], type = "l")
    if(i %% 150 == 100) plot(res_pat_sig_sd[1, ], type = "l")
    if(i %% 150 == 50 & i > burn_in) print(round(c(colMeans(res_accept[(burn_in+1):nsim,], na.rm = T), i), 3))
    if(i %% 150 == 0) plot(res_pat_sig[1, ], type = "l")
  }

  sim_use <- seq(burn_in + 1, nsim, by = thin)
  draws <- list(
    res_cuts = res_cuts[, sim_use],
    res_beta = res_beta[, sim_use],
    res_pat_sig = res_pat_sig[, sim_use],
    res_pat_sig_q = res_pat_sig_q[, sim_use],
    res_pat_sig_sd = res_pat_sig_sd[, sim_use],
    res_ar = res_ar[sim_use],
    res_sigma = res_sig[, sim_use],
    res_y = res_y[,, sim_use],
    res_pat_eff = res_pat_eff[,, sim_use],
    res_accept = res_accept[sim_use, ],
    samp_info = samp_info,
    X = X,
    Z_kron = Z_kron, z = z,
    priors = priors)
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
#' @importFrom nlme lme VarCorr lmeControl corAR1
#' @export

baseline_bmr_old <- function(formula, data, ordinal_outcome = c("y_ord"),
                         time_var = "time", patient_var = "patient_idx",
                         random_slope = F, ar_cov = TRUE, nsim = 1000,
                         burn_in = 100, thin = 10, seed = 14, verbose = TRUE,
                         sig_prior = 1000000000, sd_vec = c(0.15, 0.30),
                         N_burn_trunc = 5) {

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
    rand_form <- as.formula(paste0("~ ", time_var, "|", patient_var))
  } else {
    rand_form <- as.formula(paste0("~ 1|", patient_var))
  }

  ord_form <- reformulate(deparse(formula[[3]]), response = ordinal_outcome)
  cont_form <- reformulate(deparse(formula[[3]]), response = cont_out_var)

  ## Datasets and structures for models
  data_ord <- data[!is.na(data[, ordinal_outcome]), ]
  data_cont <- data[!is.na(data[, cont_out_var]), ]
  cor_struct <- if(ar_cov) {
    cor_struct <- corAR1(form = as.formula(~ 1 | pat_idx))
  } else {
    cor_struct <- NULL
  }

  ## Fit ordinal model
  ord_mod <- tryCatch(expr = {
    lme(ord_form, data = data_ord, random = rand_form, correlation = cor_struct)
  },
  error = function(e) {
    tryCatch(expr = {
      lme(ord_form, data = data_ord, random = rand_form,
          correlation = cor_struct, control = lmeControl(opt='optim'))
    },
    error = function(e) {
      ## Resort to no AR structure
      lme(ord_form, data = data_ord, random = rand_form,
          control = lmeControl(opt='optim'))
    })
  })

  ## Fit continuous model
  cont_mod <- tryCatch(expr = {
    lme(cont_form, data = data_cont, random = rand_form,
        correlation = cor_struct)
  },
  error = function(e) {
    tryCatch(expr = {
      lme(cont_form, data = data_cont, random = rand_form,
          correlation = cor_struct, control = lmeControl(opt='optim'))
    },
    ## Resort to no AR structure
    error = function(w) {
      ## Resort to no AR structure
      lme(cont_form, data = data_cont, random = rand_form,
          control = lmeControl(opt='optim'))
    })
  })

  ## Pass prior matrices
  priors <- as.numeric(c(VarCorr(ord_mod)[, 1][1:N_pat_effects],
                         VarCorr(cont_mod)[, 1][1:N_pat_effects]))
  prior_mat <- diag(priors) * length(priors)

  for(i in 2:nsim) {
    samp_info$num_iter <- i

    ## Regression coefficients
    vals <- bmrarm_fc_sig_beta(y, X, Z_kron, cur_draws, samp_info)
    res_beta[, i] <- cur_draws$beta <- vals$beta
    res_sig[, i] <- cur_draws$sigma <- vals$sig

    # Autoregressive parameter
    if(samp_info$ar_cov) {
      vals <- bmrarm_mh_ar(y, X, Z_kron, cur_draws, samp_info)
      cur_draws$ar <- vals$ar
      res_accept[i, 2] <- vals$accept
    }
    res_ar[i] <- cur_draws$ar

    ## Subject specific effects
    vals <- bmrarm_fc_patient(y, z, X, cur_draws, samp_info, prior_mat)
    res_pat_sig[, i] <- cur_draws$pat_sig <- vals$pat_sig
    res_pat_eff[,, i] <- cur_draws$pat_effects <- vals$pat_effects

    ## Latent values, missing values, cut points
    y_cuts <- bmrarm_fc_y_cuts(y, z, X, Z_kron, cur_draws, samp_info)
    y <-  res_y[,, i] <- y_cuts$y
    res_cuts[, i] <- cur_draws$cuts <- y_cuts$cuts
    res_accept[i, 1] <- y_cuts$accept

    ## Update missing values
    y <- res_y[,, i]<- bmrarm_fc_missing(y, z, X, Z_kron, cur_draws, samp_info)

    ## Cut points
    if(i %% 150 == 100) plot(res_cuts[4, ], type = "l")
    if(i %% 150 == 50 & i > burn_in) print(colMeans(res_accept[(burn_in+1):nsim,], na.rm = T))
    if(i %% 150 == 0) plot(res_pat_sig[1, ], type = "l")
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
#' @return bmrarm
#' @importFrom loo kfold_split_stratified
#' @export

bmr_cv <- function(formula, data, ordinal_outcome = c("y_ord"),
                         time_var = "time", patient_var = "patient_idx",
                         random_slope = F, ar_cov = TRUE, nsim = 1000,
                         burn_in = 100, thin = 10, seed = 14, verbose = TRUE,
                         sig_prior = 1000000000, sd_vec = c(0.15, 0.30)) {

  ## Long version of dataset
  full_X <- model.matrix.lm(as.formula(formula), data = data, na.action = "na.pass")
  out_vars <- setdiff(all.vars(formula), colnames(full_X))
  cont_out_var <- setdiff(out_vars, ordinal_outcome)
  data$row_num <- 1:nrow(data)
  data$pat_idx_num <- dense_rank(data$pat_idx)
  data$time <- data[, time_var]

  ## Extract matrices for multiplication of random effects
  full_Z <- matrix(rep(1, nrow(data)), ncol = 1)
  if (random_slope) {
    full_Z <- cbind(full_Z, full_X[, time_var, drop = FALSE])
  }
  full_Z_kron <- kronecker(diag(rep(1, length(out_vars))), full_Z)

  ## Max number of observations
  max_obs <- data %>%
    group_by(pat_idx) %>%
    summarise(N = n()) %>%
    ungroup() %>%
    select(N) %>%
    unlist() %>%
    max()

  ## Folds
  set.seed(seed)
  folds <- kfold_split_stratified(max(max_obs, 10), data$pat_idx)
  data$folds <- folds
  cv_mat <- matrix(NA, nrow = nrow(data), ncol = length(out_vars))
  z_true = data[, ordinal_outcome]
  y_true <- data[, cont_out_var]
  i <- 1

  for(i in 1:max(folds)) {
    data_tmp <- data
    cv_locs <- which(data$folds == i)
    data_tmp <- data[data$folds != i, ]
    samps <- baseline_bmr(formula = formula, data = data_tmp,
                          ordinal_outcome = ordinal_outcome,
                          patient_var = patient_var,
                          random_slope = random_slope,
                          time_var = time_var, ar_cov = ar_cov,
                          burn_in = burn_in, nsim = nsim, thin = thin,
                          seed = seed,
                          sd_vec = sd_vec)

    if(!ar_cov) samps$res_ar[] <- 0
    cv_mat[cv_locs, 1] <- get_loocv(samps, cv_locs, z_true, y_true, full_X,
                                    full_Z_kron, data, ar = ar_cov)

    print(i)
    print(colMeans(samps$res_accept))
    print(sum(cv_mat[, 1], na.rm = T) * 2)
  }

  list(cv_vals = cv_mat[, 1], sum_cv_vals = sum(cv_mat[, 1]))
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

baseline_bmr_test <- function(formula, data, ordinal_outcome = c("y_ord"),
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
      vals <- bmrarm_fc_patient(y, z, X, cur_draws, samp_info)
      #print("pat")
      res_pat_sig[, i] <- cur_draws$pat_sig <- vals$pat_sig
      res_pat_eff[,, i] <- cur_draws$pat_effects <- vals$pat_effects
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
    if(i %% 100 == 0) plot(res_sig[4, ], type = "l")
  }

  sim_use <- seq(burn_in + 1, nsim, by = thin)
  draws <- list(
    res_cuts = res_cuts[, sim_use],
    res_beta = res_beta[, sim_use],
    res_pat_sig = res_pat_sig,
    res_ar = res_ar[sim_use],
    res_sigma = res_sig[, sim_use],
    res_y = res_y[,, sim_use],
    res_pat_eff = res_pat_eff[,, sim_use],
    samp_info = samp_info,
    X = X,
    Z_kron = Z_kron, z = z)
  draws
}
