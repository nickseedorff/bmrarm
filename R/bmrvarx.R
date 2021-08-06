#' PX-DA MCMC routine to sample from BMRVAR model
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
#' @export


bmrvarx <- function(formula, data, ordinal_outcomes = c("y_ord", "y_bin"),
                    sig_prior = 1000000, all_draws = FALSE, nsim = 1000,
                    burn_in = 100, thin = 10, seed = 14, verbose = TRUE,
                    max_iter_rej = 500, return_y = FALSE, fast = T, old_prior_y0 = F) {

  ## Extract outcome variables, record missing values
  out_vars <- setdiff(all.vars(formula),
                      attr(terms(formula), which = "term.labels"))
  dat_out <- data[, c(ordinal_outcomes, setdiff(out_vars, ordinal_outcomes))]
  y_ord <- as.matrix(dat_out[, ordinal_outcomes, drop = FALSE])
  y_ord[is.na(y_ord)] <- -10000
  miss_mat <- matrix(as.numeric(is.na(dat_out)), ncol = length(out_vars))

  ## Starting values are linearly interpolated, initial or final values set to 0
  dat_out <- apply(dat_out, 2, function(x) na.approx(x, na.rm = F))
  dat_out[is.na(dat_out)] <- 0
  data[, colnames(dat_out)] <- dat_out

  ## Set seed and get constants
  set.seed(seed)
  N_cat <- apply(y_ord, 2, function(x) length(unique(x)) - (min(x) == -10000))
  N_ord <- ncol(y_ord)
  N_response <- ncol(dat_out)
  N_cont <- N_response - N_ord
  N_obs <- nrow(dat_out)
  covars <- model.matrix(as.formula(formula), data = data)
  N_covars <- ncol(covars)
  N_base_covars <- 0
  pat_idx <- rep(1, N_obs)

  ## Get sampling info, initialize, generate storage
  samp_info <- get_sampling_info(env = environment())
  create_storage(env = environment())

  ## Run simulation
  for(i in 2:nsim) {
    mean_mat <- covars %*% tmp_list$beta

    ## Draw latent variables
    if(!old_prior_y0) {
      y_use <- res_y[,, i] <- fc_y(
        y = t(y_use), z = y_ord, mean_mat = t(mean_mat), tmp_list = tmp_list,
        miss_mat = miss_mat, samp_info = samp_info, num_iter = i, fast = fast)
    } else {
      covars[1, ] <- 0
      y_use <- res_y[,, i] <- fc_y_old(
        y = t(y_use), z = y_ord, mean_mat = t(mean_mat), tmp_list = tmp_list,
        miss_mat = miss_mat, samp_info = samp_info, num_iter = i, fast = fast)
    }

    ## Draw expansion parameters from prior and transform data
    D_prior <- expansion_prior(cor_mat, N_ordinal = N_ord)
    w_use <- y_use %*% D_prior

    ## Covariance matrix
    sig_theta <- fc_sigma_theta_tilde(y = w_use, X = covars, y_orig = y_use,
                                      prior_precision = samp_info$prior_non_base,
                                      old_prior_y0)
    sigma_tilde <- sig_theta$sigma_tilde

    ## M and beta
    beta_tilde <- sig_theta$theta_tilde[1:N_covars, ]
    M_tilde <- t(sig_theta$theta_tilde[(N_covars + 1):(N_covars + N_response), ])

    ## Update cuts
    upper_cut_limit <- sqrt(diag(sigma_tilde)) * 10000
    cuts_tilde <- fc_cuts(y = w_use, z = y_ord, N_cat, upper_cut_limit)

    ## Values need for transformations
    diag_vals <- sqrt(diag(sigma_tilde[1:N_ord, 1:N_ord, drop = FALSE]))
    v_half <- diag(c(diag_vals, rep(1, N_cont)))
    v_half_inv <- diag(1 / diag(v_half))
    v_short_inv <- v_half_inv[1:N_ord, 1:N_ord, drop = FALSE]

    ## Store updated values
    res_M[, i] <- tmp_list$M <- v_half_inv %*% M_tilde
    res_sigma[, i] <- tmp_list$sigma <- v_half_inv %*% sigma_tilde %*% v_half_inv
    res_beta[, i] <- tmp_list$beta <- beta_tilde %*% v_half_inv
    y_use[, 1:N_ord] <- res_latent[, i, ] <-  w_use[, 1:N_ord] %*% v_short_inv
    res_y[,, i] <- y_use

    ## Store thresholds
    for(j in 1:N_ord) {
      res_cuts[1:(N_cat[j] + 1), i, j] <- tmp_list$cuts[1:(N_cat[j] + 1), , j] <-
        cuts_tilde[[j]] * v_half_inv[j, j]
    }

    ## Store correlation matrix
    cor_mat <- cov2cor(sigma_tilde)

    ## Iterations update, keep track of working parameter
    if(i %% 50 == 0) print(paste0("iteration = ", i, "; unident = ",
                                  paste(round(diag(v_half), 3), collapse = "; ")))
  }


  sim_use <- seq(burn_in + 1, nsim, by = thin)
  all <- list(res_M = res_M,
              res_sigma = res_sigma,
              res_beta = res_beta,
              res_cuts = res_cuts)

  draws <- list(res_M = res_M[, sim_use],
                res_sigma = res_sigma[, sim_use],
                res_beta = res_beta[, sim_use],
                res_cuts = res_cuts[, sim_use, ],
                res_y = res_y[,, sim_use],
                rej_accept_rate = t(rej_accept_rate))

  data_for_forecasts <- list(last_y = res_y[samp_info$last_obs_num, , sim_use])

  if(all_draws) {
    lst <- list(all = all, draws = draws, data_for_forecasts = data_for_forecasts,
                covars_used = colnames(covars))
  } else {
    lst <- list(draws = draws, data_for_forecasts = data_for_forecasts,
                covars_used = colnames(covars))
  }

  if(return_y) {
    lst[["y_draws"]] <- res_y[,, sim_use]
  }

  lst
}

#' PX-DA MCMC routine to sample from BMRVAR model
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
#' @export


bvar <- function(formula, data, sig_prior = 1000000, all_draws = FALSE,
                 nsim = 1000, burn_in = 100, thin = 10, seed = 14,
                 verbose = TRUE) {

  ## Extract outcome variables, record missing values
  out_vars <- setdiff(all.vars(formula),
                      attr(terms(formula), which = "term.labels"))
  dat_out <- data[, out_vars]
  miss_mat <- matrix(as.numeric(is.na(dat_out)), ncol = length(out_vars))

  ## Set seed and get constants
  set.seed(seed)
  N_response <- ncol(dat_out)
  N_obs <- nrow(dat_out)
  covars <- model.matrix(as.formula(formula), data = data)
  N_covars <- ncol(covars)
  pat_idx <- rep(1, N_obs)
  N_cat <- 4
  N_ord <- 0
  N_cont <- N_response - N_ord
  covars <- model.matrix(as.formula(formula), data = data)
  N_covars <- ncol(covars)
  N_base_covars <- 0
  max_iter_rej <- 0

  ## Get sampling info, initialize, generate storage
  samp_info <- get_sampling_info(env = environment())
  create_storage(env = environment())

  ## Run simulation
  for(i in 2:nsim) {

    ## Covariance matrix
    sig_theta <- fc_sigma_theta_tilde(y = y_use, X = covars, y_orig = y_use,
                                      prior_precision = samp_info$prior_non_base,
                                      old_prior_y0 = T)
    sigma_tilde <- sig_theta$sigma_tilde

    ## M and beta
    beta_tilde <- sig_theta$theta_tilde[1:N_covars, ]
    M_tilde <- t(sig_theta$theta_tilde[(N_covars + 1):(N_covars + N_response), ])

    ## Store updated values
    res_M[, i] <- tmp_list$M <- M_tilde
    res_sigma[, i] <- tmp_list$sigma <- sigma_tilde
    res_beta[, i] <- tmp_list$beta <- beta_tilde

    ## Iterations update, keep track of working parameter
    if(i %% 50 == 0) print(paste0("iteration = ", i, ";"))
  }


  sim_use <- seq(burn_in + 1, nsim, by = thin)
  all <- list(res_M = res_M,
              res_sigma = res_sigma,
              res_beta = res_beta)

  draws <- list(res_M = res_M[, sim_use],
                res_sigma = res_sigma[, sim_use],
                res_beta = res_beta[, sim_use])


  if(all_draws) {
    lst <- list(all = all, draws = draws,
                covars_used = colnames(covars))
  } else {
    lst <- list(draws = draws,
                covars_used = colnames(covars))
  }

  lst
}

