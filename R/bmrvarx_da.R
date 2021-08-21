#' PA MCMC routine to sample from BMRVAR model
#'
#' @param y matrix of multivariate observations
#' @param z vector of disease status
#' @param df dataframe with patient and time point indexes
#' @param nsim scalar, number of iterations with default of 1000
#' @param burn_in scalar, number of iterations to remove with default of 50
#' @param seed seed to default, default of 14
#' @return mcmc
#' @importFrom zoo na.approx

bmrvarx_da <- function(formula, data, ordinal_outcomes = c("y_ord", "y_bin"),
                      sig_prior = 1000000, all_draws = FALSE, nsim = 1000,
                      burn_in = 100, thin = 10, seed = 14, verbose = TRUE,
                      max_iter_rej = 500, return_y = FALSE, fast = F, old_prior_y0 = F) {

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
  N_burn_trunc <- 0

  ## Get sampling info, initialize, generate storage
  samp_info <- get_sampling_info(env = environment())
  create_storage(env = environment())
  #tmp_list$cuts[3:4, , 1] <- tmp_list$cuts[3:4,,1] + 2.5
  #if(N_ord == 2) {
  #  tmp_list$cuts[3:4, , 2] <- tmp_list$cuts[3:4, , 2] + 2.5
  #}

  ## Run simulationS
  for(i in 2:nsim) {
    mean_mat <- covars %*% tmp_list$beta

    ## Draw latent variables
    y_use <- res_y[,, i] <- fc_y(
      y = t(y_use), z = y_ord, mean_mat = t(mean_mat), tmp_list = tmp_list,
      miss_mat = miss_mat, samp_info = samp_info, num_iter = i, fast = fast)

    ## Sigma and effects
    sig_theta <- fc_sigma_theta_tilde(y = y_use, X = covars, y_orig = y_use,
                                      prior_precision = samp_info$prior_sig,
                                      old_prior_y0)
    sigma_tilde <- sig_theta$sigma_tilde

    ## M and beta
    beta_tilde <- sig_theta$theta_tilde[1:N_covars, ]
    M_tilde <- t(sig_theta$theta_tilde[(N_covars + 1):(N_covars + N_response), ])

    ## Update cuts
    cuts_tilde <- fc_cuts_da(y = y_use, z = y_ord, N_cat)

    ## Store updated values
    res_M[, i] <- tmp_list$M <- M_tilde
    res_sigma[, i] <- tmp_list$sigma <- sigma_tilde
    res_beta[, i] <- tmp_list$beta <- beta_tilde
    res_latent[, i, ] <-  y_use[, 1:N_ord]

    ## Store thresholds
    for(j in 1:N_ord) {
      tmp_list$cuts[1:(N_cat[j] + 1), , j] <- cuts_tilde[[j]]
      res_cuts[1:(N_cat[j] + 1), i, j] <- tmp_list$cuts[1:(N_cat[j] + 1), , j] <-
        cuts_tilde[[j]]
    }

    res_cuts[3, i, 1] <- tmp_list$cuts[3, , 1] <- 1.3
    res_cuts[3, i, 2] <- tmp_list$cuts[3, , 2] <- 1.2

    ## Iterations update, keep track of working parameter
    if(i %% 50 == 0) print(paste0("iteration = ", i, ""))
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
                rej_accept_rate = t(rej_accept_rate))

  ## Return final observation for future forecasts, posterior draws
  data_for_forecasts <- list(last_y = res_y[samp_info$N_obs, , sim_use])
  structure(list(draws = draws, data_for_forecasts = data_for_forecasts,
                 covars_used = colnames(covars), X = covars,
                 y = data[, c(ordinal_outcomes,
                              setdiff(out_vars, ordinal_outcomes))]),
            class = "bmrvarx")
}

#' Full conditional draws for the threshold parameters
#'
#' @param y matrix of multivariate observations
#' @param z matrix of ordinal outcomes
#' @param N_cat vector of positive integers, number of categories for each ordinal outcome
#' @return list

fc_cuts_da <- function(y, z, N_cat) {
  cuts_list <- list()
  for(i in 1:length(N_cat)) {
    if(N_cat[i] == 3) {
      cuts_list[[i]] <- c(-Inf, 0, 2.5, Inf)
    } else {
      cuts_list[[i]] <- c(-Inf, 0, 2.5, rep(NA, N_cat[i] - 3), Inf)
      ## Only update if more than 2 levels
      if(N_cat[i] >= 4) {
        y_late <- y[, i]
        for(j in 4:(length(cuts_list[[i]]) - 1)) {
          cuts_min <- max(y_late[(z[, i] == j - 1)])
          cuts_max <- min(y_late[(z[, i] == j)])
          cuts_list[[i]][j] <- runif(1, cuts_min, cuts_max)
        }
      }
    }
  }
  cuts_list
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
  N_burn_trunc <- 0

  ## Get sampling info, initialize, generate storage
  samp_info <- get_sampling_info(env = environment())
  create_storage(env = environment())

  ## Run simulation
  for(i in 2:nsim) {

    ## Covariance matrix
    sig_theta <- fc_sigma_theta_tilde_bvar(y = y_use, X = covars, y_orig = y_use,
                                           prior_precision = samp_info$prior_sig,
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

  ## Return final observation for future forecasts, posterior draws
  structure(list(draws = draws,
                 covars_used = colnames(covars), X = covars, y = y_use),
            class = "bmrvarx")
}
