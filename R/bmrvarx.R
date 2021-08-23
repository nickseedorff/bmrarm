#' PX-DA MCMC routine to implement a bmrvarx model
#'
#' @param formula an object of class "formula"; a symbolic description of the model to be fitted
#' @param data a dataframe containing outcome variables, covariates, and a patient or subject identifier
#' @param ordinal_outcomes a character string containing the names of the ordinal outcomes
#' @param sig_prior scalar, prior variance on the regression coefficients
#' @param nsim positive integer, number of iterations with default of 1000
#' @param burn_in positive integer, number of iterations to remove with default of 100
#' @param thin positive integer, specifiers the period of saving samples. Default of 10
#' @param seed positive integer, seed for random number generation
#' @param max_iter_rej maximum number of rejection algorithm attempts for multivariate truncated normal
#' @param N_burn_trunc integer, number of burn-in draws from the truncated multivariate normal Gibbs sampler
#' @importFrom zoo na.approx
#' @return bmrvarx
#' @export

bmrvarx <- function(formula, data, ordinal_outcomes,
                    sig_prior = 1000000, nsim = 1000, burn_in = 100, thin = 10,
                    seed = 14, max_iter_rej = 500, N_burn_trunc = 10) {

  ## Extract outcome variables, record missing values
  out_vars <- setdiff(all.vars(formula), labels(terms(formula)))
  dat_out <- data[, c(ordinal_outcomes, setdiff(out_vars, ordinal_outcomes))]
  y_ord <- as.matrix(dat_out[, ordinal_outcomes, drop = FALSE])

  # Stopping rules ----------------------------------------------------------

  ## Ensure ordinal outcomes are based on equally spaced integers
  ord_correct_vec <- rep(NA, length = ncol(y_ord))
  for(i in 1:length(ord_correct_vec)) {
    z <- y_ord[, i]
    ord_correct_vec[i] <- all(min(z, na.rm = T):max(z, na.rm = T) ==
                                1:length(setdiff(unique(z), NA)))
  }

  if(!all(ord_correct_vec)) {
    stop("
    The ordinal outcomes must be integer valued, start at 1, be incremented by
    1, and no integers can be missing. For example, a 5 level ordinal outcome
    must take the values 1, 2, 3, 4, or 5. There must be at least one
    observation for each value.")
  }

  # -------------------------------------------------------------------------

  ## Starting values are linearly interpolated, initial or final values set to 0
  y_ord[is.na(y_ord)] <- -10000
  miss_mat <- matrix(as.numeric(is.na(dat_out)), ncol = length(out_vars))
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

  ## Get sampling info, initialize, generate storage
  samp_info <- get_sampling_info(env = environment())
  create_storage(env = environment())
  seqq <- seq(0, nsim, length.out = 11)
  rej_sampler_mat <- matrix(NA, nrow = nrow(y_use), ncol = nsim)
  rej_vec <- rep(NA, nrow(y_use))

  ## Run simulation
  for(i in 2:nsim) {
    mean_mat <- covars %*% tmp_list$beta

    ## Draw latent variables
    y_res <- fc_y(
      y = t(y_use), z = y_ord, mean_mat = t(mean_mat), tmp_list = tmp_list,
      miss_mat = miss_mat, samp_info = samp_info, rej_vec)
    y_use <- res_y[,, i] <- y_res$y_new
    rej_sampler_mat[, i] <- y_res$rej_vec

    ## Draw expansion parameters from prior and transform data
    D_prior <- expansion_prior(cor_mat, N_ordinal = N_ord)
    w_use <- y_use %*% D_prior

    ## Covariance matrix
    sig_theta <- fc_sigma_theta_tilde(y = w_use, X = covars, y_orig = y_use,
                                      prior_precision = samp_info$prior_sig)
    sigma_tilde <- sig_theta$sigma_tilde

    ## M and beta
    beta_tilde <- sig_theta$theta_tilde[1:N_covars, ]
    M_tilde <- t(sig_theta$theta_tilde[(N_covars + 1):(N_covars + N_response),])

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
    if(i %in% seqq) cat(paste0("Iteration = ", i, "; expansion parameters = ",
                               paste(round(diag(v_half)[1:N_ord], 3),
                                     collapse = "; ")), "\n")
  }

  ## Return draws after burn in
  sim_use <- seq(burn_in + 1, nsim, by = thin)
  draws <- list(res_M = res_M[, sim_use],
                res_sigma = res_sigma[, sim_use],
                res_beta = res_beta[, sim_use],
                res_cuts = res_cuts[, sim_use, ],
                res_y = res_y[,, sim_use])

  ## Return final observation for future forecasts, posterior draws
  data_for_forecasts <- list(last_y = res_y[samp_info$N_obs, , sim_use])
  structure(list(draws = draws, data_for_forecasts = data_for_forecasts,
                 covars_used = colnames(covars), X = covars,
                 rej_sampler_tracker = t(rej_sampler_mat),
                 y = data[, c(ordinal_outcomes,
                              setdiff(out_vars, ordinal_outcomes))]),
            class = "bmrvarx")
}

