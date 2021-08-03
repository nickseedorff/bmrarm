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

bmrarm <- function(formula, data, ordinal_outcome = c("y_ord"),
                             time_var = "time", patient_var = "patient_idx",
                             random_slope = F, ar_cov = TRUE, nsim = 1000,
                             burn_in = 100, thin = 10, seed = 14, verbose = TRUE,
                             sig_prior = 1000000000, sd_vec = c(0.15, 0.30),
                             N_burn_trunc = 5, prior_siw_uni = c(0.2, 5)) {

  ## Create storage
  set.seed(seed)
  bmrarm_start(env = environment())
  cont_out_var <- setdiff(out_vars, ordinal_outcome)

  ## Starting values for ordinal outcome
  y[, 1] <- (res_cuts[z, 1] + res_cuts[z + 1, 1]) / 2
  y[is.infinite(y[, 1]) & y[, 1] < 0, 1] <- -0.5
  y[is.infinite(y[, 1]) & y[, 1] > 0, 1] <-
    max(cur_draws$cuts[samp_info$N_cat]) + 0.5
  y[is.na(y[, 1]), 1] <- 0

  ## Starting values for continuous outcomes
  df <- data.frame(patient = samp_info$pat_idx_long, y = as.numeric(y),
                   outcome = rep(1:N_outcomes, each = N_obs)) %>%
    filter(outcome != 1) %>%
    group_by(patient, outcome) %>%
    mutate(y_interp = na.approx(y, na.rm = FALSE),
           y_interp = ifelse(!is.na(y_interp), y_interp,
                             ifelse(row_number() == n(),
                                    lag(y_interp), lead(y_interp))))
  y[, 2:ncol(y)] <- matrix(df$y_interp, ncol = N_outcomes - 1)
  y[is.na(y)] <- 0
  res_accept <- matrix(NA, nsim, 6)

  for(i in 2:nsim) {
    samp_info$num_iter <- i

    ## Regression coefficients
    vals <- bmrarm_fc_sig_beta(y, X, Z_kron, cur_draws, samp_info)
    res_beta[, i] <- cur_draws$beta <- vals$beta
    res_sig[, i] <- cur_draws$sigma <- vals$sig

    ## Autoregressive parameter
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

    ## Latent values and cut points
    y_cuts <- bmrarm_fc_y_cuts(y, z, X, Z_kron, cur_draws, samp_info)
    y <-  res_y[,, i] <- y_cuts$y
    res_cuts[, i] <- cur_draws$cuts <- y_cuts$cuts
    res_accept[i, 1] <- y_cuts$accept

    ## Update missing values
    y <- res_y[,, i]<- bmrarm_fc_missing(y, z, X, Z_kron, cur_draws, samp_info)

    ## Cut points
    #if(i %% 150 == 100) plot(res_cuts[4, ], type = "l")
    #if(i %% 150 == 100) plot(res_pat_sig_sd[1, ], type = "l")
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
    Z_kron = Z_kron, z = z)
  draws
}
