#' MCMC sampler to implement a bmrarm model
#'
#' @param formula an object of class "formula"; a symbolic description of the model to be fitted
#' @param data a dataframe containing outcome variables, covariates, and a patient or subject identifier
#' @param ordinal_outcome a character string, which contains the variable name for the ordinal outcome
#' @param time_var a character string, which contains the variable name for time indexing;
#' this must be integer valued and should indicate the observation number.
#' @param patient_var a character string, which contains the variable name for patient indexing
#' @param random_Slope logical, indicating use of random slopes. The default is TRUE
#' @param ar_cov logical, indicating use of an autoregressive error term. The default is TRUE
#' @param nsim positive integer, number of iterations with default of 1000
#' @param burn_in positive integer, number of iterations to remove with default of 100
#' @param thin positive integer, specifies the period of saving samples. Default of 10
#' @param seed positive integer, seed for random number generation
#' @param sig_prior scalar, variance term for prior on the beta coefficients
#' @param sd_vec numeric vector, containing standard deviations for Metroplis Hastings proposals, either of length 4 or 6.
#' The values are associated with the MH-within-gibbs step, the ar term, and the random effects
#' @param N_burn_trunc integer, number of burn-in draws from the truncated multivariate normal Gibbs sampler
#' @param prior_siw_uni prior bounds for the uniform distribution associated with the SIW prior
#' @param prior_siw_df degrees of freedom for the SIW prior.
#' The default is 1 + number of random effects for a single person (2 or 4).
#' @param prior_siw_scale_mat scale matrix for the SIW prior. The default is an identity matrix
#' @return bmrarm
#' @importFrom zoo na.approx
#' @export

bmrarm <- function(formula, data, ordinal_outcome,
                   time_var, patient_var, random_slope = T, ar_cov = TRUE,
                   nsim = 1000, burn_in = 100, thin = 10, seed = 14,
                   sig_prior = 1000000, sd_vec = c(0.15, 0.30, rep(0.2, 4)),
                   N_burn_trunc = 10, prior_siw_uni = c(0.2, 5),
                   prior_siw_df = NULL, prior_siw_scale_mat = NULL) {

  ## Create storage
  set.seed(seed)
  bmrarm_start(env = environment())
  y_store <- y
  cont_out_var <- setdiff(out_vars, ordinal_outcome)

  # Stopping rules ----------------------------------------------------------

  if(nsim <= burn_in) {
    stop("nsim must be > than burn_in")
  }

  ## Current implementation has only tested one ordinal and one continuous
  if(length(ordinal_outcome) > 1) {
    stop("brmarm only allows one ordinal outcome")
  }

  if(length(cont_out_var) != 1) {
    stop("brmarm must be supplied one continous outcome")
  }

  if(random_slope & length(sd_vec) != 6) {
    stop("Length of sd_vec must = 6, for the 6 metropolis hastings proposals")
  }

  if(!random_slope & length(sd_vec) != 4) {
    stop("Length of sd_vec must = 4, for the 4 metropolis hastings proposals")
  }

  ## Ensure ordinal outcome is based on equally spaced integers
  if(!all(min(z, na.rm = T):max(z, na.rm = T) ==
          1:length(setdiff(unique(z), NA)))) {
    stop("
    The ordinal outcome must be integer valued, start at 1, be incremented by 1,
    and no integers can be missing. For example, a 5 level ordinal outcome must
    take the values 1, 2, 3, 4, or 5. There must be at least one observation for
    each value.")
  }

  # -------------------------------------------------------------------------

  ## SIW priors
  num_eff <- ncol(samp_info$pat_z_kron[[1]])
  if(is.null(prior_siw_df)) prior_siw_df <- num_eff + 1
  if(is.null(prior_siw_scale_mat)) prior_siw_scale_mat <- diag(num_eff)

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

  ## Storage for MH acceptance
  res_accept <- matrix(NA, nsim, 6)
  loc_accept <- 1:(4 + 2 * random_slope)
  seqq <- seq(0, nsim, length.out = 11)

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
    vals <- bmrarm_fc_patient_siw(y, z, X, cur_draws, samp_info, Z_kron,
                                  prior_siw_uni, prior_siw_df,
                                  prior_siw_scale_mat)
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

    ## Print output
    if(i %in% seqq) {
      if(i > burn_in) {
        cat("Iteration = ", i, "\n")
        vec <- round(colMeans(res_accept[(burn_in+1):nsim, loc_accept],
                              na.rm = T), 2)
        cat("MH Acceptance ratios: \u03b3 = ", vec[1], "; ",
            "\u03c1 = ", vec[2], "; ",
            "\u03be = ", paste(vec[3:length(vec)], collapse = ", "), "\n", sep = "")
      } else {
        cat("Iteration = ", i, "\n")
      }
    }
  }

  ## Return posterior draws, data, sampling info
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
    res_accept = res_accept[sim_use, loc_accept])

  data <- list(samp_info = samp_info, X = X,
               Z_kron = Z_kron, z = z, y = y_store[, 2])
  structure(list(draws = draws, data = data), class = "bmrarm")
}
