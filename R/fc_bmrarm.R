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
#' @importFrom truncnorm rtruncnorm
#' @importFrom OpenMx omxMnor
#' @return matrix
#' @export

bmrarm_fc_y_cuts <- function(y, z, X, Z_kron, cur_draws, samp_info) {

  ## Mean matrix
  mean_mat <- X %*% cur_draws$beta +
    matrix(rowSums(Z_kron * cur_draws$pat_effects[samp_info$pat_idx_long, ]),
           ncol = samp_info$N_outcomes)

  ## Generate full sigma matrix
  cuts <- cur_draws$cuts
  N_pat <- samp_info$N_pat
  N_cat <- samp_info$N_cat
  N_cuts <- N_cat + 1
  cuts_tmp <- cuts

  ## Calculate conditional mean for each subject
  cond_mean <- rep(NA, samp_info$N_obs)
  if(samp_info$ar_cov) {
    sig_list <- get_sig_list(cur_draws, samp_info)
    for(i in 1:N_pat) {
      ## Locations of vectors to sample
      locs <- samp_info$pat_locs[[i]]
      ind <- samp_info$pat_time_ind[i]
      cond_mean[locs] <- as.numeric(
        mean_mat[locs, 1] + sig_list$mean_pre_list[[ind]] %*%
          (as.numeric(y[locs, -1]) - as.numeric(mean_mat[locs, -1])))
    }
  } else {
    ## Conditional means and covariances
    sig <- cur_draws$sigma
    pre_calcs <- pre_calc_ar(sig, 1)
    cond_sd <- sqrt(as.numeric(pre_calcs$cond_cov))
    cond_mean <- mean_mat[, 1] + as.numeric(pre_calcs$mean_pre) *
      (as.numeric(y[, -1]) - as.numeric(mean_mat[, -1]))
  }

# Update cutpoints --------------------------------------------------------

  if(N_cat > 3) {
    sd_c <- samp_info$sd_c

    ## Propose new cutpoints
    for(i in 4:N_cat) {
      cuts_tmp[i] <- rtruncnorm(1, a = cuts_tmp[i - 1], b = cuts_tmp[i + 1],
                                mean = cuts_tmp[i], sd = sd_c)
    }

    ## Storage of thresholds
    cuts_low_old <- cuts[z]
    cuts_high_old <- cuts[z + 1]
    cuts_low <- cuts_tmp[z]
    cuts_high <- cuts_tmp[z + 1]

    ## Missing ordinal outcomes leads to (-Inf, Inf) for the latent value
    cuts_low_old[is.na(cuts_low_old)] <- -Inf
    cuts_low[is.na(cuts_low)] <- -Inf
    cuts_high_old[is.na(cuts_high_old)] <- Inf
    cuts_high[is.na(cuts_high)] <- Inf

    ## Calculate probabilities for cowles solution
    if(samp_info$ar_cov) {
      prob_vec <- rep(NA, samp_info$N_pats_for_probs)
      for(i in 1:samp_info$N_pats_for_probs) {
        ## Locations of vectors to sample
        pat <- samp_info$pats_for_cut_probs[i]
        ind <- samp_info$pat_time_ind[pat]
        locs <- samp_info$pat_locs[[pat]]

        ## Evaluate the multivariate normal probabilities
        prob_num <- omxMnor(lbound = cuts_low[locs],
                            ubound = cuts_high[locs], mean = cond_mean[locs],
                            covariance = sig_list$cond_cov_list[[ind]])[[1]]
        prob_denom <- omxMnor(lbound = cuts_low_old[locs],
                              ubound = cuts_high_old[locs],
                              mean = cond_mean[locs],
                              covariance = sig_list$cond_cov_list[[ind]])[[1]]
        prob_vec[i] <- prob_num / prob_denom
      }
    } else {
      prob_num <- pnorm((cuts_high - cond_mean) / cond_sd) -
        pnorm((cuts_low - cond_mean) / cond_sd)

      denom_num <- pnorm((cuts_high_old - cond_mean) / cond_sd) -
        pnorm((cuts_low_old - cond_mean) / cond_sd)

      prob_vec <- prob_num / denom_num
    }

    ## MH Step, return if nothing changed
    comp_val <- first_cut_prob(samp_info, cuts, cuts_tmp) * prod(prob_vec)

    ## Accept reject step
    if(comp_val < runif(1)) {
      return(list(y = y, cuts = cuts, accept = 0))
    }
  }

# Update latent values ----------------------------------------------------

  if(samp_info$ar_cov) {
    for(i in 1:N_pat) {
      ## Locations of vectors to sample
      locs <- samp_info$pat_locs[[i]]
      ind <- samp_info$pat_time_ind[i]

      ## Sample from multivariate truncated normal
      start_val <- pmax(pmin(cond_mean[locs], cuts_high[locs]), cuts_low[locs])
      y[locs, 1] <- rtmvnorm(
        1, mean = cond_mean[locs], H = sig_list$cond_cov_inv_list[[ind]],
        lower = cuts_low[locs], upper = cuts_high[locs], algorithm = "gibbs",
        burn.in.samples = 10, start.value = start_val)
    }
  } else {
    N_obs <- length(z)
    y[1:N_obs, 1] <- rtruncnorm(N_obs, a = cuts_low, b = cuts_high,
                                mean = cond_mean, sd = cond_sd)
  }
  list(y = y, cuts = cuts_tmp, accept = 1)
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
#' @return matrix
#' @export

bmrarm_fc_missing <- function(y, z, X, Z_kron, cur_draws, samp_info) {

  ## Mean vector
  mean_vec <- as.numeric(X %*% cur_draws$beta) +
    rowSums(Z_kron * cur_draws$pat_effects[samp_info$pat_idx_long, ])
  y_vec <- as.numeric(y)

  if(samp_info$ar_cov) {
    sig_list <- get_sig_list(cur_draws, samp_info)
  }

  ## Generate new values
  for(i in 1:samp_info$N_pat) {
    if(length(samp_info$pat_cont_miss_rank[[i]]) > 0) {

      ## Observation times and missing locations
      all_locs <- samp_info$pat_all_locs[[i]]
      miss_locs <- samp_info$pat_cont_miss_locs[[i]]
      miss_ranks <- samp_info$pat_cont_miss_rank[[i]]
      pat_y <- y_vec[all_locs]
      pat_mean <- mean_vec[all_locs]
      ind <- samp_info$pat_time_ind[i]

      ## Generate condition mean and covariances
      if(samp_info$ar_cov) {
        full_sig <- sig_list$sig_list[[ind]]
      } else {
        full_sig <- kronecker(cur_draws$sigma,
                              diag(rep(1, samp_info$pat_N_obs[i])))
      }
      pre_calcs <- pre_calc_ar(full_sig, locs = miss_ranks)

      ## Sample from multivariate normal
      tmp_mean <- as.numeric(pat_mean[miss_ranks] + pre_calcs$mean_pre %*%
                               (pat_y[-miss_ranks] - pat_mean[-miss_ranks]))
      L <- t(chol(pre_calcs$cond_cov))
      y_vec[miss_locs] <- L %*% rnorm(length(tmp_mean)) + tmp_mean
    }
  }
  matrix(y_vec, ncol = samp_info$N_outcomes)
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
#' @importFrom magic adiag
#' @return matrix
#' @export

bmrarm_fc_patient <- function(y, z, X, cur_draws, samp_info, prior_mat) {

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

  ## Covariance matrix for random effects
  pat_sig <- rinvwishart(N_pat + N_pat_eff, crossprod(res) + prior_mat)
  list(pat_effects = res, pat_sig = pat_sig)
}

#' Full conditional draws of the regression coefficients
#'
#' @param subject_effects matrix of patient specific intercepts
#' @param sigma residual covariance matrix
#' @param prior_alpha prior term for shape
#' @param prior_alpha prior term for scale
#' @importFrom LaplacesDemon rinvwishart rmatrixnorm
#' @export

bmrarm_fc_sig_beta <- function(y, X, Z_kron, cur_draws, samp_info) {

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
  prior <- diag(rep(0.00001, N_covars))
  x_inv <- chol2inv(chol(prior + cov_vals))
  beta_hat <- x_inv %*% mean_vals
  val <- resid_vals + diag(rep(1, N_outcomes)) - t(mean_vals) %*% beta_hat

  ## Draw updates
  sig <- rinvwishart(samp_info$N_obs + N_outcomes, val)
  beta <- rmatrixnorm(M = beta_hat, V = sig, U = x_inv)
  list(beta = beta, sig = sig)
}

#' Full conditional draws of the regression coefficients
#'
#' @param subject_effects matrix of patient specific intercepts
#' @param sigma residual covariance matrix
#' @param prior_alpha prior term for shape
#' @param prior_alpha prior term for scale
#' @return scalar
#' @export

bmrarm_mh_ar <- function(y, X, Z_kron, cur_draws, samp_info) {
  ## Mean vector needed for all covariance terms
  resid_mat <- y -  X %*% cur_draws$beta -
    matrix(rowSums(Z_kron * cur_draws$pat_effects[samp_info$pat_idx_long, ]),
           ncol = samp_info$N_outcomes)

  ## Propose new values
  cur_draws2 <- cur_draws
  cur_draws2$ar <- rnorm(1, cur_draws$ar, sd = samp_info$sd_ar)
  if(abs(cur_draws2$ar) >=1 ) return(list(ar = cur_draws$ar, accept = 0))

  ## Calculate comparison values
  sig_list_old <- get_sig_list(cur_draws, samp_info)
  comp_old <- dmatrix_normal_log(resid_mat, cur_draws, samp_info, sig_list_old)
  sig_list_new <- get_sig_list(cur_draws2, samp_info)
  comp_new <- dmatrix_normal_log(resid_mat, cur_draws2, samp_info, sig_list_new)

  ## Get comparison value
  compar_val <- comp_new - comp_old

  ## Accept proposal
  if(compar_val >= log(runif(1))) {
    return(list(ar = cur_draws2$ar, accept = 1))
  } else {
    return(list(ar = cur_draws$ar, accept = 0))
  }
}

#' Full conditional draws of the regression coefficients
#'
#' @param subject_effects matrix of patient specific intercepts
#' @param sigma residual covariance matrix
#' @param prior_alpha prior term for shape
#' @param prior_alpha prior term for scale
#' @return scalar
#' @export

dmatrix_normal_log <- function(resid_mat, cur_draws, samp_info, sig_list) {
  N_pat <- samp_info$N_pat
  pat_vals <- rep(NA, N_pat)
  pat_corr_mat <- list()
  sig_inv <- chol2inv(chol(cur_draws$sigma))
  for(i in 1:N_pat) {
    ## Get patient locations and differences
    pat_resid <- resid_mat[samp_info$pat_locs[[i]], ]
    time_ind <- samp_info$pat_time_ind[i]
    time_inv <- sig_list$time_inv[[time_ind]]
    time_det <- sig_list$time_det[[time_ind]]

    ## Contribution from each patient
    pat_vals[i] <- -0.5 * sum(diag(sig_inv %*% crossprod(pat_resid, time_inv)
                                   %*% pat_resid)) -
      samp_info$N_outcomes / 2 * time_det
  }
  sum(pat_vals)
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

bmrarm_fc_patient_siw <- function(y, z, X, cur_draws, samp_info, prior_list, Z_kron, prior_siw_uni) {

  ## Generate full sigma matrix
  N_pat <- samp_info$N_pat
  sig_alpha_inv <- chol2inv(chol(cur_draws$pat_sig_q))
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
    pat_Z <- samp_info$pat_z_kron[[i]] %*% diag(cur_draws$pat_sig_sd)

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
  cur_draws$pat_sig_q <- rinvwishart(N_pat + N_pat_eff, crossprod(res) + diag(rep(1, N_pat_eff)))

  ## SD parameters
  accept_vec <- rep(0, N_pat_eff)
  for(i in 1:N_pat_eff) {
    ## Propose new values
    cur_draws2 <- cur_draws
    cur_draws2$pat_sig_sd[i] <- rnorm(1, cur_draws$pat_sig_sd[i],
                                      sd = samp_info$sd_pat_sd[i])
    cur_draws2$pat_sig_sd[i] <- rtruncnorm(1, a = prior_siw_uni[1],
                                           b = prior_siw_uni[2],
                                           mean = cur_draws$pat_sig_sd[i],
                                           sd = samp_info$sd_pat_sd[i])


    resid_mat_old <- y -  X %*% cur_draws$beta -
      matrix(rowSums(Z_kron * (res[samp_info$pat_idx_long, ] %*%
                                 diag(cur_draws$pat_sig_sd))),
             ncol = samp_info$N_outcomes)

    resid_mat_new <- y -  X %*% cur_draws$beta -
      matrix(rowSums(Z_kron * (res[samp_info$pat_idx_long, ] %*%
                                 diag(cur_draws2$pat_sig_sd))),
             ncol = samp_info$N_outcomes)

    ## Calculate comparison values
    sig_list <- get_sig_list(cur_draws, samp_info)
    comp_old <- dmatrix_normal_log(resid_mat_old, cur_draws, samp_info, sig_list)
    comp_new <- dmatrix_normal_log(resid_mat_new, cur_draws2, samp_info, sig_list)
    compar_val <- comp_new - comp_old +
      log(truncnorm::dtruncnorm(cur_draws$pat_sig_sd[i], prior_siw_uni[1],
                                prior_siw_uni[2],
                                mean = cur_draws2$pat_sig_sd[i],
                                sd = samp_info$sd_pat_sd[i])) -
      log(truncnorm::dtruncnorm(cur_draws2$pat_sig_sd[i], prior_siw_uni[1],
                                prior_siw_uni[2],
                                mean = cur_draws$pat_sig_sd[i],
                                sd = samp_info$sd_pat_sd[i]))

    if(compar_val >= log(runif(1)) & cur_draws2$pat_sig_sd[i] >= -Inf &
       cur_draws2$pat_sig_sd[i] <= Inf) {
      cur_draws$pat_sig_sd[i] <- cur_draws2$pat_sig_sd[i]
      accept_vec[i] <- 1
    }
  }

  list(pat_effects = res %*% diag(cur_draws$pat_sig_sd),
       pat_sig = diag(cur_draws$pat_sig_sd) %*% cur_draws$pat_sig_q %*% diag(cur_draws$pat_sig_sd),
       pat_sig_sd = cur_draws$pat_sig_sd, accept_vec = accept_vec,
       pat_sig_q = cur_draws$pat_sig_q)
}
