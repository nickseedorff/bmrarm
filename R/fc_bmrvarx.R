#' Prior for the expansion parameters
#'
#' @param cor_mat correlation matrix from previous draw
#' @param N_ordinal number of ordinal outcomes
#' @importFrom LaplacesDemon rinvgamma
#' @return matrix

expansion_prior <- function(cor_mat, N_ordinal) {
  N_total <- ncol(cor_mat)
  N_cont <- N_total - N_ordinal
  cor_inv <- qr.solve(cor_mat)
  expansions <- vector(length = N_ordinal)

  ## Draw from inverse gamma prior
  for(i in 1:N_ordinal) {
    expansions[i] <- sqrt(rinvgamma(1, (N_total + 1) / 2, (cor_inv[i, i]) / 2))
  }

  ## Return diagonal matrix
  diag(c(expansions, rep(1, N_cont)))
}

#' Function to draw regression coefficients and covariance matrix
#'
#' @param y matrix of multivariate continuous observations
#' @param X design matrix
#' @param prior_precision prior precision matrix
#' @param y_orig unscaled version of y
#' @importFrom LaplacesDemon rinvwishart rmatrixnorm
#' @import dplyr
#' @return list

fc_sigma_theta_tilde <- function(y, X, prior_precision, y_orig) {
  N_outcome <- ncol(y)
  w_tmp <- y
  X_tilde <- cbind(X, rbind(rep(0, N_outcome), y_orig[-nrow(y), ]))

  ## Find theta hat
  x_inv <- chol2inv(chol(prior_precision + crossprod(X_tilde)))
  theta_hat <- x_inv %*% t(X_tilde) %*% w_tmp

  ## sigma_draw draw
  val <- crossprod(w_tmp) + diag(rep(1, N_outcome)) -
    t(w_tmp) %*% X_tilde %*% theta_hat
  sig_draw <- rinvwishart(nrow(w_tmp) + N_outcome + 1, val)

  ## effects_draw
  theta <- rmatrixnorm(M = theta_hat, V = sig_draw, U = x_inv)
  list(sigma_tilde = sig_draw, theta_tilde = theta)
}

#' Function to draw cutpoint parameters
#'
#' @param y matrix of multivariate observations
#' @param z matrix of ordinal outcomes
#' @param N_cat vector of positive integers, number of categories for each ordinal outcome
#' @param upper_cut_limit maximum value for the cutpoint priors (delta)
#' @return list

fc_cuts <- function(y, z, N_cat, upper_cut_limit) {
  cuts_list <- list()
  for(i in 1:length(N_cat)) {
    if(N_cat[i] == 2) {
      cuts_list[[i]] <- c(-Inf, 0, Inf)
    } else {
      cuts_list[[i]] <- c(-Inf, 0, rep(NA, N_cat[i] - 2), Inf)
      ## Only update if more than 2 levels
      if(N_cat[i] >= 3) {
        y_late <- y[, i]
        for(j in 3:(length(cuts_list[[i]]) - 1)) {
          cuts_min <- max(y_late[(z[, i] == j - 1)])
          cuts_max <- min(min(y_late[(z[, i] == j)]), upper_cut_limit[i])
          cuts_list[[i]][j] <- runif(1, cuts_min, cuts_max)
        }
      }
    }
  }
  cuts_list
}

#' Function to draw latent continuous values
#'
#' @param y matrix of continuous observations
#' @param z matrix of ordinal outcomes
#' @param tmp_list list of current parameter values
#' @param miss_mat matrix of location of missing values
#' @param samp_info list of internal information used for sampling
#' @param rej_vec vector, values track the attempts by the rejection sampler
#' @return list

fc_y <- function(y, z, mean_mat, tmp_list, miss_mat, samp_info, rej_vec) {

  ## Current parameter values
  sig <- tmp_list$sigma
  M <- tmp_list$M
  cuts <- tmp_list$cuts

  ## Constants of interest
  num_ord <- samp_info$num_ord
  rej_iters_vec <- vector(length = samp_info$N_obs)

  ## Commonly used inverses
  sig_inv <- solve(sig)
  ms_cross <- crossprod(M, sig_inv)
  msm <- ms_cross %*% M

  ## Conditional normal when t = 1 to T-1
  c_mat_def <- solve(sig_inv + msm)

  ## Pre calculation based on locations that need to be sampled
  baseline_mats <- pre_calc(c_mat_def, samp_info$samp_locs)
  final_mats <- pre_calc(sig, samp_info$samp_locs)
  default_mats <- pre_calc(c_mat_def, samp_info$samp_locs)

  for(i in 1:samp_info$N_obs) {
    ## Locations needed for sampling
    iter_samp_type <- samp_info$samp_type[[i]]
    iter_locs <- samp_info$samp_locs[[iter_samp_type]]
    iter_length <- samp_info$samp_length[[iter_samp_type]]

    ## Storage for thresholds
    cuts_low <- rep(-Inf, length = iter_length)
    cuts_high <- rep(Inf, length = iter_length)

    ## Find threshold vectors
    for(j in 1:num_ord) {
      if(miss_mat[i, j] == 0) {
        cuts_low[j] <- cuts[, , j][z[i, j]]
        cuts_high[j] <- cuts[, , j][z[i, j] + 1]
      }
    }

    ## Mean and covariance depend on if the obs is the first, last, or other
    if (i == 1) {
      d_vec <- c_mat_def %*% (sig_inv %*% mean_mat[, i] +
                                ms_cross %*% (y[, i + 1] - mean_mat[, i + 1]))
      pre_calcs <- baseline_mats[[iter_samp_type]]
    } else if (i == samp_info$N_obs) {
      d_vec <- mean_mat[, i] + M %*% y[, i - 1]
      pre_calcs <- final_mats[[iter_samp_type]]
    } else {
      d_vec <- c_mat_def %*% (sig_inv %*% (mean_mat[, i] + M %*% y[, i - 1]) +
                                ms_cross %*% (y[, i + 1] - mean_mat[, i + 1]))
      pre_calcs <- default_mats[[iter_samp_type]]
    }

    ## Store results
    tmp <- tmvn_gibbs_rej(
      y_current = y[, i], mean = d_vec, lower = cuts_low,
      upper = cuts_high, locs = iter_locs, loc_length = iter_length,
      pre_calcs = pre_calcs, max_iter = samp_info$max_iter,
      N_burn_trunc = samp_info$N_burn_trunc)
    y[iter_locs, i] <- tmp$res
    rej_vec[i] <- tmp$rej_iter
  }
  list(y_new = t(y), rej_vec = rej_vec)
}

#' Gibbs step for truncated multivariate normal
#' @param y_current continuous outcome values
#' @param mean mean vector of multivariate normal
#' @param lower scalar, lower truncation bound
#' @param upper scalar, upper truncation bound
#' @param locs vector of locations to sample
#' @param loc_length length of locs
#' @param pre_calcs list of conditional means and covariance matrices
#' @param max_iter integer, maximum number of rejection sampler attempts
#' @param N_burn_trunc integer, number of Gibbs burn-in iterations
#' @importFrom truncnorm rtruncnorm
#' @return list

tmvn_gibbs_rej <- function(y_current, mean, lower, upper, locs, loc_length,
                           pre_calcs, max_iter, N_burn_trunc) {

  ## If missing a single ordinal outcome then sample from truncated normal
  if(loc_length == 1) {
    res <- rtruncnorm(
      n = 1, a = lower[1], b = upper[1], sd = pre_calcs$cond_cov,
      mean = cond_mean_part(y_current, mean, pre_calcs$mean_pre, 1))
    list(res = res, rej_iter = 0)
  } else {
    ## Otherwise use a combination approach
    res <- NA
    iter <- 0
    cond_mean <- cond_mean_part(y_current, mean, pre_calcs$mean_pre, locs)
    while(is.na(res[1]) & iter <= max_iter) {
      y_tmp <- pre_calcs$cond_chol_mat %*% rnorm(loc_length) + cond_mean
      if(all(y_tmp >= lower & y_tmp <= upper)) {
        res <- y_tmp
      }
      iter <- iter + 1
    }

    if(iter > max_iter) {
      start_val <- pmax(pmin(y_current[locs], upper), lower)
      mean_vec <- as.vector(cond_mean_part(y_current, mean, pre_calcs$mean_pre,
                                           locs))
      res <- rtmvnorm(
        1, mean = mean_vec, H = pre_calcs$cond_chol_inv, lower = lower,
        upper = upper, algorithm = "gibbs", burn.in.samples = N_burn_trunc,
        start.value = start_val)
    }
    list(res = res, rej_iter = iter)
  }
}
