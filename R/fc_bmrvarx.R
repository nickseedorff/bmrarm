#' Full conditional draws of the latent continuous values
#'
#' @param cor_mat correlation matrix from previous draw
#' @param N_ordinal number of ordinal outcomes
#' @return matrix
#' @importFrom LaplacesDemon rinvgamma

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

#' Full conditional draws of the regression coefficients
#'
#' @param y matrix of multivariate observations
#' @param X design matrix
#' @param prior_precision prior precision matrix
#' @return matrix
#' @importFrom LaplacesDemon rinvwishart rmatrixnorm
#' @import dplyr
#' @export

fc_sigma_theta_tilde <- function(y, X, prior_precision, y_orig, old_prior_y0) {
  N_outcome <- ncol(y)
  w_tmp <- y
  X_tilde <- cbind(X, rbind(rep(0, N_outcome), y_orig[-nrow(y), ]))

  if(old_prior_y0) {
    w_tmp <- y[-1, ]
    X_tilde <- cbind(X, rbind(rep(0, N_outcome), y_orig[-nrow(y), ]))[-1, ]
  }


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

#' Full conditional draws for the threshold parameters
#'
#' @param y matrix of multivariate observations
#' @param z matrix of ordinal outcomes
#' @param N_cat vector of positive integers, number of categories for each ordinal outcome
#' @return list
#' @export

fc_cuts <- function(y, z, N_cat) {
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
          cuts_max <- min(y_late[(z[, i] == j)])
          cuts_list[[i]][j] <- runif(1, cuts_min, cuts_max)
        }
      }
    }
  }
  cuts_list
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

fc_y <- function(y, z, mean_mat, tmp_list, miss_mat, samp_info, num_iter, fast) {

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
    if (fast) {
      y[iter_locs, i] <- tmvn_gibbs_rej_fast(
        y_current = y[, i], mean = d_vec, lower = cuts_low,
        upper = cuts_high, locs = iter_locs, loc_length = iter_length,
        pre_calcs = pre_calcs, max_iter = samp_info$max_iter, N_ord = num_ord,
        burn_in = samp_info$burn_in, num_iter = num_iter, num_obs = i)
    } else {
      y[iter_locs, i] <- tmvn_gibbs_rej(
        y_current = y[, i], mean = d_vec, lower = cuts_low,
        upper = cuts_high, locs = iter_locs, loc_length = iter_length,
        pre_calcs = pre_calcs, max_iter = samp_info$max_iter, N_ord = num_ord,
        burn_in = samp_info$burn_in, num_iter = num_iter, num_obs = i)
    }
    #rej_iters_vec[i] <- rej_samp_iters
  }
  #rej_iters_tmp <<- c(median(rej_iters_vec), max(rej_iters_vec))
  t(y)
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

fc_y_old <- function(y, z, mean_mat, tmp_list, miss_mat, samp_info, num_iter, fast) {

  ## Current parameter values
  N_ord <- samp_info$num_ord
  N_cont <- samp_info$N_response - N_ord
  ord_loc <- 1:N_ord
  sig <- tmp_list$sigma
  M <- tmp_list$M
  cuts <- tmp_list$cuts

  ## Constants of interest
  num_ord <- samp_info$num_ord
  rej_iters_vec <- vector(length = samp_info$N_obs)

  ## Commonly used inverses
  sig_inv <- chol2inv(chol(sig))
  ms_cross <- crossprod(M, sig_inv)
  msm <- ms_cross %*% M

  ## Conditional normal when t = 1 to T-1
  c_mat_def <- chol2inv(chol(sig_inv + msm))
  ms0_cross <- ms_cross[ord_loc, , drop = FALSE]
  c_mat0 <- chol2inv(chol(diag(N_ord) + ms0_cross %*% M[, ord_loc, drop = FALSE]))

  ## Pre calculation based on locations that need to be sampled
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
      d_vec <- c(c_mat0 %*% ms0_cross %*% (
        y[, i + 1] - mean_mat[, i + 1] -
          M[, -(ord_loc), drop = FALSE]  %*% y[-(ord_loc), i, drop = F]
      ), rep(0, N_cont))
      #pre_calcs <- baseline_mats[[iter_samp_type]]
      pre_calcs <- list(
        mean_pre = matrix(0, nrow = N_ord, ncol = N_cont),
        cond_cov = c_mat0,
        cond_chol_mat = t(chol(c_mat0)),
        cond_chol_inv = chol2inv(chol(c_mat0)))
      if(N_ord == 1) {
        pre_calcs$cond_cov <- sqrt(c_mat0)
      }
    } else if (i == samp_info$N_obs) {
      d_vec <- mean_mat[, i] + M %*% y[, i - 1]
      pre_calcs <- final_mats[[iter_samp_type]]
    } else {
      d_vec <- c_mat_def %*% (sig_inv %*% (mean_mat[, i] + M %*% y[, i - 1]) +
                                ms_cross %*% (y[, i + 1] - mean_mat[, i + 1]))
      pre_calcs <- default_mats[[iter_samp_type]]
    }

    ## Store results
    if (fast) {
      y[iter_locs, i] <- tmvn_gibbs_rej_fast(
        y_current = y[, i], mean = d_vec, lower = cuts_low,
        upper = cuts_high, locs = iter_locs, loc_length = iter_length,
        pre_calcs = pre_calcs, max_iter = samp_info$max_iter, N_ord = num_ord,
        burn_in = samp_info$burn_in, num_iter = num_iter, num_obs = i)
    } else {
      y[iter_locs, i] <- tmvn_gibbs_rej(
        y_current = y[, i], mean = d_vec, lower = cuts_low,
        upper = cuts_high, locs = iter_locs, loc_length = iter_length,
        pre_calcs = pre_calcs, max_iter = samp_info$max_iter, N_ord = num_ord,
        burn_in = samp_info$burn_in, num_iter = num_iter, num_obs = i)
    }
  }
  t(y)
}

#' Gibbs step for truncated multivariate normal
#' @param y_current continous outcome values
#' @param mean mean vector of multivariate normal
#' @param sigma covariance matrix
#' @param lower vector of lower thresholds
#' @param upper vector of upper thresholds
#' @param locs locations to sample
#' @param pre_cals pre calculations, including the conditional variance
#' @return vector
#' @importFrom truncnorm rtruncnorm
#' @importFrom MASS mvrnorm

tmvn_gibbs_rej <- function(y_current, mean, lower, upper, locs, loc_length,
                           pre_calcs, max_iter, N_ord, burn_in, num_iter, num_obs) {

  ## Limit max iter the first 100 iterations
  if(num_iter <= 100) max_iter <- 100

  ## If missing a single ordinal outcome then sample from truncated normal
  if(loc_length == 1) {
    res <- rtruncnorm(
      n = 1, a = lower[1], b = upper[1], sd = pre_calcs$cond_cov,
      mean = cond_mean_part(y_current, mean, pre_calcs$mean_pre, 1))
    rej_samp_iters <<- 0
    ## Otherwise use a naive rejection sampling approach
  } else {
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
      if(num_iter > burn_in) {
        print(paste0("Num_obs to break is ", num_obs, "; iteration, ", num_iter, "\\n"))
        stop("Sampling from the truncated multivariate normal was inefficient and could not be completed. Please increase max_iter.")
      } else if (num_iter > 100) {
        print(paste0("Num_obs to break is ", num_obs, "; iteration, ", num_iter, "\\n"))
        warning("Sampling from the truncated multivariate normal could not be completed during burn-in.")
      }
      res <- y_tmp
      res[1:N_ord] <- ifelse(is.infinite(upper[1:N_ord]), lower[1:N_ord], upper[1:N_ord])
    }
    rej_samp_iters <<- iter
  }
  res
}

#' Gibbs step for truncated multivariate normal
#' @param y_current continous outcome values
#' @param mean mean vector of multivariate normal
#' @param sigma covariance matrix
#' @param lower vector of lower thresholds
#' @param upper vector of upper thresholds
#' @param locs locations to sample
#' @param pre_cals pre calculations, including the conditional variance
#' @return vector
#' @importFrom truncnorm rtruncnorm
#' @importFrom MASS mvrnorm

tmvn_gibbs_rej_fast <- function(y_current, mean, lower, upper, locs, loc_length,
                           pre_calcs, max_iter, N_ord, burn_in, num_iter, num_obs) {

  ## Limit max iter the first 100 iterations
  if(num_iter <= 100) max_iter <- 100

  ## If missing a single ordinal outcome then sample from truncated normal
  if(loc_length == 1) {
    res <- rtruncnorm(
      n = 1, a = lower[1], b = upper[1], sd = pre_calcs$cond_cov,
      mean = cond_mean_part(y_current, mean, pre_calcs$mean_pre, 1))
    rej_samp_iters <<- 0
    ## Otherwise use a naive rejection sampling approach
  } else {
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
      mean_vec <- as.vector(cond_mean_part(y_current, mean, pre_calcs$mean_pre, locs))
      res <- rtmvnorm(
        1, mean = mean_vec, H = pre_calcs$cond_chol_inv, lower = lower,
        upper = upper, algorithm = "gibbs", burn.in.samples = 10, start.value = start_val)
    }
  }
  res
}
