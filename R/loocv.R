#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @importFrom mvtnorm dmvnorm
#' @return scalar

get_loocv <- function(samps, cv_locs, z_true, y_true, full_X, full_Z_kron,
                      data, ar = FALSE) {

  N_outcomes <- samps$samp_info$N_outcomes
  if(!ar) {
    z <- z_true[cv_locs]
    y <- y_true[cv_locs]
    miss_mat <- is.na(cbind(z_true, y_true))[cv_locs, ]
    pat_idx <- data$pat_idx[cv_locs]
    pat_idx_long <- c(pat_idx, pat_idx)
    X <- full_X[cv_locs, ]
    Z_kron <- full_Z_kron[c(cv_locs, cv_locs + nrow(data)), ]
  } else {
    z <- z_true
    miss_mat <- is.na(cbind(z_true, y_true))
    pat_idx <- data$pat_idx
    pat_idx_long <- c(pat_idx, pat_idx)
    X <- full_X
    Z_kron <- full_Z_kron
    y_full <- matrix(NA, ncol = N_outcomes, nrow = nrow(data))
    y_full[cv_locs, 2] <- y_true[cv_locs]
  }

  ## Calculate mean deviance
  all_lppd <- sapply(1:ncol(samps$res_beta), function(x){
    beta <- matrix(samps$res_beta[, x], ncol = N_outcomes)
    pat_eff <- samps$res_pat_eff[,, x]
    cuts <- samps$res_cuts[, x]
    sigma <- matrix(samps$res_sigma[, x], ncol = N_outcomes)

    if(!ar) {
      get_lppd(
        y, X, z, Z_kron, pat_idx, pat_idx_long,
        beta, sigma, cuts, pat_eff = pat_eff, miss_mat, cv_locs)
    } else {
      y_full[-cv_locs,] <- samps$res_y[,, x]
      ar_val <- samps$res_ar[x]

      get_lppd_ar(
        y = y_full, X, z, Z_kron, pat_idx, pat_idx_long,
        beta, sigma, cuts, pat_eff, miss_mat,
        samp_info = samps$samp_info, cv_locs, ar_val, data)
    }
  })
  log(rowMeans(all_lppd))
}

#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @importFrom mvtnorm dmvnorm
#' @return scalar

get_lppd <- function(y, X, z, Z_kron, pat_idx, pat_idx_long, beta, sigma,
                        cuts, pat_eff, miss_mat, samp_info, cv_locs) {

  mean_vec <- as.numeric(X %*% beta) + rowSums(Z_kron * pat_eff[pat_idx_long, ])
  mean_mat <- matrix(mean_vec, ncol = 2)

  ## Calculate cuts
  cuts_low <- cuts[z]
  cuts_low[is.na(cuts_low)] <- -Inf
  cuts_high <- cuts[c(z + 1)]
  cuts_high[is.na(cuts_high)] <- Inf

  ## Conditional means and covariances
  pre_calcs <- pre_calc_ar(sigma, 1)
  cond_sd <- sqrt(as.numeric(pre_calcs$cond_cov))
  cond_mean <- mean_mat[, 1] + as.numeric(pre_calcs$mean_pre) *
    (y - as.numeric(mean_mat[, -1]))

  ## Probabilities when continuous y is observed
  probs <- pnorm((cuts_high - cond_mean) / cond_sd) -
    pnorm((cuts_low - cond_mean) / cond_sd)

  ## Probabilities when continuous y is missing
  probs2 <- pnorm((cuts_high - mean_mat[, 1]) / sqrt(sigma[1, 1])) -
    pnorm((cuts_low - mean_mat[, 1]) / sqrt(sigma[1, 1]))

  ## Probabilities to use for the latent outcomes
  probs <- coalesce(probs, probs2)
  cont_vals <- dnorm(y, mean = mean_mat[, 2], sd = sqrt(sigma[2, 2]))
  cont_vals[is.na(cont_vals)] <- 1
  probs * cont_vals
}


#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @importFrom mvtnorm dmvnorm
#' @return scalar

get_lppd_ar <- function(y, X, z, Z_kron, pat_idx, pat_idx_long, beta, sigma,
                        cuts, pat_eff, miss_mat, samp_info, cv_locs, ar_val,
                        data) {

  mean_vec <- as.numeric(X %*% beta) + rowSums(Z_kron * pat_eff[pat_idx_long, ])
  mean_mat <- matrix(mean_vec, ncol = ncol(y))

  ## Calculate cuts
  cuts_low <- cuts[z]
  cuts_low[is.na(cuts_low)] <- -Inf
  cuts_high <- cuts[z + 1]
  cuts_high[is.na(cuts_high)] <- Inf

  ## storage for likelihood evaluations
  probs <- matrix(NA, nrow = length(cv_locs), ncol = 2)
  res_idx <- 1
  sig_small_inv <- chol2inv(chol(sigma))
  time_old <- 0
  for(i in cv_locs) {
    ## Locations of vectors to sample
    pat <- pat_idx[i]
    ord_locs <- which(pat_idx == pat)
    all_locs <- which(pat_idx_long == pat)
    times <- data$time[ord_locs]
    ord_idx <- which(ord_locs == i)
    cont_idx <- ord_idx + length(ord_locs)

    ## Full covariance matrix
    dist_mat <- as.matrix(dist(times, diag = T, upper = T))
    sig_kron <- kronecker(sigma, ar_val ^ dist_mat)

# Likelihood contribution from continuous outcome ------------------------

    ## Covariance and mean
    sig_inv <- chol2inv(chol(sig_kron[-ord_idx, -ord_idx]))
    g <- sig_inv %*% (y[ord_locs, ][-ord_idx] - mean_mat[ord_locs, ][-ord_idx])
    var_use <- 1 / sig_inv[cont_idx - 1, cont_idx - 1]
    sd_use <- sqrt(var_use)
    tmp_mean <- y[ord_locs, ][cont_idx] - g[cont_idx - 1] * var_use

    if(!miss_mat[i, 2]) {
      probs[res_idx, 2] <- dnorm(x = y[i, 2], mean = tmp_mean, sd = sd_use)
    } else {
      probs[res_idx, 2] <- 1
    }

# Likelihood contribution from Ordinal outcome ------------------------

    ## Covariance and mean
    if(miss_mat[i, 2]) {
      pre_calcs <- pre_calc_ar(sig_kron[-cont_idx, -cont_idx], locs = ord_idx)
      cond_sd <- sqrt(as.numeric(pre_calcs$cond_cov))
      cond_mean <- mean_mat[ord_locs, ][ord_idx] +
        as.numeric(pre_calcs$mean_pre) %*%
        (as.numeric(y[setdiff(ord_locs, i), ]) -
           as.numeric(mean_mat[setdiff(ord_locs, i), ]))
    } else {
      #pre_calcs <- pre_calc_ar(sig_kron, locs = ord_idx)
      #cond_sd <- sqrt(as.numeric(pre_calcs$cond_cov))
      #cond_mean <- mean_mat[ord_locs, ][ord_idx] +
      #  as.numeric(pre_calcs$mean_pre) %*% (y[ord_locs, ][-ord_idx] -
      #                                        mean_mat[ord_locs, ][-ord_idx])
      dist_inv <- chol2inv(chol(ar_val ^ dist_mat))
      sig_inv_ord <- kronecker(sig_small_inv, dist_inv)
      y[ord_locs, ][ord_idx] <- 1000
      g <- sig_inv_ord %*% (as.numeric(y[ord_locs, ]) - as.numeric(mean_mat[ord_locs, ]))
      var_use <- 1 / sig_inv_ord[ord_idx, ord_idx]
      cond_sd <- sqrt(var_use)
      cond_mean <- y[ord_locs, ][ord_idx] - g[ord_idx] * var_use
    }

    ## Evaluate ordinal outcome
    if(!is.na(z[i])) {
      probs[res_idx, 1] <- pnorm((cuts_high[i] - cond_mean) / cond_sd) -
        pnorm((cuts_low[i] - cond_mean) / cond_sd)
    } else {
      probs[res_idx, 1] <- 1
    }

    ## Update indexing
    res_idx <- res_idx + 1
    time_old <- times
  }

  probs[, 1] * probs[, 2]
}
