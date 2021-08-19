#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @return vector
#' @export

get_DIC <- function(samps, marginal = TRUE) {

  N_outcomes <- samps$samp_info$N_outcomes
  z <- samps$z
  miss_mat <- samps$samp_info$miss_mat
  pat_idx <- samps$samp_info$pat_idx_long
  pat_idx_long <- samps$samp_info$pat_idx_long
  X <- samps$X
  Z_kron <- samps$Z_kron

  ## Calculate mean deviance
  mean_dev <- lapply(1:ncol(samps$res_beta), function(x){
    beta_tmp <- matrix(samps$res_beta[, x], ncol = N_outcomes)
    sig_tmp <- matrix(samps$res_sigma[, x], ncol = N_outcomes)
    pat_eff <- samps$res_pat_eff[,, x]
    y_use <- samps$res_y[,, x]
    cuts_tmp <- samps$res_cuts[, x]
    ar_tmp <- samps$res_ar[x]
    if(!marginal) {
      cur_draws <- list(ar = ar_tmp, sigma = sig_tmp)
      sig_list <- get_sig_list(cur_draws, samps$samp_info)
    } else {
      pat_sig <- matrix(samps$res_pat_sig[, x], ncol = ncol(pat_eff))
      cur_draws <- list(ar = ar_tmp, sigma = sig_tmp, pat_sig = pat_sig)
      sig_list <- get_sig_list_marg(cur_draws, samps$samp_info)
      pat_eff[] <- 0
    }

    for(j in 1:length(sig_list$marg_cov_list)) {
      sig_list$marg_cov_inv_list[[j]] <- chol2inv(chol(sig_list$marg_cov_list[[j]]))
    }

    dev_vals <- get_dev(
      y = y_use, X, z, Z_kron, pat_idx, pat_idx_long,
      beta = beta_tmp, sig_list, cuts = cuts_tmp, pat_eff = pat_eff,
      miss_mat = miss_mat, samps$samp_info, marginal)
  }) %>%
    unlist() %>%
    mean()

  ## Calculate deviance of the mean
  beta_tmp <- matrix(rowMeans(samps$res_beta), ncol = N_outcomes)
  sig_tmp <- matrix(rowMeans(samps$res_sigma), ncol = N_outcomes)
  pat_eff <- apply(samps$res_pat_eff, c(1, 2), mean)
  y_use <- apply(samps$res_y, c(1, 2), mean)
  cuts_tmp <- rowMeans(samps$res_cuts)
  ar_tmp <- mean(samps$res_ar)
  cur_draws <- list(ar = ar_tmp, sigma = sig_tmp)
  if(!marginal) {
    cur_draws <- list(ar = ar_tmp, sigma = sig_tmp)
    sig_list <- get_sig_list(cur_draws, samps$samp_info)
  } else {
    pat_sig <- matrix(rowMeans(samps$res_pat_sig), ncol = ncol(pat_eff))
    cur_draws <- list(ar = ar_tmp, sigma = sig_tmp, pat_sig = pat_sig)
    sig_list <- get_sig_list_marg(cur_draws, samps$samp_info)
    pat_eff[] <- 0
  }

  for(j in 1:length(sig_list$marg_cov_list)) {
    sig_list$marg_cov_inv_list[[j]] <- chol2inv(chol(sig_list$marg_cov_list[[j]]))
  }

  dev_of_means <- get_dev(
    y = y_use, X, z, Z_kron, pat_idx, pat_idx_long, beta = beta_tmp,
    sig_list, cuts = cuts_tmp, pat_eff, miss_mat, samps$samp_info, marginal)

  ## Output DIC, D, pD
  pd <- mean_dev - dev_of_means
  c(DIC = pd + mean_dev, D = mean_dev, pd = pd, dev_of_means = dev_of_means)
}

#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @importFrom mvtnorm dmvnorm
#' @return scalar

get_dev <- function(y, X, z, Z_kron, pat_idx, pat_idx_long, beta, sig_list,
                       cuts, pat_eff, miss_mat, samp_info, marginal) {

  mean_vec <- as.numeric(X %*% beta) + rowSums(Z_kron * pat_eff[pat_idx_long, ])
  mean_mat <- matrix(mean_vec, ncol = ncol(y))
  y_mat <- matrix(y, ncol = 2)

  ## Calculate cuts
  cuts_low <- cuts[z]
  cuts_low[is.na(cuts_low)] <- -Inf
  cuts_high <- cuts[z + 1]
  cuts_high[is.na(cuts_high)] <- Inf

  ## storage for likelihood evaluations
  N_pat <- length(unique(pat_idx))
  probs <- probs_check1 <- probs_check2 <- cont_vals <- rep(NA, N_pat)
  for(i in 1:N_pat) {
    ## Locations of vectors to sample
    locs <- samp_info$pat_locs[[i]]
    miss_locs <- samp_info$pat_cont_miss_rank[[i]]


    ## Conditional mean
    ind <- samp_info$pat_time_ind[i]
    #if(marginal) ind <- i
    tmp_mean <- as.numeric(
      mean_mat[locs, 1] + sig_list$mean_pre_list[[ind]] %*%
        (as.numeric(y_mat[locs, -1]) - as.numeric(mean_mat[locs, -1])))

    ## Evaluate the multivariate normal probabilities
    probs[i] <- log(omxMnor(
      lbound = cuts_low[locs], ubound = cuts_high[locs],
      mean = tmp_mean, covariance = sig_list$cond_cov_list[[ind]])[[1]])

    ## Evaulate multivariate continuous outcomes
    if(length(miss_locs) != 0) {
      y_cont <- as.numeric(y_mat[locs, -1])
      pat_cont_mean <- as.numeric(mean_mat[locs, -1])
      cont_use <- setdiff(1:length(locs), miss_locs - length(locs))
      pre_calcs <- pre_calc_ar(sig_list$marg_cov_list[[ind]], cont_use)
      cond_mean <- pat_cont_mean[cont_use] + pre_calcs$mean_pre %*%
        (y_cont[-cont_use] - pat_cont_mean[- cont_use])
      cont_vals[i] <- dmvnorm(x = y_cont[cont_use], mean = cond_mean,
                              sigma = pre_calcs$cond_cov, log = T)
    } else {
      y_cont <- as.numeric(y_mat[locs, -1])
      pat_cont_mean <- as.numeric(mean_mat[locs, -1])
      cont_vals[i] <- LaplacesDemon::dmvnp(x = y_cont, mu = pat_cont_mean,
                                           Omega = sig_list$marg_cov_inv_list[[ind]], log = T)
    }
  }

  ## Calculate deviance
  dev_val <- -2 * (sum(probs) + sum(cont_vals))
  dev_ord <- -2 * sum(probs)
  dev_cont <- -2 * sum(cont_vals)
  dev_val = dev_val
  list(dev_val = dev_val, dev_ord = dev_ord, dev_cont = dev_cont)
  dev_val
}
