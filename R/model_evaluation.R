#' Get DIC
#'
#' @param y_current current continuous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @return scalar

get_DIC <- function(samps) {

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

    dev_vals <- get_dev(
      y = y_use, X, z, Z_kron, pat_idx, pat_idx_long, beta = beta_tmp,
      sigma = sig_tmp, cuts = cuts_tmp, pat_eff, miss_mat)
  }) %>%
    unlist() %>%
    mean()

  ## Calculate deviance of the mean
  beta_tmp <- matrix(rowMeans(samps$res_beta), ncol = N_outcomes)
  sig_tmp <- matrix(rowMeans(samps$res_sigma), ncol = N_outcomes)
  pat_eff <- apply(samps$res_pat_eff, c(1, 2), mean)
  y_use <- apply(samps$res_y, c(1, 2), mean)
  cuts_tmp <- rowMeans(samps$res_cuts)

  dev_of_means <- get_dev(
    y = y_use, X, z, Z_kron, pat_idx, pat_idx_long, beta = beta_tmp,
    sigma = sig_tmp, cuts = cuts_tmp, pat_eff, miss_mat)

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
#' @return scalar

get_waic_ar <- function(samps) {

  N_outcomes <- samps$samp_info$N_outcomes
  z <- samps$z
  miss_mat <- samps$samp_info$miss_mat
  pat_idx <- samps$samp_info$pat_idx_long
  pat_idx_long <- samps$samp_info$pat_idx_long
  X <- samps$X
  Z_kron <- samps$Z_kron

  ## Calculate mean deviance
  mean_dev <- sapply(1:ncol(samps$res_beta), function(x){
    beta <- matrix(samps$res_beta[, x], ncol = N_outcomes)
    sigma <- matrix(samps$res_sigma[, x], ncol = N_outcomes)
    pat_eff <- samps$res_pat_eff[,, x]
    y_use <- samps$res_y[,, x]
    cuts<- samps$res_cuts[, x]
    #ar_tmp <- samps$res_ar[x]
    ar_tmp <- 0
    cur_draws <- list(ar = ar_tmp, sigma = sigma)
    sig_list <- get_sig_list(cur_draws, samps$samp_info)

    for(j in 1:length(sig_list$marg_cov_list)) {
      sig_list$marg_cov_inv_list[[j]] <- chol2inv(chol(sig_list$marg_cov_list[[j]]))
    }

    get_dev(
      y = y_use, X, z, Z_kron, pat_idx, pat_idx_long, beta,
      sigma, cuts, pat_eff, miss_mat, all_prob = T)
  })
  #-2 * (sum(log(rowMeans(mean_dev_obs))) - sum(apply(log(mean_dev_obs), 1, var)))
  #loo(t(log(mean_dev_obs)), is_method = "tis")
  print(loo(t(log(mean_dev))))

  ## Second piece
  sum(rowMeans(log(mean_dev))) * 2
  sum(log(rowMeans(mean_dev))) * 2

  pw <- sum(log(rowMeans(mean_dev)) - rowMeans(log(mean_dev))) * 2
  pw2 <- sum(apply(log(mean_dev), 1, var))
  lppd <- sum(log(rowMeans(mean_dev)))
  (waic <- -2 * lppd + 2 * pw)
  (waic2 <- -2 * lppd + 2 * pw2)

  list(estim =   c(waic2 = waic2, waic = waic, pw = pw, pw2 = pw2,
                   dev1 = sum(rowMeans(log(mean_dev))) * 2,
                   dev2 = sum(log(rowMeans(mean_dev))) * 2),
       l_row_mean_dev = log(rowMeans(mean_dev)))
}

#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @return scalar

get_dev <- function(y, X, z, Z_kron, pat_idx, pat_idx_long, beta, sigma, cuts,
                    pat_eff, miss_mat, all_prob = FALSE) {

  mean_vec <- as.numeric(X %*% beta) + rowSums(Z_kron * pat_eff[pat_idx_long, ])
  mean_mat <- matrix(mean_vec, ncol = ncol(y))

  ## Conditional means and covariances
  pre_calcs <- pre_calc_ar(sigma, 1)
  cond_sd <- sqrt(as.numeric(pre_calcs$cond_cov))
  cond_mean <- mean_mat[, 1] + as.numeric(pre_calcs$mean_pre) *
    (as.numeric(y[, -1]) - as.numeric(mean_mat[, -1]))

  ## Calculate cuts
  cuts_low <- cuts[z]
  cuts_low[is.na(cuts_low)] <- -Inf
  cuts_high <- cuts[z + 1]
  cuts_high[is.na(cuts_high)] <- Inf

  ## Integrate over the latent outcomes
  probs <- log(pnorm((cuts_high - cond_mean) / cond_sd) -
                 pnorm((cuts_low - cond_mean) / cond_sd))

  cont_vals <- dnorm(y[, 2], mean = mean_mat[, 2],
                     sd = sqrt(sigma[2, 2]), log = T)

  ## Calculate deviance
  dev_val <- -2 * (sum(probs) + sum(cont_vals[!miss_mat[, 2]]))
  dev_ord <- -2 * sum(probs)
  dev_cont <- -2 * sum(cont_vals[!miss_mat[, 2]])
  dev_val = dev_val
  #tail(cbind(probs, cont_vals, y, mean_mat, sqrt(sigma[2, 2]), cond_sd, cond_mean, cuts_low, cuts_high), 7)
  #list(dev_val = dev_val, dev_ord = dev_ord, dev_cont = dev_cont)
  if(all_prob) {
    exp_cont_vals <- exp(cont_vals)
    exp_cont_vals[miss_mat[, 2]] <- 1
    exp_probs <- exp(probs)
    exp_cont_vals * exp_probs
  } else {
    dev_val
  }
}

#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @return scalar

get_DIC_ar <- function(samps) {

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
    cur_draws <- list(ar = ar_tmp, sigma = sig_tmp)
    sig_list <- get_sig_list(cur_draws, samps$samp_info)

    for(j in 1:length(sig_list$marg_cov_list)) {
      sig_list$marg_cov_inv_list[[j]] <- chol2inv(chol(sig_list$marg_cov_list[[j]]))
    }

    dev_vals <- get_dev_ar(
      y = y_use, X, z, Z_kron, pat_idx, pat_idx_long,
      beta = beta_tmp, sig_list, cuts = cuts_tmp, pat_eff = pat_eff,
      miss_mat = miss_mat, samps$samp_info)
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
  sig_list <- get_sig_list(cur_draws, samps$samp_info)

  for(j in 1:length(sig_list$marg_cov_list)) {
    sig_list$marg_cov_inv_list[[j]] <- chol2inv(chol(sig_list$marg_cov_list[[j]]))
  }

  dev_of_means <- get_dev_ar(
    y = y_use, X, z, Z_kron, pat_idx, pat_idx_long, beta = beta_tmp,
    sig_list, cuts = cuts_tmp, pat_eff, miss_mat, samps$samp_info)

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

get_dev_ar <- function(y, X, z, Z_kron, pat_idx, pat_idx_long, beta, sig_list,
                       cuts, pat_eff, miss_mat, samp_info) {

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
    #probs_check1[i] <- sum(probs2[locs])
    #probs_check2[i] <- sum(cont_vals2[locs])
  }

  ## Calculate deviance
  dev_val <- -2 * (sum(probs) + sum(cont_vals))
  dev_ord <- -2 * sum(probs)
  dev_cont <- -2 * sum(cont_vals)
  dev_val = dev_val
  list(dev_val = dev_val, dev_ord = dev_ord, dev_cont = dev_cont)
  dev_val
}




#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @importFrom mvtnorm dmvnorm
#' @return scalar

get_dev_ar2 <- function(y, X, z, Z_kron, pat_idx, pat_idx_long, beta, sig_list,
                       cuts, pat_eff, miss_mat, samp_info) {

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
    locs <- samps$samp_info$pat_locs[[i]]
    miss_locs <- samps$samp_info$pat_cont_miss_rank[[i]]

    ## Conditional mean
    ind <- samps$samp_info$pat_time_ind[i]
    tmp_mean <- as.numeric(
      mean_mat[locs, 1] + sig_list$mean_pre_list[[ind]] %*%
        (as.numeric(y_mat[locs, -1]) - as.numeric(mean_mat[locs, -1])))

    ## Evaluate the multivariate normal probabilities
    probs[i] <- omxMnor(
      lbound = cuts_low[locs], ubound = cuts_high[locs],
      mean = tmp_mean, covariance = sig_list$cond_cov_list[[ind]])[[1]]

    ## Evaulate multivariate continuous outcomes
    if(length(miss_locs) != 0) {
      y_cont <- as.numeric(y_mat[locs, -1])
      pat_cont_mean <- as.numeric(mean_mat[locs, -1])
      cont_use <- setdiff(1:length(locs), miss_locs - length(locs))
      pre_calcs <- pre_calc_ar(sig_list$marg_cov_list[[ind]], cont_use)
      cond_mean <- pat_cont_mean[cont_use] + pre_calcs$mean_pre %*%
        (y_cont[-cont_use] - pat_cont_mean[- cont_use])
      cont_vals[i] <- dmvnorm(x = y_cont[cont_use], mean = cond_mean,
                              sigma = pre_calcs$cond_cov)
    } else {
      y_cont <- as.numeric(y_mat[locs, -1])
      pat_cont_mean <- as.numeric(mean_mat[locs, -1])
      cont_vals[i] <- LaplacesDemon::dmvnp(x = y_cont, mu = pat_cont_mean,
                                           Omega = sig_list$marg_cov_inv_list[[ind]])
    }
    #probs_check1[i] <- sum(probs2[locs])
    #probs_check2[i] <- sum(cont_vals2[locs])
  }

  ## Calculate deviance
  #probs + cont_vals
  probs
}

#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @importFrom mvtnorm dmvnorm
#' @return scalar

get_dev_ar3 <- function(y, X, z, Z_kron, pat_idx, pat_idx_long, beta, sig_list,
                        cuts, pat_eff, miss_mat, samp_info) {

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
    ind <- samps$samp_info$pat_time_ind[i]
    ## Locations of vectors to sample
    locs <- samps$samp_info$pat_locs[[i]]
    probs[i] <-
      LaplacesDemon::dmvnp(x = as.numeric(y_mat[locs, ]),
                           mu = as.numeric(mean_mat[locs, ]),
                           Omega = sig_list$sig_inv_list[[ind]])

  }

  ## Calculate deviance
  #probs + cont_vals
  probs
}


#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @importFrom mvtnorm dmvnorm
#' @return scalar

get_dev_ar4 <- function(y, X, z, Z_kron, pat_idx, pat_idx_long, beta, sig_list,
                        cuts, pat_eff, miss_mat, samp_info) {

  mean_vec <- as.numeric(X %*% beta) + rowSums(Z_kron * pat_eff[pat_idx_long, ])
  mean_mat <- matrix(mean_vec, ncol = ncol(y))
  y_mat <- matrix(y, ncol = 2)

  ## Calculate cuts
  cuts_low <- cuts[c(z, rep(NA, length(z)))]
  cuts_low[is.na(cuts_low)] <- -Inf
  cuts_high <- cuts[c(z + 1, rep(NA, length(z)))]
  cuts_high[is.na(cuts_high)] <- Inf

  ## storage for likelihood evaluations
  N_pat <- length(unique(pat_idx))
  probs_o <- probs_c <- rep(NA, samp_info$N_obs)
  iter_o <- iter_c <- 1
  for(i in 1:N_pat) {
    ## Locations of vectors to sample
    ind <- samp_info$pat_time_ind[i]
    locs <- c(samp_info$pat_all_locs[[i]])

    ## Meam and covariance
    sig_inv <- sig_list$sig_inv_list[[ind]]
    g <- sig_inv %*% (y[locs] - mean_vec[locs])

    for(j in 1:length(locs)) {
      var_use <- 1 / sig_inv[j, j]
      sd_use <- sqrt(var_use)
      tmp_mean <- y[locs][j] - g[j] * var_use
      if(j <= length(locs) / 2) {
        probs_o[iter_o] <- pnorm((cuts_high[locs][j] - tmp_mean) / sd_use) -
          pnorm((cuts_low[locs][j] - tmp_mean) / sd_use)
        iter_o <- iter_o + 1
      } else {
        probs_c[iter_c] <- dnorm(x = y[locs][j], mean = tmp_mean, sd = sd_use)
        iter_c <- iter_c + 1
      }
    }
  }

  ## Calculate deviance
  #probs + cont_vals
  c(probs_o, probs_c)
}
