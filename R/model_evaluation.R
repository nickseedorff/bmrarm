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

get_DIC_ar <- function(samps, marginal = FALSE) {

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

    dev_vals <- get_dev_ar(
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

  dev_of_means <- get_dev_ar(
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

get_dev_ar <- function(y, X, z, Z_kron, pat_idx, pat_idx_long, beta, sig_list,
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


#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @return scalar

get_pred <- function(samps, new_X = NULL) {

  N_outcomes <- samps$samp_info$N_outcomes
  z <- samps$z
  miss_mat <- samps$samp_info$miss_mat
  pat_idx <- samps$samp_info$pat_idx_long
  pat_idx_long <- samps$samp_info$pat_idx_long
  if(is.null(new_X)) {
    X <- samps$X
  } else {
    X <- new_X
  }
  Z_kron <- samps$Z_kron

  ## Calculate mean deviance
  mean_dev <- lapply(1:ncol(samps$res_beta), function(x){
    beta <- matrix(samps$res_beta[, x], ncol = N_outcomes)
    sig_tmp <- matrix(samps$res_sigma[, x], ncol = N_outcomes)
    pat_eff <- samps$res_pat_eff[,, x]
    y_use <- samps$res_y[,, x]
    cuts <- samps$res_cuts[, x]
    ar_tmp <- samps$res_ar[x]
    cur_draws <- list(ar = ar_tmp, sigma = sig_tmp)
    sig_list <- get_sig_list(cur_draws, samps$samp_info)

    for(j in 1:length(sig_list$marg_cov_list)) {
      sig_list$marg_cov_inv_list[[j]] <- chol2inv(chol(sig_list$marg_cov_list[[j]]))
    }

    mean_vec <- as.numeric(X %*% beta) + rowSums(Z_kron * pat_eff[pat_idx_long, ])
    mean_mat <- matrix(mean_vec, ncol = ncol(y_use))
    y_mat <- y_use
    y_mat[] <- NA
    y_ord <- as.vector(y_mat[, 1])

    ## storage for likelihood evaluations
    N_pat <- length(unique(pat_idx))
    y_last <- matrix(NA, nrow = N_pat, 2)
    for(i in 1:N_pat) {
      ## Locations of vectors to sample
      locs <- samps$samp_info$pat_locs[[i]]
      times <- samps$samp_info$pat_times[[i]]
      dist_mat <- ar_tmp ^ as.matrix(dist(times, diag = T, upper = T))
      chol_mat <- kronecker(chol(sig_tmp), chol(dist_mat))
      y_mat[locs,] <- LaplacesDemon::rmvnc(1, as.vector(mean_mat[locs, ]), chol_mat)
      #y_mat[locs,] <- as.vector(mean_mat[locs, ]) + chol_mat %*% rnorm(length(locs) * 2)

      ## Predict final locations
      locs_use <- c(length(locs), 2 * length(locs))
      locs_other <- setdiff(1:(2 * length(locs)), locs_use)
      full_sig <- kronecker(sig_tmp, dist_mat)
      pre_calcs <- pre_calc_ar(full_sig, locs_use)
      mean_val <- mean_mat[locs, ][locs_use] +
        pre_calcs$mean_pre %*% (y_use[locs, ][-locs_use] - mean_mat[locs, ][-locs_use])
      y_last[i, ] <- MASS::mvrnorm(1, mean_val, Sigma = pre_calcs$cond_cov)
    }

    ## Discretize
    for(k in 1:length(y_ord)) {
      for(l in 1:(length(cuts) - 1)) {
        if(y_mat[k, 1] > cuts[l] & y_mat[k, 1] <= cuts[l + 1]) {
          y_ord[k] <- l
        }
      }
    }

    y_last_tmp <- y_last
    for(k in 1:nrow(y_last)) {
      for(l in 1:(length(cuts) - 1)) {
        if(y_last_tmp[k, 1] > cuts[l] & y_last_tmp[k, 1] <= cuts[l + 1]) {
          y_last[k, 1] <- l
        }
      }
    }

    print(x)
    list(all_preds = cbind(y_ord, y_mat[, 2]),
         last_preds = y_last)
  })
}

#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @return scalar

get_forecasts <- function(samps) {

  N_outcomes <- samps$samp_info$N_outcomes
  z <- samps$z
  miss_mat <- samps$samp_info$miss_mat
  pat_idx <- samps$samp_info$pat_idx_long
  pat_idx_long <- samps$samp_info$pat_idx_long
  X <- samps$X
  Z_kron <- samps$Z_kron
  X_tmp <- as.matrix(samps$X_for)
  Z_kron_for <- samps$Z_kron_for
  pat_idx_for <- samps$pat_long_for
  pat_idx_long_for <- rep(samps$pat_long_for, 2)

  ## Calculate mean deviance
  mean_dev <- lapply(1:ncol(samps$res_beta), function(x){
    beta <- matrix(samps$res_beta[, x], ncol = N_outcomes)
    sig_tmp <- matrix(samps$res_sigma[, x], ncol = N_outcomes)
    pat_eff <- samps$res_pat_eff[,, x]
    y_use <- samps$res_y[,, x]
    cuts <- samps$res_cuts[, x]
    ar_tmp <- samps$res_ar[x]

    mean_vec_obs <- as.numeric(X %*% beta) + rowSums(Z_kron * pat_eff[pat_idx_long, ])
    mean_mat_obs <- matrix(mean_vec_obs, ncol = 2)

    mean_vec <- as.numeric(X_tmp %*% beta) + rowSums(Z_kron_for * pat_eff[pat_idx_long_for, ])
    mean_mat <- matrix(mean_vec, ncol = 2)

    ## storage for likelihood evaluations
    N_pat <- length(unique(pat_idx))
    y_last <- array(NA, c(N_pat, 2, 4))
    for(i in 1:N_pat) {
      ## Locations of vectors to sample
      times <- X_tmp[pat_idx_for == i, 3]
      locs <- samps$samp_info$pat_locs[[i]]
      locs_for <- which(pat_idx_for == i)
      length_locs <- length(locs)
      max_loc <- max(locs_for) - 4
      dist_mat <- ar_tmp ^ as.matrix(dist(times, diag = T, upper = T))

      for(j in 1:4) {
        full_sig <- kronecker(sig_tmp, dist_mat[c(1:length_locs, length_locs + j),
                                                c(1:length_locs, length_locs + j)])
        pre_calcs <- pre_calc_ar(full_sig, c(length_locs + 1, length_locs * 2 + 2))
        mean_val <- mean_mat[max_loc + j, ] +
          pre_calcs$mean_pre %*% (as.vector(y_use[locs, ]) - as.vector(mean_mat_obs[locs, ]))
        y_last[i,,j] <- MASS::mvrnorm(1, mean_val, Sigma = pre_calcs$cond_cov)

        ## Discretize
        tmp_val <- y_last[i,1,j]
        for(l in 1:(length(cuts) - 1)) {
          if(tmp_val > cuts[l] & tmp_val <= cuts[l + 1]) {
            y_last[i,1,j]  <- l
          }
        }

      }
    }

    print(x)
    list(y_last)
  })
}
