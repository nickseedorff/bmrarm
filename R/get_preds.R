#' Generate iterative forecasts for bmrvarx object
#'
#' @param bmrvarx_obj a bmrvarx object
#' @param steps_ahead number of observations to forcast, done recursively
#' @param X matrix of covariates to make forecasts with, number of rows must
#' equal steps_ahead
#' @param get_ord_preds logical, discretize latent values into ordinal predictions
#' @importFrom MASS mvrnorm
#' @return list
#' @export

get_preds_bmrvarx <- function(bmrvarx_obj, steps_ahead = 5, X, seed = 15, get_ord_preds = T) {

  if(steps_ahead != nrow(X)) {
    stop("nrow(X) must equal steps_ahead")
  }

  set.seed(seed)
  ## Relevant data
  latent_draws <- t(as.matrix(bmrvarx_obj$data_for_forecasts[[1]][1:2, ]))
  y_cont <- bmrvarx_obj$data_for_forecasts[[1]][3, 1]
  covar_names <- bmrvarx_obj$covars_used

  ## Build design matrix
  if("(Intercept)" %in% covar_names) {
    X_use <- cbind(1, X)
  } else {
    X_use <- X
  }

  ## Number of simulations and outcomes
  N_samp <- nrow(latent_draws)
  N_outcomes <- length(y_cont) + ncol(latent_draws)
  N_cont <- length(y_cont)
  N_ord <- ncol(latent_draws)

  ## Relevant parameters
  M_draws <- bmrvarx_obj$draws$res_M
  sigma_draws <- bmrvarx_obj$draws$res_sigma
  beta_draws <- bmrvarx_obj$draws$res_beta
  cut_draws <- bmrvarx_obj$draws$res_cuts

  ## Store Predictions
  ord_preds <- cont_preds <- array(NA, dim = c(N_samp, N_outcomes, steps_ahead))

  for(i in 1:N_samp) {
    M <- matrix(M_draws[, i], ncol = N_outcomes)
    sigma <- matrix(sigma_draws[, i], ncol = N_outcomes)
    beta <- matrix(beta_draws[, i], ncol = N_outcomes)
    for(j in 1:steps_ahead) {
      if(j == 1) {
        mean_val <- crossprod(beta, X_use[j, ]) + M %*% c(latent_draws[i, ], y_cont)
      } else {
        mean_val <- crossprod(beta, X_use[j, ]) + M %*% cont_preds[i, , j - 1]
      }
      cont_preds[i, , j] <- mvrnorm(1, mean_val, Sigma = sigma)
    }

    if(get_ord_preds) {
      cuts <- matrix(cut_draws[, i, ], ncol = N_ord)
      for(j in 1:N_ord) {
        for(k in 1:steps_ahead) {
          ord_preds[i,j, k] <- cont_to_cat(cont_preds[i, j, k], cuts[, j])
        }
      }
    }
  }

  list(cont_preds = cont_preds,
       ord_preds = ord_preds[, 1:N_ord, ])
}

#' Discretize a latent value using cutpoints
#' @param y scalar continuous value
#' @param cut_points vector of cutpoints

cont_to_cat <- function(y, cut_points) {
  for(i in 1:(length(cut_points) - 1)) {
    if(y > cut_points[i] & y <= cut_points[i + 1]) {
      return(i)
    }
  }
}

#' Generate predictions for bmrarm object. Not intended for usage.
#'
#' @param samps bmrarm object with additional structure for new time points.
#' @return list

get_preds_bmrarm <- function(samps) {

  N_outcomes <- samps$data$samp_info$N_outcomes
  z <- samps$z
  miss_mat <- samps$data$samp_info$miss_mat
  pat_idx <- samps$data$samp_info$pat_idx_long
  pat_idx_long <- samps$data$samp_info$pat_idx_long
  X <- samps$data$X
  Z_kron <- samps$data$Z_kron
  X_tmp <- as.matrix(samps$data$X_for)
  X_which_time <- which(colnames(samps$data$X_for) == "time")
  Z_kron_for <- samps$data$Z_kron_for
  pat_idx_for <- samps$data$pat_long_for
  pat_idx_long_for <- rep(samps$data$pat_long_for, 2)

  ## Calculate mean deviance
  mean_dev <- lapply(1:ncol(samps$draws$res_beta), function(x){
    beta <- matrix(samps$draws$res_beta[, x], ncol = N_outcomes)
    sig_tmp <- matrix(samps$draws$res_sigma[, x], ncol = N_outcomes)
    pat_eff <- samps$draws$res_pat_eff[,, x]
    y_use <- samps$draws$res_y[,, x]
    cuts <- samps$draws$res_cuts[, x]
    ar_tmp <- samps$draws$res_ar[x]

    mean_vec_obs <- as.numeric(X %*% beta) + rowSums(Z_kron * pat_eff[pat_idx_long, ])
    mean_mat_obs <- matrix(mean_vec_obs, ncol = 2)

    mean_vec <- as.numeric(X_tmp %*% beta) + rowSums(Z_kron_for * pat_eff[pat_idx_long_for, ])
    mean_mat <- matrix(mean_vec, ncol = 2)

    ## storage for likelihood evaluations
    N_pat <- length(unique(pat_idx))
    y_last <- array(NA, c(N_pat, 2, 4))
    for(i in 1:N_pat) {
      ## Locations of vectors to sample
      times <- X_tmp[pat_idx_for == i, X_which_time]
      locs <- samps$data$samp_info$pat_locs[[i]]
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
    list(y_last)
  })
}

#' Generate predictions for bmrarm object. Not intended for usage.
#'
#' @param samps bmrarm object with additional structure for new time points.
#' @importFrom verification rps
#' @return list

summary_preds_bmrarm <- function(preds, truth) {

  N_iter <- length(preds)
  N_preds <- dim(preds[[1]][[1]])[1] * dim(preds[[1]][[1]])[3]
  pred_ord <- pred_cont <- matrix(NA, ncol = N_iter, nrow = N_preds)

  for(i in 1:N_iter) {
    ## Matrix of predictions
    mat <- rbind(preds[[i]][[1]][,, 1],
                 preds[[i]][[1]][,, 2],
                 preds[[i]][[1]][,, 3],
                 preds[[i]][[1]][,, 4])

    ## Store preds
    pred_ord[, i] <- mat[, 1]
    pred_cont[, i] <- mat[, 2]
  }

  ## Prediction summaries
  res <- as.data.frame(matrix(NA, nrow = N_preds, ncol = 10))
  colnames(res) <- c(paste0("prob", 1:5), "mean", "lower95", "upper95",
                     "pat", "num_ahead")
  res$pat <- rep(1:48, 4)
  res$num_ahead <- rep(1:4, each = 48)

  ## Ordinal probabilities
  res[, 1:5] <- apply(pred_ord, 1, function(x) {
    c(mean(x == 1), mean(x == 2), mean(x == 3), mean(x == 4), mean(x == 5))
  }) %>% t()

  ## Posterior mean and CI for continuous outcome
  res[, 6:8] <- apply(pred_cont, 1, function(x) {
    c(mean(x), quantile(x, probs = c(0.025, 0.975)))
  }) %>% t()

  ## Compare to truth
  full <- arrange(res, pat, num_ahead) %>%
    mutate(pat_idx = truth$pat_idx, time = truth$time,
           y_ord = truth$y_ord, y_cont = truth$y2,
           sq_err = (y_cont - mean) ^ 2,
           in_ci = y_cont >= lower95 & y_cont <= upper95)

  ## Summary metrics by time pint

  metrics <- sapply(1:4, function(x) {
    tmp_full <- full[full$num_ahead == x, ]
    c(rmse = sqrt(mean(tmp_full$sq_err)), in_ci = mean(tmp_full$in_ci),
      rps = rps(tmp_full$y_ord, as.matrix(tmp_full[, 1:5]))$rps)
  }) %>% t()

  list(predictions = full, metrics = metrics)
}


#' Generate posterior predictive draws for bmrarm object. Not intended for public use.
#'
#' @param samps bmrarm object with additional structure for new time points
#' @param new_X new design matrix
#' @return list

get_pred <- function(samps, new_X = NULL) {

  N_outcomes <- samps$data$samp_info$N_outcomes
  z <- samps$data$z
  miss_mat <- samps$data$samp_info$miss_mat
  pat_idx <- samps$data$samp_info$pat_idx_long
  pat_idx_long <- samps$data$samp_info$pat_idx_long
  if(is.null(new_X)) {
    X <- samps$data$X
  } else {
    X <- new_X
  }
  Z_kron <- samps$data$Z_kron

  ## Calculate mean deviance
  mean_dev <- lapply(1:ncol(samps$draws$res_beta), function(x){
    beta <- matrix(samps$draws$res_beta[, x], ncol = N_outcomes)
    sig_tmp <- matrix(samps$draws$res_sigma[, x], ncol = N_outcomes)
    pat_eff <- samps$draws$res_pat_eff[,, x]
    y_use <- samps$draws$res_y[,, x]
    cuts <- samps$draws$res_cuts[, x]
    ar_tmp <- samps$draws$res_ar[x]
    cur_draws <- list(ar = ar_tmp, sigma = sig_tmp)
    sig_list <- get_sig_list(cur_draws, samps$data$samp_info)

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
      locs <- samps$data$samp_info$pat_locs[[i]]
      times <- samps$data$samp_info$pat_times[[i]]
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
