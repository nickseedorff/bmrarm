#' Full conditional draws of the latent continuous values
#'
#' @param bmrvarx_obj an object returned from bmrvarx()
#' @param steps_ahead number of observations to forcast, done recursively
#' @param X matrix of covariates to make forecasts with, number of rows must
#' equal stes_ahead
#' @return matrix
#' @importFrom MASS mvrnorm
#' @export

get_preds_bmrvarx <- function(bmrvarx_obj, steps_ahead = 5, X, seed = 15, get_ord_preds = T) {
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

#' Get cuts list
#' @param y continuous value to be discretized

cont_to_cat <- function(y, cut_points) {
  for(i in 1:(length(cut_points) - 1)) {
    if(y > cut_points[i] & y <= cut_points[i + 1]) {
      return(i)
    }
  }
}


#' Get DIC
#'
#' @param y_current current continous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc number of ordinal outcomes
#' @return scalar

get_preds_bmrarm <- function(samps) {

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

