#' Full conditional draws of the latent continuous values
#'
#' @param bmrvarx_obj an object returned from bmrvarx()
#' @param steps_ahead number of observations to forcast, done recursively
#' @param X matrix of covariates to make forecasts with, number of rows must
#' equal stes_ahead
#' @return matrix
#' @importFrom MASS mvrnorm
#' @export

get_preds <- function(bmrvarx_obj, steps_ahead = 5, X, seed = 15) {
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
  cont_preds <- array(NA, dim = c(N_samp, N_outcomes, steps_ahead))

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
  }

  list(cont_preds_last = cont_preds[, 3, ])
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
