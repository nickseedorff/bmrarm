#' Extract conditional mean corresponding to ordinal components
#'
#' @param y_current current continuous outcome values
#' @param mean_vec vector of mean values
#' @param pre_calc_mat a pre-calculated matrix
#' @param ord_loc locations of the ordinal outcomes
#' @return numeric vector

cond_mean_part <- function(y_current, mean_vec, pre_calc_mat, ord_loc) {
  mean_vec[ord_loc] + pre_calc_mat %*%
    (y_current[-(ord_loc)] - mean_vec[-(ord_loc)])
}


#' Extract conditional means and variances for
#'
#' @param cov_mat covariance matrix
#' @param samp_locs list of unique locations to be sampled
#' @return scalar

pre_calc <- function(cov_mat, samp_locs) {
  lst_use <- list()
  for(i in 1:length(samp_locs)){
    locs <- samp_locs[[i]]

    ## Pre calculation for mean vector
    mean_pre <- cov_mat[locs, -(locs)] %*% qr.solve(cov_mat[-(locs), -(locs)])

    ## Conditional covariance
    cond_cov <- cov_mat[locs, locs] - cov_mat[locs, -(locs)] %*%
      qr.solve(cov_mat[-(locs), -(locs)]) %*% cov_mat[-(locs), locs]

    if(length(locs) == 1) {
      cond_cov <- sqrt(cond_cov)
    }
    lst_use[[i]] <- list(mean_pre = mean_pre, cond_cov = cond_cov,
                         cond_chol_mat = t(chol(cond_cov)),
                         cond_chol_inv = chol2inv(chol(cond_cov)))
  }
  lst_use
}


#' Create list of sampling info to pass to functions
#'
#' @param env parent environment

get_sampling_info <- function(env) {
  list2env(as.list(env, all.names = T), envir = environment())

  ## Info about what outcomes to sample
  samp_type <- vector(length = N_obs)
  for(i in 1:N_obs) {
    missing_locations <- which(miss_mat[i, ] == 1)
    tmp_locs <- c(1:N_ord, setdiff(missing_locations, 1:N_ord))
    locs_char <- paste0(tmp_locs, collapse = ", ")
    if(i == 1) {
      samp_unique <- locs_char
      samp_locs <- list(tmp_locs) ## Unique combinations of values to sample
      samp_length <- length(tmp_locs)
    } else if (!locs_char %in% samp_unique){
      samp_unique <- c(samp_unique, locs_char)
      samp_locs[[length(samp_unique)]] <- tmp_locs
      samp_length <- c(samp_length, length(tmp_locs))
    }
    ## Used to index which combinations of values to sample
    samp_type[i] <- which(locs_char == samp_unique)
  }

  list(samp_type = samp_type, samp_locs = samp_locs, samp_length = samp_length,
       num_ord = N_ord, max_iter = max_iter_rej, N_response = N_response,
       burn_in = burn_in, N_obs = N_obs, N_covars = N_covars,
       N_burn_trunc = N_burn_trunc, nsim = nsim,
       prior_sig = diag(rep(1 / sig_prior, N_covars + N_response)))
}

#' Create lists, vectors, and matrices for storage and passing to functions
#'
#' @param env parent environment

create_storage <- function(env) {
  e <- environment()
  list2env(as.list(env, all.names = T), envir = e)

  ## Store values
  res_M <- res_sigma <- matrix(NA, nrow = N_response ^ 2, ncol = nsim)
  res_beta <- matrix(NA, nrow = N_response * N_covars, ncol = nsim)
  res_cuts <- array(NA, c(max(N_cat) + 1, nsim, N_ord))
  res_latent <- array(NA, c(N_obs, nsim, N_ord))
  res_y <- array(NA, c(N_obs, N_response, nsim))

  ## Initialize
  res_cuts[, 1, ] <- start_cuts(N_cat)
  res_beta[, 1] <- beta_tmp <- matrix(0, ncol = N_response, nrow = N_covars)
  M_tmp <- sig_tmp <- res_M[, 1] <- res_sigma[, 1] <- diag(rep(1, N_response))
  cor_mat <- diag(rep(1, N_response))
  y_use <- as.matrix(dat_out)

  ## Initial values for unidentifiable parameters
  cuts_tmp <- res_cuts[, 1, ,drop = FALSE]
  mean_mat <- matrix(NA, nrow = N_obs, ncol = N_response)
  y_use <- as.matrix(dat_out)

  ## List of current parameter values
  tmp_list <- list(cuts = cuts_tmp, sigma = sig_tmp,
                   beta = beta_tmp, M = M_tmp)
  rm(cuts_tmp, sig_tmp, beta_tmp, M_tmp)

  ## Return all objects to calling enviroment
  list2env(as.list(e, all.names = T), envir = env)
}
