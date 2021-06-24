## Generate full sigma matrix

for(j in 1:500) {
  N_pat <- samp_info$N_pat
  sig_alpha_inv <- chol2inv(chol(cur_draws$pat_sig))
  sig_inv <- chol2inv(chol(cur_draws$sigma))
  resid_mat <- y - X %*% cur_draws$beta
  N_pat_eff <- ncol(samp_info$pat_z_kron[[1]])

  ## Sigma inverses only needed is autoregressive covariance matrix
  if(samp_info$ar_cov) {
    sig_list <- get_sig_list(cur_draws, samp_info)
  }

  ## Patient effects
  res <- matrix(NA, nrow = N_pat, ncol = N_pat_eff)

  for(i in 1:N_pat) {
    ## Get locations and time matrix for patient
    locs <- samp_info$pat_locs[[i]]
    time_ind <- samp_info$pat_time_ind[i]
    resid_vec <- as.numeric(resid_mat[locs, ])
    pat_Z <- samp_info$pat_z_kron[[i]]

    ## Patient specific covariance matrix
    if(samp_info$ar_cov) {
      pat_sig_inv <- sig_list$sig_inv_list[[time_ind]]
    } else {
      pat_sig_inv <- kronecker(sig_inv, diag(rep(1, samp_info$pat_N_obs[[i]])))
    }

    ## Cross products, covariance, alpha hat
    Z_sig_prod <- crossprod(pat_Z, pat_sig_inv)
    post_cov <- chol2inv(chol(Z_sig_prod %*% pat_Z + sig_alpha_inv))
    post_mean <- post_cov %*% Z_sig_prod %*% resid_vec
    L <- t(chol(post_cov))
    res[i, ] <- L %*% rnorm(length(post_mean)) + post_mean
  }
  pat_sig <- rinvwishart(N_pat + 4,
                         crossprod(res) + diag(rep(1, 4)))

  vals <- list(pat_effects = res, pat_sig = pat_sig)
  res_pat_eff[,, j] <- cur_draws$pat_effects <- vals$pat_effects
  res_pat_sig[, j] <- cur_draws$pat_sig <- vals$pat_sig
}

plot(res_pat_sig[1, ], type = "l")
plot(res_pat_sig[6, ], type = "l")
plot(res_pat_sig[11, ], type = "l")
plot(res_pat_sig[16, ], type = "l")



###############################################################################
## Generate full sigma matrix
###############################################################################

for(j in 1:500) {
  N_pat <- samp_info$N_pat
  sig_alpha_inv <- chol2inv(chol(cur_draws$pat_sig))
  sig_inv <- chol2inv(chol(cur_draws$sigma))
  resid_mat <- y - X %*% cur_draws$beta
  N_pat_eff <- ncol(samp_info$pat_z_kron[[1]])

  ## Sigma inverses only needed is autoregressive covariance matrix
  if(samp_info$ar_cov) {
    sig_list <- get_sig_list(cur_draws, samp_info)
  }

  ## Patient effects
  res <- matrix(NA, nrow = N_pat, ncol = N_pat_eff)

  for(i in 1:N_pat) {
    ## Get locations and time matrix for patient
    locs <- samp_info$pat_locs[[i]]
    time_ind <- samp_info$pat_time_ind[i]
    resid_vec <- as.numeric(resid_mat[locs, ])
    pat_Z <- samp_info$pat_z_kron[[i]]

    ## Patient specific covariance matrix
    if(samp_info$ar_cov) {
      pat_sig_inv <- sig_list$sig_inv_list[[time_ind]]
    } else {
      pat_sig_inv <- kronecker(sig_inv, diag(rep(1, samp_info$pat_N_obs[[i]])))
    }

    ## Cross products, covariance, alpha hat
    Z_sig_prod <- crossprod(pat_Z, pat_sig_inv)
    post_cov <- chol2inv(chol(Z_sig_prod %*% pat_Z + sig_alpha_inv))
    post_mean <- post_cov %*% Z_sig_prod %*% resid_vec
    L <- t(chol(post_cov))
    res[i, ] <- L %*% rnorm(length(post_mean)) + post_mean
  }

  ## Update Correlation matrix
  #res <- sim_data$alpha
  sd_mat <- diag(cur_draws$pat_sig_sd)
  #sd_mat <- diag(c(0.29225371, 0.07295176, 0.23017149, 0.05783218))
  pat_sig <- rinvwishart(N_pat + N_pat_eff + 1, crossprod(res) + 4 * sd_mat)
  #pat_sig <- sim_data$sig_alpha

  ## Update SDs
  pat_sig_inv <- chol2inv(chol(pat_sig))
  pat_sig_sd <- 1 / LaplacesDemon::rinvgamma(
    N_pat_eff, shape = (2 + 4) / 2,
    scale = 2 * diag(pat_sig_inv) + 1 / 100000
  )

  vals <- list(pat_effects = res, pat_sig = pat_sig, pat_sig_sd = pat_sig_sd)
res_pat_eff[,, j] <- cur_draws$pat_effects <- vals$pat_effects
res_pat_sig[, j] <- cur_draws$pat_sig <- vals$pat_sig
res_pat_sig_sd[, j] <- cur_draws$pat_sig_sd <- vals$pat_sig_sd
}

plot(res_pat_sig[1, ], type = "l")
plot(res_pat_sig[6, ], type = "l")
plot(res_pat_sig[11, ], type = "l")
plot(res_pat_sig[16, ], type = "l")
plot(res_pat_sig_sd[1, ], type = "l")
plot(res_pat_sig_sd[2, ], type = "l")


###############################################################################
# Scaled inverse wishart --------------------------------------------------
###############################################################################


## Generate full sigma matrix

for(j in 1:500) {
  N_pat <- samp_info$N_pat
  sig_alpha_inv <- chol2inv(chol(cur_draws$pat_sig))
  sig_inv <- chol2inv(chol(cur_draws$sigma))
  resid_mat <- y - X %*% cur_draws$beta
  N_pat_eff <- ncol(samp_info$pat_z_kron[[1]])

  ## Sigma inverses only needed is autoregressive covariance matrix
  if(samp_info$ar_cov) {
    sig_list <- get_sig_list(cur_draws, samp_info)
  }

  ## Patient effects
  res <- matrix(NA, nrow = N_pat, ncol = N_pat_eff)

  for(i in 1:N_pat) {
    ## Get locations and time matrix for patient
    locs <- samp_info$pat_locs[[i]]
    time_ind <- samp_info$pat_time_ind[i]
    resid_vec <- as.numeric(resid_mat[locs, ])
    pat_Z <- samp_info$pat_z_kron[[i]]

    ## Patient specific covariance matrix
    if(samp_info$ar_cov) {
      pat_sig_inv <- sig_list$sig_inv_list[[time_ind]]
    } else {
      pat_sig_inv <- kronecker(sig_inv, diag(rep(1, samp_info$pat_N_obs[[i]])))
    }

    ## Cross products, covariance, alpha hat
    Z_sig_prod <- crossprod(pat_Z, pat_sig_inv)
    post_cov <- chol2inv(chol(Z_sig_prod %*% pat_Z + sig_alpha_inv))
    post_mean <- post_cov %*% Z_sig_prod %*% resid_vec
    L <- t(chol(post_cov))
    res[i, ] <- L %*% rnorm(length(post_mean)) + post_mean
  }

  ## Correlation matrix
  #res <- sim_data$alpha
  sd_inv <- diag(1 / cur_draws$pat_sig_sd)
  pat_sig_q <- rinvwishart(N_pat + N_pat_eff + 1,
                           sd_inv %*% crossprod(res) %*% sd_inv +
                             diag(rep(1, 4)))
  #pat_sig_q <- sim_data$sig_alpha

  ## SD parameters
  accept_vec <- rep(0, N_pat_eff)
  for(i in 1:N_pat_eff) {
    ## Propose new values
    cur_draws2 <- cur_draws
    cur_draws2$pat_sig_sd[i] <- exp(rnorm(1, log(cur_draws2$pat_sig_sd[i]),
                                          sd = samp_info$sd_pat_sd))
    cur_draws2$pat_sig <- diag(cur_draws2$pat_sig_sd) %*% pat_sig_q %*% diag(cur_draws2$pat_sig_sd)

    ## Calculate comparison values
    pat_inv_old <- chol2inv(chol(cur_draws$pat_sig))
    pat_inv_new <- chol2inv(chol(cur_draws2$pat_sig))

    comp_old <- N_pat / 2 * determinant(pat_inv_old, logarithm = T)[[1]][1] -
      0.5 * sum(diag(res %*% pat_inv_old %*% t(res))) +
      dnorm(log(cur_draws$pat_sig_sd[i]), 0, sd = 10000, log = T)

    comp_new <- N_pat / 2 * determinant(pat_inv_new, logarithm = T)[[1]][1] -
      0.5 * sum(diag(res %*% pat_inv_new %*% t(res))) +
      dnorm(log(cur_draws2$pat_sig_sd[i]), 0, sd = 10000, log = T)
    compar_val <- comp_new - comp_old

    if(compar_val >= log(runif(1))) {
      cur_draws$pat_sig_sd[i] <- cur_draws2$pat_sig_sd[i]
      accept_vec[i] <- 1
    }
  }

  vals <-   list(pat_effects = res,
                 pat_sig = diag(cur_draws$pat_sig_sd) %*% pat_sig_q %*% diag(cur_draws$pat_sig_sd),
                 pat_sig_sd = cur_draws$pat_sig_sd, accept_vec = accept_vec, pat_sig_q = pat_sig_q)
  res_pat_eff[,, j] <- cur_draws$pat_effects <- vals$pat_effects
  res_pat_sig[, j] <- cur_draws$pat_sig <- vals$pat_sig
  res_pat_sig_sd[, j] <- cur_draws$pat_sig_sd <- vals$pat_sig_sd
}

plot(res_pat_sig[1, 100:500], type = "l")
plot(res_pat_sig[6, 100:500], type = "l")
plot(res_pat_sig[11, 100:500], type = "l")
plot(res_pat_sig[16, 100:500], type = "l")
plot(res_pat_sig_sd[1, ], type = "l")
plot(res_pat_sig_sd[2, ], type = "l")
plot(res_pat_sig_sd[3, ], type = "l")
plot(res_pat_sig_sd[4, ], type = "l")

