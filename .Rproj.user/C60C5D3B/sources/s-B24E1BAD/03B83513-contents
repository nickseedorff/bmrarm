#' Extract conditional mean corresponsing to ordinal components
#'
#' @param cov_mat covariance matrix based on first, last, or middle obs
#' @param locs vector of locations to be sampled
#' @return scalar

pre_calc_ar <- function(cov_mat, locs) {

  ## Pre calculation for mean vector
  mean_pre <- cov_mat[locs, -(locs)] %*% qr.solve(cov_mat[-(locs), -(locs)])

  ## Marginal covariance
  marg_cov <- cov_mat[-(locs), -(locs)]
  marg_inv_cov <- qr.solve(marg_cov)

  ## Conditional covariance
  cond_cov <- cov_mat[locs, locs] - cov_mat[locs, -(locs)] %*%
    marg_inv_cov %*% cov_mat[-(locs), locs]

  list(mean_pre = mean_pre, cond_cov = cond_cov, marg_cov = marg_cov)
}

#' Create list of covariance matrices and optionally their inverse
#'
#' @param env parent environment

get_sig_list <- function(cur_draws, samp_info) {
  chol_list <- marg_cov_list <- sig_list <- sig_inv_list <-
    cond_cov_inv_list <- mean_pre_list <- cond_cov_list <- time_inv <-
    time_det <- list()
  sigma <- cur_draws$sigma
  inv_sigma <- chol2inv(chol(sigma))
  for(i in 1:length(samp_info$uni_pat_times)) {
    ## Get full covariance matrix
    dist_mat <- cur_draws$ar ^ samp_info$uni_dist_mat[[i]]
    time_inv[[i]] <- chol2inv(chol(dist_mat))
    time_det[[i]] <- determinant(dist_mat, logarithm = T)[[1]][1]
    sig_tmp <- kronecker(sigma, dist_mat)
    sig_list[[i]] <- sig_tmp
    sig_inv_list[[i]] <- kronecker(inv_sigma, time_inv[[i]])

    ## Calculate pre-calcs and other covariance matrices
    res_tmp <-  pre_calc_ar(sig_tmp,
                            locs = 1:length(samp_info$uni_pat_times[[i]]))
    mean_pre_list[[i]] <- res_tmp$mean_pre
    cond_cov_list[[i]] <- res_tmp$cond_cov
    cond_cov_inv_list[[i]] <- chol2inv(chol(res_tmp$cond_cov))
    marg_cov_list[[i]] <- res_tmp$marg_cov
  }
  list(sig_list = sig_list,
       mean_pre_list = mean_pre_list, cond_cov_list = cond_cov_list,
       cond_cov_inv_list = cond_cov_inv_list, sig_inv_list = sig_inv_list,
       time_inv = time_inv, time_det = time_det,
       marg_cov_list = marg_cov_list)
}

#' Create list of covariance matrices and optionally their inverse
#'
#' @param env parent environment

get_sig_list2 <- function(cur_draws, samp_info) {
  chol_list <- marg_cov_list <- sig_list <- sig_inv_list <-
    cond_cov_inv_list <- mean_pre_list <- cond_cov_list <- time_inv <-
    time_det <- list()
  sigma <- cur_draws$sigma
  inv_sigma <- chol2inv(chol(sigma))
  sig_alpha_inv <- chol2inv(chol(cur_draws$pat_sig))
  for(i in 1:length(samp_info$uni_pat_times)) {
    ## Get full covariance matrix
    dist_mat <- cur_draws$ar ^ samp_info$uni_dist_mat[[i]]
    pat_ind <- min(which(samp_info$pat_time_ind == i))
    time_inv[[i]] <- chol2inv(chol(dist_mat))
    time_det[[i]] <- determinant(dist_mat, logarithm = T)[[1]][1]
    sig_tmp <- kronecker(sigma, dist_mat)  +
      samp_info$pat_z_kron[[pat_ind]] %*% cur_draws$pat_sig %*% t(samp_info$pat_z_kron[[pat_ind]])
    sig_list[[i]] <- sig_tmp
    sig_inv_list[[i]] <- chol2inv(chol(sig_tmp))

    ## Calculate pre-calcs and other covariance matrices
    res_tmp <-  pre_calc_ar(sig_tmp,
                            locs = 1:length(samp_info$uni_pat_times[[i]]))
    mean_pre_list[[i]] <- res_tmp$mean_pre
    cond_cov_list[[i]] <- res_tmp$cond_cov
    cond_cov_inv_list[[i]] <- chol2inv(chol(res_tmp$cond_cov))
    marg_cov_list[[i]] <- res_tmp$marg_cov
  }
  list(sig_list = sig_list,
       mean_pre_list = mean_pre_list, cond_cov_list = cond_cov_list,
       cond_cov_inv_list = cond_cov_inv_list, sig_inv_list = sig_inv_list,
       time_inv = time_inv, time_det = time_det,
       marg_cov_list = marg_cov_list)
}

#' Create AR matrix
#'
#' @param env parent environment

bmrarm_start <- function(env) {
  e <- environment()
  list2env(as.list(env, all.names = T), envir = e)

  ## Extract covariance matrix, outcome names, stacked covariance matrix
  X <- model.matrix.lm(as.formula(formula), data = data, na.action = "na.pass")
  time_vec <- data[, time_var]

  ## Extract outcomes
  out_vars <- setdiff(all.vars(formula), colnames(X))
  ord_loc <- which(out_vars == ordinal_outcome)
  z <- data[, out_vars[ord_loc]]
  y_cont <- data[, out_vars[-ord_loc]]
  y <- cbind(z, y_cont)
  miss_mat <- is.na(y)

  ## Counts and constants
  pat_idx <- as.numeric(as.factor(data[, patient_var]))
  N_obs <- nrow(y)
  N_outcomes <- ncol(y)
  pat_idx_long <- rep(pat_idx, N_outcomes)
  N_cat <- length(setdiff(unique(z), NA))
  N_pat <- length(unique(pat_idx))
  N_covars <- ncol(X)
  N_pat_effects <- 1 + as.numeric(random_slope)

  ## Extract matrices for multiplication of random effects
  Z <- matrix(rep(1, N_obs), ncol = 1)
  if (random_slope) {
    Z <- cbind(Z, X[, time_var, drop = FALSE])
  }
  Z_kron <- kronecker(diag(rep(1, N_outcomes)), Z)

  ## Indexing for each patient
  pat_locs <- pat_N_obs <- pat_cont_miss_locs <- pat_X <- pat_dist_mat <-
    pat_cont_miss_ranks <- pat_z_kron <- pat_times <- pat_all_locs <-
    uni_pat_time_ranks <- uni_pat_times <- uni_dist_mat <- list()
  pat_time_ind <- rep(NA, N_pat)
  for(i in 1:N_pat) {
    ## Index locations
    pat_locs[[i]] <- which(pat_idx == i)
    pat_all_locs[[i]] <- which(pat_idx_long == i)
    pat_N_obs[[i]] <- length(pat_locs[[i]])
    pat_cont_miss_locs[[i]] <- pat_locs[[i]][miss_mat[pat_locs[[i]], 2]] +
      N_obs
    pat_cont_miss_ranks[[i]] <-
      rank(pat_all_locs[[i]])[pat_all_locs[[i]] %in% pat_cont_miss_locs[[i]]]
    pat_z_kron[[i]] <- kronecker(diag(rep(1, N_outcomes)), Z[pat_locs[[i]], ])
    pat_X[[i]] <- X[pat_locs[[i]], ]
    pat_times[[i]] <- time_vec[pat_locs[[i]]]

    # Unique combinations of outcomes -----------------------------------------

    pat_time <- time_vec[pat_locs[[i]]]
    time_string <- paste0(pat_time, collapse = "")

    if(i == 1) {
      ## Unique times for all observations
      uni_string <- time_string
      uni_pat_times[[length(uni_string)]] <- pat_time
      uni_pat_time_ranks[[length(uni_string)]] <- rank(pat_time)
      pat_time_ind[i] <- 1
      uni_dist_mat[[i]] <- as.matrix(dist(pat_time, diag = T, upper = T))
    } else if (!time_string %in% uni_string) {
      ## Unique observations
      uni_string <- c(uni_string, time_string)
      uni_pat_times[[length(uni_string)]] <- pat_time
      uni_pat_time_ranks[[length(uni_string)]] <- rank(pat_time)
      pat_time_ind[i] <- length(uni_string)
      uni_dist_mat[[length(uni_string)]] <-
        as.matrix(dist(pat_time, diag = T, upper = T))
    } else {
      pat_time_ind[i] <- which(time_string == uni_string)
    }

    # Students that need to be integrated for in Cowles solution --------------

    if (max(z[pat_locs[[i]]], na.rm = T) > 2) {
      if(!exists("pats_for_cut_probs")) {
        pats_for_cut_probs <- i
      } else {
        pats_for_cut_probs <- c(pats_for_cut_probs, i)
      }
    }
  }

  samp_info <- list(pat_locs = pat_locs, pat_N_obs = pat_N_obs,
                    pat_cont_miss_locs = pat_cont_miss_locs,
                    pat_cont_miss_ranks = pat_cont_miss_ranks,
                    pat_z_kron = pat_z_kron, miss_mat = miss_mat,
                    pat_idx = pat_idx, N_obs = N_obs, N_pat = N_pat,
                    beta_sig_prior = diag(rep(sig_prior, N_covars)),
                    N_outcomes = N_outcomes, N_covars = N_covars,
                    N_pat_effects = N_pat_effects,  N_cat = N_cat,
                    sd_c = sd_vec[1], sd_ar = sd_vec[2], sd_pat_sd = sd_vec[3:6],
                    pats_for_cut_probs = pats_for_cut_probs,
                    N_pats_for_probs = length(pats_for_cut_probs),
                    slope = random_slope,
                    pat_times = pat_times, pat_X = pat_X,
                    pat_idx_long =  pat_idx_long,
                    pat_dist_mat = pat_dist_mat, pat_all_locs = pat_all_locs,
                    ar_cov = ar_cov, pat_time_ind = pat_time_ind,
                    uni_pat_time_ranks = uni_pat_time_ranks,
                    uni_pat_times = uni_pat_times,
                    uni_dist_mat = uni_dist_mat)

  ## Store values
  res_ar <- rep(NA, nsim)
  res_beta <- matrix(NA, nrow = N_covars * N_outcomes, ncol = nsim)
  res_pat_eff <- array(NA, dim = c(N_pat, N_pat_effects * N_outcomes, nsim))
  res_pat_sig <- res_pat_sig_q <- matrix(NA, nrow = (N_outcomes * N_pat_effects) ^ 2, nsim)
  res_pat_sig_sd <- matrix(NA, nrow = N_pat_effects * N_outcomes, nsim)
  res_sig_sd <- matrix(NA, nrow = N_outcomes, nsim)
  res_sig <- matrix(NA, nrow = N_outcomes ^ 2, ncol = nsim)
  res_cuts <- matrix(NA, nrow = N_cat + 1, ncol = nsim)
  res_y <- array(NA, c(nrow(y), N_outcomes, nsim))
  res_accept = matrix(NA, ncol = 2, nrow = nsim)
  colnames(res_accept) <- c("cuts", "mh_ar")

  ## Initialize
  res_beta[, 1] <- res_pat_eff[,, 1] <- 0
  res_ar[1] <- 0
  res_pat_sig[, 1] <- pat_sig_q <- diag(rep(0.1, N_pat_effects * N_outcomes))
  pat_sig_sd <- rep(1, N_pat_effects * N_outcomes)
  res_sig[, 1] <- diag(rep(1, N_outcomes))
  res_cuts[, 1] <- start_cuts(N_cat, fixed = 2)

  ## List of current parameter values
  cur_draws <- list(
    cuts = res_cuts[, 1], sigma = matrix(res_sig[, 1], ncol = 2),
    beta = matrix(res_beta[, 1], ncol = N_outcomes), ar = res_ar[1],
    pat_effects = res_pat_eff[,, 1],
    pat_sig = matrix(res_pat_sig[, 1], ncol = N_pat_effects * N_outcomes),
    pat_sig_sd = pat_sig_sd, pat_sig_q = pat_sig_q,
    sig_sd = pat_sig_sd[1:N_outcomes])

  ## Return all objects to calling environment
  list2env(as.list(e, all.names = T), envir = env)

}

#' Create AR matrix
#'
#' @param env parent environment

get_sym_mat <- function(cur_draws, dist_mat) {
  #P1 <- exp(-rho * as.matrix(dist(times, diag = T, upper = T)))
  sigma <- cur_draws$sigma
  P1 <- cur_draws$ar ^ dist_mat
  kronecker(cur_draws$sigma, cur_draws$ar ^ dist_mat)
}

#' Starting values for cut points
#'
#' @param N_cats vector, number of categories for each ordinal variable

start_cuts <- function(N_cats, fixed = 1) {
  starts <- matrix(NA, nrow = max(N_cats) + 1, ncol = length(N_cats))
  for(i in 1:length(N_cats)){
    if(N_cats[i] == 2) {
      starts[1:(N_cats[i] + 1), i] <- c(-Inf, 0, Inf)
    } else if(fixed == 1) {
      starts[1:(N_cats[i] + 1), i] <- c(-Inf, 0, (exp(rnorm(1, sd = 0.25))) *
                                          (1:(N_cats[i] - 2)), Inf)
    } else if (N_cats[i] == 3) {
      starts[1:(N_cats[i] + 1), i] <- c(-Inf, 0, 1, Inf)
    } else {
      starts[1:(N_cats[i] + 1), i] <- c(-Inf, 0, 1,
                                        (exp(rnorm(1, sd = 0.25)) + 1) *
                                          (1:(N_cats[i] - 3)), Inf)
    }
  }
  starts
}

#' Starting values for cut points
#'
#' @param N_cats vector, number of categories for each ordinal variable

first_cut_prob <- function(samp_info, cuts, cuts_tmp) {
  first_prob <- 1
  sd_c <- samp_info$sd_c
  for(i in 4:samp_info$N_cat) {
    num <- pnorm((cuts[i + 1] - cuts[i]) / sd_c) -
      pnorm((cuts_tmp[i - 1] - cuts[i]) / sd_c)
    denom <- pnorm((cuts_tmp[i + 1] - cuts_tmp[i]) / sd_c) -
      pnorm((cuts[i - 1] - cuts_tmp[i]) / sd_c)
    first_prob <- first_prob * num / denom
  }
  first_prob
}
