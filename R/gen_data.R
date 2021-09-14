#' Generate data for time series simulation
#'
#' @param N Number of timepoints
#' @param seed seed to set, default of 10
#' @importFrom MASS mvrnorm
#' @import dplyr

gen_single <- function(N = 610, seed = 10, N_param = 3, first_obs = 101) {
  set.seed(seed)

  ## Additional effects
  cont_covar <- rnorm(N)
  cov_mat <- cbind(1, cont_covar)
  cov_df <- as.data.frame(cov_mat)
  colnames(cov_df) <- c("intercept", "x1")

  beta <- matrix(c(
    0.3, -0.1, 0.2,
    0.4, -0.15, 0.75),
    byrow = T, ncol = 3)

  ## Define parameters
  M <- matrix(c(0.4, 0.15, 0.35,
                0.2, 0.35, 0.25,
                0.12, 0.1, 0.5), ncol = 3, byrow = TRUE)
  eigen(M)
  sigma <- matrix(c(1, 0.12, 0.12,
                    0.12, 1, 0.12,
                    0.12, 0.12, 0.5), ncol = 3, byrow = T)

  if(N_param == 2) {
    M <- matrix(c(0.5, 0.40,
                  0.25, 0.65), ncol = 2, byrow = TRUE)
    eigen(M)
    sigma <- sigma[c(1, 3), c(1, 3)]
    beta <- beta[, c(1, 3)]
  }

  ## Subject and time indexing
  df <- data.frame(obs_num = 0:(N - 1),
                   pat_idx = 1)
  y <- matrix(NA, nrow = N, ncol = N_param)

  ## Generate continuous values
  for(i in 1:N) {
    if(df$obs_num[i] == 0) {
      y[i, ] <- mvrnorm(1, mu = t(beta) %*% cov_mat[i, ],
                        Sigma = sigma * 2)
    } else {
      y[i, ] <- mvrnorm(1, mu = t(beta) %*% cov_mat[i, ] +
                          M %*% y[i - 1, ], Sigma = sigma)
    }
  }

  ## Create ordinal values and binary values
  full <- data.frame(df, y)
  colnames(full)[3:ncol(full)] <- paste0("y", 1:N_param)
  which_locs <- which((1:nrow(full)) >= first_obs)
  full <- full[which_locs, ]

  cuts <- c(-Inf, 0, 1.3, 2.2, Inf)
  full$y_ord <- case_when(
    full$y1 <= cuts[2] ~ 1,
    full$y1 <= cuts[3] ~ 2,
    full$y1 <= cuts[4] ~ 3,
    full$y1 > cuts[4]  ~ 4
  )

  cuts2 <- c(-Inf, 0, 1.2, Inf)
  full$y_ord2 <- case_when(
    full$y2 <= cuts2[2] ~ 1,
    full$y2 <= cuts2[3] ~ 2,
    full$y2 > cuts2[3] ~ 3
  )

  ## Check for reasonable counts
  last_obs <- N - first_obs + 1 - 10
  table(full$y_ord)
  full <- cbind(full, cov_df[which_locs, ])
  full_no_miss <- full
  full$y_ord[sample(1:last_obs, size = 25)] <- NA
  full$y_ord2[sample(1:last_obs, size = 25)] <- NA
  full$y3[sample(1:last_obs, size = 25)] <- NA

  list(data = full, M = M, sigma = sigma, beta = beta, cuts = cuts,
       cuts2 = cuts2, X = cov_mat[which_locs, ], sig0 = sigma * 2,
       data_no_miss = full_no_miss)
}

#' Generate data for longitudinal simulation
#'
#' @param N Mean time points per patient, drawn from a poisson
#' @param N_pat number of patients
#' @param seed seed to set, default of 10
#' @importFrom MASS mvrnorm
#' @import ggplot2
#' @import dplyr
#' @importFrom magic adiag

gen_ar_errors <- function(N = 7, N_pat = 48, seed = 10, unequal = FALSE,
                          slope = F, ar_cov = T) {
  set.seed(seed)

  ## Number of observations
  if(unequal) {
    obs_counts <- sample((N - 3):N, size = N_pat, replace = T,
                         prob = c(0.042, 0.042, 0.083, 0.833))
  } else {
    obs_counts <- rep(N, N_pat)
  }
  pat_idx <- rep(1:length(obs_counts), times = obs_counts)
  pat_idx <- c(pat_idx, pat_idx)
  N_total <- sum(obs_counts)

  ## Additional effects
  for(i in 1:length(obs_counts)) {
    if(i == 1) {
      time_covar <- 0:(obs_counts[i] - 1)
    } else {
      time_covar <- c(time_covar, 0:(obs_counts[i] - 1))
    }
  }
  cov_mat <- cbind(1, time_covar)
  cov_df <- as.data.frame(cov_mat)
  colnames(cov_df) <- c("intercept", "time")

  beta <- c(0.5, 0.18, 0.05, 0.1)
  sigma <- matrix(c(0.3, 0.025,
                    0.025, 0.13), ncol = 2, byrow = T)

  subject_effects <- subj_slopes <- matrix(NA, ncol = 2, nrow = N_pat)

  sig_alpha <- matrix(c(0.22, 0.075,
                        0.075, 0.18), ncol = 2, byrow = TRUE)
  slope_alpha <- matrix(c(0.05, 0.005, 0.005, 0.04), ncol = 2)
  for(i in 1:N_pat) {
    subject_effects[i, ] <- mvrnorm(1, mu = rep(0, 2), Sigma = sig_alpha)
    subj_slopes[i, ] <- mvrnorm(1, mu = rep(0, 2), Sigma = slope_alpha)
  }

  ## Subject and time indexing
  ar_val <- c(0.35)
  y <- true_means <- err_vec <- fixed_means <- rep(NA, N_total * 2)
  kron_X <- kronecker(diag(rep(1, 2)), cov_mat)

  if(slope) {
    kron_Z <- kron_X
    subject_effects <- cbind(subject_effects[, 1], subj_slopes[, 1],
                             subject_effects[, 2], subj_slopes[, 2])
  } else {
    kron_Z <- kron_X[, c(1, 3)]
  }

  cur_draws <- list(sigma = sigma, ar = ar_val)

  for(i in 1:N_pat) {
    locs <- which(pat_idx == i)
    dist_mat <- as.matrix(dist(1:(length(locs) / 2), diag = T, upper = T))
    if(ar_cov) {
    dist_mat2 <- ar_val ^ dist_mat
    } else {
      dist_mat2 <- diag(1, length(locs) / 2)
    }
    sig_full <- kronecker(sigma, dist_mat2)
    fixed_means[locs] <- kron_X[locs, ] %*% beta
    true_means[locs] <- kron_X[locs, ] %*% beta + kron_Z[locs, ] %*% subject_effects[i, ]
    err_vec[locs] <- MASS::mvrnorm(1, mu = rep(0, length(locs)), Sigma = sig_full)
    y[locs] <- true_means[locs] + err_vec[locs]
  }

  ## Create ordinal values and binary values
  full <- data.frame(pat_idx, y, outcome = rep(c(1, 2), each = N_total)) %>%
    group_by(pat_idx, outcome) %>%
    mutate(obs_num = row_number() - 1)
  cuts <- c(-Inf, 0, 1, 1.5, 1.9, Inf)
  full$y_ord <- case_when(
    full$y <= cuts[2] ~ 1,
    full$y <= cuts[3] ~ 2,
    full$y <= cuts[4] ~ 3,
    full$y <= cuts[5]  ~ 4,
    full$y > cuts[5]  ~ 5
  )

  df <- data.frame(y1 = full$y[full$outcome == 1],
                   y_ord = full$y_ord[full$outcome == 1],
                   y2 = full$y[full$outcome == 2],
                   pat_idx = full$pat_idx[full$outcome == 1],
                   time = full$obs_num[full$outcome == 1])

  if(slope) {
    sig_alpha <- magic::adiag(sig_alpha, slope_alpha)
    sig_alpha <- sig_alpha[c(1, 3, 2, 4), c(1, 3, 2, 4)]
  }

  ## Include missing values
  truth <- df
  df[df$time == 2, c("y_ord")] <- NA
  df$y_ord[sample(which(df$time > 2), size = 4)] <- NA
  df$y2[sample(1:nrow(df), size = 9)] <- NA

# Generate additional data for prediction evaluation ----------------------
  ## Data is generated conditional on the previous data
  ## This is because simulations were already run using the previous data
  ## so we want to reuse those model objects and can draw new data by
  ## conditioning on the first set


  ## Generate new covariates
  cov_tmp <- cov_df
  cov_tmp$pat_idx <- pat_idx[1:(length(pat_idx) / 2)]
  new_cov <- cov_tmp %>%
    group_by(pat_idx) %>%
    mutate(old_time = time, time = row_number() + max(old_time)) %>%
    filter(row_number() <= 4)
  pat_idx_pred <- rep(new_cov$pat_idx, 2)
  y_pred <- rep(NA, length(pat_idx_pred))
  new_cov <-  new_cov %>%
    ungroup %>%
    dplyr::select(-old_time, -pat_idx)

  ## Kronecker product matrices
  kron_X_pred <- kronecker(diag(rep(1, 2)), as.matrix(new_cov))
  if(slope) {
    kron_Z_pred <- kron_X_pred
  } else {
    kron_Z_pred <- kron_X_pred[, c(1, 3)]
  }

  for(i in 1:N_pat) {

    ## Locations of old and new outcomes
    old_locs <- which(pat_idx == i)
    l_old_locs <- length(old_locs) / 2
    old_cov_locs <- c(1:l_old_locs, (l_old_locs + 5):(2 * l_old_locs + 4))
    new_locs <- which(pat_idx_pred == i)
    dist_mat <- as.matrix(dist(1:(l_old_locs + 4), diag = T, upper = T))

    ## Covariance matrix
    if(ar_cov) {
      dist_mat2 <- ar_val ^ dist_mat
    } else {
      dist_mat2 <- diag(1, nrow(dist_mat))
    }

    sig_full <- kronecker(sigma, dist_mat2)

    ## Marginal mean of new data
    marg_mean <- kron_X_pred[new_locs, ] %*% beta +
      kron_Z_pred[new_locs, ] %*% subject_effects[i, ]

    ## Conditional mean of new data
    pre_calc <- sig_full[-old_cov_locs, old_cov_locs] %*%
      qr.solve(sig_full[old_cov_locs, old_cov_locs])
    cond_mean <- marg_mean + pre_calc %*% (y[old_locs] - true_means[old_locs])

    ## Conditional covariance of new data
    cond_cov <- sig_full[-old_cov_locs, -old_cov_locs] - pre_calc %*%
      sig_full[old_cov_locs, -old_cov_locs]

    ## Generate new data
    y_pred[new_locs] <- mvrnorm(1, mu = cond_mean, Sigma = cond_cov)
  }

  ## Discretize the ordinal outcome
  full_pred <- data.frame(pat_idx_pred, y_pred,
                          outcome = rep(c(1, 2), each = N_pat * 4)) %>%
    group_by(pat_idx_pred, outcome) %>%
    mutate(obs_num = row_number() - 1)

  full_pred$y_ord <- case_when(
    full_pred$y_pred <= cuts[2] ~ 1,
    full_pred$y_pred <= cuts[3] ~ 2,
    full_pred$y_pred <= cuts[4] ~ 3,
    full_pred$y_pred <= cuts[5]  ~ 4,
    full_pred$y_pred > cuts[5]  ~ 5
  )

  ## Dataframe of data for predictions
  df_pred <- data.frame(y1 = full_pred$y_pred[full_pred$outcome == 1],
                        y_ord = full_pred$y_ord[full_pred$outcome == 1],
                        y2 = full_pred$y_pred[full_pred$outcome == 2],
                        pat_idx = full_pred$pat_idx_pred[full_pred$outcome == 1],
                        time = new_cov$time)

# Storage -----------------------------------------------------------------

  ## Dataframe of observed data and data for prediction evaluation
  full_df <- rbind(df, df_pred) %>%
    arrange(pat_idx, time) %>%
    as.data.frame()

  full_X <- as.data.frame(cbind(1, full_df$time))
  colnames(full_X) <- c("Constant", "time")

  list(data = df, sigma = sigma, beta = beta, cuts = cuts, X = cov_mat,
       alpha = subject_effects, sig_alpha = sig_alpha, true_means = true_means,
       fixed_means = fixed_means, truth = truth, ar = ar_val, err_vec = err_vec,
       pred_data = df_pred, pred_X = new_cov, full_df = full_df,
       full_X = full_X)
}
