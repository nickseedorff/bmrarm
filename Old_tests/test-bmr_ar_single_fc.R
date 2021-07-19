test_that("bmrarm single fcs work with AR term works as expected", {

# Test regression coef, sigma, and ar ------------------------------------------
  rm(list = ls())
  sim_data <- bmrarm:::gen_ar_errors(N = 15, N_pat = 1000, slope = T, seed = 90, unequal = T)
  formula <- cbind(y_ord, y2) ~ time; data <- sim_data$data; seed = 14;
  ordinal_outcome <- "y_ord"; patient_var <- "pat_idx"; burn_in = 100;
  ar_cov = T; sig_prior = 1000000000; time_var = "time"; nsim = 2000; thin = 10;
  verbose = TRUE; random_slope <- T; true_y <- NULL; sd_vec <- c(0.15, 0.30);

  set.seed(seed)
  bmrarm_start(env = environment())
  y <- cbind(sim_data$truth$y1, sim_data$truth$y2)

  cur_draws$cuts <- sim_data$cuts
  cur_draws$sigma <- sim_data$sigma
  cur_draws$pat_effects <- sim_data$alpha
  cur_draws$ar <- sim_data$ar
  cur_draws$beta <- matrix(sim_data$beta, 2)
  cur_draws$pat_sig <- sim_data$sig_alpha

  ## Test beta
  beta_M <- rowMeans(replicate(50, as.numeric(bmrarm_fc_sig_beta(y, X, Z_kron, cur_draws, samp_info)[[1]])))
  expect_equal(beta_M, sim_data$beta, tolerance = 0.035)

  ## Test Sigma
  sig_M <- rowMeans(replicate(25, as.numeric(bmrarm_fc_sig_beta(y, X, Z_kron, cur_draws, samp_info)[[2]])))
  expect_equal(sig_M, as.numeric(sim_data$sigma), tolerance = 0.025)

  ## Pat sig
  pat_sig_M <- rowMeans(replicate(25, as.numeric(bmrarm_fc_patient(y, z, X, cur_draws, samp_info)[[2]])))
  expect_equal(pat_sig_M, as.numeric(sim_data$sig_alpha), tolerance = 0.02)

  ## Pat effects
  pat_M <- rowMeans(replicate(25, as.numeric(bmrarm_fc_patient(y, z, X, cur_draws, samp_info)[[1]])))
  expect_gt(min(diag(cor(matrix(pat_M, ncol = 4), sim_data$alpha))), 0.5)

  ## AR term
  N_iter <- 300
  ar_res <- rep(NA, N_iter)
  for(i in 1:N_iter) {
    vals <- bmrarm_mh_ar(y, X, Z_kron, cur_draws, samp_info)
    ar_res[i] <- cur_draws$ar <- vals$ar
  }
  expect_equal(mean(ar_res), sim_data$ar, tolerance = 0.02)


# Latent outcomes ---------------------------------------------------------
  for(i in 1:10) {
    rm(list = ls())
    sim_data <- bmrarm:::gen_ar_errors(N = 7, N_pat = 400, slope = T, seed = 90, unequal = T, ar_cov = F)
    formula <- cbind(y_ord, y2) ~ time; data <- sim_data$data; seed = 14;
    ordinal_outcome <- "y_ord"; patient_var <- "pat_idx"; burn_in = 100;
    ar_cov = F; sig_prior = 1000000000; time_var = "time"; nsim = 500; thin = 10;
    verbose = TRUE; random_slope <- T; true_y <- NULL; sd_vec <- c(0.05, 0.30);

    set.seed(seed)
    bmrarm_start(env = environment())
    y <- cbind(sim_data$truth$y1, sim_data$truth$y2)

    cur_draws$cuts <- sim_data$cuts
    cur_draws$sigma <- sim_data$sigma
    cur_draws$pat_effects <- sim_data$alpha
    cur_draws$ar <- sim_data$ar
    cur_draws$beta <- matrix(sim_data$beta, 2)
    cur_draws$pat_sig <- sim_data$sig_alpha

    N_iter <- 500
    for(i in 1:N_iter) {
      y_cuts <- bmrarm_fc_y_cuts(y, z, X, Z_kron, cur_draws, samp_info)
      #print("y")
      y <-  res_y[,, i] <- y_cuts$y
      res_cuts[, i] <- cur_draws$cuts <- y_cuts$cuts
      res_accept[i, 1] <- y_cuts$accept
      if(i %% 50 == 0) print(i)
    }
    plot(res_cuts[5, ], type = "l")
  }


plot(res_cuts[4, ], type = "l")
plot(res_cuts[5, ], type = "l")


})

test_that("bmrarm single fcs work without AR term works as expected", {

  # Test regression coef, sigma, and ar ------------------------------------------
  rm(list = ls())
  sim_data <- bmrarm:::gen_ar_errors(N = 15, N_pat = 1000, slope = T, seed = 90, unequal = T, ar_cov = F)
  formula <- cbind(y_ord, y2) ~ time; data <- sim_data$data; seed = 14;
  ordinal_outcome <- "y_ord"; patient_var <- "pat_idx"; burn_in = 100;
  ar_cov = F; sig_prior = 1000000000; time_var = "time"; nsim = 2000; thin = 10;
  verbose = TRUE; random_slope <- T; true_y <- NULL; sd_vec <- c(0.15, 0.30);

  set.seed(seed)
  bmrarm_start(env = environment())
  y <- cbind(sim_data$truth$y1, sim_data$truth$y2)

  cur_draws$cuts <- sim_data$cuts
  cur_draws$sigma <- sim_data$sigma
  cur_draws$pat_effects <- sim_data$alpha
  cur_draws$ar <- sim_data$ar
  cur_draws$beta <- matrix(sim_data$beta, 2)
  cur_draws$pat_sig <- sim_data$sig_alpha

  ## Test beta
  beta_M <- rowMeans(replicate(50, as.numeric(bmrarm_fc_sig_beta(y, X, Z_kron, cur_draws, samp_info)[[1]])))
  expect_equal(beta_M, sim_data$beta, tolerance = 0.025)

  ## Test Sigma
  sig_M <- rowMeans(replicate(25, as.numeric(bmrarm_fc_sig_beta(y, X, Z_kron, cur_draws, samp_info)[[2]])))
  expect_equal(sig_M, as.numeric(sim_data$sigma), tolerance = 0.025)

  ## Pat sig
  pat_sig_M <- rowMeans(replicate(25, as.numeric(bmrarm_fc_patient(y, z, X, cur_draws, samp_info)[[2]])))
  expect_equal(pat_sig_M, as.numeric(sim_data$sig_alpha), tolerance = 0.02)

  ## Pat effects
  pat_M <- rowMeans(replicate(25, as.numeric(bmrarm_fc_patient(y, z, X, cur_draws, samp_info)[[1]])))
  expect_gt(min(diag(cor(matrix(pat_M, ncol = 4), sim_data$alpha))), 0.5)
})


