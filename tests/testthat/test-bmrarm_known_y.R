test_that("bmrarm with intercept works as expected", {
  N_iter <- 10
  keep <- matrix(NA, ncol = 14, nrow = N_iter)
  i <- 2

  for(i in 1:N_iter) {
    sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = F, seed = i,
                              slope = F, ar_cov = F)
    samps <- baseline_bmr_known_y(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                          ordinal_outcome = "y_ord", patient_var = "pat_idx",
                          random_slope = F, time_var = "time", ar_cov = F,
                          burn_in = 250, nsim = 500, thin = 1, seed = 3)

    keep[i, 1:2] <- rowMeans(samps$res_cuts[4:5, ])
    keep[i, 3:6] <- rowMeans(samps$res_beta)
    keep[i, 7:10] <- rowMeans(samps$res_sigma)
    keep[i, 11:14] <- rowMeans(samps$res_pat_sig)
    print(i)
  }

  expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.02)
  expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.02)
  expect_equal(colMeans(keep)[11:14], as.numeric(sim_data$sig_alpha), tolerance = 0.04)
})

test_that("bmrarm with slope works as expected", {
  N_iter <- 40
  keep <- matrix(NA, ncol = 14, nrow = N_iter)
  i <- 2

  for(i in 1:N_iter) {
    sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = T, seed = i,
                              slope = T, ar_cov = F)
    samps <- baseline_bmr_known_y(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                          ordinal_outcome = "y_ord", patient_var = "pat_idx",
                          random_slope = T, time_var = "time", ar_cov = F,
                          burn_in = 500, nsim = 1000, thin = 1, seed = 3)

    keep[i, 1:2] <- rowMeans(samps$res_cuts[4:5, ])
    keep[i, 3:6] <- rowMeans(samps$res_beta)
    keep[i, 7:10] <- rowMeans(samps$res_sigma)
    keep[i, 11:14] <- diag(matrix(rowMeans(samps$res_pat_sig), ncol = 4))
    print(i)
  }

  expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.02)
  expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.02)
  expect_equal(colMeans(keep)[11:14],
               as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.04)
})

test_that("bmrarm with ar covariance works as expected", {
  N_iter <- 10
  keep <- matrix(NA, ncol = 14, nrow = N_iter)
  i <- 2

  for(i in 1:N_iter) {
    sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = T, seed = i,
                              slope = T, ar_cov = F)
    samps <- baseline_bmr_known_y(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                          ordinal_outcome = "y_ord", patient_var = "pat_idx",
                          random_slope = T, time_var = "time", ar_cov = T,
                          burn_in = 500, nsim = 1000, thin = 1, seed = 3)

    keep[i, 1:2] <- rowMeans(samps$res_cuts[4:5, ])
    keep[i, 3:6] <- rowMeans(samps$res_beta)
    keep[i, 7:10] <- rowMeans(samps$res_sigma)
    keep[i, 11:14] <- diag(matrix(rowMeans(samps$res_pat_sig), ncol = 4))
    print(i)
  }

  expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.05)
  expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.05)
  expect_equal(colMeans(keep)[11:14],
               as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.06)
})
