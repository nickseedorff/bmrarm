test_that("bmrarm with intercept works as expected", {
  N_iter <- 10
  keep <- matrix(NA, ncol = 14, nrow = N_iter)
  i <- 2

  for(i in 1:N_iter) {
    sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = F, seed = i,
                              slope = F, ar_cov = F)
    samps <- baseline_bmr(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                          ordinal_outcome = "y_ord", patient_var = "pat_idx",
                          random_slope = F, time_var = "time", ar_cov = F,
                          burn_in = 250, nsim = 500, thin = 1, seed = 3)

    keep[i, 1:2] <- rowMeans(samps$res_cuts[4:5, ])
    keep[i, 3:6] <- rowMeans(samps$res_beta)
    keep[i, 7:10] <- rowMeans(samps$res_sigma)
    keep[i, 11:14] <- rowMeans(samps$res_pat_sig)
    print(i)
  }

  expect_equal(colMeans(keep)[1:2], sim_data$cuts[4:5],
               tolerance = 0.05)
  expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta),
               tolerance = 0.05)
  expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma),
               tolerance = 0.05)
  expect_equal(colMeans(keep)[11:14], as.numeric(sim_data$sig_alpha),
               tolerance = 0.02)
})

test_that("bmrarm with slope works as expected", {
  N_iter <- 10
  keep <- matrix(NA, ncol = 14, nrow = N_iter)

  for(i in 1:N_iter) {
    sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = F, seed = i,
                              slope = T, ar_cov = F)
    samps <- baseline_bmr(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                          ordinal_outcome = "y_ord", patient_var = "pat_idx",
                          random_slope = T, time_var = "time", ar_cov = F,
                          burn_in = 250, nsim = 500, thin = 1, seed = 3)

    keep[i, 1:2] <- apply(samps$res_cuts[4:5, ], 1, median)
    keep[i, 3:6] <- apply(samps$res_beta, 1, median)
    keep[i, 7:10] <- apply(samps$res_sigma, 1, median)
    keep[i, 11:14] <- diag(matrix(apply(samps$res_pat_sig, 1, median), ncol = 4))
    print(i)
  }

  expect_equal(colMeans(keep)[1:2], sim_data$cuts[4:5], tolerance = 0.05)
  expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.05)
  expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.05)
  expect_equal(colMeans(keep)[11:14],
               as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.02)
})

test_that("bmrarm with slope works as expected", {
  N_iter <- 10
  keep <- matrix(NA, ncol = 15, nrow = N_iter)
  i <- 3

  for(i in 1:N_iter) {
    sim_data <- gen_ar_errors(N = 6, N_pat = 75, unequal = T, seed = i,
                              slope = T, ar_cov = T)
    samps <- baseline_bmr(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                          ordinal_outcome = "y_ord", patient_var = "pat_idx",
                          random_slope = T, time_var = "time", ar_cov = T,
                          burn_in = 250, nsim = 500, thin = 1, seed = 3,
                          sd_vec = c(0.15, 0.15))

    #f <- paste0("C:\\Users\\nicho\\Dropbox\\_Thesis\\bmrvar\\store_apps_sims\\bmr3_samps", i, ".RDS")
    #saveRDS(samps, f)

    keep[i, 1:2] <- apply(samps$res_cuts[4:5, ], 1, median)
    keep[i, 3:6] <- apply(samps$res_beta, 1, median)
    keep[i, 7:10] <- apply(samps$res_sigma, 1, median)
    keep[i, 11:14] <- diag(matrix(apply(samps$res_pat_sig, 1, median), ncol = 4))
    keep[i, 15] <- median(samps$res_ar)
    print(i)
  }

  expect_equal(colMeans(keep)[1:2], sim_data$cuts[4:5], tolerance = 0.05)
  expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.05)
  expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.05)
  expect_equal(colMeans(keep)[11:14],
               as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.06)
  expect_equal(mean(keep[, 15]), sim_data$ar, tolerance = 0.06)
})
