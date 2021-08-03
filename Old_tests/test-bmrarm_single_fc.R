test_that("bmrarm single fcs work with AR term works as expected", {

  test_that("bmrarm with slope works as expected", {
    N_iter <- 100
    keep <- matrix(NA, ncol = 14, nrow = N_iter)
    i <- 2

    for(i in 1:N_iter) {
      sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = T, seed = i,
                                slope = T, ar_cov = F)
      #sim_data$data$time <- sim_data$data$time / 10
      samps <- baseline_bmr_test(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                                    ordinal_outcome = "y_ord", patient_var = "pat_idx",
                                    random_slope = T, time_var = "time", ar_cov = F,
                                    burn_in = 100, nsim = 250, thin = 1, seed = 3,
                                 which_fc = c(T, F, F, F, F))

      keep[i, 1:2] <- apply(samps$res_cuts[4:5, ], 1, median)
      keep[i, 3:6] <- apply(samps$res_beta, 1, median)
      keep[i, 7:10] <- apply(samps$res_sigma, 1, median)
      keep[i, 11:14] <- diag(matrix(apply(samps$res_pat_sig, 1, median), ncol = 4))
      print(i)
      if(i > 1) print(colMeans(keep[1:i, ]))
    }

    expect_equal(colMeans(keep)[1:2], sim_data$cuts[4:5], tolerance = 0.05)
    expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.02)
    expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.02)
    expect_equal(colMeans(keep)[11:14],
                 as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.04)
  })
})

test_that("bmrarm single fcs work with AR term works as expected", {

  test_that("bmrarm with slope works as expected", {
    N_iter <- 100
    keep <- matrix(NA, ncol = 14, nrow = N_iter)
    i <- 2

    for(i in 1:N_iter) {
      sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = T, seed = i,
                                slope = T, ar_cov = F)
      #sim_data$data$time <- sim_data$data$time / 10
      samps <- baseline_bmr_test(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                                 ordinal_outcome = "y_ord", patient_var = "pat_idx",
                                 random_slope = T, time_var = "time", ar_cov = F,
                                 burn_in = 100, nsim = 250, thin = 1, seed = 3,
                                 which_fc = c(F, F, T, F, F))

      keep[i, 1:2] <- apply(samps$res_cuts[4:5, ], 1, median)
      keep[i, 3:6] <- apply(samps$res_beta, 1, median)
      keep[i, 7:10] <- apply(samps$res_sigma, 1, median)
      keep[i, 11:14] <- diag(matrix(apply(samps$res_pat_sig, 1, median), ncol = 4))
      print(i)
      if(i > 1) print(colMeans(keep[1:i, ]))
    }

    expect_equal(colMeans(keep)[1:2], sim_data$cuts[4:5], tolerance = 0.05)
    expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.02)
    expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.02)
    expect_equal(colMeans(keep)[11:14],
                 as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.04)
  })
})

test_that("bmrarm single fcs work with AR term works as expected", {

  test_that("bmrarm with slope works as expected", {
    N_iter <- 100
    keep <- matrix(NA, ncol = 14, nrow = N_iter)
    i <- 2

    for(i in 1:N_iter) {
      sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = T, seed = i,
                                slope = T, ar_cov = F)
      #sim_data$data$time <- sim_data$data$time / 10
      samps <- baseline_bmr_test(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                                 ordinal_outcome = "y_ord", patient_var = "pat_idx",
                                 random_slope = T, time_var = "time", ar_cov = F,
                                 burn_in = 100, nsim = 250, thin = 1, seed = 3,
                                 which_fc = c(T, F, T, F, F))

      keep[i, 1:2] <- apply(samps$res_cuts[4:5, ], 1, median)
      keep[i, 3:6] <- apply(samps$res_beta, 1, median)
      keep[i, 7:10] <- apply(samps$res_sigma, 1, median)
      keep[i, 11:14] <- diag(matrix(apply(samps$res_pat_sig, 1, median), ncol = 4))
      print(i)
      if(i > 1) print(colMeans(keep[1:i, ]))
    }

    expect_equal(colMeans(keep)[1:2], sim_data$cuts[4:5], tolerance = 0.05)
    expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.02)
    expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.02)
    expect_equal(colMeans(keep)[11:14],
                 as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.04)
  })
})

test_that("bmrarm single fcs work with AR term works as expected", {

  test_that("bmrarm with slope works as expected", {
    N_iter <- 150
    keep <- matrix(NA, ncol = 14, nrow = N_iter)
    i <- 2

    for(i in 1:N_iter) {
      sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = T, seed = i,
                                slope = T, ar_cov = F)
      #sim_data$data$time <- sim_data$data$time / 10
      samps <- baseline_bmr_test(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                                 ordinal_outcome = "y_ord", patient_var = "pat_idx",
                                 random_slope = T, time_var = "time", ar_cov = T,
                                 burn_in = 200, nsim = 400, thin = 1, seed = 3,
                                 which_fc = c(T, F, F, F, F))

      keep[i, 1:2] <- apply(samps$res_cuts[4:5, ], 1, median)
      keep[i, 3:6] <- apply(samps$res_beta, 1, median)
      keep[i, 7:10] <- apply(samps$res_sigma, 1, median)
      keep[i, 11:14] <- diag(matrix(apply(samps$res_pat_sig, 1, median), ncol = 4))
      print(i)
      if(i > 1) print(colMeans(keep[1:i, ]))
    }

    expect_equal(colMeans(keep)[1:2], sim_data$cuts[4:5], tolerance = 0.05)
    expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.02)
    expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.02)
    expect_equal(colMeans(keep)[11:14],
                 as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.04)
  })
})

test_that("bmrarm single fcs work with AR term works as expected", {

  test_that("bmrarm with slope works as expected", {
    N_iter <- 25
    keep <- matrix(NA, ncol = 14, nrow = N_iter)
    i <- 2

    for(i in 1:N_iter) {
      sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = F, seed = i,
                                slope = T, ar_cov = F)
      samps <- baseline_bmr_test(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                                 ordinal_outcome = "y_ord", patient_var = "pat_idx",
                                 random_slope = T, time_var = "time", ar_cov = F,
                                 burn_in = 250, nsim = 500, thin = 1, seed = 3,
                                 which_fc = c(T, F, F, T, F))

      keep[i, 1:2] <- apply(samps$res_cuts[4:5, ], 1, median)
      keep[i, 3:6] <- apply(samps$res_beta, 1, median)
      keep[i, 7:10] <- apply(samps$res_sigma, 1, median)
      keep[i, 11:14] <- diag(matrix(apply(samps$res_pat_sig, 1, median), ncol = 4))
      print(i)
      if(i > 1) print(colMeans(keep[1:i, ]))
    }

    expect_equal(colMeans(keep)[1:2], sim_data$cuts[4:5], tolerance = 0.05)
    expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.02)
    expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.02)
    expect_equal(colMeans(keep)[11:14],
                 as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.04)
  })
})

test_that("bmrarm single fcs work with AR term works as expected", {

  test_that("bmrarm with slope works as expected", {
    N_iter <- 25
    keep <- matrix(NA, ncol = 14, nrow = N_iter)
    i <- 2

    for(i in 1:N_iter) {
      sim_data <- gen_ar_errors(N = 20, N_pat = 48, unequal = F, seed = i,
                                slope = T, ar_cov = F)
      samps <- baseline_bmr_test(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                                 ordinal_outcome = "y_ord", patient_var = "pat_idx",
                                 random_slope = T, time_var = "time", ar_cov = F,
                                 burn_in = 500, nsim = 1000, thin = 1, seed = 3,
                                 which_fc = c(T, F, T, T, F))

      keep[i, 1:2] <- apply(samps$res_cuts[4:5, ], 1, median)
      keep[i, 3:6] <- apply(samps$res_beta, 1, median)
      keep[i, 7:10] <- apply(samps$res_sigma, 1, median)
      keep[i, 11:14] <- diag(matrix(apply(samps$res_pat_sig, 1, median), ncol = 4))
      print(i)
      if(i > 1) print(colMeans(keep[1:i, ]))
    }

    expect_equal(colMeans(keep)[1:2], sim_data$cuts[4:5], tolerance = 0.05)
    expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.02)
    expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.02)
    expect_equal(colMeans(keep)[11:14],
                 as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.04)
  })
})

test_that("bmrarm single fcs work with AR term works as expected", {

  test_that("bmrarm with slope works as expected", {
    N_iter <- 150
    keep <- matrix(NA, ncol = 14, nrow = N_iter)
    i <- 2

    for(i in 1:N_iter) {
      sim_data <- gen_ar_errors(N = 7, N_pat = 98, unequal = F, seed = i,
                                slope = T, ar_cov = F)
      #sim_data$data$time <- sim_data$data$time / 10
      samps <- baseline_bmr_test(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                                 ordinal_outcome = "y_ord", patient_var = "pat_idx",
                                 random_slope = T, time_var = "time", ar_cov = T,
                                 burn_in = 500, nsim = 1000, thin = 1, seed = 3,
                                 which_fc = c(T, T, T, F, F))

      keep[i, 1:2] <- apply(samps$res_cuts[4:5, ], 1, median)
      keep[i, 3:6] <- apply(samps$res_beta, 1, median)
      keep[i, 7:10] <- apply(samps$res_sigma, 1, median)
      keep[i, 11:14] <- diag(matrix(apply(samps$res_pat_sig, 1, median), ncol = 4))
      print(cov2cor(matrix(apply(samps$res_pat_sig, 1, median), ncol = 4)))
      cov2cor(sim_data$sig_alpha)
      print(i)
      if(i > 1) print(colMeans(keep[1:i, ]))
    }
    head(apply(samps$res_pat_eff, c(1,2), median),10)
    head(sim_data$alpha, 10)

    expect_equal(colMeans(keep)[1:2], sim_data$cuts[4:5], tolerance = 0.05)
    expect_equal(colMeans(keep)[3:6], as.numeric(sim_data$beta), tolerance = 0.02)
    expect_equal(colMeans(keep)[7:10], as.numeric(sim_data$sigma), tolerance = 0.02)
    expect_equal(colMeans(keep)[11:14],
                 as.numeric(diag(sim_data$sig_alpha)), tolerance = 0.04)
  })
})

