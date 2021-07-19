test_that("bmrarm works as anticipated", {
  sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = F, seed = 1,
                            slope = F, ar_cov = F)

# Test brm_int model against stored object ---------------------------------

  samps_int <- baseline_bmr(
    formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
    ordinal_outcome = "y_ord", patient_var = "pat_idx",
    random_slope = F, time_var = "time", ar_cov = F,
    burn_in = 250, nsim = 500, thin = 1, seed = 3,
    sd_vec = c(0.15, 0.30, rep(0.2, 4)))

# Test brm_slope model against stored object ------------------------------

  samps_slope <- baseline_bmr(
    formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
    ordinal_outcome = "y_ord", patient_var = "pat_idx",
    random_slope = T, time_var = "time", ar_cov = F,
    burn_in = 250, nsim = 500, thin = 1, seed = 3,
    sd_vec = c(0.15, 0.30, rep(0.2, 4)))

# Test brmarm model against stored object ---------------------------------

  samps_ar <- baseline_bmr(
    formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
    ordinal_outcome = "y_ord", patient_var = "pat_idx", random_slope = T,
    time_var = "time", ar_cov = T, burn_in = 250, nsim = 500, thin = 1,
    seed = 3, sd_vec = c(0.15, 0.30, rep(0.2, 4)))

# Store results for later use or updating ---------------------------------

  #saveRDS(sim_data, ".\\Test_objects/bmrarm_data_object.RDS")
  #saveRDS(samps_int, ".\\Test_objects/bmr_int_object.RDS")
  #saveRDS(samps_slope, ".\\Test_objects/bmr_slope_object.RDS")
  #saveRDS(samps_ar, ".\\Test_objects/bmrarm_object.RDS")

  print(getwd())

  str <- ".\\..\\..\\Test_objects/"
  expect_equal_to_reference(sim_data, paste0(str, "bmrarm_data_object.RDS"))
  expect_equal_to_reference(samps_int, paste0(str, "bmr_int_object.RDS"))
  expect_equal_to_reference(samps_slope, paste0(str, "bmr_slope_object.RDS"))
  expect_equal_to_reference(samps_ar, paste0(str, "bmrarm_object.RDS"))
})
