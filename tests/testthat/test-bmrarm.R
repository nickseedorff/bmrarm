test_that("bmrarm works as anticipated", {
  sim_data <- gen_ar_errors(N = 7, N_pat = 48, unequal = F, seed = 1,
                            slope = F, ar_cov = F)

# Test brm_int model against stored object ---------------------------------

  samps_int <- bmrarm(
    formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
    ordinal_outcome = "y_ord", patient_var = "pat_idx",
    random_slope = F, time_var = "time", ar_cov = F,
    burn_in = 250, nsim = 500, thin = 1, seed = 3,
    sd_vec = c(0.15, 0.30, rep(0.2, 2)))

# Test brm_slope model against stored object ------------------------------

  samps_slope <- bmrarm(
    formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
    ordinal_outcome = "y_ord", patient_var = "pat_idx",
    random_slope = T, time_var = "time", ar_cov = F,
    burn_in = 250, nsim = 500, thin = 1, seed = 3,
    sd_vec = c(0.15, 0.30, rep(0.2, 4)))

# Test brmarm model against stored object ---------------------------------

  samps_ar <- bmrarm(
    formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
    ordinal_outcome = "y_ord", patient_var = "pat_idx", random_slope = T,
    time_var = "time", ar_cov = T, burn_in = 250, nsim = 500, thin = 1,
    seed = 3, sd_vec = c(0.15, 0.30, rep(0.2, 4)))

# Test brmarm model against stored object ---------------------------------
  sim_data_single_obs <- rbind(sim_data$data, c(0.5, 4, 0.75, 49, 0))
  samps_ar_single <- bmrarm(
    formula = cbind(y_ord, y2) ~ time, data = sim_data_single_obs,
    ordinal_outcome = "y_ord", patient_var = "pat_idx", random_slope = T,
    time_var = "time", ar_cov = T, burn_in = 5, nsim = 10, thin = 1,
    seed = 3, sd_vec = c(0.15, 0.30, rep(0.2, 4)))

  expect_equal(c(1.846025, 3.399137),
               rowMeans(samps_ar_single$draws$res_cuts)[4:5], 4)

# Store results for later use or updating ---------------------------------

   # saveRDS(sim_data, ".\\Test_objects/bmrarm_data_object.RDS")
   # saveRDS(samps_int, ".\\Test_objects/bmr_int_object.RDS")
   # saveRDS(samps_slope, ".\\Test_objects/bmr_slope_object.RDS")
   # saveRDS(samps_ar, ".\\Test_objects/bmrarm_object.RDS")

  str <- ".\\..\\..\\Test_objects/"
  expect_equal_to_reference(sim_data, paste0(str, "bmrarm_data_object.RDS"))
  expect_equal_to_reference(samps_int, paste0(str, "bmr_int_object.RDS"))
  expect_equal_to_reference(samps_slope, paste0(str, "bmr_slope_object.RDS"))
  expect_equal_to_reference(samps_ar, paste0(str, "bmrarm_object.RDS"))


# Test DIC values ---------------------------------------------------------

  set.seed(1)
  expect_equal(round(get_DIC(samps_int, marginal = FALSE), 2),
               c(DIC = 1092.98, D = 1015.99, pd = 76.99,  dev_of_means = 939.00))
  expect_equal(round(get_DIC(samps_slope, marginal = FALSE), 2),
               c(DIC = 1100.97, D = 991.46, pd = 109.51,  dev_of_means = 881.96))
  expect_equal(round(get_DIC(samps_ar, marginal = FALSE), 2),
               c(DIC = 1105.06, D = 994.30, pd = 110.76,  dev_of_means = 883.54))
  expect_equal(round(get_DIC(samps_ar), 2),
               c(DIC = 1178.25, D = 1161.74, pd = 16.52,  dev_of_means = 1145.22))

})
