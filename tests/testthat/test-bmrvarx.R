test_that("bmrarm works as anticipated", {
  sim_data <- gen_single(N = 400, seed = 10, N_param = 3)

# Test bmrvarx with one ordinal outcome ---------------------------------

  samps_one <- bmrvarx(formula = cbind(y2, y3, y_ord) ~ x1,
                       data = sim_data$data, nsim = 400, burn_in = 200,
                       seed = 1,
                       ordinal_outcomes = c("y_ord"), thin = 1,
                       max_iter_rej = 10000000)


# Test brm_slope model against stored object ------------------------------

  samps_two <- bmrvarx(formula = cbind(y3, y_ord, y_ord2) ~ x1,
                       data = sim_data$data, nsim = 400, burn_in = 200,
                       seed = 1,
                       ordinal_outcomes = c("y_ord", "y_ord2"), thin = 1,
                       max_iter_rej = 10000000)

# Store results for later use or updating ---------------------------------

  # saveRDS(sim_data, ".\\Test_objects/bmrvarx_data_object.RDS")
  # saveRDS(samps_one, ".\\Test_objects/bmrvarx_one_ord.RDS")
  # saveRDS(samps_two, ".\\Test_objects/bmrvarx_two_ord.RDS")

  str <- ".\\..\\..\\Test_objects/"
  expect_equal_to_reference(sim_data, paste0(str, "bmrvarx_data_object.RDS"))
  expect_equal_to_reference(samps_one, paste0(str, "bmrvarx_one_ord.RDS"))
  expect_equal_to_reference(samps_two, paste0(str, "bmrvarx_two_ord.RDS"))
})
