library(devtools)
load_all()
i <- 1
sim_data <- gen_ar_errors(N = 7, N_pat = 75, unequal = F, seed = i,
                          slope = T, ar_cov = T)

data <- sim_data$data

check <- rbind(select(sim_data$data, pat_idx, time, value = y_ord) %>%
                 mutate(type = 1),
               select(sim_data$data, pat_idx, time, value = y2) %>%
                 mutate(type = 2)) %>%
  filter(!is.na(value)) %>%
  cbind(pt_val = obs$pointwise[, 1]) %>%
  group_by(type, pat_idx) %>%
  summarise(sum(pt_val))

samps <- baseline_bmr(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                      ordinal_outcome = "y_ord", patient_var = "pat_idx",
                      random_slope = T, time_var = "time", ar_cov = T,
                      burn_in = 200, nsim = 600, thin = 5, seed = 3,
                      sd_vec = c(0.15, 0.15))
(obs <- get_waic_ar(samps))

samps <- bmr_cv(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                      ordinal_outcome = "y_ord", patient_var = "pat_idx",
                      random_slope = T, time_var = "time", ar_cov = F,
                      burn_in = 200, nsim = 800, thin = 5, seed = 3,
                      sd_vec = c(0.15, 0.15))

samps <- bmr_cv(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                ordinal_outcome = "y_ord", patient_var = "pat_idx",
                random_slope = T, time_var = "time", ar_cov = T,
                burn_in = 200, nsim = 800, thin = 5, seed = 3,
                sd_vec = c(0.15, 0.15))

samps <- bmr_cv(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                ordinal_outcome = "y_ord", patient_var = "pat_idx",
                random_slope = F, time_var = "time", ar_cov = F,
                burn_in = 200, nsim = 800, thin = 5, seed = 3,
                sd_vec = c(0.15, 0.15))



(obs <- get_waic_ar(samps))
get_DIC_ar(samps)
get_DIC(samps)
pareto_k_ids(obs, threshold = 0.7)
plot(samps$res_ar)

samps <- baseline_bmr(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                       ordinal_outcome = "y_ord", patient_var = "pat_idx",
                       random_slope = T, time_var = "time", ar_cov = T,
                       burn_in = 500, nsim = 2000, thin = 5, seed = 3)

get_DIC_ar(samps)

samps2 <- baseline_bmr(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                      ordinal_outcome = "y_ord", patient_var = "pat_idx",
                      random_slope = T, time_var = "time", ar_cov = F,
                      burn_in = 500, nsim = 2000, thin = 5, seed = 3)

samps2$res_ar[] <- 0
(obs2 <- get_waic_ar(samps2))
get_DIC(samps2)
obs2
(pareto_k_ids(obs2, threshold = 0.7) - (6 * 65)) %% 6
get_DIC(samps)

rowMeans(samps$res_pat_sig)[c(1, 6, 11 , 16)]
rowMeans(samps2$res_pat_sig)[c(1, 6, 11 , 16)]

rowMeans(samps$res_sigma)
rowMeans(samps2$res_sigma)

rowMeans(samps$res_cuts)
rowMeans(samps2$res_cuts)

rowMeans(samps$res_beta)
rowMeans(samps2$res_beta)

hist(samps$res_beta[1, ])
hist(samps2$res_beta[1, ])

sim_data$sigma

samps3 <- baseline_bmr(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                       ordinal_outcome = "y_ord", patient_var = "pat_idx",
                       random_slope = F, time_var = "time", ar_cov = F,
                       burn_in = 400, nsim = 2000, thin = 5, seed = 3)
get_DIC(samps3)

samps3$res_ar[] <- 0
(obs3 <- get_waic_ar(samps3))

crossprod(apply(samps$res_pat_eff, c(1, 2), mean) - sim_data$alpha)
cor(apply(samps$res_pat_eff, c(1, 2), mean), sim_data$alpha)


samps2 <- baseline_bmr_test(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                           ordinal_outcome = "y_ord", patient_var = "pat_idx",
                           random_slope = T, time_var = "time", ar_cov = T,
                           burn_in = 300, nsim = 800, thin = 5, seed = 3, which_fc = c(T, F, T, F, F),
                           base0 = T)
get_DIC(samps2)
crossprod(apply(samps2$res_pat_eff, c(1, 2), mean) - sim_data$alpha)
cor(apply(samps2$res_pat_eff, c(1, 2), mean), sim_data$alpha)


rowMeans(samps$res_beta)
rowMeans(samps2$res_beta)
matrix(rowMeans(samps$res_pat_sig), ncol = 4)
matrix(rowMeans(samps2$res_pat_sig), ncol = 4)


diag(crossprod(apply(samps$res_pat_eff, c(1, 2), mean) - sim_data$alpha)) /
  diag(crossprod(apply(samps2$res_pat_eff, c(1, 2), mean) - sim_data$alpha))

get_DIC_ar(samps2)

get_DIC(samps)

get_DIC(samps)

tail(apply(samps$res_pat_eff, c(1, 2), mean), 6)
tail(apply(samps2$res_pat_eff, c(1, 2), mean), 6)

tail(apply(samps$res_pat_eff, c(1, 2), mean), 6)
tail(apply(samps2$res_pat_eff, c(1, 2), mean), 6)
tail(sim_data$alpha, 6)

samps3 <- baseline_bmr(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                       ordinal_outcome = "y_ord", patient_var = "pat_idx",
                       random_slope = T, time_var = "time", ar_cov = F,
                       burn_in = 300, nsim = 800, thin = 20, seed = 3)

get_DIC(samps3)
samps2$res_ar[] <- 0.25
get_DIC_ar(samps2)
matrix(rowMeans(samps2$res_pat_sig), 4)
matrix(rowMeans(samps$res_pat_sig), 4)
rowMeans(samps$res_beta)
rowMeans(samps2$res_beta)
mean(samps$res_ar)

samps3 <- baseline_bmr(formula = cbind(y_ord, y2) ~ time, data = sim_data$data,
                       ordinal_outcome = "y_ord", patient_var = "pat_idx",
                       random_slope = F, time_var = "time", ar_cov = F,
                       burn_in = 300, nsim = 600, thin = 1, seed = 3)

get_DIC(samps3)
