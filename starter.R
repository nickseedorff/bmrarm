library(devtools)
load_all()
i <- 1
sim_data <- gen_ar_errors(N = 6, N_pat = 55, unequal = F, seed = 4,
                          slope = T, ar_cov = T)


crossprod(sim_data$alpha) / nrow(sim_data$alpha)


crossprod(MASS::mvrnorm(100, mu = c(0, 0, 0, 0) ,Sigma = sim_data$sig_alpha)) / 100


formula = cbind(y_ord, y2) ~ time; data = sim_data$data;
ordinal_outcome = "y_ord"; patient_var = "pat_idx";
random_slope = F; time_var = "time"; ar_cov = F;
burn_in = 250; nsim = 750; thin = 1; seed = 3;
verbose = TRUE; sig_prior = 1000000000; sd_vec = c(0.15, 0.30, 0.1)
sd_vec = c(0.12, 0.30, 0.1)

## Create storage
set.seed(seed)
bmrarm_start(env = environment())

## Starting values for ordinal outcome
y[1:N_obs, 1] <- 0.5 + res_cuts[z, 1]
y[is.infinite(y[, 1]), 1] <- -0.5
y[is.na(y[, 1]), 1] <- 0

## Starting values for continuous outcomes
df <- data.frame(patient = samp_info$pat_idx_long, y = as.numeric(y),
                 outcome = rep(1:N_outcomes, each = N_obs)) %>%
  group_by(patient, outcome) %>%
  mutate(y_interp = na.approx(y, na.rm = FALSE),
         y_interp = ifelse(!is.na(y_interp), y_interp,
                           ifelse(row_number() == n(),
                                  lag(y_interp), lead(y_interp))))
y <- matrix(df$y_interp, ncol = N_outcomes)
y[is.na(y)] <- 0
i <- 2
samp_info$burn_in <- burn_in
samp_info$max_iter <- 10000000

cur_draws$beta <- matrix(sim_data$beta, 2)
cur_draws$sigma <- sim_data$sigma
cur_draws$ar <- sim_data$ar
cur_draws$pat_sig <- sim_data$sig_alpha
y <- cbind(sim_data$truth$y1, sim_data$truth$y2)
