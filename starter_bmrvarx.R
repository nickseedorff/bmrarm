library(devtools)
load_all()

sim_data <- gen_bmrvarx(N = 600, seed = 10, N_param = 3)

## Single ordinal outcome
bmrvarx(formula = cbind(y2, y3, y_ord) ~ x1,
        data = sim_data$data, nsim = 500, burn_in = 100, seed = 1,
        ordinal_outcomes = c("y_ord"), thin = 1,
        max_iter_rej = 10000000)

## Ordinal outcome and binary outcome
bmrvarx(formula = cbind(y_ord2, y3, y_ord) ~ x1,
        data = sim_data$data, nsim = 500, burn_in = 100, seed = 1,
        ordinal_outcomes = c("y_ord", "y_ord2"), thin = 1,
        max_iter_rej = 10000000)

