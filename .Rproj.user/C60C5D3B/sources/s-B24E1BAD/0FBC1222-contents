library(devtools)
load_all()

#sim_data2 <- gen_bmrvarx(N = 600, seed = 10, N_param = 3)
sim_data <- gen_single(N = 600, seed = 10, N_param = 3)
tail(sim_data$data)

dat <- sim_data$data
rm(sim_data)
formula = cbind(y3, y_ord, y_ord2) ~ x1;
data = dat;nsim = 500; burn_in = 100; seed = 1;
ordinal_outcomes = c("y_ord", "y_ord2");thin = 1;
max_iter_rej = 10000000; sig_prior = 1000000;
i <- 1

y <- cbind(dat$y1, dat$y2, dat$y3)
tmp_list$beta <- sim_data$beta
tmp_list$M <- sim_data$M
tmp_list$sigma <- sim_data$sigma
tmp_list$cuts <- sim_data$cuts
z <- dat$y_ord

## Single ordinal outcome
f <- bmrvarx(formula = cbind(y3, y_ord, y_ord2) ~ x1,
        data = sim_data$data, nsim = 400, burn_in = 200, seed = 1,
        ordinal_outcomes = c("y_ord", "y_ord2"), thin = 1,
        max_iter_rej = 10000000, fast = T)

plot(f$draws$res_cuts[3,, 1])
plot(f$draws$res_cuts[3,, 2])


f2 <- bmrvarx_da(formula = cbind(y3, y_ord, y_ord2) ~ x1,
             data = data, nsim = 500, burn_in = 200, seed = 1,
             ordinal_outcomes = c("y_ord", "y_ord2"), thin = 1,
             max_iter_rej = 10000000, fast = T)

f2 <- bmrvarx(formula = cbind(y2, y3, y_ord) ~ x1,
             data = data, nsim = 1000, burn_in = 200, seed = 1,
             ordinal_outcomes = c("y_ord"), thin = 1,
             max_iter_rej = 10000000, fast = T)


f2 <- bvar(formula = cbind(y1, y2, y3) ~ x1,
              data = filter(data, !is.na(y2), !is.na(y3)), nsim = 10000,
              burn_in = 2000, seed = 1, thin = 1)

rowMeans(f2$draws$res_M)
rowMeans(f2$draws$res_beta)
rowMeans(f2$draws$res_sigma)


rowMeans(f2$draws$res_sigma)

rowMeans(f$draws$res_cuts)
rowMeans(f2$draws$res_cuts)

rowMeans(f$draws$res_beta)
rowMeans(f2$draws$res_beta)


## Ordinal outcome and binary outcome
bmrvarx(formula = cbind(y_ord2, y3, y_ord) ~ x1,
        data = sim_data$data, nsim = 500, burn_in = 100, seed = 1,
        ordinal_outcomes = c("y_ord", "y_ord2"), thin = 1,
        max_iter_rej = 10000000)

