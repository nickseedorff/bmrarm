#' Generated datasets for bmrvarx simulations
#'
#' A list of data frames used in the simulated bmrvarx examples. Data frame 167 from bmrvarx_simulations is
#' the randomly selected dataset used to assess convergence.
#'
#' @format Two lists of data frames. One list contains data generated with an autoregressive term, the other without:
#' \describe{
#'   \item{obs_num}{Time since baseline}
#'   \item{y1}{true continuous values for the first ordinal outcome}
#'   \item{y2}{true continuous values for the second ordinal outcome}
#'   \item{y3}{observed continuous outcome}
#'   \item{y_ord}{first ordinal outcome. Obtained by discretizing y1}
#'   \item{y_ord2}{first ordinal outcome. Obtained by discretizing y2}
#'   \item{x1}{continous covariate, randomly generated from a standard normal}
#'   \item{y3_miss}{y3 with 5 percent of values set to missing}
#'   \item{y_ord_miss}{y_ord_miss with 5 percent of values set to missing}
#'   \item{y_ord2_miss}{y_ord2_miss with 5 percent of values set to missing}
#'   \item{data_type}{identifies if data was used for training or testing. The last ten observations we used to evaluate forecast accuracy.}
#' }
#' @examples
#'
#' \dontrun{
#' # Assess convergence by fitting four separate chains
#' burn <- 5000; sims <- 25000; mod_list <- list()
#'
#' for(i in 1:4) {
#'   mod_list[[i]] <- bmrvarx(formula = cbind(y_ord, y2) ~ time,
#'     data = bmrvarx_simulations$ar_data_list[[167]],
#'     ordinal_outcome = "y_ord", patient_var = "pat_idx",
#'     random_slope = T, time_var = "time", ar_cov = T,
#'     burn_in = burn, nsim = sims, thin = 5, seed = i,
#'     sd_vec = c(0.14, 0.14, 0.35, 0.1, 0.23, 0.09))
#' }
#'
#' # Fit a bmrvarx to each of the 400 AR datasets (\eqn{rho = 0.35})
#' for(i in 1:400){
#'   samps <- bmrvarx(formula = cbind(y_ord, y2) ~ time,
#'     data = bmrvarx_simulations$ar_data_list[[i]],
#'     ordinal_outcome = "y_ord", patient_var = "pat_idx",
#'     random_slope = T, time_var = "time", ar_cov = T,
#'     burn_in = burn, nsim = sims, thin = 5, seed = i,
#'     sd_vec = c(0.14, 0.15, 0.35, 0.1, 0.2, 0.09))
#'   mDIC <- get_DIC(samps)
#'   cDIC <- get_DIC(samps, marginal = FALSE)
#'   f_out <- paste0("./ar_model_ar_data", i, ".RDS")
#'   saveRDS(list(samps = samps, data = sim_data), file = f_out)
#' }
#' }
#'
#' bmrvarx_simulations
"bmrvarx_simulations"
