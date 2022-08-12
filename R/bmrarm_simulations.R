#' Generated datasets for bmrarm simulations
#'
#' A list, which contains two lists of data frames that are used in the simulated bmrarm examples.
#' Data frame 167 from bmrarm_simulations$ar_data_list is the randomly selected dataset used to assess convergence.
#'
#' @format Two lists of data frames. One list contains data generated with an autoregressive term, the other without:
#' \describe{
#'   \item{ar_data_list}{list of 400 dataframes generated using bmrarm with \eqn{rho = 0.35}}
#'   \item{slope_data_list}{list of 400 dataframes generated using bmrarm with \eqn{rho = 0}}
#'   \item{y1}{true continuous values for the ordinal outcome, unused for modeling}
#'   \item{y_ord}{ordinal outcome, obtained from discretizing y1}
#'   \item{y2}{observed continuous outcome}
#'   \item{pat_idx}{patient identifier, integer}
#'   \item{time}{time since baseline, integer}
#'   \item{data_type}{identifies if data was used for training or testing. Only applicable to ar_data_list data frames as all data was training data for the slope_data_list data frames.}
#' }
#' @examples
#'
#' \dontrun{
#' # Assess convergence by fitting four separate chains
#' burn <- 5000; sims <- 25000; mod_list <- list()
#'
#' for(i in 1:4) {
#'   mod_list[[i]] <- bmrarm(formula = cbind(y_ord, y2) ~ time,
#'     data = bmrarm_simulations$ar_data_list[[167]],
#'     ordinal_outcome = "y_ord", patient_var = "pat_idx",
#'     random_slope = T, time_var = "time", ar_cov = T,
#'     burn_in = burn, nsim = sims, thin = 5, seed = i,
#'     sd_vec = c(0.14, 0.14, 0.35, 0.1, 0.23, 0.09))
#' }
#'
#' # Fit a bmrarm to each of the 400 AR datasets (\eqn{rho = 0.35})
#' for(i in 1:400){
#'   samps <- bmrarm(formula = cbind(y_ord, y2) ~ time,
#'     data = bmrarm_simulations$ar_data_list[[i]],
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
#' bmrarm_simulations
"bmrarm_simulations"
