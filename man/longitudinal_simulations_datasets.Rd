#' Generated datasets for bmrarm simulations
#'
#' A list of datasets used in the simulated bmrarm examples.
#'
#' @format Two lists of data frames. One list contains data generated with an autoregressive term, the other without.
#' \describe{
#'   \item{y1}{true continuous values for the ordinal outcome, unused for modeled}
#'   \item{y_ord}{ordinal outcome, obtained from discretizing y1}
#'   \item{y2}{observed continuous outcome}
#'   \item{pat_idx}{patient identifier, integer}
#'   \item{time}{time since baseline, integer}
#'   \item{data_type}{identifies if data was used for training or testing. Only applicable to ar_data_list data frames as all data was training data for the slope_data_list data frames.}
#'   ...
#' }
"bmrarm_simulations"
