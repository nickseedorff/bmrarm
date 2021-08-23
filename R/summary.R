#' Summary function for a bmrarm object
#'
#' @param x A bmrarm object
#' @param digits scalar, number of decimial places to round to
#' @export

summary.bmrarm <- function(x, digits = 3) {
  N_covar <- length(colnames(x$data$X))

  ## Beta coefficients
  df_b <- data.frame(
    outcome = rep(c("Ordinal", "Continuous"), each = N_covar),
    effect = rep(colnames(x$data$X), 2),
    round(summ(x$draws$res_beta), digits))

  cat("beta estimates")
  print(knitr::kable(df_b, align = "llccccc"))

  ## Covariance matrix
  df_s <- data.frame(
    parameter = c(paste0("Sigma_", apply(expand.grid(1:2, 1:2),
                                         1, paste0, collapse = "")), "rho"),
    round(summ(rbind(x$draws$res_sigma,
                     x$draws$res_ar)), digits))

  cat("\n R and rho estimates")
  print(knitr::kable(df_s, align = "lccccc"))

  ## Covariance matrix
  n_eff <- sqrt(nrow(x$draws$res_pat_sig))
  if(n_eff == 2) {
    descrip_vec <- c("Ord intercept", "cor", "Cont intercept", "cor")
  } else {
    descrip_vec <- c("Ord intercept", rep("cor", 4),
                     "Ord slope", rep("cor", 4),
                     "Cont intercept", rep("cor", 4),
                     "Cont slope")
  }
  df_pat <- data.frame(
    parameter = paste0("Sigma_alpha_", apply(expand.grid(1:n_eff, 1:n_eff), 1,
                                             paste0, collapse = "")),
    parameter_descrip = descrip_vec,
    round(summ(x$draws$res_sigma), digits))

  cat("\n Sigma_alpha estimates (random effects covariance matrix)")
  print(knitr::kable(df_pat, align = "llccccc"))

  ## Cut points
  tmp <- x$draws$res_cuts[, 1]
  locs <- which(!is.na(tmp) & !is.infinite(tmp))
  df_cut <- data.frame(
    parameter = paste0("Cutpoint_", locs - 1),
    round(summ(x$draws$res_cuts[locs,]), digits))

  cat("\n gamma (cutpoint) estimates")
  print(knitr::kable(df_cut, align = "lccccc"))

  list(beta = df_b, cov = df_s, rand_cov = df_pat, cuts = df_cut)
}

#' Summary function for a bmrvarx object
#'
#' @param x A bmrarm object
#' @param digits scalar, number of decimial places to round to
#' @export

summary.bmrvarx <- function(x, digits = 3) {
  N_outcomes <- ncol(x$y)
  N_covar <- length(x$covars_used)
  N_ord <- ifelse(is.na(dim(x$draws$res_cuts)[3]), 1, dim(x$draws$res_cuts)[3])

  ## Beta coefficients
  df_b <- data.frame(
    outcome = rep(names(x$y), each = N_covar),
    effect = rep(x$covars_used, N_outcomes),
    round(summ(x$draws$res_beta), digits))

  cat("beta estimates")
  print(knitr::kable(df_b, align = "llccccc"))

  ## M coefficients
  df_M <- data.frame(
    outcome = rep(names(x$y), N_outcomes),
    effect = paste0("Lagged ", rep(names(x$y), each = N_outcomes)),
    round(summ(x$draws$res_M), digits))

  cat("\n M estimates")
  print(knitr::kable(df_M, align = "llccccc"))

  ## Covariance matrix
  df_s <- data.frame(
    parameter = paste0("Sigma_", apply(expand.grid(1:N_outcomes, 1:N_outcomes),
                                       1, paste0, collapse = "")),
    round(summ(x$draws$res_sigma), digits))

  cat("\n Sigma estimates (between response covariance matrix)")
  print(knitr::kable(df_s, align = "lccccc"))

  ## Cut points
  if(N_ord == 1) {
    tmp <- x$draws$res_cuts[, 1]
    locs <- which(!is.na(tmp) & !is.infinite(tmp))
    df_cut <- data.frame(
      parameter = paste0("Cutpoint_", locs - 1),
      round(summ(x$draws$res_cuts[locs,]), digits))

    cat("\n gamma (cutpoint) estimates")
    print(knitr::kable(df_cut, align = "lccccc"))
  } else {
    df_cut <- list()
    for(i in 1:N_ord) {
      tmp <- x$draws$res_cuts[, 1, i]
      locs <- which(!is.na(tmp) & !is.infinite(tmp))
      df_tmp <- data.frame(
        parameter = paste0("Cutpoint_", locs - 1),
        round(summ(x$draws$res_cuts[locs,, i]), digits))
      df_cut[[i]] <- df_tmp

      cat("\n", names(x$y)[i], "gamma (cutpoint) estimates")
      print(knitr::kable(df_tmp, align = "lccccc"))
    }
  }
  list(beta = df_b, M = df_M, cov = df_s, cuts = df_cut)
}

#' PX-DA MCMC routine to sample from HBMRVAR model
#'
#' @param x A matrix of posterior draws
#' @importFrom coda mcmc effectiveSize

## Results for the regression coefficients
summ <- function(x) {
  res <- apply(x, 1, function(x) {
    c(mean = mean(x), sd =sd(x), lower_CI = quantile(x, probs = 0.025),
      upper_CI = quantile(x, probs = 0.975),
      eff_size = effectiveSize(mcmc(x)))
  })  %>%
    t()

  colnames(res)[3:5] <- c("lower_95_CI", "upper_95_CI", "effective_size")
  res
  }
