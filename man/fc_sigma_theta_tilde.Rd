% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_bmrvarx.R
\name{fc_sigma_theta_tilde}
\alias{fc_sigma_theta_tilde}
\title{Full conditional draws of the regression coefficients}
\usage{
fc_sigma_theta_tilde(
  y,
  X,
  prior_precision,
  y_orig,
  pat_idx = NULL,
  scale_first = FALSE
)
}
\arguments{
\item{y}{matrix of multivariate observations}

\item{X}{design matrix}

\item{prior_precision}{prior precision matrix}
}
\value{
matrix
}
\description{
Full conditional draws of the regression coefficients
}