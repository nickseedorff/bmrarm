% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmrvarx_da.R
\name{bvar}
\alias{bvar}
\title{Gibbs sampler to implement a first order bvar model}
\usage{
bvar(
  formula,
  data,
  sig_prior = 1e+06,
  nsim = 1000,
  burn_in = 100,
  thin = 10,
  seed = 14
)
}
\arguments{
\item{formula}{an object of class "formula"; a symbolic description of the model to be fitted}

\item{data}{a dataframe containing outcome variables, covariates, and a patient or subject identifier}

\item{sig_prior}{scalar, prior variance on the regression coefficients}

\item{nsim}{positive integer, number of iterations with default of 1000}

\item{burn_in}{positive integer, number of iterations to remove with default of 100}

\item{thin}{positive integer, specifiers the period of saving samples. Default of 10}

\item{seed}{positive integer, seed for random number generation}
}
\value{
bmrvarx
}
\description{
Gibbs sampler to implement a first order bvar model
}
