% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_bmrvarx.R
\name{fc_y}
\alias{fc_y}
\title{Function to draw latent continuous values}
\usage{
fc_y(y, z, mean_mat, tmp_list, miss_mat, samp_info, rej_vec)
}
\arguments{
\item{y}{matrix of continuous observations}

\item{z}{matrix of ordinal outcomes}

\item{tmp_list}{list of current parameter values}

\item{miss_mat}{matrix of location of missing values}

\item{samp_info}{list of internal information used for sampling}

\item{rej_vec}{vector, values track the attempts by the rejection sampler}
}
\value{
list
}
\description{
Function to draw latent continuous values
}
