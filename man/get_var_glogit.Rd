% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variance.r
\name{get_var_glogit}
\alias{get_var_glogit}
\title{Get GEE glogit Co-Variance Matrix}
\usage{
get_var_glogit(betas, rho, dat)
}
\arguments{
\item{betas}{vector of betas}

\item{rho}{overdispersion parameter}

\item{dat}{matrix of multinomial clusters}
}
\value{
matrix of covariance
}
\description{
Calculates GEE glogit variance co-variance matrix.
}
