% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DTF_functions.R
\name{psi_prime}
\alias{psi_prime}
\title{The derivative of the bi-square psi function.}
\usage{
psi_prime(u, k = 1.96)
}
\arguments{
\item{u}{Can be a single value, vector, or matrix.}

\item{k}{The tuning parameter. Can be a scalar or the same dimension as \code{u}.}
}
\value{
The bi-square psi function.
}
\description{
If \code{abs(u) > k} , \code{psi_prime(u) = 0}. Else,
 \code{psi_prime(u) = (1 - (u/k)^2)^2 - (2u/k)^2 (1 - (u/k)^2)}.
}
\seealso{
\code{\link[robustDIF]{rho}}
}
