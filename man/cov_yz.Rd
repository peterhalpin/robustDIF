% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustDIF_functions.R
\name{cov_yz}
\alias{cov_yz}
\title{Compute the covariance matrix for \code{\link[robustDIF]{rdif_chisq_test}}.}
\usage{
cov_yz(theta.y, theta.z, mle, log = F)
}
\arguments{
\item{theta.y}{the IRT scale parameter for the item intercepts}

\item{theta.z}{the IRT scale parameter for the item slopes}

\item{mle}{the output of \code{\link[robustDIF]{get_model_parms}}}

\item{log}{logical: use log of scaling function for the slopes?}
}
\value{
A block diagonal matrix whose blocks contain the covariance matrix of the IRT scaling functions, for each item.
}
\description{
Compute the covariance matrix for \code{\link[robustDIF]{rdif_chisq_test}}.
}
\seealso{
\code{\link[robustDIF]{y_fun}}, \code{\link[robustDIF]{var_y}}
}
