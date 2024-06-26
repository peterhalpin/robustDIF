% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustDIF_functions.R
\name{y_slope}
\alias{y_slope}
\title{The R-DIF scaling function for item slopes.}
\usage{
y_slope(mle, log = F)
}
\arguments{
\item{mle}{the output of \code{\link[robustDIF]{get_model_parms}}}

\item{log}{logical: use of scaling function?}
}
\value{
A vector of Y values.
}
\description{
Computes \code{Y = a1/a0}, or its log, for each item. Used to test DIF on item slopes.
}
