% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustDIF_functions.R
\name{rdif_chisq_test}
\alias{rdif_chisq_test}
\title{The R-DIF test of all parameters for each item.}
\usage{
rdif_chisq_test(mle, log = F)
}
\arguments{
\item{mle}{the output of \code{\link[robustDIF]{get_model_parms}}}

\item{log}{logical: use log of scaling function for the slopes?}
}
\value{
A data.frame containing the value of the chi2.test, its df, and p-value.
}
\description{
Simultaneously tests for DIF in all the item intercepts/thresholds and slopes using an asymptotic chi-square test.
}
\examples{
# Using the built-in example dataset "rdif.eg"
\dontrun{
rdif_chisq_test(mle = rdif.eg)
}

}
