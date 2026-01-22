#############################################################################
# New functions to compute tests for impact using IRLS representation of RDIF
############################################################################


# -------------------------------------------------------------------
#' The bi-square psi function.
#'
#' If \code{abs(u) > k} , \code{psi(u) = 0}. Else, \code{psi(u) = u(1 - (u/k)^2)^2}.
#'
#' @param u Can be a single value, vector, or matrix.
#' @param k The tuning parameter. Can be a scalar or the same dimension as \code{u}.
#' @return The bi-square psi function.
#' @seealso \code{\link[robustDIF]{rho}}
#' @export
# -------------------------------------------------------------------

psi <- function(u, k = 1.96) {
  if (length(k) != length(u)) {k <- k[1] + u - u}
  w <- (u / k)^2
  out <- u * (1 - w)^2
  out[abs(u) > k] <- 0
  out
}

# -------------------------------------------------------------------
#' The derivative of the bi-square psi function.
#'
#' If \code{abs(u) > k} , \code{psi_prime(u) = 0}. Else,
#'  \code{psi_prime(u) = (1 - (u/k)^2)^2 - (2u/k)^2 (1 - (u/k)^2)}.
#'
#' @param u Can be a single value, vector, or matrix.
#' @param k The tuning parameter. Can be a scalar or the same dimension as \code{u}.
#' @return The bi-square psi function.
#' @seealso \code{\link[robustDIF]{rho}}
#' @export
# -------------------------------------------------------------------

psi_prime <- function(u, k = 1.96) {
  if (length(k) != length(u)) {k <- k[1] + u - u}
  w <- (u / k)^2
  out <- (1 - w)^2 - 4 * w * (1 - w)
  out[abs(u) > k] <- 0
  out
}

# -------------------------------------------------------------------
#' Weights used to computed the variance of the RDIF scaling parameter estimated via IRLS
#'
#' These are not the IRLS weights, w, but the weights used for the variance for the IRLS estimate, v = psi_prime / sum w.
#'
#' If \code{abs(u) > k} , \code{psi_prime(u) = 0}. Else,
#'  \code{psi_prime(u) = (1 - (u/k)^2)^2 - (2u/k)^2 (1 - (u/k)^2)}.
#'
#' @inheritParams y_fun
#' @param alpha the desired false positive rate for flagging items with DIF.
#' @return A vector of weights v used to compute variance of RDIF scale parameter
#' @export
# -------------------------------------------------------------------


bsq_var_weight <- function(mle, theta, par = "intercept", log = F, alpha = .05){
  y <- y_fun(mle, par, log)
  var.y <- var_y(mle, theta, par, log)
  u <- (y - theta)/var.y
  omega <- (var.y - 1/sum(1/var.y))/var.y^2
  k <- qnorm(1 - alpha/2, 0, sqrt(omega))
  numer <- psi_prime(u, k) / var.y
  denom <- sum(bsq_weight(theta, y, var.y, alpha) / var.y)
  numer/denom
}

# -------------------------------------------------------------------
#' Maximum likelihood estimate (MLE) of scaling parameters
#'
#'
#' @inheritParams y_fun
#' @return The MLE of scaling parameter. Used internally to test DTF by comparing to RDIF estimate.
#' @export
# -------------------------------------------------------------------

ml_est <- function(mle, par = "intercept", log = F) {
  y <- y_fun(mle, par = par, log = log)
  var.y <- var_y(mle, theta = NULL, par = par, log = log)
  sum(y / var.y) / sum(1 / var.y)
}


# -------------------------------------------------------------------
#' Compares RDIF estimate of impact to maximum likelihood estimate.
#'
#' Asymptotic test of difference between RDIF estimate of impact and ML estimate of impact (based on IRT scaling functions)
#'
#'
#' @inheritParams y_fun
#' @param alpha the desired false positive rate for flagging items with DIF.
#' @return A named list with details of the test
#' @export
# -------------------------------------------------------------------

delta_test <- function(mle, par = "intercept", log = F, alpha = .05) {

  # Estimates
  rdif.theta <- rdif(mle, par, log, alpha = .05)$est
  mle.theta <- ml_est(mle, par, log)
  delta = rdif.theta - mle.theta

  # Var(delta)
  var.y.mle <- var_y(mle, theta = NULL, par, log)
  v.mle <- 1/var.y.mle / sum(1/var.y.mle)
  v.bsq <- bsq_var_weight(mle, rdif.theta, par, log, alpha)
  var.delta <- sum((v.bsq - v.mle)^2 * var.y.mle)

  # Test
  se.delta = sqrt(var.delta)
  z.test = delta/sqrt(var.delta)

  # Output
   c(rdif.est = rdif.theta,
     ml.est = mle.theta,
     delta = delta,
     se.delta = se.delta,
     z.test = z.test,
     p.val = (1 - pnorm(abs(z.test))) * 2)
}

