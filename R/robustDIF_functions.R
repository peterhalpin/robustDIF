#' An example data set with five items in two groups.
#'
#' A named list containing the maximum likelihood estimates and their estimated covariance matrix, for the 2PL IRT model fitted to 5 items in two independent groups. The first item has additive bias of .5 applied to the intercept. The groups have a mean difference of .5 standard deviations on the latent trait. The variances of the latent trait are equal in each group.
#'
#' @format A named list with 3 components:
#' \describe{
#'   \item{par.names}{A named \code{list} with the internal and original labels for parameters in the model.}
#'   \item{est}{A named \code{list}, each element containing \code{data.frame} of model parameters estimates for each group}
#'   \item{vcov0}{A named \code{list}, each element containing \code{data.frame} with the variance covariance matrix of the parameter estimates. The parameters appear in the order given by rdif.eg$par.names$internal.
#'}
#'}
"rdif.eg"

# -------------------------------------------------------------------
#' The R-DIF scaling function for item intercepts / thresholds
#'
#' Computes the scaling function \code{Y = (d1 - d0)/a1} for each intercept / threshold of each item. Used to test DIF on item intercepts.
#'
#' @param mle the output of \code{\link[robustDIF]{get_model_parms}}
#' @return A vector of Y values.
# -------------------------------------------------------------------

y_intercept <- function(mle) {
  c(t((mle$est$group.2[, -1] - mle$est$group.1[, -1]) /  mle$est$group.2$a1))
}


# -------------------------------------------------------------------
#' The R-DIF scaling function for item slopes.
#'
#' Computes \code{Y = a1/a0}, or its log, for each item. Used to test DIF on item slopes.
#'
#' @inheritParams y_intercept
#' @param log logical: use of scaling function?
#'
#' @return A vector of Y values.
# -------------------------------------------------------------------

y_slope <- function(mle, log = F) {
  y <- mle$est$group.2$a1  / mle$est$group.1$a1
  if (log) {
    y <- log(y)
  }
  y
}


# -------------------------------------------------------------------
#' R-DIF scaling functions.
#'
#' Computes the R-DIF scaling function for item intercepts/thresholds or for item slopes.
#'
#' @inheritParams y_intercept
#' @param par character: use the scaling function for item intercepts or item slopes? One of \code{c("intercept", "slope")}.
#' @param log logical: use log of scaling function? Only applies if \code{par = "slope"}.
#'
#' @return A vector of Y values.
#' @seealso \code{\link[y_intercept]{y_intercept}}, \code{\link[robustDIF]{y_slope}}.
# -------------------------------------------------------------------

y_fun <- function(mle, par = "intercept", log = F) {
  if (par == "slope") {
    y <- y_slope(mle, log)
  } else {
    y <- y_intercept(mle)
  }
  y
}

# -------------------------------------------------------------------
#' The gradient of \code{\link[robustDIF]{y_intercept}}.
#'
#' The gradient is take with respect to the item parameters. The parameter vector is organized as \code{c(a0, d0, a1, d1)}, repeated over items. The parameters \code{a0} and \code{d0} are the item slope and intercept in the reference group. The parameters \code{a1} and \code{d1} are the item slope and intercept in the comparison group.
#'
#' @inheritParams y_intercept
#' @param theta the IRT scale parameter. If omitted, uses item-specific scaling functions instead.
#' @return A matrix in which the columns are the gradient vectors of \code{\link[robustDIF]{y_intercept}}, for each item.
#'
#' @seealso \code{\link[y_intercept]{y_intercept}}
# -------------------------------------------------------------------

grad_intercept <- function(mle, theta = NULL) {
  n.items <- nrow(mle$est$group.1)
  n.item.pars <- ncol(mle$est$group.1)

  if (is.null(theta)) {
    theta <- y_intercept(mle)
  }
  if (length(theta) == 1) {
    theta <- rep(theta, times = n.items)
  }

  # Make template matrix for each item (must be easier way!)
  template <- matrix(0, nrow = n.item.pars, ncol = n.item.pars - 1)
  for (j in 1:(n.item.pars - 1)) {
    template[1, j] <- template[(j + 1), j] <- 1
  }
  template <- rbind(template, template) # generalize for > 2 groups
  template[1, ] <- 0
  template[1:n.item.pars, ] <- -1 * template[1:n.item.pars, ]

  # Divide all entries by a1
  grad.list <- lapply(mle$est$group.2$a1, function (x)
    template / x)

  # Multiply n.item.pars+1 entry by -theta
  for (i in 1:n.items) {
    grad.list[[i]][(n.item.pars + 1), ] <-
      -1 * theta[i] * grad.list[[i]][(n.item.pars + 1), ]
  }

  Matrix::bdiag(grad.list)
}

# -------------------------------------------------------------------
#' The gradient of \code{\link[robustDIF]{y_slope}}.
#'
#' The gradient is take with respect to the item parameters. The parameter vector is organized as \code{c(a0, d0, a1, d1)}, repeated over items. The parameters \code{a0} and \code{d0} are the item slope and intercept in the reference group. The parameters \code{a1} and \code{d1} are the item slope and intercept in the comparison group.
#'
#' @inheritParams grad_intercept
#' @param log logical: use log of scaling function?
#' @return A matrix in which the columns are the gradient vectors of \code{\link[robustDIF]{y_slope}}, for each item.
#'
#' @seealso \code{\link[y_intercept]{y_slope}}
# -------------------------------------------------------------------

grad_slope <- function(mle, theta = NULL, log = F) {
  n.items <- nrow(mle$est$group.1)
  n.item.pars <- ncol(mle$est$group.1)
  n.groups <- length(mle$est)

  if (is.null(theta)) {
    theta <- y_slope(mle, log)
  }
  if (length(theta) == 1) {
    theta <- rep(theta, times = n.items)
  }

  a <- mle$est$group.1$a1
  if (log) {
    theta <- exp(theta)
    a <- mle$est$group.2$a1
  }

  # Make template matrix for each item gradients
  template <- rep(0, times = n.item.pars * n.groups)
  template[1] <- -1 # -theta
  template[n.item.pars + 1] <- 1

  # Divide all entries by a
  grad.list <- lapply(a, function (x)
    template / x)

  # Multiply first entry by theta
  for (i in 1:n.items) {
    grad.list[[i]][1] <- theta[i] * grad.list[[i]][1]
  }
  Matrix::bdiag(grad.list)
}

# -------------------------------------------------------------------
#' Helper function to merge VCOV matrices from different groups.
#'
#' Puts VCOV matrix of two models into a single block diagonal VCOV matrix that is conformable with the output of the \code{\link[robustDIF]{grad_intercept}} and \code{\link[robustDIF]{grad_slope}}.
#'
#' @inheritParams y_intercept
#'
#' @importFrom Matrix bdiag
#' @return A block diagonal VCOV matrix in which each block is the VCOV of an item.

# -------------------------------------------------------------------

joint_vcov <- function(mle) {
  n.items <- nrow(mle$est$group.1)
  n.item.pars <- ncol(mle$est$group.1)
  n.groups <- length(mle$est)

  # Create blocks for each item
  mat.template <-
    diag(1:n.items) %x% matrix(1, n.item.pars, n.item.pars)

  # Extract vcov for each item, ordered items within groups
  temp1 <- lapply(mle$var.cov,
                  function(x)
                    split(as.matrix(x), mat.template)[-1])

   # Re-order to groups within items
  temp2 <- lapply(1:n.items,
                  function(x)
                    lapply(temp1, function(y)
                      y[[x]]))

  # Flatten list and reformat to matrices
  vcov.list <-
    lapply(Reduce(c, temp2), function(x)
      matrix(x,  n.item.pars))

  # Convert to block diagonal
  Matrix::bdiag(vcov.list)
}

# -------------------------------------------------------------------
#' Compute the variance of the asymptotic null distribution of IRT scaling functions.
#'
#' @param theta the IRT scale parameter.
#' @inheritParams y_fun
#'
#' @return A vector that contains the variance of the IRT scaling function, for each item.
#'
#' @seealso \code{\link[robustDIF]{y_fun}}
#'
#' @importFrom Matrix diag
# -------------------------------------------------------------------

var_y <- function(mle,
                  theta = NULL,
                  par = "intercept",
                  log = F) {
  grad.y <- grad_intercept(mle, theta)
  if (par == "slope") {
    grad.y <- grad_slope(mle, theta, log)
  }
  joint.vcov <- joint_vcov(mle)
  vcov.y <- Matrix::t(grad.y) %*% joint.vcov %*% grad.y
  Matrix::diag(vcov.y)
}


# -------------------------------------------------------------------
#' Compute the covariance matrix for \code{\link[robustDIF]{rdif_chisq_test}}.
#'
#'
#' @param theta.y the IRT scale parameter for the item intercepts
#' @param theta.z the IRT scale parameter for the item slopes
#' @inheritParams y_fun
#' @param log logical: use log of scaling function for the slopes?

#' @return A block diagonal matrix whose blocks contain the covariance matrix of the IRT scaling functions, for each item.
#' @seealso \code{\link[robustDIF]{y_fun}}, \code{\link[robustDIF]{var_y}}
#'
#' @importFrom Matrix diag
# -------------------------------------------------------------------

# vcov_yz2 <- function(theta.y, theta.z, mle, log = F) {
#   n.items <- nrow(mle$est$group.1)
#   n.item.pars <- ncol(mle$est$group.1)
#   grad.y <- grad_intercept(theta.y, mle)
#   grad.z <- grad_slope(theta.z, mle, log = log)
#   grad <- cbind(grad.y, grad.z)*0
#
#   # Re-order columns by items (intercepts, then slope)
#   for(i in 1:nrow(mle$est$group.1)){
#     m <- (i-1)*n.item.pars + 1
#     n <- m + n.item.pars - 2
#     r <- (i-1)*(n.item.pars-1) + 1
#     s <- r + n.item.pars - 2
#     grad[, m:n] <- grad.y[, r:s]
#     grad[, (n+1)] <- grad.z[, i]
#   }
#   Matrix::t(grad) %*% joint_vcov(mle) %*% grad
# }

cov_yz <- function(theta.y, theta.z, mle, log = F) {
  n.items <- nrow(mle$est$group.1)
  n.item.pars <- ncol(mle$est$group.1)

  var.y <- var_y(mle, theta.y, par = "intercept")
  w <- (1 / var.y) / sum(1 / var.y)
  grad.y <- grad_intercept(mle, theta.y)
  grad.theta <- grad.y %*% matrix(w)
  grad.y.theta <-
    grad.y - outer(grad.theta, rep(1, times = ncol(grad.y)))

  var.z <- var_y(mle, theta.y, par = "slope", log = log)
  v <- (1 / var.z) / sum(1 / var.z)
  grad.z <- grad_slope(mle, theta.y, log)
  grad.sigma <- grad.z %*% matrix(v)
  grad.z.sigma <-
    grad.z - outer(grad.sigma, rep(1, times = ncol(grad.z)))


  # Re-order columns by items (intercepts, then slope)
  grad <- cbind(grad.y, grad.z) * 0
  for (i in 1:nrow(mle$est$group.1)) {
    m <- (i - 1) * n.item.pars + 1
    n <- m + n.item.pars - 2
    r <- (i - 1) * (n.item.pars - 1) + 1
    s <- r + n.item.pars - 2
    grad[, m:n] <- grad.y.theta[, r:s]
    grad[, (n + 1)] <- grad.z.sigma[, i]
  }

  vcov <- joint_vcov(mle)
  vcov.yz <- Matrix::t(grad) %*% vcov %*% grad

  # Convert to block diag so can inverted and compute chsiq for all items at once
  # (Otherwise, would be omnibus test of MI over all items)
  ones <- matrix(1, nrow = n.item.pars, ncol = n.item.pars)
  ones.list <- lapply(vector("list", n.items), function(x)
    x = ones)
  vcov.yz * Matrix::bdiag(ones.list)
}

# -------------------------------------------------------------------
#' Compute the weights of the bi-square function.
#'
#' The weights are used for estimation via iteratively re-weighted least squares, so the function is set up to take a pre-computed values of \code{theta} and \code{var.y}
#'
#' @param theta the IRT scale parameter.
#' @param y the output of \code{\link[robustDIF]{y_fun}}.
#' @param var.y the output of \code{\link[robustDIF]{var_y}}.
#' @param alpha the desired false positive rate for flagging items with DIF.
#' @return The bi-square weights, for each item.
# -------------------------------------------------------------------

bsq_weight <- function(theta, y, var.y, alpha = .05) {
  r <- y - theta
  var.theta <- 1 / sum(1 / var.y)
  omega <- (var.y - var.theta) / var.y ^ 2
  k <- qnorm(1 - alpha / 2, 0, sqrt(omega))
  bsq.cuts <- var.y * k
  w <- (1 - (r / bsq.cuts) ^ 2) ^ 2
  w[abs(r) > bsq.cuts] <- 0
  w
}

# -------------------------------------------------------------------
#' Compute staring values for \code{\link[robustDIF]{rdif}}.
#'
#' @inheritParams y_fun
#' @param alpha the desired false positive rate for flagging items with DIF.

#' @return A vector containing the median of \code{\link[robustDIF]{y_fun}}, the least trimmed squares estimate of location for \code{\link[robustDIF]{y_fun}} with 50-percent trim rate, and the minimum of \code{\link[robustDIF]{rho_fun}}.
# -------------------------------------------------------------------

get_starts <-
  function(mle,
           par = "intercept",
           log = F,
           alpha = .05) {
    y <- y_fun(mle, par, log = log)
    var_fun <- function(theta) {
      var_y(mle, theta, par, log = log)
    }

    #This combines both choices of reference group for intercepts
    # Was not re-written for new parms
    # if (par == "intercept") {
    #   y1 <- y
    #   mle2 <- mle
    #   names(mle2)[1:2] <- names(mle)[2:1]
    #   y2 <- -1 * y_fun(mle2)
    #
    #   s <- median(y1/y2)
    #   y2 <- s * y2
    #
    #   # Drop items if y1 - y2 / sd(y) > 1.5
    #   var.y1 <- var_y(median(y1), mle)
    #   var.y2 <- s^2 * var_y(median(y2), mle)
    #   drops <- abs((y1 - y2) / sqrt(min(var.y1, var.y2))) < 1.5
    #   y <- c(y1[drops], y2[drops])
    #   var_fun <- function(theta) {
    #     c(var_y(theta, mle)[drops], s^2 * var_y(theta, mle2)[drops])
    #   }
    # }

    # median residual
    s1 <- median(y)

    # lts with 50% of data
    s2 <- lts(y)

    # grid search for min of Rho function
    s3 <- rho_grid(y, var_fun, alpha = .05)

    c(s1, s2, s3)
  }

# -------------------------------------------------------------------
#' The least trimmed squares (LTS) estimate of location
#'
#' @param y a vector of data points.
#' @param p the proportion of data points to trim.
#' @return The LTS estimate of location of \code{y}.
#'
#' @seealso \code{\link[robustDIF]{get_starts}}
# -------------------------------------------------------------------

lts <- function(y, p = .5) {
  n <- length(y)
  n.per.sample <- floor(n * p)
  n.sample <- ceiling(n - n.per.sample + 1)
  Y <- matrix(sort(y), nrow = n, ncol = n.sample)
  flag <- c(rep(1, times = n.per.sample),
            rep(0, times = n.sample))
  Flag <- matrix(rep(flag, times = n.sample)[1:(n * n.sample)],
                 nrow = n,
                 ncol = n.sample)
  dat <- Y * Flag
  sum(dat[, which.min(apply(dat, 2, var))]) / n.per.sample
}

# -------------------------------------------------------------------
#' The bi-square rho function.
#'
#' If \code{abs(u) > k} , \code{rho(u) = 1}. Else, \code{rho(u) = (1 - (1 - (u/k))^3)}.
#'
#' @param u Can be a single value, vector, or matrix.
#' @param k The tuning parameter Can be a scalar or the same dimension as \code{u}.
#' @return The bi-square rho function.
#' @seealso \code{\link[robustDIF]{rho_fun}}
# -------------------------------------------------------------------

rho <- function(u, k = 1.96) {
  if (length(k) != length(u)) {
    k <- k[1] + u - u
  }
  w <- (u / k) ^ 2
  out <- 1 - (1 - w) ^ 3
  out[abs(u) > k] <- 1
  out
}

# -------------------------------------------------------------------
#' Grid search for minimum value of the bi-square Rho function.
#'
#' This function is used by \code{\link[robustDIF]{get_starts}}, so its inputs are formatted for to take internal arguments of that function. For plotting, use \code{\link[robustDIF]{rho_fun}} instead.
#'
#' @param y output of \code{\link[robustDIF]{get_starts}}
#' @param var_fun function that variance of \code{y} as a function of only the scale parameter.
#' @param alpha The desired false positive rate for flagging items with DIF.
#' @param grid.width The width of grid points.
#'
#' @return The location parameter that minimizes the bi-square Rho function.
#' @seealso \code{\link[robustDIF]{get_starts}}

# -------------------------------------------------------------------

rho_grid <- function(y,
                     var_fun,
                     alpha = .05,
                     grid.width = .05) {
  theta <- seq(
    from = max(min(y), -1.5),
    to = min(max(y), 1.5),
    by = min(grid.width, (max(y) - min(y)) / 10)
  )

  n.items <- length(y)
  n.theta <- length(theta)

  Y <- matrix(y, nrow = n.items, ncol = n.theta)
  Var.Y <- Reduce(cbind, lapply(theta, function(x)
    var_fun(x)))
  Theta <- matrix(theta,
                  nrow = n.items,
                  ncol = n.theta,
                  byrow = T)
  Var.Theta <- matrix(
    1 / apply(1 / Var.Y, 2, sum),
    nrow = n.items,
    ncol = n.theta,
    byrow = T
  )
  Omega <- (Var.Y - Var.Theta) / Var.Y ^ 2
  U <- (Y - Theta) / Var.Y
  K <-  qnorm(1 - alpha / 2, 0, sqrt(Omega))
  R <- apply(rho(U, K), 2, sum)
  # plot(theta, R)
  theta[which.min(R)]
}


# -------------------------------------------------------------------
#' The bi-square Rho function
#'
#' Computes the objective function of the bi-square minimization problem in a location parameter, theta. The theta values are obtained internally by a grid search over the range of \code{\link[robustDIF]{y_fun}}. Useful for graphically diagnosing local solutions.
#'
#' @inheritParams y_fun
#' @param alpha the desired false positive rate for flagging items with DIF.
#' @param grid.width the width of grid points.
#'
#' @return A named list with theta values and the corresponding rho values.
#' @export

# -------------------------------------------------------------------

rho_fun <-
  function(mle,
           par = "intercept",
           log = F,
           alpha = .05,
           grid.width = .05) {
    y <- y_fun(mle, par, log = log)
    theta <-
      seq(from = max(min(y), -2.5),
          to = min(max(y), 2.5),
          by = grid.width)
    var_fun <- function(theta) {
      var_y(mle, theta, par, log = log)
    }
    n.items <- length(y)
    n.theta <- length(theta)

    Y <- matrix(y, nrow = n.items, ncol = n.theta)
    Var.Y <- Reduce(cbind, lapply(theta, function(x)
      var_fun(x)))
    Theta <- matrix(theta,
                    nrow = n.items,
                    ncol = n.theta,
                    byrow = T)
    Var.Theta <- matrix(
      1 / apply(1 / Var.Y, 2, sum),
      nrow = n.items,
      ncol = n.theta,
      byrow = T
    )
    U <- (Y - Theta) / Var.Y
    Omega <- (Var.Y - Var.Theta) / Var.Y ^ 2
    K <- qnorm(1 - alpha / 2, 0, sqrt(Omega))
    r <- apply(rho(U, K), 2, sum)
    names(r) <- NULL

    #r2 <- apply(rho((Y - Theta) / sqrt(Var.Y), 1.96), 2, sum)
    list(theta = theta, rho = r)
  }

# -------------------------------------------------------------------
#' Estimate IRT scale parameters using the RDIF procedure.
#'
#' @description
#' Implements M-estimation of an IRT scale parameter using the bi-square loss function. Also returns the bi-square weights for each item. Weights with a value of zero indicate that the corresponding item was flagged as having DIF during estimation.
#'
#' Estimation can be performed using iteratively re-weighted least squares (IRLS) or Newton-Raphson (NR). Currently, only IRLS is implemented. NR can diverge with bi-square.
#'
#' @inheritParams y_fun
#' @param alpha the desired false positive rate for flagging items with DIF.
#' @param starting.value one of \code{c("med", "lts", "min_rho", "all")} or a numerical value to be used as the starting value. The default is \code{"all"} which returns the median of the other three.
#' @param tol convergence criterion for comparing subsequent values of estimate
#' @param maxit maximum number of iterations
#' @param method one of \code{c("irls", "nr")}. Currently, only IRLS is implemented.
#'
#' @return A named list containing the estimate of the IRT scale parameter, the bi-square weights, the number of iterations performed, and the value of the convergence criterion (difference of  estimate between subsequent iterations).
#'
#' @examples
#' # Item intercepts, using the built-in example dataset "rdif.eg"
#' \dontrun{rdif(mle = rdif.eg)}
#'
#' # Item slopes
#' \dontrun{rdif(mle = rdif.eg, par = "slope")}
#'
#' @export
# -------------------------------------------------------------------

rdif <-
  function(mle,
           par = "intercept",
           log = F,
           alpha = .05,
           starting.value = "all",
           tol = 1e-7,
           maxit = 100,
           method = "irls") {
    nit <- 0
    conv <- 1

    # Set up scaling function
    y <- y_fun(mle, par, log)

    # Starting value
    starts <- get_starts(mle, par, log, alpha)
    theta <- median(starts)
    if (starting.value == "med") {
      theta <- starts[1]
    }
    if (starting.value == "lts") {
      theta <- starts[2]
    }
    if (starting.value == "min_rho") {
      theta <- starts[3]
    }
    if (is.numeric(starting.value)) {
      theta <- starting.value
    }

    # IRLS loop
    while (nit < maxit & conv > tol) {
      var.y <- var_y(mle, theta, par, log)
      w <- bsq_weight(theta, y, var.y, alpha)
      new.theta <- sum(w * y / var.y) / sum(w / var.y)
      nit <- nit + 1
      conv <- abs(theta - new.theta)
      theta <- new.theta
    }
    list(
      est = new.theta,
      weights = w,
      n.iter = nit,
      epsilon = conv
    )
  }

# -------------------------------------------------------------------
#' The R-DIF test of a single item parameter.
#'
#' Tests for DIF in either the item intercepts / thresholds or slopes using the asymptotic z-test in Theorem 1 of Halpin 2022.
#'
#' @inheritParams y_fun
#'
#' @return A \code{data.frame} containing the value of the z.test and p(|z| > |z.test|), for each item parameter.
#'
#' @examples
#' # Test intercepts, using the built-in example dataset "rdif.eg"
#' \dontrun{
#' rdif_z_test(mle = rdif.eg)
#' }
#'
#' # Test slopes
#' \dontrun{
#' rdif_z_test(.mle = rdif.eg, par = "slope")
#' }
#' @export
# -------------------------------------------------------------------

rdif_z_test <- function(mle, par = "intercept", log = F) {
  theta <- rdif(mle, par, log)$est
  y <- y_fun(mle, par, log)
  var.y <- var_y(mle, theta, par, log)
  var.theta <- 1 / (sum(1 / var.y))
  z.test <- (y - theta) / sqrt(var.y - var.theta)
  p.val <- (1 - pnorm(abs(z.test))) * 2
  if (par == "intercept") {
    names <- mle$par.names$original[grep(".d", mle$par.names$internal)]
  } else {
    names <- mle$par.names$original[grep(".a", mle$par.names$internal)]
  }
  out <- data.frame(z.test = z.test, p.val = p.val)
  row.names(out) <- names
  out
}

# -------------------------------------------------------------------
#' The R-DIF test of all parameters for each item.
#'
#' Simultaneously tests for DIF in all the item intercepts/thresholds and slopes using an asymptotic chi-square test.
#'
#' @inheritParams cov_yz

#' @return A data.frame containing the value of the chi2.test, its df, and p-value.
#'
#' @importFrom Matrix diag
#' @importFrom Matrix bdiag
#' @importFrom Matrix t
#'
#' @examples
#' # Using the built-in example dataset "rdif.eg"
#' \dontrun{
#' rdif_chisq_test(mle = rdif.eg)
#'}
#'
#' @export
# -------------------------------------------------------------------

rdif_chisq_test <- function(mle, log = F) {
  theta.y <- rdif(mle)$est
  theta.z <- rdif(mle, par = "slope", log)$est
  n.items <- nrow(mle$est$group.1)
  n.item.pars <- ncol(mle$est$group.1)

  # Set up vectors of Q form
  y <- matrix(
    y_fun(mle, par = "intercept"),
    nrow = n.items,
    ncol = n.item.pars - 1,
    byrow = T
  ) - theta.y

  z <- y_fun(mle, par = "slope", log = log) - theta.z

  yz.cbind <- cbind(y, z)
  yz.list <- lapply(1:nrow(yz.cbind), function(i)
    yz.cbind[i,])
  yz.bdiag <- Matrix::bdiag(yz.list)
  vcov.yz <- cov_yz(theta.y, theta.z, mle, log)

  # Compute test
  chi.square <-
    Matrix::diag(Matrix::t(yz.bdiag) %*% solve(vcov.yz) %*% yz.bdiag)
  p.val <- 1 - pchisq(chi.square, n.item.pars)
  out <-
    data.frame(chi.square = chi.square,
               df = n.item.pars,
               p.val = p.val)
  item.names <- unique(substr(mle$par.names$original, 1,
                              unlist(
                                gregexpr(".", mle$par.names$original, fixed = T)
                              ) - 1))
  row.names(out) <- item.names
  out
}
