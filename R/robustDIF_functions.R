# To do:

#-------------------------------------------------------------------
#' Extract 2PL item parameter estimates from \link[mirt]{mirt}.
#'
#' @param mirt.fit.2pl A \link[mirt]{mirt} object (\code{SingleGroupClass}) estimated for the 2PL model.
#' @return A data frame of 2PL item parameter estimates, in slope-intercept form.
#' @export
# -------------------------------------------------------------------

get_mirt_pars <- function(mirt.fit.2pl ){
 n.items <- mirt.fit.2pl@Data$nitems
 parms <- Reduce(rbind, coef(mirt.fit.2pl, printSE = T)[1:(n.items)])[, 1:2]
 parms <- parms[row.names(parms) == "par", ]
 parms <- data.frame(parms)
 names(parms) <- c("a", "d")
 parms
}

# -------------------------------------------------------------------
#' Extract 2PL covariance matrix of item parameter estimates from \link[mirt]{mirt}.
#' @inheritParams get_mirt_pars
#' @return The covariance matrix of 2PL item parameter estimates.
#' @export
# -------------------------------------------------------------------

get_mirt_vcov <- function(mirt.fit.2pl) {
  v <- vcov(mirt.fit.2pl)
  row.names(v) <- paste0(c("a", "d"), rep(1:(nrow(v)/2), each = 2))
  colnames(v) <- row.names(v)
  v
}

# -------------------------------------------------------------------
#' Extract 2PL item parameter estimates their covariance matrix.
#'
#' Takes a list of 2PL model fits and formats the item parameter estimates  and their covariance matrix. The R-DIF procedure assumes that the estimates were obtained by maximum likelihood and the covariance is asymptotically correct.
#'
#' Note that current implementation only supports lists of length 2 (i.e., two groups) and the first fit is treated as the reference group. The only type of fit currently supported is the \code{SingleGroupClass} of the \link[mirt]{mirt} package. To use fits from other software, use the format \code{list(par0 = est[[1]], par1 = est[[2]], vcov0 = vcov[[1]], vcov1 = vcov[[2]])}
#'
#' @param fit.list A list of 2PL model fits.
#' @param type A string indicating the package that produced the fits.
#' @return A named list of item parameter estimates and covariance matrices.
#' @export
# -------------------------------------------------------------------

get_irt_pars <- function(fit.list, type = "mirt") {
 if (type == "mirt") {
  est <- lapply(fit.list, get_mirt_pars)
  v <- lapply(fit.list, get_mirt_vcov)
  out <- list(par0 = est[[1]], par1 = est[[2]],
              vcov0 = v[[1]], vcov1 = v[[2]])
  return(out)
 }
}

# -------------------------------------------------------------------
#' The R-DIF scaling function for item intercepts.
#'
#' Computes \code{Y = (d1 - d0)/a1}. Used to test DIF on item intercepts.
#'
#' @param irt.mle The output of \link[robustDIF::get_irt_pars]{get_irt_pars}.
#' @return A vector of Y values.
#' @export
# -------------------------------------------------------------------

y_intercept <- function(irt.mle) {
   (irt.mle$par1$d - irt.mle$par0$d) / irt.mle$par1$a
}

# -------------------------------------------------------------------
#' The R-DIF scaling function for item slopes.
#'
#' Computes \code{Y = a1/d0}, or its log. Used to test DIF on item slopes.
#'
#' @inheritParams y_intercept
#' @param log Logical: return log Y instead of Y?
#'
#' @return A vector of Y values.
#' @export
# -------------------------------------------------------------------

y_slope <- function(irt.mle, log = F) {
  y <- irt.mle$par1$a / irt.mle$par0$a
  if (log) { y <- log(y) }
  y
}


# -------------------------------------------------------------------
#' Computes R-DIF scaling functions.
#'
#' Computes the R-DIF scaling function for item intercepts or for item slopes.
#'
#' @inheritParams y_intercept
#' @param par Character: Use scaling function for item intercepts or item slopes? One of \code{c("intercept", "slope")}.
#' @param log Logical: return log Y instead of Y? Only applies if \code{par = "slope"}.
#'
#' @return A vector of Y values.
#' @seealso \link[robustDIF::y_intercept]{y_intercept}, \link[robustDIF::y_slope]{y_slope}.
#' @export
# -------------------------------------------------------------------

y_fun <- function(irt.mle, par = "intercept", log = F) {
  y <- y_intercept(irt.mle)
  if (par == "slope") {
    y <- y_slope(irt.mle, log)
  }
  y
}

# -------------------------------------------------------------------
#' Compute the gradient of \link[robustDIF::y_intercept]{y_intercept}.
#'
#' For each item, the gradient are organized as \code{c(a0, d0, a1, d1)}. The parameters \code{a0} and \code{d0} are the item slope and intercept in the reference group. The parameters \code{a1} and \code{d1} are the item slope and intercept in the comparison group.
#'
#' @param theta IRT scale parameter (mu / sigma).
#' @param irt.mle the output of \code{robustDIF::get_irt_mle}.
#' @return A matrix in which the columns are the 4*n.item-dimensional gradient vectors of \code{y_intercept} for each item.
#' @export
# -------------------------------------------------------------------

grad_intercept <- function(theta, irt.mle) {
  a1 <- irt.mle$par1$a
  delta.a0 <- 0 / a1
  delta.d0 <- - 1 / a1
  delta.a1 <- - theta / a1
  delta.d1 <- 1 / a1
  grad.vec <- c(rbind(delta.a0, delta.d0, delta.a1, delta.d1))
  grad_mat(grad.vec)
}

# -------------------------------------------------------------------
#' Compute the gradient of \code{y_slope}.
#'
#' For each item, the non-null elements of the gradient are organized as \code{c(a0, d0, a1, d1)}. The parameters a0 and d0 are the item slope and intercept in the reference group. The parameters a1 and d1 are the item slope and intercept in the comparison group.
#'
#' @param theta The IRT scale parameter (sigma or log sigma).
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#' @param log Logical: Use log of scaling function?.
#' @return A matrix in which the columns are the 4*n.item-dimensional gradient vectors of \code{y_slope} for each item.
#' @export
# -------------------------------------------------------------------

grad_slope <- function(theta, irt.mle, log = F) {
  a <- irt.mle$par0$a
  if (log) {
    theta <- exp(theta)
    a <- irt.mle$par1$a
  }
  delta.a0 <- - theta / a
  delta.d0 <- 0 / a
  delta.a1 <- 1 / a
  delta.d1 <- 0 / a
  grad.vec <- c(rbind(delta.a0, delta.d0, delta.a1, delta.d1))
  grad_mat(grad.vec)
}

# -------------------------------------------------------------------
#' Helper function to put gradient vectors into a matrix.
#'
#' Called internally by the \code{robustDIF} gradient functions.
#'
#' @param grad.vec A gradient vector.
#' @return A matrix in which the columns are the 4*n.item-dimensional gradient vectors for each item.
#' @export
# -------------------------------------------------------------------

grad_mat <- function(grad.vec) {
  n.items <- length(grad.vec) / 4
  Delta <- matrix(grad.vec, nrow = 4 * n.items, ncol = n.items)
  flag <- c(rep(1, times = 4), rep(0, times = 4 * n.items))
  Flag <- matrix(rep(flag, times = n.items)[1:(4 * n.items^2)],
                   nrow = 4 * n.items, ncol = n.items)
  Delta * Flag
}

# -------------------------------------------------------------------
#' Helper function to merge two VCOV matrices.
#'
#' Puts two 2PL models into a single block diagonal VCOV matrix that is conformable with the output of the \code{robustDIF} gradient functions.
#'
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#' @return A block diagonal VCOV matrix in which each block is the VCOV of an item.
#' @export
# -------------------------------------------------------------------

joint_vcov <- function(irt.mle) {
  n.items <- nrow(irt.mle$par0)
  m <- diag(1:n.items) %x% matrix(1, 2, 2)
  v0 <- lapply(split(irt.mle$vcov0, m)[-1], matrix, 2)
  v1 <- lapply(split(irt.mle$vcov1, m)[-1], matrix, 2)
  vfull <- vector("list", n.items * 2)
  vfull[1:(n.items) * 2 - 1] <- v0
  vfull[1:(n.items) * 2 ] <- v1
  Matrix::bdiag(vfull)
}

# -------------------------------------------------------------------
#' Computes the variance of asymptotic null distribution of IRT scaling functions.
#'
#' @param theta The IRT scale parameter (mu / sigma for intercepts, sigma or log sigma for slopes).
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#' @param par Use scaling function for item intercepts or item slopes? One of \code{c("intercept", "slope")}.
#' @param log Logical: Use log of scaling function? Only applies if \code{par = "slope"}.
#' @return An n.item-dimensional vector that contains the variance of the IRT scaling function, for each item.
#' @export
# -------------------------------------------------------------------

var_y <- function(theta, irt.mle, par = "intercept", log = F) {
  grad.y <- grad_intercept(theta, irt.mle)
  if (par == "slope") {
    grad.y <- grad_slope(theta, irt.mle, log)
  }
  joint.vcov <- joint_vcov(irt.mle)
  vcov.y <- t(grad.y) %*% joint.vcov %*% grad.y
  Matrix::diag(vcov.y)
}


# -------------------------------------------------------------------
#' Computes the asymptotic covariance of IRT scaling functions for the item intercepts and the item slopes, under the null hypothesis.
#'
#' Used to compute the the quadratic form test of item intercepts and slopes.
#'
#' @param theta.y The IRT scale parameter for the item intercepts (mu / sigma)
#' @param theta.z The IRT scale parameter for the item slopes (sigma or log sigma).
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#' @param log Logical: Use log of scaling function for the slope?
#' @return An n.item-dimensional vector that contains the covariance of the IRT scaling functions, for each item.
#' @export
# -------------------------------------------------------------------

cov_yz <- function(theta.y, theta.z, irt.mle, log = F) {
  grad.y <- grad_intercept(theta.y, irt.mle)
  grad.z <- grad_slope(theta.z, irt.mle, log = log)
  joint.vcov <- joint_vcov(irt.mle)
  vcov.yz <- t(grad.y) %*% joint.vcov %*% grad.z
  Matrix::diag(vcov.yz)
}

# -------------------------------------------------------------------
#' Computes the weights of the bi-square function for use with IRLS.
#'
#' @param theta The IRT scale parameter.
#' @param y The output of \code{robustDIF::y_fun}.
#' @param var.y The variance of y.
#' @param alpha The desired false positive rate for flagging items with DIF.
#' @return The bi-square weights for (y - theta) / var.y, for each item.
#' @export
# -------------------------------------------------------------------

bsq_weight <- function(theta, y, var.y, alpha = .05){
  r <- y - theta
  var.theta <- 1/sum(1/var.y)
  omega <- (var.y - var.theta)/var.y^2
  k <- qnorm(1 - alpha/2, 0, sqrt(omega))
  bsq.cuts <- var.y * k
  w <- (1 - (r / bsq.cuts)^2)^2
  w[abs(r) > bsq.cuts] <- 0
  w
}

# -------------------------------------------------------------------
#' Computes the staring values for IRLS (or NR).
#'
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#' @param par Use scaling function for item intercepts or item slopes? One of \code{c("intercept", "slope")}.
#' @param log Logical: Use log of scaling function? Only applies if \code{par = "slope"}.
#' @return A vector with elements med(Y), LTS(Y), and min Rho.
#' @export
# -------------------------------------------------------------------

get_starts <- function(irt.mle, par = "intercept", log = F, alpha = .05){

  y <- y_fun(irt.mle, par, log = log)
  var_fun <- function(theta) {var_y(theta, irt.mle, par, log = log)}

  #This is combines both choices of reference group for intercepts
  if (par == "intercept") {
    y1 <- y
    irt.mle2 <- irt.mle
    names(irt.mle2)[1:2] <- names(irt.mle)[2:1]
    y2 <- -1 * y_fun(irt.mle2)

    s <- median(y1/y2)
    y2 <- s * y2

    # Drop items if y1 - y2 / sd(y) > 1.5
    var.y1 <- var_y(median(y1), irt.mle)
    var.y2 <- s^2 * var_y(median(y2), irt.mle)
    drops <- abs((y1 - y2) / sqrt(min(var.y1, var.y2))) < 1.5
    y <- c(y1[drops], y2[drops])
    var_fun <- function(theta) {
      c(var_y(theta, irt.mle)[drops], s^2 * var_y(theta, irt.mle2)[drops])
    }
  }

  # median residual
  s1 <- median(y)

  # lts with 50% of data
  s2 <- lts(y)

  # grid search for min of Rho function
  s3 <- rho_grid(y, var_fun, alpha = .05)

  c(s1, s2, s3)
}

# -------------------------------------------------------------------
#' Compute the least trimmed squares (LTS) estimate of location
#'
#' @param y A vector of data points.
#' @param p The proportion of data points to trim.
#' @return The LTS estimate of location.
#' @export
# -------------------------------------------------------------------

lts <- function(y, p = .5){
  n <- length(y)
  n.per.sample <- floor(n*p)
  n.sample <- ceiling(n - n.per.sample + 1)
  Y <- matrix(sort(y), nrow = n, ncol = n.sample)
  flag <- c(rep(1, times = n.per.sample),
            rep(0, times = n.sample))
  Flag <- matrix(rep(flag, times = n.sample)[1:(n * n.sample)],
                   nrow = n, ncol = n.sample)
  dat <- Y * Flag
  sum(dat[,which.min(apply(dat, 2, var))]) / n.per.sample
}

# -------------------------------------------------------------------
#' The bi-square rho function.
#'
#' If abs(u) > k , rho(u) = 1. Else, rho(u) = (1 - (1 - (u/k))^3).
#'
#' @param u Can be a single value, vector, or matrix.
#' @param k The tuning variable. Can be a scalar or the same dimension as \code{u}.
#' @return The bi-square rho function.
#' @export
# -------------------------------------------------------------------

rho <- function(u, k = 1.96) {
  if (length(k) != length(u)) {k <- k[1] + u - u}
  w <- (u / k)^2
  out <- 1 - (1 - w)^3
  out[abs(u) > k] <- 1
  out
}

# -------------------------------------------------------------------
#' Performs grid search for minimum value of the bi-square Rho function.
#'
#' Internal function used by \code{get_starts} to get staring values based on minimizing the rho-function of the bi-square. For plotting use \code{rho_fun} instead.
#'
#' @param y Output of \code{y_fun}
#' @param var_fun Computes variance of \code{y} as a function of only the scale parameter.
#' @param alpha The desired false positive rate for flagging items with DIF.
#' @param grid.width The width of grid points.
#'
#' @return The location parameter that minimizes the bi-square Rho function.

# -------------------------------------------------------------------

rho_grid <- function(y, var_fun, alpha = .05, grid.width = .05){

  theta <- seq(from = max(min(y), -1.5),
               to = min(max(y), 1.5),
               by = min(grid.width, (max(y) - min(y)) / 10))

  n.items <- length(y)
  n.theta <- length(theta)

  Y <- matrix(y, nrow = n.items, ncol = n.theta)
  Var.Y <- Reduce(cbind, lapply(theta, function(x) var_fun(x)))
  Theta <- matrix(theta, nrow = n.items, ncol = n.theta, byrow = T)
  Var.Theta <- matrix(1/apply(1/Var.Y, 2, sum),
                      nrow = n.items, ncol = n.theta, byrow = T)
  Omega <- (Var.Y - Var.Theta) / Var.Y^2
  U <- (Y - Theta) / Var.Y
  K <-  qnorm(1 - alpha/2, 0, sqrt(Omega))
  R <- apply(rho(U, K), 2, sum)
  # plot(theta, R)
  theta[which.min(R)]
}


# -------------------------------------------------------------------
#' Computes the bi-square Rho function
#'
#'  Plotting the objective function of the minimization problem in theta. Useful for diagnosing local solutions.
#'
#' @param theta The IRT scale parameter (mu / sigma for intercepts, sigma or log sigma for slopes).
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#' @param par Use scaling function for item intercepts or item slopes? One of \code{c("intercept", "slope")}.
#' @param log Logical: Use log of scaling function? Only applies if \code{par = "slope"}.
#' @param alpha The desired false positive rate for flagging items with DIF.
#' @param grid.width The width of grid points for \code{theta}.
#'
#' @return A named list with theta and Rho values, can be passed directly to plot.

# -------------------------------------------------------------------

rho_fun <- function(irt.mle, par = "intercept", log = F, alpha = .05, grid.width = .05){

  y <- y_fun(irt.mle, par, log = log)
  theta <- seq(from = max(min(y), -1.5), to = min(max(y), 2.5), by = grid.width)

  var_fun <- function(theta) {var_y(theta, irt.mle, par, log = log)}
  n.items <- length(y)
  n.theta <- length(theta)

  Y <- matrix(y, nrow = n.items, ncol = n.theta)
  Var.Y <- Reduce(cbind, lapply(theta, function(x) var_fun(x)))
  Theta <- matrix(theta, nrow = n.items, ncol = n.theta, byrow = T)
  Var.Theta <- matrix(1/apply(1/Var.Y, 2, sum),
                      nrow = n.items, ncol = n.theta, byrow = T)
  U <- (Y - Theta) / Var.Y
  Omega <- (Var.Y - Var.Theta) / Var.Y^2
  K <- qnorm(1 - alpha/2, 0, sqrt(Omega))
  r <- apply(rho(U, K), 2, sum)

  #r2 <- apply(rho((Y - Theta) / sqrt(Var.Y), 1.96), 2, sum)
  list(theta = theta, rho = r)
}

# -------------------------------------------------------------------
#' Estimates IRT scale parameters using the RDIF procedure
#'
#' Estimation is via iteratively reweighted least squares using the bi-square loss function.
#'
#' @param irt.mle the output of \code{robustDIF::irt_mle}
#' @param par Use scaling function for item intercepts or item slopes? One of \code{c("intercept", "slope")}.
#' @param log Logical: Use log of scaling function? Only applies if \code{par = "slope"}.
#' @param alpha the desired false positive rate for flagging items with DIF
#' @param starting.value one of \code{c("med", "lts", "min_rho", "all")} or a numerical value to be used as the starting value. The default is \code{"all"} which returns the median of the other three.
#' @param tol convergence criterion for comparing subsequent values of estimate
#' @param maxit maximum number of iterations to perform
#'
#' @return A named list containing the estimate, weights, and number of iterations
#' @export
# -------------------------------------------------------------------

rdif <- function(irt.mle, par = "intercept", log = F, alpha = .05, starting.value = "all", tol = 1e-7, maxit = 100){
  nit <- 0
  conv <- 1

  # Set up scaling function
  y <- y_fun(irt.mle, par, log)
  # Starting value
  starts <- get_starts(irt.mle, par, log, alpha)
  theta <- median(starts)
  if (starting.value == "med") {theta <- starts[1]}
  if (starting.value == "lts") {theta <- starts[2]}
  if (starting.value == "min_rho") {theta <- starts[3]}
  if (is.numeric(starting.value)) {theta <- starting.value}

  # Loop
   while(nit < maxit & conv > tol) {
    var.y <- var_y(theta, irt.mle, par, log)
    w <- bsq_weight(theta, y, var.y, alpha)
    new.theta <- sum(w * y / var.y ) / sum(w / var.y)
    nit <- nit + 1
    conv <- abs(theta - new.theta)
    theta <- new.theta
  }
  list(est = new.theta, weights = w, n.iter = nit, epsilon = conv)
}

# -------------------------------------------------------------------
#' R-DIF z-test.
#'
#' Tests for DIF in either the item intercept or slopes using a asymptotic z-test.
#'
#' @param theta The IRT scale parameter.
#' @param irt.mle The output of \code{robustDIF::irt_mle}.
#' @param par Use scaling function for item intercepts or item slopes? One of \code{c("intercept", "slope")}.
#' @param log Logical: Use log of scaling function? Only applies if \code{par = "slope"}.
#'
#' @return A named list with containing the value of the z.test and p(|z| > z.test), for each item.
#' @export
# -------------------------------------------------------------------

z_test <- function(theta, irt.mle, par = "intercept", log = F) {
  y <- y_fun(irt.mle, par, log)
  var.y <- var_y(theta, irt.mle, par, log)
  var.theta <- 1/(sum(1/var.y))
  z.test <- (y - theta) / sqrt(var.y - var.theta)
  p.val <- (1 - pnorm(abs(z.test))) * 2
  list(z.test = z.test, p.val = p.val)
}

# -------------------------------------------------------------------
#' R-DIF chi-square test.
#'
#' Simultaneously tests  for DIF in both the item intercept or slopes using a asymptotic chi-square test.
#'
#'
#' @param theta.y The IRT scale parameter for the item intercepts (mu / sigma)
#' @param theta.z The IRT scale parameter for the item slopes (sigma or log sigma).
#' @param irt.mle The output of \code{robustDIF::irt_mle}
#' @param log Logical: Use log of scaling function for item slopes?
#'
#' @return A named list with containing the value of the chi2.test and p(chi.square > chi2.test), for each item.
#' @export
# -------------------------------------------------------------------

chi2_test<- function(theta.y, theta.z, irt.mle, log = F) {

  # Set up vectors of Q form
  u <- y_fun(irt.mle, par = "intercept") - theta.y
  v <- y_fun(irt.mle, par = "slope", log = log) - theta.z
  uv.cbind <- cbind(u, v)
  uv.list <- lapply(1:nrow(uv.cbind), function(i) uv.cbind[i,])
  uv.bdiag <- Matrix::bdiag(uv.list)

  # Set up Sigma matrix of Q form
  var.y <- var_y(theta.y, irt.mle, par = "intercept")
  var.theta.y <- 1/sum(1/var.y)
  var.z <- var_y(theta.z, irt.mle, par = "slope", log = log)
  var.theta.z <- 1/sum(1/var.z)
  cov.yz <- cov_yz(theta.y, theta.z, irt.mle, log = log)

  Sigma11 <- var.y - var.theta.y
  Sigma22 <- var.z - var.theta.z
  Sigma12 <- cov.yz * (1 - var.theta.y/var.y - var.theta.z/var.z) +
             sum(var.theta.y/var.y * var.theta.z/var.z * cov.yz)

  Sigma.cbind <- cbind(Sigma11, Sigma12, Sigma12, Sigma22)
  Sigma.list <- lapply(1:nrow(Sigma.cbind),
                       function(i) matrix(Sigma.cbind[i,], nrow = 2, ncol = 2))
  Sigma.bdiag <- Matrix::bdiag(Sigma.list)

  # Compute test
  chi.square <- Matrix::diag(Matrix::t(uv.bdiag) %*% solve(Sigma.bdiag) %*% uv.bdiag)
  p.val <- 1 - pchisq(chi.square, 2)
  list(chi.square = chi.square, p.val = p.val)
}

