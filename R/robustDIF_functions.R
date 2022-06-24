# To do
# Revert the starts and grid function back to their ugly versions.
# Get on with programming the quadratic form test!


#-------------------------------------------------------------------
#' Extract and format 2PL item parms from mirt
#'
#' @param mirt.fit.2pl a mirt object (SingleGroupClass) for the 2PL model
#' @return a n.items X 2 data.frame of 2PL item parameter MLEs
#' @export
# -------------------------------------------------------------------
get_mirt_mle <- function(mirt.fit.2pl ){
 n.items <- mirt.fit.2pl@Data$nitems
 parms <- Reduce(rbind, coef(mirt.fit.2pl, printSE = T)[1:(n.items)])[, 1:2]
 parms <- parms[row.names(parms) == "par", ]
 parms <- data.frame(parms)
 names(parms) <- c("a", "d")
 parms
}

# -------------------------------------------------------------------
#' Extract and format 2PL VCOV from mirt
#'
#' @param mirt.fit.2pl a mirt object (SingleGroupClass) for the 2PL model
#' @return a n.items X n.items covariance matrix for item parms
#' @export
# -------------------------------------------------------------------

get_mirt_vcov <- function(mirt.fit.2pl) {
  v <- vcov(mirt.fit.2pl)
  row.names(v) <- paste0(c("a", "d"), rep(1:(nrow(v)/2), each = 2))
  colnames(v) <- row.names(v)
  v
}

# -------------------------------------------------------------------
#' Extract and format 2PL MLEs and their asymptotic VCOV matrix
#'
#' Takes a list of 2PL model fits and formats the MLEs and their VCOV matrix for use with \code{robustDIF}. Current implementation only supports lists of length 2 (i.e., two groups) and the first fit is treated as the reference group. The only type of fit currently supported is the \code{SingleGroupClass} of the \code{mirt} package. To use fits from other software, the user can manually format the MLEs and VCOVs in the same format as this function.
#'
#' @param fit.list a list of model fits.
#' @param type a string indicating the package that produced the fits
#' @return a named list of MLEs and VCOVs
#' @export
# -------------------------------------------------------------------

get_irt_mle <- function(fit.list, type = "mirt") {
 if (type == "mirt") {
  est <- lapply(fit.list, get_mirt_mle)
  v <- lapply(fit.list, get_mirt_vcov)
  out <- list(par0 = est[[1]], par1 = est[[2]],
              vcov0 = v[[1]], vcov1 = v[[2]])
  return(out)
 }
}

# -------------------------------------------------------------------
#' The scaling function for item intercepts.
#'
#' Computes Y = (d1 - d0)/a1, used to test DIF on it em intercepts.
#'
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#'
#' @return A n.items vector of Y values.
#' @export
# -------------------------------------------------------------------

y_intercept <- function(irt.mle) {
   (irt.mle$par1$d - irt.mle$par0$d) / irt.mle$par1$a
}


# -------------------------------------------------------------------
#' The scaling function for item slopes.
#'
#' Computes Y = a1/d0, or its log; used to test DIF on item slopes.
#'
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#' @param log Logical: return log Y instead of Y?
#'
#' @return A n.items vector of Y values.
#' @export
# -------------------------------------------------------------------

y_slope <- function(irt.mle, log = F) {
  y <- irt.mle$par1$a / irt.mle$par0$a
  if (log) { y <- log(y) }
  y
}


# -------------------------------------------------------------------
#' Computes IRT scaling functions.
#'
#' Computes the scaling function for item intercepts, Y = (d1 - d0)/a1, or for item slopes, Y = a1/d0.
#'
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#' @param par Assess DIF on item intercepts or item slopes? One of \code{c("intercept", "slope")}
#' @param log Logical: Return log Y instead of Y? Only applies if \code{par = "slope"}
#'
#' @return A n.items vector of Y values.
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
#' Compute the gradient of \code{y_intercept}
#'
#' For each item, the non-null elements of the gradient are organized as \code{c(a0, d0, a1, d1)}. The parameters a0 and d0 are the item slope and intercept in the reference group. The parameters a1 and d1 are the item slope and intercept in the comparison group.
#'
#' @param theta IRT scale parameter (mu / sigma)
#' @param irt.mle the output of \code{robustDIF::get_irt_mle}
#' @return A matrix in which the columns are the 4*n.item-dimensional gradient vectors of \code{y_intercept} for each item
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
#' @param theta The IRT scale parameter (sigma).
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#' @param log Return grad for log Y instead (Logical)?
#' @return A matrix in which the columns are the 4*n.item-dimensional gradient vectors of \code{y_slope} for each item.
#' @export
# -------------------------------------------------------------------

grad_slope <- function(theta, irt.mle, log = F) {
  a0 <- irt.mle$par0$a
  delta.a0 <- - theta / a0
  delta.d0 <- 0 / a0
  delta.a1 <- 1 / a0
  delta.d1 <- 0 / a0
  grad.vec <- c(rbind(delta.a0, delta.d0, delta.a1, delta.d1))

  if (log) {
    grad.vec <- c(rbind(delta.a0, delta.d0, delta.a1, delta.d1)) * a0 / irt.mle$par1$a
  }
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
  Delta <- matrix(grad.vec, nrow = 4*n.items, ncol = n.items)
  flag <- c(rep(1, times = 4), rep(0, times = 4*n.items))
  Flag <- matrix(rep(flag, times = n.items)[1:(4*n.items^2)],
                   nrow = 4*n.items, ncol = n.items)
  Delta*Flag
}

# -------------------------------------------------------------------
#' Helper function to merge two VCOV matrices.
#'
#' Puts two 2PL models into a single block diagonal VCOV matrix that is conformable with the output of the \code{robustDIF} gradient functions.
#'
#' @param irt.mle the output of \code{robustDIF::get_irt_mle}
#' @return a block diagonal VCOV matrix in which each block is the VCOV of an item
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
#' @param theta The IRT scale parameter (mu / sigma for intercepts, sigma for slopes).
#' @param irt.mle The output of \code{robustDIF::get_irt_mle}.
#' @param par Assess DIF on item intercepts or item slopes? One of \code{c("intercept", "slope")}
#' @param log Logical: return grad for log Y instead? Only applies if \code{par = "slope"}
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
#' Computes the weights of the bi-square function for use with IRLS
#'
#' @param theta IRT scale parameter
#' @param y the output of \code{robustDIF::y_fun}
#' @param var.y the variance of y
#' @param alpha the desired false positive rate for flagging items with DIF
#' @return bi-square weights for (y - theta) / var.y, for each item
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
#' Computes the staring values for IRLS (or NR)
#'
#' @param irt.mle the output of \code{robustDIF::get_irt_mle}
#' @param par Assess DIF on item intercepts or item slopes? One of \code{c("intercept", "slope")}
#' @param log Logical: return grad for log Y instead? Only applies if \code{par = "slope"}
#' @param drops logical: omit items if y changes substantially when changing ref group?
#' @return a vector with elements med(Y), LTS(Y), and min Rho
#' @export
# -------------------------------------------------------------------

get_starts <- function(irt.mle, par = "intercepts", log = F, drops = T){


# ALL this needs to be conditional on par == intercepts

    y1 <- y_fun(irt.mle)
    var.y1 <- var_y(0, irt.mle)

    # again, switching the reference group
    irt.mle2 <- irt.mle
    names(irt.mle2)[1:2] <- names(irt.mle)[2:1]
    y2 <- -1 * y_fun(irt.mle2)
    var.y2 <- var_y(0, irt.mle)

    # rescale by variance of latent trait
    s <- median(y1/y2)
    y2 <- s*y2
    var.y2 <- s^2*var.y2

    if (drops) {
      drops  <- abs((y1 - y2) / sqrt(min(var.y1, var.y2))) < 1.5
    } else {
      drops <- T
    }
    y <- c(y1[drops], y2[drops])

    # median residual
    s1 <- median(y)

    # lts with 50% of data
    s2 <- lts(y)

    # grid search for min of Rho function
    s3 <- grid(irt.mle, intercepts = intercepts, log = log, drops = drops)

    c(s1, s2, s3)
}

# -------------------------------------------------------------------
#' Compute the least trimmed squares (LTS) estimate of location
#'
#' @param y a vector of data points
#' @param p is the proportion of data points to trim
#' @return The LTS estimate of location
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
#' The bi-square rho function
#'
#' @param u can be a single value, vector, or matrix
#' @param k the tuning parameter of the bi-square
#' @return The bi-square rho function
#' @export
# -------------------------------------------------------------------

rho <- function(u, k = 1.96) {
  w <- (u / k)^2
  out <- k^2 / 6 * (1 - (1 - w)^3)
  out[abs(u) > k] <- k^2 / 6
  return(out)
}

# -------------------------------------------------------------------
#' Performs grid search for minimum value of the bi-square Rho function
#'
#' @param irt.mle the output of \code{robustDIF::get_irt_mle}
#' @param par Assess DIF on item intercepts or item slopes? One of \code{c("intercept", "slope")}
#' @param log Logical: return grad for log Y instead? Only applies if \code{par = "slope"}
#' @param alpha the desired false positive rate for flagging items with DIF
#' @param grid.width the width of grid points
#' @param drops passed from \code{get_starts}; leave at default value (\code{TRUE}) to keep all items.
#' @return The location parameter that minimizes the bi-square Rho function
#' @export
# -------------------------------------------------------------------

grid <- function(irt.mle, par = "intercepts", log = F, alpha = .05, grid.width = .05, drops = T){

# ALL this needs to be conditional on par == intercepts

  y1 <- y_fun(irt.mle)
  theta <- seq(from = max(min(y1), -1.5),
               to = min(max(y1), 1.5),
               by = grid.width)

  n.items <- length(y1)
  n.theta <- length(theta)

  Theta <- matrix(theta, nrow = n.items, ncol = n.theta, byrow = T)
  Y1 <- matrix(y1, nrow = n.items, ncol = n.theta)
  SE.Y1 <- Reduce(cbind, lapply(theta, function(x) sqrt(var_y(x, irt.mle))))
  U <- Y1 - Theta / SE.Y1

  # Again, switching reference group
  irt.mle2 <- irt.mle
  names(irt.mle2)[1:2] <- names(irt.mle)[2:1]
  y2 <- -1 * y_fun(irt.mle2)
  s <- median(y1 / y2)
  Y2 <- matrix(s * y2, nrow = n.items, ncol = n.theta)
  SE.Y2 <- Reduce(cbind, lapply(theta, function(x) s * sqrt(var_y(x, irt.mle2))))
  U2 <- Y2 - Theta / SE.Y2

  R <- apply(rho(U)*drops + rho(U2)*drops, 2, sum)
  # plot(theta, R)
  theta[which.min(R)]
}

# -------------------------------------------------------------------
#' IRLS for an IRT scale parameter using the bi-square function
#'
#' @param irt.mle the output of \code{robustDIF::irt_mle}
#' @param par Assess DIF on item intercepts or item slopes? One of \code{c("intercept", "slope")}
#' @param log Logical: return grad for log Y instead? Only applies if \code{par = "slope"}
#' @param alpha the desired false positive rate for flagging items with DIF
#' @param starting.value one of \code{c("med", "lts", "min_rho", "all")} or a numerical value to be used as the starting value. The default is \code{"all"} which returns the median of the other three.
#' @param tol convergence criterion for comparing subsequent values of estimate
#' @param maxit maximum number of iterations to perform
#'
#' @return A named list containing the estimate, weights, and number of iterations
#' @export
# -------------------------------------------------------------------

irls <- function(irt.mle, par = "intercepts", log = F, alpha = .05, starting.value = "all", tol = 1e-5, maxit = 100){
  nit <- 0
  conv <- 1

  # Set up scaling function
  y <- y_fun(irt.mle, intercepts, log)

  # Starting value
  starts <- get_starts(irt.mle, par, log)
  theta <- median(starts)
  if (starting.value == "med") {theta <- starts[1]}
  if (starting.value == "lts") {theta <- starts[2]}
  if (starting.value == "min_rho") {theta <- starts[3]}
  if (is.numeric(starting.value)) {theta <- starting.value}


  # Loop
   while(nit < maxit & conv > tol) {
    var.y <- var_y(theta, irt.mle, intercepts, log)
    w <- bsq_weight(theta, y, var.y, alpha)
    new.theta <- sum(w * y / var.y ) / sum(w / var.y)
    nit <- nit + 1
    conv <- abs(theta - new.theta)
    theta <- new.theta
  }
  list(theta = new.theta, weights = w, n.iter = nit, epsilon = conv)
}

# -------------------------------------------------------------------
#' R-DIF test for item intercepts
#'
#' @param theta IRT scale parameter
#' @param irt.mle the output of \code{robustDIF::irt_mle}
#' @param alpha the desired false positive rate for flagging items with DIF
#'
#' @return A n.items-dimensional vector of p-values
#' @export
# -------------------------------------------------------------------

d_test <- function(theta, irt.mle, alpha = .05) {
  y <- y_fun(irt.mle)
  var.y <- var_y(scale, irt.mle)
  var.theta <- 1/(sum(1/var.y))
  crit <- qnorm(1 - alpha/2)
  abs(y - theta) / sqrt(var.y - var.theta) > crit
}

