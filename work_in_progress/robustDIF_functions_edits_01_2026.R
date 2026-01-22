# To do
# Diagnostic plot method for rdif
# Multi-parameter extension of rdif_z_test (formerly rdif_chisq_test)

#' ------------------------------------------------------------------- 
#' Example data set for R-DIF functions.
#'
#' A named list containing the maximum likelihood estimates and their estimated covariance matrix, for the 2PL IRT model fitted to 5 items in two independent groups. The first item has additive bias of .5 applied to the intercept only. The groups have a mean difference of .5 standard deviations on the latent trait. The variances of the latent trait are equal in each group.
#'
#' @format A named list with 4 components:
#' \describe{
#'   \item{par0}{A \code{data.frame} or named \code{list} with \code{par0$a} containing the item slopes and \code{par0$d} containing the item intercepts, for the reference group.}
#'   \item{par1}{The item parameter estimates of the comparison groups. See \code{par0} for formatting.}
#'   \item{vcov0}{The covariance matrix of \code{par0}, formatted as either a \code{data.frame} or \code{matrix}. The parameters should be organized by item, with the slope parameter coming first and the intercept parameter coming second (e.g., \code{a.item1, d.item1, a.item2, d.item2, ...}).}
#'    \item{vcov1}{The covariance matrix of \code{par1}. See \code{vcov0} for formatting.}
#'}
"rdif.eg"

# -------------------------------------------------------------------
#' The R-DIF scaling functions for item intercepts / thresholds.
#'
#' Computes \code{(d2 - d1)/a} for each threshold (d) of each item in groups g = {1, 2}. The parameter 'a' depends on the value of \code{type}: 
#'   
#'
#' @param mle the output of \code{\link[robustDIF]{get_mle}}
#' @param type a number in \code{1:3} indicating which version of delta to compute. See description for details.
#' @return The vector of scaling function values.

#' @description 
#'  Computes the scaling function to be used from the item thresholds (d) in groups = {1, 2}. The options are:
#'  \itemize{
#'  \item \code{type = 1}: computes \code{d.fun1 = (d2 - d1)/a1}
#'  \item \code{type = 2}: computes \code{d.fun2 = (d2 - d1)/a2}
#'  \item \code{type = 3}: computes \code{d.fun3 = (d2 - d1)/sqrt{(a1^2 + a2^2)/2}}
#'  }    

#' @export
# -------------------------------------------------------------------

d_fun <- function(mle, type = 3) {
  numerator <- mle$est$group.2[, -1] - mle$est$group.1[, -1]
  
  if (type == 1) {
    y <- numerator /  mle$est$group.1$a1 
  } else if (type == 2) {
    y <- numerator /  mle$est$group.2$a1 
  } else if (type == 3) {
    y <- numerator /  sqrt((mle$est$group.1$a1^2 + mle$est$group.2$a1^2) / 2)
  } else {
    stop("type must be one a number in {1:3}")
  }  
  y <- c(t(as.matrix(y)))
  names(y) <- paste0(
    rep(rownames(mle$est$group.1), each = ncol(mle$est$group.1) - 1),
    "_d",
    rep(1:(ncol(mle$est$group.1) - 1), times = nrow(mle$est$group.1))
  )
  y
}    

# -------------------------------------------------------------------
#' The R-DIF scaling function for item slopes.
#'
#' Computes the scaling function \code{a2/a1} for item slopes (a) in groups g = {1, 2}
#'
#' @param mle the output of \code{\link[robustDIF]{get_mle}}
#' @param log logical: return of log(a2/a1)?
#'
#' @return The vector of scaling function values.
#' @export
# -------------------------------------------------------------------

a_fun <- function(mle, log = F) {
  y <- mle$est$group.2$a1  / mle$est$group.1$a1
  if (log) { y <- log(y) }
  names(y) <- paste0(rownames(mle$est$group.1),"_a")
  y
}

# -------------------------------------------------------------------
#' R-DIF scaling functions for item intercepts/thresholds and slopes.
#'
#' Computes the scaling function specified by \code{fun}, for each item. 
#'
#' @param mle the output of \code{\link[robustDIF]{get_mle}}
#' @param fun one of \code{c("a_fun1" "a_fun2", d_fun1", "d_fun2", "d_fun3")}. See description for details.
#'
#' @return A vector of scaling function values.
#' @description 
#' Computes the a scaling function using the item thresholds (d) and slopes (a) in groups = {0, 1}. The options are:
#' \itemize{
#' \item \code{"a_fun1"}: computes \code{a2/a1}
#' \item \code{"a_fun2"}: computes \code{log(a2/a1)}
#' \item \code{"d_fun1"}: computes \code{(d2 - d1)/a1} for each threshold
#' \item \code{"d_fun2"}: computes \code{(d2 - d1)/a2} for each threshold
#' \item \code{"d_fun3"}: computes \code{(d2 - d1)/sqrt{(a1^2 + a1^2)/2} for each threshold
#' }  
#' @export
# -------------------------------------------------------------------

y_fun <- function(mle, fun = "d_fun3") {
  if (fun == "a_fun1") {y <- a_fun(mle, log = F)}
  if (fun == "a_fun2") {y <- a_fun(mle, log = T)}
  if (fun == "d_fun1") {y <- d_fun(mle, type = 1)}
  if (fun == "d_fun2") {y <- d_fun(mle, type = 2)}
  if (fun == "d_fun3") {y <- d_fun(mle, type = 3)}
  if (fun %in% c("a_fun1", "a_fun2", "d_fun1", "d_fun2", "d_fun3")) {
    y
  } else {
    stop("fun must be one of c('a_fun1', 'a_fun2', 'd_fun1', 'd_fun2', 'd_fun3')")
  }
}  

# -------------------------------------------------------------------
#' The gradient matrix of \code{\link[robustDIF]{d_fun}}.
#'
#' The gradient is taken with respect to the item parameters and organized to be conformable with \code{Matrix::bdiag(mle$var.cov)}. When evaluating the gradient under the null hypothesis of no DIF, the optional argument \code{theta} can be provided. It replaces the item-specific values of d_fun in the gradient computation.
#' 
#' @inheritParams d_fun
#' @param theta (optional) the scaling parameter. Replaces item-specific values of d_fun if provided. 
#' @return A matrix in which the columns are the gradient vectors of \code{\link[robustDIF]{d_fun}}, for each item and threshold. 
#' @seealso \code{\link[robustDIF]{d_fun}}
#' @export
# -------------------------------------------------------------------

grad_d <- function(mle, theta = NULL, type = 3) {
  n.items <- nrow(mle$est$group.1)
  n.thresholds <- ncol(mle$est$group.1) - 1
  n.groups <- length(mle$est)
  n.pars <- n.items * (n.thresholds + 1) * n.groups
  
  # a1, a2, and theta are ordered to repeat thresholds within items 
  a1 <- rep(mle$est$group.1$a1, each = n.thresholds)
  a2 <- rep(mle$est$group.2$a1, each = n.thresholds)
  
  if (is.null(theta)) {
    theta <- d_fun(mle, type)
  } else if (length(theta) == 1) {
    theta <- rep(theta, times = n.items*n.thresholds)
  } else {
    stop("theta must be null or a numeric of length 1")
  }
  
  # Components of gradient, each col ordered as (a1, d1, a2, d2)
  grad.mat <- matrix(0, nrow = 4, ncol = n.items*n.thresholds)
  if (type == 1) {
    grad.mat[1, ] <-  -1/a1 * theta
    grad.mat[2, ] <-  -1/a1   
    grad.mat[3, ] <-  0   
    grad.mat[4, ] <-  1/a1
  } else if (type == 2) { 
    grad.mat[1, ] <-  0
    grad.mat[2, ] <-  -1/a2
    grad.mat[3, ] <-  -1/a2 * theta   
    grad.mat[4, ] <-  1/a2
  } else if (type == 3) {
    grad.mat[1, ] <-  -a1 / (a1^2 + a2^2) * theta   
    grad.mat[2, ] <-  -((a1^2 + a2^2)/2)^-0.5   
    grad.mat[3, ] <-  -a2 / (a1^2 + a2^2) *theta   
    grad.mat[4, ] <-  ((a1^2 + a2^2)/2)^-0.5      
  } else {
    stop("type must be one a number in {1:3}")
  }  
  
  # Convert each column of grad.mat to an n.pars by 1 vector conformable with VCOV
  grad.list <- list()
  template <- rep(0, times = n.pars)
  k <- 1
  for(i in 1:n.items){
    for(j in 1:n.thresholds){
      ind1 <- (i-1) * (n.thresholds + 1) + 1
      ind2 <- ind1 + j
      ind3 <- ind1 + n.pars/2
      ind4 <- ind2 + n.pars/2
      template [c(ind1, ind2, ind3, ind4)] <- grad.mat[, k] 
      grad.list[[k]] <- template 
      template <- template * 0
      k <- k + 1
    }
  }
  Reduce(cbind, grad.list)
}

# -------------------------------------------------------------------
#' The gradient matrix of \code{\link[robustDIF]{a_fun}}.
#'
#' The gradient is taken with respect to the item parameters and organized to be conformable with \code{Matrix::bdiag(mle$var.cov)}. When evaluating the gradient under the null hypothesis of no DIF, the optional argument \code{theta} can be provided. It replaces the item-specific values of a_fun in the gradient computation.
#' 
#' @inheritParams a_fun
#' @param theta (optional) the scaling parameter. Replaces item-specific values of alpha if provided.
#' @return A matrix in which the columns are the gradient vectors of \code{\link[robustDIF]{alpha}}, for each item. 
#' @seealso \code{\link[robustDIF]{alpha}}
#' @export
# -------------------------------------------------------------------

grad_a <- function(mle, theta = NULL, log = F) {
  n.items <- nrow(mle$est$group.1)
  n.thresholds <- ncol(mle$est$group.1) - 1
  n.groups <- length(mle$est)
  n.pars <- n.items * (n.thresholds + 1) * n.groups
  
  # a1, a2, and theta are ordered for one scaling function per items 
  a1 <- mle$est$group.1$a1
  a2 <- mle$est$group.2$a1
  
  if (is.null(theta)) {
    theta <- c(a_fun(mle, log))
  } else if(length(theta) == 1) {
    theta <- rep(theta, times = n.items)
  } else {
    stop("theta must be null or a numeric of length 1")
  }
  
  # Components of gradient, each col ordered as (a1, a2)
  grad.mat <- matrix(0, nrow = 2, ncol = n.items)
  if (log) {
    grad.mat[1, ] <-  -1/a1 
    grad.mat[2, ] <-  1/a2  
  } else {
    grad.mat[1, ] <-  -1/a1 * theta
    grad.mat[2, ] <-  1/a1
  }
  
  # Convert each column of grad.mat to an n.pars by 1 vector conformable with VCOV
  grad.list <- list()
  template <- rep(0, times = n.pars)
  for(i in 1:n.items){
    ind1 <- (i-1) * (n.thresholds + 1) + 1
    ind2 <- ind1 + n.pars/2
    template [c(ind1, ind2)] <- grad.mat[, i] 
    grad.list[[i]] <- template 
    template <- template * 0
  }
  Reduce(cbind, grad.list)
}

# -------------------------------------------------------------------
#' The covariance matrix of IRT scaling functions.
#'
#' When evaluating the covariance matrix under the null hypothesis of no DIF, the optional argument \code{theta} can be provided. It replaces the item-specific scaling functions in the gradient computation. Type should be the same as used in \code{\link[robustDIF]{y_fun}}.
#'
#' @inheritParams y_fun
#' @param theta (optional) the scaling parameter. Replaces item-specific scaling functions if provided.
#
#' @return The covariance matrix of \code{y_fun}. 
#'
#' @seealso \code{\link[robustDIF]{y_fun}}
#'
#' @importFrom Matrix t bdiag  
# -------------------------------------------------------------------

vcov_y <- function(mle, theta = NULL, fun = "d_fun3") {
  if (fun == "a_fun1") {grad <- grad_a(mle, theta, log = F)}
  if (fun == "a_fun2") {grad <- grad_a(mle, theta, log = T)}
  if (fun == "d_fun1") {grad <- grad_d(mle, theta, type = 1)}
  if (fun == "d_fun2") {grad <- grad_d(mle, theta, type = 2)}
  if (fun == "d_fun3") {grad <- grad_d(mle, theta, type = 3)}
  if (fun %in% c("a_fun1", "a_fun2", "d_fun1", "d_fun2", "d_fun3")) {
    Matrix::t(grad)%*%Matrix::bdiag(mle$var.cov)%*%grad
  } else {
    stop("fun must be one of c('a_fun1', 'a_fun2', 'd_fun1', 'd_fun2', 'd_fun3')")
  }
}  

# -------------------------------------------------------------------
#' The bi-square rho function.
#'
#' If \code{abs(u) > k} , \code{rho(u) = 1}. Else, \code{psi(u) = 1 - (1 - (u/k)^2)^3}.
#'
#' @param u Can be a single value, vector, or matrix.
#' @param k The tuning parameter. Can be a scalar or the same dimension as \code{u}.
#' @return The bi-square rho function.
#' @export
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
#' The bi-square psi function.
#'
#' If \code{abs(u) > k} , \code{psi(u) = 0}. Else, \code{psi(u) = u(1 - (u/k)^2)^2}.
#'
#' @param u Can be a single value, vector, or matrix.
#' @param k The tuning parameter. Can be a scalar or the same dimension as \code{u}.
#' @return The bi-square psi function.
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
#' @return The bi-square psi-prime function.
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
#' The bi-square weight function.
#'
#' If \code{abs(u) > k} , \code{psi(u) = 0}. Else, \code{psi(u) = (1 - (u/k)^2)^2}.
#'
#' @param u Can be a single value, vector, or matrix.
#' @param k The tuning parameter. Can be a scalar or the same dimension as \code{u}.
#' @return The bi-square psi-prime function.
#' @export
# -------------------------------------------------------------------

bsq_weight <- function(u, k = 1.96) {
 if (length(k) != length(u)) {k <- k[1] + u - u}
  w <- (u / k)^2
  out <- (1 - w)^2
  out[abs(u) > k] <- 0
  out
}

# -------------------------------------------------------------------
#' Estimate IRT scale parameters using the robust DIF procedure.
#'
#' Implements M-estimation of an IRT scale parameter using the bi-square loss function. Also returns the bi-square weights for each item. #'
#' @inheritParams y_fun
#' @param alpha the desired false positive rate for flagging items with DIF.
#' @param starting.value one of \code{c("med", "lts", "min_rho", "all")} or a numerical value to be used as the starting value. See description for details.
#' @param tol convergence criterion for comparing subsequent values of estimate
#' @param maxit maximum number of iterations
#' @param method one of \code{c("irls", "newton")}. Currently, only IRLS is implemented.
#'
#' @return A named list containing the estimate of the IRT scale parameter, the bi-square weights, the number of iterations performed, and the value of the convergence criterion (difference of  estimate between subsequent iterations). If multiple solutions were found, the one with the lowest value of the bi-square objective function is returned and the other solutions are appended to the list as \code{other.solutions}.
#'
#' @description
#' Estimation can be performed using iteratively re-weighted least squares (IRLS) or Newton-Raphson (NR). Currently, only IRLS is implemented. If \code{starting.value = "all"}, three starting values are computed: the median of \code{\link[robustDIF]{y_fun}}, the least trimmed squares estimate of location for \code{\link[robustDIF]{y_fun}} with 50-percent trim rate, and the minimum of \code{\link[robustDIF]{rho_fun}}. The estimate is computed from each starting value, and the solution with the lowest value of the bi-square objective function is returned. If there are multiple solutions, the other solutions are appended to the output list as \code{other.solutions}.

#'
#' @examples
#' # Item intercepts, using the built-in example dataset "rdif.eg"
#' \dontrun{rdif(mle = rdif.eg, fun = "d_fun3)}
#'
#' # Item slopes
#' \dontrun{rdif(mle = rdif.eg, fun = "a_fun1")}
#'
#' @export

# -------------------------------------------------------------------

rdif <- function(mle, 
                 fun = "d_fun3", 
                 alpha = .05, 
                 starting.value = "all", 
                 tol = 1e-7, 
                 maxit = 100, 
                 method = "irls") {
  nit <- 0
  conv <- 1

  # Item-level scaling functions
  y <- y_fun(mle, fun)
  
  # Starting values
  starts <- get_starts(mle, fun, alpha)
  
  # # Tuning parameter
  k <- qnorm(1 - alpha/2, 0, 1)
  # if (abs(mean(y) - starts[3]) <= 1.5*abs(starts[1] - starts[3])) {
  #    k <- 4.685
  # } else {
  #    k <- qnorm(1 - alpha/2, 0, 1)
  # }
  # 
  if (starting.value == "med") {
    starts <- starts[1]
  }
  if (starting.value == "lts") {
    starts <- starts[2]
  }
  if (starting.value == "min_rho") {
    starts <- starts[3]
  }
  if (is.numeric(starting.value)) {
    starts <- starting.value
  }
  n.starts <- length(starts)
  sols <- vector("list", n.starts)
  
  # Loop over starting values
  for (i in 1:n.starts) { 
    theta <- starts[i]
    
  # IRLS 
    while(nit < maxit & conv > tol) {
      var.y <- Matrix::diag(vcov_y(mle, theta, fun))
      u <- (y - theta) / sqrt(var.y)
      w.star <- bsq_weight(u, k) 
      w <- (w_star / var.y) / sum(w.star / var.y)
      new.theta <- sum(w * y)
      
      # Convergence
      nit <- nit + 1
      conv <- abs(theta - new.theta)
      
      # Update
      theta <- new.theta
    } 
    
    # Compute final rho value
    u <- (y - theta) / sqrt(Matrix::diag(vcov_y(mle, theta, fun)))  
    rho.value <- sum(rho(u, k))
    
    # Store solution
    sols[[i]] <- list(
                  est = theta,
                  weights = w.star,
                  rho.value = rho.value,
                  n.iter = nit,
                  epsilon = conv,
                  k = k)
    
    # Reset for next starting value
    nit <- 0
    conv <- 1
  }   
  
  # Return solution(s)
  rho.values <- sapply(sols, function(x) x$rho.value)
  min.sol <- which.min(rho.values)[1]
  out <- sols[[min.sol]]
  out$multiple.solutions <- length(unique(round(rho.values, 3))) > 1
  if (out$multiple.solutions) {
    out$other.solutions <- sols[ -min.sol]
  }  
 out
}


# -------------------------------------------------------------------
#' Compute staring values for \code{\link[robustDIF]{rdif}}.
#'
#' @inheritParams y_fun
#' @param alpha the desired false positive rate for flagging items with DIF.

#' @return A vector containing the median of \code{\link[robustDIF]{y_fun}}, the least trimmed squares estimate of location for \code{\link[robustDIF]{y_fun}} with 50-percent trim rate, and the minimum of \code{\link[robustDIF]{rho_fun}}.
# -------------------------------------------------------------------

get_starts <- function(mle, fun = "d_fun3", alpha = .05){
  y <- y_fun(mle, fun)

  s1 <- median(y)
  s2 <- lts(y)

  rho.grid <- rho_grid(mle, fun, alpha)
  s3 <- rho.grid$theta[which.min(rho.grid$rho)]
  
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

lts <- function(y, p = 0.5) {
  n <- length(y)
  
  # If p is "proportion trimmed", convert to number kept:
  # h = number of observations used in the LTS fit.
  h <- floor(n * (1 - p))
  if (h < 1L) stop("p too large: no observations left after trimming.")
  if (h >= n) return(mean(y))  
  
  y_sorted <- sort(y)
  n.sample <- n - h + 1 
  
  # Compute variance for all contiguous windows of length h
  var.y <- numeric(n.sample)
  for (j in 1:n.sample) {
    subset_j <- y_sorted[j:(j + h - 1)]
    var.y[j] <- var(subset_j)
  }
  
  # Compute mean of window with minimal variance
  j.min <- which.min(var.y)
  best_subset <- y_sorted[j.min:(j.min + h - 1)]
  mean(best_subset) 
}

# -------------------------------------------------------------------
#' Compute a grid of bi-square Rho values
#'
#' Computes the objective function of the bi-square minimization problem in a location parameter, theta. The theta values are obtained internally by a grid search over the range of \code{\link[robustDIF]{y_fun}}. Used for starting values and graphically diagnosing local solutions.
#'
#' @inheritParams y_fun
#' @param alpha the desired false positive rate for flagging items with DIF.
#' @param grid.width the width of grid points.
#'
#' @return A named list with theta values and the corresponding Rho values.
#' @export

rho_grid <- function(mle, fun = "d_fun3", alpha = .05, grid.width = .01){

  y <- y_fun(mle, fun)
  theta <- seq(from = max(min(y), -2), to = min(max(y), 2), by = grid.width)
  n.items <- length(y)
  n.theta <- length(theta)
  k <- qnorm(1 - alpha/2, 0, 1)
  
  Y <- matrix(y, nrow = n.items, ncol = n.theta)
  Theta <- matrix(theta, nrow = n.items, ncol = n.theta, byrow = T)
  
  var_fun <- function(theta) {Matrix::diag(vcov_y(mle, theta, fun))}
  Var.Y <- Reduce(cbind, lapply(theta, function(x) var_fun(x)))
  
  U <- (Y - Theta) / sqrt(Var.Y)
  r <- apply(rho(U, k), 2, sum)
  names(r) <- NULL

  list(theta = theta, rho = r)
}

# -------------------------------------------------------------------
#' Wald tests of DIF for individual item parameters.
#'
#' Tests for DIF in each value of \code{y_fun} using asymptotic z-test. If \code{theta} is not provided,  \code{\link[robustDIF]{rdif}} is called internally using default values. 
#'
#' @inheritParams y_fun
#' @param theta (optional) the scaling parameter. 
#' @return A data.frame whose rows containing the results of the test for each item parameter. 
#'
#' @examples
#' # Test thresholds, using the built-in example dataset "rdif.eg"
#' \dontrun{
#' rdif.d <- rdif(mle = rdif.eg, fun = "d_fun3")
#' rdif_z_test(mle = rdif.eg, theta = rdif.d$est, fun = "d_fun3")
#' }
#'
#' # Test slopes
#' \dontrun{
#' rdif.a <- rdif(mle = rdif.eg, par = "a_fun1")
#' rdif_z_test(mle = rdif.eg, theta = rdif.d$est, fun = "a_fun1")
#' }
#' @export
# -------------------------------------------------------------------

dif_test <- function(mle, theta = NULL, fun = "d_fun3") {
  if (is.null(theta)) {
    rdif.out <- rdif(mle, fun)
    theta <- rdif.out$est
  }
  y <- y_fun(mle, fun)
  numerator <- y - theta
  
  vcov.y <- vcov_y(mle, theta, fun)
  var.y <- Matrix::diag(vcov.y)
  I <- diag(1, length(y))
  P <- matrix((1/var.y) / sum(1/var.y), nrow = length(y), ncol = length(y)) 
  demoninator <- sqrt(Matrix::diag(Matrix::t(I - P)%*%vcov.y%*%(I - P)))
  z.test <- numerator / demoninator
  p.val <- (1 - pnorm(abs(z.test))) * 2
  
  out <- data.frame(
    delta = numerator, 
    se = demoninator,
    z.test = z.test, 
    p.val = p.val)
  
  row.names(out) <- names(y)
  out
}

# -------------------------------------------------------------------
#' Wald tests of overall DIF in IRT scale parameters.
#'
#' Tests for overall DIF in the IRT scale parameter using several variations of the bi-square loss function.
#' @inheritParams y_fun
#' @param alpha the desired false positive rate for flagging items with DIF.
#' @return A data.frame whose rows containing the results of the test for each variation of the bi-square loss function. 
#' @examples
#' # Test thresholds, using the built-in example dataset "rdif.eg"
#' \dontrun{delta_test(mle = rdif.eg, fun = "d_fun3")}
#' @export
#' -------------------------------------------------------------------

delta_test <- function(mle, fun = "d_fun3", alpha = 0.05)
{
  # Set up
  y <- y_fun(mle, fun)
  n <- length(y)
  y.bar <- mean(y)
  rdif.out <- rdif(mle, fun, alpha)
  rdif.theta <- rdif.out$est
  delta <- y.bar - rdif.theta 
  
  vcov.y <- vcov_y(mle, theta = NULL, fun) # for sandwich
  var.y <- Matrix::diag(vcov_y(mle, theta = rdif.theta, fun)) # for bsq
  u <- (y - rdif.theta) / sqrt(var.y)
  k <-  rdif.out$k
  psi.prime <- psi_prime(u, k) 
  
  # Variance weights
  v.bar <- rep(1/n, n)
  v.psi.prime <- (pmax(psi.prime, 0)  / var.y) / sum(pmax(psi.prime, 0) / var.y)  
  
  # SEs
  se.y.bar <- sqrt(t(v.bar)%*%vcov.y%*%(v.bar))[1,1]
  se.rdif.theta <- sqrt(t(v.psi.prime)%*%vcov.y%*%(v.psi.prime))[1,1]
  se.delta <- sqrt(t(v.bar - v.psi.prime)%*%vcov.y%*%(v.bar - v.psi.prime))[1,1]
  
  # Delta test
  z <- delta / se.delta
  p.val <- (1 - pnorm(abs(z))) * 2
  
  #Output
  data.frame(
    naive.est = y.bar,
    naive.se = se.y.bar,
    rdif.est = rdif.theta,
    rdif.se = se.rdif.theta,
    delta = delta, 
    delta.se = se.delta,
    z.test = z,  
    p.val = p.val)
}

delta_test_from_dif <- function(dif.items, mle, fun = "d_fun3", alpha = 0.05)
{
  # Set up
  y <- y_fun(mle, fun)
  vcov.y <- vcov_y(mle, fun = fun)
  n <- length(y)
  n0 <- n - length(dif.items)
  
  # Weights
  w.bar <- rep(1/n, n)
  w.dif <- rep(1/n0, n)
  w.dif[dif.items] <- 0
  
  # Delta
  y.bar <- sum(w.bar*y)
  y.dif <- sum(w.dif*y)
  delta <- sum((w.bar - w.dif)*y)
  
  # SEs
  se.y.bar <- sqrt(t(w.bar)%*%vcov.y%*%(w.bar))[1,1]
  se.dif <- sqrt(t(w.dif)%*%vcov.y%*%(w.dif))[1,1]
  se.delta <- sqrt(t(w.bar - w.dif)%*%vcov.y%*%(w.bar - w.dif))[1,1]
  
  # Delta test
  z <- delta / se.delta
  p.val <- (1 - pnorm(abs(z))) * 2
  
  #Output
  data.frame(
    naive.est = y.bar,
    naive.se = se.y.bar,
    dif.est = y.dif,
    dif.se = se.dif,
    delta = delta, 
    delta.se = se.delta,
    z.test = z,  
    p.val = p.val)
}

    


