
###############################################################
# New functions to compute effects and tests for impact
# Used for IMPS 2023
# Note: Depends on source packages below but also overwrites
# some of those functions, so need to call source and
# then run the functions at the top of this file
source("robustDIF_functions_old.R")
source("sim_functions.R")
###############################################################

# -------------------------------------------------------------------
#' The gradient of \code{\link[robustDIF]{y_intercept}}.
#'
#' The gradient is take with respect to the item parameters. The parameter vector is organized as \code{c(a0, d0, a1, d1)}, repeated over items. The parameters \code{a0} and \code{d0} are the item slope and intercept in the reference group. The parameters \code{a1} and \code{d1} are the item slope and intercept in the comparison group.
#'
#' @param theta the IRT scale parameter. If omitted, uses item-specific scaling functions instead.
#' @inheritParams y_intercept
#' @return A matrix in which the columns are the gradient vectors of \code{\link[robustDIF]{y_intercept}}, for each item.
#'
#' @seealso \code{\link[y_intercept]{y_intercept}}
#' @export
# -------------------------------------------------------------------

grad_intercept <- function(theta = NULL, irt.mle) {
  if (is.null(theta)) {theta <- y_intercept(irt.mle)}
  a1 <- irt.mle$par1$a
  delta.a0 <- 0 / a1
  delta.d0 <- - 1 / a1
  delta.a1 <- - theta / a1
  delta.d1 <- 1 / a1
  grad.vec <- c(rbind(delta.a0, delta.d0, delta.a1, delta.d1))
  grad_mat(grad.vec)
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
#' @export
# -------------------------------------------------------------------

grad_slope <- function(theta = NULL, irt.mle, log = F) {
  a <- irt.mle$par0$a
  if (is.null(theta)) {theta <- y_slope(irt.mle, log)}
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
#' Compute the asymptotic variance of the IRT scaling functions.
#'
#'Computations use the delta method (see Equation 13 of Halpin, 2022). If \code{theta} is not provided, the variance is computed as usual (line 1 of Equation 13). If \code{theta} is provided, it is used to compute a robust estimate of the variance, under the null hypothesis of no DIF (line 2 of Equation 13).
#'
#' @param theta the IRT scale parameter. If omitted, uses item-specific scaling functions instead.
#' @inheritParams y_fun
#'
#' @return A vector that contains the variance of the IRT scaling function, for each item.
#'
#' @seealso \code{\link[robustDIF]{y_fun}}
#'
#' @importFrom Matrix diag
#' @export
# -------------------------------------------------------------------

var_y <- function(theta = NULL, irt.mle, par = "intercept", log = F) {
  grad.y <- grad_intercept(theta, irt.mle)
  if (par == "slope") {grad.y <- grad_slope(theta, irt.mle, log)}
  joint.vcov <- joint_vcov(irt.mle)
  vcov.y <- t(grad.y) %*% joint.vcov %*% grad.y
  Matrix::diag(vcov.y)
}


# -------------------------------------------------------------------
#' The bi-square rho function.
#'
#' If \code{abs(u) > k} , \code{rho(u) = 1}. Else, \code{rho(u) = (1 - (1 - (u/k)^2)^3)}.
#'
#' @param u Can be a single value, vector, or matrix.
#' @param k The tuning parameter. Can be a scalar or the same dimension as \code{u}.
#' @return The bi-square rho function.
#' @seealso \code{\link[robustDIF]{rho_fun}}
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

bsq_var_weight <- function(theta, irt.mle, par = "intercept", log = F, alpha = .05){
  y <- y_fun(irt.mle, par, log)
  var.y <- var_y(theta, irt.mle, par, log)
  u <- (y - theta)/var.y
  omega <- (var.y - 1/sum(1/var.y))/var.y^2
  k <- qnorm(1 - alpha/2, 0, sqrt(omega))
  numer <- psi_prime(u, k) / var.y
  denom <- sum(bsq_weight(theta, y, var.y, alpha) / var.y)
  numer/denom
}

ml_est <- function(irt.mle, par = "intercept", log = F, alpha = .05) {
  y <- y_fun(irt.mle, par = par, log = log)
  var.y <- var_y(theta = NULL, irt.mle, par = par, log = log)
  sum(y / var.y) / sum(1 / var.y)
}

var_delta <- function(theta, irt.mle, par = "intercept", log = F, alpha = .05) {
  var.y.mle <- var_y(theta = NULL, irt.mle, par, log)
  w.mle <- 1/var.y.mle / sum(1/var.y.mle)
  w.bsq <- bsq_var_weight(theta, irt.mle, par, log, alpha)
  sum((w.bsq - w.mle)^2 * var.y.mle)
}

sim_study_delta <- function(n.reps = 100, n.persons = 500, n.items = 15, n.biased = 0, bias = 0, impact = c(0, 1)){

  # Item hyper-parms
  a.lower <- .9
  a.upper <- 2.5
  b.lim <- 1.5

  # Sim loop for parallelization via mclapply
  loop <- function(i){

    # DGP
    a0 <- runif(n.items, a.lower, a.upper)
    b0 <- sort(runif(n.items, -b.lim, b.lim))
    b1 <- apply_bias(b0, n.biased, bias)
    d0 <- b0*a0
    d1 <- b1*a0
    x0 <- rnorm(n.persons)
    x1 <- impact[2]*rnorm(n.persons) + impact[1]

    # sample standarized.
    # x0 <- rnorm(n.persons)
    # x0 <- (x0 - mean(x0)) / sd(x0)
    # x1 <- rnorm(n.persons)
    # x1 <- (x1 - mean(x1)) / sd(x1)
    # x1 <- impact[2]*x1 + impact[1]

    # Data gen
    dat0 <- simdata(a0, d0, n.persons, '2PL', Theta = matrix(x0))
    dat1 <- simdata(a0, d1, n.persons, '2PL', Theta = matrix(x1))

    # Fit IRT models
    fit0 <- mirt(dat0, 1, SE = T)
    fit1 <- mirt(dat1, 1, SE = T)

    # DIF Procedures
    irt.mle <- get_irt_pars(list(fit0, fit1))
    rdif.theta <- rdif(irt.mle)$est
    mle.theta <- ml_est(irt.mle)
    delta <- rdif.theta - mle.theta
    var.delta <- var_delta(rdif.theta, irt.mle)
    z.delta <- delta / sqrt(var.delta)
    p  <- (1 - pnorm(abs(z.delta))) * 2

    list(true.theta = (mean(x1) - mean(x0)) / sd(x1),
         rdif.theta = rdif.theta,
         mle.theta = mle.theta,
         delta = delta,
         z.delta = z.delta,
         p = p)
  }
  parallel::mclapply(1:n.reps, loop)
}

### Sim study

# Data gen
n.reps = 500
n.persons = 500
bias = .5
impact = c(.5,1)
n.items = 10
n.biased = 1

# Null Distribution
null.250 <- sim_study_delta(n.reps, 250, n.items = 15, n.biased = 0, bias, impact)
null.500 <- sim_study_delta(n.reps, 500, n.items = 15, n.biased = 0, bias, impact)
null.1000 <- sim_study_delta(n.reps, 1000, n.items = 15, n.biased = 0, bias, impact)
null.5000 <- sim_study_delta(n.reps, 5000, n.items = 15, n.biased = 0, bias, impact)
ds.null <- c(null.250, null.500[1:500], null.1000, null.5000)
ds <- Reduce(rbind, ds.null)
ds <- as.data.frame(ds)
ds <- as.data.frame(lapply(ds, as.numeric))
ds$n = rep(c(250, 500, 1000, 5000), each = n.reps)

head(ds)
library(xtable)
xtable(t(as.matrix(tapply(ds$p, ds$n, function(x) mean(x < .05)))))

library(ggplot2)

# Hist of p values
ggplot(ds, aes(pnorm(z.delta))) +
  geom_histogram(binwidth = .1,
                 color = "white",
                 fill= "#A6CEE3",
                 right = TRUE,
                 origin = 0) +
  facet_wrap(~n)

# Hist of z-test
ggplot(ds, aes(z.delta)) +
  geom_histogram(binwidth = .2,
                 color = "white",
                 fill = "#A6CEE3",
                 right = TRUE,
                 origin = 0) +
  facet_wrap(~n)

# QQ plots of z
ggplot(ds, aes(sample = z.delta)) +
  stat_qq(color = "#A6CEE3") +
  geom_abline(intercept = 0, slope = 1) + facet_wrap(~n)


# Power
power.0  <- sim_study_delta(n.reps, 500, n.items = 15, n.biased = 0, bias, impact)
power.1  <- sim_study_delta(n.reps, 500, n.items = 15, n.biased = 1, bias, impact)
power.2  <- sim_study_delta(n.reps, 500, n.items = 15, n.biased = 2, bias, impact)
power.3  <- sim_study_delta(n.reps, 500, n.items = 15, n.biased = 3, bias, impact)
power.4  <- sim_study_delta(n.reps, 500, n.items = 15, n.biased = 4, bias, impact)
power.5  <- sim_study_delta(n.reps, 500, n.items = 15, n.biased = 5, bias, impact)
power.6  <- sim_study_delta(n.reps, 500, n.items = 15, n.biased = 6, bias, impact)
power.7  <- sim_study_delta(n.reps, 500, n.items = 15, n.biased = 7, bias, impact)
power.8  <- sim_study_delta(n.reps, 500, n.items = 15, n.biased = 8, bias, impact)

ds.power <- c(power.0, power.1, power.2, power.3, power.4, power.5, power.6, power.7, power.8)

ds <- Reduce(rbind, ds.power)
ds <- as.data.frame(ds)
ds <- as.data.frame(lapply(ds, as.numeric))
ds$N.DIF = rep(paste0("N.DIF = ", 0:8), each = n.reps)
head(ds)

round((0:8 * .5) /15, 3)
temp <- t(as.matrix(tapply(ds$p, ds$N.DIF, function(x) mean(x < .05))))
colnames(temp) <- 0:8
xtable::xtable(rbind(temp, round((0:8 * .5) /15, 3))[c(2,1),])

ds$`p < .05` <- ds$p < .05

# Hist of z-test
library(ggplot2)
ggplot(ds, aes(z.delta, group = `p < .05`)) +
  geom_histogram(aes(fill = `p < .05`),
                 color = "grey60",
                 right = TRUE,
                 origin = 0) +
  scale_fill_manual(values = c("#A6CEE3", "#002f5b")) +
  #scale_fill_manual(values = c("blue", "red")) +
  theme_bw() +
  xlab("d / SE(d)") +
  scale_x_continuous(breaks= c(-5, 0, 5)) +
  facet_wrap(. ~ N.DIF, nrow = 3)



# boxplots of estimators
gg.data <- data.frame(
  Value = c(ds$rdif.theta, ds$mle.theta),
  Estimator = rep(c("Robust", "Naive"), each = nrow(ds)),
  N.DIF = as.factor(rep(0:8, each = n.reps)))

ggplot(gg.data, aes(N.DIF, Value, fill = Estimator)) +
  geom_boxplot(color = "grey60")  +
  scale_fill_manual(values = c("#A6CEE3", "#002f5b")) +
  #scale_fill_manual(values = c("red", "blue")) +
  theme_bw() +
  xlab ("Number of items with DIF (out of 15)") +
  geom_abline(slope = 0, intercept = .5, lty = 2)

library(xtable)

delta <-round((0:8 * .5) /15, 3)
p <- xtable(t(as.matrix(tapply(ds$p, ds$N.DIF, function(x) mean(x < .05)))), digits = 3)
temp <- tapply(gg.data$Value, list(gg.data$Estimator, gg.data$N.DIF), mean)
d <- temp[2, ] - temp[1, ]
head(ds)

round(tapply(ds$true.theta - ds$mle.theta, list(ds$N.DIF), mean), 3)

round(tapply(ds$rdif.theta - ds$mle.theta, list(ds$N.DIF), mean), 3)


round(tapply(ds$delta, list(ds$N.DIF), mean), 3)
xtable::xtable(rbind(delta, d, p), digits = 3)


round(.5 * 1:8 / 15, 3)
