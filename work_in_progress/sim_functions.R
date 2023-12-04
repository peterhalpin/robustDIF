# data sim functions

apply_bias <- function(b0, n.biased = 1, max.bias = .5, worst.case = T, multiplicative = F) {
  n.items <- length(b0)
  biased.items <- sample(n.items, n.biased)
  bias.vector <- 0 * b0
  if (worst.case) {
    bias.vector[biased.items] <- max.bias
  } else {
    bias.vector[biased.items] <- runif(n.biased, max.bias/2, max.bias)
  }
  if (!multiplicative) {out <- b0 + bias.vector}
  if (multiplicative) {out <- b0 + bias.vector *  b0}
  return(out)
}

sim_dat <- function(n.persons = 500, n.items = 15, n.biased = 0, bias = 0, impact = c(0, 1)){

  # Item hyper-parms
  a.lower <- .9
  a.upper <- 2.5
  b.lim <- 1.5
  a0 <- runif(n.items, a.lower, a.upper)
  b0 <- sort(runif(n.items, -b.lim, b.lim))
  b1 <- apply_bias(b0, n.biased, bias)
  d0 <- b0*a0
  d1 <- b1*a0
  x0 <- rnorm(n.persons)
  x1 <- impact[2]*rnorm(n.persons) + impact[1]

  # Graded
  # dat0 <- simdata(a0, matrix(c(d0, d0-.5*a0), n.items, 2), n.persons, "graded", Theta = matrix(x0))
  # dat1 <- simdata(a0, matrix(c(d1, d1-.5*a0), n.items, 2), n.persons, "graded", Theta = matrix(x1))

  # 2PL
  dat0 <- simdata(a0, d0, n.persons, '2PL', Theta = matrix(x0))
  dat1 <- simdata(a0, d1, n.persons, '2PL', Theta = matrix(x1))
  return(list(dat0 = as.data.frame(dat0), dat1 = as.data.frame(dat1)))
 }
