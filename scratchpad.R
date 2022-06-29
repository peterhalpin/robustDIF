

# weighted area between ICCs to get DIF.
# sum (abs(ICC1 * N1 - ICC2 * N2)).
#
# use delta = .5 and N1 = N2 as benchmark.
# then use the to find the equiv values for other choices of N2.
# So we are going to want a version of this we can solve for in terms of item parameters



# Testing out sims
library(mirt)

# Stuff
n.reps = 500
n.persons = 500
n.items = 15
n.biased = 5
bias = .5
impact = c(.5,1)

a.lower <- .8
a.upper <- 2
b.lim <- 1.5

  # Sim loop for parallelization via mclapply
  loop <- function(i){

    # DGP
    a0 <- runif(n.items, a.lower, a.upper)
    a1 <- a0 # apply_bias(a0, n.biased, bias, multiplicative = T)
    b0 <- sort(runif(n.items, -b.lim, b.lim))
    b1 <- apply_bias(b0, n.biased, bias)
    d0 <- b0*a0
    d1 <- b1*a0
    x0 <- rnorm(n.persons)
    x1 <- impact[2]*rnorm(n.persons) + impact[1]

    # Data gen
    dat0 <- simdata(a0, d0, n.persons, '2PL', Theta = matrix(x0))
    dat1 <- simdata(a1, d1, n.persons, '2PL', Theta = matrix(x1))

    # Fit IRT models
    fit0 <- mirt(dat0, 1, SE = T, SE.type = 'Oakes')
    fit1 <- mirt(dat1, 1, SE = T, SE.type = 'Oakes')

    # DIF Procedures
    irt.mle <- get_irt_mle(list(fit0, fit1))
    par(mfrow = c(1,2))
    rho.theta <- rho_fun(irt.mle)
    plot(rho.theta$theta, rho.theta$rho, type = "l")
    #get_starts(irt.mle)
    rho.sigma <- rho_fun(irt.mle, par = "slope", log = F)
    plot(rho.sigma$theta, rho.sigma$rho, type = "l")
    #get_starts(irt.mle, par = "slope", log = F)

    bsq.theta.out <- irls(irt.mle)
    true.theta.out <- bsq_weight(impact[1]/impact[2], y_fun(irt.mle), var_y(impact[1]/impact[2], irt.mle))
    bsq.sigma.out <- irls(irt.mle, par = "slope", log = log)
    true.sigma.out <- bsq_weight((impact[2]),
                                 y_fun(irt.mle, par = "slope", log = log),
                                 var_y(impact[2], irt.mle, par = "slope", log = log))

    z.theta <- z_test(bsq.theta.out$est, irt.mle)
    z.sigma <- z_test(bsq.sigma.out$est, irt.mle, par = "slope", log = log)
    chi2.test <- chi2_test(bsq.theta.out$est, bsq.sigma.out$est, irt.mle, log = log)
    lr.out <- lr(dat0, dat1)
    mh.out <- mh(dat0, dat1)
    lasso.out <- lasso(dat0, dat1)

    # Format output
    dif <- data.frame(a0 = a0,
                      d0 = d0,
                      dgp.a = a1-a0 != 0,
                      dgp.d = d1-d0 != 0,
                      rdif.theta.true = true.theta.out,
                      rdif.theta = bsq.theta.out$weights,
                      z.theta = z.theta$p.val < .05,
                      rdif.sigma.true = true.sigma.out,
                      rdif.sigma = bsq.sigma.out$weights,
                      z.sigma = z.sigma$p.val < .05,
                      chi2 = chi2.test$p.val < .05,
                      lr = lr.out$p < .05,
                      mh = mh.out$p < .05,
                      lasso = lasso.out$dif > 0
                      )

    scale <- data.frame(theta = bsq.theta.out$est,
                        sigma = bsq.sigma.out$est
    #                    lr = lr.out$scale.parms[1],
    #                    mh = mh.out$mu,
    #                    lasso = lasso.out$mu
                        )

    list(dif = dif, scale = scale) #, irt.mle = irt.mle)
  }
