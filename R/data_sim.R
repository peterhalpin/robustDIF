require(mirt)
require(GPCMlasso)
require(difR)

# -------------------------------------------------------------------
#' Apply bias to a random subset of item intercepts
#'
#' Applies additive bias (parm + bias) or multiplicative bias (parm + bias*parm) to a vector of item parameters. Use additive for item difficulties (not intercepts) and multiplicative for item slopes.
#'
#' @param b0 Vector of item parameters in the reference group
#' @param n.biased Number of biased items
#' @param max.bias Maximum amount of bias to apply
#' @param worst.case Logical: apply worst case bias? Default is \code{TRUE}. If set to \code{FALSE}, bias is sampled from \code{runif(max.bias/2, max.bias)}
#' @param type String: one of \code{c("add", "mult")}. See comments
#'
#' @return A vector of item intercepts, with bias added
#'
# -------------------------------------------------------------------

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

# -------------------------------------------------------------------
#' Wrapper for mirt's implementation of the likelihood ratio (LR) test of DIF for item intercept only
#'
#' @param dat0 item response data for reference group
#' @param dat1 item response data for comparison group
#'
#' @return A named list with p values from the LR test and the scale parameters estimated in the model that assumes invariance
# -------------------------------------------------------------------

lr_z <- function(dat0, dat1, which.par = "d"){
  dat <- rbind(dat0, dat1)
  group <- c(rep("0", nrow(dat0)), rep("1", nrow(dat1)))
  if (which.par == "d") {invariance = c("intercepts", "free_means")}
  if (which.par == "a1") {invariance = c("slopes", "free_vars")}

  mg.mod <- multipleGroup(dat,
              model = 1,
              group = group,
              invariance = invariance)
  lr.test <- DIF(mg.mod, which.par = which.par, scheme = "drop", Wald = F)

  #p2 <- DIF(mg.mod, which.par = parms, scheme = "drop_sequential",
  #          seq_stat = .05, max_run = 2, p.adjust = "BH", Wald = F)

  list(p = lr.test$p, scale.parms = coef(mg.mod)$`1`$GroupPars)
}

# -------------------------------------------------------------------
#' Wrapper for mirt's implementation of the likelihood ratio (LR) test of DIF for item slopes and intercepts.
#'
#' @param dat0 item response data for reference group
#' @param dat1 item response data for comparison group
#'
#' @return A named list with p values from the LR test and the scale parameters estimated in the model that assumes invariance
# -------------------------------------------------------------------

lr_chi <- function(dat0, dat1){
  dat <- rbind(dat0, dat1)
  group <- c(rep("0", nrow(dat0)), rep("1", nrow(dat1)))
  mg.mod <- multipleGroup(dat,
              model = 1,
              group = group,
              invariance = c("intercepts", "slopes", "free_means", "free_vars"))
  lr.test <- DIF(mg.mod, which.par = c("a1", "d"), scheme = "drop", Wald = F)

  #p2 <- DIF(mg.mod, which.par = parms, scheme = "drop_sequential",
  #          seq_stat = .05, max_run = 2, p.adjust = "BH", Wald = F)

  list(p = lr.test$p, scale.parms = coef(mg.mod)$`1`$GroupPars)
}


# -------------------------------------------------------------------
#' Wrapper for difR's implementation of the Mantel-Haenszel test of (uniform) DIF
#'
#' @param dat0 item response data for reference group
#' @param dat1 item response data for comparison group
#'
#' @return A named list with p values from the MH test, and the scale parameter estimated from the total scores
#'
# -------------------------------------------------------------------

mh <- function(dat0, dat1){
  dat <- rbind(dat0, dat1)
  group <- c(rep("0", nrow(dat0)), rep("1", nrow(dat1)))
  p <- difR::difMH(dat, group, focal.name = 0)$p.value
  mu <- mean(apply(dat1, 1, sum) - apply(dat0, 1, sum)) / sd(apply(dat1, 1, sum))
  list(p = p, mu = mu)
}

# -------------------------------------------------------------------
#' Wrapper for GPMClasso
#'
#' @param dat0 item response data for reference group
#' @param dat1 item response data for comparison group
#'
#' @return A named list with dif parameters for the item difficulties, and the estimated mean of .
#'
# -------------------------------------------------------------------

lasso <- function(dat0, dat1) {
  n.items <- ncol(dat0)
  dat <- data.frame(rbind(dat0, dat1))
  group <- c(rep(0, nrow(dat0)), rep(1, nrow(dat1)))
  dat.g <- data.frame(dat, group)
  mod <- as.formula(paste("cbind(",
                    paste(colnames(dat.g)[1:n.items],collapse=","),") ~ group"))
  gpcm <- GPCMlasso(mod, dat.g, model = "2PL",
                    control = ctrl_GPCMlasso(l.lambda = 25,
                                             #adaptive = T,
                                             steptol = 1e-4,
                                             gradtol = 1e-4,
                                             iterlim = 50))

  solution.index <- which.min(gpcm$BIC)
  mu <- coef(gpcm)[solution.index, ][(n.items+1)]
  dif <- abs(coef(gpcm)[solution.index, ][(n.items+2):(2*n.items+1)])
  list(dif = dif, mu = mu)
}

# -------------------------------------------------------------------
#' Runs a simulation study comparing R-DIF, the likelihood ratio test, the Mantel-Haenszel test, and  GPMClasso
#'
#' @param n.reps number of replications
#' @param n.persons number of respondents per group
#' @param n.item number of items
#' @param n.biased number of items with DIF on item intercepts
#' @param bias magnitude of bias (the same for all items with bias)
#' @param impact the mean and variance of the latent trait in the reference groups
#' @export
#' @return A n.reps length list, where each element is a list containing (a) a data.frame with the output of the DIF procedures, and (b) a vector with estimates of impact from the procedures.
#'
# -------------------------------------------------------------------

sim_study1 <- function(n.reps = 100, n.persons = 500, n.items = 15, n.biased = 0, bias = 0, impact = c(0, 1)){

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

    # Data gen
    dat0 <- simdata(a0, d0, n.persons, '2PL', Theta = matrix(x0))
    dat1 <- simdata(a0, d1, n.persons, '2PL', Theta = matrix(x1))

    # Fit IRT models
    fit0 <- mirt(dat0, 1, SE = T, SE.type = 'Oakes')
    fit1 <- mirt(dat1, 1, SE = T, SE.type = 'Oakes')

    # DIF Procedures
    irt.mle <- get_irt_mle(list(fit0, fit1))
    bsq.theta.out <- irls(irt.mle)
    true.theta.out <- bsq_weight(impact[1], y_fun(irt.mle), var_y(impact[1], irt.mle))
    bsq.sigma.out <- irls(irt.mle, par = "slope")
    chi2.test <- chi2_test(bsq.theta.out$est, bsq.sigma.out$est, irt.mle)
    lr.out <- lr_intercept(dat0, dat1)
    mh.out <- mh(dat0, dat1)
    lasso.out <- lasso(dat0, dat1)

    # Format output
    dif <- data.frame(a0 = a0,
                      d0 = d0,
                      dgp = d1-d0 != 0,
                      rdif.true = true.theta.out,
                      rdif.flag = bsq.theta.out$weights,
                      rdif.chi2 = chi2.test$p.val,
                      lr = lr.out$p,
                      mh = mh.out$p,
                      lasso = lasso.out$dif
                      )

    scale <- data.frame(theta = bsq.theta.out$est,
                        sigma = bsq.sigma.out$est
    #                    lr = lr.out$scale.parms[1],
    #                    mh = mh.out$mu,
    #                    lasso = lasso.out$mu
                        )

    list(dif = dif, scale = scale) #, irt.mle = irt.mle)
  }
  parallel::mclapply(1:n.reps, loop)
}


# -------------------------------------------------------------------
#' Runs a simulation study comparing R-DIF and the likelihood ratio test
#'
#' @param n.reps number of replications
#' @param n.persons number of respondents per group
#' @param n.item number of items
#' @param n.biased number of items with DIF on item intercepts
#' @param bias magnitude of bias (the same for all items with bias)
#' @param impact the mean and variance of the latent trait in the reference groups
#' @export
#' @return A n.reps length list, where each element is a list containing (a) a data.frame with the output of the DIF procedures, and (b) a vector with estimates of impact from the procedures.
#'
# -------------------------------------------------------------------

sim_study2 <- function(n.reps = 100, n.persons = 500, n.items = 15, bias = c(.5, 1)){

  # Item hyper-parms
  a.lower <- .9
  a.upper <- 2.5
  b.lim <- 1.5
  mu.lim <- .5
  sigma.lower <- sqrt(.5)
  sigma.upper <- sqrt(2)

  # Sim loop for parallelization via mclapply
  loop <- function(i){

    # DGP
    mu <- runif(1, -mu.lim, mu.lim)
    sigma <- runif(1, sigma.lower, sigma.upper)
    theta <- mu/sigma
    biased.item <- sample(n.items, 1)

    # Bias on slope is proportional
    a0 <- a1 <- runif(n.items, a.lower, a.upper)
    a1[biased.item] <-  a1[biased.item] +  bias[2] * a1[biased.item]

    # Bias on intercept is additive; make y_intercept = bias[1] even if bias[2] != 0 ?
    b0 <- b1 <- sort(runif(n.items, -b.lim, b.lim))
    b1[biased.item] <-  b1[biased.item] + bias[1]+ bias[1] * bias[2]
    d0 <- a0*b0
    d1 <- a0*b1
    # (d1 - d0) / a1

    x0 <- rnorm(n.persons)
    x1 <- sigma * rnorm(n.persons) + mu

    # Data gen
    dat0 <- simdata(a0, d0, n.persons, '2PL', Theta = matrix(x0))
    dat1 <- simdata(a1, d1, n.persons, '2PL', Theta = matrix(x1))

    # Fit IRT models
    fit0 <- mirt(dat0, 1, SE = T, SE.type = 'Oakes')
    fit1 <- mirt(dat1, 1, SE = T, SE.type = 'Oakes')

    # DIF Procedures
    irt.mle <- get_irt_mle(list(fit0, fit1))

    # par(mfrow = c(1,2))
    # rho.theta <- rho_fun(irt.mle)
    # plot(rho.theta$theta, rho.theta$rho)
    #
    # rho.sigma <- rho_fun(irt.mle, par = "slope")
    # plot(rho.sigma$theta, rho.sigma$rho)

    true.theta.out <- bsq_weight(theta, y_fun(irt.mle),
                                 var_y(theta, irt.mle))
    rdif.theta.out <- irls(irt.mle)
    z.test.theta <- z_test(rdif.theta.out$est, irt.mle)

    true.sigma.out <- bsq_weight(sigma, y_fun(irt.mle, par = "slope"),
                                 var_y(sigma, irt.mle, par = "slope"))
    rdif.sigma.out <- irls(irt.mle, par = "slope")
    z.test.sigma <- z_test(rdif.sigma.out$est, irt.mle , par = "slope")

    chi2 <- chi2_test(rdif.theta.out$est, rdif.sigma.out$est, irt.mle)
    chi2.true <- chi2_test(theta, sigma, irt.mle)


    lr.int <- lr_z(dat0, dat1, which.par = "d")
    lr.slope <- lr_z(dat0, dat1, which.par = "a1")
    lr.chi <- lr(dat0, dat1)

    # Format output
    dif <- data.frame(a0 = a0,
                      d0 = d0,
                      dgp = d1-d0 != 0 | a1-a0 != 0 ,
                      theta.true = true.theta.out,
                      theta.flag = rdif.theta.out$weights,
                      theta.test = z.test.theta$p.val,
                      sigma.true = true.theta.out,
                      sigma.flag = rdif.sigma.out$weights,
                      sigma.test = z.test.sigma$p.val,
                      rdif.chi2.true = chi2.true$p.val,
                      rdif.chi2 = chi2$p.val,
                      lr.intercept = lr.int$p,
                      lr.slope = lr.slope$p,
                      lr.chi = lr.chi$p
                      )

    scale <- data.frame(theta = theta,
                        sigma = sigma,
                        rdif.theta = rdif.theta.out$est,
                        rdif.sigma = rdif.sigma.out$est
    #                    lr = lr.out$scale.parms[1],
    #                    mh = mh.out$mu,
    #                    lasso = lasso.out$mu
                        )

    list(dif = dif, scale = scale) #, irt.mle = irt.mle)
  }
  parallel::mclapply(1:n.reps, loop)
}

# -------------------------------------------------------------------
#' Post processing for sim studies, to compute decision errors
#'
#' @param dif The DIF table output from one for the \code{robustDiF} simulation study functions.
#' @param test.names The col.names \code{dif} that correspond to DIF method.
#' @param cut.offs  Values used to make a decision about each DIF method in \code{test.names}. Can a single value or a vector of same length as \code{test.names}.
#'
#' @return A data.frame with false and true positives for each method length
#'
# -------------------------------------------------------------------

decision_errors <-function(ds, test.names, cut.offs){
  decisions <- data.frame(t(t(ds[,test.names]) < cut.offs))
  if(length(grep("lasso", test.names)) > 0) {
    decisions$lasso <- decisions$lasso == F
  }
  tp <- apply(decisions[ds$dgp == T, ], 2, mean)
  fp <- apply(decisions[ds$dgp == F, ], 2, mean)
  out <- data.frame(cbind(fp, tp))
  out$method <- rownames(out)
  rownames(out) <- NULL
  out
}
