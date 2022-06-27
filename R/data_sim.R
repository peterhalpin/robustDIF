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
#' Wrapper for mirt's implementation of the likelihood ratio (LR) test of DIF
#'
#' @param dat0 item response data for reference group
#' @param dat1 item response data for comparison group
#' @param parms which parameters to tests for DIF
#'
#' @return A named list with p values from the LR test and the scale parameters estimated in the model that assumes invariance
# -------------------------------------------------------------------

lr <- function(dat0, dat1, parms = c("d")){
  dat <- rbind(dat0, dat1)
  group <- c(rep("0", nrow(dat0)), rep("1", nrow(dat1)))
  mg.mod <- multipleGroup(dat,
              model = 1,
              group = group,
              invariance = c("intercepts",  "free_means"))
  lr.test <- DIF(mg.mod, which.par = parms, scheme = "drop", Wald = F)

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
                    control = ctrl_GPCMlasso(l.lambda = 20,
                                             adaptive = T,
                                             steptol = 1e-4,
                                             gradtol = 1e-4,
                                             iterlim = 50))
  mu <- coef(gpcm)[which.min(gpcm$BIC), ][(n.items+1)]
  dif <- abs(coef(gpcm)[which.min(gpcm$BIC), ][(n.items+2):(2*n.items+1)])
  list(dif = dif, mu = mu)
}


# -------------------------------------------------------------------
#' Runs a simulation study comparing R-DIF, the likelihood ratio test, the Mantel-Haenszel test, and the output of GPMClasso
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
  a.lower <- .8
  a.upper <- 2
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
    lr.out <- lr(dat0, dat1)
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
#' Post processing for \code{sim_study} to compute decision errors
#'
#' @param dif The DIF table output from  \code{sim_study}
#' @export
#' @return A data.frame with false and true positives for each method length
#'
# -------------------------------------------------------------------

decision_errors <-function(dif){
  test_names <- c("rdif.theta.true", "rdif.theta", "lr", "mh", "lasso")
  cuts <- c(1e-6, 1e-6, .05, .05, 1e-6)
  decisions <- data.frame(t(t(ds[,test_names]) < cuts))
  decisions$lasso <- decisions$lasso == F
  tp <- apply(decisions[ds$dgp == T, ], 2, mean)
  fp <- apply(decisions[ds$dgp == F, ], 2, mean)
  out <- data.frame(cbind(fp, tp))
  out$method <- rownames(out)
  rownames(out) <- NULL
  out
}
