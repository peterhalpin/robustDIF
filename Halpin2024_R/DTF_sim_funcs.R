##############################################################################
# Helper functions for DTF sim studies
# To be added to robustDIF pacakge for sim-based power analysis
##############################################################################

library(parallel)
library(mirt)

apply_dif <- function(n.items, n.dif, dif.effect.size) {
  out <- rep(0, n.items)
  dif.items <- sample(n.items, n.dif)
  out[dif.items] <- dif.effect.size
  out
}

data_sim <- function(
  n.reps = 500,
  n.persons = 500, 
  n.items = 20,
  n.dif = 0,
  dif.effect.size = 0.5) {
  
  # Item hyper-parms
  a.lower <- 0.5
  a.upper <- 2
  b.lower <- -1.5
  b.upper <- 1.5
  mu.lower <- -0.25
  mu.upper <- 0.25
  
  loop <- function(i) {
    a0 <- runif(n.items, a.lower, a.upper)
    b0 <- sort(runif(n.items, b.lower, b.upper)) #difficulties, not intercepts
    dif.vector <- apply_dif(n.items, n.dif, dif.effect.size)
    d0 <- -b0*a0
    d1 <- -(b0 - dif.vector)*a0
    #(d1 - d0) / a0
    # Generate latent trait (scaling due to model identification)
    mu <- runif(1, mu.lower, mu.upper)
    theta0 <- rnorm(n.persons, mean = 0, sd = 1)
    theta1 <- rnorm(n.persons, mean = mu, sd = 1) 
  
    # Genterate response data
    dat0 <- simdata(a0, d0, n.persons, '2PL', Theta = matrix(theta0))
    dat1 <- simdata(a0, d1, n.persons, '2PL', Theta = matrix(theta1))
  
    # Fit IRT models
    fit0 <- mirt(dat0, 1, SE = T)
    fit1 <- mirt(dat1, 1, SE = T)  
    
    # RDIF
    mle <- get_model_parms(list(fit0, fit1))
    rdif.delta <- delta_test(mle, fun = "d_fun3", alpha = .05)
     # LRT
    lrt.dif <- lrt_dif(dat0, dat1)$dif.items
    if (length(lrt.dif) > 0){
      lrt.dif.items <- which(colnames(dat0)%in%lrt.dif)
      lrt.delta <- delta_test_from_dif(
                      lrt.dif.items, 
                      mle, 
                      fun = "d_fun3", 
                      alpha = 0.05)
    } else {
      lrt.delta <- rep(NA, times = length(rdif.delta))
      names(lrt.delta) <- names(rdif.delta)
      
    }
    
    # Using true weights
    true.dif.items <- which(dif.vector > 0)
    if (length(true.dif.items) > 0){
      true.delta.hat <- delta_test_from_dif(
                          true.dif.items, 
                          mle, 
                          fun = "d_fun3", 
                          alpha = 0.05)
    } else {
      true.delta.hat <- rep(NA, times = length(rdif.delta))
      names(true.delta.hat) <- names(rdif.delta)
      
    }
    
    list(rdif = rdif.delta, 
         lrt = lrt.delta, 
         true.delta.hat = true.delta.hat, 
         true.scale = mu, 
         lrt.dif= lrt.dif)
  }  
  
  parallel::mclapply(1:n.reps, loop, mc.cores = detectCores()-1)
}  

lrt_dif <- function(dat0, dat1){
  dat.mg <- rbind(dat0, dat1)
  groups <- c(rep("g1", times = nrow(dat0)), 
              rep("g2", times = nrow(dat1))) 
  purify.const <- c("free_mean", "free_var", "slopes", "intercepts")
  fit.purify <- multipleGroup(
              dat.mg, 
              model = 1, 
              group = groups, 
              invariance = purify.const) 
  
  dif.purify <- DIF(fit.purify, 
             which.par = "d", 
             scheme = "drop") 
  dif.purify$p[is.nan(dif.purify$p)] <- 1
  anchor <- rownames(dif.purify)[dif.purify$p > .05]
  
  if(length(anchor) < ncol(dat.mg)){
     not.anchor <- colnames(dat.mg)[!colnames(dat.mg)%in%anchor]
     
     # New base model
     anchor.const <- c("free_mean", "free_var", "slopes", anchor)
     logL.baseline <- multipleGroup(
              dat.mg, 
              model = 1, 
              group = groups, 
              invariance = anchor.const)@Fit$logLik
         
     # Fit models that release constraints on each anchor item
     logL.anchor <-
       lapply(1:length(anchor),
          function(x){
            constraints <- c("free_mean", "free_var", "slopes", anchor[-x])
            multipleGroup(
              dat.mg,
              model = 1,
              group = groups,
              invariance = constraints)@Fit$logLik
          }
        )
     
     # Fit models that release constraints on each anchor item
     logL.not.anchor <- 
       lapply(1:length(not.anchor), 
          function(x){
            constraints <- c("free_mean", "free_var", "slopes", anchor, not.anchor[x])
            multipleGroup(
              dat.mg, 
              model = 1, 
              group = groups, 
              invariance = constraints)@Fit$logLik
          }
        )
      LR <- -2 * c((logL.baseline - unlist(logL.anchor)), 
                   (unlist(logL.not.anchor) - logL.baseline))
    
      LR[LR < 0] <- 0 
      p.vals <- 1 - pchisq(LR, 1)
      dif.table <- data.frame(
        items = c(anchor, not.anchor), 
        chi.square = LR, 
        p.vals = p.vals)
      list(dif.items = c(anchor, not.anchor)[
              p.adjust(p.vals, method = "BH") < .05], 
           dif.table = dif.table[order(dif.table$items), ])
  } else {
      list(dif.items = NULL, dif.table = NULL)
  }
}
