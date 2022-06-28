# Sim study 1 June 27 - 28, 2022

n.reps = 500
n.persons = 500
n.items = 15
bias = .5
n.biased = 7
impact = c(.5,1)


  # Sim loop for parallelization via mclapply
#   loop <- function(i){
#
#     # DGP
#     a0 <- runif(n.items, a.lower, a.upper)
#     b0 <- sort(runif(n.items, -b.lim, b.lim))
#     b1 <- apply_bias(b0, n.biased, bias)
#     d0 <- b0*a0
#     d1 <- b1*a0
#     x0 <- rnorm(n.persons)
#     x1 <- impact[2]*rnorm(n.persons) + impact[1]
#
#     # Data gen
#     dat0 <- simdata(a0, d0, n.persons, '2PL', Theta = matrix(x0))
#     dat1 <- simdata(a0, d1, n.persons, '2PL', Theta = matrix(x1))
#
#     # Fit IRT models
#     fit0 <- mirt(dat0, 1, SE = T, SE.type = 'Oakes')
#     fit1 <- mirt(dat1, 1, SE = T, SE.type = 'Oakes')
#
#     # DIF Procedures
#     irt.mle <- get_irt_mle(list(fit0, fit1))
#     bsq.theta.out <- irls(irt.mle)
#     true.theta.out <- bsq_weight(impact[1], y_fun(irt.mle), var_y(impact[1], irt.mle))
#     bsq.sigma.out <- irls(irt.mle, par = "slope")
#     chi2.test <- chi2_test(bsq.theta.out$est, bsq.sigma.out$est, irt.mle)
#     #lr.out <- lr(dat0, dat1)
#     #mh.out <- mh(dat0, dat1)
#     #lasso.out <- lasso(dat0, dat1)
#
#     # Format output
#     dif <- data.frame(a0 = a0,
#                       d0 = d0,
#                       dgp = d1-d0 != 0,
#                       rdif.true = true.theta.out,
#                       rdif.flag = bsq.theta.out$weights,
#                       rdif.chi2 = chi2.test$p.val
#                       #lr = lr.out$p,
#                       #mh = mh.out$p,
#                       #lasso = lasso.out$dif > 0
#                       )
#
#     scale <- data.frame(theta = bsq.theta.out$est,
#                         sigma = bsq.sigma.out$est
#     #                    lr = lr.out$scale.parms[1],
#     #                    mh = mh.out$mu,
#     #                    lasso = lasso.out$mu
#                         )
#
#     list(dif = dif, scale = scale)
# }

ds0 <- sim_study1(n.reps, n.persons, n.items, n.biased = 0, bias, impact)
ds1 <- sim_study1(n.reps, n.persons, n.items, n.biased = 1, bias, impact)
ds2 <- sim_study1(n.reps, n.persons, n.items, n.biased = 2, bias, impact)
ds3 <- sim_study1(n.reps, n.persons, n.items, n.biased = 3, bias, impact)
ds4 <- sim_study1(n.reps, n.persons, n.items, n.biased = 4, bias, impact)
ds5 <- sim_study1(n.reps, n.persons, n.items, n.biased = 5, bias, impact)
ds6 <- sim_study1(n.reps, n.persons, n.items, n.biased = 6, bias, impact)
ds7 <- sim_study1(n.reps, n.persons, n.items, n.biased = 7, bias, impact)
ds8 <- sim_study1(n.reps, n.persons, n.items, n.biased = 8, bias, impact)

sim.study2 <- list(ds0=ds0, ds1=ds1, ds2=ds2, ds3=ds3, ds4=ds4, ds5=ds5, ds6=ds6, ds7=ds7, ds8=ds8)

sim.study2.path <- "~/Dropbox/Academic/Manuscripts/DIF_via_scaling/data_analyses/sim1.june28.2022.RData"
save(sim.study2, file = sim.study2.path)

load(file = sim.study2.path)

# Pullng everything out...sheesh
temp.dif <- lapply(sim.study1, function(x) {lapply(x, function(y) y$dif)})
ds.dif <- Reduce(rbind, lapply(temp.dif, function(x) decision_errors(Reduce(rbind, x))))
ds.dif$n.biased <- as.factor(rep(0:8, each = 6))
ds.dif.long <- ds.dif %>% tidyr::gather("fp", "tp", key = decision, value = Value)

temp.scale <- lapply(ds, function(x) {lapply(x, function(y) y$scale)})
ds.scale <-  Reduce(rbind, lapply(temp.scale, function(x) Reduce(rbind, x)))
ds.scale$n.biased <- as.factor(rep(0:8, each = n.reps))
ds.scale.long <- ds.scale %>% gather("bsq", "lr", "mh", "lasso", key = Method, value = Value)

### Plots

# Error rates
ds.dif.long$method <- ordered(ds.dif.long$method, unique(ds.dif.long$method))
ds.dif.long$decision[ds.dif.long$decision == "fp"] <- "False positive rate"
ds.dif.long$decision[ds.dif.long$decision == "tp"] <- "True positive rate"
names(ds.dif.long)[1] <- "Method"
p1 <- ggplot(ds.dif.long, aes(y = Value, x = n.biased, group = Method)) +
            geom_point(aes(color = Method), size = 2.5) +
            geom_line(aes(color = Method), size = 1, linetype = 1) +
            geom_hline(yintercept = .85, col = 'grey25', linetype = 2) +
            geom_hline(yintercept = .05, col = 'grey65', linetype = 2) +
            ylab("Value") +
            xlab("Number of biased items (out of 15)") +
            theme(text = element_text(size=15)) +
            scale_colour_brewer(palette = "Paired")
p1 + facet_wrap(~ decision, nrow = 2)

# Scale bias
my_colors <- RColorBrewer::brewer.pal(5, "Paired")[2:5]
ds.scale.long$Method <- ordered(ds.scale.long$Method, unique(ds.scale.long$Method))
ggplot(ds.scale.long, aes(y = Value, x = n.biased)) +
            geom_boxplot(aes(fill = Method)) +
            geom_hline(yintercept = .5, col = 1, linetype = 2) +
            ylab("IRT scale (mu)") +
            xlab("Number of biased items (out of 15)") +
            theme(text = element_text(size=15)) +
            scale_fill_manual(values = my_colors,
                              labels = c("bsq", "lr", "mh", "lasso"))


# Scale distribution
facet.labels <- paste0("N.biased = ", 0:8,  " out of 15")
facet_labeller <- function(variable, value){
  return(facet.labels[value])
}

p2 <- ggplot(ds.scale.long[ds.scale.long$Method == "bsq", ], aes(x = Value)) +
            geom_histogram(col = "white", fill = my_colors[1]) +
            ylab("Count") +
            xlab("IRT scale (mu)") +
            theme(text = element_text(size=15))
p2 + facet_wrap(~ n.biased, nrow = 3, labeller = facet_labeller)


# per sample fp rate plotted as a function of the scale parm
bsq.scale <- ds.scale.long[ds.scale.long$Method == "bsq", ]
temp.dif2 <- lapply(temp.dif, function(x) {lapply(x, function(y) mean(y$bsq[y$dgp == F] < .00001))})
bsq.dif <- unlist(temp.dif2)
length(bsq.dif)
dim(bsq.scale)
bsq.scale$dif <- bsq.dif
bsq.dif <- Reduce(c, lapply(temp.dif2, function (x) Reduce(c, x)))
tapply(bsq.scale$Value, bsq.scale$n.biased, mean)

plot(bsq.scale$Value, bsq.scale$dif)


