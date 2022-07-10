# Sim study 2 June 30, 2022

ds1  <- sim_study2(n.reps, 500, n.items = 10, bias = c(.5, 0))
ds2  <- sim_study2(n.reps, 500, n.items = 10, bias = c(0, 1))
ds3  <- sim_study2(n.reps, 500, n.items = 10, bias = c(.35, .5))
ds.prev <- ds
ds <- list(ds1 = ds1, ds2 = ds2, ds3 = ds3)

# sim.study2.path <- "~/Dropbox/Academic/Manuscripts/DIF_via_scaling/data_analyses/sim1.june28.2022.RData"
# save(sim.study2, file = sim.study2.path)
# load(file = sim.study2.path)
# Pullng everything out...sheesh

test.names <- names(ds[[1]][[1]]$dif)[-c(1:3)] #,lr")
cut.offs <- c(1e-6, 1e-6, .05, 1e-6, 1e-6, .05, .05, .05, .05, .05, .05)

de <- function(x){
  decision_errors(x, test.names, cut.offs)
}

# not trying to correct effect sizes
temp.dif <- lapply(ds, function(x) {lapply(x, function(y) y$dif)})
Reduce(rbind, lapply(temp.dif, function(x) de(Reduce(rbind, x))))

# large / large
ds.ll200  <- sim_study2(n.reps, 200, n.items, bias = c(.5, 1))
ds.ll350  <- sim_study2(n.reps, 350, n.items, bias = c(.5, 1))
ds.ll500  <- sim_study2(n.reps, 500, n.items, bias = c(.5, 1))

# large / small
ds.ls200  <- sim_study2(n.reps, 200, n.items, bias = c(.5, .5))
ds.ls350  <- sim_study2(n.reps, 350, n.items, bias = c(.5, .5))
ds.ls500  <- sim_study2(n.reps, 500, n.items, bias = c(.5, .5))

# small / large
ds.sl200  <- sim_study2(n.reps, 200, n.items, bias = c(.25, 1))
ds.sl350  <- sim_study2(n.reps, 350, n.items, bias = c(.25, 1))
ds.sl500  <- sim_study2(n.reps, 500, n.items, bias = c(.25, 1))

# small / small
ds.ss200  <- sim_study2(n.reps, 200, n.items, bias = c(.25, .5))
ds.ss350  <- sim_study2(n.reps, 350, n.items, bias = c(.25, .5))
ds.ss500  <- sim_study2(n.reps, 500, n.items, bias = c(.25, .5))

ds <- list(ds.ll500 = ds.ll500)
ds1 <- ds
# sim.study2.path <- "~/Dropbox/Academic/Manuscripts/DIF_via_scaling/data_analyses/sim1.june28.2022.RData"
# save(sim.study2, file = sim.study2.path)
# load(file = sim.study2.path)
# Pullng everything out...sheesh

test.names <- names(ds[[1]][[1]]$dif)[-c(1:3)] #,lr")
cut.offs <- c(1e-6, 1e-6, .05, 1e-6, 1e-6, .05, .05, .05)

de <- function(x){
  decision_errors(x, test.names, cut.offs)
}


temp.dif <- lapply(ds, function(x) {lapply(x, function(y) y$dif)})
ds.dif <- Reduce(rbind, lapply(temp.dif, function(x) de(Reduce(rbind, x))))
ds.dif$n.biased <- as.factor(rep(0:8, each = 6))
ds.dif.long <- ds.dif %>% tidyr::gather("fp", "tp", key = decision, value = Value)

temp.scale <- lapply(sim.study2, function(x) {lapply(x, function(y) y$scale)})
ds.scale <-  Reduce(rbind, lapply(temp.scale, function(x) Reduce(rbind, x)))
ds.scale$n.biased <- as.factor(rep(0:8, each = n.reps))

### Plots

# Error rates
library(ggplot2)
ds.dif.long <- ds.dif.long[ds.dif.long$method != "rdif.chi2", ]
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
            theme(text = element_text(size=20)) +
            scale_colour_brewer(palette = "Paired")
p1 + facet_wrap(~ decision, nrow = 2)

# Scale distribution
facet.labels <- paste0("N.DIF: ", 0:8,  " out of 15")
facet_labeller <- function(variable, value){
  return(facet.labels[value])
}

p2 <- ggplot(ds.scale, aes(x = theta)) +
            geom_histogram(col = "white", fill = 'grey65') +
            ylab("Count") +
            xlab("R-DIF estimate of the IRT scale parameter") +
            theme(text = element_text(size=20))
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


