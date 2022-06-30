# Sim study 1 June 28-29, 2022
library(ggplot2)
library(dplyr)

### Data gen

n.reps = 500
n.persons = 500
n.items = 15
bias = .5
impact = c(.5,1)

ds0 <- sim_study1(n.reps, n.persons, n.items, n.biased = 0, bias, impact)
ds1 <- sim_study1(n.reps, n.persons, n.items, n.biased = 1, bias, impact)
ds2 <- sim_study1(n.reps, n.persons, n.items, n.biased = 2, bias, impact)
ds3 <- sim_study1(n.reps, n.persons, n.items, n.biased = 3, bias, impact)
ds4 <- sim_study1(n.reps, n.persons, n.items, n.biased = 4, bias, impact)
ds5 <- sim_study1(n.reps, n.persons, n.items, n.biased = 5, bias, impact)
ds6 <- sim_study1(n.reps, n.persons, n.items, n.biased = 6, bias, impact)
ds7 <- sim_study1(n.reps, n.persons, n.items, n.biased = 7, bias, impact)
ds8 <- sim_study1(n.reps, n.persons, n.items, n.biased = 8, bias, impact)

sim.study1 <- list(ds0=ds0, ds1=ds1, ds2=ds2, ds3=ds3, ds4=ds4, ds5=ds5, ds6=ds6, ds7=ds7, ds8=ds8)
sim.study1.path <- "~/Dropbox/Academic/Manuscripts/DIF_via_scaling/data_analyses/sim1.june29.2022.RData"
#save(sim.study1, file = sim.study1.path)
load(file = sim.study1.path)



### Plots

# Error rates
test.names <- c("rdif.true", "rdif.flag", "rdif.chi2", "lr", "mh", "lasso")
cut.offs <- c(1e-6, 1e-6, .05, .05, .05, 1e-6)
d_e <- function(x){
  decision_errors(x, test.names, cut.offs)
}

temp.dif <- lapply(sim.study1, function(x) {lapply(x, function(y) y$dif)})
ds.dif <- Reduce(rbind, lapply(temp.dif, function(x) d_e(Reduce(rbind, x))))
ds.dif$n.biased <- as.factor(rep(0:8, each = 6))
ds.dif.long <- ds.dif %>% tidyr::gather("fp", "tp", key = decision, value = Value)

#ds.dif.long <- ds.dif.long[ds.dif.long$method != "rdif.chi2", ] # Not necessary
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
temp.scale <- lapply(sim.study2, function(x) {lapply(x, function(y) y$scale)})
ds.scale <-  Reduce(rbind, lapply(temp.scale, function(x) Reduce(rbind, x)))
ds.scale$n.biased <- as.factor(rep(0:8, each = n.reps))

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


