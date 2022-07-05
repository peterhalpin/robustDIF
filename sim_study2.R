# Sim study 2 July 6, 2022
library(ggplot2)
library(dplyr)

### Data gen

n.reps = 500
n.items =  10

# sample size = 200
n.persons = 200
ds.a200  <- sim_study2(n.reps, n.persons, n.items, bias = c(.5, 0))
ds.b200  <- sim_study2(n.reps, n.persons, n.items, bias = c(0, 1))
ds.c200  <- sim_study2(n.reps, n.persons, n.items, bias = c(.35, .5))

# sample size = 250
n.persons = 350
ds.a350  <- sim_study2(n.reps, n.persons, n.items, bias = c(.5, 0))
ds.b350  <- sim_study2(n.reps, n.persons, n.items, bias = c(0, 1))
ds.c350  <- sim_study2(n.reps, n.persons, n.items, bias = c(.35, .5))

# sample size = 500
n.persons = 500
ds.a500  <- sim_study2(n.reps, n.persons, n.items, bias = c(.5, 0))
ds.b500  <- sim_study2(n.reps, n.persons, n.items, bias = c(0, 1))
ds.c500  <- sim_study2(n.reps, n.persons, n.items, bias = c(.35, .5))

sim.study2 <- list(ds.a200=ds.a200, ds.b200=ds.b200, ds.c200=ds.c200,
                   ds.a350=ds.a350, ds.b350=ds.b350, ds.c350=ds.c350,
                   ds.a500=ds.a500, ds.b500=ds.b500, ds.c500=ds.c500)


sim.study2.path <- "~/Dropbox/Academic/Manuscripts/DIF_via_scaling/data_analyses/sim2.july5.2022.RData"
 save(sim.study2, file = sim.study2.path)
# load(file = sim.study2.path)


 ### Plots

# Error rates
test.names <- names(sim.study2[[1]][[1]]$dif)[-c(1:3)] #,lr")
cut.offs <- c(1e-6, 1e-6, .05, 1e-6, 1e-6, .05, .05, .05, .05,.05, .05)

de <- function(x){
  decision_errors(x, test.names, cut.offs)
}

temp.dif <- lapply(sim.study2, function(x) {lapply(x, function(y) y$dif)})
ds.dif <- Reduce(rbind, lapply(temp.dif, function(x) de(Reduce(rbind, x))))
ds.dif$n.biased <- as.factor(rep(0:8, each = 6))
ds.dif.long <- ds.dif %>% tidyr::gather("fp", "tp", key = decision, value = Value)

temp.scale <- lapply(sim.study2, function(x) {lapply(x, function(y) y$scale)})
ds.scale <-  Reduce(rbind, lapply(temp.scale, function(x) Reduce(rbind, x)))
ds.scale$n.biased <- as.factor(rep(0:8, each = n.reps))


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


