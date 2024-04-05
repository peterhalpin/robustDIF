##############################################################################
# Numerical examples for: Halpin (2022) Differential Item Functioning
# Via Robust Scaling.
##############################################################################

### Required packages
library(mirt)
library(GPCMlasso)
library(difR)
library(dplyr)
library(ggplot2)
library(ggpattern)
library(xtable)

## Version of the R-DIF package used
source("robustDIF_functions_Halpin2022.R")

# Helper functions for running simulations
source("data_sim_functions.R")

##############################################################################
# Sim Study 1: Breakdown
##############################################################################

# Data gen
# n.reps = 500
# n.persons = 500
# n.items = 16
# bias = .5
# impact = c(.5,1)
#
# ds0 <- sim_study1(n.reps, n.persons, n.items, n.biased = 0, bias, impact)
# ds1 <- sim_study1(n.reps, n.persons, n.items, n.biased = 1, bias, impact)
# ds2 <- sim_study1(n.reps, n.persons, n.items, n.biased = 2, bias, impact)
# ds3 <- sim_study1(n.reps, n.persons, n.items, n.biased = 3, bias, impact)
# ds4 <- sim_study1(n.reps, n.persons, n.items, n.biased = 4, bias, impact)
# ds5 <- sim_study1(n.reps, n.persons, n.items, n.biased = 5, bias, impact)
# ds6 <- sim_study1(n.reps, n.persons, n.items, n.biased = 6, bias, impact)
# ds7 <- sim_study1(n.reps, n.persons, n.items, n.biased = 7, bias, impact)
# ds8 <- sim_study1(n.reps, n.persons, n.items, n.biased = 8, bias, impact)
#
# sim.study1 <- list(ds0=ds0, ds1=ds1, ds2=ds2,
#                    ds3=ds3, ds4=ds4, ds5=ds5,
#                    ds6=ds6, ds7=ds7, ds8=ds8)
# sim.study1.path <- "sim1.Nov7.2023.RData"
# save(sim.study1, file = sim.study1.path)

# Data processing
sim.study1.path <- "sim1.Nov7.2023.RData"
load(file = sim.study1.path)

test.names <- c("rdif.true", "rdif.flag", "rdif.chi2", "lr", "mh", "lasso")
cut.offs <- c(1e-6, 1e-6, .05, .05, .05, 1e-6)
d_e <- function(x){
  decision_errors(x, test.names, cut.offs)
}

temp.dif <- lapply(sim.study1, function(x) {lapply(x, function(y) y$dif)})
ds.dif <- Reduce(rbind, lapply(temp.dif, function(x) d_e(Reduce(rbind, x))))
ds.dif$n.biased <- as.factor(rep(0:8, each = 6))
ds.dif.long <- ds.dif %>% tidyr::gather("fp", "tp", key = decision, value = Value)

ds.dif.long <- ds.dif.long[ds.dif.long$method != "rdif.chi2", ] # Omitted
ds.dif.long$method <- as.factor(ds.dif.long$method)
levels(ds.dif.long$method) <- c("Lasso", "LR", "MH", "RDIF.flag", "RDIF.true")
ds.dif.long$method <- ordered(ds.dif.long$method, unique(ds.dif.long$method))
ds.dif.long$decision[ds.dif.long$decision == "fp"] <- "False Positive Rate"
ds.dif.long$decision[ds.dif.long$decision == "tp"] <- "True Positive Rate"
names(ds.dif.long)[1] <- "Method"

# Figure 1 of Paper
p1 <- ggplot(ds.dif.long, aes(y = Value, x = n.biased, group = Method)) +
            geom_hline(yintercept = 0, col = 'grey90', linetype = 2) +
            geom_point(aes(color = Method, shape = Method, fill = Method),
                       size = 2.5) +
            geom_line(aes(color = Method), size = .5, linetype = 1) +
            ylab("Value") +
            xlab("Number of biased items (out of 16)") +
            scale_shape_manual(values=c(21:25)) +
            scale_colour_brewer(palette = "Paired") +
            scale_fill_brewer(palette = "Paired") +
            #scale_color_grey() +
            #scale_fill_grey() +
            theme_bw(base_size =  16)

p1 + facet_wrap(~ decision, nrow = 2, scales = "free")

# Figure 2 of paper
n.reps <- 500
temp.scale <- lapply(sim.study1, function(x) {lapply(x, function(y) y$scale)})
ds.scale <-  Reduce(rbind, lapply(temp.scale, function(x) Reduce(rbind, x)))
ds.scale$n.biased <- as.factor(rep(0:8, each = n.reps))

facet.labels <- paste0("N.DIF: ", 0:8,  " out of 16")
facet_labeller <- function(variable, value){
  return(facet.labels[value])
}

p2 <- ggplot(ds.scale, aes(x = theta)) +
            geom_histogram(col = "white", fill = "#A6CEE3") +
            #geom_histogram(col = "white", fill = 'grey65') +
            theme_bw(base_size =  16) +
            ylab("Count") +
            xlab("R-DIF estimate of the IRT scale parameter")

p2 + facet_wrap(~ n.biased, nrow = 3, labeller = facet_labeller)

##############################################################################
# Sim study 2: Power
##############################################################################

# Data gen
n.reps = 500
n.items =  10

# sample size = 200
n.persons = 200
ds.a200  <- sim_study2(n.reps, n.persons, n.items, bias = c(.5, 0))
ds.b200  <- sim_study2(n.reps, n.persons, n.items, bias = c(0, 1))
ds.c200  <- sim_study2(250, n.persons, n.items, bias = c(.35, .5))

# sample size = 350
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

sim.study2.path <- "sim2.july5.2022.RData"
#save(sim.study2, file = sim.study2.path)

# Data processing
load(file = sim.study2.path)

# Compute decision errors
test.names <- names(sim.study2[[1]][[1]]$dif)[-c(1:3)] #,lr")
cut.offs <- c(1e-6, 1e-6, .05, 1e-6, 1e-6, .05, .05, .05, .05,.05, .05)

de <- function(x){
  decision_errors(x, test.names, cut.offs)
}

temp.dif <- lapply(sim.study2, function(x) {lapply(x, function(y) y$dif)})
names(temp.dif)
ds.dif <- Reduce(rbind, lapply(temp.dif, function(x) de(Reduce(rbind, x))))
ds.dif$n <- rep(c(200, 350, 500), each = length(test.names)*3)
ds.dif$type <- rep(c("Intercept Only", "Slope Only", "Intercept + Slope"), each = length(test.names))

# Averaging error rates for non-tested parameters
slope.methods <- unique(ds.dif$method)[c(4:6, 10)]
temp.fp <- ds.dif[ds.dif$type =="Intercept Only" & ds.dif$method%in%slope.methods, ]
temp.fp$fp <- temp.fp$fp * .9 + temp.fp$tp * .1
temp.fp$tp <- 0
temp.fp -> ds.dif[ds.dif$type =="Intercept Only" & ds.dif$method%in%slope.methods, ]

int.methods <- unique(ds.dif$method)[c(1:3, 9)]
temp.fp <- ds.dif[ds.dif$type =="Slope Only" & ds.dif$method%in%int.methods, ]
temp.fp$fp <- temp.fp$fp * .9 + temp.fp$tp * .1
temp.fp$tp <- 0
temp.fp -> ds.dif[ds.dif$type =="Slope Only" & ds.dif$method%in%int.methods, ]

# Reformat for plots
ds.dif.plot <- ds.dif %>% tidyr::gather("fp", "tp", key = decision, value = Value)
omit.methods <- unique(ds.dif$method)[c(1, 3, 4, 6, 7)]
ds.dif.plot <- ds.dif.plot[!ds.dif.plot$method%in%omit.methods, ]
ds.dif.plot$method[ds.dif.plot$method == "theta.flag"] <- "rdif.intercept"
ds.dif.plot$method[ds.dif.plot$method == "sigma.flag"] <- "rdif.slope"
ds.dif.plot$method[ds.dif.plot$method == "lr.chi"] <- "lr.chi2"
ds.dif.plot$Test <- "intercept"
ds.dif.plot$Test[ds.dif.plot$method%in%c("lr.chi2", "rdif.chi2")] <- "both"
ds.dif.plot$Test[ds.dif.plot$method%in%c("lr.slope", "rdif.slope")] <- "slope"
ds.dif.plot$Method <- "RDIF"
ds.dif.plot$Method[ds.dif.plot$method%in% c("lr.chi2", "lr.slope", "lr.intercept")]  <- "LRT"

ds.dif.plot$decision[ds.dif.plot$decision == "fp"] <- "False Positive Rate"
ds.dif.plot$decision[ds.dif.plot$decision == "tp"] <- "True Rositive Rate"
ds.dif.plot$n <- factor(ds.dif.plot$n)
method.order <- unique(ds.dif.plot$Method)[c(1, 2, 3, 4, 5, 6)]
ds.dif.plot$Method <- ordered(ds.dif.plot$Method, method.order)
type.order <- unique(ds.dif.plot$type)
ds.dif.plot$type <- ordered(ds.dif.plot$type, type.order)

# Figure 3 of paper
p1 <- ggplot(ds.dif.plot, aes(y = Value, x = n, fill = Test, pattern = Method)) +
            geom_bar_pattern(
               position = position_dodge(preserve = "single"),
               stat = "Identity",
               color = "black",
               pattern_fill = "White",
               pattern_angle = 45,
               pattern_density = 0.1,
               pattern_spacing = 0.025,
               pattern_key_scale_factor = 0.6) +
            ylab("Value") +
            xlab("Sample size") +
            theme_bw(base_size = 16) +
            #scale_fill_manual(
            #  values = c("grey40", "grey65", "grey85")) +
            #theme(text = element_text(size=18)) +
            scale_fill_manual(
              values = c("#A6CEE3","#B2DF8A","#FB9A99")) +
            scale_pattern_manual(
              values = c(RDIF = "stripe", LRT = "none")) # +
            #scale_y_continuous(
            #  breaks = c(.05, .1, .5, .7, .9),
            #  minor_breaks = c(.025, .075, .6, .8))

p1 +  guides(pattern = guide_legend(override.aes = list(fill = "grey80")),
             fill = guide_legend(override.aes = list(pattern = "none"))) +
             facet_grid(vars(decision), vars(type), scales = "free")


##############################################################################
# ECDI Example
# Data obtained from https://mics.unicef.org/surveys
# Analytical data set: "ECDI_Example_Halpin2022.RData"
##############################################################################

load("ECDI_Example_Halpin2022.RData")

# Rename items
names(data.fiji)[-c(1:3)] <- names(data.vietnam)[-c(1:3)]  <- paste0("ECD", 1:10)

# Normalize sampling weights
data.fiji$chweight <- data.fiji$chweight / sum(data.fiji$chweight) * dim(data.fiji)[1]
data.vietnam$chweight <- data.vietnam$chweight / sum(data.vietnam$chweight) * dim(data.vietnam)[1]

# Fit IRT models
fit.fiji <- mirt(data = data.fiji[, -c(1:3)],
                SE = T,
                survey.weights = data.fiji$chweight)

fit.vietnam <- mirt(data = data.vietnam[, -c(1:3)],
                SE = T,
                survey.weights = data.vietnam$chweight)


irt.mle <- get_irt_pars(list(fit.fiji, fit.vietnam))

# Figure 4 of paper
rho.int <- rho_fun(irt.mle, par = "intercept", alpha = .05, grid.width = .01)
rho.slope <- rho_fun(irt.mle, par = "slope", alpha = .05, grid.width = .01)
gg.data <- data.frame(rho = c(rho.int$rho, rho.slope$rho),
                      theta = c(rho.int$theta, rho.slope$theta),
                      par = c(rep("Intercept", times = length(rho.int$theta)),
                              rep("Slope", times = length(rho.slope$theta))))

p1 <- ggplot(gg.data, aes(x = theta, y = rho)) +
          geom_point(size = 1, color = "#A6CEE3", fill = "#A6CEE3") +
          geom_line(size = 1, linetype = 1, color = "#A6CEE3") +
          theme(text = element_text(size=18))+
          theme_bw(base_size =  18) +
          ylab("Rho") +
          xlab("IRT scale parameter")

p1 + facet_wrap(~ par, ncol = 2, scales = "free")

# Analysis for item intercepts
rdif.theta <- rdif(irt.mle)
rdif.theta.test <- z_test(rdif.theta$est, irt.mle)

# Analysis for item slopes
rdif.sigma <- rdif(irt.mle, par = "slope", starting.value = 1.2)
rdif.sigma.test <- z_test(rdif.sigma$est, irt.mle, par = "slope")

# Analysis for both item parameters
Q <- chi2_test(rdif.theta$est, rdif.sigma$est, irt.mle)

# Table 1 of paper

output <- data.frame(
           intercept.z = rdif.theta.test$z.test,
           intercept.p = rdif.theta.test$p.val,
           slope.z = rdif.sigma.test$z.test,
           slope.p = rdif.sigma.test$p.val,
           both.chi = Q$chi.square,
           both.p = Q$p.val)
xtable(output, digits = 2)

# Compare to LR test.
# Note the mirt::DIF() breaks when using survey.weights
data.full <- rbind(data.fiji, data.vietnam)
invariance = c("slopes", "intercepts", "free_mean", "free_var")

fit.full <- multipleGroup(data = data.full[, -c(1:3)],
                          group = data.full$Country,
                          invariance = invariance,
                          SE = T)
                          # survey.weights = data.full$chweight)

# DIF on intercepts
DIF(MGmodel = fit.full,
    which.par = c("d"),
    scheme = "drop_sequential",
    seq_stat = .05,
    max_run = 2)

# DIF on slopes
DIF(MGmodel = fit.full,
    which.par = c("a1"),
    scheme = "drop_sequential",
    seq_stat = .05,
    max_run = 2)

# DIF on both
DIF(MGmodel = fit.full,
    which.par = c("a1", "d"),
    scheme = "drop_sequential",
    seq_stat = .05,
    max_run = 2)
