################################################################################
# Numerical examples for: Halpin (2022) Differential Item Functioning
# Via Robust Scaling.
################################################################################

### Required packages
library(robustDIF)
library(mirt)
library(GPCMlasso)
library(difR)
library(dplyr)
library(ggplot2)
#devtools::install_github("coolbutuseless/ggpattern")
library(ggpattern)

# Helper functions for running simulations
source("Halpin2022_R/data_sim_functions.R")

###################################################
# Sim Study 1
# First run: June 28-29, 2022
# Second run: July 5-6, 2022 (new bsq_weights)
###################################################

# Data gen
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
sim.study1.path <- "Halpin2022_R/sim1.july6.2022.RData"
#save(sim.study1, file = sim.study1.path)

# Data processing
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

ds.dif.long <- ds.dif.long[ds.dif.long$method != "rdif.chi2", ] # Not necessary
ds.dif.long$method <- ordered(ds.dif.long$method, unique(ds.dif.long$method))
ds.dif.long$decision[ds.dif.long$decision == "fp"] <- "False Positive Rate"
ds.dif.long$decision[ds.dif.long$decision == "tp"] <- "True Positive Rate"
names(ds.dif.long)[1] <- "Method"

# Plot of decision errors (bw and color)
p1 <- ggplot(ds.dif.long, aes(y = Value, x = n.biased, group = Method)) +
            geom_point(aes(color = Method, shape = Method, fill = Method),
                       size = 2.5) +
            geom_line(aes(color = Method), size = 1, linetype = 1) +
            geom_hline(yintercept = .875, col = 'grey25', linetype = 2) +
            geom_hline(yintercept = .05, col = 'grey65', linetype = 2) +
            ylab("Value") +
            xlab("Number of biased items (out of 15)") +
            scale_shape_manual(values=c(21:25)) +
            # scale_colour_brewer(palette = "Paired") +
            # scale_fill_brewer(palette = "Paired") +
            # theme(text = element_text(size=18))
            scale_color_grey() +
            scale_fill_grey() +
            theme_bw(base_size =  18)

p1 + facet_wrap(~ decision, nrow = 2)

# Plot of scale distribution (bw and color)
temp.scale <- lapply(sim.study1, function(x) {lapply(x, function(y) y$scale)})
ds.scale <-  Reduce(rbind, lapply(temp.scale, function(x) Reduce(rbind, x)))
ds.scale$n.biased <- as.factor(rep(0:8, each = n.reps))

facet.labels <- paste0("N.DIF: ", 0:8,  " out of 15")
facet_labeller <- function(variable, value){
  return(facet.labels[value])
}

p2 <- ggplot(ds.scale, aes(x = theta)) +
            #geom_histogram(col = "white", fill = "#A6CEE3") +
            #theme(text = element_text(size=18)) +
            geom_histogram(col = "white", fill = 'grey65') +
            theme_bw(base_size =  18)
            ylab("Count") +
            xlab("R-DIF estimate of the IRT scale parameter")

p2 + facet_wrap(~ n.biased, nrow = 3, labeller = facet_labeller)

###################################################
# Sim study 2
# First run July 6, 2022
###################################################

# Data gen
n.reps = 500
n.items =  10

# sample size = 200
n.persons = 200
ds.a200  <- sim_study2(n.reps, n.persons, n.items, bias = c(.5, 0))
ds.b200  <- sim_study2(n.reps, n.persons, n.items, bias = c(0, 1))
ds.c200.r  <- sim_study2(250, n.persons, n.items, bias = c(.35, .5))

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


sim.study2.path <- "Halpin2022_R/sim2.july5.2022.RData"
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
ds.dif.plot$Test[ds.dif.plot$method%in%c("lr.chi2", "rdif.chi2")] <- "chi.square"
ds.dif.plot$Test[ds.dif.plot$method%in%c("lr.slope", "rdif.slope")] <- "slope"
ds.dif.plot$Method <- "R.DIF"
ds.dif.plot$Method[ds.dif.plot$method%in% c("lr.chi2", "lr.slope", "lr.intercept")]  <- "LRT"

ds.dif.plot$decision[ds.dif.plot$decision == "fp"] <- "False Positive Rate"
ds.dif.plot$decision[ds.dif.plot$decision == "tp"] <- "True Rositive Rate"
ds.dif.plot$n <- factor(ds.dif.plot$n)
method.order <- unique(ds.dif.plot$Method)[c(1, 2, 3, 4, 5, 6)]
ds.dif.plot$Method <- ordered(ds.dif.plot$Method, method.order)
type.order <- unique(ds.dif.plot$type)
ds.dif.plot$type <- ordered(ds.dif.plot$type, type.order)

# Plot of decision errors (bw and color). Uses ggpattern. Takes a while to run
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
            theme_bw(base_size = 18) +
            scale_fill_manual(
              values = c("grey40", "grey65", "grey85")) +
            #theme(text = element_text(size=18)) +
            #scale_fill_manual(
            #values = c("#A6CEE3","#B2DF8A","#FB9A99")) +
            scale_pattern_manual(
              values = c(R.DIF = "stripe", LRT = "none")) +
            scale_y_continuous(
              breaks = c(.05, .1, .5, .7, .9),
              minor_breaks = c(.025, .075, .6, .8))

p1 +  guides(pattern = guide_legend(override.aes = list(fill = "grey80")),
             fill = guide_legend(override.aes = list(pattern = "none"))) +
      facet_grid(vars(decision), vars(type), scales = "fixed")

###################################################
# ECDI Example
# Data not available for dissemination, but code presented below
###################################################

# Load and format data
mex.file <- "mexico_data.csv"

pal.file <- "palestine_data.csv"

mexico.data <- read.csv(mex.file)
palestine.data <- read.csv(pal.file)
mexico.data$Country <- "Mexico"
palestine.data$Country <- "Palestine"

ECDI.names <- c("new_ECD1", "new_ECD2",  "new_ECD4", "new_ECD7","new_ECD11", "new_ECD12", "new_ECD13", "new_ECD14", "new_ECD18", "new_ECD20", "new_ECD21", "new_ECD28", "R_ECD35", "new_ECD41", "new_ECD43", "new_ECD44", "R_ECD50", "new_ECD53", "new_ECD54", "new_ECD57")

var.names <- c("Country", ECDI.names)

data <- rbind(mexico.data[names(mexico.data)%in%var.names],
              palestine.data[names(palestine.data)%in%var.names])

learning.items <- c(1:4,6:11)
learning.data <- data[, ECDI.names[learning.items]]

# Fit irt models
fit.mex <- mirt(data = learning.data[data$Country == "Mexico", ],
                SE = T,
                survey.weights = data$chweight[data$Country == "Mexico"])

fit.pal <- mirt(data = learning.data[data$Country == "Palestine", ],
                SE = T,
                survey.weights = data$chweight[data$Country == "Palestine"])

irt.mle <- get_irt_pars(list(fit.mex, fit.pal))

# Check distribution of scaling functions for assymmetry
hist(y_fun(irt.mle, par = "intercept"))
hist(y_fun(irt.mle, par = "slope", log = F), breaks = 10)

# Check rho function for local solutions
rho.int <- rho_fun(irt.mle, par = "intercept", alpha = .05, grid.width = .01)
rho.slope <- rho_fun(irt.mle, par = "slope", alpha = .05, grid.width = .01)
gg.data <- data.frame(rho = c(rho.int$rho, rho.slope$rho)/10,
                      theta = c(rho.int$theta, rho.slope$theta),
                      par = c(rep("Intercept", times = length(rho.int$theta)),
                              rep("Slope", times = length(rho.slope$theta))))

p1 <- ggplot(gg.data, aes(x = theta, y = rho)) +
          geom_point(size = 1.2, color = "#A6CEE3", fill = "#A6CEE3") +
          geom_line(size = 1, linetype = 1, color = "#A6CEE3") +
          theme(text = element_text(size=18))+
          # geom_point(size = 1.2, color = "grey60", fill = "grey60") +
          # geom_line(size = 1, linetype = 1, color = "grey60") +
          # theme_bw(base_size =  18) +
          ylab("Rho") +
          xlab("IRT scale parameter") +
          scale_x_continuous(
              breaks = c(-.25, 0, .25, .50, 1, 1.25, 1.50, 1.75, 2))

p1 + facet_wrap(~ par, ncol = 2, scales = "free_x")

# Analysis for item intercepts
rdif.theta <- rdif(irt.mle)
rdif.theta.test <- z_test(rdif.theta$est, irt.mle)

# Check direction of biased
(y_fun(irt.mle) - rdif.theta$est)[rdif.theta.test$p.val < .05]

# Items without DIF
ECDI.names[learning.items][rdif.theta.test$p.val > .05]
# Result: [1] "new_ECD1"  "new_ECD12" "new_ECD13" "new_ECD18" "new_ECD20"

# Analysis for item slopes
rdif.sigma <- irls(irt.mle, par = "slope")
rdif.sigma.test <- z_test(rdif.sigma$est, irt.mle, par = "slope")

# Check direction of biased
(y_fun(irt.mle, par = "slope") - rdif.sigma$est)[rdif.sigma.test$p.val < .05]

# Items without DIF
ECDI.names[learning.items][rdif.sigma.test$p.val > .05]
# Result: [1] "new_ECD1"  "new_ECD4"  "new_ECD18" "new_ECD20"

# Analysis for both item parameters
Q <- chi2_test(rdif.theta$est, rdif.sigma$est, irt.mle)
ECDI.names[learning.items][Q$p.val > .05]
#Result: [1] "new_ECD1"  "new_ECD18" "new_ECD20"


