# Sim study 2 July 6, 2022
library(ggplot2)
library(dplyr)
#devtools::install_github("coolbutuseless/ggpattern")
library(ggpattern)

### Data gen

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


sim.study2.path <- "~/Dropbox/Academic/Manuscripts/DIF_via_scaling/data_analyses/sim2.july5.2022.RData"
 #save(sim.study2, file = sim.study2.path)
load(file = sim.study2.path)


### Data processing

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
head(ds.dif.plot)

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
            #theme(text = element_text(size=18)) +
            theme_bw(base_size = 18) +
            scale_fill_manual(
             values = c("grey40", "grey65", "grey85")) +
            #values = c("#A6CEE3","#B2DF8A","#FB9A99")) +
            scale_pattern_manual(
              values = c(R.DIF = "stripe", LRT = "none")) +
            scale_y_continuous(
              breaks = c(.05, .1, .5, .7, .9),
              minor_breaks = c(.025, .075, .6, .8))

p1 +  guides(pattern = guide_legend(override.aes = list(fill = "grey80")),
             fill = guide_legend(override.aes = list(pattern = "none"))) +
      facet_grid(vars(decision), vars(type), scales = "fixed")
