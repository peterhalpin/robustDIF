# Numerical example with ECDI2030
library(mirt)
library(ggplot2)
library(robustDIF)

# Load and format data
mex.file <- "~/Dropbox/FILES PSYCHOMETRICS PAPER/Psychometric article/Data/mexico_data.csv"

pal.file <- "~/Dropbox/FILES PSYCHOMETRICS PAPER/Psychometric article/Data/palestine_data.csv"

mexico.data <- read.csv(mex.file)
palestine.data <- read.csv(pal.file)
mexico.data$Country <- "Mexico"
palestine.data$Country <- "Palestine"

# Palestine weights do not sum to sample size
palestine.data$chweight <- palestine.data$chweight / sum(palestine.data$chweight) * nrow(palestine.data)

ECDI.names <- c("new_ECD1", "new_ECD2",  "new_ECD4", "new_ECD7","new_ECD11", "new_ECD12", "new_ECD13", "new_ECD14", "new_ECD18", "new_ECD20", "new_ECD21", "new_ECD28", "R_ECD35", "new_ECD41", "new_ECD43", "new_ECD44", "R_ECD50", "new_ECD53", "new_ECD54", "new_ECD57")

#ECDI.names <- c("new_ECD1", "new_ECD2",  "new_ECD4", "new_ECD7","new_ECD11", "new_ECD12", "new_ECD13", "new_ECD14", "new_ECD16", "new_ECD17", "new_ECD18", "new_ECD20", "new_ECD21", "new_ECD22", "new_ECD23", "new_ECD24", "new_ECD26", "new_ECD27", "new_ECD28", "new_ECD31", "R_ECD35", "R_ECD37", "new_ECD41", "new_ECD43", "new_ECD44", "R_ECD45", "R_ECD46", "new_ECD48", "R_ECD50", "new_ECD53", "new_ECD54", "new_ECD57")

var.names <- c("Country", ECDI.names)

data <- rbind(mexico.data[names(mexico.data)%in%var.names],
              palestine.data[names(palestine.data)%in%var.names])
learning.items <- c(1:4,6:11)
learning.data <- data[, ECDI.names[learning.items]]
fit.mex <- mirt(data = learning.data[data$Country == "Mexico", ],
                SE = T,
                survey.weights = data$chweight[data$Country == "Mexico"])

fit.pal <- mirt(data = learning.data[data$Country == "Palestine", ],
                SE = T,
                survey.weights = data$chweight[data$Country == "Palestine"])

irt.mle <- get_irt_pars(list(fit.mex, fit.pal))

# Check distribution of scaling functions for symmetry
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
get_starts(irt.mle)

### Analysis for item intercepts
rdif.theta <- rdif(irt.mle)
rdif.theta.test <- z_test(rdif.theta$est, irt.mle)

# Check direction of biased
(y_fun(irt.mle) - rdif.theta$est)[rdif.theta.test$p.val < .05]

# Items without DIF
ECDI.names[learning.items][rdif.theta.test$p.val > .05]
# Result: [1] "new_ECD1"  "new_ECD12" "new_ECD13" "new_ECD18" "new_ECD20"

### Analysis for item slopes
rdif.sigma <- irls(irt.mle, par = "slope")
rdif.sigma.test <- z_test(rdif.sigma$est, irt.mle, par = "slope")

(y_fun(irt.mle, par = "slope") - rdif.sigma$est)[rdif.sigma.test$p.val < .05]

ECDI.names[learning.items][rdif.sigma.test$p.val > .05]
# Result: [1] "new_ECD1"  "new_ECD4"  "new_ECD18" "new_ECD20"

### Analysis for both item parameters
Q <- chi2_test(rdif.theta$est, rdif.sigma$est, irt.mle)
ECDI.names[learning.items][Q$p.val > .05]

#[1] "new_ECD1"  "new_ECD18" "new_ECD20"



### LR tests: All tests omit all items
mg.mod.int <- multipleGroup(
                data = learning.data,
                model = 1,
                group = data$Country,
                invariance = c("intercepts", "free_means"),
                survey.weights = data$chweight)

lr.int <- DIF(mg.mod.int, which.par = "d", scheme = "drop", Wald = F)

mg.mod.slope <- multipleGroup(
                data = learning.data,
                model = 1,
                group = data$Country,
                invariance = c("slopes", "free_vars"),
                survey.weights = data$chweight)

lr.slope <- DIF(mg.mod.slope, which.par = "a1", scheme = "drop", Wald = F)

mg.mod <- multipleGroup(
                data = learning.data,
                model = 1,
                group = data$Country,
                invariance = c("intercepts", "free_means", "slopes", "free_vars"),
                survey.weights = data$chweight)

lr.chi <- DIF(mg.mod, which.par = c("a1", "d"), scheme = "drop", Wald = F)


