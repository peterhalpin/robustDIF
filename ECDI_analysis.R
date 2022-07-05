# Numerical example with ECDI

library(mirt)
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

var.names <- c("HL4", "windex5", "ece", "ind65", "ind62", "funcdiff", "chweight", "Country", "CAGE", ECDI.pool.names)

data <- rbind(mexico.data[names(mexico.data)%in%var.names],
              palestine.data[names(palestine.data)%in%var.names])

ECDI.names[1:11]
learning.data <- data[, ECDI.names[c(1:4,6:11)]]

fit.mex <- mirt(data = learning.data[data$Country == "Mexico", ],
                model = 1,
                itemtype = "2PL",
                SE = T,
                survey.weights = data$chweight[data$Country == "Mexico"])

fit.pal <- mirt(data = learning.data[data$Country == "Palestine", ],
                model = 1,
                itemtype = "2PL",
                SE = T,
                survey.weights = data$chweight[data$Country == "Palestine"])

irt.mle <- get_irt_mle(list(fit.mex, fit.pal))

hist(y_fun(irt.mle, par = "intercept"))
hist(y_fun(irt.mle, par = "slope", log = F), breaks = 10)

rho.int <- rho_fun(irt.mle, par = "intercept", alpha = .05)
rho.slope <- rho_fun(irt.mle, par = "slope", log = F, alpha = .05)
par(mfrow = c(1, 2))
plot(rho.int$theta, rho.int$rho, type = "l")
plot(rho.slope$theta, rho.slope$rho, type = "l")

get_starts(irt.mle, par = "intercept")
get_starts(irt.mle, par = "slope")

rdif.theta <- irls(irt.mle, alpha = .05)
y_fun(irt.mle) - rdif.theta$est
ECDI.names[c(1:4,6:11)][z_test(rdif.theta$est, irt.mle)$p.val > .05]

#[1] "new_ECD2"  "new_ECD4"  "new_ECD13" "new_ECD18" "new_ECD20"

rdif.sigma <- irls(irt.mle, par = "slope")
y_fun(irt.mle, par = "slope") - rdif.sigma$est
ECDI.names[c(1:4,6:11)][z_test(rdif.sigma$est, irt.mle, par = "slope")$p.val > .05]
#[1] "new_ECD1"  "new_ECD2"  "new_ECD4"  "new_ECD13" "new_ECD20"

q <- chi2_test(rdif.theta$est, rdif.sigma$est, irt.mle)
ECDI.names[c(1:4,6:11)][q$p.val > .05]

#[1] "new_ECD2"  "new_ECD4"  "new_ECD13" "new_ECD20"

# LR test omits all items for slopes, int, and chi
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


