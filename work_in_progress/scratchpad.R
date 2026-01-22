

library(mirt)
library(lavaan)
#source("work_in_progress/read_functions.R")
#source("R/robustDIF_functions.R")
#source("R/DTF_functions.R")
source("work_in_progress/sim_functions.R")


# Data gen
n.persons = 500
bias = .5
impact = c(.5,1)
n.items = 5
n.biased = 1



sim <- sim_dat(n.persons = n.persons,
          n.items = n.items,
          n.biased = n.biased,
          bias = bias,
          impact = impact)

mg.mirt.object<- multipleGroup(rbind(sim$dat0, sim$dat1),
                            group = factor(rep(c(0, 1), each = n.persons)),
                            itemtype = "graded",
                            SE = T)

mle <- get_model_parms(mg.mirt.object)
delta_test(mle)
rdif_chisq_test((mle))
rdif_z_test((mle))
rdif_z_test(mle, par = "slope")

rdif.eg <- mle
save(rdif.eg, file = "data/rdif.eg.rda")
delta_test(rdif.eg)


theta <- NULL

object <- mirt.object.list <- list(fit0 = mirt(sim$dat0, SE = T),
                         fit1 = mirt(sim$dat1, SE = T))


mle <- get_mle(mirt.object.list)
rdif(mle)

# checking against old code
irt.mle <- list(par0 = mle$est$group.1,
                par1 = mle$est$group.2,
                vcov0 = mle$var.cov$group.1,
                vcov1 = mle$var.cov$group.2)
names(irt.mle$par0) <- names(irt.mle$par1) <- c("a", "d")
row.names(irt.mle$vcov0) <- row.names(irt.mle$vcov1) <-
colnames(irt.mle$vcov0) <- colnames(irt.mle$vcov1) <-
  paste0(rep(c("a", "d"), times = n.items), rep(1:n.items, each = 2))

rdif_old(irt.mle)

y_fun(mle)
y_fun(mle, par = "slope")

grad_intercept(mle)
grad_intercept_old(.5, irt.mle)

grad_slope(mle, 1)
grad_slope_old(.5, irt.mle)

joint_vcov(mle)

var_y(mle, theta)
var_y_old(theta, irt.mle)

var_y(theta, mle, par = "slope", log = T)
var_y_old(theta, irt.mle, par = "slope", log = T)


cov_yz(theta, 1, mle)
cov_yz_old(theta, 1, irt.mle)

get_starts(mle)
get_starts_old(irt.mle)

rho_fun(mle)
rho_fun_old(irt.mle)

rdif(mle)
rdif_old(irt.mle)
z_test(.5, mle)
z_test_old(.5, irt.mle)
z_test_old(1, irt.mle, par = "slope")
2.3^2 + .8^2
chisq_test_old(.5, 1, irt.mle)
chisq_test(.5, 1, mle)

joint_vcov(mle)
joint_vcov_old(irt.mle)

## Checking other reads (not mplus)
object <- mirt.object.multigroup <- multipleGroup(Reduce(rbind, sim),
                               group = rep(c("0", "1"), each = n.persons),
                               SE = T)

get_pars(mirt.object.multigroup)

model <- paste0("f1 =~ ", paste0("Item_", 1:n.items, collapse = "+" ))

object <- lavaan.object.list <- list(fit0 = cfa(model = model,
                                     data = data.frame(sim$dat0),
                                     ordered = names(sim$dat0),
                                     meanstructure = T,
                                     std.lv = T),
                          fit1 = cfa(model = model,
                                     data = data.frame(sim$dat1),
                                     ordered = names(sim$dat1),
                                     meanstructure = T,
                                     std.lv = T))

get_pars(lavaan.object.list)

dat <- Reduce(rbind, sim)
dat$group <- rep(c("0", "1"), each = n.persons)
object <- lavaan.object.multigroup <- cfa(model = model,
                                     data = dat,
                                     ordered = names(dat),
                                     group = "group",
                                     meanstructure = T,
                                     std.lv = T)
get_pars(lavaan.object.list)

