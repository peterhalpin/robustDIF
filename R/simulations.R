# Sim study 1 June 27, 2022

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
save(sim.study1, file = "/Users/halpin/Dropbox/Academic/Manuscripts/DIF_via_scaling/data_analyses/sim1.june27.2022.RData")

