rm(list=ls())

source("4-init_params.R")

var = list(t0 = "t0", t1 = "t1", ENV1 = "ENV1", ENV2 = "ENV2", lik = "predicted", 
EB = "EB", ET = "ET", EM = "EM", Ha = "Ha", Hv = "Hv")


# Maximum likelihood estimation
source("fit_fonctions/analyze_function.R")
source("fit_fonctions/likdisplay_ibou.R")
source("fit_fonctions/support_limits_ibou.R")
source("fit_fonctions/likeli.R")


source("3-transition_model.R")
source("fit_fonctions/anneal_ibou.R")

test1= anneal(model = model, par = par, var = var, source_data = data, 
dep_var = "t1", pdf = PDF, par_initStep = par_initStep, par_testFct = testBounds, max_iter = 10000, initial_temp = 0.5, note = "", progress = TRUE, display=FALSE, support = FALSE)

#
par = test1$best_pars
save(test1, file="../estimated_par/test1")
write.table(test1$best_pars,"../estimated_par/par_herbivores.txt")
#
#


