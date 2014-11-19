rm(list=ls())

source("4-init_params.R")

var = list(st0 = "st0", st1 = "st1", ENV1 = "ENV1", ENV2 = "ENV2", lik = "predicted", 
EB = "EB", ET = "ET", EM = "EM", Ha = "Ha", Hv = "Hv")


# Maximum likelihood estimation
source("fit_fonctions/analyze_function.R")
source("fit_fonctions/likdisplay_ibou.R")
source("fit_fonctions/support_limits_ibou.R")
source("fit_fonctions/likeli.R")


source("3-transition_model.R")
source("fit_fonctions/anneal_ibou.R")


coarse= anneal(model = model, par = par, var = var, source_data = data, 
dep_var = "st1", pdf = PDF, par_initStep = par_initStep, par_testFct = testBounds, max_iter = 10000, initial_temp = 3, note = "", progress = TRUE, display=FALSE, support = FALSE, min_change = 1, min_drops = 10)

#
save(coarse, file="../estimated_params/coarse_m3")

source("fit_fonctions/anneal_simple.R")

par = coarse$best_pars
par_lo = lapply(par, function(x){x-5*abs(x)})
par_hi = lapply(par, function(x){x+5*abs(x)})
coarse= anneal(model = model, par = par, var = var, source_data = data, 
dep_var = "st1", pdf = PDF, par_lo = par_lo, par_hi=par_hi, max_iter = 10000, initial_temp = 3, min_change = 1, min_drops = 10)

save(coarse, file="../estimated_params/coarse_m3_step2") 
save(par_lo, par_hi, file="../estimated_params/coarse_m3_step2_lim")

par = coarse$best_pars
par_lo = lapply(par, function(x){x-5})
par_hi = lapply(par, function(x){x+5})
coarse= anneal(model = model, par = par, var = var, source_data = data, 
dep_var = "st1", pdf = PDF, par_lo = par_lo, par_hi=par_hi, max_iter = 10000, initial_temp = 3, min_change = 1, min_drops = 10)

save(coarse, file="../estimated_params/coarse_m3_step3") 
save(par_lo, par_hi, file="../estimated_params/coarse_m3_step3_lim")

par = coarse$best_pars
par_lo = lapply(par, function(x){x-5})
par_hi = lapply(par, function(x){x+5})
coarse= anneal(model = model, par = par, var = var, source_data = data, 
dep_var = "st1", pdf = PDF, par_lo = par_lo, par_hi=par_hi, max_iter = 10000, initial_temp = 3, min_change = 1, min_drops = 10)

save(coarse, file="../estimated_params/coarse_m3_step4")
save(par_lo, par_hi, file="../estimated_params/coarse_m3_step4_lim")


#write.table(coarse$best_pars,"../estimated_params/par_herbivores_m3.txt")
#
#
#
#
