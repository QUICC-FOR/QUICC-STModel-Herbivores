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
save(coarse, file="../estimated_params/coarse_m3_epsCor")
write.table(coarse$best_pars,"../estimated_params/par_herbivores_m3_epsCor.txt")



data$deer = rep(0, nrow(data))
data$moose = rep(0, nrow(data))
coarse_veget= anneal(model = model, par = par, var = var, source_data = data, 
dep_var = "st1", pdf = PDF, par_initStep = par_initStep, par_testFct = testBounds, max_iter = 10000, initial_temp = 3, note = "", progress = TRUE, display=FALSE, support = FALSE, min_change = 1, min_drops = 10)

#
#par = coarse_veget$best_pars
save(coarse_veget, file="../estimated_params/coarse_veget_m3_epsCor")
write.table(coarse_veget$best_pars,"../estimated_params/par_veget_m3_epsCor.txt")


#load("../estimated_params/coarse_rf3")
#par = coarse$best_pars
#
#fine= anneal(model = model, par = par, var = var, source_data = data, 
#dep_var = "st1", pdf = PDF, par_initStep = par_initStep, par_testFct = testBounds, max_iter = 10000, initial_temp = 3, note = "", progress = TRUE, display=FALSE, support = FALSE, min_change = 0.01, min_drops = 50)
#
##
#save(fine, file="../estimated_params/fine_rf3")
#write.table(coarse$best_pars,"../estimated_params/par_herbivores_rf3.txt")
#
#
#load("../estimated_params/coarse_veget_rf3")
#par = coarse$best_pars
#
#data$deer = rep(0, nrow(data))
#data$moose = rep(0, nrow(data))
#fine_veget= anneal(model = model, par = par, var = var, source_data = data, 
#dep_var = "st1", pdf = PDF, par_initStep = par_initStep, par_testFct = testBounds, max_iter = 10000, initial_temp = 3, note = "", progress = TRUE, display=FALSE, support = FALSE, min_change = 0.01, min_drops = 50)
#
##
##par = coarse_veget$best_pars
#save(fine_veget, file="../estimated_params/fine_veget_rf3")
#write.table(fine_veget$best_pars,"../estimated_params/par_veget_rf3.txt")
#

---------------------------------------------------------------

   # Calculate support limits
   load("../estimated_params/coarse_m3")
   load("../estimated_params/coarse_veget_m3")

    best_par = coarse$best_par
    for (i in 1:length(parnames)) {
      par[[parnames[i]]] <- best_par[[parnames[i]]]
    }
    limits<-support_limits(model = model, par = par, var = var, source_data = source_data, pdf = pdf, delta = delta, slimit = slimit)
    upper_limit<-limits$upper_limits
    lower_limit<-limits$lower_limits
    
    cbind(unlist(best_par), unlist(upper_limit), unlist(lowe_limit))





