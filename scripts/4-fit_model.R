# Source packages
source("fit_fonctions/anneal_custom.R")
source("fit_fonctions/likeli.R")
source("fit_fonctions/analyze_function.R")
source("fit_fonctions/likdisplay.R")

var = list(t0 = "t0", t1 = "t1", ENV = "ENV", lik = "predicted", 
EB = "EB", ET = "ET", EM = "EM", Ha = "Ha", Hv = "Hv")

# state transitions
data = as.data.frame(read.table("../data/data_categorical_RBTM.txt"))
# neighborhood
pred = read.table("../data/data_pred_states_RBTM.txt")

# remove transitions D->C and C->D
test = numeric(nrow(data))
test[data$t0 == "T" & data$t1 == "B"] = 1
test[data$t0 == "B" & data$t1 == "T"] = 1
data = subset(data, test!=1)
pred = subset(pred, test!=1)

data$ENV = data$av_annual_mean_tp
data$EB = pred$B
data$ET = pred$T
data$EM = pred$M 
data$Hv = data$deer
data$Ha = data$moose

# Evaluate initial parameter values
transitions = paste(data$t0,data$t1,sep = "")
sum_transitions = table(transitions)
tot_transitions = table(data$t0)

# estimates for initial parameters
eps_mn = # proba x->R given x is T B or M
(sum_transitions["BR"]+sum_transitions["MR"]+sum_transitions["TR"])/(tot_transitions["M"] + tot_transitions["B"] + tot_transitions["T"])

thetab_mn = # proba M-> B 
sum_transitions["MB"]/tot_transitions["M"]

thetat_mn = # proba M->T
sum_transitions["MT"]/tot_transitions["M"]

betab_mn = # proba T->M
sum_transitions["TM"]*(tot_transitions["B"]+tot_transitions["M"])/tot_transitions["M"]/sum(tot_transitions)

betat_mn = # proba B->M
sum_transitions["BM"]*(tot_transitions["T"]+tot_transitions["M"])/tot_transitions["M"]/sum(tot_transitions)

trRT = sum_transitions["RT"]/tot_transitions["R"]
trRB = sum_transitions["RB"]/tot_transitions["R"]
trRM = sum_transitions["RB"]/tot_transitions["R"]

alphat_mn = trRM/((tot_transitions["M"]+tot_transitions["T"])*(trRB+trRM))
alphab_mn = (trRB+trRM)/(tot_transitions["M"] + tot_transitions["B"])

# transform into logits
logit_eps_mn = log(eps_mn/(1-eps_mn))

logit_thetab_mn = log(thetab_mn/(1-thetab_mn))
logit_thetat_mn = log(thetat_mn/(1-thetat_mn))
 
logit_betab_mn = log(betab_mn/(1-betab_mn))
logit_betat_mn = log(betat_mn/(1-betat_mn))

logit_alphab_mn = log(alphab_mn/(1-alphab_mn))
logit_alphat_mn = log(alphat_mn/(1-alphat_mn))



# List initial parameters
par = list(ab0 =logit_alphab_mn,
ab1 = 0,
ab2 = 0, 
at0 = logit_alphat_mn,
at1 = 0,
at2 = 0,
bb0 = logit_betab_mn,
bb1 = 0,
bb2 = 0,
bt0 = logit_betat_mn,
bt1 = 0,
bt2 = 0,
tt0 = logit_thetat_mn,
tt1 = 0,
tt2 = 0,
tb0 = logit_thetab_mn,
tb1 = 0,
tb2 = 0,
e0 = logit_eps_mn,
e1 = 0,
e2 = 0)

par_lo = list(ab0 = -10,
ab1 = -20,
ab2 = -20, 
at0 = -10,
at1 = -10,
at2 = -20,
bb0 = -10,
bb1 = -10,
bb2 = -25,
bt0 = -30,
bt1 = -10,
bt2 = -10,
tt0 = -10,
tt1 = -10,
tt2 = -25,
tb0 = -10,
tb1 = -10,
tb2 = -20,
e0 = -10,
e1 = -10,
e2 = -10)

par_hi = list(ab0 = 50,
ab1 = 15,
ab2 = 30, 
at0 = 10,
at1 = 10,
at2 = 10,
bb0 = 10,
bb1 = 10,
bb2 = 20,
bt0 = 10,
bt1 = 10,
bt2 = 50,
tt0 = 10,
tt1 = 10,
tt2 = 20,
tb0 = 10,
tb1 = 10,
tb2 = 20,
e0 = 10,
e1 = 10,
e2 = 10)

# Maximum likelihood estimation
test = anneal(model = model, par = par, var = var, source_data = data, 
par_lo = par_lo, par_hi = par_hi, dep_var = "t1", pdf = PDF, 
max_iter = 25000, hessian = FALSE, initial_temp = 0.5)
par = test$best_pars
write.table(test$best_pars,"../estimated_par/par.txt")




