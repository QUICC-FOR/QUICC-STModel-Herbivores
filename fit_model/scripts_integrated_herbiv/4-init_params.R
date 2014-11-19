# state transitions
data = as.data.frame(read.table("../data/data_reshaped_RBTM_herbivores.txt"))
# neighborhood

# remove transitions D->C and C->D
test = numeric(nrow(data))
test[data$st0 == "T" & data$st1 == "B"] = 1
test[data$st0 == "B" & data$st1 == "T"] = 1
data = subset(data, test!=1)

data$ENV1 = data$annual_mean_temp
data$ENV2 = data$annual_pp/1000

data$Hv = data$deer*75/100 ## kg/hectare
data$Ha = data$moose*350/1000 ## kg/hectare
#data$Hv = rep(0, nrow(data))
#data$Ha = rep(0, nrow(data))

head(data)

# Evaluate initial parameter values
transitions = paste(data$st0,data$st1,sep = "")
sum_transitions = table(transitions)
tot_transitions = table(data$st0)

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

sumStates = table(c(data$st0,data$st1))/2
names(sumStates) = c("B", "M", "R", "T")
ET_mn = sumStates["T"]/nrow(data)
EB_mn = sumStates["B"]/nrow(data)
EM_mn = sumStates["M"]/nrow(data)
ER_mn = sumStates["R"]/nrow(data)

# transform into logits
logit_eps_mn = log(eps_mn/(1-eps_mn))

logit_thetab_mn = log(thetab_mn/(1-thetab_mn))
logit_thetat_mn = log(thetat_mn/(1-thetat_mn))
 
logit_betab_mn = log(betab_mn/(1-betab_mn))
logit_betat_mn = log(betat_mn/(1-betat_mn))

logit_alphab_mn = log(alphab_mn/(1-alphab_mn))
logit_alphat_mn = log(alphat_mn/(1-alphat_mn))

logit_ET_mn = log(ET_mn/(1-ET_mn))
logit_EB_mn = log(EB_mn/(1-EB_mn))
logit_EM_mn = log(EM_mn/(1-EM_mn))
logit_ER_mn = log(ER_mn/(1-ER_mn))

# List initial parameters
pars = list(ab0 =logit_alphab_mn, ab1 = 0, ab2 = 0, ab3=0, ab4=0, ab5=0, ab6=0,
at0 = logit_alphat_mn, at1 = 0, at2 = 0, at3=0, at4=0, at5=0, at6=0,
bb0 = logit_betab_mn, bb1 = 0, bb2 = 0, bb3=0, bb4=0, bb5=0, bb6=0,
bt0 = logit_betat_mn, bt1 = 0, bt2 = 0, bt3=0, bt4=0, bt5=0, bt6=0,
tt0 = logit_thetat_mn, tt1 = 0, tt2 = 0, tt3 =0, tt4=0, tt5=0, tt6=0,
tb0 = logit_thetab_mn, tb1 = 0, tb2 = 0, tb3=0, tb4=0, tb5=0, tb6=0,
e0 = logit_eps_mn, e1 = 0, e2 = 0, e3=0, e4=0, e5=0, e6=0
)

# coeff variation
cvar = 2

par_lo = list(ab0 = logit_alphab_mn - abs(logit_alphab_mn*cvar),
ab1 = logit_alphab_mn/max(abs(data$ENV1)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV1))),
ab2 = logit_alphab_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV1^2))), 
ab3 = logit_alphab_mn/max(abs(data$ENV2)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV2))), 
ab4= logit_alphab_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV2^2))),
ab5= logit_alphab_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV1^3))),
ab6= logit_alphab_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV2^3))),
at0 = logit_alphat_mn - abs(logit_alphat_mn*cvar),
at1 = logit_alphat_mn/max(abs(data$ENV1)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV1))),
at2 = logit_alphat_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV1^2))),
at3 = logit_alphat_mn/max(abs(data$ENV2)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV2))),
at4 = logit_alphat_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV2^2))),
at5= logit_alphat_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV1^3))),
at6= logit_alphat_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV2^3))),
bb0 = logit_betab_mn - abs(logit_betab_mn*cvar),
bb1 = logit_betab_mn/max(abs(data$ENV1)) - abs(cvar*logit_betab_mn/max(abs(data$ENV1))),
bb2 = logit_betab_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_betab_mn/max(abs(data$ENV1^2))),
bb3 = logit_betab_mn/max(abs(data$ENV2)) - abs(cvar*logit_betab_mn/max(abs(data$ENV2))),
bb4 = logit_betab_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_betab_mn/max(abs(data$ENV2^2))),
bb5= logit_betab_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_betab_mn/max(abs(data$ENV1^3))),
bb6= logit_betab_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_betab_mn/max(abs(data$ENV2^3))),
bt0 = logit_betat_mn - abs(logit_betat_mn*cvar),
bt1 = logit_betat_mn/max(abs(data$ENV1)) - abs(cvar*logit_betat_mn/max(abs(data$ENV1))),
bt2 = logit_betat_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_betat_mn/max(abs(data$ENV1^2))),
bt3 =logit_betat_mn/max(abs(data$ENV2)) - abs(cvar*logit_betat_mn/max(abs(data$ENV2))),
bt4 = logit_betat_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_betat_mn/max(abs(data$ENV2^2))),
bt5= logit_betat_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_betat_mn/max(abs(data$ENV1^3))),
bt6= logit_betat_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_betat_mn/max(abs(data$ENV2^3))),
tt0 = logit_thetat_mn - abs(logit_thetat_mn*cvar),
tt1 = logit_thetat_mn/max(abs(data$ENV1)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV1))),
tt2 = logit_thetat_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV1^2))),
tt3 = logit_thetat_mn/max(abs(data$ENV2)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV2))),
tt4 = logit_thetat_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV2^2))),
tt5= logit_thetat_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV1^3))),
tt6= logit_thetat_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV2^3))),
tb0 = logit_thetab_mn - abs(logit_thetab_mn*cvar),
tb1 = logit_thetab_mn/max(abs(data$ENV1)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV1))),
tb2 = logit_thetab_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV1^2))),
tb3 = logit_thetab_mn/max(abs(data$ENV2)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV2))),
tb4 = logit_thetab_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV2^2))),
tb5= logit_thetab_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV1^3))),
tb6= logit_thetab_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV2^3))),
e0 = logit_eps_mn - abs(logit_eps_mn*cvar),
e1 = logit_eps_mn/max(abs(data$ENV1)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1))),
e2 = logit_eps_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1^2))),
e3 = logit_eps_mn/max(abs(data$ENV2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2))), 
e4 = logit_eps_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2^2))), 
e5= logit_eps_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1^3))),
e6= logit_eps_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2^3)))
)

par_hi = list(ab0 = logit_alphab_mn + abs(logit_alphab_mn*cvar),
ab1 = logit_alphab_mn/max(abs(data$ENV1)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV1))),
ab2 = logit_alphab_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV1^2))), 
ab3 = logit_alphab_mn/max(abs(data$ENV2)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV2))), 
ab4= logit_alphab_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV2^2))),
ab5= logit_alphab_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV1^3))),
ab6= logit_alphab_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV2^3))),
at0 = logit_alphat_mn + abs(logit_alphat_mn*cvar),
at1 = logit_alphat_mn/max(abs(data$ENV1)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV1))),
at2 = logit_alphat_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV1^2))),
at3 = logit_alphat_mn/max(abs(data$ENV2)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV2))),
at4 = logit_alphat_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV2^2))),
at5= logit_alphat_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV1^3))),
at6= logit_alphat_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV2^3))),
bb0 = logit_betab_mn + abs(logit_betab_mn*cvar),
bb1 = logit_betab_mn/max(abs(data$ENV1)) + abs(cvar*logit_betab_mn/max(abs(data$ENV1))),
bb2 = logit_betab_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_betab_mn/max(abs(data$ENV1^2))),
bb3 = logit_betab_mn/max(abs(data$ENV2)) + abs(cvar*logit_betab_mn/max(abs(data$ENV2))),
bb4 = logit_betab_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_betab_mn/max(abs(data$ENV2^2))),
bb5= logit_betab_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_betab_mn/max(abs(data$ENV1^3))),
bb6= logit_betab_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_betab_mn/max(abs(data$ENV2^3))),
bt0 = logit_betat_mn + abs(logit_betat_mn*cvar),
bt1 = logit_betat_mn/max(abs(data$ENV1)) + abs(cvar*logit_betat_mn/max(abs(data$ENV1))),
bt2 = logit_betat_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_betat_mn/max(abs(data$ENV1^2))),
bt3 =logit_betat_mn/max(abs(data$ENV2)) + abs(cvar*logit_betat_mn/max(abs(data$ENV2))),
bt4 = logit_betat_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_betat_mn/max(abs(data$ENV2^2))),
bt5= logit_betat_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_betat_mn/max(abs(data$ENV1^3))),
bt6= logit_betat_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_betat_mn/max(abs(data$ENV2^3))),
tt0 = logit_thetat_mn + abs(logit_thetat_mn*cvar),
tt1 = logit_thetat_mn/max(abs(data$ENV1)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV1))),
tt2 = logit_thetat_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV1^2))),
tt3 = logit_thetat_mn/max(abs(data$ENV2)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV2))),
tt4 = logit_thetat_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV2^2))),
tt5= logit_thetat_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV1^3))),
tt6= logit_thetat_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV2^3))),
tb0 = logit_thetab_mn + abs(logit_thetab_mn*cvar),
tb1 = logit_thetab_mn/max(abs(data$ENV1)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV1))),
tb2 = logit_thetab_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV1^2))),
tb3 = logit_thetab_mn/max(abs(data$ENV2)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV2))),
tb4 = logit_thetab_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV2^2))),
tb5= logit_thetab_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV1^3))),
tb6= logit_thetab_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV2^3))),
e0 = logit_eps_mn + abs(logit_eps_mn*cvar),
e1 = logit_eps_mn/max(abs(data$ENV1)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1))),
e2 = logit_eps_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1^2))),
e3 = logit_eps_mn/max(abs(data$ENV2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2))), 
e4 = logit_eps_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2^2))), 
e5= logit_eps_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1^3))),
e6= logit_eps_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2^3)))
)

par_initStep =lapply(names(pars), function(x){
max(abs(as.numeric(pars[x])-as.numeric(par_lo[x])), abs(as.numeric(pars[x])-as.numeric(par_hi[x])))
})
names(par_initStep) = names(pars)



#ENV1 = data$ENV1
#ENV2 =  data$ENV2
#save(ENV1, ENV2, file ="ENV.RData")

testBounds =
function(pars, data)
{
    res = TRUE
    ENV1 = data$ENV1
    ENV2 = data$ENV2

 
    logit_alphab= pars$ab0 + pars$ab1*ENV1 + pars$ab2*ENV1^2 + pars$ab3*ENV2 + pars$ab4*ENV2^2 + pars$ab5*ENV1^3 + pars$ab6*ENV2^3
    logit_alphat= pars$at0 + pars$at1*ENV1 + pars$at2*ENV1^2 + pars$at3*ENV2 + pars$at4*ENV2^2 + pars$at5*ENV1^3 + pars$at6*ENV2^3
    logit_betab= pars$bb0 + pars$bb1*ENV1 + pars$bb2*ENV1^2 + pars$bb3*ENV2 + pars$bb4*ENV2^2 + pars$bb5*ENV1^3 + pars$bb6*ENV2^3
    logit_betat= pars$bt0 + pars$bt1*ENV1 + pars$bt2*ENV1^2 + pars$bt3*ENV2 + pars$bt4*ENV2^2 + pars$bt5*ENV1^3 + pars$bt6*ENV2^3
    logit_thetab = pars$tb0 + pars$tb1*ENV1 + pars$tb2*ENV1^2 + pars$tb3*ENV2 + pars$tb4*ENV2^2 + pars$tb5*ENV1^3 + pars$tb6*ENV2^3
    logit_thetat= pars$tt0 + pars$tt1*ENV1 + pars$tt2*ENV1^2 + pars$tt3*ENV2 + pars$tt4*ENV2^2 + pars$tt5*ENV1^3 +pars$ tt6*ENV2^3
    logit_eps= pars$e0  + pars$e1*ENV1  + pars$e2*ENV1^2 + pars$e3*ENV2 + pars$e4*ENV2^2 + pars$e5*ENV1^3 + pars$e6*ENV2^3
   
    alphab = exp(logit_alphab)/(1+exp(logit_alphab))
    alphat = exp(logit_alphat)/(1+exp(logit_alphat))
    betab = exp(logit_betab)/(1+exp(logit_betab))
    betat = exp(logit_betat)/(1+exp(logit_betat))
    thetab = exp(logit_thetab)/(1+exp(logit_thetab))
    thetat = exp(logit_thetat)/(1+exp(logit_thetat))
    eps = exp(logit_eps)/(1+exp(logit_eps))

    
    if(sum(alphab==0)>0 | sum(alphab==1)>0 ) res=FALSE
    if(sum(alphat==0)>0 | sum(alphat==1)>0 ) res=FALSE
    if(sum(betab==0)>0 | sum(betab==1)>0 ) res=FALSE
    if(sum(betat==0)>0 | sum(betat==1)>0 ) res=FALSE
    if(sum(thetab==0)>0 | sum(thetab==1)>0 ) res=FALSE
    if(sum(thetat==0)>0 | sum(thetat==1)>0 ) res=FALSE
    if(sum(eps==0)>0 | sum(eps==1)>0 ) res=FALSE    
    
return(res)
}





