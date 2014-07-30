rm(list=ls())
data = as.data.frame(read.table("../data/data_categorical.txt"))
pred = read.table("../data/data_pred_states.txt")

var = list()
var$t0 = "t0"
var$t1 = "t1"
var$ENV = "ENV"
var$lik = "predicted"
var$EC = "EC"
var$ED = "ED"
var$EM = "EM"

test = numeric(nrow(data))
test[data$t0 == "D" & data$t1 == "C"] = 1
test[data$t0 == "C" & data$t1 == "D"] = 1

data = subset(data, test!=1)
pred = subset(pred, test!=1)
data$ENV = data$av_annual_mean_tp
data$EC = pred[,1]
data$ED = pred[,2]
data$EM = pred[,3] 


# Evaluate initial parameter values
transitions = paste(data$t0,data$t1,sep = "")
sum_transitions = table(transitions)
tot_transitions = table(data$t0)

eps_mn = (sum_transitions[3]+sum_transitions[9])/2/length(transitions)

thetac_mn = sum_transitions[6]/tot_transitions[3]
thetad_mn = sum_transitions[7]/tot_transitions[2]

betac_mn = sum_transitions[7]*(tot_transitions[1]+tot_transitions[3])/tot_transitions[3]/sum(tot_transitions)
betad_mn = sum_transitions[7]*(tot_transitions[3]+tot_transitions[3])/tot_transitions[3]/sum(tot_transitions)

phic_mn = sum_transitions[10]/tot_transitions[4]
phid_mn = sum_transitions[11]/tot_transitions[4]

logit_eps_mn = log(eps_mn/(1-eps_mn))

logit_thetac_mn = log(thetac_mn/(1-thetac_mn))
logit_thetad_mn = log(thetad_mn/(1-thetad_mn))
 
logit_betac_mn = log(betac_mn/(1-betac_mn))
logit_betad_mn = log(betad_mn/(1-betad_mn))

# List initial parameters
par = list(ac0 =-5,
ac1 = 0,
ac2 = 0, 
ad0 = -5,
ad1 = 0,
ad2 = 0,
bc0 = logit_betac_mn,
bc1 = 0,
bc2 = 0,
bd0 = logit_betad_mn,
bd1 = 0,
bd2 = 0,
td0 = logit_thetad_mn,
td1 = 0,
td2 = 0,
tc0 = logit_thetac_mn,
tc1 = 0,
tc2 = 0,
e0 = logit_eps_mn,
e1 = 0,
e2 = 0)

par_initStep =lapply(par, function(x){x+0.01})


model = function(t0,t1,ED,EC,EM,ENV,ac0,ac1,ac2,ad0,ad1,ad2,bc0,bc1,bc2,bd0,bd1,bd2,tc0,tc1,tc2,td0,td1,td2,e0,e1,e2) 
{
	lik = numeric(length(t0))

	# Compute the logit
	logit_alphac 	= ac0 + ac1*ENV + ac2*ENV^2
	logit_alphad 	= ad0 + ad1*ENV + ad2*ENV^2
	logit_betac 	= bc0 + bc1*ENV + bc2*ENV^2
	logit_betad 	= bd0 + bd1*ENV + bd2*ENV^2
	logit_thetac	= tc0 + tc1*ENV + tc2*ENV^2
	logit_thetad	= td0 + td1*ENV + td2*ENV^2
	logit_eps 	= e0  + e1*ENV  + e2*ENV^2

	# Back transform into probabilities
	alphac = exp(logit_alphac)/(1+exp(logit_alphac))
	alphad = exp(logit_alphad)/(1+exp(logit_alphad))
	betac = exp(logit_betac)/(1+exp(logit_betac))*(EC+EM)
	betad = exp(logit_betad)/(1+exp(logit_betad))*(ED+EM)
	thetac = exp(logit_thetac)/(1+exp(logit_thetac))
	thetad = exp(logit_thetad)/(1+exp(logit_thetad))
	eps = exp(logit_eps)/(1+exp(logit_eps))
	phic = alphac*(EM + EC)*(1-alphad*(ED+EM))
	phid = alphad*(EM + ED)*(1-alphac*(EC+EM))
	phim = phic*phid

	# Compute the likelihood of observations
	lik[t0 == "C" & t1 == "M"] = betad[t0 == "C" & t1 == "M"] 
	lik[t0 == "C" & t1 == "T"] = eps[t0 == "C" & t1 == "T"] 	
	lik[t0 == "C" & t1 == "C"] = (1 - eps - betad)[t0 == "C" & t1 == "C"]

	lik[t0 == "D" & t1 == "D"] = (1 - eps - betac)[t0 == "D" & t1 == "D"] 	
	lik[t0 == "D" & t1 == "M"] = betac[t0 == "D" & t1 == "M"] 			
	lik[t0 == "D" & t1 == "T"] = eps[t0 == "D" & t1 == "T"] 		
	
	lik[t0 == "M" & t1 == "C"] = thetac[t0 == "M" & t1 == "C"]	
	lik[t0 == "M" & t1 == "D"] = thetad[t0 == "M" & t1 == "D"] 	
	lik[t0 == "M" & t1 == "M"] = (1 - eps - thetac - thetad)[t0 == "M" & t1 == "M"] 			
	lik[t0 == "M" & t1 == "T"] = eps[t0 == "M" & t1 == "T"] 
	
	lik[t0 == "T" & t1 == "C"] = phic[t0 == "T" & t1 == "C"] 	
	lik[t0 == "T" & t1 == "D"] = phid[t0 == "T" & t1 == "D"]	
	lik[t0 == "T" & t1 == "M"] = phim[t0 == "T" & t1 == "M"] 			
	lik[t0 == "T" & t1 == "T"] = (1 - phic - phid - phim)[t0 == "T" & t1 == "T"] 

	lik[lik==0] = NA

	return(lik)
}

PDF = function(t1,lik) log(lik)


#---

testBounds =
function(par, ENV_min = min(data$ENV), ENV_max = max(data$ENV), interval=0.2)
{
    res = TRUE
    
    for (ENV in seq(ENV_min, ENV_max, interval))
    {
    logit_alphab	= par$ac0 + par$ac1*ENV + par$ac2*ENV^2
    logit_alphat 	= par$ad0 + par$ad1*ENV + par$ad2*ENV^2
    logit_betab 	= par$bc0 + par$bc1*ENV + par$bc2*ENV^2
    logit_betat 	= par$bd0 + par$bd1*ENV + par$bd2*ENV^2
    logit_thetat	= par$td0 + par$td1*ENV + par$td2*ENV^2
    logit_thetab	= par$tc0 + par$tc1*ENV + par$tc2*ENV^2
    logit_eps 	= par$e0 + par$e1*ENV + par$e2*ENV^2

    alphab = exp(logit_alphab)/(1+exp(logit_alphab))
    alphat = exp(logit_alphat)/(1+exp(logit_alphat))
    betab = exp(logit_betab)/(1+exp(logit_betab))
    betat = exp(logit_betat)/(1+exp(logit_betat))
    thetab = exp(logit_thetab)/(1+exp(logit_thetab))
    thetat = exp(logit_thetat)/(1+exp(logit_thetat))
    eps = exp(logit_eps)/(1+exp(logit_eps))

    
    if(alphab==0 | alphab==1) res=FALSE
    if(alphat==0 | alphat==1) res=FALSE
    if(betab==0 | betab==1) res=FALSE
    if(betat==0 | betat==1) res=FALSE
    if(thetab==0 | thetab==1) res=FALSE
    if(thetat==0 | thetat==1) res=FALSE
    if(eps==0 | eps==1) res=FALSE    
    }
    
return(res)
}




#----


source("fit_fonctions/analyze_function.R")
source("fit_fonctions/likdisplay_ibou.R")
source("fit_fonctions/support_limits_ibou.R")
source("fit_fonctions/likeli.R")



source("fit_fonctions/anneal_ibou.R")

test1= anneal(model = model, par = par, var = var, source_data = data, 
dep_var = "t1", pdf = PDF, par_initStep = par_initStep, par_testFct = testBounds, max_iter = 2500, initial_temp = 0.5, note = "", progress = TRUE, display=FALSE, support = FALSE)
#temp_red = 0.95, ns = 10, nt = 50, min_change = 0.01, min_drops = 20, delta = 100, slimit = 2, 


# List initial parameters

par_lo = list(ac0 = -10,
ac1 = -20,
ac2 = -20, 
ad0 = -10,
ad1 = -10,
ad2 = -20,
bc0 = -10,
bc1 = -10,
bc2 = -25,
bd0 = -30,
bd1 = -10,
bd2 = -10,
td0 = -10,
td1 = -10,
td2 = -25,
tc0 = -10,
tc1 = -10,
tc2 = -20,
e0 = -10,
e1 = -10,
e2 = -10)

par_hi = list(ac0 = 50,
ac1 = 15,
ac2 = 30, 
ad0 = 10,
ad1 = 10,
ad2 = 10,
bc0 = 10,
bc1 = 10,
bc2 = 20,
bd0 = 10,
bd1 = 10,
bd2 = 50,
td0 = 10,
td1 = 10,
td2 = 20,
tc0 = 10,
tc1 = 10,
tc2 = 20,
e0 = 10,
e1 = 10,
e2 = 10)

source("fit_fonctions/anneal.R")


# Maximum likelihood estimation
test2 = anneal(model = model, par = par, var = var, source_data = data, 
par_lo = par_lo, par_hi = par_hi, dep_var = "t1", pdf = PDF, 
max_iter = 2500, hessian = FALSE, initial_temp = 0.5)

cbind(test1$max_likeli, test2$max_likeli)
cbind(unlist(test1$best_pars), unlist(test2$best_pars))
#----
