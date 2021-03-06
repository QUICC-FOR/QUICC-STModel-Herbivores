# optimisation
# prerun the model with glm to find some parameter values that are not jointly constrained
# use logits


library(deSolve)
library(rootSolve)


model = function(st0,st1, # vegetation states
ENV1, ENV2,
at0,at1,at2,at3,at4,at5,at6,
ab0,ab1,ab2,ab3,ab4,ab5,ab6,
bt0,bt1,bt2,bt3,bt4,bt5,bt6,
bb0,bb1,bb2,bb3,bb4,bb5,bb6,
tt0,tt1,tt2,tt3,tt4,tt5,tt6,
tb0,tb1,tb2,tb3,tb4,tb5,tb6,
e0,e1,e2,e3,e4,e5,e6) 
{
	lik = numeric(length(st0))
    
    logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV1^2 + ab3*ENV2 + ab4*ENV2^2 + ab5*ENV1^3 + ab6*ENV2^3
    logit_alphat 	= at0 + at1*ENV1 + at2*ENV1^2 + at3*ENV2 + at4*ENV2^2 + at5*ENV1^3 + at6*ENV2^3
    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV1^2 + bb3*ENV2 + bb4*ENV2^2 + bb5*ENV1^3 + bb6*ENV2^3
    logit_betat = bt0 + bt1*ENV1 + bt2*ENV1^2 + bt3*ENV2 + bt4*ENV2^2 + bt5*ENV1^3 + bt6*ENV2^3
    logit_thetab= tb0 + tb1*ENV1 + tb2*ENV1^2 + tb3*ENV2 + tb4*ENV2^2 + tb5*ENV1^3 + tb6*ENV2^3
    logit_thetat= tt0 + tt1*ENV1 + tt2*ENV1^2 + tt3*ENV2 + tt4*ENV2^2 + tt5*ENV1^3 + tt6*ENV2^3
    logit_eps 	= e0  + e1*ENV1  + e2*ENV1^2 + e3*ENV2 + e4*ENV2^2 + e5*ENV1^3 + e6*ENV2^3


    alphab = exp(logit_alphab)/(1+exp(logit_alphab))
    alphat = exp(logit_alphat)/(1+exp(logit_alphat))
    betab = exp(logit_betab)/(1+exp(logit_betab))
    betat = exp(logit_betat)/(1+exp(logit_betat))
    thetab = exp(logit_thetab)/(1+exp(logit_thetab))
    thetat = exp(logit_thetat)/(1+exp(logit_thetat))
    eps = exp(logit_eps)/(1+exp(logit_eps))
    
##
library(parallel)
E_TBMR = mclapply(1:length(st0), function(i){

#pars
pars_intern = list(alphati=alphat[i],
alphabi=alphab[i],
betati = betat[i],
betabi = betab[i],
thetati = thetat[i], 
thetabi = thetab[i],
epsi = eps[i])

##model
model_veget = function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
	
		R=1-TT-B-M
        
	    # vegetation dynamic
	    R_T = alphati*(M+TT)*(1 - alphabi*(M+B))
	    R_B = alphabi*(M+B)*(1 - alphati*(M+TT))
	    R_M = alphati*(M+TT)*alphabi*(M+B)
	    T_M = betabi*(B+M)
	    B_M = betati*(TT+M)
        dTT = R*R_T + M*thetati - TT*(epsi + T_M)
        dB = R*R_B + M*thetabi - B*(epsi + B_M)
        dM = R*R_M + TT*T_M + B*B_M - M*(epsi + thetati + thetabi)
                
                    
		return(list(c(dTT, dB, dM)))
	})
}

##init
T0 = c(TT=0.08, B=0.34, M=.32)

##stode
out = stode(func=model_veget, y=T0, parms = pars_intern, positive=TRUE, maxiter=5000)#, method="iteration")

##

ET = out$y["TT"]
EB = out$y["B"]
EM = out$y["M"]

return(list(ET, EB, EM))
})

steady = do.call(rbind, lapply(E_TBMR, unlist))

ET = steady[,"TT"]
EB = steady[,"B"]
EM = steady[,"M"]

##    
	# Compute the likelihood of observations
	lik[st0 == "B" & st1 == "M"] = (betat*(ET+EM))[st0 == "B" & st1 == "M"] 
	lik[st0 == "B" & st1 == "R"] = eps[st0 == "B" & st1 == "R"] 	
	lik[st0 == "B" & st1 == "B"] = (1 - eps - betat*(ET+EM))[st0 == "B" & st1 == "B"]

	lik[st0 == "T" & st1 == "T"] = (1 - eps - betab*(EB+EM))[st0 == "T" & st1 == "T"] 	
	lik[st0 == "T" & st1 == "M"] = (betab*(EB+EM))[st0 == "T" & st1 == "M"] 			
	lik[st0 == "T" & st1 == "R"] = eps[st0 == "T" & st1 == "R"] 		
	
	lik[st0 == "M" & st1 == "B"] = thetab[st0 == "M" & st1 == "B"]	
	lik[st0 == "M" & st1 == "T"] = thetat[st0 == "M" & st1 == "T"] 	
	lik[st0 == "M" & st1 == "M"] = (1 - eps - thetab - thetat)[st0 == "M" & st1 == "M"] 			
	lik[st0 == "M" & st1 == "R"] = eps[st0 == "M" & st1 == "R"] 
	
	phib = alphab*(EM + EB)*(1-alphat*(ET+EM))
	phit = alphat*(EM + ET)*(1-alphab*(EB+EM))
	phim = alphab*(EM + EB)*alphat*(EM + ET)
	lik[st0 == "R" & st1 == "B"] = phib[st0 == "R" & st1 == "B"] 	
	lik[st0 == "R" & st1 == "T"] = phit[st0 == "R" & st1 == "T"]	
	lik[st0 == "R" & st1 == "M"] = phim[st0 == "R" & st1 == "M"] 			
	lik[st0 == "R" & st1 == "R"] = (1 - phib - phit - phim)[st0 == "R" & st1 == "R"] 
#    lik[lik==0]=NA
#	if(sum(lik==0)>0) cat("lik=0\n")
#	lik[lik==0] = 1/1e99
    
	return(lik)
}

#-----------------------------------------------------------------------------

PDF = function(st1,lik) log(lik)

#-----------------------------------------------------------------------------



