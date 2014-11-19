# optimisation
# prerun the model with glm to find some parameter values that are not jointly constrained
# use logits


model = function(st0,st1, # vegetation states
ET,EB,EM,
ENV1, ENV2,
vegetPars,
Ha, Hv # herbivore densities at t0
) 
{
	lik = numeric(length(st0))

    logitModel = transition_model(vegetPars)
    logit_alphab 	= logitModel$logit_alphab
    logit_alphat 	= logitModel$logit_alphat
    logit_betab 	= logitModel$logit_betab
    logit_betat 	= logitModel$logit_betat
    logit_thetab	= logitModel$logit_thetab
    logit_thetat	= logitModel$logit_thetat
    logit_eps 	= logitModel$logit_eps

    alphab = exp(logit_alphab)/(1+exp(logit_alphab))
    alphat = exp(logit_alphat)/(1+exp(logit_alphat))
    betab = exp(logit_betab)/(1+exp(logit_betab))
    betat = exp(logit_betat)/(1+exp(logit_betat))
    thetab = exp(logit_thetab)/(1+exp(logit_thetab))
    thetat = exp(logit_thetat)/(1+exp(logit_thetat))
    eps = exp(logit_eps)/(1+exp(logit_eps))

    par_h = herbivory(alphat=alphat, alphab=alphab, betat=betat, betab=betab, thetat=thetat, thetab=thetab, ET=ET, EB=EB,EM=EM, Hv=Hv, Ha=Ha)
#    par_h =list(alphat_h = alphat, alphab_h = alphab, betat_h = betat, betab_h = betab, thetat_h = thetat, thetab_h =thetab)
    
	# Compute the likelihood of observations
	lik[st0 == "B" & st1 == "M"] = (par_h$betat_h*(ET+EM))[st0 == "B" & st1 == "M"] 
	lik[st0 == "B" & st1 == "R"] = eps[st0 == "B" & st1 == "R"] 	
	lik[st0 == "B" & st1 == "B"] = (1 - eps - par_h$betat_h*(ET+EM))[st0 == "B" & st1 == "B"]

	lik[st0 == "T" & st1 == "T"] = (1 - eps - par_h$betab_h*(EB+EM))[st0 == "T" & st1 == "T"] 	
	lik[st0 == "T" & st1 == "M"] = (par_h$betab_h*(EB+EM))[st0 == "T" & st1 == "M"] 			
	lik[st0 == "T" & st1 == "R"] = eps[st0 == "T" & st1 == "R"] 		
	
	lik[st0 == "M" & st1 == "B"] = par_h$thetab_h[st0 == "M" & st1 == "B"]	
	lik[st0 == "M" & st1 == "T"] = par_h$thetat_h[st0 == "M" & st1 == "T"] 	
	lik[st0 == "M" & st1 == "M"] = (1 - eps - par_h$thetab_h - par_h$thetat_h)[st0 == "M" & st1 == "M"] 			
	lik[st0 == "M" & st1 == "R"] = eps[st0 == "M" & st1 == "R"] 
	
	phib = par_h$alphab_h*(EM + EB)*(1-par_h$alphat_h*(ET+EM))
	phit = par_h$alphat_h*(EM + ET)*(1-par_h$alphab_h*(EB+EM))
	phim = par_h$alphab_h*(EM + EB)*par_h$alphat_h*(EM + ET)
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

    herbivory <- function(alphat, alphab, betat, betab,thetat, thetab, ET, EB, EM, Hv, Ha)
    {
    
#    Hv -> 1km2
#    Ha -> 10km2
#    
	#params fixes
	k0 = 0.5
	u = 0.5*1000; v=200; w = 100; x=200

    omegatHa = 0.3
    omegabHa = 0.35

    omegatHv = 0.325
    omegabHv = 0.35

    kappata = 0.25
    kappaba = 0.07

    kappatv = 0.265
    kappabv = 0.175
    
    # moose
    dva = 50000
    Ha = Ha*dva
    ua = u*dva ; va = u*dva ;  wa= u*dva ; xa= u*dva

    za = 1
    ca = 0.79
    ma0 = 0.25
    mas = 0.25
    taua = 13
    mua = (taua/0.07)*(dva/100)/sqrt(314)
    rhoa = 6
    nua = (rhoa/0.15)*(dva/100)/sqrt(314)
    phia = 0.5
    pa = 1.7
    ra = 5

    #deer
    dvv = 1000 # ha
    Hv = Hv*dvv
    uv = u*dvv ; vv = u*dvv ;  wv= u*dvv ; xv= u*dvv

    zv = 1
    cv = 0.5
    mv0 = 0.2
    mvs = 0.5
    tauv = 17.5
    muv = (tauv/0.07)*(dvv/100)/sqrt(518)
    rhov = 4.5
    nuv = (rhov/0.15)*(dvv/100)/sqrt(518)
    pv = 1.87
    phiv = 0.57
    rv = 5

	# impact of the herbivore
    # herbivores competition
	k = exp(-log(1/k0)*(dva/dvv)*Hv/Ha)
	
	# neighborgh t0
	ER = 1-ET -EB-EM
	
		
	# shelter effect
	ma_B = mas/(1+(mas/ma0-1)*(EB+EM))
	mv_B = mvs/(1+(mvs/mv0-1)*(EB+EM))
	
	
	# intakes
    Fa = ua*k*ER/Ha
    Ga = k*(va*ET+wa*EB+xa*EM)/Ha
    Fa[is.nan(Fa)] = 0
    Ga[is.nan(Ga)] = 0
    Fa[is.infinite(Fa)] = 0
    Ga[is.infinite(Ga)] = 0
    
    Fv = uv*(1-k)*ER/Hv
    Gv = (1-k)*(vv*ET+wv*EB+xv*EM)/Hv
    Fv[is.nan(Fv)] = 0
    Gv[is.nan(Gv)] = 0
    Fv[is.infinite(Fv)] = 0
    Gv[is.infinite(Gv)] = 0
        
    Ia1 = taua*Fa/(mua+Fa)
    Ia2 = (rhoa*Ga/(nua + Ga) )* (phia + (1-phia)/( 1 + exp(ra*( taua*Fa/(mua + Fa) -pa) ) ) )
    Iv1 = tauv*Fv/(muv+Fv)
    Iv2 = (rhov*Gv/(nuv + Gv) )* (phiv + (1-phiv)/( 1 + exp(rv*( tauv*Fv/(muv + Fv) -pv) ) ) )        
        
    # herbivore pressure
        
    PRB = Ha*Ia1*kappaba/(ua*ER) + Hv*Iv1*kappabv/(uv*ER)
    PRT = Ha*Ia1*kappaba/(ua*ER) + Hv*Iv1*kappabv/(uv*ER)
    PRB[is.nan(PRB)] = 0
    PRT[is.nan(PRT)] = 0
    PRB[is.infinite(PRB)] = 0
    PRT[is.infinite(PRT)] = 0

        
    sumOHa = omegatHa * ET + omegabHa * EB + (1-omegatHa-omegabHa) * EM
    OtHa = omegatHa*ET /sumOHa
    ObHa = omegabHa*EB /sumOHa
    OmHa = (1-omegatHa-omegabHa)*EM /sumOHa

    sumOHv = omegatHv * ET + omegabHv * EB + (1-omegatHv-omegabHv) * EM
    OtHv = omegatHv*ET /sumOHv
    ObHv = omegabHv*EB /sumOHv
    OmHv = (1-omegatHv-omegabHv)*EM /sumOHv
    
    PTB = (OtHa*Ha*Ia2*kappaba)/(va*ET) + OtHv*Hv*Iv2*kappabv/(vv*ET)
    PTT = (OtHa*Ha*Ia2*kappata)/(va*ET) + OtHv*Hv*Iv2*kappatv/(vv*ET)

    PBB = (ObHa*Ha*Ia2*kappaba)/(wa*EB) + ObHv*Hv*Ia2*kappabv/(wv*EB) 
    PBT = (ObHa*Ha*Ia2*kappata)/(wa*EB) + ObHv*Hv*Ia2*kappatv/(wv*EB) 

    PMB = (OmHa*Ha*Ia2*kappaba)/(xa*EM) + OmHv*Hv*Iv2*kappabv/(xv*EM) 
    PMT = (OmHa*Ha*Ia2*kappata)/(xa*EM) + OmHv*Hv*Iv2*kappatv/(xv*EM) 

    PTB[is.nan(PTB)] = 0
    PTT[is.nan(PTT)] = 0
    PBB[is.nan(PBB)] = 0
    PMB[is.nan(PMB)] = 0
    PBT[is.nan(PBT)] = 0
    PMT[is.nan(PMT)] = 0
    
    PTB[is.infinite(PTB)] = 0
    PTT[is.infinite(PTT)] = 0
    PBB[is.infinite(PBB)] = 0
    PMB[is.infinite(PMB)] = 0
    PBT[is.infinite(PBT)] = 0
    PMT[is.infinite(PMT)] = 0


   
    # compute modified alpha, theta and beta
    alphab_h = alphab*(1-PRB)
    alphat_h = alphat*(1-PRT)
    
    thetab_h = thetab*(1-PMB)
    thetat_h = thetat*(1-PMT)
    
    betab_h = (1/2)*(1 + betab * cos(pi*PTB) + (1-betab)*cos(pi*(1+PTT)) ) 
    betat_h = (1/2)*(1 + betat * cos(pi*PBT) + (1-betat)*cos(pi*(1+PBB)) )

    return(list(alphat_h=alphat_h, alphab_h=alphab_h,betat_h=betat_h, betab_h=betab_h, thetat_h=thetat_h, thetab_h=thetab_h))
    }




