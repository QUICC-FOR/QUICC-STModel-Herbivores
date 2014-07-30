# optimisation
# prerun the model with glm to find some parameter values that are not jointly constrained
# use logits


model = function(t0,t1, # vegetation states
ET,EB,EM,
ENV1, ENV2,
at0,at1,at2,at3,at4,
ab0,ab1,ab2,ab3,ab4,
bt0,bt1,bt2,bt3,bt4,
bb0,bb1,bb2,bb3,bb4,
tt0,tt1,tt2,tt3,tt4,
tb0,tb1,tb2,tb3,tb4,
e0,e1,e2,e3,e4,
Ha, Hv # herbivore densities at t0
) 
{
	lik = numeric(length(t0))


    logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV1^2 + ab3*ENV2 + ab4*ENV2^2
    logit_alphat 	= at0 + at1*ENV1 + at2*ENV1^2 + at3*ENV2 + at4*ENV2^2
    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV1^2 + bb3*ENV2 + bb4*ENV2^2
    logit_betat 	= bt0 + bt1*ENV1 + bt2*ENV1^2 + bt3*ENV2 + bt4*ENV2^2
    logit_thetab	= tb0 + tb1*ENV1 + tb2*ENV1^2 + tb3*ENV2 + tb4*ENV2^2
    logit_thetat	= tt0 + tt1*ENV1 + tt2*ENV1^2 + tt3*ENV2 + tt4*ENV2^2
    logit_eps 	= e0  + e1*ENV1  + e2*ENV1^2 + e3*ENV2 + e4*ENV2^2

    alphab = exp(logit_alphab)/(1+exp(logit_alphab))
    alphat = exp(logit_alphat)/(1+exp(logit_alphat))
    betab = exp(logit_betab)/(1+exp(logit_betab))
    betat = exp(logit_betat)/(1+exp(logit_betat))
    thetab = exp(logit_thetab)/(1+exp(logit_thetab))
    thetat = exp(logit_thetat)/(1+exp(logit_thetat))
    eps = exp(logit_eps)/(1+exp(logit_eps))

	#params fixes
	k0 = 0.5

    u = 1000 ; v = 500 ;  w= 200 ; x= 300
    omegatHa = 0.3
    omegabHa = 0.3

    omegatHv = 0.5
    omegabHv = 0

    kappata = 0.5
    kappaba = 0.3

    kappatv = 0.265
    kappabv = 0.175
    # moose
    za = 1
    ca = 0.7
    ma0 = 0.05
    mas = 0.2
    taua = 17
    pa = 6
    mua = taua
    nua = taua
    phia = 0.3
    ra = 5

    #deer
    zv = 1
    cv = 0.7
    mv = 0.2
    tauv = 8.75
    pv = 5
    muv = tauv
    nuv = tauv
    phiv = 0.2
    rv = 5


	# impact of the herbivore
    # herbivores competition
	k = exp(-log(1/k0)*Hv/Ha)
	
	
	# neighborgh t0
	R = 1-ET -EB-EM
	TT = ET
	B = EB
	M = EM

		
	# moose and B
	ma_B = mas/(1+(mas/ma0-1)*B)
	
	# intakes
    Fa = u*k*R/Ha
    Ga = k*(v*TT+w*B+x*M)/Ha
    Fa[is.nan(Fa)] = 0
    Ga[is.nan(Ga)] = 0
    Fa[is.infinite(Fa)] = 0
    Ga[is.infinite(Ga)] = 0
    
    Fv = u*(1-k)*R/Hv
    Gv = (1-k)*(v*TT+w*B+x*M)/Hv
    Fv[is.nan(Fv)] = 0
    Gv[is.nan(Gv)] = 0
    Fv[is.infinite(Fv)] = 0
    Gv[is.infinite(Gv)] = 0
        
    Ia1 = taua*Fa/(mua+Fa)
    Ia2 = (taua*Ga/(nua + Ga) )* (phia + (1-phia)/( 1 + exp(ra*( taua*Fa/(mua + Fa) -pa) ) ) )
    Iv1 = tauv*Fv/(muv+Fv)
    Iv2 = (tauv*Gv/(nuv + Gv) )* (phiv + (1-phiv)/( 1 + exp(rv*( tauv*Fv/(muv + Fv) -pv) ) ) )        
        
    # herbivore pressure
        
    PRB = (Ha*Ia1*kappaba + Hv*Iv1*kappabv)/(u*R)
    PRT = (Ha*Ia1*kappaba + Hv*Iv1*kappabv)/(u*R)
    PRB[is.nan(PRB)] = 0
    PRT[is.nan(PRT)] = 0
    PRB[is.infinite(PRB)] = 0
    PRT[is.infinite(PRT)] = 0

        
    sumOHa = omegatHa * T + omegabHa * B + (1-omegatHa-omegabHa) * M
    OtHa = omegatHa*T /sumOHa
    ObHa = omegabHa*B /sumOHa
    OmHa = (1-omegatHa-omegabHa)*M /sumOHa

    sumOHv = omegatHv * T + omegabHv * B + (1-omegatHv-omegabHv) * M
    OtHv = omegatHv*T /sumOHv
    ObHv = omegabHv*B /sumOHv
    OmHv = (1-omegatHv-omegabHv)*M /sumOHv
    
    PTB = (OtHa*Ha*Ia2*kappaba + OtHv*Hv*Iv2*kappabv)/(v*TT)
    PTT = (OtHa*Ha*Ia2*kappata + OtHv*Hv*Iv2*kappatv)/(v*TT)

    PBB = (ObHa*Ha*Ia2*kappaba + ObHv*Hv*Ia2*kappabv)/(w*B) 
    PBT = (ObHa*Ha*Ia2*kappata + ObHv*Hv*Ia2*kappatv)/(w*B) 

    PMB = (OmHa*Ha*Ia2*kappaba + OmHv*Hv*Iv2*kappabv)/(x*M) 
    PMT = (OmHa*Ha*Ia2*kappata + OmHv*Hv*Iv2*kappatv)/(x*M) 

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

    
	# Compute the likelihood of observations
	lik[t0 == "B" & t1 == "M"] = (betat_h*(ET+EM))[t0 == "B" & t1 == "M"] 
	lik[t0 == "B" & t1 == "R"] = eps[t0 == "B" & t1 == "R"] 	
	lik[t0 == "B" & t1 == "B"] = (1 - eps - betat_h*(ET+EM))[t0 == "B" & t1 == "B"]

	lik[t0 == "TT" & t1 == "TT"] = (1 - eps - betab_h*(EB+EM))[t0 == "TT" & t1 == "TT"] 	
	lik[t0 == "TT" & t1 == "M"] = (betab_h*(ET+EM))[t0 == "TT" & t1 == "M"] 			
	lik[t0 == "TT" & t1 == "R"] = eps[t0 == "TT" & t1 == "R"] 		
	
	lik[t0 == "M" & t1 == "B"] = thetab_h[t0 == "M" & t1 == "B"]	
	lik[t0 == "M" & t1 == "TT"] = thetat_h[t0 == "M" & t1 == "TT"] 	
	lik[t0 == "M" & t1 == "M"] = (1 - eps - thetab_h - thetat_h)[t0 == "M" & t1 == "M"] 			
	lik[t0 == "M" & t1 == "R"] = eps[t0 == "M" & t1 == "R"] 
	
	lik[t0 == "R" & t1 == "B"] = (alphab_h*(EM + EB)*(1-alphat_h*(ET+EM)))[t0 == "R" & t1 == "B"] 	
	lik[t0 == "R" & t1 == "TT"] = (alphat*(EM + ET)*(1-alphab*(EB+EM)))[t0 == "R" & t1 == "TT"]	
	lik[t0 == "R" & t1 == "M"] = (alphab_h*(EM + EB)*alphat*(EM + ET))[t0 == "R" & t1 == "M"] 			
	
	phib = alphab_h*(EM + EB)*(1-alphat_h*(ET+EM))
	phit = alphat*(EM + ET)*(1-alphab*(EB+EM))
	phim = alphab_h*(EM + EB)*alphat*(EM + ET)
	lik[t0 == "R" & t1 == "R"] = (1 - phib - phit - phim)[t0 == "R" & t1 == "R"] 

#	lik[lik==0] = NA
#	
#    print(lik)
	return(lik)
}

PDF = function(t1,lik) log(lik+1)





