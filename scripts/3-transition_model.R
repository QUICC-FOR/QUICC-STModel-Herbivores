# optimisation
# prerun the model with glm to find some parameter values that are not jointly constrained
# use logits


model = function(t0,t1, # vegetation states
ET,EB,EM,
ENV,
at0,at1,at2,at3,
ab0,ab1,ab2,ab3,
bt0,bt1,bt2,bt3,
bb0,bb1,bb2,bb3,
tt0,tt1,tt2,tt3,
tb0,tb1,tb2,tb3,
e0,e1,e2,e3,
Ha, Hv # herbivore densities at t0
) 
{
	lik = numeric(length(t0))


    logit_alphab 	= ab0 + ab1*ENV + ab2*ENV^2 + ab3*ENV^3
    logit_alphat 	= at0 + at1*ENV + at2*ENV^2 + at3*ENV^3
    logit_betab 	= bb0 + bb1*ENV + bb2*ENV^2 + bb3*ENV^3
    logit_betat 	= bt0 + bt1*ENV + bt2*ENV^2 + bt3*ENV^3
    logit_thetab	= tb0 + tb1*ENV + tb2*ENV^2 + tb3*ENV^3
    logit_thetat	= tt0 + tt1*ENV + tt2*ENV^2 + tt3*ENV^3
    logit_eps 	= e0  + e1*ENV  + e2*ENV^2 + e3*ENV^3

    alphab = exp(logit_alphab)/(1+exp(logit_alphab))
    alphat = exp(logit_alphat)/(1+exp(logit_alphat))
    betab = exp(logit_betab)/(1+exp(logit_betab))
    betat = exp(logit_betat)/(1+exp(logit_betat))
    thetab = exp(logit_thetab)/(1+exp(logit_thetab))
    thetat = exp(logit_thetat)/(1+exp(logit_thetat))
    eps = exp(logit_eps)/(1+exp(logit_eps))

	
	# impact of the herbivore
    # herbivores competition
	k = exp(-log(1/k0)*Hv/Ha)
	
	# state t0
	R = TT = B = M = 0
	if(t0=="R") R=1
	if(t0=="T") TT=1
	if(t0=="B") B=1
	if(t0=="M") M=1
		
	# moose and B
	ma_B = mas/(1+(mas/ma0-1)*B)
	
	# intakes
    if(Ha!=0) {
    Fa = u*k*R/Ha
    Ga = k*(v*TT+w*B+x*M)/Ha
    }else{
    Fa=0
    Ga=0
    }
    if(Hv!=0) {
    Fv = u*(1-k)*R/Hv
    Gv = (1-k)*(v*TT+w*B+x*M)/Hv
    }else{
    Fv=0
    Gv=0
    } 
        
    Ia1 = taua*Fa/(mua+Fa)
    Ia2 = (taua*Ga/(nua + Ga) )* (phia + (1-phia)/( 1 + exp(ra*( taua*Fa/(mua + Fa) -pa) ) ) )
    Iv1 = tauv*Fv/(muv+Fv)
    Iv2 = (tauv*Gv/(nuv + Gv) )* (phiv + (1-phiv)/( 1 + exp(rv*( tauv*Fv/(muv + Fv) -pv) ) ) )        
        
    # herbivore pressure
        
    if(R!=0) {
    PRB = (Ha*Ia1*kappaba + Hv*Iv1*kappabv)/(u*R)
    PRT = (Ha*Ia1*kappaba + Hv*Iv1*kappabv)/(u*R)
    }else{PRB = PRT =0}
        
    if(R!=1) 
    {
    sumOHa = omegatHa * T + omegabHa * B + (1-omegatHa-omegabHa) * M
    OtHa = omegatHa*T /sumOHa
    ObHa = omegabHa*B /sumOHa
    OmHa = (1-omegatHa-omegabHa)*M /sumOHa

    sumOHv = omegatHv * T + omegabHv * B + (1-omegatHv-omegabHv) * M
    OtHv = omegatHv*T /sumOHv
    ObHv = omegabHv*B /sumOHv
    OmHv = (1-omegatHv-omegabHv)*M /sumOHv
    
    if(TT!=0) {
    PTB = (OtHa*Ha*Ia2*kappaba + OtHv*Hv*Iv2*kappabv)/(v*TT)
    PTT = (OtHa*Ha*Ia2*kappata + OtHv*Hv*Iv2*kappatv)/(v*TT)
    }else{PTB = PTT =0}
    if(B!=0) {
    PBB = (ObHa*Ha*Ia2*kappaba + ObHv*Hv*Ia2*kappabv)/(w*B) 
    PBT = (ObHa*Ha*Ia2*kappata + ObHv*Hv*Ia2*kappatv)/(w*B) 
    }else{PBB = PBT =0}
    if(M!=0) {
    PMB = (OmHa*Ha*Ia2*kappaba + OmHv*Hv*Iv2*kappabv)/(x*M) 
    PMT = (OmHa*Ha*Ia2*kappata + OmHv*Hv*Iv2*kappatv)/(x*M) 
    }else{PMB = PMT =0}

    }else{
    PTB = PTT =0
    PBB = PBT =0
    PMB = PMT =0
    }
   
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
	lik[t0 == "R" & t1 == "TT"] = (alphad*(EM + ED)*(1-alphac*(EC+EM)))[t0 == "R" & t1 == "TT"]	
	lik[t0 == "R" & t1 == "M"] = (alphab_h*(EM + EB)*alphad*(EM + ED))[t0 == "R" & t1 == "M"] 			
	lik[t0 == "R" & t1 == "R"] = (1 - phic - phid - phim)[t0 == "R" & t1 == "R"] 

	lik[lik==0] = NA
	

	return(lik)
}

PDF = function(t1,lik) log(lik)





