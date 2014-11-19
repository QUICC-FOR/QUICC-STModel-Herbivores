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
e0,e1,e2,e3,e4,e5,e6, 
Ha, Hv) 
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
    
    source("../../model/params_herbivores.r")


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
epsi = eps[i], 
Hai = Ha[i], 
Hvi = Hv[i], 
k0 = k0, 
ua=ua,
uv =uv,
va=va,
vv=vv,
wa=wa,
wv=wv,
xa=xa,
xv=xv, 
omegatHa = omegatHa, 
omegatHv = omegatHv, 
kappata = kappata, 
kappaba = kappaba, 
kappatv = kappatv, 
kappabv = kappabv, 
dva = dva, 
dvv = dvv, 
za = za, 
ca =ca, 
ma0 = ma0, 
mas = mas, 
taua = taua, 
mua = mua, 
rhoa = rhoa, 
nua = nua, 
phia = phia, 
pa = pa, 
ra = ra, 
zv = zv, 
cv = cv, 
mv0 = mv0, 
mvs = mvs, 
tauv = tauv, 
muv = muv, 
rhov = rhov, 
nuv = nuv, 
pv = pv, 
phiv = phiv, 
rv=rv
)

##model
model_herbi = function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
	
		R=1-TT-B-M
		
			# impact of the herbivore

    # herbivores competition
	k = exp(-log(1/k0)*(dva/dvv)*Hvi/Hai)
	
		
	# shelter effect
	ma_B = mas/(1+(mas/ma0-1)*(B+M))
	mv_B = mvs/(1+(mvs/mv0-1)*(B+M))
	
	
	# intakes
    Fa = ua*k*R/Hai
    Ga = k*(va*TT+wa*B+xa*M)/Hai
    Fa[is.nan(Fa)] = 0
    Ga[is.nan(Ga)] = 0
    Fa[is.infinite(Fa)] = 0
    Ga[is.infinite(Ga)] = 0
    
    Fv = uv*(1-k)*R/Hvi
    Gv = (1-k)*(vv*TT+wv*B+xv*M)/Hvi
    Fv[is.nan(Fv)] = 0
    Gv[is.nan(Gv)] = 0
    Fv[is.infinite(Fv)] = 0
    Gv[is.infinite(Gv)] = 0
        
    Ia1 = taua*Fa/(mua+Fa)
    Ia2 = (rhoa*Ga/(nua + Ga) )* (phia + (1-phia)/( 1 + exp(ra*( taua*Fa/(mua + Fa) -pa) ) ) )
    Iv1 = tauv*Fv/(muv+Fv)
    Iv2 = (rhov*Gv/(nuv + Gv) )* (phiv + (1-phiv)/( 1 + exp(rv*( tauv*Fv/(muv + Fv) -pv) ) ) )        
        
    # herbivore pressure
        
    PRB = Hai*Ia1*kappaba/(ua*R) + Hvi*Iv1*kappabv/(uv*R)
    PRT = Hai*Ia1*kappaba/(ua*R) + Hvi*Iv1*kappabv/(uv*R)
    PRB[is.nan(PRB)] = 0
    PRT[is.nan(PRT)] = 0
    PRB[is.infinite(PRB)] = 0
    PRT[is.infinite(PRT)] = 0

        
    sumOHa = omegatHa * TT + omegabHa * B + (1-omegatHa-omegabHa) * M
    OtHa = omegatHa*TT /sumOHa
    ObHa = omegabHa*B /sumOHa
    OmHa = (1-omegatHa-omegabHa)*M /sumOHa

    sumOHv = omegatHv * TT + omegabHv * B + (1-omegatHv-omegabHv) * M
    OtHv = omegatHv*TT /sumOHv
    ObHv = omegabHv*B /sumOHv
    OmHv = (1-omegatHv-omegabHv)*M /sumOHv
    
    PTB = (OtHa*Hai*Ia2*kappaba)/(va*TT) + OtHv*Hvi*Iv2*kappabv/(vv*TT)
    PTT = (OtHa*Hai*Ia2*kappata)/(va*TT) + OtHv*Hvi*Iv2*kappatv/(vv*TT)

    PBB = (ObHa*Hai*Ia2*kappaba)/(wa*B) + ObHv*Hvi*Ia2*kappabv/(wv*B) 
    PBT = (ObHa*Hai*Ia2*kappata)/(wa*B) + ObHv*Hvi*Ia2*kappatv/(wv*B) 

    PMB = (OmHa*Hai*Ia2*kappaba)/(xa*M) + OmHv*Hvi*Iv2*kappabv/(xv*M) 
    PMT = (OmHa*Hai*Ia2*kappata)/(xa*M) + OmHv*Hvi*Iv2*kappatv/(xv*M) 

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
    alphab_h = alphabi*(1-PRB)
    alphat_h = alphati*(1-PRT)
    
    thetab_h = thetabi*(1-PMB)
    thetat_h = thetati*(1-PMT)
    
    betab_h = (1/2)*(1 + betabi * cos(pi*PTB) + (1-betabi)*cos(pi*(1+PTT)) ) 
    betat_h = (1/2)*(1 + betati * cos(pi*PBT) + (1-betati)*cos(pi*(1+PBB)) )

        
	    # vegetation dynamic
	    R_T = alphat_h*(M+TT)*(1 - alphab_h*(M+B))
	    R_B = alphab_h*(M+B)*(1 - alphat_h*(M+TT))
	    R_M = alphat_h*(M+TT)*alphab_h*(M+B)
	    T_M = betab_h*(B+M)
	    B_M = betat_h*(TT+M)
        dTT = R*R_T + M*thetat_h - TT*(epsi + T_M)
        dB = R*R_B + M*thetab_h - B*(epsi + B_M)
        dM = R*R_M + TT*T_M + B*B_M - M*(epsi + thetat_h + thetab_h)
                
                    
		return(list(c(dTT, dB, dM)))
	})
}

##init
T0 = c(TT=0.08, B=0.34, M=.32)

##stode
out = stode(func=model_herbi, y=T0, parms = pars_intern, positive=TRUE, maxiter=5000)#, method="iteration")

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



