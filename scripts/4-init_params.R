# state transitions
data = as.data.frame(read.table("../data/data_reshaped_RBTM_herbivores.txt"))
# neighborhood
pred = read.table("../data/data_pred_states_randomForest.txt")

# remove transitions D->C and C->D
test = numeric(nrow(data))
test[data$t0 == "T" & data$t1 == "B"] = 1
test[data$t0 == "B" & data$t1 == "T"] = 1
data = subset(data, test!=1)
pred = subset(pred, test!=1)

data$ENV1 = data$av_annual_mean_tp * 100
data$ENV2 = data$av_annual_pp/10
data$EB = pred$B
data$ET = pred$T
data$EM = pred$M 
data$Hv = data$deer
data$Ha = data$moose

head(data)

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
par = list(ab0 =logit_alphab_mn, ab1 = 0, ab2 = 0, ab3=0, ab4=0,
at0 = logit_alphat_mn, at1 = 0, at2 = 0, at3=0, at4=0,
bb0 = logit_betab_mn, bb1 = 0, bb2 = 0, bb3=0, bb4=0,
bt0 = logit_betat_mn, bt1 = 0, bt2 = 0, bt3=0, bt4=0,
tt0 = logit_thetat_mn, tt1 = 0, tt2 = 0, tt3 =0, tt4=0,
tb0 = logit_thetab_mn, tb1 = 0, tb2 = 0, tb3=0, tb4=0,
e0 = logit_eps_mn, e1 = 0, e2 = 0, e3=0, e4=0
)

# coeff variation
cvar = 2

par_lo = list(ab0 = logit_alphab_mn - abs(logit_alphab_mn*cvar),
ab1 = logit_alphab_mn/max(abs(data$ENV1)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV1))),
ab2 = logit_alphab_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV1^2))), 
ab3 = logit_alphab_mn/max(abs(data$ENV2)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV2))), 
ab4= logit_alphab_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV2^2))),
at0 = logit_alphat_mn - abs(logit_alphat_mn*cvar),
at1 = logit_alphat_mn/max(abs(data$ENV1)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV1))),
at2 = logit_alphat_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV1^2))),
at3 = logit_alphat_mn/max(abs(data$ENV2)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV2))),
at4 = logit_alphat_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV2^2))),
bb0 = logit_betab_mn - abs(logit_betab_mn*cvar),
bb1 = logit_betab_mn/max(abs(data$ENV1)) - abs(cvar*logit_betab_mn/max(abs(data$ENV1))),
bb2 = logit_betab_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_betab_mn/max(abs(data$ENV1^2))),
bb3 = logit_betab_mn/max(abs(data$ENV2)) - abs(cvar*logit_betab_mn/max(abs(data$ENV2))),
bb4 = logit_betab_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_betab_mn/max(abs(data$ENV2^2))),
bt0 = logit_betat_mn - abs(logit_betat_mn*cvar),
bt1 = logit_betat_mn/max(abs(data$ENV1)) - abs(cvar*logit_betat_mn/max(abs(data$ENV1))),
bt2 = logit_betat_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_betat_mn/max(abs(data$ENV1^2))),
bt3 =logit_betat_mn/max(abs(data$ENV2)) - abs(cvar*logit_betat_mn/max(abs(data$ENV2))),
bt4 = logit_betat_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_betat_mn/max(abs(data$ENV2^2))),
tt0 = logit_thetat_mn - abs(logit_thetat_mn*cvar),
tt1 = logit_thetat_mn/max(abs(data$ENV1)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV1))),
tt2 = logit_thetat_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV1^2))),
tt3 = logit_thetat_mn/max(abs(data$ENV2)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV2))),
tt4 = logit_thetat_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV2^2))),
tb0 = logit_thetab_mn - abs(logit_thetab_mn*cvar),
tb1 = logit_thetab_mn/max(abs(data$ENV1)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV1))),
tb2 = logit_thetab_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV1^2))),
tb3 = logit_thetab_mn/max(abs(data$ENV2)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV2))),
tb4 = logit_thetab_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV2^2))),
e0 = logit_eps_mn - abs(logit_eps_mn*cvar),
e1 = logit_eps_mn/max(abs(data$ENV1)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1))),
e2 = logit_eps_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1^2))),
e3 = logit_eps_mn/max(abs(data$ENV2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2))), 
e4 = logit_eps_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2^2))))

par_hi = list(ab0 = logit_alphab_mn + abs(logit_alphab_mn*cvar),
ab1 = logit_alphab_mn/max(abs(data$ENV1)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV1))),
ab2 = logit_alphab_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV1^2))), 
ab3 = logit_alphab_mn/max(abs(data$ENV2)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV2))), 
ab4= logit_alphab_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV2^2))),
at0 = logit_alphat_mn - abs(logit_alphat_mn*cvar),
at1 = logit_alphat_mn/max(abs(data$ENV1)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV1))),
at2 = logit_alphat_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV1^2))),
at3 = logit_alphat_mn/max(abs(data$ENV2)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV2))),
at4 = logit_alphat_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV2^2))),
bb0 = logit_betab_mn - abs(logit_betab_mn*cvar),
bb1 = logit_betab_mn/max(abs(data$ENV1)) + abs(cvar*logit_betab_mn/max(abs(data$ENV1))),
bb2 = logit_betab_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_betab_mn/max(abs(data$ENV1^2))),
bb3 = logit_betab_mn/max(abs(data$ENV2)) + abs(cvar*logit_betab_mn/max(abs(data$ENV2))),
bb4 = logit_betab_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_betab_mn/max(abs(data$ENV2^2))),
bt0 = logit_betat_mn - abs(logit_betat_mn*cvar),
bt1 = logit_betat_mn/max(abs(data$ENV1)) + abs(cvar*logit_betat_mn/max(abs(data$ENV1))),
bt2 = logit_betat_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_betat_mn/max(abs(data$ENV1^2))),
bt3 =logit_betat_mn/max(abs(data$ENV2)) + abs(cvar*logit_betat_mn/max(abs(data$ENV2))),
bt4 = logit_betat_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_betat_mn/max(abs(data$ENV2^2))),
tt0 = logit_thetat_mn - abs(logit_thetat_mn*cvar),
tt1 = logit_thetat_mn/max(abs(data$ENV1)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV1))),
tt2 = logit_thetat_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV1^2))),
tt3 = logit_thetat_mn/max(abs(data$ENV2)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV2))),
tt4 = logit_thetat_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV2^2))),
tb0 = logit_thetab_mn - abs(logit_thetab_mn*cvar),
tb1 = logit_thetab_mn/max(abs(data$ENV1)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV1))),
tb2 = logit_thetab_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV1^2))),
tb3 = logit_thetab_mn/max(abs(data$ENV2)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV2))),
tb4 = logit_thetab_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV2^2))),
e0 = logit_eps_mn - abs(logit_eps_mn*cvar),
e1 = logit_eps_mn/max(abs(data$ENV1)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1))),
e2 = logit_eps_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1^2))),
e3 = logit_eps_mn/max(abs(data$ENV2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2))), 
e4 = logit_eps_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2^2)))
)

par_initStep =lapply(names(par), function(x){
max(abs(as.numeric(par[x])-as.numeric(par_lo[x])), abs(as.numeric(par[x])-as.numeric(par_hi[x])))
})
names(par_initStep) = names(par)



#ENV1 = data$ENV1
#ENV2 =  data$ENV2
#save(ENV1, ENV2, file ="ENV.RData")

testBounds =
function(par, data)
{
    res = TRUE
    ENV1 = data$ENV1
    ENV2 = data$ENV2
    Hv = data$Hv
    Ha = data$Ha
    ET = data$ET
    EB = data$EB
    EM = data$EM
    

#    load("ENV.RData")
    logit_alphab 	= par$ab0 + par$ab1*ENV1 + par$ab2*ENV1^2 + par$ab3*ENV2 + par$ab4*ENV2^2
    logit_alphat 	= par$at0 + par$at1*ENV1 + par$at2*ENV1^2 + par$at3*ENV2 + par$at4*ENV2^2
    logit_betab 	= par$bb0 + par$bb1*ENV1 + par$bb2*ENV1^2 + par$bb3*ENV2 + par$bb4*ENV2^2
    logit_betat 	= par$bt0 + par$bt1*ENV1 + par$bt2*ENV1^2 + par$bt3*ENV2 + par$bt4*ENV2^2
    logit_thetab	= par$tb0 + par$tb1*ENV1 + par$tb2*ENV1^2 + par$tb3*ENV2 + par$tb4*ENV2^2
    logit_thetat	= par$tt0 + par$tt1*ENV1 + par$tt2*ENV1^2 + par$tt3*ENV2 + par$tt4*ENV2^2
    logit_eps 	= par$e0  + par$e1*ENV1  + par$e2*ENV1^2 + par$e3*ENV2 + par$e4*ENV2^2
    
    
    alphab = exp(logit_alphab)/(1+exp(logit_alphab))
    alphat = exp(logit_alphat)/(1+exp(logit_alphat))
    betab = exp(logit_betab)/(1+exp(logit_betab))
    betat = exp(logit_betat)/(1+exp(logit_betat))
    thetab = exp(logit_thetab)/(1+exp(logit_thetab))
    thetat = exp(logit_thetat)/(1+exp(logit_thetat))
    eps = exp(logit_eps)/(1+exp(logit_eps))

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

    if(sum(alphab_h==0)>0 | sum(alphab_h==1)>0 ) res=FALSE
    if(sum(alphat_h==0)>0 | sum(alphat_h==1)>0 ) res=FALSE
    if(sum(betab_h==0)>0 | sum(betab_h==1)>0 ) res=FALSE
    if(sum(betat_h==0)>0 | sum(betat_h==1)>0 ) res=FALSE
    if(sum(thetab_h==0)>0 | sum(thetab_h==1)>0 ) res=FALSE
    if(sum(thetat_h==0)>0 | sum(thetat_h==1)>0 ) res=FALSE
    if(sum(eps==0)>0 | sum(eps==1)>0 ) res=FALSE    
    
    
return(res)
}





