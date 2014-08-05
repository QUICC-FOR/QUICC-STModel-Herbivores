rm(list=ls())
# Range kutta 
#------------------------------------------
library(deSolve)

#------------------------------------------

## parameters 
# 

#------------------------------------------

## parameters vegetation model
# Av Temp = 3
#------------------------------------------

ab0 = 1.2
ab1 = -1
ab2 = 4.2
at0 = .8
at1 = 0
at2 = .08
bb0 = -5.23
bb1 = 0
bb2 = 0
bt0 = -.6
bt1 = 0
bt2 = 0
tt0 = -.4
tt1 = 0
tt2 = 0
tb0 = -3.8
tb1 = 0
tb2 = 0
e0 = -.4
e1 = 0
e2 = 0

ENV = 5


logit_alphab 	= ab0 + ab1*ENV + ab2*ENV^2
logit_alphat 	= at0 + at1*ENV + at2*ENV^2
logit_betab 	= bb0 + bb1*ENV + bb2*ENV^2
logit_betat 	= bt0 + bt1*ENV + bt2*ENV^2
logit_thetab	= tb0 + tb1*ENV + tb2*ENV^2
logit_thetat	= tt0 + tt1*ENV + tt2*ENV^2
logit_eps 	= e0  + e1*ENV  + e2*ENV^2

alphab = exp(logit_alphab)/(1+exp(logit_alphab))
alphat = exp(logit_alphat)/(1+exp(logit_alphat))
betab = exp(logit_betab)/(1+exp(logit_betab))#*(EC+EM)
betat = exp(logit_betat)/(1+exp(logit_betat))#*(ED+EM)
thetab = exp(logit_thetab)/(1+exp(logit_thetab))
thetat = exp(logit_thetat)/(1+exp(logit_thetat))
eps = exp(logit_eps)/(1+exp(logit_eps))
epsilont = eps
epsilonb = eps
epsilonm = eps


# moose
za = 1
ca = 0.79
ma = 0.25
taua = 13
mua = 3*taua
rhoa = 6
nua = rhoa
phia = 0.5
pa = 1.7
ra = 5

#deer
zv = 1
cv = 0.5
mv0 = 0.2
mvs = 0.5
tauv = 17.5
muv = tauv
rhov = 4.5
nuv = tauv
pv = 1.87
phiv = 0.57
rv = 5


#interactions

size = 1000 # ha

k0=0.5; 
u = 1000*size ; v = 200*size ;  w= 200*size ; x= 200*size
#(dry matter content)



omegatHa = 0.3
omegabHa = 0.35

omegatHv = 0.325
omegabHv = 0.35

kappata = 0.25
kappaba = 0.07

kappatv = 0.265
kappabv = 0.175

#------------------------------------------

##  coupled model

#------------------------------------------

n=100

# initial conditions
T0 = c(TT=0.25, B=0.25, M=.25, Ha = 350, Hv=75)

# simulation conditions
times = seq(0, n, by=1)

# parameters
#------------------------------------------

pars = c(alphat=alphat,
alphab=alphab,
betat = betat,
betab = betab,
thetat = thetat, 
thetab = thetab,
epsilont = epsilont,
epsilonb = epsilonb, 
epsilonm = epsilonm,
phia = phia, 
phiv = phiv, 
taua =taua, 
tauv = tauv,
rhoa = rhoa, 
rhov = rhov,
mua = mua,
muv = muv,
nua=nua, 
nuv = nuv,
ca = ca, 
cv = cv,
ra = ra, 
rv = rv,
pa = pa, 
pv = pv, 
mv0 = mv0,
mvs = mvs,
ma = ma,  
za = za, 
zv = zv,
u = u, 
v  = v,
w = w, 
x = x, 
k0 = k0, 
omegatHa = omegatHa, 
omegatHv = omegatHv, 
omegabHa = omegabHa, 
omegabHv = omegabHv, 
kappata = kappata, 
kappaba = kappaba,
kappatv = kappatv, 
kappabv = kappabv
)




# model
#------------------------------------------
model = function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
	
	    #test
	       
#        TT = 0.2856791 
#        B=0.2057428 
#        M=0.3304669 
#        Ha=1590.129 
#        Hv=123.212

        cat(TT, B, M, Ha, Hv, "\n")

		R=1-TT-B-M
 #       cat("\n R; TT; B; M\n")
#        cat(R, TT, B, M)
		   
		# herbivores competition
		k = exp(-log(1/k0)*Hv/Ha)
		
		# moose and B
		mv_B = mvs/(1+(mvs/mv0-1)*(M+B))
		
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
        Ia2 = (rhoa*Ga/(nua + Ga) )* (phia + (1-phia)/( 1 + exp(ra*( taua*Fa/(mua + Fa) -pa) ) ) )
        Iv1 = tauv*Fv/(muv+Fv)
        Iv2 = (rhov*Gv/(nuv + Gv) )* (phiv + (1-phiv)/( 1 + exp(rv*( tauv*Fv/(muv + Fv) -pv) ) ) )        
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
        
       
        # compute new alpha, theta and beta
        alphab_h = alphab*(1-PRB)
        alphat_h = alphat*(1-PRT)
        
        thetab_h = thetab*(1-PMB)
        thetat_h = thetat*(1-PMT)
        
        betab_h = (1/2)*(1 + betab * cos(pi*PTB) + (1-betab)*cos(pi*(1+PTT)) ) 
        betat_h = (1/2)*(1 + betat * cos(pi*PBT) + (1-betat)*cos(pi*(1+PBB)) )
        
	    # vegetation dynamic
	    R_T = alphat_h*(M+TT)*(1 - alphab_h*(M+B))
	    R_B = alphab_h*(M+B)*(1 - alphat_h*(M+TT))
	    R_M = alphat_h*(M+TT)*alphab_h*(M+B)
	    T_M = betab_h*(B+M)
	    B_M = betat_h*(TT+M)
        dTT = R*R_T + M*thetat_h - TT*(epsilont + T_M)
        dB = R*R_B + M*thetab_h - B*(epsilonb + B_M)
        dM = R*R_M + TT*T_M + B*B_M - M*(epsilonm + thetat_h + thetab_h)
                
        # herbivore dynamic
        dHa = Ha*( (1-ma)*((Ia1+Ia2)*ca-(pa+exp(-za*(Ia1 + Ia2)*ca))) - ma)
        dHv = Hv*( (1-mv_B) *((Iv1+Iv2)*cv-(pv+exp(-zv*(Iv1 + Iv2)*cv))) - mv_B)
        
        # negative values 
        if(dHa<=(-Ha)){ dHa = round(-Ha)}
        if(dHv<=(-Hv)) dHv = -Hv
               

                      
		list(c(dTT, dB, dM, dHa, dHv))
        
#        Ha_new = Ha+dHa
#        Hv_new = Hv+dHv
#        TT_new = TT+dTT
#        B_new = B + dB
#        M_new = M + dM
        
#        cat("\n dB; B; R; M; R_B; B_M\n")
#        cat(dB, B, R, M, R_B, B_M)

#        if(Ha_new<0){Ha_new=0}
#        if(Hv_new<0){Hv_new=0}
#        if(TT_new<0){TT_new=0}
#        if(B_new<0){B_new=0}
#        if(M_new<0){M_new=0}
        
#        list(c(TT_new, B_new, M_new, Ha_new, Hv_new))
	})
}

#------------------------------------------
# run simulation
out = ode(func=model, y=T0, parms = pars, times =times)#, method="iteration")
# voir run steady
summary(out)

times = times[1:nrow(out)]

R_serie = 1 - out[,"TT"] - out[, "B"] - out[, "M"]
T_serie = out[,"TT"]
B_serie = out[, "B"]
M_serie = out[, "M"]
Ha_serie = out[, "Ha"]
Hv_serie = out[, "Hv"]


colo = c("orange", "lightgreen", "darkgreen", "violet", "brown", "blue")


par(mfrow = c(1,2))
plot(times, R_serie, type="n", ylim = c(0,1), ylab="vegetation proportion")
polygon(c(times, rev(times)), c(R_serie, rep(0, length(times)) ), col=colo[1], border=NA)
polygon(c(times, rev(times)), c(R_serie+T_serie, rev(R_serie)), col=colo[2], border=NA)
polygon(c(times, rev(times)), c(R_serie+T_serie+B_serie, rev(R_serie+T_serie)), col=colo[3], border=NA)
polygon(c(times, rev(times)), c(R_serie+T_serie+B_serie+M_serie, rev(R_serie+T_serie+B_serie)), col=colo[4], border=NA)

plot(times, Ha_serie, type="l", ylim = c(0,max(out[,c("Ha", "Hv")])), col=colo[5], lwd=2, ylab="herbivores abundances")
lines(times, Hv_serie, col=colo[6], lwd=2)
legend("topleft", legend = c("moose", "deer"), col = colo[5:6], lwd=2, bty="n")

