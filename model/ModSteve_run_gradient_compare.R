rm(list=ls())
# Range kutta 
#------------------------------------------
library(deSolve)
library(rootSolve)
#------------------------------------------

## parameters 
# 

#------------------------------------------
#dat = read.table("../fit_model/data/data_reshaped_RBTM_herbivores.txt")
## parameters vegetation model
# Av Temp = 3
# Av precip = 
#------------------------------------------
veget_pars = read.table("../fit_model_m3_ok/estimated_params/par_herbivores_m3.txt")
attach(veget_pars)

source("params_Herbivores.r")
source("model_fct.r")

#
#tpresent = 0
#tfutur = 3
#ppresent = .9
#pfutur = .9
#

tpresent = 3
tfutur = 6
ppresent = 1
pfutur = 1

tpresent = 2
tfutur = 4
ppresent = .8
pfutur = .8



#---------------
ENV1 = tpresent
ENV2 = ppresent

 
 logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV1^2 + ab3*ENV2 + ab4*ENV2^2 + ab5*ENV1^3 + ab6*ENV2^3
    logit_alphat 	= at0 + at1*ENV1 + at2*ENV1^2 + at3*ENV2 + at4*ENV2^2 + at5*ENV1^3 + at6*ENV2^3
    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV1^2 + bb3*ENV2 + bb4*ENV2^2 + bb5*ENV1^3 + bb6*ENV2^3
    logit_betat 	= bt0 + bt1*ENV1 + bt2*ENV1^2 + bt3*ENV2 + bt4*ENV2^2 + bt5*ENV1^3 + bt6*ENV2^3
    logit_thetab	= tb0 + tb1*ENV1 + tb2*ENV1^2 + tb3*ENV2 + tb4*ENV2^2 + tb5*ENV1^3 + tb6*ENV2^3
    logit_thetat	= tt0 + tt1*ENV1 + tt2*ENV1^2 + tt3*ENV2 + tt4*ENV2^2 + tt5*ENV1^3 + tt6*ENV2^3
    logit_eps 	= e0  + e1*ENV1  + e2*ENV1^2 + e3*ENV2 + e4*ENV2^2 + e5*ENV1^3 + e6*ENV2^3

alphab = exp(logit_alphab)/(1+exp(logit_alphab))
alphat = exp(logit_alphat)/(1+exp(logit_alphat))
betab = exp(logit_betab)/(1+exp(logit_betab))
betat = exp(logit_betat)/(1+exp(logit_betat))
thetab = exp(logit_thetab)/(1+exp(logit_thetab))
thetat = exp(logit_thetat)/(1+exp(logit_thetat))
eps = exp(logit_eps)/(1+exp(logit_eps))


# parameter list
#------------------------------------------

pars = list(alphat=alphat,
alphab=alphab,
betat = betat,
betab = betab,
thetat = thetat, 
thetab = thetab,
eps = eps, 
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

#------------------------------------------
# run simulation

#--- WITH HERBIVORES
T0 = c(TT=0.08, B=0.34, M=.32, Ha = dva, Hv=dvv)

out = stode(func=model, y=T0, parms = pars, positive=TRUE, maxiter=5000)#, 
out$y

#--- WITHOUT HERBIVORES
T0 = c(TT=0.08, B=0.34, M=.32, Ha = 0, Hv=0)

out_noH = stode(func=model, y=T0, parms = pars, positive=TRUE, maxiter=5000)#, 
out_noH$y


#-------FUTURE-----------------------------------

ENV1 = tfutur
ENV2 = pfutur

 logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV1^2 + ab3*ENV2 + ab4*ENV2^2 + ab5*ENV1^3 + ab6*ENV2^3
    logit_alphat 	= at0 + at1*ENV1 + at2*ENV1^2 + at3*ENV2 + at4*ENV2^2 + at5*ENV1^3 + at6*ENV2^3
    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV1^2 + bb3*ENV2 + bb4*ENV2^2 + bb5*ENV1^3 + bb6*ENV2^3
    logit_betat 	= bt0 + bt1*ENV1 + bt2*ENV1^2 + bt3*ENV2 + bt4*ENV2^2 + bt5*ENV1^3 + bt6*ENV2^3
    logit_thetab	= tb0 + tb1*ENV1 + tb2*ENV1^2 + tb3*ENV2 + tb4*ENV2^2 + tb5*ENV1^3 + tb6*ENV2^3
    logit_thetat	= tt0 + tt1*ENV1 + tt2*ENV1^2 + tt3*ENV2 + tt4*ENV2^2 + tt5*ENV1^3 + tt6*ENV2^3
    logit_eps 	= e0  + e1*ENV1  + e2*ENV1^2 + e3*ENV2 + e4*ENV2^2 + e5*ENV1^3 + e6*ENV2^3

alphab = exp(logit_alphab)/(1+exp(logit_alphab))
alphat = exp(logit_alphat)/(1+exp(logit_alphat))
betab = exp(logit_betab)/(1+exp(logit_betab))
betat = exp(logit_betat)/(1+exp(logit_betat))
thetab = exp(logit_thetab)/(1+exp(logit_thetab))
thetat = exp(logit_thetat)/(1+exp(logit_thetat))
eps = exp(logit_eps)/(1+exp(logit_eps))


pars = list(alphat=alphat,
alphab=alphab,
betat = betat,
betab = betab,
thetat = thetat, 
thetab = thetab,
eps = eps, 
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


#--- WITHOUT HERBIVORES
evol = ode(func=model, y=out$y, parms = pars, times = 1:500)#, method="iteration")


#--- WITHOUT HERBIVORES

evol_noH = ode(func=model, y=out_noH$y, parms = pars, times = 1:500)#, 


##-- GRAPHS ESA

colo = c("orange", "lightgreen", "darkgreen", "violet", "brown", "blue")

par(mfrow = c(1,2), lwd = 2)
plot(1:500, evol[,"B"], type = "l", ylim = c(0,1), xlab = "time", ylab = "Boreal proportion",col = colo[3])
lines(1:500, evol_noH[,"B"], lty = 2, col =colo[3])
arrows(x0 = which(evol_noH[, "B"]<0.41 & evol_noH[, "B"]>0.39),  y0 =0.4, x1 = which(evol[, "B"]<0.41 & evol[, "B"]>0.39) , y1 = 0.4, col = colo[6])
legend("topright" , legend = "without herbivores", lty = 2, bty = "n", col = colo[3])

plot(1:500, evol[,"TT"], type = "l", ylim = c(0,1), xlab = "time", ylab = "Temperate proportion",col = colo[2])
lines(1:500, evol_noH[,"TT"], lty = 2, col = colo[2])
arrows(x0 = which(evol_noH[, "B"]<0.21 & evol_noH[, "B"]>0.19),  y0 =0.2, x1 = which(evol[, "B"]<0.21 & evol[, "B"]>0.19) , y1 = 0.2, col = colo[6])

legend("topright" , legend = "without herbivores", lty = 2, bty = "n", col = colo[2])

#
#plot(1:500, evol[,"M"], type = "l", ylim = c(0,1), xlab = "time", ylab = "Mixed proportion",col = colo[4])
#lines(1:500, evol_noH[,"M"], lty = 2)
#
#plot(1:500, 1-evol[,"B"]-evol[,"TT"]-evol[,"M"], type = "l", ylim = c(0,1), xlab = "time", ylab = "Post-dist proportion",col = colo[1])
#lines(1:500, 1-evol_noH[,"B"]-evol_noH[,"TT"]-evol_noH[,"M"], lty = 2)


#
dev.copy2pdf(file = "../graphs/climate_change.pdf")



