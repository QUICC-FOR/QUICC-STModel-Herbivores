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

ENV1 = 2
ENV2 = 1000/1000

ENV = expand.grid(TP = seq(-1,5,0.2), PP = seq(0.75,1.5,0.1))
TP = ENV$TP
PP = ENV$PP

res_env = mclapply(1:nrow(ENV), function(i)
{

ENV1 = TP[i]
ENV2 = PP[i]

 
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
source("model_fct.r")


# initial conditions

T0 = c(TT=0.08, B=0.34, M=.32, Ha = dva, Hv=dvv)

#------------------------------------------
# run simulation
out = stode(func=model, y=T0, parms = pars, positive=TRUE, maxiter=5000)#, method="iteration")

#res = c(out$y,R = 1-out$y["TT"]-out$y["B"]-out$y["M"])
#res = res[c(6, 1:5)]
#names(res) = c("R", "T", "B", "M", "moose", "deer")

return(out$y)

})

tab_env = t(data.frame(res_env))
rownames(tab_env) = 1:nrow(tab_env)
head(tab_env)
tab_env = data.frame(tab_env)
tab_env$R = 1 - tab_env$TT - tab_env$B - tab_env$M
    
sum(tab_env<0)
tab_env[tab_env<0]=NA

selec = which(!is.na(apply(tab_env, 1, sum)))

par(mfrow = c(2,3))
plot(lowess(tab_env[selec,"B"]~ENV$TP[selec]), type = "l", ylim = c(0,1), xlab = "T°", ylab = "Boreal proportion")
plot(lowess(tab_env[selec,"TT"]~ENV$TP[selec]), type = "l", ylim = c(0,1), xlab = "T°", ylab = "Temperate forest proportion")
plot(lowess(tab_env[selec,"M"]~ENV$TP[selec]), type = "l", ylim = c(0,1), xlab = "T°", ylab = "Mixed forest proportion")
plot(lowess(tab_env[selec,"R"]~ENV$TP[selec]), type = "l", ylim = c(0,1), xlab = "T°", ylab = "Young forest proportion")
plot(lowess(tab_env[selec,"Hv"]~ENV$TP[selec]), type = "l", xlab = "T°", ylab = "Deer biomass")
plot(lowess(tab_env[selec,"Ha"]~ENV$TP[selec]), type = "l", xlab = "T°", ylab = "Moose biomass")

tabAll = cbind(tab_env, ENV)
head(tabAll)

colo = c("orange", "lightgreen", "darkgreen", "violet", "brown", "blue")
#barplot(res, col = colo, ylim = c(0,1))




