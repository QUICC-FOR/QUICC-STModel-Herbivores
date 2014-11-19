rm(list=ls())
# Range kutta 
#------------------------------------------
library(deSolve)
library(rootSolve)
library(parallel)
#------------------------------------------

## parameters 
# 

#------------------------------------------
#dat = read.table("../fit_model/data/data_reshaped_RBTM_herbivores.txt")
## parameters vegetation model
# Av Temp = 3
# Av precip = 
#------------------------------------------
veget_pars = read.table("../fit_model/estimated_params/par_herbivores_m3_ok.txt")
attach(veget_pars)

source("params_Herbivores.r")

TP = seq(0,5, length.out = 20)
PP = seq(0.8, 1.4, length.out = 20)

ENV = expand.grid(TP, PP)

#--- WITH HERBIVORES

res_env = mclapply(1:nrow(ENV), function(i)
{

ENV1 = ENV[i,1]
ENV2 = ENV[i,2]

 
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


#--- WITHOUT HERBIVORES


res_env_noH = mclapply(1:nrow(ENV), function(i)
{

ENV1 = ENV[i,1]
ENV2 = ENV[i,2]

 
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

T0 = c(TT=0.08, B=0.34, M=.32, Ha = 0, Hv=0)

#------------------------------------------
# run simulation
out = stode(func=model, y=T0, parms = pars, positive=TRUE, maxiter=5000)#, method="iteration")

#res = c(out$y,R = 1-out$y["TT"]-out$y["B"]-out$y["M"])
#res = res[c(6, 1:5)]
#names(res) = c("R", "T", "B", "M", "moose", "deer")

return(out$y)

})

tab_env_noH = t(data.frame(res_env_noH))
rownames(tab_env_noH) = 1:nrow(tab_env_noH)
head(tab_env_noH)
tab_env_noH = data.frame(tab_env_noH)
tab_env_noH$R = 1 - tab_env_noH$TT - tab_env_noH$B - tab_env_noH$M
    
sum(tab_env_noH<0)
tab_env[tab_env_noH<0]=NA

selec_noH = which(!is.na(apply(tab_env_noH, 1, sum)))

##------
par(mfrow = c(2,1))
plot(lowess(cbind(ENV$TP,tab_env$TT)[selec,]), type = "l", xlab = "Temperature", ylab = "Vegetation proportions", col = "lightgreen")
lines(lowess(cbind(ENV$TP,tab_env$B)[selec,]), col = "darkgreen")
lines(lowess(cbind(ENV$TP,tab_env$M)[selec,]), col = "blue")
lines(lowess(cbind(ENV$TP,tab_env$R)[selec,], f = 1/2), col = "orange")

#plot(lowess(cbind(ENV$TP,tab_env_noH$TT)[selec_noH,]), type = "l", xlab = "Temperature", ylab = "Vegetation proportions", col = "lightgreen")
lines(lowess(cbind(ENV$TP,tab_env_noH$TT)[selec_noH,]), col = "lightgreen", lty=2)
lines(lowess(cbind(ENV$TP,tab_env_noH$B)[selec_noH,]), col = "darkgreen", lty = 2)
lines(lowess(cbind(ENV$TP,tab_env_noH$M)[selec_noH,]), col = "blue", lty = 2)
lines(lowess(cbind(ENV$TP,tab_env_noH$R)[selec_noH,], f = 1/2), col = "orange", lty =2)

par(mfrow = c(2,1))
plot((cbind(ENV$TP,tab_env$TT)[selec,]), type = "l", xlab = "Temperature", ylab = "Vegetation proportions", col = "lightgreen")
lines((cbind(ENV$TP,tab_env$B)[selec,]), col = "darkgreen")
lines((cbind(ENV$TP,tab_env$M)[selec,]), col = "blue")
lines((cbind(ENV$TP,tab_env$R)[selec,]), col = "orange")

plot((cbind(ENV$TP,tab_env_noH$TT)[selec_noH,]), type = "l", xlab = "Temperature", ylab = "Vegetation proportions", col = "lightgreen")
lines((cbind(ENV$TP,tab_env_noH$TT)[selec_noH,]), col = "lightgreen", lty=2)
lines((cbind(ENV$TP,tab_env_noH$B)[selec_noH,]), col = "darkgreen", lty = 2)
lines((cbind(ENV$TP,tab_env_noH$M)[selec_noH,]), col = "blue", lty = 2)
lines((cbind(ENV$TP,tab_env_noH$R)[selec_noH,]), col = "orange", lty =2)

##-----
coexistence <- function(x)
{
res = 1
if((x[1]>0.0001) & (x[2]<0.0001)) res = 1
if((x[1]<0.0001) & (x[2]>0.0001)) res = 2
if((x[1]>0.0001) & (x[2]>0.0001)) res = 3
return(res)
}
land.with = apply(tab_env[,c("TT", "B", "M")], 1, coexistence)
land.without = apply(tab_env_noH[,c("TT", "B", "M", "R")], 1, which.max)


table(land.with)
landscape = land.with
##-- GRAPHS 2d
Z = matrix(landscape,nr = length(TP), nc = length(PP))
quartz(width = 6, height = 6)
colo = c("lightgreen", "darkgreen", "orange", "red")
layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#title(title,cex=2)
#legend("center",legend = c("AltSS","Boreal Wins","Temperate Wins","Coexistence"),fill = colo,bty = "n",horiz = TRUE,cex = 0.8)
par(mar=c(5,5,0,2))
image(TP,PP*1000,Z,xlab = "Mean annual temperature", ylab = "Annual precipitation (mm)", cex.lab = 1.5, cex.axis = 1.25, col = colo, breaks = c(0:4))#grey(c(0:3)/3))

dev.copy2pdf(file = "../graphs/Coexistence_area_herbivores.pdf")

##
land.with = apply(tab_env[which(ENV[,2]>1 & ENV[,2]<1.05),c("TT", "B", "M")], 1, coexistence)

TpGrad = ENV[which(ENV[,2]>1 & ENV[,2]<1.05),1]
herbi = tab_env[which(ENV[,2]>1 & ENV[,2]<1.05),c("TT", "B", "M")]

plot(herbi$TT~TpGrad, type = "l", col = "lightgreen", lwd = 2, ylim = c(0,1))
lines(herbi$B~TpGrad, col = "darkgreen", lwd =2)


##-- GRAPHS ESA

colo = c("orange", "lightgreen", "darkgreen", "violet", "brown", "blue")

par(mfrow = c(2,2), lwd = 2)
plot(tab_env[selec,"B"]~ENV$TP[selec], type = "l", ylim = c(0,1), xlab = "T°", ylab = "Boreal proportion",col = colo[3])
axis(3, at = seq(0, 5, length.out = 6), seq(800, 1400, length.out = 6))
mtext("precipitation", 3, line = 3, cex = .8)
lines(tab_env_noH[selec_noH,"B"]~ENV$TP[selec_noH], lty = 2)

plot(tab_env[selec,"TT"]~ENV$TP[selec], type = "l", ylim = c(0,1), xlab = "T°", ylab = "Temperate proportion", col = colo[2])
axis(3, at = seq(0, 5, length.out = 6), seq(800, 1400, length.out = 6))
mtext("precipitation", 3, line = 3, cex = .8)
lines(tab_env_noH[selec_noH,"TT"]~ENV$TP[selec_noH], lty = 2)


plot(tab_env[selec,"M"]~ENV$TP[selec], type = "l", ylim = c(0,1), xlab = "T°", ylab = "Mixed proportion", col = colo[4])
axis(3, at = seq(0, 5, length.out = 6), seq(800, 1400, length.out = 6))
mtext("precipitation", 3, line = 3, cex = .8)
lines(tab_env_noH[selec_noH,"M"]~ENV$TP[selec_noH], lty = 2)


plot(tab_env[selec,"R"]~ENV$TP[selec], type = "l", ylim = c(0,1), xlab = "T°", ylab = "Post-disturb. proportion", col = colo[1])
axis(3, at = seq(0, 5, length.out = 6), seq(800, 1400, length.out = 6))
mtext("precipitation", 3, line = 3, cex = .8)
lines(tab_env_noH[selec_noH,"R"]~ENV$TP[selec_noH], lty = 2)


#plot(tab_env[selec,"Hv"]~ENV$TP[selec], type = "l", xlab = "T°", ylab = "Deer biomass")
#plot(tab_env[selec,"Ha"]~ENV$TP[selec], type = "l", xlab = "T°", ylab = "Moose biomass")
#





