rm(list=ls())
#------------------------------------------
# coexistence
#------------------------------------------
library(pracma) #lambertWp
library(rootSolve)


# both
#----------



k_fct = function(Ha, Hv)
{
k = exp(-log(1/k0)*Hv/Ha)
return(k)
}



#--------------------------------------------------------------------

obs = read.table(file="../fit_model/data/data_allyears_RBTM.txt")
obsMoose = read.table("../fit_model/data/densities_moose.txt", h=T, sep="\t")
obsMooseMerge = merge(obs, obsMoose, by = "id_plot")

mooseList = split(obsMooseMerge$class_final, obsMooseMerge$NO_ZONE)

obs_Moose = unlist(lapply(split(obsMooseMerge$densites_10ha, obsMooseMerge$NO_ZONE), mean))


rtbmList = lapply(mooseList, function(x){c(sum(x=="R")/length(x), sum(x=="T")/length(x), sum(x=="B")/length(x), sum(x=="M")/length(x))})



#--------------------------------------------------------------------

#--------
# moose
size = 5000

za = 1
ca = 0.79
ma0 = 0.2
mas = 0.25
taua = 13
mua = (taua/0.07)*(size/100)/sqrt(314)
rhoa = 6
nua = (rhoa/0.15)*(size/100)/sqrt(314)
phia = 0.5
pa = 1.7
ra = 5


k0=0.5; 
u = 0.5*1000*size ; v = 200*size ;  w= 100*size ; x= 200*size

Ia_fct = function(Ha, Hv, R, T, B, M)
{

    k = k_fct(Ha, Hv)
    
if(Ha>0)
    {
    F = k*u*R/Ha
    G = k*(v*T + w*B + x*M)/Ha
    Ia1 = taua*F/(mua + F)
    Ia2 = (rhoa*G/(nua + G) )* (phia + (1-phia)/( 1 + exp(ra*( taua*F/(mua + F) -pa) ) ) )
    }else
    {
    Ia1 = Ia2 = 0
    print("Ha is null")
    }

return(c(Ia1=Ia1, Ia2=Ia2))
}

# calc dHadt
Ha_eq_fct = function(Ha, Hv, R, T, B, M)
{
    ma = mas/(1+(mas/ma0-1)*(B+M))
    Ia = Ia_fct(Ha, Hv, R, T, B, M)
    Ia_tot = Ia["Ia1"]+Ia["Ia2"]
    dHadt = Ha * ( ca*Ia_tot - (pa + ma + exp(-za*ca*Ia_tot) ))
    if(dHadt<=(-Ha)) dHadt = -Ha
    return(dHadt)
}


# find Ha for which dHadt = 0 given RTBM and Hv
Ha_eq_solve = function(RTBMH)
{

#RTBMH = c(.5, .2, 0.2, 0.1,100); # debug
R = as.numeric(RTBMH[1])
T = as.numeric(RTBMH[2])
B = as.numeric(RTBMH[3])
M = as.numeric(RTBMH[4])

#Hv = as.numeric(RTBMH[5])
Hv=0

# 1 # find Ha for which dHadt = 0
# cherche solution positive sinon 0 est une solution quelque-soit Hv
hamin = as.numeric(Ha_eq_fct(0.0001, Hv=Hv, R=R, T=T, B=B, M=M))
hamax = as.numeric(Ha_eq_fct(1000000, Hv=Hv, R=R, T=T, B=B, M=M))
if(hamin*hamax<0) # if it is possible to cross the zero line
{
Ha_eq = uniroot(Ha_eq_fct, Hv=Hv, R=R, T=T, B=B, M=M, interval=c(0.0001, 1000000), f.lower = hamin, f.upper = hamax, tol = 0.001)$root
} else Ha_eq = sign(hamin);

return(Ha_eq)
}


#--------------------------------------------------------------------

Ha_eq_obs = lapply(rtbmList,  Ha_eq_solve) 
pred_K_Moose = unlist(Ha_eq_obs)/5/350


##
#--------------------------------------------------------------------

#--------------------------------------------------------------------


obsDeer = read.table("../fit_model/data/densities_deer.txt", h=T, sep="\t")
obsDeerMerge = merge(obs, obsDeer, by = "id_plot")

deerList = split(obsDeerMerge$class_final, obsDeerMerge$NO_ZONE)

obs_Deer= unlist(lapply(split(obsDeerMerge$densites_10ha, obsDeerMerge$NO_ZONE), mean))


rtbmList = lapply(deerList, function(x){c(sum(x=="R")/length(x), sum(x=="T")/length(x), sum(x=="B")/length(x), sum(x=="M")/length(x))})


#--------------------------------------------------------------------

#--------
#deer
size = 1000 # ha

zv = 1
cv = 0.5
mv0 = 0.2
mvs = 0.5
tauv = 17.5
muv = (tauv/0.07)*(size/100)/sqrt(518)
rhov = 4.5
nuv = (rhov/0.15)*(size/100)/sqrt(518)
pv = 1.87
phiv = 0.57
rv = 5


k0=0.5; 
u = 0.5*1000*size ; v = 200*size ;  w= 100*size ; x= 200*size


Iv_fct = function(Ha, Hv, R, T, B, M)
{

    k = k_fct(Ha, Hv)
   
if(Hv>0)
{
    F = u*(1-k)*R/Hv
    G = (1-k) * (v*T + w*B + x*M)/Hv
    Iv1 = tauv*F/(muv + F)
    Iv2 = (rhov*G/(nuv + G) )* (phiv + (1-phiv)/( 1 + exp(rv*( tauv*F/(muv + F) -pv) ) ) )
    }else
{
Iv1 = Iv2 = 0
print("Hv is null")
}

    return(c(Iv1=Iv1, Iv2=Iv2))
 
}


# calc dHvdt
Hv_eq_fct = function(Hv, Ha, R, T, B, M)
{
    mv = mvs/(1+(mvs/mv0-1)*(B+M))
    Iv = Iv_fct(Ha, Hv, R, T, B, M)
    Iv_tot = Iv["Iv1"]+Iv["Iv2"]
    dHvdt = Hv * (1-mv)*( cv*Iv_tot - (pv + exp(-zv*cv*Iv_tot)) -mv)
    if(dHvdt<=(-Hv)) dHvdt = -Hv
    return(dHvdt)
}


# find Hv for which dHvdt = 0 given RTBM and Ha
Hv_eq_solve = function(RTBMH)
{

#RTBMH = c(.5, .2, 0.2, 0.1,100); # debug
R = as.numeric(RTBMH[1])
T = as.numeric(RTBMH[2])
B = as.numeric(RTBMH[3])
M = as.numeric(RTBMH[4])

#Ha = as.numeric(RTBMH[5])
Ha = 0

# 1 # find Ha for which dHadt = 0
# cherche solution positive sinon 0 est une solution quelque-soit Hv
hvmin = as.numeric(Hv_eq_fct(0.0001, Ha=Ha, R=R, T=T, B=B, M=M))
hvmax = as.numeric(Hv_eq_fct(1000000, Ha=Ha, R=R, T=T, B=B, M=M))
if(hvmin*hvmax<0) # if it is possible to cross the zero line
{
Hv_eq = uniroot(Hv_eq_fct, Ha=Ha, R=R, T=T, B=B, M=M, interval=c(0.0001, 1000000), f.lower = hvmin, f.upper = hvmax, tol = 0.001)$root
} else Hv_eq = sign(hvmin);

return(Hv_eq)
}

#--------------------------------------------------------------------

Hv_eq_obs = lapply(rtbmList,  Hv_eq_solve) 
pred_K_Deer = unlist(Hv_eq_obs)/10/75



#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
quartz(width = 8, height = 4)

par(mfrow = c(1,2))
plot(y=pred_K_Moose, x=obs_Moose, xlim = c(0,10), ylim = c(0,10), pch=20, xlab = "observed densities (10km2)", ylab = "predicted carrying capacity", main = "Moose")
abline(0,1,lty=2)

plot(y=pred_K_Deer, x=obs_Deer, xlim = c(0,10), ylim = c(0,10), pch=20, xlab = "observed densities (km2)", ylab = "predicted carrying capacity", main = "White-tailed deer")
#text(y=pred_K_Deer, x=obs_Deer,  labels = names(Hv_eq_obs), cex = .5)
abline(0,1,lty=2)

dev.copy2pdf(file = "../graphs/validationK.pdf")
