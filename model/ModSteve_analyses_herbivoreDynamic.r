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


# params 
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

Ha = as.numeric(RTBMH[5])

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



# params 
#--------
# moose
size = 20000

za = 1
ca = 0.79
ma0 = 0.25
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

Hv = as.numeric(RTBMH[5])

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




# RTBM space
#-------------

xx= seq(0,1,0.1)
xx2 = seq(0,0,0)

RTBH = expand.grid(list(R=xx, T=xx, B=xx, H=xx2), KEEP.OUT.ATTRS = FALSE)

RTBH= cbind(RTBH, sum= apply(RTBH[,1:3], 1, sum))

RTBH = RTBH[which(RTBH[,5]<=1),]

RTBMH = RTBH[,1:4]
RTBMH$M = 1- apply(RTBMH[,1:3], 1, sum)

RTBMH = RTBMH[, c(1:3, 5, 4)]
RTBMH = rbind(RTBMH, c(.26, .08, .34, .32, 0))
#testIndex = sample(1:nrow(RTBMH), 200)
#
#RTBMH[testIndex,]


# find equilibrium functions
#-------------


# -   equilibrium RESULTS -- #
#----------------------------



#Hstars = cbind(RTBMH[,], Ha_eq, Hv_eq)
#
#
#tab = Hstars[order(Hstars$B),]
#tab[tab==-1]=NA
#
#meanHa = unlist(lapply(split(tab$Ha_eq, tab$B), max, na.rm=T))
#meanHv = unlist(lapply(split(tab$Hv_eq, tab$B), max, na.rm=T))
#Bax = as.numeric(as.character(names(meanHa)))


#pdf("../graphs/herbivore_eq_biomasses.pdf")
#
#plot(Bax , meanHa/344, main="Herbivores biomasses at equilibrium", xlab = "Boreal proportion", ylab = "#nb indiv", type ="b", ylim = c(0, max(meanHa/344,meanHv/75, na.rm=T)))
#
#points(Bax, meanHv/75, col=2, type = "b", pch = 4)
#
#legend("topright", legend = c("Moose", "Deer"), col = c(1, 2), lwd = 1)

#dev.off()

Ha_eq = apply(RTBMH[,], 1,  Ha_eq_solve) 


par(mfrow = c(1, 2))
plot(Ha_eq/20/350~RTBMH[,"B"], main="Moose biomasses at equilibrium", xlab = "Boreal proportion", ylab = "#nb indiv (350kg) per 10km2", pch=4, cex=.5)
points(Ha_eq[nrow(RTBMH)]/20/350~RTBMH[nrow(RTBMH),"B"], col="red")
plot(Ha_eq/20/350~RTBMH[,"R"], main="Moose biomasses at equilibrium", xlab = "Immature forest proportion", ylab = "#nb indiv (350kg) per 10km2", pch=4, cex=.5)
points(Ha_eq[nrow(RTBMH)]/20/350~RTBMH[nrow(RTBMH),"R"], col="red")

dev.copy2pdf(file = "../../graphs/validation_moose.pdf")


Hv_eq = apply(RTBMH[,], 1,  Hv_eq_solve) 

par(mfrow = c(1, 2))
plot(Hv_eq/10/75~RTBMH[,"B"], main="Deer biomasses at equilibrium", xlab = "Boreal proportion", ylab = "#nb indiv (75kg) per km2", pch=4, cex=.5)
points(Hv_eq[nrow(RTBMH)]/10/75~RTBMH[nrow(RTBMH),"B"], col="red")
plot(Hv_eq/10/75~RTBMH[,"R"], main="Deer biomasses at equilibrium", xlab = "Immature forest proportion", ylab = "#nb indiv (75kg) per km2", pch=4, cex=.5)
points(Hv_eq[nrow(RTBMH)]/10/75~RTBMH[nrow(RTBMH),"R"], col="red")

dev.copy2pdf(file = "../../graphs/validation_deer.pdf")



#
#dat = read.table("../data/data_reshaped_RBTM_herbivores.txt")
#table(dat$st0)/nrow(dat)
#table(dat$st1)/nrow(dat)
#

# solve equilibrium main function
#eq_H = function(RTBM)
#{
#print(RTBM)
#RTBM = c(.2, .2, 0.4, 0.2)  # debug
#R = as.numeric(RTBM[1])
#T = as.numeric(RTBM[2])
#B = as.numeric(RTBM[3])
#M = as.numeric(RTBM[4])
#
# first tries to test if there is a interval including zero (a solution)
#Hvx = exp(seq(-10,9,0.2))
#res = unlist(lapply(Hvx, eq_fct, RTBM=as.numeric(RTBM[1:4])))
#plot(Hvx, res, type="b")
#abline(h=0)
#
#res[res==-999] = NA
#
#if(sum(is.na(res))!=length(res)) # check if there is no special case problems
#{
#print(RTBM)
#
#if(sum(res==0)>1) 
# res == 0 pour infinite Hv>0, which means Ha=0 is one solution
#{
#Ha_eq= 0
#
# find Hv back pour Ha_eq=0
#
#hvmin = as.numeric(Iv_Ikv(0.0001, Ha=Ha_eq, R=R, T=T, B=B, M=M, Ikv=Ikv))
#hvmax = as.numeric(Iv_Ikv(10000, Ha=Ha_eq, R=R, T=T, B=B, M=M, Ikv=Ikv))
#if(hvmin*hvmax<0) # if it is possible to cross the zero line
#{
#Hv_eq = uniroot(Iv_Ikv, Ha=Ha_eq, R=R, T=T, B=B, M=M, Ikv=Ikv, interval=c(0.0001, 10000), f.lower = hvmin, f.upper = hvmax, tol = 0.001)$root
#} else Hv_eq = 0;
#
#
#
#if(max(res, na.rm=TRUE)>0 & min(res, na.rm=TRUE)<0 ) 
# check if there is a unique solution
#{
#    imin = Hvx[which.min(res)]
#    imax = Hvx[which.max(res)]
#    hvmin = eq_fct(imin, as.numeric(RTBM[1:4]))
#    hvmax = eq_fct(imax, as.numeric(RTBM[1:4]))
#    if(hvmin*hvmax<0) 
#    {
#    Hv_eq = uniroot(eq_fct, RTBM=as.numeric(RTBM[1:4]), interval=c(imin, imax))$root
#    } else { #no coexistence
#    if(eq_fct(0, RTBM=as.numeric(RTBM[1:4]))==0) {
#    Hv_eq = 0
#    }else Hv_eq=NA
#    #}
#}else{Hv_eq = NA}
#
# find Ha back pour Hv_eq
#if(!is.na(Hv_eq))
#{
#hamin = as.numeric(Ia_Ika(0.0001, Hv=Hv_eq, R=R, T=T, B=B, M=M))
#hamax = as.numeric(Ia_Ika(10000, Hv=Hv_eq, R=R, T=T, B=B, M=M))
#if(hamin*hamax<0) # if it is possible to cross the zero line
#{
#Ha_eq = uniroot(Ia_Ika, Hv=Hv_eq, R=R, T=T, B=B, M=M, interval=c(0.0001, 10000), f.lower = hamin, f.upper = hamax, tol = 0.001)$root
#} else Ha_eq = 0;
#
#}else {
#Ha_eq = NA
#}
#
#}else {Ha_eq = Hv_eq=NA}
#
#
#return(c(Ha_eq, Hv_eq))
#}
#
# apply function to find Hb that solve the equation


#Hstars = t(apply(RTBM[testIndex,], 1,  eq_H)); colnames(Hstars) = c("Ha", "Hv");
#Hstars = cbind(RTBM[testIndex,], Hstars)
#

# plot
#---------
#summary(Hstars)
#
#par(mfrow = c(2, 2))
#
#hist(Hstars$Ha, main="Moose")
#hist(Hstars$Hv, main="Deer")
#
#tab = Hstars[order(Hstars$B),]
#
#meanHa = lapply(split(tab$Ha, tab$B), mean)
#meanHv = lapply(split(tab$Hv, tab$B), mean)
#Bax = as.numeric(as.character(names(meanHa)))
#
#
#pdf("../graphs/herbivore_coexistence.pdf")
#
#plot(Bax , meanHa, main="Herbivores biomasses at equilibrium", xlab = "Boreal proportion", ylab = "Herbivore biomass", type ="b", ylim = c(0, max(Hstars[,c("Ha","Hv")])))
#
#points(Bax, meanHv, col=2, type = "b", pch = 4)
#
#legend("topright", legend = c("Moose", "Deer"), col = c(1, 2), lwd = 1)
#
#dev.off()
#
#
