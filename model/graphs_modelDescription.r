#------------------------------------------
# Herbivore mortality
#------------------------------------------
Intake = seq(0,15,0.01)
p = 5.84; z = 1; cc=0.7
M = p+exp(-z*cc*Intake)

pdf("../graphs/herbivore_mortality.pdf", height=3, width=3)

par(mar=c(4,4,0.1,0.1), cex=.9)
plot(M~Intake, type ="l", yaxt="n",ylab = "biomass loss (kg/kg/year)", xlab = "Biomass gain (kg/kg/year)", bty="n")
abline(h=p, lty=2)
axis(2, at = p, labels="p", las = 1, cex.axis =0.7)
abline(h=1+p, lty=2)
axis(2, at = 1+p, labels="1+p", las = 1, cex.axis =0.7)

abline(h=0)
abline(v=0)

dev.off()

## mortality en fonction de B

BB = seq(0,1,0.01)
m0 = 0.1
ms = 0.5
ma = ms/(1+(ms/m0-1)*BB)

pdf("../graphs/shelter_effect.pdf", height=4, width=4)

plot(ma~BB, type = "l", ylim = c(0,1), xlab = c("Boreal forest proportion"), ylab = "m (for Moose)")
abline(h=m0, lty=2)
abline(h=ms, lty=2)
axis(2, at = m0, labels="m0", las = 1, cex.axis = 0.7, line = -0.8)
axis(2, at = ms, labels="ms", las = 1, cex.axis = 0.7, line = -0.8)

dev.off()
#------------------------------------------
# Intake rates
#------------------------------------------
F1 = seq(0,20,1)
tau = 17; mu = tau; 
I1 = tau*F1/(tau + F1)

pdf("../graphs/intakes.pdf", height=4, width=7)

par(mfrow = c(1,2))
plot(I1~F1, type ="l", ylim = c(0,tau*1.1), xaxt="n", yaxt="n", ylab = "Intake rate", xlab="Available resource", bty="n")
abline(0,1, col = 2)
abline(v=mu, lty=2)
axis(1, at = mu, labels="mu", las = 1, cex.axis =0.7)
abline(h=tau/2, lty=2)
axis(2, at = tau/2, labels="tau/2)", las = 1, cex.axis =0.7)
abline(h=tau, lty=2)
axis(2, at = tau, labels="tau", las = 1, cex.axis =0.7)

abline(h=0)
abline(v=0)
axis(1, at = 0, labels="0", cex.axis =0.7, tick=F)
axis(2, at = 0, labels="0",las = 1,  cex.axis =0.7, tick=F)

#------------------------------------------
r=5; 
I1 = seq(0,10,0.01)
phi = 0.1
#I2 = 1/(1+exp(r*(I1-p)))
#I2 = (phi*(exp(p)-1)*exp(-p))/((exp(p)-2)*exp(-p) +exp(r*(I1-p)))
I2 = phi + (1-phi)/(1+exp(r*(I1-p)))

plot(I2~I1, type ="l", xaxt="n", yaxt="n", ylab = "Coefficient for non-preferred resource", xlab="Intake of preferred resource (kg/kg/year)", bty="n", col="purple")
abline(v=p, lty=2)

r=1
I2 = phi + (1-phi)/(1+exp(r*(I1-p)))
lines(I2~I1, col="red")

r=25
I2 = phi + (1-phi)/(1+exp(r*(I1-p)))
lines(I2~I1, col="dodgerblue2")

axis(2, at = 1, labels="1", las = 1, cex.axis =0.7)
axis(2, at = phi, labels="phi", las = 1, cex.axis =0.7)
axis(2, at = (1+phi)/2, labels="(1+phi)/2", las = 1, cex.axis =0.7)

axis(1, at = p, labels="p", cex.axis =0.7)

abline(h=1, lty=2)
abline(h=phi, lty=2)
abline(h=(1+phi)/2, lty=2)

abline(v=0)
axis(1, at = 0, labels="0", cex.axis =0.7, tick=F)
axis(2, at = 0, labels="0", las = 1, cex.axis =0.7, tick=F)


dev.off()


#------------------------------------------
# competition
#------------------------------------------
Ha = Hv = seq(0,1,0.1)
Htable = expand.grid(Ha, Hv)
kk = unlist(apply(Htable, 1, function(x){exp(-log(1/0.5)*x[2]/x[1])}))

image(Ha, Hv, matrix(kk, ncol=length(Ha)), col=rev(heat.colors(10)))


#------------------------------------------
# effects on beta
#------------------------------------------
PTT = PTB = seq(0,1,0.1)
PTTable = expand.grid(PTT, PTB)

beta0 = 0.5

kk = unlist(apply(PTTable, 1, function(x){
exp(-log(1/beta0)*((1-x[1])/(1-x[2])))
}))

kk = unlist(apply(PTTable, 1, function(x){
(1/2)*((beta0)*cos(x[2]*pi) + ((1-beta0))*cos(pi*(1+x[1])) + 1 )
}))

sel = PTTable[,1] ==0
plot(PTTable[sel,2], kk[sel])

sel = PTTable[,2] ==0
plot(PTTable[sel,1], kk[sel])

#------------------------------------------
# effects on c and alpha
#------------------------------------------
#UG = US = seq(0,1, 0.1)
#c0 = 0.5; 
#cb = c0 * (1-US)
#cg = UG*(1-c0) + c0
#
#pdf("../graphs/effectsOnVeg.pdf", height=4, width=10)
#
#cexax = 1
#par(mfrow = c(1,3), cex.lab=1.4)
#plot(cb~US, type ="l", ylim = c(0,1), xaxt="n", yaxt="n", ylab = "c", main = "S/G competition \nwhen the grazer is absent", xlab="Herbivore pressure on S")
#abline(h=c0, lty=2)
#axis(2, at = c0, labels="c0", las = 1, cex.axis =cexax)
#axis(1, at = 0, labels="0", cex.axis =cexax, tick=F)
#axis(2, at = 0, labels="0", las = 1, cex.axis =cexax, tick=F)
#axis(1, at = 1, labels="1", cex.axis =cexax, tick=F)
#axis(2, at = 1, labels="1", las = 1, cex.axis =cexax, tick=F)
#
#plot(cg~UG, type ="l", ylim = c(0,1), xaxt="n", yaxt="n", ylab = "c", main = "S/G competition \nwhen the browser is absent", xlab="Herbivore pressure on G")
#abline(h=c0, lty=2)
#axis(2, at = c0, labels="c0", las = 1, cex.axis =cexax)
#axis(1, at = 0, labels="0", cex.axis =cexax, tick=F)
#axis(2, at = 0, labels="0", las = 1, cex.axis =cexax, tick=F)
#axis(1, at = 1, labels="1", cex.axis =cexax, tick=F)
#axis(2, at = 1, labels="1", las = 1, cex.axis =cexax, tick=F)
#
#plot(cb~US, type ="l", ylim = c(0,1), xaxt="n", yaxt="n", ylab = "alpha", main = "Succession rate", xlab="Herbivore pressure on S")
#abline(h=c0, lty=2)
#axis(2, at = c0, labels="alpha0", las = 3, cex.axis =cexax)
#
#axis(1, at = 0, labels="0", cex.axis =cexax, tick=F)
#axis(2, at = 0, labels="0", las = 1, cex.axis =cexax, tick=F)
#axis(1, at = 1, labels="1", cex.axis =cexax, tick=F)
#axis(2, at = 1, labels="1", las = 1, cex.axis =cexax, tick=F)
#
#dev.off()
#
#------------------------------------------
# effects of environment
#------------------------------------------
#env = seq(-3, 3, by=.1)
#
#comp0 = 0.8
#comp = comp0*(1-pnorm(env))
#
#k = pnorm(env)
#
#uvw = .6 - .8*env/(max(env)-min(env))
#
#pdf("../graphs/environmental_gradient.pdf", height=4, width=4)
#cexax = 1
#par(mfrow = c(1,1), cex.lab=1.4)
#
#plot(uvw~env, type ="l", ylim = c(0,1), xaxt="n",  ylab = "biomass factor", main = "", xlab="environmental gradient")
#abline(h=0.5, lty=2)
#
#
#dev.off()
#
#
#
