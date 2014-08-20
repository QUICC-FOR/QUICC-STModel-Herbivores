rm(list=ls())
# Range kutta 
#------------------------------------------
library(deSolve)

source("Mod1_functions.R")

#-----------------------------------------------------------
#
# graphs -> environmental gradient & coupled model
#
#----------------------------------------------------------
#-----------------------------------------------------------
#
# runs
#
#----------------------------------------------------------
env_seq = seq(-3, 3, by=.1)
comp0 = 1
comp_seq = comp0*(1-pnorm(env_seq))
#comp_seq = rep(comp0, length(env_seq))
#uvw = .6 - .8*env_seq/(max(env_seq)-min(env_seq))
uvw = rep(1, length(env_seq))
env_list = cbind(comp_seq, uvw)

plot(comp_seq~env_seq, type="l", xlab = "Env gradient (T° decreases)", ylab = "colonisation of seedlings (c)")
#plot(uvw~env_seq)
dev.copy2pdf(file="../graphs/climate.pdf")



## runs
resHb = t(apply(env_list, 1, runHb))
resHg = t(apply(env_list, 1, runHg))
resHH = t(apply(env_list, 1, runHH))


## effect pov veget
resNoH = t(apply(env_list, 1, runNoH))


effHb = log(resHb/resNoH)
effHg = log(resHg/resNoH)
effHH = log(resHH/resNoH)


# effect pov herbi


HH_base_HH = t(apply(resHH[,c(3,1,2)],1,  coexistence))  ## effectOfVegeOnHerb
HH_base_Hb = t(apply(resHb[,c(3,1,2)],1,  coexistence))  ## effectOfVegeOnHerb
HH_base_Hg = t(apply(resHg[,c(3,1,2)],1,  coexistence))  ## effectOfVegeOnHerb
colnames(HH_base_Hb) = colnames(HH_base_Hg) = colnames(HH_base_HH)= c("Hb", "Hg")


save.image("runs_env.RData")

#-----------------------------------------------------------
#
# graphs
#
#----------------------------------------------------------
load("runs_env.RData")



#--------------------------
# graphs presentation ESA
#---------------------------
colo = c("darkgreen", "lightgreen", "orange", "blue", "red")
par(lwd = 2)
#par(mfrow = c(2,2), mar = c(5,5,5,5))

plot(resNoH[,"T"], type='l', xlab="Env gradient (T° decreases)", xaxt="n", ylab ="Proportions", ylim = c(0,1), main = "", col=colo[1])
lines(resNoH[,"S"], col=colo[2])
lines(resNoH[,"G"], col=colo[3])
axis(1, at = c(0,100), labels = c(0,1))
legend("topleft", bty="n", legend = c("Trees", "Seedlings", "Grasses"), col=colo[1:3], lwd=1, cex=.8)

dev.copy2pdf(file = "../graphs/graphs_vegetDyn.pdf")

plot(1, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt="n", bty="n")
legend("center", bty="n", legend = c("Trees", "Seedlings", "Grasses", "Grazer", "Browser"), col=colo[1:5], lty=c(1,1,1,2,2), cex=1)

#par(mfrow = c(2,2), mar = c(5,5,5,5))

plot(resHb[,"T"], type='l', xlab="decreasing temperature", xaxt="n", ylab ="Proportions", ylim = c(0,1), main = "With Browser", col=colo[1])
lines(resHb[,"S"], col=colo[2])
lines(resHb[,"G"], col=colo[3])
axis(1, at = c(0,100), labels = c(0,1))
#legend("topleft", bty="n", legend = c("Trees", "Seedlings", "Grasses"), col=colo[1:3], lwd=1, cex=.8)

lines(resHb[,"Hb"]/max(resHb[,"Hb"]), col = colo[5], lty=2)
axis(4, at = c(0,0.5,1), labels = c(0, round(max(resHb[,"Hb"])/2), round(max(resHb[,"Hb"]))), las = 3)
mtext("Herbivore biomass", 4, line =2, cex=.8) 


plot(resHg[,"T"], type='l', xlab="decreasing temperature", xaxt="n", ylab ="Proportions", ylim = c(0,1), main = "With Grazer", col=colo[1])
lines(resHg[,"S"], col=colo[2])
lines(resHg[,"G"], col=colo[3])
axis(1, at = c(0,100), labels = c(0,1))
#legend("topleft", bty="n", legend = c("Trees", "Seedlings", "Grasses", "Grazer"), col=colo[1:4], lwd=1, cex=.8)

lines(resHg[,"Hg"]/max(resHg[,"Hg"]), col = colo[4], lty=2)
axis(4, at = c(0,0.5,1), labels = c(0, round(max(resHg[,"Hg"])/2), round(max(resHg[,"Hg"]))), las = 3)
mtext("Herbivore biomass", 4, line =2, cex=.8) 

dev.copy2pdf(file="../graphs/effectOnVeget_compare.pdf")


#----------------------------------------------------------

par(mfrow = c(1,2), mar = c(5,5,5,1))

plot(sign(log((resHb[,"T"]+resHb[,"S"])/(resNoH[,"T"]+resNoH[,"S"]))), type='l', xlab="decreasing seedling performance", xaxt="n", ylab ="Effect on closed vegetation (S+T)", main = "Browser", col=colo[5])
abline(h=0, lty=2)

plot(sign(log((resHg[,"T"]+resHg[,"S"])/(resNoH[,"T"]+resNoH[,"S"]))), type='l', xlab="decreasing seedling performance", xaxt="n", ylab ="Effect on closed vegetation (S+T)", main = "Grazer", col=colo[4])
abline(h=0, lty=2)

dev.copy2pdf(file="../graphs/effectOnVeget_sign_single.pdf")

par(mfrow = c(1,2), mar = c(5,5,5,1))

plot(sign(log((resHH[,"T"]+resHH[,"S"])/(resNoH[,"T"]+resNoH[,"S"]))), type='l', xlab="decreasing seedling performance", xaxt="n", ylab ="Effect on closed vegetation (S+T)", main = "Browser x Grazer", col=1)
abline(h=0, lty=2)

dev.copy2pdf(file="../graphs/effectOnVeget_sign_both.pdf")

#----------------------------------------------------------

par(mfrow = c(1,2), mar = c(5,5,5,1))

plot(HH_base_HH[,"Hb"], type="l", col = colo[5],xlab="decreasing seedling performance", xaxt="n", ylab = "Herbivore biomass", ylim = c(0, 1000), main = "Browser")
lines(resHH[,"Hb"], type="l", col = colo[5], lty = 2)
legend("topleft", bty="n", legend = c("No feedback on vegetation", "Coupled model"), lty = c(1, 2), lwd=1.2, cex=1, col = colo[5])

plot(HH_base_HH[,"Hg"], type="l", col = colo[4],xlab="decreasing seedling performance", xaxt="n", ylab = "Herbivore biomass", ylim = c(0, 1500), main = "Grazer")
lines(resHH[,"Hg"], type="l", col = colo[4], lty = 2)
legend("topleft", bty="n", legend = c("No feedback on vegetation", "Coupled model"), lty = c(1, 2), lwd=1.2, cex=1, col = colo[4])

dev.copy2pdf(file="../graphs/feedback_importance.pdf")

#----------------------------------------------------------

par(mfrow = c(2,2), mar = c(5,5,5,1))

plot(resHH[,"T"], type='l', xlab="decreasing temperature", xaxt="n", ylab ="Proportions", ylim = c(0,1), main = "Both herbivores", col=colo[1])
lines(resHH[,"S"], col=colo[2])
lines(resHH[,"G"], col=colo[3])
axis(1, at = c(0,100), labels = c(0,1))
#legend("topleft", bty="n", legend = c("Trees", "Seedlings", "Grasses"), col=colo[1:3], lwd=1, cex=.8)

#ymax = max(resHH[,"Hg"],resHH[,"Hb"])
#lines(resHH[,"Hg"]/ymax, col = colo[4], lty=2)
#
#lines(resHH[,"Hb"]/ymax, col = colo[5], lty=2)
#axis(4, at = c(0,0.5,1), labels = c(0, round(ymax/2), round(ymax)), las = 3)

#legend("top", bty="n", legend = c("Grazer", "Browser"), col=colo[4:5], lwd=1, lty=2,cex=.8)


plot(1, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt="n", bty="n")
legend("center", bty="n", legend = c("Trees", "Seedlings", "Grasses", "Grazer", "Browser"), col=colo[1:5], lty=c(1,1,1,2,2), cex=1)



plot(resHb[,"Hb"], type="l", col = colo[5],xlab="decreasing seedling performance", xaxt="n", ylab = "Browser biomass", ylim = c(0, 1500))
lines(resHH[,"Hb"], type="l", col = colo[5], lty = 2)
legend("topleft", bty="n", legend = c("Alone", "Interaction"), lty = c(1, 2), lwd=1.2, col = colo[5], cex=1)

plot(resHg[,"Hg"], type="l", col = colo[4],xlab="decreasing seedling performance", xaxt="n", ylab = "Grazer biomass", ylim = c(0, 1500))
lines(resHH[,"Hg"], type="l", col = colo[4], lty = 2)
legend("topleft", bty="n", legend = c("Alone", "Interaction"), lty = c(1, 2), lwd=1.2, col = colo[4], cex=1)

dev.copy2pdf(file="../graphs/competition_facilitation.pdf")



