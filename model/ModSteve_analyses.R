rm(list=ls())
# Range kutta 
#------------------------------------------
library(deSolve)

#------------------------------------------
# no herbivore
#------------------------------------------
env_seq = seq(-3, 3, by=.1)
comp0 = 0.8
comp_seq = (comp0-0.2)*(1-pnorm(env_seq))+0.2
k_seq = 0.5 + env_seq/(max(env_seq)-min(env_seq))

pdf("../graphs/env_gradient2.pdf", height=4, width=10)
cexax = 1
par(mfrow = c(1,2), cex.lab=1.4)
plot(comp_seq~env_seq, type ="l", ylim = c(0,1), xaxt="n",  ylab = "seedling competitiveness (c)", main = "", xlab="temperate->boreal")
abline(h=comp0, lty=2)
abline(h=0.2, lty=2)

plot(k_seq~env_seq, type ="l", ylim = c(0,1), xaxt="n",  ylab = "conifer competitiveness (k)", main = "", xlab="temperate->boreal")
abline(h=0.5, lty=2)

dev.off()

resNoH = t(apply(cbind(k_seq, comp_seq), 1, runNoH))

colnames(resNoH) = c("C", "D", "S", "G")

pdf("../graphs/Mod2_climate.pdf", width=4, height=4)

res = resNoH
colo = c("darkgreen", "darkviolet", "brown", "orange", "blue", "darkred")
plot(res[,"C"], type='l', xlab="temperate->boreal", xaxt="n", ylab ="Proportions", ylim = c(0,1), main = "Vegetation", col=colo[1])
lines(res[,"D"], col=colo[2])
lines(res[,"S"], col=colo[3])
lines(res[,"G"], col=colo[4])
axis(1, at = c(0,100), labels = c(0,1))
legend("topleft", bty="n", legend = c("Conifers", "Deciduous", "Seedlings", "Grasses"), col=colo[1:4], lwd=1, cex=.6)


dev.off()


#-----------------------------------------------------------
#
# Herbivore / constant (no dynamic)
#
#----------------------------------------------------------
resHb = t(data.frame(lapply(seq(0, 5000, by=10), run)))
colnames(resHb) = c("C", "D", "S", "G")
rownames(resHb) = seq(0, 5000, by=10)
resHb

res0 = resHb[1,]


pdf("../graphs/Mod2_herbivore.pdf", width=8, height=4)

par(mfrow = c(1, 2), lwd=1.2)
colo = c(1, "darkgreen", "saddlebrown", "darkorange")

res=resHb
plot(rownames(res), res[,"C"], type='l', xlab="Browser abundance", ylab ="Proportions", ylim = c(0,1), main = "", col=colo[1])
lines(rownames(res),res[,"D"], col=colo[2])
lines(rownames(res),res[,"S"], col=colo[3])
lines(rownames(res),res[,"G"], col=colo[4])
legend("topright", bty="n", legend = c("Boreal", "Temperate", "Seedlings", "Grasses"), col=colo[1:4], lwd=1, cex=1)


dat = rbind(run0, apply(resHb, 2, mean))
lo = rbind(rep(NA,4), apply(resHb, 2, function(x){mean(x)-sd(x)}))
up = rbind(rep(NA,4), apply(resHb, 2, function(x){mean(x)+sd(x)}))

par(mar = c(7,4,1,1))
barplot2(as.matrix(dat), col = rep(colo[1:4], each=2), beside=TRUE, border=NA, ylab = "proportions", ylim = c(0,1), space = c(0.2, 1), plot.grid = TRUE, grid.lty = 1, grid.lwd = 1, grid.col = "white", bty='n', grid.inc = 10,  prcol = "lightgrey", plot.ci=TRUE, ci.l=lo, ci.u = up, names.arg = rep(c("No herbivore", "Browser"), time=4), ci.width = .2, las = 3 )
legend("topright", bty="n", legend = c("Boreal", "Temperate", "Seedlings", "Grasses"), col=colo[1:4], pch=15, cex=1)

dev.off()


#-----------------------------------------------------------
#
#  Combined model for a set of parameters
#
#----------------------------------------------------------


resHb = runHb(c(0.8, 1, 0.5))


resNoH = runNoH(0.5) ## cf Mod2_NoHerbivore_Climate

HH_base = herbi (as.vector(resNoH[4:1]))  ## effectOfVegeOnHerb


pdf("../graphs/Mod2H_combinated.pdf", width=8, height=7)

par(mfrow = c(2, 1), mar = c(4,4,1,1), lwd=1.2)
colo = c("black", "dodgerblue", "orangered", "purple")

dat = rbind(c(resNoH, HH_base), resHb)
eff = apply(dat, 2, function(x){log(x[2]/x[1])})


barplot2(as.matrix(dat[,1:4]), col = rep(colo[1:2], 4), beside=TRUE, border=NA, ylab = "proportions", ylim = c(0,.8), space = c(0.2, 1), plot.grid = TRUE, grid.lty = 1, grid.lwd = 1, grid.col = "white", bty='n', grid.inc = 10,  prcol = "lightgrey", names.arg = c("Boreal trees", "Temperate trees", "Seedling", "Grasses"))
legend("topright", bty="n", legend = c("No feedback", "Browser"), col=colo[1:2], pch=15, cex=1)


barplot2(as.matrix(eff), col = colo[2], beside=TRUE, border=NA, ylab = "herbivore effect", ylim = c(-.50,.5), space = c(0.2, 1), plot.grid = TRUE, grid.lty = 1, grid.lwd = 1, grid.col = "white", bty='n', grid.inc = 10,  prcol = "lightgrey", names.arg = c("Boreal trees", "Temperate trees", "Seedling", "Grasses", "Browser"))
abline(h=0, lwd=1)


dev.off()


#-----------------------------------------------------------
#
#  Runs and graphs -> environmental gradient & coupled model
#
#----------------------------------------------------------
env_seq = seq(-3, 3, by=.1)
comp0 = 0.8
comp_seq = comp0*(1-pnorm(env_seq))
uvw = .6 - .8*env_seq/(max(env_seq)-min(env_seq))
k_seq = 0.5-env_seq/(max(env_seq)-min(env_seq))

env_list = cbind(comp_seq, uvw, k_seq)


## runs
resHb = t(apply(env_list, 1, runHb))


## effect pov veget

resNoH = t(apply(cbind(k_seq, comp_seq), 1, runNoH))

effHb = log(resHb[,1:4]/resNoH)


## effect pov herbi


HH_base_Hb = t(apply(resHb[,4:1],1,  herbi))  ## effectOfVegeOnHerb

effHb2 = log(resHb[,4]/HH_base_Hb)

effHb2[effHb2==max(effHb2)]=NA

pdf("../graphs/Mod2H_combinated_env.pdf", width=8, height=4)

par(mfrow = c(1, 2), mar = c(4,4,4,1), lwd=1.2)
colo = c(1, "darkgreen", "brown", "orange", "blue", "darkred")

ploteff = function(eff, name, colo){
if(is.null(dim(eff))) 
{
eff1 = eff
} else eff1 =eff[,1]
plot(env_seq, eff1, ylim = c(min(eff, na.rm=T), max(eff, na.rm=T)), type ="l", col = colo[1], xaxt="n", xlab = "hot and wet --> cold and dry", ylab = "Feedback effect", main=name)
abline(h=0, lty=2, col ="grey", lwd=1)
try(lines(env_seq, eff[,2], col = colo[2]))
try(lines(env_seq, eff[,3], col = colo[3]))
try(lines(env_seq, eff[,4], col = colo[4]))

}

ploteff(effHb, "Vegetation", colo = colo[1:4])
legend("bottomright", bty="n", legend = c("Boreal trees", "Temperate trees", "Sedlings", "Grassland"), col=colo[1:4], lwd=1.2, cex=1)


plot(env_seq, effHb2, type ="l", col = colo[5], xaxt="n", xlab = "hot and wet --> cold and dry", ylab = "Feedback effect", main="Browser")
abline(h=0, lty=2, col ="grey", lwd=1)

dev.off()