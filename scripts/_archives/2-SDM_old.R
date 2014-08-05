rm(list=ls())
# Load data
data = as.data.frame(read.table("../data/data_categorical_RBTM.txt"))
data$E = data$av_annual_mean_tp
data$P = data$av_annual_pp
data$state = data$t1
str(data)


data = as.data.frame(read.table("../data/data_allyears_RBTM.txt"))
data$E = data$annual_mean_temp
data$P = data$annual_pp
data$state = data$pred_class
data = data[-which(data$state=="Unclass"),]
data$state = as.factor(as.character(data$state))
str(data)

# cross validation
source("BoulangeatEcoLet2012_fct.r")
sampl = sample(1:nrow(data), 2*nrow(data)/3)
calib = data #data[sampl,] 
valid = data #data[-sampl,]

# Run the models
#library(nnet)
#SDM1 = multinom(state ~ poly(E,3,raw=TRUE) + poly(P,3,raw=TRUE) + E:P , data = calib)
#SDM0 = multinom(state ~ 1, calib[,c("state", "E", "P")])
#
#pred1 = predict(SDM1, new=valid,"class")
#(HK1 = HK(pred1, valid$state)) 


library(randomForest)
set.seed(23)
SDM2 = randomForest(state ~ . , data = data[, c("state", "E", "P","lat", "lon")], ntree = 500)

set.seed(23)
pred2 = predict(SDM2,new=valid,"response", OOB=TRUE)
(HK2 = HK(pred2, valid$state)) 

set.seed(23)
pred_proba = predict(SDM2,new=valid,"prob", OOB=TRUE)

par(mfrow = c(2,2))
plot(lowess(cbind(valid$E,pred_proba[,"B"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Boreal probability", ylim = c(0,1))
plot(lowess(cbind(valid$E,pred_proba[,"T"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Temperate probability", ylim = c(0,1))
plot(lowess(cbind(valid$E,pred_proba[,"M"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Mixed probability", ylim = c(0,1))
plot(lowess(cbind(valid$E,pred_proba[,"R"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Early succession probability", ylim = c(0,1))

dev.copy2pdf(file = "../figures/SDM_temperature_randomForest.pdf")


#library(party) # plus long que randomForest
#set.seed(23)
#SDM3 = cforest(state ~  ., data = calib[, c("state", "E", "P")], control = cforest_unbiased(mtry = 2))
#set.seed(23)
#pred3 = predict(SDM3, newdata=valid[, c("state", "E", "P")], OOB = TRUE)
#
#(HK3 = HK(pred3, valid$state)) 


#---------------------------------------------------------------------
# Get predictions
newdat = expand.grid(E=seq(-1.7,7,0.1),P=seq(700, 1600, 10) )

#pred_real = predict(SDM,new=data.frame(E=data$E, P=data$P),"probs")
#pred_gradient = predict(SDM,new=data.frame(newdat),"probs")

set.seed(23)
pred_real = predict(SDM,new=data.frame(E=data$E, P=data$P),"prob")
set.seed(23)
pred_gradient = predict(SDM,new=data.frame(newdat),"prob")

write.table(pred_real,"../data/data_pred_states_RBTM.txt")
write.table(cbind(newdat,pred_gradient),"../data/data_pred_gradient_RBTM.txt")

# Plot the results
colo = c("darkgreen", "lightgreen", "orange", "blue")

plotfct <- function(env_var, leg)
{
Cmax = lapply(split(pred_gradient[,"B"], as.factor(env_var)), max)
Dmax = lapply(split(pred_gradient[,"T"], as.factor(env_var)), max)
Mmax= lapply(split(pred_gradient[,"M"], as.factor(env_var)), max)
Tmax = lapply(split(pred_gradient[,"R"], as.factor(env_var)), max)
pred_gradient_max = cbind(Cmax, Dmax, Mmax, Tmax)

par(mar=c(5,5,2,1))
xscale = as.numeric(rownames(pred_gradient_max))
plot(xscale,pred_gradient_max[,1],type = "l",ylim=c(0,1),cex.axis = 1.25, cex.lab = 1.25, xlab = leg, ylab = "Proportion",lwd = 2, col=colo[1])
lines(xscale,pred_gradient_max[,2],col = colo[2],lwd = 2)
lines(xscale,pred_gradient_max[,3],col = colo[3],lwd = 2)
lines(xscale,pred_gradient_max[,4],col = colo[4],lwd = 2)
abline(v = 7,lty = 3)
legend("topright",bty = "n", col = colo,legend = c("B","T","M","R"),lty=1)
}

plotfct (newdat$E ,leg = "Température moyenne annuelle")
dev.copy2pdf(file = "../figures/SDM_temp_RBTM.pdf")

plotfct (env_var = newdat$P , leg = "Précipitations moyennes annuelles")
dev.copy2pdf(file = "../figures/SDM_precip_RBTM.pdf")



