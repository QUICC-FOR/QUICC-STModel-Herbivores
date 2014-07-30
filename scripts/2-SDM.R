rm(list=ls())
# Load data
dataProj = as.data.frame(read.table("../data/data_reshaped_RBTM_herbivores.txt"))
dataProj$E = dataProj$av_annual_mean_tp
dataProj$P = dataProj$av_annual_pp

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
library(nnet)
SDM1 = multinom(state ~ poly(E,3,raw=TRUE) + poly(P,3,raw=TRUE) + E:P +lon + lat, data = calib)

pred1 = predict(SDM1, new=valid,"class")
(HK1 = HK(pred1, valid$state)) 
pred_proba1 = predict(SDM1,new=valid,"prob")

#--
library(randomForest)
set.seed(23)
SDM2 = randomForest(state ~ . , data = data[, c("state", "E", "P","lat", "lon")], ntree = 500)

set.seed(23)
pred2 = predict(SDM2,new=valid,"response", OOB=TRUE)
(HK2 = HK(pred2, valid$state)) 

set.seed(23)
pred_proba = predict(SDM2,new=valid,"prob", OOB=TRUE)

##-- graphs
par(mfrow = c(2,2))
plot(lowess(cbind(pred_proba[,"B"], pred_proba1[,"B"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Boreal")
abline(0,1,lty=2)
plot(lowess(cbind(pred_proba[,"T"], pred_proba1[,"T"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Temperate")
abline(0,1,lty=2)
plot(lowess(cbind(pred_proba[,"M"], pred_proba1[,"M"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Mixed")
abline(0,1,lty=2)
plot(lowess(cbind(pred_proba[,"R"], pred_proba1[,"R"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Early succession")
abline(0,1,lty=2)

pred_proba = pred_proba1 # pour graphs multinom

par(mfrow = c(2,2))
plot(lowess(cbind(valid$E,pred_proba[,"B"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Boreal probability", ylim = c(0,1))
plot(lowess(cbind(valid$E,pred_proba[,"T"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Temperate probability", ylim = c(0,1))
plot(lowess(cbind(valid$E,pred_proba[,"M"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Mixed probability", ylim = c(0,1))
plot(lowess(cbind(valid$E,pred_proba[,"R"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Early succession probability", ylim = c(0,1))

dev.copy2pdf(file = "../figures/SDM_temperature_randomForest.pdf")

par(mfrow = c(2,2))
boxplot(split(pred_proba[,"B"], ifelse(data$state=="B", "dominant", "non-dominant")), main="Boreal", ylab = "predicted probability")
boxplot(split(pred_proba[,"T"], ifelse(data$state=="T", "dominant", "non-dominant")), main="Temperate", ylab = "predicted probability")
boxplot(split(pred_proba[,"M"], ifelse(data$state=="M", "dominant", "non-dominant")), main="Mixed", ylab = "predicted probability")
boxplot(split(pred_proba[,"R"], ifelse(data$state=="R", "dominant", "non-dominant")), main="Early succession", ylab = "predicted probability")

dev.copy2pdf(file = "../figures/SDM_validation_randomForest.pdf")




#---------------------------------------------------------------------
# Get predictions
#newdat = expand.grid(E=seq(-1.7,7,0.1),P=seq(700, 1600, 10) )
SDM = SDM2 #rf
SDM = SDM1 #multinom

set.seed(23)
pred_real = predict(SDM,new=dataProj,"prob")

write.table(pred_real,"../data/data_pred_states_randomForest.txt")
write.table(pred_real,"../data/data_pred_states_multinom.txt")


