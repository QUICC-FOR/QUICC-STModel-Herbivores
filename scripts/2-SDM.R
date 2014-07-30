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

par(mfrow = c(2,2))
boxplot(split(pred_proba[,"B"], ifelse(data$state=="B", "dominant", "non-dominant")), main="Boreal", ylab = "predicted probability")
boxplot(split(pred_proba[,"T"], ifelse(data$state=="T", "dominant", "non-dominant")), main="Temperate", ylab = "predicted probability")
boxplot(split(pred_proba[,"M"], ifelse(data$state=="M", "dominant", "non-dominant")), main="Mixed", ylab = "predicted probability")
boxplot(split(pred_proba[,"R"], ifelse(data$state=="R", "dominant", "non-dominant")), main="Early succession", ylab = "predicted probability")

dev.copy2pdf(file = "../figures/SDM_validation_randomForest.pdf")




#---------------------------------------------------------------------
# Get predictions
#newdat = expand.grid(E=seq(-1.7,7,0.1),P=seq(700, 1600, 10) )
SDM = SDM2

set.seed(23)
pred_real = predict(SDM,new=dataProj,"prob")
#set.seed(23)
#pred_gradient = predict(SDM,new=data.frame(newdat),"prob")

write.table(pred_real,"../data/data_pred_states_randomForest.txt")
#write.table(cbind(newdat,pred_gradient),"../data/data_pred_gradient_RBTM.txt")

