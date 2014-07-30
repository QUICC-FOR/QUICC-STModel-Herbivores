rm(list=ls())
# Load data
data = as.data.frame(read.table("../data/data_categorical_RBTM.txt"))
#data$E = data$av_annual_mean_tp
#data$P = data$av_annual_pp
#data$state = data$t1
#str(data)


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


par(mfrow = c(2,2))
plot(lowess(cbind(valid$E,pred_proba1[,"B"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Boreal probability", ylim = c(0,1))
plot(lowess(cbind(valid$E,pred_proba1[,"T"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Temperate probability", ylim = c(0,1))
plot(lowess(cbind(valid$E,pred_proba1[,"M"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Mixed probability", ylim = c(0,1))
plot(lowess(cbind(valid$E,pred_proba1[,"R"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Early succession probability", ylim = c(0,1))

dev.copy2pdf(file = "../figures/SDM_temperature_multinom.pdf")

par(mfrow = c(2,2))
boxplot(split(pred_proba1[,"B"], ifelse(data$state=="B", "dominant", "non-dominant")), main="Boreal", ylab = "predicted probability")
boxplot(split(pred_proba1[,"T"], ifelse(data$state=="T", "dominant", "non-dominant")), main="Temperate", ylab = "predicted probability")
boxplot(split(pred_proba1[,"M"], ifelse(data$state=="M", "dominant", "non-dominant")), main="Mixed", ylab = "predicted probability")
boxplot(split(pred_proba1[,"R"], ifelse(data$state=="R", "dominant", "non-dominant")), main="Early succession", ylab = "predicted probability")

dev.copy2pdf(file = "../figures/SDM_validation_multinom.pdf")


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


#library(party) # plus long que randomForest
#set.seed(23)
#SDM3 = cforest(state ~  ., data = calib[, c("state", "E", "P")], control = cforest_unbiased(mtry = 2))
#set.seed(23)
#pred3 = predict(SDM3, newdata=valid[, c("state", "E", "P")], OOB = TRUE)
#
#(HK3 = HK(pred3, valid$state)) 


#---------------------------------------------------------------------
# Get predictions
#newdat = expand.grid(E=seq(-1.7,7,0.1),P=seq(700, 1600, 10) )
SDM = SDM2

set.seed(23)
pred_real = predict(SDM,new=data,"prob")
#set.seed(23)
#pred_gradient = predict(SDM,new=data.frame(newdat),"prob")

write.table(pred_real,"../data/data_pred_states_randomForest.txt")
#write.table(cbind(newdat,pred_gradient),"../data/data_pred_gradient_RBTM.txt")

SDM = SDM1

pred_real = predict(SDM,new=data,"prob")

write.table(pred_real,"../data/data_pred_states_multinom.txt")
