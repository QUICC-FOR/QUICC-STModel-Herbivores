library(nnet)
# Load data
data = as.data.frame(read.table("../data/data_categorical_RBTM.txt"))
data$E = data$av_annual_mean_tp
data$P = data$av_annual_pp

# Run the model
#SDM = multinom(t1 ~ E + I(E^2) + I(E^3) + P + I(P^2) + I(P^3), data)
library(randomForest)
set.seed(23)
SDM = randomForest(t1 ~ . , data = data[, c("t1", "E", "P")])

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



