###################################################################
#####                        Clean data                     #######
###################################################################

dat = read.table("../data/data_BA.txt",header = TRUE, sep = ";")

## Rm all disturbed plots
dat  <- dat[is.na(dat$disturbance),]
# BRP = brulis partiel
# CAM = coupe d'amélioration
# CB = coupe par bandes
# CD = coupes en damier
# CDL = coupe à diamètre limité
# CE = coupe partielle et épidémie légère
# CHP = chablis partiel
# CJ = coupe de jardinage
# CP = coupe partielle
# DLD = Coupe à diamètre limite avec dégagement des arbres d'avenir
# DP  = Dépérissement partiel du feuillu
# EL = Épidémie légère
# EPC = Éclaircie précommerciale
# VEP = Verglas partiel

## Order plot
dat  <- dat[order(dat$id_plot,dat$yr),]

## Rm all NA in cover type (R,M,F)
dat  <- dat[complete.cases(dat$cover_type),]

## Rm all plots with no climatic data associated
dat  <- dat[which(!is.na(dat$annual_pp)),]

## Conserve all plots with drainage 20,30,40
dat$drainage  <- as.numeric(dat$drainage)
test  = numeric(nrow(dat))
test[dat$drainage >= 20 & dat$drainage <= 41] = 1
dat = dat[test==1,]

### Get BA in hectares
dat[,7:62] <- dat[,7:62]*10000/400

## Rm plots with only one measurement
N  <- rowSums(table(dat$id_plot,dat$id_mes))
N  <- names(N[N==1])
dat  <- dat[!dat$id_plot %in% N,]

###################################################################
#####                  Classify plots                       #######
###################################################################

# List of interests species
#T_sp  <- c("boj","bop","chr","peb","peg","pet","pib","pir","prp","sal","soa")
#C_sp  <- c("epn","epb","epr","mel","pig","sab", "tho", "pru") 
#D_sp  <- c("err","ers","fra","frn","heg","osv","til","cet")

R_sp  <- c("boj","bop","chr","peb","peg","pet","pib","pir","prp","sal","soa")
B_sp  <- c("epn","epb","epr",  "mel","pig","sab") 
T_sp  <- c("err","ers","fra","frn","heg","osv","til","cet", "pru","tho")

# Subset species BA and cover type observed
class_dat  <- dat[,c(7:62,65)]
class_dat$sum_tot <- rowSums(dat[,7:62],na.rm=TRUE) 

#class_dat$T_prop  <- rowSums(dat[which(names(dat) %in% T_sp)],na.rm=T)/class_dat$sum_tot
#class_dat$C_prop  <- rowSums(dat[which(names(dat) %in% C_sp)],na.rm=T)/class_dat$sum_tot
#class_dat$D_prop  <- rowSums(dat[which(names(dat) %in% D_sp)],na.rm=T)/class_dat$sum_tot
class_dat$R_prop  <- rowSums(dat[which(names(dat) %in% R_sp)],na.rm=T)/class_dat$sum_tot
class_dat$B_prop  <- rowSums(dat[which(names(dat) %in% B_sp)],na.rm=T)/class_dat$sum_tot
class_dat$T_prop  <- rowSums(dat[which(names(dat) %in% T_sp)],na.rm=T)/class_dat$sum_tot

# Subset proportion by species class
class_dat  <- class_dat[,57:61]

# Class into state types

class_fn = function(x) { 
  classPlot  <- NULL
  if(min(complete.cases(x))==1) {
  if(sum(x[1:3]) < 1/3) {classPlot="Unclass"} 	
  else if(x[1] > 2/3) {classPlot="R"}
  else if(x[2] > 2/3) {classPlot="B"}
  else if(x[3] > 2/3) {classPlot="T"}
  else {classPlot="M"}
  } else {classPlot="Unclass"}
return(classPlot)
}

class_final  <- as.vector(apply(class_dat[,3:5],1,class_fn))

# Rename columns
class_dat  <- cbind(class_dat,class_final)
names(class_dat)[c(1,6)]  <- c("obs","pred")

# Format factors levels
class_dat$obs  <- as.factor(class_dat$obs) 
class_dat$obs  <- factor(class_dat$obs,levels=c(levels(class_dat$obs)[3],levels(class_dat$obs)[1],levels(class_dat$obs)[2],"R","Unclass"),labels=c("B","T","M","R","Unclass"))

###################################################################
#####               Reshape and export data                 #######
###################################################################

## Function to match remeasurements on the same line 
prwise  <- function(x,clim=FALSE){
 
  if (clim == FALSE){
    if (class(x)=='numeric') df  <-  data.frame(col1=numeric(length(x)-1),col2=numeric(length(x)-1))
    if (class(x)=='character') df  <-  data.frame(col1=character(length(x)-1),col2=character(length(x)-1),stringsAsFactors=FALSE)
    colnames(df)  <- c("t0","t1")
    for(i in 1:dim(df)[1]){
      df[i,1]  <-  x[i]
      df[i,2]  <- x[i+1]
    }
  }
  
  if (clim == TRUE){
    df  <- numeric(length(x)-1)
    for(i in 2:length(x)){
      df[i-1]  <-  mean(c(x[i-1],x[i]))
    }
  }
  return(df)
}

# Subset only columns needed and remove unclass plots
res_dat  <- cbind(dat[,c(2:6,68:71)],pred_class=class_dat$pred)
res_dat$pred_class  <- as.character(res_dat$pred_class,stringsAsFactors=FALSE)

write.table(res_dat,file="../data/data_allyears_RBTM.txt")


## Split by id_plot
res_dat  <- split(res_dat,res_dat$id_plot)

## Pair
pair <- function(x) { res  <- cbind(id_plot=x[-1,"id_plot"],
                                av_annual_pp=prwise(x[,"annual_pp"],clim=TRUE),
                                av_annual_mean_tp=prwise(x[,"annual_mean_temp"],clim=TRUE),
                                av_annual_min_tp=prwise(x[,"annual_min_temp"],clim=TRUE),
                                av_annual_max_tp=prwise(x[,"annual_max_temp"],clim=TRUE),
                                int=diff(x[,"yr_measured"]),
                                prwise(x[,"pred_class"],clim=FALSE))
                      return(res)
                     }

## Transpose states and compute mean
reshape_dat <- lapply(res_dat,pair)

## Final reshaping and export
reshape_dat  <- do.call(rbind,reshape_dat)
reshape_dat  <- data.frame(reshape_dat,row.names=NULL)
table(reshape_dat$t0)

reshape_dat = reshape_dat[reshape_dat$t0 != "Unclass" & reshape_dat$t1 != "Unclass",]

head(reshape_dat)

# ------------------- write

write.table(reshape_dat,file="../data/data_reshaped_RBTM.txt")

# ---------------- herbivores
reshape_dat = read.table("../data/data_reshaped_RBTM.txt")
head(reshape_dat)
dim(reshape_dat)

moose = read.table("../data/densities_moose.txt", h=T, sep="\t")
colnames(moose) = c("NO_ZONE"     ,  "id_plot", "moose")
head(moose)
dim(moose)

deer = read.table("../data/densities_deer.txt", h=T, sep="\t")
colnames(deer) = c("NO_ZONE"     ,  "id_plot", "deer")
head(deer)
dim(deer)

herbivores = merge(moose[,-1], deer[,-1], by = 'id_plot')
dim(herbivores)
head(herbivores)
herbivores[is.na(herbivores$moose),]
herbivores[is.na(herbivores$deer),]

finalTab = merge(reshape_dat, herbivores, by = "id_plot")
head(finalTab)
dim(finalTab)

library(foreign)
coords = read.csv("../data/plot_coords.csv")
head(coords)
coords$lat[which(coords$lat==0.0)]=NA
coords = na.omit(coords)
summary(coords)

fint = merge(finalTab, coords, by = 'id_plot', all.x=TRUE, all.y = FALSE)
head(unique(fint))

write.table(fint,file="../data/data_reshaped_RBTM_herbivores.txt")


