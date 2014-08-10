
library(sp)
library(maptools) 
library(rgdal)
library(RColorBrewer) # creates nice color schemes
library(classInt) # finds class intervals for continuous variables


## ----- read shapefile

wd = getwd()
setwd("../data/Zones_chasse/")
zones = readOGR(dsn=".", layer="Zone_chasse_2014")
setwd(wd)

##-----

head(zones@data)

## tableaux

tableMoose = data.frame(
zone = c('01','02','03','04','05','06','07','08','09','10','11 E','11 W','12','13','14','15','16','17','18','19 S','22','26','27','28','29', 'other'), 
densites_10ha = c(7.90, 6.80,5.95 ,1.48, 0, 2.41, 2.70, 0, 1.08, 2.40, 0.98, 2.39, 
3.20, 3.08, 1.80, 1.70, 1.11, 0.45, 1.26, 0.44, 0.26, 2.32,3.20,0.54,0.44, 0)
)


tableDeer = data.frame(
zone = c('01', '02', '02 W', '03', '03 W', '04', '05', '05 W', '06 N', '06', 
'07', '07 S', '08 E', '08', '08 S', '09', '09 W', '10', '10 W', '11 W', '12', 
'13', '26 E', '27', '28', '15 W', 'other'), 
densites_10ha = c(0.37, 0.25, 0.54, 0.70, 2.37, 4.54,6.17, 9.10, 8.71,3.94, 3.86, 9.07, 8.20, 5.58, 7.53, 1.78, 2.01,  2.57, 3.12, 2.22, 0, 0.39, 0.9, 1.59, 
0, 2.22, 0)
)


## ----- points

library(foreign)
coords = read.csv("../data/plot_coords.csv")
head(coords)
coords$lat[which(coords$lat==0.0)]=NA
coords = na.omit(coords)

library(sp)
pts = SpatialPointsDataFrame(coords[,2:3], coords, proj4string=CRS("+proj=longlat +datum=WGS84"))
ptsProj <- spTransform(pts, CRS(proj4string(zones)))

## --- croisement

newPts = over(ptsProj, zones[, "NO_ZONE"])
table(newPts)


zonesMoose = data.frame(pts, newPts)
levels(zonesMoose$NO_ZONE) = c("01" ,   "02"  ,  "02",  "03" ,   "03" , "04" ,   "04"  , "05" ,   "05" ,"06" ,   "06" , "07"  ,  "07" , "08" ,   "08" , "08" , "09" ,   "09" ,"10"  ,  "10"  ,"10" ,"10" , "11"  ,  "11 E" , "11 W" , "12" ,   "13" ,  "13" ,"14" ,   "15" ,   "15" , "16" ,   "17" ,   "18" ,   "other" ,   "19 S" , "other", "other" ,"other" ,   "other" ,   "other" ,   "22"  ,  "22" , "22" , "other" ,   "other" , "other" , "other" ,   "26"  ,  "26" , "27" ,   "27" , "27" , "28" ,  "29")

densMoose = merge(zonesMoose, tableMoose, by.x = 'NO_ZONE', by.y = 'zone', all.x=TRUE)
dim(pts)
dim(densMoose)
head(densMoose)

densMoose = unique(densMoose[, c("NO_ZONE", "id_plot", "densites_1000ha")])
densMoose$densites_10ha[is.na(densMoose$densites_10ha)] = 0
write.table(densMoose, "../data/densities_moose.txt", quote=F, row.names=F, sep="\t")


zonesDeer = data.frame(pts, newPts)
levels(zonesDeer$NO_ZONE) = c("01", "02" ,    "02 W" ,  "03"  ,   "03 W" ,  "04" ,   "04" ,   "05"  ,   "05 W",  "06"   ,  "06 N" ,  "07"  ,   "07 S" ,  "08" ,    "08 E" ,  "08 S"  , "09" ,    "09 W",  "other"  ,   "10" , "other" , "10 W"  , "other",     "other" , "11 W" ,  "12" ,    "13" ,   "13" , "other"  ,   "other" ,    "15 W"  , "other"  ,   "other"  ,   "other"  ,   "other"  ,   "other" , 
"other" , "other" , "other" ,    "other"  ,  "other"  ,   "other" ,    "other" ,  "other" ,  "other"  ,  "other" ,  "other" ,  "other"  ,   "other"  ,   "26 E",   "27"  ,   "27"  ,"27" ,  "28",     "other")

densDeer = merge(zonesDeer, tableDeer, by.x = 'NO_ZONE', by.y = 'zone', all.x=TRUE)
dim(pts)
dim(densDeer)
head(densDeer)

densDeer = unique(densDeer[, c("NO_ZONE", "id_plot", "densites_100ha")])
densDeer$densites_10ha[is.na(densDeer$densites_10ha)] = 0
write.table(densDeer, "../data/densities_deer.txt", quote=F, row.names=F, sep="\t")

#plot(zones)
#points(ptsProj, pch = ".")
