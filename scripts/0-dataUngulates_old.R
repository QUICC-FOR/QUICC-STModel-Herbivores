path = "~/Documents/RESEARCH/DATA_UQAR/SIG/"

library(sp)
library(maptools) 
library(rgdal)
library(RColorBrewer) # creates nice color schemes
library(classInt) # finds class intervals for continuous variables


## ----- read shapefile

wd = getwd()
setwd(path)
decoup = readOGR(dsn=".", layer="regions_admin_Qc")
setwd(wd)

plot(decoup)


## ----- add densities data


#locator(1) you choose the place using the pointer
#nclr <- 22
#plotclr <- brewer.pal(nclr,"BuPu")
#class <- classIntervals(plotvar, nclr, style="quantile",dataPrecision=0)
#colcode <- findColours(class, plotclr)
#title(main="your title", sub="Quantile (Equal-Frequency) Class Intervals")
#legend(locator(1), legend=names(attr(colcode, "table")), fill=attr(colcode, "palette"), cex=0.6, bty="n")

plotvar <- decoup@data$REGIO_S_ID
#plotclr = colors()[c(1:22)*4]
plotclr = c(1, 2, 2, 2, 2, 3, 4, rainbow(16))

plot(decoup)
plot(decoup, col=plotclr, add=T)
legend("topright", legend=decoup@data$TOPONYME, fill=plotclr, cex=0.6, bty="n")
legend("topleft", legend=decoup@data$CODE, fill=plotclr, cex=0.6, bty="n")

decoup@data[,-c(1:2)]

orignal = list('10' = c(16,17,22), '09' = c(18,19,23), '02' = c(28,29), 
'11' = c(1, 20, 21), 

## ----- points

library(foreign)
coords = read.csv("../data/plot_coords.csv")
head(coords)
coords$lat[which(coords$lat==0.0)]=NA
coords = na.omit(coords)

library(sp)
pts = SpatialPointsDataFrame(coords[,2:3], coords, proj4string=CRS("+proj=longlat +datum=WGS84"))
ptsProj <- spTransform(pts, CRS(proj4string(decoup)))

points(ptsProj, pch = ".")
