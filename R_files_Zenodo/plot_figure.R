library(sp)
library(spdep)
library(INLA)
library(raster)
library(ggplot2)
library(SpatialEpi)
library(maptools)
library(spatialEco)
library(rgdal)
library(reshape2)
library(plyr)
library(ROCR)
library(OptimalCutpoints)
library(viridis)
library(ggplot2)
library(ggthemes)
library(rasterVis)

#requires the file "parks.Rdata" in the same folder in additon to the shapefiles below
#take the previous modeling data
load("parks.Rdata")

#rasterize the 2017 model
mod2017<-rasterize(dnare, rast, "y17", fun = median)

#import the raster with forest cover data for 2018
def18 <- raster("defdeg20181.tif")

#change extent
extent(def18) <- extent(parks)

#replace values
values(def18) <- ifelse(values(def18) ==3, 2, ifelse(values(def18) ==1,NA, ifelse(values(def18) ==2, 1,0))) 

#import the raster with forest cover data for 2017
loss17 <- raster("3PA_Hansen_v1.5_loss2017.tif")

#eliminate the issue of high values introducing noise
loss17[loss17 == 128] <- NA

#eliminate previously converted pixels
loss17[loss17 == 1] <- NA
loss17[loss17 == 17] <- 1

#print the figure
#tiff(file = "figure2.tiff", width = 4800, height = 3200, units = "px", res = 300)
#par(mfrow=c(2,2))
#plot(mod2017, col=rev(inferno(256)), main="Predicted 2017")
#plot(loss17, breaks = c(0,0.5,1), col=rev(viridis(3)), main="Observed 2017", #axis.args=list(at=c(0,1,0), labels=c("Forest","Loss","")), horizontal=T)
#plot(delta1718, col=rev(inferno(256)), main="Predicted change in probability 2017-2018")
#plot(def18, breaks = c(0,0.99,1.99,3), col=rev(viridis(3)), main="Observed 2018", #axis.args=list(at=c(0,1.5,0,2.99), labels=c("Forest","Degradation","","Loss")), #horizontal=T)
#dev.off()

#print the figure as pdf
pdf(file = "figure2.pdf", width = 12, height = 8)
par(mfrow=c(2,2))
plot(mod2017, col=rev(inferno(256)), main="Predicted 2017")
plot(loss17, breaks = c(0,0.5,1), col=rev(viridis(3)), main="Observed 2017", axis.args=list(at=c(0,1,0), labels=c("Forest","Loss","")), horizontal=T)
plot(delta1718, col=rev(inferno(256)), main="Predicted change in probability 2017-2018")
plot(def18, breaks = c(0,0.99,1.99,3), col=rev(viridis(3)), main="Observed 2018", axis.args=list(at=c(0,1.5,0,2.99), labels=c("Forest","Degradation","","Loss")), horizontal=T)
dev.off()

