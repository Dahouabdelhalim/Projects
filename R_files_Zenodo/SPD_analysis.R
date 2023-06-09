
#Rcode for the SPD and sedimentary charcoal time series analysed in the manuscript:

#Sánchez-García, C., Revelles, J., Burjachs, F., Euba, I., Expósito, I., Ibáñez, J., 
#Schulte, L., Fernández-López de Pablo, J. What burned the forest? Wildfires, climate 
#change and human activity during the Mesolithic – Neolithic transition in SE 
#Iberian Peninsula (submitted to Catena). 

#Author: Javier Fernández-López de Pablo (https://orcid.org/0000-0002-6953-7004)
#Funding: European Comission-H2020-ERC (Ref.683018)

library(readr)
library(rcarbon)
#Loading the radiocarbon data
Vinalopo <- read.csv("~/Documents/Rdata/Charcoal/Vinalopo.csv")
View(Vinalopo)

#Loading the microcharcoal data set
Char <- read.csv("~/Documents/Rdata/Charcoal/Char.csv")

#Loading and plotting the ngrip data set
ngrip<-read.csv("~/Documents/Rdata/Charcoal/ngrip.csv")
ngrip
plot(-1*ngrip$Age_b2k, ngrip$d18O , type = 'l', xlim = c(-9500, -5500))

#Building simple SPDs
Vinalopo.caldates=calibrate(x=Vinalopo$C14Age, errors = Vinalopo$C14SD, calCurves = 'intcal20')
Vinalopo.spd=spd(Vinalopo.caldates, timeRange = c(9500, 5500))
plot(Vinalopo.spd)
plot(Vinalopo.spd,runm = 100, add = TRUE, type = "simple", col="red", lwd=1.5, lty=2)

#Binning process and binned SPDs

Vinalopo.bins=binPrep(sites = Vinalopo$SiteID, ages = Vinalopo$C14Age, h=100)
Vinalopo.spd.bins=spd(Vinalopo.caldates, bins=Vinalopo.bins, timeRange = c(9500, 5500))
plot(Vinalopo.spd.bins)
plot(Vinalopo.spd.bins, runm = 100, add = TRUE, type = "simple", col="red", lwd=1.5, lty=2)

#Sensivity analysis of bins

binsense(x=Vinalopo.caldates, y=Vinalopo$SiteID, h=seq(0, 500, 100), timeRange = c(10000, 5500))
Vinalopo.bins.med=binMed(x=Vinalopo.caldates, bins=Vinalopo.bins)
plot(Vinalopo.spd.bins, runm = 100)
barCodes(Vinalopo.bins.med, yrng=c(0,0.002))
### Building a multi-plot with the NGrip, the microcarcoal and the radiocarbon record
par(mfrow=c(3,1))
par(mar=c(3,3,0.3,0.3))
xl<-seq(-10000,-5500,500)

plot(-1*ngrip$Age_b2k, ngrip$d18O , type = 'l', col='blue', xlim = c(-10000, -5500), ylim = c(-36, -33), las=2, frame.plot=FALSE,xaxt='n')
plot(Char$Age, Char$Char, type='l', col='black', xlim = c(-10000, -5500), ylim = c(0, 3), las=2, frame.plot=FALSE,xaxt='n')
plot(Vinalopo.spd.bins, runm = 100, add = TRUE, type = "simple", col="red", lwd=1.5, lty=2)





