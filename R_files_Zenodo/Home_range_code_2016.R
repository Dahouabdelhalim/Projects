###########Home range analyses for Vole 2016 data##################
##goals:  1.  get HR values using Kernel Density Utilization Distributions
#2.  determine HR kernel overlap for voles in same enclosure
#3.  export shapefiles of HR for use in ArcGIS
rm(list=ls()) 

library(adehabitatHR)
library(sp)
library(plyr)
library(reshape2)
library(maptools)

voles0 = read.csv("GPS data_edit.csv")
###remove voles w/less than 4 captures/period
captures <- ddply(voles0,.(HR.Period,
                          Enclosure,
                          Vole.ID), 
                  summarise, Tot=length(Vole.ID))
head(captures)

mergecap <- merge(voles0, captures, 
                  by = c("HR.Period", "Enclosure", "Vole.ID"),
                  all.x= TRUE); head(mergecap)
voles = subset(mergecap, subset = Tot >4); head(voles)

#########Subsetting Data###########  

voles2 = voles[voles$Encl=="A",]#subset to one enclosure- do one at a time

df = subset(voles2, subset = HR.Period == 1)
  #OR
df = subset(voles2, subset = HR.Period == 2)
  #OR
df = subset(voles2, subset = HR.Period == 3)


#########create Spatial data frame to use for HR analyses#######
df = droplevels(df)
attributes <- data.frame(df$HR.Period, df$Enclosure, df$Vole.ID ) #attribute IDs each animal
xy <- cbind(df$Easting,df$Northing)
voleloc <- SpatialPointsDataFrame(xy, attributes)

###############Calculating  Kernels## --------------------------
kud <- kernelUD(voleloc[,3], h = 2, grid = 100, extent = 2, same4all= TRUE)

KHR = kernel.area(kud, percent=75, unin = "m",unout="m2")####give actual areas at 50%
melt(KHR) ##make it df long format
df <- ldply (KHR, data.frame) 
df = rename(df, c(".id"="Vole.ID", "X..i.."="HR"))

######Kernel overlap#######
Over = kerneloverlap(voleloc[,3], 
                      grid=100, h = 2, extent = 2, 
                      percent = 75,
                      meth="HR") #method for estimation of overlap
Over
Over2 = melt(Over)
######Shapefile of HR#####
shape <- getverticeshr(kud, percent = 75)
plot(shape)
writePolyShape(shape, "P3")

df
Over2

