#Upload library
library(lattice)
library(FME)
library(deSolve)
library(growthrates)

#Upload data or import
R.growthrate = read.table("RG", header=TRUE)

#Run model
rrs=all_splines(value ~ time|strain, data=R.growthrate, spar=0.5)
par(mfrow= c(2,4))
plot(rrs)
coef(rrs)

#Results
R.growthrate.summary=summary(rrs)
summary(rrs)
