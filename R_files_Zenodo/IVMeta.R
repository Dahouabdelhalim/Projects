setwd("H:/BigLake/Paper/Regressions")

library(doBy)
library(tidyr)
library(dplyr)
library(metafor)

secchi <- read.csv("SecchiBootstrap.csv",stringsAsFactors = F)

result <- rma(secchi, secchise, data=secchi)
summary(result)

resultunweight <- rma(secchi, secchise, weighted=FALSE, data=secchi)
summary(resultunweight)



result <- rma(maxdepth, maxdepthse, data=secchi)
summary(result)

resultunweight <- rma(maxdepth, maxdepthse, weighted=FALSE, data=secchi)
summary(resultunweight)



result <- rma(boatlaunch, boatlaunchse, data=secchi)
summary(result)

resultunweight <- rma(boatlaunch, boatlaunchse, weighted=FALSE, data=secchi)
summary(resultunweight)



result <- rma(beach, beachse, data=secchi)
summary(result)

resultunweight <- rma(beach, beachse, weighted=FALSE, data=secchi)
summary(resultunweight)



result <- rma(hotels, hotelsse, data=secchi)
summary(result)

resultunweight <- rma(hotels, hotelsse, weighted=FALSE, data=secchi)
summary(resultunweight)




result <- rma(shelter, shelterse, data=secchi)
summary(result)

resultunweight <- rma(shelter, shelterse, weighted=FALSE, data=secchi)
summary(resultunweight)



result <- rma(toilets, toiletsse, data=secchi)
summary(result)

resultunweight <- rma(toilets, toiletsse, weighted=FALSE, data=secchi)
summary(resultunweight)




result <- rma(picnic, picnicse, data=secchi)
summary(result)

resultunweight <- rma(picnic, picnicse, weighted=FALSE, data=secchi)
summary(resultunweight)




result <- rma(bbq, bbqse, data=secchi)
summary(result)

resultunweight <- rma(bbq, bbqse, weighted=FALSE, data=secchi)
summary(resultunweight)




result <- rma(marina, marinase, data=secchi)
summary(result)

resultunweight <- rma(marina, marinase, weighted=FALSE, data=secchi)
summary(resultunweight)



result <- rma(Developed, Developedse, data=secchi)
summary(result)

resultunweight <- rma(Developed, Developedse, weighted=FALSE, data=secchi)
summary(resultunweight)



result <- rma(wsecchi, wsecchise, data=secchi)
summary(result)

resultunweight <- rma(wsecchi, wsecchise, weighted=FALSE, data=secchi)
summary(resultunweight)




result <- rma(d5lakes, d5lakesse, data=secchi)
summary(result)

resultunweight <- rma(d5lakes, d5lakesse, weighted=FALSE, data=secchi)
summary(resultunweight)



result <- rma(bach, bachse, data=secchi)
summary(result)

resultunweight <- rma(bach, bachse, weighted=FALSE, data=secchi)
summary(resultunweight)




result <- rma(medinc, medincse, data=secchi)
summary(result)

resultunweight <- rma(medinc, medincse, weighted=FALSE, data=secchi)
summary(resultunweight)





result <- rma(lagy, lagyse, data=secchi)
summary(result)

resultunweight <- rma(lagy, lagyse, weighted=FALSE, data=secchi)
summary(resultunweight)





result <- rma(lagsecchi, lagsecchise, data=secchi)
summary(result)

resultunweight <- rma(lagsecchi, lagsecchise, weighted=FALSE, data=secchi)
summary(resultunweight)




result <- rma(distcbsa, distcbsase, data=secchi)
summary(result)

resultunweight <- rma(distcbsa, distcbsase, weighted=FALSE, data=secchi)
summary(resultunweight)




result <- rma(popden, popdense, data=secchi)
summary(result)

resultunweight <- rma(popden, popdense, weighted=FALSE, data=secchi)
summary(resultunweight)

