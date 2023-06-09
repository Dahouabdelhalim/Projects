#setwd("...")

#1. The modelling of diagnostic sensitivity (dSE)
#modelling_dSe.R

source("Functions_modelling_dSe.R")
require(fitdistrplus)

###############################################################################
#1) Assuming time since infection are randomly drawn from time line of 2 years
############################################################################### 

#Simulated distribution of dSe
dSE_obex<-SeAdults[,"dSE_obex"]   #dSe when testing obex of adults
dSE_obexY<-SeYearlings[,"dSE_obex"]  #dSe when testing obex of yearlings
#dSE_obexC<-SeCalves[,"dSE_obex"]

dSE_RLNobex<-SeAdults[,"dSE_RLN"]    #dSe when testing RLN&obex of adults
dSE_RLNobexY<-SeYearlings[,"dSE_RLN"]   #dSe when testing RLN&obex of yearlings
dSE_RLNobexC<-SeCalves[,"dSE_RLN"]   #dSe when testing RLN&obex of calves

nsim=length(dSE_obex)     #10000

#Proportion of infected being non-detectable
Prob0_obex<-round(table(dSE_obex==0)[[2]]/nsim,2)
Prob0_obexY<-round(table(dSE_obexY==0)[[2]]/nsim,2)
Prob0_obexC<-round(table(SeCalves[,"dSE_obex"]==0)/nsim,2)[[1]]  
Prob0_RLNobex<-round(table(dSE_RLNobex==0)[[2]]/nsim,2)
Prob0_RLNobexY<-round(table(dSE_RLNobexY==0)[[2]]/nsim,2)
Prob0_RLNobexC<-round(table(dSE_RLNobexC==0)[[2]]/nsim,2)

#Proportion of infected with dSe_RLNobex>0 and dSe_obex=0
Pr_obex0.RLNposC<-round(table(SeCalves[SeCalves[,"dSE_RLN"]>0,"dSE_obex"]==0)/nsim,3)[[1]]  
Pr_obex0.RLNposY<-round(table(SeYearlings[SeYearlings[,"dSE_RLN"]>0,"dSE_obex"]==0)[[2]]/nsim,3)  
Pr_obex0.RLNpos<-round(table(SeAdults[SeAdults[,"dSE_RLN"]>0,"dSE_obex"]==0)[[2]]/nsim,3)   

#Fitting non-zero part of the simulated dSe by a beta distribution
fitbeta_dRLNobex <- fitdist(dSE_RLNobex[dSE_RLNobex>0], "beta")
fitbeta_dRLNobexY <- fitdist(dSE_RLNobexY[dSE_RLNobexY>0], "beta")
fitbeta_dRLNobexC <- fitdist(dSE_RLNobexC[dSE_RLNobexC>0], "beta")

fitbeta_dObex <- fitdist(dSE_obex[dSE_obex>0], "beta")
fitbeta_dObexY <- fitdist(dSE_obexY[dSE_obexY>0], "beta")

dSE_betadist<-rbind(c(fitbeta_dRLNobexC$estimate[1],fitbeta_dRLNobexC$estimate[2],Prob0_RLNobexC,NA),
c(fitbeta_dRLNobexY$estimate[1],fitbeta_dRLNobexY$estimate[2],Prob0_RLNobexY,NA),
c(fitbeta_dRLNobex$estimate[1],fitbeta_dRLNobex$estimate[2],Prob0_RLNobex,NA),
c(NA, NA,Prob0_obexC,Pr_obex0.RLNposC),
c(fitbeta_dObexY$estimate[1], fitbeta_dObexY$estimate[2],Prob0_obexY,Pr_obex0.RLNposY),
c(fitbeta_dObex$estimate[1],fitbeta_dObex$estimate[2],Prob0_obex,Pr_obex0.RLNpos))

rownames(dSE_betadist)<-c("dRLNobexC","dRLNobexY","dRLNobex",
"dObexC","dObexY","dObex")
colnames(dSE_betadist)<-c("shape1","shape2","PrZeros","Pr_Obex0.RLNpos")
dSE_betadist<-round(dSE_betadist,2)
dSE_betadist

###############################################################################
#2) Assuming time since infection are randomly drawn from time line of 3 years
###############################################################################
#Time step for drawing random time since infection, assuming lenght of disease progression is 3 years

#Simulated distribution of dSe
dSE_obex_longI3<-Se3yrI_Adults[,"dSE_obex"]   #dSe when testing obex of adults
dSE_obexY_longI3<-Se3yrI_Yearlings[,"dSE_obex"] #dSe when testing obex of yearlings
dSE_obexC_longI3<-Se3yrI_Calves[,"dSE_obex"]

dSE_RLNobex_longI3<-Se3yrI_Adults[,"dSE_RLN"] #dSe when testing RLN&obex of adults
dSE_RLNobexY_longI3<-Se3yrI_Yearlings[,"dSE_RLN"] #dSe when testing RLN&obex of yearlings
dSE_RLNobexC_longI3<-Se3yrI_Calves[,"dSE_RLN"]

nsim=length(dSE_obex_longI3)

#Proportion of infected being non-detectable
Prob0_RLNobexC_longI3<-round(table(dSE_RLNobexC_longI3==0)[[1]]/nsim,3)
Prob0_RLNobexY_longI3<-round(table(dSE_RLNobexY_longI3==0)[[2]]/nsim,3)
Prob0_RLNobex_longI3<-round(table(dSE_RLNobex_longI3==0)[[2]]/nsim,3)

Prob0_obexC_longI3<-round(table(dSE_obexC_longI3==0)[[1]]/nsim,3)
Prob0_obexY_longI3<-round(table(dSE_obexY_longI3==0)[[2]]/nsim,3)
Prob0_obex_longI3<-round(table(dSE_obex_longI3==0)[[2]]/nsim,3)

#Proportion of infected with dSe_RLNobex>0 and dSe_obex=0
 
Pr_obex0.RLNposY_longI3<-round(table(Se3yrI_Yearlings[Se3yrI_Yearlings[,"dSE_RLN"]>0,"dSE_obex"]==0)[[2]]/nsim,3)  
Pr_obex0.RLNpos_longI3<-round(table(Se3yrI_Adults[Se3yrI_Adults[,"dSE_RLN"]>0,"dSE_obex"]==0)[[2]]/nsim,3)   


#Fitting non-zero part of the simulated dSe by a beta distribution
fitbeta_dRLNobex_longI3 <- fitdist(dSE_RLNobex_longI3[dSE_RLNobex_longI3>0], "beta")
fitbeta_dRLNobexY_longI3 <- fitdist(dSE_RLNobexY_longI3[dSE_RLNobexY_longI3>0], "beta")
fitbeta_dObex_longI3 <- fitdist(dSE_obex_longI3[dSE_obex_longI3>0], "beta")
fitbeta_dObexY_longI3 <- fitdist(dSE_obexY_longI3[dSE_obexY_longI3>0], "beta")

dSE_betadist_longI3<-rbind(c(NA,NA,1,NA),
c(fitbeta_dRLNobexY_longI3$estimate[1],fitbeta_dRLNobexY_longI3$estimate[2],Prob0_RLNobexY_longI3,NA),
c(fitbeta_dRLNobex_longI3$estimate[1],fitbeta_dRLNobex_longI3$estimate[2],Prob0_RLNobex_longI3,NA),
c(NA, NA,1,0),
c(fitbeta_dObexY_longI3$estimate[1], fitbeta_dObexY_longI3$estimate[2],Prob0_obexY_longI3,Pr_obex0.RLNposY_longI3),
c(fitbeta_dObex_longI3$estimate[1],fitbeta_dObex_longI3$estimate[2],Prob0_obex_longI3,Pr_obex0.RLNpos_longI3))


rownames(dSE_betadist_longI3)<-c("dRLNobexC","dRLNobexY","dRLNobex",
"dObexC","dObexY","dObex")
colnames(dSE_betadist_longI3)<-c("shape1","shape2","PrZeros","Pr_Obex0.RLNpos")
dSE_betadist_longI3<-round(dSE_betadist_longI3,2)
dSE_betadist_longI3

##############################################################################
#3) Assuming an exponential time decay distribution of time since infection 
##############################################################################

#Assume mean time since infection = 1 year
#We are deleting any time since infection > 4 years, hence mean time since infection corresponding to the simulated values will be slightly lower

dSEexp_obex<-SeExpAdults[,"dSE_obex"]
dSEexp_obexY<-SeExpYearlings[,"dSE_obex"]
dSEexp_obexC<-SeExpCalves[,"dSE_obex"]

dSEexp_RLNobex<-SeExpAdults[,"dSE_RLN"]
dSEexp_RLNobexY<-SeExpYearlings[,"dSE_RLN"]
dSEexp_RLNobexC<-SeExpCalves[,"dSE_RLN"]

nsim=length(dSEexp_obex)

#Proportion of infected being non-detectable
Prob0Exp_RLNobexC<-round(table(dSEexp_RLNobexC==0)[[2]]/nsim,2)
Prob0Exp_RLNobexY<-round(table(dSEexp_RLNobexY==0)[[2]]/nsim,2)
Prob0Exp_RLNobex<-round(table(dSEexp_RLNobex==0)[[2]]/nsim,2)

Prob0Exp_obexC<-round(table(dSEexp_obexC==0)[[1]]/nsim,2) 
Prob0Exp_obexY<-round(table(dSEexp_obexY==0)[[2]]/nsim,2)
Prob0Exp_obex<-round(table(dSEexp_obex==0)[[2]]/nsim,2)

#Proportion of infected with dSe_RLNobex>0 and dSe_obex=0
PrExp_obex0.RLNposC<-round(table(SeExpCalves[SeExpCalves[,"dSE_RLN"]>0,"dSE_obex"]==0)/nsim,2)[[1]] 
PrExp_obex0.RLNposY<-round(table(SeExpYearlings[SeExpYearlings[,"dSE_RLN"]>0,"dSE_obex"]==0)[[2]]/nsim,2)  
PrExp_obex0.RLNpos<-round(table(SeExpAdults[SeExpAdults[,"dSE_RLN"]>0,"dSE_obex"]==0)[[2]]/nsim,2)   


#Fitting non-zero part of the simulated dSe by a beta distribution
fitbeta_dObexExp <- fitdist(dSEexp_obex[dSEexp_obex>0], "beta")
fitbeta_dRLNobexExp <- fitdist(dSEexp_RLNobex[dSEexp_RLNobex>0], "beta")
fitbeta_dObexYExp <- fitdist(dSEexp_obexY[dSEexp_obexY>0], "beta")
fitbeta_dRLNobexYExp <- fitdist(dSEexp_RLNobexY[dSEexp_RLNobexY>0], "beta")
fitbeta_dRLNobexCExp <- fitdist(dSEexp_RLNobexC[dSEexp_RLNobexC>0], "beta")

dSE_betadistExp24<-rbind(c(fitbeta_dRLNobexCExp$estimate[1],fitbeta_dRLNobexCExp$estimate[2],Prob0Exp_RLNobexC,NA),
c(fitbeta_dRLNobexYExp$estimate[1],fitbeta_dRLNobexYExp$estimate[2],Prob0Exp_RLNobexY,NA),
c(fitbeta_dRLNobexExp$estimate[1],fitbeta_dRLNobexExp$estimate[2],Prob0Exp_RLNobex,NA),c(NA, NA,Prob0Exp_obexC,PrExp_obex0.RLNposC),
c(fitbeta_dObexYExp$estimate[1], fitbeta_dObexYExp$estimate[2],Prob0Exp_obexY,PrExp_obex0.RLNposY),
c(fitbeta_dObexExp$estimate[1],fitbeta_dObexExp$estimate[2],Prob0Exp_obex,PrExp_obex0.RLNpos))

rownames(dSE_betadistExp24)<-c("dRLNobexC","dRLNobexY","dRLNobex",
"dObexC","dObexY","dObex")
colnames(dSE_betadistExp24)<-c("shape1","shape2","PrZeros","Pr_Obex0.RLNpos")
dSE_betadistExp24<-round(dSE_betadistExp24,2)
dSE_betadistExp24

