setwd("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Cornell/bee_pathogens/2015-FINAL/science/verified")

wf<-read.csv("visitation_2015UR.csv")
pathogens<-read.csv("pathogens_2015UR.csv")

library(vegan)
library(bipartite)
library(plyr)
library(car)

#exclude date and flower species information to fit bipartite requirements
net<- wf[c(-1,-2)]

#most simple network: all interactions, at all sites, throughout summer
WFsum <- aggregate(. ~ Species, data=net, FUN=sum)
row.names(WFsum) <- WFsum[,1]
WFsum <-WFsum[,c(-1)]

# separate by site 
CB<- subset(wf, site=="CB")
CB<- CB[c(-1,-2)]
CO<-subset(wf, site=="CO")
CO<- CO[c(-1,-2)]
CP<-subset(wf, site=="CP")
CP<- CP[c(-1,-2)]
R4<-subset(wf, site=="R4")
R4<- R4[c(-1,-2)]
RJ<- subset(wf, site=="RJ")
RJ<- RJ[c(-1,-2)]
RL<-subset(wf, site=="RL")
RL<- RL[c(-1,-2)]
DW<-subset(wf, site=="DW")
DW<- DW[c(-1,-2)]
EI<- subset(wf, site=="EI")
EI<- EI[c(-1,-2)]
KP<-subset(wf, site=="KP")
KP<- KP[c(-1,-2)]
FC<-subset(wf, site=="FC")
FC<- FC[c(-1,-2)]
FO<- subset(wf, site=="FO")
FO<- FO[c(-1,-2)]

CB <- aggregate(. ~ Species, data=CB, FUN=sum)
row.names(CB) <- CB[,1]
CB <-CB[,c(-1)]
CO <- aggregate(. ~ Species, data=CO, FUN=sum)
row.names(CO) <- CO[,1]
CO <-CO[,c(-1)]
CP <- aggregate(. ~ Species, data=CP, FUN=sum)
row.names(CP) <- CP[,1]
CP <-CP[,c(-1)]
R4 <- aggregate(. ~ Species, data=R4, FUN=sum)
row.names(R4) <- R4[,1]
R4 <-R4[,c(-1)]
RL <- aggregate(. ~ Species, data=RL, FUN=sum)
row.names(RL) <- RL[,1]
RL <-RL[,c(-1)]
RJ <- aggregate(. ~ Species, data=RJ, FUN=sum)
row.names(RJ) <- RJ[,1]
RJ <-RJ[,c(-1)]
DW <- aggregate(. ~ Species, data=DW, FUN=sum)
row.names(DW) <- DW[,1]
DW <-DW[,c(-1)]
EI <- aggregate(. ~ Species, data=EI, FUN=sum)
row.names(EI) <- EI[,1]
EI <-EI[,c(-1)]
KP <- aggregate(. ~ Species, data=KP, FUN=sum)
row.names(KP) <- KP[,1]
KP <-KP[,c(-1)]
FC <- aggregate(. ~ Species, data=FC, FUN=sum)
row.names(FC) <- FC[,1]
FC <-FC[,c(-1)]
FO <- aggregate(. ~ Species, data=FO, FUN=sum)
row.names(FO) <- FO[,1]
FO <-FO[,c(-1)]

##################### H2 z-score ##########################
obs <- unlist(networklevel(CB, index="H2"))
nulls <- nullmodel(CB, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
CBH2Z<-as.data.frame(((obs-mean(null))/sd(null)))
CBH2Z$site<-"CB"
names(CBH2Z)<-c("H2Z","site")

obs <- unlist(networklevel(CO, index="H2"))
nulls <- nullmodel(CO, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
COH2Z<-as.data.frame(((obs-mean(null))/sd(null)))
COH2Z$site<-"CO"
names(COH2Z)<-c("H2Z","site")

obs <- unlist(networklevel(CP, index="H2"))
nulls <- nullmodel(CP, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
CPH2Z<-as.data.frame(((obs-mean(null))/sd(null)))
CPH2Z$site<-"CP"
names(CPH2Z)<-c("H2Z","site")

obs <- unlist(networklevel(DW, index="H2"))
nulls <- nullmodel(DW, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
DWH2Z<-as.data.frame(((obs-mean(null))/sd(null)))
DWH2Z$site<-"DW"
names(DWH2Z)<-c("H2Z","site")

obs <- unlist(networklevel(EI, index="H2"))
nulls <- nullmodel(EI, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
EIH2Z<-as.data.frame(((obs-mean(null))/sd(null)))
EIH2Z$site<-"EI"
names(EIH2Z)<-c("H2Z","site")

obs <- unlist(networklevel(FC, index="H2"))
nulls <- nullmodel(FC, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
FCH2Z<-as.data.frame(((obs-mean(null))/sd(null)))
FCH2Z$site<-"FC"
names(FCH2Z)<-c("H2Z","site")

obs <- unlist(networklevel(FO, index="H2"))
nulls <- nullmodel(FO, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
FOH2Z<-as.data.frame(((obs-mean(null))/sd(null)))
FOH2Z$site<-"FO"
names(FOH2Z)<-c("H2Z","site")

obs <- unlist(networklevel(KP, index="H2"))
nulls <- nullmodel(KP, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
KPH2Z<-as.data.frame(((obs-mean(null))/sd(null)))
KPH2Z$site<-"KP"
names(KPH2Z)<-c("H2Z","site")

obs <- unlist(networklevel(R4, index="H2"))
nulls <- nullmodel(R4, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
R4H2Z<-as.data.frame(((obs-mean(null))/sd(null)))
R4H2Z$site<-"R4"
names(R4H2Z)<-c("H2Z","site")

obs <- unlist(networklevel(RJ, index="H2"))
nulls <- nullmodel(RJ, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
RJH2Z<-as.data.frame(((obs-mean(null))/sd(null)))
RJH2Z$site<-"RJ"
names(RJH2Z)<-c("H2Z","site")

obs <- unlist(networklevel(RL, index="H2"))
nulls <- nullmodel(RL, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="H2")) 
RLH2Z<-as.data.frame(((obs-mean(null))/sd(null)))
RLH2Z$site<-"RL"
names(RLH2Z)<-c("H2Z","site")

H2<-rbind(CBH2Z,COH2Z,CPH2Z,DWH2Z,EIH2Z,FCH2Z,FOH2Z,KPH2Z,R4H2Z,RJH2Z,RLH2Z)

##################### Nestedness z-score ##########################
obs <- unlist(networklevel(CB, index="weighted NODF"))
nulls <- nullmodel(CB, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
CBWNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
CBWNODFZ$site<-"CB"
names(CBWNODFZ)<-c("WNODFZ","site")

obs <- unlist(networklevel(CO, index="weighted NODF"))
nulls <- nullmodel(CO, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
COWNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
COWNODFZ$site<-"CO"
names(COWNODFZ)<-c("WNODFZ","site")

obs <- unlist(networklevel(CP, index="weighted NODF"))
nulls <- nullmodel(CP, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
CPWNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
CPWNODFZ$site<-"CP"
names(CPWNODFZ)<-c("WNODFZ","site")

obs <- unlist(networklevel(DW, index="weighted NODF"))
nulls <- nullmodel(DW, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
DWWNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
DWWNODFZ$site<-"DW"
names(DWWNODFZ)<-c("WNODFZ","site")

obs <- unlist(networklevel(EI, index="weighted NODF"))
nulls <- nullmodel(EI, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
EIWNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
EIWNODFZ$site<-"EI"
names(EIWNODFZ)<-c("WNODFZ","site")

obs <- unlist(networklevel(FC, index="weighted NODF"))
nulls <- nullmodel(FC, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
FCWNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
FCWNODFZ$site<-"FC"
names(FCWNODFZ)<-c("WNODFZ","site")

obs <- unlist(networklevel(FO, index="weighted NODF"))
nulls <- nullmodel(FO, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
FOWNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
FOWNODFZ$site<-"FO"
names(FOWNODFZ)<-c("WNODFZ","site")

obs <- unlist(networklevel(KP, index="weighted NODF"))
nulls <- nullmodel(KP, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
KPWNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
KPWNODFZ$site<-"KP"
names(KPWNODFZ)<-c("WNODFZ","site")

obs <- unlist(networklevel(R4, index="weighted NODF"))
nulls <- nullmodel(R4, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
R4WNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
R4WNODFZ$site<-"R4"
names(R4WNODFZ)<-c("WNODFZ","site")

obs <- unlist(networklevel(RJ, index="weighted NODF"))
nulls <- nullmodel(RJ, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
RJWNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
RJWNODFZ$site<-"RJ"
names(RJWNODFZ)<-c("WNODFZ","site")

obs <- unlist(networklevel(RL, index="weighted NODF"))
nulls <- nullmodel(RL, N=1000, method=3)
null <- unlist(sapply(nulls, networklevel, index="weighted NODF")) 
RLWNODFZ<-as.data.frame(((obs-mean(null))/sd(null)))
RLWNODFZ$site<-"RL"
names(RLWNODFZ)<-c("WNODFZ","site")

nestedness<-rbind(CBWNODFZ,COWNODFZ,CPWNODFZ,DWWNODFZ,EIWNODFZ,FCWNODFZ,FOWNODFZ,KPWNODFZ,R4WNODFZ,RJWNODFZ,RLWNODFZ)


###################### Modularity z-score ####################################
mod <- computeModules(web=CB)
nulls <- nullmodel(CB, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
CBMODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
CBMODZ$site<-"CB"
names(CBMODZ)<-c("modularity.z","site")

mod <- computeModules(web=CO)
nulls <- nullmodel(CO, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
COMODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
COMODZ$site<-"CO"
names(COMODZ)<-c("modularity.z","site")

mod <- computeModules(web=CP)
nulls <- nullmodel(CP, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
CPMODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
CPMODZ$site<-"CP"
names(CPMODZ)<-c("modularity.z","site")

mod <- computeModules(web=DW)
nulls <- nullmodel(DW, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
DWMODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
DWMODZ$site<-"DW"
names(DWMODZ)<-c("modularity.z","site")

mod <- computeModules(web=EI)
nulls <- nullmodel(EI, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
EIMODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
EIMODZ$site<-"EI"
names(EIMODZ)<-c("modularity.z","site")

mod <- computeModules(web=FC)
nulls <- nullmodel(FC, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
FCMODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
FCMODZ$site<-"FC"
names(FCMODZ)<-c("modularity.z","site")

mod <- computeModules(web=FO)
nulls <- nullmodel(FO, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
FOMODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
FOMODZ$site<-"FO"
names(FOMODZ)<-c("modularity.z","site")

mod <- computeModules(web=KP)
nulls <- nullmodel(KP, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
KPMODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
KPMODZ$site<-"KP"
names(KPMODZ)<-c("modularity.z","site")

mod <- computeModules(web=R4)
nulls <- nullmodel(R4, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
R4MODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
R4MODZ$site<-"R4"
names(R4MODZ)<-c("modularity.z","site")

mod <- computeModules(web=RJ)
nulls <- nullmodel(RJ, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
RJMODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
RJMODZ$site<-"RJ"
names(RJMODZ)<-c("modularity.z","site")

mod <- computeModules(web=RL)
nulls <- nullmodel(RL, N=1000, method="vaznull")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
RLMODZ<-as.data.frame(mod@likelihood - mean(like.nulls))/sd(like.nulls)
RLMODZ$site<-"RL"
names(RLMODZ)<-c("modularity.z","site")

modularity<-rbind(CBMODZ,COMODZ,CPMODZ,DWMODZ,EIMODZ,FCMODZ,FOMODZ,KPMODZ,R4MODZ,RJMODZ,RLMODZ)


##################### Connectance value ##########################
#connectance constant in vaznull, so using original value instead of null z-score

CBCON<-as.data.frame(networklevel(CB, index="connectance"))
CBCON$site<-"CB"
names(CBCON)<-c("connectance","site")

COCON<-as.data.frame(networklevel(CO, index="connectance"))
COCON$site<-"CO"
names(COCON)<-c("connectance","site")

CPCON<-as.data.frame(networklevel(CP, index="connectance"))
CPCON$site<-"CP"
names(CPCON)<-c("connectance","site")

DWCON<-as.data.frame(networklevel(DW, index="connectance"))
DWCON$site<-"DW"
names(DWCON)<-c("connectance","site")

EICON<-as.data.frame(networklevel(EI, index="connectance"))
EICON$site<-"EI"
names(EICON)<-c("connectance","site")

FCCON<-as.data.frame(networklevel(FC, index="connectance"))
FCCON$site<-"FC"
names(FCCON)<-c("connectance","site")

FOCON<-as.data.frame(networklevel(FO, index="connectance"))
FOCON$site<-"FO"
names(FOCON)<-c("connectance","site")

KPCON<-as.data.frame(networklevel(KP, index="connectance"))
KPCON$site<-"KP"
names(KPCON)<-c("connectance","site")

R4CON<-as.data.frame(networklevel(R4, index="connectance"))
R4CON$site<-"R4"
names(R4CON)<-c("connectance","site")

RJCON<-as.data.frame(networklevel(RJ, index="connectance"))
RJCON$site<-"RJ"
names(RJCON)<-c("connectance","site")

RLCON<-as.data.frame(networklevel(RL, index="connectance"))
RLCON$site<-"RL"
names(RLCON)<-c("connectance","site")

connectance<-rbind(CBCON,COCON,CPCON,DWCON,EICON,FCCON,FOCON,KPCON,R4CON,RJCON,RLCON)
##################### Species richness value ##########################

CBRICH<-as.data.frame(networklevel(CB, index="number of species"))
CBRICH$site<-"CB"
CBRICH<-CBRICH[1,]
names(CBRICH)<-c("number of species","site")

CORICH<-as.data.frame(networklevel(CO, index="number of species"))
CORICH$site<-"CO"
CORICH<-CORICH[1,]
names(CORICH)<-c("number of species","site")

CPRICH<-as.data.frame(networklevel(CP, index="number of species"))
CPRICH$site<-"CP"
CPRICH<-CPRICH[1,]
names(CPRICH)<-c("number of species","site")

DWRICH<-as.data.frame(networklevel(DW, index="number of species"))
DWRICH$site<-"DW"
DWRICH<-DWRICH[1,]
names(DWRICH)<-c("number of species","site")

EIRICH<-as.data.frame(networklevel(EI, index="number of species"))
EIRICH$site<-"EI"
EIRICH<-EIRICH[1,]
names(EIRICH)<-c("number of species","site")

FCRICH<-as.data.frame(networklevel(FC, index="number of species"))
FCRICH$site<-"FC"
FCRICH<-FCRICH[1,]
names(FCRICH)<-c("number of species","site")

FORICH<-as.data.frame(networklevel(FO, index="number of species"))
FORICH$site<-"FO"
FORICH<-FORICH[1,]
names(FORICH)<-c("number of species","site")

KPRICH<-as.data.frame(networklevel(KP, index="number of species"))
KPRICH$site<-"KP"
KPRICH<-KPRICH[1,]
names(KPRICH)<-c("number of species","site")

R4RICH<-as.data.frame(networklevel(R4, index="number of species"))
R4RICH$site<-"R4"
R4RICH<-R4RICH[1,]
names(R4RICH)<-c("number of species","site")

RJRICH<-as.data.frame(networklevel(RJ, index="number of species"))
RJRICH$site<-"RJ"
RJRICH<-RJRICH[1,]
names(RJRICH)<-c("number of species","site")

RLRICH<-as.data.frame(networklevel(RL, index="number of species"))
RLRICH$site<-"RL"
RLRICH<-RLRICH[1,]
names(RLRICH)<-c("number of species","site")

richness<-rbind(CBRICH,CORICH,CPRICH,DWRICH,EIRICH,FCRICH,FORICH,KPRICH,R4RICH,RJRICH,RLRICH)

# now tying all of the site level metrics together
site_indices<-merge(connectance,modularity, by="site")
site_indices<-merge(site_indices,nestedness, by="site")
site_indices<-merge(site_indices,richness, by="site")
site_indices<-merge(site_indices,H2, by="site")

#now merge with site level dataset (includes pathogen prevalence, agriculture, etc.)





