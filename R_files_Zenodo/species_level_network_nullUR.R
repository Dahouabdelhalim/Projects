setwd("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Cornell/bee_pathogens/2015-FINAL/science/verified")

#these are the unresolved data (network maintains "genus" only identification possible in the field)

wf<-read.csv("visitation_2015UR.csv")
pathogens<-read.csv("pathogens_2015UR.csv")

require(bipartite)
require(sna)
library(dplyr)

#exclude date and flower species information to fit bipartite requirements
net<- wf[c(-1,-2)]

#most simple network: all interactions, at all sites, throughout summer
WFsum <- aggregate(. ~ Species, data=net, FUN=sum)
row.names(WFsum) <- WFsum[,1]
WFsum <-WFsum[,c(-1)]

# plots separated by site 
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
CB <-CB[,(-1)]
CO <- aggregate(. ~ Species, data=CO, FUN=sum)
row.names(CO) <- CO[,1]
CO <-CO[,(-1)]
CP <- aggregate(. ~ Species, data=CP, FUN=sum)
row.names(CP) <- CP[,1]
CP <-CP[,(-1)]
R4 <- aggregate(. ~ Species, data=R4, FUN=sum)
row.names(R4) <- R4[,1]
R4 <-R4[,(-1)]
RL <- aggregate(. ~ Species, data=RL, FUN=sum)
row.names(RL) <- RL[,1]
RL <-RL[,(-1)]
RJ <- aggregate(. ~ Species, data=RJ, FUN=sum)
row.names(RJ) <- RJ[,1]
RJ <-RJ[,(-1)]
DW <- aggregate(. ~ Species, data=DW, FUN=sum)
row.names(DW) <- DW[,1]
DW <-DW[,(-1)]
EI <- aggregate(. ~ Species, data=EI, FUN=sum)
row.names(EI) <- EI[,1]
EI <-EI[,(-1)]
KP <- aggregate(. ~ Species, data=KP, FUN=sum)
row.names(KP) <- KP[,1]
KP <-KP[,(-1)]
FC <- aggregate(. ~ Species, data=FC, FUN=sum)
row.names(FC) <- FC[,1]
FC <-FC[,(-1)]
FO <- aggregate(. ~ Species, data=FO, FUN=sum)
row.names(FO) <- FO[,1]
FO <-FO[,(-1)]

##################### Bee species level network indices  ##########################
CB<-CB[, colSums(CB) > 0]
indices <- specieslevel(CB, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(CB, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 12, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(CB)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matCB<-as.data.frame(z.mat)
z.matCB$site<-"CB"

CO<-CO[, colSums(CO) > 0]
indices <- specieslevel(CO, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(CO, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 18, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(CO)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matCO<-as.data.frame(z.mat)
z.matCO$site<-"CO"

CP<-CP[, colSums(CP) > 0]
indices <- specieslevel(CP, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(CP, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 13, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(CP)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matCP<-as.data.frame(z.mat)
z.matCP$site<-"CP"

DW<-DW[, colSums(DW) > 0]
indices <- specieslevel(DW, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(DW, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 11, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(DW)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matDW<-as.data.frame(z.mat)
z.matDW$site<-"DW"

EI<-EI[, colSums(EI) > 0]
indices <- specieslevel(EI, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(EI, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 7, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(EI)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matEI<-as.data.frame(z.mat)
z.matEI$site<-"EI"

FC<-FC[, colSums(FC) > 0]
indices <- specieslevel(FC, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(FC, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 11, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(FC)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matFC<-as.data.frame(z.mat)
z.matFC$site<-"FC"

FO<-FO[, colSums(FO) > 0]
indices <- specieslevel(FO, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(FO, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 12, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(FO)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matFO<-as.data.frame(z.mat)
z.matFO$site<-"FO"

KP<-KP[, colSums(KP) > 0]
indices <- specieslevel(KP, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(KP, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 13, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(KP)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matKP<-as.data.frame(z.mat)
z.matKP$site<-"KP"

R4<-R4[, colSums(R4) > 0]
indices <- specieslevel(R4, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(R4, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 12, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(R4)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matR4<-as.data.frame(z.mat)
z.matR4$site<-"R4"

RJ<-RJ[, colSums(RJ) > 0]
indices <- specieslevel(RJ, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(RJ, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 18, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(RJ)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matRJ<-as.data.frame(z.mat)
z.matRJ$site<-"RJ"

RL<-RL[, colSums(RL) > 0]
indices <- specieslevel(RL, index="ALLBUTD")$"higher level"
# generate null models and calculate their index values:
exanulls <- nullmodel(RL, N=1000, method=3)
exanull.res <- lapply(exanulls, function(x) specieslevel(x, index="ALLBUTD")$"higher level")
# calculate z-scores:
z.mat <- matrix(0, 16, 20)
colnames(z.mat) <- colnames(indices)
rownames(z.mat) <- colnames(RL)
for (i in 1:3){
  mean.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, mean, na.rm=T)
  sd.n <- apply(sapply(exanull.res, function(x) x[,i]), 1, sd, na.rm=T)
  z <- (indices[,i] - mean.n)/ifelse(sd.n==0, 1, sd.n) # sd is 0 when null model values are
  z.mat[,i] <- z                                       #constant (i.e. mean.n==index);then the z-score should be 0
}
z.matRL<-as.data.frame(z.mat)
z.matRL$site<-"RL"

full<-list(z.matCB,z.matCO,z.matCP,z.matDW,z.matEI,z.matFC,z.matFO,z.matKP,z.matR4,
           z.matRJ,z.matRL)
full<-bind_rows(lapply(full, add_rownames))
full$GenSp<-full$rowname
full$spst<- paste(full$GenSp, full$site, sep='.')

pathogens$spst<-paste(pathogens$genus, pathogens$species, pathogens$site, sep='.')

total<-merge(full,pathogens,by="spst")
#write.csv(total,"beenetspeciesUR.csv")