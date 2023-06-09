setwd("~/Dropbox/11.Data_Stats_Benoit_Dropbox/2020-12-17/")
## Libraries
library(MASS)
library(lme4)
library(vegan)
library(DHARMa)
library(RVAideMemoire)
library(emmeans)
library(pls)
library(rareNMtests)
library(SPECIES)
library(SpadeR)

## Functions
NORM <- function(m){
  qqnorm(residuals(m))
  qqline(residuals(m))
  shapiro.test(residuals(m))
}
se <- function(x){sd(x)/sqrt(length(x))}

###################### A/ FLORAL TRAITS ---------------------------------------------------
######## 1/ FLORAL SCENT -------------------------------------------------------------------------------------------------------
### 1.1/ Rosemary RO ---------------------------------------------------------------
RO3 <- read.csv("./2021-01-08_RO3.csv"); summary(RO3) # data set VOCs in Salvia rosmarinus floral samples (µg.min-1)
LRO3 <- read.csv("./2021-01-08_LRO3.csv"); summary(LRO3) # data set VOCs in Salvia rosmarinus leaf samples (µg.min-1)
tro <- dim(RO3)[2]

# Standardize by dry mass and by temperature :
RO3c <- RO3[,13:tro]*60/RO3$MassTot*exp(0.09*(303-(RO3$TempIn+273.15))) # µg.gDM-1.h-1
LRO3c <- LRO3[,13:tro]*60/LRO3$MassTot*exp(0.09*(303-(LRO3$TempIn+273.15))) # µg.gDM-1.h-1
# Transform variables
RO4c <- (RO3c)^(1/4)
LRO4c <- (LRO3c)^(1/4)
# Remove rare VOCs:
ROFLOvsLEAF <- data.frame(FLO=apply(RO4c,2,function(x)length(which(x>0))),LEAF=apply(LRO4c,2,function(x)length(which(x>0))))
ROFLOvsLEAF[which(ROFLOvsLEAF$FLO<3),]
ROFLOvsLEAF[which(ROFLOvsLEAF$LEAF==0),]

RO5c <- RO4c[,-which(ROFLOvsLEAF$FLO<3)]
LRO5c <- LRO4c[,-which(ROFLOvsLEAF$FLO<3)]

ROtab <- data.frame(VOC=sort(names(RO3c)),Mean=round(as.matrix(apply(RO3c,2,mean)),digits=4)[order(names(RO3c))],
                    SE=round(as.matrix(apply(RO3c,2,se)),digits=4)[order(names(RO3c))],
                    Nsup0=c(as.matrix(apply(RO3c,2,function(x)length(which(x>0)))))[order(names(RO3c))])


#### 1. Drought impact on floral scent multivariate analysis:
## 1.1. All VOCs
MAC <- RO5c; Treat <- RO3$Treatment
plsda1c <- cppls(dummy(Treat)~as.matrix(MAC),ncomp=3,scale=T,center=T)
cp <- MVA.synt(plsda1c); cp
# Test (run 4 times, 999 permutations each)
MVA.test(MAC,Treat,ncomp=3,model="PPLS-DA",cmv=TRUE)
par(mfrow=c(1,1),bty="l")
MVA.plot(plsda1c,fac=Treat,col=c("#B69BDE","mediumpurple4"),points=F,xax=1,yax=2,fac.lab="",fac.cex=0.001,pch=c(1,2))
MVA.plot(plsda1c,fac=Treat,col=c("#B69BDE","mediumpurple4"),points=F,xax=1,yax=3,fac.lab="",fac.cex=0.001,pch=c(1,2))
# Select VOCs contributing most to group separation (highest contribution along barycentres' axis):
PLRO5c <- data.frame(VOC=names(RO5c),plsda1c$projection)
names(PLRO5c)[2:4] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,3],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLRO5c$xprim <- (PLRO5c$NC1*X0[2]+PLRO5c$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLRO5c$xsec <- (PLRO5c$xprim*Xp0[2]+PLRO5c$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2)
PLRO5c$Percxx2 <- round(abs(PLRO5c$xsec)/sum(abs(PLRO5c$xsec))*100,digits=2) # change of coordinates (Xprim,Z) > Xsec
PLRO5c[order(PLRO5c$Percxx2),]
t1xx <- which(PLRO5c$Percxx2>4); length(t1xx); sum(PLRO5c$Percxx2[t1xx]) # 11 VOCs=48% of VOCs; 65.2% contrib axis Xp0
RO6c <- RO5c[,t1xx]

## 1.2 Selected VOCs RO6c final matrix:
MAC <- RO6c; Treat <- RO3$Treatment
plsda1c <- cppls(dummy(RO3$Treatment)~as.matrix(MAC),ncomp=3,scale=T,center=T)
cp <- MVA.synt(plsda1c); cp
MVA.test(MAC,Treat,ncomp=3,model="PPLS-DA",cmv=TRUE)
#save(plsda1c,file="2021-01-11_RO6c_plsda1c.Rdata")
PLRO6c <- data.frame(VOC=names(RO6c),plsda1c$projection)
names(PLRO6c)[2:4] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,3],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLRO6c$xprim <- (PLRO6c$NC1*X0[2]+PLRO6c$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLRO6c$xsec <- (PLRO6c$xprim*Xp0[2]+PLRO6c$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLRO6c$Percxx2 <- round(abs(PLRO6c$xsec)/sum(abs(PLRO6c$xsec))*100,digits=2)
PLRO6c[6:7] # Contribution of the 18 VOCs to group separation, and distinction of typical "control" VOCs (xsec<0) vs. "dry" VOCs (xsec>0)


#### 2. Total emissions: 
RO3S <- apply(RO3c,1,sum)
RO3D <- apply(RO3c,1,function(x)length(which(x>0)))
RO4S <- (RO3S)^(1/4)
t.test(RO4S~RO3$Treatment) # t=-0.214, df=20.6, P=0.83
plot(RO4S~RO3$Treatment)
tapply(RO3S,RO3$Treatment,mean)
tapply(RO3S,RO3$Treatment,se)
t.test(RO3D~RO3$Treatment) # t=-0.615 df=20.5 P=0.54
plot(RO3D~RO3$Treatment)
tapply(RO3D,RO3$Treatment,mean)
tapply(RO3D,RO3$Treatment,se)

# Define chemical family for each VOC:
ROFam <- data.frame(Family1=as.factor(c("MONOalc","MONOcet","MONOalc","MONOalc","ARalc","SESQ","SESQ","FADalc","AR","MONO",
                                        "MONO","MONOalc","MONO","MONOalc","MONOalc","MONOest","MONOest","MONO","MONOalc","MONOoxd",
                                        "MONO","MONO","MONOcet","SESQoxd","MONO","MONO","MONOcet","MONOalc","MONOalc","MONOcet",
                                        "MONOcet","MONOald","FADaci","SESQ","FADaci","MONOald","MONOcet")))
ROFam$Family2 <- ROFam$Family1
levels(ROFam$Family2) <- c("AR","ARox",rep("FADox",2),"MONO",rep("MONOox",5),"SESQ","SESQox")
ROFam$Family3 <- ROFam$Family2
levels(ROFam$Family3) <- c("AR","AR","FAD","MONO","MONO","SESQ","SESQ")
ROFam$Family4 <- ROFam$Family3
levels(ROFam$Family4) <- c("AR","FAD","TERP","TERP")
ROMoy <- data.frame(MONOsimp=c(apply(RO3c[,which(ROFam$Family1=="MONO")],1,sum)),
                    MONOox=c(apply(RO3c[,which(ROFam$Family2=="MONOox")],1,sum)),
                    FADox=c(apply(RO3c[,which(ROFam$Family2=="FADox")],1,sum)),
                    SESQ=c(apply(RO3c[,which(ROFam$Family1=="SESQ")],1,sum)))
ROMoy2 <- ROMoy^(1/4)
ROMoy2$Ech2 <- RO3$Ech2
ROMoy2$Treatment <- RO3$Treatment
plot(ROMoy2$MONOsimp~ROMoy2$Treatment)
t.test(ROMoy2$MONOsimp~ROMoy2$Treatment) # t=0.573, df=13.0, P=0.58
tapply(ROMoy$MONOsimp,ROMoy2$Treatment,mean); tapply(ROMoy$MONOsimp,ROMoy2$Treatment,se)
plot(ROMoy2$MONOox~ROMoy2$Treatment)
t.test(ROMoy2$MONOox~ROMoy2$Treatment) # t=0.235, df=20.5, P=0.82
tapply(ROMoy$MONOox,ROMoy2$Treatment,mean); tapply(ROMoy$MONOox,ROMoy2$Treatment,se)
plot(ROMoy2$FADox~ROMoy2$Treatment)
t.test(ROMoy2$FADox~ROMoy2$Treatment) # t=-4.39, df=18.3, P=0.00034 ***
tapply(ROMoy$FADox,ROMoy2$Treatment,mean); tapply(ROMoy$FADox,ROMoy2$Treatment,se)
plot(ROMoy2$SESQ~ROMoy2$Treatment)
t.test(ROMoy2$SESQ~ROMoy2$Treatment) # t=-1.06, df=21.0, P=0.30
tapply(ROMoy$SESQ,ROMoy2$Treatment,mean); tapply(ROMoy$SESQ,ROMoy2$Treatment,se)


#### 3. Leaf VOCs vs. Flower VOCs mutlivariate analysis
MAC <- rbind(RO5c,LRO5c); Treat <- as.factor(c(as.character(RO3$Type),as.character(LRO3$Type)))
plsda1c <- cppls(dummy(Treat)~as.matrix(MAC),ncomp=3,scale=T,center=T)
cp <- MVA.synt(plsda1c); cp
par(mfrow=c(1,1),bty="l",mar=c(4,4,1,1))
MVA.plot(plsda1c,fac=Treat,col=c("#cc99FFcc","mediumpurple4"),points=F,xax=1,yax=2,fac.lab="",fac.cex=0.001,pch=c(1,2))
# Select VOCs contributing most to group separation (highest contribution along barycentres' axis):
PLROLc <- data.frame(plsda1c$projection)
names(PLROLc)[1:3] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,3],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLROLc$xprim <- (PLROLc$NC1*X0[2]+PLROLc$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLROLc$xsec <- (PLROLc$xprim*Xp0[2]+PLROLc$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLROLc$Percxx2 <- round(abs(PLROLc$xsec)/sum(abs(PLROLc$xsec))*100,digits=2)
PLROLc[order(PLROLc$Percxx2),]
t1 <- which(PLROLc$Percxx2>3.2); t1; length(t1); sum(PLROLc[t1,]$Percxx2) # 14 VOCs selected, representing 66.8  % of contrib. to axis Xp0
# cca tests :
MAC <- rbind(RO5c,LRO5c)[,t1]; Treat <- as.factor(c(as.character(RO3$Type),as.character(LRO3$Type)))
plsda1c <- cppls(dummy(Treat)~as.matrix(MAC),ncomp=3,scale=T,center=T)#,Y.add=Ya)
cp <- MVA.synt(plsda1c); cp
PLROLc2 <- data.frame(plsda1c$projection)
names(PLROLc2)[1:3] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,3],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLROLc2$xprim <- (PLROLc2$NC1*X0[2]+PLROLc2$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLROLc2$xsec <- (PLROLc2$xprim*Xp0[2]+PLROLc2$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLROLc2$Percxx2 <- round(abs(PLROLc2$xsec)/sum(abs(PLROLc2$xsec))*100,digits=2)
PLROLc2[5:6] # Contribution of the 17 VOCs to group separation, and distinction of typical leafy VOCs (xsec<0) vs. floral VOCs (xsec>0)
aaaa=cca(MAC~Treat)
anova.cca(aaaa) # df=1,25, F=2.51, P=0.028 * (!!! because of cross-validation & permutation, values change at each run)

# Total emissions:
LRO3S <- apply(LRO3c,1,sum)
t.test((RO3S)^(1/4),LRO3S^(1/4)) # t = 2.01 df =3.48, P=0.13
mean(RO3S); mean(LRO3S)
se(RO3S); se(LRO3S)
LRO3D <- apply(LRO3c,1,function(x)length(which(x>0)))
t.test(RO3D,LRO3D) # t = 2.90 df = 5.94, P=0.028 *
mean(RO3D); mean(LRO3D) # 19 / 13
se(RO3D); se(LRO3D) # 1 / 2


### 1.2/ Cistus CI ---------------------------------------------------------------
CI3 <- read.csv("./2021-01-08_CI3.csv"); summary(CI3) # data set VOCs in Cistus albidus floral samples (µg.min-1)
LCI3 <- read.csv("./2021-01-08_LCI3.csv"); summary(LCI3) # data set VOCs in Cistus albidus leaf samples (µg.min-1)
tCI <- dim(CI3)[2]
# Standardize by dry mass and by temperature :
CI3c <- CI3[,13:tCI]*60/CI3$MassTot*exp(0.09*(303-(CI3$TempIn+273.15))) # µg.gDM-1.h-1
LCI3c <- LCI3[,13:tCI]*60/CI3$MassTot*exp(0.09*(303-(LCI3$TempIn+273.15))) # µg.gDM-1.h-1
# Transform variables
CI4c <- (CI3c)^(1/4)
LCI4c <- (LCI3c)^(1/4)
# Remove rare VOCs:
CIFLOvsLEAF <- data.frame(FLO=apply(CI4c,2,function(x)length(which(x>0))),LEAF=apply(LCI4c,2,function(x)length(which(x>0))))
CIFLOvsLEAF[which(CIFLOvsLEAF$FLO<3),]
CIFLOvsLEAF[which(CIFLOvsLEAF$LEAF==0),]

CI5c <- CI4c[,-which(CIFLOvsLEAF$FLO<3)]
LCI5c <- LCI4c[,-which(CIFLOvsLEAF$FLO<3)]

CItab <- data.frame(VOC=sort(names(CI3c)),Mean=round(as.matrix(apply(CI3c,2,mean)),digits=3)[order(names(CI3c))],
                    SE=round(as.matrix(apply(CI3c,2,se)),digits=3)[order(names(CI3c))],
                    Nsup0=c(as.matrix(apply(CI3c,2,function(x)length(which(x>0)))))[order(names(CI3c))])



#### 1. Drought impact of floral scent multivariate analysis:
## 1.1. All VOCs
MAC <- CI5c; Treat <- CI3$Treatment
plsda1c <- cppls(dummy(Treat)~as.matrix(MAC),ncomp=3,center=T,scale=T,y=T)
cp <- MVA.synt(plsda1c); cp
# Test (run 4 times, 999 permutations each)
MVA.test(as.matrix(CI5c),CI3$Treatment,model="PPLS-DA",cmv=TRUE,ncomp=3)
par(mfrow=c(1,1))
MVA.plot(plsda1c,fac=CI3$Treatment,col=c("lightpink1","deeppink4"),points=FALSE,xax=1,yax=2,fac.lab="",fac.cex=0.001)
# Select VOCs contributing most to group separation (highest contribution along barycentres' axis):
PLCI5c <- data.frame(VOC=names(CI5c),plsda1c$projection)
names(PLCI5c)[2:4] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,2],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLCI5c$xprim <- (PLCI5c$NC1*X0[2]+PLCI5c$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLCI5c$yprim <- (PLCI5c$NC2*Y0[2]+PLCI5c$NC3*Z0[2])/sqrt(Y0[2]^2+Z0[2]^2)
PLCI5c$xsec <- (PLCI5c$xprim*Xp0[2]+PLCI5c$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLCI5c$Percxx2 <- round(abs(PLCI5c$xsec)/sum(abs(PLCI5c$xsec))*100,digits=2)
PLCI5c[order(PLCI5c$Percxx2),]
t1xx4 <- which(PLCI5c$Percxx2>5); length(t1xx4); sum(PLCI5c$Percxx2[t1xx4]) # 7 VOCs selected (37.5%) representing 65.88% of contrib. to axis
CI6c <- CI5c[,t1xx4]

## 1.2 Selected VOCs CI6c final matrix:
MAC <- CI6c; Treat <- CI3$Treatment
plsda1c <- cppls(dummy(CI3$Treatment)~as.matrix(MAC),ncomp=3,center=T,scale=T)
cp <- MVA.synt(plsda1c); cp
#save(plsda1c,file="2021-01-11_CI63cplsda1c.Rdata")
MVA.cmv(MAC,CI3$Treatment,model="PPLS-DA",crit.inn="NMC",ncomp=3) # classification error rate 33% (6 individuals)
MVA.test(as.matrix(CI6c),CI3$Treatment,model="PPLS-DA",cmv=TRUE,ncomp=3)
par(mfrow=c(1,1),xpd=F,las=0,mar=c(4,4,1,1),bty="l")
MVA.plot(plsda1c,fac=CI3$Treatment,col=c("plum","magenta4"),points=FALSE,xax=1,yax=2,fac.lab="",fac.cex=0.001)
PLCI6c <- data.frame(VOC=names(CI6c),plsda1c$projection)
names(PLCI6c)[2:4] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,2],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLCI6c$xprim <- (PLCI6c$NC1*X0[2]+PLCI6c$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLCI6c$xsec <- (PLCI6c$xprim*Xp0[2]+PLCI6c$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLCI6c$Percxx2 <- round(abs(PLCI6c$xsec)/sum(abs(PLCI6c$xsec))*100,digits=2)
PLCI6c[6:7] # Contribution of the 18 VOCs to group separation, and distinction of typical "control" VOCs (xsec<0) vs. "dry" VOCs (xsec>0)


#### 2. Total emissions: 
CI3S <- apply(CI3c,1,sum)
CI3D <- apply(CI3c,1,function(x)length(which(x>0)))
CI4S <- (CI3S)^(1/4)
t.test(CI4S~CI3$Treatment) # t=0.527, df=20.0, P=0.60
plot(CI4S~CI3$Treatment)
tapply(CI3S,CI3$Treatment,mean)
tapply(CI3S,CI3$Treatment,se)
t.test(CI3D~CI3$Treatment) # t=1.76 df=19.72 P=0.093 .
tapply(CI3D,CI3$Treatment,mean)
tapply(CI3D,CI3$Treatment,se)
CIFam <- data.frame(Family1=as.factor(c(rep("SESQ",6),"FADalc","FADest","SESQ","SESQ",
                                        "MONO",rep("SESQ",5),"FADest","MONOcet","MONOest","SESQ",
                                        "MONOest",rep("SESQ",4),"SESQcet","SESQalc","SESQalc","SESQest")))
cbind(CIFam,names(CI3c))
CIFam$Family2 <- CIFam$Family1
levels(CIFam$Family2) <- c(rep("FADox",2),"MONO",rep("MONOox",2),"SESQ",rep("SESQox",3))
CIFam$Family3 <- CIFam$Family2
levels(CIFam$Family3) <- c("FAD","MONO","MONO","SESQ","SESQ")
CIFam$Family4 <- CIFam$Family3
levels(CIFam$Family4) <- c("FAD","TERP","TERP")
CIMoy <- data.frame(SESQsimp=c(apply(CI3c[,which(CIFam$Family1=="SESQ")],1,sum)),
                    SESQox=c(apply(CI3c[,which(CIFam$Family2=="SESQox")],1,sum)),
                    MONOox=c(apply(CI3c[,which(CIFam$Family2=="MONOox")],1,sum)),
                    FADox=c(apply(CI3c[,which(CIFam$Family2=="FADox")],1,sum)))
CIMoy2 <- CIMoy^(1/4)
CIMoy2$Ech2 <- CI3$Ech2
CIMoy2$Treatment <- CI3$Treatment
plot(CIMoy2$MONOox~CIMoy2$Treatment)
t.test(CIMoy2$MONOox~CIMoy2$Treatment) # t=-0.179, df=19.9, P=0.86
tapply(CIMoy$MONOox,CIMoy2$Treatment,mean); tapply(CIMoy$MONOox,CIMoy2$Treatment,se)
plot(CIMoy2$SESQsimp~CIMoy2$Treatment)
t.test(CIMoy2$SESQsimp~CIMoy2$Treatment) # t=1.41, df=19.9, P=0.17
tapply(CIMoy$SESQsimp,CIMoy2$Treatment,mean); tapply(CIMoy$SESQsimp,CIMoy2$Treatment,se)
plot(CIMoy2$SESQox~CIMoy2$Treatment)
t.test(CIMoy2$SESQox~CIMoy2$Treatment) # t=0.661, df=17.9, P0.52
tapply(CIMoy$SESQox,CIMoy2$Treatment,mean); tapply(CIMoy$SESQox,CIMoy2$Treatment,se)
plot(CIMoy2$FADox~CIMoy2$Treatment)
t.test(CIMoy2$FADox~CIMoy2$Treatment) # t=-0.448, df=16.2, P=0.66
tapply(CIMoy$FADox,CIMoy2$Treatment,mean); tapply(CIMoy$FADox,CIMoy2$Treatment,se)


#### 3. Leaf VOCs vs. Flower VOCs mutlivariate analysis
MAC <- rbind(CI5c,LCI5c); Treat <- as.factor(c(as.character(CI3$Type),as.character(LCI3$Type)))
plsda1c <- cppls(dummy(Treat)~as.matrix(MAC),ncomp=3,scale=T,center=T)#,Y.add=Ya)
cp <- MVA.synt(plsda1c); cp
par(mfrow=c(1,1),bty="l",mar=c(4,4,1,1),las=1)
MVA.plot(plsda1c,fac=Treat,col=c("deeppink1","deeppink4"),points=F,xax=1,yax=2,fac.lab="",fac.cex=0.001,pch=c(1,2))
# Select VOCs contributing most to group separation (highest contribution along barycentres' axis):
PLCILc <- data.frame(plsda1c$projection)
names(PLCILc)[1:3] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,3],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLCILc$xprim <- (PLCILc$NC1*X0[2]+PLCILc$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLCILc$xsec <- (PLCILc$xprim*Xp0[2]+PLCILc$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLCILc$Percxx2 <- round(abs(PLCILc$xsec)/sum(abs(PLCILc$xsec))*100,digits=2)
PLCILc[order(PLCILc$Percxx2),]
t0 <- which(PLCILc$Percxx2>7); PLCILc[t0,]
length(c(t0,2)); sum(PLCILc$Percxx2[c(t0,2)]) # 6+1 VOCs selected, representing 70.3% of contrib. axis
# (2=E-Caryophyllene, added because other in t0 have empty rows.)

# cca tests :
MAC <- rbind(CI5c,LCI5c)[,c(t0,2)]; Treat <- as.factor(c(as.character(CI3$Type),as.character(LCI3$Type)))
plsda1c <- cppls(dummy(Treat)~as.matrix(MAC),ncomp=3,scale=T,center=T)#,Y.add=Ya)
cp <- MVA.synt(plsda1c); cp
aaaa=cca(MAC~Treat) # CCA > AFCs (tableau contingence)
anova.cca(aaaa) # df=1,24, F=3.51, P=0.027 *
PLCILc2 <- data.frame(plsda1c$projection)
names(PLCILc2)[1:3] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,3],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLCILc2$xprim <- (PLCILc2$NC1*X0[2]+PLCILc2$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLCILc2$xsec <- (PLCILc2$xprim*Xp0[2]+PLCILc2$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLCILc2$Percxx2 <- round(abs(PLCILc2$xsec)/sum(abs(PLCILc2$xsec))*100,digits=2)
PLCILc2[5:6] # Contribution of the 17 VOCs to group separation, and distinction of typical leafy VOCs (xsec<0) vs. floral VOCs (xsec>0)

# Total emissions:
LCI3S <- apply(LCI3c,1,sum)
mean(CI3S); mean(LCI3S)
se(CI3S); se(LCI3S)
t.test(CI3S^(1/4),LCI3S^(1/4)) # t=-4.00 df =9.23 P=0.0029 **
LCI3D <- apply(LCI3c,1,function(x)length(which(x>0)))
t.test(CI3D,LCI3D) # t=-1.24, df=4.33 P=0.28
mean(CI3D); mean(LCI3D)
se(CI3D); se(LCI3D)


### 1.3/ Thyme TH ---------------------------------------------------------------
TH3 <- read.csv("./2021-01-08_TH3.csv"); summary(TH3)
LTH3 <- read.csv("./2021-01-08_LTH3.csv"); summary(LTH3)
tTH <- dim(TH3)[2]
TH3$Treatment <- as.factor(as.character(TH3$Treatment))
# Standardize by dry mass and by temperature :
TH3c <- TH3[,17:tTH]*60/TH3$MassTot*exp(0.09*(303-(TH3$TempIn+273.15))) # µg.gDM-1.h-1
LTH3c <- LTH3[,17:tTH]*60/LTH3$MassTot*exp(0.09*(303-(LTH3$TempIn+273.15))) # µg.gDM-1.h-1
# Transform variables
TH4c <- (TH3c)^(1/4)
LTH4c <- (LTH3c)^(1/4)
# Remove rare VOCs:
THFLOvsLEAF <- data.frame(FLO=apply(TH4c,2,function(x)length(which(x>0))),LEAF=apply(LTH4c,2,function(x)length(which(x>0))))
THFLOvsLEAF[which(THFLOvsLEAF$FLO<3),]
THFLOvsLEAF[which(THFLOvsLEAF$LEAF==0),]
THFLOvsLEAF[which(THFLOvsLEAF$FLO==0),]

TH5c <- TH4c[,-which(THFLOvsLEAF$FLO<3)]
LTH5c <- LTH4c[,-which(THFLOvsLEAF$FLO<3)]

THtab <- data.frame(VOC=sort(names(TH3c)),Mean=round(as.matrix(apply(TH3c,2,mean)),digits=4)[order(names(TH3c))],
                    SE=round(as.matrix(apply(TH3c,2,se)),digits=4)[order(names(TH3c))],
                    Nsup0=c(as.matrix(apply(TH3c,2,function(x)length(which(x>0)))))[order(names(TH3c))])


#### 1. Drought impact of floral scent multivariate analysis:
## 1.1. All VOCs
MAC <- TH5c; Treat <- TH3$Treatment; Ya <- dummy(TH3$Sex)
plsda1c <- cppls(dummy(Treat)~as.matrix(MAC),ncomp=3,Y.add=Ya)
cp <- MVA.synt(plsda1c); cp
# Test (run 4 times, 999 permutations each) :
MVA.test(MAC,Treat,ncomp=3,model="PPLS-DA",cmv=TRUE,Y.add=Ya)
par(mfrow=c(1,1))
MVA.plot(plsda1c,fac=TH3$Treatment,col=c("#ADD6FF","#5A70AD"),points=FALSE,fac.lab="",xax=2,yax=1,fac.cex=0.01)
# Select VOCs contributing most to group separation (highest contribution along barycentres' axis):
PLTH5c <- data.frame(VOC=names(TH5c),plsda1c$projection)
names(PLTH5c)[2:4] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,2],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLTH5c$xprim <- (PLTH5c$NC1*X0[2]+PLTH5c$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLTH5c$xsec <- (PLTH5c$xprim*Xp0[2]+PLTH5c$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLTH5c$Percxx2 <- round(abs(PLTH5c$xsec)/sum(abs(PLTH5c$xsec))*100,digits=2)
PLTH5c[order(PLTH5c$Percxx2),]
t1xx <- which(PLTH5c$Percxx2>5.2); PLTH5c[t1xx,]
length(t1xx); sum(PLTH5c$Percxx2[t1xx]) # 8 VOCs selected,representing 66.3 % of contrib. to axis Xp0
TH6c <- TH5c[,t1xx]

## 1.2 Selected VOCs RO6c final matrix:
MAC <- TH6c; Treat<- TH3$Treatment; Ya <- dummy(TH3$Sex)
plsda1c <- cppls(dummy(TH3$Treatment)~as.matrix(MAC),ncomp=3,Y.add=Ya)
cp <- MVA.synt(plsda1c); cp
MVA.test(MAC,Treat,ncomp=3,model="PPLS-DA",cmv=TRUE,Y.add=Ya)
#save(plsda1c,file="./2021-01-11_TH6c_plsda1c.Rdata")
par(mfrow=c(1,1),las=0,bty="l",mar=c(4,4,1,1))
MVA.plot(plsda1c,fac=TH3$Treatment,col=c("skyblue2","skyblue4"),points=FALSE,xax=2,yax=1,fac.lab="",fac.cex=0.001)
PLTH6c <- data.frame(VOC=names(TH6c),plsda1c$projection)
names(PLTH6c)[2:4] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,2],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLTH6c$xprim <- (PLTH6c$NC1*X0[2]+PLTH6c$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLTH6c$xsec <- (PLTH6c$xprim*Xp0[2]+PLTH6c$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLTH6c$Percxx2 <- round(abs(PLTH6c$xsec)/sum(abs(PLTH6c$xsec))*100,digits=2)
PLTH6c[,6:7] # Contribution of the 10 VOCs to group separation, and distinction of typical "control" VOCs (xsec<0) vs. "dry" VOCs (xsec>0)


#### 2. Total emissions: 
TH3S <- apply(TH3c,1,sum)
TH3D <- apply(TH3c*60,1,function(x)length(which(x>0)))
TH4S <- (TH3S)^(1/4)
plot(TH3S~TH3$Treatment)
t.test(TH4S~TH3$Treatment) # t=-1.01, df=7.53, P=0.34
tapply(TH3S,TH3$Treatment,mean)
tapply(TH3S,TH3$Treatment,se)
t.test(TH3D~TH3$Treatment) # t=0.880 df=8.88 P=0.40
tapply(TH3D,TH3$Treatment,mean)
tapply(TH3D,TH3$Treatment,se)

THFam <- data.frame(Family1=as.factor(c(rep("MONOalc",3),"ARalc",rep("SESQ",4),"FADalc","FADest",
                                        "AR","MONO",rep("MONOalc",3),"SESQoxd","MONO","ARoxd","MONOest","FADalc",
                                        "ARalc","ARoxd","ARald","ARcet","ARest","ARest","MONOest","SESQ")))
cbind(THFam,names(TH3c))
THFam$Family2 <- THFam$Family1
levels(THFam$Family2) <- c("AR",rep("ARox",5),rep("FADox",2),"MONO",rep("MONOox",2),"SESQ",rep("SESQox",2))
THFam$Family3 <- THFam$Family2
levels(THFam$Family3) <- c("AR","AR","FAD","MONO","MONO","SESQ","SESQ")
THFam$Family4 <- THFam$Family3
levels(THFam$Family4) <- c("AR","FAD","TERP","TERP")
THMoy <- data.frame(ARox=c(apply(TH3c[,which(THFam$Family2=="ARox")],1,sum)),
                    MONOox=c(apply(TH3c[,which(THFam$Family2=="MONOox")],1,sum)),
                    SESQsimp=c(apply(TH3c[,which(THFam$Family1=="SESQ")],1,sum)),
                    FADox=c(apply(TH3c[,which(THFam$Family2=="FADox")],1,sum)))
THMoy2 <- THMoy^(1/4)
THMoy2$Ech2 <- TH3$Ech2
THMoy2$Treatment <- TH3$Treatment
plot(THMoy2$ARox~TH3$Treatment)
t.test(THMoy2$ARox~TH3$Treatment) # t=-0.951, df=8.04, P=0.37
tapply(THMoy$ARox,THMoy2$Treatment,mean); tapply(THMoy$ARox,THMoy2$Treatment,se)
plot(THMoy2$MONOox~TH3$Treatment)
t.test(THMoy2$MONOox~TH3$Treatment) # t=-1.30, df=8.57, P=0.23
tapply(THMoy$MONOox,THMoy2$Treatment,mean); tapply(THMoy$MONOox,THMoy2$Treatment,se)
plot(THMoy2$SESQsimp~TH3$Treatment)
t.test(THMoy2$SESQsimp~TH3$Treatment) # t=-0.802, df=7.53, P=0.45
tapply(THMoy$SESQsimp,THMoy2$Treatment,mean); tapply(THMoy$SESQsimp,THMoy2$Treatment,se)
plot(THMoy2$FADox~TH3$Treatment)
t.test(THMoy2$FADox~THMoy2$Treatment) # t=-0.487, df=8.25, P=0.64
tapply(THMoy$FADox,THMoy2$Treatment,mean); tapply(THMoy$FADox,THMoy2$Treatment,se)



#### 3. Leaf VOCs vs. Flower VOCs mutlivariate analysis
MAC <- rbind(TH5c,LTH5c); Treat <- as.factor(c(as.character(TH3$Type),as.character(LTH3$Type)))
plsda1c <- cppls(dummy(Treat)~as.matrix(MAC),ncomp=3,scale=T,center=T)#,Y.add=Ya)
cp <- MVA.synt(plsda1c); cp
MVA.plot(plsda1c,fac=Treat,col=c("#ADD6FF","#5A70AD"),points=FALSE,xax=1,yax=2,fac.lab="",fac.cex=0.001)
# Select VOCs contributing most to group separation (highest contribution along barycentres' axis):
PLTHLc <- data.frame(plsda1c$projection)
names(PLTHLc)[1:3] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,3],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLTHLc$xprim <- (PLTHLc$NC1*X0[2]+PLTHLc$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLTHLc$xsec <- (PLTHLc$xprim*Xp0[2]+PLTHLc$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLTHLc$Percxx2 <- round(abs(PLTHLc$xsec)/sum(abs(PLTHLc$xsec))*100,digits=2)
PLTHLc[order(PLTHLc$Percxx2),]
t1 <- which(PLTHLc$Percxx2>5); PLTHLc[t1,]
length(t1); sum(PLTHLc$Percxx2[t1]) # 13 VOCs selected, representing 74.6% of contrib. to axis
# cca tests :
MAC <- rbind(TH5c,LTH5c)[,t1]; Treat <- as.factor(c(as.character(TH3$Type),as.character(LTH3$Type)))
aaaa=cca(MAC~Treat) # CCA > AFCs (tableau contingence)
anova.cca(aaaa) # df=1,21, F=2.50, P=0.036 *
plsda1c <- cppls(dummy(Treat)~as.matrix(MAC),ncomp=3,scale=T,center=T)#,Y.add=Ya)
cp <- MVA.synt(plsda1c); cp
PLTHLc2 <- data.frame(plsda1c$projection)
names(PLTHLc2)[1:3] <- c("NC1","NC2","NC3")
X0 <- tapply(plsda1c$scores[,1],Treat,mean)
Y0 <- tapply(plsda1c$scores[,2],Treat,mean)
Z0 <- tapply(plsda1c$scores[,3],Treat,mean)
Xprim <- c(as.matrix(plsda1c$scores[,1]*X0[2]+plsda1c$scores[,2]*Y0[2]/sqrt(X0[2]^2+Y0[2]^2))) # change of coordinates (X,Y) > Xprim
Xp0 <- tapply(Xprim,Treat,mean)
PLTHLc2$xprim <- (PLTHLc2$NC1*X0[2]+PLTHLc2$NC2*Y0[2])/sqrt(X0[2]^2+Y0[2]^2)
PLTHLc2$xsec <- (PLTHLc2$xprim*Xp0[2]+PLTHLc2$NC3*Z0[2])/sqrt(Xp0[2]^2+Z0[2]^2) # change of coordinates (Xprim,Z) > Xsec
PLTHLc2$Percxx2 <- round(abs(PLTHLc2$xsec)/sum(abs(PLTHLc2$xsec))*100,digits=2)
PLTHLc2[,5:6] # Contribution of the 17 VOCs to group separation, and distinction of typical leafy VOCs (xsec<0) vs. floral VOCs (xsec>0)

# Total emissions:
LTH3S <- apply(LTH3c,1,sum)
t.test(TH3S^(1/4),LTH3S^(1/4)) # t=1.95, df=6.83, P=0.093 .
mean(TH3S); se(TH3S)
mean(LTH3S); se(LTH3S)
LTH3D <- apply(LTH3c*60,1,function(x)length(which(x>0)))
t.test(TH3D,LTH3D) # t=2.68, df=21.0 P=0.014 *
mean(TH3D); se(TH3D)
mean(LTH3D); se(LTH3D)


######## 2/ Flowering phenology, number of Flowers and number of Fruits --------------------------------------------------------------
### 2.1/ DATA -------------------------------------------------------------
pheno <- read.table("./2020-12-17_data_Pheno_CLIMED.csv",h=T)
pheno$CalendarWeek <- pheno$CalendarDay/7
## Create synthetic table with 1L per plant individual:
pheno2 <- aggregate(pheno$Nflower,by=list(pheno$NumPlant),function(x){sum(na.omit(x))}) # calculate total number of flowers (TotFlo)
names(pheno2) <- c("NumPlant","TotFlo")
pheno2$Species <- as.factor(substr(pheno2$NumPlant,6,7))
pheno2$Treatment <- as.factor(substr(pheno2$NumPlant,1,2))
levels(pheno2$Treatment) <- c("control","drier")
pheno2$Sex <- as.factor(substr(pheno2$NumPlant,9,9))
pheno2$Sex[which(pheno2$Sex=="")] <- NA
# Define max number of fruits, max number of flowers, time of peak flowering, number of observations, Duration of flowering (number of weeks with Nflower > maxFlo/2):
pheno2$Area <- c(); pheno2$MaxFruit <- c(); pheno2$MaxFlo <- c(); pheno2$PeakFlo <- c(); pheno2$NTimes <- c()
pheno2$TimeHalf1 <- c(); pheno2$TimeHalf2 <- c()
for(i in 1:dim(pheno2)[1]){
  t0 <- which(pheno$NumPlant==as.character(pheno2$NumPlant[i]))
  pheno2$Area[i] <- pheno$Area_m2[t0][1]
  pheno2$MaxFruit[i] <- max(na.omit(pheno$Nfruit[t0]))
  pheno2$MaxFlo[i] <- max(na.omit(pheno$Nflower[t0]))
  pheno2$PeakFlo[i] <- pheno$CalendarDay[t0[which.max(pheno$Nflower[t0])]]
  pheno2$NTime[i] <- length(pheno$Nflower[t0])
  t1 <- pheno$CalendarDay[t0[which(na.omit(pheno$Nflower[t0])<=pheno2$MaxFlo[i]/2)]]
  pheno2$Timehalf1[i] <- max(t1[which(t1<pheno2$PeakFlo[i])]) # day at which number of flowers = 1/2 max flo
  pheno2$Timehalf2[i] <- min(t1[which(t1>pheno2$PeakFlo[i])]) # day at which number of flowers = 1/2 max flo
  }
pheno2[which(pheno2$Timehalf1<0),]$Timehalf1[1] <- 6*7 # OK
pheno2[which(pheno2$Timehalf2>1000),]$Timehalf2 <- 151+7
pheno2$TotFlo <- pheno2$TotFlo/pheno2$Area/0.16 # Total number of flowers per m2 ("Area=1" <> 40cm x 40 cm ==0.16m2)
pheno2$MaxFruit <- pheno2$MaxFruit/pheno2$Area/0.16 # max number of fruits per m2
pheno2$PeakFlo <- round(pheno2$PeakFlo/7) # Peak flowering time in week
pheno2$DurationFlo <- round((pheno2$Timehalf2-pheno2$Timehalf1)/7) # Duration of flowering ~ Number of weeks with at least 1/2 max number of flowers
pheno2$Plot <- as.factor(substr(pheno2$NumPlant,1,4))
summary(pheno2)

# Define sub-tables for each plant species:
RO <- pheno2[which(pheno2$Species=="RO"),] # Salvia rosmarinus
CI <- pheno2[which(pheno2$Species=="CI"),] # Cistus albidus
TH <- pheno2[which(pheno2$Species=="TH"),] # Thymus vulgaris

### 2.2/ STATS ----------------------------------------------------------------
par(mfrow=c(1,1),mar=c(4,6,1,1),bty="l",las=1)
## RO
# 1. Total number of flowers
plot(RO$TotFlo~RO$Treatment,ylab="Total number of flowers\\ncounted throughout the season and per m2",col="mediumpurple",xlab="")
mro <- glmer.nb(round(TotFlo)~Treatment+(1|Plot),data=RO)
summary(mro)
anova(mro,update(mro,~.-Treatment)) # df=1 Chisq=0.0948 P=0.76
t0 <- simulateResiduals(mro); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

# 2. Time of flowering peak
plot(RO$PeakFlo~RO$Treatm,ylab="Time of observed flowering peak\\n(calendar weeks)",col="mediumpurple",xlab="")
wilcox.test(RO$PeakFlo~RO$Treatm) # W=209 P=0.82

# 3. Duration of flowering
plot(RO$DurationFlo~RO$Treatment,ylab="Number of weeks with at least 1/2 of max. number of flowers",col="mediumpurple",xlab="")
wilcox.test(RO$DurationFlo~RO$Treatment) # W=197 P=0.93

# 4. Number of fruits
plot(RO$MaxFruit~RO$Treatment,ylab="Number of fruits",col="mediumpurple",xlab="") # W=147 P=0.15
RO$Plot <- as.factor(substr(RO$NumPlant,1,4))
m0 <- glmer.nb(as.integer(MaxFruit)~Treatment+(1|Plot),data=RO)
summary(m0)
anova(m0,update(m0,~.-Treatment),test="Chisq") # df=1 Chisq=0.369 P=0.54
t0 <- simulateResiduals(m0); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok


## CI
# 1. Total number of flowers
plot(CI$TotFlo~CI$Treatment,ylab="Total number of flowers\\ncounted throughout the season and per m2",col="deeppink",xlab="")
mCI <- glmer.nb(round(TotFlo)~Treatment+(1|Plot),data=CI)
summary(mci)
anova(mCI,update(mCI,~.-Treatment)) # df=1 Chisq=0.171 P=0.68
t0 <- simulateResiduals(mCI); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

# 2. Time of flowering peak
plot(CI$PeakFlo~CI$Treatm,ylab="Time of observed flowering peak\\n(calendar weeks)",col="deeppink",xlab="")
wilcox.test(CI$PeakFlo~CI$Treatm) # W=170 P=0.35

# 3. Duration of flowering
plot(CI$DurationFlo~CI$Treatm,ylab="Number of weeks with at least 1/2 of max. number of flowers",col="deeppink",xlab="")
wilcox.test(CI$DurationFlo~CI$Treatm) # W=157 P=0.17

# 4. Number of fruits
plot(CI$MaxFruit~CI$Treatment,ylab="Number of fruits",col="deeppink",xlab="") # W=147 P=0.15
m1 <- glmer.nb(as.integer(MaxFruit)~Treatment+(1|Plot),data=CI)
summary(m1)
anova(m1,update(m1,~.-Treatment),test="Chisq") # df=1 Chisq=0.902 P=0.34
t0 <- simulateResiduals(m1); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok


## TH
# 1. Total number of flowers
plot(TH$TotFlo~TH$Treatment,ylab="Total number of flowers\\ncounted throughout the season and per m2",col="skyblue2",xlab="")
TH$Sex <- as.factor(as.character(TH$Sex))
mth1 <- glmer.nb(round(TotFlo)~Treatment+Sex+(1|Plot),data=TH)
summary(mth1)
anova(mth1,update(mth1,~.-Treatment)) # df=1 Chisq=2.98 P=0.084 .
t0 <- simulateResiduals(mth1); plot(t0) # deviation OK
testDispersion(t0) # NO dispersion ok

# 2. Time of flowering peak
plot(TH$PeakFlo~TH$Treatment,ylab="Time of observed flowering peak\\n(calendar weeks)",col="Skyblue2",xlab="")
wilcox.test(TH$PeakFlo~TH$Treatm) # W=57.5 P=0.54

# 3. Duration of flowering
plot(TH$DurationFlo~TH$Treatment,ylab="Number of weeks with at least 1/2 of max. number of flowers",col="Skyblue2",xlab="")
wilcox.test(TH$DurationFlo[which(TH$DurationFlo<10)]~TH$Treatment[which(TH$DurationFlo<10)]) # W=41.5 P=0.86


### 2.4/ Figure Pheno --------------------------------------------------------------
RO2 <- pheno[which(pheno$PlantSpecies=="RO"),]
CI2 <- pheno[which(pheno$PlantSpecies=="CI"),]
TH2 <- pheno[which(pheno$PlantSpecies=="TH"),]

## RO2 flowering phenology
RO2$CalendarWeek <- round(RO2$CalendarDay/7)
RO2C <- RO2[which(RO2$Treatment=="control"),]
RO2D <- RO2[which(RO2$Treatment=="drier"),]
ROCM <- tapply(RO2C$Nflower/RO2C$Area/0.16,RO2C$CalendarWeek,mean)
ROCSE <- tapply(RO2C$Nflower/RO2C$Area/0.16,RO2C$CalendarWeek,se)
RODM <- tapply(RO2D$Nflower/RO2D$Area/0.16,RO2D$CalendarWeek,mean)
RODSE <- tapply(RO2D$Nflower/RO2D$Area/0.16,RO2D$CalendarWeek,se)
ROday <- as.integer(levels(as.factor(RO2$CalendarWeek)))

## CI2 flowering phenology
CI2$CalendarWeek <- round(CI2$CalendarDay/7)
CI2C <- CI2[which(CI2$Treatment=="control"),]
CI2D <- CI2[which(CI2$Treatment=="drier"),]
CICM <- tapply(CI2C$Nflower/CI2C$Area/0.16,CI2C$CalendarWeek,mean)
CICSE <- tapply(CI2C$Nflower/CI2C$Area/0.16,CI2C$CalendarWeek,se)
CIDM <- tapply(CI2D$Nflower/CI2D$Area/0.16,CI2D$CalendarWeek,mean)
CIDSE <- tapply(CI2D$Nflower/CI2D$Area/0.16,CI2D$CalendarWeek,se)
CIday <- as.integer(levels(as.factor(CI2$CalendarWeek)))

## TH2 flowering phenology
TH2$CalendarWeek <- round(TH2$CalendarDay/7)
TH2C <- TH2[which(TH2$Treatment=="control"),]
TH2D <- TH2[which(TH2$Treatment=="drier"),]
THCM <- tapply(TH2C$Nflower/TH2C$Area/0.16,TH2C$CalendarWeek,function(x)mean(na.omit(x)))
THCSE <- tapply(TH2C$Nflower/TH2C$Area/0.16,TH2C$CalendarWeek,function(x)se(na.omit(x)))
THDM <- tapply(TH2D$Nflower/TH2D$Area/0.16,TH2D$CalendarWeek,function(x)mean(na.omit(x)))
THDSE <- tapply(TH2D$Nflower/TH2D$Area/0.16,TH2D$CalendarWeek,function(x)se(na.omit(x)))
THday <- as.integer(levels(as.factor(TH2$CalendarWeek)))

## FIG
par(mfrow=c(3,1),bty="l",las=1,mar=c(5,4,1,1))
plot(ROCM~ROday,col="#B69BDE",lty=2,pch=21,type="o",ylim=c(0,2500),
     ylab="Mean number of flowers per m2",xlab="Calendar weeks",cex=1.2,xlim=c(6,22))
arrows(ROday,ROCM-ROCSE,ROday,ROCM+ROCSE,code=3,length=0.05,angle=90,col="#B69BDE")
points(RODM~ROday,col="mediumpurple4",type="o",lty=1,pch=19,cex=0.7)
arrows(ROday,RODM-RODSE,ROday,RODM+RODSE,code=3,length=0.05,angle=90,col="mediumpurple4")
text(5.6,2500,labels=expression(bold("(A)")),adj=0)
text(6.5,2500,labels=expression(italic("S. rosmarinus")),adj=0)
legend(5.6,2500,legend=c("control","drought"),col=c("#B69BDE","mediumpurple4"),lty=c(2,1),pch=c(21,19),
       pt.cex=c(1.2,0.7),bty="n",x.intersp=0.1,y.intersp=0.3,adj=0)
plot(CICM~CIday,col="#DEB2CC",type="o",lty=2,pch=21,ylim=c(0,100),ylab="Mean number of flowers per m2",
     xlab="Calendar weeks",cex=1.2,xlim=c(6,22))
arrows(CIday,CICM-CICSE,CIday,CICM+CICSE,code=3,length=0.05,angle=90,col="#DEB2CC")
points(CIDM~CIday,col="#AB538E",type="o",lty=1,pch=19,cex=0.7)
arrows(CIday,CIDM-CIDSE,CIday,CIDM+CIDSE,code=3,length=0.05,angle=90,col="#AB538E")
text(5.6,100,labels=expression(bold("(B)")),adj=0)
text(6.5,100,labels=expression(italic("C. albidus")),adj=0)
legend(5.6,100,legend=c("control","drought"),col=c("#DEB2CC","#AB538E"),lty=c(2,1),pch=c(21,19),pt.cex=c(1.2,0.7),
       bty="n",x.intersp=0.1,y.intersp=0.3)
plot(THCM~THday,col="skyblue2",type="o",lty=2,pch=21,ylim=c(0,16000),ylab="Mean number of flowers per m2",
     xlab="Calendar weeks",cex=1.2,xlim=c(6,22))
arrows(THday,THCM-THCSE,THday,THCM+THCSE,code=3,length=0.05,angle=90,col="skyblue2")
points(THDM~THday,col="skyblue4",type="o",lty=1,pch=19,cex=0.7)
arrows(THday,THDM-THDSE,THday,THDM+THDSE,code=3,length=0.05,angle=90,col="skyblue4")
text(5.6,16000,labels=expression(bold("(C)")),adj=0)
text(6.5,16000,labels=expression(italic("T. vulgaris")),adj=0)
legend(5.6,16000,legend=c("control","drought"),col=c("skyblue2","skyblue4"),lty=c(2,1),pch=c(21,19),pt.cex=c(1.2,0.7),
       bty="n",x.intersp=0.1,y.intersp=0.3)

######## 3/ NECTAR Production, Flower size--------------------------------------------------------------------------------------
### 3.1/ DATA -----------------------------------
NEC <- read.csv("2018-10-21_Nectar_Data_Climed.csv") # Raw data, nectar volumes and concentrations calculated
NEC$Sex <- NEC$Species
levels(NEC$Sex) <- c(NA,NA,"F","H")
levels(NEC$Species)[3:4] <- "TH"
summary(NEC)

RO <- NEC[which(NEC$Species=="RO"),]
RO$PlanteNum <- as.factor(as.character(RO$PlanteNum))
RO$Plot <- as.factor(as.character(RO$Plot))
summary(RO)
CI <- NEC[which(NEC$Species=="CI"),]
CI$PlanteNum <- as.factor(as.character(CI$PlanteNum))
CI$Plot <- as.factor(as.character(CI$Plot))
summary(CI)
TH <- NEC[which(NEC$Species=="TH"),]
TH$PlanteNum <- as.factor(as.character(TH$PlanteNum))
summary(TH)


### 3.2/ Rosemary RO ---------------------------------------------------
# 3.2.1/ Proportion of flowers with nectar
# Define dataframe with 1L per plant individual:
ROfl <- data.frame(Plant=levels(RO$PlanteNum),Treatment=substr(levels(RO$PlanteNum),2,2),
                   Plot=substr(levels(RO$PlanteNum),1,4))
levels(ROfl$Treatment) <- c("control","drier")
ROfl$Sampled <- tapply(RO$NectarVoluL,RO$PlanteNum,function(x){length(x)}) # Number of flowers sampled
ROfl$PropNec <- tapply(RO$NectarVoluL,RO$PlanteNum,function(x){length(which(x>0))/length(x)}) # Proportion of flowers with nectar (nectar volume>0)
ROfl$Full <- tapply(RO$NectarVoluL,RO$PlanteNum,function(x){length(which(x>0))}) # number of flowers with nectar
ROfl$Empty <- ROfl$Sampled - ROfl$Full # Number of flowers without nectar
summary(ROfl)
plot(ROfl$PropNec~ROfl$Treatment,ylim=c(0,1))

mro1 <- glmer(cbind(Empty,Full)~Treatment+(1|Plot),data=ROfl,family=binomial(link="logit"))
summary(mro1)
anova(mro1,update(mro1,~.-Treatment),test="Chisq") # df=1 Chisq=0.0442 P=0.83
t0 <- simulateResiduals(mro1); plot(t0) # OK
testDispersion(t0) # OK

# 3.2.2/ Volume in nL
RONec <- RO[which(RO$NectarVoluL>0),]
summary(RONec)
plot(RONec$NectarVoluL~RONec$Treatment)

mro2 <- glmer(round(NectarVoluL*100)~Treatment+(1|Plot/PlanteNum),data=RONec,family=poisson(link="log"))
summary(mro2)
anova(mro2,update(mro2,~.-Treatment)) # df=1 Chisq=0.171 P=0.68
t0 <- simulateResiduals(mro2); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok


# 3.2.3/ Sugar content
mro3 <- glmer.nb(round(Sugar_ug)~Treatment+(1|Plot/PlanteNum),data=RONec)
summary(mro3)
mro3b <- update(mro3,~.-Treatment)
anova(mro3,mro3b) # Chisq=0.0431 df=1 P=0.84
t0 <- simulateResiduals(mro3b); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

# 3.2.4/ Flower size
plot(RO$Size1~RO$Treatment)
mro <- lmer(Size1~Treatment+(1|Plot/PlanteNum),data=RO)
anova(mro,update(mro,~.-Treatment)) # df=1 Chisq=0.144 P=0.70
t0 <- simulateResiduals(mro); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

### 3.3/ Cistus CI ---------------------------------------------------
# 3.3.1/ Proportion of flowers with nectar
# Define dataframe with 1L per plant individual:
CIfl <- data.frame(Plant=levels(CI$PlanteNum),Treatment=substr(levels(CI$PlanteNum),2,2),
                   Plot=substr(levels(CI$PlanteNum),1,4))
levels(CIfl$Treatment) <- c("control","drier")
CIfl$PropNec <- tapply(CI$NectarVoluL,CI$PlanteNum,function(x){length(which(x>0))/length(x)}) # Proportion of flowers with nectar (nectar volume>0)
CIfl$Sampled <- tapply(CI$NectarVoluL,CI$PlanteNum,function(x){length(x)}) # Number of flowers sampled
CIfl$Full <- tapply(CI$NectarVoluL,CI$PlanteNum,function(x){length(which(x>0))}) # number of flowers with nectar
CIfl$Empty <- CIfl$Sampled - CIfl$Full # Number of flowers without nectar
summary(CIfl)

mci2 <- glmer(cbind(Full,Empty)~Treatment+(1|Plot),data=CIfl,family=binomial(link="logit"))
summary(mci2)
anova(mci2,update(mci2,~.-Treatment),test="Chisq") # df=1 Chisq=1.49 P=0.22
t0 <- simulateResiduals(mci2); plot(t0) # OK
testDispersion(t0) ~ OK

# 3.3.2/ Volume in nL
CINec <- CI[which(CI$NectarVoluL>0),]
CINec$PlanteNum <- as.factor(as.character(CINec$PlanteNum))

mci2 <- glmer(round(NectarVoluL*100)~Treatment+(1|Plot/PlanteNum),data=CINec,family=poisson(link="log"))
summary(mci2)
anova(mci2,update(mci2,~.-Treatment)) # df=1 Chisq=0.804 P=0.37
t0 <- simulateResiduals(mci2); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

# 3.3.3/ Sugar content
mci3 <- glmer.nb(round(Sugar_ug)~Treatment+(1|Plot/PlanteNum),data=CINec[which(is.na(CINec$Sugar_ug)==F),])
summary(mci3)
mci3b <- update(mci3,~.-Treatment)
anova(mci3,mci3b) # Chisq=0.730 df=1 P=0.40
t0 <- simulateResiduals(mci3); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok


# 3.3.4/ Flower size
par(mfrow=c(1,1))
CI$MSizePink <- apply(cbind(CI$Size1,CI$Size3),1,function(x){mean(na.omit(x))})
CI$SESizePink <- apply(cbind(CI$Size1,CI$Size3),1,function(x){sd(na.omit(x))/length(sqrt(na.omit(x)))})
CI$MSizeYellow <- apply(cbind(CI$Size2,CI$Size4),1,function(x){mean(na.omit(x))})
CI$SESizeYellow <- apply(cbind(CI$Size2,CI$Size4),1,function(x){sd(na.omit(x))/length(sqrt(na.omit(x)))})
summary(CI)
par(mfrow=c(2,2),mar=c(4,4,1,1),bty="l",las=1)
plot(CI$MSizePink~CI$Treatment,col="magenta")
plot(CI$MSizeYellow~CI$Treatment,col="gold")
plot(CI$SESizePink~CI$Treatment,col="magenta")
plot(CI$SESizeYellow~CI$Treatment,col="gold")

mci1 <- lmer(MSizePink~Treatment+(1|Plot/PlanteNum),data=CI)
anova(mci1,update(mci1,~.-Treatment)) # Chisq=0.812 df=1 P=0.37
t0 <- simulateResiduals(mci1); plot(t0) # OK
testDispersion(t0) # OK
mci2 <- lmer(MSizeYellow~Treatment+(1|Plot/PlanteNum),data=CI)
anova(mci2,update(mci2,~.-Treatment)) # Chisq=1.77 df=1 P=0.18
t0 <- simulateResiduals(mci2); plot(t0) # OK
testDispersion(t0) # OK


### 3.4/ Thyme TH ---------------------------------------------------
# 3.4.1/ Proportion of flowers with nectar
# Define dataframe with 1L per plant individual:
THfl <- data.frame(Plant=levels(TH$PlanteNum),Treatment=substr(levels(TH$PlanteNum),2,2),
                   Plot=substr(levels(TH$PlanteNum),1,4))
levels(THfl$Treatment) <- c("control","drier")
THfl$Sex <- NA; for(i in 1:length(THfl$Plant)){THfl$Sex[i] <- as.character(TH$Sex[which(TH$PlanteNum==as.character(THfl$Plant[i]))][1])}
THfl$Sex <- as.factor(THfl$Sex)
THfl$PropNec <- tapply(TH$NectarVoluL,TH$PlanteNum,function(x){length(which(x>0))/length(x)}) # Proportion of flowers with nectar (nectar volume>0)
THfl$Sampled <- tapply(TH$NectarVoluL,TH$PlanteNum,function(x){length(x)}) # Number of flowers sampled
THfl$Full <- tapply(TH$NectarVoluL,TH$PlanteNum,function(x){length(which(x>0))}) # number of flowers with nectar
THfl$Empty <- THfl$Sampled - THfl$Full # Number of flowers without nectar
summary(THfl)

tapply(THfl$PropNec,THfl$Treat,summary)
par(mfrow=c(1,1)); plot(THfl$PropNec~THfl$Treat,ylim=c(0,1))

mth2 <- glmer(cbind(Full,Empty)~Treatment+Sex+(1|Plot),data=THfl,family=binomial(link="logit"))
summary(mth2)
anova(mth2,update(mth2,~.-Treatment),test="Chisq") # df=1 Chisq=3.58 P=0.059
t0 <- simulateResiduals(mth2); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok
emmeans(mth2,specs=pairwise~Treatment)
tapply(THfl$PropNec,THfl$Treatment,summary)

# 3.4.2/ Volume
THNec <- TH[which(TH$NectarVoluL>0),]
plot(THNec$NectarVoluL~THNec$Treatment)

mth2 <- glmer(round(NectarVoluL*100)~Treatment+Sex+(1|Plot/PlanteNum),data=THNec,family=poisson(link="log"))
summary(mth2)
anova(mth2,update(mth2,~.-Treatment)) # df=1 Chisq=1.13 df=1 P=0.29
t0 <- simulateResiduals(mth2); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

# 3.4.3/ Sugar content
mth3 <- glmer.nb(round(Sugar_ug)~Treatment+Sex+(1|Plot/PlanteNum),data=THNec)
summary(mth3)
anova(mth3,update(mth3,~.-Treatment)) # Chisq=6.19 df=1 P=0.013 *
t0 <- simulateResiduals(mth3); plot(t0) # KS not great
testDispersion(t0) # NO dispersion ok ~~~~
tapply(THNec$Sugar_ug,THNec$Treatment,function(x){mean(na.omit(x))})
tapply(THNec$Sugar_ug,THNec$Treatment,function(x){se(na.omit(x))})

# 3.4.4/ Flower size
par(mfrow=c(1,1))
plot(TH$Size1~TH$Treatment)

mth1 <- lmer(Size1~Treatment+Sex+(1|Plot/PlanteNum),data=TH)
anova(mth1,update(mth1,~.-Treatment)) # Chisq=0.0815 df=1 P=0.78
t0 <- simulateResiduals(mth1); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

TH$Diam <- apply(cbind(TH$Size2,TH$Size3),1,function(x){mean(na.omit(x))})
mth3 <- lmer(Diam~Treatment+Sex+(1|Plot/PlanteNum),data=TH)
anova(mth3,update(mth3,~.-Treatment)) # Chisq=0.545 df=1 P=0.46
t0 <- simulateResiduals(mth3); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

### 3.5/ Figure -----------------------------------------------------------------
levels(ROfl$Treatment)[2] <- "drought"
levels(CIfl$Treatment)[2] <- "drought"
levels(THfl$Treatment)[2] <- "drought"
levels(RONec$Treatment)[2] <- "drought"
levels(CINec$Treatment)[2] <- "drought"
levels(THNec$Treatment)[2] <- "drought"

par(mfrow=c(3,3),las=1,mar=c(4,4,3,1),bty="l",xpd=TRUE)

# RO
plot(ROfl$PropNec~ROfl$Treat,xlab=NULL,ylab="Proportion of flowers with nectar",col="#B69BDE",ylim=c(0,1),outline=F,border="grey40")
r0 <- rnorm(27,0,0.1); t0 <- c(1,2)[ROfl$Treatment]+r0
points(ROfl$PropNec~t0,pch=19,cex=0.6)
text(0.1,1.3,labels=expression(bold("(A)")),adj=0)
text(0.5,1.3,labels=expression(italic("S. rosmarinus")),adj=0)
plot(RONec$NectarVoluL~RONec$Treatment,xlab=NULL,ylab="Volume of nectar (µL per flower)",col="#B69BDE",outline=F,ylim=c(0,0.8),border="grey40")
r0 <- rnorm(109,0,0.1); t0 <- c(1,2)[RONec$Treatment]+r0
points(RONec$NectarVoluL~t0,pch=19,cex=0.6)
plot(RONec$Sugar_ug~RONec$Treatment,xlab=NULL,ylab="Sugar content (µg per flower)",col="#B69BDE",outline=F,ylim=c(0,300),border="grey40")
points(RONec$Sugar_ug~t0,pch=19,cex=0.6)
# CI
plot(CIfl$PropNec~CIfl$Treat,xlab=NULL,ylab="Proportion of flowers with nectar",col="#DEB2CC",ylim=c(0,1),outline=F,border="grey40")
r0 <- rnorm(24,0,0.1); t0 <- c(1,2)[CIfl$Treatment]+r0
points(CIfl$PropNec~t0,pch=19,cex=0.6)
text(0.1,1.3,labels=expression(bold("(B)")),adj=0)
text(0.5,1.3,labels=expression(italic("C. albidus")),adj=0)
plot(CINec$NectarVoluL~CINec$Treatment,xlab=NULL,ylab="Volume of nectar (µL per flower)",col="#DEB2CC",outline=F,ylim=c(0,1.5),border="grey30")
r0 <- rnorm(68,0,0.1); t0 <- c(1,2)[CINec$Treatment]+r0
points(CINec$NectarVoluL~t0,pch=19,cex=0.6)
plot(CINec$Sugar_ug~CINec$Treatment,xlab=NULL,ylab="Quantity of sugar (µg per flower)",col="#DEB2CC",outline=F,ylim=c(0,350),border="grey30")
points(CINec$Sugar_ug~t0,pch=19,cex=0.6)
# TH
plot(THfl$PropNec~THfl$Treat,xlab=NULL,ylab="Proportion of flowers with nectar",col="skyblue2",ylim=c(0,1),border="grey40")
r0 <- rnorm(19,0,0.1); t0 <- c(1,2)[THfl$Treatment]+r0
points(THfl$PropNec~t0,pch=19,cex=0.6)
text(0.1,1.3,labels=expression(bold("(C)")),adj=0)
text(0.5,1.3,labels=expression(italic("T. vulgaris")),adj=0)
arrows(c(1,1,2),1.15,c(2,1,2),c(1.15,1.1,1),code=0,lwd=c(2,1,1))
text(1.5,1.2,labels=expression(bold(italic("P = 0.059"))),cex=0.9)
plot(THNec$NectarVoluL~THNec$Treatment,xlab=NULL,ylab="Volume of nectar (µL per flower)",col="skyblue2",outline=F,ylim=c(0,0.12),border="grey40")
r0 <- rnorm(60,0,0.1); t0 <- c(1,2)[THNec$Treatment]+r0
points(THNec$NectarVoluL~t0,pch=19,cex=0.6)
plot(THNec$Sugar_ug~THNec$Treatment,xlab=NULL,ylab="Quantity of sugar (µg per flower)",col="skyblue2",outline=F,ylim=c(0,90),border="grey40")
points(THNec$Sugar_ug~t0,pch=19,cex=0.6)
arrows(c(1,1,2),95,c(2,1,2),c(90,85,75)+5,code=0,lwd=c(2,1,1))
text(1.5,100,labels=expression(bold("*")),cex=1.1)


######## 4/ Flower colour ------------------------------------------------------------------------------------------------------
### 4.1/ DATA -----------------------------------------------------------
Spectro <- read.csv("./2019-03-22_Spectro_colour_data.csv")
summary(Spectro)

# Define species-specific dataframes:
RO <- Spectro[which(Spectro$Species=="RO"),] # rosemary
RO$Plot <- as.factor(as.character(RO$Plot))
RO$PlanteNum <- as.factor(as.character(RO$PlanteNum))
CI <- Spectro[which(Spectro$Species=="CI"),] # cistus
CI$Plot <- as.factor(as.character(CI$Plot))
CI$PlanteNum <- as.factor(as.character(CI$PlanteNum))
TH <- Spectro[which(Spectro$Species=="THH" | Spectro$Species=="THF"),] # Thyme
TH$Sex <- as.factor(as.character(TH$Species))
TH$Plot <- as.factor(as.character(TH$Plot))
TH$PlanteNum <- as.factor(as.character(TH$PlanteNum))


### 4.2/ STATS ------------------------------------------------------------------
## RO
par(mfrow=c(1,2),mar=c(4,4,1,1),bty="l",las=1,xpd=F)
plot(RO$R1~RO$Treatment,col="mediumpurple")
plot(RO$Theta1~RO$Treatment,col="mediumpurple")

mr <- lmer(R1~Treatment+(1|Plot/PlanteNum),data=RO)
summary(mr)
anova(mr,update(mr,~.-Treatment)) # Chisq=0.415 df=1 P=0.52
t0 <- simulateResiduals(mr); plot(t0) # KS not great
testDispersion(t0) # OK
mt <- lmer(Theta1~Treatment+(1|Plot/PlanteNum),data=RO)
anova(mt,update(mt,~.-Treatment)) # Chisq=0.458 df=1 P=0.50
t0 <- simulateResiduals(mt); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

## CI
par(mfrow=c(2,2),mar=c(4,4,1,1),bty="l",las=1)
plot(CI$R1~CI$Treatment,col="deeppink")
plot(CI$Theta1~CI$Treatment,col="deeppink")
plot(CI$R2~CI$Treatment,col="gold")
plot(CI$Theta2~CI$Treatment,col="gold")

mr1 <- lmer(R1~Treatment+(1|Plot/PlanteNum),data=CI)
anova(mr1,update(mr1,~.-Treatment)) # Chisq=0.275 df=1 P=0.60
t0 <- simulateResiduals(mr1); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok
mr2 <- lmer(R2~Treatment+(1|Plot/PlanteNum),data=CI)
summary(mr2)
anova(mr2,update(mr2,~.-Treatment)) # Chisq=0.426 df=1 P=0.51
t0 <- simulateResiduals(mr2); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

mt1 <- lmer(Theta1~Treatment+(1|Plot/PlanteNum),data=CI)
anova(mt1,update(mt1,~.-Treatment)) # Chisq=1.36 df=1 P=0.24
t0 <- simulateResiduals(mt1); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok
mt2 <- lmer(Theta2~Treatment+(1|Plot/PlanteNum),data=CI)
anova(mt2,update(mt2,~.-Treatment)) # Chisq=3.40 df=1 P=0.065 .
t0 <- simulateResiduals(mt2); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok
tapply(CI$Theta2,CI$Treatment,mean)
tapply(CI$Theta2,CI$Treatment,se)

## TH
par(mfrow=c(1,2),mar=c(4,4,1,1),bty="l",las=1)
plot(TH$R1~TH$Treatment,col="skyblue2")
plot(TH$Theta1~TH$Treatment,col="skyblue2")

mr <- lmer(R1~Treatment+Sex+(1|Plot/PlanteNum),data=TH)
anova(mr,update(mr,~.-Treatment)) # Chisq=0.627 df=1 P=0.43
t0 <- simulateResiduals(mr); plot(t0) # KS borderline...
testDispersion(t0) # NO dispersion ok

mt <- lmer(Theta1~Treatment+Sex+(1|Plot/PlanteNum),data=TH)
anova(mt,update(mt,~.-Treatment)) # df=1 Chisq=0.0025 P=0.96
t0 <- simulateResiduals(mt); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

### 4.3/ Figure --------------------------------------------------------------
x0=c(-sqrt(3)/2,sqrt(3)/2,-sqrt(3)/2,-sqrt(3)/2,0,0)*2/3;y0=c(-1/2,-1/2,1/2,-1/2,1,-1)*2/3 # hexagon, x-coordinates
x1=c(-sqrt(3)/2,sqrt(3)/2,0,0,sqrt(3)/2,sqrt(3)/2)*2/3;y1=c(1/2,1/2,1,-1,1/2,-1/2)*2/3 # hexagon, y-coordinates
par(mfrow=c(1,1),bty='n',xpd=TRUE,pty="s",mar=c(0.6,0.9,3,0.9),xaxt="n",yaxt="n",family="Arial")

Spectro$Colour <- as.factor(paste(Spectro$Species,Spectro$Treatment,sep="_"))
levels(Spectro$Colour) <- c("deeppink1","deeppink4","mediumpurple1","mediumpurple4",
                            "grey60","grey10","cyan1","cyan4")
Spectro$Colour <- as.character(Spectro$Colour)  
RO <- Spectro[which(Spectro$Species=="RO"),]
CI <- Spectro[which(Spectro$Species=="CI"),]
TH <- Spectro[which(Spectro$Species=="THH" | Spectro$Species=="THF"),]

par(mfrow=c(1,1))
## RO
plot(RO$X1,RO$Y1,xlim=c(-0.6,0.6),ylim=c(-0.6,0.6),xlab="",ylab="",
     pch=4,col=RO$Colour) # valeurs moyennes
arrows(x0,y0,x1,y1,length=0,lwd=1.5) # draws hexagon
text(c(-sqrt(3)/2-0.1,sqrt(3)/2+0.1,0)*2/3,c(-1/2-0.1,-1/2-0.1,1+0.08)*2/3,labels=c("E(U)","E(G)","E(B)")) # labels 3 récepteurs sensoriels
arrows(c(0,0,0),c(0,0,0),c(0,-sqrt(3)/2,sqrt(3)/2)*2/3,c(1,-1/2,-1/2)*2/3,lty=1,length=0.1,lwd=1)
legend(-0.7,0.75,legend=c("control","drought"),pch=4,cex=1,
       col=c("mediumpurple1","mediumpurple4"),
       bty="n",x.intersp=0.5,y.intersp=0.7)
text(-0.75,0.8,labels=expression(bold("(A)")),cex=1,adj=0)
text(-0.65,0.8,labels=expression(italic("S. rosmarinus")),cex=1,adj=0)
## CI
plot(CI$X1,CI$Y1,xlim=c(-0.6,0.6),ylim=c(-0.6,0.6),xlab="",ylab="",
     pch=4,col=CI$Colour) # valeurs moyennes
points(CI$X2,CI$Y2,
       pch=4,col=ifelse(Spectro$Treatment=="control","gold","orange2")) # valeurs moyennes
arrows(x0,y0,x1,y1,length=0,lwd=1.5) # draws hexagon
text(c(-sqrt(3)/2-0.1,sqrt(3)/2+0.1,0)*2/3,c(-1/2-0.1,-1/2-0.1,1+0.08)*2/3,labels=c("E(U)","E(G)","E(B)")) # labels 3 récepteurs sensoriels
arrows(c(0,0,0),c(0,0,0),c(0,-sqrt(3)/2,sqrt(3)/2)*2/3,c(1,-1/2,-1/2)*2/3,lty=1,length=0.1,lwd=1)
legend(-0.7,0.65,legend=c("","","",""),pch=4,cex=1,ncol=2,col=c("deeppink1","deeppink4","gold","orange2"),
       bty="n",x.intersp=0.3,y.intersp=0.7)
text(-0.5,c(0.52,0.58),labels=c("control","drought"),adj=0)
text(c(-0.6,-0.38),0.65,labels=c("tip","centre"),adj=c(1,0))
text(-0.75,0.8,labels=expression(bold("(B)")),cex=1,adj=0)
text(-0.65,0.8,labels=expression(italic("C. albidus")),cex=1,adj=0)
## TH
plot(TH$X1,TH$Y1,xlim=c(-0.6,0.6),ylim=c(-0.6,0.6),xlab="",ylab="",
     pch=4,col=TH$Colour) # valeurs moyennes
arrows(x0,y0,x1,y1,length=0,lwd=1.5) # draws hexagon
text(c(-sqrt(3)/2-0.1,sqrt(3)/2+0.1,0)*2/3,c(-1/2-0.1,-1/2-0.1,1+0.08)*2/3,labels=c("E(U)","E(G)","E(B)")) # labels 3 récepteurs sensoriels
arrows(c(0,0,0),c(0,0,0),c(0,-sqrt(3)/2,sqrt(3)/2)*2/3,c(1,-1/2,-1/2)*2/3,lty=1,length=0.1,lwd=1)
legend(-0.7,0.65,legend=c("","","",""),pch=4,cex=1,ncol=2,col=c("grey60","grey10","cyan1","cyan4"),
       bty="n",x.intersp=0.3,y.intersp=0.7)
text(-0.5,c(0.52,0.58),labels=c("control","drought"),adj=0)
text(c(-0.65,-0.57),0.65,labels=c("F","H"),adj=0)
text(-0.75,0.8,labels=expression(bold("(C)")),cex=1,adj=0)
text(-0.65,0.8,labels=expression(italic("T. vulgaris")),cex=1,adj=0)




# Save CLIM_colours_3species.eps 800*600

#points(CLIMRCT$X1[CLIMRCT$Species=="THF" & CLIMRCT$Treatment=="exclu"],
#      CLIMRCT$Y1[CLIMRCT$Species=="THF" & CLIMRCT$Treatment=="exclu"],
#       pch=4,col="grey20") # valeurs moyennes



######################## B/ POLLINATOR VISITS ------------------------------------------------
######## 5/ Pollinator visits - morphogroups -----------------------------------------------------------------------------------
### 5.1/ DATA -----------------------------------------------------------------
CLIM0 <- read.csv("2021-02-04_Pollinator_observations.csv") # 1L per Plant species, Plot and Date
summary(CLIM0)
CLIM0$Plot_Plant <- as.factor(paste(CLIM0$Plot,CLIM0$Plant,sep="_"))

## Exclude thyme (not enough visits)
CLIM1 <- CLIM0[which(CLIM0$Plant=="RO" | CLIM0$Plant=="CI"),]
CLIM1$Plant <- as.factor(as.character(CLIM1$Plant))
summary(CLIM1)

## CLIM2 Data frame 1L per Plot and Plant species:
CLIM2 <- aggregate(CLIM0[,-c(1:5,7,14)],by=list(CLIM0$Plot_Plant),sum)
names(CLIM2)[1] <- "Plot_Plant"
CLIM2$Plot <- as.factor(substr(as.character(CLIM2$Plot_Plant),1,4))
CLIM2$Plant <- as.factor(substr(as.character(CLIM2$Plot_Plant),6,7))
CLIM2$Sex <- as.factor(substr(as.character(CLIM2$Plot_Plant),8,8))
CLIM2$Treatment <- as.factor(substr(as.character(CLIM2$Plot),1,2))
levels(CLIM2$Treatment) <- c("control","drier")
CLIM2$NVis <- apply(CLIM2[,3:8],1,sum)
summary(CLIM2)

## Reorganize (transpose) CLIM1 data frame:
CLIM1t <- data.frame(Treatment=rep(CLIM1$Treatment,6),Plot=rep(CLIM1$Plot,6),Date=rep(CLIM1$Date,6),
                     Plant=rep(CLIM1$Plant,6),Flowers=rep(CLIM1$Flowers,6),Time=rep(CLIM1$Time,6),
                     Morphog=c(rep("SBee",dim(CLIM1)[1]),rep("BBee",dim(CLIM1)[1]),rep("Apis",dim(CLIM1)[1]),
                               rep("BumB",dim(CLIM1)[1]),rep("Coleo",dim(CLIM1)[1]),rep("Diptera",dim(CLIM1)[1])),
                     NVis=c(CLIM1$SBee,CLIM1$BBee,CLIM1$Apis,CLIM1$BumB,CLIM1$Coleo,CLIM1$Diptera))
summary(CLIM1t)

### 5.2/ Check that flower number and phenology at the plot level are not affected by drought ----------------------------------------------------------
for(i in 1:dim(CLIM2)[1]){
  t0 <- which(CLIM0$Plot_Plant==as.character(CLIM2$Plot_Plant[i]))
  CLIM2$PeakFlo[i] <- CLIM0$CalWeek[t0[which.max(CLIM0$Flowers[t0])]]
  CLIM2$NTime[i] <- length(CLIM0$Flowers[t0])
  t1 <- CLIM0$CalWeek[t0[which(CLIM0$Flowers[t0]<=max(CLIM0$Flowers[t0]/2))]]
  CLIM2$Timehalf1[i] <- max(t1[which(t1<CLIM2$PeakFlo[i])]) # day at which number of flowers = 1/2 max flo
  CLIM2$Timehalf2[i] <- min(t1[which(t1>CLIM2$PeakFlo[i])]) # day at which number of flowers = 1/2 max flo
}
CLIM2[which(CLIM2$Timehalf1<0),]$Timehalf1 <- CLIM2[which(CLIM2$Timehalf1<0),]$PeakFlo-1
CLIM2[which(CLIM2$Timehalf2>100),]$Timehalf2 <- CLIM2[which(CLIM2$Timehalf2>100),]$PeakFlo+1
CLIM2$DurationFlo <- CLIM2$Timehalf2-CLIM2$Timehalf1 # Duration of flowering ~ Number of weeks with at least 1/2 max number of flowers
summary(CLIM2)

RO <- CLIM2[which(CLIM2$Plant=="RO"),]
CI <- CLIM2[which(CLIM2$Plant=="CI"),]
TH <- CLIM2[which(CLIM2$Plant=="TH"),]

### 1. Total flowers
par(mfrow=c(1,1),mar=c(4,4,1,1),bty="l",las=1,xpd=F)
## Total flowers RO
plot(Flowers~Treatment,data=RO,col="mediumpurple")
mro <- lm(round(Flowers)~Treatment,data=RO)
summary(mro)
anova(mro,update(mro,~.-Treatment)) # df=1 F=0.0016 P=0.97
t0 <- simulateResiduals(mro); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

## Total flowers CI
plot(Flowers~Treatment,data=CI,col="deeppink")
mci <- lm(round(Flowers)~Treatment,data=CI)
summary(mci)
anova(mci,update(mci,~.-Treatment)) # df=1 F=0.136 P=0.72
t0 <- simulateResiduals(mci); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

## Total flowers TH
plot(Flowers~Treatment,data=TH,col="skyblue2")
mth <- lm(round(Flowers)~Treatment*Sex,data=TH)
mth1 <- lm(round(Flowers)~Treatment+Sex,data=TH)
summary(mth1)
anova(mth,mth1) # df=1 F=0.397 Chisq=0.544 
anova(mth1,update(mth1,~.-Treatment)) # df=1 F=4.91 P=0.044 *
anova(mth1,update(mth1,~.-Sex)) # df=1 F=0.756 P=0.40
t0 <- simulateResiduals(mth1); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok


### 2. Flowering Peak
## RO
plot(PeakFlo~Treatment,data=RO,col="mediumpurple")
wilcox.test(PeakFlo~Treatment,data=RO) # W=48.5 P=0.93
## CI
plot(PeakFlo~Treatment,data=CI,col="deeppink")
wilcox.test(PeakFlo~Treatment,data=CI) # W=45 P=0.58
## TH
plot(PeakFlo~Treatment,data=TH,col="skyblue2")
wilcox.test(PeakFlo~Treatment,data=TH) # W=22 P=0.94

### 3. Flowering Duration
## RO
plot(DurationFlo~Treatment,data=RO,col="mediumpurple")
wilcox.test(DurationFlo~Treatment,data=RO) # W=45.5 P=0.76
## CI
plot(DurationFlo~Treatment,data=CI,col="deeppink")
wilcox.test(DurationFlo~Treatment,data=CI) # W=43 P=0.54
## TH
plot(DurationFlo~Treatment,data=TH,col="skyblue2")
wilcox.test(DurationFlo~Treatment,data=TH) # W=19 P=0.82


### 5.3/ Total number of visits per plot and plant species -----------------------------------------
par(mfrow=c(1,2),mar=c(4,5,2,1),bty="l")
plot(NVis~Treatment,data=CLIM2[which(CLIM2$Plant=="RO"),],ylab="Total number of visits per plot",xlab="",col="mediumpurple")
points(NVis~Treatment,data=CLIM2[which(CLIM2$Plant=="RO"),])
plot(NVis~Treatment,data=CLIM2[which(CLIM2$Plant=="CI"),],ylab="Total number of visits per plot",xlab="",col="deeppink")
points(NVis~Treatment,data=CLIM2[which(CLIM2$Plant=="CI"),])

## STATS
mv <- glmer(NVis~Treatment*Plant+(1|Plot),data=CLIM2[which(CLIM2$Plant!="TH"),],family=poisson(link="log"))
summary(mv)
mv2 <- update(mv,~.-Treatment:Plant)
summary(mv2)
anova(mv,mv2) # df=2 Chisq=0.668 P=0.41
anova(mv2,update(mv2,~.-Treatment)) # df=1,Chisq=1.24 P=0.26
t4 <- simulateResiduals(mv2); plot(t4) # deviation NS
testDispersion(t4) # dispersion OK

### 5.3.0/ Figure Total number of visits --------------------------------
levels(CLIM2$Treatment)[2] <- "drought"
par(mfrow=c(1,2),mar=c(4,5,2,1),bty="l",las=1)
r0 <- rnorm(20,0,0.1); t0 <- c(rep(1,10),rep(2,10))+r0
plot(NVis~Treatment,data=CLIM2[which(CLIM2$Plant=="RO"),],ylim=c(0,950),ylab="Total number of visits per plot",xlab="",col="#B69BDE")
points(NVis~t0,data=CLIM2[which(CLIM2$Plant=="RO"),],pch=19,cex=0.8)
text(0.5,960,labels=expression(bold("(A)")),adj=0)
text(0.8,960,labels=expression(italic("S. rosmarinus")),adj=0)
plot(NVis~Treatment,data=CLIM2[which(CLIM2$Plant=="CI"),],outline=F,ylim=c(0,53),ylab="Total number of visits per plot",xlab="",col="#DEB2CC")
points(NVis~t0,data=CLIM2[which(CLIM2$Plant=="CI"),],pch=19,cex=0.8)
text(0.5,53.5,labels=expression(bold("(B)")),adj=0)
text(0.8,53.5,labels=expression(italic("C. albidus")),adj=0)


### 5.4/ Week-to-week visits by pollinator morphogroup - Hernel model ------------------------------------------------------------
## 5.3.1/ Presence/absence of visits
CLIM1t$NVis2 <- CLIM1t$NVis
CLIM1t$NVis2[which(CLIM1t$NVis2>0)] <- 1
summary(CLIM1t)
m0a <- glmer.nb(NVis~Treatment*Morphog*Plant+scale(Flowers)+scale(Time)+(1|Plot),data=CLIM1t)
summary(m0a)
anova(m0a,update(m0a,~.-Treatment:Morphog:Plant)) # df=5 hisq=6.876 P=0.23
t0 <- simulateResiduals(m0a); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok


m1a <- glmer(NVis2~Treatment*Morphog*Plant+scale(Flowers)+scale(Time)+(1|Plot),family=binomial(link="logit"),data=CLIM1t)
summary(m1a)
m1b <- glmer(NVis2~(Treatment+Morphog+Plant)^2+scale(Flowers)+scale(Time)+(1|Plot),family=binomial(link="logit"),data=CLIM1t)
summary(m1b)
anova(m1a,m1b) # df=5 Chisq=1.83 P=0.87
m1c <- glmer(NVis2~(Treatment+Morphog)*Plant+scale(Flowers)+scale(Time)+(1|Plot),family=binomial(link="logit"),data=CLIM1t)
summary(m1c)
anova(m1b,m1c) # df=5 Chisq=5.43 P=0.37
m1d <- glmer(NVis2~Treatment+Morphog*Plant+scale(Flowers)+scale(Time)+(1|Plot),family=binomial(link="logit"),data=CLIM1t)
summary(m1d)
anova(m1c,m1d) # df=5 Chisq=2.39 P=0.12
anova(m1d,update(m1d,~.-Morphog:Plant)) # df=5, Chisq=156, P<0.001 ***
m1e <- glmer(NVis2~Morphog*Plant+scale(Flowers)+scale(Time)+(1|Plot),family=binomial(link="logit"),data=CLIM1t)
anova(m1d,m1e) # df=1, Chisq=1.11, P=0.29
t0 <- simulateResiduals(m1e); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok

# 5.3.2/ Intensiy of presence = Number of visits
CLIM1ta <- CLIM1t[which(CLIM1t$NVis>0),]
m2a <- glmer.nb(NVis~Treatment*Morphog*Plant+scale(Flowers)+scale(Time)+(1|Plot),data=CLIM1ta)
summary(m2a)
m2b <- glmer.nb(NVis~(Treatment+Morphog+Plant)^2+scale(Flowers)+scale(Time)+(1|Plot),data=CLIM1ta)
summary(m2b)
anova(m2a,m2b) # df=4 Chisq=10.3 P=0.035 *
emmeans(m2a,specs=pairwise ~Treatment | Morphog*Plant,type="response")
## SBeeRO C<D P=0.0015 ** / BumBRO C>D P=0.016 * / ApisRO C>D P=0.011 *
t0 <- simulateResiduals(m2a); plot(t0) # deviation NS
testDispersion(t0) # NO dispersion ok
#save(m2a, file="2021-02-04_PollMorph_m2a.Rdata")


### 5.4.0/ Figure ggplot2 -----------------------------------------------------------
library(ggplot2)
## Barplot
CLIM0$Treatment_Plant <- as.factor(paste(CLIM0$Treatment,CLIM0$Plant,sep="_"))
levels(CLIM0$Treatment_Plant)[c(3,7,4,8)] <- c("control_TH","drier_TH","control_TH","drier_TH")
CLIM3m <- aggregate(CLIM0[,-c(1:5,7,14:15)],by=list(CLIM0$Treatment_Plant),mean)
CLIM3sum <- aggregate(CLIM0[,-c(1:5,7,14:15)],by=list(CLIM0$Treatment_Plant),sum)
CLIM3se <- aggregate(CLIM0[,-c(1:5,7,14:15)],by=list(CLIM0$Treatment_Plant),se)
names(CLIM3m)[1] <- "Treatment_Plant"
names(CLIM3sum)[1] <- "Treatment_Plant"
names(CLIM3se)[1] <- "Treatment_Plant"
CLIM3m$Treatment <- strsplit(as.character(CLIM3m$Treatment_Plant),"_")[[1]][1]
CLIM3m$Plant <- strsplit(as.character(CLIM3m$Treatment_Plant),"_")[[1]][1]
for(i in 1:dim(CLIM3m)[1]){
  CLIM3m$Treatment[i] <- strsplit(as.character(CLIM3m$Treatment_Plant[i]),"_")[[1]][1]
  CLIM3m$Plant[i] <- strsplit(as.character(CLIM3m$Treatment_Plant[i]),"_")[[1]][2]
}
CLIM3m$Treatment <- as.factor(CLIM3m$Treatment)
CLIM3m$Plant <- as.factor(CLIM3m$Plant)

t0 <- order(CLIM3m$Plant)
CLIM3m <- CLIM3m[t0,]
CLIM3se <- CLIM3se[t0,]
CLIM3m <- CLIM3m[,c(5:3,6:8,1,2,9:10)]
CLIM3se <- CLIM3se[,c(5:3,6:8,1,2)]

CLIMb <- data.frame(M0=c(c(as.matrix(CLIM3m[c(3:4),1:6])),c(as.matrix(CLIM3m[c(1:2),1:6])),c(as.matrix(CLIM3m[c(5:6),1:6]))),
                    SE0=c(c(as.matrix(CLIM3se[c(3:4),1:6])),c(as.matrix(CLIM3se[c(1:2),1:6])),c(as.matrix(CLIM3se[c(5:6),1:6]))),
                    Treatment_Plant=CLIM3m$Treatment_Plant[c(rep(3:4,6),rep(1:2,6),rep(5:6,6))],
                    Treatment=CLIM3m$Treatment[c(rep(3:4,6),rep(1:2,6),rep(5:6,6))],
                    Plant=CLIM3m$Plant[c(rep(3:4,6),rep(1:2,6),rep(5:6,6))],
                    Morphog=names(CLIM3m)[sort(rep(1:6,2))])
CLIMb$Treatment_Plant_Morphog <- as.factor(paste(CLIMb$Treatment_Plant,CLIMb$Morphog,sep="_"))
CLIMb$TwoSide <- sort(rep(c(1:(36/2)),2))

CLIMb2 <- CLIMb
CLIMb2[1:2,1:2] <- CLIMb2[1:2,1:2]/2

CLIMc <- CLIMb[1:24,]
aa <- (CLIMc$M0-CLIMc$SE0)*c(-1,1)
ab <- (CLIMc$M0+CLIMc$SE0)*c(-1,1)
ac <- CLIMc$TwoSide
levels(CLIMc$Treatment)[2] <- "drought"

p <- ggplot(data=CLIMc, aes(x=TwoSide, y=M0, fill=Treatment_Plant)) +
  scale_y_continuous(labels=abs,breaks=seq(-50,45,5),trans="reverse",position = "right") +
  labs(x="",y="Mean number of visits per plot, day and plant species\\n") +
  #  theme(plot.margin = unit(c(1,2,1,1), "lines"),panel.grid.major.y=element_line(colour="grey")) +
  geom_vline(xintercept=seq(7,48/3,6)-0.5,colour="grey",size=0.2) +
  geom_hline(yintercept=0,colour="grey30",size=0.2) +
  scale_fill_manual(values=c("#DEB2CC","#B69BDE","#AB538E","mediumpurple4"),guide=FALSE) +
 # geom_segment(x=CLIMb$TwoSide,y=(CLIMb$M0-CLIMb$SE0)*c(-1,1),
  #             xend=CLIMb$TwoSide,yend=(CLIMb$M0+CLIMb$SE0)*c(-1,1),size=20) +
  geom_segment(x=-ac,y=aa,
               xend=-ac,yend=ab,size=1) +
  theme(axis.text=element_text(size=10,colour="black"),
        #axis.text.x=element_text(hjust=1),
        axis.title=element_text(size=12,family="Helvetica"),#,margin = margin(t = 0, r = 20, b = 0, l = 40)),
        axis.line=element_line(size=0.1,colour="black"),
        panel.background=element_rect(fill=NA,colour=NA),
       # legend.position=c(0.95,0.85),legend.background=element_rect(fill=NA,colour=NA),
        legend.key=element_rect(fill=NA,colour=NA),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank()) +
  annotate("text",x=c(1:12),y=54,size=4,hjust=0,
           label=rep(c("honeybees","large wild bees","small wild bees","bumblebees","Coleoptera","Diptera"),2)) +
  annotate("text",x=c(1,4,3)+0.2,y=c(-25,-0.5,-4),size=5,hjust=0,
           label=c("*","*","**")) +
  annotate("text",x=13,y=c(20,-15),size=5,parse=T,
           label=c(expression(bold("control")),expression(bold("drought")))) +
  annotate("text",x=c(3.5,3.5+6)+0.2,y=-10,size=4,hjust=0,parse=F,
           label=c(expression(italic("S. rosmarinus")),expression(italic("C. albidus")))) +
  geom_bar(data=subset(CLIMc,Treatment=="control"),stat="identity",size=0.5,width=0.9,na.rm=FALSE) +
  geom_bar(data=subset(CLIMc,Treatment=="drought"),stat="identity",position="identity",size=0.5,width=0.9,
           mapping=aes(y=-M0))# +
p + coord_flip() + scale_x_reverse()

### 5.4.0/ Figure ggplot2 RO/CI >0 -----------------------------------------------------------
library(ggplot2)
## Barplot
CLIM1ta$Treatment_Plant <- as.factor(paste(CLIM1ta$Treatment,CLIM1ta$Plant,sep="_"))
CLIM3m <- data.frame(Apis=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="Apis")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="Apis")],mean),
                     BBee=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="BBee")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="BBee")],mean),
                     BumB=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="BumB")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="BumB")],mean),
                     SBee=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="SBee")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="SBee")],mean),
                     Coleo=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="Coleo")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="Coleo")],mean),
                     Diptera=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="Diptera")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="Diptera")],mean))
CLIM3size <- data.frame(Apis=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="Apis")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="Apis")],length),
                     BBee=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="BBee")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="BBee")],length),
                     BumB=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="BumB")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="BumB")],length),
                     SBee=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="SBee")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="SBee")],length),
                     Coleo=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="Coleo")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="Coleo")],length),
                     Diptera=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="Diptera")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="Diptera")],length))
CLIM3se <- data.frame(Apis=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="Apis")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="Apis")],se),
                     BBee=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="BBee")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="BBee")],se),
                     BumB=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="BumB")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="BumB")],se),
                     SBee=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="SBee")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="SBee")],se),
                     Coleo=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="Coleo")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="Coleo")],se),
                     Diptera=tapply(CLIM1ta$NVis[which(CLIM1ta$Morphog=="Diptera")],CLIM1ta$Treatment_Plant[which(CLIM1ta$Morphog=="Diptera")],se))
CLIM3m$Treatment_Plant <- row.names(CLIM3m)
CLIM3m$Treatment <- as.factor(c("control","control","drought","drought"))
CLIM3m$Plant <- as.factor(c("CI","RO","CI","RO"))

t0 <- order(CLIM3m$Plant)
CLIM3m <- CLIM3m[t0,]
CLIM3se <- CLIM3se[t0,]
CLIM3size <- CLIM3size[t0,]

CLIM3m <- CLIM3m[,c(1:2,4,3,5:6,8:9)]
CLIM3se <- CLIM3se[,c(1:2,4,3,5:6,8:9)]
CLIM3se$Coleo[which(is.na(CLIM3se$Coleo))] <- 0
CLIM3se$Diptera[which(is.na(CLIM3se$Diptera))] <- 0
CLIM3size$Diptera[which(is.na(CLIM3size$Diptera))] <- 0

CLIMb <- data.frame(M0=c(c(as.matrix(CLIM3m[c(3:4),1:6])),c(as.matrix(CLIM3m[c(1:2),1:6]))),
                    SE0=c(c(as.matrix(CLIM3se[c(3:4),1:6])),c(as.matrix(CLIM3se[c(1:2),1:6]))),
                    SIZE=c(c(as.matrix(CLIM3size[c(3:4),1:6])),c(as.matrix(CLIM3size[c(1:2),1:6]))),
                    Treatment=CLIM3m$Treatment[c(rep(3:4,6),rep(1:2,6))],
                    Plant=CLIM3m$Plant[c(rep(3:4,6),rep(1:2,6))],
                    Morphog=names(CLIM3m)[sort(rep(1:6,2))])
CLIMb$Treatment_Plant <- as.factor(paste(CLIMb$Treatment,CLIMb$Plant,sep="_"))
CLIMb$Treatment_Plant_Morphog <- as.factor(paste(CLIMb$Treatment_Plant,CLIMb$Morphog,sep="_"))
CLIMb$TwoSide <- sort(rep(c(1:(24/2)),2))

CLIMc <- CLIMb[1:24,]
aa <- (CLIMc$M0-CLIMc$SE0)*c(-1,1)
ab <- (CLIMc$M0+CLIMc$SE0)*c(-1,1)
ac <- CLIMc$TwoSide
levels(CLIMc$Treatment)[2] <- "drought"

p <- ggplot(data=CLIMc, aes(x=TwoSide, y=M0, fill=Treatment_Plant)) +
  scale_y_continuous(labels=abs,breaks=seq(-50,80,10),trans="reverse",position = "right") +
  labs(x="",y="Mean number of visits >0 per plot, day and plant species\\n") +
  #  theme(plot.margin = unit(c(1,2,1,1), "lines"),panel.grid.major.y=element_line(colour="grey")) +
  geom_vline(xintercept=seq(7,48/3,6)-0.5,colour="grey",size=0.2) +
  geom_hline(yintercept=0,colour="grey30",size=0.2) +
  scale_fill_manual(values=c("#DEB2CC","#B69BDE","#AB538E","mediumpurple4"),guide=FALSE) +
  # geom_segment(x=CLIMb$TwoSide,y=(CLIMb$M0-CLIMb$SE0)*c(-1,1),
  #             xend=CLIMb$TwoSide,yend=(CLIMb$M0+CLIMb$SE0)*c(-1,1),size=20) +
  geom_segment(x=-ac,y=aa,
               xend=-ac,yend=ab,size=1) +
  theme(axis.text=element_text(size=10,colour="black"),
        #axis.text.x=element_text(hjust=1),
        axis.title=element_text(size=12,family="Helvetica"),#,margin = margin(t = 0, r = 20, b = 0, l = 40)),
        axis.line=element_line(size=0.1,colour="black"),
        panel.background=element_rect(fill=NA,colour=NA),
        # legend.position=c(0.95,0.85),legend.background=element_rect(fill=NA,colour=NA),
        legend.key=element_rect(fill=NA,colour=NA),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank()) +
  annotate("text",x=c(1:12),y=115,size=4,hjust=0,
           label=rep(c("honeybees","large wild bees","small wild bees","bumblebees","Coleoptera","Diptera"),2)) +
  annotate("text",x=c(1:12),y=(CLIMc$M0[which(CLIMc$Treatment=="control")]+CLIMc$SE0[which(CLIMc$Treatment=="control")]+5),size=3,hjust=0,
           label=CLIMc$SIZE[which(CLIMc$Treatment=="control")],adj=0) +
  annotate("text",x=c(1:12),y=(CLIMc$M0[which(CLIMc$Treatment=="drought")]+CLIMc$SE0[which(CLIMc$Treatment=="drought")]+2)*-1,size=3,hjust=0,
           label=CLIMc$SIZE[which(CLIMc$Treatment=="drought")],adj=0) +
  annotate("text",x=c(1,4,3)+0.1,y=c(-55,-17,-17),size=5,hjust=0,
           label=c("*","*","**")) +
  annotate("text",x=13,y=c(20,-20),size=5,parse=T,
           label=c(expression(bold("control")),expression(bold("drought")))) +
  annotate("text",x=c(3.5,3.5+6)+0.2,y=-30,size=4,hjust=0,parse=F,
           label=c(expression(italic("S. rosmarinus")),expression(italic("C. albidus")))) +
  geom_bar(data=subset(CLIMc,Treatment=="control"),stat="identity",size=0.5,width=0.9,na.rm=FALSE) +
  geom_bar(data=subset(CLIMc,Treatment=="drought"),stat="identity",position="identity",size=0.5,width=0.9,
           mapping=aes(y=-M0))# +
p + coord_flip() + scale_x_reverse()

######## 6/ Pollinator captures -----------------------------------------------------------------------------------
### 6.1/ DATA ----------------------------------------------------------------
Cap <- read.csv("./2021-02-07_Pollinator_captures_identification.csv")
summary(Cap)

### 6.2/ PCAs : Community composition ------------------------------------------------------------
# Remove singletons, and select only specimens identified to the species level
CAPCCA <- Cap[,c(3:38)[-c(13,28,31)]][which(apply(Cap[,c(3:38)[-c(13,28,31)]],2,sum)>1)]
aaaa=cca(CAPCCA~Cap$Treatment) # CCA > AFCs (tableau contingence)
anova.cca(aaaa) # F=0.670 df=1,18 P=0.91
plot(aaaa)

### 6.3/ Abundance & species richness ####### 
# In Abundance : select all
# In species richness : select only specimens identified to the species level
RIC <- data.frame(Plot=Cap$Plot,Treatment=Cap$Treatment,ABUND=apply(Cap[,c(3:38)],1,sum),
                  RICHNESS=apply(Cap[,c(3:38)[-c(13,28,31)]],1,function(x)length(which(x>0))))
RIC$ABUNDNOApis <- RIC$ABUND - Cap$Apis_mellifera

## Figure:
levels(RIC$Treatment)[2] <- "drought"
par(mfrow=c(1,2),bty="l",las=1,mar=c(3,4,1,1))
r0 <- rnorm(20,0,0.1); t0 <- c(rep(1,10),rep(2,10))+r0
plot(RIC$ABUND~RIC$Treat,ylab="Number of captures",xlab="",ylim=c(0,30),col="grey")
points(RIC$ABUND~t0,pch=19,cex=0.8)
text(0.7,30.4,labels=expression(bold("(A)")))
plot(RIC$RICHNESS~RIC$Treat,ylab="Species richness",xlab="",ylim=c(0,12),col="grey")
points(RIC$RICHNESS~t0,pch=19,cex=0.8)
text(0.7,12.1,labels=expression(bold("(B)")))

plot(RIC$ABUNDNOApis~RIC$Treat,ylab="Number of captures without honeybee speciments",xlab="",ylim=c(0,30),col="grey")
points(RIC$ABUNDNOApis~t0,pch=19,cex=0.8)


t.test(RIC$ABUND~RIC$Treat) # t=0.738 df=17.7 P=0.47
NORM(lm(RIC$ABUND~RIC$Treat)) # P=0.79 OK

t.test(RIC$ABUNDNOApis~RIC$Treat) # t=1.10 df=18.0 P=0.29
NORM(lm(RIC$ABUNDNOApis~RIC$Treat)) # P=0.088 OK

t.test(RIC$RICHNESS~RIC$Treat) # t=1.04 df=17.3 P=0.31
NORM(lm(RIC$RICHNESS~RIC$Treat)) # OK


######################## C/ PLANT REPRODUCTION ---------------------------------------------------
### 7.1/ Number of seeds per fruit --------------------------------------------
seeds <- read.csv("./2020-12-18_Seeds_per_fruit.csv") # Table number of seeds per fruit
seeds$Treatment <- as.factor(substr(seeds$PlantID,1,2))
levels(seeds$Treatment) <- c("control", "drier")
seeds$Plant <- as.factor(substr(seeds$PlantID,6,7))
summary(seeds)

# Define summary table with 1L per plant individual:
SEED <- data.frame(PlantID=levels(seeds$PlantID),SEEDM=tapply(seeds$NSeeds,seeds$PlantID,mean),
                   SEEDVAR=tapply(seeds$NSeeds,seeds$PlantID,function(x){var(x)}))
SEED$Treatment <- as.factor(substr(SEED$PlantID,1,2))
levels(SEED$Treatment) <- c("control","drier")
SEED$Species <- as.factor(substr(SEED$PlantID,6,7))
SEED$Plot <- as.factor(substr(SEED$PlantID,1,4))
summary(SEED)
RO <- SEED[SEED$Species=="RO",]
CI <- SEED[SEED$Species=="CI",]

#### STATS
## RO
plot(RO$Treatment,RO$SEEDM,col="mediumpurple",ylab="mean number of seeds per fruits")
plot(RO$Treatment,RO$SEEDVAR,col="mediumpurple",ylab="Variance in number of seeds per fruits")

# Mean seed number
mro <- lmer(SEEDM~Treatment+(1|Plot),data=RO) 
summary(mro)
anova(mro,update(mro,~.-Treatment)) # df=1,38 Chisq=1.54 P=0.22
t0 <- simulateResiduals(mro); plot(t0) # deviation NS
testDispersion(t0) # OK
# Variance in seed number
mro <- glmer.nb(round(SEEDVAR*10)~Treatment+(1|Plot),data=RO) 
summary(mro)
anova(mro,update(mro,~.-Treatment)) # df=1 Chisq=1.84 P=0.17
t0 <- simulateResiduals(mro); plot(t0) # deviation NS
testDispersion(t0) # OK


## CI
plot(CI$Treatment,CI$SEEDM,col="deeppink",ylab="Mean number of seeds per fruits")
plot(CI$Treatment,CI$SEEDVAR,col="deeppink",ylab="Variance in number of seeds per fruits")

# Mean seed number
mci <- lmer(SEEDM~Treatment+(1|Plot),data=CI) 
summary(mci)
anova(mci,update(mci,~.-Treatment)) # df=1?38 Chisq=0.0408 P=0.84
t0 <- simulateResiduals(mci); plot(t0) # deviation NS
testDispersion(t0) # OK

# Variance in seed number
mci <- glmer.nb(round(SEEDVAR*10)~Treatment+(1|Plot),data=CI)
summary(mci)
anova(mci,update(mci,~.-Treatment)) # df=1 CHisq=1.23 P=0.27
t0 <- simulateResiduals(mci); plot(t0) # deviation NS
testDispersion(t0) # OK


### 7.2/ Mean seed mass -------------------------------------------------------
seeds2 <- read.csv("./2020-12-18_Seed_mass.csv")
summary(seeds2)
RO <- seeds2[which(seeds2$Species=="RO"),]
CI <- seeds2[which(seeds2$Species=="CI"),]

## RO 
plot(MeanSeedMass~Treatment,data=RO,col="mediumpurple")
mro <- lmer(MeanSeedMass~Treatment+(1|Plot),data=RO)
summary(mro)
anova(mro,update(mro,~.-Treatment)) # df=1 Chisq=0.0423 P=0.84
t0 <- simulateResiduals(mro); plot(t0) # OK
testDispersion(t0) # dispersion ok

## CI 
plot(MeanSeedMass~Treatment,data=seeds2[which(seeds2$Species=="CI"),],)
plot(MeanSeedMass~Treatment,data=CI,col="deeppink")
mci <- lmer(MeanSeedMass~Treatment+(1|Plot),data=CI)
summary(mci)
anova(mci,update(mci,~.-Treatment)) # df=1 Chisq=0.194 P=0.66
t0 <- simulateResiduals(mci); plot(t0) # OK
testDispersion(t0) # dispersion ok


### 7.0/ Figure 1 ------------------------------------------------------------------
par(mfrow=c(2,4),las=1,bty="l",mar=c(3,5,2,2),xpd=TRUE)
levels(RO$Treatment)[2] <- "drought"
levels(CI$Treatment)[2] <- "drought"
levels(SEED$Treatment)[2] <- "drought"
# 1/ RO
r0 <- rnorm(40,0,0.1); t0 <- c(rep(1,20),rep(2,20))+r0
plot(RO$MaxFruit/1000~RO$Treatment,ylab="Number of fruits per m2 (/1000)",col="#B69BDE",xlab="",ylim=c(0,30),outline=F)
points(RO$MaxFruit/1000~t0,pch=19,cex=0.8)
text(0.5,32,label=expression(bold("(A)")),cex=1.1,adj=0)
text(0.9,32,label=expression(italic("S. rosmarinus")),cex=1.1,adj=0)
plot(SEED$Treatment[SEED$Species=="RO"],SEED$SEEDM[SEED$Species=="RO"],outline=F,col="#B69BDE",ylab="mean number of seeds\\nper fruits",xlab="",ylim=c(0,4))
points(t0,SEED$SEEDM[SEED$Species=="RO"],pch=19,cex=0.8)
plot(SEED$Treatment[SEED$Species=="RO"],SEED$SEEDVAR[SEED$Species=="RO"],outline=F,col="#B69BDE",ylab="Variance in number of seeds\\nper fruits",xlab="",ylim=c(0,5))
points(t0,SEED$SEEDVAR[SEED$Species=="RO"],pch=19,cex=0.8)
plot(seeds2$MeanSeedMass[seeds2$Species=="RO"]~seeds2$Treatment[seeds2$Species=="RO"],col="#B69BDE",outline=F,
        ylab="Mean seed mass (µg)",xlab="",ylim=c(0,1))
points(t0,seeds2$MeanSeedMass[seeds2$Species=="RO"],pch=19,cex=0.8)

# 2/ CI
plot(CI$MaxFruit~RO$Treatment,ylab="Number of fruits per m2",col="#DEB2CC",xlab="",outline=F)
points(CI$MaxFruit~t0,pch=19,cex=0.8)
text(0.5,920,label=expression(bold("(B)")),cex=1.1,adj=0)
text(0.9,920,label=expression(italic("C. albidus")),cex=1.1,adj=0)
plot(SEED$Treatment[SEED$Species=="CI"],SEED$SEEDM[SEED$Species=="CI"],outline=F,col="#DEB2CC",ylab="mean number of seeds\\nper fruits",xlab="")
points(t0,SEED$SEEDM[SEED$Species=="CI"],pch=19,cex=0.8)
plot(SEED$Treatment[SEED$Species=="CI"],SEED$SEEDVAR[SEED$Species=="CI"],outline=F,col="#DEB2CC",ylab="Variance in number of seeds\\nper fruits",xlab="")
points(t0,SEED$SEEDVAR[SEED$Species=="CI"],pch=19,cex=0.8)
plot(seeds2$MeanSeedMass[seeds2$Species=="CI"]~seeds2$Treatment[seeds2$Species=="CI"],col="#DEB2CC",outline=F,
     ylab="Mean seed mass (µg)",xlab="")
points(t0,seeds2$MeanSeedMass[seeds2$Species=="CI"],pch=19,cex=0.8)


### 7.3/ Seed experiment ------------------------------------------------------
### DATA
Seedexp <- read.csv("./2021-02-01_Seed_experiment_data.csv",h=T)
summary(Seedexp)
# Define rosemary and cistus tables :
RO <- Seedexp[which(Seedexp$Plant=="RO"),]
RO$PlantID <- as.factor(as.character(RO$PlantID))
CI <- Seedexp[which(Seedexp$Plant=="CI"),]
CI$PlantID <- as.factor(as.character(CI$PlantID))


### 1/ Viability rate
# RO
mvro <- glmer(cbind(Viable,Unviable)~MatTreat*SeedTreat+(1|PlantID),family=binomial(link="logit"),data=RO)
summary(mvro)
mvro1 <- glmer(cbind(Viable,Unviable)~MatTreat+SeedTreat+(1|PlantID),family=binomial(link="logit"),data=RO)
summary(mvro1)
anova(mvro,mvro1) # df=1 Chisq=0.266 P=0.61
anova(mvro1,update(mvro1,~.-MatTreat)) # df=1 Chisq=0.0031 P=0.96
anova(mvro1,update(mvro1,~.-SeedTreat)) # df=1 Chisq=1.74 P=0.19
t0 <- simulateResiduals(mvro1); plot(t0) # OK
testDispersion(t0) # OK
# CI
mvci <- glmer(cbind(Viable,Unviable)~MatTreat*SeedTreat+(1|PlantID),family=binomial(link="logit"),data=CI)
summary(mvci)
mvci1 <- glmer(cbind(Viable,Unviable)~MatTreat+SeedTreat+(1|PlantID),family=binomial(link="logit"),data=CI)
summary(mvci1)
anova(mvci,mvci1) # df=1 Chisq=0.20 P=0.66
anova(mvci1,update(mvci1,~.-MatTreat)) # df=1 Chisq=3.79 P=0.052 .
anova(mvci1,update(mvci1,~.-SeedTreat)) # df=1 Chisq=84.8 P<0.001 ***
t0 <- simulateResiduals(mvci1); plot(t0) # OK
testDispersion(t0) # OK


### 2/ Germination rate
# RO
mgro <- glmer(cbind(Germ,ViableNotGerm)~MatTreat*SeedTreat+(1|PlantID),family=binomial(link="logit"),data=RO)
summary(mgro)
mgro1 <- glmer(cbind(Germ,ViableNotGerm)~MatTreat+SeedTreat+(1|PlantID),family=binomial(link="logit"),data=RO)
summary(mgro1)
anova(mgro,mgro1) # df=1 Chisq=1.51 P=0.22
anova(mgro1,update(mgro1,~.-MatTreat)) # df=1 Chisq=0.854 P=0.36
anova(mgro1,update(mgro1,~.-SeedTreat)) # df=1 Chisq=2.70 P=0.10
t0 <- simulateResiduals(mgro1); plot(t0) # OK
testDispersion(t0) # OK

mgci <- glmer(cbind(Germ,ViableNotGerm)~MatTreat*SeedTreat+(1|PlantID),family=binomial(link="logit"),data=CI)
summary(mgci)
mgci1 <- glmer(cbind(Germ,ViableNotGerm)~MatTreat+SeedTreat+(1|PlantID),family=binomial(link="logit"),data=CI)
summary(mgci1)
anova(mgci,mgci1) # df=1 Chisq=1.19 P=0.28
anova(mgci1,update(mgci1,~.-MatTreat)) # df=1 Chisq=3.70 P=0.054 .
anova(mgci1,update(mgci1,~.-SeedTreat)) # df=1 Chisq=47.8 P<0.001 ***
emmeans(mgci1,specs=pairwise ~ MatTreat | SeedTreat,type="response") # C>D P=0.049 *
t0 <- simulateResiduals(mgci1); plot(t0) # OK
testDispersion(t0) # OK


### 3/ Survival rate
# RO
msro <- glmer(cbind(GermDead,GermSurv)~MatTreat*SeedTreat+(1|PlantID),family=binomial(link="logit"),data=RO)
summary(msro)
msro1 <- glmer(cbind(GermDead,GermSurv)~MatTreat+SeedTreat+(1|PlantID),family=binomial(link="logit"),data=RO)
summary(msro1)
anova(msro,msro1) # df=1 Chisq=2.59 P=0.11
anova(msro1,update(msro1,~.-MatTreat)) # df=1 Chisq=0.0214 P=0.88
anova(msro1,update(msro1,~.-SeedTreat)) # df=1 Chisq=0.149 P=0.70
t0 <- simulateResiduals(msro1); plot(t0) # ~OK
testDispersion(t0) # OK
# CI
msci <- glmer(cbind(GermDead,GermSurv)~MatTreat*SeedTreat+(1|PlantID),family=binomial(link="logit"),data=CI)
summary(msci)
msci1 <- glmer(cbind(GermDead,GermSurv)~MatTreat+SeedTreat+(1|PlantID),family=binomial(link="logit"),data=CI)
summary(msci1)
anova(msci,msci1) # df=1 Chisq=0.192 P=0.66
anova(msci1,update(msci1,~.-MatTreat)) # df=1 Chisq=0.0115 P=0.91
anova(msci1,update(msci1,~.-SeedTreat)) # df=1 Chisq=0.540 P=0.46
t0 <- simulateResiduals(msci1); plot(t0) # OK
testDispersion(t0) # OK


### 4/ Time to germination
# RO
mtro <- lmer(TimeGerm~MatTreat*SeedTreat+(1|PlantID),data=RO)
summary(mtro)
mtro1 <- lmer(TimeGerm~MatTreat+SeedTreat+(1|PlantID),data=RO)
summary(mtro1)
anova(mtro,mtro1) # df=1 Chisq=0.0904 P=0.76
anova(mtro1,update(mtro1,~.-SeedTreat)) # df=1 Chisq=2.87 P=0.090 .
mtro2 <- lmer(TimeGerm~MatTreat+(1|PlantID),data=RO)
anova(mtro1,mtro2)
anova(mtro2,update(mtro2,~.-MatTreat)) # df=1 Chisq=4.24 P=0.039 *
emmeans(mtro1,specs=pairwise ~ MatTreat | SeedTreat,type="response") # C<D P=0.038 *
tapply(RO$TimeGerm,RO$MatTreat,function(x){mean(na.omit(x))}) # C=14.5 D=20.9
tapply(RO$TimeGerm,RO$MatTreat,function(x){se(na.omit(x))}) # C=0.95 D=2.46
tapply(RO$TimeGerm,RO$MatTreat,function(x){length(which(na.omit(x)>0))}) # C=14 D=15
t0 <- simulateResiduals(mtro2); plot(t0) # ~OK
testDispersion(t0) # OK

mtci <- lmer(TimeGerm~MatTreat*SeedTreat+(1|PlantID),data=CI)
summary(mtci)
mtci1 <- lmer(TimeGerm~MatTreat+SeedTreat+(1|PlantID),data=CI)
summary(mtci1)
anova(mtci,mtci1) # df=1 Chisq=0.305 P=0.58
anova(mtci1,update(mtci1,~.-MatTreat)) # df=1 Chisq=0.672 P=0.41
anova(mtci1,update(mtci1,~.-SeedTreat)) # df=1 Chisq=1.82 P=0.18
t0 <- simulateResiduals(mtci1); plot(t0) # ~OK
testDispersion(t0) # OK


### 5/ Time from germination to cotyledon stage
# RO
mtsro <- lmer(TimeCoty~MatTreat*SeedTreat+(1|PlantID),data=RO)
summary(mtsro)
mtsro1 <- lmer(TimeCoty~MatTreat+SeedTreat+(1|PlantID),data=RO)
summary(mtsro1)
anova(mtsro,mtsro1) # df=1 Chisq=0.151 P=0.70
anova(mtsro1,update(mtsro1,~.-MatTreat)) # df=1 Chisq=2.19 P=0.14
anova(mtsro1,update(mtsro1,~.-SeedTreat)) # df=1 Chisq=0.312 P=0.58
tapply(RO$TimeCoty,RO$MatTreat,function(x){length(na.omit(x))}) # C=12 D=15
t0 <- simulateResiduals(mtsro1); plot(t0) # ~OK
testDispersion(t0) # OK

mtsci <- lmer(TimeCoty~MatTreat*SeedTreat+(1|PlantID),data=CI)
summary(mtsci)
mtsci1 <- lmer(TimeCoty~MatTreat+SeedTreat+(1|PlantID),data=CI)
summary(mtsci1)
anova(mtsci,mtsci1) # df=1 Chisq=1.89 P=0.17
anova(mtsci1,update(mtsci1,~.-MatTreat)) # df=1 Chisq=0.0885 P=0.77
anova(mtsci1,update(mtsci1,~.-SeedTreat)) # df=1 Chisq=3.19 P=0.074
mtsci2 <- lmer(TimeCoty~MatTreat+(1|PlantID),data=CI)
anova(mtsci1,mtsci2) # df=1 Chisq=3.19 P=0.074 .
tapply(CI$TimeCoty,CI$MatTreat,function(x){length(na.omit(x))}) # C=32 D=27
t0 <- simulateResiduals(mtsci1); plot(t0) # ~OK
testDispersion(t0) # OK


# Figure 2 Seed experiment ------------------------------------------------------------------
names(RO)[16] <- "Treat3"
names(CI)[16] <- "Treat3"
RO$Treat3 <- as.factor(as.character(RO$Treat3))
levels(RO$Treat3) <- c("CC","CH","DC","DH")
CI$Treat3 <- as.factor(as.character(CI$Treat3))
levels(CI$Treat3) <- c("CC","CH","DC","DH")
RO$Rviable <- RO$Viable/RO$Nseeds
RO$RgermCor <- RO$Germ/RO$Viable
RO$RSurv <- RO$GermSurv/RO$Germ
CI$Rviable <- CI$Viable/CI$Nseeds
CI$RgermCor <- CI$Germ/CI$Viable
CI$RSurv <- CI$GermSurv/CI$Germ
r0 <- rnorm(75,0,0.15); t0 <- c(1:4)[RO$Treat3]+r0[1:53]
t1 <- c(1:4)[CI$Treat3]+r0

par(mfrow=c(2,5),las=1,bty="l",mar=c(3,5,3,1),xpd=TRUE)
# 1/ RO
plot(RO$Rviable~RO$Treat3,col="#B69BDE",ylab="Viability rate",ylim=c(0,1),outline=F)
points(RO$Rviable~t0,pch=19,cex=0.8)
text(0.5,1.2,label=expression(italic("S. rosmarinus")),cex=1.1,adj=0)
text(0.5,1.1,label=expression(bold("A")),cex=1,adj=0)
plot(RO$RgermCor~RO$Treat3,col="#B69BDE",ylab="Corrected germination rate",ylim=c(0,1))
points(RO$RgermCor~t0,pch=19,cex=0.8)
text(0.5,1.1,label=expression(bold("B")),cex=1,adj=0)
plot(RO$RSurv~RO$Treat3,col="#B69BDE",ylab="Survival rate\\nto cotyledone stage",ylim=c(0,1),outline=F)
points(RO$RSurv~t0,pch=19,cex=0.8)
text(0.5,1.1,label=expression(bold("C")),cex=1,adj=0)
plot(RO$TimeGerm~RO$Treat3,col="#B69BDE",ylab="Mean time to germination (days)",ylim=c(0,52))
points(RO$TimeGerm~t0,pch=19,cex=0.8)
text(0.5,56,label=expression(bold("D")),cex=1,adj=0)
arrows(c(1,3),c(30,47),c(2,4),c(30,47),length=0);arrows(1:4,c(30,30,47,47),1:4,c(30,30,47,47)-2,length=0)
arrows(1.5,50,3.5,50,length=0,lwd=2);arrows(c(1.5,3.5),50,c(1.5,3.5),c(50-20,50-3),length=0)
text(2.5,53,labels="P = 0.039",cex=0.9)
plot(RO$TimeCoty~RO$Treat3,col="#B69BDE",ylab="Mean time from germination\\nto cotyledon stage (days)",ylim=c(0,15))
points(RO$TimeCoty~t0,pch=19,cex=0.8)
text(0.5,16,label=expression(bold("E")),cex=1,adj=0)

# 2/ CI
plot(CI$Rviable~CI$Treat3,col="#DEB2CC",ylab="Viability rate",ylim=c(0,1),outline=F)
points(CI$Rviable~t1,pch=19,cex=0.8)
text(0.5,1.25,label=expression(italic("C. albidus")),cex=1.1,adj=0)
text(0.5,1.1,label=expression(bold("F")),cex=1,adj=0)
text(0.5,56,label=expression(bold("D")),cex=1,adj=0)
arrows(c(1,3),1.05,c(2,4),1.05,length=0)
arrows(1.5,1.1,3.5,1.1,length=0,lwd=2);arrows(c(1.5,3.5),1.1,c(1.5,3.5),1.05,length=0)
text(2.5,1.14,labels="P = 0.052",cex=0.9)
plot(CI$RgermCor~CI$Treat3,col="#DEB2CC",ylab="Corrected germination rate",ylim=c(0,1))
points(CI$RgermCor~t1,pch=19,cex=0.8)
arrows(c(1,3),1.05,c(2,4),1.05,length=0)
arrows(1.5,1.1,3.5,1.1,length=0,lwd=2);arrows(c(1.5,3.5),1.1,c(1.5,3.5),1.05,length=0)
text(2.5,1.14,labels="P = 0.054",cex=0.9)
plot(CI$RSurv~CI$Treat3,col="#DEB2CC",ylab="Survival rate\\nto cotyledone stage",ylim=c(0,1))
points(CI$RSurv~t1,pch=19,cex=0.8)
text(0.5,1.1,label=expression(bold("H")),cex=1,adj=0)
plot(CI$TimeGerm~CI$Treat3,col="#DEB2CC",ylab="Mean time to germination\\n(days)",ylim=c(0,52),outline=F)
points(CI$TimeGerm~t1,pch=19,cex=0.8)
text(0.5,56,label=expression(bold("I")),cex=1,adj=0)
plot(CI$TimeCoty~CI$Treat3,col="#DEB2CC",ylab="Mean time from germination\\nto cotyledone (days)",ylim=c(0,15),outline=F)
points(CI$TimeCoty~t1,pch=19,cex=0.8)
text(0.5,16,label=expression(bold("J")),cex=1,adj=0)


