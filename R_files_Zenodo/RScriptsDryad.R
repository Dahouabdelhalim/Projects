#################################################All Species#########################################################

LF2Data=read.csv('ALLKinResultsMaxLF.csv')
attach(LF2Data)
lm3 = lm(MaxProt~PreyVel)
plot(MaxProt~PreyVel)
abline(lm3)
summary(lm3)

MZ2Data=read.csv('ALLKinResultsMaxMZ.csv')
attach(MZ2Data)
lm4 = lm(MaxProt~PreyVel)
plot(MaxProt~PreyVel)
abline(lm4)
summary(lm4)

HybData=read.csv('ALLKinResultsMaxHyb.csv')
attach(HybData)
lm5 = lm(MaxProt~PreyVel)
plot(MaxProt~PreyVel)
abline(lm5)
summary(lm5)

AllKinMaxData=read.csv('ALLKinResultsMax.csv')
plot(MaxProt~Species,data=AllKinMaxData)


zebData=read.csv('zeb_and_original_results.csv')
plot(MaxProtTime~Species,data=zebData)
plot(MaxProt~Species,data=zebData)
modzeb=lm(MaxProtTime~Species,data=zebData)
summary(modzeb)

#The rest

allData=read.csv('AllSpeciesMax.csv')
allDataExclude=read.csv('AllSpeciesMaxnoblanks.csv')
attach(allData)

library(ppcor)

pcor.test(allDataExclude$PC1All,allDataExclude$Max.Vel.Resid,allDataExclude$Prot.Resid)
pcor.test(allDataExclude$Prot.Resid,allDataExclude$Max.Vel.Resid,allDataExclude$PC1All)

AllLF=subset(allData,allData$Species=='LF')
mean(AllLF$Standard.Length, na.rm=TRUE)
sd(AllLF$Standard.Length, na.rm=TRUE)
AllMZ=subset(allData,allData$Species=='MZ')
mean(AllMZ$Standard.Length, na.rm=TRUE)
sd(AllMZ$Standard.Length, na.rm=TRUE)
AllHyb=subset(allData,allData$Species=='Hybrid')
mean(AllHyb$Standard.Length, na.rm=TRUE)
sd(AllHyb$Standard.Length, na.rm=TRUE)


AllLF=subset(allData,allData$Species=='LF')
plot(Max.Vel.Resid~PC1All,data=AllLF,col='blue',cex=2.2,pch=16)
linmod=lm(Max.Vel.Resid~PC1All,data=AllLF)
linmodoutl=lm(Max.Vel.Resid~PC1AllNoOutl,data=AllLF)
abline(linmod,lwd=3,col='blue')
abline(linmodoutl,lwd=3,col='lightcyan4')
summary(linmod)
summary(linmodoutl)

AllMZ=subset(allData,allData$Species=='MZ')
plot(Max.Vel.Resid~PC1All,data=AllMZ,col='red',cex=2.2,pch=16)
linmod=lm(Max.Vel.Resid~PC1All,data=AllMZ)
abline(linmod,lwd=3,col='red')
summary(linmod)

Species = factor(Species, levels = c("LF", "Hybrid", "MZ"))
boxplot(Max.Vel.Resid~Species, ylab='Prey Velocity Residuals',ylim=c(-.23,.55),labels=FALSE,lwd=2)
x = aov(Max.Vel.Resid~Species,)
summary(x)
TukeyHSD(x)

boxplot(Prot.Resid~Species, ylab='Protrusion Residuals',ylim=c(-.0013,.0018),lwd=2)
z = aov(Prot.Resid~Species)
summary(z)
TukeyHSD(z)

boxplot(PC1All~Species, ylab='PC1',ylim=c(-.123,.11),lwd=2)
zz = aov(PC1All~Species)
summary(zz)
TukeyHSD(zz)

mod1=lm(Max.Vel.Resid~Species)
summary(mod1)

leg.txt=c("LF","Hybrid","MZ")

plot(Max.Vel.Resid[Species=='Hybrid']~PC1All[Species=='Hybrid'],ylim=c(-.25,.49),pch=16,cex=2.2,col='violetred',xlim = c(-.14,.09), ,xlab='',ylab='')
legend(.0633,.35,leg.txt,col=c('blue','violetred','red'),pch=16,cex=1.4)
points(Max.Vel.Resid[Species=='MZ']~PC1All[Species=='MZ'],pch=16,cex=2.2,col='red')
points(Max.Vel.Resid[Species=='LF']~PC1All[Species=='LF'],pch=16,cex=2.2,col='blue')
modHybn = lm(Max.Vel.Resid~PC1All)
polymodPC = lm(Max.Vel.Resid~poly(PC1All,2,raw=TRUE))
summary(modHybn)
summary(polymodPC)
abline(modHybn,lwd=3)

modHybn2 = lm(Max.Vel.Resid~PC2All)
summary(modHybn2)
modHybn3 = lm(Max.Vel.Resid~PC3All)
summary(modHybn3)
modHybn4 = lm(Max.Vel.Resid~PC4All)
summary(modHybn4)
modHybn5 = lm(Max.Vel.Resid~PC5All)
summary(modHybn5)


plot(Max.Vel.Resid~Prot.Resid)

plot(Max.Vel.Resid[Species=='Hybrid']~Prot.Resid[Species=='Hybrid'],pch=16,cex=2.2,col='violetred',xlim = c(-.0018,.0019), ylim=c(-.25,.49),xlab='Protrusion Residuals', ylab='Prey Velocity Residuals')
legend(.00147,.451,leg.txt,col=c('blue','violetred','red'),pch=16,cex=1.4)
points(Max.Vel.Resid[Species=='MZ']~Prot.Resid[Species=='MZ'],pch=16,cex=2.2,col='red')
points(Max.Vel.Resid[Species=='LF']~Prot.Resid[Species=='LF'],pch=16,cex=2.2,col='blue')
protMod1 = lm(Max.Vel.Resid~Prot.Resid)
abline(protMod1,lwd=3)
summary(protMod1)


attach(allDataExclude)

plot(Prot.Resid[Species=='Hybrid']~PC1All[Species=='Hybrid'],pch=16,cex=2.2,col='violetred',xlim = c(-.14,.09), ylim = c(-.0018,.0019),xlab='Protrusion Residuals', ylab='PC1 All')
legend(.0633,-.00095,leg.txt,col=c('blue','violetred','red'),pch=16,cex=1.4)
points(Prot.Resid[Species=='MZ']~PC1All[Species=='MZ'],pch=16,cex=2.2,col='red')
points(Prot.Resid[Species=='LF']~PC1All[Species=='LF'],pch=16,cex=2.2,col='blue')
protMod11 = lm(Prot.Resid~PC1All)
polymod = lm(Prot.Resid~poly(PC1All,2,raw=TRUE))
abline(protMod11,lwd=3)
summary(polymod)
summary(protMod11)


plot(Max.Vel.Resid[Species=='MZ']~Prot.Resid[Species=='MZ'],pch=16,cex=2.2)
moddddd=lm(Max.Vel.Resid[Species=='MZ']~Prot.Resid[Species=='MZ'])
summary(moddddd)
modddd=lm(Prot.Resid[Species=='MZ']~Max.Vel.Resid[Species=='MZ'])

summary(modddd)

plot(Max.Vel.Resid[Species=='LF']~Prot.Resid[Species=='LF'],pch=16,cex=2.2)
modddddd=lm(Max.Vel.Resid[Species=='LF']~Prot.Resid[Species=='LF'])
summary(modddddd)




############################################Size Correction###############################################

sizeCorr = read.csv('AllSpeciesMax.csv')
attach(sizeCorr)

# Prey velocity, all species.
corrModX=lm(Max.Vel~Standard.Length)
corrModX
summary(corrModX)
plot(Max.Vel[Species=='Hybrid']~Standard.Length[Species=='Hybrid'],col='green')
points(Max.Vel[Species=='MZ']~Standard.Length[Species=='MZ'],col='red')
points(Max.Vel[Species=='LF']~Standard.Length[Species=='LF'],col='blue')
abline(corrModX)
plot(corrModX$residuals)
write.csv(corrModX$residuals,file="AllSpecies_Max.Vel_residuals.csv")


# protrusion correction, all species.
corrMod2=lm(Prot~Standard.Length)
corrMod2
summary(corrMod2)
plot(Prot[Species=='Hybrid']~Standard.Length[Species=='Hybrid'],col='green')
points(Prot[Species=='MZ']~Standard.Length[Species=='MZ'],col='red')
points(Prot[Species=='LF']~Standard.Length[Species=='LF'],col='blue')
abline(corrMod2)
plot(corrMod2$residuals)
write.csv(corrMod2$residuals,file="AllSpecies_Prot_residuals.csv")

moddddd = lm(Prot~Standard.Length)
moddddd$residuals





#################################################Hybrids##################################################

allData=read.csv('AllSpecGen.csv')
attach(allData)

###PC1hyb###

plot(Max.Vel~PC1hyb,pch=16,cex=2,xlim=c(-.07,.043),ylim=c(.35,.835))
modbroad=lm(Max.Vel~PC1hyb)
abline(modbroad,lwd=3)
summary(modbroad)

#RUNX2#

unique(Runx2)
color=c('violetred','blue','red','magenta','blue')
color[unclass(Runx2)]

plot(Max.Vel[Runx2=='LF']~PC1hyb[Runx2=='LF'])

RunxHybHet=subset(allData,allData$Runx2=='Het')
modHHet = lm(Max.Vel~PC1hyb,data=RunxHybHet)
plot(Max.Vel~PC1hyb,col=color[unclass(Runx2)],pch=16,cex=2.2, data=RunxHybHet,main='Runx2 Hets')
summary(modHHet)
abline(modHHet)

RunxHybLF=subset(allData,allData$Runx2=='LF')
modHLF = lm(Max.Vel~PC1hyb,data=RunxHybLF)
plot(Max.Vel~PC1hyb,col=color[unclass(Runx2)],pch=16,cex=2.2, data=RunxHybLF,main='Runx2 LF')
summary(modHLF)
abline(modHLF)

RunxHybHet=subset(allData,allData$Runx2=='Het')
RunxHybLF=subset(allData,allData$Runx2=='LF')
modHHet = lm(Max.Vel~PC1hyb,data=RunxHybHet)
modHLF = lm(Max.Vel~PC1hyb,data=RunxHybLF)
plot(Max.Vel~PC1hyb,col=color[unclass(Runx2)],pch=16,cex=2.2, data=RunxHybHet,main='Runx2',xlim=c(-.07,.043),ylim=c(.35,.835))
points(Max.Vel~PC1hyb,col=color[unclass(Runx2)],pch=16,cex=2.2, data=RunxHybLF)
abline(modHHet,col='violetred',lwd=3)
abline(modHLF,col='blue',lwd=3)
summary(modHHet)
summary(modHLF)



#HSP47#

unique(HSP47)
color=c('violetred','red','red','blue')
color[unclass(HSP47)]

HSP47HybHet=subset(allData,allData$HSP47=='Het')
modHHet = lm(Max.Vel~PC1hyb,data=HSP47HybHet)
plot(Max.Vel~PC1hyb,col=color[unclass(HSP47)],pch=16,cex=2, data=HSP47HybHet,main='Hets')
summary(modHHet)
abline(modHHet)

HSP47HybLF=subset(allData,allData$HSP47=='MZ')
modHLF = lm(Max.Vel~PC1hyb,data=HSP47HybLF)
plot(Max.Vel~PC1hyb,col=color[unclass(HSP47)],pch=16,cex=2, data=HSP47HybLF,main='MZ')
summary(modHLF)
abline(modHLF)

HSPHybHet=subset(allData,allData$HSP47=='Het')
HSPHybMZ=subset(allData,allData$HSP47=='MZ')
modHHet = lm(Max.Vel~PC1hyb,data=HSPHybHet)
modHMZ = lm(Max.Vel~PC1hyb,data=HSPHybMZ)
plot(Max.Vel~PC1hyb,col=color[unclass(HSP47)],pch=16,cex=2.2, data=HSPHybHet,main='HSP47',xlim=c(-.07,.043),ylim=c(.35,.835))
points(Max.Vel~PC1hyb,col=color[unclass(Runx2)],pch=16,cex=2.2, data=HSPHybMZ)
abline(modHHet,col='violetred',lwd=3)
abline(modHMZ,col='blue',lwd=3)
summary(modHHet)
summary(modHMZ)


#Sox9b#

unique(Sox9b)
color=c('violetred','blue','red','magenta','blue')
color[unclass(Sox9b)]


Sox9bHybHet=subset(allData,allData$Sox9b=='Het')
modHHet = lm(Max.Vel~PC1hyb,data=Sox9bHybHet)
plot(Max.Vel~PC1hyb,col=color[unclass(Sox9b)],pch=16,cex=2, data=Sox9bHybHet,main='Sox9b Hets')
summary(modHHet)
abline(modHHet)

Sox9bHybLF=subset(allData,allData$Sox9b=='LF')
modHLF = lm(Max.Vel~PC1hyb,data=Sox9bHybLF)
plot(Max.Vel~PC1hyb,col=color[unclass(Sox9b)],pch=16,cex=2, data=Sox9bHybLF,main='Sox9b LF')
summary(modHLF)
abline(modHLF)


Sox9bHybHet=subset(allData,allData$Sox9b=='Het')
Sox9bHybLF=subset(allData,allData$Sox9b=='LF')
modHHet = lm(Max.Vel~PC1hyb,data=Sox9bHybHet)
modHLF = lm(Max.Vel~PC1hyb,data=Sox9bHybLF)
plot(Max.Vel~PC1hyb,col=color[unclass(Sox9b)],pch=16,cex=2, data=Sox9bHybHet,main='Sox9b',xlim=c(-.07,.043),ylim=c(.35,.835))
points(Max.Vel~PC1hyb,col=color[unclass(Sox9b)],pch=16,cex=2, data=Sox9bHybLF)
abline(modHHet,col='violetred',lwd=3)
abline(modHLF,col='blue',lwd=3)
summary(modHHet)
summary(modHLF)









#########################Permutation Analysis#########################

####First input all the data, then remove all NaN values, as well as non-hybrids
allData=read.csv('AllSpecGen.csv')
abbrevData = allData[c(4,17,29:30)]
abbrevDataHyb = subset(abbrevData,Runx2=='LF'|Runx2=='Het')
abbrevDataHybComp = abbrevDataHyb[complete.cases(abbrevDataHyb),]
attach(abbrevDataHybComp)

#######Create function to return p value from linear regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


#create vectors to hold R^2 and p values from models
allR2Het = NULL
allpHet = NULL
allR2LF = NULL
allpLF = NULL

for (i in 1:10000){
  
  ####Randomize Runx2 values
  Runx2Rand = sample(abbrevDataHybComp$Runx2)
  abbrevDataHybComp$Runx2Rand = Runx2Rand
  #attach(abbrevDataHybComp)
  
  ###Rerun graphs and models for randomized values
  RunxHybHet=subset(abbrevDataHybComp,abbrevDataHybComp$Runx2Rand=='Het')
  RunxHybLF=subset(abbrevDataHybComp,abbrevDataHybComp$Runx2Rand=='LF')
  modHHet = lm(Max.Vel.Resid~PC1hyb,data=RunxHybHet)
  modHLF = lm(Max.Vel.Resid~PC1hyb,data=RunxHybLF)
  #plot(Max.Vel.Resid~PC1hyb,col='violetred',pch=16,cex=2.2, data=RunxHybHet,main='Runx2',xlim=c(-.07,.043),ylim=c(.35,.77))
  #points(Max.Vel.Resid~PC1hyb,col='blue',pch=16,cex=2.2, data=RunxHybLF)
  #abline(modHHet,col='violetred',lwd=3)
  #abline(modHLF,col='blue',lwd=3)
  #abline(modbroad,lwd=3)
  summary(modHHet)
  summary(modHLF)
  
  #get r^2 and p and add to vector
  rHet = summary(modHHet)$r.squared
  pHet = lmp(modHHet)
  rLF = summary(modHLF)$r.squared
  pLF = lmp(modHLF)
  allR2Het = append(allR2Het,rHet)
  allR2LF = append(allR2LF,rLF)
  allpHet = append(allpHet,pHet)
  allpLF = append(allpLF,pLF)
}


histBreak = seq(0,1,by=.05)
histBreakpval = seq(0,1,by=.005)

hist(allR2Het,breaks=histBreakpval,main = "Het R^2 Values")
hist(allR2LF,breaks=histBreakpval,main = "LF Homozygote R^2 Values")

#### Calculate percentile of R^2 of Sox9b and Runx2 measured values
ecdf(allR2LF)(.3755)
ecdf(allR2LF)(.4595)
ecdf(allR2Het)(.01775)
ecdf(allR2Het)(.0002528)

hist(allpLF)