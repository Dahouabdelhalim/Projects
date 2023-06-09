library("car")
library("multcomp")
library("lme4")
library("ggplot2")
library("DHARMa")


## GEN1
GEN1<-read.csv('[pathname]/GEN1ants.csv',h=T)

m1i<-lmer(RMSD^2~indVar*P_M+(1|boxID),data=GEN1) ##  model w/ ind. attribute "indVar" (e.g., A vs. B) and social env. "P_M" (pure vs. mixed)
simulateResiduals(m1i, n = 1000, refit = F, integerResponse = NULL, plot = T, seed = 123)
drop1(m1i,test="Chisq")

m1<-lmer(RMSD^2~cat_treat+(1|boxID),data=GEN1) ## equivalent model w/ one 4-level treatment variable
drop1(m1,test="Chisq")

cM1 <- contrMat(n=table(GEN1$cat_treat), type='Tukey') # the full contrast matrix
lines <- c(1,3:4,6) # comparisons of interest
cM2 <- cM1[lines,]   
summary(glht(m1, linfct=mcp(cat_treat=cM2)),test = adjusted("holm")) 

## GEN2
GEN2<-read.csv('[pathname]/GEN2ants.csv',h=T)

m2b<-lmer(RMSD^(0.6)~indVar*P_M+(1|boxID),data=GEN2) 
simulateResiduals(m2b, n = 1000, refit = F, integerResponse = NULL, plot = T, seed = 123)
drop1(m2b,test="Chisq")

m2<-lmer(RMSD^(0.6)~cat_treat+(1|boxID),data=GEN2) 
drop1(m2,test="Chisq")

cM1 <- contrMat(n=table(GEN2$cat_treat), type="Tukey") 
lines <- c(1,3:4,6)
cM2 <- cM1[lines,]   
summary(glht(m2, linfct=mcp(cat_treat=cM2)),test = adjusted("holm")) 

## AGE
AGE<-read.csv('[pathname]/AGEants.csv',h=T)

m3b<-lmer(RMSD^2~indVar*P_M+(1|boxID),data=AGE) 
simulateResiduals(m3b, n = 1000, refit = F, integerResponse = NULL, plot = T, seed = 123)
drop1(m3b,test="Chisq")

m3<-lmer(RMSD^2~cat_treat+(1|boxID),data=AGE) 
drop1(m3,test="Chisq")

cM1 <- contrMat(n=table(AGE$cat_treat), type='Tukey') 
lines <- c(1,2,5:6)
cM2 <- cM1[lines,]  
summary(glht(m3, linfct=mcp(cat_treat=cM2)),test = adjusted("holm")) 

## IW
IW<-read.csv('[pathname]/IWants.csv',h=T)

m4b<-lmer(RMSD~indVar*P_M+(1|boxID),data=IW)
simulateResiduals(m4b, n = 1000, refit = F, integerResponse = NULL, plot = T, seed = 123)
drop1(m4b,test="Chisq")

m4<-lmer(RMSD~cat_treat+(1|boxID),data=IW)
drop1(m4,test="Chisq")

cM1 <- contrMat(n=table(IW$cat_treat), type='Tukey') 
lines <- c(1,3:4,6)
cM2 <- cM1[lines,]   
summary(glht(m4, linfct=mcp(cat_treat=cM2)),test = adjusted("holm")) 

