
rm(list=ls())

###Multihost code
###Phil Trans Paper
###Stuart K. J. R. Auld
data1<-read.csv(file.choose())
attach(data1)
names(data1)


###Updated 15/07/2016

spp<-factor(hostspp)
infected<-factor(infected)
line<-factor(clone)
ind<-factor(sheetno)

##load relevant libraries
library(lme4)
library(nlme)
library(survival)
library(coxme)


#Expt1 Analysis of proportion infected. Files: Expt1PastFull.csv; Expt1MetschFull.csv
m0<-glmer(infected~1+(1|spp/line/ind),subset=(treat=="1"),family=binomial)
m1<-glmer(infected~spp+(1|spp/line/ind),subset=(treat=="1"),family=binomial)
m2<-glmer(infected~spp+(1|spp/ind),subset=(treat=="1"),family=binomial)
summary(m2)
anova(m0,m1,test="LRT")
anova(m1,m2,test="LRT")

#Expt1 Analysis of spore burdens Files: Expt1PastFull.csv; Expt1MetschFull.csv
s0<-glmer(as.integer(spd)~1+(1|spp/line/ind),subset=(infected=="1"),family=poisson)
s1<-glmer(as.integer(spd)~spp+(1|spp/line/ind),subset=(infected=="1"),family=poisson)
s2<-glmer(as.integer(spd)~spp+(1|spp/ind),subset=(infected=="1"),family=poisson)
summary(s1)
anova(s0,s1,test="LRT")
anova(s1,s2,test="LRT")


#Expt1 Analysis of host fecundity. Files: Expt1PastFull.csv; Expt1MetschFull.csv
f0<-glmer(as.integer(totaloff)~1+(1|spp/line/ind),subset=(treat=="1"),family=poisson)
f1<-glmer(as.integer(totaloff)~infected*spp+(1|spp/line/ind),subset=(treat=="1"),family=poisson)
f2<-glmer(as.integer(totaloff)~infected+spp+(1|spp/line/ind),subset=(treat=="1"),family=poisson)
f3<-glmer(as.integer(totaloff)~infected+(1|spp/line/ind),subset=(treat=="1"),family=poisson)
f4<-glmer(as.integer(totaloff)~spp+(1|spp/line/ind),subset=(treat=="1"),family=poisson)
f5<-glmer(as.integer(totaloff)~infected*spp+(1|spp/ind),subset=(treat=="1"),family=poisson)
summary(f1)
anova(f3,f2,test="LRT")###testing for a spp effect
anova(f1,f2,test="LRT")###testing the inf x spp interaction
anova(f4,f2,test="LRT")###testing for an inf effect


#Expt1 Welch's t.test of transmission potential (Beta(sigma/tau)) by species. Files: Expt1PastBST.csv; Expt1MetschBST.csv
BSD<-Beta*sigdeath
t.test(BSD~spp)

#Expt1 Analysis of survival. Files: Expt1PastFull.csv; Expt1MetschFull.csv
status<-1*(death<60)
l0 <- coxme(Surv(death, status) ~ 1 + (1|spp/line/ind))
l1 <- coxme(Surv(death, status) ~ infected*spp + (1|spp/line/ind))
l2 <- coxme(Surv(death, status) ~ infected+spp + (1|spp/line/ind))
l3 <- coxme(Surv(death, status) ~ infected + (1|spp/line/ind))
l4 <- coxme(Surv(death, status) ~ spp + (1|spp/line/ind))
l5 <- coxme(Surv(death, status) ~ spp + (1|spp/ind))
anova(l3,l2,test="LRT")###testing for a spp effect
anova(l1,l2,test="LRT")###testing the inf x spp interaction
anova(l4,l2,test="LRT")###testing for an inf effect
anova(l4,l2,test="LRT")###testing for an inf effect

library(visreg)
visreg(l1,day,by="spp")

#Expt1 Test relationship between parasite growth rate and host reproductive rate (Fungus only). Files: Expt1MetschFull.csv

rparagrowth<-sqrt(spd/death)
rhostfecrate<-sqrt(totaloff/death)
m1<-lme(rparagrowth ~ rhostfecrate*spp,random=~1|spp/line,na=na.omit)

#Expt1 Test relationship between parasite burden and host day of death (Fungus only). Files: Expt1MetschFull.csv
rootspores<-sqrt(spd)
m1<-lme(rootspores ~ death*spp,random=~1|spp/line,na=na.omit)


#Expt2 Analysis of transmission rate (Beta) with first host species (host1) as a fixed effect and second host line (host 2) as a random effect. Files: Expt2PastBST.csv; Expt2MetschBST.csv
 
m1<-lme(beta~host1,random=~1|host2,method="ML")

#Expt2 Analysis of overall transmission potential (Beta/(sifma/tau)) with first host species (host1) as a fixed effect and second host line (host 2) as a random effect. Files: Expt2PastBST.csv; Expt2MetschBST.csv

m2<-lme(bs.t~host1,random=~1|host2,method="ML")

