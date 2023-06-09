#This script conducts male mate recognition (MMR) analyses
require(arm)
require(glmmADMB)

mmr<-read.csv("Dataset3.csv")

#######################################################################
##### translocation experiment 1: analyses of titia male MR sympatry (CV) vs. allopatry (BF) ##### 
#######################################################################


mmr_titia<-subset(mmr,male.sp=="titia")
m1<-bayesglm(clasp~factor(female.treatment)*sympatric,data=mmr_titia,family="binomial")

#######################################################################
##### translocation experiment 2: analyses of occisa male MR sympatry (PA2, OT) vs. allopatry (ES, UOT) ##### 
#######################################################################


mmr_occisa<-subset(mmr,male.sp=="occisa")

mmr_occisa_translocation<-subset(mmr_occisa,translocation==1)

m2<-bayesglm(clasp~factor(female.treatment)*sympatric, data=mmr_occisa_translocation, family = "binomial")

#######################################################################
##### wing darkening experiment 1: analyses of americana male MR sympatry (CV) vs. allopatry (LM) ##### 
#######################################################################


mmr_americana<-subset(mmr,male.sp=="americana")

m3<-glmmadmb(clasp~factor(female.treatment,levels=c("dw","am"))*sympatric+(1|id),data=mmr_americana,family="binomial")

#######################################################################
##### wing darkening experiment 2: analyses of occisa male MR sympatry (PA2) vs. allopatry (ES, CT) ##### 
#######################################################################


mmr_occisa_wing.darkening<-subset(mmr_occisa,wing.darkening==1)

m4<-glmmadmb(clasp~female.treatment*factor(sympatric)+(1|id),data=mmr_occisa_wing.darkening,family="binomial")
