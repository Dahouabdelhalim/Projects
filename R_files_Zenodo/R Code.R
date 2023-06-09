#Libraries to upload
library(gtools)
library(ggplot2)
library(gplots)
library(lme4)
library(lattice)
library(lmerTest)
library(lsmeans)
library(scales)
library(blmeco)

#Upload Data

#Sperm precedence data: 
#SB1C (data file --> SB1C.csv)
 
#Female remating propensity data (all females): 
#SB1C.p (data file --> Remating.csv)

#Key differences between females remating data: 

#Data for analyses if IMI 0  - for all females in experiment: 
#SBfemrem0 (data file --> Remate0.csv)

#Data for analysis of IMI 24 - for all females in experiment that didn't remate at IMI 0:
#SBfemrem24 (data file --> Remate24.csv

#Data for analysis of IMI 48 - for all females in experiment that didn't remate at IMI 0 or 24: 
#SBfemrem48 (data file --> Remate48.csv)

#Supplementary data
#See Correction Factors for Sterile Male Technique.xlsx


table(SB1C$MatingOrder)
table(SB1C$SterileOrder)
table(SB1C$IMInterval)

se = function(x) sd(x)/sqrt(length(x))
mean(SB1C$P2[SB1C$IMInterval =="0"])
se(SB1C$P2[SB1C$IMInterval =="0"])
mean(SB1C$P2[SB1C$IMInterval =="24"])
se(SB1C$P2[SB1C$IMInterval =="24"])
mean(SB1C$P2[SB1C$IMInterval =="48"])
se(SB1C$P2[SB1C$IMInterval =="48"])

mean(SB1C$P2[SB1C$MatingOrder=="PM"])
se(SB1C$P2[SB1C$MatingOrder=="PM"])
mean(SB1C$P2[SB1C$MatingOrder=="MP"])
se(SB1C$P2[SB1C$MatingOrder=="MP"])

mean(SB1C$P2[SB1C$SterileOrder =="RN"])
se(SB1C$P2[SB1C$SterileOrder =="RN"])
mean(SB1C$P2[SB1C$SterileOrder =="NR"])
se(SB1C$P2[SB1C$SterileOrder =="NR"])

IMI0<-SB1C$P2[SB1C$IMInterval=="0"]
IMI24<-SB1C$P2[SB1C$IMInterval=="24"]
IMI48<-SB1C$P2[SB1C$IMInterval=="48"]
summary(IMI0)
summary(IMI24)
summary(IMI48)

#Age data for focal females and males
mean(SB1C.p$M1FEMAge,na.rm=TRUE)
#length=111
l=111
sd=sd(SB1C.p$M1FEMAge,na.rm=TRUE)
se=(sd/sqrt(l))
se

mean(SB1C.p$M1MALEAge,na.rm=TRUE)
l=111
sd=sd(SB1C.p$M1MALEAge,na.rm=TRUE)
se=(sd/sqrt(l))
se

mean(SB1C.p$M2MaleAge0,na.rm=TRUE)
l=111
sd=sd(SB1C.p$M2MaleAge0,na.rm=TRUE)
se=(sd/sqrt(l))
se

#Egg data for females

#Number of eggs laid before 2nd mating:
mean(SB1C$PrimarySeedEggs[SB1C$IMInterval =="24"],na.rm=TRUE)
l=18
sd=sd(SB1C$PrimarySeedEggs[SB1C$IMInterval =="24"],na.rm=TRUE)
se=(sd/sqrt(l))
se

mean(SB1C$PrimarySeedEggs[SB1C$IMInterval =="48"],na.rm=TRUE)
l=14
sd=sd(SB1C$PrimarySeedEggs[SB1C$IMInterval =="48"],na.rm=TRUE)
se=(sd/sqrt(l))
se 

#Number of eggs laid after 2nd mating: 
mean(SB1C$Fecundity[SB1C$IMInterval =="0"])
l=25
sd=sd(SB1C$Fecundity[SB1C$IMInterval =="0"],na.rm=TRUE)
se=(sd/sqrt(l))
se

mean(SB1C$Fecundity[SB1C$IMInterval =="24"])
l=24
sd=sd(SB1C$Fecundity[SB1C$IMInterval =="24"],na.rm=TRUE)
se=(sd/sqrt(l))
se

mean(SB1C$Fecundity[SB1C$IMInterval =="48"])
l=17
sd=sd(SB1C$Fecundity[SB1C$IMInterval =="48"],na.rm=TRUE)
se=(sd/sqrt(l))
se




############################ SPERM PRECEDENCE ###############################

#Figure 1

#violin plot
IMI0<-SB1C$P2[SB1C$IMInterval =="0"]
IMI24<-SB1C$P2[SB1C$IMInterval =="24"]
IMI48<-SB1C$P2[SB1C$IMInterval =="48"]
vioplot(IMI0,IMI24,IMI48,names=c("0","24","48"),col="gainsboro")

#GLMM

SB1C$OLRE = as.factor(1:dim(SB1C)[1])
#Determine random effects for model, unnested and nested
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemGen)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | FemGen: FemPatriline)+(1 | M1MaleGen)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline)+(1 | OLRE),data= SB1C,family=binomial)
summary(P2.glmer) #AIC 492.6
#,control=glmerControl(optimizer="bobyqa")
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemGen)+(1 | OLRE),data= SB1C,family=binomial)
summary(P2.glmer) #AIC 490.6,488.6,486.6,484.6,482.6,480.6,478.6,476.6,474.6, 472.6, 472, 
#removing OLRE

P2.noOLRE.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemGen),data= SB1C,family=binomial) 
summary(P2.noOLRE.glmer)  #AIC 855.1
#Tests if you have overdispersion
overdisp_fun <- function(model) {
    ## number of variance parameters in an n-by-n variance-covariance matrix
    vpars <- function(m) {
        nrow(m) * (nrow(m) + 1)/2
    }
    # The next two lines calculate the residual degrees of freedom
    model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
    rdf <- nrow(model.frame(model)) - model.df
    # extracts the Pearson residuals
    rp <- residuals(model, type = "pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    # Generates a p-value. If less than 0.05, the data are overdispersed.
    pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
    c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}
overdisp_fun(P2.noOLRE.glmer)  #YES I DO! So keep OLRE in!


#Now run bivariate analyses

SB1C$M1MatingLatency2<- scale(SB1C$M1MatingLatency, center = TRUE, scale = TRUE)
SB1C$M1KickingLatency2<- scale(SB1C$M1KickingLatency, center = TRUE, scale = TRUE)
SB1C$M1RemovalLatency2<- scale(SB1C$M1RemovalLatency, center = TRUE, scale = TRUE)
SB1C$M1CopulationDuration2 <- scale(SB1C$M1CopulationDuration,center = TRUE, scale = TRUE)
SB1C$M2MatingLatency2<- scale(SB1C$M2MatingLatency, center = TRUE, scale = TRUE)
SB1C$M2KickingLatency2<- scale(SB1C$M2KickingLatency, center = TRUE, scale = TRUE)
SB1C$M2RemovalLatency2<- scale(SB1C$M2RemovalLatency, center = TRUE, scale = TRUE)
SB1C$M2CopulationDuration2 <- scale(SB1C$M2CopulationDuration, center = TRUE, scale = TRUE)
SB1C$M2EjacDiff2 <- scale(SB1C$M2EjacDiff, center = TRUE, scale = TRUE)
SB1C$M2CopDurDiff2 <- scale(SB1C$M2CopDurDiff, center = TRUE, scale = TRUE)
SB1C$M2MatingLatencyDiff2 <- scale(SB1C$M2MatingLatencyDiff, center = TRUE, scale = TRUE)
SB1C$M2KickingLatencyDiff2 <- scale(SB1C$M2KickingLatencyDiff, center = TRUE, scale = TRUE)
SB1C$M2KickingDurationDiff2 <- scale(SB1C$M2KickingDurationDiff, center = TRUE, scale = TRUE)

SB1C$OLRE = as.factor(1:dim(SB1C)[1])
#,control=glmerControl(optimizer="bobyqa")
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~ MatingOrder +(1 | FemGen)+(1 | OLRE),data= SB1C,family=binomial)
summary(P2.glmer)

#M1FEMAge           0.1531        #collinear w/IMI
#M1FEMWeight1    			0.764
#M2MaleAgeDiff		0.0528		   #collinear w/IMI
#M1MALEAge					0.7461
#M1MALEWeight1				0.818  #collinear w/m1age
#M2MaleWeightDiff			0.33555  #coll w/m2maleagdiff
#M2EjacDiff 				0.23961
#M2MatingLatencyDiff2 		0.552620
#M2KickingLatencyDiff2		0.629831  #coll w/copdur
#M2KickingDurationDiff2		0.593630  #coll w/copdur
#M2CopDurDiff2				0.429677
#as.factor(IMInterval)  0.0479*, 0.0203*
#SterileOrder		 0.0906
#MatingOrder     	0.1654

#check collinearities
#MatingOrder+as.factor(IMInterval)+SterileOrder +M1FEMAge+ M2MaleAgeDiff+ M1MALEAge+ M2MaleWeightDiff+ M2EjacDiff+ M2CopDurDiff2+

lm1=lm(M1MALEAge ~ M2CopDurDiff2,data=SB1C)
summary(lm1)
#remember male ages (F1,64=71.5, P=5.09e-12***) and weights (F1,64=13.66,P=0.000455***) are highly collinear so not considered here, though their differences are due to variation in when the second males actually mated

#collinear
#M2MaleAgeDiff ~ M1MALEAge
#M2MaleAgeDiff ~ as.factor(IMInterval)
#M1FEMAge ~ as.factor(IMInterval)
#M2MaleWeightDiff ~ M2EjacDiff

#MatingOrder+as.factor(IMInterval)+SterileOrder+ M1MALEAge+ M2EjacDiff+ M2CopDurDiff2+

#final model selection
#,control=glmerControl(optimizer="bobyqa")
P2.glmer1=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~ MatingOrder+as.factor(IMInterval)+SterileOrder+ M1MALEAge+ M2EjacDiff+ M2CopDurDiff2+(1 | FemGen)+(1 | OLRE),data= SB1C,family=binomial)
summary(P2.glmer1) #AIC 400.2
#remove M1MALEAge
P2.glmer2=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~ MatingOrder+as.factor(IMInterval)+SterileOrder+ M2EjacDiff+ M2CopDurDiff2+(1 | FemGen)+(1 | OLRE),data= SB1C,family=binomial)
summary(P2.glmer2) #AIC 398.2
#remove SterileOrder,398.4
#remove M2EjacDiff, 402.3
#remove M2CopDurDiff2, 465.4

#FINAL MODEL IS P2.glmer2

SB1C$OLRE = as.factor(1:dim(SB1C)[1])
SB1C$M2CopDurDiff2 <- scale(SB1C$M2CopDurDiff, center = TRUE, scale = TRUE)
P2.glmer2=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~ MatingOrder+as.factor(IMInterval)+SterileOrder+ M2EjacDiff+ M2CopDurDiff2+(1 | FemGen)+(1 | OLRE),data= SB1C,family=binomial)
summary(P2.glmer2)


lsmeans(P2.glmer2, pairwise ~ as.factor(IMInterval),type="response")

#Fitted values
x<-0.5541284   #0
xse<-0.05788450
y<-0.7398948  #24 
yse<-0.04483346
z<-0.7415564   #48 
zse<-0.05350834


#Table 2

#(Intercept)              -0.2902     0.2687  -1.080  0.28026
exp(-0.2902)/(1+exp(-0.2902))
exp(-0.2902-1.96* 0.2687)/(1+exp(-0.2902-1.96* 0.2687))
exp(-0.2902+1.96* 0.2687)/(1+exp(-0.2902+1.96* 0.2687))

#as.factor(IMInterval)24 0.8281     0.3072   2.695  0.00703**
exp(0.8281)/(1+exp(0.8281))
exp(0.8281-1.96* 0.3072)/(1+exp(0.8281-1.96* 0.3072))
exp(0.8281 +1.96* 0.3072)/(1+exp(0.8281 +1.96* 0.3072))

#as.factor(IMInterval)48  0.8367     0.3588   2.332  0.01970*
exp(0.8367)/(1+exp(0.8367))
exp(0.8367-1.96* 0.3588)/(1+exp(0.8367-1.96* 0.3588))
exp(0.8367 +1.96* 0.3588)/(1+exp(0.8367 +1.96* 0.3588))

#MatingOrderPM            0.5338     0.2530   2.110  0.03489*
exp(0.5338)/(1+exp(0.5338))
exp(0.5338-1.96* 0.2530)/(1+exp(0.5338-1.96* 0.2530))
exp(0.5338 +1.96* 0.2530)/(1+exp(0.5338 +1.96* 0.2530))

#SterileOrderRN            0.3838     0.2580   1.488  0.13688
exp(0.3838)/(1+exp(0.3838))
exp(0.3838-1.96* 0.2580)/(1+exp(0.3838-1.96* 0.2580))
exp(0.3838 +1.96* 0.2580)/(1+exp(0.3838 +1.96* 0.2580))

#M2EjacDiff               -0.5715     0.6491  -0.880  0.37866
exp(-0.5715)/(1+exp(-0.5715))
exp(-0.5715-1.96* 0.6491)/(1+exp(-0.5715-1.96* 0.6491))
exp(-0.5715+1.96* 0.6491)/(1+exp(-0.5715+1.96* 0.6491))

#M2CopDurDiff2            -0.0427     0.1549  -0.276  0.78287
exp(-0.0427)/(1+exp(-0.0427))
exp(-0.0427-1.96* 0.1549)/(1+exp(-0.0427-1.96* 0.1549))
exp(-0.0427+1.96* 0.1549)/(1+exp(-0.0427+1.96* 0.1549))





############################ FEMALE REMATING PROPENSITY #########################


#Figure 2

table(SB1C.p$Remated,SB1C.p$IMInterval)

#IMI 0 = 26 out of 113
#IMI 24 = 25 out of 87
#IMI 48 = 18 out of 62
#NONE = 44 remaining

#Proportion of females remating
PropRem=c(0.230,0.221,0.159)
barplot(PropRem,beside=TRUE,col=c("gainsboro"),ylim=c(0,0.3),cex.lab=1.25,cex.main=1.25,cex.axis=1.25)
title(xlab="Inter-mating Interval (Hours)",ylab="Proportion of Females Remating",cex.lab=1.25,cex.main=1.25,cex.axis=1.25)
#multiple proportions test
Remated<-c(26, 25, 18)
Females<-c(113, 87, 62)
prop.test(Remated,Females)

mean(SB1C$P2[SB1C$IMInterval =="0"])
mean(SB1C$P2[SB1C$IMInterval =="24"])
mean(SB1C$P2[SB1C$IMInterval =="48"])

se = function(x) sd(x)/sqrt(length(x))

se(SB1C$P2[SB1C$IMInterval =="0"])
se(SB1C$P2[SB1C$IMInterval =="24"])
se(SB1C$P2[SB1C$IMInterval =="48"])

#multiple proportions test on raw data (summed not averaged)
Hatched<-c(992, 538, 268)
Unhatched<-c(1570, 969, 487)
prop.test(Hatched,Unhatched)





############### KEY DIFFERENCES BETWEEN FEMALES THAT DID AND DIDN'T REMATE ###############

#NOTES: 
#Remated - YES = 1
#Remated - NO = 0
#M1RemovalLatency = Mating 1 Kicking Duration



################### IMI 0 #####################

#Determine random factors
SBfemrem0$OLRE = as.factor(1:dim(SBfemrem0)[1])
#,control=glmerControl(optimizer="bobyqa")

SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | FemGen: FemPatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline)+(1 | Group)+(1 | OLRE),data= SBfemrem0, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem0.glmer) #AIC 149.6
#remove +(1 | OLRE)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | FemGen: FemPatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline)+(1 | Group),data= SBfemrem0, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem0.glmer) #AIC 147.6
#test for overdispersion
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}   
overdisp_fun(SBrem0.glmer)  #No dispersion! p = 0.9837536
#remove  +(1 | FemGen: FemPatriline)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline)+(1 | Group),data= SBfemrem0, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem0.glmer) #AIC 145.6
#remove +(1 | M2MaleGen: M2MalePatriline)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | Group),data= SBfemrem0, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem0.glmer) #AIC 143.6
#remove +(1 | M2MaleGen: M2MaleMatriline2)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | Group),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 141.6
#remove +(1 | M2MaleMatriline)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleGen: M2MaleMatriline)+(1 | Group),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 139.6
#remove +(1 | M2MaleGen: M2MaleMatriline)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | Group),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 137.6 
#remove +(1 | FemGen: FemMatriline2)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | Group),data= SBfemrem0, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem0.glmer) #AIC 135.6
#remove +(1 | FemMatriline)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen)+(1 | FemGen: FemMatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | Group),data= SBfemrem0, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem0.glmer) #AIC 133.6 
#remove +(1 | M1MaleGen: M1MaleMatriline)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen)+(1 | FemGen: FemMatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | Group),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 131.6
#remove +(1 | FemGen)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen: FemMatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | Group),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 129.6
#remove +(1 | M1MaleGen)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen: FemMatriline)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | Group),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 127.6
#remove +(1 | M1MaleMatriline)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen: FemMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | Group),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 125.6
#remove +(1 | Group)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen: FemMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 123.6 
#remove +(1 | M1MaleGen: M1MaleMatriline2)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen: FemMatriline)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 121.6
#remove +(1 | M1MaleGen: M1MalePatriline)
SBrem0.glmer=glmer(Remated0 ~ +(1 | FemGen: FemMatriline)+(1 | M2MaleGen),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 120.1
#remove (1 | M2MaleGen), 120.8
#remove +(1 | FemGen: FemMatriline),
SBrem0.glmer=glmer(Remated0 ~(1 | M2MaleGen) ,data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer) #AIC 118.7
overdisp_fun(SBrem0.glmer)  #No dispersion! p = 0.7064434

#Bivariate analyses
SBfemrem0$M1MatingLatency2 <- scale(SBfemrem0 $M1MatingLatency, center = TRUE, scale = TRUE)
SBfemrem0$M1KickingLatency2 <- scale(SBfemrem0 $M1KickingLatency, center = TRUE, scale = TRUE)
SBfemrem0$M1RemovalLatency2 <- scale(SBfemrem0 $M1RemovalLatency, center = TRUE, scale = TRUE)
SBfemrem0$M1CopulationDuration2 <-scale(SBfemrem0$M1CopulationDuration, center = TRUE, scale = TRUE)

SBrem0.glmer=glmer(Remated0 ~ MatingOrder + (1 | M2MaleGen) ,data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer)

#M1FEMAge					0.2864
#M1FEMWeight1	0.1490
#M1FEMWeightGain			0.2710
#M1MALEAge					0.988
#M1MALEWeight1				0.370
#M1MALERawEjac		0.0344 *
#M1MatingLatency2			0.2064
#M1KickingLatency2	0.0449 *
#M1RemovalLatency2	0.011223 * 
#M1CopulationDuration2	0.13401
#M2Male2Weight0		0.148
#M2AgeDiff0					0.481
#M2MaleWeightDiff 			0.2326
#SterileOrder				0.6387
#MatingOrder				0.8920

#M1FEMWeight1+ M1MALERawEjac+M1KickingLatency2+ M1RemovalLatency2+ M1CopulationDuration2+ M2Male2Weight0

#Check collinearities
lm1=lm(M1RemovalLatency2 ~ M1CopulationDuration2,data= SBfemrem0)
summary(lm1)
plot(M1MALERawEjac ~ M1CopulationDuration,data= SBfemrem0)
abline(coef(lm1),col=2)

#Collinear
#M1MALERawEjac (keep) ~ M1KickingLatency2
#M1MALERawEjac ~ M1CopulationDuration  (slight correlation, P = 0.0504)
#M1MALERawEjac (keep) ~ M2Male2Weight0
#M1RemovalLatency2(keep) ~ M1CopulationDuration2

#Final predictors for consideration
#MatingOrder + SterileOrder + M1FEMWeight1+ M1MALERawEjac+ M1RemovalLatency2

#Final model selection
SBrem0.glmer1=glmer(Remated0 ~ MatingOrder + SterileOrder + M1FEMWeight1+ M1MALERawEjac+ M1RemovalLatency2 + (1 | M2MaleGen),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer1) #AIC 107.7
#remove SterileOrder
SBrem0.glmer2=glmer(Remated0 ~ MatingOrder + M1FEMWeight1+ M1MALERawEjac+ M1RemovalLatency2 + (1 | M2MaleGen),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer2) #AIC 106.1
#remove MatingOrder
SBrem0.glmer3=glmer(Remated0 ~ M1FEMWeight1+ M1MALERawEjac+ M1RemovalLatency2 + (1 | M2MaleGen),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer3) #AIC 104.2
#remove M1FEMWeight1
SBrem0.glmer4=glmer(Remated0 ~ M1MALERawEjac+ M1RemovalLatency2 + (1 | M2MaleGen),data= SBfemrem0, family=binomial) 
summary(SBrem0.glmer4) #AIC 103.6
#remove M1MALERawEjac, AIC 105

#final model is SBrem0.glmer4



#TABLE 3: 

#(Intercept)        -0.5188     0.7707  -0.673   0.5009
exp(est)/(1+exp(est))
exp(est-1.96*se)/(1+exp(est-1.96*se))
exp(est+1.96*se)/(1+exp(est+1.96*se))
  
#M1MALERawEjac      -3.2024     1.8380  -1.742   0.0815
exp(-3.2024)/(1+exp(-3.2024))
exp(-3.2024-1.96*1.8380)/(1+exp(-3.2024-1.96*1.8380))
exp(-3.2024+1.96*1.8380)/(1+exp(-3.2024+1.96*1.8380))

#M1RemovalLatency2  -1.3443     0.5481  -2.453   0.0142 *
exp(-1.3443 )/(1+exp(-1.3443 ))
exp(-1.3443 -1.96*0.5481)/(1+exp(-1.3443 -1.96*0.5481))
exp(-1.3443 +1.96*0.5481)/(1+exp(-1.3443 +1.96*0.5481))

#mean M1RemovalLatency (i.e. kicking duration) for yes and no
tapply(SBfemrem0$M1RemovalLatency,SBfemrem0$Remated0,mean,na.rm=TRUE)
#verified in Excel data sheet

#calculation of length using R is incorrect, comes out to 87, but is in fact 84
#calculating SE manually for remated = NO 
sdn<-sd(SBfemrem0$M1RemovalLatency[SBfemrem0$Remated0=="0"],na.rm=TRUE)
ln<-84
seno<-sdn/sqrt(ln)
seno  #

#calculation of length using R is incorrect, comes out to 26, but is in fact 24
#calculating SE manually for remated = YES 
sdy<-sd(SBfemrem0$M1RemovalLatency[SBfemrem0$Remated0=="1"],na.rm=TRUE)
ly<-24
seyes<-sdy/sqrt(ly)
seyes  #

#Checking up on male ejaculate sizes
boxplot2(M1MALERawEjac ~Remated0,xlab="Female Remated Immediately",ylab="First Male Ejaculate Size (mg)",col=c("gainsboro"),data= SBfemrem0,ylim=c(-0.1,0.8),cex.lab=1.25,cex.main=1.25,cex.axis=1.25)
text(1,-0.1,"(87)")
text(2,-0.1,"(26)")

tapply(SBfemrem0$M1MALERawEjac,SBfemrem0$Remated0,mean,na.rm=TRUE)
#calculating SE manually for remated = NO 
sdn<-sd(SBfemrem0$M1MALERawEjac[SBfemrem0$Remated0=="0"],na.rm=TRUE)
ln<-length(SBfemrem0$M1MALERawEjac[SBfemrem0$Remated0=="0"])
seno<-sdn/sqrt(ln)
seno  
#calculating SE manually for remated = YES 
sdy<-sd(SBfemrem0$M1MALERawEjac[SBfemrem0$Remated0=="1"],na.rm=TRUE)
ly<-length(SBfemrem0$M1MALERawEjac[SBfemrem0$Remated0=="1"])
seyes<-sdy/sqrt(ly)
seyes  






################### IMI 24 #####################

#Determine random factors
SBfemrem24$OLRE = as.factor(1:dim(SBfemrem24)[1])
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | FemGen: FemPatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline)+(1 | Group)+(1 | OLRE),data= SBfemrem24, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem24.glmer) #AIC 139.9
#remove +(1 | OLRE)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | FemGen: FemPatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline)+(1 | Group),data= SBfemrem24, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem24.glmer) #AIC  137.9

##test for overdispersion
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}   
overdisp_fun(SBrem24.glmer)  #NO overdispersion, P=0.3261874

#remove +(1 | Group)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | FemGen: FemPatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline),data= SBfemrem24, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem24.glmer) #AIC 135.9
#remove +(1 | FemGen: FemPatriline)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline),data= SBfemrem24, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem24.glmer) #AIC 133.9
#remove +(1 | M1MaleGen: M1MalePatriline)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline),data= SBfemrem24, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem24.glmer) #AIC 131.9
#remove +(1 | M2MaleGen: M2MalePatriline)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2),data= SBfemrem24, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem24.glmer) #AIC 129.9
#remove +(1 | M2MaleGen: M2MaleMatriline2)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline),data= SBfemrem24, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem24.glmer) #AIC 127.9
#remove +(1 | M1MaleGen: M1MaleMatriline2)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline),data= SBfemrem24, family=binomial) 
summary(SBrem24.glmer) #AIC 125.9
#remove +(1 | FemGen: FemMatriline2)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline),data= SBfemrem24, family=binomial) 
summary(SBrem24.glmer) #AIC 123.9
#remove +(1 | M2MaleMatriline)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | M2MaleGen: M2MaleMatriline),data= SBfemrem24, family=binomial) 
summary(SBrem24.glmer) #AIC 121.9
#remove +(1 | FemGen: FemMatriline)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | M2MaleGen: M2MaleMatriline),data= SBfemrem24, family=binomial) 
summary(SBrem24.glmer) #AIC 119.9
#remove +(1 | FemMatriline)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | M2MaleGen: M2MaleMatriline),data= SBfemrem24, family=binomial) 
summary(SBrem24.glmer) #AIC 117.9
#remove +(1 | M2MaleGen)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline),data= SBfemrem24, family=binomial) 
summary(SBrem24.glmer) #AIC 115.9
#remove +(1 | M1MaleGen: M1MaleMatriline)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline),data= SBfemrem24, family=binomial) 
summary(SBrem24.glmer) #AIC 113.9
#remove +(1 | M2MaleGen: M2MaleMatriline)
SBrem24.glmer=glmer(Remated24 ~ +(1 | FemGen)+(1 | M1MaleGen)+(1 | M1MaleMatriline),data= SBfemrem24, family=binomial) 
summary(SBrem24.glmer) #AIC 111.9
#remove +(1 | FemGen)
SBrem24.glmer=glmer(Remated24 ~ +(1 | M1MaleGen)+(1 | M1MaleMatriline),data= SBfemrem24, family=binomial) 
summary(SBrem24.glmer) #AIC 109.9
#remove +(1 | M1MaleGen)
SBrem24.glmer=glmer(Remated24 ~ +(1 | M1MaleMatriline),data= SBfemrem24, family=binomial) 
summary(SBrem24.glmer) #AIC 107.9 

#bivariate analyses
SBfemrem24$M1MatingLatency2 <- scale(SBfemrem24 $M1MatingLatency, center = TRUE, scale = TRUE)
SBfemrem24 $M1KickingLatency2 <- scale(SBfemrem24 $M1KickingLatency, center = TRUE, scale = TRUE)
SBfemrem24 $M1RemovalLatency2 <- scale(SBfemrem24 $M1RemovalLatency, center = TRUE, scale = TRUE)
SBfemrem24 $M1CopulationDuration2 <-scale(SBfemrem24 $M1CopulationDuration, center = TRUE, scale = TRUE)

SBrem24.glmer=glmer(Remated24 ~ M2MaleAge24 +(1 | M1MaleMatriline),data= SBfemrem24,family=binomial) 
summary(SBrem24.glmer) 

#M1FEMAge,  0.439
#M1FEMWeight1,  0.774
#M1FEMWeightGain,  0.701

#M1MALEAge, 0.616 (need to use ,control=glmerControl(optimizer="bobyqa"))
#M1MALEWeight1, 0.278
#M1MALERawEjac,	0.5924
#M1MatingLatency2, 0.32443
#M1KickingLatency2, 0.52116
#M1RemovalLatency2, 		0.120284
#M1CopulationDuration2,		0.16289

#M2FemaleAge24, 0.313
#M2FemWeight24, 			0.119

#M2MaleAge24, 0.934
#M2MaleWeight24, 0.663

#M2AgeDiff24, 0.405
#M2MaleWeightDiff (time of mating only),		0.18256
#M2FemWeightLoss0to24,		0.04711 *
#PrimarySeedEggs,			0.0124 *
#PrimarySeedEggsPROP,		 0.00490 **
#SterileOrder, 0.73255
#MatingOrder, 0.3144

#MatingOrder + SterileOrder+ M1RemovalLatency2 + M1CopulationDuration2+ M2MaleWeightDiff+ M2FemWeight24 + M2FemWeightLoss0to24+ PrimarySeedEggsPROP

#Check interactions for collinearities
lm1=lm(M2FemWeight24 ~ PrimarySeedEggs,data= SBfemrem24)
summary(lm1)

#collinear
#M1RemovalLatency2 (keep) ~ M1CopulationDuration2
#M1RemovalLatency2 ~ M2FemWeightLoss0to24 (keep)
#M1CopulationDuration2 ~ M2FemWeightLoss0to24 (keep)
#M2FemWeight24 ~ M2FemWeightLoss0to24 (keep)
#PrimarySeedEggsPROP (keep) ~ M2FemWeightLoss0to24 #so add M1RemovalLatency2 & M2FemWeight24  back


#Predictors to include 
#MatingOrder + SterileOrder+ M1RemovalLatency2+ M2FemWeight24 + M2MaleWeightDiff+ PrimarySeedEggsPROP

#,control=glmerControl(optimizer="bobyqa")
SBrem24.glmer1=glmer(Remated24 ~ MatingOrder + SterileOrder+ M1RemovalLatency2+ M2FemWeight24 + M2MaleWeightDiff+ PrimarySeedEggsPROP +(1 | M1MaleMatriline),data= SBfemrem24,family=binomial) 
summary(SBrem24.glmer1) #AIC 33.1 

#turn into a glm because m1malematriline has no effect on response variable

#Bivariate analyses of predictors using GLM: 
SBrem24.glm=glm(Remated24 ~ MatingOrder,data= SBfemrem24,family=binomial) 
summary(SBrem24.glm) 

#M1FEMAge,  0.303
#M1FEMWeight1,  0.634
#M1FEMWeightGain,  0.716

#M1MALEAge, 0.646 (need to use ,control=glmerControl(optimizer="bobyqa"))
#M1MALEWeight1, 0.304
#M1MALERawEjac,	0.5476
#M1MatingLatency2, 0.302
#M1KickingLatency2, 0.663816
#M1RemovalLatency2, 		0.114307
#M1CopulationDuration2,		0.153863    

#M2FemaleAge24, 0.331
#M2FemWeight24, 			0.116

#M2MaleAge24, 0.946
#M2MaleWeight24, 0.626

#M2AgeDiff24, 0.40791
#M2MaleWeightDiff (time of mating only),		0.1899
#M2FemWeightLoss0to24,							0.04417 *
#PrimarySeedEggs,			0.0124 *
#PrimarySeedEggsPROP,		 0.00490 **
#SterileOrder, 0.76043
#MatingOrder, 0.2664

#MatingOrder + SterileOrder+ M1RemovalLatency2 + M1CopulationDuration2+ M2MaleWeightDiff+ M2FemWeight24 + M2FemWeightLoss0to24+ PrimarySeedEggs +PrimarySeedEggsPROP

#Check interactions for collinearities
lm1=lm(M1RemovalLatency2 ~ M2FemWeightLoss0to24,data= SBfemrem24)
summary(lm1)

#collinear
#PrimarySeedEggsPROP (keep) ~ PrimarySeedEggs
#M1RemovalLatency2 (keep) ~ M1CopulationDuration2
#M1RemovalLatency2 ~ M2FemWeightLoss0to24 (keep)
#M1CopulationDuration2 ~ M2FemWeightLoss0to24 (keep)
#M2FemWeight24 ~ M2FemWeightLoss0to24 (keep)
#PrimarySeedEggsPROP (keep) ~ M2FemWeightLoss0to24 #so add M1RemovalLatency2 & M2FemWeight24  back

#so final factors:
#MatingOrder + SterileOrder+ M1RemovalLatency2 + M2FemWeight24 +M2MaleWeightDiff+ PrimarySeedEggsPROP

SBrem24.glm1=glm(Remated24 ~ MatingOrder + SterileOrder+ M1RemovalLatency2+ M2FemWeight24+ M2MaleWeightDiff+ PrimarySeedEggsPROP,data= SBfemrem24,family=binomial) 
summary(SBrem24.glm1) #AIC 31.066
#remove M1RemovalLatency2
SBrem24.glm2=glm(Remated24 ~ MatingOrder + SterileOrder+ M2FemWeight24 + M2MaleWeightDiff+ PrimarySeedEggsPROP,data= SBfemrem24,family=binomial) 
summary(SBrem24.glm2) #AIC 29.494 
#remove MatingOrder
SBrem24.glm3=glm(Remated24 ~ SterileOrder+ M2FemWeight24 +M2MaleWeightDiff+ PrimarySeedEggsPROP,data= SBfemrem24,family=binomial) 
summary(SBrem24.glm3) #AIC 27.628
#remove SterileOrder
SBrem24.glm4=glm(Remated24 ~ M2MaleWeightDiff+ M2FemWeight24 + PrimarySeedEggsPROP,data= SBfemrem24,family=binomial) 
summary(SBrem24.glm4) #AIC 27.385
#remove M2MaleWeightDiff
SBrem24.glm5=glm(Remated24 ~ M2FemWeight24 +PrimarySeedEggsPROP,data= SBfemrem24,family=binomial) 
summary(SBrem24.glm5) #AIC 25.583

#Final model is SBrem24.glm5 


#TABLE 3: 

#(Intercept)           38.219     15.723   2.431  0.01507 * 

#PrimarySeedEggsPROP  -19.261      7.344  -2.623  0.00873 **
exp(-19.261)/(1+exp(-19.261))
exp(-19.261-1.96* 7.344)/(1+exp(-19.261-1.96* 7.344))
exp(-19.261+1.96* 7.344)/(1+exp(-19.261+1.96* 7.344))


#M2FemWeight24         -5.199      2.240  -2.321  0.02031 *
exp(-5.199)/(1+exp(-5.199))
exp(-5.199-1.96* 2.240)/(1+exp(-5.199-1.96* 2.240))
exp(-5.199 +1.96* 2.240)/(1+exp(-5.199+1.96* 2.240))
 

#calculate mean for Primary Seed Eggs Proportion for remated and unremated fems
tapply(SBfemrem24 $PrimarySeedEggsPROP, SBfemrem24 $Remated24,mean,na.rm=TRUE)

#R incorrectly calculates length of 62, but should be 14
#calculating SE manually for remated = NO 
sdn<-sd(SBfemrem24$PrimarySeedEggsPROP[SBfemrem24$Remated24=="0"],na.rm=TRUE)
ln<-14
seno<-sdn/sqrt(ln)
seno  

#R incorrectly calculates length of 25, but should be 18
#calculating SE manually for remated = YES 
sdy<-sd(SBfemrem24$PrimarySeedEggsPROP[SBfemrem24$Remated24=="1"],na.rm=TRUE)
ly<-18
seyes<-sdy/sqrt(ly)
seyes

#calculate mean weight of females prior to second mating for remated and unremated fems
tapply(SBfemrem24$M2FemWeight24, SBfemrem24 $Remated24,mean,na.rm=TRUE)

#R incorrectly calculates length of 62, but should be 49
#calculating SE manually for remated = NO 
sdn<-sd(SBfemrem24$M2FemWeight24[SBfemrem24$Remated24=="0"],na.rm=TRUE)
ln<-49
seno<-sdn/sqrt(ln)
seno  

#R incorrectly calculates length of 25, but should be 19
#calculating SE manually for remated = YES 
sdy<-sd(SBfemrem24$M2FemWeight24[SBfemrem24$Remated24=="1"],na.rm=TRUE)
ly<-19
seye<-sdy/sqrt(ly)
seye



################### IMI 48 #####################

#Determine random factors
SBfemrem48$OLRE = as.factor(1:dim(SBfemrem48)[1])
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | FemGen: FemPatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline)+(1 | Group)+(1 | OLRE),data= SBfemrem48, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem48.glmer) #AIC 110.7
#remove +(1 | OLRE)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | FemGen: FemPatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline)+(1 | Group),data= SBfemrem48, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem48.glmer) #AIC 108.7
#test for overdispersion
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}   
overdisp_fun(SBrem48.glmer)  #YES dispersion! p = 0.04704377

#so keep OLRE!
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | FemGen: FemPatriline)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline)+(1 | Group)+(1 | OLRE),data= SBfemrem48, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem48.glmer) #AIC 110.7
#remove +(1 | FemGen: FemPatriline)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | M2MaleGen: M2MalePatriline)+(1 | Group)+(1 | OLRE),data= SBfemrem48, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem48.glmer) #AIC 108.7
#remove +(1 | M2MaleGen: M2MalePatriline)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M1MaleGen: M1MalePatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | Group)+(1 | OLRE),data= SBfemrem48, family=binomial,control=glmerControl(optimizer="bobyqa")) 
summary(SBrem48.glmer) #AIC 106.7
#remove +(1 | M1MaleGen: M1MalePatriline)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline2)+(1 | Group)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 104.7
#remove +(1 | M2MaleGen: M2MaleMatriline2)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline2)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | Group)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 102.7
#remove +(1 | M1MaleGen: M1MaleMatriline2)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | M2MaleMatriline)+(1 | M2MaleGen: M2MaleMatriline)+(1 | Group)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 100.7 
#remove +(1 | M2MaleMatriline)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleGen)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | M2MaleGen: M2MaleMatriline)+(1 | Group)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 98.7
#remove +(1 | M1MaleGen)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | FemGen: FemMatriline2)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | M2MaleGen: M2MaleMatriline)+(1 | Group)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 96.7
#remove +(1 | FemGen: FemMatriline2)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | M2MaleGen: M2MaleMatriline)+(1 | Group)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 94.7
#remove +(1 | Group)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | M2MaleGen: M2MaleMatriline)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 92.7
#remove +(1 | M2MaleGen: M2MaleMatriline)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | M1MaleMatriline)+(1 | M1MaleGen: M1MaleMatriline)+(1 | M2MaleGen)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 90.7
#remove +(1 | M1MaleGen: M1MaleMatriline)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemMatriline)+(1 | FemGen: FemMatriline)+(1 | M1MaleMatriline)+(1 | M2MaleGen)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 88.7
#remove +(1 | FemMatriline)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemGen: FemMatriline)+(1 | M1MaleMatriline)+(1 | M2MaleGen)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 
#remove +(1 | M2MaleGen)  86.7
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | FemGen: FemMatriline)+(1 | M1MaleMatriline)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 84.7
#remove +(1 | FemGen: FemMatriline)
SBrem48.glmer=glmer(Remated48 ~ +(1 | FemGen)+(1 | M1MaleMatriline)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 82.7
#remove  +(1 | FemGen)
SBrem48.glmer=glmer(Remated48 ~ +(1 | M1MaleMatriline)+(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 80.7
#remove +(1 | M1MaleMatriline)
SBrem48.glmer=glmer(Remated48 ~ +(1 | OLRE),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer) #AIC 78.7 
#No effect of OLRE on glmer, but overdispersion detected earlier, so try replacing with another random effect and testing for overdispersion again.
SBrem48.glmer1=glmer(Remated48 ~ +(1 | M1MaleMatriline),data= SBfemrem48, family=binomial) 
summary(SBrem48.glmer1)
overdisp_fun(SBrem48.glmer1) #no overdisperion detected, P = 0.4046522
so need to collapse to a glm

#bivariate analyses
SBfemrem48$M1MatingLatency2 <- scale(SBfemrem48 $M1MatingLatency, center = TRUE, scale = TRUE)
SBfemrem48 $M1KickingLatency2 <- scale(SBfemrem48 $M1KickingLatency, center = TRUE, scale = TRUE)
SBfemrem48 $M1RemovalLatency2 <- scale(SBfemrem48 $M1RemovalLatency, center = TRUE, scale = TRUE)
SBfemrem48 $M1CopulationDuration2 <-scale(SBfemrem48 $M1CopulationDuration, center = TRUE, scale = TRUE)

SBrem48.glm=glm(Remated48 ~ MatingOrder,data=SBfemrem48, family="binomial")
summary(SBrem48.glm) 

#M1FEMAge, 							0.117
#M1FEMWeight1, 						0.192
#M1FEMWeightGain, 0.680
#M2FemaleAge48, 					0.052
#M2FemWeight48, 0.861
#M2FemWeightLoss0to48,				0.01947 * 
#PrimarySeedEggs, don't have enough data for these! Only were counted for females that did remate, not for females that did not remate at this time point

#M1MALEAge, 0.252
#M1MALEWeight1, 0.519
#M1MALERawEjac, 0.2356
#M1MatingLatency2, 0.89121
#M1KickingLatency2, 0.64122
#M1RemovalLatency2, 0.209739
#M1CopulationDuration2, 0.33760

#M2MaleAge48, 						0.0304 *
#M2MaleWeight48,					0.0582

#M2AgeDiff48, 0.210
#M2MaleWeightDiff (time of mating only), 0.24957

#SterileOrder, 						0.07049
#MatingOrder, 0.8141


#Final predictors to consider:
#MatingOrder + SterileOrder +M1FEMAge+ M1FEMWeight1+ M2FemaleAge48+ M2FemWeightLoss0to48 + M2MaleAge48+ M2Male2Weight48+

#check collinearities
lm1=lm(M2FemaleAge48 ~ M2FemWeightLoss0to48,data= SBfemrem48)
summary(lm1)
#plot(M2FemaleAge48 ~ M2FemWeightLoss0to48,data= SBfemrem48)
#abline(coef(lm1),col=2)

#collinear
#M1FEMAge ~ M2FemaleAge48 (keep)
#M1FEMWeight1 ~ M2MaleAge48 (keep)
#M2MaleAge48 (keep) ~ M2Male2Weight48
#M2MaleAge48 ~ M2FemWeightLoss0to48 (keep)
#M2Male2Weight48 ~ M2FemWeightLoss0to48 (keep)
#M2FemaleAge48 ~ M2FemWeightLoss0to48 (borderline collinear, keep femweightloss)

#Analysis with M1FEMWeight1
#MatingOrder + SterileOrder + M2FemaleAge48 + M2FemWeightLoss0to48 + M1FEMWeight1
SBrem48.glm1=glm(Remated48 ~ MatingOrder + SterileOrder + M2FemaleAge48 + M2FemWeightLoss0to48 + M1FEMWeight1,data=SBfemrem48, family="binomial")
summary(SBrem48.glm1) #AIC 45.161
#remove MatingOrder
SBrem48.glm2=glm(Remated48 ~ SterileOrder + M2FemaleAge48 + M2FemWeightLoss0to48 + M1FEMWeight1,data=SBfemrem48, family="binomial")
summary(SBrem48.glm2) #AIC 43.659
#Remove SterileOrder, AIC is 45.146, so worse

boxplot(M1FEMWeight1~Remated48,data=SBfemrem48)

#Final model is SBrem48.glm2!


#TABLE 3: 

#(Intercept)          -12.8530     6.4143  -2.004	0.0451 *
exp(est)/(1+exp(est))
exp(est-1.96*se)/(1+exp(est-1.96*se))
exp(est+1.96*se)/(1+exp(est+1.96*se))

#SterileOrderRN         1.6234     0.9229   1.759	0.0786
exp(1.6234)/(1+exp(1.6234))
exp(1.6234-1.96* 0.9229)/(1+exp(1.6234-1.96* 0.9229))
exp(1.6234 +1.96* 0.9229)/(1+exp(1.6234 +1.96* 0.9229))

#M2FemaleAge48         -1.2683     0.5095  -2.489	0.0128 *
exp(-1.2683)/(1+exp(-1.2683))
exp(-1.2683-1.96* 0.5095)/(1+exp(-1.2683-1.96* 0.5095))
exp(-1.2683+1.96* 0.5095)/(1+exp(-1.2683+1.96* 0.5095))

#M1FEMWeight1           1.7869     0.8498   2.103	0.0355 *
exp(1.7869)/(1+exp(1.7869))
exp(1.7869-1.96* 0.8498)/(1+exp(1.7869-1.96* 0.8498))
exp(1.7869 +1.96* 0.8498)/(1+exp(1.7869 +1.96* 0.8498))

#M2FemWeightLoss0to48   3.4913     1.4409   2.423	0.0154 *
exp(3.4913)/(1+exp(3.4913))
exp(3.4913-1.96* 1.4409)/(1+exp(3.4913-1.96* 1.4409))
exp(3.4913 +1.96* 1.4409)/(1+exp(3.4913 +1.96* 1.4409))


#means for female age at time of second mating for remated and unremated fems
tapply(SBfemrem48$M2FemaleAge48,SBfemrem48$Remated48,mean,na.rm=TRUE)
#R correctly calculates the length as 44
#Calculating SE manually for remated = NO 
sdn<-sd(SBfemrem48$M2FemaleAge48[SBfemrem48$Remated48=="0"],na.rm=TRUE)
ln<-length(SBfemrem48$M2FemaleAge48[SBfemrem48$Remated48=="0"])
seno<-sdn/sqrt(ln)
seno  
#R incorrectly calculates the length as 18, when it is 16
#calculating SE manually for remated = YES 
sdy<-sd(SBfemrem48$M2FemaleAge48[SBfemrem48$Remated48=="1"],na.rm=TRUE)
ly<-16
seye<-sdy/sqrt(ly)
seye

#means for female initial weight at time of first mating for remated and unremated females
tapply(SBfemrem48$M1FEMWeight1,SBfemrem48$Remated48,mean,na.rm=TRUE)
#calculating SE manually for remated = NO 
sd<-sd(SBfemrem48$M1FEMWeight1[SBfemrem48$Remated48=="0"],na.rm=TRUE)
ly<-44  #Length for "0" / unremated = 44
se1<-sd/sqrt(ly)
se1  #0.08816876
#calculating SE manually for remated = YES 
sd<-sd(SBfemrem48$M1FEMWeight1[SBfemrem48$Remated48=="1"],na.rm=TRUE)
ly<-18  #Length for "1" / remated = 18
se1<-sd/sqrt(ly)
se1  #0.2192489


#means for female weight loss at time of second mating for remated and unremated fems
tapply(SBfemrem48$M2FemWeightLoss0to48,SBfemrem48$Remated48,mean,na.rm=TRUE)
#calculating SE manually for remated = NO 
#R incorrectly calculates the length as 44, when it is 34
sd<-sd(SBfemrem48$M2FemWeightLoss0to48[SBfemrem48$Remated48=="0"],na.rm=TRUE)
ln<-34
se1<-sd/sqrt(ln)
se1  
#calculating SE manually for remated = YES 
#R incorrectly calculates the length as 18, when it is 14
sd<-sd(SBfemrem48$M2FemWeightLoss0to48[SBfemrem48$Remated48=="1"],na.rm=TRUE)
ly<-14
se1<-sd/sqrt(ly)
se1  #0.0765799



########################          POST-HOC ANALYSES:       ########################

#Do females that lose more weight between matings lay more eggs? #YES
#(Need to use data from IMI 24 for comparison with missing egg data for IMI 48)
lm1=lm(M2FemWeightLoss0to24 ~ PrimarySeedEggs,data= SBfemrem24)
summary(lm1)
plot(M2FemWeightLoss0to24 ~ PrimarySeedEggs,data= SBfemrem24)
abline(coef(lm1),col=2)

#What is the relationship between female size (i.e. weight) and fecundity? #Larger females more fecund, borderline significant
lm1=lm(M1FEMWeight1 ~ FecundityTOTAL,data= SB1C)
summary(lm1)
plot(M1FEMWeight1 ~ FecundityTOTAL,data= SB1C)
abline(coef(lm1),col=2)

#What is the relationship between female age and fecundity?   #Younger females more fecund
lm1=lm(M1FEMAge ~ FecundityTOTAL,data= SB1C)
summary(lm1)
plot(M1FEMAge ~ FecundityTOTAL,data= SB1C)
abline(coef(lm1),col=2)


#What is the relationship between female size (i.e. weight) and male ejaculate size?  #None
lm1=lm(M1FEMWeight1 ~ M1MALERawEjac,data= SB1C)
summary(lm1)
plot(M1FEMWeight1 ~ M1MALERawEjac,data= SB1C)
abline(coef(lm1),col=2)

#What is the relationship between female size (i.e. weight) and the number of eggs laid between matings?   #None
lm1=lm(M1FEMWeight1 ~ PrimarySeedEggs,data= SB1C)
summary(lm1)
plot(M1FEMWeight1 ~ PrimarySeedEggs,data= SB1C)
abline(coef(lm1),col=2)
