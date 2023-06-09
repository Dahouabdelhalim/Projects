##### R Code for: Waller et al. (2018) 
### "Trait responses to AM fungi are stronger and more consistent than 
### fixed differences among populations of Asclepias speciosa." 
### Published in: American Journal of Botany.

m3 <- read.csv("Waller_AJB_Data.csv", header=T)  ## load full dataset
head(m3)
# load packages
library(lattice)
library(lme4)
library(lmerTest)
library(MuMIn)
library(reshape2)
library(lsmeans)

#######################################################################
######  RESHAPE FOR MULTIVARIATE TRAIT ANALYSIS #######################
### reshape centered trait values to long form ####################

melt3 <- melt(m3, measure.vars=c("TotBMs","Hgtdays","RSs","slas",
                                 "latexs","trics"), 
              id.vars=c("Pop","Fam","AMF","Seed_g"), 
              variable.name='Trait')

### examine histograms and boxplots
histogram(~value|Trait, data=melt3)
bwplot(value~AMF|Trait, data=melt3)


#########################################################
### full univariate models to answer Q1 

um1 <- lmer(value~AMF+scale(Seed_g)+(1|Pop)+(1|Pop:AMF), data=melt3, subset=Trait=="TotBMs")
um1a <- lmer(value~AMF+scale(Seed_g)+(1|Pop)+(1|Pop:Fam)+(1|Pop:AMF), 
             data=melt3, subset=Trait=="TotBMs") #### test for Family effect
anova(um1,um1a)
summary(um1)
summary(um1a)
rand(um1a)

anova(um1)
summary(um1)
rand(um1)
r.squaredGLMM(um1)

um2 <- lmer(value~AMF+scale(Seed_g)+(1|Pop)+(1|Pop:AMF), data=melt3, subset=Trait=="Hgtdays")
anova(um2)
summary(um2)
rand(um2)
r.squaredGLMM(um2)

um3 <- lmer(value~AMF+scale(Seed_g)+(1|Pop)+(1|Pop:AMF), data=melt3, subset=Trait=="RSs")
anova(um3)
summary(um3)
rand(um3)
r.squaredGLMM(um3)


um4 <- lmer(value~AMF+scale(Seed_g)+(1|Pop)+(1|Pop:AMF), data=melt3, subset=Trait=="slas")
anova(um4)
summary(um4)
rand(um4)
r.squaredGLMM(um4)


um5 <- lmer(value~AMF+scale(Seed_g)+(1|Pop)+(1|Pop:AMF), data=melt3, subset=Trait=="latexs")
anova(um5)
summary(um5)
rand(um5)
r.squaredGLMM(um5)


um6 <- lmer(value~AMF+scale(Seed_g)+(1|Pop)+(1|Pop:AMF), data=melt3, subset=Trait=="trics")
anova(um6)
summary(um6)
rand(um6)
r.squaredGLMM(um6)

####################################################################
### calculate AMF responses to answer Q2
### use raw trait values
####################################################################
wide1 <- dcast(m3, Pop+Fam~AMF,
               value.var="TotBM", mean)
names(wide1)[names(wide1)=="Yes"] <-"YT"
names(wide1)[names(wide1)=="No"] <-"NT"
wide1$amfresp <- log(wide1$YT/wide1$NT)
wide1$amfresp2 <- ((wide1$YT-wide1$NT)/wide1$YT)*100

wide2 <- dcast(m3, Pop+Fam~AMF,
               value.var="RS", mean)
names(wide2)[names(wide2)=="Yes"] <-"YR"
names(wide2)[names(wide2)=="No"] <-"NR"
wide2$RSresp <- log(wide2$YR/wide2$NR)
wide2$RSresp2 <- ((wide2$YR-wide2$NR)/wide2$YR)*100


wide3 <- dcast(m3, Pop+Fam~AMF,
               value.var="sla", mean)
names(wide3)[names(wide3)=="Yes"] <-"YS"
names(wide3)[names(wide3)=="No"] <-"NS"
wide3$slaresp <- log(wide3$YS/wide3$NS)
wide3$slaresp2 <- ((wide3$YS-wide3$NS)/wide3$YS)*100


wide4 <- dcast(m3, Pop+Fam~AMF,
               value.var="latex", mean)
names(wide4)[names(wide4)=="Yes"] <-"YL"
names(wide4)[names(wide4)=="No"] <-"NL"
wide4$latexresp <- log(wide4$YL/wide4$NL)
wide4$latexresp2 <- ((wide4$YL-wide4$NL)/wide4$YL)*100


wide5 <- dcast(m3, Pop+Fam~AMF,
               value.var="tric", mean)
names(wide5)[names(wide5)=="Yes"] <-"YTR"
names(wide5)[names(wide5)=="No"] <-"NTR"
wide5$tricresp <- log(wide5$YTR/wide5$NTR)
wide5$tricresp2 <- ((wide5$YTR-wide5$NTR)/wide5$YTR)*100


wide6 <- dcast(m3, Pop+Fam~AMF,
               value.var="Hgtday", mean)
names(wide6)[names(wide6)=="Yes"] <-"YH"
names(wide6)[names(wide6)=="No"] <-"NH"
wide6$Hgtresp <- log(wide6$YH/wide6$NH)
wide6$Hgtresp2 <- ((wide6$YH-wide6$NH)/wide6$YH)*100


w1 <- merge(wide1,wide2, by=c("Pop","Fam"))
w2 <- merge(w1,wide3, by=c("Pop","Fam"))
w3 <- merge(w2,wide4, by=c("Pop","Fam"))
w4 <- merge(w3,wide5, by=c("Pop","Fam"))
w5 <- merge(w4,wide6, by=c("Pop","Fam"))

w6 <- aggregate(cbind(amfresp2,RSresp2,slaresp2,latexresp2,tricresp2,Hgtresp2)~Pop, data=w5, mean)

library(Hmisc)
library(ape)
library(mvtnorm)
library(vegan)

###### compare two ways to calculate responsiveness
cor(w5$amfresp,w5$amfresp2)
cor(w5$slaresp,w5$slaresp2)

######  Calculate correlation matrix among responsiveness of traits
###### use percent responsiveness
respmat2 <- as.matrix(w5[,c("amfresp2","Hgtresp2","RSresp2","slaresp2","latexresp2","tricresp2")])
respcorr <- rcorr(respmat2)
respcorr


################################################################################
########## ANAYLSES FROM APPENDICES ############################################
################################################################################


##### correlations among traits in +AMF and -AMF treatments
fam <- aggregate(cbind(TotBMs,Hgtdays,RSs,slas,latexs,trics)~Pop+Fam+AMF, data=m3, mean)

m3Yes <- subset(fam, subset=AMF=="Yes")
m3Yes <- as.matrix(m3Yes[c("TotBMs","Hgtdays","RSs","slas","latexs","trics")])
rcorr(m3Yes)

m3No <- subset(fam, subset=AMF=="No")
m3No <- as.matrix(m3No[c("TotBMs","Hgtdays","RSs","slas","latexs","trics")])
rcorr(m3No)


#### multivariate model from Appendix 3
lm.3 <- lmer(value~AMF*Trait+(1|Pop)+(1|Pop:Trait:AMF),data=melt3,
             control=lmerControl(calc.derivs=F, optCtrl = list(maxfun=50000)))
anova(lm.3)




