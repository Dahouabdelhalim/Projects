library(ggplot2)

Data<-BehavioralAssayFinal

Data$Block<-as.factor(Data$Block)
Data$SP<-Data$TotalSexualP
Data$WP<-Data$TotalWP

Data$Foraging<-as.numeric(Data$Foraging)
Data$Aggression<-as.numeric(Data$Aggression)
Data$ExploratoryRate<-as.numeric(Data$ExploratoryRate)
Data$ForagerCoverage<-as.numeric(Data$GroupExploration)
Data$GroupExploration<-as.numeric(Data$ColonyExploration)

mass$SexRatio<-as.numeric((Data$TotalMale/(Data$TotalMale+mass$TotalGyne)))
mass$CasteRatio<-as.numeric((Data$TotalGyne/(Data$WPCasteRatio+mass$TotalGyne)))

###############################################################################
#Heritability

library(MCMCglmm)
library(nadiv)
library(lme4)
library(nlme)
pedigree<-read.csv("Pedigree.csv")
pedigree<-Pedigree

mass<-read.csv("Behavioral Assay Data.csv")
mass<-BehavioralAssayFinal

mass$Block<-as.factor(mass$Block)
mass$Wolbachia<-as.factor(mass$Predicted.Wolbachia)

#Foraging
hist(mass$Foraging)

mass$Foraging<-(mass$Foraging)

model<-lm(Foraging~Block+Queen+Wolbachia, data=mass)
summary(model)
drop1(model,test="Chi")#Block and Queen

mass1=merge(mass,pedigree,by.x="Genotype",by.y="animal",all.X=TRUE)
mass1$animal=interaction(mass1$Colony,mass1$Genotype)
pedigreeH5A=mass1[c("animal","dam","sire","sex")]
pedigree<-pedigree[c("animal","dam","sire","sex")]

#bind the portion of the pedigree from before H500 with the H500 portion where we have mass data
pedigree2=rbind(pedigreeH5A,pedigree)
#create the A matrix 
Mat_A <- makeS(prepPed(pedigree2),heterogametic="M",returnS=F)

prior <- list(R = list(V=1, nu=0.002), G = list(G1 = list(V=1, nu=0.002) )) 
MCMC <- MCMCglmm(Foraging~ Block+Wolbachia, random=~ animal,  ginverse=list(animal=Mat_A$Sinv),
                 nitt = 1e8, burnin = 1e4, thin = 5e2, family = c(rep("gaussian",1)), pl = T , pr = T, data = mass1,prior = prior,singular.ok=TRUE)

HPDinterval  (MCMC$VCV[, "animal"]/(MCMC$VCV[, "animal"] + MCMC$VCV[, "units"]))   
traceplot    (MCMC$VCV[, "animal"]/(MCMC$VCV[, "animal"] + MCMC$VCV[, "units"]))
densplot     (MCMC$VCV[, "animal"]/(MCMC$VCV[, "animal"] + MCMC$VCV[, "units"]))
effectiveSize(MCMC$VCV[, "animal"]/(MCMC$VCV[, "animal"] + MCMC$VCV[, "units"]))
autocorr.diag(MCMC$VCV[, "animal"]/(MCMC$VCV[, "animal"] + MCMC$VCV[, "units"]))

posterior.mode(MCMC$VCV[, "animal"]/(MCMC$VCV[, "animal"] + MCMC$VCV[, "units"]))
mean(MCMC$VCV[,1]/rowSums(MCMC$VCV))#Heritability
HPDinterval(MCMC$VCV[,1]/rowSums(MCMC$VCV))#95% CIs


###Repeat for each phenotype

#########################################################################################
#Phenotypic correlations
Data<-read.csv("Behavioral Assay Data.csv")
Data<-BehavioralAssayFinal

Foraging<-as.numeric(Data$Foraging)
Aggression<-as.numeric(Data$Aggression)
ExploratoryRate<-as.numeric(Data$ExploratoryRate)
ColonyExploration<-as.numeric(Data$ColonyExploration)
GroupExploration<-as.numeric(Data$GroupExploration)


cor.test(Foraging,Aggression,method="spearman")
cor.test(Foraging,ExploratoryRate,method="spearman")
cor.test(Foraging,GroupExploration,method="spearman")
cor.test(Foraging,ColonyExploration,method="spearman")

cor.test(Aggression,ExploratoryRate,method="spearman")
cor.test(Aggression,GroupExploration,method="spearman")
cor.test(Aggression,ColonyExploration,method="spearman")

cor.test(ExploratoryRate,GroupExploration,method="spearman")
cor.test(ExploratoryRate,ColonyExploration,method="spearman")

cor.test(GroupExploration,ColonyExploration,method="spearman")

p.adjust(pvalues, method = "fdr")

###############################################################################
#Genetic correlations

#Genetic correlations

mass<-read.csv("Behavioral Assay Data.csv")
pedigree<-read.csv("Pedigree.csv")

mass1=merge(mass,pedigree,by.x="Genotype",by.y="animal",all.X=TRUE)
mass1$animal=interaction(mass1$Colony,mass1$Genotype)

pedigreeH5A=mass1[c("animal","dam","sire","sex")]
pedigree<-pedigree[c("animal","dam","sire","sex")]

#bind the portion of the pedigree from before H500 with the H500 portion where we have mass data
pedigree2=rbind(pedigreeH5A,pedigree)
#create the A matrix 
Mat_A <- makeS(prepPed(pedigree2),heterogametic="M",returnS=F)

prior<-
  list(R=list(V=diag(2)/2,nu=2),
       G=list(G1=list(V=diag(2)/2,nu=2)))

mass1$CasteRatio<-(mass1$TotalGyne/(mass1$WPCasteRatio+mass1$TotalGyne))
mass1$SexRatio<-(mass1$TotalGyne)/(mass1$TotalGyne+mass1$TotalMale)

mass1$phen1<-as.numeric(mass1$Foraging)
mass1$phen2<-as.numeric(mass1$Aggression)

hist(mass1$phen1)
hist(mass1$phen2)

MCMC<-MCMCglmm(cbind(phen1,phen2)~trait-1,random=~us(trait):animal,
               rcov=~us(trait):units,
               ginverse=list(animal=Mat_A$Sinv),
               nitt = 1e8, burnin = 1e4, thin =100, family = c(rep("gaussian",1),rep("gaussian",1)), pl = T , pr = T, data = mass1,prior = prior)
summary(MCMC)

herit1<-MCMC$VCV[,'traitphen1:traitphen1.animal']/
  (MCMC$VCV[,'traitphen1:traitphen1.animal']+MCMC$VCV[,'traitphen1:traitphen1.units'])
herit2<-MCMC$VCV[,'traitphen2:traitphen2.animal']/
  (MCMC$VCV[,'traitphen2:traitphen2.animal']+MCMC$VCV[,'traitphen2:traitphen2.units'])
mean(herit1)
mean(herit2)

corr.gen<-MCMC$VCV[,'traitphen1:traitphen2.animal']/
  (sqrt(MCMC$VCV[,'traitphen1:traitphen1.animal']*MCMC$VCV[,'traitphen2:traitphen2.animal'])) 
mean(corr.gen) #Correlation
posterior.mode(corr.gen) 

HPDinterval(MCMC$VCV[,'traitphen1:traitphen2.animal']/
              (sqrt(MCMC$VCV[,'traitphen1:traitphen1.animal']*MCMC$VCV[,'traitphen2:traitphen2.animal'])))#95% CIs

#Repeat for each trait combination

######################################################################################
#Selection

library(mgcv)
library(gsg)
library(AER)

library(mgcv)
library(gsg)
library(AER)
library(car)
Data<-read.csv("BehavioralAssayFinal.csv")
#Use raw counts for sexual and worker pupae (no need to standardize)
Data$Block<-as.factor(Data$Block)
Data$SP<-as.numeric(Data$TotalSexualP)
Data$WP<-as.numeric(Data$TotalWP)

Data$Foraging<-as.numeric((Data$Foraging))
Data$Aggression<-as.numeric(Data$Aggression)
Data$ExRate<-as.numeric(Data$ExploratoryRate)
Data$GroupEx<-as.numeric(Data$GroupExploration)
Data$ColonyEx<-as.numeric(Data$GroupExploration)

#Remove rows with NAs
Data<-Data[!is.na(Data$Foraging),]
Data<-Data[!is.na(Data$Aggression),]
Data<-Data[!is.na(Data$ExRate),]
Data<-Data[!is.na(Data$GroupEx),]
Data<-Data[!is.na(Data$ColonyEx),]

#Phenotypic data are mean-centred and variance-standardized 

Data$z1<-as.numeric(Data$Foraging)
Data$z1<-(Data$z1-mean(Data$z1,na.rm=T))/sd(Data$z1)

Data$z2<-as.numeric(Data$Aggression)
Data$z2<-(Data$z2-mean(Data$z2,na.rm=T))/sd(Data$z2,na.rm=T)

Data$z3<-as.numeric(Data$ExRate)
Data$z3<-(Data$z3-mean(Data$z3,na.rm=T))/sd(Data$z3,na.rm=T)

Data$z4<-as.numeric(Data$GroupEx)
Data$z4<-(Data$z4-mean(Data$z4,na.rm=T))/sd(Data$z4,na.rm=T)

Data$z5<-as.numeric(Data$ColonyEx)
Data$z5<-(Data$z5-mean(Data$z5,na.rm=T))/sd(Data$z5,na.rm=T)
############################################
#Data$CasteRatio<-(Data$TotalGyne/(Data$WPCasteRatio+Data$TotalGyne))
#Data$SexRatio<-(Data$TotalGyne)/(Data$TotalGyne+Data$TotalMale)


#### MODEL FOR WORKER PUPAE ########
### use negative binomial family ####
ff<-gam(WP~s(z1)+s(z2)+s(z3)+s(z4)+s(z5)+Block,family="nb",data=Data)
#diagnostic (check that the edf value reported is not close to K'. The more the p-value for k-index is equal to 1, the better, but not always the case)
gam.check(ff) 
#### SELECTION GRADIENTS FOR WORKER PUPAE ###
fit<-gam.gradients(mod=ff,phenotype=c("z1","z2","z3","z4","z5"), covariates=c("Block"),standardized =F, 
                   se.method = "posterior",n.boot=1000,
                   refit.smooth = T)
fit

####################################
#### MODEL FOR SEXUAL PUPAE ########
### use tweedie family ####
ff<-gam(Data$SP~s(z1)+s(z2)+s(z3)+s(z4)+s(z5)+s(Block,bs="re"),family="tw",data=Data)
#adjust the par to plot the four diagnostic plots in one page

#diagnostic (check that the edf value reported is not close to K'. The more the p-value for k-index is equal to 1, the better, but not always the case)
gam.check(ff) 
#### SELECTION GRADIENTS FOR SEXUAL PUPAE ###
fit2<-gam.gradients(mod=ff,phenotype=c("z1","z2","z3","z4","z5"), covariates=c("Block"), standardized =F, 
                    se.method = "posterior",n.boot=1000,
                    refit.smooth = T)

fit2

