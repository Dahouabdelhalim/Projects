#Libraries to upload
library(gtools)
library(ggplot2)
library(gplots)
library(lme4)
library(lmerTest)
library(lsmeans)
library(scales)
library(vioplot)
library(multcomp)
library(pbkrtest)

### Data Files to Upload

#Part B - mothers' data (maternal treatments)
#Upload data 'SexySperm.PartB', name data SB1B

#Part C - sons' data (competitive double matings of sons) -> filled in variables for model
#Upload data 'SexySperm.Partc', name data SB1C.model

#Part C - sons' data (competitive double matings of sons) -> original data, includes missing variables 
#Upload data 'SexySperm.Partc', name data SB1C

#Part C - long-form male data (for post-hoc analyses of male weights)
#Upload data 'SexySperm.PartC.Males', name data SB1C.m


#########################################   SUMMARY DATA   ######################################## 

#Check Sample Sizes
table(SB1C$MatingOrder)
table(SB1C$SterilizationOrder)
table(SB1C$IntermatingInterval)
table(SB1C$MatingOrder,SB1C$SterilizationOrder)

#Are male ages and weights collinear? 

ages=lm(M1.MaleAge~M2.MaleAge,data=SB1C)
summary(ages)

weights=lm(M1.MaleWeight1~M2.MaleWeight1,data=SB1C)
summary(weights)


#Are males in part C size- and age-matched?
#Need to consider ages and sizes at the time point 0 (prior to mating), otherwise, differences based on when the second male actually mated (which is problematic for 24 and 48 hour inter-matings intervals after time elapse); to this end, I created two new columns to use for male 2: M2.MaleAge.Original and M2.MaleWeight1.Original  

#MALE SIZES
#paired t-tests
t.test(SB1C.model$M1.MaleWeight1, SB1C.model$M2.MaleWeight1.Original, paired=TRUE) #p = 0.7754

#MALE AGES
#paired t-tests
t.test(SB1C.model$M1.MaleAge, SB1C.model$M2.MaleAge.Original, paired=TRUE) #p = 0.7086


#########################################   PATERNITY SUCCESS   ######################################## 


### Table 1: Manually extract data within "SB1C.model" Excel data sheet  

#Statistics

#First, create an observation level random effect (OLRE) to deal with overdispersion
SB1C.model$OLRE = as.factor(1:dim(SB1C.model)[1])

#Second, determine random factors for model

#If get an error message when running the model that "...model failed to converge...", then include the following after family=binomial: ,control=glmerControl(optimizer="bobyqa")

P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | Sisters)+(1 | FemaleGeneration)+(1 | FemaleGeneration: FemaleMatriline)+(1 | M1.Brothers)+(1 | M1.MaleGeneration)+(1 | M1.MaleGeneration: M1.MaleMatriline)+(1 | M2.Brothers)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | Group)+(1 | OLRE),data= SB1C.model,family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(P2.glmer) #AIC 463.3
#Remove +(1 | M1.MaleGeneration: M1.MaleMatriline)
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | Sisters)+(1 | FemaleGeneration)+(1 | FemaleGeneration: FemaleMatriline)+(1 | M1.Brothers)+(1 | M1.MaleGeneration)+(1 | M2.Brothers)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | Group)+(1 | OLRE),data= SB1C.model,family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(P2.glmer) #AIC 461.3
#Remove +(1 | FemaleGeneration: FemaleMatriline)
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | Sisters)+(1 | FemaleGeneration)+(1 | M1.Brothers)+(1 | M1.MaleGeneration)+(1 | M2.Brothers)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | Group)+(1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer) #AIC 459.3
#Remove +(1 | Sisters)
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemaleGeneration)+(1 | M1.Brothers)+(1 | M1.MaleGeneration)+(1 | M2.Brothers)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | Group)+(1 | OLRE),data= SB1C.model,family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(P2.glmer) #AIC 457.3
#Remove +(1 | M1.MaleGeneration)
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemaleGeneration)+(1 | M1.Brothers)+(1 | M2.Brothers)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | Group)+(1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer) #AIC 455.3
#Remove +(1 | M2.MaleGeneration)
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemaleGeneration)+(1 | M1.Brothers)+(1 | M2.Brothers)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | Group)+(1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer) #AIC 453.3
#Remove +(1 | Group)
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemaleGeneration)+(1 | M1.Brothers)+(1 | M2.Brothers)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer) #AIC 451.3
#Remove +(1 | M2.MaleGeneration: M2.MaleMatriline)
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemaleGeneration)+(1 | M1.Brothers)+(1 | M2.Brothers)+(1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer) #AIC 449.3
#Remove +(1 | M1.Brothers)
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemaleGeneration)+(1 | M2.Brothers)+(1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer) #AIC 447.3
#Remove +(1 | M2.Brothers)
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemaleGeneration)+(1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer) #AIC 446.0
#Remove +(1 | OLRE)
P2.noOLRE.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemaleGeneration),data= SB1C.model,family=binomial)
summary(P2.noOLRE.glmer) #AIC 700.3

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
overdisp_fun(P2.noOLRE.glmer)  #YES, there is overdispersion, so keep OLRE in as a random factor

#Remove +(1 | FemaleGeneration)
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer) #AIC 449.1 (worse)

#Final random effects for model: 
P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  +(1 | FemaleGeneration)+(1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer)  

#Run bivariate analyses

#To do so, rescale the following factors:
SB1C.model$M2.EjaculateSizeDiff2 <- scale(SB1C.model$M2.EjaculateSizeDiff, center = TRUE, scale = TRUE)
SB1C.model$M2.CopulationLatencyDiff2 <- scale(SB1C.model$M2.CopulationLatencyDiff, center = TRUE, scale = TRUE)
SB1C.model$M2.KickingLatencyDiff2 <- scale(SB1C.model$M2.KickingLatencyDiff, center = TRUE, scale = TRUE)
SB1C.model$M2.KickingDurationDiff2 <- scale(SB1C.model$M2.KickingDurationDiff, center = TRUE, scale = TRUE)
SB1C.model$M2.CopulationDurationDiff2 <- scale(SB1C.model$M2.CopulationDurationDiff, center = TRUE, scale = TRUE)

P2.glmer=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  Maternal.IMI.Comparison +(1 | FemaleGeneration)+(1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer)  

#MatingOrder, 0.0778
#as.factor(IntermatingInterval), 0.0034**, 0.0538
#SterilizationOrder, 						0.45904
#M1.FemaleAge, 								0.20470   	
#M1.FemaleWeight1    						0.665 
#M2.MaleAgeDiff, 0.0887		   
#M2.MaleWeightDiff							0.930122  
#M2.EjaculateSizeDiff2, 					0.929247
#M2.CopulationLatencyDiff2 					0.491040
#M2.KickingLatencyDiff2						0.732027  
#M2.KickingDurationDiff2					0.729446  
#M2.CopulationDurationDiff2, 				0.509822
#Maternal.IMI.Comparison					0.57505
#Maternal.Mating.Number.Comparison, 		0.3870
#Maternal.Mating.Status.Comparison 			0.62474
#as.factor(Maternal.Treatment.Comparison)   warning message (too many variables), but all above 0.20


#Predictors to consider for model: 
#MatingOrder + as.factor(IntermatingInterval) + SterilizationOrder + M2.MaleAgeDiff 

#Check for collinearities for these variables: 
lm1=lm(M2.MaleAgeDiff~as.factor(IntermatingInterval),data=SB1C.model)
summary(lm1)

#The following variables are collinear:
#M2.MaleAgeDiff~as.factor(IntermatingInterval) (so only include IntermatingInterval) 

#Final model selection (AIC criteria <2)
P2.glmer1=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  MatingOrder + as.factor(IntermatingInterval) + SterilizationOrder +(1 | FemaleGeneration)+ (1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer1) #AIC 440.6
#Remove SterilizationOrder
P2.glmer2=glmer(cbind(Male2EggsRounded, Male1EggsRounded) ~  MatingOrder + as.factor(IntermatingInterval) +(1 | FemaleGeneration)+ (1 | OLRE),data= SB1C.model,family=binomial)
summary(P2.glmer2) #AIC 439.0
#Remove +(1 | FemaleGeneration), AIC 441.1 so worse

#P2.glmer2 is final model!

#Least squares means adjusted values: 
lsmeans(P2.glmer2, pairwise ~ MatingOrder,type="response")

#Calculate mean P2 by summing eggs laid across second males for each PM and MP

### FIGURE 2
PM<-SB1C.model$P2[SB1C.model$MatingOrder=="PM"]
MP<-SB1C.model$P2[SB1C.model$MatingOrder=="MP"]
vioplot(PM,MP,names=c("PM","MP"),col="gainsboro")
#Use "gainsboro" for MP 

#####Table 2 Results
#(Intercept)                       0.08548    0.24778   0.345  0.73012
#MatingOrderPM                     0.46890    0.22320   2.101  0.03566 *
#as.factor(IntermatingInterval)24  0.81943    0.25540   3.208  0.00133 **
#as.factor(IntermatingInterval)48  0.51352    0.29132   1.763  0.07794

#Mean +/- SE P2 based on mating order
se = function(x) sd(x)/sqrt(length(x))
#MP
mean(SB1C.model$P2[SB1C.model$MatingOrder=="MP"])
se(SB1C.model$P2[SB1C.model$MatingOrder=="MP"])
#PM
mean(SB1C.model$P2[SB1C.model$MatingOrder=="PM"])
se(SB1C.model$P2[SB1C.model$MatingOrder=="PM"])



##############################   Post-hoc Analyses for Mechanisms of Paternity Bias   #######################

######TABLE 3 DATA - Statistics for proportion remating females
#Manually look at data within 'SexySperm.PropRemating.csv' 

#IMI 0 
yes = c(13,13)
no = c(44,43)
 tab = as.data.frame(rbind(yes, no))
names(tab) = c("mp","pm")
tab
chisq.test(tab)
#IMI 24 
yes = c(15,10)
no = c(29,33)
 tab = as.data.frame(rbind(yes, no))
names(tab) = c("mp","pm")
tab
chisq.test(tab)
#IMI 48 
yes = c(8,10)
no = c(21,23)
 tab = as.data.frame(rbind(yes, no))
names(tab) = c("mp","pm")
tab
chisq.test(tab)

#TOTAL
yes = c(36,33)
no = c(21,23)
tab = as.data.frame(rbind(yes, no))
names(tab) = c("mono","poly")
tab
chisq.test(tab)


#####TABLE 4 DATA

######### Statistics for Mating 1: 


######### Ejaculate Size

#Check data distribution
hist(SB1C$M1.MaleEjacSize.Raw) #normally distributed

#Determine random factors for model
ejac1.lmer1=lmer(M1.MaleEjacSize.Raw ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | Group),data= SB1C)
summary(ejac1.lmer1)
#remove +(1 | Group) 
ejac1.lmer2=lmer(M1.MaleEjacSize.Raw ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline),data= SB1C)  
anova(ejac1.lmer1,ejac1.lmer2) #Not significantly different
summary(ejac1.lmer2)
#remove +(1 | FemaleGeneration:FemaleMatriline)
ejac1.lmer3=lmer(M1.MaleEjacSize.Raw ~ +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline),data= SB1C)
anova(ejac1.lmer3,ejac1.lmer2)  #Not significantly different
summary(ejac1.lmer3)
#remove +(1 | M1.MaleGeneration:M1.MaleMatriline)
ejac1.lmer4=lmer(M1.MaleEjacSize.Raw ~ +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline),data= SB1C)
anova(ejac1.lmer3,ejac1.lmer4)  #Not significantly different
summary(ejac1.lmer4)
#remove +(1 | M1.MaleMatriline)
ejac1.lmer5=lmer(M1.MaleEjacSize.Raw ~ +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration),data= SB1C)
anova(ejac1.lmer5,ejac1.lmer4)  #Not significantly different
summary(ejac1.lmer5)
#remove +(1 | M1.MaleGeneration)
ejac1.lmer6=lmer(M1.MaleEjacSize.Raw ~ +(1 | FemaleGeneration)+(1 | FemaleMatriline),data= SB1C)
anova(ejac1.lmer5,ejac1.lmer6)  #Not significantly different
summary(ejac1.lmer6)
#remove +(1 | FemaleGeneration)
ejac1.lmer7=lmer(M1.MaleEjacSize.Raw ~ +(1 | FemaleMatriline),data= SB1C)
anova(ejac1.lmer7,ejac1.lmer6)  #Not significantly different
summary(ejac1.lmer7)
#Try switching out FemaleMatriline for FemaleGeneration to see if they differ; if they don't, then revert to a LM instead of a LMER
ejac1.lmer8=lmer(M1.MaleEjacSize.Raw ~ +(1 | FemaleGeneration),data= SB1C)
anova(ejac1.lmer7,ejac1.lmer8)  #Not significantly different
summary(ejac1.lmer8)

#Thus, revert to a simple LM instead of a LMER
 
#Run bivariate analyses
ejac1.lm9=lm(M1.MaleEjacSize.Raw ~ MatingOrder,data=SB1C)
summary(ejac1.lm9) 

#MatingOrder, 0.145 
#SterilizationOrder   			0.896
#M1.FemaleWeight1   			0.8027
#M1.FemaleAge       			0.445
#M1.MaleWeight1, 0.000304 ***	
#M1.CopulationLatency   		0.743
#M1.KickingLatency 				0.384
#M1.KickingDuration  			0.619
#M1.CopulationDuration  		0.712

#Predictors to consider: 
#MatingOrder+M1.MaleWeight1

#Check for collinearities for these variables:
lm1=lm(M1.MaleWeight1 ~ MatingOrder, data=SB1C)
summary(lm1)

#The following variables are collinear:
#M1.MaleWeight1 ~ MatingOrder (keep Mating Order)
#boxplot(M1.MaleWeight1 ~ MatingOrder, data=SB1C)


#Final model selection (AIC criteria <2)
ejac1.lm9=lm(M1.MaleEjacSize.Raw ~ MatingOrder,data=SB1C)
summary(ejac1.lm9) 

#Check normality and homogeneity of residuals
plot(ejac1.lm9) #normal and homogeneous

#Mean +/- SE for absolute ejaculate size for males in first mating role
se = function(x) sd(x)/sqrt(length(x))
#MP
mean(SB1C$M1.MaleEjacSize.Raw[SB1C$MatingOrder=="MP"])
se(SB1C$M1.MaleEjacSize.Raw[SB1C$MatingOrder=="MP"])
#PM
mean(SB1C$M1.MaleEjacSize.Raw[SB1C$MatingOrder=="PM"])
se(SB1C$M1.MaleEjacSize.Raw[SB1C$MatingOrder=="PM"])


######### Copulation Latency 

#Check data distribution
hist(SB1C$M1.CopulationLatency) #NOT normally distributed

#Try a log transformation
hist(log(SB1C$M1.CopulationLatency)) #normally distributed

#Determine random factors for model
CL1.lmer1=lmer(log(M1.CopulationLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | Group),data= SB1C)
summary(CL1.lmer1)   
#remove +(1 | M1.MaleGeneration:M1.MaleMatriline)
CL1.lmer2=lmer(log(M1.CopulationLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | Group),data= SB1C)
anova(CL1.lmer1, CL1.lmer2)   #Not significantly different
summary(CL1.lmer2)
#remove +(1 | Group)
CL1.lmer3=lmer(log(M1.CopulationLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline),data= SB1C)
anova(CL1.lmer3, CL1.lmer2)   #Not significantly different
summary(CL1.lmer3)
#remove  +(1 | M1.MaleGeneration)
CL1.lmer4=lmer(log(M1.CopulationLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleMatriline),data= SB1C)
anova(CL1.lmer3, CL1.lmer4)   #Not significantly different
summary(CL1.lmer4)
#remove +(1 | M1.MaleMatriline) 
CL1.lmer5=lmer(log(M1.CopulationLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline),data= SB1C)
anova(CL1.lmer5, CL1.lmer4)   #Not significantly different
summary(CL1.lmer5)
#remove +(1 | FemaleMatriline)
CL1.lmer6=lmer(log(M1.CopulationLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration),data= SB1C)
anova(CL1.lmer5, CL1.lmer6)   #Not significantly different
summary(CL1.lmer6)
#remove +(1 | FemaleGeneration:FemaleMatriline)
CL1.lmer7=lmer(log(M1.CopulationLatency) ~  +(1 | FemaleGeneration),data= SB1C)
anova(CL1.lmer7, CL1.lmer6)   #Not significantly different
summary(CL1.lmer7)
#Try switching out FemaleGeneration for FemaleGeneration:FemaleMatriline to see if they differ; if they don't, then revert to a LM instead of a LMER
CL1.lmer8=lmer(log(M1.CopulationLatency) ~  +(1 | FemaleGeneration:FemaleMatriline),data= SB1C)
anova(CL1.lmer7, CL1.lmer8)   #Not significantly different
summary(CL1.lmer8)
#Thus, revert to a simple LM instead of a LMER

#Run bivariate analyses
CL1.lm8=lm(log(M1.CopulationLatency) ~ M1.MaleWeight1,data= SB1C)
summary(CL1.lm8)

#MatingOrder      					0.514
#SterilizationOrder, 0.163 
#M1.FemaleAge, 0.00907 **  
#M1.FemaleWeight1     				0.361        
#M1.MaleAge        					0.256    
#M1.MaleWeight1, 0.105  
 
#Predictors to consider: 
#MatingOrder + SterilizationOrder + M1.FemaleAge + M1.MaleWeight1

#Check for collinearities for these variables:
lm1=lm(M1.MaleWeight1 ~ SterilizationOrder, data=SB1C)
summary(lm1)

#The following variables are collinear: None

#Final model selection (AIC criteria <2)
CL1.lm8=lm(log(M1.CopulationLatency) ~ MatingOrder + SterilizationOrder + M1.FemaleAge + M1.MaleWeight1,data= SB1C)
summary(CL1.lm8)
#Remove M1.MaleWeight1
CL1.lm9=lm(log(M1.CopulationLatency) ~ MatingOrder + SterilizationOrder + M1.FemaleAge,data= SB1C)
anova(CL1.lm8, CL1.lm9)   #Not significantly different
summary(CL1.lm9)
#Remove SterilizationOrder
CL1.lm10=lm(log(M1.CopulationLatency) ~ MatingOrder + M1.FemaleAge,data= SB1C)
anova(CL1.lm10, CL1.lm9)   #Not significantly different
summary(CL1.lm10)

#Check model residuals
plot(CL1.lm10) #homogeneous

#MP: 
#Mean +/- SE for copulation latency for males in first mating role
mean(SB1C$M1.CopulationLatency[SB1C$MatingOrder=="MP"],na.rm=TRUE)
#Alternative calculation for SE: 
l=31
sd=sd(SB1C$M1.CopulationLatency[SB1C$MatingOrder=="MP"],na.rm=TRUE)
se=(sd/sqrt(l))
se

#PM:
mean(SB1C$M1.CopulationLatency[SB1C$MatingOrder=="PM"],na.rm=TRUE)
#Alternative calculation for SE:
l=29
sd=sd(SB1C$M1.CopulationLatency[SB1C$MatingOrder=="PM"],na.rm=TRUE)
se=(sd/sqrt(l))
se


######### Kicking Latency

#Check data distribution
hist(SB1C$M1.KickingLatency) #NOT normally distributed
#Try log transformation
hist(log(SB1C$M1.KickingLatency)) #NOT normally distributed
#Try log transformation
hist(log2(SB1C$M1.KickingLatency)) #Normally distributed

#Determine random factors for model
KL1.lmer1=lmer(log2(M1.KickingLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | Group),data= SB1C)
summary(KL1.lmer1)
#remove +(1 | Group)
KL1.lmer2=lmer(log2(M1.KickingLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline),data= SB1C)
anova(KL1.lmer1,KL1.lmer2) #Not significantly different
summary(KL1.lmer2)
#remove +(1 | M1.MaleGeneration)
KL1.lmer3=lmer(log2(M1.KickingLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleMatriline),data= SB1C)
anova(KL1.lmer3,KL1.lmer2) #Not significantly different
summary(KL1.lmer3)
#remove +(1 | M1.MaleGeneration:M1.MaleMatriline)
KL1.lmer4=lmer(log2(M1.KickingLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleMatriline),data=SB1C)
anova(KL1.lmer3,KL1.lmer4) #Not significantly different
summary(KL1.lmer4)
#remove +(1 | M1.MaleMatriline) 
KL1.lmer5=lmer(log2(M1.KickingLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline),data=SB1C)
anova(KL1.lmer5,KL1.lmer4) #Not significantly different
summary(KL1.lmer5)
#remove +(1 | FemaleMatriline)
KL1.lmer6=lmer(log2(M1.KickingLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration),data=SB1C)
anova(KL1.lmer5,KL1.lmer6) #Not significantly different
summary(KL1.lmer6)
#remove +(1 | FemaleGeneration:FemaleMatriline) 
KL1.lmer7=lmer(log2(M1.KickingLatency) ~ +(1 | FemaleGeneration),data=SB1C)
anova(KL1.lmer7,KL1.lmer6) #Not significantly different
summary(KL1.lmer7)
#Try switching out FemaleGeneration for FemaleGeneration:FemaleMatriline to see if they differ; if they don't, then revert to a LM instead of a LMER
KL1.lmer8=lmer(log2(M1.KickingLatency) ~ +(1 | FemaleGeneration:FemaleMatriline),data=SB1C)
anova(KL1.lmer7,KL1.lmer8) #Significantly different, so maintain as a LMER with FemaleGeneration only

#Run bivariate analyses
KL1.lmer7=lmer(log2(M1.KickingLatency) ~ M1.MaleWeight1 +(1 | FemaleGeneration),data= SB1C)
summary(KL1.lmer7)

#MatingOrder, 0.106
#SterilizationOrder       		0.596
#M1.FemaleAge           		0.378   
#M1.FemaleWeight1      		 	0.815      
#M1.MaleAge          			0.924  
#M1.MaleWeight1      			0.526

#Predictors to consider: 
#MatingOrder

#Final model
KL1.lmer7=lmer(log2(M1.KickingLatency) ~ MatingOrder +(1 | FemaleGeneration),data= SB1C)
summary(KL1.lmer7)

#Check model residuals
plot(KL1.lmer7) #homogeneous


#Mean +/- SE kicking latency for males in first mating role
#MP:
mean(SB1C$M1.KickingLatency[SB1C$MatingOrder=="MP"],na.rm=TRUE)
#Alternative calculation for SE:  
l=34
sd=sd(SB1C$M1.KickingLatency[SB1C$MatingOrder=="MP"],na.rm=TRUE)
se=(sd/sqrt(l))
se
#PM:
mean(SB1C$M1.KickingLatency[SB1C$MatingOrder=="PM"],na.rm=TRUE)
#Alternative calculation for SE:
l=32
sd=sd(SB1C$M1.KickingLatency[SB1C$MatingOrder=="PM"],na.rm=TRUE)
se=(sd/sqrt(l))
se



######### Kicking Duration

#Check data distribution
hist(SB1C$M1.KickingDuration) #NOT normally distributed
#Try log transformation
hist(log(SB1C$M1.KickingDuration)) #normally distributed

#Determine random factors for model
KD1.lmer1=lmer(log(M1.KickingDuration) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | Group),data= SB1C)
summary(KD1.lmer1)
#remove +(1 | M1.MaleMatriline)
KD1.lmer2=lmer(log(M1.KickingDuration) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | Group),data= SB1C)
anova(KD1.lmer1, KD1.lmer2)  #Not significantly different
summary(KD1.lmer2)
#remove +(1 | FemaleGeneration:FemaleMatriline)
KD1.lmer3=lmer(log(M1.KickingDuration) ~ +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | Group),data= SB1C)
anova(KD1.lmer3, KD1.lmer2)  #Not significantly different
summary(KD1.lmer3)
#remove +(1 | M1.MaleGeneration)
KD1.lmer4=lmer(log(M1.KickingDuration) ~ +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | Group),data= SB1C)
anova(KD1.lmer3, KD1.lmer4)  #Not significantly different
summary(KD1.lmer4)
#remove +(1 | FemaleMatriline)
KD1.lmer5=lmer(log(M1.KickingDuration) ~ +(1 | FemaleGeneration)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | Group),data= SB1C)
anova(KD1.lmer5, KD1.lmer4)  #Not significantly different
summary(KD1.lmer5)
#remove +(1 | FemaleGeneration)
KD1.lmer6=lmer(log(M1.KickingDuration) ~ +(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | Group),data= SB1C)
anova(KD1.lmer5, KD1.lmer6)  #Not significantly different
summary(KD1.lmer6)
#remove +(1 | Group)
KD1.lmer7=lmer(log(M1.KickingDuration) ~ +(1 | M1.MaleGeneration:M1.MaleMatriline),data= SB1C)
anova(KD1.lmer7, KD1.lmer6)  #Not significantly different
summary(KD1.lmer7)
#remove nested random factor
KD1.lmer8=lmer(log(M1.KickingDuration) ~ +(1 | M1.MaleGeneration) +(1 | M1.MaleMatriline),data= SB1C)
anova(KD1.lmer7, KD1.lmer8)  #Not significantly different
summary(KD1.lmer8)
#remove +(1 | M1.MaleGeneration)
KD1.lmer9=lmer(log(M1.KickingDuration) ~  +(1 | M1.MaleMatriline),data= SB1C)
anova(KD1.lmer9, KD1.lmer8)  #Not significantly different
summary(KD1.lmer9)
#Try switching out M1.MaleMatriline for M1.MaleGeneration to see if they differ; if they don't, then revert to a LM instead of a LMER
KD1.lmer10=lmer(log(M1.KickingDuration) ~  +(1 | M1.MaleGeneration),data= SB1C)
anova(KD1.lmer9, KD1.lmer10)  #Yes, significantly different, so maintain LMER with M1.MaleMatriline only

#Run bivariate analyses
KD1.lmer9=lmer(log(M1.KickingDuration) ~ M1.MaleWeight1 +(1 | M1.MaleMatriline),data= SB1C)
summary(KD1.lmer9)

#MatingOrder        			0.309 
#SterilizationOrder, 0.113
#M1.FemaleAge, 0.0397 *        			   		   
#M1.FemaleWeight1, 0.0088 **   
#M1.MaleAge          		  	0.757
#M1.MaleWeight1      			0.982

#Predictors to consider: 
#MatingOrder + SterilizationOrder + M1.FemaleAge + M1.FemaleWeight1

#Check for collinearities for these variables: 
lm1=lm(M1.FemaleAge ~ MatingOrder,data=SB1C)
summary(lm1)

#The following variables are collinear: NONE

#Final model selection (AIC criteria <2)
KD1.lmer9=lmer(log(M1.KickingDuration) ~ MatingOrder + SterilizationOrder + M1.FemaleAge + M1.FemaleWeight1 +(1 | M1.MaleMatriline),data= SB1C)
summary(KD1.lmer9)

#Check model residuals
plot(KD1.lmer9)  #homogeneous

#Mean +/- SE kicking duration for males in first mating role
#MP:
mean(SB1C$M1.KickingDuration[SB1C$MatingOrder=="MP"],na.rm=TRUE)
#Alternative calculation for SE:  
l=32
sd=sd(SB1C$M1.KickingDuration[SB1C$MatingOrder=="MP"],na.rm=TRUE)
se=(sd/sqrt(l))
se
#PM: 
mean(SB1C$M1.KickingDuration[SB1C$MatingOrder=="PM"],na.rm=TRUE)
l=32
sd=sd(SB1C$M1.KickingDuration[SB1C$MatingOrder=="PM"],na.rm=TRUE)
se=(sd/sqrt(l))
se



######### Copulation Duration

#Check data distribution
hist(SB1C$M1.CopulationDuration) #NOT normally distributed
#Try log transformation
hist(log(SB1C$M1.CopulationDuration)) #Normally distributed

#Determine random factors for model
CD1.lmer1=lmer(log(M1.CopulationDuration) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | Group),data=SB1C)
summary(CD1.lmer1)
#remove +(1 | M1.MaleMatriline)
CD1.lmer2=lmer(log(M1.CopulationDuration) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | Group),data=SB1C)   
anova(CD1.lmer1,CD1.lmer2) #Not significantly different
summary(CD1.lmer2)
#remove +(1 | FemaleGeneration:FemaleMatriline)
CD1.lmer3=lmer(log(M1.CopulationDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | Group),data=SB1C)  
anova(CD1.lmer3,CD1.lmer2)  #Not significantly different
summary(CD1.lmer3)
#remove +(1 | M1.MaleGeneration)
CD1.lmer4=lmer(log(M1.CopulationDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | Group),data=SB1C)  
anova(CD1.lmer3,CD1.lmer4)  #Not significantly different
summary(CD1.lmer4)
#remove +(1 | FemaleGeneration)
CD1.lmer5=lmer(log(M1.CopulationDuration) ~  +(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | Group),data=SB1C)  
anova(CD1.lmer5,CD1.lmer4)  #Not significantly different
summary(CD1.lmer5)
#remove +(1 | FemaleMatriline)
CD1.lmer6=lmer(log(M1.CopulationDuration) ~ +(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | Group),data=SB1C)  
anova(CD1.lmer5,CD1.lmer6)  #Not significantly different 
summary(CD1.lmer6)
#remove +(1 | Group)
CD1.lmer7=lmer(log(M1.CopulationDuration) ~ +(1 | M1.MaleGeneration:M1.MaleMatriline),data=SB1C)  
anova(CD1.lmer7,CD1.lmer6)  #Not significantly different 
summary(CD1.lmer7)
#remove nested random factor
CD1.lmer8=lmer(log(M1.CopulationDuration) ~ +(1 | M1.MaleGeneration) +(1 | M1.MaleMatriline),data=SB1C)  
anova(CD1.lmer7,CD1.lmer8)  #Not significantly different 
summary(CD1.lmer8)
#remove +(1 | M1.MaleGeneration)
CD1.lmer9=lmer(log(M1.CopulationDuration) ~  +(1 | M1.MaleMatriline),data=SB1C)  
anova(CD1.lmer9,CD1.lmer8)  #Not significantly different 
summary(CD1.lmer9)
#Try switching out M1.MaleMatriline for M1.MaleGeneration to see if they differ; if they don't, then revert to a LM instead of a LMER 
CD1.lmer10=lmer(log(M1.CopulationDuration) ~  +(1 | M1.MaleGeneration),data=SB1C)  
anova(CD1.lmer9,CD1.lmer10)  #Yes, significantly different, so keep lmer with M1.MaleMatriline only 

#Run bivariate analyses
CD1.lmer9=lmer(log(M1.CopulationDuration) ~ M1.MaleWeight1 +(1 | M1.MaleMatriline),data=SB1C)  
summary(CD1.lmer9)

#MatingOrder, 0.0537
#SterilizationOrder     		0.374 
#M1.FemaleAge, 0.0187 * 
#M1.FemaleWeight1, 0.084   
#M1.MaleAge        				0.927       
#M1.MaleWeight1    				0.357       

#Predictors to consider: 
#MatingOrder + M1.FemaleAge + M1.FemaleWeight1

#Check for collinearities for these variables: 
lm1=lm(M1.FemaleAge ~ M1.FemaleWeight1,data=SB1C)
summary(lm1)

#The following variables are collinear: NONE 

#Final model selection (AIC criteria <2)
CD1.lmer9=lmer(log(M1.CopulationDuration) ~ MatingOrder + M1.FemaleAge + M1.FemaleWeight1 +(1 | M1.MaleMatriline),data=SB1C)  
summary(CD1.lmer9)
#Remove M1.FemaleWeight1
CD1.lmer10=lmer(log(M1.CopulationDuration) ~ MatingOrder + M1.FemaleAge +(1 | M1.MaleMatriline),data=SB1C)
anova(CD1.lmer9,CD1.lmer10)  #Not significantly different 
summary(CD1.lmer10)

#Check model residuals
plot(CD1.lmer10) #homogeneous

#Mean +/- SE Copulation Duration for males in first mating role
#MP: 
mean(SB1C$M1.CopulationDuration[SB1C$MatingOrder=="MP"],na.rm=TRUE)
#Alternative calculation for SE: 
l=34
sd=sd(SB1C$M1.CopulationDuration[SB1C$MatingOrder=="MP"],na.rm=TRUE)
se=(sd/sqrt(l))
se

#PM:
mean(SB1C$M1.CopulationDuration[SB1C$MatingOrder=="PM"],na.rm=TRUE)
l=32
sd=sd(SB1C$M1.CopulationDuration[SB1C$MatingOrder=="PM"],na.rm=TRUE)
se=(sd/sqrt(l))
se




######### Statistics for Mating 2: 

######### Ejaculate size

#Check data distribution
hist(SB1C$M2.MaleEjacSize.Raw) #Not normally distributed
#Try a log transformation
hist(sqrt(SB1C$M2.MaleEjacSize.Raw)) #Normally distributed

#Determine random factors for model
ejac2.lmer1=lmer(sqrt(M2.MaleEjacSize.Raw) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data=SB1C) 
summary(ejac2.lmer1)  
#Remove +(1 | FemaleGeneration:FemaleMatriline)
ejac2.lmer2=lmer(sqrt(M2.MaleEjacSize.Raw) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data=SB1C)
anova(ejac2.lmer1,ejac2.lmer2)  #Not significantly different
summary(ejac2.lmer2) 
#Remove +(1 | M1.MaleGeneration:M1.MaleMatriline)
ejac2.lmer3=lmer(sqrt(M2.MaleEjacSize.Raw) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data=SB1C)
anova(ejac2.lmer3,ejac2.lmer2)  #Not significantly different
summary(ejac2.lmer3) 
#Remove +(1 | M2.MaleGeneration:M2.MaleMatriline)
ejac2.lmer4=lmer(sqrt(M2.MaleEjacSize.Raw) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data=SB1C)
anova(ejac2.lmer3,ejac2.lmer4)  #Not significantly different
summary(ejac2.lmer4)
#Remove +(1 | Group)
ejac2.lmer5=lmer(sqrt(M2.MaleEjacSize.Raw) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline),data=SB1C)
anova(ejac2.lmer5,ejac2.lmer4)  #Not significantly different
summary(ejac2.lmer5)
#Remove +(1 | FemaleGeneration)
ejac2.lmer6=lmer(sqrt(M2.MaleEjacSize.Raw) ~  +(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline),data=SB1C)
anova(ejac2.lmer5,ejac2.lmer6)  #Not significantly different
summary(ejac2.lmer6)
#Remove +(1 | M1.MaleGeneration)
ejac2.lmer7=lmer(sqrt(M2.MaleEjacSize.Raw) ~  +(1 | FemaleMatriline)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline),data=SB1C)
anova(ejac2.lmer7,ejac2.lmer6)  #Not significantly different
summary(ejac2.lmer7)
#Remove +(1 | M2.MaleMatriline)
ejac2.lmer8=lmer(sqrt(M2.MaleEjacSize.Raw) ~  +(1 | FemaleMatriline)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration),data=SB1C)
anova(ejac2.lmer7,ejac2.lmer8)  #Not significantly different
summary(ejac2.lmer8)
#Remove +(1 | FemaleMatriline)
ejac2.lmer9=lmer(sqrt(M2.MaleEjacSize.Raw) ~  +(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration),data=SB1C)
anova(ejac2.lmer9,ejac2.lmer8)  #Not significantly different
summary(ejac2.lmer9)
#Remove +(1 | M1.MaleMatriline)
ejac2.lmer10=lmer(sqrt(M2.MaleEjacSize.Raw) ~  +(1 | M2.MaleGeneration),data=SB1C)
anova(ejac2.lmer9,ejac2.lmer10)  #Not significantly different
summary(ejac2.lmer10)
#Try switching out M2.MaleGeneration for M1.MaleMatriline to see if they differ; if they don't, then revert to a LM instead of a LMER 
ejac2.lmer11=lmer(sqrt(M2.MaleEjacSize.Raw) ~  +(1 | M1.MaleMatriline),data=SB1C)
anova(ejac2.lmer11,ejac2.lmer10)  #Yes, significantly different, so maintain LMER with M2.MaleGeneration 

#Run bivariate analyses
ejac2.lmer10=lmer(sqrt(M2.MaleEjacSize.Raw) ~  M2.CopulationDurationDiff + (1 | M2.MaleGeneration),data=SB1C)
summary(ejac2.lmer10) 

#MatingOrder, 0.0616
#SterilizationOrder   					0.50365
#as.factor(IntermatingInterval)  		0.456, 0.740
#M1.FemaleWeight1, 0.182
#M1.FemaleAge, 0.252 
#M1.FemaleWeightGain   					0.991
#M1.MaleEjacSize.Raw  					0.455
#M1.CopulationLatency   				0.814785
#M1.KickingLatency 						0.247
#M1.KickingDuration 					0.386
#M1.CopulationDuration  				0.769
#M2.FemaleWeight1, 0.13134
#M2.MaleWeight1      					0.3147
#M2.MaleAge								0.967
#M2.CopulationLatency, 0.017 * 
#M2.KickingLatency						0.305
#M2.KickingDuration						0.845
#M2.CopulationDuration, 				0.238
#M2.CopulationDurationDiff				0.607062  

#Predictors to consider for model: 
#MatingOrder + M1.FemaleWeight1 + M1.FemaleAge + M2.FemaleWeight1 + M2.CopulationLatency 

#Check for collinearities for these variables: 
lm1=lm(M1.FemaleAge ~ M2.FemaleWeight1,data=SB1C)
summary(lm1)

#The following variables are collinear: 
#M1.FemaleWeight1~M2.FemaleWeight1 (keep M2.FemaleWeight1)

#Final model selection (AIC criteria <2)
ejac2.lmer10=lmer(sqrt(M2.MaleEjacSize.Raw) ~  MatingOrder + M1.FemaleAge + M2.FemaleWeight1 + M2.CopulationLatency + (1 | M2.MaleGeneration),data=SB1C)
summary(ejac2.lmer10) 
#revert to a LM (1 | M2.MaleGeneration) has no effect on the response variable

ejac2.lm10=lm(sqrt(M2.MaleEjacSize.Raw) ~  MatingOrder + M1.FemaleAge + M2.FemaleWeight1 + M2.CopulationLatency,data=SB1C)
summary(ejac2.lm10)
#Remove + M2.CopulationLatency
ejac2.lm11=lm(sqrt(M2.MaleEjacSize.Raw) ~  MatingOrder + M1.FemaleAge + M2.FemaleWeight1,data=SB1C)
anova(ejac2.lm10,ejac2.lm11) #cannot compare, different data set (sample sizes) 
summary(ejac2.lm11)
#remove + M1.FemaleAge
ejac2.lm12=lm(sqrt(M2.MaleEjacSize.Raw) ~  MatingOrder  + M2.FemaleWeight1,data=SB1C)
anova(ejac2.lm12,ejac2.lm11) #Not significantly different 
summary(ejac2.lm12)
#remove + M2.FemaleWeight1
ejac2.lm13=lm(sqrt(M2.MaleEjacSize.Raw) ~  MatingOrder,data=SB1C)
anova(ejac2.lm12,ejac2.lm13) #cannot compare, different data set (sample sizes)
summary(ejac2.lm13)


#Check model residuals
plot(ejac2.lm13) #homogeneous

#Mean +/- SE Ejaculate size for males in second mating role
#PM
mean(SB1C$M2.MaleEjacSize.Raw[SB1C$MatingOrder=="PM"],na.rm=TRUE)
l=32
sd=sd(SB1C$M2.MaleEjacSize.Raw[SB1C$MatingOrder=="PM"],na.rm=TRUE)
se=(sd/sqrt(l))
se

#MP
mean(SB1C$M2.MaleEjacSize.Raw[SB1C$MatingOrder=="MP"],na.rm=TRUE)
#Alternative calculation for SE:  
l=35
sd=sd(SB1C$M2.MaleEjacSize.Raw[SB1C$MatingOrder=="MP"],na.rm=TRUE)
se=(sd/sqrt(l))
se


######### Copulation Latency 

#Check data distribution
hist(SB1C$M2.CopulationLatency) #Not normally distributed
#Try log
hist(log(SB1C$M2.CopulationLatency)) #More normally distributed

#Determine random factors for model
CL2.lmer1=lmer(log(M2.CopulationLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
summary(CL2.lmer1)
#Remove +(1 | FemaleGeneration:FemaleMatriline)
CL2.lmer2=lmer(log(M2.CopulationLatency) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C) 
anova(CL2.lmer1, CL2.lmer2)  #Not significantly different
summary(CL2.lmer2) 
#Remove +(1 | M2.MaleGeneration:M2.MaleMatriline)
CL2.lmer3=lmer(log(M2.CopulationLatency) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C) 
anova(CL2.lmer3, CL2.lmer2)  #Not significantly different
summary(CL2.lmer3) 
#Remove +(1 | FemaleMatriline)
CL2.lmer4=lmer(log(M2.CopulationLatency) ~  +(1 | FemaleGeneration)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C) 
anova(CL2.lmer3, CL2.lmer4)  #Not significantly different
summary(CL2.lmer4) 
#Remove +(1 | Group)
CL2.lmer5=lmer(log(M2.CopulationLatency) ~  +(1 | FemaleGeneration)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline),data= SB1C) 
anova(CL2.lmer5, CL2.lmer4)  #Not significantly different
summary(CL2.lmer5) 
#Remove +(1 | FemaleGeneration)
CL2.lmer6=lmer(log(M2.CopulationLatency) ~  +(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline),data= SB1C) 
anova(CL2.lmer5, CL2.lmer6)  #Not significantly different
summary(CL2.lmer6) 
#Remove +(1 | M1.MaleGeneration)
CL2.lmer7=lmer(log(M2.CopulationLatency) ~  +(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline),data= SB1C) 
anova(CL2.lmer7, CL2.lmer6)  #Not significantly different
summary(CL2.lmer7) 
#Remove +(1 | M2.MaleGeneration)
CL2.lmer8=lmer(log(M2.CopulationLatency) ~  +(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleMatriline)+(1 | M2.MaleMatriline),data= SB1C) 
anova(CL2.lmer7, CL2.lmer8)  #Not significantly different
summary(CL2.lmer8) 
#Remove +(1 | M1.MaleMatriline)
CL2.lmer9=lmer(log(M2.CopulationLatency) ~  +(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M2.MaleMatriline),data= SB1C) 
anova(CL2.lmer9, CL2.lmer8)  #Not significantly different
summary(CL2.lmer9) 
#Remove +(1 | M1.MaleGeneration:M1.MaleMatriline)
CL2.lmer10=lmer(log(M2.CopulationLatency) ~  +(1 | M2.MaleMatriline),data= SB1C)
anova(CL2.lmer9, CL2.lmer10)  #Not significantly different
summary(CL2.lmer10) 
#Try switching out M2.MaleMatriline for M1.MaleGeneration:M1.MaleMatriline to see if they differ; if they don't, then revert to a LM instead of a LMER 
CL2.lmer11=lmer(log(M2.CopulationLatency) ~  +(1 | M1.MaleGeneration:M1.MaleMatriline),data= SB1C)
anova(CL2.lmer11, CL2.lmer10)  #Yes, significantly different, so maintain LMER with M2.MaleMatriline 

#Run bivariate analyses
CL2.lmer10=lmer(log(M2.CopulationLatency) ~ M2.MaleWeight1 +(1 | M2.MaleMatriline),data= SB1C) 
summary(CL2.lmer10)

#MatingOrder						0.265
#as.factor(IntermatingInterval), 0.2772, 0.0379 * 
#SterilizationOrder       			0.494
#M1.FemaleAge, 0.0474 *
#M1.FemaleWeight1          			0.399
#M1.FemaleWeightGain       			0.884
#M1.MaleAge, 0.0389 *
#M1.MaleWeight1         			0.999
#M1.MaleEjacSize.Raw         		0.891
#M1.CopulationLatency       		0.37
#M1.KickingLatency      			0.754
#M1.KickingDuration      			0.26
#M1.CopulationDuration, 0.167
#M2.FemaleWeightGain       			0.853 
#M2.MaleAge             			0.61  
#M2.MaleWeight1         			0.71   

#Predictors to consider: 
#MatingOrder + as.factor(IntermatingInterval) + M1.FemaleAge + M1.MaleAge + M1.CopulationDuration 

#Check for collinearities for these variables: 
lm1<-lm(M1.MaleAge ~ M1.CopulationDuration,data=SB1C)
summary(lm1)

#The following variables are collinear: 
#M1.FemaleAge ~ M1.CopulationDuration (keep M1.FemaleAge)

#Final model selection (AIC criteria <2)
CL2.lmer10=lmer(log(M2.CopulationLatency) ~ MatingOrder + as.factor(IntermatingInterval) + M1.FemaleAge + M1.MaleAge +(1 | M2.MaleMatriline),data= SB1C) 
summary(CL2.lmer10)
#remove + M1.FemaleAge
CL2.lmer11=lmer(log(M2.CopulationLatency) ~ MatingOrder + as.factor(IntermatingInterval)  + M1.MaleAge +(1 | M2.MaleMatriline),data= SB1C)
anova(CL2.lmer10,CL2.lmer11) #Error message, not significantly different 
summary(CL2.lmer11)
#Random factor +(1 | M2.MaleMatriline) has no effect on the response variable, so revert to a LM: 
CL2.lm12=lm(log(M2.CopulationLatency) ~ MatingOrder + as.factor(IntermatingInterval)  + M1.MaleAge,data= SB1C)
summary(CL2.lm12)

#Check model residuals
plot(CL2.lm12) #homogeneous

#Mean +/- SE XX Copulation Latency for males in second mating role
#PM:
mean(SB1C$M2.CopulationLatency[SB1C$MatingOrder=="PM"],na.rm=TRUE)
#Alternative calculation for SE: 
l=27
sd=sd(SB1C$M2.CopulationLatency[SB1C$MatingOrder=="PM"],na.rm=TRUE)
se=(sd/sqrt(l))
se

#MP:
mean(SB1C$M2.CopulationLatency[SB1C$MatingOrder=="MP"],na.rm=TRUE)
#Alternative calculation for SE: 
#MP: 
l=33
sd=sd(SB1C$M2.CopulationLatency[SB1C$MatingOrder=="MP"],na.rm=TRUE)
se=(sd/sqrt(l))
se



######### Kicking Latency 

#Check data distribution
hist(SB1C$M2.KickingLatency) #Not normally distributed
#Try log transformation
hist(log(SB1C$M2.KickingLatency)) #Normally distributed

#Determine random factors for model

KL2.lmer1=lmer(log(M2.KickingLatency) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
summary(KL2.lmer1)
#Remove +(1 | FemaleGeneration:FemaleMatriline)
KL2.lmer2=lmer(log(M2.KickingLatency)~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(KL2.lmer1,KL2.lmer2) #Not significantly different
summary(KL2.lmer2)
#Remove +(1 | M1.MaleGeneration:M1.MaleMatriline)
KL2.lmer3=lmer(log(M2.KickingLatency) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(KL2.lmer3,KL2.lmer2) #Not significantly different
summary(KL2.lmer3)
#Remove +(1 | M2.MaleGeneration:M2.MaleMatriline) 
KL2.lmer4=lmer(log(M2.KickingLatency) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(KL2.lmer3,KL2.lmer4) #Not significantly different
summary(KL2.lmer4)
#Remove +(1 | M1.MaleMatriline)
KL2.lmer5=lmer(log(M2.KickingLatency) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(KL2.lmer5,KL2.lmer4) #Not significantly different
summary(KL2.lmer5)
#Remove +(1 | FemaleGeneration)
KL2.lmer6=lmer(log(M2.KickingLatency) ~  +(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(KL2.lmer5,KL2.lmer6) #Not significantly different
summary(KL2.lmer6)
#Remove +(1 | M1.MaleGeneration)
KL2.lmer7=lmer(log(M2.KickingLatency) ~  +(1 | FemaleMatriline)+(1 | M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | Group),data= SB1C)
anova(KL2.lmer7,KL2.lmer6) #Not significantly different
summary(KL2.lmer7)
#Remove +(1 | Group)
KL2.lmer8=lmer(log(M2.KickingLatency) ~  +(1 | FemaleMatriline)+(1 | M2.MaleMatriline)+(1 | M2.MaleGeneration),data= SB1C)
anova(KL2.lmer7,KL2.lmer8) #Not significantly different
summary(KL2.lmer8)
#Remove +(1 | M2.MaleMatriline)
KL2.lmer9=lmer(log(M2.KickingLatency) ~  +(1 | FemaleMatriline)+(1 | M2.MaleGeneration),data= SB1C)
anova(KL2.lmer9,KL2.lmer8) #Not significantly different
summary(KL2.lmer9)
#Remove +(1 | FemaleMatriline)
KL2.lmer10=lmer(log(M2.KickingLatency) ~  +(1 | M2.MaleGeneration),data= SB1C)
anova(KL2.lmer9,KL2.lmer10) #Not significantly different
summary(KL2.lmer10)
#Try switching out M2.MaleGeneration for FemaleMatriline to see if they differ; if they don't, then revert to a LM instead of a LMER   
KL2.lmer11=lmer(log(M2.KickingLatency) ~  +(1 | FemaleMatriline),data= SB1C)
anova(KL2.lmer11,KL2.lmer10) #Not significantly different
summary(KL2.lmer11)

#Thus, revert to a lm

#Run bivariate analyses
KL2.lm11=lm(log(M2.KickingLatency) ~ M2.MaleWeight1,data= SB1C)
summary(KL2.lm11) 

#MatingOrder						0.926
#as.factor(IntermatingInterval), 0.000896 ***, 3.65e-07 ***
#SterilizationOrder, 0.00913 **  
#M1.FemaleAge           			0.531
#M1.FemaleWeight1           		0.913
#M1.FemaleWeightGain      			0.527
#M1.MaleAge            				0.268
#M1.MaleWeight1, 0.181
#M1.MaleEjacSize.Raw, 0.122 
#M1.CopulationLatency    			0.282
#M1.KickingLatency, 0.164 
#M1.KickingDuration, 0.121
#M1.CopulationDuration,0.0146 *  
#M2.FemaleWeightGain        		0.805
#M2.MaleAge, 0.00506 **
#M2.MaleWeight1, 0.0134 * 

#Predictors to consider:  
#MatingOrder + as.factor(IntermatingInterval) + SterilizationOrder + M1.MaleWeight1 + M1.MaleEjacSize.Raw + M1.KickingLatency + M1.KickingDuration + M1.CopulationDuration + M2.MaleAge + M2.MaleWeight1 +

#Check for collinearities for these variables:
lm1=lm(M1.CopulationDuration ~ M2.MaleWeight1,data=SB1C)
summary(lm1)

#The following variables are collinear:
#M1.MaleWeight1 ~ MatingOrder (Keep Mating Order, the main response of interest)
#M1.MaleEjacSize.Raw ~ as.factor(IntermatingInterval) (keep as.factor(IntermatingInterval))
#M1.KickingDuration ~ as.factor(IntermatingInterval) (keep as.factor(IntermatingInterval))
#M2.MaleAge ~ as.factor(IntermatingInterval) (keep as.factor(IntermatingInterval))
#M1.KickingLatency~ M1.CopulationDuration (keep M1.CopulationDuration)
 
#MatingOrder + as.factor(IntermatingInterval) + SterilizationOrder + M1.CopulationDuration + M2.MaleWeight1 +


#Final model selection (AIC criteria <2)
KL2.lm11=lm(log(M2.KickingLatency) ~ MatingOrder + as.factor(IntermatingInterval) + SterilizationOrder + M1.CopulationDuration + M2.MaleWeight1,data= SB1C)
summary(KL2.lm11)  
#Remove + SterilizationOrder
KL2.lm12=lm(log(M2.KickingLatency) ~ MatingOrder + as.factor(IntermatingInterval) + M1.CopulationDuration + M2.MaleWeight1,data= SB1C)
anova(KL2.lm11,KL2.lm12) #Not significantly different
summary(KL2.lm12) 

#Check model residuals
plot(KL2.lm12) #homogeneous

#Mean +/- SE Kicking Latency for males in second mating role
#PM: 
mean(SB1C$M2.KickingLatency[SB1C$MatingOrder=="PM"],na.rm=TRUE)
#Alternative calculation for SE: 
l=29
sd=sd(SB1C$M2.KickingLatency[SB1C$MatingOrder=="PM"],na.rm=TRUE)
se=(sd/sqrt(l))
se

#MP
mean(SB1C$M2.KickingLatency[SB1C$MatingOrder=="MP"],na.rm=TRUE)
#Alternative calculation for SE: 
l=34
sd=sd(SB1C$M2.KickingLatency[SB1C$MatingOrder=="MP"],na.rm=TRUE)
se=(sd/sqrt(l))
se



######### Kicking Duration

#Check data distribution
hist(SB1C$M2.KickingDuration) #Not normally distributed
#Try a log10 transformation
hist(log10(SB1C$M2.KickingDuration)) #Normally distributed

#Determine random factors for model

KD2.lmer1=lmer(log10(M2.KickingDuration) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
summary(KD2.lmer1)
#Remove +(1 | FemaleGeneration:FemaleMatriline)
KD2.lmer2=lmer(log10(M2.KickingDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(KD2.lmer1,KD2.lmer2) #Not significantly different
summary(KD2.lmer2)
#Remove +(1 | M1.MaleGeneration:M1.MaleMatriline)
KD2.lmer3=lmer(log10(M2.KickingDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(KD2.lmer3,KD2.lmer2) #Not significantly different
summary(KD2.lmer3)
#Remove +(1 | M2.MaleGeneration:M2.MaleMatriline)
KD2.lmer4=lmer(log10(M2.KickingDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(KD2.lmer3,KD2.lmer4) #Not significantly different
summary(KD2.lmer4)
#Remove +(1 | M2.MaleMatriline)
KD2.lmer5=lmer(log10(M2.KickingDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | Group),data= SB1C)
anova(KD2.lmer5,KD2.lmer4) #Not significantly different
summary(KD2.lmer5)
#Remove +(1 | M1.MaleMatriline)
KD2.lmer6=lmer(log10(M2.KickingDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M2.MaleGeneration)+(1 | Group),data= SB1C)
anova(KD2.lmer5,KD2.lmer6) #Not significantly different
summary(KD2.lmer6)
#Remove +(1 | FemaleMatriline)
KD2.lmer7=lmer(log10(M2.KickingDuration) ~  +(1 | FemaleGeneration)+(1 | M1.MaleGeneration)+(1 | M2.MaleGeneration)+(1 | Group),data= SB1C)
anova(KD2.lmer7,KD2.lmer6) #Not significantly different
summary(KD2.lmer7)
#Remove +(1 | Group)
KD2.lmer8=lmer(log10(M2.KickingDuration) ~  +(1 | FemaleGeneration)+(1 | M1.MaleGeneration)+(1 | M2.MaleGeneration),data= SB1C)
anova(KD2.lmer7,KD2.lmer8) #Not significantly different
summary(KD2.lmer8)
#Remove +(1 | M1.MaleGeneration)
KD2.lmer9=lmer(log10(M2.KickingDuration) ~  +(1 | FemaleGeneration)+(1 | M2.MaleGeneration),data= SB1C)
anova(KD2.lmer9,KD2.lmer8) #Not significantly different
summary(KD2.lmer9)
#Remove +(1 | FemaleGeneration)
KD2.lmer10=lmer(log10(M2.KickingDuration) ~  +(1 | M2.MaleGeneration),data= SB1C)
anova(KD2.lmer9,KD2.lmer10) #Not significantly different
summary(KD2.lmer10)
#Try switching out M2.MaleGeneration for FemaleGeneration to see if they differ; if they don't, then revert to a LM instead of a LMER
KD2.lmer11=lmer(log10(M2.KickingDuration) ~  +(1 | FemaleGeneration),data= SB1C)
anova(KD2.lmer11,KD2.lmer10) #Not significantly different

#Thus, revert to a lm

#Run bivariate analyses

#Predictors to consider: 
KD2.lm10=lm(log10(M2.KickingDuration) ~  M2.MaleWeight1,data= SB1C)
summary(KD2.lm10)    

#MatingOrder						0.452
#as.factor(IntermatingInterval)  	0.574, 0.307  
#SterilizationOrder       			0.225
#M1.FemaleAge, 0.126   
#M1.FemaleWeight1, 0.116    
#M1.FemaleWeightGain   				0.802   
#M1.MaleAge          				0.574
#M1.MaleWeight1, 0.0159 *
#M1.MaleEjacSize.Raw, 0.125   
#M1.CopulationLatency, 0.0406 *
#M1.KickingLatency, 0.0451 *
#M1.KickingDuration					0.216
#M1.CopulationDuration, 0.00834 **
#M2.FemaleWeightGain         		0.962
#M2.MaleAge          				0.984
#M2.MaleWeight1, 0.033343 *

#Predictors to consider: 
#MatingOrder + M1.FemaleAge + M1.FemaleWeight1 + M1.MaleWeight1 + M1.MaleEjacSize.Raw + M1.CopulationLatency + M1.KickingLatency + M1.CopulationDuration + M2.MaleWeight1

#Check for collinearities for these variables: 
lm1=lm(M1.MaleEjacSize.Raw ~ MatingOrder,data=SB1C)
summary(lm1)

#The following variables are collinear:
#M1.CopulationDuration~M1.FemaleAge (keep M1.CopulationDuration) 
#M1.CopulationDuration~M1.CopulationLatency (keep M1.CopulationDuration)
#M1.CopulationDuration ~ M1.KickingLatency (Keep M1.CopulationDuration)
#M2.MaleWeight1 ~ M1.MaleWeight1 (Keep M1.MaleWeight1)
#M1.MaleWeight1 ~ MatingOrder (Keep MatingOrder)


#Final model selection (AIC criteria <2)
KD2.lm10=lm(log10(M2.KickingDuration) ~  MatingOrder + M1.FemaleWeight1 + M1.MaleEjacSize.Raw + M1.CopulationDuration,data= SB1C)
summary(KD2.lm10)
#remove M1.MaleEjacSize.Raw
KD2.lm11=lm(log10(M2.KickingDuration) ~  MatingOrder + M1.FemaleWeight1 + M1.CopulationDuration,data= SB1C)
anova(KD2.lm10,KD2.lm11) #Not significantly different
summary(KD2.lm11)
#remove M1.FemaleWeight1
KD2.lm12=lm(log10(M2.KickingDuration) ~  MatingOrder + M1.CopulationDuration,data= SB1C)
anova(KD2.lm12,KD2.lm11) #Not significantly different
summary(KD2.lm12) 

#Check model residuals
plot(KD2.lm12) #homogeneous

#Mean +/- SE Kicking Duration for males in second mating role
#PM:
mean(SB1C$M2.KickingDuration[SB1C$MatingOrder=="PM"],na.rm=TRUE)
#Alternative calculation for SE: 
l=29
sd=sd(SB1C$M2.KickingDuration[SB1C$MatingOrder=="PM"],na.rm=TRUE)
se=(sd/sqrt(l))
se
#MP: 
mean(SB1C$M2.KickingDuration[SB1C$MatingOrder=="MP"],na.rm=TRUE)
l=33
sd=sd(SB1C$M2.KickingDuration[SB1C$MatingOrder=="MP"],na.rm=TRUE)
se=(sd/sqrt(l))
se


######### Copulation Duration

#Check data distribution
hist(SB1C$M2.CopulationDuration) #Not normally distributed
#Try a log transformation
hist(log2(SB1C$M2.CopulationDuration)) #Normally distributed

#Determine random factors for model
CD2.lmer1=lmer(log2(M2.CopulationDuration) ~ +(1 | FemaleGeneration:FemaleMatriline) +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
summary(CD2.lmer1)
#Remove +(1 | FemaleGeneration:FemaleMatriline)
CD2.lmer2=lmer(log2(M2.CopulationDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration:M2.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(CD2.lmer1,CD2.lmer2)  #Not significantly different
summary(CD2.lmer2)
#Remove +(1 | M2.MaleGeneration:M2.MaleMatriline)
CD2.lmer3=lmer(log2(M2.CopulationDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration:M1.MaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(CD2.lmer3,CD2.lmer2)  #Not significantly different
summary(CD2.lmer3)
#Remove +(1 | M1.MaleGeneration:M1.MaleMatriline)
CD2.lmer4=lmer(log2(M2.CopulationDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleMatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(CD2.lmer3,CD2.lmer4)  #Not significantly different
summary(CD2.lmer4)
#Remove +(1 | M1.MaleMatriline)
CD2.lmer5=lmer(log2(M2.CopulationDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline)+(1 | Group),data= SB1C)
anova(CD2.lmer5,CD2.lmer4)  #Not significantly different
summary(CD2.lmer5)
#Remove +(1 | Group)
CD2.lmer6=lmer(log2(M2.CopulationDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M2.MaleGeneration)+(1 | M2.MaleMatriline),data= SB1C)
anova(CD2.lmer5,CD2.lmer6)  #Not significantly different
summary(CD2.lmer6)
#Remove +(1 | M2.MaleGeneration)
CD2.lmer7=lmer(log2(M2.CopulationDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M1.MaleGeneration)+(1 | M2.MaleMatriline),data= SB1C)
anova(CD2.lmer7,CD2.lmer6)  #Not significantly different
summary(CD2.lmer7)
#Remove +(1 | M1.MaleGeneration)
CD2.lmer8=lmer(log2(M2.CopulationDuration) ~  +(1 | FemaleGeneration)+(1 | FemaleMatriline)+(1 | M2.MaleMatriline),data=SB1C)
anova(CD2.lmer7,CD2.lmer8)  #Not significantly different
summary(CD2.lmer8)
#Remove +(1 | FemaleGeneration)
CD2.lmer9=lmer(log2(M2.CopulationDuration) ~  +(1 | FemaleMatriline)+(1 | M2.MaleMatriline),data=SB1C)
anova(CD2.lmer9,CD2.lmer8)  #Not significantly different
summary(CD2.lmer9)
#Remove +(1 | FemaleMatriline)
CD2.lmer10=lmer(log2(M2.CopulationDuration) ~  +(1 | M2.MaleMatriline),data=SB1C)
anova(CD2.lmer9,CD2.lmer10)  #Not significantly different
summary(CD2.lmer10)
#Try switching out M2.MaleMatriline for FemaleMatriline to see if they differ; if they don't, then revert to a LM instead of a LMER 
CD2.lmer11=lmer(log2(M2.CopulationDuration) ~  +(1 | FemaleMatriline),data=SB1C)
anova(CD2.lmer11,CD2.lmer10)  #Yes, significantly different, so keep lmer with M2.MaleMatriline 

#Run bivariate analyses

CD2.lmer11=lmer(log2(M2.CopulationDuration) ~ M2.MaleWeight1 + (1 | M2.MaleMatriline),data=SB1C)
summary(CD2.lmer11)

#MatingOrder, 0.199
#as.factor(IntermatingInterval), 0.0395 *, 8.01e-05 ***
#SterilizationOrder, 0.101
#M1.FemaleAge           				0.731
#M1.FemaleWeight1         				0.213
#M1.FemaleWeightGain      				0.719
#M1.MaleAge            					0.821
#M1.MaleWeight1, 0.0843
#M1.MaleEjacSize.Raw, 0.0122 *
#M1.CopulationLatency    				0.393
#M1.KickingLatency  					0.548
#M1.KickingDuration   					0.733
#M1.CopulationDuration  	 			0.49
#M2.FemaleWeightGain        			0.985
#M2.MaleAge,0.0178 *
#M2.MaleWeight1           				0.941

#Predictors to consider: 
#MatingOrder + as.factor(IntermatingInterval) + SterilizationOrder + M1.MaleWeight1 + M1.MaleEjacSize.Raw + M2.MaleAge 

#Check for collinearities for these variables: 
lm1=lm(M2.MaleAge ~ as.factor(IntermatingInterval),data=SB1C)
summary(lm1)

#The following variables are collinear: 
#M1.MaleEjacSize.Raw ~ as.factor(IntermatingInterval) (keep as.factor(IntermatingInterval)) 
#M2.MaleAge~as.factor(IntermatingInterval) (keep as.factor(IntermatingInterval)) 
#M1.MaleWeight1~ MatingOrder (keep MatingOrder)

#MatingOrder + as.factor(IntermatingInterval) + SterilizationOrder  

#Final model selection (AIC criteria <2)
CD2.lmer11=lmer(log2(M2.CopulationDuration) ~ MatingOrder + as.factor(IntermatingInterval) + SterilizationOrder+ (1 | M2.MaleMatriline),data=SB1C)
summary(CD2.lmer11)
#remove SterilizationOrder
CD2.lmer12=lmer(log2(M2.CopulationDuration) ~ MatingOrder + as.factor(IntermatingInterval) + (1 | M2.MaleMatriline),data=SB1C) 
anova(CD2.lmer11,CD2.lmer12) #Not significantly different
summary(CD2.lmer12)

#Check model residuals
plot(CD2.lmer12) #homogeneous

#Mean +/- SE Copulation Duration for males in second mating role
#PM:
mean(SB1C$M2.CopulationDuration[SB1C$MatingOrder=="PM"],na.rm=TRUE)
#Alternative calculation for SE:  
l=28
sd=sd(SB1C$M2.CopulationDuration[SB1C$MatingOrder=="PM"],na.rm=TRUE)
se=(sd/sqrt(l))
se
#MP:
mean(SB1C$M2.CopulationDuration[SB1C$MatingOrder=="MP"],na.rm=TRUE) 
l=34
sd=sd(SB1C$M2.CopulationDuration[SB1C$MatingOrder=="MP"],na.rm=TRUE)
se=(sd/sqrt(l))
se



####FIGURE 3

w<-0.37   #   PM  Mating 1
wse<-0.03
x<-0.32   #   MP  Mating 1
xse<-0.02
y<-0.31   #   PM Mating 2
yse<-0.03
z<-0.23   #   MP Mating 2
zse<-0.02

Ejac1<-c(w,x)
Ejac2<-c(y,z)
Ejacs<-cbind(Ejac1, Ejac2)
bp<-barplot(Ejacs,beside=TRUE,col=c("white","gainsboro","white","gainsboro"),ylim=c(0,0.5),xaxt="n",cex.lab=1.25,cex.main=1.25,cex.axis=1.25,las=1)
title(ylab=expression(bold("Ejaculate Size  (mg)")),cex.lab=1.25,cex.main=1.25,cex.axis=1.25)
#axis(2,at=seq(0,1,0.2))
legend("topright", c("Monandrous","Polyandrous"),title="Maternal Treatment",lty=c(1,1),lwd=c(2.5,2.5),horiz=TRUE,col=c("white","gainsboro"),bty="o",box.lty=1, box.col=par(2),cex=0.8)
#Drawing standard error bars
ejac1se<-c(wse,xse)
ejac2se<-c(yse,zse)
se<-cbind(ejac1se,ejac2se)
arrows(bp, Ejacs,bp, Ejacs +se,lwd=1.5,angle=90,length=0.1)
box()


# Paired t-tests of ejaculate sizes for males 1 and 2 based on mating order (WITHIN TREATMENTS)

#Subset of data for MP matings
SBmp =SB1C[which(SB1C$MatingOrder=="MP"),]

#paired t-test
t.test(SBmp$M1.MaleEjacSize.Raw,SBmp$M2.MaleEjacSize.Raw, paired=TRUE). #significantly different
mean(SBmp$M1.MaleEjacSize.Raw)
mean(SBmp$M2.MaleEjacSize.Raw,na.rm=TRUE)

#Subset of data for PM matings
SBpm =SB1C[which(SB1C$MatingOrder=="PM"),]

#paired t-test
t.test(SBpm$M1.MaleEjacSize.Raw,SBpm$M2.MaleEjacSize.Raw, paired=TRUE). #not significantly different
mean(SBpm$M1.MaleEjacSize.Raw)
mean(SBpm$M2.MaleEjacSize.Raw,na.rm=TRUE)




##################################  Post-Hoc Analyses: MALE WEIGHTS  #################################### 

#Check data distribution
hist(SB1C.m$MaleWeight) #Normally distributed

#Determine random factors for model

#Determine random factors
weight.lmer1=lmer(MaleWeight ~ +(1 | MaleGeneration: MaleMatriline)+(1 | MaleGeneration)+(1 | MaleMatriline)+(1 | MalePatriline)+(1 | MaleGeneration: MalePatriline)+(1 | Group),data= SB1C.m)
summary(weight.lmer1)
#Remove +(1 | MalePatriline)
weight.lmer2=lmer(MaleWeight ~ +(1 | MaleGeneration: MaleMatriline)+(1 | MaleGeneration)+(1 | MaleMatriline)+(1 | MaleGeneration: MalePatriline)+(1 | Group),data= SB1C.m)
anova(weight.lmer1, weight.lmer2)  #Not significantly different
summary(weight.lmer2)
#Remove +(1 | MaleMatriline)
weight.lmer3=lmer(MaleWeight ~ +(1 | MaleGeneration: MaleMatriline)+(1 | MaleGeneration)+(1 | MaleGeneration: MalePatriline)+(1 | Group),data= SB1C.m)
anova(weight.lmer3, weight.lmer2)  #Not significantly different
summary(weight.lmer3)
#Remove +(1 | MaleGeneration)
weight.lmer4=lmer(MaleWeight ~ +(1 | MaleGeneration: MaleMatriline)+(1 | MaleGeneration: MalePatriline)+(1 | Group),data= SB1C.m)
anova(weight.lmer3, weight.lmer4)  #Not significantly different
summary(weight.lmer4)
#Remove +(1 | MaleGeneration: MalePatriline)
weight.lmer5=lmer(MaleWeight ~ +(1 | MaleGeneration: MaleMatriline)+(1 | Group),data= SB1C.m)
anova(weight.lmer5, weight.lmer4)  #Not significantly different
summary(weight.lmer5)
#Remove +(1 | MaleGeneration: MaleMatriline)
weight.lmer6=lmer(MaleWeight ~ +(1 | Group),data= SB1C.m)
anova(weight.lmer5, weight.lmer6)  #Not significantly different
summary(weight.lmer6)

#Check normality of residuals
qqnorm(resid(weight.lmer6))
qqline(resid(weight.lmer6)) #normal

#Run bivariate analyses
weight.lmer6=lmer(MaleWeight ~ as.factor(IntermatingInterval) +(1 | Group),data= SB1C.m)
summary(weight.lmer6)

#MaternalTreatment, 0.0345 * 
#MaternalTreatmentNumber, 			0.671
#as.factor(IntermatingInterval)		0.720, 0.662
#SterilizationOrder, 0.168
#MaleAge,3.15e-06 ***
#MaleDevelopmentTime				0.80613
#MaleSterile						0.942

#Predictors to consider: 
#MaternalTreatment + SterilizationOrder + MaleAge

#Check for collinearities for these variables: 
lm1=lm(MaleAge ~ SterilizationOrder,data=SB1C.m)
summary(lm1)

#The following variables are collinear: NONE

#Final model selection (AIC criteria <2)
weight.lmer6=lmer(MaleWeight ~ MaternalTreatment + SterilizationOrder + MaleAge +(1 | Group),data= SB1C.m)
summary(weight.lmer6)
#Remove + SterilizationOrder
weight.lmer7=lmer(MaleWeight ~ MaternalTreatment  + MaleAge +(1 | Group),data= SB1C.m)
anova(weight.lmer6,weight.lmer7) #Not significantly different
summary(weight.lmer7)

#Check model residuals
plot(weight.lmer7) #homogeneous

#Mean +/- SE XX for males in second mating role
mean(SB1C.m$MaleWeight[SB1C.m$MaternalTreatment =="Monandrous"])
#Alternative calculation for SE:  
l=70
sd=sd(SB1C.m$MaleWeight[SB1C.m$MaternalTreatment =="Monandrous"],na.rm=TRUE)
se=(sd/sqrt(l))
se
se = function(x) sd(x)/sqrt(length(x))
mean(SB1C.m$MaleWeight[SB1C.m$MaternalTreatment =="Polyandrous"])
se(SB1C.m$MaleWeight[SB1C.m$MaternalTreatment =="Polyandrous"])

#Least squares means adjusted values: 
lsmeans(weight.lmer7, pairwise ~ MaternalTreatment,type="response")

###FIGURE 3
t <- ggplot(SB1C.m, aes(MaternalTreatment, MaleWeight))
t + geom_boxplot(fill=c("white","gainsboro")) + geom_jitter(width = 0.2,size=2) + theme_classic() 



######################################   Post-Hoc Analyses (Discussion)   ###################################

#A closer look at the mothers from Part B of the experiment

#Do monandrous mothers weigh more?  NO 
boxplot(FemalePreMatingWeight ~ FemaleTreatment, xlab="Female Treatment", ylab="Female Weight Before First Mating (mg))", data=SB1B)
qqnorm(SB1B$FemalePreMatingWeight)
qqline(SB1B$FemalePreMatingWeight) #normal
var.test(SB1B$FemalePreMatingWeight ~SB1B$FemaleTreatment) #equal
t.test(SB1B$FemalePreMatingWeight ~SB1B$FemaleTreatment)

#Are monandrous mothers younger?  NO 
boxplot(FemaleAgeatFirstMating ~ FemaleTreatment, xlab="Female Treatment",ylab="Female Weight Before First Mating (mg))",data=SB1B)
qqnorm(SB1B$FemaleAgeatFirstMating)
qqline(SB1B$FemaleAgeatFirstMating) #normal
var.test(SB1B$FemaleAgeatFirstMating ~SB1B$FemaleTreatment) #equal
t.test(SB1B$FemaleAgeatFirstMating ~SB1B$FemaleTreatment)

#Did mothers lay different numbers of eggs?   NO

#TotalFecundity includes eggs laid during the inter-mating interval
boxplot(TotalFecundity ~ FemaleTreatment, xlab="Female Treatment", ylab="Number of Eggs Laid in Total (Fecundity)",ylim=c(0,90), data=SB1B)
qqnorm(SB1B$TotalFecundity)
qqline(SB1B$TotalFecundity) #normal?
var.test(SB1B$TotalFecundity ~SB1B$FemaleTreatment) #equal
t.test(SB1B$TotalFecundity ~SB1B$FemaleTreatment)

#Fecundity does not include eggs laid during the inter-mating interval
boxplot(Fecundity ~ FemaleTreatment, xlab="Female Treatment", ylab="Number of Eggs Laid After Second Mating (Fecundity)",ylim=c(0,90), data=SB1B)
qqnorm(SB1B$Fecundity)
qqline(SB1B$Fecundity) #normal
var.test(SB1B$Fecundity~SB1B$FemaleTreatment) #equal
t.test(SB1B$Fecundity~SB1B$FemaleTreatment)

#Do monandrous mothers have shorter lifespans?  #NO
boxplot(LifeSpan ~ FemaleTreatment, xlab="Female Treatment", ylab="Number of Eggs Laid in Total (Fecundity)", data=SB1B)
qqnorm(SB1B$LifeSpan)
qqline(SB1B$LifeSpan) #normal?
var.test(SB1B$LifeSpan ~SB1B$FemaleTreatment) #equal
t.test(SB1B$LifeSpan ~SB1B$FemaleTreatment)



#A closer look at focal females and focal males (i.e., sons) from Part C of the experiment


#Did focal females lay different raw numbers of inter-mating interval eggs?   NO
boxplot(PrimarySeedEggs.1~ M1.MomTreatment,data=SB1C)
qqnorm(SB1C$PrimarySeedEggs.1)
qqline(SB1C$PrimarySeedEggs.1) #Normal
var.test(SB1C$PrimarySeedEggs.1 ~SB1C$M1.MomTreatment) #equal variances
t.test(SB1C$PrimarySeedEggs ~SB1C$M1.MomTreatment)

#What about as a proportion? 
boxplot(PrimarySeedEggs.ProportionFecundity ~ M1.MomTreatment,data=SB1C)
qqnorm(SB1C$PrimarySeedEggs.ProportionFecundity)
qqline(SB1C$PrimarySeedEggs.ProportionFecundity) #not Normal
#try log transformation
logProp<-log(SB1C$PrimarySeedEggs.ProportionFecundity)
qqnorm(logProp)
qqline(logProp)
#try arcsin trans 
arcProp<-asin(sqrt(SB1C$PrimarySeedEggs.ProportionFecundity))
qqnorm(arcProp)
qqline(arcProp)  #better
var.test(arcProp~SB1C$M1.MomTreatment) #equal variances
t.test(arcProp~SB1C$M1.MomTreatment)
t.test(SB1C$PrimarySeedEggs.ProportionFecundity ~ SB1C$M1.MomTreatment)

#Did focal females lay different raw numbers of eggs after their first mating (Fecundity)?   NO
boxplot(Fecundity ~ MatingOrder,data=SB1C)
qqnorm(SB1C$Fecundity)
qqline(SB1C$Fecundity) #Normal
var.test(SB1C$Fecundity ~SB1C$MatingOrder) #equal variances
t.test(SB1C$Fecundity ~SB1C$MatingOrder)


#Are monandrous sons larger? NO, they are sig. smaller
#Looking only within first mating only since males differed in ages for second mating (and thus also weights depending on inter-mating interval) 
boxplot2(M1.MaleWeight1~ M1.MomTreatment,data=SB1C,xlab="Mom's Treatment",ylab="Male Pre-Mating Weight (mg)",ylim=c(2,6))
text(1,2,"(39)")
text(2,2,"(33)")
qqnorm(SB1C$M1.MaleWeight1)
qqline(SB1C$M1.MaleWeight1) #normal
var.test(SB1C$M1.MaleWeight1~SB1C$M1.MomTreatment)
t.test(SB1C$M1.MaleWeight1~SB1C$M1.MomTreatment)
#do differ in age? NO
boxplot(M1.MaleAge ~ M1.MomTreatment,data= SB1C)
qqnorm(SB1C$M1.MaleAge)
qqline(SB1C$M1.MaleAge) #normal?
var.test(SB1C$M1.MaleAge ~SB1C$M1.MomTreatment) #equal
t.test(SB1C$M1.MaleAge ~SB1C$M1.MomTreatment)


#Did first mates that were polyandrous transfer smaller ejacs?  NO
boxplot(M1.MaleEjacSize.Raw ~ M1.MomTreatment,data= SB1C)
qqnorm(SB1C$M1.MaleEjacSize.Raw)
qqline(SB1C$M1.MaleEjacSize.Raw) #normal
var.test(SB1C$M1.MaleEjacSize.Raw ~SB1C$M1.MomTreatment) #equal
t.test(SB1C$M1.MaleEjacSize.Raw ~SB1C$M1.MomTreatment)

#Did second mates that were monandrous transfer larger ejacs? Yes
boxplot(M2.MaleEjacSize.Raw ~ M2.MomTreatment,data= SB1C, xlab="Mom's Treatment",ylab="Second Male Ejaculate Size (mg)",ylim=c(-0.1,0.8),col="gainsboro")
text(1,-0.1,"(32)")
text(2,-0.1,"(38)")
qqnorm(SB1C$M2.MaleEjacSize.Raw)
qqline(SB1C$M2.MaleEjacSize.Raw) #normal
var.test(SB1C$M2.MaleEjacSize.Raw ~SB1C$M2.MomTreatment) #equal
t.test(SB1C$M2.MaleEjacSize.Raw ~SB1C$M2.MomTreatment)


#Did sons differ in their development time?  #No
boxplot(MaleDEVTime ~ MaleMom,data= SB1C.m)
qqnorm(log10(SB1C.m$MaleDEVTime))
qqline(log10(SB1C.m$MaleDEVTime)) # not normal
var.test(SB1C.m$MaleDEVTime ~ SB1C.m$MaleMom)  #equal
t.test(SB1C.m$MaleDEVTime ~ SB1C.m$MaleMom)  - NS different P = 0.75


#FEMALEWEIGHTGAIN 

#Did focal females gain less weight in first mating for polyandrous males? NO
boxplot(M1.FemaleWeightGain ~ M1.MomTreatment,data= SB1C)
qqnorm(SB1C$M1.FemaleWeightGain)
qqline(SB1C$M1.FemaleWeightGain) #normal
var.test(SB1C$M1.FemaleWeightGain~SB1C$M1.MomTreatment) #equal
t.test(SB1C$M1.FemaleWeightGain~SB1C$M1.MomTreatment)

#Did focal females gain more weight in second mating for monandrous males? NO
boxplot(M2.FemaleWeightGain ~ M2.MomTreatment,data= SB1C)
qqnorm(SB1C$M2.FemaleWeightGain)
qqline(SB1C$M2.FemaleWeightGain) #normal
var.test(SB1C$M2.FemaleWeightGain~SB1C$M2.MomTreatment) #equal
t.test(SB1C$M2.FemaleWeightGain~SB1C$M2.MomTreatment)


#Boxplot 
boxplot2(M2.FemaleWeightGain ~ MatingOrder, xlab="Treatment", ylab="Proportion of Eggs Sired by Second Male (P2)",col=c("gainsboro","white"),ylim=c(0,1.1),data= SB1C)

#Determine random effects for model, unnested and nested
Fem.lmer=lmer(M2.FemaleWeightGain ~ +(1 | FemaleGeneration)+(1 | FemaleGeneration: FemaleMatriline)+(1 | FemaleGeneration: FemaleMatriline2)+(1 | FemaleGeneration: FemPatriline)+(1 | M1.MaleGeneration)+(1 | M1.MaleGeneration: M1.MaleMatriline)+(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M1.MaleGeneration: M1.MalePatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | M2.MaleGeneration: M2.MaleMatriline2)+(1 | M2.MaleGeneration: M2.MalePatriline),data= SB1C)
summary(Fem.lmer) 
Fem.lmer1=lmer(M2.FemaleWeightGain ~ +(1 | FemaleGeneration)+(1 | FemaleGeneration: FemaleMatriline)+(1 | FemaleGeneration: FemaleMatriline2)+(1 | M1.MaleGeneration)+(1 | M1.MaleGeneration: M1.MaleMatriline)+(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M1.MaleGeneration: M1.MalePatriline)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | M2.MaleGeneration: M2.MaleMatriline2)+(1 | M2.MaleGeneration: M2.MalePatriline),data= SB1C)
anova(Fem.lmer,Fem.lmer1) #NS
summary(Fem.lmer1) 
Fem.lmer2=lmer(M2.FemaleWeightGain ~ +(1 | FemaleGeneration)+(1 | FemaleGeneration: FemaleMatriline)+(1 | FemaleGeneration: FemaleMatriline2)+(1 | M1.MaleGeneration)+(1 | M1.MaleGeneration: M1.MaleMatriline)+(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | M2.MaleGeneration: M2.MaleMatriline2)+(1 | M2.MaleGeneration: M2.MalePatriline),data= SB1C)
anova(Fem.lmer2,Fem.lmer1) #NS
summary(Fem.lmer2) 
Fem.lmer3=lmer(M2.FemaleWeightGain ~ +(1 | FemaleGeneration)+(1 | FemaleGeneration: FemaleMatriline)+(1 | FemaleGeneration: FemaleMatriline2)+(1 | M1.MaleGeneration)+(1 | M1.MaleGeneration: M1.MaleMatriline)+(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | M2.MaleGeneration: M2.MaleMatriline2),data= SB1C)
anova(Fem.lmer2,Fem.lmer3) #NS
summary(Fem.lmer3)
Fem.lmer4=lmer(M2.FemaleWeightGain ~ +(1 | FemaleGeneration)+(1 | FemaleGeneration: FemaleMatriline2)+(1 | M1.MaleGeneration)+(1 | M1.MaleGeneration: M1.MaleMatriline)+(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | M2.MaleGeneration: M2.MaleMatriline2),data= SB1C)
anova(Fem.lmer4,Fem.lmer3) #NS
summary(Fem.lmer4)
Fem.lmer5=lmer(M2.FemaleWeightGain ~ +(1 | FemaleGeneration)+(1 | FemaleGeneration: FemaleMatriline2)+(1 | M1.MaleGeneration)+(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline)+(1 | M2.MaleGeneration: M2.MaleMatriline2),data= SB1C)
anova(Fem.lmer4,Fem.lmer5) #NS
summary(Fem.lmer5)
Fem.lmer6=lmer(M2.FemaleWeightGain ~ +(1 | FemaleGeneration)+(1 | FemaleGeneration: FemaleMatriline2)+(1 | M1.MaleGeneration)+(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline2),data= SB1C)
anova(Fem.lmer6,Fem.lmer5) #NS
summary(Fem.lmer6)
Fem.lmer7=lmer(M2.FemaleWeightGain ~ +(1 | FemaleGeneration: FemaleMatriline2)+(1 | M1.MaleGeneration)+(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline2),data= SB1C)
anova(Fem.lmer6,Fem.lmer7) #NS
summary(Fem.lmer7)
Fem.lmer8=lmer(M2.FemaleWeightGain ~ +(1 | FemaleGeneration: FemaleMatriline2)+(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M2.MaleGeneration)+(1 | M2.MaleGeneration: M2.MaleMatriline2),data= SB1C)
anova(Fem.lmer8,Fem.lmer7) #NS
summary(Fem.lmer8)
Fem.lmer9=lmer(M2.FemaleWeightGain ~ +(1 | FemaleGeneration: FemaleMatriline2)+(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M2.MaleGeneration: M2.MaleMatriline2),data= SB1C)
anova(Fem.lmer8,Fem.lmer9) #NS
summary(Fem.lmer9)
Fem.lmer10=lmer(M2.FemaleWeightGain ~ +(1 | M1.MaleGeneration: M1.MaleMatriline2)+(1 | M2.MaleGeneration: M2.MaleMatriline2),data= SB1C)
anova(Fem.lmer10,Fem.lmer9) #NS
summary(Fem.lmer10)
Fem.lmer11=lmer(M2.FemaleWeightGain ~ +(1 | M2.MaleGeneration: M2.MaleMatriline2),data= SB1C)
anova(Fem.lmer10,Fem.lmer11) #NS
summary(Fem.lmer11)

#NO LMER NEEDED, go with a LM

SB1C$M2.EjaculateSizeDiff2 <- scale(SB1C$M2.EjaculateSizeDiff, center = TRUE, scale = TRUE)
SB1C$M2.CopulationDurationDiff2 <- scale(SB1C$M2.CopulationDurationDiff, center = TRUE, scale = TRUE)
SB1C$M2.CopulationLatencyDiff2 <- scale(SB1C$M2.CopulationLatencyDiff, center = TRUE, scale = TRUE)
SB1C$M2.KickingLatencyDiff2 <- scale(SB1C$M2.KickingLatencyDiff, center = TRUE, scale = TRUE)
SB1C$M2KickingDurationDiff2 <- scale(SB1C$M2KickingDurationDiff, center = TRUE, scale = TRUE)

Fem.lm=lm(M2.FemaleWeightGain ~ MatingOrder,data= SB1C)
summary(Fem.lm)

#M1.FemaleAge           		0.372
#M1.FemaleWeight1    			0.950
#M1.FemaleWeightGain			0.606
#M2.MaleAgeDiff				0.783		   
#M1.MaleAge				0.553874
#M1.MaleWeight1				0.259
#M2.MaleWeightDiff			0.411
#M2.EjaculateSizeDiff 			0.16
#M2.KickingLatencyDiff		0.41	
#M2KickingDurationDiff		0.643  
#M2.CopulationDurationDiff				0.626
#as.factor(IntermatingInterval)    0.0389 * ,0.3151
#SterilizationOrder		 	0.184
#MatingOrder     			0.665

#check collinearities
lm1=lm(M2.EjaculateSizeDiff ~ as.factor(IntermatingInterval),data=SB1C)
summary(lm1)

Fem.lm1=lm(M2.FemaleWeightGain ~ MatingOrder+as.factor(IntermatingInterval)+ SterilizationOrder +M2.EjaculateSizeDiff,data= SB1C)
summary(Fem.lm1)
options(contrasts=c("contr.sum","contr.poly"))
#gives type 3 tests to test effects of everything in your model 
drop1(Fem.lm1,.~.,test="F")   #NO EFFECTS ON M2.FemaleWTGAIN



