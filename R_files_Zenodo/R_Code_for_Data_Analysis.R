
#Upload libraries
library(gtools)
library(ggplot2)
library(gplots)
library(lme4)
library(lmerTest)
library(lsmeans)
library(plyr)

#######  Upload data files

#Testes measurements and sperm counts for maniculatus (BW) and polionotus (PO) males
#Upload "Testes_Data_File.csv" file and name it "Testes"
#Testes
#names(Testes)

#Testes restricted to just BW
BW = Testes[which(Testes$Species=="BW"),]

#Testes restricted to just PO
PO = Testes[which(Testes$Species=="PO"),]

#Testes measurements and litter composition data for F2 hybrids
#Upload "Litter_Data_File.csv" file and name it "Litter"
#Litter
#names(Litter)



#######  Data analysis for litter composition study

#Species differences in Testes Measurements
#Mean +/- SE Testes Data

##### TABLE 1: 
#Mean +/-SE Body Mass (mg)
#calculate mean
Me<- tapply(Testes$Body.Mass..mg., Testes$Species, mean,na.rm=TRUE)
Me/1000
#Standard error function
se = function(x) sd(x)/sqrt(length(x))
se2 = tapply(Testes$Body.Mass..mg.,Testes$Species,se)
se2/1000

#Mean +/-SE Testis Mass (mg)
#calculate mean
Me<- tapply(Testes$Left.Testis.Mass..mg., Testes$Species, mean,na.rm=TRUE)
Me
#Standard error function
se = function(x) sd(x)/sqrt(length(x))
se2 = tapply(Testes$Left.Testis.Mass..mg., Testes$Species,se)
se2

#Mean +/-SE Testis Area (mm2)
#calculate mean
Me<- tapply(Testes$Left.Testis.Area..mm2., Testes$Species, mean,na.rm=TRUE)
Me
#Standard error function
se = function(x) sd(x)/sqrt(length(x))
se2 = tapply(Testes$Left.Testis.Area..mm2., Testes$Species,se)
se2
#alternative calculation for BW
l=57
sd=8.01
se=(sd/sqrt(l))
se

## Does testis mass differ between BW and PO? 

Weight.lm=lm(Combined.Testis.Mass..mg. ~ Species + Body.Mass..mg.,data=Testes)
summary(Weight.lm) 
# Check for heteroscedasticity
plot(resid(Weight.lm)~fitted(Weight.lm)) # Not heteroscedastic
# Check Normality 
qqnorm(residuals(Weight.lm),pch=16,cex=1.2)
qqline(residuals(Weight.lm),lwd=4,col=2) # Mostly normal
#Try log transformation of variables
Weight.lm2=lm(log(Combined.Testis.Mass..mg.) ~ Species + log(Body.Mass..mg.),data=Testes)
summary(Weight.lm2)
# Check for heteroscedasticity
plot(resid(Weight.lm2)~fitted(Weight.lm2)) # More heteroscedastic
# Check Normality 
qqnorm(residuals(Weight.lm2),pch=16,cex=1.2)
qqline(residuals(Weight.lm2),lwd=4,col=2) # Mostly normal


## Does testis area differ between BW and PO?  

Area.lm=lm(Combined.Testis.Area..mm2. ~ Species + Body.Mass..mg.,data=Testes)
summary(Area.lm)
# Check for heteroscedasticity
plot(resid(Area.lm)~fitted(Area.lm)) # Not heteroscedastic
# Check Normality 
qqnorm(residuals(Area.lm),pch=16,cex=1.2)
qqline(residuals(Area.lm),lwd=4,col=2) # Not quite normal
#Try log transformation of variables
Area.lm2=lm(log(Combined.Testis.Area..mm2.) ~ Species + log(Body.Mass..mg.),data=Testes)
summary(Area.lm2)
# Check for heteroscedasticity
plot(resid(Area.lm2)~fitted(Area.lm2)) # Somewhat more heteroscedastic
# Check Normality 
qqnorm(residuals(Area.lm2),pch=16,cex=1.2)
qqline(residuals(Area.lm2),lwd=4,col=2) # More normal


#Does testis area correlate with testis weight?
  
#BW
Size.lm1=lm(Combined.Testis.Area..mm2. ~ Combined.Testis.Mass..mg. + Age + Body.Mass..mg.,data=BW)
summary(Size.lm1)
#Remove Body.Mass..mg.
Size.lm2=lm(Combined.Testis.Area..mm2. ~ Combined.Testis.Mass..mg. + Age,data= BW)
anova(Size.lm1, Size.lm2) #Not significantly different
summary(Size.lm2)
#Remove Age
Size.lm3=lm(Combined.Testis.Area..mm2. ~ Combined.Testis.Mass..mg.,data= BW)
anova(Size.lm3, Size.lm2) #Error message, but result the same whether control for age or not
summary(Size.lm3)

#PO
Size.lm1=lm(Combined.Testis.Area..mm2. ~ Combined.Testis.Mass..mg. + Age + Body.Mass..mg.,data=PO)
summary(Size.lm1)
#Remove Body.Mass..mg.
Size.lm2=lm(Combined.Testis.Area..mm2. ~ Combined.Testis.Mass..mg. + Age,data=PO)
anova(Size.lm1, Size.lm2) #Not significantly different
summary(Size.lm2)
#Remove Age
Size.lm3=lm(Combined.Testis.Area..mm2. ~ Combined.Testis.Mass..mg.,data=PO)
anova(Size.lm3, Size.lm2) #Not significantly different
summary(Size.lm3)

#Both together
Size.lm1=lm(Combined.Testis.Area..mm2. ~ Combined.Testis.Mass..mg. + Age + Body.Mass..mg.,data=Testes)
summary(Size.lm1)
#Remove Body.Mass..mg.
Size.lm2=lm(Combined.Testis.Area..mm2. ~ Combined.Testis.Mass..mg. + Age,data=Testes)
anova(Size.lm1, Size.lm2) #Not significantly different
summary(Size.lm2)
#Remove Age
Size.lm3=lm(Combined.Testis.Area..mm2. ~ Combined.Testis.Mass..mg. ,data=Testes)
anova(Size.lm3, Size.lm2) #Error message, but result the same whether control for age or not
summary(Size.lm3)



# Does social composition of the rearing environment impact testis size? 

#Check normality
qqnorm(Litter$Mean.Testis.Area..mm2.)
qqline(Litter$Mean.Testis.Area..mm2.) #normal

#Run linear mixed model (LMM)

#Determine random factors 
#Random effects to consider:
#Grandparents
#Sire.ID
Testis.lmer1=lmer(Mean.Testis.Area..mm2.~ + (1 | Grandparents)+ (1 | Sire.ID)+ (1 | Grandparents: Sire.ID),data=Litter)
summary(Testis.lmer1) 
#Remove + (1 | Sire.ID)
Testis.lmer2=lmer(Mean.Testis.Area..mm2. ~ + (1 | Grandparents)+ (1 | Grandparents: Sire.ID),data=Litter)
anova(Testis.lmer1, Testis.lmer2) 	#Not significantly different
summary(Testis.lmer2)
#Remove + (1 | Grandparents: Sire.ID)
Testis.lmer3=lmer(Mean.Testis.Area..mm2. ~ + (1 | Grandparents),data=Litter)
anova(Testis.lmer3, Testis.lmer2)	#Not significantly different
summary(Testis.lmer3)

#Fixed effects to consider: 
#Male.Litter.Mates = number of brothers within litter
#Female.Litter.Mates = number of sisters within litter
#Proportion.Males.in.Litter = proportion of males to females within focal males' litter
#Litter.Order = which litter the male was born into (ranges from 0 to 21)
#Litter.Size = total number of offspring in the males' litter
#Age = age of male at processing
#Pairing.Duration..days. = duration of time that the male was paired with a female mate
#Reproductive.Success = reproductive success of focal male, 0 = male never produced a litter, 1 = male successfully produced a litter
#Weight..g. = male's body weight
#Mate.ID = female mate ID for pairing

#Run bivariate analysis:
Testis.lmer=lmer(Mean.Testis.Area..mm2. ~ Fertility + (1 | Grandparents),data=Litter)
summary(Testis.lmer) 

#No..Male.Littermates, 5.58e-06 ***
#No..Female.Littermates, 					0.30186
#Proportion.Males.in.Litter, 0.0374 *  
#Litter.Order, 							0.907770
#Litter.Size, 0.00222 ** 
#Age, 									0.698          ## Age had no effect on testis size
#Pairing.Duration..days., 				0.474182		   ## Pairing duration had no effect on testis size
#Reproductive.Success, 0.06152
#Weight..g., 							0.238

#Predictors to consider for the full model: 
No..Male.Littermates + Proportion.Males.in.Litter + Litter.Size + Reproductive.Success

#Check collinearities between these predictors:
Count=lm(No..Male.Littermates ~ Reproductive.Success, data= Testis)
summary(Count) 

#The following variables are collinear:
#No..Male.Littermates ~ Proportion.Males.in.Litter (Use former)
#No..Male.Littermates ~ Litter.Size (Use former)

#Predictors to use within the full model: 
No..Male.Littermates + Reproductive.Success

#Model selection
Testis.lmer1=lmer(Mean.Testis.Area..mm2. ~ No..Male.Littermates + Reproductive.Success + (1 | Grandparents),data=Litter)
summary(Testis.lmer1) 

#Check assumptions
# Check for heteroscedasticity
plot(resid(Testis.lmer1)~fitted(Testis.lmer1)) # Heteroscedastic
# Check Normality 
qqnorm(residuals(Testis.lmer1),pch=16,cex=1.2)
qqline(residuals(Testis.lmer1),lwd=4,col=2) # Normal


## Table 2: 
#(Intercept)  				26.9679     0.7650   5.2600  35.254 1.88e-07 ***
#Male.Litter.Mates         	1.4851     0.3087 295.5200   4.811 2.40e-06 ***
#Reproductive.Success      	1.8403     0.8093 296.9100   2.274   0.0237 * 


# Did age affect fertility?
Age.lmer=lmer(Fertility ~ Age + (1 | Grandparents),data=Litter)
summary(Age.lmer) 
# Grandparents does not impact response variable, converting to a LM: 
Age.lm=lm(Fertility ~ Age,data=Litter)
summary(Age.lm)


# Did pairing status affect fertility?
Pairing.lmer=lmer(Fertility ~ Pairing.Duration..days. + (1 | Grandparents),data=Litter)
summary(Pairing.lmer) 
# Grandparents does not impact response variable, converting to a LM: 
Pairing.lm=lm(Fertility ~ Pairing.Duration..days.,data=Litter)
summary(Pairing.lm)


# Does testis size associate with maternal identity?
Maternal.lmer=lmer(Mean.Testis.Area..mm2.~ Dam.ID + (1 | Grandparents),data=Litter)
summary(Maternal.lmer) 


# Does testis size associate with litter order?
LitterOrder.lmer=lmer(Mean.Testis.Area..mm2.~ Litter.Order+ (1 | Grandparents),data=Litter)
summary(LitterOrder.lmer) 




