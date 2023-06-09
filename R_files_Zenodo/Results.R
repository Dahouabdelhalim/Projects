
# Upload libraries

library(lme4)
library(lmerTest)
library(lsmeans)
library(RVAideMemoire) 
library(ape)
library(geiger)
library(caper)

# Upload data


# Upload 'Aggregation by Species.no outlier.csv' file - name it "Pero.Species"

# Upload 'Aggregation by Male.csv' - name it "Pero.Male2"

# Manually remove the IS outlier from the Peromyscus by male file: 
Pero.Male= Pero.Male2[which(Pero.Male2$ID !="8997"),]

# Make subsets of this data file by species: 
BW = Pero.Male[which(Pero.Male$Species=="BW"),] # This is Peromyscus maniculatus
EP = Pero.Male[which(Pero.Male$Species=="EP"),] # This is Peromyscus eremicus
GO = Pero.Male[which(Pero.Male$Species=="GO"),] # This is Peromyscus gossypinus
IS = Pero.Male[which(Pero.Male$Species=="IS"),] # This is Peromyscus californicus
LL = Pero.Male[which(Pero.Male$Species=="LL"),] # This is Peromyscus leucopus
PO = Pero.Male[which(Pero.Male$Species=="PO"),] # This is Peromyscus polionotus

IS2 = Pero.Male2[which(Pero.Male2$Species=="IS"),] # This is Peromyscus californicus and includes the outlier male (IS8997)

# Upload the Peromyscus phylogenetic tree using read.nexus() rather than read.table():
# cytb_cleaned_mcc.trees - name it "Pero.Tree" 
Pero.Tree<-read.nexus("PATHWAY/cytb_cleaned_mcc.trees")

## Determine if the tree is ultrametric
#is.ultrametric(Pero.tree) # True (it is)

# Combine phylogenetic tree and species data set 
Pero.PGLS <- comparative.data(phy = Pero.Tree, data = Pero.Species, names.col = tip.label, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE) #Ignore warning message

# Notes: 
# Low viscosity (i.e., control media), High viscosity (i.e., experimental media with 0.75% methylcellulose unless marked otherwise)



#############################################    Methods    ############################################# 

## Sample size per species
summary(Pero.Male2$Species)

## What is the mean and standard error for the post-harvest time (i.e. time between dissection and video recording)?
mean(Pero.Male$Post.Harvest.Time..min.) 	## 59.27612
#calculate standard error
se = function(x) sd(x)/sqrt(length(x))
se(Pero.Male$Post.Harvest.Time..min.) 	## 1.42431



#############################################      RESULTS       ############################################# 



##############    AGGREGATE SIZE     ################


### Do species vary in the sizes of sperm aggregates?

# Determine random factors
AggSize.lmer=lmer(Mean.Agg.Size ~ + (1 | Family),data= Pero.Male)
summary(AggSize.lmer)

#Run bivariate analysis of traits

# Rescale Variable
Pero.Male$Age2 <- scale(Pero.Male$Age, center = TRUE, scale = TRUE)
Pero.Male$Time <- scale(Pero.Male$Post.Harvest.Time..min., center = TRUE, scale = TRUE)

AggSize.lmer=lmer(Mean.Agg.Size ~ as.factor(Species) + (1 | Family),data= Pero.Male)
summary(AggSize.lmer)

#Age2,  0.04 *
#as.factor(Paired),  0.0252 *
#Time, 0.227
#Total.Sperm.Cells, 1.95e-08 ***
#Sperm.Cells.per.Video, 0.0134 *
#Number.of.CASA.Videos, 0.000381 ***
#as.factor(Species), all *** except PO 0.97

# Fixed factors to consider in final model
#Age2 + as.factor(Paired) + Total.Sperm.Cells + Sperm.Cells.per.Video + Number.of.CASA.Videos

# Check collinearities using a linear model
lm1=lm(Time ~ as.factor(Species),data= Pero.Male)
summary(lm1)

# The following factors are collinear: 
# Sperm.Cells.per.Video ~ Total.Sperm.Cells  	* use latter
# as.factor(Paired) ~ Total.Sperm.Cells 		* use latter
# Number.of.CASA.Videos ~ as.factor(Species) 	* use latter
# Age2  ~ as.factor(Species)					* use latter

# Final model selection
AggSize.lmer1=lmer(Mean.Agg.Size ~ Total.Sperm.Cells + as.factor(Species) + (1 | Family),data = Pero.Male)
summary(AggSize.lmer1)

# Family has a minimal effect on response variable, so converting to a linear model

#Run bivariate analysis of traits
AggSize.lm1=lm(Mean.Agg.Size ~ as.factor(Species), data= Pero.Male)
summary(AggSize.lm1)

#Age2,  0.00375 ** 
#as.factor(Paired),  0.00357 **
#Time, 0.0046 **
#Total.Sperm.Cells, 0.0087 ** 
#Sperm.Cells.per.Video,  0.057
#Number.of.CASA.Videos, 1.92e-10 ***
#as.factor(Species), all *** except PO 0.964

# Check collinearities using a linear model
lm1=lm(Total.Sperm.Cells ~ as.factor(Paired),data= Pero.Male)
summary(lm1)

# The following factors are collinear: 
# Age2  ~ as.factor(Species)					* use latter
# Time ~ as.factor(Species)						* use latter
# Sperm.Cells.per.Video ~ as.factor(Species)	* use latter
# Number.of.CASA.Videos ~ as.factor(Species) 	* use latter
# as.factor(Paired) ~ Total.Sperm.Cells 		* use former

# Final model selection
AggSize.lm1=lm(Mean.Agg.Size ~ as.factor(Paired) + as.factor(Species), data= Pero.Male)
summary(AggSize.lm1)
# Remove as.factor(Paired)
AggSize.lm2=lm(Mean.Agg.Size ~ as.factor(Species), data= Pero.Male)
anova(AggSize.lm1, AggSize.lm2) # models not fit to same dataset so cannot compare
summary(AggSize.lm2) # but result is no different, and model 1 drops 5 individuals, so using model 1

#Check model assumptions
# Check for heteroscedasticity
plot(resid(AggSize.lm2)~fitted(AggSize.lm2)) # Mostly heteroscedastic
# Check Normality 
qqnorm(residuals(AggSize.lm2),pch=16,cex=1.2)
qqline(residuals(AggSize.lm2),lwd=4,col=2) # Not normal

# Are the data normal? 
hist(Pero.Male$Mean.Agg.Size)
shapiro.test(Pero.Male$Mean.Agg.Size) # No, W = 0.83859, p-value = 9.12e-11
# Try a log transformation of data
hist(log(Pero.Male$Mean.Agg.Size))
shapiro.test(log(Pero.Male$Mean.Agg.Size)) # No, but more normal, W = 0.90353, p-value = 9.142e-08

# Rerun with transformed data
AggSize.lm2b=lm(log(Mean.Agg.Size) ~ as.factor(Species),data= Pero.Male)
summary(AggSize.lm2b)

#Check model assumptions
# Check for heteroscedasticity
plot(resid(AggSize.lm2b)~fitted(AggSize.lm2b)) # More heteroscedastic
# Check Normality 
qqnorm(residuals(AggSize.lm2b),pch=16,cex=1.2)
qqline(residuals(AggSize.lm2b),lwd=4,col=2) # Normal with a few outliers

# Final model is AggSize.lm2b


## control for total sperm cells
AggSize.lm3=lm(log(Mean.Agg.Size) ~ Total.Sperm.Cells + as.factor(Species), data= Pero.Male)
anova(AggSize.lm2b, AggSize.lm3) # models are significantly different, though model 3 has one less individual in the dataset
summary(AggSize.lm3) # But R value indicates a better fit

# Are the data normal? 
hist(Pero.Male$Total.Sperm.Cells)
shapiro.test(Pero.Male$Total.Sperm.Cells) # Yes, W = 0.99012, p-value = 0.4612

#Check model assumptions
# Check for heteroscedasticity
plot(resid(AggSize.lm3)~fitted(AggSize.lm3)) # Mostly heteroscedastic
# Check Normality 
qqnorm(residuals(AggSize.lm3),pch=16,cex=1.2)
qqline(residuals(AggSize.lm3),lwd=4,col=2) # Somewhat normal, with outlier tails

# Reporting values from AggSize.lm3


# FOR TABLE 1:    
#                       Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           1.3640475  0.0658074  20.728  < 2e-16 ***

#Total.Sperm.Cells     0.0008065  0.0001718   4.694 6.89e-06 ***
exp(0.0008065)/(1+exp(0.0008065))
exp(0.0008065-1.96* 0.0001718)/(1+exp(0.0008065-1.96* 0.0001718))
exp(0.0008065 +1.96* 0.0001718)/(1+exp(0.0008065 +1.96* 0.0001718)) 

#as.factor(Species)EP -0.5270055  0.0623208  -8.456 5.85e-14 ***
exp(-0.5270055)/(1+exp(-0.5270055))
exp(-0.5270055-1.96* 0.0623208)/(1+exp(-0.5270055-1.96* 0.0623208))
exp(-0.5270055 +1.96* 0.0623208)/(1+exp(-0.5270055 +1.96* 0.0623208))

#as.factor(Species)GO -0.7913487  0.0618754 -12.789  < 2e-16 ***
exp(-0.7913487)/(1+exp(-0.7913487))
exp(-0.7913487-1.96* 0.0618754)/(1+exp(-0.7913487-1.96* 0.0618754))
exp(-0.7913487 +1.96* 0.0618754)/(1+exp(-0.7913487 +1.96* 0.0618754))

#as.factor(Species)IS -0.5834052  0.0582571 -10.014  < 2e-16 ***
exp(-0.5834052)/(1+exp(-0.5834052))
exp(-0.5834052-1.96* 0.0582571)/(1+exp(-0.5834052-1.96* 0.0582571))
exp(-0.5834052 +1.96* 0.0582571)/(1+exp(-0.5834052 +1.96* 0.0582571))

#as.factor(Species)LL -0.7565332  0.0619109 -12.220  < 2e-16 ***
exp(-0.7565332)/(1+exp(-0.7565332))
exp(-0.7565332-1.96* 0.0619109)/(1+exp(-0.7565332-1.96* 0.0619109))
exp(-0.7565332 +1.96* 0.0619109)/(1+exp(-0.7565332 +1.96* 0.0619109))

#as.factor(Species)PO -0.0180014  0.0600892  -0.300    0.765 
exp(-0.0180014)/(1+exp(-0.0180014))
exp(-0.0180014-1.96* 0.0600892)/(1+exp(-0.0180014-1.96* 0.0600892))
exp(-0.0180014 +1.96* 0.0600892)/(1+exp(-0.0180014 +1.96* 0.0600892))

### Is there more variance across than within species?  YES

var(Pero.Male$Mean.Agg.Size,na.rm=TRUE) #1.961014

var(IS$Mean.Agg.Size,na.rm=TRUE) 	#0.6184152 (2.900569 with outlier)
var(EP$Mean.Agg.Size) 				#0.3059087
var(PO$Mean.Agg.Size,na.rm=TRUE) 	#2.150459    
var(BW$Mean.Agg.Size,na.rm=TRUE) 	#0.8228839
var(LL$Mean.Agg.Size,na.rm=TRUE) 	#0.121562
var(GO$Mean.Agg.Size) 			 	#0.02637317  

# Post-hoc comparisons
lsmeans(AggSize.lm3, pairwise ~ Species,type="response")



### Does the number of cells aggregated differ by relative testes size?

# Linear model:
AggSize.lm<-lm(Mean.Agg.Size ~ Testis.Mass..mg. + Body.Mass..mg., data= Pero.Species) 
summary(AggSize.lm)  

# Controlling for phylogeny using a PGLS:

AggSize.pgls1<-pgls(Mean.Agg.Size ~ Testis.Mass..mg. + Body.Mass..mg., data= Pero.PGLS, lambda="ML") 
summary(AggSize.pgls1)  


### Does the CV for the number of cells aggregated differ by relative testes size? 

# Linear model:

AggSizeCV.lm<-lm(Agg.Size.CV ~ Testis.Mass..mg. + Body.Mass..mg., data= Pero.Species) 
summary(AggSizeCV.lm)  

# Controlling for phylogeny using a PGLS: 

AggSizeCV.pgls1<-pgls(Agg.Size.CV ~ Testis.Mass..mg. + Body.Mass..mg., data= Pero.PGLS, lambda="ML") 
summary(AggSizeCV.pgls1)  


# What is the CV for aggregate size within species? 

tapply(Pero.Species$Agg.Size.CV, Pero.Species$Species, sum,na.rm=TRUE) 




##############    AGGREGATE MOTILITY     ################


##########  Does sperm aggregation confer motility benefits within species? [By species - MOTILE CELLS ONLY]


######## VCL 


###############    LOW VISCOSITY MEDIA

se = function(x) sd(x)/sqrt(length(x)) # calculates standard error

tapply(Pero.Male2$Single.VCL..low.viscosity., Pero.Male2$Species, mean)
tapply(Pero.Male2$Single.VCL..low.viscosity., Pero.Male2$Species, se)  

tapply(Pero.Male2$Agg.VCL..low.viscosity., Pero.Male2$Species, mean,na.rm=TRUE)
tapply(Pero.Male2$Agg.VCL..low.viscosity., Pero.Male2$Species, se)

# Calculate SE for P. leucopus
l=21
sd=sd(Pero.Male$Agg.VCL..low.viscosity.[Pero.Male$Species=="LL"], na.rm=TRUE)
se2=(sd/sqrt(l))
se2  

# Paired t-tests
t.test(IS2$Single.VCL..low.viscosity., IS2$Agg.VCL..low.viscosity., paired=TRUE) # NS (t = -1.0545, df = 28, p-value = 0.3007)
t.test(EP$Single.VCL..low.viscosity., EP$Agg.VCL..low.viscosity., paired=TRUE) # NS (t = 1.7308, df = 20, p-value = 0.09889)
t.test(PO$Single.VCL..low.viscosity., PO$Agg.VCL..low.viscosity., paired=TRUE) # Sig *** (t = 9.4575, df = 23, p-value = 2.166e-09)
t.test(BW$Single.VCL..low.viscosity., BW$Agg.VCL..low.viscosity., paired=TRUE) # Sig * (t = -2.2482, df = 17, p-value = 0.03812)
t.test(LL$Single.VCL..low.viscosity., LL$Agg.VCL..low.viscosity., paired=TRUE) # NS (t = 0.24811, df = 20, p-value = 0.8066)
t.test(GO$Single.VCL..low.viscosity., GO$Agg.VCL..low.viscosity., paired=TRUE) # Sig *** (t = 5.4048, df = 20, p-value = 2.73e-05)



###############    HIGH VISCOSITY MEDIA (0.75% Methylcellulose solution)

tapply(Pero.Male2$Single.VCL..high.viscosity., Pero.Male2$Species, mean,na.rm=TRUE)
tapply(Pero.Male2$Single.VCL..high.viscosity., Pero.Male2$Species, sd,na.rm=TRUE)  

# Calculate SE manually for each

# BW l = 12, sd = 15.02989
# EP l = 10, sd = 12.19113
# GO l = 11, sd = 15.38463
# IS l = 10, sd = 25.23098
# LL l = 10, sd = 10.84781
# PO l = 14, sd = 24.60915

se2=(24.60915/sqrt(14))
se2

tapply(Pero.Male2$Agg.VCL..high.viscosity., Pero.Male2$Species, mean,na.rm=TRUE)
tapply(Pero.Male2$Agg.VCL..high.viscosity., Pero.Male2$Species, sd, na.rm=TRUE)

# Calculate SE manually for each

# BW l = 12, sd = 19.41075
# EP l = 10, sd = 18.54995
# GO l = 11, sd = 19.35749
# IS l = 10, sd = 25.41022
# LL l = 10, sd = 15.66611
# PO l = 14, sd = 19.48304

se2=(19.48304/sqrt(14))
se2

 
# Paired t-tests
t.test(IS2$Single.VCL..high.viscosity., IS2$Agg.VCL..high.viscosity., paired=TRUE) # Sig diff ** (t = -3.9495, df = 9, p-value = 0.003357)
t.test(EP$Single.VCL..high.viscosity., EP$Agg.VCL..high.viscosity., paired=TRUE) # NS (t = -0.71739, df = 9, p-value = 0.4913)
t.test(PO$Single.VCL..high.viscosity., PO$Agg.VCL..high.viscosity., paired=TRUE) # NS (t = 2.0349, df = 13, p-value = 0.06278)
t.test(BW$Single.VCL..high.viscosity., BW$Agg.VCL..high.viscosity., paired=TRUE) # Sig diff * (t = -2.9397, df = 11, p-value = 0.01345)
t.test(LL$Single.VCL..high.viscosity., LL$Agg.VCL..high.viscosity., paired=TRUE) # NS (t = 0.32429, df = 9, p-value = 0.7531)
t.test(GO$Single.VCL..high.viscosity., GO$Agg.VCL..high.viscosity., paired=TRUE) # NS (t = -0.63835, df = 10, p-value = 0.5376)




##############    AGGREGATE COMPOSITION     ################


# Do species differ in the proportion of sperm aggregates they produce that are 'abnormal'? 

# Using a multiple proportions test
Abnormal<-c(187,266, 273, 107,128,58) # Number of abnormal aggregates by species - IS, EP, PO, BW, LL, GO
Aggregates<-c(858, 513, 950, 822, 361, 254) # Number of total aggregates by species - IS, EP, PO, BW, LL, GO
prop.test(Abnormal, Aggregates)

#pairwise comparisons
pairwise.prop.test(x = Abnormal, Aggregates)


## Does aggregate composition impact aggregate speed?  YES

# Paired t-tests
t.test(IS$NormalAgg.VCL, IS$AbnormalAgg.VCL, paired=TRUE) # Sig *** (t = 3.9666, df = 27, p-value = 0.0004835)
t.test(EP$NormalAgg.VCL, EP$AbnormalAgg.VCL, paired=TRUE) # Sig * (t = 2.8157, df = 19, p-value = 0.01104)
t.test(PO$NormalAgg.VCL, PO$AbnormalAgg.VCL, paired=TRUE) # Sig *** (t = 4.9961, df = 20, p-value = 6.935e-05)
t.test(BW$NormalAgg.VCL, BW$AbnormalAgg.VCL, paired=TRUE) # Sig *** (t = 5.3539, df = 15, p-value = 8.039e-05)
t.test(LL$NormalAgg.VCL, LL$AbnormalAgg.VCL, paired=TRUE) # NS (t = 0.8074, df = 14, p-value = 0.4329)
t.test(GO$NormalAgg.VCL, GO$AbnormalAgg.VCL, paired=TRUE) # NS (t = 1.1651, df = 17, p-value = 0.2601)








#########################################       Supplementary Materials     #########################################  


# Upload data in which one outlier P. californicus (IS) male is present (ID is 8997)

# Upload the "Aggregation by Species.csv" file - name it "Pero.Species"

# Upload the "Aggregation by Male.csv" file - name it "Pero.Male"

# Make species subset for IS (this contains the outlier male IS8997)
IS = Pero.Male2[which(Pero.Male2$Species=="IS"),]
   
# Combine phylogenetic tree and species data set with the outlier IS male
Pero.PGLS <- comparative.data(phy = Pero.Tree, data = Pero.Species, names.col = tip.label, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE) #Ignore warning message

# Upload the "Viscosity 1.5% by Male.csv" file - name it Viscosity2

# Make subsets of this data file by species: 
BW2 = Viscosity2[which(Viscosity2$Species=="BW"),] # This is Peromyscus maniculatus
EP2 = Viscosity2[which(Viscosity2$Species=="EP"),] # This is Peromyscus eremicus
GO2 = Viscosity2[which(Viscosity2$Species=="GO"),] # This is Peromyscus gossypinus
IS2 = Viscosity2[which(Viscosity2$Species=="IS"),] # This is Peromyscus californicus
LL2 = Viscosity2[which(Viscosity2$Species=="LL"),] # This is Peromyscus leucopus
PO2 = Viscosity2[which(Viscosity2$Species=="PO"),] # This is Peromyscus polionotus


# Upload the "Viscosity 2.25% by Male.csv" file - name it Viscosity3  

# Make subsets of this data file by species: 
BW3 = Viscosity3[which(Viscosity3$Species=="BW"),] # This is Peromyscus maniculatus
EP3 = Viscosity3[which(Viscosity3$Species=="EP"),] # This is Peromyscus eremicus
GO3 = Viscosity3[which(Viscosity3$Species=="GO"),] # This is Peromyscus gossypinus
IS3 = Viscosity3[which(Viscosity3$Species=="IS"),] # This is Peromyscus californicus
LL3 = Viscosity3[which(Viscosity3$Species=="LL"),] # This is Peromyscus leucopus
PO3 = Viscosity3[which(Viscosity3$Species=="PO"),] # This is Peromyscus polionotus


#############################################      METHODS       ############################################# 


## Age and standard deviation of males used in study
mean(Pero.Male$Age, na.rm=TRUE) 	## 193.4815
median(Pero.Male $Age, na.rm=TRUE) 	## 120
sd(Pero.Male $Age, na.rm=TRUE) 	## 138.3578
SE = function(x) sd(x)/sqrt(length(x))
SE(Pero.Male$Age) ## 11.90795

summary(Pero$Age<=120) # 70 are less than or equal to 120 days of age 
    

## Does age affect sperm production? # No with outliers (Pero), yes with them (HemoCounts); values reported for latter

Age.lm=lm(Hemo.Cell.Density.Estimate..cells.per.uL.~ Age, data = Pero.Male)
plot(Hemo.Cell.Density.Estimate..cells.per.uL.~ Age, data = Pero.Male)
abline(coef(Age.lm),col=2)
summary(Age.lm)
# Check assumptions
plot(resid(Age.lm)~fitted(Age.lm)) # Not heteroscedastic
qqnorm(residuals(Age.lm),pch=16,cex=1.2)
qqline(residuals(Age.lm),lwd=4,col=2) # Not normal

# Check normality of variables
hist(Pero.Male$Hemo.Cell.Density.Estimate..cells.per.uL.)
shapiro.test(Pero.Male$Hemo.Cell.Density.Estimate..cells.per.uL.) # not normal, W = 0.64274, p-value = 7.162e-16

# How many males had counts over 30,000
length(which(Pero.Male$Hemo.Cell.Density.Estimate..cells.per.uL. > 30000)) # n = 4

# Remove four outliers
HemoCounts = Pero.Male[which(Pero.Male$Hemo.Cell.Density.Estimate..cells.per.uL.< 30000),]
hist(HemoCounts$Hemo.Cell.Density.Estimate..cells.per.uL.) # normal
shapiro.test(HemoCounts$Hemo.Cell.Density.Estimate..cells.per.uL.) # normal, W = 0.99619, p-value = 0.9882

#Re-run linear model
Age.lm2=lm(Hemo.Cell.Density.Estimate..cells.per.uL.~ Age, data = HemoCounts) # p = 0.0165 * 
summary(Age.lm2)
# Check assumptions
plot(resid(Age.lm2)~fitted(Age.lm2)) # More heteroscedastic
qqnorm(residuals(Age.lm2),pch=16,cex=1.2)
qqline(residuals(Age.lm2),lwd=4,col=2) # Normal


## How many males paired in study (0 = no, 1 = yes)
summary(Pero.Male$Paired=="0") #Paired:  Yes = 32, No = 99, Unknown = 4


## Does pairing affect sperm production? # No, with or without outliers; reported values for latter

Paired.lm=lm(Hemo.Cell.Density.Estimate..cells.per.uL.~ as.factor(Paired), data = Pero.Male)
plot(Hemo.Cell.Density.Estimate..cells.per.uL.~ as.factor(Paired), data = Pero.Male)
abline(coef(Paired.lm),col=2)
summary(Paired.lm)
# Check assumptions
plot(resid(Paired.lm)~fitted(Paired.lm)) # Not heteroscedastic
qqnorm(residuals(Paired.lm),pch=16,cex=1.2) 
qqline(residuals(Paired.lm),lwd=4,col=2) # Not normal

#Re-run linear model without hemocytometer outliers
Paired.lm2=lm(Hemo.Cell.Density.Estimate..cells.per.uL.~ as.factor(Paired), data = HemoCounts) # p = 0.314 
summary(Paired.lm2)
# Check assumptions
plot(resid(Paired.lm2)~fitted(Paired.lm2)) # Heteroscedastic
qqnorm(residuals(Paired.lm2),pch=16,cex=1.2)
qqline(residuals(Paired.lm2),lwd=4,col=2) # Normal



## Are estimated sperm counts generated by CASA correlated with those collected manually from the hemocytometer? YES, p = 2.12e-05 ***

# Linear model
SpermDensity.lm=lm(CASA.Leja.Estimate..cells.per.mL.~Hemo.Cell.Density.Estimate..cells.per.uL.,data= Pero.Male)
summary(SpermDensity.lm) 
# Check assumptions
plot(resid(SpermDensity.lm)~fitted(SpermDensity.lm)) # Not heteroscedastic
qqnorm(residuals(SpermDensity.lm),pch=16,cex=1.2)
qqline(residuals(SpermDensity.lm),lwd=4,col=2) # Not normal, some outliers

#Re-run linear model without hemocytometer outliers
SpermDensity.lm2=lm(CASA.Leja.Estimate..cells.per.mL.~Hemo.Cell.Density.Estimate..cells.per.uL., data = HemoCounts) # p = 1.00e-07 *** 
summary(SpermDensity.lm2)
# Check assumptions
plot(resid(SpermDensity.lm2)~fitted(SpermDensity.lm2)) # Heteroscedastic
qqnorm(residuals(SpermDensity.lm2),pch=16,cex=1.2)
qqline(residuals(SpermDensity.lm2),lwd=4,col=2) # Normal

names(Pero.Male)


## Does sperm cell density significantly differ by species BEFORE dilutions? # Yes, with and without outliers (reported w/o outliers)

boxplot2(Sperm.Cell.Density.Before.Dilutions..cells.per.ul.~Species, col=c("red"), data= Pero.Male)  # seven outliers

# Linear model
PreDilutions.lm=lm(Sperm.Cell.Density.Before.Dilutions..cells.per.ul.~Species,data= Pero.Male) # Yes, between BW, GO, LL
summary(PreDilutions.lm) 
# Check assumptions
plot(resid(PreDilutions.lm)~fitted(PreDilutions.lm)) # Not heteroscedastic
qqnorm(residuals(PreDilutions.lm),pch=16,cex=1.2)
qqline(residuals(PreDilutions.lm),lwd=4,col=2) # Not normal

# Check variable for normality
hist(Pero.Male$Sperm.Cell.Density.Before.Dilutions..cells.per.ul.) # Not normally distributed, some outliers

# How many males had counts over 300,000?
length(which(Pero.Male$Sperm.Cell.Density.Before.Dilutions..cells.per.ul. > 300000)) # n = 4

# Remove four outliers
Dilutions = Pero.Male[which(Pero.Male$Sperm.Cell.Density.Before.Dilutions..cells.per.ul.< 300000),]
hist(Dilutions$Sperm.Cell.Density.Before.Dilutions..cells.per.ul.) # not normal
shapiro.test(Dilutions$Sperm.Cell.Density.Before.Dilutions..cells.per.ul.) # not normal, W = 0.81055, p-value = 5.687e-11


# Try a transformation
hist(log(Dilutions$Sperm.Cell.Density.Before.Dilutions..cells.per.ul.)) # more normal
shapiro.test(log(Dilutions$Sperm.Cell.Density.Before.Dilutions..cells.per.ul.)) # normal, W = 0.97891, p-value = 0.06228

# Linear model with transformation and removing outliers
PreDilutions.lm2=lm(log(Sperm.Cell.Density.Before.Dilutions..cells.per.ul.)~Species,data= Dilutions)
summary(PreDilutions.lm2) 
# Check assumptions
plot(resid(PreDilutions.lm2)~fitted(PreDilutions.lm2)) # Somewhat heteroscedastic
qqnorm(residuals(PreDilutions.lm2),pch=16,cex=1.2)
qqline(residuals(PreDilutions.lm2),lwd=4,col=2) # Also not normal, though more so




## Does sperm cell density significantly differ by species AFTER dilutions? 

boxplot2(Sperm.Cell.Density.After.Dilutions..cells.per.ul.~Species, col=c("turquoise"), data= Pero.Male) # four outliers

# Linear model
PostDilutions.lm=lm(Sperm.Cell.Density.After.Dilutions..cells.per.ul.~Species,data= Pero.Male)
summary(PostDilutions.lm)  
# Check assumptions
plot(resid(PostDilutions.lm)~fitted(PostDilutions.lm)) # not heteroscedastic
qqnorm(residuals(PostDilutions.lm),pch=16,cex=1.2)
qqline(residuals(PostDilutions.lm),lwd=4,col=2) # Not normal

# Check variable for normality/outliers
hist(Pero.Male$Sperm.Cell.Density.After.Dilutions..cells.per.ul.) # Not normally distributed, some outliers
shapiro.test(Pero.Male$Sperm.Cell.Density.After.Dilutions..cells.per.ul.) # Not normal, W = 0.65197, p-value = 1.151e-15

# How many males had counts over 30,000
length(which(Pero.Male$Sperm.Cell.Density.After.Dilutions..cells.per.ul. > 30000)) # n = 3

# Remove three outliers
Dilutions2 = Pero.Male[which(Pero.Male$Sperm.Cell.Density.After.Dilutions..cells.per.ul.< 30000),]
hist(Dilutions2$Sperm.Cell.Density.After.Dilutions..cells.per.ul.) # normal
shapiro.test(Dilutions2$Sperm.Cell.Density.After.Dilutions..cells.per.ul.) # normal, W = 0.99641, p-value = 0.9921

# Linear model without outliers
PostDilutions.lm2=lm(Sperm.Cell.Density.After.Dilutions..cells.per.ul.~Species,data=Dilutions2)
summary(PostDilutions.lm2)  
# Check assumptions
plot(resid(PostDilutions.lm2)~fitted(PostDilutions.lm2)) # More heteroscedastic
qqnorm(residuals(PostDilutions.lm2),pch=16,cex=1.2)
qqline(residuals(PostDilutions.lm2),lwd=4,col=2) # Normal


# Post-hoc comparisons
lsmeans(PostDilutions.lm2, pairwise ~ Species,type="response")



## How many samples fell above and below the ideal sperm concentration?
summary(Pero.Male$CASA.Leja.Count..Sum.<"300") # n = 40
summary(Pero.Male$CASA.Leja.Count..Sum.>"400") # n = 21





#############################################      RESULTS       ############################################# 



##############    AGGREGATE SIZE     ################


### Do species vary in the sizes of sperm aggregates?

# Determine random factors
AggSize.lmer=lmer(Mean.Agg.Size ~ + (1 | Family),data= Pero.Male)
summary(AggSize.lmer)

#Run bivariate analysis of traits

# Rescale Variable
Pero.Male$Age2 <- scale(Pero.Male$Age, center = TRUE, scale = TRUE)
Pero.Male$Time <- scale(Pero.Male$Post.Harvest.Time..min., center = TRUE, scale = TRUE)

AggSize.lmer=lmer(Mean.Agg.Size ~ Age2 + (1 | Family),data= Pero.Male)
summary(AggSize.lmer)

#Age2,  0.0372 *
#as.factor(Paired),  0.00226 **
#Time, 0.0241 *
#Total.Sperm.Cells, 3.37e-06 ***
#Sperm.Cells.per.Video, 0.0826
#Number.of.CASA.Videos, 0.000674 ***
#as.factor(Species), all *** except PO 0.981090

# Fixed factors to consider in final model
#Age2 + as.factor(Paired) + Time + Total.Sperm.Cells + Sperm.Cells.per.Video + Number.of.CASA.Videos

# Check collinearities using a linear model
lm1=lm(Time ~ as.factor(Species),data= Pero.Male)
summary(lm1)

# The following factors are collinear: 
# Age2  ~ as.factor(Species)					* use latter
# Time  ~ as.factor(Species)					* use latter
# Sperm.Cells.per.Video ~ Total.Sperm.Cells  	* use latter
# as.factor(Paired) ~ Total.Sperm.Cells 		* use latter
# Number.of.CASA.Videos ~ as.factor(Species) 	* use latter


# Final model selection

AggSize.lmer1=lmer(Mean.Agg.Size ~ Total.Sperm.Cells + as.factor(Species) + (1 | Family),data = Pero.Male)
summary(AggSize.lmer1)

# Family has a minimal effect on response variable, converting to a linear model
AggSize.lm1=lm(Mean.Agg.Size ~ Total.Sperm.Cells + as.factor(Species),data= Pero.Male)
summary(AggSize.lm1)

#Check model assumptions
# Check for heteroscedasticity
plot(resid(AggSize.lm1)~fitted(AggSize.lm1)) # Not heteroscedastic
# Check Normality 
qqnorm(residuals(AggSize.lm1),pch=16,cex=1.2)
qqline(residuals(AggSize.lm1),lwd=4,col=2) # Not normal, giant outlier

# Are the data normal? 
hist(Pero.Male$Mean.Agg.Size)
shapiro.test(Pero.Male$Mean.Agg.Size) # No, W = 0.80965, p-value = 6.708e-12
# Try a log transformation of data
hist(log(Pero.Male$Mean.Agg.Size))
shapiro.test(log(Pero.Male$Mean.Agg.Size)) # No, but more normal, W = 0.90399, p-value = 8.868e-08

# Are the data normal? 
hist(Pero.Male$Total.Sperm.Cells)
shapiro.test(Pero.Male$Total.Sperm.Cells) # Yes, W = 0.99045, p-value = 0.4857

# Rerun with transformed data
AggSize.lm2=lm(log(Mean.Agg.Size) ~ Total.Sperm.Cells + as.factor(Species),data= Pero.Male)
summary(AggSize.lm2)

#Check model assumptions
# Check for heteroscedasticity
plot(resid(AggSize.lm2)~fitted(AggSize.lm2)) # More heteroscedastic
# Check Normality 
qqnorm(residuals(AggSize.lm2),pch=16,cex=1.2)
qqline(residuals(AggSize.lm2),lwd=4,col=2) # Still not normal, giant outlier


### Is there more variance across than within species?  YES

var(Pero.Male$Mean.Agg.Size,na.rm=TRUE) #1.961014

var(IS$Mean.Agg.Size,na.rm=TRUE) 	#2.900569
var(EP$Mean.Agg.Size) 				#0.3059087
var(PO$Mean.Agg.Size,na.rm=TRUE) 	#2.150459    
var(BW$Mean.Agg.Size,na.rm=TRUE) 	#0.8228839
var(LL$Mean.Agg.Size,na.rm=TRUE) 	#0.121562
var(GO$Mean.Agg.Size) 			 	#0.02637317  

# Post-hoc comparisons
lsmeans(AggSize.lm1, pairwise ~ Species,type="response")


### Does the number of cells aggregated differ by relative testes size when controlling for phylogeny using a PGLS?
AggSize.pgls1<-pgls(Mean.Agg.Size ~ Testis.Mass..mg. + Body.Mass..mg., data= Pero.PGLS, lambda="ML") 
summary(AggSize.pgls1) 


### Does the CV for the number of cells aggregated differ by relative testes size when controlling for phylogeny using a PGLS? 
AggSizeCV.pgls1<-pgls(Agg.Size.CV ~ Testis.Mass..mg. + Body.Mass..mg., data= Pero.PGLS, lambda="ML") 
summary(AggSizeCV.pgls1)



# What is the CV for aggregate size within species? 

tapply(Pero.Species$Agg.Size.CV, Pero.Species$Species, sum,na.rm=TRUE) 



##############    AGGREGATE MOTILITY     ################


##########  Does sperm aggregation confer motility benefits within species? [By species - MOTILE CELLS ONLY]



###############    HIGH VISCOSITY MEDIA (1.5% Methylcellulose solution)

se = function(x) sd(x)/sqrt(length(x)) # calculates standard error

tapply(Viscosity2$Single.VCL..high.viscosity.1.5.., Viscosity2$Species, mean)
tapply(Viscosity2$Single.VCL..high.viscosity.1.5.., Viscosity2$Species, se)  

tapply(Viscosity2$Agg.VCL..high.viscosity.1.5.., Viscosity2$Species,mean)
tapply(Viscosity2$Agg.VCL..high.viscosity.1.5.., Viscosity2$Species,se)

# Paired t-tests
t.test(IS2$Single.VCL..high.viscosity.1.5.., IS2$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # Sig diff * (t = -3.1371, df = 5, p-value = 0.02575)
t.test(EP2$Single.VCL..high.viscosity.1.5.., EP2$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # Sig diff * (t = -2.4215, df = 7, p-value = 0.04599)
t.test(PO2$Single.VCL..high.viscosity.1.5.., PO2$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # NS (t = 4.8314, df = 13, p-value = 0.0003279)
t.test(BW2$Single.VCL..high.viscosity.1.5.., BW2$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # Sig diff *** (t = -5.2835, df = 11, p-value = 0.0002589)
t.test(LL2$Single.VCL..high.viscosity.1.5.., LL2$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # NS (t = -0.50583, df = 8, p-value = 0.6266)
t.test(GO2$Single.VCL..high.viscosity.1.5.., GO2$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # Sig diff *  (t = -2.7834, df = 5, p-value = 0.03875)



###############    HIGH VISCOSITY MEDIA (2.25% Methylcellulose solution)

se = function(x) sd(x)/sqrt(length(x)) # calculates standard error

tapply(Viscosity3$Single.VCL..high.viscosity.2.25.., Viscosity3$Species, mean)
tapply(Viscosity3$Single.VCL..high.viscosity.2.25.., Viscosity3$Species, se)  

tapply(Viscosity3$Agg.VCL..motile., Viscosity3$Species,mean, na.rm=TRUE)
tapply(Viscosity3$Agg.VCL..motile., Viscosity3$Species,sd, na.rm=TRUE)

# Calculate SE manually for each

# EP l = 8, sd = 11.95143 (1 did not produce aggregates)
# GO l = 6, sd = 11.19377 (3 did not produce aggregates)
# IS l = 12, sd = 17.59711 (3 did not produce aggregates)
# LL l = 8, sd = 23.96119 (6 did not produce aggregates)

se2=(23.96119/sqrt(8))
se2

# Paired t-tests
t.test(IS3$Single.VCL..high.viscosity.2.25.., IS3$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # Sig diff ** (t = -3.164, df = 11, p-value = 0.009015)
t.test(EP3$Single.VCL..high.viscosity.2.25.., EP3$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # NS (t = -1.179, df = 7, p-value = 0.2769)
t.test(PO3$Single.VCL..high.viscosity.2.25.., PO3$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # NS (t = -2.156, df = 12, p-value = 0.05208)
t.test(BW3$Single.VCL..high.viscosity.2.25.., BW3$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # Sig diff *** (t = -5.3274, df = 11, p-value = 0.000242)
t.test(LL3$Single.VCL..high.viscosity.2.25.., LL3$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # Sig diff * (t = -2.4498, df = 7, p-value = 0.04412)
t.test(GO3$Single.VCL..high.viscosity.2.25.., GO3$Agg.VCL..high.viscosity.1.5.., paired=TRUE) # NS (t = 0.69638, df = 5, p-value = 0.5172)




