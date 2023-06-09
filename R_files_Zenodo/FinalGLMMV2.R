#### Habitat use GLMM for Bottlenose and Humpback dolphin 
### Feb 2022
### Jono B and Ale V

# Load packages
library(lme4)
library(lmerTest)
library(MuMIn)
library(usdm)

############################################## Bottlenose dolphins ##############################################################
#Set working directory
#setwd("C:/Users/Jonathan/Desktop/Ale-Resource Partitioning/Models/Final GLMM and data sets")
#setwd("~/Jonathan/5.WorkingFiles/Vargus Dolphin Resource Partitioning/Models/2022-Updated")
setwd("C:/Users/botha/Desktop/Ale models")
rm(list=ls()) # clear workspace

# Read in the data
Data<-read.csv(file ="ModelData.csv", head = T)
head(Data)
str(Data)

Data$AMPM <- ifelse(Data$timeBroad < 12,"AM","PM")

# Set the data up 
head(Data)
Data$Season <- as.factor(Data$Season)
Data$Poly <- as.factor(Data$Polygon)
#Data$Time.h <- as.factor(Data$Time.h)
#Data$Time.2h <- as.factor(Data$Time.2h)
#Data$Time.Bin <- as.factor(Data$Time.Bin)
Data$AMPM <- as.factor(Data$AMPM)
Data$MPA <- as.factor(Data$MPA)
Data$Hab <- as.factor(Data$Substratum.inshore)
Data$AReef <- as.factor(Data$Agulhas.reef)
Data$Mosaic <- as.factor(Data$Mosaic.substratum)
Data$EST <- as.factor(Data$ESTUARY)
Data$Temp <- as.numeric(Data$Temp)
Data$EncY.N <- as.factor(Data$EncY.N)

# Check the data
head(Data)
str(Data)

Data$BNDEnc <- as.factor(ifelse(Data$BND.HBD=="bnd",1,0))
Data$HBDEnc <- as.factor(ifelse(Data$BND.HBD=="hbd",1,0))

# Check correlations of fixed variables
head(Data)
viftable <- as.data.frame(cbind(Data$Hab,Data$AMPM,Data$MPA,Data$Season,Data$AReef,Data$Mosaic,Data$EST,Data$Temp,Data$Poly))
#simplify headers
colnames(viftable) <- c("Habitat","Time","mpa","Season","AReef","Mosaic","Estuary","Temp","Poly")
vifcor(viftable)

#########################################################################################################################
############################################# Bottlenose Dolphin Model ##################################################
# The full model including all predictor effects
head(Data)
BDM1=glmer(BNDEnc ~ Hab*AMPM + AReef + Mosaic + EST + MPA + Season + Temp + (1|Poly), Data,family="binomial")
summary(BDM1)

# Model selection
BDM2 <- refitML(BDM1) # refit the model using maximum likelihood estimation
summary(BDM2)

# Model selection - AIC scores
options(na.action=na.fail)
AICBD<- dredge(BDM2,rank=AICc)
summary(AICBD)
head(AICBD)
AICBD
options(na.action = na.omit)
# Hav a look at the first model
top_model <- get.models(AICBD, subset = 1)[[1]]
top_model
summary(top_model)

# Model averaging
AICBDavg <- model.avg(AICBD, subset = delta <= 2)
summary(AICBDavg)
confint(AICBDavg)
as.data.frame(importance(full.model.avg))

summary(model.avg(AICBD, subset = delta <= 2))

#as a 95% confidence set:
subs <- subset(AICBD,cumsum(weight) <= .95)
avg.coeffs <- model.avg(subs,fit=T) # get averaged coefficients
avg.coeffs
summary(avg.coeffs)
#'Best' model
summary(get.models(AIC, 1)[[1]]) ######## no interaction is best model - but interaction retained in 95% set
summary(avg.coeffs)
confint(avg.coeffs)
# Most parsimonious model following model averaging
qqnorm(residuals(BDM1)) # Check normality of residuals
shapiro.test(residuals(BDM1)) 
plot(residuals(BDM1)~fitted(BDM1))
r.squaredGLMM(BDM1)

#########################################################################################################################
############################################### Humpback Dolphin Model ##################################################
# The full model including all predictor effects
head(Data)
HDM1=glmer(HBDEnc ~ Hab*AMPM + AReef + Mosaic + EST + MPA + Season + Temp + (1|Poly), Data,family="binomial")
summary(HDM1)

# Model selection
HDM2 <- refitML(HDM1) # refit the model using maximum likelihood estimation
summary(HDM2)

# Model selection - AIC scores
options(na.action=na.fail)
AICHD<- dredge(HDM2,rank=AICc)
summary(AICHD)
head(AICHD)
AICHD
options(na.action = na.omit)
# Hav a look at the first model
top_model <- get.models(AICHD, subset = 1)[[1]]
top_model
summary(top_model)

# Model averaging
AICHDavg <- model.avg(AICHD, subset = delta <= 2)
summary(AICHDavg)
confint(AICHDavg)

summary(model.avg(AICHD, subset = delta <= 2))

#as a 95% confidence set:
subs <- subset(AICHD,cumsum(weight) <= .95)
avg.coeffs <- model.avg(subs,fit=T) # get averaged coefficients
avg.coeffs
summary(avg.coeffs)
#'Best' model
summary(get.models(AICHD, 1)[[1]]) ######## no interaction is best model - but interaction retained in 95% set
summary(avg.coeffs)
confint(avg.coeffs)
# Most parsimonious model following model averaging
qqnorm(residuals(HDM1)) # Check normality of residuals
shapiro.test(residuals(HDM1)) 
plot(residuals(HDM1)~fitted(HDM1))
r.squaredGLMM(HDM1)


