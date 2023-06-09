#### MFPS Severity Field data 2015 ######

## Updated 3-20-2020

data<-read.csv("MFPS_sev_3-18-20.csv")
head(data)
length(data$depth) # 509 with MFPS

hist(data$mfps_sev) # censored right tail, cap at 50

#### Data classes 

# Numeric
data$mfps_sev<-as.numeric(data$mfps_sev) 
class(data$mfps_sev)
data$temp<-as.numeric(data$temp) 
class(data$temp) 
class(data$temp_sd) # numeric
data$depth<-as.numeric(data$depth)
class(data$depth) # numeric
class(data$coralcov) #numeric
data$sf_num<-as.numeric(data$sf_num)
class(data$sf_num) # 
data$size<-as.numeric(data$size) # 
class(data$size) # numeric

# Factors
class(data$site) # factor
class(data$cu_all) # factor 
class(data$mfps_fan) 

# CHECKING Data
plot(fan_num~site,data=data)
plot(transect~site,data=data) # 3 per
plot(temp~site,data=data) # 1
plot(temp_sd~site,data=data) # 1
plot(sf_num~site,data=data) # 1
plot(cu_all~site,data=data) # 1
plot(coralcov~transect,data=data) # 1 each
plot(depth~transect,data=data) # 1 each
plot(size~fan_num,data=data) # 1 each

## TEST FOR COLLINEARITY 
head(data)
round(cor(data[,c(6,7,8,9,11,12)]),2)  # only temp metrics, as expected

#### Rescale numeric predictors
library(arm)

data$temp<-rescale(data$temp, binary.inputs="center")
data$temp_sd<-rescale(data$temp_sd, binary.inputs="center")
data$coralcov<-rescale(data$coralcov, binary.inputs="center")
data$depth<-rescale(data$depth, binary.inputs="center")
data$size<-rescale(data$size, binary.inputs="center")
data$sf_num<-rescale(data$sf_num, binary.inputs="center")


### Contrasts - factors ## 
contrasts(data$site) 
contrasts(data$cu_all) 
contrasts(data$cu_all) <- cbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1)) # comparing all to level 1


###############
#### MODELS ###
###############

# Make dataframe
logspots<-log(data$mfps_sev)
dataf<-data.frame(logspots,data$site, data$size,data$mfps_sev,data$temp, data$temp_sd, data$cu_all,data$coralcov,data$depth,data$sf_num)
head(dataf)

library(plm) 
head(dataf)

## Plan for analysis ######
# Censreg is meant to deal with truncated data
# (1) First run with data not transformed
# (2) Run same models with lmer and see if it's overdispersed; log transform if so to show that fixes overdipersion 
# (3) Log transform response variable and re-run models in censReg (tobit) if it is overdispersed

library(censReg) #


######### 1: SPOTS AS DEPENDENT VARIABLE, no log ####

# Predictors -  size, coralcov, temp, cu_all, depth, sf_num; separate with site

# Null
sev_null<-censReg(data.mfps_sev ~ 1, data = dataf,left = 0, right = 50, method = "BHHH" )

# All factors, additive
sev1<-censReg(data.mfps_sev ~ data.size+data.coralcov+data.sf_num+data.temp+data.cu_all+data.depth, data=dataf, left = 0, right = 50, method = "BHHH" ) # BHHH is ML method; also estimates random effects
drop1(sev1) # drop depth
sev1b<-censReg(data.mfps_sev ~ data.size+data.coralcov+data.sf_num+data.temp+data.cu_all, data=dataf, left = 0, right = 50, method = "BHHH" ) # BHHH is ML method; also estimates random effects
drop1(sev1b) # drop none

# Overdispersion? 
library(lme4)      
model<-glm(mfps_sev~size+coralcov+sf_num+temp+cu_all, family="poisson", data= data)
summary(model) # overdispersed based on residual deviance

# does log transform help?
model2<-glm(logspots~data.size+data.coralcov+data.sf_num+data.cu_all, family="poisson", data= dataf)
summary(model2) # yes.

## SO: Log transformation accomplishes the same goal as the quasipoisson as censReg does not take the 'quasipoisson'

############## 2: RUN MODELS WITH LMER ######## 

# See if it's overdispersed; log transform to show that fixes overdipersion and has similar standard errors as using quasipoisson correction ###
library(lme4)

# 1: SE for poisson 
spot_mod1<-glm(mfps_sev~size+coralcov+sf_num+cu_all,data=data,family="poisson") # 
summary(spot_mod1) # Size SE = 0.0209; overdispersed

# 2: SE for Quasi 
spot_mod1b<-glm(mfps_sev~size+coralcov+sf_num+cu_all,data=data,family="quasipoisson") 
summary(spot_mod1b) # Size SE = 0.0806

# 3: SE for log 
spot_mod1c<-lm(logspots~size+coralcov+sf_num+cu_all,data=data) 
summary(spot_mod1c) # Size SE 0.098 (similar to quasi- good)

# SO: log transformation approximately reproduces the SE of the quasi models 

# Homogeneity of variance, logspots response
plot(logspots~site, data=data)
plot(logspots~cu_all, data=data)
plot(logspots~size, data=data)
plot(logspots~temp, data=data)
plot(logspots~temp_sd, data=data)
plot(logspots~depth, data=data)
plot(logspots~coralcov, data=data)
plot(logspots~sf_num, data=data)

########## SECTION 3: LOG-SPOTS AS DEPENDENT VARIABLE, CensReg ###########

# Predictors -  size, coralcov, temp, temp_sd, cu_all, depth, sf_num; separate with site

# Right limit of censReg is now 3.92, not 50, because on log scale
library(censReg)
log(50) # 3.92 bounded
lognull<-censReg(logspots~1,data=dataf, left=-Inf, right=3.92,method="BHHH")

### Additive ###

# Temp
log1<-censReg(logspots~data.size+data.coralcov+data.temp+data.cu_all+data.depth+data.sf_num,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log1) # drop sf_num
log1b<-censReg(logspots~data.size+data.coralcov+data.temp+data.cu_all+data.depth,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log1b) # drop depth
log1c<-censReg(logspots~data.size+data.coralcov+data.temp+data.cu_all,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log1c) # drop coralcov
log1d<-censReg(logspots~data.size+data.temp+data.cu_all,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log1d) # drop none

# Temp_sd
log2<-censReg(logspots~data.size+data.coralcov+data.temp_sd+data.cu_all+data.depth+data.sf_num,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log2) # drop sf_num
log2b<-censReg(logspots~data.size+data.coralcov+data.temp_sd+data.cu_all+data.depth,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log2b) # drop depth
log2c<-censReg(logspots~data.size+data.coralcov+data.temp_sd+data.cu_all,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log2c) # drop temp_sd
log2d<-censReg(logspots~data.size+data.coralcov+data.cu_all,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log2d) # drop none

# Site (separate because interferes with convergence)
log7<-censReg(logspots~data.coralcov+data.site,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log7) # drop none

log8<-censReg(logspots~data.cu_all+data.site,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log8) # doesn't converge

log9<-censReg(logspots~data.sf_num+data.site,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log9)# doesn't converge

log10<-censReg(logspots~data.size+data.site,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log10) # drop none

log11<-censReg(logspots~data.temp+data.site,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log11) # doesn't converge

log12<-censReg(logspots~data.temp_sd+data.site,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log12)# doesn't converge

log13<-censReg(logspots~data.depth+data.site,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log13) # drop depth
log13b<-censReg(logspots~data.site,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log13b) 

####### QUADRATICS ########

# Size temp 
log4<-censReg(logspots~data.sf_num+data.size+I(data.size^2)+data.coralcov+data.temp+data.cu_all+data.depth,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log4) # drop sf_num
log4b<-censReg(logspots~data.size+I(data.size^2)+data.coralcov+data.temp+data.cu_all+data.depth,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log4b) # drop depth
log4c<-censReg(logspots~data.size+I(data.size^2)+data.coralcov+data.temp+data.cu_all,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log4c) # drop coralcov
log4d<-censReg(logspots~data.size+I(data.size^2)+data.temp+data.cu_all,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log4d) # drop none

# Size temp_sd
log5<-censReg(logspots~data.size+I(data.size^2)+data.coralcov+data.temp_sd+data.cu_all+data.depth+data.sf_num,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log5) # drop depth, sf_num
log5b<-censReg(logspots~data.size+I(data.size^2)+data.coralcov+data.temp_sd+data.cu_all,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log5b) # drop temp
log5c<-censReg(logspots~data.size+I(data.size^2)+data.coralcov+data.cu_all,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log5c) # drop none

# Size site
log6<-censReg(logspots~data.size+I(data.size^2)+data.site,data=dataf, left=-Inf, right=3.92,method="BHHH")
drop1(log6) # none


#### AIC ####
AIC(lognull) # 1536.735 (df=2)
AIC(log1d) # 1500.07 (df=7)
AIC(log2d) # 1504.327 (df=7)
AIC(log7) #1512.025 (df=17)
AIC(log10) # 1484.102 (df=17)
AIC(log13b) # 1514.757 (df=16)
AIC(log4d) #  1497.401 (df=8)
AIC(log5c) #  1501.405 (df=8)
AIC(log6) # 1479.695 (df=18)


# LRT - can use because top models are nested
library(lmtest)
lrtest(log6,log10) # log 6 merits quadratic term
lrtest(log6,lognull) # log 6 is best

### BEST MODEL - log6
summary(log6) # 

# Conf. intervals
library(MuMIn)
confint(log6,level=0.95, method="Wald") 

