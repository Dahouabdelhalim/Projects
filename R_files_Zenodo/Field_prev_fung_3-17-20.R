### Field prevalence 2015 field: fungus ###

# Fungal hyphae PREVALENCE from histology slides ##
# 90 data points (6 colonies * 15 sites)

data<-read.csv("Field_prev_90points_3-17-20.csv")
head(data)
length(data$depth) # 90 points

# Normality/ distribution - binomial response
summary(data$fung_fan) 
12/(90) # 13.3% prevalence of fungus

### For Chi square analyses ###
table(data$fung_fan,data$cop_fan)
table(data$fung_fan,data$mfps_yn)

#### REGRESSORS SET UP ########

# Numeric
data$temp_mean<-as.numeric(data$temp_mean) 
class(data$temp_mean) 
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
class(data$sex)
class(data$fung_fan)

# CHECKING Data
plot(slide~site,data=data)
plot(colony~site,data=data)
plot(temp_mean~site,data=data)
plot(temp_sd~site,data=data)
plot(sf_num~site,data=data)
plot(cu_all~site,data=data)
plot(coralcov~transect,data=data)
plot(depth~transect,data=data)
plot(size~colony,data=data)
plot(sex~colony,data=data)
plot(fung_fan~colony,data=data)
plot(sex~site,data=data)

# EXPLORATORY plots 
plot(fung_fan~site,data=data)  # some sites have no fungi
plot(fung_fan~cu_all,data=data) 
plot(fung_fan~sex,data=data) 
plot(fung_fan~size,data=data) 
plot(fung_fan~coralcov,data=data) 
plot(fung_fan~sf_num,data=data) 
plot(fung_fan~depth,data=data)
plot(fung_fan~temp_mean,data=data)
plot(fung_fan~temp_sd,data=data)

## TEST FOR COLLINEARITY 
head(data)
round(cor(data[,c(8,9,10,11,12,15)]),2)  

#### Rescale numeric predictors
library(arm)

data$temp_mean<-rescale(data$temp_mean, binary.inputs="center")
data$temp_sd<-rescale(data$temp_sd, binary.inputs="center")
data$coralcov<-rescale(data$coralcov, binary.inputs="center")
data$depth<-rescale(data$depth, binary.inputs="center")
data$size<-rescale(data$size, binary.inputs="center")
data$sf_num<-rescale(data$sf_num, binary.inputs="center")


### Contrasts - factors ## 
contrasts(data$site) 
contrasts(data$cu_all) 
contrasts(data$cu_all) <- cbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1)) # comparing all to level 1
contrasts(data$sex) <- cbind(c(-1,1,0),c(1,1,-2)) 
contrasts(data$sex) 

contrasts(data$fung_fan) # Y is 1

###############
#### MODELS ###
###############

# Response: fung_fan (binary)
# Predictors: site, temp_mean, temp_sd, depth, coralcov, sf_num, cu_all, size, sex, cop_fan

library(lme4)
library(AICcmodavg)

null<-glm(fung_fan~1,data=data,family='binomial')

#### Additive
mod1<-glm(fung_fan~sex+temp_mean+site+cop_fan+depth+coralcov+temp_sd+sf_num+cu_all+size,data=data,family='binomial')
drop1(mod1) # drop site, warning on fit for only this model 
mod1b<-glm(fung_fan~sex+temp_mean+cop_fan+depth+coralcov+temp_sd+sf_num+cu_all+size,data=data,family='binomial')
drop1(mod1b) # drop cu
mod1c<-glm(fung_fan~sex+temp_mean+cop_fan+depth+coralcov+temp_sd+sf_num+size,data=data,family='binomial')
drop1(mod1c) # drop sex
mod1d<-glm(fung_fan~temp_mean+cop_fan+depth+coralcov+temp_sd+sf_num+size,data=data,family='binomial')
drop1(mod1d) # drop size
mod1e<-glm(fung_fan~temp_mean+cop_fan+depth+coralcov+temp_sd+sf_num,data=data,family='binomial')
drop1(mod1e) # drop depth
mod1f<-glm(fung_fan~temp_mean+cop_fan+coralcov+temp_sd+sf_num,data=data,family='binomial')
drop1(mod1f) # drop sf_num
mod1g<-glm(fung_fan~temp_mean+cop_fan+coralcov+temp_sd,data=data,family='binomial')
drop1(mod1g) # drop temp_mean
mod1h<-glm(fung_fan~cop_fan+coralcov+temp_sd,data=data,family='binomial')
drop1(mod1h) # drop coralcov
mod1i<-glm(fung_fan~cop_fan+temp_sd,data=data,family='binomial')
drop1(mod1i) # drop temp_sd
mod1j<-glm(fung_fan~cop_fan,data=data,family='binomial')
drop1(mod1j) # drop none

# LRT to test vs null
library(lmtest)
lrtest(null,mod1j) # p = 0.134, Chi = 2.24; NS

# Null is the best model.

#### Additive with mfps instead of copepod 
mod2<-glm(fung_fan~sex+temp_mean+site+mfps_yn+depth+coralcov+temp_sd+sf_num+cu_all+size,data=data,family='binomial')
drop1(mod2) # drop site, warnings for this model only 
mod2b<-glm(fung_fan~sex+temp_mean+mfps_yn+depth+coralcov+temp_sd+sf_num+cu_all+size,data=data,family='binomial')
drop1(mod2b) # drop cu_all
mod2c<-glm(fung_fan~sex+temp_mean+mfps_yn+depth+coralcov+temp_sd+sf_num+size,data=data,family='binomial')
drop1(mod2c) # drop sex
mod2d<-glm(fung_fan~temp_mean+mfps_yn+depth+coralcov+temp_sd+sf_num+size,data=data,family='binomial')
drop1(mod2d) # drop size
mod2e<-glm(fung_fan~temp_mean+mfps_yn+depth+coralcov+temp_sd+sf_num,data=data,family='binomial')
drop1(mod2e) # drop depth
mod2f<-glm(fung_fan~temp_mean+mfps_yn+coralcov+temp_sd+sf_num,data=data,family='binomial')
drop1(mod2f) # drop mfps
mod2g<-glm(fung_fan~temp_mean+coralcov+temp_sd+sf_num,data=data,family='binomial')
drop1(mod2g) # drop sf_num
mod2h<-glm(fung_fan~temp_mean+coralcov+temp_sd,data=data,family='binomial')
drop1(mod2h) # drop temp_mean
mod2i<-glm(fung_fan~coralcov+temp_sd,data=data,family='binomial')
drop1(mod2i) # drop coralcov
mod2j<-glm(fung_fan~temp_sd,data=data,family='binomial')
drop1(mod2j)# drop none

# LRT to test vs null
library(lmtest)
lrtest(null,mod2j) # p = 0.1572, NS

# Null is the best model.
