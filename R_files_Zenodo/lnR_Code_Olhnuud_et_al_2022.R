
# # # # # # # # # # # # # # # # # 
# # # Pollination deficit # # # # 
#   --- Aruhan, Wopke --- #
#      - 2021.02.02 -    #


# clean the environment #
rm(list=ls())

# set work space #
setwd("F:/R/RStudio/R/r_work/Meta-analysis_Aruhan")

# install packages #
install.packages("nlme")
install.packages("car")
install.packages("emmeans")
install.packages("ggplot2")
install.packages("gridExtra")

# load libraries #
library(nlme) # for fit linear mixed models
library(ggplot2) # to plot forest plot
library(emmeans) # to compare treatment means in analysis with moderator
library(car) # to plot qqPlot
library(multcomp) # to do multiple comparison
library(gridExtra) 

# load dataset #
PL_ALL<- read.csv("Pollination deficits_DATA_Olhnuud et al.csv")
# We would extract required columns from raw database.
names(PL_ALL)
PL_data <- PL_ALL[,c(2,3,8,9,10,11,12,13,14,15,16,17,18,19,20)]
names(PL_data)

# Calculate effect size lnRPD
PL_data$R = PL_data$Xopen/PL_data$Xartificial
PL_data$lnR = log(PL_data$R)
range(PL_data$lnR) # -2.094, 1.174
names(PL_data)

# Create a separate dataset for fruit set data
Data_fruit = PL_data[!(PL_data$Measurement == "seed set"),]
Data_fruit = droplevels(Data_fruit)
dim(Data_fruit) # 145 obs with 17 columns

# Create a separate dataset for seed set data
Data_seed = PL_data[!(PL_data$Measurement == "fruit set"),]
Data_seed = droplevels(Data_seed)
dim(Data_seed) # 50 obs with 17 columns

# Assign variables for Data_fruit ----
class(Data_fruit$Study)
Data_fruit$Study<-factor(Data_fruit$Study)
length(levels(Data_fruit$Study)) #24 studies.

class(Data_fruit$Experiment)
Data_fruit$Experiment<-factor(Data_fruit$Experiment) # to change it's class as factor.
length(levels(Data_fruit$Experiment)) # 34 experiment.

class(Data_fruit$Continent)
Data_fruit$Continent<-factor(Data_fruit$Continent)
levels(Data_fruit$Continent) # "Asia","Europe","North America","South America" 4 levels

class(Data_fruit$Country)
Data_fruit$Country<-factor(Data_fruit$Country)
levels(Data_fruit$Country) # Australia, Belgium, Canada, China, India, Netherland
# Spain, Sweden, UK, USA, Yogoslavia 11 levels

class(Data_fruit$Climate)
Data_fruit$Climate<-factor(Data_fruit$Climate)
levels(Data_fruit$Climate) # Arid, Continental, Temperature 3 levels

class(Data_fruit$Cultivar)
Data_fruit$Cultivar<-factor(Data_fruit$Cultivar)
levels(Data_fruit$Cultivar) # 24 cultivars
table(Data_fruit$Cultivar)
sum(table(Data_fruit$Cultivar)) 

class(Data_fruit$Compatibility)
Data_fruit$Compatibility<-factor(Data_fruit$Compatibility)
levels(Data_fruit$Compatibility) # Self-incompatible, Semi-compatible 2 levels

class(Data_fruit$Pollen_type)
Data_fruit$Pollen_type<- factor(Data_fruit$Pollen_type)
levels(Data_fruit$Pollen_type) # Mixed, Single 2 levels

class(Data_fruit$Manipulation_level)
Data_fruit$Manipulation_level<-factor(Data_fruit$Manipulation_level)
levels(Data_fruit$Manipulation_level) # Branch, Flowers, Inflorescence, Plant 4 levels

# Data exploration---------------------------------------------------------------------------------------
# frequency distribution plot by histogram
# fruit set data
par(mfrow = c(1,2),mar=c(4.5,4.5,4.5,2.5))
freq_lnR_fruit = hist(Data_fruit$lnR,
                      main="Fruit set",
                      xlab=expression(ln(R[PD])),
                      xlim=c(-2.5,1.5),
                      breaks = 20,lwd=1.5,
                      cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)
axis(1, at=seq(-2.5,1.5, by=1), labels=seq(-2.5,1.5, by=1),cex.axis = 1.5)
segments(x0 = -2.5,x1= 1.5,y0 = 0,y1 = 0)
text(1,20,"(A)",cex = 1.5)

# seed set data
freq_lnR_seed = hist(Data_seed$lnR,
                     main="Seed set",
                     xlab=expression(ln(R[PD])),
                     xlim=c(-2,0.5),
                     breaks = 20,
                     cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)
segments(x0 = -2,x1= 0.5,y0 = 0,y1 = 0)
text(0.45,10,"(B)",cex = 1.5)

# explore patterns to know about data-----------------------------------------------------------------------#
# exploration for effect of continent on cultivars 
dim(Data_fruit)
colours = c("red","orange","green","blue")
par(mfrow = c(1,1),mar=c(3,8,3,2))
with(Data_fruit, plot(lnR,as.numeric(Cultivar),pch=16,col=colours[as.numeric(Data_fruit$Continent)],
                      yaxt = "n", ylab = ""))
cultivar.names = unique(Data_fruit$Cultivar)

axis(side=2,at=1:24,labels=sort(cultivar.names),hadj=1,las=1)
for (i in c(1:24)){
  lines(c(-2.2,1.2),c(i,i),lty=2,col="grey",lwd=1)
}
legend(-2.1,23,c("Asia","Europe","N America","S America"),pch=16,col=colours,bty="y")

# exploration for effect of manipulation level on cultivars
# draw a figure with for each variety the experiments coded by colours
#  and manipulation level by symbols

Data_fruit_mani = Data_fruit[!is.na(Data_fruit$Manipulation_level),]
dim(Data_fruit_mani) # 113 
Data_fruit_mani = droplevels(Data_fruit_mani)
Data_fruit_mani$Manipulation_level
Data_fruit_mani$Cultivar
cul.names = unique(Data_fruit_mani$Cultivar)

plot.symbols = 15:18
with(Data_fruit_mani, plot(lnR,as.numeric(Cultivar), col=Data_fruit_mani$Experiment, 
                           pch=plot.symbols[as.numeric(Data_fruit_mani$Manipulation_level)],
                           yaxt = "n", ylab = "",cex=1.2))
axis(side=2,at=1:19,labels=sort(cul.names),hadj=1,las=1)
for (i in c(1:19)){
  lines(c(-2.2,1.2),c(i,i),lty=3,col="grey",lwd=0.5)
}
legend(-2.2,19,c("Branch","Flower","Inflorescence","Plant"),cex=1.2,pch=15:18,col="black",bty="y")
# there are only two cultivar tested in more than one manipulation level,
# McIntosh tested in three manipulation level,
# Gala tested in both branch and infloresence level, 
# thus, this database could not reflect any interaction or confounding effect 
# of cultivar and manipulation level.

# exploration for effect of compatibility on cultivars
Data_fruit_compa = Data_fruit[!is.na(Data_fruit$Compatibility),]
dim(Data_fruit_compa) # 77 records
Data_fruit_compa = droplevels(Data_fruit_compa)
levels(Data_fruit_compa$Compatibility)
compa.cul.names = unique(Data_fruit_compa$Cultivar)

plot.symbols = 15:16
with(Data_fruit_compa, plot(lnR,as.numeric(Cultivar), col=Data_fruit_compa$Experiment, 
                            pch=plot.symbols[as.numeric(Data_fruit_compa$Compatibility)],
                            yaxt = "n", ylab = "",cex=1.2))
axis(side=2,at=1:11,labels=sort(compa.cul.names),hadj=1,las=1)
for (i in c(1:11)){
  lines(c(-2.2,1.2),c(i,i),lty=3,col="grey",lwd=0.5)
}
legend(-2.2,10,c("Self-incompatible","Semi-compatible"),cex=1.2,pch=15:16,col="black",bty="y")
# Here varieties (compatibility is specified) classified into two groups
# there is no self-compatible cultivar. 
# However, there is only one variety is semi-compatible, 
# I doubt this database could not support for comparison between compatibility type.

#### State for model selection ####
# To choose best random effect component to explain between-study variability 
# compare AIC (or BIC) of models fitted by "ML", and refit final model by REML
# start with estimation of means effect size,  
# then try to explain more variability of data according to add meaningful moderators


# 1.overall effect of pollination deficit----
overall_fruit0<-gls(lnR ~ 1,data = Data_fruit,method = "ML")
overall_fruit1<-lme(lnR ~ 1, random = ~1|Study, data = Data_fruit,method = "ML")
overall_fruit2<-lme(lnR ~ 1, random = ~1|Experiment, data = Data_fruit,method = "ML")
overall_fruit3<-lme(lnR ~ 1, random = ~1|Study/Experiment, data = Data_fruit,method = "ML")
anova(overall_fruit0,overall_fruit1,overall_fruit2,overall_fruit3)
# the best model is overall_fruit2, experiment as random effect component.

# Refit model by REML
overall_fruit <- lme(lnR ~ 1, random = ~1|Experiment, data = Data_fruit)
summary(overall_fruit)

intervals(overall_fruit)$fixed
# Fixed effects:
#                 lower       est.      upper
# (Intercept) -0.72881 -0.5267442 -0.3246784

lnR_FRUIT = intervals(overall_fruit)$fixed
PD = 1- exp(lnR_FRUIT)
PD
#               lower      est.     upper
# (Intercept) 0.5175172 0.4094755 0.2772403

# the mean effect size is -0.53(CI:-0.73,-0.32) for fruit set,
# that is, apple fruit set in natural condition is 41%(CI:28%, 52%)
# lower than artificial pollination.

# Test for efect of fruit set measurement-----
overall_fruit_measurement <- lme(lnR ~ Fruit_set_measurement, random = ~1|Experiment, data = Data_fruit)
summary(overall_fruit_measurement)
anova(overall_fruit_measurement) # p = 0.128 > 0.05, 
# indicating fruit set measurement is no effect, 
# lnR is not different between initial fruit,final fruit and unknown fruit set

# excluded records with Not state in Fruit set measurement
Data_fruit1 = Data_fruit[!(Data_fruit$Fruit_set_measurement == "Not state"),]
dim(Data_fruit1) # 129 records
overall_fruit_measurement1 <- lme(lnR ~ Fruit_set_measurement, random = ~1|Experiment, data = Data_fruit1)
summary(overall_fruit_measurement1) 
anova(overall_fruit_measurement1) # p = 0.3117
# same as previous result, lnR is not different between intial fruit set and final fruit set
# whatever using whole data or subset data 

# 2. what extent the pollination deficit across continents?-----
mcontinent0<-gls(lnR ~ Continent, data = Data_fruit, method = "ML")
mcontinent1<-lme(lnR ~ Continent, random = ~1|Study,data = Data_fruit, method = "ML")
mcontinent2<-lme(lnR ~ Continent, random = ~1|Experiment,data = Data_fruit, method = "ML")
mcontinent3<-lme(lnR ~ Continent, random = ~1|Study/Experiment,data = Data_fruit, method = "ML")
anova(mcontinent0,mcontinent1,mcontinent2,mcontinent3)
# the mcontinent2 is the best, the model with experiment as random effect component 

# Refit model by REML #
mcontinent<-lme(lnR ~ Continent, random = ~1|Experiment,data = Data_fruit)
summary(mcontinent)
anova(mcontinent) # p=0.0068, significnat difference between continents

# make multiple comparison #
mod_continent<-lsmeans(mcontinent, ~ Continent) # results same as model summary 
pairs(mod_continent,adjust = "tukey")
multcomp::cld(mod_continent,Letters = c(letters))

# model validation #
hist(resid(mcontinent),freq = F) 
qqPlot(resid(mcontinent))
plot(fitted(mcontinent),resid(mcontinent))
# residuals seems ok, close to normal

# 3. what extent the pollination deficit among cultivars?----
mcultivar0<-gls(lnR ~ Cultivar,data = Data_fruit,method = "ML")
mcultivar1<-lme(lnR ~ Cultivar, random = ~1|Study,data = Data_fruit,method = "ML")
mcultivar2<-lme(lnR ~ Cultivar, random = ~1|Experiment,data = Data_fruit,method = "ML")
mcultivar3<-lme(lnR ~ Cultivar, random = ~1|Study/Experiment,data = Data_fruit,method = "ML")
anova(mcultivar0,mcultivar1,mcultivar2,mcultivar3)
# model mcultivar2 is the best, the model with experiment as random effect conponent

# Refit model by REML
mcultivar<-lme(lnR ~ Cultivar, random = ~1|Experiment,data = Data_fruit)
summary(mcultivar)
anova(mcultivar)# p = 0.0063, significant difference between cultivars

# make multiple comparison #
mod_cul<-lsmeans(mcultivar, ~ Cultivar)
multcomp::cld(mod_cul,Letters = c(letters))
# letters are all same, here, indicating all cultivars are in same homogeneity group
# because anova said there is signficant difference between cultivars,so it is strange
# tukey is available to test difference between treatmens 
# when the number of treatment levels more than 10
# We doubt this is because of no enough data points under each cultivar,
# so compare means with different reference cultivar mean 

# to check any difference between control and other cultivars
contrast(mod_cul,method = "trt.vs.ctrl",ref = 8)
# contrast                 estimate    SE df t.ratio p.value
# Fuji - Eva               -1.43222 0.363 25 -3.945  0.0098
# Princesa - Eva           -0.73532 0.213 96 -3.454  0.0153 

contrast(mod_cul,method = "trt.vs.ctrl",ref = 9)
# Cox - Fuji                  1.0942 0.324 25 3.380   0.0368 
# Eva - Fuji                  1.4322 0.363 25 3.945   0.0098
# Red Delicious - Fuji        1.4125 0.395 25 3.577   0.0234  

# calculate all the differences between each cultivar and general mean
contrast(mod_cul,method = "eff",adjust = "sidak")
# there is no significant effect of any cultivars

# conclude: lnR for Fuji is significant lower than Eva and Red Delicious,
# lnR for Princesa is significantly lower than Eva.

# This results are not so much robust because of p-value should be corrected,
# hence, show effect sizes of all cultivars.
mcultivar_effect<-lme(lnR ~ Cultivar-1, random = ~1|Experiment,data = Data_fruit)
summary(mcultivar_effect)

# extract mean,confidence intervals for each cultivar
mean_cultivar = intervals(mcultivar_effect)$fixed[,2]
upper_cultivar = intervals(mcultivar_effect)$fixed[,3]
lower_cultivar = intervals(mcultivar_effect)$fixed[,1]

# reshape data set
cultivar.names = unique(Data_fruit$Cultivar)
Data_cultivar = data.frame(cultr = sort(cultivar.names),mean_cultr = mean_cultivar,
                           lower_CI = lower_cultivar, upper_CI = upper_cultivar)
names(Data_cultivar)

#  plot #
DATA_CULTIVAR = Data_cultivar[order(Data_cultivar$mean_cultr),]
par(mfrow = c(1,1),mar=c(5,5,5,2.5))
dotchart(DATA_CULTIVAR$mean_cultr,labels = DATA_CULTIVAR$cultr,
         xlim = c(-2.1,2),cex = 1.3,xlab = expression(ln(R[PD])),pch = 16,main = "Fruit set")
for (i in 1:nrow(DATA_CULTIVAR)) {
  lines(x=c(DATA_CULTIVAR[i,3],DATA_CULTIVAR[i,4]), y=c(i,i))
}
abline(v=0,lty = "dashed",col = "grey19")

# model validation #
hist(resid(mcultivar),freq = F) 
qqPlot(resid(mcultivar))
plot(fitted(mcultivar),resid(mcultivar))
# residuals seems ok, close to normal

#### Analysis in other factor(s) as explanatory variable ####
# 4. Does the pollination deficit is differ between countries?-----
mcountry_fruit0 <- gls(lnR ~ Country,data = Data_fruit,method = "ML") 
mcountry_fruit1 <- lme(lnR ~ Country, random = ~1|Study, data = Data_fruit,method = "ML") 
mcountry_fruit2 <- lme(lnR ~ Country, random = ~1|Experiment, data = Data_fruit,method = "ML")
mcountry_fruit3 <- lme(lnR ~ Country, random = ~1|Study/Experiment, data = Data_fruit,method = "ML")
anova(mcountry_fruit0,mcountry_fruit1,mcountry_fruit2,mcountry_fruit3)
# the model mcountry_fruit2 is the best, the model with experiment as radom effect component

# Refit model  by REML
mcountry_fruit <- lme(lnR ~ Country, random = ~1|Experiment, data = Data_fruit) 
summary(mcountry_fruit)
anova(mcountry_fruit) # p = 0.2268, no difference between countries

# model validation #
hist(resid(mcountry_fruit),freq = F) 
qqPlot(resid(mcountry_fruit))
plot(fitted(mcountry_fruit),resid(mcountry_fruit))
# residuals seems ok, close to normal

# 5. Does the pollination deficit is associated with climate?-----
mclimate_fruit0 <- gls(lnR ~ Climate, data = Data_fruit,method = "ML") 
mclimate_fruit1 <- lme(lnR ~ Climate, random = ~1|Study, data = Data_fruit,method = "ML") 
mclimate_fruit2 <- lme(lnR ~ Climate, random = ~1|Experiment, data = Data_fruit,method = "ML") 
mclimate_fruit3 <- lme(lnR ~ Climate, random = ~1|Study/Experiment, data = Data_fruit,method = "ML") 
anova(mclimate_fruit0,mclimate_fruit1,mclimate_fruit2,mclimate_fruit3)
# model mclimate_fruit2 is the best, the model with experiment as random effect component

# Refit model by REML
mclimate_fruit <- lme(lnR ~ Climate, random = ~1|Experiment, data = Data_fruit) 
summary(mclimate_fruit)
anova(mclimate_fruit) # p = 0.9028, no difference between climate zones

# model validation #
hist(resid(mclimate_fruit),freq = F) # residuals seems ok, close to normal.
qqPlot(resid(mclimate_fruit))
plot(fitted(mclimate_fruit),resid(mclimate_fruit))
# residuals looks good

# 6. Does the pollination deficit is associated with compatibility?-----
mcompa0<-gls(lnR ~ Compatibility,data = Data_fruit,na.action = na.exclude,method = "ML")
mcompa1<-lme(lnR ~ Compatibility,random = ~1|Study,data = Data_fruit,na.action = na.exclude,method = "ML")
mcompa2<-lme(lnR ~ Compatibility,random = ~1|Experiment,data = Data_fruit,na.action = na.exclude,method = "ML")
mcompa3<-lme(lnR ~ Compatibility,random = ~1|Study/Experiment,data = Data_fruit,na.action = na.exclude,method = "ML")
anova(mcompa0,mcompa1,mcompa2,mcompa3)
# model mcompa2 is the best, the model with experiment as random effect

# Refit model by REML
mcompa<-lme(lnR ~ Compatibility,random = ~1|Experiment,data = Data_fruit,na.action = na.exclude)
summary(mcompa)
anova(mcompa) # p = 0.6723,there is no difference between compatibilities

# model validation #
hist(resid(mcompa),freq = F) 
qqPlot(resid(mcompa))
plot(fitted(mcompa),resid(mcompa))
# residual looks good 

# 7. Does the pollination deficit is associated with manipulation level?-----
manipulation_fruit0 <- gls(lnR ~ Manipulation_level, data = Data_fruit,na.action = na.exclude,method = "ML") 
manipulation_fruit1 <- lme(lnR ~ Manipulation_level, random = ~1|Study, data = Data_fruit,na.action = na.exclude,method = "ML") 
manipulation_fruit2 <- lme(lnR ~ Manipulation_level, random = ~1|Experiment, data = Data_fruit,na.action = na.exclude,method = "ML") 
manipulation_fruit3 <- lme(lnR ~ Manipulation_level, random = ~1|Study/Experiment, data = Data_fruit,na.action = na.exclude,method = "ML") 
anova(manipulation_fruit0,manipulation_fruit1,manipulation_fruit2,manipulation_fruit3)
# model manipulation_fruit2 is the best, the model with experiment as random effect component

# Refit model by REML
manipulation_fruit <- lme(lnR ~ Manipulation_level, random = ~1|Experiment, data = Data_fruit,na.action = na.exclude) 
summary(manipulation_fruit)
anova(manipulation_fruit) # p = 0.7591, no difference between difference manipulation level

# model validation #
hist(resid(manipulation_fruit),freq = F) 
qqPlot(resid(manipulation_fruit))
plot(fitted(manipulation_fruit),resid(manipulation_fruit))
# residual looks good 

# 8. Does the pollination deficit is associated with pollen type?-----
mpollen_fruit0 <- gls(lnR ~ Pollen_type,data = Data_fruit,na.action = na.exclude,method = "ML") 
mpollen_fruit1 <- lme(lnR ~ Pollen_type, random = ~1|Study, data = Data_fruit,na.action = na.exclude,method = "ML") 
mpollen_fruit2 <- lme(lnR ~ Pollen_type, random = ~1|Experiment, data = Data_fruit,na.action = na.exclude,method = "ML") 
mpollen_fruit3 <- lme(lnR ~ Pollen_type, random = ~1|Study/Experiment, data = Data_fruit,na.action = na.exclude,method = "ML") 
anova(mpollen_fruit0,mpollen_fruit1,mpollen_fruit2,mpollen_fruit3)
# model mpollen_fruit2 is the best, the model with experiment as random effect component

# Refit model by REML
mpollen_fruit <- lme(lnR ~ Pollen_type, random = ~1|Experiment, data = Data_fruit,na.action = na.exclude) 
summary(mpollen_fruit)
anova(mpollen_fruit) # p = 0.9103, no difference between pollen types

# model validation #
hist(resid(mpollen_fruit),freq = F) 
qqPlot(resid(mpollen_fruit))
plot(fitted(mpollen_fruit),resid(mpollen_fruit))
# residual looks good

#### Analysis in seed set data #########################################################################
# Assign variables for Data_seed ----
class(Data_seed$Study)
Data_seed$Study<-factor(Data_seed$Study) # to change it's class as factor
length(levels(Data_seed$Study)) # 12 studies

class(Data_seed$Experiment)
Data_seed$Experiment<-factor(Data_seed$Experiment) 
length(levels(Data_seed$Experiment)) # 21 experiments

class(Data_seed$Continent)
Data_seed$Continent<-factor(Data_seed$Continent)
levels(Data_seed$Continent) # Asia, Europe, North America, South America, Oceania 5 levels

class(Data_seed$Cultivar)
Data_seed$Cultivar<-factor(Data_seed$Cultivar)
levels(Data_seed$Cultivar) # 18 cultivars

class(Data_seed$Country)
Data_seed$Country<-factor(Data_seed$Country)
levels(Data_seed$Country) # Australia, Belgium, China, India, Netherlands
                          # Sweden, UK, USA  9 countries

class(Data_seed$Climate)
Data_seed$Climate<-factor(Data_seed$Climate)
levels(Data_seed$Climate) # Arid, Continental, Temperature 3 levels

class(Data_seed$Pollen_type)
Data_seed$Pollen_type<- factor(Data_seed$Pollen_type)
levels(Data_seed$Pollen_type) # Mixed, Single 2 levels

# 1. Overall pollination deficit on seed set----
overall_seed0<-gls(lnR ~ 1,data = Data_seed,method = "ML")
overall_seed1<-lme(lnR ~ 1, random = ~1|Study, data = Data_seed,method = "ML")
overall_seed2<-lme(lnR ~ 1, random = ~1|Experiment, data = Data_seed,method = "ML")
overall_seed3<-lme(lnR ~ 1, random = ~1|Study/Experiment, data = Data_seed,method = "ML")
anova(overall_seed0,overall_seed1,overall_seed2,overall_seed3)
# model overall_seed2 is the best, the model with experiment as random effect

# Refit model by REML
overall_seed<-lme(lnR ~ 1,random = ~1| Experiment,data = Data_seed)
summary(overall_seed)
# lnR = -0.228, PD = 1-exp(-0.228) = 0.20, 
# means seed number would be set 20% fewer than hand pollinaiton

intervals(overall_seed)
# Fixed effects:
#   lower       est.       upper
# (Intercept) -0.3858704 -0.2281743 -0.07047824

lnR_SEED = intervals(overall_seed)$fixed
PD_SEED = 1-exp(lnR_SEED)
PD_SEED
# 
#               lower      est.      upper
# (Intercept) 0.3201414 0.2040145 0.06805198

# model validation #
hist(resid(overall_seed),freq = F) 
qqPlot(resid(overall_seed))
plot(fitted(overall_seed),resid(overall_seed))
# residuals seems ok, close to normal 

# 2. what extent the pollination deficit across continent?----
mcontinent_seed0<-gls(lnR ~ Continent,data = Data_seed,method = "ML")
mcontinent_seed1<-lme(lnR ~ Continent,random = ~1|Study,data = Data_seed,method = "ML")
mcontinent_seed2<-lme(lnR ~ Continent,random = ~1|Experiment,data = Data_seed,method = "ML")
mcontinent_seed3<-lme(lnR ~ Continent,random = ~1|Study/Experiment,data = Data_seed,method = "ML")
anova(mcontinent_seed0,mcontinent_seed1,mcontinent_seed2,mcontinent_seed3)
# model mcontinent_seed2 is the best, the model with experiment as random effect 

# Refit model by REML 
mcontinent_seed<-lme(lnR ~ Continent,random = ~1|Experiment,data = Data_seed)
summary(mcontinent_seed)
anova(mcontinent_seed)  # p = 0.1456 > 0.05, no difference between continents

# to see estimate for each continent
seed_conti<-lsmeans(mcontinent_seed, ~ Continent) # results same as model summary 
multcomp::cld(seed_conti,Letters = c(letters))

# model validation #
hist(resid(mcontinent_seed),freq = F) 
qqPlot(resid(mcontinent_seed))
plot(fitted(mcontinent_seed),resid(mcontinent_seed))
# residuals seems ok, close to normal

# 3. what extent the pollination deficit among cultivars?----
mcultivar_seed0<-gls(lnR ~ Cultivar,data = Data_seed,method = "ML")
mcultivar_seed1<-lme(lnR ~ Cultivar,random = ~1|Study,data = Data_seed,method = "ML")
mcultivar_seed2<-lme(lnR ~ Cultivar,random = ~1|Experiment,data = Data_seed,method = "ML")
mcultivar_seed3<-lme(lnR ~ Cultivar,random = ~1|Study/Experiment,data = Data_seed,method = "ML")
anova(mcultivar_seed0,mcultivar_seed1,mcultivar_seed2,mcultivar_seed3)
# model mcultivar_seed2 is the best, the model with experiment as random effect component

# Refit model by REML
mcultivar_seed<-lme(lnR ~ Cultivar,random = ~1|Experiment,data = Data_seed)
summary(mcultivar_seed)
anova(mcultivar_seed) # p = 0.4315>0.05, no difference between cultivars

# extract mean,confidence intervals for each cultivar
mcultivar_seed_plot<-lme(lnR ~ Cultivar-1,random = ~1|Experiment,data = Data_seed)
mean_cultivar_seed = intervals(mcultivar_seed_plot)$fixed[,2]
upper_cultivar_seed = intervals(mcultivar_seed_plot)$fixed[,3]
lower_cultivar_seed = intervals(mcultivar_seed_plot)$fixed[,1]

# reshape data set
cultivar.names_seed = unique(Data_seed$Cultivar)
Data_cultivar_seed = data.frame(cultr_seed = sort(cultivar.names_seed),mean_cultr_seed = mean_cultivar_seed,
                                lower_CI_seed = lower_cultivar_seed, upper_CI_seed = upper_cultivar_seed)
names(Data_cultivar_seed)

#  plot #
Cultivar_seed = Data_cultivar_seed[order(Data_cultivar_seed$mean_cultr_seed),]
par(mfrow = c(1,1),mar=c(5,5,5,2.5))
dotchart(Cultivar_seed$mean_cultr_seed,labels = Cultivar_seed$cultr_seed,
         xlim = c(-1.5,1.5),cex = 1.3,xlab = expression(ln(R[PD])),pch = 16,main = "Seed set")
for (i in 1:nrow(Cultivar_seed)) {
  lines(x=c(Cultivar_seed[i,3],Cultivar_seed[i,4]), y=c(i,i))
}
abline(v=0,lty = "dashed",col = "grey19")

# model validation #
hist(resid(mcultivar_seed),freq = F) 
qqPlot(resid(mcultivar_seed))
plot(fitted(mcultivar_seed),resid(mcultivar_seed))
# residuals seems ok 

# analysis in other factor(s) as explanatory variable #
# 4. does the degree of pollination deficit differ between countries?----
mcountry_seed0 <- gls(lnR ~ Country,data = Data_seed,method = "ML") 
mcountry_seed1 <- lme(lnR ~ Country, random = ~1|Study, data = Data_seed,method = "ML") 
mcountry_seed2 <- lme(lnR ~ Country, random = ~1|Experiment, data = Data_seed,method = "ML") 
mcountry_seed3 <- lme(lnR ~ Country, random = ~1|Study/Experiment, data = Data_seed,method = "ML") 
anova(mcountry_seed0,mcountry_seed1,mcountry_seed2,mcountry_seed3)
# model mcountry_seed0 and mcountry_seed2 is similar,check results of random effect model 

# Refit model by REML
mcountry_seed <- lme(lnR ~ Country, random = ~1|Experiment, data = Data_seed) 
summary(mcountry_seed)
# choose random effect model since there is evidence/contribution of random effect 
anova(mcountry_seed) # p = 0.2741 > 0.05, no difference between countries

# to see estimate for each country
seed_country<-lsmeans(mcountry_seed, ~ Country)

# model validation #
hist(resid(mcountry_seed),freq = F) 
qqPlot(resid(mcountry_seed))
plot(fitted(mcountry_seed),resid(mcountry_seed))
# residuals seems ok 

# 5. does the degree of pollination deficit differ across climate zones?----
mclimate_seed0 <- gls(lnR ~ Climate,data = Data_seed,method = "ML") 
mclimate_seed1 <- lme(lnR ~ Climate, random = ~1|Study, data = Data_seed,method = "ML") 
mclimate_seed2 <- lme(lnR ~ Climate, random = ~1|Experiment, data = Data_seed,method = "ML") 
mclimate_seed3 <- lme(lnR ~ Climate, random = ~1|Study/Experiment, data = Data_seed,method = "ML") 
anova(mclimate_seed0,mclimate_seed1,mclimate_seed2,mclimate_seed3)
# model mclimate_seed2 is the best, the model with experiment as random effect component

# Refit modle by REML
mclimate_seed <- lme(lnR ~ Climate, random = ~1|Experiment, data = Data_seed) 
summary(mclimate_seed)
anova(mclimate_seed) # p = 0.2099 > 0.05, no difference between climate zones

# model validation #
hist(resid(mclimate_seed),freq = F) 
qqPlot(resid(mclimate_seed))
plot(fitted(mclimate_seed),resid(mclimate_seed))
# residuals seems ok 

# 6. does the degree of pollination deficit differ between compatibility?----
mcompatibility_seed0 <- gls(lnR ~ Compatibility, data = Data_seed,na.action = na.exclude,method = "ML") 
mcompatibility_seed1 <- lme(lnR ~ Compatibility, random = ~1|Study, data = Data_seed,na.action = na.exclude,method = "ML") 
mcompatibility_seed2 <- lme(lnR ~ Compatibility, random = ~1|Experiment, data = Data_seed,na.action = na.exclude,method = "ML") 
mcompatibility_seed3 <- lme(lnR ~ Compatibility, random = ~1|Study/Experiment, data = Data_seed,na.action = na.exclude,method = "ML") 
anova(mcompatibility_seed0,mcompatibility_seed1,mcompatibility_seed2,mcompatibility_seed3)
# no difference between fixed effect model and random effect model 

# Refit random effect model to check results
mcompatibility_seed <- lme(lnR ~ Compatibility, random = ~1|Experiment, data = Data_seed,na.action = na.exclude) 
summary(mcompatibility_seed) # there is evidence of between-study(experiment) variability, so choose random effect model 
anova(mcompatibility_seed) # p = 0.4260 > 0.05, no difference between different compatibility

# model validation #
hist(resid(mcompatibility_seed),freq = F) 
qqPlot(resid(mcompatibility_seed))
plot(fitted(mcompatibility_seed),resid(mcompatibility_seed))
# residuals seems ok 

# 7. does the degree of pollination deficit differ across manipulation level ?----
manipulation_seed0 <- gls(lnR ~ Manipulation_level,data = Data_seed,na.action = na.exclude,method = "ML") 
manipulation_seed1 <- lme(lnR ~ Manipulation_level, random = ~1|Study, data = Data_seed,na.action = na.exclude,method = "ML") 
manipulation_seed2 <- lme(lnR ~ Manipulation_level, random = ~1|Experiment, data = Data_seed,na.action = na.exclude,method = "ML") 
manipulation_seed3 <- lme(lnR ~ Manipulation_level, random = ~1|Study/Experiment, data = Data_seed,na.action = na.exclude,method = "ML") 
anova(manipulation_seed0,manipulation_seed1,manipulation_seed2,manipulation_seed3)
# model manipulation_seed2 is the best, the model with experiment as random effect component

# Refit model by REML
manipulation_seed <- lme(lnR ~ Manipulation_level, random = ~1|Experiment, data = Data_seed,na.action = na.exclude) 
summary(manipulation_seed)
anova(manipulation_seed) # p = 0.92 > 0.05, no difference between difference manipulation level

# model validation #
hist(resid(manipulation_seed),freq = F) 
qqPlot(resid(manipulation_seed))
plot(fitted(manipulation_seed),resid(manipulation_seed))
# residuals seems ok 

# 8. does the degree of pollination deficit differ across pollne type ?----
mpollen_seed0 <- gls(lnR ~ Pollen_type, data = Data_seed,na.action = na.exclude,method = "ML") 
mpollen_seed1 <- lme(lnR ~ Pollen_type, random = ~1|Study, data = Data_seed,na.action = na.exclude,method = "ML") 
mpollen_seed2 <- lme(lnR ~ Pollen_type, random = ~1|Experiment, data = Data_seed,na.action = na.exclude,method = "ML") 
mpollen_seed3 <- lme(lnR ~ Pollen_type, random = ~1|Study/Experiment, data = Data_seed,na.action = na.exclude,method = "ML") 
anova(mpollen_seed0,mpollen_seed1,mpollen_seed2,mpollen_seed3)
# model mpollen_seed2 is the best, the model with experiment as random effect component

# Refit model by REML
mpollen_seed <- lme(lnR ~ Pollen_type, random = ~1|Experiment, data = Data_seed,na.action = na.exclude) 
summary(mpollen_seed)
anova(mpollen_seed) # p = 0.318 > 0.05, no difference between pollen type

# model validation #
hist(resid(mpollen_seed),freq = F) 
qqPlot(resid(mpollen_seed))
plot(fitted(mpollen_seed),resid(mpollen_seed))
# residual seems ok

# Test difference between fruit set and seed set on lnR----
Overall_fruit_seed = lme(lnR ~ Measurement, random = ~1|Experiment,data = PL_data,method = "REML")
anova(Overall_fruit_seed)
# p = 0.016 < 0.05, there is significant difference of lnR between fruit set and seed set  

#### Plotting ###############
# plot overall pollination deficit for fruit set and seed set----
# overall effect for fruit set: 
#                 lower       est.      upper
# (Intercept)   -0.72881  -0.5267442  -0.3246784


# overall effect for seed set: 
#                 lower       est.       upper
# (Intercept) -0.3858704 -0.2281743 -0.07047824

overall_fruit_seed = data.frame(response = 
                                  c("fruit set","seed set"),
                                estimate = c(-0.526,-0.228),
                                Lower = c(-0.728,-0.385),
                                Upper = c(-0.324,-0.070))

par(mar=c(4.5,4.5,4.5,2.5))
overall_PD <- ggplot(PL_data, aes(x = Measurement, y = lnR))+scale_y_continuous(limits = c(-2.2,1.2))
overall_PD <- overall_PD + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "",y=expression(lnR[PD]), x = "",
       axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) 
overall_PD<-overall_PD + geom_pointrange(overall_fruit_seed,
                                         mapping = aes(x=response, y=estimate, ymin=Lower, ymax=Upper), 
                                         size=1.3,color = c("springgreen4","blue4"),position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
overall_PD


# Plot - Continent ---------------------------------------------------------------------------------
mod_continent # estimate from model for fruit set data
seed_conti  # estimate from model for seed set data

table(Data_fruit$Continent)
table(Data_seed$Continent)

# create dataset for showing lnR across continent
Continent_fruit = data.frame(continent = 
                               c("Asia","Europe","North America","South America"),
                             estimate = c(-0.995,-0.363,-0.346,-0.111),
                             Lower = c(-1.303,-0.627,-0.904,-0.592),
                             Upper = c(-0.687,-0.1,0.213,0.369))

Continent_seed = data.frame(continent = 
                              c("Asia","South America","North America","Europe"),
                            estimate = c(-0.6309,-0.5086,-0.2283,-0.1301),
                            Lower = c(-1.082,-0.886,-1.058,-0.310),
                            Upper = c(-0.1802,-0.1308,0.6019,0.0494))

# Plot #
# fruit set data for continent
PD_fruit_continent <- ggplot(Data_fruit, aes(x = Continent, y = lnR))+scale_y_continuous(limits = c(-2.2,1.2)) # original data points
PD_fruit_continent <- PD_fruit_continent + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") + # make points spread around
  theme_bw() + theme(panel.grid = element_blank()) + # make backgroud white 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),# locate the position of legend
        legend.background=element_blank(),
        legend.text = element_text(size = 14), # make bacground of legend transparent
        axis.text=element_text(size=14, color = "black")) + # adjust the axis text 
  labs(title = "fruit set",y=expression(ln(R[PD])), x = "",
       axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) 
PD_fruit_continent<-PD_fruit_continent + geom_pointrange(Continent_fruit, # add the estimate mean value and CI
                                                         mapping = aes(x=continent, y=estimate, ymin=Lower, ymax=Upper), 
                                                         size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
fruit<-PD_fruit_continent

# seed set data for continent
Data_removeOC_seed <- Data_seed[!(Data_seed$Continent == "Oceania"),]
Data_removeOC_seed$Continent = factor(Data_removeOC_seed$Continent,levels = 
                                        c("Asia","Europe","North America","South America"))
PD_seed_continent <- ggplot(Data_removeOC_seed, aes(x = Continent, y = lnR))+ 
  scale_y_continuous(limits = c(-2.2,1.2),position = "right")
PD_seed_continent <- PD_seed_continent + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "seed set",y=expression(ln(R[PD])), 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PD_seed_continent<-PD_seed_continent + 
  geom_pointrange(Continent_seed,
                  mapping = aes(x=continent, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
seed<-PD_seed_continent

# combine plots 
require(gridExtra) 
grid.arrange(fruit, seed, ncol=2)


#### Plot - country ####
# fruit set data #
mod_country1 <- emmeans(mcountry_fruit, "Country") # results same with those from model summary 
multcomp::cld(mod_country1,Letters = c(letters))
# seed set data #
mod_country2 <- emmeans(mcountry_seed, "Country")
multcomp::cld(mod_country2,Letters = c(letters))

table(Data_fruit$Country)
table(Data_seed$Country)

# create dataset for showing lnR across continent
Country_fruit = data.frame(country= 
                             c("India","China",
                               "Netherlands","UK","Belgium","Spain","Sweden","Yogoslavia",
                               "Canada","USA","Argentina"),
                           estimate = c(-1.171,-0.933,
                                        -0.564,-0.332,-0.427,-0.427,-0.175,0.147,
                                        -0.504,-0.210,-0.110),
                           Lower = c(-1.898,-1.320,
                                     -1.135,-0.831,-1.043,-1.704,-1.316,-0.894,
                                     -1.406,-1.029,-0.642),
                           Upper = c(-0.444,-0.545,0.007,0.167, 0.188,0.849,0.966,1.187,
                                     0.399,0.608,0.421))

Country_seed = data.frame(country = 
                            c("China",
                              "UK","Belgium","Sweden","Serbia","Netherlands",
                              "USA","Argentina","New Zealand"),
                          estimate = c(-0.6315,
                                       -0.4951,-0.2103,-0.0426,-0.0223,-0.0065,
                                       -0.2283,-0.5086,0.0889),
                          Lower = c(-1.101,
                                    -1.031,-0.548,-0.898,-0.419,-0.332,
                                    -1.084,-0.879,-0.766),
                          Upper = c(-0.1617,
                                    0.0406,0.1274,0.8129,0.3743,0.3192,
                                    0.6272,-0.1386,0.9444))

# Plot #
# fruit set data for country
Data_fruit$Country = factor(Data_fruit$Country,levels = c("India","China",
                                                          "Netherlands","UK","Belgium","Spain","Sweden","Yogoslavia",
                                                          "Canada","USA","Argentina"))
# plot #
par(mar=c(4.5,4.5,4.5,2.5))
PD_fruit_country <- ggplot(Data_fruit, aes(x = Country, y = lnR))+scale_y_continuous(limits = c(-2.2,1.3)) # original data points
PD_fruit_country <- PD_fruit_country + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") + # make points spread around
  theme_bw() + theme(panel.grid = element_blank()) + # make backgroud white 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),# locate the position of legend
        legend.background=element_blank(),
        legend.text = element_text(size = 14), # make bacground of legend transparent
        axis.text=element_text(size=14, color = "black")) + # adjust the axis text 
  labs(title = "Fruit set",y=expression(ln(R[PD])), x = "",
       axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(angle =- 50, vjust = 0.5))
PD_fruit_country<-PD_fruit_country + geom_pointrange(Country_fruit, # add the estimate mean value and CI
                                                     mapping = aes(x=country, y=estimate, ymin=Lower, ymax=Upper), 
                                                     size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

PD_fruit_country

# seed set data for country
Data_seed$Country = factor(Data_seed$Country,levels = c("China",
                                                        "UK","Belgium","Sweden","Serbia","Netherlands",
                                                        "USA","Argentina","New Zealand"))
# plot #
par(mar=c(4.5,4.5,4.5,2.5))
PD_seed_country <- ggplot(Data_seed, aes(x = Country, y = lnR))+ 
  scale_y_continuous(limits = c(-2.2,1.3))
PD_seed_country <- PD_seed_country + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "Seed set",y=expression(ln(R[PD])), 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(axis.text.x=element_text(angle =- 50, vjust = 0.5))
PD_seed_country<-PD_seed_country + 
  geom_pointrange(Country_seed,
                  mapping = aes(x=country, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
PD_seed_country


#### Plot - climate ####
# fruit set data #
mod_climate1 <- emmeans(mclimate_fruit, "Climate") # results are same as results from model summary
multcomp::cld(mod_climate1,Letters = c(letters))

# seed set data #
mod_climate2 <- emmeans(mclimate_seed, "Climate")
multcomp::cld(mod_climate2,Letters = c(letters))

# check data
table(Data_fruit$Climate)
table(Data_seed$Climate)

# create dataset for showing lnR across climate zone
Climate_fruit = data.frame(climate= 
                             c("Arid","Temperate","Continental"),
                           estimate = c(-0.679,-0.510,-0.523),
                           Lower = c(-1.403,-0.751,-1.093),
                           Upper = c(0.043,-0.268,0.046))

Climate_seed = data.frame(climate = 
                            c("Arid","Temperate","Continental"),
                          estimate = c(-0.696,-0.185,-0.483),
                          Lower = c(-1.287,-0.351,-1.338),
                          Upper = c( -0.1046,-0.0201,0.3723))

# Plot #
par(mar=c(4.5,4.5,4.5,2.5))
# fruit set data for climate
Data_fruit$Climate = factor(Data_fruit$Climate,levels = c("Arid","Temperate","Continental"))
PD_fruit_climate <- ggplot(Data_fruit, aes(x = Climate, y = lnR))+scale_y_continuous(limits = c(-2.2,1.3)) # original data points
PD_fruit_climate <- PD_fruit_climate + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") + # make points spread around
  theme_bw() + theme(panel.grid = element_blank()) + # make backgroud white 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),# locate the position of legend
        legend.background=element_blank(),
        legend.text = element_text(size = 14), # make bacground of legend transparent
        axis.text=element_text(size=14, color = "black")) + # adjust the axis text 
  labs(title = "fruit set",y=expression(ln(R[PD])), x = "",
       axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PD_fruit_climate<-PD_fruit_climate + geom_pointrange(Climate_fruit, # add the estimate mean value and CI
                                                     mapping = aes(x=climate, y=estimate, ymin=Lower, ymax=Upper), 
                                                     size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

fruit_climate<-PD_fruit_climate

# seed set data for climate
Data_seed$Climate = factor(Data_seed$Climate,levels = c("Arid","Temperate","Continental"))

PD_seed_climate <- ggplot(Data_seed, aes(x = Climate, y = lnR))+ 
  scale_y_continuous(limits = c(-2.2,1.3),position = "right")
PD_seed_climate <- PD_seed_climate + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "seed set",y=expression(ln(R[PD])), 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PD_seed_climate<-PD_seed_climate + 
  geom_pointrange(Climate_seed,
                  mapping = aes(x=climate, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
seed_climate<-PD_seed_climate

climate_plot<-grid.arrange(fruit_climate,seed_climate,ncol =2)


#### Plot - compatibility ####
# fruit set data #
mod_compa1 <- emmeans(mcompa, "Compatibility")
multcomp::cld(mod_compa1,Letters = c(letters))

# seed set data #
mod_compa2 <- emmeans(mcompatibility_seed, "Compatibility")
multcomp::cld(mod_compa2,Letters = c(letters))

# create dataset for showing lnR across compatibility
Compa_fruit = data.frame(compatibility= 
                           c("Self-incompatible","Semi-compatible"),
                         estimate = c(-0.601,-0.487),
                         Lower = c(-0.924,-1.073),
                         Upper = c(-0.278, 0.100))

Compa_seed = data.frame(compatibility = 
                          c("Self-incompatible","Semi-compatible"),
                        estimate = c(-0.399,-0.171),
                        Lower = c(-0.656,-0.782),
                        Upper = c(-0.143,0.440))

# Plot #
# fruit set data for compatibility
Data_compa_fruit = Data_fruit[!is.na(Data_fruit$Compatibility),] # 77 obs
Data_compa_fruit$Compatibility = factor(Data_compa_fruit$Compatibility,
                                        levels = c("Self-incompatible","Semi-compatible"))
PD_fruit_compa <- ggplot(Data_compa_fruit, aes(x = Compatibility, y = lnR))+scale_y_continuous(limits = c(-2.2,1.3)) # original data points
PD_fruit_compa <- PD_fruit_compa + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") + # make points spread around
  theme_bw() + theme(panel.grid = element_blank()) + # make backgroud white 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),# locate the position of legend
        legend.background=element_blank(),
        legend.text = element_text(size = 14), # make bacground of legend transparent
        axis.text=element_text(size=14, color = "black")) + # adjust the axis text 
  labs(title = "fruit set",y=expression(ln(R[PD])), x = "",
       axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PD_fruit_compa<-PD_fruit_compa + geom_pointrange(Compa_fruit, # add the estimate mean value and CI
                                                 mapping = aes(x=compatibility, y=estimate, ymin=Lower, ymax=Upper), 
                                                 size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

fruit_compa<-PD_fruit_compa

# seed set data for compatibility 
Data_compa_seed = Data_seed[!is.na(Data_seed$Compatibility),] # 19 obs
Data_compa_seed$Compatibility = factor(Data_compa_seed$Compatibility,
                                       levels = c("Self-incompatible","Semi-compatible"))

PD_seed_compa <- ggplot(Data_compa_seed, aes(x = Compatibility, y = lnR))+ 
  scale_y_continuous(limits = c(-2.2,1.3),position = "right")
PD_seed_compa <- PD_seed_compa + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "seed set",y=expression(ln(R[PD])), 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PD_seed_compa<-PD_seed_compa + 
  geom_pointrange(Compa_seed,
                  mapping = aes(x=compatibility, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
seed_compa<-PD_seed_compa

require(gridExtra)
compa_plot<-grid.arrange(fruit_compa,seed_compa,ncol =2)

#### Plot - manipulation level ####
# fruit set data #
mod_mani1 <- emmeans(manipulation_fruit, "Manipulation_level")
multcomp::cld(mod_mani1,Letters = c(letters))
# seed set data #
mod_mani2 <- emmeans(manipulation_seed, "Manipulation_level")
multcomp::cld(mod_mani2,Letters = c(letters))

# create dataset for showing lnR across manipulation level
Manipulation_fruit = data.frame(Manipulation= 
                                  c("Plant","Branch","Inflorescence","Flower"),
                                estimate = c(-0.626,-0.641,-0.661,-0.173),
                                Lower = c(-1.476,-1.054,-0.954,-1.056),
                                Upper = c(0.223,-0.229,-0.367, 0.709))

Manipulation_seed = data.frame(Manipulation = 
                                 c("Branch","Inflorescence","Flower"),
                               estimate = c(-0.210,-0.306,-0.263),
                               Lower = c( -0.632, -0.563,-0.965),
                               Upper = c(0.2118,-0.0485,0.4400))

# Plot #
# fruit set data for manipulation level
Data_mani_fruit = Data_fruit[!is.na(Data_fruit$Manipulation_level),] # 113 obs
Data_mani_fruit$Manipulation_level = factor(Data_mani_fruit$Manipulation_level,
                                            levels = c("Plant","Branch","Inflorescence","Flower"))
PD_fruit_mani <- ggplot(Data_mani_fruit, aes(x = Manipulation_level, y = lnR))+scale_y_continuous(limits = c(-2.2,1.1)) # original data points
PD_fruit_mani <- PD_fruit_mani + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") + # make points spread around
  theme_bw() + theme(panel.grid = element_blank()) + # make backgroud white 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),# locate the position of legend
        legend.background=element_blank(),
        legend.text = element_text(size = 14), # make bacground of legend transparent
        axis.text=element_text(size=14, color = "black")) + # adjust the axis text 
  labs(title = "fruit set",y=expression(ln(R[PD])), x = "",
       axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PD_fruit_mani<-PD_fruit_mani + geom_pointrange(Manipulation_fruit, # add the estimate mean value and CI
                                               mapping = aes(x=Manipulation, y=estimate, ymin=Lower, ymax=Upper), 
                                               size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

fruit_manipulation<-PD_fruit_mani

# seed set data for manipulation level
Data_mani_seed = Data_seed[!is.na(Data_seed$Manipulation_level),] # 42 obs
Data_mani_seed$Manipulation_level = factor(Data_mani_seed$Manipulation_level,
                                           levels = c("Branch","Inflorescence","Flower"))

PD_seed_mani <- ggplot(Data_mani_seed, aes(x = Manipulation_level, y = lnR))+ 
  scale_y_continuous(limits = c(-2.2,1.1),position = "right")
PD_seed_mani <- PD_seed_mani + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "seed set",y=expression(ln(R[PD])), 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PD_seed_mani<-PD_seed_mani + 
  geom_pointrange(Manipulation_seed,
                  mapping = aes(x=Manipulation, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
seed_manipulation<-PD_seed_mani

grid.arrange(fruit_manipulation,seed_manipulation,ncol =2)

#### Plot - pollen type ####
# fruit set data #
mod_pollen1 <- emmeans(mpollen_fruit, "Pollen_type")
multcomp::cld(mod_pollen1,Letters = c(letters))

# seed set data #
mod_pollen2 <- emmeans(mpollen_seed, "Pollen_type")
multcomp::cld(mod_pollen2,Letters = c(letters))

# create dataset for showing lnR for pollen type
Pollen_fruit = data.frame(pollen= c("Mixed","Single"),
                          estimate = c(-0.562,-0.542),
                          Lower = c(-0.912,-0.790),
                          Upper = c(-0.213,-0.294))

Pollen_seed = data.frame(pollen = 
                           c("Mixed","Single"),
                         estimate = c(-0.088,-0.261),
                         Lower = c(-0.376,-0.468),
                         Upper = c(0.199,-0.054))

# Plot #
# fruit set data for pollen type
Data_pollen_fruit = Data_fruit[!is.na(Data_fruit$Pollen_type),] # 113 obs
Data_pollen_fruit$Pollen_type = factor(Data_pollen_fruit$Pollen_type,
                                       levels = c("Mixed","Single"))
PD_fruit_pollen <- ggplot(Data_pollen_fruit, aes(x = Pollen_type, y = lnR))+scale_y_continuous(limits = c(-2.2,1.3)) # original data points
PD_fruit_pollen <- PD_fruit_pollen + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") + # make points spread around
  theme_bw() + theme(panel.grid = element_blank()) + # make backgroud white 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),# locate the position of legend
        legend.background=element_blank(),
        legend.text = element_text(size = 14), # make bacground of legend transparent
        axis.text=element_text(size=14, color = "black")) + # adjust the axis text 
  labs(title = "fruit set",y=expression(ln(R[PD])), x = "",
       axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PD_fruit_pollen<-PD_fruit_pollen + geom_pointrange(Pollen_fruit, # add the estimate mean value and CI
                                                   mapping = aes(x=pollen, y=estimate, ymin=Lower, ymax=Upper), 
                                                   size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

fruit_pollen<-PD_fruit_pollen

# seed set data for pollen type
Data_pollen_seed = Data_seed[!is.na(Data_seed$Pollen_type),] # 47 obs
Data_pollen_seed$Pollen_type = factor(Data_pollen_seed$Pollen_type,
                                      levels = c("Mixed","Single"))

PD_seed_pollen <- ggplot(Data_pollen_seed, aes(x = Pollen_type, y = lnR))+ 
  scale_y_continuous(limits = c(-2.2,1.3),position = "right")
PD_seed_pollen <- PD_seed_pollen + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "seed set",y=expression(ln(R[PD])), 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PD_seed_pollen<-PD_seed_pollen + 
  geom_pointrange(Pollen_seed,
                  mapping = aes(x=pollen, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
seed_pollen<-PD_seed_pollen

grid.arrange(fruit_pollen,seed_pollen,ncol =2)



# to create figure of frequency dixtribution of effext size

# now, grouping four frequency plot into one figure
par(mfrow = c(1,2),mar=c(4.5,4.5,4.5,4.5))
# fruit set data for lnR
hist(Data_fruit$lnR,
     main="Fruit set",
     xlab=expression(ln(R[PD])),
     xlim=c(-2.5,1.5),
     breaks = 20,lwd=1.5,
     cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)
axis(1, at=seq(-2.5,1.5, by=1), labels=seq(-2.5,1.5, by=1),cex.axis = 1.5)
segments(x0 = -2.5,x1= 1.5,y0 = 0,y1 = 0)
text(-2.4,20,"(A)",cex = 1.5)

# seed set data for lnR
hist(Data_seed$lnR,
     main="Seed set",
     xlab=expression(ln(R[PD])),
     xlim=c(-2,0.5),
     breaks = 20,
     cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)
segments(x0 = -2,x1= 0.5,y0 = 0,y1 = 0)
text(-1.95,10,"(B)",cex = 1.5)

# funnel plots of lnRPD for fruit set and seed set -----------------------------------
# fruit set data #
length(unique(Data_fruit$Study)) # 24 

# pick study ID, in which No. of replicates was reported
ID_lnR_fruit<-unique(Data_fruit$Study)
ID_lnR_fruit

# a vector for mean PC of each study
lnR_ave_fruit<-c()

# a vector for study size
SS_lnR_fruit<-c()

# calculate average of PL of each study 
# calculate study size, study size = No. of replicates * No. of treatments

for (i in 2:27) {
  dat_temp_fruit<-Data_fruit[Data_fruit$Study==ID_lnR_fruit[i],]
  lnR_ave_fruit<-c(lnR_ave_fruit,mean(dat_temp_fruit$lnR))
  SS_lnR_fruit<-c(SS_lnR_fruit,sum(dat_temp_fruit$Nopen))
}

# mean PL is -0.53
# check how many PC are above -0.53
length(which(lnR_ave_fruit>-0.53))
# ones below -0.53
length(which(lnR_ave_fruit<=-0.53)) 

# funnel plot with average PL and Study Size of each study in my dataset
par(mfrow = c(1,2),mar=c(4.5,4.5,4.5,4.5))
plot(lnR_ave_fruit,SS_lnR_fruit,xlim=c(-2.1,1),xlab=expression(ln(R[PD])),
     ylim = c(0,150),ylab="Study size",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,tck=-0.02)

# add lines to the funnel plot
# add a vertical line representing the mean of PL 
lines(x=rep(-0.53,5),y=seq(0,150,length=5))
# draw a symmetric triangle base on mean PL and minimum PL
# add a line to the left of the vertical line according to the left end of lnR and the maximum study size
d<--0.53-(-2.1)
# set the min and max of X for the line
X<-c(-0.53-d,-0.53)
# set min and max of Y for the line
Y<-c(0,150)
# calculate the slope for the line
lm(Y~X)
# draw the line
X<-seq(-0.53-d,-0.53,length=5)
Y<-200.64  +  95.54  *X
lines(x=X,y=Y)

# add a line to the right of the vertical line according to the right end of lnR and the maximum study size
X<-c(-0.53+d,-0.53)
Y<-c(0,150)
lm(Y~X)
X<-seq(-0.53,-0.53+d,length=5)
Y<-99.36 -95.54 *X
lines(x=X,y=Y)
# add a horizontal line to close the triangle
lines(x=seq(-0.53-d,-0.53+d,length=5),y=rep(0,5))

text(-2.0,150,"(A)",cex=1.5,font=1.5)
#------------------------------------------------------------------------------------------#

# seed set #
length(unique(Data_seed$Study)) # 12 

# pick study ID, in which No. of replicates was reported
ID_lnR_seed<-unique(Data_seed$Study)
ID_lnR_seed

# a vector for mean PC of each study
lnR_ave_seed<-c()

# a vector for study size
SS_lnR_seed<-c()

# calculate average of lnR of each study 
# calculate study size, study size = No. of replicates * No. of treatments

for (i in 1:24) {
  dat_temp<-Data_seed[Data_seed$Study==ID_lnR_seed[i],]
  lnR_ave_seed<-c(lnR_ave_seed,mean(dat_temp$lnR))
  SS_lnR_seed<-c(SS_lnR_seed,sum(dat_temp$Nopen))
}

# mean lnR is -0.228
# check how many lnR are above -0.228
length(which(lnR_ave_seed>-0.228)) # 7
# ones below -0.42
length(which(lnR_ave_seed<=-0.228)) # 5

# funnel plot with average lnR and Study Size of each study in my dataset
plot(lnR_ave_seed,SS_lnR_seed,xlim=c(-1,0.5),xlab=expression(ln(R[PD])),
     ylim = c(0,150),ylab="Study size",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,tck=-0.02)

# add lines to the funnel plot
# add a vertical line representing the mean of PL 
lines(x=rep(-0.228,5),y=seq(0,150,length=5))
# draw a symmetric triangle base on mean PL and minimum PL
# add a line to the left of the vertical line according to the left end of lnR and the maximum study size
d<--0.228-(-1)
# set the min and max of X for the line
X<-c(-0.228-d,-0.228)
# set min and max of Y for the line
Y<-c(0,150)
# calculate the slope for the line
lm(Y~X)
# draw the line
X<-seq(-0.228-d,-0.228,length=5)
Y<- 194.3    + 194.3  *X
lines(x=X,y=Y)

# add a line to the right of the vertical line according to the right end of lnR and the maximum study size
X<-c(-0.228+d,-0.228)
Y<-c(0,150)
lm(Y~X)
X<-seq(-0.228,-0.228+d,length=5)
Y<-105.7 -194.3  *X
lines(x=X,y=Y)
# add a horizontal line to close the triangle
lines(x=seq(-0.228-d,-0.228+d,length=5),y=rep(0,5))

text(-0.95,150,"(B)",cex=1.5,font=1.5)

# ---------------------- end ---------------------------#