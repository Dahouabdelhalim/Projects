
# # # # # # # # # # # # # # # 
# pollination contribution #
#      --- Aruhan, Wopke #
#    - 2021.02.02 -    #

# clean the environment #
rm(list=ls())

# set work space #
setwd("F:/R/RStudio/R/r_work/Meta-analysis_Aruhan")

# install packages #
# install.packages("nlme")
# install.packages("car")
# install.packages("emmeans")
# install.packages("ggplot2")
# install.packages("gridExtra")

# load libraries #
library(nlme) # for fit linear mixed models
library(car) # # enables the qqPlot function to check residuals
library(emmeans) # to do multiple comparison and extract lsmeans
library(ggplot2) # to plot forest plot 

# load dataset---- 
PC_ALL<- read.csv("Pollinators contribution_DATA_Olhnuud et al.csv",header = T)
names(PC_ALL) # check the column name
dim(PC_ALL) # 164 observations with 19 variables. 

PC_data = PC_ALL[,c(2,3,8,9,10,11,12,13,14,15,16,17,18,19)]
names(PC_data) 

# calulate effect size (Xopen-Xexclusion)/Xopen for whole data set
PC_data$PC = (PC_data$Xopen-PC_data$Xexclusion)/PC_data$Xopen

# check how many PC equals to 1 
length(which(PC_data$PC==1)) # 25 records

# assign variables----
PC_data$Study <- factor(PC_data$Study)
levels(PC_data$Study)
length(levels(PC_data$Study)) # 29 studies

PC_data$Experiment<-factor(PC_data$Experiment)
levels(PC_data$Experiment)
length(levels(PC_data$Experiment)) # 54 experiments

PC_data$Continent<-factor(PC_data$Continent)
levels(PC_data$Continent)  #  Asia, Europe, North America, South America, Oceania 5 levels

PC_data$Country<-factor(PC_data$Country)
levels(PC_data$Country) # 12 countries

PC_data$Climate<-factor(PC_data$Climate)  
levels(PC_data$Climate) # Arid, Continental, Temperate, Tropical 4 levels

PC_data$Cultivar<-factor(PC_data$Cultivar)
length(levels(PC_data$Cultivar)) # 29 cultivars
table(PC_data$Cultivar)
sum(table(PC_data$Cultivar)) # 158 records, cultivars are unknown in 6 records

PC_data$Compatibility<-factor(PC_data$Compatibility)
levels(PC_data$Compatibility) # "Self-incompatible" "Semi-compatible" "self-compatible" 
table(PC_data$Compatibility)
sum(table(PC_data$Compatibility)) # 143 records, compatibility in 27 records were not defined

PC_data$Manipulation_level<-factor(PC_data$Manipulation_level)
levels(PC_data$Manipulation_level)#  Branch, Inflorescence, Plant 3 levels

# Create a subgroup dataset for fruit-set data.
Data_fruit = PC_data[!(PC_data$Measurement == "seed set"),] # 100 records
Data_seed = PC_data[!(PC_data$Measurement == "fruit set"),] # 64 records

# It is necessary to use droplevels function
# if we do not use this, this will bother the analysis,
# there is still 54 levels (equal to experiment levels from whole dataset)
# for experiments after separating data 
# by measurement into fruit set data and seed set data
Data_seed<-droplevels(Data_seed)
Data_fruit<-droplevels(Data_fruit)

# Data exploration-------------------------------------------------------------------------------------
# fruit set data
# histogram with added parameters
hist(Data_fruit$PC,
     main="Fruit set",
     xlab="PC",
     xlim=c(0,1),
     breaks = 20,cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)

hist(Data_seed$PC,
     main="Seed set",
     xlab="PC",
     xlim=c(0,1),
     breaks = 20,cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)

# to explore the correlation between study or correlation between experiment
plot(as.numeric(Data_fruit$Study), Data_fruit$PC)
plot(as.numeric(Data_fruit$Experiment), Data_fruit$PC)
# there is clear correlation between study, Study may could be considered as random effec component

# explore patterns to know about data #
datx = Data_fruit[-c(4,6,15),]
dim(datx)
colours = c("red","orange","green","blue")
op=par()
par(mar=c(5.1,9,4.1,2.1))
with(datx, plot(PC,as.numeric(Cultivar),pch=16,col=colours[as.numeric(datx$Continent)],
                yaxt = "n", ylab = ""))
cultivar.names = unique(datx$Cultivar)
cultivar.nrs = as.numeric(datx$Cultivar)
cultivar.names[cultivar.nrs]
axis(side=2,at=1:23,labels=sort(cultivar.names),hadj=1,las=1)
for (i in c(1:23)){
  lines(c(0,1),c(i,i),lty=3,col="grey",lwd=0.5)
}
legend(0.35,23,c("Asia","Europe","N America","S America"),pch=16,col=colours,bty="y")
# There are two varieties (9 and 10: Gala and Golden Delicious) that are tested 
# in more than one continent
datx[as.numeric(datx$Cultivar)==9,]$Cultivar
datx[as.numeric(datx$Cultivar)==10,]$Continent

# Explore whether interaction/confounding effect on cultivar and continent-----------------------------------------#
# A figure with for each variety the experiments coded
# by colours and the continents by symbols
plot.symbols = 15:18
with(datx, plot(PC,as.numeric(Cultivar), col=datx$Experiment, pch=plot.symbols[as.numeric(datx$Continent)],
                yaxt = "n", ylab = "",cex=1.2))
axis(side=2,at=1:23,labels=sort(cultivar.names),hadj=1,las=1)
for (i in c(1:23)){
  lines(c(0,1),c(i,i),lty=3,col="grey",lwd=0.5)
}
legend(0.35,23,c("Asia","Europe","N America","S America"),cex=1.2,pch=15:18,col="black",bty="y")
# This figure shows that all several independent experiments in Asia show low PC in Red Gold. 
# It also shows that independent experiments in Asia and Europe 
# for Golden Delicious find higher pollination in Europe than in Asia. 
# But for Gala there is a high value in Asia,similar to European data. 
# There is no support for a model with continent and cultivar
# as additive effects in the paper, 
# because there is not enough data to really understand the continent effect. 


# Explore whether there is interaction/confounding effect between cultivar and manipulation level
# draw a figure with for each variety the experiments coded by colours and
# manipulation level by symbols
colours = c("red","green","blue")
with(datx, plot(PC,as.numeric(Cultivar), pch=16,col=colours[as.numeric(datx$Manipulation_level)],
                yaxt = "n", ylab = "",cex=1.2))
axis(side=2,at=1:23,labels=sort(cultivar.names),hadj=1,las=1)
for (i in c(1:23)){
  lines(c(0,1),c(i,i),lty=3,col="grey",lwd=0.5)
}
legend(0.35,23,c("Branch","Inflorescence","Plant"),cex=1.2,colours,bty="y")
# This figure shows pollination contribution for each cultivar, 
# manipulation level coded by colors. 
# In fruit set data, 
# 71 out of 100 points are from experiment done on manipulation level for branch, 
# 10 points are on inflorescence level, and 19 points are on plant level. 
# The pollination treatment is done mostly on branch level, 
# in same experiment, pollination contribution varied between different cultivars. 
# Three manipulation levels appeared only in experiment of Gala and Fuji. 
# To conclude, I don¡¯t think effect of manipulation level is strong. 
# The main effect would be only from cultivars. 

# Explore effect of compatibility of cultivars here------------------------------------------------------------------#
# draw a plot for each variety in continent coded by colours and 
# compatibility by symbols
plot.symbols = 15:19
with(datx, plot(PC,as.numeric(Cultivar), col=datx$Continent, pch=plot.symbols[as.numeric(datx$Compatibility)],
                yaxt = "n", ylab = "",cex=1.2))
axis(side=2,at=1:23,labels=sort(cultivar.names),hadj=1,las=1)
for (i in c(1:23)){
  lines(c(0,1),c(i,i),lty=3,col="grey",lwd=0.5)
}
legend(0.35,23,c("self-compatible","self-incompatible","semi-compatible"),cex=1.2,pch=15:18,col="black",bty="y")
# self-compatible cultivar Red Gold indeed showing lower PC than other cultivars, 
# and semi-compatible cultivar Cox showing about 60% PC, 
# and self-compatible cultivars showing higher value, 
# beside few points in Honeycrisp, Gala, Fuji and Elstar. 
# In addition, Borkhausen, Aroma, and Amorosa could not be classified since no evidence. 
# I am not sure there is clear pattern on effect of compatibility on PC,
# and there is not enough sample for self-compatible and semi-compatible group,
# so, the effect of compatibility is uncertain.

# Explore complex model by practical way-------------------------------------------------------------------------------------------#
# more complex model to test: the interaction for continent and cultivar effect
B1<-lm(PC ~ Cultivar*Continent,data = Data_fruit,na.action = na.exclude)
summary(B1)
# there is only interaction for Europe and Gala, and so many empty(NA) in interaction
# the interaction is neither available nor suitable

B2<-lm(PC ~ Cultivar + Continent,data = Data_fruit,na.action = na.exclude)
summary(B2)
# there is also no strong support for confounding effect of cultivar and continent,
# as not all cultivar appeared all continent, only very few.

B3<-lme(PC ~ Cultivar + Continent,random = ~1|Study,data = Data_fruit,na.action = na.exclude)
summary(B3)
# there is also no strong support for confounding effect of cultivar and continent,
# as not all cultivar appeared all continent, only very few.

# The Random Intercept model or The Random Intercept and Slope Model?
B4<-lme(PC ~ Cultivar + Continent,random = ~1 + Continent|Experiment,method = "REML",data = Data_fruit,na.action = na.exclude)
summary(B4)
# there is also no strong support for confounding effect of cultivar and continent,
# as not all cultivar appeared all continent, only very few.

B5<-lm(PC ~ Cultivar*Country,data = Data_fruit,na.action = na.exclude)
summary(B5)
# so many empty in interaction, also no suppot for complex model

#### State for model selection ####
# To choose best random effect component to explain between-study variability 
# compare AIC (or BIC) of models fitted by "ML", and refit final model by REML
# start with estimation of means effect size,  
# then try to explain more variability of data according to add meaningful moderators

# 1.overall effect of pollinator contribution----
M0<- gls(PC~ 1, data = Data_fruit,method = "ML")
M1 <- lme(PC~ 1, random = ~1|Study, data=Data_fruit,method = "ML")
M2 <- lme(PC~ 1, random = ~1|Experiment, data=Data_fruit,method = "ML")
M3 <- lme(PC~ 1, random = ~1|Study/Experiment, data=Data_fruit,method = "ML")
anova(M0,M1,M2,M3) 
summary(M1)
#  M1 is better than M2, M1 is similar with M0 and M3
#  supposed to choose simpler model, fixed-effects model
#  according to the possibility of correlation between study,
#  however, study effect could be estimated though between-study variablility is small.
#  sensitivity analysis would be appllied in this case, 
#  compare fixed-effects model and mixed-effect model.

# Refit models by REML to compare results (mean effect size PC)-----
# fixed-effects model 
M0_REML <- gls(PC~ 1,data=Data_fruit)
summary(M0_REML) 
# Approximate 95% confidence intervals
# 
# Coefficients:
#   lower      est.     upper
# (Intercept) 0.632551 0.6920163 0.7514817
# attr(,"label")
# [1] "Coefficients:"
# 
# Residual standard error:
#   lower      est.     upper 
# 0.2631314 0.2996918 0.3481445

intervals(M0_REML)

# mixed-effects model 
M1_REML <- lme(PC~ 1,random = ~1|Study,data=Data_fruit,method = "REML")
summary(M1_REML) 
# Linear mixed-effects model fit by REML
# Data: Data_fruit 
# AIC      BIC   logLik
# 49.94081 57.72617 -21.9704
# 
# Random effects:
#   Formula: ~1 | Study
#         (Intercept)  Residual
# StdDev:   0.1071035 0.2815502
# 
# Fixed effects: PC ~ 1 
#                Value  Std.Error  DF t-value p-value
# (Intercept) 0.7117252 0.03928537 77 18.1168       0
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -2.2005520 -0.5005003  0.2322757  0.6032405  1.3762939 
# 
# Number of Observations: 100
# Number of Groups: 23 

intervals(M1_REML)

# model validation for fixed-effects model #
par(mfrow = c(2,2),mar=c(6,6,4,4))
hist(resid(M0_REML),freq = F,
     main = "Histogram of resid of fixed-effects model",
     cex.main=0.8)
qqPlot(resid(M0_REML))

# model validation for random effect model #
hist(resid(M1_REML),freq = F,
     main = "Histogram of resid of mixed-effects model",
     cex.main=0.8)
qqPlot(resid(M1_REML))
# compare two models, we prefer mixed-effects model, so choose mixed effects model to present.
# since between study variability was estimated, and model validation has better performance. 

# Test effect of measurement for fruit set (initial fruit set or final fruit)----
M1_REML_measurement <- lme(PC~ Fruit_set_measurement,random = ~1|Study,data=Data_fruit,method = "REML")
summary(M1_REML_measurement)
anova(M1_REML_measurement) # p = 0.746 > 0.05, 
# indicating fruit_set_measurement is no effect, PC is not different between initial fruit and final fruit. 

# 2. Does the magnitude of pollinator contribution differ between cultivar? ---- 
Mfruit_cultivar0 <- gls(PC~ Cultivar,data=Data_fruit,na.action = na.exclude,method = "ML")
Mfruit_cultivar1 <- lme(PC~ Cultivar, random = ~1|Study, data=Data_fruit,na.action = na.exclude,method = "ML")
Mfruit_cultivar2 <- lme(PC~ Cultivar, random = ~1|Experiment, data=Data_fruit,na.action = na.exclude,method = "ML")
Mfruit_cultivar3 <- lme(PC~ Cultivar, random = ~1|Study/Experiment, data=Data_fruit,na.action = na.exclude,method = "ML")
anova(Mfruit_cultivar0,Mfruit_cultivar1,Mfruit_cultivar2,Mfruit_cultivar3)
# Mfruit_cultivar3 is the best, model with study/Experiment as random effect 
# however, just this one case (there is no big difference between AICs of Mfruit_cultivar1 and Mfruit_cultivar2)
# so choose model with study as random effect 

# Refit models by REML to compare results-----
# fixed-effects model 
Mfruit_cultivar0.1 <- gls(PC~ Cultivar,data=Data_fruit,na.action = na.exclude)
summary(Mfruit_cultivar0.1)
anova(Mfruit_cultivar0.1) # p <.0001, PC is varies between cultivars 

mod_cul = lsmeans(Mfruit_cultivar0.1,~Cultivar)
multcomp::cld(mod_cul,Letters = c(letters))

# mixed-effects model by REML
Mfruit_cultivar1.1 <- lme(PC~ Cultivar,random = ~1|Study,data=Data_fruit,na.action = na.exclude,method = "REML")
summary(Mfruit_cultivar1.1)
# Random effects:
#   Formula: ~1 | Study
#         (Intercept)  Residual
# StdDev:   0.1524527 0.1471874
# between-study variability is estimated, 
# between-study SD is equal to within-study SD in this case

anova(Mfruit_cultivar1.1) # p <.0001, PC is different between cultivars.

mod_cul_random = lsmeans(Mfruit_cultivar1.1, ~ Cultivar)
multcomp::cld(mod_cul_random,Letters = c(letters))
# confidence intervals became wider when using mixed-effects model

# plot #
Mfruit_cultivar_plot<- lme(PC~ Cultivar-1,random = ~1|Study,data=Data_fruit,na.action = na.exclude,method = "REML")
summary(Mfruit_cultivar_plot)

# extract mean,confidence intervals for each cultivar
Mean_cultivar = intervals(Mfruit_cultivar_plot)$fixed[,2]
Upper_cultivar = intervals(Mfruit_cultivar_plot)$fixed[,3]
Lower_cultivar = intervals(Mfruit_cultivar_plot)$fixed[,1]

# reshape data set
cultivar_names = unique(Data_fruit$Cultivar)
Data_cul = data.frame(cultr = sort(cultivar_names),mean_cultr = Mean_cultivar,
                      lower_CI = Lower_cultivar, upper_CI = Upper_cultivar)
names(Data_cul)

#  plot #
par(mfrow = c(1,1),mar=c(4.7,4,4,2.5))
Data_Cultivar_order = Data_cul[order(Data_cul$mean_cultr),]

dotchart(Data_Cultivar_order$mean_cultr,labels = Data_Cultivar_order$cultr,
         xlim = c(0,1.5),cex = 1.19,xlab = "PC",pch = 16,main = "Fruit set")

for (i in 1:nrow(Data_Cultivar_order)) {
  lines(x=c(Data_Cultivar_order[i,3],Data_Cultivar_order[i,4]), y=c(i,i))
}
abline(v=0,lty = "dashed",col = "grey19")

#### model validation for fixed-effects model ###
par(mfrow = c(2,2),mar=c(6,6,4,4))
hist(resid(Mfruit_cultivar0.1),freq = F)
qqPlot(resid(Mfruit_cultivar0.1))
plot(fitted(Mfruit_cultivar0.1),resid(Mfruit_cultivar0.1))

#### model validation for mixed-effects model ###
par(mfrow = c(2,2),mar=c(6,6,4,4))
hist(resid(Mfruit_cultivar1.1),freq = F)
qqPlot(resid(Mfruit_cultivar1.1))
plot(fitted(Mfruit_cultivar1.1),resid(Mfruit_cultivar1.1))

# comparing qqplot of residuals from above fitted models, mixed-effects model performed better
# as most points fall approximately along this reference line, we can assume normality.

# 3. Does the magnitude of pollinator contribution differ between continents?---- 
Mfruit_continent0 <- gls(PC~ Continent,data=Data_fruit,method = "ML")
Mfruit_continent1 <- lme(PC~ Continent, random = ~1|Study, data=Data_fruit,method = "ML")
Mfruit_continent2 <- lme(PC~ Continent, random = ~1|Experiment, data=Data_fruit,method = "ML")
Mfruit_continent3 <- lme(PC~ Continent, random = ~1|Study/Experiment, data=Data_fruit,method = "ML")
anova(Mfruit_continent0,Mfruit_continent1,Mfruit_continent2,Mfruit_continent3)
# here, fixed-effects model has better performance, but still have to compare fixed- and mixed- effects model

# Refit three models to compare the results---- 
Mfruit_continent0.1 <- gls(PC~ Continent,data=Data_fruit)
Mfruit_continent1.1 <- lme(PC~ Continent, random = ~1|Study, data=Data_fruit)
Mfruit_continent2.1 <- lme(PC~ Continent, random = ~1|Experiment, data=Data_fruit)
summary(Mfruit_continent0.1)
summary(Mfruit_continent1.1) 
# there is contribution of between study variability in model with study effect as random effect
summary(Mfruit_continent2.1) 
# no between-experiment variablity was estimated (almost zero), so drop out. 

anova(Mfruit_continent0.1) # p = 0.2426 > 0.05
anova(Mfruit_continent1.1) # p = 0.5584 > 0.05
# to conclude: no matter what model using, PC is not different between continents. 

# to extract estimate from fixed-effects model
mod_fruit_conti = lsmeans(Mfruit_continent0.1,~Continent)
# to extract estimate from mixed-effects model
mod_fruit_conti_random = lsmeans(Mfruit_continent1.1,~ Continent)
# Results from two model were very similar.

#### model validation for fixed-effects model ###
par(mfrow = c(2,3),mar=c(6,6,4,4))
hist(resid(Mfruit_continent0.1),freq = F)
qqPlot(resid(Mfruit_continent0.1)) # residuals 
plot(fitted(Mfruit_continent0.1),resid(Mfruit_continent0.1))

#### model validation for mixed-effects model ###
hist(resid(Mfruit_continent1.1),freq = F)
qqPlot(resid(Mfruit_continent1.1))
plot(fitted(Mfruit_continent1.1),resid(Mfruit_continent1.1))
# prefer mixed effect model

# 4. Does the magnitude of pollinator contribution is differ between countries?----
Mcountry0 <- gls(PC~ Country,data = Data_fruit,method = "ML")
Mcountry1 <- lme(PC~ Country, random = ~1|Study,data = Data_fruit,method = "ML")
Mcountry2 <- lme(PC~ Country, random = ~1|Experiment,data = Data_fruit,method = "ML")
Mcountry3 <- lme(PC~ Country, random = ~1|Study/Experiment,data = Data_fruit,method = "ML")
anova(Mcountry0,Mcountry1,Mcountry2,Mcountry3)
# no difference between first three models, do not consider final complex one

# Refit the first three models by REML to compare results-------
Mcountry0.1 <- gls(PC~ Country,data = Data_fruit)
Mcountry1.1 <- lme(PC~ Country, random = ~1|Study,data = Data_fruit)
Mcountry2.1 <- lme(PC~ Country, random = ~1|Experiment,data = Data_fruit)

summary(Mcountry0.1)
anova(Mcountry0.1) # p = 0.0319*, significant difference between countries

summary(Mcountry1.1)
anova(Mcountry1.1) # p = 0.0913, no difference between countries, differ from fixed-effects model
# Random effects:
#   Formula: ~1 | Study
# (Intercept)  Residual
# StdDev:  0.00692139 0.2857536

# estimates of fixed-effects model and mixed-effects (study) model are very similar,
# the between study variability is very small

summary(Mcountry2.1) # there is no between experiment variability estimated, so drop this model out

# to see estimate from fixed-effects model
mod_country<-lsmeans(Mcountry0.1, ~ Country)
multcomp::cld(mod_country,Letters = c(letters))

# to see estimate from mixed-effects model 
mod_country_random<-lsmeans(Mcountry1.1, ~ Country)
# comparision indicatin no markable difference between results of fixed-effects model and mixed-effects model

# 5. does the pollinator dependence of cultivar is different between climate zones?----
Mclimate0 <- gls(PC~ Climate,data = Data_fruit,method = "ML")
Mclimate1 <- lme(PC~ Climate, random = ~1|Study,data = Data_fruit,method = "ML")
Mclimate2 <- lme(PC~ Climate, random = ~1|Experiment,data = Data_fruit,method = "ML")
Mclimate3 <- lme(PC~ Climate, random = ~1|Study/Experiment,data = Data_fruit,method = "ML")
anova(Mclimate0,Mclimate1,Mclimate2,Mclimate3)
# no difference fixed effect model and mixed-effects model

# Refit two models by REML to compare results-----
Mclimate0.1 <- gls(PC~ Climate,data = Data_fruit)
Mclimate1.1 <- lme(PC~ Climate, random = ~1|Study,data = Data_fruit)
Mclimate2.1 <- lme(PC~ Climate, random = ~1|Experiment,data = Data_fruit)
anova(Mclimate0.1,Mclimate1.1,Mclimate2.1) # no difference between three models

summary(Mclimate0.1)
anova(Mclimate0.1) # p = 0.318, no
mod_fruit_climate0.1 = lsmeans(Mclimate0.1, ~ Climate)

summary(Mclimate1.1)
# Random effects:
#   Formula: ~1 | Study
#         (Intercept)  Residual
# StdDev:   0.1062987 0.2835932

anova(Mclimate1.1) # p = 0.6208, no 
mod_fruit_climate_random = lsmeans(Mclimate1.1, ~ Climate)

summary(Mclimate2.1)
# second model is better, since explain variablity more than third model with experiment as random effect
# betweeen study variability was estimated (not equal to zero)

# conclusion: no matter what model using, PC is not differ across climate zones

#### model validation for fixed-effects model ###
par(mfrow = c(2,3),mar=c(6,6,4,4))
hist(resid(Mclimate0.1),freq = F)
qqPlot(resid(Mclimate0.1))
plot(fitted(Mclimate0.1),resid(Mclimate0.1))

#### model validation for mixed-effects model ###
hist(resid(Mclimate1.1),freq = F)
qqPlot(resid(Mclimate1.1))
plot(fitted(Mclimate1.1),resid(Mclimate1.1))
# residuals looks good

# 6. Does the magnitude of pollinator contribution is differ between compatibility?----
Mcompa0 <- gls(PC~ Compatibility,data = Data_fruit,na.action = na.exclude,method = "ML")
Mcompa1 <- lme(PC~ Compatibility, random = ~1|Study,data = Data_fruit,na.action = na.exclude,method = "ML")
Mcompa2 <- lme(PC~ Compatibility, random = ~1|Experiment,data = Data_fruit,na.action = na.exclude,method = "ML")
Mcompa3 <- lme(PC~ Compatibility, random = ~1|Study/Experiment,data = Data_fruit,na.action = na.exclude,method = "ML")
anova(Mcompa0,Mcompa1,Mcompa2,Mcompa3)
# no difference between the first three models, omit final complex one 

# Refit the first two models by REML to compare results----
Mcompa0.1 <- gls(PC~ Compatibility,data = Data_fruit,na.action = na.exclude)
Mcompa1.1 <- lme(PC~ Compatibility, random = ~1|Study,data = Data_fruit,na.action = na.exclude)
Mcompa2.1 <- lme(PC~ Compatibility, random = ~1|Experiment,data = Data_fruit,na.action = na.exclude)

summary(Mcompa0.1)
anova(Mcompa0.1) # p <.0001, yes
mod_fruit_compa = lsmeans(Mcompa0.1, ~ Compatibility)
multcomp::cld(mod_fruit_compa,Letters = c(letters))

summary(Mcompa1.1)
anova(Mcompa1.1) # p <.0001, yes
mod_fruit_compa_random = lsmeans(Mcompa1.1, ~ Compatibility)
multcomp::cld(mod_fruit_compa_random,Letters = c(letters))

summary(Mcompa2.1) # between-expeirment variablity is almost zero, so drop out

# PC is differ between different compatibility type
# However, not all cultivar specified with compatibility 
# (11 missing data records in this variable, 89 observation records out of 100 have value of compatibility)
# results are unstable, hence, we do not stress it in discussion

#### model validation for fixed-effects model ###
par(mfrow = c(2,3),mar=c(6,6,4,4))
hist(resid(Mcompa0.1),freq = F)
qqPlot(resid(Mcompa0.1))
plot(fitted(Mcompa0.1),resid(Mcompa0.1))

#### model validation for mixed-effects model ###
hist(resid(Mcompa1.1),freq = F)
qqPlot(resid(Mcompa1.1))
plot(fitted(Mcompa1.1),resid(Mcompa1.1))

# 7. Does the magnitude of pollinator contribution is associated with manipulation level?----
Mfruit_mani0 <- gls(PC~ Manipulation_level,data=Data_fruit,method = "ML")
Mfruit_mani1 <- lme(PC~ Manipulation_level, random = ~1|Study, data=Data_fruit,method = "ML")
Mfruit_mani2 <- lme(PC~ Manipulation_level, random = ~1|Experiment, data=Data_fruit,method = "ML")
Mfruit_mani3 <- lme(PC~ Manipulation_level, random = ~1|Study/Experiment, data=Data_fruit,method = "ML")
anova(Mfruit_mani0,Mfruit_mani1,Mfruit_mani2,Mfruit_mani3) 
# no difference between model1 and model2, 
# no difference between model3and model4,
# Mfruit_mani1 (study as random) is better than Mfruit_mani2 (experiment as random)
# in this case, choose simpler model in principal

# Refit the model by REML to compare results-----
Mfruit_mani0.1 <- gls(PC~ Manipulation_level,data=Data_fruit)
summary(Mfruit_mani0.1)
anova(Mfruit_mani0.1) # p = 0.3887, no

Mfruit_mani1.1 <- lme(PC~ Manipulation_level, random = ~1|Study, data=Data_fruit)
summary(Mfruit_mani1.1)
# Linear mixed-effects model fit by REML
# Data: Data_fruit 
#    AIC      BIC    logLik
# 56.61813 69.49169 -23.30907
# 
# Random effects:
#   Formula: ~1 | Study
#           (Intercept)  Residual
# StdDev:   0.1148331 0.2792932
# between-study variability is estimated in this model
anova(Mfruit_mani1.1) # p = 0.3092, no

# to see estimate from different models
mod_fruit_mani<-lsmeans(Mfruit_mani0.1, ~ Manipulation_level)
mod_fruit_mani_random<-lsmeans(Mfruit_mani1.1, ~ Manipulation_level)

#### model validation for fixed-effects model ###
par(mfrow = c(2,3),mar=c(6,6,4,4))
hist(resid(Mfruit_mani0.1),freq = F)
qqPlot(resid(Mfruit_mani0.1))
plot(fitted(Mfruit_mani0.1),resid(Mfruit_mani0.1))

#### model validation for mixed-effects model ###
hist(resid(Mfruit_mani1.1),freq = F)
qqPlot(resid(Mfruit_mani1.1))
plot(fitted(Mfruit_mani1.1),resid(Mfruit_mani1.1))
# I prefer mixed-effects model because of performace of histagram

# Analysis process of pollination contribution for seed set -----------------------------------------------------------------------------
# 1.overall effect of pollinator contribution for seed set----
M0_seed <- gls(PC~ 1,data=Data_seed,method = "ML")
M1_seed <- lme(PC~ 1,random = ~1|Study,data=Data_seed,method = "ML")
M2_seed <- lme(PC~ 1,random = ~1|Experiment,data=Data_seed,method = "ML")
M3_seed <- lme(PC~ 1,random = ~1|Study/Experiment,data=Data_seed,method = "ML")
anova(M0_seed,M1_seed,M2_seed,M3_seed) 
# no difference between fixed- and mixed- effects (study or experiment) models

# Refit model by REML to compare results of models------
M0_seed_REML <- gls(PC~ 1,data=Data_seed)
M1_seed_REML <- lme(PC~ 1,random = ~1|Study,data=Data_seed)
summary(M0_seed_REML)
summary(M1_seed_REML) # between study variability is very close to zero, 
# so choose fixed-effects model

# print confidence intervals
intervals(M0_seed_REML) # mean = 0.62, CI(0.55, 0.69) 


# 2. Does the magnitude of pollinator contribution is differ across cultivars?----
Mseed_cultivar0<-gls(PC~ Cultivar,data=Data_seed,na.action = na.exclude,method = "ML")
Mseed_cultivar1<-lme(PC~ Cultivar,random = ~1|Study,data=Data_seed,na.action = na.exclude,method = "ML")
Mseed_cultivar2<-lme(PC~ Cultivar,random = ~1|Experiment,data=Data_seed,na.action = na.exclude,method = "ML")
Mseed_cultivar3<-lme(PC~ Cultivar,random = ~1|Study/Experiment,data=Data_seed,na.action = na.exclude,method = "ML")
anova(Mseed_cultivar0,Mseed_cultivar1,Mseed_cultivar2,Mseed_cultivar3)

# Refit model by REML to compare results of models-----
Mseed_cultivar0.1<-gls(PC~ Cultivar,data=Data_seed,na.action = na.exclude)
Mseed_cultivar1.1<-lme(PC~ Cultivar,random = ~1|Study,data=Data_seed,na.action = na.exclude)
Mseed_cultivar2.1<-lme(PC~ Cultivar,random = ~1|Experiment,data=Data_seed,na.action = na.exclude)
summary(Mseed_cultivar0.1)
anova(Mseed_cultivar0.1) # p <.0001
# yes, PC for seed set is different between cultivars

summary(Mseed_cultivar1.1)
# Linear mixed-effects model fit by REML
# Data: Data_seed 
#    AIC      BIC   logLik
# 10.87243 49.71665 17.56379
# 
# Random effects:
#   Formula: ~1 | Study
#         (Intercept)  Residual
# StdDev:   0.1381623 0.1111431
anova(Mseed_cultivar1.1) # p <.0001
# yes, PC for seed set is different between cultivars

summary(Mseed_cultivar2.1)
# the between-study variability (experiment) is extimated sd also zero,
# so drop out

# multiple comparison using fixed-effects model
mod_seed<-lsmeans(Mseed_cultivar0.1, ~ Cultivar)
multcomp::cld(mod_seed,Letters = c(letters))

# to see estimate for each cultivars 
Mseed_cultivar_plot<-gls(PC~ Cultivar-1,data=Data_seed,na.action = na.exclude)
summary(Mseed_cultivar_plot)

# extract mean,confidence intervals for each cultivar 
# (results are same when using function lsmean)
Mean_cultivar_seed = intervals(Mseed_cultivar_plot)$coef[,2]
Upper_cultivar_seed = intervals(Mseed_cultivar_plot)$coef[,3]
Lower_cultivar_seed = intervals(Mseed_cultivar_plot)$coef[,1]

# plot #
# reshape data set
cultivar_names_seed = unique(Data_seed$Cultivar)
Data_cul_seed = data.frame(cultr_seed = sort(cultivar_names_seed),mean_cultr_seed = Mean_cultivar_seed,
                           lower_CI_seed = Lower_cultivar_seed, upper_CI_seed = Upper_cultivar_seed)
names(Data_cul_seed)

#  plot #
Cultivar_seed_order = Data_cul_seed[order(Data_cul_seed$mean_cultr_seed),]
par(mfrow = c(1,1),mar=c(4.5,4,4,2.5))
par(bg = "white", fg = "black",col.axis = "black",col.lab = "black")
dotchart(Cultivar_seed_order$mean_cultr_seed,labels = Cultivar_seed_order$cultr_seed,
         xlim = c(0,1.5),pch = 16,cex = 1.19,xlab = "PC",main = "Seed set")

for (i in 1:nrow(Cultivar_seed_order)) {
  lines(x=c(Cultivar_seed_order[i,3],Cultivar_seed_order[i,4]), y=c(i,i))
}
abline(v=0,lty = "dashed",col = "grey19")

# estimates from random effect model
mod_seed_random<-lsmeans(Mseed_cultivar1.1, ~ Cultivar)
multcomp::cld(mod_seed_random,Letters = c(letters))

# model validation for fixed-effects model #
par(mfrow = c(1,1),mar=c(6,6,4,4))
hist(resid(Mseed_cultivar0.1),freq = F) # looks ok
qqPlot(resid(Mseed_cultivar0.1)) # looks terrible
plot(fitted(Mseed_cultivar0.1),resid(Mseed_cultivar0.1))

# model validation for mixed-effects model #
hist(resid(Mseed_cultivar1.1),freq = F) # look ok
qqPlot(resid(Mseed_cultivar1.1)) # looks terrible
plot(fitted(Mseed_cultivar1.1),resid(Mseed_cultivar1.1)) # this looks ok

# 3. Does pollinator contribution differ between continents?---- 
Mseed_continent0 <- gls(PC~ Continent,data=Data_seed,method = "ML")
Mseed_continent1 <- lme(PC~ Continent, random = ~1|Study, data=Data_seed,method = "ML")
Mseed_continent2 <- lme(PC~ Continent, random = ~1|Experiment, data=Data_seed,method = "ML")
Mseed_continent3 <- lme(PC~ Continent, random = ~1|Study/Experiment, data=Data_seed,method = "ML")
anova(Mseed_continent0,Mseed_continent1,Mseed_continent2,Mseed_continent3)
# no difference between fixed-effects and mixed-effects model

# Remit model by REML to compare results -----
Mseed_continent0.1 <- gls(PC~ Continent,data=Data_seed)
summary(Mseed_continent0.1)

Mseed_continent1.1 <- lme(PC~ Continent, random = ~1|Study, data=Data_seed)
summary(Mseed_continent1.1) # between study variability is almost zero

Mseed_continent2.1 <- lme(PC~ Continent, random = ~1|Experiment, data=Data_seed)
summary(Mseed_continent2.1) # between study(experiment) variability is almost zero

# so choose fixed effect model
anova(Mseed_continent0.1) # p = 0.7055, no

# to extract estimate values
mod_seed_continent<-lsmeans(Mseed_continent0.1, ~ Continent)

# model validation for fixed-effects model
par(mfrow = c(1,1),mar=c(6,6,4,4))
hist(resid(Mseed_continent0.1),freq = F)
qqPlot(resid(Mseed_continent0.1))
plot(fitted(Mseed_continent0.1),resid(Mseed_continent0.1))
# histagram and qqplot looks ok

# 4. Does the pollinator contribution is differ between countries?----
Mseed_country0 <- gls(PC~ Country,data=Data_seed,method = "ML")
Mseed_country1 <- lme(PC~ Country, random = ~1|Study, data=Data_seed,method = "ML")
Mseed_country2 <- lme(PC~ Country, random = ~1|Experiment, data=Data_seed,method = "ML")
Mseed_country3<- lme(PC~ Country, random = ~1|Study/Experiment, data=Data_seed,method = "ML")
anova(Mseed_country0,Mseed_country1,Mseed_country2,Mseed_country3) 
# no difference between fixe-effects model and mixed-effects model

# refit model by REML to compare results
Mseed_country0.1 <- gls(PC~ Country,data=Data_seed)
summary(Mseed_country0.1)

Mseed_country1.1 <- lme(PC~ Country,random = ~1|Study,data=Data_seed)
summary(Mseed_country1.1) # between-study variablity is almost zero

# so choose fixed-effects model 
anova(Mseed_country0.1) # p =0.5693, no

# to extract estimate values
mod_seed_country<-lsmeans(Mseed_country0.1, ~ Country) # results same with summary(model)

# model validation
par(mfrow = c(1,1),mar=c(6,6,4,4))
hist(resid(Mseed_country0.1),freq = F)
qqPlot(resid(Mseed_country0.1))
plot(fitted(Mseed_country0.1),resid(Mseed_country0.1))
# histagram and qqplot looks good

# 5. Does the pollinator contribution is differ between climate zones?----
Mclimate_seed0 <- gls(PC~ Climate,data = Data_seed,method = "ML")
Mclimate_seed1 <- lme(PC~ Climate, random = ~1|Study,data = Data_seed,method = "ML")
Mclimate_seed2 <- lme(PC~ Climate, random = ~1|Experiment,data = Data_seed,method = "ML")
Mclimate_seed3 <- lme(PC~ Climate, random = ~1|Study/Experiment,data = Data_seed,method = "ML")
anova(Mclimate_seed0,Mclimate_seed1,Mclimate_seed2,Mclimate_seed3)
# no difference between fixed-effects model and mixed-effects model
# so should choose simpler model, fixed-effects model

# Retit model by REML to compare results----
Mclimate_seed0.1 <- gls(PC~ Climate,data = Data_seed)
summary(Mclimate_seed0.1)

anova(Mclimate_seed0.1) # p = 0.8372
Mclimate_seed1.1 <- lme(PC~ Climate, random = ~1|Study,data = Data_seed)
summary(Mclimate_seed1.1) # between-study variabibility is almost zero

# to extract estimate values from fixed effect model
mod_seed_climate<-lsmeans(Mclimate_seed0.1, ~ Climate)

# to extract estimate values from random effect model
mod_seed_climate_random<-lsmeans(Mclimate_seed1.1, ~ Climate)

# model validation for fixed-effects model
par(mfrow = c(1,1),mar=c(6,6,4,4))
hist(resid(Mclimate_seed0.1),freq = F)
qqPlot(resid(Mclimate_seed0.1)) # looks ok
plot(fitted(Mclimate_seed0.1),resid(Mclimate_seed0.1))

# 6. Does the pollinator contribution is differ between compatibility?----
Mcompa_seed0 <- gls(PC~ Compatibility,data = Data_seed,na.action = na.exclude,method = "ML")
Mcompa_seed1 <- lme(PC~ Compatibility, random = ~1|Study,data = Data_seed,na.action = na.exclude,method = "ML")
Mcompa_seed2 <- lme(PC~ Compatibility, random = ~1|Experiment,data = Data_seed,na.action = na.exclude,method = "ML")
Mcompa_seed3 <- lme(PC~ Compatibility, random = ~1|Study/Experiment,data = Data_seed,na.action = na.exclude,method = "ML")
anova(Mcompa_seed0,Mcompa_seed1,Mcompa_seed2,Mcompa_seed3)
# no difference between fixed-effects model and mixed-effects model
# so should choose fixed-effects model

# Refit model by REML to compare results-----
Mcompa_seed0.1 <- gls(PC~ Compatibility,data = Data_seed,na.action = na.exclude)
summary(Mcompa_seed0.1)
anova(Mcompa_seed0.1) # p = 2e-04, yes, PC is different between compatibility
                      # however, not all cultivars are specified with compatibility type, 
                      # hence we may not stress this result

Mcompa_seed1.1 <- lme(PC~ Compatibility, random = ~1|Study,data = Data_seed,na.action = na.exclude)
summary(Mcompa_seed1.1) # between-study variability is almost zero.

Mcompa_seed2.1 <- lme(PC~ Compatibility, random = ~1|Experiment,data = Data_seed,na.action = na.exclude)
summary(Mcompa_seed2.1) # between-study(experiment) variability is almost zero.

# to extract estimate values from different models
mod_seed_compa<-lsmeans(Mcompa_seed0.1, ~ Compatibility)
multcomp::cld(mod_seed_compa,Letters = c(letters))

mod_seed_compa_random<-lsmeans(Mcompa_seed1.1, ~ Compatibility)
multcomp::cld(mod_seed_compa_random,Letters = c(letters))

# model validation
hist(resid(Mcompa_seed0.1),freq = F)
qqPlot(resid(Mcompa_seed0.1))
plot(fitted(Mcompa_seed0.1),resid(Mcompa_seed0.1))

# 7. Does the pollinator contribution is differ between manipulation levels?----
Mseed_mani0 <- gls(PC~ Manipulation_level,data=Data_seed,method = "ML")
Mseed_mani1 <- lme(PC~ Manipulation_level, random = ~1|Study, data=Data_seed,method = "ML")
Mseed_mani2 <- lme(PC~ Manipulation_level, random = ~1|Experiment, data=Data_seed,method = "ML")
Mseed_mani3 <- lme(PC~ Manipulation_level, random = ~1|Study/Experiment, data=Data_seed,method = "ML")
anova(Mseed_mani0,Mseed_mani1,Mseed_mani2,Mseed_mani3) 
# no difference between fixed-effects model and mixed-effects model
# should choose simple model, fixed-effects model

# Refit the model by REML to compare the results-----
Mseed_mani0.1 <- gls(PC~ Manipulation_level,data=Data_seed)
summary(Mseed_mani0.1)
anova(Mseed_mani0.1) # p = 0.9796, no 

Mseed_mani1.1 <- lme(PC~ Manipulation_level, random = ~1|Study, data=Data_seed)
summary(Mseed_mani1.1) # no between-study vairability is estimated

# to extract estimate values
mod_seed_mani<-lsmeans(Mseed_mani0.1, ~ Manipulation_level)

# model validation
par(mfrow = c(1,1),mar=c(6,6,4,4))
hist(resid(Mseed_mani0.1),freq = F)
qqPlot(resid(Mseed_mani0.1))
plot(fitted(Mseed_mani0.1),resid(Mseed_mani0.1))

#### Interprete the results from mixed-effects model for fruit set data ####
# and results from fixed-effects model for seed set data #
# plot - 1. overall effect----------------------------------------------------------------------------
Mfruit_overall <- lme(PC~ 1,data=Data_fruit,random = ~1|Study,method = "REML")
Mseed_overall <- gls(PC~ 1,data=Data_seed,method = "REML")

summary(Mfruit_overall)  
summary(Mseed_overall) 
intervals(Mfruit_overall) # 0.633  0.712  0.789
intervals(Mseed_overall)  # 0.546  0.619  0.692

# plot overall pollination deficit for fruit set and seed set
PC_overall_fruitseed = data.frame(response = 
                                    c("fruit set","seed set"),
                                  estimate = c(0.712,0.619),
                                  Lower = c(0.633,0.546),
                                  Upper = c(0.789,0.692))

# plot #
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
overall_PC <- ggplot(PC_data, aes(x = Measurement, y = PC))+scale_y_continuous(limits = c(0,1))
overall_PC <- overall_PC + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size =1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),legend.direction = "horizontal",legend.title=element_blank(),
        legend.background=element_blank(),legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "",y="PC", x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5)) 
overall_PC<-overall_PC + geom_pointrange(PC_overall_fruitseed,
                                         mapping = aes(x=response, y=estimate, ymin=Lower, ymax=Upper), 
                                         size=1.3, color=c("springgreen4","blue4"),position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
overall_PC

# Plot- 3.continent -----------------------------------------------------------------------------------------
# Not interprete South amaerica and Oceania beacause of few data points
# exclude observations from two PC database
Data_fruit_continent = Data_fruit[!Data_fruit$Continent == "South America",]
dim(Data_fruit_continent) # 98 observations
Data_seed_continent = Data_seed[!Data_seed$Continent == "Oceania",]
dim(Data_seed_continent) # 61 observations

mod_fruit_conti_random # estimate of PC for fruit set in each continent 
mod_seed_continent

# create dataset for showing PC of continents, 
Continent_fruit_random = data.frame(continent = 
                                      c("Asia","Europe","North-America"),
                                    estimate = c(0.670,0.718,0.741),
                                    Lower = c(0.538,0.596,0.521),
                                    Upper = c(0.801,0.841,0.962))

Continent_seed = data.frame(continent = 
                              c("Asia","Europe","North-America"),
                            estimate = c(0.581,0.665,0.599),
                            Lower = c(0.478,0.544,0.303),
                            Upper = c(0.684,0.786,0.895))

# fruit set data for continent 
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
Data_fruit_continent$Continent = factor(Data_fruit_continent$Continent,
                                        levels = c("Asia","Europe","North-America"))
PC_fruit_continent<- ggplot(Data_fruit_continent, aes(x = Continent, y = PC))+
  scale_y_continuous(limits = c(0,1))
PC_fruit_continent <- PC_fruit_continent + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "fruit set",y="PC", 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PC_fruit_continent<-PC_fruit_continent + 
  geom_pointrange(Continent_fruit_random,
                  mapping = aes(x=continent, y=estimate, ymin=Lower, ymax=Upper),
                  size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

fruit_continent<-PC_fruit_continent

# seed set data for continent
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
Data_seed_continent$Continent = factor(Data_seed_continent$Continent,
                                       levels = c("Asia","Europe","North-America"))
PC_seed_continent <- ggplot(Data_seed_continent, aes(x = Continent, y = PC))+
  scale_y_continuous(limits = c(0,1),position = "right")
PC_seed_continent <- PC_seed_continent + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "seed set",y="PC", 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PC_seed_continent<-PC_seed_continent + 
  geom_pointrange(Continent_seed,
                  mapping = aes(x=continent, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

seed_continent<-PC_seed_continent

# Combine two plot
require(gridExtra)
grid.arrange(fruit_continent, seed_continent, ncol=2)

# Plot- 4.country -----------------------------------------------------------------------------------------
# create dataset for showing PC of countries 
Country_fruit = data.frame(country = 
                             c("India","UK","Belgium","Italy","China","USA","Sweden","Hungary","Argentina"),
                           estimate = c(0.557,0.575,0.630,0.638,0.735,0.736,0.833,0.875,0.998),
                           Lower = c(0.441,0.396,0.466,0.237,0.608,0.565,0.549,0.729,0.597),
                           Upper = c(0.672,0.755,0.793,1.040,0.862,0.908,1.117,1.022,1.40))

Country_seed = data.frame(country = 
                            c("India","UK","Belgium","Italy","China",
                              "USA","Sweden","Pakistan","Turkey","New Zealand"),
                          estimate = c(0.598,0.839,0.584,0.508,0.515,
                                       0.599,0.910,0.565,0.488,0.703),
                          Lower = c(0.477,0.5747,0.4131,0.2120,-0.0769,
                                    0.3030, 0.5682,0.3008,0.147,0.3618),
                          Upper = c(0.719,1.104,0.755,0.804,1.106,
                                    0.895,1.251,0.830,0.830, 1.045))

# fruit set data for cultivar 
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
Data_fruit$Country = factor(Data_fruit$Country,
                            levels = c("India","UK","Belgium","Italy","China","USA","Sweden","Hungary","Argentina"))
PC_fruit_country <- ggplot(Data_fruit, aes(x = Country, y = PC))+
  scale_y_continuous(limits = c(-0.1,1.5))
PC_fruit_country <- PC_fruit_country + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "fruit set",y="PC", 
       x = "",axis.text=element_text(size=11, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(axis.text.x=element_text(angle = -50, vjust = 0.5))
PC_fruit_country<-PC_fruit_country + 
  geom_pointrange(Country_fruit,
                  mapping = aes(x=country, y=estimate, ymin=Lower, ymax=Upper),
                  size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

PC_fruit_country

# seed set data for cultivar
Data_seed$Country = factor(Data_seed$Country ,
                           levels = c("India","UK","Belgium","Italy","China",
                                      "USA","Sweden","Pakistan","Turkey","New Zealand"))
# plot #
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
PC_seed_country <- ggplot(Data_seed, aes(x = Country, y = PC))+
  scale_y_continuous(limits = c(-0.1,1.5))
PC_seed_country <- PC_seed_country + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "seed set",y="PC", 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(axis.text.x=element_text(angle =- 50, vjust = 0.5))
PC_seed_country<-PC_seed_country + 
  geom_pointrange(Country_seed,
                  mapping = aes(x=country, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

PC_seed_country


# Plot- 5.climate -----------------------------------------------------------------------------------------
mod_fruit_climate_random # estimate of PC for fruit set in each climate zone 
mod_seed_climate # estimate of PC for seed set in each climate zone 

# create dataset for showing PC in different climate zone , 
Climate_fruit_random = data.frame(climate = 
                                    c("Tropical","Temperate","Continental","Arid"),
                                  estimate = c(0.623,0.700,0.730,0.831),
                                  Lower = c(0.409,0.603,0.520,0.587),
                                  Upper = c(0.838,0.797,0.939,1.075))

Climate_seed = data.frame(climate = 
                            c("Tropical","Temperate","Continental","Arid"),
                          estimate = c(0.617,0.633,0.524,0.655),
                          Lower = c(0.445,0.541,0.300,0.311),
                          Upper = c(0.788,0.725,0.749,0.998))

# fruit set data for cultivar 

Data_fruit$Climate = factor(Data_fruit$Climate,
                            levels = c("Tropical","Temperate","Continental","Arid"))

# plot #
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
PC_fruit_climate<- ggplot(Data_fruit, aes(x = Climate, y = PC))+
  scale_y_continuous(limits = c(0,1.075))
PC_fruit_climate <- PC_fruit_climate + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "fruit set",y="PC", 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PC_fruit_climate<-PC_fruit_climate + 
  geom_pointrange(Climate_fruit_random,
                  mapping = aes(x=climate, y=estimate, ymin=Lower, ymax=Upper),
                  size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

fruit_climate<-PC_fruit_climate

# seed set data for cultivar
Data_seed$Climate = factor(Data_seed$Climate,
                           levels = c("Tropical","Temperate","Continental","Arid"))
# plot #
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
PC_seed_climate <- ggplot(Data_seed, aes(x = Climate, y = PC))+
  scale_y_continuous(limits = c(0,1.075),position = "right")
PC_seed_climate <- PC_seed_climate + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "seed set",y="PC", 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PC_seed_climate<-PC_seed_climate + 
  geom_pointrange(Climate_seed,
                  mapping = aes(x=climate, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

seed_climate<-PC_seed_climate

# Combine two plots
require(gridExtra)
grid.arrange(fruit_climate, seed_climate, ncol=2)

# Plot- 6.compatibility -----------------------------------------------------------------------------------------
# there are severals observations with no compatibility specified 
# so, create subgroup data set for plotting all original observations
Data_fruit_compa = Data_fruit[!is.na(Data_fruit$Compatibility),]
dim(Data_fruit_compa) # 89 observations
Data_seed_compa = Data_seed[!is.na(Data_seed$Compatibility),]
dim(Data_seed_compa) # 54 observations

mod_fruit_compa_random # estimate of PC for fruit set in each compatibility type
mod_seed_compa # estimate of PC for seed set in each compatibility type 

# create dataset for showing PC of continents
Combatibility_fruit_random = data.frame(compa = 
                                          c("self-compatible","self-incompatible","semi-compatible"),
                                        estimate = c(0.109,0.732,0.628),
                                        Lower = c(-0.136,0.649,0.286),
                                        Upper = c(0.355,0.815,0.970))

Combatibility_seed = data.frame(compa = 
                                  c("self-compatible","self-incompatible","semi-compatible"),
                                estimate = c(0.157,0.661,0.689),
                                Lower = c(-0.057,0.582,0.385),
                                Upper = c(0.371,0.739,0.993))

# fruit set data for cuomaptibility

Data_fruit_compa$Compatibility = factor(Data_fruit_compa$Compatibility,
                                        levels = c("self-compatible","semi-compatible","self-incompatible"))

# plot #
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
PC_fruit_compa<- ggplot(Data_fruit_compa, aes(x = Compatibility, y = PC))+
  scale_y_continuous(limits = c(-0.15,1.05))
PC_fruit_compa <- PC_fruit_compa + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "fruit set",y="PC", 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PC_fruit_compa<-PC_fruit_compa + 
  geom_pointrange(Combatibility_fruit_random,
                  mapping = aes(x=compa, y=estimate, ymin=Lower, ymax=Upper),
                  size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

fruit_compatibility<-PC_fruit_compa

# seed set data for compatibility
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
Data_seed_compa$Compatibility = factor(Data_seed_compa$Compatibility,
                                       levels = c("self-compatible","semi-compatible","self-incompatible"))
PC_seed_compa <- ggplot(Data_seed_compa, aes(x = Compatibility, y = PC))+
  scale_y_continuous(limits = c(-0.15,1.05),position = "right")
PC_seed_compa <- PC_seed_compa + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "seed set",y="PC", 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PC_seed_compa<-PC_seed_compa + 
  geom_pointrange(Combatibility_seed,
                  mapping = aes(x=compa, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

seed_compatibility<-PC_seed_compa

require(gridExtra)
grid.arrange(fruit_compatibility, seed_compatibility, ncol=2)

# Plot- 7.manipulation level -----------------------------------------------------------------------------------------
mod_fruit_mani_random # estimate of PC for fruit set in each manipulation level 
mod_seed_mani # estimate of PC for seed set in each manipulation level 

# create dataset for showing PC of continents
Manipulation_fruit_random = data.frame(manipulation = 
                                         c("Branch","Inflorescence","Plant"),
                                       estimate = c(0.743,0.743,0.578),
                                       Lower = c( 0.642,0.522,0.381),
                                       Upper = c(0.845,0.964,0.776))

Manipulation_seed = data.frame(manipulation = 
                                 c("Branch","Inflorescence","Plant"),
                               estimate = c(0.624,0.610, 0.607),
                               Lower = c(0.536,0.313,0.454 ),
                               Upper = c(0.712,0.906,0.760))

# fruit set data for cuomaptibility

Data_fruit$Manipulation_level = factor(Data_fruit$Manipulation_level,
                                       levels = c("Plant","Branch","Inflorescence"))
# plot #
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
PC_fruit_manipulation<- ggplot(Data_fruit, aes(x = Manipulation_level, y = PC))+
  scale_y_continuous(limits = c(0,1))
PC_fruit_manipulation <- PC_fruit_manipulation + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "fruit set",y="PC", 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PC_fruit_manipulation<-PC_fruit_manipulation + 
  geom_pointrange(Manipulation_fruit_random,
                  mapping = aes(x=manipulation, y=estimate, ymin=Lower, ymax=Upper),
                  size=1.3, color="springgreen4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

fruit_manipulation<-PC_fruit_manipulation

# seed set data for compatibility
Data_seed$Manipulation_level = factor(Data_seed$Manipulation_level,
                                      levels = c("Plant","Branch","Inflorescence"))

# plot #
par(mfrow = c(1,1),mar=c(4.5,4.5,4.5,2.5))
PC_seed_manipulation <- ggplot(Data_seed, aes(x = Manipulation_level, y = PC))+
  scale_y_continuous(limits = c(0,1),position = "right")
PC_seed_manipulation <- PC_seed_manipulation + 
  geom_point(position = position_jitter(w = 0.35, h = 0),size = 1.5,pch = 16,color = "grey68") +
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(legend.position = c(0.2,0.93),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.text = element_text(size = 14),# locate the position of legend
        axis.text=element_text(size=14, color = "black")) + 
  labs(title = "seed set",y="PC", 
       x = "",axis.text=element_text(size=14, color = "black")) + 
  theme(plot.title = element_text(hjust = 0.5))
PC_seed_manipulation<-PC_seed_manipulation + 
  geom_pointrange(Manipulation_seed,
                  mapping = aes(x=manipulation, y=estimate, ymin=Lower, ymax=Upper), 
                  size=1.3, color="blue4",position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0.00),linetype="dashed")+ # add dashed refline 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

seed_manipulation<-PC_seed_manipulation

require(gridExtra)
grid.arrange(fruit_manipulation, seed_manipulation, ncol=2)


# Frequency distribution of PC for fruit set and seed set ----------
par(mfrow = c(1,2),mar=c(4.5,4.5,2.5,2.5))
# fruit set data for PC
hist(Data_fruit$PC,
     main="Fruit set",
     xlab="PC",
     xlim=c(0,1),
     breaks = 20,cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)
text(0.02,20,"(C)",cex = 1.5)

# seed set data for PC
hist(Data_seed$PC,
     main="Seed set",
     xlab="PC",
     xlim=c(0,1),
     breaks = 20,cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)
text(0.02,15,"(D)",cex = 1.5)

#funnel plots of PC for fruit set data and seed set data -----------------------------------
# fruit set data #
# check whether there is any study not reporting number of replicates
length(unique(Data_fruit$Study)) # 23
ID_PC_fruit<-unique(Data_fruit$Study)
ID_PC_fruit

# a vector for mean PC of each study
PC_ave_fruit<-c()

# a vector for study size
SS_PC_fruit<-c()

for (i in 1:28) {
  PC_temp_fruit<-Data_fruit[Data_fruit$Study==ID_PC_fruit[i],]
  PC_ave_fruit<-c(PC_ave_fruit,mean(PC_temp_fruit$PC))
  SS_PC_fruit<-c(SS_PC_fruit,sum(PC_temp_fruit$Nopen))
}

# mean value of PC from mixed effects model is 0.71
# check how many PC above 0.71 
length(which(PC_ave_fruit<0.71))    # 11 values
# ones below 0.71
length(which(PC_ave_fruit>=0.71))    # 12 values

# funnel plot with average PC and Study Size of each study in my seed set dataset
par(mfrow = c(1,2),mar=c(4.5,4.5,4.5,4.5))
plot(PC_ave_fruit,SS_PC_fruit,xlim=c(0.3,1.1),xlab=list("PC (fruit set)",cex=1.5),
     ylim = c(0,150),ylab="Study size",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,tck=-0.02)

# add lines to the funnel plot
# add a vertical line representing the mean of PC 
lines(x=rep(0.71,5),y=seq(0,150,length=5))
# draw a symmetric triangle base on mean PC and maximum PC
# add a line to the left of the vertical line according to the left end of PC and the maximum study size
d<-0.71-0.3
# set the min and max of X for the line
X<-c(0.71-d,0.71)
# set min and max of Y for the line
Y<-c(0,150)
# calculate the slope for the line
lm(Y~X)
# draw the line
X<-seq(0.71-d,0.71,length=5)
Y<--109.8  +  365.9 *X
lines(x=X,y=Y)

# add a line to the right of the vertical line according to the right end of PC and the maximum study size
X<-c(0.71+d,0.71)
Y<-c(0,150)
lm(Y~X)
X<-seq(0.71,0.71+d,length=5)
Y<-  409.8 -365.9 *X
lines(x=X,y=Y)
# add a horizontal line to close the triangle
lines(x=seq(0.71-d,0.71+d,length=5),y=rep(0,5))
text(0.325,150,"(C)",cex=1.5,font=1.5)


# seed set data #
# check whether there is any study not reporting number of replicates
length(unique(Data_seed$Study)) # 17
ID_PC_seed<-unique(Data_seed$Study)
ID_PC_seed

# a vector for mean PC of each study
PC_ave_seed<-c()

# a vector for study size
SS_PC_seed<-c()

for (i in 1:29) {
  PC_temp_seed<-Data_seed[Data_seed$Study==ID_PC_seed[i],]
  PC_ave_seed<-c(PC_ave_seed,mean(PC_temp_seed$PC))
  SS_PC_seed<-c(SS_PC_seed,sum(PC_temp_seed$Nopen))
}

# mean value of PC from fixed-effects model is 0.62
# check how many PC above 0.62 
length(which(PC_ave_seed<0.62))    # 9 values
# ones below 0.62
length(which(PC_ave_seed>=0.62))    # 8 values

# funnel plot with average PC and Study Size of each study in my seed set dataset
plot(PC_ave_seed,SS_PC_seed,xlim=c(0.2,1.05),xlab=list("PC (seed set)",cex=1.5),
     ylim = c(0,150),ylab="Study size",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,tck=-0.02)

# add lines to the funnel plot
# add a vertical line representing the mean of PC 
lines(x=rep(0.62,5),y=seq(0,150,length=5))
# draw a symmetric triangle base on mean PC and maximum PC
# add a line to the left of the vertical line according to the left end of PC and the maximum study size
d<-0.62-0.2
# set the min and max of X for the line
X<-c(0.62-d,0.62)
# set min and max of Y for the line
Y<-c(0,150)
# calculate the slope for the line
lm(Y~X)
# draw the line
X<-seq(0.62-d,0.62,length=5)
Y<--71.43  +  357.14 *X
lines(x=X,y=Y)

# add a line to the right of the vertical line according to the right end of PC and the maximum study size
X<-c(0.62+d,0.62)
Y<-c(0,150)
lm(Y~X)
X<-seq(0.62,0.62+d,length=5)
Y<-  371.4 -357.1 *X
lines(x=X,y=Y)
# add a horizontal line to close the triangle
lines(x=seq(0.62-d,0.62+d,length=5),y=rep(0,5))

text(0.225,150,"(D)",cex=1.5,font=1.5)

# -------------------- end ----------------------------- #