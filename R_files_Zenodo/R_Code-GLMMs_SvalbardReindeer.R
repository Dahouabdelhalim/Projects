####### "Isotope models GLMMS _Svalbard Reindeer - February 2022"
#### Tamara Hiltuen et al

#Import the data 

library(tidyverse)

# import the data. 
mydata <- read.csv("Reindeer_data_schubert_correction.csv", 
                   header = TRUE, 
                   stringsAsFactors = FALSE)

# verify that our data looks correct by printing the first1 10 lines
# to screen
head(mydata)

#mydata$Year.ch <- as.character(mydata$year)
names(mydata)

## Plots and correlation tests
par(mar=c(5,5,3,3))
plot(mydata$year, mydata$d15N, xlab="Year", ylab=(expression(paste(delta^{15}, "N (\\u2030)"))))
cor.test(mydata$year, mydata$d15N)

#Pearson's product-moment correlation
# data:  mydata$year and mydata$d15N
#t = 7.8877, df = 230, p-value = 1.237e-13
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.3536428 0.5571131
#sample estimates:
#      cor 0.4614242 

boxplot(d15N~factor(year, levels=1995:2012), data=mydata, ylab=expression(paste(delta^{15}, "N (\\u2030)")), xlab="", las=1 )
boxplot(d13C~factor(year, levels=1995:2012), data=mydata, ylab=expression(paste(delta^{13}, "C (\\u2030)")), xlab="", las=1 )

ggplot(mydata, aes(x = year, y = d15N)) +
  geom_point() +
  stat_smooth(method = lm)

year.mean<-aggregate(d15N~year, data=mydata, FUN=mean)
cor.test(as.numeric(year.mean$year),year.mean$d15N )

#Pearson's product-moment correlation
#data:  as.numeric(year.mean$year) and year.mean$d15N
#t = 2.6331, df = 14, p-value = 0.01967
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.1116504 0.8334441
#sample estimates:
#     cor 0.575504 

year.mean$d13C<-aggregate(d13C~year, data=mydata, FUN=mean)$d13C
cor.test(as.numeric(year.mean$year),year.mean$d13C )

#Pearson's product-moment correlation
# data:  as.numeric(year.mean$year) and year.mean$d13C
#t = -4.2068, df = 14, p-value = 0.0008789
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.9069734 -0.3994556
#sample estimates:
#      cor -0.7472051 

plot(year.mean$d13C, year.mean$d15N, pch=15, cex=2)
cor.test(year.mean$d13C, year.mean$d15N )

#Pearson's product-moment correlation
#
#data:  year.mean$d13C and year.mean$d15N
#t = -3.1782, df = 14, p-value = 0.006705
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.8653798 -0.2233591
#sample estimates:
#       cor -0.6473857 


#Test for correlations between Year and the Body Mass

plot(mydata$year, mydata$lweight, xlab="Year", ylab= "Body mass (kg)")
cor.test(mydata$year, mydata$lweight)

#Pearson's product-moment correlation
#data:  mydata$year and mydata$lweight
#t = 1.0955, df = 230, p-value = 0.2744
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.05728003  0.19900214
#sample estimates:
#      cor 0.0720502 


#Test for correlations of all variables


### MODELS Nitrogen - d15N

#Load in required packages 

library(mgcv)
library(lme4)
library(cAIC4)
library(lattice)
library(ggplot2)

## Mixed models with id and year as random effect

hist(table(mydata$id))

#some reindeer were capture more than once, some even 4 times.
#Population, ROS, JAT and other variables are year dependant thus it should be a random effect

dN15_lme0<-lm(d15N~factor(year)+lweight, data=mydata)
summary(dN15_lme0)

dN15_lme1<-lmer(d15N~factor(pdprog)+lweight+(1|year)+(1|ID), data=mydata)
summary(dN15_lme1)


dN15_lme1a<-lmer(d15N~population +factor(pdprog)+lweight++(1|year)+(1|ID), data=mydata)
summary(dN15_lme1a)

dN15_lme2<-lmer(d15N~ROS.mm+July.Temp.Average+factor(pdprog)+lweight+(1|year)+(1|ID), data=mydata)
summary(dN15_lme2)


dN15_lme2a<-lmer(d15N~ROS.mm+factor(pdprog)+lweight+(1|year)+(1|ID), data=mydata)
summary(dN15_lme2a)

dN15_lme2b<-lmer(d15N~July.Temp.Average+factor(pdprog)+lweight+(1|year)+(1|ID), data=mydata)
summary(dN15_lme2b)

dN15_lme3<-lmer(d15N~ROS.mm+July.Temp.Average+population +factor(pdprog)+lweight+(1|year)+(1|ID), data=mydata)
summary(dN15_lme3)


###AIC, cAIC and ANOVA best model dEtermination

AIC(dN15_lme1, dN15_lme1a, dN15_lme2, dN15_lme2a, dN15_lme2b, dN15_lme3)
#Best model appears to be s dN15_lme2b then dN15_lme2
# But kest sue the conditional AIC as this is better

cAIC(dN15_lme1)
cAIC(dN15_lme1a)
cAIC(dN15_lme2)
cAIC(dN15_lme2a)
cAIC(dN15_lme2b)
cAIC(dN15_lme3)

##SELECTED MODEL FULL MODEL AS paramater estimates reasonable and the AIC is low >2.

###Checking selected model assumptions

#Tukey-Anscombe plot
plot(dN15_lme3, type = c("p","smooth"), col.line ="black")
#smoother is almost flat on zero (look at scale), hence there is no indication that the model equation may be incorrect.
#The variability seems to be constránt however there is a slight slope


#Scale-location plot
xyplot(sqrt(abs(resid(dN15_lme3))) ~ fitted(dN15_lme3), type = c("p", "smooth", "g"),
       col.line = "black")
#The variance of the residuals does not seem to increase (or decrease) with the fitted values.

#Quantile-Quantile plots
qqnorm(resid(dN15_lme3))
#This QQ-plot shows that there is no clear evidence that the data does not follow a normal distribution.

#QQ-plot for the random effects.
qqmath(ranef(dN15_lme3))

#Although not perfect, this plot does not show a clear violation of the normality assumption for the randomeffects.

## RESIDUALS R2 NITROGEN}

library(MuMIn)
library(insight)
library(performance)
library(piecewiseSEM)

# The following lines:

r.squaredGLMM(dN15_lme1)
r.squaredGLMM(dN15_lme1a)
r.squaredGLMM(dN15_lme2)
r.squaredGLMM(dN15_lme2a)
r.squaredGLMM(dN15_lme2b)
r.squaredGLMM(dN15_lme3)


# The following lines PiecewiseSEM (Nakagawa et al 2017):

rsquared(lme4::lmer(d15N~ROS.mm+July.Temp.Average+population+factor(pdprog)+lweight+(1|year)+(1|ID), data=mydata))


#Response   family     link method  Marginal Conditional
#1     d15N gaussian identity   none 0.3393651   0.7768613

#Estimate marginal and conditional R-squared:
# Marginal = The fixed effects alone explain X% of the variance.
# Conditional = The combined fixed and random effects explain X% of the variance.

#Correlation of residuals N
dN15_lme3<-lmer(d15N~ROS.mm+July.Temp.Average+population+factor(pdprog)+lweight+(1|year)+(1|ID), data=mydata)
summary(dN15_lme3)

dN15_lm1<-lm(d15N~factor(year)+factor(pdprog)+lweight, data=mydata)
summary(dN15_lm1)
anova(dN15_lm1)

dN15_lme4<-lmer(d15N~factor(pdprog)+scale(lweight, scale=F)+(scale(lweight, scale=F)|year), data=mydata)
summary(dN15_lme4)


year.mean$ROS.mm<-aggregate(ROS.mm~year, data=mydata, FUN=mean)$ROS.mm
year.mean$July.Temp.Average<-aggregate(July.Temp.Average~year, data=mydata, FUN=mean)$July.Temp.Average
year.mean$population<-aggregate(population~year, data=mydata, FUN=mean)$population.year
year.mean$lweight<-aggregate(lweight~year, data=mydata, FUN=mean)$lweight

year.mean$log.ROS<-aggregate(log(ROS.mm+1)~year, data=mydata, FUN=mean)[,2]
year.mean$log.population<-aggregate(log(population+1)~year, data=mydata, FUN=mean)[,2]

year.mean$predN<-aggregate(predict(dN15_lme3,re.form=NA)~year, data=mydata, FUN=mean)[,2]

plot(d15N~predN, data=year.mean)
cor.test(year.mean$d15N, year.mean$predN)

cor.test(year.mean$year, year.mean$d15N)
cor.test(year.mean$year, year.mean$predN)

cor.test(year.mean$year, year.mean$d15N-year.mean$predN)

plot(predict(dN15_lme3,re.form =NA)~factor(year, levels=1995:2012), data=mydata)

plot(year.mean$year, ranef(dN15_lme3)$year[,1])

cor.test(year.mean$year, ranef(dN15_lme3)$year[,1])

36+37


### d15N INCREASES WITH BODY MASS, DECREASES IN PREGNANT AND BETWEEN YEAR VARIATION LARGELY (MORE THAN HALF OF VARIATION) EXPLAINED BY ROS AND JULY TEMPERATURE, BOTH HAVING POSITIVE EFFECTS. Population has a very small positive effect.


#Carbon - d13C

#Mixed models with ID random effect D13C


dC13_lm<-lm(d13C~factor(year)+lweight, data=mydata)
summary(dC13_lm)

dC13_lme0<-lmer(d13C~(1|year)+(1|ID), data=mydata)
summary(dC13_lme0)

dC13_lme1<-lmer(d13C~lweight+(1|year)+(1|ID), data=mydata)
summary(dC13_lme1)

dC13_lme2<-lmer(d13C~I(log(ROS.mm+1))+July.Temp.Average+lweight+(1|year)+(1|ID), data=mydata)
summary(dC13_lme2)


dC13_lme2a<-lmer(d13C~I(log(ROS.mm+1))+lweight+(1|year)+(1|ID), data=mydata)
summary(dC13_lme2a)


dC13_lme2b<-lmer(d13C~July.Temp.Average+lweight+(1|year)+(1|ID), data=mydata)
summary(dC13_lme2b)


dC13_lme3<-lmer(d13C~I(log(ROS.mm+1))+July.Temp.Average+I(log(population+1))+lweight+(1|year)+(1|ID), data=mydata)
summary(dC13_lme3)

dC13_lme3a<-lmer(d13C~I(log(population+1))+lweight+(1|year)+(1|ID), data=mydata)
summary(dC13_lme3a)

###AIC, cAIC and ANOVA best model dtermination

AIC(dC13_lme1,dC13_lme2,dC13_lme2a,dC13_lme2b,dC13_lme3,dC13_lme3a)

#Based on the above dC13_lme2a explains more than dC13_lme2b and dC13_lme2.
#But dC13_lme2a explains the most deviance in the data.<


cAIC(dC13_lme1)
cAIC(dC13_lme2)
cAIC(dC13_lme2a)
cAIC(dC13_lme2b)
cAIC(dC13_lme3)
cAIC(dC13_lme3a)

###Checking selected model assumptions

#Tukey-Anscombe plot
plot(dC13_lme3 , type = c("p","smooth"), col.line ="black")
#smoother is almost flat on zero (look at scale), hence there is no indication that the model equation may be incorrect.
#The variability seems to be constránt however there is a slight slope at the lower end


#Scale-location plot
xyplot(sqrt(abs(resid(dC13_lme3))) ~ fitted(dC13_lme3), type = c("p", "smooth", "g"),
       col.line = "black")
#The variance of the residuals does not seem to increase (or decrease) with the fitted values.


#Quantile-Quantile plots
qqnorm(resid(dC13_lme3))

#This QQ-plot shows that there is no clear evidence that the data does not follow a normal distribution.

#QQ-plot for the random effects.
qqmath(ranef(dC13_lme3))

#Although not perfect, this plot does not show a clear violation of the normality assumption for the randomeffects.


##RESIDUALS R2 CARBON}

# The following lines MuMin (Johnson, 2014):

r.squaredGLMM(dC13_lme0)
r.squaredGLMM(dC13_lme1)
r.squaredGLMM(dC13_lme2)
r.squaredGLMM(dC13_lme2a)
r.squaredGLMM(dC13_lme2b)
r.squaredGLMM(dC13_lme3)
r.squaredGLMM(dC13_lme3a)


# The following lines PiecewiseSEM (Nakagawa et al 2017):


rsquared(lme4::lmer(d13C~I(log(ROS.mm+1))+July.Temp.Average+I(log(population+1))+lweight+(1|year)+(1|ID), data=mydata))

#Response   family     link method  Marginal Conditional
#1     d13C gaussian identity   none 0.2202718   0.6456771

# Correlation of residuals c

dC13_lme3<-lmer(d13C~I(log(ROS.mm+1))+July.Temp.Average+I(log(population+1))+lweight+(1|year)+(1|ID), data=mydata)
summary(dC13_lme3)


dC13_lm1<-lmer(d13C~factor(year)+lweight+(1|ID), data=mydata)
summary(dC13_lm1)
anova(dC13_lm1)


dC13_lme4<-lmer(d13C~ scale(lweight, scale=F) + (scale(lweight, scale = F)|year), data=mydata)
summary(dC13_lme4)

year.mean$ROS.mm<-aggregate(ROS.mm~year, data=mydata, FUN=mean)$ROS.mm
year.mean$July.Temp.Average<-aggregate(July.Temp.Average~year, data=mydata, FUN=mean)$July.Temp.Average
year.mean$population.year<-aggregate(population~year, data=mydata, FUN=mean)$population.year
year.mean$lweight<-aggregate(lweight~year, data=mydata, FUN=mean)$lweight

year.mean$log.ROS<-aggregate(log(ROS.mm+1)~year, data=mydata, FUN=mean)[,2]
year.mean$log.population.year<-aggregate(log(population+1)~year, data=mydata, FUN=mean)[,2]

year.mean$predC<-aggregate(predict(dC13_lme3,re.form=NA)~year, data=mydata, FUN=mean)[,2]

plot(d13C~predC, data=year.mean)
cor.test(year.mean$d13C, year.mean$predC)

cor.test(year.mean$year, year.mean$d13C)
cor.test(year.mean$year, year.mean$predC)

cor.test(year.mean$year, year.mean$d13C-year.mean$predC)

plot(predict(dC13_lme3,re.form =NA)~factor(year, levels=1995:2012), data=mydata)

plot(year.mean$year, ranef(dC13_lme3)$year[,1])

cor.test(year.mean$year, ranef(dC13_lme3)$year[,1])

45+45

## d13C is decreasing with body mass and between year variation to some degree explained by log(ROS+1) (13 % of variation) with a negative effect
# However, still a negative temporal trend in d13C after controlling for ROS.