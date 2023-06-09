##################################################################################################################

.libPaths("E:/PhD Brown hare/PhD/Statistics/R/R/win-library")


## Start up
#Open project
install.packages("lme4")
library(lme4)
library(colorspace)
library(ggplot2)
library(lattice)
library(MuMIn)
library(AICcmodavg)
library(bbmle)
library(checkmate)
library(Hmisc)
library(lmtest)
library(car)
library(corrgram)
library(arm)
library(zoo)
library(sjPlot)
library(dplyr)
library(stats)
library(coefplot)
library(interplot)



##############################################################################################################
################################## Loading and exploring the data ############################################
##############################################################################################################

## Read in table and explore data: presence absence data
data <- spss.get("E:/PhD Brown hare/PhD/Article Castricum camera trapping/Dryad_repository/Implications of shared predation_database.sav", use.value.labels=TRUE, to.data.frame = TRUE)
attach(data)
summary(data)

## prepare weights for LMER
data$BH.w <- sqrt(BH.cases)
data$R.w <- sqrt(R.cases)
data$RF.w <- sqrt(RF.cases)

# select dataset
#data <- subset(data, BH.time.avg.corr > 0 & RF.time.avg.corr == 0 & R.time.avg.corr == 0)   # set BH 
#data <- subset(data, BH.time.avg.corr > 0 & RF.time.avg.corr == 0 & R.time.avg.corr > 0)    # set BH & R
data <- subset(data, BH.time.avg.corr > 0 & RF.time.avg.corr > 0 & R.time.avg.corr == 0)   # set BH & RF
#data <- subset(data, BH.time.avg.corr > 0 & RF.time.avg.corr > 0 & R.time.avg.corr > 0)     # set BH, R & RF

#data <- subset(data, BH.time.avg.corr == 0 & RF.time.avg.corr == 0 & R.time.avg.corr > 0)   # set R 
#data <- subset(data, BH.time.avg.corr > 0 & RF.time.avg.corr == 0 & R.time.avg.corr > 0)    # set R & BH
#data <- subset(data, BH.time.avg.corr == 0 & RF.time.avg.corr > 0 & R.time.avg.corr > 0)   # set R & RF
#data <- subset(data, BH.time.avg.corr > 0 & RF.time.avg.corr > 0 & R.time.avg.corr > 0)     # set R, BH & RF

#data <- subset(data, Species.present == "BH & RF" | Species.present == "R & RF" | Species.present == "BH & R & RF")    # set RF (with prey)
#data <- subset(data, BH.time.avg.corr == 0 & RF.time.avg.corr > 0 & R.time.avg.corr == 0)   # set RF (without prey)
#data <- subset(data, Species.present == "BH & RF" | Species.present == "BH & R & RF")    # set RF (with hare prey)



## make Area.nr a factor
data$f.Area <- factor(data$Area.nr)
# data$f.Strata <- factor(data$Strata)


## centre session, open.half.open, core.edge
data$Open.Half.open <- as.numeric(data$Open.Half.open)
data$Core.edge <- as.numeric(data$Core.edge)
data$Session <- as.numeric(data$Session.num)

data$c.Open.Half.open <- data$Open.Half.open - mean(data$Open.Half.open)
data$c.Core.edge <- data$Core.edge - mean(data$Core.edge)
data$c.Session <- data$Session - mean(data$Session)








# transformation, centring & z-score of all continous variables
data$l.shrub <- (data$shrub.mean.height)
data$z.shrub <- ((data$l.shrub) - mean(data$l.shrub))/I(2*sd(data$l.shrub))


# average time in plot
data$l.BH.avg <- log10(data$BH.time.avg.corr+0.0000001)
data$z.BH.avg <- ((data$l.BH.avg) - mean(data$l.BH.avg))/I(2*sd(data$l.BH.avg))
data$l.R.avg <- log10(data$R.time.avg.corr+0.0000001)
data$z.R.avg <- ((data$l.R.avg) - mean(data$l.R.avg))/I(2*sd(data$l.R.avg))
data$l.RF.avg <- log10(data$RF.time.avg.corr+0.0000001)
data$z.RF.avg <- ((data$l.RF.avg) - mean(data$l.RF.avg))/I(2*sd(data$l.RF.avg))


## check for distribution & outliers
dotchart(data$z.shrub)
dotchart(data$z.BH.avg)
dotchart(data$z.R.avg)
dotchart(data$z.RF.avg)



# check for (multi)collinearity
## take script from Zuur et al. 2010
source("E:/PhD Brown hare/PhD/Statistics/R/Analysis/Camera_trapping_14_1_2015/Zero_inflated/HighstatLibV6.txt")
## create a dataframe with parameters in the model: individual dataset
par <- data.frame(BH = data$z.BH.avg, R = data$z.R.avg, RF = data$z.RF.avg, shrub = data$z.shrub)


### create correlation tables ratio data
cor(par)

### calculate Variance inflation factor (VIF) for each parameter
corvif(par) 
summary(par)

#############################################################################################################
######################################## LM models #################################################
#############################################################################################################

# LMM with activity = average time in plot (log transformed and z-score) of Brown hare

### full model, no repeated effects, random effect for session
full.model <- lmer(z.BH.avg ~ z.RF.avg + c.Open.Half.open + z.shrub + c.Core.edge + (1|c.Session), data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 



# BH only
model2a <- lmer(z.BH.avg ~ z.shrub + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4a <- lmer(z.BH.avg ~ c.Open.Half.open + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model5a <- lmer(z.BH.avg ~ c.Core.edge + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 


# BH & R only
model2a <- lmer(z.BH.avg ~ z.shrub + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model2b <- lmer(z.BH.avg ~ z.shrub + z.R.avg + z.R.avg*z.shrub + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4a <- lmer(z.BH.avg ~ c.Open.Half.open + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4b <- lmer(z.BH.avg ~ z.R.avg + c.Open.Half.open + z.R.avg*c.Open.Half.open + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model5a <- lmer(z.BH.avg ~ c.Core.edge + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model6a <- lmer(z.BH.avg ~ z.R.avg + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 


# BH & RF only
model2a <- lmer(z.BH.avg ~ z.shrub + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model2c <- lmer(z.BH.avg ~ z.shrub + z.RF.avg + z.RF.avg*z.shrub + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4a <- lmer(z.BH.avg ~ c.Open.Half.open + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4c <- lmer(z.BH.avg ~ z.RF.avg + c.Open.Half.open + z.RF.avg*c.Open.Half.open + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model5a <- lmer(z.BH.avg ~ c.Core.edge + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model7a <- lmer(z.BH.avg ~ z.RF.avg + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 


# BH, R & RF only
model2a <- lmer(z.BH.avg ~ z.shrub + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model2b <- lmer(z.BH.avg ~ z.shrub + z.R.avg + z.R.avg*z.shrub + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model2c <- lmer(z.BH.avg ~ z.shrub + z.RF.avg + z.RF.avg*z.shrub + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4a <- lmer(z.BH.avg ~ c.Open.Half.open + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4b <- lmer(z.BH.avg ~ c.Open.Half.open + z.R.avg + z.R.avg*c.Open.Half.open + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4c <- lmer(z.BH.avg ~ c.Open.Half.open + z.RF.avg + z.RF.avg*c.Open.Half.open + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model5a <- lmer(z.BH.avg ~ c.Core.edge + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model6a <- lmer(z.BH.avg ~ z.R.avg + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model7a <- lmer(z.BH.avg ~ z.RF.avg + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model8a <- lmer(z.BH.avg ~ z.R.avg + z.RF.avg + z.R.avg*z.RF.avg + (1|c.Session), weights=BH.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 


# R only
model2a <- lmer(z.R.avg ~ z.shrub + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4a <- lmer(z.R.avg ~ c.Open.Half.open + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model5a <- lmer(z.R.avg ~ c.Core.edge + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 


# R & BH only
model2a <- lmer(z.R.avg ~ z.shrub + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model2b <- lmer(z.R.avg ~ z.BH.avg + z.shrub + z.BH.avg*z.shrub + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4a <- lmer(z.R.avg ~ c.Open.Half.open + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4b <- lmer(z.R.avg ~ c.Open.Half.open + z.BH.avg + z.BH.avg*c.Open.Half.open + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model5a <- lmer(z.R.avg ~ c.Core.edge + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model6a <- lmer(z.R.avg ~ z.BH.avg + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 


# R & RF only
model2a <- lmer(z.R.avg ~ z.shrub + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model2c <- lmer(z.R.avg ~ z.RF.avg + z.shrub + z.RF.avg*z.shrub + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4a <- lmer(z.R.avg ~ c.Open.Half.open + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4c <- lmer(z.R.avg ~ z.RF.avg + c.Open.Half.open + z.RF.avg*c.Open.Half.open + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model5a <- lmer(z.R.avg ~ c.Core.edge + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model7a <- lmer(z.R.avg ~ z.RF.avg + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 


# R, BH & RF only
model2a <- lmer(z.R.avg ~ z.shrub + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model2b <- lmer(z.R.avg ~ z.BH.avg + z.shrub + z.BH.avg*z.shrub + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model2c <- lmer(z.R.avg ~ z.RF.avg + z.shrub + z.RF.avg*z.shrub + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4a <- lmer(z.R.avg ~ c.Open.Half.open + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4b <- lmer(z.R.avg ~ c.Open.Half.open + z.BH.avg +  z.BH.avg*c.Open.Half.open + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model4c <- lmer(z.R.avg ~ c.Open.Half.open + z.RF.avg + z.RF.avg*c.Open.Half.open + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model5a <- lmer(z.R.avg ~ c.Core.edge + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model6a <- lmer(z.R.avg ~ z.BH.avg + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model7a <- lmer(z.R.avg ~ z.RF.avg + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 
model8a <- lmer(z.R.avg ~ z.RF.avg + z.BH.avg + z.BH.avg*z.RF.avg + (1|c.Session), weights=R.w, data=data, REML=FALSE, control = lmerControl(optimizer = "bobyqa")) 




modelxx <- model4c
confint(modelxx, method = "boot")
summary(modelxx)



# Plots
interplot(m = modelxx, var1 = "c.Open.Half.open", var2 = "z.RF.avg", plot = TRUE, steps = NULL, ci = 0.95, hist = TRUE) +
  # Add labels for X and Y axes
  xlab("Log10 fox average residence time") +
  ylab("Estimated Coefficient for hare average residence time") +
  # Change the background
  theme_bw() +
  # Add the title
  ggtitle("Estimated Coefficient of hare average residence time on average residence time by red fox") +
  theme(plot.title = element_text(face="bold")) +
  # increase size labels x axis
  theme(axis.text.x  = element_text(size=12)) +
  # increase size labels y axis
  theme(axis.text.y = element_text(size=12)) +
  # Add a horizontal line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed")


# barplot z.BH and OHO
library(gplots)
time <- tapply(data$z.BH.avg, data$c.Open.Half.open, mean)
lower <- tapply(data$z.BH.avg, data$c.Open.Half.open, function(v) t.test(v)$conf.int[1])
upper <- tapply(data$z.BH.avg, data$c.Open.Half.open, function(v) t.test(v)$conf.int[2])

barplot2(time, plot.ci = TRUE, ci.l=lower, ci.u=upper,
         ylim=c(-1.5, 2.5), xpd=FALSE,
         main = "",
         names.arg=c("open", "half-open"),
         ylab="Log10 Hare average residence time spent")
abline(h=0) 



# assumptions residuals
#get standardized predicted and residual values
residuals(modelxx)
lres <- residuals(modelxx)
plot(lres)
qqnorm(lres, main = "Normal Q-Q Plot-residuals", plot.it = TRUE, datax = FALSE)
qqline(lres)

predict(modelxx)
lpred <- predict(modelxx) # predicted not in scale response binary variable, see predictions below
plot(lpred)
qqnorm(lpred, main = "Normal Q-Q Plot-predicted", plot.it = TRUE, datax = FALSE)
qqline(lpred)
predictions <- predict(modelxx, type = "response")
data$predictions <- predict(modelxx, type = "response")
plot(predictions)


#create standardized residuals plot
plot(lpred, lres, main= "standardized Residuals Plot", xlab = "Standardized Predicted Values", ylab = "Standardized Residuals") # scatterplot
abline(0,0) # add horizontal line
# above line is underpredicted, below line is overpredicted, close to line is predicted well
# The residuals plot can also be used to test the homogeneity of variance (homoscedasticity ) assumption.
# Look at the vertical scatter at a given point along the x-axis. 
# Now look at the vertical scatter across all points along the x-axis. 
# The homogeneity of variance assumption is supported to the extent that the vertical scatter is the same across all x values.

#create residuals histogram and add normal curve
hist(lres, freq = FALSE, main = "Histogram residuals")
curve(dnorm, add = TRUE)
#create predicted histogram and add normal curve
hist(lpred, freq = FALSE, main = "Histogram predicted")
curve(dnorm, add = TRUE)



predict(full.model)
lpred <- predict(modelxx) # predicted not in scale response binary variable, see predictions below
plot(lpred, data$z.BH.avg, main= "observed vs fitted Plot", xlab = "Fitted Values", ylab = "Observed values") # scatterplot
abline(0,1) # add perfect fit line
rl <- lm(data$z.BH.avg ~ lpred)

y=predict(rl,se=TRUE)
segments(lpred,y$fit+2*y$se.fit,lpred,y$fit-2*y$se.fit,col="blue")
summary(rl)

display(modelxx)
VarCorr(modelxx)
summary(modelxx)
confint(modelxx)
print(modelxx)
ranef(modelxx)
qqmath(ranef(modelxx)) ## plotting random effects
coef(modelxx) # conditional average model
fixef(modelxx) # conditional average model

r.squaredGLMM(modelxx) # R2m is marginal R2 (only fixed effects taken into account), R2c is conditional R2 (fixed and random effects taken into account)
model.set2 <- model.sel(rl)
cor(data$z.BH.avg, lpred, method = c("pearson"))
summary(rl)$r.squared
summary(rl)$adj.r.squared