##### Immunity and co-infection: Field amoebocyte data #######

# Date: 3-11-20
# Allison Tracy

data<-read.csv("Field_amoeb_data_3-11-20.csv")
head(data)


######## Distribution #######
hist(data$amoeb_perc)
shapiro.test(data$amoeb_perc) # not normal

# Log and sqrt transformations don't help 
qqnorm(log(data$amoeb_perc))
qqline(log(data$amoeb_perc)) 

qqnorm(sqrt(data$amoeb_perc))
qqline(sqrt(data$amoeb_perc)) 
shapiro.test(sqrt(data$amoeb_perc))
hist(sqrt(data$amoeb_perc))

### Use Poisson GLM and test for overdispersion 

## Data classes
class(data$colony)
class(data$amoeb_perc) # numeric
class(data$site)  # factor
class(data$temp_mean)  # numeric
class(data$temp_sd) # numeric
data$depth<-as.numeric(data$depth)
class(data$depth) # numeric
class(data$coralcov) # numeric
data$sf_num<-as.numeric(data$sf_num)
class(data$sf_num) # numeric
data$sf_num_sum<-as.numeric(data$sf_num_sum)
class(data$sf_num_sum) # numeric
class(data$cu_all) # factor
data$size<-as.numeric(data$size)
class(data$size) # numeric
class(data$sex) # numeric

# use fan level disease status
class(data$fung_fan) 
class(data$cop_fan) 


## CHECKING data with plots
plot(colony~slide,data=data) # 1 each
plot(transect~site,data=data) # numbering is correct 
plot(temp_mean~site,data=data) # 1 each
plot(temp_sd~temp_mean,data=data) # variance increases with mean
plot(depth~site,data=data) # 1 each - only test depth as a predictor in models without "B" site to evaluate it 
plot(depth~transect,data=data) # 1 each
plot(coralcov~transect,data=data) # 1 each
plot(sf_num~transect,data=data) # 1 each
plot(sf_num_sum~site,data=data) # 1 each
plot(size~colony,data=data) # 1 each
plot(sex~colony,data=data)

plot(fung_fan~slide,data=data)
plot(cop_fan~slide,data=data)

summary(data$fung_fan) # 123
123/918 # 13.4%

### EXPLORATORY plots
# homoscedasticity

plot(amoeb_perc~temp_mean,data=data,ylim=c(0,50)) 
plot(amoeb_perc~temp_sd,data=data,ylim=c(0,50)) # similar to mean so use mean
plot(amoeb_perc~sf_num,data=data)
plot(amoeb_perc~sf_num,data=data,ylim=c(0,50))
plot(amoeb_perc~size,data=data,ylim=c(0,50)) 
plot(amoeb_perc~site,data=data,ylim=c(0,50)) 
plot(amoeb_perc~depth,data=data,ylim=c(0,50)) 
plot(amoeb_perc~coralcov,data=data,ylim=c(0,50))
plot(amoeb_perc~cu_all,data=data,ylim=c(0,50)) 
plot(amoeb_perc~sex,data=data,ylim=c(0,50))
plot(amoeb_perc~fung_fan,data=data,ylim=c(0,50))
plot(amoeb_perc~cop_fan,data=data,ylim=c(0,50))

### Test for collinearity: continuous variables ### 
head(data)
round(cor(data[,c(11,13,14,15,19)]),2) # none over 0.8 - good 

#### RESCALE numeric variables ####
library(arm)
data$temp_mean<-rescale(data$temp_mean, binary.inputs="center")
class(data$temp_mean)
data$depth<-rescale(data$depth, binary.inputs="center")
class(data$depth)
data$coralcov<-rescale(data$coralcov, binary.inputs="center")
class(data$coralcov)
data$sf_num<-rescale(data$sf_num, binary.inputs="center")
class(data$sf_num)
data$size<-rescale(data$size, binary.inputs="center")
class(data$size)

# Contrasts for categorical variables
contrasts(data$cu_all) 
contrasts(data$cu_all) <- cbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1)) # comparing all to level 1
contrasts(data$sex) <- cbind(c(-1,1,0),c(1,1,-2)) # 1st contrast really matters: F vs. M; F+M vs. N contrasts investment 
contrasts(data$sex)
contrasts(data$fung_fan) <-c(-1,1)
contrasts(data$fung_fan) 
contrasts(data$cop_fan) <-c(-1,1)
contrasts(data$cop_fan)


##################
######### Linear models ######### 
##################
# Fixed effects: temp_mean, depth, coral cover, cu_all, sf_num, size, sex, fung_fan, cop_fan 
# Random effect: site
library(lme4) 

null<-glm(amoeb_perc_int~1,data=data,family=poisson)

# additive 
mod1<-glm(amoeb_perc_int~fung_fan+cop_fan+cu_all+depth+sf_num+temp_mean+coralcov+size+sex+site,data=data,family=poisson) # 
drop1(mod1) # drop none

summary(mod1) # very overdispersed- use neg binom
 

###############################
#### Negative binomial 
##############################

library(lme4) 

null<-glm.nb(amoeb_perc_int~1,data=data)

# additive with all 9 factors
newmod1<-glm.nb(amoeb_perc_int~fung_fan+cop_fan+cu_all+sf_num+site+temp_mean+coralcov+size+sex,data=data) # depth not included because tested later in subset (see below)
drop1(newmod1) # drop fungus, copper, sf_num, size
newmod1b<-glm.nb(amoeb_perc_int~cop_fan+site+temp_mean+coralcov+sex,data=data) # not depth 
drop1(newmod1b) # drop temperature
newmod1c<-glm.nb(amoeb_perc_int~cop_fan+site+coralcov+sex,data=data) # not depth 
drop1(newmod1c) # drop none

# test fung * cop interaction: co-infection 
newmod2<-glm.nb(amoeb_perc_int~fung_fan*cop_fan+cu_all+sf_num+site+temp_mean+coralcov+size+sex,data=data) # not depth 
drop1(newmod2) # drop copper, size
newmod2b<-glm.nb(amoeb_perc_int~fung_fan*cop_fan+sf_num+site+temp_mean+coralcov+sex,data=data) # not depth 
drop1(newmod2b) # drop interaction - done

### AIC Table ###
library(AICcmodavg)

AIC(null,newmod1c,newmod2b)

AIC(null) # 18899.54
AIC(newmod1c) # 18767.01 - Best
AIC(newmod2b) # 18771.39

library(lmtest)
lrtest(newmod1c,null) #better than null

# newmod1c is best

# PLOT and summary
summary(newmod1c) # copepod (est = 0.087); sites (greater influence); coralcov (NEW; est=0.18); reproductive - (est= 0.05) 


### VIF for lmer: Collinearity for GLMs ### 
library(rms)
vif(newmod1c) # all acceptable 


# Conf. intervals
library(MuMIn)
confint(newmod1c,level=0.95, method="Wald") 

plot(newmod1c) # good


# Site multcomp - more informative
library(multcomp)
newmod1c_mc<- glht(newmod1c, linfct = mcp(site= "Tukey")) 
summary(newmod1c_mc) 

plot(amoeb_perc_int~fung_fan,data=data)
tapply(data$amoeb_perc,data$fung_fan,mean) 
tapply(data$amoeb_perc,data$cop_fan,mean) # 17.2 without, 19.6 with - yes.
(19.56-17.22)/17.22 # 13.6


################################
######## SUBSET data to test depth #########
################################

### Depth only in models without Buoy: test hypothesis that depth matters in subset of data

no_buoy<-subset(data,site!="B")
summary(no_buoy$site)

length(no_buoy$slide) # 854 vs. 920

######## Distribution #######
hist(no_buoy$amoeb_perc)
shapiro.test(no_buoy$amoeb_perc) # not normal 

# Log and square root not helpful 
qqnorm(log(no_buoy$amoeb_perc))
qqline(log(no_buoy$amoeb_perc)) 

qqnorm(sqrt(no_buoy$amoeb_perc))
qqline(sqrt(no_buoy$amoeb_perc)) 

shapiro.test(sqrt(data$amoeb_perc))
hist(sqrt(data$amoeb_perc))

### Try Poisson GLM

## Data classes
class(no_buoy$colony)
class(no_buoy$amoeb_perc) # numeric
class(no_buoy$site)  # factor
class(no_buoy$temp_mean)  # numeric
class(no_buoy$temp_sd) # numeric
no_buoy$depth<-as.numeric(no_buoy$depth)
class(no_buoy$depth) # numeric
class(no_buoy$coralcov) # numeric
no_buoy$sf_num<-as.numeric(no_buoy$sf_num)
class(no_buoy$sf_num) # numeric
no_buoy$sf_num_sum<-as.numeric(no_buoy$sf_num_sum)
class(no_buoy$sf_num_sum) # numeric
class(no_buoy$cu_all) # factor
no_buoy$size<-as.numeric(no_buoy$size)
class(no_buoy$size) # numeric
class(no_buoy$sex) # numeric

# use fan level disease status
class(no_buoy$fung_fan) 
class(no_buoy$cop_fan) 



### EXPLORATORY plots
# homoscedasticity
plot(amoeb_perc~temp_mean,data=no_buoy,ylim=c(0,50)) 
plot(amoeb_perc~sf_num,data=no_buoy)
plot(amoeb_perc~sf_num,data=no_buoy,ylim=c(0,50))
plot(amoeb_perc~size,data=no_buoy,ylim=c(0,50)) 
plot(amoeb_perc~site,data=no_buoy,ylim=c(0,50)) 
plot(amoeb_perc~depth,data=no_buoy,ylim=c(0,50)) # no longer a major outlier
plot(amoeb_perc~coralcov,data=no_buoy,ylim=c(0,50))
plot(amoeb_perc~cu_all,data=no_buoy,ylim=c(0,50))
plot(amoeb_perc~sex,data=no_buoy,ylim=c(0,50))
plot(amoeb_perc~fung_fan,data=no_buoy,ylim=c(0,50))
plot(amoeb_perc~cop_fan,data=no_buoy,ylim=c(0,50))

### Test for collinearity with correlations: continuous variables ### 
head(no_buoy)
round(cor(no_buoy[,c(11,13,14,15,19)]),2) # - good 

#### RESCALE numeric variables ####
library(arm)
no_buoy$temp_mean<-rescale(no_buoy$temp_mean, binary.inputs="center")
class(no_buoy$temp_mean)
no_buoy$depth<-rescale(no_buoy$depth, binary.inputs="center")
class(no_buoy$depth)
no_buoy$coralcov<-rescale(no_buoy$coralcov, binary.inputs="center")
class(no_buoy$coralcov)
no_buoy$sf_num_sum<-rescale(no_buoy$sf_num_sum, binary.inputs="center")
class(no_buoy$sf_num_sum)
no_buoy$size<-rescale(no_buoy$size, binary.inputs="center")
class(no_buoy$size)


###### LINEAR MODELS ######

# Fixed effects: temp_mean, depth, coral cover, cu_all, sf_num_sum, size, sex, fung_fan, cop_fan (ignore labyweb_fan for now b/c too many groups)
# Random effect: site (not colony b/c already combined)

# Distribution - try poisson first and if overdispersed use negative binomial 

# Test of Poisson
b_depth<-glm(amoeb_perc_int~depth,data=no_buoy,family=poisson) 
summary(b_depth) # way overdispersed


######################
### Negative binomial for overdispersion ## 
######################

library(MASS)
library(lme4) 

bnull<-glm.nb(amoeb_perc_int~1,data=no_buoy)

# additive
bmod1<-glm.nb(amoeb_perc_int~fung_fan+cop_fan+site+coralcov+depth+temp_mean+coralcov+size+sex,data=no_buoy) 
drop1(bmod1) # drop size and fung
bmod1b<-glm.nb(amoeb_perc_int~cop_fan+site+coralcov+depth+temp_mean+coralcov+sex,data=no_buoy) 
drop1(bmod1b) # drop temp
bmod1c<-glm.nb(amoeb_perc_int~cop_fan+site+coralcov+depth+coralcov+sex,data=no_buoy) 
drop1(bmod1c) # drop none

# test fung * cop interaction: co-infection 
bmod2<-glm.nb(amoeb_perc_int~fung_fan*cop_fan+site+coralcov+depth+temp_mean+coralcov+size+sex,data=no_buoy) 
drop1(bmod2) # drop intxn, size - done 


# AIC
AIC(bnull,bmod1d,bmod2) # 1d is best 

# Best model: mod 1d
lrtest(bnull,bmod1d) 
lrtest(bmod4b,bmod4) # additional parameters in 4 don't justify

summary(bmod1d) # depth is still not significant. 
plot(bmod1d) # good



###############
## BAR PLOT 
############


# summarySE - http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/ .....5-23-17
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#### COPEPOD PRESENCE ###
data_summarize1<-summarySE(data,measurevar="amoeb_perc", groupvars=c("cop_fan"))


## Point & LINE ###

pd <- position_dodge(0.5) # if it's smaller, they are more closely grouped by pairs

ggplot(data_summarize1, aes(x=cop_fan, y=amoeb_perc)) + 
  geom_errorbar(aes(ymin=amoeb_perc-1*se, ymax=amoeb_perc+1*se),width=.2, position=position_dodge(0.5)) +
  geom_point(position=pd,size=3) +
  xlab("Copepod presence")+
  ylab("Amoebocyte density (% area)")+
  ylim(c(14.5,26))



### 4 conditions copepod ###

## ADD coinfect column to be able to plot this

data$coinfect <- factor(data$coinfect, levels=c("N","F","C","FC")) # re-orders

data_summarize1<-summarySE(data,measurevar="amoeb_perc", groupvars=c("coinfect"))

pd <- position_dodge(0.5) # if it's smaller, they are more closely grouped by pairs

ggplot(data_summarize1, aes(x=coinfect, y=amoeb_perc)) + 
  geom_errorbar(aes(ymin=amoeb_perc-2*se, ymax=amoeb_perc+2*se),width=.2, position=position_dodge(0.5)) +
  geom_point(position=pd,size=3) +
  xlab("Infection status")+
  ylab("Amoebocyte density (% area)")


### FUNGUS ###

data_summarize1<-summarySE(data,measurevar="amoeb_perc", groupvars=c("fung_fan"))

pd <- position_dodge(0.5) # if it's smaller, they are more closely grouped by pairs

ggplot(data_summarize1, aes(x=fung_fan, y=amoeb_perc)) + 
  geom_errorbar(aes(ymin=amoeb_perc-1*se, ymax=amoeb_perc+1*se),width=.2, position=position_dodge(0.5)) +
  geom_point(position=pd,size=3) +
  xlab("Fungal presence")+
  ylab("Amoebocyte density (% area)")+ 
  ylim(c(14.5,26))


### SEX ###

data_summarize1<-summarySE(data,measurevar="amoeb_perc", groupvars=c("sex"))

pd <- position_dodge(0.5) # if it's smaller, they are more closely grouped by pairs

ggplot(data_summarize1, aes(x=sex, y=amoeb_perc)) + 
  geom_errorbar(aes(ymin=amoeb_perc-1*se, ymax=amoeb_perc+1*se),width=.2, position=position_dodge(0.5)) +
  geom_point(position=pd,size=3) +
  xlab("Repro status")+
  ylab("Amoebocyte density (% area)")+
  ylim(c(14.5,26))


### SITE ###

data_summarize1<-summarySE(data,measurevar="amoeb_perc", groupvars=c("site"))

pd <- position_dodge(0.5) # if it's smaller, they are more closely grouped by pairs

ggplot(data_summarize1, aes(x=site, y=amoeb_perc)) + 
  geom_errorbar(aes(ymin=amoeb_perc-1*se, ymax=amoeb_perc+1*se),width=.2, position=position_dodge(0.5)) +
  geom_point(position=pd,size=3) +
  xlab("Site")+
  ylab("Amoebocyte density (% area)")+
  ylim(c(14.5,26))


