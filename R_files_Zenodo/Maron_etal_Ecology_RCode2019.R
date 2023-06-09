## This file contains the R Code to conduct the main analyses presented in 
 ## Maron et al. (2019, Ecology) paper entitled:
 ## "Seedling recruitment correlates with seed input across seed sizes:
 ##   implications for coexistence".
 ## Please contact Phil Hahn (hahnp@ufl.edu) with any questions. 


library(nlstools)
library(bbmle)
library(nlme)
library(FSA)
library(lmtest)
library(lme4)
library(lmerTest)
library(emmeans)
library(tidyverse)
library(minpack.lm)

rec <-read_csv("Maron_etal_Ecology_DATA2019.csv")

#look at data
ggplot(rec %>% filter(Year==2017), aes(x=TotalSeedlings)) +
  geom_histogram() +
  facet_wrap(~Species, scales="free_x")


#######################################################################################
### H1 SKELLAM RECRUITMENT CURVES #####################################################
#######################################################################################

##################################################################################
#### SKELLAM EXAMPLE 1 ###########################################################
### Code below is an example for one species and treatment combination only.
### Example is for a model that does not propertly estimate the 'n' parameter and refit with linear model.

### SUBSET DATA BY SPECIES, RODENT =/-, AND YEAR TO BUILD EACH MODEL SEPERATELY
ls1 <- filter(rec, Species=='ACMI', Rodents=='-', Competition=='+', Year==2017) # %>% filter(!seedinput==0|TotalSeedlings<40)

#PLOT DATA
plot(TotalSeedlings ~ seedinput, data=ls1)

## SET VAGUE STARTING PARAMETERS AND ESTIMATE PARAMETERS
## USE 'nlsLM' from package:minpack.lm TO ESTIMATE STARTING PARAMETERS
b_start=.01
n_start=10
ls1.skA <- nlsLM(TotalSeedlings ~ b*n*(1-(exp(-seedinput/n))),
              data=ls1, start=list(b=b_start,n=n_start))
summary(ls1.skA)



## UPDATE STARTING PARAMETERS FROM ABOVE AND ESTIMATE PARAMETERS 
b_start=.039
n_start=50000000
ls1.sk <- nls(TotalSeedlings ~ b*n*(1-(exp(-seedinput/n))),
                data=ls1, start=list(b=b_start,n=n_start))
summary(ls1.sk) ## model cannot estimate 'n' parameter

## REDUCE MODEL TO LINEAR
b_start1=.039
ls1.ln <- nls(TotalSeedlings ~ b*seedinput,
              data=ls1, start=list(b=b_start1))
summary(ls1.ln) ### estimate parameters (see Table 1 for parameters for all species,treatment combos)

plot(nlsResiduals(ls1.ln)) ## examine residuals. Most look okay. Variance increases at higher fitted values,
                            # but residuals typically centered around zero so models seem to estimate reasonable parameters

plot(TotalSeedlings ~ seedinput, data=ls1, main=" 2017 +C-R")
plotfit(ls1.ln, smooth=T) ## plot fitted model
ax <- seq(0,40000,1)
ay1 <- coef(ls1.ln)[1]*ax ## predict linear function
ay3 <- coef(ls1.sk)[1]*coef(ls1.sk)[2]*(1-(exp(-ax/coef(ls1.sk)[2])) ) ## predict skellam function, NA here
lines(ax,ay1, lty=1,lwd=2,col="red") ## plot linear function (thick red)
lines(ax,ay3, lty=1,lwd=2,col="blue") ## plot skellam function (blue), NA here

###############################################################################################
####### SKELLAM EXAMPLE 2 #####################################################################
### Code below is an example for one species and treatment combination only.
### Example is for a model that does estimate the 'n' parameter.

### SUBSET DATA BY SPECIES, RODENT =/-, AND YEAR TO BUILD EACH MODEL SEPERATELY
ls2 <- filter(rec, Species=='LUSE', Rodents=='-', Competition=='+', Year==2017) # %>% filter(!seedinput==0|TotalSeedlings<40)

#PLOT DATA
plot(TotalSeedlings ~ seedinput, data=ls2)

## SET VAGUE STARTING PARAMETERS AND ESTIMATE PARAMETERS
## USE 'nlsLM' from package:minpack.lm TO ESTIMATE STARTING PARAMETERS
b_start=.01
n_start=10
ls2.skA <- nlsLM(TotalSeedlings ~ b*n*(1-(exp(-seedinput/n))),
                 data=ls2, start=list(b=b_start,n=n_start))
summary(ls2.skA)


## UPDATE STARTING PARAMETERS FROM ABOVE AND ESTIMATE PARAMETERS 
b_start=.66
n_start=51
ls2.sk <- nls(TotalSeedlings ~ b*n*(1-(exp(-seedinput/n))),
              data=ls2, start=list(b=b_start,n=n_start))
summary(ls2.sk) ## estimate parameter (see Table 1 for all parameters)


plot(nlsResiduals(ls2.sk)) ## examine residuals. Most look okay. Variance increases at higher fitted values,
                # but residuals typically centered around zero so models seem to estimate reasonable parameters

summary(nlsBoot(ls2.sk)) ## check estimate using nlsJack and nlsBoot from package:nlstools
                        ## parameter estimates should be close to 'summary()'


## REDUCE MODEL TO LINEAR, just for comparison
b_start1=.66
ls2.ln <- nls(TotalSeedlings ~ b*seedinput,
              data=ls2, start=list(b=b_start1))
summary(ls2.ln)


## PLOT DATA AND FIT LINEAR AND SKELLAM FOR COMPARISON
plot(TotalSeedlings ~ seedinput, data=ls2, main=" 2017 +C-R")
bx <- seq(0,400,1)
by1 <- coef(ls2.ln)[1]*bx
by3 <- coef(ls2.sk)[1]*coef(ls2.sk)[2]*(1-(exp(-bx/coef(ls2.sk)[2])) )
lines(bx,by1, lty=1,lwd=2,col="red") ## plot linear fit (red)
lines(bx,by3, lty=1,lwd=2,col="blue") ## plot saturating fit (blue)



###########################################################################
###### HYPOTHESIS 2: PCA ANALYSES ##########################################
############################################################################

## subset and summarize data for testing hypothesis 2
loop1 <- rec %>% filter(Competition=='+'&Rodents=='-', SeedAddition>.1, Species!='LIRU') %>% 
  group_by(Species) %>% mutate(proprec1=scale(log(proprec+.001)), proprec2=scale(log((proprec+.001)/(1-(proprec+.001)))))  %>% 
  group_by(Site,Species,Year) %>%
  summarize(proprec1=mean(proprec1, na.rm=T), PC1=mean(PC1), PC2=mean(PC2), proprec=mean(proprec, na.rm=T), proprec2=mean(proprec2, na.rm=T))
hist(loop1$proprec1) ## check different transformations, use proprec1
hist(loop1$proprec)
hist(loop1$proprec2)


lpc1 <- lmer(proprec1~Species*PC1+Species*PC2+(1|Site/Species), data=loop1 )
#lpc1 <- lmer(proprec1~Species*PC1*as.factor(Year)+Species*PC2+(1|Site/Species), data=loop1 )

summary(lpc1)
anova(lpc1)
emtrends(lpc1, ~Species, var="PC1")
emtrends(lpc1, ~Species, var="PC2")

ggplot(loop1, aes(x=PC1, y=proprec1))+geom_point()+geom_smooth(method='lm')+facet_wrap(~Species)
ggplot(loop1, aes(x=PC2, y=proprec1))+geom_point()+geom_smooth(method='lm')+facet_wrap(~Species)



################################################################################
#### HYPOTHESIS 3: SITE YEAR RODENT SITE ########################################
################################################################################
loop2 <- rec %>% filter(Competition=='+', SeedAddition>.1, Species!='LIRU') %>% 
  group_by(Species) %>% mutate(proprec1=log(proprec+.001))  %>% group_by(Site,Species,Year,Rodents) %>%
  summarize(proprec1=mean(proprec1, na.rm=T))

aov1 <- lmer(proprec1~as.factor(Year)*Species*Rodents+
               (1|Site/Rodents/Species),
             data=loop2)
anova(aov1)
ranova(aov1)
summary(aov1)
plot(aov1)

emmeans(aov1, pairwise~Year|Species)
emmeans(aov1, pairwise~Rodents|Species)

