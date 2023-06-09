##############################################
######## nectar bioassay #####################
##############################################

##### Package list ####
library(lme4)
library(nlme)
library(devtools)
require(ggeffects)
library(MASS)
require(shiny)
library(stats)
library(emmeans)
library(blmeco)
library(lattice)
library(car)
library(ggplot2)
library(VGAM)
library(pscl)
library(statmod)
library(plyr)
library(survival)


setwd("/Users/luis/Documents/Research/Tobacco/Data Files")
nect.bio <- read.csv("nectar.bioassay.csv",
                     header = T)
str(nect.bio)


##### collapse Treaments and Period #####
# procedure to collapse treatments to use in analysis with WF included

for(i in 1:length(nect.bio$Treat)){
  if (nect.bio$Treat[i] == "DAM" & nect.bio$Period[i] == "L"){
    nect.bio$collapse[i] <- "DAM-L"
  } else if (nect.bio$Treat[i] == "DAM" & nect.bio$Period[i] == "PRE"){
    nect.bio$collapse[i] <- "DAM-PRE"
  } else if (nect.bio$Treat[i] == "DAM" & nect.bio$Period[i] == "POST"){
    nect.bio$collapse[i] <- "DAM-POST"
  } else if (nect.bio$Treat[i] == "CON" & nect.bio$Period[i] == "L"){
    nect.bio$collapse[i] <- "CON-L"
  } else if (nect.bio$Treat[i] == "CON" & nect.bio$Period[i] == "PRE"){
    nect.bio$collapse[i] <- "CON-PRE"
  } else if (nect.bio$Treat[i] == "CON" & nect.bio$Period[i] == "POST"){
    nect.bio$collapse[i] <- "CON-POST"
  } else if (nect.bio$Treat[i] == "WF"){
    nect.bio$collapse[i] <- "WF"
  }
}


# Change names of period to more intuitive identifiers
for(i in 1:length(nect.bio$Period)){
  if(nect.bio$Period[i]=="PRE"){
    nect.bio$n.Period[i]<-"EARLY"
  } else if (nect.bio$Period[i]=="POST"){
    nect.bio$n.Period[i]<-"MID"
  } else if (nect.bio$Period[i]=="L"){
    nect.bio$n.Period[i]<-"LATE"
  } else if (nect.bio$Period[i]=="WF"){
    nect.bio$n.Period[i]<-"WF"
  }
}


nect.bio$colla <- as.factor(nect.bio$collapse)
nect.bio$n.Period <- as.factor(nect.bio$n.Period)
str(nect.bio)

###### relevel colla factor
nect.bio$colla <- relevel(nect.bio$colla, 
                          "DAM-L")
nect.bio$colla <- relevel(nect.bio$colla, 
                          "CON-L")
nect.bio$colla <- relevel(nect.bio$colla, 
                          "DAM-POST")
nect.bio$colla <- relevel(nect.bio$colla, 
                          "CON-POST")
nect.bio$colla <- relevel(nect.bio$colla, 
                          "DAM-PRE")
nect.bio$colla <- relevel(nect.bio$colla, 
                          "CON-PRE")
nect.bio$colla <- relevel(nect.bio$colla, 
                          "WF")
levels(nect.bio$colla)

##### relevel period factor
nect.bio$n.Period <- relevel(nect.bio$n.Period, 
                             "LATE")
nect.bio$n.Period <- relevel(nect.bio$n.Period, 
                             "MID")
nect.bio$n.Period <- relevel(nect.bio$n.Period, 
                             "EARLY")
levels(nect.bio$n.Period)


par(mfrow=c(1,1))
boxplot(Crithidia~colla, 
        data = nect.bio)
boxplot(Crithidia~Treat, 
        data = nect.bio)
str(nect.bio)

##### relevel Treat factor 
nect.bio$Treat <- relevel(nect.bio$Treat, 
                          "WF")
levels(nect.bio$Treat)

ggplot(nect.bio, aes(Treat, 
                     Crithidia, 
                     colour = n.Period)) + 
  geom_point() +  geom_jitter()

#### Consumption exploration #####
hist(nect.bio$CorrectedConsumption, 
     breaks = 20)
hist(nect.bio$mmwing, 
     breaks = 10)
plot(nect.bio$CorrectedConsumption)
boxplot(nect.bio$CorrectedConsumption~nect.bio$Treat)
boxplot(nect.bio$CorrectedConsumption~nect.bio$Period)
boxplot(nect.bio$CorrectedConsumption~nect.bio$collapse)
x <- ggplot(nect.bio, aes(mmwing, 
                          CorrectedConsumption))
x + geom_point() + geom_smooth(model=lm)

#### Consumption model selection####
full <- lmer(CorrectedConsumption ~ Treat * n.Period * mmwing + (1| Date) +(1|Colony),
           data = nect.bio,
           na.action = na.exclude)
plot(full) ## diagnostic plots look decent, however notice the tails in qqplot
hist(resid(full))
qqPlot(resid(full))
summary(full)

m1a <- lmer(CorrectedConsumption ~ Treat * n.Period * mmwing + (1|Colony),
          data = nect.bio,
          na.action = na.exclude)
m1b <- lmer(CorrectedConsumption ~ Treat * n.Period * mmwing + (1|Date),
            data = nect.bio,
            na.action = na.exclude)

AIC(full, m1a, m1b) ## full best

m2a <- lmer(CorrectedConsumption ~ Treat * n.Period + mmwing + (1|Colony) + (1|Date),
             data = nect.bio,
             na.action = na.exclude)

m2b <- lmer(CorrectedConsumption ~ Treat + n.Period * mmwing + (1|Colony) + (1|Date),
            data = nect.bio,
            na.action = na.exclude)
AIC(full, m2a, m2b) ## m2a best model
summary(m2a)

m3a <- lmer(CorrectedConsumption ~ Treat + n.Period + mmwing + (1|Colony) + (1|Date),
            data = nect.bio,
            na.action = na.exclude)
AIC(m2a, m3a) ## m3a
summary(m3a)
Anova(m3a) ### Get P-values for 
anova(m3a)

m4a <- lmer(CorrectedConsumption ~ n.Period + mmwing + (1|Colony) + (1|Date),
            data = nect.bio,
            na.action = na.exclude)
m4b <- lmer(CorrectedConsumption ~ Treat  + mmwing + (1|Colony) + (1|Date),
            data = nect.bio,
            na.action = na.exclude)
summary(m4b)
m4c <- lmer(CorrectedConsumption ~ Treat + n.Period + (1|Colony) + (1|Date),
            data = nect.bio,
            na.action = na.exclude)
AIC(m3a, m4a, m4b, m4c) # m4b best
summary(m4b)

m5a <- lmer(CorrectedConsumption ~ mmwing + (1|Colony) + (1|Date),
            data = nect.bio,
            na.action = na.exclude)
m5b <- lmer(CorrectedConsumption ~ Treat + (1|Colony) + (1|Date),
                   data = nect.bio,
                   na.action = na.exclude)
AIC(m4b, m5a, m5b) # m5a best
anova(m4b, m5a, m5b)
summary(m5a)

m6 <- lmer(CorrectedConsumption ~ 1 + (1|Colony) + (1|Date),
            data = nect.bio,
            na.action = na.exclude)
AIC(m5a, m6)
summary(m5a) # best model



##### Model validation ####
plot(m5a)
hist(resid(m5a))
qqPlot(resid(m5a))

pred1 <- ggpredict(model = m5a, 
                   terms = "mmwing", 
                   type = "fe") ### 95% confidence interval
f <- ggplot(pred1, 
            aes(x, 
                predicted)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = conf.low, 
                  ymax = conf.high), 
              alpha = .1) + 
  theme_classic() + 
  labs(y = expression(paste("Nectar Consumed (",
                            mu, 
                            "l)")),
       x = "Wing size (mm)")

cor.test(nect.bio$CorrectedConsumption, 
         nect.bio$mmwing, 
         na.action=na.exclude) ### There's correlation between wing and consumption

########################################
#####  Nectar Bioassays Exploration ####
########################################

###### Exploratory plots ########
par(mfrow=c(1,1))
hist(nect.bio$Crithidia, 
     breaks = 10)
plot(table(nect.bio$Crithidia)) # another way to look at it
plot(nect.bio$Crithidia~nect.bio$mmwing)
plot(nect.bio$Crithidia~nect.bio$CorrectedConsumption)

##### Modelling no WF #####
m.pois <- glmer(Crithidia ~ Treat * n.Period + mmwing + (1|Colony) + (1|Date),
              data = nect.bio,
              na.action = na.exclude,
              family = poisson, 
              control = glmerControl(optimizer = "bobyqa", 
                                     optCtrl = list(maxfun=2e5)),
              subset = Treat!="WF")

m.nbfull <- glmer.nb(Crithidia ~ Treat * n.Period + mmwing + (1|Colony) + (1|Date),
                     data = nect.bio,
                     na.action = na.exclude,
                     family = poisson, 
                     control = glmerControl(optimizer = "bobyqa", 
                                            optCtrl = list(maxfun=2e5)),
                     subset = Treat!="WF")

AIC(m.pois, m.nbfull)
summary(m.nbfull)



m.nb1 <- glmer.nb(Crithidia ~ Treat * n.Period + mmwing + (1|Colony),
                  data = nect.bio,
                  na.action = na.exclude,
                  family = poisson, 
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun=2e5)),
                  subset = Treat!="WF")

m.nb2 <- glmer.nb(Crithidia ~ Treat * n.Period + mmwing + (1|Date),
                  data = nect.bio,
                  na.action = na.exclude,
                  family = poisson, 
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun=2e5)),
                  subset = Treat!="WF")

m.nb3 <- glm.nb(Crithidia ~ Treat * n.Period + mmwing,
                  data = nect.bio,
                  na.action = na.exclude,
                  subset = Treat!="WF")

AIC(m.nbfull, m.nb1, m.nb2, m.nb3)
summary(m.nbfull)

m.nb4 <- glmer.nb(Crithidia ~ Treat + n.Period + mmwing + (1|Colony) + (1|Date),
                  data = nect.bio,
                  na.action = na.exclude,
                  family = poisson, 
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun=2e5)),
                  subset = Treat!="WF")

AIC(m.nbfull, m.nb4)

m.nb5 <-  glmer.nb(Crithidia ~ Treat + n.Period + (1|Colony) + (1|Date),
                           data = nect.bio,
                           na.action = na.exclude,
                           family = poisson, 
                           control = glmerControl(optimizer = "bobyqa", 
                                                  optCtrl = list(maxfun=2e5)),
                           subset = Treat!="WF")

m.nb6 <- glmer.nb(Crithidia ~ Treat + mmwing + (1|Colony) + (1|Date),
                  data = nect.bio,
                  na.action = na.exclude,
                  family = poisson, 
                  control = glmerControl(optimizer = "bobyqa", 
                                         optCtrl = list(maxfun=2e5)),
                  subset = Treat!="WF")

m.nb7 <- glmer.nb(Crithidia ~ n.Period + mmwing + (1|Colony) + (1|Date),
                           data = nect.bio,
                           na.action = na.exclude,
                           family = poisson, 
                           control = glmerControl(optimizer = "bobyqa",
                                                  optCtrl = list(maxfun=2e5)),
                           subset = Treat!="WF")
AIC(m.nbfull, m.nb5, m.nb7)

m.nb8 <- glmer.nb(Crithidia ~ n.Period + (1|Colony) + (1|Date),
                           data = nect.bio,
                           na.action = na.exclude,
                           family = poisson, 
                           control = glmerControl(optimizer = "bobyqa",
                                                  optCtrl = list(maxfun=2e5)),
                           subset = Treat!="WF")

m.nb9 <- glmer.nb(Crithidia ~ mmwing + (1|Colony) + (1|Date),
                  data = nect.bio,
                  na.action = na.exclude,
                  family = poisson, 
                  control = glmerControl(optimizer = "bobyqa", 
                                         optCtrl = list(maxfun=2e5)),
                  subset = Treat!="WF")

AIC(m.nb7, m.nb8, m.nb9)
summary(m.nbfull)
summary(m.nb7)

ps <-summary(m.nb7)$coefficients[,4]
p.adjust(ps, 
         method = "fdr") # multiple comparisons adjusted p-value


##### Predictive Plot #####

pred2 <- ggpredict(m.nb7, 
                   terms = c("mmwing", 
                             "n.Period"))

plot(pred2)
summary(pred2)

ggplot(pred2, 
       aes(x, 
           predicted, 
           colour = group), 
       linetype = "solid") + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, 
                  ymax= conf.high), 
              alpha= .2, 
              linetype="dashed") + 
  theme_classic() + 
  labs(y = expression(paste("Crithidia cells/0.02 ", 
                            mu, 
                            "l")), 
       x = "Wing size (mm)") + 
  scale_color_discrete(name="Collection Period") + 
  theme(legend.background = element_rect(size=0.5, 
                                         linetype="solid", 
                                         colour ="black"), 
        axis.text = element_text(size=14), 
        axis.title = element_text(size =18), 
        legend.text = element_text(size=14))


#### Survival analysis ####
smod <- coxph(Surv(d.time, status) ~ Treat + n.Period ,
              subset = Treat!= "WF" & n.Period != "WF",
              data = nect.bio,
              na.action = na.exclude)
summary(smod)

length(nect.bio$Crithidia[nect.bio$Treat=="CON" & nect.bio$Period == "L"])


##############################################################
###### Supplemental analysis of nectar biossays ##############
##############################################################


# This analysis included WF treatment groups
supplement.c <- glmer.nb(Crithidia ~ colla + mmwing + (1|Colony) + (1|Date),
                  data = nect.bio,
                  na.action = na.exclude,
                  family = poisson, 
                  control = glmerControl(optimizer = "bobyqa", 
                                         optCtrl = list(maxfun=2e5)))

hist(resid(supplement.c))
qqPlot(resid(supplement.c))
p.adjust(ps, "fdr")

