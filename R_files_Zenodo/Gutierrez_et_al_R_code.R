# R Code accompanying Gutiérrez et al.
# Royal Society Open Science 2020

# required R packages:
# lme4, car, coxme, survival and lavaan

# if needed, install these by running the following lines:
#install.packages("lme4")
#install.packages("car")
#install.packages("coxme")
#install.packages("survival")
#install.packages("lavaan")

# load data
setwd("C:\\\\...")
survival.analysis=read.table("survival analysis.txt",h=T,sep="\\t",na.strings = "NA")
development.time=read.table("development time.txt",h=T,sep="\\t",na.strings = "NA")
weight.adults=read.table("weight adults.txt",h=T,sep="\\t",na.strings = "NA")
consumption=read.table("food consumption.txt",h=T,sep="\\t",na.strings = "NA")
cann.prop.tot=read.table("cannibalism proportions.total.txt",h=T,sep="\\t",na.strings = "NA")
protein=read.table("protein - bradford.txt",h=T,sep="\\t",na.strings = "NA")
lipids=read.table("lipids.total.txt",h=T,sep="\\t",na.strings = "NA")
egg.production=read.table("eggs.time.temp.txt",h=T,sep="\\t",na.strings = "NA")
total.eggs=read.table("eggs.total.temp.txt",h=T,sep="\\t",na.strings = "NA")
lifespan=read.table("lifespan.txt",h=T,sep="\\t",na.strings = "NA")
femdata=read.table("SEMfemales.txt",h=T,sep="\\t")
maledata=read.table("SEMmales.txt",h=T,sep="\\t")


# load R packages
library(nlme)
library(lme4)
library(car)
library(coxme)
library(survival)
library(lavaan)


## Mixed-effects survival analysis
model1 <- coxme(Surv(day,status)~group*diet+(1|exp_unit),data=survival.analysis)
summary(model1)
Anova(model1)

## GLMM for development time
model2 <- glmer(adulthood ~ group*diet*sex + (1|exp_unit), family = poisson,  data = development.time)
summary(model2)
Anova(model2)

## LME for individual body mass (weight)
# sexes combined
model3 <- lme(weight ~ diet*group*sex, random=~1|exp_unit, data = weight.adults)
summary(model3)
Anova(model3)

# males only
adult_split=split(weight.adults,weight.adults$sex)

model3.males <- lme(weight ~ diet*group, random=~1|exp_unit, data = adult_split$male)
summary(model3.males)
Anova(model3.males)

# females only
model3.females <- lme(weight ~ diet*group, random=~1|exp_unit, data = adult_split$female)
summary(model3.females)
Anova(model3.females)

## GLM for food consumption
# model comparing: single_fem, single_male and grouped(male+fem)
model4.a <- glm(food_per_capita ~ diet*sex, data = consumption, family=Gamma(link ="inverse"))
summary(model4.a)
Anova(model4.a)

# model comparing: solitary(male+fem) and grouped(male+fem)
model4.b <- glm(food_per_capita ~ diet*group, data = consumption, family=Gamma(link ="inverse"))
summary(model4.b)
Anova(model4.b)

## GLM for cannibalism
model5 <- glm(proportion~diet, family = quasibinomial(link="logit"), data=cann.prop.tot)
summary(model5)
Anova(model5)

## GLM for protein content
# total protein
model6.a <- glm(prot_total_mg~sex*diet*group,data=protein, family = Gamma(link="inverse"))
summary(model6.a)
Anova(model6.a)

# relative protein content
model6.b <- glm(protconc_mgmg~sex*diet*group,data=protein, family=quasibinomial(link="logit"))
summary(model6.b)
Anova(model6.b)

## GLM for lipid content
# total lipids
model7.a=glm(total_lipids~sex*diet*group,data=lipids, family = Gamma(link="inverse"))
summary(model7.a)
Anova(model7.a)

# relative lipid content
model7.b=glm(relativ_lipid~sex*diet*group,data=lipids, family=quasibinomial(link="logit"))
summary(model7.b)
Anova(model7.b)

## GLMM for weekly egg production
model8 <- glmer(eggs~ poly(week,3)*diet*treatment + (1|code), family= negative.binomial(1), data=egg.production)
summary(model8)
Anova(model8)

## GLM for total egg production
model9 <- glm(total_eggs~diet*treatment, data=total.eggs, family=negative.binomial(1))
summary(model9)
Anova(model9)

## survival analysis for lifespan
model10 = survfit(Surv(lifespan.w, status) ~ diet, data = lifespan)
summary(model10)
survdiff(Surv(lifespan.w, status) ~ diet, data = lifespan)

## structural equation models (SEM)

females="food_per_capita ~ group 
weight_adult ~ group + diet + food_per_capita
adulthood ~ food_per_capita + diet
weight_adult ~~ adulthood"
females1=lavaan(females, fixed.x=T, std.ov=T, int.ov.free=T, auto.var=T, auto.cov.y=T, test="standard", estimator="MLMVS", data=femdata)
summary(females1)

males="food_per_capita ~ group
weight_adult ~ group + food_per_capita
adulthood ~ diet
weight_adult ~~ adulthood"
males1=lavaan(males, fixed.x=T, std.ov=T, int.ov.free=T, auto.var=T, auto.cov.y=T, test="standard", estimator="MLMVS", data=maledata)
summary(males1)




