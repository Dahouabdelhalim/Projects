# Article: "Evidence for Individual discrimination and numerical assessment in collective antipredator behaviour in jackdaws (Corvus monedula)"
# Authors: Jenny R. Coomes, Guillam E. McIvor, Alex Thornton
# Script: Recruits to playbacks


#install.packages("ggeffects")
#install.packages("piecewiseSEM")
#install.packages("blmeco")
#install.packages("influence.ME")
# install.packages("ggplot2")
# install.packages("ICC")
library(lme4)
library(MuMIn)
library(influence.ME)
library(ggplot2)
library(pastecs)
library(tidyverse)
library(tidyr)
library(dplyr)
library(knitr)
library(ICC)
library(MCMCglmm)

jd.rec<-read.csv("PB.recruits.csv",header=T,na.strings="NA")
head(jd.rec)
dim(jd.rec)
sapply(jd.rec, class)

noscold<-subset(jd.rec,scold=="no")
dim(noscold)

# Analyses
Fullmod<-glmer(recruits ~ treatment + scold + order + num.date + wind + time 
               + (1|colony/location), family = "poisson", control=glmerControl(optimizer="bobyqa"), data = jd.rec)

# Model checks
DR <- residuals(Fullmod,type="deviance")
plot(DR~fitted(Fullmod)) # Not wonderful
qqnorm(DR)
qqline(DR)
hist(DR,10)
plot(jd.rec$wind, DR, xlab="wind", ylab="Deviance Residuals")
plot(jd.rec$time, DR, xlab="time", ylab="Deviance Residuals")
plot(jd.rec$treatment, DR, xlab="treatment", ylab="Deviance Residuals")
plot(jd.rec$order, DR, xlab="order", ylab="Deviance Residuals")
plot(jd.rec$scold, DR, xlab="scold", ylab="Deviance Residuals")
plot(jd.rec$colony, DR, xlab="site", ylab="Deviance Residuals")
plot(jd.rec$location, DR, xlab="location", ylab="Deviance Residuals") # Generally ok
# Examine Cooks distances to identify cases with high influence
Fullmod.infl <- influence(Fullmod, obs = TRUE)
cooks.distance(Fullmod.infl)
plot(Fullmod.infl, which = "cook")
# we do have two highly influential cases (Cooks distance >1), rows 34 and 44. Re-run the top models without these rows
# to ensure they are not the driving force behind the relationships our analysis identifies

# Check for overdispersion - aim for a value of 1
(Fullmod <- deviance(Fullmod) / df.residual(Fullmod)) # = 4.07 so overdispersed

# We will take and IT approach to model selection, following the methods of Richards et al in the papers listed below
# https://www.jstor.org/stable/20143968  |  https://link.springer.com/article/10.1007/s00265-010-1035-8
# We will use QAICc and a quasi correction, as suggested by Ben Bolker  https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#fitting-models-with-overdispersion
# to account for overdispersion in the data. 

# We also have explored the alternative option of fitting observation level random effects as advocated by Harrison 2015
# (https://peerj.com/articles/1114/). This results in an underdispersed model, and so problems of a different kind that 
# must be accounted for (https://doi.org/10.1080/03610926.2017.1291976). Applying the same quasi-correction suggested by Bolker
# (https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#underdispersion) we end up with fairly identical models as we
# obtain with the corrected overdispersed models. However comparison of model plots to observation level random effects being
# a less suitable approach for our data

# Maximum number of terms in the model limited to 4 to avoid over-parameterisation. All models must include two terms 
# (treatment and scold) that are key to the experimental predictions. Models ranked by QAICc due to overdispersion
chat <- deviance(Fullmod) / df.residual(Fullmod)

Fullmod.dredge<-dredge(Fullmod, m.lim=c(2,4), chat=chat, fixed=c("treatment","scold"), rank="QAICc",k=2) 

Fullmod.dredge
Fullmod.nest<-subset(Fullmod.dredge,df<=7)
Fullmod.nest
Fullmod.dredge.d6<-subset(Fullmod.nest,delta<6) # call only those in the top set for final model weights calculation
Fullmod.dredge.d6

# Write both tables
#write.table(Fullmod.nest, "~/Datasets/Coomes et al - model selection table post nested.csv", sep=",")
#write.table(Fullmod.dredge.d6, "~/Datasets/Coomes et al - final model weights.csv", sep=",")



# Obtain the summaries for each of the retained models
# R-squared values can only be obtained prior to the quasi correction, but this should not be a problem given that they
# are calculated based on how well the model fits the data, and the estimates do not change as a result of the correction,
# only the SE, Z scores, and p Values. 
model3<-(get.models(Fullmod.dredge.d6, 1)[[1]])
summary(model3)

# overdispersion function from (https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion)
# Required for all the upcoming overdispersion corrections. Does not need to be re-run for each model
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)}

r.squaredGLMM(model3)
cc3 <- coef(summary(model3))
phi <- overdisp_fun(model3)["ratio"]
cc3 <- within(as.data.frame(cc3),
             {   `Std. Error` <- `Std. Error`*sqrt(phi)
             `z value` <- Estimate/`Std. Error`
             `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
             })
printCoefmat(cc3,digits=3)
# This makes our results more conservative, but the key results to the manuscript remain

# Also Call the model summaries of model1
model1<-(get.models(Fullmod.dredge.d6, 2)[[1]])
summary(model1)
r.squaredGLMM(model1)
cc1 <- coef(summary(model1))
phi <- overdisp_fun(model1)["ratio"]
cc1 <- within(as.data.frame(cc1),
             {   `Std. Error` <- `Std. Error`*sqrt(phi)
             `z value` <- Estimate/`Std. Error`
             `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
             })
printCoefmat(cc1,digits=3)

# Also Call the model summaries of model2

model2<-(get.models(Fullmod.dredge.d6, 3)[[1]])
summary(model2)
r.squaredGLMM(model2)

cc2 <- coef(summary(model2))
phi <- overdisp_fun(model2)["ratio"]
cc2 <- within(as.data.frame(cc2),
             {   `Std. Error` <- `Std. Error`*sqrt(phi)
             `z value` <- Estimate/`Std. Error`
             `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
             })
printCoefmat(cc2,digits=3)

# Finally, rerun model 3 but with the cases highlighted as influential above (rows 34 and 44) removed
cooks.model3<-glmer(recruits ~ treatment + scold + order  
                 + (1|colony/location), family = "poisson", control=glmerControl(optimizer="bobyqa"), data = jd.rec[-c(34,44),])
summary(cooks.model3)

cooks.mod3.correct <- coef(summary(cooks.model3))
phi <- overdisp_fun(cooks.model3)["ratio"]
cooks.mod3.correct <- within(as.data.frame(cooks.mod3.correct),
                            {   `Std. Error` <- `Std. Error`*sqrt(phi)
                            `z value` <- Estimate/`Std. Error`
                            `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                            })
printCoefmat(cooks.mod3.correct,digits=3)
# The outcome is essentially the same - removing the potentially influential cases does not alter our findings



### Test using only cases without responsive scolding, to determine how robust the experimental treatment effects are
nsc.mod<-glmer(recruits ~ treatment + order 
               + (1|colony/location), family = "poisson", control=glmerControl(optimizer="bobyqa"), data = noscold)
summary(nsc.mod)
r.squaredGLMM(nsc.mod)

nosco<- coef(summary(nsc.mod))
phi <- overdisp_fun(nsc.mod)["ratio"]
nosco <- within(as.data.frame(nosco),
              {   `Std. Error` <- `Std. Error`*sqrt(phi)
              `z value` <- Estimate/`Std. Error`
              `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
              })
printCoefmat(nosco,digits=3)


# Model checks
DR.nsc <- residuals(nsc.mod,type="deviance")
plot(DR.nsc~fitted(nsc.mod)) # Heteroscedasticity, probably due to few playbacks with large numbers?
qqnorm(DR.nsc)
qqline(DR.nsc)
hist(DR.nsc,10)

nsc.infl <- influence(nsc.mod, obs = TRUE) # calculate cooks distances to identify potentially influential cases
cooks.distance(nsc.infl)
plot(nsc.infl, which = "cook")  # row 22 is potentially influential

# re-run with row 22 excluded
nsc.cook<-glmer(recruits ~ treatment + order 
               + (1|colony/location), family = "poisson", control=glmerControl(optimizer="bobyqa"), data = noscold[-c(22),])
summary(nsc.cook)

nsc.cook.corr<- coef(summary(nsc.cook))
phi <- overdisp_fun(nsc.cook)["ratio"]
nsc.cook.corr <- within(as.data.frame(nsc.cook.corr),
                {   `Std. Error` <- `Std. Error`*sqrt(phi)
                `z value` <- Estimate/`Std. Error`
                `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                })
printCoefmat(nsc.cook.corr,digits=3)
# Removal of influential case does have some influence on model, but general predictions hold true
# treatment effect is stronger when cases of responsive scolding are removed, suggesting a masking effect of responsive 
# scolding with regards our experimental predictions. This to to be expected, as responsive scolding can substantially 
# alter the number of callers, and is not possible to consistently quantify


########      Post hoc tests    ########
# With responsive scolding
treatmentGS3GS5sco<- subset(jd.rec,treatment == "GS3" | treatment == "GS5")
dim(treatmentGS3GS5sco)

GS3GS5<-glmer(recruits ~ treatment + scold + order + 
               + (1|colony/location), family = "poisson", control=glmerControl(optimizer="bobyqa"), data = treatmentGS3GS5sco)
summary(GS3GS5)
(GS3GS5.disp <- deviance(GS3GS5) / df.residual(GS3GS5)) # = 2.87 so overdispersed
GS3GS5.quasi <- coef(summary(GS3GS5))
phi <- overdisp_fun(GS3GS5)["ratio"]
GS3GS5.quasi <- within(as.data.frame(GS3GS5.quasi),
                          {   `Std. Error` <- `Std. Error`*sqrt(phi)
                          `z value` <- Estimate/`Std. Error`
                          `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                          })
printCoefmat(GS3GS5.quasi,digits=3) # No difference between GS3 and GS5 when cases of responsive scolding are included

# Check for influential cases
GS3GS5.infl <- influence(GS3GS5, obs = TRUE) # calculate cooks distances to identify potentially influential cases
cooks.distance(GS3GS5.infl)
plot(GS3GS5.infl, which = "cook")  # 4 influential cases

# re-run with rows 5,6,29,30 excluded
GS3GS5.cook<-glmer(recruits ~ treatment + order + scold
                + (1|colony/location), family = "poisson", control=glmerControl(optimizer="bobyqa"), data = treatmentGS3GS5sco[-c(5,6,29,30),])
summary(GS3GS5.cook)
GS3GS5.sco.cook <- coef(summary(GS3GS5.cook))
phi <- overdisp_fun(GS3GS5.cook)["ratio"]
GS3GS5.sco.cook <- within(as.data.frame(GS3GS5.sco.cook),
                      {   `Std. Error` <- `Std. Error`*sqrt(phi)
                      `z value` <- Estimate/`Std. Error`
                      `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                      })
printCoefmat(GS3GS5.sco.cook,digits=3) # results consistent even when influential cases removed
# No difference between GS3 and GS5 when cases of responsive scolding are included



### Post-hoc GS3 GS5 Without repsonsive scolding
treatmentGS3GS5nosco<- subset(noscold, treatment == "GS3" | treatment == "GS5")
dim(treatmentGS3GS5nosco)

GS3GS5nosco<-glmer(recruits ~ treatment + order + 
                + (1|colony/location), family = "poisson", control=glmerControl(optimizer="bobyqa"), data = treatmentGS3GS5nosco)
summary(GS3GS5nosco)
(GS3GS5nosco.disp <- deviance(GS3GS5nosco) / df.residual(GS3GS5nosco)) # 2.50 so overdispersed
# Overdispersion correction to model summary
GS3GS5quasi <- coef(summary(GS3GS5nosco))
phi <- overdisp_fun(GS3GS5nosco)["ratio"]
GS3GS5quasi <- within(as.data.frame(GS3GS5quasi),
                      {   `Std. Error` <- `Std. Error`*sqrt(phi)
                      `z value` <- Estimate/`Std. Error`
                      `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                      })
printCoefmat(GS3GS5quasi,digits=3)

GS3GS5nosco.infl <- influence(GS3GS5nosco, obs = TRUE) # calculate cooks distances to identify potentially influential cases
cooks.distance(GS3GS5nosco.infl)
plot(GS3GS5nosco.infl, which = "cook")  # No influential data points - is fine as is



# Obtain decriptive stats 
GS1<-subset(jd.rec,treatment=="GS1")
stat.desc(GS1) # median = 2.5, mean (s.e.) = 11.31 (4.42) n =16
GS3<-subset(jd.rec,treatment=="GS3")
stat.desc(GS3) # median = 6.5 mean (s.e.) = 10.31 (2.43) n =16
GS5<-subset(jd.rec,treatment=="GS5")
stat.desc(GS5) # median = 5  mean (s.e.) = 14.25 (4.42) n =16

GS1.nosco<-subset(noscold,treatment=="GS1")
stat.desc(GS1.nosco) # median = 1.5, mean (s.e.) = 2.5 (0.93)  , n = 12
GS3.nosco<-subset(noscold,treatment=="GS3")
stat.desc(GS3.nosco) # median = 5, mean (s.e.) = 4.73 (1.14), n =11
GS5.nosco<-subset(noscold,treatment=="GS5")
stat.desc(GS5.nosco) # median = 3, mean(s.e.) = 7.17 (2.61), n = 12

### In response to the comments of R5 - use MCMCglmm to take into account the possibility that some individuals may elicit 
# a greater response. As they state: At present the model does not take into account the possibility that some individual 
# jackdaws elicit a greater response. The issue with this is that as you include more individuals you are more likely to 
# include the jackdaw(s) that elicit this response. An additional posthoc test that  you might consider in order to estimate
# and control for such an effect would be to employ a multi-membership random term, so that you then estimate how the 
# response varies depending on whether different individuals are present or absent. As far as I am aware a multi-membership
# design is not implemented in lme4, so I think you would need to use MCMCglmm. The syntax would be something like this and 
# you would need a column in your dataset for each individual bird and whether it was present or absent in a recording.

# Run it for the possible top models. Use the default settings for priors

model1.MCMC<-MCMCglmm(recruits~scold+treatment, random=~colony+location+idv(~y29.ind+y11.ind+y04.ind+y07.ind+y19.ind+
                                                                              y06.ind+y10.ind+y04u.ind+y01.ind+z30.ind+z15.ind+z26.ind+z20.ind+z14.ind+z22.ind+z18.ind+z19.ind+
                                                                              z28.ind+z43.ind+z02.ind+z24.ind),family="poisson",data=jd.rec,nitt=200000,burnin=20000)
summary(model1.MCMC)

model1.without<-MCMCglmm(recruits~scold+treatment, random=~colony+location+idv,family="poisson",data=jd.rec,nitt=200000,burnin=20000)
summary(model1.without)


model2.MCMC<-MCMCglmm(recruits~scold+treatment+num.date, random=~colony+location+idv(~y29.ind+y11.ind+y04.ind+y07.ind+y19.ind+
                                                                                       y06.ind+y10.ind+y04u.ind+y01.ind+z30.ind+z15.ind+z26.ind+z20.ind+z14.ind+z22.ind+z18.ind+z19.ind+
                                                                                       z28.ind+z43.ind+z02.ind+z24.ind),family="poisson",data=jd.rec,nitt=200000,burnin=20000)
summary(model2.MCMC)

model2.without<-MCMCglmm(recruits~scold+treatment+num.date, random=~colony+location,family="poisson",data=jd.rec,nitt=200000,burnin=20000)
summary(model2.without)

model3.MCMC<-MCMCglmm(recruits~scold+treatment+order, random=~colony+location+idv(~y29.ind+y11.ind+y04.ind+y07.ind+y19.ind+
                                                                                    y06.ind+y10.ind+y04u.ind+y01.ind+z30.ind+z15.ind+z26.ind+z20.ind+z14.ind+z22.ind+z18.ind+z19.ind+
                                                                                    z28.ind+z43.ind+z02.ind+z24.ind),family="poisson",data=jd.rec,nitt=200000,burnin=20000)
summary(model3.MCMC)

model3.without<-MCMCglmm(recruits~scold+treatment+date, random=~colony+location,family="poisson",data=jd.rec,nitt=200000,burnin=20000)
summary(model3.without)

# Finally considering only trials with no responsive scolding
nsc.MCMC1<-MCMCglmm(recruits~treatment+order, random=~colony+location+idv(~y29.ind+y11.ind+y04.ind+y07.ind+y19.ind+
                                                                            y06.ind+y10.ind+y04u.ind+y01.ind+z30.ind+z15.ind+z26.ind+z20.ind+z14.ind+z22.ind+z18.ind+z19.ind+
                                                                            z28.ind+z43.ind+z02.ind+z24.ind) ,family="poisson", data=noscold,nitt=200000,burnin=20000)
summary(nsc.MCMC1)
nsc.without<-MCMCglmm(recruits~treatment+order, random=~colony+location ,family="poisson", data=noscold,nitt=200000,burnin=20000)
summary(nsc.without)

# No substantial effect of the idv multimembership random term, all are not sig.dif from zero. 
# Therefore there seems to be no effect of individual caller on the number of recruits




### Intercoder reliability ###
library(ICC)
intercod<-read.csv("2019Intercoder.reliability.csv",header=T)
ICCest(sd.recruits,jc.recruits,data=intercod)
icc(intercod, model=c("twoway"), type=c("agreement"), unit=c("average"))



### Finally for the Supp Mat - reanalyse excluding cases where the same calls were used multiple times
no.rep<-read.csv("~/Datasets/2019.Withoutrepeatcalling.csv",header=T,na.strings="NA")
head(no.rep)
dim(no.rep)

n.repMod<-glmer(recruits ~ treatment + scold + order + num.date + wind + time 
                + (1|colony/location), family = "poisson", control=glmerControl(optimizer="bobyqa"), data = no.rep)
summary(n.repMod)

(chat <- deviance(n.repMod) / df.residual(n.repMod)) # 3.64 -> overdispersed
nrep.dredge<-dredge(n.repMod, m.lim=c(2,4), chat=chat, fixed=c("treatment","scold"), rank="QAICc",k=2) 

nrep.dredge
nrep.nest<-subset(nrep.dredge,df<=7)
nrep.nest
# All models retained once the nesting rule is applied have delta QAICc <6
# The two best supported models are identical to those from our main analysis. 
# Two additional models also attract weak support - one containing 'wind' as an explanatory variable, and 
# one containing time of day that the experiment was performed.

# Call the top two model summaries to ensure the results are still consitent after applying the quasi correction
modelR3<-(get.models(nrep.nest, 1)[[1]])
summary(modelR3)
r.squaredGLMM(modelR3)

ccR3 <- coef(summary(modelR3))
phi <- overdisp_fun(modelR3)["ratio"]
ccR3 <- within(as.data.frame(ccR3),
               {   `Std. Error` <- `Std. Error`*sqrt(phi)
               `z value` <- Estimate/`Std. Error`
               `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
               })
printCoefmat(ccR3,digits=3)


# And for model 1
modelR1<-(get.models(nrep.nest, 2)[[1]])
summary(modelR1)
r.squaredGLMM(modelR1)

ccR1 <- coef(summary(modelR1))
phi <- overdisp_fun(modelR1)["ratio"]
ccR1 <- within(as.data.frame(ccR1),
               {   `Std. Error` <- `Std. Error`*sqrt(phi)
               `z value` <- Estimate/`Std. Error`
               `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
               })
printCoefmat(ccR1,digits=3)

# And finally consider the subset without responsive scolding
nosc.nrep<-subset(no.rep,scold=="no")
dim(nosc.nrep)
nosc.nrep.MOD<-glmer(recruits ~ treatment + order + 
                       + (1|colony/location), family = "poisson", control=glmerControl(optimizer="bobyqa"), data = nosc.nrep)
summary(nosc.nrep.MOD)
r.squaredGLMM(nosc.nrep.MOD)
cc.nosc.rep.1 <- coef(summary(nosc.nrep.MOD))
phi <- overdisp_fun(nosc.nrep.MOD)["ratio"]
cc.nosc.rep.1<- within(as.data.frame(cc.nosc.rep.1),
                       {   `Std. Error` <- `Std. Error`*sqrt(phi)
                       `z value` <- Estimate/`Std. Error`
                       `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                       })
printCoefmat(cc.nosc.rep.1,digits=3)

# So pretty consistent, order drops out. Suggests that rather than what the reviewers thought that by chance some
# of the calls repeated may have been particularly stimulating, if anything they are less stimulating and/or birds may
# remember them from last time and be less likely to recruit


GS3GS5.nrep<-subset(nosc.nrep,treatment == "GS3" | treatment == "GS5")
dim(GS3GS5.nrep)
final.MOD<-glmer(recruits ~ treatment + order + 
                   + (1|colony/location), family = "poisson", control=glmerControl(optimizer="bobyqa"), data = GS3GS5.nrep)
summary(final.MOD) # only 20 observations - we are fast running out of data
r.squaredGLMM(final.MOD)
cc.final.MOD<- coef(summary(final.MOD))
phi <- overdisp_fun(final.MOD)["ratio"]
cc.final.MOD<- within(as.data.frame(cc.final.MOD),
                      {   `Std. Error` <- `Std. Error`*sqrt(phi)
                      `z value` <- Estimate/`Std. Error`
                      `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                      })
printCoefmat(cc.final.MOD,digits=3)
# But no change in results still, but could easily be because we have run out of statistical power





#### Graphs ####  

### New Figure 1a
# Use housemouse code to make graphs clearer: http://rpubs.com/tomhouslay/beyond-bar-line-graphs
jd.rec$treatment<-as.factor(jd.rec$treatment)
jd.rec$recruits<-as.numeric(jd.rec$recruits)
jd.rec$GS1[jd.rec$treatment=="GS1"]<- "1"
jd.rec$GS1[jd.rec$treatment=="GS3"]<- "3"
jd.rec$GS1[jd.rec$treatment=="GS5"]<- "5"
jd.rec$GS1<-as.factor(jd.rec$GS1)
std_err <- function(x){ sd(x) / sqrt(length(x))}
df_ind_sum <- jd.rec %>% 
  group_by(GS1) %>% 
  summarise(Grp_mean = mean(recruits),
            Grp_se = std_err(recruits))
Coomes.Fig1a<-ggplot(jd.rec, aes(x = GS1, y = recruits)) +
  labs(x="Number of Callers",y="Recruits")+  ylim(-1,61)+
  geom_point(size = 3, position = position_jitter(width = 0.2),
             alpha = 0.7) +
  geom_pointrange(data = df_ind_sum,
                  aes(y = Grp_mean,
                      ymin = Grp_mean - Grp_se, 
                      ymax = Grp_mean + Grp_se),
                  colour = "red",
                  alpha = 0.7,
                  size = 1) + theme_bw() +
  theme(text = element_text(size=18,face= "bold",colour="black"))+
  theme(axis.text=element_text(size=16,colour="black"))
Coomes.Fig1a



### New Figure 1b
noscold<-subset(jd.rec,scold=="no")
dim(noscold)
df_ind_sum <- noscold %>% 
  group_by(GS1) %>% 
  summarise(Grp_mean = mean(recruits),
            Grp_se = std_err(recruits))
Coomes.Fig1b<-ggplot(noscold, aes(x = GS1, y = recruits)) +
  labs(x="Number of Callers",y="Recruits")+  ylim(-1,61)+
  geom_point(size = 3, position = position_jitter(width = 0.2),
             alpha = 0.7) +
  geom_pointrange(data = df_ind_sum,
                  aes(y = Grp_mean,
                      ymin = Grp_mean - Grp_se, 
                      ymax = Grp_mean + Grp_se),
                  colour = "red",
                  alpha = 0.7,
                  size = 1) + theme_bw() +
  theme(text = element_text(size=18,face= "bold",colour="black"))+
  theme(axis.text=element_text(size=16,colour="black"))
Coomes.Fig1b



### Supp Mat Fig1

indeff<-read.csv("~/Datasets/indiveffect.csv",header=T,na.strings="NA")
head(indeff)
std_err <- function(x){ sd(x) / sqrt(length(x))}

# Plot variation by individual
df_ind_sum <- indeff %>% 
  group_by(indiv) %>% 
  summarise(Grp_mean = mean(recruits),
            Grp_se = std_err(recruits)) 

palette<-c("red2","mediumseagreen","dodgerblue")
Supp.Mat.S1<-ggplot(indeff, aes(x = indiv, y = recruits,colour=treatment)) +
  scale_fill_manual(values=palette)+scale_colour_manual(values=palette,name="Treatment")+
  labs(x="Individual",y="Recruits")+  ylim(-1,60)+
  geom_point(size = 3, position = position_jitter(width = 0.2),
             alpha = 0.7) +
  geom_pointrange(data = df_ind_sum,
                  aes(y = Grp_mean,
                      ymin = Grp_mean - Grp_se, 
                      ymax = Grp_mean + Grp_se),
                  colour = "magenta",
                  alpha = 0.7,
                  size = 1) + theme_bw() +
  theme(text = element_text(size=18,face= "bold",colour="black"))+
  theme(axis.text=element_text(size=16,colour="black"))
Supp.Mat.S1












