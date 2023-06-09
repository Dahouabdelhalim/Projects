rm(list=ls())

# load required packages
library(readr)
library(plyr)
library(lme4)
library(multcomp)
library(glmmADMB)
library(bbmle)
library(lme4)
library(phia)
library(sinaplot)
library(glmmTMB)

# function to ceck for overdispersion 
overdisp_fun <- function(model) {
  ## number of variance parameters in
  ## an n-by-n variance-covariance matrix
  vpars <- function(m)
  { nrow(m)*(nrow(m)+1)/2 }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  Pearson.chisq <- sum(residuals(model, type="pearson")^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(Pearson_dev=Pearson.chisq,rdf=rdf,ratio=prat,p=pval)
}

#drop1
LRT_fun<-function(model1,model2,df=1)
  
{LRT<-2*(summary(model1)$loglik-summary(model2)$loglik)
p_val<-pchisq(LRT,df,lower.tail=F)
return(list(LRT=LRT,p=p_val))

}

#Load data
setwd("")
predation_final <- read.csv("Data.csv", header=T, sep=";", na.strings=c("NA",""))

########################################################################################
##################overall: Model 1

final_treat <- with(predation_final, data.frame(
  territory = territory,
  activity_TF = scale(activity_TF,scale=F)[,1],
  treatment = treatment,
  aggression = aggression,
  harem = harem,
  stimulus_id = stimulus_id,
  Time_attack=Time_attack/60,
  B_H = B_H,
  id=id))


final_treat <- na.omit(final_treat)
str(final_treat)
final_treat$territory<-as.factor(final_treat$territory)
final_treat$aggression<-as.numeric(final_treat$aggression)
final_treat$harem<-as.factor(final_treat$harem)

final_treat$Time_attack[final_treat$Time_attack<0] <- 0

#Run models with different distributions
treat_poi0 <- glmmadmb(aggression ~ B_H*treatment+activity_TF+offset(log(Time_attack))+(1|stimulus_id)+
                         (1|harem), data=final_treat, family="poisson", 
                       zeroInflation=T)

treat_poi <- glmmadmb(aggression ~ B_H*treatment+activity_TF+offset(log(Time_attack))+(1|stimulus_id)+
                        (1|harem), data=final_treat, family="poisson", 
                      zeroInflation=F)

treat_nbi <- glmmadmb(aggression ~ B_H*treatment+activity_TF+offset(log(Time_attack))+(1|stimulus_id)+
                        (1|harem), data=final_treat, family="nbinom", 
                      zeroInflation=F)

###############
#Best model
treat_nbi0 <- glmmadmb(aggression ~ B_H*treatment+activity_TF+offset(log(Time_attack))+(1|stimulus_id)+
                         (1|harem), data=final_treat, family="nbinom", 
                       zeroInflation=T)



AICtab(treat_poi0,treat_poi,treat_nbi,treat_nbi0)

overdisp_fun(treat_nbi0)

summary(treat_nbi0)
drop1(treat_nbi0, test="Chisq")

#post-hoc test of interaction
testInteractions(treat_nbi0, pairwise="B_H", fixed="treatment", adjustment="none")
testInteractions(treat_nbi0, pairwise="treatment", fixed="B_H",adjustment="none")

#sample sizes
new<-data.frame(ddply(final_treat, .(treatment,B_H), summarize,N=length(aggression), size=mean(aggression),
                      sd=sd(aggression),se=sd/sqrt(N)))

#################################################################################################
#Figure 1a
social <- with(predation_final, data.frame(
  territory = territory,
  treatment = treatment,
  bre_hel = bre_hel,
  aggression = aggression,
  young_yn = young_yn,
  activity_TF=scale(activity_TF, scale=F)[,1],
  stimulus_id=stimulus_id,
  B_H=B_H,
  helper_yn = helper_yn,
  gro_tre=gro_tre,
  harem = harem,
  id=id))

social<-na.omit(social)
str(social)
social$territory<-as.factor(social$territory)
social$harem<-as.factor(social$harem)
social$young_yn<-as.factor(social$young_yn)


new<-data.frame(ddply(social, .(treatment,B_H), summarize,N=length(id), agg=mean(aggression),
                      sd=sd(aggression),se=sd/sqrt(N)))

sochelper<-subset(social, B_H=="H")


new3<-data.frame(ddply(sochelper, .(gro_tre), summarize,  agg=sum(aggression)))
finalH<-merge(new3, sochelper[!duplicated(sochelper$gro_tre),], by=("gro_tre"), all.x=T)

breeder<-subset(social, B_H=="B")
new3<-data.frame(ddply(breeder, .(gro_tre), summarize,  agg=sum(aggression)))
finalB<-merge(new3, breeder[!duplicated(breeder$gro_tre),], by=("gro_tre"), all.x=T)

#combine breeder and helper files
final<-rbind(finalB, finalH)

finalB<-subset(final, B_H=="B")
str(finalB)

#Plot 1a
finalB$pred_treat<-revalue(finalB$treatment, c("pfe"="2","elo"="1","vit"="3",
                                               "spi"="4"))
finalB$pred_treat<-as.character(finalB$pred_treat)
boxplot(agg~pred_treat, data=finalB,xlab="Type of predator", 
        ylab="Number of attacks by breeders",cex.lab=1.2, las=1, cex.axis = 2,
        xlim=c(0.8,4.2),ylim=c(0,320), yaxt="n", xaxt="n", boxwex=0.2, outline=F,
        whisklty = 1, whisklwd=2,boxfill="lightblue", medcol="red", staplelwd=2)


legend(3.5, 335, legend ="(a)", bty="n", cex=1.5)

axis(side = 1, at = seq(1, 4, by = 1),labels=c("Fish","Young", "Eggs", "Control"),
     cex.axis=1.2)
axis(side = 2, at = seq(0, 300, by = 50),cex.axis=1.2, las=1)
legend(0.6, 355, legend = "a", bty = "n", cex = 1.3, text.col="red")
legend(1.6, 160, legend = "b", bty = "n", cex = 1.3, text.col="red")
legend(2.6, 207, legend = "b", bty = "n", cex = 1.3, text.col="red")
legend(3.6, 102, legend = "c", bty = "n", cex = 1.3, text.col="red")

test<-sinaplot(finalB$agg~finalB$pred_treat, method="density", maxwidth=0.6,scale=F, alpha=0.4, plot=F)

points(y~scaled, data=test,lwd=1)
#Export 5.24x4.63

#Figure 1b

finalH<-subset(final, B_H=="H")

finalH$pred_treat<-revalue(finalH$treatment, c("pfe"="2","elo"="1","vit"="3",
                                               "spi"="4"))
finalH$pred_treat<-as.character(finalH$pred_treat)

boxplot(agg~pred_treat, data=finalH,xlab="Type of predator", 
        ylab="Number of attacks by helpers",cex.lab=1.2, las=1, cex.axis = 2,
        xlim=c(0.8,4.2),ylim=c(0,150), yaxt="n", xaxt="n", boxwex=0.2, outline=F,
        whisklty = 1, whisklwd=2,boxfill="lightblue", medcol="red", staplelwd=2)


legend(3.5, 165, legend ="(b)", bty="n", cex=1.5)

axis(side = 1, at = seq(1, 4, by = 1),labels=c("Fish","Young", "Eggs", "Control"),
     cex.axis=1.2)
axis(side = 2, at = seq(0, 140, by = 20),cex.axis=1.2, las=1)
legend(0.6, 45, legend = "a", bty = "n", cex = 1.3, text.col="red")
legend(1.5, 80, legend = "a,b", bty = "n", cex = 1.3, text.col="red")
legend(2.6, 160, legend = "b", bty = "n", cex = 1.3, text.col="red")
legend(3.6, 45, legend = "c", bty = "n", cex = 1.3, text.col="red")

test<-sinaplot(finalH$agg~finalH$pred_treat, method="density", maxwidth=0.6,scale=F, alpha=0.4, plot=F)

points(y~scaled, data=test,lwd=1)
#Export 5.24x4.63


#####################################################################
#Breeder (model 2)
######################################################################

breeder <- subset(predation_final,B_H=="B")

final_treat <- with(breeder, data.frame(
  territory = territory,
  activity_TF = scale(activity_TF,scale=F)[,1],
  treatment = treatment,
  sex = social_rank,
  aggression = aggression,
  helper_yn = helper_yn,
  young_yn = young_yn,
  harem = harem,
  stimulus_id = stimulus_id,
  fish_size=fish_sl,
  Time_attack=Time_attack/60))


final_treat <- na.omit(final_treat)
str(final_treat)


final_treat$territory<-as.factor(final_treat$territory)
final_treat$treatment<-as.factor(final_treat$treatment)
final_treat$sex<-as.factor(final_treat$sex)
final_treat$young_yn<-as.factor(final_treat$young_yn)
final_treat$helper_yn<-as.factor(final_treat$helper_yn)
final_treat$harem<-as.factor(final_treat$harem)
final_treat$aggression<-as.numeric(final_treat$aggression)


#size range of breeders within sex is very small --> Exlusion

new<-data.frame(ddply(final_treat, .(sex), summarize,N=length(fish_size), size=mean(fish_size),
                      sd=sd(fish_size),se=sd/sqrt(N)))


new<-data.frame(ddply(final_treat, .(treatment,helper_yn), summarize,N=length(aggression), size=mean(aggression),
                      sd=sd(aggression),se=sd/sqrt(N)))


#Run models to compare distribution

treat2_poiZi <- glmmadmb(aggression ~ treatment*helper_yn+young_yn+sex+activity_TF+offset(log(Time_attack))+
                          (1|harem)+(1|stimulus_id), data=final_treat, 
                        family="poisson", zeroInflation = T)

########################

treat2_poi <- glmmadmb(aggression ~ treatment*helper_yn+young_yn+sex+activity_TF+offset(log(Time_attack))+
                        (1|harem)+(1|stimulus_id), data=final_treat, 
                      family="poisson", zeroInflation = F)
########################

treat2_nbi <- glmmadmb(aggression ~ treatment*helper_yn+young_yn+sex+activity_TF+offset(log(Time_attack))+
                        (1|harem)+(1|stimulus_id), data=final_treat, 
                      family="nbinom", zeroInflation = F)

treat2_nbiZi <- glmmadmb(aggression ~ treatment*helper_yn+young_yn+activity_TF+sex+offset(log(Time_attack))+
                          (1|harem)+(1|stimulus_id), data=final_treat, 
                        family="nbinom", zeroInflation = T)


AICtab(treat2_poiZi,treat2_poi,treat2_nbi,treat2_nbiZi)

#treat_nbiZi has best fit


overdisp_fun(treat_nbiZi)
summary(treat_nbiZi)
Anova(treat_nbiZi,type=3)


testInteractions(treat_nbiZi, pairwise="helper_yn", fixed="treatment", adjustment="none")

#####################

#test for the interction between treatment and young
treat_nbiZi_young <- glmmadmb(aggression ~helper_yn+treatment*young_yn+activity_TF+sex+offset(log(Time_attack))+
                          (1|harem)+(1|stimulus_id), data=final_treat, 
                        family="nbinom", zeroInflation = F)

drop1(treat_nbiZi_young, test="Chisq")


##########################################################################
#Plot Figure 2
#Plot added values of breeder females and breeder males

boxplot(agg~young_yn, data=finalB,xlab="Presence of young", 
        ylab="Number of attacks by breeders",cex.lab=1.2, las=1, cex.axis = 2,
        xlim=c(0.5,2.5),ylim=c(0,300), yaxt="n", xaxt="n", boxwex=0.2, outline=F,
        whisklty = 1, whisklwd=2,boxfill="lightblue", medcol="red", staplelwd=2)

legend(2, 330, legend ="(a)", bty="n", cex=1.5)

axis(side = 1, at = seq(1, 2, by = 1),labels=c("No", "Yes"),
     cex.axis=1.2)
axis(side = 2, at = seq(0, 300, by = 50),cex.axis=1.2, las=1)
legend(1.15, 250, legend = "( )", bty = "n", cex = 1.5, text.col="black")
legend(1.22, 242, legend = "*", bty = "n", cex = 1.5, text.col="black")


test<-sinaplot(finalB$agg~finalB$young_yn, method="density", maxwidth=0.6,scale=F, alpha=0.4, plot=F)

points(y~scaled, data=test,lwd=1)


#save 3.5x5.63

#Figure 2B

elo<-subset(finalB, treatment=="elo")

boxplot(agg~helper_yn, data=elo,xlab="Presence of helpers", 
        ylab="Number of attacks by breeders",cex.lab=1.2, las=1, cex.axis = 2,
        xlim=c(0.5,2.5),ylim=c(0,300), yaxt="n", xaxt="n", boxwex=0.2, outline=F,
        whisklty = 1, whisklwd=2,boxfill="lightblue", medcol="red", staplelwd=2)


legend(2, 330, legend ="(b)", bty="n", cex=1.5)

axis(side = 1, at = seq(1, 2, by = 1),labels=c("No", "Yes"),
     cex.axis=1.2)
axis(side = 2, at = seq(0, 300, by = 50),cex.axis=1.2, las=1)


legend(1.21, 242, legend = "**", bty = "n", cex = 1.5, text.col="black")

legend(0.25,320, legend =c("Predator of fish", expression(italic("(L. elongatus)"))),
       bty="n", cex=1.2, y.intersp = 0.8)

test<-sinaplot(elo$agg~elo$helper_yn, method="density", maxwidth=0.6,scale=F, alpha=0.4, plot=F)

points(y~scaled, data=test,lwd=1)

box()
#safe 3.5 x 5.63

##########################################################################
#Helper investment (model 3)
##########################################################################

helper <- subset(predation_final,B_H=="H")

final_treat_H <- with(helper, data.frame(
  territory = territory,
  activity_TF = scale(activity_TF,scale=F)[,1],
  treatment = treatment,
  social_rank = social_rank,
  aggression = aggression,
  helper_yn = helper_yn,
  young_yn = young_yn,
  helper_num = scale(helper_num,scale=F)[,1],
  aggr_bino = aggr_bino,
  harem = harem,
  fish_sl = scale(fish_sl,scale=F)[,1],
  fish = fish_sl,
  stimulus_id = stimulus_id,
  id=id))

str(final_treat_H)
final_treat_H <- na.omit(final_treat_H)

final_treat_H$territory<-as.factor(final_treat_H$territory)
final_treat_H$treatment<-as.factor(final_treat_H$treatment)
final_treat_H$social_rank<-as.factor(final_treat_H$social_rank)
final_treat_H$young_yn<-as.factor(final_treat_H$young_yn)
final_treat_H$helper_yn<-as.factor(final_treat_H$helper_yn)
final_treat_H$harem<-as.factor(final_treat_H$harem)
final_treat_H$aggression<-as.numeric(final_treat_H$aggression)



treat_nbi01 <- glmer(aggr_bino ~ treatment*young_yn+fish_sl+ helper_num+activity_TF+
                      (1|harem/territory), data=final_treat_H, family="binomial",
                    glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))


overdisp_fun(treat_nbi01)
summary(treat_nbi01)
drop1(treat_nbi01, test="Chisq")



final_treat_H$pred<- predict(treat_nbi01, type="response")

new<-data.frame(ddply(final_treat_H, .(treatment,young_yn), summarize,N=length(id), probability=mean(pred),
                      sd=sd(pred),se=sd/sqrt(N)))
####
#Best interaction
testInteractions(treat_nbi01, pairwise="young_yn", fixed="treatment", adjustment = "none")


##########################################################################
#Plot 3a

final_treat_H$aggr_bino<-as.character(final_treat_H$aggr_bino)
plot(aggr_bino~fish, data=final_treat_H, type="n", ylab="Attack probability of helpers", 
     xaxt = "n",yaxt="n",bty="l",cex.lab=1.2, ylim=c(0,1.08),cex.axis=2,
     xlim=c(15,56), xlab="Body size [mm]")

#line=4.5 to change distance to axis

axis(side=2, at=seq(0,1, by = 0.2), cex.axis=1.2, las=1)
axis(side=1, at=seq(15,55, by = 10), cex.axis=1.2)

legend(47.8,1.18, legend ="(a)", bty="n", cex=1.5)


final_treat_H$aggr_bino<-as.factor(final_treat_H$aggr_bino)
final_treat_H$aggr_bino<-as.character(final_treat_H$aggr_bino)
final_treat_H$aggr_bino<-as.numeric(final_treat_H$aggr_bino)
str(final_treat_H)

#plotting without interaction and treatment as random effect
d <- glmer(aggr_bino ~ young_yn+fish+ helper_num+activity_TF+
                     (1|territory)+(1|treatment), data=final_treat_H, family="binomial",
                    glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))



newdat <- expand.grid(
  fish=seq(16,54,length.out=100),activity_TF=0,aggr_bino=0,
  helper_num=0,young_yn=0, treatment=0)



mm <- model.matrix(terms(d),newdat) 
newdat$aggr_bino <- mm %*% fixef(d)
pvar1 <- diag(mm %*% tcrossprod(vcov(d),mm))
tvar1 <- pvar1+VarCorr(d)$territory[1] +VarCorr(d)$treatment[1]



newdat <- data.frame(
  newdat
  , plo = 1/(1+exp(-newdat$aggr_bino-1.96*sqrt(pvar1)))
  , phi = 1/(1+exp(-newdat$aggr_bino+1.96*sqrt(pvar1)))
  , tlo = 1/(1+exp(-newdat$aggr_bino-1.96*sqrt(tvar1)))
  , thi = 1/(1+exp(-newdat$aggr_bino+1.96*sqrt(tvar1)))
)

newdat$aggr_bino <- 1/(1+exp(-newdat$aggr_bino))

xes<-newdat$fish

coltrans<-"#00FF2288"

polygon(x=c(xes,xes[length(xes):1],xes[1]),
        y=c(newdat$plo,newdat$phi[length(xes):1],newdat$plo[1]),
        border=NA,
        col="lightgrey")

lines(newdat$fish,newdat$plo,lty=2)
lines(newdat$fish,newdat$phi,lty=2)

lines(newdat$fish,newdat$aggr_bino, lwd=2)


box()

#export 5.27 x 5.63
########################
#Figure 3b
final_treat_H$pred<-predict(treat_nbi01, type="response")
pfeH<-subset(final_treat_H, treatment=="pfe")

pfeH$young_yn<-as.character(pfeH$young_yn)


boxplot(pred~young_yn, data=pfeH, type="n", ylab="Attack probability of helpers", 
        xaxt = "n",yaxt="n",bty="l",cex.lab=1.2, ylim=c(0,1.08),cex.axis=2,
        xlim=c(0.75,2.25), xlab="Presence of young", las=1, boxwex=0.2, outline=F,
        whisklty = 1, whisklwd=2,boxfill="lightblue", medcol="red", staplelwd=2)

test<-sinaplot(pfeH$pred~pfeH$young_yn, method="density", maxwidth=0.6,scale=F, alpha=0.4, plot=F)
points(y~scaled, data=test,lwd=1)

legend(1.85, 1.185, legend ="(b)", bty="n", cex=1.5)

axis(side = 1, at = seq(1, 2, by = 1),labels=c("No", "Yes"),
     cex.axis=1.2)
axis(side = 2, at = seq(0, 1, by = 0.2),cex.axis=1.2, las=1)


legend(1.3, 0.9, legend = "**", bty = "n", cex = 1.5, text.col="black")

legend(0.65,1.15, legend =c("Predator of young", expression(italic("(G. pfefferi)"))),
       bty="n", cex=1.2, y.intersp = 0.8)

box()

#export 3.5 x 5.63
##########################################################################
#########################################################################################################################################
# Time in shelter of helpers (model 4)

predator <- subset(predation_final,B_H=="H")
time_pfe <- subset(predator, treatment=="pfe")


final_time_pfe <- with(time_pfe, data.frame(
  territory = territory,
  treatment = treatment,
  young_yn = young_yn,
  fish_sl = scale(fish_sl,scale=F)[,1],
  fish = fish_sl,
  shelter = shelter/60,
  activity_TF =scale(activity_TF,scale=F)[,1],
  id=id,
  harem = harem))

final_time_pfe<-na.omit(final_time_pfe)
str(final_time_pfe)
final_time_pfe$territory<-as.factor(final_time_pfe$territory)
final_time_pfe$young_yn<-as.factor(final_time_pfe$young_yn)

test<-lmer(shelter~young_yn+fish_sl+activity_TF+(1|harem/territory), 
           data=final_time_pfe)

qqnorm(resid(test))
qqline(resid(test))

summary(test)
drop1(test, test="Chisq")

new<-data.frame(ddply(final_time_pfe, .(young_yn), summarize,N=length(id), probability=mean(shelter),
                      sd=sd(shelter),se=sd/sqrt(N)))

######################################################################
# Figure 3c: time spend in shelter (helpers) ~ presence of young

#Figure 3c

predator <- subset(predation_final,B_H=="H")
time_pfe <- subset(predator, treatment=="pfe")


final_time_pfe <- with(time_pfe, data.frame(
  territory = territory,
  treatment = treatment,
  young_yn = young_yn,
  fish_sl = scale(fish_sl,scale=F)[,1],
  fish = fish_sl,
  shelter = shelter,
  activity_TF =scale(activity_TF,scale=F)[,1]))


final_time_pfe$young_yn<-as.character(final_time_pfe$young_yn)

boxplot(shelter~young_yn, data=final_time_pfe,xlab="Presence of young", 
        ylab="Time in breeding chamber [s]",cex.lab=1.2, las=1, cex.axis = 2,
        xlim=c(0.5,2.5),ylim=c(0,349), yaxt="n", xaxt="n", boxwex=0.2, outline=F,
        whisklty = 1, whisklwd=2,boxfill="lightblue", medcol="red", staplelwd=2)


legend(2, 380, legend ="(c)", bty="n", cex=1.5)

axis(side = 1, at = seq(1, 2, by = 1),labels=c("No", "Yes"),
     cex.axis=1.2)
axis(side = 2, at = seq(0, 300, by = 50),cex.axis=1.2, las=1)



legend(1.21, 300, legend = "*", bty = "n", cex = 1.5, text.col="black")

legend(1.07,316, legend = "( )", bty = "n", cex = 1.8)

legend(0.25,370, legend =c("Predator of young", expression(italic("(G. pfefferi)"))),
       bty="n", cex=1.2, y.intersp = 0.8)

test<-sinaplot(final_time_pfe$shelter~final_time_pfe$young_yn, method="density", maxwidth=0.6,scale=F, alpha=0.4, plot=F)

points(y~scaled, data=test,lwd=1)
box()

#Export 3.5x5.63

###########################################################################################################################################
# Male investment to the different females (model 5)

dm<-subset(predation_final, social_rank=="DM")


new<-data.frame(ddply(dm, .(territory,treatment,harem,young_yn, helper_yn), summarize,N=length(aggression), mean=mean(aggression),
                      aggression=sum(aggression), activity_TF=mean(activity_TF), SL = mean(fish_sl), Time_attack=mean(Time_attack)))

new <- subset(new, harem!="7")
new <- subset(new, harem!="8")
new <- subset(new, harem!="17")

new$territory<-as.factor(new$territory)
new$treatment<-as.factor(new$treatment)
new$harem<-as.factor(new$harem)

###############################################
Harem.data<-with(new, data.frame(
  territory = territory,
  treatment = treatment,
  young_yn = young_yn,
  SL = scale(SL,scale=F)[,1],
  activity_TF =scale(activity_TF,scale=F)[,1],
  aggression = aggression,
  harem = harem,
  helper_yn = helper_yn,
  Time_attack =Time_attack))
str(Harem.data)
Harem.data$territory<-as.factor(Harem.data$territory)
Harem.data$young_yn<-as.factor(Harem.data$young_yn)
Harem.data$harem<-as.factor(Harem.data$harem)
Harem.data<-na.omit(Harem.data)


dm_poi0 <- glmmadmb(aggression ~ young_yn+helper_yn+SL+activity_TF+(1|treatment)+(1|harem), 
                    data=Harem.data, family="nbinom", zeroInflation=F)

dm_poi <- update(dm_poi0, zeroInflation=T)
dm_poi2 <- update(dm_poi0, family="poisson",zeroInflation=F)
dm_poi2zi <- update(dm_poi0, family="poisson",zeroInflation=T)
AICtab(dm_poi, dm_poi0,dm_poi2)

summary(dm_poi)
drop1(dm_poi, test="Chisq")

