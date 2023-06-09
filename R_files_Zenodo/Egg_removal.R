rm(list=ls())

library(plotrix)
library(car)
library(multcomp)
library(MASS)
library(nlme)
library(lme4)
library(lmerTest)




# Survival #
############

library(coxme)
data = read.csv2("survival_egg_rem.csv", header=T)
ColID= as.factor(data$ColID)
Treat= as.factor(data$Treat)
Week = as.numeric(data$Week.during.treatment)
Death = as.numeric(data$Dead)

coxmod1 <- coxme(Surv(time = Week,event= Death) ~  Treat + (1|ColID), data=data)
anova(coxmod1)

# fertility #
#############

library(plotrix)
library(car)
library(multcomp)
library(MASS)
library(nlme)
library(lme4)
library(lmerTest)

data = read.csv2("Ovary egg rem.csv", header=T)
ColID= as.factor(data$ColID)
Frag= as.factor(data$Frag)
Treat= as.factor(data$Treat)
Ov = as.numeric(data$Ovmean)
Egg = as.numeric(data$white_egg)

# ovary length
lmm1 = lmer(Ov~Treat+(1|ColID/Frag), data=data)
Anova(lmm1)

# number of white eggs

lmm1 = glmer(Egg~Treat+(1|ColID/Frag), family = "poisson", data=data)
Anova(lmm1)




# ovary length graph

moyenne <- aggregate(Ov~Treat,FUN=mean,data=data);names(moyenne)[names(moyenne)=="Egg"]<-"mean"
stderror <- aggregate(Ov~Treat,FUN=std.error,data=data);names(stderror)[names(stderror)=="Egg"]<-"std.error"
summary_data <-merge(moyenne,stderror,all.x=T,all.y=T)
summary_data
mea.Wtot = c( 3.00,3.21)
se.Wtot = c(0.11,0.12)
patric= c("Control","Treatment")
Wegg.tr = barplot(mea.Wtot, ylab="Ovary length (mm)" ,names = patric, cex.lab=0.8, cex.names=0.8, col=c("white","darkgrey"), space=0, ylim=c(0,5))
arrows(Wegg.tr,mea.Wtot-se.Wtot,Wegg.tr, mea.Wtot+se.Wtot, lwd=1.3, angle=90,length=0.1,code=3)

# White egg graph

moyenne <- aggregate(Egg~Treat,FUN=mean,data=data);names(moyenne)[names(moyenne)=="Egg"]<-"mean"
stderror <- aggregate(Egg~Treat,FUN=std.error,data=data);names(stderror)[names(stderror)=="Egg"]<-"std.error"
summary_data <-merge(moyenne,stderror,all.x=T,all.y=T)
summary_data
mea.Wtot = c(7.385965,10.837209)
se.Wtot = c(0.6047817,0.8587574)
patric= c("Control","Treatment")
Wegg.tr = barplot(mea.Wtot, ylab="Egg number" ,names = patric, cex.lab=0.8, cex.names=0.8, col=c("white","darkgrey"), space=0, ylim=c(0,15))
arrows(Wegg.tr,mea.Wtot-se.Wtot,Wegg.tr, mea.Wtot+se.Wtot, lwd=1.3, angle=90,length=0.1,code=3)
text(x=1, y=14, labels="**",font=1, ps=0.1)