######################################################################################################################################
## Script by Arianne F. Messerman
## 13 April, 2021
## This script examines differences in body size metrics between juvenile spotted and marbled salamanders (Ambystoma maculatum and A. opacum)
##   immediately prior to initial release into semi-natural outdoor experimental enclosures.
## Please see the associated manuscript for full study details: 
##   Messerman, A.F. and M. Leal. Submitted. The contributions of individual traits to survival among terrestrial juvenile 
##   pond-breeding salamanders. Ecology.
######################################################################################################################################

#setwd()

dat<-read.csv("Metamorph_2017_Release_Data.csv", header=TRUE)
str(dat)
AMMA<-subset(dat, dat$Species=="AMMA")
AMOP<-subset(dat, dat$Species=="AMOP")

mod.mass<-lm(log(Mass)~Species, dat)
summary(mod.mass) #p=0.261
anova(mod.mass)
mean(AMMA$Mass)#1.317
sd(AMMA$Mass)#0.393
mean(AMOP$Mass)#1.226
sd(AMOP$Mass)#0.253
mean(dat$Mass)#1.271
sd(dat$Mass)#0.333
mean(AMMA$Mass)/mean(AMOP$Mass)*100 #7.4995% AMMA greater than AMOP

mod.svl<-lm(log(SVL)~Species, dat)
summary(mod.svl)#p=0.895
anova(mod.svl)
mean(AMMA$SVL)#29.619
sd(AMMA$SVL)#2.986
mean(AMOP$SVL)#29.488
sd(AMOP$SVL)#2.091
mean(dat$SVL)#29.554
sd(dat$SVL)#2.571

mod.tl<-lm(log(TL)~Species, dat)
summary(mod.tl)#p<0.001
anova(mod.tl)
mean(AMMA$TL)#58.012
sd(AMMA$TL)#6.856
mean(AMOP$TL)#54.286
sd(AMOP$TL)#4.425
mean(dat$TL)#56.149
sd(dat$TL)#6.048
mean(AMMA$TL)/mean(AMOP$TL)*100 #6.864% AMMA greater than AMOP
