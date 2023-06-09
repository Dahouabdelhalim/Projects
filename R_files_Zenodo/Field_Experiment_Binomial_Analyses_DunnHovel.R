###Density-dependence analysis, purple urchins 2014 & 2017
###Combine 60 min and 24 h sampling events, analyze with generalized mixed effects model (binomial dist)
###then compare density levels using ANOVA with Tukey post hoc
###Based on reviewer 2 comments from Ecology submission (August 2018)

library(lme4)
library(nlme)
library(ggplot2)
library(plyr)
library(dplyr)
library(multcomp)
library(car)
library(bbmle)

setwd("~/PhD/Field D-D Exp")
bad<-read.table("BinomialAnalysisData.txt", header=T, sep="\\t", 
                colClasses=c("Year"="factor","Treatment"="factor", "RemainCount" = "integer",
                             "Trial" = "factor", "Hour"="factor","Reef"="factor"))
str(bad)
 
bad$InitDense<-bad$InitCount/(.75*.75)  #add density column to use as predictor,
bad$PropMort<-1-(bad$RemainCount/bad$InitCount) #add Proportional Mortality Column for reference
str(bad)
##subset into PurpleOnly & Both datasets 
badpurp<-subset(bad, Treatment=="Purple")
badboth<-subset(bad, Treatment=="Both")
str(badpurp);str(badboth)

###################################################
##Generalized linear mixed effects model###########
##Binomial Distribution############################
##Full model to explicitly test red urchin effect##

####################################################
########## with init density as numeric##########
####################################################

fullmod<-glmer(1-(RemainCount/InitCount) ~ 
                 InitDense + Treatment + Hour + InitDense*Treatment +
                 Site +  (1|Trial/Reef), data=bad,
               family="binomial", weights = InitCount)#, glmerControl(optimizer="bobyqa"))
summary(fullmod)  #Interaction N-S

noInt<-glmer(1-(RemainCount/InitCount) ~
               InitDense + Treatment +  Hour +
               (1|Site) + (1|Trial/Reef), data=bad,
             family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa"))
summary(noInt) #Treatment Significant here

noInt2<-glmer(1-(RemainCount/InitCount) ~
                InitDense + Treatment +  Hour +
                Site + (1|Trial/Reef), data=bad,
              family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa"))
summary(noInt2)  #Treatment, Hour, Site all significant

noInt3<-glmer(1-(RemainCount/InitCount) ~  
                InitDense + Treatment +  Hour +
                (1|Trial/Reef), data=bad,
              family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa"))
summary(noInt3)

AICtab(fullmod, noInt, noInt2, noInt3, weights=TRUE)#most weight (55%) behind noInt2, so separate purple... 
AICctab(fullmod, noInt, noInt2, noInt3, weights=TRUE, nobs=202)  #... and red data and run analyses on each.
summary(noInt2)                             


###########################################
#####PURPLE TRIALS ONLY####################
###########################################
plot(PropMort~jitter(InitDense), data=subset(badpurp, Hour==1), cex=2, pch="o", col="darkorchid4")
points(PropMort~jitter(InitDense), data=subset(badpurp, Hour==24), cex=2, pch="o", col="blue")

fullpurp<-glmer((InitCount-RemainCount)/InitCount ~
                  InitDense + Year +  Site + Hour + 
                  InitDense*Year + (1|Trial/Reef), data=badpurp,
                family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
summary(fullpurp) # Interaction not significant, but this model gets the most AIC (bc it includes Hour)

noIntPurp<-glmer(1-(RemainCount/InitCount) ~ 
                   InitDense + Year + Hour + (1|Site) + (1|Trial/Reef), data=badpurp,
                 family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa"))
summary(noIntPurp)
noIntPurp2<-glmer(1-(RemainCount/InitCount) ~    
                    InitDense + Year + Site + Hour + (1|Trial/Reef), data=badpurp,
                  family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa",
                                                                     optCtrl=list(maxfun=2e6)))
summary(noIntPurp2) # this doesn't converge
noIntPurp3<-glmer(1-(RemainCount/InitCount) ~   
                    InitDense + Year + Site + (1|Trial/Reef), data=badpurp,
                  family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa")) 
summary(noIntPurp3)
noIntPurp4<-glmer(1-(RemainCount/InitCount) ~   
                    InitDense + Year + Hour + (1|Trial/Reef), data=badpurp,
                  family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
summary(noIntPurp4)
noIntPurp5<-glmer(1-(RemainCount/InitCount) ~    
                    InitDense +  Hour + Site + (1|Year) + (1|Trial/Reef), data=badpurp,
                  family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa"))
summary(noIntPurp5) #doesn't converge
noIntPurp6<-glmer((InitCount-RemainCount)/InitCount ~    
                    Hour + Site + Year + (1|Trial/Reef), data=badpurp,
                  family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa"))
summary(noIntPurp6)
print(noIntPurp6)
AICtab( noIntPurp,  noIntPurp3, noIntPurp4, noIntPurp6 , weights=T)
AICctab(noIntPurp,  noIntPurp3, noIntPurp4, noIntPurp6 , weights=T)
#most weight (90%) to noIntPurp6, the model WITHOUT URCHIN DENSITY AS A PREDICTOR
summary(noIntPurp6)
anova(noIntPurp6)
summary(pairwise<-glht(noIntPurp6,linfct=mcp(Site="Tukey")))  #post-hoc comparison of mean mortalities
summary(pairwise<-glht(noIntPurp6,linfct=mcp(Year="Tukey")))  #post-hoc comparison of mean mortalities

#########################################
########Purple & Red Together############
#########################################
plot(PropMort~jitter(InitDense), data=subset(badboth, Hour==1), cex=2, pch="o", col="firebrick")
points(PropMort~jitter(InitDense), data=subset(badboth, Hour==24), cex=2, pch="o", col="pink")

fullboth<-glmer(1-(RemainCount/InitCount) ~ 
                  InitDense + Site + Hour + InitDense*Site  + 
                  (1|Trial/Reef), data=badboth,
                family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa"))
summary(fullboth) #throws a warning, be cautious with this model

noIntboth<-glmer(1-(RemainCount/InitCount) ~ 
                   InitDense + Site + Hour + (1|Trial/Reef), data=badboth,
                 family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
summary(noIntboth) #

noIntboth2<-glmer(1-(RemainCount/InitCount) ~ 
                    InitDense +  Hour + (1|Site) + (1|Trial/Reef), data=badboth,
                  family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
summary(noIntboth2) #

noIntboth3<-glmer(1-(RemainCount/InitCount) ~ 
                    Hour + (1|Site) + (1|Trial/Reef), data=badboth,
                  family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
summary(noIntboth3) #

AICtab(fullboth, noIntboth, noIntboth2, noIntboth3, weights=T)
AICctab(fullboth, noIntboth, noIntboth2, noIntboth3, weights=T)#noIntboth gets 70% support

summary(noIntboth)
Anova(noIntboth)
summary(pairwise<-glht(noIntboth,linfct=mcp(Site="Tukey")))  #post-hoc comparison of mean mortalities
summary(pairwise<-glht(noIntboth,linfct=mcp(Hour="Tukey")))  #post-hoc comparison of mean mortalities


###########################################
####Separate High and Low densities##
###############################
#Subset data into high and low, by urchin treatment (separate purple only from purple + red experiments)
badpurplo<-subset(badpurp,InitDense < 12)
badbothlo<-subset(badboth, InitDense <12)
badpurphi<-subset(badpurp,InitDense > 12)
badbothhi<-subset(badboth, InitDense > 12)

##########################
#Purple Only
##Low density
plot(PropMort~jitter(InitDense), data=subset(badpurplo, Hour==1), xlim=c(0,40), ylim=c(0,1), 
     las=1, bty="L", xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality", 
     cex.lab=1.5, cex.axis=1.2, cex=2, pch="o", col="darkorchid4")
points(PropMort~jitter(InitDense), data=subset(badpurplo, Hour==24), cex=2, pch=2, col="darkorchid4")
lopurp<-glmer(1-(RemainCount/InitCount) ~    
                    InitDense +  Site + Year + Hour + (1|Trial/Reef), data=badpurplo,
                  family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
#summary(lopurp)  
lopurp2<-glmer((InitCount-RemainCount)/InitCount ~      ###USE THIS MODEL###
                InitDense +  Year + Hour + (1|Trial/Reef), data=badpurplo,
              family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
summary(lopurp2)  #Density p=0.02
anova(lopurp2) 
lopurp3<- glmer(1-(RemainCount/InitCount) ~    
                  InitDense +  Site +  Hour + (1|Trial/Reef), data=badpurplo,
                family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
lopurp4<-glmer(1-(RemainCount/InitCount) ~    
                  Site + Year + Hour + (1|Trial/Reef), data=badpurplo,
               family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
lopurp5<-glmer(1-(RemainCount/InitCount) ~    
                  Site +  Hour + (1|Trial/Reef), data=badpurplo,
               family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
lopurp6<- glmer(1-(RemainCount/InitCount) ~    
                  InitDense + Site + Hour + InitDense*Site + (1|Trial/Reef), data=badpurplo,
                family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
AICctab(lopurp, lopurp2,lopurp3, lopurp4,lopurp5 ,lopurp6, weights=T)    #lopurp2 gets 71% weight
AICtab(lopurp, lopurp2,lopurp3, lopurp4,lopurp5 , lopurp6, weights=T)    #lopurp2 gets 67% weight
#summary(pairwise<-glht(lopurp2,linfct=mcp(InitDense="Tukey")))  #can't run post hoc when treating dens as numeric 
summary(pairwise<-glht(lopurp2,linfct=mcp(Year="Tukey")))  #post-hoc of mean 

####To plot predicted lines for lopurp2 model#####
#1 h example line
loxinit<-seq(0,12,0.01)
hixinit<-seq(12,40,0.01)

newdata<-data.frame(InitDense=loxinit, Year=factor("2014", levels = levels(badpurplo$Year)),
                    Trial=factor("10", levels = levels(badpurplo$Trial)),
                    Reef=factor("6", levels= levels(badpurplo$Reef)),
                    Hour=factor("1", levels=levels(badpurplo$Hour)))
predicted<-predict(lopurp2, newdata, type='response')

###24 h 
newdatab<-data.frame(InitDense=loxinit, Year=factor("2014", levels = levels(badpurplo$Year)),
                    Trial=factor("10", levels = levels(badpurplo$Trial)),
                    Reef=factor("6", levels= levels(badpurplo$Reef)),
                    Hour=factor("24", levels=levels(badpurplo$Hour)))
predictedb<-predict(lopurp2, newdatab, type='response')

##########################
#Purple Only Hi Density
plot(PropMort~jitter(InitDense, amount=.5), data=subset(badpurphi, Hour==1), xlim=c(0,40), ylim=c(0,1), 
     las=1, bty="L", xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality",
     cex.lab=1.5, cex.axis=1.2, cex=2, pch="o", col="darkorchid4")
points(PropMort~jitter(InitDense, amount=.5), data=subset(badpurphi, Hour==24), cex=2, pch=2, col="darkorchid4")
hipurp<-glmer(1-(RemainCount/InitCount) ~    
                InitDense +  Site + Year + Hour + (1|Trial/Reef), data=badpurphi,
              family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa"))
summary(hipurp)   
hipurp2<-glmer(1-(RemainCount/InitCount) ~    
                 InitDense +  Year + Hour + (1|Trial/Reef), data=badpurphi,
               family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
summary(hipurp2)  #Density p=0.02
hipurp3<- glmer(1-(RemainCount/InitCount) ~    
                  InitDense +  Site + Hour + (1|Trial/Reef), data=badpurphi,
                family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
hipurp4<- glmer(1-(RemainCount/InitCount) ~    
                  Year + Hour + (1|Trial/Reef), data=badpurphi,
                family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
hipurp5<- glmer(1-(RemainCount/InitCount) ~    
                  Site + Hour + (1|Trial/Reef), data=badpurphi,
                family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
hipurp6<- glmer(1-(RemainCount/InitCount) ~    
                   Hour + (1|Trial/Reef), data=badpurphi,
                family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
hipurp7<- glmer(1-(RemainCount/InitCount) ~    #Use this model
                 Site + Year + Hour + (1|Trial/Reef), data=badpurphi,
                family="binomial", weights=InitCount, glmerControl(optimizer="bobyqa"))
AICtab(hipurp,hipurp2,hipurp3,hipurp4,hipurp5,hipurp6,hipurp7, weights=TRUE)
AICctab(hipurp,hipurp2,hipurp3,hipurp4,hipurp5,hipurp6,hipurp7, weights=TRUE)  #use hipurp7, 70+% of weight
anova(hipurp7)
summary(hipurp7)
#summary(pairwise<-glht(hipurp,linfct=mcp(InitDense="Tukey")))  #can't run post hoc when treating dens as numeric 
summary(pairwise<-glht(hipurp7,linfct=mcp(Year="Tukey")))  #post-hoc of mean 

####To plot predicted lines for hipurp7 model#####
#1 h
plot(PropMort~jitter(InitDense, amount=.5), data=subset(badpurphi, Hour==1), xlim=c(0,40), ylim=c(0,1), 
     las=1, bty="L", xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality",
     cex.lab=1.5, cex.axis=1.2, cex=2, pch="o", col="darkorchid4")
points(PropMort~jitter(InitDense, amount=.5), data=subset(badpurphi, Hour==24), cex=2, pch=2, col="darkorchid4")

newdatahi<-data.frame(InitDense=hixinit, Year=factor("2014", levels = levels(badpurphi$Year)),
                    Trial=factor("2", levels = levels(badpurphi$Trial)),
                    Site = factor("SLJ", levels=levels(badpurphi$Site)),
                    Reef=factor("6", levels= levels(badpurphi$Reef)),
                    Hour=factor("1", levels=levels(badpurphi$Hour)))
predictedhi<-predict(hipurp7, newdatahi, type='response')
lines(hixinit, predictedhi, lwd=3)
points(PropMort ~ jitter(InitDense, amount=.5), data=subset(badpurphi, Hour==1), 
       cex=2, pch="o", col="darkorchid4")

#24 h
newdatahib<-data.frame(InitDense=hixinit, Year=factor("2014", levels = levels(badpurphi$Year)),
                     Trial=factor("2", levels = levels(badpurphi$Trial)),
                     Site = factor("SLJ", levels=levels(badpurphi$Site)),
                     Reef=factor("6", levels= levels(badpurphi$Reef)),
                     Hour=factor("24", levels=levels(badpurphi$Hour)))
predictedhib<-predict(hipurp7, newdatahib, type='response')


###Plot for Figure 2A###   Final Version
tiff("Figure2A.tif", width=3, height=2, units="in", res=800, pointsize=5, compression="lzw")
par(mar=c(5,5,2,2))
plot(PropMort~jitter(InitDense, amount=.5), data=subset(badpurplo, Hour==1), xlim=c(0,40), ylim=c(0,1), 
     las=1, bty="L", xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality", 
     cex.lab=1.9, cex.axis=1.7, cex=2, pch="o", col="darkorchid4")
lines(loxinit, predicted, lwd=2)
points(PropMort~jitter(InitDense, amount=.5), data=subset(badpurphi, Hour==1),cex=2, pch="o", col="darkorchid4")
lines(hixinit,predictedhi, lwd=2)
text(38.5,.925, labels="A", cex=2.5)
dev.off()

###Plot for Figure 2B###  Final version
tiff("Figure2B.tif", width=3, height=2, units="in", res=800, pointsize=5, compression="lzw")
par(mar=c(5,5,2,2))
plot(PropMort~jitter(InitDense, amount=.5), data=subset(badpurplo, Hour==24), xlim=c(0,40), ylim=c(0,1), 
     las=1, bty="L", xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality", 
     cex.lab=1.9, cex.axis=1.7, cex=2, pch="o", col="darkorchid4")
lines(loxinit, predictedb, lwd=2, col="black")
points(PropMort ~ jitter(InitDense, amount=.5), data=subset(badpurphi, Hour==24), 
       cex=2, pch="o", col="darkorchid4")
lines(hixinit, predictedhib, lwd=2, col="black")
text(38.5,.925, labels="B", cex=2.5)
dev.off()

###############################
###############################
#Purple & Red together  Low Density
plot(PropMort~jitter(InitDense, amount=.5), data=subset(badbothlo, Hour==1), xlim=c(0,40), ylim=c(0,1), 
     las=1, bty="L", xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality",
     cex.lab=1.5, cex.axis=1.2, cex=2, pch="o", col="firebrick")
points(PropMort ~ jitter(InitDense, amount=.5), data=subset(badbothhi, Hour==1), 
       cex=2, pch="o", col="firebrick")

plot(PropMort~jitter(InitDense, amount=.5), data=subset(badbothlo, Hour==24), xlim=c(0,40), ylim=c(0,1), 
     las=1, bty="L", xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality",
     cex.lab=1.5, cex.axis=1.2,  cex=2, pch="o", col="firebrick")
points(PropMort ~ jitter(InitDense, amount=.5), data=subset(badbothhi, Hour==24), 
       cex=2, pch="o", col="firebrick")

loboth<-glmer(1-(RemainCount/InitCount) ~    
                InitDense +  Site +  Hour + (1|Trial/Reef), data=badbothlo,
              family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
summary(loboth)  #Density N-S
anova(loboth)
summary(pairwise<-glht(loboth,linfct=mcp(Site="Tukey")))  #post-hoc of mean 

loboth2<-glmer(1-(RemainCount/InitCount) ~    
                 InitDense +   Hour + (1|Trial/Reef), data=badbothlo,
               family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
loboth3<-glmer(1-(RemainCount/InitCount) ~    
                 Site +  Hour + (1|Trial/Reef), data=badbothlo,
               family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
loboth4<-glmer(1-(RemainCount/InitCount) ~    
                 Hour + (1|Trial/Reef), data=badbothlo,
               family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))

AICctab(loboth, loboth2, loboth3, loboth4, weights=T)   #loboth3 gets 79% of weight
summary(loboth3)
#1 h prediction curve
newdata2<-data.frame(InitDense=loxinit, Trial=factor("B2", levels = levels(badbothlo$Trial)),
                     Site=factor("SLJ", levels=levels(badbothlo$Site)),
                     Reef=factor("4", levels= levels(badbothlo$Reef)),
                     Hour=factor("1", levels=levels(badbothlo$Hour)))
predicted2<-predict(loboth3, newdata2, type='response')

#24 h prediction curve
newdata5<-data.frame(InitDense=loxinit, Trial=factor("B2", levels = levels(badbothlo$Trial)),
                     Site=factor("SLJ", levels=levels(badbothlo$Site)),
                     Reef=factor("4", levels= levels(badbothlo$Reef)),
                     Hour=factor("24", levels=levels(badbothlo$Hour)))
predicted5<-predict(loboth3, newdata5, type='response')


###############################
#Purple & Red together  High Density
plot(PropMort~jitter(InitDense, amount=.5), data=subset(badbothhi, Hour==1), xlim=c(0,40), ylim=c(0,1), 
     las=1, bty="L", xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality",
     cex.lab=1.5, cex.axis=1.2, cex=2, pch="o", col="firebrick")
points(PropMort~jitter(InitDense, amount=.5), data=subset(badbothhi, Hour==24), cex=2, pch=2, col="firebrick")

hiboth<-glmer(1-(RemainCount/InitCount) ~    
                InitDense +  Site +  Hour + (1|Trial/Reef), data=badbothhi,
              family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))

hiboth2<-glmer(1-(RemainCount/InitCount) ~    
                InitDense +  Hour + (1|Trial/Reef), data=badbothhi,
              family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
hiboth3<-glmer(1-(RemainCount/InitCount) ~    
                Site +  Hour + (1|Trial/Reef), data=badbothhi,
              family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
hiboth4<-glmer(1-(RemainCount/InitCount) ~    
                Hour + (1|Trial/Reef), data=badbothhi,
              family="binomial", weights=InitCount)#, glmerControl(optimizer="bobyqa"))
AICtab(hiboth, hiboth2, hiboth3, hiboth4, weights=T) #hiboth3 gets 65%
AICctab(hiboth, hiboth2, hiboth3, hiboth4, weights=T) #hiboth3 gets ~70% of weight
summary(hiboth3)

#1 h prediction curve
newdatahi2<-data.frame(InitDense=hixinit, Trial=factor("B2", levels = levels(badbothhi$Trial)),
                     Site=factor("SLJ", levels=levels(badbothhi$Site)),
                     Reef=factor("3", levels= levels(badbothhi$Reef)),
                     Hour=factor("1", levels=levels(badbothhi$Hour)))
predictedhi2<-predict(hiboth3, newdatahi2, type='response')

#24 h prediction curve
newdatahi5<-data.frame(InitDense=hixinit, Trial=factor("B2", levels = levels(badbothhi$Trial)),
                     Site=factor("SLJ", levels=levels(badbothhi$Site)),
                     Reef=factor("3", levels= levels(badbothhi$Reef)),
                     Hour=factor("24", levels=levels(badbothhi$Hour)))
predictedhi5<-predict(hiboth3, newdatahi5, type='response')


###Plot for Figure 2C###   Version final
tiff("Figure2C.tif", width=3, height=2, units="in", res=800, pointsize=5, compression="lzw")
par(mar=c(5,5,2,2))
plot(PropMort~jitter(InitDense, amount=.5), data=subset(badbothlo, Hour==1), xlim=c(0,40), ylim=c(0,1), 
     las=1, bty="L", xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality", 
     cex.lab=1.9, cex.axis=1.7, cex=2, pch="o", col="firebrick")
lines(loxinit, predicted2, lwd=2)
points(PropMort~jitter(InitDense, amount=.5), data=subset(badbothhi, Hour==1),cex=2, pch="o", col="firebrick")
lines(hixinit, predictedhi2, lwd=2)
text(38.5,.925, labels="C", cex=2.5)
dev.off()

###Plot for Figure 2D###   Final Version
tiff("Figure2D.tif", width=3, height=2, units="in", res=800, pointsize=5, compression="lzw")
par(mar=c(5,5,2,2))
plot(PropMort~jitter(InitDense, amount=.5), data=subset(badbothlo, Hour==24), xlim=c(0,40), ylim=c(0,1), 
     las=1, bty="L", xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality", 
     cex.lab=1.9, cex.axis=1.7, cex=2, pch="o", col="firebrick")
lines(loxinit, predicted5, lwd=2, col="black")
points(PropMort ~ jitter(InitDense, amount=.5), data=subset(badbothhi, Hour==24), 
       cex=2, pch="o", col="firebrick")
lines(hixinit, predictedhi5, lwd=2, col="black")
text(38.5,.925, labels="D", cex=2.5)
dev.off()
