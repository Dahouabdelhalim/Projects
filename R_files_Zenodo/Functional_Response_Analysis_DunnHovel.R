##Lobster functional response fitting with package "frair"
library(frair)
library(readr)
library(bbmle)
library(emdbook)

setwd("~/PhD/Functional Response")
###########Functional Response Analyses############
###################################################

#load Both red & purple datasets: (LobsFuncRespDataBoth2.txt & LobsFuncRespData.txt)
Both<-read.table("LobsFuncRespDataBoth2.txt")
colnames(Both)<-c("Killed", "PropMort", "Initial")
Both$Density<-Both$Initial/(((22/7)*(2.18/2)^2)/2) #add data column with initial density

#Set up purple only dataset
Purp<-read.table("LobsFuncRespData.txt")
colnames(Purp)<-c("Killed", "PropMort", "Initial")
Purp$Density<-Purp$Initial/(((22/7)*(2.18/2)^2)/2)

##fitting a type II curve to gammarid data (an example of how to operate the FRAIR package)
data(gammarus)
frair_responses()
outII <- frair_fit(eaten~density, data=gammarus, response=
'rogersII', start=list(a = 1.2, h = 0.015), fixed=list(T=40/24))
plot(outII)
outIIb<-frair_boot(outII)
confint(outIIb)
plot(outIIb, type='n')
lines(outIIb, all_lines=T)
points(outIIb, pch=20)
plot(outIIb, type='n', main="Empirical 95% CI")
drawpoly(outIIb, col=rgb(0,0.5,0))
points(outIIb, pch=20)
lines(outIIb, all_lines=F)

##########################
###My Analysis############
##########################

##Purple urchins only####
purptest<-frair_test(Killed~Initial, data=Purp) #phenomonological test of a type II or type III curve
purptest
##Units for handling time = day
outI<-frair_fit(Killed~Initial, data=Purp, response='rogersII', start=list(a=.2, h=0.07), fixed=list(T=48/24)) #fit Rogers
outIb<-frair_boot(outI) #bootstrapping
confint(outIb)  #extract confidence intervals for parameter values

#Plot for Fig 1A ###FINAL VERSION###
tiff("Figure1A.tif", width=3, height=2, units="in", res=800, pointsize=5, compression="lzw")
par(mar=c(5,5,2,2))
plot(outIb, type="n", ylim=c(0,6), bty="L", xlab="Urchins offered", ylab="Urchins eaten", 
     cex.lab=2, cex.axis=1.75, las=1)
lines(outIb, all_lines=T)
lines(outIb, all_lines=F, lwd=3, col="grey")
points(Killed~jitter(Initial), data=Purp, pch="o", cex=2, col="darkorchid4")
text(22,5.5,labels="A", cex=2.5)  #if exporting via point-click, use dimensions = 550 x 400
dev.off()

#Didn't use this plot
outIb
plot(outIb, type="n", ylim=c(0,6))
points(jitter(Killed)~Initial, data=Purp, cex=1.5)
drawpoly(outIb, col="purple")
lines(outIb, all_lines=F)

###########################
#Both purple & red together
###########################
bothtest<-frair_test(Killed~Initial, data=Both) #phenomonological test of a type II or type III curve
bothtest

outII<-frair_fit(Killed~Initial, data=Both, response='rogersII', start=list(a=.2, h=0.07), fixed=list(T=48/24))#fit rogers
outIIb<-frair_boot(outII)#bootstrapping
confint(outIIb) #extract confidence intervals for parameter values

#Plot for Fig 1B ###FINAL VERSION###
tiff("Figure1B.tif", width=3, height=2, units="in", res=800, pointsize=5, compression="lzw")
par(mar=c(5,5,2,2))
plot(outIIb, type="n",ylim=c(0,6), bty="L", xlab="Urchins offered", ylab="Urchins eaten", 
     cex.lab=2, cex.axis=1.75, las=1)
lines(outIIb, all_lines=T)
lines(outIIb, all_lines=F, lwd=3, col="grey")
points(Killed~jitter(Initial), data=Both, pch="o", cex=2, col="firebrick")
text(22,5.5, labels="B", cex=2.5) #if exporting via point-click, use dimensions = 550 x 400
dev.off()
#didn't use this plot
outIIb
plot(outIIb, type="n", ylim=c(0,6))
points(jitter(Killed)~Initial, data=Both, cex=1.5)
drawpoly(outIIb, col="red")
lines(outIIb, all_lines=F)

#################################################################
####Alternate test for type II vs III explicitly testing with AIC
#Purple only
p_flex <- frair_fit(Killed ~ Initial, data=Purp,
                    response='flexpnr',
                    start=list(b = 1, q = 0, h = 0.07),
                    fixed=list(T = 48/24))
# Fit a model where q is fixed to zero:
p_II <- frair_fit(Killed ~ Initial, data=Purp,
                  response='flexpnr',
                  start=list(b = 1, h = 0.07),
                  fixed=list(T = 48/24, q = 0))
summary(p_flex$fit) # q = 0 : Type II preferred
AIC(p_flex$fit, p_II$fit)
AICtab(p_flex$fit, p_II$fit, weights=T)

#Purple & Red 
b_flex <- frair_fit(Killed ~ Initial, data=Both,
                    response='flexpnr',
                    start=list(b = 1, q = 0, h = 0.07),
                    fixed=list(T = 48/24))
# Fit a model where q is fixed to zero:
b_II <- frair_fit(Killed ~ Initial, data=Both,
                  response='flexpnr',
                  start=list(b = 1, h = 0.07),
                  fixed=list(T = 48/24, q = 0))
summary(b_flex$fit) # q = 0 : Type II preferred
AICtab(b_flex$fit, b_II$fit, weights=T) # The model without q is preferred


####################################
#Proportional Mortality & LOESS#####
####################################

##Purple Only
propmortlmpurp<-lm(PropMort ~ Density, data=Purp, contrasts=T)#fit lm of purp. prop. mort
summary(propmortlmpurp) #Statistical output for lm of purple urchin proportional mortality
xinitdense=seq(0,15,.01) #set up sequence to predict across
predictedlmpurp<-predict(propmortlmpurp, data.frame(Density=xinitdense), se=T) #predict from lm

#Plot for Fig. 1C##  #FINAL VERSION
tiff("Figure1C.tif", width=3, height=2, units="in", res=800, pointsize=5, compression="lzw")
par(mar=c(5,5,2,2))
plot(PropMort ~ jitter(Density), data=Purp, bty="L", las=1, pch="o", 
     cex.lab=1.9, cex.axis=1.7, cex=2, ylim=c(0,1.0), 
     xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality", col="darkorchid4")
lines(xinitdense, predictedlmpurp$fit, lwd=2)
text(11.75, 0.95, labels="C", cex=2.5)
dev.off()


###Both urchin species
propmortlmboth<-lm(PropMort~Density, data=Both, contrasts=T)
summary(propmortlmboth)
predictedlmboth<-predict(propmortlmboth,data.frame(Density=xinitdense), type="response", se=T)

#Plot for Figure 1D ##FINAL VERSION##
tiff("Figure1D.tif", width=3, height=2, units="in", res=800, pointsize=5, compression="lzw")
par(mar=c(5,5,2,2))
plot(PropMort ~ jitter(Density), data=Both, bty="L", las=1, pch="o", 
     cex.lab=1.9, cex.axis=1.7, cex=2, ylim=c(0,1.0), 
     xlab=expression(" Urchin density" ~ (m^{-2})), ylab="Proportional mortality",col="firebrick")
lines(xinitdense,predictedlmboth$fit, col="black", lwd=2)
text(11.75, 0.95, labels="D", cex=2.5)
dev.off()
