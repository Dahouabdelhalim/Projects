##Field experimental reefs LOESS to visualize for density-mortality relationships
#load FieldDDdata.txt from Field D-D folder
library(AICcmodavg)
library(ggplot2)
library(car)
library(lme4)
library(nlme)
library(multcomp)
setwd("~/PhD/Field D-D Exp")


###############################################################
############LOESS for Purple Only: 2014 and 2017 data##########
###############################################################
#load data
a<-read.table("FieldDDdataPurpleAll.txt", head=T) 
str(a)
plot(jitter(PropMort60, amount=.05) ~ InitDense, data=a, cex=1.5)

range(a$InitDense)
xinit<-seq(0,40, 0.01)

DDloess60<-loess(PropMort60 ~ InitDense, data=a, span=0.75, degree=2)#, control=loess.control(surface="direct"))
predictedLOESS60<-predict(DDloess60, data.frame(InitDense=xinit), se=T)
lines(xinit,predictedLOESS60$fit, lty=2, lwd=3, col="blue")

#To plot Fig. S2A   #Final Version
par(mar=c(5,5,.5,.5))
plot(PropMort60 ~ jitter(InitDense, amount=.5),data=a, xlab=expression(" Urchin density" ~ (m^{-2})),
     ylab="Proportional mortality", las=1, bty="L", pch="o", 
     cex=2, cex.axis=1.2, cex.lab=1.5, ylim=c(0,1), col="darkorchid4")
lines(xinit, predictedLOESS60$fit, lwd=3)
lines(xinit, predictedLOESS60$fit+predictedLOESS60$se.fit, lwd=2, lty=2, col="grey")
lines(xinit, predictedLOESS60$fit-predictedLOESS60$se.fit, lwd=2, lty=2, col="grey")
text(32, 0.925, labels="A", cex=3)

plot(DDloess60$residuals ~ InitDense, data=a) #check residuals for homogenous variance

##############################################
###Same thing for 24 h mortality (Purple only: 2014 & 2017)
##############################################
plot(jitter(PropMort24h, amount=.05) ~ InitDense, data=a, cex=1.5)

##LOESS
DDloess24<-loess(PropMort24h ~ InitDense, data=a, span=0.75, degree=2)#, control=loess.control(surface="direct"))
predictedLOESS24<-predict(DDloess24, data.frame(InitDense=xinit), se=T)
lines(xinit,predictedLOESS24$fit, lty=2, lwd=3, col="blue")
summary(DDloess24)

plot(DDloess24$residuals ~ InitDense, data=a) #check residuals for homogenous variance

#To plot Figure S2B    #Version final
par(mar=c(5,5,.5,.5))
plot(PropMort24h ~ jitter(InitDense, amount=0.5),data=a, xlab=expression(" Urchin density" ~ (m^{-2})),
     ylab="Proportional mortality", las=1, bty="L", cex=2, pch="o", 
     cex.axis=1.2, cex.lab=1.5, col="darkorchid4")
lines(xinit, predictedLOESS24$fit, lwd=3)
lines(xinit, predictedLOESS24$fit+predictedLOESS24$se.fit, lwd=2, lty=2, col="grey")
lines(xinit, predictedLOESS24$fit-predictedLOESS24$se.fit, lwd=2, lty=2, col="grey")
text(32, 0.925, labels="B", cex=3)

#############################################################

###############################################################
############LOESS for Purple & Red together: 2014 and 2017 data
###############################################################
#load data
d<-read.table("FieldDDdataBoth.txt", head=T)
d<-as.data.frame(d)

################################################
### 1 h mortality analysis Both urchin species
################################################
plot(jitter(PropMort60, amount=.05) ~ InitDense, data=d, cex=1.5)

range(d$InitDense)
xinit<-seq(0,40, 0.01)

DDloess60<-loess(PropMort60 ~ InitDense, data=d, span=0.75, degree=2)#, control=loess.control(surface="direct"))
predictedLOESS60<-predict(DDloess60, data.frame(InitDense=xinit), se=T)
lines(xinit,predictedLOESS60$fit, lty=2, lwd=3, col="blue")

#To plot Fig. S2C  #Version final
par(mar=c(5,5,.5,.5))
plot(PropMort60 ~ jitter(InitDense, amount=.5),data=d, xlab=expression(" Urchin density" ~ (m^{-2})),
     ylab="Proportional mortality", las=1, bty="L", pch="o",
     cex=2, cex.axis=1.2, cex.lab=1.5, ylim=c(0,1), col="firebrick")
lines(xinit, predictedLOESS60$fit, lwd=3)
lines(xinit, predictedLOESS60$fit+predictedLOESS60$se.fit, lwd=2, lty=2, col="grey")
lines(xinit, predictedLOESS60$fit-predictedLOESS60$se.fit, lwd=2, lty=2, col="grey")
text(32,.925, labels="C", cex=3)

##################################
###Same thing for 24 h mortality
#####################################
plot(jitter(PropMort24, amount=.05) ~ InitDense, data=d, cex=1.5)

##LOESS
DDloess<-loess(PropMort24 ~ InitDense, data=d, span=0.75, degree=2)#, control=loess.control(surface="direct"))
predictedLOESS<-predict(DDloess, data.frame(InitDense=xinit), se=T)
lines(xinit,predictedLOESS$fit, lty=2, lwd=3, col="blue")
summary(DDloess)

#To plot Figure S2D,    Final version
par(mar=c(5,5,.5,.5))
plot(PropMort24 ~ jitter(InitDense, amount=.5) , data=d, xlab=expression(" Urchin density" ~ (m^{-2})),
     ylab="Proportional mortality", las=1, bty="L", pch="o",
     cex=2, cex.axis=1.2, cex.lab=1.5, col="firebrick")
lines(xinit, predictedLOESS$fit, lwd=3)
lines(xinit, predictedLOESS$fit+predictedLOESS$se.fit, lwd=2, lty=2, col="grey")
lines(xinit, predictedLOESS$fit-predictedLOESS$se.fit, lwd=2, lty=2, col="grey")
text(32,.925, labels="D", cex=3)