# Here, we test if we can get the observed pattern solely due to regression to a mean 

# we know observed data mean is about 0.65 for the fast ones and 0.45 for the slow ones 

#install.packages("tidyverse")
library(tidyverse)
library(ggplot2)
library(plyr) 
library(plotrix)
library(nlme)
library(lme4)
library(reshape2)
library(MuMIn)
library(AICcmodavg)
library(devtools)
library(nls.multstart) 
library(reshape2)
library(cowplot)
library(TeamPhytoplankton)
library(nlsLoop)

rm(list=ls())




#############################################################
####### we have two populations being exactly the same in terms of their normal distribution of growth rates - apart from the distribution having a certain shape, this is PURE RANDOM ##### THIS IS THE ONE SINEAD HAS BEEN PLAYING WITH
###########################################################

#our first population X
X.1<-runif(100000, max = 1.0, min =0.4)
hist(X.1)
X.1<-data.frame(X.1)
X.1$measure<-ifelse(X.1$X.1<0.45, "below","above") # we make a crude split between below and above average values 

#can sample from that distribution, but can also work with the full distribution. not sure it makes a huge differences. random samples taken from a normal distribution often enough will still yield a normal distribution, so probably no point
samps<-sample(X.1$X.1, 100)
hist(samps) 

head(X.1)

#our second population Y

Y.1<- runif(100000, max = 2.1, min = 1.9)
Y.1<-data.frame(Y.1)
Y.1$measure2<-ifelse(Y.1$Y.1<0.45, "below","above") # we make a crude split between below and above average values 

XY.1<-cbind(X.1,Y.1)
head(XY.1)
XY.1$xy<-XY.1$X.1/XY.1$Y.1
XY.1$xymeas<-ifelse(XY.1$xy<1,"higher","lower")

plot1<-qplot(X.1,Y.1, data=XY.1,alpha=0.5, colour=xymeas, xlab="Growth rate in mono-culture", ylab="Growth rate in co culture") +geom_abline(slope = 1, intercept = 0, lty = 2)+theme_classic()+scale_colour_manual(values=c("darkgreen","blue"))+theme(legend.position = "top")
plot1 # is entirely random and occupies the entire space

XY.1$foldchange<-XY.1$Y.1/XY.1$X.1

#lets make a fold change plot from this 

plot2<-qplot(X.1, foldchange, data=subset(XY.1, foldchange<10), colour=xymeas, alpha=0.5, ylab="Fold change growth compared to mono culture",xlab="Growth rate in mono-culture")+theme_classic()+scale_colour_manual("Growth rate lower or higher than 1:1",values=c("darkgreen","blue"))+theme(legend.position = "top")
plot2 #results in vaguely L-shaped distribution


##################################
###### random with slightly higher mean growth rates in co-culture - uses mean data as per the experiments - otherwise STILL RANDOM#####
###################################


X.1<-rnorm(100000,mean=0.45, sd=(1/sqrt(100))) # the larger we make sd , the larger the spread around the mean 
hist(X.1)
X.1<-data.frame(X.1)
X.1$measure<-ifelse(X.1$X.1<0.45, "below","above") # we make a crude split between below and above average values 

Y.2<- rnorm(100000,mean=0.65, sd=(1/sqrt(100))) # we assume that mean growth rate and spread are the same 

Y.2<-data.frame(Y.2)
Y.2$measure2<-ifelse(Y.2$Y.2<0.65, "higher","lower") # we make a crude split between below and above average values 

XY.2<-cbind(X.1,Y.2)
head(XY.2)
XY.2$xy<-XY.2$X.1/XY.2$Y.2
XY.2$xymeas<-ifelse(XY.2$xy<1,"higher","lower")

plot3<-qplot(X.1,Y.2, data=XY.2,alpha=0.5, colour=xymeas, xlab="Growth in mono-culture", ylab="Growth in co-culture") +geom_abline(slope = 1, intercept = 0, lty = 2)+theme_classic()+scale_colour_manual(values=c("darkgreen","blue"))+theme(legend.position = "top")
plot3 #more 'natural' distribution - looks much more like our data

XY.2$foldchange<-XY.2$Y.2/XY.2$X.1
plot4<-qplot(X.1, foldchange, data=subset(XY.2, foldchange<15), colour=xymeas,ylab="Fold change growth compared to mono-culture",xlab="growth alone", alpha=0.5)+theme_classic()+scale_colour_manual("growth rate lower or higher than 1:1",values=c("darkgreen","blue"))+theme(legend.position = "top")
plot4 #

plot4.1<-qplot(X.1, foldchange, data=subset(XY.2, foldchange<15&foldchange>.1), colour=xymeas,ylab="Fold change growth compared to mono-culture",xlab="growth alone", alpha=0.5)+theme_classic()+scale_colour_manual("growth rate lower or higher than 1:1",values=c("darkgreen","blue"))+theme(legend.position = "top")
plot4.1

#can also check whether growth rate in mono-culture was above or below average
qplot(X.1, XY.2$Y.2/XY.2$X.1, data=XY.2, colour=measure, alpha=0.5)+theme_classic()+scale_colour_manual("growth rate below or above average in monoculture",values=c("darkred","orange")) # this just looks like a rotting banana

library(cowplot)
plot_grid(plot1,plot2,plot3,plot4.1, ncol=2)



##########################################################################
###################### now do the same, many times  #############################
##########################################################################

mean.growth.mono<- 0.45 # we don't need a more complicated way, e.g. mean(seq(0.25, 0.65,0.001))  as we create a normal distribution around the mean anyways.. 

increased.growth<- matrix(data=rnorm(100*100000,mean=0.65, sd=1/sqrt(100)), nrow = 10^5, ncol = 100)  # we can increase the mean by on average 0.5, and the samples in co-culture get to grab from that distribution randomly, regardless how fast or how slow they are growing. 

#now we make 100 different new normal distributions with the 100 possible increased growth rates ... with normal distriubtions that can take 100000 different values. This should be enough, no? Also, this WAS a loop. It took half an hour to run. It is now vectorized, take two lines of code ,and is much faster

Monodata = matrix(data=rnorm(100*100000,mean=mean.growth.mono, sd=1/sqrt(100)), nrow = 10^5, ncol = 100) #the ncols are our 100 different distributions
Codata <- Monodata+increased.growth # done permutating! The rest is for making things pretty

### make things pretty and plottable 
Mono.frame<-data.frame(Monodata)
Co.frame<-data.frame(Codata)
#we want all X in Co.frame to be Y
colnames(Co.frame) = gsub(pattern = "X", replacement = "Y", x = colnames(Co.frame))

library(reshape2)
Mono.melt<-melt(Mono.frame) # will give warning, but is alright, we WANT This to be using all measures as variables 
head(Mono.melt)
Co.melt<-melt(Co.frame) # will give warning, but is alright, we WANT This to be using all measures as variables 
head(Co.melt)
colnames(Co.melt)<-c("varCo","valueCo")

Monoco<-cbind(Mono.melt,Co.melt)
#plot all in the same (oh dear)

Monoco$xy<-Monoco$value/Monoco$valueCo
Monoco$xymeas<-ifelse(Monoco$xy<1,"below","above")

Monoco$runnumber<-as.numeric(substr(Monoco$variable,2,3))

# this takes exceedingly long, so we just plot the first 25 out of 100 now 

Monoco<-ddply(Monoco, .(runnumber), mutate, id = seq_along(value))

avMon <- ddply(Monoco, .(id), summarise, avVal =  mean(value, na.rm=T), avCO=mean(valueCo, na.rm=T))

avMon$xy<-avMon$avCO/avMon$avVal
avMon$xymeans<-ifelse(avMon$xy>1,"above","below")

#make plot from averaged dataframe 

bigplot3<-qplot(avVal,avCO, data=avMon, alpha=0.5, colour=xymeans) +geom_abline(slope = 1, intercept = 0, lty = 2)+theme_classic()+scale_colour_manual(values=c("darkgreen","blue"))
bigplot3

#make fold change plot ...also from averages, takes too long otherwise 

bigplot4<-qplot(avVal,avCO/avVal, data=avMon,alpha=0.5, colour=xymeans) +geom_abline(slope = 0, intercept = 1, lty = 2)+theme_classic()+scale_colour_manual(values=c("darkgreen","blue"))
bigplot4


#now need one with the predicted values only,and the MEAN predicted value overlayed on top....hello exponential decay,my old friend 
#note we did this with a non linear mixed effects model for our real data. But now we are fitting 100s of curves. We don't want to do that in an nlme frame, so we use a multi starter instead

exp_dec <- function(a, b, growth) {
  y <- a * exp(growth * b) # a and be being the intercept and steepness of slope 
  return(y)
}

avMon$foldchange<-avMon$avCO/avMon$avVal
Monoco$foldchange<-Monoco$valueCo/Monoco$value

#we trial this for the first 100 curves. increase to 10 000 later 

Monoco1<-subset(Monoco,id<100 & foldchange<10)
id<-unique(Monoco1$id)

exp.fits <- nlsLoop::nlsLoop(Monoco1,
                              model = foldchange ~ exp_dec(a, b, growth=value),
                              id='id', # right now we are just fitting the first 100 out of 10 000 runs  - will need to re run when not on laptop 
                              tries = 100, #should increase this to 1000 later - but am on slow computer just now. 
                              param_bds = c(0.01, 10, -10, -0.0001),
                              na.action = na.omit,
                             lower = c(a=0.01, b=-10),
                             upper = c(a=10, b=-0.0001))

exp.fits.params<-exp.fits$params
exp.fits.preds<-exp.fits$predictions

### now plot with predicted slope 
exp.fits.params$runnumber <- substr(exp.fits.params$id, 1, 2)
exp.fits.preds$runnumber<- substr(exp.fits.preds$id, 1, 2)


predplot_01<-ggplot(Monoco1) +
  geom_point(aes(x=value, y=foldchange,alpha =0.5, colour=xymeas)) + scale_colour_manual(values=c('darkgreen','blue'))+labs(x=expression(Growth~alone), y=expression(Fold~change~growth)) + geom_point(data=exp.fits.preds,aes(x=value, y=foldchange), colour='red', size=0.4) +theme_classic()+theme(legend.position='top')

predplot_01



##########################################################################
###################### run many times with upper hard limit #############################
##########################################################################

#####now there is an upper hard limit to how much they can increase
#so we use the same dataframe for growth alone, but in co-cultured growth, we don't alow the increase to go beyond a total of 1.0 

increased.growth2<- matrix(data=rnorm(100*100000,mean=0.65, sd=1/sqrt(1000)), nrow = 10^5, ncol = 100) # making sd smaller limits the distribution of samples to pick from in the first place  - this creates the upper limit from which we can sample, while allowing the mean to be shifted. We could of course also shift the mean in increments, but I thought it may make sense to use the new mean as we find it in the ThinCert experiment. 

Codata2 <-increased.growth2 # or can try Monodata+ increased.growth2 when we just increase by a mean increment.. 


Co.frame2<-data.frame(Codata2)
#we want all X in Co.frame to be Y
colnames(Co.frame2) = gsub(pattern = "X", replacement = "Y", x = colnames(Co.frame2))

Co.melt2<-melt(Co.frame2) # will give warning, but is alright, we WANT This to be using all measures as variables 
head(Co.melt2)
colnames(Co.melt2)<-c("varCo","valueCo")

Monoco2<-cbind(Mono.melt,Co.melt2)

Monoco2$xy<-Monoco2$value/Monoco2$valueCo
Monoco2$xymeas<-ifelse(Monoco2$xy<1,"above","below")

Monoco2$runnumber<-as.numeric(substr(Monoco2$variable,2,3))
Monoco2$valueCo2<-ifelse(Monoco2$valueCo>1,mean(Monoco2$valueCo),Monoco2$valueCo)

#make averaged dataframe. May need seq_along within each runnumber, and then make means by that 
Monoco2<-ddply(Monoco2, .(runnumber), mutate, id = seq_along(value))

avMon2 <- ddply(Monoco2, .(id), summarise, avVal =  mean(value, na.rm=T), avCO=mean(valueCo, na.rm=T))

avMon2$xy<-avMon2$avCO/avMon2$avVal
avMon2$xymeans<-ifelse(avMon2$xy>1,"above","below")

#make plot from averaged dataframe 

bigplot3.2<-qplot(avVal,avCO, data=subset(avMon2,avCO/avVal<10), alpha=0.5, colour=xymeans) +geom_abline(slope = 1, intercept = 0, lty = 2)+theme_classic()+scale_colour_manual(values=c("darkgreen","blue"))
bigplot3.2

#make fold change plot ...also from averages, takes too long otherwise 

bigplot4.2<-qplot(avVal,avCO/avVal, data=subset(avMon2,avCO/avVal<10),alpha=0.5, colour=xymeans) +geom_abline(slope = 0, intercept = 1, lty = 2)+theme_classic()+scale_colour_manual(values=c("darkgreen","blue"))
bigplot4.2

#now need one with the predicted values only,and the MEAN predicted value overlayed on top....hello exponential decay,my old friend 
#note we did this with a non linear mixed effects model for our real data. But now we are fitting 100s of curves. We don't want to do that in an nlme frame, so we use a multi starter instead

exp_dec <- function(a, b, growth) {
  y <- a * exp(growth * b) # a and be being the intercept and steepness of slope 
  return(y)
}

avMon2$foldchange<-avMon2$avCO/avMon2$avVal
Monoco2$foldchange<-Monoco2$valueCo/Monoco2$value

#we trial this for the first 100 curves. increase to 10 000 later 

Monoco1.2<-subset(Monoco2,id<=100 & foldchange<10)
id<-unique(Monoco1.2$id)

exp.fits2 <- nlsLoop::nlsLoop(Monoco1.2,
                             model = foldchange ~ exp_dec(a, b, growth=value),
                             id='id', # right now we are just fitting the first 100 out of 10 000 runs  - will need to re run when not on laptop 
                             tries = 100, #should increase this to 1000 later - but am on slow computer just now. 
                             param_bds = c(0.01, 10, -10, -0.0001),
                             na.action = na.omit,
                             lower = c(a=0.01, b=-10),
                             upper = c(a=10, b=-0.0001))

exp.fits.params2<-exp.fits2$params
exp.fits.preds2<-exp.fits2$predictions

### now plot with predicted slope 
exp.fits.params2$runnumber <- substr(exp.fits.params2$id, 1, 2)
exp.fits.preds2$runnumber<- substr(exp.fits.preds2$id, 1, 2)


predplot_02<-ggplot(Monoco1.2) +
  geom_point(aes(x=value, y=foldchange,alpha =0.5, colour=xymeas)) + scale_colour_manual(values=c('darkgreen','blue'))+labs(x=expression(Growth~alone), y=expression(Fold~change~growth)) + geom_point(data=exp.fits.preds2,aes(x=value, y=foldchange), colour='red', size=0.4) +theme_classic()+theme(legend.position='top')+geom_abline(slope = 0, intercept = 1, lty = 2)

predplot_02

