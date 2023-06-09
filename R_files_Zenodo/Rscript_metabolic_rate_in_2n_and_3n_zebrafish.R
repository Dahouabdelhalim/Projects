rm(list = ls())
# load packages
require(nlme)
require(car)
library(lme4)
library(AICcmodavg)
library(rsquared)
require(piecewiseSEM)
require(MuMIn)

# read data files
surv <- read.csv("Survival_rate.csv",header=T)
d<- read.csv("Metabolic_data.csv",header=T)
e<-subset(d,d$Length!="x")
e$Length<-as.numeric(e$Length)
data<-read.csv("O2 measurements corrected for blanks.csv")
t<-data[which(as.character(round(data$O2))=="20"),]

# Preliminary analyses

# survival
s1<-lme(dead.~ploidy+rtemp,random = ~ 1 | batch,data=surv)
Anova(s1, type="III")
summary(s1) # no effect of rearing temperature, but higher mortality in 3n

# length
l1<-lmer(Length~ploidy+relevel(as.factor(rearT),ref="26.5") +(1 | batch),data=e)
anova(l1) # length differs with rearing temperature, but not between 2n and 3n
r.squaredGLMM(l1) # 9.7% of the variation explained

# analyses and models

mVmax<-lme(Vmax_2~testT*rearT*efficiency,random = ~ 1 | batch,data=d) # model on 422 obs without body size
summary(mVmax) # Table S1
Anova(mVmax, type="III") # Table 1

mVmaxb<-lme(Vmax_2~I(Length^3)+testT*rearT*efficiency,random = ~ 1 | batch,data=e) # model on the 410 obs that do have body size
Anova(mVmaxb, type="III")
summary(mVmaxb) # Table S2
summary(backupl)

mKm<-lme(Km_2~testT+rearT+efficiency,random = ~ 1 | batch,data=d) # best fitted model for P50 value
summary(mKm) # Table S3
Anova(mKm, type="III") # Table 3

# Ea calculations for Table 2
E1<-lme(log(Vmax_2)~I(1/(testT+273.15))*rearT*efficiency,random = ~ 1 | batch,data=d)
summary(E1)
-0.000086173324* (5927.769 +23.5*-469.801) #2n reared at 23.5
-0.000086173324* (5927.769 +26.5*-469.801) #2n reared at 26.5
-0.000086173324* (5927.769 +29.5*-469.801) #2n reared at 29.5
-0.000086173324* (5927.769 +23.5*(-469.801+100*6.843)+100*-162.723)  #3n reared at 23.5
-0.000086173324* (5927.769 +26.5*(-469.801+100*6.843)+100*-162.723)  #3n reared at 26.5
-0.000086173324* (5927.769 +29.5*(-469.801+100*6.843)+100*-162.723)  #3n reared at 29.5
# Q10 calculations
param <- expand.grid(testT=c(23.5,29.5),rearT=c(23.5,26.5,29.5),efficiency=c(0,100),batch=c(1:8))
param$pred<-predict(mVmax,newdata=param)
# example for Q10 calculation for 23 reared at 23.5
low<-mean(param$pred[which(param$testT==23.5&param$rearT==23.5&param$efficiency==100)])
high<-mean(param$pred[which(param$testT==29.5&param$rearT==23.5&param$efficiency==100)])
(high/low)^(10/6) #1.8008

# visualisations

# Fig 1
par(mfrow=c(1,3))
plot(breakppoint~Km_1,data=d,pch=19,col="#00000020",cex=2,cex.axis=2,cex.lab=2,xlab="",ylab="Pcrit")  
abline(a=0,b=1,col="red")
plot(breakppoint~Km_2,data=d,pch=19,col="#00000020",cex=2,cex.axis=2,cex.lab=2,xlab="P50",ylab="")  # H = 2 correlates best with the breakpoint
abline(a=0,b=1,col="red")
plot(breakppoint~Km_3,data=d,pch=19,col="#00000020",cex=2,cex.axis=2,xlab="",ylab="")  
abline(a=0,b=1,col="red")

plot(Vmax_1~meanMO2,data=d,pch=19,col="#00000020",cex=2,cex.axis=2,cex.lab=2,xlab="",ylab="MO2routine")
abline(a=0,b=1,col="red")
plot(Vmax_2~meanMO2,data=d,pch=19,col="#00000020",cex=2,cex.axis=2,cex.lab=2,xlab="P50",ylab="") # H = 2 correlates best with mean MO2
abline(a=0,b=1,col="red")
plot(Vmax_3~meanMO2,data=d,pch=19,col="#00000020",cex=2,cex.axis=2,xlab="",ylab="") 
abline(a=0,b=1,col="red")

#Fig 2
m<-lmer(Vmax_2~testT*rearT*efficiency+ (1 | batch),data=d,REML = FALSE)
plotT=c(23.5,26.5,29.5)
par(mfrow=c(1,3))
par(mar=c(5,4,2,1))
for(i in 1:3){
  plot(Vmax_2~jitter(testT,0.6),main=plotT[i],data=d[which(d$rearT==plotT[i]),],ylim=c(5,20),pch=c(19,17)[as.factor(d$ploidy[which(d$rearT==plotT[i])])],col=c("#FF000070","#FFA50070")[as.factor(d$ploidy[which(d$rearT==plotT[i])])],xlab="",ylab="",xaxt="n",cex=2,cex.axis=2)
  axis(1, at=c(23.5,26.5,29.5),cex.axis=2)
  if(i==3){legend(23,20, legend = c("2n","3n"), col=c("#FF0000","#FFA500"), pch = c(19,17), cex = 2, bty = "n")}
  if(i==2){title(xlab="Test Temperature", cex.lab=2)}
  if(i==1){title(ylab="Oxygen Consumption", cex.lab=2)}
  pred<-predict (m,newdata=data.frame(testT=c(23:30),rearT=plotT[i],efficiency=100),re.form=NA)
  lines(c(23:30),pred,lwd=3,col="#FFA500")
  pred<-predict (m,newdata=data.frame(testT=c(23:30),rearT=plotT[i],efficiency=0),re.form=NA)
  lines(c(23:30),pred,lwd=3,col="#FF0000")
}

# Figure 3
par(mfrow=c(1,1))
plot(Km_2~jitter((testT-0.25),0.5),col="#FF000035",pch=19,data=d[which(d$ploidy=="2n"),],xaxt='n',ylim=c(0,38),xlim=c(23,30),xlab="Test temperature",ylab="P50")
points(Km_2~jitter((testT+0.25),0.5),col="#FFA50045",pch=17,data=d[which(d$ploidy=="3n"),])
axis(side=1,at=c(23.5,26.5,29.5))
legend(27,39, legend = c("2n","3n"), col=c("#FF0000","#FFA500"), pch = c(19,17), cex = 1, bty = "n")
abline(a=-18.632037+(-0.172867*26.5),b=1.339687,col="#FF0000",lwd=3)
abline(a=-18.632037+1.3064+(-0.172867*26.5),b=1.339687,col="#FFA500",lwd=3)

# Figure S1
# show individual traces of oxygen saturation corrected for blanks
par(mfrow=c(1,2))
data2n<-data[which(data$ploidy=="2n"),]
data3n<-data[which(data$ploidy=="3n"),]
plot(O2~absTime.s.,data=data2n,xlab="time(h)",ylab=expression("O"[2]*" (%)"),xlim=c(0,10),ylim=c(0,105),col=c("#0000FF15","#00CC0015","#FF000015")[as.factor(data2n$testT)],pch=92,cex=0.6,main="Diploid")
plot(O2~absTime.s.,data=data3n,xlab="time(h)",ylab=expression("O"[2]*" (%)"),xlim=c(0,10),ylim=c(0,105),col=c("#0000FF15","#00CC0015","#FF000015")[as.factor(data3n$testT)],pch=92,cex=0.6,main="Triploid")
legend(0.5,100, c(expression(paste("23.5",degree,"C")),expression(paste("26.5",degree,"C")),expression(paste("29.5",degree,"C"))), col=c("#0000FF","#00CC00","#FF0000"), pch=19, cex=1,bty = "n")

# Figure S2 
# show body size as a function of incubation time, incubation temperature
boxplot(Length~ploidy*rearT,data=e,col=c("red","orange"),cex.axis=1.5)

# Figure S3
# show mean time until 20% saturation is reached
par(mfrow=c(1,1))
boxplot(absTime.s.~testT,data=t,xlab="test temperature",ylab="time to reach 20% saturation (h)",col=c("blue","green", "red"))
# calculate the Ea for the rate (i.e. the inverse of time) at which oxygen is depleted to 20% saturation
m1<-lm(log(1/absTime.s.)~I(1/(testT+273.15)),data=t)
summary(m1)
summary(m1)$coefficients[2, 1]*-0.000086173324 #Ea = 0.5385
(-6249.0432+179.6108)*-0.000086173324 # Ea 0.5230232-0.5539785
