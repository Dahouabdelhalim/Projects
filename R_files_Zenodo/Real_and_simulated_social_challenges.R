
#libraries
library(ggplot2)
library(magrittr)
library(tidyverse)
library(lubridate)
library(glmmTMB)
library(DHARMa)
library(Rmisc)
library(lme4)
library(lmerTest)
library(Hmisc)


#input Experiment 1 data

#set wd

#set file name
fn<- "Exp_1_data.csv"
#import data
STIdata<- read.csv(fn, na.strings=c("","NA"))


#test that time of day does not affect T

cor.test(STIdata$capture_time_of_day_h, STIdata$logT, method="spearman")

# test that latency from capture to bleed does not affect T

cor.test(STIdata$bleed_latency_h, STIdata$logT, method="spearman")

#test of decoy and playback IDs on T in the experimental group (supplmental table 1)

STI<- STIdata %>% filter(Treatment=="STI")

m<-lm(logT~STI.decoy, data=STI)
anova(m)

m<-lm(logT~STI.audio, data=STI)
anova(m)

m<-lm(STI.prop.min.agg~STI.decoy, data=STI)
anova(m)

m<-lm(STI.prop.min.agg~STI.audio, data=STI)
anova(m)

## Statistical test for results in Table 1, Model 1


m<-lm(logT~Treatment+Mass_g+Age, data=STIdata)
anova(m)


## Statistical test for results in Table 1, Model 2

m<-lm(logT~STI.prop.min.agg+Captime.since.30min.end +Mass_g+Age, data=STI)
anova(m)

## code for making Figure 2a


sum<-summarySE(STIdata, measurevar = "logT", groupvars = "Treatment", na.rm=T)
STI <- STIdata[ which(STIdata$Treatment=='STI'), ]
CTL <- STIdata[ which(STIdata$Treatment=='Control'), ]

#create T plot
const=0.25
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(2.5,1.2,0), oma=c(0,0,0,0)) #margins
par(bty="l")
set.seed(5)
ep=0.5
cp= 1.5
stripchart(logT~Treatment,data=STIdata, tck=-0.02, bty="l",vertical=T, las=1, xaxt="n",
           cex=0, cex.lab=1.2, bg="black", xlim=c(0,2), ylim=c(-3.5,0.5),
          axes=T,  ylab="log Testosterone (ng/mL plasma)", xlab="Treatment",
           group.names=c("Control","STI")) 
axis(side=1, tck=-0.02, cex.lab=2.5,  at= c(ep, cp), labels=c("Control \\n","Experimental \\n (STI)"))


pcol=adjustcolor("darkgray",alpha.f=0.5)
mcol= adjustcolor("darkgray", alpha.f=0.5)
rect(ep-const,sum[1,"logT"]-sum[1,"se"],ep+const,sum[1,"logT"]+sum[1,"se"],col=pcol,border=NA)
rect(cp-const,sum[2,"logT"]-sum[2,"se"],cp+const,sum[2,"logT"]+sum[2,"se"],col=mcol,border=NA)

lines(c(ep-const,ep+const),rep(sum[1,"logT"],2),col="black",lwd=2)
lines(c(cp-const,cp+const),rep(sum[2,"logT"],2),col="black",lwd=2)

stripchart(CTL$logT,  tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.1, pch=1,  cex=1, 
           add=T, at=ep,
           axes=T, bty="l")
stripchart(STI$logT, tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.1, pch=19, cex=1, 
           add=T, at= cp,
           axes=T, bty="l")


### making Figure 2b

ggplot(data=STI, aes(x=STI.prop.min.agg, y=logT)) + 
  geom_smooth(method='lm',formula=y~x, color="black", se=TRUE)+ 
  xlab("Proportion of minutes with aggression")+ ylab("Log testosterone (ng/mL plasma)")+ 
  theme_classic()+ theme(text=element_text(size=13))+  geom_point(size=2)

############## Exp 2: effects of treatment on rates of aggression

#load Y1 observational data

fn<-"Exp2_obs_data.csv"
Exp2_obs<- read.csv(fn, na.strings=c("","NA"))

Exp2_Y1_obs <- Exp2_obs %>% filter(Year=="Y1")

Exp2_Y1_obs$exp_day<-as.factor(Exp2_Y1_obs$exp_day)

#### model to obtain results in Table 2

m<-glmmTMB(physagg_permin_perpair~Treatment*exp_day+start_time_h+(1|Site.block),
           ziformula=~1,family=ziGamma (link="log"), data=Exp2_Y1_obs)
summary(m)

#check model residuals
sim_res<-simulateResiduals(m)
plot(sim_res)

### make figure 3a

EXP <- Exp2_Y1_obs[ which(Exp2_Y1_obs$Treatment=="EXP"), ]
CTL <- Exp2_Y1_obs[ which(Exp2_Y1_obs$Treatment=="CTL"), ]


const=0.25
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(3, 1, 0),oma=c(0,0,0,0)) #margins
par(bty="l")
set.seed(11)
cp=0.5
ep= 1.5
stripchart(physagg_permin_perpair~Treatment, data=Exp2_Y1_obs, tck=-0.02, bty="l",vertical=T, las=1, xaxt="n",
           cex=0, cex.lab=1.2, bg="black", xlim=c(0,2), ylim=c(0,0.26),
           ylab="rate of physical aggression",
           xlab="Y1 Treatment", axes=T,  
           group.names=c("Control","Experimental")) 
axis(side=1, tck=-0.02, cex.lab=2.5,  at= c(cp, ep), labels=c("",""))

boxplot(physagg_permin_perpair~Treatment, axes = FALSE, outline=FALSE, add=TRUE, at= c(cp, ep),col="white", data=Exp2_Y1_obs)

stripchart(CTL$physagg_permin_perpair,  tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.3, pch=1,  cex=1, 
           add=T, at=cp,
           axes=T, bty="l", col=alpha("black", 1))
stripchart(EXP$physagg_permin_perpair, tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.3, pch=19, cex=1, 
           add=T, at= ep,
           axes=T, bty="l", col=alpha("black", 1))

#load Y2 obs data


Exp2_Y2_obs<- Exp2_obs %>% filter(Year=="Y2")

#model for Table 2 results, Y2
m2<-glmmTMB(physagg_permin_perpair~Treatment+start_time_h+ (1|Site.block), ziformula=~1,family=ziGamma (link="log"),data=Exp2_Y2_obs)
summary(m2)

#check residuals
sim_resm2<-simulateResiduals(m2)
plot(sim_resm2)

# make figure 3b


EXP <- Exp2_Y2_obs[ which(Exp2_Y2_obs$Treatment=="EXP"), ]
CTL <- Exp2_Y2_obs[ which(Exp2_Y2_obs$Treatment=="CTL"), ]


const=0.25
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(3, 1, 0),oma=c(0,0,0,0)) #margins
par(bty="l")
set.seed(11)
cp=0.5
ep= 1.5
stripchart(physagg_permin_perpair~Treatment, data=Exp2_Y2_obs, tck=-0.02, bty="l",vertical=T, las=1, xaxt="n",
           cex=0, cex.lab=1.2, bg="black", xlim=c(0,2), ylim=c(0,0.26),
           ylab="rate of physical aggression",
           xlab="Y2 Treatment", axes=T,  
           group.names=c("Control","Experimental")) 
axis(side=1, tck=-0.02, cex.lab=2.5,  at= c(cp, ep), labels=c("",""))

boxplot(physagg_permin_perpair~Treatment, axes = FALSE, outline=FALSE, add=TRUE, at= c(cp, ep),col="white", data=Exp2_Y2_obs)

stripchart(CTL$physagg_permin_perpair,  tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.3, pch=1,  cex=1, 
           add=T, at=cp,
           axes=T, bty="l", col=alpha("black", 1))
stripchart(EXP$physagg_permin_perpair, tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.3, pch=19, cex=1, 
           add=T, at= ep,
           axes=T, bty="l", col=alpha("black", 1))

######################

#load RFID data

fn<-"Exp2_Y2_RFIDsums.csv"

RFID<-read.csv(fn, na.strings=c("","NA"))

### stats results for # intrusion frequency by treatment in Y1, Experiment 2

m<-glmmTMB((F_intrud_ev_N)~Treatment +(1|Exp.group), data=RFID, ziformula=~1, family=poisson(link="log"))
summary(m) 


##########################

# Experiment 2 testosterone 

##Y1 

#load data

fn<-"Exp2_T.csv"
T_data<-read.csv(fn, na.strings=c("","NA"))

Y1_T_data <- T_data %>% filter (Year=="Y1")

Y1_T_data$exp_day<-as.factor(Y1_T_data$exp_day)

## stats results for Table 3, Y1

m1<-lmer(logT~Treatment*exp_day+Age+Mass_g+ (1|Site.block), data=Y1_T_data)
summary(m1)
anova(m1)


# testing for effects of day and nestbox category (supplemental table 2), Y1
EXP <- Y1_T_data[ which(Y1_T_data$Treatment=='EXP'), ]
CTL <- Y1_T_data[ which(Y1_T_data$Treatment=='CTL'), ]


m<-lmer(logT~Status+exp_day+ (1|Site.block), data=EXP)
summary(m)
anova(m)

m<-lmer(logT~Status+exp_day+ (1|Site.block), data=CTL)
summary(m)
anova(m)

#make Figure 4a

sum<- summarySE(data=Y1_T_data, measurevar="logT", groupvars="Treatment")

EXP <- Y1_T_data[ which(Y1_T_data$Treatment=='EXP'), ]
CTL <- Y1_T_data[ which(Y1_T_data$Treatment=='CTL'), ]


new.exp<-EXP[ which(EXP$Status=='new'), ]
neighbor.exp<-EXP[ which(EXP$Status=='neighbor'), ]
floater.exp<-EXP[ which(EXP$Status=='floater'), ]
neighbor.ctl<-CTL[which(CTL$Status=='neighbor'), ]
emptied.ctl<-CTL[which(CTL$Status=='emptied'), ]


const=0.25
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(3, 1, 0),oma=c(0,0,0,0)) #margins
par(bty="l")
set.seed(6)
cp=0.5
ep= 1.5
stripchart(logT~Treatment,data=Y1_T_data, tck=-0.02, bty="l",vertical=T, las=1, xaxt="n",
           cex=0, cex.lab=1.2, bg="black", xlim=c(0,2), ylim=c(-3.5,0),
           ylab="log Testosterone (ng/mL plasma)",
           xlab="Y1 Treatment", axes=T,  
           group.names=c("Control","Experimental")) 
axis(side=1, tck=-0.02, cex.lab=2.5,  at= c(cp, ep), labels=c("",""))
pcol=adjustcolor("darkgray",alpha.f=0.5)
mcol= adjustcolor("darkgray", alpha.f=0.5)
rect(cp-const,sum[1,"logT"]-sum[1,"se"],cp+const,sum[1,"logT"]+sum[1,"se"],col=pcol,border=NA)
rect(ep-const,sum[2,"logT"]-sum[2,"se"],ep+const,sum[2,"logT"]+sum[2,"se"],col=mcol,border=NA)

lines(c(ep-const,ep+const),rep(sum[2,"logT"],2),col="black",lwd=2)
lines(c(cp-const,cp+const),rep(sum[1,"logT"],2),col="black",lwd=2)

stripchart(neighbor.ctl$logT,  tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.18, pch=1,  cex=1, 
           add=T, at=cp,
           axes=T, bty="l")

stripchart(emptied.ctl$logT,  tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.18, pch=2,  cex=1, 
           add=T, at=cp,
           axes=T, bty="l")


stripchart(new.exp$logT, tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.18, pch=17, cex=1, 
           add=T, at= ep,
           axes=T, bty="l")

stripchart(neighbor.exp$logT, tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.18, pch=19, cex=1, 
           add=T, at= ep,
           axes=T, bty="l")

stripchart(floater.exp$logT, tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.18, pch=15, cex=1, 
           add=T, at= ep,
           axes=T, bty="l")

legend("bottomleft", c("neighbor", "emptied", "neighbor", "new", "floater"),
       pch=c(1, 2, 19, 17, 15), bty="n", cex=1)


#### Y2

#load data

Y2_T_data<-T_data %>% filter (Year=="Y2")

# stats results for Table 3, Y2

m1<-lmer(logT~Treatment+Age+Mass_g+(1|Site.block), data=Y2_T_data)
summary(m1)
anova(m1)


# testing for effects of nestbox category (supplemental table 2), Y2
EXP<-subset(Y2_T_data, Treatment=="EXP")

m<-lmer(logT~Status+(1|Site.block), data=EXP)
summary(m)
anova(m)

#make Figure 4b

EXP <- Y2_T_data[ which(Y2_T_data$Treatment=='EXP'), ]
CTL <- Y2_T_data[ which(Y2_T_data$Treatment=='CTL'), ]

new.exp<-EXP[ which(EXP$Status=='new'), ]
neighbor.exp<-EXP[ which(EXP$Status=='neighbor'), ]
floater.exp<-EXP[ which(EXP$Status=='floater'), ]
neighbor.ctl<-CTL[which(CTL$Status=='neighbor'), ]

sum<- summarySE(Y2_T_data, measurevar="logT", groupvars="Treatment")

#create T plot
const=0.25
par(mar=c(4.5,4.5,0.5,0.5), mgp=c(3, 1, 0),oma=c(0,0,0,0)) #margins
par(bty="l")
set.seed(5)
cp=0.5
ep= 1.5
windowsFonts(A = windowsFont("Calibri"))
stripchart(logT~Treatment,data=Y2_T_data, tck=-0.02, bty="l",vertical=T, las=1, xaxt="n",
           cex=0, cex.lab=1.2, bg="black", xlim=c(0,2), ylim=c(-3.5,0),
           ylab="log Testosterone (ng/mL plasma)", 
           xlab="Y2 Treatment", axes=T,  
           group.names=c("Control","Experimental")) 
axis(side=1, tck=-0.02, cex.lab=2.5,  at= c(cp, ep), labels=c("",""))
#title(main="Experiment 2", line=-1)
pcol=adjustcolor("darkgray",alpha.f=0.5)
mcol= adjustcolor("darkgray", alpha.f=0.5)
rect(cp-const,sum[1,"logT"]-sum[1,"se"],cp+const,sum[1,"logT"]+sum[1,"se"],col=pcol,border=NA)
rect(ep-const,sum[2,"logT"]-sum[2,"se"],ep+const,sum[2,"logT"]+sum[2,"se"],col=mcol,border=NA)

lines(c(ep-const,ep+const),rep(sum[2,"logT"],2),col="black",lwd=2)
lines(c(cp-const,cp+const),rep(sum[1,"logT"],2),col="black",lwd=2)

stripchart(neighbor.ctl$logT,  tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.1, pch=1,  cex=1, 
           add=T, at=cp,
           axes=T, bty="l")

stripchart(new.exp$logT, tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.1, pch=17, cex=1, 
           add=T, at= ep,
           axes=T, bty="l")

stripchart(neighbor.exp$logT, tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.1, pch=19, cex=1, 
           add=T, at= ep,
           axes=T, bty="l")

stripchart(floater.exp$logT, tck=-0.01,las=1,vertical=T,
           method="jitter",jitter=0.1, pch=15, cex=1, 
           add=T, at= ep,
           axes=T, bty="l")

legend("topright",
       c("Control neighbor","Exp neighbor", "Exp new", "Exp floater"),
       pch=c(1, 19, 17, 15), cex=0.7)




