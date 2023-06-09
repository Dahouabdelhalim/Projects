
setwd("C:/Users/samue/OneDrive/PhD/PhD OneDrive/Chapter Two Familiarity or Intensity/Data V2")

library(lme4)
library(glmmML)
library(pscl)
library(AER)
library(MASS)
library(RVAideMemoire)
library(sjPlot)

opt.expon.mod <- function(model, expon=seq(4/100, 4, length=100)) {
  # Function that establishes to which exponent (power) a variable 
  # should be raised to maximise the probability of the supplied
  # model's residuals testing as being normally distributed using 
  # the Shapiro-Wilks test.
  #
  # The function takes two arguments: the base model and the
  # range of exponents to be tested. This range is set by default
  # to 100 values between 0 (but not including 0) and 4.
  #
  # Author: Emile van Lieshout
  # Date: 1/8/2013
  
  
  # Could use the commented-out code below to test whether specific
  # exponents < 1 equivalent to even roots are inappropriate if the
  # response contains negative numbers.
  #	form <- formula(model)
  #	dep <- as.character(form[2])
  #	rest <- as.character(c(form[c(1,3)]))
  
  start <- proc.time()
  pvals <- numeric(length(expon))
  opt_expon <- 0
  resid.best <- NULL
  
  for(i in 1:length(expon)) {
    form.upd <- as.formula(paste(".^", expon[i], " ~ .", sep=''))
    model.upd <- update(model, form.upd)
    
    pvals[i] <- shapiro.test(residuals(model.upd))$p.value
    
    if(pvals[i] >= max(pvals)) {
      opt_expon <- expon[i]
      resid.best <- residuals(model.upd)
    }
  }
  
  cat('Maximum observed Shapiro-Wilks p =', max(pvals), 'at exponent', opt_expon, '\\n')
  
  layout(matrix(1:2, ncol=1))
  plot(expon, pvals, 
       type='l', 
       xlab="Dependent variable exponent used", 
       ylab="Observed Shapiro-Wilk test p")
  hist(resid.best)
  layout(1)
  
  cat("Calculated in", (proc.time() - start)[3], "seconds\\n")
  
}

###EJACULATE WEIGHT###

ed=read.csv("EjacData3F.csv")
ed$Treatment=as.factor(ed$Treatment)
ed$FW=scale(ed$F.Weight)
ed$MW=scale(ed$W1)
names(ed)

#Test patterns in male weight
mw1=lmer(W1~Total+Sim*Fam+(1|Eb.ID),data=ed)
plot(mw1)
mw2=lm(W1~Total+Sim*Fam,data=ed)
plot(mw2)
anova(mw1)
anova(mw2)
mw3=lmer(W1~Total+Sim+Fam+(1|Eb.ID),data=ed)
anova(mw1,mw3)

###############


EW1=lmer(Ejac~(Sim+Total)*Fam*MW+(1|Eb.ID),data=ed)
EW1b=lm(Ejac~(Sim+Total)*Fam*MW,data=ed)
plot(EW1)
plot(EW1b)
shapiro.test(EW1b$residuals)
opt.expon.mod(EW1)
opt.expon.mod(EW1b)
#Can't correct
#Remove 131, 149, 182

ed2=read.csv("EjacData3Fb.csv")
ed2$Treatment=as.factor(ed2$Treatment)
ed2$FW=scale(ed2$F.Weight)
ed2$MW=scale(ed2$W1)
names(ed2)

EW1=lmer(Ejac~(Sim+Total)*Fam*MW+(1|Eb.ID),data=ed2)
EW1b=lm(Ejac~(Sim+Total)*Fam*MW,data=ed2)
plot(EW1)
plot(EW1b)
shapiro.test(EW1b$residuals)
opt.expon.mod(EW1)
opt.expon.mod(EW1b)
#Can't correct
#Remove 31, 131, 152

ed3=read.csv("EjacData3Fc.csv")
ed3$Treatment=as.factor(ed3$Treatment)
ed3$FW=scale(ed3$F.Weight)
ed3$MW=scale(ed3$W1)
names(ed3)

T1=subset(ed3,Treatment=="1")
length(T1$Ejac)
T2=subset(ed3,Treatment=="2")
length(T2$Ejac)
T3=subset(ed3,Treatment=="3")
length(T3$Ejac)
T4=subset(ed3,Treatment=="4")
length(T4$Ejac)
T5=subset(ed3,Treatment=="5")
length(T5$Ejac)

EW1=lmer(Ejac~(Sim+Total)*Fam*MW+(1|Eb.ID),data=ed3)
EW1b=lm(Ejac~(Sim+Total)*Fam*MW,data=ed3)
plot(EW1)
plot(EW1b)
shapiro.test(EW1b$residuals)
opt.expon.mod(EW1)
opt.expon.mod(EW1b)
#Now it is normal (p=0.071)

###Test whether MW differs among treatments###
mw1=lm(W1~Sim+Total+Fam,data=ed3)
mw2=lm(W1~Total+Fam,data=ed3)
anova(mw1,mw2)
#No effect of Sim
mw3=lm(W1~Sim+Fam,data=ed3)
anova(mw1,mw3)
#Marginally non-sig effect of total
mw4=lm(W1~Sim+Total,data=ed3)
anova(mw1,mw4)
#Marginally non sig effect of Fam
#Include MW as a covariate

#Comparison of T1 and T4, to determine whether absolute exposure time is a problem.
T1T4=rbind(T1,T4)
time1=lmer(Ejac~Treatment*MW+(1|Eb.ID),data=T1T4)
time2=lmer(Ejac~Treatment+MW+(1|Eb.ID),data=T1T4)
anova(time1,time2)
#Remove MW
time3=lmer(Ejac~MW+(1|Eb.ID),data=T1T4)
anova(time2,time3)
#No sig effect of treatment

#Calculate number of groups
levels(ed3$Eb.ID)

###Now test for ejaculate weight across all combinations
ew1=lmer(Ejac~Sim+Total+Fam+MW
         +Total*Fam
         +Sim*MW
         +Total*MW
         +Fam*MW
         +(1|Eb.ID),
           data=ed3)
ew2=lm(Ejac~Sim+Total+Fam+MW
         +Total*Fam
         +Sim*MW
         +Total*MW
         +Fam*MW,
         data=ed3)
anova(ew1,ew2)
#Keep Eb.ID
ew3=lmer(Ejac~Sim+Total+Fam+MW
         +Total*Fam
         +Sim*MW
         +Total*MW
         +Fam+MW
         +(1|Eb.ID),
         data=ed3)
anova(ew1,ew3)
#Fam by MW non sig
ew4=lmer(Ejac~Sim+Total+Fam+MW
         +Total*Fam
         +Sim*MW
         +Total+MW
         +Fam*MW
         +(1|Eb.ID),
         data=ed3)
anova(ew1,ew4)
#Total by MW non sig
ew5=lmer(Ejac~Sim+Total+Fam+MW
     +Total*Fam
     +Sim+MW
     +Total*MW
     +Fam*MW
     +(1|Eb.ID),
     data=ed3)
anova(ew1,ew5)
#Sim by MW non sig
ew6=lmer(Ejac~Sim+Total+Fam+MW
     +Total*Fam
     +(1|Eb.ID),
     data=ed3)
ew7=lmer(Ejac~Sim+Total+Fam+MW
         +(1|Eb.ID),
         data=ed3)
anova(ew6,ew7)
#Fam by Total non sig
ew8=lmer(Ejac~Sim+Total+Fam
         +(1|Eb.ID),
         data=ed3)
anova(ew7,ew8)
#Keep MW

set_theme(base = theme_sjplot2(),panel.major.gridcol = "white", panel.minor.gridcol = "white",
          panel.gridcol = "white")
MWEjacplot=plot_model(ew7,type=c("pred"),terms=c("MW"),pretty=FALSE,title="",
                      axis.title=c("Scaled Male Weight (mg)","Ejaculate Weight (mg)"),
                      show.data=TRUE)
MWEjacplot

ew9=lmer(Ejac~Sim+Total+MW
         +(1|Eb.ID),
         data=ed3)
anova(ew7,ew9)
#Fam not sig
ew10=lmer(Ejac~Total+Fam+MW
          +(1|Eb.ID),
          data=ed3)
anova(ew7,ew10)
#Sim not sig
ew11=lmer(Ejac~Sim+Fam+MW
          +(1|Eb.ID),
          data=ed3)
anova(ew7,ew11)
#Total sig

###Separate Plots for 1st Order Effects###

Sim2=subset(ed3,Sim=="Two")
Sim4=subset(ed3,Sim=="Four")
MeanSim2=mean(Sim2$Ejac)
seSim2=sd(Sim2$Ejac/sqrt(length(Sim2$Ejac)))
MeanSim4=mean(Sim4$Ejac)
seSim4=sd(Sim4$Ejac/sqrt(length(Sim4$Ejac)))
SimMeans=c(MeanSim2,MeanSim4)
Simse=c(seSim2,seSim4)
SimMeans
Simse

Total2=subset(ed3,Total=="Two")
Total4=subset(ed3,Total=="Four")
MeanTotal2=mean(Total2$Ejac)
seTotal2=sd(Total2$Ejac/sqrt(length(Total2$Ejac)))
MeanTotal4=mean(Total4$Ejac)
seTotal4=sd(Total4$Ejac/sqrt(length(Total4$Ejac)))
TotalMeans=c(MeanTotal2,MeanTotal4)
Totalse=c(seTotal2,seTotal4)
TotalMeans
Totalse

FamF=subset(ed3,Fam=="F")
FamUF=subset(ed3,Fam=="UF")
MeanFamF=mean(FamF$Ejac)
seFamF=sd(FamF$Ejac/sqrt(length(FamF$Ejac)))
MeanFamUF=mean(FamUF$Ejac)
seFamUF=sd(FamUF$Ejac/sqrt(length(FamUF$Ejac)))
FamMeans=c(MeanFamF,MeanFamUF)
Famse=c(seFamF,seFamUF)
FamMeans
Famse

Treatments=c("two","four")
TreatmentsF=c("familiar","unfamiliar")

par(mfrow = c(1,3))
par(mar=c(8,6,5,2))

SimPlot=barplot(SimMeans,names.arg=Treatments,width=1,space=0.2, col=c("skyblue","seagreen"),
                ylim=c(0.15,0.225), beside=TRUE, xpd = FALSE,
                cex.names=2.5,
                cex.axis=2.5,
                cex.lab=2.5)
title(ylab="Mean Ejaculate Weight", line=4,cex.lab=2.5)
title(xlab="Number of Competitor
Males Present At Once", 
      line=6.5, cex.lab=2.5)
box(bty="l")
arrows(SimPlot, SimMeans - Simse, SimPlot,
       SimMeans + Simse, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

TotalPlot=barplot(TotalMeans,names.arg=Treatments,width=1,space=0.2,col=c("skyblue","seagreen"),
                  ylim=c(0.15,0.225), beside=TRUE, xpd = FALSE,
                  cex.names=2.5,
                  cex.axis=2.5,
                  cex.lab=2.5)
title(ylab="Mean Ejaculate Weight", line=4,cex.lab=2.5)
title(xlab="Number of Competitor
Males Present in Total", 
      line=6.5, cex.lab=2.5)
box(bty="l")
arrows(TotalPlot, TotalMeans - Totalse, TotalPlot,
       TotalMeans + Totalse, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
y = 0.2225
x=TotalPlot
offset = 0.002
lines(x[1:2],c(y, y))
lines(x[c(1,1)],c(y, y-offset))
lines(x[c(2,2)],c(y, y-offset))
text(x[1]+((x[2]-x[1])/2),y+offset,"*",cex=3)

FamPlot=barplot(FamMeans,width=1,space=0.2, col=c("skyblue","seagreen"),
                ylim=c(0.15,0.225), beside=TRUE, xpd = FALSE,
                cex.names=2.5,
                cex.axis=2.5,
                cex.lab=2.5,
                names.arg=c("familiar","unfamiliar"))
title(ylab="Mean Ejaculate Weight", line=4,cex.lab=2.5)
title(xlab="Familiarity of Final
Competitor Males", 
      line=6.5, cex.lab=2.5)
box(bty="l")
arrows(FamPlot, FamMeans - Famse, FamPlot,
       FamMeans + Famse, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

###PDF Plot for Submission###

Sim2=subset(ed3,Sim=="Two")
Sim4=subset(ed3,Sim=="Four")
MeanSim2=mean(Sim2$Ejac)
seSim2=sd(Sim2$Ejac)/sqrt(length(Sim2$Ejac))
MeanSim4=mean(Sim4$Ejac)
seSim4=sd(Sim4$Ejac)/sqrt(length(Sim4$Ejac))
SimMeans=c(MeanSim2,MeanSim4)
Simse=c(seSim2,seSim4)
SimMeans
Simse
length(Sim2$Ejac)
length(Sim4$Ejac)

Total2=subset(ed3,Total=="Two")
Total4=subset(ed3,Total=="Four")
MeanTotal2=mean(Total2$Ejac)
seTotal2=sd(Total2$Ejac)/sqrt(length(Total2$Ejac))
MeanTotal4=mean(Total4$Ejac)
seTotal4=sd(Total4$Ejac)/sqrt(length(Total4$Ejac))
TotalMeans=c(MeanTotal2,MeanTotal4)
Totalse=c(seTotal2,seTotal4)
TotalMeans
Totalse
length(Total2$Ejac)
length(Total4$Ejac)

FamF=subset(ed3,Fam=="F")
FamUF=subset(ed3,Fam=="UF")
MeanFamF=mean(FamF$Ejac)
seFamF=sd(FamF$Ejac)/sqrt(length(FamF$Ejac))
MeanFamUF=mean(FamUF$Ejac)
seFamUF=sd(FamUF$Ejac)/sqrt(length(FamUF$Ejac))
FamMeans=c(MeanFamF,MeanFamUF)
Famse=c(seFamF,seFamUF)
FamMeans
Famse
length(FamF$Ejac)
length(FamUF$Ejac)

pdf(file="Figure 1.pdf",width=8,height=4,family="Helvetica")
par(family="serif",font=1)

par(mfrow = c(1,3))
par(mar=c(6.5,5,5,2.5))

SimPlot=barplot(SimMeans,width=1,space=0.2,names.arg=c("Two","Four"),
                col=c("skyblue","seagreen"),
                ylab=c("Mean ejaculate weight (mg)"),
                ylim=c(0.15,0.23), beside=TRUE, xpd = FALSE,
                cex.names=1.6,
                cex.axis=1.6,
                cex.lab=1.6)
title(xlab="Number of competitor males
present at once", 
      line=5, cex.lab=1.8)
box(bty="l")
arrows(SimPlot, SimMeans - Simse, SimPlot,
       SimMeans + Simse, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
mtext("A",3, line=0.5, adj=-0.1, cex=1.6)

TotalPlot=barplot(TotalMeans,width=1,space=0.2,names.arg=c("Two","Four"), 
                  col=c("skyblue","seagreen"),
                  ylab=c("Mean ejaculate weight (mg)"),
                  ylim=c(0.15,0.23), beside=TRUE, xpd = FALSE,
                  cex.names=1.6,
                  cex.axis=1.6,
                  cex.lab=1.6)
title(xlab="Number of competitor males
encountered in total", 
      line=5, cex.lab=1.8)
box(bty="l")
arrows(TotalPlot, TotalMeans - Totalse, TotalPlot,
       TotalMeans + Totalse, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
y = 0.2225
x=TotalPlot
offset = 0.002
lines(x[1:2],c(y, y))
text(x[1]+((x[2]-x[1])/2),y+0.005,"*",cex=2.5)
mtext("B",3, line=0.5, adj=-0.1, cex=1.6)

FamPlot=barplot(FamMeans,width=1,space=0.2,names.arg=c("Familiar","Unfamiliar"),
                col=c("skyblue","seagreen"),
                ylab=c("Mean ejaculate weight (mg)"),
                ylim=c(0.15,0.23), beside=TRUE, xpd = FALSE,
                cex.names=1.6,
                cex.axis=1.6,
                cex.lab=1.6)
title(xlab="Familiarity of final
competitor males", 
      line=5, cex.lab=1.6)
box(bty="l")
arrows(FamPlot, FamMeans - Famse, FamPlot,
       FamMeans + Famse, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
mtext("C",3, line=0.5, adj=-0.1, cex=1.6)

dev.off()


#New plots using the effect function

library(effects)
library(ggplot2)

#Sim
ew7=lmer(Ejac~Sim+Total+Fam+MW
         +(1|Eb.ID),
         data=ed3)

par(mfrow = c(1,3))

ef1=effect("Sim",ew7)
summary(ef1)
sim=as.data.frame(ef1)
sim
ggplot(sim,aes(Sim,fit))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0.4)+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(family="serif",size=rel(1.5),angle=90),
        axis.text=element_text(family="serif",size=rel(1.5),colour="black"))+
  labs(x="",y="mean ejaculate weight (mg)")+
  geom_line(size=1.2)

ef1=effect("Total",ew7)
summary(ef1)
total=as.data.frame(ef1)
total
ggplot(total,aes(Total,fit))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0.4)+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(family="serif",size=rel(1.5),angle=90),
        axis.text=element_text(family="serif",size=rel(1.5),colour="black"))+
  labs(x="",y="mean ejaculate weight (mg)")+
  geom_line(size=1.2)

ef1=effect("Fam",ew7)
summary(ef1)
fam=as.data.frame(ef1)
fam
ggplot(fam,aes(Fam,fit))+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0.2)+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(family="serif",size=rel(1.5),angle=90),
        axis.text=element_text(family="serif",size=rel(1.5),colour="black"))+
  labs(x="",y="mean ejaculate weight (mg)")+
  geom_line(size=1.2)
?theme

####
#geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
 #             position=position_dodge(0.05))

