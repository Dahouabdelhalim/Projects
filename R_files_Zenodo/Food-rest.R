rm(list=ls())

library(plotrix)
library(car)
library(multcomp)
library(MASS)
library(nlme)
library(lme4)
library(lmerTest)



data.rest = read.csv2("Ovaries_Foodrest.csv", header=T)
data.rest$Treatment = as.factor(data.rest$Treat)
data.rest$Ovmean = as.numeric(data.rest$Ovmean)
data.rest$white_egg = as.numeric(data.rest$white_egg)

# ovary length
###############
wilcox.test(Ovmean~Treat, data=data.rest)

moyenne <- aggregate(Ovmean~Treat,FUN="mean",data=data.rest);names(moyenne)[names(moyenne)=="moy"]<-"mean"
std.error <- aggregate(Ovmean~Treat,FUN="std.error",data=data.rest);names(std.error)[names(std.error)=="se"]<-"std.error"
summary_data <-merge(moyenne,std.error,all.x=T,all.y=T)
summary_data

mea.Wtot = c(3.87,3.49)
se.Wtot = c(0.127,0.18)
patric= c("Ctrl","Treat")
Wegg.tr = barplot(mea.Wtot, ylab="Ovary length (mm)" ,names = patric, cex.lab=1.1, cex.names=1.1, col=c("grey","white"), space=0, ylim=c(0,6))
arrows(Wegg.tr,mea.Wtot-se.Wtot,Wegg.tr, mea.Wtot+se.Wtot, lwd=2, angle=90,length=0.1,code=3)
text(x=1, y=4.5, labels="*",font=3, ps=0.5,cex = 2)
# number of developping white egg
##################################

lmm1 = glm(white_egg~Treat, family = "poisson",data=data.rest)
Anova(lmm1)

# number of laid egg 
#######################

data.rest = read.csv2("Egg_Foodrest.csv", header=T)
data.rest$Treatment = as.factor(data.rest$Treat)
data.rest$Egg = as.numeric(data.rest$Egg)
data.rest$Week = as.numeric(data.rest$Week)

moyenne <- aggregate(Egg~Week*Treat,FUN=mean,data=data.rest);names(moyenne)[names(moyenne)=="Egg"]<-"mean"
stderror <- aggregate(Egg~Week*Treat,FUN=std.error,data=data.rest);names(stderror)[names(stderror)=="Egg"]<-"std.error"
summary_data <-merge(moyenne,stderror,all.x=T,all.y=T)


ymin <- min(summary_data$mean-summary_data$std.error,na.rm=T) - 0.1* (max(summary_data$mean+summary_data$std.error,na.rm=T)-min(summary_data$mean-summary_data$std.error,na.rm=T))
ymax <- max(summary_data$mean+summary_data$std.error,na.rm=T) + 0.1 *(max(summary_data$mean+summary_data$std.error,na.rm=T)-min(summary_data$mean-summary_data$std.error,na.rm=T))


plot(mean~Week,data=summary_data[summary_data$Treat=="Restriction",],type="p",col="black",bg="black",ylim=c(ymin,ymax),pch=22,xlab="Time (week)",ylab="Egg number",cex.lab=0.9,cex.main=0.9,bty="l")
arrows(summary_data[summary_data$Treat=="Restriction","Week"],summary_data[summary_data$Treat=="Restriction","mean"]-summary_data[summary_data$Treat=="Restriction","std.error"],summary_data[summary_data$Treat=="Restriction","Week"],summary_data[summary_data$Treat=="Restriction","mean"]+summary_data[summary_data$Treat=="Restriction","std.error"],col="black",angle=90,code=3,length=0.05)
arrows(summary_data[summary_data$Treat=="Crontrol","Week"],summary_data[summary_data$Treat=="Crontrol","mean"]-summary_data[summary_data$Treat=="Crontrol","std.error"],summary_data[summary_data$Treat=="Crontrol","Week"],summary_data[summary_data$Treat=="Crontrol","mean"]+summary_data[summary_data$Treat=="Crontrol","std.error"],col="black",angle=90,code=3,length=0.05)
points(mean~Week,data=summary_data[summary_data$Treat=="Crontrol",],type="p",col="black",bg="grey",pch=21)


# Final number of laid Egg
##########################

data.final= data.rest[data.rest$Week == 12,] 
lmm1 = glm(Egg~Treat, family = "poisson",data=data.final)
Anova(lmm1)


# survival
############

data.rest = read.csv2("Surv_Foodrest.csv", header=T)

data.rest
data.rest$Treatment = as.factor(data.rest$Treat)
data.rest$time = as.numeric(data.rest$time)
data.rest$Death = as.numeric(data.rest$Death)

library(survival)


survival_model<- survreg((Surv(time=time,event=Death)~Treatment),data=data.rest)
summary(survival_model)
anova(survival_model)

