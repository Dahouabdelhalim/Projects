## Heritability of all secondary traits with graphs
## Prepared by Mokhles Rahman
## March 23, 2020
## mrahman@ksu.edu
## Purpose to document five years data anlysis
rm(list=ls()) 
cat("\\f")
## Setting working directory
setwd("~/Documents/BHEARD_documents/Dissertation_research/Data_Analysis")
### loading packages
require(pastecs) 
require("lme4")
require(ggplot2)
library(leaps)
library(caret) 
library(dplyr)

### Loading data and setting entry, rep, range, and trial as factors
data_2016 <- read.csv("2016_Data.csv")
data_2016$entry <- as.factor(data_2016$entry)
data_2016$rep <- as.factor(data_2016$rep)
data_2016$range <- as.factor(data_2016$range)
data_2016$trial <- as.factor(data_2016$trial) #note when you change this to a factor it effects the next part

## Estimating heritability considering all trials together
# require(lme4)

heritab2016=data.frame(Traits=factor(),Heritability=double())
for(q in 7:ncol(data_2016)){
        mod = lmer(data_2016[,q]~ 1 + (1|trial) + (1|trial:entry) + (1|trial:rep) + (1|trial:rep:range), data=data_2016)
        e=data.frame(VarCorr(mod))[1,4]
        r=data.frame(VarCorr(mod))[5,4]
        H2=e/(e+(r/2))
        H2=round(H2,digits=2)
        T = names(data_2016[q])
        heritab2016=rbind(heritab2016,data.frame(Traits=T, Heritability=H2))
}
heritab2016
# write.csv(heritab2016,file="Heritability_allTrials_together_2016.csv",row.names=FALSE,quote=FALSE)

## Barplot
# heritab2016$Heritability=as.numeric(as.character(heritab2016$Heritability))
# heritab2016$Traits=as.character(heritab2016$Traits)
# par(mgp=c(3,0.5,0))
# par(mar=c(8,2,1.5,0))
# barplot(heritab2016$Heritability, names.arg = heritab2016$Traits,las=2,ylim=c(0,1.0),main="Heritability 2016")
# abline(v=7,col="blue",lty=1,lwd=5)

# dev.off()
# par("mar")
# par("mgp")
par(mar=c(3,3,1,1))
par(mgp=c(1.5,0.5,0))
par(mfrow=c(3,2))

aa=heritab2016[1:8,]
aa=as.matrix(aa)
aa=aa[,2]
aa=as.numeric(aa)
DAS=c("59","68","78","90","95","97","104","110")
DAS=as.numeric(DAS)
plot(DAS,aa,ylim=c(0,1),xlim=c(55,115),ylab="Heritability",xlab="Days after sowing",type="p",col="red")
lines(lowess(DAS,aa,f=2/3),col="red")
abline(h=heritab2016$Heritability[heritab2016$Traits=="GRYLD"],lty=3)
abline(v=70,lty=5)
abline(v=104,lty=3)
title("A",adj=0)
par(new=TRUE)

ab=heritab2016[9:17,]
ab=as.matrix(ab)
ab=ab[,2]
ab=as.numeric(ab)
DAS=c("57","61","66","74","90","95","98","105","110")
DAS=as.numeric(DAS)
plot(DAS,ab,ylim=c(0,1),xlim=c(55,115),ylab="Heritability",xlab="Days after sowing",type="p",col="blue")
lines(lowess(DAS,ab,f=2/3),col="blue")
legend("topleft",c("CT","NDVI"),fill=c("red","blue"),bty="n")

# ac=data.frame(heritab2016[18:26,])
# ac$Heritability=as.numeric(as.character(ac$Heritability))
# barplot(ac$Heritability, names.arg = ac$Traits, las=2, ylim=c(0,1.0), ylab = "Heritability",bty="n", col = "green4")

## Estimating heritability by trial basis
heritabbytrial2016=data.frame()
for(trial in 1:10){
        data.per.trial=data_2016[as.character(data_2016$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        for(q in 7:ncol(data.per.trial)){
                mod = lmer(data.per.trial[,q]~ 1 + (1|entry) + (1|rep:range), data=data.per.trial)
                e=data.frame(VarCorr(mod))[1,4]
                r=data.frame(VarCorr(mod))[3,4]
                Heritability=e/(e+(r/2))
                Heritability=round(Heritability,digits=2)
                Traits=names(data.per.trial[q])
                Trial=data.per.trial$trial[trial]
                herit=cbind(Trial,Traits,Heritability)
                heritabbytrial2016=rbind(heritabbytrial2016,herit)
        }
}
# write.csv(heritabbytrial2016,file="Heritability_byTrial_2016_long.csv",row.names=FALSE,quote=FALSE)

H=heritabbytrial2016$Heritability
H=matrix(H,nrow = 26,ncol = 10)
Traits=colnames(BLUE2016)[3:ncol(BLUE2016)]
H=data.frame(Traits,H)
colnames(H)=c("Traits","Trial_1","Trial_2","Trial_3","Trial_4","Trial_5","Trial_6","Trial_7","Trial_8","Trial_9","Trial_10")
write.csv(H,file="Heritability_byTrial_2016.csv",row.names=FALSE,quote=FALSE)


### 2017 data
data_2017 <- read.csv("2017_Data.csv")
data_2017$entry <- as.factor(data_2017$entry)
data_2017$rep <- as.factor(data_2017$rep)
data_2017$range <- as.factor(data_2017$range)
data_2017$trial <- as.factor(data_2017$trial) #note when you change this to a factor it effects the next part

## Estimating heritability considering all trials together
heritab2017=data.frame(Traits=factor(),Heritability=double())
for(q in 7:ncol(data_2017)){
        mod = lmer(data_2017[,q]~ 1 + (1|trial) + (1|trial:entry) + (1|trial:rep) + (1|trial:rep:range), data=data_2017)
        e=data.frame(VarCorr(mod))[1,4]
        r=data.frame(VarCorr(mod))[5,4]
        H2=e/(e+(r/2))
        H2=round(H2,digits=2)
        T = names(data_2017[q])
        heritab2017=rbind(heritab2017,data.frame(Traits=T, Heritability=H2))
}
heritab2017
# write.csv(heritab2017,file="Heritability_allTrials_together_2017.csv",row.names=FALSE,quote=FALSE)

# heritab2017$Heritability=as.numeric(as.character(heritab2017$Heritability))
# heritab2017$Traits=as.character(heritab2017$Traits)
# par(mgp=c(3,0.5,0))
# par(mar=c(8,2,1.5,0))
# barplot(heritab2017$Heritability,names.arg = heritab2017$Traits,las=2,ylim=c(0,1.0),main="Heritability 2017")
# dev.off()

aa=heritab2017[1:14,]
aa=as.matrix(aa)
aa=aa[,2]
aa=as.numeric(aa)
DAS=c("38","43","48","54","59","65","70","75","80","86","90","95","100","106")
DAS=as.numeric(DAS)
plot(DAS,aa,ylim=c(0,1),xlim=c(35,105),ylab="Heritability",xlab="Days after sowing",type="p",col="red")
# lines(spline(DAS,aa),col="red")
lines(lowess(DAS,aa,f=2/3),col="red")
abline(h=heritab2017$Heritability[heritab2017$Traits=="GRYLD"],lty=3)
abline(v=70,lty=5)
abline(v=105,lty=3)
title("B",adj=0)
par(new=TRUE)

ab=heritab2017[15:28,]
ab=as.matrix(ab)
ab=ab[,2]
ab=as.numeric(ab)
DAS=c("37","42","48","54","59","65","70","75","80","85","90","95","100","106")
DAS=as.numeric(DAS)
plot(DAS,ab,ylim=c(0,1),xlim=c(35,105),ylab="Heritability",xlab="Days after sowing",type="p",col="blue")
lines(lowess(DAS,ab,f=2/3),col="blue")
legend("topleft",c("CT","NDVI"),fill=c("red","blue"),bty="n")

# ac=data.frame(heritab2017[29:37,])
# ac$Heritability=as.numeric(as.character(ac$Heritability))
# barplot(ac$Heritability, names.arg = ac$Traits, las=2, ylim=c(0,1.0), ylab = "Heritability", col = "green4")

## Estimating heritability by trial basis
heritabbytrial2017=data.frame()
for(trial in 1:11){
        data.per.trial=data_2017[as.character(data_2017$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        for(q in 7:ncol(data.per.trial)){
                mod = lmer(data.per.trial[,q]~ 1 + (1|entry) + (1|rep:range), data=data.per.trial)
                e=data.frame(VarCorr(mod))[1,4]
                r=data.frame(VarCorr(mod))[3,4]
                Heritability=e/(e+(r/2))
                Heritability=round(Heritability,digits=2)
                Traits=names(data.per.trial[q])
                Trial=data.per.trial$trial[trial]
                herit=cbind(Trial,Traits,Heritability)
                heritabbytrial2017=rbind(heritabbytrial2017,herit)
        }
}
# write.csv(heritabbytrial2017,file="Heritability_byTrial_2017.csv",row.names=FALSE,quote=FALSE)
H2017=heritabbytrial2017$Heritability
H2017=matrix(H2017,nrow = 37,ncol = 11)
Traits=colnames(BLUE2017)[3:ncol(BLUE2017)]
H2017=data.frame(Traits,H2017)
colnames(H2017)=c("Traits","Trial_1","Trial_2","Trial_3","Trial_4","Trial_5","Trial_6","Trial_7","Trial_8","Trial_9","Trial_10","Trial_11")
write.csv(H2017,file="Heritability_byTrial_2017.csv",row.names=FALSE,quote=FALSE)

### 2018 data
data_2018 <- read.csv("2018_Data.csv")
data_2018$entry <- as.factor(data_2018$entry)
data_2018$rep <- as.factor(data_2018$rep)
data_2018$range <- as.factor(data_2018$range)
data_2018$trial <- as.factor(data_2018$trial) #note when you change this to a factor it effects the next part


## Estimating heritability considering all trials together
# require(lme4)

heritab2018=data.frame(Traits=factor(),Heritability=double())
for(q in 7:ncol(data_2018)){
        mod = lmer(data_2018[,q]~ 1 + (1|trial) + (1|trial:entry) + (1|trial:rep) + (1|trial:rep:range), data=data_2018)
        e=data.frame(VarCorr(mod))[1,4]
        r=data.frame(VarCorr(mod))[5,4]
        H2=e/(e+(r/2))
        H2=round(H2,digits=2)
        T = names(data_2018[q])
        heritab2018=rbind(heritab2018,data.frame(Traits=T, Heritability=H2))
}
# write.csv(heritab2018,file="Heritability_allTrials_together_2018.csv",row.names=FALSE,quote=FALSE)

# heritab2018$Heritability=as.numeric(as.character(heritab2018$Heritability))
# heritab2018$Traits=as.character(heritab2018$Traits)
# pdf(file="Heritability 2018 latest.pdf",width=9,height=4)
# par(mgp=c(3,0.5,0))
# par(mar=c(8,2,1.5,0))
# barplot(heritab2018$Heritability,names.arg = heritab2018$Traits,las=2,ylim=c(0,1.0),main="Heritability 2018")
# dev.off()

aa=heritab2018[1:12,]
aa=as.matrix(aa)
aa=aa[,2]
aa=as.numeric(aa)
DAS=c("58","63","68","73","77","82","88","92","96","101","106","111")
DAS=as.numeric(DAS)
plot(DAS,aa,ylim=c(0,1),xlim=c(55,115),ylab="Heritability",xlab="Days after sowing",type="p",col="red")
lines(lowess(DAS,aa,f=2/3),col="red")
abline(h = heritab2018$Heritability[heritab2018$Traits=="GRYLD"], lty=3)
abline(v=72,lty=5)
abline(v=105,lty=3)
title("C",adj=0)
par(new=TRUE)


ab=heritab2017[13:24,]
ab=as.matrix(ab)
ab=ab[,2]
ab=as.numeric(ab)
DAS=c("58","63","68","72","77","83","88","92","96","101","106","111")
DAS=as.numeric(DAS)
plot(DAS,ab,ylim=c(0,1),xlim=c(55,115),ylab="Heritability",xlab="Days after sowing",type="p",col="blue")
lines(lowess(DAS,ab,f=2/3),col="blue")
legend("topleft",c("CT","NDVI"),fill=c("red","blue"),bty="n")

# ac=data.frame(heritab2018[25:33,])
# ac$Heritability=as.numeric(as.character(ac$Heritability))
# barplot(ac$Heritability, names.arg = ac$Traits, las=2, ylim=c(0,1.0), ylab = "Heritability", col = "green4")

## Estimating heritability by trial basis
heritabbytrial2018=data.frame()
for(trial in 1:11){
        data.per.trial=data_2018[as.character(data_2018$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        for(q in 7:ncol(data.per.trial)){
                mod = lmer(data.per.trial[,q]~ 1 + (1|entry) + (1|rep:range), data=data.per.trial)
                e=data.frame(VarCorr(mod))[1,4]
                r=data.frame(VarCorr(mod))[3,4]
                Heritability=e/(e+(r/2))
                Heritability=round(Heritability,digits=2)
                Traits=names(data.per.trial[q])
                Trial=data.per.trial$trial[trial]
                herit=cbind(Trial,Traits,Heritability)
                heritabbytrial2018=rbind(heritabbytrial2018,herit)
        }
}
# write.csv(heritabbytrial2018,file="Heritability_byTrial_2018.csv",row.names=FALSE,quote=FALSE)
H2018=heritabbytrial2018$Heritability
H2018=matrix(H2018,nrow = 33,ncol = 11)
Traits=colnames(BLUE2018)[3:ncol(BLUE2018)]
H2018=data.frame(Traits,H2018)
colnames(H2018)=c("Traits","Trial_1","Trial_2","Trial_3","Trial_4","Trial_5","Trial_6","Trial_7","Trial_8","Trial_9","Trial_10","Trial_11")
write.csv(H2018,file="Heritability_byTrial_2018.csv",row.names=FALSE,quote=FALSE)


### 2019 data
data_2019 <- read.csv("2019_Data.csv")
data_2019$entry <- as.factor(data_2019$entry)
data_2019$rep <- as.factor(data_2019$rep)
data_2019$range <- as.factor(data_2019$range)
data_2019$trial <- as.factor(data_2019$trial) #note when you change this to a factor it effects the next part


## Estimating heritability considering all trials together
# require(lme4)
heritab2019=data.frame(Traits=factor(),Heritability=double())
for(q in 7:ncol(data_2019)){
        mod = lmer(data_2019[,q]~ 1 + (1|trial) + (1|trial:entry) + (1|trial:rep) + (1|trial:rep:range), data=data_2019)
        e=data.frame(VarCorr(mod))[1,4]
        r=data.frame(VarCorr(mod))[5,4]
        H2=e/(e+(r/2))
        H2=round(H2,digits=2)
        T = names(data_2019[q])
        heritab2019=rbind(heritab2019,data.frame(Traits=T, Heritability=H2))
}
# write.csv(heritab2019,file="Heritability_allTrials_together_2019.csv",row.names=FALSE,quote=FALSE)

# heritab2019$Heritability=as.numeric(as.character(heritab2019$Heritability))
# heritab2019$Traits=as.character(heritab2019$Traits)
# # pdf(file="Heritability 2019 latest.pdf",width=9,height=4)
# par(mgp=c(3,0.5,0))
# par(mar=c(8,2,1.5,0))
# barplot(heritab2019$Heritability,names.arg = heritab2019$Traits,las=2,ylim=c(0,1.0),main="Heritability 2019")
# # dev.off()

aa=heritab2019[1:13,]
aa=as.matrix(aa)
aa=aa[,2]
aa=as.numeric(aa)
DAS=c("56","60","64","69","75","82","87","93","97","103","108","112","117")
DAS=as.numeric(DAS)
plot(DAS,aa,ylim=c(0,1),xlim=c(55,120),ylab="Heritability",xlab="Days after sowing",type="p",col="red")
# lines(spline(DAS,aa),col="red")
lines(lowess(DAS,aa,f=2/3),col="red")
abline(h = heritab2019$Heritability[heritab2019$Traits=="GRYLD"], lty=3)
abline(v=77,lty=5)
abline(v=110,lty=3)
title("D",adj=0)
par(new=TRUE)


ab=heritab2019[14:26,]
ab=as.matrix(ab)
ab=ab[,2]
ab=as.numeric(ab)
DAS=c("54","60","64","69","75","82","86","92","97","103","107","112","117")
DAS=as.numeric(DAS)
plot(DAS,ab,ylim=c(0,1),xlim=c(55,120),ylab="Heritability",xlab="Days after sowing",type="p",col="blue")
lines(lowess(DAS,ab,f=2/3),col="blue")
legend("topleft",c("CT","NDVI"),fill=c("red","blue"),bty="n")

# ac=data.frame(heritab2019[27:35,])
# ac$Heritability=as.numeric(as.character(ac$Heritability))
# barplot(ac$Heritability, names.arg = ac$Traits, las=2, ylim=c(0,1.0), ylab = "Heritability", col = "green4")

## Estimating heritability by trial basis
heritabbytrial2019=data.frame()
for(trial in 1:10){
        data.per.trial=data_2019[as.character(data_2019$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        for(q in 7:ncol(data.per.trial)){
                mod = lmer(data.per.trial[,q]~ 1 + (1|entry) + (1|rep:range), data=data.per.trial)
                e=data.frame(VarCorr(mod))[1,4]
                r=data.frame(VarCorr(mod))[3,4]
                Heritability=e/(e+(r/2))
                Heritability=round(Heritability,digits=2)
                Traits=names(data.per.trial[q])
                Trial=data.per.trial$trial[trial]
                herit=cbind(Trial,Traits,Heritability)
                heritabbytrial2019=rbind(heritabbytrial2019,herit)
        }
}
# write.csv(heritabbytrial2019,file="Heritability_byTrial_2019.csv",row.names=FALSE,quote=FALSE)
H2019=heritabbytrial2019$Heritability
H2019=matrix(H2019,nrow = 35,ncol = 10)
Traits=colnames(BLUE2019)[3:ncol(BLUE2019)]
H2019=data.frame(Traits,H2019)
colnames(H2019)=c("Traits","Trial_1","Trial_2","Trial_3","Trial_4","Trial_5","Trial_6","Trial_7","Trial_8","Trial_9","Trial_10")
write.csv(H2019,file="Heritability_byTrial_2019.csv",row.names=FALSE,quote=FALSE)


### 2020 data
data_2020 <- read.csv("2020_Data.csv")
data_2020 <- data_2020[,-c(39,40)]
data_2020$entry <- as.factor(data_2020$entry)
data_2020$rep <- as.factor(data_2020$rep)
data_2020$range <- as.factor(data_2020$range)
data_2020$trial <- as.factor(data_2020$trial) #note when you change this to a factor it effects the next part


## Estimating heritability considering all trials together
# require(lme4)
heritab2020=data.frame(Traits=factor(),Heritability=double())
for(q in 7:ncol(data_2020)){
        mod = lmer(data_2020[,q]~ 1 + (1|trial) + (1|trial:entry) + (1|trial:rep) + (1|trial:rep:range), data=data_2020)
        e=data.frame(VarCorr(mod))[1,4]
        r=data.frame(VarCorr(mod))[5,4]
        H2=e/(e+(r/2))
        H2=round(H2,digits=2)
        T = names(data_2020[q])
        heritab2020=rbind(heritab2020,data.frame(Traits=T, Heritability=H2))
}
# write.csv(heritab2020,file="Heritability_allTrials_together_2020.csv",row.names=FALSE,quote=FALSE)

# heritab2020$Heritability=as.numeric(as.character(heritab2020$Heritability))
# heritab2020$Traits=as.character(heritab2020$Traits)
# pdf(file="Heritability 2020 latest.pdf",width=9,height=4)
# par(mgp=c(3,0.5,0))
# par(mar=c(8,2,1.5,0))
# barplot(heritab2020$Heritability,names.arg = heritab2020$Traits,las=2,ylim=c(0,1.0),main="Heritability 2020")
# dev.off()

aa=heritab2020[1:15,]
aa=as.matrix(aa)
aa=aa[,2]
aa=as.numeric(aa)
DAS=c("42","46","51","56","60","65","70","75","80","86","91","97","102","107","112")
DAS=as.numeric(DAS)
plot(DAS,aa,ylim=c(0,1),xlim=c(40,115),ylab="Heritability",xlab="Days after sowing",type="p",col="red")
# lines(spline(DAS,aa),col="red")
lines(lowess(DAS,aa,f=2/3),col="red")
abline(h = heritab2020$Heritability[heritab2020$Traits=="GRYLD"], lty=3)
abline(v=73,lty=5)
abline(v=104,lty=3)
title("E",adj=0)
par(new=TRUE)


ab=heritab2020[16:30,]
ab=as.matrix(ab)
ab=ab[,2]
ab=as.numeric(ab)
DAS=c("42","46","51","56","60","65","70","75","80","86","91","97","102","107","112")
DAS=as.numeric(DAS)
plot(DAS,ab,ylim=c(0,1),xlim=c(40,115),ylab="Heritability",xlab="Days after sowing",type="p",col="blue")
lines(lowess(DAS,ab,f=2/3),col="blue")
legend("topleft",c("CT","NDVI"),fill=c("red","blue"),bty="n")

# ac=data.frame(heritab2020[31:43,])
# ac$Heritability=as.numeric(as.character(ac$Heritability))
# barplot(ac$Heritability, names.arg = ac$Traits, las=2, ylim=c(0,1.0), ylab = "Heritability", col = "green4")

## Estimating heritability by trial basis
heritabbytrial2020=data.frame()
for(trial in 1:11){
        data.per.trial=data_2020[as.character(data_2020$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        for(q in 7:ncol(data.per.trial)){
                mod = lmer(data.per.trial[,q]~ 1 + (1|entry) + (1|rep:range), data=data.per.trial)
                e=data.frame(VarCorr(mod))[1,4]
                r=data.frame(VarCorr(mod))[3,4]
                Heritability=e/(e+(r/2))
                Heritability=round(Heritability,digits=2)
                Traits=names(data.per.trial[q])
                Trial=data.per.trial$trial[trial]
                herit=cbind(Trial,Traits,Heritability)
                heritabbytrial2020=rbind(heritabbytrial2020,herit)
        }
}
# write.csv(heritabbytrial2020,file="Heritability_byTrial_2020.csv",row.names=FALSE,quote=FALSE)

H2020=heritabbytrial2020$Heritability
H2020=matrix(H2020,nrow = 42,ncol = 11)
Traits=colnames(BLUE2020)[3:ncol(BLUE2020)]
H2020=data.frame(Traits,H2020)
colnames(H2020)=c("Traits","Trial_1","Trial_2","Trial_3","Trial_4","Trial_5","Trial_6","Trial_7","Trial_8","Trial_9","Trial_10","Trial_11")
write.csv(H2020,file="Heritability_byTrial_2020.csv",row.names=FALSE,quote=FALSE)

heritab.allyears.alltrials=merge(heritab2016,heritab2017,heritab2018,heritab2019,heritab2020)
heritab.allyears.alltrials


# par("mar")
# pdf(file = 'Figure 1.pdf', height = 18, width = 12)
# par(mfrow=c(5,1))
# par(mgp=c(3,0.5,0))
# par(mar=c(8,2,1.5,0))
# dev.off()
# barplot(heritab2016$Heritability,names.arg = heritab2016$Traits,las=2,ylim=c(0,1.0))
# title("A",adj=0)
# barplot(heritab2017$Heritability,names.arg = heritab2017$Traits,las=2,ylim=c(0,1.0))
# title("B",adj=0)
# barplot(heritab2018$Heritability,names.arg = heritab2018$Traits,las=2,ylim=c(0,1.0))
# title("C",adj=0)
# barplot(heritab2019$Heritability,names.arg = heritab2019$Traits,las=2,ylim=c(0,1.0))
# title("D",adj=0)
# barplot(heritab2020$Heritability,names.arg = heritab2020$Traits,las=2,ylim=c(0,1.0))
# title("E",adj=0)
# dev.off()
