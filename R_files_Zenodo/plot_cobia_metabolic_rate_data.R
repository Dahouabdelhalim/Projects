setwd("~/Documents/PhD Project/Respirometry/Respirometry")
#sumdat<-read.csv("MO2_methods5_cobia_20180914.csv")
sumdat<-read.csv("MO2_methods10_cobia_20180914.csv") #issues are that C30_28 and C34_28 crash super high cause there are low points throughout hypoxia that screw it

#assign NAs to Ccrit and Scrit for 30_28 and 34_28 because the automation method to calculate ccrit from SMR was not accurate for those two fish
sumdat$Ccrit[which(sumdat$Animal_ID=="C30"&sumdat$Test_Temp=="28")]<-NA
sumdat$Scrit[which(sumdat$Animal_ID=="C30"&sumdat$Test_Temp=="28")]<-NA
sumdat$Ccrit[which(sumdat$Animal_ID=="C34"&sumdat$Test_Temp=="28")]<-NA
sumdat$Scrit[which(sumdat$Animal_ID=="C34"&sumdat$Test_Temp=="28")]<-NA


########### Skip to line 275 to start
#unless you want to plot raw data of MMR and mRMR of all sharks

sumdat<-subset(sumdat, select=c(Animal_ID,Test_Temp,Animal_Mass, MMR,SMR,AS,FS,Scrit,Ccrit))

#get rid of two trials at ESL
#sumdat<-sumdat[-which(sumdat$Animal_ID=="C01" | sumdat$Animal_ID=="C02"),]

#get rid of trials where fish died during trial after chase
#those fish were often outliers for MMR, SMR, AS, and FS
#trials C31_28, C24_32, C28_32, C32_32
sumdat<-sumdat[-which(sumdat$Animal_ID=="C24" & sumdat$Test_Temp=="32" | 
                        sumdat$Animal_ID=="C28" & sumdat$Test_Temp=="32" |
                        sumdat$Animal_ID=="C32" & sumdat$Test_Temp=="32" |
                 sumdat$Animal_ID=="C31" & sumdat$Test_Temp=="28"),]


#take means and SE of each metric
meansAS<-data.frame(Temp=c(24,28,32))
meansAS$MMR<-tapply(sumdat$MMR,sumdat$Test_Temp,mean)
meansAS$SE_MMR<-tapply(sumdat$MMR,sumdat$Test_Temp,sd)/sqrt(tapply(sumdat$MMR,sumdat$Test_Temp,length))
meansAS$SMR<-tapply(sumdat$SMR,sumdat$Test_Temp,mean)
meansAS$SE_SMR<-tapply(sumdat$SMR,sumdat$Test_Temp,sd)/sqrt(tapply(sumdat$SMR,sumdat$Test_Temp,length))
meansAS$AS<-tapply(sumdat$AS,sumdat$Test_Temp,mean)
meansAS$SE_AS<-tapply(sumdat$AS,sumdat$Test_Temp,sd)/sqrt(tapply(sumdat$AS,sumdat$Test_Temp,length))
meansAS$FS<-tapply(sumdat$FS,sumdat$Test_Temp,mean)
meansAS$SE_FS<-tapply(sumdat$FS,sumdat$Test_Temp,sd)/sqrt(tapply(sumdat$FS,sumdat$Test_Temp,length))
meansAS$Scrit<-tapply(sumdat$Scrit,sumdat$Test_Temp,mean,na.rm=T)
meansAS$SE_Scrit<-tapply(sumdat$Scrit,sumdat$Test_Temp,sd,na.rm=T)/sqrt(tapply(sumdat$Scrit,sumdat$Test_Temp,length))
meansAS$Ccrit<-tapply(sumdat$Ccrit,sumdat$Test_Temp,mean,na.rm=T)
meansAS$SE_Ccrit<-tapply(sumdat$Ccrit,sumdat$Test_Temp,sd,na.rm=T)/sqrt(tapply(sumdat$Ccrit,sumdat$Test_Temp,length))


#MMR & RMR plot
par(mfrow=c(1,2))
temps<-c(24,28,32)
plot(MMR~temps, data=meansAS, type="b", ylim=c(50,350),ylab="Metabolic rate (mg/kg/h)",xlab="Temperature (°C)",col="red")
lines(SMR~temps,data=meansAS, type="b",col="blue")
arrows(y0=(meansAS$MMR-meansAS$SE_MMR),y1=(meansAS$MMR+meansAS$SE_MMR), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1, col="red")
arrows(y0=(meansAS$SMR-meansAS$SE_SMR),y1=(meansAS$SMR+meansAS$SE_SMR), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1, col="blue")
legend("topleft",legend=c("MMR","RMR"), fill=c("red","blue"))

#MMR SMR boxplots on same plot
boxplot(MMR~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Metabolic rate (mg/kg/h)",ylim=c(75,400))
n <- 3
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
y0s <- c(mean(sumdat$MMR[which(sumdat$Test_Temp=="24")], na.rm=T), mean(sumdat$MMR[which(sumdat$Test_Temp=="28")], na.rm=T),mean(sumdat$MMR[which(sumdat$Test_Temp=="32")], na.rm=T))
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red",lwd=2)
par(new=T)
boxplot(SMR~Test_Temp, data=sumdat, xlab="", ylab="",ylim=c(75,300),xaxt="n",yaxt="n", col="cyan3",outcol="cyan3",whiskcol="cyan3", staplecol="cyan3",medcol="dark blue")
n <- 3
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
y0s <- c(mean(sumdat$SMR[which(sumdat$Test_Temp=="24")], na.rm=T), mean(sumdat$SMR[which(sumdat$Test_Temp=="28")], na.rm=T),mean(sumdat$SMR[which(sumdat$Test_Temp=="32")], na.rm=T))
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red",lwd=2)
legend("topleft",legend=c("MMR","SMR"),fill=c("white","cyan3"))


#MMR
par(mfrow=c(1,2))
boxplot(MMR~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Metabolic rate (mg/kg/h)",ylim=c(150,400))
n <- 3
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
y0s <- c(mean(sumdat$MMR[which(sumdat$Test_Temp=="24")], na.rm=T), mean(sumdat$MMR[which(sumdat$Test_Temp=="28")], na.rm=T),mean(sumdat$MMR[which(sumdat$Test_Temp=="32")], na.rm=T))
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red",lwd=2)

barplot<-barplot(meansAS$MMR, names.arg=meansAS$Temp, ylim=c(0,400), xlab="Temperature (°C)")
arrows(y0=meansAS$MMR+meansAS$SE_MMR,y1=meansAS$MMR-meansAS$SE_MMR,
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)

#SMR
par(mfrow=c(1,2))
boxplot(SMR~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Metabolic rate (mg/kg/h)",ylim=c(50,200))
n <- 3
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
y0s <- c(mean(sumdat$SMR[which(sumdat$Test_Temp=="24")], na.rm=T), mean(sumdat$SMR[which(sumdat$Test_Temp=="28")], na.rm=T),mean(sumdat$SMR[which(sumdat$Test_Temp=="32")], na.rm=T))
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red",lwd=2)

barplot<-barplot(meansAS$SMR, names.arg=meansAS$Temp, ylim=c(0,200), xlab="Temperature (°C)")
arrows(y0=meansAS$SMR+meansAS$SE_SMR,y1=meansAS$SMR-meansAS$SE_SMR,
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)


#Aerobic Scope plot
par(mfrow=c(1,2))
boxplot(AS~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Metabolic rate (mg/kg/h)",ylim=c(50,250))
n <- 3
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
y0s <- c(mean(sumdat$AS[which(sumdat$Test_Temp=="24")], na.rm=T), mean(sumdat$AS[which(sumdat$Test_Temp=="28")], na.rm=T),mean(sumdat$AS[which(sumdat$Test_Temp=="32")], na.rm=T))
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red",lwd=2)

barplot<-barplot(meansAS$AS, names.arg=meansAS$Temp, ylim=c(0,200), xlab="Temperature (°C)")
arrows(y0=meansAS$AS+meansAS$SE_AS,y1=meansAS$AS-meansAS$SE_AS,
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)

#Factorial Scope plot
par(mfrow=c(1,2))
boxplot(FS~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Factorial scope",ylim=c(1,3))
n <- 3
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
y0s <- c(mean(sumdat$FS[which(sumdat$Test_Temp=="24")], na.rm=T), mean(sumdat$FS[which(sumdat$Test_Temp=="28")], na.rm=T),mean(sumdat$FS[which(sumdat$Test_Temp=="32")], na.rm=T))
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red",lwd=2)

barplot<-barplot(meansAS$FS, names.arg=meansAS$Temp, ylim=c(0,3), xlab="Temperature (°C)")
arrows(y0=meansAS$FS+meansAS$SE_FS,y1=meansAS$FS-meansAS$SE_FS,
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)

#Scrit plots
boxplot(Scrit~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Critical Oxygen Saturation (%)", ylim=c(15,60))
n <- 3
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
y0s <- c(mean(sumdat$Scrit[which(sumdat$Test_Temp=="24")], na.rm=T), mean(sumdat$Scrit[which(sumdat$Test_Temp=="28")], na.rm=T),mean(sumdat$Scrit[which(sumdat$Test_Temp=="32")], na.rm=T))
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red",lwd=2)

barplot<-barplot(meansAS$Scrit, names.arg=meansAS$Temp, ylim=c(0,45), xlab="Temperature (°C)")
arrows(y0=meansAS$Scrit+meansAS$SE_Scrit,y1=meansAS$Scrit-meansAS$SE_Scrit,
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)

#Ccrit plots
boxplot(Ccrit~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Critical Oxygen Concentration (mg/l)", ylim=c(1,4))
n <- 3
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
y0s <- c(mean(sumdat$Ccrit[which(sumdat$Test_Temp=="24")], na.rm=T), mean(sumdat$Ccrit[which(sumdat$Test_Temp=="28")], na.rm=T),mean(sumdat$Ccrit[which(sumdat$Test_Temp=="32")], na.rm=T))
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red",lwd=2)

barplot<-barplot(meansAS$Ccrit, names.arg=meansAS$Temp, ylim=c(0,3), xlab="Temperature (°C)")
arrows(y0=meansAS$Ccrit+meansAS$SE_Ccrit,y1=meansAS$Ccrit-meansAS$SE_Ccrit,
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)








#############
#import backtransformed LSMEANS (done in excel) from multivariate model done in SAS
#plot LSMEANS, decided to go with simple method with highest pts and 5% low pts
ls10<-read.csv("LSMEANS_cobia_no_dead10.csv") #eliminating dead individuals didn't change significance but plotting only alive ind looked a lot better and justification for this is MMR,SMR,and AS were lower on dead fish 

#low res figure for manuscript
par(mfrow=c(1,3))
temps<-c(24,28,32)
par(mar=c(5,4.6,4,2))
barplot<-barplot(height = ls10$EstimateC[which(ls10$VAR=="MMR")], names.arg = temps, ylim=c(0,400),ylab=expression(paste("MMR (mg kg"^"-1","h"^"-1",")")),cex.lab=1.25,cex.axis=1.25,cex.names=1.25)
arrows(y0=ls10$UpperC[which(ls10$VAR=="MMR")],y1=ls10$LowerC[which(ls10$VAR=="MMR")],
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)
text(0.6873683,251, "a",cex=1.25)
text(1.9009137,287, "a",cex=1.25)
text(3.1001821,380, "b",cex=1.25)
text(0.4,394, "A", cex=2)
text(1.9,394, "MMR", cex=2)

barplot<-barplot(height = ls10$EstimateC[which(ls10$VAR=="SMR")], names.arg = temps, ylim=c(0,200),xlab="Temperature (°C)",ylab=expression(paste("SMR (mg kg"^"-1","h"^"-1",")")),cex.lab=1.25,cex.axis=1.25,cex.names=1.25)
arrows(y0=ls10$UpperC[which(ls10$VAR=="SMR")],y1=ls10$LowerC[which(ls10$VAR=="SMR")],
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)
text(0.6873683,114, "a",cex=1.25)
text(1.9009137,147, "b",cex=1.25)
text(3.1001821,185, "c",cex=1.25)
text(0.4,197, "B", cex=2)
text(1.9,197, "SMR", cex=2)

barplot<-barplot(height = ls10$EstimateC[which(ls10$VAR=="AS")], names.arg = temps, ylim=c(0,205),ylab=expression(paste("AS (mg kg"^"-1","h"^"-1",")")),cex.lab=1.25,cex.axis=1.25,cex.names=1.25)
arrows(y0=ls10$UpperC[which(ls10$VAR=="AS")],y1=ls10$LowerC[which(ls10$VAR=="AS")],
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)
text(0.6873683,141, "a",cex=1.25)
text(1.9009137,143, "a",cex=1.25)
text(3.1001821,195, "b",cex=1.25)
text(0.4,201.5, "C", cex=2)
text(1.9,201.5, "AS", cex=2)


#convert Pcrit back to Ccrit
pcrits<-ls10[which(ls10$VAR=="Pcrit"),]
pcrits$mean_O2_pO2<-c(0.04467,0.04196,0.03994)
pcrits$mean_100O2content<-c(6.91,6.45,6.07)
pcrits$Ccrit<-pcrits$EstimateC*pcrits$mean_O2_pO2
pcrits$CcritL<-pcrits$LowerC*pcrits$mean_O2_pO2
pcrits$CcritU<-pcrits$UpperC*pcrits$mean_O2_pO2
pcrits$Scrit<-pcrits$Ccrit/pcrits$mean_100O2content*100
pcrits$ScritL<-pcrits$CcritL/pcrits$mean_100O2content*100
pcrits$ScritU<-pcrits$CcritU/pcrits$mean_100O2content*100

#plot oxygen metrics together
par(mfrow=c(1,3))
barplot<-barplot(height = ls10$EstimateC[which(ls10$VAR=="Pcrit")], names.arg = temps,ylim=c(0,70),ylab="Critical Oxygen Partial Pressure (mmHg)",cex.lab=1.25,cex.axis=1.25,cex.names=1.25)
arrows(y0=ls10$UpperC[which(ls10$VAR=="Pcrit")],y1=ls10$LowerC[which(ls10$VAR=="Pcrit")],
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)
text(0.6873683,47, "a",cex=1.25)
text(1.9009137,49, "a",cex=1.25)
text(3.1001821,68.5, "b",cex=1.25)
text(0.4,68, "A", cex=2)
text(1.9,68, expression(paste("P"[crit])), cex=2)
barplot<-barplot(height = pcrits$Ccrit, names.arg = temps,ylim=c(0,3),ylab="Critical Oxygen Concentration (mg/L)", xlab="Temperature (°C)",cex.lab=1.25,cex.axis=1.25,cex.names=1.25)
arrows(y0=pcrits$CcritU,y1=pcrits$CcritL,
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)
text(0.4,2.93, "B", cex=2)
text(1.9,2.93, expression(paste("C"[crit])), cex=2)
barplot<-barplot(height = pcrits$Scrit, names.arg = temps,ylim=c(0,50),ylab="Critical Oxygen Saturation (%)",cex.lab=1.25,cex.axis=1.25,cex.names=1.25)
arrows(y0=pcrits$ScritU,y1=pcrits$ScritL,
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)
text(0.4,49, "C", cex=2)
text(1.9,49, expression(paste("S"[crit])), cex=2)



