## Correlation of all secondary traits with graphs
## Prepared by Mokhles Rahman
## March 23, 2020
## mrahman@ksu.edu
## Purpose to document five years data anlysis
rm(list=ls())
cat("\\f")
### Running correlations between grain yield and the traits 2016
### Estimation of correlation considering all trials together 2016
setwd("~/Documents/BHEARD_documents/Dissertation_research/Data_Analysis")
BLUE2016=read.csv("BLUE_2016.csv")
cor.for.figure2016=data.frame()
cor.result2016=data.frame()
for(p in 3:ncol(BLUE2016)){
        correlation=cor(BLUE2016[,p],BLUE2016$GRYLD)
        correlation=round(correlation,digits=2)
        trait=names(BLUE2016[p])
        corr=cbind(trait,correlation)
        cor.for.figure2016=rbind(cor.for.figure2016,corr)
        aa=cor.test(BLUE2016[,p],BLUE2016$GRYLD)
        Sig=aa[["p.value"]]
        Sig=if(Sig>0.05){print("ns")} else if(Sig>=0.01){print("*")} else if(Sig>=0.001){print("**")} else{print("***")}
        r=paste(correlation,Sig)
        traits=names(BLUE2016[p])
        corr=cbind(traits,r)
        cor.result2016=rbind(cor.result2016,corr)
}
# write.csv(cor.result2016,file="Correlation_2016.csv",row.names=FALSE,quote=FALSE)
# write.csv(cor.for.figure,file="Correlation for figure 2016",row.names=FALSE,quote=FALSE)
dev.off()
par(mar=c(3,3,1,1))
par(mgp=c(1.5,0.5,0))
par(mfrow=c(3,2))

aa=cor.for.figure2016[1:8,]
aa=as.matrix(aa)
aa=aa[,2]
aa=as.numeric(aa)
DAS=c("59","68","78","90","95","97","104","110")
DAS=as.numeric(DAS)
plot(DAS,aa,ylim=c(-0.75,0.75),xlim=c(35,120),ylab="Correlation",xlab="Days after sowing",type="p",col="red",at=c("40","50","60","70","80","90","100","110","120"))
lines(lowess(DAS,aa,f=2/3),col="red")
abline(h=0,lty=3)
abline(v=70,lty=5)
abline(v=104,lty=3)
title("A",adj=0)
par(new=TRUE)

ab=cor.for.figure2016[9:17,]
ab=as.matrix(ab)
ab=ab[,2]
ab=as.numeric(ab)
DAS=c("57","61","66","74","90","95","98","105","110")
DAS=as.numeric(DAS)
plot(DAS,ab,ylim=c(-0.75,0.75),xlim=c(35,120),ylab="Correlation",xlab="Days after sowing",type="p",col="blue")
lines(lowess(DAS,ab,f=2/3),col="blue")
legend("topleft",c("CT","NDVI"),fill=c("red","blue"),bty="n",ncol=2)


### Estimation of correlation by trial basis 2016
cor.bytrial2016=data.frame()
cor.bytrial.forfigure2016=data.frame()
for(trial in 1:10){
        per.trial.blue=BLUE2016[as.character(BLUE2016$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        for(q in 3:ncol(per.trial.blue)){
                correlation=cor(per.trial.blue[,q],per.trial.blue$GRYLD)
                correlation=round(correlation,digits=2)
                trait=names(BLUE2016[q])
                Trial=trial
                corr=cbind(Trial,trait,correlation)
                cor.bytrial.forfigure2016=rbind(cor.bytrial.forfigure2016,corr)
                aa=cor.test(per.trial.blue[,q],per.trial.blue$GRYLD)
                Sig=aa[["p.value"]]
                Sig=if(Sig>0.05){print("ns")} else if(Sig>=0.01){print("*")} else if(Sig>=0.001){print("**")} else{print("***")}
                r=paste(correlation,Sig)
                trait=names(per.trial.blue[q])
                Trial=trial
                correl=cbind(Trial,trait,r)
                cor.bytrial2016=rbind(cor.bytrial2016,correl)
        }
}
# write.csv(cor.bytrial2016,file="Correlation by trial 2016.csv",row.names = FALSE,quote = FALSE)
# write.csv(cor.bytrial.forfigure2016,file="Correlation for figure by trial 2016.csv",row.names = FALSE,quote = FALSE)
cor2016=cor.bytrial2016$r
cor2016=matrix(cor2016,nrow = 26,ncol = 10)
Traits=colnames(BLUE2016)[3:ncol(BLUE2016)]
cor2016=data.frame(Traits,cor2016)
colnames(cor2016)=c("Traits","Trial_1","Trial_2","Trial_3","Trial_4","Trial_5","Trial_6","Trial_7","Trial_8","Trial_9","Trial_10")
cor2016=cor2016[-26,]
write.csv(cor2016,file="Correlation_byTrial_2016.csv",row.names=FALSE,quote=FALSE)

### Running correlations between grain yield and the traits 2017
### Estimation of correlation considering all trials together 2017
BLUE2017=read.csv("BLUE_2017.csv")
cor.for.figure2017=data.frame()
cor.result2017=data.frame()
for(p in 3:ncol(BLUE2017)){
        correlation=cor(BLUE2017[,p],BLUE2017$GRYLD)
        correlation=round(correlation,digits=2)
        trait=names(BLUE2017[p])
        corr=cbind(trait,correlation)
        cor.for.figure2017=rbind(cor.for.figure2017,corr)
        aa=cor.test(BLUE2017[,p],BLUE2017$GRYLD)
        Sig=aa[["p.value"]]
        Sig=if(Sig>0.05){print("ns")} else if(Sig>=0.01){print("*")} else if(Sig>=0.001){print("**")} else{print("***")}
        r=paste(correlation,Sig)
        traits=names(BLUE2017[p])
        corr=cbind(traits,r)
        cor.result2017=rbind(cor.result2017,corr)
}
# write.csv(cor.result2017,file="Correlation_2017.csv",row.names=FALSE,quote=FALSE)
# write.csv(cor.for.figure,file="Correlation for figure 2017",row.names=FALSE,quote=FALSE)

aa=cor.for.figure2017[1:14,]
aa=as.matrix(aa)
aa=aa[,2]
aa=as.numeric(aa)
DAS=c("38","43","48","54","59","65","70","75","80","86","90","95","100","106")
DAS=as.numeric(DAS)
plot(DAS,aa,ylim=c(-0.75,0.75),xlim=c(35,120),ylab="Correlation",xlab="Days after sowing",type="p",col="red",at=c("40","50","60","70","80","90","100","110","120"))
lines(lowess(DAS,aa,f=2/3),col="red")
abline(h=0,lty=3)
abline(v=70,lty=5)
abline(v=105,lty=3)
title("B",adj=0)
par(new=TRUE)

ab=cor.for.figure2017[15:28,]
ab=as.matrix(ab)
ab=ab[,2]
ab=as.numeric(ab)
DAS=c("37","42","48","54","59","65","70","75","80","85","90","95","100","106")
DAS=as.numeric(DAS)
plot(DAS,ab,ylim=c(-0.75,0.75),xlim=c(35,120),ylab="Correlation",xlab="Days after sowing",type="p",col="blue")
lines(lowess(DAS,ab,f=2/3),col="blue")
legend("topleft",c("CT","NDVI"),fill=c("red","blue"),bty="n",ncol=2)

### Estimation of correlation by trial basis 2017
cor.bytrial2017=data.frame()
cor.bytrial.forfigure2017=data.frame()
for(trial in 1:11){
        per.trial.blue=BLUE2017[as.character(BLUE2017$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        for(q in 3:ncol(per.trial.blue)){
                correlation=cor(per.trial.blue[,q],per.trial.blue$GRYLD)
                correlation=round(correlation,digits=2)
                trait=names(BLUE2017[q])
                Trial=trial
                corr=cbind(Trial,trait,correlation)
                cor.bytrial.forfigure2017=rbind(cor.bytrial.forfigure2017,corr)
                aa=cor.test(per.trial.blue[,q],per.trial.blue$GRYLD)
                Sig=aa[["p.value"]]
                Sig=if(Sig>0.05){print("ns")} else if(Sig>=0.01){print("*")} else if(Sig>=0.001){print("**")} else{print("***")}
                r=paste(correlation,Sig)
                trait=names(per.trial.blue[q])
                Trial=trial
                correl=cbind(Trial,trait,r)
                cor.bytrial2017=rbind(cor.bytrial2017,correl)
        }
}
# write.csv(cor.bytrial2017,file="Correlation by trial 2017.csv",row.names = FALSE,quote = FALSE)
# write.csv(cor.bytrial.forfigure2017,file="Correlation for figure by trial 2017.csv",row.names = FALSE,quote = FALSE)
cor2017=cor.bytrial2017$r
cor2017=matrix(cor2017,nrow = 37,ncol = 11)
Traits=colnames(BLUE2017)[3:ncol(BLUE2017)]
cor2017=data.frame(Traits,cor2017)
colnames(cor2017)=c("Traits","Trial_1","Trial_2","Trial_3","Trial_4","Trial_5","Trial_6","Trial_7","Trial_8","Trial_9","Trial_10","Trial_11")
cor2017=cor2017[-37,]
write.csv(cor2017,file="Correlation_byTrial_2017.csv",row.names=FALSE,quote=FALSE)



### Running correlations between grain yield and the traits 2018
### Estimation of correlation considering all trials together 2018
BLUE2018=read.csv("BLUE_2018.csv")
cor.for.figure2018=data.frame()
cor.result2018=data.frame()
for(p in 3:ncol(BLUE2018)){
        correlation=cor(BLUE2018[,p],BLUE2018$GRYLD)
        correlation=round(correlation,digits=2)
        trait=names(BLUE2018[p])
        corr=cbind(trait,correlation)
        cor.for.figure2018=rbind(cor.for.figure2018,corr)
        aa=cor.test(BLUE2018[,p],BLUE2018$GRYLD)
        Sig=aa[["p.value"]]
        Sig=if(Sig>0.05){print("ns")} else if(Sig>=0.01){print("*")} else if(Sig>=0.001){print("**")} else{print("***")}
        r=paste(correlation,Sig)
        traits=names(BLUE2018[p])
        corr=cbind(traits,r)
        cor.result2018=rbind(cor.result2018,corr)
}
# write.csv(cor.result2018,file="Correlation_2018.csv",row.names=FALSE,quote=FALSE)
# write.csv(cor.for.figure,file="Correlation for figure 2018",row.names=FALSE,quote=FALSE)
aa=cor.for.figure2018[1:12,]
aa=as.matrix(aa)
aa=aa[,2]
aa=as.numeric(aa)
DAS=c("58","63","68","73","77","82","88","92","96","101","106","111")
DAS=as.numeric(DAS)
plot(DAS,aa,ylim=c(-0.75,0.75),xlim=c(35,120),ylab="Correlation",xlab="Days after sowing",type="p",col="red",at=c("40","50","60","70","80","90","100","110","120"))
lines(lowess(DAS,aa,f=2/3),col="red")
abline(h = 0, lty=3)
abline(v=72,lty=5)
abline(v=105,lty=3)
title("C",adj=0)
par(new=TRUE)


ab=cor.for.figure2018[13:24,]
ab=as.matrix(ab)
ab=ab[,2]
ab=as.numeric(ab)
DAS=c("58","63","68","72","77","83","88","92","96","101","106","111")
DAS=as.numeric(DAS)
plot(DAS,ab,ylim=c(-0.75,0.75),xlim=c(35,120),ylab="Correlation",xlab="Days after sowing",type="p",col="blue")
lines(lowess(DAS,ab,f=2/3),col="blue")
legend("topleft",c("CT","NDVI"),fill=c("red","blue"),bty="n",ncol=2)

### Estimation of correlation by trial basis 2018
cor.bytrial2018=data.frame()
cor.bytrial.forfigure2018=data.frame()
for(trial in 1:11){
        per.trial.blue=BLUE2018[as.character(BLUE2018$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        for(q in 3:ncol(per.trial.blue)){
                correlation=cor(per.trial.blue[,q],per.trial.blue$GRYLD)
                correlation=round(correlation,digits=2)
                trait=names(BLUE2018[q])
                Trial=trial
                corr=cbind(Trial,trait,correlation)
                cor.bytrial.forfigure2018=rbind(cor.bytrial.forfigure2018,corr)
                aa=cor.test(per.trial.blue[,q],per.trial.blue$GRYLD)
                Sig=aa[["p.value"]]
                Sig=if(Sig>0.05){print("ns")} else if(Sig>=0.01){print("*")} else if(Sig>=0.001){print("**")} else{print("***")}
                r=paste(correlation,Sig)
                trait=names(per.trial.blue[q])
                Trial=trial
                correl=cbind(Trial,trait,r)
                cor.bytrial2018=rbind(cor.bytrial2018,correl)
        }
}
# write.csv(cor.bytrial2018,file="Correlation by trial 2018.csv",row.names = FALSE,quote = FALSE)
# write.csv(cor.bytrial.forfigure2018,file="Correlation for figure by trial 2018.csv",row.names = FALSE,quote = FALSE)
cor2018=cor.bytrial2018$r
cor2018=matrix(cor2018,nrow = 33,ncol = 11)
Traits=colnames(BLUE2018)[3:ncol(BLUE2018)]
cor2018=data.frame(Traits,cor2018)
colnames(cor2018)=c("Traits","Trial_1","Trial_2","Trial_3","Trial_4","Trial_5","Trial_6","Trial_7","Trial_8","Trial_9","Trial_10","Trial_11")
cor2018=cor2018[-33,]
write.csv(cor2018,file="Correlation_byTrial_2018.csv",row.names=FALSE,quote=FALSE)



### Running correlations between grain yield and the traits 2019
### Estimation of correlation considering all trials together 2019
BLUE2019=read.csv("BLUE_2019.csv")
cor.for.figure2019=data.frame()
cor.result2019=data.frame()
for(p in 3:ncol(BLUE2019)){
        correlation=cor(BLUE2019[,p],BLUE2019$GRYLD)
        correlation=round(correlation,digits=2)
        trait=names(BLUE2019[p])
        corr=cbind(trait,correlation)
        cor.for.figure2019=rbind(cor.for.figure2019,corr)
        aa=cor.test(BLUE2019[,p],BLUE2019$GRYLD)
        Sig=aa[["p.value"]]
        Sig=if(Sig>0.05){print("ns")} else if(Sig>=0.01){print("*")} else if(Sig>=0.001){print("**")} else{print("***")}
        r=paste(correlation,Sig)
        traits=names(BLUE2019[p])
        corr=cbind(traits,r)
        cor.result2019=rbind(cor.result2019,corr)
}
# write.csv(cor.result2019,file="Correlation_2019.csv",row.names=FALSE,quote=FALSE)
# write.csv(cor.for.figure,file="Correlation for figure 2019",row.names=FALSE,quote=FALSE)

aa=cor.for.figure2019[1:13,]
aa=as.matrix(aa)
aa=aa[,2]
aa=as.numeric(aa)
DAS=c("56","60","64","69","75","82","87","93","97","103","108","112","117")
DAS=as.numeric(DAS)
plot(DAS,aa,ylim=c(-0.75,0.75),xlim=c(35,120),ylab="Correlation",xlab="Days after sowing",type="p",col="red",at=c("40","50","60","70","80","90","100","110","120"))
lines(lowess(DAS,aa,f=2/3),col="red")
abline(h = 0, lty=3)
abline(v=77,lty=5)
abline(v=110,lty=3)
title("D",adj=0)
par(new=TRUE)

ab=cor.for.figure2019[14:26,]
ab=as.matrix(ab)
ab=ab[,2]
ab=as.numeric(ab)
DAS=c("54","60","64","69","75","82","86","92","97","103","107","112","117")
DAS=as.numeric(DAS)
plot(DAS,ab,ylim=c(-0.75,0.75),xlim=c(35,120),ylab="Correlation",xlab="Days after sowing",type="p",col="blue")
lines(lowess(DAS,ab,f=2/3),col="blue")
legend("topleft",c("CT","NDVI"),fill=c("red","blue"),bty="n",ncol=2)

### Estimation of correlation by trial basis 2019
cor.bytrial2019=data.frame()
cor.bytrial.forfigure2019=data.frame()
for(trial in 1:10){
        per.trial.blue=BLUE2019[as.character(BLUE2019$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        for(q in 3:ncol(per.trial.blue)){
                correlation=cor(per.trial.blue[,q],per.trial.blue$GRYLD)
                correlation=round(correlation,digits=2)
                trait=names(BLUE2019[q])
                Trial=trial
                corr=cbind(Trial,trait,correlation)
                cor.bytrial.forfigure2019=rbind(cor.bytrial.forfigure2019,corr)
                aa=cor.test(per.trial.blue[,q],per.trial.blue$GRYLD)
                Sig=aa[["p.value"]]
                Sig=if(Sig>0.05){print("ns")} else if(Sig>=0.01){print("*")} else if(Sig>=0.001){print("**")} else{print("***")}
                r=paste(correlation,Sig)
                trait=names(per.trial.blue[q])
                Trial=trial
                correl=cbind(Trial,trait,r)
                cor.bytrial2019=rbind(cor.bytrial2019,correl)
        }
}
# write.csv(cor.bytrial2019,file="Correlation by trial 2019.csv",row.names = FALSE,quote = FALSE)
# write.csv(cor.bytrial.forfigure2019,file="Correlation for figure by trial 2019.csv",row.names = FALSE,quote = FALSE)
cor2019=cor.bytrial2019$r
cor2019=matrix(cor2019,nrow = 35,ncol = 10)
Traits=colnames(BLUE2019)[3:ncol(BLUE2019)]
cor2019=data.frame(Traits,cor2019)
colnames(cor2019)=c("Traits","Trial_1","Trial_2","Trial_3","Trial_4","Trial_5","Trial_6","Trial_7","Trial_8","Trial_9","Trial_10")
cor2019=cor2019[-35,]
write.csv(cor2019,file="Correlaiton_byTrial_2019.csv",row.names=FALSE,quote=FALSE)


### Running correlations between grain yield and the traits 2010
### Estimation of correlation considering all trials together 2020
BLUE2020=read.csv("BLUE_2020.csv")
cor.for.figure2020=data.frame()
cor.result2020=data.frame()
for(p in 3:ncol(BLUE2020)){
        correlation=cor(BLUE2020[,p],BLUE2020$GRYLD)
        correlation=round(correlation,digits=2)
        trait=names(BLUE2020[p])
        corr=cbind(trait,correlation)
        cor.for.figure2020=rbind(cor.for.figure2020,corr)
        aa=cor.test(BLUE2020[,p],BLUE2020$GRYLD)
        Sig=aa[["p.value"]]
        Sig=if(Sig>0.05){print("ns")} else if(Sig>=0.01){print("*")} else if(Sig>=0.001){print("**")} else{print("***")}
        r=paste(correlation,Sig)
        traits=names(BLUE2020[p])
        corr=cbind(traits,r)
        cor.result2020=rbind(cor.result2020,corr)
}
# write.csv(cor.result2020,file="Correlation_2020.csv",row.names=FALSE,quote=FALSE)
# write.csv(cor.for.figure,file="Correlation for figure 2020",row.names=FALSE,quote=FALSE)

aa=cor.for.figure2020[1:15,]
aa=as.matrix(aa)
aa=aa[,2]
aa=as.numeric(aa)
DAS=c("42","46","51","56","60","65","70","75","80","86","91","97","102","107","112")
DAS=as.numeric(DAS)
plot(DAS,aa,ylim=c(-0.75,0.75),xlim=c(35,120),ylab="Correlation",xlab="Days after sowing",type="p",col="red",at=c("40","50","60","70","80","90","100","110","120"))
lines(lowess(DAS,aa,f=2/3),col="red")
abline(h = 0, lty=3)
abline(v=73,lty=5)
abline(v=104,lty=3)
title("E",adj=0)
par(new=TRUE)


ab=cor.for.figure2020[16:30,]
ab=as.matrix(ab)
ab=ab[,2]
ab=as.numeric(ab)
DAS=c("42","46","51","56","60","65","70","75","80","86","91","97","102","107","112")
DAS=as.numeric(DAS)
plot(DAS,ab,ylim=c(-0.75,0.75),xlim=c(35,120),ylab="Correlation",xlab="Days after sowing",type="p",col="blue")
lines(lowess(DAS,ab,f=2/3),col="blue")
legend("topleft",c("CT","NDVI"),fill=c("red","blue"),bty="n",ncol=2)


### Estimation of correlation by trial basis 2020
cor.bytrial2020=data.frame()
cor.bytrial.forfigure2020=data.frame()
for(trial in 1:11){
        per.trial.blue=BLUE2020[as.character(BLUE2020$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        for(q in 3:ncol(per.trial.blue)){
                correlation=cor(per.trial.blue[,q],per.trial.blue$GRYLD)
                correlation=round(correlation,digits=2)
                trait=names(BLUE2020[q])
                Trial=trial
                corr=cbind(Trial,trait,correlation)
                cor.bytrial.forfigure2020=rbind(cor.bytrial.forfigure2020,corr)
                aa=cor.test(per.trial.blue[,q],per.trial.blue$GRYLD)
                Sig=aa[["p.value"]]
                Sig=if(Sig>0.05){print("ns")} else if(Sig>=0.01){print("*")} else if(Sig>=0.001){print("**")} else{print("***")}
                r=paste(correlation,Sig)
                trait=names(per.trial.blue[q])
                Trial=trial
                correl=cbind(Trial,trait,r)
                cor.bytrial2020=rbind(cor.bytrial2020,correl)
        }
}
# write.csv(cor.bytrial2020,file="Correlation by trial 2020.csv",row.names = FALSE,quote = FALSE)
# write.csv(cor.bytrial.forfigure2020,file="Correlation for figure by trial 2020.csv",row.names = FALSE,quote = FALSE)
cor2020=cor.bytrial2020$r
cor2020=matrix(cor2020,nrow = 42,ncol = 11)
Traits=colnames(BLUE2020)[3:ncol(BLUE2020)]
cor2020=data.frame(Traits,cor2020)
colnames(cor2020)=c("Traits","Trial_1","Trial_2","Trial_3","Trial_4","Trial_5","Trial_6","Trial_7","Trial_8","Trial_9","Trial_10","Trial_11")
cor2020=cor2020[-42,]
write.csv(cor2020,file="Correlation_byTrial_2020.csv",row.names=FALSE,quote=FALSE)



