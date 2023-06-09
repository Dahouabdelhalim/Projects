## Correlation of all secondary traits with graphs
## Prepared by Mokhles Rahman
## March 23, 2020
## mrahman@ksu.edu
## Purpose to document five years data anlysis
rm(list=ls())
cat("\\f")
### Running univariate model yield prediction 2016
getwd()
setwd("~/Documents/BHEARD_documents/Dissertation_research/Data_Analysis")
BLUE2016=read.csv("BLUE_2016.csv")
univar.model.yield.pred2016=data.frame()
for(p in 3:ncol(BLUE2016)){
        univariate.model.2016 <- lm(GRYLD ~ BLUE2016[,p],data=BLUE2016) # writing the model
        correlation=cor(univariate.model.2016$fitted.values,BLUE2016$GRYLD)
        correlation=round(correlation,digits=2)
        traits=names(BLUE2016[p])
        corr=cbind(traits,correlation)
        univar.model.yield.pred2016=rbind(univar.model.yield.pred2016,corr)
}
### Plotting univariate model yield prediction 2016
x=univar.model.yield.pred2016
x=t(x)
x1=x[1,]
colnames(x)=c(x1[1:17],"Days to heading","Days to maturity","Plant height","Spike per sq.meter","Spike length","Spikelets per spike","Grains per spike","Thousand grain weight","Grain yield")
x=t(x)
x=data.frame(cbind(row.names(x),x),row.names=NULL)
x=x[,-2]
colnames(x)[1]=c("traits")
univar.model.yield.pred2016=x
univar.model.yield.pred2016=univar.model.yield.pred2016[-26,]
univar.model.yield.pred2016$correlation=as.numeric(as.character(univar.model.yield.pred2016$correlation))
univar.model.yield.pred2016$traits=as.character(univar.model.yield.pred2016$traits)



### Running univariate model yield prediction 2017
BLUE2017=read.csv("BLUE_2017.csv") 
univar.model.yield.pred2017=data.frame()
for(p in 3:ncol(BLUE2017)){
        univariate.model.2017 <- lm(GRYLD ~ BLUE2017[,p],data=BLUE2017) # writing the model
        correlation=cor(univariate.model.2017$fitted.values,BLUE2017$GRYLD)
        correlation=round(correlation,digits=2)
        traits=names(BLUE2017[p])
        corr=cbind(traits,correlation)
        univar.model.yield.pred2017=rbind(univar.model.yield.pred2017,corr)
}
### Plotting univariate model yield prediction 2017
x=univar.model.yield.pred2017
x=t(x)
x1=x[1,]
colnames(x)=c(x1[1:28],"Days to heading","Days to maturity","Plant height","Spike per sq.meter","Spike length","Spikelets per spike","Grains per spike","Thousand grain weight","Grain yield")
x=t(x)
x=data.frame(cbind(row.names(x),x),row.names=NULL)
x=x[,-2]
colnames(x)[1]=c("traits")
univar.model.yield.pred2017=x
univar.model.yield.pred2017=univar.model.yield.pred2017[-37,]
univar.model.yield.pred2017$correlation=as.numeric(as.character(univar.model.yield.pred2017$correlation))
univar.model.yield.pred2017$traits=as.character(univar.model.yield.pred2017$traits)

### Running univariate model yield prediction 2018
BLUE2018=read.csv("BLUE_2018.csv")
univar.model.yield.pred2018=data.frame()
for(p in 3:ncol(BLUE2018)){
        univariate.model.2018 <- lm(GRYLD ~ BLUE2018[,p],data=BLUE2018) # writing the model
        correlation=cor(univariate.model.2018$fitted.values,BLUE2018$GRYLD)
        correlation=round(correlation,digits=2)
        traits=names(BLUE2018[p])
        corr=cbind(traits,correlation)
        univar.model.yield.pred2018=rbind(univar.model.yield.pred2018,corr)
}
### Plotting univariate model yield prediction 2018
x=univar.model.yield.pred2018
x=t(x)
x1=x[1,]
colnames(x)=c(x1[1:24],"Days to heading","Days to maturity","Plant height","Spike per sq.meter","Spike length","Spikelets per spike","Grains per spike","Thousand grain weight","Grain yield")
x=t(x)
x=data.frame(cbind(row.names(x),x),row.names=NULL)
x=x[,-2]
colnames(x)[1]=c("traits")
univar.model.yield.pred2018=x
univar.model.yield.pred2018=univar.model.yield.pred2018[-33,]
univar.model.yield.pred2018$correlation=as.numeric(as.character(univar.model.yield.pred2018$correlation))
univar.model.yield.pred2018$traits=as.character(univar.model.yield.pred2018$traits)

### Running univariate model yield prediction 2019
BLUE2019=read.csv("BLUE_2019.csv")
univar.model.yield.pred2019=data.frame()
for(p in 3:ncol(BLUE2019)){
        univariate.model.2019 <- lm(GRYLD ~ BLUE2019[,p],data=BLUE2019) # writing the model
        correlation=cor(univariate.model.2019$fitted.values,BLUE2019$GRYLD)
        correlation=round(correlation,digits=2)
        traits=names(BLUE2019[p])
        corr=cbind(traits,correlation)
        univar.model.yield.pred2019=rbind(univar.model.yield.pred2019,corr)
}
### Plotting univariate model yield prediction 2019
x=univar.model.yield.pred2019
x=t(x)
x1=x[1,]
colnames(x)=c(x1[1:26],"Days to heading","Days to maturity","Plant height","Spike per sq.meter","Spike length","Spikelets per spike","Grains per spike","Thousand grain weight","Grain yield")
x=t(x)
x=data.frame(cbind(row.names(x),x),row.names=NULL)
x=x[,-2]
colnames(x)[1]=c("traits")
univar.model.yield.pred2019=x
univar.model.yield.pred2019=univar.model.yield.pred2019[-35,]
univar.model.yield.pred2019$correlation=as.numeric(as.character(univar.model.yield.pred2019$correlation))
univar.model.yield.pred2019$traits=as.character(univar.model.yield.pred2019$traits)

### Running univariate model yield prediction 
BLUE2020=read.csv("BLUE_2020.csv")
BLUE2020=BLUE2020[,-c(33:36)]
univar.model.yield.pred2020=data.frame()
for(p in 3:ncol(BLUE2020)){
        univariate.model.2020 <- lm(GRYLD ~ BLUE2020[,p],data=BLUE2020) # writing the model
        correlation=cor(univariate.model.2020$fitted.values,BLUE2020$GRYLD)
        correlation=round(correlation,digits=2)
        traits=names(BLUE2020[p])
        corr=cbind(traits,correlation)
        univar.model.yield.pred2020=rbind(univar.model.yield.pred2020,corr)
}
### Plotting univariate model yield prediction 2020
x=univar.model.yield.pred2020
x=t(x)
x1=x[1,]
colnames(x)=c(x1[1:30],"Days to heading","Days to maturity","Plant height","Spike per sq.meter","Spikelets per spike","Grains per spike","Thousand grain weight","Grain yield")
x=t(x)
x=data.frame(cbind(row.names(x),x),row.names=NULL)
x=x[,-2]
colnames(x)[1]=c("traits")
univar.model.yield.pred2020=x
univar.model.yield.pred2020=univar.model.yield.pred2020[-38,]
univar.model.yield.pred2020$correlation=as.numeric(as.character(univar.model.yield.pred2020$correlation))
univar.model.yield.pred2020$traits=as.character(univar.model.yield.pred2020$traits)

dev.off()
### Plotting univariate model yield prediciton onto panel plot
# pdf(file = 'Figure 2.pdf', height = 18, width = 12)
par(mfrow=c(3,2))
par(mgp=c(3,0.5,0))
par(mar=c(7.7,2,.75,0))
barplot(univar.model.yield.pred2016$correlation,names.arg = univar.model.yield.pred2016$traits,las=2,cex.axis=0.7,cex.names=0.7,ylim=c(0,0.8),col=c(rep('red',8), rep('blue',9), rep('green',8)))
# abline(h=0.5,lty=3)
title("A",adj=0)
legend("top",c("CT","NDVI","Agronomic"),fill=c("red","blue","green"),bty="n",ncol=3)
barplot(univar.model.yield.pred2017$correlation,names.arg = univar.model.yield.pred2017$traits,las=2,cex.axis=0.7,cex.names=0.7,ylim=c(0,0.8),col=c(rep('red',14), rep('blue',14), rep('green',8)))
# abline(h=0.5,lty=3)
title("B",adj=0)
legend("top",c("CT","NDVI","Agronomic"),fill=c("red","blue","green"),bty="n",ncol=3)
barplot(univar.model.yield.pred2018$correlation,names.arg = univar.model.yield.pred2018$traits,las=2,cex.axis=0.7,cex.names=0.7,ylim=c(0,0.8),col=c(rep('red',12), rep('blue',12), rep('green',8)))
# abline(h=0.5,lty=3)
title("C",adj=0)
legend("top",c("CT","NDVI","Agronomic"),fill=c("red","blue","green"),bty="n",ncol=3)
barplot(univar.model.yield.pred2019$correlation,names.arg = univar.model.yield.pred2019$traits,las=2,cex.axis=0.7,cex.names=0.7,ylim=c(0,0.8),col=c(rep('red',13), rep('blue',13), rep('green',8)))
# abline(h=0.5,lty=3)
title("D",adj=0)
legend("top",c("CT","NDVI","Agronomic"),fill=c("red","blue","green"),bty="n",ncol=3)
barplot(univar.model.yield.pred2020$correlation,names.arg = univar.model.yield.pred2020$traits,las=2,cex.axis=0.7,cex.names=0.7,ylim=c(0,0.8),col=c(rep('red',15), rep('blue',15), rep('green',11)))
# abline(h=0.5,lty=3)
title("E",adj=0)
legend("top",c("CT","NDVI","Agronomic"),fill=c("red","blue","green"),bty="n",ncol=3)




