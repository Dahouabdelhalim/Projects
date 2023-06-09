#An example of overfitting and cross validation using a subset of compiled data in Saenger and Evans 2019
#Note that formMgCa_env.dat.R should be run prior to running this script

setwd("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data");


source("/Volumes/GoogleDrive/My Drive/R_documents/funcs/cal.val.R")
source("/Volumes/GoogleDrive/My Drive/R_documents/funcs/cal.val.nls.R")

library(ncdf4)
library(reshape2)
library(plyr)


#coretop foram MgCa data and compiled environmental data from the closest gridbox
dat.all<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/Gbulloides.MgCa.env.dat.csv")

#empirically account for reductive cleaning following Rosenthal et al., 2004 relationship
red<-which(dat.all$reductive=="Y")
dat.all$Mg.Ca[red]<-(dat.all$Mg.Ca[red]-0.2)/0.85

dat1<-dat.all[which(dat.all$Reference=="Elderfield and Ganssen 2000" | dat.all$Reference=="Yu et al. 2008"),]
MgCa1<-dat1$Mg.Ca
T1<-dat1[,20:26]
T1.sd<-dat1[,62:68]
S1<-dat1[,13:19]
S1.sd<-dat1[,55:61]
Omega.s1<-dat1[,27:33]
Omega.s1.sd<-dat1[,69:75]

dat2<-dat.all[which(dat.all$Reference=="Cleroux et al. 2008"),]
MgCa2<-dat2$Mg.Ca
T2<-dat2[,20:26]
T2.sd<-dat2[,62:68]
S2<-dat2[,13:19]
S2.sd<-dat2[,55:61]
Omega.s2<-dat2[,27:33]
Omega.s2.sd<-dat2[,69:75]

MgCa12<-c(MgCa1,MgCa2)
T12<-rbind(T1,T2)
T12.sd<-rbind(T1.sd,T2.sd)
S12<-rbind(S1,S2)
S12.sd<-rbind(S1.sd,S2.sd)
Omega.s12<-rbind(Omega.s1,Omega.s2)
Omega.s12.sd<-rbind(Omega.s1.sd,Omega.s2.sd)

dat3<-dat.all[which(dat.all$Reference=="Quintana Krupinski 2017"),]
MgCa3<-dat3$Mg.Ca
T3<-dat3[,20:26]
T3.sd<-dat2[,62:68]
S3<-dat3[,13:19]
S3.sd<-dat3[,55:61]
Omega.s3<-dat3[,27:33]
Omega.s3.sd<-dat3[,69:75]

#***********************************************#
n.it<-1000
r2_1<-matrix(,nrow=n.it,ncol=3)
colnames(r2_1)<-c("T only","T & S", "T, S & Î©")
rmse_1<-matrix(,nrow=n.it,ncol=3)
colnames(rmse_1)<-colnames(r2_1)
rmse_2<-matrix(,nrow=n.it,ncol=3)
colnames(rmse_2)<-colnames(r2_1)
ce_2<-matrix(,nrow=n.it,ncol=3)
colnames(ce_2)<-colnames(r2_1)
r2_12<-rep(0,n.it)
rmse_12<-rep(0,n.it)
rmse_3<-rep(0,n.it)
ce_3<-rep(0,n.it)

#assume analytical uncertainty on Mg/Ca of 1.5% RSD
rsd<-.015 

#T-only calibration
MgCac.pred<-matrix(,nrow=length(MgCa1),ncol=n.it)
MgCav.pred<-matrix(,nrow=length(MgCa2),ncol=n.it)

for (i in 1:n.it) {
	#calibrate with record 1
	rand<-sample(1:5,1)
	Tc<-T1[,rand]+rnorm(nrow(T1),0,T1.sd[,rand])
	MgCac<-MgCa1+rnorm(1,0,MgCa1*rsd)
	fit<-lm(log(MgCac)~Tc)
	r2_1[i,1]<-summary(fit)$r.squared													#r2 of calibration
	rmse_1[i,1]<-sqrt(sum((exp(predict(fit))-MgCac)^2)/length(MgCac))					#rmse between predicted and measured calibration MgCa
	MgCac.pred[,i]<-exp(predict(fit))
	
	#validate with record 2
	rand<-sample(1:5,1)
	v<-data.frame(Tc=T2[,rand]+rnorm(nrow(T2),0,T2.sd[,rand]))
	MgCav<-MgCa2+rnorm(1,0,MgCa2*rsd)
	rmse_2[i,1]<-sqrt(sum((exp(predict(fit,v))-MgCav)^2)/nrow(v))					#rmse between MgCa predicted from the calibration for T2 and validation MgCa
	ce_2[i,1]<-1-(sum((exp(predict(fit,v))-MgCav)^2)/sum((mean(MgCav)-MgCav)^2))
	MgCav.pred[,i]<-exp(predict(fit,v))
	}

#T and S calibration
for (i in 1:n.it) {
	#calibrate with record 1
	rand<-sample(1:5,1)
	Tc<-T1[,rand]+rnorm(nrow(T1),0,T1.sd[,rand])
	Sc<-S1[,rand]+rnorm(nrow(S1),0,S1.sd[,rand])
	MgCac<-MgCa1+rnorm(1,0,MgCa1*rsd)
	fit<-lm(log(MgCac)~Tc+Sc)
	r2_1[i,2]<-summary(fit)$r.squared													#r2 of calibration
	rmse_1[i,2]<-sqrt(sum((exp(predict(fit))-MgCac)^2)/length(MgCac))					#rmse between predicted and measured calibration MgCa
	MgCac.pred[,i]<-exp(predict(fit))
	
	#validate with record 2
	rand<-sample(1:5,1)
	v<-data.frame(Tc=T2[,rand]+rnorm(nrow(T2),0,T2.sd[,rand]),Sc=S2[,rand]+rnorm(nrow(S2),0,S2.sd[,rand]))
	MgCav<-MgCa2+rnorm(1,0,MgCa2*rsd)
	rmse_2[i,2]<-sqrt(sum((exp(predict(fit,v))-MgCav)^2)/nrow(v))					#rmse between MgCa predicted from the calibration for T2 and validation MgCa
	ce_2[i,2]<-1-(sum((exp(predict(fit,v))-MgCav)^2)/sum((mean(MgCav)-MgCav)^2))
	MgCav.pred[,i]<-exp(predict(fit,v))
	
	#combined calibration with 1 and 2
	rand<-sample(1:5,1)
	Tc<-T12[,rand]+rnorm(nrow(T12),0,T12.sd[,rand])
	Sc<-S12[,rand]+rnorm(nrow(S12),0,S12.sd[,rand])
	MgCac<-MgCa12+rnorm(1,0,MgCa1*rsd)
	fit<-lm(log(MgCac)~Tc+Sc)
	r2_12[i]<-summary(fit)$r.squared													#r2 of calibration
	rmse_12[i]<-sqrt(sum((exp(predict(fit))-MgCac)^2)/length(MgCac))					#rmse between predicted and measured calibration MgCa
	
	#validate with record 3
	rand<-sample(1:5,1)
	v<-data.frame(Tc=T3[,rand]+rnorm(nrow(T3),0,T3.sd[,rand]),Sc=S3[,rand]+rnorm(nrow(S3),0,S3.sd[,rand]))
	MgCav<-MgCa3+rnorm(1,0,MgCa3*rsd)
	rmse_3[i]<-sqrt(sum((exp(predict(fit,v))-MgCav)^2)/nrow(v))					#rmse between MgCa predicted from the calibration for T2 and validation MgCa
	ce_3[i]<-1-(sum((exp(predict(fit,v))-MgCav)^2)/sum((mean(MgCav)-MgCav)^2))
	}

#T, S and surface omega calibration
for (i in 1:n.it) {
	#calibrate with record 1
	rand<-sample(1:5,1)
	Tc<-T1[,rand]+rnorm(nrow(T1),0,T1.sd[,rand])
	Sc<-S1[,rand]+rnorm(nrow(S1),0,S1.sd[,rand])
	Omega.sc<-Omega.s1[,rand]+rnorm(nrow(Omega.s1),0,Omega.s1.sd[,rand])
	MgCac<-MgCa1+rnorm(1,0,MgCa1*rsd)
	fit<-lm(log(MgCac)~Tc+Sc+log(Omega.sc))
	r2_1[i,3]<-summary(fit)$r.squared													#r2 of calibration
	rmse_1[i,3]<-sqrt(sum((exp(predict(fit))-MgCac)^2)/length(MgCac))					#rmse between predicted and measured calibration MgCa
	MgCac.pred[,i]<-exp(predict(fit))
	
	#validate with record 2
	rand<-sample(1:5,1)
	v<-data.frame(Tc=T2[,rand]+rnorm(nrow(T2),0,T2.sd[,rand]),Sc=S2[,rand]+rnorm(nrow(S2),0,S2.sd[,rand]),Omega.sc=Omega.s2[,rand]+rnorm(nrow(Omega.s2),0,Omega.s2.sd[,rand]))
	MgCav<-MgCa2+rnorm(1,0,MgCa2*rsd)
	rmse_2[i,3]<-sqrt(sum((exp(predict(fit,v))-MgCav)^2)/nrow(v))					#rmse between MgCa predicted from the calibration for T2 and validation MgCa
	ce_2[i,3]<-1-(sum((exp(predict(fit,v))-MgCav)^2)/sum((mean(MgCav)-MgCav)^2))
	MgCav.pred[,i]<-exp(predict(fit,v))
	}
	
colMeans(r2_1)
apply(r2_1,2,sd)/sqrt(n.it)*2.576
colMeans(rmse_1)
apply(rmse_1,2,sd)/sqrt(n.it)*2.576
colMeans(rmse_2)
apply(rmse_2,2,sd)/sqrt(n.it)*2.576
colMeans(ce_2)
apply(ce_2,2,sd)/sqrt(n.it)*2.576
mean(r2_12)
sd(r2_12)/sqrt(n.it)*2.576
mean(rmse_12)
sd(rmse_12)/sqrt(n.it)*2.576
mean(rmse_3)
sd(rmse_3)/sqrt(n.it)*2.576
mean(ce_3)
sd(ce_3)/sqrt(n.it)*2.576
