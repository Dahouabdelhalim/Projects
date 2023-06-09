#final models of foraminifera MgCa using all data

setwd("/Volumes/GoogleDrive/My Drive/R_documents/PSR");

#Resutls from formMgCa_cal.val suggest the following models: 
#N. pachyderma (N. pachy L): 	f(Omega.s,Omega.d, and size)
#N. incompta (N. pachy R): 		f(S,size)
#N. pachyderma & incompta:		f(T,size)
#G. ruber: 						f(T,Omega.d) 
#G. inflata: 					f(T,Omega.s)
#G. bulloides: 					f(T,Omega.d)
#Construct models using the variables above and ALL Mg/Ca for a species

source("/Volumes/GoogleDrive/My Drive/R_documents/funcs/cal.val.R")
source("/Volumes/GoogleDrive/My Drive/R_documents/funcs/cal.val.nls.R")

library(ncdf4)
library(leaps)
library(abind)
library(plyr)

MgCa.rsd<-0.015																#assume RSD of Mg/Ca analyses is 1.5%
n.it<-1000

#coretop foram MgCa data
MgCa.dat.Npachy<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/Npachy.MgCa.env.dat.csv")
MgCa.dat.Gruber<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/Gruber.MgCa.env.dat.csv")
MgCa.dat.Ginflata<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/Ginflata.MgCa.env.dat.csv")
MgCa.dat.Gbulloides<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/Gbulloides.MgCa.env.dat.csv")

#empirically account for reductive cleaning following Rosenthal et al., 2004 relationship
red<-which(MgCa.dat.Npachy$reductive=="Y")
MgCa.dat.Npachy$Mg.Ca[red]<-MgCa.dat.Npachy$Mg.Ca[red]/0.85									#omit intercept for N. pachy as this brings Yu, 2008 and E+G, 2000 into better agreement
red<-which(MgCa.dat.Gruber$reductive=="Y")
MgCa.dat.Gruber$Mg.Ca[red]<-(MgCa.dat.Gruber$Mg.Ca[red]-0.2)/0.85
red<-which(MgCa.dat.Ginflata$reductive=="Y")
MgCa.dat.Ginflata$Mg.Ca[red]<-(MgCa.dat.Ginflata$Mg.Ca[red]-0.2)/0.85
red<-which(MgCa.dat.Gbulloides$reductive=="Y")
MgCa.dat.Gbulloides$Mg.Ca[red]<-(MgCa.dat.Gbulloides$Mg.Ca[red]-0.2)/0.85

#subset Npachy
MgCa.dat.NpachyR<-MgCa.dat.Npachy[which(MgCa.dat.Npachy$Notes=="Dextral"),]
MgCa.dat.NpachyL<-MgCa.dat.Npachy[which(MgCa.dat.Npachy$Notes=="Sinstral"),]

#***********************************************#
#coretop calibration of all N.pachyL
temp<-MgCa.dat.NpachyL[,20:26]
temp.sd<-MgCa.dat.NpachyL[,62:68]
salt<-MgCa.dat.NpachyL[,13:19]
salt.sd<-MgCa.dat.NpachyL[,55:61]
Omega.s<-MgCa.dat.NpachyL[,27:33]
Omega.s.sd<-MgCa.dat.NpachyL[,69:75]
Omega.d<-MgCa.dat.NpachyL$GLODAP.OmegaC.deep.cut
Omega.d.sd<-MgCa.dat.NpachyL$GLODAP.OmegaC.deep.sd.1
size<-MgCa.dat.NpachyL$mean_size
MgCa<-MgCa.dat.NpachyL$Mg.Ca


tmp<-matrix(,nrow=n.it,ncol=5)
tmp.pred<-matrix(,nrow=n.it,length(MgCa))
tmp.t<-matrix(,nrow=n.it,ncol=5)
tmp.t.pred<-matrix(,nrow=n.it,length(MgCa))

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																						
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)

	#use coefficients from stepwise model as starting point
	fit<-nls(MgCa.tmp~a*Omega.s.tmp^b1+b2*Omega.d.tmp+b3*size,start=c(a=0.3,b1=0.7,b2=0.3,b3=-0.001))
	tmp[i,1]<-summary(fit)$coeff[1,1]
	tmp[i,2]<-summary(fit)$coeff[2,1]
	tmp[i,3]<-summary(fit)$coeff[3,1]
	tmp[i,4]<-summary(fit)$coeff[4,1]
	tmp[i,5]<-cor(MgCa.tmp,predict(fit))
	tmp.pred[i,]<-predict(fit)
	fit<-lm(log(MgCa.tmp)~temp.tmp)
	tmp.t[i,1]<-summary(fit)$coeff[1,1]
	tmp.t[i,2]<-summary(fit)$coeff[2,1]
	tmp.t[1,3]<-NA
	tmp.t[1,4]<-NA
	tmp.t[i,5]<-summary(fit)$r.squared
	tmp.t.pred[i,]<-predict(fit)
	}	

NpachyL.final<-matrix(,nrow=2,ncol=ncol(tmp)*2)
NpachyL.final[1,seq(1,length.out=ncol(tmp),by=2)]<-colMeans(tmp)
NpachyL.final[1,seq(2,length.out=ncol(tmp),by=2)]<-apply(tmp,2,sd)/sqrt(n.it)*2.576
NpachyL.final[2,seq(1,length.out=ncol(tmp),by=2)]<-colMeans(tmp.t)
NpachyL.final[2,seq(2,length.out=ncol(tmp),by=2)]<-apply(tmp.t,2,sd)/sqrt(n.it)*2.576
colnames(NpachyL.final)<-c("int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","r","r_99ci")
NpachyL.final<-as.data.frame(NpachyL.final,stringsAsFactors=FALSE)
NpachyL.final$variables<-c("Omega.s,Omega.d,size","T")

dev.new(width=8, height=11);
par(mfrow=c(3,2))
plot(MgCa,colMeans(tmp.pred),pch=16,xlim=c(0,1.5),ylim=c(0,1.5),main="N. pachyderma",las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",col="palegreen3")
points(MgCa,exp(colMeans(tmp.t.pred)),pch=16)
abline(0,1)
legend("topleft",c("multivariate (eq. 1)","T-only"),col=c("palegreen3","black"),pch=c(16,16),bty="n")

#***********************************************#
#coretop calibration of all N.pachyR, a.k.a N. incompta
temp<-MgCa.dat.NpachyR[,20:26]
temp.sd<-MgCa.dat.NpachyR[,62:68]
salt<-MgCa.dat.NpachyR[,13:19]
salt.sd<-MgCa.dat.NpachyR[,55:61]
MgCa<-MgCa.dat.NpachyR$Mg.Ca
size<-MgCa.dat.NpachyR$mean_size

tmp<-matrix(,nrow=n.it,ncol=4)
tmp.pred<-matrix(,nrow=n.it,length(MgCa))
tmp.t<-matrix(,nrow=n.it,ncol=4)
tmp.t.pred<-matrix(,nrow=n.it,length(MgCa))

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																	
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	#use coefficients from stepwise model as starting point
	fit<-nls(MgCa.tmp~a*exp(b1*salt.tmp)+b2*size,start=c(a=1e-6,b1=0.36,b2=0.001))
	tmp[i,1]<-summary(fit)$coeff[1,1]
	tmp[i,2]<-summary(fit)$coeff[2,1]
	tmp[i,3]<-summary(fit)$coeff[3,1]
	tmp[i,4]<-cor(MgCa.tmp,predict(fit))
	tmp.pred[i,]<-predict(fit)
	fit<-lm(log(MgCa.tmp)~temp.tmp)
	tmp.t[i,1]<-summary(fit)$coeff[1,1]
	tmp.t[i,2]<-summary(fit)$coeff[2,1]
	tmp.t[1,3]<-NA
	tmp.t[i,4]<-summary(fit)$r.squared
	tmp.t.pred[i,]<-predict(fit)
	}	

NpachyR.final<-matrix(,nrow=2,ncol=ncol(tmp)*2)
NpachyR.final[1,seq(1,length.out=ncol(tmp),by=2)]<-colMeans(tmp)
NpachyR.final[1,seq(2,length.out=ncol(tmp),by=2)]<-apply(tmp,2,sd)/sqrt(n.it)*2.576
NpachyR.final[2,seq(1,length.out=ncol(tmp),by=2)]<-colMeans(tmp.t)
NpachyR.final[2,seq(2,length.out=ncol(tmp),by=2)]<-apply(tmp.t,2,sd)/sqrt(n.it)*2.576
colnames(NpachyR.final)<-c("int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","r","r_99ci")
NpachyR.final<-as.data.frame(NpachyR.final,stringsAsFactors=FALSE)
NpachyR.final$variables<-c("S,size","T")

plot(MgCa,colMeans(tmp.pred),pch=16,xlim=c(1,2.5),ylim=c(1,2.5),main="N. incompta",las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",col="grey50")
points(MgCa,exp(colMeans(tmp.t.pred)),pch=16)
abline(0,1)
legend("topleft",c("multivariate (eq. 2)","T-only"),col=c("grey50","black"),pch=c(16,16),bty="n")	

#***********************************************#
#coretop calibration of all N.pachy
temp<-MgCa.dat.Npachy[,20:26]
temp.sd<-MgCa.dat.Npachy[,62:68]
size<-MgCa.dat.Npachy$mean_size
MgCa<-MgCa.dat.Npachy$Mg.Ca

tmp<-matrix(,nrow=n.it,ncol=4)
tmp.pred<-matrix(,nrow=n.it,length(MgCa))
tmp.t<-matrix(,nrow=n.it,ncol=4)
tmp.t.pred<-matrix(,nrow=n.it,length(MgCa))

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																	
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	#use coefficients from stepwise model as starting point
	fit<-nls(MgCa.tmp~a*exp(b1*temp.tmp)+b2*size,start=c(a=0.1,b1=0.1,b2=0.002))
	tmp[i,1]<-summary(fit)$coeff[1,1]
	tmp[i,2]<-summary(fit)$coeff[2,1]
	tmp[i,3]<-summary(fit)$coeff[3,1]
	tmp[i,4]<-cor(MgCa.tmp,predict(fit))
	tmp.pred[i,]<-predict(fit)
	fit<-lm(log(MgCa.tmp)~temp.tmp)
	tmp.t[i,1]<-summary(fit)$coeff[1,1]
	tmp.t[i,2]<-summary(fit)$coeff[2,1]
	tmp.t[1,3]<-NA
	tmp.t[i,4]<-summary(fit)$r.squared
	tmp.t.pred[i,]<-predict(fit)
	}	

Npachy.final<-matrix(,nrow=2,ncol=ncol(tmp)*2)
Npachy.final[1,seq(1,length.out=ncol(tmp),by=2)]<-colMeans(tmp)
Npachy.final[1,seq(2,length.out=ncol(tmp),by=2)]<-apply(tmp,2,sd)/sqrt(n.it)*2.576
Npachy.final[2,seq(1,length.out=ncol(tmp),by=2)]<-colMeans(tmp.t)
Npachy.final[2,seq(2,length.out=ncol(tmp),by=2)]<-apply(tmp.t,2,sd)/sqrt(n.it)*2.576
colnames(Npachy.final)<-c("int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","r","r_99ci")
Npachy.final<-as.data.frame(Npachy.final,stringsAsFactors=FALSE)
Npachy.final$variables<-c("T,size","T")

plot(MgCa,colMeans(tmp.pred),pch=16,xlim=c(0,2.5),ylim=c(0,2.5),main="N. pachyderma & incompta",las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",col="darkmagenta")
points(MgCa,exp(colMeans(tmp.t.pred)),pch=16)
abline(0,1)
legend("topleft",c("multivariate (eq. 6)","T-only"),col=c("darkmagenta","black"),pch=c(16,16),bty="n")	

#***********************************************#
#coretop calibration of all Gruber
temp<-MgCa.dat.Gruber[,18:22]
temp.sd<-MgCa.dat.Gruber[,48:52]
Omega.d<-MgCa.dat.Gruber$GLODAP.OmegaC.deep.cut
Omega.d.sd<-MgCa.dat.Gruber$GLODAP.OmegaC.deep.sd.1
MgCa<-MgCa.dat.Gruber$Mg.Ca

tmp<-matrix(,nrow=n.it,ncol=4)
tmp.pred<-matrix(,nrow=n.it,length(MgCa))
tmp.t<-matrix(,nrow=n.it,ncol=4)
tmp.t.pred<-matrix(,nrow=n.it,length(MgCa))

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	#use coefficients from stepwise model as starting point
	fit<-nls(MgCa.tmp~a*exp(b1*temp.tmp)+b2*Omega.d.tmp,start=c(a=0.5,b1=0.06,b2=1))
	tmp[i,1]<-summary(fit)$coeff[1,1]
	tmp[i,2]<-summary(fit)$coeff[2,1]
	tmp[i,3]<-summary(fit)$coeff[3,1]
	tmp[i,4]<-cor(MgCa.tmp,predict(fit))
	tmp.pred[i,]<-predict(fit)
	fit<-lm(log(MgCa.tmp)~temp.tmp)
	tmp.t[i,1]<-summary(fit)$coeff[1,1]
	tmp.t[i,2]<-summary(fit)$coeff[2,1]
	tmp.t[1,3]<-NA
	tmp.t[i,4]<-summary(fit)$r.squared
	tmp.t.pred[i,]<-predict(fit)
	}	

Gruber.final<-matrix(,nrow=2,ncol=ncol(tmp)*2)
Gruber.final[1,seq(1,length.out=ncol(tmp),by=2)]<-round(colMeans(tmp),7)
Gruber.final[1,seq(2,length.out=ncol(tmp),by=2)]<-round(apply(tmp,2,sd)/sqrt(n.it)*2.576,7)
Gruber.final[2,seq(1,length.out=ncol(tmp),by=2)]<-round(colMeans(tmp.t),3)
Gruber.final[2,seq(2,length.out=ncol(tmp),by=2)]<-round(apply(tmp.t,2,sd)/sqrt(n.it)*2.576,3)
colnames(Gruber.final)<-c("int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","r","r_99ci")
Gruber.final<-as.data.frame(Gruber.final,stringsAsFactors=FALSE)
Gruber.final$variables<-c("T,Omega.d","T")

plot(MgCa,colMeans(tmp.pred),pch=16,xlim=c(1.5,7),ylim=c(1.5,7),main="G. ruber",las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",col="cornflowerblue")
points(MgCa,exp(colMeans(tmp.t.pred)),pch=16)
abline(0,1)
legend("topleft",c("multivariate (eq. 3)","T-only"),col=c("cornflowerblue","black"),pch=c(16,16),bty="n")	

#***********************************************#
#coretop calibration of all Ginflata
temp<-MgCa.dat.Ginflata[,17:20]
temp.sd<-MgCa.dat.Ginflata[,41:44]
Omega.s<-MgCa.dat.Ginflata[,21:24]
Omega.s.sd<-MgCa.dat.Ginflata[,45:48]
MgCa<-MgCa.dat.Ginflata$Mg.Ca

tmp<-matrix(,nrow=n.it,ncol=4)
tmp.pred<-matrix(,nrow=n.it,length(MgCa))
tmp.t<-matrix(,nrow=n.it,ncol=4)
tmp.t.pred<-matrix(,nrow=n.it,length(MgCa))

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	#use coefficients from stepwise model as starting point
	fit<-nls(MgCa.tmp~a*exp(b1*temp.tmp)*Omega.s.tmp^b2,start=c(a=1.5,b1=0.06,b2=-1))
	tmp[i,1]<-summary(fit)$coeff[1,1]
	tmp[i,2]<-summary(fit)$coeff[2,1]
	tmp[i,3]<-summary(fit)$coeff[3,1]
	tmp[i,4]<-cor(MgCa.tmp,predict(fit))
	tmp.pred[i,]<-predict(fit)
	fit<-lm(log(MgCa.tmp)~temp.tmp)
	tmp.t[i,1]<-summary(fit)$coeff[1,1]
	tmp.t[i,2]<-summary(fit)$coeff[2,1]
	tmp.t[1,3]<-NA
	tmp.t[i,4]<-summary(fit)$r.squared
	tmp.t.pred[i,]<-predict(fit)
	}	

Ginflata.final<-matrix(,nrow=2,ncol=ncol(tmp)*2)
Ginflata.final[1,seq(1,length.out=ncol(tmp),by=2)]<-round(colMeans(tmp),7)
Ginflata.final[1,seq(2,length.out=ncol(tmp),by=2)]<-round(apply(tmp,2,sd)/sqrt(n.it)*2.576,7)
Ginflata.final[2,seq(1,length.out=ncol(tmp),by=2)]<-round(colMeans(tmp.t),3)
Ginflata.final[2,seq(2,length.out=ncol(tmp),by=2)]<-round(apply(tmp.t,2,sd)/sqrt(n.it)*2.576,3)
colnames(Ginflata.final)<-c("int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","r","r_99ci")
Ginflata.final<-as.data.frame(Ginflata.final,stringsAsFactors=FALSE)
Ginflata.final$variables<-c("T,Omega.s","T")

plot(MgCa,colMeans(tmp.pred),pch=16,xlim=c(0.5,3),ylim=c(0.5,3),main="G. inflata",las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",col="firebrick")
points(MgCa,exp(colMeans(tmp.t.pred)),pch=16)
abline(0,1)
legend("topleft",c("multivariate (eq. 4)","T-only"),col=c("firebrick","black"),pch=c(16,16),bty="n")	

#***********************************************#
#coretop calibration of all Gbulloides
temp<-MgCa.dat.Gbulloides[,20:26]
temp.sd<-MgCa.dat.Gbulloides[,62:68]
Omega.d<-MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.cut
Omega.d.sd<-MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.sd.1
MgCa<-MgCa.dat.Gbulloides$Mg.Ca

tmp<-matrix(,nrow=n.it,ncol=4)
tmp.pred<-matrix(,nrow=n.it,length(MgCa))
tmp.t<-matrix(,nrow=n.it,ncol=4)
tmp.t.pred<-matrix(,nrow=n.it,length(MgCa))

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	#use coefficients from stepwise model as starting point
	fit<-nls(MgCa.tmp~a*exp(b1*temp.tmp)+b2*Omega.d.tmp,start=c(a=0.5,b1=0.06,b2=1))
	tmp[i,1]<-summary(fit)$coeff[1,1]
	tmp[i,2]<-summary(fit)$coeff[2,1]
	tmp[i,3]<-summary(fit)$coeff[3,1]
	tmp[i,4]<-cor(MgCa.tmp,predict(fit))
	tmp.pred[i,]<-predict(fit)
	fit<-lm(log(MgCa.tmp)~temp.tmp)
	tmp.t[i,1]<-summary(fit)$coeff[1,1]
	tmp.t[i,2]<-summary(fit)$coeff[2,1]
	tmp.t[1,3]<-NA
	tmp.t[i,4]<-summary(fit)$r.squared
	tmp.t.pred[i,]<-predict(fit)
	}	

Gbulloides.final<-matrix(,nrow=2,ncol=ncol(tmp)*2)
Gbulloides.final[1,seq(1,length.out=ncol(tmp),by=2)]<-round(colMeans(tmp),7)
Gbulloides.final[1,seq(2,length.out=ncol(tmp),by=2)]<-round(apply(tmp,2,sd)/sqrt(n.it)*2.576,7)
Gbulloides.final[2,seq(1,length.out=ncol(tmp),by=2)]<-round(colMeans(tmp.t),3)
Gbulloides.final[2,seq(2,length.out=ncol(tmp),by=2)]<-round(apply(tmp.t,2,sd)/sqrt(n.it)*2.576,3)
colnames(Gbulloides.final)<-c("int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","r","r_99ci")
Gbulloides.final<-as.data.frame(Gbulloides.final,stringsAsFactors=FALSE)
Gbulloides.final$variables<-c("T,Omega.d","T")

plot(MgCa,colMeans(tmp.pred),pch=16,xlim=c(0.5,4.5),ylim=c(0.5,4.5),main="G. inflata",las=1,xlab="measured Mg/Ca (mmol/mol)",ylab="predicted Mg/Ca (mmol/mol)",col="gold3")
points(MgCa,exp(colMeans(tmp.t.pred)),pch=16)
abline(0,1)
legend("topleft",c("multivariate (eq. 4)","T-only"),col=c("gold3","black"),pch=c(16,16),bty="n")

write.csv(NpachyL.final,"proxy_data/NpachyL_alldata.cal.csv")
write.csv(NpachyR.final,"proxy_data/NpachyR_alldata.cal.csv")
write.csv(Npachy.final,"proxy_data/Npachy_alldata.cal.csv")
write.csv(Gruber.final,"proxy_data/Gruber_alldata.cal.csv")
write.csv(Ginflata.final,"proxy_data/Ginflata_alldata.cal.csv")
write.csv(Gbulloides.final,"proxy_data/Gbulloides_alldata.cal.csv")
