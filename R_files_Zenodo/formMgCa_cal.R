#calibrate and validate forward models of foraminifera MgCa using an interative approach with CE as metric 

setwd("/Volumes/GoogleDrive/My Drive/R_documents/PSR");

source("/Volumes/GoogleDrive/My Drive/R_documents/funcs/cal.val.R")
source("/Volumes/GoogleDrive/My Drive/R_documents/funcs/cal.val.nls.R")

library(ncdf4)
library(leaps)
library(abind)
library(plyr)

#constants
MgCa.rsd<-0.015																#assume RSD of Mg/Ca analyses is 1.5%																				
n.it<-10000																	#number of iterations
prct.cal<-0.75																#percent of data in calibration


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
#coretop calibration of all N.pachyR, a.k.a N. incompta
temp<-MgCa.dat.NpachyR[,20:26]
temp.sd<-MgCa.dat.NpachyR[,62:68]
salt<-MgCa.dat.NpachyR[,13:19]
salt.sd<-MgCa.dat.NpachyR[,55:61]
Omega.s<-MgCa.dat.NpachyR[,27:33]
Omega.s.sd<-MgCa.dat.NpachyR[,69:75]
pH.s<-MgCa.dat.NpachyR[,34:40]
pH.s.sd<-MgCa.dat.NpachyR[,76:82]
Omega.d<-MgCa.dat.NpachyR$GLODAP.OmegaC.deep.cut
Omega.d.sd<-MgCa.dat.NpachyR$GLODAP.OmegaC.deep.sd.1
pH.d<-MgCa.dat.NpachyR$GLODAP.pH.deep.cut
pH.d.sd<-MgCa.dat.NpachyR$GLODAP.pH.deep.sd.1
MgCa<-MgCa.dat.NpachyR$Mg.Ca
size<-MgCa.dat.NpachyR$mean_size

#univariate
nvar<-7																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																						
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	#cal.val.nls(salt.tmp, MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(1.5E-4,0.25))
	cal.val(salt.tmp,log(MgCa.tmp),prct.cal,1)
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	#cal.val.nls(temp.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(1,0.04))
	cal.val(temp.tmp,log(MgCa.tmp),prct.cal,1)
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*x1^b1},s=c(1,1))
	#cal.val(log(Omega.s.tmp),log(MgCa.tmp),prct.cal,1)
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(pH.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	#cal.val(pH.s.tmp,MgCa.tmp,prct.cal,1)
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	#cal.val(Omega.d.tmp,MgCa.tmp,prct.cal,1)
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(pH.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	#cal.val(pH.d.tmp,MgCa.tmp,prct.cal,1)
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
	cal.val.nls(size,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	#cal.val(size,MgCa.tmp,prct.cal,1)
	tmp[nvar*(i-1)+7,]<-c(7,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-colMeans(dat.tmp)[2:10]
	tmp.means[i,seq(2,18,2)]<-apply(dat.tmp,2,sd)[2:10]/sqrt(nrow(dat.tmp))*2.576
	f.good[i]<-nrow(dat.tmp)/n.it
}
colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyR.1<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyR.1$variables<-c("S","T","Omega.s","pH.s","Omega.d","pH.d","size")
NpachyR.1$f.good<-f.good

#bivariate, keeping salinity
nvar<-6																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																						
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)

	
	cal.val.nls(cbind(salt.tmp,temp.tmp), MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1+b2*x2)},s=c(1E-4,0.22,0.01))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(salt.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)*x2^b2},s=c(1E-4,0.23,0.2))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(salt.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(0.6,0.14,-0.4))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(salt.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1E-5,0.3,0.07))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(salt.tmp,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(-50,-0.03,2))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(salt.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1E-3,0.4,0.001))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}
colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyR.2<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyR.2$variables<-c("S,T","S,Omega.s","S,pH.s*","S,Omega.d","S,pH.d","S,size")
NpachyR.2$f.good<-f.good

#bivariate, keeping temperature
nvar<-6																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																						
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)

	
	cal.val.nls(cbind(temp.tmp,salt.tmp), MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1+b2*x2)},s=c(1E-4,0.02,0.2))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)*x2^b2},s=c(1,0.1,0.2))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(-0.1,0.15,0.1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1,0.02,0.02))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(-10,0.02,2))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1,0.02,0.001))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}
colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyR.2T<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyR.2T$variables<-c("T,S","T,Omega.s","T,pH.s*","T,Omega.d","T,pH.d","T,size")
NpachyR.2T$f.good<-f.good


#trivariate, keeping salinity and size
nvar<-5																					#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),1,)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(salt.tmp,temp.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1+b2*x2)+b3*x3},s=c(2.5,0.2,0.01,1E-4))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(salt.tmp,size,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)*x3^b3+b2*x2},s=c(1E-3,0.1,0.001,2))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(salt.tmp,size,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(1E-3,0.1,0.001,1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(salt.tmp,size,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(0.1,0.15,0.002,0.001))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(salt.tmp,size,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(1E-3,0.15,0.002,0.1))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyR.3<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyR.3$variables<-c("S,size,T*","S,size,Omega.s*","S,size,pH.s*","S,size,Omega.d*","S,size,pH.d*")
NpachyR.3$f.good<-f.good

#trivariate, keeping temperature, and salinity
nvar<-5																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(temp.tmp,salt.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1+b2*x2)*x3^b3},s=c(0.01,0.05,0.1,-0.1))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,salt.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1+b2*x2)+b3*x3},s=c(1E-5,0.13,0.25,0.15))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,salt.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1+b2*x2)+b3*x3},s=c(0.001,0.015,0.2,0.05))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,salt.tmp,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1+b2*x2)+b3*x3},s=c(-50,-0.015,-0.2,1))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,salt.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1+b2*x2)+b3*x3},s=c(1,0.1,0.01,1E-5))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyR.3T<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyR.3T$variables<-c("T,S,Omega.s","T,S,pH.s*","T,S,Omega.d*","T,S,pH.d*","T,S,size*")
NpachyR.3T$f.good<-f.good

#quadvariate, keeping temperature, salinity, Omega.s
nvar<-3																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(temp.tmp,salt.tmp,Omega.s.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,x4,a,b1,b2,b3,b4) {a*exp(b1*x1+b2*x2)*x3^b3+b4*x4},s=c(0.01,0.05,0.1,-0.1,0.05))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,salt.tmp,Omega.s.tmp,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,x4,a,b1,b2,b3,b4) {a*exp(b1*x1+b2*x2)*x3^b3+b4*x4},s=c(0.01,0.05,0.1,-0.1,0.5))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,salt.tmp,Omega.s.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,x4,a,b1,b2,b3,b4) {a*exp(b1*x1+b2*x2)*x3^b3+b4*x4},s=c(0.01,0.05,0.1,-0.1,0.001))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyR.4T<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyR.4T$variables<-c("T,S,Omega.s,Omega.d*","T,S,Omega.s,pH.d*","T,S,Omega.s,size*")
NpachyR.4T$f.good<-f.good


NpachyR.all<-rbind(NpachyR.1,NpachyR.2,NpachyR.3)
NpachyR.T.all<-rbind(NpachyR.1,NpachyR.2T,NpachyR.3T,NpachyR.4T)

#***********************************************#
#coretop calibration of all N.pachyL
temp<-MgCa.dat.NpachyL[,20:26]
temp.sd<-MgCa.dat.NpachyL[,62:68]
salt<-MgCa.dat.NpachyL[,13:19]
salt.sd<-MgCa.dat.NpachyL[,55:61]
Omega.s<-MgCa.dat.NpachyL[,27:33]
Omega.s.sd<-MgCa.dat.NpachyL[,69:75]
pH.s<-MgCa.dat.NpachyL[,34:40]
pH.s.sd<-MgCa.dat.NpachyL[,76:82]
Omega.d<-MgCa.dat.NpachyL$GLODAP.OmegaC.deep.cut
Omega.d.sd<-MgCa.dat.NpachyL$GLODAP.OmegaC.deep.sd.1
pH.d<-MgCa.dat.NpachyL$GLODAP.pH.deep.cut
pH.d.sd<-MgCa.dat.NpachyL$GLODAP.pH.deep.sd.1
MgCa<-MgCa.dat.NpachyL$Mg.Ca
size<-MgCa.dat.NpachyL$mean_size

#univariate
nvar<-7																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																						
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(salt.tmp, MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(0.001,0.2))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(temp.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(1,0.04))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*x1^b1},s=c(1,1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(pH.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(pH.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
	cal.val.nls(size,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+7,]<-c(7,round(cal.valsummary[,1],3))

}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}
colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyL.1<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyL.1$variables<-c("S","T","Omega.s","pH.s","Omega.d","pH.d","size")
NpachyL.1$f.good<-f.good

#bivariate, keeping Omega.d
nvar<-5																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(salt.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(-10,-0.01,1))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(-10,0.1,1))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(Omega.d.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a+b1*x1+b2*x2},s=c(-10,1,-0.1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(Omega.d.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*x2^b2+b1*x1},s=c(-10,1,-0.1))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(Omega.d.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a+b1*x1+b2*x2},s=c(-10,1,0.00001))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyL.2<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyL.2$variables<-c("S,Omega.d*","T,Omega.d","Omega.d,pH.s","Omega.d,Omega.s*","Omega.d,size")
NpachyL.2$f.good<-f.good

#trivariate, keeping Omega.d, and size
nvar<-4																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(salt.tmp,Omega.d.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(0.1,0.06,0.25,-0.001))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(0.5,0.05,0.5,-0.001))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(Omega.d.tmp,pH.s.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a+b1*x1+b2*x2+b3*x3},s=c(-5,0.3,0.5,-0.001))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(Omega.d.tmp,Omega.s.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*x2^b2+b1*x1+b3*x3},s=c(0.1,0.3,0.2,-0.001))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyL.3<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyL.3$variables<-c("S,Omega.d,size*","T,Omega.d,size","Omega.d,pH.s,size","Omega.d,Omega.s,size")
NpachyL.3$f.good<-f.good

#quadvariate, keeping Omega.d, Omega.s, and size
nvar<-2																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(salt.tmp,Omega.d.tmp,Omega.s.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,x4,a,b1,b2,b3,b4) {a*exp(b1*x1)*x3^b3+b2*x2+b4*x4},s=c(0.4,0.04,0.3,0.1,-0.001))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,Omega.s.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,x4,a,b1,b2,b3,b4) {a*exp(b1*x1)*x3^b3+b2*x2+b4*x4},s=c(0.4,0.05,0.4,-0.1,-0.001))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyL.4<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyL.4$variables<-c("S,Omega.d,Omega.s,size","T,Omega.d,Omega.s,size")
NpachyL.4$f.good<-f.good

NpachyL.all<-rbind(NpachyL.1,NpachyL.2,NpachyL.3,NpachyL.4)

#***********************************************#
#coretop calibration of all N.pachyL with T>3
temp<-MgCa.dat.NpachyL[,20:26]
MgCa.dat.NpachyLw<-MgCa.dat.NpachyL[which(rowMeans(temp)>3),]
temp<-MgCa.dat.NpachyLw[,20:26]
temp.sd<-MgCa.dat.NpachyLw[,62:68]
salt<-MgCa.dat.NpachyLw[,13:19]
salt.sd<-MgCa.dat.NpachyLw[,55:61]
Omega.s<-MgCa.dat.NpachyLw[,27:33]
Omega.s.sd<-MgCa.dat.NpachyLw[,69:75]
pH.s<-MgCa.dat.NpachyLw[,34:40]
pH.s.sd<-MgCa.dat.NpachyLw[,76:82]
Omega.d<-MgCa.dat.NpachyLw$GLODAP.OmegaC.deep.cut
Omega.d.sd<-MgCa.dat.NpachyLw$GLODAP.OmegaC.deep.sd.1
pH.d<-MgCa.dat.NpachyLw$GLODAP.pH.deep.cut
pH.d.sd<-MgCa.dat.NpachyLw$GLODAP.pH.deep.sd.1
MgCa<-MgCa.dat.NpachyLw$Mg.Ca
size<-MgCa.dat.NpachyLw$mean_size

#univariate
nvar<-7																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																						
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(salt.tmp, MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(0.001,0.2))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(temp.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(1,0.04))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*x1^b1},s=c(1,1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(pH.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(pH.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
	cal.val.nls(size,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+7,]<-c(7,round(cal.valsummary[,1],3))

}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}
colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
NpachyLw.1<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
NpachyLw.1$variables<-c("S","T","Omega.s","pH.s","Omega.d","pH.d","size")
NpachyLw.1$f.good<-f.good


#***********************************************#
#coretop calibration of all N.pachy
temp<-MgCa.dat.Npachy[,20:26]
temp.sd<-MgCa.dat.Npachy[,62:68]
salt<-MgCa.dat.Npachy[,13:19]
salt.sd<-MgCa.dat.Npachy[,55:61]
Omega.s<-MgCa.dat.Npachy[,27:33]
Omega.s.sd<-MgCa.dat.Npachy[,69:75]
pH.s<-MgCa.dat.Npachy[,34:40]
pH.s.sd<-MgCa.dat.Npachy[,76:82]
Omega.d<-MgCa.dat.Npachy$GLODAP.OmegaC.deep.cut
Omega.d.sd<-MgCa.dat.Npachy[,83]
pH.d<-MgCa.dat.Npachy$GLODAP.pH.deep.cut
pH.d.sd<-MgCa.dat.Npachy[,90]
MgCa<-MgCa.dat.Npachy$Mg.Ca
#Nurnberg does not provide size, assume mean of other studies
#nans<-which(is.na(MgCa.dat.Npachy$mean_size)==TRUE)
#MgCa.dat.Npachy$mean_size[nans]<-round(mean(MgCa.dat.Npachy$mean_size,na.rm=T))
size<-MgCa.dat.Npachy$mean_size

#univariate
nvar<-7																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(salt.tmp, MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(1E-3,0.5))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(temp.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(1,0.04))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*x1^b1},s=c(1,1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(pH.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(pH.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
	cal.val.nls(size,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+7,]<-c(7,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Npachy.1<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Npachy.1$variables<-c("S","T","Omega.s","pH.s","Omega.d","pH.d","size")
Npachy.1$f.good<-f.good
	
#bivariate, keeping salinity
nvar<-6																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(salt.tmp,temp.tmp), MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1+b2*x2)},s=c(0.01,0.4,0.05))
	tmp[nvar*(i-1)+1,]<-c(1,cal.valsummary[,1])
	cal.val.nls(cbind(salt.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)*x2^b2},s=c(1E-3,0.3,0.6))
	tmp[nvar*(i-1)+2,]<-c(2,cal.valsummary[,1])
	cal.val.nls(cbind(salt.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1E-3,0.2,-0.1))
	tmp[nvar*(i-1)+3,]<-c(3,cal.valsummary[,1])
	cal.val.nls(cbind(salt.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(-10,0.3,0.1))
	tmp[nvar*(i-1)+4,]<-c(4,cal.valsummary[,1],3)
	cal.val.nls(cbind(salt.tmp,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1E-3,0.2,-0.1))
	tmp[nvar*(i-1)+5,]<-c(5,cal.valsummary[,1],3)
	cal.val.nls(cbind(salt.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(-10,0.3,0.00001))
	tmp[nvar*(i-1)+6,]<-c(6,cal.valsummary[,1],3)
	}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Npachy.2<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Npachy.2$variables<-c("S,T","S,Omega.s*","S,pH.s*","S,Omega.d*","S,pH.d*","S,size*")
Npachy.2$f.good<-f.good

#bivariate, forcing temperature to be kept
nvar<-6																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(salt.tmp,temp.tmp), MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1+b2*x2)},s=c(0.01,0.4,0.05))
	tmp[nvar*(i-1)+1,]<-c(1,cal.valsummary[,1])
	cal.val.nls(cbind(temp.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)*x2^b2},s=c(1E-3,0.3,0.6))
	tmp[nvar*(i-1)+2,]<-c(2,cal.valsummary[,1])
	cal.val.nls(cbind(temp.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1E-3,0.2,-0.1))
	tmp[nvar*(i-1)+3,]<-c(3,cal.valsummary[,1])
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(-10,0.3,0.1))
	tmp[nvar*(i-1)+4,]<-c(4,cal.valsummary[,1])
	cal.val.nls(cbind(temp.tmp,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1E-3,0.2,-0.1))
	tmp[nvar*(i-1)+5,]<-c(5,cal.valsummary[,1])
	cal.val.nls(cbind(temp.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(-10,0.3,0.00001))
	tmp[nvar*(i-1)+6,]<-c(6,cal.valsummary[,1])
	}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-colMeans(dat.tmp,na.rm=T)[2:10]
	tmp.means[i,seq(2,18,2)]<-apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Npachy.2T<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Npachy.2T$variables<-c("S,T","T,Omega.s","T,pH.s*","T,Omega.d","T,pH.d","T,size")
Npachy.2T$f.good<-f.good


#trivariate, keeping temperature and size
nvar<-5																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(temp.tmp,salt.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1+b2*x2)+b3*x3},s=c(1E-3,0.05,0.4,0.001))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.s.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)*x2^b2+b3*x3},s=c(0.1,0.1,-0.1,0.001))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,pH.s.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(-1,0.1,-0.1,0.001))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(-1E-3,0.1,0.1,0.002))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,pH.d.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(-1,0.1,1.7,0.002))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Npachy.3T<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Npachy.3T$variables<-c("T,S,size*","T,Omega.s,size*","T,pH.s,size","T,Omega.d,size*","T,pH.d,size*")
Npachy.3T$f.good<-f.good

Npachy.all<-rbind(Npachy.1,Npachy.2,Npachy.3)
#***********************************************#
#coretop calibration of G.ruber
temp<-MgCa.dat.Gruber[,18:22]
temp.sd<-MgCa.dat.Gruber[,48:52]
salt<-MgCa.dat.Gruber[,13:17]
salt.sd<-MgCa.dat.Gruber[,43:47]
Omega.s<-MgCa.dat.Gruber[,23:27]
Omega.s.sd<-MgCa.dat.Gruber[,53:57]
pH.s<-MgCa.dat.Gruber[,28:32]
pH.s.sd<-MgCa.dat.Gruber[,58:62]
Omega.d<-MgCa.dat.Gruber$GLODAP.OmegaC.deep.cut
Omega.d.sd<-MgCa.dat.Gruber$GLODAP.OmegaC.deep.sd.1
pH.d<-MgCa.dat.Gruber$GLODAP.pH.deep.cut
pH.d.sd<-MgCa.dat.Gruber$GLODAP.pH.deep.sd.1
MgCa<-MgCa.dat.Gruber$Mg.Ca
size<-MgCa.dat.Gruber$mean_size

#univariate
nvar<-7																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(salt.tmp, MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(5,-0.1))
	#cal.val(salt.tmp, log(MgCa.tmp),prct.cal,1)
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(temp.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(1,0.07))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*x1^b1},s=c(1,1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(pH.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(pH.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
	cal.val.nls(size,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+7,]<-c(7,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Gruber.1<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Gruber.1$variables<-c("S","T","Omega.s","pH.s","Omega.d","pH.d","size")
Gruber.1$f.good<-f.good

#bivariate, keeping temperature
nvar<-6																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(temp.tmp,salt.tmp), MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1+b2*x2)},s=c(1,0.1,0.05))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)*x2^b2},s=c(0.01,0.1,0.1))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(0.1,0.1,-2))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(0.1,0.1,2))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(-50,0.1,-2))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1,0.1,-0.001))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
	}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Gruber.2<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Gruber.2$variables<-c("T,S","T,Omega.s","T,pH.s","T,Omega.d","T,pH.d","T,size")
Gruber.2$f.good<-f.good

#trivariate, keeping temperature and Omega.d
nvar<-4																					#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),1,)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,salt.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1+b3*x3)+b2*x2},s=c(0.1,0.1,2,0.1))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)*x3^b3+b2*x2},s=c(10,0.1,2,-0.01))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(10,0.1,2,-3))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(10,0.1,2,-0.01))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Gruber.3<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Gruber.3$variables<-c("T,Omega.d,S","T,Omega.d,Omega.s*","T,Omega.d,pH.s*","T,Omega.d,size")
Gruber.3$f.good<-f.good

Gruber.all<-rbind(Gruber.1,Gruber.2,Gruber.3)

#***********************************************#
#coretop calibration of G.inflata
temp<-MgCa.dat.Ginflata[,17:20]
temp.sd<-MgCa.dat.Ginflata[,41:44]
salt<-MgCa.dat.Ginflata[,13:16]
salt.sd<-MgCa.dat.Ginflata[,37:40]
Omega.s<-MgCa.dat.Ginflata[,21:24]
Omega.s.sd<-MgCa.dat.Ginflata[,45:48]
pH.s<-MgCa.dat.Ginflata[,25:28]
pH.s.sd<-MgCa.dat.Ginflata[,49:52]
Omega.d<-MgCa.dat.Ginflata$GLODAP.OmegaC.deep.cut
Omega.d.sd<-MgCa.dat.Ginflata$GLODAP.OmegaC.deep.sd.1
pH.d<-MgCa.dat.Ginflata$GLODAP.pH.deep.cut
pH.d.sd<-MgCa.dat.Ginflata$GLODAP.pH.deep.sd.1
MgCa<-MgCa.dat.Ginflata$Mg.Ca
size<-MgCa.dat.Ginflata$mean_size

#univariate
nvar<-7																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(salt.tmp, MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(0.01,0.1))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(temp.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(0.01,0.1))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*x1^b1},s=c(0.1,-0.1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(pH.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,-1))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(0.1,0.1))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(pH.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(-1,1))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
	cal.val.nls(size,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,0.1))
	tmp[nvar*(i-1)+7,]<-c(7,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Ginflata.1<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Ginflata.1$variables<-c("S","T","Omega.s","pH.s","Omega.d","pH.d","size")
Ginflata.1$f.good<-f.good


#bivariate, keeping temperature
nvar<-6																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(temp.tmp,salt.tmp), MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1+b2*x2)},s=c(1,0.1,-0.05))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)*x2^b2},s=c(2,0.1,-0.1))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1,0.1,-1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(0.5,0.1,0.5))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(2,0.1,-0.1))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1,0.1,0.001))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
	}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Ginflata.2<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Ginflata.2$variables<-c("T,S*","T,Omega.s","T,pH.s*","T,Omega.d","T,pH.d*","T,size")
Ginflata.2$f.good<-f.good

#trivariate, keeping temperature and Omega.s.
nvar<-4																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(temp.tmp,Omega.s.tmp,salt.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1+b3*x3)*x2^b2},s=c(50,0.1,-1,-0.1))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.s.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)*x2^b2+b3*x3},s=c(1,0.1,-1,1))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.s.tmp,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)*x2^b2+b3*x3},s=c(2,0.1,-1,-0.1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.s.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)*x2^b2+b3*x3},s=c(2,0.1,-1,0.01))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}
colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Ginflata.3<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Ginflata.3$variables<-c("T,Omega.s,S*","T,Omega.s,Omega.d*","T,Omega.s,pH.d*","T,Omega.s,size")
Ginflata.3$f.good<-f.good

Ginflata.all<-rbind(Ginflata.1,Ginflata.2,Ginflata.3)
#***********************************************#
#coretop calibration of G.bulloides
temp<-MgCa.dat.Gbulloides[,20:26]
temp.sd<-MgCa.dat.Gbulloides[,62:68]
salt<-MgCa.dat.Gbulloides[,13:19]
salt.sd<-MgCa.dat.Gbulloides[,55:61]
Omega.s<-MgCa.dat.Gbulloides[,27:33]
Omega.s.sd<-MgCa.dat.Gbulloides[,69:75]
pH.s<-MgCa.dat.Gbulloides[,34:40]
pH.s.sd<-MgCa.dat.Gbulloides[,76:82]
Omega.d<-MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.cut
Omega.d.sd<-MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.sd.1
pH.d<-MgCa.dat.Gbulloides$GLODAP.pH.deep.cut
pH.d.sd<-MgCa.dat.Gbulloides$GLODAP.pH.deep.sd.1
MgCa<-MgCa.dat.Gbulloides$Mg.Ca
size<-MgCa.dat.Gbulloides$mean_size

#univariate
nvar<-7																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(salt.tmp, MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(10,-0.1))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(temp.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*exp(b1*x1)},s=c(1,0.1))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a*x1^b1},s=c(0.1,2))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(pH.s.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(10,-1))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(Omega.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,1))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(pH.d.tmp,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(10,1))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
	cal.val.nls(size,MgCa.tmp,prct.cal,1,f=function(x1,a,b1) {a+b1*x1},s=c(1,-0.01))
	tmp[nvar*(i-1)+7,]<-c(7,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Gbulloides.1<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Gbulloides.1$variables<-c("S","T","Omega.s","pH.s","Omega.d","pH.d","size")
Gbulloides.1$f.good<-f.good

#bivariate, keeping temperature
nvar<-6																						#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),nrow(temp),replace=T)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(temp.tmp,salt.tmp), MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1+b2*x2)},s=c(1,0.1,-0.05))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)*x2^b2},s=c(2,0.1,-0.1))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1,0.1,-1))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(0.5,0.1,0.5))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,pH.d.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(2,0.1,-0.1))
	tmp[nvar*(i-1)+5,]<-c(5,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,a,b1,b2) {a*exp(b1*x1)+b2*x2},s=c(1,0.1,0.001))
	tmp[nvar*(i-1)+6,]<-c(6,round(cal.valsummary[,1],3))
	}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Gbulloides.2<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Gbulloides.2$variables<-c("T,S","T,Omega.s","T,pH.s","T,Omega.d","T,pH.d","T,size")
Gbulloides.2$f.good<-f.good

#trivariate, keeping temperature and Omega.d
nvar<-4																					#number of variables
tmp<-matrix(,nrow=n.it*nvar,ncol=10)

for (i in 1:n.it) {
	#random sample of depth habitat. Assumes calibration/validation share the same depth
	rand<-sample(1:ncol(temp),1,)																
	temp.tmp<-temp[cbind(seq_along(rand),rand)]+rnorm(nrow(temp),0,temp.sd[cbind(seq_along(rand),rand)])
	salt.tmp<-salt[cbind(seq_along(rand),rand)]+rnorm(nrow(salt),0,salt.sd[cbind(seq_along(rand),rand)])
	Omega.s.tmp<-Omega.s[cbind(seq_along(rand),rand)]+rnorm(nrow(Omega.s),0,Omega.s.sd[cbind(seq_along(rand),rand)])
	pH.s.tmp<-pH.s[cbind(seq_along(rand),rand)]+rnorm(nrow(pH.s),0,pH.s.sd[cbind(seq_along(rand),rand)])
	Omega.d.tmp<-Omega.d+rnorm(length(Omega.d),0,Omega.d.sd)
	pH.d.tmp<-pH.d+rnorm(length(pH.d),0,pH.d.sd)
	MgCa.tmp<-MgCa+rnorm(length(MgCa),0,mean(MgCa)*MgCa.rsd)
	
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,salt.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1+b3*x3)+b2*x2},s=c(1,0.1,1,0.1))
	tmp[nvar*(i-1)+1,]<-c(1,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,Omega.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)*x3^b3+b2*x2},s=c(10,0.1,0.5,-1))
	tmp[nvar*(i-1)+2,]<-c(2,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,pH.s.tmp),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(10,0.1,0.5,-3))
	tmp[nvar*(i-1)+3,]<-c(3,round(cal.valsummary[,1],3))
	cal.val.nls(cbind(temp.tmp,Omega.d.tmp,size),MgCa.tmp,prct.cal,1,f=function(x1,x2,x3,a,b1,b2,b3) {a*exp(b1*x1)+b2*x2+b3*x3},s=c(10,0.1,0.5,-0.01))
	tmp[nvar*(i-1)+4,]<-c(4,round(cal.valsummary[,1],3))
}	

tmp.means<-matrix(,nrow=nvar,ncol=18)
tmp.good<-tmp[which(is.na(tmp[,2])==FALSE),]
f.good<-vector(,length=nvar)
for (i in 1:nvar) {
	dat.tmp<-tmp.good[which(tmp.good[,1]==i),]
	tmp.means[i,seq(1,17,2)]<-round(colMeans(dat.tmp,na.rm=T)[2:10],3)
	tmp.means[i,seq(2,18,2)]<-round(apply(dat.tmp,2,sd,na.rm=T)[2:10]/sqrt(nrow(dat.tmp))*2.576,3)
	f.good[i]<-nrow(dat.tmp)/n.it
}

colnames(tmp.means)<-c("r2","r2_99ci","BIC","BIC_99ci","int","int_99ci","coeff1","coeff1_99ci","coeff2","coeff2_99ci","coeff3","coeff3_99ci","coeff4","coeff4_99ci","RMSE","RMSE_99ci","CE","CE_99ci")
Gbulloides.3<-as.data.frame(tmp.means,stringsAsFactors=FALSE)
Gbulloides.3$variables<-c("T,Omega.d,S*","T,Omega.d,Omega.s*","T,Omega.d,pH.s*","T,Omega.d,size")
Gbulloides.3$f.good<-f.good

Gbulloides.all<-rbind(Gbulloides.1,Gbulloides.2,Gbulloides.3)


write.csv(NpachyL.all,"proxy_data/NpachyL_summary.csv")
write.csv(NpachyR.all,"proxy_data/NpachyR_summary.csv")
write.csv(Npachy.all,"proxy_data/Npachy_summary.csv")
write.csv(Gruber.all,"proxy_data/Gruber_summary.csv")
write.csv(Ginflata.all,"proxy_data/Ginflata_summary.csv")
write.csv(Gbulloides.all,"proxy_data/Gbulloides_summary.csv")