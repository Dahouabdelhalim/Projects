rm(list = ls())
library(graphics); library(compositions); library(mgcv);
#
setwd("~/Dropbox/__CoDa_JMP")

filename<-"Sanitation.csv"; 
data<-read.csv(file=filename, header=TRUE, sep = ";", dec = ","); data0S<-data; levels(data0S$country)
filename<-"Water.csv"; 
data<-read.csv(file=filename, header=TRUE, sep = ";", dec = ","); data0W<-data; levels(data0W$country)

xlimgraf <-c(1990,2015); 
deltazero <- 0.5E-5

plotpdf<-TRUE
#plotpdf<-FALSE

################################
# FIGURA_1 & TABLE 1
################################
if(plotpdf){pdf(file=paste("plots_FIG1.pdf"),width=9, height=9)}

ylimgraf <-c(0.3,0.8);   
ylimgraf2<-c(0.5,1.); 
ylimgraf3<-c(0.1,0.6);   
ylimgraf4<-c(0.,0.5); 

grafico<-function(xx,year,ylimgraf,xlimgraf){
  yearini <- year[1]
  year0 <- year - yearini
  
  dataexp<-xx[,1]
  regr_lin<-lm(dataexp ~ year0)

  yyear<-c((-2):(max(year0))) 
  yyyear<-c((max(year0)+1):(max(year0)+2)) 
  yyyyear<-c(yyear,yyyear)
  
  A<-array(cbind(1.,yyear),dim = c(length(yyear),2))
  B<-array(regr_lin$coefficients,dim = c(2,1))
  yy <- A %*% B
  C<-array(cbind(1.,yyyear),dim = c(length(yyyear),2))
  D<-array(regr_lin$coefficients,dim = c(2,1))
  zz <- C %*% D 
  yyy<-c(yy,zz)
  
  plot(year0+yearini,dataexp, lwd=c(3), col="black",xlim=xlimgraf,ylim=ylimgraf,xlab="",ylab="", cex.axis=1.5)
  title(xlab="YEAR",ylab="Proportion",cex.lab=1.5)
  lines(yyyyear+yearini,yyy,type="l", lwd=c(3), col="black")
  
  gam <- gam(dataexp ~ s(year0,k=4), family=gaussian); gam 
  gam0.pred <- predict.gam(gam, data.frame(year0=yyyyear))
  lines(yyyyear+yearini,gam0.pred,type="l",lwd=c(3),col="blue")
  
  dd<-unclass(ilr(xx))
  comp_linmodel <- lm(dd ~ year0)
  
  Ac<-array(cbind(1.,yyear),dim = c(length(yyear),2))
  Bc<-array(comp_linmodel$coefficients,dim = c(2,(dim(dd)[2])))
  yyc <- Ac %*% Bc
  Cc<-array(cbind(1.,yyyear),dim = c(length(yyyear),2))
  Dc<-array(comp_linmodel$coefficients,dim = c(2,(dim(dd)[2])))
  zzc <- Cc %*% Dc 
  yyyc<-rbind(yyc,zzc)
  zzz0<-ilrInv(yyyc)
  lines(yyyyear+yearini,zzz0[,1],type="l",lwd=c(3),col="black",lty=3)
  
  gam <- gam(dd ~ s(year0,k=4), family=gaussian)
  gam.pred <- predict.gam(gam, data.frame(year0=yyyyear))
  zzz0<-ilrInv(array(gam.pred,dim = c(length(gam.pred),1)))
  lines(yyyyear+yearini,zzz0[,1],type="l",lwd=c(3),col="blue",lty=3)
}
errores<-function(xx,year){
  yearini <- year[1]
  year0 <- year - yearini

  dataexp<-xx[,1]
  regr_lin<-lm(dataexp ~ year0)
  A<-array(cbind(1.,year0),dim = c(length(year0),2))
  B<-array(regr_lin$coefficients,dim = c(2,1))
  yy <- A %*% B

  cc1<-1./length(dataexp)
  ccc<-mean(dataexp)
  cc2<-1./sum((dataexp-ccc)^2)
  aux11<-summary(lm(yy~dataexp))$r.squared
  aux12<-sqrt(sum((yy-dataexp)^2)*cc1)
  aux13<-1. - sum((yy-dataexp)^2)*cc2

  gamx <- gam(dataexp ~ s(year0,k=4), family=gaussian);
  aux31<-(summary(lm(gamx$fitted.values~dataexp))$r.squared)
  aux32<-sqrt(sum((gamx$fitted.values-dataexp)^2)*cc1)
  aux33<-1. - sum((gamx$fitted.values-dataexp)^2)*cc2

  dd<-unclass(ilr(xx))
  comp_linmodel <- lm(dd ~ year0)
  Ac<-array(cbind(1.,year0),dim = c(length(year),2))
  Bc<-array(comp_linmodel$coefficients,dim = c(2,(dim(dd)[2])))
  yyc <- Ac %*% Bc
  zzz0<-ilrInv(yyc)

  aux21<-summary(lm(zzz0[,1]~dataexp))$r.squared
  aux22<-sqrt(sum((zzz0[,1]-dataexp)^2)*cc1)
  aux23<-1. - sum((zzz0[,1]-dataexp)^2)*cc2

  gamx <- gam(dd ~ s(year0,k=4), family=gaussian)
  zzz0<-ilrInv(array(gamx$fitted.values,dim = c(length(gamx$fitted.values),1)))
  aux41<-summary(lm(zzz0[,1]~dataexp))$r.squared
  aux42<-sqrt(sum((zzz0[,1]-dataexp)^2)*cc1)
  aux43<-1. - sum((zzz0[,1]-dataexp)^2)*cc2
  
  resul<-matrix(c(aux11,aux12,aux13,aux21,aux22,aux23,aux31,aux32,aux33,aux41,aux42,aux43),c(4,3),byrow=TRUE)
  resul<-data.frame(resul)
  names(resul)<-c("r2","RMSE", "NSE")
  return(resul)
  }

data0 <- data0W; iipais<-1; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,5]);  xx<-data[aux,c(5:6)]; year<-data[aux,4]; #RURAL
grafico(xx,year,ylimgraf,xlimgraf)
title(paste("a.1) Improved WATER Rural - ",pais), cex.main=1.5)
legend(2007,.45, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))
paste("a.1) Improved WATER Rural - ",pais)
errores(xx,year)

data0 <- data0W; iipais<-5; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,5]);  xx<-data[aux,c(5:6)]; year<-data[aux,4]; #RURAL
grafico(xx,year,ylimgraf2,xlimgraf)
title(paste("a.2) Improved WATER Rural - ",pais), cex.main=1.5)
legend(2007,.65, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))
paste("a.2) Improved WATER Rural - ",pais)
errores(xx,year)

data0 <- data0W; iipais<-1; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,7]);  xx<-data[aux,c(7:8)]; year<-data[aux,4]; #URBAN
grafico(xx,year,ylimgraf2,xlimgraf)
title(paste("b.1) Improved WATER Urban - ",pais), cex.main=1.5)
legend(2007,.65, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))
paste("b.1) Improved WATER Urban - ",pais)
errores(xx,year)

data0 <- data0W; iipais<-5; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,7]);  xx<-data[aux,c(7:8)]; year<-data[aux,4]; #URBAN
grafico(xx,year,ylimgraf2,xlimgraf)
title(paste("b.2) Improved WATER Urban - ",pais), cex.main=1.5)
legend(2007,.65, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))
paste("b.2) Improved WATER Urban - ",pais)
errores(xx,year)

data0 <- data0S; iipais<-4; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,5]);  xx<-data[aux,c(5:6)]; year<-data[aux,4]; #RURAL
grafico(xx,year,ylimgraf4,xlimgraf)
title(paste("c.1) Improved SANITATION Rural - ",pais), cex.main=1.5)
legend(1992,.50, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))
paste("c.1) Improved SANITATION Rural - ",pais)
errores(xx,year)

data0 <- data0W; iipais<-4; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,5]);  xx<-data[aux,c(5:6)]; year<-data[aux,4]; #RURAL
grafico(xx,year,ylimgraf3,xlimgraf)
title(paste("c.2) Improved WATER Rural - ",pais), cex.main=1.5)
legend(1992,.60, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))
paste("c.2) Improved WATER Rural - ",pais)
errores(xx,year)

data0 <- data0W; iipais<-4; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,7]);  xx<-data[aux,c(7:8)]; year<-data[aux,4]; #URBAN
grafico(xx,year,ylimgraf2,xlimgraf)
title(paste("d.1) Improved WATER Urban - ",pais), cex.main=1.5)
legend(2007,.65, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))
paste("d.1) Improved WATER Urban - ",pais)
errores(xx,year)

data0 <- data0S; iipais<-5; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,5]);  xx<-data[aux,c(5:6)]; year<-data[aux,4]; #RURAL
grafico(xx,year,ylimgraf4,xlimgraf)
title(paste("d.2) Improved SANITATION Rural - ",pais), cex.main=1.5)
legend(1992,.50, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))
paste("d.2) Improved SANITATION Rural - ",pais)
errores(xx,year)

if(plotpdf){dev.off()}

################################
# FIGURA 2
################################
if(plotpdf){pdf(file=paste("plots_FIG2.pdf"),width=9, height=9)}

ylimgraf <-c(0.,0.2); 
ylimgraf2<-c(0.5,1.); 
ylimgraf3<-c(0.1,0.6); 
ylimgraf4<-c(0.,0.1); 
ylimgraf5<-c(0.4,0.8); 
ylimgraf6<-c(0.1,0.4); 

kk<-4
graficoVV<-function(xx,year,ylimgraf,xlimgraf,VV,kk){
  dataexp<-xx[,1]
  yearini <- year[1]
  tt <- year - yearini
  ttt<-c((-2):(max(tt)+2)) 
  plot(tt+yearini,dataexp, lwd=c(3), col="black",xlim=xlimgraf,ylim=ylimgraf,xlab="",ylab="", cex.axis=2.0)
  ols  <- lm(dataexp ~ tt); yyy <- predict(ols, data.frame(tt=ttt))
  gam  <- gam(dataexp ~ s(tt,k=kk), family=gaussian) 
  gam0.pred <- predict.gam(gam, data.frame(tt=ttt))
  lines(ttt+yearini,yyy,type="l",lwd=c(3),col="black")
  lines(ttt+yearini,gam0.pred,type="l",lwd=c(3),col="blue")
  dd   <- unclass(ilr(xx,VV)); 
  ols  <- lm(dd ~ tt);  yyyc <- predict(ols, data.frame(tt=ttt)); 
  zzz0 <- ilrInv(yyyc,VV)
  lines(ttt+yearini,zzz0[,1],type="l",lwd=c(3),col="black",lty=3)
  gam1 <- gam(dd[,1] ~ s(tt,k=kk), family=gaussian);     gam1.pred <- predict.gam(gam1, data.frame(tt=ttt))
  gam2 <- gam(dd[,2] ~ s(tt,k=kk), family=gaussian);     gam2.pred <- predict.gam(gam2, data.frame(tt=ttt))
  gam3 <- gam(dd[,3] ~ s(tt,k=kk), family=gaussian);     gam3.pred <- predict.gam(gam3, data.frame(tt=ttt))
  zzz1<-ilrInv(array(cbind(gam1.pred,gam2.pred,gam3.pred),dim = c(length(gam1.pred),3)),VV)
  lines(ttt+yearini,zzz1[,1],type="l",lwd=c(3),col="blue",lty=3)
}

data0 <- data0W; 
iipais<-4; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,9]);  xx<-data[aux,c(9:12)]; year<-data[aux,4]; #RURAL
xx[xx[,1]==0.,1]<-deltazero;
signs <- rbind (c( 1, 1, -1, -1),c(1, -1, 0, 0),c(0, 0, 1, -1)); VV=gsi.buildilrBase(t(signs))
graficoVV(xx,year,ylimgraf4,xlimgraf,VV,kk)
title(paste("a.1) PIPED WATER Rural - ",pais), cex.main=1.5)
legend(1992,.10, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=7, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

yy<-cbind(xx[,2:4],xx[,1]);names(yy)[4]<-names(xx)[1];xx<-yy;
yy<-rbind(VV[2:4,],VV[1,]);VV<-yy;
graficoVV(xx,year,ylimgraf3,xlimgraf,VV,kk)
title(paste("a.2) OTHER IMPROVED Rural - ",pais), cex.main=1.5)
legend(1992,.60, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

yy<-cbind(xx[,2:4],xx[,1]);names(yy)[4]<-names(xx)[1];xx<-yy;
yy<-rbind(VV[2:4,],VV[1,]);VV<-yy;
graficoVV(xx,year,ylimgraf3,xlimgraf,VV,kk)
title(paste("a.4) SURFACE WATER Rural - ",pais), cex.main=1.5)
legend(1992,.25, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

yy<-cbind(xx[,2:4],xx[,1]);names(yy)[4]<-names(xx)[1];xx<-yy;
yy<-rbind(VV[2:4,],VV[1,]);VV<-yy;
graficoVV(xx,year,ylimgraf3,xlimgraf,VV,kk)
title(paste("a.3) OTHER UNIMPROVED Rural - ",pais), cex.main=1.5)
legend(1992,.25, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

#
data0 <- data0S; 
iipais<-5; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,9]);  xx<-data[aux,c(9:12)]; year<-data[aux,4]; #RURAL
signs <- rbind (c( 1, 1, -1, -1),c(1, -1, 0, 0),c(0, 0, 1, -1)); VV=gsi.buildilrBase(t(signs))
graficoVV(xx,year,ylimgraf,xlimgraf,VV,kk)
title(paste("b.1) IMPROVED SHARED Rural - ",pais), cex.main=1.5)
legend(1992,.05, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

yy<-cbind(xx[,2:4],xx[,1]);names(yy)[4]<-names(xx)[1];xx<-yy;
yy<-rbind(VV[2:4,],VV[1,]);VV<-yy;
graficoVV(xx,year,ylimgraf,xlimgraf,VV,kk)
title(paste("b.2) IMPROVED NON-SHARED Rural - ",pais), cex.main=1.5)
legend(2007,.2, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

yy<-cbind(xx[,2:4],xx[,1]);names(yy)[4]<-names(xx)[1];xx<-yy;
yy<-rbind(VV[2:4,],VV[1,]);VV<-yy;
graficoVV(xx,year,ylimgraf5,xlimgraf,VV,kk)
title(paste("b.4) OPEN DEFECATION Rural - ",pais), cex.main=1.5)
legend(2007,.8, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

yy<-cbind(xx[,2:4],xx[,1]);names(yy)[4]<-names(xx)[1];xx<-yy;
yy<-rbind(VV[2:4,],VV[1,]);VV<-yy;
graficoVV(xx,year,ylimgraf6,xlimgraf,VV,kk)
title(paste("b.3) OTHER UNIMPROVED Rural - ",pais), cex.main=1.5)
legend(1992,.4, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

#
if(plotpdf){dev.off()}

################################
# FIGURA 3
################################
if(plotpdf){pdf(file=paste("plots_FIG3.pdf"),width=9, height=9)}

sumagam4_ols<-function(xx,tt,ttt){
  ols  <- lm(xx[,1] ~ tt)
  olsp <- predict(ols, data.frame(tt=ttt))
  ols  <- lm(xx[,2] ~ tt)
  olsp2 <- predict(ols, data.frame(tt=ttt))
  ols  <- lm(xx[,3] ~ tt)
  olsp3 <- predict(ols, data.frame(tt=ttt))
  ols  <- lm(xx[,4] ~ tt)
  olsp4 <- predict(ols, data.frame(tt=ttt))
  return(olsp+olsp2+olsp3+olsp4)
}
sumagam4_cs <-function(xx,tt,ttt){
  dataexp<-xx[,1]
  gam <- gam(dataexp ~ s(tt,bs="cr"), family=gaussian); gam 
  gam1.pred <- predict.gam(gam, data.frame(tt=ttt))
  dataexp<-xx[,2]
  gam2 <- gam(dataexp ~ s(tt,bs="cr"), family=gaussian); gam2 
  gam2.pred <- predict.gam(gam2, data.frame(tt=ttt))
  dataexp<-xx[,3]
  gam3 <- gam(dataexp ~ s(tt,bs="cr"), family=gaussian); gam3 
  gam3.pred <- predict.gam(gam3, data.frame(tt=ttt))
  dataexp<-xx[,4]
  gam4 <- gam(dataexp ~ s(tt,bs="cr"), family=gaussian); gam4 
  gam4.pred <- predict.gam(gam4, data.frame(tt=ttt))
  return(gam1.pred+gam2.pred+gam3.pred+gam4.pred)
}
sumagam4_sk <-function(xx,tt,ttt,kk){
  dataexp<-xx[,1]
  gam <- gam(dataexp ~ s(tt,k=kk), family=gaussian); gam 
  gam1.pred <- predict.gam(gam, data.frame(tt=ttt))
  dataexp<-xx[,2]
  gam2 <- gam(dataexp ~ s(tt,k=kk), family=gaussian); gam2 
  gam2.pred <- predict.gam(gam2, data.frame(tt=ttt))
  dataexp<-xx[,3]
  gam3 <- gam(dataexp ~ s(tt,k=kk), family=gaussian); gam3 
  gam3.pred <- predict.gam(gam3, data.frame(tt=ttt))
  dataexp<-xx[,4]
  gam4 <- gam(dataexp ~ s(tt,k=kk), family=gaussian); gam4 
  gam4.pred <- predict.gam(gam4, data.frame(tt=ttt))
  return(gam1.pred+gam2.pred+gam3.pred+gam4.pred)
}
sumagam3_cs <-function(xx,tt,ttt){
  signs <- rbind (c( 1, 1, -1, -1),c(1, -1, 0, 0),c(0, 0, 1, -1)); VV=gsi.buildilrBase(t(signs))
  dd<-unclass(ilr(xx,VV))
  dataexp<-dd[,1]
  gam1 <- gam(dataexp ~ s(tt,bs="cr"), family=gaussian); gam1 
  gam1.pred <- predict.gam(gam1, data.frame(tt=ttt))
  dataexp<-dd[,2]
  gam2 <- gam(dataexp ~ s(tt,bs="cr"), family=gaussian); gam2 
  gam2.pred <- predict.gam(gam2, data.frame(tt=ttt))
  dataexp<-dd[,3]
  gam3 <- gam(dataexp ~ s(tt,bs="cr"), family=gaussian); gam3 
  gam3.pred <- predict.gam(gam3, data.frame(tt=ttt))
  aux<-ilrInv(array(cbind(gam1.pred,gam2.pred,gam3.pred),dim = c(length(gam1.pred),3)),VV)
  return(aux[,1]+aux[,2]+aux[,3]+aux[,4])
}
sumagam3_sk <-function(xx,tt,ttt,kk){
  signs <- rbind (c( 1, 1, -1, -1),c(1, -1, 0, 0),c(0, 0, 1, -1)); VV=gsi.buildilrBase(t(signs))
  dd<-unclass(ilr(xx,VV))
  dataexp<-dd[,1]
  gam1 <- gam(dataexp ~ s(tt,k=kk), family=gaussian); gam1 
  gam1.pred <- predict.gam(gam1, data.frame(tt=ttt))
  dataexp<-dd[,2]
  gam2 <- gam(dataexp ~ s(tt,k=kk), family=gaussian); gam2 
  gam2.pred <- predict.gam(gam2, data.frame(tt=ttt))
  dataexp<-dd[,3]
  gam3 <- gam(dataexp ~ s(tt,k=kk), family=gaussian); gam3 
  gam3.pred <- predict.gam(gam3, data.frame(tt=ttt))
  aux<-ilrInv(array(cbind(gam1.pred,gam2.pred,gam3.pred),dim = c(length(gam1.pred),3)),VV)
  return(aux[,1]+aux[,2]+aux[,3]+aux[,4])
}

data0 <- data0S; iipais<-5; # SANITATION 
#data0 <- data0W; iipais<-5; # WATER 

pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,9]);  xx<-data[aux,c(9:12)]; year<-data[aux,4]; #RURAL
#aux<-!is.na(data[,13]);  xx<-data[aux,c(13:16)]; year<-data[aux,4]; #URBAN

yearini <- year[1]
tt <- year - yearini
ttt<-c((-2):(2020-yearini)) 

ylimgraf<-c(0.94,1.06);

kk<-4; aux<- sumagam4_sk(xx,tt,ttt,kk)
plot(ttt+yearini,aux,type="l", lwd=c(3), col="blue",xlim=xlimgraf,ylim=ylimgraf,xlab="",ylab="", cex.axis=1.5)
title(xlab="YEAR",ylab="SUM PROP.",cex.lab=1.5)
kk<-4; aux<- sumagam3_sk(xx,tt,ttt,kk)
lines(ttt+yearini,aux,type="l",lwd=c(3),col="blue",lty=3)
legend(2005,.97, c("GAM x 4","GAM(ilr) x 3"),text.width=5.5, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3), lwd=c(3,3),col=c("blue","blue"))
title(paste("GAM tps-4 - San. Rural - ",pais), cex.main=1.5)

kk<-6; aux<- sumagam4_sk(xx,tt,ttt,kk)
plot(ttt+yearini,aux,type="l", lwd=c(3), col="blue",xlim=xlimgraf,ylim=ylimgraf,xlab="",ylab="", cex.axis=1.5)
title(xlab="YEAR",ylab="SUM PROP.",cex.lab=1.5)
kk<-6; aux<- sumagam3_sk(xx,tt,ttt,kk)
lines(ttt+yearini,aux,type="l",lwd=c(3),col="blue",lty=3)
legend(2005,.97, c("GAM x 4","GAM(ilr) x 3"),text.width=5.5, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3), lwd=c(3,3),col=c("blue","blue"))
title(paste("GAM tps-6 - San. Rural - ",pais), cex.main=1.5)

kk<-8; aux<- sumagam4_sk(xx,tt,ttt,kk)
plot(ttt+yearini,aux,type="l", lwd=c(3), col="blue",xlim=xlimgraf,ylim=ylimgraf,xlab="",ylab="", cex.axis=1.5)
title(xlab="YEAR",ylab="SUM PROP.",cex.lab=1.5)
kk<-8; aux<- sumagam3_sk(xx,tt,ttt,kk)
lines(ttt+yearini,aux,type="l",lwd=c(3),col="blue",lty=3)
legend(2005,.97, c("GAM x 4","GAM(ilr) x 3"),text.width=5.5, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3), lwd=c(3,3),col=c("blue","blue"))
title(paste("GAM tps-8 - San. Rural - ",pais), cex.main=1.5)

aux<- sumagam4_cs(xx,tt,ttt)
plot(ttt+yearini,aux,type="l", lwd=c(3), col="blue",xlim=xlimgraf,ylim=ylimgraf,xlab="",ylab="", cex.axis=1.5)
title(xlab="YEAR",ylab="SUM PROP.",cex.lab=1.5)
aux<- sumagam3_cs(xx,tt,ttt)
lines(ttt+yearini,aux,type="l",lwd=c(3),col="blue",lty=3)
legend(2005,.97, c("GAM x 4","GAM(ilr) x 3"),text.width=5.5, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3), lwd=c(3,3),col=c("blue","blue"))
title(paste("GAM cs - San. Rural - ",pais), cex.main=1.5)

if(plotpdf){dev.off()}

################################
# FIGURA 4
################################
#kk<-4 GAM dof 
kk<-8

if(plotpdf){pdf(file=paste("plots_FIG4.pdf"),width=9, height=9)}

data0 <- data0S; iipais<-5; pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; 
aux<-!is.na(data[,13]);  xx<-data[aux,c(13:16)]; year<-data[aux,4]; #URBAN
signs <- rbind (c( 1, 1, -1, -1),c(1, -1, 0, 0),c(0, 0, 1, -1)); VV=gsi.buildilrBase(t(signs))
graficoVV(xx,year,ylimgraf6,xlimgraf,VV,kk)
title(paste("b) SHARED SANITATION Urban - ",pais), cex.main=1.5)
legend(1992,.4, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

yy<-cbind(xx[,2:4],xx[,1]);names(yy)[4]<-names(xx)[1];xx<-yy;
yy<-rbind(VV[2:4,],VV[1,]);VV<-yy;
graficoVV(xx,year,ylimgraf6,xlimgraf,VV,kk)
title(paste("a) IMPROVED SANITATION Urban - ",pais), cex.main=1.5)
legend(1992,.4, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

yy<-cbind(xx[,2:4],xx[,1]);names(yy)[4]<-names(xx)[1];xx<-yy;
yy<-rbind(VV[2:4,],VV[1,]);VV<-yy;
graficoVV(xx,year,ylimgraf6,xlimgraf,VV,kk)
title(paste("d) OPEN DEFECATION Urban - ",pais), cex.main=1.5)
legend(1992,.18, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

yy<-cbind(xx[,2:4],xx[,1]);names(yy)[4]<-names(xx)[1];xx<-yy;
yy<-rbind(VV[2:4,],VV[1,]);VV<-yy;
graficoVV(xx,year,ylimgraf3,xlimgraf,VV,kk)
title(paste("c) OTHER UNIMPROVED Urban - ",pais), cex.main=1.5)
legend(1992,.25, c("OLS","OLS(ilr)","GAM","GAM(ilr)"),text.width=4, seg.len=2.5,y.intersp=1.5,cex=1.3,lty=c(1,3,1,3), lwd=c(3,3,3,3),col=c("black","black","blue","blue"))

if(plotpdf){dev.off()}

################################
# TABLE 2
################################
{
  data0 <- data0S; iipais<-2; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,13]);  xx<-data[aux,c(13:16)]; year<-data[aux,4]; 
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2)
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  c(NA,NA) -> aux4 ; 
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("SAN URB - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-aux
  data0 <- data0S; iipais<-4; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,13]);  xx<-data[aux,c(13:16)]; year<-data[aux,4]; 
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2)
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  c(NA,NA) -> aux4 ; 
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("SAN URB - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux)
  data0 <- data0S; iipais<-5; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,13]);  xx<-data[aux,c(13:16)]; year<-data[aux,4]; 
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2)
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  sumagam4_cs(xx,tt,ttt2) -> aux4 ;
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("SAN URB - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux) 
}
{
  data0 <- data0S; iipais<-2; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,9]);  xx<-data[aux,c(9:12)]; year<-data[aux,4]; 
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2)
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  c(NA,NA) -> aux4 ; 
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("SAN RUR - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux) 
  data0 <- data0S; iipais<-4; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,9]);  xx<-data[aux,c(9:12)]; year<-data[aux,4];
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2);
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  c(NA,NA) -> aux4 ;
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("SAN RUR - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux) 
  data0 <- data0S; iipais<-5; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,9]);  xx<-data[aux,c(9:12)]; year<-data[aux,4]; 
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2)
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  sumagam4_cs(xx,tt,ttt2) -> aux4 ;
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("SAN RUR - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux) 
}
{
  data0 <- data0W; iipais<-1; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,13]); xx<-data[aux,c(13:16)]; year<-data[aux,4]; 
  xx[xx[,3]==0.,3]<-deltazero;
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2); 
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  sumagam4_cs(xx,tt,ttt2) -> aux4 ;
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("WAT URB - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux)  
  data0 <- data0W; iipais<-4; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,13]); xx<-data[aux,c(13:16)]; year<-data[aux,4]; 
  xx[xx[,3]==0.,3]<-deltazero;
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2); 
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  c(NA,NA) -> aux4 ;
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("WAT URB - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux)    
  data0 <- data0W; iipais<-5; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,13]); xx<-data[aux,c(13:16)]; year<-data[aux,4]; 
  xx[xx[,3]==0.,3]<-deltazero;
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2); 
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  sumagam4_cs(xx,tt,ttt2) -> aux4 ;
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("WAT URB - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux)  
}
{
  data0 <- data0W; iipais<-1; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,9]); xx<-data[aux,c(9:12)]; year<-data[aux,4]; 
  xx[xx[,1]==0.,1]<-deltazero;
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2);
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  sumagam4_cs(xx,tt,ttt2) -> aux4 ;
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("WAT RUR - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux)   
  data0 <- data0W; iipais<-4; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,9]); xx<-data[aux,c(9:12)]; year<-data[aux,4]; 
  xx[xx[,1]==0.,1]<-deltazero;
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2);
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  sumagam4_cs(xx,tt,ttt2) -> aux4 ;
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("WAT RUR - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux)  
  data0 <- data0W; iipais<-5; 
  pais<-levels(data0$country)[iipais]; data<-data0[data0$country==pais,]; aux<-!is.na(data[,9]); xx<-data[aux,c(9:12)]; year<-data[aux,4]; 
  xx[xx[,1]==0.,1]<-deltazero;
  yearini <- year[1]; tt <- year - yearini; ttt2<-c(-2, max(tt) + 2);
  kk<-4; sumagam4_sk(xx,tt,ttt2,kk) -> aux1 ; 
  kk<-6; sumagam4_sk(xx,tt,ttt2,kk) -> aux2 ;
  kk<-8; sumagam4_sk(xx,tt,ttt2,kk) -> aux3 ; 
  sumagam4_cs(xx,tt,ttt2) -> aux4 ;
  aux<-c(aux1,aux2,aux3,aux4);aux<-data.frame(t(aux));row.names(aux)<-paste("WAT RUR - ",pais); colnames(aux)<-c("1g0","1g1","2g0","2g1","3g0","3g1","4g0","4g1"); auxx<-rbind(auxx,aux)  
}
zzz<-auxx-1.;zzz