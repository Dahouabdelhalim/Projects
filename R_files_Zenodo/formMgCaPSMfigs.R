#Figure for foram Mg/Ca-environmental data calibration

setwd("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data");

library(ncdf4)
library(maps);library(mapproj); library(mapdata);

#coretop foram MgCa data
MgCa.dat.Npachy<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/Npachy.MgCa.env.dat.csv")
MgCa.dat.Gruber<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/Gruber.MgCa.env.dat.csv")
MgCa.dat.Ginflata<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/Ginflata.MgCa.env.dat.csv")
MgCa.dat.Gbulloides<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/Gbulloides.MgCa.env.dat.csv")

#subset Npachy
MgCa.dat.NpachyR<-MgCa.dat.Npachy[which(MgCa.dat.Npachy$Notes=="Dextral"),]
MgCa.dat.NpachyL<-MgCa.dat.Npachy[which(MgCa.dat.Npachy$Notes=="Sinstral"),]

#results from iterative regressions
CE.NpachyL<-read.csv("NpachyL_summary.csv")
CE.NpachyR<-read.csv("NpachyR_summary.csv")
CE.Npachy<-read.csv("Npachy_summary.csv")
CE.Gruber<-read.csv("Gruber_summary.csv")
CE.Ginflata<-read.csv("Ginflata_summary.csv")
CE.Gbulloides<-read.csv("Gbulloides_summary.csv")

#good results from iterative regressions
CE.NpachyL.good<-CE.NpachyL[which(CE.NpachyL$f.good>=0.9),]
CE.NpachyR.good<-CE.NpachyR[which(CE.NpachyR$f.good>=0.9),]
CE.Npachy.good<-CE.Npachy[which(CE.Npachy$f.good>=0.9),]
CE.Gruber.good<-CE.Gruber[which(CE.Gruber$f.good>=0.9),]
CE.Ginflata.good<-CE.Ginflata[which(CE.Ginflata$f.good>=0.9),]
CE.Gbulloides.good<-CE.Gbulloides[which(CE.Gbulloides$f.good>=0.9),]

#best models
NpachyL.final<-read.csv("NpachyL_alldata.cal.csv")
NpachyR.final<-read.csv("NpachyR_alldata.cal.csv")
Npachy.final<-read.csv("Npachy_alldata.cal.csv")
Gruber.final<-read.csv("Gruber_alldata.cal.csv")
Ginflata.final<-read.csv("Ginflata_alldata.cal.csv")
Gbulloides.final<-read.csv("Gbulloides_alldata.cal.csv")

#Figure 1 map of spatial distribution 
dev.new(width=12, height=6);
world<-map('worldHires',interior=FALSE,ylim=c(-88,88),xlim=c(-180,180))
points(MgCa.dat.Gruber$lon,MgCa.dat.Gruber$lat,cex=0.8,col="cornflowerblue")
points(MgCa.dat.Ginflata$lon,MgCa.dat.Ginflata$lat,cex=1.1,col="firebrick")
points(MgCa.dat.Gbulloides$lon,MgCa.dat.Gbulloides$lat,cex=1.4,col="gold3")
points(MgCa.dat.NpachyL$lon,MgCa.dat.NpachyL$lat,cex=0.5,col="palegreen3")
points(MgCa.dat.NpachyR$lon,MgCa.dat.NpachyR$lat,cex=0.5,col="grey50")
legend(30,88,legend=c("N. pachyderma","N. incompta","G. ruber","G. inflata","G. bulloides"), pch=c(1,1,1,1,1), col=c("palegreen3","grey50","cornflowerblue","firebrick","gold3"),bg="white")
dev.print(pdf,"MgCa_map.pdf");

#Figure 2, comparison of d18O and WOA T, generated in formMgCa_env.dat.R

#Figure 3, culture data, generated in formMgCa_culture

#Figure 4 summary of CE
dev.new(width=8, height=11);
par(mfrow=c(3,2))
par(mar=c(10,4,4,4))
plot(CE.NpachyL.good$CE,xaxt="n",xlab="",las=1,ylim=c(-0.2,0.8),ylab="CE",main="N. pachyderma",pch=16,cex=1.2)
segments(seq(1:nrow(CE.NpachyL)),CE.NpachyL.good$CE-CE.NpachyL.good$CE_99ci,seq(1:nrow(CE.NpachyL)),CE.NpachyL.good$CE+CE.NpachyL.good$CE_99ci)
points(which(CE.NpachyL.good$CE==max(CE.NpachyL.good$CE)),max(CE.NpachyL.good$CE),col="firebrick",cex=2)
axis(1,at=1:nrow(CE.NpachyL.good),labels=CE.NpachyL.good$variables,las=2,cex=0.5)
abline(0,0)

plot(CE.NpachyR.good$CE,xaxt="n",xlab="",las=1,ylim=c(-0.2,0.8),ylab="CE",main="N. incompta",pch=16,cex=1.2)
segments(seq(1:nrow(CE.NpachyR)),CE.NpachyR.good$CE-CE.NpachyR.good$CE_99ci,seq(1:nrow(CE.NpachyR)),CE.NpachyR.good$CE+CE.NpachyR.good$CE_99ci)
points(which(CE.NpachyR.good$CE==max(CE.NpachyR.good$CE)),max(CE.NpachyR.good$CE),col="firebrick",cex=2)
axis(1,at=1:nrow(CE.NpachyR.good),labels=CE.NpachyR.good$variables,las=2,cex=0.5)
abline(0,0)

#plot(CE.Npachy.good$CE,xaxt="n",xlab="",las=1,ylim=c(-0.2,0.8),ylab="CE",main="N. pachyderma & incompta",pch=16,cex=1.2)
#segments(seq(1:nrow(CE.Npachy)),CE.Npachy.good$CE-CE.Npachy.good$CE_99ci,seq(1:nrow(CE.Npachy)),CE.Npachy.good$CE+CE.Npachy.good$CE_99ci)
#points(which(CE.Npachy.good$CE==max(CE.Npachy.good$CE)),max(CE.Npachy.good$CE),col="firebrick",cex=2)
#axis(1,at=1:nrow(CE.Npachy.good),labels=CE.Npachy.good$variables,las=2,cex=0.5)
#abline(0,0)

plot(CE.Gruber.good$CE,xaxt="n",xlab="",las=1,ylim=c(-0.2,0.8),ylab="CE",main="G. ruber",pch=16,cex=1.2)
segments(seq(1:nrow(CE.Gruber)),CE.Gruber.good$CE-CE.Gruber.good$CE_99ci,seq(1:nrow(CE.Gruber)),CE.Gruber.good$CE+CE.Gruber.good$CE_99ci)
points(which(CE.Gruber.good$CE==max(CE.Gruber.good$CE)),max(CE.Gruber.good$CE),col="firebrick",cex=2)
axis(1,at=1:nrow(CE.Gruber.good),labels=CE.Gruber.good$variables,las=2,cex=0.5)
abline(0,0)

plot(CE.Ginflata.good$CE,xaxt="n",xlab="",las=1,ylim=c(-0.2,0.8),ylab="CE",main="G. inflata",pch=16,cex=1.2)
segments(seq(1:nrow(CE.Ginflata)),CE.Ginflata.good$CE-CE.Ginflata.good$CE_99ci,seq(1:nrow(CE.Ginflata)),CE.Ginflata.good$CE+CE.Ginflata.good$CE_99ci)
points(which(CE.Ginflata.good$CE==max(CE.Ginflata.good$CE)),max(CE.Ginflata.good$CE),col="firebrick",cex=2)
axis(1,at=1:nrow(CE.Ginflata.good),labels=CE.Ginflata.good$variables,las=2,cex=0.5)
abline(0,0)

plot(CE.Gbulloides.good$CE,xaxt="n",xlab="",las=1,ylim=c(-0.2,0.8),ylab="CE",main="G. bulloides",pch=16,cex=1.2)
segments(seq(1:nrow(CE.Gbulloides)),CE.Gbulloides.good$CE-CE.Gbulloides.good$CE_99ci,seq(1:nrow(CE.Gbulloides)),CE.Gbulloides.good$CE+CE.Gbulloides.good$CE_99ci)
points(which(CE.Gbulloides.good$CE==max(CE.Gbulloides.good$CE)),max(CE.Gbulloides.good$CE),col="firebrick",cex=2)
axis(1,at=1:nrow(CE.Gbulloides.good),labels=CE.Gbulloides.good$variables,las=2,cex=0.5)
abline(0,0)

dev.copy(pdf,"CEsummary.pdf")
dev.off()

#extra Figure summary of RMSE
dev.new(width=8, height=11);
par(mfrow=c(3,2))
par(mar=c(10,4,4,4))
plot(CE.NpachyL.good.good$RMSE,xaxt="n",xlab="",las=1,ylim=c(0.1,0.3),ylab="RMSE (mmol/mol)",main="N. pachyderma",pch=16,cex=1.2)
segments(seq(1:nrow(CE.NpachyL)),CE.NpachyL.good.good$RMSE-CE.NpachyL.good.good$RMSE_99ci,seq(1:nrow(CE.NpachyL)),CE.NpachyL.good.good$RMSE+CE.NpachyL.good$RMSE_99ci)
axis(1,at=1:nrow(CE.NpachyL),labels=CE.NpachyL$variables,las=2,cex=0.5)
abline(0,0)

plot(CE.NpachyR.good$RMSE,xaxt="n",xlab="",las=1,ylim=c(0.1,0.3),ylab="RMSE (mmol/mol)",main="N. incompta",pch=16,cex=1.2)
segments(seq(1:nrow(CE.NpachyR)),CE.NpachyR.good$RMSE-CE.NpachyR.good$RMSE_99ci,seq(1:nrow(CE.NpachyR)),CE.NpachyR.good$RMSE+CE.NpachyR.good$RMSE_99ci)
axis(1,at=1:nrow(CE.NpachyR),labels=CE.NpachyR$variables,las=2,cex=0.5)
abline(0,0)

#plot(CE.Npachy.good$RMSE,xaxt="n",xlab="",las=1,ylim=c(0.2,0.5),ylab="RMSE (mmol/mol)",main="N. pachyderma & incompta",pch=16,cex=1.2)
#segments(seq(1:nrow(CE.Npachy)),CE.Npachy.good$RMSE-CE.Npachy.good$RMSE_99ci,seq(1:nrow(CE.Npachy)),CE.Npachy.good$RMSE+CE.Npachy.good$RMSE_99ci)
#axis(1,at=1:nrow(CE.Npachy),labels=CE.Npachy$variables,las=2,cex=0.5)
#abline(0,0)

plot(CE.Gruber.good$RMSE,xaxt="n",xlab="",las=1,ylim=c(0.7,1.2),ylab="RMSE (mmol/mol)",main="G. ruber",pch=16,cex=1.2)
segments(seq(1:nrow(CE.Gruber)),CE.Gruber.good$RMSE-CE.Gruber.good$RMSE_99ci,seq(1:nrow(CE.Gruber)),CE.Gruber.good$RMSE+CE.Gruber.good$RMSE_99ci)
axis(1,at=1:nrow(CE.Gruber),labels=CE.Gruber$variables,las=2,cex=0.5)
abline(0,0)

plot(CE.Ginflata.good$RMSE,xaxt="n",xlab="",las=1,ylim=c(0,0.5),ylab="RMSE (mmol/mol)",main="G. inflata",pch=16,cex=1.2)
segments(seq(1:nrow(CE.Ginflata)),CE.Ginflata.good$RMSE-CE.Ginflata.good$RMSE_99ci,seq(1:nrow(CE.Ginflata)),CE.Ginflata.good$RMSE+CE.Ginflata.good$RMSE_99ci)
axis(1,at=1:nrow(CE.Ginflata),labels=CE.Ginflata$variables,las=2,cex=0.5)
abline(0,0)

plot(CE.Gbulloides.good$RMSE,xaxt="n",xlab="",las=1,ylim=c(0.8,2),ylab="RMSE (mmol/mol)",main="G. bulloides",pch=16,cex=1.2)
segments(seq(1:nrow(CE.Gbulloides)),CE.Gbulloides.good$RMSE-CE.Gbulloides.good$RMSE_99ci,seq(1:nrow(CE.Gbulloides)),CE.Gbulloides.good$RMSE+CE.Gbulloides.good$RMSE_99ci)
axis(1,at=1:nrow(CE.Gbulloides),labels=CE.Gbulloides$variables,las=2,cex=0.5)
abline(0,0)

dev.copy(pdf,"RMSEsummary.pdf")
dev.off()

#extra Figure summary of r. Note this is mislabelled as r2 in .good files
dev.new(width=8, height=11);
par(mfrow=c(3,2))
par(mar=c(10,4,4,4))
plot(CE.NpachyL.good$r2,xaxt="n",xlab="",las=1,ylim=c(0.1,0.9),ylab="r",main="N. pachyderma",pch=16,cex=1.2)
segments(seq(1:nrow(CE.NpachyL)),CE.NpachyL.good$r2-CE.NpachyL.good$r2_99ci,seq(1:nrow(CE.NpachyL)),CE.NpachyL.good$r2+CE.NpachyL.good$r2_99ci)
axis(1,at=1:nrow(CE.NpachyL),labels=CE.NpachyL$variables,las=2,cex=0.5)
abline(0,0)

plot(CE.NpachyR.good$r2,xaxt="n",xlab="",las=1,ylim=c(0.1,0.9),ylab="r",main="N. incompta",pch=16,cex=1.2)
segments(seq(1:nrow(CE.NpachyR)),CE.NpachyR.good$r2-CE.NpachyR.good$r2_99ci,seq(1:nrow(CE.NpachyR)),CE.NpachyR.good$r2+CE.NpachyR.good$r2_99ci)
axis(1,at=1:nrow(CE.NpachyR),labels=CE.NpachyR$variables,las=2,cex=0.5)
abline(0,0)

#plot(CE.Npachy.good$r2,xaxt="n",xlab="",las=1,ylim=c(0.1,0.9),ylab="r",main="N. pachyderma & incompta",pch=16,cex=1.2)
#segments(seq(1:nrow(CE.Npachy)),CE.Npachy.good$r2-CE.Npachy.good$r2_99ci,seq(1:nrow(CE.Npachy)),CE.Npachy.good$r2+CE.Npachy.good$r2_99ci)
#axis(1,at=1:nrow(CE.Npachy),labels=CE.Npachy$variables,las=2,cex=0.5)
#abline(0,0)

plot(CE.Gruber.good$r2,xaxt="n",xlab="",las=1,ylim=c(0.1,0.9),ylab="r",main="G. ruber",pch=16,cex=1.2)
segments(seq(1:nrow(CE.Gruber)),CE.Gruber.good$r2-CE.Gruber.good$r2_99ci,seq(1:nrow(CE.Gruber)),CE.Gruber.good$r2+CE.Gruber.good$r2_99ci)
axis(1,at=1:nrow(CE.Gruber),labels=CE.Gruber$variables,las=2,cex=0.5)
abline(0,0)

plot(CE.Ginflata.good$r2,xaxt="n",xlab="",las=1,ylim=c(0.1,0.9),ylab="r",main="G. inflata",pch=16,cex=1.2)
segments(seq(1:nrow(CE.Ginflata)),CE.Ginflata.good$r2-CE.Ginflata.good$r2_99ci,seq(1:nrow(CE.Ginflata)),CE.Ginflata.good$r2+CE.Ginflata.good$r2_99ci)
axis(1,at=1:nrow(CE.Ginflata),labels=CE.Ginflata$variables,las=2,cex=0.5)
abline(0,0)

plot(CE.Gbulloides.good$r2,xaxt="n",xlab="",las=1,ylim=c(0.1,0.9),ylab="r2",main="G. bulloides",pch=16,cex=1.2)
segments(seq(1:nrow(CE.Gbulloides)),CE.Gbulloides.good$r2-CE.Gbulloides.good$r2_99ci,seq(1:nrow(CE.Gbulloides)),CE.Gbulloides.good$r2+CE.Gbulloides.good$r2_99ci)
axis(1,at=1:nrow(CE.Gbulloides),labels=CE.Gbulloides$variables,las=2,cex=0.5)
abline(0,0)

dev.copy(pdf,"r_summary.pdf")
dev.off()

#Figure 5, comparison of predicted vs. observed. generated in formMgCa_all.R 

#Figure 6 (formerly Figure 9), comparison with Ocean2K originally published temperatures
DT<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/O2K_DT.csv")
DT.stda.bin<-matrix(,nrow=1,ncol=4)
dev.new(width=8, height=5);
par(mfrow=c(1,2))
plot(DT$temp[which(DT$Species=="G.ruber")],DT$T_uni[which(DT$Species=="G.ruber")],xlab="T (ºC, original)",ylab="T (ºC, this study)", xlim=c(0,32),ylim=c(0,32),las=1,col="cornflowerblue",main="T-only",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="G.bulloides")],DT$T_uni[which(DT$Species=="G.bulloides")], col="gold3",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="G.inflata")],DT$T_uni[which(DT$Species=="G.inflata")], col="firebrick",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="N.pachyderma")],DT$T_uni[which(DT$Species=="N.pachyderma")], col="palegreen3",pch=20,cex=0.4)
abline(0,1)
pal<-c("gold3","firebrick","cornflowerblue","palegreen3")

plot(DT$temp[which(DT$Species=="G.ruber")],DT$T_mult[which(DT$Species=="G.ruber")],xlab="T (ºC, original)",ylab="T (ºC, this study)", xlim=c(0,32),ylim=c(0,32),las=1,col="cornflowerblue",main="multivariate",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="G.bulloides")],DT$T_mult[which(DT$Species=="G.bulloides")], col="gold3",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="G.inflata")],DT$T_mult[which(DT$Species=="G.inflata")], col="firebrick",pch=20,cex=0.4)
points(DT$temp[which(DT$Species=="N.pachyderma")],DT$T_mult[which(DT$Species=="N.pachyderma")], col="palegreen3",pch=20,cex=0.4)
abline(0,1)
pal<-c("gold3","firebrick","cornflowerblue","palegreen3")
legend("topleft",legend=c("G. bulloides","G. inflata","G. ruber","N. pachyderam & incompta"),text.col=pal[seq(1,4)],bty="n",cex=0.8)

#Figure 7, (formerly figure 10) redo of McGregor bins
dev.new(width=6, height=6);
plot(0,type="n",xlab="mean change in T (ºC, this study-original)",ylab="change in std. dev. (ºC, this study-original)",xlim=c(-5,5),ylim=c(-1,2),las=1)
abline(h=0,col="grey50");abline(v=0,col="grey50")
for (i in 1:length(unique(DT$index))) {
	tmp<-DT[which(DT$index==unique(DT$index)[i]),]
	var<-sd(tmp$T_mult)-sd(tmp$temp)
	Dtemp<-mean(tmp$T_mult-tmp$temp)
	text(Dtemp,var,labels=tmp$index,col=pal[as.numeric(tmp$Species[1])])
	temp.sa<-(tmp$temp-mean(tmp$temp))/sd(tmp$temp)											#original T standard anomaly
	T_mult.sa<-(tmp$T_mult-mean(tmp$T_mult))/sd(tmp$T_mult)							#this study T standard anomaly					
	temp.sa.bin<-as.numeric(tapply(temp.sa,cut(tmp$year,seq(0,2000,by=200)),mean))			#bin every 200 years
	T_mult.sa.bin<-as.numeric(tapply(T_mult.sa,cut(tmp$year,seq(0,2000,by=200)),mean))
	tmp1<-cbind(rep(tmp$index[1],10),seq(100,1900,by=200),temp.sa.bin,T_mult.sa.bin)
	DT.stda.bin<-rbind(DT.stda.bin,tmp1)
	}
legend("bottomleft",legend=c("G. bulloides","G. inflata","G. ruber","N. pachyderam & incompta"),text.col=pal[seq(1,4)],bty="n")
DT.stda.bin<-DT.stda.bin[-1,]
dev.new(width=7, height=6);
par(mfrow=c(2,1),mar=c(5,5,1,1))
boxplot(DT.stda.bin[,3]~DT.stda.bin[,2],names=seq(100,1900,200),las=1,xlab="",ylab="standardized T anomaly",ylim=c(-2.5,2), main="original")
boxplot(DT.stda.bin[,4]~DT.stda.bin[,2],names=seq(100,1900,200),las=1,xlab="center year of bin (A.D.)",ylab="standardized T anomaly",ylim=c(-2.5,2),main="this study")

