

setwd("C:\\\\Users\\\\vkara\\\\Dropbox\\\\Projects\\\\quaggaGL\\\\BYONI IMPACTS")
#remotes::install_github("gavinsimpson/tsgam") #tsgam not available on CRAN
library(mgcv); library(rio); library(car); library(tsgam); 


#core functions
normcent=function(x,xNuse=1:length(x)) (x-mean(x[xNuse],na.rm=TRUE))/sd(x[xNuse],na.rm=TRUE)
normcentLev=function(x,Lvl,xNuse=rep(TRUE,length(x)),xout=1[0]){ for(i in unique(Lvl)) xout=c(xout,normcent(x[Lvl==i],xNuse[Lvl==i])); return(xout); }
kfindB=function(vr,Dati,Ks=4:12,zmdom=FALSE){ 
	Bs=1[-1]; if(!zmdom) for(k in Ks) Bs=c(Bs,tryCatch(BIC(gam(get(vr)~s(YSI,Wt,k=k),data=Dati)),error=function(e) NA)); 
	if(zmdom) for(k in Ks) Bs=c(Bs,tryCatch(BIC(gam(get(vr)~s(YSI,k=k),data=Dati)),error=function(e) NA))
	Bmin=min(Bs,na.rm=TRUE); Ks[which((Bs-Bmin)<2)[1]]
}
serialAgg=function(x, AggCats, AggTarg = (1:ncol(x))[-AggCats], CatsCode=TRUE, FUN=function(x) mean(x,na.rm=TRUE)){
    ncat = length(AggCats)
    if (ncat == 1) Cats = as.character(x[, AggCats])
    else Cats = codify(x[, AggCats])
    agged = as.matrix(aggregate(x[, AggTarg], by = list(Cats), FUN = FUN))
    if (!CatsCode & ncat > 1) agged = cbind(matrix(unlist(strsplit(agged[, 1], "_")), ncol = ncat, byrow = T), agged[, -1])
    return(agged)
}




#Import data:
dat0=import("byoniDat12-30-22.csv")

#Omit Lake Balaton mussel data in years where the littoral was scraped by ice, and omit the highly turbid Lake Balaton western basin and Lake Batorino.
dat0[dat0$Lake=="BalatonE" & dat0$Year%in%c(2013:2015,2018),c("QMB","ZM_Biom","TD_Biom")]=NA
dat2=dat0[!dat0$Lake%in%c("Batorino","BalatonW"),]

#Standardize time since invasion variable in Lake Balaton to that of other lakes (see Methods)
yb0=14; dat2$YSI[dat2$Lake=="BalatonE"]=yb0:(yb0-1+length(dat2$YSI[dat2$Lake=="BalatonE"]));

#Calculate lake- and variable-specific z-scores
dat2[,-(1:6)]=apply(as.matrix(dat2[,-(1:6)]),2,normcentLev,Lvl=dat2[,"Lake"])

#As we omit Lake Balaton time series before quagga invasion, we need to standardize Lake Balaton z-scores
#to z-scores of other lakes during the quagga mussel invasion period
Lk=as.factor(dat2[,"Lake"]); Cols=rep(c("deepskyblue","red","magenta","blue","darkorange","green","darkolivegreen4","red")[1:7],summary(Lk)[c(2,4,3,5,6,7,1)]);
LB=Lk=="BalatonE"; LQM=Lk%in%c("Oneida","Veluwe","Eem") & dat2$YSI>13; 
for(V in c("TD_Biom","Visibility","Chlorophyll")) dat2[LB,V]=dat2[LB,V]*mean(as.numeric(serialAgg(dat2[LQM,],"Lake",V,FUN=function(y) sd(y,na.rm=TRUE))[,2]),na.rm=TRUE) + mean(dat2[LQM,V],na.rm=TRUE)

#Assign lake colors for plotting and variable names
Cols=rep(c("deepskyblue","red","magenta","blue","darkorange","green","darkolivegreen4","red")[1:7],summary(Lk)[c(2,4,3,5,6,1,7)]);
VARS=c("TD_Biom","Visibility","Chlorophyll","PhytoplanktonBiom","ZoobenthosBiom","ZooplanktonBiom","MactophytesPct","TP")



#Constructs Fig 2. Setting QUAGGA=TRUE, TRENDS=FALSE constructs Fig. 3, and setting QUAGGA=FALSE and TRENDS=TRUE gives Appendix Fig S1
QUAGGA=TRUE; TRENDS=FALSE; 
MaxYrPlot=25+5*(TRENDS); MaxYrUse=30; derifLenSignif=4; par(mfrow=c(2,(4-QUAGGA))); for(i in 1:(8-2*QUAGGA)){
	vr=VARS[i]; Incl=dat2[,"YSI"]<MaxYrUse & !is.na(dat2[,vr]) & (dat2$Lake!="Lukomskoe" | dat2$YSI<21)
	Dat=dat2[Incl,]; Dat$Wt=Dat$ZM_Biom/(Dat$QMB+Dat$ZM_Biom); Dat$Wt[is.na(Dat$Wt)]=(1-(Dat$YSI2>0))[is.na(Dat$Wt)]; Dat$Wt2=Dat$Wt;
	#Correct quagga mussel under-detection in Oneida Lake in 2006-2008 for a slightly more gradual shift in dominance.
	#This slightly reduces quagga mussel impacts but does not affect results significance
	Dat[Dat$Lake=="Oneida" & Dat$YSI%in%(16:19),c("Wt","Wt2")]=cbind(c(0.88,0.7,0.64,0.21),seq(0.88,0.21,len=4)); #

	ki=as.numeric(kfindB(vr,Dat)); xft=gam(get(vr)~s(YSI,Wt,k=ki), data=Dat); #vis.gam(xft,plot.type="contour",color="topo")
	DZ=Dat[Dat$YSI2==0,]; xftderiv=gam(get(vr)~s(YSI,k=ki), data=DZ); QME=BIC(gam(get(vr)~s(YSI,k=ki[1]),data=Dat))-BIC(xft); 
	if(!QUAGGA & !TRENDS){
		#Prop trends explained = 1-sum(SSRs_local)/SSR_global. From this metric exclude lakes with <5 observations.
		ZMdom=Dat$Wt>0.6; SRglob=resid(xft)^2; SRglob[!ZMdom]=0; SSRloc=0; for(lktmp in unique(Dat[,1])){ dtmp=Dat[,1]==lktmp & ZMdom; if(sum(dtmp)<5) SRglob[dtmp]=0 else SSRloc=SSRloc+sum(resid(gam(get(vr)~s(YSI,k=kfindB(vr,Dat[dtmp,],zmdom=TRUE)),data=Dat[dtmp,]))^2); }
		print(vr); SSRglob=sum(SRglob); print(round(c(SSRloc/SSRglob,length(unique(Dat[SRglob!=0,1]))),2));
	}
	
	yrs=0:MaxYrUse; DFN=data.frame(YSI=yrs,YSI2=0,Wt=1); xftdt=fderiv(xftderiv,newdata=DFN); mu=xftdt$derivatives$YSI$deriv; se2=1.96*xftdt$derivatives$YSI$se.deriv;
	dSignif=1+0*yrs; dSignif[(mu+se2)<0]=2; dSignif[(mu-se2)>0]=4; dSignif[yrs>MaxYrPlot]=1;
	dSrle=rle(dSignif); dSignif[rep(dSrle$lengths,dSrle$lengths)<derifLenSignif]=1;
	Prd=predict.gam(xft,newdata=DFN,se.fit=TRUE); Bnd=c(Prd$fit+1.96*Prd$se.fit,rev(Prd$fit-1.96*Prd$se.fit));

	Ttl=paste0(vr,", R2=",round(summary(xft)$r.sq,2),", k=",paste0(ki,collapse="_"),", QMdBIC= ",round(QME,1)) #QMpval=",round(summary(xft)$s.pv[2],3)
	plot(yrs,Prd$fit,col=0,las=1,xlim=c(0,MaxYrPlot),ylim=range(Bnd)*c(1+0.5*(QUAGGA & vr=="Chlorophyll"),1+1*((QUAGGA | TRENDS) & vr%in%c("TD_Biom","ZoobenthosBiom"))+1.25*((QUAGGA | TRENDS) & vr=="Visibility")),xlab="",ylab="",main=Ttl,cex.main=1); 
	polygon(c(yrs,rev(yrs)),Bnd,col="gray60",border=0); box(); if(!QUAGGA & !TRENDS) points(Dat[ZMdom,c("YSI",vr)],col=Cols[Incl][ZMdom],pch=16); 
	rect(col=rgb(1,1,1,alpha=0.7),par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4]);
	if(!TRENDS){
		segs=c(0,cumsum(diff(dSignif)!=0)); for(si in unique(segs)){ sel=segs==si; lines(yrs[sel],Prd$fit[sel],lwd=4-2*QUAGGA,col=c(dSignif[sel][1],1)[1+(QUAGGA)]); }
		if(i==1 & QUAGGA) legend("bottomright",lty=1,col=c("black","indianred"),lwd=3,seg.len=2,c("Mean ZM effect (6 lakes)","Mean QM effect (4 lakes)"),box.col=0)
		if(i==1 & !QUAGGA) legend("bottomright",lty=1,col=c("blue","red","black"),lwd=4,seg.len=2,c("Variable increasing","Variable decreasing","No signficant change"),box.col=0)
	} else {
		for(L in unique(Lk)){ tmp=Dat[Dat[,"Lake"]==L,]; if(sum(!is.na(tmp[,2]))>2){ smoo=supsmu(tmp$YSI,tmp[,vr],span="cv"); qmp=smoo$x>12; X=smoo$x[qmp]; Y=smoo$y[qmp]-100*(max(tmp$YSI2)==0); 
		polygon(c(X,rev(X)),c(0.1+Y,rev(Y-0.1)),border=0,col=rgb(0.8,0.8,0,alpha=0.75)); lines(lwd=2,col=Cols[Incl][Dat[,"Lake"]==L][1],smoo); }; };
		if(i==1 & !QUAGGA) legend("bottomright",seg.len=2.2,box.col=0,lwd=2,col=unique(Cols),c("Lukomskoe","Naroch","Myastro","Oneida","Veluwe","Eem","E. Balaton"))
	}
	if(i==2 & !TRENDS & !QUAGGA) legend("bottomright",seg.len=0.5,box.col=0,lwd=4,col=unique(Cols)[-6],c("Lukomskoe","Naroch","Myastro","Oneida (ZM yrs)","Veluwe (ZM yrs)","Eem (ZM yrs)"));
	if(QUAGGA){
		DQ=Dat[Dat$Lake=="Oneida" & Dat$YSI%in%(13:26),]; DQ$Wt=DQ$Wt2; Prd2=predict.gam(xft,newdata=DQ,se.fit=TRUE); lines(DQ$YSI,Prd2$fit,lwd=3,col="indianred");
		Bnd2=c(Prd2$fit+1.96*Prd2$se.fit,rev(Prd2$fit-1.96*Prd2$se.fit)); polygon(c(DQ$YSI,rev(DQ$YSI)),Bnd2,col=rgb(.8,.36,.36,alpha=0.25),border=0); 
	}
}


#Figure S2: trends and data separately by variable and lake
varsNm=c("Mussels","Secchi Depth","Chlorophyll","Phytoplankton","Zoobenthos","Zooplankton","Macrophytes %","Total Phosphorus")
par(mfcol=c(7,8)); for(i in 1:8){
	vr=VARS[i]; Incl=!is.na(dat2[,vr]) & (dat2$Lake!="Lukomskoe" | dat2$YSI<21)
	Dat=dat2[Incl,]; YLM=quantile(Dat[,vr],c(0.02,0.98),na.rm=TRUE);
	for(L in unique(Lk)){
		MAR=c(3,6,1.5,0.1); VN=varsNm[i]; if(L!="Lukomskoe"){ MAR[3]=0.1; VN=""; }
		LN=as.character(L); if(i!=1){ MAR[2]=3; LN=""; }; par(mar=MAR);		
		tmp=Dat[Dat[,"Lake"]==L,];  if(length(tmp[,vr])<2){ plot(0,col=0,axes=FALSE,ylab=LN,xlab="",main=VN,cex.main=2,cex.lab=2); next; } 
		smoo=supsmu(tmp$YSI,tmp[,vr],span="cv"); qmp=smoo$x>12; X=smoo$x[qmp]; Y=smoo$y[qmp]#-100*(max(tmp$YSI2)==0);
		plot(0:MaxYrUse,bty="n",xlim=c(0,32),ylim=YLM,col=0,main=VN,ylab=LN,xlab="",cex.main=2,cex.lab=2);
		if(max(tmp$YSI2)>0) lines(X,Y,col="lightgray",lwd=15)
		# polygon(c(X,rev(X)),c(0.2+Y,rev(Y-0.2)),border=0,col=rgb(0.3,0.3,0.3,alpha=0.5));
		lines(lwd=2,col=Cols[Incl][Dat[,"Lake"]==L][1],smoo);
		points(tmp[,"YSI"],tmp[,vr],pch=21,col=1,bg=Cols[Incl][Dat[,"Lake"]==L][1])
	}
}





#Figure 4 results - ZM VERSUS QM EFFECTS APPROXIMATION
VarsHires0=import("Oneida Seasonal 2021 short.csv")[,-1]; 
VarsHires=VarsHires0[which(VarsHires0[,1]%in%(1975:2017) & VarsHires0[,2]%in%(3:11)),]
VarsAnns1=cbind(data.frame(apply(serialAgg(VarsHires,1),2,as.numeric)),ZMB=0,QMB=0)
qd=dat0[dat0$Lake=="Oneida",]; VarsAnns1[VarsAnns1[,1]>1991,c("ZMB","QMB")]=qd[-(1:2),c("ZM_Biom","QMB")]; VarsAnns1[VarsAnns1[,1]==1991,"ZMB"]=NA;
VarsAnns4=VarsAnns1[-(1:11),]; VarsAnns4[,3:6]=apply(VarsAnns4[,3:6],2,normcent)


VARS2=c("AvgOfSecchi","AvgOfChla","AvgOfTP","Phytoplankton")
STATS=matrix(0,nrow=4,ncol=8); rownames(STATS)=VARS2; 
colnames(STATS)=c("qmEffect","qmSE","qmPvalue","zmEffect","zmSE","zmPvalue","totalR-squared","EffectDifferencePvalue");
for(i in 1:4){ LM=lm(get(VARS2[i])~QMB+ZMB,data=VarsAnns4); lms=summary(LM); 
	STATS[i,]=c(as.vector(t(lms$coefficients[-1,-3])), lms$r.sq, linearHypothesis(LM,"ZMB = QMB")[[6]][2]); }
write.csv(STATS,"DreissenaSppEffectsTest3zscores1986.csv")






#Invasion statistics
setwd("C:\\\\Users\\\\vkara\\\\Dropbox\\\\Projects\\\\quaggaGL\\\\BYONI IMPACTS\\\\Data files")
zm=import("ZMinvasions.csv"); zmempt=qm$Locality=="[No locality description provided]"; zm$Locality[zmempt]=as.character(zm[zmempt,16])
qm=import("QMinvasions.csv"); MIsta=nchar(qm$LocalityFix)>0; qm$Locality[MIsta]=qm$LocalityFix[MIsta]; qmempt=qm$Locality=="[No locality description provided]"; qm$Locality[qmempt]=as.character(qm[qmempt,16])
zmu=zm[!duplicated(substr(zm$Locality,1,10)),]; qmu=qm[!duplicated(substr(qm$Locality,1,10)),]
two=qmu$Locality[qmu$Locality%in%zmu$Locality]; c(nrow(zmu),nrow(qmu),length(two))
twozmy=apply(cbind(zmu$Year[zmu$Locality%in%two],qmu$Year[qmu$Locality%in%two]),1,min)
tbco=table(twozmy); tbco; #sum(tbco[as.numeric(names(tbco))<2010]); sum(tbco[as.numeric(names(tbco))>2010]); #plot(as.numeric(names(tbco)),cumsum(tbco))
# library(stringdist); ClosestMatch2=function(string,stringVector,MD=3) cbind(string,stringVector[amatch(string, stringVector, maxDist=MD)])
# Matched=ClosestMatch2(qmu$Locality,zmu$Locality); sum(!is.na(Matched[,2]))


qpd=qpdn=import("QMpd.csv")
qpdn[,-(1:3)]=apply(as.matrix(qpdn[,-(1:3)]),2,normcentLev,Lvl=qpdn[,"Lake"]);
par(mfrow=c(2,3)); for(V in 4:8) gampred(qpd[,c(3,V)],1,2,k=3:9,pch=16,col=as.numeric(as.factor(qpd[,1])),seCol="gray",xlab="Year since invasion",ylab=colnames(qpd)[V]); plot(1,col=0); legend("center",col=1:3,lwd=8,seg.len=0.5,cex=2,c("Huron","Michigan","Ontario"))
par(mfrow=c(2,3)); for(V in 4:8) gampred(qpdn[,c(3,V)],1,2,k=3:9,pch=16,col=as.numeric(as.factor(qpdn[,1])),seCol="gray",xlab="Year since invasion",ylab=colnames(qpd)[V],ylim=c(-1.2,2)); plot(1); legend("center",col=1:3,lwd=8,seg.len=0.5,cex=2,c("Huron","Michigan","Ontario"))





