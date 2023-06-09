#extract environmental data for Saenger and Evans, 2019 global coretop calibration
  

setwd("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data");

library(ncdf4)
library(abind)

#constants
OmegaCcut<-1.8																				#threshold below which calcite saturation has an effect

#coretop foram MgCa data
MgCa.dat.Npachy<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Npachys.csv")
MgCa.dat.Gruber<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Gruber.csv")
MgCa.dat.Ginflata<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Ginflata.csv")
MgCa.dat.Gbulloides<-read.csv("/Volumes/GoogleDrive/My Drive/R_documents/PSR/proxy_data/ForamMg_Gbulloides.csv")

#envrionmental data
#salinity
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_s00_01v2.nc")
s.ann<-ncvar_get(nc,"s_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))	
s.ann.sd<-ncvar_get(nc,"s_sd",start=c(1,1,1,1),count=c(-1,-1,40,-1))
depth<-nc$dim$depth$vals[1:40]
lat<-nc$dim$lat$vals
lon<-nc$dim$lon$vals
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_s13_01v2.nc")
s.winter.tmp<-ncvar_get(nc,"s_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.winter.tmp.sd<-ncvar_get(nc,"s_sd",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.winterNH<-s.winter.tmp[,which(lat>=0),]
s.summerSH<-s.winter.tmp[,which(lat<0),]
s.winterNH.sd<-s.winter.tmp.sd[,which(lat>=0),]
s.summerSH.sd<-s.winter.tmp.sd[,which(lat<0),]
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_s15_01v2.nc")
s.summer.tmp<-ncvar_get(nc,"s_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.summer.tmp.sd<-ncvar_get(nc,"s_sd",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.summerNH<-s.summer.tmp[,which(lat>=0),]
s.winterSH<-s.summer.tmp[,which(lat<0),]
s.summerNH.sd<-s.summer.tmp.sd[,which(lat>=0),]
s.winterSH.sd<-s.summer.tmp.sd[,which(lat<0),]
rm(nc)
s.winter<-abind(s.winterSH,s.winterNH,along=2)
s.winter.sd<-abind(s.winterSH.sd,s.winterNH.sd,along=2)
s.summer<-abind(s.summerSH,s.summerNH,along=2)
s.summer.sd<-abind(s.summerSH.sd,s.summerNH.sd,along=2)

nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_s14_01v2.nc")
s.spring.tmp<-ncvar_get(nc,"s_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.spring.tmp.sd<-ncvar_get(nc,"s_sd",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.springNH<-s.spring.tmp[,which(lat>=0),]
s.springNH.sd<-s.spring.tmp.sd[,which(lat>=0),]
s.fallSH<-s.spring.tmp[,which(lat<0),]
s.fallSH.sd<-s.spring.tmp.sd[,which(lat<0),]
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_s16_01v2.nc")
s.fall.tmp<-ncvar_get(nc,"s_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.fall.tmp.sd<-ncvar_get(nc,"s_sd",start=c(1,1,1,1),count=c(-1,-1,40,-1))
s.fallNH<-s.fall.tmp[,which(lat>=0),]
s.fallNH.sd<-s.fall.tmp.sd[,which(lat>=0),]
s.springSH<-s.fall.tmp[,which(lat<0),]
s.springSH.sd<-s.fall.tmp.sd[,which(lat<0),]
rm(nc)
s.spring<-abind(s.springSH,s.springNH,along=2)
s.spring.sd<-abind(s.springSH.sd,s.springNH.sd,along=2)
s.fall<-abind(s.fallSH,s.fallNH,along=2)
s.fall.sd<-abind(s.fallSH.sd,s.fallNH.sd,along=2)

#temperature
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_t00_01v2.nc")
t.ann<-ncvar_get(nc,"t_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))	
t.ann.sd<-ncvar_get(nc,"t_sd",start=c(1,1,1,1),count=c(-1,-1,40,-1))
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_t13_01v2.nc")
t.winter.tmp<-ncvar_get(nc,"t_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.winter.tmp.sd<-ncvar_get(nc,"t_sd",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.winterNH<-t.winter.tmp[,which(lat>=0),]
t.summerSH<-t.winter.tmp[,which(lat<0),]
t.winterNH.sd<-t.winter.tmp.sd[,which(lat>=0),]
t.summerSH.sd<-t.winter.tmp.sd[,which(lat<0),]
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_t15_01v2.nc")
t.summer.tmp<-ncvar_get(nc,"t_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.summer.tmp.sd<-ncvar_get(nc,"t_sd",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.summerNH<-t.summer.tmp[,which(lat>=0),]
t.winterSH<-t.summer.tmp[,which(lat<0),]
t.summerNH.sd<-t.summer.tmp.sd[,which(lat>=0),]
t.winterSH.sd<-t.summer.tmp.sd[,which(lat<0),]
rm(nc)
t.winter<-abind(t.winterSH,t.winterNH,along=2)
t.winter.sd<-abind(t.winterSH.sd,t.winterNH.sd,along=2)
t.summer<-abind(t.summerSH,t.summerNH,along=2)
t.summer.sd<-abind(t.summerSH.sd,t.summerNH.sd,along=2)

nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_t14_01v2.nc")
t.spring.tmp<-ncvar_get(nc,"t_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.spring.tmp.sd<-ncvar_get(nc,"t_sd",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.springNH<-t.spring.tmp[,which(lat>=0),]
t.springNH.sd<-t.spring.tmp.sd[,which(lat>=0),]
t.fallSH<-t.spring.tmp[,which(lat<0),]
t.fallSH.sd<-t.spring.tmp.sd[,which(lat<0),]
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/WOA/woa13_decav_t16_01v2.nc")
t.fall.tmp<-ncvar_get(nc,"t_an",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.fall.tmp.sd<-ncvar_get(nc,"t_sd",start=c(1,1,1,1),count=c(-1,-1,40,-1))
t.fallNH<-t.fall.tmp[,which(lat>=0),]
t.fallNH.sd<-t.fall.tmp.sd[,which(lat>=0),]
t.springSH<-t.fall.tmp[,which(lat<0),]
t.springSH.sd<-t.fall.tmp.sd[,which(lat<0),]
rm(nc)
t.spring<-abind(t.springSH,t.springNH,along=2)
t.spring.sd<-abind(t.springSH.sd,t.springNH.sd,along=2)
t.fall<-abind(t.fallSH,t.fallNH,along=2)
t.fall.sd<-abind(t.fallSH.sd,t.fallNH.sd,along=2)

#saturation state (NOTE: GLODAP has same latitude, but different longitude)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/GLODAP/GLODAPv2.2016b.OmegaC.nc")
OmegaC<-ncvar_get(nc,"OmegaC",start=c(1,1,1),count=c(-1,-1,-1))	
OmegaC.sd<-ncvar_get(nc,"OmegaC_error",start=c(1,1,1),count=c(-1,-1,-1))
lon.glodap<-nc$dim$lon$vals
lon.glodap[which(lon.glodap>180)]<-lon.glodap[which(lon.glodap>180)]-360
depth.glodap<-ncvar_get(nc,"Depth",start=c(1),count=c(-1))	
rm(nc)
nc<-nc_open("/Volumes/GoogleDrive/My Drive/R_documents/PSR/GLODAP/GLODAPv2.2016b.pHtsinsitutp.nc")
pH<-ncvar_get(nc,"pHtsinsitutp",start=c(1,1,1),count=c(-1,-1,-1))	
pH.sd<-ncvar_get(nc,"pHtsinsitutp_error",start=c(1,1,1),count=c(-1,-1,-1))
rm(nc)

#***********************************************#

#get WOA data at appropriate depth and seasonality for N. pachy & N. incompta. Do not average depth
Npachy.d<-c(0,100)																									#assigned depth habitat
d1<-which(depth.glodap<=max(Npachy.d) & depth.glodap>=min(Npachy.d))
d<-match(depth.glodap[d1],depth)
Npachy.ta<-t.ann[,,d]
Npachy.ta.sd<-t.ann.sd[,,d]
Npachy.ta.mean<-apply(Npachy.ta,c(1,2),mean)
Npachy.t1<-t.summer[,,d]
Npachy.t1.sd<-t.summer.sd[,,d]
Npachy.t2<-t.spring[,,d]
Npachy.t2.sd<-t.spring.sd[,,d]
Npachy.s1<-s.summer[,,d]
Npachy.s1.sd<-s.summer.sd[,,d]
Npachy.s2<-s.spring[,,d]
Npachy.s2.sd<-s.spring.sd[,,d]

#Start with spring, but if mean annual T<=10, replace with summer based on Jonkers 15
Npachy.t<-Npachy.t2
Npachy.t.sd<-Npachy.t2.sd
tmp<-which(Npachy.ta.mean<=10)																			#indicies of gridboxes with mean T<10ºC
Npachy.t[tmp]<-Npachy.t1[tmp]
Npachy.t.sd[tmp]<-Npachy.t1.sd[tmp]
Npachy.s<-Npachy.s2
Npachy.s.sd<-Npachy.s2.sd
Npachy.s[tmp]<-Npachy.s1[tmp]
Npachy.s.sd[tmp]<-Npachy.s1.sd[tmp]

saveRDS(Npachy.s,"Npachy.s")
saveRDS(Npachy.s.sd,"Npachy.s.sd")

Npachy.omega.s<-OmegaC[,,d1]
Npachy.omega.s.sd<-OmegaC.sd[,,d1]
Npachy.pH.s<-pH[,,d1]
Npachy.pH.s.sd<-pH.sd[,,d1]

Npachy.WOA.s<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.WOA.s.sd<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.WOA.t<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.WOA.t.sd<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.glodap.omega.s<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.glodap.omega.s.sd<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.glodap.omega<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.glodap.omega.sd<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.glodap.pH.s<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.glodap.pH.s.sd<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.glodap.pH<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))
Npachy.glodap.pH.sd<-matrix(,nrow(MgCa.dat.Npachy),ncol=length(d))

for(i in 1:nrow(MgCa.dat.Npachy)) {
	#WOA T and S
	lon.tmp<-which(abs(lon-MgCa.dat.Npachy$lon[i])==min(abs(lon-MgCa.dat.Npachy$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Npachy$lat[i])==min(abs(lat-MgCa.dat.Npachy$lat[i])))[1]
	if (is.na(Npachy.s[lon.tmp,lat.tmp,1])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																		#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Npachy.WOA.s[i,]<-apply(Npachy.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)									
			Npachy.WOA.t[i,]<-apply(Npachy.t[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Npachy.WOA.s.sd[i,]<-apply(Npachy.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Npachy.WOA.t.sd[i,]<-apply(Npachy.t.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
		} else {
			Npachy.WOA.s[i,]<-Npachy.s[lon.tmp,lat.tmp,]
			Npachy.WOA.t[i,]<-Npachy.t[lon.tmp,lat.tmp,]
			Npachy.WOA.s.sd[i,]<-Npachy.s.sd[lon.tmp,lat.tmp,]
			Npachy.WOA.t.sd[i,]<-Npachy.t.sd[lon.tmp,lat.tmp,]
		}
	
	#GLODAP Omega and pH
	lon.tmp<-which(abs(lon.glodap-MgCa.dat.Npachy$lon[i])==min(abs(lon.glodap-MgCa.dat.Npachy$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Npachy$lat[i])==min(abs(lat-MgCa.dat.Npachy$lat[i])))[1]
	tmp<-depth.glodap-MgCa.dat.Npachy$depth[i]
	depth.tmp<-which(tmp==max(tmp[tmp<0]))[1]																					#find the glodap depth closest (without going over) to the core depth
	if (is.na(Npachy.omega.s[lon.tmp,lat.tmp,1])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																				#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Npachy.glodap.omega.s[i,]<-apply(Npachy.omega.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Npachy.glodap.omega.s.sd[i,]<-apply(Npachy.omega.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Npachy.glodap.pH.s[i,]<-apply(Npachy.pH.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Npachy.glodap.pH.s.sd[i,]<-apply(Npachy.pH.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Npachy.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			#if there's no carbonate data at this depth move up a bit
			if (all(is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp]))==TRUE) {
			depth.tmp<-ifelse(depth.tmp<25,depth.tmp-2,depth.tmp-1)
			Npachy.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))		
			} 
		} else {
			Npachy.glodap.omega.s[i,]<-Npachy.omega.s[lon.tmp,lat.tmp,]
			Npachy.glodap.omega.s.sd[i,]<-Npachy.omega.s.sd[lon.tmp,lat.tmp,]
			Npachy.glodap.pH.s[i,]<-Npachy.pH.s[lon.tmp,lat.tmp,]
			Npachy.glodap.pH.s.sd[i,]<-Npachy.pH.s.sd[lon.tmp,lat.tmp,]
			Npachy.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			if (is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp])==TRUE) {
			depth.tmp<-ifelse(depth.tmp<25,depth.tmp-2,depth.tmp-1)				
			Npachy.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Npachy.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			}
		}
	}


MgCa.dat.Npachy$WOA.s<-Npachy.WOA.s
MgCa.dat.Npachy$WOA.t<-Npachy.WOA.t
MgCa.dat.Npachy$GLODAP.OmegaC.surf<-Npachy.glodap.omega.s
MgCa.dat.Npachy$GLODAP.pH.surf<-Npachy.glodap.pH.s
MgCa.dat.Npachy$GLODAP.OmegaC.deep<-Npachy.glodap.omega
MgCa.dat.Npachy$GLODAP.pH.deep<-Npachy.glodap.pH

MgCa.dat.Npachy$WOA.s.sd<-Npachy.WOA.s.sd
MgCa.dat.Npachy$WOA.s.sd[which(is.na(rowMeans(Npachy.WOA.s.sd))==TRUE | rowMeans(Npachy.WOA.s.sd)==0),]<-colMeans(Npachy.WOA.s.sd,na.rm=T)
MgCa.dat.Npachy$WOA.t.sd<-Npachy.WOA.t.sd
MgCa.dat.Npachy$WOA.t.sd[which(is.na(rowMeans(Npachy.WOA.t.sd))==TRUE | rowMeans(Npachy.WOA.t.sd)==0),]<-colMeans(Npachy.WOA.t.sd,na.rm=T)
MgCa.dat.Npachy$GLODAP.OmegaC.surf.sd<-Npachy.glodap.omega.s.sd
MgCa.dat.Npachy$GLODAP.OmegaC.surf.sd[which(is.na(rowMeans(Npachy.glodap.omega.s.sd))==TRUE | rowMeans(Npachy.glodap.omega.s.sd)==0),]<-colMeans(Npachy.glodap.omega.s.sd,na.rm=T)
MgCa.dat.Npachy$GLODAP.pH.surf.sd<-Npachy.glodap.pH.s.sd
MgCa.dat.Npachy$GLODAP.pH.surf.sd[which(is.na(rowMeans(Npachy.glodap.pH.s.sd))==TRUE | rowMeans(Npachy.glodap.pH.s.sd)==0),]<-colMeans(Npachy.glodap.pH.s.sd,na.rm=T)
MgCa.dat.Npachy$GLODAP.OmegaC.deep.sd<-Npachy.glodap.omega.sd
MgCa.dat.Npachy$GLODAP.OmegaC.deep.sd[which(is.na(rowMeans(Npachy.glodap.omega.sd))==TRUE | rowMeans(Npachy.glodap.omega.sd)==0),]<-colMeans(Npachy.glodap.omega.sd,na.rm=T)
MgCa.dat.Npachy$GLODAP.pH.deep.sd<-Npachy.glodap.pH.sd
MgCa.dat.Npachy$GLODAP.pH.deep.sd[which(is.na(rowMeans(Npachy.glodap.pH.sd))==TRUE | rowMeans(Npachy.glodap.pH.sd)==0),]<-colMeans(Npachy.glodap.pH.sd,na.rm=T)

#set deep values to the critical value when they exceed the critical value. 
cut<-which(Npachy.glodap.omega[,1]<=OmegaCcut)
MgCa.dat.Npachy$GLODAP.OmegaC.deep.cut<-rep(OmegaCcut,length.out=nrow(MgCa.dat.Npachy))
MgCa.dat.Npachy$GLODAP.OmegaC.deep.cut[cut]<-MgCa.dat.Npachy$GLODAP.OmegaC.deep[cut]
MgCa.dat.Npachy$GLODAP.OmegaC.deep.cut[which(is.na(Npachy.glodap.omega[,1])==TRUE & MgCa.dat.Npachy$depth>2000)]<-NA
MgCa.dat.Npachy$GLODAP.pH.deep.cut<-rep(mean(MgCa.dat.Npachy$GLODAP.pH.deep[-cut],na.rm=T),length.out=nrow(MgCa.dat.Npachy))
MgCa.dat.Npachy$GLODAP.pH.deep.cut[cut]<-MgCa.dat.Npachy$GLODAP.pH.deep[cut]
MgCa.dat.Npachy$GLODAP.pH.deep.cut[which(is.na(Npachy.glodap.pH[,1])==TRUE & MgCa.dat.Npachy$depth>2000)]<-NA

#clean up to remove rows without environmental data
#This excludes 9 records
no.data<-which(is.na(rowMeans(Npachy.WOA.s))==TRUE | is.na(rowMeans(Npachy.WOA.t))==TRUE | is.na(rowMeans(Npachy.glodap.omega.s))==TRUE | is.na(rowMeans(Npachy.glodap.pH.s))==TRUE | is.na(MgCa.dat.Npachy$GLODAP.OmegaC.deep.cut)==TRUE)
MgCa.dat.Npachy<-MgCa.dat.Npachy[-no.data,]

#compute a correlation matrix
Npachy.subset<-data.frame(s=as.vector(MgCa.dat.Npachy[,12]),																	#salinity
			t=as.vector(MgCa.dat.Npachy[,13]),																					#temperature
			OmegaC.surf=as.vector(MgCa.dat.Npachy[,14]),																			#omega surface
			pH.surf=as.vector(MgCa.dat.Npachy[,15]),																				#pH surface
			OmegaC.deep=as.vector(matrix(MgCa.dat.Npachy$GLODAP.OmegaC.deep.cut,nrow=nrow(MgCa.dat.Npachy),ncol=length(d))),	#omega deep after cutoff
			pH.deep=as.vector(matrix(MgCa.dat.Npachy$GLODAP.pH.deep.cut,nrow=nrow(MgCa.dat.Npachy),ncol=length(d)))				#pH deep after cutoff	
			)

Npachy.cor.mat<-round(cor(Npachy.subset),2)

#subset Npachy
MgCa.dat.NpachyR<-MgCa.dat.Npachy[which(MgCa.dat.Npachy$Notes=="Dextral"),]
MgCa.dat.NpachyL<-MgCa.dat.Npachy[which(MgCa.dat.Npachy$Notes=="Sinstral"),]

#compute a correlation matrix for NpachyR
NpachyR.subset<-data.frame(s=as.vector(MgCa.dat.NpachyR[,12]),																	#salinity
			t=as.vector(MgCa.dat.NpachyR[,13]),																					#temperature
			OmegaC.surf=as.vector(MgCa.dat.NpachyR[,14]),																			#omega surface
			pH.surf=as.vector(MgCa.dat.NpachyR[,15]),																				#pH surface
			OmegaC.deep=as.vector(matrix(MgCa.dat.NpachyR$GLODAP.OmegaC.deep.cut,nrow=nrow(MgCa.dat.NpachyR),ncol=length(d))),	#omega deep after cutoff
			pH.deep=as.vector(matrix(MgCa.dat.NpachyR$GLODAP.pH.deep.cut,nrow=nrow(MgCa.dat.NpachyR),ncol=length(d)))				#pH deep after cutoff	
			)

NpachyR.cor.mat<-round(cor(NpachyR.subset),2)

#compute a correlation matrix for NpachyL
NpachyL.subset<-data.frame(s=as.vector(MgCa.dat.NpachyL[,12]),																	#salinity
			t=as.vector(MgCa.dat.NpachyL[,13]),																					#temperature
			OmegaC.surf=as.vector(MgCa.dat.NpachyL[,14]),																			#omega surface
			pH.surf=as.vector(MgCa.dat.NpachyL[,15]),																				#pH surface
			OmegaC.deep=as.vector(matrix(MgCa.dat.NpachyL$GLODAP.OmegaC.deep.cut,nrow=nrow(MgCa.dat.NpachyL),ncol=length(d))),	#omega deep after cutoff
			pH.deep=as.vector(matrix(MgCa.dat.NpachyL$GLODAP.pH.deep.cut,nrow=nrow(MgCa.dat.NpachyL),ncol=length(d)))				#pH deep after cutoff	
			)

NpachyL.cor.mat<-round(cor(NpachyL.subset),2)

#save environmental data and correlation matrix

write.csv(MgCa.dat.Npachy,file="Npachy.MgCa.env.dat.csv")
write.csv(Npachy.cor.mat,file="Npachy.cor.mat.csv")
Npachy95sig.cor<-2*nrow(MgCa.dat.Npachy)^-0.5

write.csv(NpachyR.cor.mat,file="NpachyR.cor.mat.csv")
write.csv(NpachyL.cor.mat,file="NpachyL.cor.mat.csv")

dev.new(width=6,height=8)
par(mfrow=c(3,2))
plot(NpachyR.subset$s,NpachyR.subset$t,xlab="salinity",ylab="temperature",main="N. pachyderma; N. incompta: 0-100m",las=1,xlim=c(33.5,36.5),ylim=c(-5,25))
points(NpachyL.subset$s,NpachyL.subset$t,pch=5,col="grey50")
plot(NpachyR.subset$s,NpachyR.subset$OmegaC.surf,xlab="salinity",ylab="surface saturation state",las=1,xlim=c(33.5,36.5),ylim=c(2,5.5))
points(NpachyL.subset$s,NpachyL.subset$OmegaC.surf,pch=5,col="grey50")
plot(NpachyR.subset$s,NpachyR.subset$pH.surf,xlab="salinity",ylab="pH",las=1,xlim=c(33.5,36.5),ylim=c(8,8.3))
points(NpachyL.subset$s,NpachyL.subset$pH.surf,pch=5,col="grey50")
plot(NpachyR.subset$t,NpachyR.subset$OmegaC.surf,xlab="temperature",ylab="surface saturation state",las=1,xlim=c(-5,25),ylim=c(2,5.5))
points(NpachyL.subset$t,NpachyL.subset$OmegaC.surf,pch=5,col="grey50")
plot(NpachyR.subset$t,NpachyR.subset$pH.surf,xlab="temperature",ylab="pH",las=1,xlim=c(-5,25),ylim=c(8,8.3))
points(NpachyL.subset$t,NpachyL.subset$pH.surf,pch=5,col="grey50")
plot(NpachyR.subset$pH.surf,NpachyR.subset$OmegaC.surf,xlab="pH",ylab="surface saturation state",las=1,xlim=c(8,8.3),ylim=c(2,5.5))
points(NpachyL.subset$pH.surf,NpachyL.subset$OmegaC.surf,pch=5,col="grey50")

#***********************************************#
#coretop calibration of G.ruber
Gruber.d<-c(0,50)
d1<-which(depth.glodap<=max(Gruber.d) & depth.glodap>=min(Gruber.d))
d<-match(depth.glodap[d1],depth)
Gruber.ta<-t.ann[,,d]
Gruber.ta.sd<-t.ann.sd[,,d]
Gruber.ta.mean<-apply(Gruber.ta,c(1,2),mean)
Gruber.t1<-Gruber.ta
Gruber.t1.sd<-Gruber.ta.sd
Gruber.t2<-t.fall[,,d]
Gruber.t2.sd<-t.fall.sd[,,d]
Gruber.s1<-s.ann[,,d]
Gruber.s1.sd<-s.ann.sd[,,d]
Gruber.s2<-s.fall[,,d]
Gruber.s2.sd<-s.fall.sd[,,d]

#start with mean annual, but if mean annual <=25, replace with fall based on Jonkers 15
Gruber.t<-Gruber.t1
Gruber.t.sd<-Gruber.t1.sd
tmp<-which(Gruber.ta.mean<=25)																			#indicies of gridboxes with mean T<25ºC
Gruber.t[tmp]<-Gruber.t2[tmp]
Gruber.t.sd[tmp]<-Gruber.t2.sd[tmp]
Gruber.s<-Gruber.s1
Gruber.s.sd<-Gruber.s1.sd
Gruber.s[tmp]<-Gruber.s2[tmp]
Gruber.s.sd[tmp]<-Gruber.s2.sd[tmp]

Gruber.omega.s<-OmegaC[,,d1]
Gruber.omega.s.sd<-OmegaC.sd[,,d1]
Gruber.pH.s<-pH[,,d1]
Gruber.pH.s.sd<-pH.sd[,,d1]

Gruber.WOA.s<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.WOA.s.sd<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.WOA.t<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.WOA.t.sd<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.glodap.omega.s<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.glodap.omega.s.sd<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.glodap.omega<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.glodap.omega.sd<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.glodap.pH.s<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.glodap.pH.s.sd<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.glodap.pH<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))
Gruber.glodap.pH.sd<-matrix(,nrow(MgCa.dat.Gruber),ncol=length(d))

for(i in 1:nrow(MgCa.dat.Gruber)) {
	#WOA T and S
	lon.tmp<-which(abs(lon-MgCa.dat.Gruber$lon[i])==min(abs(lon-MgCa.dat.Gruber$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Gruber$lat[i])==min(abs(lat-MgCa.dat.Gruber$lat[i])))[1]
	if (is.na(Gruber.s[lon.tmp,lat.tmp,1])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																		#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Gruber.WOA.s[i,]<-apply(Gruber.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gruber.WOA.t[i,]<-apply(Gruber.t[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gruber.WOA.s.sd[i,]<-apply(Gruber.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gruber.WOA.t.sd[i,]<-apply(Gruber.t.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
		} else {
			Gruber.WOA.s[i,]<-Gruber.s[lon.tmp,lat.tmp,]
			Gruber.WOA.t[i,]<-Gruber.t[lon.tmp,lat.tmp,]
			Gruber.WOA.s.sd[i,]<-Gruber.s.sd[lon.tmp,lat.tmp,]
			Gruber.WOA.t.sd[i,]<-Gruber.t.sd[lon.tmp,lat.tmp,]
		}
	#GLODAP Omega and pH
	lon.tmp<-which(abs(lon.glodap-MgCa.dat.Gruber$lon[i])==min(abs(lon.glodap-MgCa.dat.Gruber$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Gruber$lat[i])==min(abs(lat-MgCa.dat.Gruber$lat[i])))[1]
	tmp<-depth.glodap-MgCa.dat.Gruber$depth[i]
	depth.tmp<-which(tmp==max(tmp[tmp<0]))[1]																					#find the glodap depth closest (without going over) to the core depth
	if (is.na(Gruber.omega.s[lon.tmp,lat.tmp,1])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																				#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Gruber.glodap.omega.s[i,]<-apply(Gruber.omega.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gruber.glodap.omega.s.sd[i,]<-apply(Gruber.omega.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gruber.glodap.pH.s[i,]<-apply(Gruber.pH.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gruber.glodap.pH.s.sd[i,]<-apply(Gruber.pH.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gruber.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			#if there's no carbonate data at this depth move up a bit
			if (all(is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp]))==TRUE) {
			depth.tmp<-ifelse(depth.tmp<25,depth.tmp-2,depth.tmp-1)
			Gruber.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))	
			Gruber.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))	
			} 
		} else {
			Gruber.glodap.omega.s[i,]<-Gruber.omega.s[lon.tmp,lat.tmp,]
			Gruber.glodap.omega.s.sd[i,]<-Gruber.omega.s.sd[lon.tmp,lat.tmp,]
			Gruber.glodap.pH.s[i,]<-Gruber.pH.s[lon.tmp,lat.tmp,]
			Gruber.glodap.pH.s.sd[i,]<-Gruber.pH.s.sd[lon.tmp,lat.tmp,]
			Gruber.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			if (is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp])==TRUE) {
			depth.tmp<-ifelse(depth.tmp<25,depth.tmp-2,depth.tmp-1)				
			Gruber.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gruber.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			}
		}
	}

MgCa.dat.Gruber$WOA.s<-Gruber.WOA.s
MgCa.dat.Gruber$WOA.t<-Gruber.WOA.t
MgCa.dat.Gruber$GLODAP.OmegaC.surf<-Gruber.glodap.omega.s
MgCa.dat.Gruber$GLODAP.pH.surf<-Gruber.glodap.pH.s
MgCa.dat.Gruber$GLODAP.OmegaC.deep<-Gruber.glodap.omega
MgCa.dat.Gruber$GLODAP.pH.deep<-Gruber.glodap.pH

MgCa.dat.Gruber$WOA.s.sd<-Gruber.WOA.s.sd
MgCa.dat.Gruber$WOA.s.sd[which(is.na(rowMeans(Gruber.WOA.s.sd))==TRUE | rowMeans(Gruber.WOA.s.sd)==0),]<-colMeans(Gruber.WOA.s.sd,na.rm=T)
MgCa.dat.Gruber$WOA.t.sd<-Gruber.WOA.t.sd
MgCa.dat.Gruber$WOA.t.sd[which(is.na(rowMeans(Gruber.WOA.t.sd))==TRUE | rowMeans(Gruber.WOA.t.sd)==0),]<-colMeans(Gruber.WOA.t.sd,na.rm=T)
MgCa.dat.Gruber$GLODAP.OmegaC.surf.sd<-Gruber.glodap.omega.s.sd
MgCa.dat.Gruber$GLODAP.OmegaC.surf.sd[which(is.na(rowMeans(Gruber.glodap.omega.s.sd))==TRUE | rowMeans(Gruber.glodap.omega.s.sd)==0),]<-colMeans(Gruber.glodap.omega.s.sd,na.rm=T)
MgCa.dat.Gruber$GLODAP.pH.surf.sd<-Gruber.glodap.pH.s.sd
MgCa.dat.Gruber$GLODAP.pH.surf.sd[which(is.na(rowMeans(Gruber.glodap.pH.s.sd))==TRUE | rowMeans(Gruber.glodap.pH.s.sd)==0),]<-colMeans(Gruber.glodap.pH.s.sd,na.rm=T)
MgCa.dat.Gruber$GLODAP.OmegaC.deep.sd<-Gruber.glodap.omega.sd
MgCa.dat.Gruber$GLODAP.OmegaC.deep.sd[which(is.na(rowMeans(Gruber.glodap.omega.sd))==TRUE | rowMeans(Gruber.glodap.omega.sd)==0),]<-colMeans(Gruber.glodap.omega.sd,na.rm=T)
MgCa.dat.Gruber$GLODAP.pH.deep.sd<-Gruber.glodap.pH.sd
MgCa.dat.Gruber$GLODAP.pH.deep.sd[which(is.na(rowMeans(Gruber.glodap.pH.sd))==TRUE | rowMeans(Gruber.glodap.pH.sd)==0),]<-colMeans(Gruber.glodap.pH.sd,na.rm=T)

#set deep values to the critical value when they exceed the critical value. 
cut<-which(Gruber.glodap.omega[,1]<=OmegaCcut)
MgCa.dat.Gruber$GLODAP.OmegaC.deep.cut<-rep(OmegaCcut,length.out=nrow(MgCa.dat.Gruber))
MgCa.dat.Gruber$GLODAP.OmegaC.deep.cut[cut]<-MgCa.dat.Gruber$GLODAP.OmegaC.deep[cut]
MgCa.dat.Gruber$GLODAP.OmegaC.deep.cut[which(is.na(Gruber.glodap.omega[,1])==TRUE & MgCa.dat.Gruber$depth>2000)]<-NA
MgCa.dat.Gruber$GLODAP.pH.deep.cut<-rep(mean(MgCa.dat.Gruber$GLODAP.pH.deep[-cut],na.rm=T),length.out=nrow(MgCa.dat.Gruber))
MgCa.dat.Gruber$GLODAP.pH.deep.cut[cut]<-MgCa.dat.Gruber$GLODAP.pH.deep[cut]
MgCa.dat.Gruber$GLODAP.pH.deep.cut[which(is.na(Gruber.glodap.pH[,1])==TRUE & MgCa.dat.Gruber$depth>2000)]<-NA

#clean up to remove rows without environmental data
#This removes 58 records, primarily because of carbonate parameters. Maybe consider putting some back if T and S are only significant predictors?
no.data<-which(is.na(rowMeans(Gruber.WOA.s))==TRUE | is.na(rowMeans(Gruber.WOA.t))==TRUE | is.na(rowMeans(Gruber.glodap.omega.s))==TRUE | is.na(rowMeans(Gruber.glodap.pH.s))==TRUE | is.na(MgCa.dat.Gruber$GLODAP.OmegaC.deep.cut)==TRUE)
MgCa.dat.Gruber<-MgCa.dat.Gruber[-no.data,]

#compute a correlation matrix
Gruber.subset<-data.frame(s=as.vector(MgCa.dat.Gruber[,12]),																	#salinity
			t=as.vector(MgCa.dat.Gruber[,13]),																					#temperature
			OmegaC.surf=as.vector(MgCa.dat.Gruber[,14]),																			#omega surface
			pH.surf=as.vector(MgCa.dat.Gruber[,15]),																				#pH surface
			OmegaC.deep=as.vector(matrix(MgCa.dat.Gruber$GLODAP.OmegaC.deep.cut,nrow=nrow(MgCa.dat.Gruber),ncol=length(d))),	#omega deep after cutoff
			pH.deep=as.vector(matrix(MgCa.dat.Gruber$GLODAP.pH.deep.cut,nrow=nrow(MgCa.dat.Gruber),ncol=length(d)))				#pH deep after cutoff	
			)

Gruber.cor.mat<-round(cor(Gruber.subset),2)

#save environmental data and correlation matrix

write.csv(MgCa.dat.Gruber,file="Gruber.MgCa.env.dat.csv")
write.csv(Gruber.cor.mat,file="Gruber.cor.mat.csv")
Gruber95sig.cor<-2*nrow(MgCa.dat.Gruber)^-0.5

dev.new(width=6,height=8)
par(mfrow=c(3,2))
plot(Gruber.WOA.s,Gruber.WOA.t,xlab="salinity",ylab="temperature",main="G. ruber: 0-50m",las=1)
plot(Gruber.WOA.s,Gruber.glodap.omega.s,xlab="salinity",ylab="surface saturation state",las=1)
plot(Gruber.WOA.s,Gruber.glodap.pH.s,xlab="salinity",ylab="pH",las=1)
plot(Gruber.WOA.t,Gruber.glodap.omega.s,xlab="temperature",ylab="surface saturation state",las=1)
plot(Gruber.WOA.t,Gruber.glodap.pH.s,xlab="temperature",ylab="pH",las=1)
plot(Gruber.glodap.pH.s,Gruber.glodap.omega.s,xlab="pH",ylab="surface saturation state",las=1)

#***********************************************#
#coretop calibration of G.inflata
#75-150m north of 35; 150-350 south of 35º (including S. Atl). See Cleroux 07; Groenveld11.
#The deeper depth is expanded to have the same number of gridpoints as the shallower depth.
Ginflata.d1<-c(74,151)											
Ginflata.d2<-c(149,351)
d11<-which(depth.glodap<max(Ginflata.d1) & depth.glodap>min(Ginflata.d1))
d21<-which(depth.glodap<max(Ginflata.d2) & depth.glodap>min(Ginflata.d2))
d1<-match(depth.glodap[d11],depth)
d2<-match(depth.glodap[d21],depth)

Ginflata.s1<-s.spring[,,d1]
Ginflata.s1.sd<-s.spring.sd[,,d1]
Ginflata.t1<-t.spring[,,d1]
Ginflata.t1.sd<-t.spring.sd[,,d1]
Ginflata.s2<-s.spring[,,d2]
Ginflata.s2.sd<-s.spring.sd[,,d2]
Ginflata.t2<-t.spring[,,d2]
Ginflata.t2.sd<-t.spring.sd[,,d2]
Ginflata.omega1.s<-OmegaC[,,d11]
Ginflata.omega1.s.sd<-OmegaC.sd[,,d11]
Ginflata.pH1.s<-pH[,,d11]
Ginflata.pH1.s.sd<-pH.sd[,,d11]
Ginflata.omega2.s<-OmegaC[,,d21]
Ginflata.omega2.s.sd<-OmegaC.sd[,,d21]
Ginflata.pH2.s<-pH[,,d21]
Ginflata.pH2.s.sd<-pH.sd[,,d21]

Ginflata.WOA.s<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.WOA.s.sd<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.WOA.t<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.WOA.t.sd<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.glodap.omega.s<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.glodap.omega.s.sd<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.glodap.omega<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.glodap.omega.sd<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.glodap.pH.s<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.glodap.pH.s.sd<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.glodap.pH<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))
Ginflata.glodap.pH.sd<-matrix(,nrow(MgCa.dat.Ginflata),ncol=length(d1))

for(i in 1:nrow(MgCa.dat.Ginflata)) {
	#WOA T and S
	lon.tmp<-which(abs(lon-MgCa.dat.Ginflata$lon[i])==min(abs(lon-MgCa.dat.Ginflata$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Ginflata$lat[i])==min(abs(lat-MgCa.dat.Ginflata$lat[i])))[1]
	if (lat[lat.tmp]>35) {
		Ginflata.s<-Ginflata.s1
		Ginflata.t<-Ginflata.t1
		Ginflata.pH.s<-Ginflata.pH1.s
		Ginflata.omega.s<-Ginflata.omega1.s
		Ginflata.s.sd<-Ginflata.s1.sd
		Ginflata.t.sd<-Ginflata.t1.sd
		Ginflata.pH.s.sd<-Ginflata.pH1.s.sd
		Ginflata.omega.s.sd<-Ginflata.omega1.s.sd
		} else {
		Ginflata.s<-Ginflata.s2
		Ginflata.t<-Ginflata.t2
		Ginflata.pH.s<-Ginflata.pH2.s
		Ginflata.omega.s<-Ginflata.omega2.s
		Ginflata.s.sd<-Ginflata.s2.sd
		Ginflata.t.sd<-Ginflata.t2.sd
		Ginflata.pH.s.sd<-Ginflata.pH2.s.sd
		Ginflata.omega.s.sd<-Ginflata.omega2.s.sd
		}
	if (is.na(Ginflata.s[lon.tmp,lat.tmp,1])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																		#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Ginflata.WOA.s[i,]<-apply(Ginflata.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Ginflata.WOA.t[i,]<-apply(Ginflata.t[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Ginflata.WOA.s.sd[i,]<-apply(Ginflata.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Ginflata.WOA.t.sd[i,]<-apply(Ginflata.t.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
		} else {
			Ginflata.WOA.s[i,]<-Ginflata.s[lon.tmp,lat.tmp,]
			Ginflata.WOA.t[i,]<-Ginflata.t[lon.tmp,lat.tmp,]
			Ginflata.WOA.s.sd[i,]<-Ginflata.s.sd[lon.tmp,lat.tmp,]
			Ginflata.WOA.t.sd[i,]<-Ginflata.t.sd[lon.tmp,lat.tmp,]
		}
	#GLODAP Omega and pH
	lon.tmp<-which(abs(lon.glodap-MgCa.dat.Ginflata$lon[i])==min(abs(lon.glodap-MgCa.dat.Ginflata$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Ginflata$lat[i])==min(abs(lat-MgCa.dat.Ginflata$lat[i])))[1]
	tmp<-depth.glodap-MgCa.dat.Ginflata$depth[i]
	depth.tmp<-which(tmp==max(tmp[tmp<0]))[1]																					#find the glodap depth closest (without going over) to the core depth
	if (is.na(Ginflata.omega.s[lon.tmp,lat.tmp,1])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																				#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Ginflata.glodap.omega.s[i,]<-apply(Ginflata.omega.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Ginflata.glodap.omega.s.sd[i,]<-apply(Ginflata.omega.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Ginflata.glodap.pH.s[i,]<-apply(Ginflata.pH.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Ginflata.glodap.pH.s.sd[i,]<-apply(Ginflata.pH.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Ginflata.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			#if there's no carbonate data at this depth move up a bit
			if (all(is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp]))==TRUE) {
			depth.tmp<-ifelse(depth.tmp<25,depth.tmp-2,depth.tmp-1)
			Ginflata.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))		
			} 
		} else {
			Ginflata.glodap.omega.s[i,]<-Ginflata.omega.s[lon.tmp,lat.tmp,]
			Ginflata.glodap.omega.s.sd[i,]<-Ginflata.omega.s.sd[lon.tmp,lat.tmp,]
			Ginflata.glodap.pH.s[i,]<-Ginflata.pH.s[lon.tmp,lat.tmp,]
			Ginflata.glodap.pH.s.sd[i,]<-Ginflata.pH.s.sd[lon.tmp,lat.tmp,]
			Ginflata.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			if (is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp])==TRUE) {
			depth.tmp<-ifelse(depth.tmp<25,depth.tmp-2,depth.tmp-1)				
			Ginflata.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			Ginflata.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d1))
			}
		}
	}
	
MgCa.dat.Ginflata$WOA.s<-Ginflata.WOA.s
MgCa.dat.Ginflata$WOA.t<-Ginflata.WOA.t
MgCa.dat.Ginflata$GLODAP.OmegaC.surf<-Ginflata.glodap.omega.s
MgCa.dat.Ginflata$GLODAP.pH.surf<-Ginflata.glodap.pH.s
MgCa.dat.Ginflata$GLODAP.OmegaC.deep<-Ginflata.glodap.omega
MgCa.dat.Ginflata$GLODAP.pH.deep<-Ginflata.glodap.pH

MgCa.dat.Ginflata$WOA.s.sd<-Ginflata.WOA.s.sd
MgCa.dat.Ginflata$WOA.s.sd[which(is.na(rowMeans(Ginflata.WOA.s.sd))==TRUE | rowMeans(Ginflata.WOA.s.sd)==0),]<-colMeans(Ginflata.WOA.s.sd,na.rm=T)
MgCa.dat.Ginflata$WOA.t.sd<-Ginflata.WOA.t.sd
MgCa.dat.Ginflata$WOA.t.sd[which(is.na(rowMeans(Ginflata.WOA.t.sd))==TRUE | rowMeans(Ginflata.WOA.t.sd)==0),]<-colMeans(Ginflata.WOA.t.sd,na.rm=T)
MgCa.dat.Ginflata$GLODAP.OmegaC.surf.sd<-Ginflata.glodap.omega.s.sd
MgCa.dat.Ginflata$GLODAP.OmegaC.surf.sd[which(is.na(rowMeans(Ginflata.glodap.omega.s.sd))==TRUE | rowMeans(Ginflata.glodap.omega.s.sd)==0),]<-colMeans(Ginflata.glodap.omega.s.sd,na.rm=T)
MgCa.dat.Ginflata$GLODAP.pH.surf.sd<-Ginflata.glodap.pH.s.sd
MgCa.dat.Ginflata$GLODAP.pH.surf.sd[which(is.na(rowMeans(Ginflata.glodap.pH.s.sd))==TRUE | rowMeans(Ginflata.glodap.pH.s.sd)==0),]<-colMeans(Ginflata.glodap.pH.s.sd,na.rm=T)
MgCa.dat.Ginflata$GLODAP.OmegaC.deep.sd<-Ginflata.glodap.omega.sd
MgCa.dat.Ginflata$GLODAP.OmegaC.deep.sd[which(is.na(rowMeans(Ginflata.glodap.omega.sd))==TRUE | rowMeans(Ginflata.glodap.omega.sd)==0),]<-colMeans(Ginflata.glodap.omega.sd,na.rm=T)
MgCa.dat.Ginflata$GLODAP.pH.deep.sd<-Ginflata.glodap.pH.sd
MgCa.dat.Ginflata$GLODAP.pH.deep.sd[which(is.na(rowMeans(Ginflata.glodap.pH.sd))==TRUE | rowMeans(Ginflata.glodap.pH.sd)==0),]<-colMeans(Ginflata.glodap.pH.sd,na.rm=T)

#set deep values to the critical value when they exceed the critical value. 
cut<-which(Ginflata.glodap.omega[,1]<=OmegaCcut)
MgCa.dat.Ginflata$GLODAP.OmegaC.deep.cut<-rep(OmegaCcut,length.out=nrow(MgCa.dat.Ginflata))
MgCa.dat.Ginflata$GLODAP.OmegaC.deep.cut[cut]<-MgCa.dat.Ginflata$GLODAP.OmegaC.deep[cut]
MgCa.dat.Ginflata$GLODAP.OmegaC.deep.cut[which(is.na(Ginflata.glodap.omega[,1])==TRUE & MgCa.dat.Ginflata$depth>2000)]<-NA
MgCa.dat.Ginflata$GLODAP.pH.deep.cut<-rep(mean(MgCa.dat.Ginflata$GLODAP.pH.deep[-cut],na.rm=T),length.out=nrow(MgCa.dat.Ginflata))
MgCa.dat.Ginflata$GLODAP.pH.deep.cut[cut]<-MgCa.dat.Ginflata$GLODAP.pH.deep[cut]
MgCa.dat.Ginflata$GLODAP.pH.deep.cut[which(is.na(Ginflata.glodap.pH[,1])==TRUE & MgCa.dat.Ginflata$depth>2000)]<-NA

#clean up to remove rows without environmental data
#This removes 21 records
no.data<-which(is.na(rowMeans(Ginflata.WOA.s))==TRUE | is.na(rowMeans(Ginflata.WOA.t))==TRUE | is.na(rowMeans(Ginflata.glodap.omega.s))==TRUE | is.na(rowMeans(Ginflata.glodap.pH.s))==TRUE | is.na(MgCa.dat.Ginflata$GLODAP.OmegaC.deep.cut)==TRUE)
MgCa.dat.Ginflata<-MgCa.dat.Ginflata[-no.data,]

#compute a correlation matrix
Ginflata.subset<-data.frame(s=as.vector(MgCa.dat.Ginflata[,12]),																	#salinity
			t=as.vector(MgCa.dat.Ginflata[,13]),																					#temperature
			OmegaC.surf=as.vector(MgCa.dat.Ginflata[,14]),																			#omega surface
			pH.surf=as.vector(MgCa.dat.Ginflata[,15]),																				#pH surface
			OmegaC.deep=as.vector(matrix(MgCa.dat.Ginflata$GLODAP.OmegaC.deep.cut,nrow=nrow(MgCa.dat.Ginflata),ncol=length(d1))),	#omega deep after cutoff
			pH.deep=as.vector(matrix(MgCa.dat.Ginflata$GLODAP.pH.deep.cut,nrow=nrow(MgCa.dat.Ginflata),ncol=length(d1)))				#pH deep after cutoff	
			)

Ginflata.cor.mat<-round(cor(Ginflata.subset),2)

#save environmental data and correlation matrix

write.csv(MgCa.dat.Ginflata,file="Ginflata.MgCa.env.dat.csv")
write.csv(Ginflata.cor.mat,file="Ginflata.cor.mat.csv")
Ginflata95sig.cor<-2*nrow(MgCa.dat.Ginflata)^-0.5

dev.new(width=6,height=8)
par(mfrow=c(3,2))
plot(Ginflata.WOA.s,Ginflata.WOA.t,xlab="salinity",ylab="temperature",main="G. inflata",las=1)
plot(Ginflata.WOA.s,Ginflata.glodap.omega.s,xlab="salinity",ylab="surface saturation state",las=1,main="75-150m lat>35ºN; 150-350 lat<35ºN")
plot(Ginflata.WOA.s,Ginflata.glodap.pH.s,xlab="salinity",ylab="pH",las=1)
plot(Ginflata.WOA.t,Ginflata.glodap.omega.s,xlab="temperature",ylab="surface saturation state",las=1)
plot(Ginflata.WOA.t,Ginflata.glodap.pH.s,xlab="temperature",ylab="pH",las=1)
plot(Ginflata.glodap.pH.s,Ginflata.glodap.omega.s,xlab="pH",ylab="surface saturation state",las=1)

#***********************************************#
#coretop calibration of G.bulloides
Gbulloides.d<-c(0,100)
d1<-which(depth.glodap<=max(Gbulloides.d) & depth.glodap>=min(Gbulloides.d))
d<-match(depth.glodap[d1],depth)
Gbulloides.ta<-t.ann[,,d]
Gbulloides.ta.sd<-t.ann.sd[,,d]
Gbulloides.ta.mean<-apply(Gbulloides.ta,c(1,2),mean)
Gbulloides.t1<-t.summer[,,d]
Gbulloides.t1.sd<-t.summer.sd[,,d]
Gbulloides.t2<-t.spring[,,d]
Gbulloides.t2.sd<-t.spring.sd[,,d]
Gbulloides.s1<-s.summer[,,d]
Gbulloides.s1.sd<-s.summer.sd[,,d]
Gbulloides.s2<-s.spring[,,d]
Gbulloides.s2.sd<-s.spring.sd[,,d]

#Start with spring, but if mean annual T<=10, replace with summer based on Jonkers 15
Gbulloides.t<-Gbulloides.t2
Gbulloides.t.sd<-Gbulloides.t2.sd
tmp<-which(Gbulloides.ta.mean<=10)																			#indicies of gridboxes with mean T<10ºC
Gbulloides.t[tmp]<-Gbulloides.t1[tmp]
Gbulloides.t.sd[tmp]<-Gbulloides.t1.sd[tmp]
Gbulloides.s<-Gbulloides.s2
Gbulloides.s.sd<-Gbulloides.s2.sd
Gbulloides.s[tmp]<-Gbulloides.s1[tmp]
Gbulloides.s.sd[tmp]<-Gbulloides.s1.sd[tmp]

Gbulloides.omega.s<-OmegaC[,,d1]
Gbulloides.omega.s.sd<-OmegaC.sd[,,d1]
Gbulloides.pH.s<-pH[,,d1]
Gbulloides.pH.s.sd<-pH.sd[,,d1]

Gbulloides.WOA.s<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.WOA.s.sd<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.WOA.t<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.WOA.t.sd<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.glodap.omega.s<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.glodap.omega.s.sd<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.glodap.omega<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.glodap.omega.sd<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.glodap.pH.s<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.glodap.pH.s.sd<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.glodap.pH<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))
Gbulloides.glodap.pH.sd<-matrix(,nrow(MgCa.dat.Gbulloides),ncol=length(d))

for(i in 1:nrow(MgCa.dat.Gbulloides)) {
	#WOA T and S
	lon.tmp<-which(abs(lon-MgCa.dat.Gbulloides$lon[i])==min(abs(lon-MgCa.dat.Gbulloides$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Gbulloides$lat[i])==min(abs(lat-MgCa.dat.Gbulloides$lat[i])))[1]
	if (is.na(Gbulloides.s[lon.tmp,lat.tmp,1])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																		#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Gbulloides.WOA.s[i,]<-apply(Gbulloides.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gbulloides.WOA.t[i,]<-apply(Gbulloides.t[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gbulloides.WOA.s.sd[i,]<-apply(Gbulloides.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gbulloides.WOA.t.sd[i,]<-apply(Gbulloides.t.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
		} else {
			Gbulloides.WOA.s[i,]<-Gbulloides.s[lon.tmp,lat.tmp,]
			Gbulloides.WOA.t[i,]<-Gbulloides.t[lon.tmp,lat.tmp,]
			Gbulloides.WOA.s.sd[i,]<-Gbulloides.s.sd[lon.tmp,lat.tmp,]
			Gbulloides.WOA.t.sd[i,]<-Gbulloides.t.sd[lon.tmp,lat.tmp,]
		}
	
	#GLODAP Omega and pH
	lon.tmp<-which(abs(lon.glodap-MgCa.dat.Gbulloides$lon[i])==min(abs(lon.glodap-MgCa.dat.Gbulloides$lon[i])))[1]
	lat.tmp<-which(abs(lat-MgCa.dat.Gbulloides$lat[i])==min(abs(lat-MgCa.dat.Gbulloides$lat[i])))[1]
	tmp<-depth.glodap-MgCa.dat.Gbulloides$depth[i]
	depth.tmp<-which(tmp==max(tmp[tmp<0]))[1]																					#find the glodap depth closest (without going over) to the core depth
	if (is.na(Gbulloides.omega.s[lon.tmp,lat.tmp,1])==TRUE) {																	#if there's no model data at the proxy site
			lon.tmp<-seq(lon.tmp-1,lon.tmp+1)																				#take the surrounding gridboxes
			lat.tmp<-seq(lat.tmp-1,lat.tmp+1)
			Gbulloides.glodap.omega.s[i,]<-apply(Gbulloides.omega.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gbulloides.glodap.omega.s.sd[i,]<-apply(Gbulloides.omega.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gbulloides.glodap.pH.s[i,]<-apply(Gbulloides.pH.s[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gbulloides.glodap.pH.s.sd[i,]<-apply(Gbulloides.pH.s.sd[lon.tmp,lat.tmp,],3,mean,na.rm=TRUE)
			Gbulloides.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			#if there's no carbonate data at this depth move up a bit
			if (all(is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp]))==TRUE) {
			depth.tmp<-ifelse(depth.tmp<25,depth.tmp-2,depth.tmp-1)
			Gbulloides.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))			
			} 
		} else {
			Gbulloides.glodap.omega.s[i,]<-Gbulloides.omega.s[lon.tmp,lat.tmp,]
			Gbulloides.glodap.omega.s.sd[i,]<-Gbulloides.omega.s.sd[lon.tmp,lat.tmp,]
			Gbulloides.glodap.pH.s[i,]<-Gbulloides.pH.s[lon.tmp,lat.tmp,]
			Gbulloides.glodap.pH.s.sd[i,]<-Gbulloides.pH.s.sd[lon.tmp,lat.tmp,]
			Gbulloides.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			if (is.na(OmegaC[lon.tmp,lat.tmp,depth.tmp])==TRUE) {
			depth.tmp<-ifelse(depth.tmp<25,depth.tmp-2,depth.tmp-1)				
			Gbulloides.glodap.omega[i,]<-rep(mean(OmegaC[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.pH[i,]<-rep(mean(pH[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.omega.sd[i,]<-rep(mean(OmegaC.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			Gbulloides.glodap.pH.sd[i,]<-rep(mean(pH.sd[lon.tmp,lat.tmp,depth.tmp],na.rm=TRUE),length(d))
			}
		}
	}

MgCa.dat.Gbulloides$WOA.s<-Gbulloides.WOA.s
MgCa.dat.Gbulloides$WOA.t<-Gbulloides.WOA.t
MgCa.dat.Gbulloides$GLODAP.OmegaC.surf<-Gbulloides.glodap.omega.s
MgCa.dat.Gbulloides$GLODAP.pH.surf<-Gbulloides.glodap.pH.s
MgCa.dat.Gbulloides$GLODAP.OmegaC.deep<-Gbulloides.glodap.omega
MgCa.dat.Gbulloides$GLODAP.pH.deep<-Gbulloides.glodap.pH

MgCa.dat.Gbulloides$WOA.s.sd<-Gbulloides.WOA.s.sd
MgCa.dat.Gbulloides$WOA.s.sd[which(is.na(rowMeans(Gbulloides.WOA.s.sd))==TRUE | rowMeans(Gbulloides.WOA.s.sd)==0),]<-colMeans(Gbulloides.WOA.s.sd,na.rm=T)
MgCa.dat.Gbulloides$WOA.t.sd<-Gbulloides.WOA.t.sd
MgCa.dat.Gbulloides$WOA.t.sd[which(is.na(rowMeans(Gbulloides.WOA.t.sd))==TRUE | rowMeans(Gbulloides.WOA.t.sd)==0),]<-colMeans(Gbulloides.WOA.t.sd,na.rm=T)
MgCa.dat.Gbulloides$GLODAP.OmegaC.surf.sd<-Gbulloides.glodap.omega.s.sd
MgCa.dat.Gbulloides$GLODAP.OmegaC.surf.sd[which(is.na(rowMeans(Gbulloides.glodap.omega.s.sd))==TRUE | rowMeans(Gbulloides.glodap.omega.s.sd)==0),]<-colMeans(Gbulloides.glodap.omega.s.sd,na.rm=T)
MgCa.dat.Gbulloides$GLODAP.pH.surf.sd<-Gbulloides.glodap.pH.s.sd
MgCa.dat.Gbulloides$GLODAP.pH.surf.sd[which(is.na(rowMeans(Gbulloides.glodap.pH.s.sd))==TRUE | rowMeans(Gbulloides.glodap.pH.s.sd)==0),]<-colMeans(Gbulloides.glodap.pH.s.sd,na.rm=T)
MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.sd<-Gbulloides.glodap.omega.sd
MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.sd[which(is.na(rowMeans(Gbulloides.glodap.omega.sd))==TRUE | rowMeans(Gbulloides.glodap.omega.sd)==0),]<-colMeans(Gbulloides.glodap.omega.sd,na.rm=T)
MgCa.dat.Gbulloides$GLODAP.pH.deep.sd<-Gbulloides.glodap.pH.sd
MgCa.dat.Gbulloides$GLODAP.pH.deep.sd[which(is.na(rowMeans(Gbulloides.glodap.pH.sd))==TRUE | rowMeans(Gbulloides.glodap.pH.sd)==0),]<-colMeans(Gbulloides.glodap.pH.sd,na.rm=T)

#set deep values to the critical value when they exceed the critical value. 
cut<-which(Gbulloides.glodap.omega[,1]<=OmegaCcut)
MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.cut<-rep(OmegaCcut,length.out=nrow(MgCa.dat.Gbulloides))
MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.cut[cut]<-MgCa.dat.Gbulloides$GLODAP.OmegaC.deep[cut]
MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.cut[which(is.na(Gbulloides.glodap.omega[,1])==TRUE & MgCa.dat.Gbulloides$depth>2000)]<-NA
MgCa.dat.Gbulloides$GLODAP.pH.deep.cut<-rep(mean(MgCa.dat.Gbulloides$GLODAP.pH.deep[-cut],na.rm=T),length.out=nrow(MgCa.dat.Gbulloides))
MgCa.dat.Gbulloides$GLODAP.pH.deep.cut[cut]<-MgCa.dat.Gbulloides$GLODAP.pH.deep[cut]
MgCa.dat.Gbulloides$GLODAP.pH.deep.cut[which(is.na(Gbulloides.glodap.pH[,1])==TRUE & MgCa.dat.Gbulloides$depth>2000)]<-NA

#clean up to remove rows without environmental data
#This removes 17 records
no.data<-which(is.na(rowMeans(Gbulloides.WOA.s))==TRUE | is.na(rowMeans(Gbulloides.WOA.t))==TRUE | is.na(rowMeans(Gbulloides.glodap.omega.s))==TRUE | is.na(rowMeans(Gbulloides.glodap.pH.s))==TRUE | is.na(MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.cut)==TRUE)
MgCa.dat.Gbulloides<-MgCa.dat.Gbulloides[-no.data,]

#compute a correlation matrix
Gbulloides.subset<-data.frame(s=as.vector(MgCa.dat.Gbulloides[,12]),																	#salinity
			t=as.vector(MgCa.dat.Gbulloides[,13]),																					#temperature
			OmegaC.surf=as.vector(MgCa.dat.Gbulloides[,14]),																			#omega surface
			pH.surf=as.vector(MgCa.dat.Gbulloides[,15]),																				#pH surface
			OmegaC.deep=as.vector(matrix(MgCa.dat.Gbulloides$GLODAP.OmegaC.deep.cut,nrow=nrow(MgCa.dat.Gbulloides),ncol=length(d))),	#omega deep after cutoff
			pH.deep=as.vector(matrix(MgCa.dat.Gbulloides$GLODAP.pH.deep.cut,nrow=nrow(MgCa.dat.Gbulloides),ncol=length(d)))				#pH deep after cutoff	
			)

Gbulloides.cor.mat<-round(cor(Gbulloides.subset),2)

#save environmental data and correlation matrix

write.csv(MgCa.dat.Gbulloides,file="Gbulloides.MgCa.env.dat.csv")
write.csv(Gbulloides.cor.mat,file="Gbulloides.cor.mat.csv")
Gbulloides95sig.cor<-2*nrow(MgCa.dat.Gbulloides)^-0.5

dev.new(width=6,height=8)
par(mfrow=c(3,2))
plot(Gbulloides.WOA.s,Gbulloides.WOA.t,xlab="salinity",ylab="temperature",main="G.bulloides: 0-100m",las=1)
plot(Gbulloides.WOA.s,Gbulloides.glodap.omega.s,xlab="salinity",ylab="surface saturation state",las=1)
plot(Gbulloides.WOA.s,Gbulloides.glodap.pH.s,xlab="salinity",ylab="pH",las=1)
plot(Gbulloides.WOA.t,Gbulloides.glodap.omega.s,xlab="temperature",ylab="surface saturation state",las=1)
plot(Gbulloides.WOA.t,Gbulloides.glodap.pH.s,xlab="temperature",ylab="pH",las=1)
plot(Gbulloides.glodap.pH.s,Gbulloides.glodap.omega.s,xlab="pH",ylab="surface saturation state",las=1)

#***********************************************#
#Compare d18O-based temperatures to assigned temperatures
dev.new(width=6,height=6)
par(mfrow=c(2,2))
#subset Npachy data
MgCa.dat.NpachyR<-MgCa.dat.Npachy[which(MgCa.dat.Npachy$Notes=="Dextral"),]
MgCa.dat.NpachyL<-MgCa.dat.Npachy[which(MgCa.dat.Npachy$Notes=="Sinstral"),]
plot(rowMeans(MgCa.dat.NpachyR$WOA.t),MgCa.dat.NpachyR$O18T,xlab="T (ºC, this study)",ylab="T (ºC, d18O-based)",main="N. pachyderma & N. incompta",las=1,xlim=c(-5,20),ylim=c(-5,20),col="grey50")
points(rowMeans(MgCa.dat.NpachyL$WOA.t),MgCa.dat.NpachyL$O18T,col="palegreen3")
legend("topleft",legend=c("N. pachyderma","N. incompta"),pch=c(1,1),col=c("palegreen3","grey50"),cex=0.8,bty="n")
abline(0,1)
plot(rowMeans(MgCa.dat.Gruber$WOA.t),MgCa.dat.Gruber$O18T,xlab="T (ºC, this study)",ylab="T (ºC, d18O-based)",main="G. ruber",las=1,col="cornflowerblue")
abline(0,1)
plot(rowMeans(MgCa.dat.Ginflata$WOA.t),MgCa.dat.Ginflata$O18T,xlab="T (ºC, this study)",ylab="T (ºC, d18O-based)",main="G. inflata",las=1,col="firebrick")
abline(0,1)
plot(rowMeans(MgCa.dat.Gbulloides$WOA.t),MgCa.dat.Gbulloides$O18T,xlab="T (ºC, this study)",ylab="T (ºC, d18O-based)",main="G. bulloides",las=1,col="gold3")
abline(0,1)
