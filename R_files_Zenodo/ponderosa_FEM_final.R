
library(DHARMa)
library(car)
library(performance)
library(visreg)
library(lme4)
library(glmmTMB)
library(MuMIn)
library(data.table)
library(psych)
library(rgdal)
library(raster)
library(ggcorrplot)
library(car)
library(patchwork)
library(MuMIn)
library(ggeffects)
library(sjPlot)
library(viridis)
library(partR2)

library(sf)

library(PNWColors)
library(nationalparkcolors)

library(scico)
library(patchwork)
library(ggplot2)

library(boot)

stdErr <- function(x) sqrt(var(x, na.rm = T)/length(na.exclude(x)))


setwd("C:/Users/awion/Documents/PhData/PhD Data/Cones/PIPO/DATA/Core_ecoapps/")

climate<-read.csv("clim_df_FULL_min.csv")
climate<-as.data.table(climate)
climate<-climate[,-1]

cones<-read.csv("full_cones_3.csv")
cones<-as.data.table(cones)
cones<-cones[!is.na(cones$cone_total),]
cones$std<-cones$cone_total/cones$leadershoots



CWD<-raster('CWD.tif')
age<-read.csv("tree_age.csv")
AET<-raster('AET.tif')

sites<-readOGR("ponderosa_sites_10_27.shp")
names<-sites$site

pipo<-raster("PIPO_1km.tif")

bound<-readOGR("cb_2018_us_state_20m.shp")
bound<-spTransform(bound,CRS=crs(pipo))
bound$NAME
#map stuff
fullwest<-bound[c(6,11,24,28,39,43,44,45,46,51),]
west<-bound[c(11,24,28,39,43,44,46,51),]
west_df<-fortify(west)
fullwest_df<-fortify(fullwest)


lat<-crs(CWD)
aea<-crs(pipo)

CWD_AEA<-projectRaster(CWD,crs=aea)
AET_AEA<-projectRaster(AET,crs=aea)

western_pipo<-crop(pipo,west)
plot(western_pipo)
western_pipo[western_pipo>1]<-1
western_pipo[western_pipo<1]<-NA
plot(western_pipo)
extent(western_pipo)
extent(CWD_AEA)

western_pipo<-crop(western_pipo,CWD_AEA)
plot(western_pipo)

pipo_cwd_pr<-projectRaster(western_pipo,CWD_AEA)
pipo_aet_pr<-projectRaster(western_pipo,AET_AEA)

pipo_cwdmask<-mask(CWD_AEA,pipo_cwd_pr)
pipo_aetmask<-mask(AET_AEA,pipo_aet_pr)

aea_sites<-spTransform(sites,CRS=crs(pipo))
aea_sites
plot(aea_sites)
plot(pipo)
s_df<-as.data.frame(aea_sites,xy=TRUE)
p_df<-as.data.frame(pipo,xy=TRUE)
buff_crop<-c(-1500000,-500000,1350000,2500000)

pipo_crop<-crop(pipo,west)
pc_df<-as.data.frame(pipo_crop,xy=TRUE)
pc_df$layer<-ifelse(pc_df$PIPO_1km==0,NA,pc_df$PIPO_1km)
extent(pipo_crop)
plot(pipo_crop)


#setting up the climate data
#drop tdmean bc I only use it to calc VPD
#melt it
climate_stacked<-melt(climate, id.vars=c("year","site","month"))
levels(climate_stacked$variable)

#combine climate variable and month together
climate_stacked$monthvar<-paste(climate_stacked$variable,climate_stacked$month,sep="")
#take out most recent years of climate data
#climate_stacked<-subset(climate_stacked,year<2017)
tail(climate_stacked)

#ghosts of climate past
#important to note here that the years have not been adjusted from the raw scar data. in otherwords, year 0 (and year column in cone data) corresponds to year of cone pollination, or the year that the cone scar on the branch is present on

climate_minus<-climate_stacked
climate_minus$year<-climate_minus$year
climate_minus$monthvar<-paste("yr0",climate_minus$monthvar,sep="")
#climate_minus1<-subset(climate_minus1,year>1999)

climate_minus1<-climate_stacked
climate_minus1$year<-climate_minus1$year+1
climate_minus1$monthvar<-paste("yr1",climate_minus1$monthvar,sep="")
#climate_minus1<-subset(climate_minus1,year>1999)

climate_minus2<-climate_stacked
climate_minus2$year<-climate_stacked$year+2
climate_minus2$monthvar<-paste("yr2",climate_minus2$monthvar,sep="")
#climate_minus2<-subset(climate_minus2,year>1999)

#combine
climate_stacked_allyrs<-rbind(climate_minus,climate_minus1,climate_minus2)
#select relevant years
climate_stacked_allyrs<-subset(climate_stacked_allyrs,year>1999)
climate_stacked_allyrs<-subset(climate_stacked_allyrs,year<2020)
#cast it
climate_unstacked <- dcast(climate_stacked_allyrs,site+year~monthvar,measure.var=c("value"))
#return to climate name
climate<-climate_unstacked


season<-climate[,]

#####3 month

season<-climate[,yr2pptjfm:=(yr2ppt1+yr2ppt2+yr2ppt3)]
season<-climate[,yr2pptfma:=(yr2ppt2+yr2ppt3+yr2ppt4)]
season<-climate[,yr2pptmam:=(yr2ppt3+yr2ppt4+yr2ppt5)]
season<-climate[,yr2pptamj:=(yr2ppt4+yr2ppt5+yr2ppt6)]
season<-climate[,yr2pptmjj:=(yr2ppt5+yr2ppt6+yr2ppt7)]
season<-climate[,yr2pptjja:=(yr2ppt6+yr2ppt7+yr2ppt8)]
season<-climate[,yr2pptjas:=(yr2ppt7+yr2ppt8+yr2ppt9)]
season<-climate[,yr2pptaso:=(yr2ppt8+yr2ppt9+yr2ppt10)]
season<-climate[,yr2pptson:=(yr2ppt9+yr2ppt10+yr2ppt11)]
season<-climate[,yr2pptond:=(yr2ppt10+yr2ppt11+yr2ppt12)]
season<-climate[,yr2pptndj:=(yr2ppt11+yr2ppt12+yr1ppt1)]
season<-climate[,yr2pptdjf:=(yr2ppt12+yr1ppt1+yr1ppt2)]


season<-climate[,yr1pptjfm:=(yr1ppt1+yr1ppt2+yr1ppt3)]
season<-climate[,yr1pptfma:=(yr1ppt2+yr1ppt3+yr1ppt4)]
season<-climate[,yr1pptmam:=(yr1ppt3+yr1ppt4+yr1ppt5)]
season<-climate[,yr1pptamj:=(yr1ppt4+yr1ppt5+yr1ppt6)]
season<-climate[,yr1pptmjj:=(yr1ppt5+yr1ppt6+yr1ppt7)]
season<-climate[,yr1pptjja:=(yr1ppt6+yr1ppt7+yr1ppt8)]
season<-climate[,yr1pptjas:=(yr1ppt7+yr1ppt8+yr1ppt9)]
season<-climate[,yr1pptaso:=(yr1ppt8+yr1ppt9+yr1ppt10)]
season<-climate[,yr1pptson:=(yr1ppt9+yr1ppt10+yr1ppt11)]
season<-climate[,yr1pptond:=(yr1ppt10+yr1ppt11+yr1ppt12)]
season<-climate[,yr1pptndj:=(yr1ppt11+yr1ppt12+yr1ppt1)]
season<-climate[,yr1pptdjf:=(yr1ppt12+yr0ppt1+yr0ppt2)]

season<-climate[,yr0pptjfm:=(yr0ppt1+yr0ppt2+yr0ppt3)]
season<-climate[,yr0pptfma:=(yr0ppt2+yr0ppt3+yr0ppt4)]
season<-climate[,yr0pptmam:=(yr0ppt3+yr0ppt4+yr0ppt5)]
season<-climate[,yr0pptamj:=(yr0ppt4+yr0ppt5+yr0ppt6)]
season<-climate[,yr0pptmjj:=(yr0ppt5+yr0ppt6+yr0ppt7)]
season<-climate[,yr0pptjja:=(yr0ppt6+yr0ppt7+yr0ppt8)]
season<-climate[,yr0pptjas:=(yr0ppt7+yr0ppt8+yr0ppt9)]
season<-climate[,yr0pptaso:=(yr0ppt8+yr0ppt9+yr0ppt10)]
season<-climate[,yr0pptson:=(yr0ppt9+yr0ppt10+yr0ppt11)]
season<-climate[,yr0pptond:=(yr0ppt10+yr0ppt11+yr0ppt12)]



season<-climate[,yr2vpdjfm:=(yr2vpd1+yr2vpd2+yr2vpd3)/3]
season<-climate[,yr2vpdfma:=(yr2vpd2+yr2vpd3+yr2vpd4)/3]
season<-climate[,yr2vpdmam:=(yr2vpd3+yr2vpd4+yr2vpd5)/3]
season<-climate[,yr2vpdamj:=(yr2vpd4+yr2vpd5+yr2vpd6)/3]
season<-climate[,yr2vpdmjj:=(yr2vpd5+yr2vpd6+yr2vpd7)/3]
season<-climate[,yr2vpdjja:=(yr2vpd6+yr2vpd7+yr2vpd8)/3]
season<-climate[,yr2vpdjas:=(yr2vpd7+yr2vpd8+yr2vpd9)/3]
season<-climate[,yr2vpdaso:=(yr2vpd8+yr2vpd9+yr2vpd10)/3]
season<-climate[,yr2vpdson:=(yr2vpd9+yr2vpd10+yr2vpd11)/3]
season<-climate[,yr2vpdond:=(yr2vpd10+yr2vpd11+yr2vpd12)/3]
season<-climate[,yr2vpdndj:=(yr2vpd11+yr2vpd12+yr1vpd1)/3]
season<-climate[,yr2vpddjf:=(yr2vpd12+yr1vpd1+yr1vpd2)/3]


season<-climate[,yr1vpdjfm:=(yr1vpd1+yr1vpd2+yr1vpd3)/3]
season<-climate[,yr1vpdfma:=(yr1vpd2+yr1vpd3+yr1vpd4)/3]
season<-climate[,yr1vpdmam:=(yr1vpd3+yr1vpd4+yr1vpd5)/3]
season<-climate[,yr1vpdamj:=(yr1vpd4+yr1vpd5+yr1vpd6)/3]
season<-climate[,yr1vpdmjj:=(yr1vpd5+yr1vpd6+yr1vpd7)/3]
season<-climate[,yr1vpdjja:=(yr1vpd6+yr1vpd7+yr1vpd8)/3]
season<-climate[,yr1vpdjas:=(yr1vpd7+yr1vpd8+yr1vpd9)/3]
season<-climate[,yr1vpdaso:=(yr1vpd8+yr1vpd9+yr1vpd10)/3]
season<-climate[,yr1vpdson:=(yr1vpd9+yr1vpd10+yr1vpd11)/3]
season<-climate[,yr1vpdond:=(yr1vpd10+yr1vpd11+yr1vpd12)/3]
season<-climate[,yr1vpdndj:=(yr1vpd11+yr1vpd12+yr0vpd1)/3]
season<-climate[,yr1vpddjf:=(yr1vpd12+yr0vpd1+yr0vpd2)/3]

season<-climate[,yr0vpdjfm:=(yr0vpd1+yr0vpd2+yr0vpd3)/3]
season<-climate[,yr0vpdfma:=(yr0vpd2+yr0vpd3+yr0vpd4)/3]
season<-climate[,yr0vpdmam:=(yr0vpd3+yr0vpd4+yr0vpd5)/3]
season<-climate[,yr0vpdamj:=(yr0vpd4+yr0vpd5+yr0vpd6)/3]
season<-climate[,yr0vpdmjj:=(yr0vpd5+yr0vpd6+yr0vpd7)/3]
season<-climate[,yr0vpdjja:=(yr0vpd6+yr0vpd7+yr0vpd8)/3]
season<-climate[,yr0vpdjas:=(yr0vpd7+yr0vpd8+yr0vpd9)/3]
season<-climate[,yr0vpdaso:=(yr0vpd8+yr0vpd9+yr0vpd10)/3]
season<-climate[,yr0vpdson:=(yr0vpd9+yr0vpd10+yr0vpd11)/3]
season<-climate[,yr0vpdond:=(yr0vpd10+yr0vpd11+yr0vpd12)/3]

season<-climate[,yr1tmeanjja:=(yr1tmean6+yr1tmean7+yr1tmean8)/3]
season<-climate[,yr2tmeanjja:=(yr2tmean6+yr2tmean7+yr2tmean8)/3]


#####

#season<-climate[,yr0vpd:=(yr0vpd5+yr0vpd6+yr0vpd7+yr0vpd8+yr0vpd9)]
#season<-climate[,yr1vpd:=(yr1vpd5+yr1vpd6+yr1vpd7+yr1vpd8+yr1vpd9)]
#season<-climate[,yr2vpd:=(yr2vpd5+yr2vpd6+yr2vpd7+yr2vpd8+yr2vpd9)]
#
#season<-climate[,yr0ppt:=(yr0ppt5+yr0ppt6+yr0ppt7+yr0ppt8+yr0ppt9)]
#season<-climate[,yr1ppt:=(yr1ppt5+yr1ppt6+yr1ppt7+yr1ppt8+yr1ppt9)]
#season<-climate[,yr2ppt:=(yr2ppt5+yr2ppt6+yr2ppt7+yr2ppt8+yr2ppt9)]

#season<-climate[,deltaP:=(yr1pptja-yr2pptja)]
season<-climate[,deltaT:=(yr1tmeanjja-yr2tmeanjja)]
season<-climate[,deltaV:=(yr1vpdjja-yr2vpdjja)]
season<-climate[,deltaP:=(yr1pptjja-yr2pptjja)]


climatetable1<-melt(climate,id.vars=c("site","year"),na.rm=TRUE)
climatetablesdmean1<-climatetable1[,list(std=sd(value),mean=mean(value)),by=list(site,variable)]
joined<-merge(climatetable1,climatetablesdmean1,by=c("site","variable"))
joined[,zscore:=(value-mean)/std]

climate_all<-melt(joined, id=c("year","site","variable"))
CLIMATE<-subset(climate_all, variable.1=="zscore")
cast<-dcast(CLIMATE, site+year~variable, fun.aggregate = mean)

z_season<-cast

names(season)
ssn<-z_season[,-c(3:218)]
names(ssn)
#merge

#extract long term climate variables 

new_aet<-extract(AET,sites)
new_aet<-cbind.data.frame(new_aet,names)
names(new_aet)<-c("AET","site")

new_cwd<-extract(CWD,sites)
new_cwd<-cbind.data.frame(new_cwd,names)
names(new_cwd)<-c("CWD","site")


#the great merge

full_cones<-merge(cones,ssn,by=c("site","year"))

full_cones<-merge(full_cones,new_cwd,by="site")

full_cones<-merge(full_cones,new_aet,by="site")


#creating a column of previous year's cone production

treeyr_cones<-full_cones[,.(cones=mean(cone_total)),by=c("tree","year")]
test<-treeyr_cones
test$year2<-test$year+1
names(test)<-c("tree","year2","AR1","year")

AR<-test

ar<-AR[,c(1,3,4)]

testin<-merge(ar,full_cones,by=c("tree","year"))

full_cones<-testin

full_cones<-as.data.table(full_cones)




#removing stray NA's

comp_cones<-full_cones[!is.na(full_cones$cone_total),]
comp_cones<-comp_cones[!is.na(comp_cones$BA5),]


ggplot


#aggregating data

tree_cones<-comp_cones[,.(mean_cones=mean(cone_total,na.rm=TRUE),CVi=(sd(cone_total,na.rm=TRUE))/(mean(cone_total,na.rn=TRUE)),
                          CVi=(sd(cone_total,na.rm=TRUE))/(mean(cone_total,na.rn=TRUE)),
                          CWD=mean(CWD),AET=mean(AET),lat=mean(lat),long=mean(long)),
                       by=c("tree","site","leadershoots","dbh","cbh","BA5","BA20","X7m_count","crown","height")]



treeyr_cones<-comp_cones[,.(cone_total=mean(cone_total),AR1=mean(AR1),leadershoots=mean(leadershoots)),by=c("tree","year","site","CWD","AET","BA5","dbh")]
tree_cones$BA5<-(tree_cones$BA5/4.356)

#tree_cones$size<-(tree_cones$height-tree_cones$cbh)*(tree_cones$crown)



site_cones<-comp_cones[,.(mean_cones=mean(cone_total,na.rm=TRUE),CVp=(sd(cone_total,na.rm=TRUE))/(mean(cone_total,na.rn=TRUE)),
                          BA5=mean(BA5,na.rm=TRUE),BA20=mean(BA20),height=mean(height),dbh=mean(dbh),lat=mean(lat),lon=mean(long),
                          CWD=mean(CWD),AET=mean(AET),leadershoots=mean(leadershoots)),by=c("site")]


siteyr_cones<-comp_cones[,.(cone_total=mean(cone_total),branchcones=mean(cone_total),AR1=mean(AR1),AET=mean(AET),CWD=mean(CWD),
                            BA5=mean(BA5,na.rm=TRUE),dbh=mean(dbh),leadershoots=mean(leadershoots)),by=c("site","year")]



t_propzero<-treeyr_cones[,.(prop=(sum(cone_total==0)/.N)),by=c("tree")]

s_propzero<-siteyr_cones[,.(prop=(sum(cone_total==0)/.N)),by=c("site")]


#site_lats<-site_cones[,c(1,8)]
#write.csv(site_lats,"lat.csv")


###tree level synchrony 

dat<-vector()
dat2<-vector()
full_data<-comp_cones[,.(cone_total=mean(cone_total)),by=c("tree","site","year")]
seriesuniq<-unique(full_data$site)

#cone totals are heavy tailed at each site, using spearman cors below

for(i in 1:length(seriesuniq)){
  sub<-subset(full_data,site==seriesuniq[i])
  tcast<-dcast(sub,year~tree,value.var = "cone_total")
  tcast<-tcast[,year:=NULL]
  trial<-(as.matrix(tcast))
  ctrial<-cor(trial,method='pearson',use="complete.obs")
  diag(ctrial)<-NA
  cmean<-as.matrix(mean(ctrial,na.rm=TRUE))
  dat<-rbind(dat,cmean)
  ctrial<-data.table(ctrial)
  rowm<-ctrial[,.(tree_s=(colMeans(ctrial,na.rm=TRUE)),by=c(colnames(ctrial)))]
  rowm<-as.matrix(rowm)
  dat2<-rbind(dat2,rowm)
}
#View(dat)
#View(dat2)


dat<-as.numeric(dat)

dat2.1<-as.numeric(dat2[,1])

(seriesuniq)

site_r<-cbind.data.frame((seriesuniq),dat)
names(site_r)<-c("site","site_r")

tree_r<-cbind.data.frame(dat2[,2],dat2.1)
names(tree_r)<-c("tree","tree_r")







tree_data<-merge(tree_r,tree_cones,by=c("tree"))
site_data<-merge(site_r,site_cones,by="site")

tree_data<-merge(tree_data,t_propzero,by="tree")
site_data<-merge(site_data,s_propzero,by="site")

treeyr_data<-merge(treeyr_cones,ssn,by=c("site","year"))
siteyr_data<-merge(siteyr_cones,ssn,by=c("site","year"))

siteyr_data$log_cones<-(siteyr_data$cone_total+1)
treeyr_data$log_cones<-(treeyr_data$cone_total+1)

siteyr_data$cone_total<-round(siteyr_data$cone_total,0)
treeyr_data$cone_total<-round(treeyr_data$cone_total,0)


#treeyr_data$masting_index<-scale((treeyr_data$deltaV7)*(0.29)+(treeyr_data$yr1pptja*(0.32)))
#siteyr_data$masting_index<-scale((siteyr_data$deltaV7)*(0.29)+(siteyr_data$yr1pptja*(0.32)))
d_matrix<-vector()
p_matrix<-vector()

seriesuniq<-unique(siteyr_data$site)

seriesuniq

for(i in 1:length(seriesuniq)){
  sub<-subset(siteyr_data,site==seriesuniq[i])
  cone<-(sub$cone_total)
  #drop cone production data, all site data and so all htat is left are climate variables
  clim<-sub[,-(1:9)]
  mat<-corr.test(cone,clim,use="pairwise",method="spearman")
  mat_r<-mat$r
  pmat_s<-mat$p
  p_matrix<-rbind(p_matrix,pmat_s)
  d_matrix<-rbind(d_matrix,mat_r)
}

rownames(d_matrix)<-seriesuniq


p_matrix<-as.data.frame(p_matrix)
rownames(p_matrix)<-seriesuniq
names(p_matrix)
p_matrix



#h7<-as.data.frame(pmat[pmat<0.05])

#cone_climate_medians_log<-as.data.frame(apply(d_matrix,2,sd,na.rm=TRUE))
d_matrix
d_matrix<-d_matrix[,-c(1,75)]

data<-d_matrix

set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
data <- d_matrix
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


?kmeans
clust<-kmeans(d_matrix,centers=3,nstart=1)

clust
d_matrix




dsher<-aggregate(d_matrix, by=list(cluster=clust$cluster), mean)
ddd <- cbind(site_cones[,1], cluster = as.factor(clust$cluster))
dd <- cbind(ddd,d_matrix)
dddd<-cbind(site_data, cluster = as.factor(clust$cluster))

dsher

site_data$site
seriesuniq

do<-merge(ddd,treeyr_data,by="site")
di<-merge(ddd,siteyr_data,by="site")

colnames(dsher)[m.col(dsher[,-1], ties.method = "first")]


climcor<-as.data.frame(dd)
climcor$site<-rownames(climcor)
melt_clim<-melt(climcor)

melt_clim<-as.data.table(melt_clim)



new_clim<-melt_clim[,.(q50=median(value),q25=quantile(value)[2],q75=quantile(value)[4]),by=c("variable","cluster")]

new_clim$times<-as.numeric(new_clim$variable)

clim_cor1<-subset(new_clim,cluster==1)
clim_cor2<-subset(new_clim,cluster==2)
clim_cor3<-subset(new_clim,cluster==3)

ppt_box_cor<-melt_clim[grep("ppt",melt_clim$variable)]
vpd_box_cor<-melt_clim[grep("vpd",melt_clim$variable)]
dlt_box_cor<-melt_clim[grep("delta",melt_clim$variable)]

ppt_cor<-new_clim[grep("ppt",new_clim$variable)]
vpd_cor<-new_clim[grep("vpd",new_clim$variable)]
dlt_cor<-new_clim[grep("delta",new_clim$variable)]

ppt_cor.1<-subset(ppt_cor,cluster==1)
ppt_cor.2<-subset(ppt_cor,cluster==2)
ppt_cor.3<-subset(ppt_cor,cluster==3)

vpd_cor.1<-subset(vpd_cor,cluster==1)
vpd_cor.2<-subset(vpd_cor,cluster==2)
vpd_cor.3<-subset(vpd_cor,cluster==3)


ppt_plot<-ggplot()+
  geom_boxplot(ppt_box_cor,mapping=aes(y=value,x=variable,col=cluster))
vpd_plot<-ggplot()+
  geom_boxplot(vpd_box_cor,mapping=aes(y=value,x=variable,col=cluster))
dlt_plot<-ggplot()+
  geom_boxplot(dlt_box_cor,mapping=aes(y=value,x=variable,col=cluster))



vpd_line_plot<-ggplot()+
  geom_line(vpd_cor,mapping=aes(y=q50,x=times,col=cluster))+
  geom_ribbon(vpd_cor,mapping=aes(ymax=q75,ymin=q25,x=times,fill=cluster),alpha=0.15)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=47,lty=1,alpha=0.25)+
  geom_vline(xintercept=59,lty=1,alpha=0.25)+
  geom_vline(xintercept=52,lty=2,alpha=0.6)+
  geom_vline(xintercept=61,lty=2,alpha=0.6)+
  labs(title="VPD",y="Correlation coefficient", x = "Month")+
  facet_wrap(~cluster)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")



vpd_line_plot1<-ggplot()+
  geom_line(vpd_cor.1,mapping=aes(y=q50,x=times),col="red")+
  geom_ribbon(vpd_cor.1,mapping=aes(ymax=q75,ymin=q25,x=times),fill="red",alpha=0.15)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=47,lty=1,alpha=0.25)+
  geom_vline(xintercept=59,lty=1,alpha=0.25)+
  geom_vline(xintercept=52,lty=2,alpha=0.6)+
  geom_vline(xintercept=61,lty=2,alpha=0.6)+
  scale_x_continuous(breaks=NULL)+
  annotate("text", x = 45, y = -.375, label = "yr-3", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 58, y = -.375, label = "yr-2", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 69, y = -.375, label = "yr-1", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 51, y = .5, label = "JJA", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  annotate("text", x = 60, y = .5, label = "MAM", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  labs(title="VPD",y="Correlation", x = "")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),legend.position = "none")


vpd_line_plot2<-ggplot()+
  geom_line(vpd_cor.2,mapping=aes(y=q50,x=times),col="green")+
  geom_ribbon(vpd_cor.2,mapping=aes(ymax=q75,ymin=q25,x=times,fill=cluster),fill="green",alpha=0.15)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=47,lty=1,alpha=0.25)+
  geom_vline(xintercept=59,lty=1,alpha=0.25)+
  geom_vline(xintercept=52,lty=2,alpha=0.6)+
  geom_vline(xintercept=61,lty=2,alpha=0.6)+
  scale_x_continuous(breaks=NULL)+
  annotate("text", x = 45, y = -.2, label = "yr-3", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 58, y = -.2, label = "yr-2", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 69, y = -.2, label = "yr-1", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 51, y = .5, label = "JJA", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  annotate("text", x = 60, y = .5, label = "MAM", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  labs(title="VPD",y="Correlation", x = "Month prior to seed maturation")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),legend.position = "none")



vpd_line_plot3<-ggplot()+
  geom_line(vpd_cor.3,mapping=aes(y=q50,x=times),col="blue")+
  geom_ribbon(vpd_cor.3,mapping=aes(ymax=q75,ymin=q25,x=times,fill=cluster),fill="blue",alpha=0.15)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=47,lty=1,alpha=0.25)+
  geom_vline(xintercept=59,lty=1,alpha=0.25)+
  geom_vline(xintercept=52,lty=2,alpha=0.6)+
  geom_vline(xintercept=61,lty=2,alpha=0.6)+
  scale_x_continuous(breaks=NULL)+
  annotate("text", x = 45, y = -.375, label = "yr-3", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 58, y = -.375, label = "yr-2", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 69, y = -.375, label = "yr-1", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 51, y = .5, label = "JJA", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  annotate("text", x = 60, y = .5, label = "MAM", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  labs(title="VPD",y="Correlation", x = "")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),legend.position = "none")


vpd_line_plot1/ppt_line_plot1|vpd_line_plot2/ppt_line_plot2|vpd_line_plot3/ppt_line_plot3


vpd_line_plot1|vpd_line_plot2|vpd_line_plot3

ppt_var<-c("Yr-2", "JJA","Yr-1","MAM")
ppt_var


ppt_line_plot<-ggplot()+
  geom_line(ppt_cor,mapping=aes(y=q50,x=times,col=cluster))+
  geom_ribbon(ppt_cor,mapping=aes(ymax=q75,ymin=q25,x=times,fill=cluster),alpha=0.15)+
  geom_hline(yintercept=0,alpha=0.75)+
  facet_wrap(~cluster)+
  geom_vline(xintercept=13,lty=1,alpha=0.25)+
  geom_vline(xintercept=25,lty=1,alpha=0.25)+
  geom_vline(xintercept=18,lty=2,alpha=0.6)+
  geom_vline(xintercept=27,lty=2,alpha=0.6)+
  scale_x_continuous(breaks=NULL)+
  labs(title="PPT",y="Correlation", x = "")+
  annotate("text", x = 11, y = .5, label = "yr-3", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 24, y = .5, label = "yr-2", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 35, y = .5, label = "yr-1", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 17, y = -.25, label = "JJA", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  annotate("text", x = 26, y = -.25, label = "MAM", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  theme_minimal()+
  facet_wrap(~cluster)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")

ppt_line_plot1<-ggplot()+
  geom_line(ppt_cor.1,mapping=aes(y=q50,x=times),col="red")+
  geom_ribbon(ppt_cor.1,mapping=aes(ymax=q75,ymin=q25,x=times),fill="red",alpha=0.15)+
  geom_hline(yintercept=0,alpha=0.75)+
  geom_vline(xintercept=13,lty=1,alpha=0.25)+
  geom_vline(xintercept=25,lty=1,alpha=0.25)+
  geom_vline(xintercept=18,lty=2,alpha=0.6)+
  geom_vline(xintercept=27,lty=2,alpha=0.6)+
  scale_x_continuous(breaks=NULL)+
  labs(title="PPT",y="Correlation", x = "")+
  annotate("text", x = 11, y = .5, label = "yr-3", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 24, y = .5, label = "yr-2", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 35, y = .5, label = "yr-1", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 17, y = -.25, label = "JJA", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  annotate("text", x = 26, y = -.25, label = "MAM", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")



ppt_line_plot2<-ggplot()+
  geom_line(ppt_cor.2,mapping=aes(y=q50,x=times),col="green")+
  geom_ribbon(ppt_cor.2,mapping=aes(ymax=q75,ymin=q25,x=times),fill="green",alpha=0.15)+
  geom_hline(yintercept=0,alpha=0.75)+
  geom_vline(xintercept=13,lty=1,alpha=0.25)+
  geom_vline(xintercept=25,lty=1,alpha=0.25)+
  geom_vline(xintercept=18,lty=2,alpha=0.6)+
  geom_vline(xintercept=27,lty=2,alpha=0.6)+
  scale_x_continuous(breaks=NULL)+
  labs(title="PPT",y="Correlation", x = "")+
  annotate("text", x = 11, y = .5, label = "yr-3", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 24, y = .5, label = "yr-2", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 35, y = .5, label = "yr-1", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 17, y = -.25, label = "JJA", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  annotate("text", x = 26, y = -.25, label = "MAM", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")




ppt_line_plot3<-ggplot()+
  geom_line(ppt_cor.3,mapping=aes(y=q50,x=times),col="blue")+
  geom_ribbon(ppt_cor.3,mapping=aes(ymax=q75,ymin=q25,x=times),fill='blue',alpha=0.15)+
  geom_hline(yintercept=0,alpha=0.75)+
  geom_vline(xintercept=13,lty=1,alpha=0.25)+
  geom_vline(xintercept=25,lty=1,alpha=0.25)+
  geom_vline(xintercept=18,lty=2,alpha=0.6)+
  geom_vline(xintercept=27,lty=2,alpha=0.6)+
  scale_x_continuous(breaks=NULL)+
  labs(title="PPT",y="Correlation", x = "")+
  annotate("text", x = 11, y = .5, label = "yr-3", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 24, y = .5, label = "yr-2", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 35, y = .5, label = "yr-1", angle = 0,vjust = 1, hjust=1,alpha=0.3)+
  annotate("text", x = 17, y = -.25, label = "JJA", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  annotate("text", x = 26, y = -.25, label = "MAM", angle = 45,vjust = 1, hjust=1,alpha=0.6)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")

vpd_line_plot1/ppt_line_plot1|vpd_line_plot2/ppt_line_plot2|vpd_line_plot3/ppt_line_plot3


vpd_line_plot/ppt_line_plot

ggplot(dddd)+
  geom_density(mapping=aes(group=cluster,x=AET,col=cluster))

ggplot(dddd)+
  geom_boxplot(mapping=aes(group=cluster,y=AET,col=cluster))

ggplot(dddd)+
  geom_point(mapping=aes(x=lon,y=lat,col=cluster))

vpd_plot/ppt_plot

?IQR

chuck<-pairwise.t.test(dddd$CWD,dddd$cluster)
chuck
summary(chuck)
summary(model)

str(dddd)


ppt_plot/vpd_plot
dlt_plot



new_clim$times<-as.numeric(new_clim$variable)


clime<-siteyr_data[,-(1:9)]
suh<-clime[,c(77,78,61,40,53,52,18,19)]

plot(suh)

names(p_matrix)
cormat<-cbind.data.frame(seriesuniq,d_matrix)
colnames(cormat[1])<-"tree"
#ps<-as.data.frame(p_matrix[,c(19,60,104,106)])
#ts<-as.data.frame(d_matrix[,c(19,60,104,106)])




#creating a rough binary category of cones no cones
siteyr_data$mast<-ifelse(siteyr_data$cone_total>0,1,0)
treeyr_data$mast<-ifelse(treeyr_data$cone_total>0,1,0)

#for running models
options(na.action = "na.omit")

#distribution assumptions


mod_l2<-glmmTMB(log_cones~1+(1|site),data=siteyr_data)
mod_p1<-glmmTMB(cone_total~1+(1|site),data=siteyr_data,family=poisson(link="log"))
mod_pz<-glmmTMB(cone_total~1+(1|site),data=siteyr_data,ziformula = ~.,family=poisson(link="log"))
mod_n1<-glmmTMB(cone_total~1+(1|site),data=siteyr_data,family=nbinom2(link="log"))
mod_nz<-glmmTMB(cone_total~1+(1|site),data=siteyr_data,ziformula = ~.,family=nbinom2(link="log"))


AIC(mod_l2)
AIC(mod_p1)
AIC(mod_pz)
AIC(mod_n1)
AIC(mod_nz)



#here's ya null
null_mod<-glmmTMB(cone_total~1+(1|site/tree),data=treeyr_data,ziformula = ~.,family=nbinom2(link="log"))
AIC(null_mod)


mod0<-glmmTMB(cone_total~scale(AR1)+scale(deltaT)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mod1<-glmmTMB(cone_total~scale(AR1)+scale(deltaV)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mod2<-glmmTMB(cone_total~scale(AR1)+scale(yr1pptjja)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mod3<-glmmTMB(cone_total~scale(AR1)+scale(yr1vpdjas)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mod4<-glmmTMB(cone_total~scale(AR1)+scale(yr2vpdjja)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))


mo00<-glmmTMB(cone_total~scale(AR1)+yr1pptjja+yr2vpdjja+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo11<-glmmTMB(cone_total~scale(AR1)+yr1vpdjas+yr2vpdjja+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo12<-glmmTMB(cone_total~scale(AR1)+yr1vpdjas+yr1pptjja+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))

mo01<-glmmTMB(cone_total~scale(AR1)+scale(yr1pptjja)*scale(CWD)+scale(yr2vpdjja)*scale(AET)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo02<-glmmTMB(cone_total~scale(AR1)+scale(yr1pptjja)*scale(CWD)+scale(yr2vpdjja)*scale(CWD)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo03<-glmmTMB(cone_total~scale(AR1)+scale(yr1pptjja)*scale(AET)+scale(yr2vpdjja)*scale(AET)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo04<-glmmTMB(cone_total~scale(AR1)+scale(yr1pptjja)*scale(AET)+scale(yr2vpdjja)*scale(CWD)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))


AIC(mod0)
AIC(mod1)-AIC(mo01)
AIC(mod2)
AIC(mod3)-AIC(mo01)
AIC(mod4)

AIC(mo00)
AIC(mo11)-AIC(mo01)
AIC(mo12)-AIC(mo01)

AIC(mo01)-AIC(mo01)
AIC(mo02)-AIC(mo01)
AIC(mo03)-AIC(mo01)
AIC(mo04)-AIC(mo01)

AIC(mo01)-AIC(mo00)
AIC(mo02)-AIC(mo00)
AIC(mo03)-AIC(mo00)
AIC(mo04)-AIC(mo00)

#just in case

mo15<-glmmTMB(cone_total~scale(AR1)+deltaV* (AET) + (1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo16<-glmmTMB(cone_total~scale(AR1)+deltaV* (CWD) + (1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo17<-glmmTMB(cone_total~scale(AR1)+deltaP* (AET) + (1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo18<-glmmTMB(cone_total~scale(AR1)+deltaP* (CWD) + (1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo19<-glmmTMB(cone_total~scale(AR1)+deltaT* (AET) + (1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo20<-glmmTMB(cone_total~scale(AR1)+deltaT* (CWD) + (1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo21<-glmmTMB(cone_total~scale(AR1)+yr1pptjja*(AET)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo22<-glmmTMB(cone_total~scale(AR1)+yr1pptjja*(CWD)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo23<-glmmTMB(cone_total~scale(AR1)+yr1vpdjja*(AET)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo24<-glmmTMB(cone_total~scale(AR1)+yr1vpdjja*(CWD)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo25<-glmmTMB(cone_total~scale(AR1)+yr1vpdmam*(AET)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo26<-glmmTMB(cone_total~scale(AR1)+yr1vpdmam*(CWD)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo27<-glmmTMB(cone_total~scale(AR1)+yr1pptmam*(AET)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
mo28<-glmmTMB(cone_total~scale(AR1)+yr1pptmam*(CWD)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))


summary(mo15)
summary(mo16)
summary(mo17)
summary(mo18)
summary(mo19)
summary(mo20)
summary(mo21)
summary(mo22)
summary(mo23)
summary(mo24)
summary(mo25)
summary(mo26)
summary(mo27)
summary(mo28)


smo15<-glmmTMB(cone_total~scale(AR1)+deltaV* (AET) + (1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo16<-glmmTMB(cone_total~scale(AR1)+deltaV* (CWD) + (1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo17<-glmmTMB(cone_total~scale(AR1)+deltaP* (AET) + (1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo18<-glmmTMB(cone_total~scale(AR1)+deltaP* (CWD) + (1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo19<-glmmTMB(cone_total~scale(AR1)+deltaT* (AET) + (1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo20<-glmmTMB(cone_total~scale(AR1)+deltaT* (CWD) + (1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo21<-glmmTMB(cone_total~scale(AR1)+yr1pptjja*(AET)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo22<-glmmTMB(cone_total~scale(AR1)+yr1pptjja*(CWD)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo23<-glmmTMB(cone_total~scale(AR1)+yr1vpdjja*(AET)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo24<-glmmTMB(cone_total~scale(AR1)+yr1vpdjja*(CWD)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo25<-glmmTMB(cone_total~scale(AR1)+yr1vpdmam*(AET)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo26<-glmmTMB(cone_total~scale(AR1)+yr1vpdmam*(CWD)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo27<-glmmTMB(cone_total~scale(AR1)+yr1pptmam*(AET)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo28<-glmmTMB(cone_total~scale(AR1)+yr1pptmam*(CWD)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))


smo29<-glmmTMB(cone_total~scale(AR1)+yr1pptjja*CWD+yr1pptmam*(CWD)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
smo29<-glmmTMB(cone_total~scale(AR1)+yr1pptjja*CWD+yr1pptmam*(AET)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))
smo29<-glmmTMB(cone_total~scale(AR1)+yr1pptjja*CWD+yr1vpdmam*(AET)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))


smo30<-glmmTMB(cone_total~scale(AR1)+deltaV*CWD+yr1pptmam*(CWD)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo31<-glmmTMB(cone_total~scale(AR1)+yr1pptmam*(CWD)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))
smo32<-glmmTMB(cone_total~scale(AR1)+yr1pptmam*(CWD)+(1|site),ziformula= ~.,data=siteyr_data,family=nbinom2(link="log"))



AIC(smo15)
AIC(smo16)
AIC(smo17)
AIC(smo18)
AIC(smo19)
AIC(smo20)
AIC(smo21)
AIC(smo22)
AIC(smo23)
AIC(smo24)
AIC(smo25)
AIC(smo26)
AIC(smo27)
AIC(smo28)


summary(smo15)
summary(smo16)
summary(smo17)
summary(smo18)
summary(smo19)
summary(smo20)
summary(smo21)
summary(smo22)
summary(smo23)
summary(smo24)
summary(smo25)
summary(smo26)
summary(smo27)
summary(smo28)

summary(smo29)
summary(smo30)

AIC(mo15)
AIC(mo16)
AIC(mo17)
AIC(mo18)



set.seed(7)
performance::r2(mo00)
summary(mo00)


plot(predict(mo01,type="response"),treeyr_data$cone_total)
abline(a=0, b=1)


fit3<-glmmTMB(cone_total~scale(AR1)+yr2vpdjja*scale(AET)+yr1pptjja*scale(CWD)+(1|site/tree),ziformula= ~.,data=treeyr_data,family=nbinom2(link="log"))

x1a_range<-seq(from=min(siteyr_data$yr2vpdjja),to=max(siteyr_data$yr2vpdjja),by=.1)
x1b_range<-seq(from=min(siteyr_data$yr1pptjja),to=max(siteyr_data$yr1pptjja),by=.1)

quantile(siteyr_data$AET,probs=c(.10,.90))
quantile(siteyr_data$CWD,probs=c(.10,.90))

x2a_l <- 300
x2a_h <- 500

x2b_l <- 390
x2b_h <- 740

x3 <- -1

x4a= 0

x5a<-mean(siteyr_data$AET)
x5b<-mean(siteyr_data$CWD)

#visreg(mo174,"yr1pptjja",by="CWD",scale="response")

generated_data_a <- as.data.frame(expand.grid(yr2vpdjja=x1a_range, yr1pptjja=x4a, CWD= x5b, AET=c(x2a_l, x2a_h), AR1=x3, site="a", tree="AA")) 
generated_data_b <- as.data.frame(expand.grid(yr1pptjja=x1b_range, yr2vpdjja=x4a, AET= x5a, CWD=c(x2b_l, x2b_h), AR1=x3, site="a", tree="AA")) 



generated_data_a$z <- predict(fit3, newdata=generated_data_a, type = 'zprob',re.form=NA)
generated_data_a$c <- predict(fit3, newdata=generated_data_a, type = 'response',re.form=NA)
generated_data_a$aetlevel <- factor(generated_data_a$AET, labels=c("300 mm", "500 mm"), ordered=T)

generated_data_b$z <- predict(fit3, newdata=generated_data_b, type = 'zprob',re.form=NA)
generated_data_b$c <- predict(fit3, newdata=generated_data_b, type = 'response',re.form=NA)
generated_data_b$cwdlevel <- factor(generated_data_b$CWD, labels=c("390 mm", "740 mm"), ordered=T)

log.mu=predict(fit3, newdata=generated_data_a, type = 'link',re.form=NA,se=TRUE)
logit.z=predict(fit3, newdata=generated_data_a, type = 'zlink',re.form=NA,se=TRUE)


log.mub=predict(fit3, newdata=generated_data_b, type = 'link',re.form=NA,se=TRUE)
logit.zb=predict(fit3, newdata=generated_data_b, type = 'zlink',re.form=NA,se=TRUE)


pred=exp(log.mu$fit)*(1-inv.logit(logit.z$fit))
ci.lo=exp(log.mu$fit-1.96*log.mu$se.fit)*(1-boot::inv.logit(logit.z$fit+1.96*logit.z$se.fit))
ci.hi=exp(log.mu$fit+1.96*log.mu$se.fit)*(1-boot::inv.logit(logit.z$fit-1.96*logit.z$se.fit))

predb=exp(log.mub$fit)*(1-inv.logit(logit.zb$fit))
ci.lob=exp(log.mub$fit-1.96*log.mub$se.fit)*(1-boot::inv.logit(logit.zb$fit+1.96*logit.zb$se.fit))
ci.hib=exp(log.mub$fit+1.96*log.mub$se.fit)*(1-boot::inv.logit(logit.zb$fit-1.96*logit.zb$se.fit))


yuh<-cbind.data.frame(pred,ci.lo,ci.hi)
yuhb<-cbind.data.frame(predb,ci.lob,ci.hib)

###NOTE HERE sign is flipped for ZINF predictions.  So officially predictiong a mast year, not failure, 6/13/22

zprob=(1-boot::inv.logit(logit.z$fit))


ci.hi=(1-boot::inv.logit(logit.z$fit-1.96*logit.z$se.fit))
ci.lo=(1-boot::inv.logit(logit.z$fit+1.96*logit.z$se.fit))


zprobb=(1-boot::inv.logit(logit.zb$fit))

ci.hib=(1-boot::inv.logit(logit.zb$fit-1.96*logit.zb$se.fit))
ci.lob=(1-boot::inv.logit(logit.zb$fit+1.96*logit.zb$se.fit))

zpro<-cbind.data.frame(zprob,ci.lo,ci.hi)
zpro_b<-cbind.data.frame(zprobb,ci.lob,ci.hib)

tost<-cbind(yuh,generated_data_a)
testt<-cbind(zpro,generated_data_a)

tostb<-cbind(yuhb,generated_data_b)
testtb<-cbind(zpro_b,generated_data_b)


#treeyr_data$cone_prs<-ifelse(treeyr_data$cone_total>0,treeyr_data$cone_total,NA)


int1<-ggplot()+
  geom_line(testt,mapping=aes(x=yr2vpdjja,y=zprob,col=aetlevel),size=1.5)+
  geom_ribbon(testt,mapping=aes(x=yr2vpdjja,ymax=ci.hi,ymin=ci.lo,fill=aetlevel),alpha=0.2,show.legend=F)+
  theme_classic()+
  scale_color_scico_d(palette = "vik",begin=0.3, direction=-1, end=0.725,"30-yr AET")+
  scale_fill_scico_d(palette = "vik",begin=0.3, direction=-1, end=0.725,"30-yr AET")+
  annotate("text",label=expression(italic("P")~" = 0.03"),x=1.5,y=0.85,parse=TRUE,size=3)+
  labs(x="VPD, June-August (T-3)",y=expression(theta~"Cone Presence (T)"))
  #theme(legend.position="none")


int2<-ggplot()+
  #geom_point(treeyr_data,mapping=aes(x=deltaV,y=cone_prs),alpha=0.05)+
  geom_line(tost,mapping=aes(x=yr2vpdjja,y=pred,col=aetlevel),size=1.5)+
  geom_ribbon(tost,mapping=aes(x=yr2vpdjja,ymax=ci.hi,ymin=ci.lo,fill=aetlevel),alpha=0.2,show.legend=F)+
  #  coord_cartesian(ylim=c(0,155))+
  scale_color_scico_d(palette = "vik",begin=0.3, direction=-1,end=0.725,"30-yr AET")+
  scale_fill_scico_d(palette = "vik",begin=0.3, direction=-1,end=0.725,"30-yr AET")+
  labs(x="VPD, June-August (T-3)",y='Cone Abundance (T)')+
  annotate("text",label=expression(italic("P")~" = 0.30"),x=1.5,y=175,parse=TRUE,size=3)+
  theme_classic()


int3<-ggplot()+
  geom_line(testtb,mapping=aes(x=yr1pptjja,y=zprobb,col=cwdlevel),size=1.5)+
  geom_ribbon(testtb,mapping=aes(x=yr1pptjja,ymax=ci.hib,ymin=ci.lob,fill=cwdlevel),alpha=0.2,show.legend=F)+
  theme_classic()+
  scale_color_scico_d(palette = "vikO",begin=0.3, end=0.725,"30-yr CWD")+
  scale_fill_scico_d(palette = "vikO",begin=0.3, end=0.725,"30-yr CWD")+
  annotate("text",label=expression(italic("P")~" < 0.001"),x=2,y=0.875,parse=TRUE,size=3)+
  labs(x="PPT, June-August (T-2)",y=expression(theta~"Cone Presence (T)"))
  #theme(legend.position = "none")

int4<-ggplot()+
  #geom_point(treeyr_data,mapping=aes(x=deltaV,y=cone_prs),alpha=0.05)+
  geom_line(tostb,mapping=aes(x=yr1pptjja,y=predb,col=cwdlevel),size=1.5)+
  geom_ribbon(tostb,mapping=aes(x=yr1pptjja,ymax=ci.hib,ymin=ci.lob,fill=cwdlevel),alpha=0.2,show.legend=F)+
  #  coord_cartesian(ylim=c(0,155))+
  scale_color_scico_d(palette = "vikO",begin=0.3, end=0.725,"30-yr CWD")+
  scale_fill_scico_d(palette = "vikO",begin=0.3, end=0.725,"30-yr CWD")+
  annotate("text",label=expression(italic("P")~" = 0.42"),x=2,y=140,parse=TRUE,size=3)+
  labs(x="PPT, June-August (T-2)",y='Cone Abundance (T)')+
  theme_classic()

int4

#png(file="C:/Users/Andreas/Google Drive/Cones_Ponderosa/figs/int.png",width=20, height=17.5, bg="white", units="cm", res= 300)
(int1|int3)/(int2|int4)
dev.off()





scico(100,palette = 'vikO')
#

























####AGE MODELS
options(na.action = "na.omit")

aged<-merge(age,tree_data,by="tree")

aged<-aged[!is.na(aged$tree_r),]
aged<-aged[!is.na(aged$CVi),]
aged<-aged[!is.na(aged$BA5),]
aged<-aged[!is.na(aged$prop),]

#oaged<-aged
#
aged$log_mean_cones<-log(aged$mean_cones)
#
#oaged$log_mean_cones<-log(oaged$mean_cones)
#oaged$log_CVi<-log(oaged$CVi)
#
#aged<-subset(aged,tree!="MON2T5")
#
#
#oa_mod_mean<-glmmTMB(log_mean_cones~scale(dbh)+scale(BA5)+scale(AET)+scale(CWD)+scale(age)+(1|site),data=oaged)
#oa_mod_r <- glmmTMB(tree_r~scale(dbh)+scale(BA5)+scale(AET)+scale(CWD)+scale(age)+(1|site),data=oaged)
#oa_mod_cv<-glmmTMB(CVi~scale(dbh)+scale(BA5)+scale(AET)+scale(CWD)+scale(age)+(1|site),data=oaged)
#oa_mod_ar<-glmmTMB(prop~scale(dbh)+scale(BA5)+scale(AET)+scale(CWD)+scale(age)+(1|site),data=oaged,family="binomial")


a_mod_mean<-glmmTMB(log_mean_cones~scale(poly(dbh,2))+scale(poly(BA5,2))+scale(AET)+scale(CWD)+scale(age)+(1|site),data=aged)
a_mod_r <-  glmmTMB(tree_r~(poly(dbh,2))+(poly(BA5,2))+(AET)+(CWD)+(age)+(1|site),data=aged)
a_mod_cv<-  glmmTMB(CVi~scale(poly(dbh,2))+(poly(BA5,2))+scale(AET)+scale(CWD)+scale(age)+(1|site),data=aged)


a_mod_mean2<-glmmTMB(log_mean_cones~scale(poly(dbh,2))+scale(poly(BA5,2))+scale(AET)+scale(CWD)+scale(poly(age,2))+(1|site),data=aged)
a_mod_r2 <-  glmmTMB(tree_r~scale(poly(dbh,2))+scale(poly(BA5,2))+scale(AET)+scale(CWD)+scale(poly(age,2))+(1|site),data=aged)
a_mod_cv2<-  glmmTMB(CVi~scale(poly(dbh,2))+(poly(BA5,2))+scale(AET)+scale(CWD)+scale(poly(age,2))+(1|site),data=aged)


a_mod_mean3<-glmmTMB(log_mean_cones~(poly(dbh,2))+(poly(BA5,2))+(AET)+(CWD)+(poly(age,2))+(1|site),data=aged)
a_mod_r3 <-  glmmTMB(tree_r~(poly(dbh,2))+(poly(BA5,2))+(AET)+(CWD)+(poly(age,2))+(1|site),data=aged)
a_mod_cv3<-  glmmTMB(CVi~(poly(dbh,2))+(poly(BA5,2))+(AET)+(CWD)+(poly(age,2))+(1|site),data=aged)

check_collinearity(a_mod_cv2)


AIC(a_mod_mean)
AIC(a_mod_mean2)
AIC(a_mod_r )
AIC(a_mod_r2)
AIC(a_mod_cv)
AIC(a_mod_cv2)


options(na.action = "na.fail")

amd<-dredge(a_mod_mean2,rank="AIC")
amr<-dredge(a_mod_r2,rank="AIC")
amc<-dredge(a_mod_cv2,rank="AIC")

testUniformity(simulateResiduals(a_mod_mean2))

importance(amd)
importance(amr)
importance(amc)

amd
amr
amc

r.squaredGLMM(a_mod_mean3)
r.squaredGLMM(a_mod_r3)
r.squaredGLMM(a_mod_cv3)


samd<-summary(model.avg(amd))
samr<-summary(model.avg(amr))
samc<-summary(model.avg(amc))


samd
samr
samc


plot(aged$age,aged$log_mean_cones)

scico(palette = "vikO",100)

col_blu<-"#6E98B8"
col_red<-"#B25F3B"











mean_ba_preds<-ggpredict(a_mod_mean3,"BA5[all]",type="random")

mean_ba<-ggplot() +
  geom_ribbon(mean_ba_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill=col_red,alpha=0.33)+
  geom_line(mean_ba_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_cv_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=BA5,y=log_mean_cones),alpha=0.2)+
  labs(x="",y="")+
  annotate("text",label=expression(italic("P")~"= 0.03"),x=15,y=5.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


mean_ba

mean_dbh_preds<-ggpredict(a_mod_mean3,"dbh[all]",type="random")

mean_dbh<-ggplot() +
  geom_ribbon(mean_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill=col_blu,alpha=0.33)+
  geom_line(mean_dbh_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_cv_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=dbh,y=log_mean_cones),alpha=0.2)+
  labs(x="",y="")+
  annotate("text",label=expression(italic("P")~"< 0.001"),x=50,y=6,parse=TRUE, hjust=0,size=3)+
  theme_classic()


mean_dbh


mean_age_preds<-ggpredict(a_mod_mean3,"age[all]",type="random")

mean_age<-ggplot() +
  geom_ribbon(mean_age_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill=col_blu,alpha=0.33)+
  geom_line(mean_age_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_cv_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=age,y=log_mean_cones),alpha=0.2)+
  labs(x="",y="log (Mean cones)")+
  annotate("text",label=expression(italic("P")~"= 0.004"),x=175,y=5.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


mean_age



cv_ba_preds<-ggpredict(a_mod_cv3,"BA5[all]",type="random")

cv_ba<-ggplot() +
  geom_ribbon(cv_ba_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill=col_blu,alpha=0.33)+
  geom_line(cv_ba_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_cv_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_cv_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=BA5,y=CVi),alpha=0.2)+
  labs(x="",y="")+
  annotate("text",label=expression(italic("P")~"= 0.02"),x=15,y=3.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


cv_ba

cv_dbh_preds<-ggpredict(a_mod_cv3,"dbh[all]",type="random")

cv_dbh<-ggplot() +
  geom_ribbon(cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill=col_red,alpha=0.33)+
  geom_line(cv_dbh_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_cv_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=dbh,y=CVi),alpha=0.2)+
  labs(x="",y="")+
  annotate("text",label=expression(italic("P")~"< 0.001"),x=50,y=3.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


cv_dbh


cv_age_preds<-ggpredict(a_mod_cv3,"age[all]",type="random")

cv_age<-ggplot() +
  geom_ribbon(cv_age_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill=col_red,alpha=0.33)+
  geom_line(cv_age_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_cv_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=age,y=CVi),alpha=0.2)+
  labs(x="",y="CVi")+
  annotate("text",label=expression(italic("P")~"= 0.03"),x=175,y=3.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


cv_age






r_ba_preds<-ggpredict(a_mod_r3,"BA5[all]",type="random")

r_ba<-ggplot() +
  geom_ribbon(r_ba_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill=col_red,alpha=0.33)+
  geom_line(r_ba_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_r_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=BA5,y=tree_r),alpha=0.2)+
  labs(x=expression("Basal Area (m"^{"2"}~~"/ ha)"),y="")+
  annotate("text",label=expression(italic("P")~"= 0.01"),x=15,y=1.0,parse=TRUE, hjust=0,size=3)+
  theme_classic()


r_ba

r_dbh_preds<-ggpredict(a_mod_r3,"dbh[all]",type="random")

r_dbh<-ggplot() +
  geom_ribbon(r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill=col_red,alpha=0.33)+
  geom_line(r_dbh_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=dbh,y=tree_r),alpha=0.2)+
  labs(x="DBH (cm)",y="")+
  annotate("text",label=expression(italic("P")~"= 0.04"),x=50,y=1.0,parse=TRUE, hjust=0,size=3)+
  theme_classic()


r_dbh

#r_age_preds<-ggpredict(a_mod_r3,"age",type="random")
r_age_predsl<-ggpredict(a_mod_r,"age[all]",type="random")

r_age<-ggplot() +
#  geom_ribbon(r_age_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90",alpha=0.5)+
#  geom_line(r_age_preds, mapping=aes(x, predicted),lty=2) +
  geom_ribbon(r_age_predsl,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill=col_red,alpha=0.33)+
  geom_line(r_age_predsl, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=age,y=tree_r),alpha=0.2)+
  labs(x="Age",y="Synchrony")+
  annotate("text",label=expression(italic("P")~"= 0.01"),x=175,y=1,parse=TRUE, hjust=0,size=3)+
#  annotate("text",label=expression(Age^{"2"}~~italic("P")~"= 0.25"),x=140,y=0.9,parse=TRUE, hjust=0,size=3)+
  theme_classic()


r_age


(mean_age|mean_dbh|mean_ba)/(cv_age|cv_dbh|cv_ba)/(r_age|r_dbh|r_ba)



#png(file="C:/Users/Andreas/Google Drive/Cones_Ponderosa/figs/trees.png",width=20, height=20, bg="white", units="cm", res= 300)

(mean_age|mean_dbh|mean_ba)/(cv_age|cv_dbh|cv_ba)/(r_age|r_dbh|r_ba)
dev.off()







r_cwd_preds<-ggpredict(a_mod_r3,"CWD[all]",type="random")

r_cwd<-ggplot() +
  geom_ribbon(r_cwd_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(r_cwd_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=CWD,y=tree_r),alpha=0.2)+
  labs(x="CWD",y="Synchrony")+
  annotate("text",label=expression(italic("P")~"= 0.22"),x=650,y=0.65,parse=TRUE, hjust=0,size=3)+
  theme_classic()


r_cwd





r_aet_preds<-ggpredict(a_mod_r3,"AET[all]",type="random")

r_aet<-ggplot() +
  geom_ribbon(r_aet_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(r_aet_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=AET,y=tree_r),alpha=0.2)+
  labs(x="AET",y="Synchrony")+
  annotate("text",label=expression(italic("P")~"= 0.35"),x=500,y=0.65,parse=TRUE, hjust=0,size=3)+
  theme_classic()


r_aet


cv_cwd_preds<-ggpredict(a_mod_cv3,"CWD[all]",type="random")

cv_cwd<-ggplot() +
  geom_ribbon(cv_cwd_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(cv_cwd_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=CWD,y=CVi),alpha=0.2)+
  labs(x="CWD",y="CVi")+
  annotate("text",label=expression(italic("P")~"= 0.22"),x=650,y=3.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


cv_cwd

cv_aet_preds<-ggpredict(a_mod_cv3,"AET[all]",type="random")

cv_aet<-ggplot() +
  geom_ribbon(cv_aet_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(cv_aet_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=AET,y=CVi),alpha=0.2)+
  labs(x="AET",y="CVi")+
  annotate("text",label=expression(italic("P")~"= 0.52"),x=500,y=3.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


cv_aet




mean_cwd_preds<-ggpredict(a_mod_cv3,"CWD[all]",type="random")

mean_cwd<-ggplot() +
  geom_ribbon(mean_cwd_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(mean_cwd_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=CWD,y=log_mean_cones),alpha=0.2)+
  labs(x="CWD",y="log (Mean cones)")+
  annotate("text",label=expression(italic("P")~"= 0.34"),x=650,y=5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


mean_cwd

mean_aet_preds<-ggpredict(a_mod_cv3,"AET[all]",type="random")

mean_aet<-ggplot() +
  geom_ribbon(mean_aet_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(mean_aet_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(aged,mapping=aes(x=AET,y=log_mean_cones),alpha=0.2)+
  labs(x="AET",y="log (Mean cones)")+
  annotate("text",label=expression(italic("P")~"= 0.61"),x=500,y=5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


mean_aet





#png(file="C:/Users/Andreas/Google Drive/Cones_Ponderosa/figs/climate.png",width=20, height=12.5, bg="white", units="cm", res= 300)

(mean_aet/mean_cwd)|(cv_aet/cv_cwd)|(r_aet/r_cwd)
dev.off()







tree_data_<-tree_data[!is.na(tree_data$tree_r),]
tree_data_<-tree_data_[!is.na(tree_data_$CVi),]
tree_data_<-tree_data_[!is.na(tree_data_$BA5),]
tree_data_<-tree_data_[!is.na(tree_data_$prop),]

tree_data_$log_mean_cones<-log(tree_data_$mean_cones)



t_mod_mean<-glmmTMB(log_mean_cones~scale(poly(dbh,2))+scale(poly(BA5,2))+scale(AET)+scale(CWD)+(1|site),data=tree_data_)
t_mod_cv<-glmmTMB(CVi~scale(poly(dbh,2))+scale(poly(BA5,2))+scale(AET)+scale(CWD)+(1|site),data=tree_data_)
t_mod_r<-glmmTMB(tree_r~scale(poly(dbh,2))+scale(poly(BA5,2))+scale(AET)+scale(CWD)+(1|site),data=tree_data_)




t_mod_mean3<-glmmTMB(log_mean_cones~(poly(dbh,2))+(poly(BA5,2))+(AET)+(CWD)+(1|site),data=tree_data_)
t_mod_cv3<-glmmTMB(CVi~(poly(dbh,2))+(poly(BA5,2))+(AET)+(CWD)+(1|site),data=tree_data_)
t_mod_r3<-glmmTMB(tree_r~(poly(dbh,2))+(poly(BA5,2))+(AET)+(CWD)+(1|site),data=tree_data_)



options(na.action = "na.fail")

tmd<-dredge(t_mod_mean,rank="AIC")
tmr<-dredge(t_mod_r,rank="AIC")
tmc<-dredge(t_mod_cv,rank="AIC")

importance(tmd)
importance(amd)

importance(tmr)
importance(amr)

importance(tmc)
importance(amc)



tmd
tmr
tmc

r.squaredGLMM(t_mod_mean)
r.squaredGLMM(t_mod_r)
r.squaredGLMM(t_mod_cv)

summary(model.avg(tmd))
summary(model.avg(tmr))
summary(model.avg(tmc))



stmd<-summary(model.avg(tmd))
#tmd_tab<-as.data.frame(stmd$coefmat.subset)
stmr<-summary(model.avg(tmr))
#tmr_tab<-as.data.frame(stmr$coefmat.subset)
stmc<-summary(model.avg(tmc))
#tmc_tab<-as.data.frame(stmc$coefmat.subset)

stmd




tmean_ba_preds<-ggpredict(t_mod_mean3,"BA5[all]",type="random")

tmean_ba<-ggplot() +
  geom_ribbon(tmean_ba_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(tmean_ba_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_cv_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(tree_data_,mapping=aes(x=BA5,y=log_mean_cones),alpha=0.2)+
  labs(x="Basal Area",y="Mean Cone Production")+
  annotate("text",label=expression(italic("P")~"= 0.08"),x=15,y=5.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


tmean_ba

tmean_dbh_preds<-ggpredict(t_mod_mean3,"dbh[all]",type="random")

tmean_dbh<-ggplot() +
  geom_ribbon(tmean_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(tmean_dbh_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_cv_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(tree_data_,mapping=aes(x=dbh,y=log_mean_cones),alpha=0.2)+
  labs(x="DBH",y="Mean Cone Production")+
  annotate("text",label=expression(italic("P")~"= 0.08"),x=70,y=5.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


tmean_dbh



tcv_ba_preds<-ggpredict(t_mod_cv3,"BA5[all]",type="random")

tcv_ba<-ggplot() +
  geom_ribbon(tcv_ba_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(tcv_ba_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_cv_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_cv_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(tree_data_,mapping=aes(x=BA5,y=CVi),alpha=0.2)+
  labs(x="Basal Area",y="CVi")+
  annotate("text",label=expression(italic("P")~"= 0.08"),x=15,y=5.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


tcv_ba

tcv_dbh_preds<-ggpredict(t_mod_cv3,"dbh[all]",type="random")

tcv_dbh<-ggplot() +
  geom_ribbon(tcv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(tcv_dbh_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_cv_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_cv_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(tree_data_,mapping=aes(x=dbh,y=CVi),alpha=0.2)+
  labs(x="DBH",y="CVi")+
  annotate("text",label=expression(italic("P")~"= 0.08"),x=70,y=5.5,parse=TRUE, hjust=0,size=3)+
  theme_classic()


tcv_dbh




tr_ba_preds<-ggpredict(t_mod_r3,"BA5[all]",type="random")

tr_ba<-ggplot() +
  geom_ribbon(tr_ba_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(tr_ba_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_r_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(tree_data_,mapping=aes(x=BA5,y=tree_r),alpha=0.2)+
  labs(x="Basal Area",y="Synchrony")+
  annotate("text",label=expression(italic("P")~"= 0.08"),x=15,y=1.0,parse=TRUE, hjust=0,size=3)+
  theme_classic()


tr_ba

tr_dbh_preds<-ggpredict(t_mod_r3,"dbh[all]",type="random")

tr_dbh<-ggplot() +
  geom_ribbon(tr_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high),fill="grey90")+
  geom_line(tr_dbh_preds, mapping=aes(x, predicted)) +
  #geom_line(ot_r_dbh_preds, mapping=aes(x, predicted),lty=2,alpha=0.5) +
  #geom_ribbon(t_mean_ba_preds_quad,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .1)+
  #geom_ribbon(ot_r_dbh_preds,mapping=aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high), alpha = .05)+
  geom_point(tree_data_,mapping=aes(x=dbh,y=tree_r),alpha=0.2)+
  labs(x="DBH",y="Synchrony")+
  annotate("text",label=expression(italic("P")~"= 0.08"),x=70,y=1.0,parse=TRUE, hjust=0,size=3)+
  theme_classic()


tr_dbh


(mean_dbh|tmean_ba)/(tcv_dbh|tcv_ba)/(tr_dbh|tr_ba)














#map

pcc_df<-as.data.frame(western_pipo,xy=T)

red_sf<-as.data.frame(s_df[,c(1,24,25)])
names(red_sf)<-c('site','x','y')

#library(ggplot)

extent(western_pipo)

west_lat<-spTransform(west,CRS=crs(sites))
west_latc<-crop(west_lat,sites)
westlat_df<-fortify(west_latc)

plot(west_lat)


cwd.df<-as.data.frame(pipo_cwdmask,xy=T)
aet.df<-as.data.frame(pipo_aetmask,xy=T)

extent(west)
extent(western_pipo)


ext<-c(-800000.6,-770505.2,1803827,1843242)

pipo_crop<-crop(western_pipo,ext)
pipodf<-as.data.frame(pipo_crop,xy=TRUE)

plot(pipo_crop)

thumap<-ggplot()+
  #geom_path(rd_c,mapping=aes(x=long,y=lat,group=group),size=1.1,alpha=0.5)+
  geom_raster(elev.df,mapping=aes(x=x,y=y,fill=west_elev),alpha=0.5)+
  scale_fill_scico(palette = 'vik',na.value=NA,direction=-1)+
  #geom_point(hm.df,mapping=aes(x=coords.x1,y=coords.x2),size=3)+
  geom_text(hm.df,mapping=aes(x=coords.x1,y=coords.x2,label=site))+
  theme_void()+
  theme(legend.position="none")

hm1<-subset(treeyr_data,site=="HM1")
hm2<-subset(treeyr_data,site=="HM2")
hm3<-subset(treeyr_data,site=="HM3")
hm4<-subset(treeyr_data,site=="HM4")
hm5<-subset(treeyr_data,site=="HM5")

pal

hm1gg<-ggplot()+
  geom_line(hm1,mapping=aes(x=year,y=cone_total/leadershoots,col=tree),size=1.25,alpha=0.5)+
  annotate("text",label=expression(("B")),x=2003,y=2.5,parse=TRUE, hjust=0,size=8)+
  geom_hline(yintercept=2,lty=2,alpha=0.75,col="gray50",size=.5)+
  geom_hline(yintercept=1,lty=2,alpha=0.33,col="gray50",size=.5)+
  annotate("text",label=expression(bold("HM1")),x=2017,y=3,parse=TRUE, hjust=0,size=5)+
  annotate("text",label=expression(italic("r")~"= 0.36"),x=2017,y=1.25,parse=TRUE, hjust=0,size=3.5)+
  annotate("text",label=expression(CV[p]~"= 1.40"),x=2017,y=1.7, hjust=0,size=3.5)+  
  scale_color_manual(values=pall)+
  #geom_vline(xintercept=2005,lty=2,alpha=0.5,col="blue",size=1.25)+
  #geom_vline(xintercept=2010,lty=2,alpha=0.5,col="blue",size=1.25)+
  xlim(c(2003,2019))+
  theme_classic()+
  labs(x='',y='')+
  theme(legend.position="none", axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

hm2gg<-ggplot()+
  geom_line(hm2,mapping=aes(x=year,y=cone_total/leadershoots,col=tree),size=1.25,alpha=0.5)+
  xlim(c(2003,2019))+
  geom_hline(yintercept=2,lty=2,alpha=0.75,col="gray50",size=.5)+
  geom_hline(yintercept=1,lty=2,alpha=0.33,col="gray50",size=.5)+
  #geom_vline(xintercept=2005,lty=2,alpha=0.5,col="blue",size=1.25)+
  #geom_vline(xintercept=2010,lty=2,alpha=0.5,col="blue",size=1.25)+
  annotate("text",label=expression(bold("HM2")),x=2017,y=3,parse=TRUE, hjust=0,size=5)+
  annotate("text",label=expression(italic("r")~"= 0.17"),x=2017,y=1.2,parse=TRUE, hjust=0,size=3.5)+
  annotate("text",label=expression(CV[p]~"= 1.34"),x=2017,y=1.75, hjust=0,size=3.5)+    
  scale_color_manual(values=pall)+
  theme_classic()+
  labs(x='',y='')+
  theme(legend.position="none", axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

hm3gg<-ggplot()+
  geom_line(hm3,mapping=aes(x=year,y=cone_total/leadershoots,col=tree),size=1.25,alpha=0.5)+
  xlim(c(2003,2019))+
  geom_hline(yintercept=2,lty=2,alpha=0.75,col="gray50",size=.5)+
  geom_hline(yintercept=1,lty=2,alpha=0.33,col="gray50",size=.5)+
  #geom_vline(xintercept=2005,lty=2,alpha=0.5,col="blue",size=1.25)+
  #geom_vline(xintercept=2010,lty=2,alpha=0.5,col="blue",size=1.25)+
  annotate("text",label=expression(bold("HM3")),x=2017,y=3,parse=TRUE, hjust=0,size=5)+
  annotate("text",label=expression(italic("r")~"= 0.59"),x=2017,y=1.2,parse=TRUE, hjust=0,size=3.5)+
  annotate("text",label=expression(CV[p]~"= 2.66"),x=2017,y=1.75, hjust=0,size=3.5)+  
  scale_color_manual(values=pall)+
  theme_classic()+
  labs(x='',y='Std. Cone Production')+
  theme(legend.position="none", axis.title.y=element_text(size=12),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())


hm4gg<-ggplot()+
  geom_line(hm4,mapping=aes(x=year,y=cone_total/leadershoots,col=tree),size=1.25,alpha=0.5)+
  xlim(c(2003,2019))+
  geom_hline(yintercept=2,lty=2,alpha=0.75,col="gray50",size=.5)+
  geom_hline(yintercept=1,lty=2,alpha=0.33,col="gray50",size=.5)+
  #geom_vline(xintercept=2005,lty=2,alpha=0.5,col="blue",size=1.25)+
  #geom_vline(xintercept=2010,lty=2,alpha=0.5,col="blue",size=1.25)+
  annotate("text",label=expression(bold("HM4")),x=2017,y=3,parse=TRUE, hjust=0,size=5)+
  annotate("text",label=expression(italic("r")~"= 0.55"),x=2017,y=1.2,parse=TRUE, hjust=0,size=3.5)+
  annotate("text",label=expression(CV[p]~"= 1.45"),x=2017,y=1.75, hjust=0,size=3.5)+
  scale_color_manual(values=pall)+
  theme_classic()+
  labs(x='Year',y='')+
  theme(legend.position="none")

hm5gg<-ggplot()+
  geom_line(hm5,mapping=aes(x=year,y=cone_total/leadershoots,col=tree),size=1.25,alpha=0.5)+
  xlim(c(2003,2019))+
  geom_hline(yintercept=2,lty=2,alpha=0.75,col="gray50",size=.5)+
  geom_hline(yintercept=1,lty=2,alpha=0.33,col="gray50",size=.5)+
  #geom_vline(xintercept=2005,lty=2,alpha=0.5,col="blue",size=1.25)+
  #geom_vline(xintercept=2010,lty=2,alpha=0.5,col="blue",size=1.25)+
  annotate("text",label=expression(bold("HM5")),x=2017,y=3.5,parse=TRUE, hjust=0,size=5)+
  annotate("text",label=expression(italic("r")~"= 0.53"),x=2017,y=1.2,parse=TRUE, hjust=0,size=3.5)+
  annotate("text",label=expression(CV[p]~"= 1.92"),x=2017,y=1.75, hjust=0,size=3.5)+
  scale_color_manual(values=(pall))+
  theme_classic()+
  labs(x='Year',y='')+
  theme(legend.position="none")


hm1gg/hm2gg/hm3gg/hm4gg/hm5gg




library(ggspatial)

western_pipo.df<-as.data.frame(western_pipo,xy=T)

pcc_df

cwd_aet<-ggplot()+
  geom_raster(pcc_df,mapping=aes(x=x,y=y,fill=PIPO_1km),alpha=0.7)+
  geom_path(west_df,mapping=aes(x=long,y=lat,group=id),size=0.75,alpha=0.5)+
  #geom_rect(mapping=aes(xmin=-750000.6 ,xmax= -800505.2 ,ymin=1773242.0,ymax=1853827.0),lty=2,size=1,fill=NA,col='black')+
  scale_fill_scico(palette = "tofino",begin=0.8,end=0.8,na.value=NA)+
  #coord_cartesian(xlim=c(-1418833 ,-591704.6 ),ylim=c(1447632 ,2430060 ))+
  coord_cartesian(xlim=c(-1476717,-507717),ylim=c(1051792,2489792))+
  #annotate("text",label=expression("(A)"),x=-1476717,y=2489792,parse=TRUE, hjust=0,size=8)+
  labs(x="",y="")+
  geom_point(red_sf,mapping=aes(x=x,y=y,alpha=0.66),position=position_jitter(w=3,h=3),size=4,shape=16)+
  theme_classic()+
  theme(axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),legend.position = "none")+
  annotation_scale(line_width = 1,height = unit(0.5,"cm"))

cwd_aet

#png(file="C:/Users/Andreas/Google Drive/Cones_Ponderosa/figs/map.png",width=10, height=12.5, bg="white", units="cm", res= 300)
cwd_aet
dev.off()

#png(file="C:/Users/Andreas/Google Drive/Cones_Ponderosa/figs/map with mast.png",width=15, height=15, bg="white", units="cm", res= 300)
hm2gg/hm3gg/hm4gg
dev.off()



###




















(map|cwd_aet|cwd_aet)
dev.off()
plot(western_pipo)

crs(aea_sites)

us
library(cowplot)
library(sf)
library(sp)

aea_sites
box<-st_as_sfc(st_bbox(extent(aea_sites),crs=5070))
data("us_states", package = "spData")
us = st_transform(us_states, crs = 5070)

#a_crs<-crs(aea_sites)
crs(us)
crs(box)



ggm1 = ggplot() + 
  geom_sf(data = us, fill = "lightgray",) + 
  geom_sf(data = box, fill = NA, color = "red", size = 1,lty=1) +
  theme_void()

ggm1

gg_inset_map1 = ggdraw() +
  draw_plot(cwd_aet) +
  draw_plot(ggm1, x = 0.15, y = 0.66, width = 0.4, height = 0.4)

#png(file="C:/Users/Andreas/Google Drive/Cones_Ponderosa/figs/map.png",width=11.5, height=15, bg="white", units="cm", res= 300)
gg_inset_map1
dev.off()

map<-ggplot(red_sf,mapping=aes(x=x,y=y),col="black",size=3,shape=20)+
  coord_cartesian(xlim=c(-1418833 ,-591704.6 ),ylim=c(1447632 ,2430060 ))+
  geom_path(west_df,mapping=aes(x=long,y=lat,group=id),size=0.75,alpha=0.5)+
  geom_jitter(width = 100, height=100)

map



map<-ggplot()+
  geom_path(west_df,mapping=aes(x=long,y=lat,group=id),size=0.75,alpha=0.5)+
  coord_cartesian(xlim=c(-1706717,-357717),ylim=c(1400000,2419792))+
  geom_raster(pcc_df,mapping=aes(x=x,y=y,fill=PIPO_1km),alpha=0.75,show.legend=FALSE)+
  scale_fill_scico(palette = "imola",na.value=NA,"")+
  geom_point(red_sf,mapping=aes(x=x,y=y),size=3,shape=20)+
  geom_jitter()+
  theme(legend.position = "none")+
  labs(x="x",y="y")+
  theme_void()



mapp<-ggplot()+
  geom_raster(pcc_df,mapping=aes(x=x,y=y,fill=PIPO_1km))+
  theme(legend.position = "none")+
  theme_classic()


mapp


#png(file="C:/Users/Andreas/Google Drive/Cones_Ponderosa/figs/map.png",width=15, height=15, bg="white", units="cm", res= 300)
map
dev.off()
pipo




map<-ggplot()+
  geom_point(red_sf,mapping=aes(x=x,y=y),position=position_jitter(w=.25,h=.25),size=3,shape=20)+
  geom_path(west_df,mapping=aes(x=long,y=lat,group=id),size=0.75,alpha=0.5)+
  coord_cartesian(xlim=c(-1506717,-357717),ylim=c(900000,2719792))+
  geom_raster(pcc_df,mapping=aes(x=x,y=y,fill=PIPO_1km),alpha=0.75,show.legend=FALSE)+
  scale_fill_scico(palette = "imola",na.value=NA,"")+
  geom_jitter()+
  theme(legend.position = "none")+
  labs(x="x",y="y")+
  theme_void()
map


cwd_map<-ggplot()+
  geom_raster(cwd.df,mapping=aes(x=x,y=y,fill=CWD),alpha=0.7)+
  geom_path(west_df,mapping=aes(x=long,y=lat,group=id),size=0.75,alpha=0.5)+
  scale_fill_scico(palette = "vikO",na.value=NA,direction=1)+
  coord_cartesian(xlim=c(-1506717,-607717),ylim=c(1400000,2419792))+
  geom_point(red_sf,mapping=aes(x=x,y=y),size=3,shape=3)+
  #coord_cartesian(xlim=c(-1706717,-357717),ylim=c(991792,2419792))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x="",y="")+
  theme_classic()

cwd_aet<-ggplot()+
  geom_raster(aet.df,mapping=aes(x=x,y=y,fill=AET),alpha=0.7)+
  geom_path(west_df,mapping=aes(x=long,y=lat,group=id),size=0.75,alpha=0.5)+
  scale_fill_scico(palette = "vikO",na.value=NA,direction=-1)+
  coord_cartesian(xlim=c(-1506717,-607717),ylim=c(1400000,2419792))+
  geom_point(red_sf,mapping=aes(x=x,y=y),size=2,shape=16)+
  #coord_cartesian(xlim=c(-1706717,-357717),ylim=c(991792,2419792))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(x="",y="")+
  theme_classic()

cwd_aet

#png(file="C:/Users/Andreas/Google Drive/Cones_Ponderosa/figs/map.png",width=25, height=10, bg="white", units="cm", res= 300)
(cwd_map|cwd_aet)
dev.off()


##







mast1<-lm(CVi~tree_r,tree_data)









