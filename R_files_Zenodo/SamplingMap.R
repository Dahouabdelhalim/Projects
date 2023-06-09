#### Create sampling map for ECP species delimitation project #####
### Libraries ###
require(raster)
require(rgdal)
require(sp)
require(xlsx)
require(maps)
require(rgeos)
require(adegenet)
require(plotrix)

### Read in meta data for sampling sites ###
metadata<-read.xlsx("~/Dropbox/InsectTaxa/metadata/insecttaxa_metadata.xlsx",sheetName="Andesiops peruvianus (Ape)",stringsAsFactors=F)
head(metadata)
ECP_data<-metadata[grep("ECP",metadata$FULL),]
ECA_data<-metadata[grep("ECA",metadata$FULL),]

ECP_df<-data.frame(lat=ECP_data$LAT,lon=ECP_data$LON,row.names=ECP_data$FULL,stringsAsFactors=F)
ECP_df<-ECP_df[!duplicated(ECP_df[,1]),]
rownames(ECP_df)<-sapply(strsplit(rownames(ECP_df),"_"),function(x) x[3])

ECA_df<-data.frame(lat= ECA_data $LAT,lon= ECA_data $LON,row.names= ECA_data $FULL,stringsAsFactors=F)
ECA_df<-ECA_df[!duplicated(ECA_df[,1]),]
rownames(ECA_df)<-sapply(strsplit(rownames(ECA_df),"_"),function(x) x[3])

coords_df<-rbind(ECP_df,ECA_df[1,])
#coords_df<-coords_df[-c(1,2,12,13),]

### Read in K2 structure q values ###
ecp<-read.genepop("~/Dropbox/InsectTaxa/genepop/ApeECP.gen")
ecp_k2<-readLines(file("~/Dropbox/InsectTaxa/dapc_and_structure/Ape_ECP_Struc/StructureHarvesterOutput 2/K2_CLUMPP.OUT"))

q<-sapply(strsplit(ecp_k2,":"),function(x) x[length(x)])
q_df<-data.frame(lapply(strsplit(q," "),function(x) x[3:4]),stringsAsFactors=F)
colnames(q_df)<-NULL
ecp_k2_q_df<-as.matrix(q_df)
Q_vec<-as.numeric(ecp_k2_q_df[1,])
names(Q_vec)<-row.names(ecp@tab)

q_df<-data.frame(q=Q_vec,site=sapply(strsplit(names(Q_vec),"_"),function(x) x[3]),stringsAsFactors=F)
q_list<-split(q_df,q_df$site)

mean_q<-sapply(q_list,function(x) mean(x$q))
mean_q2<-data.frame(q=mean_q,q2=1-mean_q,row.names=names(mean_q),stringsAsFactors=F)
mean_q2<-mean_q2[-c(1,2,12,13),]

### Download elevational data ###
ecu_ele1<-raster("~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/Data/GIS/ASTGTM2_S01W078/ASTGTM2_S01W078_dem.tif")
ecu_ele2<-raster("~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/Data/GIS/ASTGTM2_S01W079/ASTGTM2_S01W079_dem.tif")
ecu_ele<-merge(ecu_ele1,ecu_ele2)

### Read in water shapefile ### LARGE FILE TAKES A LONG TIME TO READ IN
ecu_wat<-shapefile("~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/Data/GIS/hydrosheds-f2bbc7371b8c08308119/sa_riv_15s/sa_riv_15s.shp")

site_bbox<-as(extent(bbox(as.matrix(coords_df[c(2,1)]))),"SpatialPolygons")
proj4string(site_bbox)<-proj4string(ecu_wat)

ecu_wat_crop<-crop(ecu_wat,site_bbox)

### Create map ###
#quartz(width=6.5,height=6.5*0.7)
pdf(file="~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/Figures/samplingmap_v1.pdf")
#layout(matrix(c(1:2),nrow=1),widths=c(0.7,0.3))
#par(mar=c(4,4,1,1))
plot(coords_df[,2],coords_df[,1],xlab="Latitude",ylab="Longitude")
plot(ecu_ele,add=T,col=colorRampPalette(c("black","white"))(100))
plot(ecu_wat_crop,col="blue",add=T,lwd=2)

for(i in 1:nrow(mean_q2)){
	floating.pie(xpos=coords_df[rownames(mean_q2),2][i],ypos= coords_df[rownames(mean_q2),1][i],x=as.numeric(mean_q2[i,]),radius=0.01,col=c("red","blue"))
}

points(coords_df[nrow(coords_df),2],coords_df[nrow(coords_df),1],pch=21,bg="orange",cex=4)

dev.off()

0.008998329 #degrees to km at this latitude

map(xlim=c(-84,-35),ylim=c(-60,15))

png(file="~/Desktop/PeruWithBox.png",res=500,width=4,height=4,units="in")
par(xpd=NA)
map(regions="Ecuador",exact=T,mar=c(0,0,0,0))
rect(bbox(as.matrix(coords_df[c(2,1)]))[1],bbox(as.matrix(coords_df[c(2,1)]))[2],bbox(as.matrix(coords_df[c(2,1)]))[3],bbox(as.matrix(coords_df[c(2,1)]))[4],border="red",lwd=1,col="transparent")
dev.off()

sa_bbox<-locator(2)

### Create south america map with lil square for region of interest ###
map(xlim=c(-90,-20),ylim=c())