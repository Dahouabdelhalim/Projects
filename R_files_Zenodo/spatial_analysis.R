#spatial analyses

#open packages
library(dplyr)
library(geoR)
library(raster)
library(ggplot2)
library(SpatialEpi)
library(extrafont)
loadfonts(device="win")


#remove old data
rm(list=ls())

#import data
#import lightning frequency data from ENTLN (0.05 x 0.05 degree cells) as "fulldata_meanFRD"
str(fulldata_meanFRD)


#####################################################
#drop water bodies, ice/snow, and needleleaf forests
logical1<-fulldata_meanFRD$Type.0.Coverage+fulldata_meanFRD$Type.1.Coverage+fulldata_meanFRD$Type.3.Coverage+fulldata_meanFRD$Type.15.Coverage==100
trimdata<-fulldata_meanFRD[!logical1,]
str(trimdata)

#estimate x-y position based on cell area rather than just degrees
latlon_coords<-as.matrix(data.frame(trimdata[,c(2,1)]))
latlon_coords_distance<-latlong2grid(latlon_coords)
#add the disance values to the dataframe
trimdata$longitude_km<-latlon_coords_distance$x
trimdata$latitude_km<-latlon_coords_distance$y

###################################################
#Americas

#create a variogram for the largest continuous area in the Americas
#subset to only include the americas
Americas_data<-trimdata[trimdata$Continent=="Americas",]
str(Americas_data)

#create a version of the coordinates that are at the ca. 1 degree scale
#create 1 degree bins
breaks_lat_1deg_cut<-seq(-23.5,23.5,1)
Americas_data$latitude_1degbins<-cut(Americas_data$Latitude,breaks_lat_1deg_cut,labels = seq(-23,23,1))
breaks_long_1deg_cut<-seq(-180,180,1)
Americas_data$longitude_1degbins<-cut(Americas_data$Longitude,breaks_long_1deg_cut, labels = seq(-179.5,179.5,1))
Americas_data$bins_1deg<-paste(Americas_data$longitude_1degbins,Americas_data$latitude_1degbins,sep=".")
length(levels(as.factor(Americas_data$bins_1deg)))

#make the new lat-lon values numeric vectors
Americas_data$longitude_1degbins<-as.numeric(as.character(Americas_data$longitude_1degbins))
Americas_data$latitude_1degbins<-as.numeric(as.character(Americas_data$latitude_1degbins))

#then aggregate based on that scale of the coordinates
aggregated_Americas<-Americas_data%>%
  group_by(bins_1deg)%>%
  summarise(Longitude = mean(longitude_1degbins),
            Latitude = mean(latitude_1degbins),
            Longitude_km = mean(longitude_km),
            Latitude_km = mean(latitude_km),
            meanFRD = mean(CG...10kA.FRD))

str(aggregated_Americas)

#create a raster object from coordinates and lightning values
Americas_XYZ<-aggregated_Americas[,c(2,3,4,5,6)]
Americas_XYZ_raster<-rasterFromXYZ(Americas_XYZ)#this centers each cell on the x-y coordinates (so minimum Latitude is now -23.475 instead of -23.5)
plot(Americas_XYZ_raster)

#create a polygon for sampling 
Americas_coords<-as.matrix(data.frame(c(-60,-36,-60,-80),c(7,-9,-23.5,-5)))#minimum distance along an edge is 23.3 degrees
Americas_poly<-spPolygons(Americas_coords)
plot(Americas_XYZ_raster)
plot(Americas_poly, add=T)

#now mask all values outside the polygon
Americas_raster_extracted<-mask(Americas_XYZ_raster,Americas_poly)
plot(Americas_raster_extracted)
plot(Americas_poly, add=T)

#turn the extracted data into a dataframe
Americas_square_full<-as.data.frame(Americas_raster_extracted, xy = TRUE)
Americas_square<-Americas_square_full[!is.na(Americas_square_full$meanFRD),]
Americas_square

#variogram bins - max distance is ca. 2500 km, so we will use half of that
deg_vec = seq(100,1300,100)

#create a GEOR object using the "americas_square" data
str(Americas_square)
Americas_geoR_dat<-as.geodata(Americas_square,coords.col=3:4,data.col=5)

#create semivariogram
Americas_geoR_variog<-variog(Americas_geoR_dat,max.dist = 1300,uvec = deg_vec)

#plot the semivariogram
plot(Americas_geoR_variog)

#create a df for ploting
df_Americas<-data.frame("distance_km"=Americas_geoR_variog$u)
df_Americas$semivariance<-Americas_geoR_variog$v
df_Americas[14,]<-c(0,0)#create a point at zero for plotting

#now plot these values
#create the theme
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text=element_text(family = "Arial",colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour="black"))

#now creat figure
Variogram_americas<-ggplot(data = df_Americas, aes(x = distance_km, y = semivariance))+
  geom_point(size = 2)+
  geom_smooth(method="loess",se=FALSE, color = "black")+
  scale_x_continuous(name = "Distance (km)",breaks = seq(0,1300,100),labels = c("0",rep("",4),"500",rep("",4),"1000",rep("",3)))+
  scale_y_continuous(name = "Semivariance",breaks = seq(0,1.5,0.1),labels = c("0",rep("",4),"0.5",rep("",4),"1.0",rep("",4),"1.5"))+
  theme_basis

Variogram_americas

ggsave("Variogram_americas.tiff",Variogram_americas,dpi = 600, width = 3.5,height = 2.25,scale = 1.8, compression = "lzw")

###################################################
#now create a similar figure at a small spatial scale (<100km)
str(Americas_data)

#create a raster object from coordinates and lightning values
Americas_XYZ<-Americas_data[,c(2,1,28,29,25)]
Americas_XYZ_raster<-rasterFromXYZ(Americas_XYZ)#this centers each cell on the x-y coordinates (so minimum Latitude is now -23.475 instead of -23.5)
plot(Americas_XYZ_raster)

#create a polygon for sampling 
Americas_coords<-as.matrix(data.frame(c(-60,-36,-60,-80),c(7,-9,-23.5,-5)))#minimum distance along an edge is 23.3 degrees
Americas_poly<-spPolygons(Americas_coords)
plot(Americas_XYZ_raster)
#plot(Americas_poly, add=T) # this code can be used to creat polygon of selected area if only 1 layer is selected (e.g.,only CG lightning)

#now mask all values outside the polygon
Americas_raster_extracted<-mask(Americas_XYZ_raster,Americas_poly)
plot(Americas_raster_extracted)

#turn the extracted data into a dataframe
Americas_square_full<-as.data.frame(Americas_raster_extracted, xy = TRUE)
Americas_square<-Americas_square_full[!is.na(Americas_square_full$CG...10kA.FRD),]
Americas_square

#variogram for loop to estimate confidence with manageable dataset sizes
deg_vec = seq(0,100,5)

df_americas_local<-data.frame("distance_degrees" = rep(NA,10*20),"semivariance" = rep(NA,10*20),"iteration" = rep(NA,10*20))
for(i in 1:1000){
  #create a for loop to bootstrap variogram parameters
  #randomly subset americas_local data to a size and object type that we can analyze
  americas_local_square_subset<-sample_n(Americas_square,10000,replace=FALSE)
  americas_local_geoR_dat<-as.geodata(americas_local_square_subset,coords.col=3:4,data.col=5)
  
  #create variogram
  americas_local_geoR_variog<-variog(americas_local_geoR_dat,max.dist = 100,uvec = deg_vec)
  
  #variogram locations for saving data
  rowvals<-seq(i*20-19,i*20,1)
  
  #save variogram characteristics
  df_americas_local[rowvals,"distance_degrees"]<-americas_local_geoR_variog$u
  df_americas_local[rowvals,"semivariance"]<-americas_local_geoR_variog$v
  df_americas_local[rowvals,"iteration"]<-rep(paste("iter",i,sep="_"),20)
}

#View(df_americas_local)
#write.csv(df_americas_local,"americas_local_semivariance_output.csv")

#calculate total area we are measuring over in the americas_local
Americas_square$cell_ID<-paste(round(Americas_square$x,digits = 2),round(Americas_square$y, digits = 2), sep = ".")
Americas_data$cell_ID<-paste(Americas_data$longitude_1degbins,Americas_data$latitude_1degbins, sep = ".")

#combine datasets
americas_local_area_calc<-left_join(Americas_square,Americas_data[,c("cell_ID","bin.size.km2","Continent")], by = "cell_ID")
sum(americas_local_area_calc$bin.size.km2,na.rm=TRUE)
nrow(americas_local_area_calc)

#now calculate mean, hiCI and loCI
#create dataframe to save values
americas_local_semi_df<-data.frame("distance_km" = rep(NA,20),"mean" = rep(NA,20),"lo_CI" = rep(NA,20),"hi_CI" = rep(NA,20))
#now use for loop
for(i in 1:20){
  semi_vec<-df_americas_local[df_americas_local$distance_degrees==(i*5),"semivariance"]
  americas_local_semi_df[i,"distance_degrees"]<-i*5-2.5
  americas_local_semi_df[i,"mean"]<-mean(semi_vec)
  americas_local_semi_df[i,"lo_CI"]<-quantile(semi_vec, 0.025)
  americas_local_semi_df[i,"hi_CI"]<-quantile(semi_vec, 0.975)
}
americas_local_semi_df
str(americas_local_semi_df)
#rename distance_degrees as distance_km - it is currently misnamed
americas_local_semi_df$distance_km<-americas_local_semi_df$distance_degrees

#add a zero-zero point for plotting
americas_local_semi_df[21,]<-c(0,0,0,0,0)

#now plot these values
#create the theme
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text=element_text(family = "Arial",colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour="black"))

#now creat figure
Variogram_americas_local<-ggplot(data = americas_local_semi_df, aes(x = distance_km, y = mean))+
  geom_point(size = 2)+
  geom_errorbar(data = americas_local_semi_df, aes(ymin = lo_CI,ymax = hi_CI))+
  geom_smooth(method= "loess",se=FALSE, color = "black")+
  scale_x_continuous(name = "Distance (km)",breaks = seq(0,100,5),labels = c("0",rep("",4),"25",rep("",4),"50",rep("",4),"75",rep("",4),"100"))+
  scale_y_continuous(name = "Semivariance",breaks = seq(0,1.5,0.25),labels = c("0","","0.5","","1.0","","1.5"))+
  theme_basis

Variogram_americas_local

ggsave("Variogram_americas_local.tiff",Variogram_americas_local,dpi = 600, width = 3.5,height = 2.25,scale = 1.8, compression = "lzw")

###################################################
#Africa - all code identical to above, except within a different set of coordinates

#create a variogram for the largest continuous area in the Africa
#subset to only include the africa
Africa_data<-trimdata[trimdata$Continent=="Africa",]
str(Africa_data)

#create a version of the coordinates that are at the 1 degree scale
#also create 1 degree bins
breaks_lat_1deg_cut<-seq(-23.5,23.5,1)
Africa_data$latitude_1degbins<-cut(Africa_data$Latitude,breaks_lat_1deg_cut,labels = seq(-23,23,1))
breaks_long_1deg_cut<-seq(-180,180,1)
Africa_data$longitude_1degbins<-cut(Africa_data$Longitude,breaks_long_1deg_cut, labels = seq(-179.5,179.5,1))
Africa_data$bins_1deg<-paste(Africa_data$longitude_1degbins,Africa_data$latitude_1degbins,sep=".")
length(levels(as.factor(Africa_data$bins_1deg)))

#make the new lat-lon values numeric vectors
Africa_data$longitude_1degbins<-as.numeric(as.character(Africa_data$longitude_1degbins))
Africa_data$latitude_1degbins<-as.numeric(as.character(Africa_data$latitude_1degbins))

#then aggregate based on that scale of the coordinates
aggregated_Africa<-Africa_data%>%
  group_by(bins_1deg)%>%
  summarise(Longitude = mean(longitude_1degbins),
            Latitude = mean(latitude_1degbins),
            Longitude_km = mean(longitude_km),
            Latitude_km = mean(latitude_km),
            meanFRD = mean(CG...10kA.FRD))

str(aggregated_Africa)

#create a raster object from coordinates and lightning values
Africa_XYZ<-aggregated_Africa[,c(2,3,4,5,6)]
Africa_XYZ_raster<-rasterFromXYZ(Africa_XYZ)#this centers each cell on the x-y coordinates (so minimum Latitude is now -23.475 instead of -23.5)
plot(Africa_XYZ_raster)

#create a polygon for sampling 
Africa_coords<-as.matrix(data.frame(c(39,15,14,38),c(-17,-17,16,16)))#minimum distance along an edge is 24 degrees
Africa_poly<-spPolygons(Africa_coords)
plot(Africa_XYZ_raster)
plot(Africa_poly, add=T)

#now mask all values outside the polygon
Africa_raster_extracted<-mask(Africa_XYZ_raster,Africa_poly)
plot(Africa_raster_extracted)
plot(Africa_poly, add=T)

#turn the extracted data into a dataframe
Africa_square_full<-as.data.frame(Africa_raster_extracted, xy = TRUE)
Africa_square<-Africa_square_full[!is.na(Africa_square_full$meanFRD),]
Africa_square

#variogram bins - max distance on any edge is ca. 2500 km, so we will use half of that
deg_vec = seq(100,1300,100)

#create a GEOR object using the "africa_square" data
str(Africa_square)
Africa_geoR_dat<-as.geodata(Africa_square,coords.col=3:4,data.col=5)

#create semivariogram
Africa_geoR_variog<-variog(Africa_geoR_dat,max.dist = 1300,uvec = deg_vec)

#plot the semivariogram
plot(Africa_geoR_variog)

#create a df for ploting
df_Africa<-data.frame("distance_km"=Africa_geoR_variog$u)
df_Africa$semivariance<-Africa_geoR_variog$v
df_Africa[14,]<-c(0,0)#create 0-0 coordinate for plotting

#now plot these values
#create the theme
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text=element_text(family = "Arial",colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour="black"))

#now creat figure
Variogram_africa<-ggplot(data = df_Africa, aes(x = distance_km, y = semivariance))+
  geom_point(size = 2)+
  geom_smooth(method="loess",se=FALSE, color = "black")+
  scale_x_continuous(name = "Distance (km)",breaks = seq(0,1300,100),labels = c("0",rep("",4),"500",rep("",4),"1000",rep("",3)))+
  scale_y_continuous(name = "Semivariance",breaks = seq(0,8,1),labels = c("0","","2","","4","","6","","8"))+
  theme_basis

Variogram_africa

ggsave("Variogram_africa.tiff",Variogram_africa,dpi = 600, width = 3.5,height = 2.25,scale = 1.8, compression = "lzw")

###################################################
#now create the same figure at a small spatial scale (<100km)
str(Africa_data)

#create a raster object from coordinates and lightning values
Africa_XYZ<-Africa_data[,c(2,1,28,29,25)]
Africa_XYZ_raster<-rasterFromXYZ(Africa_XYZ)#this centers each cell on the x-y coordinates (so minimum Latitude is now -23.475 instead of -23.5)
plot(Africa_XYZ_raster)

#create a polygon for sampling 
Africa_coords<-as.matrix(data.frame(c(39,15,14,38),c(-17,-17,16,16)))#minimum distance along an edge is 24 degrees
Africa_poly<-spPolygons(Africa_coords)
plot(Africa_XYZ_raster)
#plot(Africa_poly, add=T) # this code can be used to creat polygon of selected area if only 1 layer is selected (e.g.,only CG lightning)

#now mask all values outside the polygon
Africa_raster_extracted<-mask(Africa_XYZ_raster,Africa_poly)
plot(Africa_raster_extracted)

#turn the extracted data into a dataframe
Africa_square_full<-as.data.frame(Africa_raster_extracted, xy = TRUE)
Africa_square<-Africa_square_full[!is.na(Africa_square_full$CG...10kA.FRD),]
Africa_square

#variogram for loop to estimate confidence with manageable dataset sizes
deg_vec = seq(0,100,5)

df_africa_local<-data.frame("distance_degrees" = rep(NA,10*20),"semivariance" = rep(NA,10*20),"iteration" = rep(NA,10*20))
for(i in 1:1000){
  #create a for loop to bootstrap variogram parameters
  #randomly subset africa_local data to a size and object type that we can analyze
  africa_local_square_subset<-sample_n(Africa_square,10000,replace=FALSE)
  africa_local_geoR_dat<-as.geodata(africa_local_square_subset,coords.col=3:4,data.col=5)
  
  #create variogram
  africa_local_geoR_variog<-variog(africa_local_geoR_dat,max.dist = 100,uvec = deg_vec)
  
  #variogram locations for saving data
  rowvals<-seq(i*20-19,i*20,1)
  
  #save variogram characteristics
  df_africa_local[rowvals,"distance_degrees"]<-africa_local_geoR_variog$u
  df_africa_local[rowvals,"semivariance"]<-africa_local_geoR_variog$v
  df_africa_local[rowvals,"iteration"]<-rep(paste("iter",i,sep="_"),20)
}

View(df_africa_local)
#write.csv(df_africa_local,"africa_local_semivariance_output.csv")

#calculate total area we are measuring over in the africa_local
Africa_square$cell_ID<-paste(round(Africa_square$x,digits = 2),round(Africa_square$y, digits = 2), sep = ".")
Africa_data$cell_ID<-paste(Africa_data$Longitude,Africa_data$Latitude, sep = ".")

#combine datasets
africa_local_area_calc<-left_join(Africa_square,Africa_data[,c("cell_ID","bin.size.km2","Continent")])
sum(africa_local_area_calc$bin.size.km2,na.rm=TRUE)
nrow(africa_local_area_calc)

#now calculate mean, hiCI and loCI
#create dataframe to save values
africa_local_semi_df<-data.frame("distance_degrees" = rep(NA,20),"mean" = rep(NA,20),"lo_CI" = rep(NA,20),"hi_CI" = rep(NA,20))
#now use for loop
for(i in 1:20){
  semi_vec<-df_africa_local[df_africa_local$distance_degrees==(i*5),"semivariance"]
  africa_local_semi_df[i,"distance_degrees"]<-i*5-2.5
  africa_local_semi_df[i,"mean"]<-mean(semi_vec)
  africa_local_semi_df[i,"lo_CI"]<-quantile(semi_vec, 0.025)
  africa_local_semi_df[i,"hi_CI"]<-quantile(semi_vec, 0.975)
}
africa_local_semi_df
str(africa_local_semi_df)

#rename distance_degrees
colnames(africa_local_semi_df)[1]<-"distance_km"

#add a zero-zero point for plotting
africa_local_semi_df[21,]<-c(0,0,0,0)

#now plot these values
#create the theme
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text=element_text(family = "Arial",colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour="black"))

#now creat figure
Variogram_africa_local<-ggplot(data = africa_local_semi_df, aes(x = distance_km, y = mean))+
  geom_point(size = 2)+
  geom_errorbar(data = africa_local_semi_df, aes(ymin = lo_CI,ymax = hi_CI))+
  geom_smooth(method = "loess",color = "black", se = F)+
  scale_x_continuous(name = "Distance (km)",breaks = seq(0,100,5),labels = c("0",rep("",4),"25",rep("",4),"50",rep("",4),"75",rep("",4),"100"))+
  scale_y_continuous(name = "Semivariance",breaks = seq(0,.75,0.05),labels = c("0.0",rep("",3),"0.2",rep("",3),"0.4",rep("",3),"0.6",rep("",3)))+
  theme_basis

Variogram_africa_local

ggsave("Variogram_africa_local.tiff",Variogram_africa_local,dpi = 600, width = 3.5,height = 2.25,scale = 1.8, compression = "lzw")



