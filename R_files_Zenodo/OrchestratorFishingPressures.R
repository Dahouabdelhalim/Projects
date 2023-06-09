rm(list=ls(all=TRUE))
#cat("v19 02 2021\\n")
cat("v26 07 2021\\n")

suppressMessages(library(raster))
#suppressMessages(library(cowplot))
#suppressMessages(library(ggplot2))
#suppressMessages(library(DT))


res<-1
##input data
xcolumn<-"x"
ycolumn<-"y"
speedcolumn<-"speed"
vesselidcolumn<-"vesselid"
datetimecolumn<-"datetime"
inputTable<-"NAFO_VMS_SAMPLE_2011_corrected.csv"
occurrenceTimeRange<-8

yearsSelec<-"ALL"
seasonSelec<-"ALL"
sensitivity<-100
projectionExtension<-"POLYGON((-180 -90,-180 90,180 90,180 -90,-180 -90))"

inputFAOList<-"ASFIS_sp_2019_geom.txt"
minimumNumberOfRecords<-5

## load data
dVessel<-read.csv(inputTable,header=T,sep=",")
cat("column names:",names(dVessel),"\\n")
cat("selected parameters:\\n",
    "resolution:",res,"\\n",
    "xcolumn:",xcolumn,"\\n",
    "ycolumn:",ycolumn,"\\n",
    "speedcolumn:",speedcolumn,"\\n",
    "vesselidcolumn:",vesselidcolumn,"\\n",
    "datetimecolumn:",datetimecolumn,"\\n",
    "occurrenceTimeRange:",occurrenceTimeRange,"\\n"
    ,"\\n")

names(dVessel)[names(dVessel) == xcolumn] <- "x"
names(dVessel)[names(dVessel) == ycolumn] <- "y"
names(dVessel)[names(dVessel) == speedcolumn] <- "speed"
names(dVessel)[names(dVessel) == vesselidcolumn] <- "vesselid"
names(dVessel)[names(dVessel) == datetimecolumn] <- "datetime"

dVessel<-dVessel[,c("x","y","speed","vesselid","datetime")]

cat("altered column names:",names(dVessel),"\\n")

years<<-yearsSelec
seasons<<-seasonSelec

cat("User-selected years:",years, "\\n")
cat("User-selected months:",seasons, "\\n")
cat("User-selected min records:",minimumNumberOfRecords, "\\n")
cat("User-selected resolution:",res, "\\n")
cat("User-selected sensitivity:",sensitivity, "\\n")


cat("you are working with a resolution of:",res,"\\n")
extent_x<-180
extent_y<-90
n<-extent_y*2/res
m<-extent_x*2/res

cat("1. calculating the centroids of X and Y for the selected resolution\\n")
kmaxX<-ceiling((dVessel$x+extent_x)/res)
KminX<-kmaxX-1

dVessel$Xcentroid<--extent_x+(KminX*res)+res/2

kmaxY<-ceiling((dVessel$y+extent_y)/res)
KminY<-kmaxY-1

dVessel$Ycentroid<--extent_y+(KminY*res)+res/2

# exclude variables from dataset
listvardf <- names(dVessel) %in% c("x", "y")
dVesselNew <- dVessel[!listvardf]

cat("2. adding the bathymetry information for each X and Y \\n")





#fileBathy <- "http://thredds.d4science.org/thredds/dodsC/public/netcdf/depth_b2f62dfb-7b4b-428e-8601-4d1089308e14.nc"
fileBathy <- "depth_b2f62dfb-7b4b-428e-8601-4d1089308e14.nc"

if (!file.exists(fileBathy)){
  download.file("https://thredds.d4science.org/thredds/fileServer/public/netcdf/depth_b2f62dfb-7b4b-428e-8601-4d1089308e14.nc", "depth_b2f62dfb-7b4b-428e-8601-4d1089308e14.nc", "curl", quiet = FALSE, mode = "w",cacheOK = TRUE)
}  

cat("\\taccessing bathymetry remote file\\n")
dat.multi<-suppressWarnings(brick(fileBathy))
cat("\\tretrieving centroid columns\\n")
pt <-dVesselNew[,(ncol(dVesselNew)-1):ncol(dVesselNew)] #take the last two columns

cat("\\textracting bathymetry values for points\\n")
dVesselNew$depth<-raster::extract(dat.multi[[1]], pt)

rm(pt)

cat("3. classifying fishing activity by speed: Hauling (speed <= 2kn); 2) Fishing (2kn < speed <= 5kn); 3) Steaming (speed > 5kn)\\n")

for (j in 1:nrow(dVesselNew)){
  if(dVesselNew$speed[j]<=2){
    dVesselNew$fishing_activity[j]<-"Hauling"
    
  }else if (dVesselNew$speed[j]>2 && dVesselNew$speed[j]<=5){
    dVesselNew$fishing_activity[j]<-"Fishing"
    
  }else {
    dVesselNew$fishing_activity[j]<-"Steaming"
    
  }
}
cat("4. classifying fishing activity by speed and bathymetry\\n")

for (i in 1:nrow(dVesselNew)){
  
  if(dVesselNew$speed[i]<=2){
    dVesselNew$fishing_activity_bath[i]<-"Hauling"
    
  } else if (dVesselNew$speed[i]>2 && dVesselNew$speed[i]<=5){
    if (is.na(dVesselNew$depth[i]))
      dVesselNew$fishing_activity_bath[i]<-"Trawling"
    else if(dVesselNew$depth[i]>=500){
      dVesselNew$fishing_activity_bath[i]<-"Trawling"
    } else{
      dVesselNew$fishing_activity_bath[i]<-"Midwater trawling"}
  }else{
    dVesselNew$fishing_activity_bath[i]<-"Steaming"
    
  }
}

cat("5. calculating the fishing activity in hours: fahs\\n")
allvesselsids<-unique(dVesselNew$vesselid)
dVesselProcessed<-NA

for (vid in allvesselsids){
  
  dVesselNewOrd<-dVesselNew[dVesselNew$vesselid==vid,]
  dVesselNewOrd$datetime<-as.character(dVesselNewOrd$datetime)
  # print(vid)
  formattedDates<-c()
  for (datepoint in dVesselNewOrd$datetime){
    # cat("datepoint\\n")
    # print(datepoint)
    
    transformedDate<-as.POSIXlt(strptime(datepoint,format = c("%m/%d/%Y %I:%M:%S %p",
                                                              "%m/%d/%Y %H:%M:%S",
                                                              "%m/%d/%Y",
                                                              "%Y-%m-%dT%H:%M:%SZ",
                                                              "%Y-%m-%d %H:%M:%OS",
                                                              "%Y/%m/%d %H:%M:%OS",
                                                              "%Y-%m-%d %H:%M",
                                                              "%Y/%m/%d %H:%M",
                                                              "%Y-%m-%dT%H:%M:%OS",
                                                              "%Y/%m/%dT%H:%OS",
                                                              "%Y-%m-%dT%H:%M",
                                                              "%Y/%m/%dT%H:%M",
                                                              "%Y-%m-%d",
                                                              "%Y/%m/%d"), tz = "GMT"))
    
    transformedDate<-transformedDate[which(!is.na(transformedDate))][1]
    
    transformedDatef<-format(transformedDate,"%m/%d/%Y %H:%M:%S")
    # cat("transformed\\n")
    # print(transformedDatef)
    formattedDates<-c(formattedDates,transformedDatef)
  }
  
  cat("date transformation complete\\n")
  dVesselNewOrd$datetime<- as.POSIXlt(strptime(formattedDates,format = "%m/%d/%Y %H:%M:%S",tz="GMT"))
  
  dVesselNewOrd<-dVesselNewOrd[order(dVesselNewOrd$datetime),]
  
  dVesselNewOrd$tdiff <- unlist(tapply(dVesselNewOrd$datetime, INDEX = dVesselNewOrd$vesselid,
                                       FUN = function(x) c(0, `units<-`(diff(x), "hours"))))
  
  cat("6. assigning the fishing classes according to fahs for vessel",vid,"\\n")
  
  if (length(dVesselNewOrd$tdiff)==0 || length(which(is.na(dVesselNewOrd$tdiff)))>0)
    cat("!!! Wrong or inconsistent dates in the date column - Please use the  m/d/Y H:M:S format for time\\n")
  
  for (z in 1:nrow(dVesselNewOrd)){
    
    if (dVesselNewOrd$fishing_activity_bath[z]=="Midwater trawling" || dVesselNewOrd$fishing_activity_bath[z]=="Trawling"){
      
      if(dVesselNewOrd$tdiff[z]>4.0){
        dVesselNewOrd$fahs[z]<-0
      }else{
        dVesselNewOrd$fahs[z]<-dVesselNewOrd$tdiff[z]
      }
      
    }else{
      dVesselNewOrd$fahs[z]<-0
    }
  }
  
  if(is.na(dVesselProcessed)[1]){
    dVesselProcessed<-dVesselNewOrd
  } else{
    dVesselProcessed<-rbind(dVesselProcessed,dVesselNewOrd)
  }
  
}

# adding time 
time<-dVesselProcessed[,3]
seas<-dVesselProcessed[,3]
time<-as.Date(time)
seas<-as.Date(seas)
time<-as.numeric(format(time, "%Y"))
seas<-as.numeric(format(seas, "%m"))

dVesselProcessed$years<- time
dVesselProcessed$months<- seas

dataVessel <- dVesselProcessed

# exclude variables from dataset
exclude <- names(dataVessel) %in% c("speed", "vesselid","datetime","depth","fishing_activity","fishing_activity_bath","tdiff")
dataVessel <- dataVessel[!exclude]

##### source the other script
source("fishingGoogleDataSAI.R")


#########################################
cat("Map of vessels trajectories\\n")
# crop the raster to the area of interst
xminM<-max(min(dVessel$x)-10,-180)
xmaxM<-min(max(dVessel$x)+10,180)
yminM<-max(min(dVessel$y)-10,-90)
ymaxM<-min(max(dVessel$y)+10,90)

dat.multi2 <- crop(dat.multi[[1]], extent(xminM,xmaxM,yminM,ymaxM))
dat.multi2<-dat.multi2* (-1) # bathymetry to negative
dat.multi2_df <- as.data.frame(dat.multi2, xy = TRUE)# raster as data frame to deal better with ggplot2

## need to have as date format the column for the plot, plus x and y and vessel id
dVesselOrd <- dVessel[order(dVessel$vesselid),]
dVesselOrd$datetime<-as.character(dVesselOrd$datetime)

formattedDates<-c()
for (datepoint in dVesselOrd$datetime){
  transformedDate<-as.POSIXlt(strptime(datepoint,format = c("%m/%d/%Y %I:%M:%S %p",
                                                            "%m/%d/%Y",
                                                            "%Y-%m-%dT%H:%M:%SZ",
                                                            "%Y-%m-%d %H:%M:%OS",
                                                            "%Y/%m/%d %H:%M:%OS",
                                                            "%Y-%m-%d %H:%M",
                                                            "%Y/%m/%d %H:%M",
                                                            "%Y-%m-%dT%H:%M:%OS",
                                                            "%Y/%m/%dT%H:%OS",
                                                            "%Y-%m-%dT%H:%M",
                                                            "%Y/%m/%dT%H:%M",
                                                            "%Y-%m-%d",
                                                            "%Y/%m/%d"), tz = "GMT"))
  
  transformedDate<-transformedDate[which(!is.na(transformedDate))][1]
  
  transformedDatef<-format(transformedDate,"%m/%d/%Y %H:%M:%S")
  formattedDates<-c(formattedDates,transformedDatef)
}
dVesselOrd$datetime<- as.POSIXlt(strptime(formattedDates,format = "%m/%d/%Y %H:%M:%S",tz="GMT"))

dVesselOrd<-dVesselOrd[order(dVesselOrd$datetime),]


traj_vess<-ggplot() +
  geom_raster(data = dat.multi2_df, aes(x = x,y = y,fill=layer)) +
  geom_contour(data = dat.multi2_df, aes(x = x,y = y,z=layer),colour = "white", bins = 2)+
  geom_point(data = dVesselOrd, aes(x = x, y = y))+
  labs(x ="Longitude", y = "Latitude", fill= "Depth (m)",caption = "Vessels' trajectories. Colors represent different vessels")+
  theme_map()+
  theme(axis.title.x = element_text(size=40),#10),
        axis.title.y = element_text(size=40),
        axis.text=element_text(size=40),
        legend.title = element_text(size=40),
        legend.text = element_text(size=30),plot.caption=element_text(face = "italic",size=30))+
  scale_fill_continuous( guide = guide_colourbar())
traj_vess<-traj_vess+guides(fill = guide_colourbar(barwidth = 2, barheight = 10))
traj_vess2<-traj_vess+geom_line(data = dVesselOrd,aes(x = x, y = y,colour=as.factor(vesselid)))
traj_vess2<-traj_vess+geom_path(data= dVesselOrd, aes(x = x, y  = y, colour=as.factor(vesselid)))+  scale_colour_grey() 
traj_vess2<-traj_vess2+guides(color = FALSE)


##########################################

cat("saving results \\n")
destfileFileAreas<-"Exploited fishing areas.csv"
destfileFilesSpecies<-"Bar plots.csv"
destfileFilesSummaryTab<-"Summary table.html"
destfileFilesFullTab<-"Full table.csv"

destPlotTrajectory<-"Vessels trajectories.jpeg"
destPlotExpAreas<-"Exploited fishing areas.jpeg"
destPlotBarPlot<-"Bar plots.jpeg"

cat("Exporting images in: \\n")
outputImgTrj<-ggsave(destPlotTrajectory,traj_vess2,dpi=300,width=30,height=30, units="in")
outputImgExA<-ggsave(destPlotExpAreas,expl_loc,dpi=300,width=10,height=10)
outputImgBP<-ggsave(destPlotBarPlot,BarPlots,dpi=300,width=10,height=13)

## csv 
cat("Saving information on fishing areas and exploited species (.csv and .html format)")
captions<-paste("High medium and low fishing efforts at",res,"degree resolution")
FishingAreaExport<-dataTableAll_sf
st_geometry(FishingAreaExport) <- NULL

fileheader<-seq(1,ncol(dataTableAll_sf),by=1)
fileheader[1:length(fileheader)]<-","
fileheader[1]<-captions[1]

write(paste(fileheader,collapse=""),file = destfileFileAreas) 
outputFileExA<-write.table(FishingAreaExport[,-4],destfileFileAreas,append=TRUE,sep=",",col.names = c("Xcentroid","Ycentroid","Fishing hours","Fishing effort"),quote=TRUE,row.names = F,fileEncoding = "UTF-8")

captions<-"List of species fishing activity and pressures in the fishing areas."
fileheader<-seq(1,ncol(SppFA),by=1)
fileheader[1:length(fileheader)]<-","
fileheader[1]<-captions[1]

write(paste(fileheader,collapse=""),file = destfileFilesSpecies) 
outputFileSpp<-write.table(SppFA[,-6],destfileFilesSpecies,append=TRUE,col.names = c("Scientific name","Estimated Fishing Impact","Is Threatened","Fishing hours","N. of Occurrences"),sep=",",quote=TRUE,row.names = F,fileEncoding = "UTF-8")

captions<-"Full table."
fileheader<-seq(1,ncol(all_fao_obis_species3),by=1)
fileheader[1:length(fileheader)]<-","
fileheader[1]<-captions[1]

write(paste(fileheader,collapse=""),file = destfileFilesFullTab) 
outputFileFTab<-write.table(all_fao_obis_species3,destfileFilesFullTab,append=TRUE,col.names = c("Scientific name","Class","N. of Occurrences","Geometry","Fishing hours","Normalized fishing hours","Is a stock","Is High Fishing Area","Target Score","target Fishery","Estimated Fishing Impact Score","Estimated Fishing Impact","Is Threatened"),quote=TRUE, sep=",",row.names = F,fileEncoding = "UTF-8")

file.create(destfileFilesSummaryTab)
outputFileHtml<-DT::saveWidget(summary_html, normalizePath(destfileFilesSummaryTab), selfcontained = TRUE, libdir = "lib")


