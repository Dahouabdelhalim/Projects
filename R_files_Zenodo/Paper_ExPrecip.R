#Author: Joeri Reinders
#Purpose: Script which constructs daily precipitation totals in each Brazos, Trinity and Neches watershed (watershed shapefiles are constructed in a GIS). 
# Constructs a dataset with the ten biggest precipitation events that let to a >2year flood event and the ten biggest precipitatoin events that let to discharges not rising more than the 1.5year flood

rm(list=ls())
x11()

library(ncdf4)
library(ggplot2)
library(gridExtra)
library(raster)
library(sf)
library(tmap)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(R.matlab)
library(lubridate)
library(metR)
library(cowplot)
library(RColorBrewer)
library(dplR)
library(extRemes)
library(multiApply)
library(patchwork)
library(boot)

rm(list=ls())

drive_e <- 'D:/'
drive_l <- 'C:/'
region <- "Texas"

directory <- "Data/CPC/Precip/"  #change this into the folder you need it to be
file_name <- paste('All.',region,'.nc',sep="")
var <- "precip"

nc <- nc_open(paste(drive_e,directory,file_name, sep=""))
attributes(nc$var)

lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat", verbose = F)
time <- ncvar_get(nc, "time")

ndvi.array <- ncvar_get(nc, var) # store the data in a 3-dimensional array
dim(ndvi.array) 

fillvalue <- ncatt_get(nc, var, "_FillValue")
ndvi.array[ndvi.array == fillvalue$value] <- NA

dim(ndvi.array)
ndvi.array[ndvi.array < 0] <- 0
ndvi.array[is.na(ndvi.array)] <- 0
lon_ad <- (180-(lon-180))*-1
            
            
            #watershed_brazos <- st_read(paste(drive_e,directory,"Brazos_Richmond_Watersheduf.shp", sep=""))
            #watershed_trinity <- st_read(paste(drive_e,directory,"Trinity_Romayor_Watershed.shp", sep=""))
            #watershed_neches <- st_read(paste(drive_e,directory,"Neches_Evadale_Watershed.shp", sep=""))
            
            directory <- "QGIS/HoustonRiver_wsheddel/"
            file_brazos <- "Brazos_Richmond_Watershed_ras.tif"
            file_trinity <- "Trinity_Romayor_Watershed_ras.tif"
            file_neches <- "Neches_Evadale_Watershed_ras.tif"
            
            watershed_brazos <- raster(paste(drive_e,directory,file_brazos, sep=""))
            watershed_trinity <- raster(paste(drive_e,directory,file_trinity, sep=""))
            watershed_neches <- raster(paste(drive_e,directory,file_neches, sep=""))
            
            watershed_brazos_pol <- rasterToPolygons(watershed_brazos)
            watershed_trinity_pol <- rasterToPolygons(watershed_trinity)
            watershed_neches_pol <- rasterToPolygons(watershed_neches)
            
            test <- ndvi.array[,,1]
            r <- raster(t(test), xmn=lon_ad[1], xmx=lon_ad[80], ymn=lat[1], ymx=lat[80],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
            r <- flip(r,direction="y")
            plot(r,col = rev(terrain.colors(50)))
            
            watershed_brazos_masked <- mask(r,watershed_brazos_pol)
            watershed_trinity_masked <- mask(r,watershed_trinity_pol)
            watershed_neches_masked <- mask(r,watershed_neches_pol)
            
            # Loop - extracts daily precipitation levels from the three watersheds and adds them to a vector! 
            
            precip_brazos <- 1:length(time)
            precip_trinity <- 1:length(time)
            precip_neches <- 1:length(time)
            
            for(day in 1:length(time)){
              test <- ndvi.array[,,day]
              r <- raster(t(test), xmn=lon_ad[1], xmx=lon_ad[80], ymn=lat[1], ymx=lat[80],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
              r <- flip(r,direction="y")
             
              watershed_brazos_masked <- mask(r,watershed_brazos_pol)
              watershed_trinity_masked <- mask(r,watershed_trinity_pol)
              watershed_neches_masked <- mask(r,watershed_neches_pol)
              
              sum_brazos <- raster::as.data.frame(watershed_brazos_masked, xy = TRUE)
              precip_brazos[day] <- sum(sum_brazos[,3],na.rm=T)
              
              sum_trinity <- raster::as.data.frame(watershed_trinity_masked, xy = TRUE)
              precip_trinity[day] <- sum(sum_trinity[,3],na.rm=T)
              
              sum_neches <- raster::as.data.frame(watershed_neches_masked, xy = TRUE)
              precip_neches[day] <- sum(sum_neches[,3],na.rm=T)
            }
            
            Precip_Rivers <- data.frame(Neches = precip_neches,Trinity=precip_trinity,Brazos=precip_brazos) #dataframe with the above three daily precip records
            
            
            
            directory <- "Users/joeri/OneDrive/Bureaublad/Northeastern University/Research Chapters/Chapter 4/"
            setwd(paste(drive_l,directory,"Discharg_data/",sep=""))
            getwd()
            
            #write.csv(Precip_Rivers,"Precip_RiversWatersheds.csv") #precipitation csv

            
########################################            
## Extract floods with precip record! ##
########################################            
              
  directory <- "Users/joeri/OneDrive/Bureaublad/Northeastern University/Research Chapters/Chapter 4/"
  precipdata <- read.csv(paste(drive_l,directory,"Discharg_data/Precip_RiversWatersheds.csv",sep=""))
  precipdata$Date <- seq(as.Date("1948-01-01"), as.Date("2020-12-31"), by="days")
  precipdata$Month <- month(precipdata$Date)
  precipdata$Day <- day(precipdata$Date)
  precipdata$Year <- year(precipdata$Date)
  precipdata$Date  <- as.Date(precipdata$Date , "%m/%d/%y")
  
  spring_precip_neches <- subset(precipdata, Month==4|Month==5|Month==6,select=c(Date, Day, Month, Year, Neches))
  spring_precip_trinity <- subset(precipdata, Month==4|Month==5|Month==6,select=c(Date, Day, Month, Year, Trinity))
  spring_precip_brazos <- subset(precipdata, Month==4|Month==5|Month==6,select=c(Date, Day, Month, Year, Brazos))
  
  SS_precip_Neches <- spring_precip_neches[order(spring_precip_neches$Neches, decreasing = T),]
  SS_precip_Trinity <- spring_precip_trinity[order(spring_precip_trinity$Trinity, decreasing = T),]
  SS_precip_Brazos <- spring_precip_brazos[order(spring_precip_brazos$Brazos, decreasing = T),]
  
  #flood data
    drive <- 'C:/'
    directory <- "Users/joeri/OneDrive/Bureaublad/Northeastern University/Research Chapters/Chapter 4/"
    file_neches <- paste(drive,directory,"Discharg_data/Daily_Neches@Evadale.csv",sep="")
    file_trinity <- paste(drive,directory,"Discharg_data/Daily_Trinity@Romayor.csv",sep="")
    file_brazos <- paste(drive,directory,"Discharg_data/Daily_Brazos@Richmond.csv",sep="")
    orgNeches <- read.csv(file_neches)
    orgTrinity <- read.csv(file_trinity)
    orgBrazos <- read.csv(file_brazos)
    
    spring_neches <- subset(orgNeches, Month==4|Month==5|Month==6)
    spring_trinity <- subset(orgTrinity, Month==4|Month==5|Month==6)
    spring_brazos <- subset(orgBrazos, Month==4|Month==5|Month==6)
  
    Peak_neches <- as.numeric(tapply(spring_neches$Discharge,spring_neches$Year,max))
    Peak_trinity <- as.numeric(tapply(spring_trinity$Discharge,spring_trinity$Year,max))
    Peak_brazos <- as.numeric(tapply(spring_brazos$Discharge,spring_brazos$Year,max))
    
    Neches2yF <- as.numeric(return.level.fevd(fevd(Peak_neches)))[1]
    Trinity2yF <- as.numeric(return.level.fevd(fevd(Peak_trinity)))[1]
    Brazos2yF <- as.numeric(return.level.fevd(fevd(Peak_brazos)))[1]
    
    Neches1.5yF <- as.numeric(return.level(fevd(Peak_neches), return.period =c(1.5)))[1]
    Trinity1.5yF <- as.numeric(return.level(fevd(Peak_trinity), return.period =c(1.5)))[1]
    Brazos1.5yF <- as.numeric(return.level(fevd(Peak_brazos), return.period =c(1.5)))[1]
    
  
  ###################
  ## Top 10 Neches ##  
  ###################  
    
  top10_neches_precipflood <- rep(NA, 10)
  top10_neches_precipnoflood <- rep(NA, 10)
  Index <- 1
  
  while(anyNA(c(top10_neches_precipnoflood,top10_neches_precipflood))==T){
    
    print(Index)
    d1 <- as.Date(((SS_precip_Neches$Date[Index])), format="%m/%d/%Y")
    floodin <- which(as.Date(((spring_neches$Date)), format="%m/%d/%Y")==d1)
    discharge_level <- spring_neches$Discharge[floodin:(floodin+7)]
    
    #higher
    if(any(discharge_level>Neches2yF)){
      if(sum(is.na(top10_neches_precipflood))==0){Index <- Index+1
      print("check")}
      else{
        if(sum(is.na(top10_neches_precipflood))==10){top10_neches_precipflood[1] <- as.character(d1)}
        
        vol <- abs(sum(is.na(top10_neches_precipflood))-11)
        lengte <- vol - 1 
        dif <- rep(NA,lengte)
        
        for(index in 1:lengte){
          d2 <- (as.Date((as.character(top10_neches_precipflood[index]))))
          dif[index] <- abs(as.numeric(difftime(d2,d1 , units="day")))
          }
          
          if(any(dif<16)){Index <- Index+1}
          else{ top10_neches_precipflood[vol] <- as.character(d1)
                Index <- Index+1}
          }
    }
    
    #lower
    if(length(which((discharge_level<(Neches1.5yF))==T))==8){
      
      if(sum(is.na(top10_neches_precipnoflood))==0){Index <- Index+1}
      else{
        if(sum(is.na(top10_neches_precipnoflood))==10){top10_neches_precipnoflood[1] <- as.character(d1)}
        
        vol <- abs(sum(is.na(top10_neches_precipnoflood))-11)
        lengte <- vol - 1 
        dif <- rep(NA,lengte)
        
        for(index in 1:lengte){
          d2 <- (as.Date((as.character(top10_neches_precipnoflood[index]))))
          dif[index] <- abs(as.numeric(difftime(d2,d1 , units="day")))
        }
        
        if(any(dif<16)){Index <- Index+1}
        else{ top10_neches_precipnoflood[vol] <- as.character(d1)
              Index <- Index+1}
        }
    }
    
    #inbetween
    if((length(which((discharge_level<(Neches1.5yF))==T))!=8)&& (discharge_level<Neches2yF)){Index <- Index+1
    print("test")}
      
    }
    
  
  
  
  
  ####################
  ## Top 10 Trinity ##  
  ####################  
  
  top10_trinity_precipflood <- rep(NA, 10)
  top10_trinity_precipnoflood <- rep(NA, 10)
  Index <- 1
  
  while(anyNA(c(top10_trinity_precipnoflood,top10_trinity_precipflood))==T){
    
    print(Index)
    d1 <- as.Date(((SS_precip_Trinity$Date[Index])), format="%m/%d/%Y")
    floodin <- which(as.Date(((spring_trinity$Date)), format="%m/%d/%Y")==d1)
    discharge_level <- spring_trinity$Discharge[floodin:(floodin+7)]
    
    #higher
    if(any(discharge_level>Trinity2yF)){
      if(sum(is.na(top10_trinity_precipflood))==0){Index <- Index+1
      print("check")}
      else{
        if(sum(is.na(top10_trinity_precipflood))==10){top10_trinity_precipflood[1] <- as.character(d1)}
        
        vol <- abs(sum(is.na(top10_trinity_precipflood))-11)
        lengte <- vol - 1 
        dif <- rep(NA,lengte)
        
        for(index in 1:lengte){
          d2 <- (as.Date((as.character(top10_trinity_precipflood[index]))))
          dif[index] <- abs(as.numeric(difftime(d2,d1 , units="day")))
        }
        
        if(any(dif<16)){Index <- Index+1}
        else{ top10_trinity_precipflood[vol] <- as.character(d1)
        Index <- Index+1}
      }
    }
    
    #lower
    if(length(which((discharge_level<(Trinity1.5yF))==T))==8){
      
      if(sum(is.na(top10_trinity_precipnoflood))==0){Index <- Index+1}
      else{
        if(sum(is.na(top10_trinity_precipnoflood))==10){top10_trinity_precipnoflood[1] <- as.character(d1)}
        
        vol <- abs(sum(is.na(top10_trinity_precipnoflood))-11)
        lengte <- vol - 1 
        dif <- rep(NA,lengte)
        
        for(index in 1:lengte){
          d2 <- (as.Date((as.character(top10_trinity_precipnoflood[index]))))
          dif[index] <- abs(as.numeric(difftime(d2,d1 , units="day")))
        }
        
        if(any(dif<16)){Index <- Index+1}
        else{ top10_trinity_precipnoflood[vol] <- as.character(d1)
        Index <- Index+1}
      }
    }
    
    #inbetween
    if((length(which((discharge_level<(Trinity1.5yF))==T))!=8)&& (discharge_level<Trinity2yF)){Index <- Index+1
    print("test")}
    
  }
  
  
  
  
  
  ###################
  ## Top 10 Brazos ##  
  ###################  
  
  top10_brazos_precipflood <- rep(NA, 10)
  top10_brazos_precipnoflood <- rep(NA, 10)
  Index <- 1
  
  while(anyNA(c(top10_brazos_precipnoflood,top10_brazos_precipflood))==T){
    
    print(Index)
    d1 <- as.Date(((SS_precip_Brazos$Date[Index])), format="%m/%d/%Y")
    floodin <- which(as.Date(((spring_brazos$Date)), format="%m/%d/%Y")==d1)
    discharge_level <- spring_brazos$Discharge[floodin:(floodin+7)]
    
    #higher
    if(any(discharge_level>Brazos2yF)){
      if(sum(is.na(top10_brazos_precipflood))==0){Index <- Index+1
      print("check")}
      else{
        if(sum(is.na(top10_brazos_precipflood))==10){top10_brazos_precipflood[1] <- as.character(d1)}
        
        vol <- abs(sum(is.na(top10_brazos_precipflood))-11)
        lengte <- vol - 1 
        dif <- rep(NA,lengte)
        
        for(index in 1:lengte){
          d2 <- (as.Date((as.character(top10_brazos_precipflood[index]))))
          dif[index] <- abs(as.numeric(difftime(d2,d1 , units="day")))
        }
        
        if(any(dif<16)){Index <- Index+1}
        else{ top10_brazos_precipflood[vol] <- as.character(d1)
        Index <- Index+1}
      }
    }
    
    #lower
    if(length(which((discharge_level<(Brazos1.5yF))==T))==8){
      
      if(sum(is.na(top10_brazos_precipnoflood))==0){Index <- Index+1}
      else{
        if(sum(is.na(top10_brazos_precipnoflood))==10){top10_brazos_precipnoflood[1] <- as.character(d1)}
        
        vol <- abs(sum(is.na(top10_brazos_precipnoflood))-11)
        lengte <- vol - 1 
        dif <- rep(NA,lengte)
        
        for(index in 1:lengte){
          d2 <- (as.Date((as.character(top10_brazos_precipnoflood[index]))))
          dif[index] <- abs(as.numeric(difftime(d2,d1 , units="day")))
        }
        
        if(any(dif<16)){Index <- Index+1}
        else{ top10_brazos_precipnoflood[vol] <- as.character(d1)
        Index <- Index+1}
      }
    }
    
    #inbetween
    if((length(which((discharge_level<(Brazos1.5yF))==T))!=8)&& (discharge_level<Brazos2yF)){Index <- Index+1
    print("test")}
    
  }
  
  
  
  
  
########  
  

csvFloodPrecip <- data.frame( NechesFlood = top10_neches_precipflood, NechesNoFlood = top10_neches_precipnoflood,
                              TrinityFlood = top10_trinity_precipflood, TrinityNoFlood = top10_trinity_precipnoflood,
                              BrazosFlood = top10_brazos_precipflood, BrazosNoFlood = top10_brazos_precipnoflood)

directory <- "Users/joeri/OneDrive/Bureaublad/Northeastern University/Research Chapters/Chapter 4/"
setwd(paste(drive_l,directory,"Discharg_data/",sep=""))
getwd()

write.csv(csvFloodPrecip,"precip_peak_floods.csv")


plot_list = list()
discharge_files <- list(spring_neches,spring_trinity,spring_brazos) 
plot_title <- c("Neches", "Trinity", "Brazos")
letters1 <- c("(f)  ","(e)  ","(d)  ")
letters2 <- c("(c)  ","(b)  ","(a)  ")
line <- list(Neches2yF,Trinity2yF,Brazos2yF)
index <- c(1,3,5)
test <- data.frame(spring_precip_neches$Neches,spring_precip_trinity$Trinity,spring_precip_brazos$Brazos)

for(plot in 1:3){
  river <- index[plot]
  
  simpel_fl <- matrix(NA,nrow=10,ncol=2)
  for(l in 1:10){
    floodin1 <- which(as.Date((as.character(discharge_files[[plot]]$Date)), format="%m/%d/%Y")==as.character(csvFloodPrecip[l,river]))
    floodin2 <- which(as.Date((as.character(discharge_files[[plot]]$Date)), format="%m/%d/%Y")==as.character(csvFloodPrecip[l,river+1]))
    discharge_level1 <- max(discharge_files[[plot]]$Discharge[floodin1:(floodin1+7)])
    discharge_level2 <- max(discharge_files[[plot]]$Discharge[floodin2:(floodin2+7)])
    simpel_fl[l,1] <- discharge_level1
    simpel_fl[l,2] <- discharge_level2
    
  }
  
  simpel_pr <- matrix(NA,nrow=10,ncol=2)
  for(f in 1:10){
    floodin1 <- which(as.character(spring_precip_brazos$Date)==as.character(csvFloodPrecip[f,river]))
    floodin2 <- which(as.character(spring_precip_brazos$Date)==as.character(csvFloodPrecip[f,river+1]))
    simpel_pr[f,1] <- test[floodin1, plot]
    simpel_pr[f,2] <- test[floodin2, plot]
  }

  level <- c(simpel_fl[,2], simpel_fl[,1])
  name <- c(rep("No Flooding",10),rep("Flooding",10))
  histogram_dataframe_fl <- data.frame("level"=level, "name"=name)
  
  level2 <- c(simpel_pr[,2], simpel_pr[,1])
  name2 <- c(rep("noflood",10),rep("flood",10))
  histogram_dataframe_pr <- data.frame("level2"=level2, "name2"=name2)
  
  j <- ggplot(histogram_dataframe_fl, aes(x=level, fill=name))+
        geom_histogram(color="black", alpha=0.5, position="identity", bins=20)+
        geom_vline(xintercept=line[[plot]],color="red", size=1)+
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal()+ 
        ylim(0, 4)+
        labs(x="Discharge (ft3/s)", y="Count",fill=paste("", sep=""), title=paste(letters1[plot],"  ",plot_title[plot], sep=""))+
        theme(plot.title = element_text(hjust = 0, face="bold", size=14))+
        theme(axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 13),axis.title.y = element_text(size = 14),axis.text.y = element_text(size = 13))+
        theme(legend.position="bottom")+
        theme(legend.key.size = unit(0.8, 'cm'))+
        theme(legend.text = element_text(size=14))
  
  k <- ggplot(histogram_dataframe_pr, aes(x=level2, fill=name2))+
        geom_histogram(color="black", alpha=0.5, position="identity", bins=20)+
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_minimal()+ 
        ylim(0,8)+
        labs(x="Precipitation (mm)", y="Count",fill=paste("", sep=""), title=paste(letters2[plot],"  ",plot_title[plot], sep=""))+
        theme(plot.title = element_text(hjust = 0, face="bold", size=14))+
        theme(axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 13),axis.title.y = element_text(size = 14),axis.text.y = element_text(size = 13))+
        theme(legend.position="none")+
        theme(legend.key.size = unit(0, 'cm'))+
        theme(legend.text = element_text(size=0))

  plot_list[[plot]] <- j
  plot_list[[plot+3]] <- k
}


combineda <- plot_list[[3]] + plot_list[[2]] + plot_list[[1]] 
combinedb <- plot_list[[6]] + plot_list[[5]] + plot_list[[4]] & theme(legend.position = "none")
combinedb / combineda + plot_layout(guides = "collect") & theme(legend.position = "bottom")


