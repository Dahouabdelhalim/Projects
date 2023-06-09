#Author: Joeri Reinders
#Purpose: Script which creates composite maps of meridional winds reanalysis data for days of flooding [given by an csv document] - figure 4 in paper.
#Data: contact reinders.j@northeastern.edu for orginal data files

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
library(patchwork) 
library(viridis)

rm(list=ls())

drive_e <- 'D:/'
drive_l <- 'C:/'
level <- 850
recordlength <- 26664

mean_fun <- function(data){
  resultaat <- mean(data,na.rm=T)
  return(resultaat)
} #
sd_fun <- function(data){
  resultaat <- sd(data,na.rm=T)
  return(resultaat)
} #
houston <- data.frame(x=c(-95.391), y=(c=29.5787))
plot_list = list()
myPallette1 <- c(rev(brewer.pal(11, "RdBu")))
myPallette2 <- c((brewer.pal(9, "GnBu"))[2:7],brewer.pal(9,"BuPu")[5:7])

directory <- "Users/joeri/OneDrive/Bureaublad/Northeastern University/Research Chapters/Chapter 4/"
file_10 <- paste(drive_l,directory,"Discharg_data/Top10Spring.csv",sep="")
Spring_Top10 <- read.csv(file_10)


##################
## Wind Vectors ##     Produces a dataframe wind_data which contains wind directions and speed
##################

w.level <- 850 

directory_uwind <- "Data/Reanalysis/NCEP-NCAR/Uwinds/daily_data/"  #change this into the folder you need it to be
file_name_uwind <- paste('All.',w.level,'.nc',sep="")
var_uwind <- "uwnd"

directory_vwind <- "Data/Reanalysis/NCEP-NCAR/Vwinds/daily_data/"  #change this into the folder you need it to be
file_name_vwind <- paste('All.',w.level,'.nc',sep="")
var_vwind <- "vwnd"

nc_uwind <- nc_open(paste(drive_e,directory_uwind,file_name_uwind, sep=""))
nc_vwind <- nc_open(paste(drive_e,directory_vwind,file_name_vwind, sep=""))

ulon <- ncvar_get(nc_uwind, "lon")
ulat <- ncvar_get(nc_uwind, "lat", verbose = F)
utime <- ncvar_get(nc_uwind, "time")

vlon <- ncvar_get(nc_vwind, "lon")
vlat <- ncvar_get(nc_vwind, "lat", verbose = F)
vtime <- ncvar_get(nc_vwind, "time")

u.ndvi.array <- ncvar_get(nc_uwind, var_uwind) # store the data in a 3-dimensional array
fillvalue <- ncatt_get(nc_uwind, var_uwind, "_FillValue")
u.ndvi.array[u.ndvi.array == fillvalue$value] <- NA
v.ndvi.array <- ncvar_get(nc_vwind, var_vwind) # store the data in a 3-dimensional array
fillvalue <- ncatt_get(nc_vwind, var_vwind, "_FillValue")
v.ndvi.array[v.ndvi.array == fillvalue$value] <- NA


naam_plot4A <- c("(a)   Normal Wind","(b)   Flood Wind")
for(plot4A in 1:2){
   
  if(plot4A==1){
    dates <- seq(as.Date("1948-01-01"), as.Date("2020-12-31"), by = "days") 
    months <- format(dates, "%m")
    dates_lente <- which((months=="04")|(months=="05")|(months=="06"))
    
    general_map.index <- sample(dates_lente,1000,replace=TRUE)
    }
  
  if(plot4A==2){
    date_plot <- c(as.character(Spring_Top10[,5]),as.character(Spring_Top10[,6]),as.character(Spring_Top10[,7])) #make sure to put this as character
    #combines the dates for all three rivers (5,6,7 are columns in the csv)
    vector_dates <- seq( as.Date("1948-01-01"), by=1, len=recordlength) # a vector with dates to determine the posistion of a certain day in the nc file
    dates_of_interest <- date_plot # a vector with dates in which we are interested
    general_map.index <- 1:length(dates_of_interest) #string for index numbers of the picked dates
    for(i in 1:length(general_map.index)){
      general_map.index[i] <- which(vector_dates==dates_of_interest[i])
    } #retrieves the indexes for the dates for which we are interested
    }
  
  u.general_map.data <- u.ndvi.array[,,general_map.index]
  u.general_map.slice <- apply(u.general_map.data,1:2, mean_fun)
  v.general_map.data <- v.ndvi.array[,,general_map.index]
  v.general_map.slice <- apply(v.general_map.data,1:2, mean_fun)

  u.general_map.r <- raster(t(u.general_map.slice), xmn=ulon[1], xmx=ulon[144], ymn=ulat[1], ymx=ulat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  v.general_map.r <- raster(t(v.general_map.slice), xmn=vlon[1], xmx=vlon[144], ymn=vlat[1], ymx=vlat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  df.v.gm <- raster::as.data.frame(v.general_map.r, xy = TRUE)
  df.u.gm <- raster::as.data.frame(u.general_map.r, xy = TRUE)
  wind_data.gm <- data.frame(x=df.u.gm[,1],y=df.u.gm[,2],u=df.u.gm[,3],v=df.v.gm[,3])
  wind_data.gm$velocity <- sqrt((wind_data.gm$u^2)+(wind_data.gm$v^2))
  wind_data.gm$selectie <- rep(c(0,1,2,3),(length(wind_data.gm$velocity)/4))
  wind_data.gm <- subset(wind_data.gm,selectie==1 )
  
  cells <- length(wind_data.gm$x)
  for(i in 1:cells){ 
    if(wind_data.gm$x[i]<180){wind_data.gm$x[i]<-wind_data.gm$x[i]-((360/144)/2)}
    if(wind_data.gm$x[i]>180){wind_data.gm$x[i]<-wind_data.gm$x[i]-360+((360/144)/2)}
  }

    j <- ggplot()+
          metR::geom_contour_fill(data=wind_data.gm, aes(x=x, y=y, z=velocity),binwidth = 1)+
          geom_segment(data=wind_data.gm, aes(x=x, xend = x+u, y = y, yend = y+v), arrow = arrow(length = unit(0.15, "cm")))+
          geom_point(data=houston, aes(x=x, y=y), color="red")+
          borders('world', xlim=range(wind_data.gm$x), ylim=range(wind_data.gm$y), colour='black', size=.4)+
          borders('state', colour="black", size=0.1)+
          theme_minimal()+ 
          theme(panel.ontop=TRUE, panel.background=element_blank()) +
          theme(axis.title.x = element_text(size = 13),axis.text.x = element_text(size = 12),axis.title.y = element_text(size = 13),axis.text.y = element_text(size = 12))+
          scale_fill_gradientn(colors =myPallette2, limits=c(0,16), breaks=c(0,2,4,6,8,10,12,14,16),
                               guide = guide_colorsteps(frame.colour = "black", ticks.colour = "black", barwidth=20, barheight=0.8, title.position="top",title="Wind speed (m/s)", 
                                                        title.hjust=0.5, title.theme=element_text(size=13),label.theme=element_text(size=12)))+
          coord_quickmap(xlim=c(-120,-40), ylim=c(10,60))+
          labs(x="Longitude", y="Latitude",fill=paste("", sep=""))+
          theme(legend.position = "bottom")+
          theme(plot.title = element_text(hjust = 0, face="bold", size=13))+
          theme(legend.margin=margin(20, 0, 0, 0))+
          labs(title=paste(naam_plot4A[plot4A]))
  plot_list[[plot4A+6]] <- j
}





    ###################
    ## Daily Anomoly ## DONT RUN THIS IF YOU WANT ABSOLUTE VALUES
    ###################
    
    #U wind
    u.ndvi.normal.m <- array(NA, dim=c(144,73,366))
    u.ndvi.normal.sd <- array(NA, dim=c(144,73,366))
    
    dates <- seq(as.Date("1948-01-01"), as.Date("2020-12-31"), by = "days") 
    index.max <- which(dates=="2020-12-31")
    date_dates <- format(dates, "%m-%d")
    
    dates_anom <- seq(as.Date("1981-01-01"), as.Date("2010-12-31"), by = "days")
    index_anom.min <- which((dates==as.Date("1981-01-01")))
    date_dates_anom <- format(dates_anom, "%m-%d")
    
    dates_ref <- format((seq(as.Date("1948-01-01"), as.Date("1948-12-31"), by = "days")),"%m-%d") 
    
    for(day in 1:366){
      date_index <- (which(dates_ref[day]==date_dates_anom)+index_anom.min)
      ndvi.normal <- u.ndvi.array[,,date_index]
      u.ndvi.normal.m[,,day] <- apply(ndvi.normal,1:2, mean_fun)
      u.ndvi.normal.sd[,,day] <- apply(ndvi.normal,1:2, sd_fun)
    }
    
    for(days in 1:366){
      date_index <- which(dates_ref[days]==date_dates)
      u.ndvi.array[,,date_index] <- sweep(u.ndvi.array[,,date_index],1:2,u.ndvi.normal.m[,,days])/replicate(length(date_index), u.ndvi.normal.sd[,,days], simplify="array")
    }
    
    
    #V wind
    v.ndvi.normal.m <- array(NA, dim=c(144,73,366))
    v.ndvi.normal.sd <- array(NA, dim=c(144,73,366))
    
    dates <- seq(as.Date("1948-01-01"), as.Date("2020-12-31"), by = "days") 
    index.max <- which(dates=="2020-12-31")
    date_dates <- format(dates, "%m-%d")
    
    dates_anom <- seq(as.Date("1981-01-01"), as.Date("2010-12-31"), by = "days")
    index_anom.min <- which((dates==as.Date("1981-01-01")))
    date_dates_anom <- format(dates_anom, "%m-%d")
    
    dates_ref <- format((seq(as.Date("1948-01-01"), as.Date("1948-12-31"), by = "days")),"%m-%d") 
    
    for(day in 1:366){
      date_index <- (which(dates_ref[day]==date_dates_anom)+index_anom.min)
      ndvi.normal <- v.ndvi.array[,,date_index]
      v.ndvi.normal.m[,,day] <- apply(ndvi.normal,1:2, mean_fun)
      v.ndvi.normal.sd[,,day] <- apply(ndvi.normal,1:2, sd_fun)
    }
    
    for(days in 1:366){
      date_index <- which(dates_ref[days]==date_dates)
      v.ndvi.array[,,date_index] <- sweep(v.ndvi.array[,,date_index],1:2,v.ndvi.normal.m[,,days])/replicate(length(date_index), v.ndvi.normal.sd[,,days], simplify="array")
    }

    
    
    
    
mybreak <- seq(-1,1,by=0.2)
tijd_naam <- c("Day", "Week", "Two Weeks")

bootstrap_function <- function(B,n,data,ref,delay,alpha){
  
  dates <- seq(as.Date("1948-01-01"), as.Date("2020-12-31"), by = "days") 
  months <- format(dates, "%m")
  dates_lente <- which((months=="04")|(months=="05")|(months=="06"))
  
  results.bootstrap <- array(data=NA, dim=c(144,73,B))
  
  for(bt in 1:B){
    dates_sampled <- sample(dates_lente,n,replace=TRUE)
    
    if(delay>1){
      c <- 1
      index <- rep(NA,(delay*length(dates_sampled)))
      
      for(a in 1:length(dates_sampled)){
        for(b in 1:delay){
          index[c] <- dates_sampled[a]-b+1
          c <- c+1
        }
      }
      
    }else{index <- dates_sampled}
    
    ndvi.slices.bootstrap <- data[,,index]
    results.bootstrap[,,bt] <- apply(ndvi.slices.bootstrap,1:2, mean_fun)
  }
  
  
  results.bootstrap2 <- array(data=NA, dim=c(144,73))
  for(lon.b in 1:144){
    for(lat.b in 1:73){
      bina <- sum(abs(results.bootstrap[lon.b,lat.b,])>abs(ref[lon.b,lat.b]))/B  
      if(bina<alpha){results.bootstrap2[lon.b,lat.b] <- 1}
      else{results.bootstrap2[lon.b,lat.b] <- 0}
    }
  }
  
  return(results.bootstrap2)
}
mean_fun <- function(data){
  resultaat <- mean(data,na.rm=T)
  return(resultaat)
} #
alpha <- 0.01
B <- 10000

letter1 <- c("(c)", "(d)", "(e)")
letter2 <- c("(f)", "(g)", "(h)")

  for(tijd in 1:3){
      
    date_plot <- c(as.character(Spring_Top10[,5]),as.character(Spring_Top10[,6]),as.character(Spring_Top10[,7])) #make sure to put this as character
    #combines the dates for all three rivers (5,6,7 are columns in the csv)
      
    vector_dates <- seq( as.Date("1948-01-01"), by=1, len=recordlength) # a vector with dates to determine the posistion of a certain day in the nc file
    dates_of_interest <- date_plot # a vector with dates in which we are interested
      
      if(tijd==1){
        #Day before flood:
        index_day <- 1:length(dates_of_interest) #string for index numbers of the picked dates
        for(i in 1:length(index_day)){
          index_day[i] <- which(vector_dates==dates_of_interest[i])
        } #retrieves the indexes for the dates for which we are interested

        u.ndvi.slices <- u.ndvi.array[,,index_day]
        u.ndvi.slice <- apply(u.ndvi.slices,1:2, mean_fun)
        v.ndvi.slices <- v.ndvi.array[,,index_day]
        v.ndvi.slice <- apply(v.ndvi.slices,1:2, mean_fun)
        
        u.results.bootstrap2 <- bootstrap_function(B,30,u.ndvi.array,u.ndvi.slice,1,alpha)
        v.results.bootstrap2 <- bootstrap_function(B,30,v.ndvi.array,v.ndvi.slice,1,alpha)
        
        
        #results.bootstrap2 <- bootstrap_function(1000,30,ndvi.array,ndvi.slice,1,alpha)
        
      }
      
      if(tijd==2){
        #Week before flood:
        c <- 1
        index_week <- rep(NA,(7*length(dates_of_interest)))
        
        for(a in 1:length(dates_of_interest)){
          for(b in 1:7){
            date <- (as.Date(dates_of_interest[a])-b+1)
            index_week[c] <- which(vector_dates==date)
            c <- c+1
          }
        }

        u.ndvi.slices <- u.ndvi.array[,,index_week]
        u.ndvi.slice <- apply(u.ndvi.slices,1:2, mean_fun)
        v.ndvi.slices <- v.ndvi.array[,,index_week]
        v.ndvi.slice <- apply(v.ndvi.slices,1:2, mean_fun)
        
        u.results.bootstrap2 <- bootstrap_function(B,30,u.ndvi.array,u.ndvi.slice,1,alpha)
        v.results.bootstrap2 <- bootstrap_function(B,30,v.ndvi.array,v.ndvi.slice,1,alpha)
        
      
      }
      
      if(tijd==3){
        #Three weeks before flood!:
        c <- 1
        index_3week <- rep(NA,(21*length(dates_of_interest)))
        for(a in 1:length(dates_of_interest)){
          for(b in 1:21){
            date <- (as.Date(dates_of_interest[a])-b+1)
            index_3week[c] <- which(vector_dates==date)
            c <- c+1
          }
        }
        u.ndvi.slices <- u.ndvi.array[,,index_3week]
        u.ndvi.slice <- apply(u.ndvi.slices,1:2, mean_fun)
        v.ndvi.slices <- v.ndvi.array[,,index_3week]
        v.ndvi.slice <- apply(v.ndvi.slices,1:2, mean_fun)
        
        u.results.bootstrap2 <- bootstrap_function(B,30,u.ndvi.array,u.ndvi.slice,1,alpha)
        v.results.bootstrap2 <- bootstrap_function(B,30,v.ndvi.array,v.ndvi.slice,1,alpha)
        

      }
      
      #Prepares composistion map for ggplot
    
      r.u.wind <- raster(t(u.ndvi.slice), xmn=vlon[1], xmx=vlon[144], ymn=vlat[1], ymx=vlat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      r.v.wind <- raster(t(v.ndvi.slice), xmn=vlon[1], xmx=vlon[144], ymn=vlat[1], ymx=vlat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      
      # first the raster is molded into a dataframe which can be plotted in ggplot and the lot is adjusted for the world map
      dfv <- raster::as.data.frame(r.v.wind, xy = TRUE)
      dfu <- raster::as.data.frame(r.u.wind, xy = TRUE)
      wind_data <- data.frame(x=dfu[,1],y=dfu[,2],u=dfu[,3],v=dfv[,3])
      wind_data$velocity <- sqrt((wind_data$u^2)+(wind_data$v^2))
      
      
      cells <- length(wind_data$x)
      
      for(i in 1:cells){ 
        if(wind_data$x[i]<180){wind_data$x[i]<-wind_data$x[i]-((360/144)/2)}
        if(wind_data$x[i]>180){wind_data$x[i]<-wind_data$x[i]-360+((360/144)/2)}
      }
      
      #Prepares the bootstrap results in a raster and dataframe format for ggplot
      u.r.b <- raster(t(u.results.bootstrap2), xmn=vlon[1], xmx=vlon[144], ymn=vlat[1], ymx=vlat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      u.flood.b <- raster::as.data.frame(u.r.b, xy = TRUE)
      cells <- length(u.flood.b$x)
      for(i in 1:cells){ 
        if(u.flood.b$x[i]<180){u.flood.b$x[i]<-u.flood.b$x[i]-((360/144)/2)}
        if(u.flood.b$x[i]>180){u.flood.b$x[i]<-u.flood.b$x[i]-360+((360/144)/2)}
      }
      #u.flood.b <- subset(u.flood.b, layer==1)
      index.b <- which(u.flood.b$layer==1)
      wind.data.bu <- wind_data[index.b,]
      
      #Prepares the bootstrap results in a raster and dataframe format for ggplot
      v.r.b <- raster(t(v.results.bootstrap2), xmn=vlon[1], xmx=vlon[144], ymn=vlat[1], ymx=vlat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      v.flood.b <- raster::as.data.frame(v.r.b, xy = TRUE)
      cells <- length(v.flood.b$x)
      for(i in 1:cells){ 
        if(v.flood.b$x[i]<180){v.flood.b$x[i]<-v.flood.b$x[i]-((360/144)/2)}
        if(v.flood.b$x[i]>180){v.flood.b$x[i]<-v.flood.b$x[i]-360+((360/144)/2)}
      }
      #v.flood.b <- subset(v.flood.b, layer==1)
      index.b <- which(v.flood.b$layer==1)
      wind.data.bv <- wind_data[index.b,]
      
      q <- ggplot()+
        metR::geom_contour_fill(data=wind_data, aes(x=x, y=y, z=u),binwidth = 0.01)+
        #geom_point(data=u.flood.b, aes(x=x, y=y), color="black", size=0.5)+
        geom_segment(data=wind.data.bu, aes(x=x, xend = x+(u*5), y = y, yend = y), arrow = arrow(length = unit(0.15, "cm")))+
        geom_point(data=houston, aes(x=x, y=y), color="red")+
        borders('world', xlim=range(wind_data$x), ylim=range(wind_data$y), colour='black', size=.2)+
        borders('state', colour="black", size=0.1)+
        theme_minimal()+ 
        theme(axis.title.x = element_text(size = 13),axis.text.x = element_text(size = 12),axis.title.y = element_text(size = 13),axis.text.y = element_text(size = 12))+
        theme(panel.ontop=TRUE, panel.background=element_blank()) +
        scale_fill_gradientn(colors =myPallette1, limits=c(-1.2,1.2), breaks=mybreak,
                             guide = guide_colorsteps(frame.colour = "black", ticks.colour = "black", barwidth=30, barheight=0.8, title.position="top",title="Standardized wind speed anomoly", 
                                                      title.hjust=0.5, title.theme=element_text(size=13),label.theme=element_text(size=12)))+
        coord_quickmap(xlim=c(-120,-40), ylim=c(10,60))+
        labs(x="Longitude", y="Latitude",fill=paste("", sep=""))+
        theme(legend.position = "bottom")+
        theme(plot.title = element_text(hjust = 0, face="bold", size=13))+
        theme(legend.margin=margin(20, 0, 0, 0))+
        labs(title=paste(letter1[tijd], "   Zonal wind - ",tijd_naam[tijd], sep=""))
      
      p <- ggplot()+
        metR::geom_contour_fill(data=wind_data, aes(x=x, y=y, z=v),binwidth = 0.01)+
        #geom_point(data=v.flood.b, aes(x=x, y=y), color="black", size=0.5)+
        geom_segment(data=wind.data.bv, aes(x=x, xend = x, y = y, yend = y+(v*5)), arrow = arrow(length = unit(0.15, "cm")))+
        geom_point(data=houston, aes(x=x, y=y), color="red")+
        borders('world', xlim=range(wind_data$x), ylim=range(wind_data$y), colour='black', size=.2)+
        borders('state', colour="black", size=0.1)+
        theme_minimal()+ 
        theme(axis.title.x = element_text(size = 13),axis.text.x = element_text(size = 12),axis.title.y = element_text(size = 13),axis.text.y = element_text(size = 12))+
        theme(panel.ontop=TRUE, panel.background=element_blank()) +
        scale_fill_gradientn(colors =myPallette1, limits=c(-1.2,1.2), breaks=mybreak,
                             guide = guide_colorsteps(frame.colour = "black", ticks.colour = "black", barwidth=30, barheight=0.8, title.position="top",title="Standardized wind speed anomoly", 
                                                      title.hjust=0.5, title.theme=element_text(size=13),label.theme=element_text(size=12)))+
        coord_quickmap(xlim=c(-120,-40), ylim=c(10,60))+
        labs(x="Longitude", y="Latitude",fill=paste("", sep=""))+
        theme(legend.position = "bottom")+
        theme(plot.title = element_text(hjust = 0, face="bold", size=13))+
        theme(legend.margin=margin(20, 0, 0, 0))+
        labs(title=paste(letter2[tijd],"   Meridional wind - ",tijd_naam[tijd], sep=""))
      
      plot_list[[tijd]] <- q
      plot_list[[tijd+3]] <- p
      
    }
   

combinedA <- plot_list[[7]] / plot_list[[8]]  + plot_layout(ncol = 1) & theme(legend.position = "bottom")
combinedB <- plot_list[[1]] + plot_list[[4]] + plot_list[[2]] + plot_list[[5]] + plot_list[[3]] + plot_list[[6]] + plot_layout(ncol = 2) & theme(legend.position = "bottom")
FA <- combinedA / plot_layout(guides = "collect")
FB <- combinedB + plot_layout(guides = "collect")


#combinedC <-  plot_list[[7]] + plot_spacer() + plot_list[[1]] + plot_list[[2]] + plot_list[[3]] +
#              plot_list[[8]] + plot_spacer() + plot_list[[4]] + plot_list[[5]] + plot_list[[6]] +
#              guide_area()   + plot_spacer() + plot_spacer()  + guide_area()   + plot_spacer()  + plot_layout(ncol = 5) & theme(legend.position = "bottom")
#combinedC + plot_layout(guides = "collect", widths = c(1,0.1,0.5,0.5,0.5))
#combinedD <- (plot_list[[7]] / plot_list[[8]]) | ((plot_list[[1]] + plot_list[[2]] + plot_list[[3]])/ (plot_list[[4]] + plot_list[[5]] + plot_list[[6]])) & theme(legend.position = "bottom")
#combinedD + plot_layout(guides = "collect")

combinedE <- FA | FB  + plot_layout(widths = c(0.9,1))
combinedE 

