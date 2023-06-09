#Author: Joeri Reinders
#Purpose: Script which creates composite maps of geopotentialheight reanalysis data for days of flooding [given by an csv document] - figure 3 in paper.
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
myPallette <- c(rev(brewer.pal(11, "RdBu")))
myPallette2 <- c((brewer.pal(9, "GnBu"))[2:7],brewer.pal(9,"BuPu")[5:7])

##############################
## Geopotential Height load ##
##############################

directory <- "Data/Reanalysis/NCEP-NCAR/GeopotentialHeight/daily_data/"  #change this into the folder you need it to be
file_name <- paste('All.',level,'.nc',sep="")
var <- "hgt"

nc <- nc_open(paste(drive_e,directory,file_name, sep=""))
attributes(nc$var)

lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat", verbose = F)
time <- ncvar_get(nc, "time")

ndvi.array <- ncvar_get(nc, var) # store the data in a 3-dimensional array
dim(ndvi.array) 

fillvalue <- ncatt_get(nc, var, "_FillValue")
ndvi.array[ndvi.array == fillvalue$value] <- NA


maand_titel <- c("(d)  April", "(e) May", "(f)  June")
maand <- c("04","05","06")

for(m in 1:3){ #Loop generates 3 composistion maps of geopotential heigh during 1000 random spring days in April, May, June 
dates <- seq(as.Date("1948-01-01"), as.Date("2020-12-31"), by = "days") 
months <- format(dates, "%m")
dates_lente <- which((months==maand[m]))
general_map.index <- sample(dates_lente,1000,replace=TRUE)
general_map.data <- ndvi.array[,,general_map.index]
general_map.slice <- apply(general_map.data,1:2, mean_fun)

general_map.r <- raster(t(general_map.slice), xmn=lon[1], xmx=lon[144], ymn=lat[1], ymx=lat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
general_map.flood <- raster::as.data.frame(general_map.r, xy = TRUE)
cells <- length(general_map.flood$x)
for(i in 1:cells){ 
  if(general_map.flood$x[i]<180){general_map.flood$x[i]<-general_map.flood$x[i]-((360/144)/2)}
  if(general_map.flood$x[i]>180){general_map.flood$x[i]<-general_map.flood$x[i]-360+((360/144)/2)}
}

plot_list[[3+m]] <- ggplot()+
  metR::geom_contour_fill(data=general_map.flood, aes(x=x, y=y, z=layer),binwidth = 10)+
  geom_point(data=houston, aes(x=x, y=y), color="red")+
  borders('world', xlim=range(general_map.flood$x), ylim=range(general_map.flood$y), colour='black', size=.2)+
  borders('state', colour="black", size=0.1)+
  theme_minimal()+ 
  theme(panel.ontop=TRUE, panel.background=element_blank()) +
  theme(axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 13),axis.title.y = element_text(size = 14),axis.text.y = element_text(size = 13))+
  scale_fill_gradientn(colors =myPallette2, limits=c(1300,1700), breaks=c(1300,1350,1400,1450,1500,1550,1600,1650,1700),
                       guide = guide_colorsteps(frame.colour = "black", ticks.colour = "black", barwidth=40, barheight=0.8, title.position="top",title="850 hPa geopotential height (m)",
                                                title.hjust=0.5, title.theme=element_text(size=14),label.theme=element_text(size=13)))+
  coord_quickmap(xlim=c(-160,-20), ylim=c(0,75))+
  labs(x="Longitude", y="Latitude",fill=paste("", sep=""))+
  theme(legend.position = "bottom")+
  theme(plot.title = element_text(hjust = 0, face="bold", size=14))+
  theme(legend.margin=margin(20, 0, 0, 0))+
  labs(title=maand_titel[m])

}

    ###################
    ## Daily Anomoly ## Computes daily anomolies with 1981-2010 as climate normal
    ###################
    
    ndvi.normal.m <- array(NA, dim=c(144,73,366))
    ndvi.normal.sd <- array(NA, dim=c(144,73,366))
    
    dates <- seq(as.Date("1948-01-01"), as.Date("2020-12-31"), by = "days") 
    index.max <- which(dates=="2020-12-31")
    date_dates <- format(dates, "%m-%d")
    
    dates_anom <- seq(as.Date("1981-01-01"), as.Date("2010-12-31"), by = "days")
    index_anom.min <- which((dates==as.Date("1981-01-01")))
    date_dates_anom <- format(dates_anom, "%m-%d")
    
    dates_ref <- format((seq(as.Date("1948-01-01"), as.Date("1948-12-31"), by = "days")),"%m-%d") 
    
    for(day in 1:366){
      date_index <- (which(dates_ref[day]==date_dates_anom)+index_anom.min)
      ndvi.normal <- ndvi.array[,,date_index]
      ndvi.normal.m[,,day] <- apply(ndvi.normal,1:2, mean_fun)
      ndvi.normal.sd[,,day] <- apply(ndvi.normal,1:2, sd_fun)
    }
    
    for(days in 1:366){
      date_index <- which(dates_ref[days]==date_dates)
      ndvi.array[,,date_index] <- sweep(ndvi.array[,,date_index],1:2,ndvi.normal.m[,,days])/replicate(length(date_index), ndvi.normal.sd[,,days], simplify="array")
    }



mybreak <- seq(-1,1,by=0.2)
tijd_naam <- c("(a)  Day", "(b)  Week", "(c)  Three Weeks")

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
alpha <- 0.01

for(tijd in 1:3){
  
  directory <- "Users/joeri/OneDrive/Bureaublad/Northeastern University/Research Chapters/Chapter 4/"
  file_10 <- paste(drive_l,directory,"Discharg_data/Top10Spring.csv",sep="")
  Spring_Top10 <- read.csv(file_10)
  
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
    ndvi.slices <- ndvi.array[,,index_day]
    ndvi.slice <- apply(ndvi.slices,1:2, mean_fun)
    
    results.bootstrap2 <- bootstrap_function(10000,30,ndvi.array,ndvi.slice,1,alpha)
    
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
    ndvi.slices <- ndvi.array[,,index_week]
    ndvi.slice <- apply(ndvi.slices,1:2, mean_fun)
    
    results.bootstrap2 <- bootstrap_function(10000,30,ndvi.array,ndvi.slice,7,alpha)
    
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
    
    ndvi.slices <- ndvi.array[,,index_3week] #netcdfs of the picked days
    ndvi.slice <- apply(ndvi.slices,1:2, mean_fun)
    
    results.bootstrap2 <- bootstrap_function(10000,30,ndvi.array,ndvi.slice,21,alpha)
    
  }
  
  #Prepares composistion map for ggplot
  r <- raster(t(ndvi.slice), xmn=lon[1], xmx=lon[144], ymn=lat[1], ymx=lat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  # first the raster is molded into a dataframe which can be plotted in ggplot and the lot is adjusted for the world map
  flood <- raster::as.data.frame(r, xy = TRUE)
  cells <- length(flood$x)
  for(i in 1:cells){ 
    if(flood$x[i]<180){flood$x[i]<-flood$x[i]-((360/144)/2)}
    if(flood$x[i]>180){flood$x[i]<-flood$x[i]-360+((360/144)/2)}
  }
  
  #Prepares the bootstrap results in a raster and dataframe format for ggplot
  r.b <- raster(t(results.bootstrap2), xmn=lon[1], xmx=lon[144], ymn=lat[1], ymx=lat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  flood.b <- raster::as.data.frame(r.b, xy = TRUE)
  cells <- length(flood.b$x)
  for(i in 1:cells){ 
    if(flood.b$x[i]<180){flood.b$x[i]<-flood.b$x[i]-((360/144)/2)}
    if(flood.b$x[i]>180){flood.b$x[i]<-flood.b$x[i]-360+((360/144)/2)}
  }
  flood.b <- subset(flood.b, layer==1)

  q <- ggplot()+
    metR::geom_contour_fill(data=flood, aes(x=x, y=y, z=layer),binwidth = 0.05)+
    geom_point(data=flood.b, aes(x=x, y=y), color="black", size=0.5)+
    geom_point(data=houston, aes(x=x, y=y), color="red")+
    borders('world', xlim=range(flood$x), ylim=range(flood$y), colour='black', size=.2)+
    borders('state', colour="black", size=0.1)+
    theme_minimal()+ 
    theme(panel.ontop=TRUE, panel.background=element_blank()) +
    theme(axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 13),axis.title.y = element_text(size = 14),axis.text.y = element_text(size = 13))+
    scale_fill_gradientn(colors =myPallette, limits=c(-1.2,1.2), breaks=mybreak,
                         guide = guide_colorsteps(frame.colour = "black", ticks.colour = "black", barwidth=40, barheight=0.8, title.position="top",title="Standardized 850 hPa geopotential height anomoly", 
                                                  title.hjust=0.5, title.theme=element_text(size=14),label.theme=element_text(size=13)))+
    coord_quickmap(xlim=c(-160,-20), ylim=c(0,75))+
    labs(x="Longitude", y="Latitude",fill=paste("", sep=""))+
    theme(legend.position = "bottom")+
    theme(plot.title = element_text(hjust = 0, face="bold", size=14))+
    theme(legend.margin=margin(20, 0, 0, 0))+
    labs(title=tijd_naam[tijd])
  
  plot_list[[tijd]] <- q
  
}


library(patchwork)

combineda <- plot_list[[4]] + plot_list[[5]] + plot_list[[6]] & theme(legend.position = "bottom")
combinedb <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] & theme(legend.position = "bottom")
(combinedb + plot_layout(guides = "collect")) / (combineda + plot_layout(guides = "collect")) 

