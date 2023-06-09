### Formatting weather data to use in analysis
#
# I'm just going to use billy barr's weather station data for all sites. It's the most complete
# It would have been ideal to use Jen Rudgers' interpolation code, but that only goes to 2015 >:(
#
# Michael Stemkovski 2018-11-27
#
###

setwd("/home/michael/Documents/Grad School/Research Projects/Bee phenology")

library(lubridate)
library(readxl)
library(measurements)

read.excel <- function(...) return(as.data.frame(read_excel(...)))

#weather_raw <- read.csv("Raw data/RMBL weather data/billy_bar_wrcc.csv",stringsAsFactors = F)
weather_raw <- read.excel("Raw data/RMBL weather data/aawx Michael.xlsx") #this came from billy barr's weather station
bare_ground <- read.csv("Raw data/RMBL weather data/billy_barr_bare_ground.csv") #this came from billy barr's observations too

# basic cleaning
bare_ground$doy <- yday(mdy(bare_ground$date_bare_ground))
weather <- weather_raw[3:nrow(weather_raw),]
colnames(weather) <- c("year","month","day","temp_min","temp_max","snow_cm","snow_water_in","snow_depth_cm","rain_in")
weather$snow_water_in <- conv_unit(as.numeric(weather$snow_water_in),"inch","cm")
weather$rain_in <- conv_unit(as.numeric(weather$rain_in),"inch","cm")
colnames(weather) <- c("water_year","wy_month","day","temp_min","temp_max","snow_cm","snow_water_cm","snow_depth_cm","rain_cm")
weather$rain_cm <- sapply(weather$rain_cm,function(x) ifelse(is.na(x), return(0), return(x)))

# date reformatting
get.year <- function(water_year,wy_month){
  if(wy_month >= 9 & wy_month <= 12){
    year <- paste("20",substr(water_year,nchar(water_year)-3,nchar(water_year)-2),sep="")
  } else {
    year <- paste("20",substr(water_year,nchar(water_year)-1,nchar(water_year)),sep="")
  }
  return(as.numeric(year))
}

head(weather,20)
str(weather)

month_letters <- toupper(letters[1:12])
month_nums <- c(9:12,1:8)

month_letters
month_nums

weather$month <- sapply(weather$wy_month, function(x) return(month_nums[which(month_letters == x)]), USE.NAMES = F)
weather$year <- apply(weather[,c("water_year","month")], 1, function(x) get.year(x[1], as.numeric(x[2])))
weather$date <- ymd(apply(weather[,c("year","month","day")], 1, function(x) paste(x[1],x[2],x[3],sep="-")))

years <- 2009:2018
water_years <- unique(weather$water_year)[3:12]
weather_summary <- data.frame(year = years,
                              snow = rep(NA,length(years)),
                              snow_water = rep(NA,length(years)),
                              rain = rep(NA,length(years)),
                              bare_ground = rep(NA,length(years)),
                              temp_min = rep(NA,length(years)),
                              temp_max = rep(NA,length(years)),
                              stringsAsFactors = F)

for (i in 1:length(years)){
  year <- years[i]
  w_year <- water_years[i]
  
  water_year_sub <- weather[which(weather$water_year == w_year),]
  weather_summary$snow[i] <- sum(as.numeric(water_year_sub$snow_cm))
  weather_summary$snow_water[i] <- sum(as.numeric(water_year_sub$snow_water_cm))
  #weather_summary$rain[i] <- sum(as.numeric(water_year_sub$rain_cm))
  
  summer_sub <- weather[which(weather$year == year & weather$month %in% 4:9),]
  #weather_summary$temp_min[i] <- mean(as.numeric(water_year_sub$temp_min))
  #weather_summary$temp_max[i] <- mean(as.numeric(water_year_sub$temp_max))
  weather_summary$temp_min[i] <- mean(as.numeric(summer_sub$temp_min))
  weather_summary$temp_max[i] <- mean(as.numeric(summer_sub$temp_max))
  weather_summary$rain[i] <- sum(as.numeric(summer_sub$rain_cm))
  
  weather_summary$bare_ground[i] <- bare_ground[which(bare_ground$year == year),"doy"]
  
}

plot(weather_summary)

write.csv(weather_summary,"Edited data/weather.csv", row.names = F)

# par(mar=c(4.5, 4, 3, 2))
# plot(weather$temp_min ~ weather$temp_max)
# plot(weather$snow_water_cm ~ weather$snow_cm)
# plot(weather$snow_depth_cm ~ weather$temp_min,pch=20,col=rgb(0.2,0.5,0.9,0.2))
# plot(weather$snow_depth_cm ~ weather$temp_max)
# plot(weather$snow_water_cm ~ weather$temp_min)
# plot(weather$snow_cm ~ weather$temp_min)
# plot(weather$rain_cm ~ weather$temp_max)

# colnames(weather_raw) <- 1:ncol(weather_raw)
# weather_raw[1:4,]
# nrow <- nrow(weather_raw)
# 
# # dealing with the annoying thing where precip switches columns halways though the dataset
# get.precip <- function(x){
#   x <- sapply(x,interpret.na)
#   precip <- sum(x[1],x[2],na.rm = T) # this is a placeholder
#   return(precip)
# }
# 
# interpret.na <- function(x){
#   if(x == -9999) x <- NA
#   return(as.numeric(x))
# }
# sum(2,NA)
# 
# weather_raw[10400:10402,19:20]
# as.matrix(weather_raw[10400:10402,19:20])
# apply(as.matrix(weather_raw[10400:10402,19:20]),c(1,2),interpret.na)
# apply(as.matrix(weather_raw[10400,19:20]),c(1,2),interpret.na)
# apply(as.matrix(weather_raw[10400:10402,19:20]),1,get.precip) # problem
# #precip = sapply(lapply(weather_raw[4:nrow,19:20],as.numeric),get.precip),
# 
# weather <- data.frame(date = weather_raw[4:nrow,1],
#                       min_temp = as.numeric(weather_raw[4:nrow,13]),
#                       max_temp = as.numeric(weather_raw[4:nrow,12]),
#                       avg_temp = as.numeric(weather_raw[4:nrow,14]),
#                       snow_depth = as.numeric(weather_raw[4:nrow,21]),
#                       precip = apply(as.matrix(weather_raw[4:nrow,19:20]),1,get.precip),
#                       stringsAsFactors = F)
# 
# weather <- data.frame(date = weather_raw[4:nrow,1],
#                       min_temp = lapply(weather_raw[4:nrow,13],as.numeric),
#                       max_temp = lapply(weather_raw[4:nrow,12],as.numeric),
#                       avg_temp = lapply(weather_raw[4:nrow,14],as.numeric),
#                       stringsAsFactors = F)
# 
# years <- 2009:2017
# 
# weather_summary <- data.frame(year = years,
#                               )

###### need to make sure that the days dont have unequal numbers of observations... this may skew averages if say morning measures are missing
###### also, looks like precip and acumm precip switched columns at some point >:(