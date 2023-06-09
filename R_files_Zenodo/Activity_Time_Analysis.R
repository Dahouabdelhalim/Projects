library(ggplot2)
library(modeest)
library(tools)
library(dplyr)
library(lubridate)
library(suncalc)

setwd("****************")

#Read in all file names in the folder
file.names <- list.files(full.names=FALSE, pattern = ".csv")

file <- "2A00000008331453_072919_630033.csv"

#(file.names,function(y){
#work.file <- read.csv(y, header = TRUE)
#work.file_name <- substr(y, nchar(y)-9, nchar(y)-4)

#Create moonrise and moonset data
start <- as.Date("2019-03-01", origin = "1960-10-01")
end <- as.Date("2019-08-01", origin = "1960-10-01")

moontimes<- getMoonTimes(date = seq.Date(start, end, by = 1),
                         keep = c("rise", "set"),
                         lat=34.33, lon=-106.70, tz="Etc/GMT+7")

#Create sunrise and sunset data
start <- as.Date("2019-03-01", origin = "1960-10-01")
end <- as.Date("2019-08-01", origin = "1960-10-01")

suntimes<- getSunlightTimes(date = seq.Date(start, end, by = 1),
                         keep = c("sunset", "sunrise", "night"),
                         lat=34.33, lon=-106.70, tz="Etc/GMT+7")


write.csv(moontimes, file=paste("./Calculated Time Active/", "moontimes.csv", sep=""))
write.csv(suntimes, file=paste("./Calculated Time Active/", "suntimes.csv", sep=""))


work.file <- read.csv(file)
work.file_name <- substr(file, nchar(file)-9, nchar(file)-4)



#Format Excel date as an R date format
work.file$Date_Time <- as.POSIXct(work.file$Date_Time, format = "%m/%d/%Y %H:%M")


#Round off temperature data for calculating mode below. If there are too many decimals, mode calculation can be odd.
work.file$Tb_round <- round(work.file$Temp, digits = 2)

#Pulls out the real time (not offset time) to subset data by active period
work.file$Real_time <- strftime(work.file$Date_Time, format = "%H:%M:%S")

#Set time backwards or forwards a certain number of hours (12 would set the break between days at noon). This serves to get each active period within one day.
offset_hours<- 12

work.file$Offset_Time <- work.file$Date_Time - (3600 * offset_hours)

#Pull date from the offset date_time column for analysis
work.file$Date <- as.Date(format(work.file$Offset_Time, "%Y-%m-%d"))

#Used to calculate higher mode if data are bimodal. If so, set the cutoff value higher than the lower mode. If the data are unimodal, set this value lower than the lowest temp.ONLY USE THIS IF THE TEMPERATURES DO NOT FOLLOW A DAILY CYCLE. IF THEY DO, USE THE NEXT LINES TO CUT DATA TO ONLY ACTIVE PERIODS!
cutoff_Temp <- 15
high_mode <- subset(work.file, work.file$Temp > cutoff_Temp)

#Calculate overall mode for the body temperature dataset
mode_Tb <- mlv(high_mode$Tb_round, method = "mfv")

#Calculate mode for active period
active_period <- with(high_mode, high_mode[hour(Date_Time) >=23 | hour(Date_Time) <= 3,])
active_mode <- mlv(active_period$Tb_round, method = "mfv")


##########################################################
##########################################################
###########Loop through days to calculate Activity times##
##########################################################
##########################################################

#creating an empty data frame called timesub with the above column names, so that "timesub" exists for the for loop below
timesub <- data.frame(matrix(ncol=12, nrow=0))
x <- c("Animal", "Date", "Overall_active_mode", "Cutoff", "Activity_Start", "Activity_End", "Activity_Period", "Sunset", "Sunrise", "Night", "Moonset", "Moonrise")
colnames(timesub) <- x

temp <- data.frame(matrix(ncol=12, nrow=1))
x <- c("Animal", "Date", "Overall_active_mode", "Cutoff", "Activity_Start", "Activity_End", "Activity_Period", "Sunset", "Sunrise", "Night", "Moonset", "Moonrise")
colnames(temp) <- x

maxdate<- (max(high_mode$Date))
mindate<- (min(high_mode$Date))
Days_measured<-maxdate-mindate

cutoff <- .5


for(i in mindate:maxdate){
  temp_Tb_data <- subset(high_mode, high_mode$Date == i)
 # temp_Tb_data_active <- subset(temp_Tb_data_warm, temp_Tb_data_warm[hour(Date_Time) >=18 | hour(Date_Time) <= 6,])
  temp_Tb_data_warm <- subset(temp_Tb_data, temp_Tb_data$Temp >= (active_mode - cutoff))
 # daily_active_period <- with(temp_Tb_data_warm, temp_Tb_data_warm[hour(Date_Time) >=18 | hour(Date_Time) <= 6,])
  suntime1 <- subset(suntimes,suntimes$date == i)
  suntime2 <- subset(suntimes,suntimes$date == i+1)
  ms <- subset(moontimes,moontimes$set>suntime1$sunset & moontimes$set<suntime2$sunrise)
  mr <- subset(moontimes,moontimes$rise>suntime1$sunset & moontimes$rise<suntime2$sunrise)
  
  temp$Animal <- work.file_name
  temp$Date <- as.Date(i, format = "%m-%d-%Y", origin = "01-01-1970") 
  temp$Overall_active_mode <- active_mode
  temp$Cutoff <- cutoff
  temp$Activity_Start <- min(temp_Tb_data_warm$Date_Time)
  temp$Activity_End <- max(temp_Tb_data_warm$Date_Time)
  temp$Activity_Period <- difftime(temp$Activity_End, temp$Activity_Start, units = c("hours"))
  temp$Sunset<- as.POSIXct(suntime1$sunset, format = "%m/%d/%Y %H:%M", origin = "01-01-1970")
  temp$Sunrise<- as.POSIXct(suntime2$sunrise, format = "%m/%d/%Y %H:%M", origin = "01-01-1970")
  temp$Night <- as.POSIXct(suntime1$night, format = "%m/%d/%Y %H:%M", origin = "01-01-1970")
  if(nrow(ms)==1){temp$Moonset <- ms$set}else{temp$Moonset <- as.POSIXct(NA, format = "%m/%d/%Y %H:%M", origin = "01-01-1970")}
  if(nrow(mr)==1){temp$Moonrise <- mr$rise}else{temp$Moonrise <- as.POSIXct(NA,format = "%m/%d/%Y %H:%M", origin = "01-01-1970")}
  temp$Moonrise <- temp$Moonrise - 7200
  
  if(i == mindate){timesub <- temp}else{timesub <- rbind(timesub,temp)}
  
}

head(timesub)
str(timesub) 
print(work.file_name)
write.csv(timesub, file=paste("./Calculated Time Active/", work.file_name, ".csv", sep=""))


Time_graph <- ggplot(timesub, aes(Date)) +
  geom_point(aes(Date,Activity_Period), color = "blue")

print(Time_graph)

#})