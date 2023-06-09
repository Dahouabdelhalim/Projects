library(ggplot2)
library(modeest)
library(tools)
library(dplyr)
library(lubridate)

setwd("****************")

#Read in all file names in the folder
file.names <- list.files(full.names=FALSE, pattern = ".csv")

#Use only to work on one file at a time
file <- "F500000004512253_072919_629965.csv"

work.file <- read.csv(file)
work.file_name <- substr(file, nchar(file)-9, nchar(file)-4)

#Use to calculate HI for all files in the folder
#sapply(file.names,function(y){
#  work.file <- read.csv(y, header = TRUE)
#  work.file_name <- substr(y, nchar(y)-9, nchar(y)-4)

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

#Plots raw data, active mode, and overall mode values to verify that mode makes sense
ggplot(work.file, aes(Date_Time)) +
  geom_line(aes(Date_Time,Temp), color = "gray") +
  geom_hline(yintercept = mode_Tb, color = "blue", linetype = "dotted") +
  geom_hline(yintercept = active_mode, color = "red")


#Plots a histogram and the mode value to verify that mode makes sense
ggplot(high_mode, aes(high_mode$Temp)) + 
  geom_histogram(binwidth = 0.1) + 
  geom_vline(xintercept = mode_Tb, color = "blue") +
  geom_vline(xintercept = active_mode, color = "red")



##########################################################
##########################################################
###########Loop through days to calculate HI##############
##########################################################
##########################################################

#creating an empty data frame called HIsub with the above column names, so that "HIsub" exists for the for loop below
HIsub <- data.frame(matrix(ncol=11, nrow=0))
x <- c("Animal", "Date","Overall_mode", "Overall_active_mode", "Daily_mode", "Daily_active_mode", "Max_Tb", "Min_Tb", "Mean_Tb", "HI_Overall_Mode", "HI_Daily_Mode")
colnames(HIsub) <- x

#creating an empty data frame called temp with the above column names, so that "temp" exists for the for loop below
temp <- data.frame(matrix(ncol=11, nrow=1))
x <- c("Animal", "Date","Overall_mode", "Overall_active_mode", "Daily_mode", "Daily_active_mode", "Max_Tb", "Min_Tb", "Mean_Tb", "HI_Overall_Mode", "HI_Daily_Mode")
colnames(temp) <- x

maxdate<- (max(high_mode$Date))
mindate<- (min(high_mode$Date))
Days_measured<-maxdate-mindate


for(i in mindate:maxdate){
  temp_Tb_data <- subset(high_mode, high_mode$Date == i)
  daily_mode_Tb <- max(mlv(temp_Tb_data$Tb_round, method = "mfv"))
  daily_active_period <- with(temp_Tb_data, temp_Tb_data[hour(Date_Time) >=23 | hour(Date_Time) <= 3,])
  daily_active_mode <- max(mlv(daily_active_period$Tb_round, method = "mfv"))
  
  temp_Tb_data$diff_overall <- active_mode - temp_Tb_data$Temp
  sum_overall <- sum(temp_Tb_data$diff_overall^2)
  temp_Tb_data$diff_daily <- daily_active_mode - temp_Tb_data$Temp
  sum_daily <- sum(temp_Tb_data$diff_daily^2)
  count <- nrow(temp_Tb_data)
  
  temp$Animal <- work.file_name
  temp$Date <- as.Date(i, format = "%m-%d-%Y", origin = "01-01-1970") 
  temp$Overall_mode <- mode_Tb
  temp$Overall_active_mode <- active_mode
  temp$Daily_mode <- daily_mode_Tb
  temp$Daily_active_mode <- daily_active_mode
  temp$Max_Tb <- max(temp_Tb_data$Temp)
  temp$Min_Tb <- min(temp_Tb_data$Temp)
  temp$Mean_Tb <- mean(temp_Tb_data$Temp)
  temp$HI_Overall_Mode <- sqrt((sum_overall/(nrow(temp_Tb_data)-1))) 
  temp$HI_Daily_Mode <- sqrt((sum_daily/(nrow(temp_Tb_data)-1)))
  
  
  if(i == mindate){HIsub <- temp}else{HIsub <- rbind(HIsub,temp)}
  
}


head(HIsub)
str(HIsub) 
print(work.file_name)
write.csv(HIsub, file=paste("./Calculated HI/", work.file_name, ".csv", sep=""))


HI_graph <- ggplot(HIsub, aes(Date)) +
  geom_line(aes(Date,HI_Overall_Mode), color = "blue") +
  geom_line(aes(Date,HI_Daily_Mode), color = "red") 
  
  print(HI_graph)