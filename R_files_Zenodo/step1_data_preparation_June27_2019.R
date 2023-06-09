lapply(c("dplyr","stringr","geosphere","ggmap"), require, character.only = TRUE)

data_13_14<-read.csv("2013-2014 Balearics for JP_June27_2019.csv")
data_13_14[data_13_14==""] <- NA 

clean_13_14<- data_13_14[,c(1:9,11:12,14)]

clean_13_14$Distance_m[is.na(clean_13_14$Distance_m)] <- 0
clean_13_14$Angle[is.na(clean_13_14$Angle)] <- 0

clean13<-subset(clean_13_14,clean_13_14$Year=="2013")
clean14<-subset(clean_13_14,clean_13_14$Year=="2014")

clean13$Time<-paste(str_pad(clean13$Time, 5, pad = "0"),"00",sep=":")

clean14<-subset(clean14,!is.na(clean14$Time))
clean14$Time<-format(strptime(clean14$Time, "%I:%M:%S %p"), format="%H:%M:%S")

all_data<-read.csv("All 2015-2017 Balearics for JP jessy edits.csv")
all_data[all_data==""] <- NA 
all_data<-subset(all_data,is.na(all_data$Duplicate))
names(all_data)[1] <- "ESASorML"

clean_all<- all_data[,c(1:9,11:12,14)]

clean15<-subset(clean_all,clean_all$Year=="2015")
clean16<-subset(clean_all,clean_all$Year=="2016")
clean17<-subset(clean_all,clean_all$Year=="2017")


clean15$Time<-ifelse(nchar(as.character(clean15$Time))<6, paste(clean15$Time,":00",sep=""), paste(clean15$Time))
clean16$Time<-ifelse(nchar(as.character(clean16$Time))<6, paste(clean16$Time,":00",sep=""), 
                     format(strptime(clean16$Time, "%I:%M:%S %p"), format="%H:%M:%S")) 
clean17$Time<-paste(str_pad(clean17$Time, 5, pad = "0"),"00",sep=":")
clean17$Date<-paste(sub('/([^/]*)$', '', clean17$Date),"/17",sep="")

clean17<-subset(clean17,as.Date( clean17$Date, "%m/%d/%y")< as.Date("2017-10-25"))


all_new_data <- rbind(clean13,clean14,clean15,clean16,clean17)
all_new_data$datetime<-as.POSIXct(paste(all_new_data$Date, all_new_data$Time, sep = " "),
                             format = "%m/%d/%y %H:%M:%S", tz = "gmt")


all_gps<-c()
for (i in 1:length(all_new_data[,1])){
  one_gps<-destPoint(c(all_new_data[i,7],all_new_data[i,6]), all_new_data[i,10], all_new_data[i,11]) 
  all_gps<-rbind(all_gps,one_gps)
}

all_new_data$bird_lon<-all_gps[,1]
all_new_data$bird_lat<-all_gps[,2]

all_clean<-subset(all_new_data,!is.na(all_new_data$bird_lat)) 

all_sighting_nums<- all_clean[c(13,1,2,4,12,14:15,6:7,11)]
names(all_sighting_nums)[8:9] <- c("boat_lat","boat_lon")

bird_on<-subset(all_sighting_nums,all_sighting_nums$Effort=="Effort")

setwd("C:/Users/Jessica/Desktop/BS analysis none out 1km")
saveRDS(bird_on, file = "0.0 all_sighting_on_effort_June3 2020.rds")


