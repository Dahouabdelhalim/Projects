library(sp)
library(rgeos)
library(dplyr)

setwd("C:/Users/Mohammad Fahmy/Desktop/Honours 2016-2017/Operation Kinglet/Localisation Script/Dryad") # set directory

####Localisation####
buffer <- 10 # set buffer width
error_limit <- 0.1 # the maximum error value to use when excluding locs
error_limit_micdist <- 0.06 # maximum error value to use when excluding locs for structural analysis song selection
maxmicdist <- 15 # max distance from a microphone that is acceptable for structural analysis

miclist <- read.csv("Microphone_Coordinates.csv") # read master mic coordinate file; 1=SW1, 2=NW1 or NC1, 3=NE1, 4=SE1
miclist <- miclist[order(miclist$Array, miclist$Microphone), ] # sort miclist by array & mic

data <- read.csv("RCKI_MasterFile.csv") # open data
data <- data[order(data$clipname), ] # sort data by clipname
data <- data[data$is.song == "Y",] # remove non-songs

data.array <- data.frame(matrix(nrow = max(data$array), ncol = 33)) # initialize dataframe with 1 row per array
colnames(data.array) <- c("array", "day", "week", "year" ,"area", "Ncalls", "Ninarray", "Nlowerror", "Nfinal", "hr1", "hr2", "hr3", "hr4", "hr5", "hr6", "hr7", "hr8", "hr9", "hr10", "hr11", "hr12", "hr13", "hr14", "hr15", "hr16", "hr17", "hr18", "hr19", "hr20", "hr21", "hr22", "hr23", "hr24") # hour values are for songs inside array with low error

mic.dist <- data.frame() # initialize dataframe to store locs that are within maxmicdist

for (i in 1:max(data$array)) { # run loop for each array

	data.temp <- data[data$array == i, ] # extract all detections for array i
	if(nrow(data.temp) == 0) next # skip iteration if there are no songs in data.temp
	
	mics <- miclist[miclist$Array == i, ] # extract mic coordinates for detection
	mics2 <- as.matrix(mics[, 5:4]) # make duplicate of mic coordinates for calculating loc to microphone distances
	poly <- as.matrix(mics[, 5:4]) # convert xy coordinates to matrix
	poly <- rbind(poly, poly[1,]) # add first row to end of matrix to create polygon
	poly <- Polygon(poly) # specify as polygon
	poly <- Polygons(list(poly),1) # add wrapper to polygon class
	poly <- SpatialPolygons(list(poly)) # add wrapper??
	poly <- gBuffer(poly, width = buffer, joinStyle="ROUND") # add buffer

	data.array$array[i] <- data.temp$array[1] # add array id to data.array
	data.array$day[i] <- data.temp$julien.day[1]
	data.array$week[i] <- data.temp$week[1]
	data.array$year[i] <- data.temp$year[1]
	data.array$area[i] <- gArea(poly) # add area of polygon to data.array
	data.array$Ncalls[i] <- nrow(data.temp) # print number of songs before filtering to data.array

	
	for (j in 1:nrow(data.temp)) { # run new loop for each song in data.temp
		
		coordinates <- cbind(data.temp$East_BF[j], data.temp$North_BF[j]) # specify coordinates for loc
		point <- SpatialPointsDataFrame(coordinates, data.temp[j, ]) # convert coordinates for loc to spatial point

		data.temp$in_array[j] <- gCovers(poly, point) # is point contained in polygon
		time <- data.temp$parent.start.time[j] # extract start time of song's parent file
		time <- strptime(time, "%H:%M:%S") # convert start time of song's parent file to time format
		time <- time + data.temp$song.start.time[j] # add start time of song within parent file to parent file start time
		data.temp$hour1[j] <- as.character(time, format = "%H") # print hour in which song was produced to data.temp
		
		mics3 <- rbind(coordinates, mics2) # make matrix of x,y coordinates with loc at top and mics below
		mics3 <- as.matrix(dist(mics3)) # calculate distance matrix from mics3 and then convert to class matrix
		data.temp$minmicdist[j] <- min(mics3[2:nrow(mics3), 1]) # extract min distance between point and each microphone
		data.temp$minmic[j] <- which(grepl(min(mics3[2:nrow(mics), 1]), mics3[,1])) # identify which mic has min dist to loc
		
	}

	# print all rows where loc was within maxmicdist of a mic, and SumAbsDiff_BF <= 0.05
	mic.dist <- rbind(mic.dist, data.temp[data.temp$minmicdist <= maxmicdist & data.temp$SumAbsDiff_BF <= error_limit_micdist & data.temp$is.solo == "Y", c("clipname", "array", "minmicdist", "minmic")])

	data.array$Ninarray[i] <- nrow(data.temp[data.temp$in_array == "TRUE", ]) # print number of songs produced inside array to data.array
	data.array$Nlowerror[i] <- nrow(data.temp[data.temp$SumAbsDiff_BF <= error_limit, ]) # print number of songs with acceptable error to data.array
	data.array$Nfinal[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit, ]) # print number of songs that are in array and with acceptable error
	
	# print number of songs inside array with acceptable error that were produced during each hour of day
	data.array$hr1[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "00", ])
	data.array$hr2[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "01", ])
	data.array$hr3[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "02", ])
	data.array$hr4[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "03", ])
	data.array$hr5[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "04", ])
	data.array$hr6[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "05", ])
	data.array$hr7[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "06", ])
	data.array$hr8[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "07", ])
	data.array$hr9[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "08", ])
	data.array$hr10[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "09", ])
	data.array$hr11[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "10", ])
	data.array$hr12[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "11", ])
	data.array$hr13[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "12", ])
	data.array$hr14[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "13", ])
	data.array$hr15[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "14", ])
	data.array$hr16[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "15", ])
	data.array$hr17[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "16", ])
	data.array$hr18[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "17", ])
	data.array$hr19[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "18", ])
	data.array$hr20[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "19", ])
	data.array$hr21[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "20", ])
	data.array$hr22[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "21", ])
	data.array$hr23[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "22", ])
	data.array$hr24[i] <- nrow(data.temp[data.temp$in_array == "TRUE" & data.temp$SumAbsDiff_BF <= error_limit & data.temp$hour1 == "23", ])

}


data.array[, 6:33] <- data.array[, 6:33] * 4 # correct song #s for 1 in 4 subsampling
data.array[, 5] <- data.array[, 5]/10000 #convert m^2 to ha
data.array[, 6:33] <- data.array[, 6:33] / data.array[, 5] # convert song #s to # songs per hectare

###Generating Plots#### 

#load required packages
library(ggplot2)
library(plotrix) # to use std.error function
library(ggpubr)  # to arrange plots next to each other



###Remove arrays that do not contain songs
df<-data.array
df2<-df%>%
  filter(Nfinal!=0)

####Seasonal pattern####

###Weekly pattern 2017
year2 <- df2 %>% 
  filter(year == 2017)%>%
  group_by(week)%>%
  summarize(songs.per.ha=mean(Nfinal, na.rm = TRUE),se= std.error(Nfinal, na.rm = TRUE))
year2$n <- c(5,6,5,7,4)


yr2<-2017
wksong2017<-cbind(year2, yr=yr2)

#Graph weekly pattern 2017
plot_2017 <- ggplot(year2, aes(x = week, y=songs.per.ha)) +
              xlab("Week") + ylab("Songs  per hectare per day")+
                geom_point(se=FALSE)+
                geom_errorbar(data = year2,aes(ymin=songs.per.ha-se, ymax=songs.per.ha+se),
                              width=0.1,position=position_dodge(1))+ 
                scale_y_continuous(name = "", limits = c(0,2000))+
                geom_ribbon(data = year2, aes(ymin=songs.per.ha-se, ymax=songs.per.ha+se), alpha = 0.2)+
                geom_text(x=21, y=1290, label = "5")+
                geom_text(x=22, y=958, label = "6")+
                geom_text(x=24 , y=247 , label = "5")+
                geom_text(x=25 , y=356 , label = "7")+
                geom_text(x=26 , y=230 , label = "4")+
                theme_bw(15) +
                theme_classic()+
                ggtitle("2017")


###Weekly pattern 2016
year1 <-df2%>%
  filter(year == 2016)%>%
  group_by(week)%>%
  summarize(songs.per.ha=mean(Nfinal, na.rm = TRUE),se= std.error(Nfinal, na.rm = TRUE))
year1$n <-c(9,9,8,9,5,9,4,1)

yr1<-2016
wksong2016<-cbind(year1, yr=yr1)


###Graph weekly pattern 2016
plot_2016 <- ggplot(year1, aes(x = week, y=songs.per.ha)) +
  geom_point(se=FALSE)+
  geom_errorbar(data = year1,aes(ymin=songs.per.ha-se, ymax=songs.per.ha+se),
                width=0.1,position=position_dodge(1))+ 
  geom_ribbon(data = year1, aes(ymin=songs.per.ha-se, ymax=songs.per.ha+se), alpha = 0.2)+
  scale_x_continuous(name = "Week", breaks = c(21,22,23,24,25,26,26,27,28))+
  scale_y_continuous(name = "Songs per hectare per day", limits = c(0,2000))+
  geom_text(x=21, y=689, label = "9")+
  geom_text(x=22, y=409, label = "9")+
  geom_text(x=23 , y=1076 , label = "8")+
  geom_text(x=24 , y=1631 , label = "9")+
  geom_text(x=25 , y=680 , label = "5")+
  geom_text(x=26 , y=745 , label = "9")+
  geom_text(x=27 , y=246 , label = "4")+
  geom_text(x=28 , y=101 , label = "1")+
  theme_bw(15) +
  theme_classic()+ 
  ggtitle("2016")


### Arrange existing plots next to each other using ggarrange

pdf(file = "C:/Users/Mohammad Fahmy/Desktop/Honours 2016-2017/Operation Kinglet/Localisation Script/Dryad/Fig.2 Seasonal pattern plot.pdf", width=4.7, height=4)
ggarrange(plot_2016 , plot_2017)
dev.off()



#### Diel pattern with both years combined####

byhrcombo<- df2[,4:33]

data.array.combo <- data.frame(matrix(nrow = 24, ncol = 3))
colnames(data.array.combo)<- c("hr", "mean.songs", "error")
data.array.combo$hr <- factor(c(1:24))
data.array.combo$hr <- factor(data.array.combo$hr,labels = c("01","02","03","04","05","06","07","08","09","10","11","12",
                                                             "13","14","15","16","17","18","19","20","21","22","23","00"))

data.array.combo$hr <- as.integer(data.array.combo$hr)

byhrcombo<-byhrcombo[,-1:-6]
data.array.combo$mean.songs<-sapply(byhrcombo, mean)
data.array.combo$error <-sapply(byhrcombo, std.error)
data.array.combo<- data.array.combo[c(24,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),]

pdf(file = "C:/Users/Mohammad Fahmy/Desktop/Honours 2016-2017/Operation Kinglet/Localisation Script/Dryad/Fig.3 Diel pattern.pdf", width= 6.5, height= 4)

ggplot(data.array.combo, aes(x = hr ,y = mean.songs))+
  xlab("Hour of the day") + ylab("Number of songs  per hectare")+
  geom_point()+
  scale_shape_manual(values=16)+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24))+
  geom_errorbar(data = data.array.combo,aes(ymin=mean.songs-error, ymax=mean.songs+error),
                width=0.1,position=position_dodge(1))+ 
  geom_vline(xintercept=c(3.6,4,19.9,20.6), linetype="dotted")+
  geom_vline(xintercept=c(2.7,3.3,20.7,21.4), linetype = "longdash")+
  geom_ribbon(data = data.array.combo, aes(ymin = mean.songs-error, ymax = mean.songs+error), alpha = 0.2)+
  theme_bw(15) +
  theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))+
  theme_classic()

dev.off()

