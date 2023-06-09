#Author: Joeri Reinders
#Purpose: Script which extracts the ten largest daily discharge levels from the Evadale, Romayor and Richmond USGS gages on respectively the Neches, Trinity and Brazos Rivers.

rm(list=ls())

#packages:
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
library(ggplotify)
library(cowplot)

x11()

rm(list=ls())

drive <- 'C:/'
directory <- "Users/joeri/OneDrive/Bureaublad/Northeastern University/Research Chapters/Chapter 4/"
file_neches <- paste(drive,directory,"Discharg_data/Daily_Neches@Evadale.csv",sep="")
file_trinity <- paste(drive,directory,"Discharg_data/Daily_Trinity@Romayor.csv",sep="")
file_brazos <- paste(drive,directory,"Discharg_data/Daily_Brazos@Richmond.csv",sep="")
orgNeches <- read.csv(file_neches)
orgTrinity <- read.csv(file_trinity)
orgBrazos <- read.csv(file_brazos)

which(orgNeches$Date=="1/1/1949")

spring_neches <- subset(orgNeches, Month==4|Month==5|Month==6)
spring_trinity <- subset(orgTrinity, Month==4|Month==5|Month==6)
spring_brazos <- subset(orgBrazos, Month==4|Month==5|Month==6)

fall_neches <- subset(orgNeches, Month==8|Month==9|Month==10)
fall_trinity <- subset(orgTrinity, Month==8|Month==9|Month==10)
fall_brazos <- subset(orgBrazos, Month==8|Month==9|Month==10)

SS_Neches <- spring_neches[order(spring_neches$Discharge, decreasing = T),]
SS_Trinity <- spring_trinity[order(spring_trinity$Discharge, decreasing = T),]
SS_Brazos <- spring_brazos[order(spring_brazos$Discharge, decreasing = T),]

FF_Neches <- fall_neches[order(fall_neches$Discharge, decreasing = T),]
FF_Trinity <- fall_trinity[order(fall_trinity$Discharge, decreasing = T),]
FF_Brazos <- fall_brazos[order(fall_brazos$Discharge, decreasing = T),]


## Check whether you need fall (FF) or spring (SS) dates
top10_neches <- rep(NA, 10)
Index <- 1

while(anyNA(top10_neches)==T){
  
  print(Index)
  d1 <- as.Date((as.character(SS_Neches$Date[Index])), format="%m/%d/%Y")
  
  if(sum(is.na(top10_neches))==10){top10_neches[1] <- as.character(d1)}
  
  vol <- abs(sum(is.na(top10_neches))-11)
  lengte <- vol - 1 
  dif <- rep(NA,lengte)
  
  for(index in 1:lengte){
    d2 <- (as.Date((as.character(top10_neches[index]))))
    dif[index] <- abs(as.numeric(difftime(d2,d1 , units="day")))
  }
  
  if(any(dif<16)){Index <- Index+1}
  else{top10_neches[vol] <- as.character(d1)
       Index <- Index+1
       }
  
}

top10_trinity <- rep(NA, 10)
Index <- 1

while(anyNA(top10_trinity)==T){
  print(Index)
  d1 <- as.Date((as.character(SS_Trinity$Date[Index])), format="%m/%d/%Y")
  if(sum(is.na(top10_trinity))==10){top10_trinity[1] <- as.character(d1)}
  
  vol <- abs(sum(is.na(top10_trinity))-11)
  lengte <- vol - 1 
  dif <- rep(NA,lengte)
  
  for(index in 1:lengte){
    d2 <- (as.Date((as.character(top10_trinity[index]))))
    dif[index] <- abs(as.numeric(difftime(d2,d1 , units="day")))
  }
  
  if(any(dif<16)){Index <- Index+1}
  else{top10_trinity[vol] <- as.character(d1)
  Index <- Index+1
  }
  
}

top10_brazos <- rep(NA, 10)
Index <- 1

while(anyNA(top10_brazos)==T){
  print(Index)
  d1 <- as.Date((as.character(SS_Brazos$Date[Index])), format="%m/%d/%Y")
  if(sum(is.na(top10_brazos))==10){top10_brazos[1] <- as.character(d1)}
  
  vol <- abs(sum(is.na(top10_brazos))-11)
  lengte <- vol - 1 
  dif <- rep(NA,lengte)
  
  for(index in 1:lengte){
    d2 <- (as.Date((as.character(top10_brazos[index]))))
    dif[index] <- abs(as.numeric(difftime(d2,d1 , units="day")))
  }
  
  if(any(dif<16)){Index <- Index+1}
  else{top10_brazos[vol] <- as.character(d1)
  Index <- Index+1
  }
  
}


print(top10_neches)
print(top10_trinity)
print(top10_brazos)

Top10Spring <- data.frame(Neches = top10_neches,Trinity=top10_trinity,Brazos=top10_brazos)

setwd(paste(drive,directory,"Discharg_data/",sep=""))
getwd()

write.csv(Top10Spring,"Top10Spring.csv")


Neches1973 <- subset(orgNeches, Year==1973)
plot(as.Date((as.character(Neches1973$Date)), format="%m/%d/%Y"),Neches1973$Discharge, type ="l")

dev.off()