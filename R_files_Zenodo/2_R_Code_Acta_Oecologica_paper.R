## Here we are not taking nest sites as source location for seed kernels.
#also we are looking the proportion of seeds falling in nest,roost and other sites

library(lubridate)
library(sp)
library(geosphere)
library(plotrix)

##Import movement data from csv file
data<-read.csv("a_data_GH1Br_sdk.csv",header = TRUE)
attach(data)
##Remove garbage data
data<-data[data$Long != 0 & data$Lat != 0,]
##Import gut passage data
gpt<-read.csv("regur2.csv",header=TRUE)
gpt
attach(gpt)
##Activity budget file - we will use it for sampling during bootstrapping
act=read.csv("proportions_foraging.csv",header=TRUE)

##Only consider data from days when more than 30 locations are available (to avoid underestimationof dispersal)
temp <- table(data$Date)
temp<-names(temp[temp[] > 30])
data <- subset(data,Date %in% temp)
##dectale dataframe to store the results
final = data.frame(0)

###We have to make use of nest and roost locations (change this parameter according to the focal individual)
nest_lat = 26.938917 #26.936533 #
nest_long = 92.973317 #92.967067 #


## Remove the data points 
#new data frame 
#data1 <- data

data$nest = numeric(length=nrow(data))
data$roost =numeric(length=nrow(data))
k = 1
i = 1

#Find out nest locations and remove them
for(i in 1:nrow(data))
{
  lat = data[i,"Lat"]
  long = data[i,"Long"]
  dist_nest = distHaversine(c(long,lat), c(nest_long,nest_lat))
  if(dist_nest < 20)
    data$nest[i] = 1
  #dist_roost = distHaversine(c(long,lat), c(roost_long,roost_lat))
  #if(dist_roost < 50)
  #  data$roost[i] = 1
}
data1 <- data[data$nest != 1,]
if(nrow(data1) == 0)
{
  message <- "All points lie within nest area"
  message
}

##Create a vector to represent each row of data frame as many times as activity level in that hour
times=numeric(length=nrow(data1))
for(i in 1:nrow(data1)){
  times[i]=round(act$Breeding[max(which((strptime(act$Time.of.the.day,format="%H:%M:%S") <= strptime(data1$Time[i],format="%H:%M:%S"))))])
  
}
times[is.na(times)] = 0
samp_vec = rep(1:nrow(data1),times=times)

##calculation for every time bin
#final$nest = 0
#final$roost = 0
##Loop over all gpts
for(i in 1:length(gpt$time))
{
  #Select first bin from GPT
  tim = gpt$bin[i]
  j = 1
  while(j <= 10)
  {
    #Sampling according to activity levels
    ran <- sample(samp_vec,1)
    ##Select source location and time
    latx = data1[ran,"Lat"]
    lonx = data1[ran,"Long"]  
    day = data1[ran,"Date"] 
    reftime = strptime(as.character(data1$Stamp[ran]),format = "%d/%m/%y %H:%M")
    reftime = format(reftime,"%H:%M")
    if(reftime>=format(strptime("17:00","%H:%M"),"%H:%M"))
      next
    t = 0
    temp = which(data$Location == data1$Location[ran])
    if(temp <= nrow(data))
      while(data[temp,"Date"]  == day)
      {
        temp = temp+1
        if(temp >= nrow(data))
          break
        
        dif = as.duration(strptime(data$Time[temp],format = "%H:%M") - strptime(reftime,format="%H:%M")) 
        dif = as.numeric(dif/60)
        if(!is.na(dif) & length(dif) != 0 & dif <= tim & dif > tim-1)
        {
          t = temp
          break
        }
        
        
      }
    if(t == 0)                        ##take last value if destination was not found within given interval
      t = temp-1
    
    laty = data[t,"Lat"]
    lony = data[t,"Long"]
    #dist = distHaversine(c(lony,laty), c(lonx,latx))
    
    dist = distHaversine(c(lony,laty), c(lonx,latx))
    dist_nest = distHaversine(c(lony,laty), c(nest_long,nest_lat))
    if(dist_nest < 20)
      final[k,"nest"] = dist
    
    
    #dist_roost = distHaversine(c(lony,laty), c(roost_long,roost_lat))
    #if(dist_roost < 50)
    #final[k,"roost"] = dist
    
    
    if(dist_nest >= 20)
      final[k,"other"] = dist
    
    final[k,"lat"] = latx
    final[k,"long"] = lonx
    #  final[k,"dist"] = dist
    # final[k,"time"] = gpt$time
    # final[k,"prop"] = gpt$Prop[u]
    #final[k,"flag"] = flag
    k = k + 1
    j=j+1
  }
}
write.csv(final$nest,"sdk_output_nest_bill_17Sep18_whgh_10.csv")
write.csv(final$other,"sdk_output_other_bill_17Sep18_whgh_10.csv")
final$nest<-cut(final$nest,breaks=c(0,20,150,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),dig.lab=4)
final$other<-cut(final$other,breaks=c(0,20,150,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),dig.lab=4)
tiff(file = "Dispersal_prob_BILL_17Sep18_whgh_10.tiff",width = 1420, height = 840,pointsize = 18)
x = barplot(rbind(table(final$nest)/nrow(final),table(final$other)/nrow(final)),beside = FALSE,space = 2,legend = c("Nest","Other locations"),col = c("black","gray"),xlab = "Distance(m)",ylab = "Probability of seed arrival")
dev.off()


####GH3Br#######

library(lubridate)
library(sp)
library(geosphere)
library(plotrix)

##Import movement data from csv file
data<-read.csv("a_data_GH3Br_sdk.csv",header = TRUE)
attach(data)
##Remove garbage data
data<-data[data$Long != 0 & data$Lat != 0,]
##Import gut passage data
gpt<-read.csv("regur2.csv",header=TRUE)
gpt
attach(gpt)
##Activity budget file - we will use it for sampling during bootstrapping
act=read.csv("proportions_foraging.csv",header=TRUE)

##Only consider data from days when more than 30 locations are available (to avoid underestimationof dispersal)
temp <- table(data$Date)
temp<-names(temp[temp[] > 30])
data <- subset(data,Date %in% temp)
##dectale dataframe to store the results
final = data.frame(0)

###We have to make use of nest and roost locations (change this parameter according to the focal individual)
nest_lat = 26.95423
nest_long = 92.94208

## Remove the data points 
#new data frame 
#data1 <- data

data$nest = numeric(length=nrow(data))
data$roost =numeric(length=nrow(data))
k = 1
i = 1

#Find out nest locations and remove them
for(i in 1:nrow(data))
{
  lat = data[i,"Lat"]
  long = data[i,"Long"]
  dist_nest = distHaversine(c(long,lat), c(nest_long,nest_lat))
  if(dist_nest < 20)
    data$nest[i] = 1
  #dist_roost = distHaversine(c(long,lat), c(roost_long,roost_lat))
  #if(dist_roost < 50)
  #  data$roost[i] = 1
}
data1 <- data[data$nest != 1,]
if(nrow(data1) == 0)
{
  message <- "All points lie within nest area"
  message
}

##Create a vector to represent each row of data frame as many times as activity level in that hour
times=numeric(length=nrow(data1))
for(i in 1:nrow(data1)){
  times[i]=round(act$Breeding[max(which((strptime(act$Time.of.the.day,format="%H:%M:%S") <= strptime(data1$Time[i],format="%H:%M:%S"))))])
  
}
times[is.na(times)] = 0
samp_vec = rep(1:nrow(data1),times=times)

##calculation for every time bin
#final$nest = 0
#final$roost = 0
##Loop over all gpts
for(i in 1:length(gpt$time))
{
  #Select first bin from GPT
  tim = gpt$bin[i]
  j = 1
  while(j <= 10)
  {
    #Sampling according to activity levels
    ran <- sample(samp_vec,1)
    ##Select source location and time
    latx = data1[ran,"Lat"]
    lonx = data1[ran,"Long"]  
    day = data1[ran,"Date"] 
    reftime = strptime(as.character(data1$Stamp[ran]),format = "%d/%m/%y %H:%M")
    reftime = format(reftime,"%H:%M")
    if(reftime>=format(strptime("17:00","%H:%M"),"%H:%M"))
      next
    t = 0
    temp = which(data$Location == data1$Location[ran])
    if(temp <= nrow(data))
      while(data[temp,"Date"]  == day)
      {
        temp = temp+1
        if(temp >= nrow(data))
          break
        
        dif = as.duration(strptime(data$Time[temp],format = "%H:%M") - strptime(reftime,format="%H:%M")) 
        dif = as.numeric(dif/60)
        if(!is.na(dif) & length(dif) != 0 & dif <= tim & dif > tim-1)
        {
          t = temp
          break
        }
        
      }
    if(t == 0)                        ##take last value if destination was not found within given interval
      t = temp-1
    
    laty = data[t,"Lat"]
    lony = data[t,"Long"]
    #dist = distHaversine(c(lony,laty), c(lonx,latx))
    
    dist = distHaversine(c(lony,laty), c(lonx,latx))
    dist_nest = distHaversine(c(lony,laty), c(nest_long,nest_lat))
    if(dist_nest < 20)
      final[k,"nest"] = dist
    
    
    #dist_roost = distHaversine(c(lony,laty), c(roost_long,roost_lat))
    #if(dist_roost < 50)
    #final[k,"roost"] = dist
    
    
    if(dist_nest >= 20)
      final[k,"other"] = dist
    
    final[k,"lat"] = latx
    final[k,"long"] = lonx
    #  final[k,"dist"] = dist
    # final[k,"time"] = gpt$time
    # final[k,"prop"] = gpt$Prop[u]
    #final[k,"flag"] = flag
    k = k + 1
    
    j=j+1
  }
  
}
write.csv(final$nest,"sdk_output_nest_godfather_17Sep18_whgh_10.csv")
write.csv(final$other,"sdk_output_other_godfather_17Sep18_whgh_10.csv")
final$nest<-cut(final$nest,breaks=c(0,20,150,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),dig.lab=4)
final$other<-cut(final$other,breaks=c(0,20,150,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),dig.lab=4)
tiff(file = "Dispersal_prob_GODFATHER_17Sep18_whgh_10.tiff",width = 1420, height = 840,pointsize = 18)
x = barplot(rbind(table(final$nest)/nrow(final),table(final$other)/nrow(final)),beside = FALSE,space = 2,legend = c("Nest","Other locations"),col = c("black","gray"),xlab = "Distance(m)",ylab = "Probability of seed arrival")
dev.off()

####GH4Br (Code name: Rifle) #######

library(lubridate)
library(sp)
library(geosphere)
library(plotrix)

##Import movement data from csv file
data<-read.csv("a_data_GH4Br_sdk.csv",header = TRUE)
attach(data)
##Remove garbage data
data<-data[data$Long != 0 & data$Lat != 0,]
##Import gut passage data
gpt<-read.csv("regur2.csv",header=TRUE)
gpt
attach(gpt)
##Activity budget file - we will use it for sampling during bootstrapping
act=read.csv("proportions_foraging.csv",header=TRUE)

##Only consider data from days when more than 30 locations are available (to avoid underestimationof dispersal)
temp <- table(data$Date)
temp<-names(temp[temp[] > 30])
data <- subset(data,Date %in% temp)
##dectale dataframe to store the results
final = data.frame(0)

###We have to make use of nest locations (change this parameter according to the focal individual)
nest_lat = 26.94431	
nest_long = 92.95374

## Remove the data points 
#new data frame 
#data1 <- data

data$nest = numeric(length=nrow(data))
data$roost =numeric(length=nrow(data))
k = 1
i = 1

#Find out nest locations and remove them
for(i in 1:nrow(data))
{
  lat = data[i,"Lat"]
  long = data[i,"Long"]
  dist_nest = distHaversine(c(long,lat), c(nest_long,nest_lat))
  if(dist_nest < 20)
    data$nest[i] = 1
  #dist_roost = distHaversine(c(long,lat), c(roost_long,roost_lat))
  #if(dist_roost < 50)
  #  data$roost[i] = 1
}
data1 <- data[data$nest != 1,]
if(nrow(data1) == 0)
{
  message <- "All points lie within nest area"
  message
}

##Create a vector to represent each row of data frame as many times as activity level in that hour
times=numeric(length=nrow(data1))
for(i in 1:nrow(data1)){
  times[i]=round(act$Breeding[max(which((strptime(act$Time.of.the.day,format="%H:%M:%S") <= strptime(data1$Time[i],format="%H:%M:%S"))))])
  
}
times[is.na(times)] = 0
samp_vec = rep(1:nrow(data1),times=times)

##calculation for every time bin
#final$nest = 0
#final$roost = 0
##Loop over all gpts
for(i in 1:length(gpt$time))
{
  #Select first bin from GPT
  tim = gpt$bin[i]
  j = 1
  while(j <= 10)
  {
    #Sampling according to activity levels
    ran <- sample(samp_vec,1)
    ##Select source location and time
    latx = data1[ran,"Lat"]
    lonx = data1[ran,"Long"]  
    day = data1[ran,"Date"] 
    reftime = strptime(as.character(data1$Stamp[ran]),format = "%d/%m/%y %H:%M")
    reftime = format(reftime,"%H:%M")
    if(reftime>=format(strptime("17:00","%H:%M"),"%H:%M"))
      next
    t = 0
    temp = which(data$Location == data1$Location[ran])
    if(temp <= nrow(data))
      while(data[temp,"Date"]  == day)
      {
        temp = temp+1
        if(temp >= nrow(data))
          break
        
        dif = as.duration(strptime(data$Time[temp],format = "%H:%M") - strptime(reftime,format="%H:%M")) 
        dif = as.numeric(dif/60)
        if(!is.na(dif) & length(dif) != 0 & dif <= tim & dif > tim-1)
        {
          t = temp
          break
        }
        
      }
    if(t == 0)                        ##take last value if destination was not found within given interval
      t = temp-1
    
    laty = data[t,"Lat"]
    lony = data[t,"Long"]
    #dist = distHaversine(c(lony,laty), c(lonx,latx))
    
    dist = distHaversine(c(lony,laty), c(lonx,latx))
    dist_nest = distHaversine(c(lony,laty), c(nest_long,nest_lat))
    if(dist_nest < 20)
      final[k,"nest"] = dist
    
    
    #dist_roost = distHaversine(c(lony,laty), c(roost_long,roost_lat))
    #if(dist_roost < 50)
    #final[k,"roost"] = dist
    
    
    if(dist_nest >= 20)
      final[k,"other"] = dist
    
    final[k,"lat"] = latx
    final[k,"long"] = lonx
    #  final[k,"dist"] = dist
    # final[k,"time"] = gpt$time
    # final[k,"prop"] = gpt$Prop[u]
    #final[k,"flag"] = flag
    k = k + 1
    
    j=j+1
  }
  
}
write.csv(final$nest,"sdk_output_nest_rifle_17Sep18_whgh_10.csv")
write.csv(final$other,"sdk_output_other_rifle_17Sep18_whgh_10.csv")
final$nest<-cut(final$nest,breaks=c(0,20,150,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),dig.lab=4)
final$other<-cut(final$other,breaks=c(0,20,150,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),dig.lab=4)
tiff(file = "Dispersal_prob_RIFLE_17Sep18_whgh_10.tiff",width = 1420, height = 840,pointsize = 18)
x = barplot(rbind(table(final$nest)/nrow(final),table(final$other)/nrow(final)),beside = FALSE,space = 2,legend = c("Nest","Other locations"),col = c("black","gray"),xlab = "Distance(m)",ylab = "Probability of seed arrival")
dev.off()


####WH1Br#######
library(lubridate)
library(sp)
library(geosphere)
library(plotrix)

##Import movement data from csv file
data<-read.csv("a_data_WH1Br_sdk.csv",header = TRUE)
attach(data)
##Remove garbage data
data<-data[data$Long != 0 & data$Lat != 0,]
##Import gut passage data
gpt<-read.csv("regur2.csv",header=TRUE)
gpt
attach(gpt)
##Activity budget file - we will use it for sampling during bootstrapping
act=read.csv("proportions_foraging.csv",header=TRUE)

##Only consider data from days when more than 30 locations are available (to avoid underestimationof dispersal)
temp <- table(data$Date)
temp<-names(temp[temp[] > 30])
data <- subset(data,Date %in% temp)
##dectale dataframe to store the results
final = data.frame(0)

###We have to make use of nest and roost locations (change this parameter according to the focal individual)
nest_lat = 26.936533	
nest_long = 92.967067


## Remove the data points 
#new data frame 
#data1 <- data

data$nest = numeric(length=nrow(data))
data$roost =numeric(length=nrow(data))
k = 1
i = 1

#Find out nest locations and remove them
for(i in 1:nrow(data))
{
  lat = data[i,"Lat"]
  long = data[i,"Long"]
  dist_nest = distHaversine(c(long,lat), c(nest_long,nest_lat))
  if(dist_nest < 20)
    data$nest[i] = 1
  #dist_roost = distHaversine(c(long,lat), c(roost_long,roost_lat))
  #if(dist_roost < 50)
  #  data$roost[i] = 1
}
data1 <- data[data$nest != 1,]
if(nrow(data1) == 0)
{
  message <- "All points lie within nest area"
  message
}

##Create a vector to represent each row of data frame as many times as activity level in that hour
times=numeric(length=nrow(data1))
for(i in 1:nrow(data1)){
  times[i]=round(act$Breeding[max(which((strptime(act$Time.of.the.day,format="%H:%M:%S") <= strptime(data1$Time[i],format="%H:%M:%S"))))])
  
}
times[is.na(times)] = 0
samp_vec = rep(1:nrow(data1),times=times)

##calculation for every time bin
#final$nest = 0
#final$roost = 0
##Loop over all gpts
for(i in 1:length(gpt$time))
{
  #Select first bin from GPT
  tim = gpt$bin[i]
  j = 1
  while(j <= 10)
  {
    #Sampling according to activity levels
    ran <- sample(samp_vec,1)
    ##Select source location and time
    latx = data1[ran,"Lat"]
    lonx = data1[ran,"Long"]  
    day = data1[ran,"Date"] 
    reftime = strptime(as.character(data1$Stamp[ran]),format = "%d/%m/%y %H:%M")
    reftime = format(reftime,"%H:%M")
    if(reftime>=format(strptime("17:00","%H:%M"),"%H:%M"))
      next
    t = 0
    temp = which(data$Location == data1$Location[ran])
    if(temp <= nrow(data))
      while(data[temp,"Date"]  == day)
      {
        temp = temp+1
        if(temp >= nrow(data))
          break
        
        dif = as.duration(strptime(data$Time[temp],format = "%H:%M") - strptime(reftime,format="%H:%M")) 
        dif = as.numeric(dif/60)
        if(!is.na(dif) & length(dif) != 0 & dif <= tim & dif > tim-1)
        {
          t = temp
          break
        }
        
      }
    if(t == 0)                        ##take last value if destination was not found within given interval
      t = temp-1
    
    laty = data[t,"Lat"]
    lony = data[t,"Long"]
    #dist = distHaversine(c(lony,laty), c(lonx,latx))
    
    dist = distHaversine(c(lony,laty), c(lonx,latx))
    dist_nest = distHaversine(c(lony,laty), c(nest_long,nest_lat))
    if(dist_nest < 20)
      final[k,"nest"] = dist
    
    
    #dist_roost = distHaversine(c(lony,laty), c(roost_long,roost_lat))
    #if(dist_roost < 50)
    #final[k,"roost"] = dist
    
    
    if(dist_nest >= 20)
      final[k,"other"] = dist
    
    final[k,"lat"] = latx
    final[k,"long"] = lonx
    #  final[k,"dist"] = dist
    # final[k,"time"] = gpt$time
    # final[k,"prop"] = gpt$Prop[u]
    #final[k,"flag"] = flag
    k = k + 1
    
    j=j+1
  }
  
}
write.csv(final$nest,"sdk_output_nest_gabbar_17Sep18_whgh_10.csv")
write.csv(final$other,"sdk_output_other_gabbar_17Sep18_whgh_10.csv")
final$nest<-cut(final$nest,breaks=c(0,20,150,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),dig.lab=4)
final$other<-cut(final$other,breaks=c(0,20,150,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),dig.lab=4)
tiff(file = "Dispersal_prob_GABBAR_17Sep18_whgh_10.tiff",width = 1420, height = 840,pointsize = 18)
x = barplot(rbind(table(final$nest)/nrow(final),table(final$other)/nrow(final)),beside = FALSE,space = 2,legend = c("Nest","Other locations"),col = c("black","gray"),xlab = "Distance(m)",ylab = "Probability of seed arrival")
dev.off()


####GH2NBr#######
library(lubridate)
library(sp)
library(geosphere)
library(plotrix)

##Import movement data from csv file
data<-read.csv("a_data_GH2NBr_sdk.csv",header = TRUE)
attach(data)
##Remove garbage data
data<-data[data$Long != 0 & data$Lat != 0,]
##Import gut passage data
gpt<-read.csv("regur2.csv",header=TRUE)
gpt
attach(gpt)
##Activity budget file - we will use it for sampling during bootstrapping
act=read.csv("proportions_foraging.csv",header=TRUE)

##Only consider data from days when more than 30 locations are available (to avoid underestimationof dispersal)
temp <- table(data$Date)
temp<-names(temp[temp[] > 30])
data <- subset(data,Date %in% temp)
##dectale dataframe to store the results
final = data.frame(0)

###We have to make use of nest and roost locations (change this parameter according to the focal individual)
nest_lat = 1.1111
nest_long = 1.1111


## Remove the data points 
#new data frame 
#data1 <- data

data$nest = numeric(length=nrow(data))
data$roost =numeric(length=nrow(data))
k = 1
i = 1

#Find out nest locations and remove them
for(i in 1:nrow(data))
{
  lat = data[i,"Lat"]
  long = data[i,"Long"]
  dist_nest = distHaversine(c(long,lat), c(nest_long,nest_lat))
  if(dist_nest < 20)
    data$nest[i] = 1
  #dist_roost = distHaversine(c(long,lat), c(roost_long,roost_lat))
  #if(dist_roost < 50)
  #  data$roost[i] = 1
}
data1 <- data[data$nest != 1,]
if(nrow(data1) == 0)
{
  message <- "All points lie within nest area"
  message
}

##Create a vector to represent each row of data frame as many times as activity level in that hour
times=numeric(length=nrow(data1))
for(i in 1:nrow(data1)){
  times[i]=round(act$Nonbreeding[max(which((strptime(act$Time.of.the.day,format="%H:%M:%S") <= strptime(data1$Time[i],format="%H:%M:%S"))))])
  
}
times[is.na(times)] = 0
samp_vec = rep(1:nrow(data1),times=times)

##calculation for every time bin
#final$nest = 0
#final$roost = 0
##Loop over all gpts
for(i in 1:length(gpt$time))
{
  #Select first bin from GPT
  tim = gpt$bin[i]
  j = 1
  while(j <= 10)
  {
    #Sampling according to activity levels
    ran <- sample(samp_vec,1)
    ##Select source location and time
    latx = data1[ran,"Lat"]
    lonx = data1[ran,"Long"]  
    day = data1[ran,"Date"] 
    reftime = strptime(as.character(data1$Stamp[ran]),format = "%d/%m/%y %H:%M")
    reftime = format(reftime,"%H:%M")
    if(reftime>=format(strptime("17:00","%H:%M"),"%H:%M"))
      next
    t = 0
    temp = which(data$Location == data1$Location[ran])
    if(temp <= nrow(data))
      while(data[temp,"Date"]  == day)
      {
        temp = temp+1
        if(temp >= nrow(data))
          break
        
        dif = as.duration(strptime(data$Time[temp],format = "%H:%M") - strptime(reftime,format="%H:%M")) 
        dif = as.numeric(dif/60)
        if(!is.na(dif) & length(dif) != 0 & dif <= tim & dif > tim-1)
        {
          t = temp
          break
        }
        
      }
    if(t == 0)                        ##take last value if destination was not found within given interval
      t = temp-1
    
    laty = data[t,"Lat"]
    lony = data[t,"Long"]
    #dist = distHaversine(c(lony,laty), c(lonx,latx))
    
    dist = distHaversine(c(lony,laty), c(lonx,latx))
    dist_nest = distHaversine(c(lony,laty), c(nest_long,nest_lat))
    if(dist_nest < 20)
      final[k,"nest"] = dist
    
    
    #dist_roost = distHaversine(c(lony,laty), c(roost_long,roost_lat))
    #if(dist_roost < 50)
    #final[k,"roost"] = dist
    
    
    if(dist_nest >= 20)
      final[k,"other"] = dist
    
    final[k,"lat"] = latx
    final[k,"long"] = lonx
    #  final[k,"dist"] = dist
    # final[k,"time"] = gpt$time
    # final[k,"prop"] = gpt$Prop[u]
    #final[k,"flag"] = flag
    k = k + 1
    
    j=j+1
  }
  
}
write.csv(final$nest,"sdk_output_nest_mogambo_17Sep18_whgh_10.csv")
write.csv(final$other,"sdk_output_other_mogambo_17Sep18_whgh_10.csv")
final$other<-cut(final$other,breaks=c(0,20,150,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),dig.lab=4)
tiff(file = "Dispersal_prob_MOGAMBO_17Sep18_whgh_10.tiff",width = 1420, height = 840,pointsize = 18)
x = barplot(table(final$other)/nrow(final),beside = FALSE,space = 2,legend = "Other locations",col = "gray",xlab = "Distance(m)",ylab = "Probability of seed arrival")
dev.off()


####TKBHAI#######
library(lubridate)
library(sp)
library(geosphere)
library(plotrix)

##Import movement data from csv file
data<-read.csv("a_data_GH5NBr_sdk.csv",header = TRUE)
attach(data)
##Remove garbage data
data<-data[data$Long != 0 & data$Lat != 0,]
##Import gut passage data
gpt<-read.csv("regur2.csv",header=TRUE)
gpt
attach(gpt)
##Activity budget file - we will use it for sampling during bootstrapping
act=read.csv("proportions_foraging.csv",header=TRUE)

##Only consider data from days when more than 30 locations are available (to avoid underestimationof dispersal)
temp <- table(data$Date)
temp<-names(temp[temp[] > 30])
data <- subset(data,Date %in% temp)
##dectale dataframe to store the results
final = data.frame(0)

###We have to make use of nest locations (change this parameter according to the focal individual)
nest_lat = 1.1111
nest_long = 1.1111

## Remove the data points 
#new data frame 
#data1 <- data

data$nest = numeric(length=nrow(data))
data$roost =numeric(length=nrow(data))
k = 1
i = 1

#Find out nest locations and remove them
for(i in 1:nrow(data))
{
  lat = data[i,"Lat"]
  long = data[i,"Long"]
  dist_nest = distHaversine(c(long,lat), c(nest_long,nest_lat))
  if(dist_nest < 20)
    data$nest[i] = 1
  #dist_roost = distHaversine(c(long,lat), c(roost_long,roost_lat))
  #if(dist_roost < 50)
  #  data$roost[i] = 1
}
data1 <- data[data$nest != 1,]
if(nrow(data1) == 0)
{
  message <- "All points lie within nest area"
  message
}

##Create a vector to represent each row of data frame as many times as activity level in that hour
times=numeric(length=nrow(data1))
for(i in 1:nrow(data1)){
  times[i]=round(act$Nonbreeding[max(which((strptime(act$Time.of.the.day,format="%H:%M:%S") <= strptime(data1$Time[i],format="%H:%M:%S"))))])
  
}
times[is.na(times)] = 0
samp_vec = rep(1:nrow(data1),times=times)

##calculation for every time bin
#final$nest = 0
#final$roost = 0
##Loop over all gpts
for(i in 1:length(gpt$time))
{
  #Select first bin from GPT
  tim = gpt$bin[i]
  j = 1
  while(j <= 10)
  {
    #Sampling according to activity levels
    ran <- sample(samp_vec,1)
    ##Select source location and time
    latx = data1[ran,"Lat"]
    lonx = data1[ran,"Long"]  
    day = data1[ran,"Date"] 
    reftime = strptime(as.character(data1$Stamp[ran]),format = "%d/%m/%y %H:%M")
    reftime = format(reftime,"%H:%M")
    if(reftime>=format(strptime("17:00","%H:%M"),"%H:%M"))
      next
    t = 0
    temp = which(data$Location == data1$Location[ran])
    if(temp <= nrow(data))
      while(data[temp,"Date"]  == day)
      {
        temp = temp+1
        if(temp >= nrow(data))
          break
        
        dif = as.duration(strptime(data$Time[temp],format = "%H:%M") - strptime(reftime,format="%H:%M")) 
        dif = as.numeric(dif/60)
        if(!is.na(dif) & length(dif) != 0 & dif <= tim & dif > tim-1)
        {
          t = temp
          break
        }
        
      }
    if(t == 0)                        ##take last value if destination was not found within given interval
      t = temp-1
    
    laty = data[t,"Lat"]
    lony = data[t,"Long"]
    #dist = distHaversine(c(lony,laty), c(lonx,latx))
    
    dist = distHaversine(c(lony,laty), c(lonx,latx))
    dist_nest = distHaversine(c(lony,laty), c(nest_long,nest_lat))
    if(dist_nest < 20)
      final[k,"nest"] = dist
    
    
    #dist_roost = distHaversine(c(lony,laty), c(roost_long,roost_lat))
    #if(dist_roost < 50)
    #final[k,"roost"] = dist
    
    
    if(dist_nest >= 20)
      final[k,"other"] = dist
    
    final[k,"lat"] = latx
    final[k,"long"] = lonx
    #  final[k,"dist"] = dist
    # final[k,"time"] = gpt$time
    # final[k,"prop"] = gpt$Prop[u]
    #final[k,"flag"] = flag
    k = k + 1
    
    j=j+1
  }
  
}
write.csv(final$nest,"sdk_output_nest_tkbhai_17Sep18_whgh_10.csv")
write.csv(final$other,"sdk_output_other_tkbhai_17Sep18_whgh_10.csv")
final$other<-cut(final$other,breaks=c(0,20,150,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),dig.lab=4)
tiff(file = "Dispersal_prob_TKBHAI_17Sep18_whgh_10.tiff",width = 1420, height = 840,pointsize = 18)
x = barplot(table(final$other)/nrow(final),beside = FALSE,space = 2,legend = "Other locations",col = "gray",xlab = "Distance(m)",ylab = "Probability of seed arrival")
dev.off()


###chi-square breeding time to test if there are differences in visitation patterns in one hour intervals during the breeding season
dat<-read.csv("dat_breed_time_chi.csv",header=T)
dat
head(dat)
x<-table(dat$cat,dat$time)
x
chisq.test(x)


###chi-square non-breeding time. to test if there are differences in visitation patterns in one hour intervals during the non-breeding season
dat<-read.csv("dat_nonbreed_time_chi.csv",header=T)
dat
head(dat)
x<-table(dat$cat,dat$time)
x
chisq.test(x)

###chi-square breeding and non-breeding time. to test if there are differences in relative proportions of sightings in each hour interval between breeding and non-breeding season
dat<-read.csv("dat_breed_nonbreed_time_chi.csv",header=T)
dat
head(dat)
x<-table(dat$cat,dat$time)
x
chisq.test(x)


###chi-square to test for differences in relative proportion of seeds dispersed in nests vs. non-nest sites by Great Hornbills
dat<-read.csv("dat_breed_gh_nest_other_chi.csv",header=T)
dat
head(dat)
x<-table(dat$category)
x
chisq.test(x)

###chi-square to test for differences in relative proportion of seeds dispersed in nest vs non-nest sites by each of the three Great Hornbill breeding individuals
dat<-read.csv("dat_breed_gh_nest_other_chi.csv",header=T)
dat
head(dat)
x<-table(dat$code,dat$category)
x
chisq.test(x)


###chi-square to test for differences in relative proportion of seeds dispersed in nests vs. non-nest sites by Wreathed Hornbills
dat<-read.csv("dat_breed_wh_nest_other_chi.csv",header=T)
dat
head(dat)
x<-table(dat$category)
x
chisq.test(x)
