##Home Range Script

#Do this to prevent graphics error
dev.off()

####IF GET LOTS OF ERRORS - close project and re-open. 
#Dendrograms script uses plyr before dplyr and that messes up some of the dplyr commands in here. 

library(here)

####Load Data####

#Load in waypoints
GPS16 = read.csv(here::here("Input files","All Non-breeding GPS 2016.csv"))
GPS17 = read.csv(here::here("Input files","All Non-breeding GPS 2017.csv"))
GPS18 = read.csv(here::here("Input files","All Non-breeding GPS 2018.csv"))
GPS19 = read.csv(here::here("Input files","All Non-breeding GPS 2019.csv"))

#Load in bird dataframes
bird16 = read.csv(here::here("Input files","bird16.csv"))
bird17R = read.csv(here::here("Input files","bird17R.csv"))
bird18 = read.csv(here::here("Input files","bird18.csv"))
bird19 = read.csv(here::here("Input files","bird19.csv"))

#Load birdlists with birds in communities - birds in a group on their own have been removed
birdlist16noKsubc = read.csv(here::here("Output files","birdlist16noKsubc.csv"))
birdlist17noKsubc = read.csv(here::here("Output files","birdlist17noKsubc.csv"))
birdlist18noKsubc = read.csv(here::here("Output files","birdlist18noKsubc.csv"))
birdlist19noKsubc = read.csv(here::here("Output files","birdlist19noKsubc.csv"))


##Match waypoints to birds
#Can use (all.x = T) in the merge commands below to check that all waypoint numbers in sightings match up
#to an actual waypoint. In 2016 a couple JFW waypoints missing - not in waypoint file either. Must've 
#forgotten to take them. 
#In 2017 missing 20 SAD GPS points from her last two days in the season - must not have downloaded them

#2016 didn't include initials in waypoint
library(tidyr)
bird16.gps = bird16
bird16.gps=unite(bird16.gps, WP2, c(Initials, WP), remove=FALSE, sep="") #Match initials to waypoint
bird16.gps = merge(bird16.gps,GPS16,by.x="WP2",by.y="Waypoint")
#2017
bird17.gps = bird17R
bird17.gps = merge(bird17.gps,GPS17,by.x="WP",by.y="Waypoint")
#2018
bird18.gps = bird18
bird18.gps = merge(bird18.gps,GPS18,by.x="WP",by.y="Waypoint")
#2019
bird19.gps = bird19
bird19.gps = merge(bird19.gps,GPS19,by.x="WP",by.y="Waypoint")

#Get Date as date
bird16.gps$Date = as.Date(bird16.gps$Date,"%m/%d/%y")
bird17.gps$Date = as.Date(bird17.gps$Date,"%m/%d/%y")
bird18.gps$Date = as.Date(bird18.gps$Date,"%m/%d/%y")
bird19.gps$Date = as.Date(bird19.gps$Date,"%m/%d/%y")

#Get time as time
bird16.gps$Time = as.POSIXct(bird16.gps$Time,format = "%m/%d/%y %H:%M")
bird17.gps$Time = as.POSIXct(bird17.gps$Time,format = "%m/%d/%y %H:%M")
bird18.gps$Time = as.POSIXct(bird18.gps$Time,format = "%m/%d/%y %H:%M")
bird19.gps$Time = as.POSIXct(bird19.gps$Time,format = "%m/%d/%y %H:%M")

#Get time in hours and minutes
bird16.gps$Time2 <- format(bird16.gps$Time,"%H:%M")
bird17.gps$Time2 <- format(bird17.gps$Time,"%H:%M")
bird18.gps$Time2 <- format(bird18.gps$Time,"%H:%M")
bird19.gps$Time2 <- format(bird19.gps$Time,"%H:%M")

#Combine date and time into a timestamp
bird16.gps$timestamp = as.POSIXct(paste(bird16.gps$Date,bird16.gps$Time2), format="%Y-%m-%d %H:%M", tz="Australia/Brisbane")
bird17.gps$timestamp = as.POSIXct(paste(bird17.gps$Date,bird17.gps$Time2), format="%Y-%m-%d %H:%M", tz="Australia/Brisbane")
bird18.gps$timestamp = as.POSIXct(paste(bird18.gps$Date,bird18.gps$Time2), format="%Y-%m-%d %H:%M", tz="Australia/Brisbane")
bird19.gps$timestamp = as.POSIXct(paste(bird19.gps$Date,bird19.gps$Time2), format="%Y-%m-%d %H:%M", tz="Australia/Brisbane")


#Match GPS data to birds that were assigned groups and were not in a group on their own
bird16.gps = bird16.gps[bird16.gps$Bird %in% birdlist16noKsubc$Bird,]
bird16.gps = droplevels(bird16.gps)
bird17.gps = bird17.gps[bird17.gps$Bird %in% birdlist17noKsubc$Bird,]
bird17.gps = droplevels(bird17.gps)
bird18.gps = bird18.gps[bird18.gps$Bird %in% birdlist18noKsubc$Bird,]
bird18.gps = droplevels(bird18.gps)
bird19.gps = bird19.gps[bird19.gps$Bird %in% birdlist19noKsubc$Bird,]
bird19.gps = droplevels(bird19.gps)


#Calculate how many waypoints per bird
wpbird16 = data.frame(table(bird16.gps$Bird))
colnames(wpbird16) = c("Bird","WPfreq")
hist(wpbird16$WPfreq)
wpbird17 = data.frame(table(bird17.gps$Bird))
colnames(wpbird17) = c("Bird","WPfreq")
hist(wpbird17$WPfreq)
wpbird18 = data.frame(table(bird18.gps$Bird))
colnames(wpbird18) = c("Bird","WPfreq")
hist(wpbird18$WPfreq)
wpbird19 = data.frame(table(bird19.gps$Bird))
colnames(wpbird19) = c("Bird","WPfreq")
hist(wpbird19$WPfreq)

#Calculate how many observations per bird
library(dplyr)
obird16.gps = bird16.gps %>% count(Bird,Observation) %>% count(Bird)
hist(obird16.gps$n)
colnames(obird16.gps)[2] = "obs"
obird17.gps = bird17.gps %>% count(Bird,Observation) %>% count(Bird)
hist(obird17.gps$n)
colnames(obird17.gps)[2] = "obs"
obird18.gps = bird18.gps %>% count(Bird,Observation) %>% count(Bird)
hist(obird18.gps$n)
colnames(obird18.gps)[2] = "obs"
obird19.gps = bird19.gps %>% count(Bird,Observation) %>% count(Bird)
hist(obird19.gps$n)
colnames(obird19.gps)[2] = "obs"


#Are number of observations and number of waypoints correlated?
obird16.gps = merge(obird16.gps,wpbird16,by="Bird")
obird17.gps = merge(obird17.gps,wpbird17,by="Bird")
obird18.gps = merge(obird18.gps,wpbird18,by="Bird")
obird19.gps = merge(obird19.gps,wpbird19,by="Bird")

plot(obird16.gps$obs,obird16.gps$WPfreq,pch=19)
plot(obird17.gps$obs,obird17.gps$WPfreq,pch=19)
plot(obird18.gps$obs,obird18.gps$WPfreq,pch=19)
plot(obird19.gps$obs,obird19.gps$WPfreq,pch=19)
#Definitely, as expected and as should be! 



####Calculate AKDEs####

library(ctmm)

#Get bird.gps column names in Movebank format
colnames(bird16.gps)[c(9,11,12)] = c("individual.local.identifier","location.lat","location.long")
colnames(bird17.gps)[c(8,10,11)] = c("individual.local.identifier","location.lat","location.long")
colnames(bird18.gps)[c(8,10,11)] = c("individual.local.identifier","location.lat","location.long")
colnames(bird19.gps)[c(8,10,11)] = c("individual.local.identifier","location.lat","location.long")

#Order bird.gps by timestamp
bird16.gps = bird16.gps[order(bird16.gps$timestamp),]
bird17.gps = bird17.gps[order(bird17.gps$timestamp),]
bird18.gps = bird18.gps[order(bird18.gps$timestamp),]
bird19.gps = bird19.gps[order(bird19.gps$timestamp),]

#Get a blank column for akde area
birdlist16noKsubc$akde = NA
birdlist17noKsubc$akde = NA
birdlist18noKsubc$akde = NA
birdlist19noKsubc$akde = NA

#Get starting list
akde16.listb = list()
akde17.listb = list()
akde18.listb = list()
akde19.listb = list()

#Function to calculate akde for each individual and then put akde in a list
akde.rbfw = function(birdlist,bird.gps,UD,akde.list) {
  par(mfrow=c(5,5))
  for (i in 1:nrow(birdlist)) {
    #Get bird
    bird = birdlist$Bird[i]
    #Subset bird.gps to that bird
    bird.gps2 = bird.gps[which(bird.gps$individual.local.identifier==bird),]
    #Get into ctmm object
    HR <- as.telemetry(bird.gps2, timeformat="%Y-%m-%d $H:$M:$S", timezone="Australia/Brisbane",projection = sp::CRS("+init=epsg:32756"))
    #Get variogram
    var.HR = variogram(HR)
    #Plot variogram
    plot(var.HR,main=bird)
    #Get the fit parameters
    GUESS = variogram.fit(var.HR,interactive=F,name="GUESS")
    #Fit model
    mod.HR <- ctmm.select(HR, CTMM=GUESS, verbose=F)
    #Get kernel
    akde.HR<- akde(HR, CTMM=mod.HR)
    #Put akde estimate into birdlist
    birdlist$akde[i] = summary(akde.HR, level.UD=UD)$CI[2]
    #Make akde.HR into a spatial polygons dataframe
    akde.HR.sp = SpatialPolygonsDataFrame.UD(akde.HR,level.UD =0.85)[2,] #use just the estimate polygon
    #If it's the first row in the birdlist, then don't rbind it to the previous versions
    if(i==1) {akde.HR.sp2 = akde.HR.sp} else {akde.HR.sp2 = rbind(akde.HR.sp2,akde.HR.sp)}
  }
  par(mfrow=c(1,1))
  #Add birdlist to the list as well at the end of the list - only way to return multiple items
  akde.list[[1]] = birdlist
  akde.list[[2]] = akde.HR.sp2
  return(akde.list)
}

#Do this to prevent graphics error
#dev.off()
#par(mar=c(1,1,1,1))

##Get akde's for every year
#Check the plots - they shouldn't increase with time
#This output shows the minimum sampling interval for each individual - went through each of these that were lower than 5 minutes
#and made sure they were in separate observations. Left them in if they were in separate observations. 
# #2016
# akde16.list = akde.rbfw(birdlist=birdlist16noKsubc,bird.gps=bird16.gps,UD=0.85,akde.list = akde16.listb)
# birdlist16akde = akde16.list[[1]] #Birdlist is first in list
# akde16.sp = akde16.list[[2]] #Spatial points dataframe is second in list
# saveRDS(here::here("Output files/Homerange output files",akde16.list,"akde16_list.rds"))
# #2017
# akde17.list = akde.rbfw(birdlist=birdlist17noKsubc,bird.gps=bird17.gps,UD=0.85,akde.list = akde17.listb)
# birdlist17akde = akde17.list[[1]] #Birdlist is first in list
# akde17.sp = akde17.list[[2]] #Spatial points dataframe is second in list
# saveRDS(here::here("Output files/Homerange output files",akde17.list,"akde17_list.rds"))
# #2018
# akde18.list = akde.rbfw(birdlist=birdlist18noKsubc,bird.gps=bird18.gps,UD=0.85,akde.list = akde18.listb)
# birdlist18akde = akde18.list[[1]] #Birdlist is first in list
# akde18.sp = akde18.list[[2]] #Spatial points dataframe is second in list
# saveRDS(here::here("Output files/Homerange output files",akde18.list,"akde18_list.rds"))
# #2019
# akde19.list = akde.rbfw(birdlist=birdlist19noKsubc,bird.gps=bird19.gps,UD=0.85,akde.list = akde19.listb)
# birdlist19akde = akde19.list[[1]] #Birdlist is first in list
# akde19.sp = akde19.list[[2]] #Spatial points dataframe is second in list
# saveRDS(here::here("Output files/Homerange output files",akde19.list,"akde19_list.rds"))


####Just plot variograms####
variogram.rbfw = function(birdlist,bird.gps) {
  par(mfrow=c(5,5))
  for (i in 1:nrow(birdlist)) {
    #Get bird
    bird = birdlist$Bird[i]
    #Subset bird.gps to that bird
    bird.gps2 = bird.gps[which(bird.gps$individual.local.identifier==bird),]
    #Get into ctmm object
    HR <- as.telemetry(bird.gps2, timeformat="%Y-%m-%d $H:$M:$S", timezone="Australia/Brisbane",projection = sp::CRS("+init=epsg:32756"))
    #Get variogram
    var.HR = variogram(HR)
    #Plot variogram
    plot(var.HR,main=bird)
  }
  par(mfrow=c(1,1))
}


#Plot
#variogram.rbfw(birdlist=birdlist16noKsubc,bird.gps=bird16.gps)
# #Look at YHZ, RIG, RYH
# variogram.rbfw(birdlist=birdlist17noKsubc,bird.gps=bird17.gps)
# #Look at VLL and IRW variograms for 2017 - may have moved homeranges
# variogram.rbfw(birdlist=birdlist18noKsubc,bird.gps=bird18.gps)
# #Look at GGG, LRB, LZI, RVR, RVY, LGY, VWW, ZZB, HYH, ILB
# variogram.rbfw(birdlist=birdlist19noKsubc,bird.gps=bird19.gps)
# #Look at BLB, BZR, LHI, RHY, RLR, VLR, WWH, WYY, YHR, ZYB



####Run a single individual through to AKDE####
#This section shows the basics of how the ctmm code works

#Get a test bird - can only run the model for one bird at a time
ind = bird17.gps[which(bird17.gps$individual.local.identifier=="GRG"),]#Need to order by timestamp
ind = ind[order(ind$timestamp),]

#Get time difference between timestamp
ind$timediff = NA
for (i in 2:nrow(ind)) {
  ind$timediff[i] = ind$timestamp[i]-ind$timestamp[i-1]
}

#Get into ctmm format - can set a projection here if you want. For KDEs though best not to use an ESPG
#projection, object gets projected for each individual, centered on it's territory. 
#projection = sp::CRS("+init=epsg:32756") - if you want a projection for plotting purposes
HR <- as.telemetry(ind, timeformat="%Y-%m-%d $H:$M:$S", timezone="Australia/Brisbane",projection = sp::CRS("+init=epsg:32756"))

#Plot points
plot(HR)

#Get variogram and plot the fit -"The variogram represents the average square distance traveled 
#(vertical axis) within some time lag (horizontal axis)." - from ctmm vingette.
var.HR <- variogram(HR,fast=F,dt=10000)
variogram.fit(var.HR)

#Get the fit parameters
GUESS = variogram.fit(var.HR,interactive=F,name="GUESS")

#Fit multiple models
mod.HR <- ctmm.select(HR, CTMM=GUESS, verbose=TRUE)
summary(mod.HR)

#Get the best model
iid.HR <- mod.HR[[1]]
summary(iid.HR)

#Calculate autocorrelated Kernel Density Estimate
akde.HR<- akde(HR, CTMM=iid.HR)
summary(akde.HR)
plot(HR, UD=akde.HR, level.UD=0.85)

#Get akde at a specific confidence level - 85% seems pretty good to keep territories from going over water
#but need to check this with more birds first. 
summary(akde.HR, level.UD=0.85)


####Read in AKDE lists####
akde16.list = readRDS(here::here("Output files/Homerange output files","akde16_list.rds"))
akde17.list = readRDS(here::here("Output files/Homerange output files","akde17_list.rds"))
akde18.list = readRDS(here::here("Output files/Homerange output files","akde18_list.rds"))
akde19.list = readRDS(here::here("Output files/Homerange output files","akde19_list.rds"))

#Get birdlists with adke
birdlist16akde = akde16.list[[1]] #Birdlist is first in list
birdlist17akde = akde17.list[[1]] #Birdlist is first in list
birdlist18akde = akde18.list[[1]] #Birdlist is first in list
birdlist19akde = akde19.list[[1]] #Birdlist is first in list


##Birds to remove - variograms and waypoint plots are weird - these birds clearly moved home ranges or do 
#not have enough points. 
#Removing YLI here because her akde range goes way out over the water which seems unlikely. Can take 
#her out of the bad17 list to see what her homerange looks like. 
#2016: None
#2017: 
bad17 = c("VLL", "IRW", "YLI")
#2018: 
bad18 = c("GGG", "LRB", "RVY", "LGY", "ZZB", "HYH" )
#2019: 
bad19 = c("LHI", "ZYB")

#Get birdlists without the bad birds
birdlist17akde = birdlist17akde[!birdlist17akde$Bird %in% bad17,]
birdlist18akde = birdlist18akde[!birdlist18akde$Bird %in% bad18,]
birdlist19akde = birdlist19akde[!birdlist19akde$Bird %in% bad19,]

#Merge birdlist and obird.gps to get number of observations and waypoints per bird
birdlist16akde = merge(birdlist16akde,obird16.gps,by="Bird")
birdlist17akde = merge(birdlist17akde,obird17.gps,by="Bird")
birdlist18akde = merge(birdlist18akde,obird18.gps,by="Bird")
birdlist19akde = merge(birdlist19akde,obird19.gps,by="Bird")


#Get homerange size/akde into hectares
for (i in 1:nrow(birdlist16akde)) {
  if(birdlist16akde$akde[i]>1000) {birdlist16akde$akde[i]=(birdlist16akde$akde[i]/10000)}
}
for (i in 1:nrow(birdlist17akde)) {
  if(birdlist17akde$akde[i]>1000) {birdlist17akde$akde[i]=(birdlist17akde$akde[i]/10000)}
}
for (i in 1:nrow(birdlist18akde)) {
  if(birdlist18akde$akde[i]>1000) {birdlist18akde$akde[i]=(birdlist18akde$akde[i]/10000)}
}
for (i in 1:nrow(birdlist19akde)) {
  if(birdlist19akde$akde[i]>1000) {birdlist19akde$akde[i]=(birdlist19akde$akde[i]/10000)}
}

#Does number of observations relate to home range size? 
plot(birdlist16akde$obs,birdlist16akde$akde,pch=19)
cor.test(birdlist16akde$obs,birdlist16akde$akde)
plot(birdlist17akde$obs,birdlist17akde$akde,pch=19)
cor.test(birdlist17akde$obs,birdlist17akde$akde)
plot(birdlist18akde$obs,birdlist18akde$akde,pch=19)
cor.test(birdlist18akde$obs,birdlist18akde$akde)
plot(birdlist19akde$obs,birdlist19akde$akde,pch=19)
cor.test(birdlist19akde$obs,birdlist19akde$akde)
#Maybe? Could cut off birds that have fewer than 10 observations in 2018


####How many observations needed to calculate homerange size?####

#Think that home range size decreases as the season goes on, but before I can break the season up into segments
#I need to figure out how many observations are required to get a good estimate of home range size. Using
#observations because sightings within observations are autocorrelated. 

#Try looking at 5, 10, 15 points at first

#Get birds in each year that were in at least 15 observations
birdlist16akde15 = birdlist16akde[which(birdlist16akde$obs>=15),]
birdlist17akde15 = birdlist17akde[which(birdlist17akde$obs>=15),]
birdlist18akde15 = birdlist18akde[which(birdlist18akde$obs>=15),]
birdlist19akde15 = birdlist19akde[which(birdlist19akde$obs>=15),]

###Function for simulations
#Sample.obs is the number of observations to use for each individual, runs is the number of randomizations/simulations 
akde.rbfw.sim = function(birdlist,bird.gps,UD,sample.obs,runs) {
  sim.data = data.frame(matrix(NA, nrow = runs, ncol = 1))
  colnames(sim.data) = "mean"
  sim.list = list()
  for (h in 1:nrow(sim.data)) {
  birdlist$akde.rand = NA
  for (i in 1:nrow(birdlist)) try({
    #Get bird
    bird = birdlist$Bird[i]
    #Subset bird.gps to that bird
    bird.gps2 = bird.gps[which(bird.gps$individual.local.identifier==bird),]
    #Get list of observations
    bird.obs = bird.gps2$Observation
    bird.obs = bird.obs[!duplicated(bird.obs)]
    #Sample observations 
    sample.list = sample(bird.obs,sample.obs,replace = F)
    #Get only those observations
    bird.gps3 = bird.gps2[bird.gps2$Observation %in% sample.list,]
    #Get into ctmm object
    HR <- as.telemetry(bird.gps3, timeformat="%Y-%m-%d $H:$M:$S", timezone="Australia/Brisbane",projection = sp::CRS("+init=epsg:32756"))
    #Get variogram
    var.HR = variogram(HR)
    #Get the fit parameters
    GUESS = variogram.fit(var.HR,interactive=F,name="GUESS")
    #Fit model
    mod.HR <- ctmm.select(HR, CTMM=GUESS, verbose=F)
    #Get kernel
    akde.HR<- akde(HR, CTMM=mod.HR)
    #Put akde estimate into birdlist
    birdlist$akde.rand[i] = summary(akde.HR, level.UD=UD)$CI[2]
    #Convert to hectares if in meters
    if(birdlist$akde.rand[i]>1000) {birdlist$akde.rand[i] = (birdlist$akde.rand[i]/10000)}})
  sim.list[[h]]=birdlist
  sim.data$mean[h] = mean(birdlist$akde.rand,na.rm=T)
  print(c("round",h))} 
  sim.list[[h+1]] = data.frame(sim.data)
  return(sim.list)
  }


#Run it at 5, 10, 15 obs per bird
#2016
# akde16.rand.5obs = akde.rbfw.sim(birdlist=birdlist16akde15,bird.gps=bird16.gps,UD=0.85,sample.obs=5,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde16.rand.5obs,"akde16.rand.5obs.rds"))
# akde16.rand.7obs = akde.rbfw.sim(birdlist=birdlist16akde15,bird.gps=bird16.gps,UD=0.85,sample.obs=7,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde16.rand.7obs,"akde16.rand.7obs.rds"))
# akde16.rand.10obs = akde.rbfw.sim(birdlist=birdlist16akde15,bird.gps=bird16.gps,UD=0.85,sample.obs=10,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde16.rand.10obs,"akde16.rand.10obs.rds"))
#2017
# akde17.rand.5obs = akde.rbfw.sim(birdlist=birdlist17akde15,bird.gps=bird17.gps,UD=0.85,sample.obs=5,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde17.rand.5obs,"akde17.rand.5obs.rds"))
# akde17.rand.7obs = akde.rbfw.sim(birdlist=birdlist17akde15,bird.gps=bird17.gps,UD=0.85,sample.obs=7,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde17.rand.7obs,"akde17.rand.7obs.rds"))
# akde17.rand.10obs = akde.rbfw.sim(birdlist=birdlist17akde15,bird.gps=bird17.gps,UD=0.85,sample.obs=10,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde17.rand.10obs,"akde17.rand.10obs.rds"))
# akde17.rand.15obs = akde.rbfw.sim(birdlist=birdlist17akde15,bird.gps=bird17.gps,UD=0.85,sample.obs=15,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde17.rand.15obs,"akde17.rand.15obs.rds"))
#2018
# akde18.rand.5obs = akde.rbfw.sim(birdlist=birdlist18akde15,bird.gps=bird18.gps,UD=0.85,sample.obs=5,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde18.rand.5obs,"akde18.rand.5obs.rds"))
# akde18.rand.7obs = akde.rbfw.sim(birdlist=birdlist18akde15,bird.gps=bird18.gps,UD=0.85,sample.obs=7,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde18.rand.7obs,"akde18.rand.7obs.rds"))
# akde18.rand.10obs = akde.rbfw.sim(birdlist=birdlist18akde15,bird.gps=bird18.gps,UD=0.85,sample.obs=10,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde18.rand.10obs,"akde18.rand.10obs.rds"))
#2019
# akde19.rand.5obs = akde.rbfw.sim(birdlist=birdlist19akde15,bird.gps=bird19.gps,UD=0.85,sample.obs=5,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde19.rand.5obs,"akde19.rand.5obs.rds"))
# akde19.rand.7obs = akde.rbfw.sim(birdlist=birdlist19akde15,bird.gps=bird19.gps,UD=0.85,sample.obs=7,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde19.rand.7obs,"akde19.rand.7obs.rds"))
# akde19.rand.10obs = akde.rbfw.sim(birdlist=birdlist19akde15,bird.gps=bird19.gps,UD=0.85,sample.obs=10,runs=10)
# saveRDS(here::here("Output files/Homerange output files",akde19.rand.10obs,"akde19.rand.10obs.rds"))


#Read in saved randomization files 
akde16.rand.5obs = readRDS(here::here("Output files/Homerange output files","akde16.rand.5obs.rds"))
akde16.rand.7obs = readRDS(here::here("Output files/Homerange output files","akde16.rand.7obs.rds"))
akde16.rand.10obs = readRDS(here::here("Output files/Homerange output files","akde16.rand.10obs.rds"))
akde17.rand.5obs = readRDS(here::here("Output files/Homerange output files","akde17.rand.5obs.rds"))
akde17.rand.7obs = readRDS(here::here("Output files/Homerange output files","akde17.rand.7obs.rds"))
akde17.rand.10obs = readRDS(here::here("Output files/Homerange output files","akde17.rand.10obs.rds"))
akde17.rand.15obs = readRDS(here::here("Output files/Homerange output files","akde17.rand.15obs.rds"))
akde18.rand.5obs = readRDS(here::here("Output files/Homerange output files","akde18.rand.5obs.rds"))
akde18.rand.7obs = readRDS(here::here("Output files/Homerange output files","akde18.rand.7obs.rds"))
akde18.rand.10obs = readRDS(here::here("Output files/Homerange output files","akde18.rand.10obs.rds"))
akde19.rand.5obs = readRDS(here::here("Output files/Homerange output files","akde19.rand.5obs.rds"))
akde19.rand.7obs = readRDS(here::here("Output files/Homerange output files","akde19.rand.7obs.rds"))
akde19.rand.10obs = readRDS(here::here("Output files/Homerange output files","akde19.rand.10obs.rds"))


###Histogram comparisons 
#These histograms show that 5 observations do a decent job of getting the home-range size, so using periods 
#where birds have at least 5 observations as periods for comparisons
library(ggplot2)
library(cowplot)
#2016
akde16.rand.5obs.means = akde16.rand.5obs[[length(akde16.rand.5obs)]]
akde16.rand.7obs.means = akde16.rand.7obs[[length(akde16.rand.7obs)]]
akde16.rand.10obs.means = akde16.rand.10obs[[length(akde16.rand.10obs)]]
akde16.rand.all = cbind(akde16.rand.5obs.means,akde16.rand.7obs.means,akde16.rand.10obs.means)
colnames(akde16.rand.all) = c("means5","means7","means10")
ggplot(data=akde16.rand.all) +  
  geom_density(aes(x=means5),fill="light blue",alpha=0.7) + 
  geom_density(aes(x=means7),fill="dark red",alpha=0.7) + 
  geom_density(aes(x=means10),fill="yellow", alpha=0.7) + theme_cowplot() +
  geom_vline(aes(xintercept=mean(birdlist16akde15$akde)),size=1.5,linetype=2) + xlab("Mean aKDE") +
  ggtitle("2016")
#2017
akde17.rand.5obs.means = akde17.rand.5obs[[length(akde17.rand.5obs)]]
akde17.rand.7obs.means = akde17.rand.7obs[[length(akde17.rand.7obs)]]
akde17.rand.10obs.means = akde17.rand.10obs[[length(akde17.rand.10obs)]]
akde17.rand.15obs.means = akde17.rand.15obs[[length(akde17.rand.15obs)]]
akde17.rand.all = cbind(akde17.rand.5obs.means,akde17.rand.7obs.means,akde17.rand.10obs.means,akde17.rand.15obs.means)
colnames(akde17.rand.all) = c("means5","means7","means10","means15")
ggplot(data=akde17.rand.all) +  
  geom_density(aes(x=means5),fill="light blue",alpha=0.7) + 
  geom_density(aes(x=means7),fill="dark red",alpha=0.7) + 
  geom_density(aes(x=means10),fill="yellow", alpha=0.7) + theme_cowplot() +
  geom_density(aes(x=means15),fill="pink", alpha=0.7) +
  geom_vline(aes(xintercept=mean(birdlist17akde15$akde)),size=1.5,linetype=2) + xlab("Mean aKDE") +
  ggtitle("2017")
#2018
akde18.rand.5obs.means = akde18.rand.5obs[[length(akde18.rand.5obs)]]
akde18.rand.7obs.means = akde18.rand.7obs[[length(akde18.rand.7obs)]]
akde18.rand.10obs.means = akde18.rand.10obs[[length(akde18.rand.10obs)]]
akde18.rand.all = cbind(akde18.rand.5obs.means,akde18.rand.7obs.means,akde18.rand.10obs.means)
colnames(akde18.rand.all) = c("means5","means7","means10")
ggplot(data=akde18.rand.all) +  
  geom_density(aes(x=means5),fill="light blue",alpha=0.7) + 
  geom_density(aes(x=means7),fill="dark red",alpha=0.7) + 
  geom_density(aes(x=means10),fill="yellow", alpha=0.7) + theme_cowplot() +
  geom_vline(aes(xintercept=mean(birdlist18akde15$akde)),size=1.5,linetype=2) + xlab("Mean aKDE") +
  ggtitle("2018")
#2019
akde19.rand.5obs.means = akde19.rand.5obs[[length(akde19.rand.5obs)]]
akde19.rand.7obs.means = akde19.rand.7obs[[length(akde19.rand.7obs)]]
akde19.rand.10obs.means = akde19.rand.10obs[[length(akde19.rand.10obs)]]
akde19.rand.all = cbind(akde19.rand.5obs.means,akde19.rand.7obs.means,akde19.rand.10obs.means)
colnames(akde19.rand.all) = c("means5","means7","means10")
ggplot(data=akde19.rand.all) +  
  geom_density(aes(x=means5),fill="light blue",alpha=0.7) + 
  geom_density(aes(x=means7),fill="dark red",alpha=0.7) + 
  geom_density(aes(x=means10),fill="yellow", alpha=0.7) + theme_cowplot() +
  geom_vline(aes(xintercept=mean(birdlist19akde15$akde)),size=1.5,linetype=2) + xlab("Mean aKDE") +
  ggtitle("2019")


###Plot simulations with model means
#Testing each year separately since slightly different time periods are needed to maximize the number of 
#individuals per period
#This is similar to randomized means comparisons above, but here using the model to get the mean. 
#Model accounts for social class and matches the model that will be used to calculate the mean homerange size
#for each time period for each year

#Example model:
library(lme4)
example.model = lmer(akde~1 + (1|class), data=birdlist16akde15)
summary(example.model)$coefficients[1] #coefficient of model is mean but now class is taken into account
mean(birdlist16akde15$akde) #slightly different from mean without class accounted for

#Get observed model means for each year
observed16 = summary(lmer(akde~1 + (1|class), data=birdlist16akde15))$coefficients[1]
observed17 = summary(lmer(akde~1 + (1|class), data=birdlist17akde15))$coefficients[1]
observed18 = summary(lmer(akde~1 + (1|class), data=birdlist18akde15))$coefficients[1]
observed19 = summary(lmer(akde~1 + (1|class), data=birdlist19akde15))$coefficients[1]

#Drop previous means from akde.rand list files
akde16.rand.5obs.r = akde16.rand.5obs[-11]
akde17.rand.5obs.r = akde17.rand.5obs[-11]
akde18.rand.5obs.r = akde18.rand.5obs[-11]
akde19.rand.5obs.r = akde19.rand.5obs[-11]

test = akde16.rand.5obs.r[[1]]

model.rand = function(x) {
  lm1 = lmer(akde.rand~1 + (1|class),data=x)
  lm2 = summary(lm1)$coefficients[1]
  return(lm2)
}

#Get means using models
akde16.rand.5obs.r.means = data.frame(unlist(lapply(akde16.rand.5obs.r,model.rand)))
colnames(akde16.rand.5obs.r.means) = "means"
akde17.rand.5obs.r.means = data.frame(unlist(lapply(akde17.rand.5obs.r,model.rand)))
colnames(akde17.rand.5obs.r.means) = "means"
akde18.rand.5obs.r.means = data.frame(unlist(lapply(akde18.rand.5obs.r,model.rand)))
colnames(akde18.rand.5obs.r.means) = "means"
akde19.rand.5obs.r.means = data.frame(unlist(lapply(akde19.rand.5obs.r,model.rand)))
colnames(akde19.rand.5obs.r.means) = "means"

#Plot histograms
ggplot(data=akde16.rand.5obs.r.means,aes(x=means)) + geom_density(fill="light blue",alpha=0.7) +
  theme_cowplot() + geom_vline(aes(xintercept=observed16),size=1.5,linetype=2) + xlab("Mean aKDE") +
  ggtitle("2016")
ggplot(data=akde17.rand.5obs.r.means,aes(x=means)) + geom_density(fill="light blue",alpha=0.7) +
  theme_cowplot() + geom_vline(aes(xintercept=observed17),size=1.5,linetype=2) + xlab("Mean aKDE") +
  ggtitle("2017")
ggplot(data=akde18.rand.5obs.r.means,aes(x=means)) + geom_density(fill="light blue",alpha=0.7) +
  theme_cowplot() + geom_vline(aes(xintercept=observed18),size=1.5,linetype=2) + xlab("Mean aKDE") +
  ggtitle("2018")
ggplot(data=akde19.rand.5obs.r.means,aes(x=means)) + geom_density(fill="light blue",alpha=0.7) +
  theme_cowplot() + geom_vline(aes(xintercept=observed19),size=1.5,linetype=2) + xlab("Mean aKDE") +
  ggtitle("2019")

#Result: results are similar to simulation results above, 5 observations does a decent job of 
#getting near the observed mean. But since there is some variation will want to compare to randomizations
#in final tests among periods to make sure that observed difference seen is more than expected by random



####Figure out where to cut season into comparison periods####

##First need to figure out when observation dates were
#Subset bird.gps to birds with enough points - at least 15
bird16.gps15 = bird16.gps[bird16.gps$individual.local.identifier %in% birdlist16akde15$Bird,]
bird17.gps15 = bird17.gps[bird17.gps$individual.local.identifier %in% birdlist17akde15$Bird,]
bird18.gps15 = bird18.gps[bird18.gps$individual.local.identifier %in% birdlist18akde15$Bird,]
bird19.gps15 = bird19.gps[bird19.gps$individual.local.identifier %in% birdlist19akde15$Bird,]

#Then get a list of observations for each bird and add a "n" column to sum in the next step
bird16.gps15.list = bird16.gps15 %>% distinct(individual.local.identifier,Observation,n=1,.keep_all=T)
bird17.gps15.list = bird17.gps15 %>% distinct(individual.local.identifier,Observation,n=1,.keep_all=T)
bird18.gps15.list = bird18.gps15 %>% distinct(individual.local.identifier,Observation,n=1,.keep_all=T)
bird19.gps15.list = bird19.gps15 %>% distinct(individual.local.identifier,Observation,n=1,.keep_all=T)

#Add a count column
bird16.gps15.list = bird16.gps15.list %>% group_by(individual.local.identifier) %>% mutate(count = cumsum(n))
bird17.gps15.list = bird17.gps15.list %>% group_by(individual.local.identifier) %>% mutate(count = cumsum(n))
bird18.gps15.list = bird18.gps15.list %>% group_by(individual.local.identifier) %>% mutate(count = cumsum(n))
bird19.gps15.list = bird19.gps15.list %>% group_by(individual.local.identifier) %>% mutate(count = cumsum(n))

#Get julian dates 
library(lubridate)
bird16.gps15.list$jdate = yday(bird16.gps15.list$Date)
bird17.gps15.list$jdate = yday(bird17.gps15.list$Date)
bird18.gps15.list$jdate = yday(bird18.gps15.list$Date)
bird19.gps15.list$jdate = yday(bird19.gps15.list$Date)



###Look at when 5th, 10th, 15th observations occur
#2016
ggplot() + 
  geom_histogram(data=bird16.gps15.list[which(bird16.gps15.list$count==5),], aes(x=jdate),binwidth = 1,
                 fill="dark red",alpha=0.8) + theme_cowplot() + xlim(170,230) +
  geom_histogram(data=bird16.gps15.list[which(bird16.gps15.list$count==10),], aes(x=jdate),binwidth = 1,
                 fill="light blue",alpha=0.8) +
  geom_histogram(data=bird16.gps15.list[which(bird16.gps15.list$count==15),], aes(x=jdate),binwidth = 1,
                 fill="orange",alpha=0.8) + ggtitle("2016") +
  scale_x_continuous(breaks=seq(170,250,1)) + theme(axis.text.x = element_text(angle = 90)) + 
  geom_vline(aes(xintercept=192)) + geom_vline(aes(xintercept=205))

#2017
ggplot() + 
  geom_histogram(data=bird17.gps15.list[which(bird17.gps15.list$count==10),], aes(x=jdate),binwidth = 1,
                 fill="dark red",alpha=0.8) + theme_cowplot() + xlim(170,230) +
  geom_histogram(data=bird17.gps15.list[which(bird17.gps15.list$count==20),], aes(x=jdate),binwidth = 1,
                 fill="light blue",alpha=0.8) +
  geom_histogram(data=bird17.gps15.list[which(bird17.gps15.list$count==30),], aes(x=jdate),binwidth = 1,
                 fill="orange",alpha=0.8) + ggtitle("2017") +
  scale_x_continuous(breaks=seq(170,250,1)) + theme(axis.text.x = element_text(angle = 90)) + 
  geom_vline(aes(xintercept=195)) + geom_vline(aes(xintercept=211))
#2018
ggplot() + 
  geom_histogram(data=bird18.gps15.list[which(bird18.gps15.list$count==5),], aes(x=jdate),binwidth = 1,
                 fill="dark red",alpha=0.8) + theme_cowplot() + xlim(170,230) +
  geom_histogram(data=bird18.gps15.list[which(bird18.gps15.list$count==10),], aes(x=jdate),binwidth = 1,
                 fill="light blue",alpha=0.8) +
  geom_histogram(data=bird18.gps15.list[which(bird18.gps15.list$count==15),], aes(x=jdate),binwidth = 1,
                 fill="orange",alpha=0.8) + ggtitle("2018") +
  scale_x_continuous(breaks=seq(170,250,1)) + theme(axis.text.x = element_text(angle = 90)) + 
  geom_vline(aes(xintercept=201)) + geom_vline(aes(xintercept=210))
#2019
ggplot() + 
  geom_histogram(data=bird19.gps15.list[which(bird19.gps15.list$count==5),], aes(x=jdate),binwidth = 1,
                 fill="dark red",alpha=0.8) + theme_cowplot() + xlim(170,230) +
  geom_histogram(data=bird19.gps15.list[which(bird19.gps15.list$count==10),], aes(x=jdate),binwidth = 1,
                 fill="light blue",alpha=0.8) +
  geom_histogram(data=bird19.gps15.list[which(bird19.gps15.list$count==15),], aes(x=jdate),binwidth = 1,
                 fill="orange",alpha=0.8) + ggtitle("2019") +
  scale_x_continuous(breaks=seq(170,250,1)) + theme(axis.text.x = element_text(angle = 90)) + 
  geom_vline(aes(xintercept=197)) + geom_vline(aes(xintercept=210))


##Break up seasons into three time periods
#2016
bird16.gps15.list.1 = bird16.gps15.list %>% filter(jdate<192)
bird16.gps15.list.1L = bird16.gps15.list.1 %>% count(individual.local.identifier) %>% filter(n>=5)
bird16.gps15.list.2 = bird16.gps15.list %>% filter(jdate>=192,jdate<205)
bird16.gps15.list.2L = bird16.gps15.list.2 %>% count(individual.local.identifier) %>% filter(n>=5)
bird16.gps15.list.3 = bird16.gps15.list %>% filter(jdate>=205)
bird16.gps15.list.3L = bird16.gps15.list.3 %>% count(individual.local.identifier) %>% filter(n>=5)
#2017
bird17.gps15.list.1 = bird17.gps15.list %>% filter(jdate<195)
bird17.gps15.list.1L = bird17.gps15.list.1 %>% count(individual.local.identifier) %>% filter(n>=10)
bird17.gps15.list.2 = bird17.gps15.list %>% filter(jdate>=195,jdate<211)
bird17.gps15.list.2L = bird17.gps15.list.2 %>% count(individual.local.identifier) %>% filter(n>=10)
bird17.gps15.list.3 = bird17.gps15.list %>% filter(jdate>=211)
bird17.gps15.list.3L = bird17.gps15.list.3 %>% count(individual.local.identifier) %>% filter(n>=10)
#2018
bird18.gps15.list.1 = bird18.gps15.list %>% filter(jdate<201)
bird18.gps15.list.1L = bird18.gps15.list.1 %>% count(individual.local.identifier) %>% filter(n>=5)
bird18.gps15.list.2 = bird18.gps15.list %>% filter(jdate>=201,jdate<210)
bird18.gps15.list.2L = bird18.gps15.list.2 %>% count(individual.local.identifier) %>% filter(n>=5)
bird18.gps15.list.3 = bird18.gps15.list %>% filter(jdate>=210)
bird18.gps15.list.3L = bird18.gps15.list.3 %>% count(individual.local.identifier) %>% filter(n>=5)
#2019
bird19.gps15.list.1 = bird19.gps15.list %>% filter(jdate<197)
bird19.gps15.list.1L = bird19.gps15.list.1 %>% count(individual.local.identifier) %>% filter(n>=5)
bird19.gps15.list.2 = bird19.gps15.list %>% filter(jdate>=197,jdate<210)
bird19.gps15.list.2L = bird19.gps15.list.2 %>% count(individual.local.identifier) %>% filter(n>=5)
bird19.gps15.list.3 = bird19.gps15.list %>% filter(jdate>=210)
bird19.gps15.list.3L = bird19.gps15.list.3 %>% count(individual.local.identifier) %>% filter(n>=5)



####Plot home ranges and write shapefiles####

#This code is in this location because I caught a bird that has a really weird homerange by plotting
#see details below but taking it out of following analyses. 
#Shapefiles for Google Earth Engine or any other plotting of ranges

#Get spatial points dataframes out of akde.lists
sp16 = akde16.list[[2]]
sp17 = akde17.list[[2]]
sp18 = akde18.list[[2]]
sp19 = akde19.list[[2]]

#Get list of social groups
sg16 = akde16.list[[1]]$Social.Group
sg17 = akde17.list[[1]]$Social.Group
sg18 = akde18.list[[1]]$Social.Group
sg19 = akde19.list[[1]]$Social.Group

##Grab birds that have over 15 obs - birds used in homerange change across season analysis
#First add a column with #% est
splist16 = birdlist16akde15 %>% mutate(extra="85% est") %>% mutate(spformat=paste(Bird,extra))
splist17 = birdlist17akde15 %>% mutate(extra="85% est") %>% mutate(spformat=paste(Bird,extra))
splist18 = birdlist18akde15 %>% mutate(extra="85% est") %>% mutate(spformat=paste(Bird,extra))
splist19 = birdlist19akde15 %>% mutate(extra="85% est") %>% mutate(spformat=paste(Bird,extra))

#Then subset spatial polygons dataframe
sp16s = sp16[sp16@data$name %in% splist16$spformat,]
sp17s = sp17[sp17@data$name %in% splist17$spformat,]
sp18s = sp18[sp18@data$name %in% splist18$spformat,]
sp19s = sp19[sp19@data$name %in% splist19$spformat,]

#Plot spatial points basic
sp::plot(sp16s)
sp::plot(sp17s)
sp::plot(sp18s)
sp::plot(sp19s)

#Get datapoints of spatial polygon dataframes 
dpts16 <- ggplot2::fortify(sp16s)
dpts17 <- ggplot2::fortify(sp17s)
dpts18 <- ggplot2::fortify(sp18s)
dpts19 <- ggplot2::fortify(sp19s)

#Get bird ID separated
dpts16 = dpts16 %>% separate(id,into=c("Bird","2","3"))
dpts17 = dpts17 %>% separate(id,into=c("Bird","2","3"))
dpts18 = dpts18 %>% separate(id,into=c("Bird","2","3"))
dpts19 = dpts19 %>% separate(id,into=c("Bird","2","3"))

#Get groups
groups16 = birdlist16akde %>% select(Bird,Social.Group)
groups17 = birdlist17akde %>% select(Bird,Social.Group)
groups18 = birdlist18akde %>% select(Bird,Social.Group)
groups19 = birdlist19akde %>% select(Bird,Social.Group)

#Add group to datapoints dataframe
dpts16 = merge(dpts16,groups16,by="Bird")
dpts16$Social.Group = as.factor(dpts16$Social.Group)
dpts17 = merge(dpts17,groups17,by="Bird")
dpts17$Social.Group = as.factor(dpts17$Social.Group)
dpts18 = merge(dpts18,groups18,by="Bird")
dpts18$Social.Group = as.factor(dpts18$Social.Group)
dpts19 = merge(dpts19,groups19,by="Bird")
dpts19$Social.Group = as.factor(dpts19$Social.Group)

#Transform to lat long for plotting
library(rgdal)
sp16s2 = spTransform(sp16s, CRS("+proj=longlat"))
sp17s2 = spTransform(sp17s, CRS("+proj=longlat"))
sp18s2 = spTransform(sp18s, CRS("+proj=longlat"))
sp19s2 = spTransform(sp19s, CRS("+proj=longlat"))

#Plot
ggplot() + geom_polygon(data=sp16s2, aes(x=long,y=lat,group=id,color=dpts16$Social.Group),fill=NA) + 
  theme_cowplot() + scale_fill_identity() + xlab("Longitude") + ylab("Latitude") + labs(color="Social Group") + ggtitle("2016")
ggplot() + geom_polygon(data=sp17s2, aes(x=long,y=lat,group=id,color=dpts17$Social.Group),fill=NA) + 
  theme_cowplot() + scale_fill_identity() + xlab("Longitude") + ylab("Latitude") + labs(color="Social Group") + ggtitle("2017")
ggplot() + geom_polygon(data=sp18s2, aes(x=long,y=lat,group=id,color=dpts18$Social.Group),fill=NA) + 
  theme_cowplot() + scale_fill_identity() + xlab("Longitude") + ylab("Latitude") + labs(color="Social Group") + ggtitle("2018")
ggplot() + geom_polygon(data=sp19s2, aes(x=long,y=lat,group=id,color=dpts19$Social.Group),fill=NA) + 
  theme_cowplot() + scale_fill_identity() + xlab("Longitude") + ylab("Latitude") + labs(color="Social Group") + ggtitle("2019")


#In 2017, YLI has a giant homerange that extends way over the water. Some other home ranges encounter
#water slighltly, but sometimes they do forage on vegetation in the water and in some years the water was
#quite receded. But YLI is extreme - look at shapefile in Google Earth Engine. So going to remove her from 
#analyses. YLI gets removed at the bad17 step - can take her out of that step to see what her range looks like
# spYLI = sp17s[sp17s@data$name=="YLI 85% est",]
# sp::plot(spYLI)

# #Remove YLI from spatial polygons dataframe
# sp17s = sp17s[sp17s@data$name != "YLI 85% est",]
# sp::plot(sp17s) #Looks better

###Write shapefiles for Google Earth Engine 
library(rgdal)
# setwd("/Users/Joe/Documents/Research/RBFW Data/2016 Brisbane/NDVI 2016")
# writeOGR(sp16s, dsn = '.', layer = 'akde_poly16', driver = "ESRI Shapefile")
# setwd("/Users/Joe/Documents/Research/RBFW Data/2017 Brisbane/NDVI 2017")
# writeOGR(sp17s, dsn = '.', layer = 'adke_poly17', driver = "ESRI Shapefile")
# setwd("/Users/Joe/Documents/Research/RBFW Data/2018 Brisbane/NDVI 2018")
# writeOGR(sp18s, dsn = '.', layer = 'adke_poly18', driver = "ESRI Shapefile")
# setwd("/Users/Joe/Documents/Research/RBFW Data/2019 Brisbane/NDVI 2019")
# writeOGR(sp19s, dsn = '.', layer = 'adke_poly19', driver = "ESRI Shapefile")

#DID NOT END UP USING THESE HOME RANGES FOR NDVI MEASUREMENTS - these were more likely to have a little bit of 
#overlap with vegetation in the water which would've messed up NDVI. Used 95% KDE areas instead for NDVI. 

#Shapefiles then get imported as assets into GEE as zipped files. They're then used to specify 
#area of interest where NDVI is calculated. Download last 20 years data for each year and will 
#bring back into R to subset for dates of interest for each year to calculate NDVI leading up to 
#each year. GEE converts projection to EPSG:3857 (WGS84) when shapefiles are imported. 
#Have to project at  as.telemetry step to use in Google Earth Engine. GEE projects shapefile into ESPG:4236
#when importing the shape file. If it's not projected, an error comes up saying it couldn't change it from
#the flat projection



####Compare home range size across time periods####


#Using bird.gps15 lists - break up that list into time periods then subset for birds seen at least 5/10 times 
#in each time period
#First get julian dates in bird.gps15
bird16.gps15$jdate = yday(bird16.gps15$Date)
bird17.gps15$jdate = yday(bird17.gps15$Date)
bird18.gps15$jdate = yday(bird18.gps15$Date)
bird19.gps15$jdate = yday(bird19.gps15$Date)

##Get bird.gps into periods that were selected above 
#2016
bird16.gps15.1 = bird16.gps15 %>% filter(jdate<192)
bird16.gps15.2 = bird16.gps15 %>% filter(jdate>=192,jdate<205)
bird16.gps15.3 = bird16.gps15 %>% filter(jdate>=205)
#2017
bird17.gps15.1 = bird17.gps15 %>% filter(jdate<195)
bird17.gps15.2 = bird17.gps15 %>% filter(jdate>=195,jdate<211)
bird17.gps15.3 = bird17.gps15 %>% filter(jdate>=211)
#2018
bird18.gps15.1 = bird18.gps15 %>% filter(jdate<201)
bird18.gps15.2 = bird18.gps15 %>% filter(jdate>=201,jdate<210)
bird18.gps15.3 = bird18.gps15 %>% filter(jdate>=210)
#2019
bird19.gps15.1 = bird19.gps15 %>% filter(jdate<197)
bird19.gps15.2 = bird19.gps15 %>% filter(jdate>=197,jdate<210)
bird19.gps15.3 = bird19.gps15 %>% filter(jdate>=210)

##Subset bird.gps birds for birds that have at least 5/10 observations in each period
#2016
bird16.gps15.1 = bird16.gps15.1[bird16.gps15.1$individual.local.identifier %in% bird16.gps15.list.1L$individual.local.identifier,]
bird16.gps15.2 = bird16.gps15.2[bird16.gps15.2$individual.local.identifier %in% bird16.gps15.list.2L$individual.local.identifier,]
bird16.gps15.3 = bird16.gps15.3[bird16.gps15.3$individual.local.identifier %in% bird16.gps15.list.3L$individual.local.identifier,]
#2017
bird17.gps15.1 = bird17.gps15.1[bird17.gps15.1$individual.local.identifier %in% bird17.gps15.list.1L$individual.local.identifier,]
bird17.gps15.2 = bird17.gps15.2[bird17.gps15.2$individual.local.identifier %in% bird17.gps15.list.2L$individual.local.identifier,]
bird17.gps15.3 = bird17.gps15.3[bird17.gps15.3$individual.local.identifier %in% bird17.gps15.list.3L$individual.local.identifier,]
#2018
bird18.gps15.1 = bird18.gps15.1[bird18.gps15.1$individual.local.identifier %in% bird18.gps15.list.1L$individual.local.identifier,]
bird18.gps15.2 = bird18.gps15.2[bird18.gps15.2$individual.local.identifier %in% bird18.gps15.list.2L$individual.local.identifier,]
bird18.gps15.3 = bird18.gps15.3[bird18.gps15.3$individual.local.identifier %in% bird18.gps15.list.3L$individual.local.identifier,]
#2019
bird19.gps15.1 = bird19.gps15.1[bird19.gps15.1$individual.local.identifier %in% bird19.gps15.list.1L$individual.local.identifier,]
bird19.gps15.2 = bird19.gps15.2[bird19.gps15.2$individual.local.identifier %in% bird19.gps15.list.2L$individual.local.identifier,]
bird19.gps15.3 = bird19.gps15.3[bird19.gps15.3$individual.local.identifier %in% bird19.gps15.list.3L$individual.local.identifier,]

#Get a new akde.period column in each birdlist
birdlist16akde15$akde.period = NA
birdlist17akde15$akde.period = NA
birdlist18akde15$akde.period = NA
birdlist19akde15$akde.period = NA

#Subset the birdlists for list of birds in each time period
#2016
birdlist16akde15.1 = birdlist16akde15[birdlist16akde15$Bird %in% bird16.gps15.list.1L$individual.local.identifier,]
birdlist16akde15.2 = birdlist16akde15[birdlist16akde15$Bird %in% bird16.gps15.list.2L$individual.local.identifier,]
birdlist16akde15.3 = birdlist16akde15[birdlist16akde15$Bird %in% bird16.gps15.list.3L$individual.local.identifier,]
#2017
birdlist17akde15.1 = birdlist17akde15[birdlist17akde15$Bird %in% bird17.gps15.list.1L$individual.local.identifier,]
birdlist17akde15.2 = birdlist17akde15[birdlist17akde15$Bird %in% bird17.gps15.list.2L$individual.local.identifier,]
birdlist17akde15.3 = birdlist17akde15[birdlist17akde15$Bird %in% bird17.gps15.list.3L$individual.local.identifier,]
#2018
birdlist18akde15.1 = birdlist18akde15[birdlist18akde15$Bird %in% bird18.gps15.list.1L$individual.local.identifier,]
birdlist18akde15.2 = birdlist18akde15[birdlist18akde15$Bird %in% bird18.gps15.list.2L$individual.local.identifier,]
birdlist18akde15.3 = birdlist18akde15[birdlist18akde15$Bird %in% bird18.gps15.list.3L$individual.local.identifier,]
#2019
birdlist19akde15.1 = birdlist19akde15[birdlist19akde15$Bird %in% bird19.gps15.list.1L$individual.local.identifier,]
birdlist19akde15.2 = birdlist19akde15[birdlist19akde15$Bird %in% bird19.gps15.list.2L$individual.local.identifier,]
birdlist19akde15.3 = birdlist19akde15[birdlist19akde15$Bird %in% bird19.gps15.list.3L$individual.local.identifier,]


#Function to calculate akde for each individual and then put akde in a list
akde.rbfw.period = function(birdlist,bird.gps,UD,akde.list) {
  par(mfrow=c(5,5))
  for (i in 1:nrow(birdlist)) {
    #Get bird
    bird = birdlist$Bird[i]
    #Subset bird.gps to that bird
    bird.gps2 = bird.gps[which(bird.gps$individual.local.identifier==bird),]
    #Get into ctmm object
    HR <- as.telemetry(bird.gps2, timeformat="%Y-%m-%d $H:$M:$S", timezone="Australia/Brisbane",projection = sp::CRS("+init=epsg:32756"))
    #Get variogram
    var.HR = variogram(HR)
    #Plot variogram
    plot(var.HR,main=bird)
    #Get the fit parameters
    GUESS = variogram.fit(var.HR,interactive=F,name="GUESS")
    #Fit model
    mod.HR <- ctmm.select(HR, CTMM=GUESS, verbose=F)
    #Get kernel
    akde.HR<- akde(HR, CTMM=mod.HR)
    #Put akde estimate into birdlist
    birdlist$akde.period[i] = summary(akde.HR, level.UD=UD)$CI[2]
    #Make akde.HR into a spatial polygons dataframe
    akde.HR.sp = SpatialPolygonsDataFrame.UD(akde.HR,level.UD =0.85)[2,] #use just the estimate polygon
    #If it's the first row in the birdlist, then don't rbind it to the previous versions
    if(i==1) {akde.HR.sp2 = akde.HR.sp} else {akde.HR.sp2 = rbind(akde.HR.sp2,akde.HR.sp)}
  }
  par(mfrow=c(1,1))
  #Add birdlist to the list as well - only way to return multiple items
  akde.list[[1]] = birdlist
  akde.list[[2]] = akde.HR.sp2
  return(akde.list)
}


#Sometimes need to run these again
dev.off()
par(mar=c(1,1,1,1))


##Calculate akde for each period
#2016
# akde16.list.1 = akde.rbfw.period(birdlist=birdlist16akde15.1,bird.gps=bird16.gps15.1,UD=0.85,akde.list = akde16.listb)
# saveRDS(here::here("Output files/Homerange output files",akde16.list.1,"akde16_list_1.rds"))
# akde16.list.2 = akde.rbfw.period(birdlist=birdlist16akde15.2,bird.gps=bird16.gps15.2,UD=0.85,akde.list = akde16.listb)
# saveRDS(here::here("Output files/Homerange output files",akde16.list.2,"akde16_list_2.rds"))
# akde16.list.3 = akde.rbfw.period(birdlist=birdlist16akde15.3,bird.gps=bird16.gps15.3,UD=0.85,akde.list = akde16.listb)
# saveRDS(here::here("Output files/Homerange output files",akde16.list.3,"akde16_list_3.rds"))
#2017
# akde17.list.1 = akde.rbfw.period(birdlist=birdlist17akde15.1,bird.gps=bird17.gps15.1,UD=0.85,akde.list = akde17.listb)
# saveRDS(here::here("Output files/Homerange output files",akde17.list.1,"akde17_list_1.rds"))
# akde17.list.2 = akde.rbfw.period(birdlist=birdlist17akde15.2,bird.gps=bird17.gps15.2,UD=0.85,akde.list = akde17.listb)
# saveRDS(here::here("Output files/Homerange output files",akde17.list.2,"akde17_list_2.rds"))
# akde17.list.3 = akde.rbfw.period(birdlist=birdlist17akde15.3,bird.gps=bird17.gps15.3,UD=0.85,akde.list = akde17.listb)
# saveRDS(here::here("Output files/Homerange output files",akde17.list.3,"akde17_list_3.rds"))
#2018
# akde18.list.1 = akde.rbfw.period(birdlist=birdlist18akde15.1,bird.gps=bird18.gps15.1,UD=0.85,akde.list = akde18.listb)
# saveRDS(here::here("Output files/Homerange output files",akde18.list.1,"akde18_list_1.rds"))
# akde18.list.2 = akde.rbfw.period(birdlist=birdlist18akde15.2,bird.gps=bird18.gps15.2,UD=0.85,akde.list = akde18.listb)
# saveRDS(here::here("Output files/Homerange output files",akde18.list.2,"akde18_list_2.rds"))
# akde18.list.3 = akde.rbfw.period(birdlist=birdlist18akde15.3,bird.gps=bird18.gps15.3,UD=0.85,akde.list = akde18.listb)
# saveRDS(here::here("Output files/Homerange output files",akde18.list.3,"akde18_list_3.rds"))
#2019
# akde19.list.1 = akde.rbfw.period(birdlist=birdlist19akde15.1,bird.gps=bird19.gps15.1,UD=0.85,akde.list = akde19.listb)
# saveRDS(here::here("Output files/Homerange output files",akde19.list.1,"akde19_list_1.rds"))
# akde19.list.2 = akde.rbfw.period(birdlist=birdlist19akde15.2,bird.gps=bird19.gps15.2,UD=0.85,akde.list = akde19.listb)
# saveRDS(here::here("Output files/Homerange output files",akde19.list.2,"akde19_list_2.rds"))
# akde19.list.3 = akde.rbfw.period(birdlist=birdlist19akde15.3,bird.gps=bird19.gps15.3,UD=0.85,akde.list = akde19.listb)
# saveRDS(here::here("Output files/Homerange output files",akde19.list.3,"akde19_list_3.rds"))

##Read in adke.list.# lists
#2016
akde16.list.1 = readRDS(here::here("Output files/Homerange output files","akde16_list_1.rds"))
akde16.list.2 = readRDS(here::here("Output files/Homerange output files","akde16_list_2.rds"))
akde16.list.3 = readRDS(here::here("Output files/Homerange output files","akde16_list_3.rds"))
#2017
akde17.list.1 = readRDS(here::here("Output files/Homerange output files","akde17_list_1.rds"))
akde17.list.2 = readRDS(here::here("Output files/Homerange output files","akde17_list_2.rds"))
akde17.list.3 = readRDS(here::here("Output files/Homerange output files","akde17_list_3.rds"))
#2018
akde18.list.1 = readRDS(here::here("Output files/Homerange output files","akde18_list_1.rds"))
akde18.list.2 = readRDS(here::here("Output files/Homerange output files","akde18_list_2.rds"))
akde18.list.3 = readRDS(here::here("Output files/Homerange output files","akde18_list_3.rds"))
#2019
akde19.list.1 = readRDS(here::here("Output files/Homerange output files","akde19_list_1.rds"))
akde19.list.2 = readRDS(here::here("Output files/Homerange output files","akde19_list_2.rds"))
akde19.list.3 = readRDS(here::here("Output files/Homerange output files","akde19_list_3.rds"))


##Get birdlists from akde.list.# lists
#2016
birdlist16akde15.1 = akde16.list.1[[1]] #Birdlist is first in list
birdlist16akde15.2 = akde16.list.2[[1]] #Birdlist is first in list
birdlist16akde15.3 = akde16.list.3[[1]] #Birdlist is first in list
#2017
birdlist17akde15.1 = akde17.list.1[[1]] #Birdlist is first in list
birdlist17akde15.2 = akde17.list.2[[1]] #Birdlist is first in list
birdlist17akde15.3 = akde17.list.3[[1]] #Birdlist is first in list
#2018
birdlist18akde15.1 = akde18.list.1[[1]] #Birdlist is first in list
birdlist18akde15.2 = akde18.list.2[[1]] #Birdlist is first in list
birdlist18akde15.3 = akde18.list.3[[1]] #Birdlist is first in list
#2019
birdlist19akde15.1 = akde19.list.1[[1]] #Birdlist is first in list
birdlist19akde15.2 = akde19.list.2[[1]] #Birdlist is first in list
birdlist19akde15.3 = akde19.list.3[[1]] #Birdlist is first in list

##Add a period column in each period birdlist
#2016
birdlist16akde15.1$period = as.factor("1")
birdlist16akde15.2$period = as.factor("2")
birdlist16akde15.3$period = as.factor("3")
#2017
birdlist17akde15.1$period = as.factor("1")
birdlist17akde15.2$period = as.factor("2")
birdlist17akde15.3$period = as.factor("3")
#2018
birdlist18akde15.1$period = as.factor("1")
birdlist18akde15.2$period = as.factor("2")
birdlist18akde15.3$period = as.factor("3")
#2019
birdlist19akde15.1$period = as.factor("1")
birdlist19akde15.2$period = as.factor("2")
birdlist19akde15.3$period = as.factor("3")

##Combine periods into one dataframe for each year
birdlist16akde15.all = rbind(birdlist16akde15.1,birdlist16akde15.2,birdlist16akde15.3)
birdlist17akde15.all = rbind(birdlist17akde15.1,birdlist17akde15.2,birdlist17akde15.3)
birdlist18akde15.all = rbind(birdlist18akde15.1,birdlist18akde15.2,birdlist18akde15.3)
birdlist19akde15.all = rbind(birdlist19akde15.1,birdlist19akde15.2,birdlist19akde15.3)

##If akde.period is in meters, convert to hectares
for (i in 1:nrow(birdlist16akde15.all)) {
  if(birdlist16akde15.all$akde.period[i]>1000) 
  {birdlist16akde15.all$akde.period[i]=(birdlist16akde15.all$akde.period[i]/10000)}
  }
for (i in 1:nrow(birdlist17akde15.all)) {
  if(birdlist17akde15.all$akde.period[i]>1000) 
  {birdlist17akde15.all$akde.period[i]=(birdlist17akde15.all$akde.period[i]/10000)}
}
for (i in 1:nrow(birdlist18akde15.all)) {
  if(birdlist18akde15.all$akde.period[i]>1000) 
  {birdlist18akde15.all$akde.period[i]=(birdlist18akde15.all$akde.period[i]/10000)}
}
for (i in 1:nrow(birdlist19akde15.all)) {
  if(birdlist19akde15.all$akde.period[i]>1000) 
  {birdlist19akde15.all$akde.period[i]=(birdlist19akde15.all$akde.period[i]/10000)}
}

#Get social group as a factor
birdlist16akde15.all$Social.Group = as.factor(birdlist16akde15.all$Social.Group)
birdlist17akde15.all$Social.Group = as.factor(birdlist17akde15.all$Social.Group)
birdlist18akde15.all$Social.Group = as.factor(birdlist18akde15.all$Social.Group)
birdlist19akde15.all$Social.Group = as.factor(birdlist19akde15.all$Social.Group)



####Model differences in mean akde across periods####
model.akde16 = lmer(akde.period~period + (1|Bird) + (1|class) + (1|Social.Group),data = birdlist16akde15.all)
summary(model.akde16)
model.akde17 = lmer(akde.period~period + (1|Bird) + (1|class) + (1|Social.Group),data = birdlist17akde15.all)
summary(model.akde17)
model.akde18 = lmer(akde.period~period + (1|Bird) + (1|class) + (1|Social.Group),data = birdlist18akde15.all)
summary(model.akde18)
model.akde19 = lmer(akde.period~period + (1|Bird) + (1|class) + (1|Social.Group),data = birdlist19akde15.all)
summary(model.akde19)


#emmeans comparisons
library(emmeans)
#2016
model.akde16.e = emmeans(model.akde16,specs="period")
model.akde16.e
model.akde16.e.df = data.frame(model.akde16.e)
pairs(model.akde16.e)
model.akde16.e.1.2 = data.frame(pairs(model.akde16.e))$estimate[1]
model.akde16.e.1.3 = data.frame(pairs(model.akde16.e))$estimate[2]
model.akde16.e.2.3 = data.frame(pairs(model.akde16.e))$estimate[3]
#2017
model.akde17.e = emmeans(model.akde17,specs="period")
model.akde17.e
model.akde17.e.df = data.frame(model.akde17.e)
pairs(model.akde17.e)
model.akde17.e.1.2 = data.frame(pairs(model.akde17.e))$estimate[1]
model.akde17.e.1.3 = data.frame(pairs(model.akde17.e))$estimate[2]
model.akde17.e.2.3 = data.frame(pairs(model.akde17.e))$estimate[3]
#2018
model.akde18.e = emmeans(model.akde18,specs="period")
model.akde18.e
model.akde18.e.df = data.frame(model.akde18.e)
pairs(model.akde18.e)
model.akde18.e.1.2 = data.frame(pairs(model.akde18.e))$estimate[1]
model.akde18.e.1.3 = data.frame(pairs(model.akde18.e))$estimate[2]
model.akde18.e.2.3 = data.frame(pairs(model.akde18.e))$estimate[3]
#2019
model.akde19.e = emmeans(model.akde19,specs="period")
model.akde19.e
model.akde19.e.df = data.frame(model.akde19.e)
pairs(model.akde19.e)
model.akde19.e.1.2 = data.frame(pairs(model.akde19.e))$estimate[1]
model.akde19.e.1.3 = data.frame(pairs(model.akde19.e))$estimate[2]
model.akde19.e.2.3 = data.frame(pairs(model.akde19.e))$estimate[3]

#Those are the comparisons based off of normal emmeans, but I should compare to randomizations to see if 
#observed differences between time periods are more than expected by chance. Randomizations are a better way
#to go becuase 5 observations are not necessarily a good approximation of an individual's home range. 
#In randomizations, use same individuals and same number of observations of each individual, but randomly 
#sample observations from across the entire season to see if observed difference is more than could 
#have happend by random. 

####Write .csv files for cluster randomizations####
##2016
# write.csv(bird16.gps15,here::here("Output files/Homerange output files","bird16_gps15.csv"))
# write.csv(bird16.gps15.list.1L,here::here("Output files/Homerange output files","bird16_gps15_list_1L.csv"))
# write.csv(bird16.gps15.list.2L,here::here("Output files/Homerange output files","bird16_gps15_list_2L.csv"))
# write.csv(bird16.gps15.list.3L,here::here("Output files/Homerange output files","bird16_gps15_list_3L.csv"))
# write.csv(birdlist16akde15.1,here::here("Output files/Homerange output files","birdlist16akde15_1.csv"))
# write.csv(birdlist16akde15.2,here::here("Output files/Homerange output files","birdlist16akde15_2.csv"))
# write.csv(birdlist16akde15.3,here::here("Output files/Homerange output files","birdlist16akde15_3.csv"))
##2017
# write.csv(bird17.gps15,"bird17_gps15.csv"))
# write.csv(bird17.gps15.list.1L,here::here("Output files/Homerange output files","bird17_gps15_list_1L.csv"))
# write.csv(bird17.gps15.list.2L,here::here("Output files/Homerange output files","bird17_gps15_list_2L.csv"))
# write.csv(bird17.gps15.list.3L,here::here("Output files/Homerange output files","bird17_gps15_list_3L.csv"))
# write.csv(birdlist17akde15.1,here::here("Output files/Homerange output files","birdlist17akde15_1.csv"))
# write.csv(birdlist17akde15.2,here::here("Output files/Homerange output files","birdlist17akde15_2.csv"))
# write.csv(birdlist17akde15.3,here::here("Output files/Homerange output files","birdlist17akde15_3.csv"))
##2018
# write.csv(bird18.gps15,"bird18_gps15.csv"))
# write.csv(bird18.gps15.list.1L,here::here("Output files/Homerange output files","bird18_gps15_list_1L.csv"))
# write.csv(bird18.gps15.list.2L,here::here("Output files/Homerange output files","bird18_gps15_list_2L.csv"))
# write.csv(bird18.gps15.list.3L,here::here("Output files/Homerange output files","bird18_gps15_list_3L.csv"))
# write.csv(birdlist18akde15.1,here::here("Output files/Homerange output files","birdlist18akde15_1.csv"))
# write.csv(birdlist18akde15.2,here::here("Output files/Homerange output files","birdlist18akde15_2.csv"))
# write.csv(birdlist18akde15.3,here::here("Output files/Homerange output files","birdlist18akde15_3.csv"))
##2019
# write.csv(bird19.gps15,here::here("Output files/Homerange output files","bird19_gps15.csv"))
# write.csv(bird19.gps15.list.1L,here::here("Output files/Homerange output files","bird19_gps15_list_1L.csv"))
# write.csv(bird19.gps15.list.2L,here::here("Output files/Homerange output files","bird19_gps15_list_2L.csv"))
# write.csv(bird19.gps15.list.3L,here::here("Output files/Homerange output files","bird19_gps15_list_3L.csv"))
# write.csv(birdlist19akde15.1,here::here("Output files/Homerange output files","birdlist19akde15_1.csv"))
# write.csv(birdlist19akde15.2,here::here("Output files/Homerange output files","birdlist19akde15_2.csv"))
# write.csv(birdlist19akde15.3,here::here("Output files/Homerange output files","birdlist19akde15_3.csv"))



####Determine significance across time periods using randomizations####

#Load randomization outputs
akde.rand16 = read.csv(here::here("Output files","akde_rand16.csv"))
akde.rand17 = read.csv(here::here("Output files","akde_rand17.csv"))
akde.rand18 = read.csv(here::here("Output files","akde_rand18.csv")) 
akde.rand19 = read.csv(here::here("Output files","akde_rand19.csv"))

##Get p-values by comparing observed emmeans pair comparison estimates to randomized emmeans 
#pair comparison estimates

#2016 
1 - sum(model.akde16.e.1.2>akde.rand16$X1.2)/1000 #1 vs 2
1 - sum(model.akde16.e.1.3>akde.rand16$X1.3)/1000 #1 vs 3
1 - sum(model.akde16.e.2.3>akde.rand16$X2.3)/1000 #2 vs 3
#2017 
1 - sum(model.akde17.e.1.2>akde.rand17$X1.2)/1000 #1 vs 2
1 - sum(model.akde17.e.1.3>akde.rand17$X1.3)/1000 #1 vs 3
1 - sum(model.akde17.e.2.3>akde.rand17$X2.3)/1000 #2 vs 3
#2018 
1 - sum(model.akde18.e.1.2>akde.rand18$X1.2)/1000 #1 vs 2
1 - sum(model.akde18.e.1.3>akde.rand18$X1.3)/1000 #1 vs 3
1 - sum(model.akde18.e.2.3>akde.rand18$X2.3)/1000 #2 vs 3
#2019 
1 - sum(model.akde19.e.1.2>akde.rand19$X1.2)/1000 #1 vs 2
1 - sum(model.akde19.e.1.3>akde.rand19$X1.3)/1000 #1 vs 3
1 - sum(model.akde19.e.2.3>akde.rand19$X2.3)/1000 #2 vs 3


####Plot homerange differences across time periods####

#Get middle dates of each time period for plotting - round up
#2016
first16 = round((191-bird16.gps15.1$jdate[1])/2 + bird16.gps15.1$jdate[1]) 
middle16 = round((204-192)/2 + 192) 
third16 = round((bird16.gps15.3$jdate[nrow(bird16.gps15.3)]-205)/2 + 205)
#2017
first17 = round((194-bird17.gps15.1$jdate[1])/2 + bird17.gps15.1$jdate[1]) 
middle17 = round((210-195)/2 + 195) 
third17 = round((bird17.gps15.3$jdate[nrow(bird17.gps15.3)]-211)/2 + 211)
#2018
first18 = round((200-bird18.gps15.1$jdate[1])/2 + bird18.gps15.1$jdate[1]) 
middle18 = round((209-201)/2 + 201) 
third18 = round((bird18.gps15.3$jdate[nrow(bird18.gps15.3)]-210)/2 + 210)
#2019
first19 = round((196-bird19.gps15.1$jdate[1])/2 + bird19.gps15.1$jdate[1]) 
middle19 = round((209-197)/2 + 197)
third19 = round((bird19.gps15.3$jdate[nrow(bird19.gps15.3)]-210)/2 + 210)

#Add dates to model.akde.e.df - use middle dates of each of the time periods
model.akde16.e.df$jdate = c(first16,middle16,third16)
model.akde17.e.df$jdate = c(first17,middle17,third17)
model.akde18.e.df$jdate = c(first18,middle18,third18)
model.akde19.e.df$jdate = c(first19,middle19,third19)

#Add Year to model.akde.e.dfs
model.akde16.e.df$Year="2016"
model.akde17.e.df$Year="2017"
model.akde18.e.df$Year="2018"
model.akde19.e.df$Year="2019"

#Combine model.akde.e.dfs to plot
model.akde.e.dfall = rbind(model.akde16.e.df,model.akde17.e.df,model.akde18.e.df,model.akde19.e.df)

#Write 2017 data for 2017 plot
#write.csv(model.akde17.e.df,here::here("Output files","2017 akde HR summarized.csv"),row.names=F,quote=F)

#Plot
ggplot() +
  geom_line(data=model.akde.e.dfall,aes(x=jdate,y=emmean,group=Year,color=Year),size=1) +
  geom_errorbar(data=model.akde.e.dfall,aes(x=jdate,y=emmean,ymax=emmean+SE,ymin=emmean-SE,width=2,group=Year,color=Year),size=1) +
  geom_point(data=model.akde.e.dfall,aes(x=jdate,y=emmean,group=Year,color=Year),size=3) +
  theme_cowplot() + scale_colour_manual(name="Season",values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  xlab("Julian date") + ylab("Mean aKDE (hectares)") + ylim(0,12.5) 
  
#Add significance (a,b,c) by color (year) in pdf. Each period will get a letter - 1st will be a|a|a|a then
#next is b|a|a|b   next is b|b|a|c  where letters are all colored to their year

#Will want to include oberserved coefficients, random coefficients (maybe mean and 95% range) and significance 
#values in a table for supplemental. 


##Plot by date with points different shapes for printing in black and white
ggplot() +
  geom_line(data=model.akde.e.dfall,aes(x=jdate,y=emmean,group=Year,color=Year),size=1.3) +
  geom_errorbar(data=model.akde.e.dfall,aes(x=jdate,y=emmean,ymax=emmean+SE,ymin=emmean-SE,width=2,group=Year,color=Year),size=1) +
  geom_point(data=model.akde.e.dfall,aes(x=jdate,y=emmean,group=Year,color=Year,shape=Year),size=4) +
  theme_cowplot() + scale_colour_manual(name="Season",values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  scale_shape_manual(name="Season",values=c(19,15,17,18)) +
  xlab("Julian date") + ylab("Mean aKDE (hectares)") + ylim(0,12.5) 

ggplot() +
  geom_line(data=model.akde.e.dfall,aes(x=jdate,y=emmean,group=Year,color=Year),size=1.3) +
  geom_errorbar(data=model.akde.e.dfall,aes(x=jdate,y=emmean,ymax=emmean+SE,ymin=emmean-SE,width=2,group=Year,color=Year),size=1) +
  geom_point(data=model.akde.e.dfall,aes(x=jdate,y=emmean,group=Year,fill=Year,shape=Year),size=4) +
  theme_cowplot() + scale_colour_manual(name="Season",values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  scale_shape_manual(name="Season",values=c(21,22,24,23)) + 
  scale_fill_manual(name="Season",values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  xlab("Julian date") + ylab("Mean aKDE (hectares)") + ylim(0,12.5) 






##Plot histograms of observed vs random for home range size differences
A = ggplot(data=akde.rand16, aes(x=X1.2)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde16.e.1.2,color="red",size=1.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2016: Period 1 vs Period 2; p=0.030") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
B = ggplot(data=akde.rand16, aes(x=X1.3)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde16.e.1.3,color="red",size=1.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2016: Period 1 vs Period 3; p=0.002") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
C = ggplot(data=akde.rand16, aes(x=X2.3)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde16.e.2.3,color="red",size=1.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2016: Period 2 vs Period 3; p=0.073") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
D = ggplot(data=akde.rand17, aes(x=X1.2)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde17.e.1.2,color="red",size=1.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2017: Period 1 vs Period 2; p=0.278") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
E = ggplot(data=akde.rand17, aes(x=X1.3)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde17.e.1.3,color="red",size=1.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2017: Period 1 vs Period 3; p=0.006") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
G = ggplot(data=akde.rand17, aes(x=X2.3)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde17.e.2.3,color="red",size=1.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2017: Period 2 vs Period 3; p=0.055") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
H = ggplot(data=akde.rand18, aes(x=X1.2)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde18.e.1.2,color="red",size=0.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2018: Period 1 vs Period 2; p=0.076") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
I = ggplot(data=akde.rand18, aes(x=X1.3)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde18.e.1.3,color="red",size=0.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2018: Period 1 vs Period 3; p=0.068") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
J = ggplot(data=akde.rand18, aes(x=X2.3)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde18.e.2.3,color="red",size=0.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2018: Period 2 vs Period 3; p=0.543") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
K = ggplot(data=akde.rand19, aes(x=X1.2)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde19.e.1.2,color="red",size=1.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2019: Period 1 vs Period 2; p=0.807") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
L = ggplot(data=akde.rand19, aes(x=X1.3)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde19.e.1.3,color="red",size=1.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2019: Period 1 vs Period 3; p=0.037") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
M = ggplot(data=akde.rand19, aes(x=X2.3)) + geom_histogram(binwidth = 0.5,fill="gray",color="black") + 
  geom_vline(xintercept = model.akde19.e.2.3,color="red",size=1.5) + xlab("aKDE difference (hectares)") +
  ggtitle("2019: Period 2 vs Period 3; p=0.011") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))


library(ggpubr)
ggarrange(A,B,C,D,E,G,H,I,J,K,L,M,ncol=3,nrow=4)





