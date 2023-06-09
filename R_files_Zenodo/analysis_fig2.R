library(magrittr)
library(scales)
library(lubridate)
library(viridis)
library(geosphere)
library(maps)
library(maptools)


#### IMPORTING AND PROCESSING TEMPERATURE DATA------------------------------------------------------

raw_Sth <- read.delim("input_T_erken05-16.txt")
raw_Hgl <- read.delim("input_T_tomtabacken05-16.txt")
raw_Got <- read.delim("input_T_gotlandhemse05-16.txt")
raw_Ssk <- read.delim("input_T_lund05-16.txt")

raw_Sth$day <- yday(raw_Sth$datum)
raw_Hgl$day <- yday(raw_Hgl$datum)
raw_Got$day <- yday(raw_Got$datum)
raw_Ssk$day <- yday(raw_Ssk$datum)

raw_Sth$year <- year(raw_Sth$datum)
raw_Hgl$year <- year(raw_Hgl$datum)
raw_Got$year <- year(raw_Got$datum)
raw_Ssk$year <- year(raw_Ssk$datum)

# What years are shared?
#plot(raw_Sth$day ~ raw_Sth$year, main="Stockholm")
#plot(raw_Hgl$day ~ raw_Hgl$year, main="Höglandet")
#plot(raw_Got$day ~ raw_Got$year, main="Gotland")
#plot(raw_Ssk$day ~ raw_Ssk$year, main="Skåne")

# Filter out 2014, because Stockholm had no data at all then;
# also 2016 because some populations lack the end of that year
raw_Sth <- subset(raw_Sth, year!=2014 & year!=2016)
raw_Hgl <- subset(raw_Hgl, year!=2014 & year!=2016)
raw_Got <- subset(raw_Got, year!=2014 & year!=2016)
raw_Ssk <- subset(raw_Ssk, year!=2014 & year!=2016)

# Substitute 2005 for 2014 so that the years form an unbroken row (for easier plotting)
raw_Sth[raw_Sth$year==2005,]$year <- 2014
raw_Hgl[raw_Hgl$year==2005,]$year <- 2014
raw_Got[raw_Got$year==2005,]$year <- 2014
raw_Ssk[raw_Ssk$year==2005,]$year <- 2014


standardcurve <- function(input) {
  
  # Parse and average temperature series
  meancurve <- tapply(input$temp,input$day,mean)
  smoothcurve <- filter(meancurve, rep(1/21,21), sides=2, circular=TRUE) # Smoothed average curve

  detrend <- NULL
  # Temperature variance relative to the between-year mean is higher in winter than in summer.
  # As the model focuses on simulating the growing season, only April-September variance is used as input.
  summer <- subset(input, day>91 & day<274)
  for (i in 1:nrow(summer)) { detrend[i] <- summer[i,]$temp - smoothcurve[summer[i,]$day] }
  stdev <- sd(detrend)
  
    list(smoothcurve=smoothcurve, stdev=stdev)
}


ptempcurves <- list(Sth=standardcurve(raw_Sth),
                    Hgl=standardcurve(raw_Hgl),
                    Got=standardcurve(raw_Got),
                    Ssk=standardcurve(raw_Ssk))




#### IMPORTING AND PROCESSING OBSERVATION DATA -----------------------------------------------------------

paemap <- read.delim("input_artportalen_2018.txt")
paemap$duration <- as.numeric(ymd(paemap$enddate))-as.numeric(ymd(paemap$startdate))
paemap$year <- year(paemap$startdate)
paemap$day <- yday(paemap$startdate) + floor(paemap$duration*0.5) # Take middle day of reported duration (round down) for obs. lasting 1-10 days
paemap <- subset(paemap, duration < 10 & day<300)

# Function for averaging observation densitites across the chosen years for a given geographic subset
observations <- function(dataset, from, to, bandwidth) {
  tabulation <- data.frame(density(subset(dataset, year == 2006)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2007)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2008)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2009)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2010)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2011)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2012)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2013)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2014)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2015)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2016)$day, bw=bandwidth, from=from, to=to)$y,
                           density(subset(dataset, year == 2017)$day, bw=bandwidth, from=from, to=to)$y)
  
  tabulation$average <- NA
  for (i in 1:nrow(tabulation)) { tabulation$average[i] <- mean(as.numeric(tabulation[i,1:12])) }
  return(tabulation$average) }

# Function for extracting geographic subset (100 km radius from sampled population)
popradius <- function(loc){
  totset <- paemap[c(2,3,10,11)]
  totset$distance <- 0
  for(r in 1:nrow(totset)) {
    totset[r,]$distance <- distVincentyEllipsoid(c(totset[r,]$lon,totset[r,]$lat),loc)
  }
  return (subset(totset, distance<=100000)) }

locations <- list(Sth=c(18.53,59.63), Hgl=c(14.14,57.52), Got=c(18.53,57.40), Ska=c(13.4,55.60))
popsets <- lapply(locations, popradius)

# Assumed first-gen observations, for model input
obs.firstgen <- list(Sth=subset(popsets$Sth, day<=195)$day,
                     Hgl=subset(popsets$Hgl, day<=195)$day,
                     Got=subset(popsets$Got, day<=180)$day)

# Assumed second-gen observations, for model input
obs.secndgen <- list(Got=subset(popsets$Got, day>180)$day,
                     Ska=subset(popsets$Ska, day>180)$day)

# Density curves of observation data, for plotting
obs.curve <- list(Sth=observations(popsets$Sth, from=80, to=300, bandwidth=14),
                  Hgl=observations(popsets$Hgl, from=80, to=300, bandwidth=14),
                  Got=observations(popsets$Got, from=80, to=300, bandwidth=14),
                  Ska=observations(popsets$Ska, from=80, to=300, bandwidth=14))



#### SETTING UP OTHER MODEL PARAMETERS -----------------------------------------------------------------

# Thermal reaction norm for larval development (when daylength > LDT)
Tmodel <- lm(1/c(42.5,30,20,25.45,30.25,21.42)~c(13,17,21,18,16,20))
Tmin <- as.numeric(-coef(Tmodel)[1]/coef(Tmodel)[2])
Ticpt <- as.numeric(coef(Tmodel)[1])
Tslope <- as.numeric(coef(Tmodel)[2])
rm(Tmodel)

# Photoperiodic reaction norm for larval development (when daylength < LDT)
Dicpt <- -106
Dslope <- 10
Dcap <- 30 # Minimum development time under daylength regulation
Dmin <- (Dicpt-Dcap)/-Dslope # Fall daylength at which development stops speeding up

# Fixed CDs (from Lindestad et al 2019)
CDs <- as.list(read.delim("input_CDestimates.txt"))
latitudes <- list(Vag=57.52, Kal=mean(c(56.93,57.31)), Ola=56.62, Got=57.40, Nsk=56.29, Sth=59.63, Ssk=55.60)

# Reaction norm for temperature-sensitive critical daylength
CDmodel <- lm(c(CDs$Got, CDs$Ola, CDs$Vag, CDs$Kal)~c(16,20,16,20,16,20,16,20)+c("B","B","O","O","V","V","K","K"))
CDicpt_is <- as.numeric(coef(CDmodel)[1])
CDicpt_hi <- as.numeric(coef(CDmodel)[1]+coef(CDmodel)[5])
CDicpt_ka <- as.numeric(coef(CDmodel)[1]+coef(CDmodel)[3])
CDslope <- as.numeric(coef(CDmodel)[2])
rm(CDmodel)

# Extrapolated CD slopes for Skåne & Stockholm
CDicpt_is-((CDicpt_is+CDslope*18)-CDs$Ssk[1]) -> CDicpt_ss
CDicpt_is-((CDicpt_is+CDslope*18)-CDs$Sth[1]) -> CDicpt_st



#### DEFINING MODEL FUNCTIONS -----------------------------------------------------------------------

# Daylength calculation according to Forsythe et al (1995)
daylength <- function(day, p=1.5, L) {
  24-((24/pi)*acos((sin(p*pi/180)+sin(L*pi/180)*
                      (asin(0.39795*cos((0.2163108 + 2*atan(0.9671396*tan(
                        0.00860*(day-186))))))))/(cos(L*pi/180)*
                                                    cos((asin(0.39795*cos((0.2163108 + 2*atan(0.9671396*tan(0.00860*(day-186)))))))))))}

# Lindestad et al 2019 (Ecology) model function (version 3; see original paper for details)
makelarva <- function(temps, startday, latitude, smoothtemp=6, CD, LDT=CD+0.5, stochtemp=TRUE, CDicpt=NULL, years=10) {
  
  DLprofile <- daylength(1:365, L=latitude)
  output <- data.frame(year=NA, startday=NA, endday=NA, endtemp=NA, enddaylength=NA)
  basetemp <- temps$smoothcurve
  tempvariation <- temps$stdev
  
  for (y in 1:years) {
    
    tempseries <- basetemp+rnorm(366,0,stochtemp*tempvariation) # Simulate a unique thermal year (unless stochtemp=FALSE)
    
    # For each larva within that year:
    for (i in startday) {
      
      dev <- Ticpt + Tslope*tempseries # Development rate according to daily temperature
      dev <- dev*1*(tempseries>Tmin)    # No development when T < Tmin
      dev[1:i] <- 0  # Development begins at larval hatch date
      devsum <- cumsum(dev)
      
      # At molt to second instar, photoperiod can affect development
      instar234 <- length(devsum[devsum<0.25])+1
      if (instar234>364) print("TOO LATE!", quote=FALSE)
      
      # For each day after second instar begins, check LDT and let photoperiod regulate development where appropriate
      for (d in instar234:365) {
        if (DLprofile[d]<LDT) {dev[d] <- 1/(Dicpt+Dslope*DLprofile[d]) }
        if (DLprofile[d]<Dmin) {dev[d] <- 1/Dcap }
      }
      
      # Maturity occurs when devsum = 1
      devsum <- cumsum(dev)
      
      # Results for this individual:
      enddate <- length(devsum[devsum<1])+1
      howhotthatweek <- mean(tempseries[(enddate-smoothtemp):enddate])
      howlongthatday <- DLprofile[enddate]
      output <- rbind(output,c(y, i, enddate, howhotthatweek, howlongthatday))
    }
  }  
  
  output$diapause <- output$enddaylength < CD # Diapause decision for each individual (i.e. model version 2)
  if(!is.null(CDicpt)) { output$diapause_T <- output$enddaylength < output$endtemp*CDslope + CDicpt } # Temperature-sensitive diapause decision (model version 3)
  
  return(output[-1,])
}  




#### RUNNING LIFE CYCLE MODEL ------------------------------------------------------------------------

set.seed(413)
nitt=100

# Generate simulated start dates (stdev cut in half to omit extreme autumn dates the model was not built to handle)
#start.Sth  <- round(rnorm(100,mean(obs.firstgen$Sth),sd(obs.firstgen$Sth)/2) )
#start.Hgl  <- round(rnorm(100,mean(obs.firstgen$Hgl),sd(obs.firstgen$Hgl)/2) )
#start.Got1 <- round(rnorm(100,mean(obs.firstgen$Got),sd(obs.firstgen$Got)/2) )
#start.Got2 <- round(rnorm(100,mean(obs.secndgen$Got),sd(obs.secndgen$Got)/2) )
#start.Ska  <- round(rnorm(100,mean(obs.secndgen$Ska),sd(obs.secndgen$Ska)/2) )

# Run life cycle model (100 simulated years; N=100 for each year)
#run.gen1.Sth <- makelarva(ptempcurves$Sth, startday=start.Sth+10, latitude=latitudes$Sth, CD=CDs$Sth[1], CDicpt=CDicpt_st, years=nitt)
#run.gen1.Hgl <- makelarva(ptempcurves$Hgl, startday=start.Hgl+10, latitude=latitudes$Vag, CD=CDs$Vag[1], CDicpt=CDicpt_hi, years=nitt)
#run.gen1.Got <- makelarva(ptempcurves$Got, startday=start.Got1+10,latitude=latitudes$Got, CD=CDs$Got[1], CDicpt=CDicpt_is, years=nitt)
#run.gen2.Got <- makelarva(ptempcurves$Got, startday=start.Got2+10,latitude=latitudes$Got, CD=CDs$Got[1], CDicpt=CDicpt_is, years=nitt)
#run.gen2.Ska <- makelarva(ptempcurves$Ssk, startday=start.Ska+10, latitude=latitudes$Ssk, CD=CDs$Ssk[1], CDicpt=CDicpt_ss, years=nitt)


#### EXPORTING MODEL RUNS ----------------------------------------------------------------------------

#write.table(run.gen1.Sth, file="run.gen1.Sth", quote=FALSE, sep="\\t", row.names=FALSE)
#write.table(run.gen1.Hgl, file="run.gen1.Hgl", quote=FALSE, sep="\\t", row.names=FALSE)
#write.table(run.gen1.Got, file="run.gen1.Got", quote=FALSE, sep="\\t", row.names=FALSE)
#write.table(run.gen2.Got, file="run.gen2.Got", quote=FALSE, sep="\\t", row.names=FALSE)
#write.table(run.gen2.Ska, file="run.gen2.Ska", quote=FALSE, sep="\\t", row.names=FALSE)


#### IMPORTING MODEL RUNS ------------------------------------------------------------------------------

run.gen1.Sth <- read.delim("run.gen1.Sth")
run.gen1.Hgl <- read.delim("run.gen1.Hgl")
run.gen1.Got <- read.delim("run.gen1.Got")
run.gen2.Got <- read.delim("run.gen2.Got")
run.gen2.Ska <- read.delim("run.gen2.Ska")


#### EXTRACTING DIAPAUSE STARTING TIMES FROM MODEL OUTPUT

diapause.starts <- list(Sth= quantile(subset(run.gen1.Sth, diapause_T==TRUE)$endday, c(0.25,0.5,0.75)),
                        Hgl= quantile(subset(run.gen1.Hgl, diapause_T==TRUE)$endday, c(0.25,0.5,0.75)),
                        Got1=quantile(subset(run.gen1.Got, diapause_T==TRUE)$endday, c(0.25,0.5,0.75)),
                        Got2=quantile(subset(run.gen2.Got, diapause_T==TRUE)$endday, c(0.25,0.5,0.75)),
                        Ska= quantile(subset(run.gen2.Ska, diapause_T==TRUE)$endday, c(0.25,0.5,0.75)))

diapause.25th <- unlist(lapply(diapause.starts, function(x) x[1]))
diapause.medians <- unlist(lapply(diapause.starts, function(x) x[2]))
diapause.75th <- unlist(lapply(diapause.starts, function(x) x[3]))



#### DRAWING / EXPORTING FIGURE 2 ---------------------------------------------------------


pdf(file="Fig2.pdf", width=10, height=10, pointsize=20, onefile=TRUE)

# Set up plotting space
plot(0, pch="", xlim=c(1,366), ylim=c(0,4), xlab="Ordinal date", ylab="", yaxt="n", xaxt="n", bty="n")
axis(side=2, at=c(0.5:3.5), labels=c("Skå","Got","Hgl","Sth"), las=1, tick=FALSE)
axis(side=1, at=seq(0,360,60), pos=-0.3)
midmonths<-c(15,75,136,197,259,320)
monthnames<-c("Jan","Mar","May","Jul","Sep","Nov")
axis(side=1, at=midmonths, labels=monthnames, pos=0.23, tick=FALSE)

# What days fall below 60?
subset(raw_Sth, temp<6) %>% {rect(xleft=.$day-1, xright=.$day, ytop= 3.1 + ((.$year-2005)/10)*0.9,
                                  ybottom=3.1 + ((.$year-2005)/10)*0.9 - (1/10), col="#8f8fff", border=NA)}
subset(raw_Hgl, temp<6) %>% {rect(xleft=.$day-1, xright=.$day, ytop= 2.1 + ((.$year-2005)/10)*0.9,
                                  ybottom=2.1 + ((.$year-2005)/10)*0.9 - (1/10), col="#8f8fff", border=NA)}
subset(raw_Got, temp<6) %>% {rect(xleft=.$day-1, xright=.$day, ytop= 1.1 + ((.$year-2005)/10)*0.9,
                                  ybottom=1.1 + ((.$year-2005)/10)*0.9 - (1/10), col="#8f8fff", border=NA)}
subset(raw_Ssk, temp<6) %>% {rect(xleft=.$day-1, xright=.$day, ytop= 0.1 + ((.$year-2005)/10)*0.9,
                                  ybottom=0.1 + ((.$year-2005)/10)*0.9 - (1/10), col="#8f8fff", border=NA)}

# Point where temperature dips below 6 on average
winter <- c(291,288,306,316)
summer <- c(111,109,109,95)
rect(xleft=winter, xright=366, ybottom=3.1:0.1, ytop=4:1, col=alpha("#8f8fff",0.6), border=NA)
rect(xleft=0, xright=summer, ybottom=3.1:0.1, ytop=4:1, col=alpha("#8f8fff",0.6), border=NA)

# 1st and 3rd quartiles of inferred diapause start times
rect(xleft=diapause.25th, xright=diapause.75th, ybottom=c(3.1,2.1,1.1,1.1,0.1), ytop=c(4,3,2,2,1), col="grey95", border=NA)

# Average temperatures for each day of the year
heatmap <- magma(25)
rect(xleft=0:365, xright=1:366, ybottom=3, ytop=3.1, col=heatmap[round(ptempcurves$Sth$smoothcurve)+5], border=NA)
rect(xleft=0:365, xright=1:366, ybottom=2, ytop=2.1, col=heatmap[round(ptempcurves$Hgl$smoothcurve)+5], border=NA)
rect(xleft=0:365, xright=1:366, ybottom=1, ytop=1.1, col=heatmap[round(ptempcurves$Got$smoothcurve)+5], border=NA)
rect(xleft=0:365, xright=1:366, ybottom=0, ytop=0.1, col=heatmap[round(ptempcurves$Ssk$smoothcurve)+5], border=NA)

# Timing of flight peaks
timepoints <- seq(80, 300, length.out=512) # 512 sampling points; default of density() function
points(timepoints, 3.105 + (obs.curve$Sth/max(obs.curve$Sth))*0.85, type="l", lwd=2, col=alpha("black",0.5))
points(timepoints, 2.105 + (obs.curve$Hgl/max(obs.curve$Hgl))*0.85, type="l", lwd=2, col=alpha("black",0.5))
points(timepoints, 1.105 + (obs.curve$Got/max(obs.curve$Got))*0.85, type="l", lwd=2, col=alpha("black",0.5))
points(timepoints, 0.105 + (obs.curve$Ska/max(obs.curve$Ska))*0.85, type="l", lwd=2, col=alpha("black",0.5))

# Median inferred diapause start times
segments(x0=diapause.medians, x1=diapause.medians, y0=c(4,3,2,2,1), y1=c(3.1,2.1,1.1,1.1,0.1),
         col=c("black","black",heatmap[15],"black","black"), lty="dotted")

a.length <- 10
a.height <- 0.025
a.mid.y <- c(3.55,2.55,1.7,1.4,0.55)
arrow.x <- matrix( c(diapause.medians, winter[c(1,2,3,3,4)]-a.length, winter[c(1,2,3,3,4)]-a.length, winter[c(1,2,3,3,4)],
              winter[c(1,2,3,3,4)]-a.length, winter[c(1,2,3,3,4)]-a.length, diapause.medians), nrow=5)
arrow.y <- matrix( c( a.mid.y+a.height, a.mid.y+a.height, a.mid.y+a.height*2, a.mid.y,
              a.mid.y-a.height*2, a.mid.y-a.height, a.mid.y-a.height), nrow=5)
polygon(arrow.x[1,], arrow.y[1,], col="gray90")
polygon(arrow.x[2,], arrow.y[2,], col="gray90")
polygon(arrow.x[3,], arrow.y[3,], col=heatmap[19])
polygon(arrow.x[4,], arrow.y[4,], col="gray90")
polygon(arrow.x[5,], arrow.y[5,], col="gray90")



# Create map object
swedmap <- maps::map("world",xlim=c(11,19),ylim=c(55,60.5), fill=T, col="white", region=c("Sweden","Denmark"), plot=FALSE)
polyg <- map2SpatialPolygons(swedmap, IDs=swedmap$names,
                             proj4string=CRS("+proj=longlat +datum=WGS84"))
mapdata <- data.frame(seq_len(length(polyg)), row.names=names(polyg))
spdf <- SpatialPolygonsDataFrame(polyg, data=mapdata)
plot(spdf, ylim=c(55,61), lwd=0.3, ylab="Latitude", xlab="Longitude")
points(c(18.53,14.14,18.53,13.55), c(59.63,57.52,57.4,55.65), pch=".")
axis(side=1, at=seq(11,19,2), pos=54.6)
axis(side=2, at=seq(55,61,2), pos=10)


# Add heatmap legend
heatstrip <- seq(59.5,60.5,length.out=25)
rect(14.4,heatstrip,14,heatstrip+0.045, col=heatmap, border=NA)
lines(c(14.4,14.5), c(60.5, 60.5))
lines(c(14.4,14.5), c(60, 60))
lines(c(14.4,14.5), c(59.5, 59.5))
text(14.4, c(59.5, 60, 60.5), c("-4","8","20"), pos=4, cex=0.5)


dev.off()


