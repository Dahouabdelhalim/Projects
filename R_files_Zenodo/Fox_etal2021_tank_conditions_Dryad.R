######## Fox et al. 2021 - Differential resistance and acclimation of two coral species to chronic nutrient
########                   enrichment reflect life history traits 

####### Code file 1 of 4
####### Required data files: 1) Tank_irradiance.csv
#                            2) Tank_temperatures.csv

####### Summary: This file analyzes the light and temperature data over the course of the experiment.

####### Figures produced: Figs S1

library(dplyr) 
library(ggplot2)
library(lubridate)
library(splitstackshape)
library(scales)
library(cowplot)

### LICOR LIGHT MEASUREMENTS FOR DAILY IRRADIANCE TO CHAIN TANKS

### This is the file to import, plot, and ultimately deal with the LiCOR data from the CHAIN tank
##  The logger was first deployed from (0600-1900) on 10/27 until 11/19

## Missing data were back filled using correlations between our daily irrdiance records and those from HIMB
## Specifically we used a dual model approach to account for variation in afternoon light down at the tanks
## vs. up on the hill. Times from 6-12 were based off one linear regression and the afternoon times from 
## 1 to 19 were based off of a second. 

###  load file and add column names -- first column = 100 if data is bad so use as a check to fill w/NAs

D1<-read.csv(file="Tank_irradiance.csv",header=TRUE)

D1$PAR[D1$PAR <= 0] <- NA


#average daily range 
daily_range<-D1 %>% group_by(Day) %>% summarize(
                    N=sum(!is.na(PAR)),
                    Par_max=max(PAR,na.rm=TRUE),
                    Par_min=min(PAR,na.rm=TRUE))

#calculate diel range +/- SD.

mean(daily_range$Par_max) # 374.27
sd(daily_range$Par_max)   # +- 111.77

#create plot date
D1$Year<-2015
D1$Min<-0
D1$Sec<-0

D1$Date<-mdy(paste(D1$Month,D1$Day,D1$Year,sep="-"))
D1$Time<-paste(D1$Hour,D1$Min,D1$Sec,sep=":")
D1$test<-paste(D1$Date,D1$Time,sep=" ")
D1$plot<-as.POSIXct(D1$test) 

#define limits of experiment
#stretch to match light timeseries
start<-as.POSIXct(mdy("10-19-2015"))
end<-as.POSIXct(mdy("11-24-2015"))


#timeseries plot
a<-ggplot(D1,aes(x=plot,y=PAR))+geom_path()+
  scale_x_datetime(breaks=date_breaks("7 days"),labels=date_format("%b-%d"),limits=c(start,end))+
  theme_classic()+theme(text = element_text(size=14))+
  labs(x="",y='PAR ('*mu*'E'~m^-2~s^-1*')')+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))
a

#average irradiance by hour throughout the day 
hourly_means<-D1 %>% group_by(Hour) %>% summarize(
                       N=sum(!is.na(PAR)),
                       Par_hour_mean=mean(PAR,na.rm=TRUE),
                       sd_hour_mean=sd(PAR,na.rm=TRUE),
                        se_hour=sd_hour_mean/sqrt(N))      ##generates the mean PAR per hour from all days logged

#####################
#This plot shows the SD around each hour of the day for the CHAIN Experiment from 10.27-11.09
#Using this plot we can estimate the range of irridance at the peak of each day 

b<-ggplot(hourly_means,aes(x=Hour, y=Par_hour_mean))+
  geom_ribbon(aes(ymin=Par_hour_mean-sd_hour_mean,ymax=Par_hour_mean+sd_hour_mean),alpha=.5)+
  geom_line(size=.75)+
  scale_x_continuous(breaks=seq(6,19,1))+ylab("Mean PAR")+
  scale_y_continuous(limits=c(-10,600))+
  labs(x="Time (Hrs)",y='PAR ('*mu*'E'~m^-2~s^-1*')')+
  theme_classic()+theme(text = element_text(size=14))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))
b


############################################## Tank Temperature 

#missing first few days of temp because HOBO logger had filled at stopped recording
temp<-read.csv(file="Tank_temperatures.csv",header=T)

temp2<-cSplit(temp,"Date"," ")
temp2<-cSplit(temp2,"Date_2",":")
temp2<-cSplit(temp2,"Date_1","/")

colnames(temp2)<-c("Temp","Hour","Min","Month","Day","Year")
#create plotable date
temp2$Date<-mdy(paste(temp2$Month,temp2$Day,temp2$Year,sep="-"))
temp2$Time<-paste(temp2$Hour,temp2$Min,temp2$Sec,sep=":")
temp2$test<-paste(temp2$Date,temp2$Time,sep=" ")

temp2$plot<-as.POSIXct(temp2$test) 

#get rid of bad data points on 11.10
temp2<-subset(temp2,Temp>24.8)

#timeseries plot
c<-ggplot(temp2,aes(x=plot,y=Temp))+geom_path()+
  scale_x_datetime(breaks=date_breaks("7 days"),labels=date_format("%b-%d"),limits=c(start,end))+
  theme_classic()+theme(text = element_text(size=14))+
  labs(x="Date",y="Temperature (°C)")+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))
c

#hourly plot
temp3<-temp2 %>% group_by(Hour) %>% summarize(
                N=sum(!is.na(Temp)),
                Temp_mean=mean(Temp,na.rm=TRUE),
                sd=sd(Temp,na.rm=TRUE),
                se=sd/sqrt(N))   

d<-ggplot(temp3,aes(x=Hour,y=Temp_mean))+
  ylim(25,29)+
  labs(x="Time (Hrs)",y="Temperature (°C)")+
  geom_ribbon(aes(ymin=Temp_mean-sd,ymax=Temp_mean+sd),alpha=.5)+
  geom_line()
d

#combine the plots 
fig.S1<-plot_grid(a,b,c,d, labels = "AUTO", ncol = 2, align = 'v',axis='1',
          scale=c(1,1,1,1))
fig.S1
