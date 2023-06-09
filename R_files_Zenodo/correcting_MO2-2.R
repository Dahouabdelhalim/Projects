setwd("~/Documents/PhD Project/Respirometry/Respirometry")

respo<-read.csv("C26_32data.csv") #read in respo data based on what trial you want
bg_start<-read.table("background_startC26_32.txt", skip=16) # read in start background data based on what trial you want
bg_end<-read.table("background_endC26_32.txt", skip=16) #read in end background data based on what trial you want
bg_start<-bg_start[,c(1:4,8)]
bg_end<-bg_end[,c(1:4,8)]
names(bg_start)<-c("Date","Time","Seconds","Oxygen","Temp")
names(bg_end)<-c("Date","Time","Seconds","Oxygen","Temp")

bg_start$Seconds<-round(bg_start$Seconds, digits=0)
bg_end$Seconds<-round(bg_end$Seconds,digits=0)

plot(Oxygen~Seconds,data=bg_start)
plot(Temp~Seconds,data=bg_start)

head(bg_start)
tail(bg_start)

#isolate 3 hours of decreasing O2 and replot bg_start
plot(Oxygen~Seconds,data=bg_start)



plot(Oxygen~Seconds,data=bg_end)
plot(Temp~Seconds,data=bg_end)

head(bg_end)
tail(bg_end)


#isolate 3 hours of decreasing O2 and replot bg_end
plot(Oxygen~Seconds,data=bg_end)



#put respo data in correct format
respo$Date_Time<-as.factor(paste(respo$Date, respo$time, sep=" "))
respo$Date_Time<-as.POSIXlt(respo$Date_Time, format="%m/%d/%y %H:%M:%S")

respo$Date_Timen<-as.numeric(respo$Date_Time)
respo$Time_Diff<-respo$Date_Timen-respo$Date_Timen[1]

#regression on background data
mstart<-lm(Oxygen~Seconds, data=bg_start)
slope1<-coef(mstart)[2]
mend<-lm(Oxygen~Seconds, data=bg_end)
slope2<-coef(mend)[2]

#make regression between two slopes of background
time<-c(0,respo$Time_Diff[nrow(respo)])
slopes<-c(slope1, slope2)

bg<-lm(slopes~time)
slope3<-coef(bg)[2]
intercept3<-coef(bg)[1]

y<-(slope3*respo$Time_Diff) + intercept3
plot(y~respo$Time_Diff)

#subtracting background from slope to get adjusted slopes for each point 
respo$newslope<-respo$slope-y
#turn umol O2/l/s into umol/kg/hr
respo$newMO2<--1*respo$newslope*(respo$chamber.volume-respo$animal.mass)*3600/respo$animal.mass
#turn umol O2/kg/hr into mg O2/kg/hr
respo$finalMO2<-respo$newMO2*32/1000

#write new csv
setwd("~/Documents/PhD Project/Respirometry/Respirometry/Cobia Corrected Data")








#check confidence intervals of slope
test<-lm(Oxygen~Seconds, data=bg_start)
summary(test)
confint(test,"Seconds",level=0.95)

new <- data.frame(Seconds = bg_start$Seconds)
p<-predict(test, newdata=new,se.fit=TRUE)
plot(Oxygen~Seconds,data=bg_start)
lines(bg_start$Seconds,p$fit)
