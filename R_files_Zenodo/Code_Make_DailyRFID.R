
## Load library
	library(plyr)

## Read data
	d0<-read.csv("Input_AllVisitsMerged.csv")
	
	d<-read.delim("Input_ties_by_row.txt")
		
## Make a new column with unique ID for each RFID/day combo
	d0$Visited.ubd<-paste(d0$UnitBox,d0$JDate,sep="_")

## Make a data frame with this unique code
	d0$Visited.ubd<-gsub("-","_",d0$Visited.ubd)
	lists<-as.data.frame(unique(d0$Visited.ubd))
	colnames(lists)<-"Visited.ubd"
	d2<-join(lists,d0,"Visited.ubd",match="first")
	d2<-d2[,c("Visited.ubd","JDate","Year","Offset","UnitBox","FemaleRFID","MaleRFID")]
	d2$fRFID_Day<-as.character(paste(d2$FemaleRFID,"_",d2$JDate,sep=""))
	d2$mRFID_Day<-as.character(paste(d2$MaleRFID,"_",d2$JDate,sep=""))
	
	d$Visited.ubd<-paste(d$VisitedUnitBox,d$JDate,sep="_")
	
	d$fRFID_Day<-as.character(paste(d$VisitorRFID,"_",d$JDate,sep=""))
	d$mRFID_Day<-as.character(paste(d$VisitorRFID,"_",d$JDate,sep=""))

## For each day, count number of male and female visits total and unique
	for(i in 1:nrow(d2)){
		sub.m<-subset(d,d$Visited.ubd==d2$Visited.ubd[i]&d$VisitorSex=="M")
		sub.f<-subset(d,d$Visited.ubd==d2$Visited.ubd[i]&d$VisitorSex=="F")
		sub.t<-subset(d,d$fRFID_Day==d2$fRFID_Day[i])
		sub.t.male<-subset(d,d$mRFID_Day==d2$mRFID_Day[i])
		d2$All.m[i]<-nrow(sub.m)
		d2$Uni.m[i]<-length(unique(sub.m$VisitorBand))
		d2$All.f[i]<-nrow(sub.f)
		d2$Uni.f[i]<-length(unique(sub.f$VisitorBand))
		d2$All.t[i]<-nrow(sub.t)
		d2$Uni.t[i]<-length(unique(sub.t$VisitedUnitBox))
		d2$All.t.m[i]<-nrow(sub.t.male)
		d2$All.t.m[i]<-length(unique(sub.t.male$VisitedUnitBox))
	}

	d2$All.v<-d2$All.m+d2$All.f
	d2$All.u<-d2$Uni.m+d2$Uni.f
	
## Write the output to a new file to be used in the other code
	write.csv(d2,"Input_Visits_By_Day.csv")
	
## Make a new object that is the averge number of visitors across all days 
	# with sd and sample size.
	
	d3<-as.data.frame(unique(d2$UnitBox))
	colnames(d3)<-"UnitBox"
	for(i in 1:nrow(d3)){
		sub<-subset(d2,d2$UnitBox==d3$UnitBox[i])
		d3$n[i]<-nrow(sub)
		d3$Avg.All.m[i]<-mean(sub$All.m)
		d3$SD.All.m[i]<-sd(sub$All.m)
		d3$Avg.Uni.m[i]<-mean(sub$Uni.m)
		d3$SD.Uni.m[i]<-sd(sub$Uni.m)
		d3$Avg.All.f[i]<-mean(sub$All.f)
		d3$SD.All.f[i]<-sd(sub$All.f)
		d3$Avg.Uni.f[i]<-mean(sub$Uni.f)
		d3$SD.Uni.f[i]<-sd(sub$Uni.f)
		d3$Avg.All.v[i]<-mean(sub$All.v)
		d3$SD.All.v[i]<-sd(sub$All.v)
		d3$Avg.Uni.v[i]<-mean(sub$All.u)
		d3$SD.Uni.v[i]<-sd(sub$All.u)
		d3$Avg.All.t[i]<-mean(sub$All.t)
		d3$SD.All.t[i]<-sd(sub$All.t)
		d3$Avg.Uni.t[i]<-mean(sub$Uni.t)
		d3$SD.Uni.t[i]<-sd(sub$Uni.t)
	}
	
## Write the output to a new file to be used in other script
	write.csv(d3,"Input_Avg_Daily_Visits.csv")