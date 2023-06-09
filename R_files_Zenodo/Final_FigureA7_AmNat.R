#Clear environment
rm(list=ls())

#Load necessary libraries
library(dplyr)
library(pracma)
library(lme4)
library(calibrate)
library(diptest)
library(MASS)
library(DHARMa)
library(glmmTMB)
library(fitdistrplus)
library(betareg)

#Set the working directory correctly on your own machine to read in the data.

#Read in the data
Prev_2009_2015<-read.csv("FieldPrevalence_2009-2015.csv",stringsAsFactors=F)
Prev_2016<-read.csv("FieldPrevalence_2016.csv",stringsAsFactors=F)

#Organize the data into the correct columns with correct column labels
Prev_neat_2009_2015=Prev_2009_2015[,c(2,3,1,5)]
Prev_neat_2016=cbind(rep(2016,dim(Prev_2016)[1]),Prev_2016[,c(1,2,4)])
names(Prev_neat_2009_2015)=c("Year","Jday","Lake","Infection prevalence")
names(Prev_neat_2016)=c("Year","Jday","Lake","Infection prevalence")

#All_prev is now an object that has year in the first column,
#Julian day in the second column, Lake in the third column,
#and infection prevalence, on a scale from 0-1, in the fourth column.
#Each row represents a single measurement of infection prevalence in
#a single lake.
All_prev=rbind(Prev_neat_2009_2015,Prev_neat_2016)


#Also read in the total phosphorus data
TP_2009=read.csv("2009_TP.csv")
TP_2010=read.csv("2010_TP.csv")
TP_2011=read.csv("2011_TP.csv")
TP_2012=read.csv("2012_TP.csv")
TP_2014=read.csv("2014_TP.csv")
TP_2015=read.csv("2015_TP.csv")

#Organizing lake name, year, jday, TP
All_TP1<-rbind(TP_2009,TP_2010,TP_2011,TP_2012,TP_2014,TP_2015)

Year=c(rep(2009,dim(TP_2009)[1]),rep(2010,dim(TP_2010)[1]),rep(2011,dim(TP_2011)[1]),
        rep(2012,dim(TP_2012)[1]),rep(2014,dim(TP_2014)[1]),rep(2015,dim(TP_2015)[1]))

All_TP_years=cbind(All_TP1,Year)

#Reorder columns to be structured like All_prev
All_TP<-cbind(All_TP_years$Year,All_TP_years$Jday,All_TP_years$Lake,All_TP_years$TP)

#Now for each lake year, calculate max prevalence and area under the curve
Lake_names=unique(c(All_prev[,3]))
Years=c(2009,2010,2011,2012,2013,2014,2015,2016)

#Create arrays of Lake_names and Years that can be easily indexed later
Lake_array<-array(NA,dim=c(length(Years),length(Lake_names)))
Year_array<-array(NA,dim=c(length(Years),length(Lake_names)))
for(i in 1:length(Lake_names)){
  Year_array[,i]<-Years  
}
for(j in 1:length(Years)){
  Lake_array[j,]<-Lake_names
}


#Make an array that holds key information for each lake-year.
#Each row corresponds to a different year, each column to a different lake

#Layer 1 is maximum prevalence, layer 2 is mean, layer 3 is area under the curve, layer 4 is area under the curve/time
Prev_summaries=array(NA,dim=c(length(Years),length(Lake_names),4))

#Here, layer 1 is total area under the curve of total phosphorus during the epidemic.
#Layer 2 is that total area divided by the duration of the epidemic.
TP_summaries=array(NA,dim=c(length(Years),length(Lake_names),2))

#This array captures the time range of the epidemic for each lake year.
#Measures for each lake year are only taken over this time range.
Time_range=array(NA,dim=c(length(Years),length(Lake_names),2))

#Set Julian days corresponding to the epidemic season for starting and ending
Start_time=238
End_time=312

#Loop through each lake-year combination
for(i in 1:length(Lake_names)){
  for(j in 1:length(Years)){
  #Skip if not sampled that lake year with at least two dates showing infection
  if(length(which(All_prev[,1]==Years[j]&All_prev[,3]==Lake_names[i]&All_prev[,2]<End_time&All_prev[,2]>Start_time&All_prev[,4]>0))<2) next
  
  #Get maximum prevalence  
  Prev_summaries[j,i,1]=max(All_prev[All_prev[,1]==Years[j]&All_prev[,3]==Lake_names[i]&All_prev[,2]<End_time&All_prev[,2]>Start_time,4],na.rm=T)
  
  #Get mean prevalence, weighting all dates the same
  Prev_summaries[j,i,2]=mean(All_prev[All_prev[,1]==Years[j]&All_prev[,3]==Lake_names[i]&All_prev[,2]<366&All_prev[,2]>0,4],na.rm=T)
  
  #Get the time range of the epidemic for this lake-year
  Time_range[j,i,]=c(max(All_prev[All_prev[,1]==Years[j]&All_prev[,3]==Lake_names[i]&All_prev[,4]>0,2],na.rm=T),min(All_prev[All_prev[,1]==Years[j]&All_prev[,3]==Lake_names[i]&All_prev[,4]>0,2],na.rm=T))
  
  #Get area under the curve prevalence for this lake-year
  Prev_summaries[j,i,3]=trapz(All_prev[All_prev[,1]==Years[j]&All_prev[,3]==Lake_names[i]&All_prev[,2]<=Time_range[j,i,1]&All_prev[,2]>=Time_range[j,i,2]&!is.na(All_prev[,4]),2],All_prev[All_prev[,1]==Years[j]&All_prev[,3]==Lake_names[i]&All_prev[,2]<=Time_range[j,i,1]&All_prev[,2]>=Time_range[j,i,2]&!is.na(All_prev[,4]),4])
  
  #Assuming that the time range makes sense, get the area under the curve
  #prevalence divided by the duration of the epidemic
  if(0<Time_range[j,i,2]&Time_range[j,i,1]<365) Prev_summaries[j,i,4]=ifelse(-diff(Time_range[j,i,])>0,Prev_summaries[j,i,3]/-diff(Time_range[j,i,]),NA)
  
  #Also find area under the curve prevalence for years excluding 2016
  if(j<length(Years)){
    #Get the TP summaries as well.
    TP_summaries[j,i,1]<-trapz(as.numeric(All_TP[All_TP[,1]==Years[j]&All_TP[,3]==Lake_names[i]&All_TP[,2]<=Time_range[j,i,1]&All_TP[,2]>=Time_range[j,i,2]&!is.na(All_TP[,4]),2]),as.numeric(All_TP[All_TP[,1]==Years[j]&All_TP[,3]==Lake_names[i]&All_TP[,2]<=Time_range[j,i,1]&All_TP[,2]>=Time_range[j,i,2]&!is.na(All_TP[,4]),4]))
    if(0<Time_range[j,i,2]&Time_range[j,i,1]<365) TP_summaries[j,i,2]=ifelse(-diff(Time_range[j,i,])>0,TP_summaries[j,i,1]/-diff(Time_range[j,i,]),NA)    
  }

  }
}


#set a minimum of 0 threshold for prevalence.
thresh=0

#Get mean phospohorus vectors for lake-years that had max prevalence exceeding a
#certain threshold, either 0.1 or 0.001.
field_TP_mean_0.1=c(TP_summaries[,,2])[which(c(Prev_summaries[,,1])>0.1)]
field_TP_mean_0.01=c(TP_summaries[,,2])[which(c(Prev_summaries[,,1])>0.01)]

#Get vectors of lakes and years corresponding to the field_TP_mean_0.1
#and field_TP_mean_0.01 vectors.
Lakes_final_0.1=c(Lake_array)[which(c(Prev_summaries[,,1])>0.1)]
Years_final_0.1=c(Year_array)[which(c(Prev_summaries[,,1])>0.1)]
Lakes_final_0.01=c(Lake_array)[which(c(Prev_summaries[,,1])>0.01)]
Years_final_0.01=c(Year_array)[which(c(Prev_summaries[,,1])>0.01)]

#Organize vectors of prevalence corresponding to meeting a condition of
#max prevalence greater than 0, 0.01, or 0.1
field_prev_max_0=c(Prev_summaries[,,1])[which(c(Prev_summaries[,,1])>0)]
field_prev_max_0.01=c(Prev_summaries[,,1])[which(c(Prev_summaries[,,1])>0.01)]
field_prev_mean_0.1=c(Prev_summaries[,,4])[which(c(Prev_summaries[,,1])>0.1)]

#Get a vector of mean prevalence corresponding to lake-years with max prevalence
#of at least 0.01.
field_prev_mean_0.01=c(Prev_summaries[,,4])[which(c(Prev_summaries[,,1])>0.01)]

#I exported some of these objects for use in another R script.
#saveRDS(field_prev_max_0,"field_prev_max_0.RDS")
#saveRDS(field_prev_mean_0.1,"field_prev_mean_0.1.RDS")
#saveRDS(field_prev_mean_0.01,"field_prev_mean_0.01.RDS")

#Get a vector of lake-years max prevalence where max prevalence was >0.1
field_prev_max_0.1=field_prev_max_0[which(field_prev_max_0>.1)]

#Make a beta regression of mean prevalence where a 0.01 threshold was met
#vs max prevalence meeting the same condition.
summary(simple_curve0.01<-betareg(field_prev_mean_0.01~field_prev_max_0.01))
#Check these diagnostic plots. They seem good enough.
plot(simple_curve0.01)

#Do the same for mean vs max prevalence with a 0.1 threshold instead of 0.01.
summary(simple_curve<-betareg(field_prev_mean_0.1~field_prev_max_0.1))
#Check diagnostic plots again.
plot(simple_curve)

#Make some fake data for the purpose of plotting the beta regression curve.
example_max<-linspace(0,1,1e4)
new_max=data.frame(field_prev_max_0.1=example_max)
new_mean<-predict(simple_curve,new_max,type="response")

#DO a beta regression of max prevalence with a 0.1 cutoff threshold
#upon mean total phosphorus with lake and year as random effects.
beta_random_max<-glmmTMB(field_prev_max_0.1~field_TP_mean_0.1+(1|Lakes_final_0.1)+(1|Years_final_0.1),family=beta_family())
summary(beta_random_max)
diagnostic<-simulateResiduals(beta_random_max)
#QQ plot confirms good beta distribution. Resiual vs predicted indicates no problems with difference in variance or problem with linearity assumption
plot(diagnostic)

#Check for a lower prevalence threshold
beta_random_max0.01<-glmmTMB(field_prev_max_0.01~field_TP_mean_0.01+(1|Lakes_final_0.01)+(1|Years_final_0.01),family=beta_family())
summary(beta_random_max0.01)
diagnostic<-simulateResiduals(beta_random_max0.01)
#QQ plot confirms good beta distribution. Resiual vs predicted indicates no problems with difference in variance or problem with linearity assumption
plot(diagnostic)

#Make fake data for plotting purposes.
example_TP<-linspace(0,40,1e4)
new_TP=data.frame(field_TP_mean_0.1=example_TP,Lakes_final_0.1=rep(NA,1e4),Years_final_0.1=rep(NA,1e4))
new_prevs<-predict(beta_random_max,new_TP,type="response")

#Make Figure A7 of the manuscript
{
  png("C:/Users/Jason/OneDrive - University of Pittsburgh/Documents/Thesis Chapter 2/Appendix_figure_phosphorus.png",width=8,height=5,units="in",res=600)
  cex_smallest_text=.5
  cex_minor_text=1
  cex_major_text=1
  lwd_minor=1
  lwd_major=1.5
  par(mfrow=c(1,2),mar=c(5,4,1,.5),mgp=c(1.8,.4,0),cex.lab=cex_major_text,oma=c(0,0,0,0),cex.axis=cex_minor_text,bg="white",cex=1)
  
  plot(field_TP_mean_0.1,field_prev_max_0.1,pch=20,xlab="",ylab=expression("Max prevalence:"~italic(p)~"(unitless)"),cex=cex_minor_text)
  mtext(side=1,expression(atop("Average total phosphorus:","("~mu~"g P"~"L"^"-1"~")")),padj=1.5,cex=1)
  points(example_TP,new_prevs,lwd=lwd_major,type="l")
  text(27.5,.4,"P = 0.0174",cex=cex_minor_text)
  text((par("usr")[2]-par("usr")[1])*.1+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_minor_text)
  
  
  plot(field_prev_max_0.1,field_prev_mean_0.1,xlab="",ylab=expression("Average prevalence:"~italic(p)~"(unitless)"),cex=cex_minor_text,pch=20)
  mtext(side=1,expression(atop("Max prevalence:",italic(p)~"(unitless)")),padj=1.5,cex=1)
  points(example_max,new_mean,type="l",lwd=lwd_major)  
  text(.22,.18,expression("Pseudo R"^"2"~"= 0.506"),cex=cex_minor_text)
  text((par("usr")[2]-par("usr")[1])*.1+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_minor_text)
  dev.off()
}

#Finish fixing up numbers on the panels and put this into the main text.
#Check other numbers in the main text. Get next iteration of GA running.