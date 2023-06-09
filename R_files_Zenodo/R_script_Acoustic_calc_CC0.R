#Acoustic signal calculations - summarize data from exported Raven text files

##Before starting script: acoustic recordings extracted (Raven) and edited to a csv file with copied formulas. 
##Make sure SumStats.csv is in the working directory

#Two scripts are included; one calculates waveform variables, the other call variables

#Clear the workspace
rm(list=ls(all=TRUE))


##################################
#Calculations for WAVEFORMS

# Make sure to be in the appropriate working directory
#setwd("/Users/Lex/Documents/HYLA_ANDERSONII/ANLYS_LAB/3_Acoustic/NewCalls_Labeled_05_14_13/WF")

files <- list.files(path="./", pattern=".txt", all.files=T, full.names=T)

length(files)
head(files)

#Import the final data table and double check it's uploaded
Stat <- read.csv("WFStats.csv")
head(Stat)
dim(Stat)

#For loop to calculate mean, sd, cv for each call variable. 
for(i in 1:length(files)){
	W <- read.table(files[i],header=T, sep= "\\t", fill=T)
	wf <- (max(W$Selection) - 1)
	dur <- rep(1,wf) 
	
	for(j in 1:wf){
		dur[j]<- (W[(j+1),6]-W[j,6])
		#Oldest version
		#dur[j]<- (W[(j+1),4]-W[j,4])

	}

	Stat$File[i] <- files[i]
	Stat$numWF[i] <- wf
	Stat$meanDur[i] <- mean(dur)
	Stat$sdDur[i] <- sd(dur)
	Stat$cvDur <- Stat$sdDur/Stat$meanDur*100
	Stat$minDur[i] <- min(dur)
	Stat$maxDur[i] <- max(dur)
	Stat$sumDur[i] <- sum(dur)
	Stat$durDIVwf <- Stat$sumDur/Stat$numWF
	

	cutW <- W[1:wf+1, 7]


#Oldest version (don't use amp)
	Stat$meanAmp[i] <- mean(cutW)
	Stat$sdAmp[i] <- sd(cutW)
	Stat$cvAmp <- Stat$sdAmp/Stat$meanAmp*100
	Stat$minAmp[i] <- min(cutW)
	Stat$maxAmp[i] <- max(cutW)

}

head(Stat)

#After all rows done: 
write.csv(Stat,file="allWF.csv")

###After writing them to file:
###Add column with number of calls analyzed
###Take average of multiple waveforms values
###need to split the File name
##Add additional 8 columns to the file after "File" and 1 before
#Use "Data" -> "Text to Columns" and delimit by "_"
#Region, State, Field, PopCode, Month, Day, Year
#Find replace ".//" and correct any filter listing

##Convert MM.DD.YY to _ (find replace)
##Need to collapse the A, B, C portion to field number

##(use the skip feature to avoid copying the rest of the info)

#######################################

#Calculations for CALLS

#Clear the workspace
rm(list=ls(all=TRUE))

#Change the working directory manually (command D) or automatically

#Import the data and double check it's uploaded
Stat <- read.csv("SumStats2.csv")
head(Stat)
dim(Stat)


#List all files in a particular directory (make sure the calls end in calls or analysis)
files <- list.files(path="./", pattern="s.csv", all.files=T, full.names=T)


#Double check all files made it
length(files)
head(files)


for(i in 1:length(files)){
	D <- read.csv(files[i],header=TRUE)
	First<-1
	Stat$File[i] <- files[i]
	Stat$BeginFile[i] <- paste(D$BeginFile[First])
	Stat$BeginPath[i] <- paste(D$BeginPath[First])

	Calls<-max(D$Selection)
	First<-1
	Stat$numCalls[i] <- max(D$Selection)

	Stat$DeltaTime_mean[i] <- mean(D$DeltaTime)
	Stat$DeltaTime_sd[i] <- sd(D$DeltaTime)
	Stat$DeltaTime_cv <- Stat$DeltaTime_sd/Stat$DeltaTime_mean*100

	Stat$Dur90_mean[i] <- mean(D$Dur90)
	Stat$Dur90_sd[i] <- sd(D$Dur90)
	Stat$Dur90_cv <- Stat$Dur90_sd/Stat$Dur90_mean*100

	Stat$Freq5_mean[i] <- mean(D$Freq5)
	Stat$Freq5_sd[i] <- sd(D$Freq5)
	Stat$Freq5_cv <- Stat$Freq5_sd/Stat$Freq5_mean*100
	
	Stat$Q1Freq_mean[i] <- mean(D$Q1Freq)
	Stat$Q1Freq_sd[i] <- sd(D$Q1Freq)
	Stat$Q1Freq_cv <- Stat$Q1Freq_sd/Stat$Q1Freq_mean*100

	Stat$CenterFreq_mean[i] <- mean(D$CenterFreq)
	Stat$CenterFreq_sd[i] <- sd(D$CenterFreq)
	Stat$CenterFreq_cv <- Stat$CenterFreq_sd/Stat$CenterFreq_mean*100

	Stat$MaxFreq_mean[i] <- mean(D$MaxFreq)
	Stat$MaxFreq_sd[i] <- sd(D$MaxFreq)
	Stat$MaxFreq_cv <- Stat$MaxFreq_sd/Stat$MaxFreq_mean*100

	Stat$Q3Freq_mean[i] <- mean(D$Q3Freq)
	Stat$Q3Freq_sd[i] <- sd(D$Q3Freq)
	Stat$Q3Freq_cv <- Stat$Q3Freq_sd/Stat$Q3Freq_mean*100

	Stat$Freq95_mean[i] <- mean(D$Freq95)
	Stat$Freq95_sd[i] <- sd(D$Freq95)
	Stat$Freq95_cv <- Stat$Freq95_sd/Stat$Freq95_mean*100

	Stat$MaxPower_mean[i] <- mean(D$MaxPower)
	Stat$Energy_mean[i] <- mean(D$Energy)
	Stat$Length_mean[i] <- mean(D$Length)
	
	Stat$MaxAmp_mean[i] <- mean(D$MaxAmp)
	Stat$MaxAmp_sd[i] <- sd(D$MaxAmp)
	Stat$MaxAmp_cv <- Stat$MaxAmp_sd/Stat$MaxAmp_mean*100
	
	Stat$Interval_mean[i] <- mean(D$Interval,na.rm=T)
	Stat$Interval_sd[i] <- sd(D$Interval,na.rm=T)
	Stat$Interval_cv <- Stat$Interval_sd/Stat$Interval_mean*100

	Stat$Period1_mean[i] <- mean(D$Period1,na.rm=T)

	Stat$RiseTime_mean[i] <- mean(D$RiseTime)
	Stat$RiseTime_sd[i] <- sd(D$RiseTime)
	Stat$RiseTime_cv <- Stat$RiseTime_sd/Stat$RiseTime_mean*100

	Stat$FallTime_mean[i] <- mean(D$FallTime)
	Stat$FallTime_sd[i] <- sd(D$FallTime)
	Stat$FallTime_cv <- Stat$FallTime_sd/Stat$FallTime_mean*100

	Stat$Diff50Max_mean[i] <- mean(D$Diff50Max)

	Stat$Ratio5_25.5_50_mean[i] <- mean(D$Ratio5_25.5_50)
	Stat$Ratio5_25.5_50_sd[i] <- sd(D$Ratio5_25.5_50)
	Stat$Ratio5_25.5_50_cv <- Stat$Ratio5_25.5_50_sd/Stat$Ratio5_25.5_50_mean*100

	Stat$Ratio75_95.50_95_mean[i] <- mean(D$Ratio75_95.50_95)
	Stat$Ratio75_95.50_95_sd[i] <- sd(D$Ratio75_95.50_95)
	Stat$Ratio75_95.50_95_cv <- Stat$Ratio75_95.50_95_sd/Stat$Ratio75_95.50_95_mean*100

	Stat$Ratio_Rise.Fall_mean[i] <- mean(D$Ratio_Rise.Fall)
	Stat$Ratio_Rise.Fall_sd[i] <- sd(D$Ratio_Rise.Fall)
	Stat$Ratio_Rise.Fall_cv <- Stat$Ratio_Rise.Fall_sd/Stat$Ratio_Rise.Fall_mean*100

	Stat$Ratio_SlopeOn.Off_mean[i] <- mean(D$Ratio_SlopeOn.Off)
	Stat$Ratio_SlopeOn.Off_sd[i] <- sd(D$Ratio_SlopeOn.Off)
	Stat$Ratio_SlopeOn.Off_cv <- Stat$Ratio_SlopeOn.Off_sd/Stat$Ratio_SlopeOn.Off_mean*100

	Stat$Prop5_mean[i] <- mean(D$Prop5)
	Stat$Prop5_sd[i] <- sd(D$Prop5)
	Stat$Prop5_cv <- Stat$Prop5_sd/Stat$Prop5_mean*100
	
	Stat$Prop25_mean[i] <- mean(D$Prop25)
	Stat$Prop25_sd[i] <- sd(D$Prop25)
	Stat$Prop25_cv <- Stat$Prop25_sd/Stat$Prop25_mean*100

	Stat$Prop50_mean[i] <- mean(D$Prop50)
	Stat$Prop50_sd[i] <- sd(D$Prop50)
	Stat$Prop50_cv <- Stat$Prop50_sd/Stat$Prop50_mean*100

	Stat$PropMax_mean[i] <- mean(D$PropMax)
	Stat$PropMax_sd[i] <- sd(D$PropMax)
	Stat$PropMax_cv <- Stat$PropMax_sd/Stat$PropMax_mean*100

	Stat$Prop75_mean[i] <- mean(D$Prop75)
	Stat$Prop75_sd[i] <- sd(D$Prop75)
	Stat$Prop75_cv <- Stat$Prop75_sd/Stat$Prop75_mean*100

	Stat$Prop95_mean[i] <- mean(D$Prop95)
	Stat$Prop95_sd[i] <- sd(D$Prop95)
	Stat$Prop95_cv <- Stat$Prop95_sd/Stat$Prop95_mean*100

	Stat$BW90_mean[i] <- mean(D$BW90)
	Stat$BW90_sd[i] <- sd(D$BW90)
	Stat$BW90_cv <- Stat$BW90_sd/Stat$BW90_mean*100

	Stat$IQRBW_mean[i] <- mean(D$IQRBW)
	Stat$IQRBW_sd[i] <- sd(D$IQRBW)
	Stat$IQRBW_cv <- Stat$IQRBW_sd/Stat$IQRBW_mean*100
	
	Stat$TotTime[i] <- D$EndTime[Calls]-D$BeginTime[First]
	Stat$Rate[i] <- Calls/Stat$TotTime[i]
	Stat$DutyCycle <- Stat$Rate*Stat$DeltaTime_mean
	Stat$Period2[i] <- Stat$TotTime[i]/Calls
	Stat$DiffPeriod <- abs(Stat$Period2-Stat$Period1_mean)
	Stat$PropDiffPeriod <- Stat$DiffPeriod/Stat$Period2*100
}


#After all rows done: 
write.csv(Stat,file="allCalls.csv")


###After writing them to file, need to split the BeginPath
##Convert MM.DD.YY to _ (find replace)
##Need to collapse the A, B, C portion to field number
##Add additional 7 columns to the file
#Use "Data" -> "Text to Columns" to list the following 7 columns:
#Region, State, Field, PopCode, Month, Day, Year
##(use the skip feature to avoid copying the rest of the info)
