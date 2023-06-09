##Amanda Northrop
##University of Vermont, Department of Biology
##anorthro@gmail.com
##August 2017
##requires the following files: raw_data_long_format.csv, std_curve_data.csv, classes_by_treatment.csv, PIEcontrasts.csv, pilotdata2015.csv, and hysteresis_functions.R
##requires the ggplot2, dplyr, and benthos packages
##returns hysteresis indices for each replicate in each treatment
##returns hysteresis loop plots and statistical analyses and figures for publication

########################################################################################################################
################### READ IN THE RAW DATA IN LONG FORMAT AND SOURCE FUNCTIONS FILE ######################################
########################################################################################################################

data_long<-read.csv(file = "raw_data_long_format.csv",header=TRUE, na.strings = "NA")
std_curve_long<-read.csv(file = "standard_curve_data.csv",header=TRUE, stringsAsFactors = FALSE)  #note that the absorbance for 32 mg/ml column was only measured in two curves (values= 1.549 and 1.400)
source(file = "hysteresis_functions.R")
library(ggplot2)
library(dplyr)
library(benthos)

########################################################################################################################
################### FIT STANDARD CURVE #################################################################################
########################################################################################################################

#Reshape data to wide format
std_curve_wide<-reshape(data = std_curve_long, timevar = "Curve_Number", idvar = "Concentration_ugPERml",direction = "wide")

#Take an average of standard curve data to create a single curve, omitting NA values
std_curve_wide$Average<-(rowMeans(std_curve_wide[,-1],na.rm = TRUE))
std_curve_wide$SD<-apply(MARGIN = 1,X = std_curve_wide[,2:7],FUN = 'sd')

#Subtract the average value for the blank from the standard values
std_curve_wide[,2:7]<-std_curve_wide[,2:7]-std_curve_wide$Average[1]

#Fit a polynomial curve to the standards
linmod2<-lm(Concentration_ugPERml~poly(Average,2,raw=TRUE), std_curve_wide)


########################################################################################################################
########## PREDICT CONCENTRATION VALUES FROM ABSORBANCE VALUES (SUBTRACT STD CURVE AND CORRECT FOR DILUTION FACTOR)#######
##########################################################################################################################

#Add idvar column to long format data frame temporarily
data_long$SamplingTime<-paste(data_long$Day,data_long$Time.of.Day,sep="_")

#Reshape absorbance data to wide format to calculate BSA concentration from absorbance values 
ABS_wide<-reshape(data = data_long[which(data_long$Year==2015),],
                  drop = c("D.O.","Time.of.Day","Year","Day","BSA_Concentration"),
                  idvar = c("SamplingTime"),
                  timevar = "Pitcher.Number",
                  direction = "wide",
                  new.row.names = data_long$SamplingTime[1:72])

std_curve_wide[,2:7]<-std_curve_wide[,2:7]-std_curve_wide$Average[1]                                                                                                                                                                                                   

#Remove SamplingTime column and subtract std. curve from
ABS_wide<-ABS_wide[,-1]
ABS_wide_norm<-ABS_wide-std_curve_wide$Average[1] 

#Calculate BSA concentration values by passing absorbance values through linmod2 function in hysteresis_functions.R
BSA_wide<-as.data.frame(t(apply(ABS_wide_norm,1,linfun,
                                int_coeff = linmod2$coefficients[[1]],
                                poly1_coeff = linmod2$coefficients[[2]],
                                poly2_coeff = linmod2$coefficients[[3]])))
colnames(BSA_wide)<-sub("Absorbance","BSA_Calculated",colnames(BSA_wide))

#Multiply by dilution factor (x4 except for pitcher 75 which is x2) 
for(i in (1:dim(BSA_wide)[2])){
  if(colnames(BSA_wide)[i]=="BSA_Calculated.75"){
    BSA_wide[,i]<-(BSA_wide[,i]*2)/5}
  else{
    BSA_wide[,i]<-(BSA_wide[,i]*4)/5}
}

#Get the BSA data into the long format data frame "data_long"
data_long$BSA_Calculated<-c(reshape(BSA_wide,varying = names(BSA_wide),
                                    direction = "long")$BSA_Calculated,
                            data_long$Absorbance[which(data_long$Year==2016)])

#Make data_long into wide format
DO_wide<-reshape(data = data_long,
                 drop = c("Absorbance","Time.of.Day","Year","Day","BSA_Concentration","BSA_Calculated"),
                 idvar = c("SamplingTime"),
                 timevar = "Pitcher.Number",
                 direction = "wide",
                 new.row.names = data_long$SamplingTime[1:72])


Low<-c(DO_wide[,2],DO_wide[,3],DO_wide[,4],DO_wide[,5],DO_wide[,6],DO_wide[,7])
Intermediate<-c(DO_wide[,13],DO_wide[,14],DO_wide[,15],DO_wide[,16],DO_wide[,19],DO_wide[,20])
High<-c(DO_wide[,7],DO_wide[,8],DO_wide[,9],DO_wide[,10],DO_wide[,11])

                                    

BSA_wide<-cbind(BSA_wide,reshape(data = data_long,
                                 drop = c("Absorbance","Time.of.Day","Year","Day","BSA_Concentration","D.O."),
                                 idvar = c("SamplingTime"),
                                 timevar = "Pitcher.Number",
                                 direction = "wide",
                                 new.row.names = NULL)[,13:20])

#normalize the DO data
fud<-function(vector){
  max<-max(vector,na.rm=TRUE)
  min<-min(vector,na.rm=TRUE)
  y=(vector-min)/(max-min)
  return(y)
}

DO_wide[,2:20]<-apply(DO_wide[,2:20],MARGIN = 2,FUN = fud)
BSA_wide<-apply(BSA_wide,MARGIN = 2,FUN = fud)

Low<-c(DO_wide[,2],DO_wide[,3],DO_wide[,4],DO_wide[,5],DO_wide[,6],DO_wide[,7])
Intermediate<-c(DO_wide[,13],DO_wide[,14],DO_wide[,15],DO_wide[,16],DO_wide[,19],DO_wide[,20])
High<-c(DO_wide[,7],DO_wide[,8],DO_wide[,9],DO_wide[,10],DO_wide[,11])



new_data_long<-data.frame(SamplingTime=rep(DO_wide$SamplingTime,17),
                          D.O.=c(High,Low,Intermediate),
                          BSA_Calculated=c(c(BSA_wide[,1],BSA_wide[,2],BSA_wide[,3],BSA_wide[,4],BSA_wide[,5], BSA_wide[,6]),
                                           c(BSA_wide[,7],BSA_wide[,8],BSA_wide[,9],BSA_wide[,10],BSA_wide[,11]),
                                           c(BSA_wide[,12],BSA_wide[,13],BSA_wide[,14],BSA_wide[,15],BSA_wide[,18],BSA_wide[,19])),
                          BSA_Concentration=c(rep("High",432),rep("Low",360),rep("Intermediate",432)),
                          trip=rep(c(rep("Enrichment",8),rep("Recovery",64)),17))




########################################################################################################################
############################################### #FIGURE 3 ##############################################################
#####################################################################################################################
#Create panels for a multi-frame plot
P1<-
  ggplot(new_data_long[which(new_data_long$BSA_Concentration=="Low"),], 
         aes(x=BSA_Calculated, y=D.O., color=factor(trip))) +
  labs(x="\\n  ",y="Normalized DO")+
  ggtitle("Low")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title=element_text(size=16)) +
  theme(axis.title = element_text(size=16)) +
  theme(axis.text = element_text(size = 14, colour = "black")) +
  geom_point() +
  scale_color_manual(values=c("red","black")) +
  theme(legend.position = 'none')+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank(), axis.line = element_line())+
  stat_smooth(aes(y=D.O., fill=trip), 
              new_data_long[which(new_data_long$BSA_Concentration=="Low"),], 
              method="loess",na.rm = TRUE,level = .95)+
  scale_fill_manual(values = c("brown","black")) +
  ylim(0,1)

P2<-ggplot(new_data_long[which(new_data_long$BSA_Concentration=="Intermediate"),], 
           aes(x=BSA_Calculated, y=D.O., color=factor(trip))) +
  labs(x="Normalized BSA\\n",y="")+
  ggtitle("Intermediate")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title=element_text(size=16)) +
  theme(axis.title = element_text(size=16)) +
  theme(axis.text = element_text(size = 14, colour = "black")) +
  geom_point() +
  scale_color_manual(values=c("red","black")) +
  theme(legend.position = 'none')+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank(), axis.line = element_line())+
  stat_smooth(aes(y=D.O., fill=trip), 
              new_data_long[which(new_data_long$BSA_Concentration=="Intermediate"),], 
              method="loess",na.rm = TRUE,level = .95)+
  scale_fill_manual(values = c("brown","black")) +
  ylim(0,1)

P3<-ggplot(new_data_long[which(new_data_long$BSA_Concentration=="High"),], 
           aes(x=BSA_Calculated, y=D.O., color=factor(trip))) +
  labs(x="\\n  ",y="")+
  ggtitle("High")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title=element_text(size=16)) +
  theme(axis.title = element_text(size=16)) +
  theme(axis.text = element_text(size = 14, colour = "black")) +
  scale_color_manual(values=c("red","black")) +
  theme(legend.position = c(.8, .8))+
  theme(legend.title = element_blank())+
  theme(panel.background = element_blank())+
  theme(panel.border = element_blank(), axis.line = element_line())+
  stat_smooth(aes(y=D.O., fill=trip), 
              new_data_long[which(new_data_long$BSA_Concentration=="High"),], 
              method="loess",na.rm = TRUE,level = .95)+
  scale_fill_manual(values = c("brown","black"),guide = 'none') +
  geom_point() +
  ylim(0,1)

multiplot(P1,P2,P3,cols = 3)

#Write a function to normalize all the dataframes (From Lloyd et al. 2016)
normalize_fctn<-function(datfram){
  tempdat<-datfram
  for(i in 1:dim(datfram)[2]){ 
    for(y in 1:dim(datfram)[1]){ 
      tempdat[y,i]<-(datfram[y,i]-min(datfram[,i],na.rm = TRUE))/(max(datfram[,i],na.rm = TRUE)-min(datfram[,i],na.rm = TRUE))
    }
  }
  colnames(tempdat)<-colnames(datfram)
  return(tempdat)
}

tempnorm<-new_data_long[,2:3]
#Normalize the data frames (Lloyd et al. 2016)
new_data_long_norm<-normalize_fctn(tempnorm)
new_data_long_norm$trip<-new_data_long$trip
new_data_long_norm$Treatment<-new_data_long$BSA_Concentration


new_data_long_norm<-new_data_long_norm[order(new_data_long_norm$trip,decreasing = FALSE),]

hystmat<-matrix(nrow = (1/.05)+1,ncol = 2)
colnames(hystmat)<-colnames(new_data_long_norm)[1:2]

low_all_norm<-new_data_long_norm[which(new_data_long_norm$Treatment=="Low"),]
int_all_norm<-new_data_long_norm[which(new_data_long_norm$Treatment=="Intermediate"),]
high_all_norm<-new_data_long_norm[which(new_data_long_norm$Treatment=="High"),]



k=.05
BSAvals<-seq(min(low_all_norm$BSA_Calculated,na.rm = TRUE),max(low_all_norm$D.O.,na.rm=TRUE),k*(max(low_all_norm$BSA_Calculated,na.rm=TRUE)-min(low_all_norm$BSA_Calculated,na.rm=TRUE))+min(low_all_norm$BSA_Calculated,na.rm=TRUE))
loresFL<-loess(low_all_norm$D.O.[1:136]~low_all_norm$BSA_Calculated[1:136],span = .75,na.action=na.exclude)
loresRL<-loess(low_all_norm$D.O.[137:dim(low_all_norm)[1]]~low_all_norm$BSA_Calculated[137:dim(low_all_norm)[1]],span=.75,na.action=na.exclude)
valsFL<-predict(loresFL,BSAvals)
valsRL<-predict(loresRL,BSAvals)
hystvec<-valsRL-valsFL

low_himean<-mean(hystvec,na.rm = TRUE)

BSAvals<-seq(min(int_all_norm$BSA_Calculated,na.rm = TRUE),max(int_all_norm$D.O.,na.rm=TRUE),k*(max(int_all_norm$BSA_Calculated,na.rm=TRUE)-min(int_all_norm$BSA_Calculated,na.rm=TRUE))+min(int_all_norm$BSA_Calculated,na.rm=TRUE))
loresFL<-loess(int_all_norm$D.O.[1:136]~int_all_norm$BSA_Calculated[1:136],span = .75,na.action=na.exclude)
loresRL<-loess(int_all_norm$D.O.[137:dim(int_all_norm)[1]]~int_all_norm$BSA_Calculated[137:dim(int_all_norm)[1]],span=.75,na.action=na.exclude)
valsFL<-predict(loresFL,BSAvals)
valsRL<-predict(loresRL,BSAvals)
hystvec<-valsRL-valsFL

int_himean<-mean(hystvec,na.rm = TRUE)

BSAvals<-seq(min(high_all_norm$BSA_Calculated,na.rm = TRUE),max(high_all_norm$D.O.,na.rm=TRUE),k*(max(high_all_norm$BSA_Calculated,na.rm=TRUE)-min(high_all_norm$BSA_Calculated,na.rm=TRUE))+min(high_all_norm$BSA_Calculated,na.rm=TRUE))
loresFL<-loess(high_all_norm$D.O.[1:136]~high_all_norm$BSA_Calculated[1:136],span = .75,na.action=na.exclude)
loresRL<-loess(high_all_norm$D.O.[137:dim(high_all_norm)[1]]~high_all_norm$BSA_Calculated[137:dim(high_all_norm)[1]],span=.75,na.action=na.exclude)
valsFL<-predict(loresFL,BSAvals)
valsRL<-predict(loresRL,BSAvals)
hystvec<-valsRL-valsFL

high_himean<-mean(hystvec,na.rm = TRUE)

low_himean
int_himean
high_himean

