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


########################################################################################################################
############################################ PLOT 02 TRACES ##############################################################
########################################################################################################################

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

#Calculate the mean and standard deviation of the O2 time series for each treatment
DO_wide$Mean05<-rowMeans(DO_wide[,8:12],na.rm = TRUE)
DO_wide$Mean20<-rowMeans(DO_wide[,c(13:16,19:20)],na.rm = TRUE)
DO_wide$Mean50<-rowMeans(DO_wide[,2:7],na.rm = TRUE)
DO_wide$SD05<-apply(DO_wide[,8:12],1,FUN = sd,na.rm = TRUE)
DO_wide$SD20<-apply(DO_wide[,c(13:16,19:20)],1,FUN = sd,na.rm = TRUE)
DO_wide$SD50<-apply(DO_wide[,2:7],1,FUN = sd,na.rm = TRUE)


########################################################################################################################
############################################### #FIGURE 2 ##############################################################
########################################################################################################################

par(mfrow=c(3,1),mgp=c(3,1,0),mar=c(6.0,5.0,1.0,4.0))

plot(DO_wide$Mean05~seq(0,35.5,by=.5), type="n",ylim=c(0,25),
     xlab="",ylab="Dissolved Oxygen (%)", cex = 1.5, 
     cex.axis = 1.5, cex.lab = 1.5)
rect(-1.3,-.25,36.8,2, col="gray90",border="gray90")
points(DO_wide$Mean05~seq(0,35.5,by= .5),pch=19, cex=1,col="black")
arrows(x0 = c(0,1,2,3),y1 = c(20,20,20,20),x1=c(0,1,2,3),y0=c(25,25,25,25),length = .1)
arrows(seq(0,35.5,by=.5), DO_wide$Mean05-DO_wide$SD05, seq(0,35.5,by=.5), 
       DO_wide$Mean05+DO_wide$SD05, col="blue",length=0.05, angle=90, code=3)
text(33,24,labels="Low",cex=1.5)
text(33, .9, labels = "n=5",cex=1.5)

plot(DO_wide$Mean20~seq(0,35.5,by=.5), type="n",ylim=c(0,25),
     xlab="",ylab="Dissolved Oxygen (%)", cex = 1.5, 
     cex.axis = 1.5, cex.lab = 1.5)
rect(-1.3,-.25,36.8,2, col="gray90",border="gray90")
points(DO_wide$Mean20~seq(0,35.5,by=.5),pch=19, cex=1,col="black")
arrows(x0 = c(0,1,2,3),y1 = c(20,20,20,20),x1=c(0,1,2,3),y0=c(25,25,25,25),length = .1)
arrows(seq(0,35.5,by=.5), DO_wide$Mean20-DO_wide$SD20, seq(0,35.5,by=.5), 
       DO_wide$Mean20+DO_wide$SD20, col="green",length=0.05, angle=90, code=3)
text(33,24,labels="Intermediate",cex=1.5)
text(33, .9, labels = "n=6",cex=1.5)

plot(DO_wide$Mean50~seq(0,35.5,by=.5), type="n",ylim=c(0,25),
     xlab="Days",ylab="Dissolved Oxygen (%)", cex = 1.5, 
     cex.axis = 1.5, cex.lab = 1.5)
rect(-1.3,-.25,36.8,2, col="gray90",border="gray90")
points(DO_wide$Mean50~seq(0,35.5,by= .5),pch=19, cex=1,col="black")
arrows(x0 = c(0,1,2,3),y1 = c(20,20,20,20),x1=c(0,1,2,3),y0=c(25,25,25,25),length = .1)
arrows(seq(0,35.5,by=.5), DO_wide$Mean50-DO_wide$SD50, seq(0,35.5,by=.5), 
       DO_wide$Mean50+DO_wide$SD50, col="brown",length=0.05, angle=90, code=3)
text(33,23,labels="High",cex=1.5)
text(33, .9, labels = "n=6",cex=1.5)

########################################################################################################################
############################################### #FIGURE S3 ##############################################################
########################################################################################################################

#Plot the individual DO traces for each replicate in eaach treatment
par(mfrow=c(3,6),mar=c(4,4,1.2,1.2))

for(i in 8:12){
  plot(DO_wide[,i]~seq(0,35.5,by=.5),
       col = "blue",
       ylab = "Dissolved Oxygen (%)", 
       xlab = "Days\\n",
       cex = 1,
       cex.lab=1,cex.axis=1,
       xlim = c(0,35),
       ylim = c(0,25),
       main = paste("Pitcher",strsplit(colnames(DO_wide)[2:20], split = "[.][.]")[[i-1]][2],sep=" "))
  abline(h=19,lty=2,col="grey",lwd=1.5)
}

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =c('Low', 'Intermediate', 'High'), 
       pch=16, pt.cex=1, cex=1, bty='n',
       col = c('blue', 'green', 'brown'))

for(i in c(13:16,19:20)){
  plot(DO_wide[,i]~seq(0,35.5,by=.5),
       col = "green",
       ylab = "Dissolved Oxygen (%)", 
       xlab = "Days\\n",
       cex=1,
       cex.lab=1,cex.axis=1,
       xlim = c(0,35),
       ylim = c(0,30),
       main = paste("Pitcher",strsplit(colnames(DO_wide)[2:20], split = "[.][.]")[[i-1]][2],sep=" "))
  abline(h=19,lty=2,col="grey",lwd=1.5)
}

for(i in 2:7){
  plot(DO_wide[,i]~seq(0,35.5,by=.5),
       col = "brown",
       ylab ="Dissolved Oxygen (%)", 
       xlab = "Days\\n",
       cex=1,
       cex.lab=1,cex.axis=1,
       xlim = c(0,35),
       ylim = c(0,25),
       main = paste("Pitcher",strsplit(colnames(DO_wide)[2:20], split = "[.][.]")[[i-1]][2],sep=" "))
  abline(h=19,lty=2,col="grey",lwd=1.5)
}



########################################################################################################################
############################################### PLOT STATE-SPACE FIGURES ###############################################
########################################################################################################################

#Add the 2016 BSA data to the wide format data frame for BSA concentrations
BSA_wide<-cbind(BSA_wide,reshape(data = data_long,
                                 drop = c("Absorbance","Time.of.Day","Year","Day","BSA_Concentration","D.O."),
                                 idvar = c("SamplingTime"),
                                 timevar = "Pitcher.Number",
                                 direction = "wide",
                                 new.row.names = NULL)[,13:20])

state_space_df<-data.frame(BSA_Calculated = c(rowMeans(BSA_wide[1:6],na.rm = TRUE),
                                              rowMeans(BSA_wide[7:11],na.rm = TRUE),
                                              rowMeans(BSA_wide[c(12:15,18:19)],na.rm = TRUE)),
                           D.O. = c(rowMeans(DO_wide[2:7],na.rm = TRUE),
                                    rowMeans(DO_wide[8:12],na.rm = TRUE),
                                    rowMeans(DO_wide[c(13:16,19:20)],na.rm = TRUE)),
                           BSA_Concentration  = c(rep("5",72),rep(".5",72),rep("2",72)),
                           trip = rep(c(rep("Enrichment",8),rep("Recovery",64)),3))


########################################################################################################################
############################################### #FIGURE 3 ##############################################################
########################################################################################################################

#Create panels for a multi-frame plot
P1<-
  ggplot(state_space_df[which(state_space_df$BSA_Concentration==".5"),], 
         aes(x=BSA_Calculated, y=D.O., color=factor(trip))) +
  labs(x="\\n  ",y="Dissolved Oxygen (%)")+
  xlim(0,25)+
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
              state_space_df[which(state_space_df$BSA_Concentration==".5"),], 
              method="loess",na.rm = TRUE,level = .95)+
  scale_fill_manual(values = c("brown","black")) +
  ylim(-8,25)

P2<-ggplot(state_space_df[which(state_space_df$BSA_Concentration=="2"),], 
           aes(x=BSA_Calculated, y=D.O., color=factor(trip))) +
  labs(x="\\nBSA (mg/mL)",y="")+
  xlim(0,25)+
  
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
              state_space_df[which(state_space_df$BSA_Concentration=="2"),], 
              method="loess",na.rm = TRUE,level = .95)+
  scale_fill_manual(values = c("brown","black")) +
  ylim(-8,25)

P3<-ggplot(state_space_df[which(state_space_df$BSA_Concentration=="5"),], 
           aes(x=BSA_Calculated, y=D.O., color=factor(trip))) +
  labs(x="\\n  ",y="")+
  xlim(0,25)+
  
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
              state_space_df[which(state_space_df$BSA_Concentration=="5"),], 
              method="loess",na.rm = TRUE,level = .95)+
  scale_fill_manual(values = c("brown","black"),guide = 'none') +
  geom_point() +
  ylim(-8,25)

multiplot(P1,P2,P3,cols = 3)

########################################################################################################################
############################################### Fig S4  ##############################################################
########################################################################################################################

#Calculate the mean and standard deviation of the BSA time series for each treatment
BSA_wide$Mean05<-rowMeans(BSA_wide[,7:11],na.rm = TRUE)
BSA_wide$Mean20<-rowMeans(BSA_wide[,c(12:15,18:19)],na.rm = TRUE)
BSA_wide$Mean50<-rowMeans(BSA_wide[,1:6],na.rm = TRUE)
BSA_wide$SD05<-apply(BSA_wide[,7:11],1,FUN = sd,na.rm = TRUE)
BSA_wide$SD20<-apply(BSA_wide[,c(12:15,18:19)],1,FUN = sd,na.rm = TRUE)
BSA_wide$SD50<-apply(BSA_wide[,1:6],1,FUN = sd,na.rm = TRUE)

#plot
par(mar=c(4.2,4.2,2.2,2.2),mfrow=c(1,1))
plot(BSA_wide$Mean05,xaxt="n",ylim=c(0,30),cex=1,col="blue",pch=16,
     ylab="BSA (mg/ml)",xlab="\\nTime (days)")
axis(1, at = seq(0, dim(BSA_wide)[1], by = 1),
     label = seq(0,36,by =.5), las=2,cex.axis=.75)
arrows(x0 = 1:dim(BSA_wide)[1],y0 = BSA_wide$Mean05-BSA_wide$SD05,
       x1 = 1:dim(BSA_wide)[1],y1 = BSA_wide$Mean05+BSA_wide$SD05,lwd=1.5,
       col = adjustcolor("blue",alpha.f = .40),angle = 90,length = .05,code = 3)

points(BSA_wide$Mean50,col="brown",pch=16,cex=1)
arrows(x0 = 1:dim(BSA_wide)[1],y0 = BSA_wide$Mean50-BSA_wide$SD50,
       x1 = 1:dim(BSA_wide)[1],y1 = BSA_wide$Mean50+BSA_wide$SD50,lwd=1.5,
       col = adjustcolor("brown",alpha.f = .40),angle = 90,length = .05,code = 3)


points(BSA_wide$Mean20,col="green",pch=16,cex=1)
arrows(x0 = 1:dim(BSA_wide)[1],y0 = BSA_wide$Mean20-BSA_wide$SD20,
       x1 = 1:dim(BSA_wide)[1],y1 = BSA_wide$Mean20+BSA_wide$SD20,lwd=1.5,
       col = adjustcolor("green",alpha.f = .60),angle = 90,length = .05,code = 3)

legend(x="topright", legend=c("Low","Intermediate","High"),col=c("blue","green","brown"),
       pch=c(16,16,16),bty = "n")


########################################################################################################################
############################################## HYSTERESIS INIDCES #####################################################
########################################################################################################################


#Write a function to normalize all the dataframes (From Lloyd et al. 2016)
normalize_fctn<-function(datfram){
  tempdat<-datfram[,-1]
  for(i in 1:dim(datfram)[2]){ 
    for(y in 1:dim(datfram)[1]){ 
      tempdat[y,i]<-(datfram[y,i]-min(datfram[,i],na.rm = TRUE))/(max(datfram[,i],na.rm = TRUE)-min(datfram[,i],na.rm = TRUE))
    }
  }
  colnames(tempdat)<-colnames(datfram)
  return(tempdat)
}

#Normalize the data frames (Lloyd et al. 2016)
BSA_wide_norm<-normalize_fctn(BSA_wide[,1:19])

#Remove unnecessary columns from DO_wide
DO_wide_clean<-DO_wide[,2:20]
DO_wide_norm<-normalize_fctn(DO_wide_clean)



#For each of the replicates, calculate the intervals on which you're going to measure HI 
#k is the width of the intervals
#Then calculate the hysteresis index - Trlnorm-Tflnorm
#Then calculate the mean HI over the interval. 
calculateHI<-function(do_datfram,conc_datfram,k){
  hystmat<-matrix(nrow = (1/k)+1,ncol = dim(conc_datfram)[2])
  colnames(hystmat)<-colnames(do_datfram)
  for(i in 1:dim(conc_datfram)[2]){
    #calculate the range of values of BSA over which to determine HI (using 5% intervals per Lloyd et al. 2016)
    BSAvals<-seq(min(conc_datfram[,i],na.rm = TRUE),max(conc_datfram[,i],na.rm=TRUE),k*(max(conc_datfram[,i],na.rm=TRUE)-min(conc_datfram[,i],na.rm=TRUE))+min(conc_datfram[,i],na.rm=TRUE))
    loresFL<-loess(do_datfram[1:10,i]~conc_datfram[1:10,i],span = .75,na.action=na.exclude)
    loresRL<-loess(do_datfram[11:dim(do_datfram)[1],i]~conc_datfram[11:dim(conc_datfram)[1],i],span=.75,na.action=na.exclude)
    valsFL<-predict(loresFL,BSAvals)
    valsRL<-predict(loresRL,BSAvals)
    hystmat[,i]<-valsRL-valsFL
  }
  return(hystmat)
}

HIqi<-calculateHI(DO_wide_norm, BSA_wide_norm,.05)

MeanLow<-rowMeans(na.rm = TRUE,HIqi[,c(7:11)])[1:20]
MeanHigh<-rowMeans(na.rm = TRUE,HIqi[,c(1:6)])[2:19]
MeanInt<-rowMeans(na.rm=TRUE,HIqi[,c(12:15,18:19)])[1:20]


#calculate mean HImean
HIMEANS<-colMeans(HIqi,na.rm = TRUE)


########################################################################################################################
##ANOVA test to see if mean HI differs among treatments
########################################################################################################################
HI_anova<-data.frame(Treatment = c(rep("Low",5), rep("Intermediate",6) ,rep("High",6)),
                     HI=c(HIMEANS[7:11],HIMEANS[c(12:15,18:19)],HIMEANS[1:6]),
                     Color = c(rep("blue",5),rep("green",6),rep("brown",6)))

aov_res<-aov(HI ~ Treatment, data = HI_anova)       
summary(aov_res)
TukeyHSD(aov_res)

########################################################################################################################
##################################### FIGURE 4 ########################################################
########################################################################################################################

par(mar=c(5.2,5.2,2.2,2.2),mfrow=c(1,1))
HIfram<-data.frame(Mean=tapply(HI_anova$HI, HI_anova$Treatment, mean),
                   SD=tapply(HI_anova$HI, HI_anova$Treatment, sd))
HIfram<-HIfram[c(3,2,1),]

o <- ordered(HI_anova$Treatment, levels = c("Low", "Intermediate", "High"))
boxplot(HI~o,data = HI_anova,outline = FALSE,names = c("Low","Intermediate","High"),
        boxlty=0,whisklty=0,staplelty=0,ylim=c(-.5,.5),
        xlab="Enrichment Treatment",ylab=expression("HI"["MEAN"]),cex.axis=1.5,cex=2,cex.lab=2)



stripchart(HI~o, data = HI_anova, 
           vertical = TRUE,
           cex=1.5,
           pch = 19, 
           col= c("blue","green","brown"),
           add = TRUE,method="jitter")
segments(x0 = 0.38,x1 = 3.62,y0 = 0,y1 = 0,lty = 2,col="grey",lwd=1.5)


########################################################################################################################
##################################### FIGURE S5 ###################################################################################
########################################################################################################################
#plot the individual hysteresis plots for each replicate
par(mar=c(5,4,2,2))
#add the enrichment phase to the data_long data frame
data_long$Phase<-as.factor(rep(c(rep("Enrichment",8),rep("Recovery",64)),19))

#clean the bsa data
BSA_wide_clean<-BSA_wide[,1:19]

plotlist<-vector("list",20)

data_long$ColorBSA<-NA
#add custom colors for plot to data_long
data_long$ColorBSA[which(data_long$BSA_Concentration==5)]<-"brown"
data_long$ColorBSA[which(data_long$BSA_Concentration==2)]<-"green"
data_long$ColorBSA[which(data_long$BSA_Concentration==.5)]<-"blue"

pitchernum<-unique(data_long$Pitcher.Number)

for(i in 1:6){
  plotlist[[i]]<-ggplot(data_long[which(data_long$Pitcher.Number==pitchernum[i]),], 
                        aes(x=BSA_Calculated, y=D.O., color=factor(Phase))) +
    labs(x="BSA (mg/ml)",y="Dissolved Oxygen (%)")+
    ggtitle(paste("Pitcher ",pitchernum[i]))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title=element_text(size=10)) +
    theme(axis.title = element_text(size=10)) +
    theme(axis.text = element_text(size = 10, colour = "black")) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    geom_point() +
    scale_color_manual(values=c("brown","black")) +
    theme(legend.position = 'none')+
    theme(panel.background = element_blank())+
    theme(panel.border = element_blank(), axis.line = element_line())+
    stat_smooth(aes(y=D.O., fill=Phase), 
                data_long[which(data_long$Pitcher.Number==pitchernum[i]),], 
                method="loess",na.rm = TRUE,level = .95)+
    scale_fill_manual(values = c("brown","black")) +
    ylim(-8,25)
}

for(i in 7:11){
  plotlist[[i]]<-ggplot(data_long[which(data_long$Pitcher.Number==pitchernum[i]),], 
                        aes(x=BSA_Calculated, y=D.O., color=factor(Phase))) +
    labs(x="BSA (mg/ml)",y="Dissolved Oxygen (%)")+
    ggtitle(paste("Pitcher ",pitchernum[i]))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title=element_text(size=10)) +
    theme(axis.title = element_text(size=10)) +
    theme(axis.text = element_text(size = 10, colour = "black")) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    geom_point() +
    scale_color_manual(values=c("blue","black")) +
    theme(legend.position = 'none')+
    theme(panel.background = element_blank())+
    theme(panel.border = element_blank(), axis.line = element_line())+
    stat_smooth(aes(y=D.O., fill=Phase), 
                data_long[which(data_long$Pitcher.Number==pitchernum[i]),], 
                method="loess",na.rm = TRUE,level = .95)+
    scale_fill_manual(values = c("blue","black")) +
    ylim(-8,25)
}

for(i in c(12:15,18:19)){
  plotlist[[i]]<-ggplot(data_long[which(data_long$Pitcher.Number==pitchernum[i]),], 
                        aes(x=BSA_Calculated, y=D.O., color=factor(Phase))) +
    labs(x="BSA (mg/ml)",y="Dissolved Oxygen (%)")+
    ggtitle(paste("Pitcher ",pitchernum[i]))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title=element_text(size=10)) +
    theme(axis.title = element_text(size=10)) +
    theme(axis.text = element_text(size = 10, colour = "black")) +
    theme(plot.margin = unit(c(0,0,0,0), "pt"))+
    geom_point() +
    scale_color_manual(values=c("green","black")) +
    theme(legend.position = 'none')+
    theme(panel.background = element_blank())+
    theme(panel.border = element_blank(), axis.line = element_line())+
    stat_smooth(aes(y=D.O., fill=Phase), 
                data_long[which(data_long$Pitcher.Number==pitchernum[i]),], 
                method="loess",na.rm = TRUE,level = .95)+
    scale_fill_manual(values = c("green","black")) +
    ylim(-8,25)
}

# testlist<-list(plotlist[[7]],plotlist[[12]],plotlist[[17]],
#                plotlist[[3]],plotlist[[8]],plotlist[[13]],
#                plotlist[[18]],plotlist[[4]],plotlist[[9]],
#                plotlist[[14]],plotlist[[19]],plotlist[[5]],
#                plotlist[[10]],plotlist[[15]],plotlist[[1]],
#                plotlist[[6]],plotlist[[11]],plotlist[[16]],
#                plotlist[[2]])

testlist<-list(plotlist[[7]],plotlist[[12]],plotlist[[1]],
               plotlist[[8]],plotlist[[13]],plotlist[[2]],
               plotlist[[9]],plotlist[[14]],plotlist[[3]],
               plotlist[[10]],plotlist[[15]],plotlist[[4]],
               plotlist[[11]],plotlist[[18]],plotlist[[5]],
               plotlist[[17]],plotlist[[19]],plotlist[[6]])


testlist[[20]]<-ggplot()+
  ylim(0,1)+
  xlim(0,1)+
  theme(axis.line=element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks=element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        legend.position = "none")+
  gg_circle(r=.025, xc=0.93, yc=0.60, color="blue", fill="blue")+
  gg_circle(r=.025, xc=0.93, yc=0.50, color="green", fill="green")+
  gg_circle(r=.025, xc=0.93, yc=0.38, color="brown", fill="brown")+
  gg_circle(r=.025, xc=0.93, yc=0.28, color="black", fill="black")+
  annotate(geom="text", x=.3, y=.5,fontface=2,
           label="Enrichment:\\n   Low\\n                Intermediate\\n    High\\nRecovery",color="black",cex=3)+
  theme(legend.position = 'none')+
  theme(legend.title = element_blank())+
  theme(panel.background = element_blank())+
  theme(axis.line = element_line())

multiplot(plotlist = testlist,cols=6)


########################################################################################################################
############################################## FIGURE S1 ################################################
########################################################################################################################

#read in the 2015 pilot experiment DO data
pilot_2015<-read.csv("pilotdata2015.csv",header=TRUE,row.names = 1)

#subset the data for the relevant replicates
BSA01<-pilot_2015[which(rownames(pilot_2015)%in% c(42,144,39)),]
BSA01_mean<-colMeans(BSA01)
BSA01_sd<-apply(BSA01,2,sd)

BSA05<-pilot_2015[which(rownames(pilot_2015)%in% c(28,72,66)),]
BSA05_mean<-colMeans(BSA05)
BSA05_sd<-apply(BSA05,2,sd)

BSA10<-pilot_2015[which(rownames(pilot_2015)%in% c(13,20, 117)),]
BSA10_mean<-colMeans(BSA10)
BSA10_sd<-apply(BSA10,2,sd)

wasp<-pilot_2015[which(rownames(pilot_2015)%in% c(115,103,14)),]
wasp_mean<-colMeans(wasp)
wasp_sd<-apply(wasp,2,sd)

#create a vector to represent sampling times (days)
pilotdays<-seq(0,10,.5)


par(mar=c(5.2,4.2,2.1,2.1))
plot(y=wasp_mean,x=pilotdays,type="n",col="black",xlab="Time (days)", ylab = "Dissolved Oxygen (%)",ylim=c(0,25))
lines(predict(loess(wasp_mean~pilotdays))~pilotdays, col = "black", lwd = 2)
lines(predict(loess(BSA01_mean~pilotdays))~pilotdays, col = "blue", lwd = 2)
lines(predict(loess(BSA05_mean~pilotdays))~pilotdays, col = "green", lwd = 2)
lines(predict(loess(BSA10_mean~pilotdays))~pilotdays, col = "brown", lwd = 2)

legend("topright",c("Wasp (.5 mg/ml/day)","BSA cocktail (0.1 mg/ml/day)",
                    "BSA cocktail (0.5 mg/ml/day)", "BSA cocktail (1.0 mg/ml/day)"),
       col=c("black","blue","green","brown"),bty="n",lty = 1,lwd=2)

