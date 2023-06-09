#code developed by B. Weinstein based on Fridley 2009 for A. Luxbacher's PhD research and edited by M. Lyons for current and future scenarios
Answer<- 0
elev<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_dem.asc", header=FALSE,sep="",na.strings="-9999",dec=".") #elevation ascii
totrad<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_tot.asc", header=FALSE,sep="",na.strings="-9999",dec=".") #annual radiative heating
rad015<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_015.asc", header=FALSE,sep="",na.strings="-9999",dec=".")#radiative heating for January 15th
rad046<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_046.asc", header=FALSE,sep="",na.strings="-9999",dec=".")#radiative heating for Feb 15th
rad074<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_074.asc", header=FALSE,sep="",na.strings="-9999",dec=".")
rad105<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_105.asc", header=FALSE,sep="",na.strings="-9999",dec=".")
rad135<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_135.asc", header=FALSE,sep="",na.strings="-9999",dec=".")
rad166<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_166.asc", header=FALSE,sep="",na.strings="-9999",dec=".")
rad196<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_196.asc", header=FALSE,sep="",na.strings="-9999",dec=".")
rad227<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_227.asc", header=FALSE,sep="",na.strings="-9999",dec=".")
rad258<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_258.asc", header=FALSE,sep="",na.strings="-9999",dec=".")
rad288<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_288.asc", header=FALSE,sep="",na.strings="-9999",dec=".")
rad319<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_319.asc", header=FALSE,sep="",na.strings="-9999",dec=".")
rad349<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_349.asc", header=FALSE,sep="",na.strings="-9999",dec=".")
tci<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_tci.asc", header=FALSE,sep="",na.strings="-9999",dec=".")#topographic convergence index
strdst<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/Asc/jordUTM/jord_str.asc", header=FALSE,sep="",na.strings="-9999",dec=".") #stream distance
PTempmin<-matrix(nrow=608440,ncol=12)
PTempmax<-matrix(nrow=608440,ncol=12)
#force the matrices into columns, so they can be put together in a data frame. This process reads down the 1st column, down the 2nd column...ending at the bottom of final column)
elev1<-(cbind(unmatrix(as.matrix(elev))))
colnames(elev1)<- c("elev")
strdst1<-(cbind(unmatrix(as.matrix(strdst))))
colnames(strdst1)<- c("strdst")
tci1<-(cbind(unmatrix(as.matrix(tci))))
colnames(tci)<- c("tci")
totrad1<-(cbind(unmatrix(as.matrix(totrad))))
rad0151<-(cbind(unmatrix(as.matrix(rad015))))
rad0461<-(cbind(unmatrix(as.matrix(rad046))))
rad0741<-(cbind(unmatrix(as.matrix(rad074))))
rad1051<-(cbind(unmatrix(as.matrix(rad105))))
rad1351<-(cbind(unmatrix(as.matrix(rad135))))
rad1661<-(cbind(unmatrix(as.matrix(rad166))))
rad1961<-(cbind(unmatrix(as.matrix(rad196))))
rad2271<-(cbind(unmatrix(as.matrix(rad227))))
rad2581<-(cbind(unmatrix(as.matrix(rad258))))
rad2881<-(cbind(unmatrix(as.matrix(rad288))))
rad3191<-(cbind(unmatrix(as.matrix(rad319))))
rad3491<-(cbind(unmatrix(as.matrix(rad349))))
rf <- data.frame(elev1,tci1,strdst1,totrad1,rad0151,rad0461,rad0741,rad1051,rad1351,rad1661,rad1961,rad2271,rad2581,rad2881,rad3191,rad3491)
names(rf) <- c("ELEV","TCI","LOG.STRDST","TOTRAD","January","February","March","April","May","June","July","August","September","October","November","December")
mcoef <- read.csv("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/InputValues/mincoef_out.csv",row.names=1) #coefficients based on orignal Fridley 2009 paper
mxcoef <- read.csv("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/InputValues/maxcoef_out.csv",row.names=1)#coefficients based on orignal Fridley 2009 paper
dayP<-read.table("C:/Users/cfanskozaklab/Desktop/Marta/Niche Modeling/4Fridley/InputValues/mapping_predictors_jord_acc_45_55.csv",header=TRUE, sep=",") #lapse rates specific to the climate scenario calculated by rate of temperature change with elevation
rd<-dayP$Rad.day
rd<-format(rd,width=3); rd<-sub(" ","0",rd);rd<-sub(" ","0",rd)
dayP$Rad.day<-rd
dayP$max.lps.int<-dayP$max.lps.int-Answer
dayP$lapse.intercept<-dayP$lapse.intercept-Answer
for(i in 1:12) {
  
  #select a date (format eg "9/4/2006")
  date <- dayP$Date[i]
  print(date)
  scenario = dayP$scenario[i]
  
  #take day-based parameter values from mapping_predictors file
  #JDATE<-dayP$Julian[dayP$Date==date]
  JDATE<-dayP$Julian[i]
  #rd<-dayP$Rad.day[dayP$Date==date]
  rd<-dayP$Rad.day[i]
  rd<-format(rd,width=3); rd<-sub(" ","0",rd);rd<-sub(" ","0",rd)#this corrects for blank spaces 
  max.i<-dayP$max.lps.int[i]
  max.s<-dayP$max.lps.slope[i]
  min.i<-dayP$lapse.intercept[i]
  min.s<-dayP$lapse.slope[i]
  
  #create minSYN and maxSYN values
  rf$minSYN <- min.i + rf$ELEV*min.s
  rf$maxSYN <- max.i + rf$ELEV*max.s
  
  #Calculation of predicted temps for each pixel
  
  #MINIMUM temp
  rf$minT <- mcoef[1,1] + mcoef[2,1]*rf$minSYN + mcoef[3,1]*rf[,4+i] + mcoef[7,1]*rf$LOG.STRDST + mcoef[6,1]*rf$ELEV + mcoef[9,1]*rf$minSYN*rf[,4+i] + mcoef[4,1]*cos(.0172*JDATE) + mcoef[5,1]*sin(.0172*JDATE) + mcoef[10,1]*rf$ELEV*cos(.0172*JDATE) + mcoef[11,1]*rf$ELEV*sin(.0172*JDATE) + mcoef[12,1]*rf$LOG.STRDST*cos(.0172*JDATE) + mcoef[13,1]*rf$LOG.STRDST*sin(.0172*JDATE) + mcoef[8,1]*log(rf$TCI) + mcoef[16,1]*rf$ELEV*rf$minSYN + mcoef[17,1]*rf$LOG.STRDST*rf$minSYN + mcoef[19,1]*rf[,4+i]*rf$ELEV + mcoef[20,1]*rf[,4+i]*rf$LOG.STRDST + mcoef[18,1]*rf$minSYN*log(rf$TCI) + mcoef[15,1]*sin(0.0172*JDATE)*log(rf$TCI) + mcoef[14,1]*cos(0.0172*JDATE)*log(rf$TCI)
  
  #MAX temp
  rf$maxT <- mxcoef[1,1] + mxcoef[2,1]*rf$maxSYN + mxcoef[3,1]*rf[,4+i] + mxcoef[8,1]*rf$LOG.STRDST + mxcoef[7,1]*rf$ELEV + mxcoef[10,1]*rf$maxSYN*rf[,4+i] + mxcoef[4,1]*cos(.0172*JDATE) + mxcoef[5,1]*sin(.0172*JDATE) + mxcoef[13,1]*rf$ELEV*cos(.0172*JDATE) + mxcoef[14,1]*rf$ELEV*sin(.0172*JDATE) + mxcoef[15,1]*rf$LOG.STRDST*cos(.0172*JDATE) + mxcoef[16,1]*rf$LOG.STRDST*sin(.0172*JDATE) + mxcoef[9,1]*log(rf$TCI) + mxcoef[19,1]*rf$ELEV*rf$maxSYN + mxcoef[21,1]*rf$LOG.STRDST*rf$maxSYN + mxcoef[22,1]*rf[,4+i]*rf$ELEV + mxcoef[18,1]*sin(0.0172*JDATE)*log(rf$TCI) + mxcoef[17,1]*cos(0.0172*JDATE)*log(rf$TCI) + mxcoef[12,1]*sin(0.0172*JDATE)*rf$TOTRAD + mxcoef[11,1]*cos(0.0172*JDATE)*rf$TOTRAD + mxcoef[6,1]*rf$TOTRAD + mxcoef[20,1]*rf$maxSYN*rf$TOTRAD + mxcoef[23,1]*rf[,4+i]*log(rf$TCI)
  
  #inm temp (derived)
  rf$meanT <- (rf$minT+rf$maxT)/2
  PTempmin[,i]<-cbind(rf$minT)
  PTempmax[,i]<-cbind(rf$maxT)
}
colnames(PTempmin) <- c("January","February","March","April","May","June","July","August","September","October","November","December")
colnames(PTempmax) <- c("January","February","March","April","May","June","July","August","September","October","November","December")
PTempmax[is.na(PTempmax)]<- -9999
PTempmin[is.na(PTempmin)]<- -9999
JanMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,1])
FebMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,2])
MarchMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,3])
AprilMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,4])
MayMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,5])
JuneMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,6])
JulyMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,7])
AugMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,8])
SeptMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,9])
OctMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,10])
NovMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,11])
DecMin<-matrix(nrow=574,ncol=1060,data=PTempmin[,12])
Janmax<-matrix(nrow=574,ncol=1060,data=PTempmax[,1])
Febmax<-matrix(nrow=574,ncol=1060,data=PTempmax[,2])
Marchmax<-matrix(nrow=574,ncol=1060,data=PTempmax[,3])
Aprilmax<-matrix(nrow=574,ncol=1060,data=PTempmax[,4])
Maymax<-matrix(nrow=574,ncol=1060,data=PTempmax[,5])
Junemax<-matrix(nrow=574,ncol=1060,data=PTempmax[,6])
Julymax<-matrix(nrow=574,ncol=1060,data=PTempmax[,7])
Augmax<-matrix(nrow=574,ncol=1060,data=PTempmax[,8])
Sepmax<-matrix(nrow=574,ncol=1060,data=PTempmax[,9])
Octmax<-matrix(nrow=574,ncol=1060,data=PTempmax[,10])
Novmax<-matrix(nrow=574,ncol=1060,data=PTempmax[,11])
Decmax<-matrix(nrow=574,ncol=1060,data=PTempmax[,12])
write.table(JanMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/015Min.asc",row.names=F,col.names=F) #output minimum average temperature for 2055 under ACCESS GCM for RCP4.5
write.table(FebMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/046Min.asc",row.names=F,col.names=F)
write.table(MarchMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/074Min.asc",row.names=F,col.names=F)
write.table(AprilMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/105Min.asc",row.names=F,col.names=F)
write.table(MayMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/135Min.asc",row.names=F,col.names=F)
write.table(JuneMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/166Min.asc",row.names=F,col.names=F)
write.table(JulyMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/196Min.asc",row.names=F,col.names=F)
write.table(AugMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/227Min.asc",row.names=F,col.names=F)
write.table(SeptMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/258Min.asc",row.names=F,col.names=F)
write.table(OctMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/288Min.asc",row.names=F,col.names=F)
write.table(NovMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/319Min.asc",row.names=F,col.names=F)
write.table(DecMin,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/349Min.asc",row.names=F,col.names=F)
write.table(Janmax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/015Max.asc",row.names=F,col.names=F)
write.table(Febmax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/046Max.asc",row.names=F,col.names=F)
write.table(Marchmax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/074Max.asc",row.names=F,col.names=F)
write.table(Aprilmax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/105Max.asc",row.names=F,col.names=F)
write.table(Maymax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/135Max.asc",row.names=F,col.names=F)
write.table(Junemax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/166Max.asc",row.names=F,col.names=F)
write.table(Julymax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/196Max.asc",row.names=F,col.names=F)
write.table(Augmax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/227Max.asc",row.names=F,col.names=F)
write.table(Sepmax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/258Max.asc",row.names=F,col.names=F)
write.table(Octmax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/288Max.asc",row.names=F,col.names=F)
write.table(Novmax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/319Max.asc",row.names=F,col.names=F)
write.table(Decmax,"K:/NicheModeling/FridLayers/outputs/jordani/ClimateNW/acc_45_2055/349Max.asc",row.names=F,col.names=F)