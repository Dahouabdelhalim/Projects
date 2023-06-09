

# This code allows to test the predictive power of the climate models on bins of 10º width 
# (from 180º West longitude to 180º East longitude) and represent the predicted geographical 
# distribution of the six trophic structures.


#----- Libraries Loading
library(maps)
library(ggplot2)
library(evtree)
library(irr)
library(randomForest)

rm(list=ls(all=TRUE)) # clears workspace
set.seed(1)  # Set random seed

# ----------- Data loading and preparation

# Open Data_III and copy all


read.excel <- function(header=TRUE,...) {
  read.table("clipboard",sep="\\t",header=header,...)
}

Data=read.excel()

Data$TS<-as.factor(Data$TS)

#------- Validación cruzada por bloques

Preds <- data.frame()

for (i in seq(from = -180, to = 170, by = 10)){
  
  traindata <- subset(Data, lon < i|lon>= i+10)
  testdata <- subset(Data,(lon>= i & lon< i+10))
  
  
  RFfit<- randomForest(TS~. -lon -lat,data=traindata,mtry=5,importance=F,ntree=300) 
  print(i)
  
  prd<-predict(RFfit, testdata, type="response")
  testdata$prd<-prd
  Preds<-rbind(Preds,testdata)
}

# windows(); plot(RFfit)

kappatest<-kappa2(Preds[,c(23,24)], "equal")
round(kappatest$value, 2)

tab<-table(Preds$prd,Preds$TS)
tab
round(100*(prop.table(tab,2)),0)

a<-table(Preds$prd==Preds$TS)
paste("CC:",round((a[2]*100)/(a[1]+a[2]),0),"%")


# to represent the distribution of TS predicted

NCDS<-Preds
world<-map_data('world')

windows();ggplot(legend=FALSE) + 
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks = element_blank()) +
  xlab("") + ylab("")+
  geom_point(data=NCDS,aes(x=lon,y=lat,color=prd),size=1.25) + # PTP o FTP
  scale_color_manual(values=c("grey45","navy","skyblue","gold","green3","darkgreen")) + 
  geom_path( data=world, aes(x=long, y=lat,group=group)) +
  labs(title = "Trophic Patterns - Full Data")+
  guides(colour = guide_legend(override.aes = list(size=5)))


