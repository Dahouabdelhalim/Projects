#Author: Joeri Reinders
#Purpose: Densityplot of the months in which peakflows (on the Neches, Trinity and Brazos Rivers) occur throughout the year- figure 3 in paper.
#Data: contact reinders.j@northeastern.edu for orginal data files

rm(list=ls())
x11()

#packages
library(dplyr)
library(ggplot2)
library(gridExtra)

rm(list=ls())

directory <- 'D:/Data/USGS_Peakflow/Processed/'
file_name <- 'seasonality_peakflow_Ahouston.csv'

peakflow_houston <- read.csv(paste(directory,file_name, sep=""))

peaks <- data.frame(Months= rep(0,each=9), Rivers=rep(c("Brazos","Trinity","Neches"), c(3,3,3)), Top=rep(c("Highest","Second","Third")))
cell <- c(1,4,7)

for(t in 1:3){

highestvalues <- data.frame(Discharge = peakflow_houston[,2*t])
limiet <- as.numeric(min(top_n(highestvalues, n=3)))

highestvalues <- data.frame(Discharge = peakflow_houston[,2*t], month = peakflow_houston[,(2*t)+1])
highestvalues <- subset(highestvalues, Discharge>limiet-1 )
highestvalues <- highestvalues[order(highestvalues[,1]),]
peaks[cell[t]:(cell[t]+2),1] <- highestvalues$month

}



Brazos <- data.frame(Month=peakflow_houston$Brazos_Month, River=rep("Brazos", each=length(peakflow_houston$Brazos_Month)))
Brazos <- na.omit( Brazos )
Trinity <- data.frame(Month=peakflow_houston$Trinity_Month, River=rep("Trinity", each=length(peakflow_houston$Trinity_Month)))
Trinity <- na.omit( Trinity )
Neches <- data.frame(Month=peakflow_houston$Neches_month, River=rep("Neches", each=length(peakflow_houston$Neches_month)))
Neches <- na.omit( Neches )

months <- data.frame(Months= c(Brazos$Month,Trinity$Month, Neches$Month), Rivers=c(as.character(Brazos$River),as.character(Trinity$River),as.character(Neches$River)))
peaks$y <- 0

for(i in 1:9){
  if(i==1){peaks[i,4] <- 0.25}
  if(i>1){  
    
    check <- peaks[i,1]
    list <- peaks[1:(i-1),1]
    
    if(check %in% list == T){
      if(length(which(check == list))>1){peaks[i,4] <- 0.22}
      else{peaks[i,4] <- 0.235}
      
      }
    else{peaks[i,4] <- 0.25}
  }
}

maanden <- c("J","F","M","A","M","J","J","A","S","O","N","D")

ggplot(data=months, aes(x=Months,colour=Rivers, fill=Rivers))+
  geom_density(alpha=0.3, size=1.5)+ 
  scale_colour_manual(values=c("#333652", "#90ADC6","#FAD02C"))+
  scale_fill_manual(values=c("#333652","#90ADC6","#FAD02C"))+
  geom_point(data=peaks, aes(x=Months, y=y, shape=Top), size=4)+
  scale_shape_manual(values=c(17,19,15))+
  scale_x_continuous(breaks=seq(1,12,1), limits=c(1,12), labels = maanden)+
  ylim(0,0.25)+
  ylab("Density")+
  theme_classic()+
  theme(plot.title= element_text(hjust=0, size=18,vjust =10))+
  theme(axis.title.x=element_text(size=17, hjust=0.5,margin = margin(t = 20)))+
  theme(axis.title.y=element_text(size=17, hjust=0.5,margin = margin(r = 20)))+
  theme(axis.text=element_text(size=14))+
  theme(legend.position = c(.91, .80),)+
  theme(legend.box.background=element_rect(),legend.box.margin = margin(6, 6, 6, 6))+
  theme(legend.text=element_text(size=14), legend.title=element_text(size=17, face="bold"))+
  theme(plot.margin = margin(3, 4, 0.5, 1, "cm"))


#("#122620","#B68D40",  "#F4EBD0") noga
#("#74BDCB","#FFA384",  "#EFE7BC") zalm
#("#2F5061","#E57F84",  "#F4EAE6") SF
#("#0A4158","#FF9636",  "#E4D7D0") orange
#("#333652","#90ADC6","#FAD02C") sneaker


