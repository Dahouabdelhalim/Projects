library(ggplot2)
require(dplyr)
library(quantreg)

### Optimal temperature estimation ###############################
### Model upper and lower limits of Moto distribution and find the
### temperature that maximise the gap between upper and lower limits.

All<-read.csv("Moto_re_stage_3_R3.csv",header = T)

n_th=4 # use temperature bins that include >=4 individuals

# Calculate average Moto and temperature for each life stages 
# (larva, early juvenile and later juvenile). 
# Data with missing data within each life stage were excluded
# Outliers were detected and also removed.

# Because Age.range "A" corresponds 30 days interval and "B-G" to 15 days 
# intervals for JP sardine, Age.range "A" of JP sardine was duplicated 
# before averaging for weighting.
JS1=subset(All, All$Age.range=="A" & All$Region.ID=="JP") 
All=rbind(All,JS1) # Duplicate Age.range A of JP
JS1=subset(All, All$Stage=="S1" & All$Region.ID=="JP") 
n = JS1 %>% count(Fish.ID)
n2 = n$Fish.ID[n$n==3] 
JS1=filter(JS1, JS1$Fish.ID %in% n2) # remove individuals with missing data
JS1=aggregate(JS1[,10:11], list(JS1$Fish.ID), mean)
JS1$Stage='S1'
JS1$Region.ID='JP'
JS1$TG=round(JS1$Estimated.Temperature) # Temperature bins 

boxplot(JS1$Moto.mean, plot=FALSE)$out # Detect Moto outliers
outliers <- boxplot(JS1$Moto.mean, plot=FALSE)$out
JS1<- JS1[-which(JS1$Moto.mean %in% outliers),] # remove outliers
n = JS1 %>% count(TG)
n2 = n$TG[n$n>=n_th] # remove temperature bins that includes < 4 data
JS1=filter(JS1, JS1$TG %in% n2)


JS2=subset(All, All$Stage=="S2" & All$Region.ID=="JP")
n = JS2 %>% count(Fish.ID)
n2 = n$Fish.ID[n$n==2]
JS2=filter(JS2, JS2$Fish.ID %in% n2)
JS2=aggregate(JS2[,10:11], list(JS2$Fish.ID), mean)
JS2$Stage='S2'
JS2$Region.ID='JP'

boxplot(JS2$Moto.mean, plot=FALSE)$out
outliers <- boxplot(JS2$Moto.mean, plot=FALSE)$out
JS2<- JS2[-which(JS2$Moto.mean %in% outliers),]
JS2$TG=round(JS2$Estimated.Temperature)
n = JS2 %>% count(TG)
n2 = n$TG[n$n>=n_th]
JS2=filter(JS2, JS2$TG %in% n2)


JS3=subset(All, All$Stage=="S3" & All$Region.ID=="JP")
n = JS3 %>% count(Fish.ID)
n2 = n$Fish.ID[n$n==2]
JS3=filter(JS3, JS3$Fish.ID %in% n2)
JS3=aggregate(JS3[,10:11], list(JS3$Fish.ID), mean)
JS3$Stage='S3'
JS3$Region.ID='JP'

boxplot(JS3$Moto.mean, plot=FALSE)$out
outliers <- boxplot(JS3$Moto.mean, plot=FALSE)$out
JS3<- JS3[-which(JS3$Moto.mean %in% outliers),]
JS3$TG=round(JS3$Estimated.Temperature)
n = JS3 %>% count(TG)
n2 = n$TG[n$n>=n_th]
JS3=filter(JS3, JS3$TG %in% n2)



CS1=subset(All, All$Stage=="S1" & All$Region.ID=="CA")
n = CS1 %>% count(Fish.ID)
n2 = n$Fish.ID[n$n==2]
CS1=filter(CS1, CS1$Fish.ID %in% n2)
CS1=aggregate(CS1[,10:11], list(CS1$Fish.ID), mean)
CS1$Stage='S1'
CS1$Region.ID='CA'

boxplot(CS1$Moto.mean, plot=FALSE)$out
#outliers <- boxplot(CS1$Moto.mean, plot=FALSE)$out # No outliers detected
#CS1<- CS1[-which(CS1$Moto.mean %in% outliers),]
CS1$TG=round(CS1$Estimated.Temperature)
n = CS1 %>% count(TG)
n2 = n$TG[n$n>=n_th]
CS1=filter(CS1, CS1$TG %in% n2)



CS2=subset(All, All$Stage=="S2" & All$Region.ID=="CA")
CS2=aggregate(CS2[,10:11], list(CS2$Fish.ID), mean)
CS2$Stage='S2'
CS2$Region.ID='CA'

boxplot(CS2$Moto.mean, plot=FALSE)$out
#outliers <- boxplot(CS2$Moto.mean, plot=FALSE)$out # No outliers detected
#CS2<- CS2[-which(CS2$Moto.mean %in% outliers),]
CS2$TG=round(CS2$Estimated.Temperature)
n = CS2 %>% count(TG)
n2 = n$TG[n$n>=n_th]
CS2=filter(CS2, CS2$TG %in% n2)



CS3=subset(All, All$Stage=="S3" & All$Region.ID=="CA")
CS3=aggregate(CS3[,10:11], list(CS3$Fish.ID), mean)
CS3$Stage='S3'
CS3$Region.ID='CA'
boxplot(CS3$Moto.mean, plot=FALSE)$out
outliers <- boxplot(CS3$Moto.mean, plot=FALSE)$out
CS3<- CS3[-which(CS3$Moto.mean %in% outliers),]
CS3$TG=round(CS3$Estimated.Temperature)
n = CS3 %>% count(TG)
n2 = n$TG[n$n>=n_th]
CS3=filter(CS3, CS3$TG %in% n2)

All2=rbind(JS1,JS2,JS3,CS1,CS2,CS3) # Restore filtered data 

# Create boxplot for overall view 

ggplot(All2, aes(x=Estimated.Temperature, y=Moto.mean))+
  geom_point()+
  stat_smooth(method="lm")+
  facet_grid(Stage~Region.ID)+
  theme_bw( )+  
  theme(axis.title = element_text(size = rel(1.5)),axis.text=element_text(size=rel(1.3),face="bold", colour="black"),
        strip.text=element_text(size = rel(1.3)),legend.text =element_text(size = rel(1.0),face="bold"),
        legend.title =element_text(size = rel(1),face="bold"))+
  xlab("Temperature")+
  ylab("Moto")

ggplot(data=All2,aes(x=factor(TG),y= Moto.mean))+
  geom_boxplot(position=position_dodge(1))+
  geom_jitter(size=0.4, alpha=0.9)+
  facet_grid(Stage~Region.ID)+
  theme_bw( )+  
  theme(axis.title = element_text(size = rel(1.5)),axis.text=element_text(size=rel(1.0),face="bold", colour="black"),
        strip.text=element_text(size = rel(1.3)),legend.text =element_text(size = rel(1.0),face="bold"),
        legend.title =element_text(size = rel(1),face="bold"))+
  xlab("Temperature")+
  ylab("Moto")


#model upper and lower limit of data distribution per Region and stage
b<-aggregate(All2$Moto.mean, list(All2$Region.ID,All2$Stage, All2$TG), FUN = function(x) quantile(x, probs = 0.95, na.rm=T))
c<-aggregate(All2$Moto.mean, list(All2$Region.ID, All2$Stage, All2$TG), FUN = function(x) quantile(x, probs = 0.05,na.rm=T))


colnames(b) <- c("Region.ID","Stage", "Temp", "M95")
b$M5<-NA
b$M5<-c$x


J1<-b[(b$Region.ID=="JP"&b$Stage=="S1"),]# JP S1 model
J1M1 <- lm(M95 ~ poly(Temp,2,raw=TRUE), data=J1)
J1M2<-glm(M5~Temp,data=J1, family = gaussian(link = 'log'))


J1x<-seq(14, 22, by=0.01)
J195<- J1M1$coefficients[1]+ J1M1$coefficients[2]*J1x+ J1M1$coefficients[3]*(J1x^2)
J105<-exp(J1M2$coefficients[1]+J1M2$coefficients[2]*J1x)
J1d<-J195-J105

J2<-b[(b$Region.ID=="JP"&b$Stage=="S2"),]# JP S2 model
J2M1 <- lm(M95 ~ poly(Temp,2,raw=TRUE), data=J2)
J2M2<-glm(M5~Temp,data=J2, family = gaussian(link = 'log'))
JS2=JS2[(JS2$Estimated.Temperature>=13.5&JS2$Estimated.Temperature<=22.5),]

J2x<-seq(14, 22, by=0.01)
J295<- J2M1$coefficients[1]+ J2M1$coefficients[2]*J2x+ J2M1$coefficients[3]*(J2x^2)
J205<-exp(J2M2$coefficients[1]+J2M2$coefficients[2]*J2x)
J2d<-J295-J205

J3<-b[(b$Region.ID=="JP"&b$Stage=="S3"),]# JP S3 model
J3M1 <- lm(M95 ~ poly(Temp,2,raw=TRUE), data=J3)
J3M2<-glm(M5~Temp,data=J3, family = gaussian(link = 'log'))
JS3=JS3[(JS3$Estimated.Temperature>=11.5&JS3$Estimated.Temperature<22.5),]


J3x<-seq(12, 22, by=0.01)
J395<- J3M1$coefficients[1]+ J3M1$coefficients[2]*J3x+ J3M1$coefficients[3]*(J3x^2)
J305<-exp(J3M2$coefficients[1]+J3M2$coefficients[2]*J3x)
J3d<-J395-J305


C1<-b[(b$Region.ID=="CA"&b$Stage=="S1"),]# CA S1 model
C1M1 <- lm(M95 ~ poly(Temp,2,raw=TRUE), data=C1)
C1M2<-glm(M5~Temp,data=C1, family = gaussian(link = 'log'))
CS1=CS1[(CS1$Estimated.Temperature>=10.5&CS1$Estimated.Temperature<19.5),]


C1x<-seq(12, 19, by=0.01)
C195<- C1M1$coefficients[1]+ C1M1$coefficients[2]*C1x+ C1M1$coefficients[3]*(C1x^2)
C105<-exp(C1M2$coefficients[1]+C1M2$coefficients[2]*C1x)
C1d<-C195-C105

C2<-b[(b$Region.ID=="CA"&b$Stage=="S2"),]# CA S2 model
C2M1 <- lm(M95 ~ poly(Temp,2,raw=TRUE), data=C2)
C2M2<-glm(M5~Temp,data=C2, family = gaussian(link = 'log'))
CS2=CS2[(CS2$Estimated.Temperature>=11.5&CS2$Estimated.Temperature<20.5),]


C2x<-seq(12, 20, by=0.01)
C295<- C2M1$coefficients[1]+ C2M1$coefficients[2]*C2x+ C2M1$coefficients[3]*(C2x^2)
C205<-exp(C2M2$coefficients[1]+C2M2$coefficients[2]*C2x)
C2d<-C295-C205

C3<-b[(b$Region.ID=="CA"&b$Stage=="S3"),]# CA S3 model
C3M1 <- lm(M95 ~ poly(Temp,2,raw=TRUE), data=C3)
C3M2<-glm(M5~Temp,data=C3, family = gaussian(link = 'log'))
CS3=CS3[(CS3$Estimated.Temperature>=11.5&CS3$Estimated.Temperature<20.5),]


C3x<-seq(13, 18, by=0.01)
C395<- C3M1$coefficients[1]+ C3M1$coefficients[2]*C3x+ C3M1$coefficients[3]*(C3x^2)
C305<-exp(C3M2$coefficients[1]+C3M2$coefficients[2]*C3x)
C3d<-C395-C305


# Output model coefficients
J1M1$coefficients[1:3]
J1M2$coefficients[1:2]
J2M1$coefficients[1:3]
J2M2$coefficients[1:2]
J3M1$coefficients[1:3]
J3M2$coefficients[1:2]
C1M1$coefficients[1:3]
C1M2$coefficients[1:2]
C2M1$coefficients[1:3]
C2M2$coefficients[1:2]
C3M1$coefficients[1:3]
C3M2$coefficients[1:2]


# Make figures 
par(mfrow=c(4,4),mai=c(0.3,0.3,0.3,0.3))
plot(x=J1x,y=J195,ylim=c(0.1,0.5),col="red")
points(x=J1x,y=J105, col="blue")
plot(x=J1x,y=J1d)
abline(v = c(J1x[which.max(J1d)]), col="black", lwd=3, lty=2)

plot(x=C1x,y=C195,ylim=c(0.1,0.5),col="red")
points(x=C1x,y=C105, col="blue")
plot(x=C1x,y=C1d)
abline(v = c(C1x[which.max(C1d)]), col="black", lwd=3, lty=2)

plot(x=J2x,y=J295,ylim=c(0.1,0.5),col="red")
points(x=J2x,y=J205, col="blue")
plot(x=J2x,y=J2d)
abline(v = c(J2x[which.max(J2d)]), col="black", lwd=3, lty=2)


plot(x=C2x,y=C295,ylim=c(0.1,0.5),xlim=c(9,18),col="red")
points(x=C2x,y=C205, col="blue")
plot(x=C2x,y=C2d)
abline(v = c(C2x[which.max(C2d)]), col="black", lwd=3, lty=2)

plot(x=J3x,y=J395,ylim=c(0.1,0.5),col="red")
points(x=J3x,y=J305, col="blue")
plot(x=J3x,y=J3d)
abline(v = c(J3x[which.max(J3d)]), col="black", lwd=3, lty=2)

plot(x=C3x,y=C395,ylim=c(0.1,0.5),xlim=c(9,18),col="red")
points(x=C3x,y=C305, col="blue")
plot(x=C3x,y=C3d)
abline(v = c(C3x[which.max(C3d)]), col="black", lwd=3, lty=2)

#output optimal temperature
J1x[which.max(J1d)]
J2x[which.max(J2d)]
J3x[which.max(J3d)]
C1x[which.max(C1d)]
C2x[which.max(C2d)]
C3x[which.max(C3d)]