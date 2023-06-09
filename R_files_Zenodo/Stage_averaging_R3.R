library(ggplot2)
require(dplyr)

## Exclude Moto outliers for statistical tests.

All<-read.csv("Moto_re_stage_3_R3.csv",header = T)

# Calculate average Moto and temperature for each life stages 
# (larva, early juvenile and later juvenile). 
# Data with missing data within each life stage were excluded

# Because Age.range "A" corresponds 30 days interval and "B-G" to 15 days 
# intervals for JP sardine, Age.range "A" of JP sardine was duplicated 
# before averaging for weighting.

JS1=subset(All, All$Age.range=="A" & All$Region.ID=="JP")
All=rbind(All,JS1)
JS1=subset(All, All$Stage=="S1" & All$Region.ID=="JP")
n = JS1 %>% count(Fish.ID)
n2 = n$Fish.ID[n$n==3]
JS1=filter(JS1, JS1$Fish.ID %in% n2)
JS1=aggregate(JS1[,8:11], list(JS1$Fish.ID), mean)
JS1$Stage='S1'
JS1$Region.ID='JP'


JS2=subset(All, All$Stage=="S2" & All$Region.ID=="JP")
n = JS2 %>% count(Fish.ID)
n2 = n$Fish.ID[n$n==2]
JS2=filter(JS2, JS2$Fish.ID %in% n2)
JS2=aggregate(JS2[,8:11], list(JS2$Fish.ID), mean)
JS2$Stage='S2'
JS2$Region.ID='JP'


JS3=subset(All, All$Stage=="S3" & All$Region.ID=="JP")
n = JS3 %>% count(Fish.ID)
n2 = n$Fish.ID[n$n==2]
JS3=filter(JS3, JS3$Fish.ID %in% n2)
JS3=aggregate(JS3[,8:11], list(JS3$Fish.ID), mean)
JS3$Stage='S3'
JS3$Region.ID='JP'


CS1=subset(All, All$Stage=="S1" & All$Region.ID=="CA")
n = CS1 %>% count(Fish.ID)
n2 = n$Fish.ID[n$n==2]
CS1=filter(CS1, CS1$Fish.ID %in% n2)
CS1=aggregate(CS1[,8:11], list(CS1$Fish.ID), mean)
CS1$Stage='S1'
CS1$Region.ID='CA'


CS2=subset(All, All$Stage=="S2" & All$Region.ID=="CA")
CS2=aggregate(CS2[,8:11], list(CS2$Fish.ID), mean)
CS2$Stage='S2'
CS2$Region.ID='CA'


CS3=subset(All, All$Stage=="S3" & All$Region.ID=="CA")
CS3=aggregate(CS3[,8:11], list(CS3$Fish.ID), mean)
CS3$Stage='S3'
CS3$Region.ID='CA'


All2=rbind(JS1,JS2,JS3,CS1,CS2,CS3) # Restore data
All2$TG=round(All2$Estimated.Temperature)

# write.csv(All2,file='Moto_T_bin_re_R3.csv') 

# Outliers were detected and removed for Moto analyses.

boxplot(JS1$Moto.mean, plot=FALSE)$out
outliers <- boxplot(JS1$Moto.mean, plot=FALSE)$out
JS1<- JS1[-which(JS1$Moto.mean %in% outliers),]

boxplot(JS2$Moto.mean, plot=FALSE)$out
outliers <- boxplot(JS2$Moto.mean, plot=FALSE)$out
JS2<- JS2[-which(JS2$Moto.mean %in% outliers),]

boxplot(JS3$Moto.mean, plot=FALSE)$out
outliers <- boxplot(JS3$Moto.mean, plot=FALSE)$out
JS3<- JS3[-which(JS3$Moto.mean %in% outliers),]

boxplot(CS1$Moto.mean, plot=FALSE)$out
#outliers <- boxplot(CS1$Moto.mean, plot=FALSE)$out
#CS1<- CS1[-which(CS1$Moto.mean %in% outliers),]


boxplot(CS2$Moto.mean, plot=FALSE)$out
#outliers <- boxplot(CS2$Moto.mean, plot=FALSE)$out
#CS2<- CS2[-which(CS2$Moto.mean %in% outliers),]

boxplot(CS3$Moto.mean, plot=FALSE)$out
outliers <- boxplot(CS3$Moto.mean, plot=FALSE)$out
CS3<- CS3[-which(CS3$Moto.mean %in% outliers),]


All3=rbind(JS1,JS2,JS3,CS1,CS2,CS3) # Restore data
All3$TG=round(All3$Estimated.Temperature)

write.csv(All3,file='Moto_T_bin_re_outed_R3.csv')