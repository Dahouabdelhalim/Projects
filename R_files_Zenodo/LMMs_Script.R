######## 1) Comparison of the 1) total distance traveled and 2) track duration and 3) the speed between outgoing and incoming tracks. 

#The script uses the following libraries, make sure you install them beforehand.

library(nlme)
library(lme4)

#1.1) Load the file

x <- read.table(file.choose(), header=T, sep=",")

#The file should be composed of a column with :
#1) the id of the penguin (ID_Penguin) (as.factor)
#2) the total distance traveled in metres (Total_distance_traveled_m) (as.numeric)
#3) the track duration in minutes (Track_duration_min) (as.numeric)
#4) the speed (Speed) computed by dividing the total distance traveled by the track duration (in seconds) in order to have the speed in m/s. This column should be as.numeric
#5) the track segment (OutgoingvsIncoming) (as.factor)
 

#1.2) Run the LMM models with the id of the penguin as a random effect

#1.2.1 ) distance comparison 

model1.1 <- lme(Total_distance_traveled_m ~ OutgoingvsIncoming, random=~ 1|ID_Penguin, data=x, method="ML")
summary(model1.2)  

model1.2 <- lme(Total_distance_traveled_m ~ 1, random=~ 1|ID_Penguin, data=x1,method="ML")
summary(model1.3) 

anova(model1.1,model1.2) #comparison between both models to determine if the distance traveled between both segments differed. 

#1.2.2) duration comparison

model2.1 <- lme(Track_duration_min ~ OutgoingvsIncoming,random=~ 1|ID_Penguin, data=x,method="ML")
summary(model2.2)

model2.2 <- lme(Track_duration_min ~ 1,random=~ 1|ID_Penguin, data=x, method="ML" )
summary(model2.3)

anova(model2.1, model2.2) #comparison between both models to determine if track duration of both segments differed. 

#1.2.3) speed comparison

model3.1 <- lme(Speed ~ OutgoingvsIncoming,random=~ 1|ID_Penguin, data=x,method="ML")
summary(model2.2)

model3.2 <- lme(Speed ~ 1,random=~ 1|ID_Penguin, data=x, method="ML" )
summary(model2.3)

anova(model3.1, model3.2) #comparison between both models to determine if the speed of both segments differed. 

######## 2) Only for Y tracks comparison of the speed between the I-segment part and the non I-segment part of the outgoing and incoming tracks

#2.1) Load the file 

y <- read.table(file.choose(), header=T, sep=",")

#The file should be composed of a column with :
#1) the id of the penguin (ID_Penguin) (as.factor)
#2) the track segment (OutgoingvsIncoming) (as.factor)
#3) if the section of the path was in or out of the I-segment (as.factor)
#4) the speed in m/s. This column should be (as.numeric)

#2.2) Run the LMM models with the id of the penguin as a random effect

#2.2.1 ) Comparison of the outgoing and incoming speed while birds are on the I-segment

#Select from the file only the speed data while birds were in the I-segment

idx=which(y$Out_or_In_Isegment="In")
y_Iseg=y[idx,]

model4.1 <- lme(Speed ~ Outgoing_Incoming, random=~ 1|ID_Penguin, data=y_Iseg, method="ML")
summary(model4.1)  

model4.2 <- lme(Speed ~ 1, random=~ 1|ID_Penguin, data=y_Iseg,method="ML")
summary(model4.2) 

anova(model4.1,model4.2) #comparison between both models to determine if the speed of the I-segment during outbound and inbound tracks differed 

#2.2.2) Comparison of the incoming speed while birds were in or off the I-segment

#Select from the file only the speed data of incoming path

idx1=which(y$Outgoing_Incoming="Incoming")
y_Incoming=y[idx1,]

model5.1 <- lme(Speed ~ Out_or_In_Isegment, random=~ 1|ID_Penguin, data=y_Incoming, method="ML")
summary(model5.1)  

model5.2 <- lme(Speed ~ 1, random=~ 1|ID_Penguin, data=y_Incoming,method="ML")
summary(model5.2) 

anova(model5.1,model5.2) #comparison between both models to determine if the speed of the inbound path in or off the I-segment differed.  





 





