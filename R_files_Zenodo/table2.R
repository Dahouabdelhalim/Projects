#
# computes and prints models for table 2 of Schomaker, J., Walper, D., Wittmann, B.C., & Einhäuser, W. (2017). Attention in natural scenes: Affective-motivational factors guide gaze independently of visual salience. Vision Research, 133, 161-175..
#

# load library
library(lme4)

#clear work space
rm(list = ls()) 

#load data
load("dataForExperiment2.Rdata")

# center and scale the predictors 
data$arousal = scale(data$arousal,center=TRUE,scale=TRUE)
data$valence = scale(data$valence,center=TRUE,scale=TRUE)
data$valence2 = scale(data$valence2,center=TRUE,scale=TRUE)
data$motivation = scale(data$motivation,center=TRUE,scale=TRUE)
data$motivation2 = scale(data$motivation2,center=TRUE,scale=TRUE)

data$objSize = scale(data$objSize,center=TRUE,scale=TRUE)
data$objEcc = scale(data$objEcc,center=TRUE,scale=TRUE)
data$peakSal = scale(data$peakSal,center=TRUE,scale=TRUE)


model1<-glmer((isFixated) ~ arousal + valence + valence2 + motivation + motivation2 + objSize + objEcc + peakSal + (0 + arousal | subID) + (0 + valence | subID) + (0 + valence2 | subID) +  (0 + motivation | subID) +  (0 + motivation2 | subID)  +  (0 + peakSal | subID) + (0 + objSize | subID) + (0 + objEcc | subID) + (0 + arousal | imgID) + (0 + valence | imgID) +  (0 + valence2 | imgID) + (0 + motivation | imgID) +  (0 + motivation2 | imgID) +  (0 + peakSal | imgID) + (0 + objSize | imgID) + (0 + objEcc | imgID)+ (1 | subID) + (1  | imgID), verbose=2, control=glmerControl(optimizer="bobyqa"), data=data, family=binomial)

print(summary(model1))


model2<-lmer((totalDwell) ~ arousal + valence  + valence2 + motivation + motivation2 + objSize + objEcc + peakSal + (0 + arousal | subID) + (0 + valence | subID) + (0 + valence2 | subID) +  (0 + motivation | subID) +  (0 + motivation2 | subID)  +  (0 + peakSal | subID) + (0 + objSize | subID) + (0 + objEcc | subID) + (0 + arousal | imgID) + (0 + valence | imgID) +(0 + valence2 | imgID) +  (0 + motivation | imgID) +  (0 + motivation2 | imgID) +  (0 + peakSal | imgID) + (0 + objSize | imgID) + (0 + objEcc | imgID)+ (1 | subID) + (1  | imgID),verbose=2, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e7)), data=data)


print(summary(model2))



ix<-which(is.finite(data$tFirstFix))
data<-data[ix,]    

model3<-lmer((tFirstFix) ~ arousal + valence + valence2 + motivation + motivation2 + objSize + objEcc + peakSal + (0 + arousal | subID) + (0 + valence | subID) +  (0 + valence2 | subID) +  (0 + motivation | subID) +  (0 + motivation2 | subID)  +  (0 + peakSal | subID) + (0 + objSize | subID) + (0 + objEcc | subID) + (0 + arousal | imgID) + (0 + valence | imgID) + (0 + valence2 | imgID) +  (0 + motivation | imgID) +  (0 + motivation2 | imgID) +  (0 + peakSal | imgID) + (0 + objSize | imgID) + (0 + objEcc | imgID)+ (1 | subID) + (1  | imgID),verbose=2, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e7)), data=data)


print(summary(model3))
