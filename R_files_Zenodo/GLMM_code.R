#Linear Mixed Effects Models
#Jack Weiss gives the best examples and lecture of this material here http://www.unc.edu/courses/2010fall/ecol/563/001/docs/lectures/lecture27.htm

#Set your wd
setwd("/Users/Cody/Desktop/Crossbill_Projects/SH_Call_Project")

#Install nlme and lmer
#install.packages("nlme", dependencies=T)
#install.packages("lme4", dependencies=T)
library(nlme)
data<- read.csv("playback_experiment_data.csv", header=T)
data

#Make sure nlme isn't operating concurrently with lme4
detach(package:nlme)
library(lme4)



call_type<-as.factor(data$Call_Type)


#GLMM with SVL=responsevar, Treatment=FE, Pool=RE nested within treatment
model1<- glmer(Landed~call_type+similarity_to_type2+call_type:similarity_to_type2+(1|Site),data=data, family=binomial(link="logit"))
null<-glmer(Landed~call_type+similarity_to_type2+(1|Site),data=data,family=binomial(link="logit"))
sapply(list(null, model1), AIC)
anova(null,model1)
summary(model1)
summary(null)
