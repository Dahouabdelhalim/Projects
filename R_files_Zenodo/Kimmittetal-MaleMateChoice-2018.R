setwd("/Users/akimmitt/Documents/MaleChoice/MaleChoice-Paper/Data")

d<-read.csv("Kimmittetal_DataFile_MaleChoice2016-ESM.csv",header=TRUE)

######## LOAD PACKAGES ########

require(car)
require(MASS)
require(lme4)
require(emmeans) #replacement for lsmeans
require(ggplot2)


# add 1 to get rid of zeroes
d$SRS <-d$SRS+1 
d$TS <-d$TS+1
d$PT <-d$PT+1
d$WIN5 <-d$WIN5+1
d$LRS <-d$LRS+1



#Box-Cox Test of Transformation
summary(powerTransform(lm(SRS ~ LureSubsp, data=d))) 
# advocates for the log transform 
summary(powerTransform(lm(TS ~ LureSubsp, data=d))) 
# advocates for the log transform (estimated lambda is close to 0)
summary(powerTransform(lm(PT ~ LureSubsp, data=d))) 
# estimated lambda is close to 0.5, sqrt transformation is better
summary(powerTransform(lm(WIN5 ~ LureSubsp, data=d))) 
# estimated lambda is between 0 and 0.5, so log or sqrt possible; based on histogram, sqrt will be best transformation
summary(powerTransform(lm(LRS ~ LureSubsp, data=d))) 
# estimated lambda is close to 0, so log transform


###### Models ######

########################
######### SRS ##########
########################

srs.m1<-glmer(SRS~(1|LureID)+ LureSubsp, data=d, family="gaussian"(link="log"))

Anova(srs.m1)


emm.srs<-emmeans(srs.m1,"LureSubsp")
emm.srs



### Boxplot Code-SRS ### 

srs.box1<-ggplot(d, aes(x=LureSubsp,y=SRS))+geom_boxplot()+geom_point()
srs.box1
srs.box2<-srs.box1+ theme(axis.title.x = element_text(vjust=-.3, face="bold", size=25),axis.text.x  = element_text(vjust=0.5, size=25)) +
  theme(axis.title.y = element_text(vjust=1.1, face="bold", size=25),axis.text.y  = element_text(vjust=.5, size=25))
srs.box2
srs.box3<-srs.box2+ylab("Short-range song (seconds)") + xlab("Female Migratory Strategy") + scale_color_manual(values = cols,cols <- c("blue3","red1"))
srs.box3

########################
########## TS ##########
########################

ts.m1<-glmer(TS~(1|LureID)+ LureSubsp, data=d,family="gaussian"(link="log"))

Anova(ts.m1)



emm.ts<-emmeans(ts.m1, ~LureSubsp)
emm.ts



### Boxplot Code- TS ### 

ts.box1<-ggplot(d, aes(x=LureSubsp,y=TS))+geom_boxplot()+geom_point()
ts.box1
ts.box2<-ts.box1+ theme(axis.title.x = element_text(vjust=-.3, face="bold", size=25),axis.text.x  = element_text(vjust=0.5, size=25)) +
  theme(axis.title.y = element_text(vjust=1.1, face="bold", size=25),axis.text.y  = element_text(vjust=.5, size=25))
ts.box2
ts.box3<-ts.box2+ylab("Tail Spread (seconds)") + xlab("Female Migratory Strategy") + scale_color_manual(values = cols,cols <- c("blue3","red1"))
ts.box3


########################
########## PT ##########
########################

pt.m1<-glmer(d$PT~(1|LureID)+ LureSubsp, data=d,family="gaussian"(link="sqrt"))

Anova(pt.m1)

emm.pt<- emmeans(pt.m1,~LureSubsp)
emm.pt


### Boxplot Code- PT ### 

pt.box1<-ggplot(d, aes(x=LureSubsp,y=PT))+geom_boxplot()+geom_point()
pt.box1
pt.box2<-pt.box1+ theme(axis.title.x = element_text(vjust=-.3, face="bold", size=25),axis.text.x  = element_text(vjust=0.5, size=25)) +
  theme(axis.title.y = element_text(vjust=1.1, face="bold", size=25),axis.text.y  = element_text(vjust=.5, size=25))
pt.box2
pt.box3<-pt.box2+ylab("Ptiloerection (seconds)") + xlab("Female Migratory Strategy") + scale_color_manual(values = cols,cols <- c("blue3","red1"))
pt.box3

########################
######## WIN5 ##########
########################

win5.m1<-glmer(WIN5~(1|LureID)+ LureSubsp, data=d, family="gaussian"(link="sqrt"))

Anova(win5.m1)

emm.win5<-emmeans(win5.m1, ~LureSubsp)
emm.win5


### Boxplot Code-Win5 ### 

win.box1<-ggplot(d, aes(x=LureSubsp,y=WIN5))+geom_boxplot()+geom_point()
win.box1
win.box2<-win.box1+ theme(axis.title.x = element_text(vjust=-.3, face="bold", size=25),axis.text.x  = element_text(vjust=0.5, size=25)) +
  theme(axis.title.y = element_text(vjust=1.1, face="bold", size=25),axis.text.y  = element_text(vjust=.5, size=25))
win.box2
win.box3<-win.box2+ylab("Time Spent within 5 meters (seconds)") + xlab("Female Migratory Strategy") + scale_color_manual(values = cols,cols <- c("blue3","red1"))
win.box3

########################
########## LRS ##########
########################

lrs.m1<-glmer(LRS~(1|LureID)+ LureSubsp, data=d,family="gaussian"(link="log"))

Anova(lrs.m1)

emm.lrs<-emmeans(lrs.m1, ~LureSubsp)
emm.lrs


### Boxplot Code- LRS ### 

lrs.box1<-ggplot(d, aes(x=LureSubsp,y=LRS))+geom_boxplot()+geom_point()
lrs.box1
lrs.box2<-lrs.box1+ theme(axis.title.x = element_text(vjust=-.3, face="bold", size=25),axis.text.x  = element_text(vjust=0.5, size=25)) +
  theme(axis.title.y = element_text(vjust=1.1, face="bold", size=25),axis.text.y  = element_text(vjust=.5, size=25))
lrs.box2
lrs.box3<-lrs.box2+ylab("Long-range song (count)") + xlab("Female Migratory Strategy") + scale_color_manual(values = cols,cols <- c("blue3","red1"))
lrs.box3
