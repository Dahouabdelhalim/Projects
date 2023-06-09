
library(weights)
library(lmerTest)

# Load datafile CTI_Longterm.csv
NAERCTI=read.csv(file.choose())

###
### Testing connection between mean CTI and mean temperature of regions 

#mTempSt=NAERCTI$mTemp-mean(NAERCTI$mTemp)

# Testing seasonal interaction (see Fig. 1)
ctislopem1=lmer(mCTI ~ mTempSt*Season +(1|Continent/Region), data = NAERCTI,REML=TRUE,control = lmerControl(optCtrl = list(maxfun = 20000)))
summary(ctislopem1)

# Only breeding season (see Fig. 1)
ctislopeBr=lm(NAERCTI$mCTI[which(NAERCTI$Season=="B")]~NAERCTI$mTempSt[which(NAERCTI$Season=="B")])
summary(ctislopeBr)

# Only winter season (see Fig. 1)
ctislopeWi=lm(NAERCTI$mCTI[which(NAERCTI$Season=="W")]~NAERCTI$mTempSt[which(NAERCTI$Season=="W")])
summary(ctislopeWi)

# Long-term analyses
# Calculating the weighting

NAERCTI$W2=1/(NAERCTI$CTI_se^2)
NAERCTI$Weight=(NAERCTI$W2/max(NAERCTI$W2,na.rm=TRUE))

### Comparison of observed and expected CTI changes during breeding and winter seasons (see Fig. 1C-D)
#breeding

ctitemp1=NAERCTI$dCTI[which(NAERCTI$Season=="B")]
ctitemp1=as.data.frame(ctitemp1)
ctitemp1$type="CTI"
ctitemp1$cont=NAERCTI$Continent[which(NAERCTI$Season=="B")]
ctitemp1$weight=NAERCTI$Weight[which(NAERCTI$Season=="B")]


ctitemp2=NAERCTI$CTITempPred[which(NAERCTI$Season=="B")]
ctitemp2=as.data.frame(ctitemp2)
ctitemp2$type="ExpCTI"
ctitemp2$cont=NAERCTI$Continent[which(NAERCTI$Season=="B")]
ctitemp2$weight=NAERCTI$Weight[which(NAERCTI$Season=="B")]
colnames(ctitemp2)=c("ctitemp1","type","cont","weight")

ctitempB=rbind(ctitemp1,ctitemp2)

obspreB=lm(ctitemp1~type*cont,weight=weight, data=ctitempB)
summary(obspreB)

#winter

ctitemp1=NAERCTI$dCTI[which(NAERCTI$Season=="W")]
ctitemp1=as.data.frame(ctitemp1)
ctitemp1$type="CTI"
ctitemp1$cont=NAERCTI$Continent[which(NAERCTI$Season=="W")]
ctitemp1$weight=NAERCTI$Weight[which(NAERCTI$Season=="W")]


ctitemp2=NAERCTI$CTITempPred[which(NAERCTI$Season=="W")]
ctitemp2=as.data.frame(ctitemp2)
ctitemp2$type="ExpCTI"
ctitemp2$cont=NAERCTI$Continent[which(NAERCTI$Season=="W")]
ctitemp2$weight=NAERCTI$Weight[which(NAERCTI$Season=="W")]
colnames(ctitemp2)=c("ctitemp1","type","cont","weight")

ctitempW=rbind(ctitemp1,ctitemp2)

obspreW=lm(ctitemp1~type*cont,weight=weight, data=ctitempW)

summary(obspreW)

std <- function(x) sd(x)/sqrt(length(x))
mean(ctitempW$ctitemp1[which(ctitempW$type=="CTI" & ctitempW$cont=="N")])
std(ctitempW$ctitemp1[which(ctitempW$type=="CTI" & ctitempW$cont=="N")])
mean(ctitempW$ctitemp1[which(ctitempW$type=="ExpCTI" & ctitempW$cont=="N")])
std(ctitempW$ctitemp1[which(ctitempW$type=="ExpCTI" & ctitempW$cont=="N")])
mean(ctitempW$ctitemp1[which(ctitempW$type=="CTI" & ctitempW$cont=="E")])
std(ctitempW$ctitemp1[which(ctitempW$type=="CTI" & ctitempW$cont=="E")])
mean(ctitempW$ctitemp1[which(ctitempW$type=="ExpCTI" & ctitempW$cont=="E")])
std(ctitempW$ctitemp1[which(ctitempW$type=="ExpCTI" & ctitempW$cont=="E")])

mean(ctitempB$ctitemp1[which(ctitempW$type=="CTI" & ctitempW$cont=="N")])
std(ctitempB$ctitemp1[which(ctitempW$type=="CTI" & ctitempW$cont=="N")])
mean(ctitempB$ctitemp1[which(ctitempW$type=="ExpCTI" & ctitempW$cont=="N")])
std(ctitempB$ctitemp1[which(ctitempW$type=="ExpCTI" & ctitempW$cont=="N")])
mean(ctitempB$ctitemp1[which(ctitempW$type=="CTI" & ctitempW$cont=="E")])
std(ctitempB$ctitemp1[which(ctitempW$type=="CTI" & ctitempW$cont=="E")])
mean(ctitempB$ctitemp1[which(ctitempW$type=="ExpCTI" & ctitempW$cont=="E")])
std(ctitempB$ctitemp1[which(ctitempW$type=="ExpCTI" & ctitempW$cont=="E")])




# Long-term analyses temperature, season, continent

twm00<- blmer(dCTI ~ 1+ Continent+ (1|Region) , weights=Weight,data = NAERCTI,REML=FALSE)
twm0<- blmer(dCTI ~ Temp+Continent+ (1|Region) , weights=Weight,data = NAERCTI,REML=FALSE)
twm1<- blmer(dCTI ~ Temp+Season+Continent+ (1|Region) , weights=Weight,data = NAERCTI,REML=FALSE)
twm2<- blmer(dCTI ~ Temp*Season+Continent+(1|Region) , weights=Weight,data = NAERCTI,REML=FALSE)
twm3<- blmer(dCTI ~ Temp+Season*Continent+ (1|Region) , weights=Weight,data = NAERCTI,REML=FALSE)
twm4<- blmer(dCTI ~ Temp*Season+Season*Continent+ (1|Region) , weights=Weight,data = NAERCTI,REML=FALSE)

library(MuMIn)
library(parameters)
 
outputctiall<-model.sel(twm00,twm0,twm1,twm2,twm3,twm4)
outputctiall

summary(twm1)
parameters::p_value(twm1_final)

# Season specific analyses


twW<- lm(dCTI ~ Temp, weights=Weight,data = NAERCTI[which(NAERCTI$Season=="W"),])
summary(twW)

twB<- lm(dCTI ~ Temp, weights=Weight,data = NAERCTI[which(NAERCTI$Season=="B"),])
summary(twB)

