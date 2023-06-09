#Load Library
library(car)
library(lme4)
library(lmerTest)
library("plotrix")
library(MuMIn)
options(na.action = "na.fail") 
library(polynom)

#Load data file
BSTP<-read.csv("...Brood Schedule TP Combined.csv",header=TRUE)
BSTP<-subset(BSTP,Day>109) #Remove data before April because it was too cold for extranidal activity

#Define calculated variables
BSTP$PropWPup<-(BSTP$Wpupae/BSTP$TotB)
BSTP$NonWBrood<-(BSTP$TotB-BSTP$Wpupae)
BSTP$PropL<-(BSTP$L/BSTP$TotB)
BSTP<-subset(BSTP,Day =!is.na(PropWPup))


#Original Pam v Tlong peak in brood production comparison
#=========================================================
QuadP<-lm(PropWPup~Day+I(Day^2),data=BSTP[BSTP$Sp=="P",],na.action="na.omit")
summary(QuadP)


#Will need to subset the Tlong data to exclude TLN488 with only one brood item dissected in Feb
#plot(PropWPup~Day,data=BSTP[BSTP$Sp=="T",])
BSTP<-subset(BSTP,ID!="TLN488")
QuadT<-lm(PropWPup~Day+I(Day^2),data=BSTP[BSTP$Sp=="T",],na.action="na.omit")
summary(QuadT)
plot(QuadT)

CP1<-summary(QuadP)$coefficients["Day","Estimate"]
CP2<-summary(QuadP)$coefficients["I(Day^2)","Estimate"]
IntP<-summary(QuadP)$coefficients["(Intercept)","Estimate"]

CT1<-summary(QuadT)$coefficients["Day","Estimate"]
CT2<-summary(QuadT)$coefficients["I(Day^2)","Estimate"]
IntT<-summary(QuadT)$coefficients["(Intercept)","Estimate"]

PolyP<-polynomial(c(IntP,CP1,CP2))
dPolyP <- deriv(PolyP)
pxm <- solve(dPolyP)
d2PolyP <- deriv(dPolyP) #Run 2nd derivative test
Pmax=pxm[predict(d2PolyP, pxm) < 0]

PolyT<-polynomial(c(IntT,CT1,CT2))
dPolyT <- deriv(PolyT)
txm <- solve(dPolyT)
d2PolyT <- deriv(dPolyT) 
Tmax=txm[predict(d2PolyT, txm) < 0]

Diff<-Tmax-Pmax #238.596-167.706
Diff #70.89


#RANDOMIZATION LOOP - with quadratic functions
#====================
#Load Library
library(car)
library(lme4)
library(lmerTest)
library("plotrix")
library(MuMIn)
options(na.action = "na.fail") 
library(polynom)


BSTP<-read.csv("D:\\\\Documents\\\\Research Fun\\\\Slave Makers\\\\Data\\\\Brood Schedule\\\\Brood Schedule TP combined.csv")
BSTP<-subset(BSTP,ID!="TLN488")#Outlier of Tlong with 1 brood item and it was a pupa
BSTP<-subset(BSTP,Day>109) #Remove data before April because it was too cold for extranidal activity

#Caluclate some new values
BSTP$PropWPup<-(BSTP$Wpupae/BSTP$TotB) #Proportion of brood that are worker pupae
BSTP$NonWBrood<-(BSTP$TotB-BSTP$Wpupae) #Proportion of brood that are NOT worker pupae
BSTP$PropL<-(BSTP$L/BSTP$TotB) #Proportion of brood that are larvae

#Create new data table to hold permutation results
ResultsTable <- data.frame(matrix(vector(), 1500, 1, dimnames=list(c(), c("Difference"))), stringsAsFactors=F)

#Begining of Loop
for(i in 1:1500){
  #Randomize BSTP data with respect to species (PAM v Tlong)
  #Create random number vector to add to data table
  rand_N_vector <- runif(length(BSTP$Sp),0,1) #make X random numbers between 0 and 1 where x is the number of original observations
  #Add random vector to table
  R_BSTP <- cbind(rand_N_vector,BSTP)
  #Sort species based on this random vector column (holding brood in place)
  R_BSTP <- R_BSTP[order(rand_N_vector),]
  #Re-associate randomized species with brood proportions in original table
  N_BSTP <- BSTP #take original data frame and rename with N in front for "new"
  N_BSTP$Sp <- R_BSTP$Sp #re-assign old species ids with randomized spedcies ids
  
  #Define linear models for randomized Tlong and Pam 
  QuadP<-lm(PropWPup~Day+I(Day^2),data=N_BSTP[N_BSTP$Sp=="P",],na.action="na.omit")
  QuadT<-lm(PropWPup~Day+I(Day^2),data=N_BSTP[N_BSTP$Sp=="T",],na.action="na.omit")
  
  #Extract model coefficients and intercepts
  CP1<-summary(QuadP)$coefficients["Day","Estimate"]
  CP2<-summary(QuadP)$coefficients["I(Day^2)","Estimate"]
  IntP<-summary(QuadP)$coefficients["(Intercept)","Estimate"]
  CT1<-summary(QuadT)$coefficients["Day","Estimate"]
  CT2<-summary(QuadT)$coefficients["I(Day^2)","Estimate"]
  IntT<-summary(QuadT)$coefficients["(Intercept)","Estimate"]
  
  #Pull out day for which polynomial function is maximized
  PolyP<-polynomial(c(IntP,CP1,CP2)) #polynomial equation
  dPolyP <- deriv(PolyP) #first derivative
  pxm <- solve(dPolyP) #solve for day where PolyP is maximum
  if (pxm>=0 & pxm<=365) {
    Pmax=pxm
  } else { #Only do this if not a actual maxima within the range of days
    newdata=seq(0,365) 
    PP<-predict(PolyP,newdata=newdata)
    maxp=max(PP)
    PMaxDay=solve(PolyP,maxp) #the day where PolyP is maximized
    PMaxDay_i<-as.integer(PMaxDay) #Make integer so following greater than equal to is possible
    Pmax=PMaxDay[PMaxDay_i>=0 & PMaxDay_i<=365] #select the value that falls between 0 and 365
  }
 
  
   # d2PolyP <- deriv(dPolyP)  #find second derivative to identify local max
  # Pmax=pxm[predict(d2PolyP, pxm) < 0] #calculate day at this maximum for Pam
  # 
 
  PolyT<-polynomial(c(IntT,CT1,CT2)) #repeat above for Tlong
  dPolyT <- deriv(PolyT) #first derivative
  txm <- solve(dPolyT) #solve for day where PolyP is maximum
  if (txm>=0 & txm<=365) {
    Tmax=txm
  } else { #Only do this if not a actual maxima within the range of days
    newdata=seq(0,365) 
    PT<-predict(PolyT,newdata=newdata)
    maxt=max(PT)
    TMaxDay=solve(PolyT,maxt)
    TMaxDay_i<-as.integer(TMaxDay)
    Tmax=TMaxDay[TMaxDay_i>=0 & TMaxDay_i<=365]
  }
  

  
  #calculate difference between these values
  Diff<-Tmax-Pmax 
  
  #Add difference value for current permutation to results data table
  if ((length(Tmax)>0)&(length(Pmax)>0)){
    ResultsTable[i,1] <- Diff}
  else {i <- i-1}
}

write.csv(ResultsTable, "D://Documents//Research Fun//Slave Makers//Data//Brood Schedule//Permutation Results.csv")

#Calculate P-value
#Original Diff = 70.89
hist(ResultsTable$Difference)
greater_than_mydiff<-ResultsTable[ResultsTable$Difference>70.89,]
num_greater<-length(greater_than_mydiff)
p_value<-num_greater/1500
#p=0.02333333



#MODEL SELECTION

#Load data file
BSTP<-read.csv("....Brood Schedule TP Combined.csv",header=TRUE)
BSTP<-subset(BSTP,Day>109) #Remove data before April because it was too cold for extranidal activity

BSTP$PropWPup<-(BSTP$Wpupae/BSTP$TotB)
BSTP$NonWBrood<-(BSTP$TotB-BSTP$Wpupae)
BSTP$PropL<-(BSTP$L/BSTP$TotB)
BSTP<-subset(BSTP,Day>109)

#GLM with percent pupae as the response 
library(car)
library(lme4)
library(lmerTest)
library("plotrix")
library(RVAideMemoire)
library(MuMIn)
options(na.action = "na.fail") 

BSTP$Year<-as.factor(BSTP$Year)
BSTP_NoNA<-subset(BSTP,!is.na(Wpupae))
glm_bs<-glmer(cbind(Wpupae,NonWBrood)~scale(Day)+Sp+scale(Day)*Sp+(1|Year)+(1|Loc),family=binomial(link = "logit"),data = BSTP_NoNA)
glm_bs
summary(glm_bs)
# #Interaction effect means that the effect of day is species specific
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: cbind(Wpupae, NonWBrood) ~ scale(Day) + Sp + scale(Day) * Sp +      (1 | Year) + (1 | Loc)
# Data: BSTP_NoNA
# 
# AIC      BIC   logLik deviance df.resid 
# 8838.4   8864.3  -4413.2   8826.4      549 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -16.8260  -1.5229  -0.5078   1.5269  13.2251 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# Year   (Intercept) 3.564    1.888   
# Loc    (Intercept) 1.899    1.378   
# Number of obs: 555, groups:  Year, 4; Loc, 3
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -1.96966    1.23446  -1.596    0.111    
# scale(Day)     -0.36980    0.05575  -6.633 3.28e-11 ***
#   SpT             0.38604    0.03659  10.551  < 2e-16 ***
#   scale(Day):SpT  1.40838    0.06306  22.333  < 2e-16 ***

model.sel.bs<-dredge(glm_bs)
model.sel.bs
get.models(model.sel.bs,subset=1)
#Interaction effect means that the effect of day is species specific

summary(glm_bs)
#standardizing the beta coefficients
#Day
-0.36980*(0.05575*sqrtn)
#-0.4848131
#Species
0.38604*(0.03659*sqrtn) 
#0.3321676
#SpxDay Interaction
1.40838*(0.06306*sqrtn)
#2.088509