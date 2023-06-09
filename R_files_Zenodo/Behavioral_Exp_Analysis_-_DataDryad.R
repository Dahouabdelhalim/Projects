#Load R packages
library(plyr)
library(dplyr)
library(plotrix)
library(ggplot2)
library(data.table)
library(lme4)
library(lsmeans)
library(lmerTest)
library(survival)
library(survminer)
library(psa)
library(MANOVA.RM)

#Read in and filter raw data
Beetles<-read.csv("Behavioral Exp Data.csv",na.strings = "")
Beetles$Date.Paired<-as.POSIXct(strptime(paste(Beetles$Date.Paired), format = "%m/%d/%Y"))
#Remove any trials that failed due to innappropriately timed adding of larvae, missexing, or death of non-focal individual
Beetles<-Beetles[-which(Beetles$Reason.Failed=="User Error"),]

#Calculate sample Sizes
nrow(Beetles) #365 total trials
nrow(Beetles[which(!is.na(Beetles$Condition)),]) #165 that made it through behavioral observation stage
#Hot BP
nrow(Beetles[which(Beetles$Treatment=="HOT"&Beetles$Condition=="BP"),])
#Hot UPM
nrow(Beetles[which(Beetles$Treatment=="HOT"&Beetles$Condition=="UPM"),])
#Hot UPF
nrow(Beetles[which(Beetles$Treatment=="HOT"&Beetles$Condition=="UPF"),])
#Cold BP
nrow(Beetles[which(Beetles$Treatment=="COLD"&Beetles$Condition=="BP"),])
#Cold UPM
nrow(Beetles[which(Beetles$Treatment=="COLD"&Beetles$Condition=="UPM"),])
#Cold UPF
nrow(Beetles[which(Beetles$Treatment=="COLD"&Beetles$Condition=="UPF"),])
#Number of repeated measure individuals
nrow(Beetles[which(Beetles$Treatment=="COLD"&Beetles$M.Index==2),])
nrow(Beetles[which(Beetles$Treatment=="COLD"&Beetles$F.Index==2),])
nrow(Beetles[which(Beetles$Treatment=="HOT"&Beetles$M.Index==2),])
nrow(Beetles[which(Beetles$Treatment=="HOT"&Beetles$F.Index==2),])

#Flag repeated measures to control for breeding history (i.e., parents already bred successfully prior to focal trial)
#First flag any breeding attempts for individuals that follow successful breeding attempts
Beetles <- Beetles[order(as.Date(Beetles$Date.Paired)),] 
Beetles <- Beetles %>% 
  group_by(F.ID) %>% 
  mutate(F.BreedingHistoryInclusive = cumsum(Condition!="UPM"&!is.na(Condition)) > 0)
Beetles <- Beetles %>% 
  group_by(M.ID) %>% 
  mutate(M.BreedingHistoryInclusive = cumsum(Condition!="UPF"&!is.na(Condition)) > 0)
#Create Index column for male and female successful breeding attempts
Beetles$M.Index<-NA_integer_
Beetles[which(!is.na(Beetles$Condition)&is.na(Beetles$Reason.Failed)&Beetles$Condition=="UPM"|Beetles$Condition=="BP"),]$M.Index<-rowid(Beetles[which(!is.na(Beetles$Condition)&is.na(Beetles$Reason.Failed)&Beetles$Condition=="UPM"|Beetles$Condition=="BP"),]$M.ID)
Beetles$F.Index<-NA_integer_
Beetles[which(!is.na(Beetles$Condition)&is.na(Beetles$Reason.Failed)&Beetles$Condition=="UPF"|Beetles$Condition=="BP"),]$F.Index<-rowid(Beetles[which(!is.na(Beetles$Condition)&is.na(Beetles$Reason.Failed)&Beetles$Condition=="UPF"|Beetles$Condition=="BP"),]$F.ID)
#Now combine this information so that only attempts after the initial success are flagged
Beetles$BreedingHistory<-0
Beetles[which(Beetles$F.BreedingHistoryInclusive==TRUE&is.na(Beetles$F.Index)|Beetles$F.BreedingHistoryInclusive==TRUE&Beetles$F.Index>1),]$BreedingHistory<-1
Beetles[which(Beetles$M.BreedingHistoryInclusive==TRUE&is.na(Beetles$M.Index)|Beetles$M.BreedingHistoryInclusive==TRUE&Beetles$M.Index>1),]$BreedingHistory<-1
Beetles$BreedingHistory<-as.factor(Beetles$BreedingHistory)


####1 Fitness ##################################
#Survival Analysis
#merge eclosion and death dates
dates<-read.csv("names.csv")
dates$Eclosion.Date<-as.POSIXct(strptime(paste(dates$Eclosion.Date), format = "%m/%d/%Y"))
dates$Death.Date..if.applicable.<-as.POSIXct(strptime(paste(dates$Death.Date..if.applicable.), format = "%m/%d/%Y"))
dates$time<-NA
dates[which(!is.na(dates$Death.Date..if.applicable.)),]$time<-dates[which(!is.na(dates$Death.Date..if.applicable.)),]$Death.Date..if.applicable.-dates[which(!is.na(dates$Death.Date..if.applicable.)),]$Eclosion.Date
dates[which(is.na(dates$Death.Date..if.applicable.)),]$time<-as.POSIXct("2/15/21", format="%m/%d/%Y")-dates[which(is.na(dates$Death.Date..if.applicable.)),]$Eclosion.Date
dates$status<-NA
dates[which(!is.na(dates$Death.Date..if.applicable.)),]$status<-2
dates[which(is.na(dates$Death.Date..if.applicable.)),]$status<-1
survM<-dates
colnames(survM)<-c("M.ID", "M.Eclosion.Date","M.Death.Date","time","status")
survM<-merge(survM,Beetles[which(Beetles$Condition!="UPF"),c("M.ID","Treatment","Condition","Date.Paired", "M.Index","M.Direct.Care","M.Indirect.Care","M.Duration.on.carcass")],all.x=T)
survM<-survM[-which(is.na(survM$Treatment)),]
survM$age<-survM$Date.Paired-survM$M.Eclosion.Date
survM$M.TotalCare<-survM$M.Direct.Care+survM$M.Indirect.Care
survF<-dates
colnames(survF)<-c("F.ID", "F.Eclosion.Date","F.Death.Date","time","status")
survF<-merge(survF,Beetles[which(Beetles$Condition!="UPM"),c("F.ID","Treatment","Condition","Date.Paired", "F.Index","F.Direct.Care","F.Indirect.Care","F.Duration.on.carcass..days.")],all.x=T)
survF<-survF[-which(is.na(survF$Treatment)),]
survF$age<-survF$Date.Paired-survF$F.Eclosion.Date
survF$F.TotalCare<-survF$F.Direct.Care+survF$F.Indirect.Care
survM$Sex<-"M"
survF$Sex<-"F"
colnames(survF)<-c("ID","Eclosion.Date", "Death.Date","time","status","Treatment", "Condition","Date.Paired","Index","Direct.Care","Indirect.Care","Duration.Care","age","Total.Care","Sex")
colnames(survM)<-c("ID","Eclosion.Date", "Death.Date","time","status","Treatment", "Condition","Date.Paired","Index","Direct.Care","Indirect.Care","Duration.Care","age","Total.Care","Sex")
surv<-rbind(survF,survM)
surv<-surv[which(surv$Index==1),] #Keep only individuals' first breeding attempt
fit<-coxph(Surv(time,status)~Treatment+Sex,data=surv)
summary(fit) #Hazard ratios (HR) 
#exp(coef) = 2.0075, meaning likelihood of death in harsh is ~2x that in benign
#exp(coef) = 0.498, meaning beetles in harsh are ~50% more likely to die
confint(fit) #coef CIs
exp(confint(fit)) #HR CIs
#Males live shorter


#Next look at differences in potential and realized reproductive success.
#Potential Reproductive Output
#compare fecundity
Beetles$HatchingSuccess<-Beetles$No..Hatched.Larvae/Beetles$No..of.Eggs
t.test(Beetles[which(Beetles$Treatment=="HOT"),]$No..of.Eggs,Beetles[which(Beetles$Treatment=="COLD"),]$No..of.Eggs)
t.test(Beetles[which(Beetles$Treatment=="HOT"),]$No..Hatched.Larvae,Beetles[which(Beetles$Treatment=="COLD"),]$No..Hatched.Larvae)
#And fertility
t.test(Beetles[which(Beetles$Treatment=="HOT"),]$HatchingSuccess,Beetles[which(Beetles$Treatment=="COLD"),]$HatchingSuccess)


#For realized reproductive success, remove all trials that failed before the observation stage.
Beetles<-Beetles[which(!is.na(Beetles$Condition)&is.na(Beetles$Reason.Failed)),]

#Realized output (as a measure of no. larvae dispersded)
anova(lm(No..of.Larvae.Dispersed~Treatment+BreedingHistory,data=Beetles))
#As a measure of mean larval mass
anova(lm(data=Beetles, Mean.Larval.Mass~Treatment+BreedingHistory))



###2.Behavior##########################################

#Selection Gradients: are different behaviors implicated in success between treatments?
#Calculate metrics for cumulative care 
Beetles$Total.Direct.Care<-NA_integer_
Beetles[which(Beetles$Condition=="BP"),]$Total.Direct.Care<-Beetles[which(Beetles$Condition=="BP"),]$F.Direct.Care+Beetles[which(Beetles$Condition=="BP"),]$M.Direct.Care
Beetles[which(Beetles$Condition=="UPF"),]$Total.Direct.Care<-Beetles[which(Beetles$Condition=="UPF"),]$F.Direct.Care
Beetles[which(Beetles$Condition=="UPM"),]$Total.Direct.Care<-Beetles[which(Beetles$Condition=="UPM"),]$M.Direct.Care
Beetles$Total.Indirect.Care<-NA_integer_
Beetles[which(Beetles$Condition=="BP"),]$Total.Indirect.Care<-Beetles[which(Beetles$Condition=="BP"),]$F.Indirect.Care+Beetles[which(Beetles$Condition=="BP"),]$M.Indirect.Care
Beetles[which(Beetles$Condition=="UPF"),]$Total.Indirect.Care<-Beetles[which(Beetles$Condition=="UPF"),]$F.Indirect.Care
Beetles[which(Beetles$Condition=="UPM"),]$Total.Indirect.Care<-Beetles[which(Beetles$Condition=="UPM"),]$M.Indirect.Care
Beetles$M.TotalCare<-Beetles$M.Direct.Care+Beetles$M.Indirect.Care
Beetles$F.TotalCare<-Beetles$F.Direct.Care+Beetles$F.Indirect.Care
Beetles$Total.Care<-NA_integer_
Beetles[which(Beetles$Condition=="BP"),]$Total.Care<-Beetles[which(Beetles$Condition=="BP"),]$F.TotalCare+Beetles[which(Beetles$Condition=="BP"),]$M.TotalCare
Beetles[which(Beetles$Condition=="UPF"),]$Total.Care<-Beetles[which(Beetles$Condition=="UPF"),]$F.TotalCare
Beetles[which(Beetles$Condition=="UPM"),]$Total.Care<-Beetles[which(Beetles$Condition=="UPM"),]$M.TotalCare
#Days attended by >=1 parent:
Beetles$MaxDuration<-pmax(Beetles$F.Duration.on.carcass..days.,Beetles$M.Duration.on.carcass,na.rm=TRUE)

#To calculate selection gradients, we want to use linear regression of relative fitness measures
#First calculate relative fitness measures (Lande and Arnold 1983)
Beetles$No..of.Larvae.Dispersed<-factor2number(Beetles$No..of.Larvae.Dispersed)
NoLHot<-mean(Beetles[which(Beetles$Treatment=="HOT"),]$No..of.Larvae.Dispersed,na.rm=T)
NoLCold<-mean(Beetles[which(Beetles$Treatment=="COLD"),]$No..of.Larvae.Dispersed,na.rm=T)
Beetles$Rel.NoLarvae<-NA_real_
Beetles[which(Beetles$Treatment=="HOT"),]$Rel.NoLarvae<-Beetles[which(Beetles$Treatment=="HOT"),]$No..of.Larvae.Dispersed/NoLHot
Beetles[which(Beetles$Treatment=="COLD"),]$Rel.NoLarvae<-Beetles[which(Beetles$Treatment=="COLD"),]$No..of.Larvae.Dispersed/NoLCold
LMassHot<-mean(Beetles[which(Beetles$Treatment=="HOT"),]$Mean.Larval.Mass,na.rm=T)
LMassCold<-mean(Beetles[which(Beetles$Treatment=="COLD"),]$Mean.Larval.Mass,na.rm=T)
Beetles$Rel.LarvalMass<-NA_real_
Beetles[which(Beetles$Treatment=="HOT"),]$Rel.LarvalMass<-Beetles[which(Beetles$Treatment=="HOT"),]$Mean.Larval.Mass/LMassHot
Beetles[which(Beetles$Treatment=="COLD"),]$Rel.LarvalMass<-Beetles[which(Beetles$Treatment=="COLD"),]$Mean.Larval.Mass/LMassCold

#Now run a multiple regression of ALL behaviors on the trait. Use standardized behaviors.
#We will look first at the linear relationship, then the quadratic.
#Calculate squares of behaviors for quadratic model. We will need to divide coefficients by two.
Beetles$MaxDuration2<-(Beetles$MaxDuration)^2
Beetles$Total.Care2<-(Beetles$Total.Care)^2
Beetles$Total.Direct.Care2<-(Beetles$Total.Direct.Care)^2
Beetles$Total.Indirect.Care2<-(Beetles$Total.Indirect.Care)^2
Beetles$F.DaysAttending2<-(Beetles$F.Duration.on.carcass..days.)^2
Beetles$F.TotalCare2<-(Beetles$F.TotalCare)^2
Beetles$F.Direct.Care2<-(Beetles$F.Direct.Care)^2
Beetles$F.Indirect.Care2<-(Beetles$F.Indirect.Care)^2
Beetles$M.PropDaysAttending2<-(Beetles$M.Duration.on.carcass)^2
Beetles$M.TotalCare2<-(Beetles$M.TotalCare)^2
Beetles$M.Direct.Care2<-(Beetles$M.Direct.Care)^2
Beetles$M.Indirect.Care2<-(Beetles$M.Indirect.Care)^2


#For Hot v Cold
####No. of larvae#####
#Hot linear
summary(lm(Rel.NoLarvae~scale(MaxDuration)+scale(Total.Indirect.Care)+scale(Total.Direct.Care),data=Beetles[which(Beetles$Treatment=="HOT"),])) 
#Hot nonlinear
summary(lm(Rel.NoLarvae~scale(MaxDuration)+scale(MaxDuration2)+scale(Total.Indirect.Care)+scale(Total.Indirect.Care2)+scale(Total.Direct.Care)+scale(Total.Direct.Care2),data=Beetles[which(Beetles$Treatment=="HOT"),])) 
#Cold linear
summary(lm(Rel.NoLarvae~scale(MaxDuration)+scale(Total.Indirect.Care)+scale(Total.Direct.Care),data=Beetles[which(Beetles$Treatment=="COLD"),])) 
#Cold nonlinear
summary(lm(Rel.NoLarvae~scale(MaxDuration)+scale(MaxDuration2)+scale(Total.Indirect.Care)+scale(Total.Indirect.Care2)+scale(Total.Direct.Care)+scale(Total.Direct.Care2),data=Beetles[which(Beetles$Treatment=="COLD"),])) 
####Larval mass####
#Hot linear
summary(lm(Rel.LarvalMass~scale(MaxDuration)+scale(Total.Indirect.Care)+scale(Total.Direct.Care),data=Beetles[which(Beetles$Treatment=="HOT"),])) 
#Hot nonlinear
summary(lm(Rel.LarvalMass~scale(MaxDuration)+scale(MaxDuration2)+scale(Total.Indirect.Care)+scale(Total.Indirect.Care2)+scale(Total.Direct.Care)+scale(Total.Direct.Care2),data=Beetles[which(Beetles$Treatment=="HOT"),])) 
#Cold linear
summary(lm(Rel.LarvalMass~scale(MaxDuration)+scale(Total.Indirect.Care)+scale(Total.Direct.Care),data=Beetles[which(Beetles$Treatment=="COLD"),])) 
#Cold nonlinear
ColdCLQ<-summary(lm(Rel.LarvalMass~scale(MaxDuration)+scale(MaxDuration2)+scale(Total.Indirect.Care)+scale(Total.Indirect.Care2)+scale(Total.Direct.Care)+scale(Total.Direct.Care2),data=Beetles[which(Beetles$Treatment=="COLD"),])) 


#Examine effect of social condition and perform a priori comparisons of uniparental to biparental
TotalMaintenanceHot<-lmer(data=Beetles[which(Beetles$Treatment=="HOT"),], Total.Indirect.Care~Condition+(1|F.ID)+(1|M.ID))
anova(TotalMaintenanceHot)
lsmeans(TotalMaintenanceHot,pairwise~Condition)
TotalMaintenanceCold<-lmer(data=Beetles[which(Beetles$Treatment=="COLD"),], Total.Indirect.Care~Condition+(1|F.ID)+(1|M.ID))
anova(TotalMaintenanceCold)
lsmeans(TotalMaintenanceCold,pairwise~Condition)
TotalMaintenance<-lmer(data=Beetles, Total.Indirect.Care~Condition*Treatment+(1|F.ID)+(1|M.ID))
lsmeans(TotalMaintenance,pairwise~Treatment|Condition)

TotalProvisionHot<-lmer(data=Beetles[which(Beetles$Treatment=="HOT"),], Total.Direct.Care~Condition+(1|F.ID)+(1|M.ID))
anova(TotalProvisionHot)
lsmeans(TotalProvisionHot,pairwise~Condition)
TotalProvisionCold<-lmer(data=Beetles[which(Beetles$Treatment=="COLD"),], Total.Direct.Care~Condition+(1|F.ID)+(1|M.ID))
anova(TotalProvisionCold)
lsmeans(TotalProvisionCold,pairwise~Condition)
TotalProvision<-lmer(data=Beetles, Total.Direct.Care~Condition*Treatment+(1|F.ID)+(1|M.ID))
lsmeans(TotalProvision,pairwise~Treatment|Condition)

DaysAttendingHot<-lmer(data=Beetles[which(Beetles$Treatment=="HOT"),], MaxDuration~Condition+(1|F.ID)+(1|M.ID))
anova(DaysAttendingHot)
lsmeans(DaysAttendingHot,pairwise~Condition)
DaysAttendingCold<-lmer(data=Beetles[which(Beetles$Treatment=="COLD"),], MaxDuration~Condition+(1|F.ID)+(1|M.ID))
anova(DaysAttendingCold)
lsmeans(DaysAttendingCold,pairwise~Condition)
DaysAttending<-lmer(data=Beetles, MaxDuration~Condition*Treatment+(1|F.ID)+(1|M.ID))
lsmeans(DaysAttending,pairwise~Treatment|Condition)



#####3. Plasticity###########################################################

#Split the dataset by sex and retain only individuals w/ more than one observation
Males2x<-Beetles[which(Beetles$M.Index==2),]$M.ID
MaleRM<-as.data.frame(Beetles[Beetles$M.ID %in% Males2x & Beetles$M.Index==1 |Beetles$M.ID %in% Males2x & Beetles$M.Index==2,])
MaleRM<-MaleRM[which(!is.na(MaleRM$M.ID)),]

Females2x<-Beetles[which(Beetles$F.Index==2),]$F.ID
FemaleRM<-as.data.frame(Beetles[Beetles$F.ID %in% Females2x & Beetles$F.Index==1 |Beetles$F.ID %in% Females2x & Beetles$F.Index==2,])
FemaleRM<-FemaleRM[which(!is.na(FemaleRM$F.ID)),]

#Repeated-measures ANOVAs: perform separate repeated-measures ANOVAs for each sex and social condition
#HOT males
summary(aov(M.Indirect.Care~Condition+Error(M.ID/M.Index),data=MaleRM[which(MaleRM$Treatment=="HOT"),]))
summary(aov(M.Direct.Care~Condition+Error(M.ID/M.Index),data=MaleRM[which(MaleRM$Treatment=="HOT"),]))
summary(aov(M.Duration.on.carcass~Condition+Error(M.ID/M.Index),data=MaleRM[which(MaleRM$Treatment=="HOT"),]))

#COLD males
summary(aov(M.Indirect.Care~Condition+Error(M.ID/M.Index),data=MaleRM[which(MaleRM$Treatment=="COLD"),]))
summary(aov(M.Direct.Care~Condition+Error(M.ID/M.Index),data=MaleRM[which(MaleRM$Treatment=="COLD"),]))
summary(aov(M.Duration.on.carcass~Condition+Error(M.ID/M.Index),data=MaleRM[which(MaleRM$Treatment=="COLD"),]))

#HOT females
summary(aov(F.Indirect.Care~Condition+Error(F.ID/F.Index),data=FemaleRM[which(FemaleRM$Treatment=="HOT"),]))
summary(aov(F.Direct.Care~Condition+Error(F.ID/F.Index),data=FemaleRM[which(FemaleRM$Treatment=="HOT"),]))
summary(aov(F.Duration.on.carcass..days.~Condition+Error(F.ID/F.Index),data=FemaleRM[which(FemaleRM$Treatment=="HOT"),]))

#COLD females
summary(aov(F.Indirect.Care~Condition+Error(F.ID/F.Index),data=FemaleRM[which(FemaleRM$Treatment=="COLD"),]))
summary(aov(F.Direct.Care~Condition+Error(F.ID/F.Index),data=FemaleRM[which(FemaleRM$Treatment=="COLD"),]))
summary(aov(F.Duration.on.carcass..days.~Condition+Error(F.ID/F.Index),data=FemaleRM[which(FemaleRM$Treatment=="COLD"),]))


#Run a multivariate ANOVA (MANOVA) to quantify these effects between-treatments 
#Cannot perform Mauchly's test for sphericity because we only have two levels.
#anova_test(data = MaleRM, dv = M.Indirect.Care, wid = M.ID, within = Condition)
#We want to be able to directly compare behaviors between treatments, so we should scale duration of care for development time
#Calculate mean development times in benign v. harsh
Beetles$Date.Time.Larvae.Added<-as.POSIXct(strptime(paste(Beetles$Date.Time.Larvae.Added), format = "%m/%d/%Y"))
Beetles$Dispersal.Date<-as.POSIXct(strptime(paste(Beetles$Dispersal.Date), format = "%m/%d/%Y"))
Beetles$DevelopmentTime<-NA_integer_
Beetles[which(Beetles$No..of.Larvae.Dispersed>0),]$DevelopmentTime<-as.integer(Beetles[which(Beetles$No..of.Larvae.Dispersed>0),]$Dispersal.Date-Beetles[which(Beetles$No..of.Larvae.Dispersed>0),]$Date.Time.Larvae.Added)

MaleRM$M.ScaledDuration<-NA_real_
MaleRM[which(MaleRM$Treatment=="HOT"),]$M.ScaledDuration<- MaleRM[which(MaleRM$Treatment=="HOT"),]$M.Duration.on.carcass/mean(Beetles[which(Beetles$Treatment=="HOT"),]$DevelopmentTime,na.rm=T)
MaleRM[which(MaleRM$Treatment=="COLD"),]$M.ScaledDuration<- MaleRM[which(MaleRM$Treatment=="COLD"),]$M.Duration.on.carcass/mean(Beetles[which(Beetles$Treatment=="COLD"),]$DevelopmentTime,na.rm=T)
FemaleRM$F.ScaledDuration<-NA_real_
FemaleRM[which(FemaleRM$Treatment=="HOT"),]$F.ScaledDuration<- FemaleRM[which(FemaleRM$Treatment=="HOT"),]$F.Duration.on.carcass..days./mean(Beetles[which(Beetles$Treatment=="HOT"),]$DevelopmentTime,na.rm=T)
FemaleRM[which(FemaleRM$Treatment=="COLD"),]$F.ScaledDuration<- FemaleRM[which(FemaleRM$Treatment=="COLD"),]$F.Duration.on.carcass..days./mean(Beetles[which(Beetles$Treatment=="COLD"),]$DevelopmentTime,na.rm=T)

if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}
fitM<-multRM(cbind(M.Indirect.Care,M.Direct.Care,M.ScaledDuration)~Treatment*Condition,subject="M.ID",within="Condition",iter=200, data=MaleRM)
summary(fitM)
simCI(fitM,contrast="pairwise",type="Tukey")
fitF<-multRM(cbind(F.Indirect.Care,F.Direct.Care,F.ScaledDuration)~Treatment*Condition,subject="F.ID",within="Condition",iter=200, data=FemaleRM)
summary(fitF)
simCI(fitF,contrast="pairwise",type="Tukey")

