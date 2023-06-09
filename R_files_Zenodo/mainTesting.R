#################################################################### #
#
# This script is part of a school-testing project for Belgium.
#
# The simulation model is based on Torneri et al. (2020 BMC Medicine)
# 
# copyright: Torneri et al. SIMID 2021
#################################################################### #


## LOAD DEPENDENCIES  ----
################################## #

#load packages
library(foreach)
#library(doRNG)

# load functions
source("SimFunctionsTesting.R")
source("nCov_simulator_testing.R")


## COMMAND LINE ARGUMENTS ----
################################## #

args <- commandArgs(trailingOnly = TRUE)

# option to use testthat and/or debug arguments
if(length(args)==0 && exists('args_testthat')){
  args <- args_testthat
}

#working dir
out = args[1] 
#number of cores to run in parallel
cores = as.numeric(args[2]) 
#test.sens.as: asymptomatic test sensitivity
test.sens.as = as.numeric(args[3]) 
#test.delay.as: asymptomatic test delay (to obtain the result)
test.delay.as = strtoi(args[4])
#test.sens.s: symptomatic test sensitivity
test.sens.s = as.numeric(args[5])
#test.delay.s: symptomatic test delay (to obtain the result)
test.delay.s = strtoi(args[6])
#delay between case detection and start of screening
screening.delay = as.numeric(args[7]) 
#number of tests to be conduced per week, when screening
n.test.week = as.numeric(args[8])
#inter-class contact freq == proportion of intra-class contact freq
prop.lambda.w =  as.numeric(args[9])
#protocol
strategy = args[10] 
#timewindow in which to consider cases that could trigger a closure (of class/school)
timewindow.closure = as.numeric(args[11])
#number of cases (in the time window) that result in closing the SCHOOL
threshold.school = as.numeric(args[12]) 
#number of cases (in the time window) that result in closing a CLASS 
threshold.class = as.numeric(args[13]) 
#Basic reproduction number for symptomatic adults 
R.s = as.numeric(args[14])
#Proportion of children that are immune
prop.immune.ch = as.numeric(args[15])
#Proportion of adults that are immune
prop.immune.ad = as.numeric(args[16])

#Asymptomatic infectiousness, relative to sympto
alpha.as = as.numeric(args[17])

#child seed
child.seed.arg <- as.numeric(args[18])


# Compliance
compliance = as.numeric(args[19])

#Nreas
nrs<- as.numeric(args[20])

#Seeding Time
seeding.time<-as.numeric(args[21])

#Probability that a child is symptomatic
rho.ch<-as.numeric(args[22])

#Variant
variant<-as.character(args[23])

#School Size
school.size<-as.numeric(args[24])

#Printing Network info
netw<-as.numeric(args[25])


# Number of stochastic realisations
pos <- 26
nSim <- ifelse(length(args)>=pos,as.numeric(args[pos]),100)
# Random number generator seed
pos <- 27
rnSeed<- ifelse(length(args)>=pos,as.numeric(args[pos]),14022020)

# To log all contacts
pos <- 28 
bool.log.contacts <- ifelse(length(args)>=pos,as.logical(args[pos]),FALSE)



# cat(",out=",out)
# cat(",cores=",cores)
# cat(",test.sens.as=",test.sens.as)
# cat(",test.delay.as=",test.delay.as)
# cat(",test.sens.s=",test.sens.s)
# cat(",test.delay.s=",test.delay.s)
# cat(",screening.delay=",screening.delay)
# cat(",n.test.week=",n.test.week)
# cat(",prop.lambda.w=",prop.lambda.w)
# cat(",strategy=",strategy)
# cat(",timewindow.closure=",timewindow.closure)
# cat(",threshold.school=",threshold.school)
# cat(",threshold.class=",threshold.class)
# #TODO: add nSim and rnSeed?

#print(paste0(out,cores,test.sens.as,test.delay.as,test.sens.s,test.delay.s,screening.delay,n.test.week,strategy))


## DEFAULT PARAMETERS           ----
################################## #

#pmf inferred from data
#data from: https://www.agodi.be/nieuwe-omkadering-basisonderwijs
#starts at class size 1
#!!! primary class sizes, best proxy for secondary class sizes, I'm afraid
zero = .Machine$double.eps
class_size_pmf = c(zero, 0.000452284034373587, zero, 0.000452284034373587, zero, zero, zero, zero, 0.000904568068747173, 0.00135685210312076, 0.00180913613749435, 
0.00180913613749435, 0.00814111261872456, 0.0235187697874265, 
0.0398009950248756, 0.052464947987336, 0.0687471732247852, 0.106286748077793, 
0.180009045680687, 0.26232473993668, 0.251922207146088) 

#n.teachers: number of teachers
#compute as proportion of the pupils, from:
#https://www.vlaanderen.be/publicaties/vlaams-onderwijs-in-cijfers
teacher.pupil.ratio <- 1/9

#TODO: simple seeding for now
#TODO: use Vittoria's seed function
seeds.ad = 1
seeds.ch = child.seed.arg 

#Repetitive seeding
reseeding.ch <- child.seed.arg
reseeding.ad <- 1

#lambda.b: contact frequency between classes
#data from: SOCRATES - Survey Belgium 2010 - contacts for children [6,12] years old (ALL CONTACTS EXCLUDED HOLIDAYS) - within: more than 1 hr - between total contacts - within
lambda.w = 6.62 # contact frequency within classes
lambda.b.tot <- 2.5 # contact frequency between classes in a non-pandemic situation
lambda.b = lambda.b.tot * prop.lambda.w #prop.lamba is the proportion of between contacts respect to the non-pandemic scenario

#
R.ctcs<-12.16 #(all cts excluded holidays) 

#rho.ch: prob. that child is symptomatic
#rho.ch = 1-.8
#rho.ad: prob. that adult is symptomatic
rho.ad = 1-.31

#suscep.ch: child susceptibility, relative to adult susceptibility
suscep.ch = .5

#pdiagn.ch: probability to TEST a child, when symptomatic
pdiagn.ch = .3
#pdiagn.ad: probability to TEST a adult, when symptomatic
pdiagn.ad = .5

#t.detection: time it takes before an indvidual becomes detectable
#2 days, same as in the antiviral paper
#TODO: check whether this is indeed similar to Vittoria's exposure time
t.detection = 2
#t.isolation: duration of the isolation
t.isolation = 10
#t.stop: duration of the simulation
t.stop<- 100



total.seeds<-(seeds.ch+seeds.ad)+floor(t.stop/seeding.time)*(reseeding.ad+reseeding.ch)


## PRE-PROCESSING    ----
################################## #

# set RNG stream
set.seed(rnSeed)

# sample population
classes.comp = c()
pupils = 0
while(pupils < school.size) {
  class = sample(x=1:length(class_size_pmf), prob=class_size_pmf, size=1)
  classes.comp = c(classes.comp, class)
  pupils = sum(classes.comp)
}

#n.teachers: number of teachers
n.teachers <- round(pupils*teacher.pupil.ratio, 0)

#total.considered.population
n<-pupils+n.teachers


##################   TEMPORARY   ########### #
## Parameters to test the reproduction number
# seeds.ad = 1
# seeds.ch = 10
# print(paste0('infected seeds ch/ad: ',seeds.ch,'/',seeds.ad))

#Repetitive seeding
#reseeding.ch <- 0
#reseeding.ad <- 0
#print('no additional seeding')

#suscep.ch = 1
#prop.immune <- 0
#print('no immunity')
##################   TEMPORARY   ########### #


## RUN SIMULATION    ----
################################## #

  if (strategy=="PI"){
    screening.delay <- Inf
    threshold.class <- Inf
    n.test.week     <- 0 
    reg.screening.type = NA
  }
  
  if (strategy=="RS"){
    reg.screening.type = NA
    n.reascr<-nrs
  }

if (strategy=="SI"){
  reg.screening.type = NA
}

  
  if (strategy=="RS_A"){
    reg.screening.type = "A"
  }
  
  if (strategy=="RS_B"){
    reg.screening.type = "B"
  }
  
  if (strategy=="RS_C"){
    reg.screening.type = "C"
  }

# initialise variable to aggregate output
epi.outbreak<-list()

for(i_sim in 1:nSim) {
  #print(i_sim)
  epi.outbreak[[i_sim]]<-nCov.simulator.Testing(n.teachers = n.teachers, lambda.b = lambda.b, lambda.w = lambda.w, rho.ch = rho.ch, rho.ad = rho.ad, alpha.as = alpha.as, R.s = R.s, classes.comp = classes.comp, seeds.ad = seeds.ad,seeds.ch = seeds.ch,test.sens.as = test.sens.as, test.sens.s = test.sens.s,test.delay.s = test.delay.s,test.delay.as = test.delay.as ,suscep.ch = suscep.ch,pdiagn.ch = pdiagn.ch,pdiagn.ad = pdiagn.ad, t.isolation = t.isolation, t.detection = t.detection, timewindow.closure = timewindow.closure, threshold.school = threshold.school,  prop.immune.ch=prop.immune.ch,  prop.immune.ad=prop.immune.ad, reseeding.ch=reseeding.ch, reseeding.ad=reseeding.ad, t.stop=t.stop, screening.delay = screening.delay, threshold.class = threshold.class, n.test.week = n.test.week,reg.screening.type = reg.screening.type,strategy=strategy, bool.log.contacts = bool.log.contacts,R.ctcs=R.ctcs, compliance=compliance, nrs=nrs, seeding.time=seeding.time)
} # end for-loop: nSim

finalSize<-NULL
finalSizeProp<-NULL
not.extinct<-NULL
n.tests<-NULL

n.schooldayslost <- NULL
n.classclosure <- NULL
n.schoolclosure <- NULL
n.detectedpositive <- NULL

n.immune      <- NULL
n.susceptible <- NULL
n.severe1     <- NULL 
n.severe2     <- NULL

for (i_sim in 1:nSim){
  total.seeds<-length(which(epi.outbreak[[i_sim]]$status.matrix$IndexCase==1))
  finalSize[i_sim] <-length(c(which(epi.outbreak[[i_sim]]$status.matrix$infected==-1),which(epi.outbreak[[i_sim]]$status.matrix$infected==1)))
  finalSizeProp[i_sim]<-(finalSize[i_sim]-total.seeds)/n
  n.tests      <-c(n.tests, epi.outbreak[[i_sim]]$n.tests)
  n.schooldayslost <- c(n.schooldayslost, epi.outbreak[[i_sim]]$school.days.losts/n)
  n.detectedpositive<- c(n.detectedpositive, epi.outbreak[[i_sim]]$n.detectedpositive)
  n.schoolclosure<- c(n.schoolclosure, epi.outbreak[[i_sim]]$n.schoolclosure)
  n.classclosure<- c(n.classclosure, epi.outbreak[[i_sim]]$n.classclosure)
  if (finalSize[i_sim]>round(n*0.1)){not.extinct<-c(not.extinct,i_sim)}
  n.immune      <- rbind(n.immune,sum(epi.outbreak[[i_sim]]$status.matrix[,1]==-2))
  n.susceptible <- rbind(n.susceptible,sum(epi.outbreak[[i_sim]]$status.matrix[,1]==0))
  n.severe1     <- c(n.severe1,sum(epi.outbreak[[i_sim]]$time.events[,2]==1,na.rm=T))
  n.severe2     <- c(n.severe2,sum(epi.outbreak[[i_sim]]$time.events[,2]==2,na.rm=T))
}
FinSize<-finalSize
FinalSizeProp<-finalSizeProp
n.tests<-n.tests

df <- data.frame(FinSize   = FinSize,
                 FinSizeProp = FinalSizeProp,
                 n.tests   = n.tests,
                 n.immune  = n.immune, 
                 n.susceptible = n.susceptible,
                 n.severe1 = n.severe1,
                 n.severe2 = n.severe2,
                 n.schooldayslost= n.schooldayslost,
                 n.detectedpositive = n.detectedpositive,
                 n.classclosure = n.classclosure,
                 n.schoolclosure = n.schoolclosure)
df
write.csv(df, paste0(out,"out.csv"))

if (netw==1){
  for (i in 1:nSim) {
    ctc.betw<-rep(NA,n)
    temp.ctc.betw<-epi.outbreak[[i]]$between.classes.contacts
    for (j in 1:length(temp.ctc.betw)){
      if (length(temp.ctc.betw[[j]][-1])>0){
        ctc.betw[j]<-paste(temp.ctc.betw[[j]][-1],collapse = ",")
      }
    }
    temp.status.matrix<-epi.outbreak[[i]]$status.matrix
    dfn<-data.frame(ID=1:n,
                    StudentClass=temp.status.matrix$Category,
                    TeacherClass=temp.status.matrix$TeacherAssign,
                    IsIndex=temp.status.matrix$IndexCase,
                    InfectionTime=temp.status.matrix$time.of.infection,
                    InfectorID=temp.status.matrix$infector,
                    SymptOnsDate=temp.status.matrix$TimeSymptomOnset,
                    DayPosTest=temp.status.matrix$TimePosTest,
                    Compliance=temp.status.matrix$Compliance,
                    BetweenClassesCtc=ctc.betw)
    write.csv(dfn, paste0(out,"InfNetw_Sim",i,".csv"))
  }
  
}



# save all
name<-paste(out,"epi_outbreak.rds", sep = "")
saveRDS(epi.outbreak, file = name)

