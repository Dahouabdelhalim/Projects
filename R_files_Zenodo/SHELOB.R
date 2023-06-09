#Load data
#Rodents.csv, AntsNum.csv, WaspsNum.csv, Pollinators.csv
Rodents <- read.csv("Rodents.csv")
Ants <- read.csv("Ants.csv")
Wasps <- read.csv("Wasps.csv")
Pollinators <- read.csv("Pollinators.csv")

#Data setup
#separate frequencies from treatment number
AntsNum <- Ants[,1]
WaspsNum <- Wasps[,1]
TreatmentA <- Ants[,2]
TreatmentW <- Wasps[,2]

#transform rats, mice, ants, and wasps into binomial
RatsNum <- table(Rodents$Treatment, Rodents$Rat)
RatsNum <- cbind(RatsNum, RatsNum[,1]+RatsNum[,2])
RatsNum <- RatsNum[,-1]
colnames(RatsNum) <- c("y", "N")

MiceNum <- table(Rodents$Treatment, Rodents$Mouse)
MiceNum <- cbind(MiceNum, MiceNum[,1]+MiceNum[,2])
MiceNum <- MiceNum[,-1]
colnames(MiceNum) <- c("y", "N")

TreatmentM <- 1:36
TreatmentR <- 1:36

#length of the SHELOBator datasets
nR <- length(RatsNum) 
nA <- length(AntsNum)
nW <- length(WaspsNum)

#Transform Ants and Wasps into binomial
AntsNum <- cbind(AntsNum, rep(194, nA)) #194 is max number of ants observed
WaspsNum <- cbind(WaspsNum, rep(57, nW)) #57 is max number of yellowjackets observed

#Set up the pollinator data
PlantSpp <- factor(Pollinators[,2])
PlantSpp <- as.numeric(PlantSpp)

Trials2 <- as.numeric(Pollinators[,3])
Trials <- mat.or.vec(nr=60, nc=6)
Trials[,1] <- Trials2[1:60]
Trials[,2] <- Trials2[61:120]
Trials[,3] <- Trials2[121:180]
Trials[,4] <- Trials2[181:240]
Trials[,5] <- Trials2[241:300]
Trials[,6] <- Trials2[301:360]
#convert NAs to 0's
Trials[which(is.na(Trials[,1])),1] <- 0
Trials[which(is.na(Trials[,2])),2] <- 0
Trials[which(is.na(Trials[,3])),3] <- 0
Trials[which(is.na(Trials[,4])),4] <- 0
Trials[which(is.na(Trials[,5])),5] <- 0
Trials[which(is.na(Trials[,6])),6] <- 0
#then back to NAs
Trials[which(Trials[,1]==0),1] <- NA
Trials[which(Trials[,2]==0),2] <- NA
Trials[which(Trials[,3]==0),3] <- NA
Trials[which(Trials[,4]==0),4] <- NA
Trials[which(Trials[,5]==0),5] <- NA
Trials[which(Trials[,6]==0),6] <- NA
Trials3 <- Trials[,-5]


SocialBeesTemp <- as.numeric(Pollinators[,4])
SocialBees <- mat.or.vec(nr=60, nc=6)
SocialBees[,1] <- SocialBeesTemp[1:60]
SocialBees[,2] <- SocialBeesTemp[61:120]
SocialBees[,3] <- SocialBeesTemp[121:180]
SocialBees[,4] <- SocialBeesTemp[181:240]
SocialBees[,5] <- SocialBeesTemp[241:300]
SocialBees[,6] <- SocialBeesTemp[301:360]
for(j in 1:6){
  for(i in 1:60){
    SocialBees[i,j] <- ifelse(is.na(SocialBees[i,j])==TRUE, 0, SocialBees[i,j])
  }
}
LepidopTemp <- as.numeric(Pollinators[,5])
Lepidop <- mat.or.vec(nr=60, nc=6)
Lepidop[,1] <- LepidopTemp[1:60]
Lepidop[,2] <- LepidopTemp[61:120]
Lepidop[,3] <- LepidopTemp[121:180]
Lepidop[,4] <- LepidopTemp[181:240]
Lepidop[,5] <- LepidopTemp[241:300]
Lepidop[,6] <- LepidopTemp[301:360]
for(j in 1:6){
  for(i in 1:60){
    Lepidop[i,j] <- ifelse(is.na(Lepidop[i,j])==TRUE, 0, Lepidop[i,j])
  }
}
NonSDipteraTemp <- as.numeric(Pollinators[,6])
NonSDiptera <- mat.or.vec(nr=60, nc=6)
NonSDiptera[,1] <- NonSDipteraTemp[1:60]
NonSDiptera[,2] <- NonSDipteraTemp[61:120]
NonSDiptera[,3] <- NonSDipteraTemp[121:180]
NonSDiptera[,4] <- NonSDipteraTemp[181:240]
NonSDiptera[,5] <- NonSDipteraTemp[241:300]
NonSDiptera[,6] <- NonSDipteraTemp[301:360]
for(j in 1:6){
  for(i in 1:60){
    NonSDiptera[i,j] <- ifelse(is.na(NonSDiptera[i,j])==TRUE, 0, NonSDiptera[i,j])
  }
}
NonVespWaspsTemp <- as.numeric(Pollinators[,7])
NonVespWasps <- mat.or.vec(nr=60, nc=6)
NonVespWasps[,1] <- NonVespWaspsTemp[1:60]
NonVespWasps[,2] <- NonVespWaspsTemp[61:120]
NonVespWasps[,3] <- NonVespWaspsTemp[121:180]
NonVespWasps[,4] <- NonVespWaspsTemp[181:240]
NonVespWasps[,5] <- NonVespWaspsTemp[241:300]
NonVespWasps[,6] <- NonVespWaspsTemp[301:360]
for(j in 1:6){
  for(i in 1:60){
    NonVespWasps[i,j] <- ifelse(is.na(NonVespWasps[i,j])==TRUE, 0, NonVespWasps[i,j])
  }
}
SolitaryBeesTemp <- as.numeric(Pollinators[,8])
SolitaryBees <- mat.or.vec(nr=60, nc=6)
SolitaryBees[,1] <- SolitaryBeesTemp[1:60]
SolitaryBees[,2] <- SolitaryBeesTemp[61:120]
SolitaryBees[,3] <- SolitaryBeesTemp[121:180]
SolitaryBees[,4] <- SolitaryBeesTemp[181:240]
SolitaryBees[,5] <- SolitaryBeesTemp[241:300]
SolitaryBees[,6] <- SolitaryBeesTemp[301:360]
for(j in 1:6){
  for(i in 1:60){
    SolitaryBees[i,j] <- ifelse(is.na(SolitaryBees[i,j])==TRUE, 0, SolitaryBees[i,j])
  }
}
SyrphidTemp <- as.numeric(Pollinators[,9])
Syrphid <- mat.or.vec(nr=60, nc=6)
Syrphid[,1] <- SyrphidTemp[1:60]
Syrphid[,2] <- SyrphidTemp[61:120]
Syrphid[,3] <- SyrphidTemp[121:180]
Syrphid[,4] <- SyrphidTemp[181:240]
Syrphid[,5] <- SyrphidTemp[241:300]
Syrphid[,6] <- SyrphidTemp[301:360]
for(j in 1:6){
  for(i in 1:60){
    Syrphid[i,j] <- ifelse(is.na(Syrphid[i,j])==TRUE, 0, Syrphid[i,j])
  }
}

SocialBees <- SocialBees[,-5]
NonSDiptera <- NonSDiptera[,-5]
NonVespWasps <- NonVespWasps[,-5]

T1 <- 1:60
T1 <- T1[-which(is.na(Trials[,1]))]

T2 <- 1:60
T2 <- T2[-which(is.na(Trials[,2]))]

T3 <- 1:60
T3 <- T3[-which(is.na(Trials[,3]))]

T4 <- 1:60
T4 <- T4[-which(is.na(Trials[,4]))]

T5 <- 1:60
T5 <- T5[-which(is.na(Trials[,5]))]

T52 <- 1:60
T52 <- T52[-which(is.na(Trials3[,5]))]

T6 <- 1:60
T6 <- T6[-which(is.na(Trials[,6]))]


#initial values for all thetas - note that mice and rats are different from ants and wasps
BR <- rep(0.6, 36) #mice
IR <- rep(0.6, 36) #mice
HR <- rep(0.6, 36) #mice
KR <- rep(0.15, 36) #rats
OR <- rep(0.15, 36) #rats
NR <- rep(0.15, 36) #rats
B <- rep(0.05, 60) #ants
I <- rep(0.05, 60) #ants
H <- rep(0.05, 60) #ants
K <- rep(0.1, 60) #wasps
O <- rep(0.1, 60) #wasps
N <- rep(0.1, 60) #wasps
#initial values for all slopes in the pollinator level
S <- rep(0, 6)
E <- rep(0, 6)
R <- rep(0, 6)
S5 <- rep(0, 5)
E5 <- rep(0, 5)
R5 <- rep(0, 5)



library(R2jags)
library(R2WinBUGS)

#FULL MODEL


#start the model here
#start the clock
ptm <- proc.time()
SHELOB <- function(){
  
  #SHELOBATOR LEVEL
  #Predators
  #Rats and Mice
  
  #loop through each plot/block combination
  for(h in 1:36){
    IntmeanM[h] ~ dunif(0.00000001, 1)
    IntmeanR[h] ~ dunif(0.00000001, 1)
    
    MiceNum[h,1] ~ dbin(IntmeanM[TreatmentM[h]], MiceNum[h,2])
    RatsNum[h,1] ~ dbin(IntmeanR[TreatmentR[h]], RatsNum[h,2])
  }
  
  #Ants and Wasps
  for(k in 1:nA){
    AntsNum[k,1] ~ dbin(meanA[TreatmentA[k]], AntsNum[k,2])
  }
  
  for(m in 1:nW){
    WaspsNum[m,1] ~ dbin(meanW[TreatmentW[m]], WaspsNum[m,2])
  }
  
  #loop through each treatment
  for(i in 1:60){
    meanA[i] ~ dunif(0.00000001, 1)
    
    meanW[i] ~ dunif(0.00000001, 1)
  }
  
  #use C values for A and V plots for Rats and Mice - ALL starts at 13, C starts at 25, R starts at 37, V starts at 49
  meanM[13:48] <- IntmeanM
  meanR[13:48] <- IntmeanR
  for(m in 1:12){
    meanM[m] <- IntmeanM[m+12]
    meanR[m] <- IntmeanR[m+12]
    meanM[m+48] <- IntmeanM[m+12]
    meanR[m+48] <- IntmeanR[m+12]
  }
  
  
  #POLLINATOR LEVEL
  for(j in 1:6){
    B0[j] ~ dnorm(0, 0.001)
    D0[j] ~ dnorm(0, 0.001)
    F0[j] ~ dnorm(0, 0.001)
    
    B1[j] ~ dnorm(0, 0.001)
    B2[j] ~ dnorm(0, 0.001)
    B3[j] ~ dnorm(0, 0.001)
    B4[j] ~ dnorm(0, 0.001)
    
    D1[j] ~ dnorm(0, 0.001)
    D2[j] ~ dnorm(0, 0.001)
    D3[j] ~ dnorm(0, 0.001)
    D4[j] ~ dnorm(0, 0.001)
    
    F1[j] ~ dnorm(0, 0.001)
    F2[j] ~ dnorm(0, 0.001)
    F3[j] ~ dnorm(0, 0.001)
    F4[j] ~ dnorm(0, 0.001)
  }
  
  #loop through plot/block combination)
  for(i in 1:length(T1)){
    Lepidop[T1[i],1] ~ dbin(meanLep[i,1], Trials[T1[i],1])
    SolitaryBees[T1[i],1] ~ dbin(meanSolBees[i,1], Trials[T1[i],1])
    Syrphid[T1[i],1] ~ dbin(meanSyrph[i,1], Trials[T1[i],1])
    
    logit(meanLep[i,1]) <- D0[1] + D1[1]*meanM[T1[i]] + D2[1]*meanR[T1[i]] + D3[1]*meanA[T1[i]] + D4[1]*meanW[T1[i]]
    logit(meanSolBees[i,1]) <- B0[1] + B1[1]*meanM[T1[i]] + B2[1]*meanR[T1[i]] + B3[1]*meanA[T1[i]] + B4[1]*meanW[T1[i]] 
    logit(meanSyrph[i,1]) <- F0[1] + F1[1]*meanM[T1[i]] + F2[1]*meanR[T1[i]] + F3[1]*meanA[T1[i]] + F4[1]*meanW[T1[i]] 
  }
  for(i in 1:length(T2)){
    Lepidop[T2[i],2] ~ dbin(meanLep[i,2], Trials[T2[i],2])
    SolitaryBees[T2[i],2] ~ dbin(meanSolBees[i,2], Trials[T2[i],2])
    Syrphid[T2[i],2] ~ dbin(meanSyrph[i,2], Trials[T2[i],2])
    
    logit(meanLep[i,2]) <- D0[2] + D1[2]*meanM[T2[i]] + D2[2]*meanR[T2[i]] + D3[2]*meanA[T2[i]] + D4[2]*meanW[T2[i]]
    logit(meanSolBees[i,2]) <- B0[2] + B1[2]*meanM[T2[i]] + B2[2]*meanR[T2[i]] + B3[2]*meanA[T2[i]] + B4[2]*meanW[T2[i]] 
    logit(meanSyrph[i,2]) <- F0[2] + F1[2]*meanM[T2[i]] + F2[2]*meanR[T2[i]] + F3[2]*meanA[T2[i]] + F4[2]*meanW[T2[i]] 
  }
  for(i in 1:length(T3)){
    Lepidop[T3[i],3] ~ dbin(meanLep[i,3], Trials[T3[i],3])
    SolitaryBees[T3[i],3] ~ dbin(meanSolBees[i,3], Trials[T3[i],3])
    Syrphid[T3[i],3] ~ dbin(meanSyrph[i,3], Trials[T3[i],3])
    
    logit(meanLep[i,3]) <- D0[3] + D1[3]*meanM[T3[i]] + D2[3]*meanR[T3[i]] + D3[3]*meanA[T3[i]] + D4[3]*meanW[T3[i]]
    logit(meanSolBees[i,3]) <- B0[3] + B1[3]*meanM[T3[i]] + B2[3]*meanR[T3[i]] + B3[3]*meanA[T3[i]] + B4[3]*meanW[T3[i]] 
    logit(meanSyrph[i,3]) <- F0[3] + F1[3]*meanM[T3[i]] + F2[3]*meanR[T3[i]] + F3[3]*meanA[T3[i]] + F4[3]*meanW[T3[i]] 
  }
  
  for(i in 1:length(T4)){
    Lepidop[T4[i],4] ~ dbin(meanLep[i,4], Trials[T4[i],4])
    SolitaryBees[T4[i],4] ~ dbin(meanSolBees[i,4], Trials[T4[i],4])
    Syrphid[T4[i],4] ~ dbin(meanSyrph[i,4], Trials[T4[i],4])
    
    logit(meanLep[i,4]) <- D0[4] + D1[4]*meanM[T4[i]] + D2[4]*meanR[T4[i]] + D3[4]*meanA[T4[i]] + D4[4]*meanW[T4[i]]
    logit(meanSolBees[i,4]) <- B0[4] + B1[4]*meanM[T4[i]] + B2[4]*meanR[T4[i]] + B3[4]*meanA[T4[i]] + B4[4]*meanW[T4[i]] 
    logit(meanSyrph[i,4]) <- F0[4] + F1[4]*meanM[T4[i]] + F2[4]*meanR[T4[i]] + F3[4]*meanA[T4[i]] + F4[4]*meanW[T4[i]] 
  }
  
  for(i in 1:length(T5)){
    Lepidop[T5[i],5] ~ dbin(meanLep[i,5], Trials[T5[i],5])
    SolitaryBees[T5[i],5] ~ dbin(meanSolBees[i,5], Trials[T5[i],5])
    Syrphid[T5[i],5] ~ dbin(meanSyrph[i,5], Trials[T5[i],5])
    
    logit(meanLep[i,5]) <- D0[5] + D1[5]*meanM[T5[i]] + D2[5]*meanR[T5[i]] + D3[5]*meanA[T5[i]] + D4[5]*meanW[T5[i]]
    logit(meanSolBees[i,5]) <- B0[5] + B1[5]*meanM[T5[i]] + B2[5]*meanR[T5[i]] + B3[5]*meanA[T5[i]] + B4[5]*meanW[T5[i]] 
    logit(meanSyrph[i,5]) <- F0[5] + F1[5]*meanM[T5[i]] + F2[5]*meanR[T5[i]] + F3[5]*meanA[T5[i]] + F4[5]*meanW[T5[i]] 
  }
  
  for(i in 1:length(T6)){
    Lepidop[T6[i],6] ~ dbin(meanLep[i,6], Trials[T6[i],6])
    SolitaryBees[T6[i],6] ~ dbin(meanSolBees[i,6], Trials[T6[i],6])
    Syrphid[T6[i],6] ~ dbin(meanSyrph[i,6], Trials[T6[i],6])
    
    logit(meanLep[i,6]) <- D0[6] + D1[6]*meanM[T6[i]] + D2[6]*meanR[T6[i]] + D3[6]*meanA[T6[i]] + D4[6]*meanW[T6[i]]
    logit(meanSolBees[i,6]) <- B0[6] + B1[6]*meanM[T6[i]] + B2[6]*meanR[T6[i]] + B3[6]*meanA[T6[i]] + B4[6]*meanW[T6[i]] 
    logit(meanSyrph[i,6]) <- F0[6] + F1[6]*meanM[T6[i]] + F2[6]*meanR[T6[i]] + F3[6]*meanA[T6[i]] + F4[6]*meanW[T6[i]] 
  }
  
  
  for(j in 1:5){
    C0[j] ~ dnorm(0, 0.001)
    E0[j] ~ dnorm(0, 0.001)
    G0[j] ~ dnorm(0, 0.001)
    
    C1[j] ~ dnorm(0, 0.001)
    C2[j] ~ dnorm(0, 0.001)
    C3[j] ~ dnorm(0, 0.001)
    C4[j] ~ dnorm(0, 0.001)
    
    E1[j] ~ dnorm(0, 0.001)
    E2[j] ~ dnorm(0, 0.001)
    E3[j] ~ dnorm(0, 0.001)
    E4[j] ~ dnorm(0, 0.001)
    
    G1[j] ~ dnorm(0, 0.001)
    G2[j] ~ dnorm(0, 0.001)
    G3[j] ~ dnorm(0, 0.001)
    G4[j] ~ dnorm(0, 0.001)
  }
  
  for(i in 1:length(T1)){
    SocialBees[T1[i],1] ~ dbin(meanSocBees[i,1], Trials3[T1[i],1])
    NonSDiptera[T1[i],1] ~ dbin(meanNonSDip[i,1], Trials3[T1[i],1])
    NonVespWasps[T1[i],1] ~ dbin(meanNonVWasp[i,1], Trials3[T1[i],1])
    
    logit(meanSocBees[i,1]) <- C0[1] + C1[1]*meanM[T1[i]] + C2[1]*meanR[T1[i]] + C3[1]*meanA[T1[i]] + C4[1]*meanW[T1[i]] 
    logit(meanNonSDip[i,1]) <- E0[1] + E1[1]*meanM[T1[i]] + E2[1]*meanR[T1[i]] + E3[1]*meanA[T1[i]] + E4[1]*meanW[T1[i]] 
    logit(meanNonVWasp[i,1]) <- G0[1] + G1[1]*meanM[T1[i]] + G2[1]*meanR[T1[i]]+ G3[1]*meanA[T1[i]] + G4[1]*meanW[T1[i]] 
  }
  for(i in 1:length(T2)){
    SocialBees[T2[i],2] ~ dbin(meanSocBees[i,2], Trials3[T2[i],2])
    NonSDiptera[T2[i],2] ~ dbin(meanNonSDip[i,2], Trials3[T2[i],2])
    NonVespWasps[T2[i],2] ~ dbin(meanNonVWasp[i,2], Trials3[T2[i],2])
    
    logit(meanSocBees[i,2]) <- C0[2] + C1[2]*meanM[T2[i]] + C2[2]*meanR[T2[i]] + C3[2]*meanA[T2[i]] + C4[2]*meanW[T2[i]] 
    logit(meanNonSDip[i,2]) <- E0[2] + E1[2]*meanM[T2[i]] + E2[2]*meanR[T2[i]] + E3[2]*meanA[T2[i]] + E4[2]*meanW[T2[i]] 
    logit(meanNonVWasp[i,2]) <- G0[2] + G1[2]*meanM[T2[i]] + G2[2]*meanR[T2[i]]+ G3[2]*meanA[T2[i]] + G4[2]*meanW[T2[i]] 
  }
  for(i in 1:length(T3)){
    SocialBees[T3[i],3] ~ dbin(meanSocBees[i,3], Trials3[T3[i],3])
    NonSDiptera[T3[i],3] ~ dbin(meanNonSDip[i,3], Trials3[T3[i],3])
    NonVespWasps[T3[i],3] ~ dbin(meanNonVWasp[i,3], Trials3[T3[i],3])
    
    logit(meanSocBees[i,3]) <- C0[3] + C1[3]*meanM[T3[i]] + C2[3]*meanR[T3[i]] + C3[3]*meanA[T3[i]] + C4[3]*meanW[T3[i]] 
    logit(meanNonSDip[i,3]) <- E0[3] + E1[3]*meanM[T3[i]] + E2[3]*meanR[T3[i]] + E3[3]*meanA[T3[i]] + E4[3]*meanW[T3[i]] 
    logit(meanNonVWasp[i,3]) <- G0[3] + G1[3]*meanM[T3[i]] + G2[3]*meanR[T3[i]]+ G3[3]*meanA[T3[i]] + G4[3]*meanW[T3[i]] 
  }
  
  for(i in 1:length(T4)){
    SocialBees[T4[i],4] ~ dbin(meanSocBees[i,4], Trials3[T4[i],4])
    NonSDiptera[T4[i],4] ~ dbin(meanNonSDip[i,4], Trials3[T4[i],4])
    NonVespWasps[T4[i],4] ~ dbin(meanNonVWasp[i,4], Trials3[T4[i],4])
    
    logit(meanSocBees[i,4]) <- C0[4] + C1[4]*meanM[T4[i]] + C2[4]*meanR[T4[i]] + C3[4]*meanA[T4[i]] + C4[4]*meanW[T4[i]] 
    logit(meanNonSDip[i,4]) <- E0[4] + E1[4]*meanM[T4[i]] + E2[4]*meanR[T4[i]] + E3[4]*meanA[T4[i]] + E4[4]*meanW[T4[i]] 
    logit(meanNonVWasp[i,4]) <- G0[4] + G1[4]*meanM[T4[i]] + G2[4]*meanR[T4[i]]+ G3[4]*meanA[T4[i]] + G4[4]*meanW[T4[i]] 
  }
  
  for(i in 1:length(T52)){
    SocialBees[T52[i],5] ~ dbin(meanSocBees[i,5], Trials3[T52[i],5])
    NonSDiptera[T52[i],5] ~ dbin(meanNonSDip[i,5], Trials3[T52[i],5])
    NonVespWasps[T52[i],5] ~ dbin(meanNonVWasp[i,5], Trials3[T52[i],5])
    
    logit(meanSocBees[i,5]) <- C0[5] + C1[5]*meanM[T52[i]] + C2[5]*meanR[T52[i]] + C3[5]*meanA[T52[i]] + C4[5]*meanW[T52[i]] 
    logit(meanNonSDip[i,5]) <- E0[5] + E1[5]*meanM[T52[i]] + E2[5]*meanR[T52[i]] + E3[5]*meanA[T52[i]] + E4[5]*meanW[T52[i]] 
    logit(meanNonVWasp[i,5]) <- G0[5] + G1[5]*meanM[T52[i]] + G2[5]*meanR[T52[i]]+ G3[5]*meanA[T52[i]] + G4[5]*meanW[T52[i]] 
  }
  
}


if(is.R()){
  filename <- file.path(tempdir(), "SHELOB.bug")}
write.model(SHELOB, filename)

inits <- list(list(IntmeanM=BR, IntmeanR=KR, meanA=B, meanW=B, B0=S, C0=S5, D0=S, E0=S5, F0=S, G0=S5, B1=S, B2=S, B3=S, B4=S, C1=S5, C2=S5, C3=S5, C4=S5, D1=S, D2=S, D3=S, D4=S, E1=S5, E2=S5, E3=S5, E4=S5, F1=S, F2=S, F3=S, F4=S, G1=S5, G2=S5, G3=S5, G4=S5), list(IntmeanM=IR, IntmeanR=OR, meanA=I, meanW=O, B0=E, C0=E5, D0=E, E0=E5, F0=E, G0=E5, B1=E, B2=E, B3=E, B4=E, C1=E5, C2=E5, C3=E5, C4=E5, D1=E, D2=E, D3=E, D4=E, E1=E5, E2=E5, E3=E5, E4=E5, F1=E, F2=E, F3=E, F4=E, G1=E5, G2=E5, G3=E5, G4=E5), list(IntmeanM=HR, IntmeanR=NR, meanA=H, meanW=N, B0=R, C0=R5, D0=R, E0=R5, F0=R, G0=R5, B1=R, B2=R, B3=R, B4=R, C1=R5, C2=R5, C3=R5, C4=R5, D1=R, D2=R, D3=R, D4=R, E1=R5, E2=R5, E3=R5, E4=R5, F1=R, F2=R, F3=R, F4=R, G1=R5, G2=R5, G3=R5, G4=R5)) 

data <- list("AntsNum", "TreatmentA", "WaspsNum", "TreatmentW", "nA", "nW", "MiceNum", "RatsNum", "TreatmentM", "TreatmentR", "SolitaryBees", "SocialBees", "Lepidop", "NonSDiptera", "Syrphid", "NonVespWasps", "Trials", "Trials3", "T1", "T2", "T3", "T4", "T5", "T52", "T6") 

parameters <- c("meanM", "meanR", "meanA", "meanW", "B0", "C0", "D0", "E0", "F0", "G0", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4", "E1", "E2", "E3", "E4", "F1", "F2", "F3", "F4", "G1", "G2", "G3", "G4") 

SHELOB <- jags(data=data, inits = inits, parameters.to.save=parameters, filename, n.burnin=475000, n.iter=500000, n.thin=1, n.chains=3) 
# Stop the clock
proc.time() - ptm

print(SHELOB)



SHELOB22.mcmc <- as.mcmc(SHELOB) #this saves all three runs together 