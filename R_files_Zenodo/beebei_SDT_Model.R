
#Signal detection theory model of parameter combinations that would result in observed discrimination thresholds
#based on within and among individual variation in calls of golden rocket frog calls (Anomaloglossus beebei)

#for paper:
##Tumulty, Fouilloux, Goyes Vallejos, & Bee. Predicting and measuring decision rules for social recognition in a Neotropical frog


##############################
### Load packages and data ###

#Load packages
library(ggplot2)
library(patchwork)

#Load data
data <- read.csv("beebei_call_properties.csv")

#Remove variables related to pulses 5 and 6 because of missing data
#Also remove pulse number variable
removenames <- c("Pulse3_Interval", "Pulse3_Period", "Pulse4_Duration", "Pulse4_Rise", "Pulse4_Fall",  
                 "Pulse4_HalfRise", "Pulse4_HalfFall", "Pulse4_Interval", "Pulse4_Period",
                 "Pulse5_Duration", "Pulse5_Rise", "Pulse5_Fall", "Pulse5_HalfRise", "Pulse5_HalfFall", 
                 "Pulse5_Interval", "Pulse5_Period", "Pulse6_Duration",
                 "Pulse6_Rise", "Pulse6_Fall", "Pulse6_HalfRise", "Pulse6_HalfFall",
                 "Pulse4_Max_Freq", "Pulse5_Max_Freq", "Pulse6_Max_Freq",
                 "Num_Pulses")

data <- data[,!(names(data) %in% removenames)]


##############################################################
### Compute temperature corrected values for each property ###

#New data frame to hold adjusted values
datatemp <- data[,1:4]
mean(data$Temp) #Mean = 23.5 degrees

#Run the loop
for (i in 5:ncol(data)){ # for each variable,
  fit <- lm(data[,i] ~ Temp, data=data, na.action=na.exclude) #Fit lm with temp
  intercept <- coefficients(fit)[[1]] #extract the intercept
  slope <- coefficients(fit)[[2]] #extract the slope
  pred <- intercept + slope*24 #predict value at temp of 24
  resids <- resid(fit, na.action=na.exclude) #Get model residuals, exclude NAs
  datatemp[,i] <- pred + resids #Adjust residuals by predicted value into new data frame
}

#Get the right column names
colnames(datatemp) <- colnames(data)


###############################################
### Within and Among Individual Differences ###

#Step 1) Compute pairwise percentage differences between each call and every other call in the data set
#Step 2) Label those percentage differences as "Within" individual or "Between" individual

# Function to compute Percentage Differences
percdiff <- function(x,y) {
  ifelse(x>y, (x-y)/y*100, (y-x)/x*100)
}

# Compute Pairwise Percentage Differences
out <- array(NA, dim=c(79800,ncol(datatemp)))
colnames(out) <- colnames(datatemp)
out <- as.data.frame(out)

IDdiffs <- outer(datatemp$ID, datatemp$ID, '==')
out$ID <- as.factor(IDdiffs[lower.tri(IDdiffs)])
out$ID <- as.factor(out$ID)

levels(out$ID)[levels(out$ID)=="TRUE"] <- "Within"
levels(out$ID)[levels(out$ID)=="FALSE"] <- "Between"

#Loop for all variables
for(i in 5:ncol(datatemp)){
  pdif <- outer(datatemp[,i], datatemp[,i], FUN=percdiff)
  out[,i] <- pdif[lower.tri(pdif)]
}

#Remove columns 2:5
out <- out[,-c(2:4)]


###########################################################
### Signal Detection Theory - Receiver Utility equation ###

#Utility equation
Utility <- function(pS) {
  pS*(pCD*d + (1-pCD)*m) + (1-pS)*(pFA*f + (1-pFA)*r)
}

#Save pS as range of values from 0 to 1
pS <- seq(0, 1, 0.0001)

#Thresholds, vector from 0 to 50
t <- seq(0, 50, 0.1)


########################################################################################################
############################################### Pulse Duration #########################################

#Pool data for different pulses
PD1 <- out[names(out) %in% c("ID", "Pulse1_Duration")]
PD1$Pulse <- "Pulse1"
colnames(PD1) <- c("ID", "Pulse_Duration", "Pulse")

PD2 <- out[names(out) %in% c("ID", "Pulse2_Duration")]
PD2$Pulse <- "Pulse2"
colnames(PD2) <- c("ID", "Pulse_Duration", "Pulse")

PD3 <- out[names(out) %in% c("ID", "Pulse3_Duration")]
PD3$Pulse <- "Pulse3"
colnames(PD3) <- c("ID", "Pulse_Duration", "Pulse")

PD <- rbind(PD1, PD2, PD3)

#Empty vectors to hold probabilities of false alarm (pFA) and correct detection (pCD)
pFA <- vector()
pCD <- vector()

#Compute pFA and pCD for all values of t
for(i in 1:length(t)){
  Nwithin <- length(na.omit(PD$Pulse_Duration[PD$ID=="Within"]))
  Nbetween <- length(na.omit(PD$Pulse_Duration[PD$ID=="Between"]))
  pFA[i] <- length(na.omit(PD$Pulse_Duration[PD$ID=="Within" & PD$Pulse_Duration>=t[i]]))/Nwithin
  pCD[i] <- length(na.omit(PD$Pulse_Duration[PD$ID=="Between" & PD$Pulse_Duration>=t[i]]))/Nbetween
}


### Loop for all possible values of pS (probability of stranger) and rel payoff (r-f)/(d-m)

# First modify numerator (r-f) from 1 to 20

#Utility equation by r - this one has pSi instead, so it works within loop
UtilityR <- function(r) {
  pSi*(pCD*d + (1-pCD)*m) + (1-pSi)*(pFA*f + (1-pFA)*r)
}

#f, d, and m stay the same, so that (d-m) = 1
f <- 1; d <- 2; m <- 1

#increase r to increase numerator
r <- seq(2, 21, 0.05)

#Matrix for results
t.opt.PD.r <- matrix(NA, length(r), length(pS))

#Compute optimal t for each value of pS and r
for(i in 1:length(pS)){
  pSi <- pS[i]
  for(j in 1:length(r)){
    t.opt.PD.r[j,i] <- t[which.max(UtilityR(r[j]))]
  }
}

colnames(t.opt.PD.r) <- pS
relpayoff.r <- (r-f)/(d-m)


# Then modify denominator (d-m) from 1 to 20

#Utility equation by d
UtilityD <- function(d) {
  pSi*(pCD*d + (1-pCD)*m) + (1-pSi)*(pFA*f + (1-pFA)*r)
}

#r, f, and m stay the same, so that (r-f) = 1
r <- 2; f <- 1; m <- 1

#test for range of d's (but stop at 2.05 so we can bind with the range of r's from above)
d <- seq(20, 2.05, by = -0.05)

#Matrix for results
t.opt.PD.d <- matrix(NA, length(d), length(pS))

#Compute optimal t for each value of pS and d
for(i in 1:length(pS)){
  pSi <- pS[i]
  for(j in 1:length(d)){
    t.opt.PD.d[j,i] <- t[which.max(UtilityD(d[j]))]
  }
}

colnames(t.opt.PD.d) <- pS
relpayoff.d <- (r-f)/(d-m)


#combine matrices modifying numerator and denominator
relpayoff <- c(relpayoff.d, relpayoff.r)
t.opt.PD <- rbind(t.opt.PD.d, t.opt.PD.r)
t.opt.PD <- data.frame(t.opt.PD)
t.opt.PD$relpayoff <- relpayoff


### Find conditions that result in t.opt of 9-12% ###

#Loop through all columns (all pS)
minrelpayoff <- vector()
maxrelpayoff <- vector()

for(i in 1:ncol(t.opt.PD)-1){
  minrelpayoff[i] <- min(t.opt.PD$relpayoff[t.opt.PD[,i]>=9 & t.opt.PD[,i]<=12])
  maxrelpayoff[i] <- max(t.opt.PD$relpayoff[t.opt.PD[,i]>=9 & t.opt.PD[,i]<=12])
}

PDconditions <- cbind(pS, minrelpayoff, maxrelpayoff)
PDconditions <- data.frame(PDconditions)
is.na(PDconditions) <- sapply(PDconditions, is.infinite)

p1 <- ggplot(PDconditions, aes(x=pS, y=minrelpayoff))+
  geom_line(aes(y=minrelpayoff))+
  geom_line(aes(y=maxrelpayoff))+
  geom_ribbon(aes(ymin=minrelpayoff,ymax=maxrelpayoff), fill="orange", alpha=0.5)+
  geom_hline(yintercept=1, lty=2)+
  ggtitle("Pulse Duration")+
  ylab("relative payoff ((r-f)/(d-m))")+
  xlab("stranger encounter rate (pS)")+
  scale_y_log10(breaks = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16), 
                labels = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16),
                limits=c(0.05, 20))+
  geom_point(x=0.5, y=log10(1), size=5)+ #apparently have to separately state log here for y 
  theme_bw()+
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
        text=element_text(size=12))

p1


########################################################################################################
############################################### Pulse Interval #########################################

#Pool data for different pulses
PI1 <- out[names(out) %in% c("ID", "Pulse1_Interval")]
PI1$Pulse <- "Pulse1"
colnames(PI1) <- c("ID", "Pulse_Interval", "Pulse")

PI2 <- out[names(out) %in% c("ID", "Pulse2_Interval")]
PI2$Pulse <- "Pulse2"
colnames(PI2) <- c("ID", "Pulse_Interval", "Pulse")

PI <- rbind(PI1, PI2)

#Empty vectors to hold probabilities of false alarm (pFA) and correct detection (pCD)
pFA <- vector()
pCD <- vector()

#Compute pFA and pCD for all values of t
for(i in 1:length(t)){
  Nwithin <- length(na.omit(PI$Pulse_Interval[PI$ID=="Within"]))
  Nbetween <- length(na.omit(PI$Pulse_Interval[PI$ID=="Between"]))
  pFA[i] <- length(na.omit(PI$Pulse_Interval[PI$ID=="Within" & PI$Pulse_Interval>t[i]]))/Nwithin
  pCD[i] <- length(na.omit(PI$Pulse_Interval[PI$ID=="Between" & PI$Pulse_Interval>t[i]]))/Nbetween
}


### Loop for all possible values of pS and rel payoff ###

#First modify numerator (r-f) from 1 to 20

#Utility equation by r - this one has pSi instead, so it works within loop
UtilityR <- function(r) {
  pSi*(pCD*d + (1-pCD)*m) + (1-pSi)*(pFA*f + (1-pFA)*r)
}

#f, d, and m stay the same, so that (d-m) = 1
f <- 1; d <- 2; m <- 1

#increase r to increase numerator
r <- seq(2, 21, 0.05)

#Matrix for results
t.opt.PI.r <- matrix(NA, length(r), length(pS))

#Compute optimal t for each value of pS and r
for(i in 1:length(pS)){
  pSi <- pS[i]
  for(j in 1:length(r)){
    t.opt.PI.r[j,i] <- t[which.max(UtilityR(r[j]))]
  }
}

colnames(t.opt.PI.r) <- pS


# Then modify denominator (d-m) from 1 to 20

#Utility equation by d
UtilityD <- function(d) {
  pSi*(pCD*d + (1-pCD)*m) + (1-pSi)*(pFA*f + (1-pFA)*r)
}

#r, f, and m stay the same, so that (r-f) = 1
r <- 2; f <- 1; m <- 1

#test for range of d's (but stop at 2.05 so we can bind with the range of r's from above)
d <- seq(20, 2.05, by = -0.05)

#Matrix for results
t.opt.PI.d <- matrix(NA, length(d), length(pS))

#Compute optimal t for each value of pS and d
for(i in 1:length(pS)){
  pSi <- pS[i]
  for(j in 1:length(d)){
    t.opt.PI.d[j,i] <- t[which.max(UtilityD(d[j]))]
  }
}

colnames(t.opt.PI.d) <- pS

#combine matrices modifying numerator and denominator
t.opt.PI <- rbind(t.opt.PI.d, t.opt.PI.r)
t.opt.PI <- data.frame(t.opt.PI)
t.opt.PI$relpayoff <- relpayoff


### Find conditions that result in t.opt of 9-12% ###

#Loop through all columns (all pS)
minrelpayoff <- vector()
maxrelpayoff <- vector()

for(i in 1:ncol(t.opt.PI)-1){
  minrelpayoff[i] <- min(t.opt.PI$relpayoff[t.opt.PI[,i]>=9 & t.opt.PI[,i]<=12])
  maxrelpayoff[i] <- max(t.opt.PI$relpayoff[t.opt.PI[,i]>=9 & t.opt.PI[,i]<=12])
}

PIconditions <- cbind(pS, minrelpayoff, maxrelpayoff)
PIconditions <- data.frame(PIconditions)
is.na(PIconditions) <- sapply(PIconditions, is.infinite)

p2 <- ggplot(PIconditions, aes(x=pS, y=minrelpayoff))+
  geom_line(aes(y=minrelpayoff))+
  geom_line(aes(y=maxrelpayoff))+
  geom_ribbon(aes(ymin=minrelpayoff,ymax=maxrelpayoff), fill="orange", alpha=0.5)+
  geom_hline(yintercept=1, lty=2)+
  ylab("relative payoff (r-f)/(d-m)")+
  xlab("stranger encounter rate (pS)")+
  ggtitle("Pulse Interval")+
  scale_y_log10(breaks = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16), 
                labels = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16),
                limits=c(0.05, 20))+
  geom_point(x=0.5, y=log10(1), size=5)+ #apparently have to separately state log here for y 
  theme_bw()+
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
        text=element_text(size=12))

p2


########################################################################################################
################################################ Pulse Period ##########################################

#Pool data for different pulses
PP1 <- out[names(out) %in% c("ID", "Pulse1_Period")]
PP1$Pulse <- "Pulse1"
colnames(PP1) <- c("ID", "Pulse_Period", "Pulse")

PP2 <- out[names(out) %in% c("ID", "Pulse2_Period")]
PP2$Pulse <- "Pulse2"
colnames(PP2) <- c("ID", "Pulse_Period", "Pulse")

PP <- rbind(PP1, PP2)

#Empty vectors to hold probabilities false alarm (pFA) and correct detection (pCD)
pFA <- vector()
pCD <- vector()

#Compute pFA and pCD for all values of t
for(i in 1:length(t)){
  Nwithin <- length(na.omit(PP$Pulse_Period[PP$ID=="Within"]))
  Nbetween <- length(na.omit(PP$Pulse_Period[PP$ID=="Between"]))
  pFA[i] <- length(na.omit(PP$Pulse_Period[PP$ID=="Within" & PP$Pulse_Period>t[i]]))/Nwithin
  pCD[i] <- length(na.omit(PP$Pulse_Period[PP$ID=="Between" & PP$Pulse_Period>t[i]]))/Nbetween
}


### Loop for all possible values of pS and rel payoff ###

#First modify numerator (r-f) from 1 to 20

#Utility equation by r - this one has pSi instead, so it works within loop
UtilityR <- function(r) {
  pSi*(pCD*d + (1-pCD)*m) + (1-pSi)*(pFA*f + (1-pFA)*r)
}

#f, d, and m stay the same, so that (d-m) = 1
f <- 1; d <- 2; m <- 1

#increase r to increase numerator
r <- seq(2, 21, 0.05)

#Matrix for results
t.opt.PP.r <- matrix(NA, length(r), length(pS))

#Compute optimal t for each value of pS and r
for(i in 1:length(pS)){
  pSi <- pS[i]
  for(j in 1:length(r)){
    t.opt.PP.r[j,i] <- t[which.max(UtilityR(r[j]))]
  }
}

colnames(t.opt.PP.r) <- pS


# Then modify denominator (d-m) from 1 to 20

#Utility equation by d
UtilityD <- function(d) {
  pSi*(pCD*d + (1-pCD)*m) + (1-pSi)*(pFA*f + (1-pFA)*r)
}

#r, f, and m stay the same, so that (r-f) = 1
r <- 2; f <- 1; m <- 1

#test for range of d's (but stop at 2.05 so we can bind with the range of r's from above)
d <- seq(20, 2.05, by = -0.05)

#Matrix for results
t.opt.PP.d <- matrix(NA, length(d), length(pS))

#Compute optimal t for each value of pS and d
for(i in 1:length(pS)){
  pSi <- pS[i]
  for(j in 1:length(d)){
    t.opt.PP.d[j,i] <- t[which.max(UtilityD(d[j]))]
  }
}

colnames(t.opt.PP.d) <- pS

#combine matrices modifying numerator and denominator
t.opt.PP <- rbind(t.opt.PP.d, t.opt.PP.r)
t.opt.PP <- data.frame(t.opt.PP)
t.opt.PP$relpayoff <- relpayoff


### Find conditions that result in t.opt of 9-12% ###

#Loop through all columns (all pS)
minrelpayoff <- vector()
maxrelpayoff <- vector()

for(i in 1:ncol(t.opt.PP)-1){
  minrelpayoff[i] <- min(t.opt.PP$relpayoff[t.opt.PP[,i]>=9 & t.opt.PP[,i]<=12])
  maxrelpayoff[i] <- max(t.opt.PP$relpayoff[t.opt.PP[,i]>=9 & t.opt.PP[,i]<=12])
}

PPconditions <- cbind(pS, minrelpayoff, maxrelpayoff)
PPconditions <- data.frame(PPconditions)
is.na(PPconditions) <- sapply(PPconditions, is.infinite)

p3 <- ggplot(PPconditions, aes(x=pS, y=minrelpayoff))+
  geom_line(aes(y=minrelpayoff))+
  geom_line(aes(y=maxrelpayoff))+
  geom_ribbon(aes(ymin=minrelpayoff,ymax=maxrelpayoff), fill="orange", alpha=0.5)+
  geom_hline(yintercept=1, lty=2)+
  #geom_vline(xintercept=0.5, lty=2)+
  ylab("relative payoff (r-f)/(d-m)")+
  xlab("stranger encounter rate (pS)")+
  ggtitle("Pulse Period")+
  scale_y_log10(breaks = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16), 
                labels = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16),
                limits=c(0.05, 20))+
  geom_point(x=0.5, y=log10(1), size=5)+ #apparently have to separately state log here for y 
  theme_bw()+
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),
        text=element_text(size=12))

p3


################
### Figure 5 ###

p1 + p2 + p3

