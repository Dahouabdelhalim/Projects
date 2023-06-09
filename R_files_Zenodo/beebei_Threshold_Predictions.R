
#Analysis of within and between individual differences in acoustic properties of golden rocket frog calls (Anomaloglossus beebei)

#for paper:
#Tumulty, Fouilloux, Goyes Vallejos, & Bee. Predicting and measuring decision rules for social recognition in a Neotropical frog


##############################
### Load packages and data ###

#Load packages
library(ggplot2)
library(DescTools)

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


##################################
### Predict Optimal Thresholds ###

#This approach maximizes the difference between the probability of Correct Detection (pCD) 
#and the probability of False Alarm (pFA)
#pCD = among individual differences greater than or equal to threshold
#pFA = within individual differences greater than or equal to threshold

#Note: takes a long time to run!
#Skip this and load in output file after loop to avoid running again

#Thresholds, vector from 0 to 1000
t <- seq(0, 1000, 0.1)

#Empty matrix to hold output
JMD <- matrix(NA, 6, length(2:ncol(out)))
colnames(JMD) <- colnames(out[,2:ncol(out)])
rownames(JMD) <- c("Discriminability", "t.opt", "pCR", "pFA", "pCD", "pMD")

#Empty vectors to hold values of pFA & pCD
pFA <- vector()
pCD <- vector()
D <- vector()

#Loop to compute t.opt for each property
for(j in 2:ncol(out)){
  for(i in 1:length(t)){
    Nwithin <- length(na.omit(out[,j][out$ID=="Within"]))
    Nbetween <- length(na.omit(out[,j][out$ID=="Between"]))
    pFA[i] <- length(na.omit(out[,j][out$ID=="Within" & out[,j]>=t[i]]))/Nwithin
    pCD[i] <- length(na.omit(out[,j][out$ID=="Between" & out[,j]>=t[i]]))/Nbetween
    D[i] <- pCD[i] - pFA[i] #difference between pCD and pFA for each treshold
    t.opt <- t[which.max(D)] #Max difference between pCD and pFA
    Discriminability <- AUC(x=pFA, y=pCD) #Area under the ROC curve for discriminability
    JMD["t.opt", j-1] <- t.opt
    JMD["Discriminability", j-1] <- Discriminability-0.5 #Subtract diagonal from area under curve
    JMD["pCR", j-1] <- (length(na.omit(out[,j][out$ID=="Within" & out[,j]<t.opt]))/Nwithin)*100
    JMD["pFA", j-1] <- (length(na.omit(out[,j][out$ID=="Within" & out[,j]>t.opt]))/Nwithin)*100
    JMD["pCD", j-1] <- (length(na.omit(out[,j][out$ID=="Between" & out[,j]>t.opt]))/Nbetween)*100
    JMD["pMD", j-1] <- (length(na.omit(out[,j][out$ID=="Between" & out[,j]<t.opt]))/Nbetween)*100
  }
}

JMD <- t(JMD)

#write.csv(JMD, file="beebei_JMDs.csv")

#Load back in saved JMD dataset to avoid running again
JMD <- read.csv("beebei_JMDs.csv")

#Table S-3 (or S-4 without temperature correction step)
JMD


#######################################################################
### Combine data for different pulses for pulse temporal properties ###

### Pulse Duration ###

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

#Thresholds, vector from 0 to 50 (can do smaller range for these properties)
t <- seq(0, 50, 0.1)

#Empty vectors to hold values of pFA & pCD
pFA <- vector()
pCD <- vector()
D <- vector()

#Run loop
for(i in 1:length(t)){
  Nwithin <- length(na.omit(PD$Pulse_Duration[PD$ID=="Within"]))
  Nbetween <- length(na.omit(PD$Pulse_Duration[PD$ID=="Between"]))
  pFA[i] <- length(na.omit(PD$Pulse_Duration[PD$ID=="Within" & PD$Pulse_Duration >= t[i]]))/Nwithin
  pCD[i] <- length(na.omit(PD$Pulse_Duration[PD$ID=="Between" & PD$Pulse_Duration >= t[i]]))/Nbetween
  D[i] <- pCD[i] - pFA[i]
  t.opt <- t[which.max(D)]
  Discriminability <- AUC(x=pFA, y=pCD)-0.5
  pCR.t <- (length(na.omit(PD$Pulse_Duration[PD$ID=="Within" & PD$Pulse_Duration < t.opt]))/Nwithin)*100
  pFA.t <- (length(na.omit(PD$Pulse_Duration[PD$ID=="Within" & PD$Pulse_Duration >= t.opt]))/Nwithin)*100
  pCD.t <- (length(na.omit(PD$Pulse_Duration[PD$ID=="Between" & PD$Pulse_Duration >= t.opt]))/Nbetween)*100
  pMD.t <- (length(na.omit(PD$Pulse_Duration[PD$ID=="Between" & PD$Pulse_Duration < t.opt]))/Nbetween)*100
}

#Wrangle results together into data frame
PDJMD <- data.frame(t.opt, Discriminability, pCR.t, pFA.t, pCD.t, pMD.t)

#Save pFA & pCD for ROC plots
PDROC <- cbind(t, pFA, pCD)


### Pulse Interval ###

PI1 <- out[names(out) %in% c("ID", "Pulse1_Interval")]
PI1$Pulse <- "Pulse1"
colnames(PI1) <- c("ID", "Pulse_Interval", "Pulse")

PI2 <- out[names(out) %in% c("ID", "Pulse2_Interval")]
PI2$Pulse <- "Pulse2"
colnames(PI2) <- c("ID", "Pulse_Interval", "Pulse")

PI <- rbind(PI1, PI2)

#Thresholds, vector from 0 to 50
t <- seq(0, 50, 0.1)

#Empty vectors to hold values of pFA & pCD
pFA <- vector()
pCD <- vector()
D <- vector()

#Run loop
for(i in 1:length(t)){
  Nwithin <- length(na.omit(PI[,2][PI$ID=="Within"]))
  Nbetween <- length(na.omit(PI[,2][PI$ID=="Between"]))
  pFA[i] <- length(na.omit(PI[,2][PI$ID=="Within" & PI[,2] >= t[i]]))/Nwithin
  pCD[i] <- length(na.omit(PI[,2][PI$ID=="Between" & PI[,2] >= t[i]]))/Nbetween
  D[i] <- pCD[i] - pFA[i]
  t.opt <- t[which.max(D)]
  Discriminability <- AUC(x=pFA, y=pCD)-0.5
  pCR.t <- (length(na.omit(PI[,2][PI$ID=="Within" & PI[,2] < t.opt]))/Nwithin)*100
  pFA.t <- (length(na.omit(PI[,2][PI$ID=="Within" & PI[,2] >= t.opt]))/Nwithin)*100
  pCD.t <- (length(na.omit(PI[,2][PI$ID=="Between" & PI[,2] >= t.opt]))/Nbetween)*100
  pMD.t <- (length(na.omit(PI[,2][PI$ID=="Between" & PI[,2] < t.opt]))/Nbetween)*100
}

#Wrangle results together into data frame
PIJMD <- data.frame(t.opt, Discriminability, pCR.t, pFA.t, pCD.t, pMD.t)

#Save pFA & pCD for ROC plot
PIROC <- cbind(t, pFA, pCD)


### Pulse Period ###

PP1 <- out[names(out) %in% c("ID", "Pulse1_Period")]
PP1$Pulse <- "Pulse1"
colnames(PP1) <- c("ID", "Pulse_Period", "Pulse")

PP2 <- out[names(out) %in% c("ID", "Pulse2_Period")]
PP2$Pulse <- "Pulse2"
colnames(PP2) <- c("ID", "Pulse_Period", "Pulse")

PP <- rbind(PP1, PP2)

#Thresholds, vector from 0 to 50
t <- seq(0, 50, 0.1)

#Empty vectors to hold values of pFA & pCD
pFA <- vector()
pCD <- vector()
D <- vector()

#Run loop
for(i in 1:length(t)){
  Nwithin <- length(na.omit(PP[,2][PP$ID=="Within"]))
  Nbetween <- length(na.omit(PP[,2][PP$ID=="Between"]))
  pFA[i] <- length(na.omit(PP[,2][PP$ID=="Within" & PP[,2] >= t[i]]))/Nwithin
  pCD[i] <- length(na.omit(PP[,2][PP$ID=="Between" & PP[,2] >= t[i]]))/Nbetween
  D[i] <- pCD[i] - pFA[i]
  t.opt <- t[which.max(D)]
  Discriminability <- AUC(x=pFA, y=pCD)-0.5
  pCR.t <- (length(na.omit(PP[,2][PP$ID=="Within" & PP[,2] < t.opt]))/Nwithin)*100
  pFA.t <- (length(na.omit(PP[,2][PP$ID=="Within" & PP[,2] >= t.opt]))/Nwithin)*100
  pCD.t <- (length(na.omit(PP[,2][PP$ID=="Between" & PP[,2] >= t.opt]))/Nbetween)*100
  pMD.t <- (length(na.omit(PP[,2][PP$ID=="Between" & PP[,2] < t.opt]))/Nbetween)*100
}

#Wrangle results together into data frame
PPJMD <- data.frame(t.opt, Discriminability, pCR.t, pFA.t, pCD.t, pMD.t)

#Save pFA & pCD for ROC plot
PPROC <- cbind(t, pFA, pCD)

#Merge properties together into data frame
PDJMD$Property <- "Pulse_Duration"
PIJMD$Property <- "Pulse_Interval"
PPJMD$Property <- "Pulse_Period"

PJMDs <- rbind(PDJMD, PIJMD, PPJMD)

#Table S-1 (or S-2 without temperature correction step)
PJMDs

#write.csv(PJMDs, file="beebei_mergedpulse_JMDs.csv")


#############
### Plots ###

### Fig. 2a. Within- and between-individual differences

#Pulse Duration
ggplot(PD, aes(Pulse_Duration, group=ID, fill=ID))+
  geom_density(alpha=0.5, adjust=3.5)+
  coord_cartesian(xlim = c(0,30), ylim=c(0,0.14))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0), breaks=seq(0,30, by=5))+
  scale_fill_manual(values=c("#E58601", "#46ACC8"), 
                    labels = c("Between-individual", "Within-individual"))+
  geom_vline(xintercept=PDJMD$t.opt, linetype=2, size=0.75)+
  ylab("Probatility")+
  xlab("% Difference in pulse period")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.title = element_blank())

#Pulse Interval
ggplot(PI, aes(Pulse_Interval, group=ID, fill=ID))+
  geom_density(alpha=0.5, adjust=3.5)+
  coord_cartesian(xlim = c(0,30), ylim=c(0,0.09))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0), breaks=seq(0,30, by=5))+
  scale_fill_manual(values=c("#E58601", "#46ACC8"), 
                    labels = c("Between-individual", "Within-individual"))+
  geom_vline(xintercept=PIJMD$t.opt, linetype=2, size=0.75)+
  ylab("Probatility")+
  xlab("% Difference in pulse period")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.title = element_blank())

#Pulse Period
ggplot(PP, aes(Pulse_Period, group=ID, fill=ID))+
  geom_density(alpha=0.5, adjust=3.5)+
  coord_cartesian(xlim = c(0,30), ylim=c(0,0.18))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0), breaks=seq(0,30, by=5))+
  scale_fill_manual(values=c("#E58601", "#46ACC8"), 
                    labels = c("Between-individual", "Within-individual"))+
  geom_vline(xintercept=PPJMD$t.opt, linetype=2, size=0.75)+
  ylab("Probatility")+
  xlab("% Difference in pulse period")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.title = element_blank())


### Fig. 2b. Receiver Operating Characteristics plots

#Pulse Duration
PDROC <- data.frame(PDROC)

#Also add point showing pFA and pCD corresponding to t.opt
which.max(PDROC$pCD-PDROC$pFA) #Row 67, corresponds to t=6.6

PDROC$pFA[67] #0.1743762
PDROC$pCD[67] #0.5999895

ggplot(PDROC, aes(x=pFA, y=pCD))+
  geom_path()+
  geom_abline(intercept=0, slope =1, linetype="dashed")+
  scale_x_continuous(limits=c(0,1), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,1), expand = c(0, 0))+
  geom_point(x=0.1743762, y=0.5999895, size=3)+
  ylab("Probability of correct detection (pD)")+
  xlab("Probability of false alarm (pF)")+
  theme_bw()+
  theme(panel.grid = element_blank())


#Pulse Interval
PIROC <- data.frame(PIROC)

#Also add point showing pFA and pCD corresponding to t.opt
which.max(PIROC$pCD-PIROC$pFA) #Row 93, corresponds to t=9.2

PIROC$pFA[93] #0.1888132
PIROC$pCD[93] #0.5682137

ggplot(PIROC, aes(x=pFA, y=pCD))+
  geom_path()+
  geom_abline(intercept=0, slope =1, linetype="dashed")+
  scale_x_continuous(limits=c(0,1), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,1), expand = c(0, 0))+
  geom_point(x=0.1888132, y=0.5682137, size=3)+
  ylab("Probability of correct detection (pD)")+
  xlab("Probability of false alarm (pF)")+
  theme_bw()+
  theme(panel.grid = element_blank())


#Pulse Period
PPROC <- data.frame(PPROC)

#Also add point showing pFA and pCD corresponding to t.opt
which.max(PPROC$pCD-PPROC$pFA) #Row 52, corresponds to t=5.1

PPROC$pFA[52] #0.1569338
PPROC$pCD[52] #0.6302965

ggplot(PPROC, aes(x=pFA, y=pCD))+
  geom_path()+
  geom_abline(intercept=0, slope =1, linetype="dashed")+
  scale_x_continuous(limits=c(0,1), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,1), expand = c(0, 0))+
  geom_point(x=0.1569338, y=0.6302965, size=3)+
  ylab("Probability of correct detection (pD)")+
  xlab("Probability of false alarm (pF)")+
  theme_bw()+
  theme(panel.grid = element_blank())

