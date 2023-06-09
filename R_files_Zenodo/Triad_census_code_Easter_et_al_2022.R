#Set working directory to where the data files are located on your computer
setwd()

#Load required packages
library(ggplot2)

#####################################################
#Load functions needed to extract normalized Z-scores
#####################################################

#This function first obtains the mean and standard deviation of counts of each triad across Uniform simulations with the specified density
#Note that only triads in which each node has at least one connection are included
#Triads are renumbered to match the ordering in: Waters JS, Fewell JH (2012) Information processing in social insect networks. PLoS ONE 7(7): e40337
#Z-scores are then obtained for each non-Uniform simulation by taking the count for a given triad, subtracting the corresponding mean count from the Uniform simulations, and dividing by the standard deviation (as described in the main text)

get_z_scores <- function(data, condition, density) {
  Nuni1 <- mean(data$countTriad4[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni1 <- sd(data$countTriad4[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni2 <- mean(data$countTriad5[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni2 <- sd(data$countTriad5[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni3 <- mean(data$countTriad6[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni3 <- sd(data$countTriad6[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni4 <- mean(data$countTriad7[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni4 <- sd(data$countTriad7[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni5 <- mean(data$countTriad8[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni5 <- sd(data$countTriad8[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni6 <- mean(data$countTriad11[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni6 <- sd(data$countTriad11[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni7 <- mean(data$countTriad9[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni7 <- sd(data$countTriad9[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni8 <- mean(data$countTriad10[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni8 <- sd(data$countTriad10[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni9 <- mean(data$countTriad12[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni9 <- sd(data$countTriad12[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni10 <- mean(data$countTriad13[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni10 <- sd(data$countTriad13[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni11 <- mean(data$countTriad14[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni11 <- sd(data$countTriad14[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni12 <- mean(data$countTriad15[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni12 <- sd(data$countTriad15[data$Condition == "Uniform" & data$networkDensity == density])
  Nuni13 <- mean(data$countTriad16[data$Condition == "Uniform" & data$networkDensity == density])
  sdUni13 <- sd(data$countTriad16[data$Condition == "Uniform" & data$networkDensity == density])
  
  z1 <- c()
  z2 <- c()
  z3 <- c()
  z4 <- c()
  z5 <- c()
  z6 <- c()
  z7 <- c()
  z8 <- c()
  z9 <- c()
  z10 <- c()
  z11 <- c()
  z12 <- c()
  z13 <- c()
  
  for(i in 1:max(data$simID)) {
    z1[i] <- (data$countTriad4[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni1)/sdUni1
    z2[i] <- (data$countTriad5[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni2)/sdUni2
    z3[i] <- (data$countTriad6[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni3)/sdUni3
    z4[i] <- (data$countTriad7[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni4)/sdUni4
    z5[i] <- (data$countTriad8[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni5)/sdUni5
    z6[i] <- (data$countTriad11[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni6)/sdUni6
    z7[i] <- (data$countTriad9[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni7)/sdUni7
    z8[i] <- (data$countTriad10[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni8)/sdUni8
    z9[i] <- (data$countTriad12[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni9)/sdUni9
    z10[i] <- (data$countTriad13[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni10)/sdUni10
    z11[i] <- (data$countTriad14[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni11)/sdUni11
    z12[i] <- (data$countTriad15[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni12)/sdUni12
    z13[i] <- (data$countTriad16[data$simID == i & data$Condition == condition & data$networkDensity == density] - Nuni13)/sdUni13
  }
  
  countsList <- list(z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13)
  return(countsList)
}

#Z-scores are normalized
normalize_z_scores <- function(data) {
  z1 <- c()
  z2 <- c()
  z3 <- c()
  z4 <- c()
  z5 <- c()
  z6 <- c()
  z7 <- c()
  z8 <- c()
  z9 <- c()
  z10 <- c()
  z11 <- c()
  z12 <- c()
  z13 <- c()
  
  for(i in 1:100) {
    z1[i] <- (data[[1]][i]) / ((sum(data[[1]] ^ 2)) ^ 0.5)
    z2[i] <- (data[[2]][i]) / ((sum(data[[2]] ^ 2)) ^ 0.5)
    z3[i] <- (data[[3]][i]) / ((sum(data[[3]] ^ 2)) ^ 0.5)
    z4[i] <- (data[[4]][i]) / ((sum(data[[4]] ^ 2)) ^ 0.5)
    z5[i] <- (data[[5]][i]) / ((sum(data[[5]] ^ 2)) ^ 0.5)
    z6[i] <- (data[[6]][i]) / ((sum(data[[6]] ^ 2)) ^ 0.5)
    z7[i] <- (data[[7]][i]) / ((sum(data[[7]] ^ 2)) ^ 0.5)
    z8[i] <- (data[[8]][i]) / ((sum(data[[8]] ^ 2)) ^ 0.5)
    z9[i] <- (data[[9]][i]) / ((sum(data[[9]] ^ 2)) ^ 0.5)
    z10[i] <- (data[[10]][i]) / ((sum(data[[10]] ^ 2)) ^ 0.5)
    z11[i] <- (data[[11]][i]) / ((sum(data[[11]] ^ 2)) ^ 0.5)
    z12[i] <- (data[[12]][i]) / ((sum(data[[12]] ^ 2)) ^ 0.5)
    z13[i] <- (data[[13]][i]) / ((sum(data[[13]] ^ 2)) ^ 0.5)
  }
  
  zList <- list(z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13)
  return(zList)
}

##Note that the following function is taken from:
#http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#####################
#To reproduce Fig. 4b
#####################

#Import data
triadData <- read.csv("Exp1_triadCensusData_networkDensity1000.csv", header = TRUE)

#Obtain normalized Z-scores for each set of non-Uniform simulations, using the Uniform sims as a null model
z_ActVar <- get_z_scores(data = triadData, condition = "ActVar", density = 1000)
z_ActVar_N <- normalize_z_scores(data = z_ActVar)
z_TurnVar <- get_z_scores(data = triadData, condition = "TurnVar", density = 1000)
z_TurnVar_N <- normalize_z_scores(data = z_TurnVar)
z_Uncorr <- get_z_scores(data = triadData, condition = "Uncorrelated", density = 1000)
z_Uncorr_N <- normalize_z_scores(data = z_Uncorr)
z_Corr <- get_z_scores(data = triadData, condition = "Correlated", density = 1000)
z_Corr_N <- normalize_z_scores(data = z_Corr)

#Combine normalized Z-scores into a dataframe
zScores <- c(unlist(z_ActVar_N), unlist(z_TurnVar_N), unlist(z_Uncorr_N), unlist(z_Corr_N))
triadID <- rep(c(rep(1,100), rep(2, 100), rep(3, 100), rep(4,100), rep(5, 100), rep(6, 100), rep(7,100), rep(8, 100), rep(9, 100), rep(10,100), rep(11, 100), rep(12, 100), rep(13, 100)), 4)
condition <- c(rep("ActVar", 1300), rep("TurnVar", 1300), rep("Uncorr", 1300), rep("Corr", 1300))
zData_N <- data.frame(Condition = condition, TriadID = triadID, Z = zScores)

#Only include data for fully connected triads (i.e., triad ID > 6)
zData_N <- zData_N[which(zData_N$TriadID > 6), ]

zData_N_Sum <- summarySE(data = zData_N, measurevar = "Z", groupvars = c("Condition", "TriadID"))

ggplot(data = zData_N_Sum, aes(y = Z, x = TriadID, color = as.factor(Condition))) + 
  geom_errorbar(aes(ymin = Z - se, ymax = Z + se), width = 0.1, size = 1.1, position = position_dodge(0.2)) +
  geom_line(size = 1.2, position = position_dodge(0.2)) +
  geom_point(show.legend = FALSE, size = 3.5, position = position_dodge(0.2)) +
  ylab("Normalized Z-score") + xlab("Triad ID")+
  scale_color_manual(name="Condition", labels=c("Activity: Variable", "Correlated", "TI: Variable", "Uncorrelated"), values=c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6")) +
  scale_x_continuous(limits = c(6,13.5), breaks = c(7, 8, 9, 10, 11, 12, 13)) + 
  theme(legend.position = c(0.8, 0.77)) + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 18), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) + 
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1)) +
  geom_hline(yintercept=0, size = 1.2)

#####################
#To reproduce Fig. S4
#####################

#Import the data
triadData <- read.csv("Exp2_triadCensusData.csv", header = TRUE)

#Restrict data according to whether interaction initiation and/or directionality were randomly assigned
triadData_RandomInit_RandomDirection <- triadData[which(triadData$Initiation == "Random" & triadData$Direction == "Random"), ]
triadData_RandomInit_FromActiveDirection <- triadData[which(triadData$Initiation == "Random" & triadData$Direction == "FromActive"), ]
triadData_ActivesInit_FromActiveDirection <- triadData[which(triadData$Initiation == "Actives" & triadData$Direction == "FromActive"), ]
triadData_ActivesInit_RandomDirection <- triadData[which(triadData$Initiation == "Actives" & triadData$Direction == "Random"), ]

#Obtain normalized Z-scores
#Within each initiation/directionality category, Act:Var sims were compared relative to Uniform sims
z_RI_RD <- get_z_scores(data = triadData_RandomInit_RandomDirection, condition = "ActVar", density = 1000)
z_RI_RD_N <- normalize_z_scores(data = z_RI_RD)
z_RI_AD <- get_z_scores(data = triadData_RandomInit_FromActiveDirection, condition = "ActVar", density = 1000)
z_RI_AD_N <- normalize_z_scores(data = z_RI_AD)
z_AI_AD <- get_z_scores(data = triadData_ActivesInit_FromActiveDirection, condition = "ActVar", density = 1000)
z_AI_AD_N <- normalize_z_scores(data = z_AI_AD)
z_AI_RD <- get_z_scores(data = triadData_ActivesInit_RandomDirection, condition = "ActVar", density = 1000)
z_AI_RD_N <- normalize_z_scores(data = z_AI_RD)

#Combine normalized Z-scores into a dataframe
zScores <- c(unlist(z_RI_RD_N), unlist(z_RI_AD_N), unlist(z_AI_AD_N), unlist(z_AI_RD_N))
triadID <- rep(c(rep(1,100), rep(2, 100), rep(3, 100), rep(4,100), rep(5, 100), rep(6, 100), rep(7,100), rep(8, 100), rep(9, 100), rep(10,100), rep(11, 100), rep(12, 100), rep(13, 100)), 4)
condition <- c(rep("RI_RD", 1300), rep("RI_AD", 1300), rep("AI_AD", 1300), rep("AI_RD", 1300))
zData_N <- data.frame(Condition = condition, TriadID = triadID, Z = zScores)

#Only include data for fully connected triads (i.e., triad ID > 6)
zData_N <- zData_N[which(zData_N$TriadID > 6), ]

zData_N_Sum <- summarySE(data = zData_N, measurevar = "Z", groupvars = c("Condition", "TriadID"))

ggplot(data = zData_N_Sum, aes(y = Z, x = TriadID, color = as.factor(Condition))) + 
  geom_errorbar(aes(ymin = Z - se, ymax = Z + se), width = 0.1, size = 1.1, position = position_dodge(0)) +
  geom_line(size = 1.2) +
  geom_point(show.legend = FALSE, size = 3.5, position = position_dodge(0)) +
  ylab("Normalized Z-score") + xlab("Triad ID")+
  scale_color_manual(name="Condition", labels=c("Init: Act., Dir.: Act.", "Init.: Act., Dir.: Ran.", "Init.: Ran., Dir: Act.", "Init.: Ran., Dir.: Ran."), values=c("#E69F00", "#009E73", "#0072B2", "#CC79A7")) +
  scale_x_continuous(limits = c(6,13.5), breaks = c(7, 8, 9, 10, 11, 12, 13)) + 
  theme(legend.position = c(0.85, 0.85)) + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 18), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) + 
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1)) +
  geom_hline(yintercept=0, size = 1.2)

#################
#Reproduce Fig. 7
#################

#Load data
triadData <- read.csv("NormalizedZScores_Fig7.csv", header = TRUE)

ggplot(data = triadData, aes(y = Z_N, x = TriadID, color = as.factor(Condition))) + 
  geom_errorbar(aes(ymin = Z_N - se, ymax = Z_N + se), width = 0.1, size = 1.4, position = position_dodge(0.4)) +
  geom_line(size = 1.5, position = position_dodge(0.4)) +
  geom_point(show.legend = FALSE, size = 3.5, position = position_dodge(0.4)) +
  ylab("Normalized Z-score") + xlab("Triad ID")+
  scale_color_manual(name="Null Model", labels=c("Activity: Variable", "Uniform"), values=c("#E69F00", "#56B4E9")) +
  scale_x_continuous(limits = c(0,13.5), breaks = c(1,2,3,4,5,6,7, 8, 9, 10, 11, 12, 13)) + 
  theme(legend.position = c(0.2, 0.2)) + 
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size = 18), legend.title = element_text(size = 18)) + 
  theme(axis.title.x = element_text(margin = margin(t = 12)), axis.title.y = element_text(margin = margin(r = 12))) + 
  theme(panel.background = element_rect(fill = "white"), legend.key = element_rect(fill = "white")) + 
  theme(axis.line.x.bottom = element_line(size = 1), axis.line.y.left = element_line(size = 1)) +
  geom_hline(yintercept=0, size = 1.2)


###############################################
