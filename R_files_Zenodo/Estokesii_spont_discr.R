#########################################################
#                     load packages                     #
#########################################################

# install package for data processing
library(devtools)
devtools::install_github("daniel1noble/cogdat")
library(cogdat)
# load packages
library(ggplot2)      # plotting
library(MCMCglmm)     # Bayesian modeling
library(smatr)        # Standardised major axis regression
library(simr)         # Power analysis for GLM by simulation

#########################################################
#               Repeatability function                  #
#########################################################

R_propA <- function(model){
  # variable extraction
  var_a  <- model$VCV[,1]
  var_e  <- model$VCV[,"units"]
  beta_0 <- model$Sol[,1]
  # inverse-logit transformation 
  P      <- exp(beta_0)/(1+exp(beta_0))
  
  # calculate repeatability
  Rpt <- (var_a*P*P/(1+exp(beta_0))^2) / ((var_a + var_e)*P*P/(1+exp(beta_0))^2 + P*(1-P))
  # estimates and CIs
  mea <- mean(Rpt) 
  med <- median(Rpt) 
  mod <- posterior.mode(Rpt) 
  CIs <- HPDinterval(Rpt)
  print(paste(c("mean ="), mea, c("median ="), med, c("mode ="),
              mod, c("CIlow ="), CIs[1], c("CIup ="), CIs[2]))
}

#########################################################
#                       load data                       #
#########################################################

raw  <- read.csv(file = "raw_spon_discr.csv", header = FALSE, colClasses = "character")

#########################################################
#                  data manipulation                     #
#########################################################

# transpose
proc      <- processDat(raw)

# remove empty column
proc$Date <- NULL
# turn variables into factors
proc      <- as.data.frame(unclass(proc))

# change variable types from factor to numeric
proc$SVL <- as.numeric(as.character(proc$SVL))
proc$weight <- as.numeric(as.character(proc$weight))


#########################################################
#                side preference test                   #
#########################################################

# remove the preference test from raw data set
proc_pref <- subset(proc, proc$stage == 'pref')
# check
str(proc_pref)

# remove data from ID2 (was removed from experiment)
proc_pref$choice[which(proc_pref$ID == 2)] <- NA
# remove NAs
proc_pref_noNA    <- proc_pref[complete.cases(proc_pref$choice), ]
proc_pref_noNA$ID <- factor(proc_pref_noNA$ID)
# check if correctly removed
str(proc_pref_noNA)

# translate left/right into numbers
proc_pref_noNA$choice_right[which(proc_pref_noNA$choice == 'r')] <- 1
proc_pref_noNA$choice_right[which(proc_pref_noNA$choice == 'l')] <- 0

# check if all lizards received 14 trials
aggregate(proc_pref_noNA$choice_right, list(proc_pref_noNA$ID), length)

# calculate the number of times each lizards chose the right
right_choice <- aggregate(proc_pref_noNA$choice_right, list(proc_pref_noNA$ID, proc_pref_noNA$experiment), sum)
colnames(right_choice) <- c("ID", "experiment", "Nright")

# does any of the lizards have a bias towards the right side?
# two tailed test because probability of success might either be significantly lower than 0.5 = left choice
# or significantly higher than 0.5 = right choice
right_choice$binom_p <- NA
for (i in 1:length(right_choice$Nright)){
  right_choice$binom_p[i] <- binom.test(right_choice$Nright[i], 14)$p.value
}

right_choice
# based on a binomial distribution, none of the lizards chose one side significantly more often


#########################################################
#                   body condition                      #
#              scaled mass index (SMI)                  #
#########################################################

# scaled mass index (SMI) for body condition (from: Peig, J., & Green, A. J. (2009). 
# New perspectives for estimating body condition from mass/length data: the scaled mass 
# index as an alternative method. Oikos, 118(12), 1883-1891)

# b_sma
b_sma <- sma(log2(weight) ~ log2(SVL), data = proc)$coef[[1]]$`coef(SMA)`[2]
# mean of SVL
L0 <- mean(proc$SVL)
# SMI for each individual
proc$SMI <- NA
proc$SMI <- proc$weight * (L0/proc$SVL)^b_sma

# check and save data frame
summary_SMI <- aggregate(proc$SMI, list(proc$ID, proc$experiment), mean)
# rename columns
colnames(summary_SMI) <- c("ID", "Experiment", "SMI")

# remove ID 2 (was removed from experiment)
summary_SMI[which(summary_SMI$ID == 2),] <- NA
summary_SMI <- summary_SMI[complete.cases(summary_SMI$ID), ]

# prepare data frame to be used in the wilcox test (data for independent groups in a separate column)
WC_SMI <- cbind(subset(summary_SMI, summary_SMI$Experiment == 'number'), subset(summary_SMI, summary_SMI$Experiment == 'size'))
colnames(WC_SMI) <- c("IDN", "ExpN", "SMI_N", "IDS", "ExpS", "SMI_S")

# compare SMI between experimental groups
wilcox.test(WC_SMI$SMI_N, WC_SMI$SMI_S, paired = FALSE, alternative = "two.sided", conf.int = TRUE)



#########################################################
#              Quantity discrimination                  #
#########################################################

##### further data manipulation to prepare the data set for analysis

# remove the preference test from raw data set
proc_no_pref <- subset(proc, proc$stage == 1)

# make trial, choice and session numeric
proc_no_pref$trial <- as.numeric(as.character(proc_no_pref$trial))
proc_no_pref$choice <- as.numeric(as.character(proc_no_pref$choice))
proc_no_pref$session <- as.numeric(as.character(proc_no_pref$session))

# change start time to a factor
proc_no_pref$start.time <- factor(proc_no_pref$start.time)
# turn time into a numeric variable with smallest number representing the earliest start time
# and the largest number representing the latest start time
proc_no_pref$time2 <- as.numeric(proc_no_pref$start.time)

# check variables
str(proc_no_pref)
summary(proc_no_pref)

# remove data from ID2 (was removed from experiment)
proc_no_pref$choice[which(proc_no_pref$ID == 2)] <- NA

# remove NAs
proc_no_pref_noNA <- proc_no_pref[complete.cases(proc_no_pref$choice), ]

# remove factor level for ID2
proc_no_pref_noNA$ID <- factor(proc_no_pref_noNA$ID)

# check if correctly removed
summary(proc_no_pref_noNA)
str(proc_no_pref_noNA)

# change variable types to factors
proc_no_pref_noNA$experiment <- factor(proc_no_pref_noNA$experiment)
proc_no_pref_noNA$stage <- factor(proc_no_pref_noNA$stage)
proc_no_pref_noNA$discr <- factor(proc_no_pref_noNA$discr)
proc_no_pref_noNA$researcher <- factor(proc_no_pref_noNA$researcher)

# change date to date format
proc_no_pref_noNA$date <- as.Date(proc_no_pref_noNA$date, "%d/%m/%y")


###### Probability of choosing the larger amount of food
#### comparison between combinations

# rename half of the combinations to just collate to 6 unique pairs
proc_no_pref_noNA$discr[which(proc_no_pref_noNA$discr == "4v1")] <- c("1v4")
proc_no_pref_noNA$discr[which(proc_no_pref_noNA$discr == "3v1")] <- c("1v3")
proc_no_pref_noNA$discr[which(proc_no_pref_noNA$discr == "2v1")] <- c("1v2")
proc_no_pref_noNA$discr[which(proc_no_pref_noNA$discr == "4v2")] <- c("2v4")
proc_no_pref_noNA$discr[which(proc_no_pref_noNA$discr == "3v2")] <- c("2v3")
proc_no_pref_noNA$discr[which(proc_no_pref_noNA$discr == "4v3")] <- c("3v4")

# record if the larger or smaller amount was chosen in each trial
# 1 = larger, 0 = smaller
proc_no_pref_noNA$larger <- NA
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "1v4" 
                               & proc_no_pref_noNA$choice == 4)] <- 1
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "1v4" 
                               & proc_no_pref_noNA$choice == 1)] <- 0
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "1v3" 
                               & proc_no_pref_noNA$choice == 3)] <- 1
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "1v3" 
                               & proc_no_pref_noNA$choice == 1)] <- 0
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "1v2" 
                               & proc_no_pref_noNA$choice == 2)] <- 1
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "1v2" 
                               & proc_no_pref_noNA$choice == 1)] <- 0
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "2v4" 
                               & proc_no_pref_noNA$choice == 4)] <- 1
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "2v4" 
                               & proc_no_pref_noNA$choice == 2)] <- 0
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "2v3" 
                               & proc_no_pref_noNA$choice == 3)] <- 1
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "2v3" 
                               & proc_no_pref_noNA$choice == 2)] <- 0
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "3v4" 
                               & proc_no_pref_noNA$choice == 4)] <- 1
proc_no_pref_noNA$larger[which(proc_no_pref_noNA$discr == "3v4" 
                               & proc_no_pref_noNA$choice == 3)] <- 0

# remove trials testing 0v1 and 1v0
proc_no_pref_noNA$choice[which(proc_no_pref_noNA$discr == "1v0" |
                               proc_no_pref_noNA$discr == "0v1")] <- NA
proc_no_pref_noNA <- proc_no_pref_noNA[complete.cases(proc_no_pref_noNA$choice), ]

# number of trials for each discrimination and each lizard
N_trials <- aggregate(proc_no_pref_noNA$larger, list(proc_no_pref_noNA$discr, 
                                                     proc_no_pref_noNA$ID, 
                                                     proc_no_pref_noNA$experiment), length)
# number of times larger amount was chosen
N_larger <- aggregate(proc_no_pref_noNA$larger, list(proc_no_pref_noNA$discr,
                                                     proc_no_pref_noNA$ID,
                                                     proc_no_pref_noNA$experiment), sum)
# rename columns
colnames(N_larger) <- c("discr", "ID", "experiment", "larger")
# calculate the number of trial the smaller amount was chosen
N_larger$N_trials <- N_trials$x
N_larger$smaller  <- N_larger$N_trials - N_larger$larger
N_larger$prop <- N_larger$larger/N_larger$N_trials


##################### main analysis #####################

# standard prior for binomial data (0 and 1)
prior <- list(R = list(V = 1, fix = 1), 
              G = list(G1 = list(V = diag(2), nu = 0.002)))

# iterations, burn-in and thinning interval for MCMCglmm
itts  <- 6000000
burn  <- 10000
thins <- 130

################# number test only #################
# subset to data from the number test only
proc_no_pref_noNA_number <- subset(proc_no_pref_noNA, proc_no_pref_noNA$experiment == "number")

# order contrasts
proc_no_pref_noNA_number$discr <- factor(proc_no_pref_noNA_number$discr, levels = c("1v4", "1v3", "2v4", "1v2", "2v3", "3v4"))

model_number_comb_comp <- MCMCglmm(larger ~ discr, random = ~us(session+1):ID,
                                   family = "categorical", data = proc_no_pref_noNA_number, 
                                   prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                   verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_number_comb_comp$VCV)
autocorr <- as.data.frame(autocorr(model_number_comb_comp$Sol))
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_number_comb_comp$Sol)
geweke.plot(model_number_comb_comp$VCV)
# diagnostic for chain length
heidel.diag(model_number_comb_comp$Sol)
# diagnostic for chain length and sampling
plot(model_number_comb_comp$Sol)

# model results
summary(model_number_comb_comp)

################# bold boxes of Table 1 (number test)
############# Was the larger quantity chosen significantly above chance?
number_estimates_intercept <- as.data.frame(1:6, c("1v4", "1v3", "2v4", "1v2", "2v3", "3v4"))
colnames(number_estimates_intercept) <- c('Intercepts')
# Intercept (logit scale) for each contrast
# 1v4 (first contrast)
number_estimates_intercept[1,1] <- mean(model_number_comb_comp$Sol[,1])
# 1v3
number_estimates_intercept[2,1] <- mean(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])
# 2v4
number_estimates_intercept[3,1] <- mean(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])
# 1v2
number_estimates_intercept[4,1] <- mean(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])
# 2v3
number_estimates_intercept[5,1] <- mean(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])
# 3v4
number_estimates_intercept[6,1] <- mean(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])
# lower CIs for each contrast
number_estimates_intercept$lowerCI <- NA
number_estimates_intercept[1,2] <- HPDinterval(model_number_comb_comp$Sol[,1])[1]
number_estimates_intercept[2,2] <- HPDinterval(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])[1]
number_estimates_intercept[3,2] <- HPDinterval(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])[1]
number_estimates_intercept[4,2] <- HPDinterval(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])[1]
number_estimates_intercept[5,2] <- HPDinterval(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])[1]
number_estimates_intercept[6,2] <- HPDinterval(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])[1]
# upper CIs for each contrast
number_estimates_intercept$upperCI <- NA
number_estimates_intercept[1,3] <- HPDinterval(model_number_comb_comp$Sol[,1])[2]
number_estimates_intercept[2,3] <- HPDinterval(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])[2]
number_estimates_intercept[3,3] <- HPDinterval(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])[2]
number_estimates_intercept[4,3] <- HPDinterval(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])[2]
number_estimates_intercept[5,3] <- HPDinterval(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])[2]
number_estimates_intercept[6,3] <- HPDinterval(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])[2]

################# bottom-left half of Table 1 (number test)
############# Was 1v4 performed better than the other contrasts?
number_estimates_1v4 <- as.data.frame(1:5, c("1v3", "2v4", "1v2", "2v3", "3v4"))
colnames(number_estimates_1v4) <- c("1v4")
# 1v4 against 1v3
number_estimates_1v4[1,1] <- mean(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))
# 1v4 against 2v4
number_estimates_1v4[2,1] <- mean(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))
# 1v4 against 1v2
number_estimates_1v4[3,1] <- mean(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))
# 1v4 against 2v3
number_estimates_1v4[4,1] <- mean(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))
# 1v4 against 3v4
number_estimates_1v4[5,1] <- mean(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))
# lower CIs for each contrast
number_estimates_1v4$lowerCI <- NA
number_estimates_1v4[1,2] <- HPDinterval(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))[1]
number_estimates_1v4[2,2] <- HPDinterval(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))[1]
number_estimates_1v4[3,2] <- HPDinterval(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))[1]
number_estimates_1v4[4,2] <- HPDinterval(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))[1]
number_estimates_1v4[5,2] <- HPDinterval(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))[1]
# upper CIs for each contrast
number_estimates_1v4$upperCI <- NA
number_estimates_1v4[1,3] <- HPDinterval(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))[2]
number_estimates_1v4[2,3] <- HPDinterval(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))[2]
number_estimates_1v4[3,3] <- HPDinterval(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))[2]
number_estimates_1v4[4,3] <- HPDinterval(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))[2]
number_estimates_1v4[5,3] <- HPDinterval(model_number_comb_comp$Sol[,1]-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))[2]

############# Was 1v3 performed better than the other contrasts?
number_estimates_1v3 <- as.data.frame(1:4, c("2v4", "1v2", "2v3", "3v4"))
colnames(number_estimates_1v3) <- c("1v3")
# 1v3 against 2v4
number_estimates_1v3[1,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))
# 1v3 against 1v2
number_estimates_1v3[2,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))
# 1v3 against 2v3
number_estimates_1v3[3,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))
# 1v3 against 3v4
number_estimates_1v3[4,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))
# lower CIs for each contrast
number_estimates_1v3$lowerCI <- NA
number_estimates_1v3[1,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))[1]
number_estimates_1v3[2,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))[1]
number_estimates_1v3[3,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))[1]
number_estimates_1v3[4,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))[1]
# upper CIs for each contrast
number_estimates_1v3$upperCI <- NA
number_estimates_1v3[1,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))[2]
number_estimates_1v3[2,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))[2]
number_estimates_1v3[3,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))[2]
number_estimates_1v3[4,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))[2]

############# Was 2v4 performed better than the other contrasts?
number_estimates_2v4 <- as.data.frame(1:3, c("1v2", "2v3", "3v4"))
colnames(number_estimates_2v4) <- c("2v4")
# 2v4 against 1v2
number_estimates_2v4[1,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))
# 2v4 against 2v3
number_estimates_2v4[2,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))
# 2v4 against 3v4
number_estimates_2v4[3,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))
# lower CIs for each contrast
number_estimates_2v4$lowerCI <- NA
number_estimates_2v4[1,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))[1]
number_estimates_2v4[2,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))[1]
number_estimates_2v4[3,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))[1]
# upper CIs for each contrast
number_estimates_2v4$upperCI <- NA
number_estimates_2v4[1,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))[2]
number_estimates_2v4[2,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))[2]
number_estimates_2v4[3,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))[2]

############# Was 1v2 performed better than the other contrasts?
number_estimates_1v2 <- as.data.frame(1:2, c("2v3", "3v4"))
colnames(number_estimates_1v2) <- c("1v2")
# 1v2 against 2v3
number_estimates_1v2[1,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))
# 1v2 against 3v4
number_estimates_1v2[2,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))
# lower CIs for each contrast
number_estimates_1v2$lowerCI <- NA
number_estimates_1v2[1,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))[1]
number_estimates_1v2[2,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))[1]
# upper CIs for each contrast
number_estimates_1v2$upperCI <- NA
number_estimates_1v2[1,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))[2]
number_estimates_1v2[2,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))[2]

############# Was 2v3 performed better than the other contrasts?
number_estimates_2v3 <- as.data.frame(1, c("3v4"))
colnames(number_estimates_2v3) <- c("2v3")
# 2v3 against 3v4
number_estimates_2v3[1,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))
# lower CIs for each contrast
number_estimates_2v3$lowerCI <- NA
number_estimates_2v3[1,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))[1]
# upper CIs for each contrast
number_estimates_2v3$upperCI <- NA
number_estimates_2v3[1,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6]))[2]

################# top-right half of Table 1 (number test) (mirrored bottom-left)
############# Was any other contrast performed better than 1v4?
number_estimates_1v42 <- as.data.frame(1:5, c("1v3", "2v4", "1v2", "2v3", "3v4"))
colnames(number_estimates_1v42) <- c("1v4")
# 1v4 against 1v3
number_estimates_1v42[1,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-model_number_comb_comp$Sol[,1])
# 1v4 against 2v4
number_estimates_1v42[2,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-model_number_comb_comp$Sol[,1])
# 1v4 against 1v2
number_estimates_1v42[3,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-model_number_comb_comp$Sol[,1])
# 1v4 against 2v3
number_estimates_1v42[4,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-model_number_comb_comp$Sol[,1])
# 1v4 against 3v4
number_estimates_1v42[5,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-model_number_comb_comp$Sol[,1])
# lower CIs for each contrast
number_estimates_1v42$lowerCI <- NA
number_estimates_1v42[1,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-model_number_comb_comp$Sol[,1])[1]
number_estimates_1v42[2,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-model_number_comb_comp$Sol[,1])[1]
number_estimates_1v42[3,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-model_number_comb_comp$Sol[,1])[1]
number_estimates_1v42[4,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-model_number_comb_comp$Sol[,1])[1]
number_estimates_1v42[5,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-model_number_comb_comp$Sol[,1])[1]
# upper CIs for each contrast
number_estimates_1v42$upperCI <- NA
number_estimates_1v42[1,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2])-model_number_comb_comp$Sol[,1])[2]
number_estimates_1v42[2,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-model_number_comb_comp$Sol[,1])[2]
number_estimates_1v42[3,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-model_number_comb_comp$Sol[,1])[2]
number_estimates_1v42[4,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-model_number_comb_comp$Sol[,1])[2]
number_estimates_1v42[5,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-model_number_comb_comp$Sol[,1])[2]

############# Was any other contrast performed better than 1v3?
number_estimates_1v32 <- as.data.frame(1:4, c("2v4", "1v2", "2v3", "3v4"))
colnames(number_estimates_1v32) <- c("1v3")
# 1v3 against 2v4
number_estimates_1v32[1,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))
# 1v3 against 1v2
number_estimates_1v32[2,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))
# 1v3 against 2v3
number_estimates_1v32[3,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))
# 1v3 against 3v4
number_estimates_1v32[4,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))
# lower CIs for each contrast
number_estimates_1v32$lowerCI <- NA
number_estimates_1v32[1,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))[1]
number_estimates_1v32[2,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))[1]
number_estimates_1v32[3,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))[1]
number_estimates_1v32[4,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))[1]
# upper CIs for each contrast
number_estimates_1v32$upperCI <- NA
number_estimates_1v32[1,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))[2]
number_estimates_1v32[2,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))[2]
number_estimates_1v32[3,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))[2]
number_estimates_1v32[4,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,2]))[2]

############# Was any other contrast performed better than 2v4?
number_estimates_2v42 <- as.data.frame(1:3, c("1v2", "2v3", "3v4"))
colnames(number_estimates_2v42) <- c("2v4")
# 2v4 against 1v2
number_estimates_2v42[1,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))
# 2v4 against 2v3
number_estimates_2v42[2,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))
# 2v4 against 3v4
number_estimates_2v42[3,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))
# lower CIs for each contrast
number_estimates_2v42$lowerCI <- NA
number_estimates_2v42[1,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))[1]
number_estimates_2v42[2,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))[1]
number_estimates_2v42[3,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))[1]
# upper CIs for each contrast
number_estimates_2v42$upperCI <- NA
number_estimates_2v42[1,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))[2]
number_estimates_2v42[2,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))[2]
number_estimates_2v42[3,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,3]))[2]

############# Was any other contrast performed better than 1v2?
number_estimates_1v22 <- as.data.frame(1:2, c("2v3", "3v4"))
colnames(number_estimates_1v22) <- c("1v2")
# 1v2 against 2v3
number_estimates_1v22[1,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))
# 1v2 against 3v4
number_estimates_1v22[2,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))
# lower CIs for each contrast
number_estimates_1v22$lowerCI <- NA
number_estimates_1v22[1,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))[1]
number_estimates_1v22[2,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))[1]
# upper CIs for each contrast
number_estimates_1v22$upperCI <- NA
number_estimates_1v22[1,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))[2]
number_estimates_1v22[2,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,4]))[2]

############# Was 2v3 performed better than the other contrasts?
number_estimates_2v32 <- as.data.frame(1, c("3v4"))
colnames(number_estimates_2v32) <- c("2v3")
# 2v3 against 3v4
number_estimates_2v32[1,1] <- mean((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))
# lower CIs for each contrast
number_estimates_2v32$lowerCI <- NA
number_estimates_2v32[1,2] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))[1]
# upper CIs for each contrast
number_estimates_2v32$upperCI <- NA
number_estimates_2v32[1,3] <- HPDinterval((model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,6])-(model_number_comb_comp$Sol[,1]+model_number_comb_comp$Sol[,5]))[2]



########## time of day ############

model_number_corr_time <- MCMCglmm(larger ~ time2, random = ~us(session+1):ID,
                                   family = "categorical", data = proc_no_pref_noNA_number, 
                                   prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                   verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_number_corr_time$VCV)
autocorr <- as.data.frame(autocorr(model_number_corr_time$Sol))
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_number_corr_time$Sol)
geweke.plot(model_number_corr_time$VCV)
# diagnostic for chain length
heidel.diag(model_number_corr_time$Sol)
# diagnostic for chain length and sampling
plot(model_number_corr_time$Sol)

# model results
summary(model_number_corr_time)





################# size test only #################

# subset to data from the size test only
proc_no_pref_noNA_size <- subset(proc_no_pref_noNA, proc_no_pref_noNA$experiment == "size")

# order contrasts
proc_no_pref_noNA_size$discr <- factor(proc_no_pref_noNA_size$discr, levels = c("1v4", "1v3", "2v4", "1v2", "2v3", "3v4"))

model_size_comb_comp <- MCMCglmm(larger ~ discr, random = ~us(session+1):ID,
                                 family = "categorical", data = proc_no_pref_noNA_size, 
                                 prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                 verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_size_comb_comp$VCV)
autocorr <- as.data.frame(autocorr(model_size_comb_comp$Sol))
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_size_comb_comp$Sol)
geweke.plot(model_size_comb_comp$VCV)
# diagnostic for chain length
heidel.diag(model_size_comb_comp$Sol)
# diagnostic for chain length and sampling
plot(model_size_comb_comp$Sol)

# model results
summary(model_size_comb_comp)


################# bold boxes of Table 1 (size test)
############# Was the larger quantity chosen significantly above chance?
size_estimates_intercept <- as.data.frame(1:6, c("1v4", "1v3", "2v4", "1v2", "2v3", "3v4"))
colnames(size_estimates_intercept) <- c('Intercepts')
# Intercept (logit scale) for each contrast
# 1v4 (first contrast)
size_estimates_intercept[1,1] <- mean(model_size_comb_comp$Sol[,1])
# 1v3
size_estimates_intercept[2,1] <- mean(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])
# 2v4
size_estimates_intercept[3,1] <- mean(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])
# 1v2
size_estimates_intercept[4,1] <- mean(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])
# 2v3
size_estimates_intercept[5,1] <- mean(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])
# 3v4
size_estimates_intercept[6,1] <- mean(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])
# lower CIs for each contrast
size_estimates_intercept$lowerCI <- NA
size_estimates_intercept[1,2] <- HPDinterval(model_size_comb_comp$Sol[,1])[1]
size_estimates_intercept[2,2] <- HPDinterval(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])[1]
size_estimates_intercept[3,2] <- HPDinterval(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])[1]
size_estimates_intercept[4,2] <- HPDinterval(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])[1]
size_estimates_intercept[5,2] <- HPDinterval(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])[1]
size_estimates_intercept[6,2] <- HPDinterval(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])[1]
# upper CIs for each contrast
size_estimates_intercept$upperCI <- NA
size_estimates_intercept[1,3] <- HPDinterval(model_size_comb_comp$Sol[,1])[2]
size_estimates_intercept[2,3] <- HPDinterval(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])[2]
size_estimates_intercept[3,3] <- HPDinterval(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])[2]
size_estimates_intercept[4,3] <- HPDinterval(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])[2]
size_estimates_intercept[5,3] <- HPDinterval(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])[2]
size_estimates_intercept[6,3] <- HPDinterval(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])[2]

################# bottom-left half of Table 1 (size test)
############# Was 1v4 performed better than the other contrasts?
size_estimates_1v4 <- as.data.frame(1:5, c("1v3", "2v4", "1v2", "2v3", "3v4"))
colnames(size_estimates_1v4) <- c("1v4")
# 1v3
size_estimates_1v4[1,1] <- mean(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))
# 2v4
size_estimates_1v4[2,1] <- mean(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))
# 1v2
size_estimates_1v4[3,1] <- mean(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))
# 2v3
size_estimates_1v4[4,1] <- mean(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))
# 3v4
size_estimates_1v4[5,1] <- mean(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))
# lower CIs for each contrast
size_estimates_1v4$lowerCI <- NA
size_estimates_1v4[1,2] <- HPDinterval(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))[1]
size_estimates_1v4[2,2] <- HPDinterval(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))[1]
size_estimates_1v4[3,2] <- HPDinterval(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))[1]
size_estimates_1v4[4,2] <- HPDinterval(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))[1]
size_estimates_1v4[5,2] <- HPDinterval(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))[1]
# upper CIs for each contrast
size_estimates_1v4$upperCI <- NA
size_estimates_1v4[1,3] <- HPDinterval(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))[2]
size_estimates_1v4[2,3] <- HPDinterval(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))[2]
size_estimates_1v4[3,3] <- HPDinterval(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))[2]
size_estimates_1v4[4,3] <- HPDinterval(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))[2]
size_estimates_1v4[5,3] <- HPDinterval(model_size_comb_comp$Sol[,1]-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))[2]

############# Was 1v3 performed better than the other contrasts?
size_estimates_1v3 <- as.data.frame(1:4, c("2v4", "1v2", "2v3", "3v4"))
colnames(size_estimates_1v3) <- c("1v3")
# 2v4
size_estimates_1v3[1,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))
# 1v2
size_estimates_1v3[2,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))
# 2v3
size_estimates_1v3[3,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))
# 3v4
size_estimates_1v3[4,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))
# lower CIs for each contrast
size_estimates_1v3$lowerCI <- NA
size_estimates_1v3[1,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))[1]
size_estimates_1v3[2,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))[1]
size_estimates_1v3[3,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))[1]
size_estimates_1v3[4,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))[1]
# upper CIs for each contrast
size_estimates_1v3$upperCI <- NA
size_estimates_1v3[1,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))[2]
size_estimates_1v3[2,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))[2]
size_estimates_1v3[3,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))[2]
size_estimates_1v3[4,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))[2]

############# Was 2v4 performed better than the other contrasts?
size_estimates_2v4 <- as.data.frame(1:3, c("1v2", "2v3", "3v4"))
colnames(size_estimates_2v4) <- c("2v4")
# 1v2
size_estimates_2v4[1,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))
# 2v3
size_estimates_2v4[2,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))
#3v4
size_estimates_2v4[3,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))
# lower CIs for each contrast
size_estimates_2v4$lowerCI <- NA
size_estimates_2v4[1,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))[1]
size_estimates_2v4[2,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))[1]
size_estimates_2v4[3,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))[1]
# upper CIs for each contrast
size_estimates_2v4$upperCI <- NA
size_estimates_2v4[1,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))[2]
size_estimates_2v4[2,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))[2]
size_estimates_2v4[3,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))[2]

############# Was 1v2 performed better than the other contrasts?
size_estimates_1v2 <- as.data.frame(1:2, c("2v3", "3v4"))
colnames(size_estimates_1v2) <- c("1v2")
# 2v3
size_estimates_1v2[1,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))
# 3v4
size_estimates_1v2[2,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))
# lower CIs for each contrast
size_estimates_1v2$lowerCI <- NA
size_estimates_1v2[1,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))[1]
size_estimates_1v2[2,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))[1]
# upper CIs for each contrast
size_estimates_1v2$upperCI <- NA
size_estimates_1v2[1,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))[2]
size_estimates_1v2[2,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))[2]

############# Was 2v3 performed better than the other contrasts?
size_estimates_2v3 <- as.data.frame(1, c("3v4"))
colnames(size_estimates_2v3) <- c("2v3")
# 3v4
size_estimates_2v3[1,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))
# lower CIs for each contrast
size_estimates_2v3$lowerCI <- NA
size_estimates_2v3[1,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))[1]
# upper CIs for each contrast
size_estimates_2v3$upperCI <- NA
size_estimates_2v3[1,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6]))[2]

################# top-right half of Table 1 (size test) (mirrored bottom-left)
############# Was any other contrast performed better than 1v4?
size_estimates_1v42 <- as.data.frame(1:5, c("1v3", "2v4", "1v2", "2v3", "3v4"))
colnames(size_estimates_1v42) <- c("1v4")
# 1v3
size_estimates_1v42[1,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-model_size_comb_comp$Sol[,1])
# 2v4
size_estimates_1v42[2,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-model_size_comb_comp$Sol[,1])
# 1v2
size_estimates_1v42[3,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-model_size_comb_comp$Sol[,1])
# 2v3
size_estimates_1v42[4,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-model_size_comb_comp$Sol[,1])
# 3v4
size_estimates_1v42[5,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-model_size_comb_comp$Sol[,1])
# lower CIs for each contrast
size_estimates_1v42$lowerCI <- NA
size_estimates_1v42[1,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-model_size_comb_comp$Sol[,1])[1]
size_estimates_1v42[2,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-model_size_comb_comp$Sol[,1])[1]
size_estimates_1v42[3,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-model_size_comb_comp$Sol[,1])[1]
size_estimates_1v42[4,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-model_size_comb_comp$Sol[,1])[1]
size_estimates_1v42[5,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-model_size_comb_comp$Sol[,1])[1]
# upper CIs for each contrast
size_estimates_1v42$upperCI <- NA
size_estimates_1v42[1,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2])-model_size_comb_comp$Sol[,1])[2]
size_estimates_1v42[2,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-model_size_comb_comp$Sol[,1])[2]
size_estimates_1v42[3,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-model_size_comb_comp$Sol[,1])[2]
size_estimates_1v42[4,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-model_size_comb_comp$Sol[,1])[2]
size_estimates_1v42[5,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-model_size_comb_comp$Sol[,1])[2]

############# Was any other contrast performed better than 1v3?
size_estimates_1v32 <- as.data.frame(1:4, c("2v4", "1v2", "2v3", "3v4"))
colnames(size_estimates_1v32) <- c("1v3")
# 2v4
size_estimates_1v32[1,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))
# 1v2
size_estimates_1v32[2,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))
# 2v3
size_estimates_1v32[3,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))
# 3v4
size_estimates_1v32[4,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))
# lower CIs for each contrast
size_estimates_1v32$lowerCI <- NA
size_estimates_1v32[1,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))[1]
size_estimates_1v32[2,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))[1]
size_estimates_1v32[3,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))[1]
size_estimates_1v32[4,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))[1]
# upper CIs for each contrast
size_estimates_1v32$upperCI <- NA
size_estimates_1v32[1,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))[2]
size_estimates_1v32[2,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))[2]
size_estimates_1v32[3,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))[2]
size_estimates_1v32[4,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,2]))[2]

############# Was any other contrast performed better than 2v4?
size_estimates_2v42 <- as.data.frame(1:3, c("1v2", "2v3", "3v4"))
colnames(size_estimates_2v42) <- c("2v4")
# 1v2
size_estimates_2v42[1,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))
# 2v3
size_estimates_2v42[2,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))
#3v4
size_estimates_2v42[3,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))
# lower CIs for each contrast
size_estimates_2v42$lowerCI <- NA
size_estimates_2v42[1,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))[1]
size_estimates_2v42[2,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))[1]
size_estimates_2v42[3,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))[1]
# upper CIs for each contrast
size_estimates_2v42$upperCI <- NA
size_estimates_2v42[1,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))[2]
size_estimates_2v42[2,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))[2]
size_estimates_2v42[3,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,3]))[2]

############# Was any other contrast performed better than 1v2?
size_estimates_1v22 <- as.data.frame(1:2, c("2v3", "3v4"))
colnames(size_estimates_1v22) <- c("1v2")
# 2v3
size_estimates_1v22[1,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))
# 3v4
size_estimates_1v22[2,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))
# lower CIs for each contrast
size_estimates_1v22$lowerCI <- NA
size_estimates_1v22[1,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))[1]
size_estimates_1v22[2,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))[1]
# upper CIs for each contrast
size_estimates_1v22$upperCI <- NA
size_estimates_1v22[1,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))[2]
size_estimates_1v22[2,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,4]))[2]

############# Was 2v3 performed better than the other contrasts?
size_estimates_2v32 <- as.data.frame(1, c("3v4"))
colnames(size_estimates_2v32) <- c("2v3")
# 3v4
size_estimates_2v32[1,1] <- mean((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))
# lower CIs for each contrast
size_estimates_2v32$lowerCI <- NA
size_estimates_2v32[1,2] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))[1]
# upper CIs for each contrast
size_estimates_2v32$upperCI <- NA
size_estimates_2v32[1,3] <- HPDinterval((model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,6])-(model_size_comb_comp$Sol[,1]+model_size_comb_comp$Sol[,5]))[2]



########## time of day ############

# turn time into a numeric variable with smallest number representing the earliest start time
# and the largest number representing the latest start time
proc_no_pref_noNA_size$time2 <- as.numeric(proc_no_pref_noNA_size$start.time)

model_size_corr_time <- MCMCglmm(larger ~ time2, random = ~us(session+1):ID,
                                   family = "categorical", data = proc_no_pref_noNA_size, 
                                   prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                   verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_size_corr_time$VCV)
autocorr <- as.data.frame(autocorr(model_size_corr_time$Sol))
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_size_corr_time$Sol)
geweke.plot(model_size_corr_time$VCV)
# diagnostic for chain length
heidel.diag(model_size_corr_time$Sol)
# diagnostic for chain length and sampling
plot(model_size_corr_time$Sol)

# model results
summary(model_size_corr_time)



###### Comparison between experiments ######

# subset data to only include the data collected in the number test
N_larger_num <- subset(N_larger, N_larger$experiment == "number")
N_larger_num$experiment <- NULL
colnames(N_larger_num) <- c("discr", "ID", "larger_num", "N_trials_num", "smaller_num", "prop_num")

# subset data to only include the data collected in the size test
N_larger_size <- subset(N_larger, N_larger$experiment == "size")
N_larger_size$experiment <- NULL
colnames(N_larger_size) <- c("discr", "ID", "larger_size", "N_trials_size", "smaller_size", "prop_size")

# combine datasets
N_larger_comp <- cbind(N_larger_num, N_larger_size)

# compare each combination between the number and size test
N_larger_comp_1v4 <- subset(N_larger_comp, N_larger_comp$discr == "1v4")
wilcox.test(N_larger_comp_1v4$prop_num, N_larger_comp_1v4$prop_size, paired = FALSE, 
            alternative = "two.sided", conf.int = TRUE)

N_larger_comp_1v3 <- subset(N_larger_comp, N_larger_comp$discr == "1v3")
wilcox.test(N_larger_comp_1v3$prop_num, N_larger_comp_1v3$prop_size, paired = FALSE, 
            alternative = "two.sided", conf.int = TRUE)

N_larger_comp_1v2 <- subset(N_larger_comp, N_larger_comp$discr == "1v2")
wilcox.test(N_larger_comp_1v2$prop_num, N_larger_comp_1v2$prop_size, paired = FALSE, 
            alternative = "two.sided", conf.int = TRUE)

N_larger_comp_2v3 <- subset(N_larger_comp, N_larger_comp$discr == "2v3")
wilcox.test(N_larger_comp_2v3$prop_num, N_larger_comp_2v3$prop_size, paired = FALSE, 
            alternative = "two.sided", conf.int = TRUE)

N_larger_comp_2v4 <- subset(N_larger_comp, N_larger_comp$discr == "2v4")
wilcox.test(N_larger_comp_2v4$prop_num, N_larger_comp_2v4$prop_size, paired = FALSE, 
            alternative = "two.sided", conf.int = TRUE)

N_larger_comp_3v4 <- subset(N_larger_comp, N_larger_comp$discr == "3v4")
wilcox.test(N_larger_comp_3v4$prop_num, N_larger_comp_3v4$prop_size, paired = FALSE, 
            alternative = "two.sided", conf.int = TRUE)



#########################################################
#                  Repeatability                        #
#########################################################

################### both test pooled ####################
################# repeatability across combinations #################

model_comb_comp_Radj <- MCMCglmm(larger ~ discr, random = ~us(session+1):ID,
                                 family = "categorical", data = proc_no_pref_noNA, 
                                 prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                 verbose = FALSE, pr = TRUE)

# model diagnostics
autocorr(model_comb_comp_Radj$VCV)
autocorr <- as.data.frame(autocorr(model_comb_comp_Radj$Sol))
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_comb_comp_Radj$Sol)
geweke.plot(model_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_comb_comp_Radj)


###################### number test ######################
################# repeatability across combinations #################

# intercept model to estimate agreement repeatability
model_number_comb_comp_Ragree <- MCMCglmm(larger ~ 1, random = ~us(session+1):ID,
                                          family = "categorical", data = proc_no_pref_noNA_number, 
                                          prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                          verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_number_comb_comp_Ragree$VCV)
autocorr <- as.data.frame(autocorr(model_number_comb_comp_Ragree$Sol))
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_number_comb_comp_Ragree$Sol)
geweke.plot(model_number_comb_comp_Ragree$VCV)
# diagnostic for chain length
heidel.diag(model_number_comb_comp_Ragree$Sol)
# diagnostic for chain length and sampling
plot(model_number_comb_comp_Ragree$Sol)

# calculate repeatability
R_propA(model_number_comb_comp_Ragree)

# adjustment repeatability based on the full model (with fixed effects) see above (model_number_comb_comp)

# calculate repeatability
R_propA(model_number_comb_comp)


################# individual repeatability within combinations #################

###### 1v4
# subset data to only look at 1v4
proc_no_pref_noNA_n1v4 <- subset(proc_no_pref_noNA_number, proc_no_pref_noNA_number$discr=='1v4')
# model to estimate adjustment repeatability across individuals
model_n1v4_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_n1v4, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_n1v4_comb_comp_Radj$VCV)
autocorr <- as.data.frame(model_n1v4_comb_comp_Radj$Sol)
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_n1v4_comb_comp_Radj$Sol)
geweke.plot(model_n1v4_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_n1v4_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_n1v4_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_n1v4_comb_comp_Radj)

###### 1v3
# subset data to only look at 1v3
proc_no_pref_noNA_n1v3 <- subset(proc_no_pref_noNA_number, proc_no_pref_noNA_number$discr=='1v3')
# intercept model to estimate adjustment repeatability across individuals
model_n1v3_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_n1v3, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_n1v3_comb_comp_Radj$VCV)
autocorr <- as.data.frame(model_n1v3_comb_comp_Radj$Sol)
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_n1v3_comb_comp_Radj$Sol)
geweke.plot(model_n1v3_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_n1v3_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_n1v3_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_n1v3_comb_comp_Radj)

###### 2v4
# subset data to only look at 2v4
proc_no_pref_noNA_n2v4 <- subset(proc_no_pref_noNA_number, proc_no_pref_noNA_number$discr=='2v4')
# intercept model to estimate adjustment repeatability across individuals
model_n2v4_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_n2v4, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_n2v4_comb_comp_Radj$VCV)
autocorr <- as.data.frame(model_n2v4_comb_comp_Radj$Sol)
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_n2v4_comb_comp_Radj$Sol)
geweke.plot(model_n2v4_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_n2v4_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_n2v4_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_n2v4_comb_comp_Radj)

###### 1v2
# subset data to only look at 1v2
proc_no_pref_noNA_n1v2 <- subset(proc_no_pref_noNA_number, proc_no_pref_noNA_number$discr=='1v2')
# intercept model to estimate adjustment repeatability across individuals
model_n1v2_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_n1v2, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_n1v2_comb_comp_Radj$VCV)
autocorr <- as.data.frame(autocorr(model_n1v2_comb_comp_Radj$Sol))
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_n1v2_comb_comp_Radj$Sol)
geweke.plot(model_n1v2_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_n1v2_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_n1v2_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_n1v2_comb_comp_Radj)

###### 2v3
# subset data to only look at 2v3
proc_no_pref_noNA_n2v3 <- subset(proc_no_pref_noNA_number, proc_no_pref_noNA_number$discr=='2v3')
# intercept model to estimate adjustment repeatability across individuals
model_n2v3_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_n2v3, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_n2v3_comb_comp_Radj$VCV)
autocorr <- as.data.frame(model_n2v3_comb_comp_Radj$Sol)
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_n2v3_comb_comp_Radj$Sol)
geweke.plot(model_n2v3_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_n2v3_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_n2v3_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_n2v3_comb_comp_Radj)

###### 3v4
# subset data to only look at 3v4
proc_no_pref_noNA_n3v4 <- subset(proc_no_pref_noNA_number, proc_no_pref_noNA_number$discr=='3v4')
# intercept model to estimate adjustment repeatability across individuals
model_n3v4_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_n3v4, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_n3v4_comb_comp_Radj$VCV)
autocorr <- as.data.frame(model_n3v4_comb_comp_Radj$Sol)
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_n3v4_comb_comp_Radj$Sol)
geweke.plot(model_n3v4_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_n3v4_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_n3v4_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_n3v4_comb_comp_Radj)



###################### size test ######################
################# repeatability across combinations #################

# intercept model to estimate agreement repeatability 
model_size_comb_comp_Ragree <- MCMCglmm(larger ~ 1, random = ~us(session+1):ID,
                                        family = "categorical", data = proc_no_pref_noNA_size, 
                                        prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                        verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_size_comb_comp_Ragree$VCV)
autocorr <- as.data.frame(model_size_comb_comp_Ragree$Sol)
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_size_comb_comp_Ragree$Sol)
geweke.plot(model_size_comb_comp_Ragree$VCV)
# diagnostic for chain length
heidel.diag(model_size_comb_comp_Ragree$Sol)
# diagnostic for chain length and sampling
plot(model_size_comb_comp_Ragree$Sol)

# calculate repeatability
R_propA(model_size_comb_comp_Ragree)


# adjustment repeatability based on the full model (with fixed effects) see above (model_size_comb_comp)
# calculate repeatability
R_propA(model_size_comb_comp)


################# individual repeatability within combinations #################

###### 1v4
# subset data to only look at 1v4
proc_no_pref_noNA_s1v4 <- subset(proc_no_pref_noNA_size, proc_no_pref_noNA_size$discr=='1v4')
# intercept model to estimate adjustment repeatability across individuals
model_s1v4_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_s1v4, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_s1v4_comb_comp_Radj$VCV)
autocorr <- as.data.frame(model_s1v4_comb_comp_Radj$Sol)
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_s1v4_comb_comp_Radj$Sol)
geweke.plot(model_s1v4_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_s1v4_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_s1v4_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_s1v4_comb_comp_Radj)

###### 1v3
# subset data to only look at 1v3
proc_no_pref_noNA_s1v3 <- subset(proc_no_pref_noNA_size, proc_no_pref_noNA_size$discr=='1v3')
# intercept model to estimate adjustment repeatability across individuals
model_s1v3_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_s1v3, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_s1v3_comb_comp_Radj$VCV)
autocorr <- as.data.frame(model_s1v3_comb_comp_Radj$Sol)
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_s1v3_comb_comp_Radj$Sol)
geweke.plot(model_s1v3_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_s1v3_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_s1v3_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_s1v3_comb_comp_Radj)

###### 2v4
# subset data to only look at 2v4
proc_no_pref_noNA_s2v4 <- subset(proc_no_pref_noNA_size, proc_no_pref_noNA_size$discr=='2v4')
# intercept model to estimate adjustment repeatability across individuals
model_s2v4_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_s2v4, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_s2v4_comb_comp_Radj$VCV)
autocorr <- as.data.frame(model_s2v4_comb_comp_Radj$Sol)
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_s2v4_comb_comp_Radj$Sol)
geweke.plot(model_s2v4_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_s2v4_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_s2v4_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_s2v4_comb_comp_Radj)

###### 1v2
# subset data to only look at 1v2
proc_no_pref_noNA_s1v2 <- subset(proc_no_pref_noNA_size, proc_no_pref_noNA_size$discr=='1v2')
# intercept model to estimate adjustment repeatability across individuals
model_s1v2_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_s1v2, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_s1v2_comb_comp_Radj$VCV)
autocorr <- as.data.frame(autocorr(model_s1v2_comb_comp_Radj$Sol))
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_s1v2_comb_comp_Radj$Sol)
geweke.plot(model_s1v2_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_s1v2_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_s1v2_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_s1v2_comb_comp_Radj)

###### 2v3
# subset data to only look at 2v3
proc_no_pref_noNA_s2v3 <- subset(proc_no_pref_noNA_size, proc_no_pref_noNA_size$discr=='2v3')
# intercept model to estimate adjustment repeatability across individuals
model_s2v3_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_s2v3, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_s2v3_comb_comp_Radj$VCV)
autocorr <- as.data.frame(autocorr(model_s2v3_comb_comp_Radj$Sol))
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_s2v3_comb_comp_Radj$Sol)
geweke.plot(model_s2v3_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_s2v3_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_s2v3_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_s2v3_comb_comp_Radj)

###### 3v4
# subset data to only look at 3v4
proc_no_pref_noNA_s3v4 <- subset(proc_no_pref_noNA_size, proc_no_pref_noNA_size$discr=='3v4')
# intercept model to estimate adjustment repeatability across individuals
model_s3v4_comb_comp_Radj <- MCMCglmm(larger ~ session, random = ~us(session+1):ID,
                                      family = "categorical", data = proc_no_pref_noNA_s3v4, 
                                      prior = prior, nitt = itts, burnin = burn, thin = thins,  
                                      verbose = FALSE, pr = TRUE)
# model diagnostics
autocorr(model_s3v4_comb_comp_Radj$VCV)
autocorr <- as.data.frame(autocorr(model_s3v4_comb_comp_Radj$Sol))
# transpose data frame to be able to automate the check if any value in 
# Lags bigger than 0 are bigger than 0.1 (indication for auto correlation)
autocorr <- as.data.frame(t(autocorr))
which(autocorr$`Lag 130` > 0.1 |
        autocorr$`Lag 650` > 0.1 |
        autocorr$`Lag 1300` > 0.1|
        autocorr$`Lag 6500` > 0.1)

# set plot margins to be able to view diagnostic plots
par(mar = c(1,1,1,1))
# diagnostic to check for sufficient mixing
geweke.plot(model_s3v4_comb_comp_Radj$Sol)
geweke.plot(model_s3v4_comb_comp_Radj$VCV)
# diagnostic for chain length
heidel.diag(model_s3v4_comb_comp_Radj$Sol)
# diagnostic for chain length and sampling
plot(model_s3v4_comb_comp_Radj$Sol)

# calculate repeatability
R_propA(model_s3v4_comb_comp_Radj)



#########################################################
#                    Power Analysis                     #
#########################################################

###################### number test

# glmer of main model
model_n_1 <- glmer(larger ~ discr + (1+session|ID), data = proc_no_pref_noNA_number, family = "binomial")

# set effect size (samll effect): logit 0.45 = proportion choices to larger quantity of 0.6106392
fixef(model_n_1) <- c(0.45, 0.45, 0.45, 0.45, 0.45, 0.45)
# extend sample size to 24 individuals
model_n_1_2 <- extend(model_n_1, along = "ID", n = 24)
powerSim(model_n_1_2) # 81.90% (79.37, 84.24), Test: Likelihood ratio
# Based on 1000 simulations, alpha = 0.05


###################### size test

# glmer of main model
model_s_1 <- glmer(larger ~ discr + (1+session|ID), data = proc_no_pref_noNA_size, family = "binomial")
# set effect size (samll effect): logit 0.45 = proportion choices to larger quantity of 0.6106392
fixef(model_s_1) <- c(0.45, 0.45, 0.45, 0.45, 0.45, 0.45)
# extend sample size to 24 individuals
model_s_1_2 <- extend(model_s_1, along = "ID", n = 24)
powerSim(model_s_1_2) # 82.80% (80.32, 85.09), Test: Likelihood ratio
# Based on 1000 simulations, alpha = 0.05



#########################################################
#                          Plot                         #
#########################################################
# Figure 2 

# add proportion of choices to the larger quantity
N_larger$binom_est <- N_larger$larger/N_larger$N_trials

# remove trials 0 verus 1 and 1 versus 0 to not be included in the pot
N_larger[which(N_larger$discr == "0v1" | N_larger$discr == "1v0"),] <- NA
N_larger   <- N_larger[complete.cases(N_larger$discr), ]

# order contrasts from smalles to largest ratio/ largest to smallest distance
N_larger$discr <- factor(N_larger$discr, levels = c("1v4", "1v3", "2v4", "1v2", "2v3", "3v4"))

ggplot(N_larger, aes(discr, binom_est, fill = experiment)) +
  # draw boxplot showing data from both experiments
  geom_boxplot(outlier.size = -1, position = position_dodge(1)) +
  # fill one set of boxplots (number test) white and the other grey (size test)
  scale_fill_manual(values=c("white", "grey"), name = "") +
  # add points for means
  stat_summary(fun.y = mean, geom = "point", colour = "black", size = 1, position = position_dodge(1)) +
  # add line to indicate chance level
  geom_hline(yintercept = 0.5, linetype = 2) +
  # rescale y-axis
  scale_y_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), labels = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), limits = c(0, 1))  +
  # add axis labels
  labs(x = "Combination", y = "Probability of choosing the larger quantity") +
  # white background
  theme_bw() +
  # no grid lines, x-axis label moved further away from tick labels and font size 15pt
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.12, 0.20),
        text = element_text(size = 15))
# save plot
ggsave(filename = "SpontDiscrAll.tiff", scale = 1, width = 15, height = 15, 
       units = "cm", dpi = 300)



