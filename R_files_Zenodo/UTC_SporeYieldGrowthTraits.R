## Marta Shocket, Indiana University & UCLA
## Code for Can hot temperatures limit disease transmission? A test of mechanisms in a zooplankton-fungus system in Functional Ecology.
##
## Purpose: Analyze spore yield, spore growth, and related traits when hosts are exposed to extreme high temperatures.
##          Spore yield data come from two experiments: (1) beta + u measurement assay (20/20 and 32/32 treatments only)
##          and (2) within-host parasite growth assay (20, 26, and 32 treatments).
##
## Table of Contents:
##      1. Set working directory, load libraries and data
##      2. Spore yield: MLE model comparisons
##      3. Spore Yield Bootstrapping
##      4. Within-host parasite growth rate
##      5. Host death rate, bootstrapping, and randomizations
##      6. Host growth rate, bootstrapping, and randomizations
##      7. Relative growth calculation
##      8. Figure 3


######
##### 1. Set working directory, load libraries and data
######

# Set wd
setwd("~/Dropbox/Research Hall Lab/Temperature/Upper Thermal Constraints/Final Code and Data")

# Load library
library(bbmle)

# Load spore yield data (from beta + u measurement assay and within-host parasite growth assay)
# Various subsets are used for:
#     - final spore yield / sigma (both assays, but only hosts that died from their infection) 
#     - within-host parasite growth rate (only within-host parasite growth assay, all hosts [sacrificed + dead from infection])
#     - death rate (all data from both assays, sacrificed hosts are censogrey40 in the death rate calc)
#     - NOTE: for all traits, we only use spore yield data from dev temp = exp temp treatments 
#       (i.e., 20/20 and 32/32; not 20/32 or 32/20 in the beta + u measurement assay)
sy.data <- read.csv("UTC_CombinedSporeYieldData.csv")

# Load host growth rate data
JGR.data <- read.csv("UTC_HostGrowthRateAssay.csv")


######
##### 2. Spore yield: MLE model comparisons
######

# Subset spore yield data to include only animals that died of infection and came from constant temperature treatments
#   (i.e., only final spore yield [vs. sacrificed animals mid-infection], and not the 20/32 and 32/20 exposure/development treatments)
sy.data.sy <- subset(sy.data, StatusB == 1 & TempChange == "No")

# Subset by temperature
sy.data.20 <- subset(sy.data.sy, ExpTemp == 20)
sy.data.26 <- subset(sy.data.sy, ExpTemp == 26)
sy.data.32 <- subset(sy.data.sy, ExpTemp == 32)

########## Models comparing means of same temps from different experiments (beta + u measurement assay & within-host parasite growth assay) 
# No significant differences, so we can combine the data from both experiments

# Fit 20C models with one sd and two means, and compare models via likelihood ratio test
sy.20.1mean <- mle2(TotalSpores/10000 ~ dnorm(mean = mean, sd = sd), start = list(mean = 6, sd = 2), data = sy.data.20)
summary(sy.20.1mean)
sy.20.2means <- mle2(TotalSpores/10000 ~ dnorm(mean = mean + ExperimentB*exp, sd = sd), start = list(mean = 6, exp = 1, sd = 2), data = sy.data.20)
summary(sy.20.2means)
anova(sy.20.1mean, sy.20.2means)
summary(lm(TotalSpores/10000 ~ ExperimentB, data = sy.data.20))
summary(t.test(TotalSpores/10000 ~ ExperimentB, data = sy.data.20))

# Fit 32C models with one sd and two means, and compare models via likelihood ratio test
sy.32.1mean <- mle2(TotalSpores/10000 ~ dnorm(mean = mean, sd = sd), start = list(mean = 6, sd = 2), data = sy.data.32)
summary(sy.32.1mean)
sy.32.2means <- mle2(TotalSpores/10000 ~ dnorm(mean = mean + ExperimentB*exp, sd = sd), start = list(mean = 6, exp = 1, sd = 2), data = sy.data.32)
summary(sy.32.2means)
anova(sy.32.1mean, sy.32.2means)

########## Models comparing means + sds for different temperatures

# full model: 3 means, 3 sd
sy.full <- mle2(TotalSpores/10000 ~ dnorm(mean = mean + is26*m26 + is32*m32, sd = sd + sd26*is26 + sd32*is32), 
                     start = list(mean = 6, m26 = 1, m32 = 1, sd = 1, sd26 = .5, sd32 = 0.5), data = sy.data.sy)
summary(sy.full)

# winning model: 32 mean, 32 sd
sy.2mean.2sd <- mle2(TotalSpores/10000 ~ dnorm(mean = mean + is32*m32, sd = sd + sd32*is32), 
                            start = list(mean = 6, m32 = 1, sd = 1, sd32 = .5), data = sy.data.sy)
summary(sy.2mean.2sd)

# 3 means, 2 sd (winning model + extra 26C mean)
sy.3mean <- mle2(TotalSpores/10000 ~ dnorm(mean = mean + is26*m26 + is32*m32, sd = sd + sd32*is32), 
                    start = list(mean = 6, m26 = 1, m32 = 1, sd = 1, sd32 = 0.5), data = sy.data.sy)
summary(sy.3mean)

# 2 means, 3 sd (winning model + extra 26C sd)
sy.3sd <- mle2(TotalSpores/10000 ~ dnorm(mean = mean + is32*m32, sd = sd + sd26*is26 + sd32*is32), 
                      start = list(mean = 6, m32 = 1, sd = 1, sd26 = .5, sd32 = 0.5), data = sy.data.sy)
summary(sy.3sd)

# null model: 1 mean, 1 sd
sy.null <- mle2(TotalSpores/10000 ~ dnorm(mean = mean, sd = sd), start = list(mean = 6, sd = 1), data = sy.data.sy)

# 2 means, 1 sd
sy.1sd <- mle2(TotalSpores/10000 ~ dnorm(mean = mean + is32*m32, sd = sd), start = list(mean = 6, m32 = 1, sd = 1), data = sy.data.sy)
summary(sy.1sd)

# 1 mean, 2 sd
sy.1mean <- mle2(TotalSpores/10000 ~ dnorm(mean = mean, sd = sd + sd32*is32), start = list(mean = 6, sd = 1, sd32 = .5), data = sy.data.sy)
summary(sy.1mean)

# Compare models via likelihood ratio test
anova(sy.2mean.2sd, sy.3mean) # extra mean (26C) is not significant
anova(sy.2mean.2sd, sy.3sd) # extra sd (26C) is not significant
anova(sy.2mean.2sd, sy.full) # extra mean (26C) + extra sd (26C) together also not significant, as expected
anova(sy.2mean.2sd, sy.1mean) # different 32C mean is significant
anova(sy.2mean.2sd, sy.1sd) # different 32C sd is significant
anova(sy.2mean.2sd, sy.null) # together 32C mean and sd are significant, as expected

# Calculate delta-AIC for pairwise model comparisons
AIC(sy.2mean.2sd) - AIC(sy.3mean)
AIC(sy.2mean.2sd) - AIC(sy.3sd)
AIC(sy.2mean.2sd) - AIC(sy.full)
AIC(sy.2mean.2sd) - AIC(sy.1mean)
AIC(sy.2mean.2sd) - AIC(sy.1sd)
AIC(sy.2mean.2sd) - AIC(sy.null)

# Calculate Aikaike Weights
AIC1 <- -1*(AIC(sy.2mean.2sd) - AIC(sy.3mean))
AIC2 <- -1*(AIC(sy.2mean.2sd) - AIC(sy.3sd))
AIC3 <- -1*(AIC(sy.2mean.2sd) - AIC(sy.full))
AIC4 <- -1*(AIC(sy.2mean.2sd) - AIC(sy.1mean))
AIC5 <- -1*(AIC(sy.2mean.2sd) - AIC(sy.1sd))
AIC6 <- -1*(AIC(sy.2mean.2sd) - AIC(sy.null))
AICSum <- exp(0) + exp(-0.5*AIC1) + exp(-0.5*AIC2) + exp(-0.5*AIC3) + exp(-0.5*AIC4) + exp(-0.5*AIC5) + exp(-0.5*AIC6)

AW0 <- exp(0)/AICSum
AW1 <- exp(-0.5*AIC1)/AICSum
AW2 <- exp(-0.5*AIC2)/AICSum
AW3 <- exp(-0.5*AIC3)/AICSum
AW4 <- exp(-0.5*AIC4)/AICSum
AW5 <- exp(-0.5*AIC5)/AICSum
AW6 <- exp(-0.5*AIC6)/AICSum

sum(AW0, AW1, AW2, AW3, AW4, AW5, AW6) # check to make sure it sums to 1


######
##### 3. Bootstrap spore yield (for R0 calculations and plotting 95% CIs)
######

################## Function - BoostrapSY
# Arguments: 
    # data.set    dataframe: spore yield data
    # nBoots      integer: number of bootstraping iterations
# Returns: 
    # sy.means    vector: bootstrapped mean spore yield values

BootstrapSY = function(data.set, nBoots) {
  
  # get data
  data <- data.set
  
  # create output vector
  sy.means <- numeric(nBoots)
  
  # Loop through iterations
  for (i in 1:nBoots){
    
    # Sample data
    sample.list <- sample(data$TotalSpores, replace = TRUE)
    
    # Calculate mean
    sy.means[i] <- mean(sample.list)
    
  }
  
  sy.means # return output
}

# Bootstrap spore yield
bs.sy.20 <- BootstrapSY(sy.data.20, 10000)
bs.sy.26 <- BootstrapSY(sy.data.26, 10000)
bs.sy.32 <- BootstrapSY(sy.data.32, 10000)

bs.sy <- data.frame(sy20 = bs.sy.20, sy26 = bs.sy.26, sy32 = bs.sy.32)
save(bs.sy, file = "bs.sy.Rsave")
load("bs.sy.Rsave")

# Create output dataframe for plotting
sy.out <- data.frame(Temp = c(20, 26, 32), Trait = c("sigma", "sigma", "sigma"), mean = numeric(3), upper = numeric(3), lower = numeric(3))

# Calculate and store putput
sy.out$mean[1] <- mean(sy.data.20$TotalSpores)
sy.out$mean[2] <- mean(sy.data.26$TotalSpores)
sy.out$mean[3] <- mean(sy.data.32$TotalSpores)
sy.out$lower[1] <- quantile(bs.sy$sy20, 0.025)
sy.out$upper[1] <- quantile(bs.sy$sy20, 0.975)
sy.out$lower[2] <- quantile(bs.sy$sy26, 0.025)
sy.out$upper[2] <- quantile(bs.sy$sy26, 0.975)
sy.out$lower[3] <- quantile(bs.sy$sy32, 0.025)
sy.out$upper[3] <- quantile(bs.sy$sy32, 0.975)


#####
#### 4. Within host parasite growth rate
#####

# Subset spore yield data to only the within-host growth rate experiment and by temperature
sy.data.WHG.20 <- subset(sy.data, Experiment == "WHG" & ExpTemp == 20)
sy.data.WHG.26 <- subset(sy.data, Experiment == "WHG" & ExpTemp == 26)
sy.data.WHG.32 <- subset(sy.data, Experiment == "WHG" & ExpTemp == 32)

# Trim the 20C data because it plateaus after 19 days and we are just trying to capture the growth period
sy.data.WHG.20  <- subset(sy.data.WHG.20, Day <= 19)

########## Linear Models

# 20C
lm.WHG.20 <- lm(TotalSpores ~ Day, data = sy.data.WHG.20)
summary(lm.WHG.20)
cfs.20 <- coef(lm.WHG.20)

# 26C
lm.WHG.26 <- lm(TotalSpores ~ Day, data = sy.data.WHG.26)
summary(lm.WHG.26)
cfs.26 <- coef(lm.WHG.26)

# 26C - fixed x-intercept = 3 
lm.WHG.26.fi <- lm(TotalSpores ~ I(Day - 3) - 1, data = sy.data.WHG.26)
summary(lm.WHG.26.fi)
cf.26.fi <- coef(lm.WHG.26.fi)

# 32C
lm.WHG.32 <- lm(TotalSpores ~ Day, data = sy.data.WHG.32)
summary(lm.WHG.32)
cfs.32 <- coef(lm.WHG.32)

AIC(lm.WHG.26.fi)
AIC(lm.WHG.26)
-cfs.26[1]/cfs.26[2]

# Function to calculate daily spore averages for plotting only (linear models fit to raw data)
CalcDailySporeAvg = function(data.set){
  data.set <- data.set
  day.list <- unique(data.set$Day)
  out <- data.frame(Day = numeric(length(day.list)), mean = numeric(length(day.list)), sd = numeric(length(day.list)), SE = numeric(length(day.list)), 
                    n = numeric(length(day.list)), upper = numeric(length(day.list)), lower = numeric(length(day.list)))
  for(i in 1:length(day.list)){
    data.sub <- subset(data.set, Day == day.list[i])
    out$Day[i] <- day.list[i]
    out$mean[i] <- mean(data.sub$TotalSpores)
    out$sd[i] <- sd(data.sub$TotalSpores)
    out$n[i] <- nrow(data.sub)
    out$SE[i] <- sd(data.sub$TotalSpores) / sqrt(nrow(data.sub))
    out$lower[i] <- mean(data.sub$TotalSpores) - sd(data.sub$TotalSpores) / sqrt(nrow(data.sub))
    out$upper[i] <- mean(data.sub$TotalSpores) + sd(data.sub$TotalSpores) / sqrt(nrow(data.sub))
  }
  out
}

# For graphical clarity, we're going to bin a few animals were the only animals to die on a given day 
#       with animals that died on an adjacent day
sy.data.WHG.20.bin <- sy.data.WHG.20
sy.data.WHG.26.bin <- sy.data.WHG.26
sy.data.WHG.32.bin <- sy.data.WHG.32
sy.data.WHG.20.bin$Day[24] <- 12 # lone animal that died on day 13 with day 12 animals
sy.data.WHG.20.bin$Day[35] <- 16 # lone animal that died on day 15 with lone animal that died on day 16
sy.data.WHG.26.bin$Day[10] <- 8 # lone animal that died on day 9 with day 8 animals
sy.data.WHG.32.bin$Day[10] <- 8 # lone animal that died on day 9 with day 8 animals

# Calculate daily averages
DSA.20 <- CalcDailySporeAvg(sy.data.WHG.20.bin)
DSA.26 <- CalcDailySporeAvg(sy.data.WHG.26.bin)
DSA.32 <- CalcDailySporeAvg(sy.data.WHG.32.bin)

save(DSA.20, file = "DSA.20.Rsave")
save(DSA.26, file = "DSA.26.Rsave")
save(DSA.32, file = "DSA.32.Rsave")

########## Function: 
BootstrapPGR = function(data.set, nBoots, fi) {
  
  # get data
  data <- data.set
  
  # create output vector
  PGR <- numeric(nBoots)
  
  # Loop through iterations
  for (i in 1:nBoots){
    
    # Sample data
    data.sim <- data[sample(nrow(data), replace = TRUE), ]
    
    if(fi == TRUE){
      # Fit lm and store PGR with fixed x-intercept at day 3 for 26C
      PGR.lm <- lm(TotalSpores ~ I(Day - 3) - 1, data = data.sim)
      PGR[i] <- coef(PGR.lm)
    }
    else{
      # Fit lm and store PGR without fixed intercept
      PGR.lm <- lm(TotalSpores ~ Day, data = data.sim)
      PGR[i] <- coef(PGR.lm)[2]
    }
    
  }
  
  PGR # return output
}

bs.PGR.20 <- BootstrapPGR(sy.data.WHG.20, 10000, fi = FALSE)
bs.PGR.26 <- BootstrapPGR(sy.data.WHG.26, 10000, fi = TRUE)
bs.PGR.32 <- BootstrapPGR(sy.data.WHG.32, 10000, fi = FALSE)
bs.PGR <- data.frame(PGR20 = bs.PGR.20, PGR26 = bs.PGR.26, PGR32 = bs.PGR.32)

# Create output dataframe
PGR.out <- data.frame(Temp = c(20, 26, 32), Trait = c("PGR", "PGR", "PGR", "PGR.i", "PGR.i", "PGR.i"), mean = numeric(6), upper = numeric(6), lower = numeric(6))

# Add PGR (slope) estimates from linear model
PGR.out$mean[1] <- cfs.20[2]
PGR.out$mean[2] <- cf.26.fi
PGR.out$mean[3] <- cfs.32[2]

# Add 95% CIs from bootstrapped linear models
PGR.out$lower[1] <- quantile(bs.PGR.20, 0.025)
PGR.out$lower[2] <- quantile(bs.PGR.26, 0.025)
PGR.out$lower[3] <- quantile(bs.PGR.32, 0.025)
PGR.out$upper[1] <- quantile(bs.PGR.20, 0.975)
PGR.out$upper[2] <- quantile(bs.PGR.26, 0.975)
PGR.out$upper[3] <- quantile(bs.PGR.32, 0.975)

# Add y-intercept estimates for 20 and 32
PGR.out$mean[4] <- cfs.20[1]
PGR.out$mean[6] <- cfs.32[1]

1 - ecdf(bs.PGR$PGR20)(PGR.out$mean[2]) # 20 vs. 26
ecdf(bs.PGR$PGR26)(PGR.out$mean[1]) # 20 vs. 26

1 - ecdf(bs.PGR$PGR20)(PGR.out$mean[3]) # 20 vs. 32
ecdf(bs.PGR$PGR32)(PGR.out$mean[1]) # 20 vs. 32

ecdf(bs.PGR$PGR26)(PGR.out$mean[3]) # 26 vs. 32
1 - ecdf(bs.PGR$PGR32)(PGR.out$mean[2]) # 32 vs. 26

save(bs.PGR, file = "bs.PGR.Rsave")


############ Expontential Growth Models for model comparison

# Fit exponential (log-linear) models
exp.WHG.20 <- lm(log(TotalSpores+1) ~ Day, data = sy.data.WHG.20)
exp.WHG.26 <- lm(log(TotalSpores+1) ~ Day, data = sy.data.WHG.26)
exp.WHG.32 <- lm(log(TotalSpores+1) ~ Day, data = sy.data.WHG.32)

# Extract coefficients
cf.exp.WHG.20 <- coef(exp.WHG.20)
cf.exp.WHG.26 <- coef(exp.WHG.26)
cf.exp.WHG.32 <- coef(exp.WHG.32)

# Get R^2 values from summaries
summary(exp.WHG.20)
summary(exp.WHG.26)
summary(exp.WHG.32)
summary(lm.WHG.20)
summary(lm.WHG.26)
summary(lm.WHG.32)

par(mfrow = c(2,3), mar = c(4.5, 4.5, 2, 1))

# Regular plots, both linear and exponential models
plot(TotalSpores ~ Day, data = sy.data.WHG.20, col = "blue", main = "20C")
abline(cfs.20[1], cfs.20[2])
curve(exp(cf.exp.WHG.20[1]+cf.exp.WHG.20[2]*x), add = TRUE, lty = 2)

plot(TotalSpores ~ Day, data = sy.data.WHG.26, col = "purple", main = "26C")
abline(cfs.26[1], cfs.26[2])
lines(x = c(3, 17), y = c(0, cf.26.fi*(17-3)), lwd = 2)
curve(exp(cf.exp.WHG.26[1]+cf.exp.WHG.26[2]*x), add = TRUE, lty = 2)
legend("topright", bty = "n", lty = c(1, 1, 2), lwd = c(1, 2, 1), legend = c("lin.", "lin. (fi)", "exp."), cex = 1.3)

plot(TotalSpores ~ Day, data = sy.data.WHG.32, col = "red", main = "32C")
abline(cfs.32[1], cfs.32[2])
curve(exp(cf.exp.WHG.32[1]+cf.exp.WHG.32[2]*x), add = TRUE, lty = 2)

# Log plots, with exponential (log-linear) models
plot(log(TotalSpores+1) ~ Day, data = sy.data.WHG.20, col = "blue")
abline(cf.exp.WHG.20[1], cf.exp.WHG.20[2], lty = 2)
plot(log(TotalSpores+1) ~ Day, data = sy.data.WHG.26, col = "purple")
abline(cf.exp.WHG.26[1], cf.exp.WHG.26[2], lty = 2)
plot(log(TotalSpores+1) ~ Day, data = sy.data.WHG.32, col = "red")
abline(cf.exp.WHG.32[1], cf.exp.WHG.32[2], lty = 2)


#####
#### 5. Host death rate
#####

# Re-subset the spore yield data to include all animals regardless of death observation (since we can censor the death rate calculation)
#    but still exclude the 20/32 and 32/20 exposure/development treatments
d.data <- subset(sy.data, TempChange == "No")

d.data.20 <- subset(d.data, ExpTemp == 20)
d.data.26 <- subset(d.data, ExpTemp == 26)
d.data.32 <- subset(d.data, ExpTemp == 32)

################## Function - d.NLL
# Calculates NLL (Negative Log Likelihood) for parameter d, the constant failure rate in negative exponential distribution
# Arguments: 
  # data.set    dataframe:       data with day of spore count and observation status (alive/sacrificed or died from infection)
  # d           numeric value:   d parameter from optimizing function
# Returns: 
  # NLL         numeric value:   negative log-likelihood for parameter d

d.NLL = function(d, data.set) {
  
  d=exp(d) # Apply exp transform to ensure positive death rate - the little.d function takes this transformation into account
  
  observed = data.set$StatusB # 1 if the animal died from infection, 0 if sacrificed while still alive
  
  day = data.set$Day # The day of the spore count (day of death or sacrifice)
  
  NLL = 0 # Set NLL to zero, will be cumulatively adding likelihood values for each observation/animal
  
  # Loop through each observation/animal and add its likelihood value to the total NLL
  # For animals that die, add value of the density function for the day they died: d*exp(-d*day)
  # For censogrey40 animals that didn't die, add the survivorship function for the final day: exp(-d*day) #(exp and log cancel out)
  for (i in 1:length(observed)) {
    ifelse(observed[i]==1, NLL <- NLL - log(d*exp(-d*day[i])), NLL <- NLL - (-d*day[i]))
  }
  
  NLL # Return
  
}

################## Function - little.d
# Uses MLE to optimize death rate (constant failure rate, d) for an exponential survival curve 
# Arguments: 
  # data.set    dataframe:       data with day of spore count and observation status (alive/sacrificed or died from infection)
# Returns: 
  # d           numeric value:   bets fit for parameter d

little.d = function(data.set) {
  
  # read in dataset (dataset = data frame from within R)
  LT = data.set
  
  # Use MLE to calculate most likely value of d, given the data. If nobody died, death rate = 0.
  # The d.NLL transforms d = exp(d), so the output here also transforms d before returning the parameter
  ifelse(mean(LT$StatusB) == 0, d <- 0, d <- exp(coef(mle2(minuslogl=d.NLL, start=list(d=log(0.05)), data=list(data.set=LT), skip.hessian=T))))
  
  d # Return
  
}

################## Function - Bootstrap_d - bootstrap estimates of d
Bootstrap_d = function (data.set, status.col, nBoots) {
  
  # create output vector
  out.vec = numeric(nBoots)
  
  # Loop through nBoots
  for (i in 1:nBoots) {
    
    # Simulate a dataset
    simulated.data.set = data.set[sample(x = nrow(data.set), size = nrow(data.set), replace = TRUE), ]
    
    # calculate and store death rate for the simulated data set
    out.vec[i] = little.d(simulated.data.set)
    
  }
  
  out.vec   # return
  
}

################## Function - Rand_d - randomization test for d
Rand_d = function (dataset, treat.1, treat.2, nBoots) {
  
  # create output vector
  out.vec = numeric(nBoots)
  
  # Loop through nBoots
  for (i in 1:nBoots) {
    
    # make temp dataset for loop
    simulated.data.set <- dataset
    
    # add randomized treatment labels
    simulated.data.set$Treat.rand <- sample(simulated.data.set$ExpTemp, replace = FALSE)    
    
    # subset by randomized treatment labels
    sim.data.1 <- subset(simulated.data.set, Treat.rand == treat.1)
    sim.data.2 <- subset(simulated.data.set, Treat.rand == treat.2)
    
    # calculate and store death rates
    death.1 = little.d(sim.data.1)
    death.2 = little.d(sim.data.2)
    
    # calculate and store death rate for the simulated data set
    out.vec[i] = abs(death.1 - death.2)
    
  }
  
  out.vec # return
  
}

# Bootstrap d
bs.d.20 <- Bootstrap_d(d.data.20, 3, 10000)
bs.d.26 <- Bootstrap_d(d.data.26, 3, 10000)
bs.d.32 <- Bootstrap_d(d.data.32, 3, 10000)

# Create output dataframe for plotting and calculate best estimate of d
d.out <- data.frame(Temp = c(20, 26, 32), Trait = c("d", "d", "d"), mean = numeric(3), upper = numeric(3), lower = numeric(3))
d.out$mean[1] <- little.d(d.data.20)
d.out$mean[2] <- little.d(d.data.26)
d.out$mean[3] <- little.d(d.data.32)

# Calculate and store quantiles in output dataframe
d.out$lower[1] <- quantile(bs.d.20, 0.025)
d.out$lower[2] <- quantile(bs.d.26, 0.025)
d.out$lower[3] <- quantile(bs.d.32, 0.025)
d.out$upper[1] <- quantile(bs.d.20, 0.975)
d.out$upper[2] <- quantile(bs.d.26, 0.975)
d.out$upper[3] <- quantile(bs.d.32, 0.975)

# Subset data for randomzation
d.data.2026 <- subset(d.data, ExpTemp == 20 | ExpTemp == 26)
d.data.2632 <- subset(d.data, ExpTemp == 26 | ExpTemp == 32)
d.data.2032 <- subset(d.data, ExpTemp == 20 | ExpTemp == 32)

# Randomization / Permutation tests for d
rand.d.2026 <- Rand_d(d.data.2026, 20, 26, 10000) 
rand.d.2632 <- Rand_d(d.data.2632, 26, 32, 10000) 
rand.d.2032 <- Rand_d(d.data.2032, 20, 32, 10000) 

# Calculate differences between best estimates of d
dd.2026 <- d.out$mean[2] - d.out$mean[1]
dd.2632 <- d.out$mean[3] - d.out$mean[2]
dd.2032 <- d.out$mean[3] - d.out$mean[1]

# Calculate p-values for d randomization tests from the cumulative probability distribution density
1 - ecdf(rand.d.2026)(dd.2026) # 0.02096147: p = 0
1 - ecdf(rand.d.2632)(dd.2632) # 0.009968539: p = 0.063
1 - ecdf(rand.d.2032)(dd.2032) # 0.03093001: p = 0

d.bs <- cbind(bs.d.20, bs.d.26, bs.d.32)
save(d.bs, file = "d.bs.Rsave")


#####
#### 6. Host growth rate (JGR = juvenile growth rate)
#####

# Subset JGR data by temperature
JGR.20 <- subset(JGR.data, Temp == 20)
JGR.26 <- subset(JGR.data, Temp == 26)
JGR.32 <- subset(JGR.data, Temp == 32)

# Create output dataframe and calculate means
JGR.out <- data.frame(Temp = c(20, 26, 32), Trait = c("JGR", "JGR", "JGR"), mean = numeric(3), upper = numeric(3), lower = numeric(3))
JGR.out$mean[1] <- mean(JGR.20$JGR)
JGR.out$mean[2] <- mean(JGR.26$JGR)
JGR.out$mean[3] <- mean(JGR.32$JGR)

################## Function - BootstrapJGR
BootstrapJGR = function(data.vec, nBoots) {
  
  # get data
  data.vec <- data.vec
  
  # create output vector
  JGR.means <- numeric(nBoots)
  
  # Loop through iterations
  for (i in 1:nBoots){
    
    # Sample data
    sample.list <- sample(data.vec, replace = TRUE)
    
    # Calculate mean
    JGR.means[i] <- mean(sample.list)
    
  }

  JGR.means   #output
  
}

# Bootstrap JGR
bs.JGR.20 <- BootstrapJGR(JGR.20$JGR, 10000)
bs.JGR.26 <- BootstrapJGR(JGR.26$JGR, 10000)
bs.JGR.32 <- BootstrapJGR(JGR.32$JGR, 10000)
bs.JGR <- data.frame(JGR20 = bs.JGR.20, JGR26 = bs.JGR.26, JGR32 = bs.JGR.32)
save(bs.JGR, file = "bs.JGR.Rsave")

# Save output
JGR.out$lower[1] <- quantile(bs.JGR.20, 0.025)
JGR.out$lower[2] <- quantile(bs.JGR.26, 0.025)
JGR.out$lower[3] <- quantile(bs.JGR.32, 0.025)
JGR.out$upper[1] <- quantile(bs.JGR.20, 0.975)
JGR.out$upper[2] <- quantile(bs.JGR.26, 0.975)
JGR.out$upper[3] <- quantile(bs.JGR.32, 0.975)


#####
#### 7. Relative growth calculation - how does spore yield scale with host growth rate? - ultimately not included in paper
#####

# Create output dataframe
RG.out <- data.frame(Temp = c(20, 26, 32), Trait = c("RG", "RG", "RG"), 
                      mean = numeric(3), upper = numeric(3), lower = numeric(3))

# Calculate means and add to output
RG.out$mean[1] <- sy.out$mean[1] / JGR.out$mean[1]
RG.out$mean[2] <- sy.out$mean[2] / JGR.out$mean[2]
RG.out$mean[3] <- sy.out$mean[3] / JGR.out$mean[3]

# Calculate bootstrapped RG from bootstrapped spore yield and JGR
bs.RG <- bs.sy / bs.JGR
colnames(bs.RG) <- c("RG.20", "RG.26", "RG.32")

# Calculate RG quantiles and add to output
RG.out$lower[1] <- quantile(bs.RG$RG.20, 0.025)
RG.out$lower[2] <- quantile(bs.RG$RG.26, 0.025)
RG.out$lower[3] <- quantile(bs.RG$RG.32, 0.025)
RG.out$upper[1] <- quantile(bs.RG$RG.20, 0.975)
RG.out$upper[2] <- quantile(bs.RG$RG.26, 0.975)
RG.out$upper[3] <- quantile(bs.RG$RG.32, 0.975)

# Calculate p-values by comparing best estimates to cumulative probability distributions of bootstrapped values of different temperatures
#    Note: For correct p-value use the actual cumulative probability density if the distribution mean is > the estimate
#          but use "1 -" the cumulative probability density if the estimate is > the distribution mean
ecdf(bs.RG$RG.20)(RG.out$mean[2]) # 20 vs. 26
1 - ecdf(bs.RG$RG.26)(RG.out$mean[1]) # 26 vs. 20

ecdf(bs.RG$RG.20)(RG.out$mean[3]) # 20 vs. 32
1 - ecdf(bs.RG$RG.32)(RG.out$mean[1]) # 32 vs. 20

ecdf(bs.RG$RG.26)(RG.out$mean[3]) # 26 vs. 32
1 - ecdf(bs.RG$RG.32)(RG.out$mean[2]) # 32 vs. 26

##### Collected all traits in a single data frame
spore.traits.out <- rbind(sy.out, PGR.out, d.out, JGR.out, RG.out)
save(spore.traits.out, file = "sporetraits.out.Rsave")


######
##### 8. Figure 3
######

# Load and subset saved results for plotting
load("sporetraits.out.Rsave")
sy.out <- subset(spore.traits.out, Trait == "sigma")
JGR.out <- subset(spore.traits.out, Trait == "JGR")
PGR.out <- subset(spore.traits.out, Trait == "PGR" | Trait == "PGR.i")
d.out <- subset(spore.traits.out, Trait == "d")

# Load daily average spore yields for plotting within-host parasite growth
load("DSA.20.Rsave")
load("DSA.26.Rsave")
load("DSA.32.Rsave")

DSA.26$Day <- DSA.26$Day + 0.3

##### plot settings - 2 2x2 grids put together in powerpoint
par(mfrow = c(3,2), las = 1, oma = c(0, 0, 0, 0), mar = c(5, 5, 1, 1), mgp = c(3.5, 1, 0))

frame()

##### Panel A: Final spore yield (sigma)
plot(mean/10000 ~ Temp, xlab = expression(paste("Maximum temperature (",degree,"C)")), cex.axis = 1,
     ylab = expression(paste("Final spore yield - ",sigma," (10"^-4,")")), ylim = c(5.5, 9.5),
     data = sy.out, type = "p", pch = 21, bg = "white", xaxt = "n", xlim = c(18,34), cex.lab = 1.25, cex.axis = 1)
axis(side = 1, at = c(20, 26, 32), cex.axis = 1)
arrows(sy.out$Temp, sy.out$lower/10000, sy.out$Temp, sy.out$upper/10000, angle = 90, length = 0, code = 3, col = "black", lwd = 1.4)
points(mean/10000 ~ Temp, data = sy.out, type = "p", pch = 21, col = "black", bg = c("white", "grey80", "grey40"), cex = 1.75)
legend("topleft", legend = "A", cex = 1.25, bty = "n", adj = 1)
text(x = 18.6, y = 7.2, labels = "a", cex = 1.1)
text(x = 24.6, y = 7.8, labels = "a", cex = 1.1)
text(x = 33.2, y = 6.33, labels = "b", cex = 1.1)
text(x = 21.75, y = 5.65, labels = "+/- 95% CIs", cex = 1.1)

##### Panel B: host growth rate (JGR /  gH)
plot(mean ~ Temp, data = JGR.out, ylim = c(0.3, 0.5), xlim = c(18,34), type = "n", lwd = 1, col = "grey80",
     ylab = expression(paste("Host growth rate - ",italic(g[h])," (day"^-1,")")), xlab = expression(paste("Maximum temperature (",degree,"C)")), 
     lty = 1, cex.lab = 1.25, cex.axis = 1, xaxt = "n")
axis(side = 1, at = c(20, 26, 32))
arrows(JGR.out$Temp, JGR.out$lower, JGR.out$Temp, JGR.out$upper, angle = 90, length = 0, code = 3, lwd = 1.4, col = "black")
points(mean ~ Temp, data = JGR.out, pch = 21, col = "black", bg = c("white", "grey80", "grey40"), cex = 1.75)
legend("topleft", legend = "B", bty = "n", cex = 1.25, adj = 1)
text(x = 21.75, y = 0.305, labels = "+/- 95% CIs", cex = 1.1)
text(x = 19, y = 0.363, labels = "a", cex = 1.1)
text(x = 25, y = 0.447, labels = "b", cex = 1.1)
text(x = 31, y = 0.47, labels = "c", cex = 1.1)

##### Panel C: death rate (d)
plot(mean ~ Temp, data = d.out, xlim = c(18, 34), ylim = c(0.02, 0.07), xlab = expression(paste("Maximum temperature (",degree,"C)")), 
     ylab = expression(paste("Host death rate - ",italic(d)," (day"^-1,")")),  cex.lab = 1.25, cex.axis = 1, xaxt = "n")
axis(side = 1, at = c(20, 26, 32), cex.axis = 1)
arrows(d.out$Temp, d.out$upper, d.out$Temp, d.out$lower, length = 0, code = 3, lwd = 1.4, col = "black")
points(mean ~ Temp, data = d.out, cex = 1.75, pch = 21, col = "black", bg = c("white", "grey80", "grey40"))
legend("topleft", legend = "C", bty = "n", cex = 1.25, adj = 1)
text(x = 19, y = .032, labels = "a", cex = 1.1)
text(x = 25, y = .0532, labels = "b", cex = 1.1)
text(x = 31, y = .0631, labels = "b", cex = 1.1)
text(x = 20.75, y = .022, labels = "+/- 95% CIs", cex = 1.1)

##### Panel D: Within host parasite growth - spores over time
plot(lower/10000 ~ Day, data=DSA.20, xlim = c(2.5, 19.5), type = "n", ylab = expression(paste("Spore load within hosts (10"^-4,")")), 
     ylim = c(0, 14), xlab = "Time post exposure (days)", cex.lab = 1.25, cex.axis = 1)
segments(x0 = -PGR.out$mean[4]/PGR.out$mean[1], x1 = 19, y0 = 0, y1 = (PGR.out$mean[4]+PGR.out$mean[1]*19)/10000, lwd = 2, lty = 1, col = "black")
segments(x0 = 3, x1 = 16, y0 = 0, y1 = PGR.out$mean[2]*(16-3)/10000, lwd = 2, lty = 2, col = "black")
segments(x0 = -PGR.out$mean[6]/PGR.out$mean[3], x1 = 14, y0 = 0, y1 = (PGR.out$mean[6]+PGR.out$mean[3]*14)/10000, lwd = 2, lty = 3, col = "black")
arrows(DSA.20$Day, DSA.20$lower/10000, DSA.20$Day, DSA.20$upper/10000, col = "black", length = 0, lwd = 1.4)
arrows(DSA.26$Day, DSA.26$lower/10000, DSA.26$Day, DSA.26$upper/10000, col = "black", length = 0, lwd = 1.4)
arrows(DSA.32$Day, DSA.32$lower/10000, DSA.32$Day, DSA.32$upper/10000, col = "black", length = 0, lwd = 1.4)
points(mean/10000 ~ Day, data = DSA.20, pch = 21, col = "black", bg = "white", cex = 1.5)
points(mean/10000 ~ Day, data = DSA.26, pch = c(rep(21,nrow(DSA.26)-1),22), col = "black", bg = "grey80", cex = 1.5)
points(mean/10000 ~ Day, data = DSA.32, pch = c(rep(21,nrow(DSA.32)-1),22), col = "black", bg = "grey40", cex = 1.5)
a <- expression(paste("20",degree,"C"))
b <- expression(paste("26",degree,"C"))
c <- expression(paste("32",degree,"C"))
# text(x = 8, y = 12.5, labels = expression(paste("High temperature (",degree,"C):")), cex = 1)
legend(x = 2, y = 13, legend = c(b, c, a), lty = c(2, 3, 1), bty = "n", lwd = 2, cex = 1.2)
legend(x = 4, y = 13, legend = c("", "", ""), pch = c(21), pt.bg = c("grey80", "grey40", "white"), bty = "n", cex = 1.2, pt.cex = 1.5)
legend("topleft", legend = "D", bty = "n", cex = 1.25, adj = 1)
text(x = 15.75, y = .5, labels = "+/- SE", cex = 1.1)

##### Panel E: Within host parasite growth - growth rate plot
plot(mean/1000 ~ Temp, data = PGR.out, xlim = c(18, 34),  ylim = c(5, 13.5), xlab = expression(paste("Maximum temperature (",degree,"C)")),  
     ylab = expression(paste("Parasite growth rate - ",italic(g[p]))),  cex.lab = 1.25, cex.axis = 1, xaxt = "n", type = "n")
axis(side = 1, at = c(20, 26, 32))
mtext(text = expression(paste("(spores/day 10"^-3,")")), line = 2.25, side = 2, cex = 0.85, las = 0)
arrows(PGR.out$Temp, PGR.out$upper/1000, PGR.out$Temp, PGR.out$lower/1000, angle = 90, length = 0, code = 3, col = "black", lwd = 1.4)
points(mean/1000 ~ Temp, data = PGR.out, pch = 21, cex = 1.75, col = "black", bg = c("white", "grey80", "grey40"))
legend("topleft", legend = "E", bty = "n", cex = 1.25, adj = 1)
text(x = 19, y = 8.55, labels = "a", cex = 1.1)
text(x = 25, y = 9.25, labels = "a", cex = 1.1)
text(x = 31, y = 9.25, labels = "a", cex = 1.1)
text(x = 21.75, y = 5.3, labels = "+/- 95% CIs", cex = 1.1)