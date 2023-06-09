library(plyr)
library(reshape2)
library(Rmisc)
library(ggplot2)
library(cowplot)

#### FECUNDITY (f) data cleaning ####

# Read fecundity data; selecting relevant columns
f.data <- read.csv("EvolExp_fecundity_data.csv")
names(f.data)[3] <- "rep" #NB: what I called 'Treat_rep' in data collection will be referred to as simply 'rep' for all analysis purposes
f.data$rep <- as.factor(f.data$rep) 
f.data$Treat_rep <- paste(f.data$Treatment,f.data$rep) 
#NB: concatenating Treatment name and 'rep' to make unique identifier for each experimental population
#NB: not to be confused with 'Treat_rep' designated during data collection & in raw data file

# Event data processing and cleaning

# Convert all clutch production event records to Julian values
# NB: "AMPM" columns differentiate at which one of the two records per day (12hr monitoring intervals) events occured
f.data$clutch_1_date <- as.Date(f.data$clutch_1_date, "%m/%d/%y")
f.data$julian_1 <- julian(f.data$clutch_1_date) # Adding column corresponding to Julian date of first clutch (origin used: 1970 Jan 01)
f.data$daynight_1 <- as.numeric(as.character(mapvalues(f.data$AMPM_1, c("AM", "PM"), c("0.0","0.5")))) # Converting AM/PM into numerical
f.data$timing_1 <- as.numeric(f.data$julian_1 + f.data$daynight_1) # Numerical value of clutch timing

f.data$clutch_2_date <- f.data$clutch_2_date <- as.Date(f.data$clutch_2_date, "%m/%d/%y")
f.data$julian_2 <- julian(f.data$clutch_2_date) 
f.data$daynight_2 <- as.numeric(as.character(mapvalues(f.data$AMPM_2, c("AM", "PM"), c("0.0","0.5")))) 
f.data$timing_2 <- as.numeric(f.data$julian_2 + f.data$daynight_2) 

f.data$clutch_3_date <- f.data$clutch_3_date <- as.Date(f.data$clutch_3_date, "%m/%d/%y")
f.data$julian_3 <- julian(f.data$clutch_3_date) 
f.data$daynight_3 <- as.numeric(as.character(mapvalues(f.data$AMPM_3, c("AM", "PM"), c("0.0","0.5")))) 
f.data$timing_3 <- as.numeric(f.data$julian_3 + f.data$daynight_3)

# Calculate intervals between 3 clutch events; get mean
f.data$first_int <- as.numeric(f.data$timing_2 - f.data$timing_1)
f.data$second_int <- as.numeric(f.data$timing_3 - f.data$timing_2)
f.data$mean_int <- rowMeans(f.data[,c("first_int", "second_int")], na.rm = TRUE)

# Convert mean interval per individual to rate of reproduction (multiplied by mean clutch size, identical to 2017 study)
f.data$f <- 1/f.data$mean_int * 47.32

# Collect and clean all replicates with measured 'f', for statistical analysis on groups
f.reps <- as.data.frame(cbind(f.data[,c(2,1,3,28)]))
# NB: "Source_pop" denotes sourced stock pop container from pre-exp rearing, not used in analysis (just used for extra measure of replicability)

# Subset replicate populations with non-NaN values
f.reps <- f.reps[complete.cases(f.reps),]

# Calculate summary statistics for all Treatments and reps
f.summary <- summarySE(f.reps, measurevar = "f", groupvars = c("Treatment", "rep"), na.rm = T)


#### MATURATION RATE (mu) data cleaning ####

# Reading maturation data; selecting relevant columns
mu.data <- read.csv("EvolExp_maturation_data.csv")
# NB: start of life (calling it 'init' here), i.e. first clutch, will be borrowed from "f" dataset and appended
init <- as.data.frame(f.data[,1:5])
mu.data <- cbind(init, mu.data[,5:9])
mu.data$rep <- as.factor(mu.data$rep)
mu.data$Treat_rep <- paste(mu.data$Treatment,mu.data$rep) 
#NB: concatenating Treatment name and 'rep' to make unique identifier for each experimental population, same as for 'f' above

# Event data processing and cleaning

# Convert all maturation event records to Julian values
# NB: "AMPM" columns differentiate at which one of the two records per day (12hr monitoring intervals) events occured
mu.data$birth <- julian(as.Date(mu.data$clutch_1_date, "%m/%d/%y")) 

mu.data$maturity_1_timing <- julian(as.Date(mu.data$Gravid_1_date, "%m/%d/%y"))
mu.data$maturity_2_timing <- julian(as.Date(mu.data$Gravid_2_date, "%m/%d/%y"))
mu.data$maturity_3_timing <- julian(as.Date(mu.data$Gravid_3_date, "%m/%d/%y"))
mu.data$maturity_4_timing <- julian(as.Date(mu.data$Gravid_4_date, "%m/%d/%y"))
mu.data$maturity_5_timing <- julian(as.Date(mu.data$Gravid_5_date, "%m/%d/%y"))

# Calculate 2 maturation times per clutch; get mean
mu.data$age_maturity1 <- as.numeric(mu.data$maturity_1_timing - mu.data$birth)
mu.data$age_maturity2 <- as.numeric(mu.data$maturity_2_timing - mu.data$birth)

mu.data$age_maturity_mean <- rowMeans(mu.data[,c("age_maturity1", "age_maturity2")], na.rm = TRUE)

# Convert mean maturation time to rate of maturation (mu)
mu.data$mu <- 1/mu.data$age_maturity_mean

# Collect and clean all replicates with measured 'mu', for statistical analysis on groups
mu.reps <- as.data.frame(cbind(mu.data[,c(1:3,21)]))
# NB: "Source_pop" denotes sourced stock pop container from pre-exp rearing, not used in analysis (just used for extra measure of replicability)

# Subset replicate populations with non-NaN values
mu.reps <- mu.reps[complete.cases(mu.reps),]

# calculate summary statistics for all Source+Treatment+Rep permutations
mu.summary <- summarySE(mu.reps, measurevar = "mu", groupvars = c("Treatment", "rep"), na.rm = T)