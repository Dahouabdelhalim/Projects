# Determine the effect of hatchery line on each phenotypic trait
# using fish from 2010 and linear or generalized linear models. We are only testing 
# 2010 because that is the most recent year for which we have large sample sizes
# for all traits. We are not testing earlier years because, if divergence has occurred,
# it would be observed in the latter generations but not the earlier generations.

# Read in phenotypic data collected at Roza for 2010 
# (same data as in Table S2a from Supporting Information)
indivs_2010 <- read.csv("Phenotype_Data_2010.csv", header=TRUE) 

# Each row provides data from a single fish

# Test differences between the hatchery lines using linear or generalized linear models 

# make sure the variables are factors before running the models
indivs_2010$Line <- as.factor(indivs_2010$Line)
indivs_2010$sex <- as.factor(indivs_2010$sex) 

# Now exclude any row with a missing value for age, hatchery line, or sex
indivs_2010_no_missing <- indivs_2010[which(indivs_2010$age!="NA" & indivs_2010$Line!="NA" & indivs_2010$sex!="NA"),]

######################################################################################################################################################
# Quantify the effect of hatchery line on age at maturity

# Fit a GLM for age at maturity 
glm_age <- glm(age~Line*sex,family=poisson(link='log'), data=indivs_2010_no_missing) 
summary(glm_age) #remember that the coefficients must be exponentiated because we used the log link function


########################################################################################################################################################
########################################################################################################################################################

# For the other traits, use only age 4 fish to remove the 
# confounding effect of age, since ~90% of the population is age 4
indivs_2010_no_missing_age4 <- indivs_2010_no_missing[which(indivs_2010_no_missing$age==4),] #data set is reduced to 1184 fish

#Examine data set and sample sizes
length(which(indivs_2010_no_missing_age4$Line=="SEG")) #200 fish in SEG line overall
length(which(indivs_2010_no_missing_age4$Line=="INT")) #984 fish in INT line overall

########################################################################################################################################################
# Quantify the effect of hatchery line on return time

# Fit a linear model for return time
lme_return_2010_age4 <- lm(ReturnDay~Line*sex,data=indivs_2010_no_missing_age4) 
summary(lme_return_2010_age4) 

######################################################################################################################################################
# Quantify the effect of hatchery line on fork length

# Fit a linear model for fork length
lme_length_2010_age4 <- lm(forklgth~Line*sex,data=indivs_2010_no_missing_age4) 
summary(lme_length_2010_age4) 

######################################################################################################################################################
# Quantify the effect of hatchery line on weight

# Fit a linear model for weight but include an allometric effect (essentially control for length; e.g. Thorson 2015 Mar Ecol Prog Ser 526:101-112)
lme_weight_2010_age4_log <- lm(log(weight)~Line*sex + log(forklgth),data=indivs_2010_no_missing_age4) 
summary(lme_weight_2010_age4_log) # The coefficients for the intercept, Line, sex, and Line*sex need to be exponentiated

######################################################################################################################################################
######################################################################################################################################################

# Effect of hatchery line on spawn (maturation) time and daily growth coefficient 
# can only be quantified using my sample data for these two traits.
# Again use only 2010 data (from Table S2b in Supporting Information)

hatchery_lines<-read.csv("CESRF_phenotype_data_2010.csv",header=TRUE,na.strings=NA)
# Note that Origin=Line

# make covariates factors
hatchery_lines$sex <- factor(hatchery_lines$sex) 
hatchery_lines$Origin <- factor(hatchery_lines$Origin)

# subset data to only include age 4 fish
hatchery_lines_age4 <- hatchery_lines[which(hatchery_lines$age==4),] #reduced to 112 fish

###############################################################################################################
# Quantify the effect of hatchery line on spawn time for age 4 fish only

# Fit a linear model for spawn time
lme_spawn_2010_age4 <- lm(spawn_day_of_year~Origin*sex,data=hatchery_lines_age4) 
summary(lme_spawn_2010_age4) 

############################################################################################################################
# # Quantify the effect of hatchery line on DGC for age 4 fish only

#Remove fish that have missing values for DGC
hatchery_lines_age4_DGC <- hatchery_lines_age4[which(hatchery_lines_age4$Roza_CESRF_DGC!="NA"),]

# Fit a linear model for DGC
lme_DGC_2010_age4 <- lm(Roza_CESRF_DGC~Origin*sex,data=hatchery_lines_age4_DGC) 
summary(lme_DGC_2010_age4) 
