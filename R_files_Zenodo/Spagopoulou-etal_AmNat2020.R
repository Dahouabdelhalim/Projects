#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# Clear working space
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

rm(list=ls())


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# Packages
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

library(lme4)
library(car)
library(DHARMa)


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# UPLOAD Data
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

data_P<-read.csv("Spagopoulou-etal_AmNat2020.csv", header=TRUE)


#Order rows according to Pair
data_P<-data_P[ order(data_P[,2]), ]


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# SPERM ANALYSIS
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

#~#~#~#~#~#~#~#~#~#~#~# Data Subsets #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#For Sperm Number
data_sperm <- subset(data_P, Trial=="A")

# For Sperm Velocity we need to use only the individuals that have more that 10 sperm cells recorded
data_sperm_tracks <- subset(data_sperm, Sperm_cells > 10)



#~#~#~#~#~#~#~#~#~#~#~# MODELS #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

##############____ Sperm Number ____###############

## Model with Interaction term
sn_int <- lm((sqrt(data_sperm$Sperm_number)) ~ Trt *scSL, data = data_sperm,  contrasts = list(Trt = contr.sum))

# Plot model
par(mfrow=c(2,2))
plot(sn_int)

# Summary statistics
summary(sn_int)
Anova(sn_int, type=3)



## Additive Model without interaction
sn_add <- lm((sqrt(data_sperm$Sperm_number)) ~ Trt + scSL, data = data_sperm, contrasts = list(Trt = contr.sum))

# Plot model
par(mfrow=c(2,2))
plot(sn_add)

# Summary statistics
summary(sn_add)
Anova(sn_add, type=3)

##############____ Sperm Velocity ____###############
#______________________ VCL _______________________#
## Model with Interaction term
svcl_int <- lm(VCL_average ~ Trt * scSL, data = data_sperm_tracks, contrasts = list(Trt = contr.sum))

# Plot model
par(mfrow=c(2,2))
plot(svcl_int)

# Summary statistics
summary(svcl_int)
Anova(svcl_int, type=3)



## Additive Model without interaction
svcl_add <- lm(VCL_average ~ Trt + scSL, data = data_sperm_tracks, contrasts = list(Trt = contr.sum))

# Plot model
par(mfrow=c(2,2))
plot(svcl_add)

# Summary statistics 
summary(svcl_add)
Anova(svcl_add, type=3)


#______________________ VAP _______________________#
## Model with Interaction term
svap_int <- lm(VAP_average ~ Trt * scSL, data = data_sperm_tracks, contrasts = list(Trt = contr.sum))

# Plot model
par(mfrow=c(2,2))
plot(svap_int)

# Summary statistics
summary(svap_int)
Anova(svap_int, type=3)



## Additive Model without interaction
svap_add <- lm(VAP_average ~ Trt + scSL, data = data_sperm_tracks, contrasts = list(Trt = contr.sum))

# Plot model
par(mfrow=c(2,2))
plot(svap_add)

# Summary statistics
summary(svap_add)
Anova(svap_add, type=3)


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# BEHAVIOURAL ANALYSIS
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

#~#~#~#~#~#~#~#~#~#~#~# Data Subsets #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# We excluded two individuals that had a problem with their air bludder after the sperm assay
data_behav<-subset(data_P, Pair!="19")
data_behav<-subset(data_behav, Pair!="42")


### Separate the dataset in each Day of behavioural observations
data_behav_A <- subset(data_behav, Trial=="A")
#table(data_behav_A$Trial)

data_behav_C <- subset(data_behav, Trial=="C")
#table(data_behav_C$Trial)


##############____ Copulation Atempts ____###############
#____________________ Rate of CA _____________________#

# Remove zeros, to investigate the rate of CA in the males that did perform this behaviour
data_CArate_A<-subset(data_behav_A, !CA=="0")
data_CArate_C<-subset(data_behav_C, !CA=="0")

#____________________ Binomial CA ____________________#

# Binomial CA response of males (Yes/No)
data_CAbinom_A<-data_behav_A
data_CAbinom_A$CA[data_CAbinom_A$CA>0] <- 1

data_CAbinom_C<-data_behav_C
data_CAbinom_C$CA[data_CAbinom_C$CA>0] <- 1


##############____ Courtship Displays ____###############
#__________________ Rate of Displays __________________#

# Remove zeros, to investigate the rate of displays in the males that did perform this behaviour
data_DisplayRate_A<-subset(data_behav_A, !Display=="0")
data_DisplayRate_C<-subset(data_behav_C, !Display=="0")

#__________________ Binomial Diplays __________________#

# Binomial Diplay response of males (Yes/No)
data_DisplayBinom_A<-data_behav_A
data_DisplayBinom_A$Display[data_DisplayBinom_A$Display>0] <- 1

data_DisplayBinom_C<-data_behav_C
data_DisplayBinom_C$Display[data_DisplayBinom_C$Display>0] <- 1




#~#~#~#~#~#~#~#~#~#~#~# MODELS #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

##############____ Copulation Atempts ____###############
#____________________ Rate of CA _____________________#
#~~~~~~~~~~~~~~~~~~~~~~ Day 1 ~~~~~~~~~~~~~~~~~~~~~~~#
## Model with Interaction term
ca.rate_A_int <- glmer(CA~scSL * Trt   
                  + (1|Pair) +(1|Obs)
                  , family= poisson ,data=data_CArate_A, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_ca.rate_A_int <- simulateResiduals(fittedModel = ca.rate_A_int, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_ca.rate_A_int) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_ca.rate_A_int) 

# Summary statistics
summary(ca.rate_A_int)
Anova(ca.rate_A_int, type=3)




## Additive Model without interaction
ca.rate_A_add <- glmer(CA~scSL + Trt   
                   + (1|Pair) +(1|Obs)
                   , family= poisson ,data=data_CArate_A, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_ca.rate_A_add <- simulateResiduals(fittedModel = ca.rate_A_add, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_ca.rate_A_add) 

#Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_ca.rate_A_add) 

# Summary statistics
summary(ca.rate_A_add)
Anova(ca.rate_A_add, type=3)


#~~~~~~~~~~~~~~~~~~~~~~ Day 3 ~~~~~~~~~~~~~~~~~~~~~~~#
## Model with Interaction term
ca.rate_C_int <- glmer(CA~scSL * Trt       
                  + (1|Pair) +(1|Obs)
                  , family= poisson ,data=data_CArate_C, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_ca.rate_C_int <- simulateResiduals(fittedModel = ca.rate_C_int, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_ca.rate_C_int) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_ca.rate_C_int) 

# Summary statistics
summary(ca.rate_C_int)
Anova(ca.rate_C_int, type=3)



## Additive Model without interaction
ca.rate_C_add <- glmer(CA~scSL + Trt       
                       + (1|Pair) +(1|Obs)
                       , family= poisson ,data=data_CArate_C, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_ca.rate_C_add <- simulateResiduals(fittedModel = ca.rate_C_add, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_ca.rate_C_add) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_ca.rate_C_add) 

# Summary statistics
summary(ca.rate_C_add)
Anova(ca.rate_C_add, type=3)


#____________________ Binomial CA _____________________#
#~~~~~~~~~~~~~~~~~~~~~~ Day 1 ~~~~~~~~~~~~~~~~~~~~~~~#
## Model with Interaction term
ca.binom_A_int <- glm(CA~scSL * Trt       
               , family= binomial ,data=data_CAbinom_A, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_ca.binom_A_int <- simulateResiduals(fittedModel = ca.binom_A_int, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_ca.binom_A_int) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_ca.binom_A_int) 

# Summary statistics
summary(ca.binom_A_int)
Anova(ca.binom_A_int,  type=3)



## Additive Model without interaction
ca.binom_A_add <- glm(CA~scSL + Trt       
                      , family= binomial ,data=data_CAbinom_A, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_ca.binom_A_add <- simulateResiduals(fittedModel = ca.binom_A_add, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_ca.binom_A_add) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_ca.binom_A_add) 

# Summary statistics
summary(ca.binom_A_add)
Anova(ca.binom_A_add, type=3)


#~~~~~~~~~~~~~~~~~~~~~~ Day 3 ~~~~~~~~~~~~~~~~~~~~~~~#
## Model with Interaction term
ca.binom_C_int <- glm(CA~scSL * Trt       
               , family= binomial ,data=data_CAbinom_C, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_ca.binom_C_int <- simulateResiduals(fittedModel = ca.binom_C_int, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_ca.binom_C_int) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_ca.binom_C_int) 

# Summary statistics
summary(ca.binom_C_int)
Anova(ca.binom_C_int, type=3)




##############____ Courtship Displays ____###############
#__________________ Rate of Displays ___________________#
#~~~~~~~~~~~~~~~~~~~~~~ Day 1 ~~~~~~~~~~~~~~~~~~~~~~~#
## Model with Interaction term
display.rate_A_int <- glmer(Display~scSL * Trt   
                   + (1|Pair) +(1|Obs)
                   , family= poisson ,data=data_DisplayRate_A, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_display.rate_A_int <- simulateResiduals(fittedModel = display.rate_A_int, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_display.rate_A_int) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_display.rate_A_int) 

# Summary statistics
summary(display.rate_A_int)
Anova(display.rate_A_int, type=3)



## Additive Model without interaction
display.rate_A_add <- glmer(Display~scSL + Trt   
                            + (1|Pair) +(1|Obs)
                            , family= poisson ,data=data_DisplayRate_A, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_display.rate_A_add <- simulateResiduals(fittedModel = display.rate_A_add, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_display.rate_A_add) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_display.rate_A_add) 

# Summary statistics
summary(display.rate_A_add)
Anova(display.rate_A_add, type=3)


#~~~~~~~~~~~~~~~~~~~~~~ Day 3 ~~~~~~~~~~~~~~~~~~~~~~~#
## Model with Interaction term
display.rate_C_int <- glmer(Display~scSL * Trt   
                            + (1|Pair) +(1|Obs)
                            , family= poisson ,data=data_DisplayRate_C, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_display.rate_C_int <- simulateResiduals(fittedModel = display.rate_C_int, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_display.rate_C_int) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_display.rate_C_int) 

# Summary statistics
summary(display.rate_C_int)
Anova(display.rate_C_int, type=3)



## Additive Model without interaction
display.rate_C_add <- glmer(Display~scSL + Trt   
                            + (1|Pair) +(1|Obs)
                            , family= poisson ,data=data_DisplayRate_C, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_display.rate_C_add <- simulateResiduals(fittedModel = display.rate_C_add, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_display.rate_C_add) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_display.rate_C_add) 

# Summary statistics
summary(display.rate_C_add)
Anova(display.rate_C_add, type=3)



#__________________ Binomial Displays ___________________#
#~~~~~~~~~~~~~~~~~~~~~~ Day 1 ~~~~~~~~~~~~~~~~~~~~~~~#
## Model with Interaction term
display.binom_A_int <- glm(Display~scSL * Trt       
                , family= binomial ,data=data_DisplayBinom_A, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_display.binom_A_int <- simulateResiduals(fittedModel = display.binom_A_int, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_display.binom_A_int) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_display.binom_A_int) 

# Summary statistics
summary(display.binom_A_int)
Anova(display.binom_A_int, type=3)



## Additive Model without interaction
display.binom_A_add <- glm(Display~scSL + Trt       
                           , family= binomial ,data=data_DisplayBinom_A, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_display.binom_A_add <- simulateResiduals(fittedModel = display.binom_A_add, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_display.binom_A_add) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_display.binom_A_add) 

# Summary statistics
summary(display.binom_A_add)
Anova(display.binom_A_add, type=3)




#~~~~~~~~~~~~~~~~~~~~~~ Day 3 ~~~~~~~~~~~~~~~~~~~~~~~#
## Model with Interaction term
display.binom_C_int <- glm(Display~scSL * Trt       
                           , family= binomial ,data=data_DisplayBinom_C, contrasts = list(Trt=contr.sum))

# Simulate data (for high percision, more than 1000 - default is 250)
simulationOutput_display.binom_C_int <- simulateResiduals(fittedModel = display.binom_C_int, n=2000)

# Plot model
par(mfrow=c(1,2))
plotSimulatedResiduals(simulationOutput = simulationOutput_display.binom_C_int) 

# Test for goodness-of-fit
testUniformity(simulationOutput = simulationOutput_display.binom_C_int) 

# Summary statistics
summary(display.binom_C_int)
Anova(display.binom_C_int, type=3)






#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# PATERNITY ANALYSIS
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

#~#~#~#~#~#~#~#~#~#~#~# Data Subsets #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# We need one line per individual (e.g. Trial A)
data_patern <- subset(data_P, Trial=="A")

#Remove the five last males in the dataset which died after the sperm trials
data_patern <- data_patern[-c(157:161),]
tail(data_patern)

# We excluded two individuals that had a problem with their air bludder after the sperm assay
data_patern<-subset(data_patern, Pair!="19")
data_patern<-subset(data_patern, Pair!="42")

#Reorder the rows from smallest to highest Pair value and then Competition to Control treatment within pair
data_patern<-data_patern[ order(data_patern[,2], data_patern[,4]), ]
head(data_patern)

# Create variables we will need in the model
for (i in 1:nrow(data_patern)) {
  
  if ( i==1 ){
    data_patern$DeltaHet[i] <- data_patern$Het[i+1] - data_patern$Het[i]
    data_patern$DeltaSL[i] <- data_patern$SL[i+1]-data_patern$SL[i]
    data_patern$scDeltaSL[i] <- data_patern$scSL[i+1] - data_patern$scSL[i]
    data_patern$DeltaRelGono[i] <- data_patern$Relative_Gonopodium[i+1]-data_patern$Relative_Gonopodium[i]
    data_patern$scDeltaRelGono[i] <- data_patern$scGono[i+1] - data_patern$scGono[i]
    
  } else if (i==nrow(data_patern)) {
    data_patern$DeltaHet[i] <- data_patern$Het[i]-data_patern$Het[i-1]
    data_patern$DeltaSL[i] <- data_patern$SL[i]-data_patern$SL[i-1]
    data_patern$scDeltaSL[i] <- data_patern$scSL[i] - data_patern$scSL[i-1]
    data_patern$DeltaRelGono[i] <- data_patern$Relative_Gonopodium[i]-data_patern$Relative_Gonopodium[i-1]
    data_patern$scDeltaRelGono[i] <- data_patern$scGono[i] - data_patern$scGono[i-1]
    
  } else if (data_patern$Pair[i]==data_patern$Pair[i-1]) {
    data_patern$DeltaHet[i] <- data_patern$DeltaHet[i-1]
    data_patern$DeltaSL[i] <- data_patern$DeltaSL[i-1]
    data_patern$scDeltaSL[i] <- data_patern$scDeltaSL[i-1]
    data_patern$DeltaRelGono[i] <- data_patern$DeltaRelGono[i-1]
    data_patern$scDeltaRelGono[i] <- data_patern$scDeltaRelGono[i-1]
    
  } else {
    data_patern$DeltaHet[i] <- data_patern$Het[i+1]-data_patern$Het[i]
    data_patern$DeltaSL[i] <- data_patern$SL[i+1]-data_patern$SL[i]
    data_patern$scDeltaSL[i] <- data_patern$scSL[i+1] - data_patern$scSL[i]
    data_patern$DeltaRelGono[i] <- data_patern$Relative_Gonopodium[i+1]-data_patern$Relative_Gonopodium[i]
    data_patern$scDeltaRelGono[i] <- data_patern$scGono[i+1] - data_patern$scGono[i]
    
  }
  
}

#Scale new variables
data_patern<- transform(data_patern,
                    DeltaSL_post_sc=scale(DeltaSL),   
                    DeltaRelGono_post_sc=scale(DeltaRelGono),
                    DeltaHet_post_sc=scale(DeltaHet))    

#head(data_patern)
#tail(data_patern)
#str(data_patern)

# Fix dummy variable
data_patern <- transform(data_patern, Obs=factor(seq(nrow(data_patern))))
data_patern$obs_num <- as.numeric(data_patern$Obs)


# Subset the dataset only for Competition or only for Control males
data_patern_Comp<-subset(data_patern, Trt=="Competition")
data_patern_Control<-subset(data_patern, Trt=="Control")



#~#~#~#~#~#~#~#~#~#~#~# MODELS #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#############____ Number of Offspring ____##############


p_comp <- glm(cbind(Offspring,(Total-Offspring)) ~ 1 + 
                  DeltaHet_post_sc + DeltaRelGono_post_sc + DeltaSL_post_sc +
                  scFemale_SL
                , family="quasibinomial", data=data_patern_Comp)

summary(p_comp)
Anova(p_comp, type=3)


#### Exclude the 2 males (from pairs 16 & 55) that had high levrage (outliers in DeltaSL) and influenced our model.
data_patern_Comp <- data_patern_Comp[-16,]
data_patern_Comp <- data_patern_Comp[-49,]

# Re-run the model without these 2 males
p_comp_new <- glm(cbind(Offspring,(Total-Offspring)) ~ 1 + 
                   DeltaHet_post_sc + DeltaRelGono_post_sc + DeltaSL_post_sc +
                   scFemale_SL
                 , family="quasibinomial", data=data_patern_Comp)


summary(p_comp_new)
Anova(p_comp_new, type=3)







