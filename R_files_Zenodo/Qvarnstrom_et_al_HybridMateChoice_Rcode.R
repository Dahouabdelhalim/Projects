#### Patterns of heterospecific matings Ficedula flycatchers #### 

require(lme4)

#install.packages("lmerTest")
require(lmerTest)

### combined data from Oland and Gotland ###

# dataset from the female perspective, one entry per breeding female, proportions calculated per year and area

matingdata<-read.csv("individual_level_data_2021.csv", header = TRUE)

str(matingdata)

nrow(matingdata)


#### Make new variable with proportion of pieds in the population ####

propPF <- (matingdata$PFFn+matingdata$PFMn)/(matingdata$PFFn+matingdata$PFMn + matingdata$CFFn + matingdata$CFMn) 

hist(propPF)


#####################################################################################################################
###### 	PART 1: INDIVIDUAL PROBABILITY TO MATE WITH CONSPECIFIC IN RELATION TO PROPORTION OF PIED FLYCATCHERS ######
#####################################################################################################################



# binomial GLMM to examine patterns in probability to mate with a heterospecific in relation to the percentage of pied flycatchers

#each species separately
# PF: pfms = pied flycatcher mating status, ppf = proportion females, pfarea: random effect of area
pfms <- matingdata$matingstatus[which(matingdata$spf2=="PF")]
ppf <- propPF[which(matingdata$spf2=="PF")]
pfarea <- matingdata$area[which(matingdata$spf2=="PF")]


summary(glmer(pfms~ppf + (1|pfarea), family=binomial, data= matingdata))


#CF
cfms <- matingdata$matingstatus[which(matingdata$spf2=="CF")]
proppf2 <- propPF[which(matingdata$spf2=="CF")]
cfarea <- matingdata$area[which(matingdata$spf2=="CF")]
ringf_cf <- matingdata$id_f[which(matingdata$spf2=="CF")] 

#model fit

summary(glmer(cfms~proppf2  + (1|cfarea), family=binomial, data= matingdata))


#####################################################################################################################
# PART 2. PROPORTION HETEROSPECIFIC MATINGS IN RELATION TO PF ABUNDANCE												#
#####################################################################################################################

# Population level data
matingPopdata <-read.csv("population_level_data_2021.csv", header = TRUE)

str(matingPopdata)

nrow(matingPopdata)


#convert percentage PF to proportion

matingPopdata$pPF <- (matingPopdata$PFFn+ matingPopdata$PFMn)/(matingPopdata$PFFn+ matingPopdata$PFMn + matingPopdata$CFFn + matingPopdata$CFMn) 


#Examine frequency of pied females mated to collared males among all pairings

pfHetpairs <- cbind(matingPopdata$nrPFHetParis,(matingPopdata$nrCFHetPairs+matingPopdata$homPairs))

pfhetmixmod1 <- glmer(pfHetpairs~pPF + I(pPF^2) + (1|area), family=binomial, data=matingPopdata, )
pfhetmixmod2 <- glmer(pfHetpairs~pPF + (1|area), family=binomial, data=matingPopdata)

anova(pfhetmixmod1,pfhetmixmod2)

summary(pfhetmixmod1)

#Examine frequency of collared females mated to pied males among all pairings

cfHetpairs <- cbind(matingPopdata$nrCFHetPairs,(matingPopdata$nrPFHetParis+matingPopdata$homPairs))


cfhetmixmod1 <- glmer(cfHetpairs~pPF + I(pPF^2) + (1|area), family=binomial, data=matingPopdata)
cfhetmixmod2 <- glmer(cfHetpairs~pPF + (1|area), family=binomial, data=matingPopdata)

anova(cfhetmixmod1,cfhetmixmod2)

summary(cfhetmixmod1)

#TOTAL hybridization frequency (both species combined)
# Create variable that is number of heterospecific matings and number of homospecific matings
hetMating <- cbind(matingPopdata$hetPairs,matingPopdata$homPairs)

#Model

mixmod1 <- glmer(hetMating~pPF + I(pPF^2) + (1|area), family=binomial, data=matingPopdata)

mixmod2 <- glmer(hetMating~pPF + (1|area), family=binomial, data=matingPopdata)

anova(mixmod1,mixmod2,test="Chisq")

summary(mixmod1)


##############################################################################################################
## 		PART 3 - ASYMMETRY PART -  PROPORTION PF/CF OF HETEROSPECIFIC MATINGS IN RELALTION TO PF ABUNDANCE 	##
##############################################################################################################

# subset matingPopdata to only include areas with some heterospecific pairs

asymmdata <- subset(matingPopdata, hetPairs>0)

nrow(asymmdata)

head (asymmdata)


#proportion pied flycatchers in mixed pairs in relation to pied flycatcher abundance
asymmpropPF <- (asymmdata$PFFn+ asymmdata$PFMn)/(asymmdata$PFFn+ asymmdata$PFMn + asymmdata$CFFn + asymmdata$CFMn)

#Number of piedflycatcher females in heterospecific matings
PFhyb <- cbind(asymmdata$nrPFHetParis, asymmdata$nrCFHetPairs)

#mixed model
PFmod <- glmer(PFhyb~ asymmpropPF + I(asymmpropPF^2) + (1|area), family=binomial, data= asymmdata)

summary(PFmod)

PFmod2 <- glmer(PFhyb~ asymmpropPF + (1|area), family=binomial, data= asymmdata)

summary(PFmod2)

PFmod3 <- glmer(PFhyb~ 1 + (1|area), family=binomial, data= asymmdata)


anova(PFmod,PFmod2, PFmod3)



#####################################################################################################################
# PART 4. ANALYSIS OF OVERALL HETEROSPECIFIC MATING PATTERNS																#
#####################################################################################################################

# data merging

hetTable <- aggregate(matingdata$matingstatus,by=list(spf=matingdata$spf,spm=matingdata$spm),length) # make a summary table over all pairwise species comparisons


# count nr of total pairs

hetTable[1,3]+hetTable[4,3]+hetTable[2,3]+hetTable[3,3]


# Proportion heterospecific pairs
(hetTable[2,3]+hetTable[3,3] )/(hetTable[1,3]+hetTable[4,3]+hetTable[2,3]+hetTable[3,3])


### Analyses of asymmetry bias ####

## binomial test for asymmetry in hybridization rates ##

totHyb <- hetTable[2,3] + hetTable[3,3] # total number of heterospecific pairings, n = 353


prop.test(c(226,127),c(353,353)) # binomial test to see if the proportion of females involved is asymmetric





