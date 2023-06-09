#R code for the multilevel selection analyses on individual strength and network density in males the manuscript "Multilevel selection on social network traits: resource and sex effects on levels of selection in experimental populations of forked fungus beetles"
#This script measures multilevel selection gradients on individual strength and population density (connectedness) in male only networks. This script first runs models with the observed dataset. This script then performs a node-level permutation to create 2000 randomizations of the observed dataset. This script then generates a null distribution of F-values from models that use the randomized datasets. After generating a null distribution of F-values, this script then compares the F-values from the model that uses observed data and calculates p-values.
#Data are available at **add dryad link**

#set working directory to where you saved the data file
setwd("")

#load libraries
library(lme4)
library(car)
library(glmmTMB)
library(DHARMa)
library(effects)
library(ggeffects)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(emmeans)
library(doParallel)

#change R's defult to contrasts 
options(contrasts=c("contr.sum", "contr.poly")) 

#upload data
ObsSNBeetleset <- read.csv("ObsStrengthDensityMalesBeetleset.csv")

#subset for males only
ObsSNBeetleset %>%
  filter(Survey_Sex == "M") -> ObsSNBeetleset

#permute data####

##set number of permutations
permutations=2000

##shuffle elytra, number of scans seen, treatment, period, strength, and density. Shuffle without replacement. 

###create empty lists of permuted data
PermSNShuffled<-rep(list(list()), length(permutations))

####permute
for(i in 1:permutations) {
  PermSNShuffled[[i]] <- ObsSNBeetleset %>% mutate(Elytra=sample(Elytra)) %>% mutate(ScansSeen=sample(ScansSeen)) %>% mutate(Treatment=sample(Treatment)) %>% mutate(Period=sample(Period)) %>% mutate(alpha=sample(alpha)) %>% mutate(StandardizedDensity=sample(StandardizedDensity)) 
}

###name each permutation sequentially
names(PermSNShuffled) <- paste("perm", 1:length(PermSNShuffled), sep="_")

#standardize individual-level metrics####

##observed dataset
ObsSNBeetleset %>%
  group_by(Condo, Period) %>%
  mutate(StandardizedElytra = ((Elytra-(mean(Elytra)))/sd(Elytra))) %>%
  mutate(StandardizedScansSeen = ((ScansSeen-(mean(ScansSeen)))/sd(ScansSeen))) %>%
  mutate(StandardizedAlpha = ((alpha-(mean(alpha)))/sd(alpha))) -> ObsSNBeetleset

##permuted datasets
###first, ungroup
PermSNShuffled<-lapply(PermSNShuffled, function(x) {ungroup(x)})
###then, group within each condo, period, and sex
PermSNShuffled<-lapply(PermSNShuffled, function(x) {group_by(x, Condo, Period)})
###calculate standardized elytra
PermSNShuffled<-lapply(PermSNShuffled, function(x) {mutate(x, StandardizedElytra = ((Elytra-mean(Elytra)))/sd(Elytra))})
###calculate standardized scans seen
PermSNShuffled<-lapply(PermSNShuffled, function(x) {mutate(x, StandardizedScansSeen = ((ScansSeen-mean(ScansSeen)))/sd(ScansSeen))})
###calculate standardized strength
PermSNShuffled<-lapply(PermSNShuffled, function(x) {mutate(x, StandardizedAlpha = ((alpha-mean(alpha)))/sd(alpha))})

#standardize fitness metrics####

##observed dataset
ObsSNBeetleset %>%
  ungroup() %>%
  group_by(Period) %>%
  mutate(GlobalStandardizedGuards = ((UniqueGuards/(mean(UniqueGuards))))) -> ObsSNBeetleset

##permuted datasets
###first, ungroup
PermSNShuffled<-lapply(PermSNShuffled, function(x) {ungroup(x)})
###then, group within each period and sex
PermSNShuffled<-lapply(PermSNShuffled, function(x) {group_by(x, Period)})
###calculate global standardized guards
PermSNShuffled<-lapply(PermSNShuffled, function(x) {mutate(x, GlobalStandardizedGuards = ((UniqueGuards/mean(UniqueGuards))))})

#observed models####

##guards~strength+density####
obs.guards.strength.density.selection.global<-glmmTMB(GlobalStandardizedGuards~StandardizedAlpha*Treatment+StandardizedDensity*Treatment+StandardizedElytra+StandardizedScansSeen+Period+(1|Condo/ID), data=ObsSNBeetleset)
summary(obs.guards.strength.density.selection.global)
Anova(obs.guards.strength.density.selection.global, type=3)
plot(allEffects(obs.guards.strength.density.selection.global))

###test assumptions with DHARMa
obs.guards.strength.density.selection.global.resid<-simulateResiduals(obs.guards.strength.density.selection.global, n=250)
plot(obs.guards.strength.density.selection.global.resid)
testDispersion(obs.guards.strength.density.selection.global.resid)
testZeroInflation(obs.guards.strength.density.selection.global.resid)

###extra F-values from Anova
obs.mls.males.Anova<-Anova(obs.guards.strength.density.selection.global, type=3)
obs.mls.males.Anova %>%
  #extract F-values
  dplyr::select("Chisq") %>%
  #tranform dataframe (switch rows and columns)
  t() -> obs.mls.males.Anova.Chi

###calculate marginal means
emmeans(obs.guards.strength.density.selection.global, "Period")
emmeans(obs.guards.strength.density.selection.global, "Treatment")
emtrends(obs.guards.strength.density.selection.global, pairwise~Treatment, var="StandardizedAlpha")
emtrends(obs.guards.strength.density.selection.global, pairwise~Treatment, var="StandardizedDensity")

#permuted models####

##function that extracts the type 3 sum of squares (written by VAF)
AnovaT3=function(x){
  setTimeLimit(cpu = Inf,  transient = FALSE)
  return(tryCatch(Anova(x, type=3), error=function(e) NULL))
}

##guards~strength+density####

###function that runs the model
perm.mls.males = function(PermSNMales){
  perm.mls.males.mod<-glmmTMB(GlobalStandardizedGuards~StandardizedAlpha*Treatment+StandardizedDensity*Treatment+StandardizedElytra+StandardizedScansSeen+Period+(1|Condo/ID), data=PermSNMales)
  return(perm.mls.males.mod)
}

###run the model on each permutation
perm.mls.males.results<-lapply(PermSNMales, FUN=perm.mls.males)

###extract type 3 sum of squares for each permutation
perm.mls.males.Anova<-lapply(perm.mls.males.results, FUN=AnovaT3)

###extract F-values (chi square values) from Anova
####first, clean up Anova output
perm.mls.males.Anova<-lapply(perm.mls.males.Anova, tibble::rownames_to_column)
####next, extract rownames for later
variable_names.male<-perm.mls.males.Anova[[1]][,1]
####then, extract F-values
perm.mls.males.Anova.Chi<-foreach(i = 1:length(perm.mls.males.Anova)) %dopar% {
  perm.mls.males.Anova[[i]] %>% select("Chisq")
}
####save F-values in a large dataframe
perm.mls.males.Anova.Chi.df<-do.call(cbind, perm.mls.males.Anova.Chi)
####add back in rownames
row.names(perm.mls.males.Anova.Chi.df)<-variable_names.male
####transform the dataframe (swap rows and columns)
trans.perm.mls.males.Anova.Chi.df<-as.data.frame(t(perm.mls.males.Anova.Chi.df))

#calculate p-values####
# FUNCTION: CWW (2015) Calculates one- and two-tailed p-values from a permuted distribution  
# Author Corlett Wolfe Wood

# One-tailed: The proportion of permuted values that are GREATER THAN the observed value
# Two-tailed: The proportion of permuted values that are MORE EXTREME than the observed

# NB: Calculation of two-tailed p-values on asymmetric null distributions
# On symmetric distributions, one-tailed p-values can be doubled to obtained two-tailed ones
# This doesn't work with asymmetric distributions because of skewness
# --> SOLUTION: This function below calculates two-tailed p-values as the proportion of
#     permuted.values that have a probability <= the probability of the observed value
# --> SOURCE: Pratt, J. W. and J. D. Gibbons. 1981. Concepts in Nonparametric Theory. 
#     Springer-Verlag, New York, NY pp. 29-32

# NB: p-values should NEVER be zero in a permutation test
# --> the observed data is one possible permutation, so at least one value
#     in the permuted data can be equal to the observed data (P>0)
# --> SOURCE: North, B. V., D. Curtis, and P. C. Sham. 2002. A note on the calculation of empirical 
#     p values from Monte Carlo procedures. Am. J. Hum. Genet. 71:439-441.
# --> SOURCE: Davison, A. C. and D. V. Hinkley. 1997. Bootstrap methods and their application. 
#     Cambridge University Press, Cambridge, United Kingdom
# --> SOLUTION: add 1 to the numerator and the denominator

pvals = function(observed.value, permuted.values){  
  # ONE-TAILED
  # Proportion of permuted values >= (if in upper tail) OR <= (if in lower tail) the observed value
  k.greater = length(which(permuted.values>=observed.value))
  k.less = length(which(permuted.values<=observed.value))
  n = length(permuted.values)
  # One-tailed p-value
  pval.onetail = (k.greater+1)/(n+1) # See NB above for explanation of +1  
  
  # TWO-TAILED
  # Proportion of permuted values with PROBABILITIES <= the observed value
  df <- data.frame(permuted.values=permuted.values, prob=NA)
  names(df) <- c("permuted.values", "prob")
  for(p in 1:nrow(df)){
    # Calculate probability of each permuted value
    # Proportion of permuted values >= (if in upper tail) OR <= (if in lower tail) each permuted value
    k.greater = length(which(df$permuted.values>=df$permuted.values[p]))
    k.less = length(which(df$permuted.values<=df$permuted.values[p]))
    n = nrow(df)
    df$prob[p] = (min(k.greater, k.less)+1)/(n+1) # See NB above for explanation of +1
  }
  
  # Proportion of permuted.values with a have a probability <= the probability of the observed value
  prob.less <- length(which(df$prob<=pval.onetail))
  pval.twotail <- (prob.less+1)/(n+1) # See NB above for explanation of +1
  
  return(c(one.tailed=pval.onetail, 
           two.tailed=pval.twotail))
}

##guards~strength+density####

###create a new df to hold pvalues
perm.mls.males.pvalues<-data.frame(one.tailed=numeric(), two.tailed=numeric())

###use pval function to calculate pvalues
for(i in 1:length(trans.perm.mls.males.Anova.Chi.df)){
  perm.mls.males.pvalues[i,]<-pvals(obs.mls.males.Anova.Chi[1,i], trans.perm.mls.males.Anova.Chi.df[1:nrow(trans.perm.mls.males.Anova.Chi.df),i])
}

###add a column with the name of the variables
perm.mls.males.pvalues$var.names<-variable_names.male