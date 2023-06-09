#load packages (assuming they are downloaded already)
library(car) #to get pvalues from model
library(Rmisc) #to get basic stats
library(FSA) #to get basic stats
library(emmeans) #pairwise analysis
library(ggplot2)#build graphs
library(Hmisc)
library(splines)
library(sandwich)
library(effects)
library(RcmdrMisc)

# for frequency/count data, using ZIP models
library(pscl)

#for duration/continuous data, using gamma dist 
#The Tweedie package allows you to fit a glm with any power function and any power link. 
#The var.power refers to exponent of the glm variance function, so that var.power=0 specifies a normal family, var.power=1 means Poisson family, var.power=2 means gamma family, var.power=3 means inverse Gaussian family and so on. 
#link.power=0 specifies a log-link. The link is specified in terms of Box-Cox transformation powers, so link.power=1 is the identity link and link.power=0 means log.
library(statmod)
library(tweedie)

#read in data.xlsx from directory
Dataset <- readXL("data.xlsx", rownames=FALSE, header=TRUE, na="", sheet="Compare ws and ds live.4grps", stringsAsFactors=TRUE)

#within same treatment (yellow/ws), 
#subset according to what variable you want to compare (between seasons, within treatment in one season... there should be different variable names for all possible comparisons in dataset?)
#we will be using total.dur for comparing duration between variables, and total.freq for comparing frequency between variables.

summary(model<-zeroinfl(total.freq~season, data=Dataset, dist="negbin")) ## negbin for overdispersed data, poisson for not overdispersed. All count/frequency data in this dataset is overdispersed so used negbin
summary(model<-glm(total.dur~season, family=tweedie(link.power=0, var.power=2),maxit=100,data=Dataset)) # var.power=2 == gamma dist

# check residuals/ overdispersion #
E2 <- resid(model, type = "pearson")
N  <- nrow(Dataset)
p  <- length(coef(model))   
sum(E2^2) / (N - p)


#courtship latency and mating duration
summary(model<-glm(latency.mated~treatment.mated*season.mated, family = gaussian, data=Dataset))
summary(model<-glm(duration.mated~treatment.mated*season.mated, family = gaussian, data=Dataset))

summary(model<-glm(latency.mated~type.mated, family = gaussian, data=Dataset))
summary(model<-glm(duration.mated~season.mated, family = gaussian, data=Dataset))

