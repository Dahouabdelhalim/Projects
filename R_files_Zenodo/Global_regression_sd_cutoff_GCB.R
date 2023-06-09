###############################
##Linear regression with only top models
##############################

##Import SIF values for 2015

SIF.models.2015 <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/SIF.models.2015a.csv")

SIF <- SIF.models.2015





####Get covariates####

##Functional dissimilarity##

Full.trait.list.bi.1.kg <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/Full trait list bi 1 kg.csv")
Species.list <- Full.trait.list.bi.1.kg
Species.list <- filter(Species.list,Species.list$Species %in% SIF$Spp1 | Species.list$Species %in% SIF$Spp2)
rownames(Species.list) <- Species.list$Species
Species.list$Species <- NULL
Species.list$Body.Mass <- log(Species.list$Body.Mass)
Species.list$Litter.Size <- log(Species.list$Litter.Size)

#Calculate functional dissimilarity
library(cluster)
k <- daisy(Species.list, metric = "gower",weights = c(15,3,3,3,3,3,15,15,5,5,5,15))
k <- as.matrix(k)
rownames(k) <- rownames(Species.list)
colnames(k) <- rownames(Species.list)
library(reshape2)


##Functional Richness##
#Get protected area communities

SIF.year <- SIF

Func.comm <-  matrix(NA, nrow = length(unique(SIF.year$Site.Code)), ncol = nrow(Species.list), dimnames = list(unique(SIF.year$Site.Code),rownames(Species.list)))  
for (i in 1:length(unique(SIF.year$Site.Code))){
  Code <- unique(SIF.year$Site.Code)[i]
  Community.Traits <- filter(SIF.year,SIF.year$Site.Code == Code)
  Comm.Species <- c(as.character(Community.Traits$Spp1),as.character(Community.Traits$Spp2))
  Community.list <- as.character(unique(Comm.Species))
  Community.pres <- matrix(NA, nrow = 1, ncol = nrow(Species.list))
  for (j in 1:nrow(Species.list)) {
    value <- ifelse(as.character(rownames(Species.list)[j]) %in% Community.list, 1, 0)
    Community.pres[,j] <- value
  }
  Func.comm[i,] <- Community.pres
  
  
  rownames(Func.comm) <- paste(rownames(Func.comm))
}




#Calculate functional richness
library(FD)
library(dplyr)
library(tidyr)
x <- dbFD(Species.list, Func.comm, w = c(15,3,3,3,3,3,15,15,5,5,5,15))
Site.FRich <- as.data.frame(x$FRic)
colnames(Site.FRich) <- "Func.Rich"
Site.FRich$Site.Code <- rownames(Site.FRich)
rownames(Site.FRich) <- NULL

##Merge functional dissimilarity data and functional richness with SIF models
diss <- reshape2::melt(as.matrix(k), varnames = c("Species 1", "Species 2"))
colnames(diss) <- c("Spp1", "Spp2", "Func.Diss")
SIF <- left_join(SIF,diss, by =c("Spp1","Spp2"))
SIF <- left_join(SIF,Site.FRich, by = c("Site.Code"))


##Site level covariates##
Site.covs.2015 <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/Site.covs.2015.csv")
SIF1 <- left_join(SIF,Site.covs.2015, by = "Site.Code")


SIF1$Spp.Site <- paste(SIF1$Spp1,SIF1$Spp2,SIF1$Site.Code, sep = "_")
SIF2 <- SIF1[,c(1,4,5,9:14)]
SIF2$logSIF <- log(SIF2$Mean.SIF)

##Set threshold for SD
SIF2 <- filter(SIF2, SIF2$sd.SIF < 0.2)
# The R code below uses the R2jags package to call JAGS
library(R2jags)

#The data file includes information for all species pairs, including the numerically-valued #covariates for all coarse- and fine-scale trait combinations, as well as the mean and variance of SIF #estimates.
data <- SIF2

data$Func.Diss <- scale(data$Func.Diss)
data$NDVI.2015 <- scale(data$NDVI.2015)
data$Hum.Dens <- scale(data$Hum.Dens)
data$Hab.div.2015 <- scale(data$Hab.div.2015)
data$Func.Rich <- scale(data$Func.Rich)

#Specify the number of chains (nchain), number of iterations (niter), burn-in period (nburn) and #thinning rate (nthin)
niter <- 20000
nthin <- 1
nburn <- 5000
nchain <- 3

#Specify the parameters to monitor
parameters <- c('alpha', 'bFuncDiss','bNDVI','bHumDens',
                'bHabHet','bFRich','bFDNDVI','bFDHum','bFDHab',
                'mu')

#Model text file for weighted regression
sink("model.weightedregression.txt")
cat(" model{
#Likelihood    
for (i in 1:N){
  y[i] ~ dnorm(mu[i], tau.sif.var[i])
  tau.sif.var[i] <- 1/(sif.var[i])
  mu[i] <- alpha + 
    bFuncDiss*X[i,1]+
    bFRich*X[i,2]+
    bNDVI*X[i,3]+
    bHumDens*X[i,4]+
    bHabHet*X[i,5]+
    bFDNDVI*X[i,1]*X[i,3]+
    bFDHum*X[i,1]*X[i,4]+
    bFDHab*X[i,1]*X[i,5]+
    eps.pairs[i]
  eps.pairs[i]~dnorm(0,tau.reg)
}

# Priors:
 alpha ~ dnorm(0, 0.01)
 bFuncDiss ~ dnorm(0,0.01)
 bNDVI ~ dnorm(0,0.01)
 bHumDens ~ dnorm(0,0.01)
 bHabHet ~ dnorm(0,0.01)
 bFRich ~ dnorm(0,0.01)
 bFDNDVI ~ dnorm(0,0.01)
 bFDHum ~ dnorm(0,0.01)
 bFDHab ~ dnorm(0,0.01)
 sigma.reg ~ dunif(0,10)
 tau.reg <- 1 / (sigma.reg * sigma.reg)
 


}", fill = TRUE)
sink()

##Specify X as predictor variables
X <- data[,c(4:8)]

#Load the data necessary for JAGS
jags.data <- list(y = data$logSIF, sif.var = (data$sd.SIF*data$sd.SIF), N = nrow(data), X = X)

#Run the weighted regression model and save the output as ‘out’
weighted.out <- jags(data = jags.data, parameters.to.save = parameters, n.iter = niter, n.burnin = nburn, n.chains = nchain, n.thin = nthin, model.file = "model.weightedregression.txt")

head(weighted.out$BUGSoutput$summary, n = 10)


###Plot regression results###
library(bayesplot)
library(rstan)
library(rstantools)
library(brms)
library(ggplot2)
library(ggdist)
library(tidybayes)
library(tidyverse)

#plot(data$Func.Diss,data$logSIF)
x <- weighted.out$BUGSoutput %>%
  spread_draws(alpha,bFuncDiss,bNDVI,bHumDens,
               bHabHet,bFRich,bFDNDVI,bFDHum,bFDHab)%>%
  gather(key = "Effect", value = "Est" ) #%>%
x <- x[-c(1:180000),]
x$Est <- as.numeric(x$Est)
x$Effect <- factor(x$Effect, levels = c('bFuncDiss','bNDVI','bHumDens',
                                        'bHabHet','bFRich','bFDNDVI','bFDHum','bFDHab'))
library(ggplot2)
ggplot(data = x, aes(y = Effect, x = Est, fill = Effect)) +
  stat_halfeye(.width = c(.95, .5))+
  geom_vline(xintercept = 0)+
  labs(title="Species Interaction Factor 2015- Top Model Subset",
       x ="Beta effect size", y = "Predictor variables")+
  # scale_y_discrete(labels = c("alpha" = "alpha","Func.Diss" = "Func.Diss",
  #     "NDVI"="NDVI", "Hum.dens"="Hum.dens",
  #      "Hab.div"="Hab.div","Hab.het"="Hab.het", "FD.NDVI"="FD.NDVI",
  #      "FD.Hum"="FD.Hum","FD.Hab"="FD.Hab"))+
  scale_fill_manual(values = c("grey85","grey85","grey85","grey85","grey85","grey85","red","grey85"))+
  theme_classic()+
  theme(legend.position = "none")

