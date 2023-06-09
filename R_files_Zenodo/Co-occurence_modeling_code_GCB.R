###Calculate species interaction factor for all species pairs at all TEAM
###sites in 2015 and run Bayesian regression on the results
library(dplyr)
library(tidyr)

##Camera trap data available for download on Wildlife Insights

##Import TEAM occurrence data
ternary_matrix_2017.05.03 <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/ternary_matrix_2017-05-03.csv")
q <- ternary_matrix_2017.05.03
q$Spp.name <- paste(q$genus ,q$species)

##Format and filter
df1 <- q
df1$site_name <- as.character(df1$site_name)
#df1[df1 == "Bukit Barisan"] <- "BBS"
df1[df1 == "Barro Colorado Nature Monument - Soberania National Park"] <- "BCI"
df1[df1 == "Bwindi Impenetrable Forest" ] <- "BIF"
df1[df1 == "Caxiuanã" ] <- "CAX"
df1[df1 == "Cocha Cashu - Manu National Park" ] <- "COU"
df1[df1 == "Central Suriname Nature Reserve"] <- "CSN"
df1[df1 == "Korup National Park" ] <- "KRP"
df1[df1 == "Nam Kading"] <- "NAK"
df1[df1 == "Nouabalé Ndoki"] <- "NNN"
df1[df1 == "Pasoh Forest Reserve"] <- "PSH"
df1[df1 == "Ranomafana"] <- "RNF"
df1[df1 == "Udzungwa"] <- "UDZ"
df1[df1 == "Volcán Barva"] <- "VB"
df1[df1 == "Yanachaga Chimillén National Park"] <- "YAN"
df1[df1 == "Yasuni"] <- "YAS"
q <- df1
q <- filter(q,q$site_name != "Manaus" & q$site_name != "")
q <- filter(q,q$site_name != "Bukit Barisan" & q$site_name != "")
##Select year of interest
q <- filter(q, q$year == 2015)

##Filter to species over 1 kg and produce file for species pair detection
##in the selected year

unique(q$site_name)
for (j in unique(q$site_name)) {
  
  r <- filter(q,q$site_name == j)
  s <- unique(r$Spp.name)
  s <- as.data.frame(s)
  #Import trait list and species list
  Full.trait.list.bi.1.kg <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/Full trait list bi 1 kg.csv")
  Community.traits <- read.csv(paste0("~/Documents/TEAM_co-occurence/TEAM_co-occurence/Site Trait and Occ data/",j,"/",j," Community traits.csv", sep = ""))
  Species.list <- filter(Full.trait.list.bi.1.kg,Full.trait.list.bi.1.kg$Species %in% Community.traits$Species)
  s <- filter(s,s %in% Species.list$Species)
  rr <- r
  
  rr <- filter(rr,rr$detected >= 0) ##Take only time points where camera is active
  t <- rr %>%
    group_by(Spp.name) %>%
    summarise(sum(detected))
  t <- filter(t,t$`sum(detected)`> 3)
  s <- filter(s,s %in% t$Spp.name)
  Spp.list <- combn(s$s,2)
  m <-  matrix(NA, nrow = 0, ncol = 7)
  colnames(m) <- c("CT_ID", "Site", "Year","Spp1","Spp1_Det", "Spp2","Spp2_Det")
  for (i in 1:ncol(Spp.list)){
    species1 <- as.character(Spp.list[1,i])
    species2 <- as.character(Spp.list[2,i])
    x <- filter(rr,rr$Spp.name== paste(species1))
    x$active <- 1
    y <- aggregate(detected ~Spp.name + ct_id + site_name, data = x, FUN = "sum")
    z <- aggregate(active ~Spp.name + ct_id + site_name, data = x, FUN = "sum")
    yz <- merge(y,z,by =c("ct_id","site_name","Spp.name"))
    colnames(yz) <- c("CT_ID","Site","Spp1","Spp1_Det","Samp.pers")
    #colnames(y) <- c("Spp1","CT_ID","Site","Spp1_Det")
    
    xx <- filter(rr,rr$Spp.name== paste(species2))
    xx$active <- 1
    yy <- aggregate(detected ~Spp.name + ct_id + site_name, data = xx, FUN = "sum")
    zz <- aggregate(active~Spp.name + ct_id + site_name, data = xx, FUN = "sum")
    yyzz <- merge(yy,zz,by =c("ct_id","site_name","Spp.name"))
    colnames(yyzz) <- c("CT_ID","Site","Spp2","Spp2_Det","Samp.pers")
    #colnames(yy) <- c("Spp2","CT_ID","Site","Spp2_Det")
    
    dat <- merge(yz,yyzz, by = c("CT_ID","Site","Samp.pers"))
    #dat <- merge(y,yy, by = c("CT_ID","Site"))
    
    m <- rbind(m,dat)
  }
  
  write.csv(m,file = paste(j,"pairs.csv",sep = "."))
}

##Import species pair detection files
BCI.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/BCI.pairs.csv")
#BBS.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/BBS.pairs.csv")
BIF.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/BIF.pairs.csv")
CAX.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/CAX.pairs.csv")
CSN.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/CSN.pairs.csv")
COU.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/COU.pairs.csv")
KRP.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/KRP.pairs.csv")
NAK.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/NAK.pairs.csv")
NNN.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/NNN.pairs.csv")
PSH.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/PSH.pairs.csv")
RNF.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/RNF.pairs.csv")
UDZ.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/UDZ.pairs.csv")
VB.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/VB.pairs.csv")
YAN.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/YAN.pairs.csv")
YAS.pairs <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/YAS.pairs.csv")


##Bind species pair detections for all sites

TEAM.data <- rbind(BCI.pairs,BIF.pairs,CAX.pairs,COU.pairs,
                   CSN.pairs,KRP.pairs,NNN.pairs,PSH.pairs,
                   RNF.pairs,UDZ.pairs,VB.pairs,YAN.pairs,YAS.pairs)

write.csv(TEAM.data,file="TEAM.det.data.csv")


###SINGLE YEAR MODEL
# The R code below uses the R2jags package to call JAGS
library(R2jags)
library(dplyr)

##Import species pair detection data
data <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/TEAM.det.data.csv")

#Specify the number of chains (nchain), number of iterations (niter), burn-in period (nburn) and #thinning rate (nthin)
niter <- 50000
nthin <- 1
nburn <- 10000
nchain <- 3

#Specify the parameters to monitor
parameters <- c('SIF')

#Specify the index used to loop through all species pairs
index <- unique(data[ , c("Site", "Spp1", "Spp2")])

#Create a matrix ‘result’ to store the results for all species pairs
result <- matrix(NA, nrow(index), 5)

#Create matrix to store full distribution of SIF values
#SIF.fulls <- matrix(NA, nrow(index),120000)

sink("model.txt")
cat("model {
#Specify uninformative priors for detection and occupancy parameters
pA ~ dunif(0, 1)
  	pB ~ dunif(0, 1)
	psiA ~ dunif(0, 1) 
  	psiBA ~ dunif(0, 1)
  	psiBa ~ dunif(0, 1)
#Derive the unconditional occupancy probability for species B 
psiB <- psiA * psiBA + (1 - psiA) * psiBa
     for (i in 1:n){	# n is the number of sites within a study area
zA[i] ~ dbern(psiA)	    # Model the true occurrence state for species A
		PSIB[i] <- zA[i] * psiBA + (1 - zA[i]) * psiBa
  	zB[i] ~ dbern(PSIB[i])    # Model the true occurrence state for species B
# Probability of observing a species given its true occurrence state and estimated 
# detection probability
    	  yA[i] ~ dbin(zA[i]*pA, N[i])  	# N is the total number of trap days
        yB[i] ~ dbin(zB[i]*pB, N[i])	
}
	# Derive the species interaction factor
SIF <-  psiA * psiBA / (psiA * (psiA * psiBA + (1 - psiA) * psiBa))
}", fill =TRUE)
sink()

#Loop through the matrix ‘index’ to run the independent analyses for all 1216 species pairs and #save the results to ‘result’

for(j in 1:nrow(index)){
  data1 <- list(yA = data[data$Spp1 == index$Spp1[j] & data$Spp2 == index$Spp2[j] & data$Site == index$Site[j], 'Spp1_Det'], 
                yB = data[data$Spp1 == index$Spp1[j] & data$Spp2 == index$Spp2[j]& data$Site == index$Site[j], 'Spp2_Det'], 
                N = data[data$Spp1 == index$Spp1[j] & data$Spp2 == index$Spp2[j]& data$Site == index$Site[j], 'Samp.pers'], 
                n = length(unique(data[data$Spp1 == index$Spp1[j] & data$Spp2 == index$Spp2[j] & data$Site == index$Site[j],'CT_ID'])))
  
  inits <-function() {list(zA = as.numeric(as.numeric(data1$yA>0)), zB = as.numeric(as.numeric(data1$yB>0)))}
  
  out <- jags(parameters.to.save = parameters, "model.txt", inits = inits, data = data1, n.chains = nchain, n.burnin = nburn, n.iter = niter, n.thin = nthin)
  
  ##Extract SIF values of interest- mean, standard deviation, lower credible interval, upper credible interval, RHat
  mean <- out$BUGSoutput$summary[1,1]
  std.dev <- out$BUGSoutput$summary[1,2]
  lci <- out$BUGSoutput$summary[1,3]
  hci <- out$BUGSoutput$summary[1,7]
  Rhat <- out$BUGSoutput$summary[1,8]
  outsav <- cbind(mean,std.dev,lci,hci,Rhat)
  result[j, ] <- outsav
  
  #SIF.fulls[j,] <- out$BUGSoutput$sims.list$SIF[,1]
}

#Bind the ‘result’ with the respective row in ‘index’
result  <-  cbind(index, result)

#Rename the columns of ‘result’
colnames(result) <-  c('Site.Code', 'Spp1', 'Spp2', 'Mean.SIF', 'sd.SIF', 'lci.SIF', 'hci.SIF', 'RHat')

max(result$RHat)

write.csv(result, file = "SIF.models.2015a.csv")




###############################
##Linear regression
##############################

##Import SIF values for 2015

SIF.models.2015 <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/SIF.models.2015a.csv")

SIF <- SIF.models.2015



####Get covariates####

##Functional dissimilarity##

library(dplyr)
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


#Effect of poaching?
#Poaching <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/Poaching.csv")
#Poaching <- Poaching[,c(4,10,12)]
#colnames(Poaching) <- c("Spp1","Hunted1","Site.Code")
#SIF <- left_join(SIF, Poaching , by =c("Spp1","Site.Code"))
#colnames(Poaching) <- c("Spp2","Hunted2","Site.Code")
#SIF <- left_join(SIF, Poaching , by =c("Spp2","Site.Code"))
#SIF$Hunted.pr <- paste(SIF$Hunted1,SIF$Hunted2, sep = "_")
#SIF$Hunted.pr[SIF$Hunted.pr == "N_U"|SIF$Hunted.pr == "U_N"|
#              SIF$Hunted.pr == "U_U"|SIF$Hunted.pr == "U_Y"|
#              SIF$Hunted.pr == "Y_U"] <- "Unknown"
#SIF$Hunted.pr[SIF$Hunted.pr == "N_Y"|SIF$Hunted.pr == "Y_N"] <- "N_Y"

#library(ggplot2)
#ggplot() + 
#  geom_boxplot(data=SIF, aes(x=Hunted.pr, y=Mean.SIF))



##Quick calculations
median(SIF$Mean.SIF)
quantile(SIF$Mean.SIF,probs=c(.025,.975))
median(SIF$sd.SIF)
quantile(SIF$sd.SIF,probs=c(.025,.975))
aggregate <- filter(SIF,SIF$lci.SIF > 1)
segregate<- filter(SIF,SIF$hci.SIF < 1)

library(dplyr)
site.med.SIF <- vector("numeric")
for (i in unique(SIF$Site.Code)) {
  SIF.site <- filter(SIF, SIF$Site.Code == paste(i))
  SIF.med <- median(SIF.site$Mean.SIF)
  site.med.SIF[i] <- SIF.med
  
}

##Histograms of SIF values (Figure S1)
library(ggplot2)
SIF.hist <- SIF
SIF.hist$sig <- ifelse(SIF.hist$lci.SIF > 1,"agg",ifelse(SIF.hist$hci.SIF < 1, "aseg","none"))
ggplot()+
  geom_histogram(data=SIF.hist, aes(x=Mean.SIF, fill = sig))+
  scale_fill_manual(values = c("red","blue","grey"), label = c("95% aggregation","95% segregation","none"))+
  theme_classic()+
  xlab("Mean species interaction factor")+
  ylab("Number of species pairs")+
  labs(fill = "Spatial association")



##Site level covariates##
Site.covs.2015 <- read.csv("~/Documents/TEAM_co-occurence/TEAM_co-occurence/Site.covs.2015.csv")
SIF1 <- left_join(SIF,Site.covs.2015, by = "Site.Code")


##Medians
Site.covs.med <- Site.covs.2015[-c(1,8),]
median(Site.covs.med$Hab.div.2015)
median(Site.covs.med$NDVI.2015)
median(Site.covs.med$Hum.Dens)
median(Site.FRich$Func.Rich)

###Format model summaries for regression###

SIF1$Spp.Site <- paste(SIF1$Spp1,SIF1$Spp2,SIF1$Site.Code, sep = "_")
SIF2 <- SIF1[,c(1,4,5,9:14)]
SIF2$logSIF <- log(SIF2$Mean.SIF)

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


weighted.out$BUGSoutput$summary
head(weighted.out$BUGSoutput$summary, 10)



##FIG 3- heat map
int.mcmc <- as.mcmc(weighted.out)
int.mcmc.mat <- as.matrix(int.mcmc)
int.mcmc.dat <- as.data.frame(int.mcmc.mat)

x2.sim <- seq(min(data$Hum.Dens), max(data$Hum.Dens), by = 0.05)
x3.sim <- seq(min(data$Func.Diss), max(data$Func.Diss), by = 0.05)

int.sim <- matrix(NA, nrow = length(x3.sim),ncol= length(x2.sim)) 
int.sim <- as.data.frame(int.sim)                  

for(i in 1:length(x2.sim)){
  for (j in 1:length(x3.sim)){
    int.sim[j, i] <- weighted.out$BUGSoutput$mean$alpha + weighted.out$BUGSoutput$mean$bFuncDiss*x3.sim[j] + weighted.out$BUGSoutput$mean$bFDHum*x3.sim[j]*x2.sim[i] + weighted.out$BUGSoutput$mean$bHumDens*x2.sim[i]
  }}

rownames(int.sim) <- seq(min(data$Func.Diss), max(data$Func.Diss), by = 0.05)
colnames(int.sim) <- seq(min(data$Hum.Dens), max(data$Hum.Dens), by = 0.05)

library(viridis)
int.sim.long <- gather(int.sim)
int.sim.long$Func.Diss <- rep.int(x3.sim,71)
colnames(int.sim.long) <- c("Hum.Dens","SIF","Func.Diss")
int.sim.long$Hum.Dens <- as.numeric(int.sim.long$Hum.Dens)
ggplot(int.sim.long)+
  geom_tile(aes(x=Hum.Dens,y=Func.Diss,fill=SIF))+
  scale_fill_viridis_b()+
  theme_classic()+
  xlab("Scaled human density")+
  ylab("Scaled ecological dissimilarity")+
  labs(fill = 
         "Species 
interaction 
factor")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))







##Posterior predictive check

mu.runs <- as.data.frame(weighted.out$BUGSoutput$sims.list$mu)
mu.set <- sample(seq(from = 1, to = 45000), size = 8, replace = FALSE)
mu.8 <- mu.runs[mu.set,]
mu.8 <- as.data.frame(t(mu.8))
colnames(mu.8) <- c('run1','run2','run3','run4','run5','run6','run7','run8')
data.fit <- cbind(data,mu.8)

library(ggplot2)
library(cowplot)
data.hist<-ggplot() +
  geom_histogram(data.fit, mapping = aes(x = logSIF, weight = 1/(sd.SIF)^2), bins = 100)+
  xlim(-0.1,0.1)+
  ylab(NULL)+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
mod1<-ggplot() +
  geom_histogram(data.fit, mapping = aes(x = run1), bins = 100)+
  xlim(-0.1,0.1)+
  ylab(NULL)+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
mod2 <- ggplot() +
  geom_histogram(data.fit, mapping = aes(x = run2), bins = 100)+
  xlim(-0.1,0.1)+
  ylab(NULL)+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
mod3 <- ggplot() +
  geom_histogram(data.fit, mapping = aes(x = run3), bins = 100)+
  xlim(-0.1,0.1)+
  ylab(NULL)+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
mod4 <- ggplot() +
  geom_histogram(data.fit, mapping = aes(x = run4), bins = 100)+
  xlim(-0.1,0.1)+
  ylab(NULL)+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
mod5 <- ggplot() +
  geom_histogram(data.fit, mapping = aes(x = run5), bins = 100)+
  xlim(-0.1,0.1)+
  ylab(NULL)+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
mod6 <- ggplot() +
  geom_histogram(data.fit, mapping = aes(x = run6), bins = 100)+
  xlim(-0.1,0.1)+
  ylab(NULL)+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
mod7 <- ggplot() +
  geom_histogram(data.fit, mapping = aes(x = run7), bins = 100)+
  xlim(-0.1,0.1)+
  ylab(NULL)+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
mod8 <- ggplot() +
  geom_histogram(data.fit, mapping = aes(x = run8), bins = 100)+
  xlim(-0.1,0.1)+
  ylab(NULL)+
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
plot_grid(data.hist,mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)




##Plot relationship (unweighted) between functional dissimilarity and log SIF at different sites

mu <- as.data.frame(weighted.out$BUGSoutput$summary[11:1226,c(1,3,7)])
colnames(mu) <- c("mod.mean","mod.2.5","mod.97.5")
data.fit <- cbind(data,mu)

library(patchwork)
library(gridExtra)
library(grid)
library(ggplot2)
library(gtable)
library(png)
library(cowplot)



p <- ggplot(data.fit, aes(x = Func.Diss, group = Site.Code)) +
  geom_line(aes(y = 0), color = "black") +
  geom_line(aes(y = mod.mean), color = "red") +
  geom_ribbon(aes(y = mod.mean, ymin = mod.2.5, ymax = mod.97.5),fill = "red", alpha = 0.3)+
  geom_point(aes(y = logSIF),alpha = 0.3, size = 1) +
  facet_wrap(~Site.Code, scales = "fixed") +
  theme_classic()+
  ylim(-0.175,0.175)+
  xlab("Functional dissimilarity")+
  ylab("Species interaction factor (SIF)")


#img.BBS <- readPNG("Human density maps/BBS.png")
img.BCI <- readPNG("Human density maps/BCI.png")
img.BIF <- readPNG("Human density maps/BIF.png")
img.CAX <- readPNG("Human density maps/CAX.png")
img.COU <- readPNG("Human density maps/COU.png")
img.CSN <- readPNG("Human density maps/CSN.png")
img.KRP <- readPNG("Human density maps/KRP.png")
img.NNN <- readPNG("Human density maps/NNN.png")
img.PSH <- readPNG("Human density maps/PSH.png")
img.RNF <- readPNG("Human density maps/RNF.png")
img.UDZ <- readPNG("Human density maps/UDZ.png")
img.VB <- readPNG("Human density maps/VB.png")
img.YAN <- readPNG("Human density maps/YAN.png")
img.YAS <- readPNG("Human density maps/YAS.png")

#Convert images to Grob (graphical objects)
#grobBBS <- rasterGrob(img.BBS)
grobBCI <- rasterGrob(img.BCI)
grobBIF <- rasterGrob(img.BIF)
grobCAX <- rasterGrob(img.CAX)
grobCOU <- rasterGrob(img.COU)
grobCSN <- rasterGrob(img.CSN)
grobKRP <- rasterGrob(img.KRP)
grobNNN <- rasterGrob(img.NNN)
grobPSH <- rasterGrob(img.PSH)
grobRNF <- rasterGrob(img.RNF)
grobUDZ <- rasterGrob(img.UDZ)
grobVB <- rasterGrob(img.VB)
grobYAN <- rasterGrob(img.YAN)
grobYAS <- rasterGrob(img.YAS)

#Convert park labels to Grob (text objects)
#grobtBBS <- textGrob("Bukit Barisan, Indonesia",
#                     gp = gpar(col = "black", fontsize = 7))
grobtBCI <- textGrob("Barro Colorado, Panama",
                     gp = gpar(col = "black", fontsize = 7))
grobtBIF <- textGrob("Bwindi, Uganda",
                     gp = gpar(col = "black", fontsize = 7))
grobtCAX <- textGrob("Caxiunã, Brazil",
                     gp = gpar(col = "black", fontsize = 7))
grobtCOU <- textGrob("Cocha Cashu, Peru",
                     gp = gpar(col = "black", fontsize = 7))
grobtCSN <- textGrob("Central Suriname, Suriname",
                     gp = gpar(col = "black", fontsize = 7))
grobtKRP <- textGrob("Korup, Cameroon",
                     gp = gpar(col = "black", fontsize = 7))
grobtNNN <- textGrob("Noubalé-Ndoki, Republic of Congo",
                     gp = gpar(col = "black", fontsize = 7))
grobtPSH <- textGrob("Pasoh, Malaysia",
                     gp = gpar(col = "black", fontsize = 7))
grobtRNF <- textGrob("Ranomafana, Madagscar",
                     gp = gpar(col = "black", fontsize = 7))
grobtUDZ <- textGrob("Udzungwa, Tanzania",
                     gp = gpar(col = "black", fontsize = 7))
grobtVB <- textGrob("Volcán Barva, Costa Rica",
                    gp = gpar(col = "black", fontsize = 7))
grobtYAN <- textGrob("Yanachanga, Peru",
                     gp = gpar(col = "black", fontsize = 7))
grobtYAS <- textGrob("Yasuní, Ecuador",
                     gp = gpar(col = "black", fontsize = 7))



# convert the plot to gtable
mytable <- ggplot_gtable(ggplot_build(p))

#Add a line ssame height as line 4 after line 3
# use gtable_show_layout(mytable) to see where you want to put your line
mytable <- gtable_add_cols(mytable,unit(1,"null"), 6)
mytable <- gtable_add_cols(mytable,unit(1,"null"), 12)
mytable <- gtable_add_cols(mytable,unit(1,"null"), 18)
mytable <- gtable_add_cols(mytable,unit(1,"null"), 22)

mytable <- gtable_remove_grobs(mytable, c("strip-t-1-4","strip-t-2-4",
                                          "strip-t-3-4","strip-t-4-4",
                                          "strip-t-1-3","strip-t-2-3",
                                          "strip-t-3-3","strip-t-4-3",
                                          "strip-t-1-2","strip-t-2-2",
                                          "strip-t-3-2","strip-t-4-2",
                                          "strip-t-1-1","strip-t-2-1",
                                          "strip-t-3-1","strip-t-4-1"))
# if needed mor row can be added 
# (for example to keep consistent spacing with other graphs)

#Insert the grob in the cells of the new line
#mytable <- gtable_add_grob(mytable,grobBBS,8,7)
mytable <- gtable_add_grob(mytable,grobBCI,8,7)
mytable <- gtable_add_grob(mytable,grobBIF,8,13)
mytable <- gtable_add_grob(mytable,grobCAX,8,19)
mytable <- gtable_add_grob(mytable,grobCOU,8,23)
mytable <- gtable_add_grob(mytable,grobCSN,13,7)
mytable <- gtable_add_grob(mytable,grobKRP,13,13)
mytable <- gtable_add_grob(mytable,grobNNN,13,19)
mytable <- gtable_add_grob(mytable,grobPSH,13,23)
mytable <- gtable_add_grob(mytable,grobRNF,18,7)
mytable <- gtable_add_grob(mytable,grobUDZ,18,13)
mytable <- gtable_add_grob(mytable,grobVB,18,19)
mytable <- gtable_add_grob(mytable,grobYAN,18,23)
mytable <- gtable_add_grob(mytable,grobYAS,23,7)


#Insert the text grob in the cells over each site's plot and map
mytable <- gtable_add_grob(mytable,grobtBCI,7,5,7,7)
mytable <- gtable_add_grob(mytable,grobtBIF,7,10,7,13)
mytable <- gtable_add_grob(mytable,grobtCAX,7,15,7,19)
mytable <- gtable_add_grob(mytable,grobtCOU,7,20,7,23)
mytable <- gtable_add_grob(mytable,grobtCSN,12,5,12,7)
mytable <- gtable_add_grob(mytable,grobtKRP,12,10,12,13)
mytable <- gtable_add_grob(mytable,grobtNNN,12,15,12,19)
mytable <- gtable_add_grob(mytable,grobtPSH,12,20,12,23)
mytable <- gtable_add_grob(mytable,grobtRNF,17,5,17,7)
mytable <- gtable_add_grob(mytable,grobtUDZ,17,10,17,13)
mytable <- gtable_add_grob(mytable,grobtVB,17,15,17,19)
mytable <- gtable_add_grob(mytable,grobtYAN,17,20,17,23)
mytable <- gtable_add_grob(mytable,grobtYAS,22,5,22,7)


#rendering
grid.draw(mytable)




###Plot regression results [Fig. 2]###
library(bayesplot)
library(rstan)
library(rstantools)
library(brms)
library(ggplot2)
library(ggdist)
library(tidybayes)
library(tidyverse)


x <- weighted.out$BUGSoutput %>%
  spread_draws(alpha,bFuncDiss,bNDVI,bHumDens,
               bHabHet,bFRich,bFDNDVI,bFDHum,bFDHab)%>%
  gather(key = "Effect", value = "Est" ) #%>%
#x <- x[-c(1:135000),]
x <- x[-c(1:180000),]
x$Est <- as.numeric(x$Est)
#x$Effect <- factor(x$Effect, levels = c('alpha', 'bFuncDiss','bNDVI','bHumDens',
#               'bHabHet','bFRich','bFDNDVI','bFDHum','bFDHab'))
x$Effect <- factor(x$Effect, levels = c('bFuncDiss','bNDVI','bHumDens',
                                        'bHabHet','bFRich','bFDNDVI','bFDHum','bFDHab'))
library(ggplot2)
ggplot(data = x, aes(y = Effect, x = Est, fill = Effect)) +
  stat_halfeye(.width = c(.95, .5))+
  geom_vline(xintercept = 0)+
  labs(title="Species Interaction Factor 2015",
       x ="Beta effect size", y = "Predictor variables")+
  # scale_y_discrete(labels = c("alpha" = "alpha","Func.Diss" = "Func.Diss",
  #     "NDVI"="NDVI", "Hum.dens"="Hum.dens",
  #      "Hab.div"="Hab.div","Hab.het"="Hab.het", "FD.NDVI"="FD.NDVI",
  #      "FD.Hum"="FD.Hum","FD.Hab"="FD.Hab"))+
  scale_fill_manual(values = c("grey85","grey85","grey85","grey85","grey85","grey85","red","grey85"))+
  theme_classic()+
  theme(legend.position = "none")






#########
###END###
#########