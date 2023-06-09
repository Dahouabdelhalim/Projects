#Script for analysing the effect of inbreeding on morphological traits in house sparrow populations using multivariate models in MCMCglmm
#Stefanie Muff & Alina Niskanen, alina.niskanen@gmail.com
#April 2020

#Packages
library(nadiv)
library(INLA)
library(lme4)
library(MASS)
library(MasterBayes)
library(MCMCglmm)
library(parallel)
library(coda)


#Import data
#Import the phenotypic data, LRS and inbreeding estimates
data_temp <- read.table("Data.txt", header = T, stringsAsFactors = F)

str(data_temp)


###
###Subsetting the data
###

#Include only the individuals that have been observed as adults on one of the 8 study islands
data_adult_temp <- data_temp[which(data_temp$laflok %in% c("20","22","23","24","26","27","28","38")),]
table(data_adult_temp$laflok)
#Only birds that hatched on one of the 8 study islands
data_hatch_temp <- data_temp[which(data_temp$fiflok %in% c("20","22","23","24","26","27","28","38")),]
table(data_hatch_temp$fiflok)

###
###Include individuals with morphology information
###

#Only those individuals that have LRS_OK estimate and have been adults on one of the study islands
data_morph <- data_adult_temp[!is.na(data_adult_temp$age1mass),] 
table(is.na(data_temp$age1totba))
table(is.na(data_morph$age1totba))
table(data_morph$laflok) 
table(data_morph$fiflok)

#Remove temp files
rm(list=ls(pattern="temp"))

###
###Prepare the variables for the dataset that has correct morpohology data for adult islands
###

#Mean-center FGRM
data_morph$c_FGRM <- scale(data_morph$FGRM, center = T, scale = F)

# make hatchyear and island factor covariates
data_morph$f_all_hatchyears <- as.factor(data_morph$all_hatchyears)
data_morph$f_laflok <- as.factor(data_morph$laflok)
data_morph$f_fiflok <- as.factor(data_morph$fiflok)
data_morph$f_laflok2 <- data_morph$f_laflok
data_morph$f_all_hatchyears2 <- data_morph$f_all_hatchyears


#Add combined variable for hatch year and island
data_morph$IY <- as.factor(paste(data_morph$all_hatchyears,data_morph$fiflok,sep="_"))
data_morph$IAY <- as.factor(paste(data_morph$all_hatchyears,data_morph$laflok,sep="_"))

###
###Import the SNP-pedigree
###
d.ped <- read.table("pedigree.txt", header = T, sep=" ")
names(d.ped) <- c("ringnr","mother","father")

# matrix inversion from the MCMCglmm package; but first need to order the pedigree using MasterBayes package
d.ped <- orderPed(d.ped)
invA <- inverseA(d.ped[,c("ringnr","mother","father")])
Cmatrix <- invA$Ainv

# use dimnames to sort and cut the matrix so that it fits with the data_morph file:
Cmatrix <- Cmatrix[data_morph$id,data_morph$id]


#############
#morphology multivariate runs for different trait combinations
#############

#Prior specification for two-variate case
priorTwo<-list(G=list(G1=list(V=diag(2), nu=0.001, alpha.mu=rep(0,2), alpha.V=diag(2)*1000),
                      G2=list(V=diag(2), nu=0.001, alpha.mu=rep(0,2), alpha.V=diag(2)*1000),
                      G3=list(V=diag(2), nu=0.001, alpha.mu=rep(0,2), alpha.V=diag(2)*1000),
                      G4=list(V=diag(2), nu=0.001, alpha.mu=rep(0,2), alpha.V=diag(2)*1000)),
               R=list(V=diag(2), nu=0.001))

#Prior specification for three-variate case
priorTri<-list(G=list(G1=list(V=diag(3), nu=0.001, alpha.mu=rep(0,3), alpha.V=diag(3)*1000),
                      G2=list(V=diag(3), nu=0.001, alpha.mu=rep(0,3), alpha.V=diag(3)*1000),
                      G3=list(V=diag(3), nu=0.001, alpha.mu=rep(0,3), alpha.V=diag(3)*1000),
                      G4=list(V=diag(3), nu=0.001, alpha.mu=rep(0,3), alpha.V=diag(3)*1000)),
               R=list(V=diag(3), nu=0.001))



###########################
##Run first model for Tarsus, wing and billD
###########################


##First run for first model for Tarsus, wing and billD


r.trivariate_1.1 <- MCMCglmm(cbind(age1tarsus,age1wing,age1billD) ~ trait -1 + trait:c_FGRM  + trait:gen_sex,
                         random =~ us(trait):f_all_hatchyears + 
                           us(trait):f_laflok + 
                           us(trait):id + 
                           us(trait):IAY, 
                         rcov=~us(trait):units,
                         prior = priorTri, 
                         data=data_morph,
                         ginverse = list(id = Cmatrix),
                         family=c("gaussian","gaussian","gaussian"),
                         verbose=TRUE,nitt=100000,burnin=20000,thin=20)


summary(r.trivariate_1.1)


#plotting the trace
par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.trivariate_1.1$Sol, auto.layout=F)

##Second run for first model for Tarsus, wing and billD

system.time(
r.trivariate_1.2 <- MCMCglmm(cbind(age1tarsus,age1wing,age1billD) ~ trait -1 + trait:c_FGRM  + trait:gen_sex,
                             random =~ us(trait):f_all_hatchyears + 
                               us(trait):f_laflok + 
                               us(trait):id + 
                               us(trait):IAY, 
                             rcov=~us(trait):units,
                             prior = priorTri, 
                             data=data_morph,
                             ginverse = list(id = Cmatrix),
                             family=c("gaussian","gaussian","gaussian"),
                             verbose=TRUE,nitt=100000,burnin=20000,thin=20)
)


summary(r.trivariate_1.2)


#plotting the trace
par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.trivariate_1.1$Sol, auto.layout=F)


system.time(
  r.trivariate_1.3 <- MCMCglmm(cbind(age1tarsus,age1wing,age1billD) ~ trait -1 + trait:c_FGRM  + trait:gen_sex,
                               random =~ us(trait):f_all_hatchyears + 
                                 us(trait):f_laflok + 
                                 us(trait):id + 
                                 us(trait):IAY, 
                               rcov=~us(trait):units,
                               prior = priorTri, 
                               data=data_morph,
                               ginverse = list(id = Cmatrix),
                               family=c("gaussian","gaussian","gaussian"),
                               verbose=TRUE,nitt=100000,burnin=20000,thin=20)
)

summary(r.trivariate_1.3)


#Combine the different runs
#First combine the different runs into a list, then call the traces of all runs
m3 <- list(r.trivariate_1.1,r.trivariate_1.2,r.trivariate_1.3)
m3 <- lapply(m3, function(m) m$Sol)
m3 <- do.call(mcmc.list, m3)


#plotting the traces and convergence
gelman.plot(m3, auto.layout=F)

plot(m3, ask=F, auto.layout=F)

summary(m3)
gelman.diag(m3)


###########################
#Mass, tarsus & bill length
###########################

##First run 2nd model: mass, tarsus and billL

r.trivariate_2.1 <- MCMCglmm(cbind(age1mass,age1tarsus,age1billL) ~ trait -1 + trait:c_FGRM  + trait:gen_sex,
                             random =~ us(trait):f_all_hatchyears + 
                               us(trait):f_laflok + 
                               us(trait):id + 
                               us(trait):IAY, 
                             rcov=~us(trait):units,
                             prior = priorTri, 
                             data=data_morph,
                             ginverse = list(id = Cmatrix),
                             family=c("gaussian","gaussian","gaussian"),
                             verbose=TRUE,nitt=100000,burnin=20000,thin=20)

summary(r.trivariate_2.1)


#plotting the trace
par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.trivariate_2.1$Sol, auto.layout=F)


##Second run 2nd model: mass, tarsus and billL

system.time(
r.trivariate_2.2 <- MCMCglmm(cbind(age1mass,age1tarsus,age1billL) ~ trait -1 + trait:c_FGRM  + trait:gen_sex,
                             random =~ us(trait):f_all_hatchyears + 
                               us(trait):f_laflok + 
                               us(trait):id + 
                               us(trait):IAY, 
                             rcov=~us(trait):units,
                             prior = priorTri, 
                             data=data_morph,
                             ginverse = list(id = Cmatrix),
                             family=c("gaussian","gaussian","gaussian"),
                             verbose=TRUE,nitt=100000,burnin=20000,thin=20)
)

summary(r.trivariate_2.2)


#plotting the trace
par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.trivariate_2.2$Sol, auto.layout=F)


##Third run 2nd model: mass, tarsus and billL

system.time(
  r.trivariate_2.3 <- MCMCglmm(cbind(age1mass,age1tarsus,age1billL) ~ trait -1 + trait:c_FGRM  + trait:gen_sex,
                               random =~ us(trait):f_all_hatchyears + 
                                 us(trait):f_laflok + 
                                 us(trait):id + 
                                 us(trait):IAY, 
                               rcov=~us(trait):units,
                               prior = priorTri, 
                               data=data_morph,
                               ginverse = list(id = Cmatrix),
                               family=c("gaussian","gaussian","gaussian"),
                               verbose=TRUE,nitt=100000,burnin=20000,thin=20)
)

summary(r.trivariate_2.3)


#plotting the trace
par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.trivariate_2.3$Sol, auto.layout=F)


#Combine the different runs
#First combine the different runs into a list, then call the traces of all runs
m3.2 <- list(r.trivariate_2.1,r.trivariate_2.2,r.trivariate_2.3)
m3.2 <- lapply(m3.2, function(m) m$Sol)
m3.2 <- do.call(mcmc.list, m3.2)


#plotting the traces and convergence
gelman.plot(m3.2, auto.layout=F)

plot(m3.2, ask=F, auto.layout=F)

summary(m3.2)

gelman.diag(m3.2)


###########################
#Tarsus, bill length & bill depth
###########################

##First run 3rd model: Tarsus, bill length & bill depth

r.trivariate_3.1 <- MCMCglmm(cbind(age1tarsus,age1billL,age1billD) ~ trait -1 + trait:c_FGRM  + trait:gen_sex,
                             random =~ us(trait):f_all_hatchyears + 
                               us(trait):f_laflok + 
                               us(trait):id + 
                               us(trait):IAY, 
                             rcov=~us(trait):units,
                             prior = priorTri, 
                             data=data_morph,
                             ginverse = list(id = Cmatrix),
                             family=c("gaussian","gaussian","gaussian"),
                             verbose=TRUE,nitt=100000,burnin=20000,thin=20)

summary(r.trivariate_3.1)


#plotting the trace
par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.trivariate_3.1$Sol, auto.layout=F)


##Second run 3rd model: Tarsus, bill length & bill depth

r.trivariate_3.2 <- MCMCglmm(cbind(age1tarsus,age1billL,age1billD) ~ trait -1 + trait:c_FGRM  + trait:gen_sex,
                             random =~ us(trait):f_all_hatchyears + 
                               us(trait):f_laflok + 
                               us(trait):id + 
                               us(trait):IAY, 
                             rcov=~us(trait):units,
                             prior = priorTri, 
                             data=data_morph,
                             ginverse = list(id = Cmatrix),
                             family=c("gaussian","gaussian","gaussian"),
                             verbose=TRUE,nitt=100000,burnin=20000,thin=20)

summary(r.trivariate_3.2)


#plotting the trace
par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.trivariate_3.2$Sol, auto.layout=F)



##Third run 3rd model: Tarsus, bill length & bill depth

system.time(
r.trivariate_3.3 <- MCMCglmm(cbind(age1tarsus,age1billL,age1billD) ~ trait -1 + trait:c_FGRM  + trait:gen_sex,
                             random =~ us(trait):f_all_hatchyears + 
                               us(trait):f_laflok + 
                               us(trait):id + 
                               us(trait):IAY, 
                             rcov=~us(trait):units,
                             prior = priorTri, 
                             data=data_morph,
                             ginverse = list(id = Cmatrix),
                             family=c("gaussian","gaussian","gaussian"),
                             verbose=TRUE,nitt=100000,burnin=20000,thin=20)
)

summary(r.trivariate_3.3)


#plotting the trace
par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.trivariate_3.3$Sol, auto.layout=F)


#Combine the different runs
#First combine the different runs into a list, then call the traces of all runs
m3.3 <- list(r.trivariate_3.1,r.trivariate_3.2,r.trivariate_3.3)
m3.3 <- lapply(m3.3, function(m) m$Sol)
m3.3 <- do.call(mcmc.list, m3.3)

#plotting the traces and convergence

gelman.plot(m3.3, auto.layout=F)

plot(m3.3, ask=F, auto.layout=F)

summary(m3.3)

gelman.diag(m3.3)



#################
#Visible badge and total badge size, only males
################

#First run of model for two variates: visible badge and total badge size 
system.time(
  r.twoMvariate_2.1 <- MCMCglmm(cbind(age1totba,age1visba) ~ trait -1 + trait:c_FGRM,
                                random =~ us(trait):f_all_hatchyears + 
                                  us(trait):f_laflok + 
                                  us(trait):id + 
                                  us(trait):IAY, 
                                rcov=~us(trait):units,
                                prior = priorTwo, 
                                data=data_males,
                                ginverse = list(id = Cmatrix),
                                family=c("gaussian","gaussian"),
                                verbose=TRUE,nitt=100000,burnin=20000,thin=20)
)


summary(r.twoMvariate_2.1)


#Plotting the convergence
summary(r.twoMvariate_2.1$Sol)
par("mar")

par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.twoMvariate_2.1$Sol, auto.layout=F)


#Second run of model for two variates: visible badge and total badge size 
system.time(
  r.twoMvariate_2.2 <- MCMCglmm(cbind(age1totba,age1visba) ~ trait -1 + trait:c_FGRM,
                                random =~ us(trait):f_all_hatchyears + 
                                  us(trait):f_laflok + 
                                  us(trait):id + 
                                  us(trait):IAY, 
                                rcov=~us(trait):units,
                                prior = priorTwo, 
                                data=data_males,
                                ginverse = list(id = Cmatrix),
                                family=c("gaussian","gaussian"),
                                verbose=TRUE,nitt=100000,burnin=20000,thin=20)
)

summary(r.twoMvariate_2.2)


#Plotting the convergence
summary(r.twoMvariate_2.2$Sol)
par("mar")

par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.twoMvariate_2.2$Sol, auto.layout=F)


#Third run of model for two variates: visible badge and total badge size 
system.time(
  r.twoMvariate_2.3 <- MCMCglmm(cbind(age1totba,age1visba) ~ trait -1 + trait:c_FGRM,
                                random =~ us(trait):f_all_hatchyears + 
                                  us(trait):f_laflok + 
                                  us(trait):id + 
                                  us(trait):IAY, 
                                rcov=~us(trait):units,
                                prior = priorTwo, 
                                data=data_males,
                                ginverse = list(id = Cmatrix),
                                family=c("gaussian","gaussian"),
                                verbose=TRUE,nitt=100000,burnin=20000,thin=20)
)

summary(r.twoMvariate_2.3)


#Plotting the convergence
summary(r.twoMvariate_2.3$Sol)
par("mar")

par(mfrow=c(6,2), mar=c(2,2,1,0))
plot(r.twoMvariate_2.3$Sol, auto.layout=F)


#Combine the different runs for bivariate badge size measurements
#First combine the different runs into a list, then call the traces of all runs
m_badge <- list(r.twoMvariate_2.1,r.twoMvariate_2.2,r.twoMvariate_2.3)
m_badge <- lapply(m_badge, function(m) m$Sol)
m_badge <- do.call(mcmc.list, m_badge)


#plotting the traces and convergence
gelman.plot(m_badge, auto.layout=F)

plot(m_badge, ask=F, auto.layout=F)

summary(m_badge)

gelman.diag(m_badge)



