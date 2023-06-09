######------ PART A: BAYESIAN HIERARCHICAL OCCUPANCY MODEL FOR JAGUAR HABITAT USE: WITH FALSE POSITIVES ------######

library(jagsUI)
setwd("C:/Users/lpetracca/Desktop/git/Jaguar")

#Read in site covariates and standardize them
as.data.frame(site_cov<-read.csv("SiteCovariates.csv", row.names=1))
head(site_cov, n=5)

Elev<-as.numeric(scale(site_cov$Elev))
DistStrictPA<-as.numeric(scale(site_cov$Dist_StrictPA))
DistSettle<-as.numeric(scale(site_cov$Dist_MajorSettlement))
Prey<-as.numeric(scale(site_cov$Prey))
PercCover<-as.numeric(site_cov$PercCover/100)
Agropast<-as.numeric(site_cov$Agropast)
CorridorTest<-as.numeric(site_cov$CorridorTest)

#Import detection history
as.data.frame(detect<-read.csv("Jaguar_Max6_FP.csv"))
detect[detect=="-"] <-NA 
detect<-apply(detect, 2, function(x){as.numeric(x)})
as.data.frame(detect)
y <- detect
y

#Import corridor number for detection probability (six columns)
as.data.frame(corridor<-read.csv("Corridor_Max6_FP.csv", row.names=1))
corridor[corridor=="-"]<-NA
corridor[(is.na(corridor))] <- 6
corridor<-apply(corridor, 2, function(x){as.numeric(x)}) 
as.data.frame(corridor)
Corridor <- as.numeric(corridor+1)
dim(Corridor)
ncorridors <- max(as.numeric(Corridor))

#Import effort data
as.data.frame(effort<-read.csv("Effort_Max6_FP.csv", row.names=1))
effort[effort=="-"]<-NA
effort<-apply(effort, 2, function(x){as.numeric(x)}) 
#apply mean to missing values
effort[(is.na(effort))] <- 0.545
as.data.frame(effort)

#Import data on whether or not there were four than four interviews collected in the cell
as.data.frame((morethanfour<-read.csv("MoreThanFour_v3.csv", row.names=1, header=TRUE)))
morethanfour[morethanfour=="-"]<-NA
morethanfour<-apply(morethanfour, 2, function(x){as.numeric(x)}) 
#apply mean to missing values
morethanfour[(is.na(morethanfour))] <- 0.67755
as.data.frame(morethanfour)

#Import field season data for detection (1 column)
as.data.frame(FieldSeason<-read.csv("FieldSeason.csv", header=FALSE))
FieldSeason<-apply(FieldSeason, 2, function(x){as.numeric(x)})
as.data.frame(FieldSeason)
NFieldSeason <- max(as.numeric(FieldSeason))

# Define model
sink("with_FPs.txt")
cat("
    model {
    
    #Priors
    #these are the random intercepts on psi for corridor
    for (i in 1:ncorridors){
    alpha.occ[i] ~ dnorm(mu.int, tau.int)
    }
    
    #these are the random intercepts on p for field season
    for (i in 1:NFieldSeason){
    alpha.p[i] ~ dnorm(mu.p, tau.p)
    }
    
    #hyperparameter on psi
    mu.int ~ dnorm(0, 0.001)
    tau.int <- 1 / (sigma.int * sigma.int)
    sigma.int ~ dunif(0, 100)
    
    #hyperparameter on p
    mu.p ~ dunif(0,10)  
    tau.p <- 1 / (sigma.p * sigma.p)
    sigma.p ~ dunif(0, 1)
    
    #priors for beta parameters
    agropast.beta ~ dnorm(0, 0.001)
    distsettle.beta ~ dnorm(0, 0.001)
    perccover.beta ~ dnorm(0, 0.001)
    prey.beta ~ dnorm(0, 0.001)
    strictpa.beta ~ dnorm(0, 0.001)
    elev.beta ~ dnorm(0, 0.001)

    effort.beta ~ dnorm(0, 0.001)
    
    #prior for probability of false positives
    alpha.fp ~ dunif(-10, 0)
    #prior for beta on false positive probability
    morethanfour.beta ~ dunif(0, 10)
    
    # Likelihood
    for (i in 1:R) { #start initial loop over the R sites
    # True state model for the partially observed true state
    logitpsi[i]<- alpha.occ[corridor[i]] + agropast.beta * Agropast[i] + distsettle.beta * DistSettle[i] + perccover.beta * PercCover[i] + prey.beta * Prey[i] + strictpa.beta * DistStrictPA[i] + elev.beta * Elev[i] 
    # Code to avoid issues with logit
    logitpsitrun[i]<-min(999,max(-999,logitpsi[i]))
    psi[i]<-1/(1+exp(-logitpsitrun[i]))
    z[i] ~ dbern(psi[i])		# True occupancy z at site i
    
    # Model for false positive probability, based on whether there were >4 interviews conducted in the unit
    logitfp[i]<- alpha.fp + morethanfour.beta * morethanfour[i]
    logitfptrun[i]<-min(999,max(-999,logitfp[i]))
    fp[i]<-1/(1+exp(-logitfptrun[i]))
    
    # Model for detection probability, based on field season and effort
    for (j in 1:T) { # nested random effects
    y[i,j] ~ dbern(eff.p[i,j])
    eff.p[i,j] <- z[i] * p[i,j] + (1-z[i]) * fp[i]
    logitp[i,j]<-alpha.p[FieldSeason[i]] + effort.beta * effort[i,j]
    # Code to avoid issues with logit
    logitptrun[i,j]<-min(999,max(-999,logitp[i,j]))
    p[i,j]<-1/(1+exp(-logitptrun[i,j]))
    
    
    }}
    
    # Derived quantities
    occ.fs <- sum(z[])			# Number of occupied sites among 150
    }
    ",fill=TRUE)
sink()

# Bundle data
win.data <- list(y=y, effort=effort, morethanfour=morethanfour, Agropast=Agropast, DistSettle=DistSettle, PercCover=PercCover, Prey=Prey, DistStrictPA = DistStrictPA, Elev=Elev,  corridor=as.numeric(Corridor), FieldSeason=as.numeric(FieldSeason),
                 ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)),R = dim(y)[1], T = dim(y)[2])

# Inits function
zst <- apply(y, 1, max)	
zst[is.na(zst)] <- 1
inits <- function(){list(z = zst, mu.int = rnorm(1,0,2), sigma.int= rlnorm(1), alpha.occ = rnorm(ncorridors, 0, 2), agropast.beta = rnorm(1, 0, 2), distsettle.beta = rnorm(1, 0, 2), perccover.beta = rnorm(1, 0, 2), prey.beta = rnorm(1, 0, 2), 
                         strictpa.beta = rnorm(1, 0, 2), elev.beta = rnorm(1, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7),alpha.p = runif(NFieldSeason,0,10), effort.beta= rnorm(1, 0, 2), alpha.fp = runif(1, -4, -1), morethanfour.beta= runif(1, 0, 3.99))}

# Parameters to estimate
params <- c("mu.int", "alpha.occ", "agropast.beta", "distsettle.beta", "perccover.beta", "prey.beta", "strictpa.beta", "elev.beta","mu.p", "alpha.p", "effort.beta", "alpha.fp", "morethanfour.beta")

# MCMC settings
nc <- 3
nb <- 40000
ni <- 400000
nt <- 6

# Start Gibbs sampler
out <- jags(win.data, inits, params, "with_FPs.txt", n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt)

# Print output
print(out, dig = 3)



######------ PART B: BAYESIAN HIERARCHICAL OCCUPANCY MODEL FOR JAGUAR HABITAT USE: NO FALSE POSITIVES ------######

# Model is same as above but omits possibility of false positives

sink("without_FPs.txt")
cat("
    model {
    
    #Priors
    for (i in 1:ncorridors){
    alpha.occ[i] ~ dnorm(mu.int, tau.int)
    }
    
    for (i in 1:NFieldSeason){
    alpha.p[i] ~ dnorm(mu.p, tau.p)
    }
    
    mu.int ~ dnorm(0, 0.001)
    tau.int <- 1 / (sigma.int * sigma.int)
    sigma.int ~ dunif(0, 100)
    
    mu.p ~ dunif(0,10)  
    tau.p <- 1 / (sigma.p * sigma.p)
    sigma.p ~ dunif(0, 1)
    
    agropast.beta ~ dnorm(0, 0.001)
    distsettle.beta ~ dnorm(0, 0.001)
    perccover.beta ~ dnorm(0, 0.001)
    prey.beta ~ dnorm(0, 0.001)
    strictpa.beta ~ dnorm(0, 0.001)
    elev.beta ~ dnorm(0, 0.001)
    
    effort.beta ~ dnorm(0, 0.001)
    
    # Likelihood
    for (i in 1:R) { #start initial loop over the R sites
    # True state model for the partially observed true state
    logitpsi[i]<- alpha.occ[corridor[i]] + agropast.beta * Agropast[i] + distsettle.beta * DistSettle[i] + perccover.beta * PercCover[i] + prey.beta * Prey[i] + strictpa.beta * DistStrictPA[i] + elev.beta * Elev[i]
    logitpsitrun[i]<-min(999,max(-999,logitpsi[i]))
    psi[i]<-1/(1+exp(-logitpsitrun[i]))
    z[i] ~ dbern(psi[i])		# True occupancy z at site i
    
    for (j in 1:T) { # nested random effects
    
    y[i,j] ~ dbern(eff.p[i,j])
    eff.p[i,j] <- z[i] * p[i,j] 
    
    logitp[i,j]<-alpha.p[FieldSeason[i]] + effort.beta * effort[i,j]
    logitptrun[i,j]<-min(999,max(-999,logitp[i,j]))
    p[i,j]<-1/(1+exp(-logitptrun[i,j]))        
    
    }}
    
    # Derived quantities
    occ.fs <- sum(z[])			# Number of occupied sites among 150
    }
    ",fill=TRUE)
sink()

# Bundle data
win.data <- list(y=y, effort=effort, Agropast=Agropast, DistSettle=DistSettle, PercCover=PercCover, Prey=Prey, DistStrictPA = DistStrictPA, Elev=Elev,  corridor=as.numeric(Corridor), FieldSeason=as.numeric(FieldSeason),
                 ncorridors = max(as.numeric(Corridor)), NFieldSeason = max(as.numeric(FieldSeason)), R = dim(y)[1], T = dim(y)[2])

# Inits function
zst <- apply(y, 1, max)	
zst[is.na(zst)] <- 1
inits <- function(){list(z = zst, mu.int = rnorm(1,0,2), sigma.int= rlnorm(1), alpha.occ = rnorm(ncorridors, 0, 2), agropast.beta = rnorm(1, 0, 2), distsettle.beta = rnorm(1, 0, 2), perccover.beta = rnorm(1, 0, 2), prey.beta = rnorm(1, 0, 2), 
                         strictpa.beta = rnorm(1, 0, 2), elev.beta = rnorm(1, 0, 2), mu.p=runif(1,0,10), sigma.p=runif(1,0,0.7),alpha.p = runif(NFieldSeason,0,10), effort.beta= rnorm(1, 0, 2))}



# Parameters to estimate
params <- c("mu.int", "alpha.occ", "agropast.beta", "distsettle.beta", "perccover.beta", "prey.beta", "strictpa.beta", "elev.beta","mu.p", "alpha.p", "effort.beta")


# MCMC settings
nc <- 3
nb <- 10000
ni <- 50000
nt <- 3

# Start Gibbs sampler
out <- jags(win.data, inits, params, "without_FPs.txt", n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt)

print(out, dig = 3)
