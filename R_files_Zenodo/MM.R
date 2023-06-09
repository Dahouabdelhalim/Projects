# Fuentes-G JA, Polly PD, Martins EP. 2019. Evolution.
# A Bayesian extension of phylogenetic generalized least squares (PGLS): incorporating uncertainty in the comparative study of trait relationships and evolutionary rates

# Script modifies codes from:
# de Villemereuil P, Wells JA, Edwards RD, Blomberg SP. 2012. Bayesian models for comparative analysis integrating phylogenetic uncertainty. BMC Evolutionary Biology 12:102.
# Kruschke JK. 2011. Doing Bayesian Data Analysis: A Tutorial with R and BUGS. Academic Press / Elsevier.

require(rjags)
library(ape)
library(phytools)

# Generate mixture model with unequal rates file for JAGS
modelstring = "
model {
#Linear regression and multivariate normal likelihood
  for (i in 1:Ndata) {
    mu[i] <- beta0 + beta1 * x[i] + beta2 * z[i] + beta3 * xz[i]
  }
  y[1:Ndata] ~ dmnorm(mu[],TAU[,])
  
#Priors for regression coefficients (betas) and rates (gammas)
  beta0 ~ dnorm( 0 , 1.0E-06 )
  beta1 ~ dnorm( 0 , 1.0E-06 )
  beta2 ~ dnorm( 0 , 1.0E-06 )
  beta3 ~ dnorm( 0 , 1.0E-06 )
  sigma1 ~ dunif( 0 , 100 )
  sigma2 ~ dunif( 0 , 100 )
  gamma1 <- pow( sigma1 , 2 )
  gamma2 <- pow( sigma2 , 2 )

#Equal vector of probability for tree sampling
  for (k in 1:Ntree) {
    p[k] <- 1/Ntree
  }

#Tree sampling and variance-covariance matrix construction
  K ~ dcat(p[])
  PS <- (II*LAMBDA)+ID
  V1 <- gamma1*A1[,,K]
  V2 <- gamma2*A2[,,K]
  TAU <- inverse(PS*(V1 + V2))
  
# Mixture model:
  LAMBDA <- lambdaModel[ modelIndex ]
  lambdaModel[1] <- 0
  lambdaModel[2] <- 1
  lambdaModel[3] ~ dunif(0.01,0.99)
# Hyperprior:
  modelIndex ~ dcat( modelProb[] )
  modelProb[1] <- 1/3
  modelProb[2] <- 1/3
  modelProb[3] <- 1/3
}
"
writeLines(modelstring,con="multi.txt")


# Load data
Data = read.csv(file='Data.csv', header=TRUE)
Data<-Data[order(Data$Species),]
nSubj = NROW( Data )
x = Data[,"x"]
y = Data[,"y"]
z = Data[,"z"]

# Re-center and standardize data to optimize MCMC sampling
xM = mean( x ) ; xSD = sd( x )
yM = mean( y ) ; ySD = sd( y )
zx = ( x - xM ) / xSD
zy = ( y - yM ) / ySD

# Load stochastic character maps
mtrees<-read.simmap(file="Maps.tre",format="phylip")

# List of phylogenetic covariance matrices, one for each mapped state
Ntree<-NROW(mtrees)
A1=array(NA, dim=c(nSubj,nSubj,Ntree))
A2=array(NA, dim=c(nSubj,nSubj,Ntree))
for(i in 1:Ntree){
  G<-multiC(mtrees[[i]])
  G1<-G[[1]]
  G2<-G[[2]]
  G1 <- G1[order(dimnames(G1)[[1]], na.last=NA) , ]
  G2 <- G2[order(dimnames(G2)[[1]], na.last=NA) , ]
  G1 <- G1[,sort(colnames(G1))]
  G2 <- G2[,sort(colnames(G2))]
  A1[,,i]<-G1
  A2[,,i]<-G2
}


# Construct list of data
J<-matrix(1,nrow=nSubj,ncol=nSubj)
diag(J)<-0
dataList = list(
  x = zx ,
  y = zy ,
  z = z ,
  xz = zx*z ,
  Ndata = nSubj,
  Ntree = Ntree,
  A1 = A1,
  A2 = A2,
  II = J,
  ID = diag(nSubj)
)

# Initialize chains based on data standardization
r = cor(x,y)
initsList = list(
    beta0 = 0 ,
    beta1 = r ,
    beta2 = 0 ,
    beta3 = 0 ,
    sigma1 = 1-r^2,
    sigma2 = 1-r^2
)

# Set up chains
parameters = c("beta0" , "beta1" , "beta2"  , "beta3"  , "sigma1" , "sigma2", "LAMBDA", "modelIndex", "K")  # Saved parameters
adaptSteps = 1000             # Tunning steps
burnInSteps = 10000           # Burn-in steps
nChains = 3                   # Number of chains
numSavedSteps=100000          # Total steps
thinSteps=1                   # Thinning
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain
# Initialize model
jagsMulti = jags.model( "multi.txt" , data=dataList , inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in
update( jagsMulti , n.iter=burnInSteps )
# Sample final MCMC chain
SamplMulti = coda.samples( jagsMulti , variable.names=parameters , 
                            n.iter=nPerChain , thin=thinSteps )

# Examine parameters under standardized data
summary(SamplMulti)

# Convert coda-object to facilitate handling, and extract chain values
mcmcChain.m = as.matrix( SamplMulti )
z0.m = mcmcChain.m[, "beta0" ]
z1.m = mcmcChain.m[, "beta1" ]
z2.m = mcmcChain.m[, "beta2" ]
z3.m = mcmcChain.m[, "beta3" ]
zs1.m = mcmcChain.m[, "sigma1" ]
zs2.m = mcmcChain.m[, "sigma2" ]
lambda.m = mcmcChain.m[, "LAMBDA" ]

# Get the posterior sample of modelIndex
modelIdxSample.m = mcmcChain.m[,"modelIndex"]
# Compute probability of each model
pIID.m = sum( modelIdxSample.m == 1 ) / length( modelIdxSample.m )
pBM.m = sum( modelIdxSample.m == 2 ) / length( modelIdxSample.m )
pLM.m = sum( modelIdxSample.m == 3 ) / length( modelIdxSample.m )

# Convert to original scale
b1.m = z1.m * ySD / xSD
b0.m = ( z0.m * ySD + yM - z1.m * ySD * xM / xSD )
b3.m = z3.m * ySD / xSD
b2.m = ( z2.m * ySD - z3.m * ySD * xM / xSD )
s1.m = (zs1.m * ySD)^2
s2.m = (zs2.m * ySD)^2

## Posterior means and 95% credible intervals under multimodel inference
# Intercept
mean(b0.m)
quantile(b0.m, probs=c(0.025, 0.5, 0.975))
# Slope
mean(b1.m)
quantile(b1.m, probs=c(0.025, 0.5, 0.975))
# Group effect
mean(b2.m)
quantile(b2.m, probs=c(0.025, 0.5, 0.975))
# Interaction
mean(b3.m)
quantile(b3.m, probs=c(0.025, 0.5, 0.975))
# Plesiomorphic rate
mean(s1.m)
quantile(s1.m, probs=c(0.025, 0.5, 0.975))
# Apomorphic rate
mean(s2.m)
quantile(s2.m, probs=c(0.025, 0.5, 0.975))
# Phylogenetic signal
mean(lambda.m)
median(lambda.m)           # Posterior median
quantile(lambda.m, probs=c(0.025, 0.5, 0.975))


################################################################################

# Analyze the same model under a single rate for comparison

modelstring = "
model {
#Linear regression and multivariate normal likelihood
  for (i in 1:Ndata) {
    mu[i] <- beta0 + beta1 * x[i] + beta2 * z[i] + beta3 * xz[i]
  }
  y[1:Ndata] ~ dmnorm(mu[],TAU[,])
  
#Priors for regression coefficients (betas) and rates (gammas)
  beta0 ~ dnorm( 0 , 1.0E-06 )
  beta1 ~ dnorm( 0 , 1.0E-06 )
  beta2 ~ dnorm( 0 , 1.0E-06 )
  beta3 ~ dnorm( 0 , 1.0E-06 )
  sigma ~ dunif( 0 , 100 )
  gamma <- pow( sigma , 2 )

#Equal vector of probability for tree sampling
  for (k in 1:Ntree) {
    p[k] <- 1/Ntree
  }

#Tree sampling and variance-covariance matrix construction
  K ~ dcat(p[])
  PS <- LAMBDA*A[,,K]+(1-LAMBDA)*ID
  TAU <- inverse(gamma*PS)

# Mixture model:
  LAMBDA <- lambdaModel[ modelIndex ]
  lambdaModel[1] <- 0
  lambdaModel[2] <- 1
  lambdaModel[3] ~ dunif(0.01,0.99)
# Hyperprior:
  modelIndex ~ dcat( modelProb[] )
  modelProb[1] <- 1/3
  modelProb[2] <- 1/3
  modelProb[3] <- 1/3
}
"
writeLines(modelstring,con="single.txt")


# Specify phylogenetic similarity matrices
G=array(NA,dim=c(nSubj,nSubj))
A=array(NA, dim=c(nSubj,nSubj,Ntree))
for(i in 1:Ntree){
  G = cophenetic(mtrees[[i]])
  G = max(G) - G 
  G = G/max(G) 
  G <- as.matrix(G) 
  G <- G[order(dimnames(G)[[1]], na.last=NA) , ]
  G <- G[,sort(colnames(G))]
  A[,,i]<-G
}

# Construct new list of data
dataList = list(
  x = zx ,
  y = zy ,
  z = z ,
  xz = zx*z ,
  Ndata = nSubj,
  Ntree = Ntree,
  A = A,
  ID = diag(nSubj)
)

# Initialize chains based on data standardization
r = cor(x,y)
initsList = list(
    beta0 = 0 ,
    beta1 = r ,
    beta2 = 0 ,
    beta3 = 0 ,
    sigma = 1-r^2
)

# Set up chains
parameters = c("beta0" , "beta1" , "beta2"  , "beta3"  , "sigma", "LAMBDA", "modelIndex", "K")  # Saved parameters
adaptSteps = 1000             # Tunning steps
burnInSteps = 10000           # Burn-in steps
nChains = 3                   # Number of chains
numSavedSteps=100000          # Total steps
thinSteps=1                   # Thinning
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain
# Initialize model
jagsSingle = jags.model( "single.txt" , data=dataList , inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in
update( jagsSingle , n.iter=burnInSteps )
# Sample final MCMC chain
SamplSingle = coda.samples( jagsSingle , variable.names=parameters , 
                            n.iter=nPerChain , thin=thinSteps )

# Examine parameters under standardized data
summary(SamplSingle)

# Convert coda-object to facilitate handling, and extract chain values
mcmcChain.s = as.matrix( SamplSingle )
z0.s = mcmcChain.s[, "beta0" ]
z1.s = mcmcChain.s[, "beta1" ]
z2.s = mcmcChain.s[, "beta2" ]
z3.s = mcmcChain.s[, "beta3" ]
zs.s = mcmcChain.s[, "sigma" ]
lambda.s = mcmcChain.s[, "LAMBDA" ]

# Get the posterior sample of modelIndex
modelIdxSample.s = mcmcChain.s[,"modelIndex"]
# Compute probability of each model
pIID.s = sum( modelIdxSample.s == 1 ) / length( modelIdxSample.s )
pBM.s = sum( modelIdxSample.s == 2 ) / length( modelIdxSample.s )
pLM.s = sum( modelIdxSample.s == 3 ) / length( modelIdxSample.s )

# Convert to original scale
b1.s = z1.s * ySD / xSD
b0.s = ( z0.s * ySD + yM - z1.s * ySD * xM / xSD )
b3.s = z3.s * ySD / xSD
b2.s = ( z2.s * ySD - z3.s * ySD * xM / xSD )
sg.s = (zs.s * ySD)^2

## Posterior means and 95% credible intervals under multimodel inference
# Intercept
mean(b0.s)
quantile(b0.s, probs=c(0.025, 0.5, 0.975))
# Slope
mean(b1.s)
quantile(b1.s, probs=c(0.025, 0.5, 0.975))
# Group effect
mean(b2.s)
quantile(b2.s, probs=c(0.025, 0.5, 0.975))
# Interaction
mean(b3.s)
quantile(b3.s, probs=c(0.025, 0.5, 0.975))
# Rate
mean(sg.s)
quantile(sg.s, probs=c(0.025, 0.5, 0.975))
# Phylogenetic signal
mean(lambda.s)
median(lambda.s)           # Posterior median
quantile(lambda.s, probs=c(0.025, 0.5, 0.975))






















