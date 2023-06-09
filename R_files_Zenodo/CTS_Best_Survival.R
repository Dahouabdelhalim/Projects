###California Tiger Salamander (Ambystoma californiense) Cormack-Jolly-Seber Model
###Stage-structured (metamorph or adult on first capture) individual based model with body mass covariate
###Best model of set
###Arianne Messerman
###March 2021

#setwd("")

# Load data
M <- read.table("metamorphs-v2.txt",header=T,colClasses="character")
head(M)
M$meta<-as.factor(M$meta)
M$adult<-as.factor(M$adult)
str(M) #4036 metamorphs

A <- read.table("adults-v2.txt",header=T,colClasses="character")
head(A)
A$meta<-as.factor(A$meta)
A$adult<-as.factor(A$adult)
str(A) #5309 adults

cts <-read.table("covariates-v2.txt",header=T,colClasses="character")
head(cts)
cts$meta<-as.factor(cts$meta)
cts$adult<-as.factor(cts$adult)
str(cts) #9345 salamanders total

####Build txt file capture history into matrix 
CH.M<-t(array(as.numeric(unlist(strsplit(M$history,""))),dim=c(nchar(M$history[1]),length(M$meta))))
CH.M[1:10,]
CH.A<-t(array(as.numeric(unlist(strsplit(A$history,""))),dim=c(nchar(A$history[1]),length(A$meta))))
CH.A[1:10,]
CH1<-t(array(as.numeric(unlist(strsplit(cts$history,""))),dim=c(nchar(cts$history[1]),length(cts$meta))))
CH1[1:10,]
str(CH1)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f.m <- apply(CH.M, 1, get.first)
str(f.m)
f.a <- apply(CH.A, 1, get.first)
unique(f.a)
str(f.a)

###Build stage class matrices
x.m <- matrix(NA, ncol=dim(CH.M)[2]-1,nrow=dim(CH.M)[1])
str(x.m)

for (i in 1:dim(CH.M)[1]){
  for(t in f.m[i]:(dim(CH.M)[2]-1)){
    x.m[i,t] <- 2
    x.m[i, f.m[i]] <- 1
  }#t
} #i
x.m[931:993,] #Check if matches expectation based on f.m
str(x.m)


x.a <- matrix(NA, ncol=dim(CH.A)[2],nrow=dim(CH.A)[1]) #Must initially have 9 columns to accommodate individuals caught in year 9
str(x.a)

for (i in 1:dim(CH.A)[1]){
  for(t in f.a[i]:(dim(CH.A)[2])){
    x.a[i,t] <- 2
  }#t
} #i
x.a <- x.a[1:nrow(x.a),(1:ncol(x.a)-1)] #Delete 9th column

x.a[249:311,]#Check if matches expectation based on f.a
str(x.a)


##Merge capture histories (CH),first marking vectors (f), and stage matrices (x)
CH <- rbind(CH.M, CH.A)
f <- c(f.m, f.a)
x <- rbind(x.m, x.a)
str(x)

all<-rbind(M,A)
str(all)
str(cts)
test1<-(as.numeric(all$distance)-as.numeric(cts$distance))
sum(test1) #=0, so cts and all dataframes match perfectly
test2<-(as.numeric(all$ln.mass)-as.numeric(cts$ln.mass))
sum(test2) #=0, so cts and all dataframes match perfectly
test3<-(as.numeric(all$history)-as.numeric(cts$history))
sum(test3) #=0, so cts and all dataframes match perfectly

tail(x)
tail(CH)
tail(f)

library(jagsUI)
library(mcmcplots)

set.seed(423467)

cts$ln.mass <- cts$cntmass <- as.numeric(cts$ln.mass)
lnmass<-cts$ln.mass
stdmass<-cntmass<-rep(NA,length(cts$ln.mass))

#Scale ln(body mass) individual covariate
for (i in 1:length(cts$ln.mass)) {
  stdmass[i] <- (cts$ln.mass[i]-mean(cts$ln.mass[]))/sd(cts$ln.mass[])
  cntmass[i] <- (cts$ln.mass[i]-mean(cts$ln.mass))
  cts$cntmass[i] <- (cts$ln.mass[i]-mean(cts$ln.mass))
}

mean(cntmass)#check = 0
x.all<-mean(cts$ln.mass)


#############################################################
# 35. Phi(t)P(t): Model with random time-dependent metamorph survival and adult recapture
# Constant adult survival and metamorph recapture
# With fixed stage effect on survival AND recapture
# Metamorph survival only estimated in mark years
# With body mass covariate on survival only
# Using mean-centered body mass
#############################################################
### Specify model
sink("cts-cjs-mixed-stage-14center.jags")
cat("
model {
    
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    phi[i,t] <- (1/(1+exp(-(eta.phi[x[i,t], t] + alpha*mass[i]))))   # Survival
    p[i,t] <- (1/(1+exp(-(eta.p[x[i,t], t]))))      # Recapture
    } #t
    } #i
    
    for (t in 1:(n.occasions-1)){
    phi1.est[t]<-1/(1+exp(-eta.phi[1,t] - alpha*-0.124))
    phi2.est[t]<-1/(1+exp(-eta.phi[2,t] - alpha*0.094))
    p1.est[t]<-1/(1+exp(-eta.p[1,t]))
    p2.est[t]<-1/(1+exp(-eta.p[2,t]))
    } #Back transformed stage-specific mean phi and p estimates using stage-specific mean mass
    
    #Back transformed stage-specific mean metamorph phi estimates using mean centered mass of years 1,2,6,7
    phi1.est1<-1/(1+exp(-eta.phi[1,1] - alpha*-0.09))
    phi1.est2<-1/(1+exp(-eta.phi[1,2] - alpha*-0.38))
    phi1.est6<-1/(1+exp(-eta.phi[1,6] - alpha*0.14))
    phi1.est7<-1/(1+exp(-eta.phi[1,7] - alpha*-0.12))
    
    alpha ~ dnorm(0, 0.001)I(-5, 5)     # Prior for body mass slope on survival
    
    for (t in 1:(n.occasions-1)){
      eta.phi[1,t] <- mu.phi[1]  + epsilon.phi[t] #Random time-dependent metamorph survival
      eta.phi[2,t] <- mu.phi[2]
      epsilon.phi[t] ~ dnorm(0, tau.phi)T(-15, 15)
      eta.p[1,t] <- mu.p[1]
      eta.p[2,t] <- mu.p[2] + epsilon.p[t]     #Random time-dependent adult recapture
      epsilon.p[t] ~ dnorm(0, tau.p)T(-15, 15)
    } #t
    
    for (u in 1:2){
    mean.phi[u] ~ dunif(0, 1)       #Priors on mean stage-specific survival
    mu.phi[u] <- log(mean.phi[u] / (1-mean.phi[u]))  #Transformed
    mean.p[u] ~ dunif(0, 1)       #Priors on mean stage-specific recapture
    mu.p[u] <- log(mean.p[u] / (1-mean.p[u]))  #Transformed
    }
    
    tau.phi ~ dgamma(1, 0.1)T(0.000001,200)
    sigma.phi <- 1/sqrt(tau.phi) #Priors for stage-specific sd on survival
    sigma2.phi <- pow(sigma.phi, 2)
    tau.p ~ dgamma(1, 0.1)T(0.000001,200)
    sigma.p <- 1/sqrt(tau.p) #Priors for stage-specific sd on survival
    sigma2.p <- pow(sigma.p, 2)
  
  
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture 
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- phi[i,t-1] * z[i,t-1]
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Function to create a matrix with information about known latent state z
known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],  z = known.state.cjs(CH), x=x, mass = cntmass)

# Initial values
inits <- function(){list(mean.phi = runif(2,0,1), alpha = runif(1, -3, 3),
                         mean.p = runif(2,0,1), z = cjs.init.z(CH,f))}

# Parameters monitored
parameters <- c("mean.phi", "mu.phi", "alpha", "sigma.phi", "mean.p", "mu.p", "sigma.p", "phi1.est", "phi2.est", 
                "phi1.est1", "phi1.est2", "phi1.est6", "phi1.est7","p1.est", "p2.est","eta.phi", "eta.p")

# MCMC settings
ni <- 700000
nt <- 5
nb <- 350000
nc <- 3

# Call JAGS from R
cts.cjs.mixed.stage.14center <- jags(jags.data, inits, parallel=TRUE, parameters, "cts-cjs-mixed-stage-14center.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

#save("/cts-cjs-mixed-stage-14center.Rdata")
#load("/cts-cjs-mixed-stage-14center.Rdata")

# Summarize posteriors, run time 3355.29 min
print(cts.cjs.mixed.stage.14center, digits = 3)

plot(cts.cjs.mixed.stage.14center)

cts.cjs.mixed.stage.14center$mean$deviance
cts.cjs.mixed.stage.14center$q2.5$deviance
cts.cjs.mixed.stage.14center$q97.5$deviance

library(dplyr)

#Calculate temporal residual variance from std deviations
sigma2.phi <- cts.cjs.mixed.stage.14center$sims.list$sigma.phi^2
mean(sigma2.phi)#6.24
cts.cjs.mixed.stage.14center$mean$sigma.phi^2#4.21
cts.cjs.mixed.stage.14center$q2.5$sigma.phi^2#0.78
cts.cjs.mixed.stage.14center$q97.5$sigma.phi^2#22.32

sigma2.p <- cts.cjs.mixed.stage.14center$sims.list$sigma.p^2
mean(sigma2.p)#0.33
cts.cjs.mixed.stage.14center$mean$sigma.p^2#0.30
cts.cjs.mixed.stage.14center$q2.5$sigma.p^2#0.09
cts.cjs.mixed.stage.14center$q97.5$sigma.p^2#1.00

#Difference in life stage-specific posterior distributions for phi and p
surv.dif<-cts.cjs.mixed.stage.14center$sims.list$mean.phi[,1]-cts.cjs.mixed.stage.14center$sims.list$mean.phi[,2]
plot(density(surv.dif))
recap.dif<-cts.cjs.mixed.stage.14center$sims.list$mean.p[,1]-cts.cjs.mixed.stage.14center$sims.list$mean.p[,2]
plot(density(recap.dif))


#Find mean metamorph survival and credible intervals for only years 1,2,6, and 7 when individuals were marked
#Using cohort-specific mean masses
phi.sims<-c(cts.cjs.mixed.stage.14center$sims.list$phi1.est1, cts.cjs.mixed.stage.14center$sims.list$phi1.est2, cts.cjs.mixed.stage.14center$sims.list$phi1.es6, cts.cjs.mixed.stage.14center$sims.list$phi1.est7)
mean(phi.sims)#0.44
median(phi.sims)#0.54
quantile(phi.sims, c(0.025, 0.975))#q2.5=0.03, q97.5=0.90

#Using mean metamorph mass overall
phi.simsx<-c(cts.cjs.mixed.stage.14center$sims.list$phi1.est[,1], cts.cjs.mixed.stage.14center$sims.list$phi1.est[,2], cts.cjs.mixed.stage.14center$sims.list$phi1.est[,6], cts.cjs.mixed.stage.14center$sims.list$phi1.est[,7])
mean(phi.simsx)#0.55
median(phi.simsx)#0.61
quantile(phi.simsx, c(0.025, 0.975))#q2.5=0.04, q97.5=0.99


#Pull 500 sets of iterations from posterior of coefficients to .CSV
mcmc.sample<-cts.cjs.mixed.stage.14center$mcmc.info$n.samples
sub.set <- sort(sample(1:mcmc.sample, size=500))

##Extranct random samples for export to .CSV file
#For survival
meta.mu.surv <- ad.mu.surv <- meta.mean.surv <- ad.mean.surv <- mass.alpha <- 
  mean.meta.eta.phi<- sigma.phi <- phi1 <- phi1.all <- phi2 <- c(1:length(sub.set))
meta.eta.surv <-ad.eta.surv<- epsilon.surv <- matrix(NA, nrow=length(sub.set), ncol=8)
meta.eta.phi1<-meta.eta.phi2<-meta.eta.phi6<-meta.eta.phi7 <-c(1:length(sub.set))
phi.1<-phi.2<-phi.3<-phi.4<-phi.5<-phi.6<-phi.7<-phi.8<-c(1:length(sub.set))
for(i in 1:length(sub.set)){
  meta.mu.surv[i] <- cts.cjs.mixed.stage.14center$sims.list$mu.phi[sub.set[i],1]
  ad.mu.surv[i] <- cts.cjs.mixed.stage.14center$sims.list$mu.phi[sub.set[i],2]
  meta.mean.surv[i] <- cts.cjs.mixed.stage.14center$sims.list$mean.phi[sub.set[i],1]
  ad.mean.surv[i] <- cts.cjs.mixed.stage.14center$sims.list$mean.phi[sub.set[i],2]
  mass.alpha[i] <- cts.cjs.mixed.stage.14center$sims.list$alpha[sub.set[i]]
  meta.eta.phi1[i]<-cts.cjs.mixed.stage.14center$sims.list$eta.phi[sub.set[i],1,1]
  meta.eta.phi2[i]<-cts.cjs.mixed.stage.14center$sims.list$eta.phi[sub.set[i],1,2]
  meta.eta.phi6[i]<-cts.cjs.mixed.stage.14center$sims.list$eta.phi[sub.set[i],1,6]
  meta.eta.phi7[i]<-cts.cjs.mixed.stage.14center$sims.list$eta.phi[sub.set[i],1,7]
  mean.meta.eta.phi[i]<-mean(c(meta.eta.phi1[i], meta.eta.phi2[i], meta.eta.phi6[i], meta.eta.phi7[i]))
  phi.1[i] <- cts.cjs.mixed.stage.14center$sims.list$phi1.est[sub.set[i],1]
  phi.2[i] <- cts.cjs.mixed.stage.14center$sims.list$phi1.est[sub.set[i],2]
  phi.3[i] <- cts.cjs.mixed.stage.14center$sims.list$phi1.est[sub.set[i],3]
  phi.4[i] <- cts.cjs.mixed.stage.14center$sims.list$phi1.est[sub.set[i],4]
  phi.5[i] <- cts.cjs.mixed.stage.14center$sims.list$phi1.est[sub.set[i],5]
  phi.6[i] <- cts.cjs.mixed.stage.14center$sims.list$phi1.est[sub.set[i],6]
  phi.7[i] <- cts.cjs.mixed.stage.14center$sims.list$phi1.est[sub.set[i],7]
  phi.8[i] <- cts.cjs.mixed.stage.14center$sims.list$phi1.est[sub.set[i],8]
  phi2[i] <- cts.cjs.mixed.stage.14center$sims.list$phi2.est[sub.set[i],1]
  mean(cts.cjs.mixed.stage.14center$mean$phi2.est)#Mean adult phi across all years at mean ln(mass)=0.52
  sigma.phi[i] <- cts.cjs.mixed.stage.14center$sims.list$sigma.phi[sub.set[i]]
}

phis<-data.frame(phi.1, phi.2, phi.3, phi.4, phi.5, phi.6, phi.7)
phis.sub<-data.frame(phi.1, phi.2, phi.6, phi.7)
meta.etas<-data.frame(meta.eta.phi1, meta.eta.phi2, meta.eta.phi6, meta.eta.phi7)
mean(as.matrix(phis))#mean across whole data frame=0.53
mean(as.matrix(phi1))#mean across whole data frame=0.55
phi1.all<-rowMeans(phis)#mean of means=0.53 across all study years for metamorphs
phi1<-rowMeans(phis.sub)#mean of means=0.55 across only years with data for metamorphs
mean.meta.eta.phi<-rowMeans(meta.etas)#mean=0.36

phi.sims1<-c(phi.1, phi.2, phi.6, phi.7)
mean(phi.sims1)#0.55
median(phi.sims)#0.54
quantile(phi.sims1, c(0.025, 0.975))#q2.5=0.04, q97.5=0.99

for (t in 1:8){
  for(i in 1:length(sub.set)){
    meta.eta.surv[i,t] <- cts.cjs.mixed.stage.14center$sims.list$eta.phi[sub.set[i],1,t]
    ad.eta.surv[i,t] <- cts.cjs.mixed.stage.14center$sims.list$eta.phi[sub.set[i],2,t]
  }
}

#For recapture
meta.mu.recap <- ad.mu.recap <- meta.mean.recap <- ad.mean.recap <- sigma.p <- c(1:length(sub.set))
meta.eta.recap <-ad.eta.recap<- epsilon.recap <- matrix(NA, nrow=length(sub.set), ncol=8)

for(i in 1:length(sub.set)){
  meta.mu.recap[i] <- cts.cjs.mixed.stage.14center$sims.list$mu.p[sub.set[i],1]
  ad.mu.recap[i] <- cts.cjs.mixed.stage.14center$sims.list$mu.p[sub.set[i],2]
  meta.mean.recap[i] <- cts.cjs.mixed.stage.14center$sims.list$mean.p[sub.set[i],1]
  ad.mean.recap[i] <- cts.cjs.mixed.stage.14center$sims.list$mean.p[sub.set[i],2]
  sigma.p[i] <- cts.cjs.mixed.stage.14center$sims.list$sigma.p[sub.set[i]]
}

for (t in 1:8){
  for(i in 1:length(sub.set)){
    meta.eta.recap[i,t] <- cts.cjs.mixed.stage.14center$sims.list$eta.p[sub.set[i],1,t]
    ad.eta.recap[i,t] <- cts.cjs.mixed.stage.14center$sims.list$eta.p[sub.set[i],2,t]
  }
}

##Merge samples into a single dataframe
ests<-data.frame(cbind(meta.mu.surv, ad.mu.surv, meta.mean.surv, ad.mean.surv, mass.alpha, mean.meta.eta.phi, meta.eta.surv, ad.eta.surv, sigma.phi, phi1, phi1.all, phi2, #epsilon.surv,
                       meta.mu.recap, ad.mu.recap, meta.mean.recap, ad.mean.recap, meta.eta.recap, ad.eta.recap, sigma.p)) #epsilon.recap))
names(ests)<-c("meta.mu.phi","ad.mu.phi", "meta.mean.phi","ad.mean.phi", "mass.alpha", "mean.meta.eta.phi","meta.eta.phi1", 
               "meta.eta.phi2", "meta.eta.phi3", "meta.eta.phi4", "meta.eta.phi5", "meta.eta.phi6", "meta.eta.phi7", 
               "meta.eta.phi8", "ad.eta.phi1", "ad.eta.phi2", "ad.eta.phi3", "ad.eta.phi4", 
               "ad.eta.phi5", "ad.eta.phi6", "ad.eta.phi7", "ad.eta.phi8", "sigma.phi", "phi1", "phi1.all", "phi2",
               "meta.mu.recap", "ad.mu.recap", "meta.mean.recap", "ad.mean.recap", "meta.eta.recap1", 
               "meta.eta.recap2", "meta.eta.recap3", "meta.eta.recap4", "meta.eta.recap5", "meta.eta.recap6", 
               "meta.eta.recap7", "meta.eta.recap8", "ad.eta.recap1", "ad.eta.recap2", "ad.eta.recap3", "ad.eta.recap4", 
               "ad.eta.recap5", "ad.eta.recap6", "ad.eta.recap7", "ad.eta.recap8", "sigma.p") #"epsilon.recap1", 
str(ests)

#Across all iterations of model output:
mean(cts.cjs.mixed.stage.14center$mean$phi1.est)#Mean metamorph phi across all years centered at mean ln(mass)=0.55
mean(cts.cjs.mixed.stage.14center$mean$phi2.est)#Mean adult phi across all years centered at mean ln(mass)=0.54
mean((1/(1+exp(-ests$mean.meta.eta.phi))))#Mean metamorph eta.phi (intercept) for years with centered data = 0.58

#Create random normal distributions of emigration rates with mean and SE from Searcy et al.
#Use SE instead of SD because estimating across individuals, not for each individual
sd.m<-0.00946*sqrt(206)#Calculate SD from SE for 206 metamorphs, SD=0.135
sd.a<-0.0114*sqrt(69)#SD=0.095 for 69 adults
emig.m<-rnorm(mcmc.sample, 0.257, 0.00946)#metamorph emigration
emig.a<-rnorm(mcmc.sample, 0.148, 0.0114)#adult emigration
hist(emig.m)
hist(emig.a)

#Draw random samples from this distribution and add to the ests dataframe
ests$emig.m<-sample(emig.m, length(sub.set), replace=TRUE)
ests$emig.a<-sample(emig.a, length(sub.set), replace=TRUE)
hist(ests$emig.m)
hist(ests$emig.a)
mean(ests$emig.m)#0.26
mean(ests$emig.a)#0.15

#Calculate site fidelity from emigration 
ests$fidel.m<-1-ests$emig.m
ests$fidel.a<-1-ests$emig.a
hist(ests$fidel.m)
hist(ests$fidel.a)
mean(ests$fidel.m)#0.74
mean(ests$fidel.a)#0.85
sd(ests$fidel.m)#0.009
sd(ests$fidel.a)#0.011

#Calculate s (true survival estimate) from phi (apparent survival) and fidelity
ests$surv.m<-ests$phi1/ests$fidel.m
ests$surv.a<-ests$phi2/ests$fidel.a
hist(ests$surv.m)
hist(ests$surv.a)

mean(phi1)#0.55
mean(ests$surv.m)#0.74
mean(phi2)#0.52
mean(ests$surv.a)#0.60


#Multiply metamorph survival function by metamorph body size in years 1,2,6,7
mean(as.numeric(A$ln.mass))#mean adult ln(mass)=2.494
mean(as.numeric(M$ln.mass))#mean metamorph ln(mass)=2.276

M$first<-as.factor(f.m)
M$ln.mass<-M$cntmass<-as.numeric(M$ln.mass)
A$first<-as.factor(f.a)
A$ln.mass<-A$cntmass<-as.numeric(A$ln.mass)

library(dplyr)

for (i in 1:length(M$ln.mass)){
  M$cntmass[i] <- (M$ln.mass[i]-mean(cts$ln.mass))
}

for (i in 1:length(A$ln.mass)){
  A$cntmass[i] <- (A$ln.mass[i]-mean(cts$ln.mass))
}

m.yr<-M %>%
  group_by(first) %>%
  summarize(means=mean(cntmass))
m.yr<-m.yr[order(m.yr$means),]
m.yr

m.1<-M[which(M$first=='1'),]
m.1<-m.1[order(m.1$cntmass),]
m.2<-M[which(M$first=='2'),]
m.2<-m.2[order(m.2$cntmass),]
m.6<-M[which(M$first=='6'),]
m.6<-m.6[order(m.6$cntmass),]
m.7<-M[which(M$first=='7'),]
m.7<-m.7[order(m.7$cntmass),]

#Vector of back-transformed estimated phi given year-specific mean-centered individual body masses in years 1, 2, 6, and 7
m.1$phim.est<-1/(1+exp(-cts.cjs.mixed.stage.14center$mean$eta.phi[1,1] - cts.cjs.mixed.stage.14center$mean$alpha*m.1$cntmass))
m.2$phim.est<-1/(1+exp(-cts.cjs.mixed.stage.14center$mean$eta.phi[1,2] - cts.cjs.mixed.stage.14center$mean$alpha*m.2$cntmass))
m.6$phim.est<-1/(1+exp(-cts.cjs.mixed.stage.14center$mean$eta.phi[1,6] - cts.cjs.mixed.stage.14center$mean$alpha*m.6$cntmass))
m.7$phim.est<-1/(1+exp(-cts.cjs.mixed.stage.14center$mean$eta.phi[1,7] - cts.cjs.mixed.stage.14center$mean$alpha*m.7$cntmass))

#Create dataset of metamorph body sizes
m.mass<-sort(sample(M$cntmass, size=length(sub.set), replace=F))
m.mass<-m.mass[order(m.mass)]
phim.avg<-1/(1+exp(-mean(mean.meta.eta.phi) - mean(mass.alpha)*m.mass))

A<-A[order(A$cntmass),]
A$phi.cnt<-1/(1+exp(-cts.cjs.mixed.stage.14center$mean$eta.phi[2,1] - cts.cjs.mixed.stage.14center$mean$alpha*A$cntmass))

mean(m.1$phim.est)#0.72
mean(m.2$phim.est)#0.06
mean(m.6$phim.est)#0.95
mean(m.7$phim.est)#0.57
mean(c(m.1$phim.est, m.2$phim.est, m.6$phim.est, m.7$phim.est))#0.55 from individual masses
mean(phim.avg)#0.57 when considering mean cohort mass
mean(A$phi.cnt)#0.51

par(mfrow=c(1,1))
plot(m.1$cntmass, m.1$phim.est, ylim=c(0,1), xlim=c(min(cntmass),max(cntmass)), typ='l', lty=4, lwd=3,
     ylab="Apparent survival probability", xlab="Mean centered ln(body mass (g))", bty = 'l')
lines(m.2$cntmass, m.2$phim.est, col=2, lty=2, lwd=3)
lines(m.6$cntmass, m.6$phim.est, col=3, lty=3, lwd=3)
lines(m.7$cntmass, m.7$phim.est, col=4, lty=1, lwd=3)
lines(m.mass, phim.avg, col=6, lty=5, lwd=3)
lines(A$cntmass, A$phi.cnt, col=5, lty=6, lwd=3)
legend(0.7,0.3, legend=c("metamorphs 2005", "metamorphs 2006", "metamorphs 2010", "metamorphs 2011", "metamorphs all years", "juveniles/adults all years"), 
       col=c(1, 2, 3, 4, 6, 5), lty=c(4, 2, 3, 1, 5, 6), lwd=3)

#Vector of back-transformed upper and lower phi credible intervalsgiven year-specific mean-centered individual body masses in years 1, 2, 6, and 7
m.1$phim.low<-1/(1+exp(-cts.cjs.mixed.stage.14center$q2.5$eta.phi[1,1] - cts.cjs.mixed.stage.14center$q2.5$alpha*m.1$cntmass))
m.2$phim.low<-1/(1+exp(-cts.cjs.mixed.stage.14center$q2.5$eta.phi[1,2] - cts.cjs.mixed.stage.14center$q2.5$alpha*m.2$cntmass))
m.6$phim.low<-1/(1+exp(-cts.cjs.mixed.stage.14center$q2.5$eta.phi[1,6] - cts.cjs.mixed.stage.14center$q2.5$alpha*m.6$cntmass))
m.7$phim.low<-1/(1+exp(-cts.cjs.mixed.stage.14center$q2.5$eta.phi[1,7] - cts.cjs.mixed.stage.14center$q2.5$alpha*m.7$cntmass))
A$phi.low<-1/(1+exp(-cts.cjs.mixed.stage.14center$q2.5$eta.phi[2,1] - cts.cjs.mixed.stage.14center$q2.5$alpha*A$cntmass))

m.1$phim.high<-1/(1+exp(-cts.cjs.mixed.stage.14center$q97.5$eta.phi[1,1] - cts.cjs.mixed.stage.14center$q97.5$alpha*m.1$cntmass))
m.2$phim.high<-1/(1+exp(-cts.cjs.mixed.stage.14center$q97.5$eta.phi[1,2] - cts.cjs.mixed.stage.14center$q97.5$alpha*m.2$cntmass))
m.6$phim.high<-1/(1+exp(-cts.cjs.mixed.stage.14center$q97.5$eta.phi[1,6] - cts.cjs.mixed.stage.14center$q97.5$alpha*m.6$cntmass))
m.7$phim.high<-1/(1+exp(-cts.cjs.mixed.stage.14center$q97.5$eta.phi[1,7] - cts.cjs.mixed.stage.14center$q97.5$alpha*m.7$cntmass))
A$phi.high<-1/(1+exp(-cts.cjs.mixed.stage.14center$q97.5$eta.phi[2,1] - cts.cjs.mixed.stage.14center$q97.5$alpha*A$cntmass))

library(ggplot2)
library(cowplot)

Class<-rep.int(c("metamorphs 2005","metamorphs 2006","metamorphs 2010","metamorphs 2011", "adults all years"), 
           times=c(nrow(m.1), nrow(m.2), nrow(m.6), nrow(m.7), nrow(A)))#"metamorphs all years", length(phim.avg)
masses<-c(m.1$ln.mass, m.2$ln.mass, m.6$ln.mass, m.7$ln.mass, A$ln.mass)#m.mass
phi.proj<-c(m.1$phim.est, m.2$phim.est, m.6$phim.est, m.7$phim.est, A$phi.cnt)#phim.avg
lower<-c(m.1$phim.low, m.2$phim.low, m.6$phim.low, m.7$phim.low, A$phi.low)
upper<-c(m.1$phim.high, m.2$phim.high, m.6$phim.high, m.7$phim.high, A$phi.high)
preds<-as.data.frame(Class)
preds$mass<-masses
preds$pred<-phi.proj
preds$low<-lower
preds$high<-upper

ggplot(preds, aes(mass, pred, ymin=low, ymax=high, fill=Class, colour=Class, linetype=Class))+
  geom_ribbon(alpha=0.2, colour=NA)+
  geom_smooth()+
  theme_cowplot()+
  labs(y = "Apparent survival probability", x = "Ln(body mass (g))")+
  theme(legend.position=c(0.8,0.2))


write.csv(ests,"survival-posterior-samples-CENTERED.csv")
