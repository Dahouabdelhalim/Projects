########################################
### Capture-recapture models for the effect of inbreeding on survival in Helgeland house sparrows
### Original code written by HÃ¥kon Holand
###
### April 2020 small changes by Alina Niskanen
########################################



library(jagsUI)

chh<-read.table("ch.csv",sep="|", header=F)
pd<-read.table("density.csv",sep="|", header=F)
is<-scan("island.csv")
fg<-scan("fgrm.csv")
fh<-scan("froh.csv")
io<-scan("island_status.csv")
sex<-scan("gen_sex.csv")

sex<-sex-1
sex[sex==0]<-2
sex

pd<-as.matrix(pd)
pd<-pd-mean(pd,na.rm = T)
pd[is.na(pd)]<-0

CH<-as.matrix(chh)
rm(chh)

get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)


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


#full model for FGRM including habitat type interaction with F-------------------------------------------------------------------
# Specify model in BUGS language
sink("cjs-temp-raneff.bug")
cat("
    model {

    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu + habitat[io[i]] + kjonn[sex[i]] + beta*fg[i] + beta.pd*pd[i,t] + beta.sex[sex[i]]*fg[i] + beta.hab[io[i]]*fg[i]+ beta.dens*pd[i,t]*fg[i] + epsilon[t]+ epsilon.is[is[i]]
    p[i,t] <- popaar[is[i],t]
    } #t
    } #i
    for (t in 1:(n.occasions-1)){
    epsilon[t] ~ dnorm(0, tau)
    }
    sigma ~ dunif(0, 10)                     # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)                  # Temporal variance
    
    mu ~ dnorm(0, 0.001)I(-5,5)              # Prior for logit of mean survival

    for (h in 1:8){
    epsilon.is[h] ~ dnorm(0, tau.is)
    }
    
    sigma.is ~ dunif(0, 10)                     # Prior for standard deviation
    tau.is <- pow(sigma.is, -2)
    sigma2.is <- pow(sigma.is, 2)                  # Temporal variance
    
    habitat[1]<-0
    habitat[2]~ dnorm(0, 0.001)I(-5,5)

    kjonn[1]<-0
    kjonn[2]~ dnorm(0, 0.001)I(-5,5)

    beta~ dnorm(0, 0.001)I(-10,10)
    beta.pd~ dnorm(0, 0.001)I(-10,10)
    beta.sex[1]<-0
    beta.sex[2]~ dnorm(0, 0.001)I(-10,10)

    beta.hab[1]<-0
    beta.hab[2]~ dnorm(0, 0.001)I(-10,10)

    beta.dens~ dnorm(0, 0.001)I(-10,10)
    
    for (u in 1:8){
    for (h in 1:(n.occasions-1)){
    popaar[u,h]~ dunif(0, 1)
    }
    }
    
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

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH),fg=fg,is=is,io=io,sex=sex,pd=pd)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), sigma = runif(1, 0, 10),sigma.is = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mu","habitat","kjonn","beta","beta.pd","beta.sex","beta.hab","beta.dens","sigma2","sigma2.is","popaar")

# MCMC settings
ni <- 150000
nt <- 8
nb <- 110000
nc <- 3


system.time(mod.1 <- jags(bugs.data, inits, parameters, "cjs-temp-raneff.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE))

# Summarize posteriors
print(mod.1, digits = 3)

save(mod.1,file="mod1rerun.rda")




#final model for FGRM-------------------------------------------------------------------
# Specify model in BUGS language
sink("cjs-temp-raneff.bug")
cat("
    model {

    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu + kjonn[sex[i]] + beta*fg[i] + epsilon[t]+ epsilon.is[is[i]]
    p[i,t] <- popaar[is[i],t]
    } #t
    } #i
    for (t in 1:(n.occasions-1)){
    epsilon[t] ~ dnorm(0, tau)
    }
    sigma ~ dunif(0, 10)                     # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)                  # Temporal variance
    
    mu ~ dnorm(0, 0.001)I(-5,5)              # Prior for logit of mean survival
    
    for (h in 1:8){
    epsilon.is[h] ~ dnorm(0, tau.is)
    }
    
    sigma.is ~ dunif(0, 10)                     # Prior for standard deviation
    tau.is <- pow(sigma.is, -2)
    sigma2.is <- pow(sigma.is, 2)                  # Temporal variance
    
    kjonn[1]<-0
    kjonn[2]~ dnorm(0, 0.001)I(-5,5)
    
    beta~ dnorm(0, 0.001)I(-10,10)
    
    for (u in 1:8){
    for (h in 1:(n.occasions-1)){
    popaar[u,h]~ dunif(0, 1)
    }
    }
    
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

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH),fg=fg,is=is,io=io,sex=sex,pd=pd)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), sigma = runif(1, 0, 10),sigma.is = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mu","kjonn","beta","sigma2","sigma2.is","popaar")

# MCMC settings
ni <- 150000
nt <- 8
nb <- 110000
nc <- 3


system.time(mod.3 <- jags(bugs.data, inits, parameters, "cjs-temp-raneff.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE))

# Summarize posteriors
print(mod.3, digits = 3)

save(mod.3,file="mod3rerun.rda")



#final FROH model-------------------------------------------------------------------
# Specify model in BUGS language
sink("cjs-temp-raneff.bug")
cat("
    model {
    
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu + kjonn[sex[i]] + beta*fh[i] + epsilon[t]+ epsilon.is[is[i]]
    p[i,t] <- popaar[is[i],t]
    } #t
    } #i
    for (t in 1:(n.occasions-1)){
    epsilon[t] ~ dnorm(0, tau)
    }
    sigma ~ dunif(0, 10)                     # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)                  # Temporal variance
    
    mu ~ dnorm(0, 0.001)I(-5,5)              # Prior for logit of mean survival
    
    for (h in 1:8){
    epsilon.is[h] ~ dnorm(0, tau.is)
    }
    
    sigma.is ~ dunif(0, 10)                     # Prior for standard deviation
    tau.is <- pow(sigma.is, -2)
    sigma2.is <- pow(sigma.is, 2)                  # Temporal variance
    
    kjonn[1]<-0
    kjonn[2]~ dnorm(0, 0.001)I(-5,5)
    
    beta~ dnorm(0, 0.001)I(-10,10)
    
    for (u in 1:8){
    for (h in 1:(n.occasions-1)){
    popaar[u,h]~ dunif(0, 1)
    }
    }
    
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

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH),fh=fh,is=is,io=io,sex=sex,pd=pd)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), sigma = runif(1, 0, 10),sigma.is = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mu","kjonn","beta","sigma2","sigma2.is","popaar")

# MCMC settings
ni <- 150000
nt <- 8
nb <- 110000
nc <- 3


system.time(mod.4 <- jags(bugs.data, inits, parameters, "cjs-temp-raneff.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE))

# Summarize posteriors
print(mod.4, digits = 3)

save(mod.4,file="mod4rerun.rda")



#island interaction FGRM-------------------------------------------------------------------
# Specify model in BUGS language
sink("cjs-temp-raneff.bug")
cat("
    model {
    
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu + island[is[i]] + beta*fg[i] + beta.is[is[i]]*fg[i]  + epsilon[t]
    p[i,t] <- popaar[is[i],t]
    } #t
    } #i
    for (t in 1:(n.occasions-1)){
    epsilon[t] ~ dnorm(0, tau)
    }
    sigma ~ dunif(0, 10)                     # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)                  # Temporal variance
    
    mu ~ dnorm(0, 0.001)I(-5,5)              # Prior for logit of mean survival

    island[1]<-0
    island[2]~ dnorm(0, 0.001)I(-5,5)
    island[3]~ dnorm(0, 0.001)I(-5,5)
    island[4]~ dnorm(0, 0.001)I(-5,5)
    island[5]~ dnorm(0, 0.001)I(-5,5)
    island[6]~ dnorm(0, 0.001)I(-5,5)
    island[7]~ dnorm(0, 0.001)I(-5,5)
    island[8]~ dnorm(0, 0.001)I(-5,5)
    
    beta~ dnorm(0, 0.001)I(-10,10)

    beta.is[1]<-0
    beta.is[2]~ dnorm(0, 0.001)I(-10,10)
    beta.is[3]~ dnorm(0, 0.001)I(-10,10)
    beta.is[4]~ dnorm(0, 0.001)I(-10,10)
    beta.is[5]~ dnorm(0, 0.001)I(-10,10)
    beta.is[6]~ dnorm(0, 0.001)I(-10,10)
    beta.is[7]~ dnorm(0, 0.001)I(-10,10)
    beta.is[8]~ dnorm(0, 0.001)I(-10,10)

    for (u in 1:8){
    for (h in 1:(n.occasions-1)){
    popaar[u,h]~ dunif(0, 1)
    }
    }
    
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

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH),fg=fg,is=is,io=io,sex=sex,pd=pd)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), sigma = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mu","island","beta","beta.is","sigma2","popaar")

# MCMC settings
ni <- 150000
nt <- 8
nb <- 110000
nc <- 3


system.time(mod.5 <- jags(bugs.data, inits, parameters, "cjs-temp-raneff.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE))

# Summarize posteriors
print(mod.5, digits = 3)

save(mod.5,file="mod5rerun.rda")

#island intercept for FGRM-------------------------------------------------------------------
# Specify model in BUGS language
sink("cjs-temp-raneff.bug")
cat("
    model {
    
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu + island[is[i]] + beta*fg[i]  + epsilon[t]
    p[i,t] <- popaar[is[i],t]
    } #t
    } #i
    for (t in 1:(n.occasions-1)){
    epsilon[t] ~ dnorm(0, tau)
    }
    sigma ~ dunif(0, 10)                     # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)                  # Temporal variance
    
    mu ~ dnorm(0, 0.001)I(-5,5)              # Prior for logit of mean survival
    
    island[1]<-0
    island[2]~ dnorm(0, 0.001)I(-5,5)
    island[3]~ dnorm(0, 0.001)I(-5,5)
    island[4]~ dnorm(0, 0.001)I(-5,5)
    island[5]~ dnorm(0, 0.001)I(-5,5)
    island[6]~ dnorm(0, 0.001)I(-5,5)
    island[7]~ dnorm(0, 0.001)I(-5,5)
    island[8]~ dnorm(0, 0.001)I(-5,5)
    
    beta~ dnorm(0, 0.001)I(-10,10)
    
    for (u in 1:8){
    for (h in 1:(n.occasions-1)){
    popaar[u,h]~ dunif(0, 1)
    }
    }
    
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

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH),fg=fg,is=is,io=io,sex=sex,pd=pd)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), sigma = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mu","island","beta","sigma2","popaar")

# MCMC settings
ni <- 150000
nt <- 8
nb <- 110000
nc <- 3


system.time(mod.5b <- jags(bugs.data, inits, parameters, "cjs-temp-raneff.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE))

# Summarize posteriors
print(mod.5b, digits = 3)

save(mod.5b,file="mod5brerun.rda")



#year ineraction model for FGRM-------------------------------------------------------------------
# Specify model in BUGS language
sink("cjs-temp-raneff.bug")
cat("
    model {
    
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu + year[t] + beta*fg[i] + beta.year[t]*fg[i] + epsilon.is[is[i]]
    p[i,t] <- popaar[is[i],t]
    } #t
    } #i

    mu ~ dnorm(0, 0.001)I(-5,5)              # Prior for logit of mean survival
    
    for (h in 1:8){
    epsilon.is[h] ~ dnorm(0, tau.is)
    }
    
    sigma.is ~ dunif(0, 10)                     # Prior for standard deviation
    tau.is <- pow(sigma.is, -2)
    sigma2.is <- pow(sigma.is, 2)                  # Temporal variance
    
    year[1]<-0
    year[2]~ dnorm(0, 0.001)I(-5,5)
    year[3]~ dnorm(0, 0.001)I(-5,5)
    year[4]~ dnorm(0, 0.001)I(-5,5)
    year[5]~ dnorm(0, 0.001)I(-5,5)
    year[6]~ dnorm(0, 0.001)I(-5,5)
    year[7]~ dnorm(0, 0.001)I(-5,5)
    year[8]~ dnorm(0, 0.001)I(-5,5)
    year[9]~ dnorm(0, 0.001)I(-5,5)
    year[10]~ dnorm(0, 0.001)I(-5,5)
    year[11]~ dnorm(0, 0.001)I(-5,5)
    year[12]~ dnorm(0, 0.001)I(-5,5)
    year[13]~ dnorm(0, 0.001)I(-5,5)
    year[14]~ dnorm(0, 0.001)I(-5,5)
    year[15]~ dnorm(0, 0.001)I(-5,5)
    year[16]~ dnorm(0, 0.001)I(-5,5)

    beta~ dnorm(0, 0.001)I(-10,10)
    beta.year[1]<-0
    beta.year[2]~ dnorm(0, 0.001)I(-10,10)
    beta.year[3]~ dnorm(0, 0.001)I(-10,10)
    beta.year[4]~ dnorm(0, 0.001)I(-10,10)
    beta.year[5]~ dnorm(0, 0.001)I(-10,10)
    beta.year[6]~ dnorm(0, 0.001)I(-10,10)
    beta.year[7]~ dnorm(0, 0.001)I(-10,10)
    beta.year[8]~ dnorm(0, 0.001)I(-10,10)
    beta.year[9]~ dnorm(0, 0.001)I(-10,10)
    beta.year[10]~ dnorm(0, 0.001)I(-10,10)
    beta.year[11]~ dnorm(0, 0.001)I(-10,10)
    beta.year[12]~ dnorm(0, 0.001)I(-10,10)
    beta.year[13]~ dnorm(0, 0.001)I(-10,10)
    beta.year[14]~ dnorm(0, 0.001)I(-10,10)
    beta.year[15]~ dnorm(0, 0.001)I(-10,10)
    beta.year[16]~ dnorm(0, 0.001)I(-10,10)

    
    for (u in 1:8){
    for (h in 1:(n.occasions-1)){
    popaar[u,h]~ dunif(0, 1)
    }
    }
    
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

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH),fg=fg,is=is,io=io,sex=sex,pd=pd)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f),sigma.is = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mu","year","beta","beta.year","sigma2.is","popaar")

# MCMC settings
ni <- 150000
nt <- 8
nb <- 110000
nc <- 3


system.time(mod.6 <- jags(bugs.data, inits, parameters, "cjs-temp-raneff.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE))

# Summarize posteriors
print(mod.6, digits = 3)

save(mod.6,file="mod6rerun.rda")


#year intercept model for FGRM-------------------------------------------------------------------
# Specify model in BUGS language
sink("cjs-temp-raneff.bug")
cat("
    model {
    
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu + year[t] + beta*fg[i] + epsilon.is[is[i]]
    p[i,t] <- popaar[is[i],t]
    } #t
    } #i
    
    mu ~ dnorm(0, 0.001)I(-5,5)              # Prior for logit of mean survival
    
    for (h in 1:8){
    epsilon.is[h] ~ dnorm(0, tau.is)
    }
    
    sigma.is ~ dunif(0, 10)                     # Prior for standard deviation
    tau.is <- pow(sigma.is, -2)
    sigma2.is <- pow(sigma.is, 2)                  # Temporal variance
    
    year[1]<-0
    year[2]~ dnorm(0, 0.001)I(-5,5)
    year[3]~ dnorm(0, 0.001)I(-5,5)
    year[4]~ dnorm(0, 0.001)I(-5,5)
    year[5]~ dnorm(0, 0.001)I(-5,5)
    year[6]~ dnorm(0, 0.001)I(-5,5)
    year[7]~ dnorm(0, 0.001)I(-5,5)
    year[8]~ dnorm(0, 0.001)I(-5,5)
    year[9]~ dnorm(0, 0.001)I(-5,5)
    year[10]~ dnorm(0, 0.001)I(-5,5)
    year[11]~ dnorm(0, 0.001)I(-5,5)
    year[12]~ dnorm(0, 0.001)I(-5,5)
    year[13]~ dnorm(0, 0.001)I(-5,5)
    year[14]~ dnorm(0, 0.001)I(-5,5)
    year[15]~ dnorm(0, 0.001)I(-5,5)
    year[16]~ dnorm(0, 0.001)I(-5,5)
    
    beta~ dnorm(0, 0.001)I(-10,10)

    for (u in 1:8){
    for (h in 1:(n.occasions-1)){
    popaar[u,h]~ dunif(0, 1)
    }
    }
    
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

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH),fg=fg,is=is,io=io,sex=sex,pd=pd)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f),sigma.is = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mu","year","beta","sigma2.is","popaar")

# MCMC settings
ni <- 150000
nt <- 8
nb <- 110000
nc <- 3


system.time(mod.6b <- jags(bugs.data, inits, parameters, "cjs-temp-raneff.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE))

# Summarize posteriors
print(mod.6b, digits = 3)

save(mod.6b,file="mod6brerun.rda")


