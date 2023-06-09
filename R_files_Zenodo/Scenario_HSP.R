################################################################################
# Schaub functions
################################################################################
# Define function to simulate mark-recovery data
simul.mr <- function(S, R, marked){
  n.occasions <- dim(S)[2]
  MR <- matrix(NA, ncol = n.occasions+1, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:n.occasions, marked)
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    MR[i, mark.occ[i]] <- 1    # Write an 1 at the release occasion
    for (t in mark.occ[i]:n.occasions){
      # Bernoulli trial: has individual survived occasion? 
      sur <- rbinom(1, 1, S[i,t])
      if (sur==1) next    # If still alive, move to next occasion 
      # Bernoulli trial: has dead individual been recovered? 
      rp <- rbinom(1, 1, R[i,t])
      if (rp==0){
        MR[i,t+1] <- 0
        break
      }
      if (rp==1){
        MR[i,t+1] <- 1
        break
      }
    } #t
  } #i
  # Replace the NA in the file by 0
  MR[which(is.na(MR))] <- 0
  return(MR)
}

marray.dead <- function(MR){
  nind <- dim(MR)[1]
  n.occasions <- dim(MR)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Create vector with occasion of marking 
  get.first <- function(x) min(which(x!=0))
  f <- apply(MR, 1, get.first)
  # Calculate the number of released individuals at each time period
  first <- as.numeric(table(f))
  for (t in 1:n.occasions){
    m.array[t,1] <- first[t]
  }
  # Fill m-array with recovered individuals
  rec.ind <- which(apply(MR, 1, sum)==2)
  rec <- numeric()
  for (i in 1:length(rec.ind)){
    d <- which(MR[rec.ind[i],(f[rec.ind[i]]+1):n.occasions]==1)
    rec[i] <- d + f[rec.ind[i]]
    m.array[f[rec.ind[i]],rec[i]] <- m.array[f[rec.ind[i]],rec[i]] + 1
  }
  # Calculate the number of individuals that are never recovered
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1]-sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

################################################################################
################################################################################
################################################################################
# MODELS
################################################################################
################################################################################
################################################################################


################################################################################
# state-space model 
################################################################################
# Specify model in BUGS language
sink("ssm.jags")
cat("
    model { 
    # Priors and constraints

    Na[1] ~ dunif(0, 3000)            # Prior for initial population size
    mean.lambda ~ dunif(0, 10)          # Prior for mean growth rate
    sigma.proc ~ dunif(0, 10)           # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)

    tau.proc <- pow(sigma.proc, -2)
    sigma.obs ~ dunif(0, 500)           # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(n.occasions-1)){
    lambda[t] ~ dnorm(mean.lambda, tau.proc) T(0,)
    Na[t+1] <- Na[t] * lambda[t] 
    }
    # Observation process
    for (t in 1:n.occasions) {
    y[t,1] ~ dnorm(Na[t], tau.obs)
    y[t,2] ~ dnorm(Na[t], tau.obs)
    }
    }
    ",fill = TRUE)
sink()



################################################################################
# capture-mark-recapture model
################################################################################
# Specify model in BUGS language
sink("cmr.jags")
cat("
    model {
    
    # Priors and constraints
    for (t in 1:n.occasions){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    rj[t] <- mean.r
    ra[t] <- mean.r
    }
    mean.sj ~ dunif(0, 1)              # Prior for mean juv. survival
    mean.sa ~ dunif(0, 1)              # Prior for mean ad. survival
    mean.r ~ dunif(0, 1)              # Prior for mean juv. recovery

    # Define the multinomial likelihoods
    # Calculate the number of birds released each year
    for (t in 1:n.occasions){
    marr.j[t,1:(n.occasions+1)] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:(n.occasions+1)] ~ dmulti(pr.a[t,], rel.a[t])
    }
    # Define the cell probabilities of the juvenile m-array
    # Main diagonal
    for (t in 1:n.occasions){
    pr.j[t,t] <- (1-sj[t])*rj[t]
    # Further above main diagonal
    for (j in (t+2):n.occasions){
    pr.j[t,j] <- sj[t]*prod(sa[(t+1):(j-1)])*(1-sa[j])*ra[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    } #j
    } #t
    for (t in 1:(n.occasions-1)){
    # One above main diagonal
    pr.j[t,t+1] <- sj[t]*(1-sa[t+1])*ra[t+1] 
    } #t
    # Last column: probability of non-recovery
    for (t in 1:n.occasions){
    pr.j[t,n.occasions+1] <- 1-sum(pr.j[t,1:n.occasions])
    } #t
    # Define the cell probabilities of the adult m-array
    # Main diagonal
    for (t in 1:n.occasions){
    pr.a[t,t] <- (1-sa[t])*ra[t]
    # Above main diagonal
    for (j in (t+1):n.occasions){
    pr.a[t,j] <- prod(sa[t:(j-1)])*(1-sa[j])*ra[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recovery
    for (t in 1:n.occasions){
    pr.a[t,n.occasions+1] <- 1-sum(pr.a[t,1:n.occasions])
    } #t
    }
    ",fill = TRUE)
sink()


################################################################################
################################################################################
# integrated population model...
################################################################################
################################################################################
# Specify model in BUGS language
sink("ipm.jags")
cat("
    model {

#-------------------------------------------------
# 1. Define the priors for the parameters
#-------------------------------------------------
# Observation error
tauy <- pow(sigma.y, -2)
sigma.y ~ dunif(0, 500)
sigma2.y <- pow(sigma.y, 2)

# Initial population sizes
n1 ~ dunif(300,800)   # 1-year
nad ~ dunif(700,1300)   # Adults
N1[1] <- round(n1)
Nad[1] <- round(nad)

########################################################
# fecundity prior
########################################################
mean.fec ~ dunif(0,10)

for (t in 1:(n.occasions-1)){
  f[t] <- mean.fec
}


########################################################
# 2. Derived parameters
########################################################
# Population growth rate
# for (t in 1:(n.occasions-1)){
#    lambda[t] <- Ntot[t+1] / Ntot[t]
# }

########################################################
# 3. The likelihoods of the single data sets
########################################################
# 3.1. Likelihood for population population count data (state-space model)
   # 3.1.1 System process
   for (t in 2:n.occasions){
      mean1[t] <- f[t-1] * 0.5 * sj[t-1] * Nad[t]
      N1[t] ~ dpois(mean1[t])
      Sad[t-1] ~ dbin(sa[t-1], Nad[t-1])
      Nad[t] <- Sad[t-1] + N1[t-1]
      }


   # 3.1.2 Observation process
   for (t in 1:n.occasions){
      y[t,1] ~ dnorm(Nad[t], tauy)
      y[t,2] ~ dnorm(Nad[t], tauy)
      }

    #############################################################
    # cmr model
    #############################################################

    # Priors and constraints
    for (t in 1:n.occasions){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    rj[t] <- mean.r
    ra[t] <- mean.r
    }

    mean.sj ~ dunif(0, 1)              # Prior for mean juv. survival
    mean.sa ~ dunif(0, 1)              # Prior for mean ad. survival
    mean.r ~ dunif(0, 1)              # Prior for mean juv. recovery
    # mean.ra ~ dunif(0, 1)              # Prior for mean ad. recovery
    
    # Define the multinomial likelihoods
    # Calculate the number of birds released each year
    for (t in 1:n.occasions){
    marr.j[t,1:(n.occasions+1)] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:(n.occasions+1)] ~ dmulti(pr.a[t,], rel.a[t])
    }
    # Define the cell probabilities of the juvenile m-array
    # Main diagonal
    for (t in 1:n.occasions){
    pr.j[t,t] <- (1-sj[t])*rj[t]
    # Further above main diagonal
    for (j in (t+2):n.occasions){
    pr.j[t,j] <- sj[t]*prod(sa[(t+1):(j-1)])*(1-sa[j])*ra[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    } #j
    } #t
    for (t in 1:(n.occasions-1)){
    # One above main diagonal
    pr.j[t,t+1] <- sj[t]*(1-sa[t+1])*ra[t+1] 
    } #t
    # Last column: probability of non-recovery
    for (t in 1:n.occasions){
    pr.j[t,n.occasions+1] <- 1-sum(pr.j[t,1:n.occasions])
    } #t
    # Define the cell probabilities of the adult m-array
    # Main diagonal
    for (t in 1:n.occasions){
    pr.a[t,t] <- (1-sa[t])*ra[t]
    # Above main diagonal
    for (j in (t+1):n.occasions){
    pr.a[t,j] <- prod(sa[t:(j-1)])*(1-sa[j])*ra[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recovery
    for (t in 1:n.occasions){
    pr.a[t,n.occasions+1] <- 1-sum(pr.a[t,1:n.occasions])
    } #t
    }
    ",fill = TRUE)
sink()









################################################################################
################################################################################
################################################################################
# SIMULATION
################################################################################
################################################################################
################################################################################
sims <- 100
npar <- 6
tru <- matrix(0, sims, npar)
est.cmr <- array(0, dim = c(sims, npar, 3))
est.ipm <- array(0, dim = c(sims, npar, 3))
est.ssm <- array(0, dim = c(sims, npar, 3))



# loop
for (ii in 1:sims){
  
  
# 8.3.2. Age-dependent parameters
n.occasions <- 20                  # Number of occasions

marked.j <- rep(200, n.occasions)   # Annual number of newly marked young
marked.a <- rep(200, n.occasions)    # Annual number of newly marked adults
# sjuv <- -0.7                         # Juvenile survival probability
sad <- 1.4                          # Adult survival probability
rjuv <- 0.2                        # Juvenile recovery probability
rad <- 0.2                         # Adult recovery probability
rj <- c(rjuv, rep(rad, n.occasions-1))
var <- 0.5

# Define matrices with survival and recovery probabilities
SJ <- matrix(0, ncol = n.occasions, nrow = sum(marked.j))
sjuv <- rnorm((n.occasions * marked.j[1]), -0.7, var)
release <- rep(seq(1,n.occasions,1), each = 200)

for (i in 1:nrow(SJ)){
  SJ[i,release[i]] <- plogis(sjuv[i])
  if (release[i] < n.occasions){
    SJ[i,(release[i]+1):n.occasions] <- plogis(sjuv[i] + 2.1)
  }
}

SA <- matrix(plogis(rnorm((marked.a * n.occasions), sad, var)), ncol = n.occasions, nrow = sum(marked.a))
RJ <- matrix(0, ncol = n.occasions, nrow = sum(marked.j))
for (i in 1:length(marked.j)){
  RJ[(sum(marked.j[1:i])-marked.j[i]+1):sum(marked.j[1:i]),i:n.occasions] <- matrix(rep(rj[1:(n.occasions-i+1)],marked.j[i]), ncol = n.occasions-i+1, byrow = TRUE)
}
RA <- matrix(rad, ncol = n.occasions, nrow = sum(marked.a))

# Execute simulation function
MRj <- simul.mr(SJ, RJ, marked.j)
MRa <- simul.mr(SA, RA, marked.a)

# Summarize data in m-arrays
marr.j <- marray.dead(MRj)
marr.a <- marray.dead(MRa)

################################################################################
# create population
################################################################################
sa <- plogis(sad)
sj <- mean(plogis(sjuv))
fec <- 1.1

Nj <- NULL
Na <- NULL

Nj[1] <- 550
Na[1] <- 1000

y <- matrix(NA, n.occasions, 2)
var <- 100
y[1,] <- rnorm(2, Na[1], var)

Sj <- NULL
Sa <- NULL
plogis(-0.7)
for (t in 2:n.occasions){
  Sj[t-1] <- rbinom(1, Nj[t-1], sj)
  Sa[t-1] <- rbinom(1, Na[t-1], sa)
  Na[t] <- Sj[t-1] + Sa[t-1]
  y[t,] <- rnorm(2, Na[t], var)

  Nj[t] <- rpois(1, Na[t] * fec * 0.5)  
}
# plot(Sj)
# plot(Na)

# warnings()



######################################
# MODELS.....
######################################
require(jagsUI)


######################################
# Chain Info
######################################
ni <- 50000
nt <- 5
nb <- 25000
nc <- 3


######################################
# SSM
######################################
jags.data <- list(y = y, n.occasions = n.occasions)
inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), Na = c(runif(1, 900, 1100), rep(NA, (n.occasions-1))))} 
parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "Na",'sigma.obs')
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


######################################
# CMR
######################################
jags.data <- list(marr.j = marr.j, marr.a = marr.a, n.occasions = dim(marr.j)[2]-1, rel.j = rowSums(marr.j), rel.a = rowSums(marr.a))
inits <- function(){list(mean.sj = runif(1, 0, 1), mean.sa = runif(1, 0, 1), mean.rj = runif(1, 0, 1), mean.ra = runif(1, 0, 1))}  
parameters <- c("mean.sj", "mean.sa", "mean.r")
mr <- jags(jags.data, inits, parameters, "cmr.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
# print(mr)


######################################
# IPM
######################################
jags.data <- list(marr.j = marr.j, marr.a = marr.a, n.occasions = dim(marr.j)[2]-1, rel.j = rowSums(marr.j), rel.a = rowSums(marr.a), y = y)
inits <- function(){list(mean.sj = runif(1, 0, 1), mean.sa = runif(1, 0, 1), mean.r = runif(1, 0, 1), mean.fec = runif(1, 0, 10), n1 = 550, nad = 1000, sigma.y = runif(1, 50, 150))}  
parameters <- c("mean.sj", "mean.sa", "mean.r", "mean.fec", "N1", "Nad", "Ntot", "sigma2.y",'sigma.y', "lambda",'f')
ipm <- jags(jags.data, inits, parameters, "ipm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# print(ipm)

# plot(ipm$sims.list$mean.sa)


est.cmr[ii,1,1] <- mr$q50$mean.sa
est.cmr[ii,2,1] <- mr$q50$mean.sj
est.cmr[ii,3,1] <- mr$q50$mean.r
est.cmr[ii,4,1] <- NA
est.cmr[ii,5,1] <- NA
est.cmr[ii,6,1] <- NA

est.cmr[ii,1,2] <- mr$q2.5$mean.sa
est.cmr[ii,2,2] <- mr$q2.5$mean.sj
est.cmr[ii,3,2] <- mr$q2.5$mean.r
est.cmr[ii,4,2] <- NA
est.cmr[ii,5,2] <- NA
est.cmr[ii,6,2] <- NA

est.cmr[ii,1,3] <- mr$q97.5$mean.sa
est.cmr[ii,2,3] <- mr$q97.5$mean.sj
est.cmr[ii,3,3] <- mr$q97.5$mean.r
est.cmr[ii,4,3] <- NA
est.cmr[ii,5,3] <- NA
est.cmr[ii,6,3] <- NA

est.ipm[ii,1,1] <- ipm$q50$mean.sa
est.ipm[ii,2,1] <- ipm$q50$mean.sj
est.ipm[ii,3,1] <- ipm$q50$mean.r
est.ipm[ii,4,1] <- ipm$Rhat$mean.fec
est.ipm[ii,5,1] <- ipm$q50$mean.fec
est.ipm[ii,6,1] <- ipm$q50$sigma.y

est.ipm[ii,1,2] <- ipm$q2.5$mean.sa
est.ipm[ii,2,2] <- ipm$q2.5$mean.sj
est.ipm[ii,3,2] <- ipm$q2.5$mean.r
est.ipm[ii,4,2] <- ipm$Rhat$mean.fec
est.ipm[ii,5,2] <- ipm$q2.5$mean.fec
est.ipm[ii,6,2] <- ipm$q2.5$sigma.y

est.ipm[ii,1,3] <- ipm$q97.5$mean.sa
est.ipm[ii,2,3] <- ipm$q97.5$mean.sj
est.ipm[ii,3,3] <- ipm$q97.5$mean.r
est.ipm[ii,4,3] <- ipm$Rhat$mean.fec
est.ipm[ii,5,3] <- ipm$q97.5$mean.fec
est.ipm[ii,6,3] <- ipm$q97.5$sigma.y

est.ssm[ii,1,1] <- NA
est.ssm[ii,2,1] <- NA
est.ssm[ii,3,1] <- NA
est.ssm[ii,4,1] <- NA
est.ssm[ii,5,1] <- NA
est.ssm[ii,6,1] <- ssm$q50$sigma.obs

est.ssm[ii,1,2] <- NA
est.ssm[ii,2,2] <- NA
est.ssm[ii,3,2] <- NA
est.ssm[ii,4,2] <- NA
est.ssm[ii,5,2] <- NA
est.ssm[ii,6,2] <- ssm$q2.5$sigma.obs

est.ssm[ii,1,3] <- NA
est.ssm[ii,2,3] <- NA
est.ssm[ii,3,3] <- NA
est.ssm[ii,4,3] <- NA
est.ssm[ii,5,3] <- NA
est.ssm[ii,6,3] <- ssm$q97.5$sigma.obs

tru[ii,1] <- sa
tru[ii,2] <- mean(plogis(sjuv))
tru[ii,3] <- 0.2
tru[ii,4] <- 0.2
tru[ii,5] <- fec
tru[ii,6] <- 100

print(ii)

}
save.image("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\Resubmission\\\\HiSP_final_results3.RData")

################################################################################
################################################################################
################################################################################
# PLOTS
################################################################################
################################################################################
################################################################################
load("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\Resubmission\\\\HiSP_final_results2.RData")
pdf("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\new_TEX\\\\Figure8.pdf", width=10, height=10)
par(family = 'sans', mar = c(3.1,5.1,2.1,2.1), mfrow = c(3,2))
boxplot(est.cmr[1:18,1], est.ipm[1:18,1], ylim = c(0.725,0.825), ylab = expression(phi['ad']), 
        col = c('royalblue1', 'yellow'), cex.lab = 2, names = c('Seber', expression(IPM['HiSP'])))
abline(plogis(1.4),0, lty = 2, lwd = 2)
boxplot(est.cmr[1:18,2], est.ipm[1:18,2], ylim = c(0.25,0.375), ylab = expression(phi['juv']), 
        col = c('royalblue1', 'yellow'), cex.lab = 2, names = c('Seber',expression(IPM['HiSP'])))
abline(plogis(-0.7),0, lty = 2, lwd = 2)
boxplot(est.cmr[1:18,3], est.ipm[1:18,3], ylim = c(0.17,0.21), ylab = expression(italic(r)), 
        col = c('royalblue1', 'yellow'), cex.lab = 2, names = c('Seber',expression(IPM['HiSP'])))
abline(0.2,0, lty = 2, lwd = 2)

boxplot(est.ssm[1:18,6], est.ipm[1:18,6], ylim = c(70,140), ylab = expression(sigma['y']),
        col = c('royalblue1', 'yellow'), names = c('SSM',expression(IPM['HiSP'])), cex.lab = 2)
abline(100,0, lty = 2, lwd = 2)
boxplot(I(est.ipm[1:18,5]/2), ylim = c(0.5,0.9), col = 'royalblue1', ylab = expression(italic(f)), 
        names = 'IPM3', cex.lab = 2)
abline(0.55, 0, lty = 2, lwd = 2)
dev.off()



# save.image("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\new_results.RData")


################################################################################
# calculate relative bias
################################################################################

############################################
# CMR and SSM models
############################################
# adult survival
mean((est.cmr[1:18,1] - plogis(1.4))/plogis(1.4)) # relative
mean(est.cmr[1:18,1] - plogis(1.4))               # absolute

# juvenile survival
mean((est.cmr[1:18,2] - plogis(-0.7))/plogis(-0.7)) # relative
mean(est.cmr[1:18,2] - plogis(-0.7))                # absolute

# recovery
mean((est.cmr[1:18,3] - 0.2)/0.2) # relative
mean(est.cmr[1:18,3] - 0.2)       # absolute

# observation error
mean((est.ssm[1:18,6] - 100)/100) # relative
mean(est.ssm[1:18,6] - 100)       # absolute





############################################
# CMR and SSM models
############################################
# adult survival
mean((est.ipm[1:18,1] - plogis(1.4))/plogis(1.4)) # relative
mean(est.ipm[1:18,1] - plogis(1.4))               # absolute

# juvenile survival
mean((est.ipm[1:18,2] - plogis(-0.7))/plogis(-0.7)) # relative
mean(est.ipm[1:18,2] - plogis(-0.7))                # absolute

# recovery
mean((est.ipm[1:18,3] - 0.2)/0.2) # relative
mean(est.ipm[1:18,3] - 0.2)       # absolute

# observation error
mean((est.ipm[1:18,6] - 100)/100) # relative
mean(est.ipm[1:18,6] - 100)       # absolute

# fecundity
mean((est.ipm[1:18,5] - 1.1)/1.1) # relative
mean(est.ipm[1:18,5] - 1.1)       # absolute

boxplot(est.ipm[1:18,5], ylim = c(0,3), col = 'royalblue1', ylab = 'Fecundity', 
        names = 'IPM3', cex.lab = 2)
abline(1.1, 0, lty = 2, lwd = 2)












