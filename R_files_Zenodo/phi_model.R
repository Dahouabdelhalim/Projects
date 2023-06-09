data <- read.csv("/Users/cacevedo/Dropbox/Projects/Brant_EggSize/survival_input.csv")
summary(data)



################################################################
# read mean clutch
################################################################
mu <- read.csv("/Users/cacevedo/Dropbox/Projects/Brant_EggSize/mean_clutch_vols.csv", header = T)
mu <- mu[,2:4]
summary(mu)
mu$year <- mu$year - 1986
require(plyr)
data <- join(data, mu, type = 'left', by = c('year','nest'))


################################################################
# z-standardize volume by position
################################################################
data$meanVOL <- scale(data$muVOL)
summary(data$meanVOL)

data$resVOL <- data$VOL - data$muVOL
data$rVOL <- scale(data$resVOL)
summary(data$resVOL)
summary(data$rVOL)
summary(data$meanVOL)
# change all NA's to 0's for rVol and meanVOL
data$rVOL[is.na<https://nam04.safelinks.protection.outlook.com/?url=http%3A%2F%2Fis.na&data=01%7C01%7Ctriecke%40cabnr.unr.edu%7Cc768503b856941976d9c08d76ed62e4d%7C523b4bfc0ebd4c03b2b96f6a17fd31d8%7C1&sdata=YFI59EFXEd28SX4BJUmBKb%2B%2FdNVmGcGTpzSnmllFr%2Bc%3D&reserved=0>(data$rVOL)] <- 0
data$meanVOL[is.na<https://nam04.safelinks.protection.outlook.com/?url=http%3A%2F%2Fis.na&data=01%7C01%7Ctriecke%40cabnr.unr.edu%7Cc768503b856941976d9c08d76ed62e4d%7C523b4bfc0ebd4c03b2b96f6a17fd31d8%7C1&sdata=YFI59EFXEd28SX4BJUmBKb%2B%2FdNVmGcGTpzSnmllFr%2Bc%3D&reserved=0>(data$meanVOL)] <- 0

###########################################################



summary(vol)


names(data)
data$brood <- paste(data$nest, data$year, sep = '') # this combines the two columns (nest and year) together
summary(data$brood)
data$brood <- as.factor(data$brood)
data$brood <- as.numeric(data$brood)

###
### sort, but don't mix up what goes to each brood through sorting
###
data <- data[order(data$broodID),]





R <- data$recap; summary(as.factor(R))
pos <- data$pos; summary(pos)
summary(as.factor(pos))
vol <- data$vol; summary(vol)
year <- data$year; summary(as.factor(year))
S <- as.numeric(data$state); summary(S)
W <- data$tagW; summary(W)
Ai <- data$tagA; summary(Ai)
resVOL <- data$resVol
rVOL <- data$rVOL
meanVOL <- data$meanVOL
# gosling age at capture
age <- data$age

names(data)
A <- data.frame(data$tagA, data$year)
A <- unique(A)
?order()
tmp <- A[order(A$data.year),]
names(A)

A <- tmp
A <- A[,1]

# gosling gender
# sex <- rep(0, nrow(data))
# for (i in 1:length(sex)){
#   if(data$sex[i] == 'M'){sex[i] <- 1}
# }



# find out how to sort by a column, sort among by year...

# this will create the brood effect

n.brood <- max(data$brood)
brood <- data$brood

# install.packages('jagsUI')
library(jagsUI)
getwd()









loop.start <- NULL
loop.stop <- NULL

for (j in 1:max(brood)){
  loop.start[j] <- min(which(brood == j))
  loop.stop[j] <- max(which(brood == j))
}
loop.start
loop.stop

pi.data <- rep(NA,length(loop.start))
summary(loop.start)
summary(loop.stop)
for (j in 1:n.brood){
  if(any(R[loop.start[j]]:R[loop.stop[j]] == 1)){
    pi.data[j] <- 1
  }
}
summary(pi.data)

z.start <- R
z.start[z.start == 0] <- NA
summary(z.start)



















################################################################
# model
################################################################
require(rjags)
sink("phi_model.jags")
cat("
    model {

    ###########################################################
    # annual temporal component # what drives the temporal variation
    ###########################################################
    alpha.mu<https://nam04.safelinks.protection.outlook.com/?url=http%3A%2F%2Falpha.mu&data=01%7C01%7Ctriecke%40cabnr.unr.edu%7Cc768503b856941976d9c08d76ed62e4d%7C523b4bfc0ebd4c03b2b96f6a17fd31d8%7C1&sdata=%2BGTi%2BNOF3G%2FsU7ITtPNlvf69ZYVYAKV%2FR1PyANE3p%2Bc%3D&reserved=0> ~ dnorm(0, 0.001)
    beta.mu.trend ~ dnorm(0, 0.001)
    beta.mu.phenology ~ dnorm(0, 0.001)
    tau.phi.t <- pow(sigma.phi.t, -2)
    sigma.phi.t ~ dunif(0,10)

    # this models the mean survival based on trend and phenology from 1987 to 2007
    for (t in 1:n.years){
      mu.phi[t] <- alpha.mu<https://nam04.safelinks.protection.outlook.com/?url=http%3A%2F%2Falpha.mu&data=01%7C01%7Ctriecke%40cabnr.unr.edu%7Cc768503b856941976d9c08d76ed62e4d%7C523b4bfc0ebd4c03b2b96f6a17fd31d8%7C1&sdata=%2BGTi%2BNOF3G%2FsU7ITtPNlvf69ZYVYAKV%2FR1PyANE3p%2Bc%3D&reserved=0> + beta.mu.trend * t + beta.mu.phenology * A[t]
      eps.phi[t] ~ dnorm(mu.phi[t], tau.phi.t)
      p[t] ~ dnorm(mu.p, tau.p)
      logit(pG[t]) <- p[t]
    }

    mu.p ~ dnorm(-1,0.1)
    sigma.p ~ dunif(0,5)
    tau.p <- pow(sigma.p, -2)


    ###########################################################
    # model what drives the individual variation
    ###########################################################
    beta.within ~ dnorm(0, 0.001)
    beta.vol ~ dnorm(0, 0.001)
    beta.meanVOL ~ dnorm(0, 0.0001)
    beta.rVOL ~  dnorm(0, 0.001)
    beta.int<https://nam04.safelinks.protection.outlook.com/?url=http%3A%2F%2Fbeta.int&data=01%7C01%7Ctriecke%40cabnr.unr.edu%7Cc768503b856941976d9c08d76ed62e4d%7C523b4bfc0ebd4c03b2b96f6a17fd31d8%7C1&sdata=cgvOv8Ti9RKB8ibE5tf751jmixqzItnLBVFQDKEin%2FQ%3D&reserved=0> ~ dnorm(0, 0.001)

    for (j in 1:14){
      lay[j] ~ dnorm(0, 0.001)
    }
    lay[15] <- 0


    ###
    ###
    ###
    for (k in 1:n.brood){
      pi[k] ~ dbern(xi[k])
      xi[k] <- ifelse(s[k] > 0, pG[year[k]], 0)
      s[k] <- sum(z[loop.start[k]:loop.stop[k]])
    }

    for (i in 1:n){
      R[i] ~ dbern(pi[brood[i]] * z[i])
      z[i] ~ dbern(mu[i])
    }

    # this add the brood effect onto survival based on laying order position within clutches and egg volume

    # recapture of each individual is a bernoulli trial with probability mu
    # individual, brood, and temporal effects on phi


    for (i in 1:n){
     logit(mu[i]) <- beta.within * W[i] + beta.meanVOL * meanVOL[i] + eps.phi[year[i]] +
                  lay[pos[i]] + beta.rVOL * rVOL[i] + beta.int<https://nam04.safelinks.protection.outlook.com/?url=http%3A%2F%2Fbeta.int&data=01%7C01%7Ctriecke%40cabnr.unr.edu%7Cc768503b856941976d9c08d76ed62e4d%7C523b4bfc0ebd4c03b2b96f6a17fd31d8%7C1&sdata=cgvOv8Ti9RKB8ibE5tf751jmixqzItnLBVFQDKEin%2FQ%3D&reserved=0> * meanVOL[i] * rVOL[i]
    }


###
### everything below here is a derived parameter
###

    for (j in 1:14){
      pred.lay[j] <- eps.phi[11] + lay[j]
    }

    for (k in 1:100){
      pred.vol[k] <- eps.phi[11] + beta.meanVOL * egg.vol[k]
    }


    for (t in 1:n.years){
      pred.t[t] <- alpha.mu<https://nam04.safelinks.protection.outlook.com/?url=http%3A%2F%2Falpha.mu&data=01%7C01%7Ctriecke%40cabnr.unr.edu%7Cc768503b856941976d9c08d76ed62e4d%7C523b4bfc0ebd4c03b2b96f6a17fd31d8%7C1&sdata=%2BGTi%2BNOF3G%2FsU7ITtPNlvf69ZYVYAKV%2FR1PyANE3p%2Bc%3D&reserved=0> + beta.mu.trend * t

      pred.annual[t] <- eps.phi[t]
    }


    for (k in 1:100){
      pred.w[k] <- eps.phi[11] + beta.within * hatch.date[k]
    }


    }

    ",fill = TRUE)
sink()

################################################################
# export and plot results
################################################################
#Bundle data
# dat <- list(m.sex = m.sex, m.vol = m.vol, m.W = m.W, m.Ai = m.Ai, mass = mass,
#             R = R, pos = pos, vol = vol, year = year, S = S, W = W, n.years = max(year), n = length(R), A = A, Ai = Ai, age = age, sex = sex)
egg.vol <- seq(min(meanVOL), max(meanVOL), length.out = 100)
hatch.date <- seq(min(W), max(W), length.out = 100)
dat <- list(R = R, A = A, n = length(vol), pos = pos,
            meanVOL = meanVOL[,1], rVOL = rVOL[,1], age = age, W = W, n.years = 21, year = year, S = S,
            brood = brood, n.brood = n.brood, egg.vol = egg.vol, hatch.date = hatch.date,
            z = z.start, pi = pi.data,loop.start=loop.start, loop.stop =loop.stop)
#mean(data$VOL, na.rm= T)
#hist(data$muVOL)


# Initial values
inits <- function() list(mean.p = 0)

# Parameters to monitor
pars <- c('alpha.mu<https://nam04.safelinks.protection.outlook.com/?url=http%3A%2F%2Falpha.mu&data=01%7C01%7Ctriecke%40cabnr.unr.edu%7Cc768503b856941976d9c08d76ed62e4d%7C523b4bfc0ebd4c03b2b96f6a17fd31d8%7C1&sdata=%2BGTi%2BNOF3G%2FsU7ITtPNlvf69ZYVYAKV%2FR1PyANE3p%2Bc%3D&reserved=0>','beta.mu.trend','beta.mu.phenology', 'eps.phi',
          'beta.vol', 'beta.within', 'lay', 'pred.lay', 'pG',
          'sigma.brood', 'beta.meanVOL', 'beta.rVOL', 'beta.int<https://nam04.safelinks.protection.outlook.com/?url=http%3A%2F%2Fbeta.int&data=01%7C01%7Ctriecke%40cabnr.unr.edu%7Cc768503b856941976d9c08d76ed62e4d%7C523b4bfc0ebd4c03b2b96f6a17fd31d8%7C1&sdata=cgvOv8Ti9RKB8ibE5tf751jmixqzItnLBVFQDKEin%2FQ%3D&reserved=0>', 'pred.w',
          'pred.annual', 'pred.t', 'pred.vol','b')


######################################################################################
# MCMC settings
######################################################################################
n.chains <- 2
n.thin <- 5
n.adapt <- 500
n.iter <- 15000
n.burnin <- 5000


######################################################################################
# compile model!
######################################################################################
m.phi <- jags(dat, inits, pars, "phi_model.jags",
              n.chains = n.chains, n.adapt = n.adapt, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
summary(m.phi)
par(mfrow = c(1,1))
hist(m.phi$mean$b)

m.phi$mean$beta.within
m.phi$q2.5$beta.within
m.phi$q97.5$beta.within
