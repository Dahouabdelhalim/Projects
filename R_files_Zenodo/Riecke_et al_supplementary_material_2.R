######################################################################
# 
# Supplementary Material 2: RD with Markovian emigration and pooled detection among individuals
# 
# In this script we simulate (lines 1-200), analayze (lines 200-325),
# and plot (lines 325-end) robust design capture, mark, recapture datasets
#
######################################################################
# true survival and transition probabilities
######################################################################
# set.seed(124)
params <- 4
sims <- 100

truth <- matrix(NA, sims, params)
estimate <- matrix(NA, sims, params)
estimateL <- matrix(NA, sims, params)
estimateH <- matrix(NA, sims, params)


for (k in 1:sims){

n.years <- 10      # number of primaries
n.sec.occ <- 4     # number of secondaries
n <- 100            # releases


######################################################################
# true survival and transition probabilities
######################################################################
s <- rep(runif(1,0.6,0.9), n.years)
gP <- rep(runif(1,0.2,0.6), n.years)
gPP <- rep(runif(1,0.1,0.4), n.years)

######################################################################
# true detection probabilities
######################################################################
p <- matrix(NA, n.years, n.sec.occ+1)
p[,1:n.sec.occ] <- runif(n.years*n.sec.occ,0.05,0.35)

for (i in 1:n.years){
  p[i,n.sec.occ+1] <- 1 - prod(1 - p[i,1:n.sec.occ])
}

######################################################################
# create true latent state matrix (release is synonymous with true first
# year of breeding)
######################################################################
true <- matrix(NA, n * n.years, n.years)
release <- rep(1:n.years, each = n)

######################################################################
# create 4-d array of transition probabilities...
######################################################################
PSI.STATE <- array(NA, dim=c(3, 3, n * n.years, n.years-1))
for (i in 1:(n*n.years)){
  for (t in 1:(n.years-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      s[t]*(1-gPP[t]),   s[t]*gPP[t],    1-s[t],
      s[t]*(1-gP[t]),  s[t]*gP[t],   1-s[t],
      0,                0,             1       ), nrow = 3, byrow = TRUE)
  }
}

######################################################################
# first year of breeding
######################################################################
for(i in 1:(n*n.years)){
  true[i,release[i]] <- 1
}

######################################################################
# fill true latent state matrix
######################################################################
for (i in 1:(n*n.years)){
  for (t in 1:(n.years)){
    if(t > release[i]){
      if(true[i,t-1] == 1){
        true[i,t] <- which(rmultinom(1,1,c(PSI.STATE[1,1,i,t-1],PSI.STATE[1,2,i,t-1],PSI.STATE[1,3,i,t-1]))==1)
      }
      if(true[i,t-1] == 2){
        true[i,t] <- which(rmultinom(1,1,c(PSI.STATE[2,1,i,t-1],PSI.STATE[2,2,i,t-1],PSI.STATE[2,3,i,t-1]))==1)
      }
      if(true[i,t-1] == 3){
        true[i,t] <- which(rmultinom(1,1,c(PSI.STATE[3,1,i,t-1],PSI.STATE[3,2,i,t-1],PSI.STATE[3,3,i,t-1]))==1)
      }
    }
  }
} 


######################################################################
# 'master' observations
######################################################################
obs <- array(NA, c(n*n.years, n.years, n.sec.occ))
for (i in 1:(n.years*n)){
  for (t in release[i]:n.years){
    if (true[i,t] == 1){
      for (j in 1:n.sec.occ){
        obs[i,t,j] <- rbinom(1,1,p[t,j])
      }
    }
  }
}
obs[is.na(obs)] <- 0


######################################################################
# format as Bayesian MS-RD input
######################################################################
ch <- matrix(NA, n*n.years, n.years)
for (i in 1:(n.years*n)){
  for (t in 1:n.years){
    ifelse(any(obs[i,t,1:n.sec.occ]) == 1, ch[i,t] <- 1, ch[i,t] <- 2)
  }
}


###########################################################################
# summarize detection by primary and secondary occasion
###########################################################################
test <- matrix(NA, n*n.years, n.years)
for (i in 1:nrow(test)){
  for (j in 1:ncol(test)){
    test[i,j] <- sum(obs[i,j,])
  }
}

seen <- array(0, c(n*n.years, n.years, n.sec.occ))
missed <- array(0, c(n*n.years, n.years, n.sec.occ))

for (i in 1:nrow(test)){
  for (t in 1:ncol(test)){
    for (j in 1:n.sec.occ){
      if(test[i,t] > 1 & obs[i,t,j] == 1){seen[i,t,j] <- 1}
      if(test[i,t] >= 1 & obs[i,t,j] == 0){missed[i,t,j] <- 1}
    }
  }
}

obs[1:20,,1]
test[1:20,]
seen[1:20,,1]
missed[1:20,,1]

yes <- matrix(NA, n.years, n.sec.occ)
no <- matrix(NA, n.years, n.sec.occ)

for (i in 1:nrow(yes)){
  for (j in 1:ncol(yes)){
    yes[i,j] <- sum(seen[,i,j])
    no[i,j] <- sum(missed[,i,j])
  }
}

total <- yes + no


################################################################
# cut individuals never released
################################################################
get.first <- function(x)min(which (x != 2))
first <- apply(ch,1,get.first); first[first == "Inf"] <- NA
ch <- subset(ch, !is.na(first)) 
true <- subset(true, !is.na(first))
obs <- obs[-c(is.na(first)),,]
first <- subset(first, !is.na(first))



################################################################
# cut individuals released in last primary occasion
################################################################
ch <- subset(ch, first != n.years)
first <- subset(first, first != n.years)


###########################################################################
# initial
###########################################################################
z.init <- matrix(NA, nrow(ch), ncol(ch))
for (i in 1:nrow(ch)){
  if(first[i] < ncol(z.init)){
    z.init[i,(first[i] + 1):ncol(z.init)] <- 1
  }
}


################################################################
# RD model
################################################################
require(jagsUI)
sink(file="ms_rd.jags")
cat("
    model {
    
    ######################################################################################
    # survival and breeding probability
    ######################################################################################
    phi ~ dbeta(1,1)    
    gP ~ dbeta(1,1) 
    gPP ~ dbeta(1,1) 
    mean.p ~ dbeta(1,1)
    
    ######################################################################################
    # Secondary occasions p's
    ######################################################################################
    for (t in 1:n.years){
    for (j in 1:n.sec[t]){
    p[t,j] <- mean.p
    yes[t,j] ~ dbin(p[t,j], total[t,j])
    }
    }   
    
    ######################################################################################
    # Primary occasions p's
    ######################################################################################    
    for (t in 1:n.years){
    pstar[t] <- 1 - prod(1 - p[t,])
    }
    
    ###########################################################################
    # matrices
    ###########################################################################
    s[1,1] <- phi * gPP
    s[1,2] <- phi * (1 - gPP)
    s[1,3] <- 1 - phi
    s[2,1] <- phi * gP
    s[2,2] <- phi * (1 - gP)
    s[2,3] <- 1 - phi
    s[3,1] <- 0
    s[3,2] <- 0
    s[3,3] <- 1

    for (t in 1:n.years){
      o[1,t,1] <- pstar[t]
      o[1,t,2] <- 1 - pstar[t]
      o[2,t,1] <- 0
      o[2,t,2] <- 1
      o[3,t,1] <- 0
      o[3,t,2] <- 1
    }


    ###########################################################################
    # likelihood
    ###########################################################################
    for (i in 1:n.ind){
    
    z[i,first[i]] <- ch[i,first[i]]

    for (t in (first[i]+1):n.years){

    z[i,t] ~ dcat(s[z[i,t-1], ])
    ch[i,t] ~ dcat(o[z[i,t], t, ])

    } 
    
    }
    
    ####################################################################################
    # end model
    ####################################################################################
    
    }
    
    ", fill=TRUE)
sink()


######################################################################################
# bundle data
######################################################################################
dat <- list(first = first, ch = ch, 
            n.sec = rep(n.sec.occ,n.years), n.years = ncol(ch), n.ind = nrow(ch),
            yes = yes, no = no, total = total)


######################################################################################
# Initial Values
######################################################################################
inits <- function(){list(z = z.init)}  


######################################################################################
# parameters to monitor
######################################################################################
pars <- c('pstar','mean.p','phi','gP','gPP')


######################################################################################
# MCMC settings
######################################################################################
n.chains <- 1
n.thin <- 2
n.adapt <- 100
n.iter <- 10000
n.burnin <- 7000


######################################################################################
# compile model! 
######################################################################################
msrd <- jags(dat, inits, pars, "ms_rd.jags", 
             n.chains = n.chains, n.adapt = n.adapt, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)

summary(msrd)


truth[k,1] <- s[1]
truth[k,2] <- gP[1]
truth[k,3] <- gPP[1]
truth[k,4] <- mean(p[,n.sec.occ+1])
estimate[k,1] <- msrd$mean$phi
estimate[k,2] <- msrd$mean$gP
estimate[k,3] <- msrd$mean$gPP
estimate[k,4] <- msrd$mean$pstar[1]
estimateL[k,1] <- msrd$q2.5$phi
estimateL[k,2] <- msrd$q2.5$gP
estimateL[k,3] <- msrd$q2.5$gPP
estimateL[k,4] <- msrd$q2.5$pstar[1]
estimateH[k,1] <- msrd$q97.5$phi
estimateH[k,2] <- msrd$q97.5$gP
estimateH[k,3] <- msrd$q97.5$gPP
estimateH[k,4] <- msrd$q97.5$pstar[1]

}

res.phi <- estimate[,1] - truth[,1]; mean(res.phi)
res.gP <- estimate[,2] - truth[,2]; mean(res.gP)
res.gPP <- estimate[,3] - truth[,3]; mean(res.gPP)
res.pstar <- estimate[,4] - truth[,4]; mean(res.pstar)


truth[,2] <- 1 - truth[,2]
estimate[,2] <- 1 - estimate[,2]
estimateL[,2] <- 1 - estimateL[,2]
estimateH[,2] <- 1 - estimateH[,2]
truth[,3] <- 1 - truth[,3]
estimate[,3] <- 1 - estimate[,3]
estimateL[,3] <- 1 - estimateL[,3]
estimateH[,3] <- 1 - estimateH[,3]

#######################################
# PLOTS
#######################################
bitmap("W:\\\\Ph_D\\\\Manuscripts\\\\robust_design\\\\MEE_submission\\\\Resubmission\\\\in_press_version\\\\Figure3.tiff",
       width = 12, height = 12, units = 'in', res = 300, type = 'tifflzw')
par(mar = c(5.1,6.1,2.1,2.1), mfrow = c(2,2))
plot(estimate[,1] ~ truth[,1], xlim = c(0.5,1), ylim = c(0.5,1),
     ylab = expression(phi['estimate']), xlab = expression(phi['truth']),
     pch = 19, cex.lab = 2, cex = 0.75)
arrows(truth[,1], estimateL[,1], truth[,1], estimateH[,1], angle = 90, code = 3, length = 0, lty = 2)
abline(0,1)

plot(estimate[,3] ~ truth[,3], xlim = c(0.5,1), ylim = c(0.5,1),
     ylab = expression(gamma[estimate]^"''"), xlab = expression(gamma['truth']),
     pch = 19, cex.lab = 2, cex = 0.75)
arrows(truth[,3], estimateL[,3], truth[,3], estimateH[,3], angle = 90, code = 3, length = 0, lty = 2)
abline(0,1)

plot(estimate[,2] ~ truth[,2], xlim = c(0,1), ylim = c(0,1),
     ylab = expression(gamma['estimate']^"'"), xlab = expression(gamma['truth']),
     pch = 19, cex.lab = 2, cex = 0.75)
arrows(truth[,2], estimateL[,2], truth[,2], estimateH[,2], angle = 90, code = 3, length = 0, lty = 2)
abline(0,1)

plot(estimate[,4] ~ truth[,4], xlim = c(0.45,0.75), ylim = c(0.45,0.75),
     ylab = expression("p*"['estimate']), xlab = expression("p*"['truth']),
     pch = 19, cex.lab = 2, cex = 0.75)
arrows(truth[,4], estimateL[,4], truth[,4], estimateH[,4], angle = 90, code = 3, length = 0, lty = 2)
abline(0,1)
dev.off()







hist(res.phi, main = NULL, col = 'royalblue1', breaks = 10, xlim = c(-0.15,0.15), xlab = expression())
hist(res.gPP, main = NULL, col = 'royalblue1', breaks = 10, xlim = c(-0.2,0.2))
hist(res.gP, main = NULL, col = 'royalblue1', breaks = 10, xlim = c(-0.6,0.6))
hist(res.pstar, main = NULL, col = 'royalblue1', breaks = 10, xlim = c(-0.1,0.1))










