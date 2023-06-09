#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#      Comparing traditional and Bayesian approaches to ecological meta-analysis
#                                   Pappalardo et al. 
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# All functions were written by Paula Pappalardo, with assistance from Kiona Ogle, Elizabeth Hamman and Jim Bence

#------------------------Functions to simulate datasets--------------------------------------
#--------------------------------------------------------------------------------------------

makeLognormalRow <- function(mu, tau, sigma, mean.n, c.mean){
  # simulates a row of lognormal data representing study i
  #
  # Args:
  #   mu: numeric value for the overall effect size
  #   tau: numeric value for the among-study variance 
  #   sigma: numeric value for the among-replicates variation
  #   mean.n: numeric value for number of replicates
  #   c.mean: numeric value for the control mean
  # Returns:
  #   Dataframe with a row of simulated data for a study in the meta-analysis
  #
  # Simulate number of replicates using a Poisson distribution
  # To avoid the 0 and 1 we do a small trick, and rescaled all by 2              
  pn <- rpois(1, lambda = mean.n - 2) # get a number
  n <- pn + 2 # add 2 to make the n for the experiment
  
  # calculate a RANDOM ERROR based on tau
  eta <- rnorm(1, mean = 0, sd = sqrt(tau))
  
  # simulate raw data for control and treatment          
  c <- rlnorm(n, meanlog = log(c.mean) - (sigma^2)/2, sdlog = sigma)   
  t <- rlnorm(n, meanlog = log(c.mean) - (sigma^2)/2 + mu + eta , sdlog = sigma)  
  
  # calculate log ratio and its variance
  lnrr <- log(mean(t)/mean(c))   # lnRR effect size
  var.lnrr <- sd(t)^2/(n*mean(t)^2) + sd(c)^2/(n*mean(c)^2)   # lnRR variance
  
  # calculate standard error to use in true Bayesian analysis
  se.t <- sd(t)/sqrt(n)
  se.c <- sd(c)/sqrt(n)
  
  # put results in dataframe
  results <- data.frame(lnrr = lnrr, var.lnrr = var.lnrr, yt = mean(t), yc = mean(c),
                        nt = n, nc = n, sd.t = sd(t), sd.c = sd(c), se.t = se.t, se.c = se.c)
  return(results)
}

makeLognormalDataset <- function(mu, tau, sigma, mean.n, c.mean, k){
  # create datasets for simulations
  #
  # Args:
  #   mu: numeric value for the overall effect size
  #   tau: numeric value for the among-study variance 
  #   sigma: numeric value for the among-replicates variation
  #   mean.n: numeric value for number of replicates
  #   c.mean: numeric value for the control mean
  #   k: number of studies
  # Returns:
  #   list of meta-analytic datasets for all factor combinations
  
  datasets <- list() # empty list to hold results
  
  # Loop through all combinations of parameters
  for (m in mu) {
    for (i in sigma) {
      for (j in tau) {
        for (n in mean.n) {
          for (z in k) {
            # set a name to put things in a list
            name <- paste("sigma", i, "tau", j, "mean.n", n, "k", z, "mu", m, "c.mean", c.mean, sep = " ")
            
            rows.list <- list() # make empty list to put each row
            
            # creating a dataset with k= z
            for (x in 1:z) {
              namek <- paste("k number", x)
              # get dataset for this particular combination of variables
              datarow <- makeLognormalRow(sigma = i, tau = j, mu = m, mean.n = n, c.mean = c.mean)
              # put row in list 
              rows.list[[namek]] <- datarow
            }
            df.results <- do.call(rbind, rows.list)
            df.results$true.mu <- m
            df.results$true.tau2 <- j
            df.results$mean.n <- n
            df.results$k <- z
            df.results$sigma <- i
            df.results$c.mean <- c.mean
            row.names(df.results) <- NULL
            datasets[[name]] <- df.results
          }
        }
      }
    }
  }
  # return list with datasets
  return(datasets)
}


#------------ ----Functions to run meta-analysis with traditional methods--------------------
#--------------------------------------------------------------------------------------------
runMetafor <- function(mydf, tau.method, ci){
  # Run a weighted random effects meta-analysis with Metafor
  #
  # Args:
  #   tau.method: defines the method to estimate the among-study variance
  #   ci: T to run Knapp-Hartung CI, F for default z-distribution CI
  # Returns:
  #   Summary results from the meta-analysis
  #
  library(metafor)
  
  # run meta-analysis
  res <- try(rma(lnrr, var.lnrr, data = mydf, knha = ci,
                 method = tau.method, control = list(maxiter = 200)))
  class(res)[1] -> error # to check if metafor gave an error
  
  # prepare dataframe in case of error
  if (error == "try-error") {
    # make dataframe to hold the data if there was an error
    res2 <- data.frame(error = 1)
    #print("this one didn't converge")
    return(res2)
  } else {
    # get confidence intervals for tau
    res.tau <- confint(res)
    
    # extract CI's for tau
    tau.ci.lb <-  res.tau$random["tau^2", "ci.lb"]
    tau.ci.ub <-  res.tau$random["tau^2", "ci.ub"]
    
    # make dataframe to hold answers in case there is no error
    res2 <- data.frame(obs.mu = res$b, mu.ci.lb = res$ci.lb, mu.ci.ub = res$ci.ub,
                       se = res$se, obs.tau = res$tau2, tau.ci.lb, tau.ci.ub,
                       true.mu = mydf$true.mu[1], sigma = mydf$sigma[1],
                       true.tau = mydf$true.tau2[1], k = mydf$k[1], mean.n = mydf$mean.n[1],
                       error = 0, i2 = res$I2)
    
    # calculate bias for overall effect and tau
    res2$bias.eff <- res2$obs.mu - res2$true.mu
    res2$bias.tau <- res2$obs.tau - res2$true.tau
    
    # calculate coverage by comparing the metafor CIs with the real answer
    coverage.mu <- ifelse(res2$mu.ci.lb <= res2$true.mu & res2$mu.ci.ub >= res2$true.mu, 1, 0)
    res2$cov.mu <- coverage.mu # add coverage to the final dataframe
    coverage.tau <- ifelse(res2$tau.ci.lb <= res2$true.tau & res2$tau.ci.ub >= res2$true.tau, 1, 0)
    res2$cov.tau <- coverage.tau # add coverage to the final dataframe
    
    # calculate the width of the CI
    res2$mu.width <- res2$mu.ci.ub - res2$mu.ci.lb
    res2$tau.width <- res2$tau.ci.ub - res2$tau.ci.lb
    
    # return results
    print("meta-analysis with metafor done")
    return(res2)
  }
} 

runMetafor_bootfunction <- function(mydf, indices){
  # function to bootstrap a meta-analysis with metafor
  #
  # Args:
  #   mydf: dataset to analyze
  #   indices: vector of indices which define the bootstrap sample
  # Returns:
  #   Summary results from the meta-analysis
  #
  library(metafor)
  
  res <- try(rma(lnrr, var.lnrr, data = mydf, method = "REML", subset = indices), silent = TRUE)  
  if (is.element("try-error", class(res))) {
    c(NA)
  } else {
    c(coef(res), res$tau2)
  }  
}

bootMetafor <- function(mydf){
  # Bootstrap the mean effect and associated CI from meta-analysis
  #
  # Args:
  #   mydf: dataset to analyze
  # 
  # Returns:
  #   Summary results from the meta-analysis
  #
  library(boot)
  library(metafor)
  
  # function for bootstrapping CI's
  res.boot <- try(boot(mydf, runMetafor_bootfunction, R = 5000))
  
  # use try in case function fails
  class(res.boot) -> error.1
  if (error.1 == "try-error") {
    res.error <- data.frame(obs.mu = NA, mu.ci.lb = NA, mu.ci.ub = NA)
    return(res.error)
  } else {
    # run the boot.ci function to get the various non-parametric bootstrap CI's
    res.muci <- try(boot.ci(res.boot, index = 1))
    class(res.muci) -> error.2
    if (error.2 == "try-error") {
      res.error <- data.frame(obs.mu = NA, mu.ci.lb = NA, mu.ci.ub = NA)
      return(res.error)
    } else {
      # select the bias corrected non parametric CIs
      mu.ci.lb <- res.muci$bca[4]
      mu.ci.ub <- res.muci$bca[5]
      obs.mu <- res.muci$t0[1]
      
      # make dataframe to hold answers in case there is no error
      res1 <- data.frame(obs.mu,  mu.ci.lb, mu.ci.ub, true.mu=mydf$true.mu[1], sigma= mydf$sigma[1], true.tau= mydf$true.tau2[1], k=mydf$k[1], mean.n= mydf$mean.n[1], error=0, i2=res.boot$t0[5])
      
      # calculate bias for overall effect and tau
      res1$bias.eff <- res1$obs.mu - res1$true.mu
      
      # calculate coverage by comparing the metafor CIs with the real answer
      
      coverage.mu <- ifelse(res1$mu.ci.lb <= res1$true.mu & res1$mu.ci.ub >= res1$true.mu,1,0)
      res1$cov.mu<- coverage.mu # add coverage to final dataframe
      
      # calculate the width of the CI
      res1$mu.width <- res1$mu.ci.ub - res1$mu.ci.lb 
      
      print("weighted meta-analysis boostrap CI with Metafor")
      return(res1)
    }
  }
}

# function for bootstrapping CI's for the among-study variance
bootMetaforTau <- function(mydf){
  # Bootstrap the among-study variance and associated CI from meta-analysis
  #
  # Args:
  #   mydf: dataset to analyze
  # 
  # Returns:
  #   Summary results from the meta-analysis
  #
  # function for bootstrapping CI's
  res.boot <- try(boot(mydf, runMetafor_bootfunction, R = 2000))
  # use try in case function fails
  class(res.boot) -> error.1
  if (error.1 == "try-error") {
    res.error <- data.frame(obs.tau = NA, tau.ci.lb = NA, tau.ci.ub = NA)
    return(res.error)
  } else {
    # run the boot.ci function to get the various non-parametric bootstrap CI's
    res.tauci <- try(boot.ci(res.boot, index = 2))
    res.muci <- try(boot.ci(res.boot, index = 1))
    class(res.tauci) -> error.2
    class(res.muci) -> error.3
    if (error.2 == "try-error" | error.3 == "try-error") {
      res.error <- data.frame(obs.tau = NA, tau.ci.lb = NA, tau.ci.ub = NA)
      return(res.error)
    } else {
      tau.ci.lb <- res.tauci$bca[4]
      tau.ci.ub <- res.tauci$bca[5]
      obs.tau <- res.tauci$t0[1]
      
      mu.ci.lb <- res.muci$bca[4]
      mu.ci.ub <- res.muci$bca[5]
      obs.mu <- res.muci$t0[1]
    }
    
    # make dataframe to hold answers in case there is no error
    res1 <- data.frame(obs.mu, mu.ci.lb, mu.ci.ub, obs.tau, tau.ci.lb, tau.ci.ub,
                       true.mu = mydf$true.mu[1], sigma = mydf$sigma[1],
                       true.tau = mydf$true.tau2[1], k = mydf$k[1], mean.n = mydf$mean.n[1],
                       error = 0, i2 = res.boot$t0[5])
    
    # calculate bias for overall effect and tau
    res1$bias.eff <- res1$obs.mu - res1$true.mu
    res1$bias.tau <- res1$obs.tau - res1$true.tau
    
    # calculate coverage by comparing the metafor CIs with the real answer
    coverage.mu <- ifelse(res1$mu.ci.lb <= res1$true.mu & res1$mu.ci.ub >= res1$true.mu,1,0)
    coverage.tau <- ifelse(res1$tau.ci.lb <= res1$true.tau & res1$tau.ci.ub >= res1$true.tau,1,0)
    res1$cov.mu <- coverage.mu # add coverage to final dataframe
    res1$cov.tau <- coverage.tau # add coverage to final dataframe
    
    # calculate the width of the CI
    res1$mu.width <- res1$mu.ci.ub - res1$mu.ci.lb 
    res1$tau.width <- res1$tau.ci.ub - res1$tau.ci.lb 
    
    print("weighted meta-analysis boostrap CI for tau2 with Metafor")
    return(res1)
  }
}


#-----------------------------------Functions To calculate CI ----------------------------------
#-----------------------------------------------------------------------------------------------

# Confidence intervals for coverage (binomial variable)

bi.95.l <- function(vec1){
  vec <- vec1[!is.na(vec1)]
  s <- length(which(vec == 1))
  runs <- length(vec)
  res <- binom.confint(x = s, n = runs, conf.level = 0.95, method = c("wilson"))
  res.low <- res$mean - res$lower
  return(res.low)
}

bi.95.u <- function(vec1){
  vec <- vec1[!is.na(vec1)]
  s <- length(which(vec == 1))
  runs <- length(vec)
  res <- binom.confint(x = s, n = runs, conf.level = 0.95, method = c("wilson"))
  res.up <- res$upper - res$mean
  return(res.up)
}

# Confidence intervals for bias

t.95CI <- function(vec){
  x <- vec[!is.na(vec)]
  n <- length(x)
  t <- qt(0.975, df = n - 1)
  se <- sd(x)/sqrt(n)
  ci <- t * se
  return(ci)
}

z.95CI <- function(vec){
  x <- vec[!is.na(vec)]
  n <- length(x)
  se <- sd(x)/sqrt(n)
  z <- qnorm(0.975) * se
  ci <- z * se
  return(ci)
}

#---------------------------Functions to run Bayesian meta-analysis--------------------------
#--------------------------------------------------------------------------------------------

# load libraries we need

library(rjags)
library(runjags)
library(coda)

# Create function to summarize data from the JAGS output
make_JAGS_summary <- function(sims, sigma, true.mu, true.tau, k, mean.n) {
  # summarize results from the JAGS function
  #
  # Args:
  #   sims: output of run.jags() function
  #   true.tau: numeric value for the true among-study variance
  #   true.mu: numeric value for the true overall effect
  #   sigma: numeric value for the among-replicates variation
  #   mean.n: numeric value for mean number of replicates
  #   k: numeric value for the number of setudies
  # Returns:
  #   Dataframe with the meta-analysis results
  #
  # Get summary of the results
  summ <- sims$summaries
  summ2 <- sims$summary
  obs.mean.mu <- summ["mu.theta", "Mean"]  # median of overall effect
  obs.mean.tau <- summ["var.theta", "Mean"]  # median of among studies
  obs.median.mu <- summ["mu.theta", "Median"]  # median of overall effect
  obs.median.tau <- summ["var.theta", "Median"]  # median of among studies variance
  obs.sd.mu <- summ["mu.theta", "SD"]  # SD overall mean
  obs.sd.tau <- summ["var.theta", "SD"] # SD among studies variance
  n.burnin <- sims$burnin  # number of samples used for burnin
  n.iter <- sims$sample  # number of iterations (without burn-in)
  n.chains <- summ2[["nchain"]]  # number of markov chains in the model
  
  # Asses convergence of chains, Kiona says it needs to be lower than 1.1
  Rhat.mu <- summ["mu.theta","psrf"]
  Rhat.tau2 <- summ["var.theta","psrf"]
  
  # get samples    
  mcmcsamp <- sims$mcmc
  mcmcmat <- do.call("rbind", mcmcsamp)
  
  # Convert simulations to MCMC object 
  mcmc.sims <- as.mcmc(mcmcmat)
  
  # Calculate 95% high density intervals
  hdi.sims <- HPDinterval(mcmc.sims, prob = 0.95)
  
  # HDI for the overall effect
  obs.mu.hdiL <- hdi.sims["mu.theta","lower"]
  obs.mu.hdiU <- hdi.sims["mu.theta","upper"]
  mu.width <- obs.mu.hdiU - obs.mu.hdiL # width of the high density interval
  
  # HDI for the among studies variance
  obs.tau.hdiL <- hdi.sims["var.theta","lower"]
  obs.tau.hdiU <- hdi.sims["var.theta","upper"]
  tau.width <- obs.tau.hdiU - obs.tau.hdiL # width of the high density interval
  
  # calculate coverage for the overall effect and tau
  cov.mu <- ifelse(obs.mu.hdiL <= true.mu & obs.mu.hdiU >= true.mu, 1, 0)
  cov.tau <- ifelse(obs.tau.hdiL <= true.tau & obs.tau.hdiU >= true.tau, 1, 0)
  
  # add bias estimation for the overall effect and for the among studies variance
  bias.mean.eff <- obs.mean.mu - true.mu
  bias.mean.tau <- obs.mean.tau - true.tau
  bias.median.eff <- obs.median.mu - true.mu
  bias.median.tau <- obs.median.tau - true.tau
  
  # Create dataframe with all the output we want
  bayes.results <- data.frame(true.mu, sigma, true.tau, k, mean.n, obs.mean.mu,
                              obs.median.mu, obs.sd.mu, obs.mu.hdiL, obs.mu.hdiU,
                              obs.mean.tau, obs.median.tau, obs.tau.hdiL, obs.tau.hdiU,
                              obs.sd.tau, bias.mean.eff, bias.mean.tau, bias.median.eff,
                              bias.median.tau,n.burnin, n.iter, n.chains, Rhat.mu,
                              Rhat.tau2, cov.mu, cov.tau,  mu.width, tau.width)
  
  # return data
  return(bayes.results)
}

#--------------------------Hybrid Bayesian model with uniform prior --------------------------
modelString = " 
model{
for (i in 1:N){
y[i] ~ dnorm (theta[i], tau.y[i])
theta[i] ~ dnorm (mu.theta, tau.theta)
}
# Setting priors
mu.theta ~ dnorm (0.0, 1.0E-3)   # prior for Overall effect
tau.theta <- pow(sigma.theta, -2)     # define precision
sigma.theta ~ dunif (0, 10)              # prior for standard deviation
var.theta <- pow(sigma.theta, 2)  # set variance to monitor
}
"
# write model file:
writeLines(modelString, con = "modelHalfBayes_uniform_JAGS.txt")

runHalfBayes_uniform_JAGS <- function(mydf){
  # runs Hybrid Bayesian meta-analysis with uniform prior
  #
  # Args:
  #   mydf: dataset with columns eff.size and var
  # Returns:
  #   One line dataframe with the meta-analysis results
  #
  # Organize data we need
  mydf <- mydf[complete.cases(mydf),]
  N <- nrow(mydf)  # number of studies
  y <- mydf$lnrr  # the observed effect sizes
  tau.y <- 1/mydf$var.lnrr  # the variance of the effect sizes
  muinit <- mean(mydf$lnrr)
  sdinit <- sd(mydf$lnrr)
  
  # Specify data in a list form
  datalist <- list("N" = N, "y" = y, "tau.y" = tau.y)
  
  # It is highly recommend to specify the initial values for the chains
  chain1 <- list(mu.theta = rnorm(1, muinit, sdinit), sigma.theta = runif(1, 0.8*sdinit, 1.2*sdinit))
  chain2 <- list(mu.theta = rnorm(1, muinit, sdinit), sigma.theta = runif(1, 0.8*sdinit, 1.2*sdinit))
  chain3 <- list(mu.theta = rnorm(1, muinit, sdinit), sigma.theta = runif(1, 0.8*sdinit, 1.2*sdinit))
  initslist <- list(chain1, chain2, chain3)
  
  # Start the MCMC simulation:
  sims <- run.jags(data = datalist, inits = initslist, model = "modelHalfBayes_uniform_JAGS.txt", 
                   n.chains = 3, burnin = 100000, sample = 100000,
                   monitor = c("mu.theta", "var.theta"), method = "rjags") 
  
  
  # Make dataframe with the results of the meta-analysis
  results <- make_JAGS_summary(sims, sigma = mydf$sigma[1], true.mu = mydf$true.mu[1],
                               true.tau = mydf$true.tau2[1], k = mydf$k[1], mean.n = mydf$mean.n[1])
  
  # Return dataframe and print "done"
  print("half bayesian meta-analysis with uniform prior done with JAGS")  
  return(results) 
  
  # remove things not in use
  rm(mydf, N, y, tau.y, muinit, sdinit, data, chain1, chain2, chain3, inits, sims, results)
}


# --------------------ALL HYBRID BAYESIAN MODELS--------------------------

# -------------------Folded N

modelString = "
model{
  for (i in 1:N){ # likelihood 
    y[i] ~ dnorm (theta[i], tau.y[i])
    theta[i] ~ dnorm (mu.theta, tau.theta)
  }
  # Setting priors
  mu.theta ~ dnorm (0.0, 1.0E-3)  # prior for overall effect
  sigma.norm ~ dnorm (0, 1) # folded N prior for standard deviation
  sigma.theta <- abs(sigma.norm) # here I folded it
  tau.theta <- pow(sigma.theta, -2) # set precision for the likelihood
  var.theta <- pow(sigma.theta, 2) #  set variance to monitor
}
"
# some temporary filename:
writeLines(modelString, con = "modelHalfBayes_folded_N.txt")

# Function to run half bayesian analysis with folded N prior
runHalfBayes_folded_N <- function(mydf){
  # runs Hybrid Bayesian meta-analysis with folded N prior
  #
  # Args:
  #   mydf: dataset with columns eff.size and var
  # Returns:
  #   One line dataframe with the meta-analysis results
  #
  # Organize data we need
  mydf <- mydf[complete.cases(mydf),]
  N <- nrow(mydf)  # number of studies
  y <- mydf$lnrr  # the observed effect sizes
  tau.y <- 1/mydf$var.lnrr  # the variance of the effect sizes
  muinit <- mean(mydf$lnrr)
  sdinit <- sd(mydf$lnrr)
  
  # Specify data in a list form
  datalist <- list("N" = N, "y" = y, "tau.y" = tau.y)
  
  # It is highly recommend to specify the initial values for the chains
  chain1 <- list(theta = rnorm(N, muinit, sdinit), mu.theta = rnorm(1, muinit, sdinit), sigma.norm = rnorm(1, 0, sdinit))
  chain2 <- list(theta = rnorm(N, muinit, sdinit), mu.theta = rnorm(1, muinit, sdinit), sigma.norm = rnorm(1, 0, sdinit))
  chain3 <- list(theta = rnorm(N, muinit, sdinit), mu.theta = rnorm(1, muinit, sdinit), sigma.norm = rnorm(1, 0, sdinit))
  initslist <- list(chain1, chain2, chain3)
  
  # Start the MCMC simulation:
  sims <- run.jags(data = datalist, inits = initslist, model = "modelHalfBayes_folded_N.txt", 
                   n.chains = 3, burnin = 100000, sample = 100000,
                   monitor = c("mu.theta", "var.theta"), method = "rjags")
  
  # Make dataframe with the results of the meta-analysis
  results <- make_JAGS_summary(sims, sigma = mydf$sigma[1], true.mu = mydf$true.mu[1],
                               true.tau = mydf$true.tau2[1], k = mydf$k[1], mean.n = mydf$mean.n[1])
  
  # Return dataframe and print "done"
  print("half bayesian meta-analysis with folded N prior done with JAGS")  
  return(results) 
  
  # remove things not in use
  rm(mydf, N, y, tau.y, muinit, sdinit, data, chain1, chain2, chain3, inits, sims, results)
}

# -------------------Cauchy

modelString = "
model{
  for (i in 1:N){ # Likelihood
    y[i] ~ dnorm (theta[i], tau.y[i])
    theta[i] <- mu.theta + alpha*eps[i]
    eps[i] ~ dnorm(0, tau.eps)
  }
  # Setting priors
  mu.theta ~ dnorm (0.0, 1.0E-3)   # prior for Overall effect
  # Cauchy prior for the standard deviation
  alpha ~ dnorm(0, 1)
  tau.eps~ dgamma(a, b)
  a <- 1/2
  b <- pow(AA, 2)/2
  AA <- 25
  sig.eps <- 1/sqrt(tau.eps)
  sigma.theta <- abs(alpha)/sig.eps
  var.theta <- pow(sigma.theta, 2) # set variance to monitor
}
"
# some temporary filename:
writeLines(modelString, con = "modelHalfBayes_Cauchy.txt")

# Function to run half bayesian analysis with uniform prior
runHalfBayes_cauchy <- function(mydf){
  # runs Hybrid Bayesian meta-analysis with Cauchy prior
  #
  # Args:
  #   mydf: dataset with columns eff.size and var
  # Returns:
  #   One line dataframe with the meta-analysis results
  #
  # Organize data we need
  mydf <- mydf[complete.cases(mydf),]
  N <- nrow(mydf)  # number of studies
  y <- mydf$lnrr  # the observed effect sizes
  tau.y <- 1/mydf$var.lnrr  # the variance of the effect sizes
  muinit <- mean(mydf$lnrr)
  sdinit <- sd(mydf$lnrr)
  
  # Specify data in a list form
  datalist <- list("N" = N, "y" = y, "tau.y" = tau.y)
  
  # It is highly recommend to specify the initial values for the chains
  chain1 <- list(alpha = rnorm(1), tau.eps = runif(1), mu.theta = rnorm(1, muinit, sdinit))
  chain2 <- list(alpha = rnorm(1), tau.eps = runif(1), mu.theta = rnorm(1, muinit, sdinit))
  chain3 <- list(alpha = rnorm(1), tau.eps = runif(1), mu.theta = rnorm(1, muinit, sdinit))
  initslist <- list(chain1, chain2, chain3)
  
  # Start the MCMC simulation:
  sims <- run.jags(data = datalist, inits = initslist, model = "modelHalfBayes_Cauchy.txt",
                   monitor = c("mu.theta", "var.theta"), n.chains = 3, method = "rjags",
                   burnin = 100000, sample = 100000)
  
  # Make dataframe with the results of the meta-analyis
  results <- make_JAGS_summary(sims, sigma = mydf$sigma[1], true.mu = mydf$true.mu[1],
                               true.tau = mydf$true.tau2[1], k = mydf$k[1], mean.n = mydf$mean.n[1])
  
  # Return dataframe and print "done"
  print("half Bayesian meta-analysis with Cauchy prior done")  
  return(results) 
  
  # remove things not in use
  rm(mydf, N, y, tau.y, muinit, sdinit, data, chain1, chain2, chain3, inits, sims, results)
}

# -------------------Uniform (0, 10)

modelString = "
model {
  for (i in 1:N){
    y[i] ~ dnorm (theta[i], tau.y[i])
    theta[i] ~ dnorm (mu.theta, tau.theta)
  }
  # Setting priors
  mu.theta ~ dnorm (0.0, 1.0E-3)   # prior for Overall effect
  tau.theta <- pow(sigma.theta, -2)     # define precision
  sigma.theta ~ dunif (0, 10)              # prior for standard deviation
  var.theta <- pow(sigma.theta, 2)  # set variance to monitor
}
"
# some temporary filename:
writeLines(modelString, con = "modelHalfBayes_uniform.txt")

# -------------------Uniform (0, 1)

modelString = "
model {
  for (i in 1:N){
    y[i] ~ dnorm (theta[i], tau.y[i])
    theta[i] ~ dnorm (mu.theta, tau.theta)
  }
  # Setting priors
  mu.theta ~ dnorm (0.0, 1.0E-3)   # prior for Overall effect
  tau.theta <- pow(sigma.theta, -2)     # define precision
  sigma.theta ~ dunif (0, 1)              # prior for standard deviation
  var.theta <- pow(sigma.theta, 2)  # set variance to monitor
}
"
# some temporary filename:
writeLines(modelString, con = "modelHalfBayes_uniform1.txt")

# ------------------Uniform (0, 100)
modelString = "
model {
  for (i in 1:N){
    y[i] ~ dnorm (theta[i], tau.y[i])
    theta[i] ~ dnorm (mu.theta, tau.theta)
  }
  # Setting priors
  mu.theta ~ dnorm (0.0, 1.0E-3)   # prior for Overall effect
  tau.theta <- pow(sigma.theta, -2)     # define precision
  sigma.theta ~ dunif (0, 100)              # prior for standard deviation
  var.theta <- pow(sigma.theta, 2)  # set variance to monitor
}
"
# some temporary filename:
writeLines(modelString, con = "modelHalfBayes_uniform100.txt")

# Function to run half bayesian analysis with uniform prior
runHalfBayes_uniform <- function(mydf, mymodel){
  # runs Hybrid Bayesian meta-analysis with Uniform prior
  #
  # Args:
  #   mydf: dataset with columns eff.size and var
  #   mymodel: name of jags model to use
  # Returns:
  #   One line dataframe with the meta-analysis results
  #
  # Organize data we need
  mydf <- mydf[complete.cases(mydf),]
  N <- nrow(mydf)  # number of studies
  y <- mydf$lnrr  # the observed effect sizes
  tau.y <- 1/mydf$var.lnrr  # the variance of the effect sizes
  muinit <- mean(mydf$lnrr)
  sdinit <- sd(mydf$lnrr)
  
  # Specify data in a list form
  datalist <- list("N" = N, "y" = y, "tau.y" = tau.y)
  
  # It is highly recommend to specify the initial values for the chains
  chain1 <- list(mu.theta = rnorm(1, muinit, sdinit), sigma.theta = runif(1, 0.8*sdinit, 1.2*sdinit))
  chain2 <- list(mu.theta = rnorm(1, muinit, sdinit), sigma.theta = runif(1, 0.8*sdinit, 1.2*sdinit))
  chain3 <- list(mu.theta = rnorm(1, muinit, sdinit), sigma.theta = runif(1, 0.8*sdinit, 1.2*sdinit))
  initslist <- list(chain1, chain2, chain3)
  
  # Start the MCMC simulation:
  sims <- run.jags(data = datalist, inits = initslist, model = mymodel,
                   monitor = c("mu.theta", "var.theta"), method = "rjags",
                   n.chains = 3, burnin = 100000, sample = 100000)
  
  # Make dataframe with the results of the meta-analyis
  results <- make_JAGS_summary(sims, sigma = mydf$sigma[1], true.mu = mydf$true.mu[1],
                               true.tau = mydf$true.tau2[1], k = mydf$k[1], mean.n = mydf$mean.n[1])
  
  # Return dataframe and print "done"
  print("half bayesian meta-analysis with uniform prior done")  
  return(results) 
  
  # remove things not in use
  rm(mydf, N, y, var, muinit, sdinit, data, chain1, chain2, chain3, inits, sims, results)
}

runHalfBayes_uniform1 <- function(mydf, mymodel){
  # runs Hybrid Bayesian meta-analysis with Uniform(0, 1)
  #
  # Args:
  #   mydf: dataset with columns eff.size and var
  #   mymodel: name of jags model to use
  # Returns:
  #   One line dataframe with the meta-analysis results
  #
  # Organize data we need
  mydf <- mydf[complete.cases(mydf),]
  N <- nrow(mydf)  # number of studies
  y <- mydf$lnrr  # the observed effect sizes
  tau.y <- 1/mydf$var.lnrr  # the variance of the effect sizes
  muinit <- mean(mydf$lnrr)
  sdinit <- sd(mydf$lnrr)
  
  # Specify data in a list form
  datalist <- list("N" = N, "y" = y, "tau.y" = tau.y)
  
  # It is highly recommend to specify the initial values for the chains
  chain1 <- list(mu.theta = rnorm(1, muinit, sdinit), sigma.theta = runif(1, 0, 1))
  chain2 <- list(mu.theta = rnorm(1, muinit, sdinit), sigma.theta = runif(1, 0, 1))
  chain3 <- list(mu.theta = rnorm(1, muinit, sdinit), sigma.theta = runif(1, 0, 1))
  initslist <- list(chain1, chain2, chain3)
  
  # Start the MCMC simulation:
  sims <- run.jags(data = datalist, inits = initslist, model = mymodel,
                   monitor = c("mu.theta", "var.theta"), method = "rjags",
                   n.chains = 3, burnin = 100000, sample = 100000)
  
  # Make dataframe with the results of the meta-analyis
  results <- make_JAGS_summary(sims, sigma = mydf$sigma[1], true.mu = mydf$true.mu[1],
                               true.tau = mydf$true.tau2[1], k = mydf$k[1], mean.n = mydf$mean.n[1])
  
  # Return dataframe and print "done"
  print("half bayesian meta-analysis with uniform prior (0, 1) done")  
  return(results) 
  
  # remove things not in use
  rm(mydf, N, y, tau.y, muinit, sdinit, data, chain1, chain2, chain3, inits, sims, results)
}

# -------------------Gamma

modelString = "
model {
  for (i in 1:N){ # likelihood
    y[i] ~ dnorm (theta[i], tau.y[i])
    theta[i] ~ dnorm (mu.theta, tau.theta)
  }
  # Setting priors
  mu.theta ~ dnorm (0.0, 1.0E-3)  # Prior for the overall effect
  tau.theta ~ dgamma(0.1, 0.1)    # Gamma prior for the precision
  sigma.theta <- 1/sqrt(tau.theta) # get the standard deviation
  var.theta <- pow(sigma.theta, 2)  # set variance to monitor
}
"
# some temporary filename:
writeLines(modelString, con = "modelHalfBayes_Gamma.txt")

# Function to run half bayesian analysis with gamma prior
runHalfBayes_gamma <- function(mydf){
  # runs Hybrid Bayesian meta-analysis with Gamma prior
  #
  # Args:
  #   mydf: dataset with columns eff.size and var
  # Returns:
  #   One line dataframe with the meta-analysis results
  #
  # Organize data we need
  mydf <- mydf[complete.cases(mydf),]
  N <- nrow(mydf)  # number of studies
  y <- mydf$lnrr  # the observed effect sizes
  tau.y <- 1/mydf$var.lnrr  # the variance of the effect sizes
  muinit <- mean(mydf$lnrr)
  sdinit <- sd(mydf$lnrr)
  tauinit <- 1/sdinit^2
  
  # Specify data in a list form
  datalist <- list("N" = N, "y" = y, "tau.y" = tau.y)
  
  # It is highly recommend to specify the initial values for the chains
  chain1 <- list(mu.theta = rnorm(1, muinit, sdinit), tau.theta = runif(1, 0.8*tauinit, 1.2*tauinit))
  chain2 <- list(mu.theta = rnorm(1, muinit, sdinit), tau.theta = runif(1, 0.8*tauinit, 1.2*tauinit))
  chain3 <- list(mu.theta = rnorm(1, muinit, sdinit), tau.theta = runif(1, 0.8*tauinit, 1.2*tauinit))
  initslist <- list(chain1, chain2, chain3)
  
  # Start the MCMC simulation:
  gamma <- run.jags(data = datalist, inits = initslist, model = "modelHalfBayes_Gamma.txt",
                    monitor = c("mu.theta", "var.theta"), n.chains = 3, burnin = 100000,
                    sample = 100000, method = "rjags")
  
  # Make dataframe with the results of the meta-analysis
  results <- make_JAGS_summary(gamma, sigma = mydf$sigma[1], true.mu = mydf$true.mu[1],
                               true.tau = mydf$true.tau2[1], k = mydf$k[1], mean.n = mydf$mean.n[1])
  
  # Return dataframe and print "done"
  print("half bayesian meta-analysis with gamma prior done")  
  return(results) 
  
  # remove things not in use
  rm(mydf, N, y, tau.y, muinit, sdinit, data, chain1, chain2, chain3, inits, sims, results)
}

# -------------------folded-t-

modelString = "
model{
  for (i in 1:N){ # likelihood
    y[i] ~ dnorm (theta[i], tau.y[i])
    theta[i] ~ dnorm (mu.theta, tau.theta)
  }
  # Setting priors
  mu.theta ~ dnorm (0.0, 1.0E-3)  # Overall effect prior
  # folded t prior
  A <- 5
  v <- 2
  B <- 1/(A*A) 
  t.theta ~ dt (0, B, v) 
  sigma.theta <- abs(t.theta)
  tau.theta <- pow(sigma.theta, -2) 
  # set precision for the likelihood and variance to monitor
  var.theta <- pow(sigma.theta, 2)
}
"
# write model file:
writeLines(modelString, con = "modelHalfBayes_folded_t.txt")

# Function to run bayesian analysis with Folded t prior
runHalfBayes_folded_t <- function(mydf){
  # runs Hybrid Bayesian meta-analysis with Folded t prior
  #
  # Args:
  #   mydf: dataset with columns eff.size and var
  # Returns:
  #   One line dataframe with the meta-analysis results
  #
  # Organize data we need
  mydf <- mydf[complete.cases(mydf),]
  N <- nrow(mydf)  # number of studies
  y <- mydf$lnrr  # the observed effect sizes
  tau.y <- 1/mydf$var.lnrr  # the variance of the effect sizes
  muinit <- mean(mydf$lnrr)
  sdinit <- sd(mydf$lnrr)
  
  # Specify data in a list form
  datalist <- list("N" = N, "y" = y, "tau.y" = tau.y)
  
  # It is highly recommend to specify the initial values for the chains
  chain1 <- list(mu.theta = rnorm(1, muinit, sdinit), t.theta = runif(1, 0.8*sdinit, 1.2*sdinit))
  chain2 <- list(mu.theta = rnorm(1, muinit, sdinit), t.theta = runif(1, 0.8*sdinit, 1.2*sdinit))
  chain3 <- list(mu.theta = rnorm(1, muinit, sdinit), t.theta = runif(1, 0.8*sdinit, 1.2*sdinit))
  initslist <- list(chain1, chain2, chain3)
  
  # Start the MCMC simulation:
  sims <- run.jags(data = datalist, inits = initslist, model = "modelHalfBayes_folded_t.txt",
                   monitor = c("mu.theta", "var.theta"), n.chains = 3, burnin = 100000,
                   sample = 100000, method = "rjags")
  
  # Make dataframe with the results of the meta-analysis
  results <- make_JAGS_summary(sims, sigma = mydf$sigma[1], true.mu = mydf$true.mu[1],
                               true.tau = mydf$true.tau2[1], k = mydf$k[1], mean.n = mydf$mean.n[1])
  
  # Return dataframe and print "done"
  print("half bayesian meta-analysis with folded t prior done")  
  return(results) 
  
  # remove things not in use
  rm(mydf, N, y, var, A, k, muinit, sdinit, data, chain1, chain2, chain3, inits, sims, results)
  
}


# -------------------------The following function was not written by the authors
# This function will allow for multiple plots that shared a legend, and is available from https://rpubs.com/sjackman/grid_arrange_shared_legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}
