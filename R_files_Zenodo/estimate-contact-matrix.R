# This routine realises MCMC samples of the mean and overdispersion parameters 
# in the negative binomial model for the entries in a 8x8 contact matrix. 
# The mean parameters are m_{jk} and the overdispersion parameters theta_{jk}. 
# Here j=1,..,8 (contactee age class) and k=1,..,8 (respondent age class). 
# Age classes are 0-9, 10-19, ..., 60-69, and 70+ year. Regarding the parameterisation 
# of the negative binomial model as used in this programme, see Auranen et al. "Social 
# distancing and SARS-CoV-2 transmission potential during early epidemic in Finland"; section 
# "Estimation of contact matrices".
#
# The main routine follows after the subroutines. MCMC samples are realised in subroutine runMCMC.R. 
# The length of the MCMC samples is determined by two parameters: n.mcmc and n.burnin (see the main
# routine. The parameter (la.theta) of the exponential prior distribution of the overdispersion 
# parameter is also set within the code below. 
#
# The program prints the posterior means of the mean and overdispersion parameters
# for each entry in the 8x8 matrix.
#
# Last updated February 15, 2021.

library(MASS)

lognegbin = function(j,k,cjk,ckj,mjk,thet.k,thet.j,Nk,Nj,paino.k,paino.j){
# Negative binomial log likelihood function
#
# INPUT cjk = observed counts from age class k to age class j 
#             a vector of the same length as there are respondents of age class k
#       ckj = observed counts from age class j to age class k, 
#             a vector of the same length as there are respondents of age class j 
#       mjk = the mean number of contacts to age class j by an individual of age class k
#       thet.k = overdispersion parameter (contacts k->j)
#       thet.j = overdispersion parameter (contacts j->k)
  
#       Nk     = population size of age class k
#       Nj     = population size of age class j
#       paino.k  = weight of individual likelihood contribution
#       paino.j  = weight of individual likelihood contribution 

  # Numbers of observations
  nobs.k = length(cjk) # respondents of class k
  nobs.j = length(ckj) # respondents of class j
  
  # Initialise the log likelihood
  loglike = 0
  
  ###############
  # Contacts k->j
  ###############
  
  # "success probability" of the neg. bin distribution in terms of
  # mean and overdispersion parameters
  prob = (1/thet.k)/((1/thet.k) + mjk)
  
  # Iterate over the individual contributions
  for (i in 1:nobs.k){
      loglike = loglike + paino.k[i]*(lgamma(cjk[i]+1/thet.k) - lgamma(cjk[i]+1) - lgamma(1/thet.k) +
                       (1/thet.k)*log(prob) + cjk[i]*log(1-prob))
  }
  
  ###############
  # Contacts j->k
  ###############
  if (j!=k){
    
    #
    prob =  (1/thet.j)/((1/thet.j) + mjk*Nk/Nj)
    
    for (i in 1:nobs.j){
      loglike = loglike + paino.j[i]*(lgamma(ckj[i]+1/thet.j) - lgamma(ckj[i]+1) - lgamma(1/thet.j) +
                        (1/thet.j)*log(prob) + ckj[i]*log(1-prob))
     }
  } 
  
  return(loglike)  
}


propose.par = function(oldpar,delta.par){
# A function to make a symmetric proposal in the Metropolis-Hastings algorithm
# INPUT oldpar    = current value of the a parameter
#       delta.par = length of the proposal window (centered at oldpar) 
    newpar = oldpar + (runif(1,0,1)-0.5)*delta.par
    return(newpar)
}


mcmc.negbin = function(j,k,cjk,ckj,Nk,Nj,n.mcmc,la.theta,paino.k,paino.j){
# The main MCMC sampling routine
#  
# INPUT cjk = observed counts from age class k to age class j 
#             a vector of the same length as there are respondents of age class k
#       ckj = observed counts from age class j to age class k, 
#             a vector of the same length as there are respondents of age class j 
#       Nk     = population size of age class k
#       Nj     = population size of age class j
#       n.mcmc = number of MCMC iterations  
#       la.theta = parameter of the exponential prior distribution of the inverse
  #               of the overdispersion parameter
#       paino.k = weights for individuals in age class k
#       paino.j = weights for individuals in age class j
  
  
  # Reserve space for the MCMC samples
  mjk    = rep(0,n.mcmc)
  thet.k = rep(0,n.mcmc)
  thet.j = rep(0,n.mcmc)
  
  # Proposal windows
  m.delta    = 0.8
  thet.delta = 0.8
  
  # Draw random initial values 
  mjk[1]    = runif(1,0.01,4)
  thet.k[1] = runif(1,0.04,10)
  thet.j[1] = runif(1,0.04,10)
  
  # Calculate the log-likelihood with the initial values
  cur.like = lognegbin(j,k,cjk,ckj,mjk[1],thet.k[1],thet.j[1],Nk,Nj,paino.k,paino.j)

  # MCMC iterations
  for (nn in 2:n.mcmc){
  
    
    # Default
    mjk[nn] = mjk[nn-1]
    
    # Propose a new value
    m.propose = propose.par(mjk[nn-1],m.delta)
    
    if (m.propose >0){
      
      # The log-likehood under the proposed parameter vector
      new.like = lognegbin(j,k,cjk,ckj,m.propose,thet.k[nn-1],thet.j[nn-1],Nk,Nj,paino.k,paino.j)
      
      # Calculate the log acceptance ratio
      AR = new.like - cur.like 
      
      # Metropolis Hastings step
      if (log(runif(1)) < AR) {
        mjk[nn]  = m.propose
        cur.like = new.like
        }
     }
    
    ###############
    # Update thet.k
    ###############
    
    # Default
    thet.k[nn] = thet.k[nn-1]
    
    # Propose a new value
    t.propose = thet.k[nn-1]*exp(thet.delta*(runif(1,0,1)-0.5))

    if (t.propose >0){
      
      # The log-likehood under the proposed parameter vector
      new.like = lognegbin(j,k,cjk,ckj,mjk[nn],t.propose,thet.j[nn-1],Nk,Nj,paino.k,paino.j)
      
      # Calculate the log acceptance ratio
      AR = new.like - cur.like -log(thet.k[nn-1]) + log(t.propose)  - la.theta*t.propose + la.theta*thet.k[nn-1]
      
      # Metropolis Hastings step
      if (log(runif(1)) < AR) {
        thet.k[nn] = t.propose
        cur.like  = new.like
      }
    }
    
    ###############
    # Update thet.j
    ###############
    
    # Default
    thet.j[nn] = thet.j[nn-1]
    
    # Propose a new value
    t.propose = thet.j[nn-1]*exp(thet.delta*(runif(1,0,1)-0.5))
    
    if (t.propose >0){
      
      # The log-likehood under the proposed parameter vector
      new.like = lognegbin(j,k,cjk,ckj,mjk[nn],thet.k[nn],t.propose,Nk,Nj,paino.k,paino.j)
      
      # Calculate the log acceptance ratio
      AR = new.like - cur.like - log(thet.j[nn-1]) + log(t.propose) - la.theta*t.propose + la.theta*thet.j[nn-1]
      
      # Metropolis Hastings step
      if (log(runif(1)) < AR) {
        thet.j[nn] = t.propose
        cur.like  = new.like
      }
    }
    
 
    
  }
  return(list(mjk,thet.k,thet.j))  
  
}

# 
runMCMC = function(contactdata,aged10,n.mcmc,n.burnin,la.theta){
# This is the main MCMC routine  

#
n.samples    = n.mcmc-n.burnin                   # number of samples to be retained
samples.mean = matrix(0,8*8,(2+n.samples))       # MCMC samples of the mean matrix elements
samples.overdispersion = matrix(0,8*8,(2+n.samples)) # MCMC samples of the overdispersion parameters
  
hh = 1
# 
for (j in 1:8){ 
  
  # print(j)

  for (k in 1:j){
  
  # print(k) 
   
  # Select the numbers of reported contacts by individuals in age class k
  k.ind = which(contactdata[,1] == k) # row indices for age class k 
  cjk = contactdata[k.ind,j+1]        # reported contacts from k->j
  weig.k = contactdata[k.ind,"weight"]
  
  j.ind = which(contactdata[,1] == j) # row indices for age class j   
  ckj = contactdata[j.ind,k+1] # reported contacts from j->k
  weig.j = contactdata[j.ind,"weight"]
     
  # Realise an MCMC sample of length n.mcmc
  param = mcmc.negbin(j,k,cjk,ckj,aged10[k],aged10[j],n.mcmc,la.theta,weig.k,weig.j)
 
  # Extract the MCMC output
  param1 = param[[1]][(n.burnin+1):n.mcmc] # mean parameter m_jk
  param2 = param[[2]][(n.burnin+1):n.mcmc] # theta_jk = overdispersion parameter
  param3 = param[[3]][(n.burnin+1):n.mcmc] # theta_kj = overdispersion parameter

  # Store samples
  samples.mean[hh,1]    = j
  samples.mean[hh,2]    = k
  samples.mean[hh,3:(2+n.samples)] = param1 
  samples.overdispersion[hh,1]    = j
  samples.overdispersion[hh,2]    = k
  samples.overdispersion[hh,3:(2+n.samples)] = param2 
  hh = hh+1
  
  if (j !=k){
    samples.mean[hh,1] = k
    samples.mean[hh,2] = j
    samples.mean[hh,3:(2+n.samples)] = param1*aged10[k]/aged10[j] 
    samples.overdispersion[hh,1]    = k
    samples.overdispersion[hh,2]    = j
    samples.overdispersion[hh,3:(2+n.samples)] = param3
    hh = hh+1
    }

 } 
  
}
 
 return(list(samples.mean,samples.overdispersion))

}

##############################
# The main routine starts here
############################## 

library(MASS)

#
n.mcmc     = 25000  # number of MCMC samples
n.burnin   = 5000   # number of burn-in samples 
n.samples  = n.mcmc-n.burnin # number of samples to retained after discarding the burn-in samples
la.theta   = 1 # Parameter of the exponential prior distribution of the overdispersion parameter

# Read the contact data (N.B. This is a simulated example dataset.)
contactdata = read.table("example-data.txt",header=TRUE)

# Number of participants
n.ind       = dim(contactdata)[1]

# Population sizes by 10-year age classes in Finland in 2019 
# (original source https://findikaattori.fi/fi/14).
# The age classes are 0-9, 10-19, 20-29, ..., 60-69 and 70+ years. 
aged10 = as.numeric(read.table("ageclass-sizes-FIN2020.txt",header=TRUE))

##################################
# Call the main estimation routine
##################################
out1 = runMCMC(contactdata,aged10,n.mcmc,n.burnin,la.theta) 

#######################################################################
# Based on the MCMCm sample, calculate the posterior means and 95% 
# posterior intervals of the mean parameters m_jk, j=1,...,8; k=1,...,8. 
# Cf. eTable1 (all contacts) and eTable3 (physical contacts).
#######################################################################
m.out = matrix(0,64,5)
for (i in 1:64){
  
  m.out[i,1] = out1[[1]][i,1] # contactee age class
  m.out[i,2] = out1[[1]][i,2] # respondent age class
  
  # Sort the MCMC samples in ascending order
  m.sample = sort(out1[[1]][i,3:(n.samples+2)])
  
  m.out[i,3] = mean(m.sample)                   # posterior mean 
  m.out[i,4] = m.sample[round(0.025*n.samples)] # 2.5% posterior quantile
  m.out[i,5] = m.sample[round(0.975*n.samples)] # 97.5% posterior quantile
}

# Print (row index j, column index k, posterior mean and 2.5% and 97.% posterior quantiles)
for (j in 1:8){
  for (k in 1:8){
    
    ii = which(m.out[,1]==j & m.out[,2]==k)
    print(round(m.out[ii,],2))
  }
  print("")
}

#############################################################################################
# Based on the MCMCm sample, calculate, calculate the posterior means and 95% intervals of the 
# overdispersion parameters theta_jk, j=1,...,8; k=1,...,8.
# Cf. eTable2  (all contacts) and eTable4 (physical contacts).
#############################################################################################
theta.out = matrix(0,64,5)
for (i in 1:64){
  
  theta.out[i,1] = out1[[2]][i,1] # contactee age class
  theta.out[i,2] = out1[[2]][i,2] # respondent age class
  
  # Sort the MCMC samples in ascending order
  theta.sample = sort(out1[[2]][i,3:(n.samples+2)])
  
  theta.out[i,3] = mean(theta.sample)                   # posterior mean
  theta.out[i,4] = theta.sample[round(0.025*n.samples)] # 2.5% posterior quantile
  theta.out[i,5] = theta.sample[round(0.975*n.samples)] # 97.5% posterior quantile
  
}    

# Print (row index j, column index k, posterior mean and 2.5% and 97.% posterior quantiles)
for (j in 1:8){
  for (k in 1:8){
    
    ii = which(theta.out[,1]==j & theta.out[,2]==k)
    print(round(theta.out[ii,],2))
  }
  print("")
}





