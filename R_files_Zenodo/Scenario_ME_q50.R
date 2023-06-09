#####################################################################################################
# Script for IPM bias paper
#####################################################################################################
#
# lines 0010-0800: M. Schaub functions
# lines 0800-1000: marker-bias function
# lines 0900-1100: models
#
# lines 1200-1300: run models and store results
# lines 1300-1500: plot results
#
#####################################################################################################
#
# Function to create a population based on demographic parameters
#
# M.Schaub, 9.3.2016
#
# Features of the population:
# - female-based (reproduction is the number of juv. females)
# - includes demographic stochasticity
# - does NOT include density-dependence (just exponential population growth)
#
# Features of the simulation tool:
# - allows for a flexible number of age classes
#
# Input parameters:
# - Ni: vector with age(stage)-specific population size in the first year [length of vector: number of age(stage) classes, but without immigrants]. Note that the first age(stage) class refers to the 1-year old individuals, since the model is a pre-breeding model. The number of age(stage) classes must correspond to the number of age(stage) classes in the matrix of the survival probabilities.
# - phi: matrix with age(stage)- and time-specific survival probabilities [dimension of matrix: number of age(stage) classes x number of years-1]. The rows correspond to the age(stage) classes (first row are the new-borns), the columns to the years. Note that the number of age(stage) classes specified here will correspond to the number of age(stage) classes of the simulated population. Moreover, the number of colums defines the number of years
# - f: matrix with the age(stage)- and time-specific fecundity rates [dimension of matrix: number of age(stage) classes x number of years]. The first row corresponds to the fecundity of individuals in the first age(stage) class.
# - Im: vector with the number of immigrants in each year [length of vector: number of years]
# - seed: to make sure that the random generation works properly
#
# Output parameters:
# - IND: array with the life history of all individuals in the population
# - Nu: summary statistics
#
#####################################################################################################

create.population <- function(Ni = c(10, 10), phi = matrix(c(rep(0.3, 5), rep(0.55, 5)), ncol = 5, byrow = TRUE), f = matrix(c(rep(1.6, 6), rep(1.6, 6)), ncol = 6, byrow = TRUE), Im = rep(0, 6), seed = NA){
  
  if (is.na(seed) == TRUE) seed = runif(1, 0, 100)
  # set.seed(seed)
  
  T <- ncol(phi)            # Number of years
  mAge <- nrow(phi)         # Maximal number of age(stage) classes
  
  # 1. Expand the vital rate matrices, such that the number of age classes corresponds to the number of years
  PHI <- matrix(0, ncol = T, nrow = mAge + T)
  PHI[1:nrow(phi),] <- phi
  u <- mAge + T - nrow(phi)
  if (u > 0){
    for (j in 1:u){
      PHI[nrow(phi)+j,] <- phi[nrow(phi),]
    } # j
  } # if
  
  F <- matrix(0, ncol = T + 1, nrow = mAge + T)
  F[1:nrow(f),] <- f
  u <- mAge + T - nrow(f)
  if (u > 0){
    for (j in 1:u){
      F[nrow(f)+j,] <- f[nrow(f),]
    } # j
  } # if
  
  Nindex <- c(0, cumsum(Ni))
  
  # 2. Create a Leslie matrix to determine approximately how many individuals will be ever alive in the population 
  N <- matrix(data = NA, nrow = mAge, ncol = T + 1)
  N[,1] <- Ni
  A <- array(0, dim = c(mAge, mAge, T))
  for (t in 1:T){
    for (j in 1:mAge){
      A[1,j,t] <- F[j,t] * PHI[1,t]    # First row in Leslie matrix
    } # j
    for (j in 2:mAge){
      A[j,j-1,t] <- PHI[j,t]           # Subdiagonal
    } # j
    A[mAge, mAge, t] <- PHI[mAge, t]
  } # t       
  for (t in 1:T){
    N[,t+1] <- A[,,t]%*%N[,t] + matrix(c(rep(0, mAge-1), Im[t]), ncol = 1)
  } # t
  no.ani <- round(sum(N)*5)              # 5 times as many individuals that were alive
  
  # 3. Define array for each individual
  ind <- array(NA, dim = c(mAge + 4, T + 1, no.ani))   # information about [1-ye, 2-ye, ..., mAge-ye, Juv, Im, Rep, Dead]
  
  # 4. Simulate the fates of individuals already present at t = 1 (in different age classes) and their reproduction   
  # 4.a: Simulate survival of the individuals present at t=1
  # Initialize
  for (a in 1:mAge){
    if (Ni[a]==0) next
    for (i in (Nindex[a]+1):Nindex[a+1]){
      ind[a,1,i] <- 1
    } # i
  } # if
  
  # Simulate survival
  z <- numeric()
  for (a in 1:mAge){
    for (i in (Nindex[a]+1):Nindex[a+1]){
      for (t in 1:T){
        z[t] <- rbinom(1,1,PHI[a+t,t])
      } # t
      Z <- max(sum(cumprod(z)))
      if (Z==0) {
        ind[mAge+4,2,i] <- 1
        next
      } # if
      for (u in 1:Z){
        if (a+u < mAge){
          ind[a+u,u+1,i] <- 1
        } # if
        else {
          ind[mAge,u+1,i] <- 1
        } # else                
      } # u
      # Record year of death (if any)
      if (sum(z)==T) next
      else {
        D <- min(which(z==0))
        ind[mAge+4,D+1,i] <- 1
      } # else          
    } # i
  } # a
  
  # 3.b: Survival of immigrants (in all years, not just of immigrants present at t = 1)
  Nimindex <- c(0, cumsum(Im)) + max(Nindex)
  for (t in 1:(T+1)){
    if (Im[t]==0) next
    for (i in (Nimindex[t]+1):Nimindex[t+1]){
      ind[mAge+2,t,i] <- 1
    } # i
  } # t
  for (t in 1:T){
    if (Im[t]==0) next
    for (i in (Nimindex[t]+1):Nimindex[t+1]){
      z <- numeric()
      for (d in t:T){      
        z[d-t+1] <- rbinom(1,1,PHI[mAge+d,d])
      } # d
      Z <- max(sum(cumprod(z)))
      if (Z==0){
        ind[mAge+4,t+1,i] <- 1
        next
      } # if
      for (u in 1:Z){
        if (t+u <= (T+1)){
          ind[mAge,u+t,i] <- 1
        } # if
        else {next}                
      }  # u
      # Record year of death (if any)
      if (sum(z)==T-t+1) next
      else {
        D <- min(which(z==0))
        ind[mAge+4,D+t,i] <- 1
      } # else          
    } # i
  } # t
  
  # 3.c: Simulate reproduction of all already existing individuals
  for (i in 1:max(Nimindex)){
    for (t in 1:(T+1)){
      g <- which(!is.na(ind[c(1:mAge, mAge+2),t,i]))
      if (length(g)==0) next
      if (g !=8){
        ind[mAge+3,t,i] <- rpois(1,F[g,t])
      } # if
      if (g==8){
        ind[mAge+3,t,i] <- rpois(1,F[mAge,t])
      } # if
    } # t
  } # i
  
  # 4. Simulate the fates of individuals born during the study       
  # - determine the number of nestlings
  # - determine their fate over time
  # - determine their reproduction
  nestl <- numeric()
  nestl[1] <- 0
  for (t in 1:(T+1)){
    # 4.a: Enumerate the number of nestlings
    nestl[t+1] <- sum(ind[mAge+3,t,], na.rm = TRUE)
    ind[mAge+1,t,(max(Nimindex)+max(cumsum(nestl[1:t]))+1):(max(Nimindex)+max(cumsum(nestl[1:(t+1)])))] <- 1
    if (t==(T+1)) break
    
    # 4.b: Model survival of these individuals
    for (i in (max(Nimindex)+max(cumsum(nestl[1:t]))+1):(max(Nimindex)+max(cumsum(nestl[1:(t+1)])))){
      z <- numeric()
      for (d in t:T){
        z[d-t+1] <- rbinom(1,1,PHI[d-t+1,d])
      } # d
      Z <- max(sum(cumprod(z)))
      if (Z==0){
        ind[mAge+4,t+1,i] <- 1
        next
      } # if
      for (u in 1:Z){
        if (u < mAge){
          ind[u,u+t,i] <- 1
        } # if
        else {
          ind[mAge,u+t,i] <- 1
        } # else                
      }  # u
      # Record year of death (if any)
      if (sum(z)==T-t+1) next
      else {
        D <- min(which(z==0))
        ind[mAge+4,D+t,i] <- 1
      } # else          
    } # i
    
    # 4.c: Model reproduction of the surviving individuals
    for (i in (max(Nimindex)+max(cumsum(nestl[1:t]))+1):(max(Nimindex)+max(cumsum(nestl[1:(t+1)])))){
      for(d in t:T+1){       
        g <- which(!is.na(ind[c(1:mAge),d,i]))
        if (length(g)==0) next
        if (g !=8){
          ind[mAge+3,d,i] <- rpois(1,F[g,d])
        } # if
      } # d
    } # i
  } # t
  
  # 5. Enumerate the total number of animals
  Ntotal <- sum(Ni) + sum(ind[mAge+1,1:(T+1),], na.rm = TRUE) + sum(Im)
  # Remove empty cells and reorder the array such that it starts with the Juv
  IND <- ind[,,1:Ntotal]
  IND[1,,] <- ind[mAge+1,,1:Ntotal]
  for (a in 1:mAge){
    IND[a+1,,] <- ind[a,,1:Ntotal]
  } # a
  rnames <- numeric()
  for (a in 1:mAge){
    rnames[a] <- paste(a,"-Year", sep="")
  } # a
  rnames <- c("Juv", rnames, "Im", "Rep", "Dead")
  rownames(IND) <- rnames
  
  # Summary statistics: Number of individuals in each class and year, plus immigration rate
  Nu <- matrix(NA, ncol = T+1, nrow = mAge + 4)
  for (t in 1:(T+1)){
    for (a in 1:(mAge+1)){
      Nu[a,t] <- sum(IND[a,t,], na.rm = TRUE)
    } # a
    Nu[mAge+2,t] <- sum(IND[mAge+2,t,], na.rm = TRUE)
    Nu[mAge+3,t] <- sum(Nu[2:(mAge+1),t]) + sum(IND[mAge+2,t,], na.rm = TRUE)
  } # t
  for (t in 2:(T+1)){
    Nu[mAge+4,t] <- sum(IND[mAge+2,t,], na.rm = TRUE) / Nu[mAge+3,t-1]
  } # t
  rnames <- numeric()
  for (a in 1:mAge){
    rnames[a] <- paste(a,"-Year", sep="")
  } # a
  rnames <- c("Juv", rnames, "Im", "Total", "Imm rate")
  rownames(Nu) <- rnames
  
  # 6. Output
  return(list(IND = IND, Nu = Nu))
}





#####################################################################################################
#
# Function to create population survey data with a binomial sampling process
#
# It is assumed that all individuals have the same probability to be counted
#
# Input variables
# - Nu: Annual number of individuals at risk of detection (usually population size)
# - psur: vector with the annual detection probabilities of the individuals at risk of counting
#
# Last-up date: 9.6.2016, M.Schaub
#
#####################################################################################################

create.survey.bin <- function(Nu, psur, seed = NA){
  
  if (is.na(seed) == TRUE) seed = runif(1, 0, 100)
  # set.seed(seed)
  T <- length(Nu)
  SUR <- numeric()
  for (t in 1:T){
    SUR[t] <- rbinom(1, Nu[t], psur[t])
  } # t
  return(SUR)
}



#####################################################################################################
#
# Function to create population survey data with a Normal sampling process
#
# It is assumed that individuals may be double counted and missed at the same rate
#
# Input variables
# - Nu: Annual number of individuals at risk of detection (usually population size)
# - sigma: vector with the annual observation error
#
# Last-up date: 9.6.2016, M.Schaub
#
#####################################################################################################

create.survey.norm <- function(Nu, sigma, seed = NA){
  
  if (is.na(seed) == TRUE) seed = runif(1, 0, 100)
  # set.seed(seed)
  T <- length(Nu)
  SUR <- numeric()
  for (t in 1:T){
    SUR[t] <- rnorm(1, Nu[t], sigma[t])
  } # t
  return(SUR)
}



#####################################################################################################
#
# Function to create data on reproductive success
#
# Input variables
# - ind: array with the population
# - prep: vector with the annual detection probabilities of broods
#
# Output variables
# - rep.ind: matrix with the individual reproductive output. The three columns give the output, the year of the brood and the age of the mother.
#  - rep.agg: matrix with the same data, but aggregated. The two columns give the year-specific total number of newborn and the year-specific number of surveyed broods
#
# Last-up date: 14.3.2016, M.Schaub
#
#####################################################################################################

create.reproduction <- function(ind, prep, seed = NA){
  
  if (is.na(seed) == TRUE) seed = runif(1, 0, 100)
  # set.seed(seed)
  r <- dim(ind)[1] - 1
  maxAge <- dim(ind)[1] - 4
  T <- dim(ind)[2]
  rep <- year <- age <- numeric()
  le <- 0
  for (t in 1:T){
    z <- which(!is.na(ind[r,t,]))
    for (i in 1:length(z)){
      j <- le + i
      h <- rbinom(1, 1, prep[t])
      if (h==1){
        rep[j] <- ind[r,t,z[i]]
        year[j] <- t
        age[j] <- which(!is.na(ind[2:(maxAge+2),t,z[i]]))
      } # if
      else {
        rep[j] <- NA
        year[j] <- NA
        age[j] <- NA
      } # else  
    } # i
    le <- length(rep)
  } # t
  age[age==(maxAge+1)] <- maxAge   # re-adjust age of immigrants to maximal age
  k <- which(!is.na(rep))
  rep.ind <- cbind(rep[k], year[k], age[k])
  colnames(rep.ind) <- c("Reproduction", "Year", "Age of mother")
  rep.agg <- matrix(NA, nrow = T, ncol = 2)
  for (t in 1:T){
    rep.agg[t,1] <- sum(rep.ind[rep.ind[,2]==t,1])
    rep.agg[t,2] <- length(rep.ind[rep.ind[,2]==t,1])
  }
  colnames(rep.agg) <- c("Juveniles", "Surveyed broods")
  
  return(list(rep.ind = rep.ind, rep.agg = rep.agg))
}





#####################################################################################################
#
# Function to generate capture histories from the population (ind) 
#
# Input variables
# - ind: array with the population
# - c: matrix with age- and time-specific capture probabilities (probability of first capture)
# - p: matrix with age- and time-specific REcapture probabilities
# - maxAge: maximal number of age classes that can be identified when the individuals are captured for the first time
#
# Output
# - ch: matrix with the capture histories
# - age: vector with the age class at first capture for each individual
#
# Last up-date: 11.3.2016, M.Schaub
#
#####################################################################################################

create.capturehistory <- function(ind, c, p, maxAge = 2, seed = NA){
  
  if (is.na(seed) == TRUE) seed = runif(1, 0, 100)
  # set.seed(seed)
  T <- dim(ind)[2]
  nind <- dim(ind)[3]
  nstage <- dim(ind)[1]
  aclasses <- nstage-3
  age <- first <- last <- numeric()
  
  for (i in 1:nind){
    g <- which(!is.na(ind[1:(aclasses+1),,i]), arr.ind = TRUE)
    age[i] <- g[1,1]
    first[i] <- g[1,2]
    h <- which(ind[1:(aclasses+1),,i]==1, arr.ind = TRUE)
    last[i] <- max(h[,2])
  } # i
  
  ch.true <- ch <- matrix(0, ncol = T, nrow = nind)
  for (i in 1:nind){
    ch.true[i,first[i]:last[i]] <- 1
  } # i
  # Recode age
  age[age > maxAge] <- maxAge
  
  # Sampling
  # Expand c and p (to higher age classes)
  C <- matrix(0, ncol = T, nrow = max(c(maxAge, nrow(c))) + T)
  C[1:nrow(c),] <- c
  u <- max(c(maxAge, nrow(c))) + T - nrow(c)
  if (u > 0){
    for (j in 1:u){
      C[nrow(c)+j,] <- c[nrow(c),]
    } # j
  } # if
  
  P <- matrix(0, ncol = T-1, nrow = max(c(maxAge, nrow(p))) + T)
  P[1:nrow(p),] <- p
  u <- max(c(maxAge, nrow(p))) + T - nrow(p)
  if (u > 0){
    for (j in 1:u){
      P[nrow(p)+j,] <- p[nrow(p),]
    } # j
  } # if
  
  for (i in 1:nind){
    # First capture
    ch[i,first[i]] <- rbinom(1, 1, C[age[i],first[i]])
    if (first[i]==last[i]) next
    # Recapture (conditional on first capture)
    for (t in (first[i]+1):last[i]){
      ch[i,t] <- rbinom(1, 1, P[(age[i]+t-first[i]),t-1]) * ch[i,first[i]]
    } # t
  } # i
  
  # Remove individuals that have never been captured/marked
  incl <- which(rowSums(ch)>=1)
  ch <- ch[incl,]
  age <- age[incl]
  return(list(ch = ch, age = age))  
}






#####################################################################################################
#
# Function to create age-dependent m-arrays 
#
# Input variables
# - ch: matrix with capture histories. Note, this is a single file including all age classes
# - age: vector with the age for each individual at first capture
# - mAge: maximal number of age classes for which m-arrays are constructed. Input is optional and only required if the age matrix has fewer age classes as we want to separate (e.g. CH contains only individuals marked as juveniles, and we want 2 age classes)
#
# Output
# - marr: 3-d array with the m-array. The third dimension is the age class. The last column of each m-array is the number of released individuals that were never recaptured. Thus, the total number of released individuals per occasion is the row sum of each m-array.
#
# Last up-date: 14.3.2016, M.Schaub
#
#####################################################################################################

marray.age <- function(ch, age, mAge = 1){
  
  # 1. Helper functions
  # 1.1. Function to create a m-array based on capture-histories (ch)
  marray <- function(ch){
    nind <- nrow(ch)
    n.occasions <- ncol(ch)
    m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
    # Calculate the number of released individuals at each time period
    m.array[,1] <- colSums(ch)
    for (i in 1:nind){
      pos <- which(ch[i,]==1)
      g <- length(pos)
      if (g==1) next
      for (z in 1:(g-1)){
        m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
      } # z
    } # i
    # Calculate the number of individuals never recaptured
    for (t in 1:n.occasions){
      m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
    } # t
    out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
    return(out)
  }
  
  # 1.2. Function to remove histories without any capture from a capture-recapture matrix
  clean.ch <- function(ch){
    incl <- which(rowSums(ch)>=1)
    ch <- ch[incl,]
    return(ch)
  }
  
  # 1.3. Function to remove the first capture in a capture-recapture matrix
  rm.first <- function(ch) {
    get.first <- function(x) min(which(x==1))
    first <- apply(ch, 1, get.first)
    for (i in 1:nrow(ch)){
      ch[i,first[i]] <- 0
    }
    return(ch)
  }
  
  # 1.4. Function to calculate the occasion of first capture
  get.first <- function(x) min(which(x==1))
  
  
  # 2. Calculations   
  if (is.matrix(ch)==FALSE) ch <- matrix(ch, nrow = 1)   
  maxAge <- max(c(max(age), mAge))
  nind <- nrow(ch)
  n.occasions <- ncol(ch)
  
  first <- apply(ch, 1, get.first)
  age.matrix <- matrix(0, ncol = n.occasions, nrow = nind)
  for (i in 1:nind){
    age.matrix[i,first[i]:n.occasions] <- 1:(n.occasions-first[i]+1)+(age[i]-1)
  }
  age.matrix[age.matrix > maxAge] <- maxAge
  
  # Recode capture history
  ch.rec <- ch
  for (i in 1:nind){
    h <- which(ch.rec[i,]==1)
    for (j in 1:length(h)){
      ch.rec[i,h[j]] <- j
    } # j
  } # i
  ch.rec[ch.rec > maxAge] <- maxAge
  
  ch.split <- array(0, dim = c(nrow(ch), ncol(ch), maxAge))
  for (a in 1:maxAge){
    for (i in 1:nind){
      j <- which(ch.rec[i,]==a | ch.rec[i,]==(a+1))
      if (length(j)==0) next
      ch.split[i,j[1:2],age.matrix[i,j[1]]] <- 1
      if (length(j)>1){
        ch.split[i,j[2:length(j)],age.matrix[i,j[2]]] <- 1
      }
    } # i
  } # a
  
  marr <- array(0, dim = c(n.occasions-1, n.occasions, maxAge))
  for (a in 1:(maxAge-1)){
    for (i in 1:nind){
      u <- which(ch.split[i,,a]==1)
      if (length(u)==0) next
      if (u[1]==n.occasions) next
      if (length(u)==1) marr[u,n.occasions,a] <- marr[u,n.occasions,a] + 1
      if (length(u)==2) marr[u[1],u[2]-1,a] <- marr[u[1],u[2]-1,a] + 1
    } # i
  } # a
  a <- maxAge
  
  if (is.matrix(ch.split[,,a])==FALSE){ 
    ch.split1 <- matrix(ch.split[,,a], nrow = 1)
    marr[,,a] <- marray(ch.split1)
  } # if
  else marr[,,a] <- marray(ch.split[,,a])      
  return(marr)
}



#####################################################################################################
#
# Function to create an m-array for a single age class 
#
# Input variables
# - ch: matrix with capture histories
#
# Output
# - out: m-array. The last column of each m-array is the number of released individuals that were never recaptured. Thus, the total number of released individuals per occasion is the row sum of each m-array
#
# Last up-date: 14.3.2016, M.Schaub
#
#####################################################################################################

marray <- function(ch){
  
  nind <- nrow(ch)
  n.occasions <- ncol(ch)
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  m.array[,1] <- colSums(ch)
  for (i in 1:nind){
    pos <- which(ch[i,]==1)
    g <- length(pos)
    if (g==1) next
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } # z
  } # i
  # Calculate the number of individuals never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  } # t
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}




#####################################################################################################
#
# Function to create a summary statistics from a subset of a population
#
# This function is only needed, if the sampling(counts) is conducted on a subset of the complete population only (i.e. to have purely independent data)
#
# Input variables
# - ind: array with the detailed population data
# - incl: vector with the number of the individuals that are included for the summary statsitics
#
# Last-up date: 15.3.2016, M.Schaub
#
#####################################################################################################

summary.pop <- function(ind, incl){
  
  T <- dim(ind)[2]
  mAge <- dim(ind)[1]-4
  
  # Summary statistics: Number of individuals in each class and year, plus immigration rate
  Nu <- matrix(0, ncol = T, nrow = mAge + 4)
  ind.sel <- ind[,,incl]
  for (t in 1:T){
    for (a in 1:(mAge+1)){
      Nu[a,t] <- sum(ind[a,t,], na.rm = TRUE)
    }
  }
  rnames <- numeric()
  for (a in 1:mAge){
    rnames[a] <- paste(a,"-Year", sep="")
  }
  rnames <- c("Juv", rnames, "Im", "Total", "Imm rate")
  rownames(Nu) <- rnames
  return(Nu)
}





#####################################################################################################
#
# Function to remove histories without any capture from a capture-recapture matrix
#
# Input variables
# - ch: matrix with capture histories
#
# Last up-date: 14.3.2016, M.Schaub
#
#####################################################################################################

clean.ch <- function(ch){
  incl <- which(rowSums(ch)>=1)
  ch <- ch[incl,]
  return(ch)
}




#####################################################################################################
#
# Function to remove the first capture in a capture-recapture matrix
#
# Input variables
# - ch: matrix with capture histories
#
# Last up-date: 14.3.2016, M.Schaub
#
#####################################################################################################

rm.first <- function(ch){
  get.first <- function(x) min(which(x==1))
  first <- apply(ch, 1, get.first)
  for (i in 1:nrow(ch)){
    ch[i,first[i]] <- 0
  }
  return(ch)
}




#####################################################################################################
#
# Function to compute the occasion of first capture
#
# Written: 2011, BPA
#
#####################################################################################################

get.first <- function(x) min(which(x==1))




#####################################################################################################
#
# Function that calculates probablities that are used to create a discrete uniform prior with the categorical distribution in BUGS
#
# Input variables
# - A, B: range of the discrete uniform prior 
#
# Last up-date: September 2014, Michael Schaub
#
#####################################################################################################

disc.unif <- function(A, B){
  pprob <- c(rep(0, A-1), rep(1/(B-A+1), (B-A+1)))
  return(pprob)
}




#####################################################################################################
#
# Function to translate mean and sd of survival into the shape parameters of a beta distribution
#
# Input variables
# - x: mean of survival
# - sd.x: sd of survival
#
# Last up-date: 14.3.2016, M.Schaub
#
#####################################################################################################

beta.params <- function(x_bar, sd.x){
  u <- x_bar*(1-x_bar)/sd.x^2-1
  alpha <- x_bar*u
  beta <- (1-x_bar)*u
  return(list(alpha = alpha, beta = beta))
}




































############################################################################################
# Triecke edit to Schaub's create.capturehistory function
############################################################################################
create.ch.marker.loss <- function(ind, c, p, maxAge = 2, seed = NA, retention.rate = 0){
  
  if (is.na(seed) == TRUE) seed = runif(1, 0, 100)
  # set.seed(seed)
  T <- dim(ind)[2]
  nind <- dim(ind)[3]
  nstage <- dim(ind)[1]
  aclasses <- nstage-3
  age <- first <- last <- numeric()
  
  for (i in 1:nind){
    g <- which(!is.na(ind[1:(aclasses+1),,i]), arr.ind = TRUE)
    age[i] <- g[1,1]
    first[i] <- g[1,2]
    h <- which(ind[1:(aclasses+1),,i]==1, arr.ind = TRUE)
    last[i] <- max(h[,2])
  } # i
  
  ch.true <- ch <- matrix(0, ncol = T, nrow = nind)
  for (i in 1:nind){
    ch.true[i,first[i]:last[i]] <- 1
  } # i
  # Recode age
  age[age > maxAge] <- maxAge
  
  # Sampling
  # Expand c and p (to higher age classes)
  C <- matrix(0, ncol = T, nrow = max(c(maxAge, nrow(c))) + T)
  C[1:nrow(c),] <- c
  u <- max(c(maxAge, nrow(c))) + T - nrow(c)
  if (u > 0){
    for (j in 1:u){
      C[nrow(c)+j,] <- c[nrow(c),]
    } # j
  } # if
  
  P <- matrix(0, ncol = T-1, nrow = max(c(maxAge, nrow(p))) + T)
  P[1:nrow(p),] <- p
  u <- max(c(maxAge, nrow(p))) + T - nrow(p)
  if (u > 0){
    for (j in 1:u){
      P[nrow(p)+j,] <- p[nrow(p),]
    } # j
  } # if
  
  
  ############################################################################################
  ############################################################################################
  ############################################################################################
  # SIMULATE MARKER SURVIVAL
  ############################################################################################
  ############################################################################################
  ############################################################################################  
  marker <- matrix(0, ncol = T, nrow = nind)
  for (i in 1:nind){
    # First capture
    ch[i,first[i]] <- rbinom(1, 1, C[age[i],first[i]])
    marker[i,first[i]] <- 1
    if (first[i]==last[i]) next
    # Recapture (conditional on first capture)
    for (t in (first[i]+1):last[i]){
      marker[i,t] <- rbinom(1,marker[i,t-1],retention.rate)
      ch[i,t] <- rbinom(1, 1*marker[i,t], P[(age[i]+t-first[i]),t-1]) * ch[i,first[i]]
    } # t
  } # i
  
  # Remove individuals that have never been captured/marked
  incl <- which(rowSums(ch)>=1)
  ch <- ch[incl,]
  age <- age[incl]
  return(list(ch = ch, age = age))  
}


######################################################################################################
######################################################################################################
# end sim
######################################################################################################
######################################################################################################





















######################################################################################################
######################################################################################################
# Riecke 6 October 2016
######################################################################################################
######################################################################################################

######################################################################################################
######################################################################################################
# Integrated population model with immigration
######################################################################################################
######################################################################################################
# Specify the model in BUGS language
cat(file = "ipm_omega.jags", "
    model { 
    
    ############################################################################
    # Survival priors and constraints
    ############################################################################
    mean.sj ~ dunif(0, 1)
    mean.sa ~ dunif(0, 1)
    mean.p ~ dunif(0, 1)
    mean.f ~ dunif(0, 5)
    
    for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
    }
    


    sigma.obs ~ dunif(0.01, 250)
    tau.obs <- pow(sigma.obs, -2)
    
    
    ############################################################################
    # System Process
    ############################################################################
    
    # Model for the initial population size: discrete uniform priors
    mean1[1] ~ dunif(20, 1250)
    adult.init ~ dunif(20, 1250)
    N1[1] ~ dpois(mean1[1])
    Nad[1] ~ dpois(adult.init)
    Ntot[1] <- N1[1] + Nad[1]

    
    # Process model over time
    imm ~ dunif(-100,500)
    for (t in 1:(n.occasions-1)){
    immigrants[t] <- round(imm)
    mean1[t+1] <- mean.f * sj[t] * Ntot[t]
    N1[t+1] ~ dpois(mean1[t+1])
    Nad[t+1] ~ dbin(sa[t], Ntot[t])
    Ntot[t+1] <- N1[t+1] + Nad[t+1] + immigrants[t]
    # lambda[t] <- (Ntot[t+1])/(Ntot[t])
    }

    # mean.lambda <- mean(lambda[])  


    ############################################################################
    # State-space model for count data    
    # Observation model
    ############################################################################
    for (t in 1:n.occasions){
      count[t,1] ~ dnorm(Ntot[t], tau.obs)
      count[t,2] ~ dnorm(Ntot[t], tau.obs)
    }


    # Poisson regression model for productivity data
    for (i in 1:n.J){
    J[i] ~ dpois(mean.f)
    }
    
    # Capture-recapture model (multinomial likelihood)
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    }
    # Define the cell probabilities of the m-arrays
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
    pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    pr.a[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
    } #t
    

    }
    ")





######################################################################################################
######################################################################################################
# Integrated population model WITHOUT immigration
######################################################################################################
######################################################################################################
cat(file = "ipm1.jags", "
    model { 
    
    ############################################################################
    # Survival priors
    ############################################################################
    # Priors and constraints
    mean.sj ~ dunif(0, 1)
    mean.sa ~ dunif(0, 1)
    mean.p ~ dunif(0, 1)
    mean.f ~ dunif(0, 5)
    
    for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
    }
    

    ############################################################################
    # count priors
    ############################################################################
    sigma.obs ~ dunif(0, 250)
    tau.obs <- pow(sigma.obs, -2)
 

    ############################################################################
    # System Process
    ############################################################################
    # Model for the initial population size: discrete uniform priors
    mean1[1] ~ dunif(20, 1250)
    adult.init ~ dunif(20,1250)
    N1[1] ~ dpois(mean1[1])
    Nad[1] ~ dpois(adult.init)

    for (t in 1:(n.occasions)){
      Ntot[t] <- N1[t] + Nad[t]
    }

    # Process model over time
    for (t in 1:(n.occasions-1)){
    mean1[t+1] <- mean.f * sj[t] * Ntot[t]
    N1[t+1] ~ dpois(mean1[t+1])
    Nad[t+1] ~ dbin(sa[t], Ntot[t])
    # lambda[t] <- (Ntot[t+1])/(Ntot[t])
    }

    # mean.lambda <- mean(lambda[])


    ############################################################################
    # State-space model for count data    
    # Observation model
    ############################################################################
    for (t in 1:n.occasions){
      count[t,1] ~ dnorm(Ntot[t], tau.obs)
      count[t,2] ~ dnorm(Ntot[t], tau.obs)
    }
    
    ############################################################################
    # Poisson regression model for productivity data
    ############################################################################
    for (i in 1:n.J){
    J[i] ~ dpois(mean.f)
    }
    
    ############################################################################
    # Capture-recapture model (multinomial likelihood)
    # Define the multinomial likelihood
    ############################################################################
    for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    }
    # Define the cell probabilities of the m-arrays
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
    pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    pr.a[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
    } #t
    ############################################################################
    ############################################################################


    }

    ")



######################################################################################################
######################################################################################################
# Capture mark recapture model
######################################################################################################
######################################################################################################
# Specify the model in BUGS language
cat(file = "cmr.jags", "
    model { 
    
    ############################################################################
    # Survival priors and constraints
    ############################################################################
    mean.sj ~ dunif(0.1, 0.5)
    mean.sa ~ dunif(0.3, 0.7)
    mean.p ~ dunif(0.4, 0.8)

    
    for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
    }
    
    
    
    # Capture-recapture model (multinomial likelihood)
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
    }
    # Define the cell probabilities of the m-arrays
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
    pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    pr.a[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
    } #t
    
    }
    ")






######################################################################################################
######################################################################################################
# State-space-Model
######################################################################################################
######################################################################################################
cat(file = "ssm.jags", "
    model { 
    
    ############################################################################
    # count priors
    ############################################################################
    sigma.obs ~ dunif(0, 1500)
    tau.obs <- pow(sigma.obs, -2)
    
    sigma.proc ~ dunif(0, 10)
    tau.proc <- pow(sigma.proc, -2)    
    ############################################################################
    # System Process
    ############################################################################
    # Model for the initial population size: discrete uniform priors
    mean.lambda ~ dunif(0,5)
    Ntot[1] ~ dunif(0,10000)    

    for (t in 1:(n.occasions-1)){
      Ntot[t+1] <- Ntot[t] * lambda[t]
      lambda[t] ~ dnorm(mean.lambda, tau.proc) T(0,)
    }
    
    ############################################################################
    # State-space model for count data    
    # Observation model
    ############################################################################
    for (t in 1:n.occasions){
      count[t,1] ~ dnorm(Ntot[t], tau.obs)
      count[t,2] ~ dnorm(Ntot[t], tau.obs)
    }
    
    }
    
    ")




######################################################################################################
######################################################################################################
# Integrated population model WITHOUT immigration
######################################################################################################
######################################################################################################
cat(file = "fec.jags", "
    model { 
    
    ############################################################################
    # Survival priors
    ############################################################################
    # Priors and constraints
    mean.f ~ dunif(0, 5)
    
    ############################################################################
    # Poisson regression model for productivity data
    ############################################################################
    for (i in 1:n.J){
    J[i] ~ dpois(mean.f)
    }
    

    }
    
    ")




















##############################################################################
# population parameters and observation parameters
#############################################################################
# Age specific survival probabilities (juv, adult)
sj <- 0.3
sa <- 0.5

# Capture and recapture probabilities
cjuv <- 0.3           # initial capture probability of juveniles
cad <- 0.3            # initial capture probability of adults
prec <- 0.6           # recapture probability 

# Observation error for the population survey
sigma <- 100

# Fecundity rate (females)
f1 <- 1.6           # productivity of 1 year old females
f2 <- f1             # productivity of females older than one year

# Initial population size per age class
Ni <- c(1000, 1000)

# Number of years
T <- 10

# Probability to find a brood whose reproductive ouput is recorded
pprod <- 0.3


# install and load jagsUI
# install.packages('jagsUI')
require(jagsUI)

# Initial values
inits <- function(){list(mean.sj = 0.3, mean.sa = 0.5)}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "sigma.obs", "imm", "Ntot","change")
phi.parameters <- c("mean.sj", "mean.sa", "mean.p")
ssm.parameters <- c("sigma.obs")
fec.parameters <- c("mean.f")

# MCMC settings
ni <- 50000; nt <- 5; nb <- 25000; nc <- 2













##############################################################################
# simulate data and store results
#############################################################################
iter<-100
retention <- seq(0.9, 1, length.out = 2)

#############################################################################
# Create arrays to hold the results
#############################################################################
sim.results.omega <- array(NA, dim = c(12, length(retention), iter, 3))
sim.results <- array(NA, dim = c(10, length(retention), iter, 3))
phi.results <- array(NA, dim = c(6, length(retention), iter, 3))
ssm.results <- array(NA, dim = c(2, length(retention), iter, 3))
fec.results <- array(NA, dim = c(2, length(retention), iter, 3))

for (i in 1:length(retention)){
  cat(i," "); flush.console() # add counter
  for (j in 1:iter){
    cat(j," "); flush.console() # add counter
    ###############################################################################################
    # simulate a new dataset
    ###############################################################################################
    # Create the population survey data
    ind <- create.population(phi = matrix(c(rep(sj, T-1), rep(sa, T-1)), ncol = T-1, byrow = TRUE), 
                             f = matrix(c(rep(f1, T), rep(f2, T)), ncol = T, byrow = TRUE), Im = rep(0, T), 
                             Ni = Ni, seed = NA)
    
    c1 <- round(create.survey.norm(ind$Nu["Total",], rep(sigma,T), seed = NA))
    c2 <- round(create.survey.norm(ind$Nu["Total",], rep(sigma,T), seed = NA))
    
    count <- cbind(c1,c2)
    # Create the capture histories and the corresponding m-arrays
    ch <- create.ch.marker.loss(ind$IND, c = matrix(c(rep(cjuv, T), rep(cad, T)), nrow = 2, byrow = TRUE), 
                                p = matrix(c(rep(prec, T-1), rep(prec, T-1)), nrow = 2, byrow = TRUE), 
                                seed = NA, retention.rate = retention[i])
    marray <- marray.age(ch$ch, ch$age)
    
    # Create productivity data
    P <- create.reproduction(ind$IND, rep(pprod, T), seed = NA)
    
    ###############################################################################################
    # provide ipm data
    ###############################################################################################
    # Bundle data
    bugs.data <- list(marr.j = marray[,,1], marr.a = marray[,,2], n.occasions = T, rel.j = rowSums(marray[,,1]), 
                      rel.a = rowSums(marray[,,2]), J = P$rep.ind[,1], n.J = nrow(P$rep.ind), count = count)
    

    ###############################################################################################
    # run both models and store results
    ###############################################################################################
    # Call JAGS from R (jagsUI)
    m.omega <- jags(bugs.data, inits, parameters, "ipm_omega.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    sim.results.omega[1,i,j,1] <- m.omega$q50$mean.sj
    sim.results.omega[2,i,j,1] <- m.omega$q50$mean.sa     
    sim.results.omega[3,i,j,1] <- m.omega$q50$mean.f 
    sim.results.omega[4,i,j,1] <- m.omega$q50$mean.p 
    sim.results.omega[5,i,j,1] <- m.omega$q50$sigma.obs
    sim.results.omega[6,i,j,1] <- m.omega$q50$imm
    sim.results.omega[7,i,j,1] <- m.omega$sd$mean.sj
    sim.results.omega[8,i,j,1] <- m.omega$sd$mean.sa     
    sim.results.omega[9,i,j,1] <- m.omega$sd$mean.f 
    sim.results.omega[10,i,j,1] <- m.omega$sd$mean.p 
    sim.results.omega[11,i,j,1] <- m.omega$sd$sigma.obs
    sim.results.omega[12,i,j,1] <- m.omega$sd$imm

    sim.results.omega[1,i,j,2] <- m.omega$q2.5$mean.sj
    sim.results.omega[2,i,j,2] <- m.omega$q2.5$mean.sa     
    sim.results.omega[3,i,j,2] <- m.omega$q2.5$mean.f 
    sim.results.omega[4,i,j,2] <- m.omega$q2.5$mean.p 
    sim.results.omega[5,i,j,2] <- m.omega$q2.5$sigma.obs
    sim.results.omega[6,i,j,2] <- m.omega$q2.5$imm

    sim.results.omega[1,i,j,3] <- m.omega$q97.5$mean.sj
    sim.results.omega[2,i,j,3] <- m.omega$q97.5$mean.sa     
    sim.results.omega[3,i,j,3] <- m.omega$q97.5$mean.f 
    sim.results.omega[4,i,j,3] <- m.omega$q97.5$mean.p 
    sim.results.omega[5,i,j,3] <- m.omega$q97.5$sigma.obs
    sim.results.omega[6,i,j,3] <- m.omega$q97.5$imm      
    
    # Call JAGS from R (jagsUI)
    m <- jags(bugs.data, inits, parameters, "ipm1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    sim.results[1,i,j,1] <- m$q50$mean.sj
    sim.results[2,i,j,1] <- m$q50$mean.sa     
    sim.results[3,i,j,1] <- m$q50$mean.f
    sim.results[4,i,j,1] <- m$q50$mean.p 
    sim.results[5,i,j,1] <- m$q50$sigma.obs
    sim.results[6,i,j,1] <- m$sd$mean.sj
    sim.results[7,i,j,1] <- m$sd$mean.sa     
    sim.results[8,i,j,1] <- m$sd$mean.f
    sim.results[9,i,j,1] <- m$sd$mean.p 
    sim.results[10,i,j,1] <- m$sd$sigma.obs
    # sim.results[6,i,j] <- m$mean$mean.lambda
    
    sim.results[1,i,j,2] <- m$q2.5$mean.sj
    sim.results[2,i,j,2] <- m$q2.5$mean.sa     
    sim.results[3,i,j,2] <- m$q2.5$mean.f 
    sim.results[4,i,j,2] <- m$q2.5$mean.p 
    sim.results[5,i,j,2] <- m$q2.5$sigma.obs

    
    sim.results[1,i,j,3] <- m$q97.5$mean.sj
    sim.results[2,i,j,3] <- m$q97.5$mean.sa     
    sim.results[3,i,j,3] <- m$q97.5$mean.f 
    sim.results[4,i,j,3] <- m$q97.5$mean.p 
    sim.results[5,i,j,3] <- m$q97.5$sigma.obs

    
    ###############################################################################################
    # CMR model
    ###############################################################################################
    phi.data <- list(marr.j = marray[,,1], marr.a = marray[,,2], n.occasions = T, rel.j = rowSums(marray[,,1]), rel.a = rowSums(marray[,,2]))
    cmr <- jags(phi.data, inits, phi.parameters, "cmr.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    phi.results[1,i,j,1] <- cmr$q50$mean.sj
    phi.results[2,i,j,1] <- cmr$q50$mean.sa     
    phi.results[3,i,j,1] <- cmr$q50$mean.p
    phi.results[4,i,j,1] <- cmr$sd$mean.sj
    phi.results[5,i,j,1] <- cmr$sd$mean.sa     
    phi.results[6,i,j,1] <- cmr$sd$mean.p  
    
    phi.results[1,i,j,2] <- cmr$q2.5$mean.sj
    phi.results[2,i,j,2] <- cmr$q2.5$mean.sa     
    phi.results[3,i,j,2] <- cmr$q2.5$mean.p
    
    phi.results[1,i,j,3] <- cmr$q97.5$mean.sj
    phi.results[2,i,j,3] <- cmr$q97.5$mean.sa     
    phi.results[3,i,j,3] <- cmr$q97.5$mean.p
    

    
    ###############################################################################################
    # SSM model
    ###############################################################################################
    ssm.data <- list(n.occasions = T, count = count)
    ssm <- jags(ssm.data, inits <- function(){list(sigma.obs = 100)}, ssm.parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    ssm.results[1,i,j,1] <- ssm$q50$sigma.obs
    ssm.results[2,i,j,1] <- ssm$sd$sigma.obs
    ssm.results[1,i,j,2] <- ssm$q2.5$sigma.obs
    ssm.results[1,i,j,3] <- ssm$q97.5$sigma.obs
    
    ###############################################################################################
    # FEC model
    ###############################################################################################
    fec.data <- list(J = P$rep.ind[,1], n.J = nrow(P$rep.ind))
    fec <- jags(fec.data, inits <- function(){list(mean.f = 1)}, fec.parameters, "fec.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    fec.results[1,i,j,1] <- fec$q50$mean.f
    fec.results[2,i,j,1] <- fec$sd$mean.f
    fec.results[1,i,j,2] <- fec$q2.5$mean.f
    fec.results[1,i,j,3] <- fec$q97.5$mean.f
    
  }
}

#
# sim.results.omega
# sim.results

save.image("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\Resubmission\\\\ME_final_q50_THREE.RData")


############################################################################################################
# PLOTS
############################################################################################################

pdf("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\new_TEX\\\\Figure2.pdf", width=10, height=10)
par(family = 'sans', mar = c(5.1,5.1,2.1,2.1), mfrow = c(1,1), oma = c(0,0,0,0))
boxplot(phi.results[2,1,], sim.results.omega[2,1,], sim.results[2,1,], 
        phi.results[2,2,], sim.results.omega[2,2,], sim.results[2,2,],
        phi.results[2,3,], sim.results.omega[2,3,], sim.results[2,3,],
        ylim = c(0.42,0.525), outline = FALSE, boxwex = 0.35, xlab = "Model",
        ylab = expression(phi[ad]), cex.lab = 2,
        names = c('CJS',expression(IPM[imm]),expression(IPM[no-imm]),
                  'CJS',expression(IPM[imm]),expression(IPM[no-imm]),
                  'CJS',expression(IPM[imm]),expression(IPM[no-imm])),
        at = c(1,1.5,2, 3,3.5,4, 5,5.5,6), xlim = c(0.75,6.25),
        col = c('royalblue1', 'yellow', 'grey',
                'royalblue1', 'yellow', 'grey',
                'royalblue1', 'yellow', 'grey'))
abline(v = 2.5, lwd = 1)
abline(v = 4.5, lwd = 1)
abline(0.5,0, lty = 2, lwd = 2)
text(1, 0.42, "0.10", cex = 2)
text(3, 0.42, "0.05", cex = 2)
text(5, 0.42, "0", cex = 2)
dev.off()


pdf("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\new_TEX\\\\Figure3.pdf", width=10, height=10)
par(family = 'sans', mar = c(5.1,5.1,2.1,2.1), mfrow = c(1,1), oma = c(0,0,0,0))
boxplot(phi.results[1,1,], sim.results.omega[1,1,], sim.results[1,1,], 
        phi.results[1,2,], sim.results.omega[1,2,], sim.results[1,2,],
        phi.results[1,3,], sim.results.omega[1,3,], sim.results[1,3,],
        ylim = c(0.24,0.325), outline = FALSE, boxwex = 0.35, xlab = "Model",
        ylab = expression(phi['juvenile']), cex.lab = 2, xlim = c(0.75,6.25),
        names = c('CJS',expression(IPM[imm]),expression(IPM[no-imm]),
                  'CJS',expression(IPM[imm]),expression(IPM[no-imm]),
                  'CJS',expression(IPM[imm]),expression(IPM[no-imm])),
        at = c(1,1.5,2, 3,3.5,4, 5,5.5,6),
        col = c('royalblue1', 'yellow', 'grey',
                'royalblue1', 'yellow', 'grey',
                'royalblue1', 'yellow', 'grey'))
abline(v = 2.5, lwd = 1)
abline(v = 4.5, lwd = 1)
abline(0.3,0, lty = 2, lwd = 2)
text(1, 0.24, "0.10", cex = 2)
text(3, 0.24, "0.05", cex = 2)
text(5, 0.24, "0", cex = 2)
dev.off()






pdf("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\new_TEX\\\\Figure4.pdf", width=10, height=10)
par(family = 'sans', mar = c(5.1,5.1,2.1,2.1), mfrow = c(1,1), oma = c(0,0,0,0))
boxplot(phi.results[3,1,], sim.results.omega[4,1,], sim.results[4,1,], 
        phi.results[3,2,], sim.results.omega[4,2,], sim.results[4,2,],
        phi.results[3,3,], sim.results.omega[4,3,], sim.results[4,3,],
        ylim = c(0.53,0.63), outline = FALSE, boxwex = 0.35, xlab = "Model",
        ylab = expression(italic(p)), cex.lab = 2, xlim = c(0.75,6.25),
        names = c('CJS',expression(IPM[imm]),expression(IPM[no-imm]),
                  'CJS',expression(IPM[imm]),expression(IPM[no-imm]),
                  'CJS',expression(IPM[imm]),expression(IPM[no-imm])),
        at = c(1,1.5,2, 3,3.5,4, 5,5.5,6),
        col = c('royalblue1', 'yellow', 'grey',
                'royalblue1', 'yellow', 'grey',
                'royalblue1', 'yellow', 'grey'))
abline(v = 2.5, lwd = 1)
abline(v = 4.5, lwd = 1)
abline(0.6,0, lty = 2, lwd = 2)
text(1, 0.53, "0.10", cex = 2)
text(3, 0.53, "0.05", cex = 2)
text(5, 0.53, "0", cex = 2)
dev.off()






pdf("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\new_TEX\\\\Figure5.pdf", width=10, height=10)
par(family = 'sans', mar = c(5.1,5.1,2.1,2.1), mfrow = c(1,1), oma = c(0,0,0,0))
boxplot(fec.results[1,1,], sim.results.omega[3,1,], sim.results[3,1,], 
        fec.results[1,2,], sim.results.omega[3,2,], sim.results[3,2,],
        fec.results[1,3,], sim.results.omega[3,3,], sim.results[3,3,],
        ylim = c(1.55,1.67), outline = FALSE, boxwex = 0.35, xlab = "Model",
        ylab = expression(italic(f)), cex.lab = 2, xlim = c(0.75,6.25),
        names = c('PR',expression(IPM[imm]),expression(IPM[no-imm]),
                  'PR',expression(IPM[imm]),expression(IPM[no-imm]),
                  'PR',expression(IPM[imm]),expression(IPM[no-imm])),
        at = c(1,1.5,2, 3,3.5,4, 5,5.5,6),
        col = c('royalblue1', 'yellow', 'grey',
                'royalblue1', 'yellow', 'grey',
                'royalblue1', 'yellow', 'grey'))
abline(v = 2.5, lwd = 1)
abline(v = 4.5, lwd = 1)
abline(1.6,0, lty = 2, lwd = 2)
text(1, 1.55, "0.10", cex = 2)
text(3, 1.55, "0.05", cex = 2)
text(5, 1.55, "0", cex = 2)
dev.off()




pdf("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\new_TEX\\\\Figure6.pdf", width=10, height=10)
par(family = 'sans', mar = c(5.1,5.1,2.1,2.1), mfrow = c(1,1), oma = c(0,0,0,0))
boxplot(ssm.results[1,1,], sim.results.omega[5,1,], sim.results[5,1,], 
        ssm.results[1,2,], sim.results.omega[5,2,], sim.results[5,2,],
        ssm.results[1,3,], sim.results.omega[5,3,], sim.results[5,3,],
        ylim = c(50,165), outline = FALSE, boxwex = 0.35, xlab = "Model",
        ylab = expression(sigma[y]), cex.lab = 2, xlim = c(0.75,6.25),
        names = c('SSM',expression(IPM[imm]),expression(IPM[no-imm]),
                  'SSM',expression(IPM[imm]),expression(IPM[no-imm]),
                  'SSM',expression(IPM[imm]),expression(IPM[no-imm])),
        at = c(1,1.5,2, 3,3.5,4, 5,5.5,6),
        col = c('royalblue1', 'yellow', 'grey',
                'royalblue1', 'yellow', 'grey',
                'royalblue1', 'yellow', 'grey'))
abline(v = 2.5, lwd = 1)
abline(v = 4.5, lwd = 1)
abline(100,0, lty = 2, lwd = 2)
text(1, 50, "0.10", cex = 2)
text(3, 50, "0.05", cex = 2)
text(5, 50, "0", cex = 2)
dev.off()



pdf("W:\\\\Ph_D\\\\Manuscripts\\\\IPM_bias\\\\MEE\\\\new_TEX\\\\Figure7.pdf", width=10, height=10)
par(family = 'sans', mar = c(5.1,5.1,2.1,2.1), mfrow = c(1,1), oma = c(0,0,0,0))
boxplot(sim.results.omega[6,1,],sim.results.omega[6,2,],sim.results.omega[6,3,],
        ylim = c(-50,250), ylab = expression(omega), cex.lab = 2, xlab = 'Marker Loss Rate',
        names = c('0.10','0.05','0'), outline = FALSE,
        col = c('yellow', 'yellow', 'yellow'))
abline(0, 0, lty = 2, lwd = 2)
dev.off()


























################################################################################
# calculate bias
################################################################################

################################################################################
# adult survival
################################################################################
############################################
# IPM w immigration
############################################
mean((sim.results.omega[2,3,] - 0.5)/0.5) # relative
mean((sim.results.omega[2,2,] - 0.5)/0.5) # relative
mean((sim.results.omega[2,1,] - 0.5)/0.5) # relative

mean(sim.results.omega[2,3,] - 0.5) # absolute
mean(sim.results.omega[2,2,] - 0.5) # absolute
mean(sim.results.omega[2,1,] - 0.5) # absolute

############################################
# IPM w/o immigration
############################################
mean((sim.results[2,3,] - 0.5)/0.5) # relative
mean((sim.results[2,2,] - 0.5)/0.5) # relative
mean((sim.results[2,1,] - 0.5)/0.5) # relative

mean(sim.results[2,3,] - 0.5) # relative
mean(sim.results[2,2,] - 0.5) # relative
mean(sim.results[2,1,] - 0.5) # relative


############################################
# cmr
############################################
mean((phi.results[2,3,] - 0.5)/0.5) # relative
mean((phi.results[2,2,] - 0.5)/0.5) # relative
mean((phi.results[2,1,] - 0.5)/0.5) # relative

mean(phi.results[2,3,] - 0.5) # relative
mean(phi.results[2,2,] - 0.5) # relative
mean(phi.results[2,1,] - 0.5) # relative




################################################################################
# juvenlie survival
################################################################################
############################################
# IPM w immigration
############################################
mean((sim.results.omega[1,3,] - 0.3)/0.3) # relative
mean((sim.results.omega[1,2,] - 0.3)/0.3) # relative
mean((sim.results.omega[1,1,] - 0.3)/0.3) # relative

mean(sim.results.omega[1,3,] - 0.3) # absolute
mean(sim.results.omega[1,2,] - 0.3) # absolute
mean(sim.results.omega[1,1,] - 0.3) # absolute

############################################
# IPM w/o immigration
############################################
mean((sim.results[1,3,] - 0.3)/0.3) # relative
mean((sim.results[1,2,] - 0.3)/0.3) # relative
mean((sim.results[1,1,] - 0.3)/0.3) # relative

mean(sim.results[1,3,] - 0.3) # relative
mean(sim.results[1,2,] - 0.3) # relative
mean(sim.results[1,1,] - 0.3) # relative


############################################
# cmr
############################################
mean((phi.results[1,3,] - 0.3)/0.3) # relative
mean((phi.results[1,2,] - 0.3)/0.3) # relative
mean((phi.results[1,1,] - 0.3)/0.3) # relative

mean(phi.results[1,3,] - 0.3) # relative
mean(phi.results[1,2,] - 0.3) # relative
mean(phi.results[1,1,] - 0.3) # relative




################################################################################
# detection
################################################################################
############################################
# IPM w immigration
############################################
mean((sim.results.omega[4,3,] - 0.6)/0.6) # relative
mean((sim.results.omega[4,2,] - 0.6)/0.6) # relative
mean((sim.results.omega[4,1,] - 0.6)/0.6) # relative

mean(sim.results.omega[4,3,] - 0.6) # absolute
mean(sim.results.omega[4,2,] - 0.6) # absolute
mean(sim.results.omega[4,1,] - 0.6) # absolute

############################################
# IPM w/o immigration
############################################
mean((sim.results[4,3,] - 0.6)/0.6) # relative
mean((sim.results[4,2,] - 0.6)/0.6) # relative
mean((sim.results[4,1,] - 0.6)/0.6) # relative

mean(sim.results[4,3,] - 0.6) # relative
mean(sim.results[4,2,] - 0.6) # relative
mean(sim.results[4,1,] - 0.6) # relative


############################################
# cmr
############################################
mean((phi.results[3,3,] - 0.6)/0.6) # relative
mean((phi.results[3,2,] - 0.6)/0.6) # relative
mean((phi.results[3,1,] - 0.6)/0.6) # relative

mean(phi.results[3,3,] - 0.6) # relative
mean(phi.results[3,2,] - 0.6) # relative
mean(phi.results[3,1,] - 0.6) # relative









################################################################################
# observation error
################################################################################
############################################
# IPM w immigration
############################################
mean((sim.results.omega[5,3,] - 100)/100) # relative
mean((sim.results.omega[5,2,] - 100)/100) # relative
mean((sim.results.omega[5,1,] - 100)/100) # relative

mean(sim.results.omega[5,3,] - 100) # absolute
mean(sim.results.omega[5,2,] - 100) # absolute
mean(sim.results.omega[5,1,] - 100) # absolute

############################################
# IPM w/o immigration
############################################
mean((sim.results[5,3,] - 100)/100) # relative
mean((sim.results[5,2,] - 100)/100) # relative
mean((sim.results[5,1,] - 100)/100) # relative

mean(sim.results[5,3,] - 100) # relative
mean(sim.results[5,2,] - 100) # relative
mean(sim.results[5,1,] - 100) # relative


############################################
# ssm
############################################
mean((ssm.results[1,3,] - 100)/100) # relative
mean((ssm.results[1,2,] - 100)/100) # relative
mean((ssm.results[1,1,] - 100)/100) # relative

mean(ssm.results[1,3,] - 100) # relative
mean(ssm.results[1,2,] - 100) # relative
mean(ssm.results[1,1,] - 100) # relative










################################################################################
# fecundity
################################################################################
############################################
# IPM w immigration
############################################
mean((sim.results.omega[3,3,] - 1.6)/1.6) # relative
mean((sim.results.omega[3,2,] - 1.6)/1.6) # relative
mean((sim.results.omega[3,1,] - 1.6)/1.6) # relative

mean(sim.results.omega[3,3,] - 1.6) # absolute
mean(sim.results.omega[3,2,] - 1.6) # absolute
mean(sim.results.omega[3,1,] - 1.6) # absolute

############################################
# IPM w/o immigration
############################################
mean((sim.results[3,3,] - 1.6)/1.6) # relative
mean((sim.results[3,2,] - 1.6)/1.6) # relative
mean((sim.results[3,1,] - 1.6)/1.6) # relative

mean(sim.results[3,3,] - 1.6) # relative
mean(sim.results[3,2,] - 1.6) # relative
mean(sim.results[3,1,] - 1.6) # relative


############################################
# fec
############################################
mean((fec.results[1,3,] - 1.6)/1.6) # relative
mean((fec.results[1,2,] - 1.6)/1.6) # relative
mean((fec.results[1,1,] - 1.6)/1.6) # relative

mean(fec.results[1,3,] - 1.6) # relative
mean(fec.results[1,2,] - 1.6) # relative
mean(fec.results[1,1,] - 1.6) # relative







################################################################################
# immigration
################################################################################
############################################
# IPM w immigration
############################################
mean(sim.results.omega[6,3,] - 0) # absolute
mean(sim.results.omega[6,2,] - 0) # absolute
mean(sim.results.omega[6,1,] - 0) # absolute

mean(sim.results.omega[6,3,] - 0)/1800 # absolute
mean(sim.results.omega[6,2,] - 0)/1800 # absolute
mean(sim.results.omega[6,1,] - 0)/1800 # absolute







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





























