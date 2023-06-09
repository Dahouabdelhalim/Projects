#Load required packages
library(ggplot2)
library(reshape)
library(gridExtra)
library(plyr)

#function to find the K value for females and males
lambda.0 = 4 
beta = 0.003
s.jf0 = 0.5
gamma= 0
alphaf = 0.01
s.f0 = 0.5
wff = 1
wfm = 1.4
s.jm0 = 0.6
wmf = 0.8
wmm = 1
s.m0 = 0.6
alpham = 0.009


func <- function(par){
  Nf <- exp(par[1])
  
  Nm <- exp(par[2])
  
  Xf <- ((lambda.0*exp(-beta*(Nf + Nm)))/2)*(s.jf0*exp(-gamma*(Nf + Nm)))
  
  Sf <- s.f0*exp(-alphaf*(wff*Nf + wfm*Nm))
  
  Xm <- (((Nf*lambda.0*exp(-beta*(Nf + Nm)))/Nm)/2)*(s.jm0*exp(-gamma*(Nf + Nm)))
  
  Sm <- s.m0*exp(-alpham*(wmf*Nf + wmm*Nm))
  
  lambda.f <- Xf + Sf
  lambda.m <- Xm + Sm
  
  sum( (lambda.f-1)^2 + (lambda.m-1)^2 ) 
}


vec <- c( log(50), log(70) ) #guess of parameters

exp(optim(vec, func)$par) #result (Nf, Nm)


####density dependence across time - males father newborns (not for genetic monogamy)####
#
gp.mod <- function(nsim = 10, timesteps = 20, Nf0 = 40, Nm0 = 40, A = 1, alphaf  = 0, alpham = 0, beta = 0, gamma = 0, eps.a = 0, eps.j = 0, eps.r = 0, wf = 1, wm = 1, wfm = 1, wmf =1, tau.f = 0, tau.m = 0, lambda.0 = 2, s.jf0 = 0.3, s.jm0 = 0.3, s.f0 = 0.6, s.m0 = 0.7, pa0 = 0.5) {
  q = 1-pa0
  u = vector(mode = "numeric")
  env.a = vector(mode = "numeric")
  env.j = vector(mode = "numeric")
  env.r = vector(mode = "numeric")
  #Vector for environmental fluctuations. Should be equal for all individuals, but not necessarily equal between processes 
  for(i in 1:timesteps) { 
    u[i] = rnorm(n = 1, 0,1)
    env.a[i] = eps.a*u[i] #For survival, the last term vanish because of the logit function
    env.j[i] = eps.j*u[i]
    env.r[i] = eps.r*u[i] - 0.5*(eps.r^2)
  } 
  #Define matrixes we need:
  Ntot = matrix(data = rep(c(Nf0+Nm0, rep(NA, timesteps-1)), nsim), nrow = timesteps) #Total number of individuals at each timestep
  Nf = matrix(data = rep(c(Nf0, rep(NA, timesteps-1)), nsim), nrow = timesteps) #number of females
  Nm = matrix(data = rep(c(Nm0, rep(NA, timesteps-1)), nsim), nrow = timesteps)
 varbf <- matrix(nrow = timesteps, ncol = nsim)
  meanbf <- matrix(nrow = timesteps, ncol = nsim)
  varbm <- matrix(nrow = timesteps, ncol = nsim)
 meanbm <- matrix(nrow = timesteps, ncol = nsim)
 Bm = matrix(data = rep(rep(NA, timesteps), nsim), nrow = timesteps)
Bf = matrix(data = rep(rep(NA, timesteps), nsim), nrow = timesteps)
  p.a <- matrix(data = rep(c(pa0, rep(NA, timesteps-1)), nsim), nrow = timesteps)
  p.af <- matrix(data = rep(c(pa0, rep(NA, timesteps-1)), nsim), nrow = timesteps)
  delta.p <- matrix(data = rep(rep(NA, timesteps), nsim), nrow = timesteps)
  p.am <- matrix(data = rep(c(pa0, rep(NA, timesteps-1)), nsim), nrow = timesteps)
  alleles2 <- 0 #this is where the resulting alleles for each j, i != 2 goes. Genotypes will be sampled from here
  alleles.juv <- 0 
  genefreqAA <- matrix(nrow = timesteps, ncol = nsim)
  genefreqaa <- matrix(nrow = timesteps, ncol = nsim)
  genefreqAa <- matrix(nrow = timesteps, ncol = nsim)
  for(j in 1:nsim){ #number of simulations
    genotypeind <- data.frame(time = numeric(),
                              id = character(),
                              a1 = character(),
                              a2 = character(),
                              n.off = numeric(),
                              surv = numeric(),
                              sex = character(),
                              z = numeric(),
                              pm = numeric())
    for(i in 2:timesteps) { #Number of timesteps
      b <- 0
      if(i == 2) {
        genotype <- matrix(nrow = Ntot[i-1, j], ncol = 8, dimnames = list(NULL, c("id", "a1", "a2", "n.off", "surv", "sex", "z", "pm")))
        genotype[,1] = paste(j, i-1, 1:Ntot[i-1, j])
        #alleles for parent generation.
        alleles2 <- c(replicate(round(2*q*Ntot[i-1, j]), "A"), replicate(2*pa0*Ntot[i-1, j], "a"))
        #print(alleles2)
        genotype[,2:3] = sample(alleles2, 2*Ntot[i-1, j])
        a <- paste(genotype[,2], genotype[,3])
        genefreqaa[i-1,j] = length(a[a == "a a"])/Ntot[i-1, j]
        genefreqAA[i-1,j] = length(a[a == "A A"])/Ntot[i-1, j]
        genefreqAa[i-1,j] = length(a[a %in% c("A a", "a A")])/Ntot[i-1, j]
        #sex need to be given
        sexc <- c(replicate(round(Nf0), "f"), replicate(Nm0, "m"))
        genotype[,6] = sample(sexc)
        #print(genotype)
        #females reproduce
        #poisson parameter: 
        lambda.pois = lambda.0*exp(-beta*(Ntot[i-1, j]/A)+env.r[i-1]) 
        #poisson-log-normal parameter:
        #print(lambda.pois)
        
        b = rpois(Nf[i-1, j], lambda = lambda.pois*exp(tau.f*rnorm(n = Nf[i-1, j], 0,1) - 0.5*(tau.f^2))) #nr of progeny produced by each female
        genotype[genotype[,6] == "f", 4] <- b 
        #print(genotype)
        meanbf[i-1,j] <- mean(b)
        varbf[i-1,j] <- var(b)
        ##males get their breeding potential
        genotype[genotype[,6] == "m", 7] <-rlnorm(Nm[i-1,j], meanlog = -0.5*tau.m^2, sdlog = tau.m)
        ##This gives a probability vector:
        genotype[, 8] <- as.numeric(genotype[,7])/sum(as.numeric(genotype[,7]), na.rm = T)
        #print(genotype)
        #print(sum(as.numeric(genotype[,8]), na.rm = T))
        ##males are given progeny according to their pm values with a sum that equals sum(b):
        genotype[genotype[,6] == "m", 4] <- rmultinom(1, sum(b), as.numeric(genotype[genotype[,6] == "m", 8]))
        varbm[i-1,j] <- var(as.numeric(genotype[genotype[,6] == "m", 4]))
        meanbm[i-1,j] <- mean(as.numeric(genotype[genotype[,6] == "m", 4]))
        Bm[i-1,j] <- length(genotype[genotype[,6] == "m" & genotype[,4] > 0, 4])
        Bf[i-1,j] <- length(genotype[genotype[,6] == "f" & genotype[,4] > 0, 4])
        #print(genotype)
        #print(Bm)
        #print(sum(b))
        #print(sum(as.numeric(genotype[genotype[,6] == "m", 4])))
        ##progeny recieve genotype 
        alleles.juvf <- as.vector(unlist(apply(genotype[genotype[,6] == "f", ], 1, function(dat){
          sample(x = dat[2:3], size = dat[4], replace = TRUE)
        }))) #female alleles
        alleles.juvm <- as.vector(unlist(apply(genotype[genotype[,6] == "m", ], 1, function(dat){
          sample(x = dat[2:3], size = dat[4], replace = TRUE) #male alleles
        }))) #male alleles
        #print(sum(b))
        #print(length(alleles.juvf))
      #print(length(alleles.juvm))
        #print(alleles.juvm)
        gtjuv <- matrix(nrow = sum(b), ncol = 8, dimnames = list(NULL, c("id", "a1", "a2", "n.off", "surv", "sex", "z", "pm"))) #here, b = number of offsprings
        gtjuv[,1] <- paste(j, i, seq(from = Ntot[i-1, j]+1, length.out = sum(b)))
        gtjuv[,2] <- sample(alleles.juvf, sum(b))
        gtjuv[,3] <- sample(alleles.juvm, sum(b))
        gtjuv[,4] <- rep(0, sum(b))
        gtjuv[,6] <- sample(c("f", "m"), sum(b),replace = T, c(0.5, 0.5)) #offspring get sex
        #print(gtjuv)
        
      ##adults survive to year i+1 or dies:
        #weigthed values for N. One for females, and one for males
        wNf = wf*Nf[i-1, j] + wfm*Nm[i-1, j] #females might be more affected by other females
      wNm = wmf*Nf[i-1, j] + wm*Nm[i-1, j]#males may be more affected by other males
        #females:
        genotype[genotype[,6] == "f",5] = rbinom(Nf[i-1, j], 1, prob = s.f0*exp(-alphaf*(wNf/A)+env.a[i-1])/(1-s.f0*exp(-alphaf*(wNf/A))+s.f0*exp(-alphaf*(wNf/A)+env.a[i-1]))) #For now, dp is different between sexes
        #males:
      genotype[genotype[,6] == "m",5] = rbinom(Nm[i-1, j], 1, prob = s.m0*exp(-alpham*(wNm/A)+env.a[i-1])/(1-s.m0*exp(-alpham*(wNm/A))+s.m0*exp(-alpham*(wNm/A)+env.a[i-1])))
      #print(genotype)
        #juveniles survive or die: 
      #females:
        gtjuv[gtjuv[,6] == "f",5] = rbinom(length(gtjuv[which(gtjuv[,6] == "f")]), 1, prob = s.jf0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])/(1-s.jf0*exp(-gamma*(Ntot[i-1, j]/A)) + s.jf0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])))
      #males:
      gtjuv[gtjuv[,6] == "m",5] = rbinom(length(gtjuv[which(gtjuv[,6] == "m")]), 1, prob = s.jm0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])/(1-s.jm0*exp(-gamma*(Ntot[i-1, j]/A)) + s.jm0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])))
        #print(gtjuv)
      #male juveniles recieve their z's:
      gtjuv[gtjuv[,6] == "m", 7] <-rlnorm(length(gtjuv[which(gtjuv[,6] == "m")]), meanlog = -0.5*tau.m^2, sdlog = tau.m)
      #print(gtjuv)
        #Gives the startpopulation present before reproduction: 
        genotypeind = rbind(genotypeind, data.frame(time = rep(i-1, Ntot[i-1, j]), genotype))
        
        #those that survive
        genotypeind = rbind(genotypeind, data.frame(time = rep(i, sum(as.numeric(genotype[,5]))),
                                                    subset(genotype, genotype[,5] == 1)),
                            data.frame(time = rep(i, sum(as.numeric(gtjuv[,5]))),
                                       subset(gtjuv, gtjuv[,5] == 1)))
        #print(genotypeind)
        #Popsize and allelfreq at i = 2 
        Ntot[i,j] = sum(as.numeric(genotype[,5])) + sum(as.numeric(gtjuv[,5]))
        Nf[i,j] = length(genotypeind$sex[genotypeind$sex == "f" & genotypeind$time == i])
      #print(Nf)
        Nm[i,j] = Ntot[i,j] - Nf[i,j]
        alive <- genotypeind[genotypeind$surv == "1" & genotypeind$time == i,]
        alivef <- genotypeind[genotypeind$sex == "f" & genotypeind$surv == "1" & genotypeind$time == i,]
        alivem <- genotypeind[genotypeind$sex == "m" & genotypeind$surv == "1" & genotypeind$time == i,]
        #print(alive)
        pool = with(alive, c(as.character(a1), as.character(a2))) #p.a is allelefreq right before reproduction
        #print(pool)
        poolf = with(alivef, c(as.character(a1), as.character(a2))) 
        poolm = with(alivem, c(as.character(a1), as.character(a2))) 
        p.a[i, j] = length(pool[pool == "a"])/length(pool)
        p.af[i, j] = length(poolf[poolf == "a"])/length(poolf)
        p.am[i, j] = length(poolm[poolm == "a"])/length(poolm)
        delta.p[i-1,j] = p.a[i,j] - p.a[i-1,j]
        a <- paste(alive$a1, alive$a2)
        genefreqaa[i,j] = length(a[a == "a a"])/Ntot[i, j]
        genefreqAA[i,j] = length(a[a == "A A"])/Ntot[i, j]
        genefreqAa[i,j] = length(a[a %in% c("A a", "a A")])/Ntot[i, j]
      
      }
      if(i > 2) {
        if(Nf[i-1, j] %in% c(0,NA) | Nm[i-1, j] %in% c(0, NA)) {break}
        if(!Nf[i-1, j] %in% c(0,NA) & !Nm[i-1, j] %in% c(0, NA)) {#
        genotype = as.matrix(genotypeind[genotypeind$time == i-1,])
        genotype[,1] = i
        #
        lambda.pois = lambda.0*exp(-beta*(Ntot[i-1,j]/A)+env.r[i-1])
        
        b = rpois(Nf[i-1, j], lambda = lambda.pois*exp(tau.f*rnorm(n = Nf[i-1, j], 0,1) - 0.5*(tau.f^2))) #
        
        genotype[genotype[,7] == "f", 5] <- b #
        
        #print(length(genotype[genotype[,7] == "f", 3]))
        meanbf[i-1,j] <- mean(b)
        varbf[i-1,j] <- var(b) 
       ##males get new prob-vector:
       genotype[, 9] <- as.numeric(genotype[,8])/sum(as.numeric(genotype[,8]), na.rm = T)
       #print(genotype)
       #print(sum(as.numeric(genotype[,8]), na.rm = T))
       ##males are given progeny according to their pm values with a sum that equals sum(b):
       genotype[genotype[,7] == "m", 5] <- rmultinom(1, sum(b), as.numeric(genotype[genotype[,7] == "m", 9]))
       varbm[i-1,j] <- var(as.numeric(genotype[genotype[,7] == "m", 5]))
       meanbm[i-1,j] <- mean(as.numeric(genotype[genotype[,7] == "m", 5]))
       Bm[i-1,j] <- length(genotype[genotype[,7] == "m" & genotype[,5] > 0, 5]) 
       Bf[i-1,j] <- length(genotype[genotype[,7] == "f" & genotype[,5] > 0, 5])
       #print(genotype)
       #
       alleles.juvf <- unlist(apply(as.data.frame(genotype)[genotype[,7] == "f", ], 1, function(dat){
         sample(x = dat[3:4], size = dat[5], replace = TRUE)})) #female alleles
        alleles.juvm <- unlist(apply(as.data.frame(genotype)[genotype[,7] == "m", ], 1, function(dat){
          sample(x = dat[3:4], size = dat[5], replace = TRUE)})) #
        gtjuv <- matrix(nrow = sum(b), ncol = 8, dimnames = list(NULL, c("id", "a1", "a2", "n.off", "surv", "sex", "z", "pm")))
        gtjuv[,1] <- paste(j, i, seq(from = Ntot[i-1, j]+1, length.out = sum(b)))
       gtjuv[,2] <- sample(alleles.juvf, sum(b))
       gtjuv[,3] <- sample(alleles.juvm, sum(b))
       gtjuv[,4] <- rep(0, sum(b))
       gtjuv[,6] <- sample(c("f", "m"), sum(b),replace = T, c(0.5, 0.5)) #offspring get sex
       #print(gtjuv)
       
       #
       #weigthed values for N
       wNf = wf*Nf[i-1, j] + wfm*Nm[i-1, j] #females might be more affected by other females
       wNm = wmf*Nf[i-1, j] + wm*Nm[i-1, j]#males may be more affected by other males
       #females:
       genotype[genotype[,7] == "f",6] = rbinom(Nf[i-1, j], 1, prob = s.f0*exp(-alphaf*(wNf/A)+env.a[i-1])/(1-s.f0*exp(-alphaf*(wNf/A))+s.f0*exp(-alphaf*(wNf/A)+env.a[i-1]))) #For now, dp is different between sexes
       #males:
       genotype[genotype[,7] == "m",6] = rbinom(Nm[i-1, j], 1, prob = s.m0*exp(-alpham*(wNm/A)+env.a[i-1])/(1-s.m0*exp(-alpham*(wNm/A))+s.m0*exp(-alpham*(wNm/A)+env.a[i-1])))
        
        #juveniles survive or die
       #females:
       gtjuv[gtjuv[,6] == "f",5] = rbinom(length(gtjuv[which(gtjuv[,6] == "f")]), 1, prob = s.jf0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])/(1-s.jf0*exp(-gamma*(Ntot[i-1, j]/A)) + s.jf0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])))
       #print(gtjuv)
       #males:
       gtjuv[gtjuv[,6] == "m",5] = rbinom(length(gtjuv[which(gtjuv[,6] == "m")]), 1, prob = s.jm0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])/(1-s.jm0*exp(-gamma*(Ntot[i-1, j]/A)) + s.jm0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])))
       #male juvs get z's
       gtjuv[gtjuv[,6] == "m", 7] <-rlnorm(length(gtjuv[which(gtjuv[,6] == "m")]), meanlog = -0.5*tau.m^2, sdlog = tau.m)
        #
        genotypeind = rbind(genotypeind, data.frame(subset(genotype, genotype[,6] == 1)),
                            data.frame(time = rep(i, sum(as.numeric(gtjuv[,5]))),
                                       subset(gtjuv, gtjuv[,5] == 1)))
        #print(genotypeind)
        #Popsize and allelfreq at i > 2 
        Ntot[i,j] = sum(as.numeric(genotype[,6])) + sum(as.numeric(gtjuv[,5]))
       Nf[i,j] = length(genotypeind$sex[genotypeind$sex == "f" & genotypeind$time == i])
        #print(Nf)
       Nm[i,j] = Ntot[i,j] - Nf[i,j]
        #print(Ntot)
        alive <- genotypeind[genotypeind$surv == "1" & genotypeind$time == i,]
        alivef <- genotypeind[genotypeind$sex == "f" & genotypeind$surv == "1" & genotypeind$time == i,]
        alivem <- genotypeind[genotypeind$sex == "m" & genotypeind$surv == "1" & genotypeind$time == i,]
        #print(alive)
        pool = with(alive, c(as.character(a1), as.character(a2)))
        poolf = with(alivef, c(as.character(a1), as.character(a2))) 
        poolm = with(alivem, c(as.character(a1), as.character(a2))) 
        p.a[i, j] = length(pool[pool == "a"])/length(pool)
        p.af[i, j] = length(poolf[poolf == "a"])/length(poolf)
        p.am[i, j] = length(poolm[poolm == "a"])/length(poolm)
        p.a[p.a <= 0] <- 0
        p.a[p.a >= 1] <- 1
        delta.p[i-1,j] = p.a[i,j] - p.a[i-1,j]
        a <- paste(alive$a1, alive$a2)
        genefreqaa[i,j] = length(a[a == "a a"])/Ntot[i, j]
        genefreqAA[i,j] = length(a[a == "A A"])/Ntot[i, j]
        genefreqAa[i,j] = length(a[a %in% c("A a", "a A")])/Ntot[i, j]
        
      }
      }
    }
  }
 
  return(list("u" = u, "env.a" = env.a, "env.j" = env.j, "env.r" = env.r, "p.a" = p.a, "Ntot" = Ntot, "Nf" = Nf, "Bm" = Bm, "Bf" = Bf,"delta.p" = delta.p, "genefreqAA" = genefreqAA, "genefreqAa" = genefreqAa, "genefreqaa" = genefreqaa, "meanbf" = meanbf, "varbf" = varbf,"meanbm" = meanbm, "varbm" = varbm,"p.af" = p.af,"p.am" = p.am))
}

gp.mod(nsim = 1, timesteps = 2, Nf0 = 43, Nm0 = 43, s.jf0 = 0.6, lambda.0 = 2, tau.m = 0.3, alphaf = 0.0101) 

#####model for monogamous populations####
monogam.mod <- function(nsim = 10, timesteps = 20, Nf0 = 40, Nm0 = 40, A = 1, alphaf  = 0, alpham = 0, beta = 0, gamma = 0, eps.a = 0, eps.j = 0, eps.r = 0, wf = 1, wm = 1, wfm = 1, wmf =1, tau.f = 0, tau.m = 0, lambda.0 = 2, s.jf0 = 0.3, s.jm0 = 0.3, s.f0 = 0.6, s.m0 = 0.7, pa0 = 0.5) {
  q = 1-pa0
  u = vector(mode = "numeric")
  env.a = vector(mode = "numeric")
  env.j = vector(mode = "numeric")
  env.r = vector(mode = "numeric")
  #vector for environmental variance. 
  for(i in 1:timesteps) { 
    u[i] = rnorm(n = 1, 0,1)
    env.a[i] = eps.a*u[i] #For survival, the last term vanish because of the logit function
    env.j[i] = eps.j*u[i]
    env.r[i] = eps.r*u[i] - 0.5*(eps.r^2)
  } 
  #Define matrixes we need:
  Ntot = matrix(data = rep(c(Nf0+Nm0, rep(NA, timesteps-1)), nsim), nrow = timesteps) #Total number of individuals at each timestep
  Nf = matrix(data = rep(c(Nf0, rep(NA, timesteps-1)), nsim), nrow = timesteps) #number of females
  Nm = matrix(data = rep(c(Nm0, rep(NA, timesteps-1)), nsim), nrow = timesteps)
  varbf <- matrix(nrow = timesteps, ncol = nsim)
  meanbf <- matrix(nrow = timesteps, ncol = nsim)
  varbm <- matrix(nrow = timesteps, ncol = nsim)
  meanbm <- matrix(nrow = timesteps, ncol = nsim)
  Bm = matrix(data = rep(rep(NA, timesteps), nsim), nrow = timesteps)
  Bf = matrix(data = rep(rep(NA, timesteps), nsim), nrow = timesteps)
  p.a <- matrix(data = rep(c(pa0, rep(NA, timesteps-1)), nsim), nrow = timesteps)
  delta.p <- matrix(data = rep(rep(NA, timesteps), nsim), nrow = timesteps)
  alleles2 <- 0 #this is where the resulting alleles for each j, i != 2 goes. Genotypes will be sampled from here
  alleles.juv <- 0 
  genefreqAA <- matrix(nrow = timesteps, ncol = nsim)
  genefreqaa <- matrix(nrow = timesteps, ncol = nsim)
  genefreqAa <- matrix(nrow = timesteps, ncol = nsim)
  for(j in 1:nsim){ #number of simulations
    genotypeind <- data.frame(time = numeric(),
                              id = character(),
                              a1 = character(),
                              a2 = character(),
                              n.off = numeric(),
                              surv = numeric(),
                              sex = character()
                              )
    for(i in 2:timesteps) { #Number of timesteps
      b <- 0
      if(i == 2) {
        genotype <- matrix(nrow = Ntot[i-1, j], ncol = 6, dimnames = list(NULL, c("id", "a1", "a2", "n.off", "surv", "sex")))
        genotype[,1] = paste(j, i-1, 1:Ntot[i-1, j])
        #alleles for parent generation.
        alleles2 <- c(replicate(round(2*q*Ntot[i-1, j]), "A"), replicate(2*pa0*Ntot[i-1, j], "a"))
        #print(alleles2)
        genotype[,2:3] = sample(alleles2, 2*Ntot[i-1, j])
        a <- paste(genotype[,2], genotype[,3])
        genefreqaa[i-1,j] = length(a[a == "a a"])/Ntot[i-1, j]
        genefreqAA[i-1,j] = length(a[a == "A A"])/Ntot[i-1, j]
        genefreqAa[i-1,j] = length(a[a %in% c("A a", "a A")])/Ntot[i-1, j]
        #sex need to be given
        sexc <- c(replicate(round(Nf0), "f"), replicate(Nm0, "m"))
        genotype[,6] = sample(sexc)
        #print(genotype)
        #females reproduce
        #poisson parameter: 
        lambda.pois = lambda.0*exp(-beta*(Ntot[i-1, j]/A)+env.r[i-1]) #
        #poisson-log-normal parameter:
       
        b = rpois(Nf[i-1, j], lambda = lambda.pois*exp(tau.f*rnorm(n = Nf[i-1, j], 0,1) - 0.5*(tau.f^2))) #
        
        #print(b)
        genotype[genotype[,6] == "f", 4] <- b 
        #print(genotype)
        geno1 <- genotype[genotype[,6] == "f" & genotype[,4] > 0,]
        
        if(length(b[b>0]) <= length(genotype[genotype[,6] == "m",1])) {
          gtjuv <- matrix(nrow = sum(b), ncol = 6, dimnames = list(NULL, c("id", "a1", "a2", "n.off", "surv", "sex")))
          m = sample(genotype[genotype[,6] == "m",1], length(b[b>0]), replace = F)
          genot <- cbind(matrix(geno1, ncol = 6), matrix(genotype[genotype[,1] %in% m,], ncol = 6)) #datafile containing reproducing parents
          gtjuv[,2] <- as.vector(unlist(apply(genot, 1, function(dat){
            sample(x = dat[2:3], size = dat[4], replace = TRUE)})))
          gtjuv[,3] <- as.vector(unlist(apply(genot, 1, function(dat){
            sample(x = dat[8:9], size = dat[4], replace = TRUE)})))
        }
        
        if(length(b[b>0]) > length(genotype[genotype[,6] == "m",1])) {
          genot <-cbind(matrix(geno1[sample(nrow(geno1), size = length(genotype[genotype[,6] == "m",1]), replace = F),],ncol = 6), matrix(genotype[genotype[,6] == "m",],ncol = 6))
          gtjuv <- matrix(nrow = sum(as.numeric(genot[,4])), ncol = 6, dimnames = list(NULL, c("id", "a1", "a2", "n.off", "surv", "sex")))
          gtjuv[,2] <- as.vector(unlist(apply(genot, 1, function(dat){
            sample(x = dat[2:3], size = dat[4], replace = TRUE)})))
          gtjuv[,3] <- as.vector(unlist(apply(genot, 1, function(dat){
            sample(x = dat[8:9], size = dat[4], replace = TRUE)})))
        }
        
        genotype[genotype[,1] %in% genot[,7],4] <- genot[,4] #the number of offspring per male
        genotype[!genotype[,1] %in% genot[,c(1,7)],4] <- 0 #females that first were given offspring, but didn't find a mate.
        meanbf[i-1,j] <- mean(as.numeric(genotype[genotype[,6] == "f",4]))
        varbf[i-1,j] <- var(as.numeric(genotype[genotype[,6] == "f",4]))
        meanbm[i-1,j] <- mean(as.numeric(genotype[genotype[,6] == "m",4]))
        varbm[i-1,j] <- var(as.numeric(genotype[genotype[,6] == "m",4]))
        
        Bm[i-1,j] <- length(genotype[genotype[,6] == "m" & genotype[,4] > 0, 4])
        Bf[i-1,j] <- length(genotype[genotype[,6] == "f" & genotype[,4] > 0, 4])
        #print(genotype)
       #print(genot)
  
        gtjuv[,1] <- paste(j, i, seq(from = Ntot[i-1, j]+1, length.out = sum(as.numeric(genot[,4]))))
        
        gtjuv[,4] <- 0
        gtjuv[,6] <- sample(c("f", "m"), sum(as.numeric(genot[,4])),replace = T, c(0.5, 0.5)) #offspring get sex
        #print(gtjuv)
        
        ##adults survive to i+1 or dies:
        #weigthed values for N. One for females, and one for males
        wNf = wf*Nf[i-1, j] + wfm*Nm[i-1, j] #females might be more affected by other females
        wNm = wmf*Nf[i-1, j] + wm*Nm[i-1, j]#males may be more affected by other males
        #females:
        genotype[genotype[,6] == "f",5] = rbinom(Nf[i-1, j], 1, prob = s.f0*exp(-alphaf*(wNf/A)+env.a[i-1])/(1-s.f0*exp(-alphaf*(wNf/A))+s.f0*exp(-alphaf*(wNf/A)+env.a[i-1]))) #For now, dp is different between sexes
        #males:
        genotype[genotype[,6] == "m",5] = rbinom(Nm[i-1, j], 1, prob = s.m0*exp(-alpham*(wNm/A)+env.a[i-1])/(1-s.m0*exp(-alpham*(wNm/A))+s.m0*exp(-alpham*(wNm/A)+env.a[i-1])))
        #print(genotype)
        #juveniles survive or die
        #females:
        gtjuv[gtjuv[,6] == "f",5] = rbinom(length(gtjuv[which(gtjuv[,6] == "f")]), 1, prob = s.jf0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])/(1-s.jf0*exp(-gamma*(Ntot[i-1, j]/A)) + s.jf0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])))
        #males:
        gtjuv[gtjuv[,6] == "m",5] = rbinom(length(gtjuv[which(gtjuv[,6] == "m")]), 1, prob = s.jm0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])/(1-s.jm0*exp(-gamma*(Ntot[i-1, j]/A)) + s.jm0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])))
        #print(gtjuv)
        
        #Gives the startpopulation present before reproduction: 
        genotypeind = rbind(genotypeind, data.frame(time = rep(i-1, Ntot[i-1, j]), genotype))
        
        #
        genotypeind = rbind(genotypeind, data.frame(time = rep(i, sum(as.numeric(genotype[,5]))),
                                                    subset(genotype, genotype[,5] == 1)),
                            data.frame(time = rep(i, sum(as.numeric(gtjuv[,5]))),
                                       subset(gtjuv, gtjuv[,5] == 1)))
        #print(genotypeind)
        #Popsize and allele freq
        Ntot[i,j] = sum(as.numeric(genotype[,5])) + sum(as.numeric(gtjuv[,5]))
        Nf[i,j] = length(genotypeind$sex[genotypeind$sex == "f" & genotypeind$time == i])
        # print(Nf)
        Nm[i,j] = Ntot[i,j] - Nf[i,j]
        alive <- genotypeind[genotypeind$surv == "1" & genotypeind$time == i,]
        #print(alive)
        pool = with(alive, c(as.character(a1), as.character(a2))) #
        #print(pool)
        p.a[i, j] = length(pool[pool == "a"])/length(pool)
        delta.p[i-1,j] = p.a[i,j] - p.a[i-1,j]
        a <- paste(alive$a1, alive$a2)
        genefreqaa[i,j] = length(a[a == "a a"])/Ntot[i, j]
        genefreqAA[i,j] = length(a[a == "A A"])/Ntot[i, j]
        genefreqAa[i,j] = length(a[a %in% c("A a", "a A")])/Ntot[i, j]
        
      }
      if(i > 2) {
        if(Nf[i-1, j] %in% c(0,NA) | Nm[i-1, j] %in% c(0, NA)) {break}
        if(!Nf[i-1, j] %in% c(0,NA) & !Nm[i-1, j] %in% c(0, NA)) {
          genotype = as.matrix(genotypeind[genotypeind$time == i-1,])
          genotype[,1] = i
          #
          lambda.pois = lambda.0*exp(-beta*(Ntot[i-1,j]/A)+env.r[i-1])
          
          b = rpois(Nf[i-1, j], lambda = lambda.pois*exp(tau.f*rnorm(n = Nf[i-1, j], 0,1) - 0.5*(tau.f^2))) #
          while(sum(b) == 0) {
            b = rpois(Nf[i-1, j], lambda = lambda.pois*exp(tau.f*rnorm(n = Nf[i-1, j], 0,1) - 0.5*(tau.f^2)))
          } #in small pops, there is a small chanse that no offspring may be born, which will break down the rest of the sim.
          genotype[genotype[,7] == "f", 5] <- b #
          
          geno1 <- genotype[genotype[,7] == "f" & genotype[,5] > 0,]
          #print(geno1)
          #print(genotype)
          if(length(b[b>0]) <= length(genotype[genotype[,7] == "m",1])) {
            gtjuv <- matrix(nrow = sum(b), ncol = 6, dimnames = list(NULL, c("id", "a1", "a2", "n.off", "surv", "sex")))
            m = sample(genotype[genotype[,7] == "m",2], length(b[b>0]), replace = F)
            genot <- cbind(matrix(geno1, ncol = 7), matrix(genotype[genotype[,2] %in% m,], ncol = 7)) #datafile containing reproducing parents
            gtjuv[,2] <- as.vector(unlist(apply(genot, 1, function(dat){
              sample(x = dat[3:4], size = dat[5], replace = TRUE)})))
            gtjuv[,3] <- as.vector(unlist(apply(genot, 1, function(dat){
              sample(x = dat[10:11], size = dat[5], replace = TRUE)})))
          }
          
          if(length(b[b>0]) > length(genotype[genotype[,7] == "m",2])) {
            genot <-cbind(matrix(geno1[sample(nrow(geno1), size = length(genotype[genotype[,7] == "m",2]), replace = F),], ncol = 7), matrix(genotype[genotype[,7] == "m",], ncol = 7))
            #print(genot)
            gtjuv <- matrix(nrow = sum(as.numeric(genot[,5])), ncol = 6, dimnames = list(NULL, c("id", "a1", "a2", "n.off", "surv", "sex")))
            gtjuv[,2] <- as.vector(unlist(apply(genot, 1, function(dat){
              sample(x = dat[3:4], size = dat[5], replace = TRUE)})))
            gtjuv[,3] <- as.vector(unlist(apply(genot, 1, function(dat){
              sample(x = dat[10:11], size = dat[5], replace = TRUE)})))
          }
          #print(gtjuv)
          genotype[genotype[,2] %in% genot[,9],5] <- genot[,5] #the number of offspring per male
          genotype[!genotype[,2] %in% genot[,c(2,9)],5] <- 0 #females that first were given offspring, but didn't find a mate.
          meanbf[i-1,j] <- mean(as.numeric(genotype[genotype[,7] == "f",5]))
          varbf[i-1,j] <- var(as.numeric(genotype[genotype[,7] == "f",5]))
          meanbm[i-1,j] <- mean(as.numeric(genotype[genotype[,7] == "m",5]))
          varbm[i-1,j] <- var(as.numeric(genotype[genotype[,7] == "m",5]))
          
          Bm[i-1,j] <- length(genotype[genotype[,7] == "m" & genotype[,5] > 0, 5]) 
          Bf[i-1,j] <- length(genotype[genotype[,7] == "f" & genotype[,5] > 0, 5])
          #print(genotype)
          
          gtjuv[,1] <- paste(j, i, seq(from = Ntot[i-1, j]+1, length.out = sum(as.numeric(genot[,5]))))
          gtjuv[,4] <- 0
          gtjuv[,6] <- sample(c("f", "m"), sum(as.numeric(genot[,5])),replace = T, c(0.5, 0.5)) #offspring get sex
          #print(gtjuv)
          
          ##
          #weigthed values for N
          wNf = wf*Nf[i-1, j] + wfm*Nm[i-1, j] #females might be more affected by other females
          wNm = wmf*Nf[i-1, j] + wm*Nm[i-1, j]#males may be more affected by other males
          #females:
          genotype[genotype[,7] == "f",6] = rbinom(Nf[i-1, j], 1, prob = s.f0*exp(-alphaf*(wNf/A)+env.a[i-1])/(1-s.f0*exp(-alphaf*(wNf/A))+s.f0*exp(-alphaf*(wNf/A)+env.a[i-1]))) #For now, dp is different between sexes
          #males:
          genotype[genotype[,7] == "m",6] = rbinom(Nm[i-1, j], 1, prob = s.m0*exp(-alpham*(wNm/A)+env.a[i-1])/(1-s.m0*exp(-alpham*(wNm/A))+s.m0*exp(-alpham*(wNm/A)+env.a[i-1])))
          
          #
          #females:
          gtjuv[gtjuv[,6] == "f",5] = rbinom(length(gtjuv[which(gtjuv[,6] == "f")]), 1, prob = s.jf0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])/(1-s.jf0*exp(-gamma*(Ntot[i-1, j]/A)) + s.jf0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])))
          #print(gtjuv)
          #males:
          gtjuv[gtjuv[,6] == "m",5] = rbinom(length(gtjuv[which(gtjuv[,6] == "m")]), 1, prob = s.jm0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])/(1-s.jm0*exp(-gamma*(Ntot[i-1, j]/A)) + s.jm0*exp(-gamma*(Ntot[i-1, j]/A)+env.j[i-1])))
          
          #
          genotypeind = rbind(genotypeind, data.frame(subset(genotype, genotype[,6] == 1)),
                              data.frame(time = rep(i, sum(as.numeric(gtjuv[,5]))),
                                         subset(gtjuv, gtjuv[,5] == 1)))
          #print(genotypeind)
          #
          Ntot[i,j] = sum(as.numeric(genotype[,6])) + sum(as.numeric(gtjuv[,5]))
          Nf[i,j] = length(genotypeind$sex[genotypeind$sex == "f" & genotypeind$time == i])
          #print(Nf)
          Nm[i,j] = Ntot[i,j] - Nf[i,j]
          #print(Ntot)
          alive <- genotypeind[genotypeind$surv == "1" & genotypeind$time == i,]
          #print(alive)
          pool = with(alive, c(as.character(a1), as.character(a2)))
          p.a[i, j] = length(pool[pool == "a"])/length(pool)
          p.a[p.a <= 0] <- 0
          p.a[p.a >= 1] <- 1
          delta.p[i-1,j] = p.a[i,j] - p.a[i-1,j]
          a <- paste(alive$a1, alive$a2)
          genefreqaa[i,j] = length(a[a == "a a"])/Ntot[i, j]
          genefreqAA[i,j] = length(a[a == "A A"])/Ntot[i, j]
          genefreqAa[i,j] = length(a[a %in% c("A a", "a A")])/Ntot[i, j]
          
        }
      }
    }
  }
  
  return(list("u" = u, "env.a" = env.a, "env.j" = env.j, "env.r" = env.r, "p.a" = p.a, "Ntot" = Ntot, "Nf" = Nf, "Bm" = Bm, "Bf" = Bf,"delta.p" = delta.p, "genefreqAA" = genefreqAA, "genefreqAa" = genefreqAa, "genefreqaa" = genefreqaa, "meanbf" = meanbf, "varbf" = varbf,"meanbm" = meanbm, "varbm" = varbm))
}



