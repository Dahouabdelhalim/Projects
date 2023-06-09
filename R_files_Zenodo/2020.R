
### Haploid, hermaphrodite

#=================#
# ALLOW EVOLUTION #
#=================#


## Define functions

pref  <- function(a,b) exp(-4*((a-b)^2))


## Define parameters

R1  <- 4 # potential reproductive rate
K1  <- 1000 # carrying capacity
R2  <- 4
K2  <- 1000


## Define initial conditions

n01 <- 1000 # initial abundance
z01 <- 0   # initial trait value

n02 <- 1000
z02 <- 0


t <- proc.time()
  
#====================#
# Run the simulation #
#====================#
  
  
## produce individuals
  
if (n01==0) {
  sp.1 <- data.frame(NULL)
} else {
  sp.1 <- data.frame(z_01=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_02=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_03=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_04=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_05=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_06=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_07=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_08=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_09=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_10=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1)
}
  
if (n02==0) {
  sp.2 <- data.frame(NULL)
} else {
  sp.2 <- data.frame(z_01=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_02=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_03=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_04=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_05=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_06=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_07=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_08=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_09=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_10=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1)
}
  
# take stats
  
if (nrow(sp.1)==0) {
  stats.1 <- c(0, NA, NA)
} else {
  stats.1 <- c(nrow(sp.1), mean(apply(sp.1[,1:10], 1, mean)), sd(apply(sp.1[,1:10], 1, mean)))
}
  
if (nrow(sp.2)==0) {
  stats.2 <- c(0, NA, NA)
} else {
  stats.2 <- c(nrow(sp.2), mean(apply(sp.2[,1:10], 1, mean)), sd(apply(sp.2[,1:10], 1, mean)))
}
  
result <- c(stats.1, stats.2)
  
# condition parameters
  
generation <- 0
extinction <- 0
divergence <- 0
  
for(i in 1:150) {
    
  ## produce next generation
    
  # mating process
    
  if(nrow(sp.1)!=0) { sp.1$sp <- 1 } # species identity
  if(nrow(sp.2)!=0) { sp.2$sp <- 2 }
    
  all <- rbind(sp.1, sp.2)
    
  all$trait <- apply(all[,1:10], 1, mean) # trait value  
  MPM <- outer(all$trait, all$trait, pref) # mate preference matrix
  diag(MPM) <- 0
    
  mate.id <- apply(MPM, 1, function(x) sample(1:nrow(all), 1, prob=x)) # choose a mate (chonsen one acts as male)
  mate.m  <- all[mate.id, ]
    
  consp.pairs <- which(all$sp*mate.m$sp!=2) # rows of conspecific mating
  rep.all  <- all[consp.pairs, ]     # conspecifically mating individuals (acting as female)
  mom.1 <- rep.all[rep.all$sp==1, ]
  mom.2 <- rep.all[rep.all$sp==2, ]
  rep.mate <- mate.m[consp.pairs, ]  # their mates
  dad.1 <- rep.mate[rep.mate$sp==1, ]
  dad.2 <- rep.mate[rep.mate$sp==2, ]
    
  # offspring
    
  if(nrow(mom.1)==0) {
    sp.1 <- data.frame(NULL)
  } else {
    z_01 <- as.vector(apply(cbind(mom.1$z_01, dad.1$z_01), 1, function(x) sample(x, replace=T, R1)))
    z_02 <- as.vector(apply(cbind(mom.1$z_02, dad.1$z_02), 1, function(x) sample(x, replace=T, R1)))
    z_03 <- as.vector(apply(cbind(mom.1$z_03, dad.1$z_03), 1, function(x) sample(x, replace=T, R1)))
    z_04 <- as.vector(apply(cbind(mom.1$z_04, dad.1$z_04), 1, function(x) sample(x, replace=T, R1)))
    z_05 <- as.vector(apply(cbind(mom.1$z_05, dad.1$z_05), 1, function(x) sample(x, replace=T, R1)))
    z_06 <- as.vector(apply(cbind(mom.1$z_06, dad.1$z_06), 1, function(x) sample(x, replace=T, R1)))
    z_07 <- as.vector(apply(cbind(mom.1$z_07, dad.1$z_07), 1, function(x) sample(x, replace=T, R1)))
    z_08 <- as.vector(apply(cbind(mom.1$z_08, dad.1$z_08), 1, function(x) sample(x, replace=T, R1)))
    z_09 <- as.vector(apply(cbind(mom.1$z_09, dad.1$z_09), 1, function(x) sample(x, replace=T, R1)))
    z_10 <- as.vector(apply(cbind(mom.1$z_10, dad.1$z_10), 1, function(x) sample(x, replace=T, R1)))
    juvenile <- data.frame(z_01, z_02, z_03, z_04, z_05, z_06, z_07, z_08, z_09, z_10)
    n.next <- (nrow(juvenile)/R1)*(R1*K1/(K1 + (nrow(juvenile)/R1)*(R1 - 1)))
      
    if(nrow(juvenile) > n.next) {
      survivors <- sample(1:nrow(juvenile), n.next, replace=F)
      sp.1 <- juvenile[survivors, ]
    } else {
      sp.1 <- juvenile
    }
  }
    
  if(nrow(mom.2)==0) {
    sp.2 <- data.frame(NULL)
  } else {
    z_01 <- as.vector(apply(cbind(mom.2$z_01, dad.2$z_01), 1, function(x) sample(x, replace=T, R2)))
    z_02 <- as.vector(apply(cbind(mom.2$z_02, dad.2$z_02), 1, function(x) sample(x, replace=T, R2)))
    z_03 <- as.vector(apply(cbind(mom.2$z_03, dad.2$z_03), 1, function(x) sample(x, replace=T, R2)))
    z_04 <- as.vector(apply(cbind(mom.2$z_04, dad.2$z_04), 1, function(x) sample(x, replace=T, R2)))
    z_05 <- as.vector(apply(cbind(mom.2$z_05, dad.2$z_05), 1, function(x) sample(x, replace=T, R2)))
    z_06 <- as.vector(apply(cbind(mom.2$z_06, dad.2$z_06), 1, function(x) sample(x, replace=T, R2)))
    z_07 <- as.vector(apply(cbind(mom.2$z_07, dad.2$z_07), 1, function(x) sample(x, replace=T, R2)))
    z_08 <- as.vector(apply(cbind(mom.2$z_08, dad.2$z_08), 1, function(x) sample(x, replace=T, R2)))
    z_09 <- as.vector(apply(cbind(mom.2$z_09, dad.2$z_09), 1, function(x) sample(x, replace=T, R2)))
    z_10 <- as.vector(apply(cbind(mom.2$z_10, dad.2$z_10), 1, function(x) sample(x, replace=T, R2)))
    juvenile <- data.frame(z_01, z_02, z_03, z_04, z_05, z_06, z_07, z_08, z_09, z_10)
    n.next <- (nrow(juvenile)/R2)*(R2*K2/(K2 + (nrow(juvenile)/R2)*(R2 - 1)))
      
    if(nrow(juvenile) > n.next) {
      survivors <- sample(1:nrow(juvenile), n.next, replace=F)
      sp.2 <- juvenile[survivors, ]
    } else {
      sp.2 <- juvenile
    }
  }
    
  ## take stats
    
  if (nrow(sp.1)==0) {
    stats.1 <- c(0, NA, NA)
  } else {
    stats.1 <- c(nrow(sp.1), mean(apply(sp.1[,1:10], 1, mean)), sd(apply(sp.1[,1:10], 1, mean)))
  }
    
  if (nrow(sp.2)==0) {
    stats.2 <- c(0, NA, NA)
  } else {
    stats.2 <- c(nrow(sp.2), mean(apply(sp.2[,1:10], 1, mean)), sd(apply(sp.2[,1:10], 1, mean)))
  }
    
  result <- rbind(result, c(stats.1, stats.2))
    
  # change the condition parameters
    
  generation <- generation + 1
  if(nrow(sp.1)*nrow(sp.2)==0) {
    extinction <- 1
  } else {
    if(abs(stats.1[2] - stats.2[2])==2) {divergence <- 1}
  }
    
}
proc.time() - t

colnames(result) <- c("sp.1_ab", "sp.1_trait", "sp.1_trait.sd",
                      "sp.2_ab", "sp.2_trait", "sp.2_trait.sd")
RCD.result <- result

plot(1:151, RCD.result[,"sp.1_ab"], ylim=c(0,1050), type="l", col="blue")
lines(1:151, RCD.result[, "sp.2_ab"], col="red")

plot(1:151, RCD.result[,"sp.1_trait"], ylim=c(-1,1), type="l", col="blue")
lines(1:151, RCD.result[, "sp.2_trait"], col="red")

#save(RCD.result, file="RCD.result")

#-----------------------------------------------------------------------------------

  
#=======================#
# *NOT* ALLOW EVOLUTION #
#=======================#

  
## Define functions
  
pref  <- function(a,b) exp(-4*((a-b)^2))


## Define parameters

R1  <- 4 # potential reproductive rate
K1  <- 1000 # carryng capacity
R2  <- 4
K2  <- 1000


## Define initial conditions

n01 <- 1000 # initial abundance
z01 <- 0   # initial trait value

n02 <- 1000
z02 <- 0

t <- proc.time()
  
#====================#
# Run the simulation #
#====================#
  
  
## produce individuals
  
if (n01==0) {
  sp.1 <- data.frame(NULL)
} else {
  sp.1 <- data.frame(z_01=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_02=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_03=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_04=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_05=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_06=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_07=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_08=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_09=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1,
                     z_10=2*rbinom(n01, 1, 0.5*z01 + 0.5) - 1)
}
  
if (n02==0) {
  sp.2 <- data.frame(NULL)
} else {
  sp.2 <- data.frame(z_01=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_02=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_03=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_04=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_05=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_06=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_07=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_08=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_09=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1,
                     z_10=2*rbinom(n02, 1, 0.5*z02 + 0.5) - 1)
}
  
# take stats
  
if (nrow(sp.1)==0) {
  stats.1 <- c(0, NA, NA)
} else {
  stats.1 <- c(nrow(sp.1), mean(apply(sp.1[,1:10], 1, mean)), sd(apply(sp.1[,1:10], 1, mean)))
}
  
if (nrow(sp.2)==0) {
  stats.2 <- c(0, NA, NA)
} else {
  stats.2 <- c(nrow(sp.2), mean(apply(sp.2[,1:10], 1, mean)), sd(apply(sp.2[,1:10], 1, mean)))
}
  
result <- c(stats.1, stats.2)
  
# condition parameters
  
generation <- 0
extinction <- 0
divergence <- 0
  
for (i in 1:150) {
    
  ## produce next generation
    
  # mating process
    
  if(nrow(sp.1)!=0) { sp.1$sp <- 1 } # species identity
  if(nrow(sp.2)!=0) { sp.2$sp <- 2 }
    
  all <- rbind(sp.1, sp.2)
    
  all$trait <- apply(all[,1:10], 1, mean) # trait value  
  MPM <- outer(all$trait, all$trait, pref) # mate preference matrix
  diag(MPM) <- 0
    
  mate.id <- apply(MPM, 1, function(x) sample(1:nrow(all), 1, prob=x)) # choose a mate (chonsen one acts as male)
  mate.m  <- all[mate.id, ]
    
  consp.pairs <- which(all$sp*mate.m$sp!=2) # rows of conspecific mating
  rep.all  <- all[consp.pairs, ]     # conspecifically mating individuals (acting as female)
  mom.1 <- rep.all[rep.all$sp==1, ]
  mom.2 <- rep.all[rep.all$sp==2, ]
  rep.mate <- mate.m[consp.pairs, ]  # their mates
  dad.1 <- rep.mate[rep.mate$sp==1, ]
  dad.2 <- rep.mate[rep.mate$sp==2, ]
    
  # offspring
    
  if(nrow(mom.1)==0) {
    sp.1 <- data.frame(NULL)
  } else {
    z_01=2*rbinom(R1*nrow(mom.1), 1, 0.5*z01 + 0.5) - 1
    z_02=2*rbinom(R1*nrow(mom.1), 1, 0.5*z01 + 0.5) - 1
    z_03=2*rbinom(R1*nrow(mom.1), 1, 0.5*z01 + 0.5) - 1
    z_04=2*rbinom(R1*nrow(mom.1), 1, 0.5*z01 + 0.5) - 1
    z_05=2*rbinom(R1*nrow(mom.1), 1, 0.5*z01 + 0.5) - 1
    z_06=2*rbinom(R1*nrow(mom.1), 1, 0.5*z01 + 0.5) - 1
    z_07=2*rbinom(R1*nrow(mom.1), 1, 0.5*z01 + 0.5) - 1
    z_08=2*rbinom(R1*nrow(mom.1), 1, 0.5*z01 + 0.5) - 1
    z_09=2*rbinom(R1*nrow(mom.1), 1, 0.5*z01 + 0.5) - 1
    z_10=2*rbinom(R1*nrow(mom.1), 1, 0.5*z01 + 0.5) - 1
    juvenile <- data.frame(z_01, z_02, z_03, z_04, z_05, z_06, z_07, z_08, z_09, z_10)
    n.next <- (nrow(juvenile)/R1)*(R1*K1/(K1 + (nrow(juvenile)/R1)*(R1 - 1)))
      
    if(nrow(juvenile) > n.next) {
      survivors <- sample(1:nrow(juvenile), n.next, replace=F)
      sp.1 <- juvenile[survivors, ]
    } else {
      sp.1 <- juvenile
    }
  }
    
  if(nrow(mom.2)==0) {
    sp.2 <- data.frame(NULL)
  } else {
    z_01=2*rbinom(R2*nrow(mom.2), 1, 0.5*z02 + 0.5) - 1
    z_02=2*rbinom(R2*nrow(mom.2), 1, 0.5*z02 + 0.5) - 1
    z_03=2*rbinom(R2*nrow(mom.2), 1, 0.5*z02 + 0.5) - 1
    z_04=2*rbinom(R2*nrow(mom.2), 1, 0.5*z02 + 0.5) - 1
    z_05=2*rbinom(R2*nrow(mom.2), 1, 0.5*z02 + 0.5) - 1
    z_06=2*rbinom(R2*nrow(mom.2), 1, 0.5*z02 + 0.5) - 1
    z_07=2*rbinom(R2*nrow(mom.2), 1, 0.5*z02 + 0.5) - 1
    z_08=2*rbinom(R2*nrow(mom.2), 1, 0.5*z02 + 0.5) - 1
    z_09=2*rbinom(R2*nrow(mom.2), 1, 0.5*z02 + 0.5) - 1
    z_10=2*rbinom(R2*nrow(mom.2), 1, 0.5*z02 + 0.5) - 1
    juvenile <- data.frame(z_01, z_02, z_03, z_04, z_05, z_06, z_07, z_08, z_09, z_10)
    n.next <- (nrow(juvenile)/R2)*(R2*K2/(K2 + (nrow(juvenile)/R2)*(R2 - 1)))
      
    if(nrow(juvenile) > n.next) {
      survivors <- sample(1:nrow(juvenile), n.next, replace=F)
      sp.2 <- juvenile[survivors, ]
    } else {
      sp.2 <- juvenile
    }
  }
    
  ## take stats
    
  if (nrow(sp.1)==0) {
    stats.1 <- c(0, NA, NA)
  } else {
    stats.1 <- c(nrow(sp.1), mean(apply(sp.1[,1:10], 1, mean)), sd(apply(sp.1[,1:10], 1, mean)))
  }
    
  if (nrow(sp.2)==0) {
    stats.2 <- c(0, NA, NA)
  } else {
    stats.2 <- c(nrow(sp.2), mean(apply(sp.2[,1:10], 1, mean)), sd(apply(sp.2[,1:10], 1, mean)))
  }
    
  result <- rbind(result, c(stats.1, stats.2))
    
  # change the condition parameters
    
  generation <- generation + 1
  if(nrow(sp.1)*nrow(sp.2)==0) {
    extinction <- 1
  } else {
    if(abs(stats.1[2] - stats.2[2])==2) {divergence <- 1}
  }

}
proc.time() - t

colnames(result) <- c("sp.1_ab", "sp.1_trait", "sp.1_trait.sd",
                      "sp.2_ab", "sp.2_trait", "sp.2_trait.sd")

extinction.result <- result

plot(1:151, extinction.result[,"sp.1_ab"], ylim=c(0,1050), type="l", col="blue")
lines(1:151, extinction.result[, "sp.2_ab"], col="red")

plot(1:151, extinction.result[,"sp.1_trait"], ylim=c(-1,1), type="l", col="blue")
lines(1:151, extinction.result[, "sp.2_trait"], col="red")

#save(extinction.result, file="extinction.result")

