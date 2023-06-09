library("PBD")
library(geiger)
library(ape)
library(diversitree)
library(lattice)
library(misc3d)
library(plot3D)

# simulate protracted speciation tree

#birth_range <- c(0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)
#convert_range <- c(0.03, 0.05, 0.07, 0.1, 0.3, 0.5, 0.7, 1.0, 5.0, 10.0)
#death_range <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
birth_range <- c(0.5, 0.55, 0.6, 0.65, 0.7)
convert_range <- c(0.01, 0.06, 0.11, 0.16, 0.21)
death_range <- c(0.25, 0.3, 0.35, 0.4, 0.45)

# how many parameters to use from each variable
n <- 5

# generate a table to document each combination
parameters_summary <- data.frame(b=numeric(n^3), x=numeric(n^3), mu=numeric(n^3), spe_rate=numeric(n^3), ext_rate=numeric(n^3))

#fill in birth rate
for (i in 1:n){
  parameters_summary$b[(1+(i-1)*n^2):(n^2+(i-1)*n^2)] <- rep(birth_range[i], n^2)
}
#fill in convert rate
xvector <- numeric(n^2)
for (i in 1:n){
  xvector[(1+(i-1)*n):(n+(i-1)*n)] <- rep(convert_range[i], n)
}
for (i in 1:n){
  parameters_summary$x[(1+(i-1)*n^2):(n^2+(i-1)*n^2)] <- xvector
}
#fill in extinction rate
for (i in 1:n){
  parameters_summary$mu[(1+(i-1)*n^2):(n^2+(i-1)*n^2)] <- rep(death_range[1:n],n)
}


#test run, 125 combinations of PBD trees

# generate empty result tables for 125 comnination of PBD trees, each combination has 5 simulations
N <- n^3 # number of combinations
m <- 100  # replications per combination

# a list that holds all combinations of replications, in test case, 8 list,  each has 10 sublists
All_trees <- replicate(n^3, list()) 
for (b in 1:N){
  All_trees[[b]] <- replicate(m, list())
}


counter <- 1
time <- 15


for (i in 1:n){
  b <- birth_range[i]
  for (j in 1:n){
    x <- convert_range[j]
    for (k in 1:n) {
      mu <- death_range[k]
      for (l in 1:m){
        simulation <- pbd_sim(c(b, x, b, mu, mu),time)
        All_trees[[counter]][[l]] <- simulation[[2]]
      }
      counter <- counter + 1
    }
  }
}

#125 * 5 simulations for time 10 takes less than 1 min
#125 * 5 simulations for time 20 takes several hours, up to 12?

spe_num <- numeric(length(All_trees))
for (o in 1:length(All_trees)){
  spe_num[o] <- length(All_trees[[o]][[4]]$tip.label)
}

hist(spe_num)
table(spe_num)


#test data 125 Birth-death estimates, five replicates, pretty fast. A few minuets.

test <- n^3
spe_rate <- numeric(test)
ext_rate <- numeric(test)

for (p in 1:test){
  spe_temp <- numeric(m)
  ext_temp <- numeric(m)
  for (q in 1:m){
    brt <- branching.times(All_trees[[p]][[q]])
    bd_estimate <- bd_ML(brts = brt, cond=1)
    spe_temp[q] <- bd_estimate$lambda0
    ext_temp[q] <- bd_estimate$mu0
  }
  spe_rate[p] <- mean(spe_temp)
  ext_rate[p] <- mean(ext_temp)
}


parameters_summary$spe_rate <- round(spe_rate, 2)
parameters_summary$ext_rate <- round(ext_rate, 2)



#contourplot(spe_rate ~ b * mu, data=parameters_summary, region = TRUE, xlab="Population initiation rate", ylab="Population extinction rate", main = "Distribution of speciation rate of 125 trees")
#contourplot(spe_rate ~ b * x, data=parameters_summary, region = TRUE, xlab="Population initiation rate", ylab="Population conversion rate", main = "Distribution of speciation rate of 125 trees")
#contourplot(spe_rate ~ mu * x, data=parameters_summary, region = TRUE, xlab="Population extinction rate", ylab="Population conversion rate", main = "Distribution of speciation rate of 125 trees")

#Speciation rate
levelplot(spe_rate ~ b * mu, data=parameters_summary, region = TRUE, xlab="Population initiation rate", ylab="Population extinction rate", main = "Distribution of speciation rate of 125 trees")
levelplot(spe_rate ~ b * x, data=parameters_summary, region = TRUE, xlab="Population initiation rate", ylab="Population conversion rate", main = "Distribution of speciation rate of 125 trees")
levelplot(spe_rate ~ mu * x, data=parameters_summary, region = TRUE, xlab="Population extinction rate", ylab="Population conversion rate", main = "Distribution of speciation rate of 125 trees")
#contourplot(spe_rate ~ mu * x, data=parameters_summary, region = TRUE, xlab="Population extinction rate", ylab="Population conversion rate", main = "Distribution of speciation rate of 125 trees")

#extinction rate
levelplot(ext_rate ~ b * mu, data=parameters_summary, region = TRUE, xlab="Population initiation rate", ylab="Population extinction rate", main = "Distribution of extinction rate of 125 trees")
levelplot(ext_rate ~ b * x, data=parameters_summary, region = TRUE, xlab="Population initiation rate", ylab="Population conversion rate", main = "Distribution of extinction rate of 125 trees")
levelplot(ext_rate ~ mu * x, data=parameters_summary, region = TRUE, xlab="Population extinction rate", ylab="Population conversion rate", main = "Distribution of extinction rate of 125 trees")


levelplot(ext_rate ~ mu * x, data=parameters_summary, at =seq(0, 0.03,length=100), region = TRUE, xlab="Population extinction rate", ylab="Population conversion rate", main = "Distribution of extinction rate of 125 trees")
contourplot(ext_rate ~ mu * x, data=parameters_summary, region = TRUE, xlab="Population extinction rate", ylab="Population conversion rate", main = "Distribution of extinction rate of 125 trees")
