#simulate bivariate normal distribution
library(MASS)
set.seed(123)

# parameters for univariate normal distributions
# Number of random samples
N <- 5
#mean, sd of 1st variable
mu1 <- 0; s1 <- 1
#mean, sd of 2nd variable
mu2 <- 0; s2 <- 1

# Parameters for bivariate normal distribution
rho <- 0.40 # correlation coefficient
mu <- c(mu1,mu2) # means
sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),2) # covariance matrix

#make empty dataframe to gather results of many empirical tests
df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df) <- c("sim", "emp.rho", "95CI_lower", "95CI_upper","emp.p")

#number of empirical tests, to determine power
sims <- 1000

for(i in 1:sims){
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X1","bvn1_X2")

  #empirical test of null hypotheses (rho = 0)
  emp <- cor.test(bvn1[,1],bvn1[,2], method=c("pearson"))
  emp.rho <- emp$estimate
  low <- emp$conf.int[1]
  high <- emp$conf.int[2]
  emp.p <- emp$p.value
  df[nrow(df)+1,] <- c(i,emp.rho,low,high,emp.p)
}

mean(df[,2], na.rm=T)
mean(df[,3])
mean(df[,4])

#fraction rejecting Ho -- power
rejects <- nrow(df[df$emp.p<=0.05,])/nrow(df)
rejects