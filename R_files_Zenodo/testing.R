# load libraries

library(diversitree)
library(phytools)

# load tree and rescale it to edge of 1.0

tree <- read.nexus("my.tree.nex")

mbt <- max(branching.times(tree))

tree$edge.length <- tree$edge.length / mbt

taxnames <- tree$tip.label

# empty lists for the loop

simtr <- list()
simtr.sd <- list()
p <- list()
p.c <- list()
lik.c <- list()
mle.c <- list()
xr <- list()
linear.x <- list()
c.cc <- list()
p.ll <- list()
lik.l <- list()
mle.l <- list()

# table to save the results

res <- data.frame(matrix(, 100, 18, dimnames = list(seq(1:100) , c("constant.lambda.1", "constant.mu.1", "constant.diffusion", "constant.lambda.2", "constant.mu.2", "constant.k", "constant.lnLik", "constant.AIC", "linear.lambda.1.c", "linear.lambda.1.m", "linear.mu.1", "linear.diffusion", "linear.lambda.2.c", "linear.lambda.2.m", "linear.mu.2", "linear.k", "linear.lnLik", "linear.AIC"))))

# constraining drift to 0 and assuming that both partitions have the same
# difussion coefficient

nodrift <- function(f)
  constrain(f, drift.1 ~ 0, drift.2 ~ 0, diffusion.2 ~ diffusion.1)

# control parameters

control <- list(parscale = 0.1, reltol = 0.001)

# the looooop

for(i in 1:100){
  
  # simulate trait under BM (this is for sigma^2 = 0.0125; use different values)
  
  simtr[[i]] <- fastBM(tree, sig2 = 0.0125)
  simtr.sd[[i]] <- 0.001 * sd(simtr[[i]])
  
  # starting point

  p[[i]] <- starting.point.quasse(tree, simtr[[i]])

  # starting point (constant speciation)
  
  p.c[[i]] <- c(p[[i]], p[[i]][1:2])

  # likelihood function for constant speciation and constant extinction
  
  lik.c[[i]] <- make.quasse.split(tree, simtr[[i]], simtr.sd[[i]], constant.x, constant.x, 204, Inf, sampling.f = 192/247)
  
  # constant model

  mle.c[[i]] <- find.mle(nodrift(lik.c[[i]]), p.c[[i]], lower = 0, control = control, verbose = 0)

  # piecewise function (flat outside the range)
  
  xr[[i]] <- range(simtr[[i]]) + c(-1, 1) * 20 * p[[i]]["diffusion"]
  
  linear.x[[i]] <- make.linear.x(xr[[i]][1], xr[[i]][2])
  
  # parameters for the linear model
  
  c.cc[[i]] <- coef(mle.c[[i]])
  p.ll[[i]] <- c(c.cc[[i]][1], 0, c.cc[[i]][2:4], 0, c.cc[[i]][5])

  # likelihood function for linear speciation and constant extinction
  
  lik.l[[i]] <- make.quasse.split(tree, simtr[[i]], simtr.sd[[i]], linear.x[[i]], constant.x, 204, Inf, sampling.f = 192/247)
  
  # linear model
  
  mle.l[[i]] <- find.mle(nodrift(lik.l[[i]]), p.ll[[i]], control = control, verbose = 0)

  # fill table
  
  res[i, ] <- c(mle.c[[i]]$par, length(mle.c[[i]]$par), mle.c[[i]]$lnLik, AIC(mle.c[[i]]), mle.l[[i]]$par, length(mle.l[[i]]$par), mle.l[[i]]$lnLik, AIC(mle.l[[i]]))

}

write.csv(res, file = "quasse.test.split_0125.csv", row.names = F)
