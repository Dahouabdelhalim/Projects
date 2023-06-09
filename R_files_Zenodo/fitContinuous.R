library(geiger)
library(ape)
phy <- read.tree("C_tree.nwk")
dat <- read.csv("Continuous.csv", header = T, row.names = 1)

BM <- fitContinuous(phy, dat, SE = 0,
              model = "BM",
              bounds= list(), control = list(method = "L-BFGS-B",
                                             niter = 1000, FAIL = 1e+200, hessian = FALSE, CI = 0.95), ncores=NULL)

OU <- fitContinuous(phy, dat, SE = 0,
                    model = "OU",
                    bounds= list(), control = list(method = "L-BFGS-B",
                                                   niter = 1000, FAIL = 1e+200, hessian = FALSE, CI = 0.95), ncores=NULL)

EB <- fitContinuous(phy, dat, SE = 0,
                    model = "EB",
                    bounds= list(), control = list(method = "L-BFGS-B",
                                                   niter = 1000, FAIL = 1e+200, hessian = FALSE, CI = 0.95), ncores=NULL)

White <- fitContinuous(phy, dat, SE = 0,
                    model = "white",
                    bounds= list(), control = list(method = "L-BFGS-B",
                                                   niter = 1000, FAIL = 1e+200, hessian = FALSE, CI = 0.95), ncores=NULL)

Pagel <- fitContinuous(phy, dat, SE = 0,
                    model = "lambda",
                    bounds= list(), control = list(method = "L-BFGS-B",
                                                   niter = 1000, FAIL = 1e+200, hessian = FALSE, CI = 0.95), ncores=NULL)

fitContinuous(phy, dat, SE = 0,
              model = c("BM","OU","EB","rate_trend","lambda","kappa","delta","mean_trend","white"),
              bounds= list(), control = list(method = c("subplex","L-BFGS-B"),
                                             niter = 100, FAIL = 1e+200, hessian = FALSE, CI = 0.95), ncores=NULL)