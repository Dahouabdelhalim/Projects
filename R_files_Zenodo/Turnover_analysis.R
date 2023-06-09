## Models to test for community turnover 

library(MCMCglmm)

## Read in full data set
full_data_long <- read.csv("OTU_matrix_by_spp_long.csv")

## Read in terpene data set
terpene_data_long <- read.csv("Terpene_OTU_long.csv")

## Set K
k<-1000

## Full dataset model with maternal family
FamilyMod <- MCMCglmm(fixed = Abundance ~ log(Tot_tips + 0.000001), 
                  random = ~OTU + Site:OTU + Family:OTU + Grid:OTU,
                  family = "poisson",
                  data = full_data_long, 
                  nitt = 2000000, burnin = 40000,
                  prior = list(B=list(V=matrix(c(1e7, 0, 0, 1e-7), ncol = 2), mu = c(0, 1)),
                               G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k),
                                      G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k),
                                      G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k),
                                      G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k))), 
                  thin = 200,
                  pr = F)

summary(FamilyMod)

## Need more complex fixed prior to account for extra fixed effects
PriorM <- diag(3)*1e7
## Set variance for log(Tot_tips) to be very small to fix parameter at 1 for offset
PriorM[2,2] <- 1e-7

## Terpene dataset model with chemodiversity
TerpModel <- MCMCglmm(fixed = Abundance ~ log(Tot_tips + 0.0000001) + Chemodiversity, 
                             random = ~idh(1 + Chemodiversity):OTU + Site:OTU + Grid:OTU,
                             family = "poisson",
                             data = terpene_data_long, 
                             nitt = 2000000, burnin = 40000,
                             prior = list(B=list(V=PriorM, mu = c(0, 1, 0)),
                                          G=list(G1=list(V=diag(2),nu=1,aplha.mu=0,alpha.V=diag(2)*k),
                                                 G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k),
                                                 G1=list(V=1,nu=1,alpha.mu=0,alpha.V=k))),
                             thin = 200,
                             pr = F)

summary(TerpModel)
