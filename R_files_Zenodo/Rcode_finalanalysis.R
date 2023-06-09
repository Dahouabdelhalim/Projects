## load library
library(MCMCglmm)

########################################################
## define names of social variables
Wsocial <- c( "indegree", 
              "outdegree", 
              "betweenness", 
              "outcloseness", 
              "incloseness", 
              "local_clustering", 
              "global_clustering", 
              "ave_shortest_path", 
              "eigenv", 
              "outstrength", 
              "instrength", 
              "outcloseness_weight", 
              "incloseness_weight", 
              "eigenv_weight"
             )
dataNames <- paste(Wsocial, rep("FWaff",length(Wsocial)),sep=".")
###############################################
## load the data
dat <- read.csv("finaldata_longevity_sociality.csv")

##Prior for MCMCglmm model
prior <- list(
  R = list (V =diag(c(1,0.00002)), nu = 1.002, fix = 2), 
  G = list( G1 = list(V = diag(2), nu = 3, alpha.mu = rep(0,2), alpha.V=diag(25^2,2)),
           G2 = list( V = 1, nu = 0.002),
           G3 = list( V = 1, nu = 0.002)
  ) 
)

##Loop running bivariate models for each SNT and longevity.
for (j in 1:length(dataNames)) {
    dat$SNT <- scale(dat[[dataNames[j]]])
    dat$loglg <- scale(log(dat$longevity))
    m1 <- MCMCglmm(c(SNT, loglg) ~ trait -1 + at.level(trait,1):log(cort_ngL) + at.level(trait,1):age + trait:valley,
               random = ~ us(trait):uid + idh(at.level(trait,2)):yrborn + idh(at.level(trait,1)):year,
               rcov = ~ idh(trait):units,
               data=dat, 
               family= c("gaussian","gaussian"),
               prior=prior,
               nitt = 2300000, thin = 200, burnin = 300000,
               pr = TRUE
               )
    summary(m1)

    assign( paste0("m_",dataNames[j]), m1, envir = .GlobalEnv )
}

##saving all the models in a R object
save( list = ls(pattern="m_"), file = "FWaff_results_scale.rda", compress = TRUE)
