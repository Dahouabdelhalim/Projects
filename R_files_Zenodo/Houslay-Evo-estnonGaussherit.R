
#### Overview ####

##
# Code for estimating heritabilities from non-Gaussian MCMCglmm models as shown in
#  'Contributions of genetic and non-genetic sources to variation in cooperative behaviour in a cooperative mammal'
#  by TM Houslay, JF Nielsen, TH Clutton-Brock and published in Evolution.
#
#  The below uses QGglmm (de Villemereuil et al 2016) with MCMCglmm model runs.
#
# Pre-compiled models have been provided, as well as the code to run these again
#  (note that running models again will cause estimates / intervals to be slightly
#  different than shown here because of the MCMC sampling procedure).
##


#### Libraries ####

library(tidyverse)
library(MCMCglmm)
library(QGglmm)


#### Babysitting ####

mcmc_BS_am_bin_1 <- readRDS("Models/mcmc_BS_am_bin_1.rds")

summary(mcmc_BS_am_bin_1)

plot((mcmc_BS_am_bin_1$VCV[,"prov_LitterName"] +
        mcmc_BS_am_bin_1$VCV[,"prov_GroupRef"] +
        mcmc_BS_am_bin_1$VCV[,"Group_Bseason"])/rowSums(mcmc_BS_am_bin_1$VCV))


plot((mcmc_BS_am_bin_1$VCV[,"units"])/rowSums(mcmc_BS_am_bin_1$VCV))


#### .h2 no FE ####

df_BS_heritnoFE <- data.frame(mu = as.vector(mcmc_BS_am_bin_1[["Sol"]][, "(Intercept)"]),
                              va = as.vector(mcmc_BS_am_bin_1[["VCV"]][, "animal"]),
                              vp = rowSums(mcmc_BS_am_bin_1[["VCV"]]))

df_BS_heritnoFE_post <- df_BS_heritnoFE$va/(df_BS_heritnoFE$vp + 1)

mean(df_BS_heritnoFE_post)
HPDinterval(as.mcmc(df_BS_heritnoFE_post))

#### .Fixed effects variance ####

# Turn this into a function
yhatdist_BS_am <- sapply(1:nrow(mcmc_BS_am_bin_1[["Sol"]]), function(i) {
  var(predict(mcmc_BS_am_bin_1, 
              type = "terms",
              it = i))			#Using Eq. 5 in the article
})


df_BS_herit <- data.frame(pr = yhatdist_BS_am,
                          va = as.vector(mcmc_BS_am_bin_1[["VCV"]][, "animal"]),
                          vp = rowSums(mcmc_BS_am_bin_1[["VCV"]]))

head(df_BS_herit)


# 
# df_BS_herit_post <- do.call("rbind", 
#                             apply(as.data.frame(df_BS_herit), 1, function(row){
#                               QGparams(predict = row[["pr"]],
#                                        var.a = row[["va"]],
#                                        var.p = row[["vp"]],
#                                        model = "Poisson.log", 
#                                        verbose = FALSE)
#                             }))

plot(as.mcmc(df_BS_herit$va/(df_BS_herit$va + df_BS_herit$pr + df_BS_herit$vp + 1)))
mean(as.mcmc(df_BS_herit$va/(df_BS_herit$va + df_BS_herit$pr + df_BS_herit$vp + 1)))
HPDinterval(as.mcmc(df_BS_herit$va/(df_BS_herit$va + df_BS_herit$pr + df_BS_herit$vp + 1)))



df_BS_replong <- data.frame(pr = yhatdist_BS_am,
                            va = as.vector(mcmc_BS_am_bin_1[["VCV"]][, "animal"] +
                                             mcmc_BS_am_bin_1[["VCV"]][, "focal_LitterRef"] +
                                             mcmc_BS_am_bin_1[["VCV"]][, "ID"] +
                                             mcmc_BS_am_bin_1[["VCV"]][, "dam"] +
                                             mcmc_BS_am_bin_1[["VCV"]][, "mother"]),
                            vp = rowSums(mcmc_BS_am_bin_1[["VCV"]]))

plot(as.mcmc(df_BS_replong$va/(df_BS_replong$va + df_BS_replong$pr + df_BS_replong$vp + 1)))
mean(as.mcmc(df_BS_replong$va/(df_BS_replong$va + df_BS_replong$pr + df_BS_replong$vp + 1)))
HPDinterval(as.mcmc(df_BS_replong$va/(df_BS_replong$va + df_BS_replong$pr + df_BS_replong$vp + 1)))



df_BS_repshort <- data.frame(pr = yhatdist_BS_am,
                             va = as.vector(mcmc_BS_am_bin_1[["VCV"]][, "animal"] +
                                              mcmc_BS_am_bin_1[["VCV"]][, "focal_LitterRef"] +
                                              mcmc_BS_am_bin_1[["VCV"]][, "ID"] +
                                              mcmc_BS_am_bin_1[["VCV"]][, "ID_Bseason"] +
                                              mcmc_BS_am_bin_1[["VCV"]][, "dam"] +
                                              mcmc_BS_am_bin_1[["VCV"]][, "mother"]),
                             vp = rowSums(mcmc_BS_am_bin_1[["VCV"]]))

plot(as.mcmc(df_BS_repshort$va/(df_BS_repshort$va + df_BS_repshort$pr + df_BS_repshort$vp + 1)))
mean(as.mcmc(df_BS_repshort$va/(df_BS_repshort$va + df_BS_repshort$pr + df_BS_repshort$vp + 1)))
HPDinterval(as.mcmc(df_BS_repshort$va/(df_BS_repshort$va + df_BS_repshort$pr + df_BS_repshort$vp + 1)))



#### Sentinel ####

mcmc_GD_am_pois_LONG <- readRDS("Models/mcmc_GD_am_pois_LONG.RDS")

summary(mcmc_GD_am_pois_LONG)

plot((mcmc_GD_am_pois_LONG$VCV[,"Group_Bseason"] + 
        mcmc_GD_am_pois_LONG$VCV[,"GroupRef"])/
       rowSums(mcmc_GD_am_pois_LONG$VCV))

plot((mcmc_GD_am_pois_LONG$VCV[,"units"])/
       rowSums(mcmc_GD_am_pois_LONG$VCV))


mean(mcmc_GD_am_pois_LONG$VCV[,"animal"])
HPDinterval(mcmc_GD_am_pois_LONG$VCV[,"animal"])

mean(as.mcmc(rowSums(mcmc_GD_am_pois_LONG[["VCV"]])))
HPDinterval(as.mcmc(rowSums(mcmc_GD_am_pois_LONG[["VCV"]])))

#### .Heritability No FE ####


df_GD_heritNOFE <- data.frame(mu = as.vector(mcmc_GD_am_pois_LONG[["Sol"]][, "(Intercept)"]),
                              va = as.vector(mcmc_GD_am_pois_LONG[["VCV"]][, "animal"]),
                              vp = rowSums(mcmc_GD_am_pois_LONG[["VCV"]]))

head(df_GD_heritNOFE)

df_GD_heritNOFE_post <- do.call("rbind", 
                                apply(df_GD_heritNOFE, 1, function(row){
                                  QGparams(mu = row[["mu"]],
                                           var.a = row[["va"]],
                                           var.p = row[["vp"]],
                                           model = "Poisson.log", 
                                           verbose = FALSE)
                                }))

plot(as.mcmc(df_GD_heritNOFE_post$h2.obs))
HPDinterval(as.mcmc(df_GD_heritNOFE_post$h2.obs), 0.95)

mean(as.mcmc(df_GD_heritNOFE_post$h2.obs))

#### .Fixed effects variance ####

# Turn this into a function
yhatdist_GD_pois <- sapply(1:nrow(mcmc_GD_am_pois_LONG[["Sol"]]), function(i) {
  var(predict(mcmc_GD_am_pois_LONG, 
              type = "terms",
              it = i))			#Using Eq. 5 in the article
})

df_GD_herit <- data.frame(pr = yhatdist_GD_pois,
                          va = as.vector(mcmc_GD_am_pois_LONG[["VCV"]][, "animal"]),
                          vp = rowSums(mcmc_GD_am_pois_LONG[["VCV"]]))

head(df_GD_herit)



df_GD_herit_post <- do.call("rbind", 
                            apply(as.data.frame(df_GD_herit), 1, function(row){
                              QGparams(predict = row[["pr"]],
                                       var.a = row[["va"]],
                                       var.p = row[["vp"]],
                                       model = "Poisson.log", 
                                       verbose = FALSE)
                            }))

plot(as.mcmc(df_GD_herit_post$h2.obs))

mean(as.mcmc(df_GD_herit_post$h2.obs))
posterior.mode(as.mcmc(df_GD_herit_post$h2.obs))
HPDinterval(as.mcmc(df_GD_herit_post$h2.obs))



#### ..Repeatability ####


df_GD_replong <- data.frame(pr = yhatdist_GD_pois,
                            va = as.vector(mcmc_GD_am_pois_LONG[["VCV"]][, "animal"] +
                                             mcmc_GD_am_pois_LONG[["VCV"]][, "focal_LitterRef"] +
                                             mcmc_GD_am_pois_LONG[["VCV"]][, "ID"] +
                                             mcmc_GD_am_pois_LONG[["VCV"]][, "dam"] +
                                             mcmc_GD_am_pois_LONG[["VCV"]][, "mother"]),
                            vp = rowSums(mcmc_GD_am_pois_LONG[["VCV"]]))

head(df_GD_replong)



df_GD_replong_post <- do.call("rbind", 
                              apply(as.data.frame(df_GD_replong), 1, function(row){
                                QGparams(predict = row[["pr"]],
                                         var.a = row[["va"]],
                                         var.p = row[["vp"]],
                                         model = "Poisson.log", 
                                         verbose = FALSE)
                              }))


plot(as.mcmc(df_GD_replong_post$h2.obs))

mean(as.mcmc(df_GD_replong_post$h2.obs))
posterior.mode(as.mcmc(df_GD_replong_post$h2.obs))
HPDinterval(as.mcmc(df_GD_replong_post$h2.obs))




df_GD_repshort <- data.frame(pr = yhatdist_pupfeed_am_pois_trial_ginv2,
                             va = as.vector(mcmc_GD_am_pois_LONG[["VCV"]][, "animal"] +
                                              mcmc_GD_am_pois_LONG[["VCV"]][, "focal_LitterRef"] +
                                              mcmc_GD_am_pois_LONG[["VCV"]][, "ID"] +
                                              mcmc_GD_am_pois_LONG[["VCV"]][, "ID_Bseason"] +
                                              mcmc_GD_am_pois_LONG[["VCV"]][, "dam"] +
                                              mcmc_GD_am_pois_LONG[["VCV"]][, "mother"]),
                             vp = rowSums(mcmc_GD_am_pois_LONG[["VCV"]]))

head(df_GD_repshort)



df_GD_repshort_post <- do.call("rbind", 
                               apply(as.data.frame(df_GD_repshort), 1, function(row){
                                 QGparams(predict = row[["pr"]],
                                          var.a = row[["va"]],
                                          var.p = row[["vp"]],
                                          model = "Poisson.log", 
                                          verbose = FALSE)
                               }))


plot(as.mcmc(df_GD_repshort_post$h2.obs))

mean(as.mcmc(df_GD_repshort_post$h2.obs))
posterior.mode(as.mcmc(df_GD_repshort_post$h2.obs))
HPDinterval(as.mcmc(df_GD_repshort_post$h2.obs))





#### Pup feeding ####


mcmc_pupfeed_am_pois_trial_ginv2 <- readRDS("Models/mcmc_pupfeed_am_pois_trial_ginv2.rds")

summary(mcmc_pupfeed_am_pois_trial_ginv2)

mean(as.mcmc(rowSums(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]])))
HPDinterval(as.mcmc(rowSums(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]])))

plot((mcmc_pupfeed_am_pois_trial_ginv2$VCV[,"prov_LitterName"])/
       rowSums(mcmc_pupfeed_am_pois_trial_ginv2$VCV))

plot((mcmc_pupfeed_am_pois_trial_ginv2$VCV[,"prov_LitterName"] +
        mcmc_pupfeed_am_pois_trial_ginv2$VCV[,"Group_Bseason"] + 
        mcmc_pupfeed_am_pois_trial_ginv2$VCV[,"prov_GroupRef"])/
       rowSums(mcmc_pupfeed_am_pois_trial_ginv2$VCV))


mcmc_pupfeed_am_pois_trial_ginv2_nodme <- readRDS("Models/mcmc_pupfeed_am_pois_trial_ginv2_nodme.rds")

summary(mcmc_pupfeed_am_pois_trial_ginv2_nodme)

plot(mcmc_pupfeed_am_pois_trial_ginv2_nodme$VCV)


#### .h2 no FE ####


df_pupfeed_heritNOFE <- data.frame(mu = as.vector(mcmc_pupfeed_am_pois_trial_ginv2[["Sol"]][, "(Intercept)"]),
                                   va = as.vector(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "animal"]),
                                   vp = rowSums(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]]))

head(df_pupfeed_heritNOFE)

df_pupfeed_heritNOFE_post <- do.call("rbind", 
                                     apply(df_pupfeed_heritNOFE, 1, function(row){
                                       QGparams(mu = row[["mu"]],
                                                var.a = row[["va"]],
                                                var.p = row[["vp"]],
                                                model = "Poisson.log", 
                                                verbose = FALSE)
                                     }))

plot(as.mcmc(df_pupfeed_heritNOFE_post$h2.obs))
HPDinterval(as.mcmc(df_pupfeed_heritNOFE_post$h2.obs), 0.95)

mean(as.mcmc(df_pupfeed_heritNOFE_post$h2.obs))

#### .Fixed effects variance ####

# Turn this into a function
yhatdist_pupfeed_am_pois_trial_ginv2 <- sapply(1:nrow(mcmc_pupfeed_am_pois_trial_ginv2[["Sol"]]), function(i) {
  var(predict(mcmc_pupfeed_am_pois_trial_ginv2, 
              type = "terms",
              it = i))			#Using Eq. 5 in the article
})


mu_pupfeed_am_trial <- posterior.mode(mcmc_pupfeed_am_pois_trial_ginv2[["Sol"]][,"(Intercept)"])
va_pupfeed_am_trial <- posterior.mode(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][,"animal"])
vp_pupfeed_am_trial <- posterior.mode(as.mcmc(rowSums(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]])))

QGparams(mu = mu_pupfeed_am_trial, 
         var.a = va_pupfeed_am_trial, 
         var.p = vp_pupfeed_am_trial,
         model = "Poisson.log")

QGparams(predict = posterior.mode(as.mcmc(yhatdist_pupfeed_am_pois_trial_ginv2)), 
         var.a = va_pupfeed_am_trial, 
         var.p = vp_pupfeed_am_trial,
         model = "Poisson.log")

df_pupfeed_trial_herit <- data.frame(pr = yhatdist_pupfeed_am_pois_trial_ginv2,
                                     va = as.vector(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "animal"]),
                                     vp = rowSums(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]]))

head(df_pupfeed_trial_herit)



df_pupfeed_trial_herit_post <- do.call("rbind", 
                                       apply(as.data.frame(df_pupfeed_trial_herit), 1, function(row){
                                         QGparams(predict = row[["pr"]],
                                                  var.a = row[["va"]],
                                                  var.p = row[["vp"]],
                                                  model = "Poisson.log", 
                                                  verbose = FALSE)
                                       }))

plot(as.mcmc(df_pupfeed_trial_herit_post$h2.obs))

mean(as.mcmc(df_pupfeed_trial_herit_post$h2.obs))
posterior.mode(as.mcmc(df_pupfeed_trial_herit_post$h2.obs))
HPDinterval(as.mcmc(df_pupfeed_trial_herit_post$h2.obs))

plot(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][,"animal"]/
       as.mcmc(rowSums(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]])))

mean(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][,"animal"]/
       as.mcmc(rowSums(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]])))

posterior.mode(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][,"animal"]/
                 as.mcmc(rowSums(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]])))


#### ..Repeatability ####


df_pupfeed_trial_replong <- data.frame(pr = yhatdist_pupfeed_am_pois_trial_ginv2,
                                       va = as.vector(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "animal"] +
                                                        mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "focal_LitterRef"] +
                                                        mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "ID"] +
                                                        mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "dam"] +
                                                        mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "mother"]),
                                       vp = rowSums(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]]))

head(df_pupfeed_trial_replong)



df_pupfeed_trial_replong_post <- do.call("rbind", 
                                         apply(as.data.frame(df_pupfeed_trial_replong), 1, function(row){
                                           QGparams(predict = row[["pr"]],
                                                    var.a = row[["va"]],
                                                    var.p = row[["vp"]],
                                                    model = "Poisson.log", 
                                                    verbose = FALSE)
                                         }))


plot(as.mcmc(df_pupfeed_trial_replong_post$h2.obs))

mean(as.mcmc(df_pupfeed_trial_replong_post$h2.obs))
posterior.mode(as.mcmc(df_pupfeed_trial_replong_post$h2.obs))
HPDinterval(as.mcmc(df_pupfeed_trial_replong_post$h2.obs))




df_pupfeed_trial_repshort <- data.frame(pr = yhatdist_pupfeed_am_pois_trial_ginv2,
                                        va = as.vector(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "animal"] +
                                                         mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "focal_LitterRef"] +
                                                         mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "ID"] +
                                                         mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "ID_Bseason"] +
                                                         mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "dam"] +
                                                         mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]][, "mother"]),
                                        vp = rowSums(mcmc_pupfeed_am_pois_trial_ginv2[["VCV"]]))

head(df_pupfeed_trial_repshort)



df_pupfeed_trial_repshort_post <- do.call("rbind", 
                                          apply(as.data.frame(df_pupfeed_trial_repshort), 1, function(row){
                                            QGparams(predict = row[["pr"]],
                                                     var.a = row[["va"]],
                                                     var.p = row[["vp"]],
                                                     model = "Poisson.log", 
                                                     verbose = FALSE)
                                          }))


plot(as.mcmc(df_pupfeed_trial_repshort_post$h2.obs))

mean(as.mcmc(df_pupfeed_trial_repshort_post$h2.obs))
posterior.mode(as.mcmc(df_pupfeed_trial_repshort_post$h2.obs))
HPDinterval(as.mcmc(df_pupfeed_trial_repshort_post$h2.obs))

