#### Overview ####

##
# Model specifications for the major analyses described in
#  'Contributions of genetic and non-genetic sources to variation in cooperative behaviour in a cooperative mammal'
#  by TM Houslay, JF Nielsen, TH Clutton-Brock and published in Evolution.
#
# All models below are also provided in a 'pre-compiled' form as well, 
#  along with a separate script describing the methods for determining
#  heritabilities etc.
#
# Note that the code given here provides model specs to enable reproducible analysis only,
#  and does not include methods for diagnostics, visualisation etc.
##

#### Libraries ####

library(tidyverse)
library(MCMCglmm)

select <- dplyr::select

#### Priors ####

prior_px_10 <- list(R = list(V = 1,
                             nu = 0.002),
                    G = list(G1 = list(V = 1,
                                       nu = 1,
                                       alpha.mu = 0,
                                       alpha.V = 25^2),
                             G2 = list(V = 1,
                                       nu = 1,
                                       alpha.mu = 0,
                                       alpha.V = 25^2),
                             G3 = list(V = 1,
                                       nu = 1,
                                       alpha.mu = 0,
                                       alpha.V = 25^2),
                             G4 = list(V = 1,
                                       nu = 1,
                                       alpha.mu = 0,
                                       alpha.V = 25^2),
                             G5 = list(V = 1,
                                       nu = 1,
                                       alpha.mu = 0,
                                       alpha.V = 25^2),
                             G6 = list(V = 1,
                                       nu = 1,
                                       alpha.mu = 0,
                                       alpha.V = 25^2),
                             G7 = list(V = 1,
                                       nu = 1,
                                       alpha.mu = 0,
                                       alpha.V = 25^2),
                             G8 = list(V = 1,
                                       nu = 1,
                                       alpha.mu = 0,
                                       alpha.V = 25^2),
                             G9 = list(V = 1,
                                       nu = 1,
                                       alpha.mu = 0,
                                       alpha.V = 25^2),
                             G10 = list(V = 1,
                                        nu = 1,
                                        alpha.mu = 0,
                                        alpha.V = 25^2)))


prior_px_10_AgeCats4 <- list(R = list(V = diag(3),
                                      nu = 2.002),
                             G = list(G1 = list(V = diag(3),
                                                nu = 3,
                                                alpha.mu = c(0,0,0),
                                                alpha.V = diag(3)*25^2),
                                      G2 = list(V = diag(3),
                                                nu = 3,
                                                alpha.mu = c(0,0,0),
                                                alpha.V = diag(3)*25^2),
                                      G3 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G4 = list(V = diag(3),
                                                nu = 3,
                                                alpha.mu = c(0,0,0),
                                                alpha.V = diag(3)*25^2),
                                      G5 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G6 = list(V = diag(3),
                                                nu = 3,
                                                alpha.mu = c(0,0,0),
                                                alpha.V = diag(3)*25^2),
                                      G7 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G8 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G9 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G10 = list(V = 1,
                                                 nu = 1,
                                                 alpha.mu = 0,
                                                 alpha.V = 25^2)))

prior_px_10_AgeTraj2 <- list(R = list(V = diag(3),
                                      nu = 2.002),
                             G = list(G1 = list(V = diag(2),
                                                nu = 2,
                                                alpha.mu = c(0,0),
                                                alpha.V = diag(2)*25^2),
                                      G2 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G3 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G4 = list(V = diag(2),
                                                nu = 2,
                                                alpha.mu = c(0,0),
                                                alpha.V = diag(2)*25^2),
                                      G5 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G6 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G7 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G8 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G9 = list(V = 1,
                                                nu = 1,
                                                alpha.mu = 0,
                                                alpha.V = 25^2),
                                      G10 = list(V = 1,
                                                 nu = 1,
                                                 alpha.mu = 0,
                                                 alpha.V = 25^2)))


#### Run variables ####

nitt_long <- 510000
burnin_long <- 10000
thin_long <- 50

#### Load data ####

##
# For the sake of simplicity there is a separate file for
#  each phenotypic trait. There is also a separate 'ainv' 
#  matrix for each trait, as we pruned the pedigree for each
#  prior to determining the ainv matrix for the model in order
#  to only include pedigree relationships of interest.
##

#### .Phenotypic data ####

df_feed_ped_rel <- read_csv("df_feed_ms.csv",
                            col_types = "nncnnnnnnnnnnncccccccccc")

df_BS_ped_rel <- read_csv("df_BS_ms.csv",
                          col_types = "nncnnnnnnnnnnncccccccccc")

df_guard_ped_rel <- read_csv("df_guard_ms.csv",
                             col_types = "nncnnnnnncccccccccc")

#### .Inverse additive relationship matrices ####

invA_feed2 <- readRDS("invA_feed2.rds")
invA_BS2 <- readRDS("invA_BS2.rds")
invA_GD2 <- readRDS("invA_GD2.rds")

#### MODELS ####

#### Pup feeding ####

#### .Main model ####

mcmc_pupfeed_am_pois_trial_ginv2 <- MCMCglmm(FeedCount ~ AdlibMinutesAtGroup_log_sc +
                                               AgeCat*AgeMo_cen + isDom + isM + 
                                               AgeCat:isDom +
                                               AgeCat:isM +
                                               isDom:isM +
                                               prov_LitterSize_cen *
                                               prov_GroupSize_cen +
                                               isMixedLitter +
                                               focal_coef_inbred_sc +
                                               isDrySeason +
                                               obs_LitterR_cen +
                                               seasonLitterSeq +
                                               damWeight_sc,
                                             random =~ animal +
                                               focal_LitterRef +
                                               ID_Bseason +
                                               ID +
                                               dam +
                                               mother +
                                               prov_LitterName +
                                               Group_Bseason +
                                               prov_GroupRef +
                                               Bseason,
                                             family = "poisson",
                                             prior = prior_px_10,
                                             ginverse= list(animal = invA_feed2,
                                                            dam = invA_feed2),
                                             nitt = nitt_long,
                                             burnin = burnin_long,
                                             thin = thin_long,
                                             data = as.data.frame(df_feed_ped_rel))

saveRDS(mcmc_pupfeed_am_pois_trial_ginv2, file = "mcmc_pupfeed_am_pois_trial_ginv2.rds")


#### .Allow age category changes ####

mcmc_pupfeed_am_pois_trial_ginv2_agecats <- MCMCglmm(FeedCount ~ AdlibMinutesAtGroup_log_sc +
                                                       AgeCat*AgeMo_cen + isDom + isM + 
                                                       AgeCat:isDom +
                                                       AgeCat:isM +
                                                       isDom:isM +
                                                       prov_LitterSize_cen *
                                                       prov_GroupSize_cen +
                                                       isMixedLitter +
                                                       focal_coef_inbred_sc +
                                                       isDrySeason +
                                                       obs_LitterR_cen +
                                                       seasonLitterSeq +
                                                       damWeight_sc,
                                                     random =~ us(AgeCat):animal +
                                                       us(AgeCat):focal_LitterRef +
                                                       ID_Bseason +
                                                       us(AgeCat):ID +
                                                       dam +
                                                       us(AgeCat):mother +
                                                       prov_LitterName +
                                                       Group_Bseason +
                                                       prov_GroupRef +
                                                       Bseason,
                                                     rcov =~ idh(AgeCat):units,
                                                     family = "poisson",
                                                     prior = prior_px_10_AgeCats4,
                                                     ginverse= list(animal = invA_feed2,
                                                                    dam = invA_feed2),
                                                     nitt = nitt_long,
                                                     burnin = burnin_long,
                                                     thin = thin_long,
                                                     data = as.data.frame(df_feed_ped_rel))

saveRDS(mcmc_pupfeed_am_pois_trial_ginv2_agecats, file = "mcmc_pupfeed_am_pois_trial_ginv2_agecats.rds")


#### .Allow age trajectories ####

mcmc_pupfeed_am_pois_trial_ginv2_agetraj <- MCMCglmm(FeedCount ~ AdlibMinutesAtGroup_log_sc +
                                                       AgeCat*AgeMo_cen + isDom + isM + 
                                                       AgeCat:isDom +
                                                       AgeCat:isM +
                                                       isDom:isM +
                                                       prov_LitterSize_cen *
                                                       prov_GroupSize_cen +
                                                       isMixedLitter +
                                                       focal_coef_inbred_sc +
                                                       isDrySeason +
                                                       obs_LitterR_cen +
                                                       seasonLitterSeq +
                                                       damWeight_sc,
                                                     random =~ us(1 + AgeMo_cen):animal +
                                                       focal_LitterRef +
                                                       ID_Bseason +
                                                       us(1 + AgeMo_cen):ID +
                                                       dam +
                                                       mother +
                                                       prov_LitterName +
                                                       Group_Bseason +
                                                       prov_GroupRef +
                                                       Bseason,
                                                     rcov =~ idh(AgeCat):units,
                                                     family = "poisson",
                                                     prior = prior_px_10_AgeTraj2,
                                                     ginverse= list(animal = invA_feed2,
                                                                    dam = invA_feed2),
                                                     nitt = nitt_long,
                                                     burnin = burnin_long,
                                                     thin = thin_long,
                                                     data = as.data.frame(df_feed_ped_rel))

saveRDS(mcmc_pupfeed_am_pois_trial_ginv2_agetraj, file = "mcmc_pupfeed_am_pois_trial_ginv2_agetraj.rds")


#### .Remove direct maternal effect ####

mcmc_pupfeed_am_pois_trial_ginv2_nodme <- MCMCglmm(FeedCount ~ AdlibMinutesAtGroup_log_sc +
                                                     AgeCat*AgeMo_cen + isDom + isM + 
                                                     AgeCat:isDom +
                                                     AgeCat:isM +
                                                     isDom:isM +
                                                     prov_LitterSize_cen *
                                                     prov_GroupSize_cen +
                                                     isMixedLitter +
                                                     focal_coef_inbred_sc +
                                                     isDrySeason +
                                                     obs_LitterR_cen +
                                                     seasonLitterSeq,
                                                   random =~ animal +
                                                     focal_LitterRef +
                                                     ID_Bseason +
                                                     ID +
                                                     dam +
                                                     mother +
                                                     prov_LitterName +
                                                     Group_Bseason +
                                                     prov_GroupRef +
                                                     Bseason,
                                                   family = "poisson",
                                                   prior = prior_px_10,
                                                   ginverse= list(animal = invA_feed2,
                                                                  dam = invA_feed2),
                                                   nitt = nitt_long,
                                                   burnin = burnin_long,
                                                   thin = thin_long,
                                                   data = as.data.frame(df_feed_ped_rel))

saveRDS(mcmc_pupfeed_am_pois_trial_ginv2_nodme, file = "mcmc_pupfeed_am_pois_trial_ginv2_nodme.rds")



#### Babysitting ####

#### .Main model ####

mcmc_BS_am_bin_1 <- MCMCglmm(cbind(HalfDaysBS, HalfDaysNoBS) ~ AgeCat*AgeMo_cen*isM + 
                               isDom +
                               AgeCat:isDom +
                               isDom:isM +
                               prov_LitterSize_cen *
                               prov_GroupSize_cen +
                               isMixedLitter +
                               focal_coef_inbred_sc +
                               isDrySeason +
                               obs_LitterR_cen +
                               seasonLitterSeq +
                               damWeight_sc,
                             random =~ animal +
                               focal_LitterRef +
                               ID_Bseason +
                               ID +
                               dam +
                               mother +
                               prov_LitterName +
                               Group_Bseason +
                               prov_GroupRef +
                               Bseason,
                             family = "multinomial2",
                             prior = prior_px_10,
                             ginverse= list(animal = invA_BS2,
                                            dam = invA_BS2),
                             nitt = nitt_long,
                             burnin = burnin_long,
                             thin = thin_long,
                             data = as.data.frame(df_BS_ped_rel))

saveRDS(mcmc_BS_am_bin_1, file = "mcmc_BS_am_bin_1.rds")





#### .Allow age category changes ####


mcmc_BS_am_bin_agecats <- MCMCglmm(cbind(HalfDaysBS, HalfDaysNoBS) ~ AgeCat*AgeMo_cen*isM + 
                                     isDom +
                                     AgeCat:isDom +
                                     isDom:isM +
                                     prov_LitterSize_cen *
                                     prov_GroupSize_cen +
                                     isMixedLitter +
                                     focal_coef_inbred_sc +
                                     isDrySeason +
                                     obs_LitterR_cen +
                                     seasonLitterSeq +
                                     damWeight_sc,
                                   random =~ us(AgeCat):animal +
                                     us(AgeCat):focal_LitterRef +
                                     ID_Bseason +
                                     us(AgeCat):ID +
                                     dam +
                                     us(AgeCat):mother +
                                     prov_LitterName +
                                     Group_Bseason +
                                     prov_GroupRef +
                                     Bseason,
                                   rcov =~ idh(AgeCat):units,
                                   family = "multinomial2",
                                   prior = prior_px_10_AgeCats4,
                                   ginverse= list(animal = invA_BS2,
                                                  dam = invA_BS2),
                                   nitt = nitt_long,
                                   burnin = burnin_long,
                                   thin = thin_long,
                                   data = as.data.frame(df_BS_ped_rel))

saveRDS(mcmc_BS_am_bin_agecats, file = "mcmc_BS_am_bin_agecats.rds")


#### .Allow age trajectories ####

mcmc_BS_am_bin_agetraj <- MCMCglmm(cbind(HalfDaysBS, HalfDaysNoBS) ~ AgeCat*AgeMo_cen*isM + 
                                     isDom +
                                     AgeCat:isDom +
                                     isDom:isM +
                                     prov_LitterSize_cen *
                                     prov_GroupSize_cen +
                                     isMixedLitter +
                                     focal_coef_inbred_sc +
                                     isDrySeason +
                                     obs_LitterR_cen +
                                     seasonLitterSeq +
                                     damWeight_sc,
                                   random =~ us(1 + AgeMo_cen):animal +
                                     focal_LitterRef +
                                     ID_Bseason +
                                     us(1 + AgeMo_cen):ID +
                                     dam +
                                     mother +
                                     prov_LitterName +
                                     Group_Bseason +
                                     prov_GroupRef +
                                     Bseason,
                                   rcov =~ idh(AgeCat):units,
                                   family = "multinomial2",
                                   prior = prior_px_10_AgeTraj2,
                                   ginverse= list(animal = invA_BS2,
                                                  dam = invA_BS2),
                                   nitt = nitt_long,
                                   burnin = burnin_long,
                                   thin = thin_long,
                                   data = as.data.frame(df_BS_ped_rel))

saveRDS(mcmc_BS_am_bin_agetraj, file = "mcmc_BS_am_bin_agetraj.rds")


#### .Remove direct maternal effect ####


mcmc_BS_am_bin_nodme <- MCMCglmm(cbind(HalfDaysBS, HalfDaysNoBS) ~ AgeCat*AgeMo_cen*isM + 
                                   isDom +
                                   AgeCat:isDom +
                                   isDom:isM +
                                   prov_LitterSize_cen *
                                   prov_GroupSize_cen +
                                   isMixedLitter +
                                   focal_coef_inbred_sc +
                                   isDrySeason +
                                   obs_LitterR_cen +
                                   seasonLitterSeq,
                                 random =~ animal +
                                   focal_LitterRef +
                                   ID_Bseason +
                                   ID +
                                   dam +
                                   mother +
                                   prov_LitterName +
                                   Group_Bseason +
                                   prov_GroupRef +
                                   Bseason,
                                 family = "multinomial2",
                                 prior = prior_px_10,
                                 ginverse= list(animal = invA_BS2,
                                                dam = invA_BS2),
                                 nitt = nitt_long,
                                 burnin = burnin_long,
                                 thin = thin_long,
                                 data = as.data.frame(df_BS_ped_rel))

saveRDS(mcmc_BS_am_bin_nodme, file = "mcmc_BS_am_bin_nodme.rds")



#### Sentinel ####

#### .Main model ####

mcmc_GD_am_pois_LONG <- MCMCglmm(GuardDuration_mins ~ AdlibMinutesAtGroup_log_sc +
                                   AgeCat*AgeMo_cen*isM + isDom +
                                   AgeCat:isDom +
                                   isDom:isM +
                                   GroupSize_cen +
                                   damWeight_sc +
                                   focal_coef_inbred_sc,
                                 random =~ animal +
                                   focal_LitterRef +
                                   ID_Bseason +
                                   ID +
                                   dam +
                                   mother +
                                   Group_Bseason +
                                   GroupRef +
                                   Quarter_fac +
                                   Bseason,
                                 family = "poisson",
                                 prior = prior_px_10,
                                 ginverse= list(animal = invA_GD2,
                                                dam = invA_GD2),
                                 nitt = nitt_long,
                                 burnin = burnin_long,
                                 thin = thin_long,
                                 data = as.data.frame(df_guard_ped_rel))


saveRDS(mcmc_GD_am_pois_LONG, file = "mcmc_GD_am_pois_LONG.RDS")





#### .Allow age category changes ####


mcmc_GD_am_pois_LONG_agecats <- MCMCglmm(GuardDuration_mins ~ AdlibMinutesAtGroup_log_sc +
                                           AgeCat*AgeMo_cen*isM + isDom +
                                           AgeCat:isDom +
                                           isDom:isM +
                                           GroupSize_cen +
                                           damWeight_sc +
                                           focal_coef_inbred_sc,
                                         random =~ us(AgeCat):animal +
                                           us(AgeCat):focal_LitterRef +
                                           ID_Bseason +
                                           us(AgeCat):ID +
                                           dam +
                                           us(AgeCat):mother +
                                           Group_Bseason +
                                           GroupRef +
                                           Quarter_fac +
                                           Bseason,
                                         rcov =~ idh(AgeCat):units,
                                         family = "poisson",
                                         prior = prior_px_10_AgeCats4,
                                         ginverse= list(animal = invA_GD2,
                                                        dam = invA_GD2),
                                         nitt = nitt_long,
                                         burnin = burnin_long,
                                         thin = thin_long,
                                         data = as.data.frame(df_guard_ped_rel))

saveRDS(mcmc_GD_am_pois_LONG_agecats, file = "mcmc_GD_am_pois_LONG_agecats.rds")


#### .Allow age trajectories ####


mcmc_GD_am_pois_LONG_agetraj <- MCMCglmm(GuardDuration_mins ~ AdlibMinutesAtGroup_log_sc +
                                           AgeCat*AgeMo_cen*isM + isDom +
                                           AgeCat:isDom +
                                           isDom:isM +
                                           GroupSize_cen +
                                           damWeight_sc +
                                           focal_coef_inbred_sc,
                                         random =~ us(1 + AgeMo_cen):animal +
                                           focal_LitterRef +
                                           ID_Bseason +
                                           us(1 + AgeMo_cen):ID +
                                           dam +
                                           mother +
                                           Group_Bseason +
                                           GroupRef +
                                           Quarter_fac +
                                           Bseason,
                                         rcov =~ idh(AgeCat):units,
                                         family = "poisson",
                                         prior = prior_px_10_AgeTraj2,
                                         ginverse= list(animal = invA_GD2,
                                                        dam = invA_GD2),
                                         nitt = nitt_long,
                                         burnin = burnin_long,
                                         thin = thin_long,
                                         data = as.data.frame(df_guard_ped_rel))

saveRDS(mcmc_GD_am_pois_LONG_agetraj, file = "mcmc_GD_am_pois_LONG_agetraj.rds")


#### .Remove direct maternal effect ####


mcmc_GD_am_pois_LONG_nodme <- MCMCglmm(GuardDuration_mins ~ AdlibMinutesAtGroup_log_sc +
                                         AgeCat*AgeMo_cen*isM + isDom +
                                         AgeCat:isDom +
                                         isDom:isM +
                                         GroupSize_cen +
                                         #damWeight_sc +
                                         focal_coef_inbred_sc,
                                       random =~ animal +
                                         focal_LitterRef +
                                         ID_Bseason +
                                         ID +
                                         dam +
                                         mother +
                                         Group_Bseason +
                                         GroupRef +
                                         Quarter_fac +
                                         Bseason,
                                       family = "poisson",
                                       prior = prior_px_10,
                                       ginverse= list(animal = invA_GD2,
                                                      dam = invA_GD2),
                                       nitt = nitt_long,
                                       burnin = burnin_long,
                                       thin = thin_long,
                                       data = as.data.frame(df_guard_ped_rel))


saveRDS(mcmc_GD_am_pois_LONG_nodme, file = "mcmc_GD_am_pois_LONG_nodme.rds")
