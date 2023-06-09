
#### Overview ####

##
# The below gives model specifications for analyses in the paper,
#  "Benefits of cooperation in captive Damaraland mole rats".
#
# Note that:
# - IDs of groups and individuals have been anonymised within data sets,
#    so will not match across CSV files
# - Each MCMC run will give slightly different results due to how the 
#    machinery works, but should be basically the same as presented 
#    in the paper. I can provide stored versions of the models used in 
#    the paper itself upon request
# - I do not include code for diagnostics, plotting etc here but again
#    am happy to respond to any requests for such code
#
# Any correspondence can be addressed to me at houslay@gmail.com
#
# -- Tom Houslay
#
##

#### Libraries ####

library(tidyverse)

library(MCMCglmm)

#### Analysis ####

#### Cumulative helper workload ####

df_cum_work <- read_csv("cumulative_work.csv")


mcmc_cum_work <- MCMCglmm(totHelper ~ group_size_cen,
                              random =~ ProdLitterRef,
                              rcov =~ units,
                              family = "poisson",
                              nitt = 260000,
                              burnin = 10000,
                              thin = 250,
                              data = as.data.frame(df_cum_work))


#### Individual behaviour ####

df_ind_behav <- read_csv("ind_behav.csv")

#### .Set up priors ####

# 4 random effects, parameter-expanded prior

prior_RI4 <-list(R =list(V = 1, nu = 0.002),
                 G =list(G1 =list(V = 1, 
                                  nu = 1,
                                  alpha.mu = 0,
                                  alpha.V= 25^2),
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
                                   alpha.V = 25^2)))

#### .Work ####

mcmc_ind_work <- MCMCglmm(totalWork ~ inFirstYear +
                            SexM +
                            Status*group_size_cen +
                            startPM,
                          random =~ AnimalID +
                            ScanRef +
                            LitterRef +
                            Colony, 
                          rcov =~ units,
                          family = "poisson",
                          prior = prior_RI4,
                          nitt = 260000,
                          burnin = 10000,
                          thin = 250,
                          data = as.data.frame(df_ind_behav))

#### .Eating and resting ####

mcmc_ind_rest <- MCMCglmm(totalRelax ~ inFirstYear +
                                     SexM +
                                     Status*group_size_cen +
                                     startPM,
                                   random =~ AnimalID +
                                     ScanRef +
                                     LitterRef +
                                     Colony, 
                                   rcov =~ units,
                                   family = "gaussian",
                                   prior = prior_RI4,
                                   nitt = 260000,
                                   burnin = 10000,
                                   thin = 250,
                                   data = as.data.frame(df_ind_behav))


#### Multivariate behaviour and gestational weight gain ####

df_mv_behav <- read_csv("mv_behav.csv",
                        col_types = c("ccddiii"))


#### .Set up prior ####

# 4-trait model, with random regression on weight at the first random effect
prior_quad_RR_px_1  <- list(R =list(V = diag(4), nu = 3.002),
                            G =list(G1 = list(V = diag(5),
                                              nu = 5,
                                              alpha.mu = rep(0,5),
                                              alpha.V = diag(25^2,5,5)),
                                    G2 =list(V = diag(4),
                                             nu = 4,
                                             alpha.mu = rep(0,4),
                                             alpha.V = diag(25^2,4,4))))

#### .Model ####

mcmc_quad_workweight_exclittersize <- MCMCglmm(cbind(scale(Weight), 
                                                     totMother, 
                                                     totHelper,
                                                     scale(totRelax)) ~ trait +
                                                 trait:(daysIntoGestation_sc),
                                               random =~ us(trait +
                                                              at.level(trait,1):daysIntoGestation_sc):ProdLitterRef +
                                                 us(trait):AnimalID,
                                               rcov =~ us(trait):units,
                                               family = c("gaussian", 
                                                          "poisson",
                                                          "poisson",
                                                          "gaussian"),
                                               prior = prior_quad_RR_px_1,
                                               nitt = 510000,
                                               burnin = 10000,
                                               thin = 50,
                                               pr = TRUE,
                                               data = as.data.frame(df_mv_behav))


#### Litter effects ####

df_litters <- read_csv("litter_details.csv")

#### .Inter-litter interval ####

#### ..Set up prior ####

# Random regression

prior_RR_1 <-list(R =list(V = 1, nu = 0.002),
                  G =list(G1 =list(V =diag(2), 
                                   nu = 2,
                                   alpha.mu =rep(0, 2),
                                   alpha.V=diag(25^2, 2, 2))))

#### ..Model ####

##
# Note that we get rid of 2 litters
#  where the interval is >2 years
##

mcmc_ILI <- MCMCglmm(daysFromLL ~ avg_groupsize_cen + 
                       Litter_seq_cen + 
                       tunnelSize,
                     random =~ us(1 + avg_groupsize_cen):AnimalID,
                     rcov =~ units,
                     family = "poisson",
                     prior = prior_RR_1,
                     nitt = 260000,
                     burnin = 10000,
                     thin = 250,
                     data = as.data.frame(df_litters[df_litters$daysFromLL < (365*2),]))


#### .Litter size ####


#### ..Set up prior ####

# R is fixed to 1 for an ordinal model
prior_ord <- list(R = list(V = 1, fix = 1), 
                    G = list(G1 = list(V = 1, nu = 0)))

#### ..Model ####

mcmc_littersize <- MCMCglmm(Number_offspring ~ avg_groupsize_cen + 
                              Litter_seq_cen +
                              tunnelSize,
                            random =~ AnimalID,
                            rcov =~ units,
                            family = "ordinal",
                            prior = prior_ord,
                            nitt = 510000,
                            burnin = 10000,
                            thin = 50,
                            pl = TRUE,
                            data = as.data.frame(df_litters))


#### Gestational weight gain ####

df_wt <- read_csv("gest_weights.csv")

#### .Set up prior ####

# 3 random effects, the first of which
#  has quadratic random regression

prior_RR_o3 <-list(R =list(V = 1, nu = 0.002),
                   G =list(G1 =list(V =diag(3), 
                                    nu = 3,
                                    alpha.mu =rep(0, 3),
                                    alpha.V=diag(25^2, 3, 3)),
                           G2 = list(V = 1,
                                     nu = 1,
                                     alpha.mu = 0,
                                     alpha.V = 25^2),
                           G3 = list(V = 1,
                                     nu = 1,
                                     alpha.mu = 0,
                                     alpha.V = 25^2)))

#### .Model ####

mcmc_matgrowth <- MCMCglmm(Weight ~ poly(daysIntoGestation,2) * 
                                         geom_avg_groupsize_cen +
                                         tunnelSize +
                                         scale(Litter_seq, scale = FALSE),
                                       random =~ us(1 + 
                                                      daysIntoGestation_sc +
                                                      I(daysIntoGestation_sc^2)):ProdLitterRef +
                                         AnimalID +
                                         Colony,
                                       rcov =~ units,
                                       prior = prior_RR_o3,
                                       nitt = 260000,
                                       burnin = 10000,
                                       thin = 250,
                                       data = as.data.frame(df_wt))
