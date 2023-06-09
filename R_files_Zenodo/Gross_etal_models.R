#A priori models of trait dispersion
##Note: a priori models are only shown for SES of mean pairwise distance for all traits based on the global pool.
###Model selection steps, including predictors included in intermediate composite models, are shown in a separate Excel file. 
library(MuMIn)

biogeo.1 <- lm(ses.traits.is.mpd.mpd$mpd.obs.z ~ 
                                Latitude,
                              data = site_data)
summary(biogeo.1)
AICc(biogeo.1) 

biogeo.2 <- lm(ses.traits.is.mpd.mpd$mpd.obs.z ~ 
                                Latitude + 
                                Basin,
                              data = site_data)
summary(biogeo.2)
AICc(biogeo.2)

biogeo.3 <- lm(ses.traits.is.mpd.mpd$mpd.obs.z ~ 
                                Latitude *
                                Basin,
                              data = site_data)
summary(biogeo.3)
AICc(ses.traits.mpd.biogeo.3)

biogeo.4 <- lm(ses.traits.is.mpd.mpd$mpd.obs.z ~ 
                                Latitude *
                                Basin + 
                                Ocean,
                              data = site_data)
summary(biogeo.4)
AICc(biogeo.4) 

biogeo.5 <- lm(ses.traits.is.mpd.mpd$mpd.obs.z ~ 
                                Latitude *
                                Basin +
                                Latitude *
                                Ocean,
                              data = site_data)
summary(biogeo.5)
AICc(biogeo.5) 

abiotic <- lm(ses.traits.is.mpd.mpd$mpd.obs.z ~ 
                               Temperature.C +
                               Salinity.ppt + 
                               Mean.Leaf.PercN,
                             data = site_data)
summary(abiotic)
AICc(ses.traits.mpd.abiotic) 

temp.1 <- lm(ses.traits.is.mpd.mpd$mpd.obs.z ~ 
                              sstmean,
                            data = site_data)
summary(temp.1)
AICc(temp.1) 

temp.2 <- lm(ses.traits.is.mpd.mpd$mpd.obs.z ~ 
                              sstrange,
                            data = site_data)
summary(temp.2)
AICc(temp.2) 

temp.3 <- lm(ses.traits.is.mpd.mpd$mpd.obs.z ~ 
                              sstmean *
                              sstrange,
                            data = site_data)
summary(temp.3)
AICc(temp.3) 

community <- lm(ses.traits.is.mpd.mpd$mpd.obs.z ~
                                 log(Mean.Std.Total.Abund.Crustaceans) +
                                 crustaceans.med.size,
                               data = site_data)
summary(community)
AICc(community) 

total.biodiv <- lm(ses.traits.is.mpd$mpd.obs.z ~
                                    log(Site.Epifaunal.Richness),
                                  data = site_data)
summary(total.biodiv)
AICc(total.biodiv) 

biodiv <- lm(ses.traits.is.mpd$mpd.obs.z ~
                              log(Site.Peracarid.Richness),
                            data = site_data)
summary(biodiv)
AICc(ses.traits.mpd.biodiv) 

habitat <- lm(ses.traits.is.mpd$mpd.obs.z ~
                                PC1 + 
                                PC2 + 
                                log(Macroalgae.g + 1),
                              data = site_data)
summary(habitat)
AICc(habitat) 

pred <- lm(ses.traits.is.mpd$mpd.obs.z ~
                            asin(Mean.Pred.Amphipod),
                          data = site_data)
summary(pred)
AICc(pred) 

resource.1 <- lm(ses.traits.is.mpd$mpd.obs.z ~
                                  log(Mean.Site.Std.Periphyton) +
                                  log(chlomean),
                                data = site_data)
summary(resource.1)
AICc(resource.1) 

resource.2 <- lm(ses.traits.is.mpd$mpd.obs.z ~
                                  sqrt(nitrate) +
                                  parmean,
                                data = site_data)
summary(resource.2)
AICc(resource.2) 