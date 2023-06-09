
library (glmmTMB)
library (ggplot2)
library (performance)
library (lattice)
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
library (ggpubr)
library (cowplot)
library(grid)
library(gridExtra)
library (dplyr)
library (ggrepel)
library (car)
library (MuMIn)
library (arm)

# import data from your source folder (where you have downloaded it into)

data <- read.csv("/Users/KrishnaBalasubramaniam/Box Sync/McCowanLab_Krishna/Projects_CNHS/CNHS_projects/Ms_KNB_Networks&Disease_Anth&SocialNetworks/Anth&SocialNetworks_Rev1/Anth&SocialNetworks_R2/JAE-2021-00351_DatasetforGLMMs.csv", na.strings = NA)
attach(data)

# ignore outlier individual with low observation time (ID: )

# Plot outcome variables

hist(HumanCointeractionStrengthPercentile)
hist(HumanCointeractionEigenvectorPercentile)


*********************Analysis for Strength Centrality*************************
  
# Models for Strength
  
mod.Str.nonetworkmeasures <- glmmTMB(HumanCointeractionStrengthPercentile ~ Sex + RankPercentile + Species + PropScansHumanProximityPercentile + (1|Group), data=data, family = beta_family())
mod.Str.groomingandaffiliation <- glmmTMB(HumanCointeractionStrengthPercentile ~ Sex + RankPercentile + Species + GroomStrengthPercentile + PropScansHumanProximityPercentile + OtherAffiliationStrengthPercentile + (1|Group), data=data, family = beta_family())
mod.Str.proximityandaffiliation <- glmmTMB(HumanCointeractionStrengthPercentile ~ Sex + RankPercentile + Species + ProximityStrengthPercentile + PropScansHumanProximityPercentile + OtherAffiliationStrengthPercentile + (1|Group), data=data, family = beta_family())
mod.Str.groomingandaffiliation.groomingbyspecies <- glmmTMB(HumanCointeractionStrengthPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + GroomStrengthPercentile*Species + OtherAffiliationStrengthPercentile + (1|Group), data=data, family = beta_family())
mod.Str.groomingandaffiliation.affiliationbyspecies <- glmmTMB(HumanCointeractionStrengthPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + GroomStrengthPercentile + OtherAffiliationStrengthPercentile*Species + (1|Group), data=data, family = beta_family())
mod.Str.proximityandaffiliation.proximitybyspecies <- glmmTMB(HumanCointeractionStrengthPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + ProximityStrengthPercentile*Species + OtherAffiliationStrengthPercentile + (1|Group), data=data, family = beta_family())
mod.Str.proximityandaffiliation.affiliationbyspecies <- glmmTMB(HumanCointeractionStrengthPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + ProximityStrengthPercentile + OtherAffiliationStrengthPercentile*Species + (1|Group), data=data, family = beta_family())

##Get the AICc scores to determine the best-fit model

AICc(mod.Str.nonetworkmeasures)
AICc(mod.Str.groomingandaffiliation)
AICc(mod.Str.proximityandaffiliation)
AICc(mod.Str.groomingandaffiliation.groomingbyspecies)
AICc(mod.Str.groomingandaffiliation.affiliationbyspecies)
AICc(mod.Str.proximityandaffiliation.proximitybyspecies)
##This model should have the lowest AICc score
AICc(mod.Str.proximityandaffiliation.affiliationbyspecies)

##Get model summary statistics

summary(mod.Str.nonetworkmeasures)
summary(mod.Str.groomingandaffiliation)
summary(mod.Str.proximityandaffiliation)
summary(mod.Str.groomingandaffiliation.groomingbyspecies)
summary(mod.Str.groomingandaffiliation.affiliationbyspecies)
summary(mod.Str.proximityandaffiliation.proximitybyspecies)
summary(mod.Str.proximityandaffiliation.affiliationbyspecies)


##check validity of best-fit####

resideStr <- resid (mod.Str.proximityandaffiliation.proximitybyspecies)
hist (resideStr)
diagnostics.plot(mod.Str.proximityandaffiliation.proximitybyspecies)

##Get model-predicted values using the "predict" function 

data$ModPredStr <- predict(mod.Str.proximityandaffiliation.proximitybyspecies, type="response")

##use relevel function to get interaction coefficients for each combination of species ####

data$Species <- relevel(data$Species, ref = "Rhesus")
mod.Str.proximityandaffiliation.proximitybyspeciesRhesus <- glmmTMB(HumanCointeractionStrengthPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + ProximityStrengthPercentile*Species + OtherAffiliationStrengthPercentile + (1|Group), data=data, family = beta_family())
summary(mod.Str.proximityandaffiliation.proximitybyspeciesRhesus)

data$Species <- relevel(data$Species, ref = "Long-tailed")
mod.Str.proximityandaffiliation.proximitybyspeciesLTM <- glmmTMB(HumanCointeractionStrengthPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + ProximityStrengthPercentile*Species + OtherAffiliationStrengthPercentile + (1|Group), data=data, family = beta_family())
summary(mod.Str.proximityandaffiliation.proximitybyspeciesLTM)


#Get model-predicted values 

data$ModPredStr <- predict(mod.Str.proximityandaffiliation.proximitybyspecies, type="response")


###Plot the significant effects

#Plot the effect of affiliation Strength on human co-interaction Strength

plotaffiliationStr <- ggplot(data, aes(x = OtherAffiliationStrengthPercentile, y = HumanCointeractionStrengthPercentile)) +
  geom_point(size = 1, alpha = 0.9) +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network Strength", x = "Short-duration affiliation network Strength", tag = "a") +
  theme_cowplot() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotaffiliationStr


#Plot the effect of affiliation Strength on model-predicted Human co-interaction Strength

plotaffiliationStrmodpred <- ggplot(data, aes(x = OtherAffiliationStrengthPercentile, y = ModPredStr)) +
  geom_point(size = 1, alpha = 0.9) +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network Strength", x = "Short-duration affiliation network Strength", tag = "a") +
  theme_cowplot() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotaffiliationStrmodpred


#Plot the effect of proximity Strength by species on human co-interaction Strength

plotproximityStrbyspecies <- ggplot(data, aes(x = ProximityStrengthPercentile, y = HumanCointeractionStrengthPercentile, color = Species)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network Strength", x = "Proximity network Strength", tag = "b") +
  theme_cowplot() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotproximityStrbyspecies


#Plot the effect of proximity Strength by species on model-predicted human co-interaction Strength

plotproximityStrbyspeciesmodpred <- ggplot(data, aes(x = ProximityStrengthPercentile, y = ModPredStr, color = Species)) +
  scale_color_manual(values=c("grey1","grey60", "blue")) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network Strength", x = "Proximity network Strength", tag = "b") +
  theme_cowplot() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotproximityStrbyspeciesmodpred


#Plot the effect of sex on human co-interaction network Strength

plotSexStr <- ggplot(data, aes(x = Sex, y = HumanCointeractionStrengthPercentile)) + 
  geom_boxplot() +
  geom_point(size = 1, alpha = 0.9) +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network Strength", x = "Macaque sex") +
  theme_cowplot() +
  labs(title = "Tau = 0.10") +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotSexStr


#Plot the effect of sex on model-predicted human co-interaction network Strength

plotSexStrModPred <- ggplot(data, aes(x = Sex, y = ModPredStr)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  geom_point(size = 1, alpha = 0.9) +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network Strength", x = "Macaque sex", tag = "c") +
  theme_cowplot() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotSexStrModPred




*********************Analysis for Eigenvector Centrality*************************

# Models for Eigenvector

mod.evc.nonetworkmeasures <- glmmTMB(HumanCointeractionEigenvectorPercentile ~ Sex + RankPercentile + Species + PropScansHumanProximityPercentile + (1|Group), data=data, family = beta_family())
mod.evc.groomingandaffiliation <- glmmTMB(HumanCointeractionEigenvectorPercentile ~ Sex + RankPercentile + Species + GroomEigenvectorPercentile + PropScansHumanProximityPercentile + OtherAffiliationEigenvectorPercentile + (1|Group), data=data, family = beta_family())
mod.evc.proximityandaffiliation <- glmmTMB(HumanCointeractionEigenvectorPercentile ~ Sex + RankPercentile + Species + ProximityEigenvectorPercentile + PropScansHumanProximityPercentile + OtherAffiliationEigenvectorPercentile + (1|Group), data=data, family = beta_family())
mod.evc.groomingandaffiliation.groomingbyspecies <- glmmTMB(HumanCointeractionEigenvectorPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + GroomEigenvectorPercentile*Species + OtherAffiliationEigenvectorPercentile + (1|Group), data=data, family = beta_family())
mod.evc.groomingandaffiliation.affiliationbyspecies <- glmmTMB(HumanCointeractionEigenvectorPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + GroomEigenvectorPercentile + OtherAffiliationEigenvectorPercentile*Species + (1|Group), data=data, family = beta_family())
mod.evc.proximityandaffiliation.proximitybyspecies <- glmmTMB(HumanCointeractionEigenvectorPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + ProximityEigenvectorPercentile*Species + OtherAffiliationEigenvectorPercentile + (1|Group), data=data, family = beta_family())
mod.evc.proximityandaffiliation.affiliationbyspecies <- glmmTMB(HumanCointeractionEigenvectorPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + ProximityEigenvectorPercentile + OtherAffiliationEigenvectorPercentile*Species + (1|Group), data=data, family = beta_family())

##Get the AICc scores to determine the best-fit model

AICc(mod.evc.nonetworkmeasures)
AICc(mod.evc.groomingandaffiliation)
AICc(mod.evc.proximityandaffiliation)
AICc(mod.evc.groomingandaffiliation.groomingbyspecies)
AICc(mod.evc.groomingandaffiliation.affiliationbyspecies)
AICc(mod.evc.proximityandaffiliation.proximitybyspecies)
##This model should have the lowest AICc score
AICc(mod.evc.proximityandaffiliation.affiliationbyspecies)

##Get model summary statistics

summary(mod.evc.nonetworkmeasures)
summary(mod.evc.groomingandaffiliation)
summary(mod.evc.proximityandaffiliation)
summary(mod.evc.groomingandaffiliation.groomingbyspecies)
summary(mod.evc.groomingandaffiliation.affiliationbyspecies)
summary(mod.evc.proximityandaffiliation.proximitybyspecies)
summary(mod.evc.proximityandaffiliation.affiliationbyspecies)


##check validity of best-fit####

resideEvc <- resid (mod.evc.proximityandaffiliation.proximitybyspecies)
hist (resideEvc)
diagnostics.plot(mod.evc.proximityandaffiliation.proximitybyspecies)

##Get model-predicted values using the "predict" function 

data$ModPredEvc <- predict(mod.evc.proximityandaffiliation.proximitybyspecies, type="response")

##use relevel function to get interaction coefficients for each combination of species ####

data$Species <- relevel(data$Species, ref = "Rhesus")
mod.evc.proximityandaffiliation.proximitybyspeciesRhesus <- glmmTMB(HumanCointeractionEigenvectorPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + ProximityEigenvectorPercentile*Species + OtherAffiliationEigenvectorPercentile + (1|Group), data=data, family = beta_family())
summary(mod.evc.proximityandaffiliation.proximitybyspeciesRhesus)

data$Species <- relevel(data$Species, ref = "Long-tailed")
mod.evc.proximityandaffiliation.proximitybyspeciesLTM <- glmmTMB(HumanCointeractionEigenvectorPercentile ~ Sex + RankPercentile + PropScansHumanProximityPercentile + ProximityEigenvectorPercentile*Species + OtherAffiliationEigenvectorPercentile + (1|Group), data=data, family = beta_family())
summary(mod.evc.proximityandaffiliation.proximitybyspeciesLTM)


#Get model-predicted values 

data$ModPredEvc <- predict(mod.evc.proximityandaffiliation.proximitybyspecies, type="response")


###Plot the significant effects

#Plot the effect of affiliation eigenvector on human co-interaction eigenvector

plotaffiliationevc <- ggplot(data, aes(x = OtherAffiliationEigenvectorPercentile, y = HumanCointeractionEigenvectorPercentile)) +
  geom_point(size = 1, alpha = 0.9) +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network eigenvector", x = "Short-duration affiliation network eigenvector", tag = "a") +
  theme_cowplot() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotaffiliationevc


#Plot the effect of affiliation eigenvector on model-predicted Human co-interaction eigenvector

plotaffiliationevcmodpred <- ggplot(data, aes(x = OtherAffiliationEigenvectorPercentile, y = ModPredEvc)) +
  geom_point(size = 1, alpha = 0.9) +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network eigenvector", x = "Short-duration affiliation network eigenvector", tag = "a") +
  theme_cowplot() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotaffiliationevcmodpred


#Plot the effect of proximity eigenvector by species on human co-interaction eigenvector

plotproximityevcbyspecies <- ggplot(data, aes(x = ProximityEigenvectorPercentile, y = HumanCointeractionEigenvectorPercentile, color = Species)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network eigenvector", x = "Proximity network eigenvector", tag = "b") +
  theme_cowplot() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotproximityevcbyspecies


#Plot the effect of proximity eigenvector by species on model-predicted human co-interaction eigenvector

plotproximityevcbyspeciesmodpred <- ggplot(data, aes(x = ProximityEigenvectorPercentile, y = ModPredEvc, color = Species)) +
  scale_color_manual(values=c("grey1","grey60", "blue")) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network eigenvector", x = "Proximity network eigenvector", tag = "b") +
  theme_cowplot() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotproximityevcbyspeciesmodpred


#Plot the effect of sex on human co-interaction network eigenvector

plotSexEvc <- ggplot(data, aes(x = Sex, y = HumanCointeractionEigenvectorPercentile)) + 
  geom_boxplot() +
  geom_point(size = 1, alpha = 0.9) +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network eigenvector", x = "Macaque sex") +
  theme_cowplot() +
  labs(title = "Tau = 0.10") +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotSexEvc


#Plot the effect of sex on model-predicted human co-interaction network eigenvector

plotSexEvcModPred <- ggplot(data, aes(x = Sex, y = ModPredEvc)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  geom_point(size = 1, alpha = 0.9) +
  geom_smooth(method = "lm") +
  labs(y = "Human co-interaction network eigenvector", x = "Macaque sex", tag = "c") +
  theme_cowplot() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  theme(#axis.title.x = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    axis.text.x = element_text(size=12, color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    strip.text.x = element_text(size = 12),
    plot.tag = element_text(size=14))

plotSexEvcModPred




