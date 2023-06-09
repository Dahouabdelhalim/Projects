#=========================================================================================
# R script for: MK Rhodes, JB Fant & KA Skogen (2017) Pollinator identity and spatial    #
# isolation influence multiple paternity in an annual plant.                             #
#                                                                                        #
# Script by MK Rhodes. Last updated 26 January 2017.                                     #
#=========================================================================================

# << DIRECTORIES >> ----------------------------------------------------------------------
# >> Set your directory here
setwd("~/")

# << POLLINATOR VISITATION ANALYSES >> ---------------------------------------------------
# >> DATA
pmObs <- read.csv("pmPollinators.csv")
amObs <- read.csv("amPollinators.csv")
obs <- read.csv("overallPollinators.csv")
bouts <- read.csv("pollinatorBouts.csv")

# >> FORAGING BOUTS PER FLOWER PER HOUR
# Comparisons for evening observations
mothPM <- pmObs[pmObs$pollinator == "moth", "abundance"]
beePM <- pmObs[pmObs$pollinator == "bee", "abundance"]
wilcox.test(mothPM,beePM)

# Comparisons for morning observations
mothAM <- amObs[amObs$pollinator == "moth", "abundance"]
beeAM <- amObs[amObs$pollinator == "bee", "abundance"]
wilcox.test(mothAM,beeAM)

# Overall comparison for evening & morning observations combined
moth <- obs[obs$pollinator == "moth", "abundance"]
bee <- obs[obs$pollinator =="bee", "abundance"]
wilcox.test(moth, bee)

# >> NUMBER OF DIFFERENT FLOWERS VISITED PER BOUT
hawkmoth <- bouts[bouts$visitor=="hawkmoth", "flowersVisited"]
bee <- bouts[bouts$visitor=="solitary bee", "flowersVisited"]
wilcox.test(hawkmoth, bee)

# >> LIKELIHOOD OF MULTIPLE-FLOWER VISITS
# Transform number of flowers visited into binary response
bouts$multipleFlowers <- ifelse(bouts$flowersVisited > 1, 1, 0)

# Model specification & selection
multi <- glm(multipleFlowers ~ visitor, quasibinomial, data = bouts)
multiNull <- glm(multipleFlowers ~ 1, quasibinomial, data = bouts)
anova(multiNull, multi, test = "Chi")
summary(multi)

# << POLLEN DISPERSAL ANALYSES >> --------------------------------------------------------
# >> DATA
dispersal <- read.csv("pollenDispersal.csv")

# >> INFLUENCE OF POLLINATOR IDENTITY ON POLLEN DISPERSAL DISTANCE
# Model specification & selection
dMod <- lm(log10(distance + 1) ~ treatment, data = dispersal)
dNull <- lm(log10(distance + 1) ~ 1, data = dispersal)
anova(dNull, dMod)
summary(dMod)

# Summary metrics of pattens in pollen dispersal by treatment
control <- dispersal[dispersal$treatment=="C", "distance"]
summary(control)

bee <- dispersal[dispersal$treatment=="NE", "distance"]
summary(bee)

hawkmoth <- dispersal[dispersal$treatment=="DE", "distance"]
summary(hawkmoth)

# << MULTIPLE PATERNITY ~ POLLINATORS AND SPATIAL ISOLATION ANALYSIS >> ------------------
# >> DATA
mates <- read.csv("multiplePaternity.csv")

# >> CREATE RESPONSE & PREDICTOR VARIABLES
# Translate correlated paternity to number of pollen donors
mates$effectiveDads <- round(1/mates$correlatedPaternity)

# Transform multiple mating into binary response
mates$multipleMating <- ifelse(mates$effectiveDads > 1, 1, 0)

# Transform spatial isolation into binary predictor
mates$isolation <- ifelse(mates$isolation20 > 1000, 1, 0)

# >> TEST FOR INFLUENCE OF MATERNAL LINE IDENTITY
# Model specification & selection
t1 <- glm(multipleMating ~ treatment*isolation*plantID, binomial, data = mates)
t2 <- glm(multipleMating ~ treatment*isolation + treatment*plantID + 
            isolation*plantID + treatment + isolation + plantID, binomial,
          data = mates)
anova(t2, t1, test = "Chi")
# three-way interaction is not significant
summary(t2)
# remove treatment*plantID

t3 <- glm(multipleMating ~ treatment*isolation + isolation*plantID + treatment
          + isolation + plantID, binomial, data = mates)
anova(t3, t2, test = "Chi")
# treatment*plantID is not significant
summary(t3)
# remove isolation*plantID

t4 <- glm(multipleMating ~ treatment*isolation + treatment + isolation + 
            plantID, binomial, data = mates)
anova(t4, t3, test = "Chi")
# isolation*plantID is not significant
summary(t4)
# remove plantID

t5 <- glm(multipleMating ~ treatment*isolation + treatment + isolation, 
          binomial, data = mates)
anova(t5, t4, test = "Chi")
# plantID is not significant
# >> Maternal line ID does not appear to be a significant confounding factor; final models
# will be specified using only pollinator identity and spatial isolation as predictors.

# >> FINAL MODELS
# Model specification & selection
m1 <- glm(multipleMating ~ treatment*isolation, binomial, data = mates)
m2 <- glm(multipleMating ~ treatment+isolation, binomial, data = mates)
anova(m2, m1, test = "Chi")
summary(m2)

m3 <- glm(multipleMating ~ isolation, binomial, data = mates)
anova(m3, m2, test = "Chi")

m4 <- glm(multipleMating ~ treatment, binomial, data = mates)
anova(m4, m2, test = "Chi")

summary(m2)
# This is the minimum adequate model

# >> EXTRACT MODEL PREDICTIONS AND APPEND TO DATAFRAME
polyandry.response.fit <- predict.glm(m2, type = "response", se.fit = TRUE)
mates$polyandry.response.fit <- polyandry.response.fit$fit
mates$polyandry.response.fit.se <- polyandry.response.fit$se.fit

# >> INTERPRET COEFFICIENTS IN MINIMUM ADEQUATE MODEL
1/exp(m2$coefficients["treatmentde"][[1]])
# The odds of open-pollinated seed families being multiply sired were 1.81x 
# the odds of hawkmoth-pollinated seed families being multiply sired.
# This difference is NOT statistically significant.

1/exp(m2$coefficients["treatmentne"][[1]])
# The odds of open-pollinated seed families being multiply sired were 10.66x
# the odds of bee-pollinated seed families being multiply sired.
# This difference IS statistically significant.

1/exp(m2$coefficients["isolation"][[1]])
# The odds of seed families from non-isolated individuals being multiply sired were 8.72x 
# the odds of seed families from isolated individuals being multiply sired.
# This difference IS statistically significant.

# >> SUMMARIZE INFLUENCE OF MULTIPLE PATERNITY IN A SPATIAL CONTEXT
(unique(mates$polyandry.response.fit[mates$treatment == "c" & mates$isolation == "0"]) -
    unique(mates$polyandry.response.fit[mates$treatment == "c" & mates$isolation == "1"]))/
  unique(mates$polyandry.response.fit[mates$treatment == "c" & mates$isolation == "0"])
# For open-pollinated flowers, the predicted incidence of multiple paternity decreases 
# ~ 30% for isolated individuals

(unique(mates$polyandry.response.fit[mates$treatment == "de" & mates$isolation == "0"]) -
    unique(mates$polyandry.response.fit[mates$treatment == "de" & mates$isolation == "1"]))/
  unique(mates$polyandry.response.fit[mates$treatment == "de" & mates$isolation == "0"])
# For hawkmoth-pollinated flowers, the predicted incidence of multiple paternity decreases 
# ~ 43% for isolated individuals

(unique(mates$polyandry.response.fit[mates$treatment == "ne" & mates$isolation == "0"]) -
    unique(mates$polyandry.response.fit[mates$treatment == "ne" & mates$isolation == "1"]))/
  unique(mates$polyandry.response.fit[mates$treatment == "ne" & mates$isolation == "0"])
# For bee-pollinated flowers, the predicted incidence of multiple paternity decreases 
# ~ 75% for isolated individuals

# << FIGURES >> --------------------------------------------------------------------------
# >> PACKAGES
library("ggplot2")

# >> FIGURE 2: POLLINATOR ABUNDANCE
# data
pa <- read.csv("obs.csv")

# plot
limits <- aes(ymax = rate + se, ymin = rate)
dodge <- position_dodge(width = 1)
fig2 <- ggplot(data=pa, aes(x=time, y=rate, fill=insect))
fig2 +
  geom_bar(position = dodge, stat = "identity", colour = "black") +
  geom_errorbar(limits, position = dodge, width = 0.1, colour = "black") +
  scale_fill_manual(values=c("grey55", "white")) +
  xlab("Observation Time") +
  ylab("Foraging Bouts per Flower per Hour") +
  scale_x_discrete(labels = c("Evening (n = 17 hours)", "Morning (n = 12 hours)")) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3), breaks = c(0,1,2,3)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.ticks.x = element_blank(),
        legend.position = "none")

# >> FIGURE 3: POLLEN DISPERSAL DISTANCE
fig3 <- ggplot(data=dispersal, aes(x = treatment, y = log10(distance + 1)))
fig3 +
  geom_boxplot(colour = "black") +
  xlab("Exclusion Treatment") +
  ylab("Pollen Dispersal Distance (m)") +
  scale_x_discrete(labels = c("Open-pollinated (C)", "Hawkmoth (DE)", "Solitary Bee (NE)")) +
  scale_y_continuous(breaks = c(0,1,2,3), labels = c(0,10,100,1000)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.ticks.x = element_blank(),
        legend.position = "none")

# >> FIGURE 4: MULTIPLE PATERNITY ~ POLLINATORS AND SPATIAL ISOLATION
limits <- aes(ymax = polyandry.response.fit + polyandry.response.fit.se,
              ymin = polyandry.response.fit - polyandry.response.fit.se)
fig4 <- ggplot(mates, aes(x = isolation, y = polyandry.response.fit, group = treatment))
fig4 +  
  geom_point(size = 3, aes(shape = factor(treatment)), position = position_dodge(width = 0.5)) +
  #geom_errorbar(limits, size = 0.5, width = 0.1, position = position_dodge(width = 0.5), colour = "black") +
  geom_pointrange(limits, size = 0.35, position = position_dodge(width = 0.5)) +
  scale_x_continuous(breaks = c(0,1), labels = c("Not Isolated", "Isolated")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1.0)) +
  xlab("Spatial Isolation") +
  ylab("Incidence of Multiple Paternity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.ticks.x = element_blank(),
        legend.position = "none")