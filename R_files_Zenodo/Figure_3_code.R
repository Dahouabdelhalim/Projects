### This code produces Figure 3 of Stewart Merrill et al. submitted

# Required packages
library(cowplot) # for consolidating figure panels
library(dplyr) # for organizing data
library(ggplot2) # for creating figures
library(glmmTMB) # for statistical analyses

##############
# PREPARATION
##############

# Import data
data <- read.csv("master_SIN_microscopic evaluations.csv")
# In 'data', individual Daphnia are the unit of replication

###########
# ANALYSES
###########

# For attacks and successful attacks, run a generalized linear mixed model
# with a negative binomial error distribution and tube as a random effect
# To facilitate model convergence, we sometimes needed to update the optimizer
# While our standard form of the negative binomial was "nbinom2", we sometimes
# needed to use "nbinom1" for model convergence

# Sub-setting to run the models separately at each dose
low  <- subset(data, data$InitialDoseTreatment == "low")
high <- subset(data, data$InitialDoseTreatment == "high")

# Sample sizes for supplement
table(low$TubeReplicate)
table(high$TubeReplicate)

# Spore attack by host density at the low initial dose
# Requires an update of the optimizer for convergence
mdl1 <- glmmTMB(AttackingSpores ~ HostDensityTreatment + 
                (1| TubeReplicate),
                family = nbinom2, 
                data = low)
mdl1 <- update(mdl1, control=glmmTMBControl(optimizer=optim,
                                                  optArgs=list(method="BFGS")))
summary(mdl1)
# We detect a significant negative fixed effect of host density

# Successful attack by host density at low initial dose
# Requires "nbinom1" and an update of the optimizer for convergence
mdl2 <- glmmTMB(SuccessfulAttacks ~ HostDensityTreatment + 
                (1 | TubeReplicate),
                family = nbinom1, 
                data = low)
mdl2 <- update(mdl2, control=glmmTMBControl(optimizer=optim,
                                            optArgs=list(method="BFGS")))
summary(mdl2)
# There is no significant fixed effect of host density

# Rerun the model as intercept only for cumulative successful attacks
mdl2.b <- glmmTMB(SuccessfulAttacks ~ 1 + 
                (1 | TubeReplicate),
                family = nbinom1, 
                data = low)
mdl2.b <- update(mdl2.b, control=glmmTMBControl(optimizer=optim,
                                            optArgs=list(method="BFGS")))

# Spore attack by host density at the high initial dose
mdl3 <- glmmTMB(AttackingSpores ~ HostDensityTreatment + 
                (1 | TubeReplicate),
                family = nbinom2, 
                data = high)
summary(mdl3)
# There is no significant fixed effect of host density

# Rerun the model as intercept only for cumulative attacks
# Requires an update of the optimizer for convergence
mdl3.b <- glmmTMB(AttackingSpores ~ 1 + 
                (1 | TubeReplicate),
                family = nbinom2, 
                data = high)
mdl3.b <- update(mdl3.b, control=glmmTMBControl(optimizer=optim,
                                                optArgs=list(method="BFGS")))

# Successful attacks by host density at the high initial dose
mdl4 <- glmmTMB(SuccessfulAttacks ~ HostDensityTreatment + 
                (1 | TubeReplicate),
                family = nbinom2, 
                data = high)
summary(mdl4)
# There is no significant fixed effect of host density

# Rerun the model as intercept only for cumulative successful attacks
mdl4.b <- glmmTMB(SuccessfulAttacks ~ 1 + 
                (1 | TubeReplicate),
                family = nbinom2, 
                data = high)

# Note: these models were originally run with poisson, negative binomial, and negative binomial with zero inflation
# The qualitative results did not change across model types, and a reviewer suggested proceeding with
# negative binomial to reduce issues of over-dispersion

# For the interested reader, we also include a full model that
# includes both doses combined (results provided in Supplement section S2)

# We need a unique tube replicate across both the doses
# (because the name of a tube replicate could be repeated in low and high doses,
# even though they are different tubes)
data$UniqueTubeRep <- paste0(data$InitialDoseTreatment, "_", data$TubeReplicate)

# The full model for spore attack
full.mdl1 <- glmmTMB(AttackingSpores ~ HostDensityTreatment*InitialDoseTreatment + (1| UniqueTubeRep),
                family = nbinom2, data = data)
car::Anova(full.mdl1)

# The full model for successful attacks
full.mdl2 <- glmmTMB(SuccessfulAttacks ~ HostDensityTreatment*InitialDoseTreatment + (1| UniqueTubeRep),
                family = nbinom2, data = data)
car::Anova(full.mdl2)

##########################################################
# GENERATING MODEL PREDICTIONS TO OVERLAY CURVES ON PLOTS
##########################################################

# Attack by density at low initial dose
newdata1 <- data.frame(HostDensityTreatment = seq(0, 16, 0.01), TubeReplicate = NA)
newdata1$pred <- predict(mdl1, newdata = newdata1, type = "response", re.form = ~0)
newdata1$cumulative <- newdata1$pred*newdata1$HostDensityTreatment
# where 'pred' is the per-capita prediction, and 'cumulative' is the cumulative prediction

# Infection by density at low initial dose
newdata2 <- data.frame(HostDensityTreatment = seq(0, 16, 0.01), TubeReplicate = NA)
newdata2$pred <- predict(mdl2.b, newdata = newdata2, type = "response", re.form = ~0)
newdata2$cumulative <- newdata2$pred*newdata2$HostDensityTreatment

# Attack by density at high initial dose
newdata3 <- data.frame(HostDensityTreatment = seq(0, 16, 0.01), TubeReplicate = NA)
newdata3$pred <- predict(mdl3.b, newdata = newdata3, type = "response", re.form = ~0)
newdata3$cumulative <- newdata3$pred*newdata3$HostDensityTreatment

# Infection by density at high initial dose
newdata4 <- data.frame(HostDensityTreatment = seq(0, 16, 0.01), TubeReplicate = NA)
newdata4$pred <- predict(mdl4.b, newdata = newdata4, type = "response", re.form = ~0)
newdata4$cumulative <- newdata4$pred*newdata4$HostDensityTreatment

#########
# FIGURE
#########

# Figure presets
presets <- theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.border=element_blank(),
        axis.line = element_line(colour = 'black', size = 0.5),
        axis.ticks=element_line(colour = 'black', size = 0.5),
        axis.title.x = element_text(family="sans", size=13, color="black", margin=margin(12,0,0,0)),
        axis.title.y = element_text(family="sans", size=13, color="black", margin=margin(0,12,0,0)),
        axis.text.x = element_text(family="sans", size=10, color="black"),
        axis.text.y = element_text(family="sans", size=10, color="black"),
        strip.background = element_blank(),
        strip.text.x = element_blank())

########################
# PER-CAPITA MAIN PLOTS
########################

# Attacking spores in the low initial dose
fig3a <- ggplot() +
    geom_line(data = newdata1, aes(x = HostDensityTreatment/50, y = pred), col = "black", lwd = 1) +
    geom_jitter(data = low, aes(x = HostDensityTreatment/50, y = AttackingSpores),
                size = 1.5, pch = 21, col = "black", fill = "cyan3", height = 0.0, width = 0.01) +
    scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
    xlab("Density treatment (hosts per ml)") +
    ylab("Attacking spores per host") +
    ylim(0,12) +
    presets

# Attacking spores in the high initial dose
fig3b <- ggplot() +
    geom_line(data = newdata3, aes(x = HostDensityTreatment/50, y = pred), col = "black", lwd = 1, lty = "dashed") +
    geom_jitter(data = high, aes(x = HostDensityTreatment/50, y = AttackingSpores),
                size = 1.5, pch = 21, col = "black", fill = "brown3", height = 0.0, width = 0.01) +
    scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
    xlab("Density treatment (hosts per ml)") +
    ylab("Attacking spores per host") +
    ylim(0,12) +
    presets

# Successful attacks in the low initial dose
fig3c <- ggplot() +
    geom_line(data = newdata2, aes(x = HostDensityTreatment/50, y = pred), col = "black", lwd = 1, lty = "dashed") +
    geom_jitter(data = low, aes(x = HostDensityTreatment/50, y = SuccessfulAttacks),
                size = 1.5, pch = 21, col = "black", fill = "cyan3", height = 0.0, width = 0.01) +
    scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
    xlab("Density treatment (hosts per ml)") +
    ylab("Successful attacks per host") +
    ylim(0,12) +
    presets

# Successful attacks in the high initial dose
fig3d <- ggplot() +
    geom_line(data = newdata4, aes(x = HostDensityTreatment/50, y = pred), col = "black", lwd = 1, lty = "dashed") +
    geom_jitter(data = low, aes(x = HostDensityTreatment/50, y = SuccessfulAttacks),
                size = 1.5, pch = 21, col = "black", fill = "brown3", height = 0.0, width = 0.01) +
    scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
    xlab("Density treatment (hosts per ml)") +
    ylab("Successful attacks per host") +
    ylim(0,12) +
    presets

plot_grid(fig3a, fig3b, fig3c, fig3d, ncol = 2)

###############################
# CUMULATIVE SUBPLOTS (INLAID)
###############################

# Attacking spores in the low intial dose
sp3a <- ggplot(data = newdata1, aes(x = HostDensityTreatment/50, y = cumulative)) +
        geom_smooth(col = "cyan3", lwd = 3) +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        ylim(0,20) +
        xlab("Density treatment (hosts per ml") +
        ylab("Cumulative attack") +
        theme(legend.position = "none") +
        presets

# Attacking spores in the high initial dose
sp3b <- ggplot(data = newdata3, aes(x = HostDensityTreatment/50, y = cumulative)) +
    geom_smooth(col = "brown3", lwd = 3) +
    scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
    ylim(0,60) +
    xlab("Density treatment (hosts per ml)") +
    ylab("Cumulative attack") +
    theme(legend.position = "none") +
    presets

# Successful attacks in the low initial dose
sp3c <- ggplot(data = newdata2, aes(x = HostDensityTreatment/50, y = cumulative)) +
        geom_smooth(col = "cyan3", lwd = 3) +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        ylim(0,40) +
        xlab("Density treatment (hosts per ml") +
        ylab("Cumulative success") +
        theme(legend.position = "none") +
        presets

# Successful attacks in the high initial dose
sp3d <- ggplot(data = newdata4, aes(x = HostDensityTreatment/50, y = cumulative)) +
        geom_smooth(col = "brown3", lwd = 3) +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        ylim(0,60) +
        xlab("Density treatment (hosts per ml") +
        ylab("Cumulative success") +
        theme(legend.position = "none") +
        presets

##############
# END OF CODE
##############