### This code produces Figure 2 of Stewart Merrill et al. submitted

# Required packages
library(cowplot) # for consolidating figure panels
library(dplyr) # for organizing data
library(ggplot2) # for creating figures
library(glmmTMB) # for statistical analyses
library(plotrix) # for calculating standard error
library(tidyr) # for data cleaning/manipulation

##############
# PREPARATION
##############

# Import prevalence data
data <- read.csv('master_SIN_lab infection prevalence.csv')
# In 'data', individual Daphnia are the unit of replication

# Remove individuals with unknown infection statuses (those that died prematurely)
data <- subset(data, data$InfectionStatus != "")

# Import spore count data
spor <- read.csv("master_SIN_lab spore counts.csv")
# In 'spor', tubes are the unit of replication, with total spores from the tube indicated
# These are only tubes that had infected Daphnia
# These data must be merged with uninfected tubes to build in the zero spore values (see Line 60)

##############################
# ORGANIZING DATA FOR FIGURES
##############################

# For prevalence figure panels 2a, 2c, 2e, and 2g:

# Quantify prevalence for each tube replicate
data$status <- as.numeric(recode(data$InfectionStatus, i = 1, u = 0))
pdata <- data %>%
	        group_by(InitialDoseTreatment, Generation, HostDensityTreatment, TubeReplicate) %>%
		summarise_at(vars("status"), mean)
pdata <- as.data.frame(pdata)
# In 'pdata', the tube is the unit of replication and 'status' is prevalence
# Note: prevalence values are calculated from surviving Daphnia, so prevalence*density will not
# return the number of infected individuals

# Now combine tube replicates to quantify average prevalence for each spore dose, generation,
# and host density treatment
prev <- pdata %>%
		group_by(InitialDoseTreatment, Generation, HostDensityTreatment) %>%
		summarise_at(vars("status"), funs(mean, std.error))
prev <- as.data.frame(prev)
# The 'prev' dataframe is now ready to create the figure

# For spore count figure panels 2b, 2d, 2f, and 2h:

# Generate spores per ml from total spore counts by dividing by 50 ml water volume
spor$Spores.per.ml <- spor$TotalSporesProduced / 50

# Reduce the dataframe to columns of interest
spor <- select(spor, c('InitialDoseTreatment', 'Generation', 'HostDensityTreatment', 'TubeReplicate', 'Spores.per.ml'))

# Merge spore totals from each tube ('spor') to prevalence values for each tube ('pdata')
spor.data <- merge(pdata, spor, by=c('InitialDoseTreatment', 'Generation', 'HostDensityTreatment', 'TubeReplicate'), all.x = T, all.y = T)

# Tubes with no spore value receive a 0
spor.data$Spores.per.ml[is.na(spor.data$Spores.per.ml)] <- 0

# Confirm that all zero prevalence tubes have zero spores
plot(spor.data$status, spor.data$Spores.per.ml) 

# Combine tube replicates to quantify average spore counts for each spore dose, generation,
# and host density treatment
sdens <- spor.data %>%
        group_by(InitialDoseTreatment, Generation, HostDensityTreatment) %>%
        summarise_at(vars("Spores.per.ml"), funs(mean, std.error))
sdens <- as.data.frame(sdens)
# The 'sdens' dataframe is now ready to create the figure

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

# Prevalence plot for low initial dose, first generation
fig2a <- ggplot(data = subset(prev, prev$InitialDoseTreatment == "low" & prev$Generation == 1), aes(x = HostDensityTreatment/50, y = mean)) +
        geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0) +
        geom_line(alpha = 0.6) +
        geom_point(size = 3, pch = 21, col = "black", fill = "cyan3") +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        scale_y_continuous(limits = c(0, 0.65)) +
        xlab("Host density (hosts per ml") +
        ylab("Infection prevalence") +
        presets

# Spore plot for low initial dose, first generation
fig2b <- ggplot(data = subset(sdens, sdens$InitialDoseTreatment == "low" & sdens$Generation == 1), aes(x = HostDensityTreatment/50, y = mean)) +
        geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0) +
        geom_line(alpha = 0.6) +
        geom_point(size = 3, pch = 21, col = "black", fill = "cyan3") +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        coord_cartesian(ylim = c(0, 6500)) +
        xlab("Density treatment (hosts per ml)") +
        ylab("Spore production") +
        presets

# Prevalence plot for low initial dose, second generation
fig2c <- ggplot(data = subset(prev, prev$InitialDoseTreatment == "low" & prev$Generation == 2), aes(x = HostDensityTreatment/50, y = mean)) +
        geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0) +
        geom_line(alpha = 0.6) +
        geom_point(size = 3, pch = 21, col = "black", fill = "cyan3") +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        scale_y_continuous(limits = c(0, 0.65)) +
        xlab("Host density (hosts per ml") +
        ylab("Infection prevalence") +
        presets

# Spore plot for low initial dose, second generation
fig2d <- ggplot(data = subset(sdens, sdens$InitialDoseTreatment == "low" & sdens$Generation == 2), aes(x = HostDensityTreatment/50, y = mean)) +
        geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0) +
        geom_line(alpha = 0.6) +
        geom_point(size = 3, pch = 21, col = "black", fill = "cyan3") +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        coord_cartesian(ylim = c(0, 6500)) +
        xlab("Density treatment (hosts per ml)") +
        ylab("Spore production") +
        presets

# Prevalence plot for high initial dose, first generation        
fig2e <- ggplot(data = subset(prev, prev$InitialDoseTreatment == "high" & prev$Generation == 1), aes(x = HostDensityTreatment/50, y = mean)) +
        geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0) +
        geom_line(alpha = 0.6) +
        geom_point(size = 3, pch = 21, col = "black", fill = "brown3") +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        scale_y_continuous(limits = c(0, 0.65)) +
        xlab("Host density (hosts per ml") +
        ylab("Infection prevalence") +
        presets

# Spore plot for high initial dose, first generation
fig2f <- ggplot(data = subset(sdens, sdens$InitialDoseTreatment == "high" & sdens$Generation == 1), aes(x = HostDensityTreatment/50, y = mean)) +
        geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0) +
        geom_line(alpha = 0.6) +
        geom_point(size = 3, pch = 21, col = "black", fill = "brown3") +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        coord_cartesian(ylim = c(0, 6500)) +
        xlab("Density treatment (hosts per ml)") +
        ylab("Spore production") +
        presets

# Prevalence plot for high initial dose, second generation 
fig2g <- ggplot(data = subset(prev, prev$InitialDoseTreatment == "high" & prev$Generation == 2), aes(x = HostDensityTreatment/50, y = mean)) +
        geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0) +
        geom_line(alpha = 0.6) +
        geom_point(size = 3, pch = 21, col = "black", fill = "brown3") +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        scale_y_continuous(limits = c(0, 0.65)) +
        xlab("Host density (hosts per ml") +
        ylab("Infection prevalence") +
        presets

# Spore plot for high initial dose, second generation
fig2h <- ggplot(data = subset(sdens, sdens$InitialDoseTreatment == "high" & sdens$Generation == 2), aes(x = HostDensityTreatment/50, y = mean)) +
        geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0) +
        geom_line(alpha = 0.6) +
        geom_point(size = 3, pch = 21, col = "black", fill = "brown3") +
        scale_x_continuous(breaks = c(0.02,0.08,0.16,0.24,0.32)) +
        coord_cartesian(ylim = c(0, 6500)) +
        xlab("Density treatment (hosts per ml)") +
        ylab("Spore production") +
        presets

plot_grid(fig2a, fig2b, fig2c, fig2d, fig2e, fig2f, fig2g, fig2h, ncol = 4)

#####################
# HOST RISK ANALYSES
#####################

# For host risk, we run a generalized linear mixed model
# with a binomial distribution and tube as a random effect
# We compare the full models to their nested variants with AIC

# For individual host risk, use 'data' which has the individual as the unit of replication
high <- subset(data, data$InitialDoseTreatment=="high") # Data for high initial dose
low <- subset(data, data$InitialDoseTreatment=="low") # Data for low initial dose

# High initial dose treatment analyses
# These results are reported in Table S1
### Full model: All fixed effects and interactions
mdl1 <- glmmTMB(status ~ HostDensityTreatment*Generation + 
                        (1 | TubeReplicate), 
                        family = binomial, 
                        data = high)
summary(mdl1)

### Nested variant: Host density and generation as fixed effects
mdl1.b <- glmmTMB(status ~ HostDensityTreatment + Generation + 
                          (1 | TubeReplicate), 
                  family = binomial, 
                  data = high)
summary(mdl1.b)

### Nested variant: Host density as sole fixed efffect
mdl1.c <- glmmTMB(status ~ HostDensityTreatment +
                          (1 | TubeReplicate), 
                  family = binomial, 
                  data = high)
summary(mdl1.c)

### Null model: intercept only
mdl1.d <- glmmTMB(status ~ 1 +
                          (1 | TubeReplicate), 
                  family = binomial, 
                  data = high)
summary(mdl1.d)

AIC(mdl1)
AIC(mdl1.b)
AIC(mdl1.c)
AIC(mdl1.d)

# Low initial dose treatment analyses
# These results are reported in Table S1
### Full model: All fixed effects and interactions
mdl2 <- glmmTMB(status ~ HostDensityTreatment*Generation +
                        (1 | TubeReplicate),
                        family = binomial, 
                        data = low)
summary(mdl2)

### Nested variant: Host density and generation as fixed effects
mdl2.b <- glmmTMB(status ~ HostDensityTreatment + Generation +
                        (1 | TubeReplicate),
                family = binomial, 
                data = low)
summary(mdl2.b)

### Nested variant: Host density as sole fixed effect
mdl2.c <- glmmTMB(status ~ HostDensityTreatment + 
                          (1 | TubeReplicate),
                  family = binomial, 
                  data = low)
summary(mdl2.c)

### Null model: intercept only
mdl2.d <- glmmTMB(status ~ 1 + 
                          (1 | TubeReplicate),
                  family = binomial, 
                  data = low)
summary(mdl2.d)

AIC(mdl2)
AIC(mdl2.b)
AIC(mdl2.c)
AIC(mdl2.d)

# For the interested reader, we also include a full model that
# includes both doses combined (results provided in Supplement section S2)

# For the full model, we need to create a new random effect that includes dose in its name
# because there can be a tube "1A" that appears in the low dose and one in the high dose.
data$UniqueTubeRep <- paste0(data$initial.dose.treatment, "_", data$tube.replicate)

# Full model combining both doses
mdlfull <- glmmTMB(status ~ HostDensityTreatment*Generation*InitialDoseTreatment +
                            (1 | UniqueTubeRep),
                    family = binomial,
                    data = data)
car::Anova(mdlfull)

############################
# SPORE PRODUCTION ANALYSES
############################

# For spore production, we run a generalized linear model for count data
# with a zero-inflated negative binomial distribution
# We compare the full models to their nested variants with AIC

# For spore production, use 'spor.data', where tube is the unit of replication
hspor <- subset(spor.data, spor.data$InitialDoseTreatment == "high") # Data for high initial dose
lspor <- subset(spor.data, spor.data$InitialDoseTreatment == "low") # Data for low initial dose
hspor$Spores.per.ml <- round(hspor$Spores.per.ml) # Rounding spores per ml for count data
lspor$Spores.per.ml <- round(lspor$Spores.per.ml) # Rounding spores per ml for count data

# High initial dose treatment analyses
# These results are reported in Table S2
### Full model: All fixed effects and interactions
# Density needed to be scaled to facilitate model convergence in nested variants
mdl3 <- glmmTMB(Spores.per.ml ~ scale(HostDensityTreatment)*Generation, 
                ziformula = ~1, 
                family = nbinom2, 
                data = hspor)
summary(mdl3)

### Nested variant: Host density and generation as fixed effects
mdl3.b <- glmmTMB(Spores.per.ml ~ scale(HostDensityTreatment) + Generation, 
                  ziformula = ~1, 
                  family = nbinom2, 
                  data = hspor)
summary(mdl3.b)

### Nested variant: Host density as sole fixed effect
mdl3.c <- glmmTMB(Spores.per.ml ~ scale(HostDensityTreatment), 
                  ziformula = ~1, 
                  family = nbinom2, 
                  data = hspor)
summary(mdl3.c)

### Null model: intercept only
mdl3.d <- glmmTMB(Spores.per.ml ~ 1, 
                  ziformula = ~1, 
                  family = nbinom2, 
                  data = hspor)
summary(mdl3.d)

AIC(mdl3)
AIC(mdl3.b)
AIC(mdl3.c)
AIC(mdl3.d)

# Low initial dose treatment analyses
# These results are reported in Table S2
### Full model: All fixed effects and interactions
# Density needed to be scaled to facilitate model convergence
mdl4 <- glmmTMB(Spores.per.ml ~ scale(HostDensityTreatment)*Generation, 
                ziformula = ~1, 
                family = nbinom2, 
                data = lspor)
summary(mdl4)

### Nested variant: Host density and generation as fixed effects
mdl4.b <- glmmTMB(Spores.per.ml ~ scale(HostDensityTreatment) + Generation, 
                  ziformula = ~1, 
                  family = nbinom2, 
                  data = lspor)
summary(mdl4.b)

### Nested variant: Host density as sole fixed effect
mdl4.c <- glmmTMB(Spores.per.ml ~ scale(HostDensityTreatment),
                  ziformula = ~1, 
                  family = nbinom2, 
                  data = lspor)
summary(mdl4.c)

### Null model: intercept only
mdl4.d <- glmmTMB(Spores.per.ml ~ 1,
                  ziformula = ~1, 
                  family = nbinom2, 
                  data = lspor)
summary(mdl4.d)

AIC(mdl4)
AIC(mdl4.b)
AIC(mdl4.c)
AIC(mdl4.d)

# For the interested reader, we also include a full model that
# includes both doses combined (results provided in Supplement section S2)

# Need to round spore counts out in the spore data set
spor.data$Spores.per.ml <- round(spor.data$Spores.per.ml)

# The full model combining both doses
s.mdl.full <- glmmTMB(Spores.per.ml ~ HostDensityTreatment*Generation*InitialDoseTreatment,
                     ziformula = ~1,
                     family = nbinom2,
                     data = spor.data)
car::Anova(s.mdl.full)

##############
# END OF CODE
##############