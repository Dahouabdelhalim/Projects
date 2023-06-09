### This code produces figure 3 of Stewart Merrill et al. submitted

# Required packages
library(cowplot) # for consolidating figure panels
library(dplyr) # for organizing data
library(ggplot2) # for creating figures
library(glmmTMB) # for statistical analyses
library(tidyr) # for organizing data

##############
# PREPARATION
##############

# Import data
data <- read.csv('master_SIN_field density and prevalence.csv')

# Make sampling event a factor for downstream analyses
data$fevent <- as.factor(data$SamplingEvent)

#######################
# GENERATING TIME LAGS
#######################

# We lag the predictor to see effects of predictor at t = i-1 on response at t = i
# (This is equivalent to effects of predictor at t = i on response at t = i+1)

# t = i-1 (first timescale lag)
lag1 <- data %>%
        group_by(Lake) %>%
        mutate(lagdens1 = lag(lnTotalArealDensity, 1, order_by = SamplingEvent))
lag1 <- as.data.frame(lag1)

# t = i-2 (second timescale lag)
lag2 <- data %>%
        group_by(Lake) %>%
        mutate(lagdens2 = lag(lnTotalArealDensity, 2, order_by = SamplingEvent))
lag2 <- as.data.frame(lag2)

# Checking to see whether lagged density terms are correlated
cor1 <- select(data, c('Lake', 'SamplingEvent', 'EpidemicPeriod', 'lnTotalArealDensity'))
cor2 <- select(lag1, c('Lake', 'SamplingEvent', 'EpidemicPeriod', 'lagdens1'))
cor3 <- select(lag2, c('Lake', 'SamplingEvent', 'EpidemicPeriod', 'lagdens2'))

# Merge the data sets with lagged density terms
corr.data <- merge(cor1, cor2, by = c('Lake', 'SamplingEvent'))
corr.data <- merge(corr.data, cor3, by = c('Lake', 'SamplingEvent'))

# Remove pre-epidemic data (these were not used in our analyses)
corr.data <- subset(corr.data, corr.data$EpidemicPeriod != "pre-epidemic")

# Correlation tests
cor.test(corr.data$lnTotalArealDensity, corr.data$lagdens1)
# t = i and t = i + 1 correlation = 0.207, p = 0.321
cor.test(corr.data$lagdens1, corr.data$lagdens2)
# t = i + 1 and t = i + 2 correlation = 0.297, p = 0.111

###########################################################
# REFORMATTING PREVALENCE TO INDIVIDUAL INFECTION STATUSES
###########################################################

# To get the data in long form (one row for each individual)
# we pivot our raw data from prevalences to individual infection statuses
stat1 <- select(data, c('Lake', 'SamplingEvent', 'InfectedN', 'UninfectedN'))
stat1 <- stat1 %>%
                pivot_longer(c(InfectedN, UninfectedN), names_to = "Status", values_to = "Count") %>%
                drop_na(Count) %>%
                uncount(Count)

# Recode "Uninfected" as 0 and "Infected" as 1
stat1$Status <- ifelse(stat1$Status=="UninfectedN",0,ifelse(stat1$Status=="InfectedN",1,NA))

# Now merge our three datasets (t = i data, lag1, and lag2) to the status data

# First starting with our t = i dataset:
long.data.supp <- select(data, c('Lake', 
                                 'SamplingEvent',
                                 'lnTotalArealDensity',
                                 'EpidemicPeriod',
                                 'fevent',
                                 'Prevalence')) # Select only columns of interest for analyses
data.b <- merge(stat1, long.data.supp, by = c('Lake', 'SamplingEvent'), all.x = T, all.y = F)
data.b <- subset(data.b, data.b$EpidemicPeriod != "pre-epidemic") # Remove pre-epidemic data prior to analyses

# First lag at t = i - 1 dataset:
long.data.supp.2 <- select(lag1, c('Lake', 
                                 'SamplingEvent',
                                 'lagdens1',
                                 'EpidemicPeriod',
                                 'fevent',
                                 'Prevalence')) # Select only columns of interest for analyses
lag1.b <- merge(stat1, long.data.supp.2, by = c('Lake', 'SamplingEvent'), all.x = T, all.y = F)
lag1.b <- subset(lag1.b, lag1.b$EpidemicPeriod != "pre-epidemic") # Remove pre-epidemic data prior to analyses

# Second lag at t = i - 2 dataset:
long.data.supp.3 <- select(lag2, c('Lake', 
                                   'SamplingEvent',
                                   'lagdens2',
                                   'EpidemicPeriod',
                                   'fevent',
                                   'Prevalence')) # Select only columns of interest for analyses
lag2.b <- merge(stat1, long.data.supp.3, by = c('Lake', 'SamplingEvent'), all.x = T, all.y = F)
lag2.b <- subset(lag2.b, lag2.b$EpidemicPeriod != "pre-epidemic") # Remove pre-epidemic data prior to analyses

###########
# ANALYSES
###########

# Prevalence by host density: t = i
mdl1 <- glmmTMB(Status ~ lnTotalArealDensity + 
                ar1(fevent + 0 | Lake), 
                family = binomial, 
                data = data.b)
summary(mdl1)

# Prevalence by host density: t = i + 1
mdl2 <- glmmTMB(Status ~ lagdens1 + 
                ar1(fevent + 0 | Lake), 
                family = binomial, 
                data = lag1.b)
summary(mdl2)

# Prevalence by host density: t = i + 2
mdl3 <- glmmTMB(Status ~ lagdens2 + 
                ar1(fevent + 0 | Lake), 
                family = binomial, 
                data = lag2.b)
summary(mdl3)

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
                 axis.title.x = element_text(family="Arial", size=13, color="black", margin=margin(12,0,0,0)),
                 axis.title.y = element_text(family="Arial", size=13, color="black", margin=margin(0,12,0,0)),
                 axis.text.x = element_text(family="Arial", size=10, color="black"),
                 axis.text.y = element_text(family="Arial", size=10, color="black"),
                 strip.background = element_blank())

# Remove pre-epidemic points from the figure
data <- subset(data, data$EpidemicPeriod != "pre-epidemic")
lag1 <- subset(lag1, lag1$EpidemicPeriod != "pre-epidemic")
lag2 <- subset(lag2, lag2$EpidemicPeriod != "pre-epidemic")

Fig4a <- ggplot(data = data, aes(x = lnTotalArealDensity, y = Prevalence)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", fullrange = T, col = "black") +
        xlim(8.3, 12.7) +
        xlab("ln(Total host density)") +
        ylab("Host risk\\n(Infection prevalence") +
        presets

Fig4b <- ggplot(data = lag1, aes(x = lagdens1, y = Prevalence)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", fullrange = T, col = "black") +
        xlim(8.3, 12.7) +
        xlab("ln(Total host density)") +
        ylab("Host risk\\n(Infection prevalence") +
        presets

Fig4c <- ggplot(data = lag2, aes(x = lagdens2, y = Prevalence)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", fullrange = T, col = "black") +
        xlim(8.3, 12.7) +
        xlab("ln(Total host density)") +
        ylab("Host risk\\n(Infection prevalence") +
        presets

plot_grid(Fig4a, Fig4b, Fig4c, ncol = 3)

##############
# END OF CODE
##############