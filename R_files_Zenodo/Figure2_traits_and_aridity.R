### SEED SIZE (FIGURE 2A) ###

### 01. import all seed size results from ImageJ particle analysis -------------
setwd("~/Desktop/Dryad/")
dat <- read.csv(file = "2A_seed_size_data.csv", header = T, stringsAsFactors = F)

nrow(dat) # 33093 seeds
length(unique(dat$maternal_line)) # 291 maternal lines
length(unique(dat$population)) # 12 populations


### 02. assess seed size size number correlation -------------------------------
cor_dat <- data.frame(
    maternal_line = names(tapply(dat$seed_area, dat$maternal_line, mean)),
    avg_seed_size = tapply(dat$seed_area, dat$maternal_line, mean),
    n_seeds = as.vector(table(dat$maternal_line)))

cor_dat$population <- sapply(strsplit(cor_dat$maternal_line, split = "_"), function(x) x[1])

cor.test(cor_dat$avg_seed_size, cor_dat$n_seeds) 
    # r = -0.13, p = 0.017, weak correlation between seed size and seed number

plot(avg_seed_size ~ n_seeds, col = "gray", pch = 16, data = cor_dat)

# take into account the effect of population
library(lmerTest)

m0 <- lmer(avg_seed_size ~ n_seeds + (1|population), data = cor_dat)
summary(m0) # p = 0.699


### 03. import climate data ----------------------------------------------------
clim_dat <- read.csv(file = "Plantago_climate_data_simple.csv", header = T, stringsAsFactors = F)

dat$population[dat$population == "Huckaby"] <- "HuckabyTrailhead"

# remove underscores from clim_dat$Population
clim_dat$Population <-  gsub(pattern = "_", replacement = "", x = clim_dat$Population)

sort(unique(dat$population)) %in% clim_dat$Population

# merge data
all_dat <- merge(dat, clim_dat, by.x = "population", by.y = "Population")


### 05. model seed size ~ aridity of source climate ----------------------------
library(ggplot2)
library(ggeffects)
library(scales)

plot(seed_area ~ jitter(AI, factor = 5), data = all_dat,
     pch = 16, cex = 0.5, col = alpha("gray", alpha = 0.2))

# Aridity Index
m1.mixed.AI <- lmer(seed_area ~ AI + (1|population) + (1|maternal_line), 
                data = all_dat)

summary(m1.mixed.AI)

AI_df <- ggpredict(m1.mixed.AI, terms = c("AI"))

p1 <- ggplot(AI_df, aes(x, predicted)) + 
    geom_jitter(data = all_dat, aes(x = AI, y = seed_area), 
        width = 0.125, pch = 16, col = "darkgray", alpha = 0.025) +
    geom_line(aes(linetype=group, color=group), size = 1) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha = 0.1) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    
    theme_light() +
    xlab("") +
    theme(axis.text.x = element_blank()) +
    ylab(expression(bold("Seed area (mm"^2*")"))) +
    theme(axis.title.y = element_text(size = 10, face = "bold")) +
    
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") 


# Moisture Index
m1.mixed.MI <- lmer(seed_area ~ MI + (1|population) + (1|maternal_line), 
                    data = all_dat)

summary(m1.mixed.MI)

MI_df <- ggpredict(m1.mixed.MI, terms = c("MI"))

ggplot(MI_df, aes(x, predicted)) + 
    geom_line(aes(linetype=group, color=group), size = 1) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha = 0.05) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    
    geom_jitter(data = all_dat, aes(x = MI, y = seed_area), width = 0.0015, pch = 16, col = "darkgray", alpha = 0.05) +
    
    theme_light() +
    xlab("Climatic Moisture Index") +
    theme(axis.title.x = element_text(size = 12, face = "bold")) +
    ylab("Seed Area (mm^2)") +
    theme(axis.title.y = element_text(size = 12, face = "bold")) +
    
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") 



### GERMINATION TIMING VARIATION (FIGURE 2B) ###
### 00. import data ------------------------------------------------------------
germ_dat <- read.csv(
    file = "2B_germination_data.csv",
    header = TRUE,
    stringsAsFactors = FALSE)

# remove doubled planted and unplanted cells
germ_dat2 <- subset(germ_dat, Remove != "Remove")

# remove herbarium populations
germ_dat3 <- germ_dat2[!grepl(patter = "herb", x = germ_dat2$Population),]
nrow(germ_dat3) # 1726 seeds planted in germination experiment

# merge with climate data

# remove underscores from population name
germ_dat3$Population <- gsub(pattern = "_", replacement = "", x = germ_dat3$Population)

g_dat <- merge(germ_dat3, clim_dat, by.x = "Population", by.y = "Population")
nrow(g_dat) # 1726


### 01. explore germination timing in relation to climatic variables -----------

# add days to germination
g_dat$days_to_germ <- g_dat$Jdate_germ - 151 # seeds planted/begin watering on May 31st.


# Aridity Index
m1.germ.mixed.AI <- lmer(days_to_germ ~ AI + 
    (1|Population) + (1|tray) + (1|Envelope_ID), 
    data = g_dat)

summary(m1.germ.mixed.AI)

germ_AI_df <- ggpredict(m1.germ.mixed.AI, terms = c("AI"))

p2 <- ggplot(germ_AI_df, aes(x, predicted)) + 
    geom_jitter(data = g_dat, aes(x = AI, y = days_to_germ), size = 1.5,
                width = 0.125, pch = 16, col = "darkgray", alpha = 0.25) +
    geom_line(aes(linetype=group, color=group), size = 1) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha = 0.25) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +

    theme_light() +
    xlab("") +
    theme(axis.text.x = element_blank()) +
    ylab("Days to germinate") +
    theme(axis.title.y = element_text(size = 10, face = "bold")) +
    
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") 


# Climate Moisture Index
m1.germ.mixed.MI <- lmer(days_to_germ ~ MI + 
                             (1|Population) + (1|tray) + (1|Envelope_ID), 
                         data = g_dat)

summary(m1.germ.mixed.MI)

germ_MI_df <- ggpredict(m1.germ.mixed.MI, terms = c("MI"))

ggplot(germ_MI_df, aes(x, predicted)) + 
    geom_line(aes(linetype=group, color=group), size = 1) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha = 0.25) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    
    geom_jitter(data = g_dat, aes(x = MI, y = days_to_germ), size = 3,
                width = 0.001, pch = 16, col = "darkgray", alpha = 0.75) +
    
    theme_light() +
    xlab("Aridity Index") +
    theme(axis.title.x = element_text(size = 12, face = "bold")) +
    ylab("Days to germinate") +
    theme(axis.title.y = element_text(size = 12, face = "bold")) +
    
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") 


### SPECIFIC LEAF AREA VARIATION (FIGURE 2C) ###

### 00. import SLA data --------------------------------------------------------
SLA_dat <- read.csv(file = "2C_SLA_data.csv", header = T, stringsAsFactors = F)

# add SLA column to data.frame
SLA_dat$SLA <- SLA_dat$leaf_area_mm_squared / SLA_dat$leaf_weight_mg

table(!is.na(SLA_dat$SLA)) # 470 total individuals

# remove resurrection plants
l_dat <- SLA_dat[!grepl(pattern = "herb", x = SLA_dat$Population), ]

table(!is.na(l_dat$SLA)) # 307 "contemporary" individuals

# merge with climate data
sort(unique(l_dat$Population))

l_dat$Population <- gsub(pattern = "_", replacement = "", x = l_dat$Population)

l_dat <- merge(l_dat, clim_dat, by.x = "Population", by.y = "Population")


# Aridity Index
m1.sla.mixed.AI <- lmer(SLA ~ AI + 
    (1|Population) + (1|tray) + (1|Envelope_ID), 
    data = l_dat)

summary(m1.sla.mixed.AI)

sla_AI_df <- ggpredict(m1.sla.mixed.AI, terms = c("AI"))

p3 <- ggplot(sla_AI_df, aes(x, predicted)) + 
    geom_jitter(data = l_dat, aes(x = AI, y = SLA), size = 1.5,
                width = 0.125, pch = 16, col = "darkgray", alpha = 0.5) +
    geom_line(aes(linetype=group, color=group), size = 1) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha = 0.25) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    
    theme_light() +

    scale_y_continuous(breaks = c(10,20,30), labels = c(10,20,30), limits = c(10,30)) +
    xlab("Aridity Index") +
    theme(axis.title.x = element_text(size = 10, face = "bold")) +
    ylab("Specific leaf area") +
    theme(axis.title.y = element_text(size = 10, face = "bold")) +
    
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") 


# Moisture Index
m1.sla.mixed.MI <- lmer(SLA ~ MI + 
    (1|Population) + (1|tray) + (1|Envelope_ID), 
    data = l_dat)

summary(m1.sla.mixed.MI)

sla_MI_df <- ggpredict(m1.sla.mixed.MI, terms = c("MI"))

ggplot(sla_MI_df, aes(x, predicted)) + 
    geom_line(aes(linetype=group, color=group), size = 1) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha = 0.25) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    
    geom_jitter(data = l_dat, aes(x = MI, y = SLA), size = 3,
                width = 0.0001, pch = 16, col = "darkgray", alpha = 0.75) +
    
    theme_light() +
    xlab("Aridity Index") +
    theme(axis.title.x = element_text(size = 12, face = "bold")) +
    ylab("Specific Leaf Area") +
    theme(axis.title.y = element_text(size = 12, face = "bold")) +
    
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") 


### 04. FIGURE 2 (build multi-panel plot) --------------------------------------
source("multiplot.R")

dev.new()
multiplot(p1,p2,p3, cols = 1)

# quartz.save(file = "Figure2.jpg", type = "jpg", dpi = 300)

