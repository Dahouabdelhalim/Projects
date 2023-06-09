### Analysis script for:

# "Asymmetrical reproductive barriers in sympatric Jewelflowers: 
# are floral isolation, genetic incompatibilities, 
# and floral trait displacement connected?â€� (Manuscript ID: BJLS-6658)

# K. Christie, J.P. Doan, W.C. McBride, S.Y. Strauss 2021

# Biological Journal of the Linnean Society; (Manuscript ID: BJLS-6658)


### Floral rewards analysis


### 00. import pollen production data ------------------------------------------
setwd("~/Desktop/Dryad_data/")
dat = read.csv(file = "3_pollen_production_data.csv", header = T, stringsAsFactors = F)


### 01. summarize data ---------------------------------------------------------
nrow(dat) # 81 measurements
length(unique(dat$GH_cone_ID)) # 41 individuals
length(unique(dat$species_by_site)) # 10 sites

# add geography code
dat$species_full <- paste(dat$species, dat$regional_geography, sep = "-")

# add a species X geography grouping
dat$species_group <- paste(dat$species, dat$regional_geography, sep = "-")


### 02. plot data --------------------------------------------------------------
library("scales")

dev.new()
par(mfrow = c(2,1))
par(font.lab = 2, cex.lab = 1.25, mar = c(1,5,1,5), cex.names = 1.5)

# boxplot
boxplot(dat$total_pollen_grains/1000 ~ dat$species_full,
    xlab = "", ylab = "Pollen grains per flower (x1000)",
    outline = F, ylim = c(25,300),
    names = c("","",""))

# add raw data
brew_allo <- subset(dat, species_group == "S_breweri-allo")
for(i in 1:nrow(brew_allo)){
    points(x = jitter(1, factor = 5), y = brew_allo$total_pollen_grains[i]/1000,
    col = alpha("blue", alpha = 0.65),
    cex = 1.5, pch = 16)
}

brew_sym <- subset(dat, species_group == "S_breweri-sym")
for(i in 1:nrow(brew_sym)){
    points(x = jitter(2, factor = 3), y = brew_sym$total_pollen_grains[i]/1000,
    col = alpha("blue", alpha = 0.65),
    cex = 1.5, pch = 16)
}

hesp_sym <- subset(dat, species_group == "S_hesperidis-sym")
for(i in 1:nrow(hesp_sym)){
    points(x = jitter(3, factor = 2), y = hesp_sym$total_pollen_grains[i]/1000,
    col = alpha("darkgreen", alpha = 0.65),
    cex = 1.5, pch = 16)
}


### 03 build mixed model -------------------------------------------------------
library("lme4")
library("lmerTest")

tapply(dat$total_pollen_grains, dat$species_group, mean)
tapply(dat$total_pollen_grains, dat$species_group, sd)

m1_mixed <- lmer(total_pollen_grains ~ species + species_group + (1|GH_cone_ID), data = dat)
summary(m1_mixed)

anova(m1_mixed)

# conduct Tukey's
library("emmeans")
emmeans(m1_mixed, list(pairwise ~ species_group), adjust = "tukey")

# add Tukey's groupings from the plot
text(1, 285, "A", font = 2, cex = 1.5)
text(2, 285, "A", font = 2, cex = 1.5)
text(3, 285, "A", font = 2, cex = 1.5)


### 04. import nectar sugar data -----------------------------------------------
dat  <-  read.csv(file = "3_nectar_sugar_data.csv", header = TRUE, stringsAsFactors = FALSE)

# summarize data
nrow(dat) # samples from 63 different individuals
table(dat$species) # 43 breweri, 20 hesperidis

table(dat$species, dat$location)
# 24 breweri from McL (sym); 19 from regional populations (allo)
# 20 hesperidis from McL (sym)


### 05. plot of nectar sugar concentration for three groups --------------------
    # allopatric S. breweri
    # sympatric S. breweri
    # S. hesperidis (all sympatric for this study)

dat$species_group <- paste(dat$species, dat$regional_geography, sep = "-")

par(mar = c(5,5,1,5))
boxplot(dat$brix ~ dat$species_group, xlab = "", 
        names = c("S. breweri (allo)", "S. breweri (sym)", "S. hesperidis"),
        ylab = "Nectar sugar (% sucrose)",
        ylim = c(10,50))

# add raw data
brew_allo <- subset(dat, species_group == "breweri-allopatric")
for(i in 1:nrow(brew_allo)){
    points(x = jitter(1, factor = 5), y = brew_allo$brix[i],
           col = alpha("blue", alpha = 0.65),
           cex = 1.5, pch = 16)
}

brew_sym <- subset(dat, species_group == "breweri-sympatric")
for(i in 1:nrow(brew_sym)){
    points(x = jitter(2, factor = 3), y = brew_sym$brix[i],
           col = alpha("blue", alpha = 0.65),
           cex = 1.5, pch = 16)
}

hesp_sym <- subset(dat, species_group == "hesperidis-sympatric")
for(i in 1:nrow(hesp_sym)){
    points(x = jitter(3, factor = 2), y = hesp_sym$brix[i],
           col = alpha("darkgreen", alpha = 0.65),
           cex = 1.5, pch = 16)
}


### 06. build models -----------------------------------------------------------
tapply(dat$brix, dat$species_group, mean)

m2 <- lm(brix ~ species + species_group, data = dat)
summary(m2)
anova(m2)

# conduct Tukey's
emmeans(m2, list(pairwise ~ species_group), adjust = "tukey")

# add Tukey's groupings from the plot
text(1, 49, "A", font = 2, cex = 1.5)
text(2, 49, "B", font = 2, cex = 1.5)
text(3, 49, "A", font = 2, cex = 1.5)


### 07. summarize nectar sugar concentration at focal field sites --------------
dat2 <- subset(dat, site == "Quarry_View" | site == "Napa_Junction")
tapply(dat2$brix, dat2$species, mean)
    # S. breweri - 31.8% (correct number in text)
    # S. hesperidis - 22.0% (correct number in text)


### 08. test for relationship between latitude and nectar sugar concentration --
breweri <- subset(dat, species == "breweri")
sites <- read.csv(file = "3_population_locations.csv", header = T, stringsAsFactors = F)
brew2 <- merge(breweri, sites, by.x = "site", by.y = "site")

# remove McLaughlin sites for the regional analysis
brew_regional = subset(brew2, 
        site != "Coyote_Hill" &
        site != "Napa_Junction" &
        site != "Quarry_View" &
        site != "Lomatium_Hill" &
        site != "Slide_Canyon")

m3  <-  lm(brew_regional$brix ~ brew_regional$lat)
summary(m3) # no relationship between latitude and nectar sugar based on northern sites
    # p = 0.27 (correct p-value in text)






