
# Exotic herbivores dominate Australian high-elevation grasslands
# https://doi.org/10.1111/csp2.601
# Data analyses and plots
# Renee Hartley
# 12/11/2021


# Data preparation

# Load the packages required for Bayesian analysis and grouping data
library(brms)
library(dplyr)
library(DHARMa)
library(glmmTMB)
library(ggplot2)
library(emmeans)
library(splines)
library(mgcv)
library(scales)
library(car)
library(lmerTest)
library(cowplot)
library(stats)
library(tidyverse)

# prepare data
GIS <- read.csv("Hartley_etal_CSP2-601.csv")

#Distribution of herbivore families
E <- subset(GIS, GIS$TinvC > 0) # 94.0% of sites with exotic herbivores present
N <- subset(GIS, GIS$TnatC > 0) # 32.8% of sites with native herbivores present
Z <- subset(GIS, GIS$TinvC == 0 & GIS$TnatC == 0) # one site where neither exotic nor native herbivores are present
NO <- subset(GIS, GIS$TinvC == 0 & GIS$TnatC > 0) # 3 sites where only native herbivores are present
EO <- subset(GIS, GIS$TinvC > 0 & GIS$TnatC == 0) # 44 sites where only exotic herbivores are present
CO <- subset(GIS, GIS$TinvC > 0 & GIS$TnatC > 0) # 19 sites where both exotic and native herbivores are present

# tabulate contribution of herbivore families
knitr::kable(x=GIS %>% summarise(mean(TtgpC), mean(TgenC), mean(ThorC), mean(TdeeC), mean(TpigC), 
                                 mean(TlepC), mean(TmacC), mean(TwomC), mean(TBTRscat)), digits = 2,
             caption = "Table 2: Contribution of each herbivore family to total herbivore activity.")


# convert % cover variables to proportions and add small value to variables with zeros
GIS$forbProp <- (GIS$forb.mean/100)
GIS$weedProp <- (GIS$weed.mean/100)+0.0001
GIS$baregroundProp <- (GIS$bareground.mean/100)+0.0001
GIS$grassProp <- (GIS$grass.mean/100)-0.001 #subtracting small value as variable contains 100% grass cover and no zeros
GIS$shrubProp <- (GIS$shrub.mean/100)+0.0001
GIS$sedrusProp <- (GIS$sedge.rush.mean/100)+0.0001

# Log-transform the predictor variables. Add small value (1) to enable log-transformation of zeros.
GIS$tgpCL <- log(GIS$TtgpC+1)
GIS$horCL <- log(GIS$ThorC+1)
GIS$deeCL <- log(GIS$TdeeC+1)
GIS$lepCL <- log(GIS$TlepC+1)
GIS$macCL <- log(GIS$TmacC+1)

# Add small value (== half of the lowest value in the variable) to distance measures with zero values to enable log-transformation of zeros.
GIS$DistWood.min[GIS$DistWood.min==0] <- 0.5
GIS$DistWtr.min[GIS$DistWtr.min==0] <- 1

# Log transform distance variables
GIS$DistRdL <- log(GIS$DistRd.min)
GIS$DistWoodL <- log(GIS$DistWood.min)
GIS$DistWtrL <- log(GIS$DistWtr.min)

# Scale the continuous predictor variables. 
GIS$eleS <- scale(GIS$ele.mean)[,1]
GIS$DistRdS <- scale(GIS$DistRdL)[,1]
GIS$DistWoodS <- scale(GIS$DistWoodL)[,1]
GIS$DistWtrS <- scale(GIS$DistWtrL)[,1]
GIS$TWIS <- scale(GIS$TWImean)[,1]
GIS$tgpS <- scale(GIS$tgpCL)[,1]
GIS$horS <- scale(GIS$horCL)[,1]
GIS$deeS <- scale(GIS$deeCL)[,1]
GIS$lepS <- scale(GIS$lepCL)[,1]
# GIS$pigS <- scale(GIS$pigCL)[,1] # excluded from analyses, too few data
GIS$macS <- scale(GIS$macCL)[,1]
# GIS$womS <- scale(GIS$womCL)[,1] #  excluded from analyses, too few data

# create new dataframe for soil variables, removing NA values 
soils <- subset(GIS, GIS$soil.mean!="NA")

#in preparation for following tables (Note 'grassheight' refers to height of groundstorey regardless of plant order/form intercepted):
grassheight.noNA <- subset(GIS, GIS$Site != c("JTJ01")) # remove sites with NA values from GIS df

sep <- GIS %>% group_by(vegtype) %>% summarise(mean(TtgpC), mean(ThorC), mean(TdeeC), mean(TpigC), 
                                              mean(TlepC), mean(TmacC), mean(TwomC), n())
csep <- GIS %>% group_by(vegtype) %>% summarise(sum(tgpPA), sum(horPA), sum(deePA), sum(pigPA), 
                                               sum(lepPA), sum(macPA), sum(womPA), n())
vegtypeveg <- grassheight.noNA %>% group_by(vegtype) %>% summarise(mean(grassheight.mean),mean(fol5.mean), mean(sedge.rush.mean), mean(shrub.mean), n())

envcor <- with(grassheight.noNA,
               cor(cbind(grassheight.mean, eleS, DistWtrS, 
                         DistRdS, DistWoodS, TWIS),use = "everything"))

# Check Spearman and Pearson coefficients 
vdistw <- ggplot(grassheight.noNA, aes(grassheight.mean, DistWoodS)) + geom_point() + labs(x="Mean vegetation height (cm)", y="Log-distance to woodland") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
vdistt <- ggplot(grassheight.noNA, aes(grassheight.mean, DistWtrS)) + geom_point() + labs(x="Mean vegetation height (cm)", y="Log-distance to waterbody") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
vdistr <- ggplot(grassheight.noNA, aes(grassheight.mean, DistRdS)) + geom_point() + labs(x="Mean vegetation height (cm)", y="Log-distance to roads") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
vele  <- ggplot(grassheight.noNA, aes(grassheight.mean, eleS)) + geom_point()  + labs(x="Mean vegetation height (cm)", y="Elevation") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
vTWI <- ggplot(grassheight.noNA, aes(grassheight.mean, TWIS)) + geom_point()  + labs(x="Mean vegetation height (cm)", y="Topographic Wetness Index") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
vgr <- ggplot(grassheight.noNA, aes(grassheight.mean, vegtype)) + geom_point()  + labs(x="Mean vegetation height (cm)", y="Vegetation type") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
#vlong <- ggplot(grassheight.noNA, aes(grassheight.mean, longGDA94)) + geom_point()  + labs(x="Mean vegetation height (mm)", y="Longitude") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
#vlat <- ggplot(grassheight.noNA, aes(grassheight.mean, latGDA94)) + geom_point() + labs(x="Mean vegetation height (mm)", y="Latitude") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
fdistw <- ggplot(grassheight.noNA, aes(fol5.mean, DistWoodS)) + geom_point() + labs(x="Mean foliage density", y="Log-distance to woodland") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
fdistt <- ggplot(grassheight.noNA, aes(fol5.mean, DistWtrS)) + geom_point() + labs(x="Mean foliage density", y="Log-distance to waterbody") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
fdistr <- ggplot(grassheight.noNA, aes(fol5.mean, DistRdS)) + geom_point() + labs(x="Mean foliage density", y="Log-distance to roads") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
fele  <- ggplot(grassheight.noNA, aes(fol5.mean, eleS)) + geom_point()  + labs(x="Mean foliage density", y="Elevation") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
fTWI <- ggplot(grassheight.noNA, aes(fol5.mean, TWIS)) + geom_point()  + labs(x="Mean foliage density", y="Topographic Wetness Index") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
fgr <- ggplot(grassheight.noNA, aes(fol5.mean, vegtype)) + geom_point()  + labs(x="Mean foliage density", y="Vegetation type") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
odistw <- ggplot(grassheight.noNA, aes(forb.mean, DistWoodS)) + geom_point() + labs(x="Mean proportion of forbs", y="Log-distance to woodland") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
odistt <- ggplot(grassheight.noNA, aes(forb.mean, DistWtrS)) + geom_point() + labs(x="Mean proportion of forbs", y="Log-distance to waterbody") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
odistr <- ggplot(grassheight.noNA, aes(forb.mean, DistRdS)) + geom_point() + labs(x="Mean proportion of forbs", y="Log-distance to roads") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
oele  <- ggplot(grassheight.noNA, aes(forb.mean, eleS)) + geom_point()  + labs(x="Mean proportion of forbs", y="Elevation") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
oTWI <- ggplot(grassheight.noNA, aes(forb.mean, TWIS)) + geom_point()  + labs(x="Mean proportion of forbs", y="Topographic Wetness Index") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
ogr <- ggplot(grassheight.noNA, aes(forb.mean, vegtype)) + geom_point()  + labs(x="Mean proportion of forbs", y="Vegetation type") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
wdistw <- ggplot(grassheight.noNA, aes(weed.mean, DistWoodS)) + geom_point() + labs(x="Mean proportion of weeds", y="Log-distance to woodland") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
wdistt <- ggplot(grassheight.noNA, aes(weed.mean, DistWtrS)) + geom_point() + labs(x="Mean proportion of weeds", y="Log-distance to waterbody") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
wdistr <- ggplot(grassheight.noNA, aes(weed.mean, DistRdS)) + geom_point() + labs(x="Mean proportion of weeds", y="Log-distance to roads") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
wele  <- ggplot(grassheight.noNA, aes(weed.mean, eleS)) + geom_point()  + labs(x="Mean proportion of weeds", y="Elevation") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
wTWI <- ggplot(grassheight.noNA, aes(weed.mean, TWIS)) + geom_point()  + labs(x="Mean proportion of weeds", y="Topographic Wetness Index") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
wgr <- ggplot(grassheight.noNA, aes(weed.mean, vegtype)) + geom_point()  + labs(x="Mean proportion of weeds", y="Vegetation type") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
bdistw <- ggplot(grassheight.noNA, aes(bareground.mean, DistWoodS)) + geom_point() + labs(x="Mean proportion of bare ground", y="Log-distance to woodland") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
bdistt <- ggplot(grassheight.noNA, aes(bareground.mean, DistWtrS)) + geom_point() + labs(x="Mean proportion of bare ground", y="Log-distance to waterbody") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
bdistr <- ggplot(grassheight.noNA, aes(bareground.mean, DistRdS)) + geom_point() + labs(x="Mean proportion of bare ground", y="Log-distance to roads") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
bele  <- ggplot(grassheight.noNA, aes(bareground.mean, eleS)) + geom_point()  + labs(x="Mean proportion of bare ground", y="Elevation") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
bTWI <- ggplot(grassheight.noNA, aes(bareground.mean, TWIS)) + geom_point()  + labs(x="Mean proportion of bare ground", y="Topographic Wetness Index") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
bgr <- ggplot(grassheight.noNA, aes(bareground.mean, vegtype)) + geom_point()  + labs(x="Mean proportion of bare ground", y="Vegetation type") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
sdistw <- ggplot(soils, aes(soil.mean, DistWoodS)) + geom_point() + labs(x="Soil compactions (kPa)", y="Log-distance to woodland") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
sdistt <- ggplot(soils, aes(soil.mean, DistWtrS)) + geom_point() + labs(x="Soil compactions (kPa)", y="Log-distance to waterbody") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
sdistr <- ggplot(soils, aes(soil.mean, DistRdS)) + geom_point() + labs(x="Soil compactions (kPa)", y="Log-distance to roads") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
sele  <- ggplot(soils, aes(soil.mean, eleS)) + geom_point()  + labs(x="Soil compactions (kPa)", y="Elevation") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
sTWI <- ggplot(soils, aes(soil.mean, TWIS)) + geom_point()  + labs(x="Soil compactions (kPa)", y="Topographic Wetness Index") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))
sgr <- ggplot(soils, aes(soil.mean, vegtype)) + geom_point()  + labs(x="Soil compactions (kPa)", y="Vegetation type") + theme(text = element_text(size = 8), axis.title.x = element_text(hjust=1))


#check correlation and significance of vegetation height in relation to landscape features
cor.test(grassheight.noNA$grassheight.mean, grassheight.noNA$DistWoodS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$grassheight.mean, grassheight.noNA$DistWtrS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$grassheight.mean, grassheight.noNA$DistRdS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$grassheight.mean, grassheight.noNA$eleS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$grassheight.mean, grassheight.noNA$TWIS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of foliage density in relation to landscape features
cor.test(grassheight.noNA$fol5.mean, grassheight.noNA$DistWoodS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$fol5.mean, grassheight.noNA$DistWtrS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$fol5.mean, grassheight.noNA$DistRdS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$fol5.mean, grassheight.noNA$eleS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$fol5.mean, grassheight.noNA$TWIS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of forb cover in relation to landscape features
cor.test(grassheight.noNA$forb.mean, grassheight.noNA$DistWoodS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$forb.mean, grassheight.noNA$DistWtrS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$forb.mean, grassheight.noNA$DistRdS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$forb.mean, grassheight.noNA$eleS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$forb.mean, grassheight.noNA$TWIS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of weed cover in relation to landscape features
cor.test(grassheight.noNA$weed.mean, grassheight.noNA$DistWoodS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$weed.mean, grassheight.noNA$DistWtrS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$weed.mean, grassheight.noNA$DistRdS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$weed.mean, grassheight.noNA$eleS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$weed.mean, grassheight.noNA$TWIS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of bare ground cover in relation to landscape features
cor.test(grassheight.noNA$bareground.mean, grassheight.noNA$DistWoodS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$bareground.mean, grassheight.noNA$DistWtrS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$bareground.mean, grassheight.noNA$DistRdS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$bareground.mean, grassheight.noNA$eleS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$bareground.mean, grassheight.noNA$TWIS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of soil compaction in relation to landscape features
cor.test(soils$soil.mean, soils$DistWoodS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(soils$soil.mean, soils$DistWtrS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(soils$soil.mean, soils$DistRdS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(soils$soil.mean, soils$eleS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(soils$soil.mean, soils$TWIS, 
         alternative = "two.sided",
         method = "pearson", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of vegetation height in relation to landscape features
cor.test(grassheight.noNA$grassheight.mean, grassheight.noNA$DistWoodS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$grassheight.mean, grassheight.noNA$DistWtrS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$grassheight.mean, grassheight.noNA$DistRdS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$grassheight.mean, grassheight.noNA$eleS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$grassheight.mean, grassheight.noNA$TWIS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of foliage density in relation to landscape features
cor.test(grassheight.noNA$fol5.mean, grassheight.noNA$DistWoodS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$fol5.mean, grassheight.noNA$DistWtrS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$fol5.mean, grassheight.noNA$DistRdS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$fol5.mean, grassheight.noNA$eleS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$fol5.mean, grassheight.noNA$TWIS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of forb cover in relation to landscape features
cor.test(grassheight.noNA$forb.mean, grassheight.noNA$DistWoodS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$forb.mean, grassheight.noNA$DistWtrS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$forb.mean, grassheight.noNA$DistRdS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$forb.mean, grassheight.noNA$eleS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$forb.mean, grassheight.noNA$TWIS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of weed cover in relation to landscape features
cor.test(grassheight.noNA$weed.mean, grassheight.noNA$DistWoodS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$weed.mean, grassheight.noNA$DistWtrS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$weed.mean, grassheight.noNA$DistRdS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$weed.mean, grassheight.noNA$eleS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$weed.mean, grassheight.noNA$TWIS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of bare ground cover in relation to landscape features
cor.test(grassheight.noNA$bareground.mean, grassheight.noNA$DistWoodS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$bareground.mean, grassheight.noNA$DistWtrS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$bareground.mean, grassheight.noNA$DistRdS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$bareground.mean, grassheight.noNA$eleS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(grassheight.noNA$bareground.mean, grassheight.noNA$TWIS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)

#check correlation and significance of soil compaction in relation to landscape features
cor.test(soils$soil.mean, soils$DistWoodS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(soils$soil.mean, soils$DistWtrS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(soils$soil.mean, soils$DistRdS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(soils$soil.mean, soils$eleS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)
cor.test(soils$soil.mean, soils$TWIS, 
         alternative = "two.sided",
         method = "spearman", exact = NULL, conf.level = 0.95,
         continuity = FALSE)


knitr::kable(x=GIS %>% summarise(sum(tgpPA), sum(horPA), sum(deePA), sum(pigPA), 
                                 sum(lepPA), sum(macPA), sum(womPA), sum(BTRPA), n()), 
             caption = "Table 1: Number of sites at which each herbivore group is present.")

knitr::kable(x=summary(as.factor(GIS$vegtype)), caption = "Table 2: Number of sites within each grassland type.")

knitr::kable(x=csep, digits = 2, caption = 
               "Table 3: Number of sites with herbivore sign present for each grassland type.")

knitr::kable(x=sep, digits = 2, caption = 
               "Table 4: Mean grazing sign frequency for each grassland type.")

# look for potential correlation between plant height and elevation
knitr::kable(x=envcor, digits = 2, caption = 
               "Table 5: Correlation matrix of ground cover vegetation height and environmental variables.")

# test with correlation matrix
sitecor <- with(GIS,
                cor(cbind(eleS, DistWtrS, DistWoodS, 
                          DistRdS, TWIS),use = "everything"))
# test with correlation matrix
grazcor <- with(GIS,
                cor(cbind(tgpS, horS, deeS, 
                          lepS, macS),use = "everything"))

# horse and total herbivore activity highly correlated (0.8). Therefore model separately
knitr::kable(x=grazcor, digits = 2, caption = 
               "Table 6: Correlation matrix of herbivore activity according to herbivore group.")

# test with correlation matrix
vegcor <- with(soils,
               cor(cbind(grassheight.mean, fol5.mean, forbProp, weedProp, baregroundProp,
                         soil.mean),use = "everything"))
# Vegetation height and foliar density highly correlated (0.91)
knitr::kable(x=vegcor, digits = 2, caption = 
               "Table 7: Correlation matrix of vegetation and soil measures.")

knitr::kable(x=read.csv("C:/Users/email/Documents/ANU_PhD/Data/Analyses&Datasets/GrazImpact/gr_spp_cooccur_tbl.csv"), digits = 3, caption = 
               "Table 8: Co-occurrence of herbivore families.") # need to work how to do this in R



########## herbivore activity ########## 

# Dependent variables: herbivore activity; total, horse, deer, leporid, macropod (count)  
# Independent variables: elevation (cts), dist wood (cts), dist water (cts), dist road (cts), TWI (cts), vegetation classification (cat)

# Fit the negative binomial using brms.
bprior <- c(prior_string("student_t(7,0,2.5)", class = "b"),
            prior_string("student_t(7,0,2.5)", class = "Intercept"))

## Total herbivore activity analysis ("tgpC", all sign, as response variable) 

brmtgpC.nb <- brm(TtgpC ~  eleS + TWIS + vegtype + DistRdS + DistWtrS + DistWoodS, 
                  family =negbinomial(), data=GIS, cores=2, 
                  prior=bprior,control=list(adapt_delta = 0.9))
summary(brmtgpC.nb) 

# Use WAIC to determine whether the Poisson or negative binomial fits best. 

# Run model selection 
brmtgpC.nb.m63<-update(brmtgpC.nb,formula. = ~ . - DistWoodS, cores=2)  
brmtgpC.nb.m62<-update(brmtgpC.nb,formula. = ~ . - DistWtrS, cores=2)  
brmtgpC.nb.m61<-update(brmtgpC.nb,formula. = ~ . - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m60<-update(brmtgpC.nb,formula. = ~ . - DistRdS, cores=2)  
brmtgpC.nb.m59<-update(brmtgpC.nb,formula. = ~ . - DistRdS - DistWoodS, cores=2)  
brmtgpC.nb.m58<-update(brmtgpC.nb,formula. = ~ . - DistRdS - DistWtrS, cores=2)  
brmtgpC.nb.m57<-update(brmtgpC.nb,formula. = ~ . - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m56<-update(brmtgpC.nb,formula. = ~ . - vegtype, cores=2)  
brmtgpC.nb.m55<-update(brmtgpC.nb,formula. = ~ . - vegtype - DistWoodS, cores=2)  
brmtgpC.nb.m54<-update(brmtgpC.nb,formula. = ~ . - vegtype - DistWtrS, cores=2)  
brmtgpC.nb.m53<-update(brmtgpC.nb,formula. = ~ . - vegtype - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m52<-update(brmtgpC.nb,formula. = ~ . - vegtype - DistRdS, cores=2)  
brmtgpC.nb.m51<-update(brmtgpC.nb,formula. = ~ . - vegtype - DistRdS - DistWoodS, cores=2)  
brmtgpC.nb.m50<-update(brmtgpC.nb,formula. = ~ . - vegtype - DistRdS - DistWtrS, cores=2)  
brmtgpC.nb.m49<-update(brmtgpC.nb,formula. = ~ . - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m48<-update(brmtgpC.nb,formula. = ~ . - TWIS , cores=2)  
brmtgpC.nb.m47<-update(brmtgpC.nb,formula. = ~ . - TWIS  - DistWoodS, cores=2)  
brmtgpC.nb.m46<-update(brmtgpC.nb,formula. = ~ . - TWIS  - DistWtrS, cores=2)  
brmtgpC.nb.m45<-update(brmtgpC.nb,formula. = ~ . - TWIS  - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m44<-update(brmtgpC.nb,formula. = ~ . - TWIS  - DistRdS, cores=2)  
brmtgpC.nb.m43<-update(brmtgpC.nb,formula. = ~ . - TWIS  - DistRdS - DistWoodS, cores=2)  
brmtgpC.nb.m42<-update(brmtgpC.nb,formula. = ~ . - TWIS  - DistRdS - DistWtrS, cores=2)  
brmtgpC.nb.m41<-update(brmtgpC.nb,formula. = ~ . - TWIS - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m40<-update(brmtgpC.nb,formula. = ~ . - TWIS - vegtype, cores=2)  
brmtgpC.nb.m39<-update(brmtgpC.nb,formula. = ~ . - TWIS - vegtype - DistWoodS, cores=2)  
brmtgpC.nb.m38<-update(brmtgpC.nb,formula. = ~ . - TWIS - vegtype - DistWtrS, cores=2)  
brmtgpC.nb.m37<-update(brmtgpC.nb,formula. = ~ . - TWIS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m36<-update(brmtgpC.nb,formula. = ~ . - TWIS - vegtype - DistRdS, cores=2)  
brmtgpC.nb.m35<-update(brmtgpC.nb,formula. = ~ . - TWIS - vegtype - DistRdS - DistWoodS, cores=2)  
brmtgpC.nb.m34<-update(brmtgpC.nb,formula. = ~ . - TWIS - vegtype - DistRdS - DistWtrS, cores=2)  
brmtgpC.nb.m33<-update(brmtgpC.nb,formula. = ~ . - TWIS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m32<-update(brmtgpC.nb,formula. = ~ . - eleS   , cores=2)  
brmtgpC.nb.m31<-update(brmtgpC.nb,formula. = ~ . - eleS -  DistWoodS, cores=2)  
brmtgpC.nb.m30<-update(brmtgpC.nb,formula. = ~ . - eleS -   DistWtrS, cores=2)  
brmtgpC.nb.m29<-update(brmtgpC.nb,formula. = ~ . - eleS -  DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m28<-update(brmtgpC.nb,formula. = ~ . - eleS -  DistRdS, cores=2)  
brmtgpC.nb.m27<-update(brmtgpC.nb,formula. = ~ . - eleS -  DistRdS - DistWoodS, cores=2)  
brmtgpC.nb.m26<-update(brmtgpC.nb,formula. = ~ . - eleS -  DistRdS - DistWtrS, cores=2)  
brmtgpC.nb.m25<-update(brmtgpC.nb,formula. = ~ . - eleS -  DistRdS - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m24<-update(brmtgpC.nb,formula. = ~ . - eleS - vegtype, cores=2)  
brmtgpC.nb.m23<-update(brmtgpC.nb,formula. = ~ . - eleS - vegtype - DistWoodS, cores=2)  
brmtgpC.nb.m22<-update(brmtgpC.nb,formula. = ~ . - eleS - vegtype - DistWtrS, cores=2)  
brmtgpC.nb.m21<-update(brmtgpC.nb,formula. = ~ . - eleS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m20<-update(brmtgpC.nb,formula. = ~ . - eleS - vegtype - DistRdS, cores=2)  
brmtgpC.nb.m19<-update(brmtgpC.nb,formula. = ~ . - eleS - vegtype - DistRdS - DistWoodS, cores=2)  
brmtgpC.nb.m18<-update(brmtgpC.nb,formula. = ~ . - eleS - vegtype - DistRdS - DistWtrS, cores=2)  
brmtgpC.nb.m17<-update(brmtgpC.nb,formula. = ~ . - eleS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m16<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS , cores=2)    
brmtgpC.nb.m15<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS  - DistWoodS, cores=2)  
brmtgpC.nb.m14<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS  - DistWtrS, cores=2)  
brmtgpC.nb.m13<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS  - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m12<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS  - DistRdS, cores=2)  
brmtgpC.nb.m11<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS  - DistRdS - DistWoodS, cores=2)  
brmtgpC.nb.m10<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS  - DistRdS - DistWtrS, cores=2)  
brmtgpC.nb.m9<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m8<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS - vegtype, cores=2)  
brmtgpC.nb.m7<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS - vegtype - DistWoodS, cores=2)  
brmtgpC.nb.m6<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS - vegtype - DistWtrS, cores=2)  
brmtgpC.nb.m5<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmtgpC.nb.m4<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS - vegtype - DistRdS, cores=2)  
brmtgpC.nb.m3<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWoodS, cores=2)  
brmtgpC.nb.m2<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWtrS, cores=2)  
brmtgpC.nb.m1<-update(brmtgpC.nb,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)

# Compute the WAIC for each of the 63 models.
waic.nb.tgp<-sapply(list(brmtgpC.nb.m1,brmtgpC.nb.m2,brmtgpC.nb.m3,brmtgpC.nb.m4,brmtgpC.nb.m5,brmtgpC.nb.m6,brmtgpC.nb.m7,brmtgpC.nb.m8,brmtgpC.nb.m9,brmtgpC.nb.m10,brmtgpC.nb.m11,brmtgpC.nb.m12,brmtgpC.nb.m13,brmtgpC.nb.m14,brmtgpC.nb.m15,brmtgpC.nb.m16,brmtgpC.nb.m17,brmtgpC.nb.m18,brmtgpC.nb.m19,brmtgpC.nb.m20,brmtgpC.nb.m21,brmtgpC.nb.m22,brmtgpC.nb.m23,brmtgpC.nb.m24,brmtgpC.nb.m25,brmtgpC.nb.m26,brmtgpC.nb.m27,brmtgpC.nb.m28,brmtgpC.nb.m29,brmtgpC.nb.m30,brmtgpC.nb.m31,brmtgpC.nb.m32,brmtgpC.nb.m33,brmtgpC.nb.m34,brmtgpC.nb.m35,brmtgpC.nb.m36,brmtgpC.nb.m37,brmtgpC.nb.m38,brmtgpC.nb.m39,brmtgpC.nb.m40,brmtgpC.nb.m41,brmtgpC.nb.m42,brmtgpC.nb.m43,brmtgpC.nb.m44,brmtgpC.nb.m45,brmtgpC.nb.m46,brmtgpC.nb.m47,brmtgpC.nb.m48,brmtgpC.nb.m49,brmtgpC.nb.m50,brmtgpC.nb.m51,brmtgpC.nb.m52,brmtgpC.nb.m53,brmtgpC.nb.m54,brmtgpC.nb.m55,brmtgpC.nb.m56,brmtgpC.nb.m57,brmtgpC.nb.m58,brmtgpC.nb.m59,brmtgpC.nb.m60, brmtgpC.nb.m61,brmtgpC.nb.m62,brmtgpC.nb.m63,brmtgpC.nb),function(x)waic(x)$estimates[3,1])

# Identify the model with the lowest WAIC.
waic.nb.tgp-min(waic.nb.tgp)

# Select the model with the lowest WAIC and most simple formula.
summary(brmtgpC.nb.m41)

emmeans(brmtgpC.nb.m41,pairwise~BaseT,type="response")

pairwise(brmtgpC.nb.m41)


#check for multicolinearity using Variance Inflation Factor 
tgp.nb.lm <- lm(TtgpC ~  eleS + TWIS + vegtype + DistRdS + DistWtrS + DistWoodS, data=GIS)
vif(tgp.nb.lm) 

## herbivore activity Model Residual Plots 
# For information on the diagnostics presented below, refer to vignette("DHARMa", package="DHARMa").

### Total 
Mtgp <- glmmTMB(TtgpC ~ eleS + vegtype, family =nbinom2(), data=GIS) # equivalent to brmtgpC.nb.m41
summary(Mtgp) 
plot(residuals(Mtgp))
testResiduals(Mtgp) 
tgpsim <- simulateResiduals(Mtgp, n=250, refit = F, plot = F, seed = 123) 
plot(tgpsim)  

## Analyses of residuals against independent variables

### Total 
plot(residuals(Mtgp) ~ GIS$eleS)
plot(residuals(Mtgp) ~ GIS$DistWtrS)
plot(residuals(Mtgp) ~ GIS$DistWoodS)
plot(residuals(Mtgp) ~ GIS$DistRdS)
plot(residuals(Mtgp) ~ GIS$TWIS)
boxplot(residuals(Mtgp) ~ GIS$vegtype, varwidth = TRUE)
subset(cbind(GIS, residual = (residuals(Mtgp)), fitted = fitted(Mtgp)), abs(residual) >45)
# no obvious error # model okay to proceed with



######### Analyses of presence or absence of herbivore groups #########

### Horse
brmhorPA.b <- brm(horPA ~  eleS + TWIS + vegtype +
                    DistRdS + DistWtrS + DistWoodS, family =bernoulli(), 
                  data=GIS, cores=2, prior=bprior,control=list(adapt_delta = 0.9))

brmhorPA.lm <- lm(horPA ~  eleS + TWIS + vegtype +
                    DistRdS + DistWtrS + DistWoodS, data=GIS)
vif(brmhorPA.lm) 

brmhorPA.b.m63<-update(brmhorPA.b,formula. = ~ . - DistWoodS, cores=2)  
brmhorPA.b.m62<-update(brmhorPA.b,formula. = ~ . - DistWtrS, cores=2)  
brmhorPA.b.m61<-update(brmhorPA.b,formula. = ~ . - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m60<-update(brmhorPA.b,formula. = ~ . - DistRdS, cores=2)  
brmhorPA.b.m59<-update(brmhorPA.b,formula. = ~ . - DistRdS - DistWoodS, cores=2)  
brmhorPA.b.m58<-update(brmhorPA.b,formula. = ~ . - DistRdS - DistWtrS, cores=2)  
brmhorPA.b.m57<-update(brmhorPA.b,formula. = ~ . - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m56<-update(brmhorPA.b,formula. = ~ . - vegtype, cores=2)  
brmhorPA.b.m55<-update(brmhorPA.b,formula. = ~ . - vegtype - DistWoodS, cores=2)  
brmhorPA.b.m54<-update(brmhorPA.b,formula. = ~ . - vegtype - DistWtrS, cores=2)  
brmhorPA.b.m53<-update(brmhorPA.b,formula. = ~ . - vegtype - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m52<-update(brmhorPA.b,formula. = ~ . - vegtype - DistRdS, cores=2)  
brmhorPA.b.m51<-update(brmhorPA.b,formula. = ~ . - vegtype - DistRdS - DistWoodS, cores=2)  
brmhorPA.b.m50<-update(brmhorPA.b,formula. = ~ . - vegtype - DistRdS - DistWtrS, cores=2)  
brmhorPA.b.m49<-update(brmhorPA.b,formula. = ~ . - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m48<-update(brmhorPA.b,formula. = ~ . - TWIS , cores=2)  
brmhorPA.b.m47<-update(brmhorPA.b,formula. = ~ . - TWIS  - DistWoodS, cores=2)  
brmhorPA.b.m46<-update(brmhorPA.b,formula. = ~ . - TWIS  - DistWtrS, cores=2)  
brmhorPA.b.m45<-update(brmhorPA.b,formula. = ~ . - TWIS  - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m44<-update(brmhorPA.b,formula. = ~ . - TWIS  - DistRdS, cores=2)  
brmhorPA.b.m43<-update(brmhorPA.b,formula. = ~ . - TWIS  - DistRdS - DistWoodS, cores=2)  
brmhorPA.b.m42<-update(brmhorPA.b,formula. = ~ . - TWIS  - DistRdS - DistWtrS, cores=2)  
brmhorPA.b.m41<-update(brmhorPA.b,formula. = ~ . - TWIS - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m40<-update(brmhorPA.b,formula. = ~ . - TWIS - vegtype, cores=2)  
brmhorPA.b.m39<-update(brmhorPA.b,formula. = ~ . - TWIS - vegtype - DistWoodS, cores=2)  
brmhorPA.b.m38<-update(brmhorPA.b,formula. = ~ . - TWIS - vegtype - DistWtrS, cores=2)  
brmhorPA.b.m37<-update(brmhorPA.b,formula. = ~ . - TWIS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m36<-update(brmhorPA.b,formula. = ~ . - TWIS - vegtype - DistRdS, cores=2)  
brmhorPA.b.m35<-update(brmhorPA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWoodS, cores=2)  
brmhorPA.b.m34<-update(brmhorPA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWtrS, cores=2)  
brmhorPA.b.m33<-update(brmhorPA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m32<-update(brmhorPA.b,formula. = ~ . - eleS   , cores=2)  
brmhorPA.b.m31<-update(brmhorPA.b,formula. = ~ . - eleS -  DistWoodS, cores=2)  
brmhorPA.b.m30<-update(brmhorPA.b,formula. = ~ . - eleS -   DistWtrS, cores=2)  
brmhorPA.b.m29<-update(brmhorPA.b,formula. = ~ . - eleS -  DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m28<-update(brmhorPA.b,formula. = ~ . - eleS -  DistRdS, cores=2)  
brmhorPA.b.m27<-update(brmhorPA.b,formula. = ~ . - eleS -  DistRdS - DistWoodS, cores=2)  
brmhorPA.b.m26<-update(brmhorPA.b,formula. = ~ . - eleS -  DistRdS - DistWtrS, cores=2)  
brmhorPA.b.m25<-update(brmhorPA.b,formula. = ~ . - eleS -  DistRdS - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m24<-update(brmhorPA.b,formula. = ~ . - eleS - vegtype, cores=2)  
brmhorPA.b.m23<-update(brmhorPA.b,formula. = ~ . - eleS - vegtype - DistWoodS, cores=2)  
brmhorPA.b.m22<-update(brmhorPA.b,formula. = ~ . - eleS - vegtype - DistWtrS, cores=2)  
brmhorPA.b.m21<-update(brmhorPA.b,formula. = ~ . - eleS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m20<-update(brmhorPA.b,formula. = ~ . - eleS - vegtype - DistRdS, cores=2)  
brmhorPA.b.m19<-update(brmhorPA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWoodS, cores=2)  
brmhorPA.b.m18<-update(brmhorPA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWtrS, cores=2)  
brmhorPA.b.m17<-update(brmhorPA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m16<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS , cores=2)    
brmhorPA.b.m15<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS  - DistWoodS, cores=2)  
brmhorPA.b.m14<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS  - DistWtrS, cores=2)  
brmhorPA.b.m13<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS  - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m12<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS  - DistRdS, cores=2)  
brmhorPA.b.m11<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS  - DistRdS - DistWoodS, cores=2)  
brmhorPA.b.m10<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS  - DistRdS - DistWtrS, cores=2)  
brmhorPA.b.m9<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m8<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS - vegtype, cores=2)  
brmhorPA.b.m7<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWoodS, cores=2)  
brmhorPA.b.m6<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWtrS, cores=2)  
brmhorPA.b.m5<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmhorPA.b.m4<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS, cores=2)  
brmhorPA.b.m3<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWoodS, cores=2)  
brmhorPA.b.m2<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWtrS, cores=2)  
brmhorPA.b.m1<-update(brmhorPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)

#Compute the WAIC for each of the 63 models.
waic.b.hor<-sapply(list(brmhorPA.b.m1,brmhorPA.b.m2,brmhorPA.b.m3,brmhorPA.b.m4,brmhorPA.b.m5,brmhorPA.b.m6,brmhorPA.b.m7,brmhorPA.b.m8,brmhorPA.b.m9,brmhorPA.b.m10,brmhorPA.b.m11,brmhorPA.b.m12,brmhorPA.b.m13,brmhorPA.b.m14,brmhorPA.b.m15,brmhorPA.b.m16,brmhorPA.b.m17,brmhorPA.b.m18,brmhorPA.b.m19,brmhorPA.b.m20,brmhorPA.b.m21,brmhorPA.b.m22,brmhorPA.b.m23,brmhorPA.b.m24,brmhorPA.b.m25,brmhorPA.b.m26,brmhorPA.b.m27,brmhorPA.b.m28,brmhorPA.b.m29,brmhorPA.b.m30,brmhorPA.b.m31,brmhorPA.b.m32,brmhorPA.b.m33,brmhorPA.b.m34,brmhorPA.b.m35,brmhorPA.b.m36,brmhorPA.b.m37,brmhorPA.b.m38,brmhorPA.b.m39,brmhorPA.b.m40,brmhorPA.b.m41,brmhorPA.b.m42,brmhorPA.b.m43,brmhorPA.b.m44,brmhorPA.b.m45,brmhorPA.b.m46,brmhorPA.b.m47,brmhorPA.b.m48,brmhorPA.b.m49,brmhorPA.b.m50,brmhorPA.b.m51,brmhorPA.b.m52,brmhorPA.b.m53,brmhorPA.b.m54,brmhorPA.b.m55,brmhorPA.b.m56,brmhorPA.b.m57,brmhorPA.b.m58,brmhorPA.b.m59,brmhorPA.b.m60,brmhorPA.b.m61,brmhorPA.b.m62,brmhorPA.b.m63,brmhorPA.b),function(x)waic(x)$estimates[3,1])

#Identify the model with the lowest WAIC.
waic.b.hor-min(waic.b.hor)

# Select the model with the lowest WAIC and most simple formula.  
summary(brmhorPA.b.m35) 


### Deer
brmdeePA.b <- brm(deePA ~  eleS + TWIS + vegtype +
                    DistRdS + DistWtrS + DistWoodS, family =bernoulli(), 
                  data=GIS, cores=2, prior=bprior,control=list(adapt_delta = 0.9))

brmdeePA.lm <- lm(deePA ~  eleS + TWIS + vegtype +
                    DistRdS + DistWtrS + DistWoodS, data=GIS)
vif(brmdeePA.lm) 


brmdeePA.b.m63<-update(brmdeePA.b,formula. = ~ . - DistWoodS, cores=2)  
brmdeePA.b.m62<-update(brmdeePA.b,formula. = ~ . - DistWtrS, cores=2)  
brmdeePA.b.m61<-update(brmdeePA.b,formula. = ~ . - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m60<-update(brmdeePA.b,formula. = ~ . - DistRdS, cores=2)  
brmdeePA.b.m59<-update(brmdeePA.b,formula. = ~ . - DistRdS - DistWoodS, cores=2)  
brmdeePA.b.m58<-update(brmdeePA.b,formula. = ~ . - DistRdS - DistWtrS, cores=2)  
brmdeePA.b.m57<-update(brmdeePA.b,formula. = ~ . - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m56<-update(brmdeePA.b,formula. = ~ . - vegtype, cores=2)  
brmdeePA.b.m55<-update(brmdeePA.b,formula. = ~ . - vegtype - DistWoodS, cores=2)  
brmdeePA.b.m54<-update(brmdeePA.b,formula. = ~ . - vegtype - DistWtrS, cores=2)  
brmdeePA.b.m53<-update(brmdeePA.b,formula. = ~ . - vegtype - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m52<-update(brmdeePA.b,formula. = ~ . - vegtype - DistRdS, cores=2)  
brmdeePA.b.m51<-update(brmdeePA.b,formula. = ~ . - vegtype - DistRdS - DistWoodS, cores=2)  
brmdeePA.b.m50<-update(brmdeePA.b,formula. = ~ . - vegtype - DistRdS - DistWtrS, cores=2)  
brmdeePA.b.m49<-update(brmdeePA.b,formula. = ~ . - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m48<-update(brmdeePA.b,formula. = ~ . - TWIS , cores=2)  
brmdeePA.b.m47<-update(brmdeePA.b,formula. = ~ . - TWIS  - DistWoodS, cores=2)  
brmdeePA.b.m46<-update(brmdeePA.b,formula. = ~ . - TWIS  - DistWtrS, cores=2)  
brmdeePA.b.m45<-update(brmdeePA.b,formula. = ~ . - TWIS  - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m44<-update(brmdeePA.b,formula. = ~ . - TWIS  - DistRdS, cores=2)  
brmdeePA.b.m43<-update(brmdeePA.b,formula. = ~ . - TWIS  - DistRdS - DistWoodS, cores=2)  
brmdeePA.b.m42<-update(brmdeePA.b,formula. = ~ . - TWIS  - DistRdS - DistWtrS, cores=2)  
brmdeePA.b.m41<-update(brmdeePA.b,formula. = ~ . - TWIS - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m40<-update(brmdeePA.b,formula. = ~ . - TWIS - vegtype, cores=2)  
brmdeePA.b.m39<-update(brmdeePA.b,formula. = ~ . - TWIS - vegtype - DistWoodS, cores=2)  
brmdeePA.b.m38<-update(brmdeePA.b,formula. = ~ . - TWIS - vegtype - DistWtrS, cores=2)  
brmdeePA.b.m37<-update(brmdeePA.b,formula. = ~ . - TWIS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m36<-update(brmdeePA.b,formula. = ~ . - TWIS - vegtype - DistRdS, cores=2)  
brmdeePA.b.m35<-update(brmdeePA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWoodS, cores=2)  
brmdeePA.b.m34<-update(brmdeePA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWtrS, cores=2)  
brmdeePA.b.m33<-update(brmdeePA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m32<-update(brmdeePA.b,formula. = ~ . - eleS   , cores=2)  
brmdeePA.b.m31<-update(brmdeePA.b,formula. = ~ . - eleS -  DistWoodS, cores=2)  
brmdeePA.b.m30<-update(brmdeePA.b,formula. = ~ . - eleS -   DistWtrS, cores=2)  
brmdeePA.b.m29<-update(brmdeePA.b,formula. = ~ . - eleS -  DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m28<-update(brmdeePA.b,formula. = ~ . - eleS -  DistRdS, cores=2)  
brmdeePA.b.m27<-update(brmdeePA.b,formula. = ~ . - eleS -  DistRdS - DistWoodS, cores=2)  
brmdeePA.b.m26<-update(brmdeePA.b,formula. = ~ . - eleS -  DistRdS - DistWtrS, cores=2)  
brmdeePA.b.m25<-update(brmdeePA.b,formula. = ~ . - eleS -  DistRdS - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m24<-update(brmdeePA.b,formula. = ~ . - eleS - vegtype, cores=2)  
brmdeePA.b.m23<-update(brmdeePA.b,formula. = ~ . - eleS - vegtype - DistWoodS, cores=2)  
brmdeePA.b.m22<-update(brmdeePA.b,formula. = ~ . - eleS - vegtype - DistWtrS, cores=2)  
brmdeePA.b.m21<-update(brmdeePA.b,formula. = ~ . - eleS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m20<-update(brmdeePA.b,formula. = ~ . - eleS - vegtype - DistRdS, cores=2)  
brmdeePA.b.m19<-update(brmdeePA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWoodS, cores=2)  
brmdeePA.b.m18<-update(brmdeePA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWtrS, cores=2)  
brmdeePA.b.m17<-update(brmdeePA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m16<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS , cores=2)    
brmdeePA.b.m15<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS  - DistWoodS, cores=2)  
brmdeePA.b.m14<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS  - DistWtrS, cores=2)  
brmdeePA.b.m13<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS  - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m12<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS  - DistRdS, cores=2)  
brmdeePA.b.m11<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS  - DistRdS - DistWoodS, cores=2)  
brmdeePA.b.m10<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS  - DistRdS - DistWtrS, cores=2)  
brmdeePA.b.m9<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m8<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS - vegtype, cores=2)  
brmdeePA.b.m7<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWoodS, cores=2)  
brmdeePA.b.m6<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWtrS, cores=2)  
brmdeePA.b.m5<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmdeePA.b.m4<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS, cores=2)  
brmdeePA.b.m3<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWoodS, cores=2)  
brmdeePA.b.m2<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWtrS, cores=2)  
brmdeePA.b.m1<-update(brmdeePA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)

#Compute the WAIC for each of the 63 models.
waic.b.dee<-sapply(list(brmdeePA.b.m1,brmdeePA.b.m2,brmdeePA.b.m3,brmdeePA.b.m4,brmdeePA.b.m5,brmdeePA.b.m6,brmdeePA.b.m7,brmdeePA.b.m8,brmdeePA.b.m9,brmdeePA.b.m10,brmdeePA.b.m11,brmdeePA.b.m12,brmdeePA.b.m13,brmdeePA.b.m14,brmdeePA.b.m15,brmdeePA.b.m16,brmdeePA.b.m17,brmdeePA.b.m18,brmdeePA.b.m19,brmdeePA.b.m20,brmdeePA.b.m21,brmdeePA.b.m22,brmdeePA.b.m23,brmdeePA.b.m24,brmdeePA.b.m25,brmdeePA.b.m26,brmdeePA.b.m27,brmdeePA.b.m28,brmdeePA.b.m29,brmdeePA.b.m30,brmdeePA.b.m31,brmdeePA.b.m32,brmdeePA.b.m33,brmdeePA.b.m34,brmdeePA.b.m35,brmdeePA.b.m36,brmdeePA.b.m37,brmdeePA.b.m38,brmdeePA.b.m39,brmdeePA.b.m40,brmdeePA.b.m41,brmdeePA.b.m42,brmdeePA.b.m43,brmdeePA.b.m44,brmdeePA.b.m45,brmdeePA.b.m46,brmdeePA.b.m47,brmdeePA.b.m48,brmdeePA.b.m49,brmdeePA.b.m50,brmdeePA.b.m51,brmdeePA.b.m52,brmdeePA.b.m53,brmdeePA.b.m54,brmdeePA.b.m55,brmdeePA.b.m56,brmdeePA.b.m57,brmdeePA.b.m58,brmdeePA.b.m59,brmdeePA.b.m60,brmdeePA.b.m61,brmdeePA.b.m62,brmdeePA.b.m63,brmdeePA.b),function(x)waic(x)$estimates[3,1])

#Identify the model with the lowest WAIC.
waic.b.dee-min(waic.b.dee)

# Select the model with the lowest WAIC and most simple formula.  
summary(brmdeePA.b.m43) 


### Leporid
brmlepPA.b <- brm(lepPA ~  eleS + TWIS + vegtype +
                    DistRdS + DistWtrS + DistWoodS, family =bernoulli(), 
                  data=GIS, cores=2, prior=bprior,control=list(adapt_delta = 0.9))

brmlepPA.lm <- lm(lepPA ~  eleS + TWIS + vegtype +
                    DistRdS + DistWtrS + DistWoodS, data=GIS)
vif(brmlepPA.lm) 

brmlepPA.b.m63<-update(brmlepPA.b,formula. = ~ . - DistWoodS, cores=2)  
brmlepPA.b.m62<-update(brmlepPA.b,formula. = ~ . - DistWtrS, cores=2)  
brmlepPA.b.m61<-update(brmlepPA.b,formula. = ~ . - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m60<-update(brmlepPA.b,formula. = ~ . - DistRdS, cores=2)  
brmlepPA.b.m59<-update(brmlepPA.b,formula. = ~ . - DistRdS - DistWoodS, cores=2)  
brmlepPA.b.m58<-update(brmlepPA.b,formula. = ~ . - DistRdS - DistWtrS, cores=2)  
brmlepPA.b.m57<-update(brmlepPA.b,formula. = ~ . - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m56<-update(brmlepPA.b,formula. = ~ . - vegtype, cores=2)  
brmlepPA.b.m55<-update(brmlepPA.b,formula. = ~ . - vegtype - DistWoodS, cores=2)  
brmlepPA.b.m54<-update(brmlepPA.b,formula. = ~ . - vegtype - DistWtrS, cores=2)  
brmlepPA.b.m53<-update(brmlepPA.b,formula. = ~ . - vegtype - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m52<-update(brmlepPA.b,formula. = ~ . - vegtype - DistRdS, cores=2)  
brmlepPA.b.m51<-update(brmlepPA.b,formula. = ~ . - vegtype - DistRdS - DistWoodS, cores=2)  
brmlepPA.b.m50<-update(brmlepPA.b,formula. = ~ . - vegtype - DistRdS - DistWtrS, cores=2)  
brmlepPA.b.m49<-update(brmlepPA.b,formula. = ~ . - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m48<-update(brmlepPA.b,formula. = ~ . - TWIS , cores=2)  
brmlepPA.b.m47<-update(brmlepPA.b,formula. = ~ . - TWIS  - DistWoodS, cores=2)  
brmlepPA.b.m46<-update(brmlepPA.b,formula. = ~ . - TWIS  - DistWtrS, cores=2)  
brmlepPA.b.m45<-update(brmlepPA.b,formula. = ~ . - TWIS  - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m44<-update(brmlepPA.b,formula. = ~ . - TWIS  - DistRdS, cores=2)  
brmlepPA.b.m43<-update(brmlepPA.b,formula. = ~ . - TWIS  - DistRdS - DistWoodS, cores=2)  
brmlepPA.b.m42<-update(brmlepPA.b,formula. = ~ . - TWIS  - DistRdS - DistWtrS, cores=2)  
brmlepPA.b.m41<-update(brmlepPA.b,formula. = ~ . - TWIS - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m40<-update(brmlepPA.b,formula. = ~ . - TWIS - vegtype, cores=2)  
brmlepPA.b.m39<-update(brmlepPA.b,formula. = ~ . - TWIS - vegtype - DistWoodS, cores=2)  
brmlepPA.b.m38<-update(brmlepPA.b,formula. = ~ . - TWIS - vegtype - DistWtrS, cores=2)  
brmlepPA.b.m37<-update(brmlepPA.b,formula. = ~ . - TWIS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m36<-update(brmlepPA.b,formula. = ~ . - TWIS - vegtype - DistRdS, cores=2)  
brmlepPA.b.m35<-update(brmlepPA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWoodS, cores=2)  
brmlepPA.b.m34<-update(brmlepPA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWtrS, cores=2)  
brmlepPA.b.m33<-update(brmlepPA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m32<-update(brmlepPA.b,formula. = ~ . - eleS   , cores=2)  
brmlepPA.b.m31<-update(brmlepPA.b,formula. = ~ . - eleS -  DistWoodS, cores=2)  
brmlepPA.b.m30<-update(brmlepPA.b,formula. = ~ . - eleS -   DistWtrS, cores=2)  
brmlepPA.b.m29<-update(brmlepPA.b,formula. = ~ . - eleS -  DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m28<-update(brmlepPA.b,formula. = ~ . - eleS -  DistRdS, cores=2)  
brmlepPA.b.m27<-update(brmlepPA.b,formula. = ~ . - eleS -  DistRdS - DistWoodS, cores=2)  
brmlepPA.b.m26<-update(brmlepPA.b,formula. = ~ . - eleS -  DistRdS - DistWtrS, cores=2)  
brmlepPA.b.m25<-update(brmlepPA.b,formula. = ~ . - eleS -  DistRdS - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m24<-update(brmlepPA.b,formula. = ~ . - eleS - vegtype, cores=2)  
brmlepPA.b.m23<-update(brmlepPA.b,formula. = ~ . - eleS - vegtype - DistWoodS, cores=2)  
brmlepPA.b.m22<-update(brmlepPA.b,formula. = ~ . - eleS - vegtype - DistWtrS, cores=2)  
brmlepPA.b.m21<-update(brmlepPA.b,formula. = ~ . - eleS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m20<-update(brmlepPA.b,formula. = ~ . - eleS - vegtype - DistRdS, cores=2)  
brmlepPA.b.m19<-update(brmlepPA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWoodS, cores=2)  
brmlepPA.b.m18<-update(brmlepPA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWtrS, cores=2)  
brmlepPA.b.m17<-update(brmlepPA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m16<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS , cores=2)    
brmlepPA.b.m15<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS  - DistWoodS, cores=2)  
brmlepPA.b.m14<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS  - DistWtrS, cores=2)  
brmlepPA.b.m13<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS  - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m12<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS  - DistRdS, cores=2)  
brmlepPA.b.m11<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS  - DistRdS - DistWoodS, cores=2)  
brmlepPA.b.m10<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS  - DistRdS - DistWtrS, cores=2)  
brmlepPA.b.m9<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m8<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS - vegtype, cores=2)  
brmlepPA.b.m7<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWoodS, cores=2)  
brmlepPA.b.m6<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWtrS, cores=2)  
brmlepPA.b.m5<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmlepPA.b.m4<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS, cores=2)  
brmlepPA.b.m3<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWoodS, cores=2)  
brmlepPA.b.m2<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWtrS, cores=2)  
brmlepPA.b.m1<-update(brmlepPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)

#Compute the WAIC for each of the 63 models.
waic.b.lep<-sapply(list(brmlepPA.b.m1,brmlepPA.b.m2,brmlepPA.b.m3,brmlepPA.b.m4,brmlepPA.b.m5,brmlepPA.b.m6,brmlepPA.b.m7,brmlepPA.b.m8,brmlepPA.b.m9,brmlepPA.b.m10,brmlepPA.b.m11,brmlepPA.b.m12,brmlepPA.b.m13,brmlepPA.b.m14,brmlepPA.b.m15,brmlepPA.b.m16,brmlepPA.b.m17,brmlepPA.b.m18,brmlepPA.b.m19,brmlepPA.b.m20,brmlepPA.b.m21,brmlepPA.b.m22,brmlepPA.b.m23,brmlepPA.b.m24,brmlepPA.b.m25,brmlepPA.b.m26,brmlepPA.b.m27,brmlepPA.b.m28,brmlepPA.b.m29,brmlepPA.b.m30,brmlepPA.b.m31,brmlepPA.b.m32,brmlepPA.b.m33,brmlepPA.b.m34,brmlepPA.b.m35,brmlepPA.b.m36,brmlepPA.b.m37,brmlepPA.b.m38,brmlepPA.b.m39,brmlepPA.b.m40,brmlepPA.b.m41,brmlepPA.b.m42,brmlepPA.b.m43,brmlepPA.b.m44,brmlepPA.b.m45,brmlepPA.b.m46,brmlepPA.b.m47,brmlepPA.b.m48,brmlepPA.b.m49,brmlepPA.b.m50,brmlepPA.b.m51,brmlepPA.b.m52,brmlepPA.b.m53,brmlepPA.b.m54,brmlepPA.b.m55,brmlepPA.b.m56,brmlepPA.b.m57,brmlepPA.b.m58,brmlepPA.b.m59,brmlepPA.b.m60,brmlepPA.b.m61,brmlepPA.b.m62,brmlepPA.b.m63,brmlepPA.b),function(x)waic(x)$estimates[3,1])

#Identify the model with the lowest WAIC.
waic.b.lep-min(waic.b.lep)

# Select the model with the lowest WAIC and most simple formula. 
summary(brmlepPA.b.m3)


### Macropod
brmmacPA.b <- brm(macPA ~  eleS + TWIS + vegtype +
                    DistRdS + DistWtrS + DistWoodS, family =bernoulli(), 
                  data=GIS, cores=2, prior=bprior,control=list(adapt_delta = 0.9))

brmmacPA.lm <- lm(macPA ~  eleS + TWIS + vegtype +
                    DistRdS + DistWtrS + DistWoodS, data=GIS)
vif(brmmacPA.lm) 


brmmacPA.b.m63<-update(brmmacPA.b,formula. = ~ . - DistWoodS, cores=2)  
brmmacPA.b.m62<-update(brmmacPA.b,formula. = ~ . - DistWtrS, cores=2)  
brmmacPA.b.m61<-update(brmmacPA.b,formula. = ~ . - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m60<-update(brmmacPA.b,formula. = ~ . - DistRdS, cores=2)  
brmmacPA.b.m59<-update(brmmacPA.b,formula. = ~ . - DistRdS - DistWoodS, cores=2)  
brmmacPA.b.m58<-update(brmmacPA.b,formula. = ~ . - DistRdS - DistWtrS, cores=2)  
brmmacPA.b.m57<-update(brmmacPA.b,formula. = ~ . - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m56<-update(brmmacPA.b,formula. = ~ . - vegtype, cores=2)  
brmmacPA.b.m55<-update(brmmacPA.b,formula. = ~ . - vegtype - DistWoodS, cores=2)  
brmmacPA.b.m54<-update(brmmacPA.b,formula. = ~ . - vegtype - DistWtrS, cores=2)  
brmmacPA.b.m53<-update(brmmacPA.b,formula. = ~ . - vegtype - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m52<-update(brmmacPA.b,formula. = ~ . - vegtype - DistRdS, cores=2)  
brmmacPA.b.m51<-update(brmmacPA.b,formula. = ~ . - vegtype - DistRdS - DistWoodS, cores=2)  
brmmacPA.b.m50<-update(brmmacPA.b,formula. = ~ . - vegtype - DistRdS - DistWtrS, cores=2)  
brmmacPA.b.m49<-update(brmmacPA.b,formula. = ~ . - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m48<-update(brmmacPA.b,formula. = ~ . - TWIS , cores=2)  
brmmacPA.b.m47<-update(brmmacPA.b,formula. = ~ . - TWIS  - DistWoodS, cores=2)  
brmmacPA.b.m46<-update(brmmacPA.b,formula. = ~ . - TWIS  - DistWtrS, cores=2)  
brmmacPA.b.m45<-update(brmmacPA.b,formula. = ~ . - TWIS  - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m44<-update(brmmacPA.b,formula. = ~ . - TWIS  - DistRdS, cores=2)  
brmmacPA.b.m43<-update(brmmacPA.b,formula. = ~ . - TWIS  - DistRdS - DistWoodS, cores=2)  
brmmacPA.b.m42<-update(brmmacPA.b,formula. = ~ . - TWIS  - DistRdS - DistWtrS, cores=2)  
brmmacPA.b.m41<-update(brmmacPA.b,formula. = ~ . - TWIS - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m40<-update(brmmacPA.b,formula. = ~ . - TWIS - vegtype, cores=2)  
brmmacPA.b.m39<-update(brmmacPA.b,formula. = ~ . - TWIS - vegtype - DistWoodS, cores=2)  
brmmacPA.b.m38<-update(brmmacPA.b,formula. = ~ . - TWIS - vegtype - DistWtrS, cores=2)  
brmmacPA.b.m37<-update(brmmacPA.b,formula. = ~ . - TWIS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m36<-update(brmmacPA.b,formula. = ~ . - TWIS - vegtype - DistRdS, cores=2)  
brmmacPA.b.m35<-update(brmmacPA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWoodS, cores=2)  
brmmacPA.b.m34<-update(brmmacPA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWtrS, cores=2)  
brmmacPA.b.m33<-update(brmmacPA.b,formula. = ~ . - TWIS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m32<-update(brmmacPA.b,formula. = ~ . - eleS   , cores=2)  
brmmacPA.b.m31<-update(brmmacPA.b,formula. = ~ . - eleS -  DistWoodS, cores=2)  
brmmacPA.b.m30<-update(brmmacPA.b,formula. = ~ . - eleS -   DistWtrS, cores=2)  
brmmacPA.b.m29<-update(brmmacPA.b,formula. = ~ . - eleS -  DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m28<-update(brmmacPA.b,formula. = ~ . - eleS -  DistRdS, cores=2)  
brmmacPA.b.m27<-update(brmmacPA.b,formula. = ~ . - eleS -  DistRdS - DistWoodS, cores=2)  
brmmacPA.b.m26<-update(brmmacPA.b,formula. = ~ . - eleS -  DistRdS - DistWtrS, cores=2)  
brmmacPA.b.m25<-update(brmmacPA.b,formula. = ~ . - eleS -  DistRdS - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m24<-update(brmmacPA.b,formula. = ~ . - eleS - vegtype, cores=2)  
brmmacPA.b.m23<-update(brmmacPA.b,formula. = ~ . - eleS - vegtype - DistWoodS, cores=2)  
brmmacPA.b.m22<-update(brmmacPA.b,formula. = ~ . - eleS - vegtype - DistWtrS, cores=2)  
brmmacPA.b.m21<-update(brmmacPA.b,formula. = ~ . - eleS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m20<-update(brmmacPA.b,formula. = ~ . - eleS - vegtype - DistRdS, cores=2)  
brmmacPA.b.m19<-update(brmmacPA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWoodS, cores=2)  
brmmacPA.b.m18<-update(brmmacPA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWtrS, cores=2)  
brmmacPA.b.m17<-update(brmmacPA.b,formula. = ~ . - eleS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m16<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS , cores=2)    
brmmacPA.b.m15<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS  - DistWoodS, cores=2)  
brmmacPA.b.m14<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS  - DistWtrS, cores=2)  
brmmacPA.b.m13<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS  - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m12<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS  - DistRdS, cores=2)  
brmmacPA.b.m11<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS  - DistRdS - DistWoodS, cores=2)  
brmmacPA.b.m10<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS  - DistRdS - DistWtrS, cores=2)  
brmmacPA.b.m9<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS - DistRdS - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m8<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS - vegtype, cores=2)  
brmmacPA.b.m7<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWoodS, cores=2)  
brmmacPA.b.m6<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWtrS, cores=2)  
brmmacPA.b.m5<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistWtrS - DistWoodS, cores=2)  
brmmacPA.b.m4<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS, cores=2)  
brmmacPA.b.m3<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWoodS, cores=2)  
brmmacPA.b.m2<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWtrS, cores=2)  
brmmacPA.b.m1<-update(brmmacPA.b,formula. = ~ . - eleS - TWIS - vegtype - DistRdS - DistWtrS - DistWoodS, cores=2)

#Compute the WAIC for each of the 63 models.
waic.b.mac<-sapply(list(brmmacPA.b.m1,brmmacPA.b.m2,brmmacPA.b.m3,brmmacPA.b.m4,brmmacPA.b.m5,brmmacPA.b.m6,brmmacPA.b.m7,brmmacPA.b.m8,brmmacPA.b.m9,brmmacPA.b.m10,brmmacPA.b.m11,brmmacPA.b.m12,brmmacPA.b.m13,brmmacPA.b.m14,brmmacPA.b.m15,brmmacPA.b.m16,brmmacPA.b.m17,brmmacPA.b.m18,brmmacPA.b.m19,brmmacPA.b.m20,brmmacPA.b.m21,brmmacPA.b.m22,brmmacPA.b.m23,brmmacPA.b.m24,brmmacPA.b.m25,brmmacPA.b.m26,brmmacPA.b.m27,brmmacPA.b.m28,brmmacPA.b.m29,brmmacPA.b.m30,brmmacPA.b.m31,brmmacPA.b.m32,brmmacPA.b.m33,brmmacPA.b.m34,brmmacPA.b.m35,brmmacPA.b.m36,brmmacPA.b.m37,brmmacPA.b.m38,brmmacPA.b.m39,brmmacPA.b.m40,brmmacPA.b.m41,brmmacPA.b.m42,brmmacPA.b.m43,brmmacPA.b.m44,brmmacPA.b.m45,brmmacPA.b.m46,brmmacPA.b.m47,brmmacPA.b.m48,brmmacPA.b.m49,brmmacPA.b.m50,brmmacPA.b.m51,brmmacPA.b.m52,brmmacPA.b.m53,brmmacPA.b.m54,brmmacPA.b.m55,brmmacPA.b.m56,brmmacPA.b.m57,brmmacPA.b.m58,brmmacPA.b.m59,brmmacPA.b.m60,brmmacPA.b.m61,brmmacPA.b.m62,brmmacPA.b.m63,brmmacPA.b),function(x)waic(x)$estimates[3,1])

#Identify the model with the lowest WAIC.
waic.b.mac-min(waic.b.mac)

# Select the model with the lowest WAIC and most simple formula.  
summary(brmmacPA.b.m41)


## herbivore Presence Model Diagnostics
# For information on the diagnostics presented below, refer to vignette("DHARMa", package="DHARMa").

### Horse
summary(brmhorPA.b.m35) 
plot(fitted(brmhorPA.b.m35)[,1],residuals(brmhorPA.b.m35)[,1]) 
hor.res <- glmmTMB(horPA ~  eleS + DistWtrS , family = binomial(link = "logit"), data=GIS)
summary(hor.res) 
plot(residuals(hor.res))
testResiduals(hor.res)
horsim <- simulateResiduals(hor.res, n=250, refit = F, plot = F, seed = 123) 
plot(horsim) 


### Deer
summary(brmdeePA.b.m43)
plot(fitted(brmdeePA.b.m43)[,1],residuals(brmdeePA.b.m43)[,1]) 
dee.res <- glmmTMB(deePA ~  eleS + vegtype + DistWtrS, family = binomial(link = "logit"), data=GIS)
summary(dee.res) 
plot(residuals(dee.res))
testResiduals(dee.res)
deesim <- simulateResiduals(dee.res, n=250, refit = F, plot = F, seed = 123) 
plot(deesim) 


### Leporid
summary(brmlepPA.b.m3)
plot(fitted(brmlepPA.b.m3)[,1],residuals(brmlepPA.b.m3)[,1]) 
lep.res <- glmmTMB(lepPA ~  DistWtrS,
                   family = binomial(link = "logit"), data=GIS)
summary(lep.res) 
plot(residuals(lep.res))
testResiduals(lep.res)
lepsim <- simulateResiduals(lep.res, n=250, refit = F, plot = F, seed = 123) 
plot(lepsim) 


### Macropod
summary(brmmacPA.b.m41) 
plot(fitted(brmmacPA.b.m41)[,1],residuals(brmmacPA.b.m41)[,1]) 
roo.res <- glmmTMB(macPA ~  eleS + vegtype,
                   family = binomial(link = "logit"), data=GIS)
summary(roo.res) 
plot(residuals(roo.res))
testResiduals(roo.res)
macSim <- simulateResiduals(roo.res, n=250, refit = F, plot = F, seed = 123) 
plot(macSim) 


###########  VEG AND SOIL ASSOCIATIONS  ######################


####### Height #######

## vs. Total herbivore activity analysis  
brm.grassheight.tgp <- brm(grassheight.mean ~  tgpS, 
                      family=Gamma(link = "log"), data=grassheight.noNA, cores=2, 
                      prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.grassheight.tgp)

## vs each herbivore family analysis  
brm.grassheight <- brm(grassheight.mean ~  horS + deeS  + lepS + macS , 
                  family=Gamma(link = "log"), data=grassheight.noNA, cores=2, 
                  prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.grassheight) 


grassheight.lm <- lm(grassheight.mean ~  horS + deeS  + lepS + macS, data=grassheight.noNA) # check multicolinearity
vif(grassheight.lm) 


brm.grassheight.m15<-update(brm.grassheight,formula. = ~ .  - macS, cores=2)  
brm.grassheight.m14<-update(brm.grassheight,formula. = ~ .  - lepS, cores=2)  
brm.grassheight.m13<-update(brm.grassheight,formula. = ~ .  - lepS - macS, cores=2)  
brm.grassheight.m12<-update(brm.grassheight,formula. = ~ . - deeS , cores=2)  
brm.grassheight.m11<-update(brm.grassheight,formula. = ~ . - deeS  - macS, cores=2)  
brm.grassheight.m10<-update(brm.grassheight,formula. = ~ . - deeS  - lepS, cores=2)  
brm.grassheight.m9<-update(brm.grassheight,formula. = ~ . - deeS  - lepS - macS, cores=2)  
brm.grassheight.m8<-update(brm.grassheight,formula. = ~ . - horS  , cores=2)  
brm.grassheight.m7<-update(brm.grassheight,formula. = ~ . - horS  - macS, cores=2)  
brm.grassheight.m6<-update(brm.grassheight,formula. = ~ . - horS  - lepS, cores=2)  
brm.grassheight.m5<-update(brm.grassheight,formula. = ~ . - horS  - lepS - macS, cores=2)  
brm.grassheight.m4<-update(brm.grassheight,formula. = ~ . - horS - deeS , cores=2)  
brm.grassheight.m3<-update(brm.grassheight,formula. = ~ . - horS - deeS  - macS, cores=2)  
brm.grassheight.m2<-update(brm.grassheight,formula. = ~ . - horS - deeS  - lepS, cores=2)  
brm.grassheight.m1<-update(brm.grassheight,formula. = ~ . - horS - deeS  - lepS - macS, cores=2)  

# Compute the WAIC for each of the 16 models.
waic.grassheight<-sapply(list(brm.grassheight.m1,brm.grassheight.m2,brm.grassheight.m3,brm.grassheight.m4,brm.grassheight.m5,brm.grassheight.m6,brm.grassheight.m7,brm.grassheight.m8,brm.grassheight.m9,brm.grassheight.m10,brm.grassheight.m11,brm.grassheight.m12,brm.grassheight.m13,brm.grassheight.m14,brm.grassheight.m15,brm.grassheight,brm.grassheight.tgp),function(x)waic(x)$estimates[3,1])

# Identify the model with the lowest WAIC.
waic.grassheight-min(waic.grassheight)

# Select the model with the lowest WAIC and most simple formula.
summary(brm.grassheight.m11)


####### Foliar density at groundstorey #######

## vs. Total herbivore activity analysis  
brm.foldens.tgp <- brm(fol5.mean ~  tgpS, 
                       family=Gamma(link = "log"), data=GIS, cores=2, 
                       prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.foldens.tgp) 


## vs each herbivore family analysis  
brm.foldens <- brm(fol5.mean ~  horS + deeS  + lepS + macS , 
                   family=Gamma(link = "log"), data=GIS, cores=2, 
                   prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.foldens) 

foldens.lm <- lm(fol5.mean ~  horS + deeS  + lepS + macS, data=GIS) # check multicolinearity
vif(foldens.lm) 


brm.foldens.m15<-update(brm.foldens,formula. = ~ .  - macS, cores=2)  
brm.foldens.m14<-update(brm.foldens,formula. = ~ .  - lepS, cores=2)  
brm.foldens.m13<-update(brm.foldens,formula. = ~ .  - lepS - macS, cores=2)  
brm.foldens.m12<-update(brm.foldens,formula. = ~ . - deeS , cores=2)  
brm.foldens.m11<-update(brm.foldens,formula. = ~ . - deeS  - macS, cores=2)  
brm.foldens.m10<-update(brm.foldens,formula. = ~ . - deeS  - lepS, cores=2)  
brm.foldens.m9<-update(brm.foldens,formula. = ~ . - deeS  - lepS - macS, cores=2)  
brm.foldens.m8<-update(brm.foldens,formula. = ~ . - horS  , cores=2)  
brm.foldens.m7<-update(brm.foldens,formula. = ~ . - horS  - macS, cores=2)  
brm.foldens.m6<-update(brm.foldens,formula. = ~ . - horS  - lepS, cores=2)  
brm.foldens.m5<-update(brm.foldens,formula. = ~ . - horS  - lepS - macS, cores=2)  
brm.foldens.m4<-update(brm.foldens,formula. = ~ . - horS - deeS , cores=2)  
brm.foldens.m3<-update(brm.foldens,formula. = ~ . - horS - deeS  - macS, cores=2)  
brm.foldens.m2<-update(brm.foldens,formula. = ~ . - horS - deeS  - lepS, cores=2)  
brm.foldens.m1<-update(brm.foldens,formula. = ~ . - horS - deeS  - lepS - macS, cores=2)  

# Compute the WAIC for each of the 16 models.
waic.foldens<-sapply(list(brm.foldens.m1,brm.foldens.m2,brm.foldens.m3,brm.foldens.m4,brm.foldens.m5,brm.foldens.m6,brm.foldens.m7,brm.foldens.m8,brm.foldens.m9,brm.foldens.m10,brm.foldens.m11,brm.foldens.m12,brm.foldens.m13,brm.foldens.m14,brm.foldens.m15,brm.foldens,brm.foldens.tgp),function(x)waic(x)$estimates[3,1])

# Identify the model with the lowest WAIC.
waic.foldens-min(waic.foldens)

# Select the model with the lowest WAIC and most simple formula.
summary(brm.foldens.m11)


####### Forb cover #######

## vs. Total herbivore activity analysis  
brm.forb.tgp <- brm(forbProp ~  tgpS, 
                    family=Beta(), data=GIS, cores=2, 
                    prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.forb.tgp) 


## vs each herbivore family analysis  
brm.forb <- brm(forbProp ~  horS + deeS  + lepS + macS , 
                family=Beta(), data=GIS, cores=2, 
                prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.forb) 

forb.lm <- lm(forbProp ~  horS + deeS  + lepS + macS, data=GIS) # check multicolinearity
vif(forb.lm) 


brm.forb.m15<-update(brm.forb,formula. = ~ .  - macS, cores=2)  
brm.forb.m14<-update(brm.forb,formula. = ~ .  - lepS, cores=2)  
brm.forb.m13<-update(brm.forb,formula. = ~ .  - lepS - macS, cores=2)  
brm.forb.m12<-update(brm.forb,formula. = ~ . - deeS , cores=2)  
brm.forb.m11<-update(brm.forb,formula. = ~ . - deeS  - macS, cores=2)  
brm.forb.m10<-update(brm.forb,formula. = ~ . - deeS  - lepS, cores=2)  
brm.forb.m9<-update(brm.forb,formula. = ~ . - deeS  - lepS - macS, cores=2)  
brm.forb.m8<-update(brm.forb,formula. = ~ . - horS  , cores=2)  
brm.forb.m7<-update(brm.forb,formula. = ~ . - horS  - macS, cores=2)  
brm.forb.m6<-update(brm.forb,formula. = ~ . - horS  - lepS, cores=2)  
brm.forb.m5<-update(brm.forb,formula. = ~ . - horS  - lepS - macS, cores=2)  
brm.forb.m4<-update(brm.forb,formula. = ~ . - horS - deeS , cores=2)  
brm.forb.m3<-update(brm.forb,formula. = ~ . - horS - deeS  - macS, cores=2)  
brm.forb.m2<-update(brm.forb,formula. = ~ . - horS - deeS  - lepS, cores=2)  
brm.forb.m1<-update(brm.forb,formula. = ~ . - horS - deeS  - lepS - macS, cores=2)  

# Compute the WAIC for each of the 16 models.
waic.forb<-sapply(list(brm.forb.m1,brm.forb.m2,brm.forb.m3,brm.forb.m4,brm.forb.m5,brm.forb.m6,brm.forb.m7,brm.forb.m8,brm.forb.m9,brm.forb.m10,brm.forb.m11,brm.forb.m12,brm.forb.m13,brm.forb.m14,brm.forb.m15,brm.forb,brm.forb.tgp),function(x)waic(x)$estimates[3,1])

# Identify the model with the lowest WAIC.
waic.forb-min(waic.forb)

# Select the model with the lowest WAIC and most simple formula.
summary(brm.forb.m9)


####### Weed cover #######

## vs. Total herbivore activity analysis  
brm.weed.tgp <- brm(weedProp ~  tgpS, 
                    family=Beta(), data=GIS, cores=2, 
                    prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.weed.tgp) 


## vs each herbivore family analysis  
brm.weed <- brm(weedProp ~  horS + deeS  + lepS + macS , 
                family=Beta(), data=GIS, cores=2, 
                prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.weed) 

weed.lm <- lm(weedProp ~  horS + deeS  + lepS + macS, data=GIS) # check multicolinearity
vif(weed.lm) 


brm.weed.m15<-update(brm.weed,formula. = ~ .  - macS, cores=2)  
brm.weed.m14<-update(brm.weed,formula. = ~ .  - lepS, cores=2)  
brm.weed.m13<-update(brm.weed,formula. = ~ .  - lepS - macS, cores=2)  
brm.weed.m12<-update(brm.weed,formula. = ~ . - deeS , cores=2)  
brm.weed.m11<-update(brm.weed,formula. = ~ . - deeS  - macS, cores=2)  
brm.weed.m10<-update(brm.weed,formula. = ~ . - deeS  - lepS, cores=2)  
brm.weed.m9<-update(brm.weed,formula. = ~ . - deeS  - lepS - macS, cores=2)  
brm.weed.m8<-update(brm.weed,formula. = ~ . - horS  , cores=2)  
brm.weed.m7<-update(brm.weed,formula. = ~ . - horS  - macS, cores=2)  
brm.weed.m6<-update(brm.weed,formula. = ~ . - horS  - lepS, cores=2)  
brm.weed.m5<-update(brm.weed,formula. = ~ . - horS  - lepS - macS, cores=2)  
brm.weed.m4<-update(brm.weed,formula. = ~ . - horS - deeS , cores=2)  
brm.weed.m3<-update(brm.weed,formula. = ~ . - horS - deeS  - macS, cores=2)  
brm.weed.m2<-update(brm.weed,formula. = ~ . - horS - deeS  - lepS, cores=2)  
brm.weed.m1<-update(brm.weed,formula. = ~ . - horS - deeS  - lepS - macS, cores=2)  

# Compute the WAIC for each of the 16 models.
waic.weed<-sapply(list(brm.weed.m1,brm.weed.m2,brm.weed.m3,brm.weed.m4,brm.weed.m5,brm.weed.m6,brm.weed.m7,brm.weed.m8,brm.weed.m9,brm.weed.m10,brm.weed.m11,brm.weed.m12,brm.weed.m13,brm.weed.m14,brm.weed.m15,brm.weed,brm.weed.tgp),function(x)waic(x)$estimates[3,1])

# Identify the model with the lowest WAIC.
waic.weed-min(waic.weed)

# Select the model with the lowest WAIC and most simple formula.
summary(brm.weed.tgp)


####### Bare ground cover #######
summary(GIS$baregroundProp)

## vs. Total herbivore activity analysis  
brm.bareground.tgp <- brm(baregroundProp ~  tgpS, 
                          family=Beta(), data=GIS, cores=2, 
                          prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.bareground.tgp) 


## vs each herbivore family analysis  
brm.bareground <- brm(baregroundProp ~  horS + deeS  + lepS + macS , 
                      family=Beta(), data=GIS, cores=2, 
                      prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.bareground) 

bareground.lm <- lm(baregroundProp ~  horS + deeS  + lepS + macS, data=GIS) # check multicolinearity
vif(bareground.lm) 



brm.bareground.m15<-update(brm.bareground,formula. = ~ .  - macS, cores=2)  
brm.bareground.m14<-update(brm.bareground,formula. = ~ .  - lepS, cores=2)  
brm.bareground.m13<-update(brm.bareground,formula. = ~ .  - lepS - macS, cores=2)  
brm.bareground.m12<-update(brm.bareground,formula. = ~ . - deeS , cores=2)  
brm.bareground.m11<-update(brm.bareground,formula. = ~ . - deeS  - macS, cores=2)  
brm.bareground.m10<-update(brm.bareground,formula. = ~ . - deeS  - lepS, cores=2)  
brm.bareground.m9<-update(brm.bareground,formula. = ~ . - deeS  - lepS - macS, cores=2)  
brm.bareground.m8<-update(brm.bareground,formula. = ~ . - horS  , cores=2)  
brm.bareground.m7<-update(brm.bareground,formula. = ~ . - horS  - macS, cores=2)  
brm.bareground.m6<-update(brm.bareground,formula. = ~ . - horS  - lepS, cores=2)  
brm.bareground.m5<-update(brm.bareground,formula. = ~ . - horS  - lepS - macS, cores=2)  
brm.bareground.m4<-update(brm.bareground,formula. = ~ . - horS - deeS , cores=2)  
brm.bareground.m3<-update(brm.bareground,formula. = ~ . - horS - deeS  - macS, cores=2)  
brm.bareground.m2<-update(brm.bareground,formula. = ~ . - horS - deeS  - lepS, cores=2)  
brm.bareground.m1<-update(brm.bareground,formula. = ~ . - horS - deeS  - lepS - macS, cores=2)  

# Compute the WAIC for each of the 16 models.
waic.bareground<-sapply(list(brm.bareground.m1,brm.bareground.m2,brm.bareground.m3,brm.bareground.m4,brm.bareground.m5,brm.bareground.m6,brm.bareground.m7,brm.bareground.m8,brm.bareground.m9,brm.bareground.m10,brm.bareground.m11,brm.bareground.m12,brm.bareground.m13,brm.bareground.m14,brm.bareground.m15,brm.bareground,brm.bareground.tgp),function(x)waic(x)$estimates[3,1])

# Identify the model with the lowest WAIC.
waic.bareground-min(waic.bareground)

# Select the model with the lowest WAIC and most simple formula.
summary(brm.bareground.m3)


####### Soil Compaction #######

## vs. Total herbivore activity analysis  
# brm.soilcomp.tgp <- brm(soil.mean ~  tgpS + SoilMoisture, 
#                        family=Gamma(link = "log"), data=soils, cores=2, 
#                        prior=bprior,control=list(adapt_delta = 0.9)) 
### model diagnostics suggest non-linear relationship
brm.soilcomp.tgp <- brm(soil.mean ~  bs(tgpS, df = 3)  + SoilMoisture, 
                        family=Gamma(link = "log"), data=soils, cores=2, 
                        prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.soilcomp.tgp) 


## vs each herbivore family analysis  
brm.soilcomp <- brm(soil.mean ~  horS + deeS  + lepS + macS  
                    + SoilMoisture, family=Gamma(link = "log"), data=soils,  
                    cores=2, prior=bprior,control=list(adapt_delta = 0.9))
summary(brm.soilcomp) 

soil.lm <- lm(soil.mean ~  horS + deeS  + lepS + macS, data=soils) # check multicolinearity
vif(soil.lm) 



brm.soilcomp.m15<-update(brm.soilcomp,formula. = ~ .  - macS, cores=2)  
brm.soilcomp.m14<-update(brm.soilcomp,formula. = ~ .  - lepS, cores=2)  
brm.soilcomp.m13<-update(brm.soilcomp,formula. = ~ .  - lepS - macS, cores=2)  
brm.soilcomp.m12<-update(brm.soilcomp,formula. = ~ . - deeS , cores=2)  
brm.soilcomp.m11<-update(brm.soilcomp,formula. = ~ . - deeS  - macS, cores=2)  
brm.soilcomp.m10<-update(brm.soilcomp,formula. = ~ . - deeS  - lepS, cores=2)  
brm.soilcomp.m9<-update(brm.soilcomp,formula. = ~ . - deeS  - lepS - macS, cores=2)  
brm.soilcomp.m8<-update(brm.soilcomp,formula. = ~ . - horS  , cores=2)  
brm.soilcomp.m7<-update(brm.soilcomp,formula. = ~ . - horS  - macS, cores=2)  
brm.soilcomp.m6<-update(brm.soilcomp,formula. = ~ . - horS  - lepS, cores=2)  
brm.soilcomp.m5<-update(brm.soilcomp,formula. = ~ . - horS  - lepS - macS, cores=2)  
brm.soilcomp.m4<-update(brm.soilcomp,formula. = ~ . - horS - deeS , cores=2)  
brm.soilcomp.m3<-update(brm.soilcomp,formula. = ~ . - horS - deeS  - macS, cores=2)  
brm.soilcomp.m2<-update(brm.soilcomp,formula. = ~ . - horS - deeS  - lepS, cores=2)  
brm.soilcomp.m1<-update(brm.soilcomp,formula. = ~ . - horS - deeS  - lepS - macS, cores=2)  

# Compute the WAIC for each of the 16 models.
waic.soilcomp<-sapply(list(brm.soilcomp.m1,brm.soilcomp.m2,brm.soilcomp.m3,brm.soilcomp.m4,brm.soilcomp.m5,brm.soilcomp.m6,brm.soilcomp.m7,brm.soilcomp.m8,brm.soilcomp.m9,brm.soilcomp.m10,brm.soilcomp.m11,brm.soilcomp.m12,brm.soilcomp.m13,brm.soilcomp.m14,brm.soilcomp.m15,brm.soilcomp,brm.soilcomp.tgp),function(x)waic(x)$estimates[3,1])

# Identify the model with the lowest WAIC.
waic.soilcomp-min(waic.soilcomp)

# Select the model with the lowest WAIC and most simple formula.
summary(brm.soilcomp.m9)


## Grazing Impacts Model Diagnostics - total herbivore activity and each herbivore family
# For information on the diagnostics presented below, refer to vignette("DHARMa", package="DHARMa").


### Groundstorey Height
Mheight.tgp <- glmmTMB(grassheight.mean ~  tgpS, 
                       family=Gamma(link = "log"), data=grassheight.noNA) # brm.height.tgp
plot(residuals(Mheight.tgp))
testResiduals(Mheight.tgp)
height.tgp.sim <- simulateResiduals(Mheight.tgp, n=250, refit = F, plot = F, seed = 123) 
plot(height.tgp.sim) 

Mheight <- glmmTMB(grassheight.mean ~  horS + lepS, 
                   family=Gamma(link = "log"), data=grassheight.noNA)# (brm.grassheight.m11)
plot(residuals(Mheight))
testResiduals(Mheight)
height.sim <- simulateResiduals(Mheight, n=250, refit = F, plot = F, seed = 123) 
plot(height.sim) 

### Foliage density 
Mfoldens.tgp <- glmmTMB(fol5.mean ~  tgpS, 
                        family=Gamma(link = "log"), data=GIS)
plot(residuals(Mfoldens.tgp))
testResiduals(Mfoldens.tgp)
foldens.tgp.sim <- simulateResiduals(Mfoldens.tgp, n=250, refit = F, plot = F, seed = 123) 
plot(foldens.tgp.sim) 

Mfoldens <- glmmTMB(fol5.mean ~  horS + lepS, 
                    family=Gamma(link = "log"), data=GIS)# (brm.foldens.m11)
plot(residuals(Mfoldens))
testResiduals(Mfoldens)
foldens.sim <- simulateResiduals(Mfoldens, n=250, refit = F, plot = F, seed = 123) 
plot(foldens.sim) 


### Forb cover 
Mforb.tgp <- glmmTMB(forbProp ~  tgpS, 
                     family=beta_family(), data=GIS)
plot(residuals(Mforb.tgp))
testResiduals(Mforb.tgp)
forb.tgp.sim <- simulateResiduals(Mforb.tgp, n=250, refit = F, plot = F, seed = 123) 
plot(forb.tgp.sim)  

Mforb <- glmmTMB(forbProp ~  horS, 
                 family=beta_family(), data=GIS) #(brm.forb.m9)
plot(residuals(Mforb))
testResiduals(Mforb)
forb.sim <- simulateResiduals(Mforb, n=250, refit = F, plot = F, seed = 123) 
plot(forb.sim)  

### Weed cover  
Mweed.tgp <- glmmTMB(weedProp ~  tgpS, 
                     family=beta_family(), data=GIS)
plot(residuals(Mweed.tgp))
testResiduals(Mweed.tgp)
weed.tgp.sim <- simulateResiduals(Mweed.tgp, n=250, refit = F, plot = F, seed = 123) 
plot(weed.tgp.sim)  

Mweed <- glmmTMB(weedProp ~  horS , 
                 family=beta_family(), data=GIS) # summary(brm.weed.m9)
plot(residuals(Mweed))
testResiduals(Mweed)
weed.sim <- simulateResiduals(Mweed, n=250, refit = F, plot = F, seed = 123) 
plot(weed.sim)  

### Bare ground cover 
Mbareground.tgp <- glmmTMB(baregroundProp ~  tgpS, 
                           family=beta_family(), data=GIS)
plot(residuals(Mbareground.tgp))
testResiduals(Mbareground.tgp)
bareground.tgp.sim <- simulateResiduals(Mbareground.tgp, n=250, refit = F, plot = F, seed = 123) 
plot(bareground.tgp.sim) 

Mbareground <- glmmTMB(baregroundProp ~  lepS, 
                       family=beta_family(), data=GIS) # summary(brm.bareground.m3)
summary(Mbareground)
plot(residuals(Mbareground))
testResiduals(Mbareground)
bareground.sim <- simulateResiduals(Mbareground, n=250, refit = F, plot = F, seed = 123) 
plot(bareground.sim) 


### Soil compaction
Msoil.tgp <- glmmTMB(soil.mean ~  bs(tgpS, df = 3)
                     + SoilMoisture, family=Gamma(link = "log"), data=soils)
plot(residuals(Msoil.tgp))
testResiduals(Msoil.tgp)
soil.tgp.sim <- simulateResiduals(Msoil.tgp, n=250, refit = F, plot = F, seed = 123) 
plot(soil.tgp.sim) 

Msoil <- glmmTMB(soil.mean ~  horS 
                 + SoilMoisture, family=Gamma(link = "log"), data=soils) # summary(brm.soilcomp.m9)
plot(residuals(Msoil))
testResiduals(Msoil)
soil.sim <- simulateResiduals(Msoil, n=250, refit = F, plot = F, seed = 123) 
plot(soil.sim) 



## Analyses of residuals against predictor variables

### Height 
plot(residuals(Mheight.tgp)  ~  grassheight.noNA$tgpS)
plot(residuals(Mheight) ~ grassheight.noNA$deeS)
plot(residuals(Mheight) ~ grassheight.noNA$lepS)

### Foliage density 
plot(residuals(Mfolden.tgps)  ~  GIS$tgpS)
plot(residuals(Mfoldens)  ~  GIS$deeS)
plot(residuals(Mfoldens)  ~  GIS$macS)

### Forb cover 
plot(residuals(Mforb.tgp)  ~  GIS$tgpS)
plot(residuals(Mforb)  ~  GIS$horS)

### Weed cover  
plot(residuals(Mweed.tgp)  ~  GIS$tgpS)
plot(residuals(Mweed)  ~  GIS$horS)

### Bare ground cover 
plot(residuals(Mbareground.tgp)  ~  GIS$tgpS)
plot(residuals(Mbareground)  ~  GIS$lepS)

### Soil compaction
plot(residuals(Msoil.tgp)  ~  soils$tgpS)
plot(residuals(Msoil)  ~  soils$horS)




########################## PLOTS ##############################################
ce.tgp<-conditional_effects(brmtgpC.nb.m41)
ce.tgp.ele <- ce.tgp$eleS
ce.tgp.vegtype <- ce.tgp$vegtype
ce.tgp.vegtype$vegtype <- as.character(ce.tgp.vegtype$vegtype)
ce.tgp.vegtype$vegtype[ce.tgp.vegtype$vegtype == "shrub.grass"] <- "sh"
ce.tgp.vegtype$vegtype[ce.tgp.vegtype$vegtype == "sedge.grass"] <- "gr"
ce.tgp.vegtype$vegtype[ce.tgp.vegtype$vegtype == "tussock.grass"] <- "tu"
ce.tgp.ele$eleM = ce.tgp.ele$eleS*sd(GIS$ele.mean)+mean(GIS$ele.mean)

TtgpC_ele <- ggplot(ce.tgp.ele, aes(x=eleM, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Elevation (m a.s.l.)", y="Total herbivore activity (index)") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
TtgpC_ele
ggsave("ele_tgp.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("ele_tgp.png", width = 35, height = 35, units = "mm", dpi=900)

TtgpC_vegtype <- ggplot(ce.tgp.vegtype, aes(x=vegtype, y=estimate__)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), stat = 'identity')  +
  labs(x="Grassland type", y="Total herbivore activity (index)") +
  theme_classic() + scale_y_continuous(limits=c(0, 60),breaks=c(20, 40, 60)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
TtgpC_vegtype
ggsave("veg_tgp.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("veg_tgp.png", width = 35, height = 35, units = "mm", dpi=900)


emmeans(brmtgpC.nb.m41,pairwise~vegtype,type="response") 


ce.hor<-conditional_effects(brmhorPA.b.m35)
ce.hor.ele <- ce.hor$eleS
ce.hor.ele$eleM = ce.hor.ele$eleS*sd(GIS$ele.mean)+mean(GIS$ele.mean)
ce.hor.DistWtr <- ce.hor$DistWtrS
ce.hor.DistWtr$DistWtrM = exp(ce.hor.DistWtr$DistWtrS*sd(GIS$DistWtrL)+mean(GIS$DistWtrL))


horPA_ele <- ggplot(ce.hor.ele, aes(x=eleM, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Elevation (m a.s.l.)", y="Probability of horse presence") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits=c(0, 1),breaks=c(0.2, 0.4, 0.6, 0.8)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
horPA_ele
ggsave("ele_horPA.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("ele_horPA.png", width = 35, height = 35, units = "mm", dpi=900)


horPA_DistWtr <- ggplot(ce.hor.DistWtr, aes(x=DistWtrM, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Nearest waterbody (m)", y="Probability of horse presence") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits=c(0, 1),breaks=c(0.2, 0.4, 0.6, 0.8)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
horPA_DistWtr
ggsave("DistWtr_horPA.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("DistWtr_horPA.png", width = 35, height = 35, units = "mm", dpi=900)


ce.dee<-conditional_effects(brmdeePA.b.m43)
ce.dee.DistWtr <- ce.dee$DistWtrS
ce.dee.DistWtr$DistWtr = exp(ce.dee.DistWtr$DistWtrS*sd(GIS$DistWtrL)+mean(GIS$DistWtrL))

ce.dee.ele <- ce.dee$eleS
ce.dee.ele$eleM = ce.dee.ele$eleS*sd(GIS$ele.mean)+mean(GIS$ele.mean)

TdeeC_ele <- ggplot(ce.dee.ele, aes(x=eleM, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Elevation (m a.s.l.)", y="Probability of deer presence") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits=c(0, 1),breaks=c(0.2, 0.4, 0.6, 0.8)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
TdeeC_ele
ggsave("ele_deePA.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("ele_deePA.png", width = 35, height = 35, units = "mm", dpi=900)


TdeeCPA_DistWtr <- ggplot(ce.dee.DistWtr, aes(x=DistWtr, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity') +
  labs(x="Nearest waterbody (m)", y="Probability of deer presence") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits=c(0, 1),breaks=c(0.2, 0.4, 0.6, 0.8)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
TdeeCPA_DistWtr

ggsave("Wtr_deePA.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("Wtr_deePA.png", width = 35, height = 35, units = "mm", dpi=900)

ce.dee.vegtype <- ce.dee$vegtype
ce.dee.vegtype$vegtype <- as.character(ce.dee.vegtype$vegtype)
ce.dee.vegtype$vegtype[ce.dee.vegtype$vegtype == "shrub.grass"] <- "sh"
ce.dee.vegtype$vegtype[ce.dee.vegtype$vegtype == "sedge.grass"] <- "gr"
ce.dee.vegtype$vegtype[ce.dee.vegtype$vegtype == "tussock.grass"] <- "tu"


TdeeC_vegtype <- ggplot(ce.dee.vegtype, aes(x=vegtype, y=estimate__)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), stat = 'identity')  +
  labs(x="Grassland type", y="Probability of deer presence") +
  theme_classic() + scale_y_continuous(expand = c(0, 0), limits=c(0, 1),breaks=c(0.2, 0.4, 0.6, 0.8)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
TdeeC_vegtype 
ggsave("veg_deePA.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("veg_deePA.png", width = 35, height = 35, units = "mm", dpi=900)

emmeans(brmdeePA.b.m43,pairwise~vegtype,type="response") 


ce.lep<-conditional_effects(brmlepPA.b.m3)
ce.lep.DistWtr <- ce.lep$DistWtrS
ce.lep.DistWtr$DistWtr = exp(ce.lep.DistWtr$DistWtrS*sd(GIS$DistWtrL)+mean(GIS$DistWtrL))

lepPA_DistWtr <- ggplot(ce.lep.DistWtr, aes(x=DistWtr, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Nearest waterbody (m)", y="Probability of rabbit/hare presence") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0),limits=c(0, 1),breaks=c(0.2, 0.4, 0.6, 0.8)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
lepPA_DistWtr 
ggsave("Wtr_lepPA.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("Wtr_lepPA.png", width = 35, height = 35, units = "mm", dpi=900)

ce.mac<-conditional_effects(brmmacPA.b.m41)
ce.mac.ele <- ce.mac$eleS
ce.mac.ele$eleM = ce.mac.ele$eleS*sd(GIS$ele.mean)+mean(GIS$ele.mean)

TmacC_ele <- ggplot(ce.mac.ele, aes(x=eleM, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Elevation (m a.s.l.)", y="Probability of kangaroo/wallaby presence") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits=c(0, 1),breaks=c(0.2, 0.4, 0.6, 0.8)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
TmacC_ele
ggsave("ele_macPA.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("ele_macPA.png", width = 35, height = 35, units = "mm", dpi=900)

ce.mac.vegtype <- ce.mac$vegtype
ce.mac.vegtype$vegtype <- as.character(ce.mac.vegtype$vegtype)
ce.mac.vegtype$vegtype[ce.mac.vegtype$vegtype == "shrub.grass"] <- "sh"
ce.mac.vegtype$vegtype[ce.mac.vegtype$vegtype == "sedge.grass"] <- "gr"
ce.mac.vegtype$vegtype[ce.mac.vegtype$vegtype == "tussock.grass"] <- "tu"

TmacC_vegtype <- ggplot(ce.mac.vegtype, aes(x=vegtype, y=estimate__)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), stat = 'identity')  +
  labs(x="Grassland type", y="Probability of kangaroo/wallaby presence") +
  theme_classic() + scale_y_continuous(expand = c(0, 0), limits=c(0, 1),breaks=c(0.2, 0.4, 0.6, 0.8)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
TmacC_vegtype
ggsave("veg_macPA.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("veg_macPA.png", width = 35, height = 35, units = "mm", dpi=900)

emmeans(brmmacPA.b.m43,pairwise~vegtype,type="response") 



### FIGURE 4
### Vegetation height

ce.height<-conditional_effects(brm.grassheight.m11)
ce.height.hor <- ce.height$horS
ce.height.lep <- ce.height$lepS
ce.height.hor$ThorC = exp(ce.height.hor$horS*sd(grassheight.noNA$horCL)+mean(grassheight.noNA$horCL))
ce.height.lep$TlepC = exp(ce.height.lep$lepS*sd(grassheight.noNA$lepCL)+mean(grassheight.noNA$lepCL))

height_hor <- ggplot(ce.height.hor, aes(x=ThorC, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Horse activity (index)", y="Vegetation height (cm)") +
  theme_classic() + scale_x_continuous(expand = c(0, 0), breaks = c(5,10,15,20), limits = c(2,22)) + scale_y_continuous(expand = c(0, 0), breaks = c(10,15,20,25), limits=c(10, 25)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
height_hor
ggsave("height_hor.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("height_hor.png", width = 35, height = 35, units = "mm", dpi=900)

height_lep <- ggplot(ce.height.lep, aes(x=TlepC, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Rabbit/hare activity (index)", y="Vegetation height (cm)") +
  theme_classic() + scale_x_continuous(expand = c(0, 0), breaks = c(5,10,15,20), limits = c(2,22)) + scale_y_continuous(expand = c(0, 0), breaks = c(10,15,20,25), limits=c(10, 25)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
height_lep
ggsave("height_lep.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("height_lep.png", width = 35, height = 35, units = "mm", dpi=900)



### Foliage density

foldens_hor <- ggplot(ce.foldens.hor, aes(x=ThorC, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Horse activity (index)", y="Foliage density <50cm") +
  theme_classic() + scale_x_continuous(expand = c(0, 0), limits = c(2,22)) + scale_y_continuous(expand = c(0, 0), breaks = c(30,40,50,60), limits=c(29,52)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
foldens_hor 
ggsave("foldens_hor.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("foldens_hor.png", width = 35, height = 35, units = "mm", dpi=900)

foldens_lep <- ggplot(ce.foldens.lep, aes(x=TlepC, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Rabbit/hare activity (index)", y="Foliage density <50cm") +
  theme_classic() + scale_x_continuous(expand = c(0, 0), limits = c(2,22)) + scale_y_continuous(expand = c(0, 0), breaks = c(30,40,50,60), limits=c(29,52)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
foldens_lep 
ggsave("foldens_lep.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("foldens_lep.png", width = 35, height = 35, units = "mm", dpi=900)


### Forb cover

ce.forb<-conditional_effects(brm.forb.m9)
ce.forb.hor <- ce.forb$horS
ce.forb.hor$ThorC = exp(ce.forb.hor$horS*sd(GIS$horCL)+mean(GIS$horCL))

forb_hor <- ggplot(ce.forb.hor, aes(x=ThorC, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Horse activity (index)", y="Proportion of forb cover") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = c(0.1, 0.2, 0.3), limits=c(0.1,0.31)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
forb_hor
ggsave("forb_hor.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("forb_hor.png", width = 35, height = 35, units = "mm", dpi=900)



### Weed cover

ce.weed.tgp<-conditional_effects(brm.weed.tgp)
ce.weed.tgpS <- ce.weed.tgp$tgpS
ce.weed.tgpS$TtgpC = exp(ce.weed.tgpS$tgpS*sd(GIS$tgpCL)+mean(GIS$tgpCL))

weed_tgp <- ggplot(ce.weed.tgpS, aes(x=TtgpC, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Total herbivore activity (index)", y="Proportion of weed cover") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = c(0,0.1, 0.2), limits=c(0,0.21)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.1,hjust=1, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
weed_tgp
ggsave("weed_tgp.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("weed_tgp.png", width = 35, height = 35, units = "mm", dpi=900)



### Bare Ground

ce.bareground<-conditional_effects(brm.bareground.m3)
ce.bareground.lep <- ce.bareground$lepS
ce.bareground.lep$TlepC = exp(ce.bareground.lep$lepS*sd(GIS$lepCL)+mean(GIS$lepCL))

bareground_lep <- ggplot(ce.bareground.lep, aes(x=TlepC, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Rabbit/hare activity (index)", y="Proportion of bare ground") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.01, 0.02), limits=c(0,0.02)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))

bareground_lep  
ggsave("bareground_lep.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("bareground_lep.png", width = 35, height = 35, units = "mm", dpi=900)


### Soil Compaction


ce.soilcomp<-conditional_effects(brm.soilcomp.m9)
ce.soilcomp.hor <- ce.soilcomp$horS
ce.soilcomp.hor$ThorC = exp(ce.soilcomp.hor$horS*sd(soils$horCL)+mean(soils$horCL))

soilcomp_hor <- ggplot(ce.soilcomp.hor, aes(x=ThorC, y=estimate__)) +
  geom_line()  +
  geom_smooth(aes(ymin=lower__, ymax=upper__), stat='identity')+
  labs(x="Horse activity (index)", y="Soil compaction (kPa)") +
  theme_classic() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = c(0,1,2,3), limits = c(0.5,2.5)) +
  theme(legend.position="none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8,vjust=0.5,margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=8,vjust=0.5,angle = 90,margin = margin(t = 0, r = 5, b = 0, l = 0)))
soilcomp_hor
ggsave("soilcomp_hor.pdf", width = 35, height = 35, units = "mm", dpi=900)
ggsave("soilcomp_hor.png", width = 35, height = 35, units = "mm", dpi=900)



############### FIGURES #########################

## Figure 2
# not made in R

## Figure 3 
plot_grid(horPA_ele, horPA_DistWtr, TdeeC_ele, TdeeC_vegtype,TdeeCPA_DistWtr, 
         lepPA_DistWtr,  TmacC_ele,TmacC_vegtype,
         TtgpC_ele,  TtgpC_vegtype, 
         ncol = 4,
          labels = "auto",label_size=10,
          label_x = 0.3, label_y = 1)

ggsave("Figure_3.pdf", width = 180, height = 140, units = "mm", dpi=600)


## Figure 4
plot_grid(height_hor, height_lep, foldens_hor, foldens_lep,  
          forb_hor, weed_tgp, bareground_lep, soilcomp_hor,
          ncol = 4,
          labels = "auto", label_size = 10,
          label_x = 0.3, label_y = 1)
ggsave("Figure_4.pdf", width = 180, height = 90, units = "mm", dpi=600)


## Supplementary Information Plots
## 1.1a

distw <- ggplot(GIS, aes(DistWtr.min)) + geom_histogram() + 
labs(x="Nearest waterbody (m)") + 
  theme_classic()  + 
  theme(legend.position="none",
        text = element_text(size=8))

distt <- ggplot(GIS, aes(DistWood.min)) + geom_histogram()+
labs(x="Nearest woodland (m)") +
  theme_classic()  + 
  theme(legend.position="none",
        text = element_text(size=8))

distr <- ggplot(GIS, aes(DistRd.min)) + geom_histogram() + 
  theme_classic() + 
labs(x="Nearest track (m)") + 
  theme(legend.position="none",
        text = element_text(size=8))

ele  <- ggplot(GIS, aes(ele.mean)) + geom_histogram() +
  labs(x="Elevation (m)") +
  theme_classic() + 
  theme(legend.position="none",
        text = element_text(size=8))

TWI <- ggplot(GIS, aes(TWImean)) + geom_histogram() +
  labs(x="Landscape position (TWI)") + 
  theme_classic()  + 
  theme(legend.position="none",
        text = element_text(size=8))

long <- ggplot(GIS, aes(longGDA94)) + geom_histogram()+ 
  labs(x="Longitude")+ theme_classic() + 
  theme(legend.position="none",
        text = element_text(size=8))

lat <- ggplot(GIS, aes(latGDA94)) + geom_histogram() +
    labs(x="Latitude")+ theme_classic()  + 
  theme(legend.position="none",
        text = element_text(size=8))

plot_grid(distw, distt, distr, ele,  TWI, NULL,lat, long, ncol=2)
ggsave("Figure_S1_1a.tiff", width = 120, height = 240, units = "mm", dpi=300)


## Supp Info 1.1b

hor  <- ggplot(GIS, aes(ThorC)) + geom_histogram() + 
  labs(x= "Horse Grazing Activity (Index)") +
  theme_classic() + 
  theme(legend.position="none",
        text = element_text(size=8))
dee <- ggplot(GIS, aes(TdeeC)) + geom_histogram() +
  labs(x= "Deer Grazing Activity (Index)") +
  theme_classic() + 
  theme(legend.position="none",
        text = element_text(size=8))
lep <- ggplot(GIS, aes(TlepC)) + geom_histogram() +
  labs(x= "Rabbit & Hare Grazing Activity (Index)") +
  theme_classic()  + 
  theme(legend.position="none",
        text = element_text(size=8),
        axis.title.x = element_text(hjust=1)) + 
  scale_y_continuous(
    labels = scales::number_format(accuracy = 1))
mac <- ggplot(GIS, aes(TmacC)) + geom_histogram() +
      labs(x= "Kangaroo & Wallaby Grazing Activity (Index)") +
  theme_classic()  + 
  theme(legend.position="none",
        text = element_text(size=8),
        axis.title.x = element_text(hjust=1))
tgp <- ggplot(GIS, aes(TtgpC)) + geom_histogram() +  
  labs(x= "Total Grazing Activity (Index)") +
  theme_classic()  + 
  theme(legend.position="none",
        text = element_text(size=8))
plot_grid(hor, dee, lep, mac, tgp, ncol=2)
ggsave("Figure_S1_1b.tiff", width = 120, height = 180, units = "mm", dpi=300)


## Supp Info 1.2

plot_grid(vdistw, vdistt, vdistr, vele, vTWI, fdistw, fdistt, fdistr, fele, fTWI,
          odistw, odistt, odistr, oele, oTWI, wdistw, wdistt, wdistr, wele, wTWI,
          bdistw, bdistt, bdistr, bele, bTWI, sdistw, sdistt, sdistr, sele, sTWI,
          ncol=5)
ggsave("Figure_S1_2.png", width = 22, height = 26.4, units = "cm", dpi=900)


############### END ####################################################
