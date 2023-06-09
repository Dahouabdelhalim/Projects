##### WENDLANDT ET AL 2022 OECOLOGIA #####
##### "Wild legumes maintain beneficial soil rhizobia populations despite decades of nitrogen deposition"
##### DOI: 10.1007/s00442-022-05116-9


##### Code for Figures and Analyses ############################

##### Notes for using this code #####

## Terms: in the manuscript, we use the terms "NonNitrogen.PC1" and "NonNitrogen.PC2", but in the code file, those terms are called "General.PC1" and "General.PC2", respectively

## SETWD AND LIBRARIES ####

setwd("C:/Users/Camille/Documents/Cami_Lab_Backup/2_Whole Soil Inoc/Dryad data and code")

library(plyr)
library(dplyr)
library(tidyr)
library(reshape)
library(lme4)
library(gmodels)
library(gplots)
library(visreg)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(emmeans)
library(DHARMa)
library(factoextra)
library(maps)
library(graphics)
library(mapplots)
library(usmap)
library(rgdal)

sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19043)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] rgdal_1.5-28     sp_1.4-6         usmap_0.5.2      mapplots_1.5.1   maps_3.4.0      
# [6] factoextra_1.0.7 DHARMa_0.4.4     emmeans_1.7.1-1  gridExtra_2.3    reshape2_1.4.4  
# [11] ggplot2_3.3.5    visreg_2.7.0     gplots_3.1.1     gmodels_2.18.1   lme4_1.1-27.1   
# [16] Matrix_1.3-4     reshape_0.8.8    tidyr_1.1.4      dplyr_1.0.7      plyr_1.8.6      



## USEFUL FUNCTIONS ####


# function to remove NA's from select columns
complete_fun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

CI <- function(x) {mean(x)+sd(x)/sqrt(length(x))}
CI2 <- function(x) {mean(x)-sd(x)/sqrt(length(x))}




###   IMPORT SOIL DATA (D1)   ############################################
    
D1 <- read.csv("WSI_soil_data.csv", header = TRUE)

# set data types
    D1$Year <- as.factor(D1$Year)
    D1$Source <- as.factor(D1$Sample)
    D1$Metric <- as.factor(D1$Metric)
    D1$Population <- as.factor(D1$Population)
    D1$Value <- as.numeric(D1$Value)
    D1$Notes <- as.factor(D1$Notes)
    
# pivot data so each population has a column of data
temp <- pivot_wider(D1, names_from = "Population", values_from = "Value")

# make a label to use for selecting particular datapoints
temp$Label <- paste(temp$Metric, temp$Year, temp$Source, sep = ".")    


# *** PCA for Nitrogen Soil Traits  ####
    
# make a small dataframe for nitrogen ordination
Nit <- temp[which(temp$Label=="NO3.N.2015.Soil" |
                  temp$Label=="NH4.N.2015.Soil" |
                  temp$Label=="TKN.2015.Soil" |
                  temp$Label=="Mineral.N.2013.Soil" |
                  temp$Label=="Total.N.2013.Soil" |
                  temp$Label=="Dry.N.Dep.2013.Air"),]
Nit <- as.data.frame(Nit)
row.names(Nit) <- Nit$Label
Nit <- select(Nit, c("Anz", "BMR", "Cla", "Gri", "UCR", "Yuc"))
head(Nit)

# transpose the nitrogen dataframe so that populations are now rows and each N measurement is a column
Nit <- as.data.frame(t(as.matrix(Nit)))

# add a  column to calculate organic nitrogen
Nit$Organic.N.2015.Soil <- Nit$TKN.2015.Soil - Nit$NH4.N.2015.Soil

# remove TKN column
Nit <- select(Nit, -"TKN.2015.Soil")


# look at dataframe "Nit", ready for PCA
Nit
#     NO3.N.2015.Soil NH4.N.2015.Soil Mineral.N.2013.Soil Total.N.2013.Soil Dry.N.Dep.2013.Air Organic.N.2015.Soil
# Anz               9             3.9                2.02              0.01               1.68               239.1
# BMR               2             4.8                4.08              0.01               0.34               122.2
# Cla              29             4.3               10.81              0.11               7.42              1457.7
# Gri              21             7.8                6.70              0.04               4.32              1396.2
# UCR              32             3.7               20.47              0.07               7.67               630.3
# Yuc               7             4.5                7.04              0.03               1.61               231.5


# do PCA
pca.nit = princomp(Nit, cor = TRUE)
summary(pca.nit)
# Importance of components:
#                           Comp.1    Comp.2    Comp.3     Comp.4       Comp.5       Comp.6
# Standard deviation     2.0192046 1.2060361 0.6198714 0.28949685 1.551178e-02 8.210072e-09
# Proportion of Variance 0.6795312 0.2424205 0.0640401 0.01396807 4.010253e-05 1.123421e-17
# Cumulative Proportion  0.6795312 0.9219517 0.9859918 0.99995990 1.000000e+00 1.000000e+00

    
loadings(pca.nit)
# the soil N measures vary positively with the PC1 score

# Loadings:
#                     Comp.1 Comp.2 Comp.3 Comp.4 Comp.5 Comp.6
# NO3.N.2015.Soil      0.487         0.139  0.554  0.602  0.270
# NH4.N.2015.Soil            -0.794  0.456 -0.221  0.223 -0.253
# Mineral.N.2013.Soil  0.401  0.335  0.661 -0.400 -0.235  0.274
# Total.N.2013.Soil    0.463        -0.485 -0.622  0.375 -0.146
# Dry.N.Dep.2013.Air   0.492                0.304 -0.357 -0.731
# Organic.N.2015.Soil  0.382 -0.502 -0.315        -0.514  0.483
# 
#                Comp.1 Comp.2 Comp.3 Comp.4 Comp.5 Comp.6
# SS loadings     1.000  1.000  1.000  1.000  1.000  1.000
# Proportion Var  0.167  0.167  0.167  0.167  0.167  0.167
# Cumulative Var  0.167  0.333  0.500  0.667  0.833  1.000

  
pca.nit$scores
#         Comp.1     Comp.2      Comp.3      Comp.4       Comp.5        Comp.6
# Anz -1.8913066  0.4824804 -0.41209568  0.53649065  0.004558677  2.220446e-16
# BMR -2.3671893  0.1524744  0.08041312 -0.24425876 -0.026328574 -3.053113e-15
# Cla  2.6835997 -0.1094637 -1.06971239 -0.12357906 -0.003959475 -5.551115e-16
# Gri  0.5659465 -2.4590993  0.49847635  0.08872139  0.002370418  9.992007e-16
# UCR  2.4086158  1.4876356  0.85689750  0.09391341 -0.003078556 -5.551115e-16
# Yuc -1.3996661  0.4459726  0.04602110 -0.35128764  0.026437509  2.498002e-15


# save Nitrogen.PC1 values in a new dataframe to merge into plant data
    nitrogenPC <- as.data.frame(pca.nit$scores[,1])
    nitrogenPC$Inoculum_source <- row.names(nitrogenPC)
    names(nitrogenPC) <- c("Nitrogen.PC1", "Inoculum_source")
    nitrogenPC
    #     Nitrogen.PC1 Inoculum_source
    # Anz   -1.8913066             Anz
    # BMR   -2.3671893             BMR
    # Cla    2.6835997             Cla
    # Gri    0.5659465             Gri
    # UCR    2.4086158             UCR
    # Yuc   -1.3996661             Yuc
    

 
# *** PCA for Non-Nitrogen Soil Traits ####

# make a small dataframe for general soil trait ordination
Gen <- temp[which(temp$Label=="pH.2015.Soil" |
                temp$Label=="CEC.2015.Soil" |
                temp$Label=="P1.Weak.Bray.2015.Soil" |
                temp$Label=="NaHCO3.P.2015.Soil" |
                temp$Label=="Potassium.2015.Soil" |
                temp$Label=="Calcium.2015.Soil" |
                temp$Label=="Magnesium.2015.Soil" |
                temp$Label=="SO4.S.2015.Soil" |
                temp$Label=="Sodium.2015.Soil"),]
Gen <- as.data.frame(Gen)
row.names(Gen) <- Gen$Label
Gen <- select(Gen, c("Anz", "BMR", "Cla", "Gri", "UCR", "Yuc"))

# transpose the dataframe so that populations are now rows and each soil measurement is a column
Gen <- as.data.frame(t(as.matrix(Gen)))

# choose just one phosphorus test based on soil pH
# see this lecture, slide 31
# http://www.agronext.iastate.edu/soilfertility/presentations/mbotest.pdf
Gen$Phos.2015.Soil <- ifelse(Gen$pH.2015.Soil <= 7.3, Gen$P1.Weak.Bray.2015.Soil, Gen$NaHCO3.P.2015.Soil)

# remove redundant phosphorus columns
Gen <- select(Gen, -c("P1.Weak.Bray.2015.Soil", "NaHCO3.P.2015.Soil"))
    
# princomp can only be used with more units (rows) than variables (columns)
# but here I have 6 units (sites) and 8 variables (soil measures)
# can get around this by using prcomp instead (Q-mode instead of R-mode)

pca.gen = prcomp(Gen, scale = TRUE)
summary(pca.gen)
# Importance of components:
#   PC1    PC2    PC3     PC4     PC5       PC6
# Standard deviation     1.8542 1.6509 1.1024 0.68407 0.39173 1.358e-16
# Proportion of Variance 0.4297 0.3407 0.1519 0.05849 0.01918 0.000e+00
# Cumulative Proportion  0.4297 0.7704 0.9223 0.98082 1.00000 1.000e+00

# PC1 explains 43.0% of the variation in soil traits

# what are the axes? eigenvectors
# this means that along the first PC, all the soil measures decrease
# so sites with high values for PC1 will have low values for all these soil traits
# the traits most dramatically associated with PC1 are CEC, Ca, and Mg

pca.gen$rotation
#                             PC1        PC2         PC3         PC4         PC5         PC6
# pH.2015.Soil        -0.01369422 -0.2319535  0.79522857 -0.42307715 -0.04911388 -0.33803209
# CEC.2015.Soil       -0.51549785  0.1637705 -0.08555482 -0.05502381  0.14014372  0.03813420
# Potassium.2015.Soil -0.21561126 -0.5097061 -0.31969291  0.11262217  0.11377475 -0.64412198
# Calcium.2015.Soil   -0.49769028 -0.1536362  0.07412820 -0.32767233  0.42101139  0.46136989
# Magnesium.2015.Soil -0.48309623  0.2598070 -0.10445262 -0.02945443  0.01484607 -0.36106933
# SO4.S.2015.Soil     -0.09848324 -0.4890164  0.28246697  0.67111315  0.21656332  0.22283076
# Sodium.2015.Soil    -0.42009036  0.2426210  0.30017244  0.39183736 -0.57938409  0.01751195
# Phos.2015.Soil      -0.14259654 -0.5200681 -0.26794015 -0.30516923 -0.63635389  0.27603172
    
pca.gen$x
# PC1        PC2        PC3         PC4          PC5           PC6
# Anz  0.3979644 -0.2118295  2.1312075  0.36345103 -0.114021494  1.665335e-16
# BMR  1.8611415  1.6822265 -0.8649859  0.35404254 -0.435732468  7.494005e-16
# Cla -0.5374248 -2.7007382 -0.5435904 -0.38171864 -0.361534025  9.992007e-16
# Gri -3.3958880  1.4634277 -0.1764831 -0.01862357  0.003434728 -8.326673e-17
# UCR  0.4218577 -0.9193805 -0.6673506  0.81249277  0.557797887 -1.452831e-15
# Yuc  1.2523491  0.6862939  0.1212024 -1.12964414  0.350055372 -9.714451e-17

# since the top loadings of PC1 and PC2 are NEGATIVELY correlated with PC1 and PC2,
# will multiply the scores by -1 so that the PC scores are POSTIVELY correlated
# with the loadings (more intuitive to discuss)
pca.gen$x.transformed <- pca.gen$x * -1

# save General.PC1 values in a new dataframe to merge into plant data
# will actually use orig values instead of transformed, as of 12-9-21
pca.gen$x.transformed
pca.gen$x
generalPC <- as.data.frame(pca.gen$x[,1:2])
generalPC$Inoculum_source <- row.names(generalPC)
names(generalPC) <- c("General.PC1", "General.PC2", "Inoculum_source")
generalPC
# General.PC1 General.PC2 Inoculum_source
# Anz   0.3979644  -0.2118295             Anz
# BMR   1.8611415   1.6822265             BMR
# Cla  -0.5374248  -2.7007382             Cla
# Gri  -3.3958880   1.4634277             Gri
# UCR   0.4218577  -0.9193805             UCR
# Yuc   1.2523491   0.6862939             Yuc

# look at raw soil data for top loadings of General.PC1
Gen <- cbind(Gen, generalPC)
plot(Gen$CEC.2015.Soil ~ Gen$General.PC1)
plot(Gen$Calcium.2015.Soil ~ Gen$General.PC1)
plot(Gen$Magnesium.2015.Soil ~ Gen$General.PC1)
 
Gen$CEC.2015.Soil
    
#####   IMPORT NOD COLOR DATA (D2)  ##########################

# Import nodule color data
D2 <- read.csv("WSI_nodcolor_data.csv", header = TRUE)

# set data types
D2$Plant_ID <- as.factor(D2$Plant_ID)

# Average nodule counts for each color bin within each unique plant
nodcolors <- D2 %>%
  group_by(Plant_ID) %>%
  summarize(Brown.black.nods = mean(Brown.black.nods, na.rm = TRUE),
            Green.nods = mean(Green.nods, na.rm = TRUE), 
            Green.red.pink.nods = mean(Green.red.pink.nods, na.rm = TRUE),
            Red.pink.nods = mean(Red.pink.nods, na.rm = TRUE),
            White.nods = mean(White.nods, na.rm = TRUE))



#####   IMPORT PLANT DATA (D3)  ##########################

# Import nodule color data
D3 <- read.csv("WSI_plant_data.csv", header = TRUE)

# set data types
D3$Plant_ID <- as.factor(D3$Plant_ID)
D3$Experiment <- as.factor(D3$Experiment)
D3$Block <- as.factor(D3$Block)
D3$Sort_group <- as.factor(D3$Sort_group)
D3$Cluster <- as.factor(D3$Cluster)
D3$Position <- as.factor(D3$Position)
D3$Host_population <- as.factor(D3$Host_population)
D3$Host_line <- as.factor(D3$Host_line)
D3$Inoculum_source <- as.factor(D3$Inoculum_source)
D3$Inoculum_type <- as.factor(D3$Inoculum_type)
D3$Notes <- as.factor(D3$Notes)
D3$Initial_true_leaf_count <- as.numeric(D3$Initial_true_leaf_count)
D3$Initial_cotyledon_count <- as.numeric(D3$Initial_cotyledon_count)
D3$Harvester_initials <- as.factor(D3$Harvester_initials)
D3$Algae <- as.factor(D3$Algae)
D3$Flowers <- as.factor(D3$Flowers)
D3$Harvest_date <- as.Date(D3$Harvest_date, format = "%m/%d/%Y")
D3$Shoots_mg <- as.numeric(D3$Shoots_mg)
D3$Roots_mg <- as.numeric(D3$Roots_mg)
D3$Total_nodule_count <- as.numeric(D3$Total_nodule_count)
D3$Total_nodule_mass_mg <- as.numeric(D3$Total_nodule_mass_mg)


# make abbreviated Inoculum type column
D3$Inoculum_type_abbrev <- ifelse(D3$Inoculum_type=="Sterilized","Ster",
                                       ifelse(D3$Inoculum_type=="Live", "Live", ""))
D3$Inoculum_type_abbrev <- as.factor(D3$Inoculum_type_abbrev)
D3$Inoculum_type_abbrev <- ordered(D3$Inoculum_type_abbrev, levels = c("Ster", "Live"))


# make Inoculum column
D3$Inoculum <- paste(D3$Inoculum_source, D3$Inoculum_type, sep = "_")

# remove rows corresponding to missing plants
D3 <- D3[which(D3$Notes!="Plant_ID_nonexistent"),]

# merge in nitrogen PC1 values
D3 <- merge(D3,nitrogenPC, by="Inoculum_source", all.x=T)

# merge in other PC1 values
D3 <- merge(D3,generalPC, by="Inoculum_source", all.x=T)

# merge in nodule color data; change NaN to NA
D3 <- merge(D3,nodcolors, by="Plant_ID", all.x=T)

D3$Brown.black.nods[is.nan(D3$Brown.black.nods)] <- NA
D3$Green.nods[is.nan(D3$Green.nods)] <- NA
D3$Green.red.pink.nods[is.nan(D3$Green.red.pink.nods)] <- NA
D3$Red.pink.nods[is.nan(D3$Red.pink.nods)] <- NA
D3$White.nods[is.nan(D3$White.nods)] <- NA

# check that sum of averaged nodule color counts always equals total nod count
# note, there are 2 plants with nodule counts but no nodule color data: plants 1238 & 1405
D3$colorsum <- D3$Brown.black.nods + D3$Green.nods +
  D3$Green.red.pink.nods + D3$Red.pink.nods +
  D3$White.nods
D3$check <- ifelse(abs(D3$Total_nodule_count-D3$colorsum)<0.01,TRUE,FALSE)
D3[which(D3$colorsum>0),]$check %>% all()
# TRUE
D3 <- select(D3, -c("check", "colorsum"))


# add control plant shoot and root biomass as columns
# in the experimental design, each live-inoculated plant is paired with one 
    # unique sterilized-inoculated plant (paired by Block, line, and inoculum source)
# some plants died, so their paired plant will be missing and I will assign them 
    # control plant values separately

# get shoot mass and root mass of control plants for each Host_line in each Block
temp <- D3[which(D3$Inoculum_type=="Sterilized"),c("Block", "Host_line", "Inoculum_source", "Shoots_mg", "Roots_mg")]
head(temp)
# Block Host_line Inoculum_source Shoots_mg Roots_mg
# <fct> <fct>     <fct>               <dbl>    <dbl>
# 1 4     Anz11_01  Anz                 16.6      13.9
# 2 3     Anz11_01  Anz                  7.6       7.8
# 3 1     Anz11_01  Anz                 12.8       7.7
# 4 2     Anz11_01  Anz                 11.6       8  
# 5 6     Anz11_01  Anz                  9.80      6.3
# 6 7     Anz11_01  Anz                 15.8       7.2
names(temp) = c("Block", "Host_line", "Inoculum_source", "Shoots_control_mg", "Roots_control_mg")

# merge control plant shoot mass and root mass into D3 
D3 <- merge(D3,temp, by=c("Block", "Host_line", "Inoculum_source"), all.x=T)

# some control plants were missing from the experiment-- need to fill in data for these
# how many live-inoculated plants have NA for Shoots_control_mg?
temp <- D3[which(D3$Inoculum_type=="Live" & is.na(D3$Shoots_control_mg)),]
temp$Host_line %>% droplevels() %>% levels()
# "Anz10_01" "Anz13_04" "Gri01_13" are the 3 host lines missing control data
missing <- temp$Plant_ID %>% droplevels()
missing %>% length()
# 17 live-inoculated plants have NA for Shoots_control_mg
missing
# 1265 1271 1273 1276 1268 1275 1027 1147 1024 1031 1307 1305 1311 1312 1303 1308 1036


# calculate average shoot and root masses for all sterilized treatments
temp <- D3[which(D3$Inoculum_type=="Sterilized"),]
groupmeans <- temp %>%
  group_by(Host_line, Inoculum_source) %>%
  summarize(MeanShoots = mean(Shoots_mg), MeanRoots = mean(Roots_mg))
groupmeans$ID <- paste(groupmeans$Host_line, groupmeans$Inoculum_source, sep = "_")
head(groupmeans)
# Host_line    Inoculum_source MeanShoots MeanRoots ID              
# <fct>        <fct>                <dbl>     <dbl> <chr>           
# 1 A_heermannii Anz                   14.6      13.8 A_heermannii_Anz
# 2 A_heermannii BMR                   15.0      12.3 A_heermannii_BMR
# 3 A_heermannii Cla                   18.0      15.5 A_heermannii_Cla
# 4 A_heermannii Gri                   16.5      14.8 A_heermannii_Gri
# 5 A_heermannii UCR                   17.9      13.9 A_heermannii_UCR
# 6 A_heermannii Yuc                   14.0      12.0 A_heermannii_Yuc


# make a temporary column to concatenate host line and inoculum source
D3$Temp <- paste(D3$Host_line, D3$Inoculum_source, sep = "_")

# look up each Plant_ID with no control values and assign it values
for (i in missing){
  # from D3, save host line and inoc source for i as ID
  ID <- D3[which(D3$Plant_ID==i),]$Temp
  # use ID to index the temp sheet and save control value
  shoot <- groupmeans[which(groupmeans$ID==ID),]$MeanShoots
  root <- groupmeans[which(groupmeans$ID==ID),]$MeanRoots
  # assign plant i the saved control value
  D3[which(D3$Plant_ID==i),]$Shoots_control_mg = shoot
  D3[which(D3$Plant_ID==i),]$Roots_control_mg = root
}

# remove temporary column
D3 <- select(D3, -"Temp")

# make other calculated columns
D3$Initial_leaf_count <- D3$Initial_cotyledon_count + D3$Initial_true_leaf_count
D3$Totmass_mg <- D3$Shoots_mg + D3$Roots_mg
D3$Totmass_control_mg <- D3$Shoots_control_mg + D3$Roots_control_mg
D3$Relgro_totmass <- D3$Totmass_mg / D3$Totmass_control_mg
D3$Root_shoot <- D3$Roots_mg / D3$Shoots_mg
D3$Nod_size_mg <- D3$Total_nodule_mass_mg / D3$Total_nodule_count
D3$Red_nod_freq <- D3$Red.pink.nods / D3$Total_nodule_count
D3$logTotmass <- log(D3$Totmass_mg)
D3$Host_type <- ifelse(D3$Host_line=="A_heermannii", "A_heermannii",
                       ifelse(D3$Host_line=="Anz13_04", "Anz13_04",
                              ifelse(D3$Host_line=="Cla12_04", "Cla12_04",
                                     "Sympatric host")))
    D3$Host_type <- ordered(D3$Host_type,
                        levels = c("Sympatric host", "Anz13_04", "Cla12_04", "A_heermannii"))



### IMPORT CULTURE DATA (D4) ##############################
  
# Import filtrate culture data
D4 <- read.csv("WSI_culture_data.csv", header = TRUE)

# set data types
D4$Inoculum_source <- as.factor(D4$Inoculum_source)
D4$Filtrate_dilution <- as.numeric(D4$Filtrate_dilution)
D4$Plate_replicate <- as.factor(D4$Plate_replicate)
D4$Colony_count <- as.numeric(D4$Colony_count)

# estimate total CFU on each plate
D4$Total_CFU <- ifelse(D4$Colony_count>=30 & D4$Colony_count<=300,
                       D4$Colony_count*10^D4$Filtrate_dilution, NA)

# estimate CFU per mL inoculum
D4$CFU_per_mL <- 10*D4$Total_CFU

# calculate mean CFU per mL for each inoculum

meanCFU <- D4 %>%
  group_by(Inoculum_source) %>%
  summarize(Mean = mean(CFU_per_mL, na.rm = TRUE))

meanCFU
# Inoculum_source      Mean
# 1 Anz               673500 
# 2 BMR              2697500 
# 3 Cla             25766667.
# 4 Gri             27333333.
# 5 UCR              5882000 
# 6 Yuc             10166667.


################# SUBSET THE PLANT DATA ###################

# check on plants that have NA for nod color data; were they all non-nodulated?
D3[which(is.na(D3$Brown.black.nods)),]$Total_nodule_count %>% unique # 0 52 27
D3[which(D3$Total_nodule_count==52 & is.na(D3$Brown.black.nods)),]$Plant_ID # 1405
D3[which(D3$Total_nodule_count==27 & is.na(D3$Brown.black.nods)),]$Plant_ID # 1238

# flag early harvest plants
D3$Early.harvest <- ifelse(D3$Notes=="Shoot_portion_collected_early" |
                           D3$Notes=="Entire_shoot_collected_early",
                           "yes", "no")
D3[which(D3$Early.harvest=="yes"),] %>% nrow
# 27 plants, correct

all.noAnz.live <- D3[which(D3$Inoculum_source!="Anz" &
                           D3$Inoculum_type!="Sterilized"),] %>% droplevels()
    all.noAnz.live %>% nrow()
    # 244
    # order inoculum source levels by Nitrogen.PC1
    all.noAnz.live$Inoculum_source <- ordered(all.noAnz.live$Inoculum_source,
                                          levels = c("BMR", "Yuc", "Gri", "UCR", "Cla"))
    


symp <- D3[which(D3$Experiment == "Sympatric"),] %>% droplevels()
    symp %>% nrow()
    # 219
    
symp.noAnz <- symp[which(symp$Inoculum_source!="Anz"),] %>% droplevels()
    # change order of inoculum type levels
    symp.noAnz$Inoculum_type <- ordered(symp.noAnz$Inoculum_type,
                                        levels = c("Sterilized", "Live"))
    # order inoculum source levels by Nitrogen.PC1
    symp.noAnz$Inoculum_source <- ordered(symp.noAnz$Inoculum_source,
                                          levels = c("BMR", "Yuc", "Gri", "UCR", "Cla"))
    symp.noAnz %>% nrow()
    # 187
    
symp.noAnz.live <- symp.noAnz[which(symp.noAnz$Inoculum_type=="Live"),] %>% droplevels()
    symp.noAnz.live %>% nrow()
    # 94
    
symp.noAnz.ster <- symp.noAnz[which(symp.noAnz$Inoculum_type=="Sterilized"),] %>% droplevels()
    symp.noAnz.ster %>% nrow()
    # 93

univ <- D3[which(D3$Experiment == "Universal"),] %>% droplevels()
    univ %>% nrow()
    # 348

univ.noAnz <- univ[which(univ$Inoculum_source!="Anz"),] %>% droplevels()
    # change order of inoculum type levels
    univ.noAnz$Inoculum_type <- ordered(univ.noAnz$Inoculum_type,
                                        levels = c("Sterilized", "Live"))
    # order inoculum source levels by Nitrogen.PC1
    univ.noAnz$Inoculum_source <- ordered(univ.noAnz$Inoculum_source,
                                               levels = c("BMR", "Yuc", "Gri", "UCR", "Cla"))
    univ.noAnz %>% nrow()
    # 290
    
univ.noAnz.live <- univ.noAnz[which(univ.noAnz$Inoculum_type=="Live"),] %>% droplevels()
    univ.noAnz.live %>% nrow()
    # 150
    
all.noAnz.live$Host_type %>% levels
all.noAnz.live$Host_type_label <- ifelse(all.noAnz.live$Host_type=="Sympatric host", "Sympatric host",
                                         ifelse(all.noAnz.live$Host_type=="Anz13_04", "Anz13.04 host",
                                                ifelse(all.noAnz.live$Host_type=="Cla12_04", "Cla12.04 host",
                                                       ifelse(all.noAnz.live$Host_type=="A_heermannii", "A. heermannii host", ""))))
all.noAnz.live$Host_type_label <- ordered(all.noAnz.live$Host_type_label, levels = c("Sympatric host", "Anz13.04 host", "Cla12.04 host", "A. heermannii host"))

    

 
#### FIGURES ###########################################################################  

# *** Fig 1: PCA of soil traits ######

# A. Nitrogen PCA
p1 <- factoextra::fviz_pca(pca.nit, repel = TRUE)
ggplot(p1 +
    theme(element_text(size = 4)))
grid.arrange(p1)

# B. Non-nitrogen PCA
p2 <- factoextra::fviz_pca(pca.gen, repel = TRUE)
grid.arrange(p2)



# *** Fig 2. Total plant mass reaction norms ####

# ****** A. Sympatric hosts

usedata <- complete_fun(symp.noAnz, "Totmass_mg")
usedata %>% nrow()
# 187

smallgroups <- usedata %>%
  group_by(Inoculum_type_abbrev, Inoculum_source) %>%
  summarize(Mean = mean(Totmass_mg), CIup = CI(Totmass_mg), CIlow = CI2(Totmass_mg))



p1 <- ggplot(smallgroups, aes(x = Inoculum_type_abbrev, y = Mean, color = Inoculum_source)) +
  ylab("Total plant mass (mg)") +
  xlab("") +
  ylim(0,225) +
  geom_line(aes(group = Inoculum_source), 
            size = 1.1,
            alpha = 0.8, 
            show.legend = TRUE) +
  scale_color_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                     values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Inoculum_source, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none")
p1 <- arrangeGrob(p1, top = textGrob("A. Sympatric host", x = unit(0, "npc"),
                                     y = unit(0.7, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=14)))
grid.arrange(p1)


# ****** B. Anz13.04 hosts

usedata <- complete_fun(univ.noAnz, "Totmass_mg")
usedata <- usedata[which(usedata$Host_line=="Anz13_04"),]
usedata %>% nrow()
# 90

smallgroups <- usedata %>%
  group_by(Inoculum_type_abbrev, Inoculum_source, Host_line) %>%
  summarize(Mean = mean(Totmass_mg), CIup = CI(Totmass_mg), CIlow = CI2(Totmass_mg))


p2 <- ggplot(smallgroups, aes(x = Inoculum_type_abbrev, y = Mean, color = Inoculum_source)) +
  ylab("Total plant mass (mg)") +
  xlab("") +
  ylim(0,225) +
  geom_line(aes(group = Inoculum_source), 
            size = 1.1,
            alpha = 0.8, 
            show.legend = TRUE) +
  scale_color_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                     values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Inoculum_source, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none")
p2 <- arrangeGrob(p2, top = textGrob("B. Anz13.04 host", x = unit(0, "npc"),
                                     y = unit(0.7, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=14)))
grid.arrange(p2)


# ****** C. Cla12.04 hosts

usedata <- complete_fun(univ.noAnz, "Totmass_mg")
usedata <- usedata[which(usedata$Host_line=="Cla12_04"),]
usedata %>% nrow()
# 100

smallgroups <- usedata %>%
  group_by(Inoculum_type_abbrev, Inoculum_source, Host_line) %>%
  summarize(Mean = mean(Totmass_mg), CIup = CI(Totmass_mg), CIlow = CI2(Totmass_mg))


p3 <- ggplot(smallgroups, aes(x = Inoculum_type_abbrev, y = Mean, color = Inoculum_source)) +
  ylab("Total plant mass (mg)") +
  xlab("") +
  ylim(0,225) +
  geom_line(aes(group = Inoculum_source), 
            size = 1.1,
            alpha = 0.8, 
            show.legend = TRUE) +
  scale_color_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                     values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Inoculum_source, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none")
p3 <- arrangeGrob(p3, top = textGrob("C. Cla12.04 host", x = unit(0, "npc"),
                                     y = unit(0.7, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=14)))
grid.arrange(p3)


# ****** D. A. heermannii hosts

usedata <- complete_fun(univ.noAnz, "Totmass_mg")
usedata <- usedata[which(usedata$Host_line=="A_heermannii"),]
usedata %>% nrow()
# 100

smallgroups <- usedata %>%
  group_by(Inoculum_type_abbrev, Inoculum_source, Host_line) %>%
  summarize(Mean = mean(Totmass_mg), CIup = CI(Totmass_mg), CIlow = CI2(Totmass_mg))

p4 <- ggplot(smallgroups, aes(x = Inoculum_type_abbrev, y = Mean, color = Inoculum_source)) +
  ylab("Total plant mass (mg)") +
  xlab("") +
  ylim(0,225) +
  geom_line(aes(group = Inoculum_source), 
            size = 1.1,
            alpha = 0.8, 
            show.legend = TRUE) +
  scale_color_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                     values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Inoculum_source, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none")
p4 <- arrangeGrob(p4, top = textGrob("D. A. heermannii host", x = unit(0, "npc"),
                                     y = unit(0.7, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=14)))
grid.arrange(p4)

grid.arrange(p1, p2, p3, p4, nrow = 2)






# *** Fig 3. Root:shoot reaction norms  ####


# ****** A. Sympatric hosts

usedata <- complete_fun(symp.noAnz, "Root_shoot")
usedata %>% nrow()
# 187


smallgroups <- usedata %>%
  group_by(Inoculum_type_abbrev, Inoculum_source) %>%
  summarize(Mean = mean(Root_shoot), CIup = CI(Root_shoot), CIlow = CI2(Root_shoot))


p1 <- ggplot(smallgroups, aes(x = Inoculum_type_abbrev, y = Mean, color = Inoculum_source)) +
  ylab("Root:shoot ratio") +
  xlab("") +
  ylim(0,1) +
  geom_line(aes(group = Inoculum_source), 
            size = 1.1,
            alpha = 0.8, 
            show.legend = TRUE) +
  scale_color_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                     values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Inoculum_source, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none")
p1 <- arrangeGrob(p1, top = textGrob("A. Sympatric host", x = unit(0, "npc"),
                                     y = unit(0.7, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=14)))
grid.arrange(p1)


# ****** B. Anz13.04 hosts

usedata <- complete_fun(univ.noAnz, "Root_shoot")
usedata <- usedata[which(usedata$Host_line=="Anz13_04"),]
usedata %>% nrow()
# 90

smallgroups <- usedata %>%
  group_by(Inoculum_type_abbrev, Inoculum_source, Host_line) %>%
  summarize(Mean = mean(Root_shoot), CIup = CI(Root_shoot), CIlow = CI2(Root_shoot))


p2 <- ggplot(smallgroups, aes(x = Inoculum_type_abbrev, y = Mean, color = Inoculum_source)) +
  ylab("Root:shoot ratio") +
  xlab("") +
  ylim(0,1) +
  geom_line(aes(group = Inoculum_source), 
            size = 1.1,
            alpha = 0.8, 
            show.legend = TRUE) +
  scale_color_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                     values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Inoculum_source, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none")
p2 <- arrangeGrob(p2, top = textGrob("B. Anz13.04 host", x = unit(0, "npc"),
                                     y = unit(0.7, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=14)))
grid.arrange(p2)


# ****** C. Cla12.04 hosts 

usedata <- complete_fun(univ.noAnz, "Root_shoot")
usedata <- usedata[which(usedata$Host_line=="Cla12_04"),]
usedata %>% nrow()
# 100

smallgroups <- usedata %>%
  group_by(Inoculum_type_abbrev, Inoculum_source, Host_line) %>%
  summarize(Mean = mean(Root_shoot), CIup = CI(Root_shoot), CIlow = CI2(Root_shoot))


p3 <- ggplot(smallgroups, aes(x = Inoculum_type_abbrev, y = Mean, color = Inoculum_source)) +
  ylab("Root:shoot ratio") +
  xlab("") +
  ylim(0,1) +
  geom_line(aes(group = Inoculum_source), 
            size = 1.1,
            alpha = 0.8, 
            show.legend = TRUE) +
  scale_color_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                     values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Inoculum_source, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none")
p3 <- arrangeGrob(p3, top = textGrob("C. Cla12.04 host", x = unit(0, "npc"),
                                     y = unit(0.7, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=14)))
grid.arrange(p3)


# ****** D. A. heermannii hosts

usedata <- complete_fun(univ.noAnz, "Root_shoot")
usedata <- usedata[which(usedata$Host_line=="A_heermannii"),]
usedata %>% nrow()
# 100

smallgroups <- usedata %>%
  group_by(Inoculum_type_abbrev, Inoculum_source, Host_line) %>%
  summarize(Mean = mean(Root_shoot), CIup = CI(Root_shoot), CIlow = CI2(Root_shoot))

p4 <- ggplot(smallgroups, aes(x = Inoculum_type_abbrev, y = Mean, color = Inoculum_source)) +
  ylab("Root:shoot ratio") +
  xlab("") +
  ylim(0,1) +
  geom_line(aes(group = Inoculum_source), 
            size = 1.1,
            alpha = 0.8, 
            show.legend = TRUE) +
  scale_color_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                     values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Inoculum_source, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none")
p4 <- arrangeGrob(p4, top = textGrob("D. A. heermannii host", x = unit(0, "npc"),
                                     y = unit(0.7, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=14)))
grid.arrange(p4)

grid.arrange(p1, p2, p3, p4, nrow = 2)



# *** Fig 4. Nodule count ####

usedata <- complete_fun(all.noAnz.live, "Total_nodule_count")
usedata %>% nrow()
# 244

smallgroups <- usedata %>%
  group_by(Inoculum_source, Host_type_label) %>%
  summarize(Mean = mean(Total_nodule_count), CIup = CI(Total_nodule_count), CIlow = CI2(Total_nodule_count))

p1 <- ggplot(smallgroups, aes(x = Inoculum_source, y = Mean, fill = Inoculum_source)) +
  ylab("Total nodule count") +
  xlab("") +
  geom_col() +
  scale_fill_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                    values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Host_type_label, scales = "free_x", nrow = 2) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none")

grid.arrange(p1)



# *** Fig S1. Map of field sites ####

# extract GPS coordinates from D1
coord <- D1[which(D1$Notes=="decimal degrees"),]

# pivot data so Lat and Long are in diff columns
coord <- pivot_wider(coord, names_from = "Metric", values_from = "Value")

# rearrange columns
coord <- select(coord, c("Longitude", "Latitude", "Population"))

# transform for usmap projection
coord.t <- usmap_transform(coord)

coord.t$Population <- ordered(coord.t$Population, levels = c("BMR", "Anz", "Yuc", "Gri", "UCR", "Cla"))

plot_usmap(include = "CA", fill = "grey", alpha = 0.25) +
  ggrepel::geom_label_repel(data = coord.t,
                            aes(x = Longitude.1, y = Latitude.1, label = Population),
                            size = 4, alpha = 0.8,
                            segment.color = "black", 
                            segment.size = 0.7,
                            min.segment.length = 0,
                            box.padding = 0.9,
                            nudge_x = 0.5,
                            nudge_y = 0.5,
                            seed = 1003) +
  geom_point(data = coord.t,
             aes(x = Longitude.1, y = Latitude.1, size = 1, color = Population)) +
  scale_color_manual(breaks = c("BMR", "Anz", "Yuc", "Gri", "UCR", "Cla"), 
                     values=c("royalblue", "cornflowerblue", "lightblue", "pink", "red", "red4")) +
  theme(legend.position = "none")



# *** Fig S4. Plant relative growth ####

usedata <- complete_fun(all.noAnz.live, "Relgro_totmass")
usedata %>% nrow()
# 244

smallgroups <- usedata %>%
  group_by(Inoculum_source, Host_type_label) %>%
  summarize(Mean = mean(Relgro_totmass), CIup = CI(Relgro_totmass), CIlow = CI2(Relgro_totmass))

p1 <- ggplot(smallgroups, aes(x = Inoculum_source, y = Mean, fill = Inoculum_source)) +
  ylab("Plant relative growth") +
  xlab("") +
  geom_col() +
  scale_fill_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                     values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Host_type_label, scales = "free_x", nrow = 2) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none")

grid.arrange(p1)




# *** Fig S5. Red nodule frequency ####

usedata <- complete_fun(all.noAnz.live, "Red_nod_freq")
usedata %>% nrow()
# 242

smallgroups <- usedata %>%
  group_by(Inoculum_source, Host_type_label) %>%
  summarize(Mean = mean(Red_nod_freq), CIup = CI(Red_nod_freq), CIlow = CI2(Red_nod_freq))

p1 <- ggplot(smallgroups, aes(x = Inoculum_source, y = Mean, fill = Inoculum_source)) +
  ylab("Red nodule frequency") +
  xlab("") +
  geom_col() +
  scale_fill_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                    values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Host_type_label, scales = "free_x", nrow = 2) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none")

grid.arrange(p1)





# *** Fig S6. Total nodule mass ####

usedata <- complete_fun(all.noAnz.live, "Total_nodule_mass_mg")
usedata %>% nrow()
# 244

smallgroups <- usedata %>%
  group_by(Inoculum_source, Host_type_label) %>%
  summarize(Mean = mean(Total_nodule_mass_mg), CIup = CI(Total_nodule_mass_mg), CIlow = CI2(Total_nodule_mass_mg))

p1 <- ggplot(smallgroups, aes(x = Inoculum_source, y = Mean, fill = Inoculum_source)) +
  ylab("Total nodule mass (mg)") +
  xlab("") +
  geom_col() +
  scale_fill_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                    values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Host_type_label, scales = "free_x", nrow = 2) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none")

grid.arrange(p1)





# *** Fig S7. Mean nodule size ####

usedata <- complete_fun(all.noAnz.live, "Nod_size_mg")
usedata %>% nrow()
# 244

smallgroups <- usedata %>%
  group_by(Inoculum_source, Host_type_label) %>%
  summarize(Mean = mean(Nod_size_mg), CIup = CI(Nod_size_mg), CIlow = CI2(Nod_size_mg))

p5 <- ggplot(smallgroups, aes(x = Inoculum_source, y = Mean, fill = Inoculum_source)) +
  ylab("Mean nodule size (mg)") +
  xlab("") +
  geom_col() +
  scale_fill_manual(breaks = c("BMR", "Yuc", "Gri", "UCR", "Cla"), 
                    values=c("royalblue", "lightblue", "pink", "red", "red4")) +
  geom_errorbar(data = smallgroups, aes(ymin=CIlow, ymax=CIup, width = 0.4)) +
  facet_wrap(~Host_type_label, scales = "free_x", nrow = 2) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none")

grid.arrange(p5)



#### ANALYSES ####

#### SYMPATRIC HOST / TOTAL PLANT MASS ####

#### *** Data exploration ####

par(mfrow = c(1,1))
dotchart(symp.noAnz$Totmass_mg, groups = symp.noAnz$Inoculum_type)
dotchart(symp.noAnz$Totmass_mg, groups = symp.noAnz$Host_line)
dotchart(symp.noAnz$Totmass_mg, groups = symp.noAnz$Nitrogen.PC1)
dotchart(symp.noAnz$Totmass_mg, groups = symp.noAnz$General.PC1)
# don't see any outliers
# sterilized inculum has much lower variance than live inoculum
# one of the Gri host lines seems to have higher variance than the other lines
# variance for PC1 variables looks pretty good

data = symp.noAnz$Totmass_mg
summary(data)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.80   14.85   31.50   63.67  104.90  257.00 

# calculate mean and SE
mean(data, na.rm = TRUE)
# 63.66524
sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
# 4.348615

## *** usedata ####

usedata <- complete_fun(symp.noAnz, "logTotmass")
usedata %>% nrow()
# 187

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 176

####  *** Nitrogen.PC1 model ####

M1 <- lmer(logTotmass ~ Nitrogen.PC1 + Inoculum_type + Nitrogen.PC1:Inoculum_type +
          (1|Host_line) + (1|Block), data = usedata, REML = FALSE)
    
# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept) 0.06610  0.2571  
# Block     (Intercept) 0.01366  0.1169  
# Residual              0.27311  0.5226  
# Number of obs: 187, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#                         Estimate Std. Error t value
# (Intercept)                   3.61094    0.09911  36.435
# Nitrogen.PC1                  0.06923    0.04447   1.557
# Inoculum_type.L               1.27039    0.05489  23.143
# Nitrogen.PC1:Inoculum_type.L  0.04503    0.02603   1.730


# test effect of interaction
drop1(M1, test = "Chisq")
#                            Df    AIC   LRT Pr(Chi)  
# <none>                        324.93                
# Nitrogen.PC1:Inoculum_type  1 325.90 2.967 0.08498 .
     
# make reduced model
M.int <- update(M1, .~. - Nitrogen.PC1:Inoculum_type)

# test main effects
drop1(M.int, test = "Chisq")
#               Df    AIC     LRT Pr(Chi)    
# <none>           325.90                    
# Nitrogen.PC1   1 326.06   2.161  0.1416    
# Inoculum_type  1 571.96 248.058  <2e-16 ***

# Totmass differs between live and sterilized treatments, but not between nitrogen regimes



# effect size of Inoculum_type

live <- usedata[which(usedata$Inoculum_type=="Live"),"Totmass_mg"] %>% mean()
live
# 109.9404 mg totmass with live inoculum

ster <- usedata[which(usedata$Inoculum_type=="Sterilized"),"Totmass_mg"] %>% mean()
ster
# 16.89247 mg totmass with sterilized inoculum

(live-ster)/ster
# 5.508249
# plants with live inoculum have 551% greater total biomass than plants with sterilized inoculum

live/ster
# 6.508249
# plants with live inoculum have 6.5x the biomass of plants with sterilized inoculum



# *** NonNitrogen.PC1 model ####    
    
M1 <- lmer(logTotmass ~ General.PC1 + Inoculum_type + General.PC1:Inoculum_type +
             (1|Host_line) + (1|Block), data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept) 0.02777  0.1666  
# Block     (Intercept) 0.01152  0.1073  
# Residual              0.27781  0.5271  
# Number of obs: 187, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                  3.63914    0.07411  49.102
# General.PC1                  0.13990    0.03758   3.722
# Inoculum_type.L              1.28322    0.05475  23.438
# General.PC1:Inoculum_type.L -0.01399    0.03251  -0.430


# test effect of interaction
drop1(M1, test = "Chisq")
# Df    AIC     LRT Pr(Chi)
# <none>                       320.97                
# General.PC1:Inoculum_type  1 319.16 0.18503  0.6671

# make reduced model
M.int <- update(M1, .~. - General.PC1:Inoculum_type)

# test main effects
drop1(M.int, test = "Chisq")
# Df    AIC     LRT   Pr(Chi)    
# <none>           319.16                      
# General.PC1    1 326.06   8.902  0.002848 ** 
# Inoculum_type  1 565.49 248.333 < 2.2e-16 ***

# interesting-- soils DO differ in the growth-promoting effects of their microbes
# but just not predicted by N, but other soil traits
# General.PC1 and TotMass are positively related

# overall conclusion from analysis of total plant mass:
    # soil physicochemical differences do explain some variation in the growth-promoting
    # ability of soil microbes for sympatric plants, but nitrogen regimes do not drive
    # these differences
    

# *** NonNitrogen.PC2 model ####    

M1 <- lmer(logTotmass ~ General.PC2 + Inoculum_type + General.PC2:Inoculum_type +
             (1|Host_line) + (1|Block), data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

# look at parameter estimates
summary(M1)

# test effect of interaction
drop1(M1, test = "Chisq")
# npar    AIC    LRT Pr(Chi)  
# <none>                         324.49                 
# General.PC2:Inoculum_type    1 328.05 5.5601 0.01837 *

# make reduced model
M.int <- update(M1, .~. - General.PC2:Inoculum_type)

# test main effects
drop1(M.int, test = "Chisq")
# npar    AIC     LRT Pr(Chi)    
# <none>             328.05                    
# General.PC2      1 326.06   0.007  0.9344    
# Inoculum_type    1 574.14 248.087  <2e-16 ***


# ****** exploratory figure ####

p1 <- ggplot(usedata, aes(x = General.PC2, y = logTotmass)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("log(Total plant mass)") +
  xlab("General.PC2") +
  facet_wrap(~Inoculum_type, scales = "free_x") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

print(p1)




#### SYMPATRIC HOST / ROOT:SHOOT RATIO  ####

#### *** Data exploration ####

par(mfrow = c(1,1))
dotchart(symp.noAnz$Root_shoot, groups = symp.noAnz$Inoculum_type)
dotchart(symp.noAnz$Root_shoot, groups = symp.noAnz$Host_line)
dotchart(symp.noAnz$Root_shoot, groups = symp.noAnz$Nitrogen.PC1)
dotchart(symp.noAnz$Root_shoot, groups = symp.noAnz$General.PC1)
# possible high outlier in sterilized treatment, PlantID = 1074
# sterilized inculum has higher variance than live inoculum

data = symp.noAnz$Root_shoot
summary(data)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1063  0.2823  0.4286  0.5288  0.7213  1.8500 

# calculate mean and SE
mean(data, na.rm = TRUE)
# 0.5287716
sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
# 0.02226785
    

## *** usedata ####

usedata <- complete_fun(symp.noAnz, "Root_shoot")
usedata %>% nrow()
# 187

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 176

usedata.o <- usedata[which(usedata$Plant_ID!=1074),]
    
    
####  *** Nitrogen.PC1 model  ####
    
M1 <- lmer(Root_shoot ~ Nitrogen.PC1 + Inoculum_type + Nitrogen.PC1:Inoculum_type +
            (1|Host_line) + (1|Block), data = usedata, REML = FALSE)
    
# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance  Std.Dev.
# Host_line (Intercept) 0.0079906 0.08939 
# Block     (Intercept) 0.0004994 0.02235 
# Residual              0.0301434 0.17362 
# Number of obs: 187, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                   0.521824   0.032483  16.064
# Nitrogen.PC1                 -0.007105   0.015346  -0.463
# Inoculum_type.L              -0.328795   0.018236 -18.030
# Nitrogen.PC1:Inoculum_type.L -0.003643   0.008647  -0.421

# test effect of interaction
drop1(M1, test = "Chisq")
# Df     AIC     LRT Pr(Chi)
# <none>                        -90.066                
# Nitrogen.PC1:Inoculum_type  1 -91.888 0.17745  0.6736

# make reduced model
M.int <- update(M1, .~. - Nitrogen.PC1:Inoculum_type)

# test main effects
drop1(M.int, test = "Chisq")
# Df     AIC     LRT Pr(Chi)    
# <none>           -91.888                    
# Nitrogen.PC1   1 -93.676   0.212  0.6453    
# Inoculum_type  1  93.278 187.167  <2e-16 ***


# the root:shoot ratio of plants shifts with the presence or absence of sympatric microbes
# plants with microbes have lower root:shoot ratio (less need to forage for nutrients)



# effect size of Inoculum_type

live <- usedata[which(usedata$Inoculum_type=="Live"),"Root_shoot"] %>% mean()
live
# 0.2957113 with live inoculum

ster <- usedata[which(usedata$Inoculum_type=="Sterilized"),"Root_shoot"] %>% mean()
ster
# 0.764338 with sterilized inoculum

(ster-live)/live
# 1.584744
# plants with sterilized inoculum have 158% higher root:shoot ratio than plants with live inoculum

    
# was potential outlier influential?    
M1o <- lmer(Root_shoot ~ Nitrogen.PC1 + Inoculum_type + Nitrogen.PC1:Inoculum_type +
          (1|Host_line) + (1|Block), data = usedata.o, REML = FALSE)
    

# test effect of interaction
drop1(M1o, test = "Chisq")
# Df     AIC    LRT Pr(Chi)
# <none>                        -128.26               
# Nitrogen.PC1:Inoculum_type  1 -128.98 1.2793   0.258

# make reduced model
Mo.int <- update(M1o, .~. - Nitrogen.PC1:Inoculum_type)

# test main effects
drop1(Mo.int, test = "Chisq")
# Df      AIC     LRT Pr(Chi)    
# <none>           -128.977                    
# Nitrogen.PC1   1 -130.922   0.055  0.8144    
# Inoculum_type  1   74.153 205.130  <2e-16 ***
   
# excluding outlier did not change significance of terms; will use M1
    
    
####  *** NonNitrogen.PC1 model  ####

M1 <- lmer(Root_shoot ~ General.PC1 + Inoculum_type + General.PC1:Inoculum_type +
            (1|Host_line) + (1|Block), data = usedata, REML = FALSE)
    
# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# ok residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance  Std.Dev.
# Host_line (Intercept) 0.0021131 0.04597 
# Block     (Intercept) 0.0004785 0.02187 
# Residual              0.0295509 0.17190 
# Number of obs: 187, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                  0.52074    0.02059  25.286
# General.PC1                 -0.04473    0.01114  -4.015
# Inoculum_type.L             -0.32721    0.01786 -18.326
# General.PC1:Inoculum_type.L  0.01796    0.01060   1.694


# test effect of interaction
drop1(M1, test = "Chisq")
# Df     AIC    LRT Pr(Chi)  
# <none>                       -102.97                 
# General.PC1:Inoculum_type  1 -102.13 2.8427 0.09179 .

# make reduced model
M.int <- update(M1, .~. - General.PC1:Inoculum_type)

# test main effects
drop1(M.int, test = "Chisq")
# Df      AIC     LRT   Pr(Chi)    
# <none>           -102.128                      
# General.PC1    1  -93.676  10.452  0.001225 ** 
# Inoculum_type  1   84.039 188.168 < 2.2e-16 ***


# Similar to TotMass, Root:Shoot ratio responds to PC1.G

# overall, soils with high levels of Ca, CEC, and Mg (e.g., Gri) tend to have plants with high total
# shoot mass, and low root:shoot ratios

# but we don't find any interaction between soil traits and microbial effects
# microbes (live treatment) always increase shoot mass and decrease root:shoot ratio,
# irrespective of soil traits

    
# was potential outlier influential?    

M1o <- lmer(Root_shoot ~ General.PC1 + Inoculum_type + General.PC1:Inoculum_type +
              (1|Host_line) + (1|Block), data = usedata.o, REML = FALSE)


# test effect of interaction
drop1(M1o, test = "Chisq")
# Df     AIC    LRT Pr(Chi)
# <none>                       -138.49               
# General.PC1:Inoculum_type  1 -138.57 1.9162  0.1663

# make reduced model
Mo.int <- update(M1o, .~. - General.PC1:Inoculum_type)

# test main effects
drop1(Mo.int, test = "Chisq")
# Df     AIC     LRT   Pr(Chi)    
# <none>           -138.57                      
# General.PC1    1 -130.92   9.647  0.001896 ** 
# Inoculum_type  1   65.72 206.290 < 2.2e-16 ***

# excluding outlier did not change significance of terms; will use M1



# *** NonNitrogen.PC2 model ####    

M1 <- lmer(Root_shoot ~ General.PC2 + Inoculum_type + General.PC2:Inoculum_type +
             (1|Host_line) + (1|Block), data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

# look at parameter estimates
summary(M1)

# test effect of interaction
drop1(M1, test = "Chisq")
# npar     AIC    LRT Pr(Chi)
# <none>                         -92.075               
# General.PC2:Inoculum_type    1 -92.353 1.7221  0.1894

# make reduced model
M.int <- update(M1, .~. - General.PC2:Inoculum_type)

# test main effects
drop1(M.int, test = "Chisq")
# npar     AIC     LRT Pr(Chi)    
# <none>             -92.353                    
# General.PC2      1 -93.676   0.676  0.4109    
# Inoculum_type    1  92.972 187.325  <2e-16 ***



  
#### SYMPATRIC HOST / PLANT RELATIVE GROWTH ####

#### *** Data exploration ####

par(mfrow = c(1,1))
dotchart(symp.noAnz.live$Relgro_totmass, groups = symp.noAnz.live$Host_line)
dotchart(symp.noAnz.live$Relgro_totmass, groups = symp.noAnz.live$Nitrogen.PC1)
dotchart(symp.noAnz.live$Relgro_totmass, groups = symp.noAnz.live$General.PC1)
# no outliers

data = symp.noAnz.live$Relgro_totmass
summary(data)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2941  4.4878  7.2390  7.6090  9.9356 22.4932 

# calculate mean and SE
mean(data, na.rm = TRUE)
# 7.609027
sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
# 0.4577726

## *** usedata ####

usedata <- complete_fun(symp.noAnz.live, "Relgro_totmass")
usedata %>% nrow()
# 94

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 84

####  *** Nitrogen.PC1 model  ####

M1 <- lmer(Relgro_totmass ~ Nitrogen.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)


# check singularity
model = M1
isSingular(model, tol = 1e-05)
# TRUE

summary(M1)
# block has variance of zero

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept)  4.808   2.193   
# Block     (Intercept)  0.000   0.000   
# Residual              14.290   3.780   
# Number of obs: 94, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)    7.3599     0.8134   9.049
# Nitrogen.PC1   0.3311     0.3932   0.842

# test main effects
drop1(M1, test = "Chisq")
# npar    AIC     LRT Pr(Chi)
# <none>            540.88                
# Nitrogen.PC1    1 539.56 0.68252  0.4087



####  *** NonNitrogen.PC1 model  ####

M1 <- lmer(Relgro_totmass ~ General.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# TRUE

summary(M1)
# block variance estimated as zero

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept)  4.984   2.232   
# Block     (Intercept)  0.000   0.000   
# Residual              14.252   3.775   
# Number of obs: 94, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   7.4795     0.8108   9.225
# General.PC1   0.3671     0.4534   0.810

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC     LRT Pr(Chi)
# <none>         540.91                
# General.PC1  1 539.56 0.64485   0.422




# *** NonNitrogen.PC2 model ####    

M1 <- lmer(Relgro_totmass ~ General.PC2 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# TRUE

summary(M1)
# block variance estimated very low

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

# look at parameter estimates
summary(M1)

# test main effect
drop1(M1, test = "Chisq")
# npar    AIC    LRT Pr(Chi)
# <none>           539.67               
# General.PC2    1 539.56 1.8884  0.1694




#### SYMPATRIC HOST / RED NODULE FREQUENCY ####

#### *** Data exploration ####

par(mfrow = c(1,1))
dotchart(symp.noAnz.live$Red_nod_freq, groups = symp.noAnz.live$Host_line)
dotchart(symp.noAnz.live$Red_nod_freq, groups = symp.noAnz.live$Nitrogen.PC1)
dotchart(symp.noAnz.live$Red_nod_freq, groups = symp.noAnz.live$General.PC1)
# no outliers

data = symp.noAnz.live$Red_nod_freq
summary(data)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  0.3684  0.5490  0.5440  0.7381  0.9928       1 

# calculate mean and SE
mean(data, na.rm = TRUE)
# 0.5439951
sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
# 0.02614036

## *** usedata ####

usedata <- complete_fun(symp.noAnz.live, "Red_nod_freq")
usedata %>% nrow()
# 93

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 83


####  *** Nitrogen.PC1 model  ####

M1 <- lmer(Red_nod_freq ~ Nitrogen.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)


# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept) 0.007009 0.08372 
# Block     (Intercept) 0.012517 0.11188 
# Residual              0.041695 0.20419 
# Number of obs: 93, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   0.52687    0.04959  10.625
# Nitrogen.PC1  0.02762    0.01668   1.656

# test main effects
drop1(M1, test = "Chisq")
# Df     AIC    LRT Pr(Chi)
# <none>          0.44085               
# Nitrogen.PC1  1 0.83938 2.3985  0.1214

# no effect of Nitrogen on Red nodule frequency


####  *** NonNitrogen.PC1 model  ####

M1 <- lmer(Red_nod_freq ~ General.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept) 0.01025  0.1013  
# Block     (Intercept) 0.01276  0.1129  
# Residual              0.04157  0.2039  
# Number of obs: 93, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  0.536146   0.052634  10.186
# General.PC1 -0.009263   0.021794  -0.425

# test main effects
drop1(M1, test = "Chisq")
# Df     AIC     LRT Pr(Chi)
# <none>         2.65923                
# General.PC1  1 0.83938 0.18015  0.6712

# no effect of General soil traits on Red nodule frequency


# *** NonNitrogen.PC2 model ####    

M1 <- lmer(Red_nod_freq ~ General.PC2 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

# look at parameter estimates
summary(M1)

# test main effect
drop1(M1, test = "Chisq")
# npar      AIC    LRT Pr(Chi)  
# <none>           -3.11595                 
# General.PC2    1  0.83938 5.9553 0.01467 *



# ****** exploratory figure ####

p1 <- ggplot(usedata, aes(x = General.PC2, y = Red_nod_freq)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Red nodule frequency)") +
  xlab("General.PC2") +
  #facet_wrap(~Inoculum_type, scales = "free_x") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

print(p1)

#### SYMPATRIC HOST / TOTAL NODULE MASS ####

#### *** Data exploration ####

par(mfrow = c(1,1))
dotchart(symp.noAnz.live$Total_nodule_mass_mg, groups = symp.noAnz.live$Host_line)
dotchart(symp.noAnz.live$Total_nodule_mass_mg, groups = symp.noAnz.live$Nitrogen.PC1)
dotchart(symp.noAnz.live$Total_nodule_mass_mg, groups = symp.noAnz.live$General.PC1)
# no outliers

data = symp.noAnz.live$Total_nodule_mass_mg
summary(data)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.700   8.125  11.550  11.319  14.775  23.400 

# calculate mean and SE
mean(data, na.rm = TRUE)
# 11.31915
sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
# 0.5300084

## *** usedata ####

usedata <- complete_fun(symp.noAnz.live, "Total_nodule_mass_mg")
usedata %>% nrow()
# 94

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 84


####  *** Nitrogen.PC1 model  ####

M1 <- lmer(Total_nodule_mass_mg ~ Nitrogen.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept)  0.8181  0.9045  
# Block     (Intercept)  0.9588  0.9792  
# Residual              22.9310  4.7886  
# Number of obs: 94, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   11.0751     0.6595   16.79
# Nitrogen.PC1   0.5833     0.2778    2.10

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)  
# <none>          577.33                 
# Nitrogen.PC1  1 578.92 3.5942 0.05798 .

# marginal effect of nitrogen on total nodule mass

# effect size of Nitrogen.PC1
# for every 1 unit increase in Nitrogen.PC1, total nodule mass increases by 0.5833 mg


####  *** NonNitrogen.PC1 model  ####

M1 <- lmer(Total_nodule_mass_mg ~ General.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept)  2.33    1.526   
# Block     (Intercept)  1.12    1.058   
# Residual              22.77    4.772   
# Number of obs: 94, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept) 11.25912    0.77102  14.603
# General.PC1 -0.08694    0.39853  -0.218

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC      LRT Pr(Chi)
# <none>         580.88                 
# General.PC1  1 578.92 0.047286  0.8279

# no effect of General soil traits on Total nodule mass


# *** NonNitrogen.PC2 model ####    

M1 <- lmer(Total_nodule_mass_mg ~ General.PC2 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

# look at parameter estimates
summary(M1)

# test main effect
drop1(M1, test = "Chisq")
# npar    AIC    LRT Pr(Chi)  
# <none>           577.68                 
# General.PC2    1 578.92 3.2463 0.07159 .


#### SYMPATRIC HOST / TOTAL NODULE COUNT ####

#### *** Data exploration ####

par(mfrow = c(1,1))
dotchart(symp.noAnz.live$Total_nodule_count, groups = symp.noAnz.live$Host_line)
dotchart(symp.noAnz.live$Total_nodule_count, groups = symp.noAnz.live$Nitrogen.PC1)
dotchart(symp.noAnz.live$Total_nodule_count, groups = symp.noAnz.live$General.PC1)
# no outliers

data = symp.noAnz.live$Total_nodule_count
summary(data)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.00   25.00   30.00   32.34   38.75   80.00 

# calculate mean and SE
mean(data, na.rm = TRUE)
# 32.34043
sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
# 1.445854

## *** usedata ####

usedata <- complete_fun(symp.noAnz.live, "Total_nodule_count")
usedata %>% nrow()
# 94

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 84

####  *** Nitrogen.PC1 model  ####

M1 <- lmer(Total_nodule_count ~ Nitrogen.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept)  25.702   5.070  
# Block     (Intercept)   7.628   2.762  
# Residual              140.410  11.849  
# Number of obs: 94, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   31.4176     2.2403  14.024
# Nitrogen.PC1   2.2110     0.9927   2.227

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)  
# <none>          755.34                 
# Nitrogen.PC1  1 757.34 4.0026 0.04543 *

# significant effect of nitrogen on total nodule count
# for every 1 unit increase in Nitrogen.PC1, nodule count increases by 2.21 nodules


#### *** Follow-up on Nitrogen model ####

# is positive effect of Nitrogen.PC1 on nodule count due to increase in intrinsic plant size along the N dep gradient?

# can test using total mass of sympatric plants in sterilized inoculum treatment

usedata <- complete_fun(symp.noAnz.ster, "logTotmass")
usedata %>% nrow()
# 93


M1 <- lmer(logTotmass ~ Nitrogen.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept) 0.15172  0.3895  
# Block     (Intercept) 0.01091  0.1044  
# Residual              0.07948  0.2819  
# Number of obs: 93, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   2.73170    0.13340  20.478
# Nitrogen.PC1  0.03829    0.06291   0.609


# test main effects
drop1(M1, test = "Chisq")
# Df    AIC     LRT Pr(Chi)
# <none>          74.760                
# Nitrogen.PC1  1 73.123 0.36365  0.5465


# no, nitrogen in sterilized soil does not predict total plant size 



####  *** NonNitrogen.PC1 model  ####

M1 <- lmer(Total_nodule_count ~ General.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept)  44.875   6.699  
# Block     (Intercept)   7.952   2.820  
# Residual              140.378  11.848  
# Number of obs: 94, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  32.2096     2.6174  12.306
# General.PC1   0.4932     1.3799   0.357


# test main effects
drop1(M1, test = "Chisq")
# Df    AIC     LRT Pr(Chi)
# <none>         759.21                
# General.PC1  1 757.34 0.12607  0.7225

# no effect of General soil traits on Total nodule count




# *** NonNitrogen.PC2 model ####    

M1 <- lmer(Total_nodule_count ~ General.PC2 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

# look at parameter estimates
summary(M1)

# test main effect
drop1(M1, test = "Chisq")
# npar    AIC    LRT  Pr(Chi)   
# <none>           751.49                   
# General.PC2    1 757.34 7.8477 0.005089 **


# ****** exploratory figure ####

p1 <- ggplot(usedata, aes(x = General.PC2, y = Total_nodule_count)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Total nodule count)") +
  xlab("General.PC2") +
  #facet_wrap(~Inoculum_type, scales = "free_x") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

print(p1)


#### SYMPATRIC HOST / MEAN NODULE SIZE  ####

#### *** Data exploration ####

par(mfrow = c(1,1))

dotchart(symp.noAnz.live$Nod_size_mg, groups = symp.noAnz.live$Host_line)
dotchart(symp.noAnz.live$Nod_size_mg, groups = symp.noAnz.live$Nitrogen.PC1)
dotchart(symp.noAnz.live$Nod_size_mg, groups = symp.noAnz.live$General.PC1)
# one potential outlier, greater than 2 mg

data = symp.noAnz.live$Nod_size_mg
summary(data)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0650  0.2312  0.3567  0.3918  0.4491  2.0571 

# calculate mean and SE
mean(data, na.rm = TRUE)
# 0.3917847
sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
# 0.02681531

## *** usedata ####

usedata <- complete_fun(symp.noAnz.live, "Nod_size_mg")
usedata %>% nrow()
# 94

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 84

usedata.o <- usedata[which(usedata$Plant_ID!="1041"),]
usedata.o %>% nrow()
# 93


####  *** Nitrogen.PC1 model  ####

M1 <- lmer(Nod_size_mg ~ Nitrogen.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept) 0.008338 0.09131 
# Block     (Intercept) 0.001825 0.04272 
# Residual              0.055562 0.23572 
# Number of obs: 94, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   0.39620    0.04092   9.683
# Nitrogen.PC1 -0.01676    0.01856  -0.903

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>          16.329               
# Nitrogen.PC1  1 15.113 0.7837   0.376

# no effect of nitrogen on nodule size
# interesting because nitrogen increased nodule count and nodule mass (marginal effect)
# but nod size does not change with nitrogen


# check if outlier was influential

M1o <- lmer(Nod_size_mg ~ Nitrogen.PC1 + (1|Host_line) + (1|Block),
           data = usedata.o, REML = FALSE)

# test main effects
drop1(M1o, test = "Chisq")
# Df     AIC     LRT Pr(Chi)
# <none>          -43.027                
# Nitrogen.PC1  1 -44.909 0.11708  0.7322

# outlier was not influential; even without it, nitrogen has no effect on nod size


####  *** NonNitrogen.PC1 model  ####

M1 <- lmer(Nod_size_mg ~ General.PC1 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups    Name        Variance Std.Dev.
# Host_line (Intercept) 0.007388 0.08595 
# Block     (Intercept) 0.002054 0.04532 
# Residual              0.055332 0.23523 
# Number of obs: 94, groups:  Host_line, 10; Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  0.38846    0.03940   9.859
# General.PC1 -0.02701    0.02095  -1.289


# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>         15.576               
# General.PC1  1 15.113 1.5374   0.215

# no effect of General soil traits on nodule size

# check if outlier was influential

M1o <- lmer(Nod_size_mg ~ General.PC1 + (1|Host_line) + (1|Block),
            data = usedata.o, REML = FALSE)

# test main effects
drop1(M1o, test = "Chisq")
# Df     AIC     LRT Pr(Chi)
# <none>         -43.697                
# General.PC1  1 -44.909 0.78744  0.3749

# outlier was not influential; even without it, general soil traits have no effect on nod size




# *** NonNitrogen.PC2 model ####    

M1 <- lmer(Nod_size_mg ~ General.PC2 + (1|Host_line) + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

# look at parameter estimates
summary(M1)

# test main effect
drop1(M1, test = "Chisq")
# npar    AIC     LRT Pr(Chi)
# <none>           16.267                
# General.PC2    1 15.113 0.84616  0.3576





#### COMMON HOST / TOTAL PLANT MASS ####
    
#### *** Data exploration ####
    
par(mfrow = c(1,1))

hist(univ.noAnz$Totmass_mg, col = "skyblue", breaks = 20)
# looks bimodal (due to live and ster inocula)

dotchart(univ.noAnz$Totmass_mg, groups = univ.noAnz$Inoculum_type)
dotchart(univ.noAnz$Totmass_mg, groups = univ.noAnz$Host_line)
# A_heermannii has more variance than Anz and Cla

summary(univ.noAnz$Totmass_mg)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.80   19.50   46.50   73.62  117.97  290.00 

data = univ.noAnz$Totmass_mg
mean(data, na.rm = TRUE)
# 73.61877
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 3.758739
    
    
#### *** usedata ####

usedata <- complete_fun(univ.noAnz, "logTotmass")
usedata %>% nrow
# 290

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 276
    
    
####  *** Nitrogen.PC1 model  ####
    

M1 <- lmer(logTotmass ~ Nitrogen.PC1 + Inoculum_type + Host_line +
             Nitrogen.PC1:Inoculum_type + Nitrogen.PC1:Host_line + Inoculum_type:Host_line + 
             Nitrogen.PC1:Inoculum_type:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)

summary(M1)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.009326 0.09657 
# Residual             0.103780 0.32215 
# Number of obs: 290, groups:  Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                              4.214277   0.044802  94.064
# PC1.N                                    0.043973   0.016039   2.742
# Inoculum_type.L                          1.214627   0.046360  26.200
# Host_lineAnz13_04                       -0.441048   0.047974  -9.193
# Host_lineCla12_04                       -0.818139   0.046360 -17.648
# PC1.N:Inoculum_type.L                   -0.001462   0.022683  -0.064
# PC1.N:Host_lineAnz13_04                 -0.035888   0.023381  -1.535
# PC1.N:Host_lineCla12_04                 -0.006184   0.022683  -0.273
# Inoculum_type.L:Host_lineAnz13_04       -0.132264   0.067845  -1.949
# Inoculum_type.L:Host_lineCla12_04        0.368104   0.065562   5.615
# PC1.N:Inoculum_type.L:Host_lineAnz13_04 -0.037481   0.033066  -1.134
# PC1.N:Inoculum_type.L:Host_lineCla12_04 -0.009546   0.032078  -0.298

      # for every 1 unit increase in PC1.N, there is a 0.043 unit increase in logTotMass
      # meaning that inocula from higher-N soils produce greater plant growth
      # but this could be due to the different baseline N levels in the inocula, not the microbes
      # because there is no PC1.N:Inoculum_type interaction

drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>                           206.81               
# PC1.N:Inoculum_type:Host_line  2 204.17 1.3648  0.5054

M1.no3way <- update(M1, .~. -Nitrogen.PC1:Inoculum_type:Host_line)

drop1(M1.no3way, test = "Chisq")
# Df    AIC    LRT   Pr(Chi)    
# <none>                     204.17                     
# PC1.N:Inoculum_type      1 203.65  1.480    0.2237    
# PC1.N:Host_line          2 203.09  2.921    0.2321    
# Inoculum_type:Host_line  2 257.25 57.079 4.032e-13 ***

M1.testN <- update(M1.no3way, .~. -Nitrogen.PC1:Inoculum_type -Nitrogen.PC1:Host_line)
M1.testtype <- update(M1.no3way, .~. -Nitrogen.PC1:Inoculum_type -Inoculum_type:Host_line)
M1.testhost <- update(M1.no3way, .~. -Nitrogen.PC1:Host_line -Inoculum_type:Host_line)

drop1(M1.testN, test = "Chisq")
# Df    AIC    LRT   Pr(Chi)    
# <none>                     202.78                     
# PC1.N                    1 210.38  9.602  0.001943 ** 
# Inoculum_type:Host_line  2 255.01 56.231 6.161e-13 ***

drop1(M1.testtype, test = "Chisq")
# Df    AIC    LRT Pr(Chi)    
# <none>             256.46                   
# Inoculum_type    1 831.53 577.07  <2e-16 ***
# PC1.N:Host_line  2 255.01   2.55  0.2791    

drop1(M1.testhost, test = "Chisq")
# Df    AIC     LRT Pr(Chi)    
# <none>                 255.63                    
# Host_line            2 437.03 185.392  <2e-16 ***
# PC1.N:Inoculum_type  1 255.01   1.377  0.2406    


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)

# test Inoculum_type:Host_Line interaction

# testing among host lines, within each inoculum type
M1 %>%
  emmeans(~ Host_line * Inoculum_type) %>%
  contrast(., method = "pairwise", by = "Inoculum_type") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                Inoculum_type estimate     SE  df t.ratio p.value
# A_heermannii - Anz13_04 Sterilized     0.35107 0.0704 295  4.989  <.0001 
# A_heermannii - Cla12_04 Sterilized     1.07821 0.0657 291 16.403  <.0001 
# Anz13_04 - Cla12_04     Sterilized     0.72714 0.0704 295 10.334  <.0001 
# A_heermannii - Anz13_04 Live           0.55817 0.0657 291  8.491  <.0001 
# A_heermannii - Cla12_04 Live           0.56274 0.0657 291  8.561  <.0001 
# Anz13_04 - Cla12_04     Live           0.00457 0.0657 291  0.070  0.9446 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 6 tests 

# testing among inoculum types, within each host line
M1 %>%
  emmeans(~ Host_line * Inoculum_type) %>%
  contrast(., method = "pairwise", by = "Host_line") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast          Host_line    estimate     SE  df t.ratio p.value
# Sterilized - Live A_heermannii    -1.72 0.0657 291 -26.120 <.0001 
# Sterilized - Live Anz13_04        -1.51 0.0704 295 -21.458 <.0001 
# Sterilized - Live Cla12_04        -2.23 0.0657 291 -33.962 <.0001 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 


####  *** NonNitrogen.PC1 model  ####

M1 <- lmer(logTotmass ~ General.PC1 + Inoculum_type + Host_line +
                 General.PC1:Inoculum_type + General.PC1:Host_line + Inoculum_type:Host_line + 
                 General.PC1:Inoculum_type:Host_line + (1|Block),
               data = usedata, REML = FALSE)
    
# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE
    
# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged
    
# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look great

summary(M1)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.009304 0.09646 
# Residual             0.104475 0.32323 
# Number of obs: 290, groups:  Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                              4.228969   0.044464  95.110
# General.PC1                                    0.024396   0.017530   1.392
# Inoculum_type.L                          1.214294   0.045754  26.540
# Host_lineAnz13_04                       -0.454506   0.047352  -9.598
# Host_lineCla12_04                       -0.821624   0.045754 -17.958
# General.PC1:Inoculum_type.L                   -0.002762   0.024791  -0.111
# General.PC1:Host_lineAnz13_04                 -0.001614   0.025554  -0.063
# General.PC1:Host_lineCla12_04                  0.014394   0.024791   0.581
# Inoculum_type.L:Host_lineAnz13_04       -0.141520   0.066966  -2.113
# Inoculum_type.L:Host_lineCla12_04        0.365322   0.064705   5.646
# General.PC1:Inoculum_type.L:Host_lineAnz13_04 -0.061633   0.036138  -1.705
# General.PC1:Inoculum_type.L:Host_lineCla12_04 -0.010424   0.035059  -0.297


drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>                           208.68               
# General.PC1:Inoculum_type:Host_line  2 207.93 3.2569  0.1962

M1.no3way <- update(M1, .~. -General.PC1:Inoculum_type:Host_line)

    drop1(M1.no3way, test = "Chisq")
# Df    AIC    LRT   Pr(Chi)    
# <none>                     207.93                     
# General.PC1:Inoculum_type      1 208.90  2.966   0.08502 .  
# General.PC1:Host_line          2 204.55  0.614   0.73581    
# Inoculum_type:Host_line  2 260.32 56.386 5.702e-13 ***
#   ---


M1.testG <- update(M1.no3way, .~. -General.PC1:Inoculum_type -General.PC1:Host_line)
M1.testtype <- update(M1.no3way, .~. -General.PC1:Inoculum_type -Inoculum_type:Host_line)
M1.testhost <- update(M1.no3way, .~. -General.PC1:Host_line -Inoculum_type:Host_line)

drop1(M1.testG, test = "Chisq")
# Df    AIC    LRT  Pr(Chi)    
# <none>                     205.61                    
# General.PC1                    1 210.38  6.773 0.009253 ** 
# Inoculum_type:Host_line  2 257.32 55.716 7.97e-13 ***

drop1(M1.testtype, test = "Chisq")
# Df    AIC    LRT Pr(Chi)    
# <none>             260.74                   
# Inoculum_type    1 832.10 573.36  <2e-16 ***
#   General.PC1:Host_line  2 257.32   0.58  0.7482    

drop1(M1.testhost, test = "Chisq")
# Df    AIC     LRT Pr(Chi)    
# <none>                 256.82                    
# Host_line            2 437.64 184.819  <2e-16 ***
#   General.PC1:Inoculum_type  1 257.32   2.504  0.1135    


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)

# test Inoculum_type:Host_Line interaction

# testing among host lines, within each inoculum type
M1 %>%
  emmeans(~ Host_line * Inoculum_type) %>%
  contrast(., method = "pairwise", by = "Inoculum_type") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                Inoculum_type estimate     SE  df t.ratio p.value
# A_heermannii - Anz13_04 Sterilized     0.35110 0.0706 295  4.973  <.0001 
# A_heermannii - Cla12_04 Sterilized     1.07821 0.0660 291 16.348  <.0001 
# Anz13_04 - Cla12_04     Sterilized     0.72712 0.0706 295 10.299  <.0001 
# A_heermannii - Anz13_04 Live           0.55817 0.0660 291  8.463  <.0001 
# A_heermannii - Cla12_04 Live           0.56274 0.0660 291  8.532  <.0001 
# Anz13_04 - Cla12_04     Live           0.00457 0.0660 291  0.069  0.9448 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 6 tests 


# testing among inoculum types, within each host line
M1 %>%
  emmeans(~ Host_line * Inoculum_type) %>%
  contrast(., method = "pairwise", by = "Host_line") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast          Host_line    estimate     SE  df t.ratio p.value
# Sterilized - Live A_heermannii    -1.72 0.0660 291 -26.033 <.0001 
# Sterilized - Live Anz13_04        -1.51 0.0706 295 -21.387 <.0001 
# Sterilized - Live Cla12_04        -2.23 0.0660 291 -33.849 <.0001 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 
    
    
####  *** NonNitrogen.PC2 model  ####
    
M1 <- lmer(logTotmass ~ General.PC2 + Inoculum_type + Host_line +
             General.PC2:Inoculum_type + General.PC2:Host_line + Inoculum_type:Host_line + 
             General.PC2:Inoculum_type:Host_line + (1|Block),
           data = usedata, REML = FALSE)
  
# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

summary(M1)

drop1(M1, test = "Chisq")
# npar    AIC     LRT Pr(Chi)
# <none>                                   216.41                
# General.PC2:Inoculum_type:Host_line    2 212.57 0.16068  0.9228


M1.no3way <- update(M1, .~. -General.PC2:Inoculum_type:Host_line)

drop1(M1.no3way, test = "Chisq")
# npar    AIC    LRT   Pr(Chi)    
# <none>                         212.57                     
# General.PC2:Inoculum_type    1 211.48  0.905    0.3415    
# General.PC2:Host_line        2 211.45  2.875    0.2375    
# Inoculum_type:Host_line      2 264.11 55.542 8.696e-13 ***

M1.testG <- update(M1.no3way, .~. -General.PC2:Inoculum_type -General.PC2:Host_line)
M1.testtype <- update(M1.no3way, .~. -General.PC2:Inoculum_type -Inoculum_type:Host_line)
M1.testhost <- update(M1.no3way, .~. -General.PC2:Host_line -Inoculum_type:Host_line)

drop1(M1.testG, test = "Chisq")
# npar    AIC    LRT   Pr(Chi)    
# <none>                       210.51                     
# General.PC2                1 210.38  1.874    0.1711    
# Inoculum_type:Host_line    2 261.34 54.834 1.239e-12 ***

drop1(M1.testtype, test = "Chisq")
# npar    AIC    LRT Pr(Chi)    
# <none>                     262.86                   
# Inoculum_type            1 832.38 571.53  <2e-16 ***
# General.PC2:Host_line    2 261.34   2.49  0.2883   

drop1(M1.testhost, test = "Chisq")
# npar    AIC     LRT Pr(Chi)    
# <none>                         262.47                    
# Host_line                    2 440.58 182.105  <2e-16 ***
# General.PC2:Inoculum_type    1 261.34   0.871  0.3506  
  
  

#### COMMON HOST / ROOT:SHOOT RATIO ####
    
#### *** Data exploration ####

par(mfrow = c(1,1))

hist(univ.noAnz$Root_shoot, col = "skyblue", breaks = 20)
# looks bimodal (due to live and ster inocula)

dotchart(univ.noAnz$Root_shoot, groups = univ.noAnz$Inoculum_type)
dotchart(univ.noAnz$Root_shoot, groups = univ.noAnz$Host_line)

summary(univ.noAnz$Root_shoot)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1179  0.2939  0.5141  0.5358  0.7289  1.3871 

data = univ.noAnz$Root_shoot
mean(data, na.rm = TRUE)
# 0.5358255
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.0152961

    
#### *** usedata ####

usedata <- complete_fun(univ.noAnz, "Root_shoot")
length(usedata$Totmass_mg)
# 290

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 276

    
####  *** Nitrogen.PC1 model  ####
    
M1 <- lmer(Root_shoot ~ Nitrogen.PC1 + Inoculum_type + Host_line +
                 Nitrogen.PC1:Inoculum_type + Nitrogen.PC1:Host_line + Inoculum_type:Host_line + 
                 Nitrogen.PC1:Inoculum_type:Host_line + (1|Block),
               data = usedata, REML = FALSE)
    
# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look ok

drop1(M1, test = "Chisq")
# Df     AIC    LRT Pr(Chi)
# <none>                           -363.34               
# Nitrogen.PC1:Inoculum_type:Host_line  2 -363.53 3.8093  0.1489

M1.no3way <- update(M1, .~. -Nitrogen.PC1:Inoculum_type:Host_line)

drop1(M1.no3way, test = "Chisq")
# Df     AIC     LRT   Pr(Chi)    
# <none>                     -363.53                      
# Nitrogen.PC1:Inoculum_type      1 -362.31  3.2206   0.07272 .  
# Nitrogen.PC1:Host_line          2 -365.61  1.9278   0.38141    
# Inoculum_type:Host_line  2 -342.63 24.9002 3.917e-06 ***

M1.testN <- update(M1.no3way, .~. -Nitrogen.PC1:Inoculum_type -Nitrogen.PC1:Host_line)
M1.testtype <- update(M1.no3way, .~. -Nitrogen.PC1:Inoculum_type -Inoculum_type:Host_line)
M1.testhost <- update(M1.no3way, .~. -Nitrogen.PC1:Host_line -Inoculum_type:Host_line)

drop1(M1.testN, test = "Chisq")
# Df     AIC     LRT   Pr(Chi)    
# <none>                     -364.47                      
# Nitrogen.PC1                    1 -364.12  2.3485    0.1254    
# Inoculum_type:Host_line  2 -344.00 24.4674 4.864e-06 ***

drop1(M1.testtype, test = "Chisq")
# Df     AIC    LRT Pr(Chi)    
# <none>             -341.69                   
# Inoculum_type    1   18.30 361.99  <2e-16 ***
# Nitrogen.PC1:Host_line  2 -344.00   1.69  0.4306    

drop1(M1.testhost, test = "Chisq")
# Df     AIC     LRT Pr(Chi)    
# <none>                 -344.87                    
# Host_line            2 -248.83 100.044 < 2e-16 ***
# Nitrogen.PC1:Inoculum_type  1 -344.00   2.868 0.09038 .  


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)

# test Inoculum_type:Host_Line interaction

# testing among host lines, within each inoculum type
M1 %>%
  emmeans(~ Host_line * Inoculum_type) %>%
  contrast(., method = "pairwise", by = "Inoculum_type") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# contrast                Inoculum_type estimate     SE  df t.ratio p.value
# A_heermannii - Anz13_04 Sterilized      0.2516 0.0267 299  9.411  <.0001 
# A_heermannii - Cla12_04 Sterilized      0.0878 0.0251 291  3.500  0.0011 
# Anz13_04 - Cla12_04     Sterilized     -0.1639 0.0267 299 -6.128  <.0001 
# A_heermannii - Anz13_04 Live            0.1617 0.0251 291  6.449  <.0001 
# A_heermannii - Cla12_04 Live            0.1826 0.0251 291  7.282  <.0001 
# Anz13_04 - Cla12_04     Live            0.0209 0.0251 291  0.833  0.4056 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 6 tests 

# testing among inoculum types, within each host line
M1 %>%
  emmeans(~ Host_line * Inoculum_type) %>%
  contrast(., method = "pairwise", by = "Host_line") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# contrast          Host_line    estimate     SE  df t.ratio p.value
# Sterilized - Live A_heermannii    0.406 0.0251 291 16.173  <.0001 
# Sterilized - Live Anz13_04        0.316 0.0267 299 11.806  <.0001 
# Sterilized - Live Cla12_04        0.500 0.0251 291 19.955  <.0001 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 


####  *** NonNitrogen.PC1 model  ####
    
M1 <- lmer(Root_shoot ~ General.PC1 + Inoculum_type + Host_line +
                 General.PC1:Inoculum_type + General.PC1:Host_line + Inoculum_type:Host_line + 
                 General.PC1:Inoculum_type:Host_line + (1|Block),
               data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)

summary(M1)
# Random effects:
#   Groups   Name        Variance  Std.Dev.
# Block    (Intercept) 7.509e-05 0.008666
# Residual             1.532e-02 0.123765
# Number of obs: 290, groups:  Block, 10
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                              0.652060   0.012688  51.394
# General.PC1                                   -0.001915   0.006712  -0.285
# Inoculum_type.L                         -0.285470   0.017519 -16.295
# Host_lineAnz13_04                       -0.207095   0.018070 -11.461
# Host_lineCla12_04                       -0.135299   0.017519  -7.723
# General.PC1:Inoculum_type.L                   -0.016853   0.009493  -1.775
# General.PC1:Host_lineAnz13_04                  0.005659   0.009785   0.578
# General.PC1:Host_lineCla12_04                  0.001294   0.009493   0.136
# Inoculum_type.L:Host_lineAnz13_04        0.063625   0.025555   2.490
# Inoculum_type.L:Host_lineCla12_04       -0.067677   0.024776  -2.732
# General.PC1:Inoculum_type.L:Host_lineAnz13_04 -0.001250   0.013838  -0.090
# General.PC1:Inoculum_type.L:Host_lineCla12_04  0.007779   0.013424   0.579

# sig General.PC1:inoculum_type effect
# with live inocula, one unit increase in General.PC1 reduces root:shoot by 0.017
# but effect does not occur for sterilized inocula


drop1(M1, test = "Chisq")
# Df     AIC     LRT Pr(Chi)
# <none>                           -359.52                
# General.PC1:Inoculum_type:Host_line  2 -363.01 0.51563  0.7727

M1.no3way <- update(M1, .~. -General.PC1:Inoculum_type:Host_line)

drop1(M1.no3way, test = "Chisq")
# Df     AIC     LRT   Pr(Chi)    
# <none>                     -363.01                      
# General.PC1:Inoculum_type      1 -358.32  6.6896  0.009698 ** 
# General.PC1:Host_line          2 -366.68  0.3255  0.849797    
# Inoculum_type:Host_line  2 -342.15 24.8545 4.008e-06 ***

M1.testG <- update(M1.no3way, .~. -General.PC1:Inoculum_type -General.PC1:Host_line)
M1.testtype <- update(M1.no3way, .~. -General.PC1:Inoculum_type -Inoculum_type:Host_line)
M1.testhost <- update(M1.no3way, .~. -General.PC1:Host_line -Inoculum_type:Host_line)

drop1(M1.testG, test = "Chisq")
# Df     AIC     LRT   Pr(Chi)    
# <none>                     -362.12                      
# General.PC1                    1 -364.12  0.0015    0.9692    
# Inoculum_type:Host_line  2 -341.85 24.2701 5.368e-06 ***

drop1(M1.testtype, test = "Chisq")
# Df     AIC    LRT Pr(Chi)    
# <none>             -338.03                   
# Inoculum_type    1   19.36 359.39  <2e-16 ***
# General.PC1:Host_line  2 -341.85   0.18  0.9162    

drop1(M1.testhost, test = "Chisq")
# Df     AIC     LRT Pr(Chi)    
# <none>                 -345.85                    
# Host_line            2 -249.52 100.330 < 2e-16 ***
#   General.PC1:Inoculum_type  1 -341.85   5.999 0.01432 *  


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)

# test Inoculum_type:Host_Line interaction

# testing among host lines, within each inoculum type
M1 %>%
  emmeans(~ Host_line * Inoculum_type) %>%
  contrast(., method = "pairwise", by = "Inoculum_type") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# NOTE: Results may be misleading due to involvement in interactions
# contrast                Inoculum_type estimate     SE  df t.ratio p.value
# A_heermannii - Anz13_04 Sterilized      0.2516 0.0269 299  9.344  <.0001 
# A_heermannii - Cla12_04 Sterilized      0.0878 0.0253 291  3.476  0.0012 
# Anz13_04 - Cla12_04     Sterilized     -0.1638 0.0269 299 -6.084  <.0001 
# A_heermannii - Anz13_04 Live            0.1617 0.0253 291  6.405  <.0001 
# A_heermannii - Cla12_04 Live            0.1826 0.0253 291  7.232  <.0001 
# Anz13_04 - Cla12_04     Live            0.0209 0.0253 291  0.827  0.4088 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 6 tests 

# testing among inoculum types, within each host line
M1 %>%
  emmeans(~ Host_line * Inoculum_type) %>%
  contrast(., method = "pairwise", by = "Host_line") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast          Host_line    estimate     SE  df t.ratio p.value
# Sterilized - Live A_heermannii    0.406 0.0253 291 16.063  <.0001 
# Sterilized - Live Anz13_04        0.316 0.0269 299 11.730  <.0001 
# Sterilized - Live Cla12_04        0.500 0.0253 291 19.819  <.0001 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 


    
####  *** NonNitrogen.PC2 model  ####

M1 <- lmer(Root_shoot ~ General.PC2 + Inoculum_type + Host_line +
             General.PC2:Inoculum_type + General.PC2:Host_line + Inoculum_type:Host_line + 
             General.PC2:Inoculum_type:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look ok

summary(M1)

drop1(M1, test = "Chisq")
# npar     AIC    LRT Pr(Chi)  
# <none>                                   -367.41                 
# General.PC2:Inoculum_type:Host_line    2 -364.33 7.0824 0.02898 *

M1.no3way <- update(M1, .~. -General.PC2:Inoculum_type:Host_line)

drop1(M1.no3way, test = "Chisq")
# npar     AIC     LRT   Pr(Chi)    
# <none>                         -364.33                      
# General.PC2:Inoculum_type    1 -365.29  1.0324    0.3096    
# General.PC2:Host_line        2 -366.69  1.6406    0.4403    
# Inoculum_type:Host_line      2 -343.36 24.9689 3.785e-06 ***

M1.testG <- update(M1.no3way, .~. -General.PC2:Inoculum_type -General.PC2:Host_line)
M1.testtype <- update(M1.no3way, .~. -General.PC2:Inoculum_type -Inoculum_type:Host_line)
M1.testhost <- update(M1.no3way, .~. -General.PC2:Host_line -Inoculum_type:Host_line)

drop1(M1.testG, test = "Chisq")
# npar     AIC     LRT   Pr(Chi)    
# <none>                       -367.59                      
# General.PC2                1 -364.12  5.4707   0.01934 *  
# Inoculum_type:Host_line    2 -346.86 24.7332 4.259e-06 ***

drop1(M1.testtype, test = "Chisq")
# npar     AIC    LRT Pr(Chi)    
# <none>                     -344.41                   
# Inoculum_type            1   17.52 363.94  <2e-16 ***
# General.PC2:Host_line    2 -346.86   1.55  0.4597 

drop1(M1.testhost, test = "Chisq")
# npar     AIC     LRT Pr(Chi)    
# <none>                         -345.86                    
# Host_line                    2 -249.53 100.332  <2e-16 ***
# General.PC2:Inoculum_type    1 -346.86   0.998  0.3177    





# ****** exploratory figure ####

p1 <- ggplot(usedata, aes(x = General.PC2, y = Root_shoot, fill = Inoculum_type)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Root:shoot ratio)") +
  xlab("General.PC2") +
  facet_wrap(~Host_line, scales = "free_x") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "right")
print(p1)



    
#### COMMON HOST / PLANT RELATIVE GROWTH ####
    
#### *** Data exploration ####

par(mfrow = c(1,1))

hist(univ.noAnz.live$Relgro_totmass, col = "skyblue")
# looks fairly right-skewed

dotchart(univ.noAnz.live$Relgro_totmass, groups = univ.noAnz.live$Host_line)
dotchart(univ.noAnz.live$Relgro_totmass, groups = univ.noAnz.live$Block)
# may be one high outlier; plant ID = 1453

summary(univ.noAnz.live$Relgro_totmass)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.017   4.594   6.056   7.231   8.812  36.250 

data = univ.noAnz.live$Relgro_totmass
mean(data, na.rm = TRUE)
# 7.231113
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.3738366

# summarize relative growth by host line

# A. heermannii
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="A_heermannii"),]
mean(usedata$Relgro_totmass, na.rm = TRUE)
# 5.987835
se <- sd(usedata$Relgro_totmass, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Relgro_totmass)))
se
# 0.3204444

# Anz13_04
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="Anz13_04"),]
mean(usedata$Relgro_totmass, na.rm = TRUE)
# 5.118389
se <- sd(usedata$Relgro_totmass, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Relgro_totmass)))
se
# 0.3590152

# Cla12_04
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="Cla12_04"),]
mean(usedata$Relgro_totmass, na.rm = TRUE)
# 10.58744
se <- sd(usedata$Relgro_totmass, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Relgro_totmass)))
se
# 0.8312132


#### *** usedata ####

length(univ.noAnz.live$Relgro_totmass)
# 150

usedata <- complete_fun(univ.noAnz.live, "Relgro_totmass")
length(usedata$Relgro_totmass)
# 150

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 139

usedata.o <- usedata[which(usedata$Plant_ID!="1453"),]


####  *** Nitrogen.PC1 model  ####

M1 <- lmer(Relgro_totmass ~ Nitrogen.PC1 + Host_line + Nitrogen.PC1:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look good

drop1(M1, test = "Chisq")
# Df    AIC   LRT Pr(Chi)
# <none>             842.97              
# Nitrogen.PC1:Host_line  2 840.08 1.112  0.5735

M1.noint <- update(M1, .~. -Nitrogen.PC1:Host_line)

drop1(M1.noint, test = "Chisq")
# Df    AIC    LRT   Pr(Chi)    
# <none>       840.08                     
# Nitrogen.PC1      1 840.79  2.712   0.09959 .  
# Host_line  2 886.99 50.913 8.798e-12 ***


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)


# no p-value adjustment
M1 %>%
  emmeans(~ Host_line) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                estimate    SE  df t.ratio p.value
# A_heermannii - Anz13_04    0.869 0.761 145   1.143  0.2551
# A_heermannii - Cla12_04   -4.600 0.761 145  -6.045  <.0001
# Anz13_04 - Cla12_04       -5.469 0.761 145  -7.187  <.0001
# 
# Degrees-of-freedom method: kenward-roger 

# Holm p-value adjustment (sequential)
M1 %>%
  emmeans(~ Host_line) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                estimate    SE  df t.ratio p.value
# A_heermannii - Anz13_04    0.869 0.761 145   1.143  0.2551
# A_heermannii - Cla12_04   -4.600 0.761 145  -6.045  <.0001
# Anz13_04 - Cla12_04       -5.469 0.761 145  -7.187  <.0001
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests


# was outlier significant?
M1o <- lmer(Relgro_totmass ~ Nitrogen.PC1 + Host_line + Nitrogen.PC1:Host_line + (1|Block),
           data = usedata.o, REML = FALSE)

drop1(M1o, test = "Chisq")
# npar    AIC    LRT Pr(Chi)
# <none>                      788.53               
# Nitrogen.PC1:Host_line    2 786.00 1.4774  0.4777

M1o.noint <- update(M1o, .~. - Nitrogen.PC1:Host_line)

drop1(M1o.noint, test = "Chisq")
# npar    AIC    LRT   Pr(Chi)    
# <none>            786.00                     
# Nitrogen.PC1    1 785.08  1.080    0.2988    
# Host_line       2 836.48 54.476 1.482e-12 ***

# outlier was not influential



####  *** NonNitrogen.PC1 model  ####

M1 <- lmer(Relgro_totmass ~ General.PC1 + Host_line + General.PC1:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)

drop1(M1, test = "Chisq")
# npar    AIC    LRT Pr(Chi)
# <none>                     840.72               
# General.PC1:Host_line    2 839.01 2.2844  0.3191

M1.noint <- update(M1, .~. -General.PC1:Host_line)

drop1(M1.noint, test = "Chisq")
# npar    AIC    LRT   Pr(Chi)    
# <none>           839.01                     
# General.PC1    1 840.82  3.814   0.05082 .  
# Host_line      2 886.24 51.238 7.477e-12 ***


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)


# no p-value adjustment
M1 %>%
  emmeans(~ Host_line) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                estimate    SE  df t.ratio p.value
# A_heermannii - Anz13_04    0.869 0.755 145   1.152  0.2513
# A_heermannii - Cla12_04   -4.600 0.755 145  -6.094  <.0001
# Anz13_04 - Cla12_04       -5.469 0.755 145  -7.246  <.0001
# 
# Degrees-of-freedom method: kenward-roger 

# Holm p-value adjustment (sequential)
M1 %>%
  emmeans(~ Host_line) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# contrast                estimate    SE  df t.ratio p.value
# A_heermannii - Anz13_04    0.869 0.755 145   1.152  0.2513
# A_heermannii - Cla12_04   -4.600 0.755 145  -6.094  <.0001
# Anz13_04 - Cla12_04       -5.469 0.755 145  -7.246  <.0001

# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 

# was outlier significant?

M1o <- lmer(Relgro_totmass ~ General.PC1 + Host_line + General.PC1:Host_line + (1|Block),
            data = usedata.o, REML = FALSE)

drop1(M1o, test = "Chisq")
# npar    AIC    LRT Pr(Chi)
# <none>                     785.96               
# General.PC1:Host_line    2 784.52 2.5607  0.2779

M1o.noint <- update(M1o, .~. - General.PC1:Host_line)

drop1(M1o.noint, test = "Chisq")
# npar    AIC    LRT   Pr(Chi)    
# <none>           784.52                     
# General.PC1    1 785.08  2.561    0.1095    
# Host_line      2 835.51 54.985 1.148e-12 ***

# outlier was not influential


    
####  *** NonNitrogen.PC2 model  ####

M1 <- lmer(Relgro_totmass ~ General.PC2 + Host_line + General.PC2:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

drop1(M1, test = "Chisq")
# npar    AIC    LRT Pr(Chi)
# <none>                     845.27               
# General.PC2:Host_line    2 841.41 0.1384  0.9331

M1.noint <- update(M1, .~. -General.PC2:Host_line)

drop1(M1.noint, test = "Chisq")
# npar    AIC    LRT   Pr(Chi)    
# <none>           841.41                     
# General.PC2    1 840.82  1.410     0.235    
# Host_line      2 887.92 50.506 1.078e-11 ***


    
#### COMMON HOST / RED NODULE FREQUENCY  ####
    
#### *** Data exploration ####

par(mfrow = c(1,1))

hist(univ.noAnz.live$Red_nod_freq, col = "skyblue")
# looks fairly normal

dotchart(univ.noAnz.live$Red_nod_freq, groups = univ.noAnz.live$Host_line)
dotchart(univ.noAnz.live$Red_nod_freq, groups = univ.noAnz.live$Block)
# variance looks pretty even; no outliers

summary(univ.noAnz.live$Red_nod_freq)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.04167 0.40741 0.57273 0.56501 0.75000 0.96667       1 

data = univ.noAnz.live$Red_nod_freq
mean(data, na.rm = TRUE)
# 0.5650105
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.01857959


# summarize red nodule freq by host line

# A. heermannii
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="A_heermannii"),]
mean(usedata$Red_nod_freq, na.rm = TRUE)
# 0.5139612
se <- sd(usedata$Red_nod_freq, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Red_nod_freq)))
se
# 0.02512109

# Anz13_04
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="Anz13_04"),]
mean(usedata$Red_nod_freq, na.rm = TRUE)
# 0.4902547
se <- sd(usedata$Red_nod_freq, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Red_nod_freq)))
se
# 0.03163941

# Cla12_04
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="Cla12_04"),]
mean(usedata$Red_nod_freq, na.rm = TRUE)
# 0.6933832
se <- sd(usedata$Red_nod_freq, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Red_nod_freq)))
se
# 0.03192025


#### *** usedata ####

usedata <- complete_fun(univ.noAnz.live, "Red_nod_freq")
length(usedata$Red_nod_freq)
# 149

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 138


# *** figure: effect of early-harvested plants ####

# full dataset
p1 <- ggplot(usedata, aes(x = General.PC1, y = Red_nod_freq)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Red nodule frequency") +
  xlab("General.PC1") +
  facet_wrap(~Host_line, scales = "free_x") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")
p1 <- arrangeGrob(p1, top = textGrob("All plants", x = unit(0, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=18)))

# excluding 27 early-harvest plants
p2 <- ggplot(usedata.e, aes(x = General.PC1, y = Red_nod_freq)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Red nodule frequency") +
  xlab("General.PC1") +
  facet_wrap(~Host_line, scales = "free_x") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")
p2 <- arrangeGrob(p2, top = textGrob("Minus early harvest plants", x = unit(0, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=18)))

grid.arrange(p1, p2, nrow = 1)


####  *** Nitrogen.PC1 model  ####
    
M1 <- lmer(Red_nod_freq ~ Nitrogen.PC1 + Host_line + Nitrogen.PC1:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

drop1(M1, test = "Chisq")
# Df     AIC   LRT Pr(Chi)
# <none>             -100.11              
# Nitrogen.PC1:Host_line  2 -100.85 3.261  0.1958

M1.noint <- update(M1, .~. -Nitrogen.PC1:Host_line)

drop1(M1.noint, test = "Chisq")
# Nitrogen.PC1      1 -102.420  0.433    0.5105    
# Host_line  2  -59.349 45.504 1.315e-10 ***


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)


# Holm p-value adjustment (sequential)
M1 %>%
  emmeans(~ Host_line) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                estimate     SE  df t.ratio p.value
# A_heermannii - Anz13_04   0.0235 0.0305 144  0.768  0.4438 
# A_heermannii - Cla12_04  -0.1818 0.0307 144 -5.920  <.0001 
# Anz13_04 - Cla12_04      -0.2053 0.0307 144 -6.683  <.0001 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 



####  *** NonNitrogen.PC1 model  ####
    
M1 <- lmer(Red_nod_freq ~ General.PC1 + Host_line + General.PC1:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

drop1(M1, test = "Chisq")
# Df     AIC    LRT Pr(Chi)  
# <none>             -102.00                 
# General.PC1:Host_line  2 -101.13 4.8653  0.0878 .

M1.noint <- update(M1, .~. -General.PC1:Host_line)

drop1(M1.noint, test = "Chisq")
# Df      AIC    LRT   Pr(Chi)    
# <none>       -101.131                     
# General.PC1      1 -102.420  0.711    0.3991    
# Host_line  2  -59.484 45.647 1.224e-10 ***

# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)


# Holm p-value adjustment (sequential)
M1 %>%
emmeans(~ Host_line) %>%
contrast(., method = "pairwise") %>%
summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                estimate     SE  df t.ratio p.value
# A_heermannii - Anz13_04   0.0236 0.0303 144  0.779  0.4370 
# A_heermannii - Cla12_04  -0.1816 0.0305 144 -5.954  <.0001 
# Anz13_04 - Cla12_04      -0.2052 0.0305 144 -6.729  <.0001 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 


####  *** NonNitrogen.PC2 model  ####

M1 <- lmer(Red_nod_freq ~ General.PC2 + Host_line + General.PC2:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

drop1(M1, test = "Chisq")
# npar      AIC   LRT Pr(Chi)
# <none>                      -97.815              
# General.PC2:Host_line    2 -100.434 1.381  0.5013

M1.noint <- update(M1, .~. -General.PC2:Host_line)

drop1(M1.noint, test = "Chisq")
# npar      AIC    LRT   Pr(Chi)    
# <none>           -100.434                     
# General.PC2    1 -102.420  0.014    0.9046    
# Host_line      2  -58.993 45.441 1.357e-10 ***


    
#### COMMON HOST / TOTAL NODULE MASS ####
    
#### *** Data exploration ####

par(mfrow = c(1,1))

hist(univ.noAnz.live$Total_nodule_mass_mg, col = "skyblue")
# looks fairly normal

dotchart(univ.noAnz.live$Total_nodule_mass_mg, groups = univ.noAnz.live$Host_line)
dotchart(univ.noAnz.live$Total_nodule_mass_mg, groups = univ.noAnz.live$Block)
# variance looks pretty even; no outliers

summary(univ.noAnz.live$Total_nodule_mass_mg)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.60    8.65   11.90   11.91   15.28   26.10 

data = univ.noAnz.live$Total_nodule_mass_mg
mean(data, na.rm = TRUE)
# 11.90867
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.3773627

# summarize data by host line

# A. heermannii
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="A_heermannii"),]
mean(usedata$Total_nodule_mass_mg, na.rm = TRUE)
# 14.068
se <- sd(usedata$Total_nodule_mass_mg, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Total_nodule_mass_mg)))
se
# 0.6699951

# Anz13_04
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="Anz13_04"),]
mean(usedata$Total_nodule_mass_mg, na.rm = TRUE)
# 10.622
se <- sd(usedata$Total_nodule_mass_mg, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Total_nodule_mass_mg)))
se
# 0.6189411

# Cla12_04
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="Cla12_04"),]
mean(usedata$Total_nodule_mass_mg, na.rm = TRUE)
# 11.036
se <- sd(usedata$Total_nodule_mass_mg, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Total_nodule_mass_mg)))
se
# 0.5679736


#### *** usedata ####

usedata <- complete_fun(univ.noAnz.live, "Total_nodule_mass_mg")
length(usedata$Total_nodule_mass_mg)
# 150

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 139

    
####  *** Nitrogen.PC1 model  ####
    
M1 <- lmer(Total_nodule_mass_mg ~ Nitrogen.PC1 + Host_line + Nitrogen.PC1:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look good

drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>             869.96               
# Nitrogen.PC1:Host_line  2 868.28 2.3141  0.3144

M1.noint <- update(M1, .~. -Nitrogen.PC1:Host_line)

drop1(M1.noint, test = "Chisq")
# Df    AIC     LRT   Pr(Chi)    
# <none>       868.28                      
# Nitrogen.PC1      1 869.43  3.1555   0.07567 .  
# Host_line  2 884.30 20.0215 4.492e-05 ***


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)


# Holm p-value adjustment (sequential)
M1 %>%
  emmeans(~ Host_line) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                estimate    SE  df t.ratio p.value
# A_heermannii - Anz13_04    3.446 0.819 145  4.207  0.0001 
# A_heermannii - Cla12_04    3.032 0.819 145  3.701  0.0006 
# Anz13_04 - Cla12_04       -0.414 0.819 145 -0.505  0.6141 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 


####  *** NonNitrogen.PC1 model  #### 

M1 <- lmer(Total_nodule_mass_mg ~ General.PC1 + Host_line + General.PC1:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>             873.50               
# General.PC1:Host_line  2 871.31 1.8009  0.4064

M1.noint <- update(M1, .~. -General.PC1:Host_line)

drop1(M1.noint, test = "Chisq")
# Df    AIC     LRT   Pr(Chi)    
# <none>       871.31                      
# General.PC1      1 869.43  0.1268    0.7218    
# Host_line  2 886.93 19.6216 5.486e-05 ***


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)


# Holm p-value adjustment (sequential)
M1 %>%
emmeans(~ Host_line) %>%
contrast(., method = "pairwise") %>%
summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                estimate   SE  df t.ratio p.value
# A_heermannii - Anz13_04    3.446 0.83 145  4.154  0.0002 
# A_heermannii - Cla12_04    3.032 0.83 145  3.655  0.0007 
# Anz13_04 - Cla12_04       -0.414 0.83 145 -0.499  0.6185 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 


####  *** NonNitrogen.PC2 model  #### 

M1 <- lmer(Total_nodule_mass_mg ~ General.PC2 + Host_line + General.PC2:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look good

drop1(M1, test = "Chisq")
# npar    AIC    LRT Pr(Chi)
# <none>                     870.06               
# General.PC2:Host_line    2 867.16 1.0906  0.5797

M1.noint <- update(M1, .~. -General.PC2:Host_line)

drop1(M1.noint, test = "Chisq")
# npar    AIC     LRT   Pr(Chi)    
# <none>           867.16                      
# General.PC2    1 869.43  4.2767   0.03864 *  
# Host_line      2 883.33 20.1714 4.167e-05 ***


# ****** exploratory figure ####

p1 <- ggplot(usedata, aes(x = General.PC2, y =Total_nodule_mass_mg)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Total nodule mass (mg))") +
  xlab("General.PC2") +
  #facet_wrap(~Host_line, scales = "free_x") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "right")

print(p1)


#### COMMON HOST / TOTAL NODULE COUNT  ####
    
#### *** Data exploration ####

par(mfrow = c(1,1))

hist(univ.noAnz.live$Total_nodule_count, col = "skyblue")
# looks fairly normal, a little right-skewed

dotchart(univ.noAnz.live$Total_nodule_count, groups = univ.noAnz.live$Host_line)
dotchart(univ.noAnz.live$Total_nodule_count, groups = univ.noAnz.live$Block)
# variance looks pretty even
# potentially one high outlier, plant ID 1515

summary(univ.noAnz.live$Total_nodule_count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.00   28.00   35.00   38.67   46.00  110.00 

data = univ.noAnz.live$Total_nodule_count
mean(data, na.rm = TRUE)
# 38.67333
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 1.292363

# summarize data by host line

# A. heermannii
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="A_heermannii"),]
mean(usedata$Total_nodule_count, na.rm = TRUE)
# 48.58
se <- sd(usedata$Total_nodule_count, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Total_nodule_count)))
se
# 2.650081

# Anz13_04
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="Anz13_04"),]
mean(usedata$Total_nodule_count, na.rm = TRUE)
# 33.1
se <- sd(usedata$Total_nodule_count, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Total_nodule_count)))
se
# 1.52349

# Cla12_04
usedata <- univ.noAnz.live[which(univ.noAnz.live$Host_line=="Cla12_04"),]
mean(usedata$Total_nodule_count, na.rm = TRUE)
# 34.34
se <- sd(usedata$Total_nodule_count, na.rm = TRUE)/sqrt(sum(!is.na(usedata$Total_nodule_count)))
se
# 1.694796
    
    
#### *** usedata ####
    
length(univ.noAnz.live$Total_nodule_count)
# 150

usedata <- complete_fun(univ.noAnz.live, "Total_nodule_count")
length(usedata$Total_nodule_count)
# 150

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 139

usedata.o <- usedata[which(usedata$Plant_ID!="1515"),]

    
####  *** Nitrogen.PC1 model  ####
    
M1 <- lmer(Total_nodule_count ~ Nitrogen.PC1 + Host_line + Nitrogen.PC1:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#TRUE

summary(M1)
# block variance estimated as zero

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>             1231.3               
# Nitrogen.PC1:Host_line  2 1229.8 2.5586  0.2782

M1.noint <- update(M1, .~. -Nitrogen.PC1:Host_line)

drop1(M1.noint, test = "Chisq")
# Df    AIC    LRT   Pr(Chi)    
# <none>       1229.8                     
# Nitrogen.PC1      1 1230.1  2.233    0.1351    
# Host_line  2 1259.4 33.582 5.102e-08 ***


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)


# Holm p-value adjustment (sequential)
M1 %>%
  emmeans(~ Host_line) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                estimate   SE  df t.ratio p.value
# A_heermannii - Anz13_04    15.48 2.83 145  5.467  <.0001 
# A_heermannii - Cla12_04    14.24 2.83 145  5.029  <.0001 
# Anz13_04 - Cla12_04        -1.24 2.83 145 -0.438  0.6621 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 


# was outlier influential?

M1o <- lmer(Total_nodule_count ~ Nitrogen.PC1 + Host_line + Nitrogen.PC1:Host_line + (1|Block),
            data = usedata.o, REML = FALSE)    

# check singularity
model = M1o
isSingular(model, tol = 1e-05)
# FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1o
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1o)
plot(simulationOutput)
# good residuals

drop1(M1o, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>             1200.0               
# Nitrogen.PC1:Host_line  2 1200.6 4.5895  0.1008

M1o.noint <- update(M1o, .~. -Nitrogen.PC1:Host_line)

drop1(M1o.noint, test = "Chisq")
# Df    AIC    LRT   Pr(Chi)    
# <none>       1200.6                     
# Nitrogen.PC1      1 1202.4  3.805   0.05108 .  
# Host_line  2 1228.6 32.014 1.118e-07 ***

# dropping plant 1515 makes Nitrogen.PC1 almost significant
# so it is sort of influential, but not enough to justify removing it


####  *** NonNitrogen.PC1 model  ####
    
M1 <- lmer(Total_nodule_count ~ General.PC1 + Host_line + General.PC1:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#TRUE

summary(M1)
# block estimated as zero

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged? output = 0

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>             1232.2               
# General.PC1:Host_line  2 1231.0 2.7237  0.2562

M1.noint <- update(M1, .~. -General.PC1:Host_line)

drop1(M1.noint, test = "Chisq")
# Df    AIC    LRT Pr(Chi)    
# <none>       1231.0                   
# General.PC1      1 1230.1  1.125  0.2889    
# Host_line  2 1260.3 33.360 5.7e-08 ***


# multiple comparisons with p-value adjustment

# see https://stats.stackexchange.com/questions/424304/when-to-correct-for-multiple-comparisons-with-specific-reference-to-emmeans-in
# the by = NULL term in summary() is important; if set to "eStrain" instead of NULL,
# the Holm correction will be done for each contrast separately (i.e., no correction will be done)


# Holm p-value adjustment (sequential)
M1 %>%
  emmeans(~ Host_line) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# contrast                estimate   SE  df t.ratio p.value
# A_heermannii - Anz13_04    15.48 2.84 145  5.449  <.0001 
# A_heermannii - Cla12_04    14.24 2.84 145  5.012  <.0001 
# Anz13_04 - Cla12_04        -1.24 2.84 145 -0.436  0.6632 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: holm method for 3 tests 


# was outlier influential?

M1o <- lmer(Total_nodule_count ~ General.PC1 + Host_line + General.PC1:Host_line + (1|Block),
            data = usedata.o, REML = FALSE)

# check singularity
model = M1o
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1o
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1o)
plot(simulationOutput)
# good residuals

drop1(M1o, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>             1204.5               
# General.PC1:Host_line  2 1202.4 1.8702  0.3925

M1o.noint <- update(M1o, .~. -General.PC1:Host_line)

drop1(M1o.noint, test = "Chisq")
# Df    AIC    LRT   Pr(Chi)    
# <none>       1202.4                     
# General.PC1      1 1202.4  2.009    0.1563    
# Host_line  2 1230.2 31.730 1.288e-07 ***

# outlier not influential

    
####  *** NonNitrogen.PC2 model  ####

M1 <- lmer(Total_nodule_count ~ General.PC2 + Host_line + General.PC2:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

drop1(M1, test = "Chisq")
# npar    AIC    LRT Pr(Chi)  
# <none>                     1227.8                 
# General.PC2:Host_line    2 1229.7 5.8409 0.05391 .

M1.noint <- update(M1, .~. -General.PC2:Host_line)

drop1(M1.noint, test = "Chisq")
# npar    AIC    LRT   Pr(Chi)    
# <none>           1229.7                     
# General.PC2    1 1230.1  2.402    0.1211    
# Host_line      2 1259.3 33.616 5.016e-08 ***

    
#### COMMON HOST / MEAN NODULE SIZE  ####
    
#### *** Data exploration ####

par(mfrow = c(1,1))

hist(univ.noAnz.live$Nod_size_mg, col = "skyblue")
# looks fairly normal, one possible outlier

dotchart(univ.noAnz.live$Nod_size_mg, groups = univ.noAnz.live$Host_line)
dotchart(univ.noAnz.live$Nod_size_mg, groups = univ.noAnz.live$Block)
# variance looks pretty even
# potentially one high outlier, plant ID 1285

summary(univ.noAnz.live$Nod_size_mg)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.08286 0.24470 0.29744 0.32776 0.38755 1.00476

data = univ.noAnz.live$Nod_size_mg
mean(data, na.rm = TRUE)
# 0.3277597
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.01047956
    
#### *** usedata ####

usedata <- complete_fun(univ.noAnz.live, "Nod_size_mg")
length(usedata$Nod_size_mg)
# 150

# check if early-harvested plants were influential
usedata.e <- usedata[which(usedata$Early.harvest=="no"),]
usedata.e %>% nrow
# 139

usedata.o <- usedata[which(usedata$Plant_ID!="1285"),]


####  *** Nitrogen.PC1 model  ####
    
M1 <- lmer(Nod_size_mg ~ Nitrogen.PC1 + Host_line + Nitrogen.PC1:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

drop1(M1, test = "Chisq")
# Df     AIC     LRT Pr(Chi)
# <none>             -180.37                
# Nitrogen.PC1:Host_line  2 -183.40 0.96708  0.6166

M1.noint <- update(M1, .~. -Nitrogen.PC1:Host_line)

drop1(M1.noint, test = "Chisq")
# Df     AIC    LRT Pr(Chi)
# <none>       -183.40               
# Nitrogen.PC1      1 -185.20 0.1989  0.6556
# Host_line  2 -184.09 3.3078  0.1913

# was outlier influential?

M1o <- lmer(Nod_size_mg ~ Nitrogen.PC1 + Host_line + Nitrogen.PC1:Host_line + (1|Block),
           data = usedata.o, REML = FALSE)

# check singularity
model = M1o
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1o
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1o)
plot(simulationOutput)
# good residuals

drop1(M1o, test = "Chisq")
# Df     AIC      LRT Pr(Chi)
# <none>             -208.46                 
# Nitrogen.PC1:Host_line  2 -212.36 0.098903  0.9518

M1o.noint <- update(M1o, .~. -Nitrogen.PC1:Host_line)

drop1(M1o.noint, test = "Chisq")
# Df     AIC     LRT Pr(Chi)
# <none>       -212.36                
# Nitrogen.PC1      1 -214.35 0.00527  0.9421
# Host_line  2 -213.49 2.87187  0.2379

# outlier was not influential

    
####  *** NonNitrogen.PC1 model  ####
    
M1 <- lmer(Nod_size_mg ~ General.PC1 + Host_line + General.PC1:Host_line + (1|Block),
               data = usedata, REML = FALSE)
    
# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

drop1(M1, test = "Chisq")
# Df     AIC    LRT Pr(Chi)  
# <none>             -188.39                 
# General.PC1:Host_line  2 -186.99 5.4036 0.06709 .

M1.noint <- update(M1, .~. -General.PC1:Host_line)

drop1(M1.noint, test = "Chisq")
# Df     AIC    LRT Pr(Chi)  
# <none>       -186.99                 
# General.PC1      1 -185.20 3.7835 0.05176 .
# Host_line  2 -187.59 3.3925 0.18337  


# was outlier influential?
    
M1o <- lmer(Nod_size_mg ~ General.PC1 + Host_line + General.PC1:Host_line + (1|Block),
                data = usedata.o, REML = FALSE)
    
# check singularity
model = M1o
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1o
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1o)
plot(simulationOutput)
# good residuals

drop1(M1o, test = "Chisq")
# Df     AIC    LRT Pr(Chi)  
# <none>             -221.29                 
# General.PC1:Host_line  2 -217.63 7.6593 0.02172 *

# outlier is influential! interaction is now significant

M1o.noint <- update(M1o, .~. -General.PC1:Host_line)

drop1(M1o.noint, test = "Chisq")
# Df     AIC    LRT Pr(Chi)  
# <none>       -217.63                 
# General.PC1      1 -214.35 5.2745 0.02164 *
# Host_line  2 -218.65 2.9764 0.22578  

# plant 1285 had 21 nodules and total nod mass of 21.1 mg; 
# not exceptional for either of those traits:
# nod count was ~ 38 on average, 11-110 nods (common host expt)
# total nod mass was 2.6-26 mg (common host expt)
# I checked harvest photos and nodule count is correct
# also, one nod was very large, so could account for large average nod size
# we will keep this plant in the dataset

    
####  *** NonNitrogen.PC2 model  ####
    

M1 <- lmer(Nod_size_mg ~ General.PC2 + Host_line + General.PC2:Host_line + (1|Block),
           data = usedata, REML = FALSE)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
#FALSE

# will check gradient, should be < 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# good residuals

drop1(M1, test = "Chisq")
# npar     AIC   LRT Pr(Chi)
# <none>                     -184.32              
# General.PC2:Host_line    2 -183.85 4.469   0.107 

M1.noint <- update(M1, .~. -General.PC2:Host_line)

drop1(M1.noint, test = "Chisq")
# npar     AIC    LRT Pr(Chi)
# <none>           -183.85               
# General.PC2    1 -185.20 0.6439  0.4223
# Host_line      2 -184.53 3.3182  0.1903
