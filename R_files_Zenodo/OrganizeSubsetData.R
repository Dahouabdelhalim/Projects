

#################################################################################################################
## PAPER: Occupancy winners in tropical protected forests: a pantropical analysis (Semper-Pascual et al. 2022) ##
#################################################################################################################

### ---- Organize and subset covariates ---- ###

library(zoo)
library(tidyverse)
library(reshape)
library(corrplot)

## Set working directory
setwd("C:/Users/SemperPascual_et_al_2022/")

## Load TEAM information (number of protected areas, years, species...)
load("InfoData_SemperPascual_et_al_2022.RData")
data.list <- data.list2

# File containing species list
species_table <- read.csv("SpeciesList_SemperPascual_et_al_2022.csv", sep = ",")

## Dimensions used in the model
# This is extracted automatically, but here I manually indicate the dimensions)
n.cont <- 3
n.areas <- c(8,5,3)
n.sites <- matrix(c(60, 60, 60, 65, 89, 60, 63, 60, 60, 60, 64, 61, 60, 0, 0, 0, 60, 60, 60, 0, 0, 0, 0, 0),
                  ncol=3, nrow=8, byrow = F)


################################################
##########       Read covariates      ##########
################################################

# -- SPECIES covariates ---
sp_cov <- read.csv("SpeciesCovariates.csv")
head(sp_cov)

# Covariate for detection: Forest strata
levels(as.factor(sp_cov$ForStrat))  # 3 classes
# Covariate for occupancy: Trophic guild
levels(as.factor(sp_cov$T.Guild))  # 4 classes

# Check which functional group is more abundant and then use it as reference level (intercept) 
table(sp_cov$T.Guild)    # Herbivores is the most common group
table(sp_cov$ForStrat)   # Ground is the most common group

# Rename with columns
names(sp_cov)[3] <- "ForStrat.tmp"

sp_cov <- sp_cov %>% 
  mutate(ForStrat = ifelse(sp_cov$ForStrat.tmp == "G", 0, 1),              # Simplify to 2 classes ("G" --> 0; "Ar" and "S" --> 1). I give G 0 (most abundant group) so it becomes the reference value (intercepts)
         Carn = ifelse(sp_cov$T.Guild == "Carnivorous", 1, 0),      # Create dummy variables for Trophic guild.
         Insec = ifelse(sp_cov$T.Guild == "Insectivorous", 1, 0),   # There are 4 trophic guilds so I create 3 dummy variables (Herbivorous as reference because its the most common)
         Omni = ifelse(sp_cov$T.Guild == "Omnivorous", 1, 0))  

## Merge species_table with the covariates (ForStrat and Hab.Breadth)
species_table_cov <- species_table %>% 
  left_join(sp_cov, by ="Species.Name", keep) %>% 
  mutate(logHab.Breadth_unSd = log(Hab.Breadth),
         logMass_unSd = log(Mass.g)) %>% 
  mutate(logHab.Breadth = scale(logHab.Breadth_unSd), # Standardize
         logMass = scale(logMass_unSd)) %>% 
  select(c("Species.Name","Species.Name2","logMass","logHab.Breadth","ForStrat","T.Guild","Carn","Insec","Omni"))

# Check for NAs
subset(species_table_cov,is.na(species_table_cov$logHab.Breadth))
subset(species_table_cov,is.na(species_table_cov$logMass))

species_table_cov2 <-species_table_cov[order(species_table_cov$Species.Name2),] # order by Species.Name2

# Select covariates for detection and occupancy
names(species_table_cov2)
covs.p <- species_table_cov2[,c("ForStrat","logMass")]
covs.psi <- species_table_cov2[,c("logMass","logHab.Breadth","Carn", "Insec", "Omni")]
head(covs.psi)

# -- SPATIAL covariates ---

buffer <- 10000

# a) Perc. forest and fragmentation
percForest_Large_cov <- read.csv(paste0("perc_forest_buff_",buffer,"m.csv"))
division_cov <- read.csv(paste0("division_buff_",buffer,"m.csv"))
# b) Human population
pop_cov <- read.csv(paste0("pop2015_buff_",buffer,"m.csv"))

spatial_cov <- do.call("cbind",list(percForest_Large_cov,
                                    division_cov,
                                    pop_cov))

head(spatial_cov)
spatial_cov <- spatial_cov[,c(1,2,3,4,8,12)]

# Remove Madagascar
spatial_cov <- spatial_cov %>% 
  filter(!(Site.Code == "RNF")) %>% 
  # and remove the camera-trap site from Manaus that only worked 6 days and detected birds
  filter(!(Sampling.Unit.Name == "CT-MAS-3-06"))

# Histograms
par(mfrow = c(1,2))
hist(log(spatial_cov$pop2015))
hist(spatial_cov$percForest_Large_mean)

# Check correlation between landscape variables
names(spatial_cov)
par(mfrow=c(1,1))
corrplot(cor(spatial_cov[,c(4:6)]), tl.col = "red", bg = "White", tl.srt = 35, 
         addCoef.col = "black") # Forest and division are correlated !!!

names(spatial_cov)
spatial_cov <- spatial_cov %>% 
  mutate(logPop2015_unSd = log(pop2015+1)) %>% # Log transform
  mutate(percForest_Large_mean = scale(percForest_Large_mean), # Standardize
         division_mean = scale(division_mean),
         logPop2015 = scale(logPop2015_unSd))

head(spatial_cov)

# Create two arrays, one for each landscape variable
# Three dimensions array (sites*pa*c)

# -- 1) Percentage of forest -- #
names(spatial_cov)
percForest_df <- spatial_cov[,c(1,3,4)] 
head(percForest_df)

PAs <- unique(percForest_df$Site.Code)

percForest <- array(NA, dim = c(max(n.sites),max(n.areas),n.cont))
for (c in 1:n.cont){
  ifelse(c == 1, Cont2 <- "Ame",ifelse(c == 2, Cont2 <- "Afr", Cont2 <- "Asia"))
  percForest_oneCont <- array(NA, dim = c(max(n.sites),max(n.areas)))
  for(pa in 1:n.areas[c]){
    PA <- data.list$paID[which(data.list$paID$Cont == Cont2),]$Site.Code[pa] # Name of protected area "pa" and continent "c"
    percForest_onePA <- subset(percForest_df, Site.Code == PA)
    percForest_oneCont[1:nrow(percForest_onePA),pa] <- percForest_onePA$percForest_Large_mean
    
  }
  percForest[,,c] <- percForest_oneCont
}

# -- 2) Division index -- #
names(spatial_cov)
Division_df <- spatial_cov[,c(1,3,5)] # Select one buffer (e.g. 10000 meters)
head(Division_df)

PAs <- unique(Division_df$Site.Code)

Division <- array(NA, dim = c(max(n.sites),max(n.areas),n.cont))
for (c in 1:n.cont){
  ifelse(c == 1, Cont2 <- "Ame",ifelse(c == 2, Cont2 <- "Afr", Cont2 <- "Asia"))
  Division_oneCont <- array(NA, dim = c(max(n.sites),max(n.areas)))
  for(pa in 1:n.areas[c]){
    PA <- data.list$paID[which(data.list$paID$Cont == Cont2),]$Site.Code[pa] # Name of protected area "pa" and continent "c"
    Division_onePA <- subset(Division_df, Site.Code == PA)
    Division_oneCont[1:nrow(Division_onePA),pa] <- Division_onePA$division_mean
    
  }
  Division[,,c] <- Division_oneCont
}

# -- 3) Population -- #
names(spatial_cov)
Pop_df <- spatial_cov[,c(1,3,8)] # Select one buffer (e.g. 10000 meters)
head(Pop_df)

PAs <- unique(Pop_df$Site.Code)

logPop <- array(NA, dim = c(max(n.sites),max(n.areas),n.cont))
for (c in 1:n.cont){
  ifelse(c == 1, Cont2 <- "Ame",ifelse(c == 2, Cont2 <- "Afr", Cont2 <- "Asia"))
  Pop_oneCont <- array(NA, dim = c(max(n.sites),max(n.areas)))
  for(pa in 1:n.areas[c]){
    PA <- data.list$paID[which(data.list$paID$Cont == Cont2),]$Site.Code[pa] # Name of protected area "pa" and continent "c"
    Pop_onePA <- subset(Pop_df, Site.Code == PA)
    Pop_oneCont[1:nrow(Pop_onePA),pa] <- Pop_onePA$logPop2015
    
  }
  logPop[,,c] <- Pop_oneCont
}



