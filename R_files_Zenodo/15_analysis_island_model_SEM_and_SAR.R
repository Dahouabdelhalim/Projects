# Renske Onstein - SEMs #
# What explain derived woody species richness on islands? #
#Load libraries
library(lavaan)
library(spdep)
library(akima)
library(raster)
library(ncf)
library(dplyr)
library(writexl)

# set seed
set.seed(2505)

# define function to select highest p-value

identify_weakest <-  function(x){
  rem <- data.frame(lhs = summary(x, stand=T, rsq=T, fit.measures=T)$pe$lhs,
                    rhs = summary(x, stand=T, rsq=T, fit.measures=T)$pe$rhs,
                    pvalue = summary(x, stand=T, rsq=T, fit.measures=T)$pe$pvalue)
  ret <- rem[which(rem$pvalue == max(rem$pvalue, na.rm = TRUE)),]
  
  return(ret)
}

#Load the trait data
#dd <- read.csv("output/final_island_dataset_for_checks.csv")
load(file = "output/data_for_island_model_final.rda")
dd <- iw_model
str(dd)

#Transforming data to obtain (close to) normality
dd$log_GIFT_list_spnum<-log(dd$GIFT_list_spnum)
dd$sqrt_phylowood_iwspnum<-sqrt(dd$phylowood_iwspnum)
dd$log_area<-log(dd$area)
dd$sqrt_dist<-sqrt(dd$dist)
dd$sqrt_mam_total<-sqrt(dd$mam_total)
dd$sqrt_mean_CHELSA_NFD<-sqrt(dd$mean_CHELSA_NFD)
dd$sqrt_mean_ai_yr<-sqrt(dd$mean_ai_yr)
dd$log_mean_b12_v0 <-log(dd$mean_b12_v0 )
dd$log_mean_b1_v0  <-log(dd$mean_b1_v0  )
dd$sqrt_mean_elev<-sqrt(dd$mean_elev+5.1)

#Scaling the continuous data
normalized<-function(y) {
  x<-y[!is.na(y)]
  x<-(x - min(x)) / (max(x) - min(x))
  y[!is.na(y)]<-x
  return(y)
}

str(dd)
#dd2<-apply(dd[c(5:8, 13:36)],2,normalized)
dd2<-apply(dd[c(3:4, 8:28)],2,normalized)
dd2<-as.data.frame(dd2)
dd2$longitude<-dd$longitude
dd2$latitude<-dd$latitude
dd2$geo_entity<-dd$geo_entity
dd2$Arch_name<-dd$Arch_name
dd2$geology_simple<-dd$geology_simple
str(dd2)
dd<-dd2
str(dd)

# Histograms of the variables
hist(dd$sqrt_phylowood_iwspnum)
hist(dd$log_GIFT_list_spnum)
hist(dd$log_area)
hist(dd$sqrt_dist)
hist(dd$sqrt_mam_total)  
hist(dd$sqrt_mean_CHELSA_NFD)
hist(dd$sqrt_mean_ai_yr)
hist(dd$log_mean_b12_v0 )
hist(dd$log_mean_b1_v0)
hist(dd$sqrt_mean_elev )
hist(dd$mean_wc2.0_bio_30s_04)
hist(dd$mean_wc2.0_bio_30s_15)
hist(dd$mean_wc2.0_bio_30s_18)

#------------------ SEM analyses --------------------------------#
# Criteria for fit measures:
# Model Chi-Square with its df and p-value: prefer p-value greater than 0.05
# Root Mean Square Error of Approximation (RMSEA): prefer lower 90%CI to be < 0.05
# Comparative Fit Index (CFI): prefer value greater than 0.90

#The dataset
colnames(dd)
str(dd)
summary(dd)

#Original data saved under name dd3
dd3<-dd

# Rename variables
TotRich <- dd$log_GIFT_list_spnum
Wood <- dd$sqrt_phylowood_iwspnum
Area <- dd$log_area
Dist <- dd$sqrt_dist
Mammal <- dd$sqrt_mam_total
Frost <- dd$sqrt_mean_CHELSA_NFD
Arid <-dd$sqrt_mean_ai_yr
Prec_vel <- dd$log_mean_b12_v0
Temp_vel <- dd$log_mean_b1_v0 
Elev <- dd$sqrt_mean_elev
Season_temp <- dd$mean_wc2.0_bio_30s_04
Season_prec <- dd$mean_wc2.0_bio_30s_15
Prec_wq <- dd$mean_wc2.0_bio_30s_18
Y <- dd$longitude
X <- dd$latitude

# Combine data in new data.frame
Woody.dat <- data.frame(TotRich, Wood, Area, Dist, Mammal, Frost,  Arid, Temp_vel, Prec_vel, Elev, Season_temp, Season_prec, Prec_wq, X, Y)
summary(Woody.dat)
varTable(Woody.dat)

#######################################################
# SEM using derived woody species richness as the response variable #
# Step-wise model selection - remove variables that are least significant in the model, until #
# only significant variables remain #

# SEM approach - Full model
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

# remove prec_wq
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec +  Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

# remove temp_vel
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec +  Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T) 
identify_weakest(sem.fit.1)

# remove prec_vel
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Frost +  Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec +  Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T) 
identify_weakest(sem.fit.1)

# remove frost
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec +  Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T) 
identify_weakest(sem.fit.1)

# remove season_prec
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid +  Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T) 
identify_weakest(sem.fit.1)

# remove arid
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist +  Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T) 
identify_weakest(sem.fit.1)

# remove elev
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Season_temp + TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist +  Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T) 
identify_weakest(sem.fit.1)

# --> FINAL MODEL

out <- data.frame(
  resp = summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$lhs,3,
  pred = summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$rhs,3,
  est_full = round(summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$est,3),
  se_full = round(summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$se,3))

out <- out[1:21,]
out_full <- out %>% 
  mutate(resp = factor(resp, levels = c("Wood", "TotRich", "Mammal"))) %>% 
  group_by(resp) %>% 
  arrange(desc(abs(est_full)), .by_group = TRUE) 

# write_csv(out_full, file = "output/sem_effects_full_model.csv")

#Check modindices and rediduals
modindices(sem.fit.1)
resid(sem.fit.1, type="standardized")

#######################################################
# SEM for all islands excl. Hawaii and Canary

#The dataset
colnames(dd3)
str(dd3)
summary(dd3)

#Exclude Canary and Hawaii
dd2<-subset(dd3, Arch_name != "Canary Islands")
dd2<-subset(dd2, Arch_name != "Hawaiian Islands")
dd<-dd2

# Rename variables
TotRich <- dd$log_GIFT_list_spnum
Wood <- dd$sqrt_phylowood_iwspnum
Area <- dd$log_area
Dist <- dd$sqrt_dist
Mammal <- dd$sqrt_mam_total
Frost <- dd$sqrt_mean_CHELSA_NFD
Arid <-dd$sqrt_mean_ai_yr
Temp_vel <- dd$log_mean_b12_v0
Prec_vel <- dd$log_mean_b1_v0 
Elev <- dd$sqrt_mean_elev
Season_temp <- dd$mean_wc2.0_bio_30s_04
Season_prec <- dd$mean_wc2.0_bio_30s_15
Prec_wq <- dd$mean_wc2.0_bio_30s_18
Y <- dd$longitude
X <- dd$latitude

# Combine data in new dataframe
Woody.dat <- data.frame(TotRich, Wood, Area, Dist, Mammal, Frost,  Arid, Temp_vel, Prec_vel, Elev, Season_temp, Season_prec, Prec_wq, X, Y)
summary(Woody.dat)
varTable(Woody.dat)

# SEM using derived woody species richness as the response - for all islands excl. Hawaii and Canary
# Step-wise model selection - remove variables that are least significant in the model, until #
# only significant variables remain #

# SEM approach - Full model
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    
identify_weakest(sem.fit.1)

# remove prec_wq
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec +Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    
identify_weakest(sem.fit.1)

# remove prec_wq
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_vel + Frost + Temp_vel + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec +Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    
identify_weakest(sem.fit.1)

# remove temp_vel
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_vel + Frost +  Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec +Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    
identify_weakest(sem.fit.1)

# remove elev
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_vel + Frost +  Season_temp + TotRich
TotRich ~ Area +  Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec +Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    
identify_weakest(sem.fit.1)

# remove frost
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_vel +  Season_temp + TotRich
TotRich ~ Area +  Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec +Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    
identify_weakest(sem.fit.1)

# remove season_prec
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_vel +  Season_temp + TotRich
TotRich ~ Area +  Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    
identify_weakest(sem.fit.1)

# remove arid
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_vel +  Season_temp + TotRich
TotRich ~ Area +  Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    
identify_weakest(sem.fit.1)

# remove prec_vel
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec +  Season_temp + TotRich
TotRich ~ Area +  Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)    
identify_weakest(sem.fit.1)

#--> FINAL MODEL
#modindices(sem.fit.1)
#resid(sem.fit.1, type="standardized")

out <- data.frame(
  resp = summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$lhs,3,
  pred = summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$rhs,3,
  est_full_no_can_haw = round(summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$est,3),
  se_full_no_can_haw = round(summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$se,3))


out <- out[1:20,]
out_full_no_can_haw <- out %>% 
  mutate(resp = factor(resp, levels = c("Wood", "TotRich", "Mammal"))) %>% 
  group_by(resp) %>% 
  arrange(desc(abs(est_full_no_can_haw)), .by_group = TRUE) 

results <- out_full %>% 
  full_join(out_full_no_can_haw)

#########################################################################
# SEM only for oceanic islands

#The dataset
colnames(dd3)
str(dd3)
summary(dd3)

#Only oceanic islands
dd2<-subset(dd3,  geology_simple == "oceanic")

dd<-dd2

# Rename variables
TotRich <- dd$log_GIFT_list_spnum
Wood <- dd$sqrt_phylowood_iwspnum
Area <- dd$log_area
Dist <- dd$sqrt_dist
Mammal <- dd$sqrt_mam_total
Frost <- dd$sqrt_mean_CHELSA_NFD
Arid <-dd$sqrt_mean_ai_yr
Temp_vel <- dd$log_mean_b12_v0
Prec_vel <- dd$log_mean_b1_v0 
Elev <- dd$sqrt_mean_elev
Season_temp <- dd$mean_wc2.0_bio_30s_04
Season_prec <- dd$mean_wc2.0_bio_30s_15
Prec_wq <- dd$mean_wc2.0_bio_30s_18
Y <- dd$longitude
X <- dd$latitude

# Combine data in new dataframe
Woody.dat <- data.frame(TotRich, Wood, Area, Dist, Mammal, Frost,  Arid, Temp_vel, Prec_vel, Elev, Season_temp, Season_prec, Prec_wq, X, Y)
summary(Woody.dat)
varTable(Woody.dat)

# SEM using derived woody species richness as the response - only for oceanic islands
# Step-wise model selection - remove variables that are least significant in the model, until #
# only significant variables remain #

# SEM approach - Full model excl. total plant sp. richness
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

# remove temp_vel
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost +Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

# remove arid
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Season_prec + Prec_wq + Prec_vel + Frost +Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

# remove arid
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Season_prec + Prec_wq + Prec_vel + Frost + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)


# remove elev
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Season_prec + Prec_wq + Prec_vel + Frost +  TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

# remove area
WOODY_Model_1 <- 'Wood ~ Elev + Dist + Mammal + Season_prec + Prec_wq + Prec_vel + Frost +  TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

# remove frost
WOODY_Model_1 <- 'Wood ~ Elev + Dist + Mammal + Season_prec + Prec_wq + Prec_vel + Frost +  TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

#remove preq_wq
WOODY_Model_1 <- 'Wood ~ Elev + Dist + Mammal + Season_prec + Prec_wq + Prec_vel + Frost + TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

#remove mammal 
WOODY_Model_1 <- 'Wood ~ Elev + Dist + Season_prec + Prec_wq + Prec_vel + Frost + TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T) 
modindices(sem.fit.1)

# add covariate --> FINAL MODEL
WOODY_Model_1 <- 'Wood ~ Elev + Dist +Season_prec + Prec_wq + Prec_vel + Frost +  TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~  Elev + Dist + Arid + Season_prec + Season_temp
TotRich ~~      Mammal
'
sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  

# --> FINAL MODEL

out <- data.frame(
  resp = summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$lhs,3,
  pred = summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$rhs,3,
  est_oceanic = round(summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$est,3),
  se_oceanic = round(summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$se,3))


out <- out[1:20,]
out_oceanic <- out %>% 
  mutate(resp = factor(resp, levels = c("Wood", "TotRich", "Mammal"))) %>% 
  group_by(resp) %>% 
  arrange(desc(abs(est_oceanic)), .by_group = TRUE) 

results <- results %>% 
  full_join(out_oceanic)

#modindices(sem.fit.1)
#resid(sem.fit.1, type="standardized")

#######################################################
# SEM only for oceanic islands excl. Hawaii and Canary

#The dataset
colnames(dd)
str(dd)
summary(dd)

#Exclude Canary and Hawaii
dd2<-subset(dd, Arch_name != "Canary Islands")
dd2<-subset(dd2, Arch_name != "Hawaiian Islands")
dd<-dd2

# Rename variables
TotRich <- dd$log_GIFT_list_spnum
Wood <- dd$sqrt_phylowood_iwspnum
Area <- dd$log_area
Dist <- dd$sqrt_dist
Mammal <- dd$sqrt_mam_total
Frost <- dd$sqrt_mean_CHELSA_NFD
Arid <-dd$sqrt_mean_ai_yr
Temp_vel <- dd$log_mean_b12_v0
Prec_vel <- dd$log_mean_b1_v0 
Elev <- dd$sqrt_mean_elev
Season_temp <- dd$mean_wc2.0_bio_30s_04
Season_prec <- dd$mean_wc2.0_bio_30s_15
Prec_wq <- dd$mean_wc2.0_bio_30s_18
Y <- dd$longitude
X <- dd$latitude

# Combine data in new dataframe
Woody.dat <- data.frame(TotRich, Wood, Area, Dist, Mammal, Frost,  Arid, Temp_vel, Prec_vel, Elev, Season_temp, Season_prec, Prec_wq, X, Y)
summary(Woody.dat)
varTable(Woody.dat)

# SEM using derived woody species richness as the response - only for oceanic islands excl. Hawaii and Canary

# SEM approach - Full model
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel + Season_temp + TotRich
TotRich ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

# Remove elev
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel + Season_temp + TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)   
identify_weakest(sem.fit.1)

# Remove season_temp
WOODY_Model_1 <- 'Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel +TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

# Remove area
WOODY_Model_1 <- 'Wood ~ Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel +TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Area + Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)   
identify_weakest(sem.fit.1)

# Remove area
WOODY_Model_1 <- 'Wood ~ Elev + Dist + Mammal + Arid + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel +TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T) 
identify_weakest(sem.fit.1)

# Remove arid
WOODY_Model_1 <- 'Wood ~ Elev + Dist + Mammal + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel +TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Elev + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

#remove frost
WOODY_Model_1 <- 'Wood ~ Elev + Dist + Mammal + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel +TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Elev + Dist + Arid + Season_prec + Prec_wq + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

#remove preq_wq
WOODY_Model_1 <- 'Wood ~ Elev + Dist + Mammal + Season_prec + Prec_wq + Prec_vel + Frost + Temp_vel +TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Elev + Dist + Arid + Season_prec + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

#remove temp_vel --> FINAL MODEL
WOODY_Model_1 <- 'Wood ~ Elev + Dist + Mammal + Season_prec + Prec_wq + Prec_vel + Frost + TotRich
TotRich ~ Area + Dist + Arid + Season_prec + Prec_wq + Frost + Season_temp
Mammal ~ Elev + Dist + Arid + Season_prec + Season_temp
'

sem.fit.1 <- sem(WOODY_Model_1, data=Woody.dat)
summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)  
identify_weakest(sem.fit.1)

#modindices(sem.fit.1)
#resid(sem.fit.1, type="standardized")

out <- data.frame(
  resp = summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$lhs,3,
  pred = summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$rhs,3,
  est_oceanic_no_haw_no_can = round(summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$est,3),
  se_oceanic_no_haw_no_can = round(summary(sem.fit.1, stand=T, rsq=T, fit.measures=T)$pe$se,3))


out <- out[1:20,]
out_oceanic_no_haw_no_can <- out %>% 
  mutate(resp = factor(resp, levels = c("Wood", "TotRich", "Mammal"))) %>% 
  group_by(resp) %>% 
  arrange(desc(abs(est_oceanic_no_haw_no_can)), .by_group = TRUE) 

results <- results %>% 
  full_join(out_oceanic_no_haw_no_can) %>% 
  dplyr::select(-contains("X")) %>% 
  arrange(resp, desc(abs(est_full)))

write_xlsx(results, path = "output/sem_effect_estimates.xlsx")


#####################    SPATIAL    ########################################

dd<-dd3

# Rename variables
TotRich <- dd$log_GIFT_list_spnum
Wood <- dd$sqrt_phylowood_iwspnum
Area <- dd$log_area
Dist <- dd$sqrt_dist
Mammal <- dd$sqrt_mam_total
Frost <- dd$sqrt_mean_CHELSA_NFD
Arid <-dd$sqrt_mean_ai_yr
Temp_vel <- dd$log_mean_b12_v0
Prec_vel <- dd$log_mean_b1_v0 
Elev <- dd$sqrt_mean_elev
Season_temp <- dd$mean_wc2.0_bio_30s_04
Season_prec <- dd$mean_wc2.0_bio_30s_15
Prec_wq <- dd$mean_wc2.0_bio_30s_18
Y <- dd$longitude
X <- dd$latitude

# Combine data in new dataframe
Woody.dat <- data.frame(TotRich, Wood, Area, Dist, Mammal, Frost,  Arid, Temp_vel, Prec_vel, Elev, Season_temp, Season_prec, Prec_wq, X, Y)
summary(Woody.dat)
varTable(Woody.dat)

# Examine the structure of the dataset
wood<-Woody.dat
names(wood)
str(wood)

# Make a OLS regression model
lm_wood<-lm(Wood ~ Area + Elev + Dist + Mammal + Arid + Season_prec + Prec_wq +  Season_temp + TotRich)
summary(lm_wood)
hist(residuals(lm_wood))

par(mfrow=c(3,3))
termplot(lm_wood, partial=TRUE, terms=1, pch=16)
termplot(lm_wood, partial=TRUE, terms=2, pch=16)
termplot(lm_wood, partial=TRUE, terms=3, pch=16)
termplot(lm_wood, partial=TRUE, terms=4, pch=16)
termplot(lm_wood, partial=TRUE, terms=5, pch=16)
termplot(lm_wood, partial=TRUE, terms=6, pch=16)
termplot(lm_wood, partial=TRUE, terms=7, pch=16)
termplot(lm_wood, partial=TRUE, terms=8, pch=16)
termplot(lm_wood, partial=TRUE, terms=9, pch=16)

#Spatial structure of residuals

# Correlograms with latlon = FALSE
cor.OBL<-correlog(wood$X, wood$Y, z=wood$Wood, na.rm=T, increment=1, resamp=1, latlon = FALSE)
cor.OBL
plot(cor.OBL$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)  

#With latlon = FALSE and increment=52.9 (minimum distance between cells), resampling 999 times
cor.OBL_1000<-correlog(wood$X, wood$Y, z=wood$Wood, na.rm=T, increment=52.9, resamp=999, latlon = FALSE)  #uses km because latlon = TRUE
plot(cor.OBL_1000$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0) 
summary(cor.OBL_1000)

#Correlogram for residuals
cor.res_1000<-correlog(wood$X, wood$Y, z=residuals(lm_wood), na.rm=T, increment=52.9, resamp=999, latlon = FALSE)

#Plot both residuals and raw data
plot(cor.OBL_1000$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)                  
points(cor.res_1000$correlation, pch=1, cex=1.2)
lines(cor.res_1000$correlation, lwd=1.5)

#Make coordinate list
coords_wood<-as.matrix(cbind(wood$X,wood$Y))
plot(coords_wood)

#Minimum distance to connect to at least one neighbor
wood_knear <- knn2nb(knearneigh(coords_wood, k=1))
summary(wood_knear)
dsts_wood<-unlist(nbdists(wood_knear, coords_wood, longlat = FALSE))
summary(dsts_wood)
max(dsts_wood)

#Calculate Moran's I values:
wood_nb_moran1000<-dnearneigh(coords_wood,0,52.9, longlat=F)
nb_moran1000<-nb2listw(wood_nb_moran1000, glist=NULL, style="W", zero.policy=TRUE)

# For raw data
moran.test(wood$Wood, listw=nb_moran1000, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")

# For residuals of OLS model
moran.test(residuals(lm_wood), listw=nb_moran1000, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")
