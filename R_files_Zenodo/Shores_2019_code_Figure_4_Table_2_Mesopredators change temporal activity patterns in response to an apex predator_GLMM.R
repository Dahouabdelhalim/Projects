#Load packages
library(mgcv)
library(bbmle)
library(multcomp)
library(glmmTMB)

#Set working directory 
setwd("C:/file path")


#Import data
data1 <- read.csv("C:/file directory.csv")

#Data specific to certain study areas
aeneas_only <- data1[ which(data1$study_area=='A'),]
bonaparte_only <- data1[ which(data1$study_area=='B'),]
strawberry_only <- data1[ which(data1$study_area=='S'),]
ncicn_only <- data1[ which(data1$study_area=='N'),]



#####Human models
human_null <- glmmTMB(human_activity ~ 
                        1, data = data1, family = poisson, ziformula = ~ 1)
#Human location model only 
#(1|location) is a random effect for camera location
#ziformula is zero-inflated formula. It can be general or set to specific parameters
#family = error distribution family
human_location <- glmmTMB(human_activity ~ 
                            (1|location), data = data1, family = poisson, ziformula = ~ 1)

#"offset" accounts for "exposure time" for different time categories
#Dusk and dawn are each 2 hours, day and night are 10 hours
#Used offset indicating fraction of a day for each time period 
human_location__time_season <- glmmTMB(human_activity ~ 
                                         (1|location) + 
                                         time_cat + offset(log(day)) + 
                                         season, data = data1, family = poisson, ziformula = ~ 1)
#Get model summary
summary(human_location)

#######Wolf models
#Ran models relevant to study hypotheses, not every possible combination of parameters
wolf_null <- glmmTMB(wolf_activity ~ 
                       1, data = data1, family = poisson, ziformula = ~ 1)

wolf_location <- glmmTMB(wolf_activity ~ 
                           (1|location), data = data1, family = poisson, ziformula = ~ 1)
wolf_location__time_season <- glmmTMB(wolf_activity ~ 
                                        time_cat + offset(log(day)) + 
                                        season + 
                                        (1|location), data = data1, family = poisson, ziformula = ~ 1)
wolf_location__time_season_hum <- glmmTMB(wolf_activity ~ 
                                            (1|location) + 
                                            time_cat + offset(log(day)) + 
                                            season + 
                                            human_activity, data = data1, family = poisson, ziformula = ~ 1)

#########Coyote models
coyote_null <- glmmTMB(coyote_activity ~ 
                         1, data = data1, family = poisson, ziformula = ~ 1)
coyote_location <- glmmTMB(coyote_activity ~ 
                             (1|location), data = data1, family = poisson, ziformula = ~ 1)
coyote_loc_sea_time <- glmmTMB(coyote_activity ~ 
                                 (1|location) + 
                                 season + 
                                 time_cat + offset(log(day)), data = data1, family = poisson, ziformula = ~ 1)

coyote_loc_sea_time_human <- glmmTMB(coyote_activity ~ 
                                 human_activity + 
                                 (1|location) + 
                                 season + 
                                 time_cat + offset(log(day)), data = data1, family = poisson, ziformula = ~ 1)

#The wolf_presence*time_cat category tests for differences in coyote activity between study areas with and without wolves
coyote_loc_sea_time_intrct_wolfpresence_wolf<- glmmTMB(coyote_activity ~ 
                                                                  (1|location) + 
                                                                  season + 
                                                                  wolf_activity +
                                                                  wolf_presence*time_cat +  + offset(log(day)), data = data1, family = poisson, ziformula = ~ 1)
coyote_loc_sea_time_intrct_wolfpresence_human_wolf<- glmmTMB(coyote_activity ~ 
                                                                             time_cat * wolf_presence + offset(log(day)) + 
                                                                             (1|location) + 
                                                                             season + 
                                                                             human_activity + 
                                                                             wolf_activity, data = data1, family = poisson, ziformula = ~ 1)


#####Bobcat models
bobcat_null <- glmmTMB(bobcat_activity ~ 
                         1, data = data1, family = poisson, ziformula = ~ 1)
bobcat_location <- glmmTMB(bobcat_activity ~ 
                             (1|location), data = data1, family = poisson, ziformula = ~ 1)

bobcat_loc_sea_time <- glmmTMB(bobcat_activity ~ 
                                 (1|location) + 
                                 season + 
                                 time_cat + offset(log(day)), family = poisson, ziformula = ~ 1, data = data1)
bobcat_hum_loc_sea_time <- glmmTMB(bobcat_activity ~ 
                                     (1|location) + 
                                     season + 
                                     time_cat+offset(log(day)) + 
                                     human_activity, family = poisson, ziformula = ~ 1, data = data1)
bobcat_coy_loc_sea_time <- glmmTMB(bobcat_activity ~ 
                                     (1|location) + 
                                     season + 
                                     time_cat + offset(log(day)) + 
                                     coyote_activity, family = poisson, ziformula = ~ 1, data = data1)
bobcat_loc_sea_time_humIntcoy <- glmmTMB((bobcat_activity ~ 
                                            (1|location) + 
                                            season + 
                                            time_cat + offset(log(day)) + 
                                            human_activity + 
                                            coyote_activity), family = poisson, ziformula = ~ 1, data = data1)
bobcat_loc_sea_wolfInttime_humIntcoy <- glmmTMB((bobcat_activity ~ 
                                                   (1|location) + 
                                                   season + 
                                                   wolf_presence*time_cat + offset(log(day)) + 
                                                   wolf_activity + 
                                                   human_activity*coyote_activity), family = poisson, ziformula = ~ 1, data = data1)
bobcat_loc_sea_wolfInttime_wolfIntcoy <- glmmTMB((bobcat_activity ~ 
                                                    (1|location) + 
                                                    season + 
                                                    wolf_presence*time_cat + 
                                                    offset(log(day)) + 
                                                    human_activity + 
                                                    wolf_activity + 
                                                    coyote_activity), family = poisson, ziformula = ~ 1, data = data1)
bobcat_loc_sea_wolfInttime_wolf <- glmmTMB((bobcat_activity ~ 
                                                    (1|location) + 
                                                    season + 
                                                    wolf_presence*time_cat + 
                                                    offset(log(day)) + 
                                                    wolf_activity), family = poisson, ziformula = ~ 1, data = data1)

bobcat_loc_sea_wolfInttime_coy_hum_wolf <- glmmTMB((bobcat_activity ~ 
                                                      (1|location) + 
                                                      season + 
                                                      wolf_presence*time_cat + offset(log(day)) + 
                                                      human_activity + 
                                                      wolf_activity + 
                                                      coyote_activity), family = poisson, ziformula = ~ 1, data = data1)

