################################################################################
#                                                                              #
#                                                                              #
#                        STATISTICAL ANALYSIS                                  #
#                                                                              #
#                                                                              #
################################################################################

#Load (and install) required packages

### Read packages 
library(sp)
library(rgdal)
library(dplyr)
library(lubridate)
library(amt)
library(tidyr)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(purrr)
library(pbapply)
library(reshape2)
library(grid)
library(RColorBrewer)
library(forcats)
library(lme4)
library(nlme)
library(piecewiseSEM)
library(ggResidpanel)
library(ggeffects)
library(MuMIn)
require(glmmTMB)
require(MuMIn)
require(afex)
require(sf)
require(raster)


################################################################################
#
#                      Home range area (KDE 95%)
#
################################################################################

##--load our df to R environment
hr_ind<-get(load("D:/Escritorio/GRIFFON_VULTURE_HOME_RANGE/DATA_ANALYSES/GV_DISTANCE_HR/Gv_hr_year.RData"))

#inspect the first rows of the df
head(hr_ind)

#remove XCM individual since it belongs to Gyps ruppell (not the focus species)
hr_ind <- hr_ind[!(hr_ind$id %in% c("XCM")),]

#load our population metadata to join it later to our df (always with UTF-8 or ISO encoding)
setwd("D:/Escritorio/GRIFFON_VULTURE_HOME_RANGE/DATA_ANALYSES")
pops<-read.csv("POPS.csv",header=T,dec=".",sep=";",encoding = "ISO-8859-1")

#rename some factor levels to avoid discordancies when joining
pops$id<-recode_factor(pops$id, 
                       "Alo\\xf1a"="Aloña", 
                       "O\\xf1ati"="Oñati") 

#join both df to keep the population identity for each bird
hr_ind_def<-left_join(hr_ind,pops,by="id")

#get only data for the adults
hr_ind_clean<-hr_ind_def %>% filter( age.x %in% c("adult"))

#recode levels (our interest is only in KDE 95)
hr_ind_clean$level<-recode_factor(hr_ind_clean$level, 
                                  "0.95"="KDE 95", 
                                  "0.75"="KDE 75",
                                  "0.5"="KDE 50"
) 

#recode factor levels again, this time for the season
hr_ind_clean$month_name<-recode_factor(hr_ind_clean$month, 
                                       "1"="January","2"="February","3"="March",
                                       "4"="April","5"="May","6"="June","7"="July","8"="August","9"="September","10"="October","11"="November","12"="December")

#create our season variable
hr_ind_clean <- hr_ind_clean %>%
  mutate(season = case_when(
    month %in% c("7","8","9") ~ "Summer",
    month %in% c("10","11","12") ~ "Autumn",
    month %in% c("1","2","3") ~ "Winter",
    month %in% c("4","5","6") ~ "Spring"))

#remove the individual from castellon, since it is the only representant for the population
hr_ind_clean <- hr_ind_clean[!(hr_ind_clean$population %in% c("Castellon")),]

#recode population names
hr_ind_clean$population<-recode_factor(hr_ind_clean$population, 
                                                "Bardenas"="Alto Ebro", 
                                                "Burgos"="Segovia",
                                                "Cadiz"="C?diz",
                                                "Castellon"="Castellon",
                                                "Cazorla"="Cazorla",
                                                "Guipuzkoa"="Alto Ebro",
                                                "Pirineo Central"="Pyrenees",
                                                "Pirineo Oriental"="Pyrenees",
                                                "Soria"="Segovia") 

#Extract KDE 95 
hr_ind_clean_95<-hr_ind_clean %>% filter( level %in% c("KDE 95"))

##-- Compute some statistics 
hr_ind_clean_95 %>% group_by(season) %>% summarise(mean=mean(area_km2),sd=sd(area_km2),range=range(area_km2))

##--Modelling

#This is for model dredging with MUMin package
options(na.action = "na.fail")

#run our general model by using glmmTMB. We will use the gaussian family in this case and assume normal distribution of the data
#Year and id are included as random terms to control for the pseudoreplication  
general_mod<-glmmTMB(area_km2~factor(sex)+factor(population)+factor(season)+(1|year/id),data=hr_ind_clean_95, ziformula=~1)

#model dredging
dredge(general_mod)

#The best model is that including all variables (the saturated model). 
best_mod<-glmmTMB(area_km2~factor(sex)+factor(population)+factor(season)+(1|year/id),data=hr_ind_clean_95, ziformula=~1)

#get summary statistics for the best model
summary(best_mod)

#get the R2 variance by using performance package and r2_nakagawa function
performance::r2_nakagawa(best_mod)

#--get the marginal means by using emmeans package for each factor

#sex
emmeans::emmeans(best_mod, specs = pairwise ~ factor(sex), type = "response")

#population
emmeans::emmeans(best_mod, specs = pairwise ~ factor(population), type = "response")

#season
emmeans::emmeans(best_mod, specs = pairwise ~ factor(season), type = "response")



################################################################################
#
#                      Cumulative distance
#
################################################################################


##--We do the same for travelled distance
dist_ind<-get(load("D:/Escritorio/GRIFFON_VULTURE_HOME_RANGE/DATA_ANALYSES/GV_DISTANCE_HR/Gv_distance_cumulated.RData"))

#Check the first lines of our cumulative distance dataframe
head(dist_ind)

#Rename id column
colnames(dist_ind)[1]<-"id"

#remove XCM individual since it belongs to Gyps ruppell (not the focus species)
dist_ind <- dist_ind[!(dist_ind$id %in% c("XCM")),]

#load our population metadata to join it later to our df (always with UTF-8 or ISO encoding)
setwd("D:/Escritorio/GRIFFON_VULTURE_HOME_RANGE/DATA_ANALYSES")
pops<-read.csv2("POPS.csv",header=T,dec=".",sep=";")

#join both df to keep the population identity for each bird
dist_ind_def<-left_join(dist_ind,pops,by="id")

#recode factor levels again, this time for the season
dist_ind_def$month_name<-recode_factor(dist_ind_def$month, 
                                       "1"="January","2"="February","3"="March",
                                       "4"="April","5"="May","6"="June","7"="July","8"="August","9"="September","10"="October","11"="November","12"="December")

#recode population names
dist_ind_def$population<-recode_factor(dist_ind_def$population, 
                                       "Bardenas"="Alto Ebro", 
                                       "Burgos"="Segovia",
                                       "Cadiz"="C?diz",
                                       "Castellon"="Castellon",
                                       "Cazorla"="Cazorla",
                                       "Guipuzkoa"="Alto Ebro",
                                       "Pirineo Central"="Pyrenees",
                                       "Pirineo Oriental"="Pyrenees",
                                       "Soria"="Segovia") 

#create our season variable
dist_ind_def <- dist_ind_def %>%
  mutate(season = case_when(
    month %in% c("7","8","9") ~ "Summer",
    month %in% c("10","11","12") ~ "Autumn",
    month %in% c("1","2","3") ~ "Winter",
    month %in% c("4","5","6") ~ "Spring"))

#remove the individual from castellon, since it is the only representant for the population
dist_ind_def <- dist_ind_def[!(dist_ind_def$population %in% c("Castellon")),]

##--Modelling

#This is for model dredging with MUMin package
options(na.action = "na.fail")

#run our general model by using glmmTMB. We will use the gaussian family in this case and assume normal distribution of the data
#Year and id are included as random terms to control for the pseudoreplication  
general_mod<-glmmTMB(cum_dist~factor(sex)+factor(population)+factor(season)+(1|year/id),data=dist_ind_def, ziformula=~1)

#model dredging
dredge(general_mod)

#The best model is that including all variables (the saturated model). 
best_mod<-glmmTMB(cum_dist~factor(sex)+factor(population)+factor(season)+(1|year/id),data=dist_ind_def, ziformula=~1)

#get summary statistics for the best model
summary(best_mod)

#get the R2 variance by using performance package and r2_nakagawa function
performance::r2_nakagawa(best_mod)

#--get the marginal means by using emmeans package for each factor

#sex
emmeans::emmeans(best_mod, specs = pairwise ~ factor(sex), type = "response")

#population
emmeans::emmeans(best_mod, specs = pairwise ~ factor(population), type = "response")

#season
emmeans::emmeans(best_mod, specs = pairwise ~ factor(season), type = "response")


################################################################################
#
#                      Individual overlap in home range area
#
################################################################################

##--We do the same for home range overlap in home range area
Gv_ind_over<-get(load("D:/Escritorio/GRIFFON_VULTURE_HOME_RANGE/DATA_ANALYSES/GV_DISTANCE_HR/Gv_overlap_ind.RData"))

#Check the first lines of our home range overlap dataframe
head(Gv_ind_over)

#get rid of the non necessary rows by selecting only the overalp values for each id and time period (months)
Gv_overlap_ind_def<-Gv_ind_over %>% filter(variable==c("overlap"))

# Make timestamp a date/time variable
Gv_overlap_ind_def$from <-as.POSIXct(Gv_overlap_ind_def$from, format="%Y-%m-%d", tz="UTC")
Gv_overlap_ind_def$to <-as.POSIXct(Gv_overlap_ind_def$to, format="%Y-%m-%d", tz="UTC")


Gv_overlap_ind_def$month_1<-strftime(Gv_overlap_ind_def$from, "%B")
Gv_overlap_ind_def$month_2<-strftime(Gv_overlap_ind_def$to, "%B")

Gv_overlap_ind_def$period<-paste(Gv_overlap_ind_def$month_1,Gv_overlap_ind_def$month_2, sep = "_")

#rename the period names 
Gv_indv_over <- Gv_overlap_ind_def %>%
  mutate(season = case_when(
    period %in% c("diciembre_enero","enero_febrero","febrero_marzo","noviembre_enero","enero_marzo","diciembre_febrero","diciembre_marzo") ~ "Winter",
    period %in% c("marzo_abril","abril_mayo","mayo_junio","marzo_mayo","marzo_junio") ~ "Spring",
    period %in% c("junio_julio","julio_agosto","agosto_septiembre","junio_septiembre") ~ "Summer",
    period %in% c("septiembre_octubre","octubre_noviembre","noviembre_diciembre","octubre_diciembre","septiembre_noviembre") ~ "Autumn"))

#extract the year from our bimonth period
Gv_indv_over$year<-year(Gv_indv_over$from)

#recode the period variable
Gv_indv_over_def<-Gv_indv_over %>%
  mutate(period = fct_relevel(period, 
                              "enero_febrero","febrero_marzo","marzo_abril","abril_mayo","mayo_junio","junio_julio","julio_agosto",
                              "agosto_septiembre","septiembre_octubre","octubre_noviembre","noviembre_diciembre","diciembre_enero")) 

#rename the recoded variable
Gv_indv_over_def$period<-recode_factor(Gv_indv_over_def$period, 
                                       "enero_febrero"="January-February","febrero_marzo"="February-March","marzo_abril"="March-April","abril_mayo"="April-May","mayo_junio"="May-June","junio_julio"="June-July","julio_agosto"="July-August",
                                       "agosto_septiembre"="August-September","septiembre_octubre"="September-October","octubre_noviembre"="October-November","noviembre_diciembre"="November-December","diciembre_enero"="December-January") 

#rename the id column
colnames(Gv_indv_over_def)[5]<-"id"

#load our population metadata to join it later to our df (always with UTF-8 or ISO encoding)
setwd("D:/Escritorio/GRIFFON_VULTURE_HOME_RANGE/DATA_ANALYSES")
pops<-read.csv2("POPS.csv",header=T,dec=".",sep=";")

#join both df to keep the population identity for each bird
Gv_indv_over_clean_def<-left_join(Gv_indv_over_def,pops,by="id")

Gv_indv_over_clean_def<-na.omit(Gv_indv_over_clean_def)

Gv_indv_over_clean_def <- Gv_indv_over_clean_def[!(Gv_indv_over_clean_def$population %in% c("Castellon")),]

Gv_indv_over_clean_def$population<-recode_factor(Gv_indv_over_clean_def$population, 
                                                 "Guipuzkoa"="Alto Ebro", 
                                                 "Bardenas"="Alto Ebro", 
                                                 "Burgos"="Segovia",
                                                 "Cadiz"="C?diz",
                                                 "Cazorla"="Cazorla",
                                                 "Pirineo Central"="Pyrenees",
                                                 "Pirineo Oriental"="Pyrenees",
                                                 "Soria"="Segovia") 


##--some statistics
Gv_indv_over_clean_def %>% group_by(season,sex) %>% summarise(mean=mean(value),sd=sd(value))


#This is for model dredging with MUMin package
options(na.action = "na.fail")

#run our general model by using glmmTMB. We will use the beta family in this case since they are proportion data (from 0 to 1)
#Year and id are included as random terms to control for the pseudoreplication  
general_mod<-glmmTMB(value~factor(sex)+factor(population)+factor(season)+(1|year/id),data=dist_ind_def, ziformula=~1,family=beta_family(link="logit"))

#model dredging
dredge(general_mod)

#The best model is that including all variables (the saturated model). 
best_mod<-glmmTMB(value~factor(sex)+factor(population)+factor(season)+(1|year/id),data=dist_ind_def, ziformula=~1,family=beta_family(link="logit"))

#get summary statistics for the best model
summary(best_mod)

#get the R2 variance by using performance package and r2_nakagawa function
performance::r2_nakagawa(best_mod)

#--get the marginal means by using emmeans package for each factor

#sex
emmeans::emmeans(best_mod, specs = pairwise ~ factor(sex), type = "response")

#population
emmeans::emmeans(best_mod, specs = pairwise ~ factor(population), type = "response")

#season
emmeans::emmeans(best_mod, specs = pairwise ~ factor(season), type = "response")


################################################################################
################################################################################









