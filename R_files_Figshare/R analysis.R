## Covid-19 will embolden autocrats - To do what?

##install.packages("arsenal", dependencies = TRUE)
##install.packages("ggplot2", dependencies = TRUE)
##install.packages("plm", dependencies = TRUE) 
##install.packages("tidyverse", dependencies = TRUE)
##install.packages("maps", dependencies = TRUE)
##install.packages("ggmap", dependencies = TRUE)
##install.packages("mapdata", dependencies = TRUE)
##install.packages(c("rnaturalearth", "rnaturalearthdata"))
##install.packages("dplyr", dependencies = TRUE)
##install.packages("viridis", dependencies = TRUE)
##install.packages("rnaturalearth", dependencies = TRUE)
##install.packages("scales", dependencies = TRUE)
##install.packages("survival", dependencies = TRUE)
##install.packages("survminer", dependencies = TRUE)
##install.packages("lmtest", dependencies = TRUE)
##install.packages("eha", dependencies = TRUE)
##install.packages("sandwich", dependencies = TRUE)
##install.packages("plyr", dependencies = TRUE)
##install.packages("betareg")
##install.packages("gridExtra")

remove.packages(all)

## update.packages(ask = FALSE, checkBuilt = TRUE)

library(arsenal)
library(ggplot2)
library(MASS)
library(tidyverse)
library(maps)
library(ggmap)
library(mapdata)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(ggmap)
library(maps)
library(stringr)
library(zoo)
library(plyr)
library(tidyr)
library(viridis)
library(rnaturalearth)
library(scales)
library(tidyr)
library(survival)
library(survminer)
library(lmtest)
library(sandwich)
library(car)
library(broom)
library(purrr)
library(texreg)
library(betareg)
library(gridExtra)
library(stargazer)

## loading V-Dem Data

V_Dem_Data <- read.csv(file.choose(), header = TRUE, 
                       stringsAsFactors = FALSE, na.strings = "NA",
                       skipNul = TRUE)

## Loading response tracker data into R

response_tracker <- subset(read.csv(file.choose(), header = TRUE, 
                             stringsAsFactors = FALSE, na.strings = "NA",
                             skipNul = TRUE), Jurisdiction == "NAT_TOTAL")

## loading further control variables 

GDP_per_cap <- read.csv(file.choose(), skip = 4, header = TRUE, 
                        stringsAsFactors = FALSE, na.strings = "NA",
                        skipNul = TRUE)

population <- read.csv(file.choose(), skip = 4, header = TRUE, 
                       stringsAsFactors = FALSE, na.strings = "NA",
                       skipNul = TRUE)

health_sector_capacity <- read.csv(file.choose(), header = TRUE, 
                                      stringsAsFactors = FALSE, na.strings = "NA",
                                      skipNul = TRUE)


## data preparation for merging and processing 
## changing column and country names to match and merge

colnames(V_Dem_Data)[1] <- "CountryName"
colnames(GDP_per_cap)[1] <- "CountryName"
colnames(population)[1] <- "CountryName"

summary(comparedf(aggregate(StringencyIndexForDisplay ~ CountryName, response_tracker, max), subset(V_Dem_Data, year == 2019, CountryName), by = "CountryName"))

response_tracker$CountryName[response_tracker$CountryName == "Congo"] <- "Republic of the Congo"
response_tracker$CountryName[response_tracker$CountryName == "Cote d'Ivoire"] <- "Ivory Coast"
response_tracker$CountryName[response_tracker$CountryName == "Democratic Republic of Congo"] <- "Democratic Republic of the Congo"
response_tracker$CountryName[response_tracker$CountryName == "Gambia"] <- "The Gambia"
response_tracker$CountryName[response_tracker$CountryName == "Myanmar"] <- "Burma/Myanmar"
response_tracker$CountryName[response_tracker$CountryName == "Slovak Republic"] <- "Slovakia"
response_tracker$CountryName[response_tracker$CountryName == "United States"] <- "United States of America"
response_tracker$CountryName[response_tracker$CountryName == "Kyrgyz Republic"] <- "Kyrgyzstan"

summary(comparedf(subset(V_Dem_Data, year == 2019, CountryName), GDP_per_cap, by = "CountryName"))

GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Myanmar"]                     <- "Burma/Myanmar"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Cabo Verde"]                  <- "Cape Verde"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Congo, Dem. Rep."]            <- "Democratic Republic of the Congo"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Egypt, Arab Rep."]            <- "Egypt"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Hong Kong SAR, China"]        <- "Hong Kong"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Iran, Islamic Rep."]          <- "Iran"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Cote d'Ivoire"]               <- "Ivory Coast"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Kyrgyz Republic"]             <- "Kyrgyzstan"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Lao PDR"]                     <- "Laos"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Korea, Dem. Peopleâ???Ts Rep."] <- "North Korea"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Congo, Rep."]                 <- "Republic of the Congo"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Russian Federation"]          <- "Russia"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Slovak Republic"]             <- "Slovakia"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Korea, Rep."]                 <- "South Korea"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Syrian Arab Republic"]        <- "Syria"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Gambia, The"]                 <- "The Gambia"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "United States"]               <- "United States of America"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Venezuela, RB"]               <- "Venezuela"
GDP_per_cap$CountryName[GDP_per_cap$CountryName == "Yemen, Rep."]                 <- "Yemen"

colnames(GDP_per_cap)[64] <- "GDP_per_capita_2019"

summary(comparedf(subset(V_Dem_Data, year == 2019, CountryName), population, by = "CountryName"))

population$CountryName[population$CountryName == "Myanmar"]                     <- "Burma/Myanmar"
population$CountryName[population$CountryName == "Cabo Verde"]                  <- "Cape Verde"
population$CountryName[population$CountryName == "Congo, Dem. Rep."]            <- "Democratic Republic of the Congo"
population$CountryName[population$CountryName == "Egypt, Arab Rep."]            <- "Egypt"
population$CountryName[population$CountryName == "Hong Kong SAR, China"]        <- "Hong Kong"
population$CountryName[population$CountryName == "Iran, Islamic Rep."]          <- "Iran"
population$CountryName[population$CountryName == "Cote d'Ivoire"]               <- "Ivory Coast"
population$CountryName[population$CountryName == "Kyrgyz Republic"]             <- "Kyrgyzstan"
population$CountryName[population$CountryName == "Lao PDR"]                     <- "Laos"
population$CountryName[population$CountryName == "Korea, Dem. Peopleâ???Ts Rep."] <- "North Korea"
population$CountryName[population$CountryName == "Congo, Rep."]                 <- "Republic of the Congo"
population$CountryName[population$CountryName == "Russian Federation"]          <- "Russia"
population$CountryName[population$CountryName == "Slovak Republic"]             <- "Slovakia"
population$CountryName[population$CountryName == "Korea, Rep."]                 <- "South Korea"
population$CountryName[population$CountryName == "Syrian Arab Republic"]        <- "Syria"
population$CountryName[population$CountryName == "Gambia, The"]                 <- "The Gambia"
population$CountryName[population$CountryName == "United States"]               <- "United States of America"
population$CountryName[population$CountryName == "Venezuela, RB"]               <- "Venezuela"
population$CountryName[population$CountryName == "Yemem"]                       <- "Yemen"

colnames(population)[64] <- "population2019"


summary(comparedf(health_sector_capacity, subset(V_Dem_Data, year == 2019, CountryName), by = "CountryName"))

health_sector_capacity$CountryName[health_sector_capacity$CountryName == "Congo (Democratic Republic)"] <- "Democratic Republic of the Congo"
health_sector_capacity$CountryName[health_sector_capacity$CountryName == "Congo (Brazzaville)"] <- "Republic of the Congo"
health_sector_capacity$CountryName[health_sector_capacity$CountryName == "eSwatini (Swaziland)"] <- "Eswatini"
health_sector_capacity$CountryName[health_sector_capacity$CountryName == "Gambia"] <- "The Gambia"
health_sector_capacity$CountryName[health_sector_capacity$CountryName == "Kyrgyz Republic"] <- "Kyrgyzstan"
health_sector_capacity$CountryName[health_sector_capacity$CountryName == "Myanmar"] <- "Burma/Myanmar"
health_sector_capacity$CountryName[health_sector_capacity$CountryName == "Cabo Verde"] <- "Cape Verde"
health_sector_capacity$CountryName[health_sector_capacity$CountryName == "CÃ´te d'Ivoire"] <- "Ivory Coast"


## tackling the issue of missing data/ LOCF
## imputing missing values in the GDP data up to 2017

GDP_per_cap$X2018 <- ifelse(is.na(GDP_per_cap$X2018),
                            GDP_per_cap$X2017,
                            GDP_per_cap$X2018)

GDP_per_cap$GDP_per_capita_2019 <- ifelse(is.na(GDP_per_cap$GDP_per_capita_2019),
                           GDP_per_cap$X2018,
                           GDP_per_cap$GDP_per_capita_2019)


## merge datasets into one to derive datasets for analysis from 

V_Dem_tracker <- subset(join(join(join(join(join(join(join(
  subset(response_tracker, ,
         c(CountryName, CountryCode, 
           Date, ConfirmedCases, 
           ConfirmedDeaths, StringencyIndex, 
           StringencyIndexForDisplay, GovernmentResponseIndex, 
           GovernmentResponseIndexForDisplay, ContainmentHealthIndex, 
           ContainmentHealthIndexForDisplay,
           EconomicSupportIndex, EconomicSupportIndexForDisplay)), 
  subset(V_Dem_Data, year == 2019, 
         c(CountryName, year, country_text_id, COWcode, 
           v2x_regime, v2svindep,
           v2xnp_client, v2xnp_client_codelow, 
           v2xnp_client_codehigh, v2xnp_client_sd,
           v2exl_legitperf,
           v2exl_legitperf_codelow, v2exl_legitperf_codehigh,
           v2exl_legitperf_mean, v2exl_legitperf_sd,
           v2exl_legitperf_nr,
           v2xps_party,
           e_regionpol_6C, v2x_regime_amb)),
  by = "CountryName", match = "all"), 
  GDP_per_cap[,c(1,64)], by = "CountryName", match = "all"),
  population[,c(1,64)], by = "CountryName", match = "all"),
  subset(V_Dem_Data, year == 2019, 
         c(CountryName, e_wbgi_gee)), 
  by = "CountryName", match = "all"),
  subset(V_Dem_Data, year == 2019, 
         c(CountryName, v2regendtype)), 
  by = "CountryName", match = "all"),
  subset(V_Dem_Data, year == 2019, 
         c(CountryName, v2regint)), 
  by = "CountryName", match = "all"),
  health_sector_capacity,
  by = "CountryName", match = "all"),
  v2x_regime <2 & v2svindep == 1 & v2regint == 1)


## rescale indicators on performance legitimation, clientelism and state capacity
## to a 1-10 scale preserving the original intervalls 

V_Dem_tracker$v2exl_legitperf <- rescale(V_Dem_tracker$v2exl_legitperf, to = c(0,10), from = c(-5, 5))

V_Dem_tracker$v2xnp_client <- rescale(V_Dem_tracker$v2xnp_client, to = c(0,10), from = c(0, 1))

V_Dem_tracker$e_wbgi_gee <- rescale(V_Dem_tracker$e_wbgi_gee, to = c(0,10), from = c(-2.5, 2.5))

V_Dem_tracker$v2xps_party <- rescale(V_Dem_tracker$v2xps_party, to = c(0,10), from = c(0, 1))


## incorporating regional dummy variables

V_Dem_tracker$East.Europe_Centr.Asia <- ifelse(V_Dem_tracker$e_regionpol_6C == 1,
                                               1, 0)

V_Dem_tracker$LA.America_Carribean <- ifelse(V_Dem_tracker$e_regionpol_6C == 2,
                                             1, 0)

V_Dem_tracker$MiddleEast_N.Africa <- ifelse(V_Dem_tracker$e_regionpol_6C == 3,
                                            1, 0)

V_Dem_tracker$Sub.Sahara_Africa <- ifelse(V_Dem_tracker$e_regionpol_6C == 4,
                                          1, 0)

V_Dem_tracker$West.Europe_N.America <- ifelse(V_Dem_tracker$e_regionpol_6C == 5,
                                              1, 0)

V_Dem_tracker$Asia_Pacific <- ifelse(V_Dem_tracker$e_regionpol_6C == 6,
                                          1, 0)

## renaming columns for later interpretation

colnames(V_Dem_tracker)[19] <- "clientelism_inf.coopt"

colnames(V_Dem_tracker)[23] <- "perf_legitimation"

colnames(V_Dem_tracker)[29] <- "party_institutionalization"

colnames(V_Dem_tracker)[34] <- "state_capacity"

## to control for new cases, the increase in confirmed cases is calculated from 
## the cumulated case reported, normalised on 100.000 citizens

## normalise Covid-19 cases and Deaths on 100.000 citizens

V_Dem_tracker$Cases_per_100.000_cit <- (V_Dem_tracker$ConfirmedCases/(V_Dem_tracker$population2019/100000))

V_Dem_tracker$Deaths_per_100.000_cit <- (V_Dem_tracker$ConfirmedDeaths/(V_Dem_tracker$population2019/100000))

V_Dem_tracker <- V_Dem_tracker %>%
                 group_by(CountryName) %>%
                 fill(Cases_per_100.000_cit,
                      .direction = "down") %>%
                 fill(Deaths_per_100.000_cit,
                      .direction = "down")

V_Dem_tracker$Cases_per_100.000_cit <- ifelse(is.na(V_Dem_tracker$Cases_per_100.000_cit) &
                                                V_Dem_tracker$CountryName != "Turkmenistan",
                                              0, V_Dem_tracker$Cases_per_100.000_cit)

V_Dem_tracker$Deaths_per_100.000_cit <- ifelse(is.na(V_Dem_tracker$Deaths_per_100.000_cit) &
                                                V_Dem_tracker$CountryName != "Turkmenistan",
                                              0, V_Dem_tracker$Deaths_per_100.000_cit)

## then calculate the day to day difference
## giving the new cases every week

diffmod <- function(x) {
  c(x[[1]], diff(x))
}

auxiliary_function02 <- V_Dem_tracker %>% 
                        group_by(CountryName) %>%
                        mutate(new_cases_on_100.000_citizens = c(0, diff(Cases_per_100.000_cit, 
                               lag = 1, differences = 1))) %>%
                        mutate(new_deaths_on_100.000_citizens = c(0, diff(Deaths_per_100.000_cit, 
                               lag = 1, differences = 1)))

## calculating the new cases the previous week normalised on 1.000 citizens
## and incorporating it into the panel dataset

V_Dem_tracker02 <- auxiliary_function02 %>%
                   group_by(CountryName) %>%
                   mutate(new_cases_prev_week = lag(zoo::rollsumr(new_cases_on_100.000_citizens, 
                                                 7, fill = 0), default = 0)) %>%
                   mutate(new_deaths_prev_week = lag(zoo::rollsumr(new_deaths_on_100.000_citizens, 
                                                 7, fill = 0), default = 0)) %>%
                   ungroup

## limiting dataset to first 100 days after first domestic case

## limit response tracker dataset to 100 days after first case

V_Dem_tracker02 <- subset(V_Dem_tracker02, ConfirmedCases > 0)
V_Dem_tracker02 <- group_by(V_Dem_tracker02, CountryName) %>% slice_head(n=100)

## create dataset incorporating the maximum index value for analysis on measure
## stringency
## line 383ff. filter the disease burden previous to implementation
## of most stringent measures

max_response <- subset(join(join(join(
                aggregate(StringencyIndexForDisplay ~ CountryName,
                          V_Dem_tracker02, max),
                aggregate(Cases_per_100.000_cit ~ CountryName,
                          V_Dem_tracker02, max),
                by = "CountryName", match = "first"),
                subset(V_Dem_tracker02, ,c(CountryName, country_text_id, COWcode,
                                           v2x_regime, v2svindep, 
                                           clientelism_inf.coopt, perf_legitimation,
                                           party_institutionalization,
                                           state_capacity, health_capacity_index, 
                                           GDP_per_capita_2019,  
                                           v2regint, v2x_regime_amb,
                                           East.Europe_Centr.Asia,
                                           LA.America_Carribean, MiddleEast_N.Africa,
                                           Sub.Sahara_Africa, West.Europe_N.America,
                                           Asia_Pacific, population2019)),
                by = "CountryName", match = "first"),
                V_Dem_tracker02 %>% 
                  group_by(CountryName) %>% 
                  slice(start = 1, end = (which.max(StringencyIndexForDisplay)-1)) %>% 
                  top_n(n = 1, wt = Cases_per_100.000_cit) %>% 
                  dplyr::select(CountryName, Cases_per_100.000_cit),
                by = "CountryName", match = "first"),
                v2x_regime <2 & v2svindep == 1 & v2regint == 1)


  
## extracting the new cases one week prior to stringest measure implementation
## and adding it to the dataset

new_cases_prev_week <- subset(V_Dem_tracker02, StringencyIndexForDisplay ==
                              ave(StringencyIndexForDisplay, CountryName,
                              FUN = function(x) max(x)),
                              c(CountryName, new_cases_prev_week, Date))

new_cases_on_100.000_citizens <- subset(V_Dem_tracker02, StringencyIndexForDisplay ==
                                        ave(StringencyIndexForDisplay, CountryName,
                                        FUN = function(x) max(x)),
                                        c(CountryName, new_cases_on_100.000_citizens, Date))

new_deaths_prev_week <- subset(V_Dem_tracker02, StringencyIndexForDisplay ==
                                ave(StringencyIndexForDisplay, CountryName,
                                    FUN = function(x) max(x)),
                              c(CountryName, new_deaths_prev_week, Date))

new_deaths_on_100.000_citizens <- subset(V_Dem_tracker02, StringencyIndexForDisplay ==
                                          ave(StringencyIndexForDisplay, CountryName,
                                              FUN = function(x) max(x)),
                                        c(CountryName, new_deaths_on_100.000_citizens, Date))

max_response <- join(join(join(join(max_response, new_cases_prev_week,
                     by = "CountryName", match = "first"),
                     new_cases_on_100.000_citizens,
                     by = "CountryName", match = "first"),
                     new_deaths_prev_week,
                     by = "CountryName", match = "first"),
                     new_deaths_on_100.000_citizens,
                     by = "CountryName", match = "first")


## changing column names and making for appropriate data

colnames(max_response)[colnames(max_response) == "Cases_per_100.000_cit"] <- "total_disease_burden"
colnames(max_response)[colnames(max_response) == "Cases_per_100.000_cit.1"] <- "disease_burden_prev"
colnames(max_response)[29] <- "Date2"

max_response$new_cases_prev_week <- max_response$new_cases_prev_week + 1

max_response$new_cases_on_100.000_citizens <- max_response$new_cases_on_100.000_citizens + 1

max_response$new_deaths_prev_week <- max_response$new_deaths_prev_week + 1

max_response$new_deaths_on_100.000_citizens <- max_response$new_deaths_on_100.000_citizens + 1

max_response$disease_burden_prev <- max_response$disease_burden_prev + 1

max_response$GDP_per_capita_2019 <- max_response$GDP_per_capita_2019/100



## with data processing mostly done, the next stage of analysis incorporates
## creation of maps and scatterplots

## incorporate and process mapping data provided in the ggplot2 package

world_map <- subset(ne_countries(scale = "medium", returnclass = "sf"), sovereignt != "Antarctica")

colnames(world_map)[4] <- "CountryName"

summary(comparedf(max_response, world_map, by = "CountryName"))

world_map$CountryName[world_map$CountryName == "Myanmar"]                     <- "Burma/Myanmar"
world_map$CountryName[world_map$CountryName == "Swaziland"]                   <- "Eswatini"
world_map$CountryName[world_map$CountryName == "Republic of Congo"]           <- "Republic of the Congo"
world_map$CountryName[world_map$CountryName == "Republic of Serbia"]          <- "Serbia"
world_map$CountryName[world_map$CountryName == "United Republic of Tanzania"] <- "Tanzania"

mapping_data <- merge(world_map, max_response, by = "CountryName",
                      all.x = TRUE)

mapping_data$v2x_regime[mapping_data$v2x_regime == 0] <- "Closed Autocracy"

mapping_data$v2x_regime[mapping_data$v2x_regime == 1] <- "Electoral Autocracy"

mapping_data$v2x_regime[is.na(mapping_data$v2x_regime)] <- "Not Covered"


## map paper coverage based on Data availability

map_coverage <- ggplot(mapping_data, aes()) +
                theme_minimal() +
                geom_sf(aes(fill = v2x_regime)) + 
                scale_fill_manual(values=c("red3", "blue3", "grey")) +
                ggtitle("Countries covered in this Paper") + 
                guides(fill=guide_legend(title="Regime Type")) +
                theme(legend.position = "bottom") 


## mapping the previous burden of countries covered

map_previous_disease_burden <- ggplot(mapping_data, aes()) +
                               theme_minimal() +
                               geom_sf(aes(fill = log(disease_burden_prev-1))) +
                               scale_fill_viridis(begin = 0.15, option = "plasma", direction = -1, na.value = "grey") +
                               ggtitle("Disease Burden by Country one day before most stringent measures") + 
                               guides(fill=guide_legend(title="Infected per 100.000 citizens logged")) +
                               theme(legend.position = "bottom")


## mapping the Stringency Index maximum values by Country

map_stirngency_index <- ggplot(mapping_data, aes()) +
                        theme_minimal() +
                        geom_sf(aes(fill = StringencyIndexForDisplay), color = "black") +
                        scale_fill_viridis(begin = 0.15, direction = -1, option = "plasma", na.value = "grey") +
                        ggtitle("Maximum Stringency Index by Country") + 
                        guides(fill=guide_legend(title="Stringency Index")) +
                        theme(legend.position="bottom")

## mapping Performance Legitimation by Country

map_perf_legit <- ggplot(mapping_data, aes()) +
                  theme_minimal() +
                  geom_sf(aes(fill = perf_legitimation), color = "black") +
                  scale_fill_viridis(begin = 0.15, direction = -1, option = "plasma") +
                  ggtitle("Performance Legitimation by Country") +
                  guides(fill=guide_legend(title="Performance Legitimation")) +
                  theme(legend.position="bottom")

## mapping Clientelism by Country

map_clientelism <- ggplot(mapping_data, aes()) +
                   theme_minimal() +
                   geom_sf(aes(fill = clientelism_inf.coopt), color = "black") +
                   scale_fill_viridis(begin = 0.15, direction = -1, option = "plasma") +
                   ggtitle("Clientelism by Country") + 
                   guides(fill=guide_legend(title="Clientelism Level")) +
                   theme(legend.position="bottom")

map_institutionalization <- ggplot(mapping_data, aes()) +
                            theme_minimal() +
                            geom_sf(aes(fill = party_institutionalization), color = "black") +
                            scale_fill_viridis(begin = 0.15, direction = -1, option = "plasma") +
                            ggtitle("Party Institutionalization by Country") + 
                            guides(fill=guide_legend(title="Party Institutionalization")) +
                            theme(legend.position="bottom")


## density plot for Stringency Index 

stringency_density <- ggplot(max_response, aes(x = StringencyIndexForDisplay)) + 
                      geom_density(fill = "grey", alpha = 0.4) + 
                      geom_vline(aes(xintercept=mean(StringencyIndexForDisplay)),
                      color="blue", linetype="dashed", size=1) +
                      theme_minimal() +
                      labs(title = "Distribution of the Stringency Index",
                           x = "Stringency Index", y = "Density") +
                      theme(legend.position = "none")


## descriptive statistics on the variables

summary_stringency <- summary(max_response$StringencyIndexForDisplay)

## for independent variables
summary_perf_legit <- summary(max_response$perf_legitimation)
summary_perf_legit_robust <- summary(max_response$perf_legitimation_solid)
summary_perf_legit_imputed <- summary(max_response$perf_legitimation[is.na(max_response$perf_legitimation_solid)])

hist(max_response$perf_legitimation_solid)

plot(density(max_response$StringencyIndexForDisplay))
high_stringency_index <- subset(max_response, StringencyIndexForDisplay > 80)

summary_clientelism <- summary(max_response$clientelism_inf.coopt)

summary_clientelism_closed <- summary(max_response$clientelism_inf.coopt[max_response$v2x_regime == 0])

summary_clientelism_electoral <- summary(max_response$clientelism_inf.coopt[max_response$v2x_regime == 1])

hist(max_response$clientelism_inf.coopt[max_response$v2x_regime == 0])

hist(max_response$clientelism_inf.coopt[max_response$v2x_regime == 1])

subset_high_client <- subset(max_response, clientelism_inf.coopt > 6.220)

subset_high_perf_legit <- subset(max_response, perf_legitimation > 5.693)

chisq.test(subset(max_response, , c(perf_legitimation, clientelism_inf.coopt)))

## for Covid-19 cases

summary_disease_burden_prev <- summary(max_response$disease_burden_prev)

summary_new_cases_prev_week <- summary(max_response$new_cases_prev_week)

plot(density(log(max_response$disease_burden_prev)[!is.na(max_response$disease_burden_prev)]))

plot(density(log(max_response$new_cases_prev_week)[!is.na(max_response$new_cases_prev_week)]))


summary(max_response$clientelism_inf.coopt[max_response$v2x_regime == 0])
summary(max_response$clientelism_inf.coopt[max_response$v2x_regime == 1])

## checking how many countries experienced first case after April

nrow(subset(max_response, first_case_after_march == 1))

## checking if 2019 estimates on performance legitimation are potentially missing non-random

summary(max_response$GDP_per_capita_2019[!is.na(max_response$perf_legitimation_solid)])
summary(max_response$GDP_per_capita_2019[is.na(max_response$perf_legitimation_solid)])

## checking correlation between health and state capacity

cor(max_response$health_capacity_index, max_response$state_capacity, method = "pearson")

## regression analysis regarding the maximum stringency of measures
## regression of Stringency Index on Performance Legitimation

stringency_on_perf_legit <- lm(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                 log(new_cases_prev_week) + GDP_per_capita_2019 +
                                 state_capacity + health_capacity_index +
                                 v2x_regime, 
                               data = max_response, na.action = na.exclude)

summary(stringency_on_perf_legit)

## regression Government Response Index on Clientelism

stringency_on_clientelism <- lm(formula = StringencyIndexForDisplay ~ clientelism_inf.coopt +
                               log(new_cases_prev_week) + GDP_per_capita_2019 +
                               state_capacity + health_capacity_index +
                               v2x_regime, 
                             data = max_response, na.action = na.exclude)

summary(stringency_on_clientelism)


## fitting a model with log(new_cases_prev_week)

stringency_on_client_perf <- lm(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                  clientelism_inf.coopt + log(new_cases_prev_week) +
                                  state_capacity + health_capacity_index +
                                  v2x_regime, 
                                data = max_response, na.action = na.exclude)

summary(stringency_on_client_perf)


## testing linear regression assumptions for this model

Res_v_Fitted <- plot(stringency_on_client_perf, 1)

## testing linearity assumption

Res_v_Fitted_log_cases <- plot(lm(StringencyIndexForDisplay ~ log(new_cases_prev_week),
                                  data = max_response, na.action = na.exclude), 1)

Res_v_Fitted_perf_legit <- plot(lm(StringencyIndexForDisplay ~ perf_legitimation,
                                   data = max_response, na.action = na.exclude), 1)

Res_v_Fitted_client <- plot(lm(StringencyIndexForDisplay ~ clientelism_inf.coopt,
                                   data = max_response, na.action = na.exclude), 1)


## testing normality of errors assumption

QQ_res_plot <- plot(stringency_on_client_perf, 2)

res_plot <- plot(residuals(stringency_on_client_perf))

shapiro_test <- shapiro.test(residuals(stringency_on_client_perf))


## testing homoscedasticity assumption

scale_location_plot <- plot(stringency_on_client_perf, 3)

white_test_stringency_on_client_perf <- bptest(stringency_on_client_perf, ~ perf_legitimation*clientelism_inf.coopt + I(perf_legitimation^2) + I(clientelism_inf.coopt^2) +
                                               clientelism_inf.coopt*log(new_cases_prev_week) + I(log(new_cases_prev_week)^2) +
                                               log(new_cases_prev_week)*GDP_per_capita_2019 + I(GDP_per_capita_2019^2) +
                                               GDP_per_capita_2019*state_capacity + I(state_capacity^2),
                                             data = max_response)

bp_test_stringency_on_client_perf <- bptest(stringency_on_client_perf, data = max_response)


## detecting possible outliners and looking for reason

outlier_test <- outlierTest(stringency_on_client_perf, cutoff = 0.5)

Cooks_dist <- plot(stringency_on_client_perf, 4, id.n = 10)

Res_v_Leverage <- plot(stringency_on_client_perf, 5)

## testing for misspecification

reset_test <- resettest(stringency_on_client_perf, power = 2:3, type = "regressor",
                        data = max_response)

## testing for collinearity in the regresors

Var_Inflation_Factor <- vif(stringency_on_client_perf)


## obtaining heteroscedasticity robust standard errors

stringency_on_client_perf_het_robust <- coeftest(stringency_on_client_perf, vcov = vcovHC(stringency_on_client_perf, method = "arellano", type = "HC3"))


## plotting fitted values

plot(fitted.values(stringency_on_client_perf))

fitted_values_baseline <- data.frame(fitted.values(stringency_on_client_perf))

plot(x = max_response$StringencyIndexForDisplay, y = fitted.values(stringency_on_client_perf))

## fitting a model with alternative specification of Covid-19 cases as
## log(new_deaths_prev_week)

stringency_on_client_perf_deaths <- lm(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                           clientelism_inf.coopt + log(new_deaths_prev_week) +
                                           state_capacity + health_capacity_index + 
                                           v2x_regime, 
                                          data = max_response, na.action = na.exclude)

summary(stringency_on_client_perf_deaths)

AIC(stringency_on_client_perf, stringency_on_client_perf_deaths)
BIC(stringency_on_client_perf, stringency_on_client_perf_deaths)


## testing linear regression assumptions for this model

## testing linearity assumption

Res_v_Fitted_deaths <- plot(stringency_on_client_perf_deaths, 1)

Res_v_Fitted_log_deaths <- plot(lm(StringencyIndexForDisplay ~ log(new_deaths_prev_week),
                                   data = max_response, na.action = na.exclude), 1)


## testing normality of errors assumption

QQ_res_plot_deaths <- plot(stringency_on_client_perf_deaths, 2)

res_plot_dis_burden <- plot(residuals(stringency_on_client_perf_dis_burden))

shapiro_test_deaths <- shapiro.test(residuals(stringency_on_client_perf_deaths))


## testing homoscedasticity assumption

scale_location_plot_dis_burden <- plot(stringency_on_client_perf_dis_burden, 3)

white_test_stringency_client_perf_deaths <- bptest(stringency_on_client_perf_dis_burden, ~ perf_legitimation*clientelism_inf.coopt + I(perf_legitimation^2) + I(clientelism_inf.coopt^2) +
                                                      clientelism_inf.coopt*log(new_deats_prev_week) + I(log(new_deats_prev_week)^2) +
                                                      log(new_deats_prev_week)*GDP_per_capita_2019 + I(GDP_per_capita_2019^2) +
                                                      GDP_per_capita_2019*state_capacity + I(state_capacity^2),
                                                     data = max_response)

bp_test_stringency_on_client_perf_deaths <- bptest(stringency_on_client_perf_deaths, data = max_response)


## detecting possible outliners and looking for reason

outlier_test_deaths <- outlierTest(stringency_on_client_perf_deaths, cutoff = 0.5)

Cooks_dist_deaths <- plot(stringency_on_client_perf_deaths, 4, id.n = 10)

Res_v_Leverage_deaths <- plot(stringency_on_client_perf_deaths, 5)

## testing for misspecification using ramsey reset test

reset_test_dis_burden <- resettest(stringency_on_client_perf_deaths, power = 2:3, type = "regressor",
                                   data = max_response)

## testing for collinearity in the regresors

Var_Inflation_Factor_dis_burden <- vif(stringency_on_client_perf_deaths)


## obtaining heteroscedasticity robust standard errors

stringency_on_client_perf_dis_burden_het_robust <- coeftest(stringency_on_client_perf_deaths, vcov = vcovHC(stringency_on_client_perf_deaths, method = "arellano", type = "HC3"))


## plotting fitted values

plot(fitted.values(stringency_on_client_perf_deaths))

## dealing with the problematic outliers by estimating a robust regression 
## model using "MM" estimator 


stringency_on_client_perf_robust <-  rlm(StringencyIndexForDisplay ~ perf_legitimation + 
                                         clientelism_inf.coopt  + log(new_cases_prev_week) +
                                         state_capacity + health_capacity_index +
                                         v2x_regime, 
                                         na.action = na.exclude, method = "MM",
                                         data = max_response)

summary(stringency_on_client_perf_robust)

## testing linearity assumption

res_v_fitted_robust <- plot(stringency_on_client_perf_robust, 1)

## plotting fitted values

plot(fitted.values(stringency_on_client_perf_robust))


## dealing with the problematic outliers by estimating a robust regression 
## model using "MM" estimator with cases operationalizes as disease burden


stringency_on_client_perf_deaths_robust <-  rlm(StringencyIndexForDisplay ~ perf_legitimation + 
                                                    clientelism_inf.coopt  + log(new_deaths_prev_week) +
                                                    state_capacity + 
                                                    health_capacity_index +
                                                    v2x_regime, 
                                                    data = max_response,
                                                    na.action = na.exclude, method = "MM")

summary(stringency_on_client_perf_deaths_robust)

## testing linearity assumption

res_v_fitted_dis_burden_robust <- plot(stringency_on_client_perf_deaths_robust, 1)

## plotting fitted values

plot(fitted.values(stringency_on_client_perf_deaths_robust))

## ensuring robustness of the results running seperate regressions of closed 
## and electoral autocracies 

## first: closed autocracies

max_response_closed <- subset(max_response, v2x_regime == 0)


stringency_on_client_perf_closed <-  lm(StringencyIndexForDisplay ~ perf_legitimation + 
                                           clientelism_inf.coopt  + log(new_cases_prev_week) +
                                           state_capacity + 
                                           health_capacity_index, 
                                           data = max_response_closed,
                                           na.action = na.exclude)

summary(stringency_on_client_perf_closed)


## estimating robust standard errors

stringency_on_client_perf_closed_het_robust <- coeftest(stringency_on_client_perf_closed, vcov = vcovHC(stringency_on_client_perf_closed, method = "arellano", type = "HC3"))


plot(fitted.values(stringency_on_client_perf_closed))


## ensuring robustness of the results running seperate regressions of closed 
## and electoral autocracies 

## second: electoral autocracies

max_response_electoral <- subset(max_response, v2x_regime == 1)


stringency_on_client_perf_electoral <-  lm(StringencyIndexForDisplay ~ perf_legitimation + 
                                           clientelism_inf.coopt  + log(new_cases_prev_week) +
                                           state_capacity + 
                                           health_capacity_index, 
                                         data = max_response_electoral,
                                         na.action = na.exclude)

summary(stringency_on_client_perf_electoral)


## testing linear regression assumptions for regression of closed autocracies

outlierTest(stringency_on_client_perf_electoral)

plot(stringency_on_client_perf_electoral, 1)

bptest(stringency_on_client_perf_electoral, data = max_response_electoral)

plot(stringency_on_client_perf_electoral, 5)

plot(stringency_on_client_perf_electoral, 4)

outlierTest(stringency_on_client_perf_electoral)

vif(stringency_on_client_perf_electoral)


## estimating robust standard errors

stringency_on_client_perf_electoral_het_robust <- coeftest(stringency_on_client_perf_electoral, vcov = vcovHC(stringency_on_client_perf_electoral, method = "arellano", type = "HC3"))

plot(fitted.values(stringency_on_client_perf_electoral))

## regression for electoral autocracies using death numbers

stringency_on_client_perf_electoral_deaths <-  lm(StringencyIndexForDisplay ~ perf_legitimation + 
                                                  clientelism_inf.coopt  + log(new_deaths_prev_week) +
                                                  state_capacity + 
                                                  health_capacity_index, 
                                                 data = max_response_electoral,
                                                 na.action = na.exclude)

summary(stringency_on_client_perf_electoral_deaths)

stringency_on_client_perf_electoral_deaths_het_robust <- coeftest(stringency_on_client_perf_electoral_deaths, vcov = vcovHC(stringency_on_client_perf_electoral_deaths, method = "arellano", type = "HC3"))


## programming beta regression 

max_response_for_beta <- max_response

max_response_for_beta$StringencyIndexForDisplay <- max_response_for_beta$StringencyIndexForDisplay/100

max_response_for_beta$StringencyIndexForDisplay <- ((max_response_for_beta$StringencyIndexForDisplay*(73-1)+0.5)/73)


stringency_on_client_perf_beta <- betareg(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                          clientelism_inf.coopt  + log(new_cases_prev_week) +
                                          state_capacity + 
                                          health_capacity_index + 
                                          v2x_regime, 
                                         data = max_response_for_beta , na.action = na.exclude, link =  "logit")

summary(stringency_on_client_perf_beta)

stringency_on_client_perf_beta_probit <- betareg(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                                          clientelism_inf.coopt  + log(new_cases_prev_week) +
                                                          state_capacity + 
                                                          health_capacity_index, 
                                                        data = max_response_for_beta , na.action = na.exclude, 
                                                        link =  "probit")

summary(stringency_on_client_perf_beta_probit)

stringency_on_client_perf_beta_loglog <- betareg(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                                          clientelism_inf.coopt  + log(new_cases_prev_week) +
                                                          state_capacity + 
                                                          health_capacity_index, 
                                                        data = max_response_for_beta , na.action = na.exclude, 
                                                        link =  "loglog")

summary(stringency_on_client_perf_beta_loglog)

stringency_on_client_perf_beta_cloglog <- betareg(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                                       clientelism_inf.coopt  + log(new_cases_prev_week) +
                                                       state_capacity + 
                                                       health_capacity_index, 
                                                     data = max_response_for_beta , na.action = na.exclude, 
                                                     link =  "cloglog")

summary(stringency_on_client_perf_beta_cloglog)


AIC(stringency_on_client_perf_beta, 
    stringency_on_client_perf_beta_probit,
    stringency_on_client_perf_beta_loglog,
    stringency_on_client_perf_beta_cloglog)

AIC(stringency_on_client_perf,
    stringency_on_client_perf_robust,
    stringency_on_client_perf_beta)


plot(fitted.values(stringency_on_client_perf_beta))

plot(stringency_on_client_perf_beta, which = 1, type = "pearson", sub.caption = "")
plot(stringency_on_client_perf_beta, which = 2, type = "pearson", sub.caption = "")
plot(stringency_on_client_perf_beta, which = 3, type = "pearson", sub.caption = "")
plot(stringency_on_client_perf_beta, which = 4, type = "pearson", sub.caption = "")


plot(stringency_on_client_perf_beta, which = 5, type = "deviance", sub.caption = "")
plot(stringency_on_client_perf_beta, which = 1, type = "deviance", sub.caption = "")

plot(cooks.distance(stringency_on_client_perf_beta))

stringency_on_client_perf_beta_ex_B_B_N <- betareg(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                                     clientelism_inf.coopt  + log(new_cases_prev_week) +
                                                     state_capacity + 
                                                     health_capacity_index + 
                                                     v2x_regime, 
                                                   data = subset(max_response_for_beta, 
                                                                 CountryName != "Belarus"&
                                                                 CountryName != "Burundi"&
                                                                 CountryName != "Nicaragua"), na.action = na.exclude, link =  "logit")

summary(stringency_on_client_perf_beta_ex_B_B_N)



stringency_on_client_perf_deaths_beta <- betareg(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                                     clientelism_inf.coopt  + log(new_deaths_prev_week) +
                                                     state_capacity + 
                                                     health_capacity_index + 
                                                     v2x_regime, 
                                                    data = max_response_for_beta , na.action = na.exclude, link =  "logit")

summary(stringency_on_client_perf_deaths_beta)

max_response_closed_for_beta <- subset(max_response_for_beta, v2x_regime == 0)

stringency_on_client_perf_closed_beta <- betareg(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                                 clientelism_inf.coopt  + log(new_cases_prev_week) +
                                                 state_capacity + 
                                                 health_capacity_index, 
                                                data = max_response_closed_for_beta , na.action = na.exclude, 
                                                link =  "logit")

summary(stringency_on_client_perf_closed_beta)

waldtest(stringency_on_client_perf_closed_beta)


plot(fitted.values(stringency_on_client_perf_closed_beta))

plot(stringency_on_client_perf_closed_beta, which = 1, type = "pearson", sub.caption = "")
plot(stringency_on_client_perf_closed_beta, which = 2, type = "pearson", sub.caption = "")
plot(stringency_on_client_perf_closed_beta, which = 3, type = "pearson", sub.caption = "")
plot(stringency_on_client_perf_closed_beta, which = 4, type = "pearson", sub.caption = "")


plot(stringency_on_client_perf_closed_beta, which = 5, type = "deviance", sub.caption = "")
plot(stringency_on_client_perf_closed_beta, which = 1, type = "deviance", sub.caption = "")


max_response_electoral_for_beta <- subset(max_response_for_beta, v2x_regime == 1)

stringency_on_client_perf_electoral_beta <- betareg(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                                   clientelism_inf.coopt  + log(new_cases_prev_week) +
                                                   state_capacity + 
                                                   health_capacity_index , 
                                                 data = max_response_electoral_for_beta , na.action = na.exclude, 
                                                 link =  "logit")

summary(stringency_on_client_perf_electoral_beta)

plot(fitted.values(stringency_on_client_perf_electoral_beta))


## visualising the findings in more intuitive tables, enabling comparison

## for the baseline models and the regression on closed and electoral regimes seperately

robust_se_stringency_on_client_perf = vcovHC(stringency_on_client_perf, type = "HC3") %>% diag() %>% sqrt()

robust_se_stringency_on_client_perf_deaths = vcovHC(stringency_on_client_perf_deaths, type = "HC3") %>% diag() %>% sqrt()

robust_se_stringency_on_client_perf_closed = vcovHC(stringency_on_client_perf_closed, type = "HC3") %>% diag() %>% sqrt()

robust_se_stringency_on_client_perf_electoral = vcovHC(stringency_on_client_perf_electoral, type = "HC3") %>% diag() %>% sqrt()


stargazer(list(stringency_on_client_perf, stringency_on_client_perf,
               stringency_on_client_perf_deaths, stringency_on_client_perf_deaths,
               stringency_on_client_perf_electoral, stringency_on_client_perf_electoral),
          type = "html",
          se = list(NULL, robust_se_stringency_on_client_perf,
                             NULL, robust_se_stringency_on_client_perf_deaths,
                             NULL, robust_se_stringency_on_client_perf_electoral),
          title = "Maximum Stringency Index as a Function of Performance Legitimation and Clientelism",
          out = "C:/Users/philipp.becker/Documents/LSE/Publication/Covid/Plots/stringency_ols_output.html",
          column.labels = c("Model 1", "Model 1 HC3",
                                 "Model 2", "Model 2 HC3",
                                 "Model 1 Electoral", "Model 1 Electoral HC3"),
          covariate.labels = c("Performance Legitimation", "Clientelism", "log(New Cases Previous Week)",
                                "log(New Deaths Previous Week)", 
                                "State Capacity", "Health Sector Capacity",
                                "Electoral Autocracy"),
          dep.var.labels = "Maximum Stringency Index (Display Version)",
          df = TRUE, style = "ajps", model.numbers = FALSE, model.names = TRUE)


stargazer(list(stringency_on_client_perf, stringency_on_client_perf_beta, 
               stringency_on_client_perf_deaths, stringency_on_client_perf_deaths_beta,
               stringency_on_client_perf_electoral, stringency_on_client_perf_electoral_beta),
          type = "html",
          title = "Maximum Stringency Index as a Function of Performance Legitimation and Clientelism",
          out = "C:/Users/philipp.becker/Documents/LSE/Publication/Covid/Plots/stringency_ols_and_beta_output.html",
          column.labels = c("Model 1", "Model 1",
                            "Model 2", "Model 2",
                            "Model 1 Electoral", "Model 1 Electoral"),
          covariate.labels = c("Performance Legitimation", "Clientelism", "log(New Cases Previous Week)",
                               "log(New Deaths Previous Week)", 
                               "State Capacity", "Health Sector Capacity",
                               "Electoral Autocracy"),
          colnames = FALSE,
          dep.var.labels = "Maximum Stringency Index (Display Version)",
          df = TRUE)

## plotting outlier treatment

stringency_on_client_perf_electoral_robust <- rlm(StringencyIndexForDisplay ~ perf_legitimation + 
                                                    clientelism_inf.coopt  + log(new_cases_prev_week) +
                                                    state_capacity + 
                                                    health_capacity_index, 
                                                    data = max_response_electoral,
                                                    na.action = na.exclude, method = "MM")

summary(stringency_on_client_perf_electoral_robust)

stringency_on_client_perf_electoral_ex_B_B_N <- lm(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                                     clientelism_inf.coopt  + log(new_cases_prev_week) +
                                                     state_capacity + 
                                                     health_capacity_index, 
                                                   data = subset(max_response_electoral, 
                                                                 CountryName != "Belarus"&
                                                                   CountryName != "Burundi"&
                                                                   CountryName != "Nicaragua"), na.action = na.exclude)

summary(stringency_on_client_perf_electoral_ex_B_B_N)


stargazer(list(stringency_on_client_perf_robust, stringency_on_client_perf_beta_ex_B_B_N,
               stringency_on_client_perf_electoral_robust, stringency_on_client_perf_electoral_ex_B_B_N),
          type = "html",
          title = "Maximum Stringency Index as a Function of Performance Legitimation and Clientelism",
          out = "C:/Users/philipp.becker/Documents/LSE/Publication/Covid/Plots/output_outlier_treatment.html",
          column.labels = c("Model 1", "Model 1 excluding Outliers",
                            "Model 1 Electoral", "Model 1 electoral excluding Outliers"),
          covariate.labels = c("Performance Legitimation", "Clientelism", "log(New Cases Previous Week)",
                               "State Capacity", "Health Sector Capacity",
                               "Electoral Autocracy"),
          colnames = FALSE,
          dep.var.labels = "Maximum Stringency Index (Display Version)",
          df = TRUE)


## constructing model with party institutionalization and disease burden

stringency_on_client_perf_dis_burden <- lm(StringencyIndexForDisplay ~ perf_legitimation + 
                                             clientelism_inf.coopt  + log(disease_burden_prev) +
                                             state_capacity + 
                                             health_capacity_index + v2x_regime, 
                                           data = max_response,
                                           na.action = na.exclude)

summary(stringency_on_client_perf_dis_burden)

stringency_on_party_inst_perf_dis_burden <- lm(StringencyIndexForDisplay ~ perf_legitimation + 
                                               party_institutionalization  + log(disease_burden_prev) +
                                               state_capacity + 
                                               health_capacity_index + v2x_regime, 
                                             data = max_response,
                                             na.action = na.exclude)

summary(stringency_on_party_inst_perf_dis_burden)

stringency_on_client_perf_dis_burden_beta <- betareg(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                                      clientelism_inf.coopt  + log(disease_burden_prev) +
                                                      state_capacity + 
                                                      health_capacity_index + v2x_regime, 
                                                    data = max_response_for_beta , na.action = na.exclude, 
                                                    link =  "logit")

summary(stringency_on_client_perf_dis_burden_beta)


stringency_on_party_inst_perf_dis_burden_beta <- betareg(formula = StringencyIndexForDisplay ~ perf_legitimation + 
                                                       party_institutionalization  + log(disease_burden_prev) +
                                                       state_capacity + 
                                                       health_capacity_index + v2x_regime, 
                                                     data = max_response_for_beta , na.action = na.exclude, 
                                                     link =  "logit")

summary(stringency_on_party_inst_perf_dis_burden_beta)


stargazer(list(stringency_on_client_perf_dis_burden, stringency_on_party_inst_perf_dis_burden,
               stringency_on_client_perf_dis_burden_beta, stringency_on_party_inst_perf_dis_burden_beta),
          type = "html",
          title = "Maximum Stringency Index as a Function of Performance Legitimation and Clientelism/Party Institutionalization",
          out = "C:/Users/philipp.becker/Documents/LSE/Publication/Covid/Plots/output_dis_burden_inst.html",
          column.labels = c("Model 3", "Model 3 Party Inst.",
                            "Model 3", "Model 3 Party Inst."),
          covariate.labels = c("Performance Legitimation", "Clientelism", "Party Institutionalization", 
                               "log(Prev. Disease Burden)",
                               "State Capacity", "Health Sector Capacity",
                               "Electoral Autocracy"),
          colnames = FALSE,
          dep.var.labels = "Maximum Stringency Index (Display Version)",
          df = TRUE)


## visualizing results for closed autocracies

ggplot(max_response_closed, aes(x = clientelism_inf.coopt, y=StringencyIndexForDisplay),) +
  geom_smooth(method = "loess") + geom_point() + geom_smooth(method = "lm", weight = 1, se = FALSE,
                                                             color = "red3") +
  theme_classic() + labs(x = "Clientelism", y = "Stringency Index", 
                         title = "Stringency Index on Clientelism in Closed Autocracies")




## for the seperate regression models on clientelism and performance legitimation


robust_se_stringency_on_perf_legit = vcovHC(stringency_on_perf_legit, type = "HC3") %>% diag() %>% sqrt()

robust_se_stringency_on_clientelism = vcovHC(stringency_on_clientelism, type = "HC3") %>% diag() %>% sqrt()


stargazer(list(stringency_on_perf_legit, stringency_on_perf_legit,
               stringency_on_clientelism, stringency_on_clientelism),
          type = "html", 
          se = list(NULL, robust_se_stringency_on_perf_legit,
                    NULL, robust_se_stringency_on_clientelism),
          title = "Estimates on a seperate Regression on Performance Legitimation and Clientelism",
          out = "C:/Users/philipp.becker/Documents/LSE/Publication/Covid/Plots/output_separate_reg.html",
          column.labels = c("Model 1 Performance Legitimation", "Model 1 Performance Legitimation HC3",
                            "Model 1 Clientelism", "Model 1 Clientelism HC3"),
          covariate.labels = c("Performance Legitimation", "Clientelism", "log(New Cases Previous Week)", 
                               "GDP per Capita", "State Capacity", "Health Sector Capacity", 
                               "Electoral Autocracy"),
          dep.var.labels = "Maximum Stringency Index (Display Version)",
          df = TRUE)


## exporting dataset for sharing

write.csv(max_response,"C:/Users/philipp.becker/Documents/LSE/Publication/Covid/max_response.csv", row.names = FALSE)


## matching analysis for closed autocracies

##install.packages("Matching", dependencies = TRUE)
library(Matching)

max_response_for_matching <- join(join(max_response_closed, subset(V_Dem_Data,
                                                              year == 2018,
                                                              c(CountryName, e_migdppcln, e_wb_pop)),
                                  by = "CountryName", match = "first"),
                                  subset(V_Dem_Data, year == 2019, c(CountryName, e_pelifeex)),
                                  by = "CountryName", match = "first")

max_response_for_matching$one_party <- ifelse(max_response_for_matching$CountryName == "China",
                                              1, ifelse(max_response_for_matching$CountryName == "Cuba",
                                                        1, ifelse(max_response_for_matching$CountryName == "Laos",
                                                                   1,ifelse(max_response_for_matching$CountryName == "Vietnam",
                                                                             1,0))))


max_response_for_matching$monarchy <- ifelse(max_response_for_matching$CountryName == "Bahrain",
                                              1, ifelse(max_response_for_matching$CountryName == "Eswatini",
                                                        1, ifelse(max_response_for_matching$CountryName == "Jordan",
                                                                  1,ifelse(max_response_for_matching$CountryName == "Kuwait",
                                                                           1, ifelse(max_response_for_matching$CountryName == "Marocco",
                                                                                     1, ifelse(max_response_for_matching$CountryName == "Saudi Arabia",
                                                                                               1, ifelse(max_response_for_matching$CountryName == "United Arab Emirated",
                                                                                                         1,0)))))))

max_response_for_matching$personal <- ifelse(max_response_for_matching$CountryName == "Eritrea",
                                              1, ifelse(max_response_for_matching$CountryName == "Uzbekistan",
                                                        1, ifelse(max_response_for_matching$CountryName == "Somalia",
                                                                  1,ifelse(max_response_for_matching$CountryName == "South Sudan",
                                                                           1, ifelse(max_response_for_matching$CountryName == "Sudan",
                                                                                     1, ifelse(max_response_for_matching$CountryName == "Syria",
                                                                                               1, ifelse(max_response_for_matching$CountryName == "Thailand",
                                                                                                         1,0)))))))


colnames(max_response_for_matching)[colnames(max_response_for_matching) == "e_migdppcln"] <- "GDP_per_cap_log_10"

colnames(max_response_for_matching)[colnames(max_response_for_matching) == "e_wb_pop"] <- "population"

colnames(max_response_for_matching)[colnames(max_response_for_matching) == "e_pelifeex"] <- "life_expect"


summary(max_response_for_matching$clientelism_inf.coopt)

max_response_for_matching$treatment_first <- ifelse(max_response_for_matching$clientelism_inf.coopt > 5.805,
                                              1, 0)

max_response_for_matching$treatment_second <- ifelse(max_response_for_matching$clientelism_inf.coopt > 7.348,
                                              1, 0)

##Matching Analysis

Match_analysis_T1 <- Match(Y = max_response_for_matching$StringencyIndexForDisplay,
                           Tr = max_response_for_matching$treatment_first,
                           X = max_response_for_matching$state_capacity &
                               max_response_for_matching$health_capacity_index,
                           BiasAdjust = TRUE, M = 1)

summary.Match(Match_analysis_T1)

Match_analysis_T2 <- Match(Y=max_response_for_matching$StringencyIndexForDisplay,
                           Tr = max_response_for_matching$treatment_second,
                           X = max_response_for_matching$state_capacity &
                             max_response_for_matching$health_capacity_index,
                           BiasAdjust = TRUE, M = 1)

summary.Match(Match_analysis_T2)


## checking match distributions

Balance1 <- MatchBalance(treatment_first ~ state_capacity + I(state_capacity^2) +
               health_capacity_index + I(health_capacity_index^2) ,
             data = max_response_for_matching, ks = TRUE, nboots = 1000)

Balance2 <- MatchBalance(treatment_second ~ state_capacity + I(state_capacity^2) +
               health_capacity_index + I(health_capacity_index^2) ,
             data = max_response_for_matching, ks = TRUE, nboots = 1000)

summary(Balance1)

summary (Balance2)


## propensity score based estimation

prop_score_T1 <- glm(treatment_first ~ state_capacity + I(state_capacity^2) + health_capacity_index + I(health_capacity_index^2),
                  data = max_response_for_matching, family = binomial)

summary(prop_score_T1)

with(summary(prop_score_T1), 1 - deviance/null.deviance)

match_prop_score_T1 <- Match(Y=max_response_for_matching$StringencyIndexForDisplay,
                             Tr = max_response_for_matching$treatment_first,
                             X = prop_score$fitted.values,
                             BiasAdjust = TRUE, M = 1)

summary.Match(match_prop_score_T1)


prop_score_T2 <- glm(treatment_second ~ state_capacity + I(state_capacity^2) + health_capacity_index + I(health_capacity_index^2),
                     data = max_response_for_matching, family = binomial)

summary(prop_score_T2)

with(summary(prop_score_T2), 1 - deviance/null.deviance)

match_prop_score_T2 <- Match(Y=max_response_for_matching$StringencyIndexForDisplay,
                             Tr = max_response_for_matching$treatment_second,
                             X = prop_score$fitted.values,
                             BiasAdjust = TRUE, M = 1)

summary.Match(match_prop_score_T2)


htmlreg(l = list(Match_analysis_T1, Match_analysis_T2,
                 match_prop_score_T1, match_prop_score_T2),
        file = "C:/Users/philipp.becker/Documents/LSE/Publication/Covid/Plots/output_matching.html")


