##Regression models fitted to Willow Warbler first egg dates, incubation status, and moult scores
library(readr)
library(dplyr)
library(glmmTMB)
library(sjPlot)
library(moult)

##First egg phenology
#read data and set up factor variables
egg_df <- readr::read_csv('dryad/wilwa_first_egg_data.csv') %>%
  mutate(year_fac = factor(year))

egg_model <- glmmTMB(feg~ nst*sptdiff+(1|year_fac),data=egg_df)

#display parameter summary (as in Table 1)
summary(egg_model)
tab_model(egg_model)



##Incubation phenology -  We use a probit GLM to estimate the population the population distribution of the timing of brood patch regression, analogous to the moult phenology model propsed by Rothery & Newton (2002): https://doi.org/10.1046/j.1474-919X.2002.00072.x

#read data and set up factor variables
incubation_df <- readr::read_csv('dryad/wilwa_incubation_data.csv') %>%
  mutate(incubation = factor(incubation),
         year_fac = factor(year_fac))
#fit probit glm
bp_binary_3way <- glmmTMB(incubation ~ yday*latst*sptdiff+ (1|year_fac), family = binomial(link = 'probit'), data = incubation_df)

#display parameter summary (as in Table 2)
summary(bp_binary_3way)#link scale
sjPlot::tab_model(bp_binary_3way, digits = 4)#risk ratios


##Moult phenology
#read data and set up factor variables
moult_df <- readr::read_csv('dryad/wilwa_primary_moult_data.csv') %>%
  mutate(sex = factor(sex))
#fit probit glm
moult_model <- moult(pp_mass ~ yday|sptdiff*nst|nst*sptdiff+sex|sex, data = moult_df, type = 5)
#display parameter summary (as in Table 1)
summary(moult_model)
