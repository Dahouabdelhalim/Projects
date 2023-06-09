

library(tidyverse)
library(lubridate)
library(lme4)       
library(lmerTest)
library(multcomp)

#            Modelling of actual vs. reverse simulated migration
#--------------------------------------------------------------------------------
# aloft data set
aloft <- read.csv("aloftData20211112.csv")
head(aloft)
# Let's separate dataset by outbound and inbound route
aloft_out <- aloft %>% filter(route=="outbound")
aloft_in <- aloft %>% filter(route=="inbound")

# 1. "wind support"-  outbound model
WS_out_mod <- lmer(WS~datatype+(1+datatype|trackId), data=aloft_out)
summary(WS_out_mod)

# 2. "wind support" - inbound model
WS_inb_mod <- lmer(WS~datatype+(1+datatype|trackId), data=aloft_in)
summary(WS_inb_mod)

# 3. "thermal" - outbound
thermal_out_mod <- lmer(thermal~datatype+(1+datatype|trackId), data=aloft_out)
summary(thermal_out_mod)

# 4. "Thermal" - inbound
thermal_in_mod <- lmer(thermal~datatype+(1+datatype|trackId), data=aloft_in)
summary(thermal_in_mod)


# on-ground dataset
ground <- read.csv("groundData20211112.csv")

# Let's separate dataset by outbound and inbound route
ground_out <- ground %>% filter(route=="outbound")
ground_in <- ground %>% filter(route=="inbound")

# 5. "temperature-outbound" model
temp_out_mod <- lmer(grTempC~datatype+(1+datatype|trackId), data=ground_out)
summary(temp_out_mod)

# 6. "temperature-inbound"
temp_in_mod <- lmer(grTempC~datatype+(1+datatype|trackId), data=ground_in)
summary(temp_in_mod)

# 7. "NDVI-outbound"
NDVI_out_mod <- lmer(ndvi~datatype+(1+datatype|trackId), data=ground_out)
summary(NDVI_out_mod)

# 8. "NDVI-inbound"
NDVI_in_mod <- lmer(ndvi~datatype+(1+datatype|trackId), data=ground_in)
summary(NDVI_in_mod)




