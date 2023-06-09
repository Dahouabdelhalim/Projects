#### Supplemental Code for Madsen, Vander Meiden & Shizuka: Social partners and temperature jointly affect morning foraging activity of small birds in winter#####
### PART 1: Temperature Analyses#####

#generates stats for analysis of feeder visits and temperature
#generates Figure 2 & Supplemental Figure 3 & 4


#######required libraries
require(tidyverse)
require(lubridate)
require(mgcv)
require(lme4)
library(rsq)
library(MuMIn)
library(arm) 
#####

##load visitation data
all_visits = read.csv("all_visits.csv") #loads 'all_visits' dataframe

##import raw weather data, then summarize it into nightly low temperature
weather_dat = read.csv("weather_raw.csv")
weather <- weather_dat %>%
  mutate(Date = mdy(Date)) %>%
  mutate(datetime = mdy_hm(date_time)) %>%
  mutate(Hour = hour(datetime)) %>%
  filter(Hour >= 19 | Hour <= 4) %>%
  mutate(Lagdate = ifelse(Hour <= 4, paste0(lag(Date)), paste0(Date)))  %>%
  group_by(Date) %>%
  summarize(nightlows = min(Temp_C)) %>%
  ungroup() %>%
  mutate(index = as.numeric(Date-min(Date)+1))

#import individual bird data
birds <- read.csv("RFID_Records_filtered.csv")
birds = birds %>% dplyr::select(-Date)

## Species Models for DOWO and WBNU
morn_visits <- all_visits %>%
  mutate(Hour = hour(Datetime)) %>%
  mutate(Date = date(Datetime)) %>%
  filter(Hour >= 6 & Hour <= 11) %>%
  group_by(RFID, Date) %>%
  summarise(sumvisits = n()) %>%
  ungroup() %>%
  left_join(weather, by = "Date") %>%
  left_join(birds, by = "RFID") %>%
  dplyr::select(Date, RFID, Species, Weight, Sex, nightlows, sumvisits) %>%
  filter(Date < "2019-03-11" & Date > "2019-01-25") %>%
  filter(Species == "DOWO" | Species == "WBNU") %>%
  dplyr::select(Date, RFID, Species, Sex, sumvisits,nightlows)

morn_visits$RFID <- as.factor(morn_visits$RFID)


md <- gamm(sumvisits ~ s(nightlows, fx = FALSE, bs = "tp") + s(RFID, bs = "re"),
           family = poisson,
           data = morn_visits %>% filter(Species == "DOWO"))
summary(md$gam)
summary(md$lme)
plot(md$gam)

mw <- gamm(sumvisits ~ s(nightlows, fx = FALSE, bs = "tp") + s(RFID, bs = "re"),
           family = poisson,
           data = morn_visits %>% filter(Species == "WBNU"))
summary(mw$gam)
summary(mw$lme)
plot(mw$gam)

### species level LMMs
test_dowo <- glmer(sumvisits ~ scale(nightlows) + (1|RFID), data = morn_visits %>% filter(Species == "DOWO"), family = "poisson")


test_wbnu <- glmer(sumvisits ~ scale(nightlows) + (1|RFID), data = morn_visits %>% filter(Species == "WBNU"), family = "poisson")


summary(test_dowo) 
summary(test_wbnu)
r.squaredGLMM(test_dowo)
r.squaredGLMM(test_wbnu)


## Figure 2
## Panel A
md_temp <- morn_visits %>%
  filter(Species == "DOWO")
md_gp <- as.data.frame(predict(md$gam, re.form = TRUE, se = TRUE, type = "response", exclude = s(RFID)))
md_pred <- cbind(md_temp, md_gp)

pana <- ggplot(md_pred) +
  geom_line(aes(nightlows, fit)) +
  geom_ribbon(data = md_pred, aes(x = nightlows, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)), linetype = 2, alpha = 0.2) +
  theme_bw() +
  labs(x = "", y = "", title = "Downy Woodpeckers") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

## Panel B
mw_temp <- morn_visits %>%
  filter(Species == "WBNU")
mw_gp <- as.data.frame(predict(mw$gam, re.form = TRUE, se = TRUE, type = "response", exclude = s(RFID)))
mw_pred <- cbind(mw_temp, mw_gp)

panb <- ggplot(mw_pred) +
  geom_line(aes(nightlows, fit)) +
  geom_ribbon(data = mw_pred, aes(x = nightlows, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)), linetype = 2, alpha = 0.2) +
  theme_bw() +
  labs(x = "", y = "", title = "White-breasted Nuthatches") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

test <- rbind(md_pred, mw_pred)
test <- test %>%
  mutate(SpeciesName = ifelse(Species == "DOWO", paste0("Downy Woodpecker"), paste0("White-breasted Nuthatch")))


ggplot(test) +
  geom_line(aes(nightlows, fit)) +
  geom_ribbon(data = test, aes(x = nightlows, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)), linetype = 2, alpha = 0.2) +
  theme_bw() +
  labs(x = "Lowest Overnight Temperature (\\u00B0C)", y = "Morning Foraging Activity (Visits/Individual)") +
  facet_wrap("SpeciesName", ncol = 2) +
  theme(strip.text.x = element_text(size = 16), strip.background = element_blank(), strip.placement = "outside", axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13))

ggsave("Figure2.tiff")

## Figure S1
dowo_df <- morn_visits %>% 
  filter(Species == "DOWO")
ggplot(dowo_df, aes(nightlows, sumvisits)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("RFID") +
  theme_bw() +
  labs(title = "Downy Woodpeckers", x = "Lowest Overnight Temperature (\\u00B0C)" , y = "Morning Foraging Activity (total visits)") +
  theme(plot.title = element_text(size = 22, face = "bold"), axis.title.x = element_text(size = 22, face = "bold"), axis.title.y = element_text(size = 22, face = "bold"), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))

wbnu_df <- morn_visits %>% 
  filter(Species == "WBNU")
ggplot(wbnu_df, aes(nightlows, sumvisits)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("RFID") +
  theme_bw() +
  labs(title = "White-breasted Nuthatches", x = "Lowest Overnight Temperature (\\u00B0C)" , y = "Morning Foraging Activity (total visits)") +
  theme(plot.title = element_text(size = 22, face = "bold"), axis.title.x = element_text(size = 22, face = "bold"), axis.title.y = element_text(size = 22, face = "bold"), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))

##