## Juvenile focal follow analyses
## Ali Galezo

library(dplyr)
library(reshape2)
library(perm)
library(ggplot2)

####################################################
## Point sample analysis
####################################################

## Load juvenile point sample data
ps_juvs <- read.csv("juvenile_pointsamples.csv")

## Calculate proportion of time focal spends on different activites.
ps_rates <-
  ps_juvs %>%
  group_by(focal, sex) %>%
  count(activity_category) %>%
  dcast(., focal + sex ~ activity_category, value.var = "n") %>%
  replace(is.na(.), 0) %>%
  merge(x = .,
        y = ps_juvs %>% distinct(focal, duration, num_points, num_obs),
        by = "focal",
        all.x = TRUE) %>%
  mutate(num_points_nounks = num_points - UNKNOWN) %>%
  mutate(SOCrate = SOCIAL/num_points_nounks) %>%
  mutate(FORrate = FORAGE/num_points_nounks) %>%
  mutate(TRVrate = TRAVEL/num_points_nounks) %>%
  mutate(RSTrate = REST/num_points_nounks) %>%
  filter(duration >= 180) %>%
  filter(num_obs >= 2)

## Check sample sizes.
ps_rates %>% count(sex)
ps_rates %>%
  summarize(num_females = sum(sex == "FEMALE"),
            num_males = sum(sex == "MALE"),
            follows = sum(num_obs))
ps_rates %>%
  distinct(focal, duration) %>%
  summarize(mins = sum(duration),
            hours = mins/60)

## Compare proportion time spent socializing
ps_summarize <- function(x){list(mean = mean(x),
                                 se = sd(x)/sqrt(length(x)))}
male_soc <- ps_rates[ps_rates$sex == "MALE",]$SOCrate
female_soc <- ps_rates[ps_rates$sex == "FEMALE",]$SOCrate
ps_summarize(male_soc)
ps_summarize(female_soc)
permTS(male_soc, female_soc, alternative = "two.sided", exact = TRUE)

## Compare proportion time spent foraging
male_for <- ps_rates[ps_rates$sex == "MALE",]$FORrate
female_for <- ps_rates[ps_rates$sex == "FEMALE",]$FORrate
ps_summarize(male_for)
ps_summarize(female_for)
permTS(male_for, female_for, alternative = "two.sided", exact = TRUE)

## Data for figure
plotdata <- data.frame(
  Sex = c("Female","Male","Female","Male"),
  Activity = c("Social","Social","Forage","Forage"),
  Mean = c(mean(female_soc),
           mean(male_soc),
           mean(female_for),
           mean(male_for)),
  SE = c(sd(female_soc)/sqrt(length(female_soc)),
         sd(male_soc)/sqrt(length(male_soc)),
         sd(female_for)/sqrt(length(female_for)),
         sd(male_for)/sqrt(length(male_for)))
)

## Significance star coordinates
significance.stars <- data.frame(Activity = c("Forage"),
                                 Mean = c(0.65),
                                 Sex = NA,
                                 SE = NA)

## Arc coordinates
r <- .3
t <- seq(0, 180, by = 1) * pi / 180
x <- r * cos(t)
y <- r*0.2 * sin(t)
y[20:162] <- y[20]
arc.df <- data.frame(Activity = x, Mean = y, Sex = NA, SE = NA)

## Plot
ggplot(
  data = plotdata,
  aes(x = Activity, y = Mean, fill = Sex, ymin = Mean - SE, ymax = Mean + SE)) +
  geom_bar(stat = "identity", position=position_dodge(), color = "black") +
  scale_fill_manual(values=c("#3d3d3d", "#FFFFFF")) +
  ylab("Proportion of Time") +
  geom_errorbar(position=position_dodge(width = 0.9), width = 0.3, color = "#000000") +
  theme_classic() +
  geom_text(data = significance.stars, label = "*", size = 9) +
  geom_line(data = arc.df, aes(Activity+1, Mean+0.62), lty = 2) +
  scale_y_continuous(limits = c(0,0.7), expand = c(0,0))

ggsave("FocalActivityBudgets.tif", plot = last_plot(), device = "tiff", width = 85, height = 85, units = "mm", dpi = 300)

rm(r,t,x,y,arc.df,significance.stars)


####################################################
## Behavioral event analysis
####################################################

## Load juvenile behavioral event data.
jevents <- read.csv("juvenile_events.csv")

## For each focal juvenile, determine rates of agg/aff/sex behaviors.
event_rates <-
  jevents %>%
  count(dolphin_id, sex, duration, num_follows, behavior_category) %>%
  dcast(., dolphin_id + sex + duration + num_follows ~ behavior_category, value.var = "n") %>%
  replace(is.na(.), 0) %>%
  filter(duration >= 180) %>%
  filter(num_follows >= 2) %>%
  mutate(AffRate = Aff/duration) %>%
  mutate(AggRate = Agg/duration) %>%
  mutate(SexRate = Sex/duration)

## Do rates of aggressive events vary by sex?
permTS(event_rates[event_rates$sex == "MALE",]$AggRate,
       event_rates[event_rates$sex == "FEMALE",]$AggRate,
       alternative = "two.sided",
       method = "exact.network")

## Do rates of sociosexual events vary by sex?
permTS(event_rates[event_rates$sex == "MALE",]$SexRate,
       event_rates[event_rates$sex == "FEMALE",]$SexRate,
       alternative = "two.sided",
       method = "exact.network")

## Do rates of affiliative events vary by sex?
permTS(event_rates[event_rates$sex == "MALE",]$AffRate,
       event_rates[event_rates$sex == "FEMALE",]$AffRate,
       alternative = "two.sided",
       method = "exact.network")

## Summary stats
ratesummary <- function(x){return(c(median(x), IQR(x)))}
ratesummary(event_rates[event_rates$sex == "MALE",]$AggRate*60) 
ratesummary(event_rates[event_rates$sex == "FEMALE",]$AggRate*60)
ratesummary(event_rates[event_rates$sex == "MALE",]$AffRate*60) 
ratesummary(event_rates[event_rates$sex == "FEMALE",]$AffRate*60)
ratesummary(event_rates[event_rates$sex == "MALE",]$SexRate*60) 
ratesummary(event_rates[event_rates$sex == "FEMALE",]$SexRate*60)


