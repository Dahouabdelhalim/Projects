## Juvenile fission-fusion rate
## Ali Galezo

options(stringsAsFactors = FALSE)
library(dplyr)

## Load list of juvenile focal follows
jfollows <- read.csv("juvenile_joinleaves.csv")

## Determine hourly fission fusion rate for each juvenile focal.
jfollows$rate <- jfollows$n_fissionfusion_events/jfollows$cumulative_duration*60

## Only include focals with at least 30 minutes of observation.
jfollows <- jfollows %>% filter(cumulative_duration >= 30)

## Mean, SD, SE juvenile fission-fusion rate
jfollows %>%
  summarize(n_juveniles = n(),
            mean = mean(rate),
            median = median(rate),
            sd = sd(rate),
            se = sd/sqrt(n_juveniles),
            every_x_minutes = 60/mean)
