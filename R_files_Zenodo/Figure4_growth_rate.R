### Analyze relative growth rate in relation to seed size 

library(lmerTest)
library(ggplot2)
library(ggeffects)


### 00. import data ------------------------------------------------------------
setwd("~/Desktop/Dryad/")

dat <- read.csv(
    file = "Plantago_seed_weights_and_growth_rates.csv" ,
    header = TRUE,
    stringsAsFactors = FALSE)

# add column for age on July 14th (Julian 195)
dat$age_Julian195 <- 195 - dat$Jdate_germ

# remove rows for which seed didn't germinate
dat <-  dat[!is.na(dat$Height_July_14_mm), ]

# isolate cohort of seeds that germinated the same day
dat <- subset(dat, Jdate_germ == 184)



### 02. model growth rate as a function of seed size ---------------------------

# Growth Rate (RGR)
    # RGR = log(X2 − X1)/(t2 − t1), 
        # where X1 = height at time 1
        # where X2 = height at time 2
        # where t1 = Julian 177
        # where t2 = Julian 195
# RGR = log(Plant Height) / time2 - time1

# remove plants that have height = 0
dat2 <- subset(dat, Height_July_14_mm > 0)

dat2$days_old <- 195 - dat2$Jdate_germ

# calculate RGR
dat2$RGR <- log(dat2$Height_July_14_mm) / (195 - dat2$Jdate_germ)

### model relative growth rate of all seeds ------------------------------------
m1 <- lmer(RGR ~ seed_weight_mg + 
        (1|Population) + (1|Maternal_line) + (1|Plate), data = dat2)
summary(m1)


# prediction data.frame
df1 <- ggpredict(m1, terms = c("seed_weight_mg[all]")) # constrains prediction interval

ggplot(df1, aes(x, predicted)) + 
      geom_point(data = dat2, aes(x = seed_weight_mg, y = RGR), size = 3, pch = 16, col = "darkgray", alpha = 0.75) +
      geom_line(aes(linetype=group, color=group), size = 1) +
      geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha = 0.25) +
      scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
      scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5), labels = c(0.25, 0.5, 0.75, 1, 1.25, 1.5), limits = c(0.25, 1.5)) +
      theme_light() +
      ylim(0, 0.3) +
      xlab("Seed mass (mg)") +
      theme(axis.title.x = element_text(size = 12, face = "bold")) +
      ylab("Growth rate (cohort of 62 seeds)") +
      theme(axis.title.y = element_text(size = 12, face = "bold")) +
      theme(legend.title = element_blank()) +
      theme(legend.position = "none") 

#quartz.save(file = "Figure4_Growth_Rate_07_27_2022.jpg", type = "jpg", dpi = 300)



