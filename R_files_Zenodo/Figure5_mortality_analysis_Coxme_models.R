### Analyze mortality data from drought experiment -----------------------------

### 00. import mortality data --------------------------------------------------
setwd("~/Desktop/Dryad/")

dat <- read.csv(
      file = "mortality_data.csv" ,
      header = TRUE,
      stringsAsFactors = FALSE)

# remove doubled planted and unplanted cells
dat <- subset(dat, Remove != "Remove")

# add column for days to death (stopped watering on June 28th = Julian 179)
dat$days_to_death <- dat$Jdate_dead - 179

# remove resurrection populations
dat2 <-  dat[!grepl(pattern = "herb", x = dat$Population),]
nrow(dat2)


### 01. import climate data ----------------------------------------------------
clim_dat <- read.csv(file = "Plantago_climate_data_simple.csv", header = T, stringsAsFactors = F)

# merge data
mort_dat <- merge(dat2, clim_dat, by.x = "Population", by.y = "Population")


### 02. conduct survival analysis ----------------------------------------------
library(coxme)
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(coxme)
library(car)


# fit mixed effects Cox survival model
cox_me <- coxme(Surv(days_to_death) ~ 
                   plant_height_mm_on_06_29 + AI + Jdate_germ +
                   (1|Population) + (1|Envelope_ID) + (1|tray), 
                data = mort_dat)
summary(cox_me)
Anova(cox_me)


# remove plant height from model
cox_me_test <- coxme(Surv(days_to_death) ~ 
                   AI + Jdate_germ +
                   (1|Population) + (1|Envelope_ID) + (1|tray), 
                data = mort_dat)
summary(cox_me_test)
Anova(cox_me_test)


### 03. plot the data ----------------------------------------------------------

# add categories for plant height
range(mort_dat$plant_height_mm_on_06_29, na.rm = T)


# 3 categories -----------------------------------------------------------------
mort_dat$height_cat3 <-  
   ifelse(mort_dat$plant_height_mm_on_06_29 < 5, "<5mm",
      ifelse(mort_dat$plant_height_mm_on_06_29 < 10, "5-10mm", ">10mm"))

height_fit3 <- survfit(Surv(days_to_death) ~ height_cat3, data = mort_dat)

# convert prediction intervals to data.frame
hf3_df <- fortify(height_fit3)
names(hf3_df)[names(hf3_df) == "strata"] <- "Height"

hf3_df$Height = factor(hf3_df$Height, levels = c("<5mm", "5-10mm", ">10mm"))

# plot data
ggplot(data = hf3_df, aes(x = time, y = surv, color = Height)) +
   geom_line(size = 1.25) +
   geom_ribbon(aes(ymin = lower, ymax = upper, fill = Height), linetype = 0, alpha = 0.25) +
   theme_light() +
   xlim(5,25) +
   theme(axis.title.x = element_text(size = 12, face = "bold")) +
   theme(axis.title.y = element_text(size = 12, face = "bold")) +
   xlab("Days of drought") +
   ylab("Survival") +
   theme(legend.position = c(0.15,0.25)) +
   theme(legend.key.size = unit(1.25, "cm")) +
   theme(legend.title = element_text(size = 12, face = "bold")) +
   theme(legend.text = element_text(size = 10)) +
   theme(legend.title.align=0.5) 

# quartz.save(file = "Figure5_survival_curves_01_26_22.jpg", type = "jpg", dpi = 300)


### 04. conduct survival analysis for plants with SLA data ---------------------
sla_dat <- read.csv(file = "SLA_data.csv", header = T, stringsAsFactors = F)

# add SLA column
sla_dat$SLA <- sla_dat$leaf_area_mm_squared / sla_dat$leaf_weight_mg

# remove records without SLA data
sla_dat <- sla_dat[!is.na(sla_dat$SLA),]

# remove resurrection populations
sla_dat <- sla_dat[!grepl(pattern = "herb", x = sla_dat$Population),]

# merge SLA data and climate data
sla_dat2 <- merge(sla_dat, clim_dat, by.x = "Population", by.y = "Population")

# add column for days to death (stopped watering on June 28th = Julian 179)
sla_dat2$days_to_death <- sla_dat2$Jdate_dead - 179

# conduct survival analysis
cox_me_SLA <- coxme(Surv(days_to_death) ~ 
                     plant_height_mm_on_06_29 + AI + Jdate_germ + SLA +
                     (1|Population) + (1|Envelope_ID) + (1|tray), 
                     data = sla_dat2)

summary(cox_me_SLA)
Anova(cox_me_SLA)




