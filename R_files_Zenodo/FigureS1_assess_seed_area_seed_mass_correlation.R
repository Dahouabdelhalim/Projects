### Assess correlation between seed mass and seed area


### 01 import data -------------------------------------------------------------
setwd("~/Desktop/Dryad/")
dat <- read.csv(file = "S1_Plantago_seed_weights.csv", header = T, stringsAsFactors = F)

### 02. assess simple correlation ----------------------------------------------
plot(seed_weight_mg ~ area, data = dat,
     pch = 21, bg = "lightgray",
     xlab = "Seed Area (mm^2)",
     ylab = "Seed mass (mg)")

cor.test(dat$seed_weight_mg, dat$area)


### 03. build linear model to account for population ---------------------------
library(lmerTest)
library(ggplot2)
library(ggeffects)

mm1 <- lmer(seed_weight_mg ~ area + (1|Population), data = dat)
summary(mm1)

# calculate R2 for linear mixed model
library(partR2)

partR2(mm1, data = dat, R2_type = "marginal", partvars = c("area"), nboot = 1000)


# build prediction data.frame
df <- ggpredict(mm1, terms = c("area"))

ggplot(df, aes(x, predicted)) + 
      geom_point(data = dat, aes(x = area, y = seed_weight_mg), size = 3, pch = 16, col = "darkgray", alpha = 0.5) +
      geom_line(aes(linetype=group, color=group), size = 1, col = "blue") +
      geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha = 0.5, fill = "blue") +
      scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
      theme_light() +
      xlab(expression(bold("Seed area (mm"^2*")"))) +
      theme(axis.title.x = element_text(size = 12, face = "bold")) +
      ylab("Seed mass (mg)") +
      theme(axis.title.y = element_text(size = 12, face = "bold")) +
      theme(legend.title = element_blank()) +
      theme(legend.position = "none") 

# quartz.save(file = "FigureS1_seed_area_mass_correlation.jpg", type = "jpg", dpi = 300)
