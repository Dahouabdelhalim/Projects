##Morphological analyses

lake_data <- read.csv("morphology_lake_data.csv", sep = ";")

#Plots and linear models

#Benthic diet

#Linear model only for solitary populations
plot(lake_data$benthicdiet_score ~ lake_data$log_surface_area, col = lake_data$color, pch = lake_data$shape)
abline(lm(lake_data$benthicdiet_score[1:21] ~ lake_data$log_surface_area[1:21]))
summary(lm(lake_data$benthicdiet_score[1:21] ~ lake_data$log_surface_area[1:21]))

#Linear model including sympatric species
plot(lake_data$benthicdiet_score ~ lake_data$log_surface_area, col = lake_data$color, pch = lake_data$shape)
abline(lm(lake_data$benthicdiet_score ~ lake_data$log_surface_area))
summary(lm(lake_data$benthicdiet_score ~ lake_data$log_surface_area))

#Adjust p values for multiple comparisons using Bonferroni correction
p.bonferroni <- p.adjust(c(1.70e-05, 0.00254), method = "bonferroni", n = 2)

#Residuals
diet_surface.lm <- lm(lake_data$benthicdiet_score ~ lake_data$log_surface_area)
diet_surface.res <- resid(diet_surface.lm)
lake_data$diet_res <- diet_surface.res

plot(lake_data$diet_res ~ lake_data$log_surface_area, col = lake_data$color, pch = lake_data$shape)
abline(h=0)

#Gill raker numbers

#Linear model only for solitary populations
plot(lake_data$gill_raker_num ~ lake_data$log_surface_area, col = lake_data$color, pch = lake_data$shape)
abline(lm(lake_data$gill_raker_num[1:21] ~ lake_data$log_surface_area[1:21]))
summary(lm(lake_data$gill_raker_num[1:21] ~ lake_data$log_surface_area[1:21]))

#Linear model including sympatric species
plot(lake_data$gill_raker_num ~ lake_data$log_surface_area, col = lake_data$color, pch = lake_data$shape)
abline(lm(lake_data$gill_raker_num ~ lake_data$log_surface_area))
summary(lm(lake_data$gill_raker_num ~ lake_data$log_surface_area))

#Adjust p values for multiple comparisons using Bonferroni correction
p.bonferroni <- p.adjust(c(0.00565, 0.03607), method = "bonferroni", n = 2)

#Residuals
gill_raker_surface.lm <- lm(lake_data$gill_raker_num ~ lake_data$log_surface_area)
gill_raker_surface.res <- resid(gill_raker_surface.lm)
lake_data$gill_raker_res <- gill_raker_surface.res

plot(lake_data$gill_raker_res ~ lake_data$log_surface_area, col = lake_data$color, pch = lake_data$shape)
abline(h=0)
