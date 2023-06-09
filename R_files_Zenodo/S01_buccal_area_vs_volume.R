# Comparing measurements of buccal cavity sagittal section area (the main
# dataset) to actual buccal cavity volume, as measured on micro-CT scans (of
# which we only have 20)

# bookkeeping ####
data <- read.csv("data/ct_scan_measurements.csv")

# convert measurements in voxels into mm^3 and mm^2 based on voxel size:
data$volume <- data$Endocast_volume_voxels * 
  data$Voxel.size..mm.^3
data$area <- data$Sagittal_section_voxels * 
  data$Voxel.size..mm.^2
cols <- c("#fccb4e", "#BFBFBF",
          "#39d69a", "#3e4fbc")

# model fit ####

# so that both are in the same units (instead of regressing mm^3 on mm^2):
v <- data$volume^(1/3)
a <- data$area^(1/2)

# fit linear model:
pfit <- lm(v ~ a)
summary(pfit)

# plotting results ####

# 1. Raw data
plot(data$area, data$volume,
     col = cols[factor(data$behavior)],
     pch = 19,
     ylab = bquote("Volume"~(mm^3)),
     xlab = bquote("Surface area"~(mm^2)),
     main = "Sagittal section vs. total buccal cavity volume\\n(raw data)")
legend(0, 800, 
       legend = c("Mouthbrooding only",
                  "Winnowing only",
                  "Both behaviors",
                  "Neither behavior"),
       fill = cols[c(1, 4, 3, 2)], cex = 0.7)

# 2. Linear fit
plot(a, v,
     col = cols[factor(data$behavior)],
     pch = 19,
     ylab = bquote(sqrt("Volume"~(mm^3), 3)),
     xlab = bquote(sqrt("Surface area"~(mm^2))),
     main = "Sagittal section vs. total buccal cavity volume")
abline(pfit, lty = 2, lwd = 2, col = "red")
text(x = 7, y = 7.2, labels = bquote(R^2~"= 0.96"))
text(x = 7, y = 6.6, labels = bquote(italic(p)~"= 2.67e-14"))

# 3. Residuals 
plot(y = (pfit$residuals)^3, 
     x = data$Standard_length_mm, 
     xlab = "Standard length (mm)",
     ylab = bquote("Volume residuals"~(mm^3)),
     main = "Buccal cavity residuals vs. standard length",
     pch = 19, 
     col = cols[factor(data$behavior)], 
     panel.first = abline(h = 0, lty = 2, lwd = 2, col = grey(0.3))) 
arrows(40, -0.05, 40, -0.25, length = 0.1, code = 2)
arrows(40, 0.05, 40, 0.25, length = 0.1, code = 2)
text(40, 0.23, pos = 4,
     label = "Volume underestimated\\nby area", cex = 0.7)
text(40, -0.23, pos = 4,
     label = "Volume overestimated\\nby area", cex = 0.7)

mean_pct_error <- round(mean(abs(pfit$residuals^3) / v^3) * 100, digits = 2)
text(110, 0.3, 
     label = paste0("Avg. error:\\n", 
                    mean_pct_error, "% of actual volume"),
     cex = 0.7)
