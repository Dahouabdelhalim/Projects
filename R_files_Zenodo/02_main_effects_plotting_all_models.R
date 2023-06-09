
library(tidyverse)
library(patchwork)
library(visreg)
library(ggplot2)

load(file = "Data/Diuron_photo_model_and_data2.RDA") # Diuron photosynthesis models and data
load(file = "Data/DIN_photo_model_and_data2.RDA") # DIN photosynthesis models and data
load(file = "Data/Diuron_growth_model_and_data2.RDA") # Diuron growth models and data
load(file = "Data/DIN_growth_model_and_data2.RDA") # Diuron growth models and data


# MAIN EFFECTS PLOTS ------------------------------------------------------

# DIN (growth) -----------------------------------------------------

# Makes dataframe of all combos of the main variables of interest - DIN
preddatDIN <- expand.grid(hours = unique(DIN_growth_data$hours),
                          DIN_num = unique(DIN_growth_data$DIN_num),
                          Light_num = unique(DIN_growth_data$Light_num),
                          t0 = 47,
                          block = 0)

# note: also need hours as a factor, so make that: 
preddatDIN$hours_fact <- factor(preddatDIN$hours)

# Predict for variables in above dataframe, ignoring the block and sample ID random effects
m1fit <- predict(m1_din_light, newdata = preddatDIN, type = "response",
                 exclude = c("block"),
                 se.fit = TRUE)
preddatDIN$fit <- m1fit$fit
preddatDIN$sefit <- m1fit$se.fit

#Olivia, to plot just select DIN values, filter the dataframe like this:
preddatDIN_lim <- filter(preddatDIN, DIN_num %in% c(0, 2.67, 26.7))
# %in% means belongs to, so it will select rows that match either
# 2.67 or 26.7

#Note use \\n to add a new line in a label
leg_lab <- expression("Light level \\n (μmol photons m"^-2*"s"^-1*")")

#Now use that to plot: 
fig1 <- ggplot(preddatDIN_lim) + 
  aes(x = hours, y = fit, color = Light_num,
      fill = Light_num,
      group = Light_num) + 
  geom_ribbon(aes(ymin = fit -sefit*1.96, ymax = fit+sefit*1.96),
              alpha = 0.5, color = NA) +
  geom_line() +
  facet_grid(.~DIN_num) +
  theme_classic()  +
  scale_color_manual(values = c("midnightblue","blue","skyblue")) +
  scale_fill_manual(values = c("midnightblue","blue","skyblue")) +
  aes(x = hours, y = fit, color = factor(Light_num),
      fill = factor(Light_num),
      group = Light_num) +
  labs(x = expression(Hours), y = expression(Log~Cell~density~(OD[685]))) + #changes the x and y axis titles
  guides(color=guide_legend(leg_lab),
         fill=guide_legend(leg_lab)) + # changes the legend title
  labs(title = "A") + #add A and B labels +
  theme(legend.position = "none") # removes legend, so A part of the figure doesn't have a legend

fig1 <- fig1 + geom_hline( yintercept = unique(log(preddatDIN_lim$t0)), linetype = "dotted")
  
# # save the plot into output folder
# ggsave("Outputs/maineffects_DINgrowth_light_IR_limited_2021-12-03.png", 
#        width = 30, height = 13, units = "cm", dpi = 300 )

# DIURON (growth) -----------------------------------------------------

# Makes dataframe of all combos of the main variables of interest - DIURON
preddatDIURON <- expand.grid(hours = unique(Diuron_growth_data$hours),
                             Diuron_num = unique(Diuron_growth_data$Diuron_num),
                             Light_num = unique(Diuron_growth_data$Light_num),
                             t0 = 47,
                             block = 0)

# note: also need hours as a factor, so make that: 
preddatDIURON$hours_fact <- factor(preddatDIURON$hours)

# Predict for variables in above dataframe, ignoring the block and sample ID random effects
m2fit <- predict(m1_diuron_light, newdata = preddatDIURON, type = "response",
                 exclude = c("block"),
                 se.fit = TRUE)
preddatDIURON$fit <- m2fit$fit
preddatDIURON$sefit <- m2fit$se.fit

#to plot only select diuron values, filter the dataframe like this:
preddatDiuron_lim <- filter(preddatDIURON, Diuron_num %in% c(0, 0.11, 0.33, 1, 3))
# %in% means belongs to, so it will select rows that match either
# 0.11, 0.33, 1, 3

#Note use \\n to add a new line in a label
leg_lab <- expression("Light level \\n (μmol photons m"^-2*"s"^-1*")")

#Now use that to plot: 
fig2 <- ggplot(preddatDiuron_lim) + 
  aes(x = hours, y = fit, color = Light_num,
      fill = Light_num,
      group = Light_num) + 
  geom_ribbon(aes(ymin = fit -sefit*1.96, ymax = fit+sefit*1.96),
              alpha = 0.5, color = NA) +
  geom_line() +
  facet_grid(.~Diuron_num) + 
  theme_classic() +
  scale_color_manual(values = c("midnightblue","blue","skyblue")) + #changes the colour of the lines
  scale_fill_manual(values = c("midnightblue","blue","skyblue")) + #changes the colour of standard error shading
  aes(x = hours, y = fit, color = factor(Light_num),
      fill = factor(Light_num),
      group = Light_num) +
  labs(x = expression(Hours), y = expression(Log~Cell~density~(OD[685]))) + #changes the x and y axis titles
  guides(color=guide_legend(leg_lab),
         fill=guide_legend(leg_lab)) + # changes the legend title
  labs(title = "A") + #add A and B labels
  theme(legend.position = "none") # removes legend, so A part of the figure doesn't have a legend

fig2 <- fig2 + geom_hline( yintercept = unique(log(preddatDiuron_lim$t0)), linetype = "dotted")

# # save the plot into output folder
# ggsave("Outputs/maineffects_Diurongrowth_light_IR_limited_2021-12-03.png", 
#        width = 30, height = 13, units = "cm", dpi = 300 )


# DIN (photosynthesis) -----------------------------------------------------

# Makes dataframe of all combos of the main variables of interest - DIN
preddatDINphoto <- expand.grid(hours = unique(DIN_photo_data$hours),
                               DIN_num = unique(DIN_photo_data$DIN_num),
                               Light_num = unique(DIN_photo_data$Light_num),
                               t0 = 0.52,
                               block = 0)

# note: also need hours as a factor, so make that: 
preddatDINphoto$hours_fact <- factor(preddatDINphoto$hours)

# Predict for variables in above dataframe, ignoring the block and sample ID random effects
m3fit <- predict(m1_DIN_yield, newdata = preddatDINphoto, type = "response",
                 exclude = c("block"),
                 se.fit = TRUE)
preddatDINphoto$fit <- m3fit$fit
preddatDINphoto$sefit <- m3fit$se.fit

preddatDINphoto_lim <- filter(preddatDINphoto, DIN_num %in% c(0, 2.67, 26.7))

#Note use \\n to add a new line in a label
leg_lab <- expression("Light level \\n (μmol photons m"^-2*"s"^-1*")")

#Now use that to plot: 
fig3 <- ggplot(preddatDINphoto_lim) + 
  aes(x = hours, y = fit, color = Light_num,
      fill = Light_num,
      group = Light_num) + 
  geom_ribbon(aes(ymin = fit -sefit*1.96, ymax = fit+sefit*1.96),
              alpha = 0.5, color = NA) +
  geom_line() +
  facet_grid(.~DIN_num) +
  theme_classic()  +
  scale_color_manual(values = c("midnightblue","blue","skyblue")) +
  scale_fill_manual(values = c("midnightblue","blue","skyblue")) +
  aes(x = hours, y = fit, color = factor(Light_num),
      fill = factor(Light_num),
      group = Light_num) +
  xlab("Hours") + ylab("Chlorophyll-a fluorescence") + #changes the x and y axis titles
  guides(color=guide_legend(leg_lab),
         fill=guide_legend(leg_lab)) + # changes the legend title
  labs(title = "B") #add A and B labels

fig3 <- fig3 + geom_hline( yintercept = unique(preddatDINphoto_lim$t0), linetype = "dotted")

# # save the plot into output folder
# ggsave("Outputs/maineffects_DINphoto_light_IR_limited_2021-12-03.png", 
#        width = 30, height = 13, units = "cm", dpi = 300 )

# DIURON (photosynthesis) -----------------------------------------------------

# Makes dataframe of all combos of the main variables of interest - DIURON
preddatDIURONphoto <- expand.grid(hours = unique(Diuron_photo_data$hours),
                                  Diuron_num = unique(Diuron_photo_data$Diuron_num),
                                  Light_num = unique(Diuron_photo_data$Light_num),
                                  t0 = 0.52,
                                  block = 0)

# note: also need hours as a factor, so make that: 
preddatDIURONphoto$hours_fact <- factor(preddatDIURONphoto$hours)

# Predict for variables in above dataframe, ignoring the block and sample ID random effects
m4fit <- predict(m1_diuron_yield, newdata = preddatDIURONphoto, type = "response",
                 exclude = c("block"),
                 se.fit = TRUE)
preddatDIURONphoto$fit <- m4fit$fit
preddatDIURONphoto$sefit <- m4fit$se.fit

#to plot only select diuron values, filter the dataframe like this:
preddatDiuronphoto_lim <- filter(preddatDIURONphoto, Diuron_num %in% c(0, 0.11, 0.33, 1, 3))
# %in% means belongs to, so it will select rows that match either
# 0.11, 0.33, 1, 3

#Note use \\n to add a new line in a label
leg_lab <- expression("Light level \\n (μmol photons m"^-2*"s"^-1*")")

#Now use that to plot: 
fig4 <- ggplot(preddatDiuronphoto_lim) + 
  aes(x = hours, y = fit, color = Light_num,
      fill = Light_num,
      group = Light_num) + 
  geom_ribbon(aes(ymin = fit -sefit*1.96, ymax = fit+sefit*1.96),
              alpha = 0.5, color = NA) +
  geom_line() +
  facet_grid(.~Diuron_num) + 
  theme_classic() +
  scale_color_manual(values = c("midnightblue","blue","skyblue")) + #changes the colour of the lines
  scale_fill_manual(values = c("midnightblue","blue","skyblue")) + #changes the colour of standard error shading
  aes(x = hours, y = fit, color = factor(Light_num),
      fill = factor(Light_num),
      group = Light_num) +
  xlab("Hours") + ylab("Chlorophyll-a fluorescence") + #changes the x and y axis titles
  guides(color=guide_legend(leg_lab),
         fill=guide_legend(leg_lab)) + # changes the legend title
  labs(title = "B") #add A and B labels

fig4 <- fig4 + geom_hline( yintercept = unique(preddatDiuronphoto_lim$t0), linetype = "dotted")

# save the plot into output folder
# ggsave("Outputs/maineffects_Diuronphoto_light_IR_limited_2021-12-03.png", 
#        width = 30, height = 13, units = "cm", dpi = 300 )

# PATCHWORK TOGETHER... ---------------------------------------------------

pw1 <- (fig1 + fig3) + #two photosynthesis figures for DIN
  plot_layout(ncol = 2, widths = c(1,1)) & 
  theme(plot.tag = element_text(size = 16))
pw1

ggsave("Outputs/maineffects_DIN_GrowthPhoto_2022-02-09.png",
       plot = pw1, 
       width = 30, height = 13, units = "cm", dpi = 300 )

pw2 <-(fig2 + fig4) #two growth figures for DIURON
plot_layout(ncol = 2, widths = c(1,1)) & 
  theme(plot.tag = element_text(size = 16))
pw2

ggsave("Outputs/maineffects_DIURON_GrowthPhoto_2021-12-03.png", 
       plot = pw2, 
       width = 30, height = 13, units = "cm", dpi = 300 )

