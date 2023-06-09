# ***************************************************************************************
# GSEARIOModel: Version 2.0 by Yaxin Zhang
# GSEARIOModel is extended based on Wang,D (2020),DaopingW/economic-impact-model: Disaster Footprint Model (Version 1.0) [Source code].https://doi.org/10.5281/zenodo.4290117.
# ***************************************************************************************

# Part4_2: Plot
# this file is used to plot figures
# called from GSEARIOModel_Main
# coding with UTF-8

# Set working directory -------------------------------------------------------------

rootpathway <- "D:/Workspace/GSEARIOModel_v2.0"
setwd(rootpathway)

# Load packages -----------------------------------------------------------

library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggsci)
library(ggThemeAssist)
library(ggalluvial)
library(gtable)
library(grid)
library(readxl)
library(gg.gap)
library(ggplot2)
library(grid)
library(cowplot)
library(showtext)
library(hrbrthemes)
library(ggbreak)
library(gg.gap)

# Figure requirement -------------------------------------------------------
wordsize <- 6
linesize <- 0.4
pointsize <- 1.5
errlinesize <- 0.15
dpisize <- 1200

# P1 The changes of economic output and labor demand under the COVID-19 pandemic ----------------------------------------------------

np1 <- pal_npg("nrc")(5)
show_col(np1)
np2 <- c("grey", "grey")
np3 <- c("grey", "grey")

DataforPlotTMC <- function(df) {
  df1 <- data.table(df)
  df1$Week <- as.Date(df1$Week)
  df1$Min <- as.numeric(df1$Min)
  df1$Max <- as.numeric(df1$Max)
  df1$Mea <- as.numeric(df1$Mea)
  df1[Sce == "L22", `:=`(Sce, "ELS")]
  df1[Sce == "L30", `:=`(Sce, "LS")]
  df1[Sce == "T30", `:=`(Sce, "TS")]
  df1[Sce == "R30", `:=`(Sce, "CS")]
  df1[Sce == "B22", `:=`(Sce, "BAU")]
  df1[Sce == "N22", `:=`(Sce, "NS")]
  df1[Sce == "B30", `:=`(Sce, "BAU")]
  df1[Sce == "N30", `:=`(Sce, "NS")]
  df1$Sce <- factor(df1$Sce, levels = c("ELS", "TS",
                                        "NS", "CS", "LS", "BAU"), ordered = T)
  return(df1)
}

df <- data.table(df_IOX_Value_sum_plot)
df <- df[Sce != "L31" & Sce != "L32" & Sce != "T31" &
           Sce != "T32" & Sce != "R31" & Sce != "R32" & Sce !=
           "R22" & Sce != "T22"]
df1 <- DataforPlotTMC(df)
df1$Min <- df1$Min/1e+06
df1$Max <- df1$Max/1e+06
df1$Mea <- df1$Mea/1e+06
df1 <- df1[c("BAU", "NS"), on = "Sce"]
P1a <- ggplot() + theme_bw() + geom_ribbon(data = df1[Cty ==
                                                        "GB"], aes(x = Week, ymin = Min, ymax = Max, fill = Sce),
                                           alpha = 0.9) + geom_line(data = df1[Cty == "GB" &
                                                                                 Sce == "BAU"], aes(Week, Min), color = "black", size = linesize,
                                                                    alpha = 0.8) + geom_line(data = df1[Cty == "GB" &
                                                                                                          Sce == "NS"], aes(Week, Mea), color = "grey50", size = linesize,
                                                                                             alpha = 0.8) + scale_fill_manual(values = np2) + theme(panel.grid.major.y = element_blank(),
                                                                                                                                                    panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
                                                                                                                                                    panel.grid.minor.x = element_blank(), axis.text.y = element_text(angle = 0,
                                                                                                                                                                                                                     size = wordsize, colour = "black"), axis.text.x = element_text(angle = 0,
                                                                                                                                                                                                                                                                                    size = wordsize, colour = "black"), axis.title.y = element_text(size = wordsize,
                                                                                                                                                                                                                                                                                                                                                    color = "black"), axis.title.x = element_text(size = wordsize,
                                                                                                                                                                                                                                                                                                                                                                                                  color = "black"), axis.text.y.right = element_blank(),
                                                                                                                                                    axis.ticks.y.right = element_blank(), legend.position = "none",
                                                                                                                                                    panel.background = element_rect(colour = "black")) +
  ylab("Changes of economic output\\n(tillion EURO)") +
  xlab("") + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max(df1$Max) * 1.05),
                     expand = c(0, 0), breaks = seq(0, max(df1$Max) *
                                                      1.05, 1), labels = seq(0, max(df1$Max) * 1.05,
                                                                             1))

df <- data.table(df_JB_Value_sum_plot)
df <- df[Sce != "L31" & Sce != "L32" & Sce != "T31" &
           Sce != "T32" & Sce != "R31" & Sce != "R32" & Sce !=
           "R22" & Sce != "T22"]
df1 <- DataforPlotTMC(df)
df1$Min <- df1$Min/1000
df1$Max <- df1$Max/1000
df1$Mea <- df1$Mea/1000
df1 <- df1[c("BAU", "NS"), on = "Sce"]
P1b <- ggplot() + theme_bw() + geom_ribbon(data = df1[Cty ==
                                                        "GB"], aes(x = Week, ymin = Min, ymax = Max, fill = Sce),
                                           alpha = 0.9) + geom_line(data = df1[Cty == "GB" &
                                                                                 Sce == "BAU"], aes(Week, Min), color = "black", size = linesize,
                                                                    alpha = 0.8) + geom_line(data = df1[Cty == "GB" &
                                                                                                          Sce == "NS"], aes(Week, Mea), color = "grey50", size = linesize,
                                                                                             alpha = 0.8) + scale_fill_manual(values = np2) + theme(panel.grid.major.y = element_blank(),
                                                                                                                                                    panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
                                                                                                                                                    panel.grid.minor.x = element_blank(), axis.text.y = element_text(angle = 0,
                                                                                                                                                                                                                     size = wordsize, colour = "black"), axis.text.x = element_text(angle = 0,
                                                                                                                                                                                                                                                                                    size = wordsize, colour = "black"), axis.title.y = element_text(size = wordsize,
                                                                                                                                                                                                                                                                                                                                                    color = "black"), axis.title.x = element_text(size = wordsize,
                                                                                                                                                                                                                                                                                                                                                                                                  color = "black"), axis.text.y.right = element_blank(),
                                                                                                                                                    axis.ticks.y.right = element_blank(), legend.position = "none",
                                                                                                                                                    panel.background = element_rect(colour = "black")) +
  ylab("Changes of labor demand loss\\n(million people)") +
  xlab("") + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max(df1$Max) * 1.05),
                     expand = c(0, 0), breaks = seq(0, max(df1$Max) *
                                                      1.05, 15), labels = seq(0, max(df1$Max) *
                                                                                1.05, 15))

CtyCode1 <- read_excel("./Data/parameter_data_v6_0527.xlsx",
                       sheet = "cty", "AC71:AC93", col_names = F)
CtyCodeOrder1 <- unlist(read_excel("./Data/parameter_data_v6_0527.xlsx",
                                   sheet = "cty", "AH71:AH83", col_names = F))
CtyCodeOrder2 <- unlist(read_excel("./Data/parameter_data_v6_0527.xlsx",
                                   sheet = "cty", "AK71:AK83", col_names = F))


pal3 <- c("#787ebc", "#ef7a7d", "#fbe135", "#add5b0",
          "#71b62b")
mycolorpie <- colorRampPalette(colors = pal3, interpolate = "line")(15)
pal1 <- colorRampPalette(colors = c("#113260", "#2D9EE0"),
                         interpolate = "line")(2)
pal2 <- colorRampPalette(colors = c("#f47834", "#fbe135"),
                         interpolate = "line")(3)
pal3 <- colorRampPalette(colors = c("#b31217", "#f03e3e"),
                         interpolate = "line")(5)
pal4 <- colorRampPalette(colors = c("#6a3093", "#a044ff"),
                         interpolate = "line")(2)
pal <- c(pal1, pal2, pal3, pal4)
show_col(pal)

P1e_1 <- ggplot() + geom_bar(data = df11, mapping = aes(x = "Cty",
                                                        y = Value, fill = Cty), alpha = 0.9, stat = "identity",
                             position = "stack", colour = "white", lwd = linesize/1.5) +
  coord_polar(theta = "y") + labs(x = "", y = "", title = "") +
  theme_bw() + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.background = element_rect(fill = "transparent",
                                                     colour = NA), plot.background = element_rect(fill = "transparent",
                                                                                                  colour = NA), plot.title = element_text(colour = "black",
                                                                                                                                          face = "italic", size = 12, vjust = -5, hjust = 0.5),
                     legend.position = "none", panel.border = element_blank()) +
  scale_fill_manual(values = pal)

P1e_2 <- ggplot() + geom_bar(data = df22, mapping = aes(x = "Cty",
                                                        y = Value, fill = Cty), alpha = 0.9, stat = "identity",
                             position = "stack", colour = "white", lwd = linesize/1.5) +
  coord_polar(theta = "y") + labs(x = "", y = "", title = "") +
  theme_bw() + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     panel.background = element_rect(fill = "transparent",
                                                     colour = NA), plot.background = element_rect(fill = "transparent",
                                                                                                  colour = NA), plot.title = element_text(colour = "black",
                                                                                                                                          face = "italic", size = 12, vjust = -5, hjust = 0.5),
                     legend.position = "none", panel.border = element_blank()) +
  scale_fill_manual(values = pal)


CtyCode <- read_excel("./Data/parameter_data_v6_0527.xlsx",
                      sheet = "cty", "AC44:AC94", col_names = F)
CtyCodeOrder <- unlist(read_excel("./Data/parameter_data_v6_0527.xlsx",
                                  sheet = "cty", "AK71:AK83", col_names = F))
CtyCode_1 <- unlist(CtyCode)
CtyCode_2 <- rep(CtyCode_1, each = 3)

df_GB_point <- data.table(errbar_JB_skill)  #point 
df_GB_bar <- data.table(errbar_JB_mean)  #bar
df_GB_point <- df_GB_point[, c("Cty", "Ski", "Value",
                               "ci")]
df_GB_bar <- df_GB_bar[, c("Cty", "Value", "ci")]

df1_point <- JB_Rate_RR_sum_ski_NS_MMM_skip[, c("Cty",
                                                "Ski", "Value", "ci")]
df1_point$Value <- as.numeric(df1_point$Value)
df1_point$ci <- as.numeric(df1_point$ci)
df2_bar <- aggregate(df1_point[, c("Value", "ci")], by = list(df1_point$Cty),
                     mean)
colnames(df2_bar)[1] <- "Cty"

df1_point <- rbind(df1_point, df_GB_point[1:3, ])
df2_bar <- rbind(df2_bar, df_GB_bar[1, ])

CtyCode_bar <- CtyCode_1[-50]
CtyCode_point <- rep(CtyCode_bar, each = 3)

df1_point$CtyCode <- CtyCode_point
df2_bar$CtyCode <- CtyCode_bar

df1_point <- aggregate(df1_point[, c("Value", "ci")],
                       by = list(df1_point$CtyCode, df1_point$Ski), mean)
df2_bar <- aggregate(df2_bar[, c("Value", "ci")], by = list(df2_bar$CtyCode),
                     mean)
colnames(df1_point)[1:2] <- c("Cty", "Ski")
colnames(df2_bar)[1] <- "Cty"

df1_point <- data.table(df1_point)
df2_bar <- data.table(df2_bar)
df1_point <- df1_point[Cty != "Drop"]
df2_bar <- df2_bar[Cty != "Drop"]

df1_point$Cty <- factor(df1_point$Cty, levels = CtyCodeOrder,
                        ordered = T)
df2_bar$Cty <- factor(df2_bar$Cty, levels = CtyCodeOrder,
                      ordered = T)

df1_point <- setorder(df1_point, Cty)
df2_bar <- setorder(df2_bar, Cty)

df1_point <- data.table(df1_point)
df2_bar <- data.table(df2_bar)

df1_point[Ski == "L", `:=`(Ski, "Low-skilled")]
df1_point[Ski == "M", `:=`(Ski, "Medium-skilled")]
df1_point[Ski == "H", `:=`(Ski, "High-skilled")]

df1_point$Value <- df1_point$Value * 100
df1_point$ci <- df1_point$ci * 100
df2_bar$Value <- df2_bar$Value * 100
df2_bar$ci <- df2_bar$ci * 100

df2_bar <- df2_bar[Cty != "The Globe"]
df1_point <- df1_point[Cty != "The Globe"]

df_GB_point <- aggregate(df1_point[, c("Value", "ci")],
                         by = list(df1_point$Ski), mean)
df_GB_bar <- apply(df2_bar[, c("Value", "ci")], 2, mean)

P2_3 <- ggplot() + geom_bar(data = df2_bar[Cty != "The Globe"],
                            aes(Cty, Value), fill = "orange", alpha = 0.7, width = 0.5,
                            stat = "identity", position = position_dodge()) +
  geom_errorbar(data = df2_bar[Cty != "The Globe"],
                aes(Cty, ymin = Value - ci, ymax = Value + ci),
                size = linesize * 0.75, width = errlinesize, color = "black") +
  geom_point(data = df1_point[Cty != "The Globe"], aes(Cty,
                                                       Value, color = Ski), size = pointsize) + geom_errorbar(data = df1_point[Cty !=
                                                                                                                                 "The Globe"], aes(Cty, ymin = Value - ci, ymax = Value +
                                                                                                                                                     ci), size = linesize * 0.75, width = errlinesize,
                                                                                                              color = "grey30") + geom_hline(aes(yintercept = df_GB_bar[1]),
                                                                                                                                             size = linesize, color = "darkblue") + geom_hline(aes(yintercept = df_GB_point[1,
                                                                                                                                                                                                                            2]), linetype = "dashed", size = linesize * 0.75,
                                                                                                                                                                                               color = np1[1]) + geom_hline(aes(yintercept = df_GB_point[2,
                                                                                                                                                                                                                                                         2]), linetype = "dashed", size = linesize * 0.75,
                                                                                                                                                                                                                            color = np1[2]) + geom_hline(aes(yintercept = df_GB_point[3,
                                                                                                                                                                                                                                                                                      2]), linetype = "dashed", size = linesize * 0.75,
                                                                                                                                                                                                                                                         color = np1[3]) + scale_color_npg() + theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.background = element_rect(fill = NA)) +
  coord_flip() + labs(x = "", y = "Labor demand loss decline during the COVID(%)") +
  scale_y_continuous(limits = c(-25.5, 0), expand = c(0,
                                                      0)) + theme(legend.position = "NULL", axis.text.y = element_text(angle = 0,
                                                                                                                       size = wordsize, colour = "black"), axis.text.x = element_text(angle = 0,
                                                                                                                                                                                      size = wordsize, colour = "black"), axis.title = element_text(angle = 0,
                                                                                                                                                                                                                                                    size = wordsize, colour = "black")) + scale_x_discrete(position = "top")

gg_a_b_e <- ggdraw() + draw_plot(P1a, 0, 0.5, 0.6, 0.5) +
  draw_plot(P1b, 0, 0, 0.6, 0.5) + draw_plot(P1e_1,
                                             0.225, 0.51, 0.42, 0.42) + draw_plot(P1e_2, 0.225,
                                                                                  0.01, 0.42, 0.42) + draw_plot(P2_3, 0.6, 0, 0.4, 1) +
  draw_plot_label(c("a", "b", "c"), c(0.07525, 0.085,
                                      0.615), c(0.987, 0.485, 0.987), size = wordsize)

P1 <- gg_a_b_e

setwd(rootpathway)
opfile <- paste("./Plot/", "P1.png", sep = "")
ht <- 11
ggsave(opfile, gg_a_b_e, height = ht, width = 1.2 * ht,
       dpi = dpisize, units = "cm")
file.show(opfile)


# P2 Income inequality change caused by the COVID-19 pandemic as measured by the Gini index and the Theil index. ----------------------------------------------
library(maps)
library(plyr)
library(maptools)
library(readxl)

load("./Output/df_allRR.RData")

P4_Cov_Gini <- ggplot() + geom_polygon(data = df_allRR,
                                       aes(x = long, y = lat, group = group, fill = (Gini_COV)),
                                       colour = "grey40") + scale_fill_gradient2(low = "#3262c3",
                                                                                 mid = "white", high = "OrangeRed", breaks = c(-10,
                                                                                                                               0, 20), limits = c(-10, 20)) + theme(panel.grid = element_blank(),
                                                                                                                                                                    panel.background = element_blank(), axis.text = element_blank(),
                                                                                                                                                                    axis.ticks = element_blank(), axis.title.y = element_blank(),
                                                                                                                                                                    legend.position = "none") + xlab("(a) Gini index") +
  theme(axis.title.x = element_text(size = 12))

P4_Cov_TheilL <- ggplot() + geom_polygon(data = df_allRR,
                                         aes(x = long, y = lat, group = group, fill = TheilL_COV),
                                         colour = "grey40") + scale_fill_gradient2(low = "#3262c3",
                                                                                   mid = "white", high = "OrangeRed", breaks = c(-10,
                                                                                                                                 0, 20), limits = c(-10, 20)) + theme(panel.grid = element_blank(),
                                                                                                                                                                      panel.background = element_blank(), axis.text = element_blank(),
                                                                                                                                                                      axis.ticks = element_blank(), axis.title.y = element_blank(),
                                                                                                                                                                      legend.position = "none") + xlab("(b) Theil index") +
  theme(axis.title.x = element_text(size = 12))

P2 <- P4_Cov_Gini/P4_Cov_TheilL

setwd(rootpathway)
opfile <- paste("./Plot/", "P2.png", sep = "")
ggsave(opfile, P2, height = 9, width = 9, dpi = dpisize,
       units = "in")
file.show(opfile)

# P3 The impact of fiscal stimuli on economic output and labor demand. X-axis is the three fiscal stimulus scenarios. -------------------------------------------------------------------

df1 <- df_IOX_Value_RR_sum_ski_sce_MMM_skip_GB
df1_0 <- df1
df <- data.table(df1)
df[Ski == "L", `:=`(Ski, "Low-skilled")]
df[Ski == "M", `:=`(Ski, "Medium-skilled")]
df[Ski == "H", `:=`(Ski, "High-skilled")]
setorder(df, Ski)
df[Sce == "B22", `:=`(Sce, "BAU")]
df[Sce == "N22", `:=`(Sce, "NS")]
df[Sce == "R30", `:=`(Sce, "CS")]
df[Sce == "T30", `:=`(Sce, "TS")]
df[Sce == "L30", `:=`(Sce, "LS")]
df[Sce == "L22", `:=`(Sce, "ELS")]
df1 <- df[Sce != "ELS" & Sce != "R22"]
df1$Sce <- factor(df1$Sce, levels = c("LS", "TS", "CS"),
                  order = TRUE)

df2 <- df_JB_Value_RR_sum_ski_sce_MMM_skip_GB
df <- data.table(df2)
df$Ski <- factor(df$Ski, levels = c("L", "M", "H"), ordered = T)
df$ci <- as.numeric(df$ci)
df$Value <- as.numeric(df$Value)
df$Value2 <- as.numeric(df$Value2)
df[Ski == "L", `:=`(Ski, "Low-skilled")]
df[Ski == "M", `:=`(Ski, "Medium-skilled")]
df[Ski == "H", `:=`(Ski, "High-skilled")]
setorder(df, Ski)
df[Sce == "B22", `:=`(Sce, "BAU")]
df[Sce == "N22", `:=`(Sce, "NS")]
df[Sce == "R30", `:=`(Sce, "CS")]
df[Sce == "T30", `:=`(Sce, "TS")]
df[Sce == "L30", `:=`(Sce, "LS")]
df[Sce == "L22", `:=`(Sce, "ELS")]
df2 <- df[Sce != "ELS" & Sce != "R22"]
df2$Sce <- factor(df2$Sce, levels = c("LS", "TS", "CS"),
                  order = TRUE)

df1[Sce == "LS", `:=`(Sce, "1")]
df1[Sce == "TS", `:=`(Sce, "2")]
df1[Sce == "CS", `:=`(Sce, "3")]
df1$Sce <- as.numeric(df1$Sce)
df2[Sce == "LS", `:=`(Sce, "1")]
df2[Sce == "TS", `:=`(Sce, "2")]
df2[Sce == "CS", `:=`(Sce, "3")]
df2$Sce <- as.numeric(df2$Sce)

df1$ci <- as.numeric(df1$ci)/1e+06
df1$Value <- as.numeric(df1$Value)/1e+06
df1$Value2 <- as.numeric(df1$Value2)/1e+06
df2$ci <- as.numeric(df2$ci) * 55/1e+06
df2$Value <- as.numeric(df2$Value) * 55/1e+06
df2$Value2 <- as.numeric(df2$Value2) * 55/1e+06


barwidth <- 0.35
windowsFonts(A = windowsFont("Arial"), B = windowsFont("Bookman Old Style"),
             C = windowsFont("Comic Sans MS"))
P3 <- ggplot() + geom_bar(data = df1, mapping = aes(x = Sce,
                                                    y = Value), fill = "orange", alpha = 0.7, width = 0.3,
                          stat = "identity", position = "stack") + geom_bar(data = df2,
                                                                            mapping = aes(x = Sce + barwidth + 0.03, y = Value,
                                                                                          fill = Ski), alpha = 0.7, width = 0.3, stat = "identity",
                                                                            position = "stack") + geom_errorbar(data = df1, mapping = aes(x = Sce,
                                                                                                                                          ymin = Value2 - ci, ymax = Value2 + ci), size = linesize *
                                                                                                                  0.6, width = errlinesize, color = "grey30") + geom_errorbar(data = df2,
                                                                                                                                                                              mapping = aes(x = Sce + barwidth + 0.03, ymin = Value2 -
                                                                                                                                                                                              ci, ymax = Value2 + ci), size = linesize * 0.6,
                                                                                                                                                                              width = errlinesize, color = "grey30") + scale_fill_npg() +
  theme_bw() + scale_x_continuous(breaks = c(4.15, 5.19,
                                             6.19), labels = c("LS", "TS", "CS")) + scale_y_continuous(name = "Economic output increment (tillion EURO)",
                                                                                                       sec.axis = sec_axis(trans = ~./55 * 1000, name = "Labor demand increment (million people)"),
                                                                                                       limits = c(0, max(df1$Value) * 1.125), expand = c(0,
                                                                                                                                                         0)) + theme(panel.grid.minor = element_blank(),
                                                                                                                                                                     panel.grid.major = element_blank(), panel.background = element_rect(fill = NA),
                                                                                                                                                                     axis.text.x = element_text(size = wordsize, colour = "black"),
                                                                                                                                                                     axis.text.y = element_text(size = wordsize, colour = "black"),
                                                                                                                                                                     axis.title.y = element_text(size = wordsize, colour = "black"),
                                                                                                                                                                     legend.text = element_text(size = wordsize, family = "A",
                                                                                                                                                                                                colour = "black")) + labs(x = " ") + theme(legend.position = "NULL",
                                                                                                                                                                                                                                           legend.title = element_blank())

setwd(rootpathway)
opfile <- paste("./Plot/", "P3.png", sep = "")
ggsave(opfile, P3, height = 6.45, width = 6.3, dpi = dpisize,
       units = "cm")
file.show(opfile)

# P4 The impacts of fiscal stimuli on labor demand and income inequality by country and skill level. --------------------------------------------------------------------

Plot2DataFunc <- function(df) {
  df$ci <- as.numeric(df$ci)
  df$Value <- as.numeric(df$Value)
  df$Value2 <- as.numeric(df$Value2)
  df$Ski <- factor(df$Ski, levels = c("L", "M", "H"),
                   ordered = T)
  df$Sce <- factor(df$Sce, levels = c("CS", "LS", "TS"),
                   ordered = T)
  df <- df[Sce != "ELS"]
  df <- ChangeCtyName(df)
  return(df)
}

dfa <- Plot2DataFunc(df_JB_Value_RR_sum_ski_sce_MMM_skip_EU_Rest)
Cty1 <- dfa[, .N, by = .(Cty)]

CtyGet <- unique(dfa$Cty)
for (j in CtyGet) {
  dfb <- dfa[Cty == j]
  P8_Cty <- ggplot() + geom_bar(data = dfb, aes(Sce,
                                                Value, fill = Ski), alpha = 0.7, width = 0.45,
                                stat = "identity", position = "stack") + geom_errorbar(data = dfb,
                                                                                       aes(x = Sce, ymin = Value2 - ci, ymax = Value2 +
                                                                                             ci), size = linesize * 0.6, width = errlinesize,
                                                                                       color = "grey30") + scale_fill_npg(alpha = 1) +
    scale_y_continuous(limits = c(0, (max(dfb$Value2) +
                                        max(dfb$ci)) * 1.05), expand = c(0, 0)) +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(), panel.background = element_rect(fill = NA),
                       plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text.y = element_text(angle = 0,
                                                                                           size = wordsize, colour = "black"), axis.text.x = element_text(angle = 0,
                                                                                                                                                          size = wordsize, colour = "black"), axis.title.y = element_text(size = wordsize,
                                                                                                                                                                                                                          color = "black"), axis.title.x = element_text(size = wordsize,
                                                                                                                                                                                                                                                                        color = "black"), legend.position = "NULL",
                       strip.text.x = element_text(size = wordsize, color = "black")) +
    facet_wrap(~Cty, scales = "free", ncol = 5, labeller = ) +
    xlab(" ") + ylab(" ")
  
  opfile <- paste("./Plot/", j, ".png", sep = "")
  ggsave(opfile, P8_Cty, height = 3.38, width = 3.15,
         dpi = dpisize, units = "cm")
}