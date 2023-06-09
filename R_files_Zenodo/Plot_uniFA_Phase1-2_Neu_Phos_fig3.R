#04 Feb 2022
#v4.1.1
#-----------


# Cleaned up for publication


### Paired with
# "FA_10d_Neu_01.csv"
# "FA_10d_Phos_01.csv"
# "FA_37d_Neu_01.csv"
# "FA_37d_Phos_01.csv"




#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(ggplot2)    # needed for plotting
library(ggpubr)     # needed to plot panels A and B together
library(cowplot)    # needed to plot panels A and B together




setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #sets WD to folder this is saved in



#...................................................................................................
### #################################### Data Wrangling #########################################
#...................................................................................................



d10_bad.neu <- read.csv("FA_10d_Neu_01.csv", header = FALSE)

d10.neu <- data.frame(t(d10_bad.neu[-1]))
colnames(d10.neu) <- d10_bad.neu[, 1]




d10_bad.phos <- read.csv("FA_10d_Phos_01.csv", header = FALSE)

d10.phos <- data.frame(t(d10_bad.phos[-1]))
colnames(d10.phos) <- d10_bad.phos[, 1]





d37_bad.neu <- read.csv("FA_37d_Neu_01.csv", header = FALSE)

d37.neu <- data.frame(t(d37_bad.neu[-1]))
colnames(d37.neu) <- d37_bad.neu[, 1]




d37_bad.phos <- read.csv("FA_37d_Phos_01.csv", header = FALSE)

d37.phos <- data.frame(t(d37_bad.phos[-1]))
colnames(d37.phos) <- d37_bad.phos[, 1]









d10.neu.uni <- d10.neu %>%
  select(Treatment, Tank, C18.1n9, C20.5, C22.6) %>%
  mutate(Phase = 1) %>%
  mutate(FA = "Neutral")



d10.phos.uni <- d10.phos %>%
  select(Treatment, Tank, C18.1n9, C20.5, C22.6) %>%
  mutate(Phase = 1) %>%
  mutate(FA = "Phospholipid")



d37.neu.uni <- d37.neu %>%
  select(Treatment, Tank, C18.1n9, C20.5, C22.6) %>%
  mutate(Phase = 2) %>%
  mutate(FA = "Neutral")


d37.phos.uni <- d37.phos %>%
  select(Treatment, Tank, C18.1n9, C20.5, C22.6) %>%
  mutate(Phase = 2) %>%
  mutate(FA = "Phospholipid")






both.uni <- bind_rows(d10.neu.uni, d37.neu.uni)
three.uni <- bind_rows(both.uni, d10.phos.uni)
all.uni <- bind_rows(three.uni, d37.phos.uni)




FA.means <- all.uni %>%
  group_by(Phase, FA, Treatment) %>%
  summarize(Avg.OA = mean(as.numeric(C18.1n9)), Avg.EPA = mean(as.numeric(C20.5)), Avg.DHA = mean(as.numeric(C22.6)), 
            SD.OA = sd(as.numeric(C18.1n9)), SD.EPA = sd(as.numeric(C20.5)), SD.DHA = sd(as.numeric(C22.6))) %>%
  mutate(OA.min = Avg.OA - SD.OA) %>%
  mutate(OA.max = Avg.OA + SD.OA) %>%
  mutate(EPA.min = Avg.EPA - SD.EPA) %>%
  mutate(EPA.max = Avg.EPA + SD.EPA) %>%
  mutate(DHA.min = Avg.DHA - SD.DHA) %>%
  mutate(DHA.max = Avg.DHA + SD.DHA)



Neu.means <- FA.means %>%
  filter(FA == "Neutral")


Phos.means <- FA.means %>%
  filter(FA == "Phospholipid")








#...................................................................................................
### #################################### Plots of individual panels #########################################
#...................................................................................................




treatment.colors <- c("light grey", "grey55", "black")




x11()
pEPA.neu <- ggplot(Neu.means,  aes(x = factor(Phase), y = Avg.EPA, color = Treatment)) +
  geom_point(size = 3, position = position_dodge(width = 0.40)) +
  geom_errorbar(aes(x = factor(Phase), ymax = (EPA.max), ymin = (EPA.min)), width = 0.25, size = 1, position = position_dodge(width = 0.40)) +
  scale_color_manual(values = treatment.colors, labels = c("Control", "PUFA", "PUFA + E")) +
  ylim(0, 15) +
  labs(x = "Phase", y = "% Neutral Lipid EPA \\u00B1 SD") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 12, color = "black"),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
pEPA.neu



x11()
pDHA.neu <- ggplot(Neu.means,  aes(x = factor(Phase), y = Avg.DHA, color = Treatment)) +
  geom_point(size = 3, position = position_dodge(width = 0.40)) +
  geom_errorbar(aes(x = factor(Phase), ymax = (DHA.max), ymin = (DHA.min)), width = 0.25, size = 1, position = position_dodge(width = 0.40)) +
  scale_color_manual(values = treatment.colors, labels = c("Control", "PUFA", "PUFA + E")) +
  ylim(0, 15) +
  labs(x = "Phase", y = "% Neutral Lipid DHA \\u00B1 SD") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 12, color = "black"),
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
pDHA.neu




x11()
pEPA.phos <- ggplot(Phos.means,  aes(x = factor(Phase), y = Avg.EPA, color = Treatment)) +
  geom_point(size = 3, position = position_dodge(width = 0.40)) +
  geom_errorbar(aes(x = factor(Phase), ymax = (EPA.max), ymin = (EPA.min)), width = 0.25, size = 1, position = position_dodge(width = 0.40)) +
  scale_color_manual(values = treatment.colors, labels = c("Control", "PUFA", "PUFA + E")) +
  ylim(0, 15) +
  labs(x = "Phase", y = "% Phospholipid EPA \\u00B1 SD") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 12, color = "black"),
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
pEPA.phos




x11()
pDHA.phos <- ggplot(Phos.means,  aes(x = factor(Phase), y = Avg.DHA, color = Treatment)) +
  geom_point(size = 3, position = position_dodge(width = 0.40)) +
  geom_errorbar(aes(x = factor(Phase), ymax = (DHA.max), ymin = (DHA.min)), width = 0.25, size = 1, position = position_dodge(width = 0.40)) +
  scale_color_manual(values = treatment.colors, labels = c("Control", "PUFA", "PUFA + E")) +
  ylim(0, 30) +
  labs(x = "Phase", y = "% Phospholipid DHA \\u00B1 SD") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 12, color = "black"),
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
pDHA.phos














#...................................................................................................
### #################################### Plotting all panels together #########################################
#...................................................................................................



x11()
pp.ms <- ggarrange(pEPA.neu, pDHA.neu, pEPA.phos, pDHA.phos, ncol = 2, nrow = 2, common.legend = TRUE, legend = "top") +
  draw_plot_label(label = c("(A)", "(B)", "(C)", "(D)"), size = 14, x = c(0.10, 0.60, 0.10, 0.60), y = c(0.95, 0.95, 0.47, 0.47))

pp.ms




ggsave("Fig3_EPA_DHA_Neu_Phos_MS.tif", plot = last_plot(), device = "tiff",
       width = 6.5, height = 6.5, units = "in", dpi = 600)


















