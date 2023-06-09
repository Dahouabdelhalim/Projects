#04 Feb 2022
#v4.1.1
#-----------


# Cleaned up to consolidate all plot code for Artemia and exclude analysis code


### Paired with
# "01_Artemia_Neu_FA.csv"
# "02_Artemia_Phos_FA.csv"


#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(forcats)    # needed to rename factors
library(ggplot2)    # needed for plotting
library(ggpubr)     # needed to plot panels A and B together
library(cowplot)    # needed to plot panels A and B together





setwd(dirname(rstudioapi::getActiveDocumentContext()$path))   #sets WD to folder this R file is saved in



# --------------------------------------- Neutral ----------------------------------------------


artemia_neu_bad <- read.csv("01_Artemia_Neu_FA.csv", header = FALSE)

artemia_neu <- data.frame(t(artemia_neu_bad[-1]))
colnames(artemia_neu) <- artemia_neu_bad[, 1]


UniArt.neu <- artemia_neu %>%
  select(Treatment, C18.1n9, EPA, DHA)



UniArtPlot.neu <- UniArt.neu %>%
  gather('C18.1n9', 'EPA', 'DHA', key = "FA", value = "Percent") %>%
  group_by(Treatment, FA) %>%
  summarize(Avg = mean(as.numeric(Percent)), 
            SD = sd(as.numeric(Percent))) %>%
  mutate(min = Avg - SD) %>%
  mutate(max = Avg + SD)


UniArtPlot.neu$FA <- fct_recode(UniArtPlot.neu$FA, OA = "C18.1n9")






treatment.colors <- c("black", "blue", "orange", "gray")

x11()
uni.neu <- ggplot(UniArtPlot.neu,  aes(x = factor(FA), y = Avg, color = factor(Treatment))) +
  geom_point(size = 3, position = position_dodge(width = 0.40)) +
  geom_errorbar(aes(x = factor(FA), ymax = (max), ymin = (min), color = factor(Treatment)), width = 0.5, size = 1.0, position = position_dodge(width = 0.40)) +
  scale_color_manual(values = treatment.colors, labels = c("Control", "PUFA", "PUFA + E", "Unenriched")) +
  # ylim(0, 20) +
  labs(x = "", y = "Mean Percent Fatty Acid \\u00B1 SD") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
uni.neu














# --------------------------------------- Phospholipid ----------------------------------------------


artemia_phos_bad <- read.csv("02_Artemia_Phos_FA.csv", header = FALSE)

artemia_phos <- data.frame(t(artemia_phos_bad[-1]))
colnames(artemia_phos) <- artemia_phos_bad[, 1]


UniArt.phos <- artemia_phos %>%
  select(Treatment, C18.1n9, EPA, DHA)



UniArtPlot.phos <- UniArt.phos %>%
  gather('C18.1n9', 'EPA', 'DHA', key = "FA", value = "Percent") %>%
  group_by(Treatment, FA) %>%
  summarize(Avg = mean(as.numeric(Percent)), 
            SD = sd(as.numeric(Percent))) %>%
  mutate(min = Avg - SD) %>%
  mutate(max = Avg + SD)


UniArtPlot.phos$FA <- fct_recode(UniArtPlot.phos$FA, OA = "C18.1n9")






treatment.colors <- c("black", "blue", "orange", "gray")

x11()
uni.phos <- ggplot(UniArtPlot.phos,  aes(x = factor(FA), y = Avg, color = factor(Treatment))) +
  geom_point(size = 3, position = position_dodge(width = 0.40)) +
  geom_errorbar(aes(x = factor(FA), ymax = (max), ymin = (min), color = factor(Treatment)), width = 0.5, size = 1.0, position = position_dodge(width = 0.40)) +
  scale_color_manual(values = treatment.colors, labels = c("Control", "PUFA", "PUFA + E", "Unenriched")) +
  # ylim(0, 20) +
  labs(x = "", y = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
uni.phos





#Plotting the 2 together



x11()
pp <- ggarrange(uni.neu, uni.phos, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right") +
  draw_plot_label(label = c("(A) Neutral lipids", "(B) Phospholipids"), size = 14, x = c(0.05, 0.49), y = c(0.99, 0.99))

pp


ggsave("02_Art_uni_Neu_Phos.tif", plot = last_plot(), device = "tiff",
       width = 6.5, height = 4.0, units = "in", dpi = 150)

















