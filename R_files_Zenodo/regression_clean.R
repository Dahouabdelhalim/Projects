### Linear regression plot for the Ammophila-Paraxenos project (A Double-Edged Sword: Parental care increases risk of offspring infection by a maternally-vectored parasite)
## Authors: Rebecca Jean Millena, Jay Rosenheim

library(tidyverse)
library(ggpubr)
library(ggplot2)

## Loading in a data subset containing only the provisioning data and parasitism rates for each species
preygraph <- read_csv("meanprey_percentparasitized.csv")

## Loading in the species as factors in order to color the points in a gradient
preygraph$Species <- factor(preygraph$species, levels = c("aberti", "pruinosa", "azteca", "harti", "urnaria", "placida", "femurrubra",  "dysmica", "juncea", "boharti", "procera", "wrightii", "marshi", "stangei", "nigricans", "zanthoptera"))
preygraph$Species

## Plotting the regression with ggplot2
ggplot(data = preygraph, mapping = aes(x = mean_prey, y = percent_parasitized)) +
  annotate("text", x = 1, y = 8, label = "y = -0.98 + 0.98x", hjust = 0, family = "Times") +
  ## method = 'lm' means that I am fitting a linear regression model; shaded area is the 95% confidence interval that signifies there is 95% confidence that the true regression line lies within the shaded region
  geom_smooth(method ='lm',color = "black") +
  xlab("Mean Prey/Nest") +
  ylab("% Parasitized") +
  theme_gray(base_family = "Times") +
  stat_cor(color = "black", family = "Times") +
  geom_jitter(alpha = 0.8, aes(color = Species)) +
  theme(text = element_text(color = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(plot.caption = element_text(hjust = 0.3)) +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  guides(shape = guide_legend(override.aes = list(size = 1))) +
  theme(legend.key.size = unit(0.001, "cm")) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(face = "italic", color = "#595959"),
        axis.ticks = element_line(color = "#595959"),
        axis.line = element_line(color = "#595959")) +
  # labs(color = "Ammophila Species") +
  ## Adding the numbers of specimens scored for each species to the legend
  scale_color_hue(labels = c("aberti  (622)", "pruinosa  (1412)", "azteca  (2465)", "harti  (32)", "urnaria  (332)", "placida  (445)", "femurrubra  (435)",  "dysmica  (206)", "juncea  (125)", "boharti  (307)", "procera  (1012)", "wrightii  (659)", "marshi  (346)", "stangei  (241)", "nigricans  (179)", "zanthoptera  (138)"))

## Saving the plot as a .tiff with ggsave
ggsave("../ammophila_strepsiptera/strepgraph_paper_2022.pdf")

## Plotting a subset tree for comparison in tandem with the phylogenetic contrast
preygraph_subset <- read_csv("meanprey_percentparasitized_subset.csv")
preygraph_subset$Species <- factor(preygraph_subset$species, levels = c("aberti", "azteca", "urnaria", "femurrubra",  "dysmica", "procera", "wrightii", "marshi", "stangei", "nigricans"))

ggplot(data = preygraph_subset, mapping = aes(x = mean_prey, y = percent_parasitized)) +
  annotate("text", x = 1, y = 8, label = "y = -1.2 + 1.1x", hjust = 0, family = "Lato") +
  ## method = 'lm' means that I am fitting a linear regression model; shaded area is the 95% confidence interval that signifies there is 95% confidence that the true regression line lies within the shaded region
  geom_smooth(method ='lm',color = "black") +
  xlab("Mean Prey/Nest") +
  ylab("% Parasitized") +
  theme_gray(base_family = "Lato") +
  stat_cor(color = "black", family = "Lato") +
  geom_jitter(alpha = 0.8, aes(color = Species)) +
  theme(text = element_text(color = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(plot.caption = element_text(hjust = 0.3)) +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  guides(shape = guide_legend(override.aes = list(size = 1))) +
  theme(legend.key.size = unit(0.001, "cm")) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(face = "italic", color = "#595959"),
        axis.ticks = element_line(color = "#595959"),
        axis.line = element_line(color = "#595959")) +
  ## Adding the numbers of specimens scored for each species to the legend
  scale_color_hue(labels = c("aberti  (622)", "azteca  (2465)", "urnaria  (332)", "femurrubra  (435)",  "dysmica  (206)", "procera  (1012)", "wrightii  (659)", "marshi  (346)", "stangei  (241)", "nigricans  (179)"))

ggsave("../ammophila_strepsiptera/strepgraph_subset_pres2022.tiff")
