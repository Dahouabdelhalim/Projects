#03 Apr 2022
#v4.1.1
#-----------

# Cleaned up for publication

### Paired with

#Tank means
# "Larval_end-point_tank_totals.csv"
# "Juvenile_end-point1_tank_totals.csv"

#Mass data
# "Larval_end-point_size.csv"
# "Juvenile_end-point1_size.csv"




#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(ggplot2)    # needed for plotting
library(ggpubr)     # needed to plot panels A and B together
library(cowplot)    # needed to plot panels A and B together



setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets WD to folder this R file is saved in




#--------------------------------------------------------------------------------------------------
### ----------------------------- Data for plotting for Phase 1 ----------------------------------
#--------------------------------------------------------------------------------------------------


#...................................................................................................
### #################################### Data Wrangling #########################################
#...................................................................................................



Larv.tanks <- read.csv("Larval_end-point_tank_totals.csv", header = TRUE)

Larv <- read.csv("Larval_end-point_size.csv", header = TRUE)


# Assigning each replicate tank to a treatment
## Treatments to match those in Larv.tanks
#     A = Control (OA), B = PUFA, C = PUFA + E
Larv.2 <- Larv %>%
  mutate(Treatment = Tank)

Larv.2$Treatment <- as.character(Larv.2$Treatment)
Larv.2$Treatment[grepl("1", Larv.2$Treatment)] <- "C"
Larv.2$Treatment[grepl("2", Larv.2$Treatment)] <- "C"
Larv.2$Treatment[grepl("3", Larv.2$Treatment)] <- "B"
Larv.2$Treatment[grepl("4", Larv.2$Treatment)] <- "A"
Larv.2$Treatment[grepl("5", Larv.2$Treatment)] <- "A"
Larv.2$Treatment[grepl("6", Larv.2$Treatment)] <- "B"
Larv.2$Treatment[grepl("7", Larv.2$Treatment)] <- "B"
Larv.2$Treatment[grepl("8", Larv.2$Treatment)] <- "A"
Larv.2$Treatment[grepl("9", Larv.2$Treatment)] <- "C"




#Calculating means and standard deviations for mass and length of inflated fish
Larv.size.inf.means <- Larv.2 %>%
  filter(SB == "1") %>%
  group_by(Tank) %>%
  summarize(Avg.larv.inf.mass = mean(Weight..mg.), Avg.larv.inf.TL = mean(length..mm.),
            SD.larv.inf.mass = sd(Weight..mg.), SD.larv.inf.TL = sd(length..mm.))

#Calculating means and standard deviations for mass and length of uninflated fish
Larv.size.uninf.means <- Larv.2 %>%
  filter(SB == "0") %>%
  group_by(Tank) %>%
  summarize(Avg.larv.uninf.mass = mean(Weight..mg.), Avg.larv.uninf.TL = mean(length..mm.),
            SD.larv.uninf.mass = sd(Weight..mg.), SD.larv.uninf.TL = sd(length..mm.))


#Calculating out beginning and end values for error bars that would be +/- 1 standard deviation
Larv.size.means.plus.dens <- full_join(Larv.tanks, Larv.size.inf.means, by = "Tank", "Treatment")
Larv.size.means.plus.dens2 <- full_join(Larv.size.means.plus.dens, Larv.size.uninf.means, by = "Tank") %>%
  mutate(mass.larv.inf.min = Avg.larv.inf.mass - SD.larv.inf.mass) %>%
  mutate(mass.larv.inf.max = Avg.larv.inf.mass + SD.larv.inf.mass) %>%
  mutate(mass.larv.uninf.min = Avg.larv.uninf.mass - SD.larv.uninf.mass) %>%
  mutate(mass.larv.uninf.max = Avg.larv.uninf.mass + SD.larv.uninf.mass) %>%
  mutate(length.larv.inf.min = Avg.larv.inf.TL - SD.larv.inf.TL) %>%
  mutate(length.larv.inf.max = Avg.larv.inf.TL + SD.larv.inf.TL) %>%
  mutate(length.larv.uninf.min = Avg.larv.uninf.TL - SD.larv.uninf.TL) %>%
  mutate(length.larv.uninf.max = Avg.larv.uninf.TL + SD.larv.uninf.TL) %>%
  mutate(Density.larv = Trial.Stocking/50) %>%
  mutate(SGR.P1 = 100 * (((Avg.larv.inf.mass/2.55)^(1/10)) - 1))




colnames(Larv.size.means.plus.dens2)[colnames(Larv.size.means.plus.dens2)=="Treatment"] <- "Treatment.larv"



#Calculating treatment averages of mass to determine how many times greater mass is from one treatment to another
avg.larv.inf.mass <- Larv.size.means.plus.dens2 %>%
  group_by(Treatment.larv) %>%
  summarize(Avg.Avg.larv.mass = mean(Avg.larv.inf.mass),
            SD.Avg.inf.mass = sd(Avg.larv.inf.mass),
            Avg.larv.SGR = mean(SGR.P1),
            SD.larv.SGR = sd(SGR.P1))

PUFA.over.control <- 20.0 / 17.3
PUFAE.over.control <- 19.8 / 17.3
PUFA.over.PUFAE <- 20.0 / 19.8

PUFA.over.control <- 22.7 / 21.0
PUFAE.over.control <- 22.6 / 21.0
PUFA.over.PUFAE <- 22.7 / 22.6





#...................................................................................................
### #################################### Figure 1a Mass Plot #########################################
#...................................................................................................


#Plot for part (A) - end of phase 1 - of mass plot

x11()
p1 <- ggplot(Larv.size.means.plus.dens2,  aes(x = Density.larv, y = Avg.larv.inf.mass, group = factor(Treatment.larv))) +
  geom_point(aes(color = factor(Treatment.larv), shape = factor(Treatment.larv)), size = 3) +
  geom_errorbar(aes(x = Density.larv, ymax = (mass.larv.inf.max), ymin = (mass.larv.inf.min), color = factor(Treatment.larv)), width = 0.5, size = 1) +
  geom_smooth(method = "lm", alpha = 0.0, size = 1, aes(color = factor(Treatment.larv))) +
  scale_color_grey(start = 0.75, end = 0.0, labels = c("Control", "PUFA", "PUFA + E")) +
  scale_shape_manual(values = c(17, 19, 18), labels = c("Control", "PUFA", "PUFA + E")) +
  scale_y_continuous(limits = c(10, 30)) +
  labs(x = " ", y = "Mass (mg) \\u00B1 1 SD") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 12, color = "black"),
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
p1 







#...................................................................................................
### #################################### Figure 1c SGR Plot #########################################
#...................................................................................................


#Plot for part (C) - end of phase 1 - of SGR plot

x11()
p3 <- ggplot(Larv.size.means.plus.dens2,  aes(x = Density.larv, y = SGR.P1, group = factor(Treatment.larv))) +
  geom_point(aes(color = factor(Treatment.larv), shape = factor(Treatment.larv)), size = 3) +
  geom_smooth(method = "lm", alpha = 0.0, size = 1, aes(color = factor(Treatment.larv))) +
  scale_color_grey(start = 0.75, end = 0.0, labels = c("Control", "PUFA", "PUFA + E")) +
  scale_shape_manual(values = c(17, 19, 18), labels = c("Control", "PUFA", "PUFA + E")) +
  scale_y_continuous(limits = c(18, 26)) +
  labs(x = "Density (no./L)", y = "SGR (%/d)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 12, color = "black"),
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
p3 








#--------------------------------------------------------------------------------------------------
### ----------------------------- Data for plotting for Phase 2 ----------------------------------
#--------------------------------------------------------------------------------------------------

#...................................................................................................
### #################################### Data Wrangling #########################################
#...................................................................................................


Juv.tanks <- read.csv("Juvenile_end-point1_tank_totals.csv", header = TRUE) %>%
  select(Tank, Treatment, Initial.density)


Juv.sizes <- read.csv("Juvenile_end-point1_size.csv", header = TRUE) %>%
  select(-Preserved_formalin, -Comments)
Juv.size <- Juv.sizes[complete.cases(Juv.sizes),]


# Assigning each replicate tank to a treatment
## Treatments to match those in Juv.tanks (which match those in Larv.tanks)
#     A = Control (OA), B = PUFA, C = PUFA + E
Juv.2 <- Juv.size %>%
  mutate(Treatment = Tank)

Juv.2$Treatment <- as.character(Juv.2$Treatment)
Juv.2$Treatment[grepl("1", Juv.2$Treatment)] <- "C"
Juv.2$Treatment[grepl("2", Juv.2$Treatment)] <- "C"
Juv.2$Treatment[grepl("3", Juv.2$Treatment)] <- "B"
Juv.2$Treatment[grepl("4", Juv.2$Treatment)] <- "A"
Juv.2$Treatment[grepl("5", Juv.2$Treatment)] <- "A"
Juv.2$Treatment[grepl("6", Juv.2$Treatment)] <- "B"
Juv.2$Treatment[grepl("7", Juv.2$Treatment)] <- "B"
Juv.2$Treatment[grepl("8", Juv.2$Treatment)] <- "A"
Juv.2$Treatment[grepl("9", Juv.2$Treatment)] <- "C"




#Calculating means and standard deviations for mass and length of inflated fish
Juv.size.inf.means <- Juv.2 %>%
  filter(Inflated. == "Inflated") %>%
  group_by(Tank) %>%
  summarize(Avg.inf.mass = mean(Wet.mass_g), Avg.inf.SL = mean(Length.standard_mm), Avg.inf.TL = mean(Length.total_mm),
            SD.inf.mass = sd(Wet.mass_g), SD.inf.SL = sd(Length.standard_mm), SD.inf.TL = sd(Length.total_mm))

#Calculating means and standard deviations for mass and length of uninflated fish
Juv.size.uninf.means <- Juv.2 %>%
  filter(Inflated. == "Uninflated") %>%
  group_by(Tank) %>%
  summarize(Avg.uninf.mass = mean(Wet.mass_g), Avg.uninf.SL = mean(Length.standard_mm), Avg.uninf.TL = mean(Length.total_mm),
            SD.uninf.mass = sd(Wet.mass_g), SD.uninf.SL = sd(Length.standard_mm), SD.uninf.TL = sd(Length.total_mm))


#Calculating out beginning and end values for error bars that would be +/- 1 standard deviation
Juv.size.means.plus.dens <- full_join(Juv.tanks, Juv.size.inf.means, by = "Tank")
Juv.size.means.plus.dens2 <- full_join(Juv.size.means.plus.dens, Juv.size.uninf.means, by = "Tank") %>%
  mutate(mass.inf.min = Avg.inf.mass - SD.inf.mass) %>%
  mutate(mass.inf.max = Avg.inf.mass + SD.inf.mass) %>%
  mutate(mass.uninf.min = Avg.uninf.mass - SD.uninf.mass) %>%
  mutate(mass.uninf.max = Avg.uninf.mass + SD.uninf.mass) %>%
  mutate(length.inf.min = Avg.inf.TL - SD.inf.TL) %>%
  mutate(length.inf.max = Avg.inf.TL + SD.inf.TL) %>%
  mutate(length.uninf.min = Avg.uninf.TL - SD.uninf.TL) %>%
  mutate(length.uninf.max = Avg.uninf.TL + SD.uninf.TL) %>%
  mutate(SGR.P2 = 100 * ((((Avg.inf.mass * 1000)/Larv.size.inf.means$Avg.larv.inf.mass)^(1/27)) - 1)) # See Eq. 3 & 4 in Crane et al. 2020




colnames(Juv.size.means.plus.dens2)[colnames(Juv.size.means.plus.dens2)=="Treatment"] <- "Treatment.juv"




#Calculating treatment averages of mass to determine how many times greater mass is from one treatment to another
avg.juv.mass <- Juv.size.means.plus.dens2 %>%
  group_by(Treatment.juv) %>%
  summarize(Avg.Avg.juv.mass = mean(Avg.inf.mass * 1000),
            SD.Avg.inf.mass = sd(Avg.inf.mass * 1000),
            Avg.juv.SGR = mean(SGR.P2),
            SD.juv.SGR = sd(SGR.P2))

PUFA.over.control <- 0.41 / 0.36
PUFAE.over.control <- 0.42 / 0.36

PUFA.over.control <- 407 / 357
PUFAE.over.control <- 423 / 357



PUFA.over.control <- 11.8 / 11.9
PUFAE.over.control <- 12.0 / 11.9





#...................................................................................................
### #################################### Figure 1b Plot #########################################
#...................................................................................................


#Plot for part (B) - end of phase 2 - of mass plot


x11()
p2 <- ggplot(Juv.size.means.plus.dens2,  aes(x = Initial.density, y = Avg.inf.mass, group = factor(Treatment.juv))) +
  geom_point(aes(color = factor(Treatment.juv), shape = factor(Treatment.juv)), size = 3) +
  geom_errorbar(aes(x = Initial.density, ymax = (mass.inf.max), ymin = (mass.inf.min), color = factor(Treatment.juv)), width = 0.5, size = 1) +
  geom_smooth(method = "lm", alpha = 0.0, size = 1, aes(color = factor(Treatment.juv))) +
  scale_color_grey(start = 0.75, end = 0.0, labels = c("Control", "PUFA", "PUFA + E")) +
  scale_shape_manual(values = c(17, 19, 18), labels = c("Control", "PUFA", "PUFA + E")) +
  scale_y_continuous(limits = c(0.15, 0.6)) +
  labs(x = " ", y = "Mass (g) \\u00B1 1 SD") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 12, color = "black"),
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
p2









#...................................................................................................
### #################################### Figure 1d Plot #########################################
#...................................................................................................


#Plot for part (D) - end of phase 2 - of SGR plot


x11()
p4 <- ggplot(Juv.size.means.plus.dens2,  aes(x = Initial.density, y = SGR.P2, group = factor(Treatment.juv))) +
  geom_point(aes(color = factor(Treatment.juv), shape = factor(Treatment.juv)), size = 3) +
  geom_smooth(method = "lm", alpha = 0.0, size = 1, aes(color = factor(Treatment.juv))) +
  scale_color_grey(start = 0.75, end = 0.0, labels = c("Control", "PUFA", "PUFA + E")) +
  scale_shape_manual(values = c(17, 19, 18), labels = c("Control", "PUFA", "PUFA + E")) +
  scale_y_continuous(limits = c(10, 14)) +
  labs(x = "Density (no./L)", y = "SGR (%/d)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 12, color = "black"),
        legend.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
p4








#...................................................................................................
### #################################### Combining the plots into Fig. 1 #########################################
#...................................................................................................



x11()
pppp <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "top") +
  draw_plot_label(label = c("(A)", "(B)", "(C)", "(D)"), size = 14, x = c(0.07, 0.58, 0.07, 0.57), y = c(0.95, 0.95, 0.48, 0.48))

pppp


ggsave("Fig1_Phase1&2_inf_mass_sgr_Manuscript.tif", plot = last_plot(), device = "tiff",
       width = 6.5, height = 5.5, units = "in", dpi = 600)



