# This script was made and run in 
# R version 4.0.2 (2020-06-22)
# Authors: Emily Cavaliere and Hilary Dugan 

# Load packages
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)
library(ggpubr) # for theme_classic2()

# Read in data
mend_chem <- read_csv("Mendota_UnderIceProfiles_2018_2019.csv") %>% 
  dplyr::mutate(winterconditions = if_else(year(DateTime) == 2018, 'Clear Ice', 'Snow and Ice'))
mend_light <- read_csv("MendotaUnderIceLight_2018_2019.csv") %>% 
  mutate(Light_lum_m2 = Light_lum_ft2/0.092903) %>% # convert from lum/ft2 to lum/m2
  dplyr::mutate(winterconditions = if_else(year(DateTime) == 2018, 'Clear Ice', 'Snow and Ice'))

# Mendota light intensity plot 
p1 = ggplot(mend_light %>% filter(month(DateTime) == 2)) + 
  geom_density(aes(x = Light_lum_m2, fill = winterconditions), 
               color = 'black', alpha = 0.7, position="identity") +
  scale_fill_manual(name = 'Winter Conditions',
                    values = c('#82d2ee','grey60')) +
  scale_x_log10() +
  xlab(expression("Lumens per"~ m^{2})) +
  ylab("Density Probability") +
  theme(legend.position = "bottom",
        legend.box = "vertical") +
  theme_classic2(base_size = 8); p1

# Mendota under-ice profile plots
temp <- ggplot(mend_chem) +
  geom_path(aes(y=Depthm, x=TempC, color=winterconditions), size=1.5) +
  labs(x = "Temperature (Â°C)", y = "Depth (m)", color = "Winter Year") +
  scale_y_continuous(trans = "reverse") +
  scale_color_manual(values = c('#82d2ee','grey60'), guide=F) +
  theme_classic2(base_size = 8); temp

xtitle <- expression(atop(paste("Specific Conductance", " (",mu, "S ", cm^-1,")")))
cond<-ggplot(mend_chem) +
  geom_path(aes(y=Depthm, x=SpConduScm, color=winterconditions), size=1.5) +
  labs( y = "Depth (m)", color = "Winter Year") +
  xlab(xtitle) +
  scale_y_continuous(trans = "reverse") +
  scale_color_manual(values = c('#82d2ee','grey60'), guide=F) +
  theme_classic2(base_size = 8); cond

ph <- ggplot(mend_chem) +
  geom_path(aes(y=Depthm, x=pH, color=winterconditions), size=1.5) +
  labs( y = "Depth (m)", color = "Winter Year") +
  xlab("pH") +
  scale_y_continuous(trans = "reverse") +
  scale_color_manual(values = c('#82d2ee','grey60'), guide=F) +
  theme_classic2(base_size = 8); ph

chl <- ggplot(mend_chem) +
  geom_path(aes(y=Depthm, x=ChlorophyllugL, color=winterconditions), size=1.5) +
  labs(y = "Depth (m)", color = "Winter Year") +
  xlab(bquote(paste('Chlorophyll ',italic("a"), ' (' *mu,'g ', L^-1*')'))) +
  scale_y_continuous(trans = "reverse") +
  scale_color_manual(values = c('#82d2ee','grey60'), guide=F) +
  theme_classic2(base_size = 8);chl

dosat <- ggplot(mend_chem) +
  geom_path(aes(y=Depthm, x=ODOsat, color=winterconditions), size=1.5) +
  labs(y = "Depth (m)", color = "Winter Year") +
  xlab("Dissolved Oxygen Saturation (%)") +
  scale_y_continuous(trans = "reverse") +
  scale_color_manual(values = c('#82d2ee','grey60'), guide=F) +
  theme_classic2(base_size = 8); dosat

# Patchwork plots together and save figure
combined <- p1 + temp + chl + dosat + ph + cond + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') & 
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        plot.tag = element_text(size = 8)); combined

# Output 
# ggsave('ME_LICC_casestudy.png', dpi = 500, unit="in", width = 6.5, height=5)

