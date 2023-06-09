library(dplyr)
library(ggplot2)

#set the working directory
setwd()

light <- read.csv("attenuation_light.csv")

#rename the population code
light$pop <- ifelse(light$pop=="glades","Everglades","Wakulla")

#make depth 38.1 as 38
light$depth <- ifelse(light$depth==38.1,38,light$depth)

#make another column with depth as a factor
light$depth2 <- as.factor(light$depth)

#calculate the average spectra for each depth at eac population
light_average <- light %>%
  group_by(depth2, pop, wavelength) %>%
  summarise(
    mean_photons= mean(photons)
  )



light_spectrum <- ggplot(light_average, aes(x=wavelength, y=mean_photons, color=depth2)) + 
  geom_line(size=0.75) + 
  scale_color_manual(values=c("gray80","gray50", "black"), name="Depth(cm)") +
  facet_grid(.~pop) + 
  xlab("wavelength (nm)") + 
  ylab("downwelling irradiance\\n(photons/cm2/s/nm)") + 
  #labs(fill = "Depth (cm)") +
  theme_bw() + 
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=13), 
        strip.text=element_text(size=14))

light_spectrum



