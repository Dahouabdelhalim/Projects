rm(list=ls())

########
# Xing et al., Extended Data Figure 06
########


##goal: Visualize skull size ranges through box plots

###SKULL SIZE

#please set working directory
setwd(path/to/your/working/directory)

#call libraries
library(ggplot2)

#load data
skull <- read.csv("skull.orbit.avg.csv", header = TRUE)
weenie <- data.frame(taxon = "HPG-15-3", ID = "HPG-15-3", skull_length = 7.1, orbit_length = 4.92, group_1 = "HPG-15-3", group_2 = "HPG-15-3")
dat <- rbind(skull, weenie)

#plot
bxplot <- ggplot(dat, aes(x = reorder(group_2, log10(skull_length),FUN=median), y = log10(skull_length), fill=group_2)) +
  geom_boxplot(alpha = 0.5) +
  theme_classic(base_size = 12)+ theme(legend.position="none") +
  labs(x = "Birds", y = "Log10 Postnasal Skull length, mm")
bxplot

#save as PDF
pdf("Ext_Data_Fig_06.pdf", 6, 4)
bxplot
dev.off()