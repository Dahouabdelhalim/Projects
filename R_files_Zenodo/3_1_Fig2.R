#3.1: Generating base figure 2 of Castagner et al. (2022)
#Written by Ariane Castagner, finalized in June 2022

#Import packages

library(ggplot2)

# Set working directory (please change to your local directory for the project)

setwd("...//ExcessIceBetaRegression//")

# Import data

sentinel <- read.csv("data//SentinelBHdata.csv") #Import Sentinel training data

# Order visible ice classes as factors

sentinel$Visible.ice.1 <- factor(sentinel$Visible.ice.1, levels=c("No visible ice", "Low", "Medium to high", "High"))

# Plot

fig2 <- ggplot(sentinel, aes(x=VI.num.1, y=EI)) +
  geom_point() +
  xlab("Logged visible ice estimate (%)") +
  ylab("Measured volumetric excess ice content") +
  scale_x_continuous(breaks=c(seq(0,100,10))) +
  scale_y_continuous(breaks=c(seq(0,1,0.1))) +
  theme(text = element_text(size = 15))

ggsave(filename = "output//figs//FinalFigure2.jpg", device="jpg", plot=fig2, scale = 1, dpi = 300)

fig2