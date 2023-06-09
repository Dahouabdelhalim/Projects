#3.2: Generating base figure 3 of Castagner et al. (2022)
#Written by Ariane Castagner, finalized in June 2022

#Import packages

library(ggplot2)
library(reshape2)

# Set working directory (please change to your local directory for the project)

setwd("...//ExcessIceBetaRegression//")

# Import data

sentinel <- read.csv("data//SentinelBHdata.csv") #Import Sentinel training data

# Setting up

sentinelVIlogs <- sentinel[which(sentinel$VI.num.2 != "NA"), ] #Subsetting the Sentinel data to only include intervals which have a reported numeric visible ice content in both logs (log 2 did not use as many samples)
breaks <- seq(0, 100, by=5)
n <- 236

VI1 <- sentinelVIlogs$VI.num.1
VI1.cut <- cut(VI1, breaks, right=FALSE)
VI1.freq <- table(VI1.cut)
VI1.cumfreq <- c(0, cumsum(VI1.freq))
VI1.df <- as.data.frame(VI1.cumfreq)
VI1.df$breaks <- breaks

VI2 <- sentinelVIlogs$VI.num.2
VI2.cut <- cut(VI2, breaks, right=FALSE)
VI2.freq <- table(VI2.cut)
VI2.cumfreq <- c(0, cumsum(VI2.freq))
VI2.df <- as.data.frame(VI2.cumfreq)
VI2.df$breaks <- breaks

VI.merged.df <- merge(VI1.df, VI2.df, by="breaks")
VI.merged.melt <- reshape2::melt(VI.merged.df, id.var="breaks")
VI.merged.melt$valueperc <- VI.merged.melt$value / n * 100

# Plot

fig3 <- ggplot(VI.merged.melt, aes(x=breaks, y=valueperc, colour=variable)) +
  geom_line(size=1) +
  geom_point(size=2) +
  xlab("Visible ice content") +
  ylab("Cumulative frequency") +
  theme(text = element_text(size = 15)) +
  theme(legend.position = "none")

ggsave(filename = "output//figs//FinalFigure3.jpg", device="jpg", plot=fig3, scale = 1, dpi = 300)

fig3