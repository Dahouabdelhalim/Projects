#3.5: Generating base figure 6 of Castagner et al. (2022)
#Written by Ariane Castagner, finalized in June 2022

#Import packages

library(ggplot2)
library(dplyr)

# Set working directory (please change to your local directory for the project)

setwd("...//ExcessIceBetaRegression//")

# Import data

datamassive <- read.csv("data//CryostratigraphyBHdataMI.csv")
pred_cum <- read.csv("output//data//pred_cum.csv") #Import predicted cumulative EI content
Profiles <- read.csv("data//FinalEIProfiles.csv") #Import final dataset derived from predicted and observed values

# Source functions

source("src//depthmeanmassive.R")
source("src//massivemeancumulative.R")
source("src//icysedimentmeancumulative.R")

# Prepare data

Profiles <- read.csv("data//FinalEIProfiles.csv") #Import final dataset derived from predicted and observed values
Profiles <- Profiles[-c(1, 9)]
ProfilesCumulative <- Profiles[c(1, 14, 15, 16, 17, 18, 19)]
ProfilesIntervals <- Profiles[c(1, 2, 3, 4, 5, 6, 7)]
ProfilesCumulativeMelted <- reshape2::melt(ProfilesCumulative, id.var = "Intervals")
ProfilesIntervalsMelted <- reshape2::melt(ProfilesIntervals, id.var = "Intervals")
pred_cum <- pred_cum[c(7, 8, 6, 4)]

# Define depth interval range

range <- seq(0.5, 10, by=0.5)

# Apply the depthmassivemean function to the massive ice data to derive mean massive ice content values per depth interval

datamassivemean <- depthmeanmassive(datamassive, range)
datamassivemean$Average.massive.ice <- datamassivemean$Average.massive.ice / 100
datamassivemean$conf.upr.massive <- datamassivemean$conf.upr.massive / 100
datamassivemean$conf.lwr.massive <- datamassivemean$conf.lwr.massive / 100

# Convert excess ice content to excess ice thickness (m) by multiplying the excess ice content by the interval thickness (0.5m)

datamassivemean$Average.massive.ice <- datamassivemean$Average.massive.ice * 0.5
datamassivemean$conf.upr.massive <- datamassivemean$conf.upr.massive * 0.5
datamassivemean$conf.lwr.massive <- datamassivemean$conf.lwr.massive * 0.5 

# Apply the massivemeancumulative function to derive the cumulative massive ice content for each interval

datamassivemeancumulative <- massivemeancumulative(datamassivemean, datamassivemean$Depth)

# Merge cumulative massive ice means with predicted icy seds for plotting

merged_data <- merge(pred_cum, datamassivemeancumulative, var.id="Depth")

# Melt dataset to use with ggplot

pred_boot_melted <- reshape2::melt(merged_data, id.var = "Depth")

# Plot

LINES <- c("boot.upr" = "blank", "estimate" = "solid", "boot.corr.lwr" = "blank", "Massive.ice" = "solid", "conf.upr.massive" = "blank", "conf.lwr.massive" = "blank")
COLOURS <- c("boot.upr" = "darkblue", "estimate" = "black", "boot.corr.lwr" = "darkblue", "Massive.ice" = "black", "conf.upr.massive" = "darkred", "conf.lwr.massive" = "darkred")
SIZES <- c("boot.upr" = 0.5, "estimate" = 1, "boot.corr.lwr" = 0.5, "Massive.ice" = 2, "conf.upr.massive" = 0.5, "conf.lwr.massive" = 0.5)

fig6 <- ggplot(pred_boot_melted, aes(x=Depth, y=value)) +
  geom_line(data=pred_boot_melted, aes(x=Depth, y=value, colour = variable, linetype=factor(variable), size=variable)) +
  coord_flip() +
  scale_x_reverse() +
  geom_ribbon(data=pred_cum, aes(y=estimate, ymin=boot.corr.lwr, ymax=boot.upr), alpha=0.3, fill="blue") +
  geom_ribbon(data=datamassivemeancumulative, aes(y=Massive.ice, ymin=conf.lwr.massive, ymax=conf.upr.massive), alpha=0.3, fill="red") +
  ylim(0, 3) +
  xlab("Depth (m)") +
  ylab("Excess ice thickness (m)") +
  scale_colour_manual(values=COLOURS) +
  scale_linetype_manual(values=LINES) +
  scale_size_manual(values=SIZES) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=12))

#Save plot

ggsave(filename = "output//figs//FinalFigure6.jpg", device="jpg", plot=fig6, scale = 1, dpi = 300)