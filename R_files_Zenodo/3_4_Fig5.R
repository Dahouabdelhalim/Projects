#3.4: Generating base figure 5 of Castagner et al. (2022)
#Written by Ariane Castagner, finalized in June 2022

#Import packages

library(ggplot2)

# Set working directory (please change to your local directory for the project)

setwd("...//ExcessIceBetaRegression//")

# Import data

datamassive <- read.csv("data//CryostratigraphyBHdataMI.csv") #Import cryostratigraphic dataset with massive ice intervals
pred_boot <- read.csv("output//data//pred_boot.csv") #Import predicted interval EI data

# Source functions

source("src//depthmeanmassive.R")
source("src//massivemeancumulative.R")
source("src//icysedimentmeancumulative.R")

# Prepare data

pred_boot <- pred_boot[c(7, 8, 6, 4)] #Subset and reorder the predicted EI data

# Define range and apply datamassivemean function to datamassive
range <- seq(0.5, 10, by=0.5)
datamassivemean <- depthmeanmassive(datamassive, range)
datamassivemean$Average.massive.ice = datamassivemean$Average.massive.ice / 100
datamassivemean$conf.upr.massive = datamassivemean$conf.upr.massive / 100
datamassivemean$conf.lwr.massive = datamassivemean$conf.lwr.massive / 100

#Merge massive ice means with predicted icy seds
merged_data <- merge(pred_boot, datamassivemean, var.id="Depth")

#Melt dataset to use with ggplot
pred_boot_melted <- reshape2::melt(merged_data, id.var = "Depth")

#Plot
LINES <- c("boot.upr" = "blank", "estimate" = "solid", "boot.corr.lwr" = "blank", "Average.massive.ice" = "solid", "conf.upr.massive" = "blank", "conf.lwr.massive" = "blank")
COLOURS <- c("boot.upr" = "darkblue", "estimate" = "black", "boot.corr.lwr" = "darkblue", "Average.massive.ice" = "black", "conf.upr.massive" = "darkred", "conf.lwr.massive" = "darkred")
SIZES <- c("boot.upr" = 0.5, "estimate" = 1, "boot.corr.lwr" = 0.5, "Average.massive.ice" = 2, "conf.upr.massive" = 0.5, "conf.lwr.massive" = 0.5)

fig5 <- ggplot(pred_boot_melted, aes(x=Depth, y=value)) +
  geom_line(data=pred_boot_melted, aes(x=Depth, y=value, colour = variable, linetype=factor(variable), size=variable)) +
  coord_flip() +
  scale_x_reverse() +
  geom_ribbon(data=pred_boot, aes(y=estimate, ymin=boot.corr.lwr, ymax=boot.upr), alpha=0.3, fill="blue") +
  geom_ribbon(data=datamassivemean, aes(y=Average.massive.ice, ymin=conf.lwr.massive, ymax=conf.upr.massive), alpha=0.3, fill="red") +
  ylim(0, 0.4) +
  xlab("Depth (m)") +
  ylab("Average excess ice content") +
  scale_colour_manual(values=COLOURS) +
  scale_linetype_manual(values=LINES) +
  scale_size_manual(values=SIZES) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=12))

#Save plot
ggsave(filename = "output//figs//FinalFigure5.tiff", device="tiff", plot=fig5, scale = 1, dpi = 300)