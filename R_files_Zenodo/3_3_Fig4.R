#3.3: Generating base figure 4 of Castagner et al. (2022)
#Written by Ariane Castagner, finalized in June 2022

#Import packages

library(ggplot2)
library(dplyr)

# Set working directory (please change to your local directory for the project)

setwd("...//ExcessIceBetaRegression//")

# Source functions

source("src//setxorder.R")
source("src//maxdepth.R")
source("src//rectplot.R")

# Import data

cryodata <- read.csv("data//Cryo566AllPredicted.csv") # This data file is different as it also includes all of the non-frozen intervals
landscapetypes <- read.csv("data//LandscapeTypes.csv")

#Prepare data

boreholes <- split(cryodata, cryodata$Borehole.ID) #List of data frames by borehole ID
cryodata <- lapply(boreholes, maxdepth) #Apply function to list of boreholes
cryodata <- bind_rows(cryodata) #Combine split borehole data
cryodata <- merge(cryodata, landscapetypes, by="Borehole.ID")

#Separate dataset by landscape types

landscapetypes <- split(cryodata, cryodata$Deposit.type.x) #Separate by landscape type

glaciolacustrine <- landscapetypes$Glaciolacustrine
glaciolacustrine <- glaciolacustrine[order(glaciolacustrine$maxdepth, glaciolacustrine$Borehole.ID),] #Sort by maxdepth, then borehole ID

moraine <- landscapetypes$Moraine
moraine <- moraine[order(moraine$maxdepth, moraine$Borehole.ID),] #Sort by maxdepth, then borehole ID

#Set up function for x axis plotting order based on maxdepth. Make sure the dataset is ordered by maxdepth, then borehole ID. Combine list of data sets, then split again once xorder is done to use in the plotting function with lapply (see below).
#Dataset is already ordered by maxdepth and borehole ID

cryodata <- cryodata[with(cryodata, order(cryodata$maxdepth, cryodata$Borehole.ID)), ]

#Apply xorder to each dataset
cryodata <- setxorder(cryodata)
glaciolacustrine <- setxorder(glaciolacustrine)
moraine <- setxorder(moraine)

#Fix the borehole IDs if needed
cryodata$Borehole.ID <- as.character(cryodata$Borehole.ID)
cryodata$Borehole.ID <- as.factor(cryodata$Borehole.ID)
glaciolacustrine$Borehole.ID <- as.character(glaciolacustrine$Borehole.ID)
glaciolacustrine$Borehole.ID <- as.factor(glaciolacustrine$Borehole.ID)
moraine$Borehole.ID <- as.character(moraine$Borehole.ID)
moraine$Borehole.ID <- as.factor(moraine$Borehole.ID)

#Split datasets by borehole ID to use in plotting function below. If the plotting function is returning blank boreholes in the plots, you may need to convert the borehole IDs to strings, then back to factors.
boreholes <- split(cryodata, cryodata$Borehole.ID)
glaciolacustrine <- split(glaciolacustrine, glaciolacustrine$Borehole.ID) #Split by borehole ID
moraine <- split(moraine, moraine$Borehole.ID) #Split by borehole ID


# PLOTTING ------------------------------------------------------------------------------

#Plot for all samples

jpeg(filename = "output\\\\figs\\\\fig4.jpg", quality = 100, bg = "white", res = 100, width = 973, height = 659)
plot(rev(c(0, 360)), -(c(0, 22)), type = "n", xlab = "Boreholes", ylab = "Depth (m)",
     main = "Stratigraphic logs of predicted excess ice content in glaciolacustrine and moraine boreholes \\nalong the ITH by campaign (n=351 boreholes)")
rectplotall <- lapply(boreholes, rectplot)
dev.off()

#Plot for only glaciolacustrine samples

jpeg(filename = "output\\\\figs\\\\RectPlotGL.jpg", quality = 100, bg = "white", res = 100, width = 973, height = 659)
plot(rev(c(0, 250)), -(c(0, 22)), type = "n", xlab = "Boreholes", ylab = "Depth (m)",
     main = "Stratigraphic logs of predicted excess ice content in glaciolacustrine boreholes \\nalong the ITH by campaign (n=241 boreholes)")
rectplotgl <- lapply(glaciolacustrine, rectplot)
dev.off()

#Plot for only moraine samples

jpeg(filename = "output\\\\figs\\\\RectPlotMN.jpg", quality = 100, bg = "white", res = 100, width = 973, height = 659)
plot(rev(c(0, 120)), -(c(0, 22)), type = "n", xlab = "Boreholes", ylab = "Depth (m)",
     main = "Stratigraphic logs of predicted excess ice content in moraine boreholes \\nalong the ITH by campaign (n=109 boreholes)")
rectplotmn <- lapply(moraine, rectplot)
dev.off()
