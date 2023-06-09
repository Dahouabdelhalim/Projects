## Marta Shocket, Indiana University & UCLA
## Code for Can hot temperatures limit disease transmission? A test of mechanisms in a zooplankton-fungus system in Functional Ecology.
##
## Purpose: Plot motivating field pattern (possible upper thermal limit for epidemic starts)
##
## Table of Contents:
##      1. Set Working Directory, Load Libraries and Data
##      2. Functions for Bootstrapping Foraging Rate as a function of temperature
##      3. Console input for Bootstrapping Foraging Rate as a function of temperature

##### Set wd
setwd("~/Dropbox/Research Hall Lab/Temperature/Upper Thermal Constraints/Final Code and Data")

# Get data (same data set from Shocket et al. 2018 American Naturalist but with maximum Epilimnion temp added)
sd.data <- read.csv("UTC_FieldDataWeightedStartTemps.csv")

##### Figure

par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(4.5, 4.5, 1, 1))

hist(sd.data$EpiTempatStart, xlab = expression(paste("Epilimnion temperature (",degree,"C)")), breaks = seq(10, 35, 1),
     cex.lab = 1.25, main = "", col = rgb(0.2, 0.2, 0.2, 0.5), freq = TRUE, ylim = c(0, 20), xlim = c(10,35))
hist(sd.data$MaxEpiTemp, breaks = seq(10, 35, 1), col = rgb(0.9, 0.9, 0.9, 0.6), freq = TRUE, add = T, cex.axis = 1)
legend(x = 12, y = 0.2, legend = c("At epidemic start", "Season maximum"), pt.bg = c(rgb(0.2, 0.2, 0.2, 0.5), rgb(0.9, 0.9, 0.9, 0.5)), pch = 22, bty = "n", cex  = 1.1)

##### Colors!

# Color scheme 1
hist(sd.data$EpiTempatStart, xlab = expression(paste("Epilimnion temperature (",degree,"C)")), breaks = seq(10, 35, 1),
     cex.lab = 1.25, main = "", col = rgb(0, 0.8, 0, 0.5), freq = TRUE, ylim = c(0, 20), xlim = c(10,35))
hist(sd.data$MaxEpiTemp, breaks = seq(10, 35, 1), col = rgb(0, 0, 0.8, 0.5), freq = TRUE, add = T, cex.axis = 1)
legend(x = 12, y = 0.2, legend = c("At Epidemic Starts", "Season Max"), pt.bg = c(rgb(0, 0.8, 0, 0.5), rgb(0, 0, 0.8, 0.5)), pch = 22, bty = "n")

# Color scheme 2
hist(sd.data$EpiTempatStart, xlab = expression(paste("Epilimnion temperature (",degree,"C)")), breaks = seq(10, 35, 1),
     cex.lab = 1.25, main = "", col = rgb(0, 0.8, 0.8, 0.75), freq = TRUE, ylim = c(0, 20), xlim = c(10,35))
hist(sd.data$MaxEpiTemp, breaks = seq(10, 35, 1), col = rgb(0.8, 0, 0, 0.5), freq = TRUE, add = T, cex.axis = 1)
legend(x = 12, y = 0.2, legend = c("At Epidemic Starts", "Season Max"), pt.bg = c(rgb(0, 0.8, 0.8, 0.5), rgb(0.8, 0, 0, 0.5)), pch = 22, bty = "n")
