# load packages
library(readr)
library(dplyr)


# Quadrant level data (5.8x5.8km) ####

# load data
quadrant_data <- read_csv("Quadrant_data.csv")
str(quadrant_data)

quadrant_data$QUADRANT <- as.factor(quadrant_data$QUADRANT)
quadrant_data$LANDUSE <- as.factor(quadrant_data$LANDUSE)

# change proportion land use (near-natural, agriculture, urban) to percentage from 0 to 100%
quadrant_data$ALL_NATURAL100 <- 100 * quadrant_data$ALL_NATURAL
quadrant_data$ALL_AGRI100 <- 100 * quadrant_data$ALL_AGRI
quadrant_data$URBAN100 <- 100 * quadrant_data$URBAN

# get average values and standard deviation across quadrants####

# temperature means and standard deviation
summary(quadrant_data$TEMPERATURE)
sd(quadrant_data$TEMPERATURE)

# precipitation
summary(quadrant_data$PRECIPITATION)
sd(quadrant_data$PRECIPITATION)

# natural-dominated
summary(quadrant_data$ALL_NATURAL100)
sd(quadrant_data$ALL_NATURAL100)

# agriculture
summary(quadrant_data$ALL_AGRI100)
sd(quadrant_data$ALL_AGRI100)

# urban
summary(quadrant_data$URBAN100) 
sd(quadrant_data$URBAN100)

par(mfrow=c(2,2))

# get correlation coefficients for relationship between temperature and land use on quadrant level
cor.test(quadrant_data$TEMPERATURE, quadrant_data$ALL_NATURAL100)
cor.test(quadrant_data$TEMPERATURE, quadrant_data$ALL_AGRI100)
cor.test(quadrant_data$TEMPERATURE, quadrant_data$URBAN100)

# plot correlations of temperature with land use and add Pearson's r - MAIN TEXT Figure 5

par(mfrow = c(2, 2), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(3, 1, 0))

plot(quadrant_data$TEMPERATURE ~ quadrant_data$ALL_NATURAL100, 
     outer=T,
     cex.lab = 2,
     cex.axis=1.5,
     ylab="Multi-annual mean temperature (°C)", xlab="",
     xlim = c(0, 100),
     xaxt="n",
     pch = 19,
     cex = 1.5,
     col="grey30",
     las=1)
text(c(98), c(9.8), c("(a)"), cex=2)
text(c(50), c(4.1), c("r = -0.24"), cex=2)

plot(quadrant_data$TEMPERATURE ~ quadrant_data$ALL_AGRI100,
     ylab="", xlab="", 
     xlim = c(0, 100),
     cex.axis=1.5,
     pch = 19,
     cex = 1.5,
     col="grey30",
     yaxt="n")
text(c(98), c(9.8), c("(b)"), cex=2)
text(c(50), c(4.1), c("r = -0.03"), cex=2)

plot(quadrant_data$TEMPERATURE ~ quadrant_data$URBAN100, 
     ylab="", xlab="", 
     xlim = c(0, 100),
     cex.axis=1.5,
     pch = 19,
     cex = 1.5,
     col="grey30",
     las=1)
text(c(98), c(9.8), c("(c)"), cex=2)
text(c(50), c(4.1), c("r = 0.33"), cex=2)

plot(quadrant_data$TEMPERATURE ~ quadrant_data$URBAN100, 
     axes = FALSE, 
     #bty= "l",
     type = "n",
     xaxt="n",
     yaxt="n"
     )


# get correlation coefficients for relationship between precipitation and land use on quadrant level
cor.test(quadrant_data$PRECIPITATION, quadrant_data$ALL_NATURAL100)
cor.test(quadrant_data$PRECIPITATION, quadrant_data$ALL_AGRI100)
cor.test(quadrant_data$PRECIPITATION, quadrant_data$URBAN100)

# get correlation coefficients for relationship between precipitation and temperature
cor.test(quadrant_data$PRECIPITATION, quadrant_data$TEMPERATURE)

# plot correlations of land use with precipitation and add Pearson's r - APPENDIX Figure S1

par(mfrow = c(2, 2), mar = c(0, 0, 0, 0), oma = c(4, 6, 0.5, 0.5), tcl = -0.25, mgp = c(4, 1, 0))

plot(quadrant_data$PRECIPITATION ~ quadrant_data$ALL_NATURAL100, 
     outer=T,
     cex.lab = 2,
     cex.axis=1.5,
     ylab="Annual precipitation (mm)", xlab="", 
     ylim = c(0, 3000),
     xlim = c(0, 100),
     xaxt="n",
     pch = 19,
     cex = 1.5,
     col="grey30",
     las=1)
text(c(98), c(2900), c("A"), cex=2)
text(c(50), c(50), c("r = 0.3"), cex=2)

plot(quadrant_data$PRECIPITATION ~ quadrant_data$ALL_AGRI100,
     ylab="", xlab="", 
     ylim = c(0, 3000),
     xlim = c(0, 100),
     cex.axis=1.5,
     pch = 19,
     cex = 1.5,
     col="grey30",
     yaxt="n")
text(c(98), c(2900), c("B"), cex=2)
text(c(50), c(50), c("r = -0.23"), cex=2)

plot(quadrant_data$PRECIPITATION ~ quadrant_data$URBAN100, 
     ylab="", xlab="", 
     ylim = c(0, 3000),
     xlim = c(0, 100),
     cex.axis=1.5,
     pch = 19,
     cex = 1.5,
     col="grey30",
     las=1)
text(c(98), c(2900), c("C"), cex=2)
text(c(50), c(50), c("r = -0.12"), cex=2)

plot(quadrant_data$PRECIPITATION ~ quadrant_data$URBAN100, # empty panel
     axes = FALSE,
     type = "n",
     xaxt="n",
     yaxt="n"
)



# Plot level data (179 plots, 1-km landscape radius) ####

# load plot data
plot_data <- read_csv("plot_data.csv") 

plot_data$PlotID <- as.factor(plot_data$PlotID)
plot_data$habitat <- as.factor(plot_data$habitat)

str(plot_data)

# change proportion land use (forest, grassland, arable, settlement) to percentage from 0 to 100%
plot_data$forest100 <- 100 * plot_data$forest
plot_data$agricultur100 <- 100 * plot_data$agriculture
plot_data$urban100 <- 100 * plot_data$urban
plot_data$grassland100 <- 100 * plot_data$grassland

# get number of plots per local land-use type
plot_data %>% 
  group_by(habitat) %>%
  summarise(no_rows = length(habitat))

# get average values and standard deviation across 179 plots ####

# temperature
summary(plot_data$temp)
sd(plot_data$temp)

# precipitation
summary(plot_data$precip)
sd(plot_data$precip)

# forest 
summary(plot_data$forest100)
sd(plot_data$forest100)

# agriculture
summary(plot_data$agricultur100)
sd(plot_data$agricultur100)

# urban
summary(plot_data$urban100) 
sd(plot_data$urban100)

# grassland
summary(plot_data$grassland100)  
sd(plot_data$grassland100)

# edge density
summary(plot_data$edge_density) 
sd(plot_data$edge_density)

# elevation
summary(plot_data$elevation)
sd(plot_data$elevation)

# natural_se
summary(100*plot_data$natural_se)
sd(100*plot_data$natural_se)

# water
summary(100*plot_data$water) 
sd(100*plot_data$water)

# get correlation coefficients for relationship between temperature and land use/edge density on plot level

cor.test(plot_data$temp, plot_data$forest100)
cor.test(plot_data$temp, plot_data$agricultur100)
cor.test(plot_data$temp, plot_data$urban100)
cor.test(plot_data$temp, plot_data$grassland100)
cor.test(plot_data$temp, plot_data$edge_density)

# get correlation coefficients for relationship between precipitation and land use/edge density on plot level

cor.test(plot_data$precip, plot_data$forest100)
cor.test(plot_data$precip, plot_data$agricultur100)
cor.test(plot_data$precip, plot_data$urban100)
cor.test(plot_data$precip, plot_data$grassland100)
cor.test(plot_data$precip, plot_data$edge_density)

# get correlation coefficients for relationship between precipitation and temperature on plot level

cor.test(plot_data$precip, plot_data$temp)

# plot correlations of land use with temperature and add Pearson's r - MAIN TEXT Figure 5
par(mfrow = c(2, 2), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(3, 1, 0))

plot(plot_data$temp ~ plot_data$forest100, 
     outer=T,
     cex.lab = 2,
     cex.axis=1.5,
     ylab="Multi-annual mean temperature (°C)", xlab="", 
     ylim = c(4, 10),
     xlim = c(0, 100),
     xaxt="n",
     pch = 19,
     cex = 1.5,
     col="grey30",
     las=1)
text(c(98), c(9.9), c("(d)"), cex=2)
text(c(50), c(4.1), c("r = -0.29"), cex=2)

plot(plot_data$temp ~ plot_data$grassland100,
     ylab="", xlab="", 
     ylim = c(4, 10),
     xlim = c(0, 100),
     cex.axis=1.5,
     pch = 19,
     cex = 1.5,
     col="grey30",
     yaxt="n",
     xaxt="n")
text(c(98), c(9.9), c("(e)"), cex=2)
text(c(50), c(4.1), c("r = -0.17"), cex=2)

plot(plot_data$temp ~ plot_data$agricultur100, 
     ylab="", xlab="", 
     ylim = c(4, 10),
     xlim = c(0, 100),
     cex.axis=1.5,
     pch = 19,
     cex = 1.5,
     col="grey30",
     las=1)
text(c(98), c(9.9), c("(f)"), cex=2)
text(c(50), c(4.1), c("r = 0.21"), cex=2)

plot(plot_data$temp ~ plot_data$urban100, 
     ylab="", xlab="", 
     ylim = c(4, 10),
     xlim = c(0, 100),
     cex.axis=1.5,
     pch = 19,
     cex = 1.5,
     col="grey30",
     las=1,
     yaxt="n"
)
text(c(98), c(9.9), c("(g)"), cex=2)
text(c(50), c(4.1), c("r = 0.26"), cex=2)

# plot correlations of land use with precipitation and add Pearson's r - APPENDIX Figure S1

par(mfrow = c(2, 2), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(3, 1, 0))

plot(plot_data$precip ~ plot_data$forest100,
     outer=T,
     cex.lab = 2,
     cex.axis=1.5,
     ylab="Annual precipitation (mm)", xlab="", 
     ylim = c(0,2500),
     xlim = c(0, 100),
     xaxt="n",
     pch = 19,
     cex = 1.5,
     col="grey30",
     las=1)
text(c(98), c(2400), c("D"), cex=2)
text(c(50), c(50), c("r = 0.27"), cex=2)

plot(plot_data$precip ~ plot_data$grassland100,
     ylab="", xlab="", 
     ylim = c(0,2500),
     xlim = c(0, 100),
     cex.axis=1.5,
     pch = 19,
     cex = 1.5,
     col="grey30",
     yaxt="n",
     xaxt="n")
text(c(98), c(2400), c("E"), cex=2)
text(c(50), c(50), c("r = 0.31"), cex=2)

plot(plot_data$precip ~ plot_data$agricultur100, 
     ylab="", xlab="", 
     ylim = c(0,2500),
     xlim = c(0, 100),
     cex.axis=1.5,
     pch = 19,
     cex = 1.5,
     col="grey30",
     las=1)
text(c(98), c(2400), c("F"), cex=2)
text(c(50), c(50), c("r = -0.46"), cex=2)

plot(plot_data$precip ~ plot_data$urban100, 
     ylab="", xlab="", 
     ylim = c(0,2500),
     xlim = c(0, 100),
     cex.axis=1.5,
     pch = 19,
     cex = 1.5,
     col="grey30",
     las=1,
     yaxt="n"
)
text(c(98), c(2400), c("G"), cex=2)
text(c(50), c(50), c("r = -0.06"), cex=2)


# Relationships between landscape composition and configuration ####

# correlations of edge density (landscape configuration) with % land use types (landscape composition)

# all plots (n=179)
cor.test(plot_data$forest, plot_data$edge_density)
cor.test(plot_data$grassland, plot_data$edge_density)
cor.test(plot_data$agriculture, plot_data$edge_density)
cor.test(plot_data$urban, plot_data$edge_density)


# habitat-specific subsets
forest <- plot_data[plot_data$habitat == "forest",]
cor.test(forest$forest, forest$edge_density)

grassland <- plot_data[plot_data$habitat == "grassland",]
cor.test(grassland$grassland, grassland$edge_density)

agriculture <- plot_data[plot_data$habitat == "agriculture",]
cor.test(agriculture$agriculture, agriculture$edge_density)

urban <- plot_data[plot_data$habitat == "urban",]
cor.test(urban$urban, urban$edge_density)

# plot relationships between edge density and land use types (%) for all 179 compared to habitat-specific subplots  - APPENDIX Figure S2
par(mfrow=c(4,2), mar=c(4,5,2,1))

plot(plot_data$edge_density ~ plot_data$forest100, ylab="Edge density (m/ha)", xlab="Forest cover (%)", main = "All plots", xlim=c(0,100), ylim=c(0,70))
plot(forest$edge_density ~ forest$forest100, ylab="Edge density (m/ha)", xlab="Forest cover (%)", main= "Habitat subsets", xlim=c(0,100), ylim=c(0,70))
plot(plot_data$edge_density ~ plot_data$grassland100, ylab="Edge density (m/ha)", xlab="Grassland cover (%)", xlim=c(0,100), ylim=c(0,70))
plot(grassland$edge_density ~ grassland$grassland100, ylab="Edge density (m/ha)", xlab="Grassland cover (%)", xlim=c(0,100), ylim=c(0,70))
plot(plot_data$edge_density ~ plot_data$agricultur100, ylab="Edge density (m/ha)", xlab="Arable  land cover (%)", xlim=c(0,100), ylim=c(0,70))
plot(agriculture$edge_density ~ agriculture$agricultur100, ylab="Edge density (m/ha)", xlab="Arable land cover (%)", xlim=c(0,100), ylim=c(0,70))
plot(plot_data$edge_density ~ plot_data$urban100, ylab="Edge density (m/ha)", xlab="Settlement cover (%)", xlim=c(0,100), ylim=c(0,70))
plot(urban$edge_density ~ urban$urban100, ylab="Edge density (m/ha)", xlab="Settlement cover (%)", xlim=c(0,100), ylim=c(0,70))


# plot histograms of correlations from 10,000 random selections

# Histogram correlations####
settlement_random <- read_csv("settlement_random.csv")
grassland_random <- read_csv("grassland_random.csv")
forest_random <- read_csv("forest_random.csv")
arable_random <- read_csv("arable_random.csv")


# plot correlations from 10,000 random relections for all land-use types and add lines of actual correlations acorss all 179 plots and across land-use-specific plots - MAIN TEXT Figure 7

par(mfrow=c(2,2), mar = c(2, 4, 0, 0), oma = c(3, 3.5, 0.5, 0.5), tcl = -0.25, mgp = c(3, 1, 0))

hist(unique(forest_random$correlation), 
     las=1,xlim = c(-0.85,-0.3), col="lightgrey", 
     ylim=c(0,1000), 
     main = "", 
     xlab = "Pearson's r coefficient", 
     ylab = "",
     cex.axis = 1.5,
     cex.lab=1.5,
     breaks = 50)
text(c(-0.84), c(970), c("A"), cex= 2)
abline(v = -0.31, col="red", lwd = 3, xpd=F)
abline(v = -0.66, col="blue", lwd = 3, xpd=F)
mtext("Frequency", side = 2, line =4.9, adj = -0.35, cex = 2)

hist(unique(grassland_random$correlation), 
     las=1, 
     xlim = c(-0.3, 0.6), 
     col="lightgrey", 
     ylim=c(0,1000), 
     main = "", 
     xlab = "Pearson's r coefficient",
     ylab = "",
     cex.axis = 1.5,
     cex.lab=1.5,
     breaks = 50)
text(c(-0.29), c(970), c("B"), cex= 2)
abline(v = 0.51, col="red", lwd = 3, xpd=F)
abline(v = 0.19, col="blue", lwd = 3, xpd=F)

hist(unique(arable_random$correlation), 
     las=1, 
     xlim = c(-0.7, 0.2), 
     col="lightgrey", 
     ylim=c(0,1000), 
     main = "", 
     xlab = "Pearson's r coefficient", 
     ylab = "",
     cex.axis = 1.5,
     cex.lab=1.5,
     breaks = 50)
text(c(-0.69), c(970), c("C"), cex= 2)
abline(v = 0.085, col="red", lwd = 3, xpd=F)
abline(v = -0.42, col="blue", lwd = 3, xpd=F)

hist(unique(settlement_random$correlation), 
     las=1, 
     col="lightgrey", 
     ylim=c(0,1000), 
     main = "", 
     xlab = "Pearson's r coefficient", 
     ylab = "",
     cex.axis = 1.5,
     cex.lab=1.5,
     breaks = 50)
abline(v = -0.076, col="red", lwd = 3, xpd=F)
abline(v = -0.35, col="blue", lwd = 3, xpd=F)
text(c(-0.64), c(970), c("D"), cex= 2)

mtext("Pearson's r coefficient", side = 1, line =3.5, adj = -2, cex = 2)
