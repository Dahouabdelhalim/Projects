### calculate association between predictors 
## and phytochemical resistance arsenal similarity

# import data
setwd("~/Desktop/DRYAD//")
list.files(pattern = ".csv")

# all import matrices, with row and column names
geo_dist <- read.csv(file = "04_geographic_distance_matrix.csv" , header = F, stringsAsFactors = F)

herb_dist <- read.csv(file = "04_herbivore_braycurtis_distance_matrix.csv" , header = F, stringsAsFactors = F)
veg_dist <- read.csv(file = "04_vegetation_braycurtis_distance_matrix.csv" , header = F, stringsAsFactors = F)

clim_dist_ALL <- read.csv(file = "04_climate_braycurtis_distance_matrix.csv", header = F, stringsAsFactors = F)
clim_dist_TEMP <- read.csv(file = "04_climate_TEMP_braycurtis_distance_matrix.csv", header = F, stringsAsFactors = F)
clim_dist_PRECIP <- read.csv(file = "04_climate_PRECIP_braycurtis_distance_matrix.csv", header = F, stringsAsFactors = F)
clim_dist_SEASONALITY <- read.csv(file = "04_climate_SEASONALITY_braycurtis_distance_matrix.csv", header = F, stringsAsFactors = F)

trait_dist <- read.csv(file = "04_phytochemical_resistance_arsenal_braycurtis_distance_matrix.csv" , header = F, stringsAsFactors = F)


# remove row and column names
geo_dist1 <- geo_dist[,-1]
geo_dist1 <- geo_dist1[-1,]

herb_dist1 <- herb_dist[,-1]
herb_dist1 <- herb_dist1[-1,]

veg_dist1 <- veg_dist[,-1]
veg_dist1 <- veg_dist1[-1,]

clim_dist_ALL_1 <- clim_dist_ALL[,-1]
clim_dist_ALL_1 <- clim_dist_ALL_1[-1,]

clim_dist_TEMP_1 <- clim_dist_TEMP[,-1]
clim_dist_TEMP_1 <- clim_dist_TEMP_1[-1,]

clim_dist_PRECIP_1 <- clim_dist_PRECIP[,-1]
clim_dist_PRECIP_1 <- clim_dist_PRECIP_1[-1,]

clim_dist_SEASONALITY_1 <- clim_dist_SEASONALITY[,-1]
clim_dist_SEASONALITY_1 <- clim_dist_SEASONALITY_1[-1,]

trait_dist <- trait_dist[,-1]
trait_dist_1 <- trait_dist[-1,]

# convert to "distance" objects in R
geo_dist2  <-  as.dist(geo_dist1)
herb_dist2  <-  as.dist(herb_dist1)
veg_dist2  <-  as.dist(veg_dist1)

clim_dist_ALL_2  <-  as.dist(clim_dist_ALL_1)
clim_dist_TEMP_2  <-  as.dist(clim_dist_TEMP_1)
clim_dist_PRECIP_2  <-  as.dist(clim_dist_PRECIP_1)
clim_dist_SEASONALITY_2  <-  as.dist(clim_dist_SEASONALITY_1)

trait_dist_2  <-  as.dist(trait_dist_1)

# check format of all matrices
str(geo_dist2)
str(herb_dist2)
str(veg_dist2)

str(clim_dist_ALL_2)
str(clim_dist_TEMP_2)
str(clim_dist_PRECIP_2)
str(clim_dist_SEASONALITY_2)

str(trait_dist_2)


### 03. conduct mantel tests ---------------------------------------------------
library(vegan)
library(scales)
library(ggeffects)


### phytochemical arsenal similarity ~ predictor variables

mantel(ydis = trait_dist_2, xdis = geo_dist2, method = "pearson", permutations = 10000)
mantel(ydis = trait_dist_2, xdis = herb_dist2, method = "pearson", permutations = 10000)
mantel(ydis = trait_dist_2, xdis = veg_dist2, method = "pearson", permutations = 10000)

mantel(ydis = trait_dist_2, xdis = clim_dist_TEMP_2, method = "pearson", permutations = 10000)
mantel(ydis = trait_dist_2, xdis = clim_dist_PRECIP_2, method = "pearson", permutations = 10000)
mantel(ydis = trait_dist_2, xdis = clim_dist_SEASONALITY_2, method = "pearson", permutations = 10000)


# plot - FIGURE 4
dev.new(width = 8, height = 5, noRStudioGD = TRUE)
par(mfrow = c(2,3))
par(mar = c(4,1,1,1))
par(oma = c(2,5,2,2))

### physical distance ----------------------------------------------------------
geo_dist3 <- geo_dist2/ (1000*1000) # 1,000s of km

mantel(ydis = trait_dist_2, xdis = geo_dist3, method = "pearson", permutations = 10000)

plot(trait_dist_2 ~ geo_dist3, pch = 16, cex = 1.5, col = alpha("gray", 0.35),
     ylab = "Phytochemical arsenal distance", xlab = "Geographic distance (1,000s of km)",
     font.lab = 2, cex.lab = 1.25)

m1 <- lm(trait_dist_2 ~ geo_dist3)
summary(m1)
abline(m1, lwd = 2, lty = 2)

text(x = 2.5, y = 0.15, "Mantel r = 0.04", cex = 1.25)
text(x = 2.5, y = 0.075, "p = 0.08", cex = 1.25)


### herbivore distance ---------------------------------------------------------
mantel(ydis = trait_dist_2, xdis = herb_dist2, method = "pearson", permutations = 10000)

plot(trait_dist_2 ~ herb_dist2, pch = 16, cex = 1.5, col = alpha("gray", 0.35), 
     ylab = "", xlab = "Herbivore community distance",
     yaxt = "n", font.lab = 2, cex.lab = 1.25)
axis(side = 2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels = F)

m2 <- lm(trait_dist_2 ~ herb_dist2)
summary(m2)
abline(m2, lwd = 2, lty = 2, col = NULL)

text(x = 0.35, y = 0.15, "Mantel r = -0.06", cex = 1.25)
text(x = 0.35, y = 0.075, "p = 0.81", cex = 1.25)


### vegetation distance --------------------------------------------------------
mantel(ydis = trait_dist_2, xdis = veg_dist2, method = "pearson", permutations = 10000)

plot(trait_dist_2 ~ veg_dist2, pch = 16, cex = 1.5, col = alpha("gray", 0.35), 
     ylab = "", xlab = "Vegetation community distance",
     yaxt = "n", font.lab = 2, cex.lab = 1.25)
axis(side = 2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels = F)

m3 <- lm(trait_dist_2 ~ veg_dist2)
summary(m3)
abline(m3, lwd = 2, lty = 2, col = NULL)

text(x = 0.45, y = 0.15, "Mantel r = 0.03", cex = 1.25)
text(x = 0.45, y = 0.075, "p = 0.36", cex = 1.25)


### climate (TEMP) distance ----------------------------------------------------
mantel(ydis = trait_dist_2, xdis = clim_dist_TEMP_2, method = "pearson", permutations = 10000)

plot(trait_dist_2 ~ clim_dist_TEMP_2, pch = 16, cex = 1.5, col = alpha("gray", 0.35), 
     ylab = "Phytochemical arsenal distance", xlab = "Climate (temperature) distance",
     font.lab = 2, cex.lab = 1.25)

m4 <- lm(trait_dist_2 ~ clim_dist_TEMP_2)
summary(m4)
abline(m4, lwd = 2, lty = 1)

text(x = 0.15, y = 0.15, "Mantel r = 0.25", cex = 1.25)
text(x = 0.15, y = 0.075, "p < 0.01", cex = 1.25)


### climate (PRECIP) distance --------------------------------------------------
mantel(ydis = trait_dist_2, xdis = clim_dist_PRECIP_2, method = "pearson", permutations = 10000)

plot(trait_dist_2 ~ clim_dist_PRECIP_2, pch = 16, cex = 1.5, col = alpha("gray", 0.35), 
     ylab = "", xlab = "Climate (precipitation) distance", yaxt = "n",
     font.lab = 2, cex.lab = 1.25)
axis(side = 2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels = F)

m5 <- lm(trait_dist_2 ~ clim_dist_PRECIP_2)
summary(m5)
abline(m5, lwd = 2, lty = 1, col = NULL)

text(x = 0.25, y = 0.15, "Mantel r = 0.01", cex = 1.25)
text(x = 0.25, y = 0.075, "p = 0.45", cex = 1.25)


### climate (SEASONALITY) distance ---------------------------------------------
mantel(ydis = trait_dist_2, xdis = clim_dist_SEASONALITY_2, method = "pearson", permutations = 10000)

plot(trait_dist_2 ~ clim_dist_SEASONALITY_2, pch = 16, cex = 1.5, col = alpha("gray", 0.35), 
     ylab = "", xlab = "Climate (seasonality) distance", yaxt = "n",
     font.lab = 2, cex.lab = 1.25)
axis(side = 2, at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), labels = F)

m6 <- lm(trait_dist_2 ~ clim_dist_SEASONALITY_2)
summary(m6)
abline(m6, lwd = 2, lty = 2)

text(x = 0.15, y = 0.15, "Mantel r = 0.07", cex = 1.25)
text(x = 0.15, y = 0.075, "p = 0.10", cex = 1.25)

# add y-axis title
mtext(text = "Phytochemical arsenal distance", side = 2, line = 2, font = 2, outer = TRUE)

# quartz.save(file = "Figure4_6panel_mantel.jpg", type = "jpg", dpi = 300)


