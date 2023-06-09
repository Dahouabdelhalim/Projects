library(rgdal)
library(dismo)
library(raster)
library(rgeos)
library(sp)
library(scales)


### 00. import shapefiles (in .Rdata format) -----------------------------------
setwd("~/Documents/Post_Doc/0_Plantago_Research/0_MANUSCRIPTS/01_seed_size_manuscript/02_PUBLICATION_DOCUMENTS/Dryad/")

# import Four Corners shapefile
FourCorners <- get(load(file = "FourCorners_shape.Rdata"))

# import Colorado Plateau shapefile
CP <- get(load(file = "Colorado_Plateau_shape.Rdata"))
CP <-  aggregate(CP, dissolve = T)

# import georeferenced records of Plantago patagonica from the Colorado Plateau from GBIF and SEINET
Ppat_cp <- get(load(file = "Ppat_records_GBIF_and_SEINET_Colorado_Plateau.Rdata"))

# import study site locations
source_dat  <-  read.csv(file = "Plantago_collection_locations.csv", header = T, stringsAsFactors = F)

pops_xy  <-  source_dat[, c("Lon", "Lat")]
pops_xy  <-  SpatialPoints(pops_xy,)

crs(pops_xy) <-  crs(Ppat_cp)


### 01. import climate data ----------------------------------------------------
Ppat_clim <- read.csv(file = "Plantago_WorldClim_climate_data.csv", header = T)


### 02. conduct PCA of Bioclim/WorldClim variables -----------------------------
library(factoextra)

# combine focal Plantago populations with all other Plantago populations
ALL_Ppat <- rbind(pops_xy, Ppat_cp)
      # first 12 records seed collection locations
      # remaining records are collections of P.patagonica from the Colorado Plateau

# pca
res.pca <- prcomp(Ppat_clim, scale = TRUE)

# plot variables
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

  # Dim1 46.1%
  # Dim2 29.6%

# isolate individual PCA plotting coordinates
res.ind <- get_pca_ind(res.pca)
coords <-  data.frame(x = res.ind$coord[,1], y = res.ind$coord[,2])

# subset source location sites
seed_coords <- coords[1:12,]

# add plotting color
seed_coords$color  <- rep("black", times = 12)


### 03. plot Figure 1 (with aridity gradient) ----------------------------------
dev.new()
par(mfrow = c(1,2))

# collection panel
plot(FourCorners)
plot(CP, add = T, border = "gray")
plot(Ppat_cp, add = T, pch = 16, cex = 0.75, col = alpha("gray", alpha = 0.45))

# define aridity colors for each population
seed_coords$Population <- source_dat$Population

# import Aridity Index for all focal populations
arid_dat <- read.csv(file = "Plantago_climate_data_simple.csv", header = T, stringsAsFactors = F)
order(arid_dat$AI)

# generate a continuous color palette and add aridity gradient
rbPal <- colorRampPalette(c("blue", "red"))

arid_dat2  <-  arid_dat[order(arid_dat$AI),]
arid_dat2$colors  <-  rbPal(12)[as.numeric(cut(arid_dat2$AI, breaks = 12))]
arid_dat3  <-  merge(arid_dat2, seed_coords, by.x = "Population", by.y = "Population")

# plot seed collection locations
plot(pops_xy, col = alpha(arid_dat3$colors, alpha = 0.95), pch = 16, cex = 1.5, add = T)

# add aridity gradient
legend_colors = rbPal(100)
legend_image <- as.raster(matrix(rev(legend_colors), ncol = 1))

rasterImage(legend_image,
            xleft = -105.5,
            xright = -104.5,
            ybottom = 32.5,
            ytop = 36.5)

text(x = -104, y = c(32.5, 36.5), labels = c(65, 90), font = 1)


# PCA panel
plot(x = coords$x, y = coords$y, col = NULL, cex = 1, pch = 16,
     xlab = "PC1 (46.1%)", ylab = "PC2 (29.6%)", 
     main = "")

abline(v = 0, col = "black", lwd = 1, lty = 2)
abline(h = 0, col = "black", lwd = 1, lty = 2)

points(x = coords$x, y = coords$y, col = alpha("gray", alpha = 0.5), pch = 16, cex = 1)
points(x = seed_coords$x, y = seed_coords$y, col = alpha(arid_dat3$colors, alpha = 0.95),  pch = 16, cex = 1.75)

rasterImage(legend_image,
            xleft = 7.5,
            xright = 9,
            ybottom = 1.5,
            ytop = 5.5)

text(x = 9.75, y = c(1.5, 5.5), labels = c(65, 90), font = 1)

quartz.save(file = "Figure1_collection_locations_and_climate_PCA.jpg", type = "jpg", dpi = 300)





