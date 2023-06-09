### NMDS and PERMANOVA for:
  # 1) phytochemical arsenals
  # 2) herbivore communities
  # 3) vegetation communities
  # 4) climate


### PHYTOCHEMICAL ARSENALS #####################################################
### 0. import chemistry data ---------------------------------------------------
setwd("~/Desktop/DRYAD/")
dat <- read.csv(file = "03_predicted_trait_values_01_07_2022.csv", header = T, stringsAsFactors = F )

# alphabetize by Population name
dat <- dat[order(dat$Population),]

# import Population plotting colors
col_df <- read.csv(file = "03_Population_plotting_colors.csv", header = T, stringsAsFactors = F)

# find disagreement between datasets
table(dat$Population == col_df$Population)

# merge data
dat2 <- merge(dat, col_df, by.x = "Population", by.y = "Population")

# remove column with population name for analyses
row.names(dat2) <- dat2$Population

# isolate chemistry data
    # unk10, calcA, conan, verb, calcB, mimulo, unk16
trait_dat <- dat2[, 
    c("UNK_10", "CALC_A", "CONAN", "VERB", "CALC_B", "MIMULO", "UNK_16")]


### 1. calculate Bray-Curtis distance ------------------------------------------
library(vegan)

# chemical arsenal distance
trait_CHEM_dist_bc <- vegdist(trait_dat, method = "bray")
trait_CHEM_dist_bc_matrix <- as.matrix(trait_CHEM_dist_bc)

# write.csv(trait_CHEM_dist_bc_matrix, file = "arsenal_CHEM_braycurtis_distance_matrix.csv")


### 01. conduct NMDS ordination ------------------------------------------------
chem_ord <- metaMDS(trait_dat, try = 10000)

stressplot(chem_ord) # non-metric R2 = 0.96

# extract ordination coordinates
chem_xy <- as.data.frame(chem_ord$points)

dev.new()
par(mfrow = c(2,2))
par(mar = c(2,5,2,2))

plot(chem_xy, pch = 21, col = NULL,
     xlab = "", ylab = "NMDS2", font.lab = 2,
     main = "Phytochemical arsenals",
     xlim = c(-0.25, 0.4), ylim = c(-0.25, 0.25))

abline(v = 0, col = "gray", lwd = 2, lty = 2)
abline(h = 0, col = "gray", lwd = 2, lty = 2)

# add hulls
ordered_group_cols <- c("#F8766D" , "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

ordihull(chem_ord, col_df$Subregion, draw = "line", label = F,
         col = ordered_group_cols, lwd = 2, lty = 1)

points(chem_xy, pch = 16, cex = 1, col = col_df$plot_color)
# text(chem_xy, labels = (dat$Population), pos = 3, font = 1, cex = 0.75)

legend("bottomright", 
      c("Coastal", "Cordilleran", "ENA", "Northern", "Southern", "UK"),
       col = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"),
       pch = 16, cex = 1, bty = "n")


### 02. conduct PERMANOVA / adonis ---------------------------------------------
    # to test whether populations differ in phytochemical composition

groups <- as.factor(dat2$Subregion)
dis <- vegdist(trait_dat, method = "bray")

# test homogeneity of group variances 
mod <- betadisper(dis, groups)
anova(mod) # p-value is not significant, therefore group dispersions are homogeneous
# TukeyHSD(mod)

# isolate group data
col_df$Population == row.names(trait_dat) # checks outs
group_df <- data.frame(Subregion = col_df$Subregion)

# conduct test
adonis2(trait_dat ~ Subregion, data = group_df, permutations = 10000)
    ### Biogeographic subregions differ significantly in chemical trait composition
    ### p = 0.04


### 03. conduct pairwise adonis ------------------------------------------------
# library(devtools)
# install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)

trait_output <- adonis.pair(vegdist(trait_dat), as.factor(col_df$Subregion))

trait_pvals  <-  matrix(nrow = 6, ncol = 6)
rownames(trait_pvals) <- as.character(sort(unique(col_df$Subregion)))
colnames(trait_pvals) <- as.character(sort(unique(col_df$Subregion)))

trait_pvals[lower.tri(trait_pvals)]  <-  trait_output$P.value
trait_pvals  <-  t(trait_pvals)

# write.csv(trait_pvals, file = "Table_SX_pairwise_phytochemical_differences.csv", row.names = T)



### HERBIVORE COMMUNITIES ######################################################
### 0. import herbivore data ---------------------------------------------------
herb_dat <- read.csv(file = "03_Rotter_Mimgut_herbivore_communities.csv", header = T, stringsAsFactors = F )

# add zero values where no richness was entered
herb_dat[is.na(herb_dat)] <- 0

# remove Population column from the data.frame
pops <- herb_dat[,1]
herb_dat <- herb_dat[, -1]

row.names(herb_dat) <- pops
dim(herb_dat) # 41 populations of Mimulus, abundance for 22 families of herbivores


### 1. calculate Bray-Curtis distance ------------------------------------------
herb_dist_bc <- vegdist(herb_dat, method = "bray")
herb_dist_bc_matrix <- as.matrix(herb_dist_bc)
row.names(herb_dist_bc_matrix) <- pops
colnames(herb_dist_bc_matrix) <- pops

# write.csv(herb_dist_bc_matrix, file = "herbivore_braycurtis_distance_matrix.csv")


### 01. conduct NMDS ordination ------------------------------------------------
herb_ord <- metaMDS(herb_dat, try = 10000) # stress = 0.18

# stressplot(herb_ord) # non-metric R2 0.97

# extract ordination coordinates
herb_xy <- as.data.frame(herb_ord$points)

par(mar = c(2,2,2,2))
plot(herb_xy, pch = 21, col = NULL,
     xlab = "", ylab = "", font.lab = 2,
     main = "Herbivore communities")

abline(v = 0, col = "gray", lwd = 2, lty = 2)
abline(h = 0, col = "gray", lwd = 2, lty = 2)

# add hulls
ordered_group_cols <- c("#F8766D" , "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

ordihull(herb_ord, col_df$Subregion, draw = "line", label = F,
         col = ordered_group_cols, lwd = 2, lty = 1)

points(herb_xy, pch = 16, cex = 1, col = col_df$plot_color)


### 02. conduct PERMANOVA / adonis ---------------------------------------------
# to test whether populations differ in phytochemical composition

groups <- as.factor(dat2$Subregion)
dis_herb <- vegdist(herb_dat, method = "bray")

# test homogeneity of group variances 
mod_herb <- betadisper(dis_herb, groups)
anova(mod_herb) # p-value is not significant, therefore group dispersions are homogeneous


# conduct PERMANOVA test
adonis2(herb_dat ~ Subregion, data = group_df, permutations = 10000)
### Biogeographic subregions differ significantly in herbivore community composition
### p < 0.001


### 03. conduct pairwise adonis ------------------------------------------------
herb_output <- adonis.pair(vegdist(herb_dat), as.factor(col_df$Subregion))

herb_pvals  <-  matrix(nrow = 6, ncol = 6)
rownames(herb_pvals) <- as.character(sort(unique(col_df$Subregion)))
colnames(herb_pvals) <- as.character(sort(unique(col_df$Subregion)))

herb_pvals[lower.tri(herb_pvals)]  <-  herb_output$P.value
herb_pvals  <-  t(herb_pvals)

# write.csv(herb_pvals, file = "Table_SX_pairwise_herbivore_differences.csv", row.names = T)



### VEGETATION COMMUNITIES #####################################################
### 0. import vegetation data ---------------------------------------------------
veg_dat <- read.csv(file = "03_Rotter_Mimgut_plant_communities.csv", header = T, stringsAsFactors = F )

# add zero values where no richness was entered
veg_dat[is.na(veg_dat)] <- 0

# remove Population column from the data.frame
veg_dat <- veg_dat[, -1]

row.names(veg_dat) <- pops
dim(veg_dat) # 41 populations of Mimulus, abundance of 50 plant families 


### 1. calculate Bray-Curtis distance ------------------------------------------
veg_dist_bc <- vegdist(veg_dat, method = "bray")
veg_dist_bc_matrix <- as.matrix(veg_dist_bc)
row.names(veg_dist_bc_matrix) <- pops
colnames(veg_dist_bc_matrix) <- pops

# write.csv(veg_dist_bc_matrix, file = "vegetation_braycurtis_distance_matrix.csv")


### 01. conduct NMDS ordination ------------------------------------------------
veg_ord <- metaMDS(veg_dat, try = 10000) # stress 0.24

stressplot(veg_ord) # non-metric R2 = 0.94

# extract ordination coordinates
veg_xy <- as.data.frame(veg_ord$points)

par(mar = c(5,5,2,2))
plot(veg_xy, pch = 21, col = NULL,
     xlab = "NMDS1", ylab = "NMDS2", font.lab = 2,
     main = "Vegetation communities")

abline(v = 0, col = "gray", lwd = 2, lty = 2)
abline(h = 0, col = "gray", lwd = 2, lty = 2)

# add hulls
ordered_group_cols <- c("#F8766D" , "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

ordihull(veg_ord, col_df$Subregion, draw = "line", label = F,
         col = ordered_group_cols, lwd = 2, lty = 1)

points(veg_xy, pch = 16, cex = 1, col = col_df$plot_color)


### 02. conduct PERMANOVA / adonis ---------------------------------------------
# to test whether populations differ in phytochemical composition

groups <- as.factor(dat2$Subregion)
dis_veg <- vegdist(veg_dat, method = "bray")

# test homogeneity of group variances 
mod_veg <- betadisper(dis_veg, groups)
anova(mod_veg) # p-value is not significant, therefore group dispersions are homogeneous


# conduct PERMANOVA test
adonis2(veg_dat ~ Subregion, data = group_df, permutations = 10000)
### Biogeographic subregions differ significantly in chemical trait composition
### p < 0.001


### 03. conduct pairwise adonis ------------------------------------------------
veg_output <- adonis.pair(vegdist(veg_dat), as.factor(col_df$Subregion))

veg_pvals  <-  matrix(nrow = 6, ncol = 6)
rownames(veg_pvals) <- as.character(sort(unique(col_df$Subregion)))
colnames(veg_pvals) <- as.character(sort(unique(col_df$Subregion)))

veg_pvals[lower.tri(veg_pvals)]  <-  veg_output$P.value
veg_pvals  <-  t(veg_pvals)

# write.csv(veg_pvals, file = "Table_SX_pairwise_vegetation_differences.csv", row.names = T)



### CLIMATE ####################################################################

# load spatial libraries
library(dismo)
library(rgdal)
library(raster)


### 01. download BIOCLIM data --------------------------------------------------
setwd("~/Desktop/DRYAD/")
mim_clim <- read.csv(file = "03_Mimgut_bioclim_data.csv", header = T, stringsAsFactors = F)

# convert first column to row.names
row.names(mim_clim) <- mim_clim$X
mim_clim <- mim_clim[, -1]

#  TEMPERATURE
# bio_1 BIO1 = Annual Mean Temperature (in units of degrees C * 10)
# bio_5 BIO5 = Max Temperature of Warmest Month
# bio_6 BIO6 = Min Temperature of Coldest Month
# bio_8 BIO8 = Mean Temperature of Wettest Quarter
# bio_9 BIO9 = Mean Temperature of Driest Quarter
# bio_10 BIO10 = Mean Temperature of Warmest Quarter
# bio_11 BIO11 = Mean Temperature of Coldest Quarter

# PRECIPITATION
# bio_12 BIO12 = Annual Precipitation (in millimeters)
# bio_13 BIO13 = Precipitation of Wettest Month
# bio_14 BIO14 = Precipitation of Driest Month
# bio_16 BIO16 = Precipitation of Wettest Quarter
# bio_17 BIO17 = Precipitation of Driest Quarter
# bio_18 BIO18 = Precipitation of Warmest Quarter
# bio_19 BIO19 = Precipitation of Coldest Quarter

# SEASONALITY
# bio_2 BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# bio_3 BIO3 = Isothermality (BIO2/BIO7) (* 100)
# bio_4 BIO4 = Temperature Seasonality (standard deviation *100)
# bio_7 BIO7 = Temperature Annual Range (BIO5-BIO6)
# bio_15 BIO15 = Precipitation Seasonality (Coefficient of Variation)

# add scalar to all bioclim variables that have negative values
apply(mim_clim2, MARGIN = 2, function(x) table(x<0)) # 1,6,8,9,11

mim_clim[,1] = mim_clim[,1] + abs(min(mim_clim[,1]))
mim_clim[,6] = mim_clim[,6] + abs(min(mim_clim[,6]))
mim_clim[,8] = mim_clim[,8] + abs(min(mim_clim[,8]))
mim_clim[,9] = mim_clim[,9] + abs(min(mim_clim[,9]))
mim_clim[,11] = mim_clim[,11] + abs(min(mim_clim[,11]))

# isolate climatic variables into groups
mim_clim_T <- mim_clim[, c("bio1", "bio5", "bio6", "bio8", "bio9", "bio10", "bio11")]
mim_clim_P <- mim_clim[, c("bio12", "bio13", "bio14", "bio16", "bio17", "bio18", "bio19")]
mim_clim_S <- mim_clim[, c("bio2", "bio3", "bio4", "bio7", "bio15")]


### 1. calculate Bray-Curtis distance ------------------------------------------

# overall climate
clim_dist_bc <- vegdist(mim_clim, method = "bray")
clim_dist_bc_matrix <- as.matrix(clim_dist_bc)
row.names(clim_dist_bc_matrix) <- pops
colnames(clim_dist_bc_matrix) <- pops

# write.csv(clim_dist_bc_matrix, file = "climate_braycurtis_distance_matrix.csv")


# temperature
clim_dist_bc_T <- vegdist(mim_clim_T, method = "bray")
clim_dist_bc_matrix_T <- as.matrix(clim_dist_bc_T)
row.names(clim_dist_bc_matrix_T) <- pops
colnames(clim_dist_bc_matrix_T) <- pops

# write.csv(clim_dist_bc_matrix_T, file = "climate_TEMP_braycurtis_distance_matrix.csv")


# precipitation
clim_dist_bc_P <- vegdist(mim_clim_P, method = "bray")
clim_dist_bc_matrix_P <- as.matrix(clim_dist_bc_P)
row.names(clim_dist_bc_matrix_P) <- pops
colnames(clim_dist_bc_matrix_P) <- pops

# write.csv(clim_dist_bc_matrix_P, file = "climate_PRECIP_braycurtis_distance_matrix.csv")


# seasonality
clim_dist_bc_S <- vegdist(mim_clim_S, method = "bray")
clim_dist_bc_matrix_S <- as.matrix(clim_dist_bc_S)
row.names(clim_dist_bc_matrix_S) <- pops
colnames(clim_dist_bc_matrix_S) <- pops

# write.csv(clim_dist_bc_matrix_S, file = "climate_SEASONALITY_braycurtis_distance_matrix.csv")


### 01. conduct NMDS ordination ------------------------------------------------
clim_ord <- metaMDS(mim_clim, try = 10000) # stress = 0.10

stressplot(clim_ord) # non-metric R2 = 0.99

# extract ordination coordinates
clim_xy <- as.data.frame(clim_ord$points)

par(mar = c(5,2,2,2))
plot(clim_xy, pch = 21, col = NULL,
     xlab = "NMDS1", ylab = "", font.lab = 2,
     main = "Climate space")

abline(v = 0, col = "gray", lwd = 2, lty = 2)
abline(h = 0, col = "gray", lwd = 2, lty = 2)

# add hulls
ordered_group_cols <- c("#F8766D" , "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

ordihull(clim_ord, col_df$Subregion, draw = "line", label = F,
         col = ordered_group_cols, lwd = 2, lty = 1)

points(clim_xy, pch = 16, cex = 1, col = col_df$plot_color)



### 02. conduct PERMANOVA / adonis ---------------------------------------------
# to test whether populations differ in phytochemical composition

groups <- as.factor(dat2$Subregion)
dis_clim <- vegdist(mim_clim, method = "bray")

# test homogeneity of group variances 
mod_clim <- betadisper(dis_clim, groups)
anova(mod_clim) # p-value is not significant, therefore group dispersions are homogeneous


# conduct PERMANOVA test
adonis2(mim_clim ~ Subregion, data = group_df, permutations = 10000)
### Biogeographic subregions differ significantly in climates
### p < 0.001


### 03. conduct pairwise adonis ------------------------------------------------
clim_output <- adonis.pair(vegdist(mim_clim), as.factor(col_df$Subregion))

clim_pvals  <-  matrix(nrow = 6, ncol = 6)
rownames(clim_pvals) <- as.character(sort(unique(col_df$Subregion)))
colnames(clim_pvals) <- as.character(sort(unique(col_df$Subregion)))

clim_pvals[lower.tri(clim_pvals)]  <-  clim_output$P.value
clim_pvals  <-  t(clim_pvals)

# write.csv(clim_pvals, file = "Table_SX_pairwise_climate_differences.csv", row.names = T)



### 04. Plot Figure 3 --------------------------------------------------------------
dev.new(width = 7.5, height = 7.5, noRStudioGD = TRUE)
par(mfrow = c(2,2))
par(mar = c(2,5,2,2))


# Phytochemical arsenals
plot(chem_xy, pch = 21, col = NULL,
     xlab = "", ylab = "NMDS2", font.lab = 2,
     main = "Phytochemical arsenals",
     xlim = c(-0.25, 0.4), ylim = c(-0.25, 0.25))

abline(v = 0, col = "gray", lwd = 2, lty = 2)
abline(h = 0, col = "gray", lwd = 2, lty = 2)

# add hulls
ordered_group_cols <- c("#F8766D" , "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

ordihull(chem_ord, col_df$Subregion, draw = "line", label = F,
         col = ordered_group_cols, lwd = 2, lty = 1)

points(chem_xy, pch = 16, cex = 1, col = col_df$plot_color)
# text(chem_xy, labels = (dat$Population), pos = 3, font = 1, cex = 0.75)

legend(x = 0.15, y = -0.075,
       c("Coastal", "Cordilleran", "ENA", "Northern", "Southern", "UK"),
       col = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"),
       pch = 16, cex = 1, bty = "n")


# Herbivores
par(mar = c(2,2,2,2))
plot(herb_xy, pch = 21, col = NULL,
     xlab = "", ylab = "", font.lab = 2,
     main = "Herbivore communities")

abline(v = 0, col = "gray", lwd = 2, lty = 2)
abline(h = 0, col = "gray", lwd = 2, lty = 2)

# add hulls
ordered_group_cols <- c("#F8766D" , "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

ordihull(herb_ord, col_df$Subregion, draw = "line", label = F,
         col = ordered_group_cols, lwd = 2, lty = 1)

points(herb_xy, pch = 16, cex = 1, col = col_df$plot_color)


# Vegetation communities
par(mar = c(5,5,2,2))
plot(veg_xy, pch = 21, col = NULL,
     xlab = "NMDS1", ylab = "NMDS2", font.lab = 2,
     main = "Vegetation communities")

abline(v = 0, col = "gray", lwd = 2, lty = 2)
abline(h = 0, col = "gray", lwd = 2, lty = 2)

# add hulls
ordered_group_cols <- c("#F8766D" , "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

ordihull(veg_ord, col_df$Subregion, draw = "line", label = F,
         col = ordered_group_cols, lwd = 2, lty = 1)

points(veg_xy, pch = 16, cex = 1, col = col_df$plot_color)


# Climate
par(mar = c(5,2,2,2))
plot(clim_xy, pch = 21, col = NULL,
     xlab = "NMDS1", ylab = "", font.lab = 2,
     main = "Climate space")

abline(v = 0, col = "gray", lwd = 2, lty = 2)
abline(h = 0, col = "gray", lwd = 2, lty = 2)

# add hulls
ordered_group_cols <- c("#F8766D" , "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

ordihull(clim_ord, col_df$Subregion, draw = "line", label = F,
         col = ordered_group_cols, lwd = 2, lty = 1)

points(clim_xy, pch = 16, cex = 1, col = col_df$plot_color)

# quartz.save(file = "Figure3_4panel_NMDS.jpg", type = "jpg", dpi = 300)










