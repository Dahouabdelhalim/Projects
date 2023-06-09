#!/usr/bin/env Rscript
#' ----------------------------------------------------------------------------
#' - NAME:          reproducing_calculations.R
#' - AUTHOR:        Deborah Morgenstern
#' - DATE:          2022-11-23
#' ----------------------------------------------------------------------------
#' - OVERVIEW:      Reproducing the core analysis for the paper
#'                  "Thunderstorm Types in Europe" (Morgenstern et al.)
#' 
#' - DESCRIPTION:   - Load and inspect the data.
#'                  - Do a regional comparison using a principal component
#'                    analysis (PCA).
#'                  - Find thunderstorm types using k-means clustering.
#'                  - Label the found clusters as thunderstorm types using
#'                    a decision tree.
#'                  - Order the found thunderstorm types.
#'                  - Plot figures 2, 3, 5, 6, and 8.
#'
#' - DEMANDS:       Takes about half a minute to run.
#'                  If you have problems with this script, inspect my
#'                  sessionInfo() at the end.
#' ----------------------------------------------------------------------------

###############################################################################
#' 
###############################################################################

library(data.table)
library(colorspace)

# You might want to set the working directory manually
# setwd("path/to/the/data")

# For simplicity the algorithm is demonstrated at one domain and the plots
# consider only k = 3. To look at other domains or k change the values here.
id <- "A"  # A - L are available
k  <- 3    # k = 2 - 6 is what we used (3 is presented in the paper)

##############################################################################
# Introduction to the data
##############################################################################

# *************
# Get an overview over the data by looking at one domain

# Open data
one_domain <- read.csv(sprintf("data_domain%s.csv", id),
                       header = TRUE, row.names = 1,
                       strip.white = TRUE)

# Get an overview over the data
head(one_domain)

# This is the data for one representative sample in domain A
# Each line represents one lightning observation, i.e. ERA5 cell-hours where
# at least one flash was observed. 
# The first columns are meta-data (lon, lat, time, season).
# The last columns are the results from k-means clustering for k = 2-6. These
# values are reproduced later in this script (thunderstorm types).
# The ERA5 values for each observation are the 25 columns in the middle. 
# The column names refer to the short names of each variable.

# Which ERA5 variables are in the data set?
era5_vars <- names(one_domain[5:29])
era5_vars

# Details on the era5-variables can be found in era5_vars.csv
era5_details <- read.csv("era5_vars.csv",
                         header = TRUE, row.names = 1,
                         strip.white = TRUE)

# These are the more speaking names of the variables
era5_details$niceName

# Quick check if the ERA5 short names in the data and in era5_vars.csv
# are the same
stopifnot(era5_details$shortName == era5_vars)

# If you want to use more speaking names for the era5 variables in the plots,
# then use the other line.
# era5_lbl <- era5_details$niceName
era5_lbl <- era5_vars

# era5_vars.csv provides even more information
names(era5_details)

# Back to our data: each season is equally represented
table(one_domain$season)


##############################################################################
# Part 1: Reproducing Sec. 4.1: Regional analysis
##############################################################################

# In the regional analysis, the domain means are compared. 
# The domain means are provided in Table 1 of the paper but also in the
# supplements as means_domains.csv.
# This section reproduces means_domains.csv, computes the PCA, and reproduces
# Figures 2 and 3.

# *************
# Load domain means (original)
domain_means_original <- read.csv("means_domains.csv",
                                  header = TRUE, row.names = 1,
                                  strip.white = TRUE)

# *************
# Start reproduction (rep) by initializing an empty data.frame
domain_means_rep <- data.frame(matrix(nrow = 25, ncol = 12),
                               row.names = era5_vars)
colnames(domain_means_rep) <- LETTERS[1:12]

# Get the mean value for each domain
# This takes a few moments as the files for each domain are opened
for (idx in LETTERS[1:12]) {
  # Open data
  dat <- read.csv(paste0("data_domain", idx, ".csv"),
                  header = TRUE, row.names = 1,
                  strip.white = TRUE)
  # Calculate domain means
  domain_means_rep[[idx]] <- colMeans(dat[, 5:29])
}

# *************
# Compare the reproduction with the original

# There might be minor differences in precision. Omit them by rounding.
domMean_original_round     <- round(domain_means_original, digits = 7)
domMean_reproduction_round <- round(domain_means_rep, digits = 7)

for (idx in LETTERS[1:12]) {
  stopifnot(identical(domMean_original_round[[idx]],
                      domMean_reproduction_round[[idx]]))
}

# *****************************************************************************
# Transform and scale data (apply equations 1 and 2)
    
domMean_t          <- t(domain_means_original)                # transpose
domMean_trafo_t    <- sign(domMean_t) * sqrt(abs(domMean_t))  # transform
domMean_scaled_t   <- scale(domMean_trafo_t)                  # scale
domMean_scaled     <- t(domMean_scaled_t)                     # transpose back

# *****************************************************************************
# Calculate PCA
pca_domains <- prcomp(domMean_scaled_t)

# Get loadings
loadings <- pca_domains$rotation[, c(1:2)]

# How much variance is explained by the first two PC?
pc1_domains <- round(summary(pca_domains)$importance[2,1] * 100, digits = 2)
pc2_domains <- round(summary(pca_domains)$importance[2,2] * 100, digits = 2)


# *****************************************************************************
# Prepare plotting: Define colors
# *****************************************************************************

cols_region <- colorspace::qualitative_hcl(n = 4, h = c(-260, 190),
                                           c = 60, l = 70)
cols_region_named <- c("Mediterranean"  = cols_region[4],
                       "Alpine-central" = cols_region[3],
                       "Coastal"        = cols_region[2],
                       "Continental"    = cols_region[1])


# *****************************************************************************
# Plotting: PCA of domain means (Fig. 2)
# *****************************************************************************

# Scaling of loadings to have loadings and PC1 & 2 at the same axis
# Formula is based on getAnywhere(biplot.prcomp)
scale_factor <- 1.4
loadingsS <- loadings
lambda1 <- pca_domains$sdev[1] * sqrt(length(pca_domains$x[,1])) ^ scale_factor
lambda2 <- pca_domains$sdev[2] * sqrt(length(pca_domains$x[,2])) ^ scale_factor
loadingsS[,1] <- loadings[,1] * lambda1
loadingsS[,2] <- loadings[,2] * lambda2

# *************
# Plot

# Save figure as pdf
pdf("rep_figure_02.pdf", width = 11, height = 6)
par(mar = c(3, 3, 3, 3) + 2, xpd = FALSE)


# Empty plot
plot(pca_domains$x[,1], pca_domains$x[,2],
     type = "n", xlim = c(-8, 10), ylim = c(-6, 6), las = 1,
     xlab = paste("PC 1 (", pc1_domains, "%)"),
     ylab = paste("PC 2 (", pc2_domains, "%)"),
     main = "Principal component analysis of domain means")

# Lines at x = 0 and y = 0
abline(v = 0, col = "darkgray", lty = 2)
abline(h = 0, col = "darkgray", lty = 2)

# Add biplot / loadings
arrows(0, 0, x1 = loadingsS[,1], y1 = loadingsS[,2], length = .1,
       col = "lightsteelblue3")

# Add labels to arrows
# Many labels are overlapping. For the paper the position of the labels are
# manually tweaked.
text(loadingsS, labels = era5_lbl, cex = .9,
     pos = 4, offset = 0.3, col = "lightsteelblue4", xpd = TRUE)

# Add colored circles
points(pca_domains$x[,1], pca_domains$x[,2],
       col = c(rep(cols_region[2], 4),   # coastal:        A, B, C, D
               rep(cols_region[1], 3),   # continental:    E, F, G
               rep(cols_region[3], 2),   # alpine-central: H, I
               rep(cols_region[4], 3)),  # mediterranean:  J, K, L
       pch = c(rep(17, 4), rep(18, 3), rep(16, 2), rep(15, 3)),
       cex = c(rep(5, 4), rep(6.5, 3), rep(5.5, 2), rep(5, 3)),
       cex.lab = 1.7, cex.axis = 1.3, las = 1, lwd = 2)

# Add labels to colored circles
text(x = pca_domains$x[,1], y = pca_domains$x[,2],
     labels = rownames(pca_domains$x), offset = 0.2, col = "black", xpd = TRUE)

# Legend
legend("bottomleft", bty = "n", pch = 16, pt.cex = 2,
       legend = c("Coastal", "Continental", "Mediterranean", "Alpine-central"),
       col = cols_region_named[c("Coastal", "Continental",
                                 "Mediterranean", "Alpine-central")])
dev.off()


# *****************************************************************************
# Plotting: Parallel coordinate plot of domain means (Fig. 3)
# *****************************************************************************

# Calculate means for each region
region_means <- data.frame(
  "Coastal"       = rowMeans(domMean_scaled[ , c("A", "B", "C", "D")]),
  "Continental"   = rowMeans(domMean_scaled[ , c("E", "F", "G")]),
  "AlpineCentral" = rowMeans(domMean_scaled[ , c("H", "I")]),
  "Mediterranean" = rowMeans(domMean_scaled[ , c("J", "K", "L")])
  )

# *************
# Plot

# Save figure as pdf
pdf("rep_figure_03.pdf", width = 10, height = 6)
par(mar = c(7, 4, 5, 6), mgp = c(2, .6, 0), xpd = TRUE)

# Empty plot
matplot(region_means, type = "n", ylim = c(-1.9, 2.5), xaxt = "n", las = 1,
        cex.lab = 1.3, cex.axis = 1.3, ylab = "Scaled value",
        main = "Parallel coordinate plot: Region means")

# xlab
text(y = -2.2, x = seq(1, 25), xpd = TRUE, labels = era5_lbl,
     las = 1, srt = 45, adj = 1, cex = 1.2)

# Gray vertical lines
abline(v = seq(1, 25), col = "gray", lty = 1, xpd = FALSE)

# Horizontal 0 line
lines(x = c(par('usr')[1], par('usr')[2]), y = c(0, 0), lty = 1, lwd = 2)

# Plot
matplot(region_means, col = cols_region_named, type = "o", cex = 1.5,
        pch = seq(15, 18), lty = 1, lwd = 2.5, xaxt = "n", add = TRUE)

# Domain annotations
text(x = 26, y = c(1.3, 0.3, -0.5, -1), pos = 4, xpd = TRUE, cex = 1.3,
     labels = c("H, I", "E, F, G", "J, K, L", "A, B, C, D"),
     col = cols_region_named[c("Alpine-central", "Continental",
                               "Mediterranean", "Coastal")])

# legend
legend("top", inset = c(0, -.1), horiz = TRUE, xpd = TRUE, bty = "n",
       text.width = 4.6, pch = seq(15, 18), pt.lwd = 3,
       col = cols_region_named, text.col = cols_region_named,
       legend = names(cols_region_named))

dev.off()

# For simplicity, the category annotations are omitted


##############################################################################
# Part 2: Reproducing Sec. 4.2, 4.3, and 4.4: Thunderstorm types
##############################################################################

# Find thunderstorm types using k-means clustering.
# This results in a table of cluster means as provided in in Table 2 of the 
# paper and in the supplements as means_clusters.csv
# This section computes k-means clustering, labels the found clusters according 
# to the decision tree, reproduces means_clusters.csv, and reproduces
# Figures 5, 6, and 8.

# *************
# Load data (Table 2, original)
cluster_means_original <- read.csv("means_clusters.csv",
                                   header = TRUE, row.names = 1,
                                   strip.white = TRUE)

# Have a look at the data. Each domain is represented by the cluster means (rows)
# of the three thunderstorm types found there (columns). 
cluster_means_original

# In the following, we only reproduce the results for one domain.
# Subset the respective columns.
cm_original <- cluster_means_original[, grep(id, names(cluster_means_original))]

# Get the names of the thunderstorm types we want to find in the reproduction.
labels_original <- substr(colnames(cm_original), 2, 20)
labels_original

# *************
# Load data for one domain
obs <- read.csv(paste0("data_domain", id, ".csv"),
                header = TRUE, row.names = 1,
                strip.white = TRUE)

# Have a look at the data. Each row is one lightning observation, i.e. an ERA5
# cell-hour where at least one lightning stroke was observed. Columns are
# meta-data or names of ERA5 variables. The last six columns are the results
# from the cluster analysis presented in the paper. In the following, we
# reproduce one of these cluster analysis.
head(obs)

# Each season is equally represented in the data
table(obs$season)

# For the reproduction, only ERA5 variables are required
# (remove meta-data and results)
obs_era5 <- as.data.table(obs[, 5:29])

# *****************************************************************************
# Transformation and scaling

obs_trafo  <- sign(obs_era5) * sqrt(abs(obs_era5))  # transform
obs_scaled <- as.data.table(scale(obs_trafo))       # scale


# In domains without sea, the land-sea mask has the same value everywhere,
# which results in problems during scaling. When this occurs, the scaled
# value for land-sea mask is manually set to zero, i.e. an average value.
if (length(unique(obs_scaled$lsm_bin)) == 1) {
  obs_scaled$lsm_bin <- 0
}

# There should be no NAs in the data.
if (any(is.na(obs_scaled))) {
  warning("  * Careful! NA in data. \\n")
}

# *****************************************************************************
# Cluster analysis

# Calculate k-means clustering
set.seed(1)
clK <- kmeans(obs_scaled, centers = k, iter.max = 100,
              nstart = 150, algorithm = "MacQueen")

# *************
# Save results

# Get cluster means
cm_scaled   <- clK$centers
cm_unscaled <- setorder(obs_era5[, lapply(.SD, mean), by = clK$cluster], clK)
# Add cluster id to obs data
cl_name <- sprintf("rep_cluster_k%02d", k)
obs[[cl_name]] <- clK$cluster


# *****************************************************************************
# Labeling the found clusters with the corresponding thunderstorm type.
# i.e. describe the clusters meteorologically.
# *****************************************************************************

# To label the clusters, first the mean in each category is calculated. 
# Based on the resulting values, the thunderstorm types are assigned using the
# decision tree (Fig. 4). Then the clusters are ordered to have well organized 
# legends and to have always the same order within the barplots later.
# After all this, the reproduction is compared to the original values.

# Get the era5 variable categories
categories <- era5_details$meteorologicalCategory
names(categories) <- era5_vars
categories

# *************
get_increased_categories <- function(cl_mns, ctgr = categories,
                                     cls = k, thresh = 0.3) {
  # cls_mns = cluster means as returned from kmeans()
  # ctgr = categories; named vector of category-variable pairs
  # cls = number of clusters
  # thresh = threshold for considering a category for naming
  
  # Calculate category means per cluster, i.e.,
  # Meteorological characterization of clusters.
  # Figure out which category means are above or below the threshold.
  # This is the basis for labeling with the decision tree later.
  # Therefore, 1) Compute category means for each cluster.
  # 2) Label each cluster with all the increased categories.
  # e.g. if wind field and surface exchange are both above threshold the 
  # cluster is labeled wf_sfxexH.
  # If none is above it results in noLabelX (x is a number)
  
  # Take absolute values in tSfc1000. This variable is just a difference.
  cl_mns[,"tSfc1000"] <- abs(cl_mns[,"tSfc1000"])
  
  # Calculate means for all values in each category
  lbl        <- data.frame(matrix(ncol = 0, nrow = cls))
  lbl$mfH    <- rowMeans(cl_mns[, names(ctgr[ctgr == "mass field"])])
  lbl$mfL    <- lbl$mfH * -1
  lbl$wf     <- rowMeans(cl_mns[, names(ctgr[ctgr == "wind field"])])
  lbl$cp     <- rowMeans(cl_mns[, names(ctgr[ctgr == "cloud physics"])])
  lbl$sfcexH <- rowMeans(cl_mns[, names(ctgr[ctgr == "surface exchange"])])
  lbl$sfcexL <- lbl$sfcexH * -1
  
  # If the mean in a category is above threshold,
  # include the category name in the label
  # (assign label names given threshold).
  lbl_labeled <- lbl  # Prepare results
  nlab_id <- 1        # Number, in case its a noLabel
  for (cl in seq(1, cls)) {  # Loop over each cluster
    suppressWarnings(cl_sort <- sort(lbl[cl, ], decreasing = TRUE))
    idx <- which(cl_sort > thresh)
    nm <- names(cl_sort)[idx]
    if (length(nm) == 0) {
      nm <- paste0("noLabel", nlab_id)
      nlab_id <- nlab_id + 1
    }
    lbl_labeled$label[cl] <- paste(nm, collapse = '_')
  }
  
  # Create return vector
  lbs <- lbl_labeled$label
  names(lbs) <- as.character(seq(1, length(lbs)))
  return(lbs)
}

# *************
decision_tree <- function(incr_cat) {
  # incr_cat = result of get_increased_categories()
  # Apply decision tree to name the clusters properly
  
  named <- incr_cat
  
  for (i in seq(1, length(named))) {
    if (grepl("noLabel", named[i])) {
      next   
    }
    
    if (grepl("wf", named[i])) {
      if (grepl("cp", named[i])) {
        named[i] <- "wind_field_cp"
      } else {
        named[i] <- "wind_field"
      }
    } else if (grepl("mfL", named[i])) {
      named[i] <- "wind_field_noMF"
    } else if (grepl("sfcexH", named[i])) {
      named[i] <- "mass_field_sx"
    } else {
      named[i] <- "mass_field"
    }
  }
  
  return(named)
}

# *************
# Apply functions
increased_categories <- get_increased_categories(cl_mns = cm_scaled)
clusters_labeled <- decision_tree(incr_cat = increased_categories)

# The clusters have two names / labels now: Their original number obtained
# from kmeans() and a label obtained from decision_tree (the label from
# get_increased_categories() is no longer required).
# Both names (number and label) are stored in the variable clusters_labeled
clusters_labeled

# *****************************************************************************
# Ordering the clusters as follows:
# Mass field SX, mass field, wind_field noMF, wind_field, wind_field CP
# Purpose: Have always the same color ordering in the barplots and in
# the legends.

# Difficulty:
# Not every cluster is present in every domain. 
# At higher k, some cluster labels occur twice.

# *************
# First, order the available cluster labels.
# If there were always all five thunderstorm types present,
# we could easily use match().
clbl <- clusters_labeled

tmpBlueDark <- clbl[which(clbl == "wind_field_cp")]
clbl <- clbl[!clbl %in% tmpBlueDark]
tmpBlueMiddle <- clbl[which(clbl == "wind_field")]
clbl <- clbl[!clbl %in% tmpBlueMiddle]
tmpBlueLight <- clbl[which(clbl == "wind_field_noMF")]
clbl <- clbl[!clbl %in% tmpBlueLight]
tmpRedDark <- clbl[which(clbl == "mass_field_sx")]
tmpRedLight <- clbl[!clbl %in% tmpRedDark]

clusters_labeled_ord <- unique(c(tmpRedDark, tmpRedLight,
                               tmpBlueLight,
                               tmpBlueMiddle, tmpBlueDark))

# Have a look at the ordered cluster labels
clusters_labeled_ord

# *************
# Second, order the original cluster names (numbers).
# If k >= 4 there are sometimes two clusters in one domain that have the 
# same label according to the decision tree. If there were no duplicates,
# we could easily use match(). But now we need this loop.
clusters_ord <- c()

for (i in seq(1, length(clusters_labeled))) {
  pos <- which(clusters_labeled_ord[i] == unname(clusters_labeled))
  if (length(pos) == 1) {
    clusters_ord <- c(clusters_ord, pos)
  } else {
    for (j in seq(1, length(clusters_labeled))) {
      # if value is not yet in clusters_ord, add it!
      if (!pos[j] %in% clusters_ord) {
        if (is.na(pos[j])) {next}
        clusters_ord <- c(clusters_ord, pos[j])
      }
    }
  }
}

# This are the ordered cluster labels
labels_ord <- clusters_labeled[clusters_ord]
labels_ord

# *************
# Rename and order the cluster-means data (scaled and unscaled)
cm_scaled_ord           <- cm_scaled[clusters_ord,]
rownames(cm_scaled_ord) <- labels_ord

cm_unscaled_ord         <- cm_unscaled[clusters_ord,]
cm_unscaled_ord$clK     <- labels_ord

# *************
# Check if the reproduction is the same as the original

# Double check, if there are the same ERA5 variables
stopifnot(rownames(cm_original) == colnames(cm_scaled_ord))
stopifnot(era5_vars == colnames(cm_scaled_ord))

# Table 2 / cm_original contains only data for k = 3
if ( k == 3 ) {
  # Same thunderstorm types?
  stopifnot(labels_original == labels_ord)
  
  # Same mean values?
  # Round to 5 digits as there are sometimes difficulties with the precision
  cm_original_round <- as.matrix(round(cm_original, digits = 5))
  cm_unscaled_round <- round(t(cm_unscaled_ord[, 2:26]), digits = 5)
  # no column / row names
  colnames(cm_original_round) <- NULL
  rownames(cm_unscaled_round) <- NULL
  stopifnot(cm_original_round == cm_unscaled_round)
}


# *****************************************************************************
# Prepare plotting: Define colors
# *****************************************************************************

make_cols_type <- function(lbs) {
  # Returns a named vector of the hcl-colors for thunderstorm types
  
  all_col_h <- c()
  all_col_c <- c()
  all_col_l <- c()
  
  # -------------
  # define noLabel colors
  for (i in seq(1, length(lbs))) {
    
    if (grepl("noLabel", lbs[i])) {
      all_col_h <- c(all_col_h, 120)
      all_col_c <- c(all_col_c, 30)
      all_col_l <- c(all_col_l, 90)
      next
    }
    
    # -------------
    if (lbs[i] == "wind_field_cp") {
      # dark blue
      col_h <- 250
      col_c <- 45
      col_l <- 40
    } else if (lbs[i] == "wind_field") {
      # middle blue
      col_h <- 250
      col_c <- 45
      col_l <- 60
    } else if (lbs[i] == "wind_field_noMF") {
      # light blue
      col_h <- 250
      col_c <- 45
      col_l <- 80
    } else if (lbs[i] == "mass_field") {
      # light red
      col_h <- 10
      col_c <- 45
      col_l <- 70
    } else if (lbs[i] == "mass_field_sx") {
      # darkred
      col_h <- 10
      col_c <- 45
      col_l <- 40
    } else {
      stop("  * Error: Cluster name does not exist \\n")
    }
    
    all_col_h <- c(all_col_h, col_h)
    all_col_c <- c(all_col_c, col_c)
    all_col_l <- c(all_col_l, col_l)
  }
  
  cols <- c(hcl(all_col_h, all_col_c, all_col_l))
  names(cols) <- lbs
  
  return(cols)
}

# *************
# Define symbols

sym_types <- c("wind_field_cp" = 1,
               "wind_field" = 2,
               "wind_field_noMF" = 3,
               "mass_field" = 4,
               "mass_field_sx" = 5)


# *****************************************************************************
# Plotting: Parallel coordinates plot of cluster means (Fig. 5)
# *****************************************************************************

pdf("rep_figure_05.pdf", width = 7, height = 3.7)
par(mar = c(.5, .5, .5, .5), oma = c(5, 4, 4, 2),
    mgp = c(2, .6, 0), xpd = TRUE)

# Empty plot
matplot(t(cm_scaled_ord), type = "n", ylim = c(-1.5, 2), xaxt = "n", las = 1,
        ylab = "")

# Gray vertical lines
abline(v = seq(1, 25), col = "gray", lty = 1, xpd = FALSE)

# Horizontal 0 line
lines(x = c(par('usr')[1], par('usr')[2]), y = c(0, 0), lty = 1, lwd = 2)

# Horizontal 0.3 lines
lines(x = c(par('usr')[1], par('usr')[2]),
      y = c(0.3, 0.3), lty = 2, lwd = 2, col = "gray")
lines(x = c(par('usr')[1], par('usr')[2]),
      y = c(-0.3, -0.3), lty = 2, lwd = 2, col = "gray")

# Add curves
matplot(t(cm_scaled_ord), lwd = 2,
        col = make_cols_type(labels_ord),
        type = "o", pch = sym_types[labels_ord],
        cex = 1.3, lty = 1, xaxt = "n", ylab = "", add = TRUE)

# ylab
text(x = -3, y = 0, srt = 90, xpd = NA, label = "Scaled Value")

# xlab
text(y = -2, x = seq(1, 25), labels = era5_lbl,
     las = 1, srt = 45, adj = 1, xpd = NA)

# Legend
legend(x = 0, y = 3,text.width = c(4.6, 3.3, 5.1, 3.3),
  xpd = NA, horiz = TRUE, bty = "n", cex = .8,
  legend = c("Wind-field CP", "Wind-field", "Wind-field noMF",
             "Mass-field", "Mass-field SX"),
  pch = sym_types, pt.cex = 1.3,  pt.lwd = 3,
  col = make_cols_type(names(sym_types)))

# Title
title(sprintf("Parallel coordinate plot: Cluster means (Domain %s)", id),
      line = 1.8, outer = TRUE)

dev.off()

# For simplicity, the category annotations are omitted


# *****************************************************************************
# Plotting: Barplot, Seasonal variation of thunderstorm types (Fig. 6)
# *****************************************************************************

# *************
# Preparation

# Calculations
cluster_name <- sprintf("cluster_k%02d_label", k)
tab <- table(obs$season, obs[[cl_name]])
tab <- tab[c("winter", "spring", "summer", "fall"), clusters_ord]

# Look at the table
tab

# Remove column names and row names because we make our own labels and legend
colnames(tab) <- rep('', length(colnames(tab)))
rownames(tab) <- rep('', length(rownames(tab)))

# *************
# Plot

pdf("rep_figure_06.pdf", width = 6, height = 5)
par(mar = c(5, 2, 2, 9), xpd = TRUE, lheight = 0.8)

# Bars
mosaicplot(tab, color = make_cols_type(labels_ord),
    xlab = "", ylab = "", cex = 1.4, off = c(8, 0), las = 2, main = "")

# x labels
text(x = c(0.2, 0.45, 0.7, 0.95), y = -0.1,
     labels = c("Winter", "Spring", "Summer", "Fall"),
     srt = 45, cex = 1.2, pos = 2)

# y annotation
text(x = 0.05, y = c(0, .25, .5, .75, 1) - .045,
     labels = paste(as.character(seq(0, 1, 0.25)), " -"),
     cex = 1, pos = 2, xpd = NA)

# Legend
legend(x = 1.04, y = .7, xpd = NA, pch = 15, pt.cex = 2, bty = "n",
       col = make_cols_type(rev(names(sym_types))), 
       legend = c("Mass-field SX", "Mass-field", "Wind-field noMF",
                  "Wind-field", "Wind-field CP"))

# Title
title(sprintf("Seasonal variation of thunderstorm types (Domain %s)", id ),
      adj = 0.01)

dev.off()


# *****************************************************************************
# Plotting: PCA of cluster means (Fig. 8)
# *****************************************************************************

# *************
# Preparation

# Transform and scale
cm_all_t      <- t(cluster_means_original)             # transpose
cm_all_trafo  <- sign(cm_all_t) * sqrt(abs(cm_all_t))  # transform
cm_all_scaled <- scale(cm_all_trafo)                   # scale

# Get colors
col_cm_all <- make_cols_type(sub('.','', colnames(cluster_means_original)))

# Calculate PCA and get loadings
pca_cluMeans <- prcomp(cm_all_scaled)
loadings <- pca_cluMeans$rotation[, c(1:2)]

# Scale loadings to have them on the same axis as PC1 & 2
# Formula based on getAnywhere(biplot.prcomp)
scale_factor <- 1.2
loadingsS <- loadings
lambda1 <- pca_cluMeans$sdev[1] * sqrt(length(pca_cluMeans$x[,1])) ^ scale_factor
lambda2 <- pca_cluMeans$sdev[2] * sqrt(length(pca_cluMeans$x[,2])) ^ scale_factor
loadingsS[,1] <- loadings[,1] * lambda1
loadingsS[,2] <- loadings[,2] * lambda2

# How much variance is explained by the first two PC?
pc1_cluMeans <- round(summary(pca_cluMeans)$importance[2, 1] * 100, digits = 2)
pc2_cluMeans <- round(summary(pca_cluMeans)$importance[2, 2] * 100, digits = 2)

# *************
# Plot

pdf("rep_figure_08.pdf", width = 11, height = 6)
par(mar = c(3, 3, 3, 3) + 2, xpd = FALSE)

# Empty plot
plot(pca_cluMeans$x[,1], pca_cluMeans$x[,2],
     type = "n", xlim = c(-13, 15), ylim = c(-7, 8), las = 1,
     xlab = paste("PC 1 (", pc1_cluMeans, "%)"),
     ylab = paste("PC 2 (", pc2_cluMeans, "%)"),
     main = "Principal component analysis on all cluster means"
)

# Add loadings (arrows)
arrows(0, 0, x1 = loadingsS[,1], y1 = loadingsS[,2], length = .1,
       col = "darkgray")

# Add lines (x = 0 , y = 0)
abline(v = 0, col = "darkgray", lty = 2)
abline(h = 0, col = "darkgray", lty = 2)

# Add colored circles
points(pca_cluMeans$x[,1], pca_cluMeans$x[,2],
       col = adjustcolor(col_cm_all, alpha.f = 0.9),
       pch = 16, cex = 5, cex.lab = 1.7, cex.axis = 1.3, las = 1, lwd = 2)

# Label colored circles
text(x = pca_cluMeans$x[,1], y = pca_cluMeans$x[,2],
     labels = rep(LETTERS[1:12], each = 3), col = "black", xpd = TRUE)

# Label loadings (arrows)
text(loadingsS, labels = era5_lbl,
     cex = .9, pos = 4, offset = 0.3, col = "lightsteelblue4", xpd = TRUE)

# Legend
legend("bottomleft", bty = "n", pch = 16, pt.cex = 2,
       legend = c("Wind-field CP", "Wind-field", "Wind-field noMF",
                  "Mass-field", "Mass-field SX"),
       col = make_cols_type(c("wind_field_cp", "wind_field", "wind_field_noMF",
                              "mass_field", "mass_field_sx")))

dev.off()


##############################################################################
# Troubleshooting
##############################################################################

# In case you have problems with this script, compare your sessionInfo() with
# the one in which this script was written:

# R version 4.2.2 (2022-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 10 (buster)
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C            LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] data.table_1.14.2
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.8.3        compiler_4.2.2      pillar_1.8.1        class_7.3-20       
# [5] tools_4.2.2         digest_0.6.29       evaluate_0.16       lifecycle_1.0.1    
# [9] tibble_3.1.8        lattice_0.20-45     pkgconfig_2.0.3     rlang_1.0.4        
# [13] DBI_1.1.2          cli_3.3.0           rstudioapi_0.13     yaml_2.3.5         
# [17] xfun_0.31          fastmap_1.1.0       e1071_1.7-9         dplyr_1.0.9        
# [21] knitr_1.39         generics_0.1.2      vctrs_0.4.1         classInt_0.4-3     
# [25] grid_4.2.2         tidyselect_1.1.2    glue_1.6.2          sf_1.0-7           
# [29] R6_2.5.1           rnaturalearth_0.1.0 fansi_1.0.3         rmarkdown_2.14     
# [33] sp_1.4-7           purrr_0.3.4         magrittr_2.0.3      htmltools_0.5.2    
# [37] units_0.8-0        assertthat_0.2.1    colorspace_2.0-3    KernSmooth_2.23-20 
# [41] utf8_1.2.2         proxy_0.4-26   

# *************
# Uncomment to get your sessionInfo.
# sessionInfo()

