#' ---
#' title: "Differentiating lightning in winter and summer with characteristics of wind-field and mass-field"
#' author: "Deborah Morgenstern"
#' date: "2021-10-12"
#' ---
#' * NAME:                morgenstern2021_data_analysis.R
#' * SHORT DESCRIPTION:   Reproducing the core analysis and Figs. 2-4
#' 
#' This script reproduces the core analysis of the paper "Differentiating
#' lightning in winter and summer with characteristics of wind-field and
#' mass-field".
#' Everything is done for 1) the data presented in the paper and
#' 2) for a full year including transitional seasons.
#' Data is loaded, sqrt-transformed, and scaled.
#' A k-means cluster analysis and a PCA is performed. 
#' Then, the figures 2-4 of the paper are reproduced.

#' ### Load data

# Data of the representative sample which is presented in the paper.
dat <- read.csv("morgenstern2021_data.csv")
# Data of one representative sample for the whole year.
dat_allYear <- read.csv("morgenstern2021_data_allYear.csv")
# Columns containing ERA5 variables:
cov <- names(dat[7:41])

#' ### Overview

#' The first three columns (time, latitude, longitude) give information on
#' the lightning observations (cell-hours) rounded to the era5 resolution of
#' 0.25 x 0.25 degree and one hour (time is rounded down to the last full hour).
#' 
#' The column "scenario" is one of winter lightning (wl),
#' winter no lightning (wnol), summer lightning (sl),
#' and summer no lightning (snol). 
#' In dat_allYear additionally: spring lightning (spl),
#' spring no lightning (spnol), fall lightning (fl), fall no lightning (fnol).
#' For wl, all observations between 2010 - 2019
#' in the observational region are included. The others are sampled.
#' 
#' The column "cluster" is the result from the k-means cluster analysis
#' presented in the paper. It is reproduced in this script.
#' 
#' The remaining columns are era5 observations.
#' More details are given in the corresponding variable description file.

head(dat)

#' ### Transformation

#' A square root transformation is applied before scaling to reduce the skewness.
#' Three variables are distributed over positive and negative values. There the
#' transformation is applied to the absolute value multiplied by its sign.
#' Other negative values are set to zero. e.g. if ERA5 contains small negative
#' precipitation values.

transforming <- function(unscaled, era5_vars) {
  
  trafo <- subset(unscaled, select = era5_vars)
  
  # Deal with the variables that have meaningful positive and negative values.
  posneg <- c("viiwd", "mvimc", "sshf_up")
  for (i in posneg) trafo[[i]] <- sign(trafo[[i]]) * sqrt(abs(trafo[[i]]))
  
  # Apply square root transformation to the others.
  # Set negative values zero.
  for (i in setdiff(names(trafo), posneg)) trafo[[i]] <- sqrt(pmax(trafo[[i]], 0)) 
  
  return(trafo)
}

#' Apply transformation
dat_trafo <- transforming(unscaled = dat, era5_vars = cov)
dat_trafo_allYear <- transforming(unscaled = dat_allYear, era5_vars = cov)

#' ### Scaling
#' Mean and standard deviation is obtained from the scenarios without lightning
#' only (snol, wnol) because they dominate naturally.

scaling <- function(unscaled, trafo, nol = c("snol", "wnol")) {
  
  ms <- scale(trafo[unscaled$scenario %in% nol,])
  light_scaled <- scale(trafo,
                        center = attr(ms, "scaled:center"),
                        scale = attr(ms, "scaled:scale"))
  light_scaled <- as.data.frame(light_scaled)
  
  return(light_scaled)
}

#' Apply scaling
dat_scaled <- scaling(unscaled = dat, trafo = dat_trafo, nol = c("snol", "wnol"))
dat_scaled_allYear <- scaling(unscaled = dat_allYear, trafo = dat_trafo_allYear,
                              nol = c("snol", "wnol", "spnol", "fnol"))

#' ### Cluster Analysis

# Perform the cluster analysis
set.seed(as.integer(23))
clK <- kmeans(dat_scaled, centers = 5, iter.max = 100,
              nstart = 150, algorithm = "MacQueen")
set.seed(as.integer(9))
clK_allYear <- kmeans(dat_scaled_allYear, centers = 5, iter.max = 100,
                      nstart = 150, algorithm = "MacQueen")

#' Add the results to our data

tab <- proportions(table(clK$cluster, dat$scenario), 1)

# The cluster names origin from the plots later in this script.
dat$cluster_new[clK$cluster == 1] <- "wind_field"
dat$cluster_new[clK$cluster == 2] <- "CAPE"
dat$cluster_new[clK$cluster == 3] <- "average"
dat$cluster_new[clK$cluster == 4] <- "cloud_physics_CAPE" 
dat$cluster_new[clK$cluster == 5] <- "cloud_physics_wind_field" 

# Same for all_Year
tab_allYear <- proportions(table(dat_allYear$scenario, clK_allYear$cluster), 1)

dat_allYear$cluster_new[clK_allYear$cluster == 1] <- "cloud_physics_wind_field"
dat_allYear$cluster_new[clK_allYear$cluster == 2] <- "average"
dat_allYear$cluster_new[clK_allYear$cluster == 3] <- "cloud_physics_CAPE"
dat_allYear$cluster_new[clK_allYear$cluster == 4] <- "wind_field" 
dat_allYear$cluster_new[clK_allYear$cluster == 5] <- "CAPE" 

#' Check if the new cluster analysis has the same results as in the manuscript
same <- identical(dat$cluster, dat$cluster_new)
if (!same) {warning(
  "  * Results differ from the results presented in the paper.")}

same_allYear <- identical(dat_allYear$cluster, dat_allYear$cluster_new)
if (!same) {warning(
  "  * Results for the allYear analysis differ.")}

#' ### Table 2 (unscaled cluster medians)
clu_medians_unscaled   <- aggregate(dat[, cov], list(dat$cluster), median)
names(clu_medians_unscaled)[1] <- "cluster"
clu_medians_unscaled_allYear   <- aggregate(dat_allYear[, cov],
                                            list(dat_allYear$cluster), median)
names(clu_medians_unscaled_allYear)[1] <- "cluster"

#' Have a look at the table. Units are described in Table 1.
clu_medians_unscaled


#' ### PCA
# Calculate PCA (pcaL = pca on lightning data).
pcaL <- prcomp(dat_scaled)
pcaL_allYear <- prcomp(dat_scaled_allYear)

# To get loadings and PCA on the same axis in the plot,
# scaling is perfomed on the loadings.
# The formula based on getAnywhere(biplot.prcomp)
calculate_loadings <- function(pcaL, scale_factor = 0.5) {
  
  # Get loadings
  loadings <- pcaL$rotation[, 1:2]
  
  # Scaling as in biplot to have loadings and PC1 & 2 at the same axis
  # Formula based on getAnywhere(biplot.prcomp)
  loadingsS <- loadings
  lambda1 <- pcaL$sdev[1] * sqrt(length(pcaL$x[,1])) ^ scale_factor
  lambda2 <- pcaL$sdev[2] * sqrt(length(pcaL$x[,2])) ^ scale_factor
  loadingsS[,1] <- loadings[,1] * lambda1
  loadingsS[,2] <- loadings[,2] * lambda2
  
  return(loadingsS)
}

#' Apply function
loadingsS <- calculate_loadings(pcaL = pcaL)
loadingsS_allYear <- calculate_loadings(pcaL = pcaL_allYear)

#' # Plotting

#' ### Define features for the plots

cluster_names <- c("cloud_physics_wind_field", "wind_field", "average",
                   "CAPE", "cloud_physics_CAPE")

# Color definition
cols <- rev(hcl(c(10, 10, 80, 250, 250), c(80, 100, 90, 100, 80),
                c(20, 60, 90, 60, 20)))
names(cols) <- cluster_names

# Color definition light
cols_mosaic <- rev(hcl(c(10, 10, 80, 250, 250), c(70, 60, 50, 60, 70),
                       c(20, 60, 90, 60, 20)))
names(cols_mosaic) <- cluster_names 

# Plotting symbols
syms <- c("cloud_physics_wind_field" = 19, "wind_field" = 19, "average" = 18,
          "CAPE" = 17, "cloud_physics_CAPE" = 17)

# Legend labels
legend_labels <- c("Cloud physics & Wind field ", "Wind field", 
                   "Average", "CAPE",
                   "Cloud physics & CAPE")
# Legend labels with line breaks
legend_labels_2lines <- legend_labels
legend_labels_2lines[c(1, 5)] <- c("Cloud physics\\n& Wind field",
                                   "Cloud physics\\n& CAPE")
names(legend_labels_2lines) <- cluster_names

# Long version of variable names. Useful for labels.
names_plot <- c("CAPE", "CIN > 0", "-10 C isotherm height agl.",
                "Mean sea level pressure", "Wind direction at 10 m",
                "Wind speed at 10 m", "Shear below cloud",
                "Cloud shear", "Boundary layer dissipation",
                "Max. vertical velocity (up)",
                "Liquids updraft around -10 C", "Liquids around -10 C",
                "Cloud liquids -10 to -20 C", "Solids around -10 C", 
                "Cloud snow -10 to -20 C", "Cloud ice -10 to -20 C",
                "Cloud snow -20 to -40 C", "Cloud ice -20 to -40 C",
                "Supercooled liquids, total", "Snow, total", "Ice, total",
                "Ice divergence", "Cloud base height agl.", "Cloud thickness",
                "Convective prcp. 1h-sum", "Large scale prcp. 1h-sum",
                "Max. precipitation rate (hour)",
                "Vapor -10 to -20 C", "Moisture convergence",
                "Vapor, total", "Dew point at 2 m",
                "Surface sensible heat (up)", "Surface latent heat (up)",
                "Surface solar radiation (down)", "Boundary layer height"
)
renaming <- data.frame(cov, names_plot)


#' ### PCA (Figure 02)

#' The manuscript uses the long names as annotations. Here the short names are
#' used to reduce overlapping of the labels.
#' For long names uncomment:
# rownames(loadingsS) <- renaming[match(rownames(loadingsS), renaming[,1]) ,2]
# rownames(loadingsS_allYear) <- renaming[match(rownames(loadingsS_allYear),
#                                               renaming[,1]) ,2]

makePcaPlot <- function(pcaL, dat, loadingsS,
                        xlim = c(-8, 25), ylim = c(-10, 12),
                        filename,
                        fig_title = "PCA, Figure 02") {
  
  pdf(filename, width = 12, height = 6)
  par(mar = c(3, 3, 3, 3) + 2, xpd = FALSE)
  
  plot(pcaL$x[,1], pcaL$x[,2],
       type = "n",
       xlab = paste("PC 1 (", round(summary(pcaL)$importance[2,1] * 100,
                                    digits = 2), "%)"),
       ylab = paste("PC 2 (", round(summary(pcaL)$importance[2,2] * 100,
                                    digits = 2), "%)"),
       xlim = xlim,
       ylim = ylim,
       las = 1,
       main = fig_title
  )
  
  # Add biplot / loadings
  arrows(0, 0,
         x1 = loadingsS[,1],
         y1 = loadingsS[,2],
         length = .1,
         col = "gray40")
  
  # Add points
  points(pcaL$x[,1], pcaL$x[,2],
         col = adjustcolor(cols[dat$cluster_new], alpha.f = 0.1),
         pch = syms[dat$cluster_new],
         cex = 1.7,
         cex.lab = 1.7,
         cex.axis = 1.3,
         las = 1,
         lwd = 2
  )
  
  # Add labels to arrows
  text(loadingsS,
       labels = rownames(loadingsS),
       cex = .7,
       offset = .2,
       col = "black",
       xpd = TRUE)
  
  # Legend
  lgnd <- legend("topright",
                 legend = rep("", 5),
                 bty = 'n',
                 pch = syms,
                 col = cols,
                 cex = 1.2,
                 pt.cex = 1.2,
                 pt.lwd = 3)
  text(lgnd$rect$left,
       lgnd$text$y,
       legend_labels,
       pos = 2)
  
  # Helper lines
  abline(v = 0, col = "darkgray", lty = 2)
  abline(h = 0, col = "darkgray", lty = 2)
  
  invisible(dev.off()) 
}

#' Create PCA plot
makePcaPlot(pcaL = pcaL, dat = dat, loadingsS = loadingsS,
            filename = "morgenstern2021_fig02.pdf")
makePcaPlot(pcaL = pcaL_allYear, dat = dat_allYear,
            loadingsS = loadingsS_allYear,
            filename = "morgenstern2021_fig02_allYear.pdf",
            fig_title = "PCA, Figure 02, whole year")


#' ### Mosaic Plot (Figure 03)

makeMosaicPlot <- function(dat, filename, dat_type,
                           fig_title = "Mosaicplot, Figure 03") {

  # Calculate table
  tab <- table(dat$scenario, dat$cluster_new)
  
  # Order rows
  if (dat_type == "paper") {
    tab <- tab[c(3, 4, 2, 1), c(4, 5, 1, 2, 3)]
  } else if (dat_type == "allYear") {
    tab <- tab[c(7, 5, 3, 1, 8, 6, 4, 2), c(5, 4, 3, 2, 1)]
  } else {
    stop("Give either 'paper' or 'allYear' as dat_type.")
  }
  
  ord <- colnames(tab)
  
  # Remove names to prevent plotting
  colnames(tab) <- rep('', length(colnames(tab)))
  rownames(tab) <- rep('', length(rownames(tab)))
  
  
  #' ### Plot Figure 03
  
  pdf(filename, width = 8, height = 6)
  
  par(mar = c(8, 14, 2, 0), xpd = TRUE, lheight = 0.8, cex.main = 1.7)
  
  # Bars
  mosaicplot(tab,
             main = "",
             color = cols_mosaic[ord],
             cex = 1.4,
             off = c(15, 0),
             las = 2)
  
  # Title
  title(main = fig_title)
  
  # x labels
  if (dat_type == "allYear") {
    text(x = seq(0.12, 1, 0.12),
         y = -0.18,
         labels = c("Winter", "Spring", "Summer", "Fall",
                    "Winter", "Spring ", "Summer", "Fall"),
         srt = 45,
         cex = 1.8,
         pos = 2)
    text(x = c(0.48, 1),
         y = -0.1,
         labels = c("With lightning",
                    "Without lightning"),
         cex = 1.8,
         pos = 2)
    
  } else {
    text(x = c(0.2, 0.45, 0.7, 0.95),
         y = -0.1,
         labels = c("Lightning   \\nin winter", "No lightning   \\nin winter",
                    "No lightning   \\nin summer", "Lightning   \\nin summer"),
         srt = 45,
         cex = 1.8,
         pos = 2)
  }
  
  # Legend symbols
  lgnd <- legend(-.1, .96,
                 col = cols_mosaic[ord],
                 pch = 15,
                 cex = 2.5,
                 legend = rep("\\n", 5),
                 bty = "n")
  
  # Legend text
  text(lgnd$text$x - .15,
       lgnd$text$y,
       legend_labels_2lines[ord],
       pos = 2,
       cex = 1.8)

 invisible(dev.off())
}

#' Create mosaic plot
makeMosaicPlot(dat = dat, filename = "morgenstern2021_fig03.pdf",
               dat_type = "paper")
makeMosaicPlot(dat = dat_allYear, dat_type = "allYear",
               filename = "morgenstern2021_fig03.pdf_allYear",
               fig_title = "Mosaicplot, Figure 03, whole year")


#' ### Matplot (Figure 04)

makeMatplot <- function(dat, dat_scaled, filename,
                        fig_title = "Matplot, Figure 04") {
  
  # Prepare Plot
  # Calculate cluster means
  clu_means <- aggregate(dat_scaled, list(dat$cluster_new), mean)
  
  # Short names to long names
  colnames(clu_means) <- renaming[match(colnames(clu_means), renaming[,1]) ,2]
  colnames(clu_means)[1] <- "cluster"
  
  #' Plot
  pdf(filename, width = 12, height = 6)
  
  par(mar = c(11, 5, 7, 2), xpd = TRUE, cex.main = 1.7)
  
  # Empty plot
  matplot(t(clu_means[, 2:ncol(clu_means)]),
          type = "n",
          ylim = c(-2, 6),
          xaxt = "n",
          las = 1,
          cex.lab = 1.3,
          cex.axis = 1.3,
          ylab = "Scaled value"
  )
  
  # Title
  title(fig_title, line = 5)
  
  # Gray vertical lines
  xtick <- seq(1, ncol(clu_means) - 1)
  abline(v = xtick, col = "gray", xpd = FALSE)
  
  # xlab
  text(y = -2.8,
       x = xtick,
       labels = names(clu_means)[2:ncol(clu_means)],
       las = 1,
       srt = 45,
       adj = 1,
       cex = 1.2)
  
  # Legend
  legend("top",
         inset = c(0, -.6),
         ncol = 5,
         legend = legend_labels_2lines,
         col = cols,
         pch = syms,
         cex = 1.5,
         pt.lwd = 3,
         bty = "n")
  
  # Horizontal 0 line
  abline(h = 0, lwd = 2, xpd = FALSE)
  
  # Annotations of physical categories
  xpos <- c(par('usr')[1], 4.5, 10.5, 27.5, 31.5, par('usr')[2])
  labs <- c("Mass field", "Wind field", "Cloud physics",
            "Moisture\\nfield", "Surface\\nexchange")
  text(xpos[-1] - diff(xpos)/2, par('usr')[4] - .25, labs,
       pos = 1, cex = 1.5, font = 3)
  abline(v = xpos, lwd = 2, xpd = FALSE)  
  
  # Lines
  matplot(t(clu_means[, 2:ncol(clu_means)]),
          lwd = 2,
          col = cols[clu_means$cluster],
          type = "o",
          pch = syms[clu_means$cluster],
          cex = 1.3,
          lty = 1,
          ylim = c(-2, 6),
          xaxt = "n",
          cex.lab = 1.3,
          cex.axis = 1.3,
          ylab = "Scaled value",
          add = TRUE
  )
  
  invisible(dev.off())
}

#' Create matplot
makeMatplot(dat = dat, dat_scaled = dat_scaled,
            filename = "morgenstern2021_fig04.pdf")
makeMatplot(dat = dat_allYear, dat_scaled = dat_scaled_allYear,
            filename = "morgenstern2021_fig04_allYear.pdf",
            fig_title = "Matplot, Figure 04, whole year")

#' ### Session information

#' In case you have problems with this script, compare your sessionInfo with
#' the one in which this script was written:
#' 
#' R version 4.1.1 (2021-08-10)
#' 
#' Platform: x86_64-pc-linux-gnu (64-bit)
#' 
#' Running under: Debian GNU/Linux 10 (buster)
#' 
#' 
#' Matrix products: default
#' 
#' BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
#' 
#' LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
#' 
#' 
#' locale:
#'  LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8     
#'  LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#'  LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C          
#'  LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C  
#'      
#' 
#' attached base packages:
#' stats     graphics  grDevices utils     datasets  methods   base     
#' 
#' loaded via a namespace (and not attached):
#' compiler_4.1.1 magrittr_2.0.1 markdown_1.1   tools_4.1.1    stringi_1.7.4
#' highr_0.9      knitr_1.34    stringr_1.4.0  xfun_0.25      evaluate_0.14
#'

# Get your sessionInfo:
# sessionInfo()


