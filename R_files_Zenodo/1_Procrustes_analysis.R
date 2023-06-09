
# --------------------------------------------------------------------------------
# Install and load packages

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

packages <- c("geomorph", "RRPP", "ggplot2")
ipak(packages)

# global options
options(scipen = 20)


# --------------------------------------------------------------------------------
# read morphologika data into R

setwd("path-to-directory") # edit file path for working directory

Fox3D <- read.morphologika("Fox_data_Morphologika.txt")

# --------------------------------------------------------------------------------
# GPA analysis in geomorph

FoxGPA <- gpagen(A = Fox3D$coords)

# make geomorph data frame with sex, population, and centroid size
FoxDF <- geomorph.data.frame(
    sex = factor(Fox3D$labels[, 1]), # labels from morphologika file
    population = factor(Fox3D$labels[, 2]), # labels from morphologika file
    coords = FoxGPA$coords, 
    csize = FoxGPA$Csize
    )


# --------------------------------------------------------------------------------
# Manova on Procrustes Coordinates (i.e., Procrustes Anova)

# shape differences between populations and sexes
shape_diff <- procD.lm(coords ~ sex + population, iter = 99999, SS.type = "III", RRPP = TRUE, data = FoxDF)
summary(shape_diff)

# Pairwise comparison between populations
pairwise_shape_diff <- pairwise(shape_diff, groups = FoxDF$population)

# Pairwise comparisons of mean differences
mean_shape_diff <- summary(pairwise_shape_diff, test.type = "dist") # pairwise differences in means 
mean_shape_diff$summary.table$p.adj <- p.adjust(mean_shape_diff$summary.table$`Pr > d`, method = "holm") # familywise error
mean_shape_diff

# Pairwise comparisons of variance differences
var_shape_diff <- summary(pairwise_shape_diff, test.type = "var") # pairwise differences in variances 
var_shape_diff$summary.table$p.adj <- p.adjust(var_shape_diff$summary.table$`Pr > d`, method = "holm") # familywise error
var_shape_diff


# --------------------------------------------------------------------------------
# Allometry

# Common allometry (parallel allometric relationships for each population, but different means)
common_allometry <- procD.lm(coords ~ csize + population + sex, data = FoxDF, iter = 99999)
summary(common_allometry)

# Unique allometry (different allometric relationships for each population)
unique_allometry <- procD.lm(coords ~ csize * population + sex, data = FoxDF, iter = 99999)
summary(unique_allometry)

# model comparison
anova(common_allometry, unique_allometry)

# calculate shape scores
common_allometry_data <- plotAllometry(common_allometry, size = FoxDF$csize, logsz = FALSE, method = "PredLine")
unique_allometry_data <- plotAllometry(unique_allometry, size = FoxDF$csize, logsz = FALSE, method = "PredLine")

# combine shape score and size data
allometry <- data.frame(
    specimen = rep(names(common_allometry_data$size.var), times = 4),
    allometry_type = c(rep(c("Common Allometry", "Unique Allometry"), each = 73), 
                       rep(c("Common Allometry", "Unique Allometry"), each = 73)),
    population = rep(FoxDF$population, times = 4), 
    shape_metric = rep(c("PredLine", "RegScore"), each = 73*2),
    shape_score = c(common_allometry_data$PredLine, unique_allometry_data$PredLine,
                    common_allometry_data$RegScore, unique_allometry_data$RegScore),
    csize = rep(common_allometry_data$size.var, times = 4)
)       


# save objects
save(FoxDF, 
     shape_diff, 
     mean_shape_diff, 
     var_shape_diff, 
     allometry, 
     file = "shape_objects.Rdata", 
     compress = "gzip")



