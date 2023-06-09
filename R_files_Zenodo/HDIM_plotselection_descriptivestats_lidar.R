### make histograms of site level CAO LiDAR data
## Feb 2014 DSG, revised for new polygons April 2014
## August 2022 getting some site data for Natalie's paper

library(sp)
library(rgdal)
library(raster)
library(maptools)
library(tidyverse)
library(ggplot2)
library(readxl)
library(here)

setwd("G:/My Drive/GrunerGIS/DimensionsHI/plot_selection3/cao")

RAST.list <- list(
  "HAVO_Olaa_Old_CAOMN.tif",
  "HAVO_Thurston_2polygons_CAOMN.tif",
  "HiloFR_Humuula_2polygons_CAOMN.tif",
  "HiloFR_Laupahoehoe_Old_CAOMN.tif",
  "NAR_Laupahoehoe_65-250K_CAOMN.tif",
  "NAR_Laupahoehoe_HIPPNET_Young_CAOMN.tif",
  "NAR_PuuOumi_Old_4polygons_CAOMN.tif",
  "NAR_PuuOumi_Young02_CAOMN.tif",
  "HAVO_EscapeRd_1973_CAOMN.tif",
  "HAVO_EscapeRd_HighStature_CAOMN.tif",
  "KauFR_Alili_2polygons_CAOMN.tif",
  "NAR_PuuMakaala_Old_2polygons_CAOMN.tif",
  "NAR_PuuMakaala_Young_2polygons_CAOMN.tif",
  "TNC_Kaiholena_Old_CAOMN.tif",
  "TNC_Kaiholena_Young_CAOMN.tif",
  "Kauai_8m_CAOMN.tif"
)


## output histograms for each raster
for (i in 1:length(RAST.list)) { 
  RAST.i <- raster(RAST.list[[i]])
  File.Name <- paste(names(RAST.i), "AUG22.jpg", sep = "")
  jpeg(file = File.Name)
  raster::hist(RAST.i,breaks=20)
  dev.off()
}

## to obtain values from each cell in each raster
for (i in 1:length(RAST.list)) { 
  RAST.i <- raster(RAST.list[[i]])
  vec.i <- raster::getValues(RAST.i)
  vec.i<-vec.i[!is.na(vec.i)]
  vec.i[complete.cases(vec.i)]  # probably redundant with above
  write.csv(vec.i, file=paste0(names(RAST.i), ".csv"), row.names = FALSE)
}

# all csv files have leading "x" as first row. Don't know why

CSV <- read.csv("NAR_PuuMakaala_Old_2polygons_CAOMN.csv")
CSV <- CSV[-1,]
summary(CSV)
hist(CSV)

# could write a loop to remove the first column, but also need to extract column names
# easier to paste in excel to create "SiteValues.xlsx"
# here is the list of files, if proceeding with the for loop approach
CSV.list <- list(
  "HAVO_Olaa_Old_CAOMN.csv",
  "HAVO_Thurston_2polygons_CAOMN.csv",
  "HiloFR_Humuula_2polygons_CAOMN.csv",
  "HiloFR_Laupahoehoe_Old_CAOMN.csv",
  "NAR_Laupahoehoe_65-250K_CAOMN.csv",
  "NAR_Laupahoehoe_HIPPNET_Young_CAOMN.csv",
  "NAR_PuuOumi_Old_4polygons_CAOMN.csv",
  "NAR_PuuOumi_Young02_CAOMN.csv",
  "HAVO_EscapeRd_1973_CAOMN.csv",
  "HAVO_EscapeRd_HighStature_CAOMN.csv",
  "KauFR_Alili_2polygons_CAOMN.csv",
  "NAR_PuuMakaala_Old_2polygons_CAOMN.csv",
  "NAR_PuuMakaala_Young_2polygons_CAOMN.csv",
  "TNC_Kaiholena_Old_CAOMN.csv",
  "TNC_Kaiholena_Young_CAOMN.csv",
  "Kauai_8m_CAOMN.csv"
)


# read in spreadsheet that contains all the CAOMN pixel values by site
df <- read_xlsx("SiteValues.xlsx", sheet = 1, col_names = TRUE, na = "")

#TMP1 <- data.frame(summary(df))

# calculate descriptive stats in for loops
MEDIAN <- vector("double", ncol(df)) 
for (i in seq_along(df)) {  
  MEDIAN[[i]] <- median(df[[i]],na.rm=TRUE)  
}
MEDIAN

MEAN <- vector("double", ncol(df))
for (i in seq_along(df)) {
  MEAN[[i]] <- mean(df[[i]],na.rm=TRUE)
}
MEAN

MIN <- vector("double", ncol(df))
for (i in seq_along(df)) {
  MIN[[i]] <- min(df[[i]],na.rm=TRUE)
}
MIN

MAX <- vector("double", ncol(df))
for (i in seq_along(df)) {
  MAX[[i]] <- max(df[[i]],na.rm=TRUE)
}
MAX

SD <- vector("double", ncol(df))
for (i in seq_along(df)) {
  SD[[i]] <- sd(df[[i]],na.rm=TRUE)
}
SD

#quantile(df[[16]], probs = c(0,0.25,0.40,0.6,0.75,1.0), na.rm = TRUE)
Q40 <- vector("double", ncol(df))
for (i in seq_along(df)) {
  Q40[[i]] <- quantile(df[[i]], probs = 0.6, na.rm = TRUE)
}
Q40

# get states for kurtosis and skewness in datawizard (see documentation)
require(datawizard)
SKEW  <- datawizard::skewness(df,na.rm=TRUE)
KURT  <- datawizard::kurtosis(df,na.rm=TRUE)

# bind columns for simple stats. Kurtosis and Skewness output separately
CBOUND <- data.frame(
  cbind(MEAN, MEDIAN, MIN, MAX, SD, Q40),
  row.names = c(
    "HAVO_Olaa_Old",
    "HAVO_Thurston",
    "HiloFR_Humuula",
    "HiloFR_Laupahoehoe_Old",
    "NAR_Laupahoehoe_65-250K",
    "NAR_Laupahoehoe_HIPPNET",
    "NAR_PuuOumi_Old",
    "NAR_PuuOumi_Young02",
    "HAVO_EscapeRd_1973",
    "HAVO_EscapeRd_HighStature",
    "KauFR_Alili",
    "NAR_PuuMakaala_Old",
    "NAR_PuuMakaala_Young",
    "TNC_Kaiholena_Old",
    "TNC_Kaiholena_Young",
    "Kauai_8m"
  )
)


write_excel_csv(CBOUND, file="summary.csv")
write_excel_csv(KURT, file="kurtosis.csv")
write_excel_csv(SKEW, file="skewness.csv")

