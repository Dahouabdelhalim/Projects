### worldclim_analyses_clean.R
### Libby Natola
### Started 17 January 2022, cleaned up 22 March 2022
### Plot the transitions of the 4 most closely associated worldclim traits (according to Natola and Burg 2018) 

### set working directory
setwd("~/Documents/UBC/Bioinformatics/rnrb/worldclim")

#library packages
# !!!NOTE!!! tidyr has a function called extract which bothers raster. Unload it first with 
#.rs.unloadPackage("tidyr")
library(raster)
library(sp)
library(geosphere)
library(ggplot2)

### get the worldclim data, use res 10 (broadest, low res), 5, 2.5, 0.5(smallest, high res) 
r <- getData("worldclim",var="bio",res=10)
r_elev <- getData("worldclim", var="alt", res=10)

### S. nuchalis with isothermality [BIO3] (58.9%) and mean annual temperature [BIO1] (23%) S. ruber with precipitation in the coldest quarter [BIO19] (40.9%) and temperature of the wettest quarter [BIO8] (25.4%)

#removed this for Shawn's data analysis
#r <- r[[c(1,3,8,19)]]
#names(r) <- c("Temp","Isothermality","TempWettestQuarter","PrecipColdestQuarter")

### load transect data
nbc_info <- as.data.frame(read.table("../sample_info_nbc.txt", header=TRUE))
man_info <- as.data.frame(read.table("../sample_info_man.txt", header=TRUE))
caor_info <- as.data.frame(read.table("../sample_info_caor.txt", header=TRUE))


########## extract data along the latitudinal transect crossing the isocline midpoint #########
# make a df for line across plot at 0.5 isocline midpoint (found in isocline_clean.R)
# nbc midpoint # -126.2092, 52.41907
# caor midpoint # -120.0276, 42.88465,
# man midpoint # -120.9136, 49.25926

caor_midpoint_longs <- c(-123.56, -117.27)
caor_midpoint_lat <- c(42.88465, 42.88465)
caor_line <- data.frame(caor_midpoint_longs, caor_midpoint_lat)
caor_trans <- Line(caor_line)


nbc_midpoint_longs <- c(-127, -119)
nbc_midpoint_lat <- c(52.41907, 52.41907)
nbc_line <- data.frame(nbc_midpoint_longs, nbc_midpoint_lat)
nbc_trans <- Line(nbc_line)


man_midpoint_longs <- c(-122.1, -119.1)
man_midpoint_lat <- c(49.25926, 49.25926)
man_line <- data.frame(man_midpoint_longs, man_midpoint_lat)
man_trans <- Line(man_line)


### id a series of equidistant points, one per worldclim grid, along the spatial line
# used the 30 seconds (~1 km2) grid  so I want one point approx every 1 km. So measure length of line in km, then divide by 1 km to decide how many points to sample

# caor
(distCosine(c(-123.56, 42.88465), c(-117.27, 42.88465)))/1000
# [1] 512.9 points lets do 513
caor_trans_pts <- spsample(caor_trans, 513, type = "regular")
plot(caor_trans_pts)
caor_trans_pts <- caor_trans_pts@coords
colnames(caor_trans_pts) <- c("Long", "Lat")
write.csv(caor_trans_pts, "caor_trans_pts.csv")

# nbc -126.2092, 52.41907 c(-127, -119)
(distCosine(c(-127, 52.41907), c(-119, 52.41907)))/1000
# [1] 542.9 points lets do 543
nbc_trans_pts <- spsample(nbc_trans, 543, type = "regular")
plot(nbc_trans_pts)
nbc_trans_pts <- nbc_trans_pts@coords
colnames(nbc_trans_pts) <- c("Long", "Lat")
write.csv(nbc_trans_pts, "nbc_trans_pts.csv")


# man -120.9136, 49.25926 c(-122.1, -119.1)
(distCosine(c(-122.1, 49.25926), c(-119.1, 49.25926)))/1000
# [1] 217.9 points lets do 218
man_trans_pts <- spsample(man_trans, 218, type = "regular")
plot(man_trans_pts)
man_trans_pts <- man_trans_pts@coords
colnames(man_trans_pts) <- c("Long", "Lat")
write.csv(man_trans_pts, "man_trans_pts.csv")


### plot the values along those lines
# get the worldclim data
r_elev <- getData("worldclim", var="alt", res=2.5)

#random forest layers
rf <- r[[c(2,4,9,12,13,15)]]
names(rf) <- c("MeanDiurnalRange", "TemperatureSeasonality", "MeanTemperatureDriestQuarter", "AnnualPrecipitation", "PrecipitationWettestMonth", "PrecipitationSeasonality")

names(r) <- c("Annual_Mean_Temperature", "Mean_Diurnal_Range", "Isothermality", "Temperature_Seasonality", "Max_Temperature_of_Warmest_Month", "Min_Temperature_of_Coldest_Month", "Temperature_Annual_Range", "Mean_Temperature_of_Wettest_Quarter", "Mean_Temperature_of_Driest_Quarter", "Mean_Temperature_of_Warmest_Quarter", "Mean_Temperature_of_Coldest_Quarter", "Annual_Precipitation", "Precipitation_of_Wettest_Month", "Precipitation_of_Driest_Month", "Precipitation_Seasonality", "Precipitation_of_Wettest_Quarter", "Precipitation_of_Driest_Quarter", "Precipitation_of_Warmest_Quarter", "Precipitation_of_Coldest_Quarter")

# make an object for the values along the lines
#caor
caor_trans_pts_values <- extract(r,caor_trans_pts)
caor_trans_pts_elev_values <- extract(r_elev,caor_trans_pts)
caor_trans_df <- cbind.data.frame(coordinates(caor_trans_pts),caor_trans_pts_values,caor_trans_pts_elev_values)
ggplot(caor_trans_df, aes(x=Long, y=Precipitation_Seasonality)) + geom_point()

#nbc
nbc_trans_pts_values <- extract(r,nbc_trans_pts)
nbc_trans_pts_elev_values <- extract(r_elev,nbc_trans_pts)
nbc_trans_df <- cbind.data.frame(coordinates(nbc_trans_pts),nbc_trans_pts_values,nbc_trans_pts_elev_values)
ggplot(nbc_trans_df, aes(x=Long, y=Temperature_Seasonality)) + geom_point()

#man
man_trans_pts_values <- extract(r,man_trans_pts)
man_trans_pts_elev_values <- extract(r_elev,man_trans_pts)
man_trans_df <- cbind.data.frame(coordinates(man_trans_pts),man_trans_pts_values,man_trans_pts_elev_values)
ggplot(man_trans_df, aes(x=Long, y=Precipitation_of_Wettest_Month)) + geom_point()

# convert the x axis into the same x-axis as the clines, so distance from the midpoint on the 0.5 isocline
#caor
# # -120.0276, 42.88465,
caor_trans_dist <- (distCosine(caor_trans_pts, c(-120.0276, 42.88465)))/1000
caor_trans_negs <- (caor_trans_pts[,1]-(-120.0276))/(abs(caor_trans_pts[,1]-(-120.0276)))
caor_trans_dist <- caor_trans_dist * caor_trans_negs

#nbc
# # -126.2092, 52.41907
nbc_trans_dist <- (distCosine(nbc_trans_pts, c(-126.2092, 52.41907)))/1000
nbc_trans_negs <- (nbc_trans_pts[,1]-(-126.2092))/(abs(nbc_trans_pts[,1]-(-126.2092)))
nbc_trans_dist <- nbc_trans_dist * nbc_trans_negs

#man
# # -120.9136, 49.25926,
man_trans_dist <- (distCosine(man_trans_pts, c(-120.9136, 49.25926)))/1000
man_trans_negs <- (man_trans_pts[,1]-(-120.9136))/(abs(man_trans_pts[,1]-(-120.9136)))
man_trans_dist <- man_trans_dist * man_trans_negs

#get isocline data
caor_loess_envirocline <- read.csv("../isoclines/caor_isocline_envirocline_ancestry.csv")
colnames(caor_loess_envirocline) <- c("X", "iso_ancestry")
nbc_loess_envirocline <- read.csv("../isoclines/nbc_isocline_envirocline_ancestry.csv")
colnames(nbc_loess_envirocline) <- c("X", "iso_ancestry")
man_loess_envirocline <- read.csv("../isoclines/man_isocline_envirocline_ancestry.csv")
colnames(man_loess_envirocline) <- c("X", "iso_ancestry")


# combine the data for the distance and the worldclim data per point and iso_ancestry in a df
#caor
caor_trans_worldclim <- cbind.data.frame(coordinates(caor_trans_pts),caor_trans_pts_values,caor_trans_pts_elev_values, caor_trans_dist, caor_loess_envirocline$iso_ancestry)
write.csv(caor_trans_worldclim, "caor_trans_worldclim.csv")

#nbc
nbc_trans_worldclim <- cbind.data.frame(coordinates(nbc_trans_pts),nbc_trans_pts_values,nbc_trans_pts_elev_values, nbc_trans_dist, nbc_loess_envirocline$iso_ancestry)
write.csv(nbc_trans_worldclim, "nbc_trans_worldclim.csv")

#man
man_trans_worldclim <- cbind.data.frame(coordinates(man_trans_pts),man_trans_pts_values,man_trans_pts_elev_values, man_trans_dist, man_loess_envirocline$iso_ancestry)
write.csv(man_trans_worldclim, "man_trans_worldclim.csv")

