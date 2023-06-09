### isocline_clean.R
### Libby Natola
### Started 29 Jan 2022 then cleaned up 20 March 2022
### Plot isoclines for 0.5 ancestry through transect, then measure distance of points to isocline to make x-axis for cline analysis

library(ggplot2)
library(sf)
library(stringr)
library(raster)
library(sp)
library(rgeos)
library(geosphere)
library(splines)
library(dplyr)
library(rebus)
library(fANCOVA)
library(maptools)

setwd("~/Documents/UBC/Bioinformatics/rnrb/isoclines/")

########## caor ##########

# read in metadata as data frame
caor_info <- read.table("../sample_info_caor.txt", header=TRUE)
as.data.frame(caor_info)

# read in admixture k=2 output 
caor_admix_values <- read.table("../admixture/caor_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.2.Q")

# read in admixture input to label sample IDs
caor_fam <- read.table("../admixture/caor_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.fam")
as.data.frame(caor_fam)

# add sample IDs to 2 clusters
caor_admix_data <- cbind(caor_fam$V2, caor_admix_values)

#label columns
colnames(caor_admix_data) <- c("ID", "p", "q")

# join with rest of metadata
caor_admix_data <- left_join(caor_admix_data, caor_info)

# make loess, fiddle with span until it looks right
caor_topo_loess <- loess(q ~ long*lat, caor_admix_data, span = 0.9, degree = 2, normalize = F)

# specify the spatial extent of the grid
caor_grid_lat_long <- list(long = seq(-123.56, -117.27, 0.01), lat = seq(40.41, 45.61, 0.01))
caor_grid <- expand.grid(caor_grid_lat_long)
caor_loess_grid <- predict(caor_topo_loess, caor_grid, se = FALSE)

#plot the loess isoclines and points on the grid
pdf("caor_isoclines_contour.pdf")
contour(caor_grid_lat_long$long, caor_grid_lat_long$lat, caor_loess_grid)
points(caor_admix_data$long, caor_admix_data$lat)
dev.off()

caor_isoclines <- contourLines(caor_grid_lat_long$long, caor_grid_lat_long$lat, caor_loess_grid)
caor_isoclines   # This shows all the things stored in this object. By looking through, I see the 3rd item (in this case only) is the 0.5 isocline. So extract that: the line 1 we want (circly guy) is [[3]] line 0 is 
caor_fifty_percent_isocline <- caor_isoclines[[3]]
caor_x_isocline <- as.numeric(unlist(caor_fifty_percent_isocline[2]))
caor_y_isocline <- as.numeric(unlist(caor_fifty_percent_isocline[3]))

#save the isocline in a file so I can grab it for other analyses
caor_fifty_percent_isocline_df <- cbind(caor_x_isocline, caor_y_isocline)
colnames(caor_fifty_percent_isocline_df) <- c("long", "lat")
write.csv(caor_fifty_percent_isocline_df, "caor_0.5_isocline.csv", row.names = FALSE)

### add 0.5 isocline as a line in sp, 
caor_isocline <- spLines(cbind(caor_x_isocline,caor_y_isocline))

# make df of points
caor_points <- SpatialPoints(coords = cbind(caor_admix_data$long, caor_admix_data$lat))

# use geosphere() to calculate distance of caor points to line with dist2Line
caor_dist <- as.data.frame(dist2Line(caor_points, caor_isocline, distfun=distGeo))
caor_dist$distance <- (caor_dist$distance)/1000

# find out if each point is west (longitude < longitude of closest spot on line) of the isocline and therefore a negative value or east of the line and therefore a positive value
caor_longitudes <- cbind(caor_admix_data$long,caor_dist$lon)
caor_negs <- caor_longitudes[,1]-caor_longitudes[,2]
unlist(caor_negs)
caor_dist$distance <- caor_dist$distance * (caor_negs/abs(caor_negs))
write.csv(caor_dist, "caor_dist.csv")

# add to the rest of the metadata
caor_admix_data$dist <- caor_dist$distance

ggplot(caor_admix_data, aes(x=dist, y=q, color=spp)) + geom_point()

# get the midpoint of the isocline line
caor_midpoint <- getSpatialLinesMidPoints(caor_isocline)
# -120.0276, 42.88465,

# make a df for line across plot at 0.5 isocline midpoint
caor_midpoint_longs <- c(-123.56, -117.27)
caor_midpoint_lat <- c(42.88465, 42.88465)
caor_line <- data.frame(caor_midpoint_longs, caor_midpoint_lat)

#plot them together
pdf("caor_isoclines_and_transect.pdf")
contour(caor_grid_lat_long$long, caor_grid_lat_long$lat, caor_loess_grid)
points(caor_admix_data$long, caor_admix_data$lat, col=rgb(33/255, 145/255, 140/255))
lines(caor_line, col=rgb(33/255, 145/255, 140/255))
dev.off()


########## nbc ##########

# read in metadata as data frame
nbc_info <- read.table("../sample_info_nbc.txt", header=TRUE)
as.data.frame(nbc_info)

# read in admixture k=2 output 
nbc_admix_values <- read.table("../admixture/nbc_max-miss30_filtered.snps_unrelated_filtered.birds_noeasties_refiltered.pruned.2.Q")

# read in admixture input to label sample IDs
nbc_fam <- read.table("../admixture/nbc_max-miss30_filtered.snps_unrelated_filtered.birds_noeasties_refiltered.pruned.fam")
as.data.frame(nbc_fam)

# add sample IDs to 2 clusters
nbc_admix_data <- cbind(nbc_fam$V2, nbc_admix_values)

#label columns
colnames(nbc_admix_data) <- c("ID", "p", "q")

# join with rest of metadata
nbc_admix_data <- left_join(nbc_admix_data, nbc_info)

# make loess, fiddle with span until it looks right
nbc_topo_loess <- loess(q ~ long*lat, nbc_admix_data, span = 0.8, degree = 2, normalize = F)
max(nbc_admix_data$lat)
min(nbc_admix_data$lat)
max(nbc_admix_data$long)
min(nbc_admix_data$long)
# specify the spatial extent of the grid
nbc_grid_lat_long <- list(long = seq(-127, -119, 0.01), lat = seq(50.5, 53, 0.01))
nbc_grid <- expand.grid(nbc_grid_lat_long)
nbc_loess_grid <- predict(nbc_topo_loess, nbc_grid, se = FALSE)

#plot the loess isoclines and points on the grid
pdf("nbc_isoclines_contour.pdf")
contour(nbc_grid_lat_long$long, nbc_grid_lat_long$lat, nbc_loess_grid)
points(nbc_admix_data$long, nbc_admix_data$lat)
dev.off()

nbc_isoclines <- contourLines(nbc_grid_lat_long$long, nbc_grid_lat_long$lat, nbc_loess_grid)
nbc_isoclines   # This shows all the things stored in this object. By looking through, I see the 14th  item (in this case only) is the 0.5 isocline. So extract that:
nbc_fifty_percent_isocline <- nbc_isoclines[[14]]
# to see what is in this:
nbc_x_isocline <- as.numeric(unlist(nbc_fifty_percent_isocline[2]))
nbc_y_isocline <- as.numeric(unlist(nbc_fifty_percent_isocline[3]))

#save the isocline in a file so I can grab it for other analyses
nbc_fifty_percent_isocline_df <- cbind(nbc_x_isocline, nbc_y_isocline)
colnames(nbc_fifty_percent_isocline_df) <- c("long", "lat")
write.csv(nbc_fifty_percent_isocline_df, "nbc_0.5_isocline.csv", row.names = FALSE)

### add 0.5 isocline as a line in sp, 
# get the values of the loess on the grid
nbc_isocline <- spLines(cbind(nbc_x_isocline,nbc_y_isocline))

# make df of points
nbc_points <- SpatialPoints(coords = cbind(nbc_admix_data$long, nbc_admix_data$lat))

# use geosphere() to calculate distance of caor points to line with dist2Line, it reports the distance in m so I divided by 1000 to make km
nbc_dist <- as.data.frame(dist2Line(nbc_points, nbc_isocline, distfun=distGeo))
nbc_dist$distance <- (nbc_dist$distance)/1000

# find out if each point is west (longitude < longitude of closest spot on line) of the isocline and therefore a negative value or east of the line and therefore a positive value
nbc_longitudes <- cbind(nbc_admix_data$long,nbc_dist$lon)
nbc_negs <- nbc_longitudes[,1]-nbc_longitudes[,2]
unlist(nbc_negs)
nbc_dist$distance <- nbc_dist$distance * (nbc_negs/abs(nbc_negs))
write.csv(nbc_dist, "nbc_dist.csv")

# add to the rest of the metadata
nbc_admix_data$dist <- nbc_dist$distance

ggplot(nbc_admix_data, aes(x=dist, y=q, color=spp)) + geom_point()

# get the midpoint of the isocline line
nbc_midpoint <- getSpatialLinesMidPoints(nbc_isocline)
# -126.2092, 52.41907

# make a df for line across plot at 0.5 isocline midpoint
nbc_midpoint_longs <- c(-127, -119)
nbc_midpoint_lat <- c(52.41907, 52.41907)
nbc_line <- data.frame(nbc_midpoint_longs, nbc_midpoint_lat)

#plot them together
pdf("nbc_isoclines_and_transect.pdf")
contour(nbc_grid_lat_long$long, nbc_grid_lat_long$lat, nbc_loess_grid)
points(nbc_admix_data$long, nbc_admix_data$lat, col=rgb(253/255, 231/255, 37/255))
lines(nbc_line, col=rgb(253/255, 231/255, 37/255))
dev.off()


########## man ##########

# read in metadata as data frame
man_info <- read.table("../sample_info_man.txt", header=TRUE)
as.data.frame(man_info)

# read in admixture k=2 output 
man_admix_values <- read.table("../admixture/man_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.2.Q")

# read in admixture input to label sample IDs
man_fam <- read.table("../admixture/man_max-miss30_filtered.snps_unrelated_filtered.birds.pruned.fam")
as.data.frame(man_fam)

# add sample IDs to 2 clusters
man_admix_data <- cbind(man_fam$V2, man_admix_values)

#label columns
colnames(man_admix_data) <- c("ID", "q", "p")

# join with rest of metadata
man_admix_data <- left_join(man_admix_data, man_info)

# make loess, fiddle with span until it looks right
man_topo_loess <- loess(q ~ long*lat, man_admix_data, span = 0.8, degree = 2, normalize = F)
max(man_admix_data$lat)
min(man_admix_data$lat)
max(man_admix_data$long)
min(man_admix_data$long)
# specify the spatial extent of the grid
man_grid_lat_long <- list(long = seq(-122.1, -119.1, 0.05), lat = seq(49, 49.5, 0.05))
man_grid <- expand.grid(man_grid_lat_long)
man_loess_grid <- predict(man_topo_loess, man_grid, se = FALSE)

#plot the loess isoclines and points on the grid
pdf("man_isoclines_contour.pdf")
contour(man_grid_lat_long$long, man_grid_lat_long$lat, man_loess_grid, nlevels=15)
points(man_admix_data$long, man_admix_data$lat)
dev.off()

man_isoclines <- contourLines(man_grid_lat_long$long, man_grid_lat_long$lat, man_loess_grid, nlevels=15)
man_isoclines   # This shows all the things stored in this object. By looking through, I see the 11th item (in this case only) is the 0.5 isocline. So extract that:
man_fifty_percent_isocline <- man_isoclines[[11]]
# to see what is in this:
man_x_isocline <- as.numeric(unlist(man_fifty_percent_isocline[2]))
man_y_isocline <- as.numeric(unlist(man_fifty_percent_isocline[3]))

#save the isocline in a file so I can grab it for other analyses
man_fifty_percent_isocline_df <- cbind(man_x_isocline, man_y_isocline)
colnames(man_fifty_percent_isocline_df) <- c("long", "lat")
write.csv(man_fifty_percent_isocline_df, "man_0.5_isocline.csv", row.names = FALSE)

### add 0.5 isocline as a line in sp, 
# get the values of the loess on the grid
man_isocline <- spLines(cbind(man_x_isocline,man_y_isocline))

# make df of points
man_points <- SpatialPoints(coords = cbind(man_admix_data$long, man_admix_data$lat))

# use geosphere() to calculate distance of caor points to line with dist2Line, it reports the distance in m so I divided by 1000 to make km
man_dist <- as.data.frame(dist2Line(man_points, man_isocline, distfun=distGeo))
man_dist$distance <- (man_dist$distance)/1000

# find out if each point is west (longitude < longitude of closest spot on line) of the isocline and therefore a negative value or east of the line and therefore a positive value
man_longitudes <- cbind(man_admix_data$long,man_dist$lon)
man_negs <- man_longitudes[,1]-man_longitudes[,2]
unlist(man_negs)
man_dist$distance <- man_dist$distance * (man_negs/abs(man_negs))
write.csv(man_dist, "man_dist.csv")


# add to the rest of the metadata
man_admix_data$dist <- man_dist$distance

ggplot(man_admix_data, aes(x=dist, y=q, color=spp)) + geom_point()

# get the midpoint of the isocline line
man_midpoint <- getSpatialLinesMidPoints(man_isocline)
# -120.9136, 49.25926

# make a df for line across plot at 0.5 isocline midpoint
man_midpoint_longs <- c(-122.1, -119.1)
man_midpoint_lat <- c(49.25926, 49.25926)
man_line <- data.frame(man_midpoint_longs, man_midpoint_lat)

#plot them together
pdf("man_isoclines_and_transect.pdf")
contour(man_grid_lat_long$long, man_grid_lat_long$lat, man_loess_grid)
points(man_admix_data$long, man_admix_data$lat, col=rgb(68/255, 1/255, 84/255))
lines(man_line, col=rgb(68/255, 1/255, 84/255))
dev.off()



#### get isocline ancestry estimate for each transect 
#read in caor envirocline
caor_enviro_trans <- read.csv("../worldclim/caor_trans_pts.csv")
colnames(caor_enviro_trans) <- c("X", "long", "lat")
caor_loess_envirocline <- as.data.frame(predict(caor_topo_loess, caor_enviro_trans, se = FALSE))
write.csv(caor_loess_envirocline, "caor_isocline_envirocline_ancestry.csv")

#read in caor envirocline
nbc_enviro_trans <- read.csv("../worldclim/nbc_trans_pts.csv")
colnames(nbc_enviro_trans) <- c("X", "long", "lat")
nbc_loess_envirocline <- as.data.frame(predict(nbc_topo_loess, nbc_enviro_trans, se = FALSE))
write.csv(nbc_loess_envirocline, "nbc_isocline_envirocline_ancestry.csv")

#read in caor envirocline
man_enviro_trans <- read.csv("../worldclim/man_trans_pts.csv")
colnames(man_enviro_trans) <- c("X", "long", "lat")
man_loess_envirocline <- as.data.frame(predict(man_topo_loess, man_enviro_trans, se = FALSE))
write.csv(man_loess_envirocline, "man_isocline_envirocline_ancestry.csv")
