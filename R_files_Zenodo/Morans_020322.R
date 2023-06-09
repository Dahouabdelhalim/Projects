# -------------------------------------------------------------------------
# Analyzing eBird data mean arrival dates and data quantity/quality
# -------------------------------------------------------------------------
# Load Libraries
library(dplyr)
library(readr)
library(stringr)
library(lme4)
library(spdep)
library(ape)

# Set Working Directory
setwd("~/Data repository")

# ---------------------------------
# Read in data
# ---------------------------------
eBird_Filtered <- read_csv("Routine 7/Merged output/MAD_Output_All_Filter_1.csv")

# ----------------------------------------
# MAPPPING CODE
# ----------------------------------------
# Extract lat long from the grids
eBird_Filtered$lat <- str_extract_all(eBird_Filtered$GridID, "([-]*[0-9]+[x])")
eBird_Filtered$lat <- as.numeric(str_remove_all(eBird_Filtered$lat, "[x]"))
eBird_Filtered$long <- as.numeric(str_extract_all(eBird_Filtered$GridID, "[0-9]+$"))

# -----------------------------------------------
# Calculate species shifts
# -----------------------------------------------
# Get unique species names
x <- unique(eBird_Filtered$Common_Name)
MAD_Shift_list <- list()
for (i in 1:length(x)) {
  # i = 1
  eBird_subset <- eBird_Filtered[eBird_Filtered$Common_Name == x[i] ,]
  eBird_subset_lmer <- lmer(MAD ~ Year + (1|GridID), data = eBird_subset)

  # Get dists.inv
  dists <- as.matrix(dist(cbind(eBird_subset$long, eBird_subset$lat)))
  dists.inv <- 1/dists
  diag(dists.inv) <- 0
  dists.inv[is.infinite(dists.inv)] <- 0
  
  # Perform Moran's I test
  Moran_result <- Moran.I(x = residuals(eBird_subset_lmer), weight = dists.inv)
  
  
  eBird_subset_ci <- confint.merMod(eBird_subset_lmer)
  MAD_Shift <- data.frame(unique(eBird_subset$Common_Name),
                          length(unique(eBird_subset$GridID)),
                          fixef(eBird_subset_lmer)[[2]],
                          eBird_subset_ci[4,1],
                          eBird_subset_ci[4,2],
                          # Moran data
                          Moran_result$observed,
                          Moran_result$expected,
                          Moran_result$sd,
                          Moran_result$p.value
                          )
  colnames(MAD_Shift) <- c("Common_Name", "n_Grids", "MAD_Shift", "MAD_Shift_2_5_pct", "MAD_Shift_97_5_pct", "Morans_I_observed", "Morans_I_expected", "Morans_I_sd", "Morans_I_P_value")
  MAD_Shift_list[[i]] <- MAD_Shift
}

Bird_Shifts <- bind_rows(MAD_Shift_list)

write_csv(Bird_Shifts, "Bird_Shifts_Morans_I.csv")

Bird_Shifts %>% 
  count(Morans_I_P_value < 0.05)
summary(Bird_Shifts$Morans_I_P_value)


# ---------------------------
# Moran's I: MAD_Shift
# ---------------------------
# Get grid-specific shifts
Grid_coefficients <- function (gridid){
  eBird_Filtered <- eBird_Filtered[eBird_Filtered$GridID == gridid ,]
  regression <- lmer(MAD ~ Year + (1|Common_Name), data = eBird_Filtered)
  MAD.Shift <- data.frame(fixef(regression))
  grid.id <- unique(eBird_Filtered$GridID)
  data.frame(grid.id, MAD.Shift[2,])
}

# Create dataframe of slopes == MAD Shifts
Grids <- unique(eBird_Filtered$GridID)
Grid_Shifts_lmer <- lapply(Grids, Grid_coefficients)   
Grid_Shifts_lmer <- data.frame(do.call(rbind, Grid_Shifts_lmer))
colnames(Grid_Shifts_lmer) <- c("GridID", "MAD_Shift_3")

# Extract lat long from the grids
eBird_Filtered$lat <- str_extract_all(eBird_Filtered$GridID, "([-]*[0-9]+[x])")
eBird_Filtered$lat <- as.numeric(str_remove_all(eBird_Filtered$lat, "[x]"))
eBird_Filtered$long <- as.numeric(str_extract_all(eBird_Filtered$GridID, "[0-9]+$"))
eBird_coords <- data.frame(x = eBird_Filtered$long, y = eBird_Filtered$lat)
# Make spatialpointsdataframe
eBird_sp <- SpatialPointsDataFrame(coords = eBird_coords, data = eBird_Filtered)
proj4string(eBird_sp) <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
# Set radius
radius = 100000
# define the plot edges based upon the plot radius.
unique_Grids <- eBird_Filtered %>%
  select(GridID, lat, long) %>%
  distinct()

yPlus <- unique_Grids$lat+radius
xPlus <- unique_Grids$long+radius
yMinus <- unique_Grids$lat-radius
xMinus <- unique_Grids$long-radius
# calculate polygon coordinates for each plot centroid.
square = cbind(xMinus, yPlus,  # NW corner
               xPlus, yPlus,  # NE corner
               xPlus, yMinus,  # SE corner
               xMinus, yMinus, # SW corner
               xMinus, yPlus)  # NW corner again - close ploygon
# Extract the plot ID information
ID = unique_Grids$GridID

# create spatial polygons from coordinates
polys <- SpatialPolygons(mapply(function(poly, id)
{
  xy <- matrix(poly, ncol=2, byrow=TRUE)
  Polygons(list(Polygon(xy)), ID=id)
},
split(square, row(square)), ID),
proj4string = CRS(as.character("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")))
# Make into spdf
polys <- SpatialPolygonsDataFrame(Sr = polys, data = Grid_Shifts_lmer, match.ID = F)

# We need to figure out which polygons are adjacent to each other
adjacent_polys <- poly2nb(polys, row.names = polys$GridID)
class(adjacent_polys)
summary(adjacent_polys)
summary(adjacent_polys)
str(adjacent_polys)

ww <- nb2listw(adjacent_polys, style = "B")
# Now we can use the moran function. Have a look at ?moran. The function is defined as ‘moran(y, ww, n, Szero(ww))’. Note the odd arguments n and S0. I think they are odd, because “ww” has that information. Anyway, we supply them and it works. There probably are cases where it makes sense to use other values.
moran(polys$MAD_Shift_3, ww, n = length(ww$neighbours), S0 = Szero(ww))
# Now we can test for significance. First analytically, using linear regression based logic and assumptions.
moran.test(polys$MAD_Shift_3, ww, randomisation = FALSE)

MC_MAD_Shift <- moran.mc(polys$MAD_Shift_3, ww, nsim = 999)
MC_MAD_Shift
# Plot the distribution (note that this is a density plot instead of a histogram)
plot(MC_MAD_Shift, main="", las=1)
