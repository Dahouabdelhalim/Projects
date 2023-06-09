# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 2
# ------------------------------------------------

options(scipen = 999)

# Load Libraries
library(dplyr)
library(readr)

# Set Working Directory
setwd("~/Data Repository")

# ------------------------------------------------
# Create grids
# ------------------------------------------------
# These are centroids for the X and Y coords
Xs <- seq(from = 100000, to = 2900000, by = 200000);nX <- length(Xs)
Ys <- seq(from = -300000, to = 2900000, by = 200000);nY <- length(Ys)
dX <- data.frame(X=Xs,dist=NA)
dY <- data.frame(Y=Ys,dist=NA)

# ------------------------------------------------
# Read in datasets
# ------------------------------------------------
# ------------------------------------------------
# Connecticut
# ------------------------------------------------
# Read in regional file
Connecticut_1 <- read_csv("Routine 1/Connecticut_1.csv")
# Get
Connecticut_1_sites <- Connecticut_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Connecticut_1 <- Connecticut_1[1:100000,]
nsites <- length(Connecticut_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Connecticut_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Connecticut_1_sites$iX <- NA
Connecticut_1_sites$iY <- NA
Connecticut_1_sites$GY <- NA
Connecticut_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Connecticut_1_sites$X[i])
  Connecticut_1_sites$iX[i] <- which.min(dX$dist)
  Connecticut_1_sites$GX[i] <- dX$X[Connecticut_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Connecticut_1_sites$Y[i])
  Connecticut_1_sites$iY[i] <- which.min(dY$dist)
  Connecticut_1_sites$GY[i] <- dY$Y[Connecticut_1_sites$iY[i]]
  p1 <- c(Connecticut_1_sites$X[i],Connecticut_1_sites$Y[i])
  p2 <- c(Connecticut_1_sites$GX[i],Connecticut_1_sites$GY[i])
  Connecticut_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Connecticut_1_sites$GridID <- paste("Grid_",Connecticut_1_sites$GY,"*",Connecticut_1_sites$GX,sep="")
# Copy sited results to gsited
gConnecticut_1_sites <- Connecticut_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uConnecticut_1_sites <- unique(gConnecticut_1_sites)
# Clean-up sited
Connecticut_1_sites$iY <- Connecticut_1_sites$iX <- Connecticut_1_sites$GY <- Connecticut_1_sites$GX <- NULL
# Join dataframes
Connecticut_2 <- left_join(Connecticut_1, Connecticut_1_sites)
# Write dataframes
write_csv(Connecticut_2, "Routine 2/Connecticut_2.csv")
write_csv(uConnecticut_1_sites, "Routine 2/Connecticut_2_UGrids.csv")
# Clean up dataframes
Connecticut_1 <- Connecticut_2 <- gConnecticut_1_sites <- Connecticut_1_sites <- uConnecticut_1_sites <- NULL

# ------------------------------------------------
# Delaware
# ------------------------------------------------
# Read in regional file
Delaware_1 <- read_csv("Routine 1/Delaware_1.csv")
# Get
Delaware_1_sites <- Delaware_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Delaware_1 <- Delaware_1[1:100000,]
nsites <- length(Delaware_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Delaware_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Delaware_1_sites$iX <- NA
Delaware_1_sites$iY <- NA
Delaware_1_sites$GY <- NA
Delaware_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Delaware_1_sites$X[i])
  Delaware_1_sites$iX[i] <- which.min(dX$dist)
  Delaware_1_sites$GX[i] <- dX$X[Delaware_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Delaware_1_sites$Y[i])
  Delaware_1_sites$iY[i] <- which.min(dY$dist)
  Delaware_1_sites$GY[i] <- dY$Y[Delaware_1_sites$iY[i]]
  p1 <- c(Delaware_1_sites$X[i],Delaware_1_sites$Y[i])
  p2 <- c(Delaware_1_sites$GX[i],Delaware_1_sites$GY[i])
  Delaware_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Delaware_1_sites$GridID <- paste("Grid_",Delaware_1_sites$GY,"*",Delaware_1_sites$GX,sep="")
# Copy sited results to gsited
gDelaware_1_sites <- Delaware_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uDelaware_1_sites <- unique(gDelaware_1_sites)
# Clean-up sited
Delaware_1_sites$iY <- Delaware_1_sites$iX <- Delaware_1_sites$GY <- Delaware_1_sites$GX <- NULL
# Join dataframes
Delaware_2 <- left_join(Delaware_1, Delaware_1_sites)
# Write dataframes
write_csv(Delaware_2, "Routine 2/Delaware_2.csv")
write_csv(uDelaware_1_sites, "Routine 2/Delaware_2_UGrids.csv")
# Clean up dataframes
Delaware_1 <- Delaware_2 <- gDelaware_1_sites <- Delaware_1_sites <- uDelaware_1_sites <- NULL

# ------------------------------------------------
# Illinois
# ------------------------------------------------
# Read in regional file
Illinois_1 <- read_csv("Routine 1/Illinois_1.csv")
# Get
Illinois_1_sites <- Illinois_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Illinois_1 <- Illinois_1[1:100000,]
nsites <- length(Illinois_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Illinois_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Illinois_1_sites$iX <- NA
Illinois_1_sites$iY <- NA
Illinois_1_sites$GY <- NA
Illinois_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Illinois_1_sites$X[i])
  Illinois_1_sites$iX[i] <- which.min(dX$dist)
  Illinois_1_sites$GX[i] <- dX$X[Illinois_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Illinois_1_sites$Y[i])
  Illinois_1_sites$iY[i] <- which.min(dY$dist)
  Illinois_1_sites$GY[i] <- dY$Y[Illinois_1_sites$iY[i]]
  p1 <- c(Illinois_1_sites$X[i],Illinois_1_sites$Y[i])
  p2 <- c(Illinois_1_sites$GX[i],Illinois_1_sites$GY[i])
  Illinois_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Illinois_1_sites$GridID <- paste("Grid_",Illinois_1_sites$GY,"*",Illinois_1_sites$GX,sep="")
# Copy sited results to gsited
gIllinois_1_sites <- Illinois_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uIllinois_1_sites <- unique(gIllinois_1_sites)
# Clean-up sited
Illinois_1_sites$iY <- Illinois_1_sites$iX <- Illinois_1_sites$GY <- Illinois_1_sites$GX <- NULL
# Join dataframes
Illinois_2 <- left_join(Illinois_1, Illinois_1_sites)
# Write dataframes
write_csv(Illinois_2, "Routine 2/Illinois_2.csv")
write_csv(uIllinois_1_sites, "Routine 2/Illinois_2_UGrids.csv")
# Clean up dataframes
Illinois_1 <- Illinois_2 <- gIllinois_1_sites <- Illinois_1_sites <- uIllinois_1_sites <- NULL

# ------------------------------------------------
# Indiana
# ------------------------------------------------
# Read in regional file
Indiana_1 <- read_csv("Routine 1/Indiana_1.csv")
# Get
Indiana_1_sites <- Indiana_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Indiana_1 <- Indiana_1[1:100000,]
nsites <- length(Indiana_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Indiana_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Indiana_1_sites$iX <- NA
Indiana_1_sites$iY <- NA
Indiana_1_sites$GY <- NA
Indiana_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Indiana_1_sites$X[i])
  Indiana_1_sites$iX[i] <- which.min(dX$dist)
  Indiana_1_sites$GX[i] <- dX$X[Indiana_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Indiana_1_sites$Y[i])
  Indiana_1_sites$iY[i] <- which.min(dY$dist)
  Indiana_1_sites$GY[i] <- dY$Y[Indiana_1_sites$iY[i]]
  p1 <- c(Indiana_1_sites$X[i],Indiana_1_sites$Y[i])
  p2 <- c(Indiana_1_sites$GX[i],Indiana_1_sites$GY[i])
  Indiana_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Indiana_1_sites$GridID <- paste("Grid_",Indiana_1_sites$GY,"*",Indiana_1_sites$GX,sep="")
# Copy sited results to gsited
gIndiana_1_sites <- Indiana_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uIndiana_1_sites <- unique(gIndiana_1_sites)
# Clean-up sited
Indiana_1_sites$iY <- Indiana_1_sites$iX <- Indiana_1_sites$GY <- Indiana_1_sites$GX <- NULL
# Join dataframes
Indiana_2 <- left_join(Indiana_1, Indiana_1_sites)
# Write dataframes
write_csv(Indiana_2, "Routine 2/Indiana_2.csv")
write_csv(uIndiana_1_sites, "Routine 2/Indiana_2_UGrids.csv")
# Clean up dataframes
Indiana_1 <- Indiana_2 <- gIndiana_1_sites <- Indiana_1_sites <- uIndiana_1_sites <- NULL

# ------------------------------------------------
# Kentucky
# ------------------------------------------------
# Read in regional file
Kentucky_1 <- read_csv("Routine 1/Kentucky_1.csv")
# Get
Kentucky_1_sites <- Kentucky_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Kentucky_1 <- Kentucky_1[1:100000,]
nsites <- length(Kentucky_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Kentucky_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Kentucky_1_sites$iX <- NA
Kentucky_1_sites$iY <- NA
Kentucky_1_sites$GY <- NA
Kentucky_1_sites$GX <- NA

for (i in 1:nsites) {
    dX$dist <- abs(dX$X - Kentucky_1_sites$X[i])
    Kentucky_1_sites$iX[i] <- which.min(dX$dist)
    Kentucky_1_sites$GX[i] <- dX$X[Kentucky_1_sites$iX[i]]
    dY$dist <- abs(dY$Y - Kentucky_1_sites$Y[i])
    Kentucky_1_sites$iY[i] <- which.min(dY$dist)
    Kentucky_1_sites$GY[i] <- dY$Y[Kentucky_1_sites$iY[i]]
    p1 <- c(Kentucky_1_sites$X[i],Kentucky_1_sites$Y[i])
    p2 <- c(Kentucky_1_sites$GX[i],Kentucky_1_sites$GY[i])
    Kentucky_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
  }

Kentucky_1_sites$GridID <- paste("Grid_",Kentucky_1_sites$GY,"*",Kentucky_1_sites$GX,sep="")
# Copy sited results to gsited
gKentucky_1_sites <- Kentucky_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uKentucky_1_sites <- unique(gKentucky_1_sites)
# Clean-up sited
Kentucky_1_sites$iY <- Kentucky_1_sites$iX <- Kentucky_1_sites$GY <- Kentucky_1_sites$GX <- NULL
# Join dataframes
Kentucky_2 <- left_join(Kentucky_1, Kentucky_1_sites)
# Write dataframes
write_csv(Kentucky_2, "Routine 2/Kentucky_2.csv")
write_csv(uKentucky_1_sites, "Routine 2/Kentucky_2_UGrids.csv")
# Clean up dataframes
Kentucky_1 <- Kentucky_2 <- gKentucky_1_sites <- Kentucky_1_sites <- uKentucky_1_sites <- NULL

# ------------------------------------------------
# Maine
# ------------------------------------------------
# Read in regional file
Maine_1 <- read_csv("Routine 1/Maine_1.csv")
# Get
Maine_1_sites <- Maine_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Maine_1 <- Maine_1[1:100000,]
nsites <- length(Maine_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Maine_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Maine_1_sites$iX <- NA
Maine_1_sites$iY <- NA
Maine_1_sites$GY <- NA
Maine_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Maine_1_sites$X[i])
  Maine_1_sites$iX[i] <- which.min(dX$dist)
  Maine_1_sites$GX[i] <- dX$X[Maine_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Maine_1_sites$Y[i])
  Maine_1_sites$iY[i] <- which.min(dY$dist)
  Maine_1_sites$GY[i] <- dY$Y[Maine_1_sites$iY[i]]
  p1 <- c(Maine_1_sites$X[i],Maine_1_sites$Y[i])
  p2 <- c(Maine_1_sites$GX[i],Maine_1_sites$GY[i])
  Maine_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Maine_1_sites$GridID <- paste("Grid_",Maine_1_sites$GY,"*",Maine_1_sites$GX,sep="")
# Copy sited results to gsited
gMaine_1_sites <- Maine_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uMaine_1_sites <- unique(gMaine_1_sites)
# Clean-up sited
Maine_1_sites$iY <- Maine_1_sites$iX <- Maine_1_sites$GY <- Maine_1_sites$GX <- NULL
# Join dataframes
Maine_2 <- left_join(Maine_1, Maine_1_sites)
# Write dataframes
write_csv(Maine_2, "Routine 2/Maine_2.csv")
write_csv(uMaine_1_sites, "Routine 2/Maine_2_UGrids.csv")
# Clean up dataframes
Maine_1 <- Maine_2 <- gMaine_1_sites <- Maine_1_sites <- uMaine_1_sites <- NULL

# ------------------------------------------------
# Maryland
# ------------------------------------------------
# Read in regional file
Maryland_1 <- read_csv("Routine 1/Maryland_1.csv")
# Get
Maryland_1_sites <- Maryland_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Maryland_1 <- Maryland_1[1:100000,]
nsites <- length(Maryland_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Maryland_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Maryland_1_sites$iX <- NA
Maryland_1_sites$iY <- NA
Maryland_1_sites$GY <- NA
Maryland_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Maryland_1_sites$X[i])
  Maryland_1_sites$iX[i] <- which.min(dX$dist)
  Maryland_1_sites$GX[i] <- dX$X[Maryland_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Maryland_1_sites$Y[i])
  Maryland_1_sites$iY[i] <- which.min(dY$dist)
  Maryland_1_sites$GY[i] <- dY$Y[Maryland_1_sites$iY[i]]
  p1 <- c(Maryland_1_sites$X[i],Maryland_1_sites$Y[i])
  p2 <- c(Maryland_1_sites$GX[i],Maryland_1_sites$GY[i])
  Maryland_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Maryland_1_sites$GridID <- paste("Grid_",Maryland_1_sites$GY,"*",Maryland_1_sites$GX,sep="")
# Copy sited results to gsited
gMaryland_1_sites <- Maryland_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uMaryland_1_sites <- unique(gMaryland_1_sites)
# Clean-up sited
Maryland_1_sites$iY <- Maryland_1_sites$iX <- Maryland_1_sites$GY <- Maryland_1_sites$GX <- NULL
# Join dataframes
Maryland_2 <- left_join(Maryland_1, Maryland_1_sites)
# Write dataframes
write_csv(Maryland_2, "Routine 2/Maryland_2.csv")
write_csv(uMaryland_1_sites, "Routine 2/Maryland_2_UGrids.csv")
# Clean up dataframes
Maryland_1 <- Maryland_2 <- gMaryland_1_sites <- Maryland_1_sites <- uMaryland_1_sites <- NULL

# ------------------------------------------------
# Massachusetts
# ------------------------------------------------
# Read in regional file
Massachusetts_1 <- read_csv("Routine 1/Massachusetts_1.csv")
# Get
Massachusetts_1_sites <- Massachusetts_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Massachusetts_1 <- Massachusetts_1[1:100000,]
nsites <- length(Massachusetts_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Massachusetts_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Massachusetts_1_sites$iX <- NA
Massachusetts_1_sites$iY <- NA
Massachusetts_1_sites$GY <- NA
Massachusetts_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Massachusetts_1_sites$X[i])
  Massachusetts_1_sites$iX[i] <- which.min(dX$dist)
  Massachusetts_1_sites$GX[i] <- dX$X[Massachusetts_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Massachusetts_1_sites$Y[i])
  Massachusetts_1_sites$iY[i] <- which.min(dY$dist)
  Massachusetts_1_sites$GY[i] <- dY$Y[Massachusetts_1_sites$iY[i]]
  p1 <- c(Massachusetts_1_sites$X[i],Massachusetts_1_sites$Y[i])
  p2 <- c(Massachusetts_1_sites$GX[i],Massachusetts_1_sites$GY[i])
  Massachusetts_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Massachusetts_1_sites$GridID <- paste("Grid_",Massachusetts_1_sites$GY,"*",Massachusetts_1_sites$GX,sep="")
# Copy sited results to gsited
gMassachusetts_1_sites <- Massachusetts_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uMassachusetts_1_sites <- unique(gMassachusetts_1_sites)
# Clean-up sited
Massachusetts_1_sites$iY <- Massachusetts_1_sites$iX <- Massachusetts_1_sites$GY <- Massachusetts_1_sites$GX <- NULL
# Join dataframes
Massachusetts_2 <- left_join(Massachusetts_1, Massachusetts_1_sites)
# Write dataframes
write_csv(Massachusetts_2, "Routine 2/Massachusetts_2.csv")
write_csv(uMassachusetts_1_sites, "Routine 2/Massachusetts_2_UGrids.csv")
# Clean up dataframes
Massachusetts_1 <- Massachusetts_2 <- gMassachusetts_1_sites <- Massachusetts_1_sites <- uMassachusetts_1_sites <- NULL

# ------------------------------------------------
# Michigan
# ------------------------------------------------
# Read in regional file
Michigan_1 <- read_csv("Routine 1/Michigan_1.csv")
# Get
Michigan_1_sites <- Michigan_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Michigan_1 <- Michigan_1[1:100000,]
nsites <- length(Michigan_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Michigan_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Michigan_1_sites$iX <- NA
Michigan_1_sites$iY <- NA
Michigan_1_sites$GY <- NA
Michigan_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Michigan_1_sites$X[i])
  Michigan_1_sites$iX[i] <- which.min(dX$dist)
  Michigan_1_sites$GX[i] <- dX$X[Michigan_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Michigan_1_sites$Y[i])
  Michigan_1_sites$iY[i] <- which.min(dY$dist)
  Michigan_1_sites$GY[i] <- dY$Y[Michigan_1_sites$iY[i]]
  p1 <- c(Michigan_1_sites$X[i],Michigan_1_sites$Y[i])
  p2 <- c(Michigan_1_sites$GX[i],Michigan_1_sites$GY[i])
  Michigan_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Michigan_1_sites$GridID <- paste("Grid_",Michigan_1_sites$GY,"*",Michigan_1_sites$GX,sep="")
# Copy sited results to gsited
gMichigan_1_sites <- Michigan_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uMichigan_1_sites <- unique(gMichigan_1_sites)
# Clean-up sited
Michigan_1_sites$iY <- Michigan_1_sites$iX <- Michigan_1_sites$GY <- Michigan_1_sites$GX <- NULL
# Join dataframes
Michigan_2 <- left_join(Michigan_1, Michigan_1_sites)
# Write dataframes
write_csv(Michigan_2, "Routine 2/Michigan_2.csv")
write_csv(uMichigan_1_sites, "Routine 2/Michigan_2_UGrids.csv")
# Clean up dataframes
Michigan_1 <- Michigan_2 <- gMichigan_1_sites <- Michigan_1_sites <- uMichigan_1_sites <- NULL

# ------------------------------------------------
# New_Brunswick
# ------------------------------------------------
# Read in regional file
New_Brunswick_1 <- read_csv("Routine 1/New_Brunswick_1.csv")
# Get
New_Brunswick_1_sites <- New_Brunswick_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# New_Brunswick_1 <- New_Brunswick_1[1:100000,]
nsites <- length(New_Brunswick_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
New_Brunswick_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
New_Brunswick_1_sites$iX <- NA
New_Brunswick_1_sites$iY <- NA
New_Brunswick_1_sites$GY <- NA
New_Brunswick_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - New_Brunswick_1_sites$X[i])
  New_Brunswick_1_sites$iX[i] <- which.min(dX$dist)
  New_Brunswick_1_sites$GX[i] <- dX$X[New_Brunswick_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - New_Brunswick_1_sites$Y[i])
  New_Brunswick_1_sites$iY[i] <- which.min(dY$dist)
  New_Brunswick_1_sites$GY[i] <- dY$Y[New_Brunswick_1_sites$iY[i]]
  p1 <- c(New_Brunswick_1_sites$X[i],New_Brunswick_1_sites$Y[i])
  p2 <- c(New_Brunswick_1_sites$GX[i],New_Brunswick_1_sites$GY[i])
  New_Brunswick_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

New_Brunswick_1_sites$GridID <- paste("Grid_",New_Brunswick_1_sites$GY,"*",New_Brunswick_1_sites$GX,sep="")
# Copy sited results to gsited
gNew_Brunswick_1_sites <- New_Brunswick_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uNew_Brunswick_1_sites <- unique(gNew_Brunswick_1_sites)
# Clean-up sited
New_Brunswick_1_sites$iY <- New_Brunswick_1_sites$iX <- New_Brunswick_1_sites$GY <- New_Brunswick_1_sites$GX <- NULL
# Join dataframes
New_Brunswick_2 <- left_join(New_Brunswick_1, New_Brunswick_1_sites)
# Write dataframes
write_csv(New_Brunswick_2, "Routine 2/New_Brunswick_2.csv")
write_csv(uNew_Brunswick_1_sites, "Routine 2/New_Brunswick_2_UGrids.csv")
# Clean up dataframes
New_Brunswick_1 <- New_Brunswick_2 <- gNew_Brunswick_1_sites <- New_Brunswick_1_sites <- uNew_Brunswick_1_sites <- NULL

# ------------------------------------------------
# New_Hampshire
# ------------------------------------------------
# Read in regional file
New_Hampshire_1 <- read_csv("Routine 1/New_Hampshire_1.csv")
# Get
New_Hampshire_1_sites <- New_Hampshire_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# New_Hampshire_1 <- New_Hampshire_1[1:100000,]
nsites <- length(New_Hampshire_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
New_Hampshire_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
New_Hampshire_1_sites$iX <- NA
New_Hampshire_1_sites$iY <- NA
New_Hampshire_1_sites$GY <- NA
New_Hampshire_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - New_Hampshire_1_sites$X[i])
  New_Hampshire_1_sites$iX[i] <- which.min(dX$dist)
  New_Hampshire_1_sites$GX[i] <- dX$X[New_Hampshire_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - New_Hampshire_1_sites$Y[i])
  New_Hampshire_1_sites$iY[i] <- which.min(dY$dist)
  New_Hampshire_1_sites$GY[i] <- dY$Y[New_Hampshire_1_sites$iY[i]]
  p1 <- c(New_Hampshire_1_sites$X[i],New_Hampshire_1_sites$Y[i])
  p2 <- c(New_Hampshire_1_sites$GX[i],New_Hampshire_1_sites$GY[i])
  New_Hampshire_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

New_Hampshire_1_sites$GridID <- paste("Grid_",New_Hampshire_1_sites$GY,"*",New_Hampshire_1_sites$GX,sep="")
# Copy sited results to gsited
gNew_Hampshire_1_sites <- New_Hampshire_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uNew_Hampshire_1_sites <- unique(gNew_Hampshire_1_sites)
# Clean-up sited
New_Hampshire_1_sites$iY <- New_Hampshire_1_sites$iX <- New_Hampshire_1_sites$GY <- New_Hampshire_1_sites$GX <- NULL
# Join dataframes
New_Hampshire_2 <- left_join(New_Hampshire_1, New_Hampshire_1_sites)
# Write dataframes
write_csv(New_Hampshire_2, "Routine 2/New_Hampshire_2.csv")
write_csv(uNew_Hampshire_1_sites, "Routine 2/New_Hampshire_2_UGrids.csv")
# Clean up dataframes
New_Hampshire_1 <- New_Hampshire_2 <- gNew_Hampshire_1_sites <- New_Hampshire_1_sites <- uNew_Hampshire_1_sites <- NULL

# ------------------------------------------------
# New_Jersey
# ------------------------------------------------
# Read in regional file
New_Jersey_1 <- read_csv("Routine 1/New_Jersey_1.csv")
# Get
New_Jersey_1_sites <- New_Jersey_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# New_Jersey_1 <- New_Jersey_1[1:100000,]
nsites <- length(New_Jersey_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
New_Jersey_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
New_Jersey_1_sites$iX <- NA
New_Jersey_1_sites$iY <- NA
New_Jersey_1_sites$GY <- NA
New_Jersey_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - New_Jersey_1_sites$X[i])
  New_Jersey_1_sites$iX[i] <- which.min(dX$dist)
  New_Jersey_1_sites$GX[i] <- dX$X[New_Jersey_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - New_Jersey_1_sites$Y[i])
  New_Jersey_1_sites$iY[i] <- which.min(dY$dist)
  New_Jersey_1_sites$GY[i] <- dY$Y[New_Jersey_1_sites$iY[i]]
  p1 <- c(New_Jersey_1_sites$X[i],New_Jersey_1_sites$Y[i])
  p2 <- c(New_Jersey_1_sites$GX[i],New_Jersey_1_sites$GY[i])
  New_Jersey_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

New_Jersey_1_sites$GridID <- paste("Grid_",New_Jersey_1_sites$GY,"*",New_Jersey_1_sites$GX,sep="")
# Copy sited results to gsited
gNew_Jersey_1_sites <- New_Jersey_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uNew_Jersey_1_sites <- unique(gNew_Jersey_1_sites)
# Clean-up sited
New_Jersey_1_sites$iY <- New_Jersey_1_sites$iX <- New_Jersey_1_sites$GY <- New_Jersey_1_sites$GX <- NULL
# Join dataframes
New_Jersey_2 <- left_join(New_Jersey_1, New_Jersey_1_sites)
# Write dataframes
write_csv(New_Jersey_2, "Routine 2/New_Jersey_2.csv")
write_csv(uNew_Jersey_1_sites, "Routine 2/New_Jersey_2_UGrids.csv")
# Clean up dataframes
New_Jersey_1 <- New_Jersey_2 <- gNew_Jersey_1_sites <- New_Jersey_1_sites <- uNew_Jersey_1_sites <- NULL

# ------------------------------------------------
# New_York
# ------------------------------------------------
# Read in regional file
New_York_1 <- read_csv("Routine 1/New_York_1.csv")
# Get
New_York_1_sites <- New_York_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# New_York_1 <- New_York_1[1:100000,]
nsites <- length(New_York_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
New_York_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
New_York_1_sites$iX <- NA
New_York_1_sites$iY <- NA
New_York_1_sites$GY <- NA
New_York_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - New_York_1_sites$X[i])
  New_York_1_sites$iX[i] <- which.min(dX$dist)
  New_York_1_sites$GX[i] <- dX$X[New_York_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - New_York_1_sites$Y[i])
  New_York_1_sites$iY[i] <- which.min(dY$dist)
  New_York_1_sites$GY[i] <- dY$Y[New_York_1_sites$iY[i]]
  p1 <- c(New_York_1_sites$X[i],New_York_1_sites$Y[i])
  p2 <- c(New_York_1_sites$GX[i],New_York_1_sites$GY[i])
  New_York_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

New_York_1_sites$GridID <- paste("Grid_",New_York_1_sites$GY,"*",New_York_1_sites$GX,sep="")
# Copy sited results to gsited
gNew_York_1_sites <- New_York_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uNew_York_1_sites <- unique(gNew_York_1_sites)
# Clean-up sited
New_York_1_sites$iY <- New_York_1_sites$iX <- New_York_1_sites$GY <- New_York_1_sites$GX <- NULL
# Join dataframes
New_York_2 <- left_join(New_York_1, New_York_1_sites)
# Write dataframes
write_csv(New_York_2, "Routine 2/New_York_2.csv")
write_csv(uNew_York_1_sites, "Routine 2/New_York_2_UGrids.csv")
# Clean up dataframes
New_York_1 <- New_York_2 <- gNew_York_1_sites <- New_York_1_sites <- uNew_York_1_sites <- NULL

# ------------------------------------------------
# Nova_Scotia
# ------------------------------------------------
# Read in regional file
Nova_Scotia_1 <- read_csv("Routine 1/Nova_Scotia_1.csv")
# Get
Nova_Scotia_1_sites <- Nova_Scotia_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Nova_Scotia_1 <- Nova_Scotia_1[1:100000,]
nsites <- length(Nova_Scotia_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Nova_Scotia_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Nova_Scotia_1_sites$iX <- NA
Nova_Scotia_1_sites$iY <- NA
Nova_Scotia_1_sites$GY <- NA
Nova_Scotia_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Nova_Scotia_1_sites$X[i])
  Nova_Scotia_1_sites$iX[i] <- which.min(dX$dist)
  Nova_Scotia_1_sites$GX[i] <- dX$X[Nova_Scotia_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Nova_Scotia_1_sites$Y[i])
  Nova_Scotia_1_sites$iY[i] <- which.min(dY$dist)
  Nova_Scotia_1_sites$GY[i] <- dY$Y[Nova_Scotia_1_sites$iY[i]]
  p1 <- c(Nova_Scotia_1_sites$X[i],Nova_Scotia_1_sites$Y[i])
  p2 <- c(Nova_Scotia_1_sites$GX[i],Nova_Scotia_1_sites$GY[i])
  Nova_Scotia_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Nova_Scotia_1_sites$GridID <- paste("Grid_",Nova_Scotia_1_sites$GY,"*",Nova_Scotia_1_sites$GX,sep="")
# Copy sited results to gsited
gNova_Scotia_1_sites <- Nova_Scotia_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uNova_Scotia_1_sites <- unique(gNova_Scotia_1_sites)
# Clean-up sited
Nova_Scotia_1_sites$iY <- Nova_Scotia_1_sites$iX <- Nova_Scotia_1_sites$GY <- Nova_Scotia_1_sites$GX <- NULL
# Join dataframes
Nova_Scotia_2 <- left_join(Nova_Scotia_1, Nova_Scotia_1_sites)
# Write dataframes
write_csv(Nova_Scotia_2, "Routine 2/Nova_Scotia_2.csv")
write_csv(uNova_Scotia_1_sites, "Routine 2/Nova_Scotia_2_UGrids.csv")
# Clean up dataframes
Nova_Scotia_1 <- Nova_Scotia_2 <- gNova_Scotia_1_sites <- Nova_Scotia_1_sites <- uNova_Scotia_1_sites <- NULL

# ------------------------------------------------
# Ohio
# ------------------------------------------------
# Read in regional file
Ohio_1 <- read_csv("Routine 1/Ohio_1.csv")
# Get
Ohio_1_sites <- Ohio_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Ohio_1 <- Ohio_1[1:100000,]
nsites <- length(Ohio_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Ohio_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Ohio_1_sites$iX <- NA
Ohio_1_sites$iY <- NA
Ohio_1_sites$GY <- NA
Ohio_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Ohio_1_sites$X[i])
  Ohio_1_sites$iX[i] <- which.min(dX$dist)
  Ohio_1_sites$GX[i] <- dX$X[Ohio_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Ohio_1_sites$Y[i])
  Ohio_1_sites$iY[i] <- which.min(dY$dist)
  Ohio_1_sites$GY[i] <- dY$Y[Ohio_1_sites$iY[i]]
  p1 <- c(Ohio_1_sites$X[i],Ohio_1_sites$Y[i])
  p2 <- c(Ohio_1_sites$GX[i],Ohio_1_sites$GY[i])
  Ohio_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Ohio_1_sites$GridID <- paste("Grid_",Ohio_1_sites$GY,"*",Ohio_1_sites$GX,sep="")
# Copy sited results to gsited
gOhio_1_sites <- Ohio_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uOhio_1_sites <- unique(gOhio_1_sites)
# Clean-up sited
Ohio_1_sites$iY <- Ohio_1_sites$iX <- Ohio_1_sites$GY <- Ohio_1_sites$GX <- NULL
# Join dataframes
Ohio_2 <- left_join(Ohio_1, Ohio_1_sites)
# Write dataframes
write_csv(Ohio_2, "Routine 2/Ohio_2.csv")
write_csv(uOhio_1_sites, "Routine 2/Ohio_2_UGrids.csv")
# Clean up dataframes
Ohio_1 <- Ohio_2 <- gOhio_1_sites <- Ohio_1_sites <- uOhio_1_sites <- NULL



# ------------------------------------------------
# Ontario
# ------------------------------------------------
# Read in regional file
Ontario_1 <- read_csv("Routine 1/Ontario_1.csv")
# Get
Ontario_1_sites <- Ontario_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Ontario_1 <- Ontario_1[1:100000,]
nsites <- length(Ontario_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Ontario_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Ontario_1_sites$iX <- NA
Ontario_1_sites$iY <- NA
Ontario_1_sites$GY <- NA
Ontario_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Ontario_1_sites$X[i])
  Ontario_1_sites$iX[i] <- which.min(dX$dist)
  Ontario_1_sites$GX[i] <- dX$X[Ontario_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Ontario_1_sites$Y[i])
  Ontario_1_sites$iY[i] <- which.min(dY$dist)
  Ontario_1_sites$GY[i] <- dY$Y[Ontario_1_sites$iY[i]]
  p1 <- c(Ontario_1_sites$X[i],Ontario_1_sites$Y[i])
  p2 <- c(Ontario_1_sites$GX[i],Ontario_1_sites$GY[i])
  Ontario_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Ontario_1_sites$GridID <- paste("Grid_",Ontario_1_sites$GY,"*",Ontario_1_sites$GX,sep="")
# Copy sited results to gsited
gOntario_1_sites <- Ontario_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uOntario_1_sites <- unique(gOntario_1_sites)
# Clean-up sited
Ontario_1_sites$iY <- Ontario_1_sites$iX <- Ontario_1_sites$GY <- Ontario_1_sites$GX <- NULL
# Join dataframes
Ontario_2 <- left_join(Ontario_1, Ontario_1_sites)
# Write dataframes
write_csv(Ontario_2, "Routine 2/Ontario_2.csv")
write_csv(uOntario_1_sites, "Routine 2/Ontario_2_UGrids.csv")
# Clean up dataframes
Ontario_1 <- Ontario_2 <- gOntario_1_sites <- Ontario_1_sites <- uOntario_1_sites <- NULL

# ------------------------------------------------
# Pennsylvania
# ------------------------------------------------
# Read in regional file
Pennsylvania_1 <- read_csv("Routine 1/Pennsylvania_1.csv")
# Get
Pennsylvania_1_sites <- Pennsylvania_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Pennsylvania_1 <- Pennsylvania_1[1:100000,]
nsites <- length(Pennsylvania_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Pennsylvania_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Pennsylvania_1_sites$iX <- NA
Pennsylvania_1_sites$iY <- NA
Pennsylvania_1_sites$GY <- NA
Pennsylvania_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Pennsylvania_1_sites$X[i])
  Pennsylvania_1_sites$iX[i] <- which.min(dX$dist)
  Pennsylvania_1_sites$GX[i] <- dX$X[Pennsylvania_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Pennsylvania_1_sites$Y[i])
  Pennsylvania_1_sites$iY[i] <- which.min(dY$dist)
  Pennsylvania_1_sites$GY[i] <- dY$Y[Pennsylvania_1_sites$iY[i]]
  p1 <- c(Pennsylvania_1_sites$X[i],Pennsylvania_1_sites$Y[i])
  p2 <- c(Pennsylvania_1_sites$GX[i],Pennsylvania_1_sites$GY[i])
  Pennsylvania_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Pennsylvania_1_sites$GridID <- paste("Grid_",Pennsylvania_1_sites$GY,"*",Pennsylvania_1_sites$GX,sep="")
# Copy sited results to gsited
gPennsylvania_1_sites <- Pennsylvania_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uPennsylvania_1_sites <- unique(gPennsylvania_1_sites)
# Clean-up sited
Pennsylvania_1_sites$iY <- Pennsylvania_1_sites$iX <- Pennsylvania_1_sites$GY <- Pennsylvania_1_sites$GX <- NULL
# Join dataframes
Pennsylvania_2 <- left_join(Pennsylvania_1, Pennsylvania_1_sites)
# Write dataframes
write_csv(Pennsylvania_2, "Routine 2/Pennsylvania_2.csv")
write_csv(uPennsylvania_1_sites, "Routine 2/Pennsylvania_2_UGrids.csv")
# Clean up dataframes
Pennsylvania_1 <- Pennsylvania_2 <- gPennsylvania_1_sites <- Pennsylvania_1_sites <- uPennsylvania_1_sites <- NULL

# ------------------------------------------------
# Quebec
# ------------------------------------------------
# Read in regional file
Quebec_1 <- read_csv("Routine 1/Quebec_1.csv")
# Get
Quebec_1_sites <- Quebec_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Quebec_1 <- Quebec_1[1:100000,]
nsites <- length(Quebec_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Quebec_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Quebec_1_sites$iX <- NA
Quebec_1_sites$iY <- NA
Quebec_1_sites$GY <- NA
Quebec_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Quebec_1_sites$X[i])
  Quebec_1_sites$iX[i] <- which.min(dX$dist)
  Quebec_1_sites$GX[i] <- dX$X[Quebec_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Quebec_1_sites$Y[i])
  Quebec_1_sites$iY[i] <- which.min(dY$dist)
  Quebec_1_sites$GY[i] <- dY$Y[Quebec_1_sites$iY[i]]
  p1 <- c(Quebec_1_sites$X[i],Quebec_1_sites$Y[i])
  p2 <- c(Quebec_1_sites$GX[i],Quebec_1_sites$GY[i])
  Quebec_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Quebec_1_sites$GridID <- paste("Grid_",Quebec_1_sites$GY,"*",Quebec_1_sites$GX,sep="")
# Copy sited results to gsited
gQuebec_1_sites <- Quebec_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uQuebec_1_sites <- unique(gQuebec_1_sites)
# Clean-up sited
Quebec_1_sites$iY <- Quebec_1_sites$iX <- Quebec_1_sites$GY <- Quebec_1_sites$GX <- NULL
# Join dataframes
Quebec_2 <- left_join(Quebec_1, Quebec_1_sites)
# Write dataframes
write_csv(Quebec_2, "Routine 2/Quebec_2.csv")
write_csv(uQuebec_1_sites, "Routine 2/Quebec_2_UGrids.csv")
# Clean up dataframes
Quebec_1 <- Quebec_2 <- gQuebec_1_sites <- Quebec_1_sites <- uQuebec_1_sites <- NULL

# ------------------------------------------------
# Rhode_Island
# ------------------------------------------------
# Read in regional file
Rhode_Island_1 <- read_csv("Routine 1/Rhode_Island_1.csv")
# Get
Rhode_Island_1_sites <- Rhode_Island_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Rhode_Island_1 <- Rhode_Island_1[1:100000,]
nsites <- length(Rhode_Island_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Rhode_Island_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Rhode_Island_1_sites$iX <- NA
Rhode_Island_1_sites$iY <- NA
Rhode_Island_1_sites$GY <- NA
Rhode_Island_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Rhode_Island_1_sites$X[i])
  Rhode_Island_1_sites$iX[i] <- which.min(dX$dist)
  Rhode_Island_1_sites$GX[i] <- dX$X[Rhode_Island_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Rhode_Island_1_sites$Y[i])
  Rhode_Island_1_sites$iY[i] <- which.min(dY$dist)
  Rhode_Island_1_sites$GY[i] <- dY$Y[Rhode_Island_1_sites$iY[i]]
  p1 <- c(Rhode_Island_1_sites$X[i],Rhode_Island_1_sites$Y[i])
  p2 <- c(Rhode_Island_1_sites$GX[i],Rhode_Island_1_sites$GY[i])
  Rhode_Island_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Rhode_Island_1_sites$GridID <- paste("Grid_",Rhode_Island_1_sites$GY,"*",Rhode_Island_1_sites$GX,sep="")
# Copy sited results to gsited
gRhode_Island_1_sites <- Rhode_Island_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uRhode_Island_1_sites <- unique(gRhode_Island_1_sites)
# Clean-up sited
Rhode_Island_1_sites$iY <- Rhode_Island_1_sites$iX <- Rhode_Island_1_sites$GY <- Rhode_Island_1_sites$GX <- NULL
# Join dataframes
Rhode_Island_2 <- left_join(Rhode_Island_1, Rhode_Island_1_sites)
# Write dataframes
write_csv(Rhode_Island_2, "Routine 2/Rhode_Island_2.csv")
write_csv(uRhode_Island_1_sites, "Routine 2/Rhode_Island_2_UGrids.csv")
# Clean up dataframes
Rhode_Island_1 <- Rhode_Island_2 <- gRhode_Island_1_sites <- Rhode_Island_1_sites <- uRhode_Island_1_sites <- NULL

# ------------------------------------------------
# Vermont
# ------------------------------------------------
# Read in regional file
Vermont_1 <- read_csv("Routine 1/Vermont_1.csv")
# Get
Vermont_1_sites <- Vermont_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Vermont_1 <- Vermont_1[1:100000,]
nsites <- length(Vermont_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Vermont_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Vermont_1_sites$iX <- NA
Vermont_1_sites$iY <- NA
Vermont_1_sites$GY <- NA
Vermont_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Vermont_1_sites$X[i])
  Vermont_1_sites$iX[i] <- which.min(dX$dist)
  Vermont_1_sites$GX[i] <- dX$X[Vermont_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Vermont_1_sites$Y[i])
  Vermont_1_sites$iY[i] <- which.min(dY$dist)
  Vermont_1_sites$GY[i] <- dY$Y[Vermont_1_sites$iY[i]]
  p1 <- c(Vermont_1_sites$X[i],Vermont_1_sites$Y[i])
  p2 <- c(Vermont_1_sites$GX[i],Vermont_1_sites$GY[i])
  Vermont_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Vermont_1_sites$GridID <- paste("Grid_",Vermont_1_sites$GY,"*",Vermont_1_sites$GX,sep="")
# Copy sited results to gsited
gVermont_1_sites <- Vermont_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uVermont_1_sites <- unique(gVermont_1_sites)
# Clean-up sited
Vermont_1_sites$iY <- Vermont_1_sites$iX <- Vermont_1_sites$GY <- Vermont_1_sites$GX <- NULL
# Join dataframes
Vermont_2 <- left_join(Vermont_1, Vermont_1_sites)
# Write dataframes
write_csv(Vermont_2, "Routine 2/Vermont_2.csv")
write_csv(uVermont_1_sites, "Routine 2/Vermont_2_UGrids.csv")
# Clean up dataframes
Vermont_1 <- Vermont_2 <- gVermont_1_sites <- Vermont_1_sites <- uVermont_1_sites <- NULL

# ------------------------------------------------
# Virginia
# ------------------------------------------------
# Read in regional file
Virginia_1 <- read_csv("Routine 1/Virginia_1.csv")
# Get
Virginia_1_sites <- Virginia_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Virginia_1 <- Virginia_1[1:100000,]
nsites <- length(Virginia_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Virginia_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Virginia_1_sites$iX <- NA
Virginia_1_sites$iY <- NA
Virginia_1_sites$GY <- NA
Virginia_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Virginia_1_sites$X[i])
  Virginia_1_sites$iX[i] <- which.min(dX$dist)
  Virginia_1_sites$GX[i] <- dX$X[Virginia_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Virginia_1_sites$Y[i])
  Virginia_1_sites$iY[i] <- which.min(dY$dist)
  Virginia_1_sites$GY[i] <- dY$Y[Virginia_1_sites$iY[i]]
  p1 <- c(Virginia_1_sites$X[i],Virginia_1_sites$Y[i])
  p2 <- c(Virginia_1_sites$GX[i],Virginia_1_sites$GY[i])
  Virginia_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Virginia_1_sites$GridID <- paste("Grid_",Virginia_1_sites$GY,"*",Virginia_1_sites$GX,sep="")
# Copy sited results to gsited
gVirginia_1_sites <- Virginia_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uVirginia_1_sites <- unique(gVirginia_1_sites)
# Clean-up sited
Virginia_1_sites$iY <- Virginia_1_sites$iX <- Virginia_1_sites$GY <- Virginia_1_sites$GX <- NULL
# Join dataframes
Virginia_2 <- left_join(Virginia_1, Virginia_1_sites)
# Write dataframes
write_csv(Virginia_2, "Routine 2/Virginia_2.csv")
write_csv(uVirginia_1_sites, "Routine 2/Virginia_2_UGrids.csv")
# Clean up dataframes
Virginia_1 <- Virginia_2 <- gVirginia_1_sites <- Virginia_1_sites <- uVirginia_1_sites <- NULL

# ------------------------------------------------
# West_Virginia
# ------------------------------------------------
# Read in regional file
West_Virginia_1 <- read_csv("Routine 1/West_Virginia_1.csv")
# Get
West_Virginia_1_sites <- West_Virginia_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# West_Virginia_1 <- West_Virginia_1[1:100000,]
nsites <- length(West_Virginia_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
West_Virginia_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
West_Virginia_1_sites$iX <- NA
West_Virginia_1_sites$iY <- NA
West_Virginia_1_sites$GY <- NA
West_Virginia_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - West_Virginia_1_sites$X[i])
  West_Virginia_1_sites$iX[i] <- which.min(dX$dist)
  West_Virginia_1_sites$GX[i] <- dX$X[West_Virginia_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - West_Virginia_1_sites$Y[i])
  West_Virginia_1_sites$iY[i] <- which.min(dY$dist)
  West_Virginia_1_sites$GY[i] <- dY$Y[West_Virginia_1_sites$iY[i]]
  p1 <- c(West_Virginia_1_sites$X[i],West_Virginia_1_sites$Y[i])
  p2 <- c(West_Virginia_1_sites$GX[i],West_Virginia_1_sites$GY[i])
  West_Virginia_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

West_Virginia_1_sites$GridID <- paste("Grid_",West_Virginia_1_sites$GY,"*",West_Virginia_1_sites$GX,sep="")
# Copy sited results to gsited
gWest_Virginia_1_sites <- West_Virginia_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uWest_Virginia_1_sites <- unique(gWest_Virginia_1_sites)
# Clean-up sited
West_Virginia_1_sites$iY <- West_Virginia_1_sites$iX <- West_Virginia_1_sites$GY <- West_Virginia_1_sites$GX <- NULL
# Join dataframes
West_Virginia_2 <- left_join(West_Virginia_1, West_Virginia_1_sites)
# Write dataframes
write_csv(West_Virginia_2, "Routine 2/West_Virginia_2.csv")
write_csv(uWest_Virginia_1_sites, "Routine 2/West_Virginia_2_UGrids.csv")
# Clean up dataframes
West_Virginia_1 <- West_Virginia_2 <- gWest_Virginia_1_sites <- West_Virginia_1_sites <- uWest_Virginia_1_sites <- NULL

# ------------------------------------------------
# Wisconsin
# ------------------------------------------------
# Read in regional file
Wisconsin_1 <- read_csv("Routine 1/Wisconsin_1.csv")
# Get
Wisconsin_1_sites <- Wisconsin_1 %>% 
  dplyr::select(`LOCALITY ID`, X, Y) %>% 
  unique()
# Wisconsin_1 <- Wisconsin_1[1:100000,]
nsites <- length(Wisconsin_1_sites$`LOCALITY ID`)
# Distance km from nearest grid centroid
Wisconsin_1_sites$Dist_m <- NA
# Compute look-Up coordinates for 200x200km grid
Wisconsin_1_sites$iX <- NA
Wisconsin_1_sites$iY <- NA
Wisconsin_1_sites$GY <- NA
Wisconsin_1_sites$GX <- NA

for (i in 1:nsites) {
  dX$dist <- abs(dX$X - Wisconsin_1_sites$X[i])
  Wisconsin_1_sites$iX[i] <- which.min(dX$dist)
  Wisconsin_1_sites$GX[i] <- dX$X[Wisconsin_1_sites$iX[i]]
  dY$dist <- abs(dY$Y - Wisconsin_1_sites$Y[i])
  Wisconsin_1_sites$iY[i] <- which.min(dY$dist)
  Wisconsin_1_sites$GY[i] <- dY$Y[Wisconsin_1_sites$iY[i]]
  p1 <- c(Wisconsin_1_sites$X[i],Wisconsin_1_sites$Y[i])
  p2 <- c(Wisconsin_1_sites$GX[i],Wisconsin_1_sites$GY[i])
  Wisconsin_1_sites$Dist_m[i]<- sqrt((p1[1] - p2[1]) ^ 2 + (p1[2] - p2[2]) ^ 2)
}

Wisconsin_1_sites$GridID <- paste("Grid_",Wisconsin_1_sites$GY,"*",Wisconsin_1_sites$GX,sep="")
# Copy sited results to gsited
gWisconsin_1_sites <- Wisconsin_1_sites[ ,c("GridID","iX","iY","GY","GX")]
uWisconsin_1_sites <- unique(gWisconsin_1_sites)
# Clean-up sited
Wisconsin_1_sites$iY <- Wisconsin_1_sites$iX <- Wisconsin_1_sites$GY <- Wisconsin_1_sites$GX <- NULL
# Join dataframes
Wisconsin_2 <- left_join(Wisconsin_1, Wisconsin_1_sites)
# Write dataframes
write_csv(Wisconsin_2, "Routine 2/Wisconsin_2.csv")
write_csv(uWisconsin_1_sites, "Routine 2/Wisconsin_2_UGrids.csv")
# Clean up dataframes
Wisconsin_1 <- Wisconsin_2 <- gWisconsin_1_sites <- Wisconsin_1_sites <- uWisconsin_1_sites <- NULL
