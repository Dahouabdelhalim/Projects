total_start <- Sys.time()
# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 7
# ------------------------------------------------
# Load Libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(tidyr)
library(stringr)
library(nls2)
library(nlstools)
library(ggpubr)
library(truncnorm)
library(EnvStats)
library(zoo)
library(lme4)
library(data.table)
library(ecodist)
library(scales)

# Set Working Directory
setwd("~/Data repository")

# ------------------------------------------------
# Start function subroutines
# ------------------------------------------------
# Arrival Date function
ArrivalDate <- function(Slice) {
  
  # Fit nls model
  Output <- nls2(Daily_Site_Occupancy ~ SSlogis(Jd, Asym, xmid, scal), data = Slice)
  
  # Get nls model summary
  Output_summary <- summary(Output)
  
  # Create summary statistic dataframe columns
  Common_Name <- Slice$Common_Name[1]
  GridID <- Slice$GridID[1]
  NJd <- length(Slice$Jd)
  MinJd <- min(Slice$Jd)
  MaxJd <- max(Slice$Jd)
  JdDetect <- nrow(Slice[Slice$Daily_Site_Sightings > 0, ])
  Sighting_Mean <- round(mean(Slice$Daily_Site_Sightings), digits = 3)
  Sighting_Sum <- sum(Slice$Daily_Site_Sightings)
  Site_Mean <- round(mean(Slice$Daily_Site_Total), digits = 3)
  Site_Sum <- sum(Slice$Daily_Site_Total)
  Percent_Site_Sighting_Sums <- round(Sighting_Sum/Site_Sum, digits = 3)
  Observed_Bird <- subset(Slice, Daily_Site_Sightings > 0)
  FAD <- min(Observed_Bird$Jd)
  # Cox and Snell (ML)
  Prsq <- round(with(Slice, 1 - sum(resid(Output)^2)/sum((Daily_Site_Occupancy - mean(Daily_Site_Occupancy))^2)),digits=3)
  
  # Get Mean Arrival Date (MAD) and proportion occupancy on that date
  scal_Estimate <- unname(round(coef(Output)[1], digits = 3))
  MAD <- unname(round(coef(Output)[2], digits = 3))
  MAD_Occ <- unname(round(coef(Output)[3], digits = 3)/2)
  CIpre <- confint2(Output)
  MAD_CI_neg <- round(CIpre[2,1], digits = 3)
  MAD_CI_pos <- round(CIpre[2,2], digits = 3)
  MAD_CI_range <- MAD_CI_pos - MAD_CI_neg
  DF_resid <- df.residual(Output)
  Resid_dev <- deviance(Output)
  Dev_over_DF <- Resid_dev/DF_resid
  P_val_Asym <- Output_summary$parameters[1,4]
  P_val_xmid <- Output_summary$parameters[2,4]
  P_val_scal <- Output_summary$parameters[3,4]
  
  # Plot the MAD, points
  plot(x = Slice$Jd, 
       y = Slice$Daily_Site_Occupancy,
       xlab = "Julian Day",
       ylab = "Daily Site Occupancy",
       ylim = c(0, 1.0),
       cex = 1
  )
  r <- range(Slice$Jd)
  xNew <- seq(r[1],r[2],length.out = 200)
  yNew <- predict(Output, list(Jd = xNew))
  lines(xNew,yNew)
  abline(v = MAD, lty=2)
  rect(xleft = MAD_CI_neg,
       ybottom = -0.1,
       xright = MAD_CI_pos,
       ytop = 1.1,
       border = NA,
       col = rgb(0,0,0,alpha=0.2))
  mtext(paste(Slice$Common_Name, " - ", Slice$Year, " - ", Slice$GridID), cex = 1)
  
  return(    
    data.frame(
      Common_Name,
      GridID,
      NJd,
      MinJd,
      MaxJd,
      JdDetect,
      Sighting_Mean,
      Sighting_Sum,
      Site_Mean,
      Site_Sum,
      Percent_Site_Sighting_Sums,
      FAD,
      Prsq,
      MAD,
      MAD_Occ,
      MAD_CI_neg,
      MAD_CI_pos,
      MAD_CI_range,
      scal_Estimate,
      DF_resid,
      Resid_dev,
      Dev_over_DF,
      P_val_Asym,
      P_val_xmid,
      P_val_scal
    )
  )
}
# addtrunc

# Get random vector with a specific sum
rand_vect <- function(N, M, sd = 1, min, max, pos.only = TRUE) {
  vec <- rnormTrunc(N, M/N, sd, min, max)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}


simcor <- function (x, ymean=0, ysd=1, correlation=0) {
  n <- length(x)
  
  y <- rnorm(n)
  z <- correlation * scale(x)[,1] + sqrt(1 - correlation^2) * 
    scale(resid(lm(y ~ x)))[,1]
  yresult <- ymean + ysd * z
  yresult
}

# ------------------------------------------------
# Finish function subroutines
# ------------------------------------------------

# ---------------------------------
# Read in data
# ---------------------------------
eBird_Filtered <- read_csv("Routine 7/Merged Output/MAD_Output_All_Filter_1.csv")
# # Filter out species that we cannot use
# eBird_Filtered <- eBird_Filtered %>% 
#   filter(!Common_Name %in% c("Common Grackle", "Brown Thrasher", "Brown-headed Cowbird", "Eastern Phoebe", "Gray Catbird", "Hermit Thrush", "Savannah Sparrow", "Yellow-rumped Warbler", "Tree Swallow", "Chimney Swift"))
# Calculate Grid-specific sampling coverage for species
Species_Grid_Total_nyrs <- eBird_Filtered %>% 
  group_by(Common_Name, GridID) %>% 
  summarize(min_Year = min(Year),
            max_Year = max(Year),
            Duration_Years = max_Year - min_Year + 1,
            nyrs = n()) %>% 
  ungroup()
# Join these summary statistics to the data.frame
eBird_Filtered <- left_join(eBird_Filtered, Species_Grid_Total_nyrs)
# Filter for species-grids that have a duration of at least 15 years
eBird_Filtered <- eBird_Filtered %>% 
  filter(Duration_Years >= 15)
# Filter for species-grids that have number of years of at least 78 years
eBird_Filtered <- eBird_Filtered %>%
  filter(nyrs >= 8)
# Look at the maximum number of species that occur in each grid cell
Species_Grid_Max_n <- eBird_Filtered %>% 
  select(Common_Name, GridID) %>% 
  distinct() %>% 
  group_by(GridID) %>% 
  count(name = "GridID_n_spp") %>% 
  ungroup()
# Remove grids that only have one species
eBird_Filtered <- left_join(eBird_Filtered, Species_Grid_Max_n)
eBird_Filtered <- eBird_Filtered %>% 
  filter(GridID_n_spp > 1)

# Calculate Grid-specific sampling coverage for species
Species_Grid_Total_nyrs <- eBird_Filtered %>% 
  group_by(Common_Name, GridID) %>% 
  summarize(min_Year = min(Year),
            max_Year = max(Year),
            Duration_Years = max_Year - min_Year + 1,
            nyrs = n()) %>% 
  ungroup()
# Summarize by year
Yearly_coverage <- eBird_Filtered %>% 
  group_by(Year) %>% 
  summarize(Year_count = n()) %>% 
  ungroup()
# Species year summary
Species_coverage <- eBird_Filtered %>% 
  group_by(Common_Name, Year) %>% 
  summarize(Common_Name_count = n()) %>% 
  ungroup()
mean(Species_coverage$Common_Name_count)
sd(Species_coverage$Common_Name_count)


# Subset for an early data.frame and a late data.frame (2002-2004) and (2017-2019)
eBird_Filtered_Early <- eBird_Filtered %>% 
  filter(Year %in% c(2002:2004))
eBird_Filtered_Late <- eBird_Filtered %>% 
  filter(Year %in% c(2017:2019))



# -------------------------------------------
# Create simulation data
# -------------------------------------------
# Number of grids, species, and years
Years = 18
Grids = 57
Species = 44
# rows
Years*Grids*Species

# Create a Years df
Years_df <- data.frame(Year = c(1:18))

# Create a Grids df
Grids_df <- data.frame(GridID = paste("Grid_", c(1:57), sep = ""))
# Some grids will have more presences than others in general, like if they are close to cities
eBird_Filtered %>% 
  group_by(GridID) %>% 
  summarize(Site_Sum_mean = mean(Site_Sum)) %>% 
  ungroup() %>% 
  ggplot(aes(Site_Sum_mean)) +
  geom_histogram(colour = "black")
# 
Grid_Summary_1 <- eBird_Filtered %>% 
  group_by(GridID) %>% 
  summarize(Site_Sum_mean = mean(Site_Sum)) %>% 
  ungroup()
mean(Grid_Summary_1$Site_Sum_mean); sd(Grid_Summary_1$Site_Sum_mean)
# 
Grid_Summary_2 <- eBird_Filtered %>% 
  group_by(GridID) %>% 
  summarize(Sighting_Sum_mean = mean(Sighting_Sum)) %>% 
  ungroup()
mean(Grid_Summary_2$Sighting_Sum_mean); sd(Grid_Summary_2$Sighting_Sum_mean)
# 
Species_Summary_1 <- eBird_Filtered %>% 
  group_by(Common_Name) %>% 
  summarize(Sighting_Sum_mean = mean(Sighting_Sum)) %>% 
  ungroup()
mean(Species_Summary_1$Sighting_Sum_mean); sd(Species_Summary_1$Sighting_Sum_mean)
# Most have moderate amounts of sites/sightings, but some have many
# Add a multiplier following a lognormal distribution to grids 1-57
test <- eBird_Filtered %>% 
  group_by(GridID) %>% 
  summarize(Site_Sum_mean = mean(Site_Sum)) %>% 
  mutate(Site_Sum_mean = (Site_Sum_mean - mean(Site_Sum_mean)) / sd(Site_Sum_mean)) %>% 
  mutate(Site_Sum_mean = Site_Sum_mean + (min(Site_Sum_mean)*-1 + median(Site_Sum_mean) + 1)) %>%
  ungroup()
set.seed(2)
ggplot(test, aes(Site_Sum_mean)) +
  # geom_histogram(fill = "skyblue", binwidth = 0.05, colour = "black") +
  geom_density(colour = "red") +
  geom_density(aes(x = rlnormTrunc(n = 57, meanlog = 0, sdlog = 0.6, min = 0.5, max = 5)),
               colour = "black")
set.seed(2)
# Add Grid_Mult_Site
Grids_df <- Grids_df %>% 
  mutate(Grid_Mult_Site = rlnormTrunc(n = 57, meanlog = 0, sdlog = 0.6, min = 0.5, max = 5))
ggplot(Grids_df, aes(Grid_Mult_Site)) +
  geom_histogram(binwidth = 0.25,
                 colour = "black")

# Sighting
eBird_Filtered %>% 
  group_by(GridID) %>% 
  summarize(Sighting_Sum_mean = mean(Sighting_Sum)) %>% 
  ungroup() %>% 
  ggplot(aes(Sighting_Sum_mean)) +
  geom_histogram(colour = "black")
# Does it differ from sites?
eBird_Filtered %>% 
  mutate(Site_Sum = (Site_Sum - mean(Site_Sum)) / sd(Site_Sum),
         Sighting_Sum = (Sighting_Sum - mean(Sighting_Sum)) / sd(Sighting_Sum)) %>% 
ggplot() +
  geom_density(aes(Site_Sum), colour = "blue") +
  geom_density(aes(Sighting_Sum), colour = "red")


# test2 <- eBird_Filtered %>% 
#   group_by(GridID) %>% 
#   summarize(Sighting_Sum_mean = mean(Sighting_Sum)) %>% 
#   mutate(Sighting_Sum_mean = (Sighting_Sum_mean - mean(Sighting_Sum_mean)) / sd(Sighting_Sum_mean)) %>% 
#   mutate(Sighting_Sum_mean = Sighting_Sum_mean + (min(Sighting_Sum_mean)*-1 + median(Sighting_Sum_mean) + 1)) %>%
#   ungroup()
# ggplot(test2, aes(Sighting_Sum_mean)) +
#   geom_histogram(fill = "skyblue", binwidth = 0.05, colour = "black") +
#   geom_vline(aes(xintercept = median(Sighting_Sum_mean)))
# 
# set.seed(2)
# ggplot(test2, aes(Sighting_Sum_mean)) +
#   # geom_histogram(fill = "skyblue", binwidth = 0.05, colour = "black") +
#   geom_density(colour = "red") +
#   geom_density(aes(x = rlnormTrunc(n = 57, meanlog = 0, sdlog = 1.1, min = 0.9, max = 4)),
#                colour = "black")
# set.seed(2)
# # Add Grid_Muly_Sighting
# Grids_df <- Grids_df %>% 
#   mutate(Grid_Mult_Sighting = rlnormTrunc(n = 57, meanlog = 0, sdlog = 1.1, min = 0.9, max = 4))
# ggplot(Grids_df, aes(Grid_Mult_Sighting)) +
#   geom_histogram(binwidth = 0.25,
#                  colour = "black")


# -----------------------------------------
# Calculate mean arrival dates
# -----------------------------------------
# Create a Years_Species df
Species_Years <- merge(data.frame(Year = c(1:18)),
                       data.frame(Common_Name = paste("Species_", c(1:44), sep = "")))
# Create a Stationary data.frame
# Species_Years_Arrivals_Stationary
Species_Years_Arrivals_Stationary <- Species_Years
# Get species in general
Species_names <- data.frame(Common_Name = unique(Species_Years$Common_Name))
# Assess the raw data for the overall spread of mean arrival dates (MAD)
MAD_stats <- eBird_Filtered %>%
  group_by(Common_Name) %>%
  summarize(MAD_mean = mean(MAD),
            MAD_sd = sd(MAD),
            MAD_min = min(MAD),
            MAD_max = max(MAD)) %>%
  ungroup()
summary(MAD_stats)
# Plot the data
ggplot(MAD_stats, aes(MAD_mean)) +
  geom_histogram(binwidth = 1)
ggplot(MAD_stats, aes(MAD_sd)) +
  geom_histogram(binwidth = 0.5)
ggplot(MAD_stats, aes(MAD_min)) +
  geom_histogram(binwidth = 2)
ggplot(MAD_stats, aes(MAD_max)) +
  geom_histogram(binwidth = 2)
# Grab the data
MAD_mean_mean <- mean(MAD_stats$MAD_mean)
MAD_mean_sd <- sd(MAD_stats$MAD_mean)
MAD_sd_mean <- mean(MAD_stats$MAD_sd)
MAD_sd_sd <- sd(MAD_stats$MAD_sd)
# Species_names with MAD
set.seed(1)
Species_names <- Species_names %>%
  mutate(MAD_Stationary_mean = round(rnorm(n = n(),
                                           mean = MAD_mean_mean,
                                           sd = MAD_mean_sd)),
         MAD_Stationary_sd = round(rnorm(n = n(),
                                         mean = MAD_sd_mean,
                                         sd = MAD_sd_sd)))

# Join to Species_Years_Arrivals_Stationary
Species_Years_Arrivals_Stationary <- left_join(Species_Years_Arrivals_Stationary, Species_names)
Species_Years_Arrivals_Stationary <- Species_Years_Arrivals_Stationary %>%
  mutate(MAD_Stationary = round(rnorm(n = n(),
                                      mean = MAD_Stationary_mean,
                                      sd = MAD_Stationary_sd)))
# Plot the data
ggplot(Species_Years_Arrivals_Stationary, aes(Year, MAD_Stationary, colour = Common_Name)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(legend.position = "none") +
  facet_wrap(~Common_Name)

# JCreate grid_years data.frames
Grid_Years <- merge(Years_df,
                    Grids_df)
Grid_Years <- Grid_Years %>% 
  mutate(MinimumJd = 80,
         MaximumJd = 185,
         ID = paste(Year, GridID, sep = "_"))

Grid_Years_2 <- Grid_Years %>% 
  rowwise() %>%
  do(data.frame(ID = .$ID, Jd = seq(.$MinimumJd, .$MaximumJd, by = 1)))
Grid_Years_2 <- left_join(Grid_Years_2, Grid_Years)

# Grid_Adj makes each grid slightly more similar to itself over time
Grid_Adj <- data.frame(GridID = unique(Grid_Years_2$GridID))
set.seed(1)
Grid_Adj <- Grid_Adj %>% 
  mutate(Grid_Adj_mean = round(rnorm(n = n(), mean = 0, sd = 1)),
         Grid_Adj_sd = round(rnorm(n = n(), mean = 0, sd = 1)))

Grid_Years_2 <- left_join(Grid_Years_2, Grid_Adj)

# Each Species arrival date with vary by GridID
Grid_Years_3 <- left_join(Grid_Years, distinct(select(Grid_Years_2, Year, GridID, Grid_Adj_mean, Grid_Adj_sd)))
Grid_Years_Species <- left_join(Grid_Years_3, Species_Years_Arrivals_Stationary)
# Create a new Grid_Adj_mean for species, and it is similar to Grid_Adj_mean
set.seed(1)
Grid_Years_Species <- Grid_Years_Species %>%
  mutate(Grid_Adj_mean_spp = round(rnorm(n = n(), mean = 0, sd = 3)))
# Adjust the MAD
set.seed(1)
Grid_Years_Species <- Grid_Years_Species %>%
  mutate(MAD_Stationary_Adj = MAD_Stationary + Grid_Adj_mean_spp)
# Plot the data
ggplot(Grid_Years_Species, aes(factor(Year), MAD_Stationary_Adj, colour = Common_Name)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  facet_wrap(~Common_Name)



# ----------------------------------------------------------
# "Shifting" dataset
# ----------------------------------------------------------
# Create a Shifting data.frame
# Species_Years_Arrivals_Shifting
Species_Years_Arrivals_Shifting <- Species_Years
# Get species in general
Species_names <- data.frame(Common_Name = unique(Species_Years$Common_Name))
# Assess the raw data for the overall spread of mean arrival dates (MAD)
# Early
MAD_stats_early <- eBird_Filtered_Early %>%
  group_by(Common_Name) %>%
  summarize(MAD_early_mean = mean(MAD),
            MAD_early_sd = sd(MAD),
            MAD_early_min = min(MAD),
            MAD_early_max = max(MAD)) %>%
  ungroup()
summary(MAD_stats_early)
# Plot the data
ggplot(MAD_stats_early, aes(MAD_early_mean)) +
  geom_histogram(binwidth = 1)
ggplot(MAD_stats_early, aes(MAD_early_sd)) +
  geom_histogram(binwidth = 0.5)
ggplot(MAD_stats_early, aes(MAD_early_min)) +
  geom_histogram(binwidth = 2)
ggplot(MAD_stats_early, aes(MAD_early_max)) +
  geom_histogram(binwidth = 2)
# Late
MAD_stats_late <- eBird_Filtered_Late %>%
  group_by(Common_Name) %>%
  summarize(MAD_late_mean = mean(MAD),
            MAD_late_sd = sd(MAD),
            MAD_late_min = min(MAD),
            MAD_late_max = max(MAD)) %>%
  ungroup()
summary(MAD_stats_late)
# Plot the data
ggplot(MAD_stats_late, aes(MAD_late_mean)) +
  geom_histogram(binwidth = 1)
ggplot(MAD_stats_late, aes(MAD_late_sd)) +
  geom_histogram(binwidth = 0.5)
ggplot(MAD_stats_late, aes(MAD_late_min)) +
  geom_histogram(binwidth = 2)
ggplot(MAD_stats_late, aes(MAD_late_max)) +
  geom_histogram(binwidth = 2)
# Join the data
MAD_stats_join <- left_join(MAD_stats_early, MAD_stats_late)
MAD_stats_join <- MAD_stats_join %>%
  mutate(MAD_diff_mean = MAD_late_mean - MAD_early_mean,
         MAD_diff_sd = MAD_late_sd - MAD_early_sd)
summary(select(MAD_stats_join, MAD_diff_mean, MAD_diff_sd))
# Plot the data
ggplot(MAD_stats_join, aes(MAD_diff_mean)) +
  geom_histogram(binwidth = 0.5)
ggplot(MAD_stats_join, aes(MAD_diff_sd)) +
  geom_histogram(binwidth = 0.5)




# Grab the data
MAD_early_mean_mean <- mean(MAD_stats_early$MAD_early_mean)
MAD_early_mean_sd <- sd(MAD_stats_early$MAD_early_mean)
MAD_early_sd_mean <- mean(MAD_stats_early$MAD_early_sd)
MAD_early_sd_sd <- sd(MAD_stats_early$MAD_early_sd)
MAD_late_mean_mean <- mean(MAD_stats_late$MAD_late_mean)
MAD_late_mean_sd <- sd(MAD_stats_late$MAD_late_mean)
MAD_late_sd_mean <- mean(MAD_stats_late$MAD_late_sd)
MAD_late_sd_sd <- sd(MAD_stats_late$MAD_late_sd)
# Get diff
MAD_diff_mean_mean <- summary(MAD_stats_join$MAD_diff_mean)[[2]]
MAD_diff_mean_sd <- sd(MAD_stats_join$MAD_diff_mean)
MAD_diff_sd_mean <- mean(MAD_stats_join$MAD_diff_sd)
MAD_diff_sd_sd <- sd(MAD_stats_join$MAD_diff_sd)

# Species_names with MAD
set.seed(1)
# Get early shifting mean
Species_names <- Species_names %>%
  mutate(MAD_early_Shifting_mean = round(rnorm(n = n(),
                                               mean = MAD_early_mean_mean,
                                               sd = MAD_early_mean_sd)),
         MAD_early_Shifting_sd = round(rnorm(n = n(),
                                             mean = MAD_early_sd_mean,
                                             sd = MAD_early_sd_sd)))
# Get late shifting mean
set.seed(1)
Species_names <- Species_names %>%
  mutate(MAD_late_Shifting_mean = MAD_early_Shifting_mean + round(rnorm(n = n(),
                                                                        mean = MAD_diff_mean_mean,
                                                                        sd = MAD_diff_mean_sd)),
         MAD_late_Shifting_sd = MAD_early_Shifting_sd + round(rnorm(n = n(),
                                                                    mean = MAD_diff_sd_mean,
                                                                    sd = MAD_diff_sd_sd)))
# Create data.frame with Year = 1 and Year = 18 with early and late MAD values on each
Species_names_1 <- Species_names %>%
  select(Common_Name) %>%
  mutate(Year = 1)
Species_names_2 <- Species_names %>%
  select(Common_Name) %>%
  mutate(Year = 18)
# Join MAD values
Species_names_1 <- left_join(Species_names_1, select(Species_names, Common_Name, MAD_early_Shifting_mean, MAD_early_Shifting_sd)) %>%
  rename("MAD_Shifting_mean" = "MAD_early_Shifting_mean", "MAD_Shifting_sd" = "MAD_early_Shifting_sd")
Species_names_2 <- left_join(Species_names_2, select(Species_names, Common_Name, MAD_late_Shifting_mean, MAD_late_Shifting_sd)) %>%
  rename("MAD_Shifting_mean" = "MAD_late_Shifting_mean", "MAD_Shifting_sd" = "MAD_late_Shifting_sd")
# Bind the data.frames
Species_names_3 <- rbind(Species_names_1, Species_names_2)


# Join to Species_Years_Arrivals_Shifting
Species_Years_Arrivals_Shifting <- left_join(Species_Years_Arrivals_Shifting, Species_names_3)
Species_Years_Arrivals_Shifting <- left_join(Species_Years_Arrivals_Shifting, Species_names)
# Impute values
Species_Years_Arrivals_Shifting <- Species_Years_Arrivals_Shifting %>%
  mutate(MAD_Shifting = round(na.approx(MAD_Shifting_mean)))
# Add some random variation to the values
set.seed(100)
Species_Years_Arrivals_Shifting <- Species_Years_Arrivals_Shifting %>%
  mutate(MAD_Shifting = MAD_Shifting + round(rnorm(n = n(), mean = 0, sd = 1)))
# # Get the values
# Species_Years_Arrivals_Shifting <- Species_Years_Arrivals_Shifting %>%
#   mutate(MAD_Shifting = round(rnorm(n = n(),
#                                       mean = MAD_Shifting_mean,
#                                       sd = MAD_Shifting_sd)))
# Plot the data
ggplot(Species_Years_Arrivals_Shifting, aes(Year, MAD_Shifting, colour = Common_Name)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(legend.position = "none") +
  facet_wrap(~Common_Name)
ggplot(Species_Years_Arrivals_Shifting, aes(Year, MAD_Shifting)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(legend.position = "none")


# Join to Grid_Years_3
# Each Species arrival date with vary by GridID
# Grid_Years_3 <- left_join(Grid_Years, distinct(select(Grid_Years_2, Year, GridID, Grid_Adj_mean, Grid_Adj_sd)))
Grid_Years_Species <- left_join(Grid_Years_Species, Species_Years_Arrivals_Shifting)
# # Create a new Grid_Adj_mean for species, and it is similar to Grid_Adj_mean
# set.seed(1)
# Grid_Years_Species <- Grid_Years_Species %>%
#   mutate(Grid_Adj_mean_spp = round(rnorm(n = n(), mean = 0, sd = 3)))
# Adjust the MAD
set.seed(1)
Grid_Years_Species <- Grid_Years_Species %>%
  mutate(MAD_Shifting_Adj = MAD_Shifting + Grid_Adj_mean_spp)
# Plot the data
ggplot(Grid_Years_Species, aes(factor(Year), MAD_Shifting_Adj, colour = Common_Name)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  facet_wrap(~Common_Name)
ggplot(Grid_Years_Species, aes(Year, MAD_Shifting_Adj)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(legend.position = "none") +
  facet_wrap(~Common_Name)
ggplot(Grid_Years_Species, aes(Year, MAD_Shifting_Adj)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme(legend.position = "none")



# ----------------------------------------------
# Site_Sum
# ----------------------------------------------
# Insert metrics for Site_Sum
mean((eBird_Filtered_Early$Site_Sum))
sd((eBird_Filtered_Early$Site_Sum))
mean((eBird_Filtered_Late$Site_Sum))
sd((eBird_Filtered_Late$Site_Sum))

# Site_Sum follows a lognormal distribution 
Site_Sum_max = max(eBird_Filtered$Site_Sum)
# mean
Early_Site_Sum_mean = mean(log(eBird_Filtered_Early$Site_Sum))
Late_Site_Sum_mean = mean(log(eBird_Filtered_Late$Site_Sum))
# sd
Early_Site_Sum_sd = sd(log(eBird_Filtered_Early$Site_Sum))
Late_Site_Sum_sd = sd(log(eBird_Filtered_Late$Site_Sum))
# slope
Site_Sum_mean_slope = (Late_Site_Sum_mean - Early_Site_Sum_mean) / Years
Site_Sum_sd_slope = (Late_Site_Sum_sd - Early_Site_Sum_sd) / Years
# Execute calculation (multiplying mean by grid mult)
set.seed(1)
Grid_Years <- Grid_Years %>%
  mutate(Site_Sum = round(rlnormTrunc(n = n(), meanlog = Early_Site_Sum_mean + (Year*Site_Sum_mean_slope) * Grid_Mult_Site, sdlog = Early_Site_Sum_sd + (Year*Site_Sum_sd_slope), min = 10, max = Site_Sum_max)))
# Plot and subset data to ensure it worked correctly
ggplot(Grid_Years, aes(Site_Sum, colour = Year, group = Year)) +
  geom_density() +
  scale_x_continuous(trans = "log")
max(Grid_Years$Site_Sum)
ggplot() +
  geom_density(data = Grid_Years, aes(Site_Sum, colour = Year, group = Year))
sd(Grid_Years[Grid_Years$Year == 18, ]$Site_Sum)
sd(Grid_Years[Grid_Years$Year == 1, ]$Site_Sum)
sd(log(Grid_Years[Grid_Years$Year == 18, ]$Site_Sum))
sd(log(Grid_Years[Grid_Years$Year == 1, ]$Site_Sum))
# Compare to the raw data
ggplot(Grid_Years, aes(Site_Sum, colour = Year, group = Year)) +
  geom_density() +
  scale_x_continuous(trans = "log")
ggplot(eBird_Filtered, aes(Site_Sum, colour = Year, group = Year)) +
  geom_density() +
  scale_x_continuous(trans = "log")

ggplot(Grid_Years, aes(Site_Sum, colour = Year, group = Year)) +
  geom_density()
ggplot(eBird_Filtered, aes(Site_Sum, colour = Year, group = Year)) +
  geom_density()


# mean
Early_Site_Sum_mean_check = mean((eBird_Filtered_Early$Sighting_Sum))
Late_Site_Sum_mean_check = mean((eBird_Filtered_Late$Sighting_Sum))
# sd
Early_Site_Sum_sd_check = sd((eBird_Filtered_Early$Sighting_Sum))
Late_Site_Sum_sd_check = sd((eBird_Filtered_Late$Sighting_Sum))

# Create a data.frame of Species
Species_df <- data.frame(Common_Name = paste("Species_", c(1:44), sep = ""))
Species_df$Common_Name <- factor(Species_df$Common_Name, levels = Species_df$Common_Name, ordered = TRUE)
# See how the overall abundance of different species is different in the raw data
ggplot(eBird_Filtered, aes(Sighting_Sum, group = Common_Name, colour = Common_Name)) +
  geom_density() +
  theme(legend.position = "none")
ggplot(eBird_Filtered, aes(Sighting_Sum, group = Common_Name, colour = Common_Name)) +
  geom_density() +
  scale_x_continuous(trans = "log")
ggplot(eBird_Filtered, aes(Year, Sighting_Sum, group = Common_Name, colour = Common_Name)) +
  geom_point() +
  scale_y_continuous(trans = "log") +
  facet_wrap(~Common_Name)
# Sighting_Sum increases in a linear fashion for each species
# There is variation among species' Sighting_Sum values
eBird_Filtered %>% 
  group_by(Common_Name) %>% 
  summarize(Sighting_Sum_mean_spp = mean(Sighting_Sum)) %>% 
  ungroup() %>% 
  arrange(desc(Sighting_Sum_mean_spp)) %>% 
  base::print(n = 44)
# Plot the values
eBird_Filtered %>% 
  group_by(Common_Name) %>% 
  summarize(Sighting_Sum_mean_spp = mean(Sighting_Sum)) %>% 
  ungroup() %>% 
  arrange(desc(Sighting_Sum_mean_spp)) %>% 
  base::print(n = 44) %>% 
  ggplot(aes(Sighting_Sum_mean_spp)) +
  geom_histogram(binwidth = 50, colour = "black")
# Min ~ 200, Max ~ 1200
# Look at the early vs. late periods
# Early
eBird_Filtered %>% 
  filter(Year %in% c(2002:2004)) %>% 
  group_by(Common_Name) %>% 
  summarize(Sighting_Sum_mean_spp = mean(Sighting_Sum)) %>% 
  ungroup() %>% 
  arrange(desc(Sighting_Sum_mean_spp)) %>% 
  base::print(n = 44)
# Plot the values
eBird_Filtered %>% 
  filter(Year %in% c(2002:2004)) %>% 
  group_by(Common_Name) %>% 
  summarize(Sighting_Sum_mean_spp = mean(Sighting_Sum)) %>% 
  ungroup() %>% 
  arrange(desc(Sighting_Sum_mean_spp)) %>% 
  base::print(n = 44) %>% 
  ggplot(aes(Sighting_Sum_mean_spp)) +
  geom_histogram(binwidth = 25, colour = "black")

# Try plotting density plots
set.seed(1)
eBird_Filtered %>% 
  filter(Year %in% c(2002:2004)) %>% 
  group_by(Common_Name) %>% 
  summarize(Sighting_Sum_mean_spp = mean(Sighting_Sum)) %>% 
  ungroup() %>% 
  arrange(desc(Sighting_Sum_mean_spp)) %>% 
  base::print(n = 44) %>% 
  ggplot(aes(Sighting_Sum_mean_spp)) +
  geom_density( colour = "black") +
  geom_density(aes(x = rlnormTrunc(n = 44, meanlog = 1.5, sdlog = 2, min = 0.4, max = 2)*100), colour = "red")


# Min ~ 40, Max ~ 250
# Late
eBird_Filtered %>% 
  filter(Year %in% c(2017:2019)) %>% 
  group_by(Common_Name) %>% 
  summarize(Sighting_Sum_mean_spp = mean(Sighting_Sum)) %>% 
  ungroup() %>% 
  arrange(desc(Sighting_Sum_mean_spp)) %>% 
  base::print(n = 44)
# Plot the values
eBird_Filtered %>% 
  filter(Year %in% c(2017:2019)) %>% 
  group_by(Common_Name) %>% 
  summarize(Sighting_Sum_mean_spp = mean(Sighting_Sum)) %>% 
  ungroup() %>% 
  arrange(desc(Sighting_Sum_mean_spp)) %>% 
  base::print(n = 44) %>% 
  ggplot(aes(Sighting_Sum_mean_spp)) +
  geom_histogram(binwidth = 100, colour = "black")
# Min ~ 400, Max ~ 3500
# If you have a species that starts at 50, and ends at 350 over 18 years, that is an increase of ~20 sightings per year
(400-40)/18
# If you have a species that starts at 250, and ends at 3500 over 18 year, thats an increase of ~175 sightings per year
(3500-363)/18

# So, we can have species that are low in abundance, those that are high in abundance, and those that fall in between
# Increment between species
# Early (we are going to have 29 species in one group, and 15 in another (departed species))
(250-40)/29
(250-40)/15


# For Species 1-29, it can go 15, 25, 35.... 115 per year increase, starting at 50, 60 (group that stays), and starting at 50, 70, 90... 250 (group that departs)
# The distribution is lognormal
set.seed(1)
summary(rlnormTrunc(n = 44, meanlog = 1.5, sdlog = 2, min = 0.4, max = 2)*100)
set.seed(1)
Early_Sighting_Sum <- c(sort(round(rlnormTrunc(n = 29, meanlog = 1.5, sdlog = 2, min = 0.4, max = 2)*100)), sort(round(rlnormTrunc(n = 15, meanlog = 1.5, sdlog = 2, min = 0.4, max = 2)*100)))

# Bind to data.frame
Species_df <- cbind(Species_df, Early_Sighting_Sum)
# Create vector of per-year increases
Sightings_Yearly_Increase <- c(seq(15, 155, 5), seq(15, 160, 10))
# Bind to data.frame
Species_df <- cbind(Species_df, Sightings_Yearly_Increase)
# Calculate the Late Sightings sum
Species_df <- Species_df %>% 
  mutate(Late_Sighting_Sum = Early_Sighting_Sum + (Years*Sightings_Yearly_Increase))

# Create a new data.frame that takes these values and stretches them across years
Species_df_2 <- Species_df %>% 
  select(Common_Name, Early_Sighting_Sum) %>% 
  mutate(Year = 1) %>% 
  rename("Sighting_Sum" = "Early_Sighting_Sum")
Species_df_3 <- Species_df %>% 
  select(Common_Name, Late_Sighting_Sum) %>% 
  mutate(Year = 18) %>% 
  rename("Sighting_Sum" = "Late_Sighting_Sum")
Species_df_4 <- rbind(Species_df_2, Species_df_3)

# Join the Sightings data to the Species_Years df
Species_df_5 <- left_join(Species_Years, Species_df_4)
# Fill in the gaps
Species_df_5 <- Species_df_5 %>% 
  group_by(Common_Name) %>% 
  mutate(Sighting_Sum = round(na.approx(Sighting_Sum))) %>% 
  ungroup()
# Denote birds as departed or not (Species 30-44 = Departed)
Species_df_5 <- Species_df_5 %>% 
  mutate(Departed_spp = ifelse(Common_Name %in% paste("Species_", c(30:44), sep = ""), TRUE, FALSE))
# Join Species_df so that we have this information for later on
Species_df_5 <- left_join(Species_df_5, Species_df)



# Insert metrics for MinJd
mean(eBird_Filtered_Early$MinJd)
mean(eBird_Filtered_Late$MinJd)
sd(eBird_Filtered_Early$MinJd)
sd(eBird_Filtered_Late$MinJd)
eBird_Filtered_Early %>% 
  count(MinJd == 80)
eBird_Filtered_Late %>% 
  count(MinJd == 80)
eBird_Filtered_Early %>% 
  filter(MinJd != 80) %>% 
  summarize(mean(MinJd),
            sd(MinJd))
eBird_Filtered_Late %>% 
  filter(MinJd != 80) %>% 
  summarize(mean(MinJd),
            sd(MinJd))
ggplot(eBird_Filtered_Early, aes(MinJd, colour = Year, group = Year)) +
  geom_density()
ggplot(eBird_Filtered_Late, aes(MinJd, colour = Year, group = Year)) +
  geom_density()
# Set it always at 80
Species_df_5 <- Species_df_5 %>% 
  mutate(MinJd = 80)

# Insert metrics for MaxJd
# Check distributions
ggplot(eBird_Filtered_Early, aes(MaxJd)) +
  geom_density()
ggplot(eBird_Filtered_Late, aes(MaxJd)) +
  geom_density()
# How many species typically have MaxJd values that are less than 180?
eBird_Filtered_Early %>% 
  filter(MaxJd < 160) %>% 
  count(Common_Name) %>% 
  filter(n > 0) %>% 
  base::print(n = 44)
eBird_Filtered_Late %>% 
  filter(MaxJd < 160) %>% 
  count(Common_Name) %>% 
  filter(n > 0) %>% 
  base::print(n = 44)
# About half of the species in the late period, most species in the early period.
# To be consistent throughout, set 2/3 of species as consistently remaining, and 1/3 as departing
# Species 1-29 will have 180 as their MaxJd, and their arrivals will follow a pnorm dist
# Species 30-44 will have ~142-160 as their MaxJd, and their arrivals will follow a dnorm
# 
# Set to 185 to account for the cutoff mechanism
# Try 185 for all species. See how cut-off mechanism works
set.seed(1)
MaxJd_mean = 5 ; MaxJd_sd = 2
Species_df_5 <- Species_df_5 %>%
  mutate(MaxJd = 185)
# Species_df_5 <- Species_df_5 %>% 
#   mutate(MaxJd = case_when(Common_Name %in% paste("Species_", c(1:20), sep = "") ~ 185,
#                            Common_Name %in% "Species_21" ~ 144,
#                            Common_Name %in% "Species_22" ~ 146,
#                            Common_Name %in% "Species_23" ~ 148,
#                            Common_Name %in% "Species_24" ~ 150,
#                            Common_Name %in% "Species_25" ~ 152,
#                            Common_Name %in% "Species_26" ~ 154,
#                            Common_Name %in% "Species_27" ~ 156,
#                            Common_Name %in% "Species_28" ~ 158,
#                            Common_Name %in% "Species_29" ~ 160))

  
# Plot and subset data to ensure it worked correctly
ggplot(Species_df_5, aes(MaxJd)) +
  geom_density()
sd(Species_df_5[Species_df_5$Year == 18, ]$MaxJd)
sd(Species_df_5[Species_df_5$Year == 1, ]$MaxJd)

# Insert NJd
Species_df_5 <- Species_df_5 %>% 
  mutate(NJd = MaxJd - MinJd + 1)
ggplot(Species_df_5, aes(NJd)) +
  geom_density()
sd(Species_df_5[Species_df_5$Year == 18, ]$NJd)
sd(Species_df_5[Species_df_5$Year == 1, ]$NJd)

# Check out grid metrics and species metrics
testing_2 <- eBird_Filtered %>% 
  group_by(GridID) %>% 
  summarise(mean_Site_Sum = mean(Site_Sum),
            mean_Sighting_Sum = mean(Sighting_Sum)) %>% 
  ungroup()
mean(testing_2$mean_Site_Sum)
mean(testing_2$mean_Sighting_Sum)
sd(testing_2$mean_Site_Sum)
sd(testing_2$mean_Sighting_Sum)

# Check out grid metrics and species metrics
testing <- eBird_Filtered %>% 
  group_by(Common_Name) %>% 
  summarise(mean_Sighting_Sum = mean(Sighting_Sum)) %>% 
  ungroup()
mean(testing$mean_Sighting_Sum)
sd(testing$mean_Sighting_Sum)


# Join Species_df_5 to Grid_Years_Species
Species_df_5 <- left_join(Species_df_5, Grid_Years_Species)
# Insert metrics for JdDetect
summary(eBird_Filtered_Early$JdDetect)
summary(eBird_Filtered_Late$JdDetect)
sd(eBird_Filtered_Early$JdDetect)
sd(eBird_Filtered_Late$JdDetect)
# Check distributions
ggplot(eBird_Filtered_Early, aes(JdDetect)) +
  geom_density()
ggplot(eBird_Filtered_Late, aes(JdDetect)) +
  geom_density()
# These are two different distributions
ggplot(eBird_Filtered, aes(factor(Year), JdDetect, colour = Year)) +
  geom_boxplot()
ggplot(eBird_Filtered, aes(Year, JdDetect, colour = Year)) +
  geom_point() +
  geom_smooth() +
  geom_smooth(method = "lm", colour = "black")
# It appears as through JdDetect increased linearly from 2002-2012
ggplot(eBird_Filtered, aes(JdDetect, NJd)) +
  geom_point() +
  geom_abline()
# JdDetect follows a normal distribution in later years, and is log normal in early years both are truncated
# mean
Early_JdDetect_mean = mean(eBird_Filtered_Early$JdDetect)
# Early_1_JdDetect_mean = mean(log(eBird_Filtered_Early$JdDetect))
# Early_2_JdDetect_mean = mean(log(eBird_Filtered_Late$JdDetect))
Late_JdDetect_mean = mean(eBird_Filtered_Late$JdDetect)
# sd
Early_JdDetect_sd = sd(eBird_Filtered_Early$JdDetect)

JdDetect_sd <- eBird_Filtered %>%
  group_by(GridID, Common_Name) %>%
  summarize(JdDetect_sd = sd(JdDetect)) %>%
  ungroup() %>% 
  summarize(JdDetect_sd_mean = mean(JdDetect_sd, na.rm = TRUE),
            JdDetect_sd_sd = sd(JdDetect_sd, na.rm = TRUE)) %>%
  ungroup()

# Early_JdDetect_sd <- eBird_Filtered_Early %>%
#   group_by(GridID) %>%
#   summarize(JdDetect_sd = sd(JdDetect)) %>%
#   summarize(JdDetect_sd_mean = mean(JdDetect_sd, na.rm = TRUE),
#             JdDetect_sd_sd = sd(JdDetect_sd, na.rm = TRUE)) %>%
#   ungroup()
# Early_1_JdDetect_sd <- eBird_Filtered_Early %>% 
#   group_by(GridID) %>% 
#   summarize(JdDetect_sd = sd(log(JdDetect))) %>% 
#   summarize(JdDetect_sd_mean = mean(JdDetect_sd, na.rm = TRUE),
#             JdDetect_sd_sd = sd(JdDetect_sd, na.rm = TRUE)) %>% 
#   ungroup()
# Early_2_JdDetect_sd <- eBird_Filtered_Late %>% 
#   group_by(GridID) %>% 
#   summarize(JdDetect_sd = sd(log(JdDetect))) %>% 
#   summarize(JdDetect_sd_mean = mean(JdDetect_sd),
#             JdDetect_sd_sd = sd(JdDetect_sd)) %>% 
#   ungroup()
# Late_JdDetect_sd <- eBird_Filtered_Late %>% 
#   group_by(GridID) %>% 
#   summarize(JdDetect_sd = sd(JdDetect)) %>% 
#   summarize(JdDetect_sd_mean = mean(JdDetect_sd),
#             JdDetect_sd_sd = sd(JdDetect_sd)) %>% 
#   ungroup()
Late_JdDetect_sd = sd(eBird_Filtered_Late$JdDetect)
# Sd stays relatively the same
JdDetect_sd = sd(eBird_Filtered$JdDetect)
# # sd
# Early_1_JdDetect_sd = sd(log(eBird_Filtered_Early$JdDetect))
# Early_2_JdDetect_sd = sd(log(eBird_Filtered_Late$JdDetect))
# Late_JdDetect_sd = sd(eBird_Filtered_Late$JdDetect)

# How does JdDetect relate to departed birds?
eBird_Filtered %>%
  group_by(MaxJd < 180) %>%
  summarize(JdDetect_mean = mean(JdDetect),
            JdDetect_sd = sd(JdDetect)) %>%
  ungroup() %>%
  base::print(n = 44)
eBird_Filtered %>% 
  count(MaxJd < 180)
eBird_Filtered <- eBird_Filtered %>% 
  mutate(departed = if_else(MaxJd < 180, TRUE, FALSE))
eBird_Filtered$departed <- as.factor(eBird_Filtered$departed)

# How does JdDetect relate to species?
eBird_Filtered %>%
  group_by(Common_Name) %>%
  summarize(JdDetect_mean = mean(JdDetect),
            JdDetect_sd = sd(JdDetect)) %>%
  ungroup() %>%
  base::print(n = 44)
# How does JdDetect relate to Sighting_Sum
ggplot(eBird_Filtered, aes(Sighting_Sum, JdDetect)) +
  geom_point() +
  scale_x_log10()
# Segment by quantiles
# Empirical data
quantile(eBird_Filtered$Sighting_Sum)
quantile(Species_df_5$Sighting_Sum)
eBird_Filtered <- eBird_Filtered %>% 
  mutate(sighting_sum_class = case_when(
    Sighting_Sum < 150 ~ "Class 1",
    Sighting_Sum >= 150 & Sighting_Sum < 361 ~ "Class 2",
    Sighting_Sum >= 361 & Sighting_Sum < 909 ~ "Class 3",
    Sighting_Sum >= 909 ~ "Class 4" 
  ))
eBird_Filtered$sighting_sum_class <- factor(eBird_Filtered$sighting_sum_class, levels = c("Class 1", "Class 2", "Class 3", "Class 4"), ordered = TRUE)
sighting_sum_class_JdDetect <- eBird_Filtered %>% 
  group_by(sighting_sum_class) %>% 
  summarize(JdDetect_mean = mean(JdDetect),
            JdDetect_sd = sd(JdDetect)) %>% 
  ungroup()
sighting_sum_class_JdDetect
# Simulated data
Species_df_5 <- Species_df_5 %>% 
  mutate(sighting_sum_class = case_when(
    Sighting_Sum < 323.75 ~ "Class 1",
    Sighting_Sum >= 323.75 & Sighting_Sum < 690.50 ~ "Class 2",
    Sighting_Sum >= 690.50 & Sighting_Sum < 1288.00 ~ "Class 3",
    Sighting_Sum >= 1288.00 ~ "Class 4" 
  ))
Species_df_5$sighting_sum_class <- factor(Species_df_5$sighting_sum_class, levels = c("Class 1", "Class 2", "Class 3", "Class 4"), ordered = TRUE)



sighting_detect_lm <- lm(JdDetect ~ log(Sighting_Sum), eBird_Filtered)
summary(sighting_detect_lm)

# What about year?
ggplot(eBird_Filtered, aes(Year, JdDetect)) +
  geom_point() +
  scale_x_log10()
detect_year_lm <- lm(JdDetect ~ Year, eBird_Filtered)
summary(detect_year_lm)

# Departed?
ggplot(eBird_Filtered, aes(departed, JdDetect)) +
  geom_boxplot()
deparated_detect_lm <- lm(JdDetect ~ departed, eBird_Filtered)
summary(deparated_detect_lm)

# All
sighting_year_departed_detect_lm <- lm(JdDetect ~ Year + log(Sighting_Sum) + departed, eBird_Filtered)
summary(sighting_year_departed_detect_lm)

# All with MAD
MAD_sighting_year_departed_detect_lm <- lm(JdDetect ~ Year + log(Sighting_Sum) + departed + MAD, eBird_Filtered)
summary(MAD_sighting_year_departed_detect_lm)
# Residual standard error:
RSE = 7.786

# Predict values
eBird_Filtered$JdDetect_pred <- predict(MAD_sighting_year_departed_detect_lm)
# calculate RMSE
hydroGOF::rmse(eBird_Filtered$JdDetect, eBird_Filtered$JdDetect_pred)
# 
# test_lm <- lm(log(JdDetect) ~ Year + log(Sighting_Sum) + departed, eBird_Filtered)
# summary(test_lm)

# Minimum increase
minimum_slope <- (15 - 10) / Years
maximum_slope <- (106 - 97) / Years


# slope
JdDetect_mean_slope = (Late_JdDetect_mean - Early_JdDetect_mean) / Years
JdDetect_sd_slope = (Late_JdDetect_sd - Early_JdDetect_sd) / Years
# JdDetect_mean_slope_Early = (Early_2_JdDetect_mean - Early_1_JdDetect_mean) / (2012-2002+1)
# JdDetect_sd_slope_Early = (Early_2_JdDetect_sd$JdDetect_sd_mean - Early_1_JdDetect_sd$JdDetect_sd_mean) / (2012-2002+1)
# No slope for Late
# # Execute calculation
# set.seed(1)
# Species_df_5$JdDetect1 <- sighting_detect_lm$coefficients[[1]] + log(Species_df_5$Sighting_Sum) * sighting_detect_lm$coefficients[[2]]
# Species_df_5$JdDetect2 <- detect_year_lm$coefficients[[1]] + (Species_df_5$Year + 2001) * detect_year_lm$coefficients[[2]]
# Species_df_5$JdDetect <- (Species_df_5$JdDetect1 + Species_df_5$JdDetect2) / 2
# Calculate using regression equations

Species_df_5$JdDetect_Stationary <- MAD_sighting_year_departed_detect_lm$coefficients[[1]] + (Species_df_5$Year + 2001) * MAD_sighting_year_departed_detect_lm$coefficients[[2]] + log(Species_df_5$Sighting_Sum) * MAD_sighting_year_departed_detect_lm$coefficients[[3]] + ifelse(Species_df_5$Departed_spp == TRUE, MAD_sighting_year_departed_detect_lm$coefficients[[4]], 0) + Species_df_5$MAD_Stationary * MAD_sighting_year_departed_detect_lm$coefficients[[5]]

Species_df_5$JdDetect_Shifting <- MAD_sighting_year_departed_detect_lm$coefficients[[1]] + (Species_df_5$Year + 2001) * MAD_sighting_year_departed_detect_lm$coefficients[[2]] + log(Species_df_5$Sighting_Sum) * MAD_sighting_year_departed_detect_lm$coefficients[[3]] + ifelse(Species_df_5$Departed_spp == TRUE, MAD_sighting_year_departed_detect_lm$coefficients[[4]], 0) + Species_df_5$MAD_Shifting * MAD_sighting_year_departed_detect_lm$coefficients[[5]]

# # Add random variation
# set.seed(1)
# Species_df_5 <- Species_df_5 %>%
#   mutate(JdDetect = round(JdDetect + rnormTrunc(n = n(),
#                                           mean = 0,
#                                           sd = JdDetect_sd$JdDetect_sd_mean,
#                                           min = 10 - JdDetect,
#                                           max = 106 - JdDetect)))
# Species_df_5 <- Species_df_5 %>% 
#   mutate(JdDetect = if_else(
#     Year <= 11,
#     round(rlnormTrunc(n = n(), 
#                       meanlog = Early_1_JdDetect_mean + (Year-1 * JdDetect_mean_slope_Early), # Minus 1 because first year 
#                       sdlog = Early_1_JdDetect_sd$JdDetect_sd_mean + rlnorm(n = n(), meanlog = 0, sdlog = Early_1_JdDetect_sd$JdDetect_sd_sd) + (Year-1 * JdDetect_sd_slope_Early),
#                       min = 10,
#                       max = 97)), # 97 is the maximum for the early period (see)
#     round(rtruncnorm(n = n(), 
#                       mean = Late_JdDetect_mean,
#                       sd = Late_JdDetect_sd$JdDetect_sd_mean + rnorm(n = n(), mean = 0, sd = Late_JdDetect_sd$JdDetect_sd_sd),
#                       a = 15, # 15 is the minimum for the late period (see summary)
#                       b = NJd))))

# Plot and subset data to ensure it worked correctly
ggplot(Species_df_5, aes(JdDetect_Stationary, colour = Year, group = Year)) +
  geom_density()
ggplot(Species_df_5, aes(JdDetect_Shifting, colour = Year, group = Year)) +
  geom_density()
sd(Species_df_5[Species_df_5$Year == 18, ]$JdDetect_Stationary)
sd(Species_df_5[Species_df_5$Year == 1, ]$JdDetect_Stationary)
# # Plot to make sure JdDetect aligns with NJd
# ggplot(Species_df_5, aes(JdDetect, NJd, colour = Year)) +
#   geom_point(alpha = 0.25) +
#   geom_abline()
# Species_df_5 %>% 
#   filter(NJd < 100) %>% 
# ggplot(aes(JdDetect, NJd, colour = Year)) +
#   geom_point(alpha = 0.25) +
#   geom_abline()
# Plot against raw data
ggplot(Species_df_5, aes(JdDetect_Stationary, colour = Year, group = Year)) +
  geom_density() +
  scale_x_continuous(limits = c(10, 106))
ggplot(Species_df_5, aes(JdDetect_Shifting, colour = Year, group = Year)) +
  geom_density() +
  scale_x_continuous(limits = c(10, 106))
ggplot(eBird_Filtered, aes(JdDetect, colour = Year, group = Year)) +
  geom_density() +
  scale_x_continuous(limits = c(10, 106))
# The distribution is slightly different, but much closer than before
# Look at the relationship between JdDetect and Year
ggplot(Species_df_5, aes(factor(Year), JdDetect_Stationary, colour = Year)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(10, 106))
ggplot(Species_df_5, aes(factor(Year), JdDetect_Shifting, colour = Year)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(10, 106))
ggplot(eBird_Filtered, aes(factor(Year), JdDetect, colour = Year)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(10, 106))

# Join Species_df_5 to Grid_Years
head(Grid_Years)
head(Species_df_5)
# Join
Blocks <- left_join(Grid_Years, Species_df_5)
Blocks <- as_tibble(Blocks)

# Correction for bias: Stationary
set.seed(1)
Blocks$JdDetect_Stationary_cor <- rnorm(n = length(Blocks$JdDetect_Stationary), mean = 0, sd = RSE/2)
# Correct
Blocks$JdDetect_Stationary <- Blocks$JdDetect_Stationary + Blocks$JdDetect_Stationary_cor
Blocks$JdDetect_Stationary_cor <- NULL


# Correction for bias: Shifting
set.seed(2)
Blocks$JdDetect_Shifting_cor <- rnorm(n = length(Blocks$JdDetect_Shifting), mean = 0, sd = RSE/2)
# Correct
Blocks$JdDetect_Shifting <- Blocks$JdDetect_Shifting + Blocks$JdDetect_Shifting_cor
Blocks$JdDetect_Shifting_cor <- NULL



# Plot against raw data
ggplot(Blocks, aes(JdDetect_Stationary, colour = Year, group = Year)) +
  geom_density() +
  scale_x_continuous(limits = c(10, 106))
ggplot(Blocks, aes(JdDetect_Shifting, colour = Year, group = Year)) +
  geom_density() +
  scale_x_continuous(limits = c(10, 106))
ggplot(eBird_Filtered, aes(JdDetect, colour = Year, group = Year)) +
  geom_density() +
  scale_x_continuous(limits = c(10, 106))
# The distribution is slightly different, but much closer than before
# Look at the relationship between JdDetect and Year
ggplot(Blocks, aes(factor(Year), JdDetect_Stationary, colour = Year)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(10, 106))
ggplot(Blocks, aes(factor(Year), JdDetect_Shifting, colour = Year)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(10, 106))
ggplot(eBird_Filtered, aes(factor(Year), JdDetect, colour = Year)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(10, 106))
# # Still missing some higher outliers, as well as some low ones
# # Calculate the median value of each year, as well as the difference between the max/min and 106/10
# yearly_median_JdDetect <- Blocks %>% 
#   group_by(Year) %>% 
#   summarize(median_JdDetect = median(JdDetect),
#             max_106_diff = 106 - max(JdDetect),
#             min_10_diff = min(JdDetect) - 10) %>% 
#   ungroup()
# # Join to the data
# Blocks <- left_join(Blocks, yearly_median_JdDetect)
# Normalize the data to these maximum and minimum values
minimum_slope <- (15 - 10) / Years
maximum_slope <- (106 - 97) / Years
# Add minimum and maximum JdDetect based on the year
Blocks <- Blocks %>% 
  mutate(minimum_JdDetect = round(10 + (Year * minimum_slope)),
         maximum_JdDetect = round(97 + (Year * maximum_slope)))

summary(select(Blocks, JdDetect_Stationary, JdDetect_Shifting))

# Adjust bottom
Blocks$JdDetect_Stationary <- if_else(Blocks$JdDetect_Stationary < Blocks$minimum_JdDetect, Blocks$minimum_JdDetect + 1, Blocks$JdDetect_Stationary)
# Adjust top
Blocks$JdDetect_Stationary <- if_else(Blocks$JdDetect_Stationary > Blocks$maximum_JdDetect, Blocks$maximum_JdDetect - 1, Blocks$JdDetect_Stationary)
# Adjust bottom
Blocks$JdDetect_Shifting <- if_else(Blocks$JdDetect_Shifting < Blocks$minimum_JdDetect, Blocks$minimum_JdDetect + 1, Blocks$JdDetect_Shifting)
# Adjust top
Blocks$JdDetect_Shifting <- if_else(Blocks$JdDetect_Shifting > Blocks$maximum_JdDetect, Blocks$maximum_JdDetect - 1, Blocks$JdDetect_Shifting)

summary(select(Blocks, JdDetect_Stationary, JdDetect_Shifting))

# # Add median JdDetect_Stationary
# Blocks <- Blocks %>% 
#   group_by(Year) %>% 
#   mutate(median_JdDetect_Stationary = median(JdDetect_Stationary),
#          maximum_meas_Stationary = max(JdDetect_Stationary),
#          minimum_meas_Stationary = min(JdDetect_Stationary)) %>% 
#   ungroup()
# # Add median JdDetect_Shifting
# Blocks <- Blocks %>% 
#   group_by(Year) %>% 
#   mutate(median_JdDetect_Shifting = median(JdDetect_Shifting),
#          maximum_meas_Shifting = max(JdDetect_Shifting),
#          minimum_meas_Shifting = min(JdDetect_Shifting)) %>% 
#   ungroup()
# # Add specific maximum and minimum diff values
# Blocks <- Blocks %>% 
#   group_by(Year) %>%
#   mutate(below_median_div_Stationary = maximum_meas_Stationary / minimum_JdDetect,
#          above_median_mult_Stationary = maximum_JdDetect / minimum_meas_Stationary) %>% 
#   ungroup()
# Blocks <- Blocks %>% 
#   group_by(Year) %>%
#   mutate(below_median_div_Shifting = maximum_meas_Shifting / minimum_JdDetect,
#          above_median_mult_Shifting = maximum_JdDetect / minimum_meas_Shifting) %>% 
#   ungroup()
# # Add adjusted below_median_div/above_median_mult depending on how far the value is from the median value
# # Get below and above obs
# Blocks_below_Stationary <- Blocks %>% 
#   filter(JdDetect_Stationary <= median_JdDetect_Stationary)
# Blocks_above_Stationary <- Blocks %>% 
#   filter(JdDetect_Stationary >= median_JdDetect_Stationary)
# Blocks_below_Shifting <- Blocks %>% 
#   filter(JdDetect_Shifting <= median_JdDetect_Shifting)
# Blocks_above_Shifting <- Blocks %>% 
#   filter(JdDetect_Shifting >= median_JdDetect_Shifting)
# 
# Blocks_orig <- Blocks
# # Arrange Df
# Blocks_below_Stationary <- Blocks_below_Stationary %>% 
#   arrange(Year, JdDetect_Stationary)
# Blocks_above_Stationary <- Blocks_above_Stationary %>% 
#   arrange(Year, JdDetect_Stationary)
# Blocks_below_Shifting <- Blocks_below_Shifting %>% 
#   arrange(Year, JdDetect_Shifting)
# Blocks_above_Shifting <- Blocks_above_Shifting %>% 
#   arrange(Year, JdDetect_Shifting)
# # Add adjustment column for lowest/higest value and medians
# Blocks_below_Stationary <- Blocks_below_Stationary %>% 
#   mutate(JdDetect_Stationary_adj = ifelse(JdDetect_Stationary == minimum_meas_Stationary, below_median_div_Stationary,
#                                ifelse(JdDetect_Stationary == median_JdDetect_Stationary, 1, 
#                                       NA)))
# Blocks_above_Stationary <- Blocks_above_Stationary %>% 
#   mutate(JdDetect_Stationary_adj = ifelse(JdDetect_Stationary == maximum_meas_Stationary, above_median_mult_Stationary,
#                                ifelse(JdDetect_Stationary == median_JdDetect_Stationary, 1, 
#                                       NA)))
# Blocks_below_Shifting <- Blocks_below_Shifting %>% 
#   mutate(JdDetect_Shifting_adj = ifelse(JdDetect_Shifting == minimum_meas_Shifting, below_median_div_Shifting,
#                                           ifelse(JdDetect_Shifting == median_JdDetect_Shifting, 1, 
#                                                  NA)))
# Blocks_above_Shifting <- Blocks_above_Shifting %>% 
#   mutate(JdDetect_Shifting_adj = ifelse(JdDetect_Shifting == maximum_meas_Shifting, above_median_mult_Shifting,
#                                           ifelse(JdDetect_Shifting == median_JdDetect_Shifting, 1, 
#                                                  NA)))
# 
# # Use group_by and na.approx to fill in the gaps
# Blocks_below_Stationary <- Blocks_below_Stationary %>% 
#   group_by(Year) %>% 
#   mutate(JdDetect_Stationary_adj = na.approx(JdDetect_Stationary_adj)) %>% 
#   ungroup()
# Blocks_above_Stationary <- Blocks_above_Stationary %>% 
#   group_by(Year) %>% 
#   mutate(JdDetect_Stationary_adj = na.approx(JdDetect_Stationary_adj)) %>% 
#   ungroup()
# Blocks_below_Shifting <- Blocks_below_Shifting %>% 
#   group_by(Year) %>% 
#   mutate(JdDetect_Shifting_adj = na.approx(JdDetect_Shifting_adj)) %>% 
#   ungroup()
# Blocks_above_Shifting <- Blocks_above_Shifting %>% 
#   group_by(Year) %>% 
#   mutate(JdDetect_Shifting_adj = na.approx(JdDetect_Shifting_adj)) %>% 
#   ungroup()
# # Remove the median values from Blocks_above
# Blocks_above_Stationary <- Blocks_above_Stationary %>% 
#   filter(JdDetect_Stationary != median_JdDetect_Stationary)
# Blocks_above_Shifting <- Blocks_above_Shifting %>% 
#   filter(JdDetect_Shifting != median_JdDetect_Shifting)
# # Rbind together
# Blocks <- rbind(Blocks_below_Stationary, Blocks_above_Stationary)
# Shifting_df <- rbind(Blocks_below_Shifting, Blocks_above_Shifting)
# Shifting_df <- Shifting_df %>% 
#   select(Year, GridID, Common_Name, JdDetect_Shifting_adj)
# Blocks <- left_join(Blocks, Shifting_df)
# 
# nrow(Blocks) == nrow(Blocks_orig)
# # Multiply the values together to get an adjusted JdDetect
# Blocks <- Blocks %>% 
#   mutate(JdDetect_Stationary_pre = JdDetect_Stationary,
#          JdDetect_Shifting_pre = JdDetect_Shifting) %>% 
#   mutate(JdDetect_Stationary = if_else(JdDetect_Stationary <= median_JdDetect_Stationary,
#                                        JdDetect_Stationary / JdDetect_Stationary_adj,
#                                        JdDetect_Stationary * JdDetect_Stationary_adj),
#          JdDetect_Shifting = if_else(JdDetect_Shifting <= median_JdDetect_Shifting,
#                                      JdDetect_Shifting / JdDetect_Shifting_adj,
#                                      JdDetect_Shifting * JdDetect_Shifting_adj))
# # Re arrange
# Blocks <- Blocks %>% 
#   arrange(GridID, Year, Common_Name)
# 
# # # Use rescale to normalize yearly JdDetect values to the yearly min/max values
# # Blocks$JdDetect_Pre <- Blocks$JdDetect
# # 
# # Blocks <- Blocks %>% 
# #   group_by(Year) %>% 
# #   mutate(JdDetect = round((maximum_JdDetect - minimum_JdDetect) * ((JdDetect - min(JdDetect)) / (max(JdDetect) - min(JdDetect))) + minimum_JdDetect)) %>% 
# #   ungroup()
# 
# # Plot against raw data
# ggplot(Blocks, aes(JdDetect_Stationary_pre, colour = Year, group = Year)) +
#   geom_density() +
#   scale_x_continuous(limits = c(10, 106))
# ggplot(Blocks, aes(JdDetect_Stationary, colour = Year, group = Year)) +
#   geom_density() +
#   scale_x_continuous(limits = c(10, 106))
# ggplot(Blocks, aes(JdDetect_Shifting_pre, colour = Year, group = Year)) +
#   geom_density() +
#   scale_x_continuous(limits = c(10, 106))
# ggplot(Blocks, aes(JdDetect_Shifting, colour = Year, group = Year)) +
#   geom_density() +
#   scale_x_continuous(limits = c(10, 106))
# 
# ggplot(eBird_Filtered, aes(JdDetect, colour = Year, group = Year)) +
#   geom_density() +
#   scale_x_continuous(limits = c(10, 106))
# # The distribution is slightly different, but much closer than before
# # Look at the relationship between JdDetect and Year
# ggplot(Blocks, aes(factor(Year), JdDetect_Stationary_pre, colour = Year)) +
#   geom_boxplot() +
#   scale_y_continuous(limits = c(10, 106))
# ggplot(Blocks, aes(factor(Year), JdDetect_Stationary, colour = Year)) +
#   geom_boxplot() +
#   scale_y_continuous(limits = c(10, 106))
# ggplot(Blocks, aes(factor(Year), JdDetect_Shifting_pre, colour = Year)) +
#   geom_boxplot() +
#   scale_y_continuous(limits = c(10, 106))
# ggplot(Blocks, aes(factor(Year), JdDetect_Shifting, colour = Year)) +
#   geom_boxplot() +
#   scale_y_continuous(limits = c(10, 106))
# ggplot(eBird_Filtered, aes(factor(Year), JdDetect, colour = Year)) +
#   geom_boxplot() +
#   scale_y_continuous(limits = c(10, 106))



# Are sightings ever exceeding sites?
Blocks <- Blocks %>% 
  mutate(Site_Sighting_Diff = Site_Sum - Sighting_Sum)
summary(Blocks$Site_Sighting_Diff)
Blocks %>% 
  count(Site_Sighting_Diff < 0)
# True in 676 cases
# Adjust sighting sum to Grid_Mult_Site
Blocks <- Blocks %>% 
  mutate(Sighting_Sum = round(Sighting_Sum*Grid_Mult_Site),
         Site_Sighting_Diff = Site_Sum - Sighting_Sum)
summary(Blocks$Site_Sighting_Diff)
Blocks %>% 
  count(Site_Sighting_Diff < 0)
# True in 423 cases

# truncate site_sum after joining by multiplying it by Njd/100 (this makes the number smaller for species that have early departures)
Blocks <- Blocks %>% 
  mutate(Site_Sum = round(Site_Sum* (NJd/100)),
         Site_Sighting_Diff = Site_Sum - Sighting_Sum)
Blocks %>% 
  count(Site_Sighting_Diff < 0)
# True in 351 cases

# Cap any sightings tally that is 65%+ of the sight tally at 65% of 
# Calculate percent site sighting sums
Blocks <- Blocks %>% 
  mutate(Percent_Site_Sighting_Sums = Sighting_Sum / Site_Sum)
Blocks %>% 
  count(Percent_Site_Sighting_Sums > 0.65)
# Create a multiplier
Blocks <- Blocks %>% 
  mutate(Sighting_Multiplier = ifelse(Percent_Site_Sighting_Sums > 0.65,
                                       rnorm(n = n(), mean = 0.65, sd = 0.05),
                                       1)) %>% 
  mutate(Sighting_Sum = ifelse(Percent_Site_Sighting_Sums > 0.65,
                               Sighting_Multiplier * Site_Sum,
                               Sighting_Sum)) 
# Round  to nearest even integer *************
# This will help with computations later on *******
Blocks <- Blocks %>% 
  mutate(Site_Sum = ceiling(Site_Sum) - ceiling(Site_Sum) %% 2,
         Sighting_Sum = ceiling(Sighting_Sum) - ceiling(Sighting_Sum) %% 2,
         JdDetect_Stationary = ceiling(JdDetect_Stationary) - ceiling(JdDetect_Stationary) %% 2,
         JdDetect_Shifting = ceiling(JdDetect_Shifting) - ceiling(JdDetect_Shifting) %% 2)

# *********************************
# Check to make sure that JdDetect is not greater than Site_Sum
Blocks %>% 
  filter(JdDetect_Stationary >= Site_Sum)
Blocks %>% 
  filter(JdDetect_Shifting >= Site_Sum)
# Blocks to make sure that JdDetect is not greater than Sighting_Sum
Blocks %>% 
  filter(JdDetect_Stationary > Sighting_Sum)
Blocks %>% 
  filter(JdDetect_Shifting > Sighting_Sum)
# Set JdDetect equal to Sighting_Sum
Blocks <- Blocks %>% 
  mutate(JdDetect_Stationary = if_else(JdDetect_Stationary > Sighting_Sum,
                            Sighting_Sum,
                            JdDetect_Stationary))
Blocks <- Blocks %>% 
  mutate(JdDetect_Shifting = if_else(JdDetect_Shifting > Sighting_Sum,
                                       Sighting_Sum,
                                     JdDetect_Shifting))
Blocks %>% 
  filter(JdDetect_Stationary > Sighting_Sum)
Blocks %>% 
  filter(JdDetect_Stationary == Sighting_Sum)
Blocks %>% 
  filter(JdDetect_Shifting > Sighting_Sum)
Blocks %>% 
  filter(JdDetect_Shifting == Sighting_Sum)
# Subtract some random variation from JdDetect so that they are all equal to Sighting sum
# Subtract an extra 1 so that they're not equal
set.seed(1)
Blocks <- Blocks %>% 
  mutate(JdDetect_Stationary = ifelse(JdDetect_Stationary != Sighting_Sum, JdDetect_Stationary,
                           JdDetect_Stationary - round(rnormTrunc(n = n(), mean = JdDetect_Stationary/10, sd = JdDetect_Stationary/20, min = 0, max = JdDetect_Stationary/5)) - 1))
Blocks %>% 
  filter(JdDetect_Stationary >= Sighting_Sum)
ggplot(Blocks, aes(JdDetect_Stationary, Sighting_Sum)) +
  geom_point() +
  geom_abline()
ggplot(Blocks, aes(JdDetect_Stationary, Sighting_Sum)) +
  geom_point() +
  scale_y_log10()

set.seed(1)
Blocks <- Blocks %>% 
  mutate(JdDetect_Shifting = ifelse(JdDetect_Shifting != Sighting_Sum, JdDetect_Shifting,
                                      JdDetect_Shifting - round(rnormTrunc(n = n(), mean = JdDetect_Shifting/10, sd = JdDetect_Shifting/20, min = 0, max = JdDetect_Shifting/5)) - 1))
Blocks %>% 
  filter(JdDetect_Shifting >= Sighting_Sum)
ggplot(Blocks, aes(JdDetect_Shifting, Sighting_Sum)) +
  geom_point() +
  geom_abline()
ggplot(Blocks, aes(JdDetect_Shifting, Sighting_Sum)) +
  geom_point() +
  scale_y_log10()



# Get Percent
Blocks <- Blocks %>% 
  mutate(Percent_Site_Sighting_Sums = Sighting_Sum / Site_Sum)
# Plot against raw data
ggplot(Blocks, aes(Percent_Site_Sighting_Sums)) +
  geom_density()
ggplot(eBird_Filtered, aes(Percent_Site_Sighting_Sums)) +
  geom_density()
summary(Blocks$Percent_Site_Sighting_Sums)
summary(eBird_Filtered$Percent_Site_Sighting_Sums)


# Check metrics to ensure they conform with raw data
colnames(Blocks)
# Site_Sum
# Simulated
ggplot(Blocks, aes(Site_Sum, group = Year, colour = Year)) +
  geom_density()
# Real
ggplot(eBird_Filtered, aes(Site_Sum, group = Year, colour = Year)) +
  geom_density()

# Sighting_Sum
# Simulated
ggplot(Blocks, aes(Sighting_Sum, group = Year, colour = Year)) +
  geom_density()
# Real
ggplot(eBird_Filtered, aes(Sighting_Sum, group = Year, colour = Year)) +
  geom_density()

# # MinJd
# # Simulated
# ggplot(Blocks, aes(MinJd, group = Year, colour = Year)) +
#   geom_density()
# # Real
# ggplot(eBird_Filtered, aes(MinJd, group = Year, colour = Year)) +
#   geom_density()

# MaxJd
# Simulated
ggplot(Blocks, aes(MaxJd, group = Year, colour = Year)) +
  geom_density()
# Real
ggplot(eBird_Filtered, aes(MaxJd, group = Year, colour = Year)) +
  geom_density()
# Simulated
Blocks %>% 
  filter(MaxJd < 180) %>% 
ggplot(aes(MaxJd, group = Year, colour = Year)) +
  geom_density()
# Real
eBird_Filtered %>% 
  filter(MaxJd < 180) %>% 
ggplot(aes(MaxJd, group = Year, colour = Year)) +
  geom_density()

# NJd
# Simulated
ggplot(Blocks, aes(NJd, group = Year, colour = Year)) +
  geom_density()
# Real
ggplot(eBird_Filtered, aes(NJd, group = Year, colour = Year)) +
  geom_density()

# JdDetect
# Simulated
ggplot(Blocks, aes(JdDetect_Stationary, group = Year, colour = Year)) +
  geom_density()
ggplot(Blocks, aes(JdDetect_Shifting, group = Year, colour = Year)) +
  geom_density()

# Real
ggplot(eBird_Filtered, aes(JdDetect, group = Year, colour = Year)) +
  geom_density()


# Percent_Site_Sighting_Sums
# Simulated
ggplot(Blocks, aes(Percent_Site_Sighting_Sums, group = Year, colour = Year)) +
  geom_density()
# Real
ggplot(eBird_Filtered, aes(Percent_Site_Sighting_Sums, group = Year, colour = Year)) +
  geom_density()

# Simulated
ggplot(Blocks, aes(y = Percent_Site_Sighting_Sums, x = factor(Year))) +
  geom_boxplot()
# Real
ggplot(eBird_Filtered, aes(y = Percent_Site_Sighting_Sums, x = factor(Year))) +
  geom_boxplot()

# ---------------------------------------------
# Create arrival data from Blocks
# ---------------------------------------------
# Create an ID column from Common_Name, GridID, and Year
Blocks <- Blocks %>% 
  mutate(ID = paste(Common_Name, GridID, Year, sep = "_"))
Blocks_2 <- Blocks %>% 
  rowwise() %>%
  do(data.frame(ID = .$ID, Jd = seq(.$MinJd, .$MaxJd, by = 1)))
# Extract from ID to create the following columns:
# Common_Name, Year, Jd, Daily_Site_Sightings, Daily_Site_Total, Daily_Site_Occupancy, GridID
Blocks_2 <- Blocks_2 %>% 
  mutate(
    Common_Name = str_extract(ID, "[A-z]*[_][0-9]+"),
    GridID = str_extract(ID, "Grid[_][0-9]*"),
    Year = str_extract(ID, "[0-9]*$")
         )
Blocks_2$Year <- as.numeric(Blocks_2$Year)
head(Blocks_2)
# Reorder columns
Blocks_2 <- Blocks_2 %>% 
  # mutate(
  #   Daily_Site_Sightings = NA,
  #        Daily_Site_Total = NA,
  #        Daily_Site_Occupancy = NA) %>% 
  select(Common_Name, Year, Jd, 
         # Daily_Site_Sightings, Daily_Site_Total, Daily_Site_Occupancy, 
         GridID)
# Join key data to Blocks_2
Blocks_2 <- left_join(Blocks_2, Blocks)

# Join to the main data.frame
Blocks_2 <- left_join(Blocks_2, select(Grid_Years_Species, Year, GridID, Common_Name, MAD_Stationary, MAD_Shifting))



# ----------------------------------------
# Read in eBird_Daily_Occupancy_all_grids
# ----------------------------------------
eBird_Daily_Occupancy_all_grids <- fread("Routine 6/Routine 6 - merged files/eBird_Daily_Occupancy_all_grids.csv")
# # Filter the dataset
# test <- eBird_Daily_Occupancy_all_grids %>% 
#   filter(GridID %in% eBird_Filtered$GridID)


# Filter the dataset
# Create ID of filtered eBird dataset
eBird_Filtered_ID <- paste(eBird_Filtered$Common_Name, eBird_Filtered$GridID, eBird_Filtered$Year, sep = "_")
# Create ID for filtering
eBird_Daily_Occupancy_all_grids <- eBird_Daily_Occupancy_all_grids %>% 
  mutate(ID = paste(eBird_Daily_Occupancy_all_grids$Common_Name, eBird_Daily_Occupancy_all_grids$GridID, eBird_Daily_Occupancy_all_grids$Year, sep = "_"))
# Filter
orig_nrow <- nrow(eBird_Daily_Occupancy_all_grids)
eBird_Daily_Occupancy_all_grids <- eBird_Daily_Occupancy_all_grids %>% 
  filter(ID %in% eBird_Filtered_ID)
rev_nrow <- nrow(eBird_Daily_Occupancy_all_grids)
rev_nrow/orig_nrow*100

# ID
ID_2 <- unique(eBird_Daily_Occupancy_all_grids$ID)
head(setdiff(eBird_Filtered_ID, ID_2))


# Jd min is less than 80 for some species
# Just look at 80-180
eBird_Daily_Occupancy_all_grids <- eBird_Daily_Occupancy_all_grids %>% 
  filter(Jd %in% c(80:180))


eBird_Daily_Occupancy_all_grids_early <- eBird_Daily_Occupancy_all_grids %>%
  filter(Year %in% c(2002:2004))
eBird_Daily_Occupancy_all_grids_late <- eBird_Daily_Occupancy_all_grids %>%
  filter(Year %in% c(2017:2019))

# Grid sizes
eBird_Daily_Occupancy_all_grids %>% 
  group_by(GridID) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  base::print(n = 57)

# Species sizes
eBird_Daily_Occupancy_all_grids %>% 
  group_by(Common_Name) %>% 
  count() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  base::print(n = 44)

# Look at American Redstart in a large grid
eBird_Daily_Occupancy_all_grids %>% 
  filter(Common_Name == "Red-eyed Vireo",
         Jd %in% c(80:180),
         GridID == "Grid_300000x700000") %>% 
ggplot(aes(x = Jd, y = Daily_Site_Total, group = Year, colour = Year)) +
  geom_line() +
  geom_smooth() +
  facet_wrap(~Year, scales = "free_y")
# What is the peak day?
eBird_Daily_Occupancy_all_grids %>% 
  filter(Common_Name == "Red-eyed Vireo",
         Jd %in% c(80:180),
         GridID == "Grid_300000x700000") %>% 
  group_by(Year) %>% 
  filter(Daily_Site_Total == max(Daily_Site_Total)) %>% 
  ungroup() %>% 
  summary()
# Early May

# Look at American Redstart in a medium grid
eBird_Daily_Occupancy_all_grids %>% 
  filter(Common_Name == "Red-eyed Vireo",
         Jd %in% c(80:180),
         GridID == "Grid_700000x1300000") %>% 
  ggplot(aes(x = Jd, y = Daily_Site_Total, group = Year, colour = Year)) +
  geom_line() +
  geom_smooth() +
  facet_wrap(~Year, scales = "free_y")
# What is the peak day?
eBird_Daily_Occupancy_all_grids %>% 
  filter(Common_Name == "Red-eyed Vireo",
         Jd %in% c(80:180),
         GridID == "Grid_700000x1300000") %>% 
  group_by(Year) %>% 
  filter(Daily_Site_Total == max(Daily_Site_Total)) %>% 
  ungroup() %>% 
  summary(Jd)
# Early May

# Calculate the mean Daily_Site_Total for each Grid for each year
Daily_Totals_mean_early <- eBird_Daily_Occupancy_all_grids_early %>% 
  group_by(Jd, GridID) %>% 
  summarize(Daily_Site_Total_mean = round(mean(Daily_Site_Total)),
            Daily_Site_Sightings_mean = round(mean(Daily_Site_Sightings)),
            Daily_Site_Occupancy_mean = round(mean(Daily_Site_Occupancy))) %>% 
  ungroup()
Daily_Totals_mean_early

Daily_Totals_mean_late <- eBird_Daily_Occupancy_all_grids_late %>% 
  group_by(Jd, GridID) %>% 
  summarize(Daily_Site_Total_mean = round(mean(Daily_Site_Total)),
            Daily_Site_Sightings_mean = round(mean(Daily_Site_Sightings)),
            Daily_Site_Occupancy_mean = round(mean(Daily_Site_Occupancy))) %>% 
  ungroup()
Daily_Totals_mean_late
# Sample 15 grids
# set.seed(1)
# fifteen_grids <- sample(eBird_Daily_Occupancy_all_grids_early$GridID, size = 15)
# # Filter for these grids
# Daily_Totals_mean_early <- Daily_Totals_mean_early %>% 
#   filter(GridID %in% fifteen_grids)
# Daily_Totals_mean_late <- Daily_Totals_mean_late %>% 
#   filter(GridID %in% fifteen_grids)

# Create a grid size rank based on sightings
# Early
Daily_Totals_mean_early <- Daily_Totals_mean_early %>% 
  group_by(GridID) %>% 
  mutate(Grid_rank = sum(Daily_Site_Sightings_mean)) %>% 
  ungroup()
Daily_Totals_mean_early <- Daily_Totals_mean_early %>% 
  arrange(desc(Grid_rank))
unique(Daily_Totals_mean_early$Grid_rank)
# Make factor
Daily_Totals_mean_early$Grid_rank <- as.character(Daily_Totals_mean_early$Grid_rank )
Daily_Totals_mean_early$Grid_rank <- factor(Daily_Totals_mean_early$Grid_rank, 
                                            levels = unique(Daily_Totals_mean_early$Grid_rank),
                                            ordered = TRUE)
# Make grid ID 2
GridID_2 <- Daily_Totals_mean_early %>% 
  select(GridID) %>% 
  distinct() %>% 
  mutate(GridID2 = row_number())

# Join
Daily_Totals_mean_early <- left_join(Daily_Totals_mean_early, GridID_2)


# Late
Daily_Totals_mean_late <- Daily_Totals_mean_late %>% 
  group_by(GridID) %>% 
  mutate(Grid_rank = sum(Daily_Site_Sightings_mean)) %>% 
  ungroup()
Daily_Totals_mean_late <- Daily_Totals_mean_late %>%
  arrange(desc(Grid_rank))
unique(Daily_Totals_mean_late$Grid_rank)
# Make factor
Daily_Totals_mean_late$Grid_rank <- as.character(Daily_Totals_mean_late$Grid_rank )
Daily_Totals_mean_late$Grid_rank <- factor(Daily_Totals_mean_late$Grid_rank, 
                                            levels = unique(Daily_Totals_mean_late$Grid_rank),
                                            ordered = TRUE)
# Make grid ID 2
GridID_2 <- Daily_Totals_mean_late %>% 
  select(GridID) %>% 
  distinct() %>% 
  mutate(GridID2 = row_number())
# Join
Daily_Totals_mean_late <- left_join(Daily_Totals_mean_late, GridID_2)




# # Try plotting
# early_grids <- unique(Daily_Totals_mean_early$GridID)
# early_grids[1:ceiling((length(early_grids)/2))]
# last <- as.numeric(length(early_grids))
# first <- ceiling((length(early_grids)/2))+1
# early_grids[first:last]
# 
# # early
# # 1
# Daily_Totals_mean_early %>% 
#   filter(GridID %in% early_grids[1:(length(early_grids)/2)]) %>% 
ggplot(Daily_Totals_mean_early, aes(x = Jd, y = Daily_Site_Total_mean)) +
  geom_line() +
  geom_smooth(se = FALSE) +
  scale_x_continuous(limits = c(80,180)) +
  expand_limits(y = 0) +
  facet_wrap(vars(GridID2, Grid_rank), scales = "free")
# # 2
# Daily_Totals_mean_early %>% 
#   filter(GridID %in% early_grids[first:last]) %>% 
#   ggplot(aes(x = Jd, y = Daily_Site_Total_mean)) +
#   geom_line() +
#   geom_smooth(se = FALSE) +
#   scale_x_continuous(limits = c(80,180)) +
#   expand_limits(y = 0) +
#   facet_wrap(vars(GridID2, Grid_rank), scales = "free")

# late
ggplot(Daily_Totals_mean_late, aes(x = Jd, y = Daily_Site_Total_mean)) +
  geom_line() +
  geom_smooth(se = FALSE) +
  scale_x_continuous(limits = c(80,180)) +
  expand_limits(y = 0) +
  facet_wrap(vars(Grid_rank, GridID), scales = "free")

# pnorm_grids <- c(1, 2, 3, 5, 6, 10, 12, 14, 23, 31, 33, 35, 41, 43)
# rnorm_grids <- c(4, 17, 19, 20, 21, 24, 28, 29, 30, 32, 34, 36, 37, 39, 44)
# dnorm_grids <- c(7, 8, 9, 11, 13, 15, 16, 18, 22, 25, 26, 27, 38, 40, 42)

# Determine the means 
# early 
# # What is the sd of Daily_Site_Total_mean?
# Daily_Totals_mean_early <- Daily_Totals_mean_early %>% 
#   group_by(GridID) %>% 
#   summarize(mean_Daily_Site_Total_mean = mean(Daily_Site_Total_mean),
#             sd_Daily_Site_Total_mean = sd(Daily_Site_Total_mean)) %>% 
#   ungroup()
# mean(Daily_Totals_mean_early$sd_Daily_Site_Total_mean)
# sd(Daily_Totals_mean_early$sd_Daily_Site_Total_mean)
# hist(Daily_Totals_mean_early$sd_Daily_Site_Total_mean)
# # Log scale
# mean(log(Daily_Totals_mean_early$sd_Daily_Site_Total_mean))
# sd(log(Daily_Totals_mean_early$sd_Daily_Site_Total_mean))
# hist(log(Daily_Totals_mean_early$sd_Daily_Site_Total_mean))


# Continue
Daily_Site_Total_uncount_early <- Daily_Totals_mean_early %>% 
  uncount(Daily_Site_Total_mean)
Daily_Site_Total_uncount_early <- left_join(Daily_Site_Total_uncount_early, GridID_2)
Daily_Totals_mean_early <- left_join(Daily_Totals_mean_early, GridID_2)
# # pnorm-esque
# # Jd
# Daily_Site_Total_uncount_early_stats_pnorm <- Daily_Site_Total_uncount_early %>% 
#   filter(GridID2 %in% pnorm_grids) %>% 
#   group_by(GridID) %>% 
#   summarize(mean_Jd = mean(Jd),
#             sd_Jd = sd(Jd)) %>% 
#   ungroup()
# # Daily_Site_Total_mean (log)
# Daily_Site_Total_uncount_early_stats_pnorm <- left_join(Daily_Site_Total_uncount_early_stats_pnorm,
#                                                         Daily_Totals_mean_early %>% 
#   filter(GridID2 %in% pnorm_grids) %>% 
#   group_by(GridID) %>% 
#   summarize(mean_log_Daily_Site_Total_mean = mean(log(Daily_Site_Total_mean)),
#             sd_log_Daily_Site_Total_mean = sd(log(Daily_Site_Total_mean))) %>% 
#   ungroup())
# 
# summary(select(Daily_Site_Total_uncount_early_stats_pnorm, mean_Jd, sd_Jd))
# sd(Daily_Site_Total_uncount_early_stats_pnorm$mean_Jd)
# sd(Daily_Site_Total_uncount_early_stats_pnorm$sd_Jd)
# rnorm-esque
# Jd
Daily_Site_Total_uncount_early_stats_rnorm <- Daily_Site_Total_uncount_early %>% 
  # filter(GridID2 %in% rnorm_grids) %>% 
  group_by(GridID) %>% 
  summarize(mean_Jd = mean(Jd),
            sd_Jd = sd(Jd)) %>% 
  ungroup()
summary(select(Daily_Site_Total_uncount_early_stats_rnorm, mean_Jd, sd_Jd))
sd(Daily_Site_Total_uncount_early_stats_rnorm$mean_Jd)
sd(Daily_Site_Total_uncount_early_stats_rnorm$sd_Jd)
# Daily_Site_Total_mean (log)
Daily_Site_Total_uncount_early_stats_rnorm <- left_join(Daily_Site_Total_uncount_early_stats_rnorm,
                                                        Daily_Totals_mean_early %>% 
                                                          # filter(GridID2 %in% rnorm_grids) %>% 
                                                          group_by(GridID) %>% 
                                                          summarize(mean_log_Daily_Site_Total_mean = mean(log(Daily_Site_Total_mean)),
                                                                    sd_log_Daily_Site_Total_mean = sd(log(Daily_Site_Total_mean))) %>% 
                                                          ungroup())

# # dnorm-esque
# Daily_Site_Total_uncount_early_stats_dnorm <- Daily_Site_Total_uncount_early %>% 
#   filter(GridID2 %in% dnorm_grids) %>% 
#   group_by(GridID) %>% 
#   summarize(mean_Jd = mean(Jd),
#             sd_Jd = sd(Jd)) %>% 
#   ungroup()
# summary(select(Daily_Site_Total_uncount_early_stats_dnorm, mean_Jd, sd_Jd))
# sd(Daily_Site_Total_uncount_early_stats_dnorm$mean_Jd)
# sd(Daily_Site_Total_uncount_early_stats_dnorm$sd_Jd)
# # Daily_Site_Total_mean (log)
# Daily_Site_Total_uncount_early_stats_dnorm <- left_join(Daily_Site_Total_uncount_early_stats_dnorm,
#                                                         Daily_Totals_mean_early %>% 
#                                                           filter(GridID2 %in% dnorm_grids) %>% 
#                                                           group_by(GridID) %>% 
#                                                           summarize(mean_log_Daily_Site_Total_mean = mean(log(Daily_Site_Total_mean)),
#                                                                     sd_log_Daily_Site_Total_mean = sd(log(Daily_Site_Total_mean))) %>% 
#                                                           ungroup())

# Define the means/sds
# Jd
# # pnorm_early
# pnorm_early_mean_mean <- round(mean(Daily_Site_Total_uncount_early_stats_pnorm$mean_Jd), digits = 1)
# pnorm_early_mean_sd <- round(sd(Daily_Site_Total_uncount_early_stats_pnorm$mean_Jd), digits = 1)
# pnorm_early_sd_mean <- round(mean(Daily_Site_Total_uncount_early_stats_pnorm$sd_Jd), digits = 1)
# pnorm_early_sd_sd <- round(sd(Daily_Site_Total_uncount_early_stats_pnorm$sd_Jd), digits = 1)
# rnorm_early
rnorm_early_mean_mean <- round(mean(Daily_Site_Total_uncount_early_stats_rnorm$mean_Jd), digits = 1)
rnorm_early_mean_sd <- round(sd(Daily_Site_Total_uncount_early_stats_rnorm$mean_Jd), digits = 1)
rnorm_early_sd_mean <- round(mean(Daily_Site_Total_uncount_early_stats_rnorm$sd_Jd), digits = 1)
rnorm_early_sd_sd <- round(sd(Daily_Site_Total_uncount_early_stats_rnorm$sd_Jd), digits = 1)
# # dnorm_early
# dnorm_early_mean_mean <- round(mean(Daily_Site_Total_uncount_early_stats_dnorm$mean_Jd), digits = 1)
# dnorm_early_mean_sd <- round(sd(Daily_Site_Total_uncount_early_stats_dnorm$mean_Jd), digits = 1)
# dnorm_early_sd_mean <- round(mean(Daily_Site_Total_uncount_early_stats_dnorm$sd_Jd), digits = 1)
# dnorm_early_sd_sd <- round(sd(Daily_Site_Total_uncount_early_stats_dnorm$sd_Jd), digits = 1)
# log_Daily_Site_Total_mean
# # pnorm_early
# pnorm_early_mean_mean_2 <- round(mean(Daily_Site_Total_uncount_early_stats_pnorm$mean_log_Daily_Site_Total_mean), digits = 1)
# pnorm_early_mean_sd_2 <- round(sd(Daily_Site_Total_uncount_early_stats_pnorm$mean_log_Daily_Site_Total_mean), digits = 1)
# pnorm_early_sd_mean_2 <- round(mean(Daily_Site_Total_uncount_early_stats_pnorm$sd_log_Daily_Site_Total_mean), digits = 1)
# pnorm_early_sd_sd_2 <- round(sd(Daily_Site_Total_uncount_early_stats_pnorm$sd_log_Daily_Site_Total_mean), digits = 1)
# rnorm_early
rnorm_early_mean_mean_2 <- round(mean(Daily_Site_Total_uncount_early_stats_rnorm$mean_log_Daily_Site_Total_mean), digits = 1)
rnorm_early_mean_sd_2 <- round(sd(Daily_Site_Total_uncount_early_stats_rnorm$mean_log_Daily_Site_Total_mean), digits = 1)
rnorm_early_sd_mean_2 <- round(mean(Daily_Site_Total_uncount_early_stats_rnorm$sd_log_Daily_Site_Total_mean), digits = 1)
rnorm_early_sd_sd_2 <- round(sd(Daily_Site_Total_uncount_early_stats_rnorm$sd_log_Daily_Site_Total_mean), digits = 1)
# # dnorm_early
# dnorm_early_mean_mean_2 <- round(mean(Daily_Site_Total_uncount_early_stats_dnorm$mean_log_Daily_Site_Total_mean), digits = 1)
# dnorm_early_mean_sd_2 <- round(sd(Daily_Site_Total_uncount_early_stats_dnorm$mean_log_Daily_Site_Total_mean), digits = 1)
# dnorm_early_sd_mean_2 <- round(mean(Daily_Site_Total_uncount_early_stats_dnorm$sd_log_Daily_Site_Total_mean), digits = 1)
# dnorm_early_sd_sd_2 <- round(sd(Daily_Site_Total_uncount_early_stats_dnorm$sd_log_Daily_Site_Total_mean), digits = 1)



# late
Daily_Site_Total_uncount_late <- Daily_Totals_mean_late %>% 
  uncount(Daily_Site_Total_mean)
# Jd
Daily_Site_Total_uncount_late_stats <- Daily_Site_Total_uncount_late %>% 
  group_by(GridID) %>% 
  summarize(mean_Jd = mean(Jd),
            sd_Jd = sd(Jd)) %>% 
  ungroup()
summary(select(Daily_Site_Total_uncount_late_stats, mean_Jd, sd_Jd))
sd(Daily_Site_Total_uncount_late_stats$mean_Jd)
sd(Daily_Site_Total_uncount_late_stats$sd_Jd)
# Daily_Site_Total_mean (log)
Daily_Site_Total_uncount_late_stats <- left_join(Daily_Site_Total_uncount_late_stats,
                                                        Daily_Totals_mean_late %>% 
                                                          group_by(GridID) %>% 
                                                          summarize(mean_log_Daily_Site_Total_mean = mean(log(Daily_Site_Total_mean)),
                                                                    sd_log_Daily_Site_Total_mean = sd(log(Daily_Site_Total_mean))) %>% 
                                                          ungroup())

# Define the means/sds
# late
# Jd
late_mean_mean <- round(mean(Daily_Site_Total_uncount_late_stats$mean_Jd), digits = 1)
late_mean_sd <- round(sd(Daily_Site_Total_uncount_late_stats$mean_Jd), digits = 1)
late_sd_mean <- round(mean(Daily_Site_Total_uncount_late_stats$sd_Jd), digits = 1)
late_sd_sd <- round(sd(Daily_Site_Total_uncount_late_stats$sd_Jd), digits = 1)
# log_Daily_Site_Total_mean
late_mean_mean_2 <- round(mean(Daily_Site_Total_uncount_late_stats$mean_log_Daily_Site_Total_mean), digits = 1)
late_mean_sd_2 <- round(sd(Daily_Site_Total_uncount_late_stats$mean_log_Daily_Site_Total_mean), digits = 1)
late_sd_mean_2 <- round(mean(Daily_Site_Total_uncount_late_stats$sd_log_Daily_Site_Total_mean), digits = 1)
late_sd_sd_2 <- round(sd(Daily_Site_Total_uncount_late_stats$sd_log_Daily_Site_Total_mean), digits = 1)

# Jd
# # pnorm_middle
# pnorm_middle_mean_mean <- round((pnorm_early_mean_mean + late_mean_mean) / 2, digits = 1)
# pnorm_middle_mean_sd <- round((pnorm_early_mean_sd + late_mean_sd) / 2, digits = 1)
# pnorm_middle_sd_mean <- round((pnorm_early_sd_mean + late_sd_mean) / 2, digits = 1)
# pnorm_middle_sd_sd <- round((pnorm_early_sd_sd + late_sd_sd) / 2, digits = 1)
# rnorm_middle
rnorm_middle_mean_mean <- round((rnorm_early_mean_mean + late_mean_mean) / 2, digits = 1)
rnorm_middle_mean_sd <- round((rnorm_early_mean_sd + late_mean_sd) / 2, digits = 1)
rnorm_middle_sd_mean <- round((rnorm_early_sd_mean + late_sd_mean) / 2, digits = 1)
rnorm_middle_sd_sd <- round((rnorm_early_sd_sd + late_sd_sd) / 2, digits = 1)
# # dnorm_middle
# dnorm_middle_mean_mean <- round((dnorm_early_mean_mean + late_mean_mean) / 2, digits = 1)
# dnorm_middle_mean_sd <- round((dnorm_early_mean_sd + late_mean_sd) / 2, digits = 1)
# dnorm_middle_sd_mean <- round((dnorm_early_sd_mean + late_sd_mean) / 2, digits = 1)
# dnorm_middle_sd_sd <- round((dnorm_early_sd_sd + late_sd_sd) / 2, digits = 1)

# log_Daily_Site_Total_mean
# # pnorm_middle
# pnorm_middle_mean_mean_2 <- round((pnorm_early_mean_mean_2 + late_mean_mean_2) / 2, digits = 1)
# pnorm_middle_mean_sd_2 <- round((pnorm_early_mean_sd_2 + late_mean_sd_2) / 2, digits = 1)
# pnorm_middle_sd_mean_2 <- round((pnorm_early_sd_mean_2 + late_sd_mean_2) / 2, digits = 1)
# pnorm_middle_sd_sd_2 <- round((pnorm_early_sd_sd_2 + late_sd_sd_2) / 2, digits = 1)
# rnorm_middle
rnorm_middle_mean_mean_2 <- round((rnorm_early_mean_mean_2 + late_mean_mean_2) / 2, digits = 1)
rnorm_middle_mean_sd_2 <- round((rnorm_early_mean_sd_2 + late_mean_sd_2) / 2, digits = 1)
rnorm_middle_sd_mean_2 <- round((rnorm_early_sd_mean_2 + late_sd_mean_2) / 2, digits = 1)
rnorm_middle_sd_sd_2 <- round((rnorm_early_sd_sd_2 + late_sd_sd_2) / 2, digits = 1)
# # dnorm_middle
# dnorm_middle_mean_mean_2 <- round((dnorm_early_mean_mean_2 + late_mean_mean_2) / 2, digits = 1)
# dnorm_middle_mean_sd_2 <- round((dnorm_early_mean_sd_2 + late_mean_sd_2) / 2, digits = 1)
# dnorm_middle_sd_mean_2 <- round((dnorm_early_sd_mean_2 + late_sd_mean_2) / 2, digits = 1)
# dnorm_middle_sd_sd_2 <- round((dnorm_early_sd_sd_2 + late_sd_sd_2) / 2, digits = 1)



# Plot some of the stats
ggplot(Daily_Site_Total_uncount_late_stats, aes(mean_Jd)) +
  geom_histogram(binwidth = 1)
ggplot(Daily_Site_Total_uncount_late_stats, aes(sd_Jd)) +
  geom_histogram(binwidth = 0.25)
  

# Summary of distributions:
# Early:
# - Some distributions are a very sharp cumulative normal distribution 
# - Some are pretty randomly distributed, could just be random samples from a normal distribution
# - Some increase up until Jd ~= 150, then decline
# - Could do: 1/3 = pnorm, 1/3 = rnorm, and 1/3 = dnorm in early
# - Then, 1/6 = pnorm, 1/6 = rnorm, and 2/3 = dnorm in mid
# - Then 100% dnorm in late

# x = seq(80, 180, 1) (minJd to maxJd)
# mean = 137 (extremes: 124-159; 1st/3rd quantiles: 130-141; sd = 10)
# sd = 26 (extremes: 15-31; 1st/3rd quantiles: 26-28; sd = 3)

# Late:
# - Pretty much all distributions increase up until 125-150, then decrease
# function: dnorm
# x = seq(80, 180, 1) (minJd to maxJd)
# mean = 130 (extremes: 123-141; 1st/3rd quantiles: 127-132; sd = 4)
# sd = 27 (extremes: 25-29; 1st/3rd quantiles: 26-27; sd = 0.8)
# # Test late statistics
# # test
# set.seed(1)
# test_vec_mean_late <- round(rnorm(n = 57, mean = 130, sd = 4))
# summary(test_vec_mean_late)
# 
# 
# set.seed(1)
# test <- data.frame(x = seq(80, 180, 1),
#            y = dnorm(x = seq(80, 180, 1),
#                      mean = 140,
#                      sd = 27))
# ggplot(test, aes(x = x, y = y)) +
#   geom_point() +
#   scale_x_continuous(limits = c(80, 180)) +
#   expand_limits(y = 0)
# # test2
# test2 <- data.frame(seq_from = rep(80, Grids),
#                     seq_to = rep(180, Grids),
#                     seq_by = rep(1, Grids),
#                     mean = round(rnorm(n = Grids, mean = 130, sd = 4)),
#                     sd = round(rnorm(n = Grids, mean = 27, sd = 0.8)))
# # test3
# round(dnorm(x = seq(test2$seq_from[1], test2$seq_to[1], test2$seq_by[1]),
#       mean = test2$mean[1],
#       sd = test2$sd[1])*1000)


# -------------------------------------------------------
# Create Grid_Years with sampling data
# -------------------------------------------------------
head(Grid_Years)

# Randomly sample grids to place in each of the sampling cohorts
# 3 sampling cohorts for early and middle late: pnorm (1), rnorm (2), and dnorm (3)
# 1 sampling cohort for late: dnorm only (1)
# set.seed(1)
# early_grid_cohorts <- split(sample(unique(Grid_Years$GridID)), 1:3)
# set.seed(1)
# middle_grid_cohorts <- split(sample(unique(Grid_Years$GridID)), 1:6)
# middle_grid_cohorts$`3` <- c(middle_grid_cohorts$`3`, middle_grid_cohorts$`4`, middle_grid_cohorts$`5`, middle_grid_cohorts$`6`)
# middle_grid_cohorts$`4` <-  middle_grid_cohorts$`5` <-  middle_grid_cohorts$`6` <- NULL
# late_grid_cohorts <- unique(Grid_Years$GridID)


# Create info about data.frames
Grid_Years <- Grid_Years %>% 
  mutate(time_cohort = case_when(Year %in% c(1:6) ~ "Early",
                                 Year %in% c(7:12) ~ "Middle",
                                 Year %in% c(13:18) ~ "Late"))
# Grid_Years <- Grid_Years %>% 
#   mutate(distribution_cohort = case_when(
#     time_cohort == "Early" & GridID %in% early_grid_cohorts$`1` ~ "pnorm",
#                                          time_cohort == "Early" & GridID %in% early_grid_cohorts$`2` ~ "rnorm",
#                                          time_cohort == "Early" & GridID %in% early_grid_cohorts$`3` ~ "dnorm",
#                                          time_cohort == "Middle" & GridID %in% middle_grid_cohorts$`1` ~ "pnorm",
#                                          time_cohort == "Middle" & GridID %in% middle_grid_cohorts$`2` ~ "rnorm",
#                                          time_cohort == "Middle" & GridID %in% middle_grid_cohorts$`3` ~ "dnorm",
#                                          time_cohort == "Late" ~ "dnorm"
#     ))

# Join these data to Blocks_2
Grid_Years_select <- Grid_Years %>% 
  select(Year, GridID, time_cohort
         # , distribution_cohort
         )
Blocks_2 <- left_join(Blocks_2, Grid_Years_select)



# ----------------------------------
# Create data in Grid_Years
# ----------------------------------
# Early:
# - Some distributions are a very sharp cumulative normal distribution 
# - Some are pretty randomly distributed, could just be random samples from a normal distribution
# - Some increase up until Jd ~= 150, then decline
# - Could do: 1/3 = pnorm, 1/3 = rnorm, and 1/3 = dnorm in early
# - Then, 1/6 = pnorm, 1/6 = rnorm, and 2/3 = dnorm in mid
# - Then 100% dnorm in late

# x = seq(80, 180, 1) (minJd to maxJd)
# mean = 137 (extremes: 124-159; 1st/3rd quantiles: 130-141; sd = 10)
# sd = 26 (extremes: 15-31; 1st/3rd quantiles: 26-28; sd = 3)

Grid_Years_2 <- left_join(Grid_Years_2,
                          distinct(select(Blocks,
                                          GridID, Year, Site_Sum)))

# Late:
# - Pretty much all distributions increase up until 125-150, then decrease
# function: dnorm
# x = seq(80, 180, 1) (minJd to maxJd)
# mean = 130 (extremes: 123-141; 1st/3rd quantiles: 127-132; sd = 4)
# sd = 27 (extremes: 25-29; 1st/3rd quantiles: 26-27; sd = 0.8)
# Create x distributions according to the number we need for each grid
# Create early grid data.frames
early_years <- c(1:6)
early_grid_cohort <- Grid_Years_2 %>%
  filter(Year %in% early_years)
# early_grid_cohort_pnorm <- Grid_Years_2 %>%
#   filter(Year %in% early_years,
#          GridID %in% early_grid_cohorts$`1`)
# early_grid_cohort_rnorm <- Grid_Years_2 %>%
#   filter(Year %in% early_years,
#          GridID %in% early_grid_cohorts$`2`)
# early_grid_cohort_dnorm <- Grid_Years_2 %>%
#   filter(Year %in% early_years,
#          GridID %in% early_grid_cohorts$`3`)

# Create middle grid data.frames
middle_years <- c(7:12)
middle_grid_cohort <- Grid_Years_2 %>%
  filter(Year %in% middle_years)
# middle_grid_cohort_pnorm <- Grid_Years_2 %>%
#   filter(Year %in% middle_years,
#          GridID %in% middle_grid_cohorts$`1`)
# middle_grid_cohort_rnorm <- Grid_Years_2 %>%
#   filter(Year %in% middle_years,
#          GridID %in% middle_grid_cohorts$`2`)
# middle_grid_cohort_dnorm <- Grid_Years_2 %>%
#   filter(Year %in% middle_years,
#          GridID %in% middle_grid_cohorts$`3`)

# Create late grid data.frames
late_years <- c(13:18)
late_grid_cohort <- Grid_Years_2 %>%
  filter(Year %in% late_years)
# late_grid_cohort_dnorm <- Grid_Years_2 %>%
#   filter(Year %in% late_years,
#          GridID %in% late_grid_cohorts)


# # early_grid_cohort_pnorm
# # Create values for Daily_Site_Total
# head(early_grid_cohort_pnorm)
# early_grid_cohort_pnorm_grids <- early_grid_cohort_pnorm %>% 
#   group_by(GridID) %>% 
#   mutate(Site_Sum_mean = mean(Site_Sum)) %>% 
#   ungroup() %>% 
#   select(GridID, Site_Sum_mean, Grid_Adj_mean, Grid_Adj_sd) %>% 
#   distinct() %>% 
#   arrange(-desc(Site_Sum_mean))
# # stats
# set.seed(1)
# early_pnorm_grids_stats <- data.frame(mean = round(rnorm(n = nrow(early_grid_cohort_pnorm_grids), mean = pnorm_early_mean_mean, sd = pnorm_early_mean_sd)),
#          sd = round(rnorm(n = nrow(early_grid_cohort_pnorm_grids), mean = pnorm_early_sd_mean, sd = pnorm_early_sd_sd)))
# early_pnorm_grids_stats <- early_pnorm_grids_stats %>% 
#   arrange(desc(mean))
# # stats_2
# set.seed(1)
# early_pnorm_grids_stats_2 <- data.frame(mean_2 = round(rlnorm(n = nrow(early_grid_cohort_pnorm_grids), meanlog = pnorm_early_mean_mean_2, sdlog = pnorm_early_mean_sd_2)),
# sd_2 = round(rlnorm(n = nrow(early_grid_cohort_pnorm_grids), meanlog = pnorm_early_sd_mean_2, sdlog = pnorm_early_sd_sd_2)))
# early_pnorm_grids_stats_2 <- early_pnorm_grids_stats_2 %>% 
#   arrange(-desc(mean_2))
# # cbind
# early_grid_cohort_pnorm_grids <- cbind(early_grid_cohort_pnorm_grids, early_pnorm_grids_stats, early_pnorm_grids_stats_2)
# # early_grid_cohort_pnorm_grids$Site_Sum_mean <- NULL
# # Add grid_adj
# early_grid_cohort_pnorm_grids <- early_grid_cohort_pnorm_grids %>%
#   mutate(mean = mean + Grid_Adj_mean,
#          sd = sd + Grid_Adj_sd)
# 
# early_grid_cohort_pnorm <- left_join(early_grid_cohort_pnorm, early_grid_cohort_pnorm_grids)
# # Create Daily_Site_Total column
# early_grid_cohort_pnorm$Daily_Site_Total <- NA
# early_grid_cohort_pnorm$Daily_Site_Total <- as.numeric(early_grid_cohort_pnorm$Daily_Site_Total)
# # Loop
# for (i in 1:length(unique(early_grid_cohort_pnorm$ID))) {
# # i = 1
# i0 = i - 1
#   ID_i <- unique(early_grid_cohort_pnorm$ID)[i]
#   df_i <- filter(early_grid_cohort_pnorm, ID == ID_i)
#   set.seed(i)
#   vec_i <- pnorm(q = seq(unique(df_i$MinimumJd), unique(df_i$MaximumJd), 1),
#                  mean = unique(df_i$mean),
#                  sd = unique(df_i$sd))*100
#   # Get vec_i sum
#   # Correct vec_i to Site_Sum
#   site_sum_i <- unique(df_i$Site_Sum)
#   correction_factor <- site_sum_i/sum(vec_i)
#   vec_i <- vec_i * correction_factor
#   # Introduce some random variability
#   set.seed(i)
#   vec_i <- vec_i + rnorm(n = length(c(unique(df_i$MinimumJd):unique(df_i$MaximumJd))), mean = 0, sd = unique(df_i$sd_2))
#   # Correct negative values to zero
#   vec_i <- ifelse(vec_i < 0, 0, vec_i)
#   # Correct 1 values to 2
#   vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
#   # How many non-zero days are there
#   length_above_zero <- length(vec_i[vec_i > 0])
#   # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
#   # Index the minimum spot
#   min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
#   vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
#   # Get vec_i sum again
#   # How many sites is this?
#   correction_factor_2 <- site_sum_i / sum(vec_i)
#   # Multiply the vector by the correction factor again
#   vec_i <- vec_i * correction_factor_2
#   # Round the vector
#   vec_i <- round(vec_i)
#   # What is the absolute difference between the site numbers?
#   vec_diff <- site_sum_i - sum(vec_i)
#   # Add/subtract this to a position at the max value positions
#   vec_i_df <- data.frame(vec_i = vec_i)
#   vec_i_df <- vec_i_df %>%
#     mutate(rowid = row_number())
# 
#   if(vec_diff < 0){
#     # Max rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_max(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Subtract from max rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i - 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else if(vec_diff > 0) {
#     # Min rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_min(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Add to min rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i + 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else {
#     vec_i <- vec_i
#   }
# 
# 
#   # Add to Daily_Site_Total column (15)
#   start_ind <- 106*i0 + 1
#   end_ind <- 106*i0 + 1 + 105
#   early_grid_cohort_pnorm[c(start_ind:end_ind), 18] <- vec_i
# }
# # Plot
# early_grid_cohort_pnorm %>%
#   filter(Year == 1) %>%
# ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID)
# # Plot
# early_grid_cohort_pnorm %>%
#   filter(Year == 1) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID, scales = "free_y")
# # Plot
# early_grid_cohort_pnorm %>%
#   filter(Year == 1) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0)

# # early_grid_cohort_rnorm
# # Create values for Daily_Site_Total
# head(early_grid_cohort_rnorm)
# early_grid_cohort_rnorm_grids <- early_grid_cohort_rnorm %>%
#   group_by(GridID) %>%
#   mutate(Site_Sum_mean = mean(Site_Sum)) %>%
#   ungroup() %>%
#   select(GridID, Site_Sum_mean, Grid_Adj_mean, Grid_Adj_sd) %>%
#   distinct() %>%
#   arrange(-desc(Site_Sum_mean))
# set.seed(1)
# # stats
# early_rnorm_grids_stats <- data.frame(mean = round(rnorm(n = nrow(early_grid_cohort_rnorm_grids), mean = rnorm_early_mean_mean, sd = rnorm_early_mean_sd)),
#                                       sd = round(rnorm(n = nrow(early_grid_cohort_rnorm_grids), mean = rnorm_early_sd_mean, sd = rnorm_early_sd_sd)))
# set.seed(1)
# early_rnorm_grids_stats <- early_rnorm_grids_stats %>%
#   arrange(desc(mean))
# # stats_2
# early_rnorm_grids_stats_2 <- data.frame(mean_2 = round(rlnorm(n = nrow(early_grid_cohort_rnorm_grids), meanlog = rnorm_early_mean_mean_2, sdlog = rnorm_early_mean_sd_2)),
#                                         sd_2 = round(rlnorm(n = nrow(early_grid_cohort_rnorm_grids), meanlog = rnorm_early_sd_mean_2, sdlog = rnorm_early_sd_sd_2)))
# early_rnorm_grids_stats_2 <- early_rnorm_grids_stats_2 %>%
#   arrange(-desc(mean_2))
# # cbind
# early_grid_cohort_rnorm_grids <- cbind(early_grid_cohort_rnorm_grids, early_rnorm_grids_stats, early_rnorm_grids_stats_2)
# # Add grid_adj
# early_grid_cohort_rnorm_grids <- early_grid_cohort_rnorm_grids %>%
#   mutate(mean = mean + Grid_Adj_mean,
#          sd = sd + Grid_Adj_sd)
# 
# early_grid_cohort_rnorm <- left_join(early_grid_cohort_rnorm, early_grid_cohort_rnorm_grids)
# # Create Daily_Site_Total column
# early_grid_cohort_rnorm$Daily_Site_Total <- NA
# early_grid_cohort_rnorm$Daily_Site_Total <- as.numeric(early_grid_cohort_rnorm$Daily_Site_Total)
# # Loop
# for (i in 1:length(unique(early_grid_cohort_rnorm$ID))) {
#   # i = 91
#   i0 = i - 1
#   ID_i <- unique(early_grid_cohort_rnorm$ID)[i]
#   df_i <- filter(early_grid_cohort_rnorm, ID == ID_i)
#   set.seed(i)
#   vec_i <- round(rnorm(n = seq(unique(df_i$MinimumJd), unique(df_i$MaximumJd), 1),
#                        mean = unique(df_i$mean_2),
#                        sd = unique(df_i$sd_2)))
#   # Get vec_i sum
#   # Correct vec_i to Site_Sum
#   site_sum_i <- unique(df_i$Site_Sum)
#   correction_factor <- site_sum_i/sum(vec_i)
#   vec_i <- vec_i * correction_factor
#   # Introduce some random variability
#   set.seed(i)
#   vec_i <- vec_i + rnorm(n = length(c(unique(df_i$MinimumJd):unique(df_i$MaximumJd))), mean = 0, sd = unique(df_i$sd_2))
#   # Correct negative values to zero
#   vec_i <- ifelse(vec_i < 0, 0, vec_i)
#   # Correct 1 values to 2
#   vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
#   # How many non-zero days are there
#   length_above_zero <- length(vec_i[vec_i > 0])
#   # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
#   # Index the minimum spot
#   min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
#   vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
#   # Get vec_i sum again
#   # How many sites is this?
#   correction_factor_2 <- site_sum_i / sum(vec_i)
#   # Multiply the vector by the correction factor again
#   vec_i <- vec_i * correction_factor_2
#   # Round the vector
#   vec_i <- round(vec_i)
#   # What is the absolute difference between the site numbers?
#   vec_diff <- site_sum_i - sum(vec_i)
#   # Add/subtract this to a position at the max value positions
#   vec_i_df <- data.frame(vec_i = vec_i)
#   vec_i_df <- vec_i_df %>%
#     mutate(rowid = row_number())
# 
#   if(vec_diff < 0){
#     # Max rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_max(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Subtract from max rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i - 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else if(vec_diff > 0) {
#     # Min rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_min(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Add to min rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i + 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else {
#     vec_i <- vec_i
#   }
# 
#   # max_pos_i <- which.max(vec_i)
#   # vec_i[max_pos_i] <- vec_i[max_pos_i] + vec_diff
#   # Add to Daily_Site_Total column (15)
#   # Add to Daily_Site_Total column (15)
#   start_ind <- 106*i0 + 1
#   end_ind <- 106*i0 + 1 + 105
#   early_grid_cohort_rnorm[c(start_ind:end_ind), 18] <- vec_i
# }
# # Plot
# early_grid_cohort_rnorm %>%
#   filter(Year == 1) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID)
# # Plot
# early_grid_cohort_rnorm %>%
#   filter(Year == 1) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID, scales = "free_y")
# # Plot
# early_grid_cohort_rnorm %>%
#   filter(Year == 1) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0)


# early_grid_cohort
# Create values for Daily_Site_Total
head(early_grid_cohort)
early_grid_cohort_grids <- early_grid_cohort %>%
  group_by(GridID) %>%
  mutate(Site_Sum_mean = mean(Site_Sum)) %>%
  ungroup() %>%
  select(GridID, Site_Sum_mean, Grid_Adj_mean, Grid_Adj_sd) %>%
  distinct() %>%
  arrange(-desc(Site_Sum_mean))
set.seed(1)
# stats
early_grids_stats <- data.frame(mean = round(rnorm(n = nrow(early_grid_cohort_grids), mean = rnorm_early_mean_mean, sd = rnorm_early_mean_sd)),
                                      sd = round(rnorm(n = nrow(early_grid_cohort_grids), mean = rnorm_early_sd_mean, sd = rnorm_early_sd_sd)))
set.seed(1)
early_grids_stats <- early_grids_stats %>%
  arrange(desc(mean))
# stats_2
early_grids_stats_2 <- data.frame(mean_2 = round(rlnorm(n = nrow(early_grid_cohort_grids), meanlog = rnorm_early_mean_mean_2, sdlog = rnorm_early_mean_sd_2)),
                                        sd_2 = round(rlnorm(n = nrow(early_grid_cohort_grids), meanlog = rnorm_early_sd_mean_2, sdlog = rnorm_early_sd_sd_2)))
early_grids_stats_2 <- early_grids_stats_2 %>%
  arrange(-desc(mean_2))
# cbind
early_grid_cohort_grids <- cbind(early_grid_cohort_grids, early_grids_stats, early_grids_stats_2)
# Add grid_adj
early_grid_cohort_grids <- early_grid_cohort_grids %>%
  mutate(mean = mean + Grid_Adj_mean,
         sd = sd + Grid_Adj_sd)

early_grid_cohort <- left_join(early_grid_cohort, early_grid_cohort_grids)
# Create Daily_Site_Total column
early_grid_cohort$Daily_Site_Total <- NA
early_grid_cohort$Daily_Site_Total <- as.numeric(early_grid_cohort$Daily_Site_Total)
# Loop
for (i in 1:length(unique(early_grid_cohort$ID))) {
  # i = 91
  i0 = i - 1
  ID_i <- unique(early_grid_cohort$ID)[i]
  df_i <- filter(early_grid_cohort, ID == ID_i)
  set.seed(i)
  vec_i <- round(rnorm(n = seq(unique(df_i$MinimumJd), unique(df_i$MaximumJd), 1),
                       mean = unique(df_i$mean_2),
                       sd = unique(df_i$sd_2)))
  # Get vec_i sum
  # Correct vec_i to Site_Sum
  site_sum_i <- unique(df_i$Site_Sum)
  correction_factor <- site_sum_i/sum(vec_i)
  vec_i <- vec_i * correction_factor
  # Introduce some random variability
  set.seed(i)
  vec_i <- vec_i + rnorm(n = length(c(unique(df_i$MinimumJd):unique(df_i$MaximumJd))), mean = 0, sd = unique(df_i$sd_2))
  # Correct negative values to zero
  vec_i <- ifelse(vec_i < 0, 0, vec_i)
  # Correct 1 values to 2
  vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
  # How many non-zero days are there
  length_above_zero <- length(vec_i[vec_i > 0])
  # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
  # Index the minimum spot
  min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
  vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
  # Get vec_i sum again
  # How many sites is this?
  correction_factor_2 <- site_sum_i / sum(vec_i)
  # Multiply the vector by the correction factor again
  vec_i <- vec_i * correction_factor_2
  # Round the vector
  vec_i <- round(vec_i)
  # What is the absolute difference between the site numbers?
  vec_diff <- site_sum_i - sum(vec_i)
  # Add/subtract this to a position at the max value positions
  vec_i_df <- data.frame(vec_i = vec_i)
  vec_i_df <- vec_i_df %>%
    mutate(rowid = row_number())
  
  if(vec_diff < 0){
    # Max rows
    vec_i_df_subset <- vec_i_df %>%
      slice_max(vec_i, n = abs(vec_diff), with_ties = FALSE)
    # Other rows
    vec_i_df_not_subset <- vec_i_df %>%
      filter(!rowid %in% vec_i_df_subset$rowid)
    # Subtract from max rows
    vec_i_df_subset <- vec_i_df_subset %>%
      mutate(vec_i = vec_i - 1)
    # Bind data.frames
    vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
    # Arrange data.frame
    vec_i_df <- vec_i_df %>%
      arrange(-desc(rowid))
    # Get vector
    vec_i <- vec_i_df$vec_i
    
  } else if(vec_diff > 0) {
    # Min rows
    vec_i_df_subset <- vec_i_df %>%
      slice_min(vec_i, n = abs(vec_diff), with_ties = FALSE)
    # Other rows
    vec_i_df_not_subset <- vec_i_df %>%
      filter(!rowid %in% vec_i_df_subset$rowid)
    # Add to min rows
    vec_i_df_subset <- vec_i_df_subset %>%
      mutate(vec_i = vec_i + 1)
    # Bind data.frames
    vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
    # Arrange data.frame
    vec_i_df <- vec_i_df %>%
      arrange(-desc(rowid))
    # Get vector
    vec_i <- vec_i_df$vec_i
    
  } else {
    vec_i <- vec_i
  }
  
  # max_pos_i <- which.max(vec_i)
  # vec_i[max_pos_i] <- vec_i[max_pos_i] + vec_diff
  # Add to Daily_Site_Total column (15)
  # Add to Daily_Site_Total column (15)
  start_ind <- 106*i0 + 1
  end_ind <- 106*i0 + 1 + 105
  early_grid_cohort[c(start_ind:end_ind), 16] <- vec_i
}
# Plot
early_grid_cohort %>%
  filter(Year == 1) %>%
  ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
  geom_line() +
  expand_limits(y = 0) +
  facet_wrap(~GridID)
# Plot
early_grid_cohort %>%
  filter(Year == 1) %>%
  ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
  geom_line() +
  expand_limits(y = 0) +
  facet_wrap(~GridID, scales = "free_y")
# Plot
early_grid_cohort %>%
  filter(Year == 1) %>%
  ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
  geom_line() +
  expand_limits(y = 0)

# # early_grid_cohort_dnorm
# # Create values for Daily_Site_Total
# head(early_grid_cohort_dnorm)
# early_grid_cohort_dnorm_grids <- early_grid_cohort_dnorm %>%
#   group_by(GridID) %>%
#   mutate(Site_Sum_mean = mean(Site_Sum)) %>%
#   ungroup() %>%
#   select(GridID, Site_Sum_mean, Grid_Adj_mean, Grid_Adj_sd) %>%
#   distinct() %>%
#   arrange(-desc(Site_Sum_mean))
# set.seed(1)
# # stats
# early_dnorm_grids_stats <- data.frame(mean = round(rnorm(n = nrow(early_grid_cohort_dnorm_grids), mean = dnorm_early_mean_mean, sd = dnorm_early_mean_sd)),
#                                       sd = round(rnorm(n = nrow(early_grid_cohort_dnorm_grids), mean = dnorm_early_sd_mean, sd = dnorm_early_sd_sd)))
# set.seed(1)
# early_dnorm_grids_stats <- early_dnorm_grids_stats %>%
#   arrange(desc(mean))
# # stats_2
# early_dnorm_grids_stats_2 <- data.frame(mean_2 = round(rlnorm(n = nrow(early_grid_cohort_dnorm_grids), meanlog = dnorm_early_mean_mean_2, sdlog = dnorm_early_mean_sd_2)),
#                                         sd_2 = round(rlnorm(n = nrow(early_grid_cohort_dnorm_grids), meanlog = dnorm_early_sd_mean_2, sdlog = dnorm_early_sd_sd_2)))
# early_dnorm_grids_stats_2 <- early_dnorm_grids_stats_2 %>%
#   arrange(-desc(mean_2))
# # cbind
# early_grid_cohort_dnorm_grids <- cbind(early_grid_cohort_dnorm_grids, early_dnorm_grids_stats, early_dnorm_grids_stats_2)
# # Add grid_adj
# early_grid_cohort_dnorm_grids <- early_grid_cohort_dnorm_grids %>%
#   mutate(mean = mean + Grid_Adj_mean,
#          sd = sd + Grid_Adj_sd)
# 
# early_grid_cohort_dnorm <- left_join(early_grid_cohort_dnorm, early_grid_cohort_dnorm_grids)
# # Create Daily_Site_Total column
# early_grid_cohort_dnorm$Daily_Site_Total <- NA
# early_grid_cohort_dnorm$Daily_Site_Total <- as.numeric(early_grid_cohort_dnorm$Daily_Site_Total)
# # Loop
# for (i in 1:length(unique(early_grid_cohort_dnorm$ID))) {
#   # i = 1
#   i0 = i - 1
#   ID_i <- unique(early_grid_cohort_dnorm$ID)[i]
#   df_i <- filter(early_grid_cohort_dnorm, ID == ID_i)
#   set.seed(i)
#   vec_i <- round(dnorm(x = seq(unique(df_i$MinimumJd), unique(df_i$MaximumJd), 1),
#                        mean = unique(df_i$mean),
#                        sd = unique(df_i$sd))*10000)
#   # Get vec_i sum
#   # Correct vec_i to Site_Sum
#   site_sum_i <- unique(df_i$Site_Sum)
#   correction_factor <- site_sum_i/sum(vec_i)
#   vec_i <- vec_i * correction_factor
#   # Introduce some random variability
#   set.seed(i)
#   vec_i <- vec_i + rnorm(n = length(c(unique(df_i$MinimumJd):unique(df_i$MaximumJd))), mean = 0, sd = unique(df_i$sd_2))
#   # Correct negative values to zero
#   vec_i <- ifelse(vec_i < 0, 0, vec_i)
#   # Correct 1 values to 2
#   vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
#   # How many non-zero days are there
#   length_above_zero <- length(vec_i[vec_i > 0])
#   # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
#   # Index the minimum spot
#   min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
#   vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
#   # Get vec_i sum again
#   # How many sites is this?
#   correction_factor_2 <- site_sum_i / sum(vec_i)
#   # Multiply the vector by the correction factor again
#   vec_i <- vec_i * correction_factor_2
#   # Round the vector
#   vec_i <- round(vec_i)
#   # What is the absolute difference between the site numbers?
#   vec_diff <- site_sum_i - sum(vec_i)
#   # # Add/subtract this to a position at the max value position
#   # max_pos_i <- which.max(vec_i)
#   # vec_i[max_pos_i] <- vec_i[max_pos_i] + vec_diff
#   # Add/subtract this to a position at the max value positions
#   vec_i_df <- data.frame(vec_i = vec_i)
#   vec_i_df <- vec_i_df %>%
#     mutate(rowid = row_number())
# 
#   if(vec_diff < 0){
#     # Max rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_max(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Subtract from max rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i - 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else if(vec_diff > 0) {
#     # Min rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_min(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Add to min rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i + 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else {
#     vec_i <- vec_i
#   }
# 
#   # Add to Daily_Site_Total column (15)
#   start_ind <- 106*i0 + 1
#   end_ind <- 106*i0 + 1 + 105
#   early_grid_cohort_dnorm[c(start_ind:end_ind), 18] <- vec_i
# }
# # Plot
# early_grid_cohort_dnorm %>%
#   filter(Year == 1) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID)
# # Plot
# early_grid_cohort_dnorm %>%
#   filter(Year == 1) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID, scales = "free_y")
# # Plot
# early_grid_cohort_dnorm %>%
#   filter(Year == 1) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0)



# # middle_grid_cohort_pnorm
# # Create values for Daily_Site_Total
# head(middle_grid_cohort_pnorm)
# middle_grid_cohort_pnorm_grids <- middle_grid_cohort_pnorm %>%
#   group_by(GridID) %>%
#   mutate(Site_Sum_mean = mean(Site_Sum)) %>%
#   ungroup() %>%
#   select(GridID, Site_Sum_mean, Grid_Adj_mean, Grid_Adj_sd) %>%
#   distinct() %>%
#   arrange(-desc(Site_Sum_mean))
# set.seed(1)
# # stats
# middle_pnorm_grids_stats <- data.frame(mean = round(rnorm(n = nrow(middle_grid_cohort_pnorm_grids), mean = pnorm_middle_mean_mean, sd = pnorm_middle_mean_sd)),
#                                       sd = round(rnorm(n = nrow(middle_grid_cohort_pnorm_grids), mean = pnorm_middle_sd_mean, sd = pnorm_middle_sd_sd)))
# set.seed(1)
# middle_pnorm_grids_stats <- middle_pnorm_grids_stats %>%
#   arrange(desc(mean))
# # stats_2
# middle_pnorm_grids_stats_2 <- data.frame(mean_2 = round(rlnorm(n = nrow(middle_grid_cohort_pnorm_grids), meanlog = pnorm_middle_mean_mean_2, sdlog = pnorm_middle_mean_sd_2)),
#                                         sd_2 = round(rlnorm(n = nrow(middle_grid_cohort_pnorm_grids), meanlog = pnorm_middle_sd_mean_2, sdlog = pnorm_middle_sd_sd_2)))
# middle_pnorm_grids_stats_2 <- middle_pnorm_grids_stats_2 %>%
#   arrange(-desc(mean_2))
# # cbind
# middle_grid_cohort_pnorm_grids <- cbind(middle_grid_cohort_pnorm_grids, middle_pnorm_grids_stats, middle_pnorm_grids_stats_2)
# # Add grid_adj
# middle_grid_cohort_pnorm_grids <- middle_grid_cohort_pnorm_grids %>%
#   mutate(mean = mean + Grid_Adj_mean,
#          sd = sd + Grid_Adj_sd)
# 
# middle_grid_cohort_pnorm <- left_join(middle_grid_cohort_pnorm, middle_grid_cohort_pnorm_grids)
# # Create Daily_Site_Total column
# middle_grid_cohort_pnorm$Daily_Site_Total <- NA
# middle_grid_cohort_pnorm$Daily_Site_Total <- as.numeric(middle_grid_cohort_pnorm$Daily_Site_Total)
# # Loop
# for (i in 1:length(unique(middle_grid_cohort_pnorm$ID))) {
#   # i = 1
#   i0 = i - 1
#   ID_i <- unique(middle_grid_cohort_pnorm$ID)[i]
#   df_i <- filter(middle_grid_cohort_pnorm, ID == ID_i)
#   set.seed(i)
#   vec_i <- pnorm(q = seq(unique(df_i$MinimumJd), unique(df_i$MaximumJd), 1),
#                  mean = unique(df_i$mean),
#                  sd = unique(df_i$sd))*100
#   # Get vec_i sum
#   # Correct vec_i to Site_Sum
#   site_sum_i <- unique(df_i$Site_Sum)
#   correction_factor <- site_sum_i/sum(vec_i)
#   vec_i <- vec_i * correction_factor
#   # Introduce some random variability
#   set.seed(i)
#   vec_i <- vec_i + rnorm(n = length(c(unique(df_i$MinimumJd):unique(df_i$MaximumJd))), mean = 0, sd = unique(df_i$sd_2))
#   # Correct negative values to zero
#   vec_i <- ifelse(vec_i < 0, 0, vec_i)
#   # Correct 1 values to 2
#   vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
#   # How many non-zero days are there
#   length_above_zero <- length(vec_i[vec_i > 0])
#   # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
#   # Index the minimum spot
#   min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
#   vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
#   # Get vec_i sum again
#   # How many sites is this?
#   correction_factor_2 <- site_sum_i / sum(vec_i)
#   # Multiply the vector by the correction factor again
#   vec_i <- vec_i * correction_factor_2
#   # Round the vector
#   vec_i <- round(vec_i)
#   # What is the absolute difference between the site numbers?
#   vec_diff <- site_sum_i - sum(vec_i)
#   # # Add/subtract this to a position at the max value position
#   # max_pos_i <- which.max(vec_i)
#   # vec_i[max_pos_i] <- vec_i[max_pos_i] + vec_diff
#   # Add/subtract this to a position at the max value positions
#   vec_i_df <- data.frame(vec_i = vec_i)
#   vec_i_df <- vec_i_df %>%
#     mutate(rowid = row_number())
# 
#   if(vec_diff < 0){
#     # Max rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_max(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Subtract from max rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i - 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else if(vec_diff > 0) {
#     # Min rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_min(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Add to min rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i + 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else {
#     vec_i <- vec_i
#   }
# 
#   # Add to Daily_Site_Total column (15)
#   start_ind <- 106*i0 + 1
#   end_ind <- 106*i0 + 1 + 105
#   middle_grid_cohort_pnorm[c(start_ind:end_ind), 18] <- vec_i
# }
# # Plot
# middle_grid_cohort_pnorm %>%
#   filter(Year == 7) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID)
# # Plot
# middle_grid_cohort_pnorm %>%
#   filter(Year == 7) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID, scales = "free_y")
# # Plot
# middle_grid_cohort_pnorm %>%
#   filter(Year == 7) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0)

# # middle_grid_cohort_rnorm
# # Create values for Daily_Site_Total
# head(middle_grid_cohort_rnorm)
# middle_grid_cohort_rnorm_grids <- middle_grid_cohort_rnorm %>%
#   group_by(GridID) %>%
#   mutate(Site_Sum_mean = mean(Site_Sum)) %>%
#   ungroup() %>%
#   select(GridID, Site_Sum_mean, Grid_Adj_mean, Grid_Adj_sd) %>%
#   distinct() %>%
#   arrange(-desc(Site_Sum_mean))
# set.seed(1)
# # stats
# middle_rnorm_grids_stats <- data.frame(mean = round(rnorm(n = nrow(middle_grid_cohort_rnorm_grids), mean = rnorm_middle_mean_mean, sd = rnorm_middle_mean_sd)),
#                                        sd = round(rnorm(n = nrow(middle_grid_cohort_rnorm_grids), mean = rnorm_middle_sd_mean, sd = rnorm_middle_sd_sd)))
# set.seed(1)
# middle_rnorm_grids_stats <- middle_rnorm_grids_stats %>%
#   arrange(desc(mean))
# # stats_2
# middle_rnorm_grids_stats_2 <- data.frame(mean_2 = round(rlnorm(n = nrow(middle_grid_cohort_rnorm_grids), meanlog = rnorm_middle_mean_mean_2, sdlog = rnorm_middle_mean_sd_2)),
#                                          sd_2 = round(rlnorm(n = nrow(middle_grid_cohort_rnorm_grids), meanlog = rnorm_middle_sd_mean_2, sdlog = rnorm_middle_sd_sd_2)))
# middle_rnorm_grids_stats_2 <- middle_rnorm_grids_stats_2 %>%
#   arrange(-desc(mean_2))
# # cbind
# middle_grid_cohort_rnorm_grids <- cbind(middle_grid_cohort_rnorm_grids, middle_rnorm_grids_stats, middle_rnorm_grids_stats_2)
# # Add grid_adj
# middle_grid_cohort_rnorm_grids <- middle_grid_cohort_rnorm_grids %>%
#   mutate(mean = mean + Grid_Adj_mean,
#          sd = sd + Grid_Adj_sd)
# 
# middle_grid_cohort_rnorm <- left_join(middle_grid_cohort_rnorm, middle_grid_cohort_rnorm_grids)
# # Create Daily_Site_Total column
# middle_grid_cohort_rnorm$Daily_Site_Total <- NA
# middle_grid_cohort_rnorm$Daily_Site_Total <- as.numeric(middle_grid_cohort_rnorm$Daily_Site_Total)
# # Loop
# for (i in 1:length(unique(middle_grid_cohort_rnorm$ID))) {
#   # i = 1
#   i0 = i - 1
#   ID_i <- unique(middle_grid_cohort_rnorm$ID)[i]
#   df_i <- filter(middle_grid_cohort_rnorm, ID == ID_i)
#   set.seed(i)
#   vec_i <- round(rnorm(n = seq(unique(df_i$MinimumJd), unique(df_i$MaximumJd), 1),
#                        mean = unique(df_i$mean_2),
#                        sd = unique(df_i$sd_2)))
#   # Get vec_i sum
#   # Correct vec_i to Site_Sum
#   site_sum_i <- unique(df_i$Site_Sum)
#   correction_factor <- site_sum_i/sum(vec_i)
#   vec_i <- vec_i * correction_factor
#   # Introduce some random variability
#   set.seed(i)
#   vec_i <- vec_i + rnorm(n = length(c(unique(df_i$MinimumJd):unique(df_i$MaximumJd))), mean = 0, sd = unique(df_i$sd_2))
#   # Correct negative values to zero
#   vec_i <- ifelse(vec_i < 0, 0, vec_i)
#   # Correct 1 values to 2
#   vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
#   # How many non-zero days are there
#   length_above_zero <- length(vec_i[vec_i > 0])
#   # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
#   # Index the minimum spot
#   min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
#   vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
#   # Get vec_i sum again
#   # How many sites is this?
#   correction_factor_2 <- site_sum_i / sum(vec_i)
#   # Multiply the vector by the correction factor again
#   vec_i <- vec_i * correction_factor_2
#   # Round the vector
#   vec_i <- round(vec_i)
#   # What is the absolute difference between the site numbers?
#   vec_diff <- site_sum_i - sum(vec_i)
#   # Add/subtract this to a position at the max value positions
#   vec_i_df <- data.frame(vec_i = vec_i)
#   vec_i_df <- vec_i_df %>%
#     mutate(rowid = row_number())
# 
#   if(vec_diff < 0){
#     # Max rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_max(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Subtract from max rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i - 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else if(vec_diff > 0) {
#     # Min rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_min(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Add to min rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i + 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else {
#     vec_i <- vec_i
#   }
# 
#   # max_pos_i <- which.max(vec_i)
#   # vec_i[max_pos_i] <- vec_i[max_pos_i] + vec_diff
#   # Add to Daily_Site_Total column (15)
#   # Add to Daily_Site_Total column (15)
#   start_ind <- 106*i0 + 1
#   end_ind <- 106*i0 + 1 + 105
#   middle_grid_cohort_rnorm[c(start_ind:end_ind), 18] <- vec_i
# }
# # Plot
# middle_grid_cohort_rnorm %>%
#   filter(Year == 7) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID)
# # Plot
# middle_grid_cohort_rnorm %>%
#   filter(Year == 7) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID, scales = "free_y")
# # Plot
# middle_grid_cohort_rnorm %>%
#   filter(Year == 7) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0)

# middle_grid_cohort
# Create values for Daily_Site_Total
head(middle_grid_cohort)
middle_grid_cohort_grids <- middle_grid_cohort %>%
  group_by(GridID) %>%
  mutate(Site_Sum_mean = mean(Site_Sum)) %>%
  ungroup() %>%
  select(GridID, Site_Sum_mean, Grid_Adj_mean, Grid_Adj_sd) %>%
  distinct() %>%
  arrange(-desc(Site_Sum_mean))
set.seed(1)
# stats
middle_grids_stats <- data.frame(mean = round(rnorm(n = nrow(middle_grid_cohort_grids), mean = rnorm_middle_mean_mean, sd = rnorm_middle_mean_sd)),
                                       sd = round(rnorm(n = nrow(middle_grid_cohort_grids), mean = rnorm_middle_sd_mean, sd = rnorm_middle_sd_sd)))
set.seed(1)
middle_grids_stats <- middle_grids_stats %>%
  arrange(desc(mean))
# stats_2
middle_grids_stats_2 <- data.frame(mean_2 = round(rlnorm(n = nrow(middle_grid_cohort_grids), meanlog = rnorm_middle_mean_mean_2, sdlog = rnorm_middle_mean_sd_2)),
                                         sd_2 = round(rlnorm(n = nrow(middle_grid_cohort_grids), meanlog = rnorm_middle_sd_mean_2, sdlog = rnorm_middle_sd_sd_2)))
middle_grids_stats_2 <- middle_grids_stats_2 %>%
  arrange(-desc(mean_2))
# cbind
middle_grid_cohort_grids <- cbind(middle_grid_cohort_grids, middle_grids_stats, middle_grids_stats_2)
# Add grid_adj
middle_grid_cohort_grids <- middle_grid_cohort_grids %>%
  mutate(mean = mean + Grid_Adj_mean,
         sd = sd + Grid_Adj_sd)

middle_grid_cohort <- left_join(middle_grid_cohort, middle_grid_cohort_grids)
# Create Daily_Site_Total column
middle_grid_cohort$Daily_Site_Total <- NA
middle_grid_cohort$Daily_Site_Total <- as.numeric(middle_grid_cohort$Daily_Site_Total)
# Loop
for (i in 1:length(unique(middle_grid_cohort$ID))) {
  # i = 1
  i0 = i - 1
  ID_i <- unique(middle_grid_cohort$ID)[i]
  df_i <- filter(middle_grid_cohort, ID == ID_i)
  set.seed(i)
  vec_i <- round(rnorm(n = seq(unique(df_i$MinimumJd), unique(df_i$MaximumJd), 1),
                       mean = unique(df_i$mean_2),
                       sd = unique(df_i$sd_2)))
  # Get vec_i sum
  # Correct vec_i to Site_Sum
  site_sum_i <- unique(df_i$Site_Sum)
  correction_factor <- site_sum_i/sum(vec_i)
  vec_i <- vec_i * correction_factor
  # Introduce some random variability
  set.seed(i)
  vec_i <- vec_i + rnorm(n = length(c(unique(df_i$MinimumJd):unique(df_i$MaximumJd))), mean = 0, sd = unique(df_i$sd_2))
  # Correct negative values to zero
  vec_i <- ifelse(vec_i < 0, 0, vec_i)
  # Correct 1 values to 2
  vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
  # How many non-zero days are there
  length_above_zero <- length(vec_i[vec_i > 0])
  # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
  # Index the minimum spot
  min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
  vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
  # Get vec_i sum again
  # How many sites is this?
  correction_factor_2 <- site_sum_i / sum(vec_i)
  # Multiply the vector by the correction factor again
  vec_i <- vec_i * correction_factor_2
  # Round the vector
  vec_i <- round(vec_i)
  # What is the absolute difference between the site numbers?
  vec_diff <- site_sum_i - sum(vec_i)
  # Add/subtract this to a position at the max value positions
  vec_i_df <- data.frame(vec_i = vec_i)
  vec_i_df <- vec_i_df %>%
    mutate(rowid = row_number())
  
  if(vec_diff < 0){
    # Max rows
    vec_i_df_subset <- vec_i_df %>%
      slice_max(vec_i, n = abs(vec_diff), with_ties = FALSE)
    # Other rows
    vec_i_df_not_subset <- vec_i_df %>%
      filter(!rowid %in% vec_i_df_subset$rowid)
    # Subtract from max rows
    vec_i_df_subset <- vec_i_df_subset %>%
      mutate(vec_i = vec_i - 1)
    # Bind data.frames
    vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
    # Arrange data.frame
    vec_i_df <- vec_i_df %>%
      arrange(-desc(rowid))
    # Get vector
    vec_i <- vec_i_df$vec_i
    
  } else if(vec_diff > 0) {
    # Min rows
    vec_i_df_subset <- vec_i_df %>%
      slice_min(vec_i, n = abs(vec_diff), with_ties = FALSE)
    # Other rows
    vec_i_df_not_subset <- vec_i_df %>%
      filter(!rowid %in% vec_i_df_subset$rowid)
    # Add to min rows
    vec_i_df_subset <- vec_i_df_subset %>%
      mutate(vec_i = vec_i + 1)
    # Bind data.frames
    vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
    # Arrange data.frame
    vec_i_df <- vec_i_df %>%
      arrange(-desc(rowid))
    # Get vector
    vec_i <- vec_i_df$vec_i
    
  } else {
    vec_i <- vec_i
  }
  
  # max_pos_i <- which.max(vec_i)
  # vec_i[max_pos_i] <- vec_i[max_pos_i] + vec_diff
  # Add to Daily_Site_Total column (15)
  # Add to Daily_Site_Total column (15)
  start_ind <- 106*i0 + 1
  end_ind <- 106*i0 + 1 + 105
  middle_grid_cohort[c(start_ind:end_ind), 16] <- vec_i
}
# Plot
middle_grid_cohort %>%
  filter(Year == 7) %>%
  ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
  geom_line() +
  expand_limits(y = 0) +
  facet_wrap(~GridID)
# Plot
middle_grid_cohort %>%
  filter(Year == 7) %>%
  ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
  geom_line() +
  expand_limits(y = 0) +
  facet_wrap(~GridID, scales = "free_y")
# Plot
middle_grid_cohort %>%
  filter(Year == 7) %>%
  ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
  geom_line() +
  expand_limits(y = 0)

# # middle_grid_cohort_dnorm
# # Create values for Daily_Site_Total
# head(middle_grid_cohort_dnorm)
# middle_grid_cohort_dnorm_grids <- middle_grid_cohort_dnorm %>%
#   group_by(GridID) %>%
#   mutate(Site_Sum_mean = mean(Site_Sum)) %>%
#   ungroup() %>%
#   select(GridID, Site_Sum_mean, Grid_Adj_mean, Grid_Adj_sd) %>%
#   distinct() %>%
#   arrange(-desc(Site_Sum_mean))
# set.seed(1)
# # stats
# middle_dnorm_grids_stats <- data.frame(mean = round(rnorm(n = nrow(middle_grid_cohort_dnorm_grids), mean = dnorm_middle_mean_mean, sd = dnorm_middle_mean_sd)),
#                                        sd = round(rnorm(n = nrow(middle_grid_cohort_dnorm_grids), mean = dnorm_middle_sd_mean, sd = dnorm_middle_sd_sd)))
# set.seed(1)
# middle_dnorm_grids_stats <- middle_dnorm_grids_stats %>%
#   arrange(desc(mean))
# # stats_2
# middle_dnorm_grids_stats_2 <- data.frame(mean_2 = round(rlnorm(n = nrow(middle_grid_cohort_dnorm_grids), meanlog = dnorm_middle_mean_mean_2, sdlog = dnorm_middle_mean_sd_2)),
#                                          sd_2 = round(rlnorm(n = nrow(middle_grid_cohort_dnorm_grids), meanlog = dnorm_middle_sd_mean_2, sdlog = dnorm_middle_sd_sd_2)))
# middle_dnorm_grids_stats_2 <- middle_dnorm_grids_stats_2 %>%
#   arrange(-desc(mean_2))
# # cbind
# middle_grid_cohort_dnorm_grids <- cbind(middle_grid_cohort_dnorm_grids, middle_dnorm_grids_stats, middle_dnorm_grids_stats_2)
# # Add grid_adj
# middle_grid_cohort_dnorm_grids <- middle_grid_cohort_dnorm_grids %>%
#   mutate(mean = mean + Grid_Adj_mean,
#          sd = sd + Grid_Adj_sd)
# 
# middle_grid_cohort_dnorm <- left_join(middle_grid_cohort_dnorm, middle_grid_cohort_dnorm_grids)
# # Create Daily_Site_Total column
# middle_grid_cohort_dnorm$Daily_Site_Total <- NA
# middle_grid_cohort_dnorm$Daily_Site_Total <- as.numeric(middle_grid_cohort_dnorm$Daily_Site_Total)
# # Loop
# for (i in 1:length(unique(middle_grid_cohort_dnorm$ID))) {
#   # i = 1
#   i0 = i - 1
#   ID_i <- unique(middle_grid_cohort_dnorm$ID)[i]
#   df_i <- filter(middle_grid_cohort_dnorm, ID == ID_i)
#   set.seed(i)
#   vec_i <- round(dnorm(x = seq(unique(df_i$MinimumJd), unique(df_i$MaximumJd), 1),
#                        mean = unique(df_i$mean),
#                        sd = unique(df_i$sd))*3000)
#   # Get vec_i sum
#   # Correct vec_i to Site_Sum
#   site_sum_i <- unique(df_i$Site_Sum)
#   correction_factor <- site_sum_i/sum(vec_i)
#   vec_i <- vec_i * correction_factor
#   # Introduce some random variability
#   set.seed(i)
#   vec_i <- vec_i + rnorm(n = length(c(unique(df_i$MinimumJd):unique(df_i$MaximumJd))), mean = 0, sd = unique(df_i$sd_2))
#   # Correct negative values to zero
#   vec_i <- ifelse(vec_i < 0, 0, vec_i)
#   # Correct 1 values to 2
#   vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
#   # How many non-zero days are there
#   length_above_zero <- length(vec_i[vec_i > 0])
#   # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
#   # Index the minimum spot
#   min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
#   vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
#   # Get vec_i sum again
#   # How many sites is this?
#   correction_factor_2 <- site_sum_i / sum(vec_i)
#   # Multiply the vector by the correction factor again
#   vec_i <- vec_i * correction_factor_2
#   # Round the vector
#   vec_i <- round(vec_i)
#   # What is the absolute difference between the site numbers?
#   vec_diff <- site_sum_i - sum(vec_i)
#   # # Add/subtract this to a position at the max value position
#   # max_pos_i <- which.max(vec_i)
#   # vec_i[max_pos_i] <- vec_i[max_pos_i] + vec_diff
#   # Add/subtract this to a position at the max value positions
#   vec_i_df <- data.frame(vec_i = vec_i)
#   vec_i_df <- vec_i_df %>%
#     mutate(rowid = row_number())
# 
#   if(vec_diff < 0){
#     # Max rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_max(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Subtract from max rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i - 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else if(vec_diff > 0) {
#     # Min rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_min(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Add to min rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i + 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else {
#     vec_i <- vec_i
#   }
# 
#   # Add to Daily_Site_Total column (15)
#   start_ind <- 106*i0 + 1
#   end_ind <- 106*i0 + 1 + 105
#   middle_grid_cohort_dnorm[c(start_ind:end_ind), 18] <- vec_i
# }
# # Plot
# middle_grid_cohort_dnorm %>%
#   filter(Year == 7) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID)
# # Plot
# middle_grid_cohort_dnorm %>%
#   filter(Year == 7) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID, scales = "free_y")
# # Plot
# middle_grid_cohort_dnorm %>%
#   filter(Year == 7) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0)



# # late_grid_cohort_dnorm
# # Create values for Daily_Site_Total
# head(late_grid_cohort_dnorm)
# late_grid_cohort_dnorm_grids <- late_grid_cohort_dnorm %>%
#   group_by(GridID) %>%
#   mutate(Site_Sum_mean = mean(Site_Sum)) %>%
#   ungroup() %>%
#   select(GridID, Site_Sum_mean, Grid_Adj_mean, Grid_Adj_sd) %>%
#   distinct() %>%
#   arrange(-desc(Site_Sum_mean))
# set.seed(1)
# # stats
# late_dnorm_grids_stats <- data.frame(mean = round(rnorm(n = nrow(late_grid_cohort_dnorm_grids), mean = late_mean_mean, sd = late_mean_sd)),
#                                        sd = round(rnorm(n = nrow(late_grid_cohort_dnorm_grids), mean = late_sd_mean, sd = late_sd_sd)))
# set.seed(1)
# late_dnorm_grids_stats <- late_dnorm_grids_stats %>%
#   arrange(desc(mean))
# # stats_2
# late_dnorm_grids_stats_2 <- data.frame(mean_2 = round(rlnorm(n = nrow(late_grid_cohort_dnorm_grids), meanlog = late_mean_mean_2, sdlog = late_mean_sd_2)),
#                                          sd_2 = round(rlnorm(n = nrow(late_grid_cohort_dnorm_grids), meanlog = late_sd_mean_2, sdlog = late_sd_sd_2)))
# late_dnorm_grids_stats_2 <- late_dnorm_grids_stats_2 %>%
#   arrange(-desc(mean_2))
# # cbind
# late_grid_cohort_dnorm_grids <- cbind(late_grid_cohort_dnorm_grids, late_dnorm_grids_stats, late_dnorm_grids_stats_2)
# # Add grid_adj
# late_grid_cohort_dnorm_grids <- late_grid_cohort_dnorm_grids %>%
#   mutate(mean = mean + Grid_Adj_mean,
#          sd = sd + Grid_Adj_sd)
# 
# late_grid_cohort_dnorm <- left_join(late_grid_cohort_dnorm, late_grid_cohort_dnorm_grids)
# # Create Daily_Site_Total column
# late_grid_cohort_dnorm$Daily_Site_Total <- NA
# late_grid_cohort_dnorm$Daily_Site_Total <- as.numeric(late_grid_cohort_dnorm$Daily_Site_Total)
# # Loop
# for (i in 1:length(unique(late_grid_cohort_dnorm$ID))) {
#   # i = 1
#   i0 = i - 1
#   ID_i <- unique(late_grid_cohort_dnorm$ID)[i]
#   df_i <- filter(late_grid_cohort_dnorm, ID == ID_i)
#   set.seed(i)
#   vec_i <- round(dnorm(x = seq(unique(df_i$MinimumJd), unique(df_i$MaximumJd), 1),
#                        mean = unique(df_i$mean),
#                        sd = unique(df_i$sd))*3000)
#   # Get vec_i sum
#   # Correct vec_i to Site_Sum
#   site_sum_i <- unique(df_i$Site_Sum)
#   correction_factor <- site_sum_i/sum(vec_i)
#   vec_i <- vec_i * correction_factor
#   # Introduce some random variability
#   set.seed(i)
#   vec_i <- vec_i + rnorm(n = length(c(unique(df_i$MinimumJd):unique(df_i$MaximumJd))), mean = 0, sd = unique(df_i$sd_2))
#   # Correct negative values to zero
#   vec_i <- ifelse(vec_i < 0, 0, vec_i)
#   # Correct 1 values to 2
#   vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
#   # How many non-zero days are there
#   length_above_zero <- length(vec_i[vec_i > 0])
#   # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
#   # Index the minimum spot
#   min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
#   vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
#   # Get vec_i sum again
#   # How many sites is this?
#   correction_factor_2 <- site_sum_i / sum(vec_i)
#   # Multiply the vector by the correction factor again
#   vec_i <- vec_i * correction_factor_2
#   # Round the vector
#   vec_i <- round(vec_i)
#   # What is the absolute difference between the site numbers?
#   vec_diff <- site_sum_i - sum(vec_i)
#   # # Add/subtract this to a position at the max value position
#   # max_pos_i <- which.max(vec_i)
#   # vec_i[max_pos_i] <- vec_i[max_pos_i] + vec_diff
#   # Add/subtract this to a position at the max value positions
#   vec_i_df <- data.frame(vec_i = vec_i)
#   vec_i_df <- vec_i_df %>%
#     mutate(rowid = row_number())
# 
#   if(vec_diff < 0){
#     # Max rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_max(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Subtract from max rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i - 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else if(vec_diff > 0) {
#     # Min rows
#     vec_i_df_subset <- vec_i_df %>%
#       slice_min(vec_i, n = abs(vec_diff), with_ties = FALSE)
#     # Other rows
#     vec_i_df_not_subset <- vec_i_df %>%
#       filter(!rowid %in% vec_i_df_subset$rowid)
#     # Add to min rows
#     vec_i_df_subset <- vec_i_df_subset %>%
#       mutate(vec_i = vec_i + 1)
#     # Bind data.frames
#     vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
#     # Arrange data.frame
#     vec_i_df <- vec_i_df %>%
#       arrange(-desc(rowid))
#     # Get vector
#     vec_i <- vec_i_df$vec_i
# 
#   } else {
#     vec_i <- vec_i
#   }
# 
#   # Add to Daily_Site_Total column (15)
#   start_ind <- 106*i0 + 1
#   end_ind <- 106*i0 + 1 + 105
#   late_grid_cohort_dnorm[c(start_ind:end_ind), 18] <- vec_i
# }
# # Plot
# late_grid_cohort_dnorm %>%
#   filter(Year == 13) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID)
# # Plot
# late_grid_cohort_dnorm %>%
#   filter(Year == 13) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0) +
#   facet_wrap(~GridID, scales = "free_y")
# # Plot
# late_grid_cohort_dnorm %>%
#   filter(Year == 13) %>%
#   ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
#   geom_line() +
#   expand_limits(y = 0)

# late_grid_cohort
# Create values for Daily_Site_Total
head(late_grid_cohort)
late_grid_cohort_grids <- late_grid_cohort %>%
  group_by(GridID) %>%
  mutate(Site_Sum_mean = mean(Site_Sum)) %>%
  ungroup() %>%
  select(GridID, Site_Sum_mean, Grid_Adj_mean, Grid_Adj_sd) %>%
  distinct() %>%
  arrange(-desc(Site_Sum_mean))
set.seed(1)
# stats
late_grids_stats <- data.frame(mean = round(rnorm(n = nrow(late_grid_cohort_grids), mean = late_mean_mean, sd = late_mean_sd)),
                                 sd = round(rnorm(n = nrow(late_grid_cohort_grids), mean = late_sd_mean, sd = late_sd_sd)))
set.seed(1)
late_grids_stats <- late_grids_stats %>%
  arrange(desc(mean))
# stats_2
late_grids_stats_2 <- data.frame(mean_2 = round(rlnorm(n = nrow(late_grid_cohort_grids), meanlog = late_mean_mean_2, sdlog = late_mean_sd_2)),
                                   sd_2 = round(rlnorm(n = nrow(late_grid_cohort_grids), meanlog = late_sd_mean_2, sdlog = late_sd_sd_2)))
late_grids_stats_2 <- late_grids_stats_2 %>%
  arrange(-desc(mean_2))
# cbind
late_grid_cohort_grids <- cbind(late_grid_cohort_grids, late_grids_stats, late_grids_stats_2)
# Add grid_adj
late_grid_cohort_grids <- late_grid_cohort_grids %>%
  mutate(mean = mean + Grid_Adj_mean,
         sd = sd + Grid_Adj_sd)

late_grid_cohort <- left_join(late_grid_cohort, late_grid_cohort_grids)
# Create Daily_Site_Total column
late_grid_cohort$Daily_Site_Total <- NA
late_grid_cohort$Daily_Site_Total <- as.numeric(late_grid_cohort$Daily_Site_Total)
# Loop
for (i in 1:length(unique(late_grid_cohort$ID))) {
  # i = 1
  i0 = i - 1
  ID_i <- unique(late_grid_cohort$ID)[i]
  df_i <- filter(late_grid_cohort, ID == ID_i)
  set.seed(i)
  vec_i <- round(rnorm(n = seq(unique(df_i$MinimumJd), unique(df_i$MaximumJd), 1),
                       mean = unique(df_i$mean_2),
                       sd = unique(df_i$sd_2)))
  # Get vec_i sum
  # Correct vec_i to Site_Sum
  site_sum_i <- unique(df_i$Site_Sum)
  correction_factor <- site_sum_i/sum(vec_i)
  vec_i <- vec_i * correction_factor
  # Introduce some random variability
  set.seed(i)
  vec_i <- vec_i + rnorm(n = length(c(unique(df_i$MinimumJd):unique(df_i$MaximumJd))), mean = 0, sd = unique(df_i$sd_2))
  # Correct negative values to zero
  vec_i <- ifelse(vec_i < 0, 0, vec_i)
  # Correct 1 values to 2
  vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
  # How many non-zero days are there
  length_above_zero <- length(vec_i[vec_i > 0])
  # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
  # Index the minimum spot
  min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
  vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
  # Get vec_i sum again
  # How many sites is this?
  correction_factor_2 <- site_sum_i / sum(vec_i)
  # Multiply the vector by the correction factor again
  vec_i <- vec_i * correction_factor_2
  # Round the vector
  vec_i <- round(vec_i)
  # What is the absolute difference between the site numbers?
  vec_diff <- site_sum_i - sum(vec_i)
  # Add/subtract this to a position at the max value positions
  vec_i_df <- data.frame(vec_i = vec_i)
  vec_i_df <- vec_i_df %>%
    mutate(rowid = row_number())
  
  if(vec_diff < 0){
    # Max rows
    vec_i_df_subset <- vec_i_df %>%
      slice_max(vec_i, n = abs(vec_diff), with_ties = FALSE)
    # Other rows
    vec_i_df_not_subset <- vec_i_df %>%
      filter(!rowid %in% vec_i_df_subset$rowid)
    # Subtract from max rows
    vec_i_df_subset <- vec_i_df_subset %>%
      mutate(vec_i = vec_i - 1)
    # Bind data.frames
    vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
    # Arrange data.frame
    vec_i_df <- vec_i_df %>%
      arrange(-desc(rowid))
    # Get vector
    vec_i <- vec_i_df$vec_i
    
  } else if(vec_diff > 0) {
    # Min rows
    vec_i_df_subset <- vec_i_df %>%
      slice_min(vec_i, n = abs(vec_diff), with_ties = FALSE)
    # Other rows
    vec_i_df_not_subset <- vec_i_df %>%
      filter(!rowid %in% vec_i_df_subset$rowid)
    # Add to min rows
    vec_i_df_subset <- vec_i_df_subset %>%
      mutate(vec_i = vec_i + 1)
    # Bind data.frames
    vec_i_df <- rbind(vec_i_df_subset, vec_i_df_not_subset)
    # Arrange data.frame
    vec_i_df <- vec_i_df %>%
      arrange(-desc(rowid))
    # Get vector
    vec_i <- vec_i_df$vec_i
    
  } else {
    vec_i <- vec_i
  }
  
  # max_pos_i <- which.max(vec_i)
  # vec_i[max_pos_i] <- vec_i[max_pos_i] + vec_diff
  # Add to Daily_Site_Total column (15)
  # Add to Daily_Site_Total column (15)
  start_ind <- 106*i0 + 1
  end_ind <- 106*i0 + 1 + 105
  late_grid_cohort[c(start_ind:end_ind), 16] <- vec_i
}
# Plot
late_grid_cohort %>%
  filter(Year == 13) %>%
  ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
  geom_line() +
  expand_limits(y = 0) +
  facet_wrap(~GridID)
# Plot
late_grid_cohort %>%
  filter(Year == 13) %>%
  ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
  geom_line() +
  expand_limits(y = 0) +
  facet_wrap(~GridID, scales = "free_y")
# Plot
late_grid_cohort %>%
  filter(Year == 13) %>%
  ggplot(aes(x = Jd, y = Daily_Site_Total, group = GridID)) +
  geom_line() +
  expand_limits(y = 0)



# -----------------------------------------------------------
# Join the grids with Daily_Site_Total to Blocks_2
# -----------------------------------------------------------
# Bind data.frames
# Daily_Site_Total_rbind <- rbind(early_grid_cohort_pnorm,
#                                 early_grid_cohort_rnorm,
#                                 early_grid_cohort_dnorm,
#                                 middle_grid_cohort_pnorm,
#                                 middle_grid_cohort_rnorm,
#                                 middle_grid_cohort_dnorm,
#                                 late_grid_cohort_dnorm)
Daily_Site_Total_rbind <- rbind(early_grid_cohort,
                                middle_grid_cohort,
                                late_grid_cohort)
# Select needed columns
Daily_Site_Total_rbind <- Daily_Site_Total_rbind %>%
  select(Jd, Year, GridID, Daily_Site_Total)
# Join to Blocks_2
Blocks_2 <- left_join(Blocks_2, Daily_Site_Total_rbind)



# ----------------------------------------------------------
# Create Daily_Site_Sightings data
# ----------------------------------------------------------
# ----------------------------------------------------------
# "Stationary" dataset
# ----------------------------------------------------------
# In this dataset, the mean values for the distributions for each species will be drawn randomly from a normal distribution, with standard deviations similar to those from species in the raw data

# If MaxJd == 180, then pnorm
# If MaxJd < 180, then dnorm

# # Create function that makes a new vector based on correlations with another one
# complement <- function(y, rho, x) {
#   if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
#   y.perp <- residuals(lm(x ~ y))
#   rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
# }

# y <- rnorm(50, 10)
# x <- 1:50 # optional
# rho <- seq(0, 1, length.out=6) * rep(c(-1,1), 3)
#
# X <- data.frame(z=as.vector(sapply(rho, function(rho) complement(y, rho))),
#                 rho=ordered(rep(signif(rho, 2), each=length(y))),
#                 y=rep(y, length(rho)))
#
# ggplot(X, aes(y,z, group=rho)) +
#   geom_smooth(method="lm", color="Black") +
#   geom_rug(sides="b") +
#   geom_point(aes(fill=rho), alpha=1/2, shape=21) +
#   facet_wrap(~ rho, scales="free")

# # Try out complement
# test <- Blocks_2 %>%
#   filter(ID == "Species_1_Grid_57_18")
# test_comp <- complement(y = test$Daily_Site_Total,
#                         rho = 0.2,
#                         x = pnorm(q = seq(81, 180, 1),
#                                   mean = 130,
#                                   sd = 5))
# test_df <- data.frame(Daily_Site_Total = test$Daily_Site_Total,
#                       Daily_Site_Sightings = test_comp,
#                       Jd = c(81:180))
# # plot
# ggplot(test_df) +
#   geom_line(aes(Jd, Daily_Site_Total), colour = "blue") +
#   geom_line(aes(Jd, Daily_Site_Sightings), colour = "red")
# # plot
# ggplot(test_df) +
#   geom_point(aes(Daily_Site_Total, Daily_Site_Sightings))






# ---------------------------------------------
# Model expectations
# ---------------------------------------------
# Get species coefficients function
# Stationary
coefficients_Stationary_1 <- function (commonname){
  Grid_Years_Species <- Grid_Years_Species[Grid_Years_Species$Common_Name == commonname ,]
  regression <- lmer(MAD_Stationary ~ Year + (1|GridID), data = Grid_Years_Species)
  MAD.Shift <- data.frame(fixef(regression))
  Common.name <- unique(Grid_Years_Species$Common_Name)
  data.frame(Common.name, MAD.Shift[2,])
}

x <- unique(Grid_Years_Species$Common_Name)
lmer_species_Stationary <- lapply(x, coefficients_Stationary_1)
lmer_species_Stationary <- data.frame(do.call(rbind, lmer_species_Stationary))
colnames(lmer_species_Stationary) <- c("Common_Name", "MAD_Shift")
# Plot
ggplot(lmer_species_Stationary, aes(MAD_Shift)) +
  geom_histogram(binwidth = 0.02, colour = "black") +
  geom_vline(xintercept = median(lmer_species_Stationary$MAD_Shift), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "solid") +
  labs(x = "Mean arrival date shift (days/year)",
       y = "# of species-grid cells") +
  theme_bw()
summary(lmer_species_Stationary$MAD_Shift)
sd(lmer_species_Stationary$MAD_Shift)

# Shifting
coefficients_Shifting_1 <- function (commonname){
  Grid_Years_Species <- Grid_Years_Species[Grid_Years_Species$Common_Name == commonname ,]
  regression <- lmer(MAD_Shifting ~ Year + (1|GridID), data = Grid_Years_Species)
  MAD.Shift <- data.frame(fixef(regression))
  Common.name <- unique(Grid_Years_Species$Common_Name)
  data.frame(Common.name, MAD.Shift[2,])
}

x <- unique(Grid_Years_Species$Common_Name)
lmer_species_Shifting <- lapply(x, coefficients_Shifting_1)
lmer_species_Shifting <- data.frame(do.call(rbind, lmer_species_Shifting))
colnames(lmer_species_Shifting) <- c("Common_Name", "MAD_Shift")
# Plot
ggplot(lmer_species_Shifting, aes(MAD_Shift)) +
  geom_histogram(binwidth = 0.02, colour = "black") +
  geom_vline(xintercept = median(lmer_species_Shifting$MAD_Shift), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "solid") +
  labs(x = "Mean arrival date shift (days/year)",
       y = "# of species-grid cells") +
  theme_bw()
summary(lmer_species_Shifting$MAD_Shift)
sd(lmer_species_Shifting$MAD_Shift)



# ---------------------------------------------
# Get Daily_Site_Sightings
# ---------------------------------------------
# Try plotting
# early
ggplot(Daily_Totals_mean_early, aes(x = Jd, y = Daily_Site_Sightings_mean)) +
  geom_line() +
  geom_smooth(se = FALSE) +
  facet_wrap(vars(Grid_rank), scales = "free")

# late
ggplot(Daily_Totals_mean_late, aes(x = Jd, y = Daily_Site_Sightings_mean)) +
  geom_line() +
  geom_smooth(se = FALSE) +
  facet_wrap(vars(Grid_rank), scales = "free")

# Determine the means
# early
Daily_Site_Sightings_uncount_early <- Daily_Totals_mean_early %>%
  uncount(Daily_Site_Sightings_mean)
Daily_Site_Sightings_uncount_early_stats <- Daily_Site_Sightings_uncount_early %>%
  group_by(GridID) %>%
  summarize(mean_Jd = mean(Jd),
            sd_Jd = sd(Jd)) %>%
  ungroup()
summary(select(Daily_Site_Sightings_uncount_early_stats, mean_Jd, sd_Jd))
sd(Daily_Site_Sightings_uncount_early_stats$mean_Jd)
sd(Daily_Site_Sightings_uncount_early_stats$sd_Jd)

# Daily_Site_Total_mean (log)
Daily_Site_Sightings_uncount_early_stats <- left_join(Daily_Site_Sightings_uncount_early_stats,
                                                Daily_Totals_mean_early %>%
                                                  group_by(GridID) %>%
                                                  summarize(mean_log_Daily_Site_Sightings_mean = mean(log(Daily_Site_Sightings_mean + 0.001)),
                                                            sd_log_Daily_Site_Sightings_mean = sd(log(Daily_Site_Sightings_mean + 0.001))) %>%
                                                          ungroup())

summary(select(Daily_Site_Sightings_uncount_early_stats, mean_Jd, sd_Jd, mean_log_Daily_Site_Sightings_mean, sd_log_Daily_Site_Sightings_mean))
sd(Daily_Site_Sightings_uncount_early_stats$mean_Jd)
sd(Daily_Site_Sightings_uncount_early_stats$sd_Jd)

# Define the means/sds
# Jd
early_mean_mean_sightings <- round(mean(Daily_Site_Sightings_uncount_early_stats$mean_Jd), digits = 1)
early_mean_sd_sightings <- round(sd(Daily_Site_Sightings_uncount_early_stats$mean_Jd), digits = 1)
early_sd_mean_sightings <- round(mean(Daily_Site_Sightings_uncount_early_stats$sd_Jd), digits = 1)
early_sd_sd_sightings <- round(sd(Daily_Site_Sightings_uncount_early_stats$sd_Jd), digits = 1)
# log_Daily_Site_Sightings_mean
early_mean_mean_sightings_2 <- round(mean(Daily_Site_Sightings_uncount_early_stats$mean_log_Daily_Site_Sightings_mean), digits = 1)
early_mean_sd_sightings_2 <- round(sd(Daily_Site_Sightings_uncount_early_stats$mean_log_Daily_Site_Sightings_mean), digits = 1)
early_sd_mean_sightings_2 <- round(mean(Daily_Site_Sightings_uncount_early_stats$sd_log_Daily_Site_Sightings_mean), digits = 1)
early_sd_sd_sightings_2 <- round(sd(Daily_Site_Sightings_uncount_early_stats$sd_log_Daily_Site_Sightings_mean), digits = 1)

# late
Daily_Site_Sightings_uncount_late <- Daily_Totals_mean_late %>%
  uncount(Daily_Site_Sightings_mean)
# Jd
Daily_Site_Sightings_uncount_late_stats <- Daily_Site_Sightings_uncount_late %>%
  group_by(GridID) %>%
  summarize(mean_Jd = mean(Jd),
            sd_Jd = sd(Jd)) %>%
  ungroup()
summary(select(Daily_Site_Sightings_uncount_late_stats, mean_Jd, sd_Jd))
sd(Daily_Site_Sightings_uncount_late_stats$mean_Jd)
sd(Daily_Site_Sightings_uncount_late_stats$sd_Jd)
# Daily_Site_Total_mean (log)
Daily_Site_Sightings_uncount_late_stats <- left_join(Daily_Site_Sightings_uncount_late_stats,
                                                Daily_Totals_mean_late %>%
                                                  group_by(GridID) %>%
                                                  summarize(mean_log_Daily_Site_Sightings_mean = mean(log(Daily_Site_Sightings_mean + 0.001)),
                                                            sd_log_Daily_Site_Sightings_mean = sd(log(Daily_Site_Sightings_mean + 0.001))) %>%
                                                  ungroup())
summary(select(Daily_Site_Sightings_uncount_late_stats, mean_Jd, sd_Jd, mean_log_Daily_Site_Sightings_mean, sd_log_Daily_Site_Sightings_mean))
sd(Daily_Site_Sightings_uncount_late_stats$mean_Jd)
sd(Daily_Site_Sightings_uncount_late_stats$sd_Jd)

# Define the means/sds
# late
# Jd
late_mean_mean_sightings <- round(mean(Daily_Site_Sightings_uncount_late_stats$mean_Jd), digits = 1)
late_mean_sd_sightings <- round(sd(Daily_Site_Sightings_uncount_late_stats$mean_Jd), digits = 1)
late_sd_mean_sightings <- round(mean(Daily_Site_Sightings_uncount_late_stats$sd_Jd), digits = 1)
late_sd_sd_sightings <- round(sd(Daily_Site_Sightings_uncount_late_stats$sd_Jd), digits = 1)
# Daily_Site_Total_mean (log)
late_mean_mean_sightings_2 <- round(mean(Daily_Site_Sightings_uncount_late_stats$mean_log_Daily_Site_Sightings_mean), digits = 1)
late_mean_sd_sightings_2 <- round(sd(Daily_Site_Sightings_uncount_late_stats$mean_log_Daily_Site_Sightings_mean), digits = 1)
late_sd_mean_sightings_2 <- round(mean(Daily_Site_Sightings_uncount_late_stats$sd_log_Daily_Site_Sightings_mean), digits = 1)
late_sd_sd_sightings_2 <- round(sd(Daily_Site_Sightings_uncount_late_stats$sd_log_Daily_Site_Sightings_mean), digits = 1)

# middle
# Jd
middle_mean_mean_sightings <- round((early_mean_mean_sightings + late_mean_mean_sightings) / 2, digits = 1)
middle_mean_sd_sightings <- round((early_mean_sd_sightings + late_mean_sd_sightings) / 2, digits = 1)
middle_sd_mean_sightings <- round((early_sd_mean_sightings + late_sd_mean_sightings) / 2, digits = 1)
middle_sd_sd_sightings <- round((early_sd_sd_sightings + late_sd_sd_sightings) / 2, digits = 1)
# Daily_Site_Total_mean (log)
middle_mean_mean_sightings_2 <- round((early_mean_mean_sightings_2 + late_mean_mean_sightings_2) / 2, digits = 1)
middle_mean_sd_sightings_2 <- round((early_mean_sd_sightings_2 + late_mean_sd_sightings_2) / 2, digits = 1)
middle_sd_mean_sightings_2 <- round((early_sd_mean_sightings_2 + late_sd_mean_sightings_2) / 2, digits = 1)
middle_sd_sd_sightings_2 <- round((early_sd_sd_sightings_2 + late_sd_sd_sightings_2) / 2, digits = 1)

# Plot the stats
ggplot(Daily_Site_Sightings_uncount_late_stats, aes(mean_Jd)) +
  geom_histogram(binwidth = 1)
ggplot(Daily_Site_Sightings_uncount_late_stats, aes(sd_Jd)) +
  geom_histogram(binwidth = 0.25)

# What is the general correlation among Daily_Site_Total and Daily_Site_Sightings?
Daily_Stats_cor <- eBird_Daily_Occupancy_all_grids %>%
  group_by(Common_Name) %>%
  summarize(Daily_cor = cor(Daily_Site_Sightings, Daily_Site_Total)) %>%
  ungroup()
summary(Daily_Stats_cor$Daily_cor)
sd(Daily_Stats_cor$Daily_cor)
Daily_Stats_cor_mean <- mean(Daily_Stats_cor$Daily_cor)
Daily_Stats_cor_sd <- sd(Daily_Stats_cor$Daily_cor)

ggplot(Daily_Stats_cor, aes(Daily_cor)) +
  geom_histogram(binwidth = 0.05)

# Get correlation coefficients per species
set.seed(1)
Species_df_cor <- data.frame(Common_Name = unique(Species_df$Common_Name),
                             daily_cor = round(rnorm(n = 44, mean = Daily_Stats_cor_mean, sd = Daily_Stats_cor_sd), digits = 3))
# Join to main data.frame
Blocks_2 <- left_join(Blocks_2, Species_df_cor)



# --------------------------------------------------
# Add sd for daily site sightings
# --------------------------------------------------
Grid_Years_Species_select <- Grid_Years_Species %>%
  select(Year, GridID, Common_Name) %>%
  distinct()
set.seed(1)
Grid_Years_Species_select <- Grid_Years_Species_select %>%
  mutate(sd_2_sightings = case_when(
    Year %in% c(1:6) ~ (rlnorm(n = n(), meanlog = early_mean_sd_sightings_2, sdlog = early_sd_sd_sightings_2)),
    Year %in% c(7:12) ~ (rlnorm(n = n(), meanlog = middle_mean_sd_sightings_2, sdlog = middle_sd_sd_sightings_2)),
    Year %in% c(13:18) ~ (rlnorm(n = n(), meanlog = late_mean_sd_sightings_2, sdlog = late_sd_sd_sightings_2 ))
    ))
# Join the Blocks_2
Blocks_2 <- left_join(Blocks_2, Grid_Years_Species_select)





# --------------------------------------------------
# Adjust JdDetect
# --------------------------------------------------
Daily_Site_Total_zero_sum_df <- Blocks_2 %>%
  group_by(ID) %>%
  filter(Daily_Site_Total == 0) %>%
  summarize(Daily_Site_Total_zero_sum = n()) %>%
  ungroup()
Daily_Site_Total_zero_sum_df$Daily_Site_Total_zero_sum <- str_replace_na(Daily_Site_Total_zero_sum_df$Daily_Site_Total_zero_sum, 0)
# Join
Blocks_2 <- left_join(Blocks_2, Daily_Site_Total_zero_sum_df)
# Where NA, put 0
Blocks_2$Daily_Site_Total_zero_sum <- str_replace_na(Blocks_2$Daily_Site_Total_zero_sum, 0)
Blocks_2$Daily_Site_Total_zero_sum <- as.numeric(Blocks_2$Daily_Site_Total_zero_sum)
# Daily_Site_Total_zero_sum_2
Blocks_2 <- Blocks_2 %>%
  mutate(Daily_Site_Total_zero_sum_2 = NJd - Daily_Site_Total_zero_sum)
# See if there are differences
Blocks_2 <- Blocks_2 %>%
  mutate(JdDetect_Stationary_flag = ifelse(JdDetect_Stationary > Daily_Site_Total_zero_sum_2 , TRUE, FALSE))
Blocks_2 <- Blocks_2 %>%
  mutate(JdDetect_Shifting_flag = ifelse(JdDetect_Shifting > Daily_Site_Total_zero_sum_2 , TRUE, FALSE))
# Also, times where there is no Daily_Site_Total == 0
# Count
Blocks_2 %>%
  count(JdDetect_Stationary_flag)
Blocks_2 %>%
  count(JdDetect_Shifting_flag)

# catch <- Blocks_2 %>%
#   filter(JdDetect_flag == TRUE)
#
# catch %>%
#   select(ID, JdDetect, Daily_Site_Total_zero_sum_2) %>%
#   distinct() %>%
# ggplot(aes(JdDetect, Daily_Site_Total_zero_sum_2)) +
#   geom_point(alpha = 0.2) +
#   geom_abline()
#
# rm(catch)

# Set JdDetect to be this number
Blocks_2 <- Blocks_2 %>%
  mutate(JdDetect_Stationary = if_else(JdDetect_Stationary_flag == TRUE, Daily_Site_Total_zero_sum_2, JdDetect_Stationary))
Blocks_2 %>%
  filter(JdDetect_Stationary > Daily_Site_Total_zero_sum_2)

Blocks_2 <- Blocks_2 %>%
  mutate(JdDetect_Shifting = if_else(JdDetect_Shifting_flag == TRUE, Daily_Site_Total_zero_sum_2, JdDetect_Shifting))
Blocks_2 %>%
  filter(JdDetect_Shifting > Daily_Site_Total_zero_sum_2)



# # Add column
# Blocks_2$Daily_Site_Sightings_Stationary <- NA
# Blocks_2$Daily_Site_Sightings_Shifting <- NA
# Blocks_2$Daily_Site_Sightings_Stationary <- as.numeric(Blocks_2$Daily_Site_Sightings_Stationary)
# Blocks_2$Daily_Site_Sightings_Shifting <- as.numeric(Blocks_2$Daily_Site_Sightings_Shifting)
Blocks_2 <- as.data.table(Blocks_2)
write_csv(Blocks_2, "Blocks_2_final.csv")
# Blocks_2 <- fread("Blocks_2.csv")
# # Get unique Blocks_2 IDs
Blocks_2_IDs <- unique(Blocks_2$ID)


# ---------------------------------------
# Test out ecodist
# ---------------------------------------
ecodist_test <- Blocks_2 %>%
  filter(ID == Blocks_2_IDs[1])

xy <- corgen(len = unique(ecodist_test$NJd),
             x = ecodist_test$Daily_Site_Total,
             r = unique(ecodist_test$daily_cor),
             epsilon = 0)
xy

xy2 <- xy
center <- attr(xy2$y, "scaled:center")
scale <- attr(xy2$y, "scaled:scale")
xy2$y <- xy2$y * scale
xy2$y <- xy2$y + center
xy2

# Plot
plot(xy$x, xy$y)
# Scaled
plot(ecodist_test$Jd, xy$x, col = "red")
points(ecodist_test$Jd, xy$y[,1])
# Unscaled
plot(ecodist_test$Jd, xy$x, col = "red")
points(ecodist_test$Jd, xy2$y[,1])
# Occupancy
plot(ecodist_test$Jd, xy2$y[,1]/xy2$x)


# ---------------------------------------
# Test out simcor
# ---------------------------------------
simcor_test <- Blocks_2 %>%
  filter(ID == Blocks_2_IDs[1])

simcor_test_unc <- simcor_test %>%
  uncount(Daily_Site_Total)

ggplot(simcor_test, aes(Jd, Daily_Site_Total)) +
  geom_point()

# ---------------------------------------
# Test out mvrnorm
# ---------------------------------------
mvrnorm_test <- Blocks_2 %>%
  filter(ID == Blocks_2_IDs[1])

mvrnorm_test_unc <- mvrnorm_test %>%
  uncount(Daily_Site_Total)


MASS::mvrnorm(n = unique(mvrnorm_test$NJd),
              mu = c(mean(mvrnorm_test_unc$Jd), unique(mvrnorm_test$MAD_Stationary)),
              Sigma = matrix(c(1, unique(mvrnorm_test_unc$daily_cor), unique(mvrnorm_test_unc$daily_cor), 1), ncol = 2),
              empirical = TRUE)


# ---------------------------------------
# Test out fabricatr
# ---------------------------------------

# test data.frame
test <- Blocks_2
set.seed(1)
# sample_blocks <- sample(Blocks_2_IDs, size = 1000)
sample_species <- sample(unique(test$Common_Name), 3)
test <- test %>%
  filter(Common_Name %in% sample_species)
sample_blocks <- unique(test$ID)


# Create list for storing results
big_list <- vector(mode = "list", length = length(Blocks_2_IDs))

# Number of iterations
imax <- c(length(unique(Blocks_2$ID)))
imax_Blocks_2 <- c(length(unique(Blocks_2$ID)))


start_time <- Sys.time()
# Initiate the bar
pb <- txtProgressBar(min = 0, max = imax, style = 3)


# Loop Stationary
for (i in 42353:length(unique(Blocks_2$ID))) {
# for (i in 25821:length(unique(Blocks_2$ID))) {
# for (i in 1:29) {
# for (i in 1:length(unique(Blocks_2$ID))) {
    # i = 1
    i0 = i - 1
    ID_i <- Blocks_2_IDs[i]
    df_i <- Blocks_2[ID == ID_i]
    # df_i <- filter(Blocks_2, ID == ID_i)
    set.seed(i)
    vec_i <- ifelse(df_i$Departed_spp == FALSE,
                    pnorm(q = seq(unique(df_i$MinJd), unique(df_i$MaxJd), 1),
                          mean = unique(df_i$MAD_Stationary),
                          sd = ifelse(
                            unique(df_i$Year) %in% c(1:6), early_mean_sd_sightings,
                            ifelse(unique(df_i$Year) %in% c(7:12), middle_mean_sd_sightings,
                                   late_mean_sd_sightings)
                          )),
                    dnorm(x = seq(unique(df_i$MinJd), unique(df_i$MaxJd), 1),
                          mean = unique(df_i$MAD_Stationary),
                          sd = ifelse(
                            unique(df_i$Year) %in% c(1:6), early_mean_sd_sightings,
                            ifelse(unique(df_i$Year) %in% c(7:12), middle_mean_sd_sightings,
                                   late_mean_sd_sightings)
                          )))
    # Multiply by Daily_Site_Total values converted to 0 and 1 (anything positive)
    site_total_1_0_i <- ifelse(df_i$Daily_Site_Total > 0, 1, 0)
    vec_i <- vec_i * site_total_1_0_i
    # Only x days can be fpositive, where x is JdDetect
    JdDetect_i <- unique(df_i$JdDetect_Stationary)
    # JdDetect_inv_i <- length(vec_i) - JdDetect_i
    vec_detect_i <- data.frame(vec_i = vec_i)
    # Get difference between days detected and days that are supposed to be detected
    vec_detect_i_diff <- length(vec_i[vec_i > 0 ]) - JdDetect_i
      # Add needed columns
      vec_detect_i <- vec_detect_i %>%
        mutate(detect_flag = if_else(vec_i > 0 , TRUE, FALSE),
               row_id = seq(unique(df_i$MinJd), unique(df_i$MaxJd), 1),
               Dothis = if_else(vec_detect_i_diff >= 0, TRUE, FALSE))
      # Join Daily_Site_Total
      vec_detect_i <- left_join(vec_detect_i, select(df_i, "row_id" = "Jd", Daily_Site_Total), by = "row_id")
      # Break data.frame into TRUE and FALSE
      vec_detect_i_TRUE <- vec_detect_i %>%
        filter(detect_flag == TRUE)
      vec_detect_i_FALSE <- vec_detect_i %>%
        filter(detect_flag == FALSE)
    if(vec_detect_i_diff >= 0){
      # Create a subtraction/addition vector to enforce a set amount of detections
      pt1 <- floor(vec_detect_i_diff/5)
      pt2 <- floor(vec_detect_i_diff/2.5)
      pt3 <- ceiling(vec_detect_i_diff/2.5)
      pttally1 <- pt1 + pt2 + pt3
      pttally1diff <- vec_detect_i_diff - pttally1
      pt3 <- pt3 + pttally1diff
      pt4 <- floor(JdDetect_i/5)
      pt5 <- ceiling(JdDetect_i/1.25)
      pttally2 <- pt4 + pt5
      pttally2diff <- JdDetect_i - pttally2
      pt5 <- pt5 + pttally2diff

      set.seed(i)
      subtract_vec_i <- c(
        sample(c(rep(0, (abs(pt1))))),
        sample(c(rep(0, (abs(pt2))), rep(1, (pt4)))),
        sample(c(rep(0, (abs(pt3))), rep(1, (pt5))))
      )

      # Bind to data.frame
      vec_detect_i_TRUE <- cbind(vec_detect_i_TRUE, subtract_vec_i)
      # Multiply the vectors together if vec_detect_i_diff > 0
      vec_detect_i_TRUE <- vec_detect_i_TRUE %>%
        mutate(vec_i = if_else(
          Dothis == TRUE,
          vec_i * subtract_vec_i,
          vec_i
        ))
      vec_detect_i_TRUE$subtract_vec_i <- NULL
      # Add 1 to the FALSE vector if vec_detect_i_diff == 0
      vec_detect_i_FALSE_length <- length(vec_detect_i_FALSE$vec_i)
      set.seed(i)
      if(vec_detect_i_FALSE_length > 0){
        vec_detect_i_FALSE_length_i <- sample(c(1:vec_detect_i_FALSE_length), size = 1)
        vec_detect_i_FALSE$vec_i[vec_detect_i_FALSE_length_i] <- ifelse(
          vec_detect_i_TRUE$Dothis[vec_detect_i_FALSE_length_i] == FALSE,
          vec_detect_i_FALSE$vec_i[vec_detect_i_FALSE_length_i] + 1,
          vec_detect_i_FALSE$vec_i[vec_detect_i_FALSE_length_i]
        )

        vec_detect_i <- rbind(vec_detect_i_TRUE, vec_detect_i_FALSE)
      } else {
        vec_detect_i <- vec_detect_i_TRUE
      }
      # bind the TRUE and FALSE data.frames together and reorder

      vec_detect_i <- vec_detect_i %>%
        arrange(-desc(row_id))
    } else {
      # print("skipped")
      vec_detect_i <- vec_detect_i
    }
    # print(ifelse(vec_detect_i_diff > 0, "vec_detect_i_diff positive", "vec_detect_i_diff negative"))

    # Get vec_i sum
    # Correct vec_i to Sighting_Sum
    sighting_sum_i <- unique(df_i$Sighting_Sum)
    correction_factor <- sighting_sum_i/sum(vec_detect_i$vec_i)
    vec_detect_i$vec_i <- vec_detect_i$vec_i * correction_factor
    # Introduce some random variability to the positive values
    # mean_df <- mean(vec_detect_i[vec_detect_i$vec_i > 0,]$vec_i)
    sd_df <- sd(vec_detect_i[vec_detect_i$vec_i > 0,]$vec_i)
    set.seed(i)
    vec_detect_i <- vec_detect_i %>%
      mutate(vec_i = if_else(
        vec_i > 0,
        vec_i + rnorm(n = length(c(unique(df_i$MinJd):unique(df_i$MaxJd))), mean = 0, sd = unique(df_i$sd_2_sightings)/4),
        vec_i
      ))
    # Correct negative values to above zero
    vec_detect_i <- vec_detect_i %>%
      mutate(vec_i = if_else(vec_i < 0, vec_i*-1, vec_i))
    # Get vec_i sum agains
    # How many sites is this?
    correction_factor_2 <- sighting_sum_i / sum(vec_detect_i$vec_i)
    # Multiply the vector by the correction factor again
    vec_detect_i$vec_i <- vec_detect_i$vec_i * correction_factor_2
    # Round the vector
    vec_detect_i$vec_i <- ceiling(vec_detect_i$vec_i)

    # Convert any values that are greater than or equal to Daily_Site_Total to be less than Daily_Site_Total
    set.seed(i)
    mult <- round(rnorm(n = 1, mean = 0.5, sd = 0.05), digits = 1)

    # Shift downwards
    vec_detect_i <- vec_detect_i %>%
      mutate(vec_i= if_else(
        vec_i >= Daily_Site_Total,
        floor(vec_i * mult),
        vec_i
      ))

    # If this resulted in JdDetect being less than the number of positive days, then fill in some extra days at which == 0
    diff_interest <- JdDetect_i - nrow(vec_detect_i[vec_detect_i$vec_i > 0,])
    position_interest <- which(vec_detect_i$vec_i == 0 & vec_detect_i$detect_flag == TRUE)
    position_interest <- position_interest[length(position_interest):(length(position_interest) - diff_interest + 1)]
    vec_detect_i$vec_i[position_interest] <- ifelse(
      nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) < JdDetect_i,
      1,
      vec_detect_i$vec_i[position_interest]
    )

    # What is the absolute difference between the sighting numbers?
    number <- sum(vec_detect_i$vec_i)
    vec_detect_i <- vec_detect_i %>%
      mutate(vec_diff = sighting_sum_i - number)
    # What is the difference between vec_i and Daily_Site_Total for each day?
    vec_detect_i <- vec_detect_i %>%
      mutate(vec_i_diff = Daily_Site_Total - vec_i)
    # How many observations are positive and have a difference greater than max vec_i
    vec_detect_i_1 <- vec_detect_i %>%
      filter(vec_i > 1 & vec_i_diff > ifelse(max(vec_detect_i$vec_i) < 5,  max(vec_detect_i$vec_i), 5))
    # # Filter for days that have vec_i > 1 and differences that are low
    # vec_detect_i_2 <- vec_detect_i %>%
    #   filter(vec_i > 1 & vec_i_diff <= quantile(vec_detect_i_1$vec_i_diff)[[4]])
    # Filter for the rest
    vec_detect_i_3 <- vec_detect_i %>%
      filter(!row_id %in% c(vec_detect_i_1$row_id
                            # , vec_detect_i_2$row_id
                            ))
    # If vec_diff is negative, we need to remove observations
    # If vec_diff is positive, we need to add observations
    rows_1 <- nrow(vec_detect_i_1)
    rows_1 <- as.numeric(rows_1)
    # rows_2 <- nrow(vec_detect_i_2)
    set.seed(i)
    if (rows_1 > 2) {
    vec_adjust_1 <- rand_vect(N = rows_1,
                M = abs(unique(vec_detect_i$vec_diff)),
                sd = sd(vec_detect_i$vec_i),
                max = max(vec_detect_i$vec_i)-1,
                min = min(vec_detect_i$vec_i),
                pos.only = TRUE)
    # Sort in increasing or decreasing, order, depending on the sign of vec_diff (increasing if negative)
    vec_adjust_1 <- sort(vec_adjust_1, decreasing = if_else(unique(vec_detect_i_1$vec_diff) > 0, TRUE, FALSE))
    } else if (rows_1 == 1){
      vec_adjust_1 <- abs(unique(vec_detect_i$vec_diff))
    } else {
    vec_adjust_1 <- 0
    }
    #
    # rnormTrunc(n = rows_1, mean = abs(unique(vec_detect_i$vec_diff)), sd = sd(vec_detect_i$vec_i), max = max(vec_detect_i$vec_i)-1, min = min(vec_detect_i$vec_i))

    # Sort data.frame so that the difference between vec_i_diff is ascending or descending depending on vec_diff
    vec_detect_i_1 <- vec_detect_i_1 %>%
      arrange(if_else(vec_diff > 0,
                      desc(vec_i),
                      -desc(vec_i)))
    # If vec_adjust_1 has any values that are equal to the maximum value of vec_i, then subtract one from these, and then place 1s where there are zeros in vec_adjust_1
    at_max_vec_i <- which(vec_adjust_1 == median(vec_detect_i_1$vec_i))
    vec_adjust_1[at_max_vec_i] <- vec_adjust_1[at_max_vec_i] - 1
    zero_places <- which(vec_adjust_1 == 0)
    zero_places <- zero_places[1:length(at_max_vec_i)]
    vec_adjust_1[zero_places] <- ifelse(length(at_max_vec_i) > 0,  vec_adjust_1[zero_places] + 1, vec_adjust_1[zero_places])

    # Add vec_adjust to adding data.frame if vec_diff is positive
    if(length(at_max_vec_i) > 0){
      vec_detect_i_1 <- vec_detect_i_1 %>%
        mutate(vec_i = if_else(
          vec_diff > 0,
          vec_i + c(vec_adjust_1),
          vec_i - c(vec_adjust_1)
        ))
    } else {
      vec_detect_i_1 <- vec_detect_i_1
    }
    # Nothing needs to happen if vec_diff is zero
    # Bind together and sort
    vec_detect_i <- rbind(vec_detect_i_1,
                          vec_detect_i_3)
    vec_detect_i <- vec_detect_i %>%
      arrange(-desc(row_id))

    # Which rows are exceeding daily_site_total?
    # Recalculate:
    # What is the difference between vec_i and Daily_Site_Total for each day?
    vec_detect_i <- vec_detect_i %>%
      mutate(vec_i_diff = Daily_Site_Total - vec_i)
    # problem_rows <- which(vec_detect_i$vec_i> 0 & vec_detect_i$vec_i > vec_detect_i$Daily_Site_Total)
    vec_detect_i_problem <- vec_detect_i %>%
      filter(vec_i_diff < 0)
    vec_detect_i_problem <- vec_detect_i_problem %>%
      mutate(vec_i = vec_i - abs(vec_i_diff))
    vec_detect_i_noprob <- vec_detect_i %>%
      filter(!row_id %in% vec_detect_i_problem$row_id)
    # Bind together and sort
    vec_detect_i <- rbind(vec_detect_i_problem,
                          vec_detect_i_noprob)
    vec_detect_i <- vec_detect_i %>%
      arrange(-desc(row_id))
    # Add to the rows that have the greaBlocks_2 difference between Daily_Site_Total and vec_i
    # and are in the latter half of the time series
    big_gap_rows <- vec_detect_i %>%
      filter(vec_i > 0,
             row_id > quantile(row_id)[[2]]) %>%
      slice_max(vec_i_diff, n = nrow(vec_detect_i_problem))
    big_gap_rows <- big_gap_rows %>%
      arrange(desc(vec_i_diff))

    if(nrow(big_gap_rows) > 1){
      big_gap_rows <- big_gap_rows[1:nrow(vec_detect_i_problem),]
    } else {
      big_gap_rows <- big_gap_rows
    }

    small_gap_rows <- vec_detect_i %>%
      filter(!row_id %in% big_gap_rows$row_id)
    big_gap_rows <- big_gap_rows %>%
      mutate(Dothis2 = ifelse(nrow(vec_detect_i_problem) > 0, TRUE, FALSE))
    # big_gap_rows$Dothis2 <- ifelse(nrow(vec_detect_i_problem) > 0, TRUE, FALSE)
    big_gap_rows <- big_gap_rows %>%
      mutate(vec_i = if_else(
        Dothis2 == TRUE,
        vec_i + c(abs(vec_detect_i_problem$vec_i_diff)),
        vec_i))
    big_gap_rows$Dothis2 <- NULL
    # Bind together and sort
    vec_detect_i <- rbind(big_gap_rows,
                          small_gap_rows)
    vec_detect_i <- vec_detect_i %>%
      arrange(-desc(row_id))

    # Make sure one last time that Daily_Site_Total is not exceeded by vec_i
    vec_detect_i <- vec_detect_i %>%
      mutate(vec_i = ifelse(vec_i > Daily_Site_Total, Daily_Site_Total, vec_i))
    # Replace NAs with 0
    if(anyNA(vec_detect_i$vec_i) == TRUE) {
      vec_detect_i <- vec_detect_i %>%
        mutate(vec_i = str_replace_na(vec_i, 0))
      vec_detect_i$vec_i <- as.numeric(vec_detect_i$vec_i)
    } else {
      vec_detect_i <- vec_detect_i
    }

    # Replace negatives with 0s
    if(any(vec_detect_i$vec_i < 0) == TRUE) {
      vec_detect_i <- vec_detect_i %>%
        mutate(vec_i = if_else(vec_i < 0, 0, vec_i))
    } else {
      vec_detect_i <- vec_detect_i
    }



    setTxtProgressBar(pb, i)
    # print(i)
    # # Recalculate:
    # vec_detect_i <- vec_detect_i %>%
    #   mutate(vec_i_diff = Daily_Site_Total - vec_i)
    #
    # # If we still have problems, add some observations
    # if(nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) < JdDetect_i){
    #   # How many days do we need to add?
    #   needed_days <- JdDetect_i - nrow(vec_detect_i[vec_detect_i$vec_i > 0,])
    #   # How many observations do we need to add?
    #   needed_obs <- unique(df_i$Sighting_Sum) - sum(vec_detect_i$vec_i)
    #   # If needed_obs is negative, turn it to 1
    #   needed_obs <- ifelse(needed_obs < 0, 1, needed_obs)
    #   # Calculate the ratio between obs/days to find an appropriate subset
    #   ratio <- ceiling(needed_obs/needed_days)
    #   # Filter for a number of eligible days
    #   set.seed(1)
    #   vec_detect_i_subin <- vec_detect_i %>%
    #     filter(detect_flag == TRUE & vec_i_diff >= ratio & vec_i == 0) %>%
    #     sample_n(needed_days)
    #   vec_detect_i_rest <- vec_detect_i %>%
    #     filter(!row_id %in% vec_detect_i_subin$row_id)
    #   # Create vector needed
    #   N1 <- floor(needed_obs/needed_days)
    #   N2 <- needed_obs%%needed_days
    #   N3 <- c(N1, N2)
    #   N3 <- N3[N3 > 0]
    #   # Add to vec_detect_i_subin
    #   vec_detect_i_subin <- vec_detect_i_subin %>%
    #     mutate(vec_i = vec_i + N3)
    #   # Bind together and sort
    #   vec_detect_i <- rbind(vec_detect_i_subin,
    #                         vec_detect_i_rest)
    #   vec_detect_i <- vec_detect_i %>%
    #     arrange(-desc(row_id))
    # } else if(nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) > JdDetect_i) {
    #   # How many days do we need to subtract?
    #   needed_days <- nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) - JdDetect_i
    #
    #   vec_detect_i <- rbind(vec_detect_i_subin,
    #                         vec_detect_i_rest)
    #   vec_detect_i <- vec_detect_i %>%
    #     arrange(-desc(row_id))
    # } else if (nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) == JdDetect_i){
    #   vec_detect_i <- vec_detect_i
    # }
    #
    # # If you have too many observations?
    # if(sum(vec_detect_i$vec_i) > unique(df_i$Sighting_Sum))

    # If you have the wrong amount of days detected?
    if(abs(nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) - JdDetect_i) != 0)
    {
      df_i <- df_i %>%
        mutate(JdDetect_Stationary_2 = nrow(vec_detect_i[vec_detect_i$vec_i > 0,]))
    } else {
      df_i <- df_i %>%
        mutate(JdDetect_Stationary_2 = JdDetect_Stationary)
    }

    # If you have the wrong amount of observations detected?
    if((abs(sum(vec_detect_i$vec_i) - unique(df_i$Sighting_Sum))) != 0)
    {
      df_i <- df_i %>%
        mutate(Sighting_Sum_2 = sum(vec_detect_i$vec_i))
    } else {
      df_i$Sighting_Sum_2 <- df_i$Sighting_Sum
    }

    if((nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) == unique(df_i$JdDetect_Stationary_2)) == FALSE)
    {
      break
    }
    # print("JdDetect Good")

    if((sum(vec_detect_i$vec_i) == unique(df_i$Sighting_Sum_2)) == FALSE)
    {
      break
    }
    # print("Sighting_Sum Good")


    head(vec_detect_i)
    vec_detect_i_join <- vec_detect_i %>%
      select("Daily_Site_Sightings_Stationary" = "vec_i", "Jd" = "row_id")
    df_i_2 <- df_i %>%
      select(Jd, ID, JdDetect_Stationary_2, Sighting_Sum_2)
    vec_detect_i_join <- left_join(df_i_2, vec_detect_i_join, by = "Jd")
    # Join
    big_list[i] <- list(vec_detect_i_join)

    # # Add to Daily_Site_Total column (26)
    # start_ind <- ifelse(i == 1, 100*i0 + 1, end_ind + 1)
    # end_ind <- ifelse(i == 1, 100*i0 + 1 + unique(df_i$MaxJd - df_i$MinJd), start_ind + unique(df_i$NJd))
    # # Blocks_2 <- Blocks_2
    # Blocks_2[c(start_ind:end_ind), 26] <- vec_detect_i$vec_i
  }

end_time <- Sys.time()
end_time - start_time

big_list_rbind <- rbindlist(big_list)
#
#
#
# # Blocks_2 <- left_join(Blocks_2, big_list_rbind)
# #
# # length(Blocks_2_IDs)/29
# #
# # Blocks_22 <- Blocks_2 %>%
# #   filter(!is.na(Daily_Site_Sightings_Stationary))
# #
# # ggplot(Blocks_22) +
# #   geom_line(aes(Jd, Daily_Site_Total), colour = "blue") +
# #   geom_line(aes(Jd, Daily_Site_Sightings_Stationary), colour = "red", linetype = "dashed") +
# #   facet_wrap(~ID)
# #
# # Blocks_22 <- Blocks_22 %>%
# #   mutate(Daily_Site_Occupancy_Stationary = Daily_Site_Sightings_Stationary / Daily_Site_Total)
# # Blocks_22$Daily_Site_Occupancy_Stationary[is.nan(Blocks_22$Daily_Site_Occupancy_Stationary)] <- 0
# #
# # ggplot(Blocks_22) +
# #   geom_line(aes(Jd, Daily_Site_Occupancy_Stationary), colour = "black") +
# #   facet_wrap(~ID)
#
# -----------------------------
# Shifting
# -----------------------------
# Create list for storing results
big_list_2 <- vector(mode = "list", length = length(Blocks_2_IDs))

# Number of iterations
imax <- c(length(unique(Blocks_2$ID)))


start_time <- Sys.time()
# Initiate the bar
pb <- txtProgressBar(min = 0, max = imax_Blocks_2, style = 3)

# Loop Shifting
for (i in 1:length(unique(Blocks_2$ID))) {
# for (i in 1:29) {
# for (i in 25821:length(unique(Blocks_2$ID))) {
  # i = 44
  i0 = i - 1
  ID_i <- Blocks_2_IDs[i]
  df_i <- Blocks_2[ID == ID_i]
  # df_i <- filter(Blocks_2, ID == ID_i)
  set.seed(i)
  vec_i <- ifelse(df_i$Departed_spp == FALSE,
                  pnorm(q = seq(unique(df_i$MinJd), unique(df_i$MaxJd), 1),
                        mean = unique(df_i$MAD_Shifting),
                        sd = ifelse(
                          unique(df_i$Year) %in% c(1:6), early_mean_sd_sightings,
                          ifelse(unique(df_i$Year) %in% c(7:12), middle_mean_sd_sightings,
                                 late_mean_sd_sightings)
                        ))*100,
                  dnorm(x = seq(unique(df_i$MinJd), unique(df_i$MaxJd), 1),
                        mean = unique(df_i$MAD_Shifting),
                        sd = ifelse(
                          unique(df_i$Year) %in% c(1:6), early_mean_sd_sightings,
                          ifelse(unique(df_i$Year) %in% c(7:12), middle_mean_sd_sightings,
                                 late_mean_sd_sightings)
                        )))
  # Multiply by Daily_Site_Total values converted to 0 and 1 (anything positive)
  site_total_1_0_i <- ifelse(df_i$Daily_Site_Total > 0, 1, 0)
  vec_i <- vec_i * site_total_1_0_i
  # Only x days can be positive, where x is JdDetect
  JdDetect_i <- unique(df_i$JdDetect_Shifting)
  # JdDetect_inv_i <- length(vec_i) - JdDetect_i
  vec_detect_i <- data.frame(vec_i = vec_i)
  # Get difference between days detected and days that are supposed to be detected
  vec_detect_i_diff <- length(vec_i[vec_i > 0 ]) - JdDetect_i
  # Add needed columns
  vec_detect_i <- vec_detect_i %>%
    mutate(detect_flag = if_else(vec_i > 0 , TRUE, FALSE),
           row_id = seq(unique(df_i$MinJd), unique(df_i$MaxJd), 1),
           Dothis = if_else(vec_detect_i_diff >= 0, TRUE, FALSE))
  # Join Daily_Site_Total
  vec_detect_i <- left_join(vec_detect_i, select(df_i, "row_id" = "Jd", Daily_Site_Total), by = "row_id")
  # Break data.frame into TRUE and FALSE
  vec_detect_i_TRUE <- vec_detect_i %>%
    filter(detect_flag == TRUE)
  vec_detect_i_FALSE <- vec_detect_i %>%
    filter(detect_flag == FALSE)
  if(vec_detect_i_diff >= 0){
    # Create a subtraction/addition vector to enforce a set amount of detections
    pt1 <- floor(vec_detect_i_diff/5)
    pt2 <- floor(vec_detect_i_diff/2.5)
    pt3 <- ceiling(vec_detect_i_diff/2.5)
    pttally1 <- pt1 + pt2 + pt3
    pttally1diff <- vec_detect_i_diff - pttally1
    pt3 <- pt3 + pttally1diff
    pt4 <- floor(JdDetect_i/5)
    pt5 <- ceiling(JdDetect_i/1.25)
    pttally2 <- pt4 + pt5
    pttally2diff <- JdDetect_i - pttally2
    pt5 <- pt5 + pttally2diff

    set.seed(i)
    subtract_vec_i <- c(
      sample(c(rep(0, (abs(pt1))))),
      sample(c(rep(0, (abs(pt2))), rep(1, (pt4)))),
      sample(c(rep(0, (abs(pt3))), rep(1, (pt5))))
    )

    # Bind to data.frame
    vec_detect_i_TRUE <- cbind(vec_detect_i_TRUE, subtract_vec_i)
    # Multiply the vectors together if vec_detect_i_diff > 0
    vec_detect_i_TRUE <- vec_detect_i_TRUE %>%
      mutate(vec_i = if_else(
        Dothis == TRUE,
        vec_i * subtract_vec_i,
        vec_i
      ))
    vec_detect_i_TRUE$subtract_vec_i <- NULL
    # Add 1 to the FALSE vector if vec_detect_i_diff == 0
    vec_detect_i_FALSE_length <- length(vec_detect_i_FALSE$vec_i)
    set.seed(i)
    if(vec_detect_i_FALSE_length > 0){
      vec_detect_i_FALSE_length_i <- sample(c(1:vec_detect_i_FALSE_length), size = 1)
      vec_detect_i_FALSE$vec_i[vec_detect_i_FALSE_length_i] <- ifelse(
        vec_detect_i_TRUE$Dothis[vec_detect_i_FALSE_length_i] == FALSE,
        vec_detect_i_FALSE$vec_i[vec_detect_i_FALSE_length_i] + 1,
        vec_detect_i_FALSE$vec_i[vec_detect_i_FALSE_length_i]
      )

      vec_detect_i <- rbind(vec_detect_i_TRUE, vec_detect_i_FALSE)
    } else {
      vec_detect_i <- vec_detect_i_TRUE
    }
    # bind the TRUE and FALSE data.frames together and reorder

    vec_detect_i <- vec_detect_i %>%
      arrange(-desc(row_id))
  } else {
    # print("skipped")
    vec_detect_i <- vec_detect_i
  }
  # print(ifelse(vec_detect_i_diff > 0, "vec_detect_i_diff positive", "vec_detect_i_diff negative"))

  # Get vec_i sum
  # Correct vec_i to Sighting_Sum
  sighting_sum_i <- unique(df_i$Sighting_Sum)
  correction_factor <- sighting_sum_i/sum(vec_detect_i$vec_i)
  vec_detect_i$vec_i <- vec_detect_i$vec_i * correction_factor
  # Introduce some random variability to the positive values
  # mean_df <- mean(vec_detect_i[vec_detect_i$vec_i > 0,]$vec_i)
  sd_df <- sd(vec_detect_i[vec_detect_i$vec_i > 0,]$vec_i)
  set.seed(i)
  vec_detect_i <- vec_detect_i %>%
    mutate(vec_i = if_else(
      vec_i > 0,
      vec_i + rnorm(n = length(c(unique(df_i$MinJd):unique(df_i$MaxJd))), mean = 0, sd = unique(df_i$sd_2_sightings)/2),
      vec_i
    ))
  # Correct negative values to above zero
  vec_detect_i <- vec_detect_i %>%
    mutate(vec_i = if_else(vec_i < 0, vec_i*-1, vec_i))
  # Get vec_i sum agains
  # How many sites is this?
  correction_factor_2 <- sighting_sum_i / sum(vec_detect_i$vec_i)
  # Multiply the vector by the correction factor again
  vec_detect_i$vec_i <- vec_detect_i$vec_i * correction_factor_2
  # Round the vector
  vec_detect_i$vec_i <- ceiling(vec_detect_i$vec_i)

  # Convert any values that are greater than or equal to Daily_Site_Total to be less than Daily_Site_Total
  set.seed(i)
  mult <- round(rnorm(n = 1, mean = 0.5, sd = 0.05), digits = 1)

  # Shift downwards
  vec_detect_i <- vec_detect_i %>%
    mutate(vec_i= if_else(
      vec_i >= Daily_Site_Total,
      floor(vec_i * mult),
      vec_i
    ))

  # If this resulted in JdDetect being less than the number of positive days, then fill in some extra days at which == 0
  diff_interest <- JdDetect_i - nrow(vec_detect_i[vec_detect_i$vec_i > 0,])
  position_interest <- which(vec_detect_i$vec_i == 0 & vec_detect_i$detect_flag == TRUE)
  position_interest <- position_interest[length(position_interest):(length(position_interest) - diff_interest + 1)]
  vec_detect_i$vec_i[position_interest] <- ifelse(
    nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) < JdDetect_i,
    1,
    vec_detect_i$vec_i[position_interest]
  )

  # What is the absolute difference between the sighting numbers?
  number <- sum(vec_detect_i$vec_i)
  vec_detect_i <- vec_detect_i %>%
    mutate(vec_diff = sighting_sum_i - number)
  # What is the difference between vec_i and Daily_Site_Total for each day?
  vec_detect_i <- vec_detect_i %>%
    mutate(vec_i_diff = Daily_Site_Total - vec_i)
  # How many observations are positive and have a difference greater than max vec_i
  vec_detect_i_1 <- vec_detect_i %>%
    filter(vec_i > 1 & vec_i_diff > ifelse(max(vec_detect_i$vec_i) < 5,  max(vec_detect_i$vec_i), 5))
  # # Filter for days that have vec_i > 1 and differences that are low
  # vec_detect_i_2 <- vec_detect_i %>%
  #   filter(vec_i > 1 & vec_i_diff <= quantile(vec_detect_i_1$vec_i_diff)[[4]])
  # Filter for the rest
  vec_detect_i_3 <- vec_detect_i %>%
    filter(!row_id %in% c(vec_detect_i_1$row_id
                          # , vec_detect_i_2$row_id
    ))
  # If vec_diff is negative, we need to remove observations
  # If vec_diff is positive, we need to add observations
  rows_1 <- nrow(vec_detect_i_1)
  rows_1 <- as.numeric(rows_1)
  # rows_2 <- nrow(vec_detect_i_2)
  set.seed(i)
  if (rows_1 > 2) {
    vec_adjust_1 <- rand_vect(N = rows_1,
                              M = abs(unique(vec_detect_i$vec_diff)),
                              sd = sd(vec_detect_i$vec_i),
                              max = max(vec_detect_i$vec_i)-1,
                              min = min(vec_detect_i$vec_i),
                              pos.only = TRUE)
    # Sort in increasing or decreasing, order, depending on the sign of vec_diff (increasing if negative)
    vec_adjust_1 <- sort(vec_adjust_1, decreasing = if_else(unique(vec_detect_i_1$vec_diff) > 0, TRUE, FALSE))
  } else if (rows_1 == 1){
    vec_adjust_1 <- abs(unique(vec_detect_i$vec_diff))
  } else {
    vec_adjust_1 <- 0
  }
  #
  # rnormTrunc(n = rows_1, mean = abs(unique(vec_detect_i$vec_diff)), sd = sd(vec_detect_i$vec_i), max = max(vec_detect_i$vec_i)-1, min = min(vec_detect_i$vec_i))

  # Sort data.frame so that the difference between vec_i_diff is ascending or descending depending on vec_diff
  vec_detect_i_1 <- vec_detect_i_1 %>%
    arrange(if_else(vec_diff > 0,
                    desc(vec_i),
                    -desc(vec_i)))
  # If vec_adjust_1 has any values that are equal to the maximum value of vec_i, then subtract one from these, and then place 1s where there are zeros in vec_adjust_1
  at_max_vec_i <- which(vec_adjust_1 == median(vec_detect_i_1$vec_i))
  vec_adjust_1[at_max_vec_i] <- vec_adjust_1[at_max_vec_i] - 1
  zero_places <- which(vec_adjust_1 == 0)
  zero_places <- zero_places[1:length(at_max_vec_i)]
  vec_adjust_1[zero_places] <- ifelse(length(at_max_vec_i) > 0,  vec_adjust_1[zero_places] + 1, vec_adjust_1[zero_places])

  # Add vec_adjust to adding data.frame if vec_diff is positive
  if(length(at_max_vec_i) > 0){
    vec_detect_i_1 <- vec_detect_i_1 %>%
      mutate(vec_i = if_else(
        vec_diff > 0,
        vec_i + c(vec_adjust_1),
        vec_i - c(vec_adjust_1)
      ))
  } else {
    vec_detect_i_1 <- vec_detect_i_1
  }
  # Nothing needs to happen if vec_diff is zero
  # Bind together and sort
  vec_detect_i <- rbind(vec_detect_i_1,
                        vec_detect_i_3)
  vec_detect_i <- vec_detect_i %>%
    arrange(-desc(row_id))

  # Which rows are exceeding daily_site_total?
  # Recalculate:
  # What is the difference between vec_i and Daily_Site_Total for each day?
  vec_detect_i <- vec_detect_i %>%
    mutate(vec_i_diff = Daily_Site_Total - vec_i)
  # problem_rows <- which(vec_detect_i$vec_i> 0 & vec_detect_i$vec_i > vec_detect_i$Daily_Site_Total)
  vec_detect_i_problem <- vec_detect_i %>%
    filter(vec_i_diff < 0)
  vec_detect_i_problem <- vec_detect_i_problem %>%
    mutate(vec_i = vec_i - abs(vec_i_diff))
  vec_detect_i_noprob <- vec_detect_i %>%
    filter(!row_id %in% vec_detect_i_problem$row_id)
  # Bind together and sort
  vec_detect_i <- rbind(vec_detect_i_problem,
                        vec_detect_i_noprob)
  vec_detect_i <- vec_detect_i %>%
    arrange(-desc(row_id))
  # Add to the rows that have the greaBlocks_2 difference between Daily_Site_Total and vec_i
  # and are in the latter half of the time series
  big_gap_rows <- vec_detect_i %>%
    filter(vec_i > 0,
           row_id > quantile(row_id)[[2]]) %>%
    slice_max(vec_i_diff, n = nrow(vec_detect_i_problem))
  big_gap_rows <- big_gap_rows %>%
    arrange(desc(vec_i_diff))

  if(nrow(big_gap_rows) > 1){
    big_gap_rows <- big_gap_rows[1:nrow(vec_detect_i_problem),]
  } else {
    big_gap_rows <- big_gap_rows
  }

  small_gap_rows <- vec_detect_i %>%
    filter(!row_id %in% big_gap_rows$row_id)
  big_gap_rows <- big_gap_rows %>%
    mutate(Dothis2 = ifelse(nrow(vec_detect_i_problem) > 0, TRUE, FALSE))
  # big_gap_rows$Dothis2 <- ifelse(nrow(vec_detect_i_problem) > 0, TRUE, FALSE)
  big_gap_rows <- big_gap_rows %>%
    mutate(vec_i = if_else(
      Dothis2 == TRUE,
      vec_i + c(abs(vec_detect_i_problem$vec_i_diff)),
      vec_i))
  big_gap_rows$Dothis2 <- NULL
  # Bind together and sort
  vec_detect_i <- rbind(big_gap_rows,
                        small_gap_rows)
  vec_detect_i <- vec_detect_i %>%
    arrange(-desc(row_id))

  # Make sure one last time that Daily_Site_Total is not exceeded by vec_i
  vec_detect_i <- vec_detect_i %>%
    mutate(vec_i = ifelse(vec_i > Daily_Site_Total, Daily_Site_Total, vec_i))
  # Replace NAs with 0
  if(anyNA(vec_detect_i$vec_i) == TRUE) {
    vec_detect_i <- vec_detect_i %>%
      mutate(vec_i = str_replace_na(vec_i, 0))
    vec_detect_i$vec_i <- as.numeric(vec_detect_i$vec_i)
  } else {
    vec_detect_i <- vec_detect_i
  }

  # Replace negatives with 0s
  if(any(vec_detect_i$vec_i < 0) == TRUE) {
    vec_detect_i <- vec_detect_i %>%
      mutate(vec_i = if_else(vec_i < 0, 0, vec_i))
  } else {
    vec_detect_i <- vec_detect_i
  }



  setTxtProgressBar(pb, i)
  # print(i)
  # # Recalculate:
  # vec_detect_i <- vec_detect_i %>%
  #   mutate(vec_i_diff = Daily_Site_Total - vec_i)
  #
  # # If we still have problems, add some observations
  # if(nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) < JdDetect_i){
  #   # How many days do we need to add?
  #   needed_days <- JdDetect_i - nrow(vec_detect_i[vec_detect_i$vec_i > 0,])
  #   # How many observations do we need to add?
  #   needed_obs <- unique(df_i$Sighting_Sum) - sum(vec_detect_i$vec_i)
  #   # If needed_obs is negative, turn it to 1
  #   needed_obs <- ifelse(needed_obs < 0, 1, needed_obs)
  #   # Calculate the ratio between obs/days to find an appropriate subset
  #   ratio <- ceiling(needed_obs/needed_days)
  #   # Filter for a number of eligible days
  #   set.seed(1)
  #   vec_detect_i_subin <- vec_detect_i %>%
  #     filter(detect_flag == TRUE & vec_i_diff >= ratio & vec_i == 0) %>%
  #     sample_n(needed_days)
  #   vec_detect_i_rest <- vec_detect_i %>%
  #     filter(!row_id %in% vec_detect_i_subin$row_id)
  #   # Create vector needed
  #   N1 <- floor(needed_obs/needed_days)
  #   N2 <- needed_obs%%needed_days
  #   N3 <- c(N1, N2)
  #   N3 <- N3[N3 > 0]
  #   # Add to vec_detect_i_subin
  #   vec_detect_i_subin <- vec_detect_i_subin %>%
  #     mutate(vec_i = vec_i + N3)
  #   # Bind together and sort
  #   vec_detect_i <- rbind(vec_detect_i_subin,
  #                         vec_detect_i_rest)
  #   vec_detect_i <- vec_detect_i %>%
  #     arrange(-desc(row_id))
  # } else if(nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) > JdDetect_i) {
  #   # How many days do we need to subtract?
  #   needed_days <- nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) - JdDetect_i
  #
  #   vec_detect_i <- rbind(vec_detect_i_subin,
  #                         vec_detect_i_rest)
  #   vec_detect_i <- vec_detect_i %>%
  #     arrange(-desc(row_id))
  # } else if (nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) == JdDetect_i){
  #   vec_detect_i <- vec_detect_i
  # }
  #
  # # If you have too many observations?
  # if(sum(vec_detect_i$vec_i) > unique(df_i$Sighting_Sum))

  # If you have the wrong amount of days detected?
  if(abs(nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) - JdDetect_i) != 0)
  {
    df_i <- df_i %>%
      mutate(JdDetect_Shifting_2 = nrow(vec_detect_i[vec_detect_i$vec_i > 0,]))
  } else {
    df_i <- df_i %>%
      mutate(JdDetect_Shifting_2 = JdDetect_Shifting)
  }

  # If you have the wrong amount of observations detected?
  if((abs(sum(vec_detect_i$vec_i) - unique(df_i$Sighting_Sum))) != 0)
  {
    df_i <- df_i %>%
      mutate(Sighting_Sum_2 = sum(vec_detect_i$vec_i))
  } else {
    df_i$Sighting_Sum_2 <- df_i$Sighting_Sum
  }

  if((nrow(vec_detect_i[vec_detect_i$vec_i > 0,]) == unique(df_i$JdDetect_Shifting_2)) == FALSE)
  {
    break
  }
  # print("JdDetect Good")

  if((sum(vec_detect_i$vec_i) == unique(df_i$Sighting_Sum_2)) == FALSE)
  {
    break
  }
  # print("Sighting_Sum Good")


  head(vec_detect_i)
  vec_detect_i_join <- vec_detect_i %>%
    select("Daily_Site_Sightings_Shifting" = "vec_i", "Jd" = "row_id")
  df_i_2 <- df_i %>%
    select(Jd, ID, JdDetect_Shifting_2, Sighting_Sum_2)
  vec_detect_i_join <- left_join(df_i_2, vec_detect_i_join, by = "Jd")
  # Join
  big_list_2[i] <- list(vec_detect_i_join)

  # # Add to Daily_Site_Total column (26)
  # start_ind <- ifelse(i == 1, 100*i0 + 1, end_ind + 1)
  # end_ind <- ifelse(i == 1, 100*i0 + 1 + unique(df_i$MaxJd - df_i$MinJd), start_ind + unique(df_i$NJd))
  # # Blocks_2 <- Blocks_2
  # Blocks_2[c(start_ind:end_ind), 26] <- vec_detect_i$vec_i
}


end_time <- Sys.time()
end_time - start_time
# Time difference of 49.25463 mins


big_list_2_rbind <- rbindlist(big_list_2)

# Write big lists
write_csv(big_list_rbind, "big_list_rbind_final.csv")
write_csv(big_list_2_rbind, "big_list_2_rbind_final.csv")
write_csv(Blocks_2, "Blocks_2_final.csv")



# ---------------------------------------------
# Join and clean datasets
# ---------------------------------------------
# # Read big lists
# big_list_rbind <- read_csv("big_list_rbind.csv")
# big_list_2_rbind <- read_csv("big_list_2_rbind.csv")
# Blocks_2 <- read_csv("Blocks_2.csv")
# 
# # Change to JdDetect_2b and Sighting_Sum_2b
# big_list_2_rbind <- big_list_2_rbind %>%
#   rename("JdDetect_Shifting_2b" = "JdDetect_Shifting_2", "Sighting_Sum_2b" = "Sighting_Sum_2")
# 
# # Join big lists to Blocks_2
# colnames(Blocks_2)
# # Join
# # Blocks_2ING
# # Blocks_2 <- Blocks_2
# 
# Blocks_2 <- left_join(Blocks_2, big_list_rbind)
# Blocks_2 <- left_join(Blocks_2, big_list_2_rbind)
# # Write Blocks_2
# # write_csv(Blocks_2, "Blocks_2_final.csv")

# Read Blocks_2
Blocks_2 <- read_csv("Blocks_2_final.csv")

Blocks_2 %>%
  count(is.na(Daily_Site_Sightings_Stationary))
Blocks_2 %>%
  count(is.na(Daily_Site_Sightings_Shifting))


# Replace NAs with 0s
Blocks_2$Daily_Site_Sightings_Shifting <- as.numeric(str_replace_na(Blocks_2$Daily_Site_Sightings_Shifting, 0))
Blocks_2$Daily_Site_Sightings_Stationary <- as.numeric(str_replace_na(Blocks_2$Daily_Site_Sightings_Stationary, 0))

# Replace negative numbers with 0 (fix this in loop if you do this again - force it at the end)
Blocks_2 %>%
  count(Daily_Site_Sightings_Shifting < 0)
Blocks_2 %>%
  count(Daily_Site_Sightings_Stationary < 0)
# # Must have slipped through
# Blocks_2 <- Blocks_2 %>%
#   mutate(Daily_Site_Sightings_Shifting = if_else(Daily_Site_Sightings_Shifting < 0, 0, Daily_Site_Sightings_Shifting))
# Blocks_2 <- Blocks_2 %>%
#   mutate(Daily_Site_Sightings_Stationary = if_else(Daily_Site_Sightings_Stationary < 0, 0, Daily_Site_Sightings_Stationary))
# quantile(Blocks_2$Daily_Site_Sightings_Stationary)
# quantile(Blocks_2$Daily_Site_Sightings_Shifting)

# Calculate proportional occupancy
Blocks_2 <- Blocks_2 %>%
  mutate(Daily_Site_Occupancy_Stationary = Daily_Site_Sightings_Stationary/Daily_Site_Total,
         Daily_Site_Occupancy_Shifting = Daily_Site_Sightings_Shifting/Daily_Site_Total)
# Replace Nan 0/0 with 0
Blocks_2$Daily_Site_Occupancy_Stationary[is.nan(Blocks_2$Daily_Site_Occupancy_Stationary)] <- 0
Blocks_2$Daily_Site_Occupancy_Shifting[is.nan(Blocks_2$Daily_Site_Occupancy_Shifting)] <- 0
quantile(Blocks_2$Daily_Site_Occupancy_Stationary)
quantile(Blocks_2$Daily_Site_Occupancy_Shifting)



# -----------------------------------------------------
# END DATE CODE: Stationary
# -----------------------------------------------------
# -----------------------------------------------------
# Filter Jd ranges for master dataframe
# -----------------------------------------------------
# speciesname <- spp_names[28]
# speciesname <- str_replace_all(speciesname, " ", "_")
# for (j in 1:15) {
# enddateadd <- j
enddateadd <- 5
# for (k in c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) {
# Jd_min <- k
Jd_min <- unique(Blocks_2$MinJd)
# for (k in c(140, 150, 160, 170, 180, 190, 200)) {
# Jd_max <- (k + enddateadd)
Jd_max <- (180 + enddateadd)
# Jd_max For writing the dataframe
Jd_max_write <- (180)
# Jd_max_write <- (k)
# for (l in c(0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25)) {
# Occupancy_Threshold <- l
Occupancy_Threshold <- 0.1
Blocks_2 <- Blocks_2 %>%
  filter(Jd >= Jd_min & Jd <= Jd_max)
# Using Jd = 185 because the last 5 days will be removed due to the later code
Blocks_2_enddate_stationary <- Blocks_2 %>%
  select(Year, Common_Name, GridID, Jd, Daily_Site_Sightings_Stationary, Daily_Site_Total, Daily_Site_Occupancy_Stationary)
# Apply end-date filter to each species to remove dates after which occupancy drop and the model becomes difficult to fit
# Add column that shows the maximum daily occupancy for each bird-year. Try a 10% of Max-Occ Threshold as a starting point
Blocks_2_enddate_stationary <- Blocks_2_enddate_stationary %>%
  group_by(GridID,
    Common_Name,
           Year) %>%
  mutate(Max_Occ = max(Daily_Site_Occupancy_Stationary)) %>%
  mutate(Max_Occ_Fraction = Max_Occ*Occupancy_Threshold) %>%
  ungroup()
# Add column that shows the date of the maximum daily occupancy for each bird-year
Blocks_2_enddate_stationary <- Blocks_2_enddate_stationary %>%
  group_by(GridID,
           Common_Name,
           Year) %>%
  mutate(Max_Occ_Jd = ifelse(Max_Occ == Daily_Site_Occupancy_Stationary, Jd, 1000)) %>%
  ungroup()
Blocks_2_enddate_stationary <- Blocks_2_enddate_stationary %>%
  group_by(GridID,
           Common_Name,
           Year) %>%
  mutate(Max_Occ_Jd_Min = min(Max_Occ_Jd)) %>%
  ungroup()
Blocks_2_enddate_stationary$Max_Occ_Jd <- NULL
# Remove observations that occur after the date of maximum occupancy and have total site occupancy metrics that are lower than the total site occupancy fraction
# Add column that shows this (labelled "Departed", with 1 representing that the bird has departed from the landscape, and 0 representing that it still persists)
Blocks_2_enddate_stationary <- Blocks_2_enddate_stationary %>%
  mutate(Departed = ifelse(Jd > Max_Occ_Jd_Min & Daily_Site_Occupancy_Stationary <= Max_Occ_Fraction, 1, 0))
# Add column that shows the rolling sum of every certain number of days (start with 5) from the departed column. This essentially shows the last day during a certain day rolling sum where the rolling sum meets the interval that has been set to calculate the sum with (i.e. if the rolling sum is 5, and the interval is 5 days, then it means the birds havve been gone for 5 days)
n <- enddateadd
Blocks_2_enddate_stationary$Departed_Roll_Sum <- c(rep_len(NA, n - 1), rowSums(embed(Blocks_2_enddate_stationary$Departed, n)))
# Remove all observations after the first instance of the maximum sum
Blocks_2_enddate_stationary <- Blocks_2_enddate_stationary %>%
  group_by(GridID,
           Common_Name,
           Year) %>%
  mutate(First_Max_Sum = min(which(Departed_Roll_Sum == n | row_number() == n()))) %>%
  filter(row_number() <= First_Max_Sum) %>%
  ungroup()
# Remove additional rows of the maximum sum
# Blocks_2 this**********
Blocks_2_enddate_stationary <- Blocks_2_enddate_stationary %>%
  group_by(GridID,
           Common_Name,
           Year) %>%
  mutate(Spp_Year_Rows = n()) %>%
  filter(Spp_Year_Rows >= 5) %>%
  slice(1:(n()-n)) %>%
  ungroup()

# Blocks_2 <- Blocks_2_enddate_stationary %>%
#   filter(ID == Blocks_2_IDs[1])

# # Plot
# # Daily_Site_Occupancy
# plot(example_enddate[example_enddate$Common_Name == "One",]$Jd,
#      example_enddate[example_enddate$Common_Name == "One",]$Daily_Site_Occupancy,
#      col = "red",
#      xlab = "Day of year",
#      ylab = "Daily site occupancy")
# points(example_enddate[example_enddate$Common_Name == "Two",]$Jd,
#        example_enddate[example_enddate$Common_Name == "Two",]$Daily_Site_Occupancy,
#        col = "blue")
# head(example_enddate)
# tail(example_enddate)
#
# # Plot
# # Daily_Site_Occupancy
# plot(example[example$Common_Name == "One",]$Jd,
#      example[example$Common_Name == "One",]$Daily_Site_Occupancy,
#      col = "red",
#      xlab = "Day of year",
#      ylab = "Daily site occupancy")
# points(example[example$Common_Name == "Two",]$Jd,
#        example[example$Common_Name == "Two",]$Daily_Site_Occupancy,
#        col = "blue")
# abline(v = max(example_enddate[example_enddate$Common_Name == "Two",]$Jd), col = "blue")
# head(example)
# tail(example)

# # -----------------------------------------------------
# # Filter for American Redstart 2013 Grid Grid_-100000x1300000
# # -----------------------------------------------------
# head(eBird_dat_enddate)
# nrow(eBird_dat_enddate)
# colnames(eBird_dat_enddate)
# eBird_dat_enddate_figure <- eBird_dat_enddate %>%
#   filter(Common_Name == "Palm Warbler",
#          Year == 2017,
#          GridID == "Grid_300000x700000")
# eBird_dat_enddate_figure$GridID <- str_replace_all(eBird_dat_enddate_figure$GridID, "_", " ")
# nrow(eBird_dat_enddate_figure)

# -----------------------------------------------------
# ANALYSIS AND PLOTTING BEGINS
# -----------------------------------------------------
## Initialize NLS routine
#
#--------------------------------------------------------------------------------
##  By year/species plot data and residual plots in pdf  1 page per species/year
##  Store results of fitting in data.frame
#--------------------------------------------------------------------------------

# -----------------------------------------------------
# WHOLE DATASET
# -----------------------------------------------------
# Tricky birds: Brown-headed Cowbird, Common Grackle, Hermit Thrush, Eastern Phoebe, Tree Swallow
# eBird_dat_enddate_Blocks_2 <- eBird_dat_enddate %>%
#   filter(Common_Name %in% c("Red-eyed Vireo", "Eastern Phoebe", "Tree Swallow"))

# Begin PDF writing
# datestub <- format(Sys.Date(), format="%Y%m%d");# Today's date as YYYYMMDD
# pdffile <- paste("Routine 7/Processing by species 2/", speciesname, "_ebird_MAD_Jd_", Jd_min, "_", Jd_max_write, "_w_enddate_", enddateadd, "_Occ_Thresh_", Occupancy_Threshold, "_2002-2019_", datestub, ".pdf", sep = "")
# pdf(file = pdffile, width = 8.5, height = 11)
datestub <- format(Sys.Date(), format="%Y%m%d");# Today's date as YYYYMMDD
pdffile <- paste("C:/Users/dpgil/Documents/1 U of T/eBird/eBird Project Processed Data/Routine 7/Simulation/Stationary/", "Simulation_Stationary", "_ebird_MAD_Jd_", Jd_min, "_", Jd_max_write, "_w_enddate_", enddateadd, "_Occ_Thresh_", Occupancy_Threshold, "_2002-2019_", datestub, ".pdf", sep = "")
pdf(file = pdffile, width = 8.5, height = 11)

# Rename
Blocks_2_enddate_stationary <- Blocks_2_enddate_stationary %>%
  rename("Daily_Site_Occupancy" = "Daily_Site_Occupancy_Stationary",
         "Daily_Site_Sightings" = "Daily_Site_Sightings_Stationary")

# arrivBlocks_2 <- Blocks_2_enddate_stationary[1:70,]

# # Try unloading package
# detach("package:EnvStats", unload=TRUE)


# Output <- nls2(Daily_Site_Occupancy ~ SSlogis(Jd, Asym, xmid, scal), data = arrivBlocks_2)

# # Try example data
# example <- read_csv("example.csv")
# example1 <- example[1:90,]
#
# Output <- nls2(Daily_Site_Occupancy ~ SSlogis(Jd, Asym, xmid, scal), data = example1)

# Designate layout()
par(mfrow = c(4,5))
layout(matrix(c(1, 2, 3,
                4, 5, 6,
                7, 8, 9,
                10, 11, 12,
                13, 14, 15,
                16, 17, 18), nrow = 6, byrow = TRUE))
layout.show(n = 19)

# Start function
ArrivDateSum <- ddply(Blocks_2_enddate_stationary,
                               c("GridID",
                                 "Common_Name",
                                 "Year"),
                      purrr::possibly(ArrivalDate, NULL))
nrow(ArrivDateSum)
ArrivDateSum$enddateadd <- factor(enddateadd)
ArrivDateSum$OccThresh <- factor(Occupancy_Threshold)
ArrivDateSum$Jd_max <- factor(Jd_max_write)

Blocks_2_enddate_stationary %>%
  select(GridID, Common_Name, Year) %>%
  distinct() %>%
  nrow()

# # Add columns to figure out if the model worked or not
# ArrivDateSum <- ArrivDateSum %>%
#   mutate(b_Est_Neg = ifelse(b_Estimate < 0, TRUE, FALSE),
#          MAD_CI_NA = ifelse(is.na(MAD_CI_pos) == TRUE, TRUE, FALSE),
#          MAD_Less_150 = ifelse(MAD < 150, TRUE, FALSE),
#          MAD_Over_80 = ifelse(MAD > 80, TRUE, FALSE),
#          Worked = ifelse(b_Est_Neg == TRUE &
#                            MAD_CI_NA == FALSE &
#                            MAD_Less_150 == TRUE &
#                            MAD_Over_80 == TRUE,
#                          TRUE, FALSE),
#          Special_Params = FALSE
#   )


# # close the device to complete the drawing to the pdf file
dev.off()
# Write out summary fit table to csv
csvfile <- paste("Routine 7/Simulation/Stationary/", "Simulation_Stationary", "_ebird_MAD_Jd_", Jd_min, "_", Jd_max_write, "_w_enddate_", enddateadd, "_Occ_Thresh_", Occupancy_Threshold, "_2002-2019_", datestub, ".csv", sep = "")
write_csv(ArrivDateSum, csvfile)



# -----------------------------------------------------
# END DATE CODE: Shifting
# -----------------------------------------------------
# -----------------------------------------------------
# Filter Jd ranges for master dataframe
# -----------------------------------------------------
# speciesname <- spp_names[28]
# speciesname <- str_replace_all(speciesname, " ", "_")
# for (j in 1:15) {
# enddateadd <- j
enddateadd <- 5
# for (k in c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) {
# Jd_min <- k
Jd_min <- unique(Blocks_2$MinJd)
# for (k in c(140, 150, 160, 170, 180, 190, 200)) {
# Jd_max <- (k + enddateadd)
Jd_max <- (180 + enddateadd)
# Jd_max For writing the dataframe
Jd_max_write <- (180)
# Jd_max_write <- (k)
# for (l in c(0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25)) {
# Occupancy_Threshold <- l
Occupancy_Threshold <- 0.1
Blocks_2 <- Blocks_2 %>%
  filter(Jd >= Jd_min & Jd <= Jd_max)
# Using Jd = 185 because the last 5 days will be removed due to the later code
Blocks_2_enddate_Shifting <- Blocks_2 %>%
  select(Year, Common_Name, GridID, Jd, Daily_Site_Sightings_Shifting, Daily_Site_Total, Daily_Site_Occupancy_Shifting)
# Apply end-date filter to each species to remove dates after which occupancy drop and the model becomes difficult to fit
# Add column that shows the maximum daily occupancy for each bird-year. Try a 10% of Max-Occ Threshold as a starting point
Blocks_2_enddate_Shifting <- Blocks_2_enddate_Shifting %>%
  group_by(GridID,
           Common_Name,
           Year) %>%
  mutate(Max_Occ = max(Daily_Site_Occupancy_Shifting)) %>%
  mutate(Max_Occ_Fraction = Max_Occ*Occupancy_Threshold) %>%
  ungroup()
# Add column that shows the date of the maximum daily occupancy for each bird-year
Blocks_2_enddate_Shifting <- Blocks_2_enddate_Shifting %>%
  group_by(GridID,
           Common_Name,
           Year) %>%
  mutate(Max_Occ_Jd = ifelse(Max_Occ == Daily_Site_Occupancy_Shifting, Jd, 1000)) %>%
  ungroup()
Blocks_2_enddate_Shifting <- Blocks_2_enddate_Shifting %>%
  group_by(GridID,
           Common_Name,
           Year) %>%
  mutate(Max_Occ_Jd_Min = min(Max_Occ_Jd)) %>%
  ungroup()
Blocks_2_enddate_Shifting$Max_Occ_Jd <- NULL
# Remove observations that occur after the date of maximum occupancy and have total site occupancy metrics that are lower than the total site occupancy fraction
# Add column that shows this (labelled "Departed", with 1 representing that the bird has departed from the landscape, and 0 representing that it still persists)
Blocks_2_enddate_Shifting <- Blocks_2_enddate_Shifting %>%
  mutate(Departed = ifelse(Jd > Max_Occ_Jd_Min & Daily_Site_Occupancy_Shifting <= Max_Occ_Fraction, 1, 0))
# Add column that shows the rolling sum of every certain number of days (start with 5) from the departed column. This essentially shows the last day during a certain day rolling sum where the rolling sum meets the interval that has been set to calculate the sum with (i.e. if the rolling sum is 5, and the interval is 5 days, then it means the birds havve been gone for 5 days)
n <- enddateadd
Blocks_2_enddate_Shifting$Departed_Roll_Sum <- c(rep_len(NA, n - 1), rowSums(embed(Blocks_2_enddate_Shifting$Departed, n)))
# Remove all observations after the first instance of the maximum sum
Blocks_2_enddate_Shifting <- Blocks_2_enddate_Shifting %>%
  group_by(GridID,
           Common_Name,
           Year) %>%
  mutate(First_Max_Sum = min(which(Departed_Roll_Sum == n | row_number() == n()))) %>%
  filter(row_number() <= First_Max_Sum) %>%
  ungroup()
# Remove additional rows of the maximum sum
# Blocks_2 this**********
Blocks_2_enddate_Shifting <- Blocks_2_enddate_Shifting %>%
  group_by(GridID,
           Common_Name,
           Year) %>%
  mutate(Spp_Year_Rows = n()) %>%
  filter(Spp_Year_Rows >= 5) %>%
  slice(1:(n()-n)) %>%
  ungroup()

# Blocks_2 <- Blocks_2_enddate_Shifting %>%
#   filter(ID == Blocks_2_IDs[1])

# # Plot
# # Daily_Site_Occupancy
# plot(example_enddate[example_enddate$Common_Name == "One",]$Jd,
#      example_enddate[example_enddate$Common_Name == "One",]$Daily_Site_Occupancy,
#      col = "red",
#      xlab = "Day of year",
#      ylab = "Daily site occupancy")
# points(example_enddate[example_enddate$Common_Name == "Two",]$Jd,
#        example_enddate[example_enddate$Common_Name == "Two",]$Daily_Site_Occupancy,
#        col = "blue")
# head(example_enddate)
# tail(example_enddate)
#
# # Plot
# # Daily_Site_Occupancy
# plot(example[example$Common_Name == "One",]$Jd,
#      example[example$Common_Name == "One",]$Daily_Site_Occupancy,
#      col = "red",
#      xlab = "Day of year",
#      ylab = "Daily site occupancy")
# points(example[example$Common_Name == "Two",]$Jd,
#        example[example$Common_Name == "Two",]$Daily_Site_Occupancy,
#        col = "blue")
# abline(v = max(example_enddate[example_enddate$Common_Name == "Two",]$Jd), col = "blue")
# head(example)
# tail(example)

# # -----------------------------------------------------
# # Filter for American Redstart 2013 Grid Grid_-100000x1300000
# # -----------------------------------------------------
# head(eBird_dat_enddate)
# nrow(eBird_dat_enddate)
# colnames(eBird_dat_enddate)
# eBird_dat_enddate_figure <- eBird_dat_enddate %>%
#   filter(Common_Name == "Palm Warbler",
#          Year == 2017,
#          GridID == "Grid_300000x700000")
# eBird_dat_enddate_figure$GridID <- str_replace_all(eBird_dat_enddate_figure$GridID, "_", " ")
# nrow(eBird_dat_enddate_figure)

# -----------------------------------------------------
# ANALYSIS AND PLOTTING BEGINS
# -----------------------------------------------------
## Initialize NLS routine
#
#--------------------------------------------------------------------------------
##  By year/species plot data and residual plots in pdf  1 page per species/year
##  Store results of fitting in data.frame
#--------------------------------------------------------------------------------

# -----------------------------------------------------
# WHOLE DATASET
# -----------------------------------------------------
# Tricky birds: Brown-headed Cowbird, Common Grackle, Hermit Thrush, Eastern Phoebe, Tree Swallow
# eBird_dat_enddate_Blocks_2 <- eBird_dat_enddate %>%
#   filter(Common_Name %in% c("Red-eyed Vireo", "Eastern Phoebe", "Tree Swallow"))

# Begin PDF writing
# datestub <- format(Sys.Date(), format="%Y%m%d");# Today's date as YYYYMMDD
# pdffile <- paste("Routine 7/Processing by species 2/", speciesname, "_ebird_MAD_Jd_", Jd_min, "_", Jd_max_write, "_w_enddate_", enddateadd, "_Occ_Thresh_", Occupancy_Threshold, "_2002-2019_", datestub, ".pdf", sep = "")
# pdf(file = pdffile, width = 8.5, height = 11)
datestub <- format(Sys.Date(), format="%Y%m%d");# Today's date as YYYYMMDD
pdffile <- paste("C:/Users/dpgil/Documents/1 U of T/eBird/eBird Project Processed Data/Routine 7/Simulation/Shifting/", "Simulation_Shifting", "_ebird_MAD_Jd_", Jd_min, "_", Jd_max_write, "_w_enddate_", enddateadd, "_Occ_Thresh_", Occupancy_Threshold, "_2002-2019_", datestub, ".pdf", sep = "")
pdf(file = pdffile, width = 8.5, height = 11)

# Rename
Blocks_2_enddate_Shifting <- Blocks_2_enddate_Shifting %>%
  rename("Daily_Site_Occupancy" = "Daily_Site_Occupancy_Shifting",
         "Daily_Site_Sightings" = "Daily_Site_Sightings_Shifting")

# arrivBlocks_2 <- Blocks_2_enddate_Shifting[1:70,]

# # Try unloading package
# detach("package:EnvStats", unload=TRUE)


# Output <- nls2(Daily_Site_Occupancy ~ SSlogis(Jd, Asym, xmid, scal), data = arrivBlocks_2)

# # Try example data
# example <- read_csv("example.csv")
# example1 <- example[1:90,]
#
# Output <- nls2(Daily_Site_Occupancy ~ SSlogis(Jd, Asym, xmid, scal), data = example1)

# Designate layout()
par(mfrow = c(4,5))
layout(matrix(c(1, 2, 3,
                4, 5, 6,
                7, 8, 9,
                10, 11, 12,
                13, 14, 15,
                16, 17, 18), nrow = 6, byrow = TRUE))
layout.show(n = 19)

# Start function
ArrivDateSum_2 <- ddply(Blocks_2_enddate_Shifting,
                      c("GridID",
                        "Common_Name",
                        "Year"),
                      purrr::possibly(ArrivalDate, NULL))
nrow(ArrivDateSum_2)
ArrivDateSum_2$enddateadd <- factor(enddateadd)
ArrivDateSum_2$OccThresh <- factor(Occupancy_Threshold)
ArrivDateSum_2$Jd_max <- factor(Jd_max_write)

# # Add columns to figure out if the model worked or not
# ArrivDateSum <- ArrivDateSum %>%
#   mutate(b_Est_Neg = ifelse(b_Estimate < 0, TRUE, FALSE),
#          MAD_CI_NA = ifelse(is.na(MAD_CI_pos) == TRUE, TRUE, FALSE),
#          MAD_Less_150 = ifelse(MAD < 150, TRUE, FALSE),
#          MAD_Over_80 = ifelse(MAD > 80, TRUE, FALSE),
#          Worked = ifelse(b_Est_Neg == TRUE &
#                            MAD_CI_NA == FALSE &
#                            MAD_Less_150 == TRUE &
#                            MAD_Over_80 == TRUE,
#                          TRUE, FALSE),
#          Special_Params = FALSE
#   )


# # close the device to complete the drawing to the pdf file
dev.off()
# Write out summary fit table to csv
csvfile <- paste("Routine 7/Simulation/Shifting/", "Simulation_Shifting", "_ebird_MAD_Jd_", Jd_min, "_", Jd_max_write, "_w_enddate_", enddateadd, "_Occ_Thresh_", Occupancy_Threshold, "_2002-2019_", datestub, ".csv", sep = "")
write_csv(ArrivDateSum_2, csvfile)












# -----------------------------------------------------
# Finish code
# -----------------------------------------------------



# # -----------------------------------------------------
# # Take a look at the median value
# # -----------------------------------------------------
# # Use uncount to expand the dates
# example_expand <- example_enddate %>%
#   uncount(Daily_Site_Sightings)
# # Calculate the median days
# median_One <- median(example_expand[example_expand$Common_Name == "One", ]$Jd)
# median_Two <- median(example_expand[example_expand$Common_Name == "Two", ]$Jd)
# # plot the dates
# plot(example_enddate[example_enddate$Common_Name == "One", ]$Jd,
#      example_enddate[example_enddate$Common_Name == "One", ]$Daily_Site_Sightings,
#      col = "red",
#      type = "h",
#      xlab = "Day of year",
#      ylab = "Daily sightings count",
#      main = "Species One")
# abline(v = median_One, col = "black", lty = 3, lwd = 3)
# abline(v = ArrivDateSum$MAD[1], col = "black", lty = 2, lwd = 3)
# plot(example_enddate[example_enddate$Common_Name == "Two", ]$Jd,
#      example_enddate[example_enddate$Common_Name == "Two", ]$Daily_Site_Sightings,
#      col = "blue",
#      type = "h",
#      xlab = "Day of year",
#      ylab = "Daily sightings count",
#      main = "Species Two")
# abline(v = median_Two, col = "black", lty = 3, lwd = 3)
# abline(v = ArrivDateSum$MAD[2], col = "black", lty = 2, lwd = 3)

# hist(example_expand[example_expand$Common_Name == "One", ]$Jd,
#      xlab = "Day of year",
#      ylab = "Daily sightings count",
#      main = "Species One")
# hist(example_expand[example_expand$Common_Name == "Two", ]$Jd,
#      xlab = "Day of year",
#      ylab = "Daily sightings count",
#      main = "Species Two")




# }
# }
# }
# }
# }

# eBird_dat_enddate_figure_2 <- eBird_Dat %>%
#   filter(Common_Name == "Palm Warbler",
#          Year == 2017,
#          GridID == "Grid_300000x700000",
#          Jd >= 80,
#          Jd <= 180)
#
# Output <- nls2(Daily_Site_Occupancy ~ SSlogis(Jd, Asym, xmid, scal), data = eBird_dat_enddate_figure)
#
# plot(x = eBird_dat_enddate_figure_2$Jd,
#      y = eBird_dat_enddate_figure_2$Daily_Site_Occupancy,
#      xlab = "Julian Day",
#      ylab = "Daily Site Occupancy",
#      ylim = c(0, 1.0),
#      cex = 1
# )
#
# r <- range(eBird_dat_enddate_figure_2$Jd)
# xNew <- seq(r[1],ArrivDateSum$MaxJd,length.out = 200)
# yNew <- predict(Output, list(Jd = xNew))
# lines(xNew,yNew)
# abline(v = ArrivDateSum$MAD, lty=2)
# rect(xleft = ArrivDateSum$MAD_CI_neg,
#      ybottom = -0.1,
#      xright = ArrivDateSum$MAD_CI_pos,
#      ytop = 1.1,
#      border = NA,
#      col = rgb(0,0,0,alpha=0.2))
#
#
# mtext(paste(eBird_dat_enddate_figure_2$Common_Name, " - ", eBird_dat_enddate_figure_2$Year, " - ", eBird_dat_enddate_figure_2$GridID), cex = 1)
# abline(v = 145, lty=1)
#
# # Redo plot
# ggplot() +
#   geom_point(data = eBird_dat_enddate_figure_2, aes(x = Jd,
#      y = Daily_Site_Occupancy), shape = 1, size = 2) +
#      xlab("Julian day") +
#      ylab ("Daily site occupancy") +
#      scale_y_continuous(limits = c(0,1.0)
#                         # , expand = c(0, 0)
#                         ) +
#   scale_x_continuous(breaks = c(80, 100, 120, 140, 160, 180)
#     # expand = c(0, 0)
#     ) +
#   geom_rect(aes(xmin = ArrivDateSum$MAD_CI_neg, xmax = ArrivDateSum$MAD_CI_pos, ymin = -Inf, ymax = Inf ), fill = "gray", alpha =0.5) +
#   geom_line(aes(x = xNew, y = yNew)) +
#   geom_vline(xintercept = ArrivDateSum$MAD, linetype = "dashed") +
#   geom_vline(xintercept = 145, linetype = "solid") +
#   # geom_segment(aes(x = 80, y = coef(Output_summary)[1], xend = 145, yend = coef(Output_summary)[1] ), linetype = "dotted", alpha = 0.5) +
#   # geom_hline(yintercept = coef(Output_summary)[1], linetype = "dotted", alpha = 0.5) +
#   theme_classic() +
#   theme(text = element_text(size = 14))
#










# # # -------------------------------------------
# # # Create simulation data
# # # -------------------------------------------
# example <- data.frame(
#   # Species One
#   Year = 2002,
#   Common_Name = c(rep("One", length(c(80:180))), rep("Two", length(c(80:180)))),
#   GridID = "example",
#   Jd = rep(c(80:180), 2),
#   Daily_Site_Sightings = c(pnorm(q = (seq(80, 180, 1)),
#                             mean = 130,
#                             sd = 8),
#                       dnorm(x = (seq(80, 180, 1)),
#                             mean = 130,
#                             sd = 5)*10),
#   Daily_Site_Total = rnorm(n = length(rep(c(80:180), 2)),
#                            mean = 200,
#                            sd = 10))
# # multiply the values by 100
# example$Daily_Site_Sightings <- example$Daily_Site_Sightings * 100
# # Round to whole numbers
# example$Daily_Site_Sightings <- round(example$Daily_Site_Sightings , digits = 0)
# example$Daily_Site_Total <- round(example$Daily_Site_Total , digits = 0)
# # Introduce some variability
# example$Daily_Site_Sightings <- example$Daily_Site_Sightings + round(rnorm(n = length(c(80:180)),
#                                                                  mean = 10,
#                                                                  sd = 10), digits = 0)
# # Subtract some from the initial days in the series
# example$Daily_Site_Sightings <- ifelse(example$Jd < 100, example$Daily_Site_Sightings - 20, example$Daily_Site_Sightings)
# # Correct negative values to 0
# example$Daily_Site_Sightings <- ifelse(example$Daily_Site_Sightings < 0, 0, example$Daily_Site_Sightings)
#
# # Calculate Daily_Site_Occupancy
# example$Daily_Site_Occupancy <- example$Daily_Site_Sightings / example$Daily_Site_Total
# summary(example$Daily_Site_Occupancy)
#
# # Plot example data
# # Daily_Site_Sightings
# plot(example[example$Common_Name == "One",]$Jd,
#      example[example$Common_Name == "One",]$Daily_Site_Sightings,
#      col = "red",
#      xlab = "Day of year",
#      ylab = "Daily sightings")
# points(example[example$Common_Name == "Two",]$Jd,
#        example[example$Common_Name == "Two",]$Daily_Site_Sightings,
#        col = "blue")
# # Daily_Site_Total
# plot(example[example$Common_Name == "One",]$Jd,
#      example[example$Common_Name == "One",]$Daily_Site_Total,
#      col = "red",
#      xlab = "Day of year",
#      ylab = "Daily sightings")
# points(example[example$Common_Name == "Two",]$Jd,
#        example[example$Common_Name == "Two",]$Daily_Site_Total,
#        col = "blue")
# # Daily_Site_Occupancy
# plot(example[example$Common_Name == "One",]$Jd,
#      example[example$Common_Name == "One",]$Daily_Site_Occupancy,
#      col = "red",
#      xlab = "Day of year",
#      ylab = "Daily site occupancy")
# points(example[example$Common_Name == "Two",]$Jd,
#        example[example$Common_Name == "Two",]$Daily_Site_Occupancy,
#        col = "blue")
# head(example)
# tail(example)

#
# # # ------------------------------------------------
# # # Read in grid datasets
# # # ------------------------------------------------
# # eBird_Dat <- read_csv("Routine 6 - Final Grids/eBird_Daily_Occupancy_final_grids.csv")
# # length(unique(eBird_Dat$GridID))
# # Grid_Blocks_2 <- eBird_Dat %>%
# #   filter(GridID == "Grid_-100000x300000")
# # # -----------------------------------------------------
# # # Go through all species
# # # -----------------------------------------------------
# # eBird_Dat_Spp <- eBird_Dat %>%
# #   filter(Common_Name == spp_names[28])
# # length(unique(eBird_Dat_Spp$GridID))

# # Blocks_2 complement function
# complement <- function(y, rho, x) {
#   if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
#   y.perp <- residuals(lm(x ~ y))
#   rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
# }

# # Add column
# Blocks_2$Daily_Site_Sightings_Stationary <- NA
# Blocks_2$Daily_Site_Sightings_Shifting <- NA
# Blocks_2$Daily_Site_Sightings_Stationary <- as.numeric(Blocks_2$Daily_Site_Sightings_Stationary)
# Blocks_2$Daily_Site_Sightings_Shifting <- as.numeric(Blocks_2$Daily_Site_Sightings_Shifting)

# # Blocks_2
# i = 1
# i0 = i - 1
# ID_i <- unique(Blocks_2$ID)[i]
# df_i <- filter(Blocks_2, ID == ID_i)
# set.seed(i)
# vec_i <- ifelse(df_i$MaxJd == 180,
#                 pnorm(q = seq(unique(df_i$MinJd), unique(df_i$MaxJd), 1),
#                             mean = unique(df_i$MAD_Stationary),
#                             sd = middle_sd_mean)*100,
#                 dnorm(x = seq(unique(df_i$MinJd), unique(df_i$MaxJd), 1),
#                             mean = unique(df_i$MAD_Stationary),
#                             sd = middle_sd_mean)*3000)
# # Multiply by Daily_Site_Total values converted to 0 and 1 (anything positive)
# site_total_1_0_i <- ifelse(df_i$Daily_Site_Total > 0, 1, 0)
# vec_i <- vec_i * site_total_1_0_i
# # Only x days can be positive, where x is JdDetect
# JdDetect_i <- unique(df_i$JdDetect)
# JdDetect_inv_i <- length(vec_i) - JdDetect_i
# vec_detect_i <- data.frame(vec_i = vec_i)
# vec_detect_i_length <- length(vec_i[vec_i > 0])
# # flag positive days
# vec_detect_i$detect_flag <- ifelse(vec_detect_i$vec_i > 0, TRUE, FALSE)
# # Add row_id
# vec_detect_i$row_id <- seq(1, nrow(vec_detect_i), 1)
# # Get difference between days detected and days that are supposed to be detected
# vec_detect_i_diff <- JdDetect_inv_i - JdDetect_i
# # Break data.frame into TRUE and FALSE
# vec_detect_i_TRUE <- vec_detect_i %>%
#   filter(detect_flag == TRUE)
# vec_detect_i_FALSE <- vec_detect_i %>%
#   filter(detect_flag == FALSE)
# # Create a subtraction vector to enforce a set amount of detections
# set.seed(i)
# subtract_vec_i <- sample(c(rep(0, vec_detect_i_diff), rep(1, JdDetect_i)))
# # Multiply the vectors together
# vec_detect_i_TRUE$vec_i <- vec_detect_i_TRUE$vec_i * subtract_vec_i
# # bind the TRUE and FALSE data.frames together and reorder
# vec_detect_i <- rbind(vec_detect_i_TRUE, vec_detect_i_FALSE)
# vec_detect_i <- vec_detect_i %>%
#   arrange(-desc(row_id))
# # Extract vec_i
# vec_i <- vec_detect_i$vec_i
# # Get vec_i sum
# # Correct vec_i to Sighting_Sum
# sighting_sum_i <- unique(df_i$Sighting_Sum)
# correction_factor <- sighting_sum_i/sum(vec_i)
# vec_i <- vec_i * correction_factor
# # Introduce some random variability to the positive values
# set.seed(i)
# vec_i <- ifelse(vec_i > 0, vec_i + rnorm(n = length(c(80:180)), mean = mean(vec_i[vec_i > 0]), sd = mean(vec_i[vec_i > 0])/10), vec_i)
# # Correct negative values to zero
# vec_i <- ifelse(vec_i < 0, 0, vec_i)
# # Correct 1 values to 2
# vec_i <- ifelse(vec_i > 0 & vec_i < 1.5, 2, vec_i)
# # How many non-zero days are there
# length_above_zero <- length(vec_i[vec_i > 0])
# # If length above zero is odd, change the Daily_Site_Total that is the lowest to zero
# # Index the minimum spot
# min_spot <- which(vec_i == min(vec_i[vec_i > 0]))
# vec_i[min_spot] <- ifelse(length_above_zero %% 2 == 0, vec_i[min_spot], 0)
# # Get vec_i sum again
# # How many sites is this?
# correction_factor_2 <- sighting_sum_i / sum(vec_i)
# # Multiply the vector by the correction factor again
# vec_i <- vec_i * correction_factor_2
# Round the vector
# vec_i <- round(vec_i)
# # What is the absolute difference between the site numbers?
# vec_diff <- sighting_sum_i - sum(vec_i)
# # Add/subtract this to a days when observations are greater than 1
# vec_detect_i <- data.frame(vec_i = vec_i)
# # flag days above 1
# vec_detect_i$above_plus_flag <- ifelse(vec_detect_i$vec_i > 1, TRUE, FALSE)
# # Add row_id
# vec_detect_i$row_id <- seq(1, nrow(vec_detect_i), 1)
# # Break data.frame into TRUE and FALSE
# vec_detect_i_TRUE <- vec_detect_i %>%
#   filter(above_plus_flag == TRUE)
# vec_detect_i_FALSE <- vec_detect_i %>%
#   filter(above_plus_flag == FALSE)
# set.seed(i)
# vec_detect_i_TRUE_2 <- vec_detect_i_TRUE %>%
#   sample_n(abs(vec_diff))
# vec_detect_i_TRUE <- vec_detect_i_TRUE %>%
#   filter(!row_id %in% vec_detect_i_TRUE_2$row_id)
# # Subtract
# vec_detect_i_TRUE_2 <- vec_detect_i_TRUE_2 %>%
#   mutate(vec_i = vec_i + vec_diff/abs(vec_diff))
# # Bind together and sort
# vec_detect_i <- rbind(vec_detect_i_TRUE,
#                       vec_detect_i_TRUE_2,
#                       vec_detect_i_FALSE)
# vec_detect_i <- vec_detect_i %>%
#   arrange(-desc(row_id))
# vec_i <- vec_detect_i$vec_i
# # Add to Daily_Site_Total column (23)
# start_ind <- 101*i0 + 1
# end_ind <- 101*i0 + 1 + 100
# Blocks_2 <- Blocks_2
# Blocks_2[c(start_ind:end_ind), 26] <- vec_i\\

# Blocks_2ing <- Blocks_2 %>%
#   mutate(Daily_round = round(Daily_Site_Total)) %>%
#   mutate(Daily_flag = if_else(Daily_round != Daily_Site_Total, TRUE, FALSE))
# Blocks_2ing %>%
#   select(ID, Daily_flag) %>%
#   distinct() %>%
#   count(Daily_flag)
# Blocks_2ing2 <- Blocks_2ing %>%
#   filter(Daily_flag == TRUE)
# nrow(Blocks_2ing2)




# looking <- Blocks_2 %>%
#   filter(ID %in% "Species_14_Grid_1_11")
# ggplot(looking) +
#   geom_line(aes(Jd, Daily_Site_Total), colour = "blue") +
#   geom_line(aes(Jd, Daily_Site_Sightings_Shifting), colour = "red") +
#   expand_limits(y = 0) +
#   scale_x_continuous(limits = c(80, 185))
# ggplot(looking) +
#   geom_line(aes(Jd, Daily_Site_Occupancy_Shifting), colour = "red") +
#   expand_limits(y = 0) +
#   scale_x_continuous(limits = c(80, 185))

total_end <- Sys.time()
total_end - total_start
