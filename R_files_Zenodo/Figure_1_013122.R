# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 7
# ------------------------------------------------
# Load Libraries
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(nls2)
library(nlstools)
library(data.table)

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
# ------------------------------------------------
# Finish function subroutines
# ------------------------------------------------

# ------------------------------------------------
# Read in grid datasets
# ------------------------------------------------
eBird_Dat <- fread("Routine 6/Routine 6 - merged files/eBird_Daily_Occupancy_all_grids.csv")
length(unique(eBird_Dat$GridID))

# -----------------------------------------------------
# Go through all species
# -----------------------------------------------------
eBird_Dat_Spp <- eBird_Dat %>% 
  filter(Common_Name == "Palm Warbler")
length(unique(eBird_Dat_Spp$GridID))
rm(eBird_Dat) ; gc()
# -----------------------------------------------------
# END DATE CODE
# -----------------------------------------------------
# -----------------------------------------------------
# Filter Jd ranges for master dataframe
# -----------------------------------------------------
speciesname <- "Palm Warbler"
speciesname <- str_replace_all(speciesname, " ", "_")
  # for (j in 1:15) {
# enddateadd <- j
enddateadd <- 5
# for (k in c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) {
# Jd_min <- k
Jd_min <- 60
# for (k in c(140, 150, 160, 170, 180, 190, 200)) {
# Jd_max <- (k + enddateadd)
Jd_max <- (180 + enddateadd)
# Jd_max For writing the dataframe
Jd_max_write <- (180)
# Jd_max_write <- (k)
# for (l in c(0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25)) {
# Occupancy_Threshold <- l
Occupancy_Threshold <- 0.1
eBird_dat_filter <- eBird_Dat_Spp %>%
  filter(Jd >= Jd_min & Jd <= Jd_max)
# Using Jd = 185 because the last 5 days will be removed due to the later code
eBird_dat_enddate <- eBird_dat_filter
# Apply end-date filter to each species to remove dates after which occupancy drop and the model becomes difficult to fit
# Add column that shows the maximum daily occupancy for each bird-year. Try a 10% of Max-Occ Threshold as a starting point
eBird_dat_enddate <- eBird_dat_enddate %>%
  group_by(GridID,
    # Common_Name,
           Year) %>%
  mutate(Max_Occ = max(Daily_Site_Occupancy)) %>%
  mutate(Max_Occ_Fraction = Max_Occ*Occupancy_Threshold) %>%
  ungroup()
# Add column that shows the date of the maximum daily occupancy for each bird-year
eBird_dat_enddate <- eBird_dat_enddate %>%
  group_by(GridID,
    # Common_Name,
    Year) %>%
  mutate(Max_Occ_Jd = ifelse(Max_Occ == Daily_Site_Occupancy, Jd, 1000)) %>%
  ungroup()
eBird_dat_enddate <- eBird_dat_enddate %>%
  group_by(GridID,
    # Common_Name,
    Year) %>%
  mutate(Max_Occ_Jd_Min = min(Max_Occ_Jd)) %>%
  ungroup()
eBird_dat_enddate$Max_Occ_Jd <- NULL
# Remove observations that occur after the date of maximum occupancy and have total site occupancy metrics that are lower than the total site occupancy fraction
# Add column that shows this (labelled "Departed", with 1 representing that the bird has departed from the landscape, and 0 representing that it still persists)
eBird_dat_enddate <- eBird_dat_enddate %>%
  mutate(Departed = ifelse(Jd > Max_Occ_Jd_Min & Daily_Site_Occupancy <= Max_Occ_Fraction, 1, 0))
# Add column that shows the rolling sum of every certain number of days (start with 5) from the departed column. This essentially shows the last day during a certain day rolling sum where the rolling sum meets the interval that has been set to calculate the sum with (i.e. if the rolling sum is 5, and the interval is 5 days, then it means the birds havve been gone for 5 days)
n <- enddateadd
eBird_dat_enddate$Departed_Roll_Sum <- c(rep_len(NA, n - 1), rowSums(embed(eBird_dat_enddate$Departed, n)))
# Remove all observations after the first instance of the maximum sum
eBird_dat_enddate <- eBird_dat_enddate %>%
  group_by(GridID,
    # Common_Name,
           Year) %>%
  mutate(First_Max_Sum = min(which(Departed_Roll_Sum == n | row_number() == n()))) %>%
  filter(row_number() <= First_Max_Sum) %>%
  ungroup()
# Remove additional rows of the maximum sum
# Test this**********
eBird_dat_enddate <- eBird_dat_enddate %>%
  group_by(GridID,
    # Common_Name,
    Year) %>%
  mutate(Spp_Year_Rows = n()) %>%
  filter(Spp_Year_Rows >= 5) %>%
  slice(1:(n()-n)) %>%
  ungroup()


# -----------------------------------------------------
# Filter for Palm Warbler 2017 Grid Grid_300000x700000
# -----------------------------------------------------
head(eBird_dat_enddate)
nrow(eBird_dat_enddate)
colnames(eBird_dat_enddate)
eBird_dat_enddate_figure <- eBird_dat_enddate %>%
  filter(Common_Name == "Palm Warbler",
         Year == 2017,
         GridID == "Grid_300000x700000")
eBird_dat_enddate_figure$GridID <- str_replace_all(eBird_dat_enddate_figure$GridID, "_", " ")
nrow(eBird_dat_enddate_figure)

# -----------------------------------------------------
# ANALYSIS AND PLOTTING BEGINS
# -----------------------------------------------------
# Initialize NLS routine

#--------------------------------------------------------------------------------
#  By year/species plot data and residual plots in pdf  1 page per species/year
#  Store results of fitting in data.frame
#--------------------------------------------------------------------------------
# Start function
ArrivDateSum <- plyr::ddply(eBird_dat_enddate_figure,
                               c("GridID",
                                 # "Common_Name",
                                 "Year"),
                            purrr::possibly(ArrivalDate, NULL))
ArrivDateSum$enddateadd <- factor(enddateadd)
ArrivDateSum$OccThresh <- factor(Occupancy_Threshold)
ArrivDateSum$Jd_max <- factor(Jd_max_write)

# close the device to complete the drawing to the pdf file
# dev.off()
# Not writing file because this is just for the figure
# datestub <- "today"
# # Write out summary fit table to csv
# csvfile <- paste("Routine 7/Processing by species 2/", speciesname, "_ebird_MAD_Jd_", Jd_min, "_", Jd_max_write, "_w_enddate_", enddateadd, "_Occ_Thresh_", Occupancy_Threshold, "_2002-2019_", datestub, ".csv", sep = "")
# write_csv(ArrivDateSum, csvfile)

# -----------------------------------------------------
# Finish code
# -----------------------------------------------------

eBird_dat_enddate_figure_2 <- eBird_Dat_Spp %>% 
  filter(Common_Name == "Palm Warbler",
         Year == 2017,
         GridID == "Grid_300000x700000",
         Jd >= 60, 
         Jd <= 180)

Output <- nls2(Daily_Site_Occupancy ~ SSlogis(Jd, Asym, xmid, scal), data = eBird_dat_enddate_figure)

# Base R plot to show in sequence how each piece fits together
plot(x = eBird_dat_enddate_figure_2$Jd, 
     y = eBird_dat_enddate_figure_2$Daily_Site_Occupancy,
     xlab = "Ordinal date",
     ylab = "Daily Site Occupancy",
     ylim = c(0, 1.0),
     cex = 1
)

r <- range(eBird_dat_enddate_figure_2$Jd)
xNew <- seq(r[1],ArrivDateSum$MaxJd,length.out = 200)
yNew <- predict(Output, list(Jd = xNew))
lines(xNew,yNew)
abline(v = ArrivDateSum$MAD, lty=2)
rect(xleft = ArrivDateSum$MAD_CI_neg,
     ybottom = -0.1,
     xright = ArrivDateSum$MAD_CI_pos,
     ytop = 1.1,
     border = NA,
     col = rgb(0,0,0,alpha=0.2))


mtext(paste(eBird_dat_enddate_figure_2$Common_Name, " - ", eBird_dat_enddate_figure_2$Year, " - ", eBird_dat_enddate_figure_2$GridID), cex = 1)
abline(v = 145, lty=1)

# ggplot plot for Figure 1
ggplot() +
  geom_point(data = eBird_dat_enddate_figure_2, aes(x = Jd, 
     y = Daily_Site_Occupancy), shape = 1, size = 2) +
     xlab("Ordinal date") +
     ylab ("Daily site occupancy") +
     scale_y_continuous(limits = c(0,1.0)
                        ) +
  scale_x_continuous(breaks = c(60, 80, 100, 120, 140, 160, 180)
    ) +
  geom_rect(aes(xmin = ArrivDateSum$MAD_CI_neg, xmax = ArrivDateSum$MAD_CI_pos, ymin = -Inf, ymax = Inf ), fill = "gray", alpha =0.5) +
  geom_line(aes(x = xNew, y = yNew)) +
  geom_vline(xintercept = ArrivDateSum$MAD, linetype = "dashed") +
  geom_vline(xintercept = 145, linetype = "solid") +
  theme_classic() +
  theme(text = element_text(size = 14))

