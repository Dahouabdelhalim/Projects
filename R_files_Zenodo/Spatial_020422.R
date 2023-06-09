# ---------------------------------------------
# Analyzing mean arrival dates
# ---------------------------------------------
# Load libraries
library(tidyverse)
library(nlme)
library(ape)
library(MuMIn)
library(adespatial)
library(vegan)
library(ade4)
library(adegraphics)
library(spdep)

# Set Working Directory
setwd("~/Data repository")


# # Additional functions
# source("C:/Users/dpgil/Documents/1 U of T/eBird/Literature/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/plot.links.R")
# source("C:/Users/dpgil/Documents/1 U of T/eBird/Literature/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/sr.value.R")
# source("C:/Users/dpgil/Documents/1 U of T/eBird/Literature/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/quickMEM.R")
# # source("C:/Users/dpgil/Documents/1 U of T/eBird/Literature/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/scalog.R")



# ---------------------------------------------
# Read in data
# ---------------------------------------------
eBird <- read_csv("Routine 7/Merged Output/MAD_Output_All_Filter_1.csv")



# ---------------------------------------------
# Transform and subset data
# ---------------------------------------------
# First, define the irregular lattice
# Set radius
radius = 100000
# Extract lat long from the grids
eBird$lat <- str_extract_all(eBird$GridID, "([-]*[0-9]+[x])")
eBird$lat <- as.numeric(str_remove_all(eBird$lat, "[x]"))
eBird$long <- as.numeric(str_extract_all(eBird$GridID, "[0-9]+$"))
# Add jitter
set.seed(456)
head(eBird$lat)
eBird$lat_orig <- eBird$lat
eBird$lat <- jitter(eBird$lat, amount = 1)
head(eBird$lat)
set.seed(456)
eBird$long_orig <- eBird$long
eBird$long <- jitter(eBird$long, amount = 01)
eBird_coords <- data.frame(x = eBird$long, y = eBird$lat)
head(eBird_coords)
# Select only columns that are needed for analysis
eBird <- eBird %>% 
  select(GridID, Year, Common_Name, MAD, lat, long, lat_orig, long_orig)
# Change GridID to a number
Grid_numbers <- data.frame(GridID = unique(eBird$GridID))
Grid_numbers$ID <- 1:nrow(Grid_numbers)
eBird <- left_join(eBird, Grid_numbers)
# eBird$GridID <- NULL

# Get a unique ID for each species
birds <- data.frame(Common_Name = unique(eBird$Common_Name))
birds$bird_int <- 1:nrow(birds)
# Get a subset for each species
eBird_1 <- eBird %>% filter(Common_Name == birds$Common_Name[1])
eBird_2 <- eBird %>% filter(Common_Name == birds$Common_Name[2])
eBird_3 <- eBird %>% filter(Common_Name == birds$Common_Name[3])
eBird_4 <- eBird %>% filter(Common_Name == birds$Common_Name[4])
eBird_5 <- eBird %>% filter(Common_Name == birds$Common_Name[5])
eBird_6 <- eBird %>% filter(Common_Name == birds$Common_Name[6])
eBird_7 <- eBird %>% filter(Common_Name == birds$Common_Name[7])
eBird_8 <- eBird %>% filter(Common_Name == birds$Common_Name[8])
eBird_9 <- eBird %>% filter(Common_Name == birds$Common_Name[9])
eBird_10 <- eBird %>% filter(Common_Name == birds$Common_Name[10])
eBird_11 <- eBird %>% filter(Common_Name == birds$Common_Name[11])
eBird_12 <- eBird %>% filter(Common_Name == birds$Common_Name[12])
eBird_13 <- eBird %>% filter(Common_Name == birds$Common_Name[13])
eBird_14 <- eBird %>% filter(Common_Name == birds$Common_Name[14])
eBird_15 <- eBird %>% filter(Common_Name == birds$Common_Name[15])
eBird_16 <- eBird %>% filter(Common_Name == birds$Common_Name[16])
eBird_17 <- eBird %>% filter(Common_Name == birds$Common_Name[17])
eBird_18 <- eBird %>% filter(Common_Name == birds$Common_Name[18])
eBird_19 <- eBird %>% filter(Common_Name == birds$Common_Name[19])
eBird_20 <- eBird %>% filter(Common_Name == birds$Common_Name[20])
eBird_21 <- eBird %>% filter(Common_Name == birds$Common_Name[21])
eBird_22 <- eBird %>% filter(Common_Name == birds$Common_Name[22])
eBird_23 <- eBird %>% filter(Common_Name == birds$Common_Name[23])
eBird_24 <- eBird %>% filter(Common_Name == birds$Common_Name[24])
eBird_25 <- eBird %>% filter(Common_Name == birds$Common_Name[25])
eBird_26 <- eBird %>% filter(Common_Name == birds$Common_Name[26])
eBird_27 <- eBird %>% filter(Common_Name == birds$Common_Name[27])
eBird_28 <- eBird %>% filter(Common_Name == birds$Common_Name[28])
eBird_29 <- eBird %>% filter(Common_Name == birds$Common_Name[29])
eBird_30 <- eBird %>% filter(Common_Name == birds$Common_Name[30])
eBird_31 <- eBird %>% filter(Common_Name == birds$Common_Name[31])
eBird_32 <- eBird %>% filter(Common_Name == birds$Common_Name[32])
eBird_33 <- eBird %>% filter(Common_Name == birds$Common_Name[33])
eBird_34 <- eBird %>% filter(Common_Name == birds$Common_Name[34])
eBird_35 <- eBird %>% filter(Common_Name == birds$Common_Name[35])
eBird_36 <- eBird %>% filter(Common_Name == birds$Common_Name[36])
eBird_37 <- eBird %>% filter(Common_Name == birds$Common_Name[37])
eBird_38 <- eBird %>% filter(Common_Name == birds$Common_Name[38])
eBird_39 <- eBird %>% filter(Common_Name == birds$Common_Name[39])
eBird_40 <- eBird %>% filter(Common_Name == birds$Common_Name[40])
eBird_41 <- eBird %>% filter(Common_Name == birds$Common_Name[41])
eBird_42 <- eBird %>% filter(Common_Name == birds$Common_Name[42])
eBird_43 <- eBird %>% filter(Common_Name == birds$Common_Name[43])
eBird_44 <- eBird %>% filter(Common_Name == birds$Common_Name[44])
# test
unique(eBird_10$Common_Name) ; birds$Common_Name[10]
unique(eBird_29$Common_Name) ; birds$Common_Name[29]



# -----------------------------------------------
# Fit mixed model
# -----------------------------------------------
# eBird_1
# Get spatial information
eBird_1
eBird_1_coords <- cbind(eBird_1$long, eBird_1$lat)
# eBird_1_coords <- unique(eBird_1_coords)
samples_dist_eBird_1 <- as.matrix(dist(eBird_1_coords))
samples_dist_inv_eBird_1 <- 1/samples_dist_eBird_1
is.na(samples_dist_inv_eBird_1) <- sapply(samples_dist_inv_eBird_1, is.infinite)
samples_dist_inv_eBird_1[is.na(samples_dist_inv_eBird_1)] <- 0
# Fit model
eBird_1_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_1)
summary(eBird_1_mm)
semivario_eBird_1 <- Variogram(eBird_1_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_1, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_1)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_1$MAD, samples_dist_inv_eBird_1, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_1_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_1)
summary(eBird_1_mm_exp)
# gaussian
eBird_1_mm_gaus <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_1)
summary(eBird_1_mm_gaus)
# spherical
eBird_1_mm_spher <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_1)
summary(eBird_1_mm_spher)
# linear
eBird_1_mm_lin <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_1)
summary(eBird_1_mm_lin)
# ratio
eBird_1_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_1)
summary(eBird_1_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_1_mm,
  eBird_1_mm_exp,
  eBird_1_mm_gaus,
  eBird_1_mm_spher,
  # eBird_1_mm_lin,
  eBird_1_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_1_mm), residuals(eBird_1_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_1_mm_spher), residuals(eBird_1_mm_spher))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_1_mm))
qqline(residuals(eBird_1_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_1_mm_spher))
# qqline(residuals(eBird_1_mm_spher))
# Semivariogram of normal model
semivario_eBird_1_b <- Variogram(eBird_1_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_1_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_1_best <- Variogram(eBird_1_mm_spher, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_1_best, smooth = TRUE)
# Summary of normal model
summary(eBird_1_mm)
# # Summary of best model
# summary(eBird_1_mm_spher)
# Get coefficients
intervals(eBird_1_mm)


coef_df <- data.frame(
  Common_Name = unique(eBird_1$Common_Name),
  # Model
  Model = c(
    "mm",
    "exp",
    "gaus",
    "rat",
    "lin",
    "spher"
  ),
  # Slope
  Slope = c(
    eBird_1_mm$coefficients$fixed[2],
    eBird_1_mm_exp$coefficients$fixed[2],
    eBird_1_mm_gaus$coefficients$fixed[2],
    eBird_1_mm_rat$coefficients$fixed[2],
    NA,
    eBird_1_mm_spher$coefficients$fixed[2]
  ),
  # SE minus
  SE_minus = c(
    eBird_1_mm$coefficients$fixed[2] - 1.96*summary(eBird_1_mm)$tTable[2,2],
    eBird_1_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_1_mm_exp)$tTable[2,2],
    eBird_1_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_1_mm_gaus)$tTable[2,2],
    eBird_1_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_1_mm_rat)$tTable[2,2],
    NA,
    eBird_1_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_1_mm_spher)$tTable[2,2]
  ),
  # SE plus
  SE_plus = c(
    eBird_1_mm$coefficients$fixed[2] + 1.96*summary(eBird_1_mm)$tTable[2,2],
    eBird_1_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_1_mm_exp)$tTable[2,2],
    eBird_1_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_1_mm_gaus)$tTable[2,2],
    eBird_1_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_1_mm_rat)$tTable[2,2],
    NA,
    eBird_1_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_1_mm_spher)$tTable[2,2]
  ),
  # p value
  p_value = c(
    summary(eBird_1_mm)$tTable[2,5],
    summary(eBird_1_mm_exp)$tTable[2,5],
    summary(eBird_1_mm_gaus)$tTable[2,5],
    summary(eBird_1_mm_rat)$tTable[2,5],
    # summary(eBird_1_mm_lin)$tTable[2,5],
    NA,
    summary(eBird_1_mm_spher)$tTable[2,5]
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))spdep



# eBird_2
# Get spatial information
eBird_2
eBird_2_coords <- cbind(eBird_2$long, eBird_2$lat)
# eBird_2_coords <- unique(eBird_2_coords)
samples_dist_eBird_2 <- as.matrix(dist(eBird_2_coords))
samples_dist_inv_eBird_2 <- 1/samples_dist_eBird_2
is.na(samples_dist_inv_eBird_2) <- sapply(samples_dist_inv_eBird_2, is.infinite)
samples_dist_inv_eBird_2[is.na(samples_dist_inv_eBird_2)] <- 0
# Fit model
eBird_2_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_2)
summary(eBird_2_mm)
semivario_eBird_2 <- Variogram(eBird_2_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_2, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_2)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_2$MAD, samples_dist_inv_eBird_2, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_2_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_2)
summary(eBird_2_mm_exp)
# gaussian
eBird_2_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_2)
summary(eBird_2_mm_gaus)
# spherical
eBird_2_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_2)
summary(eBird_2_mm_spher)
# linear
eBird_2_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_2)
summary(eBird_2_mm_lin)
# ratio
eBird_2_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_2)
summary(eBird_2_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_2_mm,
  eBird_2_mm_exp,
  eBird_2_mm_gaus,
  eBird_2_mm_spher,
  eBird_2_mm_lin,
  eBird_2_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_2_mm), residuals(eBird_2_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_2_mm_rat), residuals(eBird_2_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_2_mm))
qqline(residuals(eBird_2_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_2_mm_rat))
# qqline(residuals(eBird_2_mm_rat))
# Semivariogram of normal model
semivario_eBird_2_b <- Variogram(eBird_2_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_2_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_2_best <- Variogram(eBird_2_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_2_best, smooth = TRUE)
# Summary of normal model
summary(eBird_2_mm)
# # Summary of best model
# summary(eBird_2_mm_rat)
# Get coefficients
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_2$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_2_mm$coefficients$fixed[2],
      eBird_2_mm_exp$coefficients$fixed[2],
      eBird_2_mm_gaus$coefficients$fixed[2],
      eBird_2_mm_rat$coefficients$fixed[2],
      eBird_2_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_2_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_2_mm$coefficients$fixed[2] - 1.96*summary(eBird_2_mm)$tTable[2,2],
      eBird_2_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_2_mm_exp)$tTable[2,2],
      eBird_2_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_2_mm_gaus)$tTable[2,2],
      eBird_2_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_2_mm_rat)$tTable[2,2],
      eBird_2_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_2_mm_lin)$tTable[2,2],
      # NA,
      eBird_2_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_2_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_2_mm$coefficients$fixed[2] + 1.96*summary(eBird_2_mm)$tTable[2,2],
      eBird_2_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_2_mm_exp)$tTable[2,2],
      eBird_2_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_2_mm_gaus)$tTable[2,2],
      eBird_2_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_2_mm_rat)$tTable[2,2],
      eBird_2_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_2_mm_lin)$tTable[2,2],
      # NA,
      eBird_2_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_2_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_2_mm)$tTable[2,5],
      summary(eBird_2_mm_exp)$tTable[2,5],
      summary(eBird_2_mm_gaus)$tTable[2,5],
      summary(eBird_2_mm_rat)$tTable[2,5],
      summary(eBird_2_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_2_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_3
# Get spatial information
eBird_3
eBird_3_coords <- cbind(eBird_3$long, eBird_3$lat)
# eBird_3_coords <- unique(eBird_3_coords)
samples_dist_eBird_3 <- as.matrix(dist(eBird_3_coords))
samples_dist_inv_eBird_3 <- 1/samples_dist_eBird_3
is.na(samples_dist_inv_eBird_3) <- sapply(samples_dist_inv_eBird_3, is.infinite)
samples_dist_inv_eBird_3[is.na(samples_dist_inv_eBird_3)] <- 0
# Fit model
eBird_3_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_3)
summary(eBird_3_mm)
semivario_eBird_3 <- Variogram(eBird_3_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_3, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_3)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_3$MAD, samples_dist_inv_eBird_3, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_3_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_3)
summary(eBird_3_mm_exp)
# gaussian
eBird_3_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_3)
summary(eBird_3_mm_gaus)
# spherical
eBird_3_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_3)
summary(eBird_3_mm_spher)
# linear
eBird_3_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_3)
summary(eBird_3_mm_lin)
# ratio
eBird_3_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_3)
summary(eBird_3_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_3_mm,
  eBird_3_mm_exp,
  eBird_3_mm_gaus,
  eBird_3_mm_spher,
  # eBird_3_mm_lin,
  eBird_3_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_3_mm), residuals(eBird_3_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_3_mm_rat), residuals(eBird_3_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_3_mm))
qqline(residuals(eBird_3_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_3_mm_rat))
# qqline(residuals(eBird_3_mm_rat))
# Semivariogram of normal model
semivario_eBird_3_b <- Variogram(eBird_3_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_3_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_3_best <- Variogram(eBird_3_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_3_best, smooth = TRUE)
# Summary of normal model
summary(eBird_3_mm)
# # Summary of best model
# summary(eBird_3_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_3$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_3_mm$coefficients$fixed[2],
      eBird_3_mm_exp$coefficients$fixed[2],
      eBird_3_mm_gaus$coefficients$fixed[2],
      eBird_3_mm_rat$coefficients$fixed[2],
      # eBird_3_mm_lin$coefficients$fixed[2],
      NA,
      eBird_3_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_3_mm$coefficients$fixed[2] - 1.96*summary(eBird_3_mm)$tTable[2,2],
      eBird_3_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_3_mm_exp)$tTable[2,2],
      eBird_3_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_3_mm_gaus)$tTable[2,2],
      eBird_3_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_3_mm_rat)$tTable[2,2],
      # eBird_3_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_3_mm_lin)$tTable[2,2],
      NA,
      eBird_3_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_3_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_3_mm$coefficients$fixed[2] + 1.96*summary(eBird_3_mm)$tTable[2,2],
      eBird_3_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_3_mm_exp)$tTable[2,2],
      eBird_3_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_3_mm_gaus)$tTable[2,2],
      eBird_3_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_3_mm_rat)$tTable[2,2],
      # eBird_3_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_3_mm_lin)$tTable[2,2],
      NA,
      eBird_3_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_3_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_3_mm)$tTable[2,5],
      summary(eBird_3_mm_exp)$tTable[2,5],
      summary(eBird_3_mm_gaus)$tTable[2,5],
      summary(eBird_3_mm_rat)$tTable[2,5],
      # summary(eBird_3_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_3_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_4
# Get spatial information
eBird_4
eBird_4_coords <- cbind(eBird_4$long, eBird_4$lat)
# eBird_4_coords <- unique(eBird_4_coords)
samples_dist_eBird_4 <- as.matrix(dist(eBird_4_coords))
samples_dist_inv_eBird_4 <- 1/samples_dist_eBird_4
is.na(samples_dist_inv_eBird_4) <- sapply(samples_dist_inv_eBird_4, is.infinite)
samples_dist_inv_eBird_4[is.na(samples_dist_inv_eBird_4)] <- 0
# Fit model
eBird_4_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_4)
summary(eBird_4_mm)
semivario_eBird_4 <- Variogram(eBird_4_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_4, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_4)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_4$MAD, samples_dist_inv_eBird_4, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_4_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_4)
summary(eBird_4_mm_exp)
# gaussian
eBird_4_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_4)
summary(eBird_4_mm_gaus)
# spherical
eBird_4_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_4)
summary(eBird_4_mm_spher)
# linear
eBird_4_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_4)
summary(eBird_4_mm_lin)
# ratio
eBird_4_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_4)
summary(eBird_4_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_4_mm,
  eBird_4_mm_exp,
  eBird_4_mm_gaus,
  eBird_4_mm_spher,
  eBird_4_mm_lin,
  eBird_4_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_4_mm), residuals(eBird_4_mm))
abline(h=0,lty=3)
# Residuals of best model
plot(fitted(eBird_4_mm_rat), residuals(eBird_4_mm_rat))
abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_4_mm))
qqline(residuals(eBird_4_mm))
# Q-Q plot of best model
qqnorm(residuals(eBird_4_mm_rat))
qqline(residuals(eBird_4_mm_rat))
# Semivariogram of normal model
semivario_eBird_4_b <- Variogram(eBird_4_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_4_b, smooth = TRUE)
# Semivariogram of best model
semivario_eBird_4_best <- Variogram(eBird_4_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
plot(semivario_eBird_4_best, smooth = TRUE)
# Summary of normal model
summary(eBird_4_mm)
# Summary of best model
summary(eBird_4_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_4$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_4_mm$coefficients$fixed[2],
      eBird_4_mm_exp$coefficients$fixed[2],
      eBird_4_mm_gaus$coefficients$fixed[2],
      eBird_4_mm_rat$coefficients$fixed[2],
      eBird_4_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_4_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_4_mm$coefficients$fixed[2] - 1.96*summary(eBird_4_mm)$tTable[2,2],
      eBird_4_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_4_mm_exp)$tTable[2,2],
      eBird_4_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_4_mm_gaus)$tTable[2,2],
      eBird_4_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_4_mm_rat)$tTable[2,2],
      eBird_4_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_4_mm_lin)$tTable[2,2],
      # NA,
      eBird_4_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_4_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_4_mm$coefficients$fixed[2] + 1.96*summary(eBird_4_mm)$tTable[2,2],
      eBird_4_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_4_mm_exp)$tTable[2,2],
      eBird_4_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_4_mm_gaus)$tTable[2,2],
      eBird_4_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_4_mm_rat)$tTable[2,2],
      eBird_4_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_4_mm_lin)$tTable[2,2],
      # NA,
      eBird_4_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_4_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_4_mm)$tTable[2,5],
      summary(eBird_4_mm_exp)$tTable[2,5],
      summary(eBird_4_mm_gaus)$tTable[2,5],
      summary(eBird_4_mm_rat)$tTable[2,5],
      summary(eBird_4_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_4_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_5
# Get spatial information
eBird_5
eBird_5_coords <- cbind(eBird_5$long, eBird_5$lat)
# eBird_5_coords <- unique(eBird_5_coords)
samples_dist_eBird_5 <- as.matrix(dist(eBird_5_coords))
samples_dist_inv_eBird_5 <- 1/samples_dist_eBird_5
is.na(samples_dist_inv_eBird_5) <- sapply(samples_dist_inv_eBird_5, is.infinite)
samples_dist_inv_eBird_5[is.na(samples_dist_inv_eBird_5)] <- 0
# Fit model
eBird_5_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_5)
summary(eBird_5_mm)
semivario_eBird_5 <- Variogram(eBird_5_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_5, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_5)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_5$MAD, samples_dist_inv_eBird_5, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_5_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_5)
summary(eBird_5_mm_exp)
# gaussian
eBird_5_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_5)
summary(eBird_5_mm_gaus)
# spherical
eBird_5_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_5)
summary(eBird_5_mm_spher)
# linear
eBird_5_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_5)
summary(eBird_5_mm_lin)
# ratio
eBird_5_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_5)
summary(eBird_5_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_5_mm,
  eBird_5_mm_exp,
  eBird_5_mm_gaus,
  eBird_5_mm_spher,
  eBird_5_mm_lin,
  eBird_5_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_5_mm), residuals(eBird_5_mm))
abline(h=0,lty=3)
# Residuals of best model
plot(fitted(eBird_5_mm_exp), residuals(eBird_5_mm_exp))
abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_5_mm))
qqline(residuals(eBird_5_mm))
# Q-Q plot of best model
qqnorm(residuals(eBird_5_mm_exp))
qqline(residuals(eBird_5_mm_exp))
# Semivariogram of normal model
semivario_eBird_5_b <- Variogram(eBird_5_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_5_b, smooth = TRUE)
# Semivariogram of best model
semivario_eBird_5_best <- Variogram(eBird_5_mm_exp, form = ~ 1|ID + long + lat, resType = "normalized")
plot(semivario_eBird_5_best, smooth = TRUE)
# Summary of normal model
summary(eBird_5_mm)
# Summary of best model
summary(eBird_5_mm_exp)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_5$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_5_mm$coefficients$fixed[2],
      eBird_5_mm_exp$coefficients$fixed[2],
      eBird_5_mm_gaus$coefficients$fixed[2],
      eBird_5_mm_rat$coefficients$fixed[2],
      eBird_5_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_5_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_5_mm$coefficients$fixed[2] - 1.96*summary(eBird_5_mm)$tTable[2,2],
      eBird_5_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_5_mm_exp)$tTable[2,2],
      eBird_5_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_5_mm_gaus)$tTable[2,2],
      eBird_5_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_5_mm_rat)$tTable[2,2],
      eBird_5_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_5_mm_lin)$tTable[2,2],
      # NA,
      eBird_5_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_5_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_5_mm$coefficients$fixed[2] + 1.96*summary(eBird_5_mm)$tTable[2,2],
      eBird_5_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_5_mm_exp)$tTable[2,2],
      eBird_5_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_5_mm_gaus)$tTable[2,2],
      eBird_5_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_5_mm_rat)$tTable[2,2],
      eBird_5_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_5_mm_lin)$tTable[2,2],
      # NA,
      eBird_5_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_5_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_5_mm)$tTable[2,5],
      summary(eBird_5_mm_exp)$tTable[2,5],
      summary(eBird_5_mm_gaus)$tTable[2,5],
      summary(eBird_5_mm_rat)$tTable[2,5],
      summary(eBird_5_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_5_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_6
# Get spatial information
eBird_6
eBird_6_coords <- cbind(eBird_6$long, eBird_6$lat)
# eBird_6_coords <- unique(eBird_6_coords)
samples_dist_eBird_6 <- as.matrix(dist(eBird_6_coords))
samples_dist_inv_eBird_6 <- 1/samples_dist_eBird_6
is.na(samples_dist_inv_eBird_6) <- sapply(samples_dist_inv_eBird_6, is.infinite)
samples_dist_inv_eBird_6[is.na(samples_dist_inv_eBird_6)] <- 0
# Fit model
eBird_6_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_6)
summary(eBird_6_mm)
semivario_eBird_6 <- Variogram(eBird_6_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_6, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_6)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_6$MAD, samples_dist_inv_eBird_6, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_6_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_6)
summary(eBird_6_mm_exp)
# gaussian
eBird_6_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_6)
summary(eBird_6_mm_gaus)
# spherical
eBird_6_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_6)
summary(eBird_6_mm_spher)
# linear
eBird_6_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_6)
summary(eBird_6_mm_lin)
# ratio
eBird_6_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_6)
summary(eBird_6_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_6_mm,
  eBird_6_mm_exp,
  eBird_6_mm_gaus,
  eBird_6_mm_spher,
  # eBird_6_mm_lin,
  eBird_6_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_6_mm), residuals(eBird_6_mm))
abline(h=0,lty=3)
# Residuals of best model
plot(fitted(eBird_6_mm_exp), residuals(eBird_6_mm_exp))
abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_6_mm))
qqline(residuals(eBird_6_mm))
# Q-Q plot of best model
qqnorm(residuals(eBird_6_mm_exp))
qqline(residuals(eBird_6_mm_exp))
# Semivariogram of normal model
semivario_eBird_6_b <- Variogram(eBird_6_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_6_b, smooth = TRUE)
# Semivariogram of best model
semivario_eBird_6_best <- Variogram(eBird_6_mm_exp, form = ~ 1|ID + long + lat, resType = "normalized")
plot(semivario_eBird_6_best, smooth = TRUE)
# Summary of normal model
summary(eBird_6_mm)
# Summary of best model
summary(eBird_6_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_6$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_6_mm$coefficients$fixed[2],
      eBird_6_mm_exp$coefficients$fixed[2],
      eBird_6_mm_gaus$coefficients$fixed[2],
      eBird_6_mm_rat$coefficients$fixed[2],
      # eBird_6_mm_lin$coefficients$fixed[2],
      NA,
      eBird_6_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_6_mm$coefficients$fixed[2] - 1.96*summary(eBird_6_mm)$tTable[2,2],
      eBird_6_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_6_mm_exp)$tTable[2,2],
      eBird_6_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_6_mm_gaus)$tTable[2,2],
      eBird_6_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_6_mm_rat)$tTable[2,2],
      # eBird_6_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_6_mm_lin)$tTable[2,2],
      NA,
      eBird_6_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_6_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_6_mm$coefficients$fixed[2] + 1.96*summary(eBird_6_mm)$tTable[2,2],
      eBird_6_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_6_mm_exp)$tTable[2,2],
      eBird_6_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_6_mm_gaus)$tTable[2,2],
      eBird_6_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_6_mm_rat)$tTable[2,2],
      # eBird_6_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_6_mm_lin)$tTable[2,2],
      NA,
      eBird_6_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_6_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_6_mm)$tTable[2,5],
      summary(eBird_6_mm_exp)$tTable[2,5],
      summary(eBird_6_mm_gaus)$tTable[2,5],
      summary(eBird_6_mm_rat)$tTable[2,5],
      # summary(eBird_6_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_6_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_7
# Get spatial information
eBird_7
eBird_7_coords <- cbind(eBird_7$long, eBird_7$lat)
# eBird_7_coords <- unique(eBird_7_coords)
samples_dist_eBird_7 <- as.matrix(dist(eBird_7_coords))
samples_dist_inv_eBird_7 <- 1/samples_dist_eBird_7
is.na(samples_dist_inv_eBird_7) <- sapply(samples_dist_inv_eBird_7, is.infinite)
samples_dist_inv_eBird_7[is.na(samples_dist_inv_eBird_7)] <- 0
# Fit model
eBird_7_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_7)
summary(eBird_7_mm)
semivario_eBird_7 <- Variogram(eBird_7_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_7, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_7)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_7$MAD, samples_dist_inv_eBird_7, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_7_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_7)
summary(eBird_7_mm_exp)
# gaussian
eBird_7_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_7)
summary(eBird_7_mm_gaus)
# spherical
eBird_7_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_7)
summary(eBird_7_mm_spher)
# linear
eBird_7_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_7)
summary(eBird_7_mm_lin)
# ratio
eBird_7_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_7)
summary(eBird_7_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_7_mm,
  eBird_7_mm_exp,
  eBird_7_mm_gaus,
  eBird_7_mm_spher,
  eBird_7_mm_lin,
  eBird_7_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_7_mm), residuals(eBird_7_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_7_mm_rat), residuals(eBird_7_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_7_mm))
qqline(residuals(eBird_7_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_7_mm_rat))
# qqline(residuals(eBird_7_mm_rat))
# Semivariogram of normal model
semivario_eBird_7_b <- Variogram(eBird_7_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_7_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_7_best <- Variogram(eBird_7_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_7_best, smooth = TRUE)
# Summary of normal model
summary(eBird_7_mm)
# # Summary of best model
# summary(eBird_7_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_7$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_7_mm$coefficients$fixed[2],
      eBird_7_mm_exp$coefficients$fixed[2],
      eBird_7_mm_gaus$coefficients$fixed[2],
      eBird_7_mm_rat$coefficients$fixed[2],
      eBird_7_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_7_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_7_mm$coefficients$fixed[2] - 1.96*summary(eBird_7_mm)$tTable[2,2],
      eBird_7_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_7_mm_exp)$tTable[2,2],
      eBird_7_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_7_mm_gaus)$tTable[2,2],
      eBird_7_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_7_mm_rat)$tTable[2,2],
      eBird_7_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_7_mm_lin)$tTable[2,2],
      # NA,
      eBird_7_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_7_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_7_mm$coefficients$fixed[2] + 1.96*summary(eBird_7_mm)$tTable[2,2],
      eBird_7_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_7_mm_exp)$tTable[2,2],
      eBird_7_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_7_mm_gaus)$tTable[2,2],
      eBird_7_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_7_mm_rat)$tTable[2,2],
      eBird_7_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_7_mm_lin)$tTable[2,2],
      # NA,
      eBird_7_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_7_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_7_mm)$tTable[2,5],
      summary(eBird_7_mm_exp)$tTable[2,5],
      summary(eBird_7_mm_gaus)$tTable[2,5],
      summary(eBird_7_mm_rat)$tTable[2,5],
      summary(eBird_7_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_7_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_8
# Get spatial information
eBird_8
eBird_8_coords <- cbind(eBird_8$long, eBird_8$lat)
# eBird_8_coords <- unique(eBird_8_coords)
samples_dist_eBird_8 <- as.matrix(dist(eBird_8_coords))
samples_dist_inv_eBird_8 <- 1/samples_dist_eBird_8
is.na(samples_dist_inv_eBird_8) <- sapply(samples_dist_inv_eBird_8, is.infinite)
samples_dist_inv_eBird_8[is.na(samples_dist_inv_eBird_8)] <- 0
# Fit model
eBird_8_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_8)
summary(eBird_8_mm)
semivario_eBird_8 <- Variogram(eBird_8_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_8, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_8)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_8$MAD, samples_dist_inv_eBird_8, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_8_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_8)
summary(eBird_8_mm_exp)
# gaussian
eBird_8_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_8)
summary(eBird_8_mm_gaus)
# spherical
eBird_8_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_8)
summary(eBird_8_mm_spher)
# linear
eBird_8_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_8)
summary(eBird_8_mm_lin)
# ratio
eBird_8_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_8)
summary(eBird_8_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_8_mm,
  eBird_8_mm_exp,
  eBird_8_mm_gaus,
  eBird_8_mm_spher,
  eBird_8_mm_lin,
  eBird_8_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_8_mm), residuals(eBird_8_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_8_mm_rat), residuals(eBird_8_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_8_mm))
qqline(residuals(eBird_8_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_8_mm_rat))
# qqline(residuals(eBird_8_mm_rat))
# Semivariogram of normal model
semivario_eBird_8_b <- Variogram(eBird_8_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_8_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_8_best <- Variogram(eBird_8_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_8_best, smooth = TRUE)
# Summary of normal model
summary(eBird_8_mm)
# # Summary of best model
# summary(eBird_8_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_8$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_8_mm$coefficients$fixed[2],
      eBird_8_mm_exp$coefficients$fixed[2],
      eBird_8_mm_gaus$coefficients$fixed[2],
      eBird_8_mm_rat$coefficients$fixed[2],
      eBird_8_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_8_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_8_mm$coefficients$fixed[2] - 1.96*summary(eBird_8_mm)$tTable[2,2],
      eBird_8_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_8_mm_exp)$tTable[2,2],
      eBird_8_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_8_mm_gaus)$tTable[2,2],
      eBird_8_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_8_mm_rat)$tTable[2,2],
      eBird_8_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_8_mm_lin)$tTable[2,2],
      # NA,
      eBird_8_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_8_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_8_mm$coefficients$fixed[2] + 1.96*summary(eBird_8_mm)$tTable[2,2],
      eBird_8_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_8_mm_exp)$tTable[2,2],
      eBird_8_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_8_mm_gaus)$tTable[2,2],
      eBird_8_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_8_mm_rat)$tTable[2,2],
      eBird_8_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_8_mm_lin)$tTable[2,2],
      # NA,
      eBird_8_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_8_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_8_mm)$tTable[2,5],
      summary(eBird_8_mm_exp)$tTable[2,5],
      summary(eBird_8_mm_gaus)$tTable[2,5],
      summary(eBird_8_mm_rat)$tTable[2,5],
      summary(eBird_8_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_8_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_9
# Get spatial information
eBird_9
eBird_9_coords <- cbind(eBird_9$long, eBird_9$lat)
# eBird_9_coords <- unique(eBird_9_coords)
samples_dist_eBird_9 <- as.matrix(dist(eBird_9_coords))
samples_dist_inv_eBird_9 <- 1/samples_dist_eBird_9
is.na(samples_dist_inv_eBird_9) <- sapply(samples_dist_inv_eBird_9, is.infinite)
samples_dist_inv_eBird_9[is.na(samples_dist_inv_eBird_9)] <- 0
# Fit model
eBird_9_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_9)
summary(eBird_9_mm)
semivario_eBird_9 <- Variogram(eBird_9_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_9, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_9)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_9$MAD, samples_dist_inv_eBird_9, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_9_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_9)
summary(eBird_9_mm_exp)
# gaussian
eBird_9_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_9)
summary(eBird_9_mm_gaus)
# spherical
eBird_9_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_9)
summary(eBird_9_mm_spher)
# linear
eBird_9_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_9)
summary(eBird_9_mm_lin)
# ratio
eBird_9_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_9)
summary(eBird_9_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_9_mm,
  eBird_9_mm_exp,
  eBird_9_mm_gaus,
  eBird_9_mm_spher,
  eBird_9_mm_lin,
  eBird_9_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_9_mm), residuals(eBird_9_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_9_mm_rat), residuals(eBird_9_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_9_mm))
qqline(residuals(eBird_9_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_9_mm_rat))
# qqline(residuals(eBird_9_mm_rat))
# Semivariogram of normal model
semivario_eBird_9_b <- Variogram(eBird_9_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_9_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_9_best <- Variogram(eBird_9_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_9_best, smooth = TRUE)
# Summary of normal model
summary(eBird_9_mm)
# # Summary of best model
# summary(eBird_9_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_9$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_9_mm$coefficients$fixed[2],
      eBird_9_mm_exp$coefficients$fixed[2],
      eBird_9_mm_gaus$coefficients$fixed[2],
      eBird_9_mm_rat$coefficients$fixed[2],
      eBird_9_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_9_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_9_mm$coefficients$fixed[2] - 1.96*summary(eBird_9_mm)$tTable[2,2],
      eBird_9_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_9_mm_exp)$tTable[2,2],
      eBird_9_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_9_mm_gaus)$tTable[2,2],
      eBird_9_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_9_mm_rat)$tTable[2,2],
      eBird_9_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_9_mm_lin)$tTable[2,2],
      # NA,
      eBird_9_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_9_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_9_mm$coefficients$fixed[2] + 1.96*summary(eBird_9_mm)$tTable[2,2],
      eBird_9_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_9_mm_exp)$tTable[2,2],
      eBird_9_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_9_mm_gaus)$tTable[2,2],
      eBird_9_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_9_mm_rat)$tTable[2,2],
      eBird_9_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_9_mm_lin)$tTable[2,2],
      # NA,
      eBird_9_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_9_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_9_mm)$tTable[2,5],
      summary(eBird_9_mm_exp)$tTable[2,5],
      summary(eBird_9_mm_gaus)$tTable[2,5],
      summary(eBird_9_mm_rat)$tTable[2,5],
      summary(eBird_9_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_9_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_10
# Get spatial information
eBird_10
eBird_10_coords <- cbind(eBird_10$long, eBird_10$lat)
# eBird_10_coords <- unique(eBird_10_coords)
samples_dist_eBird_10 <- as.matrix(dist(eBird_10_coords))
samples_dist_inv_eBird_10 <- 1/samples_dist_eBird_10
is.na(samples_dist_inv_eBird_10) <- sapply(samples_dist_inv_eBird_10, is.infinite)
samples_dist_inv_eBird_10[is.na(samples_dist_inv_eBird_10)] <- 0
# Fit model
eBird_10_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_10)
summary(eBird_10_mm)
semivario_eBird_10 <- Variogram(eBird_10_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_10, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_10)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_10$MAD, samples_dist_inv_eBird_10, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_10_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_10)
summary(eBird_10_mm_exp)
# gaussian
eBird_10_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_10)
summary(eBird_10_mm_gaus)
# spherical
eBird_10_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_10)
summary(eBird_10_mm_spher)
# linear
eBird_10_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_10)
summary(eBird_10_mm_lin)
# ratio
eBird_10_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_10)
summary(eBird_10_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_10_mm,
  eBird_10_mm_exp,
  eBird_10_mm_gaus,
  eBird_10_mm_spher,
  # eBird_10_mm_lin,
  eBird_10_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_10_mm), residuals(eBird_10_mm))
abline(h=0,lty=3)
# Residuals of best model
plot(fitted(eBird_10_mm_gaus), residuals(eBird_10_mm_gaus))
abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_10_mm))
qqline(residuals(eBird_10_mm))
# Q-Q plot of best model
qqnorm(residuals(eBird_10_mm_gaus))
qqline(residuals(eBird_10_mm_gaus))
# Semivariogram of normal model
semivario_eBird_10_b <- Variogram(eBird_10_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_10_b, smooth = TRUE)
# Semivariogram of best model
semivario_eBird_10_best <- Variogram(eBird_10_mm_gaus, form = ~ 1|ID + long + lat, resType = "normalized")
plot(semivario_eBird_10_best, smooth = TRUE)
# Summary of normal model
summary(eBird_10_mm)
# Summary of best model
summary(eBird_10_mm_gaus)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_10$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_10_mm$coefficients$fixed[2],
      eBird_10_mm_exp$coefficients$fixed[2],
      eBird_10_mm_gaus$coefficients$fixed[2],
      eBird_10_mm_rat$coefficients$fixed[2],
      # eBird_10_mm_lin$coefficients$fixed[2],
      NA,
      eBird_10_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_10_mm$coefficients$fixed[2] - 1.96*summary(eBird_10_mm)$tTable[2,2],
      eBird_10_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_10_mm_exp)$tTable[2,2],
      eBird_10_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_10_mm_gaus)$tTable[2,2],
      eBird_10_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_10_mm_rat)$tTable[2,2],
      # eBird_10_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_10_mm_lin)$tTable[2,2],
      NA,
      eBird_10_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_10_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_10_mm$coefficients$fixed[2] + 1.96*summary(eBird_10_mm)$tTable[2,2],
      eBird_10_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_10_mm_exp)$tTable[2,2],
      eBird_10_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_10_mm_gaus)$tTable[2,2],
      eBird_10_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_10_mm_rat)$tTable[2,2],
      # eBird_10_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_10_mm_lin)$tTable[2,2],
      NA,
      eBird_10_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_10_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_10_mm)$tTable[2,5],
      summary(eBird_10_mm_exp)$tTable[2,5],
      summary(eBird_10_mm_gaus)$tTable[2,5],
      summary(eBird_10_mm_rat)$tTable[2,5],
      # summary(eBird_10_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_10_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_11
# Get spatial information
eBird_11
eBird_11_coords <- cbind(eBird_11$long, eBird_11$lat)
# eBird_11_coords <- unique(eBird_11_coords)
samples_dist_eBird_11 <- as.matrix(dist(eBird_11_coords))
samples_dist_inv_eBird_11 <- 1/samples_dist_eBird_11
is.na(samples_dist_inv_eBird_11) <- sapply(samples_dist_inv_eBird_11, is.infinite)
samples_dist_inv_eBird_11[is.na(samples_dist_inv_eBird_11)] <- 0
# Fit model
eBird_11_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_11)
summary(eBird_11_mm)
semivario_eBird_11 <- Variogram(eBird_11_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_11, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_11)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_11$MAD, samples_dist_inv_eBird_11, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_11_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_11)
summary(eBird_11_mm_exp)
# gaussian
eBird_11_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_11)
summary(eBird_11_mm_gaus)
# spherical
eBird_11_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_11)
summary(eBird_11_mm_spher)
# linear
eBird_11_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_11)
summary(eBird_11_mm_lin)
# ratio
eBird_11_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_11)
summary(eBird_11_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_11_mm,
  eBird_11_mm_exp,
  eBird_11_mm_gaus,
  eBird_11_mm_spher,
  eBird_11_mm_lin,
  eBird_11_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_11_mm), residuals(eBird_11_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_11_mm_rat), residuals(eBird_11_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_11_mm))
qqline(residuals(eBird_11_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_11_mm_rat))
# qqline(residuals(eBird_11_mm_rat))
# Semivariogram of normal model
semivario_eBird_11_b <- Variogram(eBird_11_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_11_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_11_best <- Variogram(eBird_11_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_11_best, smooth = TRUE)
# Summary of normal model
summary(eBird_11_mm)
# # Summary of best model
# summary(eBird_11_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_11$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_11_mm$coefficients$fixed[2],
      eBird_11_mm_exp$coefficients$fixed[2],
      eBird_11_mm_gaus$coefficients$fixed[2],
      eBird_11_mm_rat$coefficients$fixed[2],
      eBird_11_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_11_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_11_mm$coefficients$fixed[2] - 1.96*summary(eBird_11_mm)$tTable[2,2],
      eBird_11_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_11_mm_exp)$tTable[2,2],
      eBird_11_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_11_mm_gaus)$tTable[2,2],
      eBird_11_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_11_mm_rat)$tTable[2,2],
      eBird_11_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_11_mm_lin)$tTable[2,2],
      # NA,
      eBird_11_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_11_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_11_mm$coefficients$fixed[2] + 1.96*summary(eBird_11_mm)$tTable[2,2],
      eBird_11_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_11_mm_exp)$tTable[2,2],
      eBird_11_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_11_mm_gaus)$tTable[2,2],
      eBird_11_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_11_mm_rat)$tTable[2,2],
      eBird_11_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_11_mm_lin)$tTable[2,2],
      # NA,
      eBird_11_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_11_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_11_mm)$tTable[2,5],
      summary(eBird_11_mm_exp)$tTable[2,5],
      summary(eBird_11_mm_gaus)$tTable[2,5],
      summary(eBird_11_mm_rat)$tTable[2,5],
      summary(eBird_11_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_11_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_12
# Get spatial information
eBird_12
eBird_12_coords <- cbind(eBird_12$long, eBird_12$lat)
# eBird_12_coords <- unique(eBird_12_coords)
samples_dist_eBird_12 <- as.matrix(dist(eBird_12_coords))
samples_dist_inv_eBird_12 <- 1/samples_dist_eBird_12
is.na(samples_dist_inv_eBird_12) <- sapply(samples_dist_inv_eBird_12, is.infinite)
samples_dist_inv_eBird_12[is.na(samples_dist_inv_eBird_12)] <- 0
# Fit model
eBird_12_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_12)
summary(eBird_12_mm)
semivario_eBird_12 <- Variogram(eBird_12_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_12, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_12)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_12$MAD, samples_dist_inv_eBird_12, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_12_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_12)
summary(eBird_12_mm_exp)
# gaussian
eBird_12_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_12)
summary(eBird_12_mm_gaus)
# spherical
eBird_12_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_12)
summary(eBird_12_mm_spher)
# linear
eBird_12_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_12)
summary(eBird_12_mm_lin)
# ratio
eBird_12_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_12)
summary(eBird_12_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_12_mm,
  eBird_12_mm_exp,
  eBird_12_mm_gaus,
  eBird_12_mm_spher,
  eBird_12_mm_lin,
  eBird_12_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_12_mm), residuals(eBird_12_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_12_mm_rat), residuals(eBird_12_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_12_mm))
qqline(residuals(eBird_12_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_12_mm_rat))
# qqline(residuals(eBird_12_mm_rat))
# Semivariogram of normal model
semivario_eBird_12_b <- Variogram(eBird_12_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_12_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_12_best <- Variogram(eBird_12_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_12_best, smooth = TRUE)
# Summary of normal model
summary(eBird_12_mm)
# # Summary of best model
# summary(eBird_12_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_12$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_12_mm$coefficients$fixed[2],
      eBird_12_mm_exp$coefficients$fixed[2],
      eBird_12_mm_gaus$coefficients$fixed[2],
      eBird_12_mm_rat$coefficients$fixed[2],
      eBird_12_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_12_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_12_mm$coefficients$fixed[2] - 1.96*summary(eBird_12_mm)$tTable[2,2],
      eBird_12_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_12_mm_exp)$tTable[2,2],
      eBird_12_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_12_mm_gaus)$tTable[2,2],
      eBird_12_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_12_mm_rat)$tTable[2,2],
      eBird_12_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_12_mm_lin)$tTable[2,2],
      # NA,
      eBird_12_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_12_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_12_mm$coefficients$fixed[2] + 1.96*summary(eBird_12_mm)$tTable[2,2],
      eBird_12_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_12_mm_exp)$tTable[2,2],
      eBird_12_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_12_mm_gaus)$tTable[2,2],
      eBird_12_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_12_mm_rat)$tTable[2,2],
      eBird_12_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_12_mm_lin)$tTable[2,2],
      # NA,
      eBird_12_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_12_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_12_mm)$tTable[2,5],
      summary(eBird_12_mm_exp)$tTable[2,5],
      summary(eBird_12_mm_gaus)$tTable[2,5],
      summary(eBird_12_mm_rat)$tTable[2,5],
      summary(eBird_12_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_12_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_13
# Get spatial information
eBird_13
eBird_13_coords <- cbind(eBird_13$long, eBird_13$lat)
# eBird_13_coords <- unique(eBird_13_coords)
samples_dist_eBird_13 <- as.matrix(dist(eBird_13_coords))
samples_dist_inv_eBird_13 <- 1/samples_dist_eBird_13
is.na(samples_dist_inv_eBird_13) <- sapply(samples_dist_inv_eBird_13, is.infinite)
samples_dist_inv_eBird_13[is.na(samples_dist_inv_eBird_13)] <- 0
# Fit model
eBird_13_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_13)
summary(eBird_13_mm)
semivario_eBird_13 <- Variogram(eBird_13_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_13, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_13)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_13$MAD, samples_dist_inv_eBird_13, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_13_mm_exp <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_13)
summary(eBird_13_mm_exp)
# gaussian
eBird_13_mm_gaus <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_13)
summary(eBird_13_mm_gaus)
# spherical
eBird_13_mm_spher <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_13)
summary(eBird_13_mm_spher)
# linear
eBird_13_mm_lin <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_13)
summary(eBird_13_mm_lin)
# ratio
eBird_13_mm_rat <- lme(fixed = MAD ~ Year, 
                      random = ~ 1|ID,
                      correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_13)
summary(eBird_13_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_13_mm,
  eBird_13_mm_exp,
  eBird_13_mm_gaus,
  eBird_13_mm_spher,
  eBird_13_mm_lin,
  eBird_13_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_13_mm), residuals(eBird_13_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_13_mm_rat), residuals(eBird_13_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_13_mm))
qqline(residuals(eBird_13_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_13_mm_rat))
# qqline(residuals(eBird_13_mm_rat))
# Semivariogram of normal model
semivario_eBird_13_b <- Variogram(eBird_13_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_13_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_13_best <- Variogram(eBird_13_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_13_best, smooth = TRUE)
# Summary of normal model
summary(eBird_13_mm)
# # Summary of best model
# summary(eBird_13_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_13$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_13_mm$coefficients$fixed[2],
      eBird_13_mm_exp$coefficients$fixed[2],
      eBird_13_mm_gaus$coefficients$fixed[2],
      eBird_13_mm_rat$coefficients$fixed[2],
      eBird_13_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_13_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_13_mm$coefficients$fixed[2] - 1.96*summary(eBird_13_mm)$tTable[2,2],
      eBird_13_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_13_mm_exp)$tTable[2,2],
      eBird_13_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_13_mm_gaus)$tTable[2,2],
      eBird_13_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_13_mm_rat)$tTable[2,2],
      eBird_13_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_13_mm_lin)$tTable[2,2],
      # NA,
      eBird_13_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_13_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_13_mm$coefficients$fixed[2] + 1.96*summary(eBird_13_mm)$tTable[2,2],
      eBird_13_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_13_mm_exp)$tTable[2,2],
      eBird_13_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_13_mm_gaus)$tTable[2,2],
      eBird_13_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_13_mm_rat)$tTable[2,2],
      eBird_13_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_13_mm_lin)$tTable[2,2],
      # NA,
      eBird_13_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_13_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_13_mm)$tTable[2,5],
      summary(eBird_13_mm_exp)$tTable[2,5],
      summary(eBird_13_mm_gaus)$tTable[2,5],
      summary(eBird_13_mm_rat)$tTable[2,5],
      summary(eBird_13_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_13_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_14
# Get spatial information
eBird_14
eBird_14_coords <- cbind(eBird_14$long, eBird_14$lat)
# eBird_14_coords <- unique(eBird_14_coords)
samples_dist_eBird_14 <- as.matrix(dist(eBird_14_coords))
samples_dist_inv_eBird_14 <- 1/samples_dist_eBird_14
is.na(samples_dist_inv_eBird_14) <- sapply(samples_dist_inv_eBird_14, is.infinite)
samples_dist_inv_eBird_14[is.na(samples_dist_inv_eBird_14)] <- 0
# Fit model
eBird_14_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_14)
summary(eBird_14_mm)
semivario_eBird_14 <- Variogram(eBird_14_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_14, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_14)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_14$MAD, samples_dist_inv_eBird_14, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_14_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_14)
summary(eBird_14_mm_exp)
# gaussian
eBird_14_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_14)
summary(eBird_14_mm_gaus)
# spherical
eBird_14_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_14)
summary(eBird_14_mm_spher)
# linear
eBird_14_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_14)
summary(eBird_14_mm_lin)
# ratio
eBird_14_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_14)
summary(eBird_14_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_14_mm,
  eBird_14_mm_exp,
  eBird_14_mm_gaus,
  eBird_14_mm_spher,
  eBird_14_mm_lin,
  eBird_14_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_14_mm), residuals(eBird_14_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_14_mm_spher), residuals(eBird_14_mm_spher))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_14_mm))
qqline(residuals(eBird_14_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_14_mm_spher))
# qqline(residuals(eBird_14_mm_spher))
# Semivariogram of normal model
semivario_eBird_14_b <- Variogram(eBird_14_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_14_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_14_best <- Variogram(eBird_14_mm_spher, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_14_best, smooth = TRUE)
# Summary of normal model
summary(eBird_14_mm)
# # Summary of best model
# summary(eBird_14_mm_spher)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_14$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_14_mm$coefficients$fixed[2],
      eBird_14_mm_exp$coefficients$fixed[2],
      eBird_14_mm_gaus$coefficients$fixed[2],
      eBird_14_mm_rat$coefficients$fixed[2],
      eBird_14_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_14_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_14_mm$coefficients$fixed[2] - 1.96*summary(eBird_14_mm)$tTable[2,2],
      eBird_14_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_14_mm_exp)$tTable[2,2],
      eBird_14_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_14_mm_gaus)$tTable[2,2],
      eBird_14_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_14_mm_rat)$tTable[2,2],
      eBird_14_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_14_mm_lin)$tTable[2,2],
      # NA,
      eBird_14_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_14_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_14_mm$coefficients$fixed[2] + 1.96*summary(eBird_14_mm)$tTable[2,2],
      eBird_14_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_14_mm_exp)$tTable[2,2],
      eBird_14_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_14_mm_gaus)$tTable[2,2],
      eBird_14_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_14_mm_rat)$tTable[2,2],
      eBird_14_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_14_mm_lin)$tTable[2,2],
      # NA,
      eBird_14_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_14_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_14_mm)$tTable[2,5],
      summary(eBird_14_mm_exp)$tTable[2,5],
      summary(eBird_14_mm_gaus)$tTable[2,5],
      summary(eBird_14_mm_rat)$tTable[2,5],
      summary(eBird_14_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_14_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_15
# Get spatial information
eBird_15
eBird_15_coords <- cbind(eBird_15$long, eBird_15$lat)
# eBird_15_coords <- unique(eBird_15_coords)
samples_dist_eBird_15 <- as.matrix(dist(eBird_15_coords))
samples_dist_inv_eBird_15 <- 1/samples_dist_eBird_15
is.na(samples_dist_inv_eBird_15) <- sapply(samples_dist_inv_eBird_15, is.infinite)
samples_dist_inv_eBird_15[is.na(samples_dist_inv_eBird_15)] <- 0
# Fit model
eBird_15_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_15)
summary(eBird_15_mm)
semivario_eBird_15 <- Variogram(eBird_15_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_15, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_15)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_15$MAD, samples_dist_inv_eBird_15, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_15_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_15)
summary(eBird_15_mm_exp)
# gaussian
eBird_15_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_15)
summary(eBird_15_mm_gaus)
# spherical
eBird_15_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_15)
summary(eBird_15_mm_spher)
# linear
eBird_15_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_15)
summary(eBird_15_mm_lin)
# ratio
eBird_15_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_15)
summary(eBird_15_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_15_mm,
  eBird_15_mm_exp,
  eBird_15_mm_gaus,
  eBird_15_mm_spher,
  eBird_15_mm_lin,
  eBird_15_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_15_mm), residuals(eBird_15_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_15_mm_rat), residuals(eBird_15_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_15_mm))
qqline(residuals(eBird_15_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_15_mm_rat))
# qqline(residuals(eBird_15_mm_rat))
# Semivariogram of normal model
semivario_eBird_15_b <- Variogram(eBird_15_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_15_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_15_best <- Variogram(eBird_15_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_15_best, smooth = TRUE)
# Summary of normal model
summary(eBird_15_mm)
# # Summary of best model
# summary(eBird_15_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_15$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_15_mm$coefficients$fixed[2],
      eBird_15_mm_exp$coefficients$fixed[2],
      eBird_15_mm_gaus$coefficients$fixed[2],
      eBird_15_mm_rat$coefficients$fixed[2],
      eBird_15_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_15_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_15_mm$coefficients$fixed[2] - 1.96*summary(eBird_15_mm)$tTable[2,2],
      eBird_15_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_15_mm_exp)$tTable[2,2],
      eBird_15_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_15_mm_gaus)$tTable[2,2],
      eBird_15_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_15_mm_rat)$tTable[2,2],
      eBird_15_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_15_mm_lin)$tTable[2,2],
      # NA,
      eBird_15_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_15_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_15_mm$coefficients$fixed[2] + 1.96*summary(eBird_15_mm)$tTable[2,2],
      eBird_15_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_15_mm_exp)$tTable[2,2],
      eBird_15_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_15_mm_gaus)$tTable[2,2],
      eBird_15_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_15_mm_rat)$tTable[2,2],
      eBird_15_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_15_mm_lin)$tTable[2,2],
      # NA,
      eBird_15_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_15_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_15_mm)$tTable[2,5],
      summary(eBird_15_mm_exp)$tTable[2,5],
      summary(eBird_15_mm_gaus)$tTable[2,5],
      summary(eBird_15_mm_rat)$tTable[2,5],
      summary(eBird_15_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_15_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_16
# Get spatial information
eBird_16
eBird_16_coords <- cbind(eBird_16$long, eBird_16$lat)
# eBird_16_coords <- unique(eBird_16_coords)
samples_dist_eBird_16 <- as.matrix(dist(eBird_16_coords))
samples_dist_inv_eBird_16 <- 1/samples_dist_eBird_16
is.na(samples_dist_inv_eBird_16) <- sapply(samples_dist_inv_eBird_16, is.infinite)
samples_dist_inv_eBird_16[is.na(samples_dist_inv_eBird_16)] <- 0
# Fit model
eBird_16_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_16)
summary(eBird_16_mm)
semivario_eBird_16 <- Variogram(eBird_16_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_16, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_16)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_16$MAD, samples_dist_inv_eBird_16, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_16_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_16)
summary(eBird_16_mm_exp)
# gaussian
eBird_16_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_16)
summary(eBird_16_mm_gaus)
# spherical
eBird_16_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_16)
summary(eBird_16_mm_spher)
# linear
eBird_16_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_16)
summary(eBird_16_mm_lin)
# ratio
eBird_16_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_16)
summary(eBird_16_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_16_mm,
  eBird_16_mm_exp,
  eBird_16_mm_gaus,
  eBird_16_mm_spher,
  # eBird_16_mm_lin,
  eBird_16_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_16_mm), residuals(eBird_16_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_16_mm_rat), residuals(eBird_16_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_16_mm))
qqline(residuals(eBird_16_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_16_mm_rat))
# qqline(residuals(eBird_16_mm_rat))
# Semivariogram of normal model
semivario_eBird_16_b <- Variogram(eBird_16_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_16_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_16_best <- Variogram(eBird_16_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_16_best, smooth = TRUE)
# Summary of normal model
summary(eBird_16_mm)
# # Summary of best model
# summary(eBird_16_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_16$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_16_mm$coefficients$fixed[2],
      eBird_16_mm_exp$coefficients$fixed[2],
      eBird_16_mm_gaus$coefficients$fixed[2],
      eBird_16_mm_rat$coefficients$fixed[2],
      # eBird_16_mm_lin$coefficients$fixed[2],
      NA,
      eBird_16_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_16_mm$coefficients$fixed[2] - 1.96*summary(eBird_16_mm)$tTable[2,2],
      eBird_16_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_16_mm_exp)$tTable[2,2],
      eBird_16_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_16_mm_gaus)$tTable[2,2],
      eBird_16_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_16_mm_rat)$tTable[2,2],
      # eBird_16_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_16_mm_lin)$tTable[2,2],
      NA,
      eBird_16_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_16_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_16_mm$coefficients$fixed[2] + 1.96*summary(eBird_16_mm)$tTable[2,2],
      eBird_16_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_16_mm_exp)$tTable[2,2],
      eBird_16_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_16_mm_gaus)$tTable[2,2],
      eBird_16_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_16_mm_rat)$tTable[2,2],
      # eBird_16_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_16_mm_lin)$tTable[2,2],
      NA,
      eBird_16_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_16_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_16_mm)$tTable[2,5],
      summary(eBird_16_mm_exp)$tTable[2,5],
      summary(eBird_16_mm_gaus)$tTable[2,5],
      summary(eBird_16_mm_rat)$tTable[2,5],
      # summary(eBird_16_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_16_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_17
# Get spatial information
eBird_17
eBird_17_coords <- cbind(eBird_17$long, eBird_17$lat)
# eBird_17_coords <- unique(eBird_17_coords)
samples_dist_eBird_17 <- as.matrix(dist(eBird_17_coords))
samples_dist_inv_eBird_17 <- 1/samples_dist_eBird_17
is.na(samples_dist_inv_eBird_17) <- sapply(samples_dist_inv_eBird_17, is.infinite)
samples_dist_inv_eBird_17[is.na(samples_dist_inv_eBird_17)] <- 0
# Fit model
eBird_17_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_17)
summary(eBird_17_mm)
semivario_eBird_17 <- Variogram(eBird_17_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_17, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_17)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_17$MAD, samples_dist_inv_eBird_17, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_17_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_17)
summary(eBird_17_mm_exp)
# gaussian
eBird_17_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_17)
summary(eBird_17_mm_gaus)
# spherical
eBird_17_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_17)
summary(eBird_17_mm_spher)
# linear
eBird_17_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_17)
summary(eBird_17_mm_lin)
# ratio
eBird_17_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_17)
summary(eBird_17_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_17_mm,
  eBird_17_mm_exp,
  eBird_17_mm_gaus,
  eBird_17_mm_spher,
  eBird_17_mm_lin,
  eBird_17_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_17_mm), residuals(eBird_17_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_17_mm_rat), residuals(eBird_17_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_17_mm))
qqline(residuals(eBird_17_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_17_mm_rat))
# qqline(residuals(eBird_17_mm_rat))
# Semivariogram of normal model
semivario_eBird_17_b <- Variogram(eBird_17_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_17_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_17_best <- Variogram(eBird_17_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_17_best, smooth = TRUE)
# Summary of normal model
summary(eBird_17_mm)
# # Summary of best model
# summary(eBird_17_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_17$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_17_mm$coefficients$fixed[2],
      eBird_17_mm_exp$coefficients$fixed[2],
      eBird_17_mm_gaus$coefficients$fixed[2],
      eBird_17_mm_rat$coefficients$fixed[2],
      eBird_17_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_17_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_17_mm$coefficients$fixed[2] - 1.96*summary(eBird_17_mm)$tTable[2,2],
      eBird_17_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_17_mm_exp)$tTable[2,2],
      eBird_17_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_17_mm_gaus)$tTable[2,2],
      eBird_17_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_17_mm_rat)$tTable[2,2],
      eBird_17_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_17_mm_lin)$tTable[2,2],
      # NA,
      eBird_17_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_17_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_17_mm$coefficients$fixed[2] + 1.96*summary(eBird_17_mm)$tTable[2,2],
      eBird_17_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_17_mm_exp)$tTable[2,2],
      eBird_17_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_17_mm_gaus)$tTable[2,2],
      eBird_17_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_17_mm_rat)$tTable[2,2],
      eBird_17_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_17_mm_lin)$tTable[2,2],
      # NA,
      eBird_17_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_17_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_17_mm)$tTable[2,5],
      summary(eBird_17_mm_exp)$tTable[2,5],
      summary(eBird_17_mm_gaus)$tTable[2,5],
      summary(eBird_17_mm_rat)$tTable[2,5],
      summary(eBird_17_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_17_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_18
# Get spatial information
eBird_18
eBird_18_coords <- cbind(eBird_18$long, eBird_18$lat)
# eBird_18_coords <- unique(eBird_18_coords)
samples_dist_eBird_18 <- as.matrix(dist(eBird_18_coords))
samples_dist_inv_eBird_18 <- 1/samples_dist_eBird_18
is.na(samples_dist_inv_eBird_18) <- sapply(samples_dist_inv_eBird_18, is.infinite)
samples_dist_inv_eBird_18[is.na(samples_dist_inv_eBird_18)] <- 0
# Fit model
eBird_18_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_18)
summary(eBird_18_mm)
semivario_eBird_18 <- Variogram(eBird_18_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_18, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_18)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_18$MAD, samples_dist_inv_eBird_18, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_18_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_18)
summary(eBird_18_mm_exp)
# gaussian
eBird_18_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_18)
summary(eBird_18_mm_gaus)
# spherical
eBird_18_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_18)
summary(eBird_18_mm_spher)
# linear
eBird_18_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_18)
summary(eBird_18_mm_lin)
# ratio
eBird_18_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_18)
summary(eBird_18_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_18_mm,
  eBird_18_mm_exp,
  eBird_18_mm_gaus,
  eBird_18_mm_spher,
  eBird_18_mm_lin,
  eBird_18_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_18_mm), residuals(eBird_18_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_18_mm_rat), residuals(eBird_18_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_18_mm))
qqline(residuals(eBird_18_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_18_mm_rat))
# qqline(residuals(eBird_18_mm_rat))
# Semivariogram of normal model
semivario_eBird_18_b <- Variogram(eBird_18_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_18_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_18_best <- Variogram(eBird_18_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_18_best, smooth = TRUE)
# Summary of normal model
summary(eBird_18_mm)
# # Summary of best model
# summary(eBird_18_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_18$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_18_mm$coefficients$fixed[2],
      eBird_18_mm_exp$coefficients$fixed[2],
      eBird_18_mm_gaus$coefficients$fixed[2],
      eBird_18_mm_rat$coefficients$fixed[2],
      eBird_18_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_18_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_18_mm$coefficients$fixed[2] - 1.96*summary(eBird_18_mm)$tTable[2,2],
      eBird_18_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_18_mm_exp)$tTable[2,2],
      eBird_18_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_18_mm_gaus)$tTable[2,2],
      eBird_18_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_18_mm_rat)$tTable[2,2],
      eBird_18_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_18_mm_lin)$tTable[2,2],
      # NA,
      eBird_18_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_18_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_18_mm$coefficients$fixed[2] + 1.96*summary(eBird_18_mm)$tTable[2,2],
      eBird_18_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_18_mm_exp)$tTable[2,2],
      eBird_18_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_18_mm_gaus)$tTable[2,2],
      eBird_18_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_18_mm_rat)$tTable[2,2],
      eBird_18_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_18_mm_lin)$tTable[2,2],
      # NA,
      eBird_18_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_18_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_18_mm)$tTable[2,5],
      summary(eBird_18_mm_exp)$tTable[2,5],
      summary(eBird_18_mm_gaus)$tTable[2,5],
      summary(eBird_18_mm_rat)$tTable[2,5],
      summary(eBird_18_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_18_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_19
# Get spatial information
eBird_19
eBird_19_coords <- cbind(eBird_19$long, eBird_19$lat)
# eBird_19_coords <- unique(eBird_19_coords)
samples_dist_eBird_19 <- as.matrix(dist(eBird_19_coords))
samples_dist_inv_eBird_19 <- 1/samples_dist_eBird_19
is.na(samples_dist_inv_eBird_19) <- sapply(samples_dist_inv_eBird_19, is.infinite)
samples_dist_inv_eBird_19[is.na(samples_dist_inv_eBird_19)] <- 0
# Fit model
eBird_19_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_19)
summary(eBird_19_mm)
semivario_eBird_19 <- Variogram(eBird_19_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_19, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_19)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_19$MAD, samples_dist_inv_eBird_19, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_19_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_19)
summary(eBird_19_mm_exp)
# gaussian
eBird_19_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_19)
summary(eBird_19_mm_gaus)
# spherical
eBird_19_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_19)
summary(eBird_19_mm_spher)
# linear
eBird_19_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_19)
summary(eBird_19_mm_lin)
# ratio
eBird_19_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_19)
summary(eBird_19_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_19_mm,
  eBird_19_mm_exp,
  eBird_19_mm_gaus,
  eBird_19_mm_spher,
  # eBird_19_mm_lin,
  eBird_19_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_19_mm), residuals(eBird_19_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_19_mm_rat), residuals(eBird_19_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_19_mm))
qqline(residuals(eBird_19_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_19_mm_rat))
# qqline(residuals(eBird_19_mm_rat))
# Semivariogram of normal model
semivario_eBird_19_b <- Variogram(eBird_19_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_19_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_19_best <- Variogram(eBird_19_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_19_best, smooth = TRUE)
# Summary of normal model
summary(eBird_19_mm)
# # Summary of best model
# summary(eBird_19_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_19$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_19_mm$coefficients$fixed[2],
      eBird_19_mm_exp$coefficients$fixed[2],
      eBird_19_mm_gaus$coefficients$fixed[2],
      eBird_19_mm_rat$coefficients$fixed[2],
      # eBird_19_mm_lin$coefficients$fixed[2],
      NA,
      eBird_19_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_19_mm$coefficients$fixed[2] - 1.96*summary(eBird_19_mm)$tTable[2,2],
      eBird_19_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_19_mm_exp)$tTable[2,2],
      eBird_19_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_19_mm_gaus)$tTable[2,2],
      eBird_19_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_19_mm_rat)$tTable[2,2],
      # eBird_19_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_19_mm_lin)$tTable[2,2],
      NA,
      eBird_19_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_19_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_19_mm$coefficients$fixed[2] + 1.96*summary(eBird_19_mm)$tTable[2,2],
      eBird_19_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_19_mm_exp)$tTable[2,2],
      eBird_19_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_19_mm_gaus)$tTable[2,2],
      eBird_19_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_19_mm_rat)$tTable[2,2],
      # eBird_19_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_19_mm_lin)$tTable[2,2],
      NA,
      eBird_19_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_19_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_19_mm)$tTable[2,5],
      summary(eBird_19_mm_exp)$tTable[2,5],
      summary(eBird_19_mm_gaus)$tTable[2,5],
      summary(eBird_19_mm_rat)$tTable[2,5],
      # summary(eBird_19_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_19_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_20
# Get spatial information
eBird_20
eBird_20_coords <- cbind(eBird_20$long, eBird_20$lat)
# eBird_20_coords <- unique(eBird_20_coords)
samples_dist_eBird_20 <- as.matrix(dist(eBird_20_coords))
samples_dist_inv_eBird_20 <- 1/samples_dist_eBird_20
is.na(samples_dist_inv_eBird_20) <- sapply(samples_dist_inv_eBird_20, is.infinite)
samples_dist_inv_eBird_20[is.na(samples_dist_inv_eBird_20)] <- 0
# Fit model
eBird_20_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_20)
summary(eBird_20_mm)
semivario_eBird_20 <- Variogram(eBird_20_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_20, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_20)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_20$MAD, samples_dist_inv_eBird_20, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_20_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_20)
summary(eBird_20_mm_exp)
# gaussian
eBird_20_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_20)
summary(eBird_20_mm_gaus)
# spherical
eBird_20_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_20)
summary(eBird_20_mm_spher)
# linear
eBird_20_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_20)
summary(eBird_20_mm_lin)
# ratio
eBird_20_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_20)
summary(eBird_20_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_20_mm,
  eBird_20_mm_exp,
  eBird_20_mm_gaus,
  eBird_20_mm_spher,
  eBird_20_mm_lin,
  eBird_20_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_20_mm), residuals(eBird_20_mm))
abline(h=0,lty=3)
# Residuals of best model
plot(fitted(eBird_20_mm_exp), residuals(eBird_20_mm_exp))
abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_20_mm))
qqline(residuals(eBird_20_mm))
# Q-Q plot of best model
qqnorm(residuals(eBird_20_mm_exp))
qqline(residuals(eBird_20_mm_exp))
# Semivariogram of normal model
semivario_eBird_20_b <- Variogram(eBird_20_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_20_b, smooth = TRUE)
# Semivariogram of best model
semivario_eBird_20_best <- Variogram(eBird_20_mm_exp, form = ~ 1|ID + long + lat, resType = "normalized")
plot(semivario_eBird_20_best, smooth = TRUE)
# Summary of normal model
summary(eBird_20_mm)
# Summary of best model
summary(eBird_20_mm_exp)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_20$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_20_mm$coefficients$fixed[2],
      eBird_20_mm_exp$coefficients$fixed[2],
      eBird_20_mm_gaus$coefficients$fixed[2],
      eBird_20_mm_rat$coefficients$fixed[2],
      eBird_20_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_20_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_20_mm$coefficients$fixed[2] - 1.96*summary(eBird_20_mm)$tTable[2,2],
      eBird_20_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_20_mm_exp)$tTable[2,2],
      eBird_20_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_20_mm_gaus)$tTable[2,2],
      eBird_20_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_20_mm_rat)$tTable[2,2],
      eBird_20_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_20_mm_lin)$tTable[2,2],
      # NA,
      eBird_20_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_20_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_20_mm$coefficients$fixed[2] + 1.96*summary(eBird_20_mm)$tTable[2,2],
      eBird_20_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_20_mm_exp)$tTable[2,2],
      eBird_20_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_20_mm_gaus)$tTable[2,2],
      eBird_20_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_20_mm_rat)$tTable[2,2],
      eBird_20_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_20_mm_lin)$tTable[2,2],
      # NA,
      eBird_20_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_20_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_20_mm)$tTable[2,5],
      summary(eBird_20_mm_exp)$tTable[2,5],
      summary(eBird_20_mm_gaus)$tTable[2,5],
      summary(eBird_20_mm_rat)$tTable[2,5],
      summary(eBird_20_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_20_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_21
# Get spatial information
eBird_21
eBird_21_coords <- cbind(eBird_21$long, eBird_21$lat)
# eBird_21_coords <- unique(eBird_21_coords)
samples_dist_eBird_21 <- as.matrix(dist(eBird_21_coords))
samples_dist_inv_eBird_21 <- 1/samples_dist_eBird_21
is.na(samples_dist_inv_eBird_21) <- sapply(samples_dist_inv_eBird_21, is.infinite)
samples_dist_inv_eBird_21[is.na(samples_dist_inv_eBird_21)] <- 0
# Fit model
eBird_21_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_21)
summary(eBird_21_mm)
semivario_eBird_21 <- Variogram(eBird_21_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_21, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_21)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_21$MAD, samples_dist_inv_eBird_21, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_21_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_21)
summary(eBird_21_mm_exp)
# gaussian
eBird_21_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_21)
summary(eBird_21_mm_gaus)
# spherical
eBird_21_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_21)
summary(eBird_21_mm_spher)
# linear
eBird_21_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_21)
summary(eBird_21_mm_lin)
# ratio
eBird_21_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_21)
summary(eBird_21_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_21_mm,
  eBird_21_mm_exp,
  eBird_21_mm_gaus,
  eBird_21_mm_spher,
  eBird_21_mm_lin,
  eBird_21_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_21_mm), residuals(eBird_21_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_21_mm_rat), residuals(eBird_21_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_21_mm))
qqline(residuals(eBird_21_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_21_mm_rat))
# qqline(residuals(eBird_21_mm_rat))
# Semivariogram of normal model
semivario_eBird_21_b <- Variogram(eBird_21_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_21_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_21_best <- Variogram(eBird_21_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_21_best, smooth = TRUE)
# Summary of normal model
summary(eBird_21_mm)
# # Summary of best model
# summary(eBird_21_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_21$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_21_mm$coefficients$fixed[2],
      eBird_21_mm_exp$coefficients$fixed[2],
      eBird_21_mm_gaus$coefficients$fixed[2],
      eBird_21_mm_rat$coefficients$fixed[2],
      eBird_21_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_21_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_21_mm$coefficients$fixed[2] - 1.96*summary(eBird_21_mm)$tTable[2,2],
      eBird_21_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_21_mm_exp)$tTable[2,2],
      eBird_21_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_21_mm_gaus)$tTable[2,2],
      eBird_21_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_21_mm_rat)$tTable[2,2],
      eBird_21_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_21_mm_lin)$tTable[2,2],
      # NA,
      eBird_21_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_21_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_21_mm$coefficients$fixed[2] + 1.96*summary(eBird_21_mm)$tTable[2,2],
      eBird_21_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_21_mm_exp)$tTable[2,2],
      eBird_21_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_21_mm_gaus)$tTable[2,2],
      eBird_21_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_21_mm_rat)$tTable[2,2],
      eBird_21_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_21_mm_lin)$tTable[2,2],
      # NA,
      eBird_21_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_21_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_21_mm)$tTable[2,5],
      summary(eBird_21_mm_exp)$tTable[2,5],
      summary(eBird_21_mm_gaus)$tTable[2,5],
      summary(eBird_21_mm_rat)$tTable[2,5],
      summary(eBird_21_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_21_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_22
# Get spatial information
eBird_22
eBird_22_coords <- cbind(eBird_22$long, eBird_22$lat)
# eBird_22_coords <- unique(eBird_22_coords)
samples_dist_eBird_22 <- as.matrix(dist(eBird_22_coords))
samples_dist_inv_eBird_22 <- 1/samples_dist_eBird_22
is.na(samples_dist_inv_eBird_22) <- sapply(samples_dist_inv_eBird_22, is.infinite)
samples_dist_inv_eBird_22[is.na(samples_dist_inv_eBird_22)] <- 0
# Fit model
eBird_22_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_22)
summary(eBird_22_mm)
semivario_eBird_22 <- Variogram(eBird_22_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_22, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_22)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_22$MAD, samples_dist_inv_eBird_22, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_22_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_22)
summary(eBird_22_mm_exp)
# gaussian
eBird_22_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_22)
summary(eBird_22_mm_gaus)
# spherical
eBird_22_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_22)
summary(eBird_22_mm_spher)
# linear
eBird_22_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_22)
summary(eBird_22_mm_lin)
# ratio
eBird_22_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_22)
summary(eBird_22_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_22_mm,
  eBird_22_mm_exp,
  eBird_22_mm_gaus,
  eBird_22_mm_spher,
  eBird_22_mm_lin,
  eBird_22_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_22_mm), residuals(eBird_22_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_22_mm_rat), residuals(eBird_22_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_22_mm))
qqline(residuals(eBird_22_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_22_mm_rat))
# qqline(residuals(eBird_22_mm_rat))
# Semivariogram of normal model
semivario_eBird_22_b <- Variogram(eBird_22_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_22_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_22_best <- Variogram(eBird_22_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_22_best, smooth = TRUE)
# Summary of normal model
summary(eBird_22_mm)
# # Summary of best model
# summary(eBird_22_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_22$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_22_mm$coefficients$fixed[2],
      eBird_22_mm_exp$coefficients$fixed[2],
      eBird_22_mm_gaus$coefficients$fixed[2],
      eBird_22_mm_rat$coefficients$fixed[2],
      eBird_22_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_22_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_22_mm$coefficients$fixed[2] - 1.96*summary(eBird_22_mm)$tTable[2,2],
      eBird_22_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_22_mm_exp)$tTable[2,2],
      eBird_22_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_22_mm_gaus)$tTable[2,2],
      eBird_22_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_22_mm_rat)$tTable[2,2],
      eBird_22_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_22_mm_lin)$tTable[2,2],
      # NA,
      eBird_22_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_22_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_22_mm$coefficients$fixed[2] + 1.96*summary(eBird_22_mm)$tTable[2,2],
      eBird_22_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_22_mm_exp)$tTable[2,2],
      eBird_22_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_22_mm_gaus)$tTable[2,2],
      eBird_22_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_22_mm_rat)$tTable[2,2],
      eBird_22_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_22_mm_lin)$tTable[2,2],
      # NA,
      eBird_22_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_22_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_22_mm)$tTable[2,5],
      summary(eBird_22_mm_exp)$tTable[2,5],
      summary(eBird_22_mm_gaus)$tTable[2,5],
      summary(eBird_22_mm_rat)$tTable[2,5],
      summary(eBird_22_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_22_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_23
# Get spatial information
eBird_23
eBird_23_coords <- cbind(eBird_23$long, eBird_23$lat)
# eBird_23_coords <- unique(eBird_23_coords)
samples_dist_eBird_23 <- as.matrix(dist(eBird_23_coords))
samples_dist_inv_eBird_23 <- 1/samples_dist_eBird_23
is.na(samples_dist_inv_eBird_23) <- sapply(samples_dist_inv_eBird_23, is.infinite)
samples_dist_inv_eBird_23[is.na(samples_dist_inv_eBird_23)] <- 0
# Fit model
eBird_23_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_23)
summary(eBird_23_mm)
semivario_eBird_23 <- Variogram(eBird_23_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_23, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_23)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_23$MAD, samples_dist_inv_eBird_23, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_23_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_23)
summary(eBird_23_mm_exp)
# gaussian
eBird_23_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_23)
summary(eBird_23_mm_gaus)
# spherical
eBird_23_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_23)
summary(eBird_23_mm_spher)
# linear
eBird_23_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_23)
summary(eBird_23_mm_lin)
# ratio
eBird_23_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_23)
summary(eBird_23_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_23_mm,
  eBird_23_mm_exp,
  eBird_23_mm_gaus,
  eBird_23_mm_spher,
  eBird_23_mm_lin,
  eBird_23_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_23_mm), residuals(eBird_23_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_23_mm_rat), residuals(eBird_23_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_23_mm))
qqline(residuals(eBird_23_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_23_mm_rat))
# qqline(residuals(eBird_23_mm_rat))
# Semivariogram of normal model
semivario_eBird_23_b <- Variogram(eBird_23_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_23_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_23_best <- Variogram(eBird_23_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_23_best, smooth = TRUE)
# Summary of normal model
summary(eBird_23_mm)
# Summary of best model
summary(eBird_23_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_23$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_23_mm$coefficients$fixed[2],
      eBird_23_mm_exp$coefficients$fixed[2],
      eBird_23_mm_gaus$coefficients$fixed[2],
      eBird_23_mm_rat$coefficients$fixed[2],
      eBird_23_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_23_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_23_mm$coefficients$fixed[2] - 1.96*summary(eBird_23_mm)$tTable[2,2],
      eBird_23_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_23_mm_exp)$tTable[2,2],
      eBird_23_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_23_mm_gaus)$tTable[2,2],
      eBird_23_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_23_mm_rat)$tTable[2,2],
      eBird_23_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_23_mm_lin)$tTable[2,2],
      # NA,
      eBird_23_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_23_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_23_mm$coefficients$fixed[2] + 1.96*summary(eBird_23_mm)$tTable[2,2],
      eBird_23_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_23_mm_exp)$tTable[2,2],
      eBird_23_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_23_mm_gaus)$tTable[2,2],
      eBird_23_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_23_mm_rat)$tTable[2,2],
      eBird_23_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_23_mm_lin)$tTable[2,2],
      # NA,
      eBird_23_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_23_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_23_mm)$tTable[2,5],
      summary(eBird_23_mm_exp)$tTable[2,5],
      summary(eBird_23_mm_gaus)$tTable[2,5],
      summary(eBird_23_mm_rat)$tTable[2,5],
      summary(eBird_23_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_23_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_24
# Get spatial information
eBird_24
eBird_24_coords <- cbind(eBird_24$long, eBird_24$lat)
# eBird_24_coords <- unique(eBird_24_coords)
samples_dist_eBird_24 <- as.matrix(dist(eBird_24_coords))
samples_dist_inv_eBird_24 <- 1/samples_dist_eBird_24
is.na(samples_dist_inv_eBird_24) <- sapply(samples_dist_inv_eBird_24, is.infinite)
samples_dist_inv_eBird_24[is.na(samples_dist_inv_eBird_24)] <- 0
# Fit model
eBird_24_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_24)
summary(eBird_24_mm)
semivario_eBird_24 <- Variogram(eBird_24_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_24, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_24)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_24$MAD, samples_dist_inv_eBird_24, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_24_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_24)
summary(eBird_24_mm_exp)
# gaussian
eBird_24_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_24)
summary(eBird_24_mm_gaus)
# spherical
eBird_24_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_24)
summary(eBird_24_mm_spher)
# linear
eBird_24_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_24)
summary(eBird_24_mm_lin)
# ratio
eBird_24_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_24)
summary(eBird_24_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_24_mm,
  eBird_24_mm_exp,
  eBird_24_mm_gaus,
  eBird_24_mm_spher,
  eBird_24_mm_lin,
  eBird_24_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_24_mm), residuals(eBird_24_mm))
abline(h=0,lty=3)
# Residuals of best model
plot(fitted(eBird_24_mm_exp), residuals(eBird_24_mm_exp))
abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_24_mm))
qqline(residuals(eBird_24_mm))
# Q-Q plot of best model
qqnorm(residuals(eBird_24_mm_exp))
qqline(residuals(eBird_24_mm_exp))
# Semivariogram of normal model
semivario_eBird_24_b <- Variogram(eBird_24_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_24_b, smooth = TRUE)
# Semivariogram of best model
semivario_eBird_24_best <- Variogram(eBird_24_mm_exp, form = ~ 1|ID + long + lat, resType = "normalized")
plot(semivario_eBird_24_best, smooth = TRUE)
# Summary of normal model
summary(eBird_24_mm)
# Summary of best model
summary(eBird_24_mm_exp)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_24$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_24_mm$coefficients$fixed[2],
      eBird_24_mm_exp$coefficients$fixed[2],
      eBird_24_mm_gaus$coefficients$fixed[2],
      eBird_24_mm_rat$coefficients$fixed[2],
      eBird_24_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_24_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_24_mm$coefficients$fixed[2] - 1.96*summary(eBird_24_mm)$tTable[2,2],
      eBird_24_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_24_mm_exp)$tTable[2,2],
      eBird_24_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_24_mm_gaus)$tTable[2,2],
      eBird_24_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_24_mm_rat)$tTable[2,2],
      eBird_24_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_24_mm_lin)$tTable[2,2],
      # NA,
      eBird_24_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_24_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_24_mm$coefficients$fixed[2] + 1.96*summary(eBird_24_mm)$tTable[2,2],
      eBird_24_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_24_mm_exp)$tTable[2,2],
      eBird_24_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_24_mm_gaus)$tTable[2,2],
      eBird_24_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_24_mm_rat)$tTable[2,2],
      eBird_24_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_24_mm_lin)$tTable[2,2],
      # NA,
      eBird_24_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_24_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_24_mm)$tTable[2,5],
      summary(eBird_24_mm_exp)$tTable[2,5],
      summary(eBird_24_mm_gaus)$tTable[2,5],
      summary(eBird_24_mm_rat)$tTable[2,5],
      summary(eBird_24_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_24_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_25
# Get spatial information
eBird_25
eBird_25_coords <- cbind(eBird_25$long, eBird_25$lat)
# eBird_25_coords <- unique(eBird_25_coords)
samples_dist_eBird_25 <- as.matrix(dist(eBird_25_coords))
samples_dist_inv_eBird_25 <- 1/samples_dist_eBird_25
is.na(samples_dist_inv_eBird_25) <- sapply(samples_dist_inv_eBird_25, is.infinite)
samples_dist_inv_eBird_25[is.na(samples_dist_inv_eBird_25)] <- 0
# Fit model
eBird_25_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_25)
summary(eBird_25_mm)
semivario_eBird_25 <- Variogram(eBird_25_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_25, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_25)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_25$MAD, samples_dist_inv_eBird_25, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_25_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_25)
summary(eBird_25_mm_exp)
# gaussian
eBird_25_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_25)
summary(eBird_25_mm_gaus)
# spherical
eBird_25_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_25)
summary(eBird_25_mm_spher)
# linear
eBird_25_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_25)
summary(eBird_25_mm_lin)
# ratio
eBird_25_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_25)
summary(eBird_25_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_25_mm,
  eBird_25_mm_exp,
  eBird_25_mm_gaus,
  eBird_25_mm_spher,
  eBird_25_mm_lin,
  eBird_25_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_25_mm), residuals(eBird_25_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_25_mm_rat), residuals(eBird_25_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_25_mm))
qqline(residuals(eBird_25_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_25_mm_rat))
# qqline(residuals(eBird_25_mm_rat))
# Semivariogram of normal model
semivario_eBird_25_b <- Variogram(eBird_25_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_25_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_25_best <- Variogram(eBird_25_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_25_best, smooth = TRUE)
# Summary of normal model
summary(eBird_25_mm)
# # Summary of best model
# summary(eBird_25_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_25$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_25_mm$coefficients$fixed[2],
      eBird_25_mm_exp$coefficients$fixed[2],
      eBird_25_mm_gaus$coefficients$fixed[2],
      eBird_25_mm_rat$coefficients$fixed[2],
      eBird_25_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_25_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_25_mm$coefficients$fixed[2] - 1.96*summary(eBird_25_mm)$tTable[2,2],
      eBird_25_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_25_mm_exp)$tTable[2,2],
      eBird_25_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_25_mm_gaus)$tTable[2,2],
      eBird_25_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_25_mm_rat)$tTable[2,2],
      eBird_25_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_25_mm_lin)$tTable[2,2],
      # NA,
      eBird_25_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_25_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_25_mm$coefficients$fixed[2] + 1.96*summary(eBird_25_mm)$tTable[2,2],
      eBird_25_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_25_mm_exp)$tTable[2,2],
      eBird_25_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_25_mm_gaus)$tTable[2,2],
      eBird_25_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_25_mm_rat)$tTable[2,2],
      eBird_25_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_25_mm_lin)$tTable[2,2],
      # NA,
      eBird_25_mm_spher$coefficients$fixed[2] + summary(eBird_25_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_25_mm)$tTable[2,5],
      summary(eBird_25_mm_exp)$tTable[2,5],
      summary(eBird_25_mm_gaus)$tTable[2,5],
      summary(eBird_25_mm_rat)$tTable[2,5],
      summary(eBird_25_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_25_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_26
# Get spatial information
eBird_26
eBird_26_coords <- cbind(eBird_26$long, eBird_26$lat)
# eBird_26_coords <- unique(eBird_26_coords)
samples_dist_eBird_26 <- as.matrix(dist(eBird_26_coords))
samples_dist_inv_eBird_26 <- 1/samples_dist_eBird_26
is.na(samples_dist_inv_eBird_26) <- sapply(samples_dist_inv_eBird_26, is.infinite)
samples_dist_inv_eBird_26[is.na(samples_dist_inv_eBird_26)] <- 0
# Fit model
eBird_26_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_26)
summary(eBird_26_mm)
semivario_eBird_26 <- Variogram(eBird_26_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_26, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_26)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_26$MAD, samples_dist_inv_eBird_26, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_26_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_26)
summary(eBird_26_mm_exp)
# gaussian
eBird_26_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_26)
summary(eBird_26_mm_gaus)
# spherical
eBird_26_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_26)
summary(eBird_26_mm_spher)
# linear
eBird_26_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_26)
summary(eBird_26_mm_lin)
# ratio
eBird_26_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_26)
summary(eBird_26_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_26_mm,
  eBird_26_mm_exp,
  eBird_26_mm_gaus,
  eBird_26_mm_spher,
  eBird_26_mm_lin,
  eBird_26_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_26_mm), residuals(eBird_26_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_26_mm_rat), residuals(eBird_26_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_26_mm))
qqline(residuals(eBird_26_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_26_mm_rat))
# qqline(residuals(eBird_26_mm_rat))
# Semivariogram of normal model
semivario_eBird_26_b <- Variogram(eBird_26_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_26_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_26_best <- Variogram(eBird_26_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_26_best, smooth = TRUE)
# Summary of normal model
summary(eBird_26_mm)
# Summary of best model
summary(eBird_26_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_26$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_26_mm$coefficients$fixed[2],
      eBird_26_mm_exp$coefficients$fixed[2],
      eBird_26_mm_gaus$coefficients$fixed[2],
      eBird_26_mm_rat$coefficients$fixed[2],
      eBird_26_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_26_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_26_mm$coefficients$fixed[2] - 1.96*summary(eBird_26_mm)$tTable[2,2],
      eBird_26_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_26_mm_exp)$tTable[2,2],
      eBird_26_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_26_mm_gaus)$tTable[2,2],
      eBird_26_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_26_mm_rat)$tTable[2,2],
      eBird_26_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_26_mm_lin)$tTable[2,2],
      # NA,
      eBird_26_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_26_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_26_mm$coefficients$fixed[2] + 1.96*summary(eBird_26_mm)$tTable[2,2],
      eBird_26_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_26_mm_exp)$tTable[2,2],
      eBird_26_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_26_mm_gaus)$tTable[2,2],
      eBird_26_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_26_mm_rat)$tTable[2,2],
      eBird_26_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_26_mm_lin)$tTable[2,2],
      # NA,
      eBird_26_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_26_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_26_mm)$tTable[2,5],
      summary(eBird_26_mm_exp)$tTable[2,5],
      summary(eBird_26_mm_gaus)$tTable[2,5],
      summary(eBird_26_mm_rat)$tTable[2,5],
      summary(eBird_26_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_26_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_27
# Get spatial information
eBird_27
eBird_27_coords <- cbind(eBird_27$long, eBird_27$lat)
# eBird_27_coords <- unique(eBird_27_coords)
samples_dist_eBird_27 <- as.matrix(dist(eBird_27_coords))
samples_dist_inv_eBird_27 <- 1/samples_dist_eBird_27
is.na(samples_dist_inv_eBird_27) <- sapply(samples_dist_inv_eBird_27, is.infinite)
samples_dist_inv_eBird_27[is.na(samples_dist_inv_eBird_27)] <- 0
# Fit model
eBird_27_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_27)
summary(eBird_27_mm)
semivario_eBird_27 <- Variogram(eBird_27_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_27, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_27)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_27$MAD, samples_dist_inv_eBird_27, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_27_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_27)
summary(eBird_27_mm_exp)
# gaussian
eBird_27_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_27)
summary(eBird_27_mm_gaus)
# spherical
eBird_27_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_27)
summary(eBird_27_mm_spher)
# linear
eBird_27_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_27)
summary(eBird_27_mm_lin)
# ratio
eBird_27_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_27)
summary(eBird_27_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_27_mm,
  eBird_27_mm_exp,
  eBird_27_mm_gaus,
  eBird_27_mm_spher,
  eBird_27_mm_lin,
  eBird_27_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_27_mm), residuals(eBird_27_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_27_mm_spher), residuals(eBird_27_mm_spher))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_27_mm))
qqline(residuals(eBird_27_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_27_mm_spher))
# qqline(residuals(eBird_27_mm_spher))
# Semivariogram of normal model
semivario_eBird_27_b <- Variogram(eBird_27_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_27_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_27_best <- Variogram(eBird_27_mm_spher, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_27_best, smooth = TRUE)
# Summary of normal model
summary(eBird_27_mm)
# # Summary of best model
# summary(eBird_27_mm_spher)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_27$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_27_mm$coefficients$fixed[2],
      eBird_27_mm_exp$coefficients$fixed[2],
      eBird_27_mm_gaus$coefficients$fixed[2],
      eBird_27_mm_rat$coefficients$fixed[2],
      eBird_27_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_27_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_27_mm$coefficients$fixed[2] - 1.96*summary(eBird_27_mm)$tTable[2,2],
      eBird_27_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_27_mm_exp)$tTable[2,2],
      eBird_27_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_27_mm_gaus)$tTable[2,2],
      eBird_27_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_27_mm_rat)$tTable[2,2],
      eBird_27_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_27_mm_lin)$tTable[2,2],
      # NA,
      eBird_27_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_27_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_27_mm$coefficients$fixed[2] + 1.96*summary(eBird_27_mm)$tTable[2,2],
      eBird_27_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_27_mm_exp)$tTable[2,2],
      eBird_27_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_27_mm_gaus)$tTable[2,2],
      eBird_27_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_27_mm_rat)$tTable[2,2],
      eBird_27_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_27_mm_lin)$tTable[2,2],
      # NA,
      eBird_27_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_27_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_27_mm)$tTable[2,5],
      summary(eBird_27_mm_exp)$tTable[2,5],
      summary(eBird_27_mm_gaus)$tTable[2,5],
      summary(eBird_27_mm_rat)$tTable[2,5],
      summary(eBird_27_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_27_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_28
# Get spatial information
eBird_28
eBird_28_coords <- cbind(eBird_28$long, eBird_28$lat)
# eBird_28_coords <- unique(eBird_28_coords)
samples_dist_eBird_28 <- as.matrix(dist(eBird_28_coords))
samples_dist_inv_eBird_28 <- 1/samples_dist_eBird_28
is.na(samples_dist_inv_eBird_28) <- sapply(samples_dist_inv_eBird_28, is.infinite)
samples_dist_inv_eBird_28[is.na(samples_dist_inv_eBird_28)] <- 0
# Fit model
eBird_28_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_28)
summary(eBird_28_mm)
semivario_eBird_28 <- Variogram(eBird_28_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_28, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_28)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_28$MAD, samples_dist_inv_eBird_28, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_28_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_28)
summary(eBird_28_mm_exp)
# gaussian
eBird_28_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_28)
summary(eBird_28_mm_gaus)
# spherical
eBird_28_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_28)
summary(eBird_28_mm_spher)
# linear
eBird_28_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_28)
summary(eBird_28_mm_lin)
# ratio
eBird_28_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_28)
summary(eBird_28_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_28_mm,
  eBird_28_mm_exp,
  eBird_28_mm_gaus,
  eBird_28_mm_spher,
  # eBird_28_mm_lin,
  eBird_28_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_28_mm), residuals(eBird_28_mm))
abline(h=0,lty=3)
# Residuals of best model
plot(fitted(eBird_28_mm_spher), residuals(eBird_28_mm_spher))
abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_28_mm))
qqline(residuals(eBird_28_mm))
# Q-Q plot of best model
qqnorm(residuals(eBird_28_mm_spher))
qqline(residuals(eBird_28_mm_spher))
# Semivariogram of normal model
semivario_eBird_28_b <- Variogram(eBird_28_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_28_b, smooth = TRUE)
# Semivariogram of best model
semivario_eBird_28_best <- Variogram(eBird_28_mm_spher, form = ~ 1|ID + long + lat, resType = "normalized")
plot(semivario_eBird_28_best, smooth = TRUE)
# Summary of normal model
summary(eBird_28_mm)
# Summary of best model
summary(eBird_28_mm_spher)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_28$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_28_mm$coefficients$fixed[2],
      eBird_28_mm_exp$coefficients$fixed[2],
      eBird_28_mm_gaus$coefficients$fixed[2],
      eBird_28_mm_rat$coefficients$fixed[2],
      # eBird_28_mm_lin$coefficients$fixed[2],
      NA,
      eBird_28_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_28_mm$coefficients$fixed[2] - 1.96*summary(eBird_28_mm)$tTable[2,2],
      eBird_28_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_28_mm_exp)$tTable[2,2],
      eBird_28_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_28_mm_gaus)$tTable[2,2],
      eBird_28_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_28_mm_rat)$tTable[2,2],
      # eBird_28_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_28_mm_lin)$tTable[2,2],
      NA,
      eBird_28_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_28_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_28_mm$coefficients$fixed[2] + 1.96*summary(eBird_28_mm)$tTable[2,2],
      eBird_28_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_28_mm_exp)$tTable[2,2],
      eBird_28_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_28_mm_gaus)$tTable[2,2],
      eBird_28_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_28_mm_rat)$tTable[2,2],
      # eBird_28_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_28_mm_lin)$tTable[2,2],
      NA,
      eBird_28_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_28_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_28_mm)$tTable[2,5],
      summary(eBird_28_mm_exp)$tTable[2,5],
      summary(eBird_28_mm_gaus)$tTable[2,5],
      summary(eBird_28_mm_rat)$tTable[2,5],
      # summary(eBird_28_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_28_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_29
# Get spatial information
eBird_29
eBird_29_coords <- cbind(eBird_29$long, eBird_29$lat)
# eBird_29_coords <- unique(eBird_29_coords)
samples_dist_eBird_29 <- as.matrix(dist(eBird_29_coords))
samples_dist_inv_eBird_29 <- 1/samples_dist_eBird_29
is.na(samples_dist_inv_eBird_29) <- sapply(samples_dist_inv_eBird_29, is.infinite)
samples_dist_inv_eBird_29[is.na(samples_dist_inv_eBird_29)] <- 0
# Fit model
eBird_29_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_29)
summary(eBird_29_mm)
semivario_eBird_29 <- Variogram(eBird_29_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_29, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_29)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_29$MAD, samples_dist_inv_eBird_29, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_29_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_29)
summary(eBird_29_mm_exp)
# gaussian
eBird_29_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_29)
summary(eBird_29_mm_gaus)
# spherical
eBird_29_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_29)
summary(eBird_29_mm_spher)
# linear
eBird_29_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_29)
summary(eBird_29_mm_lin)
# ratio
eBird_29_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_29)
summary(eBird_29_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_29_mm,
  eBird_29_mm_exp,
  eBird_29_mm_gaus,
  eBird_29_mm_spher,
  eBird_29_mm_lin,
  eBird_29_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_29_mm), residuals(eBird_29_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_29_mm_rat), residuals(eBird_29_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_29_mm))
qqline(residuals(eBird_29_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_29_mm_rat))
# qqline(residuals(eBird_29_mm_rat))
# Semivariogram of normal model
semivario_eBird_29_b <- Variogram(eBird_29_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_29_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_29_best <- Variogram(eBird_29_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_29_best, smooth = TRUE)
# Summary of normal model
summary(eBird_29_mm)
# # Summary of best model
# summary(eBird_29_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_29$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_29_mm$coefficients$fixed[2],
      eBird_29_mm_exp$coefficients$fixed[2],
      eBird_29_mm_gaus$coefficients$fixed[2],
      eBird_29_mm_rat$coefficients$fixed[2],
      eBird_29_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_29_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_29_mm$coefficients$fixed[2] - 1.96*summary(eBird_29_mm)$tTable[2,2],
      eBird_29_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_29_mm_exp)$tTable[2,2],
      eBird_29_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_29_mm_gaus)$tTable[2,2],
      eBird_29_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_29_mm_rat)$tTable[2,2],
      eBird_29_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_29_mm_lin)$tTable[2,2],
      # NA,
      eBird_29_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_29_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_29_mm$coefficients$fixed[2] + 1.96*summary(eBird_29_mm)$tTable[2,2],
      eBird_29_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_29_mm_exp)$tTable[2,2],
      eBird_29_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_29_mm_gaus)$tTable[2,2],
      eBird_29_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_29_mm_rat)$tTable[2,2],
      eBird_29_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_29_mm_lin)$tTable[2,2],
      # NA,
      eBird_29_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_29_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_29_mm)$tTable[2,5],
      summary(eBird_29_mm_exp)$tTable[2,5],
      summary(eBird_29_mm_gaus)$tTable[2,5],
      summary(eBird_29_mm_rat)$tTable[2,5],
      summary(eBird_29_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_29_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_30
# Get spatial information
eBird_30
eBird_30_coords <- cbind(eBird_30$long, eBird_30$lat)
# eBird_30_coords <- unique(eBird_30_coords)
samples_dist_eBird_30 <- as.matrix(dist(eBird_30_coords))
samples_dist_inv_eBird_30 <- 1/samples_dist_eBird_30
is.na(samples_dist_inv_eBird_30) <- sapply(samples_dist_inv_eBird_30, is.infinite)
samples_dist_inv_eBird_30[is.na(samples_dist_inv_eBird_30)] <- 0
# Fit model
eBird_30_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_30)
summary(eBird_30_mm)
semivario_eBird_30 <- Variogram(eBird_30_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_30, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_30)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_30$MAD, samples_dist_inv_eBird_30, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_30_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_30)
summary(eBird_30_mm_exp)
# gaussian
eBird_30_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_30)
summary(eBird_30_mm_gaus)
# spherical
eBird_30_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_30)
summary(eBird_30_mm_spher)
# linear
eBird_30_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_30)
summary(eBird_30_mm_lin)
# ratio
eBird_30_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_30)
summary(eBird_30_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_30_mm,
  eBird_30_mm_exp,
  eBird_30_mm_gaus,
  eBird_30_mm_spher,
  # eBird_30_mm_lin,
  eBird_30_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_30_mm), residuals(eBird_30_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_30_mm_rat), residuals(eBird_30_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_30_mm))
qqline(residuals(eBird_30_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_30_mm_rat))
# qqline(residuals(eBird_30_mm_rat))
# Semivariogram of normal model
semivario_eBird_30_b <- Variogram(eBird_30_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_30_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_30_best <- Variogram(eBird_30_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_30_best, smooth = TRUE)
# Summary of normal model
summary(eBird_30_mm)
# # Summary of best model
# summary(eBird_30_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_30$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_30_mm$coefficients$fixed[2],
      eBird_30_mm_exp$coefficients$fixed[2],
      eBird_30_mm_gaus$coefficients$fixed[2],
      eBird_30_mm_rat$coefficients$fixed[2],
      # eBird_30_mm_lin$coefficients$fixed[2],
      NA,
      eBird_30_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_30_mm$coefficients$fixed[2] - 1.96*summary(eBird_30_mm)$tTable[2,2],
      eBird_30_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_30_mm_exp)$tTable[2,2],
      eBird_30_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_30_mm_gaus)$tTable[2,2],
      eBird_30_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_30_mm_rat)$tTable[2,2],
      # eBird_30_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_30_mm_lin)$tTable[2,2],
      NA,
      eBird_30_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_30_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_30_mm$coefficients$fixed[2] + 1.96*summary(eBird_30_mm)$tTable[2,2],
      eBird_30_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_30_mm_exp)$tTable[2,2],
      eBird_30_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_30_mm_gaus)$tTable[2,2],
      eBird_30_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_30_mm_rat)$tTable[2,2],
      # eBird_30_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_30_mm_lin)$tTable[2,2],
      NA,
      eBird_30_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_30_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_30_mm)$tTable[2,5],
      summary(eBird_30_mm_exp)$tTable[2,5],
      summary(eBird_30_mm_gaus)$tTable[2,5],
      summary(eBird_30_mm_rat)$tTable[2,5],
      # summary(eBird_30_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_30_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_31
# Get spatial information
eBird_31
eBird_31_coords <- cbind(eBird_31$long, eBird_31$lat)
# eBird_31_coords <- unique(eBird_31_coords)
samples_dist_eBird_31 <- as.matrix(dist(eBird_31_coords))
samples_dist_inv_eBird_31 <- 1/samples_dist_eBird_31
is.na(samples_dist_inv_eBird_31) <- sapply(samples_dist_inv_eBird_31, is.infinite)
samples_dist_inv_eBird_31[is.na(samples_dist_inv_eBird_31)] <- 0
# Fit model
eBird_31_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_31)
summary(eBird_31_mm)
semivario_eBird_31 <- Variogram(eBird_31_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_31, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_31)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_31$MAD, samples_dist_inv_eBird_31, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_31_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_31)
summary(eBird_31_mm_exp)
# gaussian
eBird_31_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_31)
summary(eBird_31_mm_gaus)
# spherical
eBird_31_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_31)
summary(eBird_31_mm_spher)
# linear
eBird_31_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_31)
summary(eBird_31_mm_lin)
# ratio
eBird_31_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_31)
summary(eBird_31_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_31_mm,
  eBird_31_mm_exp,
  eBird_31_mm_gaus,
  eBird_31_mm_spher,
  eBird_31_mm_lin,
  eBird_31_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_31_mm), residuals(eBird_31_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_31_mm_rat), residuals(eBird_31_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_31_mm))
qqline(residuals(eBird_31_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_31_mm_rat))
# qqline(residuals(eBird_31_mm_rat))
# Semivariogram of normal model
semivario_eBird_31_b <- Variogram(eBird_31_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_31_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_31_best <- Variogram(eBird_31_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_31_best, smooth = TRUE)
# Summary of normal model
summary(eBird_31_mm)
# # Summary of best model
# summary(eBird_31_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_31$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_31_mm$coefficients$fixed[2],
      eBird_31_mm_exp$coefficients$fixed[2],
      eBird_31_mm_gaus$coefficients$fixed[2],
      eBird_31_mm_rat$coefficients$fixed[2],
      eBird_31_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_31_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_31_mm$coefficients$fixed[2] - 1.96*summary(eBird_31_mm)$tTable[2,2],
      eBird_31_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_31_mm_exp)$tTable[2,2],
      eBird_31_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_31_mm_gaus)$tTable[2,2],
      eBird_31_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_31_mm_rat)$tTable[2,2],
      eBird_31_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_31_mm_lin)$tTable[2,2],
      # NA,
      eBird_31_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_31_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_31_mm$coefficients$fixed[2] + 1.96*summary(eBird_31_mm)$tTable[2,2],
      eBird_31_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_31_mm_exp)$tTable[2,2],
      eBird_31_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_31_mm_gaus)$tTable[2,2],
      eBird_31_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_31_mm_rat)$tTable[2,2],
      eBird_31_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_31_mm_lin)$tTable[2,2],
      # NA,
      eBird_31_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_31_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_31_mm)$tTable[2,5],
      summary(eBird_31_mm_exp)$tTable[2,5],
      summary(eBird_31_mm_gaus)$tTable[2,5],
      summary(eBird_31_mm_rat)$tTable[2,5],
      summary(eBird_31_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_31_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_32
# Get spatial information
eBird_32
eBird_32_coords <- cbind(eBird_32$long, eBird_32$lat)
# eBird_32_coords <- unique(eBird_32_coords)
samples_dist_eBird_32 <- as.matrix(dist(eBird_32_coords))
samples_dist_inv_eBird_32 <- 1/samples_dist_eBird_32
is.na(samples_dist_inv_eBird_32) <- sapply(samples_dist_inv_eBird_32, is.infinite)
samples_dist_inv_eBird_32[is.na(samples_dist_inv_eBird_32)] <- 0
# Fit model
eBird_32_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_32)
summary(eBird_32_mm)
semivario_eBird_32 <- Variogram(eBird_32_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_32, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_32)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_32$MAD, samples_dist_inv_eBird_32, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_32_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_32)
summary(eBird_32_mm_exp)
# gaussian
eBird_32_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_32)
summary(eBird_32_mm_gaus)
# spherical
eBird_32_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_32)
summary(eBird_32_mm_spher)
# linear
eBird_32_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_32)
summary(eBird_32_mm_lin)
# ratio
eBird_32_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_32)
summary(eBird_32_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_32_mm,
  eBird_32_mm_exp,
  eBird_32_mm_gaus,
  eBird_32_mm_spher,
  eBird_32_mm_lin,
  eBird_32_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_32_mm), residuals(eBird_32_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_32_mm_rat), residuals(eBird_32_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_32_mm))
qqline(residuals(eBird_32_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_32_mm_rat))
# qqline(residuals(eBird_32_mm_rat))
# Semivariogram of normal model
semivario_eBird_32_b <- Variogram(eBird_32_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_32_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_32_best <- Variogram(eBird_32_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_32_best, smooth = TRUE)
# Summary of normal model
summary(eBird_32_mm)
# # Summary of best model
# summary(eBird_32_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_32$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_32_mm$coefficients$fixed[2],
      eBird_32_mm_exp$coefficients$fixed[2],
      eBird_32_mm_gaus$coefficients$fixed[2],
      eBird_32_mm_rat$coefficients$fixed[2],
      eBird_32_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_32_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_32_mm$coefficients$fixed[2] - 1.96*summary(eBird_32_mm)$tTable[2,2],
      eBird_32_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_32_mm_exp)$tTable[2,2],
      eBird_32_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_32_mm_gaus)$tTable[2,2],
      eBird_32_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_32_mm_rat)$tTable[2,2],
      eBird_32_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_32_mm_lin)$tTable[2,2],
      # NA,
      eBird_32_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_32_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_32_mm$coefficients$fixed[2] + 1.96*summary(eBird_32_mm)$tTable[2,2],
      eBird_32_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_32_mm_exp)$tTable[2,2],
      eBird_32_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_32_mm_gaus)$tTable[2,2],
      eBird_32_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_32_mm_rat)$tTable[2,2],
      eBird_32_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_32_mm_lin)$tTable[2,2],
      # NA,
      eBird_32_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_32_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_32_mm)$tTable[2,5],
      summary(eBird_32_mm_exp)$tTable[2,5],
      summary(eBird_32_mm_gaus)$tTable[2,5],
      summary(eBird_32_mm_rat)$tTable[2,5],
      summary(eBird_32_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_32_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_33
# Get spatial information
eBird_33
eBird_33_coords <- cbind(eBird_33$long, eBird_33$lat)
# eBird_33_coords <- unique(eBird_33_coords)
samples_dist_eBird_33 <- as.matrix(dist(eBird_33_coords))
samples_dist_inv_eBird_33 <- 1/samples_dist_eBird_33
is.na(samples_dist_inv_eBird_33) <- sapply(samples_dist_inv_eBird_33, is.infinite)
samples_dist_inv_eBird_33[is.na(samples_dist_inv_eBird_33)] <- 0
# Fit model
eBird_33_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_33)
summary(eBird_33_mm)
semivario_eBird_33 <- Variogram(eBird_33_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_33, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_33)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_33$MAD, samples_dist_inv_eBird_33, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_33_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_33)
summary(eBird_33_mm_exp)
# gaussian
eBird_33_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_33)
summary(eBird_33_mm_gaus)
# spherical
eBird_33_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_33)
summary(eBird_33_mm_spher)
# linear
eBird_33_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_33)
summary(eBird_33_mm_lin)
# ratio
eBird_33_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_33)
summary(eBird_33_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_33_mm,
  eBird_33_mm_exp,
  eBird_33_mm_gaus,
  eBird_33_mm_spher,
  eBird_33_mm_lin,
  eBird_33_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_33_mm), residuals(eBird_33_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_33_mm_rat), residuals(eBird_33_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_33_mm))
qqline(residuals(eBird_33_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_33_mm_rat))
# qqline(residuals(eBird_33_mm_rat))
# Semivariogram of normal model
semivario_eBird_33_b <- Variogram(eBird_33_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_33_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_33_best <- Variogram(eBird_33_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_33_best, smooth = TRUE)
# Summary of normal model
summary(eBird_33_mm)
# # Summary of best model
# summary(eBird_33_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_33$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_33_mm$coefficients$fixed[2],
      eBird_33_mm_exp$coefficients$fixed[2],
      eBird_33_mm_gaus$coefficients$fixed[2],
      eBird_33_mm_rat$coefficients$fixed[2],
      eBird_33_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_33_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_33_mm$coefficients$fixed[2] - 1.96*summary(eBird_33_mm)$tTable[2,2],
      eBird_33_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_33_mm_exp)$tTable[2,2],
      eBird_33_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_33_mm_gaus)$tTable[2,2],
      eBird_33_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_33_mm_rat)$tTable[2,2],
      eBird_33_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_33_mm_lin)$tTable[2,2],
      # NA,
      eBird_33_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_33_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_33_mm$coefficients$fixed[2] + 1.96*summary(eBird_33_mm)$tTable[2,2],
      eBird_33_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_33_mm_exp)$tTable[2,2],
      eBird_33_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_33_mm_gaus)$tTable[2,2],
      eBird_33_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_33_mm_rat)$tTable[2,2],
      eBird_33_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_33_mm_lin)$tTable[2,2],
      # NA,
      eBird_33_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_33_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_33_mm)$tTable[2,5],
      summary(eBird_33_mm_exp)$tTable[2,5],
      summary(eBird_33_mm_gaus)$tTable[2,5],
      summary(eBird_33_mm_rat)$tTable[2,5],
      summary(eBird_33_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_33_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_34
# Get spatial information
eBird_34
eBird_34_coords <- cbind(eBird_34$long, eBird_34$lat)
# eBird_34_coords <- unique(eBird_34_coords)
samples_dist_eBird_34 <- as.matrix(dist(eBird_34_coords))
samples_dist_inv_eBird_34 <- 1/samples_dist_eBird_34
is.na(samples_dist_inv_eBird_34) <- sapply(samples_dist_inv_eBird_34, is.infinite)
samples_dist_inv_eBird_34[is.na(samples_dist_inv_eBird_34)] <- 0
# Fit model
eBird_34_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_34)
summary(eBird_34_mm)
semivario_eBird_34 <- Variogram(eBird_34_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_34, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_34)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_34$MAD, samples_dist_inv_eBird_34, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_34_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_34)
summary(eBird_34_mm_exp)
# gaussian
eBird_34_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_34)
summary(eBird_34_mm_gaus)
# spherical
eBird_34_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_34)
summary(eBird_34_mm_spher)
# linear
eBird_34_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_34)
summary(eBird_34_mm_lin)
# ratio
eBird_34_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_34)
summary(eBird_34_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_34_mm,
  eBird_34_mm_exp,
  eBird_34_mm_gaus,
  eBird_34_mm_spher,
  eBird_34_mm_lin,
  eBird_34_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_34_mm), residuals(eBird_34_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_34_mm_rat), residuals(eBird_34_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_34_mm))
qqline(residuals(eBird_34_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_34_mm_rat))
# qqline(residuals(eBird_34_mm_rat))
# Semivariogram of normal model
semivario_eBird_34_b <- Variogram(eBird_34_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_34_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_34_best <- Variogram(eBird_34_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_34_best, smooth = TRUE)
# Summary of normal model
summary(eBird_34_mm)
# # Summary of best model
# summary(eBird_34_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_34$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_34_mm$coefficients$fixed[2],
      eBird_34_mm_exp$coefficients$fixed[2],
      eBird_34_mm_gaus$coefficients$fixed[2],
      eBird_34_mm_rat$coefficients$fixed[2],
      eBird_34_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_34_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_34_mm$coefficients$fixed[2] - 1.96*summary(eBird_34_mm)$tTable[2,2],
      eBird_34_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_34_mm_exp)$tTable[2,2],
      eBird_34_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_34_mm_gaus)$tTable[2,2],
      eBird_34_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_34_mm_rat)$tTable[2,2],
      eBird_34_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_34_mm_lin)$tTable[2,2],
      # NA,
      eBird_34_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_34_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_34_mm$coefficients$fixed[2] + 1.96*summary(eBird_34_mm)$tTable[2,2],
      eBird_34_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_34_mm_exp)$tTable[2,2],
      eBird_34_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_34_mm_gaus)$tTable[2,2],
      eBird_34_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_34_mm_rat)$tTable[2,2],
      eBird_34_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_34_mm_lin)$tTable[2,2],
      # NA,
      eBird_34_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_34_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_34_mm)$tTable[2,5],
      summary(eBird_34_mm_exp)$tTable[2,5],
      summary(eBird_34_mm_gaus)$tTable[2,5],
      summary(eBird_34_mm_rat)$tTable[2,5],
      summary(eBird_34_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_34_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_35
# Get spatial information
eBird_35
eBird_35_coords <- cbind(eBird_35$long, eBird_35$lat)
# eBird_35_coords <- unique(eBird_35_coords)
samples_dist_eBird_35 <- as.matrix(dist(eBird_35_coords))
samples_dist_inv_eBird_35 <- 1/samples_dist_eBird_35
is.na(samples_dist_inv_eBird_35) <- sapply(samples_dist_inv_eBird_35, is.infinite)
samples_dist_inv_eBird_35[is.na(samples_dist_inv_eBird_35)] <- 0
# Fit model
eBird_35_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_35)
summary(eBird_35_mm)
semivario_eBird_35 <- Variogram(eBird_35_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_35, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_35)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_35$MAD, samples_dist_inv_eBird_35, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_35_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_35)
summary(eBird_35_mm_exp)
# gaussian
eBird_35_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_35)
summary(eBird_35_mm_gaus)
# spherical
eBird_35_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_35)
summary(eBird_35_mm_spher)
# linear
eBird_35_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_35)
summary(eBird_35_mm_lin)
# ratio
eBird_35_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_35)
summary(eBird_35_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_35_mm,
  eBird_35_mm_exp,
  eBird_35_mm_gaus,
  eBird_35_mm_spher,
  eBird_35_mm_lin,
  eBird_35_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_35_mm), residuals(eBird_35_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_35_mm_rat), residuals(eBird_35_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_35_mm))
qqline(residuals(eBird_35_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_35_mm_rat))
# qqline(residuals(eBird_35_mm_rat))
# Semivariogram of normal model
semivario_eBird_35_b <- Variogram(eBird_35_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_35_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_35_best <- Variogram(eBird_35_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_35_best, smooth = TRUE)
# Summary of normal model
summary(eBird_35_mm)
# # Summary of best model
# summary(eBird_35_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_35$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_35_mm$coefficients$fixed[2],
      eBird_35_mm_exp$coefficients$fixed[2],
      eBird_35_mm_gaus$coefficients$fixed[2],
      eBird_35_mm_rat$coefficients$fixed[2],
      eBird_35_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_35_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_35_mm$coefficients$fixed[2] - 1.96*summary(eBird_35_mm)$tTable[2,2],
      eBird_35_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_35_mm_exp)$tTable[2,2],
      eBird_35_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_35_mm_gaus)$tTable[2,2],
      eBird_35_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_35_mm_rat)$tTable[2,2],
      eBird_35_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_35_mm_lin)$tTable[2,2],
      # NA,
      eBird_35_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_35_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_35_mm$coefficients$fixed[2] + 1.96*summary(eBird_35_mm)$tTable[2,2],
      eBird_35_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_35_mm_exp)$tTable[2,2],
      eBird_35_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_35_mm_gaus)$tTable[2,2],
      eBird_35_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_35_mm_rat)$tTable[2,2],
      eBird_35_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_35_mm_lin)$tTable[2,2],
      # NA,
      eBird_35_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_35_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_35_mm)$tTable[2,5],
      summary(eBird_35_mm_exp)$tTable[2,5],
      summary(eBird_35_mm_gaus)$tTable[2,5],
      summary(eBird_35_mm_rat)$tTable[2,5],
      summary(eBird_35_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_35_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_36
# Get spatial information
eBird_36
eBird_36_coords <- cbind(eBird_36$long, eBird_36$lat)
# eBird_36_coords <- unique(eBird_36_coords)
samples_dist_eBird_36 <- as.matrix(dist(eBird_36_coords))
samples_dist_inv_eBird_36 <- 1/samples_dist_eBird_36
is.na(samples_dist_inv_eBird_36) <- sapply(samples_dist_inv_eBird_36, is.infinite)
samples_dist_inv_eBird_36[is.na(samples_dist_inv_eBird_36)] <- 0
# Fit model
eBird_36_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_36)
summary(eBird_36_mm)
semivario_eBird_36 <- Variogram(eBird_36_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_36, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_36)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_36$MAD, samples_dist_inv_eBird_36, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_36_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_36)
summary(eBird_36_mm_exp)
# gaussian
eBird_36_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_36)
summary(eBird_36_mm_gaus)
# spherical
eBird_36_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_36)
summary(eBird_36_mm_spher)
# linear
eBird_36_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_36)
summary(eBird_36_mm_lin)
# ratio
eBird_36_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_36)
summary(eBird_36_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_36_mm,
  eBird_36_mm_exp,
  eBird_36_mm_gaus,
  eBird_36_mm_spher,
  # eBird_36_mm_lin,
  eBird_36_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_36_mm), residuals(eBird_36_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_36_mm_rat), residuals(eBird_36_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_36_mm))
qqline(residuals(eBird_36_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_36_mm_rat))
# qqline(residuals(eBird_36_mm_rat))
# Semivariogram of normal model
semivario_eBird_36_b <- Variogram(eBird_36_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_36_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_36_best <- Variogram(eBird_36_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_36_best, smooth = TRUE)
# Summary of normal model
summary(eBird_36_mm)
# # Summary of best model
# summary(eBird_36_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_36$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_36_mm$coefficients$fixed[2],
      eBird_36_mm_exp$coefficients$fixed[2],
      eBird_36_mm_gaus$coefficients$fixed[2],
      eBird_36_mm_rat$coefficients$fixed[2],
      # eBird_36_mm_lin$coefficients$fixed[2],
      NA,
      eBird_36_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_36_mm$coefficients$fixed[2] - 1.96*summary(eBird_36_mm)$tTable[2,2],
      eBird_36_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_36_mm_exp)$tTable[2,2],
      eBird_36_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_36_mm_gaus)$tTable[2,2],
      eBird_36_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_36_mm_rat)$tTable[2,2],
      # eBird_36_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_36_mm_lin)$tTable[2,2],
      NA,
      eBird_36_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_36_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_36_mm$coefficients$fixed[2] + 1.96*summary(eBird_36_mm)$tTable[2,2],
      eBird_36_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_36_mm_exp)$tTable[2,2],
      eBird_36_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_36_mm_gaus)$tTable[2,2],
      eBird_36_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_36_mm_rat)$tTable[2,2],
      # eBird_36_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_36_mm_lin)$tTable[2,2],
      NA,
      eBird_36_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_36_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_36_mm)$tTable[2,5],
      summary(eBird_36_mm_exp)$tTable[2,5],
      summary(eBird_36_mm_gaus)$tTable[2,5],
      summary(eBird_36_mm_rat)$tTable[2,5],
      # summary(eBird_36_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_36_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_37
# Get spatial information
eBird_37
eBird_37_coords <- cbind(eBird_37$long, eBird_37$lat)
# eBird_37_coords <- unique(eBird_37_coords)
samples_dist_eBird_37 <- as.matrix(dist(eBird_37_coords))
samples_dist_inv_eBird_37 <- 1/samples_dist_eBird_37
is.na(samples_dist_inv_eBird_37) <- sapply(samples_dist_inv_eBird_37, is.infinite)
samples_dist_inv_eBird_37[is.na(samples_dist_inv_eBird_37)] <- 0
# Fit model
eBird_37_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_37)
summary(eBird_37_mm)
semivario_eBird_37 <- Variogram(eBird_37_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_37, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_37)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_37$MAD, samples_dist_inv_eBird_37, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_37_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_37)
summary(eBird_37_mm_exp)
# gaussian
eBird_37_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_37)
summary(eBird_37_mm_gaus)
# spherical
eBird_37_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_37)
summary(eBird_37_mm_spher)
# linear
eBird_37_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_37)
summary(eBird_37_mm_lin)
# ratio
eBird_37_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_37)
summary(eBird_37_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_37_mm,
  eBird_37_mm_exp,
  eBird_37_mm_gaus,
  eBird_37_mm_spher,
  eBird_37_mm_lin,
  eBird_37_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_37_mm), residuals(eBird_37_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_37_mm_rat), residuals(eBird_37_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_37_mm))
qqline(residuals(eBird_37_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_37_mm_rat))
# qqline(residuals(eBird_37_mm_rat))
# Semivariogram of normal model
semivario_eBird_37_b <- Variogram(eBird_37_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_37_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_37_best <- Variogram(eBird_37_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_37_best, smooth = TRUE)
# Summary of normal model
summary(eBird_37_mm)
# # Summary of best model
# summary(eBird_37_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_37$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_37_mm$coefficients$fixed[2],
      eBird_37_mm_exp$coefficients$fixed[2],
      eBird_37_mm_gaus$coefficients$fixed[2],
      eBird_37_mm_rat$coefficients$fixed[2],
      eBird_37_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_37_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_37_mm$coefficients$fixed[2] - 1.96*summary(eBird_37_mm)$tTable[2,2],
      eBird_37_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_37_mm_exp)$tTable[2,2],
      eBird_37_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_37_mm_gaus)$tTable[2,2],
      eBird_37_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_37_mm_rat)$tTable[2,2],
      eBird_37_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_37_mm_lin)$tTable[2,2],
      # NA,
      eBird_37_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_37_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_37_mm$coefficients$fixed[2] + 1.96*summary(eBird_37_mm)$tTable[2,2],
      eBird_37_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_37_mm_exp)$tTable[2,2],
      eBird_37_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_37_mm_gaus)$tTable[2,2],
      eBird_37_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_37_mm_rat)$tTable[2,2],
      eBird_37_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_37_mm_lin)$tTable[2,2],
      # NA,
      eBird_37_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_37_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_37_mm)$tTable[2,5],
      summary(eBird_37_mm_exp)$tTable[2,5],
      summary(eBird_37_mm_gaus)$tTable[2,5],
      summary(eBird_37_mm_rat)$tTable[2,5],
      summary(eBird_37_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_37_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_38
# Get spatial information
eBird_38
eBird_38_coords <- cbind(eBird_38$long, eBird_38$lat)
# eBird_38_coords <- unique(eBird_38_coords)
samples_dist_eBird_38 <- as.matrix(dist(eBird_38_coords))
samples_dist_inv_eBird_38 <- 1/samples_dist_eBird_38
is.na(samples_dist_inv_eBird_38) <- sapply(samples_dist_inv_eBird_38, is.infinite)
samples_dist_inv_eBird_38[is.na(samples_dist_inv_eBird_38)] <- 0
# Fit model
eBird_38_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_38)
summary(eBird_38_mm)
semivario_eBird_38 <- Variogram(eBird_38_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_38, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_38)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_38$MAD, samples_dist_inv_eBird_38, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_38_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_38)
summary(eBird_38_mm_exp)
# gaussian
eBird_38_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_38)
summary(eBird_38_mm_gaus)
# spherical
eBird_38_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_38)
summary(eBird_38_mm_spher)
# linear
eBird_38_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_38)
summary(eBird_38_mm_lin)
# ratio
eBird_38_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_38)
summary(eBird_38_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_38_mm,
  eBird_38_mm_exp,
  eBird_38_mm_gaus,
  eBird_38_mm_spher,
  eBird_38_mm_lin,
  eBird_38_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_38_mm), residuals(eBird_38_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_38_mm_rat), residuals(eBird_38_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_38_mm))
qqline(residuals(eBird_38_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_38_mm_rat))
# qqline(residuals(eBird_38_mm_rat))
# Semivariogram of normal model
semivario_eBird_38_b <- Variogram(eBird_38_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_38_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_38_best <- Variogram(eBird_38_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_38_best, smooth = TRUE)
# Summary of normal model
summary(eBird_38_mm)
# # Summary of best model
# summary(eBird_38_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_38$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_38_mm$coefficients$fixed[2],
      eBird_38_mm_exp$coefficients$fixed[2],
      eBird_38_mm_gaus$coefficients$fixed[2],
      eBird_38_mm_rat$coefficients$fixed[2],
      eBird_38_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_38_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_38_mm$coefficients$fixed[2] - 1.96*summary(eBird_38_mm)$tTable[2,2],
      eBird_38_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_38_mm_exp)$tTable[2,2],
      eBird_38_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_38_mm_gaus)$tTable[2,2],
      eBird_38_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_38_mm_rat)$tTable[2,2],
      eBird_38_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_38_mm_lin)$tTable[2,2],
      # NA,
      eBird_38_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_38_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_38_mm$coefficients$fixed[2] + 1.96*summary(eBird_38_mm)$tTable[2,2],
      eBird_38_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_38_mm_exp)$tTable[2,2],
      eBird_38_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_38_mm_gaus)$tTable[2,2],
      eBird_38_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_38_mm_rat)$tTable[2,2],
      eBird_38_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_38_mm_lin)$tTable[2,2],
      # NA,
      eBird_38_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_38_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_38_mm)$tTable[2,5],
      summary(eBird_38_mm_exp)$tTable[2,5],
      summary(eBird_38_mm_gaus)$tTable[2,5],
      summary(eBird_38_mm_rat)$tTable[2,5],
      summary(eBird_38_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_38_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_39
# Get spatial information
eBird_39
eBird_39_coords <- cbind(eBird_39$long, eBird_39$lat)
# eBird_39_coords <- unique(eBird_39_coords)
samples_dist_eBird_39 <- as.matrix(dist(eBird_39_coords))
samples_dist_inv_eBird_39 <- 1/samples_dist_eBird_39
is.na(samples_dist_inv_eBird_39) <- sapply(samples_dist_inv_eBird_39, is.infinite)
samples_dist_inv_eBird_39[is.na(samples_dist_inv_eBird_39)] <- 0
# Fit model
eBird_39_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_39)
summary(eBird_39_mm)
semivario_eBird_39 <- Variogram(eBird_39_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_39, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_39)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_39$MAD, samples_dist_inv_eBird_39, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_39_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_39)
summary(eBird_39_mm_exp)
# gaussian
eBird_39_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_39)
summary(eBird_39_mm_gaus)
# spherical
eBird_39_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_39)
summary(eBird_39_mm_spher)
# linear
eBird_39_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_39)
summary(eBird_39_mm_lin)
# ratio
eBird_39_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_39)
summary(eBird_39_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_39_mm,
  eBird_39_mm_exp,
  eBird_39_mm_gaus,
  eBird_39_mm_spher,
  # eBird_39_mm_lin,
  eBird_39_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_39_mm), residuals(eBird_39_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_39_mm_rat), residuals(eBird_39_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_39_mm))
qqline(residuals(eBird_39_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_39_mm_rat))
# qqline(residuals(eBird_39_mm_rat))
# Semivariogram of normal model
semivario_eBird_39_b <- Variogram(eBird_39_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_39_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_39_best <- Variogram(eBird_39_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_39_best, smooth = TRUE)
# Summary of normal model
summary(eBird_39_mm)
# # Summary of best model
# summary(eBird_39_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_39$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_39_mm$coefficients$fixed[2],
      eBird_39_mm_exp$coefficients$fixed[2],
      eBird_39_mm_gaus$coefficients$fixed[2],
      eBird_39_mm_rat$coefficients$fixed[2],
      # eBird_39_mm_lin$coefficients$fixed[2],
      NA,
      eBird_39_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_39_mm$coefficients$fixed[2] - 1.96*summary(eBird_39_mm)$tTable[2,2],
      eBird_39_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_39_mm_exp)$tTable[2,2],
      eBird_39_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_39_mm_gaus)$tTable[2,2],
      eBird_39_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_39_mm_rat)$tTable[2,2],
      # eBird_39_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_39_mm_lin)$tTable[2,2],
      NA,
      eBird_39_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_39_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_39_mm$coefficients$fixed[2] + 1.96*summary(eBird_39_mm)$tTable[2,2],
      eBird_39_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_39_mm_exp)$tTable[2,2],
      eBird_39_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_39_mm_gaus)$tTable[2,2],
      eBird_39_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_39_mm_rat)$tTable[2,2],
      # eBird_39_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_39_mm_lin)$tTable[2,2],
      NA,
      eBird_39_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_39_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_39_mm)$tTable[2,5],
      summary(eBird_39_mm_exp)$tTable[2,5],
      summary(eBird_39_mm_gaus)$tTable[2,5],
      summary(eBird_39_mm_rat)$tTable[2,5],
      # summary(eBird_39_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_39_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_40
# Get spatial information
eBird_40
eBird_40_coords <- cbind(eBird_40$long, eBird_40$lat)
# eBird_40_coords <- unique(eBird_40_coords)
samples_dist_eBird_40 <- as.matrix(dist(eBird_40_coords))
samples_dist_inv_eBird_40 <- 1/samples_dist_eBird_40
is.na(samples_dist_inv_eBird_40) <- sapply(samples_dist_inv_eBird_40, is.infinite)
samples_dist_inv_eBird_40[is.na(samples_dist_inv_eBird_40)] <- 0
# Fit model
eBird_40_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_40)
summary(eBird_40_mm)
semivario_eBird_40 <- Variogram(eBird_40_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_40, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_40)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_40$MAD, samples_dist_inv_eBird_40, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_40_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_40)
summary(eBird_40_mm_exp)
# gaussian
eBird_40_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_40)
summary(eBird_40_mm_gaus)
# spherical
eBird_40_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_40)
summary(eBird_40_mm_spher)
# linear
eBird_40_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_40)
summary(eBird_40_mm_lin)
# ratio
eBird_40_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_40)
summary(eBird_40_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_40_mm,
  eBird_40_mm_exp,
  eBird_40_mm_gaus,
  eBird_40_mm_spher,
  eBird_40_mm_lin,
  eBird_40_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_40_mm), residuals(eBird_40_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_40_mm_rat), residuals(eBird_40_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_40_mm))
qqline(residuals(eBird_40_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_40_mm_rat))
# qqline(residuals(eBird_40_mm_rat))
# Semivariogram of normal model
semivario_eBird_40_b <- Variogram(eBird_40_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_40_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_40_best <- Variogram(eBird_40_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_40_best, smooth = TRUE)
# Summary of normal model
summary(eBird_40_mm)
# # Summary of best model
# summary(eBird_40_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_40$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_40_mm$coefficients$fixed[2],
      eBird_40_mm_exp$coefficients$fixed[2],
      eBird_40_mm_gaus$coefficients$fixed[2],
      eBird_40_mm_rat$coefficients$fixed[2],
      eBird_40_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_40_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_40_mm$coefficients$fixed[2] - 1.96*summary(eBird_40_mm)$tTable[2,2],
      eBird_40_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_40_mm_exp)$tTable[2,2],
      eBird_40_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_40_mm_gaus)$tTable[2,2],
      eBird_40_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_40_mm_rat)$tTable[2,2],
      eBird_40_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_40_mm_lin)$tTable[2,2],
      # NA,
      eBird_40_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_40_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_40_mm$coefficients$fixed[2] + 1.96*summary(eBird_40_mm)$tTable[2,2],
      eBird_40_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_40_mm_exp)$tTable[2,2],
      eBird_40_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_40_mm_gaus)$tTable[2,2],
      eBird_40_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_40_mm_rat)$tTable[2,2],
      eBird_40_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_40_mm_lin)$tTable[2,2],
      # NA,
      eBird_40_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_40_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_40_mm)$tTable[2,5],
      summary(eBird_40_mm_exp)$tTable[2,5],
      summary(eBird_40_mm_gaus)$tTable[2,5],
      summary(eBird_40_mm_rat)$tTable[2,5],
      summary(eBird_40_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_40_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_41
# Get spatial information
eBird_41
eBird_41_coords <- cbind(eBird_41$long, eBird_41$lat)
# eBird_41_coords <- unique(eBird_41_coords)
samples_dist_eBird_41 <- as.matrix(dist(eBird_41_coords))
samples_dist_inv_eBird_41 <- 1/samples_dist_eBird_41
is.na(samples_dist_inv_eBird_41) <- sapply(samples_dist_inv_eBird_41, is.infinite)
samples_dist_inv_eBird_41[is.na(samples_dist_inv_eBird_41)] <- 0
# Fit model
eBird_41_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_41)
summary(eBird_41_mm)
semivario_eBird_41 <- Variogram(eBird_41_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_41, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_41)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_41$MAD, samples_dist_inv_eBird_41, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_41_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_41)
summary(eBird_41_mm_exp)
# gaussian
eBird_41_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_41)
summary(eBird_41_mm_gaus)
# spherical
eBird_41_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_41)
summary(eBird_41_mm_spher)
# linear
eBird_41_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_41)
summary(eBird_41_mm_lin)
# ratio
eBird_41_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_41)
summary(eBird_41_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_41_mm,
  eBird_41_mm_exp,
  eBird_41_mm_gaus,
  eBird_41_mm_spher,
  eBird_41_mm_lin,
  eBird_41_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_41_mm), residuals(eBird_41_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_41_mm_rat), residuals(eBird_41_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_41_mm))
qqline(residuals(eBird_41_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_41_mm_rat))
# qqline(residuals(eBird_41_mm_rat))
# Semivariogram of normal model
semivario_eBird_41_b <- Variogram(eBird_41_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_41_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_41_best <- Variogram(eBird_41_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_41_best, smooth = TRUE)
# Summary of normal model
summary(eBird_41_mm)
# # Summary of best model
# summary(eBird_41_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_41$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_41_mm$coefficients$fixed[2],
      eBird_41_mm_exp$coefficients$fixed[2],
      eBird_41_mm_gaus$coefficients$fixed[2],
      eBird_41_mm_rat$coefficients$fixed[2],
      eBird_41_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_41_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_41_mm$coefficients$fixed[2] - 1.96*summary(eBird_41_mm)$tTable[2,2],
      eBird_41_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_41_mm_exp)$tTable[2,2],
      eBird_41_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_41_mm_gaus)$tTable[2,2],
      eBird_41_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_41_mm_rat)$tTable[2,2],
      eBird_41_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_41_mm_lin)$tTable[2,2],
      # NA,
      eBird_41_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_41_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_41_mm$coefficients$fixed[2] + 1.96*summary(eBird_41_mm)$tTable[2,2],
      eBird_41_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_41_mm_exp)$tTable[2,2],
      eBird_41_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_41_mm_gaus)$tTable[2,2],
      eBird_41_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_41_mm_rat)$tTable[2,2],
      eBird_41_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_41_mm_lin)$tTable[2,2],
      # NA,
      eBird_41_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_41_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_41_mm)$tTable[2,5],
      summary(eBird_41_mm_exp)$tTable[2,5],
      summary(eBird_41_mm_gaus)$tTable[2,5],
      summary(eBird_41_mm_rat)$tTable[2,5],
      summary(eBird_41_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_41_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_42
# Get spatial information
eBird_42
eBird_42_coords <- cbind(eBird_42$long, eBird_42$lat)
# eBird_42_coords <- unique(eBird_42_coords)
samples_dist_eBird_42 <- as.matrix(dist(eBird_42_coords))
samples_dist_inv_eBird_42 <- 1/samples_dist_eBird_42
is.na(samples_dist_inv_eBird_42) <- sapply(samples_dist_inv_eBird_42, is.infinite)
samples_dist_inv_eBird_42[is.na(samples_dist_inv_eBird_42)] <- 0
# Fit model
eBird_42_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_42)
summary(eBird_42_mm)
semivario_eBird_42 <- Variogram(eBird_42_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_42, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_42)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_42$MAD, samples_dist_inv_eBird_42, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_42_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_42)
summary(eBird_42_mm_exp)
# gaussian
eBird_42_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_42)
summary(eBird_42_mm_gaus)
# spherical
eBird_42_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_42)
summary(eBird_42_mm_spher)
# linear
eBird_42_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_42)
summary(eBird_42_mm_lin)
# ratio
eBird_42_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_42)
summary(eBird_42_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_42_mm,
  eBird_42_mm_exp,
  eBird_42_mm_gaus,
  eBird_42_mm_spher,
  # eBird_42_mm_lin,
  eBird_42_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_42_mm), residuals(eBird_42_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_42_mm_rat), residuals(eBird_42_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_42_mm))
qqline(residuals(eBird_42_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_42_mm_rat))
# qqline(residuals(eBird_42_mm_rat))
# Semivariogram of normal model
semivario_eBird_42_b <- Variogram(eBird_42_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_42_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_42_best <- Variogram(eBird_42_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_42_best, smooth = TRUE)
# Summary of normal model
summary(eBird_42_mm)
# # Summary of best model
# summary(eBird_42_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_42$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_42_mm$coefficients$fixed[2],
      eBird_42_mm_exp$coefficients$fixed[2],
      eBird_42_mm_gaus$coefficients$fixed[2],
      eBird_42_mm_rat$coefficients$fixed[2],
      # eBird_42_mm_lin$coefficients$fixed[2],
      NA,
      eBird_42_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_42_mm$coefficients$fixed[2] - 1.96*summary(eBird_42_mm)$tTable[2,2],
      eBird_42_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_42_mm_exp)$tTable[2,2],
      eBird_42_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_42_mm_gaus)$tTable[2,2],
      eBird_42_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_42_mm_rat)$tTable[2,2],
      # eBird_42_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_42_mm_lin)$tTable[2,2],
      NA,
      eBird_42_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_42_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_42_mm$coefficients$fixed[2] + 1.96*summary(eBird_42_mm)$tTable[2,2],
      eBird_42_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_42_mm_exp)$tTable[2,2],
      eBird_42_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_42_mm_gaus)$tTable[2,2],
      eBird_42_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_42_mm_rat)$tTable[2,2],
      # eBird_42_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_42_mm_lin)$tTable[2,2],
      NA,
      eBird_42_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_42_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_42_mm)$tTable[2,5],
      summary(eBird_42_mm_exp)$tTable[2,5],
      summary(eBird_42_mm_gaus)$tTable[2,5],
      summary(eBird_42_mm_rat)$tTable[2,5],
      # summary(eBird_42_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_42_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_43
# Get spatial information
eBird_43
eBird_43_coords <- cbind(eBird_43$long, eBird_43$lat)
# eBird_43_coords <- unique(eBird_43_coords)
samples_dist_eBird_43 <- as.matrix(dist(eBird_43_coords))
samples_dist_inv_eBird_43 <- 1/samples_dist_eBird_43
is.na(samples_dist_inv_eBird_43) <- sapply(samples_dist_inv_eBird_43, is.infinite)
samples_dist_inv_eBird_43[is.na(samples_dist_inv_eBird_43)] <- 0
# Fit model
eBird_43_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_43)
summary(eBird_43_mm)
semivario_eBird_43 <- Variogram(eBird_43_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_43, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_43)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_43$MAD, samples_dist_inv_eBird_43, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_43_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_43)
summary(eBird_43_mm_exp)
# gaussian
eBird_43_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_43)
summary(eBird_43_mm_gaus)
# spherical
eBird_43_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_43)
summary(eBird_43_mm_spher)
# linear
eBird_43_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_43)
summary(eBird_43_mm_lin)
# ratio
eBird_43_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_43)
summary(eBird_43_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_43_mm,
  eBird_43_mm_exp,
  eBird_43_mm_gaus,
  eBird_43_mm_spher,
  # eBird_43_mm_lin,
  eBird_43_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_43_mm), residuals(eBird_43_mm))
abline(h=0,lty=3)
# # Residuals of best model
# plot(fitted(eBird_43_mm_rat), residuals(eBird_43_mm_rat))
# abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_43_mm))
qqline(residuals(eBird_43_mm))
# # Q-Q plot of best model
# qqnorm(residuals(eBird_43_mm_rat))
# qqline(residuals(eBird_43_mm_rat))
# Semivariogram of normal model
semivario_eBird_43_b <- Variogram(eBird_43_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_43_b, smooth = TRUE)
# # Semivariogram of best model
# semivario_eBird_43_best <- Variogram(eBird_43_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
# plot(semivario_eBird_43_best, smooth = TRUE)
# Summary of normal model
summary(eBird_43_mm)
# # Summary of best model
# summary(eBird_43_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_43$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_43_mm$coefficients$fixed[2],
      eBird_43_mm_exp$coefficients$fixed[2],
      eBird_43_mm_gaus$coefficients$fixed[2],
      eBird_43_mm_rat$coefficients$fixed[2],
      # eBird_43_mm_lin$coefficients$fixed[2],
      NA,
      eBird_43_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_43_mm$coefficients$fixed[2] - 1.96*summary(eBird_43_mm)$tTable[2,2],
      eBird_43_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_43_mm_exp)$tTable[2,2],
      eBird_43_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_43_mm_gaus)$tTable[2,2],
      eBird_43_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_43_mm_rat)$tTable[2,2],
      # eBird_43_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_43_mm_lin)$tTable[2,2],
      NA,
      eBird_43_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_43_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_43_mm$coefficients$fixed[2] + 1.96*summary(eBird_43_mm)$tTable[2,2],
      eBird_43_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_43_mm_exp)$tTable[2,2],
      eBird_43_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_43_mm_gaus)$tTable[2,2],
      eBird_43_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_43_mm_rat)$tTable[2,2],
      # eBird_43_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_43_mm_lin)$tTable[2,2],
      NA,
      eBird_43_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_43_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_43_mm)$tTable[2,5],
      summary(eBird_43_mm_exp)$tTable[2,5],
      summary(eBird_43_mm_gaus)$tTable[2,5],
      summary(eBird_43_mm_rat)$tTable[2,5],
      # summary(eBird_43_mm_lin)$tTable[2,5],
      NA,
      summary(eBird_43_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# eBird_44
# Get spatial information
eBird_44
eBird_44_coords <- cbind(eBird_44$long, eBird_44$lat)
# eBird_44_coords <- unique(eBird_44_coords)
samples_dist_eBird_44 <- as.matrix(dist(eBird_44_coords))
samples_dist_inv_eBird_44 <- 1/samples_dist_eBird_44
is.na(samples_dist_inv_eBird_44) <- sapply(samples_dist_inv_eBird_44, is.infinite)
samples_dist_inv_eBird_44[is.na(samples_dist_inv_eBird_44)] <- 0
# Fit model
eBird_44_mm <- lme(fixed = MAD ~ Year, random = ~ 1|ID, method = "ML", data = eBird_44)
summary(eBird_44_mm)
semivario_eBird_44 <- Variogram(eBird_44_mm, form = ~ 1|ID, robust = TRUE)
plot(semivario_eBird_44, smooth = TRUE)
# regression <- lme4::lmer(MAD ~ Year + (1|ID), data = eBird_44)
# summary(regression)
# Test for spatial autocorrelation
Moran.I(eBird_44$MAD, samples_dist_inv_eBird_44, alternative = "two.sided")

# Fit models with autocorrelation structures
# exponential
eBird_44_mm_exp <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "exponential", nugget = TRUE), method = "ML", data = eBird_44)
summary(eBird_44_mm_exp)
# gaussian
eBird_44_mm_gaus <- lme(fixed = MAD ~ Year, 
                        random = ~ 1|ID,
                        correlation = corSpatial(form = ~ long + lat , type = "gaussian", nugget = TRUE), method = "ML", data = eBird_44)
summary(eBird_44_mm_gaus)
# spherical
eBird_44_mm_spher <- lme(fixed = MAD ~ Year, 
                         random = ~ 1|ID,
                         correlation = corSpatial(form = ~ long + lat , type = "spherical", nugget = TRUE), method = "ML", data = eBird_44)
summary(eBird_44_mm_spher)
# linear
eBird_44_mm_lin <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "linear", nugget = TRUE), method = "ML", data = eBird_44)
summary(eBird_44_mm_lin)
# ratio
eBird_44_mm_rat <- lme(fixed = MAD ~ Year, 
                       random = ~ 1|ID,
                       correlation = corSpatial(form = ~ long + lat , type = "ratio", nugget = TRUE), method = "ML", data = eBird_44)
summary(eBird_44_mm_rat)
# Model selection, omitting ones with errors
model.sel(
  eBird_44_mm,
  eBird_44_mm_exp,
  eBird_44_mm_gaus,
  eBird_44_mm_spher,
  eBird_44_mm_lin,
  eBird_44_mm_rat
)
# Residuals of normal model
plot(fitted(eBird_44_mm), residuals(eBird_44_mm))
abline(h=0,lty=3)
# Residuals of best model
plot(fitted(eBird_44_mm_rat), residuals(eBird_44_mm_rat))
abline(h=0,lty=3)
# Q-Q plot of normal model
qqnorm(residuals(eBird_44_mm))
qqline(residuals(eBird_44_mm))
# Q-Q plot of best model
qqnorm(residuals(eBird_44_mm_rat))
qqline(residuals(eBird_44_mm_rat))
# Semivariogram of normal model
semivario_eBird_44_b <- Variogram(eBird_44_mm, form = ~ 1|ID, resType = "normalized")
plot(semivario_eBird_44_b, smooth = TRUE)
# Semivariogram of best model
semivario_eBird_44_best <- Variogram(eBird_44_mm_rat, form = ~ 1|ID + long + lat, resType = "normalized")
plot(semivario_eBird_44_best, smooth = TRUE)
# Summary of normal model
summary(eBird_44_mm)
# Summary of best model
summary(eBird_44_mm_rat)
# Get coefficients
coef_df <- rbind(
  coef_df,
  data.frame(
    Common_Name = unique(eBird_44$Common_Name),
    # Model
    Model = c(
      "mm",
      "exp",
      "gaus",
      "rat",
      "lin",
      "spher"
    ),
    # Slope
    Slope = c(
      eBird_44_mm$coefficients$fixed[2],
      eBird_44_mm_exp$coefficients$fixed[2],
      eBird_44_mm_gaus$coefficients$fixed[2],
      eBird_44_mm_rat$coefficients$fixed[2],
      eBird_44_mm_lin$coefficients$fixed[2],
      # NA,
      eBird_44_mm_spher$coefficients$fixed[2]
    ),
    # SE minus
    SE_minus = c(
      eBird_44_mm$coefficients$fixed[2] - 1.96*summary(eBird_44_mm)$tTable[2,2],
      eBird_44_mm_exp$coefficients$fixed[2] - 1.96*summary(eBird_44_mm_exp)$tTable[2,2],
      eBird_44_mm_gaus$coefficients$fixed[2] - 1.96*summary(eBird_44_mm_gaus)$tTable[2,2],
      eBird_44_mm_rat$coefficients$fixed[2] - 1.96*summary(eBird_44_mm_rat)$tTable[2,2],
      eBird_44_mm_lin$coefficients$fixed[2] - 1.96*summary(eBird_44_mm_lin)$tTable[2,2],
      # NA,
      eBird_44_mm_spher$coefficients$fixed[2] - 1.96*summary(eBird_44_mm_spher)$tTable[2,2]
    ),
    # SE plus
    SE_plus = c(
      eBird_44_mm$coefficients$fixed[2] + 1.96*summary(eBird_44_mm)$tTable[2,2],
      eBird_44_mm_exp$coefficients$fixed[2] + 1.96*summary(eBird_44_mm_exp)$tTable[2,2],
      eBird_44_mm_gaus$coefficients$fixed[2] + 1.96*summary(eBird_44_mm_gaus)$tTable[2,2],
      eBird_44_mm_rat$coefficients$fixed[2] + 1.96*summary(eBird_44_mm_rat)$tTable[2,2],
      eBird_44_mm_lin$coefficients$fixed[2] + 1.96*summary(eBird_44_mm_lin)$tTable[2,2],
      # NA,
      eBird_44_mm_spher$coefficients$fixed[2] + 1.96*summary(eBird_44_mm_spher)$tTable[2,2]
    ),
    # p value
    p_value = c(
      summary(eBird_44_mm)$tTable[2,5],
      summary(eBird_44_mm_exp)$tTable[2,5],
      summary(eBird_44_mm_gaus)$tTable[2,5],
      summary(eBird_44_mm_rat)$tTable[2,5],
      summary(eBird_44_mm_lin)$tTable[2,5],
      # NA,
      summary(eBird_44_mm_spher)$tTable[2,5]
    )
  )
)
# Plot
ggplot(coef_df) +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Common name")



# ------------------------------------
# Further analysis
# ------------------------------------
coef_df %>% 
  filter(Model == "lin") %>% 
  count(is.na(Slope))

# Subset birds with names that are with delta-AIC > 2 that isn't the mixed model
# Plot
coef_df %>%
  filter(
    Common_Name %in% c(
      unique(eBird_4$Common_Name),
      unique(eBird_5$Common_Name),
      unique(eBird_24$Common_Name)
    )
  ) %>%
ggplot() +
  geom_point(aes(x = Common_Name, y = Slope, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = SE_minus, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = SE_plus, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Common name")

# Review model selection
# Model selection, omitting ones with errors
# eBird_4
model.sel(
  eBird_4_mm,
  eBird_4_mm_exp,
  eBird_4_mm_gaus,
  eBird_4_mm_spher,
  eBird_4_mm_lin,
  eBird_4_mm_rat
)
summary(eBird_4_mm)
summary(eBird_4_mm_rat)
# 2.48
# eBird_5
model.sel(
  eBird_5_mm,
  eBird_5_mm_exp,
  eBird_5_mm_gaus,
  eBird_5_mm_spher,
  eBird_5_mm_lin,
  eBird_5_mm_rat
)
summary(eBird_5_mm)
summary(eBird_5_mm_spher)
# 2.15
# eBird_24
model.sel(
  eBird_24_mm,
  eBird_24_mm_exp,
  eBird_24_mm_gaus,
  eBird_24_mm_spher,
  eBird_24_mm_lin,
  eBird_24_mm_rat
)
summary(eBird_24_mm)
summary(eBird_24_mm_rat)
# 3.48

write_csv(coef_df, "coef_df_20200204.csv")


# Look closer at coef_df
coef_df <- coef_df %>%
  mutate(p_value_signif = ifelse(p_value < 0.05, TRUE, FALSE),
         CI_plus = 1.96*SE_plus,
         CI_minus = 1.96*SE_minus)

ggplot(coef_df, aes(x = Common_Name, fill = p_value_signif)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# Get mixed model only
coef_df_mm <- coef_df %>%
  filter(Model == "mm") %>%
  select(Common_Name, Slope, p_value,CI_plus, CI_minus) %>%
  rename("Slope_mm" = "Slope", "p_value_mm" = "p_value",
         "CI_plus_mm" = "CI_plus", "CI_minus_mm" = "CI_minus")
# Get a summary of the slopes
summary(coef_df_mm$Slope_mm) ; sd(coef_df_mm$Slope_mm) ; nrow(coef_df_mm)
# Get a summary of the absolute slopes
summary(abs(coef_df_mm$Slope_mm)) ; sd(abs(coef_df_mm$Slope_mm)) ; nrow(coef_df_mm)

# Join to the main data.frame
coef_df <- left_join(coef_df, coef_df_mm)
coef_df <- coef_df %>%
  mutate(
    Slope_difference = Slope - Slope_mm,
    p_value_difference = p_value - p_value_mm,
    CI_plus_difference = CI_plus - CI_plus_mm,
    CI_minus_difference = CI_minus - CI_minus_mm
  )

# Filter out the mixed models
coef_df_no_mm <- coef_df %>%
  filter(Model != "mm")
# Get a summary of the absolute p_value differences
summary(abs(coef_df_no_mm$p_value_difference)) ; sd(abs(coef_df_no_mm$p_value_difference), na.rm = T)
nrow(coef_df_no_mm) - summary(abs(coef_df_no_mm$p_value_difference))[[7]]

# Get a summary of the slope differences
summary(abs(coef_df_no_mm$Slope_difference)) ; sd(abs(coef_df_no_mm$Slope_difference), na.rm = T)
nrow(coef_df_no_mm) - summary(abs(coef_df_no_mm$Slope_difference))[[7]]

# Get a summary of the CI differences
summary(abs(coef_df_no_mm$CI_plus_difference)) ; sd(abs(coef_df_no_mm$CI_plus_difference), na.rm = T)
nrow(coef_df_no_mm) - summary(abs(coef_df_no_mm$CI_plus_difference))[[7]]

# Look for signs switching
coef_df_no_mm <- coef_df_no_mm %>% 
  mutate(CI_plus_switch = ifelse(sign(CI_plus) == sign(CI_plus_mm), FALSE, TRUE),
         CI_minus_switch = ifelse(sign(CI_minus) == sign(CI_minus_mm), FALSE, TRUE))
coef_df_no_mm %>% 
  count(CI_plus_switch)

# Get the differences only for the 6 significant birds
coef_df_no_mm_spatial_sig <- coef_df_no_mm %>%
  filter(
    Common_Name %in% c(
      unique(eBird_4$Common_Name),
      unique(eBird_5$Common_Name),
      unique(eBird_24$Common_Name)
    )
  )
# Get a summary of the p vaue differences
summary(abs(coef_df_no_mm_spatial_sig$p_value_difference)) ; sd(abs(coef_df_no_mm_spatial_sig$p_value_difference), na.rm = T)
nrow(coef_df_no_mm_spatial_sig)
# Signif only
coef_df_no_mm_spatial_sig_2 <- coef_df_no_mm_spatial_sig %>%
  filter(p_value_signif == TRUE)
summary(abs(coef_df_no_mm_spatial_sig_2$p_value_difference)) ; sd(abs(coef_df_no_mm_spatial_sig_2$p_value_difference), na.rm = T)
nrow(coef_df_no_mm_spatial_sig_2)
# Get a summary of the slope differences
summary(abs(coef_df_no_mm_spatial_sig$Slope_difference)) ; sd(abs(coef_df_no_mm_spatial_sig$Slope_difference), na.rm = T)
nrow(coef_df_no_mm_spatial_sig)

# Get a summary of the slopes for mm
summary((coef_df_mm$Slope_mm)) ; sd((coef_df_mm$Slope_mm), na.rm = T)
nrow(coef_df_no_mm_spatial_sig)

coef_df_no_mm <- coef_df_no_mm %>% 
  mutate(CI_range = CI_plus - CI_minus)
cor(coef_df_no_mm$CI_range, abs(coef_df_no_mm$Slope_difference), use = "pairwise.complete.obs")
plot(coef_df_no_mm$CI_range, abs(coef_df_no_mm$Slope_difference))

# ---------------------------------------------------------
# dbMEM
# ---------------------------------------------------------
# sr.value <- function (dfxy, z, xax = 1, yax = 2, method = c("bubble",
#                                                             "greylevel"), zmax = NULL, csize = 1, cpoint = 0, pch = 20,
#                       clegend = 0.75, neig = NULL, cneig = 1, xlim = NULL, ylim = NULL,
#                       grid = TRUE, addaxes = TRUE, cgrid = 0.75, include.origin = TRUE,
#                       origin = c(0, 0), sub = "", csub = 1, possub = "topleft",
#                       pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE)
#   #
#   # Slightly modified version of ade4's s.value() graphical function.
#   # Draws round instead of square bubbles in some plots when argument
#   # "bubble" is called.
#   #
#   # License: GPL-2
#   # Author of the original function s.value: Daniel Chessel
#   # Modification: Francois Gillet, 25 August 2012
#   #
# {
#   dfxy <- data.frame(dfxy)
#   if (length(z) != nrow(dfxy))
#     stop(paste("Non equal row numbers", nrow(dfxy), length(z)))
#   opar <- par(mar = par("mar"))
#   on.exit(par(opar))
#   par(mar = c(0.1, 0.1, 0.1, 0.1))
#   coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax,
#                           xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes,
#                           cgrid = cgrid, include.origin = include.origin, origin = origin,
#                           sub = sub, csub = csub, possub = possub, pixmap = pixmap,
#                           contour = contour, area = area, add.plot = add.plot)
#   if (!is.null(neig))
#   {
#     if (is.null(class(neig))) neig <- NULL
#     if (class(neig) != "neig") neig <- NULL
#     deg <- attr(neig, "degrees")
#     if (length(deg) != length(coo$x)) neig <- NULL
#   }
#   if (!is.null(neig))
#   {
#     fun <- function(x, coo)
#     {
#       segments(coo$x[x[1]], coo$y[x[1]], coo$x[x[2]], coo$y[x[2]],
#                lwd = par("lwd") * cneig)
#     }
#     apply(unclass(neig), 1, fun, coo = coo)
#   }
#   method <- method[1]
#   if (method == "greylevel")
#   {
#     br0 <- pretty(z, 6)
#     nborn <- length(br0)
#     coeff <- diff(par("usr")[1:2])/15
#     numclass <- cut.default(z, br0, include = TRUE, lab = FALSE)
#     valgris <- seq(1, 0, le = (nborn - 1))
#     h <- csize * coeff
#     for (i in 1:(nrow(dfxy)))
#     {
#       symbols(coo$x[i], coo$y[i], circles = h/2,
#               bg = gray(valgris[numclass[i]]),
#               add = TRUE, inch = FALSE)
#     }
#     scatterutil.legend.circle.grey(br0, valgris, h/2, clegend)
#     if (cpoint > 0) points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
#   }
#   else if (method == "bubble")
#   {
#     coeff <- diff(par("usr")[1:2])/15
#     sq <- sqrt(abs(z))
#     if (is.null(zmax)) zmax <- max(abs(z))
#     w1 <- sqrt(zmax)
#     sq <- csize * coeff * sq/w1
#     for (i in 1:(nrow(dfxy)))
#     {
#       if (sign(z[i]) >= 0)
#       {
#         symbols(coo$x[i], coo$y[i], circles = sq[i]/2, bg = "black",
#                 fg = "white", add = TRUE, inch = FALSE)
#       }
#       else
#       {
#         symbols(coo$x[i], coo$y[i], circles = sq[i]/2, bg = "white",
#                 fg = "black", add = TRUE, inch = FALSE)
#       }
#     }
#     br0 <- pretty(z, 4)
#     l0 <- length(br0)
#     br0 <- (br0[1:(l0 - 1)] + br0[2:l0])/2
#     sq0 <- sqrt(abs(br0))
#     sq0 <- csize * coeff * sq0/w1
#     sig0 <- sign(br0)
#     if (clegend > 0) scatterutil.legend.bw.circle(br0, sq0, sig0, clegend)
#     if (cpoint > 0) points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
#   }
#   else if (method == "circlesize") print("not yet implemented")
#   if (!add.plot) box()
#   invisible(match.call())
# }
#
#
#
# scatterutil.legend.bw.circle <- function (br0, sq0, sig0, clegend)
# {
#   br0 <- round(br0, dig = 6)
#   cha <- as.character(br0[1])
#   for (i in (2:(length(br0)))) cha <- paste(cha, br0[i], sep = " ")
#   cex0 <- par("cex") * clegend
#   yh <- max(c(strheight(cha, cex = cex0), sq0))
#   h <- strheight(cha, cex = cex0)
#   y0 <- par("usr")[3] + yh/2 + h/2
#   ltot <- strwidth(cha, cex = cex0) + sum(sq0) + h
#   rect(par("usr")[1] + h/4, y0 - yh/2 - h/4,
#        par("usr")[1] + ltot + h/4, y0 + yh/2 + h/4, col = "white")
#   x0 <- par("usr")[1] + h/2
#   for (i in (1:(length(sq0))))
#   {
#     cha <- br0[i]
#     cha <- paste(" ", cha, sep = "")
#     xh <- strwidth(cha, cex = cex0)
#     text(x0 + xh/2, y0, cha, cex = cex0)
#     z0 <- sq0[i]
#     x0 <- x0 + xh + z0/2
#     if (sig0[i] >= 0)
#       symbols(x0, y0, circles = z0/2, bg = "black", fg = "white",
#               add = TRUE, inch = FALSE)
#     else symbols(x0, y0, circles = z0/2, bg = "white", fg = "black",
#                  add = TRUE, inch = FALSE)
#     x0 <- x0 + z0/2
#   }
#   invisible()
# }
#
#
#
# scatterutil.legend.circle.grey <- function (br0, valgris, h, clegend)
# {
#   if (clegend <= 0) return(invisible())
#   br0 <- round(br0, dig = 6)
#   nborn <- length(br0)
#   cex0 <- par("cex") * clegend
#   x0 <- par("usr")[1] + h
#   x1 <- x0
#   for (i in (2:(nborn)))
#   {
#     x1 <- x1 + h
#     cha <- br0[i]
#     cha <- paste(cha, "]", sep = "")
#     xh <- strwidth(cha, cex = cex0)
#     if (i == (nborn)) break
#     x1 <- x1 + xh + h
#   }
#   yh <- max(strheight(paste(br0), cex = cex0), h)
#   y0 <- par("usr")[3] + yh/2 + h/2
#   rect(par("usr")[1] + h/4, y0 - yh/2 - h/4, x1 - h/4, y0 + yh/2 + h/4,
#        col = "white")
#   x0 <- par("usr")[1] + h
#   for (i in (2:(nborn)))
#   {
#     symbols(x0, y0, circles = h/2, bg = gray(valgris[i - 1]), add = TRUE,
#             inch = FALSE)
#     x0 <- x0 + h
#     cha <- br0[i]
#     if (cha < 1e-05) cha <- round(cha, dig = 3)
#     cha <- paste(cha, "]", sep = "")
#     xh <- strwidth(cha, cex = cex0)
#     if (i == (nborn)) break
#     text(x0 + xh/2, y0, cha, cex = cex0)
#     x0 <- x0 + xh + h
#   }
#   invisible()
# }
#
# load("mite.xy.RDa")
# load("mite.RDa")
# load("mite.env.RDa")
# ## Step 1. Construct the matrix of dbMEM variables
# mite.dbmem.tmp <- dbmem(mite.xy, silent = FALSE)
# mite.dbmem <- as.data.frame(mite.dbmem.tmp)
# # Truncation distance used above:
# (thr <- give.thresh(dist(mite.xy)))
# # Display and count the eigenvalues
# attributes(mite.dbmem.tmp)$values
# length(attributes(mite.dbmem.tmp)$values)
# 
# # Transform data
# mite.h <- decostand(mite, "hel")
# # Detrend data
# mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data = mite.xy))
# 
# ## Step 2. Run the global dbMEM analysis on the detrended
# ## Hellinger-transformed mite data
# (mite.dbmem.rda <- rda(mite.h.det ~., mite.dbmem))
# anova(mite.dbmem.rda)
# ## Step 3. Since the R-square is significant, compute the adjusted
# ## R2 and run a forward selection of the dbmem variables
# (mite.R2a <- RsquareAdj(mite.dbmem.rda)$adj.r.squared)
# (mite.dbmem.fwd <- forward.sel(mite.h.det, as.matrix(mite.dbmem),
#                                adjR2thresh = mite.R2a))
# (nb.sig.dbmem <- nrow(mite.dbmem.fwd)) # Number of signif. dbMEM
# # Identity of the significant dbMEM in increasing order
# (dbmem.sign <- sort(mite.dbmem.fwd[ ,2]))
# # Write the significant dbMEM to a new object
# dbmem.red <- mite.dbmem[ ,c(dbmem.sign)]
# ## Step 4. New dbMEM analysis with 8 significant dbMEM variables
# ## Adjusted R-square after forward selection: R2adj = 0.2418
# (mite.dbmem.rda2 <- rda(mite.h.det ~ ., data = dbmem.red))
# (mite.fwd.R2a <- RsquareAdj(mite.dbmem.rda2)$adj.r.squared)
# anova(mite.dbmem.rda2)
# (axes.test <- anova(mite.dbmem.rda2, by = "axis"))
# # Number of significant axes
# (nb.ax <- length(which(axes.test[ , ncol(axes.test)] <= 0.05)))
# ## Step 5. Plot the significant canonical axes
# mite.rda2.axes <-
#   scores(mite.dbmem.rda2,
#          choices = c(1:nb.ax),
#          display = "lc",
#          scaling = 1)
# par(mfrow = c(1,nb.ax))
# for(i in 1:nb.ax){
#   sr.value(mite.xy, mite.rda2.axes[ ,i],
#            sub = paste("RDA",i),
#            csub = 2)
# }
# 
# # Interpreting the spatial variation: regression of the significant
# # canonical axes on the environmental variables, with Shapiro-Wilk
# # normality tests of residuals
# mite.rda2.axis1.env <- lm(mite.rda2.axes[ ,1] ~ ., data = mite.env)
# shapiro.test(resid(mite.rda2.axis1.env))
# summary(mite.rda2.axis1.env)
# mite.rda2.axis2.env <- lm(mite.rda2.axes[ ,2] ~ ., data = mite.env)
# shapiro.test(resid(mite.rda2.axis2.env))
# summary(mite.rda2.axis2.env)
# mite.rda2.axis3.env <- lm(mite.rda2.axes[ ,3] ~ ., data = mite.env)
# shapiro.test(resid(mite.rda2.axis3.env))
# summary(mite.rda2.axis3.env)
# 
# 
# 
# # Quick MEM
# mite.dbmem.quick <- quickMEM(mite.h, mite.xy)
# summary(mite.dbmem.quick)
# # Eigenvalues
# mite.dbmem.quick[[2]] # OR mite.dbmem.quick$eigenvalues
# # Results of forward selection
# mite.dbmem.quick[[3]] # OR mite.dbmem.quick$fwd.sel
# 
# # Extract and plot RDA results from a quickMEM output (scaling 2)
# plot(mite.dbmem.quick$RDA, scaling = 2)
# sp.scores2 <-
#   scores(mite.dbmem.quick$RDA,
#          choices = 1:2,
#          scaling = 2,
#          display = "sp")
# arrows(0, 0,
#        sp.scores2[ ,1] * 0.9,
#        sp.scores2[ ,2] * 0.9,
#        length = 0,
#        lty = 1,
#        col = "red"
# )
# 
# 
# 
# # eBird_1 MEM
# eBird_1_coords_2 <- eBird %>%
#   filter(Common_Name == birds$Common_Name[1]) %>%
#   select(long_orig, lat_orig)
# head(eBird_1_coords_2)
# head(eBird_1_coords)
# head(eBird_1)
# head(eBird_1)[4]
# eBird_1_dbmem_quick <- quickMEM(eBird_1[4], eBird_1_coords_2)
# summary(eBird_1_dbmem_quick)
# # Eigenvalues
# eBird_1_dbmem_quick[[2]] # OR mite.dbmem.quick$eigenvalues
# # Results of forward selection
# eBird_1_dbmem_quick[[3]] # OR mite.dbmem.quick$fwd.sel
# 
# # Extract and plot RDA results from a quickMEM output (scaling 2)
# plot(eBird_1_dbmem_quick$RDA, scaling = 2)
# sp.scores2 <-
#   scores(eBird_1_dbmem_quick$RDA,
#          choices = 1:2,
#          scaling = 2,
#          display = "sp")
# arrows(0, 0,
#        sp.scores2[ ,1] * 0.9,
#        sp.scores2[ ,2] * 0.9,
#        length = 0,
#        lty = 1,
#        col = "red"
# )




# -----------------------------------------------
# Calculate species shifts, lme4 and inla, mixed effects models
# -----------------------------------------------
# Get species-specific shifts
# Do this with lme4 and compare
# Get species-specific shifts
coefficients_lme4 <- function (commonname){
  eBird_Filtered <- eBird[eBird$Common_Name == commonname ,]
  regression <- lmer(MAD ~ Year + (1|ID), data = eBird_Filtered)
  MAD.Shift <- data.frame(fixef(regression))
  Confint <- confint(regression)
  Common.name <- unique(eBird_Filtered$Common_Name)
  data.frame(Common.name,MAD.Shift[2,], Confint[4], Confint[8])
}

#Create dataframe of slopes == MAD Shifts
Bird_Shifts_lme4 <- lapply(birds$Common_Name, coefficients_lme4)
Bird_Shifts_lme4 <- data.frame(do.call(rbind, Bird_Shifts_lme4))
colnames(Bird_Shifts_lme4) <- c("Common_Name", "MAD_Shift_lme4", "MAD_Shift_lme4_0.025", "MAD_Shift_lme4_0.975")
head(Bird_Shifts_lme4)


# Do this using inla and grid as a random intercept
# coefficients <- function (commonname){
  eBirdtest <- eBird[eBird$Common_Name == birds$Common_Name[1] ,]
  inlares <- inla(MAD ~ Year + f(ID, model = "iid"),
                  data = eBirdtest,
                  control.compute = list(dic = TRUE),
                  control.predictor = list(compute = TRUE))
  DIC <- inlares$dic$dic
  MAD.Shift <- data.frame(inlares$summary.fixed[2,1])
  MAD.Shift.0.025 <- data.frame(inlares$summary.fixed[2,3])
  MAD.Shift.0.975 <- data.frame(inlares$summary.fixed[2,5])
  Common.name <- unique(eBirdtest$Common_Name)
  data.frame(Common.name, DIC, MAD.Shift, MAD.Shift.0.025, MAD.Shift.0.975)
  rm(Common.name, eBirdtest, inlares, DIC, MAD.Shift, MAD.Shift.0.025, MAD.Shift.0.975)
# }
coefficients_iid <- function (commonname){
  eBirdsub <- eBird[eBird$Common_Name == commonname ,]
  inlares <- inla(MAD ~ Year + f(ID, model = "iid"),
                  data = eBirdsub,
                  control.compute = list(dic = TRUE),
                  control.predictor = list(compute = TRUE))
  DIC <- inlares$dic$dic
  MAD.Shift <- data.frame(inlares$summary.fixed[2,1])
  MAD.Shift.0.025 <- data.frame(inlares$summary.fixed[2,3])
  MAD.Shift.0.975 <- data.frame(inlares$summary.fixed[2,5])
  Common.name <- unique(eBirdsub$Common_Name)
  data.frame(Common.name, DIC, MAD.Shift, MAD.Shift.0.025, MAD.Shift.0.975)
}
Bird_Shifts <- lapply(birds$Common_Name, coefficients_iid)
Bird_Shifts <- data.frame(do.call(rbind, Bird_Shifts))
colnames(Bird_Shifts) <- c("Common_Name", "DIC_iid", "MAD_Shift_iid_mean", "MAD_Shift_iid_0.025", "MAD_Shift_iid_0.975")
# Join the lme4 estimates
Bird_Shifts <- left_join(Bird_Shifts_lme4, Bird_Shifts)
# plot the values
ggplot(Bird_Shifts, aes(x = MAD_Shift_lme4, y = MAD_Shift_iid_mean)) +
  geom_point()
cor(Bird_Shifts$MAD_Shift_iid_mean, Bird_Shifts$MAD_Shift_lme4)
# Virtually identical estimates, which is expected.



# -----------------------------------------------
# Calculate species shifts, inla, spatial model
# -----------------------------------------------
# First, define the irregular lattice
# Set radius
radius = 100000
# Extract lat long from the grids
eBird$lat <- str_extract_all(eBird$GridID, "([-]*[0-9]+[x])")
eBird$lat <- as.numeric(str_remove_all(eBird$lat, "[x]"))
eBird$long <- as.numeric(str_extract_all(eBird$GridID, "[0-9]+$"))
eBird_coords <- data.frame(x = eBird$long, y = eBird$lat)
head(eBird_coords)
# define the plot edges based upon the plot radius.
unique_Grids <- eBird %>%
  select(ID, lat, long) %>%
  distinct() %>%
  arrange(-desc(long)) %>%
  arrange(desc(lat))
# Get the edges
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
ID = unique_Grids$ID
# create spatial polygons from coordinates
polys <- SpatialPolygons(mapply(function(poly, id)
{
  xy <- matrix(poly, ncol=2, byrow=TRUE)
  Polygons(list(Polygon(xy)), ID=id)
},
split(square, row(square)), ID),
proj4string = CRS(as.character("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")))
# plot
plot(polys)
# Make a spdf
unique_Grids_2 <- unique_Grids %>%
  select(ID)
unique_Grids_2 <- as.data.frame(unique_Grids_2)
rownames(unique_Grids_2) <- unique_Grids_2$ID
polys_df <- SpatialPolygonsDataFrame(Sr = polys, data = unique_Grids_2)
# plot
plot(polys_df)

# To compute the adjacency, poly2nd can be used
polys_adj <- poly2nb(polys_df)
# By default, it will create a binary adjacency matrix, so that two regions are neighbors only if they share at least one point in common boundary (i.e., it is a queen adjacency).
# Furthermore, binary and row-standardized adjacency matrices will be computed as they will be needed by the spatial models presented below:
W.polys <- nb2mat(polys_adj, style = "B")
W.polys.rs <- nb2mat(polys_adj, style = "W")

# First of all, we previously fir a model with i.i.d. random effects. This will provide a baseline to assess whether spatial random effects are really required when modeling these data.
# The estimates of the random effects (i.e., posterior means) will be added to the SpatialPolygonsDataFrame to be plotted later.

# Besag's
coefficients_besag <- function (commonname){
  eBirdsub <- eBird[eBird$Common_Name == commonname ,]
  inlares <- inla(MAD ~ Year + f(ID, model = "besag", graph = W.polys),
                  data = eBirdsub,
                  control.compute = list(dic = TRUE),
                  control.predictor = list(compute = TRUE))
  DIC <- inlares$dic$dic
  MAD.Shift <- data.frame(inlares$summary.fixed[2,1])
  MAD.Shift.0.025 <- data.frame(inlares$summary.fixed[2,3])
  MAD.Shift.0.975 <- data.frame(inlares$summary.fixed[2,5])
  Common.name <- unique(eBirdsub$Common_Name)
  data.frame(Common.name, DIC, MAD.Shift, MAD.Shift.0.025, MAD.Shift.0.975)
}
Bird_Shifts_besag <- lapply(birds$Common_Name, coefficients_besag)
Bird_Shifts_besag <- data.frame(do.call(rbind, Bird_Shifts_besag))
colnames(Bird_Shifts_besag) <- c("Common_Name", "DIC_besag", "MAD_Shift_besag_mean", "MAD_Shift_besag_0.025", "MAD_Shift_besag_0.975")
# Join the besag estimates
Bird_Shifts <- left_join(Bird_Shifts, Bird_Shifts_besag)

# besagproper's
coefficients_besagproper <- function (commonname){
  eBirdsub <- eBird[eBird$Common_Name == commonname ,]
  inlares <- inla(MAD ~ Year + f(ID, model = "besagproper", graph = W.polys),
                  data = eBirdsub,
                  control.compute = list(dic = TRUE),
                  control.predictor = list(compute = TRUE),
                  verbose = T)
  DIC <- inlares$dic$dic
  MAD.Shift <- data.frame(inlares$summary.fixed[2,1])
  MAD.Shift.0.025 <- data.frame(inlares$summary.fixed[2,3])
  MAD.Shift.0.975 <- data.frame(inlares$summary.fixed[2,5])
  Common.name <- unique(eBirdsub$Common_Name)
  data.frame(Common.name, DIC, MAD.Shift, MAD.Shift.0.025, MAD.Shift.0.975)
}
Bird_Shifts_besagproper <- lapply(birds$Common_Name, coefficients_besagproper)
Bird_Shifts_besagproper <- data.frame(do.call(rbind, Bird_Shifts_besagproper))
colnames(Bird_Shifts_besagproper) <- c("Common_Name", "DIC_besagproper", "MAD_Shift_besagproper_mean", "MAD_Shift_besagproper_0.025", "MAD_Shift_besagproper_0.975")
# Join the besagproper estimates
Bird_Shifts <- left_join(Bird_Shifts, Bird_Shifts_besagproper)

# bym
coefficients_bym <- function (commonname){
  eBirdsub <- eBird[eBird$Common_Name == commonname ,]
  inlares <- inla(MAD ~ Year + f(ID, model = "bym", graph = W.polys),
                  data = eBirdsub,
                  control.compute = list(dic = TRUE),
                  control.predictor = list(compute = TRUE))
  DIC <- inlares$dic$dic
  MAD.Shift <- data.frame(inlares$summary.fixed[2,1])
  MAD.Shift.0.025 <- data.frame(inlares$summary.fixed[2,3])
  MAD.Shift.0.975 <- data.frame(inlares$summary.fixed[2,5])
  Common.name <- unique(eBirdsub$Common_Name)
  data.frame(Common.name, DIC, MAD.Shift, MAD.Shift.0.025, MAD.Shift.0.975)
}
Bird_Shifts_bym <- lapply(birds$Common_Name, coefficients_bym)
Bird_Shifts_bym <- data.frame(do.call(rbind, Bird_Shifts_bym))
colnames(Bird_Shifts_bym) <- c("Common_Name", "DIC_bym", "MAD_Shift_bym_mean", "MAD_Shift_bym_0.025", "MAD_Shift_bym_0.975")
# Join the bym estimates
Bird_Shifts <- left_join(Bird_Shifts, Bird_Shifts_bym)

# # ler
# Q <- Diagonal(x = sapply(polys_adj, length))
# Q <- Q - as_dsTMatrix_listw(nb2listw(polys_adj, style = "B"))
#
# C <- Diagonal(x = 1, n = nrow(boston.tr)) - Q



# print out the DICs
DIC_summary <- Bird_Shifts %>%
  select(Common_Name, DIC_iid, DIC_besag, DIC_bym) %>%
  mutate(DIC_min = pmin(DIC_iid, DIC_besag, DIC_bym)) %>%
  mutate(DIC_iid_del = DIC_iid - DIC_min,
         DIC_besag_del = DIC_besag - DIC_min,
         DIC_bym_del = DIC_bym - DIC_min)


DIC_summary_2 <- DIC_summary %>%
  select(Common_Name,
         DIC_iid_del,
         DIC_besag_del,
         DIC_bym_del) %>%
  pivot_longer(cols = c(DIC_iid_del,
                        DIC_besag_del,
                        DIC_bym_del),
               names_to = "Model",
               values_to = "DIC_del")

# Order the DICs for processing
DIC_summary_2 <- DIC_summary_2 %>%
  group_by(Common_Name) %>%
  arrange(-desc(DIC_del), .by_group = T) %>%
  ungroup()

# Get birds where delta DIC is greater than 2
DIC_summary_3 <-DIC_summary_2 %>%
  group_by(Common_Name) %>%
  filter(DIC_del != max(DIC_del)) %>%
  filter(DIC_del > 2) %>%
  ungroup()
DIC_summary_3
write_csv(DIC_summary_3, "C:/Users/dpgil/Documents/1 U of T/eBird/Figures/Supplementary analysis/INLA/DIC_summary_3.csv")


# Correlation among MAD values
cor(select(Bird_Shifts, MAD_Shift_lme4, MAD_Shift_iid_mean, MAD_Shift_besag_mean, MAD_Shift_bym_mean))
cor(select(Bird_Shifts, MAD_Shift_lme4_0.025, MAD_Shift_iid_0.025, MAD_Shift_besag_0.025, MAD_Shift_bym_0.025))
cor(select(Bird_Shifts, MAD_Shift_lme4_0.975, MAD_Shift_iid_0.975, MAD_Shift_besag_0.975, MAD_Shift_bym_0.975))

# Plot a summary of this
# Rearrange data frame
# Bird_Shifts_tidy_1 <- Bird_Shifts %>%
#   pivot_longer(cols = c(MAD_Shift_lme4, MAD_Shift_iid_mean, MAD_Shift_besag_mean, MAD_Shift_bym_mean),
#                names_to = "Model",
#                values_to = "MAD_Shift")
# iid
Bird_Shifts_iid_2 <- Bird_Shifts %>%
  select(Common_Name, DIC_iid, MAD_Shift_iid_mean, MAD_Shift_iid_0.025, MAD_Shift_iid_0.975) %>%
  rename("DIC" = "DIC_iid", "MAD_Shift_mean" = "MAD_Shift_iid_mean", "MAD_Shift_0.025" = "MAD_Shift_iid_0.025", "MAD_Shift_0.975" = "MAD_Shift_iid_0.975") %>%
  mutate(Model = "iid")
# besag
Bird_Shifts_besag_2 <- Bird_Shifts %>%
  select(Common_Name, DIC_besag, MAD_Shift_besag_mean, MAD_Shift_besag_0.025, MAD_Shift_besag_0.975) %>%
  rename("DIC" = "DIC_besag", "MAD_Shift_mean" = "MAD_Shift_besag_mean", "MAD_Shift_0.025" = "MAD_Shift_besag_0.025", "MAD_Shift_0.975" = "MAD_Shift_besag_0.975") %>%
  mutate(Model = "besag")
# bym
Bird_Shifts_bym_2 <- Bird_Shifts %>%
  select(Common_Name, DIC_bym, MAD_Shift_bym_mean, MAD_Shift_bym_0.025, MAD_Shift_bym_0.975) %>%
  rename("DIC" = "DIC_bym", "MAD_Shift_mean" = "MAD_Shift_bym_mean", "MAD_Shift_0.025" = "MAD_Shift_bym_0.025", "MAD_Shift_0.975" = "MAD_Shift_bym_0.975") %>%
  mutate(Model = "bym")
# lme4
Bird_Shifts_lme4_2 <- Bird_Shifts_lme4 %>%
  mutate(DIC_lme4 = NA) %>%
  select(Common_Name, DIC_lme4, MAD_Shift_lme4, MAD_Shift_lme4_0.025, MAD_Shift_lme4_0.975) %>%
  rename("DIC" = "DIC_lme4", "MAD_Shift_mean" = "MAD_Shift_lme4", "MAD_Shift_0.025" = "MAD_Shift_lme4_0.025", "MAD_Shift_0.975" = "MAD_Shift_lme4_0.975") %>%
  mutate(Model = "lme4")
# bind
Bird_Shifts_tidy <- rbind(Bird_Shifts_iid_2, Bird_Shifts_besag_2, Bird_Shifts_bym_2, Bird_Shifts_lme4_2)

ggplot(Bird_Shifts_tidy) +
  geom_point(aes(x = Common_Name, y = MAD_Shift_mean, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = MAD_Shift_0.025, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = MAD_Shift_0.975, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

Bird_Shifts_tidy %>%
  filter(Common_Name %in% DIC_summary_3$Common_Name) %>%
  ggplot() +
  geom_point(aes(x = Common_Name, y = MAD_Shift_mean, colour = Model), shape = 1) +
  geom_point(aes(x = Common_Name, y = MAD_Shift_0.025, colour = Model), shape = 6) +
  geom_point(aes(x = Common_Name, y = MAD_Shift_0.975, colour = Model), shape = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

Bird_Shifts_tidy %>%
  filter(Common_Name %in% DIC_summary_3$Common_Name) %>%
  group_by(Common_Name) %>%
  arrange(desc(Common_Name), .by_group = T) %>%
  ungroup()

write_csv(Bird_Shifts, "C:/Users/dpgil/Documents/1 U of T/eBird/Figures/Supplementary analysis/INLA/Bird_Shifts.csv")
write_csv(Bird_Shifts_tidy, "C:/Users/dpgil/Documents/1 U of T/eBird/Figures/Supplementary analysis/INLA/Bird_Shifts_tidy.csv")
