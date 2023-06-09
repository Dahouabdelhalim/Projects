# ----------------------------------------------
# Analyzing results of eBird Model Testing
# ----------------------------------------------
# Load libraries
library(dplyr)
library(ggplot2) 
library(readr) 
library(tibble) 
library(tidyr) 
library(stringr) 
library(corrr)
library(GGally)

# Set Working Directory
setwd("~/Data repository")



# ----------------------------------------------
# Read in files
# ----------------------------------------------
# F is the one where you have everything filtered down, but all of the Jd_min values. Only a few observations are removed afterwards because of grids and Jd_min
eBird_data <- read_csv("Routine 7/Merged output/MAD_Output_All_F.csv")
eBird_data_grid_spp <- eBird_data %>% 
  select(GridID, Common_Name) %>% 
  distinct() %>% 
  count(GridID)
top_grids <- eBird_data_grid_spp %>% 
  filter(n == max(n))
top_grids_2 <- eBird_data %>% 
  filter(GridID %in% top_grids$GridID)
# Which one is longer?
top_grids_2 %>% count(GridID)
# Grid_300000x700000   

eBird_data <- eBird_data %>%
  filter(GridID == "Grid_300000x700000")



# ----------------------------------------------
# Data exploration and statistical analyses
# ----------------------------------------------
# Plot the data
set.seed(123)
ggplot(data = eBird_data, aes(x = Jd_min, y = MAD)) +
  geom_jitter(alpha = 0.5, aes(colour = Common_Name))
set.seed(123)
ggplot(data = eBird_data, aes(x = Jd_min, y = MAD_CI_range)) +
  geom_jitter(alpha = 0.5, aes(colour = Common_Name))
set.seed(123)
ggplot(data = eBird_data, aes(x = Jd_min, y = MAD_CI_range)) +
  geom_jitter(alpha = 0.5, aes(colour = Common_Name)) +
  scale_y_log10()
set.seed(123)
ggplot(data = eBird_data, aes(x = Jd_min, y = MAD)) +
  geom_jitter(alpha = 0.5, aes(colour = Common_Name)) +
  facet_wrap(~Common_Name) +
  theme(legend.position = "none")
ggplot(data = eBird_data, aes(x = factor(Jd_min), y = MAD)) +
  geom_boxplot() +
  # facet_wrap(~Common_Name) +
  theme(legend.position = "none")
ggplot(data = eBird_data, aes(x = factor(Jd_min), y = MAD)) +
  geom_boxplot(aes(colour = Common_Name)) +
  facet_wrap(~Common_Name) +
  theme(legend.position = "none")

# Summarize the data
eBird_data_summary <- eBird_data %>% 
  group_by(Jd_min, Common_Name) %>% 
  summarize(median_MAD = median(MAD),
            mean_MAD = mean(MAD),
            median_MAD_CI_range = median(MAD_CI_range),
            mean_MAD_CI_range = mean(MAD_CI_range),
            n_worked = n()) %>% 
  ungroup()
eBird_data_summary <- eBird_data_summary[order(eBird_data_summary$Common_Name),]

# Conduct ANOVA
ANOVA <- aov(MAD ~ Jd_min, data = eBird_data)
summary(ANOVA)
pair_test_test <- pairwise.t.test(eBird_data$MAD, eBird_data$Jd_min, p.adjust.method = "BH")
pair_test_test
pair_test_test_df <- as.data.frame(pair_test_test$p.value)
pair_test_test_df <- rownames_to_column(pair_test_test_df)
pair_test_test_df <- pair_test_test_df %>% 
  mutate(across(where(is.numeric), round, digits = 3))
write_csv(pair_test_test_df, "pairwise_t_tests_rev_spp.csv")
P_values <- pair_test_test$p.value
P_values_non_NA <- colSums(!is.na(P_values))
sum(P_values_non_NA)
pair_test_test$p.value > 0.05
(sum(pair_test_test$p.value > 0.05, na.rm = TRUE))  / sum(P_values_non_NA) * 100

# Pairwise correlations
eBird_data_cor <- eBird_data %>% 
  select(Year, Common_Name, MAD, Jd_min) %>% 
  pivot_wider(names_from = Jd_min, values_from = MAD)
correlate(eBird_data_cor[,3:13])
ggpairs(eBird_data_cor[,3:13], mapping = aes(alpha = 0.5))



# ----------------------------------------------
# NA analysis
# ----------------------------------------------
# Sum the NAs in the summary file
eBird_data_all <- eBird_data_cor %>% 
  select(Common_Name) %>% 
  distinct()
Years <- data.frame(Year = 2002:2019)
eBird_data_all <- merge(eBird_data_all, Years)
# Replace NAs with 0
eBird_data_cor <- eBird_data_cor %>% 
  mutate_all(funs(str_replace_na(., 0)))
eBird_data_cor$Year <- as.numeric(eBird_data_cor$Year)
eBird_data_cor_all <- left_join(eBird_data_all, eBird_data_cor)
eBird_zero_sum <- eBird_data_cor_all %>% 
  select(-Year, -Common_Name) %>% 
  summarise_all(funs(sum(.==0, na.rm=TRUE)))
eBird_zero_sum_long <- eBird_zero_sum %>% 
  gather(key = "Jd_min", value = "Zeroes")
# Group by spp
eBird_zero_sum_spp <- eBird_data_cor_all %>% 
  group_by(Common_Name) %>% 
  select(-Year, -Common_Name) %>% 
  summarise_all(funs(sum(.==0, na.rm=TRUE)))
eBird_zero_sum_long_spp <- eBird_zero_sum_spp %>%
  pivot_longer(cols = `1`:`90`, names_to = "Jd_min", values_to = "Zeroes")

# Get NA sum so we know how many years failed for all iterations
eBird_NA_sum <- eBird_data_cor_all %>% 
  select(-Year, -Common_Name) %>% 
  summarise_all(funs(sum(is.na(.))))
eBird_NA_sum <- data.frame(NA_sum = eBird_NA_sum$`1`)
eBird_NA_sum$Bird_Year_Sum <- nrow(eBird_data_all)
eBird_NA_sum$Bird_Year_Sum_Prop <- 1 - eBird_NA_sum$NA_sum/eBird_NA_sum$Bird_Year_Sum
eBird_NA_sum$GridID <- unique(eBird_data$GridID)



# Convert to ordered factor
eBird_zero_sum_long$Jd_min <- factor(eBird_zero_sum_long$Jd_min, levels = c("1", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), ordered = TRUE)
# Plot the 0s
ggplot(data = eBird_zero_sum_long, aes(x = Jd_min, y = Zeroes)) +
  geom_col()

# Convert to ordered factor
eBird_zero_sum_long_spp$Jd_min <- factor(eBird_zero_sum_long_spp$Jd_min, levels = c("1", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"), ordered = TRUE)
# Plot the 0s
ggplot(data = eBird_zero_sum_long_spp, aes(x = Jd_min, y = Zeroes, fill = Common_Name)) +
  geom_col(position = position_dodge())

# What is the minimum for each species?
eBird_zero_sum_long_spp_min <- eBird_zero_sum_long_spp %>% 
  group_by(Common_Name) %>% 
  slice(which.min(Zeroes))
eBird_zero_sum_long_spp_min

# Look at histograms of the MAD (to look for outliers)
ggplot(data = eBird_data, aes(x = MAD)) +
  geom_histogram(binwidth = 5) +
  facet_wrap(~Jd_min)



# ----------------------------------------------
# Mode deviation analysis
# ----------------------------------------------
# Look at the deviations from the most common MAD or Mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Create a MAD that is rounded to one decimal place
# Find modes
eBird_data <- eBird_data %>% 
  mutate(MAD_round_0 = round(MAD, digits = 0)) 
eBird_data_mode <- eBird_data %>% 
  group_by(Common_Name, Year) %>% 
  summarize(MAD_Mode = getmode(MAD),
            MAD_round_0_Mode = getmode(MAD_round_0)) %>% 
  ungroup()
# Join the mode columns to the ebird_data
eBird_data <- left_join(eBird_data, eBird_data_mode)
# Create a column that calculates the difference between the MAD for a given bird-year for a given Jd_min trial
# 0 = not different from mode
# positive = MAD_Mode was later than MAD for that bird-year trial
# negative = MAD_Mode was earlier than MAD for that bird-year trial
# Create an absolute column as well
eBird_data <- eBird_data %>% 
  mutate(MAD_Mode_Diff = MAD_Mode - MAD,
         MAD_Mode_Diff_Abs = abs(MAD_Mode_Diff),
         MAD_round_0_Mode_Diff = MAD_round_0_Mode - MAD_round_0,
         MAD_round_0_Mode_Diff_Abs = abs(MAD_round_0_Mode_Diff))

# Get the sum of zeros
eBird_zero_sum_mode <- eBird_data %>% 
  select(Jd_min, MAD_Mode_Diff, MAD_round_0_Mode_Diff) %>% 
  group_by(Jd_min) %>% 
  summarise_at(vars(MAD_Mode_Diff, MAD_round_0_Mode_Diff), ~sum(. == 0))

# Plot the data
# Histogram
ggplot(eBird_data, aes(x = MAD_Mode_Diff)) + 
  geom_histogram(binwidth = 3) +
  facet_wrap(~Jd_min)
ggplot(eBird_data, aes(x = MAD_Mode_Diff)) + 
  geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-25, 25)) +
  facet_wrap(~Jd_min)
# Sum of zeros bar plot
ggplot(eBird_zero_sum_mode, aes(x = Jd_min, y = MAD_Mode_Diff)) +
  geom_col()
ggplot(eBird_zero_sum_mode, aes(x = Jd_min, y = MAD_round_0_Mode_Diff)) +
  geom_col()
