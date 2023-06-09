# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 7G
# ------------------------------------------------
# Load Libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# Set Working Directory
setwd("~/Data repository")

# ------------------------------------------------
# Read in processed data
# ------------------------------------------------
Routine_7_files <- read_csv("Routine 7/Merged output/MAD_Output_All_F.csv")
spp_orig <- read_csv("Routine 7/Merged output/spp_orig.csv")

# -------------------------------------------------------
# Re-visit each filter (except JdDetect and CI)
# -------------------------------------------------------
# Filter for species-grids that have number of years of at least 10 years
Routine_7_files <- Routine_7_files %>% 
  select(-nyrs, -min_Year, -max_Year, -Duration_Years, -Spp_n_GridID, -GridID_n_spp, -ngrids)

Species_Grid_nyrs <- Routine_7_files %>% 
  group_by(Common_Name, GridID, Jd_min) %>% 
  summarize(nyrs = n()) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Routine_7_files <- left_join(Routine_7_files, Species_Grid_nyrs, by = c("Common_Name", "GridID", "Jd_min"))
Routine_7_files <- Routine_7_files %>%
  filter(nyrs >= 10)
nrow(Routine_7_files)
# 144740
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# 43 species removed

# Remember to group for Jd_min, because you haven't decided on that yet
Species_Grid_Duration_Years <- Routine_7_files %>% 
  group_by(Common_Name, GridID, Jd_min) %>% 
  summarize(min_Year = min(Year),
            max_Year = max(Year),
            Duration_Years = max_Year - min_Year + 1) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Routine_7_files <- left_join(Routine_7_files, Species_Grid_Duration_Years, by = c("Common_Name", "GridID", "Jd_min"))

# Filter for species-grids that have a duration of at least 15 years
Routine_7_files <- Routine_7_files %>% 
  filter(Duration_Years >= 15)
nrow(Routine_7_files)
# 144740
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# 43 species removed

# Remove species that have less than 7 grids (i.e., 1-6)
Species_Grid_n <- Routine_7_files %>% 
  select(Common_Name, GridID, Jd_min) %>%
  distinct() %>% 
  group_by(Common_Name, Jd_min) %>% 
  count(name = "Spp_n_GridID") %>% 
  ungroup()
Routine_7_files <- left_join(Routine_7_files, Species_Grid_n, by = c("Common_Name", "Jd_min"))
Routine_7_files <- Routine_7_files %>% 
  filter(Spp_n_GridID >= 6)
# Check # of rows
nrow(Routine_7_files)
# 144740
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name))) %>% print(n = 32)
# 43 species removed

# Now, we want to see how many species there are across each of the grids
# Look at the maximum number of species that occur in each grid cell
# Remember to group for Jd_min, because you haven't decided on that yet
Species_Grid_Max_n <- Routine_7_files %>% 
  select(Common_Name, GridID, Jd_min) %>% 
  distinct() %>% 
  group_by(GridID, Jd_min) %>% 
  count(name = "GridID_n_spp") %>% 
  ungroup()
# Remove grids that only have one species
Routine_7_files <- left_join(Routine_7_files, Species_Grid_Max_n, by = c("GridID", "Jd_min"))
Routine_7_files <- Routine_7_files %>% 
  filter(GridID_n_spp > 1)
# Check # of rows
nrow(Routine_7_files)
# 144740
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# 43 species removed

# Check out the number of grids that a species has for a year
Species_Years <- Routine_7_files %>% 
  group_by(Common_Name, Year, Jd_min) %>% 
  summarize(ngrids = n()) %>% 
  ungroup()
# Remove years that have less than 7 grids
Routine_7_files <- left_join(Routine_7_files, Species_Years, by = c("Common_Name", "Year", "Jd_min"))
Routine_7_files <- Routine_7_files %>% 
  filter(ngrids >= 6)
nrow(Routine_7_files)
# 144740
anti <- anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))

# # write the progress *** Do note write because there was no rows that were removed
# write_csv(Routine_7_files, "Routine 7/Merged output/MAD_Output_All_G.csv")
rm(Species_Grid_Duration_Years, Species_Years, Species_Grid_Max_n, Species_Grid_n, Species_Grid_nyrs)
gc()


# # ------------------------------------------------
# # Read in min_JD_summary 
# # ------------------------------------------------
# min_JD_summary <- read_csv("Routine 7/Merged output/min_JD_summary.csv")
# min_JD_summary <- min_JD_summary %>% 
#   rename("Jd_min_top" = "Jd_min")
# # Join to main data.frame
# Routine_7_files <- left_join(Routine_7_files, min_JD_summary)
# # # Analyze MAD_CI_range and Prsq
# # Quality_summary_pre_QAQC <- Routine_7_files %>% 
# #   group_by(Common_Name, Jd_min, Jd_min_top) %>% 
# #   summarize(MAD_CI_range_mean = mean(MAD_CI_range),
# #             MAD_CI_range_sd = sd(MAD_CI_range),
# #             Prsq_mean = mean(Prsq),
# #             Prsq_sd = sd(Prsq)) %>% 
# #   ungroup()
# # ggplot(Quality_summary_pre_QAQC, aes(x = Common_Name, y = Prsq_mean, fill = factor(Jd_min))) + 
# #   geom_col(position = position_dodge()) +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
# # MAD_CI_range
# Routine_7_files_pivot_fits <- Routine_7_files %>%
#   select(Common_Name, GridID, Year, MAD_CI_range, Jd_min) %>%
#   pivot_wider(names_from = Jd_min, values_from = MAD_CI_range, names_prefix = "MAD_CI_range_")
# 
# cor_result <- cor(Routine_7_files_pivot_fits[4:14], use = "pairwise.complete.obs")
# cor_result_0.9 <- cor_result > 0.9
# sum(cor_result_0.9, na.rm = TRUE) 
# nrow(cor_result) ; ncol(cor_result_0.9) ; nrow(cor_result_0.9) * ncol(cor_result_0.9)
# (sum(cor_result_0.9, na.rm = TRUE))  / (nrow(cor_result_0.9) * ncol(cor_result_0.9)) * 100
# # Prsq
# Routine_7_files_pivot_fits_2 <- Routine_7_files %>%
#   select(Common_Name, GridID, Year, Prsq, Jd_min) %>%
#   pivot_wider(names_from = Jd_min, values_from = Prsq, names_prefix = "Prsq_")
# cor_result_2 <- cor(Routine_7_files_pivot_fits_2[4:14], use = "pairwise.complete.obs")
# cor_result_2_0.9 <- cor_result_2 > 0.9
# sum(cor_result_2_0.9, na.rm = TRUE) 
# (sum(cor_result_2_0.9, na.rm = TRUE))  / (nrow(cor_result_2_0.9) * ncol(cor_result_2_0.9)) * 100
# # cor_result_df_2 <- data.frame(Jd_min_class = row.names(cor_result), cor_result)
# # cor_result

# ------------------------------------------------
# Assess MinJd overall and by species
# ------------------------------------------------
# **** Take into account JdDetect to break ties in the number of blocks

# How many species grid combinations exist across all Jd_min iterations? Are some blocks only around for certain Jd_min values, or is the maximum number of Jd_min values inclusive of all blocks
Blocks_summary_Jd_min <- Routine_7_files %>%
  group_by(Common_Name, GridID, Year) %>%
  summarize(block_count_all_Jd_min = n()) %>%
  ungroup()
Blocks_summary_Jd_min %>%
  group_by(block_count_all_Jd_min) %>%
  count() %>%
  ungroup() %>%
  arrange(desc(n))
Blocks_summary_Jd_min_2 <- Routine_7_files %>%
  select(Common_Name, GridID, Year) %>%
  distinct() %>%
  group_by(Common_Name) %>%
  summarize(block_sum_all_Jd_min = n()) %>%
  ungroup()
Blocks_summary_Jd_min_3 <- Routine_7_files %>%
  group_by(Common_Name, Jd_min) %>%
  summarize(block_sum_by_Jd_min = n()) %>%
  ungroup()
Blocks_summary_Jd_min_4 <- Routine_7_files %>%
  group_by(Common_Name, Jd_min) %>%
  summarize(spp_Jd_min_JdDetect = sum(JdDetect)) %>%
  ungroup()

Routine_7_files <- left_join(Routine_7_files, Blocks_summary_Jd_min) %>%
  left_join(Blocks_summary_Jd_min_2) %>%
  left_join(Blocks_summary_Jd_min_3) %>%
  left_join(Blocks_summary_Jd_min_4)


(Jd_min_summary_1 <- Routine_7_files %>%
    group_by(Jd_min) %>%
    count() %>%
    ungroup() %>%
    arrange(desc(n)) %>%
    mutate(n_pct_max = n/max(n)*100) %>%
    mutate(n_pct_diff_max = max(n_pct_max) - n_pct_max))
# 50, 70, and 60 are almost the same, followed by 40, 30, 20, and 80. 90, 10, and 1 are much worse, and 100 much more so
Jd_min_summary_2 <- Routine_7_files %>%
  mutate(block_sum_pct = block_sum_by_Jd_min / block_sum_all_Jd_min * 100) %>%
  # group_by(Common_Name, Jd_min, block_sum_all_Jd_min, block_sum_by_Jd_min, block_sum_pct) %>%
  # count() %>%
  # ungroup() %>%
  select(Common_Name, Jd_min, block_sum_all_Jd_min, block_sum_by_Jd_min, block_sum_pct, spp_Jd_min_JdDetect) %>%
  distinct() %>%
  group_by(Common_Name) %>%
  mutate(n_pct_max = block_sum_by_Jd_min/max(block_sum_by_Jd_min)*100) %>%
  ungroup() %>%
  group_by(Common_Name) %>%
  mutate(n_pct_diff_max = max(n_pct_max) - n_pct_max) %>%
  ungroup()
# Does block_sum_by_Jd_min correlate with spp_Jd_min_JdDetect?
cor(Jd_min_summary_2$block_sum_by_Jd_min, Jd_min_summary_2$spp_Jd_min_JdDetect)
# Typically, yes 0.9302545
# Make sure that there aren't any high-value block-sums that are low JdDetects
# 1-15
Jd_min_summary_2 %>%
  filter(Common_Name %in% unique(Routine_7_files$Common_Name)[1:15]) %>%
  ggplot(aes(block_sum_by_Jd_min, spp_Jd_min_JdDetect)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~Common_Name, scales = "free")
# 16-30
Jd_min_summary_2 %>%
  filter(Common_Name %in% unique(Routine_7_files$Common_Name)[16:30]) %>%
  ggplot(aes(block_sum_by_Jd_min, spp_Jd_min_JdDetect)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~Common_Name, scales = "free")
# 31-44
Jd_min_summary_2 %>%
  filter(Common_Name %in% unique(Routine_7_files$Common_Name)[31:44]) %>%
  ggplot(aes(block_sum_by_Jd_min, spp_Jd_min_JdDetect)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~Common_Name, scales = "free")
# # 76-87
# Jd_min_summary_2 %>%
#   filter(Common_Name %in% spp$common_name[76:87]) %>%
#   ggplot(aes(block_sum_by_Jd_min, spp_Jd_min_JdDetect)) +
#   geom_point(alpha = 0.6) +
#   facet_wrap(~Common_Name, scales = "free")


# Look at each species
Jd_min_summary_3 <- Jd_min_summary_2 %>%
  filter(n_pct_diff_max <= 1)
# How many are 90 or 100?
Jd_min_summary_3 %>%
  filter(Jd_min > 80) %>%
  nrow()
Jd_min_summary_3 %>%
  filter(Jd_min > 80)

# What is the correlation between MAD for each Jd_min? For each Jd_min for each species?
Routine_7_files_pivot <- Routine_7_files %>%
  select(Common_Name, GridID, Year, MAD, Jd_min) %>%
  pivot_wider(names_from = Jd_min, values_from = MAD, names_prefix = "MAD_")
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[1,]
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[2,]
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[3,]
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[4,]
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[5,]
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[6,]
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[7,]
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[8,]
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[9,]
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[10,]
cor(Routine_7_files_pivot[4:14], use = "pairwise.complete.obs")[11,]

# Filter for each species
Jd_min_summary_4 <- Jd_min_summary_2 %>%
  filter(n_pct_diff_max == 0)
nrow(Jd_min_summary_4)
# Which speices?
Jd_min_summary_4 %>%
  group_by(Common_Name) %>%
  count() %>%
  filter(n > 1)
# There are some ties, go with the one that has more JdDetect
Jd_min_summary_5 <- Jd_min_summary_4 %>%
  group_by(Common_Name) %>%
  filter(spp_Jd_min_JdDetect == max(spp_Jd_min_JdDetect)) %>%
  ungroup()
nrow(Jd_min_summary_5)
Jd_min_summary_5 %>%
  group_by(Common_Name) %>%
  count() %>%
  filter(n > 1)
# No more ties
# # There are still some ties, go with the one that has a higher Jd_min, because the reduces the likelihood of influence from early outliers
# Jd_min_summary_6 <- Jd_min_summary_5 %>%
#   group_by(Common_Name) %>%
#   filter(Jd_min == max(Jd_min)) %>%
#   ungroup()
# nrow(Jd_min_summary_6)



# ----------------------------------
# Filter the overall data for the "best" subset of Jd_min
# ----------------------------------
Jd_min_summary_5_select <- Jd_min_summary_5 %>%
  select(Common_Name, Jd_min)
Routine_7_files_top <- left_join(Jd_min_summary_5_select, Routine_7_files)
write_csv(Routine_7_files_top, "Routine 7/Merged output/MAD_Output_All_Filter_pre.csv")
