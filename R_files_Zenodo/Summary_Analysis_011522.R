# -------------------------------------------------------------------------
# Analyzing eBird data mean arrival dates and data quantity/quality
# -------------------------------------------------------------------------
# Load Libraries
library(tidyverse)
library(lme4)
library(GGally)
library(ape)

# Set Working Directory
setwd("~/Data repository")


# ---------------------------------
# Read in data
# ---------------------------------
Routine_7_files_top <- read_csv("Routine 7/Merged Output/MAD_Output_All_Filter_1.csv")

# ------------------------------------------------
# Read in species names
# ------------------------------------------------
spp <- read_csv("MorphologyParameters.csv")

# ----------------------------------
# Minimum JD summary
# ----------------------------------
min_JD_summary <- eBird_Filtered %>% 
  group_by(Common_Name, Jd_min) %>% 
  summarize(Mean_MAD = mean(MAD)) %>% 
  ungroup()
plot(min_JD_summary$Jd_min, min_JD_summary$Mean_MAD)
cor(min_JD_summary$Jd_min, min_JD_summary$Mean_MAD)
min_JD_summary %>% 
  count(Jd_min) %>% 
  arrange(desc(n))
write_csv(min_JD_summary, "min_JD_summary.csv")


# ----------------------------------
# Model fit summary
# ----------------------------------
# Summarize by year
Yearly_coverage <- eBird_Filtered %>% 
  group_by(Year) %>% 
  summarize(Year_count = n())
ggplot(Yearly_coverage, aes(Year, Year_count)) +
  geom_point() +
  scale_x_continuous(n.breaks = 10)

# Species year summary
Species_coverage <- eBird_Filtered %>% 
  group_by(Common_Name, Year) %>% 
  summarize(Common_Name_count = n()) %>% 
  ungroup()
mean(Species_coverage$Common_Name_count)
sd(Species_coverage$Common_Name_count)

# Number of days detected
mean(eBird_Filtered$JdDetect)
sd(eBird_Filtered$JdDetect)

# MAD CI range
hist(eBird_Filtered$MAD_CI_range)
quantile(eBird_Filtered$MAD_CI_range)
# below 10
MAD_CI_below_10 <- eBird_Filtered %>% 
  filter(MAD_CI_range < 10)
(nrow(MAD_CI_below_10)/nrow(eBird_Filtered))*100
mean(MAD_CI_below_10$JdDetect)
sd(MAD_CI_below_10$JdDetect)
summary(MAD_CI_below_10$JdDetect)
# Above 10
MAD_CI_above_10 <- eBird_Filtered %>% 
  filter(MAD_CI_range >= 10)
mean(MAD_CI_above_10$JdDetect)
sd(MAD_CI_above_10$JdDetect)
summary(MAD_CI_above_10$JdDetect)

# prsq
mean(eBird_Filtered$Prsq)
median(eBird_Filtered$Prsq)
sd(eBird_Filtered$Prsq)
summary(eBird_Filtered$Prsq)



# -----------------------------------------------
# Calculate species shifts
# -----------------------------------------------
# Get unique species names
x <- unique(eBird_Filtered$Common_Name)
testdf <- filter(eBird_Filtered, Site_Mean != 2.069)
MAD_Shift_list <- list()
for (i in 1:length(x)) {
# i = 1
  eBird_subset <- testdf[testdf$Common_Name == x[i] ,]
  eBird_subset_lmer <- lmer(MAD ~ Year + (1|GridID), data = eBird_subset)
  eBird_subset_ci <- confint.merMod(eBird_subset_lmer)
  MAD_Shift <- data.frame(unique(eBird_subset$Common_Name),
                          length(unique(eBird_subset$GridID)),
                          fixef(eBird_subset_lmer)[[2]],
                          eBird_subset_ci[4,1],
                          eBird_subset_ci[4,2])
  colnames(MAD_Shift) <- c("Common_Name", "n_Grids", "MAD_Shift", "MAD_Shift_2_5_pct", "MAD_Shift_97_5_pct")
  MAD_Shift_list[[i]] <- MAD_Shift
}

Bird_Shifts <- bind_rows(MAD_Shift_list)

# Join the spp data to Bird_Shifts
Bird_Shifts <- left_join(Bird_Shifts, spp, by = c("Common_Name" = "Common.name"))

summary(round(Bird_Shifts$MAD_Shift, digits = 3))





# -----------------------------------------------
# Calculate mean arrival date shifts by Grid
# -----------------------------------------------
# Get grid-specific shifts
Grid_coefficients <- function (gridid){
  eBird_Filtered <- eBird_Filtered[eBird_Filtered$GridID == gridid ,]
  regression <- lmer(MAD ~ Year + (1|Common_Name), data = eBird_Filtered)
  MAD.Shift <- data.frame(fixef(regression))
  grid.id <- unique(eBird_Filtered$GridID)
  data.frame(grid.id, MAD.Shift[2,], MAD.Shift[1,])
}

#Create dataframe of slopes == MAD Shifts
Grids <- unique(eBird_Filtered$GridID)
Grid_Shifts_lmer <- lapply(Grids, Grid_coefficients)
Grid_Shifts_lmer <- data.frame(do.call(rbind, Grid_Shifts_lmer))
colnames(Grid_Shifts_lmer) <- c("GridID", "MAD_Shift", "MAD_Intercept")
# 397.43914 - intercept for all of the data
summary(Grid_Shifts_lmer$MAD_Intercept)



# ---------------------------------------------------------------
# Identify species-Grid combinations for passing vs. staying
# ---------------------------------------------------------------
passage_grid_data <- eBird_Filtered %>%
  group_by(Common_Name, GridID) %>%
  summarize(mean_MaxJd = mean(MaxJd),
            median_MaxJd = median(MaxJd),
            max_MaxJd = max(MaxJd),
            min_MaxJd = min(MaxJd)) %>%
  ungroup()
# Plot the data by species (median)
ggplot(passage_grid_data, aes(median_MaxJd), colour = "black") +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  scale_x_continuous(n.breaks = 6) +
  facet_wrap(~Common_Name)
# Plot the data by species (mean)
ggplot(passage_grid_data, aes(mean_MaxJd), colour = "black") +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  scale_x_continuous(n.breaks = 6) +
  facet_wrap(~Common_Name)

# Try using the median MaxJd as an indicator, and set the level to 180
# Create a new column representing this
passage_grid_data <- passage_grid_data %>%
  mutate(passing_grid = if_else(median_MaxJd < 180, TRUE, FALSE))
# How many passing grids are there per species?
passage_grid_data_count <- passage_grid_data %>%
  group_by(Common_Name, passing_grid) %>%
  count(name = "n_passage_grids_per_spp") %>%
  pivot_wider(names_from = passing_grid, names_prefix = "n_passage_grids_per_spp_", values_from = n_passage_grids_per_spp) %>%
  ungroup()

# Plot the summary of grids in each category for each species
ggplot(passage_grid_data, aes(passing_grid)) +
  geom_bar(colour = "black") +
  theme_bw() +
  facet_wrap(~Common_Name)

# Join the passage grid data to the main data.frame
eBird_Filtered <- left_join(eBird_Filtered, select(passage_grid_data, Common_Name, GridID, median_MaxJd, passing_grid)) %>% left_join(passage_grid_data_count)
# Filter the dataset
# Filter for if the passing grids are true, and if there are at least 5 passing grids (or staying grids) that are true for that species.
eBird_Filtered_passing_TRUE <- eBird_Filtered %>%
  filter(n_passage_grids_per_spp_TRUE >= 5 & passing_grid == TRUE)
eBird_Filtered_passing_FALSE <- eBird_Filtered %>%
  filter(n_passage_grids_per_spp_FALSE >= 5 & passing_grid == FALSE)


# Get passing-grid-specific shifts
# Create functions
coefficients_PT <- function (commonname){
  eBird_Filtered_passing_TRUE <- eBird_Filtered_passing_TRUE[eBird_Filtered_passing_TRUE$Common_Name == commonname ,]
  regression <- lmer(MAD ~ Year + (1|GridID), data = eBird_Filtered_passing_TRUE)
  MAD.Shift <- data.frame(fixef(regression))
  Common.name <- unique(eBird_Filtered_passing_TRUE$Common_Name)
  data.frame(Common.name,MAD.Shift[2,])
}
coefficients_PF <- function (commonname){
  eBird_Filtered_passing_FALSE <- eBird_Filtered_passing_FALSE[eBird_Filtered_passing_FALSE$Common_Name == commonname ,]
  regression <- lmer(MAD~Year + (1|GridID), data = eBird_Filtered_passing_FALSE)
  MAD.Shift <- data.frame(fixef(regression))
  Common.name <- unique(eBird_Filtered_passing_FALSE$Common_Name)
  data.frame(Common.name,MAD.Shift[2,])
}

#Create dataframe of slopes == MAD Shifts
# Passing == TRUE (PT)
xPT <- unique(eBird_Filtered_passing_TRUE$Common_Name)
Bird_Shifts_PT <- lapply(xPT, coefficients_PT)
Bird_Shifts_PT <- data.frame(do.call(rbind, Bird_Shifts_PT))
colnames(Bird_Shifts_PT) <- c("Common_Name", "Passing grids")
# Passing == FALSE (PF)
xPF <- unique(eBird_Filtered_passing_FALSE$Common_Name)
Bird_Shifts_PF <- lapply(xPF, coefficients_PF)
Bird_Shifts_PF <- data.frame(do.call(rbind, Bird_Shifts_PF))
colnames(Bird_Shifts_PF) <- c("Common_Name", "Staying grids")

# Join together
Bird_Shifts_PTPF <- left_join(Bird_Shifts, Bird_Shifts_PF) %>%
  left_join(Bird_Shifts_PT)
Bird_Shifts_PTPF <- Bird_Shifts_PTPF %>%
  rename("All grids" = "MAD_Shift")
Bird_Shifts_PTPF <- Bird_Shifts_PTPF %>% 
  select(Common_Name, `All grids`, `Passing grids`, `Staying grids`)
# Plot the correlations
Bird_Shifts_PTPF %>% 
  count(!is.na(`All grids`))
Bird_Shifts_PTPF %>% 
  count(!is.na(`Staying grids`))
Bird_Shifts_PTPF %>% 
  count(!is.na(`Passing grids`))
# 
Bird_Shifts_PTPF %>% 
  filter(!is.na(`Staying grids`) & !is.na(`Passing grids`)) %>% 
  summarize(mean_All = mean(`All grids`),
            mean_Passing = mean(`Passing grids`),
            mean_Staying = mean(`Staying grids`))
# 
Bird_Shifts_PTPF %>% 
  filter(!is.na(`Staying grids`) & !is.na(`Passing grids`)) %>% 
  mutate(`Passing grids abs` = abs(`Passing grids`),
         `Staying grids abs` = abs(`Staying grids`)) %>% 
  mutate(staying_greater_mag_passing = ifelse(`Passing grids abs` < `Staying grids abs`, TRUE, FALSE)) %>% 
  count(staying_greater_mag_passing)


lowerfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_abline(linetype = "dashed") +
    geom_smooth(method = "lm", colour = "black", size = 0.6) +
    geom_point(alpha = 0.5)
}


ggpairs(Bird_Shifts_PTPF[2:4],
        lower = list(continuous = wrap(lowerfun)))  +
  theme_bw() +
  labs(x = "Mean arrival date shift (days/year)",
       y = "Mean arrival date shift (days/year)") +
  theme(panel.spacing.x = unit(4, "mm"))

ggsave(filename = "Figure_S6.png", units = "cm", dpi = 300, height = 14, width = 17)

Bird_Shifts_PTPF_noNA <- na.omit(Bird_Shifts_PTPF)
nrow(Bird_Shifts_PTPF_noNA)
Bird_Shifts_PTPF_noNA$`Staying grids` < Bird_Shifts_PTPF_noNA$`Passing grids`
mean(Bird_Shifts_PTPF_noNA$`Staying grids` )
mean(Bird_Shifts_PTPF_noNA$`Passing grids` )

# Plot
ggplot(Bird_Shifts_PTPF_noNA, aes(`Staying grids`, `Passing grids`)) +
  geom_abline(linetype = "dashed") +
  geom_smooth(method = "lm", colour = "black", size = 0.6) +
  geom_point(aes(colour = Common_Name), alpha = 0.5) +
  theme_bw()

Bird_Shifts_PTPF_noNA %>%
  count(`Staying grids` < `Passing grids`)

# Early vs. late site sighting sum
eBird_Filtered_sampling_intensity_Early <- eBird_Filtered %>%
  filter(Year %in% c(2002:2004)) %>%
  select(GridID, "Site Sum mean" = "Site_Sum", "Sighting Sum mean" = "Sighting_Sum", "Sighting sum / Site sum (%) mean" = "Percent_Site_Sighting_Sums") %>%
  group_by(GridID) %>%
  summarise(`Site Sum mean` = mean(`Site Sum mean`),
            `Sighting Sum mean` = mean(`Sighting Sum mean`),
            `Sighting sum / Site sum (%) mean` = mean(`Sighting sum / Site sum (%) mean`)) %>%
  mutate(Period = "Early") %>%
  ungroup()
  # select(GridID, "Early_Site_Sum" = "Site_Sum", "Early_Sighting_Sum" = "Sighting_Sum", "Early_Percent_Site_Sighting_Sums" = "Percent_Site_Sighting_Sums")
eBird_Filtered_sampling_intensity_Late <- eBird_Filtered %>%
  filter(Year %in% c(2017:2019)) %>%
  select(GridID, "Site Sum mean" = "Site_Sum", "Sighting Sum mean" = "Sighting_Sum", "Sighting sum / Site sum (%) mean" = "Percent_Site_Sighting_Sums") %>%
  group_by(GridID) %>%
  summarise(`Site Sum mean` = mean(`Site Sum mean`),
            `Sighting Sum mean` = mean(`Sighting Sum mean`),
            `Sighting sum / Site sum (%) mean` = mean(`Sighting sum / Site sum (%) mean`)) %>%
  mutate(Period = "Late") %>%
  ungroup()
  # select(GridID, "Late_Site_Sum" = "Site_Sum", "Late_Sighting_Sum" = "Sighting_Sum", "Late_Percent_Site_Sighting_Sums" = "Percent_Site_Sighting_Sums")
eBird_Filtered_sampling_intensity <- rbind(eBird_Filtered_sampling_intensity_Early, eBird_Filtered_sampling_intensity_Late)


Grid_Shifts_lmer_join <- left_join(Grid_Shifts_lmer, eBird_Filtered_sampling_intensity)

Grid_Shifts_lmer_join_long <- Grid_Shifts_lmer_join %>%
  pivot_longer(cols = c(`Site Sum mean`, `Sighting Sum mean`, `Sighting sum / Site sum (%) mean`), names_to = "Stat", values_to = "Value")


# PLOT
Grid_Shifts_lmer_join_long %>%
  filter(!is.na(Period)) %>%
ggplot(aes(x = Value, y = MAD_Shift)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  theme_bw() +
  scale_y_continuous(n.breaks = 10) +
  labs(x = "Sampling intensity value",
       y = "Mean arrival date shift (days/year)") +
  facet_grid(rows = vars(Period), cols = vars(Stat), scales = "free_x")
