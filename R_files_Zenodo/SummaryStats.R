## Juvenile sociality summary statistics
## Ali Galezo

library(dplyr)
library(ggplot2)
library(tidyr)

options(stringsAsFactors = FALSE)

# Load data on all surveys with at least one juvenile
surveys <- read.csv("juvenile_surveys.csv")

# FIGURE: group-size histogram for all-juvie groups
surveys %>%
  filter(all_weaned_juvies == TRUE) %>% # all-juvie groups only
  distinct(observation_id, group_size, group_sex) %>% # one row per survey
  ggplot(aes(x = group_size)) +
  geom_bar(fill = "#333333") +
  xlab("Group Size") +
  ylab("Frequency") +
  scale_y_continuous(breaks = c(0,200,400,600,800,1000,1200), expand = c(0.01,0), limits = c(0, 1205)) +
  scale_x_continuous(breaks = 1:7) +
  theme_classic()
ggsave("GroupSizeHistogram.tif", plot = last_plot(), device = "tiff", width = 85, height = 85, units = "mm", dpi = 300)

# FIGURE: group size vs. sex composition for all-juvie groups
surveys %>%
  ungroup() %>%
  filter(all_weaned_juvies == TRUE) %>% # all-juvie groups only
  distinct(observation_id, group_size, group_sex) %>% # one row per survey
  filter(group_sex != "unknowns") %>%
  mutate(group_size = as.character(group_size),
         group_size = ifelse(group_size >= 4, "4+", group_size), # pool groups sized 4+ together
         group_sex = ifelse(group_sex == "Mixed Sex", "Mixed", group_sex), # rename 'mixed sex' to 'mixed'
         group_sex = factor(group_sex, levels = c("All Female", "Mixed", "All Male"))) %>% # set order of factor levels
  count(group_size, group_sex) %>%
  ggplot(aes(x = group_size, y = n, fill = group_sex)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  xlab("Group Size") +
  ylab("Proportion of Groups") +
  guides(fill = guide_legend(title = "Group Type")) +
  theme_classic() +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0,0)) +
  scale_fill_grey(start = 0, end = 1)
ggsave("GroupSizexSex.tif", plot = last_plot(), device = "tiff", width = 85, height = 80, units = "mm", dpi = 300)

# Number of survey sightings with all weaned juvies
surveys %>%
  ungroup() %>%
  filter(all_weaned_juvies == TRUE) %>%
  summarize(total = n_distinct(observation_id))

# Proportion time alone
alone <- surveys %>%
  group_by(dolphin_id, sex, weaned_juvie) %>%
  count(group_size) %>%
  filter(sex %in% c("FEMALE","MALE")) %>% # target individual must have a known sex
  filter(weaned_juvie == TRUE) %>% # only include sightings when target individual was a juvie
  mutate(total = sum(n)) %>%
  filter(total >= 15) %>%  # animals appearing in at least 15 surveys
  mutate(proportion = n/total) %>%
  ungroup() %>%
  complete(nesting(dolphin_id, sex, weaned_juvie, total), group_size, fill = list(n = 0, proportion = 0)) %>%
  filter(group_size == 1)
alone %>%
  group_by(sex) %>%
  summarize(mean = mean(proportion),
            n = n(),
            sd = sd(proportion),
            se = sd / sqrt(n))
perm::permTS(x = alone %>% filter(sex == "MALE") %>% pull(proportion),
             y = alone %>% filter(sex == "FEMALE") %>% pull(proportion))

# Mean group size
# in any groups:
any_group <- surveys %>%
  group_by(dolphin_id, sex, weaned_juvie) %>%
  summarize(min = min(group_size),
            max = max(group_size),
            mean = mean(group_size),
            sd = sd(group_size),
            n = n(),
            se = sd / sqrt(n)) %>%
  filter(sex %in% c("FEMALE","MALE")) %>% # target individual must have a known sex
  filter(weaned_juvie == TRUE) %>% # only include sightings when target individual was a juvie
  filter(n >= 15)  # animals appearing in at least 15 surveys during juvenile period
any_group %>%
  group_by(sex) %>%
  summarize(min_mean = min(mean),
            max_mean = max(mean),
            mean_mean = mean(mean),
            n = n(),
            sd = sd(mean),
            se = sd / sqrt(n))
perm::permTS(x = any_group %>% filter(sex == "FEMALE") %>% pull(mean),
             y = any_group %>% filter(sex == "MALE") %>% pull(mean))
# in all-juvie groups:
juvie_group <- surveys %>%
  filter(all_weaned_juvies == TRUE) %>%
  group_by(dolphin_id, sex) %>%
  summarize(min = min(group_size),
            max = max(group_size),
            mean = mean(group_size),
            sd = sd(group_size),
            n = n(),
            se = sd / sqrt(n)) %>%
  filter(sex %in% c("FEMALE","MALE")) %>%
  filter(n >= 15) # animals appearing in at least 15 all-juvie groups
juvie_group %>%
  group_by(sex) %>%
  summarize(min_mean = min(mean),
            max_mean = max(mean),
            mean_mean = mean(mean),
            n = n(),
            sd = sd(mean),
            se = sd / sqrt(n))
perm::permTS(x = juvie_group %>% filter(sex == "FEMALE") %>% pull(mean),
             y = juvie_group %>% filter(sex == "MALE") %>% pull(mean))
