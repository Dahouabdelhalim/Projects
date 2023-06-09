library(tidyverse)
library(lme4)
library(lmerTest)
library(gridExtra)
setwd("~/Dropbox/EBM07 Shaker Study")

acti_data <- read_csv("activity_data.csv")
veggie_relish <- read_csv("veggie_data.csv")

activity8 <- acti_data %>% filter(obs<9) %>% group_by(subj_id, Day) %>% 
  summarize(total=sum(Activity), active_bins = sum(Activity > 0), age = min(age))
activity80 <- acti_data %>% filter(obs<81) %>% group_by(subj_id, Day) %>% 
  summarize(total=sum(Activity), active_bins = sum(Activity > 0), age = min(age))

lmer(active_bins ~ age * Day + (1|subj_id), data = activity8) %>% anova()
lmer(total ~ age * Day + (1|subj_id), data = activity8) %>% anova()

lmer(active_bins ~ age * Day + (1|subj_id), data = activity80) %>% anova()
lmer(total ~ age * Day + (1|subj_id), data = activity80) %>% anova()

# First 2 minutes: how many of 8 bins have any activity (Figure 3A)

day_names <- list("1" = "Day 1", "2" = "Day 2")
day_labeller <- function(variable, value) {return(day_names[value])}

lmer(active_bins ~ age * Day + (1|subj_id), data = activity8) %>% anova()
(plot3a <- ggplot(activity8, aes(x=age, y=active_bins, color=factor(Day))) + 
  geom_jitter(width=0.1, height=0.1, size = 3, show.legend = FALSE) + 
  geom_smooth(method="lm", show.legend = FALSE) + 
  facet_grid(. ~ Day, labeller = day_labeller) + scale_color_manual(values=c("navy", "darkgoldenrod1")) + 
  xlab("Age (years)") + ylab("Number of 15-sec bins with any activity in first 2 min") +
  theme_bw() + theme(strip.background=element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=20)) +
  labs(subtitle = "A") + theme(plot.subtitle=element_text(size=27)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")))

# First 2 minutes: total activity counts across 2 minutes (Figure 3B)

lmer(total ~ age * Day + (1|subj_id), data = activity8) %>% anova()
(plot3b <- ggplot(activity8, aes(x=age, y=total, color=factor(Day))) + 
  geom_jitter(width=0.1, height=0.1, size = 3, show.legend = FALSE) + 
  geom_smooth(method="lm", show.legend = FALSE) + 
  facet_grid(. ~ Day, labeller = day_labeller) + scale_color_manual(values=c("navy", "darkgoldenrod1")) + 
  xlab("Age (years)") + ylab("Total activity counts across 2 minutes") +
  theme_bw() + theme(strip.background=element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=20)) +
  labs(subtitle = "B") + theme(plot.subtitle=element_text(size=27)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")))

# All 20 minutes: how many of 80 bins have any activity (Figure 3C)

lmer(active_bins ~ age * Day + (1|subj_id), data = activity80) %>% anova()
(plot3c <- ggplot(activity80, aes(x=age, y=active_bins, color=factor(Day))) + 
    geom_jitter(width=0.1, height=0.1, size = 3, show.legend = FALSE) + 
    geom_smooth(method="lm", show.legend = FALSE) + 
    facet_grid(. ~ Day, labeller = day_labeller) + scale_color_manual(values=c("navy", "darkgoldenrod1")) + 
    xlab("Age (years)") + ylab("Number of 15-sec bins with any activity in entire 20 min") +
    theme_bw() + theme(strip.background=element_rect(fill="white")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(text = element_text(size=20)) +
    labs(subtitle = "C") + theme(plot.subtitle=element_text(size=27)) +
    theme(plot.margin = unit(c(1,1,1,1), "cm")))

# All 20 minutes: total activity counts across 20 minutes

lmer(total ~ age * Day + (1|subj_id), data = activity80) %>% anova()
(plot3d <- ggplot(activity80, aes(x=age, y=total, color=factor(Day))) + 
    geom_jitter(width=0.1, height=0.1, size=3, show.legend = FALSE) + 
    geom_smooth(method="lm", show.legend = FALSE) + 
    facet_grid(. ~ Day, labeller = day_labeller) + scale_color_manual(values=c("navy", "darkgoldenrod1")) + 
    xlab("Age (years)") + ylab("Total activity counts across 20 minutes") +
    theme_bw() + theme(strip.background=element_rect(fill="white")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(text = element_text(size=20)) +
    labs(subtitle = "D") + theme(plot.subtitle=element_text(size=27)) +
    theme(plot.margin = unit(c(1,1,1,1), "cm")))

png("Figure3.png", width = 1600, height = 1419)
grid.arrange(plot3a, plot3b, plot3c, plot3d, nrow = 2)
dev.off()

# Veggie relish: how much is lost after test period

lmer(loss ~ age * day + (1|subj_id), data = veggie_relish) %>% anova()
png("Figure4.png", width = 800, height = 709)
ggplot(veggie_relish, aes(x=age, y=loss, color=factor(day))) + 
  geom_jitter(width=0.25, height=0.25, size=3, show.legend = FALSE) + 
  geom_smooth(method="lm", show.legend = FALSE) + 
  facet_grid(. ~ day, labeller = day_labeller) + scale_color_manual(values=c("navy", "darkgoldenrod1")) + 
  xlab("Age (years)") + ylab("Veggie relish lost during test period (g)") +
  theme_bw() + theme(strip.background=element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=20))
dev.off()

# efficiency measure

efficient <- left_join(activity80, veggie_relish, by = c("subj_id" = "subj_id", "Day" = "day", "age" = "age")) %>% 
  mutate(efficiency = loss/active_bins) %>% select(subj_id, Day, total, active_bins, age, loss, efficiency) %>% 
  filter(!is.na(efficiency) & efficiency < 2.5)

lmer(efficiency ~ age * Day + (1|subj_id), data = efficient) %>% anova()
ggplot(efficient, aes(x=age, y=efficiency, color=factor(Day))) + 
  geom_jitter(width=0.25, height=0.25) + geom_smooth(method="lm") + 
  facet_grid(. ~ factor(Day)) + scale_color_brewer(palette="Set1", direction = -1) + 
  guides(color=guide_legend(title="Day")) + 
  ggtitle("Veggie relish lost per active bin") + 
  theme(plot.title = element_text(hjust = 0.5))

#outliers in activity80
activity80 %>% group_by(Day) %>% summarize(max(total))
activity80_no_outliers <- activity80 %>% filter(!(Day == 1 & total == 36939)) %>% 
  filter(!(Day == 2 & total == 53144))

lmer(active_bins ~ age * Day + (1|subj_id), data = activity80_no_outliers) %>% anova()
ggplot(activity80_no_outliers, aes(x=age, y=active_bins, color=factor(Day))) + 
  geom_jitter(width=0.1, height=0.1) + geom_smooth(method="lm") + 
  facet_grid(. ~ factor(Day)) + scale_color_brewer(palette="Set1", direction = -1) + 
  guides(color=guide_legend(title="Day")) + 
  ggtitle("How many of 80 bins over 20 minutes have any activity") + 
  theme(plot.title = element_text(hjust = 0.5))

lmer(total ~ age * Day + (1|subj_id), data = activity80_no_outliers) %>% anova()
ggplot(activity80_no_outliers, aes(x=age, y=total, color=factor(Day))) + 
  geom_jitter(width=0.1, height=0.1) + geom_smooth(method="lm") + 
  facet_grid(. ~ factor(Day)) + scale_color_brewer(palette="Set1", direction = -1) + 
  guides(color=guide_legend(title="Day")) + 
  ggtitle("Total activity counts across 20 minutes") + 
  theme(plot.title = element_text(hjust = 0.5))

# Figure 2: age distribution

ag <- activity8 %>% filter(Day == 2)
table(ag$age)
ggplot(ag, aes(x=age)) + geom_bar() + 
  xlab("Age (years)") + ylab("Number of monkeys per age") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,15)) + theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black"))

#monkeys over age 19 unable to retrieve reward from tube?

reward <- veggie_relish %>% mutate(retrieved = (loss>0), over19 = (age>19)) %>% 
  group_by(subj_id) %>% summarize(successes = sum(retrieved), age=min(age)) %>% 
  mutate(success = (successes>0), over19 = (age>19))
table(reward$success, reward$over19)
chisq.test(table(reward$success, reward$over19))

#monkeys over 24 do not manipulate tube?
#latency to first nonzero activity bin

latency <- acti_data %>% filter(Activity>0) %>% group_by(subj_id, Day) %>% summarize(firstactivebin = min(obs))
latency$firstactivebin[latency$firstactivebin > 80] <- 80
latency <- left_join(latency, activity80, by = c("subj_id", "Day")) %>% select(-active_bins)
latency$age[latency$firstactivebin == 80]

latency$age[latency$firstactivebin > 8 & latency$Day == 1]
mean(latency$age[latency$firstactivebin > 8 & latency$Day == 1])
length(latency$age[latency$firstactivebin > 8 & latency$Day == 1])
sort(latency$age[latency$firstactivebin > 8 & latency$Day == 1])

latency$age[latency$firstactivebin > 8 & latency$Day == 2]
mean(latency$age[latency$firstactivebin > 8 & latency$Day == 2])
length(latency$age[latency$firstactivebin > 8 & latency$Day == 2])
