# Dominance analyses

library(dplyr)
library(tidyr)
library(scales)

dat = read.csv("Ranks_CogScores.csv")

# Pre-analysis data massaging
# Get scaled ranks within each group 
dat = dat %>%
  group_by(season, group) %>%
  mutate(rank_scaled = rescale(rank, to=c(0,1))) %>%
# Remove birds that didn't complete enough cognition trials  
# Or that don't have cognition data at all
  filter(InitTrials > 19, RevTrials > 19) 

# Group birds into quartiles by rank
dat = dat %>%
  group_by(season, elevation, group) %>%
  arrange(season,elevation, group, rank_scaled) %>%
  mutate(rank_Q = ntile(rank_scaled, 4)) %>%
  mutate(rank_Q = as.factor(rank_Q))

# Reshape data in preparation for learning curve analysis 

# Initial testing scores
dat_Ierr = dat %>%
  ungroup() %>%
  select(tag, season, group, elevation, rank_scaled, rank_Q, InitTrials, RevTrials,
         InitScore_20, InitScore_10, InitScore_5, InitScore_3) %>%
  rename(Trials.03 = InitScore_3,
         Trials.05 = InitScore_5,
         Trials.10 = InitScore_10,
         Trials.20 = InitScore_20) %>%
  pivot_longer(cols = !c(tag, season, group, elevation, rank_scaled, rank_Q, InitTrials, RevTrials),
               names_to = "NumTrials",
               values_to = "MeanInitScore")

# Reversal testing scores
dat_Rerr = dat %>%
  ungroup() %>%
  select(tag, season, elevation,
         RevScore_20, RevScore_10, RevScore_5, RevScore_3) %>%
  rename(Trials.03 = RevScore_3,
         Trials.05 = RevScore_5,
         Trials.10 = RevScore_10,
         Trials.20 = RevScore_20) %>%
  pivot_longer(cols = !c(tag, season, elevation),
               names_to = "NumTrials",
               values_to = "MeanRevScore")

# Combine
dat_Qlong = full_join(dat_Ierr, dat_Rerr, by = c("tag", "season", "elevation", "NumTrials"))
rm(dat_Ierr, dat_Rerr)


#===========================================================================================================

# Regression - 20 trials score vs dominance
library(car)
library(ggfortify)

# Split data by seasons 
season1920 = subset(dat, season == "S19-20")
season2021 = subset(dat, season == "S20-21")

#--- 19-20 ------------------------------------------------------
# Initial
mI20_1920 = lm(InitScore_20 ~ elevation*rank_scaled, data = season1920)
Anova(mI20_1920, type = "III")
summary(mI20_1920)
autoplot(mI20_1920)

# Check model without outlier
mI20_1920_o = lm(InitScore_20 ~ elevation*rank_scaled, data = season1920[-95,])
Anova(mI20_1920_o, type = "III")
summary(mI20_1920_o)

# Reversal
mR20_1920 = lm(RevScore_20 ~ elevation*rank_scaled, data = season1920)
Anova(mR20_1920, type = "III")
summary(mR20_1920)
autoplot(mR20_1920)

# Check model without outlier
mR20_1920_o = lm(RevScore_20 ~ elevation*rank_scaled, data = season1920[-126,])
Anova(mR20_1920_o, type = "III")
summary(mR20_1920_o )

#--- 20-21 -----------------------------
# Initial
mI20_2021 = lm(InitScore_20 ~ elevation*rank_scaled, data = season2021)
Anova(mI20_2021, type = "III")
summary(mI20_2021)
autoplot(mI20_2021)

# Reversal
mR20_2021 = lm(RevScore_20 ~ elevation*rank_scaled, data = season2021)
Anova(mR20_2021, type = "III")
summary(mR20_2021)
autoplot(mR20_2021)

#=================================================================================================================
# Learning curves 
library(lme4)
library(emmeans)

season1920_err = subset(dat_Qlong, season == "S19-20")
season2021_err = subset(dat_Qlong, season == "S20-21")

# ---19-20 ---
# Initial  
mLearn1920_3w = lmer(MeanInitScore ~ (elevation*rank_Q*NumTrials) + (1|tag), data = season1920_err)
Anova(mLearn1920_3w, type = "III", test.statistic = "F")

mLearn1920 = lmer(MeanInitScore ~ (elevation + rank_Q + NumTrials)^2 + (1|tag), data = season1920_err)
Anova(mLearn1920, type = "III", test.statistic = "F")

anova(mLearn1920_3w, mLearn1920)

emmeans(mLearn1920, pairwise ~ elevation*rank_Q | NumTrials)
emmeans(mLearn1920, pairwise ~ rank_Q | elevation | NumTrials)

# Reversal
mLearnR1920_3w = lmer(MeanRevScore ~ (elevation*rank_Q*NumTrials) + (1|tag), data = season1920_err)
Anova(mLearnR1920_3w, type = "III", test.statistic = "F")

mLearnR1920 = lmer(MeanRevScore ~ (elevation + rank_Q + NumTrials)^2 + (1|tag), data = season1920_err)
Anova(mLearnR1920, type = "III", test.statistic = "F")

anova(mLearnR1920_3w, mLearnR1920)

emmeans(mLearnR1920, pairwise ~ elevation*rank_Q | NumTrials)
emmeans(mLearnR1920, pairwise ~ rank_Q | elevation | NumTrials)

# ---20-21 -----------------------------
# Initial  
mLearn2021_3w = lmer(MeanInitScore ~ (elevation*rank_Q*NumTrials) + (1|tag), data = season2021_err)
Anova(mLearn2021_3w, type = "III", test.statistic = "F")

mLearn2021 = lmer(MeanInitScore ~ (elevation + rank_Q + NumTrials)^2 + (1|tag), data = season2021_err)
Anova(mLearn2021, type = "III", test.statistic = "F")

anova(mLearn2021_3w, mLearn2021)

emmeans(mLearn2021, pairwise ~ elevation*rank_Q | NumTrials)
emmeans(mLearn2021, pairwise ~ rank_Q | elevation | NumTrials)


# Reversal
mLearnR2021_3w = lmer(MeanRevScore ~ (elevation*rank_Q*NumTrials) + (1|tag), data = season2021_err)
Anova(mLearnR2021_3w, type = "III", test.statistic = "F")
emmeans(mLearnR2021_3w, pairwise ~ elevation*rank_Q | NumTrials)
emmeans(mLearnR2021_3w, pairwise ~ rank_Q  | elevation | NumTrials)

mLearnR2021 = lmer(MeanRevScore ~ (elevation + rank_Q + NumTrials)^2 + (1|tag), data = season2021_err)
Anova(mLearnR2021, type = "III", test.statistic = "F")

anova(mLearnR2021_3w, mLearnR2021)

mLearnR2021_3wH = lmer(MeanRevScore ~ (rank_Q*NumTrials) + (1|tag), data = subset(season2021_err, elevation == "H"))
Anova(mLearnR2021_3wH, type = "III", test.statistic = "F")
emmeans(mLearnR2021_3wH, pairwise ~ rank_Q  | NumTrials)

#===========================================================================================================
# Ordered alternatives test
library(PMCMRplus)
cuzickTest(season2021_err$MeanInitScore, season2021_err$rank_Q)
jonckheereTest(season2021_err$MeanInitScore, season2021_err$rank_Q, alternative = "less")

# Make sure levels go the right way 
rlevels = c(1,2,3,4)
season2021_err = season2021_err %>% mutate(rank_Qo=ordered(rank_Q,levels=rlevels))
season1920_err = season1920_err %>% mutate(rank_Qo=ordered(rank_Q,levels=rlevels))

season2021_err %>% filter(NumTrials == "Trials.03") %>%
  with(.,jonckheereTest(MeanInitScore, rank_Qo,alternative = "greater"))

season2021_err %>% filter(NumTrials == "Trials.05") %>%
  with(.,jonckheereTest(MeanInitScore, rank_Qo,alternative = "greater"))

season2021_err %>% filter(NumTrials == "Trials.10") %>%
  with(.,jonckheereTest(MeanInitScore, rank_Qo,alternative = "greater"))

season2021_err %>% filter(NumTrials == "Trials.20") %>%
  with(.,jonckheereTest(MeanInitScore, rank_Qo,alternative = "greater"))


season2021_err %>% filter(NumTrials == "Trials.03") %>%
  with(.,jonckheereTest(MeanRevScore, rank_Qo, alternative = "greater"))

season2021_err %>% filter(NumTrials == "Trials.05") %>%
  with(.,jonckheereTest(MeanRevScore, rank_Qo,alternative = "greater"))

season2021_err %>% filter(NumTrials == "Trials.10") %>%
  with(.,jonckheereTest(MeanRevScore, rank_Qo,alternative = "greater"))

season2021_err %>% filter(NumTrials == "Trials.20") %>%
  with(.,jonckheereTest(MeanRevScore, rank_Qo,alternative = "greater"))


#===========================================================================================================
# Figures!
# Get colors

library(ggplot2)
library(patchwork)
colors = c("#e66101", "#1D50CE")

#=================================================================================================================
# Initial

ySeq1 = (seq(0, max(season2021$InitScore_20), by = 0.5))
ySeq2 = (seq(0, 3.5, by = 0.5))

# Initial 20-21
p1 = ggplot(subset(season2021), 
            aes(x = rank_scaled,
                y = InitScore_20)) +
  geom_jitter(aes(color = elevation),
              alpha = 0.9,
              height = 0,
              width = 0.03) +
  expand_limits(y=0) +
  scale_x_continuous(limits = c(0, 1), oob = scales::squish) +
  scale_y_continuous(limits = c(0, 3), breaks = ySeq2) +
  scale_color_manual(values = colors,
                     name = "Elevation",
                     labels = c("High", "Low")) +
  xlab("Dominance rank") +
  ylab("") +
  theme_bw() +
  ggtitle("(b)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 14))  +
  theme(axis.title = element_text(size = 16)) +
  theme(legend.justification=c(0,1),
        legend.position=c(0.03, 0.97),
        legend.background = element_rect(colour = 1, size = 0.2))

# Initial 19-20
p3 = ggplot(subset(season1920), 
            aes(x = rank_scaled,
                y = InitScoreTrue_20)) +
  geom_jitter(aes(color = elevation),
              alpha = 0.9,
              height = 0,
              width = 0.03) +
  expand_limits(y=c(0)) +
  scale_x_continuous(limits = c(0, 1), oob = scales::squish) +
  scale_y_continuous(limits = c(0, 3), breaks = ySeq2) +
  #  geom_smooth(method = "lm", se = F) +
  scale_color_manual(values = colors,
                     name = "Elevation",
                     labels = c("High", "Low")) +
  xlab("Dominance rank") +
  ylab("Mean number of location errors per trial") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 14))  +
  theme(axis.title = element_text(size = 16)) +
  ggtitle("(a)") +
  theme(legend.justification=c(0,1),
        legend.position=c(0.03, 0.97),
        legend.background = element_rect(colour = 1, size = 0.2))

pInitial = p3 + p1 + plot_layout(guides = "collect")
pInitial

# ---- Reversal -----------------------

# Reversal 20-21
p2 = ggplot(subset(season2021), 
            aes(x = rank_scaled,
                y = RevScore_20)) +
  geom_jitter(aes(color = elevation),
              alpha = 0.9,
              height = 0,
              width = 0.03) +
  expand_limits(y=0) +
  scale_x_continuous(limits = c(0, 1), oob = scales::squish) +
  scale_y_continuous(limits = c(0, 1.5), breaks = ySeq2) +
  stat_smooth(method = "lm", se = T,
              level = 0.95,
              color = "black") +
  scale_color_manual(values = colors,
                     name = "Elevation",
                     labels = c("High", "Low")) +
  xlab("Dominance rank") +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 14))  +
  theme(axis.title = element_text(size = 16)) +
  ggtitle("(b)") +
  theme(legend.justification=c(0,1),
        legend.position=c(0.03, 0.97),
        legend.background = element_rect(colour = 1, size = 0.2))


# Reversal 19-20
p4 = ggplot(subset(season1920), 
            aes(x = rank_scaled,
                y = RevScore_20)) +
  geom_jitter(aes(color = elevation),
              alpha = 0.9,
              height = 0,
              width = 0.03) +
  expand_limits(y=0) +
  scale_x_continuous(limits = c(0, 1), oob = scales::squish) +
  scale_y_continuous(limits = c(0, 1.5), breaks = ySeq2) +
  scale_color_manual(values = colors,
                     name = "Elevation",
                     labels = c("High", "Low")) +
  xlab("Dominance rank") +
  ylab("Mean number of location errors per trial") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 14))  +
  theme(axis.title = element_text(size = 16)) +
  ggtitle("(a)") +
  theme(legend.justification=c(0,1),
        legend.position=c(0.03, 0.97),
        legend.background = element_rect(colour = 1, size = 0.2))


pReversal = p4 + p2 + plot_layout(guides = "collect")
pReversal

# --- Learning curves ----------------------------------------------------


elevation.labs = c("High", "Low")
names(elevation.labs) = c("H", "L")


task.labs = c("Spatial learning and memory task", "Reversal spatial learning task")
names(task.labs) = c("MeanInitScore", "MeanRevScore")

#------------------------------------------
season2021_err_l = season2021_err %>%
  pivot_longer(cols = c(MeanInitScore, MeanRevScore),
               names_to = "Task", values_to = "Score")

season1920_err_l = season1920_err %>%
  pivot_longer(cols = c(MeanInitScore, MeanRevScore),
               names_to = "Task", values_to = "Score")

# 20-21
p10 = ggplot(subset(season2021_err_l, !is.na(rank_Q)), aes(x = NumTrials,
                                                           y = Score,
                                                           group = rank_Q,
                                                           color = rank_Q)) +
  stat_summary(geom = "point",
               position = position_dodge(width = 0.2),
               na.rm = T,
               size = 2) +
  stat_summary(geom = "errorbar",
               width = 0,
               position = position_dodge(width = 0.2),
               na.rm = T,
               show.legend = F) +
  stat_summary(geom = "line",
               position = position_dodge(width = 0.2),
               na.rm = T,
               show.legend = F) +
  scale_color_brewer(name = "Rank\\nquartile",
                     palette = "RdBu") +
  expand_limits(y=0) +
  ylab("Mean number of location errors per trial") +
  xlab("Number of trials") +
  scale_x_discrete(labels = c("3", "5", "10", "20")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(rows = vars(elevation),
             cols = vars(Task),
             labeller = labeller(elevation = elevation.labs, Task = task.labs))

p10


# 19-20

p11 = ggplot(subset(season1920_err_l, !is.na(rank_Q)), aes(x = NumTrials,
                                                           y = Score,
                                                           group = rank_Q,
                                                           color = rank_Q)) +
  stat_summary(geom = "point",
               position = position_dodge(width = 0.2),
               na.rm = T,
               size = 2) +
  stat_summary(geom = "errorbar",
               width = 0,
               position = position_dodge(width = 0.2),
               na.rm = T,
               show.legend = F) +
  stat_summary(geom = "line",
               position = position_dodge(width = 0.2),
               na.rm = T,
               show.legend = F) +
  scale_color_brewer(name = "Rank\\nquartile",
                     palette = "RdBu") +
  expand_limits(y=0) +
  ylab("Mean number of location errors per trial") +
  xlab("Number of trials") +
  scale_x_discrete(labels = c("3", "5", "10", "20")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(rows = vars(elevation),
             cols = vars(Task),
             labeller = labeller(elevation = elevation.labs, Task = task.labs))

p11
