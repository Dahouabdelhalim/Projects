# Visualisation + Data Analysis of Transmission experiment conducted in Fowlers Gap in 2018
# This is for the dataset where all vocalisations are bandpassed with bins.
# By Hugo Loning, September 2020

library(ggplot2)
library(lme4)
library(ggpubr)
library(emmeans)

# Load + prepare data
oka_dool3 <- read.csv('Fig3_Okanoya_Dooling_1987.csv')  # absolute threshold of zebra finches (audibility curve)

df <- read.csv('processed_bandpassed_vocalisations.csv', stringsAsFactors = TRUE)
summary(df)
str(df)
df <- subset(df, quality == 'include')  # remove all observations that have too much noise
df$transect <- as.factor(df$transect)
df$repeat. <- as.factor(df$repeat.)
df$kfreq <- df$freq/1000

# Stats

# Amplitude
full_model <- lmer(signal_dB~log2(dist)*kfreq*noise_dB+voc_type+(1|transect), data = df, REML = FALSE)
summary(full_model)

no_3_way <- update(full_model, ~.-log2(dist):kfreq:noise_dB)
summary(no_3_way)
anova(full_model, no_3_way) # so the model is better with 3 way interaction included, so need to keep that one in, and also the 2 way interactions between freq and doubling of dist, noise and dist and noise and frequency

emmeans(full_model, list(pairwise ~ voc_type), adjust = "tukey")

no_voc <- lmer(signal_dB~log2(dist)*kfreq*noise_dB+(1|transect), data = df, REML = FALSE)
summary(no_voc)
anova(full_model, no_voc)  # There is an effect of vocalisation

no_freq <- lmer(signal_dB~log2(dist)*noise_dB+voc_type+(1|transect), data = df, REML = FALSE)
plot(no_freq)
anova(full_model, no_freq)  # of course freq is significant

no_dist <- lmer(signal_dB~kfreq*noise_dB+voc_type+(1|transect), data = df, REML = FALSE)
plot(no_dist)
anova(full_model, no_dist)  # of course dist is significant


# Smoothed line plots
muted = rev(c('#332288','#88ccee','#44aa99','#117733','#999933','#ddcc77','#cc6677','#882255','#aa4499'))

p_song_ampl <- ggplot() + 
  geom_point(data = subset(df, voc_type == 'song'),
             aes(x = kfreq, y = signal_dB, group=dist, colour = as.factor(dist)),
             alpha = 0.05, size = 0.2) +
  geom_smooth(method='gam', span = 0.2, se = FALSE, 
              data = subset(df, voc_type == 'song'),
              aes(x = kfreq, y = signal_dB, group=dist, colour = as.factor(dist)), size = 0.7) +
  geom_smooth(se = FALSE, 
              data = df, aes(x = kfreq, y = noise_dB), colour = 'black', linetype = 3, size = 0.7) +
  geom_line(data = oka_dool3, aes(x = Freq_Hz/1000, y = SPL_dB), size = 0.7) +
  labs(x = '', y = 'Amplitude (dB)', colour = 'Distance (m)', subtitle = 'Song') + 
  scale_x_continuous(breaks = seq(0,8,1)) +
  scale_y_continuous(breaks = seq(-10, 50, 10)) +
  coord_cartesian(xlim = c(0,8), ylim = c(-8,58), expand = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 12), plot.subtitle = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = muted)

p_dcm_ampl <- ggplot() + 
  geom_point(data = subset(df, voc_type == 'dc_male'), 
             aes(x = kfreq, y = signal_dB, group=dist, colour = as.factor(dist)), 
             alpha = 0.05, size = 0.2) +
  geom_smooth(method='gam', span = 0.2, se = FALSE, 
              data = subset(df, voc_type == 'dc_male'), 
              aes(x = kfreq, y = signal_dB, group=dist, colour = as.factor(dist)), size = 0.7) +
  geom_smooth(se = FALSE, 
              data = df, aes(x = kfreq, y = noise_dB), colour = 'black', linetype = 3, size = 0.7) +
  geom_line(data = oka_dool3, aes(x = Freq_Hz/1000, y = SPL_dB), size = 0.7) +
  labs(x = 'Frequency (kHz)', y = '', colour = 'Distance (m)', subtitle = 'Male calls') + 
  scale_x_continuous(breaks = seq(0,8,1)) +
  scale_y_continuous(breaks = seq(-10, 50, 10)) +
  coord_cartesian(xlim = c(0,8), ylim = c(-8,58), expand = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 12), plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y =  element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = muted)

p_dcf_ampl <- ggplot() + 
  geom_point(data = subset(df, voc_type == 'dc_female'), 
             aes(x = kfreq, y = signal_dB, group=dist, colour = as.factor(dist)), 
             alpha = 0.05, size = 0.2) +
  geom_smooth(method='gam', span = 0.2, se = FALSE, 
              data = subset(df, voc_type == 'dc_female'), 
              aes(x = kfreq, y = signal_dB, group=dist, colour = as.factor(dist)), size = 0.7) + 
  geom_smooth(se = FALSE, 
              data = df, aes(x = kfreq, y = noise_dB), colour = 'black', linetype = 3, size = 0.7) +
  geom_line(data = oka_dool3, aes(x = Freq_Hz/1000, y = SPL_dB), size = 0.7) +
  labs(x = '', y = '', colour = 'Distance (m)', subtitle = 'Female calls') + 
  scale_x_continuous(breaks = seq(0,8,1)) +
  scale_y_continuous(breaks = seq(-10, 50, 10)) +
  coord_cartesian(xlim = c(0,8), ylim = c(-8,58), expand = FALSE) +
  theme_classic() +
  theme(text = element_text(size = 12), plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y =  element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_color_manual(values = muted)

ggarrange(
  p_song_ampl, p_dcm_ampl, p_dcf_ampl, 
  ncol=3, nrow=1, common.legend=TRUE, legend="right", labels='auto', font.label = list(size = 12, face = "plain"), widths = c(1.25,1,1), label.x = c(0.167,0,0))
ggsave('figure2.tiff', width = 177, height = 120, units = 'mm', dpi = 1200, compression = "lzw")
