library(lme4)
library(ggplot2)
library(ggpubr)
#Read file
transect_data <- read.csv("transect_dataset.csv")
transect_data <- subset(transect_data, location != "Unknown"); transect_data$location = factor(transect_data$location)  # don't use unknown location
str(transect_data)
summary(transect_data)

singing_birds <- subset(transect_data, song == 's')
# we had a total of 116 transects (not taken from this data), during 116-39= 77 of which we observed zebra finches a total of 265 times, of which 92 observations included singing birds

mean(singing_birds$height.m., na.rm = TRUE)
sd(singing_birds$height.m., na.rm = TRUE) / sqrt(sum(!is.na(singing_birds$height.m.)))  # SE

p_height <- ggplot(singing_birds, aes(x=song, y=height.m.)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis='y', stackdir='center', alpha = 1, binwidth = 0.1, dotsize = 0.55) +
  labs(y = 'Perch height (m)', x = '') +
  scale_y_continuous(expand=c(0,0), breaks = seq(0, 4, 1), limits = c(0,3.5)) +
  scale_x_discrete(breaks=c('')) +
  theme_classic() +
  theme(text = element_text(size = 12))
p_height
# Of 49 song observations we scored the height of the singing individual. In the environment of our field site, which is mostly dominated by low shrubs and trees), the mean height of singing birds was 1.6m (range was 0.3m-3m).

#distance to other individuals when singing
# of 43 song observations we scored the distance to other individuals
mean(transect_data$within_distance.m., na.rm = TRUE)
sd(transect_data$within_distance.m., na.rm = TRUE) / sqrt(sum(!is.na(transect_data$within_distance.m.)))  # SE

p_birds <- ggplot(singing_birds, aes(x=song, y=within_distance.m.)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis='y', stackdir='center', alpha = 1, binwidth = 0.1) +
  labs(y = 'Max distance between group members (m)', x = '') +
  scale_y_continuous(expand=c(0,0), breaks = seq(0, 6, 1), limits = c(0,6.5)) +
  scale_x_discrete(breaks=c('')) +
  theme_classic() +
  theme(text = element_text(size = 12))
p_birds

# secondary group distance - more anecdotal
mean(transect_data$distance_second.m., na.rm = TRUE)
sd(transect_data$distance_second.m., na.rm = TRUE) / sqrt(sum(!is.na(transect_data$distance_second.m.)))  # SE

p_secondary <- ggplot(singing_birds, aes(x=song, y=distance_second.m.)) +
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(binaxis='y', stackdir='center',alpha = 1, binwidth = 2, dotsize = 0.3) +
  labs(y = 'Distance between groups (m)', x = '') +
  scale_y_continuous(expand=c(0,0), breaks = seq(0, 40, 5), limits = c(0,40)) +
  scale_x_discrete(breaks=c('')) +
  theme_classic() +
  theme(text = element_text(size = 12))
p_secondary

ggarrange(
  p_height, p_birds, p_secondary, 
  ncol=3, nrow=1, common.legend=FALSE, labels='auto', font.label = list(size = 12, face = "plain"))
ggsave('figure3.tiff', width = 129, height = 100, units = 'mm', dpi = 1200, compression = "lzw")
