# Data analysis and figure plotting script for zebra finch song ecology data 
# By Hugo Loning, 2022

library(tidyverse)
library(patchwork)
library(glmmTMB)
library(ggpubr)

##### GANTT CHART 2016-2020 #####

data_info <- read_csv("breeding_over_the_years_2016-2020.csv")

{
  tick_marks <- data.frame(xa = c(12.5,12.5), ya = c(-100,0))
  
plot_rainfall_breeding <- ggplot(data = data_info, mapping = aes(x = month_index)) +
    geom_bar(mapping = aes(y = monthly_rainfall_mm/150, fill = "Rainfall"), stat = "identity", color = "black", size = 0.3, alpha = 0.8) +
    geom_bar(mapping = aes(y = clutches_laid/140, fill = "Clutches"), stat = "identity", color = "black", size = 0.3, alpha = 0.75) +
    geom_line(data = tick_marks, aes(x = xa-12, y = ya), size = 0.3) +
    geom_line(data = tick_marks, aes(x = xa, y = ya), size = 0.3) +
    geom_line(data = tick_marks, aes(x = xa+12, y = ya), size = 0.3) +
    geom_line(data = tick_marks, aes(x = xa+24, y = ya), size = 0.3) +
    geom_line(data = tick_marks, aes(x = xa+36, y = ya), size = 0.3) +
    scale_x_continuous(breaks = seq(6.5, 54.5, 12), 
                       labels = c("2016", "2017", "2018", "2019", "2020")) +
    scale_fill_manual(values = c("Clutches" = "#DDAA33", "Rainfall" = "#004488"), name = "Type") +
    scale_y_continuous(breaks = seq(0, 1, 1/7), labels = c("0", "20", "40", "60", "80", "100", "120", "140"),
                       sec.axis = sec_axis(~ ., 
                                           seq(0,1,1/3), 
                                           labels = c("0", "50", "100", "150"), 
                                           name = "Monthly rainfall (mm)")) +
    coord_cartesian(xlim = c(3.25, 57.75), ylim = c(0, 1)) +
    labs(x = "Year", y = "Clutches laid per month") +
    theme_classic() +
    theme(axis.title.y.right = element_text(angle=90),
          legend.position = c(0.55, 0.85),
          #legend.background = element_rect(color = "black"),
          legend.key.height= unit(0.2, 'cm'),
          legend.key.width= unit(0.2, 'cm')); plot_rainfall_breeding
}

{
tick_marks <- data.frame(xa = c(12.5,12.5), ya = c(-100,0.05))

nest_recs2016 <- data.frame(xa = c(9.5,10.5), ya = c(0.3,0.3))
playback2017 <- data.frame(xa = c(21.5,22.5), ya = c(0.1,0.1))
focal_recs2018 <- data.frame(xa = c(32.5,34.5), ya = c(0.4,0.4))
songmeters2018 <- data.frame(xa = c(33.5,46.5), ya = c(0.2,0.2))
transects_2018 <- data.frame(xa = c(32.5,36.5), ya = c(0.5,0.5))
transects_2019 <- data.frame(xa = c(44.5,48.5), ya = c(0.5,0.5))
transects_2020 <- data.frame(xa = c(58.5,59.5), ya = c(0.5,0.5))

plot_breeding2 <- ggplot(data = data_info, mapping = aes(x = month_index)) +
  geom_bin2d(mapping = aes(y = rep.int(0.5,60), fill = clutches_laid/max(clutches_laid)), stat = "identity", alpha = 1) +
  scale_fill_gradient(low="white", high = "#DDAA33", guide = "none") +
  geom_line(data = tick_marks, aes(x = xa-12, y = ya), size = 0.3) +
  geom_line(data = tick_marks, aes(x = xa, y = ya), size = 0.3) +
  geom_line(data = tick_marks, aes(x = xa+12, y = ya), size = 0.3) +
  geom_line(data = tick_marks, aes(x = xa+24, y = ya), size = 0.3) +
  geom_line(data = tick_marks, aes(x = xa+36, y = ya), size = 0.3) +
  geom_line(data = tick_marks, aes(x = xa+48, y = ya), size = 0.3) +
  geom_line(data = tick_marks, aes(x = xa+48, y = ya), size = 0.3) +
  
  geom_line(data = nest_recs2016, aes(x = xa, y = ya), size = 3) +
  geom_line(data = playback2017, aes(x = xa, y = ya), size = 3) +
  geom_line(data = focal_recs2018, aes(x = xa, y = ya), size = 3) +
  geom_line(data = songmeters2018, aes(x = xa, y = ya), size = 3) +
  geom_line(data = transects_2018, aes(x = xa, y = ya), size = 3) +
  geom_line(data = transects_2019, aes(x = xa, y = ya), size = 3) +
  geom_line(data = transects_2020, aes(x = xa, y = ya), size = 3) +
  scale_x_continuous(breaks = seq(6.5, 54.5, 12), 
                     labels = c("2016", "2017", "2018", "2019", "2020")) +
  scale_y_continuous(breaks = seq(0.1,0.5,0.1),
                     labels = c("playback exp.", "year-round rec.", "nest rec.", "focal rec.", "transects")) +
  coord_cartesian(xlim = c(3.25, 57.9), ylim = c(0.05, 0.5)) +
  labs(x = "Year", y = "Activity") +
  theme_classic(); plot_breeding2
}


##### TRANSECT OBSERVATIONS 2018-2020 #####

### Read + tidy data
data_transect <- read_csv("transects_2018-2020.csv")
data_transect <- data_transect %>% 
  filter(location != "Unknown" & # don't use unknown location
           type != "f" & # don't use flying observations since zebra finches don't sing while flying (and there were probably some mistakes in the 2020 data in this regard by slightly inexperienced observers)
           !is.na(song) # don't use observations where it is unclear whether it is song or not
         ) %>% 
  mutate(song_binary = if_else(song == 's', 1, 0), # add numerical equivalent for song for the geom_smooth logistic regression
         breeding = if_else(year == 2020, 'Breeding', 'Non-breeding'), # add binary vector for breeding
         .after = 'song'
         ) %>% 
  mutate(across(c('location', 'type', 'song', 'breeding'), # turn all relevant chr into factor
                .fns = as.factor
                )
         ) %>% 
  mutate(breeding = factor(breeding, levels = c("Non-breeding", "Breeding")) # relevel breeding factor
         )
head(data_transect)

# Of 570 observations, there is 1 data point of a group size at 80 (not singing) while the next biggest group is 43, 
# this makes all the plots very zoomed out, so I decide to omit this point

### Effect of observer?
ggplot(data = filter(data_transect, total != 80), 
       mapping = aes(x=total, y=song_binary)) + 
  geom_count(alpha = 0.3) +
  geom_smooth(colour = "black", method = 'glm', method.args = list(family=binomial), formula = y ~ x) +
  facet_wrap(~observer, nrow = 1) +
  labs(y = 'Song presence', x = 'Group size') +
  theme_classic() +
  theme(strip.background = element_rect(colour="white"))

# I will exclude James Madden's observations since I don't completely trust his judgement of what is song and not.
# Firstly, James also sometimes scored song of flying birds (this was never the case for Daniel Kovicz), I suspect that
# he sometimes scored repetitive distance calls as song. Understandable, because the song is very repetitive and soft 
# and contains the distance call elements also as loudest element. So, to a relatively untrained ear they sound similar.
# I have never heard a zebra finch sing during flight, so I don't want to include these observations.
# Additionally, Daniel Kovicz finds the same baseline as me for song of a single individual which indeed should be 
# probably below 50% since male/female ratio is probably around 50-50 or 60-40 ish at most and in any case solo males don't 
# sing that often in my experience. The base 60% for James Madden also is striking in this regard.

summary(filter(data_transect, observer != "James Madden"))

transect_plot <- ggplot(data = filter(data_transect, total != 80 & observer != "James Madden"), 
                        mapping = aes(x=total, y=song_binary)) + 
  geom_count(alpha = 0.3) +
  geom_smooth(colour = "black", method = 'glm', method.args = list(family=binomial), formula = y ~ x) +
  labs(y = 'Song presence', x = 'Group size', size = 'Observations') +
  theme_classic() +
  theme(strip.background = element_rect(colour="white"),
        legend.position = c(0.9,0.35),
        legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.4, 'cm')
        #,legend.background = element_rect(color = "black")
        )

### Group size different between the different years?
ggplot(data = filter(data_transect, observer != "James Madden"),
       mapping = aes(x=as.factor(year), y=total)) +
  geom_boxplot()
summary(glmmTMB(total~as.factor(year), 
                family = genpois, 
                data = filter(data_transect, observer != "James Madden")))
# 2018 - 2019: z = 0.33, df = 351, P = 0.74
# 2018 - 2020: z = 1.09, df = 351, P = 0.28

# by dividing 2020 by year, I reverse the order, so that 2020 is first, then 2019, then 2018
# this is helpful for making the last comparison between 2019 and 2020, which is not possible when 2018
# is in the intercept
summary(glmmTMB(total~as.factor(2020/year), 
                family = genpois, 
                data = filter(data_transect, observer != "James Madden"))) # nope, they are similar

# 2019-2020: z = 0.73, df = 351, P = 0.47

### Effect of year?
plot_transect_year <- transect_plot +
  facet_wrap(~year, nrow = 1, labeller = labeller(year = c("2018" = "2018 (Non-breeding)", 
                                                           "2019" = "2019 (Marginal breeding)", 
                                                           "2020" = "2020 (Breeding)"))); plot_transect_year

summary(glmmTMB(song~total*as.factor(year), 
                family = binomial, 
                data = filter(data_transect, observer != "James Madden")))
# interaction effect between group size and year:
# 2018 - 2019: z = 1.59, df = 349, P = 0.11
# 2018 - 2020: z = 3.23, df = 349, P = 0.001

summary(glmmTMB(song~total*as.factor(2020/year), 
                family = binomial, 
                data = filter(data_transect, observer != "James Madden")))
# of course total is significant in this model (2020 is in the intercept), and there is a significant 
# interaction with 2019 and 2018, because in these years there is not a significant effect of group size 
# on song presence
# 2019 - 2020: z = 2.02, df = 349, P = 0.04


##### NEST RECORDINGS 2016 SONG POST BEHAVIOUR #####

### Read + tidy data
data_nestbox_recs <- read_csv("nestbox_recordings_2016_song_bouts.csv")
data_nestbox_songpost <- data_nestbox_recs %>% 
  select(c(recording = origin, individual, quality)) %>% 
  filter(!is.na(individual)) %>% 
  unique() %>% 
  mutate(quality = fct_relevel(quality, c("H", "M", "L")))
head(data_nestbox_songpost)

data_nestbox_songpost_high <- data_nestbox_songpost %>% 
  filter(quality == "H") %>% 
  count(recording, name = "singing_individuals") %>% 
  count(singing_individuals, name = "High")

data_nestbox_songpost_all <- data_nestbox_songpost %>% 
  group_by(recording) %>% 
  summarise(singing_individuals = max(individual)) %>% 
  count(singing_individuals, name = "Any") %>% 
  full_join(data_nestbox_songpost_high, by = "singing_individuals") %>% 
  replace_na(list(high = 0)) %>% 
  pivot_longer(-singing_individuals, names_to = "quality", values_to = "observations") %>% 
  mutate(quality = fct_relevel(quality, c("High", "Any")))

plot_nestbox <- ggplot(data = data_nestbox_songpost_all, mapping = aes(x = singing_individuals, y = observations, fill = quality)) +
  geom_bar(stat = "identity", position = "dodge2", color = "black") +
  scale_fill_manual(values = c("#BBBBBB", "white")) +
  scale_x_continuous(breaks = 1:15) +
  labs(x = "Singing males near nest", y = "Nest sites", fill = "Song quality") +
  theme_classic() +
  theme(legend.position = c(0.8,0.7),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width= unit(0.3, 'cm')
        #legend.background = element_rect(color = "black")
        ); plot_nestbox

#calculate mean and sd nr of males heard near the nestbox per quality 

data_nestbox_songpost_all <- within(data_nestbox_songpost_all, to_divide <-  singing_individuals*observations)
# average number of individuals heard at a nest with any quality spectrograms
mean_any <- with(filter(data_nestbox_songpost_all, quality == "Any"), 
                 sum(to_divide)/sum(observations)) 

sd_any <- with(filter(data_nestbox_songpost_all, quality == "Any"), 
               sqrt(sum(((singing_individuals - mean_any)^2)*observations)/sum(observations)))

# average number of individuals heard at a nest with high quality spectrograms
mean_high <- with(filter(data_nestbox_songpost_all, quality == "High", !is.na(observations)), 
                  sum(to_divide)/sum(observations)) 

sd_high <- with(filter(data_nestbox_songpost_all, quality == "High", !is.na(observations)), 
               sqrt(sum(((singing_individuals - mean_high)^2)*observations)/sum(observations)))

##### NEST RECORDINGS 2016 SONG OVER BREEDING STAGE #####

### Read + tidy data
data_nestbox_recs <- read_csv("nestbox_recordings_2016_song_bouts.csv")
data_nestbox_breeding <- read_csv("nestbox_recordings_2016_breeding.csv")

data_nestbox_combined <- data_nestbox_breeding %>% 
  full_join(data_nestbox_recs, by = "origin") %>% # combine datasets
  mutate(duration = end-start, # calculate duration of each song bout
         nest_stage = as.factor(nest_stage)) %>% # turn nest_stage into factor
  mutate(nest_stage = fct_relevel(nest_stage, c("Nest building", "Egg laying", "Incubation", "Nestlings"))) %>% # order the nest stages
  select(nestbox, nest_stage, abandoned, owner_id, individual, observation_hours, duration) %>% # drop unneeded vars
  filter(owner_id == individual & # remove all observations not by the owner
         abandoned == "no") %>%  # remove all nests that were abandoned
  group_by(nestbox, nest_stage, observation_hours) %>% 
  summarise(total_song_duration = sum(duration), count = n()) %>% 
  mutate(sing_act = total_song_duration/observation_hours)


head(data_nestbox_combined)
length(levels(as.factor(data_nestbox_combined$nestbox)))


plot_stage <- ggplot(data = data_nestbox_combined, mapping = aes(x = nest_stage, y = sing_act)) +
  geom_boxplot(alpha = 0.66, outlier.shape = NA) +
  geom_dotplot(color = "black", fill = "white", binaxis='y', stackdir='center', alpha = 1, binwidth = 1, dotsize = 1.2, stackratio = 1.5, position = position_nudge(x = 0.013)) +
  geom_line(data = data.frame(xa = c(2,2,3,3), ya = c(76,77,77,76)), aes(x = xa, y = ya), size = 0.3) +
  annotate("text", x = 2.5, y = 78, label = "**", size = 3) +
  geom_line(data = data.frame(xa = c(2,2,4,4), ya = c(81,82,82,81)), aes(x = xa, y = ya), size = 0.3) +
  annotate("text", x = 3, y = 83, label = "**", size = 3) +
  scale_y_continuous(breaks = seq(0, 80, 10), limits = c(0,84)) +
  labs(x = "Nest stage", y = "Male nest owner song (s / hr)") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 25, hjust = 1)); plot_stage

library(car)
qqPlot(resid(aov(sing_act ~ nest_stage, data = data_nestbox_combined)))
leveneTest(sing_act ~ nest_stage, data = data_nestbox_combined)
summary(aov(sing_act ~ nest_stage, data = data_nestbox_combined))
TukeyHSD(aov(sing_act ~ nest_stage, data = data_nestbox_combined))


##### SONGMETER RECORDINGS 2018+2019 #####

### Read + tidy data
data_songmeters <- read_csv("songmeters_2018-2019.csv")
data_songmeters$busy_score <- data_songmeters$busy_score
data_songmeters <- data_songmeters %>% 
  mutate(across(c('songmeter', 'observer'), # turn all relevant chr into factor
                .fns = as.factor
                ),
  month_index = if_else(year == 2018,
                        month-9, # 10,11,12 should become 1,2,3
                        month+3 # 1,2,3,.. should become 4,5,6,..
                        ),
  month_label = str_c(month.abb[month], as.character(year), sep = ' '), # combine year and month into 1 category
  .after = month
  )
head(data_songmeters)

data_songmeters_summ <- data_songmeters %>% 
  group_by(songmeter, month_index
           ) %>%
  summarise(across(c('song', 'busy_score'),
                   .fns = list(mean = mean)
                   )
            ) %>% 
  left_join(unique(select(data_songmeters, 
                          month_index, 
                          month_label)
                   ),
            by = "month_index"
            )
head(data_songmeters_summ)

plot_songmeter_song <- ggplot(data_songmeters_summ, mapping = aes(x = month_index, y = song_mean)) +
  geom_smooth(color = "black") +
  geom_point(mapping = aes(shape = songmeter, color = songmeter), alpha = 0.8, position = position_dodge(width = 0.33)) +
  scale_x_continuous(expand=c(0.012, 0.012), 
                     breaks = 1:13, 
                     labels = unique(data_songmeters_summ$month_label)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25)) +
  scale_shape(labels = c("1","2","3","4")) +
  scale_color_manual(values = c("#004488", "#BB5566", "#DDAA33", "#000000"),
                     labels = c("1","2","3","4")
                     ) +
  labs(x = "Month", y = "Daily hours with song  (%)", shape = "Site", color = "Site") +
  theme_classic() +
  theme(strip.background = element_rect(colour="white"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.6,0.8),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width= unit(0.3, 'cm')); plot_songmeter_song

plot_songmeter_activity <- ggplot(data_songmeters_summ, mapping = aes(x = month_index, y = busy_score_mean)) +
  geom_smooth(color = "black") +
  geom_point(mapping = aes(shape = songmeter, color = songmeter), alpha = 0.8, position = position_dodge(width = 0.33)) +
  scale_x_continuous(expand=c(0.012, 0.012), 
                     breaks = 1:13, 
                     labels = unique(data_songmeters_summ$month_label)) +
  scale_y_continuous(breaks = seq(0,4,1), 
                     labels = c("0 (no activity)", "1 (quiet overall)", "2 (not busy)", "3 (quite busy) ", "4 (very busy)")) +
  scale_shape(guide = "none") +
  scale_color_manual(guide = "none",
                     values = c("#004488", "#BB5566", "#DDAA33", "#000000")
  ) +
  labs(x = "Month", y = "Daily mean hourly activity score", shape = "Site", color = "Site") +
  theme_classic() +
  theme(strip.background = element_rect(colour="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)); plot_songmeter_activity


chisq.test(as.factor(data_songmeters$song), as.factor(data_songmeters$busy_score))
chisq.test(as.factor(data_songmeters$song), as.factor(data_songmeters$calls))


###### PLAYBACK EXPERIMENT 2017 ######

### Read + tidy data
pb_data <- read_csv("playback_experiments_2017.csv") %>% 
  mutate(across(c("Set", "Treatment"),
                .fns = as.factor)) %>% 
  filter(Continuation == "no") %>% # some visits are continuation visits where the birds present at the start of the video already arrived in the previous video
  filter(Playback_working == "yes") %>% # on rare occasions something went wrong with the playback, e.g. the speaker fell over 
  filter(Overlapping_visit == "no") %>% # exclude visits that happen when there are already birds present as these may attract birds on top of any playback effect
  filter(Exclude == "no") %>% # sometimes birds come in over the ground as a foraging group, these are not clearly responding to the playback so we exclude these and one other example where bird was not clearly responding as it was very far away
  select(c(Set:Song_present, -End_of_pb_s, -Latency))
pb_data

plot_latency <- ggplot(data = pb_data, mapping = aes(y = Latency_inclusive, x = S_first)) +
  geom_violin(mapping = aes(color = S_first)) +
  geom_dotplot(color = 'black', fill = 'white',
               binaxis='y', stackdir='center', binwidth = 15, dotsize = 1, stackratio = 1.2,
               position = position_nudge(x = 0.003)) +
  geom_line(data = data.frame(xa = c(2,2,3,3), ya = c(980,990,990,980)), aes(x = xa, y = ya), size = 0.2) +
  annotate("text", x = 2.5, y = 1010, label = "***", size = 2) +
  geom_line(data = data.frame(xa = c(1,1,3,3), ya = c(1060,1070,1070,1060)), aes(x = xa, y = ya), size = 0.2) +
  annotate("text", x = 2, y = 1090, label = "***", size = 2) +
  scale_color_manual(values = c("#004488", "#BB5566", "#DDAA33"), guide = "none") +
  scale_x_discrete(label = c(ZF = "Zebra finch", N = "Nightingale", aS = "Silence")) +
  scale_y_continuous(breaks = seq(0, 1000, 250), limits = c(0,1090)) +
  labs(x = "Playback", y = "Latency since playback end (s)") +
  theme_classic(); plot_latency 

summary(glmmTMB(Latency_inclusive~ZF_first + (1|Set),
                data = pb_data,
                family = tweedie))

summary(glmmTMB(Latency_inclusive~S_first + (1|Set),
                data = pb_data,
                family = tweedie))

pb_visits_S_first <- pb_data %>% # nice for the plots to have silence first
  select(Set:File) %>%
  count(Set, S_first) %>% 
  complete(Set, S_first, fill = list(n = 0)) %>% 
  rename(Visits = n)

pb_visits_ZF_first <- pb_data %>% # for the complete glmm comparisons, also need another factor first
  select(Set:File) %>%
  count(Set, ZF_first) %>% 
  complete(Set, ZF_first, fill = list(n = 0)) %>% 
  rename(Visits = n)

plot_visits <- ggplot(data = pb_visits_S_first, mapping = aes(x = S_first, y = Visits)) +
  geom_boxplot(mapping = aes(color = S_first),
               outlier.shape = NA) +
  geom_dotplot(color = 'black', fill = 'white',
               binaxis='y', stackdir='center', binwidth = 1, dotsize = 0.75, stackratio = 1.2,
               position = position_nudge(x = 0.005)) +
  geom_line(data = data.frame(xa = c(1,1,2,2), ya = c(31.65,32,32,31.65)), aes(x = xa, y = ya), size = 0.22) +
  annotate("text", x = 1.5, y = 32.5, label = "*", size = 2.15) +
  scale_color_manual(values = c("#004488", "#BB5566", "#DDAA33"), guide = "none") +
  scale_x_discrete(label = c(ZF = "Zebra finch", N = "Nightingale", aS = "Silence")) +
  labs(x = "Playback", y = "Site visits") +
  theme_classic(); plot_visits


summary(glmmTMB(Visits~ZF_first + (1|Set), 
                data = pb_visits_ZF_first, 
                family = genpois))

summary(glmmTMB(Visits~S_first + (1|Set), 
                data = pb_visits_S_first,  
                family = genpois))

count(filter(pb_data, Treatment == "Zebra Finch"), Visit_type, .drop = FALSE)
group_by(pb_data, Treatment) %>% count(Song_present, .drop = FALSE)



###### MANUSCRIPT PLOTS ######
# For 2-column formats (such as research articles and reviews), the sizes are 85 mm (1 column), 114 mm (1.5 columns), and 174 mm (full width of the page)
font_size <- theme(text = element_text(size = 7, face = 'plain'))

(plot_rainfall_breeding + theme(axis.title.x = element_blank()) + font_size) + (plot_breeding2 + font_size) +
  plot_layout(design = 'A\\nB', tag_level = 'new') +
  plot_annotation(tag_levels = list(c('A', 'B'))) & (theme = theme(plot.tag = element_text(face = 'bold')))
ggsave('figure1.tiff', width = 85, height = 85, units = 'mm', dpi = 1000, compression = 'lzw')

(plot_transect_year + font_size) + (plot_nestbox + font_size) + (plot_stage + font_size) +
  plot_layout(design = 'AA\\nBC', tag_level = 'new') +
  plot_annotation(tag_levels = list(c('A', 'B', 'C'))) & (theme = theme(plot.tag = element_text(face = 'bold')))
ggsave('figure2.tiff', width = 174, height = 174, units = 'mm', dpi = 1000, compression = 'lzw')

(plot_songmeter_song + font_size + theme(axis.title.x = element_blank())) + (plot_songmeter_activity + font_size) +
  plot_layout(design = 'A\\nB', tag_level = 'new') +
  plot_annotation(tag_levels = list(c('A', 'B'))) & (theme = theme(plot.tag = element_text(face = 'bold')))
ggsave('figure3.tiff', width = 174, height = 114, units = 'mm', dpi = 1000, compression = 'lzw')

(plot_latency + font_size + theme(axis.title.x = element_blank())) + (plot_visits + font_size) +
  plot_layout(design = 'A\\nB', tag_level = 'new') +
  plot_annotation(tag_levels = list(c('A', 'B'))) & (theme = theme(plot.tag = element_text(face = 'bold')))
ggsave('figure4.tiff', width = 85, height = 114, units = 'mm', dpi = 1000, compression = 'lzw')
