###################################################################
###################################################################

#Code for "Milkweed plants bought at nurseries may expose monarch caterpillars 
#to harmful pesticide residues" published in Biological Conservation.

#Author: Chris Halsch
#Date: July, 27, 2022

#Below you will find the code to replicate all results and main figures from the 
#manuscript. The code is split into blocks to replicate sections of the paper and
#can be run independently. See the below table of content:

#Table of contents:
#Lines 26-142: Creation of heat map for figure 1
#Lines 144-210: Creation of map for figure 1
#Lines 212-290: Summary results from the main text
#Lines 292-545: Analysis of pesticide richness and diversity and summary figures
#Lines 547-850: Analysis of exceedances and the creation of figure 2
#Lines 852-1022: Ordination and the creation of figure 3
#Lines 1024-1056: Change in concentration over two weeks

###################################################################
###################################################################

#############################################
##### Creation of heat map for figure 1 #####
#############################################

library(tidyverse)

setwd("~/Documents/MANUSCRIPTS/pesticides_part2/data_share/")

#Import data
dat <- read.csv("pest_data1.csv")
dat <- dat[-c(16:25),] #removes the resamples

#Locations for every sample
locs <- read.csv("site_locations.csv")
site_key <- locs %>% 
  select(-lon, -lat)

#Make a location key for the figure
locs_key <- locs %>% 
  filter(Code %in% unique(dat$Code)) %>%
  mutate(Site = paste(State, City, sep = "_")) %>% 
  select(-Code) %>% 
  distinct() %>% 
  arrange(factor(State, levels = rev(c("OR", "CA", "NV", "AZ", "NM", "CO", "TX", "OK", "NE", "KS", "IA", "MN", "IN", "OH", "NH"))), lon) %>% 
  select(-Site) %>% 
  mutate(lab_x = 1:25)

#Make a key for chemicals
chems <- read.csv("chems.csv") #import chemicals (there are more listed here than we detected)
chems <- chems %>% 
  filter(Cmpd %in% colnames(dat[8:68])) %>% 
  arrange(factor(type, levels = c("synergist", "herbicide", "fungicide", "insecticide")), desc(Cmpd)) %>% 
  mutate(lab_y = 1:61)

#Known lethal effects on monarchs from Table S4
lethal <- data.frame(Cmpd = c("Imidacloprid", "Thiamethoxam", "Chlorantraniliprole"),
                   leth = c(130, 420, 1.6))

#Known sub-lethal effects on monarchs from Table S4
sub_lethal <- data.frame(Cmpd = c("Azoxystrobin", "Trifloxystrobin"),
                   sub = c(0.67, 4.48))

#Create data for heatmep
dat1 <- dat %>% 
  gather(8:68, key = "Cmpd", value = "ppb") %>% 
  mutate(ppb = ifelse(ppb == 0, 0.001, ppb),
         ppb = as.numeric(ppb)) %>% 
  left_join(site_key, by = "Code") %>%
  group_by(State, City, Cmpd) %>% 
  summarise(av_ppb = mean(ppb)) %>% 
  left_join(chems, by = "Cmpd") %>% 
  left_join(locs_key, by = c("State", "City")) %>% 
  left_join(lethal, by = "Cmpd") %>%
  left_join(sub_lethal, by = "Cmpd") %>% 
  mutate(over1 = ifelse(av_ppb > leth, 1, 0), #Lethal overages
         over1_sym = ifelse(over1 == 1, paste('\\u2022'), " "), #Lethal overage symbol
         over2 = ifelse(av_ppb > sub, 1, 0), #Sublethal overages
         over2_sym = ifelse(over2 == 1, paste('\\u2022'), " ")) %>%  #Sublethal overage symbol 
  ungroup() %>% 
  select(-lat, -lon, -State, -City)

#Custom colors
colvec <- rev(paletteer::paletteer_c("grDevices::Sunset", 50))
colvec <- c("#FFFFFF", colvec)

#levels(reorder(factor(Cmpd), desc(lab_y)))
heat <- ggplot(data = dat1) +
  geom_tile(aes(factor(lab_x), as.numeric(reorder(factor(Cmpd), desc(lab_y))), fill = av_ppb), color = "black", size = 0.1, alpha = 1) + #builds the heat map
  geom_text(aes(factor(lab_x), as.numeric(reorder(factor(Cmpd), desc(lab_y))), label = over1_sym), color = "black", vjust = 0.48, size = 7) + #adds the white dots
  geom_text(aes(factor(lab_x), as.numeric(reorder(factor(Cmpd), desc(lab_y))), label = over2_sym), color = "white", vjust = 0.48, size = 7) + #adds the white dots
  geom_segment(aes(x = 0.5, y = 0.5, xend = 0.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 1.5, y = 0.5, xend = 1.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 2.5, y = 0.5, xend = 2.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 3.5, y = 0.5, xend = 3.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 4.5, y = 0.5, xend = 4.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 5.5, y = 0.5, xend = 5.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 6.5, y = 0.5, xend = 6.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 7.5, y = 0.5, xend = 7.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 9.5, y = 0.5, xend = 9.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 10.5, y = 0.5, xend = 10.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 11.5, y = 0.5, xend = 11.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 12.5, y = 0.5, xend = 12.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 13.5, y = 0.5, xend = 13.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 14.5, y = 0.5, xend = 14.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 21.5, y = 0.5, xend = 21.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_segment(aes(x = 25.5, y = 0.5, xend = 25.5, yend = 61.5), size = 1, color = "black") + #manually make darker lines. Using .5 puts the lines between tiles, not over them.
  geom_hline(aes(yintercept = 0.5), size = 1, color = "black") + #manually make darker lines
  geom_hline(aes(yintercept = 24.5), size = 1, color = "black") + #manually make darker lines
  geom_hline(aes(yintercept = 50.5), size = 1, color = "black") + #manually make darker lines
  geom_hline(aes(yintercept = 60.5), size = 1, color = "black") + #manually make darker lines
  geom_hline(aes(yintercept = 61.5), size = 1, color = "black") + #manually make darker lines
  scale_fill_gradientn(name = "Concentration\\n        (ppb)", colors = colvec, trans = "log", na.value = "red", breaks = c(min(dat1$av_ppb), 1, 10, 100, max(dat1$av_ppb)), labels = c(0, 1, 10, 100, 1000)) + #set color of tiles
  scale_x_discrete(labels = rev(c("\\nOR", "CA", "NV", "AZ", "NM", "CO", "TX", "\\nOK", "NE", "KS", "IA", "MN", "IN", "OH", "NH")), breaks = c(1,2,3,4,5,6,7,9,10,11,12,13,14,18,24)) +
  #manually add labels. I played with the spacing alot.
  scale_y_continuous(breaks = 1:length(levels(reorder(factor(dat1$Cmpd), desc(dat1$lab_y)))),
                     labels = levels(reorder(factor(dat1$Cmpd), desc(dat1$lab_y))),
                     expand = c(0, 0.1),
                     limits = c(0.5, 61.5),
                     sec.axis = sec_axis(~., breaks = c(12, 37, 55),labels = c("Insecticides", "   Fungicides", "   Herbicides"))) +
  coord_fixed(0.3) + #adjusts the shape of the tiles
  coord_flip() +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text.x.bottom = element_text(face="bold", color="black", size=12, angle = 60, hjust = 1.05, vjust = 1.05), 
        axis.text.x.top = element_text(face="bold", color="black", size=15), 
        axis.text.y = element_text(face="bold", color="black", size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(face="bold", color="black", size=12),
        legend.title = element_text(face="bold", color="black", size=15),
        legend.position = c(0.5, -0.3),
        legend.direction="horizontal",
        legend.key.width=unit(0.9,"cm"))
heat
#ggsave("raw_values_loq.png", device = "png", height = 8, width = 12, dpi = 300)

rm(list = ls())

##########################################################
##########     Creation of map for figure 1    ###########
##########################################################

setwd("~/Documents/MANUSCRIPTS/pesticides_part2/data_share/")

library(ggmap)
library(rgdal)
library(tidyverse)

ggmap::register_google(key = "AIzaSyCq0qfwfU7HTj1xlU3ZjwHxeqfYk1R3klA")

#Import shapefile of western US
west <- readOGR("cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
west <- subset(west, NAME != "Hawaii" & NAME != "Alaska" & NAME != "Puerto Rico")

#Import monarch routes shapefile
arr <- readOGR("monarch_route/arrow1.shp")

dat <- read.csv("pest_data.csv")
dat <- dat[-c(16:25),]

locs <- read.csv("site_locations.csv")
retail_info <- read.csv("retail_sample_details.csv")

dat2 <- dat %>% 
  gather(8:68, key = "Cmpd", value = "ppb") %>% 
  mutate(ppb = ifelse(ppb == "n.d.", 0, ppb),
         ppb = as.numeric(ppb)) %>% 
  filter(ppb > 0) %>% 
  group_by(Code) %>% 
  summarise(rich = length(Cmpd)) %>% 
  left_join(retail_info, by = "Code") %>%
  left_join(locs, by = "Code") %>% 
  mutate(storeID = paste(RetailID, City, State, sep = "_")) %>% 
  select(State, City, lat, lon, storeID) %>% 
  group_by(State, City, lat, lon) %>% 
  summarise(n_stores = length(unique(storeID)))

buffer <- 1.1

map <- ggmap(get_stamenmap(bbox = c(left = -124.7258 - buffer, bottom = 24.49813 - buffer, right = -66.94989 + buffer, top = 49.38436 + 2*buffer), zoom = 5, scale = 2, maptype ="terrain-background")) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "white", alpha = 0.03, show.legend = F) +
  geom_path(data = west, aes(long, lat, group = group), show.legend = F) +
  geom_path(data = arr, aes(long, lat, group = group), size = 3, color = "#F46D43", show.legend = F, lineend = "round", linejoin="mitre", arrow =  arrow(length = unit(0.7,"cm"), type = "open")) +
  geom_point(data = dat2, aes(x = lon, y = lat, size = n_stores), pch = 21, fill = "white", color = "black", size = 4, show.legend = F) + 
  ggrepel::geom_label_repel(data = dat2, aes(x  = lon, y = lat, label = City), label.padding =  unit(0.15, "lines"), nudge_x = 0.5, nudge_y = 0.5, size = 7, fontface = "bold") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.23, 0.83),
        legend.background = element_rect(fill=alpha('white', 0.85)),
        legend.key=element_blank(),
        rect = element_rect(fill = "transparent"), # all rectangles
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
map

#ggsave("~/Documents/MANUSCRIPTS/pesticides_part2/map.png", height = 10, width = 12)

rm(list = ls())

###############################################################
##########    Summary results from the main text     ##########
###############################################################

library(tidyverse)

setwd("~/Documents/MANUSCRIPTS/pesticides_part2/data_share/")

#Import data
dat <- read.csv("pest_data1.csv")
dat <- dat[-c(16:25),] #removes the resamples

#Locations for every sample
locs <- read.csv("site_locations.csv")
site_key <- locs %>% 
  select(-lon, -lat)

#Make a location key for the figure
locs_key <- locs %>% 
  filter(Code %in% unique(dat$Code)) %>%
  mutate(Site = paste(State, City, sep = "_")) %>% 
  select(-Code) %>% 
  distinct() %>% 
  arrange(factor(State, levels = rev(c("OR", "CA", "NV", "AZ", "NM", "CO", "TX", "OK", "NE", "KS", "IA", "MN", "IN", "OH", "NH"))), lon) %>% 
  select(-Site) %>% 
  mutate(lab_x = 1:25)

#Import chemical key
chems <- read.csv("chems.csv")

#Known lethal effects on monarchs from Table S4
lethal <- data.frame(Cmpd = c("Imidacloprid", "Thiamethoxam", "Chlorantraniliprole"),
                     leth = c(130, 420, 1.6))

#Known sub-lethal effects on monarchs from Table S4
sub_lethal <- data.frame(Cmpd = c("Azoxystrobin", "Trifloxystrobin"),
                         sub = c(0.67, 4.48))

#Number of lethal and sublethal exceedances of mean city concentration
#11/25 cities contain at least one store with an exceedance
over_by_city <- dat %>% 
  gather(8:68, key = "Cmpd", value = "ppb") %>% 
  mutate(ppb = ifelse(ppb == 0, 0, ppb),
         ppb = as.numeric(ppb)) %>% 
  filter(ppb > 0) %>% 
  left_join(site_key, by = "Code") %>% 
  group_by(State, City, Cmpd) %>% 
  summarise(av_ppb = mean(ppb)) %>% 
  left_join(chems, by = "Cmpd") %>% 
  left_join(locs_key, by = c("State", "City")) %>% 
  left_join(lethal, by = "Cmpd") %>%
  left_join(sub_lethal, by = "Cmpd") %>% 
  mutate(over1 = ifelse(av_ppb > leth, 1, 0),
         over2 = ifelse(av_ppb > sub, 1, 0)) %>% 
  ungroup() %>% 
  group_by(City, State) %>% 
  summarise(n_lethal_exceed = sum(over1, na.rm = T),
            n_sub_exceed = sum(over2, na.rm = T))

over_by_city

#Number of lethal and sublethal exceedances for each compound
#73 exceedances of Azoxystrobin and 22 of Trifloxystrobin
over_by_compound <- dat %>% 
  gather(8:68, key = "Cmpd", value = "ppb") %>% 
  mutate(ppb = ifelse(ppb == 0, 0, ppb),
         ppb = as.numeric(ppb)) %>% 
  left_join(lethal, by = "Cmpd") %>%
  left_join(sub_lethal, by = "Cmpd") %>% 
  mutate(over1 = ifelse(ppb > leth, 1, 0),
         over2 = ifelse(ppb > sub, 1, 0)) %>% 
  ungroup() %>% 
  group_by(Cmpd) %>% 
  summarise(n_lethal_exceed = sum(over1, na.rm = T),
            n_sub_exceed = sum(over2, na.rm = T))
  
over_by_compound

rm(list = ls())

####################################################################################
#####     Analysis of pesticide richness and diversity and summary figures     #####
####################################################################################

library(car)
library(ggbeeswarm)
library(ggsci)
library(lme4)
library(modEvA)
library(multcomp)
library(randomForest)
library(tidyverse)
library(wesanderson)

setwd("~/Documents/MANUSCRIPTS/pesticides_part2/data_share/")

#Import raw data
dat <- read.csv("pest_data1.csv")
dat <- dat[-c(16:25),] #removes the resamples

#Import location information
locs <- read.csv("site_locations.csv")

#Import chemical information
chems <- read.csv("chems.csv")

#Import retail information
retail_info <- read.csv("retail_sample_details.csv")

#Import store size information
stores <- read.csv("stores.csv")

#Make data for richness and diversity response variables
to_hill <- dat %>% 
  gather(8:68, key = "Cmpd", value = "ppb") %>% 
  mutate(ppb = ifelse(ppb == 0, 0, ppb),
        ppb = as.numeric(ppb)) %>% 
  spread(key = Cmpd, value = ppb) %>% 
  arrange(Code) %>% 
  left_join(retail_info, by = "Code") %>% 
  left_join(locs, by = "Code") %>% 
  mutate(storeID = paste(RetailID.x, City, State, sep = "_")) %>% 
  left_join(stores) %>% 
  mutate(Store.Size.Me = ifelse(Num_stores > 100, 3, 2), #Makes high store category (> 100 stores)
         Store.Size.Me = ifelse(Num_stores == 1, 1, Store.Size.Me), #Makes low store category (1 store)
         Store.Size.Me = as.factor(Store.Size.Me),
         Species = as.factor(Species.x),
         Grower.Size = factor(Grower.Size, levels = c(3, 2, 1)), #Reorder so that larger growers have larger values
         tag = as.factor(ifelse(tag > 0, 1, 0)),
         region = ifelse(lon < -108, "west", "east")) #East and west split by lon -108 (see map)

to_hill$h0 <- hillR::hill_taxa(to_hill[,8:68], q = 0) #make a richness column
to_hill$shan <- vegan::diversity(to_hill[,8:68], index = "shannon") #make a shannon column
to_hill$exp_shan <- exp(to_hill$shan) #Exp shannon column

##########################
#######  Richness  #######
##########################

#Model
mod0 <- glmer(h0 ~ Species + Store.Size.Me + region + tag  + (1|storeID), family = "poisson", data = to_hill)
summary(mod0)
Anova(mod0)
glht(mod0)

#Variance partitioning
p <- partR2::partR2(mod0, data = to_hill, R2_type = "marginal", partvars = c( partvars = c("Species", "Store.Size.Me", "tag")))
p
p_dat <- data.frame(term = p$R2$term, est = p$R2$estimate)

#Variance partitioning plot
modEvA::varPart(A = p_dat[2,2], B =  p_dat[3,2], C =  p_dat[4,2],
                AB =  p_dat[5,2], AC =  p_dat[6,2], BC = p_dat[7,2], ABC =  p_dat[8,2])

###########################
#######  Diversity  #######
###########################

#Model
mod1 <- lmerTest::lmer(exp_shan ~ Species + Store.Size.Me + region + tag + (1|storeID), data = to_hill)
summary(mod1)
Anova(mod1)
glht(mod1)

#Variance partitioning
p <- partR2::partR2(mod1, data = to_hill, R2_type = "marginal", partvars = c("Species", "Store.Size.Me", "tag"))
p_dat <- data.frame(term = p$R2$term, est = p$R2$estimate)

#Variance partitioning plot
modEvA::varPart(A = p_dat[2,2], B =  p_dat[3,2], C =  p_dat[4,2],
                AB =  p_dat[5,2], AC =  p_dat[6,2], BC = p_dat[7,2], ABC =  p_dat[8,2])

##################################################
##########         Summary plots        ##########
##################################################

#This code generates the supplemental plots. You can toggle between richness and diversity by switching the 
#comments in the blocks.

hill_plot <- to_hill %>% gather(90:92, key = "H", value = "div")

summ_hill <- to_hill %>% group_by(Species) %>% summarise(av = mean(exp_shan), se = (sd(exp_shan)/sqrt(length(exp_shan))))

#plots richness or diversity by retailer size
p1 <- ggplot(data = to_hill) +
  #Richness
  #geom_beeswarm(aes(x=as.factor(Store.Size.Me), y=h0, fill = as.factor(Store.Size.Me)), pch = 21, size = 3, colour='black', alpha = 0.25, show.legend = F, dodge.width=0.5) +
  #stat_summary(aes(as.factor(Store.Size.Me), h0), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  #stat_summary(aes(as.factor(Store.Size.Me), h0), fun = mean, geom = "point", size=5.5, colour='black', fill = "white", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  #labs(x = "Retailer size", y = "# of compounds") +
  
  #Diversity
  geom_beeswarm(aes(x=as.factor(Store.Size.Me), y=exp_shan, fill = as.factor(Store.Size.Me)), pch = 21, size = 3, colour='black', alpha = 0.25, show.legend = F, dodge.width=0.5) +
  stat_summary(aes(as.factor(Store.Size.Me), exp_shan), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  stat_summary(aes(as.factor(Store.Size.Me), exp_shan), fun = mean, geom = "point", size=5.5, colour='black', fill = "white", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  labs(x = "Retailer size", y = "Diversity of compounds") +
  
  scale_fill_manual(values = ghibli::ghibli_palette(name = "LaputaMedium", direction = -1), label = c("Single store", "2-100 stores", "> 100 stores")) +
  scale_x_discrete(labels = c("Single store", "2-100 stores", "> 100 stores")) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 13, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm")
  )

p1

#plots richness or diversity by milkweed species
p2 <- ggplot(data = to_hill) +
  #Richness
  #geom_beeswarm(aes(x=as.factor(Species), y=h0, fill = as.factor(Species)), pch = 21, size = 3, colour='black', alpha = 0.25, show.legend = F, dodge.width=0.5) +
  #stat_summary(aes(as.factor(Species), h0), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  #stat_summary(aes(as.factor(Species), h0), fun = mean, geom = "point", size=5.5, colour='black', fill = "white", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  #labs(x = "Species", y = "# of compounds") +
  
  #Diversity
  geom_beeswarm(aes(x=as.factor(Species), y=exp_shan, fill = as.factor(Species)), pch = 21, size = 3, colour='black', alpha = 0.25, show.legend = F, dodge.width=0.5) +
  stat_summary(aes(as.factor(Species), exp_shan), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  stat_summary(aes(as.factor(Species), exp_shan), fun = mean, geom = "point", size=5.5, colour='black', fill = "white", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  labs(x = "Species", y = "Diversity of compounds") +
  
  scale_fill_manual(values = c("#6B4E71", "#F25F5C", "#FFE066", "#247BA0", "#61F2C2")) +
  scale_x_discrete(labels = c("A. cur.", "A. fas.", "A. inc.", "A. spe.", "A. tub.")) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 13, colour = "black", face = "italic"),
        axis.text.y = element_text(size = 13, colour = "black"),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 13, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm")
  )

p2

#plots richness or diversity by region
p3 <- ggplot(data = to_hill) +
  #Richness
  #geom_beeswarm(aes(x=as.factor(region), y=h0, fill = as.factor(region)), pch = 21, size = 3, colour='black', alpha = 0.25, show.legend = F, dodge.width=0.5) +
  #stat_summary(aes(as.factor(region), h0), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  #stat_summary(aes(as.factor(region), h0), fun = mean, geom = "point", size=5.5, colour='black', fill = "white", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  #labs(x = "Region", y = "# of compounds") +
  
  #Diversity
  geom_beeswarm(aes(x=as.factor(region), y=exp_shan, fill = as.factor(region)), pch = 21, size = 3, colour='black', alpha = 0.25, show.legend = F, dodge.width=0.5) +
  stat_summary(aes(as.factor(region), exp_shan), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  stat_summary(aes(as.factor(region), exp_shan), fun = mean, geom = "point", size=5.5, colour='black', fill = "white", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  labs(x = "Region", y = "Diversity of compounds") +
  
  scale_fill_manual(values = ghibli::ghibli_palette(name = "KikiMedium")[3:5]) +
  scale_x_discrete(labels = c("East", "West")) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 13, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm")
  )

p3

#plots richness or diversity by wildlife label
p4 <- ggplot(data = subset(to_hill, !is.na(tag == T))) +
  #Richness
  #geom_beeswarm(aes(x=as.factor(tag), y=h0, fill = as.factor(tag)), pch = 21, size = 3, colour='black', alpha = 0.25, show.legend = F, dodge.width=0.5) +
  #stat_summary(aes(as.factor(tag), h0), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  #stat_summary(aes(as.factor(tag), h0), fun = mean, geom = "point", size=5.5, colour='black', fill = "white", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  #labs(x = "Label", y = "# of compounds") +
  
  #Diversity
  geom_beeswarm(aes(x=as.factor(tag), y=exp_shan, fill = as.factor(tag)), pch = 21, size = 3, colour='black', alpha = 0.25, show.legend = F, dodge.width=0.5) +
  stat_summary(aes(as.factor(tag), exp_shan), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  stat_summary(aes(as.factor(tag), exp_shan), fun = mean, geom = "point", size=5.5, colour='black', fill = "white", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  labs(x = "Label", y = "Diversity of compounds") +
  
  scale_fill_manual(values = wes_palette("Darjeeling2")[2:3], label = c("No label", "Label")) +
  scale_x_discrete(label = c("No label", "Label")) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 13, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm")
  )

p4

#merge plots
ggpubr::ggarrange(p1, p2, p3, p4, 
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2,
                  legend = "bottom", widths = c(1.4, 1),
                  common.legend = T)

#ggsave("figS2_rich.png", device = "png", height = 8, width = 10, dpi = 300)
#ggsave("figS3_div.png", device = "png", height = 8, width = 10, dpi = 300)

rm(list = ls())

####################################################################
#####   Analysis of exceedances and the creation of figure 2   #####
####################################################################

library(car)
library(ggbeeswarm)
library(ggpubr)
library(lme4)
library(multcomp)
library(randomForest)
library(tidyverse)
library(wesanderson)

setwd("~/Documents/MANUSCRIPTS/pesticides_part2/data_share/")

#Import raw data
dat <- read.csv("pest_data1.csv")
dat <- dat[-c(16:25),] #removes the resamples

#Import location information
locs <- read.csv("site_locations.csv")

#Import chemical information
chems <- read.csv("chems.csv")

#Import retail information
retail_info <- read.csv("retail_sample_details.csv")

#Import store size information
stores <- read.csv("stores.csv")

#Using the same lethal and sub_lethal exceedances from Figure 1. See Table S4.
#Known lethal effects on monarchs from Table S4
lethal <- data.frame(Cmpd = c("Imidacloprid", "Thiamethoxam", "Chlorantraniliprole"),
                     leth = c(130, 420, 1.6))

#Known sub-lethal effects on monarchs from Table S4
sub_lethal <- data.frame(Cmpd = c("Azoxystrobin", "Trifloxystrobin"),
                         sub = c(0.67, 4.48))


dat1 <- dat %>% 
  gather(8:68, key = "Cmpd", value = "ppb") %>% 
  mutate(ppb = ifelse(ppb == 0, 0, ppb),
         ppb = as.numeric(ppb)) %>% 
  left_join(lethal, by = "Cmpd") %>% 
  left_join(sub_lethal, by = "Cmpd") %>% 
  mutate(over1 = ifelse(ppb > leth, 1, 0),
         over1 = ifelse(is.na(over1)==T,0, over1),
         over2 = ifelse(ppb > sub, 1, 0),
         over2 = ifelse(is.na(over2)==T,0, over2),
         over3 = over1+over2) %>% 
  group_by(Code) %>% 
  summarise(exceed = sum(over2)) %>% #changed for subs
  left_join(retail_info, by = "Code") %>%
  left_join(locs, by = "Code") %>% 
  mutate(storeID = paste(RetailID, City, State, sep = "_")) %>% 
  left_join(stores) %>% 
  mutate(Store.Size.Me = ifelse(Num_stores > 100, 3, 2),
         Store.Size.Me = ifelse(Num_stores == 1, 1, Store.Size.Me),
         Store.Size.Me = as.factor(Store.Size.Me),
         Species = as.factor(Species),
         Grower.Size = factor(Grower.Size, levels = c(3, 2, 1)),
         tag = ifelse(tag>0, 1, 0),
         tag = factor(tag),
         exceed_bin = ifelse(exceed > 0, 1, 0), #Makes exceedance a binary category
         region = ifelse(lon < -108, "west", "east"))
  
#Model
mod1 <- glmer(exceed_bin ~ Species + Store.Size.Me + region + tag + (1|storeID), 
              family = "binomial", glmerControl(optimizer = "nlminbwrap", optCtrl = list(maxfun = 2e5)), data = dat1)
summary(mod1)
Anova(mod1)
glht(mod1)

#Variance partioning
p <- partR2::partR2(mod1, data = dat1, R2_type = "marginal", partvars = c("Species", "Store.Size.Me", "tag"))
p_dat <- data.frame(term = p$R2$term, est = p$R2$estimate)

modEvA::varPart(A = p_dat[2,2], B =  p_dat[3,2], C =  p_dat[4,2],
                AB =  p_dat[5,2], AC =  p_dat[6,2], BC = p_dat[7,2], ABC =  p_dat[8,2])

#Plot prob exceedance by Retailer size
plot_dat <- dat1 %>% 
  group_by(storeID, Store.Size.Me) %>% 
  summarise(prob_contam = sum(exceed_bin)/length(exceed_bin))

pp1 <- ggplot(data = plot_dat) +
  geom_beeswarm(aes(x=as.factor(Store.Size.Me), y=prob_contam, fill = as.factor(Store.Size.Me)), pch = 21, size = 3, colour='black', alpha = 1, show.legend = F) +
  stat_summary(aes(as.factor(Store.Size.Me), prob_contam), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  stat_summary(aes(as.factor(Store.Size.Me), prob_contam), fun = mean, geom = "point", size=5.5, colour='black', fill = "gray75", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  labs(x = "Retailer size", y = "# of compounds") +
  scale_fill_manual(values = ghibli::ghibli_palette(name = "LaputaMedium", direction = -1), label = c("Single store", "2-100 stores", "> 100 stores")) +
  scale_x_discrete(labels = c("Single store", "2-100 stores", "> 100 stores")) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 15, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.6, 0.01, 0.01, 0.1), "cm")
  )

pp1

#Plot prob exceedance by Species
plot_dat <- dat1 %>% 
  group_by(storeID, Species) %>% 
  summarise(prob_contam = sum(exceed_bin)/length(exceed_bin))

pp2 <- ggplot(data = plot_dat) +
  geom_beeswarm(aes(x=as.factor(Species), y=prob_contam, fill = as.factor(Species)), pch = 21, size = 3, colour='black', alpha = 1, show.legend = F) +
  stat_summary(aes(as.factor(Species), prob_contam), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  stat_summary(aes(as.factor(Species), prob_contam), fun = mean, geom = "point", size=5.5, colour='black', fill = "gray75", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  labs(x = "Species", y = "# of compounds") +
  scale_fill_manual(values = c("#6B4E71", "#F25F5C", "#FFE066", "#247BA0", "#61F2C2")) +
  scale_x_discrete(labels = c("A. cur.", "A. fas.", "A. inc.", "A. spe.", "A. tub.")) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11, colour = "black", face="italic"),
        axis.text.y = element_text(size = 11, colour = "black"),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 15, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.6, 0.01, 0.01, 0.1), "cm")
  )

pp2

#Plot prob exceedance by region
plot_dat <- dat1 %>% 
  group_by(storeID, region) %>% 
  summarise(prob_contam = sum(exceed_bin)/length(exceed_bin))

pp3 <- ggplot(data = plot_dat) +
  geom_beeswarm(aes(x=as.factor(region), y=prob_contam, fill = as.factor(region)), pch = 21, size = 3, colour='black', alpha = 1, show.legend = F) +
  stat_summary(aes(as.factor(region), prob_contam), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  stat_summary(aes(as.factor(region), prob_contam), fun = mean, geom = "point", size=5.5, colour='black', fill = "gray75", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  labs(x = "Region", y = "# of compounds") +
  scale_fill_manual(values = ghibli::ghibli_palette(name = "KikiMedium")[3:5]) +
  scale_x_discrete(labels = c("East", "West")) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 15, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.6, 0.01, 0.01, 0.1), "cm")
  )

pp3

#Plot prob exceedance by wildlife label
plot_dat <- dat1 %>% 
  group_by(storeID, tag) %>% 
  summarise(prob_contam = sum(exceed_bin)/length(exceed_bin))

pp4 <- ggplot(data = subset(plot_dat, !is.na(tag == T))) +
  geom_beeswarm(aes(x=as.factor(tag), y=prob_contam, fill = as.factor(tag)), pch = 21, size = 3, colour='black', alpha = 1, show.legend = F) +
  stat_summary(aes(as.factor(tag), prob_contam), fun.data = mean_cl_boot, geom = "errorbar", color='black', width = 0.2, size=1, show.legend = F, position = position_dodge(width = 0.5)) + 
  stat_summary(aes(as.factor(tag), prob_contam), fun = mean, geom = "point", size=5.5, colour='black', fill = "gray75", pch = 21, show.legend = F, position = position_dodge(width = 0.5)) +
  labs(x = "Label", y = "# of compounds") +
  scale_fill_manual(values = wes_palette("Darjeeling2")[2:3], label = c("No label", "Label")) +
  scale_x_discrete(labels = c("No label", "Label")) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 15, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.6, 0.01, 0.01, 0.1), "cm")
  )

pp4

#Prepare data for random forest to identify which compounds are the most associated with labels
to_rf <- dat %>% 
  gather(8:68, key = "Cmpd", value = "ppb") %>% 
  mutate(ppb = ifelse(ppb == 0, 0, ppb),
         ppb = as.numeric(ppb),
         ppb = ifelse(ppb == 0, 0 ,1)) %>% 
  spread(key = Cmpd, value = ppb) %>% 
  left_join(dat1, by = "Code") %>%  
  dplyr::select(c(8:68, tag)) %>% 
  na.omit()

#Runs random forest. Here tag is a binary category (so classification is run) and the predictors are all detected chemicals
set.seed(545)
rfm <- randomForest(factor(to_rf$tag) ~ ., data = to_rf[,1:61], localImp = T)

#Extract and arrange most important variables (prepare for plotting)
vip <- data.frame(varImpPlot(rfm)) %>% arrange(MeanDecreaseAccuracy)
vip$names <- rownames(vip)
vip <- vip[(nrow(vip)-19):nrow(vip),]
vip$names <- factor(vip$names,  levels = vip$names)
colnames(vip)[3] <- "Cmpd"
vip <- left_join(vip, chems)
vip$Cmpd <- factor(vip$Cmpd,  levels = vip$Cmpd)
vip$type <- factor(vip$type,  levels = c("insecticide", "fungicide", "herbicide", "synergist"))

#Make random forest plot
rfp1 <- ggplot(data = vip) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = Cmpd, yend = Cmpd), size = 1, linetype = "dotted", show.legend = F) +
  geom_point(aes(MeanDecreaseAccuracy, Cmpd, shape = type), fill = "gray75", size = 5) +
  scale_shape_manual(values = c(22,23,24,25), labels = c("Insecticide", "Fungicide", "Herbicide", "Synergist")) +
  labs(x = "Increase MSE") +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, colour = "black"),
        legend.position = c(0.79, 0.14),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.6, 0.01, 0.01, 0.25), "cm")
  )

rfp1

#make additional panel for most important compounds
to_plot <- to_rf %>% 
  dplyr::select(tag, levels(vip$Cmpd)[(nrow(vip)-2):nrow(vip)]) %>% 
  gather(2:ncol(.), key = "Cmpd", value = "pres") %>% 
  group_by(tag, Cmpd) %>% 
  summarise(tot = sum(pres))
to_plot$Cmpd <- factor(to_plot$Cmpd, levels = c("Azoxystrobin", "Trifloxystrobin", "Piperonyl.butoxide"))

rfp2 <- ggplot(data = to_plot) +
  geom_col(aes(x = Cmpd, y = tot, fill = factor(tag)), position="dodge", size=3, alpha = 1, show.legend = T) +
  labs(y = "Number of detections", fill = "Presence of label") +
  scale_fill_manual(values = wes_palette("Darjeeling2")[2:3], label = c("No label", "Label")) +
  guides(fill=guide_legend(reverse = T)) +
  coord_flip() +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, colour = "black"),
        legend.position = c(0.82, 0.095),
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.6, 0.01, 0.01, 0.25), "cm")
  )
  
rfp2

#merge plots
comp1 <- annotate_figure(ggarrange(pp1, pp2, pp3, pp4, ncol = 1, nrow = 4, labels = c("A", "B", "C", "D")), 
                         left = text_grob("% of plants with exceedance", face = "bold", size = 18, rot = 90))

comp2 <- ggarrange(rfp1, rfp2, ncol = 1, nrow = 2, labels = c("E", "F"))

ggarrange(comp1, comp2,
          ncol = 2, nrow = 1, legend = "none", widths = c(1, 1.4))

#ggsave("fig2.png", device = "png", height = 12, width = 10, dpi = 300)

rm(list = ls())

#################################################################
##########   Ordination and the creation of figure 3   ##########
#################################################################

library(tidyverse)
library(wesanderson)

setwd("~/Documents/MANUSCRIPTS/pesticides_part2/data_share/")

#Import raw data
dat <- read.csv("pest_data1.csv")
dat <- dat[-c(16:25),] #removes the resamples

#Import location information
locs <- read.csv("site_locations.csv")

#Import retail information
retail_info <- read.csv("retail_sample_details.csv")

#Import store size information
stores <- read.csv("stores.csv")

dat1 <- dat %>% 
  gather(8:68, key = "Cmpd", value = "ppb") %>% 
  mutate(ppb = ifelse(ppb == 0, 0, ppb),
         ppb = as.numeric(ppb)) %>% 
  spread(key = "Cmpd", value = ppb) %>% 
  left_join(locs, by = "Code") %>% 
  left_join(retail_info, by = "Code") %>% 
  mutate(storeID = paste(RetailID.x, City, State, sep = "_")) %>% 
  left_join(stores, by = "storeID") %>% 
  mutate(Store.Size.Me = ifelse(Num_stores > 100, 3, 2),
         Store.Size.Me = ifelse(Num_stores == 1, 1, Store.Size.Me),
         tag = factor(tag),
         region = ifelse(lon < -108, "west", "east"))

#Subset to only a matrix of concentrations
to_pcoa <- dat1[,8:68]

#Assign all non-zero values a 1 (presence/absence)
to_pcoa[to_pcoa > 0] <- 1

#Jaccard distance
dis <- vegan::vegdist(to_pcoa, method = "jaccard")

#perform pcoa
pc <- ape::pcoa(dis)

#Create plotting data frame for retailer size
plot_dat <- data.frame(pc1 = pc$vectors[,1],
                       pc2 = pc$vectors[,2],
                       store.size = as.factor(dat1$Store.Size.Me),
                       species = as.factor(dat1$Species.x),
                       grower = as.factor(dat1$Grower.Size),
                       tag = as.factor(dat1$tag),
                       region = as.factor(dat1$region))

p1 <-ggplot(data = plot_dat) +
  geom_point(aes(pc1, pc2, fill = as.factor(store.size)), pch = 21, size = 3, alpha = 0.75) +
  labs(x = paste0("PCoA1 (", round(pc$values$Relative_eig[1], 2)*100, "%)"), y = paste0("PCoA2 (", round(pc$values$Relative_eig[2], 2)*100, "%)"), fill = "Retailer size") +
  scale_fill_manual(values = ghibli::ghibli_palette(name = "LaputaMedium", direction = -1), label = c("Single store", "2-100 stores", "> 100 stores")) +
  guides(fill=guide_legend(nrow=1)) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank())
p1

#Create plotting data frame for milkweed species
plot_dat <- data.frame(pc1 = pc$vectors[,1],
                       pc2 = pc$vectors[,2],
                       store.size = as.factor(dat1$Store.Size.Me),
                       species = as.factor(dat1$Species.x),
                       grower = as.factor(dat1$Grower.Size),
                       tag = as.factor(dat1$tag),
                       region = as.factor(dat1$region))

p2 <-ggplot(data = plot_dat) +
  geom_point(aes(pc1, pc2, fill = as.factor(species)), pch = 22, size = 3, alpha = 0.75) +
  labs(x = paste0("PCoA1 (", round(pc$values$Relative_eig[1], 2)*100, "%)"), y = paste0("PCoA2 (", round(pc$values$Relative_eig[2], 2)*100, "%)"), fill = "Species") +
  scale_fill_manual(values = c("#6B4E71", "#F25F5C", "#FFE066", "#247BA0", "#61F2C2"), labels = c("A. cur.", "A. fas.", "A. inc.", "A. spe.", "A. tub.")) +
  guides(fill=guide_legend(nrow=1)) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13, colour = "black", face = "italic"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank())
p2

#Create plotting data frame for region
plot_dat <- data.frame(pc1 = pc$vectors[,1],
                       pc2 = pc$vectors[,2],
                       store.size = as.factor(dat1$Store.Size.Me),
                       species = as.factor(dat1$Species.x),
                       grower = as.factor(dat1$Grower.Size),
                       tag = as.factor(dat1$tag),
                       region = as.factor(dat1$region))

p6 <-ggplot(data = plot_dat) +
  geom_point(aes(pc1, pc2, fill = as.factor(region)), pch = 23, size = 3, alpha = 0.75) +
  labs(x = paste0("PCoA1 (", round(pc$values$Relative_eig[1], 2)*100, "%)"), y = paste0("PCoA2 (", round(pc$values$Relative_eig[2], 2)*100, "%)"), fill = "Region") +
  scale_fill_manual(values = ghibli::ghibli_palette(name = "KikiMedium")[3:5], labels = c("East", "West")) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank())

p6

#Create plotting data frame for wildlife label
p8 <-ggplot(data = subset(plot_dat, !is.na(tag == T))) +
  geom_point(aes(pc1, pc2, fill = as.factor(tag)), pch = 24, size = 3, alpha = 0.75) +
  labs(x = paste0("PCoA1 (", round(pc$values$Relative_eig[1], 2)*100, "%)"), y = paste0("PCoA2 (", round(pc$values$Relative_eig[2], 2)*100, "%)"), fill = "Label") +
  scale_fill_manual(values = wes_palette("Darjeeling2")[2:3], label = c("No label", "Label")) +
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold", size = 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13, colour = "black"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.spacing.x = unit(0.1, "cm"),
        legend.background = element_rect(fill=alpha('white', 0.0)),
        legend.box.background = element_rect(colour = "black"),
        panel.grid = element_blank())

p8

#merge plots
ggpubr::ggarrange(p1, p2, p6, p8, 
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2,
                  legend = "bottom", heights = c(1.1, 1),
                  common.legend = F)

#ggsave("fig3_ord.png", device = "png", height = 10, width = 10, dpi = 300)

rm(list = ls())

##################################################
##### Change in concentration over two weeks #####
##################################################

library(tidyverse)

setwd("~/Documents/MANUSCRIPTS/pesticides_part2/data_share/")

#Import data
dat <- read.csv("pest_data1.csv")

period1 <- dat[c(1:10),]
period1$period <- 1

period2 <- dat[c(16:25),]
period2$period <- 2
period2$Code <- substr(period2$Code, 1, nchar(period2$Code)-2)

dat1 <- rbind(period1, period2)
dat1 <- dat1 %>% 
  dplyr::select(Code, period, everything()) %>% 
  gather(8:68, key = "Cmpd", value = "ppb") %>% 
  mutate(ppb = ifelse(ppb == 0, 0.001, ppb),
         ppb = as.numeric(ppb)) %>% 
  dplyr::select(Code, period, Cmpd, ppb) %>% 
  spread(key = period, value = ppb) %>% 
  mutate(delta = `2` - `1`) %>% 
  group_by(Cmpd) %>% 
  summarise(av_delta = mean(delta))

dat1

rm(list = ls())
