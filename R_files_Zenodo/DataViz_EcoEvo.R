#Graphing and Graphics
#By the awesome Chloe Fouilloux

#If you like my data vis, gimme a shout! Thanks for your interest. 

#First, let's load some packages!
library(ggplot2)
library(lubridate)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(Hmisc)
###
#For raincloud plot
library(PupillometryR)
library(rlang)
library(forcats)
#For images
library(png)

#Let's load some data
frog1<- read.csv("frogpool19_ECOEVO.csv", 
                 header = T,fileEncoding="UTF-8-BOM")

pool_2020 <- read.csv("frogpool20_ECOEVO.csv", header = T, sep= ",") 

meteo<- read.csv("PARARE_Meteo_Janv-Juin2020.csv", 
                 header = T, fileEncoding="UTF-8-BOM") #ANDRIUS

#####################################################################################
                              # Tidy data is the best data #
#####################################################################################

#Quick clean of 2019 data
frog1 <- frog1 %>%
  mutate(Total_Tads = tp_femo_hahneli + oophagus_tp + tadpoles_tinc, 
         Tad_Y_N = ifelse(Total_Tads > 0, 1, 0), 
         Tad_Y_N = as.factor(Tad_Y_N), 
         species = ifelse(tp_femo_hahneli> 0, "Femoralis", 
                          ifelse(oophagus_tp >0, "Oophagus", "Tinctorius")), 
         species = as.factor(species))

frog1$Species <- c("All")

#Assign 2020 pool types
pool_2020[pool_2020$Type %in% c("Waka_Pool", "Buttress", "Tree_Hole"), "pool_type"] <- "alive"
pool_2020[pool_2020$Type %in% c("Palm", "Fallen_Tree"), "pool_type"] <- "dead"
pool_2020[pool_2020$Type %in% c("Andrius"), "pool_type"] <- "arboreal"

###############
### SPECIES SUBSET
###########
##Now I have to think how to do this for 3 species.
tincpool<- frog1 %>%
  dplyr::select(-oophagus_tp, -tp_femo_hahneli)%>%
  mutate(Tadpoles = tadpoles_tinc,
         Species = "Tinctorius")%>%
  dplyr::select(height, Tadpoles, Species,
                pool_type, water_capacity, 
                leaflit_vol_ml, KH, hardness, no3,
                salinity, invert_diversity, predator_count,
                predator_size_sum, amphib_diversity, Total_Tads)%>%
  filter( Tadpoles > 0)

oopool<- frog1 %>%
  dplyr::select(-tadpoles_tinc, -tp_femo_hahneli) %>%
  mutate(Tadpoles = oophagus_tp,
         Species = "Oophagus")%>%
  dplyr::select(height, Tadpoles, Species,
                pool_type, water_capacity, 
                leaflit_vol_ml, KH, hardness, no3,
                salinity, invert_diversity, predator_count,
                predator_size_sum, amphib_diversity, Total_Tads)%>%
  filter( Tadpoles > 0)

femopool<- frog1 %>%
  dplyr::select(-tadpoles_tinc, -oophagus_tp)%>%
  mutate(Tadpoles = tp_femo_hahneli,
         Species = "Femoralis")%>%
  dplyr::select(height, Tadpoles, Species,
                pool_type, water_capacity, 
                leaflit_vol_ml, KH, hardness, no3,
                salinity, invert_diversity, predator_count,
                predator_size_sum, amphib_diversity, Total_Tads)%>%
  filter( Tadpoles > 0)

pools_species<- full_join(tincpool, femopool, 
                          by = c("height", "Species", "Tadpoles", 
                                 "KH", "Total_Tads", "leaflit_vol_ml", 
                                 "salinity", "predator_size_sum", 
                                 "predator_count", "hardness", "pool_type", "water_capacity"))

pools_species<- full_join(pools_species, oopool, 
                          by = c("height", "Species", "Tadpoles", 
                                 "KH", "Total_Tads", "leaflit_vol_ml",
                                 "salinity", "predator_size_sum", 
                                 "predator_count", "hardness", "pool_type","water_capacity")) 

pools_species$Species <- as.factor(pools_species$Species)

###TERRESTRIAL POOLS ----
######################################

pool_terr <- frog1 %>%
  filter(height < 220) %>%
  mutate(Total_Tads = tp_femo_hahneli + oophagus_tp + tadpoles_tinc, 
         Tad_Y_N = ifelse(Total_Tads > 0, 1, 0), 
         Tad_Y_N = as.factor(Tad_Y_N), 
         Water_Cap = ifelse(water_capacity> 5000, "Large", 
                            ifelse(water_capacity >1000, "Medium", "Small")), 
         Water_Cap = as.factor(Water_Cap),
         species = ifelse(tp_femo_hahneli> 0, "Femoralis", 
                          ifelse(oophagus_tp >0, "Oophagus", 
                                 ifelse(tadpoles_tinc > 0, "Tinctorius", "None"))),
         species = as.factor(species))

pool_terr$Water_Cap<- fct_relevel(pool_terr$Water_Cap, c("Small", "Medium", "Large"))
pool_terr$Water_Cap1<- recode(pool_terr$Water_Cap, 
                              Small = "Small pool", 
                              Medium = "Medium pool", 
                              Large = "Large pool")

pool_terr$species<- fct_relevel(pool_terr$species, c("Femoralis", "Oophagus", "Tinctorius", "None"))

##################################
#### Weather! 2020 ----

meteo$timestamp<- dmy_hm(meteo$TIMESTAMP)
class(meteo$Rain_mm_Tot)
meteo$Rain_mm <- as.numeric(as.character(meteo$Rain_mm_Tot)) #workaround, cannot go as.numeric directly
#now let's subset for the month we were there
precip_Feb <- meteo %>%
  filter(timestamp >= as.Date('2020-01-31') & timestamp <= as.Date('2020-02-27')) %>%
  tidyr::drop_na()%>%
  group_by(week = week(timestamp))%>%
  summarise (week_ave_rain = mean(Rain_mm), 
             week_total_rain = sum(Rain_mm))

precip_Feb$Week <- c("1","2", "3", "4", "5")

#So now let's plot this together with tadpole numbers
#Now let's combine data_frames
chloe_rain <- pool_2020 %>%
  mutate(Week = as.character(Week))%>%
  dplyr::select("Week", "pH", "Dt_Tadpole_Num", "pool_type")

chloe_rain_all<- dplyr::full_join(chloe_rain, precip_Feb,by = "Week")
chloe_rain_all<-chloe_rain_all %>%
  mutate(Week = as.numeric(Week))
#####################################################################################

#In order of appearance:

#FIGURE 1---- 
#Photograph

#FIGURE 2----
#Line drawing by the immensely talented Andrius Pasukonis

#FIGURE 3----
#Tadpole presence across the vertical landscape

#LETS GET PLOTTIN'!
a1<-ggplot(pools_species)+
  geom_hline(yintercept = 220, linetype = "dashed")+
  geom_point(aes(y = height, x = Species, color = Species), 
             position = position_jitter(width = .15), 
             alpha = 0.6, size = 3)+
  geom_point(data = frog1, aes(y = height, x = Species, shape = Tad_Y_N), 
             position = position_jitter(width = .3), size = 3, 
             color = "black", alpha= 0.8)+
  ylab("Height (cm)")+
  expand_limits(x = 4) +
  scale_color_manual( values = c("darkgrey", "#ffd16b", "#80e6ff"),
                      labels = c(bquote(italic("A. femoralis")),
                                 bquote(italic("O. oophagus")),
                                 bquote(italic("D. tinctorius")),
                                 "None")) +
  scale_fill_manual( values = c("darkgrey", "#ffd16b", "#80e6ff"))+
  scale_shape_manual(values = c(2, 17), 
                     name = "Tadpole presence", labels = c("No", "Yes"))+
  scale_x_discrete(labels =
                     c("All",
                       bquote(italic("A. femoralis")),
                       bquote(italic("O. oophagus")),
                       bquote(italic("D. tinctorius"))))+
  labs(tag = "A") +
  coord_cartesian(ylim = c(0, 2300))+
  theme_bw(base_line_size = 0)

a1

legend1<- get_legend(a1)

a1<- a1 + theme(legend.position = "none")


#Just terrestrial pools
b1<- ggplot(pools_species, 
            aes(y = height, x = Species))+
  geom_flat_violin(aes(fill = Species), color = "transparent", 
                   position = position_nudge(x = .1, y = 0),
                   alpha = 0.4)+
  geom_boxplot(width = .1, aes(fill = Species), 
               outlier.shape = NA, alpha = 0.3) +
  geom_point(aes(y = height, color = Species), 
             position = position_jitter(width = .1), 
             alpha = 0.6, size = 3)+
  theme_bw(base_line_size = 0)+
  ylim(0, 225)+
  ylab("Height (cm)")+
  scale_color_manual( values = c("darkgrey", "#ffd16b", "#80e6ff" ))+
  scale_fill_manual( values = c("darkgrey", "#ffd16b", "#80e6ff"))+
  scale_x_discrete(labels =
                     c(
                       bquote(italic("A. femoralis")),
                       bquote(italic("O. oophagus")),
                       bquote(italic("D. tinctorius"))))+
  theme(legend.position = "none")+
  labs(tag = "B")

b1

a1b1<- grid.arrange(a1,legend1, b1, ncol = 2, widths = c(2.5,0.45))
a1b1

#FIGURE 4----

#Terrestrial

b<- ggplot(pool_terr, 
           aes(water_capacity, as.numeric(Tad_Y_N)-1, color = species, shape = species)) +
  geom_point(position=position_jitter(height=0.03, width=0), size = 4, alpha = 0.6) +
  xlab("Water-holding capacity (mL)") + 
  ylab("Probability (Pool occupancy)")+
  facet_wrap(~ Water_Cap1, scale= "free_x")+
  scale_color_manual( values = c("darkgrey", "#ffd16b", "#80e6ff", "black"), 
                      name = "Species", 
                      labels = c(bquote(italic("A. femoralis")),
                                 bquote(italic("O. oophagus")),
                                 bquote(italic("D. tinctorius")),
                                 "None"))+   
  scale_shape_manual(values = c(16, 16, 16, 1), 
                     name = "Species",
                     labels = c(bquote(italic("A. femoralis")),
                                bquote(italic("O. oophagus")),
                                bquote(italic("D. tinctorius")),
                                "None"))+
  theme_bw(base_line_size = 0)+
  theme(panel.spacing.x = unit(9, "mm"),
        text = element_text(size=13))+
  labs(tag = "A")

b 


legend<- get_legend(b)

b<- b + theme(legend.position = "none")

c<- ggplot(pool_terr, 
           aes(x= leaflit_vol_ml, y = height, color =species, shape = species))+
  geom_point(alpha = 0.6, size = 4)+
  xlab("Leaf litter (vol)")+
  theme_bw(base_line_size = 0)+ 
  theme(legend.position = "none",
        text = element_text(size=13))+
  geom_hline(yintercept = 75, linetype = "dashed", color = "darkgrey")+
  ylab("Height (cm)")+
  facet_wrap(~ Water_Cap1)+
  scale_color_manual( values = c("darkgrey", "#ffd16b", "#80e6ff", "black"))+
  scale_shape_manual(values = c(16, 16, 16, 1))+
  labs(tag = "B")

c

grid.arrange(b,legend, c, ncol = 2, widths = c(2.5,0.45))

#FIGURE 5----

kh_high<-
  ggplot(frog1, aes(y = height, x = KH))+
  stat_smooth(method = "glm", formula = y~ (x), aes(x= KH, y = height),
              se = T, fill = "lightgrey", color = "black", fullrange = F)+
  geom_jitter(alpha = 0.8, size = 4, 
              height = 0, width = 0, aes(shape = species, color = species))+
  theme_bw(base_line_size = 0)+ 
  theme( text = element_text(size=13), legend.position = "none")+
  xlab("KH")+
  ylab("Height (cm)")+
  geom_vline( xintercept = 6.2, linetype = "dashed",
              color = "#ffd16b", size = 1)+
  geom_vline( xintercept = 6.4, linetype = "dashed",
              color = "darkgrey", size = 1)+
  geom_vline( xintercept = 20, linetype = "dashed",
              color = "#80e6ff", size = 1)+
  coord_cartesian(ylim = c(0, 2100))+
  scale_color_manual( values = c("darkgrey", "#ffd16b", "#80e6ff", "black"))+
  scale_shape_manual( values = c(16, 16,16, 1))+
  labs(tag = "A")

kh_high

leaf_pred<- 
  ggplot(frog1, aes(x= leaflit_vol_ml, y = predator_count))+
  stat_smooth(formula = y~ (x), method = "glm", aes(x= leaflit_vol_ml, y = predator_count),
              se = T, fill = "lightgrey", color = "black", fullrange = T, 
              size = 0.8)+
  geom_jitter(alpha = 0.7, size = 4, height = 0.3, width = 0,
              aes(shape = species, color = species))+
  coord_cartesian(xlim = c(0, 2500))+
  theme_bw(base_line_size = 0)+ 
  theme(text = element_text(size=13), legend.position = "none")+
  geom_vline( xintercept = 1200, linetype = "dashed",
              color = "#80e6ff", size = 1)+
  geom_vline( xintercept = 50, linetype = "dashed",
              color = "#ffd16b", size = 1)+
  geom_vline( xintercept = 2500, linetype = "dashed",
              color = "darkgrey", size = 1)+
  xlab("Leaf litter (vol/mL)")+
  ylab("Predator Count")+
  scale_color_manual( values = c("darkgrey", "#ffd16b", "#80e6ff", "black"))+
  scale_shape_manual( values = c(16, 16,16, 1))+
  labs(tag = "B")

leaf_pred

sal_hard<-
  ggplot(frog1, aes(x= salinity, y = hardness))+
  theme_bw(base_line_size = 0)+ 
  theme(text = element_text(size=13))+
  coord_cartesian(ylim = c(0, 25)) +
  xlab("Salinity (ppm)")+
  ylab("Hardness")+
  stat_smooth(method = "glm", formula = y ~ (x), aes(x= salinity, y = hardness),
              fill = "lightgrey", color = "black")+
  geom_jitter(alpha = 0.8, size = 4, height = 0.2, width = 0,
              aes(shape = species, color = species))+
  geom_vline( xintercept = 240, linetype = "dashed",
              color = "#ffd16b", size = 1)+
  geom_vline( xintercept = 53, linetype = "dashed",
              color = "darkgrey", size = 1)+
  geom_vline( xintercept = 955, linetype = "dashed",
              color = "#80e6ff", size = 1)+
  scale_color_manual( values = c("darkgrey", "#ffd16b", "#80e6ff", "black"), 
                      name = "Species", 
                      labels = c(bquote(italic("A. femoralis")),
                                 bquote(italic("O. oophagus")),
                                 bquote(italic("D. tinctorius")),
                                 "None"))+
  scale_shape_manual( values = c(16, 16,16, 1), name = "Species",
                      labels = c(bquote(italic("A. femoralis")),
                                 bquote(italic("O. oophagus")),
                                 bquote(italic("D. tinctorius")),
                                 "None"))+
  labs(tag = "C")

sal_hard

top<- grid.arrange(kh_high, leaf_pred, ncol = 2)

grid.arrange(top, sal_hard)

#FIGURE 6----
#See PCA analysis

#FIGURE 7----

#ph/time graph ----

coeff1<-10

ph<-
  ggplot(subset(chloe_rain_all, pool_type %in% c("alive", "dead")), aes(x= Week)) +
  geom_point(aes(y = pH, color = pool_type), size =3, alpha = 0.9, shape = 1)+
  xlim (1,4)+
  stat_smooth(aes(y = pH, color = pool_type), 
              method = "loess", formula = y ~ log(x), se = F, size = 2)+
  geom_line(aes(y = week_total_rain/coeff1), color = "lightblue", 
            linetype = "longdash", size = 2)+
  geom_point(aes(y = week_total_rain/coeff1), color = "blue", size = 3)+
  scale_y_continuous(name = "pH",
                     sec.axis = sec_axis(~.*coeff1, name = "Cumulative rain (mm/week)"))+
  theme_bw(base_line_size = 0)+
  theme(text = element_text(size=13), 
        legend.position = "bottom",
        axis.title.y.right =  element_text(margin = margin(t = 0, r = 0, b = 0, l = 15)))+
  scale_color_manual( values = c("#00b606", "#e6a81e"), name = "Pool substrate",
                      labels = c("Live", "Dead"))+
  labs(tag ="A")

ph


pH_Tad<-
  ggplot(subset(pool_2020, pool_type %in% c("alive", "dead")), 
         aes( x = pH, y = Dt_Tadpole_Num, color = pool_type))+
  geom_point(size =4, shape = 1)+
  stat_smooth( method = "glm", se = T, 
               formula = y ~ log(x), fullrange = T, 
               size = 2, fill = "#f7f7f7")+
  stat_smooth(aes(color = pool_type),
              method = "glm", se = T, 
              formula = y ~ log(x), fill = "transparent", geom = "ribbon",
              linetype = "dashed")+
  ylim(0, 22)+
  ylab("Tadpole count")+
  theme_bw(base_line_size = 0)+
  scale_color_manual(values = c("#00b606", "#e6a81e" ),
                     name = "Pool Substrate", 
                     labels = c("Live", "Dead"))+
  theme( text = element_text(size=13), 
         legend.position = "none")+
  labs(tag ="B")

pH_Tad
