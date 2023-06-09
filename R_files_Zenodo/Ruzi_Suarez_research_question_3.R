## Assessing the amount of damamge E. ruidum workers have on Zanthoxylum
## ekmanii seeds, and if they have the bite force needed to crush the seed coat.

## RQ3: To what extent do E. ruidum damage the seeds that they disperse?

# abbreviations: Zan = Zanthoxylum ekmanii; ER = Ectatomma ruidum

# -- load libraries ####
library(here)
library(tidyverse) # includes ggplot2

# -- set paths ####
data_path <- here::here("data")
figure_path <- here::here("figs")

# -- read in the data ####

processing_data <- read_csv(paste(data_path, "Ruzi_Suarez_data_seed_surface_area.csv", sep = "/"),
         col_types = 
           cols(
             Seed_species = col_character(),
             file_name = col_character(),
             Control_Colony = col_character(),
             groupID = col_double(),
             Chamber = col_character(),
             total_SA = col_double(),
             Darea_1 = col_double(),
             D_area2 = col_double(),
             D_area3 = col_double(),
             D_area4 = col_double(),
             D_area5 = col_double(),
             Included_good_areas = col_double(),
             Sum_damage = col_double(),
             Seed_damage_percent = col_double(),
             notes1 = col_character(),
             notes2 = col_character()))
names(processing_data)

## - t-test of damage between the two conditions ####

control_data <- processing_data %>%
  filter(Control_Colony == "Control")
nrow(control_data) # 41

cont_seeds_percent <- control_data$Seed_damage_percent


colony_data <- processing_data %>%
  filter(Control_Colony != "Control")
nrow(colony_data) #36

colony_seeds_percent <- colony_data$Seed_damage_percent

# conducting a t-test between the two
model1 <- t.test(cont_seeds_percent, colony_seeds_percent) #leaving var.equal as false
model1
#Welch Two Sample t-test
#data:  cont_seeds_percent and colony_seeds_percent
#t = -5.0415, df = 48.902, p-value = 6.736e-06
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -25.33157 -10.89195
#sample estimates:
#  mean of x mean of y 
#8.428073 26.539833 

6.736e-06 < 0.05 # TRUE so there is a significant difference in the amount of "processing" damage
6.736e-06 < 0.001 # TRUE


## - analyze ER bite force data ####
# bite force data obtained from Michael Rivera and reported in Table 2 of the manuscript
# main text

max_ER_bites_N <- c(0.078, 0.060, 0.092, 0.092)
shapiro.test(max_ER_bites_N)
#Shapiro-Wilk normality test
#data:  max_ER_bites_N
#W = 0.85603, p-value = 0.2463
#null hypothesis is that it is normally distributed

## comapring the force needed to rupture the seed reported in Zalamea et al. 2018
# 68.3 is the force reported in Zalamea et al. 2018

t.test(max_ER_bites_N, mu=68.3, alternative="two.sided") #though visually it is so much smaller
#mu = the hypothetical mean so taking that from Camilo's paper
#One Sample t-test
#data:  max_ER_bites_N
#t = -8990, df = 3, p-value = 3.035e-12
#alternative hypothesis: true mean is not equal to 68.3
#95 percent confidence interval:
#  0.05635043 0.10464957
#sample estimates:
#  mean of x 
#0.0805 
3.035e-12 < 0.001 # TRUE

mean(max_ER_bites_N)
68.3/mean(max_ER_bites_N) # 848.4472 x greater than the mean ER worker max bite force

## - to graph the seed processing data ####

# neet to set up a column that has either Control or Colony
processing_data2 <- processing_data %>%
  mutate(Grouping = Control_Colony) %>%
  mutate(Grouping = if_else(Control_Colony == "Control", "Control", Grouping)) %>%
  mutate(Grouping = if_else(Control_Colony != "Control", "Colony", Grouping))

seed_processing_plot <- ggplot(processing_data2, aes(x=Grouping, y=Seed_damage_percent))+
  geom_boxplot(fill="grey40")+
  ylab("Seed Surface Area Damaged (%)")+
  scale_y_continuous(limits=c(0,110), expand=c(0,0))+
  xlab("Treatment")+
  theme(axis.line=element_line(colour="black"),
        strip.text.y=element_text(size=16, face="italic", angle=0),
        strip.background=element_blank(),
        strip.text.x=element_text(face="bold", colour="black", angle=45, hjust=1, vjust=1, size=16),
        axis.text.y=element_text(colour="black", size=12),
        axis.title.x=element_text(face="bold", colour="black", size=12,margin = margin(20,0,0,0)),
        axis.title.y=element_text(face="bold", colour="black", size=12, margin=margin(0,20,0,0)),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        legend.title=element_text(colour="black", face="bold", size=12),
        legend.text=element_text(colour="black", size=12),
        #plot.title=element_text(family="Arial", face="bold", colour="black", size=rel(1.0)),
        panel.border=element_blank(),
        #legend.key=element_rect(colour="black", size=1),
        legend.position="top",
        panel.background=element_blank(),
        axis.line.y	=element_line(colour="black"),
        axis.line.x=element_line(colour="black"))
seed_processing_plot

## - save the figures ####

ggsave("Ruzi_Suarez_Figure3a.png", width = 4, height = 4,
       units = "in", dpi = 600, plot = seed_processing_plot,
       path = figure_path)

