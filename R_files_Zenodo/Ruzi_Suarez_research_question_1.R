## Describing the number of Zanthoylum ekmanii seeds moved by Ectatomma ruidum,
## how far they are moved, where they are taken, and if the dispersal distance is
## great enough to deposit the seed outside of the tree crown.

## RQ1: Do E. ruidum workers conduct directed dispersal depositing seeds in their
## colonies at distances that are biologically relevant (i.e. potentially far enough
## away to escape negative density dependent drivers of mortality)?

# abbreviations: Zan = Zanthoxylum ekmanii; ER = Ectatomma ruidum

# -- load libraries ####
library(here)
library(tidyverse) # includes ggplot2
library(e1071)
library(plyr)

# -- set paths ####
data_path <- here::here("data")
figure_path <- here::here("figs")

# -- read in the data ####

movement_data <- read_csv(paste(data_path, "Ruzi_Suarez_data_Ectatomma_ruidum_distance_deposition.csv", sep = "/"),
         col_types = 
           cols(
             site = col_character(),
             rep_num = col_double(),
             seed_color = col_character(),
             distance_cm = col_double(),
             destination = col_character(),
             ant_species = col_character(),
             num_foragers = col_character(),
             movement_num = col_character(),
             year = col_double(),
             ID = col_character()))
movement_data
names(movement_data)

overall_data <- read_csv(paste(data_path, "Ruzi_Suarez_data_overall_seed_removal_2hours.csv", sep = "/"),
         col_types = 
           cols(
             Plot = col_character(),
             Year = col_double(),
             Rep = col_double(),
             Side = col_character(),
             Seeds_placed = col_double(),
             Seeds_removed = col_double(),
             Seeds_remaining = col_double(),
             Date = col_character()))
nrow(overall_data)

## - to get the total number of seeds placed over the two years ####

5*5*3 # 5 sites by 5 seeds by 3 replicates for 2015 = 75 seeds
5*5*2 # 5 sites by 5 seeds by 2 replicates for 2016 = 50 seeds
75+50 #total number of seeds placed then = 125 seeds 

## - to get the total number of seeds moved after 2 hours ####
names(overall_data)
# these could have been moved by workers of any ant species over the 2 hour duration of observation

removal_data <- overall_data %>%
  # get a proportin of seeds removed to seeds placed
  mutate(prop = Seeds_removed/Seeds_placed) %>%
  # to turn that proportion into a percent
  mutate(percent = prop*100)
removal_data

se <- function(x) sqrt(var(x)/length(x))

removal_summary <- ddply(removal_data, "Year", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent,na.rm=TRUE),
                         se=se(percent),
                         total_removed=sum(Seeds_removed),
                         total_seeds_placed=sum(Seeds_placed))
removal_summary
# Year  n     Mean       se total_removed total_seeds_placed
#1 2015 15 58.66667 11.45869            44                 75
#2 2016 10 58.00000 13.80821            29                 50
 
# only needed plyr package for the above code, detach now so that will not
# interfere with remaining code
detach("package:plyr", unload = TRUE)

## - to get the number of seeds moved once, twice, and three times by ER ####
movement_data %>%
  count(movement_num)
# A tibble: 3 x 2
#  movement_num       n
#    <chr>        <int>
#  1 one             41
#  2 three            1
#  3 two              9

# percentation of seeds moved each distance
41/125 # 0.328
9/125 # 0.072
1/125 # 0.008

## - first movement and destinations ####

first_movement <- movement_data %>%
  filter(movement_num == "one")
nrow(first_movement) # 41, matches the number of seeds moved first distance by ER

# will set na.rm = TRUE to calculate mean and sd while removing seeds that were lost 
# during tracking and therefore do not have distances measured

mean(first_movement$distance_cm, na.rm = TRUE) # 99.73846
sd(first_movement$distance_cm, na.rm = TRUE) # 91.13132

# to get the maximum first straight line distance
max(first_movement$distance_cm, na.rm = TRUE) # 405.7

# to get how many seeds were lost during their first movement
first_movement %>%
  filter(is.na(distance_cm)) 
# 2 rows returned therefore 2 seeds were lost during tacking
# of the first movement distances, and 39 (41-2) have a first movement distance
# and were not lost

# to get the locations of where first movements took seeds
#counting dropped, ground, and leaf litter all as "dropped"
first_movement %>%
  count(destination)
# A tibble: 5 x 2
#   destination     n
#   <chr>       <int>
#1 colony         28
#2 dropped         5
#3 ground          2
#4 leaf litter     2
#5 lost            4

## - second movement and destinations ####

second_movement <- movement_data %>%
  filter(movement_num == "two")
nrow(second_movement) # 9, matches the number of seeds moved second distance by ER

# will set na.rm = TRUE to calculate mean and sd while removing seeds that were lost 
# during tracking and therefore do not have distances measured

mean(second_movement$distance_cm, na.rm = TRUE) # 61.50556
sd(second_movement$distance_cm, na.rm = TRUE) # 70.92034
max(second_movement$distance_cm, na.rm = TRUE) # 208.5

# to check if any seeds were lost during second movement
second_movement %>%
  filter(is.na(distance_cm))
# returns 0 rows therefore none were lost

# the destination of the seeds moved a second distance by ER
#counting dropped, ground, and leaf litter all as "dropped"
second_movement %>%
  count(destination)
# A tibble: 4 x 2
#   destination     n
#   <chr>       <int>
#1 colony          4
#2 dropped         1
#3 ground          3
#4 leaf litter     1

## - third distance movement and destination ####

third_movement <- movement_data %>%
  filter(movement_num == "three")
nrow(third_movement) # 1 row

third_movement$distance_cm # 4
third_movement$destination # "ground"

## - to describe the distance distributions ####
# using the e1071 package

## first movement

shapiro.test(first_movement$distance_cm)
#Shapiro-Wilk normality test
#data:  first_movement$Distance_cm
#W = 0.81892, p-value = 2.113e-05

2.113e-05 < 0.001 # TRUE, therefore data is not normal

skewness(first_movement$distance_cm, na.rm=TRUE)
#1.560971


kurtosis(first_movement$distance_cm, na.rm=TRUE)
#2.157506

## second movement

nrow(second_movement)
shapiro.test(second_movement$distance_cm)
#Shapiro-Wilk normality test
#data:  second_movement$Distance_cm
#W = 0.80672, p-value = 0.02436
# also not normally distributed

skewness(second_movement$distance_cm, na.rm=TRUE)
#1.015333


kurtosis(second_movement$distance_cm, na.rm=TRUE)
# -0.5660578

## - to compare the movement distance to tree crown ####
# using the crown area of Zan from Park et al. 2019 of 211 m2 (Â± 38 SD)
# and assuming the crown is in the shape of a circle, will use the area of a cirlce
# formula to calculate the radius, which would indicate the middle of the crown, where the trunk would be

# area of a circle: A = pi(r^2)
# rearranged to get r: r = sqrt(A/pi)

#the area
radius <- sqrt(211/pi) # 8.195327

# the SD
radius_sd <- sqrt(38/pi) # 3.477898

# as the distances were not normally distributed, will conduct a non-parametric test
# conducting the test on the first movement distances only
# need to remove na distances first
first_movement_na_rm <- first_movement %>%
  filter(!is.na(distance_cm))
first_movement_distances <- first_movement_na_rm$distance_cm
class(first_movement_distances)
wilcox.test(first_movement_distances, mu = radius)
# Wilcoxon signed rank test with continuity correction
# 
# data:  first_movement_distances
# V = 766, p-value = 1.571e-07
# alternative hypothesis: true location is not equal to 8.195327
# 
# Warning message:
#   In wilcox.test.default(first_movement_distances, mu = radius) :
#   cannot compute exact p-value with ties

1.571e-07 < 0.001 # TRUE


## - to make the distance figures ####
names(movement_data)

movement_data2 <- movement_data %>%
  # remove the missing distances from lost seeds
  filter(!is.na(distance_cm))

distances_figure <- ggplot(movement_data2, aes(x=distance_cm, fill=movement_num, color=movement_num))+ #color=Species, fill=Species
  geom_histogram(position="identity", binwidth=15, alpha=0.6)+ #alpha=0.6
  scale_x_continuous(breaks=seq(0,405,15), expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab("Distance (cm)")+
  ylab("Count")+
  scale_fill_manual(breaks=c("one","two","three"),
                    labels=c("One","Two","Three"),
                    values=c("grey70", "grey40","black"),
                    name="Movement\\nNumber")+
  scale_color_manual(breaks=c("one","two","three"),
                     labels=c("One","Two","Three"),
                     values=c("grey60", "grey40","black"),
                     name="Movement\\nNumber")+
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_blank(), #gets rid of square going around the entire graph
        axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="top",
        legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(face="bold", size=12, color="black"), #size of x-axis title
        axis.title.y = element_text(face="bold", size=12, color="black"), #size of y-axis title
        axis.text.x = element_text(size=11, color="black", angle=45, vjust=1, hjust=1), #size of x-axis text
        axis.text.y = element_text(size=12, color="black"),
        axis.line.x	=element_line(colour="black"),
        axis.line.y	=element_line(colour="black"))#size of y-axis text
distances_figure 

## - desposition figure ####
# first need to clean up some of the location categories

movement_data %>%
  count(destination)

nrow(movement_data)# 51 destination locations total

movement_data3 <- movement_data %>%
  # going to clean up some of the destination locations for easier graphing
  # will merge the dropped, ground, and leaf litter to the dropped category
  mutate(destination2 = if_else(destination == "colony", "Colony", destination)) %>%
  mutate(destination2 = if_else(destination == "lost", "Lost", destination2)) %>%
  mutate(destination2 = if_else(destination2 != "Colony" & destination2 != "Lost",
                                "Dropped", destination2)) %>%
  # fix the movement_num values so that they are capitalized
  mutate(movement_num = if_else(movement_num == "one", "One", movement_num)) %>%
  mutate(movement_num = if_else(movement_num == "two", "Two", movement_num)) %>%
  mutate(movement_num = if_else(movement_num == "three", "Three", movement_num)) %>%
  # to get number of seeds at each destination by movement number
  count(movement_num, destination2) %>%
  # to get the percentage of seeds at that destination based on all destination
  # locations across all 51 recorded movements and destinations
  mutate(percent = n/51*100) %>% # 100 to change from frequency to a percent
  # plotting with just this information will have the bars take up the full space 
  # if categories with 0 values are not included in the data, therefore need to add 
  # some 0 destination locations to some of the movement numbers
  add_row(movement_num = "Two", destination2 = "Lost", n = 0, percent = 0) %>%
  add_row(movement_num = "Three", destination2 = "Colony", n = 0, percent = 0) %>%
  add_row(movement_num = "Three", destination2 = "Lost", n = 0, percent = 0)

movement_data3
# movement_num destination2     n percent
# <chr>        <chr>        <dbl>   <dbl>
# 1 One          Colony          28   54.9 
# 2 One          Dropped          9   17.6 
# 3 One          Lost             4    7.84
# 4 Three        Dropped          1    1.96
# 5 Two          Colony           4    7.84
# 6 Two          Dropped          5    9.80
# 7 Two          Lost             0    0   
# 8 Three        Colony           0    0   
# 9 Three        Lost             0    0

28 + 4 # 32 into colony
9 + 1 + 5 # 15 dropped
4 # 4 lost


movement_data3$movement_num <- factor(movement_data3$movement_num, levels=c("One","Two","Three"))
levels(movement_data3$movement_num) # to change the order of the levels
movement_data3

destinations<- ggplot(movement_data3, aes(x=destination2, y=percent, fill=movement_num))+
  geom_bar(stat="identity", position="dodge",colour="black")+
  scale_y_continuous(limits=c(0,100), expand=c(0,0))+ 
  #guides(fill=guide_legend(title="Movement Number"), override.aes=list(colour=NULL))+
  #scale_fill_manual(values=c("blue","lightblue"), breaks=c("2015.one","2015.two","2016.one","2016.two"),
  #                 labels=c("2015\\nFirst","2015\\nSecond", "2016\\nFirst", "2016\\nSecond"))+
  #scale_fill_manual(values=c("thistle4", "thistle"), breaks=c("one","two"), labels=c("First", "Second")) + #guide=FALSE, not working
  ylab("Frequency (%)")+
  xlab("Destination")+
  scale_fill_manual(breaks=c("One","Two","Three"),
                    values=c("grey70", "grey40","black"),
                    name="Movement\\nNumber")+
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_blank(), #gets rid of square going around the entire graph
        axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="top",
        legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(face="bold", size=12, color="black"), #size of x-axis title
        axis.title.y = element_text(face="bold", size=12, color="black"), #size of y-axis title
        axis.text.x = element_text(size=11, color="black"), #size of x-axis text
        axis.text.y = element_text(size=12, color="black"),
        axis.line.x	=element_line(colour="black"),
        axis.line.y	=element_line(colour="black"))#size of y-axis text
destinations

## - save the figures ####


ggsave("Ruzi_Suarez_Figure1b.png", width = 12, height = 4,
       units = "in", dpi = 600, plot = distances_figure,
       path = figure_path)

ggsave("Ruzi_Suarez_Figure1c.png", width = 3.5, height = 4,
       units = "in", dpi = 600, plot = destinations,
       path = figure_path)
