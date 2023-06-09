#Load in necessary packages
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(wesanderson)
library(maps)
library(rgeos)
library(rworldmap)
library(mapproj)
library(viridis)

setwd("set your working directory here")

# Read in the most recent version of the data 

IR_GHG <- read.csv("Data_Extraction_Refs_IR_GHG_2022-01-05.csv", skip=1) #for IRES studies
str(IR_GHG) #43 obs. of  33 variables

Dam_GHG <- read.csv("Data_Extraction_Refs_Dam_GHG_2022-01-05.csv",  skip=1) #for damming studies
str(Dam_GHG) #49 obs. of  40 variables

#Plot the count of number of publications per year
Year_IR <- ggplot(IR_GHG, aes(Publication.Year)) + geom_bar() + theme_pubr() 
Year_IR

Year_Dam <- ggplot(Dam_GHG, aes(Publication.Year)) + geom_bar() + theme_pubr()
Year_Dam

#In order to have both on one figure try merging

#First subset ID and Year, then make a new column for Topic = either GHG or IR. Then merge by Topic... 
IR_GHG_sub <- IR_GHG[c(1,4)]
Dam_GHG_sub <- Dam_GHG[c(1,4)]
IR_GHG_sub$Topic <- c("IRES")
Dam_GHG_sub$Topic <- c("Dam")

combine <- rbind(Dam_GHG_sub, IR_GHG_sub)

#Now plot count of publications per year by group
#Plot the count of number of publications per year
Year_comb <- ggplot(combine, aes(Publication.Year, fill= factor(Topic))) + geom_bar() + theme_pubr() 
Year_comb #this is a little difficult to decipher with the overlap, can we make the bars side by side?

Counts <- combine %>%
  filter(Topic %in% c("Dam", "IRES")) %>%
  group_by(Publication.Year, Topic) %>%
  summarise(counts = n()) 
head(Counts, 4)
Counts

#Remove any NA's

Counts <- Counts %>% drop_na()



setEPS()
postscript("Fig2.eps")


Year_comb2 <- p <- ggplot(Counts, aes(x = Publication.Year, y = counts)) +
  geom_bar( 
    aes(color = Topic, fill = Topic),
    stat = "identity", 
    width = 0.8, 
    position = position_dodge2(preserve = "single")
  ) + theme(panel.background = element_blank()) + theme(axis.line = element_line(colour = "black")) + scale_y_continuous(expand = c(0, 0)) +  theme(legend.title = element_blank(), text = element_text(size = 12))+  labs(y= "Count", x = "Publication Year") +
  scale_fill_manual(values = c("#79969B","#C27D38" )) +  
  scale_colour_manual(values = c("#79969B","#C27D38"))
Year_comb2


ggsave("Fig2.eps", plot = Year_comb2, width = 5.5,height = 3, units = c("in") ) 


#ggsave("Publication_Year_by_Topic2.png", plot = Year_comb2, width = 5.5,height = 3, units = c("in") )

#####################################################################################################
## Make 2 figures of GHG and damming studies by country ## Use a bubble map

## Summarize number of studies per country data
## IRES studies

# Subset ID, country, Study type

IRES_Country <- IR_GHG[c(1,7,12)]

# Exclude review and model studies (only include field or mesocosm studies)
names(IRES_Country)[names(IRES_Country) == "Study.type..Field..lab..both..review..conceptual."] <- "Study_type"

IRES_Country <- IRES_Country[!grepl("Review", IRES_Country$Study_type),]
IRES_Country <- IRES_Country[!grepl("Model", IRES_Country$Study_type),]

#Remove NA's
IRES_Country <-IRES_Country %>% drop_na()

# Since we are interested in data by country, change any of the states to USA

IRES_Country$Country[IRES_Country$Country == "Rhode Island, USA"] <- "USA"

IRES_Country$Country[IRES_Country$Country == "Arizona, USA"] <- "USA"

IRES_Country$Country[IRES_Country$Country == "Indiana, USA"] <- "USA"


IRES_Country$Country[IRES_Country$Country == "Global (29 countries)"] <- "Global"

IRES_Country$Country[IRES_Country$Country == "Global (22 countries)"] <- "Global"

IRES_Country$Country[IRES_Country$Country == "Global (17 countries, all continents except Antarctica)"] <- "Global"

IRES_Country$Country[IRES_Country$Country == "Global (89 independent systems)"] <- "Global"

# Make a count of each country 

# Make country a factor
IRES_Country$Country <- as.factor(IRES_Country$Country)

IRES_Country_count <-  IRES_Country %>%
  group_by(Country) %>%
  summarise(counts = n()) 

IRES_Country_count

as.data.frame(IRES_Country_count)
str(IRES_Country_count) #Country: Factor w/ 8 levels


#Load centroid long and lat coordinates

country_centroids <- read.csv("Country_centroids.csv")

country_centroids$name[country_centroids$name == "United States"] <- "USA" #Change to match data.frame

#Merge the long and lat centroid columns with IRES_Country_count
IRES_Country_data <- merge(IRES_Country_count, country_centroids, by.x = "Country", by.y = "name")

#######################################
#Damming studies

Dam_Country <- Dam_GHG[c(1,7,12)]

# Exclude review and model studies (only include field or mesocosm studies)
names(Dam_Country)[names(Dam_Country) == "Study.type..Field..lab..both..review..conceptual."] <- "Study_type"


Dam_Country <- Dam_Country[!grepl("Review", Dam_Country$Study_type),]
Dam_Country <- Dam_Country[!grepl("Model", Dam_Country$Study_type),]

#Remove NA's
Dam_Country <-Dam_Country %>% drop_na()

# Since we are interested in data by country, change any as necessary

Dam_Country$Country[Dam_Country$Country == "Quebec, Canada"] <- "Canada"

Dam_Country$Country[Dam_Country$Country == "Ohio, USA"] <- "USA"

Dam_Country$Country[Dam_Country$Country == "Delaware, USA"] <- "USA"

Dam_Country$Country[Dam_Country$Country == "French"] <- "France"

Dam_Country$Country[Dam_Country$Country == "Africa"] <- "Zambia"

#If one study has two countries, add an additional entry with the second country
# DAM_GHG_247 Bolivia/Peru Field-estimations

Dam_Country$Country[Dam_Country$Country == "Bolivia/Peru"] <- "Bolivia"

#DAM_GHG_315 French Guiana and Panama Field
Dam_Country$Country[Dam_Country$Country == "French Guiana and Panama"] <- "French Guiana"

# DAM_GHG_408 French Guiana and Brazil Field
Dam_Country$Country[Dam_Country$Country == "French Guiana and Brazil"] <- "French Guiana"

ID <- c("DAM_GHG_247b", "DAM_GHG_315b", "DAM_GHG_408b")
Country <- c("Peru", "Panama", "Brazil")
Study_type <- c("Field-estimations", "Field", "Field")
Duplicates <- data.frame(ID, Country, Study_type) 

Dam_Country <- rbind(Dam_Country, Duplicates) 

# Make country a factor
Dam_Country$Country <- as.factor(Dam_Country$Country)

Dam_Country_count <-  Dam_Country %>%
  group_by(Country) %>%
  summarise(counts = n()) 

Dam_Country_count

as.data.frame(Dam_Country_count)
str(Dam_Country_count) #Country: Factor w/ 17 levels

Dam_Country_data <- merge(Dam_Country_count, country_centroids, by.x = "Country", by.y = "name")


#########################################
# Load world map

world <- map_data("world")

#Now make map with points

mybreaks <- c(1,2,5, 10, 13) 

IRES_Map <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey", alpha=0.3) +
  geom_point( data=IRES_Country_data, aes(x=longitude, y=latitude, size=counts), colour="red", alpha=0.5, shape=16) +
  scale_size_continuous(name="IRES Studies",  range=c(1,13) , breaks=mybreaks)+
 theme_void() +  coord_map(xlim=c(-180,180)) 
IRES_Map

Dam_Map <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey", alpha=0.3) +
  geom_point( data=Dam_Country_data, aes(x=longitude, y=latitude, size=counts), colour="blue", alpha=0.5, shape=16) +
  scale_size_continuous(name="Dam Studies")+
  theme_void() +  coord_map(xlim=c(-180,180)) 
Dam_Map

##################
## Now create a new data frame with a group column for IRES and Dam categories so you can plot the data on a single map
IRES_Country_data$Topic <- c("IRES")
Dam_Country_data$Topic <- c("Dam")

all_country_data <- rbind(IRES_Country_data, Dam_Country_data) 

mybreaks2 <- c(1,4, 8,16) 

Bubble_Map <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey", alpha=0.3) +
  geom_point( data=all_country_data, aes(x=longitude, y=latitude, size=counts, colour=Topic), alpha=0.6, shape=16) +
  scale_size_continuous(name="Studies", trans="log", range=c(3,20), breaks=mybreaks2)+
  theme_void() +  coord_map(xlim=c(-180,180)) +
  scale_color_manual(values=c("#79969B", "#C27D38")) + 
  theme(
    legend.position = c(0.1, 0.42),
    text = element_text(color = "#22211d", size=20),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    ) + guides(shape = guide_legend(override.aes = list(size = 1)))

Bubble_Map

#can't change the size of legend bubbles here, so I must do it manually in Adobe Illustrator

ggsave("Global_distribution_of_GHG_studies.png", plot = Bubble_Map , units = c("in"), width = 8, height = 7 ) #5.42 x 4.41


#Subset only the IRES study for presentation

IRES_countries <- subset(all_country_data, Topic == "IRES")
mybreaks3 <- c(1,6,12) 

Bubble_Map2 <- ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey", alpha=0.3) +
  geom_point( data=IRES_countries, aes(x=longitude, y=latitude, size=counts, colour=Topic), alpha=0.6, shape=16) +
  scale_size_continuous(name="Number of studies", trans="log", range=c(3,20), breaks=mybreaks3)+
  theme_void() +  coord_map(xlim=c(-180,180)) +
  scale_color_manual(values=c( "#C27D38")) + 
  theme(
    legend.position = c(0.1, 0.42),
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
  ) + guides(shape = guide_legend(override.aes = list(size = 0.5)))

Bubble_Map2

ggsave("Global_distribution_of_IRES_GHG_studies.png", plot = Bubble_Map2 , units = c("in"), width = 8, height = 7 ) #5.42 x 4.41


##########################################
## Make a pie chart of the number of studies per GHG
## Or to make the values relative to each type of study, could we do some kind of bubble chart of bar graph?

GHGs <-  read.csv("type_of_GHGs_per_study.csv")

## Or to make the values relative to each type of study, could we do some kind of bubble chart of bar graph?
#try using a balloon plot


tiff("GHG_per_study_type.tiff", units="in", width=4, height=4, res=1200)


GHG_plot <- ggplot(GHGs, aes(x = Topic, y = Gas, colour=Topic)) +
  geom_point(aes(size = Count), shape = 16, alpha=0.7) +
  scale_size_area(max_size = 40, guide = FALSE) +
  #geom_text(aes(
   # y = as.numeric(as.factor(Gas)) - sqrt(Count)/35, label = Count),
   # vjust = -0.5,
   # colour = "black",
   # size = 4) + 
  theme_classic() +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.line=element_blank(),axis.ticks=element_blank(), panel.border = element_blank(), panel.grid = element_blank()) +
  scale_color_manual(values=c("#79969B", "#C27D38")) 
GHG_plot

# insert ggplot code
dev.off()

ggsave("GHG_per_study_type.png", plot = GHG_plot , units = c("in"), width = 4, height = 4 ) #5.42 x 4.41


#Now make a similar plot that can be beside (e.g. Figure 4 A and B) except with the sclaes: global, network, local

#subset the scale column for IR

IR_scale_sub <- IR_GHG[c(1,12,14)]

# Exclude review, metabolism and model studies (only include field or mesocosm studies)
names(IR_scale_sub)[names(IR_scale_sub) == "Study.type..Field..lab..both..review..conceptual."] <- "Study_type"

IR_scale_sub <- IR_scale_sub[!grepl("Review", IR_scale_sub$Study_type),]
IR_scale_sub <- IR_scale_sub[!grepl("Model", IR_scale_sub$Study_type),]
IR_scale_sub <- IR_scale_sub[!grepl("metabolism", IR_scale_sub$Study_type),]


#Clean entries with extra words

IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (10 temporary tributaries)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (~40 km section of a river)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (100m reach)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Network (12 intermittent streams 30-m study reaches)"] <- "Network"
IR_scale_sub$Scale[IR_scale_sub$Scale == "2 networks Fluvia` and Muga rivers (19 sites)"] <- "Network"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Global (200 dry IRES reaches)"] <- "Global"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Global (212 stream reaches)"] <- "Global"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (33 stream reaches)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (8 stream reaches, across 2 regions)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Catchment (multiple sites along sub-catchment)"] <- "Network"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Mesocosm (16 artificial streams)"] <- "NA"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (perennial, 345 m long)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Global (196 dry inland water reaches and adjacent upland)"] <- "Global"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (one pool at a reach)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (1 intermittent and 1 perennial reach, 50m each)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (6 sites along a 345 m reach)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Two temporary saline streams 6 reaches on each"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (single reach with 12 sampling locations)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (1 intermittent and 1 perennial reach, 50m each)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (2 main study pools on 2 streams)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Network (3 streams in network, 8 sites each)"] <- "Network"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Network (2 catchments)"] <- "Network"
IR_scale_sub$Scale[IR_scale_sub$Scale == "50 m long reach"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Stream (3 reaches)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Stream (16 reaches)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Stream (2 streams)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Stream (2 channels of of a stream)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Stream (three 104-295 m long reaches)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach (226 km long, 4 sites)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Two temporary saline streams\\n6 reaches on each"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Reach\\n(10 temporary\\ntributaries)"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "Network\\n(12 intermittent streams\\n30-m study reaches)"] <- "Network"
IR_scale_sub$Scale[IR_scale_sub$Scale == "2 networks\\nFluvia` and Muga rivers (19 sites)"] <- "Network"
IR_scale_sub$Scale[IR_scale_sub$Scale == "<NA>"] <- "Reach"
IR_scale_sub$Scale[IR_scale_sub$Scale == "NA"] <- "Reach"


#Remove NA's
IRES_Country <-IRES_Country %>% drop_na()

IR_scale_sub$Scale <- as.factor(IR_scale_sub$Scale)


#subset the scale column for damming

DAM_scale_sub <- Dam_GHG[c(1,12,14)]

# Exclude review, metabolism and model studies (only include field or mesocosm studies)
names(DAM_scale_sub)[names(DAM_scale_sub) == "Study.type..Field..lab..both..review..conceptual."] <- "Study_type"

DAM_scale_sub <- DAM_scale_sub[!grepl("Review", DAM_scale_sub$Study_type),]
DAM_scale_sub <- DAM_scale_sub[!grepl("Model", DAM_scale_sub$Study_type),]
DAM_scale_sub <- DAM_scale_sub[!grepl("metabolism", DAM_scale_sub$Study_type),]

#Clean

DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Upstream and downstream of dam"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Downstream of dam"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Various"] <- "NA"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Upstream and downstream of dam, and forests"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Upstream and downstream of dam, pre and post restoration conditions"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Upstream and downstream of dams, and control/reference river"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Upstream and downstream of dams (8 dams along a river--7 on maainstem and 1 on tributary)"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Catchment"] <- "Network"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Downstream of dam (increasing distance)"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Upstream and downstream of dam (23 sites)"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Reach (regulated and unregulated river)"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Multiple catchments impacted by reservoirs"] <- "Network"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Upstream and downstream of dams"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Dam impacted network (56 sites)"] <- "Network"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "In reservoir and downstream of dam"] <- "Reach"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "Dam impacted network"] <- "Network"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "NA"] <- "<NA>"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == ""] <- "<NA>"
DAM_scale_sub$Scale[DAM_scale_sub$Scale == "<NA>"] <- "NA"


DAM_scale_sub$Scale <- as.factor(DAM_scale_sub$Scale)

DAM_scale_sub <-DAM_scale_sub %>% drop_na()

### Now combine for a figure

#add column for study topics

IR_scale_sub$Topic <- c("IRES")
DAM_scale_sub$Topic <- c("Dam")

#Now combine the IRES and Dam data for scale

all_scales <- rbind(IR_scale_sub, DAM_scale_sub)

#add a column for counts

all_scales_counts <- all_scales %>%
  filter(Topic %in% c("Dam", "IRES")) %>%
  group_by(Scale, Topic) %>%
  summarise(counts = n()) 

head(all_scales_counts, 4)
all_scales_counts

#Remove any NA's

all_scales_counts <- all_scales_counts %>% drop_na()

all_scales_counts

as.data.frame(all_scales_counts) #get rid of the blank and NA row


#reorder factor levels

as.data.frame(all_scales_counts)
all_scales_counts$Scale <- factor(all_scales_counts$Scale, levels=c("Reach", "Network", "Global"))

all_scales_counts %>%
  filter(!is.na(Scale))

#remove NA

all_scales_counts <- na.omit(all_scales_counts) 
  

tiff("scale_per_study_type_notext.tiff", units="in", width=4, height=4, res=1200)

scales_plot <- ggplot(all_scales_counts, aes(x = Topic, y = Scale, colour=Topic)) +
  geom_point(aes(size = counts), shape = 16, alpha=0.7) +
  scale_size_area(max_size = 40, guide = FALSE) +
 # geom_text(aes(
    #y = as.numeric(as.factor(Scale)) - sqrt(counts)/35, label = counts),
   # vjust = -0.5,
    #colour = "black",
    #size = 4) + 
  theme_classic() +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.line=element_blank(),axis.ticks=element_blank(), panel.border = element_blank(), panel.grid = element_blank()) +
  scale_color_manual(values=c("#79969B", "#C27D38")) 
scales_plot

# insert ggplot code
dev.off()

ggsave("scale_per_study_type.png", plot = scales_plot , units = c("in"), width = 4, height = 4 ) #5.42 x 4.41



################################

## Now make bubble plot of the main drivers per gas, per IRES/Dam to transpose into Drivers Figure

IRES_drivers <- read.csv("IRES_drivers_counts.csv")

#select just the count >1

IRES_drivers <- IRES_drivers %>% filter(Count > 1)


IRES_drivers_plot <- ggplot(IRES_drivers, aes(x = Gas, y = Factor, colour=Gas)) +
  geom_point(aes(size = Count), shape = 16, alpha=0.7) +
  scale_size_area(max_size = 40, guide = FALSE) +
  geom_text(aes(
    x = as.numeric(as.factor(Gas)) - sqrt(Count)/35, label = Count),
    colour = "black",
    size = 4
  ) + theme_classic() +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.line=element_blank(),axis.ticks=element_blank(), panel.border = element_blank(), panel.grid = element_blank()) +
  scale_color_manual(values=c("#C27D38", "#C27D38", "##C27D38")) 
IRES_drivers_plot

ggsave("IRES_drivers_plot.png", plot = IRES_drivers_plot , units = c("in"), width = 4, height = 5 ) #5.42 x 4.41



Dams_drivers <- read.csv("Dams_drivers_counts.csv")

Dams_drivers <- Dams_drivers %>% filter(Count > 2) #too many if we filter by 1


Dams_drivers_plot <- ggplot(Dams_drivers, aes(x = Gas, y = Factor, colour=Gas)) +
  geom_point(aes(size = Count), shape = 16, alpha=0.7) +
  scale_size_area(max_size = 40, guide = FALSE) +
  geom_text(aes(
    x = as.numeric(as.factor(Gas)) - sqrt(Count)/35, label = Count),
    colour = "black",
    size = 4
  ) + theme_classic() +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.line=element_blank(),axis.ticks=element_blank(), panel.border = element_blank(), panel.grid = element_blank()) +
  scale_color_manual(values=c("#79969B", "#79969B", "#79969B")) 
Dams_drivers_plot

ggsave("Dams_drivers_plot.png", plot = Dams_drivers_plot , units = c("in"), width = 5, height = 7.5 ) #5.42 x 4.41



#Now try to combine them 


all_drivers <- rbind(IRES_drivers, Dams_drivers)

#Make a new column pasting topic and gas

all_drivers$ID <-  paste(all_drivers$Topic, all_drivers$Gas, sep="_")



all_drivers_plot <- ggplot(all_drivers, aes(x = ID, y = Factor,colour=ID)) +
  geom_point(aes(size = Count), shape = 16, alpha=0.7) +
  scale_size_area(max_size = 40, guide = FALSE) +
  geom_text(aes(
    x = as.numeric(as.factor(ID)) - sqrt(Count)/35, label = Count),
    colour = "black",
    size = 4
  ) + theme_classic() +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.line=element_blank(),axis.ticks=element_blank(), panel.border = element_blank(), panel.grid = element_blank())  + scale_color_manual(values=c("#79969B", "#79969B", "#79969B", "#C27D38", "#C27D38")) 
all_drivers_plot

ggsave("all_drivers_plot.tiff", plot = all_drivers_plot , units = c("in"), width = 6.5, height = 11.5 ) #5.42 x 4.41
