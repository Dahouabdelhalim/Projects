#Figures and Stats for Pycno Manuscript
#Dr. Sarah Gravem(2021)
# Contents
   # Logistic Declines in Presence Scatterplot
  # Epidemic Timeline Bar Graph

#############################################################################################
#Creating figure 1a
#############################################################################################
#load ggplot
library(ggplot2)
library(RColorBrewer)
library(scales) # for formatting date axes later
setwd("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/FinalDataAndScripts/Data_ForGit/")

#Import Incidence data with already calculated predicted incidence from logistic regression model
library(readxl)
Inc <- read_excel("Incidence_2012-2019.xlsx")

#label vectors that will be in the graph
Date <-as.Date(Inc$date)
PredInc<-Inc$`Prob[1] By Region_Crash`
Reg <-as.factor(Inc$Region_Crash)

#change order of levelsfor graphing
Reg<-factor(Reg, levels=c("WesternAlaska"  ,"EastGulfofAlaska","SoutheastAlaska" , 
                          "BritishColumbia", "SalishSea", "Washington" , "Oregon"  ,
                          "NorthernCalifornia" , "CentralCalifornia" , "SouthernCalifornia", "Baja" ))


#full names for legend
legendlab <- c("western Alaska" ,"east Gulf of Alaska","southeast Alaska" , 
                 "British Columbia*", "Salish Sea", 
                 "Washington outer\\ncoast*" , "Oregon",
                 "northern California" , "central California" , "southern California", "Baja California" )

#Color wheel at https://chichacha.netlify.app/2018/12/09/having-bit-of-party-with-material-colour-palette/

#12 Regions
#Order is 
  #gray, pink, purple, 
  #indigo, blue, cyan, 
  #green, lime, yellow, 
  #orange, red, brown. 
cols12 <-   c("#9E9E9E","#F06292", "#9C27B0", 
              "#3F51B5", "#2196F3", "#00BCD4", 
              "#4CAF50", "#CDDC39", "#FDD835", 
              "#F99800", "#F44336", "#795548")

#11 Regions, When Aleuts and W GOA are merged (both pink): 
cols11 <-   c("#F06292", "#9C27B0", 
              "#3F51B5", "#2196F3", "#00BCD4", 
              "#4CAF50", "#CDDC39", "#FDD835", 
              "#F99800", "#F44336", "#795548")

#load data frame with crash date information for overlay
Events <- read_excel("CrashEventsForRPlot.xlsx", col_types = c("text", "text", "date", "numeric"))

#label vector for overlay points
DateEvent <- as.Date(Events$Date)
IncEvent <-Events$Incidence
RegEvent <-Events$Region_Crash
RegEvent<-factor(RegEvent, levels=c("WesternAlaska"  ,"EastGulfofAlaska","SoutheastAlaska" , 
                          "BritishColumbia", "SalishSea", "Washington outer coast*" , "Oregon"  ,
                          "NorthernCalifornia" , "CentralCalifornia" , "SouthernCalifornia", "Baja" ))

Scatter_PredInc_Date_Reg <- ggplot(Inc, aes(x=Date, y=PredInc, color = Reg))+
  geom_smooth(aes(y=PredInc), size = 1.5, se = FALSE)+ 
  geom_point(data = Events, aes( x = DateEvent, y = IncEvent, color = RegEvent), size = 3, shape = 21, fill = cols11, color = "black")+
  xlab("Date")+
  ylab(expression(paste("Predicted incidence of" , italic(" Pycnopodia")))) +
  labs (title = "a")+
  theme_classic()+
  theme(aspect.ratio = 3/3)+
  theme(axis.title = element_text(size=16)) + 
  theme(axis.text = element_text(size=16, colour="black")) + 
  theme(legend.text = element_text(size=11), legend.position = "left") + 
  theme(plot.title = element_text(size=12)) +
  scale_colour_manual(values= cols11, name="Region", labels=legendlab)+
  scale_y_continuous(limits=c(0, 1))+ 
  scale_x_date(date_breaks = "years" , date_labels = "%Y")+
  theme(axis.text.x = element_text(angle = 45 , hjust = 1))

plot(Scatter_PredInc_Date_Reg)
ggsave("Figure1a.jpeg",dpi = 300)
ggsave("Scatter_PredInc_Date_Reg.png")
ggsave("Scatter_PredInc_Date_Reg.eps")

############################################################################################################
#Outbreak Timeline Figure
#This is called a gantt chart
#resource https://www.molecularecologist.com/2019/01/03/simple-gantt-charts-in-r-with-ggplot2-and-the-tidyverse/
#needed to rearrange Regions_Timelines to have column for each region, phase, Start or End, and Date
############################################################################################################
library(tidyverse)
library(readxl)
EpidemicPhases <- read_excel("EpidemicPhases.xlsx")

#name key columns
reg <- as.factor(EpidemicPhases$Region_Full)
Phase <- as.factor(EpidemicPhases$EpidemicPhase)
StartEnd <- as.factor(EpidemicPhases$`Start-End`)
date <- as.Date(EpidemicPhases$Date)

#order columns and label regions with 
reg <- factor(reg, levels=rev(c("western Alaska^" ,"east Gulf of Alaska","southeast Alaska" , 
               "British Columbia*", "Salish Sea", 
               "Washington outer coast*" , "Oregon",
               "northern California" , "central California" , "southern California", "Baja California^" )))

Phase <- factor(Phase, levels = c("Pre-epidemic", "Emerging epidemic", "Epidemic", "Post-epidemic"))

#color the phases
colsPhase <- (c("#f9d926","#eb7556", "#9524a1", "#3d1d98"))

#make list of Ns for figure
annotations = rev(c("N = 54", "N = 236", "N = 219", "N = 1,552", "N = 10,128","N = 179","N = 503","N = 242","N = 1,114","^N = 2,758","N = 270" ))

#make a date value to place toxt at along the x axis
dateposition <- as.Date(c("2019-12-01"))

#graph
StackBar_EpiPhase_Reg <- ggplot(EpidemicPhases, aes(date, reg, color = Phase)) +
  geom_line(size = 6) +
  labs(x="Date", y="Region", title="b")+
  scale_colour_manual(values= colsPhase, name="Epidemic Phase")+
  theme_classic()+
  theme(aspect.ratio = 3/3)+
  theme(axis.title = element_text(size=16)) + 
  theme(axis.text.y = element_text(size=14, colour="black"), axis.text.x = element_text(size = 16, colour = "black")) + 
  theme(legend.position = "none") + 
  theme(plot.title = element_text(size=12)) +
  scale_x_date(date_breaks = "years" , date_labels = "%Y")+
  theme(axis.text.x = element_text(angle = 45 , hjust = 1))

plot(StackBar_EpiPhase_Reg)
ggsave("Fig1b.jpeg", dpi = 300)
ggsave("StackBar_EpiPhase_Reg_Date.png")
ggsave("StackBar_EpiPhase_Reg_Date.eps")
