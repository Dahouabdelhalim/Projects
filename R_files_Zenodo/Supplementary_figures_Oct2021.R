#African wild dog feeding priority
#Supplementary figures

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggsci)
library(ggplot2)
library(RColorBrewer)
library(ggforce)
library(wesanderson)
library(colorRamps) 
library(viridis)
library(cowplot)
library(ggsci)
library(LaplacesDemon)

theme_set(theme_minimal())

col <- "black"
colBG <- "white"
font <- "Helvetica"
fontSize <- 20

pointsize <- 3 # define size of points in the plot
lineSize <- 3 # define line size
alphaRibbon <- 0.3



#setwd to source file and read "Latency.csv")
#assuming R studio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
Lat_Feeding <- read.csv("Latency.csv", header=TRUE, sep=",")
names(Lat_Feeding)
Lat<-Lat_Feeding
             
colnames(Lat)[3] <- 'Latency'  
colnames(Lat)[4] <- 'InvLat'

Lat$Latency<-as.numeric(Lat$Latency)
Lat$InvLat<-as.numeric(Lat$InvLat)
Lat$PFQD2<-as.numeric(Lat$PFQD2)

summary(Lat)

cor(Lat$Latency, Lat$PFQD2) 
#-0.5698435
cor(Lat$Latency, Lat$PFQD2, method = "kendall") 
#-0.3517969

##labs(x = "Position in feeding queue (PFQ)", y="Latency to feed")+
g2<-ggplot(Lat, aes(x = PFQD2, y = InvLat)) +
  geom_point() +
  stat_smooth()+
  
  scale_y_continuous(name = "Latency to feed", limits = c(0,1), breaks = seq(0,1, by =0.2))+
  scale_x_continuous(name = "Position in Feeding Queue (PFQ)", limits = c(1,7), breaks = seq(0,7, by =1))+

# We use colorblind-friendly colors from the 'ggsci' package for fill colors:
scale_fill_jco()+
  
  # customize colors, font size, font, positions, etc. of axes, text, legends, ...
  theme(axis.text = element_text(colour = col, size = fontSize, family = font),
        axis.title.x = element_text(colour = col, size = fontSize, family = font, margin = margin(15,0,0,0)),
        axis.title.y = element_text(colour = col, size = fontSize, family = font, margin = margin(0,15,0,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = col),
        axis.ticks = element_line(colour = col),
        legend.position = c(0.6, 1), 
        legend.justification = c(1,1), 
        legend.box.margin=margin(c(0,20,0,20)),
        legend.text = element_text(colour = col, size = fontSize, family = font),
        legend.title = element_blank(),
        legend.key.width = unit(35,"pt"),
        plot.title = element_text(colour = col, size = fontSize, family = font, hjust = -0.068, margin = margin(0,0,10,0)),
        panel.border = element_blank(),
        plot.background = element_blank(),
  ) 

# view plot:
print(g2)

# export figure as PDF:
pdf(file = "SuppFig1.pdf", width = 15, height = 10, family = "Helvetica", fonts = "Helvetica")
g2
dev.off()



####
##Supp Figure 2
#Proportion of time spent feeding
Prop_Feeding <- read.csv("PropFeeding.csv", header=TRUE, sep=",")
names(Prop_Feeding)
PF<-Prop_Feeding

PF$PropFeed<-as.numeric(PF$PropFeed)
PF$PFQ_D2<-as.numeric(PF$PFQ_D2)

summary(PF)

g3<-ggplot(PF, aes(x = PFQ_D2, y = PropFeed)) +
  geom_point() +
  stat_smooth()+
  
  scale_y_continuous(name = "Proportion of time feeding", limits = c(0,1), breaks = seq(0,1, by =0.2))+
  scale_x_continuous(name = "Position in Feeding Queue (PFQ)", limits = c(1,7), breaks = seq(0,7, by =1))+
  
  # We use colorblind-friendly colors from the 'ggsci' package for fill colors:
  scale_fill_jco()+
  
  # customize colors, font size, font, positions, etc. of axes, text, legends, ...
  theme(axis.text = element_text(colour = col, size = fontSize, family = font),
        axis.title.x = element_text(colour = col, size = fontSize, family = font, margin = margin(15,0,0,0)),
        axis.title.y = element_text(colour = col, size = fontSize, family = font, margin = margin(0,15,0,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = col),
        axis.ticks = element_line(colour = col),
        legend.position = c(0.6, 1), 
        legend.justification = c(1,1), 
        legend.box.margin=margin(c(0,20,0,20)),
        legend.text = element_text(colour = col, size = fontSize, family = font),
        legend.title = element_blank(),
        legend.key.width = unit(35,"pt"),
        plot.title = element_text(colour = col, size = fontSize, family = font, hjust = -0.068, margin = margin(0,0,10,0)),
        panel.border = element_blank(),
        plot.background = element_blank(),
  ) 

# view plot:
print(g3)

# export figure as PDF:
pdf(file = "SuppFig2.pdf", width = 15, height = 10, family = "Helvetica", fonts = "Helvetica")
g3
dev.off()


#Supplemantary figure 3

#per capita approaches by pfq (and carcass stage)
#load packages
library(readxl)
library(ggplot2)

#Set theme
theme_set(theme_minimal())

col <- "black"
colBG <- "white"
font <- "Helvetica"
fontSize <- 20

pointsize <- 3 # define size of points in the plot
lineSize <- 3 # define line size
alphaRibbon <- 0.3


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read_excel("percapitaApproaches.xlsx")
data<-as.data.frame(data)
data$Carcass_stage<-as.factor(data$Carcass_stage)
data$PFQ<-as.factor(data$PFQ)
summary(data)
names(data)

data$Carcass_stage <- factor(data$Carcass_stage, levels = c("Early", "Mid", "Late"))

SupFig3<-ggplot(data, aes(x=PFQ, y=Mean_per_capita_approaches)) +
  geom_pointrange(
    aes(ymin = Mean_per_capita_approaches-SE, ymax = Mean_per_capita_approaches+SE, color = Carcass_stage),
    position = position_dodge(0.4), shape=19
  )+
  scale_color_manual(values = c("#2E9FDF", "#FFA500", "#B22222"))+
  theme(axis.text = element_text(colour = col, size = fontSize, family = font),
        axis.title.x = element_text(colour = col, size = fontSize, family = font, margin = margin(15,0,0,0)),
        axis.title.y = element_text(colour = col, size = fontSize, family = font, margin = margin(0,15,0,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = col),
        axis.ticks = element_line(colour = col),
        legend.position = c(1,1), 
        legend.justification = c(1,1), 
        legend.box.margin=margin(c(0,20,0,20)),
        legend.text = element_text(colour = col, size = 12, family = font),
        legend.title = element_text(colour = col, size = 12, family = font),
        legend.key.width = unit(35,"pt"),
        plot.title = element_text(colour = col, size = fontSize, family = font, hjust = -0.068, margin = margin(0,0,10,0)),
        panel.border = element_blank(),
        plot.background = element_blank(),
  ) +
  xlab("Position in feeding queue (PFQ)")+
  ylab("Mean approaches (per capita)")

# view plot:
print(SupFig3)

# export figure as PDF:
pdf(file = "SupFig3.pdf", width = 6, height = 5, family = "Helvetica", fonts = "Helvetica")
SupFig3
dev.off()



