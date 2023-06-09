#Ryan Stephens; Finalized on Aug. 19, 2020
rm(list=ls())# clears workspace
setwd("C:/Users/ryans/Documents/Isotopic_Routing_Small_Mammals/Map_Isotopes")# set working directory


#################################################################################################
#Global map of C3, C4, and mixed systems
#################################################################################################
library(raster)
library(dplyr)
library(ggplot2)

#########################
#System data 
#########################
#data on the C4 percentage of vegetation within a 1 degree grid cell (lat./long.) downloaded from
#https://daac.ornl.gov/ISLSCP_II/guides/c4_percent_1deg.html
map <- raster("C:/Users/ryans/Documents/Isotopic_Routing_Small_Mammals/Map_Isotopes/c4_percent_1d.asc")#C4 percentage data
map.p <- rasterToPoints(map)#convert the raster to points for plotting
df <- data.frame(map.p)#Make the points a dataframe for ggplot
DF<-select(df,long = x,lat = y, C4_percent  =c4_percent_1d)#rename headers

#C4 systems are 99% C4 plants or greater
C4<-filter(DF, C4_percent >= 99)
C4$C4<-"C4"
head(C4)

#Mixed systems are between 1% and 99% C4 plants
Mixed<-filter(DF, C4_percent > 1 & C4_percent < 99)
Mixed$Mixed<-"Mixed"
head(Mixed)
#########################


#########################
#Geographic feature data
#Data on coastlines, land masses, and glaciers downloaded from https://www.naturalearthdata.com/features/
#Coastline
coastline <- shapefile("C:/Users/ryans/Documents/Isotopic_Routing_Small_Mammals/Map_Isotopes/Coastlines/ne_110m_coastline.shp")
coastline_df <- fortify(coastline)# shapefile converted to a dataframe in ggplot2

#Land
land <- shapefile("C:/Users/ryans/Documents/Isotopic_Routing_Small_Mammals/Map_Isotopes/Land/ne_110m_land.shp")
land_df <- fortify(land)# shapefile converted to a dataframe in ggplot2

#Glacier
glacier <- shapefile("C:/Users/ryans/Documents/Isotopic_Routing_Small_Mammals/Map_Isotopes/Glacier/ne_110m_glaciated_areas.shp")
glacier_df <- fortify(glacier)# shapefile converted to a dataframe in ggplot2

#Field study locations (not used)
#Note that field study locations are collected from manuscripts or estimated (African sites are offset so that layered points are visible)
#Field_Study<- read.csv("Field_Study_locations.csv",header=T)
#head(Field_Study)
#########################


#########################
#Mapping
library(ggplot2)
Map<-ggplot(Mixed, aes(x=long, y=lat)) + 
#land (C3 plants) Plotted as background with mixed and C4 systems overlaid
  geom_polygon(data = land_df, aes(x = long, y = lat, group = group), size = .2, color=NA,fill = "#ABCB8C")+
#Mixed (C3 & C4 plants)
  geom_raster(aes(fill=Mixed), fill="mediumorchid3")+
#C4 plants
  geom_raster(data = C4, aes(x = long, y = lat,fill=C4),fill="#E9BA8C")+
#coastline (Mixed)
  geom_path(data = coastline_df, aes(x = long, y = lat, group = group),size = .2, color="mediumorchid3")+
#Greenland glaciers
  geom_polygon(data = glacier_df, aes(x = long, y = lat, group = group), size = .2, color=NA, fill = "white")+
#Field Locations
  #geom_point(data = Field_Study, aes(shape=Class,fill=Diet_Source),color="black",size=1.5, stroke=.25)+
  #scale_fill_manual(values =c("#ABCB8C","#E9BA8C","#96C0E8","mediumorchid3"))+
  #scale_shape_manual(values =c(24,21,22))+#different shapes for diet class
#set x and y coordinates
  scale_x_continuous(expand = c(0, 0),limits = c(-180, 180.1),breaks=c(-150,-100,-50,0,50,100,150))+
  scale_y_continuous(expand = c(0, 0),limits = c(-60, 85),breaks=c(-50,0,50))+
#background set to blue to be the ocean
  theme_bw()+
  theme(panel.background = element_rect(fill = "#96C0E8"),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+#removes grid lines, but keeps border
  theme(panel.grid.major.x = element_blank())+
#legend
  guides(fill=FALSE)+
  theme(legend.title=element_blank(),
        legend.position = c(0.08, 0.13),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.text=element_text(size=7, color ="gray50"),
        legend.spacing.x = unit(.05, 'cm'),
        legend.key.size = unit(.25, "cm"))+
#axes
  xlab("Longitude")+
  ylab("Latitude")+
  theme(axis.title.x= element_text(size=9),
        axis.text.x  = element_text(angle=0, hjust = 0.5,size=7, colour="black"))+
  theme(axis.title.y = element_text(size=9),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=7, colour="black"))
Map

ggsave("Map_Systems.jpeg", #could be tiff, png, etc (see below for saving as a pdf)
       plot = last_plot(), width =6, height = 3, units = "in",
       dpi = 600) 
#########################
#################################################################################################




#################################################################################################
#d13C and d15N isotopes (C3 plants, C4 plants, Plankton)
#################################################################################################

#########################
#Isotope data for C3 and C4 plants
#########################
#C3 and C4 plants downloaded from Plant Trait (TRY) Database https://www.try-db.org/TryWeb/Home.php
library(data.table)
Try_data<-fread("TRY_Data.txt",header=T, sep="\\t", dec =".", quote ="", data.table = T)#TRY data accessed March 23, 2020
Try_data<-as.data.frame(Try_data)
head(Try_data)

library(dplyr)
data<-select(Try_data,AccSpeciesName,TraitID,OrigValueStr,StdValue)#select needed columns
head(data)
str(data)
#Traits (these traits were all requested, but only 'pathway' and leaf d15N and leaf d13C were used for data analysis)
#22 = pathway (C3 or C4)
#146 = C/N rato
#78 = leaf d15N
#89 = leaf d13C
#14 = leaf N
#13 = leaf C

#get pathway for each plant species
pathway<-filter(data,TraitID=="22")#pathway
pathway<-select(pathway, AccSpeciesName, Pathway=OrigValueStr)#select needed columns and rename
head(pathway)

pathway<-filter(pathway,Pathway=="C3" |Pathway=="C4" )#filter to only C3 and C4 plants
pathway_disctinct<-distinct(pathway, AccSpeciesName, .keep_all = TRUE)#remove duplicate species names   
#########################


#########################
#d13C for plants
#########################
C13<-filter(data,TraitID=="89")#d13C
head(C13)
d13C<-select(C13,AccSpeciesName, d13C = StdValue)#rename values
head(d13C)

pathway_C13<-left_join(d13C, pathway_disctinct, by="AccSpeciesName")#match d13C values with pathway using species
head(pathway_C13)
str(pathway_C13)
pathway_C13<-na.omit(pathway_C13)#remove plants with unknown pathway
str(pathway_C13)

#preliminary plot of data
library(ggplot2)
ggplot(pathway_C13, aes(x = Pathway, y = d13C)) +#plot data
  geom_boxplot()

pathway_C13$Pathway<-as.factor(pathway_C13$Pathway)

#Based on manually looking at outlier points, a few plants were mislabeled
#No C3 plants had d13C values over -18 and no C4 plants had values under -20
pathway_C13<- pathway_C13 %>%  mutate(Pathway=ifelse(Pathway== "C3" & d13C >= -18 ,"C4",paste(Pathway)))#plants under -18 are C4
pathway_C13<- pathway_C13 %>%  mutate(Pathway=ifelse(Pathway== "C4" & d13C <= -20 ,"C3",paste(Pathway)))#plants over -20 are C3

#replot data
library(ggplot2)
ggplot(pathway_C13, aes(x = Pathway, y = d13C)) +
  geom_boxplot()

C3_C4<-select(pathway_C13,System=Pathway,d13C)#select needed columns
head(C3_C4)
#########################


#########################
#d13C for Marine plankton 
#########################
#Plankton isotope data kindly provided on request from K. W. McMahon:
#McMahon, K. W., L. L. Hamady, and S. R. Thorrold. 2013. Ocean ecogeochemistry: a review. 
#Oceanography and Marine Biology: An Annual Review 51:327-374.
Marine<- read.csv("Plankton_C.csv",header=T)#plankton data
Marine$System<-"Marine"#add system name
Marine<-select(Marine,System,d13C)#select needed columns
head(Marine)
#########################


#########################
#Graphing d13C
#########################
Food<-rbind(C3_C4,Marine)#bind plant d13C data to plankton d13C data
Food<-Food%>%mutate(System=ifelse(System == "C3","C[3]", ifelse(System == "C4","C[4]",paste(System))))#rename C3 and C4 for graphing
str(Food)

#get minimum and maximum for setting x axis limits
min(Food$d13C)
max(Food$d13C)

#get means and sd for points and plotting labels
Food_Summary<- Food%>% 
  group_by (System)%>%#Grouping variables
  summarise(Max = max(d13C), Min = min(d13C), Meadian = median(d13C), Mean=mean(d13C), SD=sd(d13C), n()) %>%
  as.data.frame()
Food_Summary$density<- -0.02# height to graph at
Food_Summary

#plot d13C data
library(ggplot2)
Hist_d13C<-ggplot(Food,aes(x=d13C, color=System, fill=System))+
#histogram
  geom_histogram(aes(y=..density..), size=.25, alpha=0.45, position='identity',binwidth=.5)+
  scale_color_manual(values =c("chartreuse4","darkorange3","dodgerblue3"))+#color by system
  scale_fill_manual(values =c("chartreuse4","darkorange3","dodgerblue3"))+#fill by system
#background
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white",colour = "white"))+
#legend
  theme(legend.position="none")+
#labels
  annotate("text", x = -29.5, y = .17, label = "C[3]", parse = TRUE, color="chartreuse4", size = 3)+#c3 label
  annotate("text", x = -23, y = .22, label = "Marine", color="dodgerblue3", size = 3)+#Marine label
  annotate("text", x = -13, y = .32, label = "C[4]", parse = TRUE, color="darkorange3", size = 3)+#C4 label
#means and sd
  geom_errorbarh(aes(x=Mean, y = density, xmin=Mean-SD, xmax=Mean+SD), data=Food_Summary, size=1, alpha=.75)+
  geom_point(data = Food_Summary, aes(x=Mean, y = density,fill=System),size=2)+
  geom_point(data = Food_Summary, aes(x=Mean, y = density),color="white",size=1)+
#axes
  scale_x_continuous(expand = c(0, 0),limits = c(-39,-9),breaks=c(-35,-30,-25,-20,-15,-10))+#set x axis
  scale_y_continuous(expand = c(0, 0),limits = c(-.04,.35),breaks=c(0,.1,.2,.3))+#set y axis
  xlab(expression(delta^13 * "C (\\211)"))+#x axis label
  ylab("Density")+#y axis label
  theme(axis.title.x= element_text(size=9),
        axis.text.x  = element_text(angle=0, hjust = 0.5,size=7, colour="black"))+
  theme(axis.title.y = element_text(size=9),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=7, colour="black"))

Hist_d13C

ggsave("System_d13C.jpeg", #could be tiff, png, etc (see below for saving as a pdf)
       plot = last_plot(), width =6, height = 2, units = "in",
       dpi = 600) 
#########################


#########################
#d15N for plants
#########################
N15<-filter(data,TraitID=="78")#filter out d15N data for plants
head(N15)
d15N<-select(N15,AccSpeciesName, d15N = StdValue)#select needed columns
head(d15N)

pathway_N15<-left_join(d15N, pathway_disctinct, by="AccSpeciesName")#match d13C values with pathway
head(pathway_N15)
pathway_N15<-na.omit(pathway_N15)#remove plants with unknown pathway

library(ggplot2)#preliminary plot (although no way to verify that records are correct)
ggplot(pathway_N15, aes(x = Pathway, y = d15N)) +
  geom_boxplot()

C3_C4_15N<-select(pathway_N15,System=Pathway,d15N)#select needed data
head(C3_C4_15N)
#########################


#########################
#d15N for Marine plankton 
#########################
#Plankton data kindly provided on request from K. W. McMahon:
#McMahon, K. W., L. L. Hamady, and S. R. Thorrold. 2013. Ocean ecogeochemistry: a review. 
#Oceanography and Marine Biology: An Annual Review 51:327-374.
Marine<- read.csv("Plankton_N.csv",header=T)#plankton data
Marine$System<-"Marine"#add system
Marine_15N<-select(Marine,System,d15N)#select need columns
head(Marine_15N)
#########################


#########################
#Graphing d15N
#########################
Food_15N<-rbind(C3_C4_15N,Marine_15N)#bind plant d15N data to marine d15N data
head(Food_15N)
Food_15N<-Food_15N%>%mutate(System=ifelse(System == "C3","C[3]",ifelse(System == "C4","C[4]",paste(System))))#rename for graphing

#get minimum and maximum for setting x axis limits
min(Food_15N$d15N)
max(Food_15N$d15N)

#get means and sd for points and plotting labels
Food_15N_Summary<- Food_15N%>% 
  group_by (System) %>%#Grouping variables
  summarise(Max = max(d15N), Min = min(d15N), Meadian = median(d15N), Mean=mean(d15N), SD=sd(d15N), n()) %>%
  as.data.frame()
Food_15N_Summary$density<- -0.02

#plot d15N data
library(ggplot2)
Hist_d15N<-ggplot(Food_15N,aes(x=d15N, color=System, fill=System))+
#Histogram
  geom_histogram(aes(y=..density..),size=.25, alpha=0.45, position='identity',binwidth=.5)+
  scale_color_manual(values =c("chartreuse4","darkorange3","dodgerblue3"))+#color by system
  scale_fill_manual(values =c("chartreuse4","darkorange3","dodgerblue3"))+#fill color by system
#background
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white",colour = "white"))+
#legend
  theme(legend.position="none")+
#labels
  annotate("text", x = 1, y = .17, label = "C[3]", parse = TRUE, color="chartreuse4", size = 3)+#c3
  annotate("text", x = 2, y = .2, label = "C[4]", parse = TRUE, color="darkorange3", size = 3)+ #C4
  annotate("text", x = 7.1, y = .18, label = "Marine", color="dodgerblue3", size = 3)+#Marine
#Means and sd
  geom_errorbarh(aes(x=Mean, y = density, xmin=Mean-SD, xmax=Mean+SD), data=Food_15N_Summary, height=0, size=1,alpha=.75)+
  geom_point(data = Food_15N_Summary, aes(x=Mean, y = density,fill=System),size=2)+
  geom_point(data = Food_15N_Summary, aes(x=Mean, y = density),color="white",size=1)+
#Axes
  scale_x_continuous(expand = c(0, 0),limits = c(-15,22),breaks=c(-10,-5,-0,5,10,15,20))+
  scale_y_continuous(expand = c(0, 0),limits = c(-.04,.35),breaks=c(0,.1,.2,.3))+
  xlab(expression(delta^15 * "N (\\211)"))+#x axis label
  ylab("Density")+#y axis label
  theme(axis.title.x= element_text(size=9),
        axis.text.x  = element_text(angle=0, hjust = 0.5,size=7, colour="black"))+
  theme(axis.title.y = element_text(size=9),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=7, colour="black"))

Hist_d15N

ggsave("System_d15N.jpeg", #could be tiff, png, etc (see below for saving as a pdf)
       plot = last_plot(), width =6, height = 2, units = "in",
       dpi = 600) 
#########################
#################################################################################################




#################################################################################################
#Compile graphs using Cowplot to make final figure
#################################################################################################
#Use silhouettes from http://phylopic.org/ to represent systems
library(rphylopic)
Plankton <- image_data("4924b6bd-cfb8-4d60-a32a-442d02afbe85", size = "512")[[1]]#plankton
#Plankton1 <- image_data("8eadd350-45ac-4bc2-8994-b163aee72609", size = "512")[[1]]#not used
#Salmon <- image_data("7413aa3a-d736-435a-8635-0c316ff73f26", size = "512")[[1]]#not used
#Squid <- image_data("443c7494-aa78-4a21-b2b4-cfa3639e1346", size = "512")[[1]]#not used
Tree <- image_data("1b329337-8f8b-4380-8aae-23d50e0db14f", size = "512")[[1]]#C3 plant
Corn <- image_data("cd2aa8d4-a965-44dc-875f-bbc658813f30", size = "512")[[1]]#C4 plant

#Complied figure
library(cowplot)
ggdraw() +
  #plots
  draw_plot(Map, x = 0, y = 0, width = .75, height = .6)+#Map
  draw_plot(Hist_d13C, x = 0, y = .6, width = .5, height = .4)+#d13C histogram
  draw_plot(Hist_d15N, x = .5, y = .6, width = .5, height = .4)+#d15N histogram
  #legend
  draw_label("System",           x = 0.87, y = .53, size=9)+#legend title
  draw_label(parse(text="C[3]"), x = .88, y = .44, size = 9, hjust = 0, color="chartreuse4")+#C3 label
  draw_label(parse(text="C[4]"), x = .88, y = .32, size = 9, hjust = 0, color="darkorange3")+#C4 label
  draw_label("Marine",           x = .88, y = .20, size = 9, hjust = 0, color="dodgerblue3")+#marine label
  draw_label("Mixed",            x = .88, y = .09, size = 9, hjust = 0, color="mediumorchid3")+#mixed label
  #C3 image
  add_phylopic(Tree, alpha=0.50, x=.85, y=.44, ysize=.1, color="chartreuse4")+#C3 plants
  #C4 image
  add_phylopic(Corn, alpha=0.50, x=.84, y=.32, ysize=.1, color="darkorange3")+#C4 plants
  add_phylopic(Corn, alpha=0.50, x=.85, y=.32, ysize=.1, color="darkorange3")+#C4 plants
  add_phylopic(Corn, alpha=0.50, x=.86, y=.32, ysize=.1, color="darkorange3")+#C4 plants
  #marine image
  add_phylopic(Plankton, alpha=1, x=.85, y=.205, ysize=.07, color="dodgerblue3")+#marine
  #add_phylopic(Salmon, alpha=0.50, x=.83, y=.655, ysize=.04, color="dodgerblue3")+#not used
  #draw_image(Marine,hjust = -.352, vjust = -.16,scale = .09)+#not used
  #Mixed image
  add_phylopic(Tree, alpha=0.50, x=.835, y=.11, ysize=.07, color="mediumorchid3")+#C3 plants
  add_phylopic(Corn, alpha=0.50, x=.85, y=.11, ysize=.07, color="mediumorchid3")+#C4 plants
  add_phylopic(Corn, alpha=0.50, x=.855, y=.11, ysize=.07, color="mediumorchid3")+#C4 plants
  add_phylopic(Corn, alpha=0.50, x=.86, y=.11, ysize=.07, color="mediumorchid3")+#C4 plants
  add_phylopic(Plankton, alpha=1, x=.85, y=.05, ysize=.050, color="mediumorchid3")+#marine
  #draw_image(Marine_mixed,hjust = -.362, vjust = -.04,scale = .09)+#not used
  #letters for plots
  draw_plot_label(label = c("(a)", "(b)", "(c)"), size = 9,
                  x = c(-.01, .49, -.01),
                  y = c(1, 1, .6))

ggsave("Figure_1.tiff",
       plot = last_plot(), width = 6, height = 3.8, units = "in", dpi = 600) 
#################################################################################################




#################################################################################################
#Get C/N ratio for plant species used in meta analysis
#################################################################################################
CN<-filter(data,TraitID=="146")#pathway
head(CN)
library(tidyverse)
Plant_C_N<-CN %>% filter(str_detect(AccSpeciesName, 'Cynodon'))
Plant_C_N$OrigValueStr<-as.numeric(Plant_C_N$OrigValueStr)
str(Plant_C_N)
Plant_C_N%>%
  summarise(Mean=mean(OrigValueStr), n = n())
#################################################################################################












