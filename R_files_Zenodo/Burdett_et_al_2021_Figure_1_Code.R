#Burdett et al. 2021
#Figure 1 creation code
# Dec 2019

library(tidyverse)
library(RColorBrewer)
library(wesanderson)
library(patchwork)
library(ggpubr)
library(ggpattern)
library(lemon)

#To create Figure 1: Stacked bar chart
#Load and tidy data
df<-read_csv("Burdett_et_al_2021_Figure_1_Data.csv")
df$Keyword <- factor(df$Keyword, levels=c("Public engagement", "Education", "Outreach", "Science communication"))
talk<-subset(df, Type=="Talks")
paper<-subset(df, Type=="Papers")
df$Year<-as.factor(df$Year)
df$Type<-factor(df$Type, levels=c("Papers", "Talks"))

## Generate plots with facet
#Color version
fig_col<-ggplot(df, aes(x = Count, y = reorder(Year, desc(Year)), fill = Keyword)) +
  geom_col() +theme_classic()+
  facet_rep_grid(. ~ Type) + 
  coord_capped_cart(bottom='both', left='both')+
  ylab("Year")+
  scale_fill_manual(values=c("#fdae61", "#2c7bb6", "#d7191c", "#abd9e9"))+
  labs(fill = "Key word")+
  theme(axis.title.x=element_text(size=16), 
        axis.title.y = element_text(size=16), 
        axis.text.x = element_text(size=12, vjust=0.5), 
        axis.text.y=element_text(size=12),
        strip.background = element_blank(),
        strip.text.x = element_text(size=16, face="bold", hjust=-0.02),
        legend.text=element_text(size=12),
        legend.justification=c(1,1), legend.position=c(0.98,0.98),
        legend.title = element_text(size=12, face="bold"))
fig_col

ggsave("C:/Users/owner/Desktop/SFS/Fig1_Color.tiff",
       plot = fig_col,
       width=7.5,
       height=5,
       units=c("in"),
       dpi = 600)

#Greyscale version
fig_gry<-ggplot(df, aes(x = Count, y = reorder(Year, desc(Year)), fill = Keyword)) +
  geom_col() +theme_classic()+
  facet_rep_grid(. ~ Type) + 
  coord_capped_cart(bottom='both', left='both')+
  ylab("Year")+
  scale_fill_grey()+
  labs(fill = "Key word")+
  theme(axis.title.x=element_text(size=16), 
        axis.title.y = element_text(size=16), 
        axis.text.x = element_text(size=12, vjust=0.5), 
        axis.text.y=element_text(size=12),
        strip.background = element_blank(),
        strip.text.x = element_text(size=16, face="bold", hjust=-0.02),
        legend.text=element_text(size=12),
        legend.justification=c(1,1), legend.position=c(0.98,0.98),
        legend.title = element_text(size=12, face="bold"))
fig_gry

ggsave("C:/Users/owner/Desktop/SFS/Fig1_Grey_12.12.tiff",
       plot = fig_gry,
       width=7.5,
       height=5,
       units=c("in"),
       dpi = 600)
