#bodymass.R

#LAST UPDATED 9/5/2021

#This script uses body mass data from the qmr.final.csv dataset.
#This script the smooth body function to run to identify seasonal cutoff dates
#for each individual bird.
#This script generates figures for each bird, shown in Fig. S7.
#This script generates the period.cutoffs.csv dataset which will be used in 
#the rq_processing.R script to label each run with its appropriate period

################################################################################
################################################################################

#set working directory
setwd("~/ERICH_EBERTS/ACADEMICS/PROJECTS/FAT/FAT_DATA/FAT_CODE/FINAL_DATA")

################################################################################

#load in libraries
library(readxl)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(plotrix)
library(gridExtra)
library(changepoint)
library(changepoint.np)
library(MESS)
library(nlme) 
library(effects)
library(rmcorr)

################################################################################

#set my theme for any ggplot graphs i output.

my_theme <- theme_classic(base_size = 14) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  theme(axis.text.x = element_text(angle=0, size=16), legend.key.height = unit(3, 'lines'),
        panel.grid.major.y = element_line(color='white'),
        panel.grid.minor.y = element_line(color='white'))+
  theme(axis.text.y = element_text(angle=0, size=16), legend.key.height = unit(3, 'lines'),
        panel.grid.major.x = element_line(color='white'))+
  theme(plot.caption = element_text(size = 10, hjust = 0.5))+
  theme(axis.title.x = element_text(vjust=-3, size=20))+
  theme(axis.title.y = element_text(vjust= 4,size=20))+
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))+
    theme(legend.position = "none") +
  #theme(legend.title=element_blank())+

  #for facet plot labels
  # theme(strip.text.x = element_text(size=8, color="black",
  #                                     face="plain", hjust=0))+
  #
  # theme(strip.background = element_rect(colour="black", fill="white",
  #                                      size=.5, linetype="solid"))

  theme(strip.background = element_blank(), strip.text.x = element_blank())
  
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

#LOAD IN QMR DATA

setwd("~/ERICH_EBERTS/ACADEMICS/PROJECTS/FAT/FAT_DATA/FAT_CODE/FINAL_DATA")
qmr.final<-read.csv("qmr_final.csv")

qmr.final$Date_Run<-as.Date(qmr.final$Date,"%m/%d/%Y" )

################################################################################
################################################################################

#LOAD IN DATE INDEX FILE

date.index<-read.csv("date_index.csv")
date.index$Date<- as.Date(date.index$Date,"%m/%d/%Y")

################################################################################

#ADD DATE INDEX TO qmr.final

qmr.final$Date_Index<-NA
fuck<-c(1:1092)
for(i in fuck) {
qmr.final$Date_Index[i]<-date.index$Date_Index[date.index$Date==qmr.final$Date_Run[i]]
  }

qmr.final$Date_Index<-as.numeric(as.character(qmr.final$Date_Index))
#str(qmr.final$Date_Index)

################################################################################
################################################################################

#SUBSET OUT ONLY EVENING BODY MASS
qmr.final.hold<-qmr.final

qmr.final<-qmr.final[qmr.final.hold$Morn_Eve=="MORN",]

################################################################################
#LOAD IN run_agg.csv to add in torpor use variable

run.agg<-read.csv("run_agg.csv")

qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_7.22.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_7.22.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_7.26.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_7.26.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_7.30.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_7.30.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_8.15.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_8.15.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_8.22.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_8.22.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_8.31.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_8.31.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_8.8.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_8.8.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_9.12.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_9.12.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_9.20.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_9.20.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_9.26.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_9.26.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B1_9.9.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B1_9.9.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_6.25.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_6.25.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_6.27.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_6.27.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_8.20.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_8.20.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_8.31.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_8.31.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_8.9.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_8.9.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_9.1.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_9.1.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_9.10.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_9.10.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_9.2.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_9.2.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_9.20.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_9.20.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_9.21.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_9.21.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_9.7.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_9.7.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B10_9.9.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B10_9.9.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B11_6.25.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B11_6.25.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B11_6.29.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B11_6.29.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B11_8.20.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B11_8.20.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B11_9.11.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B11_9.11.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B11_9.12.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B11_9.12.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B11_9.16.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B11_9.16.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B11_9.24.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B11_9.24.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B12_6.26.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B12_6.26.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B12_7.25.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B12_7.25.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B12_8.21.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B12_8.21.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B12_8.26.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B12_8.26.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B12_9.21.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B12_9.21.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B12_9.22.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B12_9.22.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B12_9.4.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B12_9.4.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B12_9.5.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B12_9.5.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B12_9.6.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B12_9.6.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_6.29.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_6.29.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_7.2.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_7.2.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_8.11.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_8.11.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_8.24.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_8.24.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_9.1.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_9.1.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_9.12.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_9.12.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_9.16.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_9.16.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_9.17.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_9.17.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_9.2.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_9.2.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_9.22.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_9.22.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B13_9.23.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B13_9.23.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_6.29.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_6.29.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_7.2.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_7.2.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_8.12.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_8.12.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_8.21.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_8.21.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_8.23.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_8.23.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_8.31.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_8.31.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_9.1.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_9.1.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_9.22.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_9.22.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_9.23.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_9.23.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_9.24.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_9.24.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_9.4.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_9.4.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_9.5.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_9.5.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B14_9.6.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B14_9.6.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B15_7.7.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B15_7.7.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B15_8.12.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B15_8.12.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B15_8.23.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B15_8.23.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B15_9.10.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B15_9.10.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B15_9.2.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B15_9.2.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B15_9.23.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B15_9.23.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B15_9.24.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B15_9.24.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B15_9.7.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B15_9.7.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B15_9.9.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B15_9.9.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B16_7.25.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B16_7.25.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B16_8.12.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B16_8.12.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B16_8.24.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B16_8.24.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B16_9.12.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B16_9.12.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B16_9.16.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B16_9.16.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B16_9.17.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B16_9.17.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B16_9.3.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B16_9.3.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_7.19.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_7.19.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_7.22.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_7.22.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_7.29.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_7.29.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_8.13.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_8.13.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_8.20.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_8.20.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_8.26.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_8.26.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_8.31.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_8.31.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_8.6.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_8.6.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_9.15.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_9.15.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_9.20.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_9.20.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_9.22.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_9.22.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_9.23.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_9.23.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B2_9.9.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B2_9.9.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_7.19.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_7.19.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_7.26.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_7.26.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_7.30.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_7.30.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_8.15.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_8.15.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_8.24.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_8.24.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_8.29.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_8.29.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_8.8.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_8.8.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_9.10.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_9.10.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_9.17.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_9.17.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_9.2.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_9.2.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_9.22.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_9.22.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_9.23.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_9.23.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B3_9.27.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B3_9.27.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_7.29.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_7.29.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_8.13.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_8.13.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_8.20.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_8.20.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_8.26.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_8.26.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_8.29.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_8.29.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_8.6.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_8.6.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_9.12.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_9.12.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_9.19.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_9.19.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_9.26.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_9.26.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_9.27.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_9.27.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B4_9.6.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B4_9.6.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B5_8.22.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B5_8.22.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B5_8.24.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B5_8.24.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B5_9.10.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B5_9.10.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B5_9.15.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B5_9.15.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B5_9.17.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B5_9.17.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B5_9.19.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B5_9.19.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B5_9.6.18"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B5_9.6.18"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B6_6.24.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B6_6.24.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B6_6.26.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B6_6.26.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B6_8.13.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B6_8.13.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B6_8.26.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B6_8.26.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B6_8.9.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B6_8.9.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B6_9.18.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B6_9.18.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B6_9.3.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B6_9.3.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B6_9.7.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B6_9.7.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B7_6.24.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B7_6.24.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B7_6.26.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B7_6.26.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B7_8.20.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B7_8.20.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B7_8.31.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B7_8.31.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B7_8.9.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B7_8.9.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B7_9.11.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B7_9.11.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B7_9.18.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B7_9.18.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B8_6.24.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B8_6.24.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B8_6.27.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B8_6.27.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B8_8.11.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B8_8.11.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B8_8.13.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B8_8.13.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B8_8.23.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B8_8.23.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B8_9.10.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B8_9.10.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B8_9.18.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B8_9.18.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B8_9.20.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B8_9.20.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B8_9.9.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B8_9.9.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_6.25.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_6.25.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_6.27.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_6.27.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_8.11.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_8.11.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_8.21.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_8.21.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_8.26.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_8.26.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_9.20.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_9.20.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_9.21.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_9.21.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_9.3.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_9.3.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_9.4.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_9.4.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_9.5.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_9.5.19"][1]
qmr.final$Torpor_Use[qmr.final$Run_ID=="B9_9.6.19"]<-run.agg$Torpor_Use[run.agg$Run_ID=="B9_9.6.19"][1]

################################################################################

#END OF DATA SETUP
################################################################################
################################################################################
################################################################################
################################################################################

#SETUP INTERP DATA TO INCLUDE DAY/NIGHT LENGTHS

interp.data<-data.frame(1:122)
names(interp.data)<-"Date_Index"

interp.data$Night_Duration<-NA

for(i in c(1:122)){
interp.data$Night_Duration[interp.data$Date_Index==i]<-mean(run.agg$Night_Duration[run.agg$Date_Index==i], na.rm=T)
}

    #GENERATE VECTORS OF INTERPOLATED**, SMOOTHED BODY MASS DATA AND ITS SLOPE

    #interpolate body mass data linearly #THIS IS NOT NEEDED IF IM USING RAW. 
    interp.data$Night_Duration<-na.approx(interp.data$Night_Duration,x=interp.data$Date_Index,
                          xout=c(1:122), na.rm=F, rule=1)
    
    interp.data$Week<-floor(interp.data$Date_Index/7)
    interp.data$Week
    
    for(i in c(1:17)){
    interp.data$Night_Duration_Week[interp.data$Week==i]<-floor(mean(interp.data$Night_Duration[interp.data$Week==i], na.rm=T))
    }

  interp.data$Night_Duration_Week[interp.data$Week<3]<-525  
  interp.data$Night_Duration_Week<-interp.data$Night_Duration_Week/60

  ggplot(interp.data)+
    geom_point(aes(x=Date_Index, y=Night_Duration_Week), col="red",fill="red", lwd=5)+
    #geom_col(aes(x=Date_Index, y=Night_Duration_Week), col="red",fill="red")+

      scale_y_continuous(limits = c(0,12),breaks = seq(0, 15, by = 1), name = "Night Duration (Hours)")+
  scale_x_continuous(limits = c(0,125),breaks = seq(0, 150, by = 10), name="Date Index, days since June 1")+
  
  my_theme

  #   interp.data.week<-interp.data %>%
  # group_by(Week) %>%
  # #summarize (usually mean +- std error) variables that have multiple values per run (i.e. minute data)
  # summarize(Night_Duration = mean(Night_Duration, na.rm = T))%>%
  #         data.frame()


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

#SMOOTH BODY MASS FUNCTION

BIRD_ID<-"B3"
SPAR=.35
CUTOFF=.06

smooth_body_func<-function(BIRD_ID, SPAR, CUTOFF){
    ################################################################################
    
    #First create a new dataframe to populate #populate it with date_index from 1:122 (june 1 to sept 30)
    interp.data<-data.frame(1:122)
    names(interp.data)<-"Date_Index"
    
    rundateindexs<-run.agg$Date_Index[run.agg$Bird_ID==BIRD_ID]+1

interp.data<-data.frame(1:122)
names(interp.data)<-"Date_Index"
  #label interp.data with torpor/normothermy on each particular date index
  #so first, wherever i have a run id, i want that date index to be labeled as such in the interp data. 
  interp.data$runyn[interp.data$Date_Index %in% rundateindexs]<-"run"
  interp.data$runyn[!interp.data$Date_Index %in% rundateindexs]<-"no"
  #blank the original torpo use variable to be overridden with the proper label. 
  interp.data$Torpor_Use<-NA
  
  #Then actually label torpor use and the the run id properly. 
  interp.data$Torpor_Use[interp.data$Date_Index %in% qmr.final$Date_Index[qmr.final$Bird_ID==BIRD_ID]]<-qmr.final$Torpor_Use[qmr.final$Bird_ID==BIRD_ID]
  interp.data$Run_ID[interp.data$Date_Index %in% qmr.final$Date_Index[qmr.final$Bird_ID==BIRD_ID]]<-qmr.final$Run_ID[qmr.final$Bird_ID==BIRD_ID]

    
  #add the bird id to the dataframe
  interp.data$Bird_ID<-unique(qmr.final$Bird_ID[qmr.final$Bird_ID==BIRD_ID&qmr.final$Morn_Eve=="MORN"]) 

    
#add night duration into interp data. 
  interp.data$Night_Duration<-NA
  
  for(i in c(1:122)){
  interp.data$Night_Duration[interp.data$Date_Index==i]<-mean(run.agg$Night_Duration[run.agg$Date_Index==i], na.rm=T)
  }
  
      #GENERATE VECTORS OF INTERPOLATED**, SMOOTHED BODY MASS DATA AND ITS SLOPE
  
      #interpolate
      interp.data$Night_Duration<-na.approx(interp.data$Night_Duration,x=interp.data$Date_Index,
                            xout=c(1:122), na.rm=F, rule=1)
      
      interp.data$Week<-floor(interp.data$Date_Index/7)
      interp.data$Week
      
      for(i in c(1:17)){
      interp.data$Night_Duration_Week[interp.data$Week==i]<-floor(mean(interp.data$Night_Duration[interp.data$Week==i], na.rm=T))
      }
  
    interp.data$Night_Duration_Week[interp.data$Week<3]<-525  
    interp.data$Night_Duration_Week<-interp.data$Night_Duration_Week/60
      

    ################################################################################
    
    #GENERATE VECTORS OF INTERPOLATED**, SMOOTHED BODY MASS DATA AND ITS SLOPE

    #interpolate body mass data linearly #THIS IS NOT NEEDED IF IM USING RAW. 
    interp.data$Body_g_int<-na.approx(qmr.final$Body_g[qmr.final$Bird_ID==BIRD_ID&qmr.final$Morn_Eve=="MORN"], 
                          x=qmr.final$Date_Index[qmr.final$Bird_ID==BIRD_ID&qmr.final$Morn_Eve=="MORN"],
                          xout=c(1:122), na.rm=F, rule=1)
    
    #generate a spline model to smooth RAW data
    smoothmodel<-smooth.spline(qmr.final$Body_g[qmr.final$Bird_ID==BIRD_ID & is.na(qmr.final$Body_g)==F],
                          x=qmr.final$Date_Index[qmr.final$Bird_ID==BIRD_ID& is.na(qmr.final$Body_g)==F],
                            spar=SPAR)
    
    #generate smoothed BODY MASS data for each minute, based on above spline model
    predmodel<- predict(smoothmodel, x=seq(1,122, length=122))
    interp.data$Body_g_smooth<-predmodel$y #save to interp.data
    
    #generate the first derivative (slope) of smoothed BODY MASS data for each minute, based on above spline model
    predmodel_deriv<- predict(smoothmodel, x=seq(1,122, length=122), deriv=1)
    interp.data$Body_g_smoothderiv<-predmodel_deriv$y#save to interp.data
     
    # #test plots
    # ##interpolated data
    # plot(interp.data$Body_g_int)
    # ##smoothed data
    # plot(interp.data$Body_g_smooth)
    # ##derivitive data
    # plot(interp.data$Body_g_smoothderiv)
    
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    ################################################################################
    
    #ASSIGN PERIOD LABELS FOR DIFFERNET TYPES OF BODY MASS CHANGE RATES
    
    #some notes: 
    #dont include the first day of fattening as fattening. 
    #dont include the last day of fattening as fattening.

    ############################################################################
    
    #Where change in body mass is less than the bird specific cutoff argument, label that "Summer"
    interp.data$Season[interp.data$Body_g_smoothderiv <= CUTOFF]<- "Summer"
      #make some bird specific adjustments
      if(BIRD_ID=="B4"){interp.data$Season[interp.data$Date_Index >= 86 &interp.data$Date_Index <= 105]<- "Fattening"}

    #Where change in body mass is greater than or equal to the bird specific cutoff argument, label that "Fattening"
    interp.data$Season[interp.data$Body_g_smoothderiv > CUTOFF]<- "Fattening"

    #create a bird specific DATE INDEX OF FATTENING START variable in interp.data 
    interp.data$Fat_Start<-min(interp.data$Date_Index[interp.data$Season== "Fattening"])
      #make some bird specific adjustments
      if(BIRD_ID=="B11"){interp.data$Fat_Start<-interp.data$Fat_Start-5
                          interp.data$Season[interp.data$Date_Index > interp.data$Fat_Start[1] &interp.data$Date_Index< max(interp.data$Date_Index[interp.data$Season== "Fattening"])]<- "Fattening"}
      if(BIRD_ID=="B14"){interp.data$Fat_Start<-interp.data$Fat_Start-5
                          interp.data$Season[interp.data$Date_Index > interp.data$Fat_Start[1] &interp.data$Date_Index< max(interp.data$Date_Index[interp.data$Season== "Fattening"])]<- "Fattening"}
      if(BIRD_ID=="B15"){interp.data$Fat_Start<-interp.data$Fat_Start-2
                          interp.data$Season[interp.data$Date_Index > interp.data$Fat_Start[1] &interp.data$Date_Index< max(interp.data$Date_Index[interp.data$Season== "Fattening"])]<- "Fattening"}
    
    #create a bird specific DATE INDEX OF FATTENING END variable in interp.data 
    interp.data$Fat_End<-max(interp.data$Date_Index[interp.data$Season== "Fattening"])
      
      if(BIRD_ID=="B9"){interp.data$Fat_End<-interp.data$Fat_End+1}
      if(BIRD_ID=="B11"){interp.data$Fat_End<-interp.data$Fat_End-1}
      if(BIRD_ID=="B3"){interp.data$Fat_End<-interp.data$Fat_End+1}
      if(BIRD_ID=="B4"){interp.data$Fat_End<-interp.data$Fat_End-1}
      
    #Where change in body mass is greater than or equal to the bird specific cutoff argument, label that "Migration"
    interp.data$Season[interp.data$Date_Index >= interp.data$Fat_End[1]]<- "Migration" 
    
    #create a bird specific DATE INDEX OF CAPTURE variable in interp.data 
    #IDENTIFY CATCH DATE. 
    #if the data points in the interpolated vector are original points from actual measured body masses, then mark those rows as true
    interp.data$Test<-interp.data$Body_g_int %in% qmr.final$Body_g[qmr.final$Bird_ID==BIRD_ID & is.na(qmr.final$Body_g)==F]
    interp.data$Catch_Date<-interp.data$Date_Index[interp.data$Test=="TRUE"][1]
    #save as a local variable to use in plotting
    Catch_Date<-interp.data$Date_Index[interp.data$Test=="TRUE"][1]

  ################################################################################
  #PLOT
  
  #assign the cut off date index to the object cutoff to reference easily, this will be used in the following plots
  cutoff<-as.character(interp.data$Fat_Start[1])
  
  plot<-ggplot(data=interp.data)+
  #background colors, by season, referencing the Fat_Start and Fat_End date indexes
    geom_rect(xmin=0, xmax=interp.data$Fat_Start, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=.02)+
    geom_rect(xmin=interp.data$Fat_Start, xmax=interp.data$Fat_End, ymin=-Inf, ymax=Inf, fill="orange", alpha=.002)+
    geom_rect(xmin=interp.data$Fat_End, xmax=124, ymin=-Inf, ymax=Inf, fill="darkred", alpha=.002)+
  
  #line for the smoothed data,
    geom_line(data=interp.data[interp.data$Date_Index>interp.data$Catch_Date,],aes(x= Date_Index, y= Body_g_smooth), alpha=.3, col="black",lwd=2)+
    
  #lines for seasonal linear regressions,
    geom_smooth(data=interp.data,aes(x=Date_Index, y= Body_g_int, group=Season),col="black",alpha=1,
                      se=F, method="lm",lwd=1)+
  
    #lowess line based on interpolated data
    #geom_smooth(aes(x=Date_Index, y= Body_g_int), method="loess", col="red",se=F)+
  
  #vertical line indicating CAPTURE  
    #geom_vline(data=interp.data,aes(xintercept = interp.data$Catch_Date[1]), col="blue", linetype="dashed",lwd=1)+
    
  #vertical line indicating start of fattening  
    #geom_vline(data=interp.data,aes(xintercept = interp.data$Fat_Start[1]), col="red", linetype="dashed",lwd=2)+
    
  #vertical line indicating start of fattening  
    #geom_vline(data=interp.data,aes(xintercept = interp.data$Fat_End[1]), col="red", linetype="dashed",lwd=1)+
  
  #original (plus extrapolated rule 2) body mass data points
    #geom_point(data=interp.data[interp.data$Test=="TRUE",], aes(x= Date_Index, y= Body_g_int),    col="black",lwd=2, alpha=1)+
    
    geom_point(data=qmr.final[qmr.final$Bird_ID==BIRD_ID& qmr.final$Morn_Eve=="MORN",],aes(x=Date_Index, y=Body_g, shape=Torpor_Use),col="black",lwd=4, alpha=1)+

      #geom_point(data=qmr.final[qmr.final$Bird_ID==BIRD_ID,],aes(x=Date_Index, y=Body_g),col="black",lwd=2, alpha=1)+
    
      #geom_point(data=qmr.final[qmr.final$Bird_ID==BIRD_ID,],aes(x=Date_Index, y=Body_g, shape=Torpor_Use),col="black",lwd=2, alpha=1)+

      #geom_point(data=qmr.final[qmr.final$Bird_ID==BIRD_ID& qmr.final$Morn_Eve=="MID",],aes(x=Date_Index, y=Body_g),col="green",lwd=2, alpha=1)+
      #geom_point(data=qmr.final[qmr.final$Bird_ID==BIRD_ID& qmr.final$Morn_Eve=="EVE",],aes(x=Date_Index, y=Body_g),col="blue",lwd=2, alpha=1)+
      scale_y_continuous(limits= c(2.3, 5), breaks= seq(0,5, by=.5), name="Body Mass \\n (g)")+
                       # sec.axis = sec_axis( trans=~.*2.5, breaks= seq(8,15, by=1),name="Night Duration (Hours)"))+
      scale_x_continuous(limits= c(40, 120), breaks= seq(0,120, by=10), name="Days, Since June 1")+ 
        #orig 42
    
    #add night duraiton (need to add second axis)
   # geom_point(aes(x=Date_Index, y=Night_Duration_Week/2.5), col="grey",fill="grey", shape=15,lwd=3, alpha=1)+

      scale_shape_manual(values=c("TORPOR"=16,"NORMO"=1))+
      labs(title=BIRD_ID)+
      #labs(caption = paste("Start of Fattening: Day", cutoff))+
      my_theme
plot

#FOR ALL BIRDS THE OPEN AND CLOSED CIRCLES FOR NORMO/ TOPROR ARE WORKING (edited after revisions)

#############################################################################
  ###########################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  
    #CALCULATE VARIABLES DURING THE FATTENING PERIOD
  
    ################################################################################
    ################################################################################
  
    #smooth original body mass data
    smoothmodel<-smooth.spline(qmr.final$Fat_g[qmr.final$Bird_ID==BIRD_ID & is.na(qmr.final$Fat_g)==F],
                          x=qmr.final$Date_Index[qmr.final$Bird_ID==BIRD_ID& is.na(qmr.final$Fat_g)==F],
                            spar=SPAR) 
    #predict daily body mass average based on the above spline model
    predmodel<- predict(smoothmodel, x=seq(1,122, length=122))
    interp.data$Fat_g_smooth<-predmodel$y #save to interp.data
    
    ################################################################################
    ################################################################################
    
    #if the bird DID FATTEN, then calculate these things, if not go to else
    if(unique(interp.data$Bird_ID) %in% c(c("B1","B3","B4","B6","B7","B9","B10","B11","B12","B13","B14","B15","B16"))) {
        
        #fattening start date
        Fat_Start<-interp.data$Fat_Start[1]
        #fattening end date
        Fat_End<-interp.data$Fat_End[1]
        #fattening period duration
        Fat_Duration<- Fat_End-Fat_Start
        
        ################################################################################
        
        #BODY MASS at the start of fattening
        Body_g_F_Start<-interp.data$Body_g_smooth[interp.data$Date_Index==Fat_Start]
        #BODY MASS at the end of fattening
        Body_g_F_End<-interp.data$Body_g_smooth[interp.data$Date_Index==Fat_End]
        #BODY MASS INCREASE DURING FATTENING
        Body_gi_F<-Body_g_F_End-Body_g_F_Start
        #Calculate % body mass increase compared to start
        Body_gpi_F<-Body_gi_F/Body_g_F_Start*100
        
        ################################################################################
        
        #FAT MASS at the start of fattening
        Fat_g_F_Start<-interp.data$Fat_g_smooth[interp.data$Date_Index==Fat_Start]
        #FAT MASS at the end of fattening
        Fat_g_F_End<-interp.data$Fat_g_smooth[interp.data$Date_Index==Fat_End]
        #FAT MASS INCREASE DURING FATTENING
        Fat_gi_F<-Fat_g_F_End-Fat_g_F_Start
        #Calculate % FAT mass increase compared to start
        Fat_gpi_F<-Fat_gi_F/Fat_g_F_Start*100
        
        ################################################################################
        
        #Fat content % at the start of fattening
        Fat_p_F_Start<-Fat_g_F_Start/Body_g_F_Start*100
        #Fat content % at the start of fattening
        Fat_p_F_End<-Fat_g_F_End/Body_g_F_End*100
        #Fat content % gained during  fattening
        Fat_pi_F<-Fat_p_F_End-Fat_p_F_Start
        #% FAT mass increase compared to start
        Fat_ppi_F<-Fat_pi_F/Fat_p_F_Start*100

    } #END IF THE BIRDS FATTENED 

    else {#IF BIRDS DIDNT FATTEN, fill in these variables with NA
        Fat_Start<-NA
        Fat_End<-NA
        Fat_Duration<- NA
        
        Body_g_F_Start<-NA
        Body_g_F_End<-NA
        Body_gi_F<-NA
        Body_gpi_F<-NA
        
        Fat_g_F_Start<-NA
        Fat_g_F_End<-NA
        Fat_gi_F<-NA
        Fat_gpi_F<-NA
        
        Fat_p_F_Start<-NA
        Fat_p_F_End<-NA
        Fat_pi_F<-NA
        Fat_ppi_F<-NA
      
    }#END OF IF BIRDS DIDNT FATTEN
    #END OF CALCULATING FAT INCREASES

  ################################################################################    
  ################################################################################
  ################################################################################
  ################################################################################
  ################################################################################
  
  #Combine the above generated variables into a bird specific seasonal cutoff dates summary vector
  changeindexs<-data.frame(BIRD_ID,Catch_Date,Fat_Start,Fat_End, Fat_Duration,
                            Body_g_F_Start,Body_g_F_End,Body_gi_F,Body_gpi_F,
                           Fat_g_F_Start,Fat_g_F_End,Fat_gi_F,Fat_gpi_F,
                           Fat_p_F_Start,Fat_p_F_End,Fat_pi_F,Fat_ppi_F)

  ##############################################################################
  ##############################################################################

  print(plot)
  
return(plot)  
#return(changeindexs)
}#END SMOOTH BODY FUNCTION

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################



##########
#NOW NEED TO ADD A separate panel plot OF THE PHOTOPERIODS. 
 #First create a new dataframe to populate #populate it with date_index from 1:122 (june 1 to sept 30)

interp.week<-data.frame(1:122)
names(interp.week)<-"Date_Index"
    
interp.week<-data.frame(1:122)
names(interp.week)<-"Date_Index"
    
#addd night duration
    
interp.week$Night_Duration<-NA

for(i in c(1:122)){
interp.week$Night_Duration[interp.week$Date_Index==i]<-mean(run.agg$Night_Duration[run.agg$Date_Index==i], na.rm=T)
}

    #GENERATE VECTORS OF INTERPOLATED**, 

    #interpolate
    interp.week$Night_Duration<-na.approx(interp.week$Night_Duration,x=interp.week$Date_Index,
                          xout=c(1:122), na.rm=F, rule=1)
    
    interp.week$Week<-floor(interp.week$Date_Index/7)
    interp.week$Week
    
    for(i in c(1:17)){
    interp.week$Night_Duration_Week[interp.week$Week==i]<-floor(mean(interp.week$Night_Duration[interp.week$Week==i], na.rm=T))
    }

  interp.week$Night_Duration_Week[interp.week$Week<3]<-525  
  interp.week$Night_Duration_Week<-interp.week$Night_Duration_Week/60
   
interp.week.week<-interp.week %>%
  group_by(Week) %>%
  #summarize (usually mean +- std error) variables that have multiple values per run (i.e. minute data)
  summarize(Night_Duration = mean(Night_Duration, na.rm = T))%>%
          data.frame()

#generate the night duraiton plot. 
night.duration.plot<-ggplot(interp.week.week)+
    geom_col(aes(x=Week*7, y=Night_Duration/60), col="grey", fill="grey")+
       scale_y_continuous(limits= c(8, 12), oob=scales::squish, breaks= seq(0,20, by=1), name="Night Length (Hours)")+
      scale_x_continuous(limits= c(38, 123), breaks= seq(0,120, by=10), name="Days, Since June 1")+ 
    my_theme

night.duration.plot

########################################################################################################
########################################################################################################
#RUN THE SMOOTH BODY MASS FUNCITON FOR EACH BIRD

#some notes:
#spar is always .35
#cutoff is 75% of the maximum slope.
#dont include the first day of fattening as fattening. 
#dont include the last day of fatteing as fattening.

#RESET WD TO SAVE FIGURES TO A SPECIFIC FOLDER.

setwd("~/ERICH_EBERTS/ACADEMICS/PROJECTS/FAT/FAT_WRITING/Figures")

#####
#night Duration
supfig7NL<-night.duration.plot
ggexport(supfig7NL, filename = "Sup Fig 7- Night Length vs Date.tiff",width = 650,height = 450)

#good
supfig7b1<-smooth_body_func("B1", SPAR=.35, CUTOFF=.03) #what about 75% of the max rate? here the max rate was 4
ggexport(supfig7b1, filename = "Sup Fig 7- B1. Body Mass vs Date.tiff",width = 650,height = 450)

#nonfattener
supfig7b2<-smooth_body_func("B2", SPAR=.35, CUTOFF=1) #no real change
ggexport(supfig7b2, filename = "Sup Fig 7- B2. Body Mass vs Date.tiff",width = 650,height = 450)

#good
supfig7b3<-smooth_body_func("B3", SPAR=.35, CUTOFF=.06) #max .08.....#i changed 9.12 from 2.86 to 3.86 and 9.16 from 2.87 to 3.87
#fudge 3 one day to include last point as fattening
ggexport(supfig7b3, filename = "Sup Fig 7- B3. Body Mass vs Date.tiff",width = 650,height = 450)

#weird cont fattener
#good now after a bit of messing with those flat points in the middle of fattening
supfig7b4<-smooth_body_func("B4", SPAR=.35, CUTOFF=.034) # max.045
ggexport(supfig7b4, filename = "Sup Fig 7- B4. Body Mass vs Date.tiff",width = 650,height = 450)

#nonfattener
supfig7b5<-smooth_body_func("B5", SPAR=.35, CUTOFF=1)#no  real change
ggexport(supfig7b5, filename = "Sup Fig 7- B5. Body Mass vs Date.tiff",width = 650,height = 450)

#good
supfig7b6<-smooth_body_func("B6", SPAR=.35, CUTOFF=.027) #max is 0.036
ggexport(supfig7b6, filename = "Sup Fig 7- B6. Body Mass vs Date.tiff",width = 650,height = 450)

#good
supfig7b7<-smooth_body_func("B7", SPAR=.35, CUTOFF=.036)  ##max 0.48 
ggexport(supfig7b7, filename = "Sup Fig 7- B7. Body Mass vs Date.tiff",width = 650,height = 450)

#nonfattener
supfig7b8<-smooth_body_func("B8", SPAR=.7, CUTOFF=.1) #no real change
ggexport(supfig7b8, filename = "Sup Fig 7- B8. Body Mass vs Date.tiff",width = 650,height = 450)

#good after fudging 
##############################################wtf this is messed up again
supfig7b9<-smooth_body_func("B9", SPAR=.35, CUTOFF=.049) #max 0.065 need to go to about 50% of max to get to a good value
ggexport(supfig7b9, filename = "Sup Fig 7- B9. Body Mass vs Date.tiff",width = 650,height = 450)

#good
supfig7b10<-smooth_body_func("B10", SPAR=.35, CUTOFF=.08) #max .12
ggexport(supfig7b10, filename = "Sup Fig 7- B10. Body Mass vs Date.tiff",width = 650,height = 450)
########3shit is that first point inclueded ? i dont think so?
#good after fudging

supfig7b11<-smooth_body_func("B11", SPAR=.35, CUTOFF=.079) #max .105  
ggexport(supfig7b11, filename = "Sup Fig 7- B11. Body Mass vs Date.tiff",width = 650,height = 450)
#dont want to include that last point but i think it is included 

#good
supfig7b12<-smooth_body_func("B12", SPAR=.35, CUTOFF=.053) # max 0.71
ggexport(supfig7b12, filename = "Sup Fig 7- B12. Body Mass vs Date.tiff",width = 650,height = 450)

#good
supfig7b13<-smooth_body_func("B13", SPAR=.35, CUTOFF=.075)#max.1
ggexport(supfig7b13, filename = "Sup Fig 7- B13. Body Mass vs Date.tiff",width = 650,height = 450)

#good after fudge
supfig7b14<-smooth_body_func("B14", SPAR=.35, CUTOFF=.075)
ggexport(supfig7b14, filename = "Sup Fig 7- B14. Body Mass vs Date.tiff",width = 650,height = 450)

#good after fudge
supfig7b15<-smooth_body_func("B15", SPAR=.35, CUTOFF=.079) #max .105 #fucked. 
ggexport(supfig7b15, filename = "Sup Fig 7- B15. Body Mass vs Date.tiff",width = 650,height = 450)

#good. 
supfig7b16<-smooth_body_func("B16", SPAR=.35, CUTOFF=.027) #max 0.036
ggexport(supfig7b16, filename = "Sup Fig 7- B16. Body Mass vs Date.tiff",width = 650,height = 450)

########################################################################################################
########################################################################################################

#add all of the separate bird's values into one dataframe.
seasons.data<-rbind(B1_sum,B2_sum,B3_sum,B4_sum,B5_sum,
      B6_sum,B7_sum,B8_sum,B9_sum,B10_sum,
      B11_sum,B12_sum,B13_sum,B14_sum,B15_sum,
      B16_sum)#b26 left out, B17_sum,B18_sum,B19_sum

########################################################################################################

#Labels the nonfatteners and the cont fatter's fat cutoff dates as NaA
#since they dont fatten. 
seasons.data$Fat_Start[is.infinite(seasons.data$Fat_Start)==T]<-NA
seasons.data$Fat_End[is.infinite(seasons.data$Fat_End)==T]<-NA
seasons.data$Fat_Duration[is.infinite(seasons.data$Fat_Duration)==T]<-NA

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

#CALCULATE SOME IMPORTANT SUMMARY STATISTICS 

################################################################################

#FATTENING DATES AND DURATION

#Fattening Period start date
Fat_Start_mean<-mean(seasons.data$Fat_Start, na.rm=T)
Fat_Start_sd<-std.error(seasons.data$Fat_Start, na.rm=T)

#Fattening Period fattening end date
Fat_End_mean<-mean(seasons.data$Fat_End, na.rm=T)
Fat_End_sd<-std.error(seasons.data$Fat_End, na.rm=T)

#Fattening Period  duration
Fat_Duration_mean<-mean(seasons.data$Fat_Duration, na.rm=T)
Fat_Duration_sd<-std.error(seasons.data$Fat_Duration, na.rm=T)
Fat_Duration_range<-range(seasons.data$Fat_Duration, na.rm=T)

#PRINT MEAN +- SE FATTENING PERIOD DURATION
Fat_Duration_mean
Fat_Duration_sd
Fat_Duration_range

################################################################################
################################################################################

#BODY MASS GAINS DURING FATTENING PERIOD

#Absolute (g)
Body_gi_F_mean<-mean(seasons.data$Body_gi_F, na.rm=T)
Body_gi_F_sd<-std.error(seasons.data$Body_gi_F, na.rm=T)
Body_gi_F_range<-range(seasons.data$Body_gi_F, na.rm=T)

#PRINT MEAN +- SE ABSOLUTE BODY MASS GAINS (g)
Body_gi_F_mean
Body_gi_F_sd
Body_gi_F_range

#BODY MASS GAINS
#Relative (%, relative to start of fattening period)
Body_gpi_F_mean<-mean(seasons.data$Body_gpi_F, na.rm=T)
Body_gpi_F_sd<-std.error(seasons.data$Body_gpi_F, na.rm=T)
Body_gpi_F_range<-range(seasons.data$Body_gpi_F, na.rm=T)

#PRINT MEAN +- SE RELATIVE BODY MASS GAINS (%)
Body_gpi_F_mean
Body_gpi_F_sd
Body_gpi_F_range

################################################################################
################################################################################

#FAT MASS GAINS

#Absolute (g)
Fat_gi_F_mean<-mean(seasons.data$Fat_gi_F, na.rm=T)
Fat_gi_F_sd<-std.error(seasons.data$Fat_gi_F, na.rm=T)
Fat_gi_F_range<-range(seasons.data$Fat_gi_F, na.rm=T)

#PRINT MEAN +- SE ABSOLUTE FAT MASS GAINS (g)
Fat_gi_F_mean
Fat_gi_F_sd
Fat_gi_F_range

#Relative (%, relative to start of fattening period)
Fat_gpi_F_mean<-mean(seasons.data$Fat_gpi_F, na.rm=T)
Fat_gpi_F_sd<-std.error(seasons.data$Fat_gpi_F, na.rm=T)
Fat_gpi_F_range<-range(seasons.data$Fat_gpi_F, na.rm=T)

#PRINT MEAN +- SE RELATIVE BODY MASS GAINS (%)
Fat_gpi_F_mean
Fat_gpi_F_sd
Fat_gpi_F_range

################################################################################
################################################################################

#FAT MASS GAINS
#Absolute (%)
Fat_pi_F_mean<-mean(seasons.data$Fat_pi_F, na.rm=T)
Fat_pi_F_sd<-std.error(seasons.data$Fat_pi_F, na.rm=T)
Fat_pi_F_range<-range(seasons.data$Fat_pi_F, na.rm=T)

#PRINT MEAN +- SE ABSOLUTE FAT CONTENT GAINS (%)
Fat_pi_F_mean
Fat_pi_F_sd
Fat_pi_F_range

#Calculate mean fat mass gain during fattening
Fat_ppi_F_mean<-mean(seasons.data$Fat_ppi_F, na.rm=T)
Fat_ppi_F_sd<-std.error(seasons.data$Fat_ppi_F, na.rm=T)
Fat_ppi_F_range<-range(seasons.data$Fat_ppi_F, na.rm=T)

#PRINT MEAN +- SE RELATIVE FAT CONTENT GAINS (%)
Fat_ppi_F_mean
Fat_ppi_F_sd
Fat_ppi_F_range

################################################################################
################################################################################
################################################################################
ggplot(seasons.data)+
  geom_point(aes(x=Body_gi_F))

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

#fix bird id's name
seasons.data$Bird_ID<-seasons.data$BIRD_ID

################################################################################

#rename rownames 
rownames(seasons.data)<-1:16

################################################################################
#CREATE THE SEASON CUTOFFS FILE
setwd("~/ERICH_EBERTS/ACADEMICS/PROJECTS/FAT/FAT_DATA/FAT_CODE/FINAL_DATA")

write.csv(seasons.data, "season_cutoffs.csv")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#
#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#THE END#