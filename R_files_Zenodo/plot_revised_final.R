#PLOT_FAT.R DATA

#LAST UPDATED 9/5/2021

#This script uses run_agg.csv, all_data_1min, period_cutoffs_final, food_data_plot, dataframes (csvs_) (generated in analysis.R)
#This script generates and saves the  figures for the Fat manuscript, which are finalized in powerpoint. 

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
library(MESS)
library(effects)
library(lme4)
library(lmerTest)
library(rmcorr)
library(emmeans)
library(viridis)

################################################################################

#Set the theme for any ggplot figures

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

#LOAD IN DATA

################################################################################

#LOAD IN run_agg.csv

run.agg<-read.csv("run_agg.csv")

run.agg$Season<-factor(run.agg$Season,levels=c("Summer","Fattening","Migration","Non-Fattener"))

run.agg$Bird_ID<-factor(run.agg$Bird_ID,levels=c("B1","B2","B3","B4","B5", "B6","B7","B8","B9","B10","B11","B12","B13","B14","B15","B16"))

################################################################################
################################################################################

#LOAD IN period_cutoffs_final.csv
period.cutoffs.final<-read.csv("period_cutoffs_final.csv")

################################################################################
################################################################################
#LOAD IN all_data_1min_csv

#this will be used for instantaneous fat vs time plots (FIG S6)
all.data.1min<-read.csv("all_data_1min.csv")

all.data.1min$Bird_ID<-factor(all.data.1min$Bird_ID,levels=c("B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12","B13","B14","B15","B16"))

all.data.1min$Season<-factor(all.data.1min$Season,levels=c("Summer","Fattening","Migration","Non-Fattener"))

all.data.1min$Fat_p[all.data.1min$Torpor_Use=="NORMO"]


################################################################################
################################################################################
#LOAD IN food_data_plot

food.data<-read.csv("food_data_plot.csv")

food.data$Bird_ID<-factor(food.data$Bird_ID,levels=c("B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12","B13","B14","B15","B16"))

food.data$Season<-factor(food.data$Season, levels=c("Summer","Fattening","Migration", "Non-Fattener"))

################################################################################
################################################################################
#################################################################################
#################################################################################
#################################################################################

#RESET WD TO SAVE FIGURES TO A SPECIFIC FOLDER.

setwd("~/ERICH_EBERTS/ACADEMICS/PROJECTS/FAT/FAT_WRITING/Figures")

################################################################################
################################################################################
################################################################################
#END DATA SETUP
unique(run.agg$Run_ID)
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

#MAIN FIGURES

################################################################################
################################################################################

#Figure 1
#BODY MASS- SEASONAL SLOPES VS DATE

fig1<-ggplot(data=run.agg[run.agg$Season %in% c("Summer","Fattening","Migration","Non-Fattener","Cont-Fattener"),])+

  geom_point(aes(x= Date_Index ,y=Body_Morn_g, col=Season), lwd=3, alpha=.6)+
    #geom_point(aes(x= Date_Index ,y=Body_Morn_g, col=Bird_ID), lwd=3, alpha=.6)+

  #geom_smooth(method="lm",aes(x=Date_Index,y=Body_Morn_g, col=Season),se=F)+
  
  scale_color_manual(values=c("#56B4E9", "#F0E442","red","black", "grey"))+
  
  scale_y_continuous(limits = c(2.25,4.5),breaks = seq(0, 5, by = .5), name = "Morning Body Mass (g)")+
  scale_x_continuous(limits = c(0,125),breaks = seq(0, 125, by = 30),name="Day, Since June 1")+

  facet_wrap(~Bird_ID)+

  my_theme

fig1
ggexport(fig1, filename = "Fig 1. Body Mass vs Date.tiff",width = 650,height = 600)

################################################################################
################################################################################
################################################################################
################################################################################

#FIGURE 2 A

#TORPOR DURATION- SEASONAL SLOPES VS FATUSING DURATION AS % OF NIGHT

fig2a<-ggplot(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"& run.agg$Season %in% c("Summer"),])+

  geom_smooth(method="lm",aes(x=Fat_Eve_p, y=Torpor_Duration_ET/60, col=Season),se=F)+
  geom_point(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"& run.agg$Season %in% c("Fattening","Summer", "Migration"),],aes(x=Fat_Eve_p, y=Torpor_Duration_ET/60, color=Season, shape=Torpor_Use), alpha=.6, lwd=3)+

  geom_point(data=run.agg[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Fattening","Summer", "Migration"),],aes(x=Fat_Eve_p, y=Torpor_Duration_ET/60, color=Season, shape=Torpor_Use), alpha=.6, lwd=3)+

  scale_color_manual(values=c("#F0E442","red", "#56B4E9","black", "grey"))+
  scale_shape_manual(values = c(1,16,17,18,19))+
  
  scale_y_continuous(limits = c(0,9),breaks = seq(0, 10, by = 1), name = "Torpor Duration (hrs)")+
  scale_x_continuous(limits = c(0,40.4),breaks = seq(0, 100, by = 10), name="Evening Fat Content (%)")+
  my_theme

fig2a
ggexport(fig2a, filename = "Fig 2A. Torpor Duration vs Eve Fat.tiff",width = 650,height = 450)

################################################################################

#FIGURE 2 B

#Overnight mass loss vs evening fat 

fig2b<-ggplot(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"& run.agg$Season %in% c("Summer"),])+

  geom_smooth(method="lm",aes(x=Fat_Eve_p, y=Fat_delta_g, col=Season),se=F)+
  geom_point(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"& run.agg$Season %in% c("Fattening","Summer", "Migration"),],aes(x=Fat_Eve_p, y=Fat_delta_g, color=Season, shape=Torpor_Use), alpha=.6, lwd=3)+

  geom_point(data=run.agg[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Fattening","Summer", "Migration"),],aes(x=Fat_Eve_p, y=Fat_delta_g, color=Season, shape=Torpor_Use), alpha=.6, lwd=3)+

  scale_color_manual(values=c("#F0E442","red", "#56B4E9","black", "grey"))+
  scale_shape_manual(values = c(1,16,17,18,19))+
  
  scale_y_continuous(limits = c(0,.25),breaks = seq(0, 10, by = .05), name = "Overnight Fat Mass Loss (g)")+
  scale_x_continuous(limits = c(0,40.4),breaks = seq(0, 100, by = 10), name="Evening Fat Content (%)")+
  my_theme

fig2b
ggexport(fig2b, filename = "Fig 2B. Fat Loss vs Eve Fat.tiff",width = 650,height = 450)

###################################################################
#############################################################

#normothermic fat mass losses
##
mean(run.agg$Fat_delta_g[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Summer")], na.rm=T)
std.error(run.agg$Fat_delta_g[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Summer")], na.rm=T)

mean(run.agg$Fat_delta_g[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Fattening")], na.rm=T)
std.error(run.agg$Fat_delta_g[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Fattening")], na.rm=T)

mean(run.agg$Fat_delta_g[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Migration")], na.rm=T)
std.error(run.agg$Fat_delta_g[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Migration")], na.rm=T)

################################################################################
################################################################################
################################################################################
################################################################################

#FIGURE 3A

#TORPOR ENTRY FAT LEVELS- SEASONAL SLOPES VS TIME OF ENTRY
fig3a<-ggplot(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"&run.agg$Season %in% c("Migration"),])+

  #geom_smooth(method="lm",aes(x=Crit_Time_p, y=Crit_Fat_p, group=Bird_ID),col="grey", linetype=3,se=F)+
  geom_smooth(method="lm",aes(x=Crit_Time_p, y=Crit_Fat_p, group=Season),col="red", linetype=1,se=F)+
  
  geom_point(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"&run.agg$Season %in% c("Summer","Fattening","Migration"),],aes(x=Crit_Time_p, y=Crit_Fat_p, color=Season), alpha=.6, lwd=3)+

  geom_point(data=run.agg[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Fattening","Summer", "Migration"),],aes(x=Crit_Time_p+5, y=Crit_Fat_p, color=Season, shape=Torpor_Use), alpha=.6, lwd=3)+

  scale_color_manual(values=c("#56B4E9","#F0E442", "red","black", "grey"))+
  scale_shape_manual(values = c(1,16,17,18,19))+

  scale_y_continuous(limits = c(0,41),breaks = seq(0, 100, by = 10), name = "Fat Content at Torpor Entry (%)")+
  scale_x_continuous(limits = c(0,105),breaks = seq(0, 100, by = 20), name="Time of Night at Torpor Entry (%)")+
  
  my_theme

fig3a
ggexport(fig3a, filename = "Fig 3A. Threshold vs Time.tiff",width = 650,height = 450)

################################################################################

#FIGURE 3B
#TORPOR ENTRY FAT LEVELS- SEASONAL SLOPES VS DATE_INDEX

fig3b<-ggplot(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"&run.agg$Season %in% c("Summer","Fattening","Migration"),])+
#  geom_smooth(method="lm",aes(x=Date_Index, y=Crit_Fat_p, col=Season),se=F)+
 
  geom_point(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"&run.agg$Season %in% c("Summer","Fattening","Migration"),],aes(x=Date_Index, y=Crit_Fat_p, color=Season), alpha=.6, lwd=3)+
#  geom_smooth(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"&run.agg$Season %in% c("Summer","Fattening","Migration"),],method="lm",aes(x=Date_Index, y=Crit_Fat_p, col=Season),se=F)+

  geom_point(data=run.agg[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Fattening","Summer", "Migration"),],aes(x=Date_Index, y=Crit_Fat_p, color=Season, shape=Torpor_Use), alpha=.6, lwd=3)+

  scale_color_manual(values=c("#56B4E9","#F0E442","red","black"))+  
  scale_shape_manual(values = c(1,16,17,18,19))+

  scale_y_continuous(limits = c(0,40),breaks = seq(0, 100, by = 10), name = "Fat Content at at Torpor Entry (%)")+
  scale_x_continuous(limits = c(0,125),breaks = seq(0, 200, by = 30), name="Date")+
  
  my_theme

fig3b
ggexport(fig3b, filename = "Fig 3B. Thresholds vs Date.tiff",width = 650,height = 450)

################################################################################
################################################################################
################################################################################
################################################################################

#FIGURE 4 A

#WITHIN FATTENING PERIOD (use period cutoffs final)
#BODY MASS GAINS VS TORPOR DURATION 
#AND FAT DURATOIN VS TPROR DUARTION

coef<-20 #CONVERSION FACTOR FOR SECOND SCALE
fig4a<-ggplot(data=period.cutoffs.final)+
  #body mass gains on left axis
  geom_point(aes(x= Fattening_TorDur_mean,y=Body_gi_F), lwd=3, alpha=.6)+
  geom_smooth(aes(x=Fattening_TorDur_mean,y=Body_gi_F),col="black", lwd=1, alpha=.1, method="lm",se=T)+
  #fattenign duration on right axis
  geom_point(aes(x= Fattening_TorDur_mean,y=Fat_Duration/coef), col="#56B4E9",lwd=3, alpha=.6)+  
  #geom_smooth(aes(x=Fattening_TorDur_mean_h,y=Fat_Duration/coef),col="#56B4E9", lwd=1, alpha=.1, method="lm",se=T)+
  
  scale_y_continuous(limits = c(0,1),breaks = seq(0, 1, by = .2), name = "Body Mass Gains \\n During Fattening Period (g)",
                    sec.axis = sec_axis( trans=~.*coef, name="Fattening Period Duration (days)"))+
  
  scale_x_continuous(limits = c(-.5,6),breaks = seq(0, 10, by = 1),name="Mean Torpor Duration During Fattening Period (hrs)")+
  
  my_theme+
  theme( axis.line.y.right = element_line(color = "#56B4E9"), axis.ticks.y.right = element_line(color = "#56B4E9"),
         axis.text.y.right = element_text(color = "#56B4E9"),axis.title.y.right = element_text(color = "#56B4E9"))+
  theme(axis.title.y.right = element_text(vjust= 4,size=20))

fig4a
ggexport(fig4a, filename = "Fig 4A. Fattening Gains and Duration vs Torpor.tiff",width = 650,height = 450)

################################################################################

#FIGURE 4 B

#WITHIN FATTENING PERIOD (use period cutoffs final)
#BODY MASS GAINS VS FOOD CONSUMPTION
#AND FAT DURATOIN VS FOOD CONSUMPTION

coef2<-20
fig4b<-ggplot(data=period.cutoffs.final)+
  #body mass gains
  geom_point(aes(x= Fattening_Food_mean,y=Body_gi_F),col="black", lwd=3, alpha=.6)+
  #geom_smooth(aes(x=Fattening_Food_mean,y=Body_gi_F),col="black", lwd=1, alpha=.1, method="lm",se=T)+
  
  #fattening duration
  geom_point(aes(x= Fattening_Food_mean,y=Fat_Duration/coef2), col="#56B4E9", lwd=3, alpha=.6)+
  #geom_smooth(aes(x=Fattening_Food_mean,y=Fat_Duration/coef2),col="#56B4E9", lwd=1, alpha=.1, method="lm",se=T)+
  
  scale_y_continuous(limits = c(0,1),breaks = seq(0, 1, by = .2), name = "Body Mass Gains \\n During Fattening Period (g)",
                     sec.axis = sec_axis( trans=~.*coef2, name="Fattening Period Duration (days)"))+
  
  scale_x_continuous(limits = c(14,25),breaks = seq(0, 25, by = 2),name="Mean Daily Food Consumption \\n Duration During Fattening Period (mL)")+
  
  my_theme+
  theme( axis.line.y.right = element_line(color = "#56B4E9"), axis.ticks.y.right = element_line(color = "#56B4E9"),
         axis.text.y.right = element_text(color = "#56B4E9"),axis.title.y.right = element_text(color = "#56B4E9"))+
  theme(axis.title.y.right = element_text(vjust= 4,size=20))



fig4b
ggexport(fig4b, filename = "Fig 4B. Fattening Gains and Duration vs Food.tiff",width = 650,height = 450)


################################################################################
#END MAIN FIGURES
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

#SUPPLEMENTARY FIGURES

################################################################################
################################################################################

#SUP FIGURE 1 A
#FAT VS BODY CORRELATION

supfig1a<-ggplot(data=run.agg[run.agg$Season %in% c("Summer","Fattening","Migration","Non-Fattener"),])+
  geom_smooth(method="lm",aes(x=Body_Morn_g,y=Fat_Morn_g,group=Bird_ID),col="grey", linetype=3,lwd=.5,se=F,alpha=.1)+
  geom_smooth(method="lm",aes(x=Body_Morn_g,y=Fat_Morn_g),se=F, col="black",lwd=1,alpha=.1)+
  #geom_smooth(method="lm",aes(x=Body_Morn_g,y=Fat_Morn_g,col=Season),lwd=.5,se=F,alpha=.1)+
  geom_point(aes(x= Body_Morn_g ,y=Fat_Morn_g, col=Season), lwd=3, alpha=.6, shape="circle")+
  
  scale_color_manual(values=c("#56B4E9", "#F0E442","red","black", "grey"))+
  
  scale_y_continuous(limits = c(0,1.7),breaks = seq(0, 5, by = .5), name = "Morning Fat Mass (g)")+
  scale_x_continuous(limits = c(2.4,4.4),breaks = seq(0, 5, by = .5),name="Morning Body Mass (g)")+

  my_theme 

supfig1a
ggexport(supfig1a, filename = "Sup Fig 1A. Fat vs Body.tiff",width = 650,height = 450)

################################################################################

#SUP FIGURE 1 B
#LEAN VS BODY CORRELATION
supfig1b<-ggplot(data=run.agg[run.agg$Season %in% c("Summer","Fattening","Migration","Non-Fattener"),])+
  geom_point(aes(x= Body_Morn_g ,y=Lean_Morn_g, col=Season), lwd=3, alpha=.6, shape="circle")+
  geom_smooth(method="lm",aes(x=Body_Morn_g,y=Lean_Morn_g,group=Bird_ID),col="grey", linetype=3,lwd=.5,se=F,alpha=.1)+  
  geom_smooth(method="lm",aes(x=Body_Morn_g,y=Lean_Morn_g),col="black", lwd=1,se=F,alpha=.1)+  

  scale_color_manual(values=c("#56B4E9", "#F0E442","red","black", "grey"))+  
  
  scale_y_continuous(limits = c(0,3.5),breaks = seq(0, 5, by = .5), name = "Morning Lean Mass (g)")+
  scale_x_continuous(limits = c(2.4,4.4),breaks = seq(0, 5, by = .5),name="Morning Body Mass (g)")+
  
  my_theme

supfig1b
ggexport(supfig1b, filename = "Sup Fig 1B. Lean vs Body.tiff",width = 650,height = 450)

################################################################################
################################################################################
################################################################################
################################################################################

#SUP FIGURE 2
#TORPOR USE LOGISTIC VS EVENING FAT

supfig2<-ggplot(data=run.agg[run.agg$Season %in% c("Summer","Migration"),])+ 
  #geom_point(aes(x=Fat_Eve_p, y=Torpor_Use_bin, color=Season), alpha=.6, lwd=3)+
  geom_smooth(method="glm",aes(x=Fat_Eve_p, y=Torpor_Use_bin, color=Season),method.args=list(family="binomial"),se=F,fullrange=T)+

  geom_point(data=run.agg[run.agg$Season %in% c("Summer","Fattening","Migration"),],aes(x=Fat_Eve_p, y=Torpor_Use_bin, color=Season,shape=Torpor_Use), alpha=.6, lwd=3)+
  #geom_smooth(data=run.agg[run.agg$Season %in% c("Summer","Fattening","Migration"),],method="glm",aes(x=Fat_Eve_p, y=Torpor_Use_bin, color=Season),method.args=list(family="binomial"),se=F,fullrange=T)+

  scale_color_manual(values=c("#F0E442","red", "#56B4E9","black", "grey"))+
  scale_shape_manual(values = c(1,16,17,18,19))+
  
  scale_y_continuous(limits = c(0,1),breaks = seq(0, 5, by = .25), name = "Probability of Torpor Use")+
  scale_x_continuous(limits = c(0,40),breaks = seq(0, 100, by = 10), name="Evening Fat Content (%)")+
  
  my_theme

supfig2
ggexport(supfig2, filename = "Sup Fig 2. Torpor Log vs Eve Fat.tiff",width = 650,height = 450)

################################################################################
################################################################################
################################################################################
################################################################################

#SUP FIGURE 3 A
#PRETORPOR EXPENDIUTURE VS EVENING FAT

#DECIDE IF I SHOULD INCLUDE OPEN CIRCLES FOR NORMOTHERMIC NIGHTS

supfig3a<-ggplot(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"&run.agg$Season %in% c("Summer"),])+
  geom_point(aes(x=Fat_Eve_p, y=Fat_Loss_pretorpor_kJ, color=Season), alpha=.6, lwd=3)+
  geom_smooth(method="lm",aes(x=Fat_Eve_p, y=Fat_Loss_pretorpor_kJ, col=Season),se=F)+
  
  geom_point(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"&run.agg$Season %in% c("Fattening","Migration"),],aes(x=Fat_Eve_p, y=Fat_Loss_pretorpor_kJ, color=Season), alpha=.6, lwd=3)+
  #geom_smooth(method="lm",aes(x=Fat_Eve_p, y=Fat_Loss_pretorpor_kJ, col=Season),se=F)+
  #if want to include normothermic nights turn this on.
  geom_point(data=run.agg[run.agg$Torpor_Use=="NORMO",],aes(x=Fat_Eve_p, y=Expend_Total_kJ, color=Season), alpha=.6, lwd=3, shape=1)+

  scale_color_manual(values=c("#F0E442","red", "#56B4E9", "#56B4E9","black", "grey"))+
  scale_shape_manual(values = c(16,1,17,18,19))+

  scale_y_continuous(limits = c(0,10),breaks = seq(0, 10, by = 1), name = "Pretorpor Fat Expenditure (kJ)")+
  scale_x_continuous(limits = c(0,41),breaks = seq(0, 100, by = 10), name="Evening Fat Content (%)")+
  
  #if include normothermic nights, y scale max is 10, if not max is 7
  my_theme

supfig3a
ggexport(supfig3a, filename = "Sup Fig 3A. PreExpend vs Eve Fat.tiff",width = 650,height = 450)

################################################################################
################################################################################
################################################################################
################################################################################

#SUP FIGURE 3 B
#TIME OF TORPOR ENTRY VS EVENING FAT

supfig3b<-ggplot(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"& run.agg$Season %in% c("Summer"),])+
  geom_point(aes(x=Fat_Eve_p, y=Crit_Time_p, color=Season), alpha=.6, lwd=3)+
  geom_smooth(method="lm",aes(x=Fat_Eve_p, y=Crit_Time_p, col=Season),se=F)+

  geom_point(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"& run.agg$Season %in% c("Summer","Fattening","Migration"),],aes(x=Fat_Eve_p, y=Crit_Time_p, color=Season), alpha=.6, lwd=3)+
#  geom_smooth(method="lm",aes(x=Fat_Eve_p, y=Crit_Time_p, col=Season),se=F)+
  geom_point(data=run.agg[run.agg$Torpor_Use=="NORMO",],aes(x=Fat_Eve_p, y=105, color=Season), alpha=.6, lwd=3, shape=1)+

  scale_color_manual(values=c("#F0E442","red", "#56B4E9", "#56B4E9","black", "grey"))+
  scale_shape_manual(values = c(16,1,17,18,19))+
  
  scale_y_continuous(limits = c(0,105),breaks = seq(0, 100, by = 20), name = "Torpor Entry Time (% of night)")+
  scale_x_continuous(limits = c(0,41),breaks = seq(0, 100, by = 10), name="Evening Fat Content (%)")+
  my_theme

supfig3b
ggexport(supfig3b, filename = "Sup Fig 3B. Entry Time vs Eve Fat.tiff",width = 650,height = 450)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#NORMOTHERMIC NIGHTS ARENT WORKING.

#SUP FIGURE 4
#OVERNIGHT FAT MASS LOSS VS TOPROR DURATION

supfig4<-ggplot(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"& run.agg$Season %in% c("Summer","Fattening","Migration"),])+

  geom_smooth(method="lm",aes(x=Torpor_Duration_ET/60, y=Fat_delta_g, col=Season),se=F)+
  geom_point(data=run.agg[run.agg$Torpor_Use=="TORPOR"&run.agg$Torpor_Type!="Dip"&run.agg$Torpor_Type!="Double"& run.agg$Season %in% c("Fattening","Summer", "Migration"),],aes(x=Torpor_Duration_ET/60, y=Fat_delta_g, color=Season, shape=Torpor_Use), alpha=.6, lwd=3)+
  geom_smooth(method="lm",aes(x=Torpor_Duration_ET/60, y=Fat_delta_g), col="black",se=F)+

  geom_point(data=run.agg[run.agg$Torpor_Use=="NORMO" & run.agg$Season %in% c("Summer","Fattening", "Migration"),],aes(x=Torpor_Duration_ET/60, y=Fat_delta_g, color=Season, shape=Torpor_Use), alpha=.6, lwd=3)+
  #geom_smooth(method="lm",aes(x=Torpor_Duration_ET/60, y=Fat_delta_g, col=Season),se=F)+

  scale_color_manual(values=c("#56B4E9","#F0E442","red","black", "grey"))+
  scale_shape_manual(values = c(1,16,17,18,19))+
  
  scale_y_continuous(limits = c(0,.25),breaks = seq(0, 1, by = .05), name="Overnight Fat Mass Loss (g)")+
  scale_x_continuous(limits = c(0,9),breaks = seq(0, 10, by = 1), name = "Torpor Duration (hrs)")+
  
  my_theme

supfig4
ggexport(supfig4, filename = "Sup Fig 4. Mass Loss vs Duration.tiff",width = 650,height = 450)

################################################################################
################################################################################
################################################################################
################################################################################

# #SUP FIGURE 5A
# #BODY MASS- SEASONAL SLOPES VS DATE

supfig5a<-ggplot(data=run.agg[run.agg$Season %in% c("Fattening"),])+

  geom_smooth(method="lm",aes(x=Date_Index,y=Body_Morn_g, col=Season),se=F)+

  geom_point(data=run.agg[run.agg$Season %in% c("Summer","Fattening","Migration","Non-Fattener"),],
             aes(x= Date_Index ,y=Body_Morn_g, col=Season), lwd=3, alpha=.6)+

  scale_color_manual(values=c("#F0E442","red","black","#56B4E9"))+

  scale_y_continuous(limits = c(2.25,4.5),breaks = seq(0, 5, by = .5), name = "Morning Body Mass (g)")+
  scale_x_continuous(limits = c(0,125),breaks = seq(0, 125, by = 30),name="Day, Since June 1")+

  #facet_wrap(~Bird_ID)+

  my_theme

supfig5a
ggexport(supfig5a, filename = "Sup Fig 5A. Body Mass vs Date .tiff",width = 650,height = 450)

################################################################################

#SUP FIGURE 5B
#FAT CONTENT (P)- SEASONAL SLOPES VS DATE
supfig5b<-ggplot(data=run.agg[run.agg$Season %in% c("Fattening"),])+
  geom_smooth(method="lm",aes(x=Date_Index,y=Fat_Morn_p, col=Season),se=F)+
  
  geom_point(data=run.agg[run.agg$Season %in% c("Summer","Fattening","Migration","Non-Fattener"),],
             aes(x= Date_Index ,y=Fat_Morn_p, col=Season), lwd=3, alpha=.6)+
  
  scale_color_manual(values=c("#F0E442","red","black","#56B4E9"))+
  
  scale_y_continuous(limits = c(0,40),breaks = seq(0, 100, by = 10), name = "Morning Fat Content (%)")+
  scale_x_continuous(limits = c(0,125),breaks = seq(0, 125, by = 30),name="Date")+

  #facet_wrap(~Bird_ID)+
  
  my_theme

supfig5b
ggexport(supfig5b, filename = "Sup Fig 5B. Fat Content vs Date.tiff",width = 650,height = 450)

#################################################################################

#SUP FIGURE 5C
#LEAN MASS- SEASONAL SLOPES VS DATE

supfig5c<-ggplot(data=run.agg[run.agg$Season %in% c("Migration") &run.agg$Lean_Quality=="GOOD",])+
  geom_point(aes(x= Date_Index ,y=Lean_Morn_g, col=Season), lwd=3, alpha=.6)+
  geom_smooth(method="lm",aes(x=Date_Index,y=Lean_Morn_g, col=Season),se=F)+

  geom_point(data=run.agg[run.agg$Season %in% c("Summer","Fattening","Migration","Non-Fattener") &run.agg$Lean_Quality=="GOOD",],aes(x= Date_Index ,y=Lean_Morn_g, col=Season), lwd=3, alpha=.6)+  
  
  scale_color_manual(values=c("#F0E442","red","black","#56B4E9"))+  
  
  scale_y_continuous(limits = c(0,3.5),breaks = seq(0, 5, by = .5), name = "Morning Lean Mass (g)")+
  scale_x_continuous(limits = c(0,125),breaks = seq(0, 125, by = 30),name="Date")+
  
  #facet_wrap(~Bird_ID)+
  
  my_theme
supfig5c
ggexport(supfig5c, filename = "Sup Fig 5C. Lean Mass vs Date.tiff",width = 650,height = 450)

#################################################################################

#SUP FIGURE 5D
#FOOD INTAKE- SEASONAL SLOPES VS DATE
#i wonder if i should put this into 3 day avgs or something cuz theres so many frigin points

supfig5d<-ggplot(data=food.data[food.data$Season %in% c("Migration"),])+

  geom_smooth(method="lm",aes(x=Date_Index,y=Food_Intake, col=Season),se=F)+
  
  geom_point(data=food.data[food.data$Season %in% c("Summer","Fattening","Migration","Non-Fattener"),],aes(x= Date_Index ,y=Food_Intake, col=Season), lwd=3, alpha=.2)+
  
  scale_color_manual(values=c("#F0E442","red","black","#56B4E9"))+  
  scale_y_continuous(limits = c(0,31),breaks = seq(0, 35, by = 5), name = "Daily Food Consumption (mL)")+
  scale_x_continuous(limits = c(0,125),breaks = seq(0, 125, by = 30),name="Date")+
  
  #facet_wrap(~Bird_ID)+
  
  my_theme

supfig5d
ggexport(supfig5d, filename = "Sup Fig 5D. Food Consumption vs Date.tiff",width = 650,height = 450)

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#SUP FIGURE 6
#INST FAT VS TIME OF NIGHT, BIG OLD FACET PLOT

#Set the theme for following ggplot figures

my_theme <- theme_classic(base_size = 14) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  theme(axis.text.x = element_text(angle=0, size=16), legend.key.height = unit(3, 'lines'),
        panel.grid.major.y = element_line(color='white'),
        panel.grid.minor.y = element_line(color='white'))+
  theme(axis.text.y = element_text(angle=0, size=16), legend.key.height = unit(3, 'lines'),
        panel.grid.major.x = element_line(color='white'))+
  theme(plot.caption = element_text(size = 10, hjust = 0.5))+
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(vjust=-3, size=20))+
  theme(axis.title.y = element_text(vjust= 4,size=20))+
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))+
  theme(legend.position = "none") 
  #for facet plot labels
  # theme(strip.text.x = element_text(size=8, color="black",
  #                                     face="plain", hjust=0))+
  #         
  # theme(strip.background = element_rect(colour="black", fill="white", 
  #                                      size=.5, linetype="solid"))

  #theme(strip.background = element_blank(), strip.text.x = element_blank())

################################################################################


BIRD_ID<-"B1"

plot(all.data.1min$Fat_p[all.data.1min$Torpor_Use=="NORMO"&all.data.1min$Bird_ID==BIRD_ID])
(all.data.1min$Time_p[all.data.1min$Torpor_Use=="NORMO"&all.data.1min$Bird_ID==BIRD_ID])

plotfatp<-function(BIRD_ID){

  plot<-ggplot(data=all.data.1min[all.data.1min$Bird_ID==BIRD_ID&is.na(all.data.1min$Bird_ID)==F,])+
    #geom_point(aes(x=Time_p,y=Fat_p,group=Run_ID, col=Torpor_Use),lwd=.1,alpha=.4)+
    geom_line(aes(x=Time_p,y=Fat_p, group=Run_ID, col=Torpor_Use), lwd=1)+
    scale_color_manual(values=c("NORMO"="red","TORPOR"="blue"))+
    
    geom_hline(aes(yintercept=5.5), linetype=2)+ #mean summer threshold
  
    geom_hline(aes(yintercept=5.5+.8), linetype=2, col="grey")+ #mean summer threshold
    geom_hline(aes(yintercept=5.5-.8), linetype=2, col="grey")+ #mean summer threshold
  
    geom_hline(aes(yintercept=33), linetype=2)+ #mean summer threshold
    geom_hline(aes(yintercept=33+.8), linetype=2, col="grey")+ #mean summer threshold
    geom_hline(aes(yintercept=33-.8), linetype=2, col="grey")+ #mean summer threshold
  
    facet_wrap(~Season,nrow=1,ncol=3)+
    
    scale_y_continuous(limits=c(0,45), breaks = seq(from = 0, to = 50, by = 5), name= "Instantaneous Fat Content (%)")+
    scale_x_continuous(limits=c(0,100), breaks = seq(from = 0, to = 200, by = 25), name="Time of Night (%)")+
    labs(title=BIRD_ID)+
    my_theme
  print(plot)
return(plot)
}

B1_instfat<-plotfatp("B1") 
B2_instfat<-plotfatp("B2") #nonfattener
B3_instfat<-plotfatp("B3") 
B4_instfat<-plotfatp("B4") 
B5_instfat<-plotfatp("B5") #nonfattener
B6_instfat<-plotfatp("B6") 
B7_instfat<-plotfatp("B7") 
B8_instfat<-plotfatp("B8") #nonfatter
B9_instfat<-plotfatp("B9") 
B10_instfat<-plotfatp("B10") 
B11_instfat<-plotfatp("B11") 
B12_instfat<-plotfatp("B12") 
B13_instfat<-plotfatp("B13") 
B14_instfat<-plotfatp("B14") 
B15_instfat<-plotfatp("B15") 
B16_instfat<-plotfatp("B16")

ggexport(B1_instfat, filename = "Sup Fig 6-B1. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B2_instfat, filename = "Sup Fig 6-B2. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B3_instfat, filename = "Sup Fig 6-B3. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B4_instfat, filename = "Sup Fig 6-B4. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B5_instfat, filename = "Sup Fig 6-B5. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B6_instfat, filename = "Sup Fig 6-B6. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B7_instfat, filename = "Sup Fig 6-B7. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B8_instfat, filename = "Sup Fig 6-B8. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B9_instfat, filename = "Sup Fig 6-B9. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B10_instfat, filename = "Sup Fig 6-B10. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B11_instfat, filename = "Sup Fig 6-B11. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B12_instfat, filename = "Sup Fig 6-B12. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B13_instfat, filename = "Sup Fig 6-B13. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B14_instfat, filename = "Sup Fig 6-B14. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B15_instfat, filename = "Sup Fig 6-B15. Inst Fat vs Time.tiff",width = 650,height = 450)
ggexport(B16_instfat, filename = "Sup Fig 6-B16. Inst Fat vs Time.tiff",width = 650,height = 450)



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
#END OF SUPPLEMENTARY PLOTS
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################


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