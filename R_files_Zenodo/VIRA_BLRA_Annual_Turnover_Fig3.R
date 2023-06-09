# THIS FILE GRAPHS ANNUAL COVARIATES (+1 standard error) FOR BLACK (BLRA) AND VIRGINIA (VIRA) RAILS shown in Figure 3


#Load some packages
library(ggplot2)
library(cowplot)

#setwd("") to location with RDS files


#Get BLRA data
BCOL_ALL <- readRDS("BLRA_COL_ALL.RDS")
BEXT_ALL <- readRDS("BLRA_EXT_ALL.RDS")
str(BCOL_ALL)
str(BEXT_ALL)


# Get VIRA data
VCOL_ALL <- readRDS("C:/Users/beis/Box Sync/Black_Rails/Extinction Dynamics/2021/Analysis/VIRA/VIRA_COL_ALL.RDS")
VEXT_ALL <- readRDS("C:/Users/beis/Box Sync/Black_Rails/Extinction Dynamics/2021/Analysis/VIRA/VIRA_EXT_ALL.RDS")
str(VCOL_ALL)
str(VEXT_ALL)



#VIRA Graphs

VWNVcur_Precip_Col  <- ggplot(VCOL_ALL, aes(x=stdvalue)) +
  geom_line(aes(y=Col.Precip.mean), color = "blue", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Col.Precip.meanminus1SE, ymax = Col.Precip.meanplus1SE), alpha = 0.3, colour = "blue", fill = "blue", inherit.aes = TRUE) +
  geom_line(aes(y=Col.WNVcur_Precip.plus1SD), color = "orange", linetype = "dashed", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Col.WNVcur_Precip.plus1SD.minus1SE, ymax = Col.WNVcur_Precip.plus1SD.plus1SE), alpha = 0.3, colour = "gray",  inherit.aes = TRUE) +
  geom_line(aes(y=Col.WNVcur_Precip.minus1SD), color = "orange", linetype = "dotted", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Col.WNVcur_Precip.minus1SD.minus1SE, ymax = Col.WNVcur_Precip.minus1SD.plus1SE), alpha = 0.3, colour = "gray",  inherit.aes = TRUE)+
  coord_cartesian(ylim = c(0, 0.35), xlim = c(-2,2)) + 
  xlab("Precipitation") +
  ylab("Colonization") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +   # gets rid of the grid lines
  theme(legend.position = "none")


VPrecip_Ext  <- ggplot(VEXT_ALL, aes(x=stdvalue)) +
  geom_line(aes(y=Ext.Precip.mean), color = "blue", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.Precip.minus1se, ymax =  Ext.Precip.plus1se  ), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  coord_cartesian(ylim = c(0, 0.35), xlim = c(-2,2)) + 
  xlab("Precipitation") +
  ylab("Extinction") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

VWNVcur_Ext  <- ggplot(VEXT_ALL, aes(x=stdvalue)) +
  geom_line(aes(y=Ext.WNVcur.mean), color = "orange", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.WNVcur.minus1se, ymax = Ext.WNVcur.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  coord_cartesian(ylim = c(0, 0.35), xlim = c(-2,2)) + 
  xlab("West Nile Virus") +
  ylab("Extinction") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")


#***VIRA turnover next best model*** need to read in the data for this first and then graph.
VIRA_precip_all <- readRDS("VIRA_precip_all.RDS")

VPrecip_turnover <- ggplot(VIRA_precip_all, aes(x=stdvalue)) +
  geom_line(aes(y=colprecip.mean), color = "blue", size=1) +
  geom_ribbon(aes(x = stdvalue, ymin = colprecip.minus1se, ymax = colprecip.plus1se), alpha = 0.3, colour = "grey", inherit.aes = TRUE) +
  geom_line(aes(y=extprecip.mean), color = "red", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = extprecip.minus1se, ymax = extprecip.plus1se), alpha = 0.3, colour = "grey", inherit.aes = TRUE) +
  xlim(-2.01,2.01) +
  ylim(0,0.35) +
  xlab("Precipitation") +
  ylab("Turnover") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")



#***BLRA Colonization Graphs***

BPrecip_col_shaded <- ggplot(BCOL_ALL, aes(x=stdvalue)) +
  geom_line(aes(y=col.precip.mean), color = "blue", size=1) +
  geom_ribbon(aes(x = stdvalue, ymin = col.precip.minus1se, ymax = col.precip.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  xlim(-2.0,2.0) +
  ylim(0,0.25) +
  xlab("Precipitation") +
  ylab("Colonization") + 
  #theme_classic () +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

BWNV_col_shaded <- ggplot(BCOL_ALL, aes(x=stdvalue)) +
  geom_line(aes(y=col.WNV.mean), color = "orange", size=1) +
  geom_ribbon(aes(x = stdvalue, ymin = col.WNV.minus1se, ymax = col.WNV.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  xlim(-2.0,2.0) +
  ylim(0,0.25) +
  xlab("West Nile Virus") +
  ylab("Colonization") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

BFreezedays_col_shaded <- ggplot(BCOL_ALL, aes(x=stdvalue)) +
  geom_line(aes(y=col.freezedays.mean), color = "purple", size=1) +
  geom_ribbon(aes(x = stdvalue, ymin = col.freezedays.minus1se, ymax = col.freezedays.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  xlim(-2.0,2.0) +
  ylim(0,0.25) +
  xlab("Freeze days") +
  ylab("Colonization") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")


#*** BLRA Extinction graphs***


# Precip*WNV  MINIMUM CONFIDENCE INTERVALS sometimes go below zero, which is why I used the coord_cartesian command (works well)
BPrecip_WNV_Ext  <- ggplot(BEXT_ALL, aes(x=stdvalue)) +
  geom_line(aes(y=Ext.Precip.mean), color = "blue", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.Precip.minus1se, ymax = Ext.Precip.plus1se), alpha = 0.3, colour = "blue", fill = "blue", inherit.aes = TRUE) +
  
  geom_line(aes(y=Ext.Precip_WNV.plus1SD.mean), color = "orange", linetype = "dashed", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.Precip_WNV.plus1SD.minus1se, ymax = Ext.Precip_WNV.plus1SD.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  
  geom_line(aes(y=Ext.Precip_WNV.minus1SD.mean), color = "orange", linetype = "dotted", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.Precip_WNV.minus1SD.minus1se, ymax = Ext.Precip_WNV.minus1SD.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE)+
  coord_cartesian(ylim = c(0, 0.50), xlim = c(-2,2)) + 
  xlab("Precipitation") +
  ylab("Extinction") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")


#FREEZEDAYS * WNV   
BFreezedays_WNV_Ext  <- ggplot(BEXT_ALL, aes(x=stdvalue)) +
  geom_line(aes(y=Ext.Freezedays.mean), color = "purple", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.Freezedays.minus1se, ymax = Ext.Freezedays.plus1se), alpha = 0.3, colour = "purple", fill="purple", inherit.aes = TRUE) +
  
  geom_line(aes(y=Ext.Freezedays_WNV.plus1SD.mean), color = "orange", linetype = "dashed", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.Freezedays_WNV.plus1SD.minus1se, ymax = Ext.Freezedays_WNV.plus1SD.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  
  geom_line(aes(y=Ext.Freezedays_WNV.minus1SD.mean), color = "orange", linetype = "dotted", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.Freezedays_WNV.minus1SD.minus1se, ymax = Ext.Freezedays_WNV.minus1SD.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  coord_cartesian(ylim = c(0, 0.50), xlim = c(-2,2)) +
  xlab("Freeze days") +
  ylab("Extinction") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")


#Precip * FREEZEDAYS 
BPrecip_Freezedays_Ext  <- ggplot(BEXT_ALL, aes(x=stdvalue)) +
  geom_line(aes(y=Ext.Precip.mean), color = "blue", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.Precip.minus1se, ymax = Ext.Precip.plus1se), alpha = 0.3, colour = "blue", fill = "blue", inherit.aes = TRUE) +
  
  geom_line(aes(y=Ext.Precip_Freezedays.plus1SD.mean), color = "purple", linetype = "dashed", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.Precip_Freezedays.plus1SD.minusSE, ymax = Ext.Precip_Freezedays.plus1SD.plusSE), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  
  geom_line(aes(y=Ext.Precip_Freezedays.minus1SD.mean), color = "purple", linetype = "dotted", size =1) +
  geom_ribbon(aes(x = stdvalue, ymin = Ext.Precip_Freezedays.minus1SD.minusSE, ymax = Ext.Precip_Freezedays.minus1SD.plusSE), alpha = 0.3, colour = "gray", inherit.aes = TRUE)+
  coord_cartesian(ylim = c(0, 0.50), xlim = c(-2,2)) + 
  xlab("Precipitation") +
  ylab("Extinction") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")


#PUTTING ALL THE GRAPHS TOGETHER
ALL_AnnPlots <- cowplot::plot_grid(BPrecip_col_shaded, BWNV_col_shaded, BFreezedays_col_shaded, BPrecip_WNV_Ext, BFreezedays_WNV_Ext, BPrecip_Freezedays_Ext, 
                                   VWNVcur_Precip_Col,VPrecip_Ext, VWNVcur_Ext, VPrecip_turnover,nrow = 4, labels = "auto",label_size = 10)   
ALL_AnnPlots




