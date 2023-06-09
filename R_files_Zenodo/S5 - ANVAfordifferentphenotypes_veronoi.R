## Date 9252019
## Personal: Jagdeep Singh Sidhu (jagdeeproots@gmail.com)
## Purpose: This code is to run ANOVA for different phenotypes by depth using Voronoi and Whole plot values. ANOVA-HSD Significant values are ploted using ggplot.
## Data needed: The data used here comes from "S4 - Calculating Voronoi Rustic and TOST_analysis_against_Wholeplot.R file" and it has RLD values as wholeplot, voronoi, for each crop 
##              which is also provided as a file named "RVW_combined_meth.csv".
## Note: Its basically to test whether voronoi and whole plot values can differentiate different phenotypes or not. Methods also include rustic, which is the unweighted average of all the coring locations.
## Figures/tables in "An analysis of soil coring strategies to estimate root depth in maize (Zea mays) and common bean (Phaseolus vulgaris)" from this code: Figure 1.


## Needed Packages 
library(plyr)
library(plotrix)
library("tripack")
library(tidyr)
library(ggplot2)
library(reshape2)
library(Rmisc) ## for summarySE
library(ggrepel) ## for avoiding overlapping labels
library(magrittr) ## for changing class of mulitple columns
library (agricolae) ## for HSD test

######################################################################################################################### 
############################ SET DIR: Set any directory you like and keep the "RVW_combined_meth.csv" file in it #########
setwd("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Comparingdifferenephenotypes")

#################################################################################
########################### Editable variables #################################
crop_type = c("bean", "maize")
methods = c("voronoi", "rustic", "Wholeplot")
#crp = "maize"
# meth = "voronoi"

#######################################################################################
###########################RLD #########################################################
## Core dimensions are needed to convert root length per slice to root length per cm3 

core_dia_cm = 4.4 # diameter of core
core_area = 3.14*(4.4)^2
slice_height_cm = 10 
core_slice_vol_cm3 = 3.14*core_dia_cm*slice_height_cm


################################################################################
#################### LOAD DATA FILE ###########################################
### The data used here comes from 7. Cal Voronoi Rustic and Wholeplot.R file and it has RLD values as wholeplot, voronoi, for each crop 


RVW_frANOVA <-read.csv("RVW_combined_meth.csv") %>% mutate_at(c("Depth", "type","method", "crop", "Rep"), funs(factor(.)) )
RVW_frANOVA$Lengthcmpercm3 <- RVW_frANOVA$Length.cm/core_slice_vol_cm3
summary(RVW_frANOVA)


#### To get summary for each type and crop combinations
summarization <- function (){ 
  for (crp in crop_type){
    for (meth in methods){
      whlplot1 <- subset(RVW_frANOVA, method == meth) ## subset method
      whlplot <- subset(whlplot1, crop == crp) ## subset crop
      summary (whlplot)
      
  if (crp == "bean"){
  bean_summary <- summarySE(whlplot, measurevar="Lengthcmpercm3", groupvars=c("type","Depth"))
  write.csv(bean_summary, paste(crp, meth, "summary_for_plotting_by_depth.csv", sep=""))} else {
    maize_summary <- summarySE(whlplot, measurevar="Lengthcmpercm3", groupvars=c("type","Depth"))
    write.csv(maize_summary, paste(crp, meth, "summary_for_plotting_by_depth.csv", sep=""))
    }
    }
  }
}
###################################################################################################
################################# ANOVA Function ##################################################
## Note: ANOVA to find out if different methods can differentiate different phenotypes of each crop
ANVAfrDiffPheno <- function(){ 
  
  for (crp in crop_type){
    for (meth in methods){
    
   
   whlplot1 <- subset(RVW_frANOVA, method == meth) ## subset method
   whlplot <- subset(whlplot1, crop == crp) ## subset crop
   summary (whlplot)
   
    ########### ANOVA ###########################################
    tD <- with(whlplot, interaction(type, Depth))
   
    aov1 = aov( Lengthcmpercm3 ~ type + Depth + tD, data = whlplot) ##
    print(crp)
    print(meth)
    print(summary(aov1))
    
    ########### tukey HSD ######################################
    HSD_tukey <- agricolae::HSD.test(aov1, "tD", group=TRUE) 
    write.csv(HSD_tukey$groups, paste(crp, meth, "HSD.csv", sep = ""))
    HSD_TK <- read.csv(paste(crp, meth, "HSD.csv", sep="")) %>% separate(X, into = paste0('comp', 1:2), sep = "[.]")
    names(HSD_TK)[1:2] <- c("type", "Depth")
    write.csv(HSD_tukey$groups, paste(crp, meth, "HSD_for_plotting.csv", sep = ""))
    
    ########## Converge HSD values with Summary data for plotting ########
    smry <- read.csv(paste(crp, meth, "summary_for_plotting_by_depth.csv", sep=""))
    smry$Depth <- as.numeric(as.character(bean_summary$Depth))
    HSD_withsummary <- merge(smry, HSD_TK, by =c("type", "Depth"))
    names(HSD_withsummary)[names(HSD_withsummary) == "Lengthcmpercm3.x"] <- "Lengthcmpercm3"
    write.csv(HSD_withsummary, paste(crp, meth, "HSD_withsummary_for_plotting.csv", sep = ""))
    
    ########## Plotting RLD for each crop and each method along with showing significant HSD values ##
    ########## This will generate the four panels of Figure 
    plot11 <- (ggplot(HSD_withsummary, aes(Lengthcmpercm3,Depth, color = type, shape = type)))
    Wholplot_by_depth_Root_Length  <- plot11 + 
      geom_errorbarh(aes(xmin=Lengthcmpercm3-se, xmax=Lengthcmpercm3+se), colour="black", height=.3)+
      geom_path(size = 0.5) +
      geom_point(size=3, fill="white")+
      geom_text(aes(label= as.factor(HSD_withsummary$groups_within_depth)),hjust=1.5, vjust= -0.5, size = 5, show.legend = FALSE)+
      
      scale_shape_manual(values=1:7)+
      ggtitle(paste(crp,"\\nComparing root length by depth among different Phenotypes using\\n", meth, sep = )) + 
      theme(axis.text = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 20), axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+ 
      ylab("Depth (cm)") + xlab("Root length density (cm/cm3)") +
      scale_y_continuous(breaks = c(-10,-20,-30,-40,-50,-60))+
      scale_x_continuous(limits = c(0, 2 ),expand = c(0,0)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), legend.key = element_blank()) #+scale_y_continuous(breaks=0:60*4)
    
    Wholplot_by_depth_Root_Length
    
    jpgname <- (paste(crp,meth, "_RL_byDepth_TUKEYHSD.svg"))
    ggsave(paste(Sys.Date(),jpgname, sep =""), plot = Wholplot_by_depth_Root_Length, height = 6, width = 8)
    
    ### I saved in a specific directory, you can change it accordingly ####
    ggsave(paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Figures/",Sys.Date(),jpgname, sep=""), 
           plot = Wholplot_by_depth_Root_Length, height = 6, width = 8)
    
    }
    }

}
##########################
## Running the functions##
summarization()
ANVAfrDiffPheno()

############################################
## What you should have after running this:
## 1. Summarized data (Mean, SD etc) for each crop and each method. 6 files in total
## 2. Files containing the ANOVA and Tukey HSD results. 6 files in total
## 3. Plots showing RLD of each crop represented by each of the 3 reference methods (Wholeplot, Voronoi, and rustic) along with Tukey HSD values: 6 files in total


 