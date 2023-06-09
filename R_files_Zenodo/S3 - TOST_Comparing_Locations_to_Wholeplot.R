## Date: 9-11-2019
## Personal: Jagdeep Singh Sidhu (jagdeeproots@gmail.com)
## Purpose: This code:
##          1. Runs TOST (Equivalence Test) to test (at a particular depth) if a coring location is equivalent to wholeplot average or not. 
##             Significant values over RLD for each equivalence are plotted using ggplot. The comparison presented here in the plots is only against whole plot average.
## Data needed: The data used here comes from "S2 - Simroot_fortyday&wholeplotaverage"/
## Figures/tables in "An analysis of soil coring strategies to estimate root depth in maize (Zea mays) and common bean (Phaseolus vulgaris)" 
##                from this code: Figure 8.


## Needed packages 
#install.packages("tripack")
#install.packages("tidyverse")
#install.packages("plyr")
#install.packages("plotrix")
#install.packages("Rmisc")
#install.packages("equivalence)
library(plyr)
library(plotrix)
library("tripack")
library(tidyverse)
library(ggplot2)
library(reshape2)
library(Rmisc) ## for summarySE
library(ggrepel) ## for avoiding overlapping labels
library(equivalence) ## for TOST

#######################################################################################
########################### Editable variables - List #################################
Root_type <- c("maize_shallow_", "maize_", "maize_deep_", "bean_shallow_", "bean_deep_", "bean_")

#one time Initiation of the combined results dataframe
TOST_Combined = data.frame("Depth", "Location", "comparedto", "P-value", "decision","Type", stringsAsFactors=FALSE) ### dataframe for combined results

#######################################################################################
###########################RLD #########################################################
## Core dimensions are needed to convert root length per slice to root length per cm3 
core_dia_cm = 4.4
core_area = 3.14*(4.4)^2
slice_height_cm = 10
core_slice_vol_cm3 = 3.14*core_dia_cm*slice_height_cm

###############################################################################################
############################# 1. SET directory function #######################################
setdir <- function(){
  directory <- paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/",type, sep ="")
  
  #Set Working Directory
  setwd(paste(directory))
  getwd()
  message(paste("/working directory - ",directory))
  
}

#############################################################################
############## 2. Load whole plot and locations data ########################
loaddata <- function(){
  
  whlplt <- read.csv(paste(type, "wholeplotfraov.csv")) ##### load wholeplot for each rep - 100 reps for simulations - comes from "S2 - Simroot_fortyday&wholeplotaverage.R"
  whlplt$Lengthcmpercm3 <- whlplt$Length.cm/core_slice_vol_cm3
  summary(whlplt)
  
  locswise <- read.csv(paste(type,"for_corecomparison.csv",sep="")) #### read loc wise data for all the locs  ## comes from "S2 - Simroot_fortyday&wholeplotaverage"
  locswise$Location <- as.factor(locswise$Location) ### for merging loc wise and wholeplot 
  whlpandlocs <- rbind(whlplt,locswise[,-1])
  summary(whlpandlocs)
  whlpandlocs$Location = as.factor(whlpandlocs$Location) ## change to factors 
  whlpandlocs$Depth = as.factor(whlpandlocs$Depth) ## Change depth to factor
  
  whlpandlocs <<- whlpandlocs
  
  ####### read sumarized wholeplot average ##
  wholeplot <- read.csv(paste(type,"locandwholeplot.csv")) ## comes from "S2 - Simroot_fortyday&wholeplotaverage"
  wholeplot <<- wholeplot
  
  locAndWhlFrplot <- dplyr::arrange(wholeplot, desc(Depth))
  locAndWhlFrplot <<- locAndWhlFrplot
}

##############################################################################################
####### TOST for finding locations by depth which are close to whole plot average ############
TOST_Fnc <- function(){

  
whlpandlocs$Location~as.factor(whlpandlocs$Location)

##### creating variables for loop ###
Location = c("1","2", "3", "4", "5", "6")
Horizon = c("-10", "-20", "-30", "-40", "-50", "-60")
TOST_results = data.frame("Depth", "Location", "comparedto", "P-value", "decision","Type", stringsAsFactors=FALSE) ### dataframe for saving results
i = 1 ## variable for increasing the row number in results dataframe

for (X in Horizon){
for (Loc in Location){

  Y =   subset(whlpandlocs, Depth == X & Location =="Wholeplot") ### Wholeplot at X depth
  Z = subset(whlpandlocs, Depth == X & Location ==Loc) ### RLD at X depth and Location - Loc
results = tost(Y$Lengthcmpercm3, Z$Lengthcmpercm3, epsilon = 20, paired = FALSE, var.equal = FALSE, conf.level = 0.95, alpha = 0.05) ## epsilon is the user given difference.
TOST_results[i,] = c(X, Loc, "Wholeplot", results$tost.p.value, results$result, type)
i = i+1
TOST_results <<- TOST_results
}
  }

#write csv of the results
write.csv(TOST_results, paste(type, "TOST_results.csv", sep =""))
TOST_Combined = rbind(TOST_Combined, TOST_results)
TOST_Combined <<- TOST_Combined

#### data for plot
names(TOST_results) <- c("Depth", "Location", "compareto", "p-value", "TOST_decision", "type")

# Merge the outcome values from TOST_results into locandwhlFrplot
test1 = merge(locAndWhlFrplot, TOST_results, by=c("Depth","Location"), all.x=TRUE)
test = cbind(test1, ifelse (test1$`p-value`< 0.05, "*", "")) ## * for significane
colnames(test)[13] <- "signi"
test <<- test

}


#########################################################################################################
############################# plot RLD by depth for each loc and wholeplot -TOST ########################

Plot_whole_core<- function(){  
  plot11 <- (ggplot(test, aes(Lengthcmpercm3,Depth, color = Location, shape = Location)))
  Core_loc_by_depth_Root_Length  <- plot11 + 
    geom_errorbarh(aes(xmin=Lengthcmpercm3-se, xmax=Lengthcmpercm3+se), colour="black", height=.3)+
    geom_path(size = 1) +
    geom_point(size=3, fill="white")+
    geom_text(aes(label= as.factor(test$signi)),hjust=0, vjust= 0, size = 8)+
    scale_shape_manual(values=1:7)+
    ggtitle(paste(type,"\\nRoot length by depth in 6 coring locations and the whole plot average\\n*TOST-signi from wholeplot average")) + 
    theme(plot.title = element_text(face ="italic", hjust = 0.5))+ 
    ylab("Depth (cm)") + xlab("Root Segment length (cm)") +
    scale_y_continuous(breaks = c(-10,-20,-30,-40,-50,-60))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) #+scale_y_continuous(breaks=0:60*4)
  
  Core_loc_by_depth_Root_Length
  
  
  
  jpgname <- (paste(type, "RL_by_depth_6locs&Wholeplot_TOST.jpg"))
  ggsave(jpgname, plot = Core_loc_by_depth_Root_Length, height = 6, width = 8)
  ggsave(paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Figures/",jpgname, sep=""), plot = Core_loc_by_depth_Root_Length, height = 6, width = 8)
  
}

for (type in Root_type){
  setdir()
  vor_by_rep()
  RUSTIC()
  #loaddata()
  #TOST_Fnc()
  #Plot_whole_core()
} 
write.csv(TOST_Combined, paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/TOST/TOST_results.csv", sep =""))

#############################################################################################
## After running this code, within in each phenotype (six phenotypes) folder you should have:
## 1. 2 files (type, "TOST_results_meth.csv") containing TOST results. 
## 2. 1 plot showing comparison among locations and wholeplot average along with TOST significance values (6 panels in Figure 8.)




































