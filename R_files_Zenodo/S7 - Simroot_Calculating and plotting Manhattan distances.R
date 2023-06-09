## Date: 5-20-2020
## Personal: Jagdeep Singh Sidhu (jagdeeproots@gmail.com)
## Purpose: This code:
##          1. Calculates voronoi and Rustic values based on 6 coring locations using the SoilCoreTool developed in "An analysis of soil coring strategies to estimate root depth in maize (Zea mays) and common bean (Phaseolus vulgaris)" 
##             which is available on Gitlab (https://gitlab.com/ericnord/soilCoreTools).
##          2. Then generates plots showing Manhattan distances between reference (wholeplot or voronoi) derived RLD (cm cm-3) 
##            and estimated RLD (cm cm-3) from individual and combined soil coring locations for each of the three phenotypes for each crop. 
##          
## Data needed: The data used here is comes "S2 - Simroot_fortyday&wholeplotaverage.R".

## Figures/tables in "An analysis of soil coring strategies to estimate root depth in maize (Zea mays) and common bean (Phaseolus vulgaris)" 
##              from this code: Figure 4 and 5.



## load functions

load_funcs <- function (){
source("D:/DES/Penn State/roots lab/Coring_study_Jimmy/gitlab/soilCoreTools-master/R/CheckCoreData.R")
source("D:/DES/Penn State/roots lab/Coring_study_Jimmy/gitlab/soilCoreTools-master/R/CoreLayout.R")
source("D:/DES/Penn State/roots lab/Coring_study_Jimmy/gitlab/soilCoreTools-master/R/PlotCores.R")
source("D:/DES/Penn State/roots lab/Coring_study_Jimmy/gitlab/soilCoreTools-master/R/CoreVoronoiWeights.R")
source("D:/DES/Penn State/roots lab/Coring_study_Jimmy/gitlab/soilCoreTools-master/R/VoronoiWeightedReferences.R")
source("D:/DES/Penn State/roots lab/Coring_study_Jimmy/gitlab/soilCoreTools-master/R/CompareCoreCombinations.R")
source("D:/DES/Penn State/roots lab/Coring_study_Jimmy/gitlab/soilCoreTools-master/R/sad.R")
source("D:/DES/Penn State/roots lab/Coring_study_Jimmy/gitlab/soilCoreTools-master/R/DX.R")
source("D:/DES/Penn State/roots lab/Coring_study_Jimmy/gitlab/soilCoreTools-master/R/CoreVisualizationTools.R")
#install.packages("tripack")
#install.packages("tidyverse")
#install.packages("plyr")
#install.packages("plotrix")
#install.packages("Rmisc")
library(plyr)
library(plotrix)
library("tripack")
library(tidyverse)
library(ggplot2)
library(reshape2)
library(Rmisc) ## for summarySE
}

#################Editable variables - list ##################################
Root_type <- c("maize_shallow_", "maize_", "maize_deep_", "bean_shallow_", "bean_deep_", "bean_")
methods <- c("voronoi", "Wholeplot")

###########################RLD ############################
######################Core dimensions#######################
############################################################
core_dia_cm = 4.4
core_area = 3.14*(4.4)^2
slice_height_cm = 10
core_slice_vol_cm3 = 3.14*core_dia_cm*slice_height_cm

################ Set dir and read file ###########
setdir <- function(){
  directory <- paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/",type, sep ="")
  
  #Set Working Directory
  setwd(paste(directory))
  getwd()
  message(paste("/working directory - ",directory))

  }

main<- function(){
  
  ###################################################################
  ##############################read file########################
  #file = read.csv(paste(type,"for_corecomparison.csv",sep=""))
  file<-(checkCoreData(file=paste(type,"for_corecomparison.csv",sep="")))
  file$RLD <- file$Length.cm/core_slice_vol_cm3
  
  
  
  
  # to incorporate intermediate into the plots of intermediate angle #
  if (type =="bean_"|type == "maize_"){plot_name = paste(type,"intermediate",sep = "")
  }else{plot_name = paste(type)
  }
  #plot_name <<- plot_name #to export to another function
###############################################################################################################
##################### 2. Core layout around the focal plant ###################################################
layout.simroot <- coreLayout(px=76.2,py=23,cores=1:6,x.locs="c(0.25*py,0,0.5*py,0.5*px,0.5*px,0)",
                           y.locs="c(0,0.5*py,0,0,0.5*py,0.25*py)",ids=c("L1","L2","L3","L4","L5","L6"),
                           type=c("e","c","e","c","c","e"))
layout.maize = layout.simroot

plotCores(layout.simroot,size="whole")


################################################################################################################
###################### 3. Calculate area-based weights for core locations ######################################
# all core locations
Simroot.weights<-coreVoronoiWeights(layout.simroot,plot.cores="whole")

# only some core locations
#coreVoronoiWeights(layout.beanbean,plot.cores="whole",which=c(1,3,4,6)) # only some

##################################################################################################################
##### 4.  Calculate reference values #############################################################################
voronoi<-voronoiWeightedReferenceValues(core.weights=Simroot.weights, core.data=file, data.cols=7)
colnames(voronoi)[2] = "RLD"
voronoi

####### calculate rustic averages ###
rustic = as.data.frame(tapply(file$RLD,file$Depth,mean))
colnames(rustic) = c("RLD")
rustic
Depth = c(-60,-50,-40,-30,-20,-10)
rustic = cbind(Depth,rustic)
rustic

########## Using whole plot average #######
wholeplot <- read.csv(paste(type,"wholeplot_summary.csv"))
wholeplot$RLD <- wholeplot$Length.cm/core_slice_vol_cm3
Wholeplot1 = wholeplot %>% select(Depth, "RLD")

####bind rustic, voronoi, and whole plot average by depth
rustic_and_voronoi = merge.data.frame(rustic,voronoi, by =  "Depth")
rus_vor_whole = merge.data.frame(rustic_and_voronoi,Wholeplot1,by = "Depth")
colnames(rus_vor_whole) = c("Depth", "Rustic", "Voronoi","Wholeplot")
rus_vor_whole
rus_vor_whole2 = melt(rus_vor_whole, id.vars = "Depth")
colnames(rus_vor_whole2) = c("Depth", "Method", "RL")


# plot1 <- (ggplot(rus_vor_whole2, aes(RL,Depth, color = Method)))
# Compare_ref_methods <- plot1 + geom_point(size = 3) + 
#   geom_path() + ggtitle(paste(plot_name, "reference comparison")) + theme(plot.title = element_text(face ="italic"))+ ylab("Depth (cm)") + xlab("Root Segment Length (cm)") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Compare_ref_methods
# jpgname = paste("Compare_ref_methods_",plot_name,".jpg")
# ggsave(jpgname, plot = Compare_ref_methods, height = 6, width = 8)
# ggsave(paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Figures/",jpgname, sep=""), plot = Compare_ref_methods, height = 6, width = 8)
# 

################ write csv files for refs to be read by corecompare function #########
write.csv(voronoi,paste("voronoi",plot_name,".csv",sep = ""))
write.csv(rustic,paste("rustic",plot_name,".csv",sep = ""))
write.csv(Wholeplot1,paste("Wholeplot",plot_name,".csv",sep = ""))



for (method in methods){
  refs = read.csv(paste(method,plot_name,".csv",sep = ""))
  #file = subset(file, Rep == 1) only for testing if the code is working fine 
##################################################################################################################
##### 5. Compare core combinations against the reference values ##################################################
Simroot.results<-compareCoreCombinations(x=file,layout = layout.maize, refs = refs,var = "RLD", verbose=TRUE, save.wts=TRUE)
warnings()

Simroot.results<-Simroot.results[order(Simroot.results$n,Simroot.results$SAD.m),]
 Simroot.results[1:10,1:11] # only showing means to not clutter the document

#The function `sad()` calculates the sum of the absolute differences between 2 vectors.
#The function `DX()` calculates depth to a given percentile, such as D95 or D75.


##### 6. Plotting 

###plotting by locs 
Simroot.results$n <- as.numeric(Simroot.results$n)

#To get top three combinations for n
oneloc <- filter(Simroot.results, n == '1')
morelocs <- filter(Simroot.results, n > 1)

Top3_locs_echcomb_morelocs = as.data.frame (morelocs %>% group_by(n) %>% top_n(-3, SAD.m))
Top3_locs_echcomb_1 = rbind(Top3_locs_echcomb_morelocs, oneloc)
summary(Top3_locs_echcomb_1)

## order the x axis
Top3_locs_echcomb_2<-Top3_locs_echcomb_1[order(Top3_locs_echcomb_1$n,Top3_locs_echcomb_1$SAD.m),]
Top3_locs_echcomb_2$locs<- factor(Top3_locs_echcomb_2$locs, levels=unique(Top3_locs_echcomb_2$locs))
Top3_locs_echcomb_2$n = as.factor(Top3_locs_echcomb_2$n)
Top3_locs_echcomb = Top3_locs_echcomb_2 %>% rename("No._of_locs" = n)

##########################################################################################################
### by depth ---one loc ###################################################################################
##########################################################################################################
### by depth ---one loc ###
# one_loc_bydepth1 = (oneloc[,c(1,6:11)] %>% 
#                       rename("-10" = 'Dif.0.m', "-20" = 'Dif.1.m',"-30" = 'Dif.2.m', "-40" = 'Dif.3.m',"-50" = 'Dif.4.m', "-60" = 'Dif.5.m'))
# one_loc_bydepth =melt(data = one_loc_bydepth1, id.vars = c("locs" ), value.name = "Length.cm", variable.name = "Depth")
# one_loc_bydepth$Depth = as.numeric(as.character(one_loc_bydepth$Depth))
# 
# 
# plot11 <- (ggplot(one_loc_bydepth, aes(Length.cm,Depth, color = locs)))
# loc_by_depth_Prediction_Precision_bean <- plot11 + geom_point(size = 3) + 
#   geom_path() + ggtitle(paste(type,"loc_by_depth_Prediction_Precision")) + theme(plot.title = element_text(face ="italic"))+ ylab("Depth (cm)") + xlab(paste("[",type, "- Simulated coring] Root Segment length (cm)")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# loc_by_depth_Prediction_Precision_bean
# jpgname <- paste(Sys.Date(),plot_name,"loc_by_depth_Prediction_Precision.jpg")
# ggsave(jpgname, plot = loc_by_depth_Prediction_Precision_bean, height = 6, width = 8)
# ggsave(paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Figures/",jpgname, sep=""), plot = loc_by_depth_Prediction_Precision_bean, height = 6, width = 8)


#######################################################################################################################
############################ Ploting against reference values ###########################################################
##

Simroot_against_ref = ggplot(Top3_locs_echcomb, aes(x=locs, y=SAD.m, ymin=SAD.m - SAD.sd, ymax=SAD.m + SAD.sd))+
  geom_pointrange((aes(ymin=SAD.m - SAD.sd, ymax=SAD.m + SAD.sd))) + 
  geom_linerange(aes(color = No._of_locs), size = 1, alpha = 0.5)+
  geom_point(aes(color = No._of_locs, shape = No._of_locs), size = 2)+
  xlab('Location or combination of locations')+ 
  ylab("Sum of Absolute differences \\n in RLD (cm/cm3)\\n") +
  scale_y_continuous(limits = c(0,15),expand = c(0,0)) +
  #ggtitle(paste(plot_name, " Simroot data vs ", method, " estimated whole plot average\\n Best Combined locations", sep = ""))+ theme_grey(base_size = 16) +
  theme(axis.text = element_text(colour = "black", size = 14,angle = 90, hjust = 1, vjust=0.2), axis.title = element_text(colour = "black", size = 20), axis.line = element_line(colour = "black", size = 0.2, linetype = "solid"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.key = element_blank(), legend.position = c(0.9, 0.5)) #+scale_y_continuous(breaks=0:60*4)


jpgname1 <- paste(Sys.Date(),plot_name,"Simroot_ vs",method, ".eps")
ggsave(jpgname1, plot = Simroot_against_ref, width = 12)
ggsave(paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Figures/",jpgname1, sep=""), plot = Simroot_against_ref, width = 12)

# ## plotting by n
# Simroot.results$n <- as.factor(Simroot.results$n)
# Simroot_MPC =ggplot(Simroot.results, aes(x=n, y= SAD.m)) + 
#   geom_boxplot() +
#   xlab('Number of locations pooled')+ 
#   ylab("Sum of Absolute differences [Root Length (cm)]\\n") + 
#   ggtitle(paste(plot_name, " Simroot data vs ", method, " \\nestimated whole plot average", sep = "")) + theme_grey(base_size = 16) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5))
# jpgname2 <- paste(Sys.Date(),plot_name,"Simroot_ vs",method, "MPC.jpg")
# ggsave(jpgname2, plot = Simroot_MPC, width = 8)
# ggsave(paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Figures/",jpgname2, sep=""), plot = Simroot_MPC, width = 8)
 }

}

######## Call function ##########
load_funcs()
for (type in Root_type) {
  setdir()
  main()
  }
    warnings()

    
#############################################################################################
## After running this, for each crop and year you should have:
## 1. Plot showing showing Manhattan distances between voronoi and wholeplot derived RLD (cm cm-3) from individual and combined soil coring locations

