## Date: 9-25-2019
## Personal: Jagdeep Singh Sidhu (jagdeeproots@gmail.com)
## Purpose: This code:
##          1. Calculates voronoi and Rustic values based on 6 coring locations using the SoilCoreTool developed in "An analysis of soil coring strategies to estimate root depth in maize (Zea mays) and common bean (Phaseolus vulgaris)" 
##             which is available on Gitlab (https://gitlab.com/ericnord/soilCoreTools).
##          2. Voronoi and rustic averages are then compared with wholeplot average (by depth) for each phenotype (Deep, intermediate, and shallow for bean and maize)
##             using TOST equivalence test
##          3. Comparisons are then shown on a RLD plot, where "*" indicate significance equivalence of a method to wholeplot at a particular depth
##
## Data needed: The data used here is generated using "S2. Simroot_fortyday&wholeplotaverage.R file". So make sure you run that.

## Note: Its basically to test whether voronoi estimated RLD values are equivalent to wholeplot average. Methods also include rustic, which is the unweighted average of all the coring locations.
## Figures/tables in "An analysis of soil coring strategies to estimate root depth in maize (Zea mays) and common bean (Phaseolus vulgaris)" from this code: Figure 3.

#######################################################
## Needed packages can be loaded using load functions##


load_funcs <- function (){
  ## Source files below can be found on Gitlab (https://gitlab.com/ericnord/soilCoreTools)
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
  #install.packages("dplyr")
  #install.packages("plotrix")
  #install.packages("Rmisc")
  #install.packages("equivalence)
  library(dplyr)
  library(plotrix)
  library("tripack")
  library("tidyverse")
  library(ggplot2)
  library(reshape2)
  library(Rmisc) ## for summarySE
  library(ggrepel) ## for avoiding overlapping labels
  library(equivalence) ## for TOST
}

#############################################################################
#################Editable variables - list ##################################
Root_type <- c("maize_shallow_", "maize_", "maize_deep_", "bean_shallow_", "bean_deep_", "bean_")

##########################################################
## One time initiation of the combined results dataframe
TOST_Combined_meth = data.frame("Depth", "Location", "comparedto", "P-value", "decision","Type", stringsAsFactors=FALSE) ### dataframe for combined results
RVW_combined <- data.frame("Depth", "Length.cm", "Rep", "method","type", "crop", "Lengthcmpercm3", stringsAsFactors=FALSE)
names(RVW_combined) <- c("Depth", "Length.cm", "Rep", "method","type", "crop","Lengthcmpercm3" )
RVW_combined <- RVW_combined[-1,]

#######################################################################################
###########################RLD #########################################################
## Core dimensions are needed to convert root length per slice to root length per cm3 
core_dia_cm = 4.4
core_area = 3.14*(4.4)^2
slice_height_cm = 10
core_slice_vol_cm3 = 3.14*core_dia_cm*slice_height_cm


############################# 1. SET directory function #######################################
setdir <- function(){
  directory <- paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/",type, sep ="")
  
  #Set Working Directory
  setwd(paste(directory))
  getwd()
  message(paste("/working directory - ",directory))
  file_name = read.csv(paste(type,"for_corecomparison.csv",sep="")) ### read file for voronoi and rustic averages
  file_name <<- file_name
}

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
##### 4.  Calculate Voronoi values #############################################################################

vor_by_rep <- function(){
  voronoi_by_rep = data.frame("Depth", "Length.cm", "Rep", "type", stringsAsFactors=FALSE) ### dataframe for saving results
  colnames(voronoi_by_rep) <- c("Depth", "Length.cm", "Rep", "type")
  voronoi_by_rep <- voronoi_by_rep[-1,]
  i = 1 ## variable for increasing the row number in results dataframe
  r =1
  for (r in 1:100) {
    voronoi_temp <- voronoiWeightedReferenceValues(core.weights=Simroot.weights, core.data=subset(file_name,Rep ==r), data.cols=6)
    
    rep_type <- data.frame("Rep" = r, "type" = type) # creating a data frame for rep and type
    rep_type1 <- rep_type %>% slice(rep(1:n(), each = 6))
    
    voronoi_temp <- cbind(voronoi_temp, rep_type1)
    names(voronoi_temp)[2] <- "Length.cm"
    voronoi_by_rep <- rbind(voronoi_temp, voronoi_by_rep)
    
    
    i = i+6
    r = r+1
    voronoi_by_rep <<- voronoi_by_rep
  }
  voronoi_by_rep$method <- "voronoi"
  write.csv(voronoi_by_rep, paste(type, "voronoibyrep.csv", sep =""))
  rep_type1 <<- rep_type1
}

#################################################################################################################
####### 5. calculate rustic averages ###############################################################################
RUSTIC <- function(){
  rustic1 <- as.data.frame (aggregate(Length.cm ~ Depth + Rep, data = file_name, FUN= "mean" ))
  rustic <- cbind(rustic1, rep_type1$type)
  names(rustic)[4]<- "type"
  rustic$method <- "rustic"
  write.csv(rustic, paste(type, "rusticbyrep.csv", sep =""))
}

for (type in Root_type){
  load_funcs()
  setdir()
  vor_by_rep()
  RUSTIC()
} 

####################################################################################################################
####### 6. Get wholeplot average and merge with voronoi and rustic #################################################
RVW_func <- function() {
  whlplt <- read.csv(paste(type, "wholeplotfraov.csv")) ##### load wholeplot for each rep - 100 reps for siumulations
voronoi = read.csv(paste(type,"voronoibyrep",".csv",sep = "")) 
rustic <- read.csv(paste(type,"rusticbyrep",".csv",sep = ""))
summary(whlplt)
summary(voronoi)
summary(rustic)

rus_vor <- rbind(rustic, voronoi)
rus_vor1 <- rus_vor[,-c(1)] # remove x column

whlplt1 <- whlplt[,-c(2)] # remove time column
whlplt1$type <- type
names(whlplt1)[4] <- "method"
summary(whlplt1)
summary(rus_vor1)

RVW <- rbind(whlplt1, rus_vor1)

if (type == "maize_shallow_"| type == "maize_" | type == "maize_deep_"){
  RVW$crop <- "maize"
}else {
  RVW$crop <- "bean"
}

RVW$Lengthcmpercm3 <- RVW$Length.cm/core_slice_vol_cm3 ###### create column for RLD per cm3
RVW <<-  RVW
RVW_combined <- rbind (RVW, RVW_combined)
RVW_combined <<- RVW_combined
}


###################################################################################################
################ for combining all phenotypes #####################################################
RVWcomb <- function(){
RVW_combined <- rbind (RVW, RVW_combined)
}

###################################################################################################
####### 7. TOST for finding locations by depth which are close to wholeplot average
TOST_Fnc_meth <- function(){
  
  ##### creating variables for loop ###
  method = c("rustic", "voronoi")
  Horizon = c("-10", "-20", "-30", "-40", "-50", "-60")
  TOST_results_meth = data.frame("Depth", "Location", "comparedto", "P-value", "decision","Type", stringsAsFactors=FALSE) ### dataframe for saving results
  i = 1 ## variable for increasing the row number in results dataframe
  
  for (X in Horizon){
    for (meth in method){
      
      Y =   subset(RVW, Depth == X & method =="Wholeplot") ### Wholeplot at X depth
      Z = subset(RVW, Depth == X &  method == meth) ### RLD at X depth and Location - Loc
      cut = (3/10*mean(Y$Lengthcmpercm3)) #defining cutoff
      results = tost(Y$Lengthcmpercm3, Z$Lengthcmpercm3, epsilon = cut, paired = FALSE, var.equal = FALSE, conf.level = 0.95, alpha = 0.05) ## epsilon is the user given difference.
      TOST_results_meth[i,] = c(X, meth, "Wholeplot", results$tost.p.value, results$result, type)
      i = i+1
      TOST_results_meth <<- TOST_results_meth
    }
  }
  
  #write csv of the results
  write.csv(TOST_results_meth, paste(type, "TOST_results_meth.csv", sep =""))
  TOST_Combined_meth = rbind(TOST_Combined_meth, TOST_results_meth)
  TOST_Combined_meth <<- TOST_Combined_meth
  
}

################################################################################################################
############################# generate data for plotting #######################################################
data_fr_plotting<- function (){
  
RVW_summary <- summarySE(RVW, measurevar="Lengthcmpercm3", groupvars=c("method","Depth"))
names(TOST_results_meth) <- c("Depth", "method", "compareto", "p-value", "TOST_decision", "type")

# Merge the outcome values from TOST_results into RVW_summary
test1_meth = merge(RVW_summary, TOST_results_meth, by=c("Depth","method"), all.x=TRUE)
test1_meth$`p-value` <- as.numeric(test1_meth$`p-value`)

test_meth = cbind(test1_meth, ifelse (test1_meth$`p-value`< 0.05, "*", "")) ## * for significane
colnames(test_meth)[12] <- "signi"
test_meth <<- test_meth
}


############################# plot RLD by depth for each loc and wholeplot ###
Plot_ref_comp<- function(){  
  plot11 <- (ggplot(test_meth, aes(Lengthcmpercm3,Depth, color = method, shape = method)))
  Core_loc_by_depth_Root_Length  <- plot11 + 
    geom_errorbarh(aes(xmin=Lengthcmpercm3-se, xmax=Lengthcmpercm3+se), colour="black", height=.3)+
    geom_path(size = 0.5) +
    geom_point(size=3, fill="white")+
    geom_text(aes(label= as.factor(test_meth$signi)),hjust=0, vjust= 0, size = 8, show.legend = FALSE)+
    scale_shape_manual(values=1:7)+
    ggtitle(paste(type,"\\nRoot length by depth in different reference methods\\n*TOST-signi equivalent to wholeplot average")) + 
    theme(axis.text = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 20), axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+ 
    ylab("Depth (cm)") + xlab("Root length density (cm/cm3)") +
    scale_y_continuous(breaks = c(-10,-20,-30,-40,-50,-60))+
    scale_x_continuous(limits = c(0, 2.5 ),expand = c(0,0)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), legend.key = element_blank()) #+scale_y_continuous(breaks=0:60*4)
  
  Core_loc_by_depth_Root_Length
  
  
  
  jpgname <- (paste(Sys.Date(),type, "RL_by_depth_reference_methods&Wholeplot_TOST.svg"))
  ggsave(jpgname, plot = Core_loc_by_depth_Root_Length, height = 6, width = 8)
  ggsave(paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Figures/",jpgname, sep=""), plot = Core_loc_by_depth_Root_Length, height = 6, width = 8)
  
}




for (type in Root_type){
  setdir()
  RVW_func()
  RVWcomb()
  TOST_Fnc_meth()
  data_fr_plotting()
  Plot_ref_comp()
} 

write.csv(TOST_Combined_meth, paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/TOST/TOST_results_meth.csv", sep =""))
write.csv(RVW_combined, paste("D:/DES/Penn State/roots lab/Coring_study_Jimmy/Simroot/Burridge_coring_20190508/Comparingdifferenephenotypes/RVW_combined_meth.csv", sep =""))


#############################################################################################
## After running this code, within in each phenotype (six phenotypes) folder you should have:
## 1. 1 file (type,"voronoibyrep.csv") containing Voronoi estimated RLD values by rep. 
## 2. 1 file (type, "rusticbyrep.csv") containing rustic estimated RLD values by rep.
## 3. 1 file (type, "TOST_results_meth.csv") containing TOST results.
## 4. 1 plot showing comparison among reference methods along with TOST significance values (6 panels in Figure 3.)