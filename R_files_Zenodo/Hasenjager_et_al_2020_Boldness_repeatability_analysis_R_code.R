###############################################################################
#Code for repeatability analysis presented in:
#Hasenjager, M. J., et al. (2020). Proc. R. Soc. B. doi: 10.1098/rspb.2020.1871
###############################################################################

#Set working directory to where the data files are located on your computer
setwd()

#Load data
repeatData_MainStudy<-read.csv("Hasenjager et al 2020_Boldness repeatability data_Full study.csv", header=TRUE)
repeatData_PilotStudy<-read.csv("Hasenjager et al 2020_Boldness repeatability data_Pilot study.csv", header=TRUE)

#Load packages
library(rptR)
library(ggplot2)

#Obtain repeatability estimates, SEs, and bootstrapped 95% CI (note: these can take a while to run)

#Main study
rpt(BoldnessScore~(1|FishID)+(1|Cohort), grname=c("FishID","Cohort"), data=repeatData_MainStudy, datatype="Poisson", nboot=10000)

#Pilot study
rpt(BoldnessScore~(1|FishID), grname=c("FishID"), data=repeatData_PilotStudy, datatype="Poisson", nboot=10000)

#######################
#To reproduce Figure S2
#######################

#Load data
expGroupsData<-read.csv("Hasenjager et al 2020_Guppy individual-level data.csv", header=TRUE)
allIndivsTestedData<-read.csv("Hasenjager et al 2020_Mean boldness scores_All individuals assayed.csv", header=TRUE)

ggplot(dat = expGroupsData, aes(x = MeanBoldScore)) + 
    geom_histogram(fill = "red", alpha = 1) +
    geom_histogram(dat = allIndivsTestedData, fill="blue", alpha=0.4) +
    ylab("Frequency") + xlab("Mean boldness score (sec)") + 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

