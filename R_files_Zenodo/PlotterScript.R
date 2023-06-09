# In the Name of God
# This script allows the user to plot a multiplot entity from multiple aggregated assessment data structure 


# install.packages("ggpubr")
# setwd("~/GitHub/utlp-spm-mux")
library(ggplot2) 
library(ggpubr)
library(igraph)
library(LinkPrediction)
library(lattice)

#### Loading Saved Files ####

assessment_11 = 
  readRDS(file = "AggregatedAssesment_Air-Train-Train-k=1-paper.rds")

assessment_12 = 
  readRDS(file = "AggregatedAssesment_Air-Train-Train-k=1-paper.rds")

assessment_13 = 
  readRDS(file = "AggregatedAssesment_Air-Train-Train-k=1-paper.rds")

assessment_21 = 
  readRDS(file = "AggregatedAssesment_Air-Train-Train-k=69-paper.rds")

assessment_22 = 
  readRDS(file = "AggregatedAssesment_Air-Train-Train-k=69-paper.rds")

assessment_23 = 
  readRDS(file = "AggregatedAssesment_Air-Train-Train-k=69-paper.rds")



#### Creating Plots Using Link prediction Function ####
# Metric is determined here
P11 = plot_lp_performance(assessment_11, metric = "auroc", err = "sd")
P12 = plot_lp_performance(assessment_12, metric = "recall_at_k", err = "sd")
P13 = plot_lp_performance(assessment_13, metric = "avg_prec", err = "sd")
P21 = plot_lp_performance(assessment_21, metric = "auroc", err = "sd")
P22 = plot_lp_performance(assessment_22, metric = "recall_at_k", err = "sd")
P23 = plot_lp_performance(assessment_23, metric = "avg_prec", err = "sd")

# Changing labels here is possible
# y changing y labels 
# To add k=1 \\n is used to go to the next line
P11$labels$y = "K=1 \\n AUROC"
P21$labels$y = "K=All \\n AUROC"
P12$labels$y = "Precision"
P22$labels$y = "Precision"

# omiting x label
# P11$labels$x = ""
# P12$labels$x = ""
# P13$labels$x = ""
# P21$labels$x = ""
# P22$labels$x = ""
# P23$labels$x = ""

# Adding title to the upper row
# theme(plot.title = element_text(hjust = 0.5) places the title in the middel
P11 = P11 + ggtitle("Train") +  theme(plot.title = element_text(hjust = 0.5))
P12 = P12 + ggtitle("Train") +  theme(plot.title = element_text(hjust = 0.5))
P13 = P13 + ggtitle("Train") +  theme(plot.title = element_text(hjust = 0.5))

# deactivating background grid 
P11$theme$panel.grid = element_blank()
P12$theme$panel.grid = element_blank()
P13$theme$panel.grid = element_blank()
P21$theme$panel.grid = element_blank()
P22$theme$panel.grid = element_blank()
P23$theme$panel.grid = element_blank()

# changing ticks on x axix. 
# seq(0.0, 1, 0.2) means start from 0 to 1 by 0.2 steps.
# limits defines start and end.
P11 = P11  + scale_x_continuous(breaks=seq(0.0, 1, 0.1), limits=c(0.1,0.9))
P12 = P12  + scale_x_continuous(breaks=seq(0.0, 1, 0.1), limits=c(0.1,0.9))
P13 = P13  + scale_x_continuous(breaks=seq(0.0, 1, 0.1), limits=c(0.1,0.9))
P21 = P21  + scale_x_continuous(breaks=seq(0.0, 1, 0.1), limits=c(0.1,0.9))
P22 = P22  + scale_x_continuous(breaks=seq(0.0, 1, 0.1), limits=c(0.1,0.9))
P23 = P23  + scale_x_continuous(breaks=seq(0.0, 1, 0.1), limits=c(0.1,0.9))

# changing ticks on y axix. 
# seq(0.0, 1, 0.2) means start from 0 to 1 by 0.2 steps.
# limits defines start and end.
# to change start and end, both breaks=seq(0,1,by=0.2), limits=c(0,1) should be changed
# for example breaks=seq(0.5,1,by=0.1), limits=c(0.5,1) changes limits from 0.5 to 1 by 0.1 steps 
P11 = 
  P11 + 
  scale_y_continuous(breaks=seq(0.5, 1, by=0.1), 
                     limits=c(0.5, 1)) 
P12 = 
  P12 + 
  scale_y_continuous(breaks=seq(0, 0.5, by=0.1), 
                     limits=c(0, 0.5)) 
P13 = 
  P13 + 
  scale_y_continuous(breaks=seq(0, 0.3, by=0.1), 
                     limits=c(0, 0.3)) 
P21 = 
  P21 + 
  scale_y_continuous(breaks=seq(0.5, 1, by=0.1), 
                     limits=c(0.5, 1)) 
P22 = 
  P22 + 
  scale_y_continuous(breaks=seq(0, 0.5, by=0.1), 
                     limits=c(0, 0.5)) 
P23 = 
  P23 + 
  scale_y_continuous(breaks=seq(0, 0.3, by=0.1), 
                     limits=c(0, 0.3)) 

# emits the legend title which used to be "predictor"
P11$theme$legend.title = element_blank()
P12$theme$legend.title = element_blank()
P13$theme$legend.title = element_blank()
P21$theme$legend.title = element_blank()
P22$theme$legend.title = element_blank()
P23$theme$legend.title = element_blank()

# It legend are to be omitted, the following lines should be uncommented
# P11$theme$legend.position = 'none'
# P12$theme$legend.position = 'none'
# P13$theme$legend.position = 'none'
# P21$theme$legend.position = 'none'
# P22$theme$legend.position = 'none'
# P23$theme$legend.position = 'none'


# drawing the multiline 
# ncol determines the number of columns
# nrow determines the number of rows
ggarrange(P11, P12, P13, P21, P22, P23, 
          ncol = 3,
          nrow = 2)




