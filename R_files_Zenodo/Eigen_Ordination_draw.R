#GC plots for lab meeting doc--replicates FPCA_Plotting_Stacy ggplot code_Shuxi 
#edits (which was for DSI plots)
setwd("~/Dropbox/Things we generated during Shuxi's visit/R code/Eigen_function_Ordination_draw")
library(plotrix)
library(gridExtra)
rm(list=ls())
load("plot_DSI_F_pc.RData")
###one run of plotting
################################
##0,1 Comparison
###Extract the residuals by minusing the mean value
Residuals=unlist(residuals)
process_mean=process_value_fitted
##create a list for baboon names
##Complete.data is the cleaned data
#name.list=unique(Complete.data$sname)
#Random Pick Up an Individual
random.name=sample(name.list,1)

library(dplyr)
library(ggplot2)
library(ggforce)

#Principal components summary

#plot size is 720x720 using export
#pdf("DSI_F_Eigen.pdf")
plot(age_grids,eigen_1,type='l',col="black",ylim=c(-1.5,2),
     ylab="Eigen Function",xlab="Age", cex.axis=1.5, cex.lab=1.5)
lines(age_grids,eigen_1_up,lty=2)
lines(age_grids,eigen_1_down,lty=2)
lines(age_grids,eigen_2,col="red")
lines(age_grids,eigen_2_up,lty=2)
lines(age_grids,eigen_2_down,lty=2)
legend("topright",legend=c(paste("1st PC",format(var_explain[1]*100,digits=4),"%"),paste("2nd PC",format(var_explain[2]*100,digits=4),"%")),lty=1,col=c("black","red"), 
       cex=1.2, bty = 'n')
#dev.off()

#Find same animals on graph that I did for DSI. Note that the quadrants are not the same, but this way
#you can follow the same animals through the entire section

#Need to put sname and adversities data into process data

pca_scores_gcs <- data.frame(score_1, score_2, Treatment_Ind,sname=unique(Complete.data$sname))
snames_cumadversities <- distinct(Complete.data, sname, adv_cumulative)
pca_scores_gcs_cumadversities <- data.frame(pca_scores_gcs, snames_cumadversities)
pca_scores_gcs_cumadversities$adv_cum_asfactor <- as.factor(pca_scores_gcs_cumadversities$adv_cumulative)
levels(pca_scores_gcs_cumadversities$adv_cum_asfactor) <- list(zero=(0), one=(1), twoplus=c(2,3,4,5))

#Need residuals in the data set, saved as numeric

Complete.data$resids <-Residuals
Complete.data$resids <- as.numeric(Complete.data$resids)

#Need to save process data into new data frame with sname in it

Process.data_gcs=data.frame(
  sname=rep(name.list,each=length(age_grids)),
  grids=rep(age_grids,length(name.list)),
  process.val=as.vector(process_mean))

#Now subset into individuals to be used for graphs
# c("EPI", "LAV","LAZ", "LOC")
# c("APP", "HON","LAD", "MON")
# c("NAI", "POK")
# c("CYC", "MYS","NEA", "WOI")
##
crital_value_1_l=quantile(pca_scores_gcs$score_1,0.2)
crital_value_2_l=quantile(pca_scores_gcs$score_2,0.2)
crital_value_1_u=quantile(pca_scores_gcs$score_1,0.8)
crital_value_2_u=quantile(pca_scores_gcs$score_2,0.8)

upper_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2>crital_value_2_u]
upper_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2>crital_value_2_u]
lower_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2<crital_value_2_l]
lower_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2<crital_value_2_l]
# name_upperleft=c("EPI", "LAV")
# name_upperright=c("APP", "HON")
# name_lowerright=c("NAI","POK")
# name_lowerleft=c("CYC","MYS")
set.seed(1000)
name_upperleft=sample(as.character(upper_left_name_list),2)
name_upperright=sample(as.character(upper_right_name_list),2)
name_lowerright=sample(as.character(lower_right_name_list),2)
name_lowerleft=sample(as.character(lower_left_name_list),2)

name_combined=c(name_upperleft,name_upperright,name_lowerright,name_lowerleft)
name_order=rep(1:4,each=2)
completedata_plottigindividuals_ll <- subset(Complete.data, sname %in% name_combined)
processdata_plottingindividuals_ll <- subset(Process.data_gcs, sname %in% name_combined)
##Order_data with name
completedata_plottigindividuals_ll$order_name = unlist(lapply(
  completedata_plottigindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
completedata_plottigindividuals_ll = completedata_plottigindividuals_ll[order(completedata_plottigindividuals_ll$order_name), ]

completedata_plottigindividuals_ll$position="Q1"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==2]="Q2"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==3]="Q3"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==4]="Q4"
completedata_plottigindividuals_ll$position=as.factor(completedata_plottigindividuals_ll$position)

processdata_plottingindividuals_ll$order_name=unlist(lapply(
  processdata_plottingindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
processdata_plottingindividuals_ll = processdata_plottingindividuals_ll[order(processdata_plottingindividuals_ll$order_name), ]
#QQ
processdata_plottingindividuals_ll$position="Q1"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==2]="Q2"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==3]="Q3"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==4]="Q4"
processdata_plottingindividuals_ll$position=as.factor(processdata_plottingindividuals_ll$position)


#Size is 720x700 for export
#pdf("try.pdf")
#for (pos in c("Q1","Q2","Q3","Q4"))

breaks=c(4,8,12,16)

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q1")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q1")

p_q1=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
   labs(x = "",
        y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
#  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q2")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q2")

p_q2=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q3")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q3")

p_q3=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q4")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q4")

p_q4=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
#  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")


library(grid)
pdf("../../Figures/Multiple_Individual_Trajectories/DSI_F_individual_trajectory.pdf",width=6,height=10)
grid.arrange(p_q1,p_q2,p_q3,p_q4,nrow=4,left="Social connectedness residuals",bottom="Age at sample collection")
dev.off()

# #plot size is 720x720 using exportp
# #pdf("try.pdf")
# ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
#   geom_point()+
#   stat_ellipse()+
#   scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
#   geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
#              aes(label=sname))+
#   labs(x="Score on PC 1", y="Score on PC 2", color="Cumulative adversities")+
#   theme_bw()+
#   theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
#         legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
#         legend.text = element_text(size=20), legend.title = element_text(size=20)) 
#   #dev.off()
# #This is the above plot modified for a powerpoint presentation, w/out labeled animals

#png("try.png")
#pdf("DSI_F_Ordination.pdf")
ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
  geom_point(size=2.0)+
  stat_ellipse(size=1.0)+
  scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
  # geom_label(data=data.frame(name=c("Q1","Q2","Q3","Q4"),
  #                          x=c(min(score_1)+1,max(score_1)-1,max(score_1)-1,min(score_1)+1),
  #                          y=c(2,2,-2,-2)),aes(label=name,x=x,y=y))+
  geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
             aes(label=sname),show.legend= F )+
  # geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
  #            aes(label=sname))+
  geom_vline(xintercept = median(score_1))+
  geom_hline(yintercept = median(score_2))+
  xlim(-0.4,0.4)+
  ylim(-0.6,0.6)+
  # ylim(quantile(pca_scores_gcs_cumadversities$score_2,0.01),quantile(pca_scores_gcs_cumadversities$score_2,0.99))+
  # xlim(quantile(pca_scores_gcs_cumadversities$score_1,0.01),quantile(pca_scores_gcs_cumadversities$score_1,0.99))+
  labs(x="Score on PC 1", y="Score on PC 2", color="Number of adversities")+
  theme_bw()+
  theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
        legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
        legend.text = element_text(size=20), legend.title = element_text(size=20)) + 
  theme(axis.text.x=element_text(hjust=0.75))
#dev.off()
###############################



load("plot_gc_pc.RData")
###Second Run
##############################################
##0,1 Comparison
###Extract the residuals by minusing the mean value
Residuals=unlist(residuals)
process_mean=process_value_fitted
##create a list for baboon names
##Complete.data is the cleaned data
#name.list=unique(Complete.data$sname)
#Random Pick Up an Individual
random.name=sample(name.list,1)

library(dplyr)
library(ggplot2)
library(ggforce)

#Principal components summary

#plot size is 720x720 using export
#pdf("GC_Eigenfunction.pdf")
plot(age_grids,eigen_1,type='l',col="black",ylim=c(-1.5,2),
     ylab="Eigen Function",xlab="Age", cex.axis=1.5, cex.lab=1.5)
lines(age_grids,eigen_1_up,lty=2)
lines(age_grids,eigen_1_down,lty=2)
lines(age_grids,eigen_2,col="red")
lines(age_grids,eigen_2_up,lty=2)
lines(age_grids,eigen_2_down,lty=2)
legend("topright",legend=c(paste("1st PC",format(var_explain[1]*100,digits=4),"%"),paste("2nd PC",format(var_explain[2]*100,digits=4),"%")),lty=1,col=c("black","red"), 
       cex=1.2, bty = 'n')
#dev.off()

#Find same animals on graph that I did for DSI. Note that the quadrants are not the same, but this way
#you can follow the same animals through the entire section

#Need to put sname and adversities data into process data

pca_scores_gcs <- data.frame(score_1, score_2, Treatment_Ind,sname=unique(Complete.data$sname))
snames_cumadversities <- distinct(Complete.data, sname, adv_cumulative)
pca_scores_gcs_cumadversities <- data.frame(pca_scores_gcs, snames_cumadversities)
pca_scores_gcs_cumadversities$adv_cum_asfactor <- as.factor(pca_scores_gcs_cumadversities$adv_cumulative)
levels(pca_scores_gcs_cumadversities$adv_cum_asfactor) <- list(zero=(0), one=(1), twoplus=c(2,3,4,5))

#Need residuals in the data set, saved as numeric

Complete.data$resids <-Residuals
Complete.data$resids <- as.numeric(Complete.data$resids)

#Need to save process data into new data frame with sname in it

Process.data_gcs=data.frame(
  sname=rep(name.list,each=length(age_grids)),
  grids=rep(age_grids,length(name.list)),
  process.val=as.vector(process_mean))

#Now subset into individuals to be used for graphs
# c("EPI", "LAV","LAZ", "LOC")
# c("APP", "HON","LAD", "MON")
# c("NAI", "POK")
# c("CYC", "MYS","NEA", "WOI")
# ##
# crital_value_1_l=quantile(pca_scores_gcs$score_1,0.2)
# crital_value_2_l=quantile(pca_scores_gcs$score_2,0.2)
# crital_value_1_u=quantile(pca_scores_gcs$score_1,0.8)
# crital_value_2_u=quantile(pca_scores_gcs$score_2,0.8)
# 
# upper_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2>crital_value_2_u]
# upper_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2>crital_value_2_u]
# lower_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2<crital_value_2_l]
# lower_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2<crital_value_2_l]
# # name_upperleft=c("EPI", "LAV")
# # name_upperright=c("APP", "HON")
# # name_lowerright=c("NAI","POK")
# # name_lowerleft=c("CYC","MYS")
# set.seed(1000)
# name_upperleft=sample(as.character(upper_left_name_list),2)
# name_upperright=sample(as.character(upper_right_name_list),2)
# name_lowerright=sample(as.character(lower_right_name_list),2)
# name_lowerleft=sample(as.character(lower_left_name_list),2)
# 
# name_combined=c(name_upperleft,name_upperright,name_lowerright,name_lowerleft)
# name_order=rep(1:4,each=2)
completedata_plottigindividuals_ll <- subset(Complete.data, sname %in% name_combined)
processdata_plottingindividuals_ll <- subset(Process.data_gcs, sname %in% name_combined)
##Order_data with name
completedata_plottigindividuals_ll$order_name = unlist(lapply(
  completedata_plottigindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
completedata_plottigindividuals_ll = completedata_plottigindividuals_ll[order(completedata_plottigindividuals_ll$order_name), ]

completedata_plottigindividuals_ll$position="Q1"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==2]="Q2"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==3]="Q3"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==4]="Q4"
completedata_plottigindividuals_ll$position=as.factor(completedata_plottigindividuals_ll$position)

processdata_plottingindividuals_ll$order_name=unlist(lapply(
  processdata_plottingindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
processdata_plottingindividuals_ll = processdata_plottingindividuals_ll[order(processdata_plottingindividuals_ll$order_name), ]
#QQ
processdata_plottingindividuals_ll$position="Q1"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==2]="Q2"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==3]="Q3"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==4]="Q4"
processdata_plottingindividuals_ll$position=as.factor(processdata_plottingindividuals_ll$position)


#Size is 720x700 for export
#pdf("try.pdf")
#for (pos in c("Q1","Q2","Q3","Q4"))

breaks=c(4,8,12,16)

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q1")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q1")

p_q1=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q2")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q2")

p_q2=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q3")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q3")

p_q3=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q4")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q4")

p_q4=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")


pdf("../../Figures/Multiple_Individual_Trajectories/GC_individual_trajectory.pdf",width=6,height=10)
grid.arrange(p_q1,p_q2,p_q3,p_q4,nrow=4,left="fGC residuals",bottom="Age at sample collection")
dev.off()

# #plot size is 720x720 using exportp
# #pdf("try.pdf")
# ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
#   geom_point()+
#   stat_ellipse()+
#   scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
#   geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
#              aes(label=sname))+
#   labs(x="Score on PC 1", y="Score on PC 2", color="Cumulative adversities")+
#   theme_bw()+
#   theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
#         legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
#         legend.text = element_text(size=20), legend.title = element_text(size=20)) 
#   #dev.off()
# #This is the above plot modified for a powerpoint presentation, w/out labeled animals

#png("try.png")
#pdf("GC_Ordination.pdf")
ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
  geom_point(size=2.0)+
  stat_ellipse(size=1.0)+
  scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
  # geom_label(data=data.frame(name=c("Q1","Q2","Q3","Q4"),
  #                          x=c(min(score_1)+1,max(score_1)-1,max(score_1)-1,min(score_1)+1),
  #                          y=c(2,2,-2,-2)),aes(label=name,x=x,y=y))+
  geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
             aes(label=sname),show.legend= F )+
  # geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
  #            aes(label=sname))+
  geom_vline(xintercept = median(score_1))+
  geom_hline(yintercept = median(score_2))+
  xlim(max(min(score_1),-20)-0.05,min(max(score_1),20)+0.05)+
  ylim(max(min(score_2),-20)-0.05,min(max(score_2),20)+0.05)+
  # ylim(quantile(pca_scores_gcs_cumadversities$score_2,0.01),quantile(pca_scores_gcs_cumadversities$score_2,0.99))+
  # xlim(quantile(pca_scores_gcs_cumadversities$score_1,0.01),quantile(pca_scores_gcs_cumadversities$score_1,0.99))+
  labs(x="Score on PC 1", y="Score on PC 2", color="Number of adversities")+
  theme_bw()+
  theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
        legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
        legend.text = element_text(size=20), legend.title = element_text(size=20)) + 
  theme(axis.text.x=element_text(hjust=0.75))
#dev.off()

#####################################


load("plot_DSI_M_pc.RData")
###Third Run
##############################################
##0,1 Comparison
###Extract the residuals by minusing the mean value
Residuals=unlist(residuals)
process_mean=process_value_fitted
##create a list for baboon names
##Complete.data is the cleaned data
#name.list=unique(Complete.data$sname)
#Random Pick Up an Individual
random.name=sample(name.list,1)

library(dplyr)
library(ggplot2)
library(ggforce)

#Principal components summary

#plot size is 720x720 using export
#pdf("DSI_M_Eigenfunction.pdf")
plot(age_grids,eigen_1,type='l',col="black",ylim=c(-1.5,2),
     ylab="Eigen Function",xlab="Age", cex.axis=1.5, cex.lab=1.5)
lines(age_grids,eigen_1_up,lty=2)
lines(age_grids,eigen_1_down,lty=2)
lines(age_grids,eigen_2,col="red")
lines(age_grids,eigen_2_up,lty=2)
lines(age_grids,eigen_2_down,lty=2)
legend("topright",legend=c(paste("1st PC",format(var_explain[1]*100,digits=4),"%"),paste("2nd PC",format(var_explain[2]*100,digits=4),"%")),lty=1,col=c("black","red"), 
       cex=1.2, bty = 'n')
#dev.off()

#Find same animals on graph that I did for DSI. Note that the quadrants are not the same, but this way
#you can follow the same animals through the entire section

#Need to put sname and adversities data into process data

pca_scores_gcs <- data.frame(score_1, score_2, Treatment_Ind,sname=unique(Complete.data$sname))
snames_cumadversities <- distinct(Complete.data, sname, adv_cumulative)
pca_scores_gcs_cumadversities <- data.frame(pca_scores_gcs, snames_cumadversities)
pca_scores_gcs_cumadversities$adv_cum_asfactor <- as.factor(pca_scores_gcs_cumadversities$adv_cumulative)
levels(pca_scores_gcs_cumadversities$adv_cum_asfactor) <- list(zero=(0), one=(1), twoplus=c(2,3,4,5))

#Need residuals in the data set, saved as numeric

Complete.data$resids <-Residuals
Complete.data$resids <- as.numeric(Complete.data$resids)

#Need to save process data into new data frame with sname in it

Process.data_gcs=data.frame(
  sname=rep(name.list,each=length(age_grids)),
  grids=rep(age_grids,length(name.list)),
  process.val=as.vector(process_mean))

#Now subset into individuals to be used for graphs
# c("EPI", "LAV","LAZ", "LOC")
# c("APP", "HON","LAD", "MON")
# c("NAI", "POK")
# c("CYC", "MYS","NEA", "WOI")
# ##
# crital_value_1_l=quantile(pca_scores_gcs$score_1,0.2)
# crital_value_2_l=quantile(pca_scores_gcs$score_2,0.2)
# crital_value_1_u=quantile(pca_scores_gcs$score_1,0.8)
# crital_value_2_u=quantile(pca_scores_gcs$score_2,0.8)
# 
# upper_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2>crital_value_2_u]
# upper_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2>crital_value_2_u]
# lower_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2<crital_value_2_l]
# lower_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2<crital_value_2_l]
# # name_upperleft=c("EPI", "LAV")
# # name_upperright=c("APP", "HON")
# # name_lowerright=c("NAI","POK")
# # name_lowerleft=c("CYC","MYS")
# set.seed(1000)
# name_upperleft=sample(as.character(upper_left_name_list),2)
# name_upperright=sample(as.character(upper_right_name_list),2)
# name_lowerright=sample(as.character(lower_right_name_list),2)
# name_lowerleft=sample(as.character(lower_left_name_list),2)
# 
# name_combined=c(name_upperleft,name_upperright,name_lowerright,name_lowerleft)
# name_order=rep(1:4,each=2)
completedata_plottigindividuals_ll <- subset(Complete.data, sname %in% name_combined)
processdata_plottingindividuals_ll <- subset(Process.data_gcs, sname %in% name_combined)
##Order_data with name
completedata_plottigindividuals_ll$order_name = unlist(lapply(
  completedata_plottigindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
completedata_plottigindividuals_ll = completedata_plottigindividuals_ll[order(completedata_plottigindividuals_ll$order_name), ]

completedata_plottigindividuals_ll$position="Q1"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==2]="Q2"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==3]="Q3"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==4]="Q4"
completedata_plottigindividuals_ll$position=as.factor(completedata_plottigindividuals_ll$position)

processdata_plottingindividuals_ll$order_name=unlist(lapply(
  processdata_plottingindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
processdata_plottingindividuals_ll = processdata_plottingindividuals_ll[order(processdata_plottingindividuals_ll$order_name), ]
#QQ
processdata_plottingindividuals_ll$position="Q1"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==2]="Q2"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==3]="Q3"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==4]="Q4"
processdata_plottingindividuals_ll$position=as.factor(processdata_plottingindividuals_ll$position)


#Size is 720x700 for export
#pdf("try.pdf")
#for (pos in c("Q1","Q2","Q3","Q4"))

breaks=c(4,8,12,16)

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q1")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q1")

p_q1=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q2")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q2")

p_q2=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q3")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q3")

p_q3=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q4")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q4")

p_q4=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")


pdf("../../Figures/Multiple_Individual_Trajectories/DSI_M_individual_trajectory.pdf",width=6,height=10)
grid.arrange(p_q1,p_q2,p_q3,p_q4,nrow=4,left="Social connectedness residuals",bottom="Age at sample collection")
dev.off()

# #plot size is 720x720 using exportp
# #pdf("try.pdf")
# ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
#   geom_point()+
#   stat_ellipse()+
#   scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
#   geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
#              aes(label=sname))+
#   labs(x="Score on PC 1", y="Score on PC 2", color="Cumulative adversities")+
#   theme_bw()+
#   theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
#         legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
#         legend.text = element_text(size=20), legend.title = element_text(size=20)) 
#   #dev.off()
# #This is the above plot modified for a powerpoint presentation, w/out labeled animals

#png("try.png")
#pdf("DSI_M_Ordination.pdf")
ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
  geom_point(size=2.0)+
  stat_ellipse(size=1.0)+
  scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
  # geom_label(data=data.frame(name=c("Q1","Q2","Q3","Q4"),
  #                          x=c(min(score_1)+1,max(score_1)-1,max(score_1)-1,min(score_1)+1),
  #                          y=c(2,2,-2,-2)),aes(label=name,x=x,y=y))+
  geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
             aes(label=sname),show.legend= F )+
  # geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
  #            aes(label=sname))+
  geom_vline(xintercept = median(score_1))+
  geom_hline(yintercept = median(score_2))+
  xlim(max(min(score_1),-20)-0.05,min(max(score_1),20)+0.05)+
  ylim(max(min(score_2),-20)-0.05,min(max(score_2),20)+0.05)+
  # ylim(quantile(pca_scores_gcs_cumadversities$score_2,0.01),quantile(pca_scores_gcs_cumadversities$score_2,0.99))+
  # xlim(quantile(pca_scores_gcs_cumadversities$score_1,0.01),quantile(pca_scores_gcs_cumadversities$score_1,0.99))+
  labs(x="Score on PC 1", y="Score on PC 2", color="Number of adversities")+
  theme_bw()+
  theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
        legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
        legend.text = element_text(size=20), legend.title = element_text(size=20)) + 
  theme(axis.text.x=element_text(hjust=0.75))
#dev.off()

#####################################


load("plot_SCI_F_pc.RData")
###Fourth Run
##############################################
##0,1 Comparison
###Extract the residuals by minusing the mean value
Residuals=unlist(residuals)
process_mean=process_value_fitted
##create a list for baboon names
##Complete.data is the cleaned data
#name.list=unique(Complete.data$sname)
#Random Pick Up an Individual
random.name=sample(name.list,1)

library(dplyr)
library(ggplot2)
library(ggforce)

#Principal components summary

#plot size is 720x720 using export
#pdf("SCI_F_Eigenfunction.pdf")
plot(age_grids,eigen_1,type='l',col="black",ylim=c(-1.5,2),
     ylab="Eigen Function",xlab="Age", cex.axis=1.5, cex.lab=1.5)
lines(age_grids,eigen_1_up,lty=2)
lines(age_grids,eigen_1_down,lty=2)
lines(age_grids,eigen_2,col="red")
lines(age_grids,eigen_2_up,lty=2)
lines(age_grids,eigen_2_down,lty=2)
legend("topright",legend=c(paste("1st PC",format(var_explain[1]*100,digits=4),"%"),paste("2nd PC",format(var_explain[2]*100,digits=4),"%")),lty=1,col=c("black","red"), 
       cex=1.2, bty = 'n')
#dev.off()

#Find same animals on graph that I did for DSI. Note that the quadrants are not the same, but this way
#you can follow the same animals through the entire section

#Need to put sname and adversities data into process data

pca_scores_gcs <- data.frame(score_1, score_2, Treatment_Ind,sname=unique(Complete.data$sname))
snames_cumadversities <- distinct(Complete.data, sname, adv_cumulative)
pca_scores_gcs_cumadversities <- data.frame(pca_scores_gcs, snames_cumadversities)
pca_scores_gcs_cumadversities$adv_cum_asfactor <- as.factor(pca_scores_gcs_cumadversities$adv_cumulative)
levels(pca_scores_gcs_cumadversities$adv_cum_asfactor) <- list(zero=(0), one=(1), twoplus=c(2,3,4,5))

#Need residuals in the data set, saved as numeric

Complete.data$resids <-Residuals
Complete.data$resids <- as.numeric(Complete.data$resids)

#Need to save process data into new data frame with sname in it

Process.data_gcs=data.frame(
  sname=rep(name.list,each=length(age_grids)),
  grids=rep(age_grids,length(name.list)),
  process.val=as.vector(process_mean))

#Now subset into individuals to be used for graphs
# c("EPI", "LAV","LAZ", "LOC")
# c("APP", "HON","LAD", "MON")
# c("NAI", "POK")
# c("CYC", "MYS","NEA", "WOI")
# ##
# crital_value_1_l=quantile(pca_scores_gcs$score_1,0.2)
# crital_value_2_l=quantile(pca_scores_gcs$score_2,0.2)
# crital_value_1_u=quantile(pca_scores_gcs$score_1,0.8)
# crital_value_2_u=quantile(pca_scores_gcs$score_2,0.8)
# 
# upper_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2>crital_value_2_u]
# upper_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2>crital_value_2_u]
# lower_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2<crital_value_2_l]
# lower_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2<crital_value_2_l]
# # name_upperleft=c("EPI", "LAV")
# # name_upperright=c("APP", "HON")
# # name_lowerright=c("NAI","POK")
# # name_lowerleft=c("CYC","MYS")
# set.seed(1000)
# name_upperleft=sample(as.character(upper_left_name_list),2)
# name_upperright=sample(as.character(upper_right_name_list),2)
# name_lowerright=sample(as.character(lower_right_name_list),2)
# name_lowerleft=sample(as.character(lower_left_name_list),2)
# 
# name_combined=c(name_upperleft,name_upperright,name_lowerright,name_lowerleft)
# name_order=rep(1:4,each=2)
completedata_plottigindividuals_ll <- subset(Complete.data, sname %in% name_combined)
processdata_plottingindividuals_ll <- subset(Process.data_gcs, sname %in% name_combined)
##Order_data with name
completedata_plottigindividuals_ll$order_name = unlist(lapply(
  completedata_plottigindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
completedata_plottigindividuals_ll = completedata_plottigindividuals_ll[order(completedata_plottigindividuals_ll$order_name), ]

completedata_plottigindividuals_ll$position="Q1"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==2]="Q2"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==3]="Q3"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==4]="Q4"
completedata_plottigindividuals_ll$position=as.factor(completedata_plottigindividuals_ll$position)

processdata_plottingindividuals_ll$order_name=unlist(lapply(
  processdata_plottingindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
processdata_plottingindividuals_ll = processdata_plottingindividuals_ll[order(processdata_plottingindividuals_ll$order_name), ]
#QQ
processdata_plottingindividuals_ll$position="Q1"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==2]="Q2"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==3]="Q3"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==4]="Q4"
processdata_plottingindividuals_ll$position=as.factor(processdata_plottingindividuals_ll$position)


#Size is 720x700 for export
#pdf("try.pdf")
#for (pos in c("Q1","Q2","Q3","Q4"))

breaks=c(4,8,12,16)

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q1")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q1")

p_q1=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q2")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q2")

p_q2=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q3")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q3")

p_q3=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q4")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q4")

p_q4=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")


pdf("../../Figures/Multiple_Individual_Trajectories/SCI_F_individual_trajectory.pdf",width=6,height=10)
grid.arrange(p_q1,p_q2,p_q3,p_q4,nrow=4,left="Social connectedness residuals",bottom="Age at sample collection")
dev.off()

# #plot size is 720x720 using exportp
# #pdf("try.pdf")
# ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
#   geom_point()+
#   stat_ellipse()+
#   scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
#   geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
#              aes(label=sname))+
#   labs(x="Score on PC 1", y="Score on PC 2", color="Cumulative adversities")+
#   theme_bw()+
#   theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
#         legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
#         legend.text = element_text(size=20), legend.title = element_text(size=20)) 
#   #dev.off()
# #This is the above plot modified for a powerpoint presentation, w/out labeled animals

#pdf("SCI_F_Ordination.pdf")
ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
  geom_point(size=2.0)+
  stat_ellipse(size=1.0)+
  scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
  # geom_label(data=data.frame(name=c("Q1","Q2","Q3","Q4"),
  #                          x=c(min(score_1)+1,max(score_1)-1,max(score_1)-1,min(score_1)+1),
  #                          y=c(2,2,-2,-2)),aes(label=name,x=x,y=y))+
  geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
             aes(label=sname),show.legend= F )+
  # geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
  #            aes(label=sname))+
  geom_vline(xintercept = median(score_1))+
  geom_hline(yintercept = median(score_2))+
  xlim(max(min(score_1),-20)-0.05,min(max(score_1),20)+0.05)+
  ylim(max(min(score_2),-20)-0.05,min(max(score_2),20)+0.05)+
  # ylim(quantile(pca_scores_gcs_cumadversities$score_2,0.01),quantile(pca_scores_gcs_cumadversities$score_2,0.99))+
  # xlim(quantile(pca_scores_gcs_cumadversities$score_1,0.01),quantile(pca_scores_gcs_cumadversities$score_1,0.99))+
  labs(x="Score on PC 1", y="Score on PC 2", color="Number of adversities")+
  theme_bw()+
  theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
        legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
        legend.text = element_text(size=20), legend.title = element_text(size=20)) + 
  theme(axis.text.x=element_text(hjust=0.75))
#dev.off()

#####################################


load("plot_SCI_M_pc.RData")
###Fourth Run
##############################################
##0,1 Comparison
###Extract the residuals by minusing the mean value
Residuals=unlist(residuals)
process_mean=process_value_fitted
##create a list for baboon names
##Complete.data is the cleaned data
#name.list=unique(Complete.data$sname)
#Random Pick Up an Individual
random.name=sample(name.list,1)

library(dplyr)
library(ggplot2)
library(ggforce)

#Principal components summary

#plot size is 720x720 using export
#pdf("SCI_M_Eigenfunction.pdf")
plot(age_grids,eigen_1,type='l',col="black",ylim=c(-1.5,2),
     ylab="Eigen Function",xlab="Age", cex.axis=1.5, cex.lab=1.5)
lines(age_grids,eigen_1_up,lty=2)
lines(age_grids,eigen_1_down,lty=2)
lines(age_grids,eigen_2,col="red")
lines(age_grids,eigen_2_up,lty=2)
lines(age_grids,eigen_2_down,lty=2)
legend("topright",legend=c(paste("1st PC",format(var_explain[1]*100,digits=4),"%"),paste("2nd PC",format(var_explain[2]*100,digits=4),"%")),lty=1,col=c("black","red"), 
       cex=1.2, bty = 'n')
#dev.off()

#Find same animals on graph that I did for DSI. Note that the quadrants are not the same, but this way
#you can follow the same animals through the entire section

#Need to put sname and adversities data into process data

pca_scores_gcs <- data.frame(score_1, score_2, Treatment_Ind,sname=unique(Complete.data$sname))
snames_cumadversities <- distinct(Complete.data, sname, adv_cumulative)
pca_scores_gcs_cumadversities <- data.frame(pca_scores_gcs, snames_cumadversities)
pca_scores_gcs_cumadversities$adv_cum_asfactor <- as.factor(pca_scores_gcs_cumadversities$adv_cumulative)
levels(pca_scores_gcs_cumadversities$adv_cum_asfactor) <- list(zero=(0), one=(1), twoplus=c(2,3,4,5))

#Need residuals in the data set, saved as numeric

Complete.data$resids <-Residuals
Complete.data$resids <- as.numeric(Complete.data$resids)

#Need to save process data into new data frame with sname in it

Process.data_gcs=data.frame(
  sname=rep(name.list,each=length(age_grids)),
  grids=rep(age_grids,length(name.list)),
  process.val=as.vector(process_mean))

#Now subset into individuals to be used for graphs
# c("EPI", "LAV","LAZ", "LOC")
# c("APP", "HON","LAD", "MON")
# c("NAI", "POK")
# c("CYC", "MYS","NEA", "WOI")
# ##
# crital_value_1_l=quantile(pca_scores_gcs$score_1,0.2)
# crital_value_2_l=quantile(pca_scores_gcs$score_2,0.2)
# crital_value_1_u=quantile(pca_scores_gcs$score_1,0.8)
# crital_value_2_u=quantile(pca_scores_gcs$score_2,0.8)
# 
# upper_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2>crital_value_2_u]
# upper_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2>crital_value_2_u]
# lower_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2<crital_value_2_l]
# lower_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2<crital_value_2_l]
# # name_upperleft=c("EPI", "LAV")
# # name_upperright=c("APP", "HON")
# # name_lowerright=c("NAI","POK")
# # name_lowerleft=c("CYC","MYS")
# set.seed(1000)
# name_upperleft=sample(as.character(upper_left_name_list),2)
# name_upperright=sample(as.character(upper_right_name_list),2)
# name_lowerright=sample(as.character(lower_right_name_list),2)
# name_lowerleft=sample(as.character(lower_left_name_list),2)
# 
# name_combined=c(name_upperleft,name_upperright,name_lowerright,name_lowerleft)
# name_order=rep(1:4,each=2)
completedata_plottigindividuals_ll <- subset(Complete.data, sname %in% name_combined)
processdata_plottingindividuals_ll <- subset(Process.data_gcs, sname %in% name_combined)
##Order_data with name
completedata_plottigindividuals_ll$order_name = unlist(lapply(
  completedata_plottigindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
completedata_plottigindividuals_ll = completedata_plottigindividuals_ll[order(completedata_plottigindividuals_ll$order_name), ]

completedata_plottigindividuals_ll$position="Q1"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==2]="Q2"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==3]="Q3"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==4]="Q4"
completedata_plottigindividuals_ll$position=as.factor(completedata_plottigindividuals_ll$position)

processdata_plottingindividuals_ll$order_name=unlist(lapply(
  processdata_plottingindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
processdata_plottingindividuals_ll = processdata_plottingindividuals_ll[order(processdata_plottingindividuals_ll$order_name), ]
#QQ
processdata_plottingindividuals_ll$position="Q1"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==2]="Q2"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==3]="Q3"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==4]="Q4"
processdata_plottingindividuals_ll$position=as.factor(processdata_plottingindividuals_ll$position)


#Size is 720x700 for export
#pdf("try.pdf")
#for (pos in c("Q1","Q2","Q3","Q4"))

breaks=c(4,8,12,16)

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q1")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q1")

p_q1=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q2")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q2")

p_q2=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q3")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q3")

p_q3=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q4")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q4")

p_q4=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")


pdf("../../Figures/Multiple_Individual_Trajectories/SCI_M_individual_trajectory.pdf",width=6,height=10)
grid.arrange(p_q1,p_q2,p_q3,p_q4,nrow=4,left="Social connectedness residuals",bottom="Age at sample collection")
dev.off()

# #plot size is 720x720 using exportp
# #pdf("try.pdf")
# ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
#   geom_point()+
#   stat_ellipse()+
#   scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
#   geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
#              aes(label=sname))+
#   labs(x="Score on PC 1", y="Score on PC 2", color="Cumulative adversities")+
#   theme_bw()+
#   theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
#         legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
#         legend.text = element_text(size=20), legend.title = element_text(size=20)) 
#   #dev.off()
# #This is the above plot modified for a powerpoint presentation, w/out labeled animals

#png("try.png")
#pdf("SCI_M_Ordination.pdf")
ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
  geom_point(size=2.0)+
  stat_ellipse(size=1.0)+
  scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
  # geom_label(data=data.frame(name=c("Q1","Q2","Q3","Q4"),
  #                          x=c(min(score_1)+1,max(score_1)-1,max(score_1)-1,min(score_1)+1),
  #                          y=c(2,2,-2,-2)),aes(label=name,x=x,y=y))+
  geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
             aes(label=sname),show.legend= F )+
  # geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
  #            aes(label=sname))+
  geom_vline(xintercept = median(score_1))+
  geom_hline(yintercept = median(score_2))+
  xlim(max(min(score_1),-20)-0.05,min(max(score_1),20)+0.05)+
  ylim(-2,3)+
  # ylim(quantile(pca_scores_gcs_cumadversities$score_2,0.01),quantile(pca_scores_gcs_cumadversities$score_2,0.99))+
  # xlim(quantile(pca_scores_gcs_cumadversities$score_1,0.01),quantile(pca_scores_gcs_cumadversities$score_1,0.99))+
  labs(x="Score on PC 1", y="Score on PC 2", color="Number of adversities")+
  theme_bw()+
  theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
        legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
        legend.text = element_text(size=20), legend.title = element_text(size=20)) + 
  theme(axis.text.x=element_text(hjust=0.75))
#dev.off()

#####################################


load("plot_proprank_pc.RData")
###Fifth Run
##############################################
##0,1 Comparison
###Extract the residuals by minusing the mean value
Residuals=unlist(residuals)
process_mean=process_value_fitted
##create a list for baboon names
##Complete.data is the cleaned data
#name.list=unique(Complete.data$sname)
#Random Pick Up an Individual
random.name=sample(name.list,1)

library(dplyr)
library(ggplot2)
library(ggforce)

#Principal components summary

#plot size is 720x720 using export
#pdf("Proprank_Eigenfunction.pdf")
plot(age_grids,eigen_1,type='l',col="black",ylim=c(-1.5,2),
     ylab="Eigen Function",xlab="Age", cex.axis=1.5, cex.lab=1.5)
lines(age_grids,eigen_1_up,lty=2)
lines(age_grids,eigen_1_down,lty=2)
lines(age_grids,eigen_2,col="red")
lines(age_grids,eigen_2_up,lty=2)
lines(age_grids,eigen_2_down,lty=2)
legend("topright",legend=c(paste("1st PC",format(var_explain[1]*100,digits=4),"%"),paste("2nd PC",format(var_explain[2]*100,digits=4),"%")),lty=1,col=c("black","red"), 
       cex=1.2, bty = 'n')
#dev.off()

#Find same animals on graph that I did for DSI. Note that the quadrants are not the same, but this way
#you can follow the same animals through the entire section

#Need to put sname and adversities data into process data

pca_scores_gcs <- data.frame(score_1, score_2, Treatment_Ind,sname=unique(Complete.data$sname))
snames_cumadversities <- distinct(Complete.data, sname, adv_cumulative)
pca_scores_gcs_cumadversities <- data.frame(pca_scores_gcs, snames_cumadversities)
pca_scores_gcs_cumadversities$adv_cum_asfactor <- as.factor(pca_scores_gcs_cumadversities$adv_cumulative)
levels(pca_scores_gcs_cumadversities$adv_cum_asfactor) <- list(zero=(0), one=(1), twoplus=c(2,3,4,5))

#Need residuals in the data set, saved as numeric

Complete.data$resids <-Residuals
Complete.data$resids <- as.numeric(Complete.data$resids)

#Need to save process data into new data frame with sname in it

Process.data_gcs=data.frame(
  sname=rep(name.list,each=length(age_grids)),
  grids=rep(age_grids,length(name.list)),
  process.val=as.vector(process_mean))

#Now subset into individuals to be used for graphs
# c("EPI", "LAV","LAZ", "LOC")
# c("APP", "HON","LAD", "MON")
# c("NAI", "POK")
# c("CYC", "MYS","NEA", "WOI")
# ##
# crital_value_1_l=quantile(pca_scores_gcs$score_1,0.2)
# crital_value_2_l=quantile(pca_scores_gcs$score_2,0.2)
# crital_value_1_u=quantile(pca_scores_gcs$score_1,0.8)
# crital_value_2_u=quantile(pca_scores_gcs$score_2,0.8)
# 
# upper_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2>crital_value_2_u]
# upper_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2>crital_value_2_u]
# lower_right_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1>crital_value_1_u&pca_scores_gcs$score_2<crital_value_2_l]
# lower_left_name_list= pca_scores_gcs$sname[pca_scores_gcs$score_1<crital_value_1_l&pca_scores_gcs$score_2<crital_value_2_l]
# # name_upperleft=c("EPI", "LAV")
# # name_upperright=c("APP", "HON")
# # name_lowerright=c("NAI","POK")
# # name_lowerleft=c("CYC","MYS")
# set.seed(1000)
# name_upperleft=sample(as.character(upper_left_name_list),2)
# name_upperright=sample(as.character(upper_right_name_list),2)
# name_lowerright=sample(as.character(lower_right_name_list),2)
# name_lowerleft=sample(as.character(lower_left_name_list),2)
# 
# name_combined=c(name_upperleft,name_upperright,name_lowerright,name_lowerleft)
# name_order=rep(1:4,each=2)
completedata_plottigindividuals_ll <- subset(Complete.data, sname %in% name_combined)
processdata_plottingindividuals_ll <- subset(Process.data_gcs, sname %in% name_combined)
##Order_data with name
completedata_plottigindividuals_ll$order_name = unlist(lapply(
  completedata_plottigindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
completedata_plottigindividuals_ll = completedata_plottigindividuals_ll[order(completedata_plottigindividuals_ll$order_name), ]

completedata_plottigindividuals_ll$position="Q1"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==2]="Q2"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==3]="Q3"
completedata_plottigindividuals_ll$position[completedata_plottigindividuals_ll$order_name==4]="Q4"
completedata_plottigindividuals_ll$position=as.factor(completedata_plottigindividuals_ll$position)

processdata_plottingindividuals_ll$order_name=unlist(lapply(
  processdata_plottingindividuals_ll$sname,
  FUN = function(x) {
    name_order[which(name_combined == x)]
  }
))
processdata_plottingindividuals_ll = processdata_plottingindividuals_ll[order(processdata_plottingindividuals_ll$order_name), ]
#QQ
processdata_plottingindividuals_ll$position="Q1"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==2]="Q2"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==3]="Q3"
processdata_plottingindividuals_ll$position[processdata_plottingindividuals_ll$order_name==4]="Q4"
processdata_plottingindividuals_ll$position=as.factor(processdata_plottingindividuals_ll$position)


#Size is 720x700 for export
#pdf("try.pdf")
#for (pos in c("Q1","Q2","Q3","Q4"))

breaks=c(4,8,12,16)

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q1")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q1")

p_q1=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q2")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q2")

p_q2=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q3")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q3")

p_q3=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(labels = rep("", length(breaks)), breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")

completedata_plottigindividuals_ll_sub=subset(completedata_plottigindividuals_ll,position=="Q4")
processdata_plottingindividuals_ll_sub=subset(processdata_plottingindividuals_ll,position=="Q4")

p_q4=ggplot() +
  geom_line(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_point(data=completedata_plottigindividuals_ll_sub, aes(x=age_sample, y=resids))+
  geom_line(linetype = 2, data=processdata_plottingindividuals_ll_sub,aes(x=grids, y=process.val))+
  labs(x = "",
       y = "")+ylim(min(completedata_plottigindividuals_ll$resids)-0.1,max(completedata_plottigindividuals_ll$resids)+0.1)+
  scale_x_continuous(breaks = breaks)+
  theme_linedraw()+
  theme_bw()+
  #  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))+
  facet_grid(position~sname,switch="y")


pdf("../../Figures/Multiple_Individual_Trajectories/Proprank_individual_trajectory.pdf",width=6,height=10)
grid.arrange(p_q1,p_q2,p_q3,p_q4,nrow=4,left="Prop Rank residuals",bottom="Age at sample collection")
dev.off()

# #plot size is 720x720 using exportp
# #pdf("try.pdf")
# ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
#   geom_point()+
#   stat_ellipse()+
#   scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
#   geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
#              aes(label=sname))+
#   labs(x="Score on PC 1", y="Score on PC 2", color="Cumulative adversities")+
#   theme_bw()+
#   theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
#         legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
#         legend.text = element_text(size=20), legend.title = element_text(size=20)) 
#   #dev.off()
# #This is the above plot modified for a powerpoint presentation, w/out labeled animals

#png("try.png")
#pdf("Proprank_Ordination.pdf")
ggplot(pca_scores_gcs_cumadversities, aes(score_1, score_2, color=adv_cum_asfactor))+
  geom_point(size=2.0)+
  stat_ellipse(size=1.0)+
  scale_color_brewer(palette="Dark2", labels=c("0", "1", "2+"))+
  # geom_label(data=data.frame(name=c("Q1","Q2","Q3","Q4"),
  #                          x=c(min(score_1)+1,max(score_1)-1,max(score_1)-1,min(score_1)+1),
  #                          y=c(2,2,-2,-2)),aes(label=name,x=x,y=y))+
  geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
             aes(label=sname),show.legend= F )+
  # geom_label(data=subset(pca_scores_gcs_cumadversities, sname%in%name_combined),
  #            aes(label=sname))+
  geom_vline(xintercept = median(score_1))+
  geom_hline(yintercept = median(score_2))+
  xlim(max(min(score_1),-20)-0.05,min(max(score_1),20)+0.05)+
  ylim(-1,0.5)+
  # ylim(quantile(pca_scores_gcs_cumadversities$score_2,0.01),quantile(pca_scores_gcs_cumadversities$score_2,0.99))+
  # xlim(quantile(pca_scores_gcs_cumadversities$score_1,0.01),quantile(pca_scores_gcs_cumadversities$score_1,0.99))+
  labs(x="Score on PC 1", y="Score on PC 2", color="Number of adversities")+
  theme_bw()+
  theme(legend.justification=c(0,1), legend.position=c(0.02, 0.98),
        legend.background=element_blank(), axis.text = element_text(size=20), axis.title = element_text(size=20),
        legend.text = element_text(size=20), legend.title = element_text(size=20)) + 
  theme(axis.text.x=element_text(hjust=0.75))
#dev.off()
#
#####################################