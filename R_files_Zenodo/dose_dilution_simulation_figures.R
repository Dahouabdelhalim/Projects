
##Create figures for simulation output

#load in libraries
require(here)
require(R2jags)
require(mcmcplots)
require(ggplot2)
require(plyr)
require(dplyr)
require(tidyr)
require(lubridate)
library(cowplot)


#Set wd
library(here)
no_par <- read.csv("no_par_rev.csv")
no_par <- no_par[which(no_par$r2 > 0 & no_par$r2 < 2),]
for(i in 1:length(no_par$I2)){
  if((no_par$I2[i] + no_par$S2[i]) < 0.001){
    no_par$I2[i] <- 0
    no_par$S2[i] <- 0
  }}
no_par$focdens <- no_par$S1 + no_par$I1



###REDONE Base FIGURE###
curves1 <- read.csv("base_rev.csv")
curves1 <- curves1[which(curves1$r2 > 0 & curves1$r2 < 2),]
for(i in 1:length(curves1$I2)){
  if((curves1$I2[i] + curves1$S2[i]) < 0.001){
    curves1$I2[i] <- 0
    curves1$S2[i] <- 0
  }}
for(i in 1:length(curves1$I2)){
  if((curves1$I1[i] + curves1$S1[i]) < 0.001){
    curves1$I1[i] <- 0
    curves1$S1[i] <- 0
  }}
curves1$Prevalence <- curves1$I1/(curves1$I1 + curves1$S1)
curves1$row <- "Fixed Excretion, ??=0"
curves1$compdens <- curves1$S2 + curves1$I2
curves1$focdens <- curves1$S1 + curves1$I1
curves1$focdens <- (curves1$focdens - no_par$focdens)/no_par$focdens
curves1$nopar <- no_par$focdens

curves2 <- read.csv("exc_neg_rev.csv")
curves2 <- curves2[which(curves2$r2 > 0 & curves2$r2 < 2),]
for(i in 1:length(curves2$I2)){
  if((curves2$I2[i] + curves2$S2[i]) < 0.001){
    curves2$I2[i] <- 0
    curves2$S2[i] <- 0
  }}
for(i in 1:length(curves2$I2)){
  if((curves2$I1[i] + curves2$S1[i]) < 0.001){
    curves2$I1[i] <- 0
    curves2$S1[i] <- 0
  }}
curves2$Prevalence <- curves2$I1/(curves2$I1 + curves2$S1)
curves2$row <- "Dec. Excretion, ??=-3"
curves2$compdens <- curves2$S2 + curves2$I2
curves2$focdens <- curves2$S1 + curves2$I1
curves2$focdens <- (curves2$focdens - no_par$focdens)/no_par$focdens
curves2$nopar <- no_par$focdens


curves3 <- read.csv("exc_pos_rev.csv")
curves3 <- curves3[which(curves3$r2 > 0 & curves3$r2 < 2),]
for(i in 1:length(curves3$I2)){
  if((curves3$I2[i] + curves3$S2[i]) < 0.001){
    curves3$I2[i] <- 0
    curves3$S2[i] <- 0
  }}
for(i in 1:length(curves3$I2)){
  if((curves3$I1[i] + curves3$S1[i]) < 0.001){
    curves3$I1[i] <- 0
    curves3$S1[i] <- 0
  }}
curves3$Prevalence <- curves3$I1/(curves3$I1 + curves3$S1)
curves3$row <- "Inc. Excretion, ??=0.5"
curves3$compdens <- curves3$S2 + curves3$I2
curves3$focdens <- curves3$S1 + curves3$I1
curves3$focdens <- (curves3$focdens - no_par$focdens)/no_par$focdens
curves3$nopar <- no_par$focdens

curves <- rbind(curves1,curves2,curves3)
curves$K <- as.factor(curves$K)
curves$dilution <- as.factor(curves$dilution)
colnames(curves)[5] <- "k"

#curves$row <- factor(curves$row, levels=unique(curves$row))


prev_1 <- ggplot(curves, aes(x=compdens,y=Prevalence,color = k,linetype = dilution)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=12, face='bold'), plot.title = element_text(hjust = 0.5)) +      #Size of legend Title
  geom_line(size = 2) +
  #geom_point(data = points, aes(x=(S1+I1),y=prev), size = 2) +
  #scale_linetype_manual(values=c("solid", "solid","solid","solid"))+
  facet_wrap(as.factor(curves$row)) +
  theme(strip.text = element_text(face="bold", size=15)) +
  theme(strip.background = element_blank()) +
  scale_fill_gradient(low="blue", high="yellow") +  
  labs(x="Competitor Density",y="Focal Prevalence",title="") +
  #labs(color='k') +
  theme(legend.text=element_text(size=12)) +
  xlim(0,1.8) +
  scale_colour_manual(values = c( "#1f78b4","#a6cee3","#b2df8a")) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  guides(color = FALSE, linetype = FALSE)


prop_1 <- ggplot(curves, aes(x=compdens,y=log(P + 1),color = k,linetype = dilution)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=12, face='bold'), plot.title = element_text(hjust = 0.5)) +      #Size of legend Title
  geom_line(size = 2) +
  #geom_point(data = points, aes(x=(S1+I1),y=prev), size = 2) +
  #scale_linetype_manual(values=c("solid", "solid","solid","solid"))+
  facet_wrap(as.factor(curves$row)) +
  scale_fill_gradient(low="blue", high="yellow") +  
  labs(x="Competitor Density",y="Propagule Density",title="") +
  #labs(color='k') +
  theme(legend.text=element_text(size=12)) +
  xlim(0,1.8) +
  scale_colour_manual(values = c("#1f78b4","#a6cee3","#b2df8a")) +
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(linetype = "Competitor Excretion")


  
dens_1 <- ggplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=12, face='bold'), plot.title = element_text(hjust = 0.5)) +      #Size of legend Title
  geom_line(data = curves, aes(x=compdens,y=focdens,color = k,linetype = dilution), size = 2) +
  geom_hline(yintercept = -0.1835568, size = 1) +
  facet_wrap(~as.factor(row)) +
  scale_fill_gradient(low="blue", high="yellow") +  
  labs(x="Competitor Density",y="Focal Density",title="") +
  theme(legend.text=element_text(size=12)) +
  xlim(0,1.8) +
  scale_colour_manual(values = c( "#1f78b4","#a6cee3","#b2df8a")) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  guides(color = FALSE, linetype = FALSE) 

  

tiff(filename = "base_figure_analysis.tiff", width = 3000, height = 1800, pointsize = 12,res=300)

plot_grid(prev_1, prop_1,
          ncol = 1, nrow = 2, align = 'v', axis = "lr",
          rel_heights = c(1.1,1))

dev.off()


######################################################
#######################################################
########################################################3


###REDONE mort FIGURE###
curves1 <- read.csv("mort_dec_rev.csv")
curves1 <- curves1[which(curves1$r2 > 0 & curves1$r2 < 2),]
for(i in 1:length(curves1$I2)){
  if((curves1$I2[i] + curves1$S2[i]) < 0.001){
    curves1$I2[i] <- 0
    curves1$S2[i] <- 0
  }}
for(i in 1:length(curves1$I2)){
  if((curves1$I1[i] + curves1$S1[i]) < 0.001){
    curves1$I1[i] <- 0
    curves1$S1[i] <- 0
  }}
curves1$Prevalence <- curves1$I1/(curves1$I1 + curves1$S1)
curves1$row <- "Decelerating, ??=0.5"
curves1$compdens <- curves1$S2 + curves1$I2
curves1$focdens <- curves1$S1 + curves1$I1
curves1$focdens <- (curves1$focdens - no_par$focdens)/no_par$focdens

curves2 <- read.csv("mort_lin_rev.csv")
curves2 <- curves2[which(curves2$r2 > 0 & curves2$r2 < 2),]
for(i in 1:length(curves2$I2)){
  if((curves2$I2[i] + curves2$S2[i]) < 0.001){
    curves2$I2[i] <- 0
    curves2$S2[i] <- 0
  }}
for(i in 1:length(curves1$I2)){
  if((curves2$I1[i] + curves2$S1[i]) < 0.001){
    curves2$I1[i] <- 0
    curves2$S1[i] <- 0
  }}
curves2$Prevalence <- curves2$I1/(curves2$I1 + curves2$S1)
curves2$row <- "Linear, ??=1.0"
curves2$compdens <- curves2$S2 + curves2$I2
curves2$focdens <- curves2$S1 + curves2$I1
curves2$focdens <- (curves2$focdens - no_par$focdens)/no_par$focdens

curves3 <- read.csv("mort_acc_rev.csv")
curves3 <- curves3[which(curves3$r2 > 0 & curves3$r2 < 2),]
for(i in 1:length(curves3$I2)){
  if((curves3$I2[i] + curves3$S2[i]) < 0.001){
    curves3$I2[i] <- 0
    curves3$S2[i] <- 0
  }}
for(i in 1:length(curves1$I2)){
  if((curves3$I1[i] + curves3$S1[i]) < 0.001){
    curves3$I1[i] <- 0
    curves3$S1[i] <- 0
  }}
curves3$Prevalence <- curves3$I1/(curves3$I1 + curves3$S1)
curves3$row <- "Accelerating, ??=1.5"
curves3$compdens <- curves3$S2 + curves3$I2
curves3$focdens <- curves3$S1 + curves3$I1
curves3$focdens <- (curves3$focdens - no_par$focdens)/no_par$focdens


curves <- rbind(curves1,curves2,curves3)
curves$K <- as.factor(curves$K)
curves$dilution <- as.factor(curves$dilution)
colnames(curves)[5] <- "k"

curves$row <- factor(curves$row, levels=unique(curves$row))


prev_1 <- ggplot(curves, aes(x=compdens,y=Prevalence,color = k,linetype = dilution)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=12, face='bold'), plot.title = element_text(hjust = 0.5)) +      #Size of legend Title
  geom_line(size = 2) +
  #geom_point(data = points, aes(x=(S1+I1),y=prev), size = 2) +
  #scale_linetype_manual(values=c("solid", "solid","solid","solid"))+
  facet_wrap(as.factor(curves$row)) +
  theme(strip.text = element_text(face="bold", size=15)) +
  theme(strip.background = element_blank()) +
  scale_fill_gradient(low="blue", high="yellow") +  
  labs(x="Competitor Density",y="Focal Prevalence",title="") +
  #labs(color='k') +
  theme(legend.text=element_text(size=12)) +
  xlim(0,1.8) +
  scale_colour_manual(values = c( "#1f78b4","#a6cee3","#b2df8a")) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  guides(color = FALSE, linetype = FALSE)


prop_1 <- ggplot(curves, aes(x=compdens,y=log(P + 1),color = k,linetype = dilution)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=12, face='bold'), plot.title = element_text(hjust = 0.5)) +      #Size of legend Title
  geom_line(size = 2) +
  #geom_point(data = points, aes(x=(S1+I1),y=prev), size = 2) +
  #scale_linetype_manual(values=c("solid", "solid","solid","solid"))+
  facet_wrap(as.factor(curves$row)) +
  scale_fill_gradient(low="blue", high="yellow") +  
  labs(x="Competitor Density",y="Propagule Density",title="") +
  #labs(color='k') +
  theme(legend.text=element_text(size=12)) +
  xlim(0,1.8) +
  scale_colour_manual(values = c("#1f78b4","#a6cee3","#b2df8a")) +
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(linetype = "Competitor Excretion")



dens_1 <- ggplot(curves, aes(x=compdens,y=(S1+I1),color = k,linetype = dilution)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=12, face='bold'), plot.title = element_text(hjust = 0.5)) +      #Size of legend Title
  geom_line(size = 2) +
  #geom_hline(yintercept = -0.1835568, size = 1) +
  facet_wrap(~as.factor(row)) +
  scale_fill_gradient(low="blue", high="yellow") +  
  labs(x="Competitor Density",y="Focal Density",title="") +
  theme(legend.text=element_text(size=12)) +
  xlim(0,1.8) +
  scale_colour_manual(values = c( "#1f78b4","#a6cee3","#b2df8a")) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  guides(color = FALSE, linetype = FALSE)



tiff(filename = "mort_figure_analysis.tiff", width = 3000, height = 1800, pointsize = 12,res=300)

plot_grid(prev_1, prop_1,
          ncol = 1, nrow = 2, align = 'v', axis = "lr",
          rel_heights = c(1.1,1,1))

dev.off()


######################################################
#######################################################
########################################################3



###REDONE combine FIGURE###
curves1 <- read.csv("exc_pos_mort_dec_rev.csv")
curves1 <- curves1[which(curves1$r2 > 0 & curves1$r2 < 2),]
for(i in 1:length(curves1$I2)){
  if((curves1$I2[i] + curves1$S2[i]) < 0.001){
    curves1$I2[i] <- 0
    curves1$S2[i] <- 0
  }}
for(i in 1:length(curves1$I2)){
  if((curves1$I1[i] + curves1$S1[i]) < 0.001){
    curves1$I1[i] <- 0
    curves1$S1[i] <- 0
  }}
curves1$Prevalence <- curves1$I1/(curves1$I1 + curves1$S1)
curves1$row <- "Decelerating, ??=0.5"
curves1$compdens <- curves1$S2 + curves1$I2
curves1$focdens <- curves1$S1 + curves1$I1
curves1$focdens <- (curves1$focdens - no_par$focdens)/no_par$focdens

curves2 <- read.csv("exc_pos_mort_lin_rev.csv")
curves2 <- curves2[which(curves2$r2 > 0 & curves2$r2 < 2),]
for(i in 1:length(curves2$I2)){
  if((curves2$I2[i] + curves2$S2[i]) < 0.001){
    curves2$I2[i] <- 0
    curves2$S2[i] <- 0
  }}
for(i in 1:length(curves1$I2)){
  if((curves2$I1[i] + curves2$S1[i]) < 0.001){
    curves2$I1[i] <- 0
    curves2$S1[i] <- 0
  }}
curves2$Prevalence <- curves2$I1/(curves2$I1 + curves2$S1)
curves2$row <- "Linear, ??=1.0"
curves2$compdens <- curves2$S2 + curves2$I2
curves2$focdens <- curves2$S1 + curves2$I1
curves2$focdens <- (curves2$focdens - no_par$focdens)/no_par$focdens

curves3 <- read.csv("exc_pos_mort_acc_rev.csv")
curves3 <- curves3[which(curves3$r2 > 0 & curves3$r2 < 2),]
for(i in 1:length(curves3$I2)){
  if((curves3$I2[i] + curves3$S2[i]) < 0.001){
    curves3$I2[i] <- 0
    curves3$S2[i] <- 0
  }}
for(i in 1:length(curves1$I2)){
  if((curves3$I1[i] + curves3S1[i]) < 0.001){
    curves3$I[i] <- 0
    curves3$S1[i] <- 0
  }}
curves3$Prevalence <- curves3$I1/(curves3$I1 + curves3$S1)
curves3$row <- "Accelerating, ??=1.5"
curves3$compdens <- curves3$S2 + curves3$I2
curves3$focdens <- curves3$S1 + curves3$I1
#curves3$compdens <- round(curves1$compdens, digits = 2)
curves3$focdens <- (curves3$focdens - no_par$focdens)/no_par$focdens


curves <- rbind(curves1,curves2,curves3)
curves$K <- as.factor(curves$K)
curves$dilution <- as.factor(curves$dilution)
colnames(curves)[5] <- "k"
curves <- curves[which(curves$S1 > 0.001 | curves$I1 > 0.001),]

curves$row <- factor(curves$row, levels=unique(curves$row))


prev_1 <- ggplot(curves, aes(x=compdens,y=Prevalence,color = k,linetype = dilution)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=12, face='bold'), plot.title = element_text(hjust = 0.5)) +      #Size of legend Title
  geom_line(size = 2) +
  #geom_point(data = points, aes(x=(S1+I1),y=prev), size = 2) +
  #scale_linetype_manual(values=c("solid", "solid","solid","solid"))+
  facet_wrap(as.factor(curves$row)) +
  theme(strip.text = element_text(face="bold", size=15)) +
  theme(strip.background = element_blank()) +
  scale_fill_gradient(low="blue", high="yellow") +  
  labs(x="Competitor Density",y="Focal Prevalence",title="") +
  #labs(color='k') +
  theme(legend.text=element_text(size=12)) +
  xlim(0,1.8) +
  scale_colour_manual(values = c( "#1f78b4","#a6cee3","#b2df8a")) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  guides(color = FALSE, linetype = FALSE)


prop_1 <- ggplot(curves, aes(x=compdens,y=log(P + 1),color = k,linetype = dilution)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=12, face='bold'), plot.title = element_text(hjust = 0.5)) +      #Size of legend Title
  geom_line(size = 2) +
  #geom_point(data = points, aes(x=(S1+I1),y=prev), size = 2) +
  #scale_linetype_manual(values=c("solid", "solid","solid","solid"))+
  facet_wrap(as.factor(curves$row)) +
  scale_fill_gradient(low="blue", high="yellow") +  
  labs(x="Competitor Density",y="Propagule Density",title="") +
  #labs(color='k') +
  theme(legend.text=element_text(size=12)) +
  xlim(0,1.8) +
  scale_colour_manual(values = c("#1f78b4","#a6cee3","#b2df8a")) +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(linetype = "Competitor Excretion")



dens_1 <- ggplot(curves, aes(x=compdens,y=focdens,color = k,linetype = dilution)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=12, face='bold'), plot.title = element_text(hjust = 0.5)) +      #Size of legend Title
  geom_line(size = 2) +
  geom_hline(yintercept = -0.1835568, size = 1) +
  facet_wrap(~as.factor(row)) +
  scale_fill_gradient(low="blue", high="yellow") +  
  labs(x="Competitor Density",y="Focal Density",title="") +
  theme(legend.text=element_text(size=12)) +
  xlim(0,1.8) +
  scale_colour_manual(values = c( "#1f78b4","#a6cee3","#b2df8a")) +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  guides(color = FALSE, linetype = FALSE)



tiff(filename = "mort_exc_figure_analysis.tiff", width = 3000, height = 1800, pointsize = 12,res=300)

plot_grid(prev_1, prop_1,
          ncol = 1, nrow = 2, align = 'v', axis = "lr",
          rel_heights = c(1.1,1))

dev.off()


#######################################################################
#######################################################################
###    
###           FRIENDLY COMPETITION Figure
###
#######################################################################
#######################################################################

friend1 <- read.csv("friendly_competition1.csv",header = TRUE)
colnames(friend1) <- c("draw","r2","inter","dilution", "K", "S1", "I1","S2","I2","P","rho","gam")

#calculate focal and comp densities
friend$focdens <- friend$S1 + friend$I1
friend$comdens <- friend$S2 + friend$I2
for(i in 1:length(friend$draw)){if(friend$r2[i] == 0){friend$comdens[i] = 0}}

#calculate difference in focal densities from when competitor density is 0
index = 1
for(i in 1:length(friend$draw)){
  if(friend$r2[i] == 0){index = i}
  friend$comdiff[i] <- friend$focdens[i] - friend$focdens[index] 
  if(friend$focdens[index] > 0.999){friend$comdiff[i] = 1}
}

friend$comdiff <- round(friend$comdiff,10)

#aggregate min and max

min_friend <- aggregate(friend,by = list(friend$inter,friend$dilution,friend$K,friend$rho,friend$gam),FUN = min)
max_friend <- aggregate(friend,by = list(friend$inter,friend$dilution,friend$K,friend$rho,friend$gam),FUN = max)

##Add max onto min sheet and then mark things as full negative, full positive, or mix

min_friend$min <- min_friend$comdiff
min_friend$max <- max_friend$comdiff

for(i in 1:length(min_friend$min)){
  if(min_friend$min[i] >= 0 & min_friend$max[i] >= 0){min_friend$result[i] = 1} #all friendly competition
  if(min_friend$min[i] < 0 & min_friend$max[i] > 0){min_friend$result[i] = 1} #some friendly competition
  if(min_friend$min[i] <= 0 & min_friend$max[i] <= 0){min_friend$result[i] = 0} #all negative, no friendly competition
}



##OK, now average over k to get an average for how many scenarios are diluting
min_friend <- aggregate(min_friend,by = list(min_friend$inter,min_friend$dilution,min_friend$rho,min_friend$gam),FUN = mean)

min_friend$result <- as.factor(min_friend$result)




min_friend$result <- as.factor(min_friend$result)
levels(min_friend$result) <- c("None","k=1.5","k=1.5,1.0","k=1.5,1.0,0.5")

min_friend <- min_friend[,-c(1:9)]

min_friend$rho <- as.factor(min_friend$rho)
levels(min_friend$rho) <- c("?? = 0","?? = 0.5","?? = 1.0","?? = 1.5")

min_friend$gam <- as.factor(min_friend$gam)
levels(min_friend$gam) <- c("?? = -3","?? = 0","?? = 0.5")

tiff(filename = "friendly_competition_figure.tiff", width = 2500, height = 2000, pointsize = 12,res=300)

ggplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=12, face='bold'), plot.title = element_text(hjust = 0.5)) +       #Size of legend Title
  geom_tile(data = min_friend, aes(x = dilution, y = inter, fill = result)) +
  facet_grid(rho ~ gam) +
  theme(strip.text = element_text(face="bold", size=15)) +
  theme(strip.background = element_blank()) +
  labs(y="Interspecific Competition (??12,??21)",x="Competitor Excretion (x2)",title="",fill = "Friendly Competition") +
  scale_fill_manual(values = c( "black","#b2df8a","#a6cee3","#1f78b4")) +
  theme(legend.text=element_text(size=12))
  

dev.off()


