

#Here we make meta-analysis figures from the simulation

#load in libraries
require(here)
require(R2jags)
require(mcmcplots)
require(ggplot2)
require(plyr)
require(dplyr)
require(tidyr)
require(lubridate)
library(here)


##Read in meta-results

meta_results <- read.csv("meta_results_flex_fed.csv")

for(i in 1:length(meta_results$X)){
  if(meta_results$k.min[i] > 1){meta_results$type[i] <- "Accelerating"}
  if(meta_results$k.max[i] < 1){meta_results$type[i] <- "Decelerating"}
  if(meta_results$k.max[i] > 1 & meta_results$k.min[i] < 1){meta_results$type[i] <- "Linear"}
  
}

##Make figure 2
tiff(filename = "meta_figure_flex_feed.tiff", width = 1200, height = 1200, pointsize = 12,res=300)

ggplot(meta_results, aes(x=k.mean,y=reorder(X,k.mean),color=type)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=12),axis.text.x  = element_text(size=8)) + 
  theme(axis.title.y = element_text(face="bold", size=12),axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  theme(title =element_text(size=10, face='bold'), plot.title = element_text(hjust = 0.5)) +
  geom_errorbarh(size=0.5, height=0, aes(xmin = k.min, xmax = k.max)) +
  geom_point(size=1, shape = 18) +
  #geom_point(aes(x=sig,y=reorder(X,k.mean)),color="black") +
  geom_vline(xintercept = 1, size = 1.0, linetype = 2) +
  scale_colour_manual(values = c("#b2df8a", "#1f78b4", "#a6cee3")) +
  labs(x="Dose Shape Parameter (k)",y="Host/Pathogen Pair",title="Antagonism         Synergism") +
  coord_cartesian(xlim=c(0,2))

dev.off()


#Now we are repeating this process while holding feeding rate constant, to visualize
##Whether feeding rate has an effect on values of k. 

meta_results <- read.csv("meta_results_high_f_fed.csv")

for(i in 1:length(meta_results$X)){
  if(meta_results$k.min[i] > 1){meta_results$type[i] <- "Accelerating"}
  if(meta_results$k.max[i] < 1){meta_results$type[i] <- "Decelerating"}
  if(meta_results$k.max[i] > 1 & meta_results$k.min[i] < 1){meta_results$type[i] <- "Linear"}
  
}

tiff(filename = "meta_figure_high_f_feed.tiff", width = 1200, height = 1200, pointsize = 12,res=300)

ggplot(meta_results, aes(x=k.mean,y=reorder(X,k.mean),color=type)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=12),axis.text.x  = element_text(size=8)) + 
  theme(axis.title.y = element_text(face="bold", size=12),axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  theme(title =element_text(size=10, face='bold'), plot.title = element_text(hjust = 0.5)) +
  geom_errorbarh(size=0.5, height=0, aes(xmin = k.min, xmax = k.max)) +
  geom_point(size=1, shape = 18) +
  #geom_point(aes(x=sig,y=reorder(X,k.mean)),color="black") +
  geom_vline(xintercept = 1, size = 1.0, linetype = 2) +
  scale_colour_manual(values = c("#b2df8a", "#1f78b4", "#a6cee3")) +
  labs(x="Dose Shape Parameter (k)",y="Host/Pathogen Pair",title="Antagonism         Synergism") +
  coord_cartesian(xlim=c(0,2))

dev.off()

meta_results <- read.csv("meta_results_med_f_fed.csv")

for(i in 1:length(meta_results$X)){
  if(meta_results$k.min[i] > 1){meta_results$type[i] <- "Accelerating"}
  if(meta_results$k.max[i] < 1){meta_results$type[i] <- "Decelerating"}
  if(meta_results$k.max[i] > 1 & meta_results$k.min[i] < 1){meta_results$type[i] <- "Linear"}
  
}

tiff(filename = "meta_figure_med_f_feed.tiff", width = 1200, height = 1200, pointsize = 12,res=300)

ggplot(meta_results, aes(x=k.mean,y=reorder(X,k.mean),color=type)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=12),axis.text.x  = element_text(size=8)) + 
  theme(axis.title.y = element_text(face="bold", size=12),axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  theme(title =element_text(size=10, face='bold'), plot.title = element_text(hjust = 0.5)) +
  geom_errorbarh(size=0.5, height=0, aes(xmin = k.min, xmax = k.max)) +
  geom_point(size=1, shape = 18) +
  #geom_point(aes(x=sig,y=reorder(X,k.mean)),color="black") +
  geom_vline(xintercept = 1, size = 1.0, linetype = 2) +
  scale_colour_manual(values = c("#b2df8a", "#1f78b4", "#a6cee3")) +
  labs(x="Dose Shape Parameter (k)",y="Host/Pathogen Pair",title="Antagonism         Synergism") +
  coord_cartesian(xlim=c(0,2))

dev.off()

meta_results <- read.csv("meta_results_low_f_fed.csv")

for(i in 1:length(meta_results$X)){
  if(meta_results$k.min[i] > 1){meta_results$type[i] <- "Accelerating"}
  if(meta_results$k.max[i] < 1){meta_results$type[i] <- "Decelerating"}
  if(meta_results$k.max[i] > 1 & meta_results$k.min[i] < 1){meta_results$type[i] <- "Linear"}
  
}

tiff(filename = "meta_figure_low_f_feed.tiff", width = 1200, height = 1200, pointsize = 12,res=300)

ggplot(meta_results, aes(x=k.mean,y=reorder(X,k.mean),color=type)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) + 
  theme(axis.title.x = element_text(face="bold", size=12),axis.text.x  = element_text(size=8)) + 
  theme(axis.title.y = element_text(face="bold", size=12),axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
  theme(title =element_text(size=10, face='bold'), plot.title = element_text(hjust = 0.5)) +
  geom_errorbarh(size=0.5, height=0, aes(xmin = k.min, xmax = k.max)) +
  geom_point(size=1, shape = 18) +
  #geom_point(aes(x=sig,y=reorder(X,k.mean)),color="black") +
  geom_vline(xintercept = 1, size = 1.0, linetype = 2) +
  scale_colour_manual(values = c("#b2df8a", "#1f78b4", "#a6cee3")) +
  labs(x="Dose Shape Parameter (k)",y="Host/Pathogen Pair",title="Antagonism         Synergism") +
  coord_cartesian(xlim=c(0,2))

dev.off()
