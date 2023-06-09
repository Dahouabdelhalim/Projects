## ant morphological traits

## -- load libraries ####
library(here)
library(tidyverse)
library(reshape)
library(grid)

## -- set file paths ####
raw_data_path <- here::here("data/raw_data/csv")
rdata_path <- here::here("data/rdata")
figure_path <- here::here("figs")

## -- read in data file ####

Ruzi_etal_data_ant_morphological_measurements <- read_csv(paste(raw_data_path, "Ruzi_etal_data_ant_morphological_measurements.CSV", sep = "/"),
                                                       col_types = cols(
                                                         .default = col_double(),
                                                         Species = col_character(),
                                                         Sp_ID = col_character(),
                                                         Pitfalls_highest_cutoff_met = col_character(),
                                                         Above_seed_removal = col_character(),
                                                         Experiment = col_character(),
                                                         Specimen_ID = col_character()
                                                       ))
Ruzi_etal_data_ant_morphological_measurements

## this file has all the ant species that were found in at least 10% of traps (in at least
# 4 traps) at one site with the wet and dry seasons pooled and all the ant species 
# that were in the samples associated with hourly change in seed count and observed removing
# seeds

## -- process the data ####

## need to make alter the experiment column that says whether the species was observed
# removing seeds or was "common" in pitfall traps

mydata2 <- Ruzi_etal_data_ant_morphological_measurements %>%
  select(Sp_ID, Experiment, Num_Individuals, HW_avg, HL_avg, ML_avg, SL_avg, EW_avg,
         EL_avg, HF_avg, WL_avg, EP_avg) %>%
  # to make it into only two categories instead of three
  mutate(Experiment = if_else(Experiment == "both", "removal", Experiment))
mydata2

## log transform the data
log_data <- log(mydata2[,c("HW_avg", "HL_avg", "ML_avg", "SL_avg",
                           "EW_avg", "EL_avg", "HF_avg", "WL_avg", "EP_avg")])
log_data

log_data.matrix <- as.matrix(log_data)
data_experiment <- mydata2$Experiment

## -- make initial PCA and save it to R data ####
# commented out so not as to overwrite

#log_pca <- prcomp(log_data, center = TRUE, scale.= TRUE)
##scale.=TRUE normalized the variables to have standard deviaiton equal to 1

#save(log_pca, file = paste(rdata_path, "log_pca.rda", sep = "/"))

## -- load saved PCA data ####

load(paste(rdata_path,"log_pca.rda", sep = "/"))

round(log_pca$rotation, 4) # used to make the table in the manuscript
#         PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9
#HW_avg -0.3386 -0.2386 -0.2875 -0.3433  0.3182 -0.1911  0.4853 -0.1441 -0.4835
#HL_avg -0.3370 -0.2818 -0.2452 -0.2025 -0.6062  0.3151  0.2591  0.2590  0.3240
#ML_avg -0.3311 -0.3982 -0.2681  0.1024  0.5209  0.2101 -0.4353 -0.0982  0.3672
#SL_avg -0.3385  0.0192  0.4449  0.4620  0.0585  0.5592  0.2758 -0.1878 -0.2147
#EW_avg -0.3284  0.4697 -0.1749  0.1200 -0.0990 -0.2861  0.1803 -0.5734  0.4167
#EL_avg -0.3157  0.5746 -0.3927  0.1885  0.1091  0.0958 -0.1325  0.5404 -0.2230
#HF_avg -0.3383 -0.1574  0.4539  0.1787  0.1317 -0.5535  0.1663  0.4546  0.2546
#WL_avg -0.3416 -0.2199 -0.0010  0.1971 -0.4638 -0.2935 -0.5140 -0.1964 -0.4382
#EP_avg -0.3300  0.2843  0.4433 -0.7046  0.0492  0.1542 -0.2997 -0.0344  0.0084

plot(log_pca, type="line") # most of the variance accounted for in PC1
summary(log_pca)
#Importance of components:
#                       PC1     PC2     PC3    PC4     PC5     PC6     PC7     PC8     PC9
#Standard deviation     2.8745 0.66992 0.39620 0.2510 0.17113 0.15185 0.10223 0.06653 0.03758
#Proportion of Variance 0.9181 0.04987 0.01744 0.0070 0.00325 0.00256 0.00116 0.00049 0.00016
#Cumulative Proportion  0.9181 0.96793 0.98537 0.9924 0.99563 0.99819 0.99935 0.99984 1.00000

## -- graph PCA nicer and save the graph ####

df_out <- as.data.frame(log_pca$x)
df_out$Group <- mydata2$Experiment
df_out$species <- Ruzi_etal_data_ant_morphological_measurements$Sp_ID
head(df_out)


pca_plot_sp_updated <- ggplot(df_out, aes(x=PC1, y=PC2, label=species, color=Group))+
  geom_text(size=3)+
  ylab("PC2 (4.99%)") +
  xlab("PC1 (91.81%)") +
  xlim(c(-4,4)) +
  scale_color_manual(values = c("#999999", "red"),
                     breaks = c("pitfall", "removal"),
                     labels=c("Pitfalls\\n(n = 9)", "Removal\\n(n=13)"))+
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_rect(colour="black", fill=NA),
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(size=14, color="black"), #size of x-axis title
        axis.title.y = element_text(size=14, color="black"), #size of y-axis title
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text
        axis.text.y = element_text(size=12, color="black"))
pca_plot_sp_updated

## commenting out so not to overwrite
#ggsave("Ruzi_etal_Figure4a.png", width = 4, height = 4.5,
#       units = "in", dpi = 300, plot = pca_plot_sp_updated,
#       path = figure_path) 

## -- run one-way ANOVAs on the PCs ####
fit_PC1 <- aov(PC1 ~ Group, data=df_out)
fit_PC1

summary(fit_PC1)
#             Df Sum Sq Mean Sq F value Pr(>F)
#Group        1   2.02   2.022   0.236  0.633
#Residuals   20 171.49   8.575

layout(matrix(c(1,2,3,4),2,2)) # optional layout 
plot(fit_PC1) # diagnostic plots

fit_PC2 <- aov(PC2 ~ Group, data=df_out)
fit_PC2

summary(fit_PC2)
#             Df Sum Sq Mean Sq F value Pr(>F)  
#Group        1  1.653  1.6534   4.255 0.0524 .
#Residuals   20  7.771  0.3886                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

layout(matrix(c(1,2,3,4),2,2)) # optional layout 
plot(fit_PC2) # diagnostic plots

## homogeneity of variance tests
pc1_hov <- bartlett.test(PC1 ~ Group, data=df_out)
pc1_hov
#Bartlett test of homogeneity of variances
#data:  PC1 by Group
#Bartlett's K-squared = 3.0424, df = 1, p-value = 0.08112

pc2_hov <- bartlett.test(PC2 ~ Group, data=df_out)
pc2_hov
#Bartlett test of homogeneity of variances
#data:  PC2 by Group
#Bartlett's K-squared = 3.2407, df = 1, p-value = 0.07183

## -- graph and save the PC variances ####

df_out2 <- melt(df_out)
df_out2_subset <- subset(df_out2, variable=="PC1" | variable=="PC2")
df_out2_subset


PCA_variance <- ggplot(df_out2_subset, aes(x=variable, y=value, fill=Group))+
  geom_boxplot()+
  ylab("Loadings")+
  xlab("PC variable")+
  scale_fill_manual(values=c("#999999", "red"), # values=c("#999999", "#E69F00", "#56B4E9"),
                    breaks=c("pitfall", "removal"),
                    labels=c("Pitfalls\\n(n = 9)", "Removal\\n(n=13)"))+
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_blank(), #gets rid of square going around the entire graph
        axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="top",
        legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        strip.background=element_blank(),
        legend.title=element_blank(),
        strip.text.x=element_text(colour="black", size=12),
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(size=14, color="black"), #size of x-axis title
        axis.title.y = element_text(size=14, color="black"), #size of y-axis title
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text
        axis.text.y = element_text(size=12, color="black"),
        axis.line.x	=element_line(colour="black"),
        axis.line.y	=element_line(colour="black"))#size of y-axis text
PCA_variance

#ggsave("Ruzi_etal_Figure4b.png", width = 4, height = 4.5,
#       units = "in", dpi = 300, plot = PCA_variance,
#       path = figure_path)