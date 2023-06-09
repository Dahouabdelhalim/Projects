######################################################################################################  
#### The following code will process all simulation files and output results figures
#### Sara M. Schaal

######################################################################################################
### Load Packages and Download Data Files ###
## List Packages Needed 
packages_needed <- c("scales", "vcfR", "ggplot2", "RColorBrewer",
                      "gridExtra", "akima", "fields", "ggnewscale",
                      "ash", "plotly", "stringr", "tidyverse",  
                      "viridisLite", "ggpubr", "purrr", "dplyr")

## Install packages that aren't installed already and load them
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library( packages_needed[i], character.only = TRUE)
}

g_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

folderIn <- "./results/Inversion/20220419_lastRerun/"
folderOutFig <- "./figures/20220419/ManuscriptFigs/"
df.summary <- read.table(paste0(folderIn, "outputSumData.txt"))
colnames(df.summary) <- c("seed", "Va_perc_In", "LA_final", "num_inv", "num_inv_NA", "num_inv_NS",
                          "capture_gain_p", "capture_no_gain_p", "neutral_gain_p", "neutral_no_gain_p",
                          "capture_gain_p_NA", "capture_no_gain_p_NA", "neutral_gain_p_NA", "neutral_no_gain_p_NA",
                          "capture_gain_p_NS", "capture_no_gain_p_NS", "neutral_gain_p_NS", "neutral_no_gain_p_NS",
                          "ave_start_QTNs", "ave_start_QTNs_NA", "ave_start_QTNs_NS",
                          "ave_end_QTNs", "ave_end_QTNs_NA", "ave_end_QTNs_NS", 
                          "ave_start_FST", "ave_start_FST_NA", "ave_start_FST_NS", 
                          "ave_end_FST", "ave_end_FST_NA", "ave_end_FST_NS",  
                          "ave_abV_start_qtnSelCoef", "ave_abV_start_qtnSelCoef_NA", "ave_abV_start_qtnSelCoef_NS",
                          "ave_abV_end_qtnSelCoef", "ave_abV_end_qtnSelCoef_NA", "ave_abV_end_qtnSelCoef_NS", 
                          "true_pos_pcadapt", "true_pos_outflank", "false_neg_pcadapt", "false_neg_outflank", 
                          "true_neg_pcadapt", "true_neg_outflank", "false_pos_pcadapt", "false_pos_outflank",
                          "true_neg_pcadapt_NS",  "false_pos_pcadapt_NS", 
                          "true_neg_outflank_NS", "false_pos_outflank_NS", "av.effect", "av.perc.VA", 
                          "num.adapt.overlap.invs", "num.adapt.within.invs")


df.simStats <- read.table(paste0(folderIn, "FullSet_dfparams.txt"), header = TRUE)
df.invCharFinalGen <- read.table(paste0(folderIn, "outputInvChar_finalGen.txt"))
colnames(df.invCharFinalGen) <- c("seed", "inv_id", "inv_age", "inv_length", "num_qtns_Lscaled", "adaptInv")
df.invCharAllData <- read.table(paste0(folderIn, "outputInvChar_allData.txt"), fill = TRUE)
df.invGenome <- read.table(paste0(folderIn, "outputInvGenome.txt"))

#### END Run the following code chunk once to get full data files to do further analyses on ####
######################################################################################################  

## Sanity check this should be 0
dim(df.simStats[!df.simStats$seed %in% df.summary$seed,])

######################################################################################################  
#### Unique Parameters ####
# Subset original parameters dataframe to identify all unique parameter combinations. 
# Then using those parameter values identify the simulation seeds of replicate simulations 
# for each set of parameters.
unique.params <- unique(df.simStats[c("chromNum", "mig1", "mig2", "n", "muProp", "muInv", 
                                      "r", "alpha", "sigmaK", "theta1", "theta2", "burnin", "dom", "enVar")])

unique.params$params <- NULL
for(i in 1:nrow(unique.params)){
  unique.params$params[i] <- paste(unique.params[i,1:14], collapse = " ")
}

df.simStats$params <- NULL
for(i in 1:nrow(df.simStats)){
  df.simStats$params[i] <- paste(df.simStats[i,c(14, 9, 10, 7, 2, 3, 15, 5, 4, 11, 12, 13, 16, 8)], 
                                 collapse = " ")
}

reps <- 5
for(i in 1:reps){
  seedColName <- paste("Seed", i, sep = "")
  unique.params[,seedColName] <- NA
}

for(i in 1:nrow(unique.params)){
 seeds <- df.simStats[unique.params$params[i] == df.simStats[,17], ]$seed
 unique.params[i, 16:20] <- seeds
}

######################################################################################################    


######################################################################################################
#### Q1 & 2: Plot Percent Additive Genetic Variance in Adaptive Inversions

## Join summary data together with parameters
df.sumData <- left_join(df.summary, df.simStats[,1:(ncol(df.simStats)-1)], by = "seed")
dim(df.sumData)

## Separate data into three different inversion mutation rates
df.muInv0 <- df.sumData[df.sumData$muInv == 0,]
df.muInv3 <- df.sumData[df.sumData$muInv == 1e-03,]

## Join Inv mu 0 and with both of the other Inv mutation rates to calculate local adaptation difference
df.muInv3_0 <- full_join(df.muInv3[,c(2:3, 49:50, 53, 55:57, 59:60)], df.muInv0[,c(2:3, 49:50, 53, 55:57, 59:60)], 
                         by = c("muProp", "sigmaK", "alpha", "enVar", "mig1", "rep"))
colnames(df.muInv3_0)[c(1:4,11:14)] <- c("VA_perc_In_3", "LA_final_3", "av_qtn_effect_size_3", "av_qtn_perc_VA_3", "Va_perc_In_0", 
                                           "LA_final_0", "av_qtn_effect_size_0", "av_qtn_perc_VA_0")
df.muInv3_0$LA_diff <- df.muInv3_0$LA_final_3 - df.muInv3_0$LA_final_0

df.muInv3_0[df.muInv3_0$mig1 == 0.4 & df.muInv3_0$muProp == 0.1 & df.muInv3_0$sigmaK == 1.5 & 
              df.muInv3_0$alpha == 0.002,]

## Summarize the additive genetic variation and Local adaptation difference across replicates 
df.muInv3_0_av <- aggregate(cbind(VA_perc_In_3, LA_diff, LA_final_3, LA_final_0)~muProp + sigmaK + alpha + enVar + mig1, 
                            FUN=mean, data = df.muInv3_0)
df.muInv3_0_sd <- aggregate(cbind(VA_perc_In_3, LA_diff, LA_final_3, LA_final_0)~muProp + sigmaK + alpha + enVar + mig1, 
                            FUN=sd, data = df.muInv3_0)
df.muInv3_0_av$LA_diff_sd <-  df.muInv3_0_sd$LA_diff
df.muInv3_0_av$LA_diff_upSD <- df.muInv3_0_av$LA_diff + df.muInv3_0_av$LA_diff_sd
df.muInv3_0_av$LA_diff_lowSD <- df.muInv3_0_av$LA_diff - df.muInv3_0_av$LA_diff_sd
df.muInv3_0_av$LA_final_3_sd <-  df.muInv3_0_sd$LA_final_3
df.muInv3_0_av$LA_final_3_upSD <- df.muInv3_0_av$LA_final_3 + df.muInv3_0_av$LA_final_3_sd
df.muInv3_0_av$LA_final_3_lowSD <- df.muInv3_0_av$LA_final_3 - df.muInv3_0_av$LA_final_3_sd
df.muInv3_0_av$LA_final_0_sd <-  df.muInv3_0_sd$LA_final_0
df.muInv3_0_av$LA_final_0_upSD <- df.muInv3_0_av$LA_final_0 + df.muInv3_0_av$LA_final_0_sd
df.muInv3_0_av$LA_final_0_lowSD <- df.muInv3_0_av$LA_final_0 - df.muInv3_0_av$LA_final_0_sd
df.muInv3_0_av$VA_perc_sd <-  df.muInv3_0_sd$VA_perc_In_3
df.muInv3_0_av$VA_perc_upSD <- df.muInv3_0_av$VA_perc_In_3 + df.muInv3_0_av$VA_perc_sd
df.muInv3_0_av$VA_perc_lowSD <- df.muInv3_0_av$VA_perc_In_3 - df.muInv3_0_av$VA_perc_sd

## Plot local adaptation difference
df.muInv3_0$muProp <- as.factor(df.muInv3_0$muProp)
df.muInv3_0$mig1 <- as.factor(df.muInv3_0$mig1)
df.muInv3_0$sigmaK <- as.factor(df.muInv3_0$sigmaK)
sigmaK.labels <- c("0.75" = "Strong Selection", "1.5" = "Moderate Selection", "3" = "Weak Selection")

# make all subsetting columns factors
for(i in 1:6){
  df.muInv3_0_av[,i] <- as.factor(df.muInv3_0_av[,i])
}

# recode necessary factor levels for plotting 
df.muInv3_0_av$sigmaK <- factor(df.muInv3_0_av$sigmaK, c(3, 1.5, 0.75))
  
#### Local Adaptation Plots ####
# Difference in local adaptation
# ploygenic architecture
plot.LA_diff_inv3_pgen <- ggplot(df.muInv3_0_av[df.muInv3_0_av$enVar == 0.1 & 
                                              df.muInv3_0_av$alpha == 0.002 &
                                              df.muInv3_0_av$muProp == 0.1,], 
                            aes(x = mig1, y = LA_diff, group = sigmaK)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin=LA_diff_lowSD, ymax=LA_diff_upSD), size =  0.2, width=.5) +
  geom_point(fill = viridis(4)[2], color = "navy", shape = 21, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  labs(y = " ",
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(-0.07, 0.6)

# oligogenic architecture
plot.LA_diff_inv3_ogen <- ggplot(df.muInv3_0_av[df.muInv3_0_av$enVar == 0.1 & 
                                                df.muInv3_0_av$alpha == 0.2 &
                                                df.muInv3_0_av$muProp == 0.01,], 
                               aes(x = mig1, y = LA_diff, group = sigmaK)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin=LA_diff_lowSD, ymax=LA_diff_upSD), size =  0.2, width=.5) +
  geom_point(fill = viridis(4)[4], color = "grey13", shape = 24, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  labs(x = " ", 
       y = expression("LA"[Inv]*" - LA"[noInv])) + 
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(-0.07, 0.6)

#### Amount of Additive Genetic Variance Plots ####
colnames(df.invGenome) <- c("seed", "uniqueBases", "numOverlap", "percGenome")
df.invGenomeParam <- full_join(df.invGenome, df.simStats, by = "seed")
df.invGenome_av <- aggregate(percGenome~muProp+ sigmaK + muInv + alpha + enVar + mig1 + mig2, 
                             FUN=mean, data = df.invGenomeParam)
df.invGenome_sd <- aggregate(percGenome~muProp + sigmaK + muInv + alpha + enVar + mig1 + mig2, 
                             FUN=sd, data = df.invGenomeParam)
df.invGenome_av$percGenome_sd <- df.invGenome_sd$percGenome
df.invGenome_av$percGenome_lowSD <- df.invGenome_av$percGenome - df.invGenome_av$percGenome_sd
df.invGenome_av$percGenome_upSD <- df.invGenome_av$percGenome + df.invGenome_av$percGenome_sd

for(i in 1:6){
  df.invGenome_av[,i] <- as.factor(df.invGenome_av[,i]) 
}

df.invGenome_av$sigmaK <- factor(df.invGenome_av$sigmaK, c(3, 1.5, 0.75))
df.VA <- left_join(df.muInv3_0_av, df.invGenome_av[df.invGenome_av$muInv == 0.001,], by = c("muProp", "sigmaK", "alpha", 
                                                  "enVar", "mig1"))
df.VA$VA_perc_In_3 <- as.numeric(as.character(df.VA$VA_perc_In_3))

# polygenic architecture
plot.VA_in_inv3_pgen <- ggplot(df.VA[df.VA$enVar == 0.1 & 
                                df.VA$alpha == 0.002 &
                                df.VA$muProp == 0.1,], 
                          aes(x = mig1, y = VA_perc_In_3, group = sigmaK)) + 
  geom_errorbar(aes(ymin=VA_perc_lowSD, ymax=VA_perc_upSD), size =  0.2, width=.5) +
  geom_ribbon(aes(ymin= percGenome_lowSD, ymax= percGenome_upSD), fill = "navy", alpha=0.2) +
  geom_point(fill = viridis(4)[2], color = "navy", shape = 21, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
  labs(y = " ",
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[2]) + 
  ylim(c(0, 80))

# oligogenic architecture
plot.VA_in_inv3_ogen <- ggplot(df.VA[df.VA$enVar == 0.1 & 
                                             df.VA$alpha == 0.2 &
                                             df.VA$muProp == 0.01,], 
                          aes(x = mig1, y = VA_perc_In_3, group = sigmaK)) + 
  geom_errorbar(aes(ymin=VA_perc_lowSD, ymax=VA_perc_upSD), size = 0.2, width=.5) +
  geom_ribbon(aes(ymin= percGenome_lowSD, ymax= percGenome_upSD), fill = "goldenrod", alpha=0.2) +
  geom_point(fill = viridis(4)[4], color = "grey13", shape = 24, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
  labs(y = expression("%VA"[inv]),
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[4]) + 
  ylim(c(0, 80))
         

## total LA
# polygenic architecture
plot.totalLA_pgen <- ggplot(df.muInv3_0_av[df.muInv3_0_av$enVar == 0.1 & 
                                           df.muInv3_0_av$alpha == 0.002 &
                                           df.muInv3_0_av$muProp == 0.1,], 
                            aes(x = mig1, y = LA_final_3, group = sigmaK)) + 
  geom_errorbar(aes(ymin=LA_final_3_lowSD, ymax=LA_final_3_upSD), size = 0.2, width=.5) +
  geom_ribbon(aes(ymin= LA_final_0_lowSD, ymax= LA_final_0_upSD), fill = "navy", alpha=0.2) +
  geom_point(fill = viridis(4)[2], color = "navy", shape = 21, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
  labs(title = "Highly Polygenic Architecture",
       y = " ",
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(c(-0.015,1))

plot.totalLA_ogen <- ggplot(df.muInv3_0_av[df.muInv3_0_av$enVar == 0.1 & 
                                           df.muInv3_0_av$alpha == 0.2 &
                                           df.muInv3_0_av$muProp == 0.01,], 
                            aes(x = mig1, y = LA_final_3, group = sigmaK)) + 
  geom_errorbar(aes(ymin=LA_final_3_lowSD, ymax=LA_final_3_upSD), size = 0.2, width=.5) +
  geom_ribbon(aes(ymin= LA_final_0_lowSD, ymax= LA_final_0_upSD), fill = "goldenrod", alpha=0.2) +
  geom_point(fill = viridis(4)[4], color = "grey13", shape = 24, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
  labs(title = "Polygenic Architecture",
       y = expression("LA"[inv]),
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(c(-0.015,1))

#### Number of Inversions Plots ####
numInv_av <- aggregate(num_inv~muProp + sigmaK + alpha + enVar + mig1, 
                       FUN=mean, data = df.muInv3)
numInv_sd <- aggregate(num_inv~muProp + sigmaK + alpha + enVar + mig1, 
                       FUN=sd, data = df.muInv3)

numInv_av$numInv_sd <- numInv_sd$num_inv
numInv_av$numInv_lowSD <- numInv_av$num_inv - numInv_av$numInv_sd
numInv_av$numInv_upSD <- numInv_av$num_inv + numInv_av$numInv_sd

for(i in 1:5 ){
  numInv_av[,i] <- as.factor(numInv_av[,i])
}
for(i in (ncol(df.muInv3)-14):ncol(df.muInv3)){
  df.muInv3[,i] <- as.factor(df.muInv3[,i])
}

numInv_av$sigmaK <- factor(numInv_av$sigmaK, c(3, 1.5, 0.75))
df.muInv3$sigmaK <- factor(df.muInv3$sigmaK, c(3, 1.5, 0.75))
numInv_av$numInv_lowSD_0 <- numInv_av$numInv_lowSD
numInv_av$numInv_lowSD_0[numInv_av$numInv_lowSD <  0 ] <- 0
plot.numInv.pgen <- ggplot(numInv_av[numInv_av$enVar == 0.1 & numInv_av$alpha == 0.002 & 
                                  numInv_av$muProp == 0.1 ,], 
                      aes(x = mig1, y = num_inv, group = sigmaK)) + 
  geom_errorbar(aes(ymin=numInv_lowSD_0, ymax=numInv_upSD), size =  0.2, width=.5) +
  geom_point(fill = viridis(4)[2], color = "navy", shape = 21, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) )+ 
  labs(y = " ",
       x = "Migration Rate") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  scale_color_manual(values=viridis(4)[2]) + ylim(0, 15)

plot.numInv.ogen <- ggplot(numInv_av[numInv_av$enVar == 0.1 & numInv_av$alpha == 0.2
                                           & numInv_av$muProp == 0.01,], 
                                 aes(x = mig1, y = num_inv, group = sigmaK)) + 
  geom_errorbar(aes(ymin=numInv_lowSD_0, ymax=numInv_upSD), size =  0.2, width=.5) +
  geom_point(fill = viridis(4)[4], color = "grey13", shape = 24, size = 3) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels) ) + 
  labs(y = expression(bar(N)[inv]),  x = "Migration Rate") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=18,face="bold"))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(0, 15) +
  scale_x_discrete(breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)


pdf(file = paste0(folderOutFig, "Fig2_LAplots_envar.pdf"), width = 15, height = 15)
ggarrange(plot.totalLA_ogen, plot.totalLA_pgen, plot.LA_diff_inv3_ogen, plot.LA_diff_inv3_pgen,
          plot.VA_in_inv3_ogen, plot.VA_in_inv3_pgen, plot.numInv.ogen, plot.numInv.pgen, 
            nrow = 4, ncol = 2, labels = c("A", "B", "C", "D", "E", "F", "G", "H"))
dev.off()
 
#### end Q1 & 2 plotting
######################################################################################################

######################################################################################################
#### Q3: Plot Proportion of Adaptive inversions in each evo history category

head(df.muInv3)
dim(df.muInv3)
df.muInv3$capture_gain_count <- df.muInv3$num_inv*df.muInv3$capture_gain_p
df.muInv3$capture_no_gain_count <- df.muInv3$num_inv*df.muInv3$capture_no_gain_p
df.muInv3$neutral_gain_count <- df.muInv3$num_inv*df.muInv3$neutral_gain_p
df.muInv3$neutral_no_gain_count <- df.muInv3$num_inv*df.muInv3$neutral_no_gain_p

df.ev.hist.av <- aggregate(cbind(capture_gain_count, capture_no_gain_count, 
                neutral_gain_count)~muProp + sigmaK + alpha + enVar + mig1 + mig2, 
          FUN=mean, data = df.muInv3)
df.ev.hist.sd <- aggregate(cbind(capture_gain_count, capture_no_gain_count, 
                                 neutral_gain_count)~muProp + sigmaK + alpha + enVar + mig1 + mig2, 
                           FUN=sd, data = df.muInv3)

df.ev.hist.long <- as.data.frame(pivot_longer(df.ev.hist.av, cols = c(capture_gain_count, capture_no_gain_count, 
                                     neutral_gain_count),
             names_to = "evHist", values_to = "count"))
df.ev.hist.sd.long <- as.data.frame(pivot_longer(df.ev.hist.sd, cols = c(capture_gain_count, capture_no_gain_count, 
                                                                      neutral_gain_count),
                                              names_to = "evHist", values_to = "sd"))

for(i in 1:7){
  df.ev.hist.long[,i] <- as.factor(df.ev.hist.long[,i])
  df.ev.hist.sd.long[,i] <- as.factor(df.ev.hist.sd.long[,i])
}

df.ev.hist.long <- full_join(df.ev.hist.long, df.ev.hist.sd.long, by = c("muProp", "sigmaK", "alpha", "enVar",  "mig1",  "mig2", "evHist"))
df.ev.hist.long$upSD <- df.ev.hist.long$count + df.ev.hist.long$sd
df.ev.hist.long$lowSD <-  df.ev.hist.long$count - df.ev.hist.long$sd
df.ev.hist.long$lowSD[df.ev.hist.long$lowSD<0] <- 0
df.ev.hist.long$evHist <- recode_factor(df.ev.hist.long$evHist, 'capture_gain_count' = 'Capture & Gain', 
                                        'capture_no_gain_count' = 'Capture & No Gain', 
                                        'neutral_gain_count' = 'Neutral & Gain')
df.ev.hist.long$sigmaK <- factor(df.ev.hist.long$sigmaK, c(3, 1.5, 0.75))

df.ev.hist.boxplot <- as.data.frame(pivot_longer(df.muInv3[68:70], cols = c(capture_gain_count, capture_no_gain_count, 
                                     neutral_gain_count),
             names_to = "evHist", values_to = "count"))
df.ev.hist.boxplot$evHist <- recode_factor(df.ev.hist.boxplot$evHist, 'capture_gain_count' = 'Capture & Gain', 
                                        'capture_no_gain_count' = 'Capture & No Gain', 
                                        'neutral_gain_count' = 'Neutral & Gain')

pdf(file = paste0(folderOutFig, "fig3_evoHist.pdf"), width = 5, height = 4)

ggplot(df.ev.hist.boxplot[complete.cases(df.ev.hist.boxplot),], 
       aes(x = as.factor(evHist), y = count, fill = as.factor(evHist))) + 
  geom_boxplot() + 
  scale_fill_manual(values = viridis(5)[c(4,2,1)]) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.title=element_text(size=16),
        axis.text.x = element_text(angle = 90))  +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Inversion Evolutionary History", y = "Number of Inversions", x = "") 
dev.off()

#### end Q3 plotting
######################################################################################################

######################################################################################################
#### Q4: Plot Characteristics

head(df.invCharFinalGen)
tail(df.invCharFinalGen)
for(i in 2:ncol(df.simStats)){
  df.simStats[,i] <- as.factor(as.character(df.simStats[,i]))
}
df.invCharParam <- left_join(df.invCharFinalGen, df.simStats, by = "seed")
df.invChar.muInv3 <- df.invCharParam[df.invCharParam$muInv == 0.001,]
df.invChar.muInv3$sigmaK <- factor(df.invChar.muInv3$sigmaK, c(3, 1.5, 0.75))

color_scale <- c(inferno(4)[3], "darkgrey", inferno(4)[1])
df.invChar.muInv3$adaptInv <- factor(df.invChar.muInv3$adaptInv, levels = c("Adaptive", "Nonadaptive", "No selection"))

## subset to get data for only parameter sets that contain all 3 factor levels
plot.output <- NULL
for(i in 1:length(unique(df.invChar.muInv3$params))){
  df <- df.invChar.muInv3[df.invChar.muInv3$params==unique(df.invChar.muInv3$params)[i],]
  if(length(unique(df$adaptInv)) == 3){
    plot.output <- rbind(plot.output, df)
  }
}
plot.output$inv_length_cm <- NULL
plot.output$inv_length_cm <- plot.output$inv_length*0.0001

plot.age.pgen <- ggplot(data = plot.output[plot.output$enVar == 0.1 & 
                                                      plot.output$muProp == 0.1 & 
                                                      plot.output$alpha == 0.002,], 
  aes(x = mig1, y= inv_age, fill = adaptInv, color = adaptInv)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  scale_color_manual(values = c("orange", "black", "grey21")) +
  scale_fill_manual(values = color_scale) + 
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size=18),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=30),
        title = element_text(size = 25))+
  theme(legend.text=element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Highly Polygenic Architecture", y = " ", x = " ")

plot.age.ogen <- ggplot(data = plot.output[plot.output$enVar == 0.1 &
                                                   plot.output$muProp == 0.01 &
                                                   plot.output$alpha == 0.2,], 
       aes(x = mig1, y= inv_age, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values = c("orange", "black", "grey21"))+
  scale_fill_manual(values = color_scale) + 
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size=18),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=30),
        title = element_text(size = 25))+
  theme(legend.text=element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Polygenic Architecture", y = expression("Age"[inv]*" (gen)"), x = " ") +
  scale_x_discrete(" ", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)


plot.length.pgen <- ggplot(data = plot.output[plot.output$enVar == 0.1 & 
                                                         plot.output$muProp == 0.1 & 
                                                         plot.output$alpha == 0.002,], 
       aes(x = mig1, y= inv_length_cm, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK,  labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = color_scale) + 
  scale_color_manual(values = c("orange", "black", "grey21"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme(strip.text = element_text(size=18),
        axis.text.y = element_text(size=25),
        axis.text.x = element_text(size=25),
        axis.title = element_text(size=30),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=12),
        legend.title = element_text(size=15))+
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = " ")

plot.length.ogen <- ggplot(data = plot.output[plot.output$enVar == 0.1 & 
                                                      plot.output$muProp == 0.01 & 
                                                      plot.output$alpha == 0.2,], 
       aes(x = mig1, y= inv_length_cm, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK,  labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = color_scale) + 
  scale_color_manual(values = c("orange", "black", "grey21"))+
  theme_classic() +
  theme(legend.position = "none") +
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme(strip.text = element_text(size=18),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=30),
        title = element_text(size = 25))+
  theme(legend.text=element_text(size=17)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = expression("Length"[inv]*" (cM)"), x = " ") +
  scale_x_discrete("", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)

plot.numQTN.pgen <- ggplot(data = plot.output[plot.output$enVar == 0.1 & 
                                                         plot.output$muProp == 0.1 & 
                                                         plot.output$alpha == 0.002,], 
       aes(x = mig1, y= num_qtns_Lscaled, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK,  labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = color_scale) + 
  scale_color_manual(values = c("orange", "black", "grey21")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme(strip.text = element_text(size=18),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=30),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = "Migration Rate") + 
  ylim(c(0,NA))

plot.numQTN.ogen <- ggplot(data = plot.output[plot.output$enVar == 0.1  &
                                                      plot.output$muProp == 0.01 &
                                                      plot.output$alpha == 0.2,], 
                              aes(x = mig1, y= num_qtns_Lscaled, fill = adaptInv, color = adaptInv)) +
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = color_scale) + 
  scale_color_manual(values = c("orange", "black", "grey21")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme(strip.text = element_text(size=18),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.title = element_text(size=30),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = expression(bar(N)[QTNs] / "(Length"[inv]*")")) + 
  ylim(c(0,0.012)) +
  scale_x_discrete("Migration Rate", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)


## blank graph
invChar.legend <- g_legend(plot.length.pgen)
plot.length.pgen.noLeg <- plot.length.pgen + theme(legend.position = "none")
blank <- ggplot() + theme_void()
pdf(file = paste0(folderOutFig, "SFigX_characteristics.pdf"), width = 20, height = 17)
ggarrange(plot.age.ogen, plot.age.pgen, blank, plot.length.ogen, plot.length.pgen.noLeg, 
          invChar.legend, plot.numQTN.ogen, plot.numQTN.pgen, blank, ncol = 3, nrow = 3,
          widths = c(2.3, 2.3, 0.8, 2.3, 2.3, 0.8, 2.3, 2.3, 0.8), 
          labels = c("A", "B", "", "C", "D", "", "E", "F"),
          font.label = list(size = 30))
dev.off()

plot.age.onePanel <- ggplot(data = plot.output[plot.output$enVar == 0.1,], 
                        aes(x = adaptInv, y= inv_age, fill = adaptInv, color = adaptInv)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values = c("orange", "black", "grey21")) +
  scale_fill_manual(values = color_scale) + 
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        title = element_text(size = 17))+
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "All Simulations", y = expression("Age"[inv]*" (gen)"), x = " ")


plot.length.onePanel <- ggplot(data = plot.output[plot.output$enVar == 0.1,], 
                                   aes(x = adaptInv, y= inv_length_cm, fill = adaptInv, color = adaptInv)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values = c("orange", "black", "grey21")) +
  scale_fill_manual(values = color_scale) + 
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        title = element_text(size = 17))+
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "All Simulations", y = expression("Length"[inv]*" (cM)"), x = " ")

plot.numQTNs.pgen.onePanel <- ggplot(data = plot.output[plot.output$enVar == 0.1  &
                                                          plot.output$muProp == 0.1 &
                                                          plot.output$alpha == 0.002,], 
                                    aes(x = adaptInv, y= num_qtns_Lscaled, fill = adaptInv, color = adaptInv)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values = c("orange", "black", "grey21")) +
  scale_fill_manual(values = color_scale) + 
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        title = element_text(size = 17))+
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Highly Polygenic Architecture", y = " ", x = " ")

plot.numQTNs.ogen.onePanel <- ggplot(data = plot.output[plot.output$enVar == 0.1  &
                                                          plot.output$muProp == 0.01 &
                                                          plot.output$alpha == 0.2,], 
                                     aes(x = adaptInv, y= num_qtns_Lscaled, fill = adaptInv, color = adaptInv)) +
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values = c("orange", "black", "grey21")) +
  scale_fill_manual(values = color_scale) + 
  guides(fill = guide_legend(title = "Inversion Status"),
         color = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        title = element_text(size = 17))+
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "Polygenic Architecture", y = expression(bar(N)[QTNs] / "(Length"[inv]*")"), x = " ")

pdf(file = paste0(folderOutFig, "Fig3_characteristics_envar.pdf"), width = 9, height = 8)
ggarrange(plot.age.onePanel, plot.length.onePanel, blank, plot.numQTNs.ogen.onePanel, plot.numQTNs.pgen.onePanel, invChar.legend, 
          ncol = 3, nrow = 2,
          widths = c(2.3, 2.3, 1, 2.3, 2.3, 1), 
          labels = c("A", "B"," ", "C", "D", ""),
          font.label = list(size=20))
dev.off()

### end Q4 plotting
######################################################################################################


######################################################################################################
#### Q5: Plot Genome scan results

df.outliers <- df.summary 
df.outlierSumData <- left_join(df.outliers, df.simStats[,1:(ncol(df.simStats)-1)], by = "seed")
df.muInv3_outliers <- df.outlierSumData[df.outlierSumData$muInv == 1e-03,]

df.pcadapt.av <- aggregate(cbind(true_pos_pcadapt, false_neg_pcadapt, 
                                 true_neg_pcadapt, false_pos_pcadapt,
                                 true_neg_pcadapt_NS, false_pos_pcadapt_NS)~muProp + sigmaK + alpha + enVar + mig1 + mig2, 
                           FUN=mean, data = df.muInv3_outliers)
df.pcadapt.sd <- aggregate(cbind(true_pos_pcadapt, false_neg_pcadapt, 
                                 true_neg_pcadapt, false_pos_pcadapt,
                                 true_neg_pcadapt_NS, false_pos_pcadapt_NS)~muProp + sigmaK + alpha + enVar + mig1 + mig2, 
                           FUN=sd, data = df.muInv3_outliers)
df.outflank.av <- aggregate(cbind(true_pos_outflank, false_neg_outflank, 
                                 true_neg_outflank, false_pos_outflank,
                                 true_neg_outflank_NS, false_pos_outflank_NS)~muProp + sigmaK + alpha + enVar + mig1 + mig2, 
                           FUN=mean, data = df.muInv3_outliers)
df.outflank.sd <- aggregate(cbind(true_pos_outflank, false_neg_outflank, 
                                  true_neg_outflank, false_pos_outflank,
                                  true_neg_outflank_NS, false_pos_outflank_NS)~muProp + sigmaK + alpha + enVar + mig1 + mig2, 
                            FUN=sd, data = df.muInv3_outliers)

df.pcadapt.long <- as.data.frame(pivot_longer(df.pcadapt.av, cols = c(true_pos_pcadapt, false_neg_pcadapt, 
                                                                      true_neg_pcadapt, false_pos_pcadapt, true_neg_pcadapt_NS, false_pos_pcadapt_NS),
                                              names_to = "outcome", values_to = "count"))
df.pcadapt.sd.long <- as.data.frame(pivot_longer(df.pcadapt.sd, cols =  c(true_pos_pcadapt, false_neg_pcadapt, 
                                                                          true_neg_pcadapt, false_pos_pcadapt, true_neg_pcadapt_NS, false_pos_pcadapt_NS),
                                                 names_to = "outcome", values_to = "sd"))
df.outflank.long <- as.data.frame(pivot_longer(df.outflank.av, cols = c(true_pos_outflank, false_neg_outflank, 
                                                                        true_neg_outflank, false_pos_outflank, true_neg_outflank_NS, false_pos_outflank_NS),
                                              names_to = "outcome", values_to = "count"))
df.outflank.sd.long <- as.data.frame(pivot_longer(df.outflank.sd, cols =  c(true_pos_outflank, false_neg_outflank, 
                                                                            true_neg_outflank, false_pos_outflank, true_neg_outflank_NS, false_pos_outflank_NS),
                                                 names_to = "outcome", values_to = "sd"))

for(i in 1:6 ){
  df.pcadapt.long[,i] <- as.factor(df.pcadapt.long[,i])
  df.pcadapt.sd.long[,i] <- as.factor(df.pcadapt.sd.long[,i])
  df.outflank.long[,i] <- as.factor(df.outflank.long[,i])
  df.outflank.sd.long[,i] <- as.factor(df.outflank.sd.long[,i])
}

df.pcadapt.long.plot <- full_join(df.pcadapt.long, df.pcadapt.sd.long, by = c("muProp", "sigmaK", "alpha", "enVar",  "mig1",  "mig2", "outcome"))
df.pcadapt.long.plot$upSD <- df.pcadapt.long.plot$count + df.pcadapt.long.plot$sd
df.pcadapt.long.plot$lowSD <-  df.pcadapt.long.plot$count - df.pcadapt.long.plot$sd
df.pcadapt.long.plot$lowSD[df.pcadapt.long.plot$lowSD<0] <- 0
df.outflank.long.plot <- full_join(df.outflank.long, df.outflank.sd.long, by = c("muProp", "sigmaK", "alpha", "enVar",  "mig1",  "mig2", "outcome"))
df.outflank.long.plot$upSD <- df.outflank.long.plot$count + df.outflank.long.plot$sd
df.outflank.long.plot$lowSD <-  df.outflank.long.plot$count - df.outflank.long.plot$sd
df.outflank.long.plot$lowSD[df.outflank.long.plot$lowSD<0] <- 0

df.pcadapt.long.plot$sigmaK <- factor(df.pcadapt.long.plot$sigmaK, c(3, 1.5, 0.75))
df.outflank.long.plot$sigmaK <- factor(df.outflank.long.plot$sigmaK, c(3, 1.5, 0.75))

df.pcadapt.long.plot$outcome <- as.factor(df.pcadapt.long.plot$outcome)
df.outflank.long.plot$outcome <- as.factor(df.outflank.long.plot$outcome)
df.pcadapt.long.Selplot <- df.pcadapt.long.plot[df.pcadapt.long.plot$outcome == 'true_pos_pcadapt' | df.pcadapt.long.plot$outcome == 'true_neg_pcadapt' |
                                                  df.pcadapt.long.plot$outcome == 'false_pos_pcadapt' | df.pcadapt.long.plot$outcome == 'false_neg_pcadapt',]
df.pcadapt.long.Selplot$outcome <- recode_factor(df.pcadapt.long.Selplot$outcome, 'true_pos_pcadapt' = 'Adaptive Outlier', 'true_neg_pcadapt' = 'Nonadaptive Nonoutlier', 
                                         'false_pos_pcadapt' = 'Nonadaptive Outlier', 'false_neg_pcadapt' = 'Adaptive Nonoutlier')
df.outflank.long.Selplot <- df.outflank.long.plot[df.outflank.long.plot$outcome == 'true_pos_outflank' | df.outflank.long.plot$outcome == 'true_neg_outflank' |
                                                  df.outflank.long.plot$outcome == 'false_pos_outflank' | df.outflank.long.plot$outcome == 'false_neg_outflank',]
df.outflank.long.Selplot$outcome <- as.factor(df.outflank.long.Selplot$outcome)
df.outflank.long.Selplot$outcome <- recode_factor(df.outflank.long.Selplot$outcome, 'true_pos_outflank' = 'Adaptive Outlier', 'true_neg_outflank' = 'Nonadaptive Nonoutlier', 
                                         'false_pos_outflank' = 'Nonadaptive Outlier', 'false_neg_outflank' = 'Adaptive Nonoutlier')

df.pcadapt.pgen <- df.pcadapt.long.Selplot[df.pcadapt.long.Selplot$enVar == 0.1 & df.pcadapt.long.Selplot$alpha == "0.002" & df.pcadapt.long.Selplot$muProp == 0.1,]
plot.pcadapt.pgen.stack <- ggplot(df.pcadapt.pgen[order(df.pcadapt.pgen$outcome),], 
                            aes(x = mig1, y = count, fill = outcome)) + 
  geom_bar(position="stack", stat="identity") + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = c("cornflowerblue", "navy", "lightgoldenrod1", "goldenrod")) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(strip.text = element_text(size=14),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=18),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = "Highly Polygenic Architecture\\nNumber of inversions", x = " ") + 
  ylim(c(0, 50)) + 
  scale_x_discrete("Migration Rate", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)


df.pcadapt.ogen <- df.pcadapt.long.Selplot[df.pcadapt.long.Selplot$enVar == 0.1 & df.pcadapt.long.Selplot$alpha == "0.2" & df.pcadapt.long.Selplot$muProp == 0.01,]
plot.pcadapt.ogen.stack <- ggplot(df.pcadapt.ogen[order(df.pcadapt.ogen$outcome),], 
                            aes(x = mig1, y = count, fill = outcome)) + 
  geom_bar(position = "stack", stat="identity") + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = c("cornflowerblue", "navy", "lightgoldenrod1", "goldenrod")) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(strip.text = element_text(size=14),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=18),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "pcadapt", y = "Polygenic Architecture\\nNumber of inversions", x = " ") + 
  ylim(c(0, 50)) +
  scale_x_discrete(" ", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)

df.outflank.pgen <- df.outflank.long.Selplot[df.outflank.long.Selplot$enVar == 0.1 & df.outflank.long.Selplot$alpha == "0.002" & df.outflank.long.Selplot$muProp == 0.1,]
plot.outflank.pgen.stack <- ggplot(df.outflank.pgen[order(df.outflank.pgen$outcome),], 
                             aes(x = mig1, y = count, fill = outcome)) + 
  geom_bar(position = "stack", stat="identity") + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = c("cornflowerblue", "navy", "lightgoldenrod1", "goldenrod")) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(strip.text = element_text(size=14),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=18),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = "Migration Rate") + 
  ylim(c(0, 50)) +
  scale_x_discrete("Migration Rate", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)


df.outflank.ogen <- df.outflank.long.Selplot[df.outflank.long.Selplot$enVar == 0.1 & df.outflank.long.Selplot$alpha == "0.2" & df.outflank.long.Selplot$muProp == 0.01,]
plot.outflank.ogen.stack <- ggplot(df.outflank.ogen[order(df.outflank.ogen$outcome),], 
                             aes(x = mig1, y = count, fill = outcome)) + 
  geom_bar(position = "stack", stat="identity") + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = c("cornflowerblue", "navy", "lightgoldenrod1", "goldenrod")) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(strip.text = element_text(size=14),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=18),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "OutFLANK", y = " ", x = " ") + 
  ylim(c(0, 50)) +
  scale_x_discrete(" ", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE)


eval.legend.stack <- g_legend(plot.pcadapt.pgen.stack)
plot.pcadapt.pgen.stack.noLeg <- plot.pcadapt.pgen.stack + theme(legend.position = "none")


pdf(file = paste0(folderOutFig, "SFigX_outliers_count_envar.pdf"), width = 15, height = 12)
ggarrange(plot.pcadapt.ogen.stack, plot.outflank.ogen.stack, blank, 
          plot.pcadapt.pgen.stack.noLeg, plot.outflank.pgen.stack,
          eval.legend.stack, ncol = 3, nrow = 2, widths = c(2.3,2.3,0.8,2.3,2.3,0.8),
          labels = c("A", "B", "", "C", "D"),
          font.label = list(size = 30))
dev.off()

## No Selection Sims

df.pcadapt.long.noSelplot <- df.pcadapt.long.plot[df.pcadapt.long.plot$outcome == 'true_neg_pcadapt_NS' | df.pcadapt.long.plot$outcome == 'false_pos_pcadapt_NS',]
df.pcadapt.long.noSelplot$outcome <- recode_factor(df.pcadapt.long.noSelplot$outcome, 'true_neg_pcadapt_NS' = 'Nonadaptive Nonoutlier', 
                                                   'false_pos_pcadapt_NS' = 'Nonadaptive Outlier')
df.outflank.long.noSelplot <- df.outflank.long.plot[df.outflank.long.plot$outcome == 'true_neg_outflank_NS' | df.outflank.long.plot$outcome == 'false_pos_outflank_NS',]
df.outflank.long.noSelplot$outcome <- recode_factor(df.outflank.long.noSelplot$outcome,'true_neg_outflank_NS' = 'Nonadaptive Nonoutlier', 
                                                    'false_pos_outflank_NS' = 'Nonadaptive Outlier')

df.pcadapt.long.noSelplot$outcome <- factor(df.pcadapt.long.noSelplot$outcome, levels = c('Nonadaptive Nonoutlier', 'Nonadaptive Outlier'))
df.pcadapt.NS.pgen <- df.pcadapt.long.noSelplot[df.pcadapt.long.noSelplot$enVar == 0.1 & df.pcadapt.long.noSelplot$alpha == "0.002" & df.pcadapt.long.noSelplot$muProp == 0.1,]

plot.pcadapt.NS.pgen.stack <- ggplot(df.pcadapt.NS.pgen[order(df.pcadapt.NS.pgen$outcome),], 
                               aes(x = mig1, y = count, fill = outcome)) + 
  geom_bar(position = "stack", stat="identity") + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = c("navy", "lightgoldenrod1")) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(strip.text = element_text(size=14),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=18),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = "Highly Polygenic Architecture\\nNumber of inversions", x = "Migration Rate") +
  scale_x_discrete("Migration Rate", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE) +
  ylim(0,30)

df.pcadapt.NS.ogen <- df.pcadapt.long.noSelplot[df.pcadapt.long.noSelplot$enVar == 0.1 & df.pcadapt.long.noSelplot$alpha == "0.2" & df.pcadapt.long.noSelplot$muProp == 0.01,]
plot.pcadapt.NS.ogen.stack <- ggplot(df.pcadapt.NS.ogen[order(df.pcadapt.NS.ogen$outcome),], 
                               aes(x = mig1, y = count, fill = outcome)) + 
  geom_bar(position = "stack", stat="identity") + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = c("navy", "lightgoldenrod1")) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(strip.text = element_text(size=14),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=18),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "pcadapt", y = "Polygenic Architecture\\nNumber of inversions", x = " ") +
  scale_x_discrete(" ", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE) +
  ylim(0,30)


df.outflank.long.noSelplot$outcome <- factor(df.outflank.long.noSelplot$outcome, levels = c('Nonadaptive Nonoutlier', 'Nonadaptive Outlier'))

df.outflank.NS.pgen <- df.outflank.long.noSelplot[df.outflank.long.noSelplot$enVar == 0.1 & df.outflank.long.noSelplot$alpha == "0.002" & df.outflank.long.noSelplot$muProp == 0.1,]
plot.outflank.NS.pgen.stack <- ggplot(df.outflank.NS.pgen[order(df.outflank.NS.pgen$outcome),], 
                                aes(x = mig1, y = count, fill = outcome)) + 
  geom_bar(position = "stack", stat="identity") + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = c("navy", "lightgoldenrod1")) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(strip.text = element_text(size=14),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=18),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = " ", y = " ", x = "Migration Rate") +
  scale_x_discrete("Migration Rate", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE) +
  ylim(0,30)

df.outflank.NS.ogen <- df.outflank.long.noSelplot[df.outflank.long.noSelplot$enVar == 0.1 & df.outflank.long.noSelplot$alpha == "0.2" & df.outflank.long.noSelplot$muProp == 0.01,]
plot.outflank.NS.ogen.stack <- ggplot(df.outflank.NS.ogen[order(df.outflank.NS.ogen$outcome),], 
                                aes(x = mig1, y = count, fill = outcome)) + 
  geom_bar(position = "stack", stat="identity") + 
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels))+ 
  scale_fill_manual(values = c("navy", "lightgoldenrod1")) + 
  guides(fill = guide_legend(title = "Inversion Status")) +
  theme_classic() +
  theme(strip.text = element_text(size=14),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        axis.title = element_text(size=18),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=15)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  labs(title = "OutFLANK", y = " ", x = " ") +
  scale_x_discrete(" ", breaks=factor(c(0.001,0.01,0.1,0.25,0.4,0.5)), drop=FALSE) +
  ylim(0,25)

genomeScans.legend.stack <- g_legend(plot.pcadapt.NS.pgen.stack)
plot.pcadapt.NS.pgen.stack.noLeg <- plot.pcadapt.NS.pgen.stack + theme(legend.position = "none")

pdf(file = paste0(folderOutFig, "SFigX_outliersNS_count_envar.pdf"), width = 15, height = 12)
ggarrange(plot.pcadapt.NS.ogen.stack, plot.outflank.NS.ogen.stack, blank, plot.pcadapt.NS.pgen.stack.noLeg, plot.outflank.NS.pgen.stack, 
          genomeScans.legend.stack, ncol = 3, nrow = 2, widths = c(2.3,2.3,0.8,2.3,2.3,0.8), 
          labels = c("A", "B", "", "C", "D"),
          font.label = list(size = 30))
dev.off()

#### end Q4 plotting
######################################################################################################

######################################################################################################
#### average percent VA ####

df.muInv3_qtn <- aggregate(cbind(av.effect, av.perc.VA)~muProp + sigmaK + alpha + enVar + mig1, 
                            FUN=mean, data = df.muInv3)

mean(df.muInv3[df.muInv3$muProp == 0.1 & df.muInv3$alpha == 0.002,]$av.perc.VA)
min(df.muInv3_qtn[df.muInv3$muProp == 0.1 & df.muInv3$alpha == 0.002,]$av.perc.VA)

mean(df.muInv3_0_qtn[df.muInv3_0_qtn$muProp == 0.0000001 & df.muInv3_0_qtn$alpha == 0.2,]$av_qtn_perc_VA)
min(df.muInv3_0_qtn[df.muInv3_0_qtn$muProp == 0.0000001 & df.muInv3_0_qtn$alpha == 0.2,]$av_qtn_perc_VA)

mean(abs(df.muInv3[df.muInv3$muProp == 0.01 & df.muInv3$alpha == 0.2 & df.muInv3$enVar == 0.1,]$av.perc.VA))
mean(abs(df.muInv3[df.muInv3$muProp == 0.1 & df.muInv3$alpha == 0.002 & df.muInv3$enVar == 0.1,]$av.perc.VA))
min(abs(df.muInv3[df.muInv3$muProp == 0.1 & df.muInv3$alpha == 0.002,]$av.perc.VA))

max(abs(df.muInv3[df.muInv3$muProp == 0.01 & df.muInv3$alpha == 0.2,]$av.effect))
min(abs(df.muInv3[df.muInv3$muProp == 0.01 & df.muInv3$alpha == 0.2,]$av.effect))


max(df.summary$av_qtn_effect_size)

#### end average percent VA
######################################################################################################


#######################################################################################################
#### Overlapping inversions

plot.overlap_pgen <- ggplot(df.muInv3[df.muInv3$enVar == 0.1 & 
                                        df.muInv3$alpha == 0.002 &
                                        df.muInv3$muProp == 0.1,], 
                            aes(x = mig1, y = num.adapt.overlap.invs),color =  viridis(4)[2], fill = viridis(4)[2]) + 
  geom_boxplot(color =  'black', fill = viridis(4)[2]) +
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  labs(title = "Highly Polygenic Architecture",
       y = " ",
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(c(0,9))


plot.overlap_ogen <- ggplot(df.muInv3[df.muInv3$enVar == 0.1 & 
                                        df.muInv3$alpha == 0.2 &
                                        df.muInv3$muProp == 0.01,], 
                            aes(x = mig1, y = num.adapt.overlap.invs),color =  viridis(4)[2], fill = viridis(4)[2]) + 
  geom_boxplot(color =  'black', fill = viridis(4)[4]) +
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  labs(title = "Polygenic Architecture",
       y = "Number of Overlapping Inversions",
       x = " ") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(c(0,9))


plot.within_pgen <- ggplot(df.muInv3[df.muInv3$enVar == 0.1 & 
                                       df.muInv3$alpha == 0.002 &
                                       df.muInv3$muProp == 0.1,], 
                           aes(x = mig1, y = num.adapt.within.invs),color =  viridis(4)[2], fill = viridis(4)[2]) + 
  geom_boxplot(color =  'black', fill = viridis(4)[2]) +
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  labs(title = "Highly Polygenic Architecture",
       y = " ",
       x = "Migration Rate") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
 ylim(c(0,9))



plot.within_ogen <- ggplot(df.muInv3[df.muInv3$enVar == 0.1 & 
                                       df.muInv3$alpha == 0.2 &
                                       df.muInv3$muProp == 0.01,], 
                             aes(x = mig1, y = num.adapt.within.invs)) + 
  geom_boxplot(color =  'black', fill = viridis(4)[4]) +
  facet_wrap(~sigmaK, labeller = labeller(sigmaK = sigmaK.labels)) + 
  labs(title = "Polygenic Architecture",
       y = "Number of Inversions-Within-Inversions",
       x = "Migration Rate") +
  guides(color = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)"),
         shape = guide_legend(title = "QTN Mutation Rate\\n(Ne*mu)")) +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(c(0,9))

pdf(file = paste0(folderOutFig, "SFig7_overlappingInversions.pdf"), width = 15, height = 12)
ggarrange(plot.overlap_ogen, plot.overlap_pgen, plot.within_ogen, plot.within_pgen, ncol = 2, nrow = 2,
          widths = c(2.3,2.3,2.3,2.3), 
          labels = c("A", "B", "C", "D"))
dev.off()


#### end overlapping inversions
######################################################################################################

