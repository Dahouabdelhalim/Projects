# set directories
DATADIR <- "/scratch/schaal.s/InversionSimulations/results/20220419/"
PARAMSDIR <- "/scratch/schaal.s/InversionSimulations/src/"
DATAOUT <- "/scratch/schaal.s/InversionSimulations/figures/20220429/"

# install packages
packages_needed <- c("tidyverse", "viridisLite", "ggpubr")
for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library( packages_needed[i], character.only = TRUE)
}

# load data
df.params <- read.table(paste0(PARAMSDIR, "FullSet_dfparams.txt"))

# loop through pop dynamics files and combine them
df.muts <- NULL
count <- 0
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  mutsNewFile <- read.table(paste(DATADIR, seed, "_outputMutations.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
  if(nrow(mutsNewFile) > 0){
    mutsNewFile$seed <- seed
    df.muts <-  rbind(df.muts, mutsNewFile)
  }
  count <- count + 1
  print(count)
}
write.table(df.muts, paste0(DATAOUT,"FullSet_muts.txt"), row.names = FALSE)

# df.muts <- read.table(paste0(DATAOUT, "FullSet_muts.txt"), header =TRUE)
df.muts_params <- left_join(df.muts, df.params, by ="seed")

######################


qtns <- df.muts_params[df.muts_params$type == "m2",]
df.muts_plotting <- NULL
for(i in 1:length(unique(qtns$seed))){
  subset.df <- qtns[qtns$seed == unique(qtns$seed)[i],]
  subset.df$perc_VA <- ((subset.df$selCoef^2*subset.df$freq*(1-subset.df$freq))/sum(subset.df$selCoef^2*subset.df$freq*(1-subset.df$freq)))*100
  df.muts_plotting <- rbind(df.muts_plotting, subset.df)
}

#av.effect <- mean(qtns$selCoef)
for(i in 11:20){
  df.muts_plotting[,i] <- as.factor(df.muts_plotting[,i])
}
df.muts_plotting$sigmaK <- factor(df.muts_plotting$sigmaK, c(3, 1.5, 0.75))
df.muts_plotting$sigmaK <- recode_factor(df.muts_plotting$sigmaK, '3' = 'Weak Selection', '1.5' = 'Moderate Selection','0.75' = 'Strong Selection')


plot.qtns_pgen <- ggplot(df.muts_plotting[df.muts_plotting$enVar == 0.1 & 
                                        df.muts_plotting$alpha == 0.002 &
                                        df.muts_plotting$muProp == 0.1 &
                                        df.muts_plotting$muInv == 0.001,], 
                            aes(x = mig1, y = selCoef),color =  viridis(4)[2], fill = viridis(4)[2]) + 
  geom_boxplot(color =  'black', fill = viridis(4)[2]) +
  facet_wrap(~sigmaK) + 
  labs(title = "Highly Polygenic Architecture",
       y = " ",
       x = " ") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"))
  #ylim(c(0,6))


plot.qtns_ogen <- ggplot(df.muts_plotting[df.muts_plotting$enVar == 0.1 & 
                                        df.muts_plotting$alpha == 0.2 &
                                        df.muts_plotting$muProp == 0.01 &
                                        df.muts_plotting$muInv == 0.001,], 
                            aes(x = mig1, y = selCoef),color =  viridis(4)[4], fill = viridis(4)[4]) + 
  #geom_errorbar(aes(ymin=num.adapt.overlap.invs_lowSD, ymax=num.adapt.overlap.invs_upSD), size = 0.2, width=.5) +
  geom_boxplot(color =  'black', fill = viridis(4)[4]) +
  #geom_point(fill = viridis(4)[2], color = "navy", shape = 21, size = 3) + 
  facet_wrap(~sigmaK) + 
  labs(title = "Polygenic Architecture",
       y = "QTN Effect Sizes",
       x = " ") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) 
 # ylim(c(0,6))


plot.qtns_pgen_noInv <- ggplot(df.muts_plotting[df.muts_plotting$enVar == 0.1 & 
                                       df.muts_plotting$alpha == 0.002 &
                                       df.muts_plotting$muProp == 0.1 &
                                       df.muts_plotting$muInv == 0,], 
                           aes(x = mig1, y = selCoef),color =  viridis(4)[2], fill = viridis(4)[2]) + 
  geom_boxplot(color =  'black', fill = viridis(4)[2]) +
  facet_wrap(~sigmaK) + 
  labs(title = "Highly Polygenic Architecture",
       y = " ",
       x = " ") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"))
  #ylim(c(0,5))



plot.qtns_ogen_noInv <- ggplot(df.muts_plotting[df.muts_plotting$enVar == 0.1 & 
                                       df.muts_plotting$alpha == 0.2 &
                                       df.muts_plotting$muProp == 0.01  &
                                       df.muts_plotting$muInv == 0,], 
                           aes(x = mig1, y = selCoef), color =  viridis(4)[2], fill = viridis(4)[2]) + 
  geom_boxplot(color =  'black', fill = viridis(4)[4]) +
  facet_wrap(~sigmaK) + 
  labs(title = "Polygenic Architecture",
       y = "QTN Effect Sizes",
       x = " ") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) 
  #ylim(c(0,5))
max(df.muts_plotting[df.muts_plotting$enVar == 0.1 & 
                   df.muts_plotting$alpha == 0.002 &
                   df.muts_plotting$muProp == 0.1 &
                   df.muts_plotting$muInv == 0.001,]$perc_VA)

plot.VAperc_pgen <- ggplot(df.muts_plotting[df.muts_plotting$enVar == 0.1 & 
                                            df.muts_plotting$alpha == 0.002 &
                                            df.muts_plotting$muProp == 0.1 &
                                            df.muts_plotting$muInv == 0.001,], 
                         aes(x = mig1, y = perc_VA),color =  viridis(4)[2], fill = viridis(4)[2]) + 
  geom_boxplot(color =  'black', fill = viridis(4)[2]) +
  facet_wrap(~sigmaK) + 
  labs(title = "Highly Polygenic Architecture",
       y = " ",
       x = "Migration Rate") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"))
#ylim(c(0,6))


plot.VAperc_ogen <- ggplot(df.muts_plotting[df.muts_plotting$enVar == 0.1 & 
                                            df.muts_plotting$alpha == 0.2 &
                                            df.muts_plotting$muProp == 0.01 &
                                            df.muts_plotting$muInv == 0.001,], 
                         aes(x = mig1, y = perc_VA),color =  viridis(4)[4], fill = viridis(4)[4]) + 
  #geom_errorbar(aes(ymin=num.adapt.overlap.invs_lowSD, ymax=num.adapt.overlap.invs_upSD), size = 0.2, width=.5) +
  geom_boxplot(color =  'black', fill = viridis(4)[4]) +
  #geom_point(fill = viridis(4)[2], color = "navy", shape = 21, size = 3) + 
  facet_wrap(~sigmaK) + 
  labs(title = "Polygenic Architecture",
       y = expression(bold("QTN %V"[A])),
       x = "Migration Rate") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) 
# ylim(c(0,6))


plot.VAperc_pgen_noInv <- ggplot(df.muts_plotting[df.muts_plotting$enVar == 0.1 & 
                                                  df.muts_plotting$alpha == 0.002 &
                                                  df.muts_plotting$muProp == 0.1 &
                                                  df.muts_plotting$muInv == 0,], 
                               aes(x = mig1, y = perc_VA),color =  viridis(4)[2], fill = viridis(4)[2]) + 
  geom_boxplot(color =  'black', fill = viridis(4)[2]) +
  facet_wrap(~sigmaK) + 
  labs(title = "Highly Polygenic Architecture",
       y = " ",
       x = "Migration Rate") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92"))
#ylim(c(0,5))



plot.VAperc_ogen_noInv <- ggplot(df.muts_plotting[df.muts_plotting$enVar == 0.1 & 
                                                  df.muts_plotting$alpha == 0.2 &
                                                  df.muts_plotting$muProp == 0.01  &
                                                  df.muts_plotting$muInv == 0,], 
                               aes(x = mig1, y = perc_VA), color =  viridis(4)[2], fill = viridis(4)[2]) + 
  geom_boxplot(color =  'black', fill = viridis(4)[4]) +
  facet_wrap(~sigmaK) + 
  labs(title = "Polygenic Architecture",
       y = expression(bold("QTN %V"[A])),
       x = "Migration Rate") +
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title = element_text(size=18)) +
  theme(panel.background = element_blank(), 
        strip.background = element_rect(colour = "white", fill = "grey92")) 


pdf(file = paste0(DATAOUT, "SFig3_qtnEffect_VA.pdf"), width = 15, height = 12)
ggarrange(plot.qtns_ogen, plot.qtns_pgen, plot.VAperc_ogen, plot.VAperc_pgen, ncol = 2, nrow = 2,
          widths = c(2.3,2.3,2.3,2.3), 
          labels = c("A", "B", "C", "D"))
dev.off()


pdf(file = paste0(DATAOUT, "SFig3_qtnEffect_VA_noInv.pdf"), width = 15, height = 12)
ggarrange(plot.qtns_ogen_noInv, plot.qtns_pgen_noInv, plot.VAperc_ogen_noInv, plot.VAperc_pgen_noInv, ncol = 2, nrow = 2,
          widths = c(2.3,2.3,2.3,2.3), 
          labels = c("A", "B", "C", "D"))
dev.off()


