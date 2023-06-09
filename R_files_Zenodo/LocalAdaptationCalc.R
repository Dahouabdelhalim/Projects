# set directories
DATADIR <- "/scratch/schaal.s/InversionSimulations/results/20220220/10gen/"
PARAMSDIR <- "/scratch/schaal.s/InversionSimulations/src/"
DATAOUT <- "/scratch/schaal.s/InversionSimulations/figures/"

# install packages
packages_needed <- c("tidyverse")
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

# load data
df.params <- read.table(paste0(PARAMSDIR, "FullSet_dfparams.txt"))

# loop through pop dynamics files and combine them
df.popDyn <- NULL
count <- 0
for(i in 1:nrow(df.params)){
  seed <- df.params$seed[i]
  popDynNewFile <- read.table(paste(DATADIR, seed, "_outputPopDynam.txt", sep=""), header = TRUE, stringsAsFactors = FALSE)
  if(nrow(popDynNewFile) > 0){
    popDynNewFile$seed <- seed
    df.popDyn <-  rbind(df.popDyn, popDynNewFile)
  }
  count <- count + 1
  print(count)
}
write.table(df.popDyn, paste0(DATAOUT,"FullSet_popDyn.txt"), row.names = FALSE)

df.popDyn_stats <- left_join(df.popDyn, df.params, by ="seed")

# Average data across parameter combinations 
popDyn_plotting_df <- aggregate(cbind(localAdaptSA, invVA, outVA)~sim_gen +muProp + muInv + sigmaK + alpha + rep + mig1 + enVar, FUN = mean, data = df.popDyn_stats)

# Write out the averages dataframe
write.table(popDyn_plotting_df, paste0(DATAOUT, "FullSet_popDynam.txt"))

#popDyn_plotting_df <-read.table(paste0(folder, "FullSet_popDynam.txt"))

########################################



########################################
### Manuipulate and plot data

for(i in 2:(ncol(popDyn_plotting_df)-3)){
  popDyn_plotting_df[,i] <- as.factor(popDyn_plotting_df[,i])
}
popDyn_plotting_df$sigmaK <- factor(popDyn_plotting_df$sigmaK, c(3, 1.5, 0.75))
popDyn_plotting_df$sigmaK <- recode_factor(popDyn_plotting_df$sigmaK, '3' = 'Weak Selection', '1.5' = 'Moderate Selection','0.75' = 'Strong Selection')
popDyn_plotting_df$muInv <- recode_factor(popDyn_plotting_df$muInv, '0.001' = 'Inversions', '0' = 'No Inversions')



LA.highpoly.envar <- ggplot(data = popDyn_plotting_df[popDyn_plotting_df$muProp ==0.1 &
                                   popDyn_plotting_df$alpha ==0.002 &
                                   popDyn_plotting_df$enVar == 0.1,],
       aes(x = sim_gen, y = localAdaptSA, group = interaction(mig1, rep, muInv, sigmaK))) +
  geom_line(aes(color = mig1, linetype = rep)) +
  scale_color_manual(values = viridis(6)) +
  geom_vline(xintercept =  10000, lty = "dashed") + 
  facet_grid(sigmaK~muInv) +
  labs(title = "Highly Polygenic - Environmental Variance 0.1", x = "Generation", y = " ") +
  guides(color = guide_legend(title = "Migration Rate"),
         linetype = guide_legend(title = "Replicate")) +
  theme_classic() +
  theme(strip.text = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        title = element_text(size = 14))+
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(-0.05, 1)

LA.poly.envar <- ggplot(data = popDyn_plotting_df[popDyn_plotting_df$muProp ==0.01 &
                                   popDyn_plotting_df$alpha ==0.2 &
                                   popDyn_plotting_df$enVar == 0.1,],
       aes(x = sim_gen, y = localAdaptSA, group = interaction(mig1, rep, muInv, sigmaK))) +
  geom_line(aes(color = mig1, linetype = rep)) +
  scale_color_manual(values = viridis(6)) +
  geom_vline(xintercept =  10000, lty = "dashed") + 
  facet_grid(sigmaK~muInv) +
  labs(title = "Polygenic - Environmental Variance 0.1", x = "Generation", y = "Local Adaptation") +
  guides(color = guide_legend(title = "Migration Rate"),
         linetype = guide_legend(title = "Replicate")) +
  theme_classic() +
  theme(strip.text = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        title = element_text(size = 14)) +
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  theme(legend.position = "none") +
  ylim(-0.05, 1)


LA.highpoly.Noenvar <- ggplot(data = popDyn_plotting_df[popDyn_plotting_df$muProp ==0.1 &
                                                        popDyn_plotting_df$alpha ==0.002 &
                                                        popDyn_plotting_df$enVar == 0,],
                            aes(x = sim_gen, y = localAdaptSA, group = interaction(mig1, rep, muInv, sigmaK))) +
  geom_line(aes(color = mig1, linetype = rep)) +
  geom_vline(xintercept =  10000, lty = "dashed") + 
  scale_color_manual(values = viridis(6)) +
  #facet_wrap(~sigmaK + muInv, nrow = 3) +
  facet_grid(sigmaK~muInv) +
  labs(title = "Highly Polygenic - No Environmental Variance", x = "Generation", y = " ") +
  guides(color = guide_legend(title = "Migration Rate"),
         linetype = guide_legend(title = "Replicate")) +
  theme_classic() +
  theme(strip.text = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        title = element_text(size = 14))+
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(-0.05, 1)

LA.poly.Noenvar <- ggplot(data = popDyn_plotting_df[popDyn_plotting_df$muProp ==0.01 &
                                                    popDyn_plotting_df$alpha ==0.2 &
                                                    popDyn_plotting_df$enVar == 0,],
                        aes(x = sim_gen, y = localAdaptSA, group = interaction(mig1, rep, muInv, sigmaK))) +
  geom_line(aes(color = mig1, linetype = rep)) +
  scale_color_manual(values = viridis(6)) +
  geom_vline(xintercept =  10000, lty = "dashed") + 
  facet_grid(sigmaK~muInv) +
  labs(title = "Polygenic - No Environmental Variance", x = "Generation", y = "Local Adaptation") +
  guides(color = guide_legend(title = "Migration Rate"),
         linetype = guide_legend(title = "Replicate")) +
  theme_classic() +
  theme(strip.text = element_text(size=12),
        axis.text = element_text(size=10),
        axis.title = element_text(size=15),
        title = element_text(size = 14))  +
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  theme(legend.position = "none") +
  ylim(-0.05, 1)

LA_leg <- g_legend(LA.highpoly.envar)
LA.highpoly.envar.noLeg <- LA.highpoly.envar + theme(legend.position = "none") 
LA.highpoly.Noenvar.noLeg <- LA.highpoly.Noenvar +  theme(legend.position = "none")

folderOut <- "./figures/20220419/ManuscriptFigs/"
pdf(paste0(folderOut,"SFig4_LAthroughTime_enVar.pdf"), height = 10, width = 15)
  ggarrange(LA.poly.envar, LA.highpoly.envar.noLeg,  LA_leg, ncol = 3, widths = c(2.3,2.3,0.8),
            labels = c("A", "B", ""))
dev.off()

pdf(paste0(folderOut, "SFig5_LAthroughTime_NoenVar.pdf"), height = 10, width = 15)
  ggarrange(LA.poly.Noenvar, LA.highpoly.Noenvar.noLeg, LA_leg, ncol = 3, widths = c(2.3,2.3,0.8),
            labels = c("A", "B", ""))
dev.off()

df.VA <- pivot_longer(popDyn_plotting_df[popDyn_plotting_df$muInv == "Inversions", ], cols = c(invVA, outVA), names_to = "inOut", values_to = "VA")
df.VA$inOut <- recode_factor(df.VA$inOut, 'invVA' = 'Inversions QTN', 'outVA' = 'Collinear QTN')
df.VA$VA <- as.numeric(df.VA$VA)
VA.highpoly.envar <- ggplot(data = df.VA[df.VA$muProp ==0.1 &
                                         df.VA$alpha ==0.002 &
                                         df.VA$enVar == 0.1,],
                            aes(x = sim_gen, y = VA, group = interaction(inOut, mig1, rep, sigmaK))) +
  geom_line(aes(color = inOut, linetype = rep)) +
  scale_color_manual(values = c("firebrick", "navyblue")) +
  geom_vline(xintercept =  10000, lty = "dashed") + 
  facet_grid(sigmaK~mig1) +
  labs(title = "Highly Polygenic - Environmental Variance 0.1", x = "Generation", y = "  ") +
  guides(color = guide_legend(title = "QTN Location"),
         linetype = guide_legend(title = "Replicate"),
         size = 3) +
  theme_classic() +
  theme(strip.text = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=8),
        axis.title = element_text(size=30),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(-0.05, 100)+
  scale_x_continuous(breaks = c(20000,40000,60000))


VA.poly.envar <- ggplot(data = df.VA[df.VA$muProp ==0.01 &
                                           df.VA$alpha ==0.2 &
                                           df.VA$enVar == 0.1,],
                            aes(x = sim_gen, y = VA, group = interaction(inOut, mig1, rep, sigmaK))) +
  geom_line(aes(color = inOut, linetype = rep)) +
  scale_color_manual(values = c("firebrick", "navyblue")) +
  geom_vline(xintercept =  10000, lty = "dashed") + 
  facet_grid(sigmaK~mig1) +
  labs(title = "Polygenic - Environmental Variance 0.1", x = "Generation", y = expression("%V"["A"])) +
  guides(color = guide_legend(title = "QTN Location"),
         linetype = guide_legend(title = "Replicate")) +
  theme_classic() +
  theme(strip.text = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=8),
        axis.title = element_text(size=30),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  theme(legend.position = "none") +
  ylim(-0.05, 100)+
  scale_x_continuous(breaks = c(20000,40000,60000))


VA.highpoly.Noenvar <- ggplot(data = df.VA[df.VA$muProp ==0.1 &
                                           df.VA$alpha ==0.002 &
                                           df.VA$enVar == 0,],
                            aes(x = sim_gen, y = VA, group = interaction(inOut, mig1, rep, sigmaK))) +
  geom_line(aes(color = inOut, linetype = rep)) +
  scale_color_manual(values = c("firebrick", "navyblue")) +
  geom_vline(xintercept =  10000, lty = "dashed") + 
  facet_grid(sigmaK~mig1) +
  labs(title = "Highly Polygenic - No Environmental Variance", x = "Generation", y = "  ") +
  guides(color = guide_legend(title = "QTN Location"),
         linetype = guide_legend(title = "Replicate"),
         size = 3) +
  theme_classic() +
  theme(strip.text = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=8),
        axis.title = element_text(size=30),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  ylim(-0.05, 100)+
  scale_x_continuous(breaks = c(20000,40000,60000))


VA.poly.Noenvar <- ggplot(data = df.VA[df.VA$muProp ==0.01 &
                                       df.VA$alpha ==0.2 &
                                       df.VA$enVar == 0,],
                        aes(x = sim_gen, y = VA, group = interaction(inOut, mig1, rep, sigmaK))) +
  geom_line(aes(color = inOut, linetype = rep)) +
  scale_color_manual(values = c("firebrick", "navyblue")) +
  geom_vline(xintercept =  10000, lty = "dashed") + 
  facet_grid(sigmaK~mig1) +
  labs(title = "Polygenic - No Environmental Variance", x = "Generation", y = expression("%V"["A"])) +
  guides(color = guide_legend(title = "QTN Location"),
         linetype = guide_legend(title = "Replicate")) +
  theme_classic() +
  theme(strip.text = element_text(size=20),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=8),
        axis.title = element_text(size=30),
        title = element_text(size = 20))+
  theme(legend.text=element_text(size=20)) +
  theme(panel.background = element_blank(),
        strip.background = element_rect(colour = "white", fill = "grey92")) +
  theme(legend.position = "none") +
  ylim(-0.05, 100) +
  scale_x_continuous(breaks = c(20000,40000,60000))

VA_leg <- g_legend(VA.highpoly.envar)
VA.highpoly.envar.noLeg <- VA.highpoly.envar + theme(legend.position = "none") 
VA.highpoly.Noenvar.noLeg <- VA.highpoly.Noenvar +  theme(legend.position = "none")

pdf(paste0(folderOut,"SFig6_VAthroughTime_enVar.pdf"), height = 10, width = 27)
  ggarrange(VA.poly.envar, VA.highpoly.envar.noLeg, VA_leg, ncol = 3, widths = c(2.5,2.5,0.6),
            labels = c("A", "B", ""), font.label = list(size = 25))
dev.off()
pdf(paste0(folderOut, "SFig7_VAthroughTime_NoenVar.pdf"), height = 10, width = 27)
  ggarrange(VA.poly.Noenvar,VA.highpoly.Noenvar.noLeg, VA_leg, ncol = 3, widths = c(2.5,2.5,0.6),
            labels = c("A", "B", ""), font.label = list(size = 25))
dev.off()

