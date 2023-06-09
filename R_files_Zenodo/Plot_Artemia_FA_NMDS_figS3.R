#04 Feb 2022
#v4.1.1
#-----------


# Cleaned up for publication


### Paired with
# "01_Artemia_Neu_FA.csv"
# "02_Artemia_Phos_FA.csv"


#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(vegan)      # needed for NMDS
library(MASS)       # needed for NMDS
library(ggplot2)    # needed for plotting
library(ggrepel)    # needed for geom_text_repel
library(ggpubr)     # needed to plot panels A and B together
library(cowplot)    # needed to plot panels A and B together





setwd(dirname(rstudioapi::getActiveDocumentContext()$path))   #sets WD to folder this R file is saved in



# --------------------------------------- Neutral ----------------------------------------------


artemia_neu_bad <- read.csv("01_Artemia_Neu_FA.csv", header = FALSE)

artemia_neu <- data.frame(t(artemia_neu_bad[-1]))
colnames(artemia_neu) <- artemia_neu_bad[, 1]





#...................................................NMDS..................................................................


Art.neu.nmds <- matrix(as.numeric(as.character(unlist(artemia_neu[, 3:20]))), nrow = nrow(artemia_neu))
colnames(Art.neu.nmds) <- colnames(artemia_neu)[3:20]


#Run NMDS with a number of different axes (k)
Art.neu.mds2 <- metaMDS(Art.neu.nmds, k=2)
Art.neu.mds3 <- metaMDS(Art.neu.nmds, k=3)
Art.neu.mds4 <- metaMDS(Art.neu.nmds, k=4)
Art.neu.mds5 <- metaMDS(Art.neu.nmds, k=5)
Art.neu.mds6 <- metaMDS(Art.neu.nmds, k=6)

#Calculate how much stress changes with addition of axes
##Rule of thumb, only keep axes that decrease stress by at least ~0.05
Art.neu.mds2$stress - Art.neu.mds3$stress
Art.neu.mds3$stress - Art.neu.mds4$stress

#using 3 axes here
plot(Art.neu.mds3)
Art.neu.mds3$species  #ZP scores
Art.neu.mds3$points  #Individual scores
Art.neu.mds3$stress   #Stress value = 0.022
stressplot(Art.neu.mds3)  #Diagnostic - pretty good




scrs.neu <- as.data.frame(scores(Art.neu.mds3, display = "sites")) #values for each axis
scrs.neu <- cbind(scrs.neu, Group = artemia_neu$Treatment)  #Adding the associated treatment group names





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Axes 1 & 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Only fitting vectors for the FA most driving differences between groups (from SIMPER analysis) to minimize clutter on figure, but only for axes 1 and 2
vec.sp.neu <- envfit(Art.neu.mds3 ~ EPA + DHA + C18.1n9 + C18.3 + C18.2 + C24.0, as.data.frame(Art.neu.nmds))
spp.scrs.neu <- as.data.frame(scores(vec.sp.neu, display = "vectors"))
spp.scrs.neu <- cbind(spp.scrs.neu, Species = c("EPA", "DHA", "C18:1n-9", "C18:3", "C18:2", "C24:0"))



NMDS.treatment.colors <- c("black", "dodgerblue4", "darkcyan", "gold")


x11()
p1 <- ggplot(scrs.neu, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = factor(Group)), size = 2) +
  coord_fixed() +
  geom_segment(data = spp.scrs.neu,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour="grey",
               size = 1.0) +
  geom_text_repel(data = spp.scrs.neu, aes(x = NMDS1, y = NMDS2, label = Species), size = 3,
            show.legend = FALSE, color = "black")  +
  ylim(-1,1) + #xlim(-2,3) + 
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.title = element_text(size = 12), legend.text = element_text(color = NMDS.treatment.colors, size = 12),
        legend.key.size = unit(2, "lines"))
neu.p1 <- p1 +  scale_color_manual(values = NMDS.treatment.colors,
                                  labels = c("Control", "PUFA", "PUFA + E", "Unenriched"),
                                  name = "") +
  guides(color = guide_legend(override.aes = list(size = 2.0)))
neu.p1



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Axes 1 & 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Only fitting vectors for the FA most driving differences between groups (from SIMPER analysis) to minimize clutter on figure, but only for axes 1 and 3
vec.sp.neu <- envfit(Art.neu.mds3 ~ EPA + DHA + C18.1n9 + C18.3 + C18.2 + C24.0, as.data.frame(Art.neu.nmds), perm=1000, choices = c(1,3))
spp.scrs.neu <- as.data.frame(scores(vec.sp.neu, display = "vectors"))
spp.scrs.neu <- cbind(spp.scrs.neu, Species = c("EPA", "DHA", "C18:1n-9", "C18:3", "C18:2", "C24:0"))



x11()
p2 <- ggplot(scrs.neu, aes(x = NMDS1, y = NMDS3)) +
  geom_point(aes(color = factor(Group)), size = 2) +
  coord_fixed() +
  # stat_ellipse(aes(color = factor(Group)), size = 2, type = "norm", show.legend = FALSE) +
  geom_segment(data = spp.scrs.neu,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS3),
               arrow = arrow(length = unit(0.25, "cm")), colour="grey",
               size = 1.0) +
  geom_text_repel(data = spp.scrs.neu, aes(x = NMDS1, y = NMDS3, label = Species), size = 3,
                  show.legend = FALSE, color = "black")  +
  ylim(-1,1) + #xlim(-2,3) + 
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.title = element_text(size = 12), legend.text = element_text(color = NMDS.treatment.colors, size = 12),
        legend.key.size = unit(2, "lines"))
neu.p2 <- p2 +  scale_color_manual(values = NMDS.treatment.colors,
                                  labels = c("Control", "PUFA", "PUFA + E", "Unenriched"),
                                  name = "") +
  guides(color = guide_legend(override.aes = list(size = 2.0)))
neu.p2














# --------------------------------------- Phospholipid ----------------------------------------------


artemia_phos_bad <- read.csv("02_Artemia_Phos_FA.csv", header = FALSE)

artemia_phos <- data.frame(t(artemia_phos_bad[-1]))
colnames(artemia_phos) <- artemia_phos_bad[, 1]





#...................................................NMDS..................................................................



Art.phos.nmds <- matrix(as.numeric(as.character(unlist(artemia_phos[, 3:20]))), nrow = nrow(artemia_phos))
colnames(Art.phos.nmds) <- colnames(artemia_phos)[3:20]


#Run NMDS with a number of different axes (k)
Art.phos.mds2 <- metaMDS(Art.phos.nmds, k=2)
# Warning message:
#   In metaMDS(Art.phos.nmds, k = 2) :
#   stress is (nearly) zero: you may have insufficient data














#Plotting the 2 pairs of neutral plots together


x11()
pp <- ggarrange(neu.p1, neu.p2, ncol = 1, nrow = 2, common.legend = TRUE, legend = "top") +
  draw_plot_label(label = c("(A)", "(B)"), size = 14, x = c(0.35, 0.35), y = c(0.90, 0.45))

pp




ggsave("01_Art_NMDS.tif", plot = last_plot(), device = "tiff",
       width = 120, height = 160, units = "mm", dpi = 300)
#













