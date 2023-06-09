# 03 Feb 2022
# v4.1.1
#-------------


# Cleaned up for publication

### Paired with

#Tank means
# "FA_10d_Neu_01.csv"
# ""FA_10d_Phos_01.csv""




#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(vegan)      # needed for NMDS
library(MASS)       # needed for NMDS
library(ggplot2)    # needed for plotting
# library(grid) # Is this needed?
library(ggrepel)    # needed for geom_text_repel
library(ggpubr)     # needed to plot panels A and B together
library(cowplot)    # needed to plot panels A and B together




setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets WD to folder this is saved in



#...................................................................................................
### #################################### Phase 1 Neutral FA #########################################
#...................................................................................................



Neu.bad <- read.csv("FA_10d_Neu_01.csv", header = FALSE)

Neu <- data.frame(t(Neu.bad[-1]))
colnames(Neu) <- Neu.bad[, 1]



#...................................................NMDS..................................................................


Neu.nmds <- matrix(as.numeric(as.character(unlist(Neu[, 3:18]))), nrow = nrow(Neu))
colnames(Neu.nmds) <- colnames(Neu)[3:18]


#Run NMDS with a number of different axes (k)
Neu.mds2 <- metaMDS(Neu.nmds, k=2, trymax = 200)
#Have to go with 2 - all other number of axes have stress at almost 0, indicating insufficient data
Neu.mds3 <- metaMDS(Neu.nmds, k=3, trymax = 200)
Neu.mds4 <- metaMDS(Neu.nmds, k=4, trymax = 200)
Neu.mds5 <- metaMDS(Neu.nmds, k=5, trymax = 200)
Neu.mds6 <- metaMDS(Neu.nmds, k=6, trymax = 200)

#Calculate how much stress changes with addition of axes
##Rule of thumb, only keep axes that decrease stress by at least ~0.05
###However, don't use more than 2nd axis here bc stress is nearly zero above that
Neu.mds2$stress - Neu.mds3$stress
Neu.mds3$stress - Neu.mds4$stress

#using 2 axes here
plot(Neu.mds2)
Neu.mds2$species  #ZP scores
Neu.mds2$points  #Individual scores
Neu.mds2$stress   #Stress value = 0.027
stressplot(Neu.mds2)  #Diagnostic - pretty good



scrs.neu <- as.data.frame(scores(Neu.mds2, display = "sites")) #values for each axis
scrs.neu <- cbind(scrs.neu, Group = Neu$Treatment)  #Adding the associated treatment group names



#Only fitting vectors for the FA most driving differences between groups (from SIMPER analysis) to minimize clutter on figure
vec.sp.neu <- envfit(Neu.mds2 ~ C20.0 + C20.5 + C18.3 + C22.6 + C18.1n9 + C18.2, as.data.frame(Neu.nmds))
spp.scrs.neu <- as.data.frame(scores(vec.sp.neu, display = "vectors"))
spp.scrs.neu <- cbind(spp.scrs.neu, Species = c("C20:0", "EPA", "C18:3", "DHA", "C18:1n-9", "C18:2"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Axes 1 & 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NMDS.treatment.colors <- c("black", "grey55", "black")


x11()
p1 <- ggplot(scrs.neu, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = factor(Group), shape = factor(Group)), size = 3) +
  coord_fixed() +
  geom_segment(data = spp.scrs.neu,
               aes(x = 0, xend = NMDS1 * 0.5, y = 0, yend = NMDS2 * 0.5),
               arrow = arrow(length = unit(0.25, "cm")), colour = "Black",
               size = 1.0) +
  geom_text_repel(data = spp.scrs.neu, aes(x = NMDS1 * 0.55, y = NMDS2 * 0.55, label = Species), size = 3,
                  show.legend = FALSE, color = "grey50")  +
  xlim(-0.5, 0.5) + ylim(-0.5, 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 12, color = "black"),
        legend.key.size = unit(2, "lines"))
p10.neu <- p1 +  scale_color_manual(values = NMDS.treatment.colors,
                                    labels = c("Control", "PUFA", "PUFA + E"),
                                    name = "") +
  scale_shape_manual(values = c(2, 16, 16),
                     labels = c("Control", "PUFA", "PUFA + E"),
                     name = "") +
  guides(color = guide_legend(override.aes = list(size = 5.0)))
p10.neu












#...................................................................................................
### #################################### Phase 1 Phospholipid FA #########################################
#...................................................................................................




Phos.bad <- read.csv("FA_10d_Phos_01.csv", header = FALSE)

Phos <- data.frame(t(Phos.bad[-1]))
colnames(Phos) <- Phos.bad[, 1]





#...................................................NMDS..................................................................



Phos.nmds <- matrix(as.numeric(as.character(unlist(Phos[, 3:19]))), nrow = nrow(Phos))
colnames(Phos.nmds) <- colnames(Phos)[3:19]


#Run NMDS with a number of different axes (k)
Phos.mds2 <- metaMDS(Phos.nmds, k=2, trymax = 200)
#Have to go with 2 - all other number of axes have stress at almost 0, indicating insufficient data
Phos.mds3 <- metaMDS(Phos.nmds, k=3, trymax = 200)
Phos.mds4 <- metaMDS(Phos.nmds, k=4, trymax = 200)
Phos.mds5 <- metaMDS(Phos.nmds, k=5, trymax = 200)
Phos.mds6 <- metaMDS(Phos.nmds, k=6, trymax = 200)

#Calculate how much stress changes with addition of axes
##Rule of thumb, only keep axes that decrease stress by at least ~0.05
###However, don't use more than 2nd axis here bc stress is nearly zero above that
Phos.mds2$stress - Phos.mds3$stress
Phos.mds3$stress - Phos.mds4$stress

#using 2 axes here
plot(Phos.mds2)
Phos.mds2$species  #ZP scores
Phos.mds2$points  #Individual scores
Phos.mds2$stress   #Stress value = 0.0559
stressplot(Phos.mds2)  #Diagnostic - pretty good




scrs.phos <- as.data.frame(scores(Phos.mds2, display = "sites")) #values for each axis
scrs.phos <- cbind(scrs.phos, Group = Phos$Treatment)  #Adding the associated treatment group names



#Only fitting vectors for the FA most driving differences between groups (from SIMPER analysis) to minimize clutter on figure
vec.sp.phos <- envfit(Phos.mds2 ~ C20.0 + C20.5 + C24.0 + C22.6 + C18.1n9 + C18.2 + C16.0, as.data.frame(Phos.nmds))
spp.scrs.phos <- as.data.frame(scores(vec.sp.phos, display = "vectors"))
spp.scrs.phos <- cbind(spp.scrs.phos, Species = c("C20:0", "EPA", "C24:0", "DHA", "C18:1n-9", "C18:2", "C16:0"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Axes 1 & 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NMDS.treatment.colors <- c("black", "grey55", "black")


x11()
p2 <- ggplot(scrs.phos, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = factor(Group), shape = factor(Group)), size = 3) +
  coord_fixed() +
  geom_segment(data = spp.scrs.phos,
               aes(x = 0, xend = NMDS1*0.5, y = 0, yend = NMDS2*0.5),
               arrow = arrow(length = unit(0.25, "cm")), colour="black",
               size = 1.0) +
  geom_text_repel(data = spp.scrs.phos, aes(x = NMDS1*0.5, y = NMDS2*0.5, label = Species), size = 3,
                  show.legend = FALSE, color = "grey50")  +
  xlim(-0.5, 0.5) + ylim(-0.5, 0.5) +
  theme_classic() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 12, color = "black"),
        legend.key.size = unit(2, "lines"))
p10.phos <- p2 +  scale_color_manual(values = NMDS.treatment.colors,
                                     labels = c("Control", "PUFA", "PUFA + E"),
                                     name = "") +
  scale_shape_manual(values = c(2, 16, 16),
                     labels = c("Control", "PUFA", "PUFA + E"),
                     name = "") +
  guides(color = guide_legend(override.aes = list(size = 5.0)))
p10.phos






















#...................................................................................................
### #################################### Combining the plots into Fig. 1 #########################################
#...................................................................................................



x11()
pp <- ggarrange(p10.neu, p10.phos, ncol = 2, nrow = 1, common.legend = TRUE, legend = "top") +
  draw_plot_label(label = c("(A)", "(B)"), size = 14, x = c(0.11, 0.61), y = c(0.85, 0.85)) +
  draw_plot_label(label = c("Stress = 0.027", "Stress = 0.056"), size = 12, x = c(0.11, 0.61), y = c(0.25, 0.25))

pp


ggsave("Fig2_Phase1_Neu&Phos_Manuscript2.tif", plot = last_plot(), device = "tiff",
       width = 6.5, height = 3.5, units = "in", dpi = 600)












