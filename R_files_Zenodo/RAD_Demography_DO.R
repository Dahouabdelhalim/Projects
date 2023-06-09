# Run the FUNCTIONS, LOAD and CLEAN R scripts prior to producing the outputs

source("RAD_Demography_FUNCTIONS.R")
source("RAD_Demography_LOAD.R")
source("RAD_Demography_CLEAN.R")

# GCTA relatedness --------------------------------------------------------

# Canaries
GCTA_relatedness_CAN <- ggplot(relateGCTA_val_CAN,aes(x = V4))+
  xlim(-0.05,0.7)+
  ylim(0,29)+
  xlab("")+
  ylab("Count")+
  geom_histogram(bins = 140)+
  theme_claudia()+
  facet_wrap( ~ Pop1)

GCTA_relatedness_CAN


# Maderian
GCTA_relatedness_MAD <- ggplot(relateGCTA_val_MAD,aes(x = V4))+
  xlim(-0.05,0.7)+
  ylim(0,29)+
  xlab("Relatedness")+
  ylab("Count")+
  geom_histogram(bins = 140)+
  theme_claudia()+
  facet_wrap( ~ Pop1)

GCTA_relatedness_MAD

# Canaries
GCTA_relatedness_SEL <- ggplot(relateGCTA_val_SEL,aes(x = V4))+
  xlim(-0.2,0.7)+
  ylim(0,29)+
  xlab("")+
  ylab("Count")+
  geom_histogram(bins = 140)+
  theme_claudia()

GCTA_relatedness_SEL


# table of closely related individuals 
close_relatives_CAN <- subset(relateGCTA_val_CAN, V4 >0.2)
close_relatives_CAN <- subset(close_relatives_CAN, Indiv1 != Indiv2)


# Plink relatedness. * FIGURE S1 *-------------------------------------------------------

# Canaries
Plink_relatedness_CAN_plot <- ggplot(Plink_relatedness_CAN,aes(x = rel))+
  theme_claudia()+
  xlim(-0.05,0.7)+
  ylim(0,29)+
  xlab("Relatedness")+
  ylab("Count")+
  geom_histogram(bins = 140)+
  facet_wrap( ~ POP1)

Plink_relatedness_CAN_plot

# Maderia
Plink_relatedness_MAD_plot <- ggplot(Plink_relatedness_MAD,aes(x = rel))+
  theme_claudia()+
  xlim(-0.05,0.7)+
  ylim(0,29)+
  xlab("Relatedness")+
  ylab("Count")+
  geom_histogram(bins = 140)+
  facet_wrap( ~ POP1)

Plink_relatedness_MAD_plot

# Selvagens
Plink_relatedness_SEL_plot <- ggplot(Plink_relatedness_SEL,aes(x = rel))+
  theme_claudia()+
  xlim(-0.05,0.7)+
  ylim(0,29)+
  xlab("Relatedness")+
  ylab("Count")+
  geom_histogram(bins = 140)+
  ggtitle("SG")+
  theme(plot.title = element_text(hjust = 0.5, size = 10))
Plink_relatedness_SEL_plot

Plink_relatedness_SEL_plot <- cowplot::plot_grid(Plink_relatedness_SEL_plot, NULL, NULL, nrow = 1, rel_widths = c(1.19,1,1))


### plot for 3 archipelagos

Plink_relatedness_CAN_plot <- cowplot::plot_grid(NULL, Plink_relatedness_CAN_plot, rel_widths = c(0.1,2), labels = c("A", ""))
Plink_relatedness_MAD_plot <- cowplot::plot_grid(NULL, Plink_relatedness_MAD_plot, rel_widths = c(0.1,2), labels = c("B", ""))
Plink_relatedness_SEL_plot <- cowplot::plot_grid(NULL, Plink_relatedness_SEL_plot, rel_widths = c(0.1,2), labels = c("C", ""))

cowplot::plot_grid(Plink_relatedness_CAN_plot, NULL, Plink_relatedness_MAD_plot, NULL, Plink_relatedness_SEL_plot, 
                   nrow = 5, rel_heights = c(3, 0.04, 1.15, 0.04, 1.13))


### Plot both relatedness measures 

relatedness_CAN <- cowplot::plot_grid(GCTA_relatedness_CAN, Plink_relatedness_CAN_plot, labels = c("A","B"), nrow = 2)
relatedness_CAN
cowplot::plot_grid(GCTA_relatedness_CAN, labels = c("A"), nrow = 1)
cowplot::plot_grid(Plink_relatedness_CAN_plot, labels = c("B"), nrow = 1)

relatedness_MAD <- cowplot::plot_grid(GCTA_relatedness_MAD, Plink_relatedness_MAD_plot, labels = c("A","B"), nrow = 2)
relatedness_MAD

Plink_relatedness_MAD_plot <- cowplot::plot_grid(NULL, Plink_relatedness_MAD_plot, ncol = 2, rel_widths = c(0.1,2))
Plink_relatedness_CAN_plot <- cowplot::plot_grid(NULL, Plink_relatedness_CAN_plot, ncol = 2, rel_widths = c(0.1,2))

Plink_relatedness <- cowplot::plot_grid(Plink_relatedness_CAN_plot, NULL, Plink_relatedness_MAD_plot, NULL, labels = c("A","", "B", ""), nrow = 4, rel_heights = c(3,0.2, 1.2, 0.2))
Plink_relatedness


# Relatedness correlation test. * Correlation results reported in text *  -------------------------------------------

### Plink relatedness calculations but using the --make-grm-gz comand so that the output is the same 

# Canaries

#relateGCTA_val_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_noZ_GCTA_rel.grm.gz")

#relatePlink_val_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_Plinkrel0.2.grm.gz")

#relatePlink_val_CAN <- subset(relatePlink_val_CAN,relatePlink_val_CAN$V1 != relatePlink_val_CAN$V2)
#relateGCTA_val_CAN <- subset(relateGCTA_val_CAN,relateGCTA_val_CAN$V1 != relateGCTA_val_CAN$V2)

#all(relateGCTA_val_CAN$V1 == relatePlink_val_CAN$V1)
#all(relateGCTA_val_CAN$V2 == relatePlink_val_CAN$V2)

#cor(relateGCTA_val_CAN$V4, relatePlink_val_CAN$V4)
#plot(relateGCTA_val_CAN$V4 ~ relatePlink_val_CAN$V4)



# Maderia 

#relateGCTA_val_MAD <- read.table("data/All_SNPs/Berthelots/Madeira_subset/Demography/Relatedness/GCTA_relatedness/PIP_BER_parseV-90-3flag.Allsnp.pop.nl.trimmed_Maderia_subset_maf_trim_0.03_LD_trim_0.4_noZ_GCTA_rel.grm.gz")

#relatePlink_val_MAD <- read.table("data/All_SNPs/Berthelots/Madeira_subset/Demography/Relatedness/Plink_relatedness/PIP_BER_parseV-90-3flag.Allsnp.pop.nl.trimmed_Maderia_subset_maf_trim_0.03_LD_trim_0.4_noZ_Plink_rel.grm.gz")


#relatePlink_val_MAD <- subset(relatePlink_val_MAD,relatePlink_val_MAD$V1 != relatePlink_val_MAD$V2)
#relateGCTA_val_MAD <- subset(relateGCTA_val_MAD,relateGCTA_val_MAD$V1 != relateGCTA_val_MAD$V2)

#all(relateGCTA_val_MAD$V1 == relatePlink_val_MAD$V1)
#all(relateGCTA_val_MAD$V2 == relatePlink_val_MAD$V2)

#cor(relateGCTA_val_MAD$V4, relatePlink_val_MAD$V4)
#cor.test(relateGCTA_val_MAD$V4, relatePlink_val_MAD$V4)
#plot(relateGCTA_val_MAD$V4 ~ relatePlink_val_MAD$V4)


# Selvagens

#relatePlink_val_SEL <- read.table("data/All_SNPs/Berthelots/Selvagens_subset/Relatedness/Plink/Triangle/PIP_BER_parseV-90-3flag.Allsnp.pop.nl.trimmed_SelvagensSubset_maf0.03_noZ_Plink_rel.grm.gz")


#relatePlink_val_SEL <- subset(relatePlink_val_SEL,relatePlink_val_SEL$V1 != relatePlink_val_SEL$V2)
#relateGCTA_val_SEL <- subset(relateGCTA_val_SEL,relateGCTA_val_SEL$V1 != relateGCTA_val_SEL$V2)

#all(relateGCTA_val_SEL$V1 == relatePlink_val_SEL$V1)
#all(relateGCTA_val_SEL$V2 == relatePlink_val_SEL$V2)

#cor(relateGCTA_val_SEL$V4, relatePlink_val_SEL$V4)
#cor.test(relateGCTA_val_SEL$V4, relatePlink_val_SEL$V4)
#plot(relateGCTA_val_SEL$V4 ~ relatePlink_val_SEL$V4)



# PCA * FIGURE 3 A&B * ---------------------------------------------------------------------

# Canaries 

PCA_CAN <- ggplot(pca_CAN,aes(x = PC1,y = PC2,col = Population))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.spacing = unit(0.8, "lines"),
        axis.title.x = element_text(size = 12, colour = "black", vjust = 1),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.line = element_blank(),
        legend.title = element_blank())+
  scale_shape_manual(values=c(2,17,6,19,1,10,15,0,9))+
  scale_color_manual(values=c(dark_blue, dark_blue,dark_blue, fire_orange,fire_orange,fire_orange,purple,purple,purple))+
  geom_point(aes(shape=Population), size = 2)

PCA_CAN <- PCA_CAN + theme(legend.text = element_text(color = "Black", size = 11))

PCA_CAN <- cowplot::plot_grid(NULL, PCA_CAN, rel_widths = c(0.1,2), labels = c("A", ""))
PCA_CAN

# Maderia 

PCA_MAD <- ggplot(pca_MAD,aes(x = PC1,y = PC2,col = Population))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.spacing = unit(0.8, "lines"),
        axis.title.x = element_text(size = 12, colour = "black", vjust = 1),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.line = element_blank(),
        legend.title = element_blank())+
  scale_shape_manual(values=c(10,2,5))+
  scale_color_manual(values=c(yellow,mid_green,grey))+
  geom_point(aes(shape=Population), size = 2)
PCA_MAD

PCA_MAD <- PCA_MAD + theme(legend.text = element_text(color = "Black", size = 11))


###PCA_MAD <- cowplot::plot_grid(NULL, PCA_MAD, NULL, rel_heights = c(0.3,1,0.3), nrow = 3)
PCA_MAD <- cowplot::plot_grid(NULL, PCA_MAD, rel_widths = c(0.1,2), labels = c("B", ""))


# next to eachother... the most space saving
PCA_plots <- cowplot::plot_grid(PCA_CAN, NULL, PCA_MAD,
                                nrow = 1,
                                ncol = 3,
                                rel_heights = c(1,1),
                                rel_widths = c(1,0.05,1))
# on top..
PCA_plots <- cowplot::plot_grid(PCA_CAN, NULL, PCA_MAD,
                                nrow = 3,
                                ncol = 1,
                                rel_heights = c(1, 0.05,1),
                                rel_widths = c(1,1,1))

PCA_plots



# LD by population * FIGURE 2B and S7 * --------------------------------------------------------


# Selvagens ---------------------------------------------------------------
# SG

SG <- ggplot(LD_SG,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  #ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  geom_smooth(method = "loess", colour = fire_orange, fill = fire_orange, size = 0.5)+
  ggtitle("SG")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,1.2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
SG


# Canary Islands ----------------------------------------------------------

# EH

EH <- ggplot(LD_EH,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  #ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  geom_smooth(method = "loess", fill = purple, colour = purple, size = 0.5)+
  ggtitle("EH")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
EH

# LP

LP <- ggplot(LD_LP,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  #ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  geom_smooth(method = "loess", fill = purple, colour = purple, size = 0.5)+
  ggtitle("LP")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,1.2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
LP

# GOM

GOM <- ggplot(LD_GOM,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  #ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  geom_smooth(method = "loess", fill = purple, colour = purple, size = 0.5)+
  ggtitle("GOM")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,1.2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
GOM

# TEID

TEID <- ggplot(LD_TEID,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  #ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  geom_smooth(method = "loess", fill = purple, colour = purple, size = 0.5)+
  ggtitle("TEID")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
TEID

# TF

TF <- ggplot(LD_TF,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  geom_smooth(method = "loess", size = 0.5, fill = purple, colour = purple)+
  ggtitle("TF")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))

TF


# GC

GC <- ggplot(LD_GC,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  geom_smooth(method = "loess", size = 0.5, fill = purple, colour = purple)+
  ggtitle("GC")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,1.2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))

GC


# FV

FV <- ggplot(LD_FV,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  geom_smooth(method = "loess", fill = purple, colour = purple, size = 0.5)+
  ggtitle("FV")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,1.2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
FV


# LZ

LZ <- ggplot(LD_LZ,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  geom_smooth(method = "loess", fill = purple, colour = purple, size = 0.5)+
  ggtitle("LZ")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,1.2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
LZ

# GRA

GRA <- ggplot(LD_GRA,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  geom_smooth(method = "loess", size = 0.5, fill = purple, colour = purple)+
  ggtitle("GRA")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
GRA

# put all the Canary Islands together
# three islands together 
CI_arch_LD <- rbind(LD_EH,LD_LP,LD_GOM,LD_TEID,LD_TF,LD_GC,LD_FV,LD_LZ,LD_GRA)

CI_all_on_1 <- ggplot(CI_arch_LD,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  #theme(axis.text.y = element_text(colour = "white"))+
  ylab(expression(paste("LD r" ^ "2")))+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 1, pch = 1, alpha = 0/10)+
  geom_smooth(data = LD_EH, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, colour = yellow, fill = yellow, se = F)+
  #geom_smooth(data = LD_LP, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, colour = purple, fill = purple, se = F)+
  #geom_smooth(data = LD_GOM, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, colour = purple, fill = purple, se = F)+
  geom_smooth(data = LD_TF, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, colour = purple, fill = purple, se = F)+
  #geom_smooth(data = LD_TEID, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, colour = purple, fill = purple, se = F)+
  #geom_smooth(data = LD_GC, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, colour = purple, fill = purple, se = F)+
  #geom_smooth(data = LD_FV, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, colour = purple, fill = purple, se = F)+
  #geom_smooth(data = LD_LZ, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, colour = purple, fill = purple, se = F)+
  geom_smooth(data = LD_GRA, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, colour = dark_blue, fill = dark_blue, se = F)+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))
CI_all_on_1

# Madeira -----------------------------------------------------------------


# DG

LD_DG$Island <- "DG"

DG <- ggplot(LD_DG,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  #theme(axis.text.y = element_text(colour = "white"))+
  ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  #ylim(0,1)+
  #xlim(0.01,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  #geom_smooth(method = "gam", , formula = y ~ s(x), size = 0.5, colour = purple, fill = purple)+
  geom_smooth(method = "loess", size = 0.5, fill = dark_green, colour = dark_green)+
  ggtitle("DG")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
DG

# M

LD_M$Island <- "M"

M <- ggplot(LD_M,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  #ylim(0,1)+
  #xlim(0.01,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  #geom_smooth(method = "gam", , formula = y ~ s(x), size = 0.5, colour = purple, fill = purple)+
  geom_smooth(method = "loess", size = 0.5, fill = dark_green, colour = dark_green)+
  ggtitle("M")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
M

# PS

LD_PS$Island <- "PS"

PS <- ggplot(LD_PS,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  #theme(axis.text.y = element_text(colour = "white"))+
  #ylab(expression(paste("LD r" ^ "2")))+
  ylab(NULL)+
  xlab(NULL)+
  #ylim(0,1)+
  #xlim(0,150)+
  geom_point(cex = 0.5, pch = 1, alpha = 0.5/10)+
  #geom_smooth(method = "gam", formula = y ~ s(x), size = 0.5, colour = purple, fill = purple)+
  geom_smooth(method = "loess", size = 0.5, fill = dark_green, colour = dark_green)+
  ggtitle("PS")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(limits = c(0,150), expand = c(0,2))+
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01))
PS

# three islands together 
Madeiran_arch_LD <- rbind(LD_PS,LD_DG,LD_M)

Madeira_all_on_1 <- ggplot(Madeiran_arch_LD,aes(x = BP_dis/1000, y= R2))+
  theme_bw()+
  xlab("Distance (Kb)")+
  #theme(axis.text.y = element_text(colour = "white"))+
  ylab(expression(paste("LD r" ^ "2")))+
  ylim(0,1)+
  xlim(0,150)+
  geom_point(cex = 1, pch = 1, alpha = 0/10)+
  geom_smooth(data = LD_PS, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, fill = purple, colour = purple, se = F)+
  geom_smooth(data = LD_DG, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, fill = purple, colour = purple, se = F)+
  geom_smooth(data = LD_M, aes(x= BP_dis/1000, y= R2), method = "loess", size = 0.5, fill = fire_orange, colour = fire_orange, se = F)+
  #geom_smooth(method = "loess")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))
Madeira_all_on_1


# Put plots together  -----------------------------------------------------

CI_line1 <- cowplot::plot_grid(EH, GOM, LP, labels = c("",""),
                               nrow = 1,
                               ncol = 3,
                               rel_heights = c(1,1),
                               rel_widths = c(1,1))


CI_line_2 <- cowplot::plot_grid(TEID, TF, GC, labels = c("",""),
                                nrow = 1,
                                ncol = 3,
                                rel_heights = c(1,1),
                                rel_widths = c(1,1))

M_SG_line_3 <- cowplot::plot_grid(FV, LZ, GRA, labels = c("","", "", ""),
                                  nrow = 1,
                                  ncol = 3,
                                  rel_heights = c(1,1),
                                  rel_widths = c(1,1))

DG_line_4 <- cowplot::plot_grid(M, PS, DG, SG, NULL, NULL, NULL,
                                nrow = 1,
                                ncol = 3,
                                rel_heights = c(1,1),
                                rel_widths = c(1,1))

SG_line_5 <- cowplot::plot_grid(SG, NULL, NULL,
                                nrow = 1,
                                ncol = 3,
                                rel_heights = c(1,1),
                                rel_widths = c(1,1))


library(cowplot)
library(grid)
library(gridExtra)
y.grob <- textGrob(expression(paste("LD r" ^ "2")), 
                   gp=gpar(col="black", fontsize=13), rot=90)

x.grob <- textGrob("Distance (Kb)", 
                   gp=gpar(col="black", fontsize=13))

all_islands <- cowplot::plot_grid(CI_line1,CI_line_2, M_SG_line_3, DG_line_4, SG_line_5, labels = c("","", "", "", ""),
                                  nrow = 5,
                                  ncol = 1,
                                  rel_heights = c(1,1),
                                  rel_widths = c(1,1))
all_islands


grid.arrange(arrangeGrob(all_islands, left = y.grob, bottom = x.grob))

####

Madeiran_plots <- cowplot::plot_grid(M, PS, DG, labels = c("",""),
                                     nrow = 1,
                                     ncol = 3,
                                     rel_heights = c(1,1),
                                     rel_widths = c(1,1))
Madeiran_plots

Madeiran_plots <- grid.arrange(arrangeGrob(Madeiran_plots, left = y.grob, bottom = x.grob))



CI_bottleneck_examples <- cowplot::plot_grid(EH, TF, GRA, labels = c("",""),
                                             nrow = 1,
                                             ncol = 3,
                                             rel_heights = c(1,1),
                                             rel_widths = c(1,1))
CI_bottleneck_examples

CI_bottleneck_examples <- grid.arrange(arrangeGrob(CI_bottleneck_examples, left = y.grob, bottom = x.grob))


plots <- cowplot::plot_grid(Madeiran_plots, NULL, CI_bottleneck_examples, labels = c("B","",""),
                            nrow = 3,
                            ncol = 1,
                            rel_heights = c(1,0.2,1),
                            rel_widths = c(1,1,1))
plots




# Number of SNPs in each section ------------------------------------------

length(unique(LD_SG$SNP_A))

length(unique(LD_DG$SNP_A))
length(unique(LD_PS$SNP_A))
length(unique(LD_M$SNP_A))

length(unique(LD_EH$SNP_A))
length(unique(LD_GRA$SNP_A))
length(unique(LD_LP$SNP_A))
length(unique(LD_TF$SNP_A))
length(unique(LD_TEID$SNP_A))


# TreeMix All Pipits with Tawny Root * FIGURE S2 A&B *--------------------------------------

plot_tree("All_pipits_TAW_0mig", mbar = F, disp = 0.001) # disp sets the distance that the labels are away from the line
get_f("All_pipits_TAW_0mig")

plot_tree("All_pipits_TAW_1mig", mbar = F)
get_f("All_pipits_TAW_1mig")

# you must create the "pop_order" file with one population name on each row to determine the order of comparison for residual plots
plot_resid("All_pipits_TAW_0mig", "pop_order") # no migration events are actually the best here 
plot_resid("All_pipits_TAW_1mig", "pop_order")



# TreeMix Berthelotâ€™s pipit * FIGURE 2A, S3, S4, S5 and S6 *----------------------------

# Windows of 20 SNPs
plot_tree("Berthelots_0mig_k20", ybar = 1.1, disp = 0.0004, lwd = 1, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_1mig_k20", ybar = 0.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_2mig_k20", ybar = 1.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_3mig_k20", ybar = 0.2, disp = 0.0005, arrow = 0.1, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_4mig_k20", ybar = 0.2, disp = 0.0004, arrow = 0.15, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_5mig_k20", ybar = 0.5, disp = 0.0004, arrow = 0.15, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_6mig_k20", ybar = 1.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_7mig_k20", ybar = 1.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_8mig_k20", ybar = 0.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line

get_f("Berthelots_0mig_k20") # determine the variance in relatedness explained by the models, repeat over each model

# create file called pop_order in this location. This is a file one population name per line eg. C_05TF
# create the residuals plot
plot_resid("Berthelots_0mig_k20", "pop_order")

# Windows of 50 SNPs
plot_tree("Berthelots_0mig_k50", ybar = 1.1, disp = 0.0004, lwd = 1, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_1mig_k50", ybar = 0.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_2mig_k50", ybar = 1.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_3mig_k50", ybar = 0.2, disp = 0.0005, arrow = 0.1, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_4mig_k50", ybar = 0.2, disp = 0.0004, arrow = 0.15, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_5mig_k50", ybar = 0.5, disp = 0.0004, arrow = 0.15, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_6mig_k50", ybar = 1.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_7mig_k50", ybar = 1.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_8mig_k50", ybar = 0.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line

get_f("Berthelots_0mig_k50") # determine the variance in relatedness explained by the models, repeat over each model

# create file called pop_order in this location 
# create the residuals plot
plot_resid("Berthelots_0mig_k50", "pop_order")
plot_resid("Berthelots_4mig_k50", "pop_order")

# Likelihood increase plot for k50
likelihood <- matrix(c(686.25, 692.276, 691.259, 711.571, 700.385, 711.154, 715.685, 719.219, 706.2), nrow = 9, ncol = 1)
likelihood <- as.data.frame(likelihood)


likelihood$number_migrations <- c(0,1,2,3,4,5,6,7,8)
likelihood

# Plot the log likelihood of each of the models with different migration parameters
ggplot(data=likelihood,aes(x = number_migrations, y = V1))+
  geom_line(cex=0.3)+
  ylim(650,720)+
  theme_bw()+
  geom_point()+
  ylab("Log likelihood")+
  xlab("Migration events")


# Windows of 80 SNPs
plot_tree("Berthelots_0mig_k80", ybar = 1.1, disp = 0.0004, lwd = 1, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_1mig_k80", ybar = 0.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_2mig_k80", ybar = 1.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_3mig_k80", ybar = 0.2, disp = 0.0005, arrow = 0.1, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_4mig_k80", ybar = 0.2, disp = 0.0004, arrow = 0.15, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_5mig_k80", ybar = 0.5, disp = 0.0004, arrow = 0.15, cex = 1, mbar = T, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_6mig_k80", ybar = 1.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_7mig_k80", ybar = 1.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line
plot_tree("Berthelots_8mig_k80", ybar = 0.1, disp = 0.0004, arrow = 0.15, cex = 1, mbar = F, plus = 0.01) # disp sets the distance that the labels are away from the line

get_f("Berthelots_0mig_k80") # determine the variance in relatedness explained by the models, repeat over each model

# create file called pop_order in this location 
# create the residuals plot
plot_resid("Berthelots_0mig_k80", "pop_order")
