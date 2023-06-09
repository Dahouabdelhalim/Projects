# This R script examines pollinator observation data. Specifically, 
# it tests whether the proportion of floral visits matching syndrome 
# varies by syndrome (Figure 2) using a phylogenetic ANOVA. It uses 
# the phylogenetic tree file "chronogram_v8_2.tre" and the data file 
# "Pollinator_observations.csv"

library(dplyr)
library(geiger) #for phylANOVA
library(car) #levene's test
library(ggplot2)

setwd("~/Dropbox/Costus grant/Dena and Kathleen/FloralEvolution/Dryad_Deposit")

#######
### Prune phylogeny to match pollinator data
#######

#import chronogram
phy <- read.nexus("chronogram_v8_2.tre")
plot(phy, cex=0.4)

#import species pollinator observations
phy.data <- read.csv("Pollinator_observations.csv")

#remove species without observations
phy.data <- filter(phy.data, observed == "yes")
phy.data$syndrome <- as.factor(phy.data$syndrome)

(name.check(phy, phy.data, phy.data$tip_label) -> phyOverlap)
drop.tip(phy, phyOverlap$tree_not_data) -> phyComparativeTree
plot(phyComparativeTree, cex=0.6)

attach(phy.data)

#######
### plot prop visits matching synd by syndrome
#######

mycol <- c("#1F639B","#ED553B")

pdf(file="Fig2.pdf", width=5,height=4)
prop_plot <- ggplot(phy.data, aes(x = syndrome, y = prop_matching, color = syndrome, fill = syndrome)) +
  geom_jitter(size = log(phy.data$Total_visits_observed), show.legend = FALSE, width = 0.3) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  geom_boxplot(fill = NA, show.legend = FALSE) +
  scale_color_manual(values = mycol) +
  ylab("Proportion of visits matching syndrome") +
  ylim(0.5,1)+
  xlab("Pollination syndrome")+
  theme(legend.position = "top")+
  theme_classic()
prop_plot
dev.off()

#######
### Does variance in proportion of visits matching syndrome differ by syndrome?
#######

leveneTest(phy.data$prop_matching ~ phy.data$syndrome)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value  Pr(>F)  
# group  1  6.7828 0.01502 *
#       26                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#######
### PHY ANOVA on proportion of visits matching syndrome by syndrome
#######

dat = phy.data$prop_matching
dat = as.data.frame(dat)
rownames(dat) = phy.data$tip_label
dat = na.omit(dat)
grp = as.factor(phy.data$syndrome)
names(grp)=phy.data$tip_label
t=aov.phylo(dat ~ grp, phyComparativeTree) 
# Analysis of Variance Table
# 
# Response: dat
# Df   Sum-Sq  Mean-Sq F-value   Pr(>F) Pr(>F) given phy
# group      1 0.032215 0.032215  6.3412 0.018288           0.1399
# Residuals 26 0.132085 0.005080   

#not significant with phylogeny
#significant without phylogeny
