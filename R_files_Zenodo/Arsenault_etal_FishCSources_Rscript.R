#########################################################
#####      NSF MACRO Fish EAA carbon isotopes       #####
#####      E. R. Arsenault et al.                   #####
#####      Last updated August 2022                 #####
#########################################################

# load libraries
library(tidyverse)
library(MixSIAR)
library(ggpubr)
library(rstatix)
library(ggfortify)
library(MASS)
library(lemon)
library(RColorBrewer)
library(vegan)
library(viridis)
library(ggh4x)
library(ggnewscale)
library(cowplot)
library(ggrepel)
library(hrbrthemes)
library(rjags)
library(mcmcOutput)
library(lme4)
library(performance)
library(lmerTest)

#########################################################
##### Make Figure 2 - Plot of Site Characteristics  #####
#########################################################

# read in site locations file
locs <- read.csv("SiteLocations_NSFMACRO.csv")

# change file format for plotting
sitechar.A <- locs %>% 
  pivot_longer(names_to = "name", values_to = "value", cols=c(CanCoverBank, CanCoverMid))
sitechar.B <- locs %>% 
  pivot_longer(names_to = "name", values_to = "value", cols=WettedWidth)
sitechar.C <- locs %>% 
  pivot_longer(names_to = "name", values_to = "value", cols=StreamOrder)

# plot mid-stream and bank canopy cover together (panel A)
sitechar.A$name <- factor(sitechar.A$name, levels = c("CanCoverMid", "CanCoverBank"))
sitechar.A.plot <- ggplot()+
  geom_jitter(data=sitechar.A, aes(x=name, y=value, fill=Ecoregion, shape=Country), 
              size=1, alpha=0.6, width = 0.3, height = 0, color="black")+
  scale_color_manual(guide="none")+
  scale_shape_manual(values = c(22, 24))+
  scale_fill_manual(values = c("white", "grey", "black"), labels = c("Grassland", "Mountain Steppe", "Semi-Arid Terminal Basin"))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text=element_text(size=7, color = "black"), axis.title.y=element_text(size=6))+
  xlab(label = "")+
  ylab(label = "Percent")+
  ylim(c(0, 100))+
  scale_x_discrete(labels = c('Mid Canopy','Bank Canopy'))

# plot wetted width (panel B)
sitechar.B.plot <- ggplot()+
  geom_jitter(data=sitechar.B, aes(x=name, y=value, fill=Ecoregion, shape=Country), 
              size=1, alpha=0.6, width = 0.5, height = 0, color="black")+
  scale_shape_manual(values = c(22, 24))+
  scale_fill_manual(values = c("white", "grey", "black"), guide = "legend", labels = c("Grassland", "Mountain Steppe", "Semi-Arid Terminal Basin"))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text=element_text(size=7, color = "black"), axis.title.y=element_text(size=6))+
  xlab(label = "")+
  ylab(label = "Meters")+
  scale_x_discrete(labels = c('Wetted Width'))

# plot Strahler stream order (panel C)
sitechar.C.plot <- ggplot()+
  geom_jitter(data=sitechar.C, aes(x=name, y=value, fill=Ecoregion, shape=Country), 
              size=1, alpha=0.6, width = 0.5, height = 0, color="black")+
  scale_shape_manual(values = c(22, 24))+
  scale_fill_manual(values = c("white", "grey", "black"), guide = "legend", labels = c("Grassland", "Mountain Steppe", "Semi-Arid Terminal Basin"))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text=element_text(size=7, color = "black"), axis.title.y=element_text(size=6))+
  xlab(label = "")+
  ylab(label = "")+
  ylim(c(1, 5))+
  scale_x_discrete(labels = c('Stream Order'))

# make larger plot from panels A, B, and C
plot.all.site <- ggarrange(sitechar.A.plot, sitechar.B.plot, sitechar.C.plot, ncol=3, widths = c(1, 0.5, 0.5))

# save as Figure 2
ggsave("Figure2.png", units = "mm", width = 125, height = 70, plot.all.site)

#########################################################
##### Import and Normalize d13CEAA Vals for Sources #####
#########################################################

# load source data
# NOTE: this is an aggregated dataset that contains data from the previous literature in addition to original data. see dataset and associated manuscript for citations
sources <- na.omit(read.csv("MACROMS_AggregatedFoodSources.csv", header = T, row.names = NULL))

# normalize source data
sources.rowmeans <- rowMeans(sources[,5:(ncol(sources))])
sources.rowmeans.data <- cbind(sources, sources.rowmeans)
aa <- c("Ile", "Leu", "Phe", "Thr", "Val")
sources.norm <- data.frame(matrix(nrow = nrow(sources.rowmeans.data), ncol=length(aa)))
colnames(sources.norm)<-aa
for(i in aa){
  x <- sources.rowmeans.data[,i]
  y <- sources.rowmeans.data[,"sources.rowmeans"]
  z <- (x-y)
  sources.norm[,i]<-z
}
sources.norm <- cbind(sources.rowmeans.data[,1:4], sources.norm[aa]) #normalized source data final

################################################################
##### Make Supp Fig 1 - Plot d13C Vals for Each AA & Group #####
################################################################

# reformat source dataframe for plotting
resourceaa.plot <- sources %>% 
  pivot_longer(cols = c(Ile, Leu, Phe, Thr, Val), names_to = "AminoAcid")

# calculate mean and SD d13C for each AA and each group
resourceaa.plot.sum <- resourceaa.plot %>% 
  group_by(Group, AminoAcid) %>% 
  dplyr::summarize(meanaa = mean(value), sdaa = sd(value)) %>% 
  arrange(factor(Group, levels = c("Aquatic", "Bacterial", 
                                   "Fungal", "Terrestrial")))

resourceaa.plot.sum$Group <- factor(resourceaa.plot.sum$Group, 
                                    levels = c("Aquatic", "Bacterial", 
                                               "Fungal", "Terrestrial"))

# make Supplemental Figure 1
ggplot(data = resourceaa.plot.sum, aes(x = AminoAcid, y = meanaa, fill = Group, color = Group)) +
  geom_pointrange(aes(ymin=meanaa-sdaa, ymax=meanaa+sdaa), color = "black", position = position_dodge(width=0.7))+
  geom_point(size = 4, shape = 21, stroke = 1, color = "black", position = position_dodge(width=0.7)) +
  ylab(expression(paste(delta^{13}, "C (\\u2030)"))) +
  xlab("")+
  scale_fill_manual(values=c("#0072B2",  "#E69F00", "#D55E00",  "#009E73"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial"),
                    name = "")+
  scale_x_discrete(labels=c("Ile", "Leu", "Phe", "Thr", "Val"))+
  theme(legend.position = c(0.87, 0.84), panel.background=element_rect(colour = "black", fill = NA, size=1),
        panel.grid = element_blank(), text=element_text(size=14, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank())

# save as Supplemental Figure 1
ggsave("SuppFig1.png", width = 8, height = 6)

################################################################
##### Make Supp Fig 2 - Plot PCA Sources by Group + Origin #####
################################################################

# label dataframe for plotting
pca.sources.labels <- sources.norm[c("Group", "Origin")]
row.names(sources.norm) <- paste(sources.norm$Group, sources.norm$Origin, row.names(sources.norm), sep="_")

# run principal components analysis
pca.sources <- prcomp(sources.norm[aa])
summary(pca.sources)

# prepare PCA for plotting
pca.sources.plot <- as.data.frame(pca.sources$x)
pca.sources.plot$group <- sapply(strsplit(as.character(row.names(sources.norm)), "_"), "[[", 1)
pca.sources.plot$origin <- sapply(strsplit(as.character(row.names(sources.norm)), "_"), "[[", 2)
variance.pca.sources <- (pca.sources$sdev)^2
variance.pca.sources.perc <- round(variance.pca.sources/sum(variance.pca.sources) * 100, 1)
perc.pca.sources <- paste(colnames(pca.sources.plot), "(", paste(as.character(variance.pca.sources), "%", ")", sep="") )
loadings.pca.sources <- data.frame(Variables = rownames(pca.sources$rotation), pca.sources$rotation)
pca.sources.plot$group <- factor(pca.sources.plot$group, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial"))
pca.sources.plot$origin <- factor(pca.sources.plot$origin, levels = c("Commercial", "Culture", "Wild"))

# label plot axes with proportion of variance explained
axislabels <- as.data.frame(summary(pca.sources)$importance)
x.axis <- axislabels %>% 
  slice(2) %>% 
  dplyr::select(PC1, PC2)

# plot PCA for sources labeled by group (color and shape) and origin (size)
plot.sources <- ggplot(data = pca.sources.plot, aes(x=PC1, y=PC2, fill=group, color=group, shape=group, size=origin, alpha=group, ellipse=group))+
  stat_ellipse(type = "t", linetype = 1, size=1, level = 0.95)+
  geom_point(stroke=1)+
  scale_alpha_manual(values = c(0.8, 0.8, 0.8, 0.8))+
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73"), labels=c("Aquatic", "Bacterial", "Fungal", "Terrestrial"))+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73"))+
  scale_shape_manual(values=c(21, 22, 24, 23), labels=c("Aquatic", "Bacterial", "Fungal", "Terrestrial"))+
  scale_size_manual(values = c(3, 3, 1), guide=F)+
  geom_segment(data = loadings.pca.sources, aes(x = 0, y = 0, xend = (PC1*5), yend = (PC2*5)), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "black", inherit.aes = FALSE)+
  annotate("text", x = (loadings.pca.sources$PC1*5), y = (loadings.pca.sources$PC2*5),
           label = loadings.pca.sources$Variables, vjust=1)+
  theme(legend.position = "bottom", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=16, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), plot.tag.position = c(0.25, 0.95))+
  xlab(paste("PC1 (", round((x.axis$PC1)*100, 1), "%)", sep = "")) +
  ylab(paste("PC2 (", round((x.axis$PC2)*100, 1), "%)", sep = ""))
ggsave("SuppFig2.png", height = 6, width = 8)

#################################################
##### Check Assumption of Normality for LDA #####
#################################################

# view QQ Plots
ggqqplot(sources.norm$Ile[sources.norm$Group=="Aquatic"])
ggqqplot(sources.norm$Leu[sources.norm$Group=="Aquatic"])
ggqqplot(sources.norm$Phe[sources.norm$Group=="Aquatic"])
ggqqplot(sources.norm$Thr[sources.norm$Group=="Aquatic"])
ggqqplot(sources.norm$Val[sources.norm$Group=="Aquatic"])
ggqqplot(sources.norm$Ile[sources.norm$Group=="Bacterial"])
ggqqplot(sources.norm$Leu[sources.norm$Group=="Bacterial"])
ggqqplot(sources.norm$Phe[sources.norm$Group=="Bacterial"])
ggqqplot(sources.norm$Thr[sources.norm$Group=="Bacterial"])
ggqqplot(sources.norm$Val[sources.norm$Group=="Bacterial"])
ggqqplot(sources.norm$Ile[sources.norm$Group=="Fungal"])
ggqqplot(sources.norm$Leu[sources.norm$Group=="Fungal"])
ggqqplot(sources.norm$Phe[sources.norm$Group=="Fungal"])
ggqqplot(sources.norm$Thr[sources.norm$Group=="Fungal"])
ggqqplot(sources.norm$Val[sources.norm$Group=="Fungal"])
ggqqplot(sources.norm$Ile[sources.norm$Group=="Terrestrial"])
ggqqplot(sources.norm$Leu[sources.norm$Group=="Terrestrial"])
ggqqplot(sources.norm$Phe[sources.norm$Group=="Terrestrial"])
ggqqplot(sources.norm$Thr[sources.norm$Group=="Terrestrial"])
ggqqplot(sources.norm$Val[sources.norm$Group=="Terrestrial"])

# Shapiro-Wilk tests
shapiro.test(sources.norm$Ile[sources.norm$Group=="Aquatic"])
shapiro.test(sources.norm$Leu[sources.norm$Group=="Aquatic"])
shapiro.test(sources.norm$Phe[sources.norm$Group=="Aquatic"])
shapiro.test(sources.norm$Thr[sources.norm$Group=="Aquatic"])
shapiro.test(sources.norm$Val[sources.norm$Group=="Aquatic"])
shapiro.test(sources.norm$Ile[sources.norm$Group=="Bacterial"])
shapiro.test(sources.norm$Leu[sources.norm$Group=="Bacterial"])
shapiro.test(sources.norm$Phe[sources.norm$Group=="Bacterial"])
shapiro.test(sources.norm$Thr[sources.norm$Group=="Bacterial"])
shapiro.test(sources.norm$Val[sources.norm$Group=="Bacterial"])
shapiro.test(sources.norm$Ile[sources.norm$Group=="Fungal"])
shapiro.test(sources.norm$Leu[sources.norm$Group=="Fungal"])
shapiro.test(sources.norm$Phe[sources.norm$Group=="Fungal"])
shapiro.test(sources.norm$Thr[sources.norm$Group=="Fungal"])
shapiro.test(sources.norm$Val[sources.norm$Group=="Fungal"])
shapiro.test(sources.norm$Ile[sources.norm$Group=="Terrestrial"])
shapiro.test(sources.norm$Leu[sources.norm$Group=="Terrestrial"])
shapiro.test(sources.norm$Phe[sources.norm$Group=="Terrestrial"])
shapiro.test(sources.norm$Thr[sources.norm$Group=="Terrestrial"])
shapiro.test(sources.norm$Val[sources.norm$Group=="Terrestrial"])

###############################################
##### Make Supp Fig 3 - LDA Source Groups #####
###############################################

# define group variable for LDA formula
group <- sources.norm$Group
ldaout <- lda(formula = group ~ Ile + Leu + Phe + Thr + Val, data = sources.norm)

# conduct leave-one-out cross validation to determine predictive ability of the LDA
ldaout.cv <- lda(group~Ile + Leu + Phe + Thr + Val, data = sources.norm, CV = TRUE)
table(ldaout.cv$class, group) #classify by "Group" variable
ldaout.p <- predict(ldaout, newdata = sources.norm[,c("Ile", "Leu", "Phe", "Thr", "Val")])$class
ldaout.p
ldatable <- table(ldaout.p, data = group)
accur <- sum(diag(ldatable))/sum(ldatable)*100
accur #overall accuracy of LDA based on leave-one-out cross validation
ldaout.plot <- predict(ldaout, newdata = sources.norm[,c("Ile", "Leu", "Phe", "Thr", "Val")])
ldahist(data = ldaout.plot$x[,1], g = group)

# make plot
newdata <- data.frame(type = group, lda=ldaout.plot$x)
newdata$type <- factor(newdata$type, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial"))
ggplot(newdata) + 
  stat_ellipse(aes(lda.LD1, lda.LD2, colour = type, alpha=type), type = "t", linetype = 1, size=1, level = 0.95)+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73")) + #using cbPalette,  color-blind friendly
  geom_point(aes(lda.LD1, lda.LD2, alpha=type, color=type, fill=type, shape=type, size=type)) +
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73"), labels=c("Aquatic", "Bacterial", "Fungal", "Terrestrial"))+ #using cbPalette,  color-blind friendly
  theme(legend.position = "bottom", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=16, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), plot.tag.position = c(0.25, 0.95))+
  scale_shape_manual(values=c(21, 22, 24, 23), labels=c("Aquatic", "Bacterial", "Fungal", "Terrestrial"))+
  scale_size_manual(values = c(3, 3, 3, 3))+
  scale_alpha_manual(values = c(0.8, 0.8, 0.8, 0.8))+
  xlab(paste("LD1 (", round(prop.table(ldaout$svd^2)[1]*100, 1), "%)", sep = ""))+
  ylab(paste("LD2 (", round(prop.table(ldaout$svd^2)[2]*100, 1), "%)", sep = ""))

# save as Supp Fig 3
ggsave("SuppFig3.png", height = 6, width = 8)

###########################################################
##### Import and Normalize d13CEAA Vals for Consumers #####
###########################################################

# load consumer data
consumer <- na.omit(read.csv("MACROMS_AllFishData.csv", header = T))

# normalize consumer data (value-(sum of d13C for 5 EAAs in a sample/number of EAAs))
consumer.rowmeans <- rowMeans(consumer[,9:(ncol(consumer))])
consumer.rowmeans.data <- cbind(consumer, consumer.rowmeans)
aa <- c("Ile", "Leu", "Phe", "Thr", "Val")
consumer.norm <- data.frame(matrix(nrow = nrow(consumer.rowmeans.data), ncol=length(aa)))
colnames(consumer.norm)<-aa
for(i in aa){
  x <- consumer.rowmeans.data[,i]
  y <- consumer.rowmeans.data[,"consumer.rowmeans"]
  z <- (x-y)
  consumer.norm[,i]<-z
}
consumer.norm <- cbind(consumer.rowmeans.data[,1:8], consumer.norm[aa])

#################################################################
##### Run and Plot PCAs for Source + Consumer d13CEAA Vals  #####
#################################################################

# separate consumer data by ecoregion for PCA (normalized data)
fish.USGR <- consumer.norm %>% 
  filter(EcoregionName=="USGR") %>% 
  mutate(Group = "USGR")
fish.USMT <- consumer.norm %>% 
  filter(EcoregionName=="USMT")%>% 
  mutate(Group = "USMT")
fish.USSA <- consumer.norm %>% 
  filter(EcoregionName=="USSA")%>% 
  mutate(Group = "USSA") 
fish.MNGR <- consumer.norm %>% 
  filter(EcoregionName=="MNGR")%>% 
  mutate(Group = "MNGR") 
fish.MNMT <- consumer.norm %>% 
  filter(EcoregionName=="MNMT")%>% 
  mutate(Group = "MNMT") 
fish.MNSA <- consumer.norm %>% 
  filter(EcoregionName=="MNSA")%>% 
  mutate(Group = "MNSA") 

# run PCA and prepare dataframe for plotting - USGR
pca.data.norm.USGR <- rbind(sources.norm[aa], fish.USGR[aa]) #US grassland
pca.data.norm.labels.USGR <- rbind(sources.norm["Group"], fish.USGR["Group"])
pca.data.norm.USGR <- cbind(pca.data.norm.USGR, pca.data.norm.labels.USGR)
row.names(pca.data.norm.USGR) <- paste(pca.data.norm.USGR$Group, row.names(pca.data.norm.USGR), sep="_")
pca.USGR <- prcomp(pca.data.norm.USGR[aa])
summary(pca.USGR)
pca.USGR.data <- as.data.frame(pca.USGR$x)
pca.USGR.data$group <- sapply(strsplit(as.character(row.names(pca.data.norm.USGR)), "_"), "[[", 1 )
variance.USGR <- (pca.USGR$sdev)^2
variance.USGR.perc <- round(variance.USGR/sum(variance.USGR) * 100, 1)
perc.USGR <- paste(colnames(pca.USGR.data), "(", paste(as.character(variance.USGR.perc), "%", ")", sep="") )
loadings.USGR <- data.frame(Variables = rownames(pca.USGR$rotation), pca.USGR$rotation)
pca.USGR.data$group <- factor(pca.USGR.data$group, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "USGR"))
axislabels.USGR <- as.data.frame(summary(pca.USGR)$importance)
x.axis.USGR <- axislabels.USGR %>% 
  slice(2) %>% 
  dplyr::select(PC1, PC2)

# run PCA and prepare dataframe for plotting - USMT
pca.data.norm.USMT <- rbind(sources.norm[aa], fish.USMT[aa]) #US mountain steppe
pca.data.norm.labels.USMT <- rbind(sources.norm["Group"], fish.USMT["Group"])
pca.data.norm.USMT <- cbind(pca.data.norm.USMT, pca.data.norm.labels.USMT)
row.names(pca.data.norm.USMT) <- paste(pca.data.norm.USMT$Group, row.names(pca.data.norm.USMT), sep="_")
pca.USMT <- prcomp(pca.data.norm.USMT[aa])
summary(pca.USMT) 
pca.USMT.data <- as.data.frame(pca.USMT$x)
pca.USMT.data$group <- sapply(strsplit(as.character(row.names(pca.data.norm.USMT)), "_"), "[[", 1 )
variance.USMT <- (pca.USMT$sdev)^2
variance.USMT.perc <- round(variance.USMT/sum(variance.USMT) * 100, 1)
perc.USMT <- paste(colnames(pca.USMT.data), "(", paste(as.character(variance.USMT.perc), "%", ")", sep="") )
loadings.USMT <- data.frame(Variables = rownames(pca.USMT$rotation), pca.USMT$rotation)
pca.USMT.data$group <- factor(pca.USMT.data$group, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "USMT"))
axislabels.USMT <- as.data.frame(summary(pca.USMT)$importance)
x.axis.USMT <- axislabels.USMT %>% 
  slice(2) %>% 
  dplyr::select(PC1, PC2)

# run PCA and prepare dataframe for plotting - USSA
pca.data.norm.USSA <- rbind(sources.norm[aa], fish.USSA[aa]) #US terminal basin (semi-arid)
pca.data.norm.labels.USSA <- rbind(sources.norm["Group"], fish.USSA["Group"])
pca.data.norm.USSA <- cbind(pca.data.norm.USSA, pca.data.norm.labels.USSA)
row.names(pca.data.norm.USSA) <- paste(pca.data.norm.USSA$Group, row.names(pca.data.norm.USSA), sep="_")
pca.USSA <- prcomp(pca.data.norm.USSA[aa])
summary(pca.USSA) 
pca.USSA.data <- as.data.frame(pca.USSA$x)
pca.USSA.data$group <- sapply(strsplit(as.character(row.names(pca.data.norm.USSA)), "_"), "[[", 1 )
variance.USSA <- (pca.USSA$sdev)^2
variance.USSA.perc <- round(variance.USSA/sum(variance.USSA) * 100, 1)
perc.USSA <- paste(colnames(pca.USSA.data), "(", paste(as.character(variance.USSA.perc), "%", ")", sep="") )
loadings.USSA <- data.frame(Variables = rownames(pca.USSA$rotation), pca.USSA$rotation)
pca.USSA.data$group <- factor(pca.USSA.data$group, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "USSA"))
axislabels.USSA <- as.data.frame(summary(pca.USSA)$importance)
x.axis.USSA <- axislabels.USSA %>% 
  slice(2) %>% 
  dplyr::select(PC1, PC2)

# run PCA and prepare dataframe for plotting - MNGR
pca.data.norm.MNGR <- rbind(sources.norm[aa], fish.MNGR[aa]) #Mongolia grassland
pca.data.norm.labels.MNGR <- rbind(sources.norm["Group"], fish.MNGR["Group"])
pca.data.norm.MNGR <- cbind(pca.data.norm.MNGR, pca.data.norm.labels.MNGR)
row.names(pca.data.norm.MNGR) <- paste(pca.data.norm.MNGR$Group, row.names(pca.data.norm.MNGR), sep="_")
pca.MNGR <- prcomp(pca.data.norm.MNGR[aa])
summary(pca.MNGR) 
pca.MNGR.data <- as.data.frame(pca.MNGR$x)
pca.MNGR.data$group <- sapply(strsplit(as.character(row.names(pca.data.norm.MNGR)), "_"), "[[", 1 )
variance.MNGR <- (pca.MNGR$sdev)^2
variance.MNGR.perc <- round(variance.MNGR/sum(variance.MNGR) * 100, 1)
perc.MNGR <- paste(colnames(pca.MNGR.data), "(", paste(as.character(variance.MNGR.perc), "%", ")", sep="") )
loadings.MNGR <- data.frame(Variables = rownames(pca.MNGR$rotation), pca.MNGR$rotation)
pca.MNGR.data$group <- factor(pca.MNGR.data$group, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "MNGR"))
axislabels.MNGR <- as.data.frame(summary(pca.MNGR)$importance)
x.axis.MNGR <- axislabels.MNGR %>% 
  slice(2) %>% 
  dplyr::select(PC1, PC2)

# run PCA and prepare dataframe for plotting - MNMT
pca.data.norm.MNMT <- rbind(sources.norm[aa], fish.MNMT[aa]) #Mongolia mountain steppe
pca.data.norm.labels.MNMT <- rbind(sources.norm["Group"], fish.MNMT["Group"])
pca.data.norm.MNMT <- cbind(pca.data.norm.MNMT, pca.data.norm.labels.MNMT)
row.names(pca.data.norm.MNMT) <- paste(pca.data.norm.MNMT$Group, row.names(pca.data.norm.MNMT), sep="_")
pca.MNMT <- prcomp(pca.data.norm.MNMT[aa])
summary(pca.MNMT) 
pca.MNMT.data <- as.data.frame(pca.MNMT$x)
pca.MNMT.data$group <- sapply(strsplit(as.character(row.names(pca.data.norm.MNMT)), "_"), "[[", 1 )
variance.MNMT <- (pca.MNMT$sdev)^2
variance.MNMT.perc <- round(variance.MNMT/sum(variance.MNMT) * 100, 1)
perc.MNMT <- paste(colnames(pca.MNMT.data), "(", paste(as.character(variance.MNMT.perc), "%", ")", sep="") )
loadings.MNMT <- data.frame(Variables = rownames(pca.MNMT$rotation), pca.MNMT$rotation)
pca.MNMT.data$group <- factor(pca.MNMT.data$group, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "MNMT"))
axislabels.MNMT <- as.data.frame(summary(pca.MNMT)$importance)
x.axis.MNMT <- axislabels.MNMT %>% 
  slice(2) %>% 
  dplyr::select(PC1, PC2)

# run PCA and prepare dataframe for plotting - MNSA
pca.data.norm.MNSA <- rbind(sources.norm[aa], fish.MNSA[aa]) #Mongolia terminal basin (semi-arid)
pca.data.norm.labels.MNSA <- rbind(sources.norm["Group"], fish.MNSA["Group"])
pca.data.norm.MNSA <- cbind(pca.data.norm.MNSA, pca.data.norm.labels.MNSA)
row.names(pca.data.norm.MNSA) <- paste(pca.data.norm.MNSA$Group, row.names(pca.data.norm.MNSA), sep="_")
pca.MNSA <- prcomp(pca.data.norm.MNSA[aa])
summary(pca.MNSA) 
pca.MNSA.data <- as.data.frame(pca.MNSA$x)
pca.MNSA.data$group <- sapply(strsplit(as.character(row.names(pca.MNSA.data)), "_"), "[[", 1 )
variance.MNSA <- (pca.MNSA$sdev)^2
variance.MNSA.perc <- round(variance.MNSA/sum(variance.MNSA) * 100, 1)
perc.MNSA <- paste(colnames(pca.MNSA.data), "(", paste(as.character(variance.MNSA.perc), "%", ")", sep="") )
loadings.MNSA <- data.frame(Variables = rownames(pca.MNSA$rotation), pca.MNSA$rotation)
pca.MNSA.data$group <- factor(pca.MNSA.data$group, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "MNSA"))
axislabels.MNSA <- as.data.frame(summary(pca.MNSA)$importance)
x.axis.MNSA <- axislabels.MNSA %>% 
  slice(2) %>% 
  dplyr::select(PC1, PC2)

############################################################
##### Make Supp Fig 4 - PCA Plots Sources + Consumers  #####
############################################################

# plot PCA source + consumer d13CEAA vals - USGR
plot.USGR.pca <- ggplot(pca.USGR.data, aes(x=PC1, y=PC2, fill=group, shape=group, size=group, alpha=group, color=group))+
  stat_ellipse(aes(linetype = group), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point(stroke = 1) +
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 1),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ #using cbPalette,  color-blind friendly
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23, 8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 2),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  geom_segment(data = loadings.USGR, aes(x = 0, y = 0, xend = (PC1*5), yend = (PC2*5)), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "black", inherit.aes = FALSE)+
  annotate("text", x = (loadings.USGR$PC1*5), y = (loadings.USGR$PC2*5),
           label = loadings.USGR$Variables, vjust=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=16, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), plot.tag.position = c(0.25, 0.93))+
  labs(tag="D")+
  xlab(paste("PC1 (", round((x.axis.USGR$PC1)*100, 1), "%)", sep = "")) +
  ylab(paste("PC2 (", round((x.axis.USGR$PC2)*100, 1), "%)", sep = ""))

# plot PCA source + consumer d13CEAA vals - USMT
plot.USMT.pca <- ggplot(pca.USMT.data, aes(x=PC1, y=PC2, fill=group, shape=group, size=group, alpha=group, color=group))+
  stat_ellipse(aes(linetype = group), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point(stroke = 1) +
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 1),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ 
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23, 8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 2),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  geom_segment(data = loadings.USMT, aes(x = 0, y = 0, xend = (PC1*5), yend = (PC2*5)), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "black", inherit.aes = FALSE)+
  annotate("text", x = (loadings.USMT$PC1*5), y = (loadings.USMT$PC2*5),
           label = loadings.USMT$Variables, vjust=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=16, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), plot.tag.position = c(0.25, 0.93))+
  labs(tag="E")+
  xlab(paste("PC1 (", round((x.axis.USMT$PC1)*100, 1), "%)", sep = "")) +
  ylab(paste("PC2 (", round((x.axis.USMT$PC2)*100, 1), "%)", sep = ""))

# plot PCA source + consumer d13CEAA vals - USSA
plot.USSA.pca <- ggplot(pca.USSA.data, aes(x=PC1, y=PC2, fill=group, shape=group, size=group, alpha=group, color=group))+
  stat_ellipse(aes(linetype = group), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point(stroke = 1) +
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 1),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ 
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23, 8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 2),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  geom_segment(data = loadings.USSA, aes(x = 0, y = 0, xend = (PC1*5), yend = (PC2*5)), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "black", inherit.aes = FALSE)+
  annotate("text", x = (loadings.USSA$PC1*5), y = (loadings.USSA$PC2*5),
           label = loadings.USSA$Variables, vjust=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=16, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), plot.tag.position = c(0.25, 0.93))+
  labs(tag="F")+
  xlab(paste("PC1 (", round((x.axis.USSA$PC1)*100, 1), "%)", sep = "")) +
  ylab(paste("PC2 (", round((x.axis.USSA$PC2)*100, 1), "%)", sep = ""))

# plot PCA source + consumer d13CEAA vals - MNGR
plot.MNGR.pca <- ggplot(pca.MNGR.data, aes(x=PC1, y=PC2, fill=group, shape=group, size=group, color=group, alpha=group))+
  stat_ellipse(aes(linetype = group), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point(stroke = 1) +
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 1),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ #using cbPalette,  color-blind friendly
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23, 8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 2),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  geom_segment(data = loadings.MNGR, aes(x = 0, y = 0, xend = (PC1*5), yend = (PC2*5)), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "black", inherit.aes = FALSE)+
  annotate("text", x = (loadings.MNGR$PC1*5), y = (loadings.MNGR$PC2*5),
           label = loadings.MNGR$Variables, vjust=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=16, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), plot.tag.position = c(0.25, 0.93))+
  labs(tag="A")+
  xlab(paste("PC1 (", round((x.axis.MNGR$PC1)*100, 1), "%)", sep = "")) +
  ylab(paste("PC2 (", round((x.axis.MNGR$PC2)*100, 1), "%)", sep = ""))

# plot PCA source + consumer d13CEAA vals - MNMT
plot.MNMT.pca <- ggplot(pca.MNMT.data, aes(x=PC1, y=PC2, fill=group, shape=group, size=group, color=group, alpha=group))+
  stat_ellipse(aes(linetype = group), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point(stroke = 1) +
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 1),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ #using cbPalette,  color-blind friendly
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23, 8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 2),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  geom_segment(data = loadings.MNMT, aes(x = 0, y = 0, xend = (PC1*5), yend = (PC2*5)), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "black", inherit.aes = FALSE)+
  annotate("text", x = (loadings.MNMT$PC1*5), y = (loadings.MNMT$PC2*5),
           label = loadings.MNMT$Variables, vjust=1)+
  theme(legend.position = "right", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=16, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), plot.tag.position = c(0.25, 0.93))+
  labs(tag="B")+
  xlab(paste("PC1 (", round((x.axis.MNMT$PC1)*100, 1), "%)", sep = "")) +
  ylab(paste("PC2 (", round((x.axis.MNMT$PC2)*100, 1), "%)", sep = ""))

# plot PCA source + consumer d13CEAA vals - MNSA
plot.MNSA.pca <- ggplot(pca.MNSA.data, aes(x=PC1, y=PC2, fill=group, shape=group, size=group, alpha=group, color=group))+
  stat_ellipse(aes(linetype = group), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point(stroke = 1) +
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 1),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ #using cbPalette,  color-blind friendly
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23, 8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 2),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  geom_segment(data = loadings.MNSA, aes(x = 0, y = 0, xend = (PC1*5), yend = (PC2*5)), 
               arrow = arrow(length = unit(1/2, "picas")),
               color = "black", inherit.aes = FALSE)+
  annotate("text", x = (loadings.MNSA$PC1*5), y = (loadings.MNSA$PC2*5),
           label = loadings.MNSA$Variables, vjust=1, hjust=1)+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=16, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), plot.tag.position = c(0.25, 0.93))+
  labs(tag="C")+
  xlab(paste("PC1 (", round((x.axis.MNSA$PC1)*100, 1), "%)", sep = "")) +
  ylab(paste("PC2 (", round((x.axis.MNSA$PC2)*100, 1), "%)", sep = ""))

# combine all ecoregion plots into one larger plot
plot.all <- grid_arrange_shared_legend(plot.MNGR.pca, plot.MNMT.pca, plot.MNSA.pca, 
                                       plot.USGR.pca, plot.USMT.pca, plot.USSA.pca, 
                                       ncol = 3, nrow = 2, position="bottom")

# save as Supp Fig 4
ggsave("SuppFig4.png", width = 10, height = 8, plot.all)

#####################################################################
##### Run Linear Discriminant Analysis for Diet Classification  #####
#####################################################################

# source d13CEAA vals as training set and consumer d13CEAA vals as test set for classification
trainingset <- sources.norm
grouptraining <- trainingset$Group
testset <- consumer.norm

row.names(testset) <- paste(testset$Fish_ID, testset$SpeciesName, testset$Site, testset$EcoregionName, testset$Country, testset$Feeding, sep="_")

lda.litvals.train <- lda(grouptraining ~ Ile + Leu + Phe + Thr + Val, data = trainingset) # create a training lda (source data)
lda.fishvals.predict <- predict(lda.litvals.train, testset)
summary(lda.fishvals.predict$class)
outputsources <- round(lda.fishvals.predict$posterior, 3)
write.csv(outputsources, file = "LDAposteriors_Fish.csv")
posteriors <- read.csv("LDAposteriors_Fish.csv", header = T)
par(mfrow=c(1,1))
plot(lda.fishvals.predict$x[,1], lda.fishvals.predict$class, col=as.factor(trainingset$Group))

########################################################
##### Make Fig 3 - Plot LDA (Sources + Consumers)  #####
########################################################

# format dataframes for plotting
lda.sources.plot <- data.frame(type=group, ldaout.plot$x) #ldaout.plot$class is predicted category
lda.consumers.plot <- data.frame(type=consumer.norm$EcoregionName, lda.fishvals.predict$x)
lda.all.plot <- rbind(lda.sources.plot, lda.consumers.plot)

# prepare MNGR dataframe for plotting
lda.fish.MNGR.plot <- lda.consumers.plot %>% 
  filter(type=="MNGR")
lda.plot.MNGR <- rbind(lda.sources.plot, lda.fish.MNGR.plot)
lda.plot.MNGR$type <- factor(lda.plot.MNGR$type, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "MNGR"))

# plot LDA - MNGR
plot.MNGR.lda <- ggplot(lda.plot.MNGR, aes(x=LD1, y=LD2, col=type, fill=type, size=type, shape=type, alpha=type))+
  stat_ellipse(aes(linetype = type), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point()+
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 0.8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ 
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23,  8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 1),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=10, color = "black"),
        axis.text=element_text(colour="black", size=8), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), axis.title = element_text(), plot.tag.position = c(0.28, 0.93),
        axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))+
  xlab("")+
  ylab("")+
  labs(tag="A")

# prepare USGR dataframe for plotting
lda.fish.USGR.plot <- lda.consumers.plot %>% 
  filter(type=="USGR")
lda.plot.USGR <- rbind(lda.sources.plot, lda.fish.USGR.plot)
lda.plot.USGR$type <- factor(lda.plot.USGR$type, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "USGR"))

# plot LDA - USGR
plot.USGR.lda <- ggplot(lda.plot.USGR, aes(x=LD1, y=LD2, col=type, fill=type, size=type, shape=type, alpha=type))+
  stat_ellipse(aes(linetype = type), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point()+
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 0.8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ 
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23,  8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 1),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=10, color = "black"),
        axis.text=element_text(colour="black", size=8), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), axis.title = element_text(), plot.tag.position = c(0.28, 0.93),
        axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))+
  xlab(paste("LD1 (", round(prop.table(lda.litvals.train$svd^2)[1]*100, 2), "%)", sep = ""))+
  ylab(paste("LD2 (", round(prop.table(lda.litvals.train$svd^2)[2]*100, 2), "%)", sep = ""))+
  labs(tag="D")

# prepare MNMT dataframe for plotting
lda.fish.MNMT.plot <- lda.consumers.plot %>% 
  filter(type=="MNMT")
lda.plot.MNMT <- rbind(lda.sources.plot, lda.fish.MNMT.plot)
lda.plot.MNMT$type <- factor(lda.plot.MNMT$type, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "MNMT"))

# plot LDA - MNMT
plot.MNMT.lda <- ggplot(lda.plot.MNMT, aes(x=LD1, y=LD2, col=type, size=type, fill=type, shape=type, alpha=type))+
  stat_ellipse(aes(linetype = type), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point()+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ 
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 0.8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23 ,8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 1),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=10, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), axis.title = element_text(), plot.tag.position = c(0.28, 0.93))+
  xlab("")+
  ylab("")+
  labs(tag="B")

# prepare USMT dataframe for plotting
lda.fish.USMT.plot <- lda.consumers.plot %>% 
  filter(type=="USMT")
lda.plot.USMT <- rbind(lda.sources.plot, lda.fish.USMT.plot)
lda.plot.USMT$type <- factor(lda.plot.USMT$type, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "USMT"))

# plot LDA - USMT
plot.USMT.lda <- ggplot(lda.plot.USMT, aes(x=LD1, y=LD2, col=type, fill=type, size=type, shape=type, alpha=type))+
  stat_ellipse(aes(linetype = type), level = 0.95, size = 1)+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ 
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point()+
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 0.8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23, 8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 1),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=10, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), axis.title = element_text(), plot.tag.position = c(0.28, 0.93))+
  xlab("")+
  ylab("")+
  labs(tag="E")

# prepare MTSA dataframe for plotting
lda.fish.MNSA.plot <- lda.consumers.plot %>% 
  filter(type=="MNSA")
lda.plot.MNSA <- rbind(lda.sources.plot, lda.fish.MNSA.plot)
lda.plot.MNSA$type <- factor(lda.plot.MNSA$type, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "MNSA"))

# plot LDA - MNSA
plot.MNSA.lda <- ggplot(lda.plot.MNSA, aes(x=LD1, y=LD2, col=type, fill=type, size=type, shape=type, alpha=type))+
  stat_ellipse(aes(linetype = type), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point()+
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 0.8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ 
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23, 8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 1),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=10, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), axis.title = element_text(), plot.tag.position = c(0.28, 0.93))+
  xlab("")+
  ylab("")+
  labs(tag="C")

# prepare USSA dataframe for plotting
lda.fish.USSA.plot <- lda.consumers.plot %>% 
  filter(type=="USSA")
lda.plot.USSA <- rbind(lda.sources.plot, lda.fish.USSA.plot)
lda.plot.USSA$type <- factor(lda.plot.USSA$type, levels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "USSA"))

# plot LDA - USSA
plot.USSA.lda <- ggplot(lda.plot.USSA, aes(x=LD1, y=LD2, col=type, fill=type, size=type, shape=type, alpha=type))+
  stat_ellipse(aes(linetype = type), level = 0.95, size = 1)+
  scale_linetype_manual(values=c(1, 1, 1, 1, 0), guide="none")+
  geom_point()+
  scale_alpha_manual(values = c(0.5, 0.5, 0.5, 0.5, 0.8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+ 
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#009E73", "#000000"),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_shape_manual(values=c(21, 22, 24, 23, 8),
                     labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  scale_size_manual(values = c(2, 2, 2, 2, 1),
                    labels = c("Aquatic", "Bacterial", "Fungal", "Terrestrial", "Fish Consumers"))+
  theme(legend.position = "none", panel.background=element_rect(colour = "black", size=1, fill = NA),
        panel.grid = element_blank(), text=element_text(size=10, color = "black"),
        axis.text=element_text(colour="black"), legend.key = element_rect(fill = NA),
        legend.title = element_blank(), axis.title = element_text(), plot.tag.position = c(0.28, 0.93))+
  xlab("")+
  ylab("")+
  labs(tag="F")

# combine all ecoregion plots into one figure
plot.all.lda <- grid_arrange_shared_legend(plot.MNGR.lda, plot.MNMT.lda, plot.MNSA.lda, 
                                           plot.USGR.lda, plot.USMT.lda, plot.USSA.lda, 
                                           ncol = 3, nrow = 2, position="bottom")

# save as Figure 3
ggsave("Figure3.png", units = "mm", width = 173, height = 150, plot.all.lda)

#########################################################
##### Determine dietary contributions using MixSIAR #####
#########################################################

# prepare source data file with group mean and SD for normalized d13CEAA values
sources.norm.model <- sources.norm %>%
  dplyr::group_by(Group) %>%
  dplyr::summarize(MeanIle=mean(Ile), MeanLeu=mean(Leu), MeanPhe=mean(Phe), MeanThr=mean(Thr), MeanVal=mean(Val),
                   SDIle=sd(Ile), SDLeu=sd(Leu), SDPhe=sd(Phe), SDThr=sd(Thr), SDVal=sd(Val)) %>% 
  dplyr::mutate(n=c(nrow(sources.norm[sources.norm$Group=="Aquatic",]), 
                    nrow(sources.norm[sources.norm$Group=="Bacterial",]),
                    nrow(sources.norm[sources.norm$Group=="Fungal",]),
                    nrow(sources.norm[sources.norm$Group=="Terrestrial",])))
write_csv(sources.norm.model, "sources.aggregated.norm.csv")

# prepare file of small, nonzero values for TDF
discrim.model <- data.frame(Group = c("Aquatic", "Bacterial", "Fungal", "Terrestrial"),
                            MeanIle  = c(0.1, 0.1, 0.1, 0.1),
                            MeanLeu = c(0.1, 0.1, 0.1, 0.1),
                            MeanPhe = c(0.1, 0.1, 0.1, 0.1),
                            MeanThr = c(0.1, 0.1, 0.1, 0.1),
                            MeanVal = c(0.1, 0.1, 0.1, 0.1),
                            SDIle = c(0.1, 0.1, 0.1, 0.1),
                            SDLeu = c(0.1, 0.1, 0.1, 0.1),
                            SDPhe = c(0.1, 0.1, 0.1, 0.1),
                            SDThr = c(0.1, 0.1, 0.1, 0.1),
                            SDVal = c(0.1, 0.1, 0.1, 0.1))
colnames(discrim.model)[1] <- ""
write_csv(discrim.model, "discrim.csv") 

# prepare file for normalized consumer d13CEAA values
write_csv(consumer.norm, "consumer.norm.model.csv") 

# create MixSIAR model with Site and SpeciesName as fixed, nested factors
model.macro <- load_mix_data(filename = "consumer.norm.model.csv",
                             iso_names = c("Ile", "Leu", "Phe", "Thr", "Val"),
                             factors = c("Site", "SpeciesName"), fac_random = c(FALSE, FALSE), fac_nested = c(FALSE, TRUE), cont_effects = NULL)

# load MixSIAR source file based on means and SD of aggregated resource d13CEAA values
model.sources.macro <- load_source_data(filename="sources.aggregated.norm.csv", source_factors = NULL,
                                        conc_dep=FALSE, data_type="means", model.macro)

# load trophic discrimination factor data
discr.sources.macro <- load_discr_data(filename="discrim.csv", model.macro)

# plot source and consumer data to check assumptions for mixing models (pairwise xy plot for each EAA combo)
plot_data(filename="model.test", plot_save_pdf = F, plot_save_png = F, 
          mix=model.macro, source=model.sources.macro, discr=discr.sources.macro, return_obj = T)

# set up error structure and write model
model_filename <- "model.macro"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, model.macro, model.sources.macro)

# run MixSIAR model
set.seed(2022) #set seed to make output reproducible
jags.model.macro <- run_model(run="long", model.macro, model.sources.macro, discr.sources.macro, model_filename, 
                              resid_err = resid_err, process_err = process_err) 

# set output options
output_options <- list(summary_save = TRUE,
                       summary_name = "output_model.macro",
                       sup_post = TRUE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot.macro",
                       sup_xy = FALSE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = TRUE, #FALSE
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics_model.macro",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       return_obj = TRUE,
                       diag_save_ggmcmc = TRUE)

# to address graphics error
options(expressions=100000)

# to allow printing of entire summary table
options(max.print = 10000)

# view model output plots, summary, and diagnostics
output_JAGS(jags.model.macro, model.macro, model.sources.macro, output_options) #much lower DIC than NULL model, proceed with this model

# convert MixSIAR model into mcmc object
jags.mcmc<-as.mcmc(jags.model.macro)
jags.sum<-summary(jags.mcmc)
jags.sum.stat<-jags.sum$statistics

# get MCE (Monte carlo standard error) to determine model precision
getRhat(jags.model.macro, bad = NA, sort = FALSE) #no values for Rhat <1.05

# get MCE (Monte carlo standard error) to determine model precision
getMCE(jags.model.macro, bad = NA, pc = TRUE, sort = FALSE) #no values greater than 5

#visual examination of trace and density plots (Rhat should be <1.05 and MCEpc should be less than 5, as quantified above)
diagPlot(jags.model.macro)

#############################################################
##### Make Figs 4 + 5 - Summarize Dietary Contributions #####
#############################################################

# format LDA posterior dataframe to calculate summary statistics
lda.outputs_all <- posteriors %>% 
  separate(X, c("Fish_ID", "Species", "Site", "Ecoregion", "Country", "Feeding")) %>% 
  mutate(EcoregionName=Ecoregion) %>% 
  dplyr::select(Fish_ID, Species, Country, EcoregionName, Site, Feeding, Aquatic, Bacterial, Fungal, Terrestrial)

# site-level summary statistics - LDA
lda.outputs_site <- lda.outputs_all %>% 
  dplyr::group_by(Site) %>% 
  dplyr::summarize(meanAqua = mean(Aquatic)*100, sdAqua = sd(Aquatic)*100,
                   meanBact = mean(Bacterial)*100, sdBact = sd(Bacterial)*100,
                   meanFung = mean(Fungal)*100, sdFung = sd(Fungal)*100,
                   meanTerr = mean(Terrestrial)*100, sdTerr = sd(Terrestrial)*100)

# ecoregion-level summary statistics - LDA
lda.outputs_ecoregion <- lda.outputs_all %>% 
  dplyr::group_by(EcoregionName) %>% 
  dplyr::summarize(meanAqua = mean(Aquatic)*100, sdAqua = sd(Aquatic)*100,
                   meanBact = mean(Bacterial)*100, sdBact = sd(Bacterial)*100,
                   meanFung = mean(Fungal)*100, sdFung = sd(Fungal)*100,
                   meanTerr = mean(Terrestrial)*100, sdTerr = sd(Terrestrial)*100)

# country-level summary statistics - LDA
lda.outputs_country <- lda.outputs_all %>% 
  dplyr::group_by(Country) %>% 
  dplyr::summarize(meanAqua = mean(Aquatic)*100, sdAqua = sd(Aquatic)*100,
                   meanBact = mean(Bacterial)*100, sdBact = sd(Bacterial)*100,
                   meanFung = mean(Fungal)*100, sdFung = sd(Fungal)*100,
                   meanTerr = mean(Terrestrial)*100, sdTerr = sd(Terrestrial)*100)

# feeding group-level summary statistics - LDA
lda.outputs_feeding <- lda.outputs_all %>% 
  dplyr::group_by(Feeding) %>% 
  dplyr::summarize(meanAqua = mean(Aquatic)*100, sdAqua = sd(Aquatic)*100,
                   meanBact = mean(Bacterial)*100, sdBact = sd(Bacterial)*100,
                   meanFung = mean(Fungal)*100, sdFung = sd(Fungal)*100,
                   meanTerr = mean(Terrestrial)*100, sdTerr = sd(Terrestrial)*100)

# format LDA posteriors dataframe for plotting
lda.outputs_all.plot <- lda.outputs_all %>% 
  pivot_longer(names_to = "type", values_to = "value", cols=c(Aquatic, Bacterial, Fungal, Terrestrial)) %>% 
  group_by(EcoregionName, Site, Species, type) %>% 
  summarize(meanPost = mean(value), sdPost = sd(value)) %>% 
  mutate(plot_label = case_when(EcoregionName == "MNGR" ~ "A",
                                EcoregionName == "MNMT" ~ "B",
                                EcoregionName == "MNSA" ~ "C",
                                EcoregionName == "USGR" ~ "D",
                                EcoregionName == "USMT" ~ "E",
                                EcoregionName == "USSA" ~ "F"))
lda.outputs_all.plot[is.na(lda.outputs_all.plot)] <- 0

lda.plot.all.outputs <- ggplot(lda.outputs_all.plot, aes(x=type, y=meanPost, fill=type, color=type, shape=type)) +
  geom_pointrange(data = lda.outputs_all.plot, aes(ymin=meanPost-sdPost, ymax=meanPost+sdPost), position=position_jitter(width=0.4), size = 0.5, alpha = 0.5)+
  scale_fill_manual(values = c("#0072B2",  "#E69F00", "#D55E00",  "#009E73"),
                    label = c("Aquatic", "Bacterial", "Fungal", "Terrestrial")) +
  scale_color_manual(values = c("#0072B2",  "#E69F00", "#D55E00",  "#009E73"),
                     label = c("Aquatic", "Bacterial", "Fungal", "Terrestrial")) +
  scale_shape_manual(values = c(21, 22, 24, 23),
                     label = c("Aquatic", "Bacterial", "Fungal", "Terrestrial"))+
  scale_y_continuous(breaks = seq(0.25, 1.0, 0.25))+
  theme(
    legend.position="bottom",
    plot.title = element_text(size=9), axis.title = element_blank(), axis.text.x = element_blank(), 
    legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_rect("white"), panel.border = element_rect(color="black", fill=NA, size=1), 
    plot.subtitle = element_blank(), legend.box.background = element_blank(),
    legend.key = element_rect(fill = NA), axis.text = element_text(size = 9, color="black"), 
    legend.text = element_text(size=9), strip.text = element_blank(),
    strip.background = element_blank()) +
  xlab("") +
  geom_text(data = lda.outputs_all.plot, aes(label = plot_label), x = 4.3, y = 0.97, color = "black") +
  facet_wrap(~EcoregionName, ncol = 3)
  
# save as Fig 4
ggsave("Figure4.png", units = "mm", width = 173, height = 150, lda.plot.all.outputs)

# read in MixSIAR models results summary
mix.outputs_all <- read.delim("output_model.macro.txt", sep = "", skip = 6)
mix.outputs_all$Sample <- rownames(mix.outputs_all)
mix.outputs_all <- mix.outputs_all %>% 
  separate(Sample, into = c("p", "Site", "SpeciesName", "Group"))
mix.outputs_all <- mix.outputs_all[-c(1, 2, 3, 4, 5),] %>% 
  dplyr::select(Site, SpeciesName, Group, Mean, SD)

# format MixSIAR model results summary to be able to calculate summary statistics
mix.outputs_all_wider <- mix.outputs_all %>% 
  mutate(EcoregionName = str_sub(Site, 1, 4)) %>%
  mutate(Country = str_sub(Site, 1, 2)) %>%
  pivot_wider(names_from = Group, values_from = Mean, id_cols = c("Country", "EcoregionName", "Site", "SpeciesName"), id_expand = F) %>% 
  dplyr::select(Country, EcoregionName, Site, SpeciesName, Aquatic, Bacterial, Fungal, Terrestrial)

# site-level summary statistics - MixSIAR
mix.outputs_site <- mix.outputs_all_wider %>% 
  dplyr::group_by(Site, SpeciesName) %>% 
  dplyr::summarize(meanAqua = mean(Aquatic)*100, sdAqua = sd(Aquatic)*100,
                   meanBact = mean(Bacterial)*100, sdBact = sd(Bacterial)*100,
                   meanFung = mean(Fungal)*100, sdFung = sd(Fungal)*100,
                   meanTerr = mean(Terrestrial)*100, sdTerr = sd(Terrestrial)*100) %>% 
  print(n=1000)

# ecoregion-level summary statistics - MixSIAR
mix.outputs_ecoregion <- mix.outputs_all_wider %>% 
  dplyr::group_by(EcoregionName) %>% 
  dplyr::summarize(meanAqua = mean(Aquatic)*100, sdAqua = sd(Aquatic)*100,
                   meanBact = mean(Bacterial)*100, sdBact = sd(Bacterial)*100,
                   meanFung = mean(Fungal)*100, sdFung = sd(Fungal)*100,
                   meanTerr = mean(Terrestrial)*100, sdTerr = sd(Terrestrial)*100)

# country-level summary statistics - MixSIAR
mix.outputs_country <- mix.outputs_all_wider %>% 
  dplyr::group_by(Country) %>% 
  dplyr::summarize(meanAqua = mean(Aquatic)*100, sdAqua = sd(Aquatic)*100,
                   meanBact = mean(Bacterial)*100, sdBact = sd(Bacterial)*100,
                   meanFung = mean(Fungal)*100, sdFung = sd(Fungal)*100,
                   meanTerr = mean(Terrestrial)*100, sdTerr = sd(Terrestrial)*100)

# format MixSIAR output summary for plotting
mix.outputs_all.plot <- mix.outputs_all %>% 
  mutate(EcoregionName = str_sub(Site, 1, 4)) %>%
  mutate(Country = str_sub(Site, 1, 2)) %>%
  mutate(plot_label = case_when(EcoregionName == "MNGR" ~ "A",
                                EcoregionName == "MNMT" ~ "B",
                                EcoregionName == "MNSA" ~ "C",
                                EcoregionName == "USGR" ~ "D",
                                EcoregionName == "USMT" ~ "E",
                                EcoregionName == "USSA" ~ "F"))

# plot MixSIAR diet estimates by species and ecoregion
mix.plot.all.outputs <- ggplot(data = mix.outputs_all.plot, aes(x=Group, y=Mean, fill=Group, color=Group, shape=Group)) +
  geom_pointrange(data = mix.outputs_all.plot, aes(ymin=Mean-SD, ymax=Mean+SD), position=position_jitter(width=0.4), size = 0.5, alpha = 0.5)+
  scale_fill_manual(values = c("#0072B2",  "#E69F00", "#D55E00",  "#009E73"),
                    label = c("Aquatic", "Bacterial", "Fungal", "Terrestrial")) +
  scale_color_manual(values = c("#0072B2",  "#E69F00", "#D55E00",  "#009E73"),
                     label = c("Aquatic", "Bacterial", "Fungal", "Terrestrial")) +
  scale_shape_manual(values = c(21, 22, 24, 23),
                     label = c("Aquatic", "Bacterial", "Fungal", "Terrestrial"))+
  scale_y_continuous(breaks = seq(0.25, 1.0, 0.25))+
  theme(
    legend.position="bottom",
    plot.title = element_text(size=9), axis.title = element_blank(), axis.text.x = element_blank(), 
    legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_rect("white"), panel.border = element_rect(color="black", fill=NA, size=1), 
    plot.subtitle = element_blank(), legend.box.background = element_blank(),
    legend.key = element_rect(fill = NA), axis.text = element_text(size = 9, color="black"), 
    legend.text = element_text(size=9), strip.text = element_blank(),
    strip.background = element_blank()) +
  xlab("") +
  geom_text(data = mix.outputs_all.plot, aes(label = plot_label), x = 4.3, y = 0.97, color = "black") +
  facet_wrap(~EcoregionName, ncol = 3)

# save as Fig 5
ggsave("Figure5.png", units = "mm", width = 173, height = 150, mix.plot.all.outputs)

############################################################################
##### Run Linear Mixed-Effects Models - Habitat + Source Relationships #####
############################################################################

# make single file that includes both MixSIAR model outputs and site characteristics
lmer.habitat_locs <- merge(mix.outputs_all_wider, locs, by=c("Site"))

# add Family-level information for fish consumers
lmer.habitat_locs <- lmer.habitat_locs %>% 
  mutate(FamilyName = case_when(SpeciesName == "Thymallusgrubei" | SpeciesName == "Thymallusbrevirostris" | SpeciesName == "Salvelinusfontinalis" | 
                                  SpeciesName == "Brachymystaxlenok" | SpeciesName == "Proposiumwilliamsoni" |
                                  SpeciesName == "Oncorhynchusclarki" | SpeciesName == "Oncorhynchusmykiss" | 
                                  SpeciesName =="SalmotruttaxSalvelinusfontinalis" | SpeciesName == "Salmotrutta"  ~ "Salmonidae",
                                SpeciesName == "Barbatulabarbatula" ~ "Nemacheilidae",
                                SpeciesName == "Rhynchocyprislagowskii" | SpeciesName == "Phoxinusphoxinus" |
                                  SpeciesName == "Rhinichthyscataractae" | SpeciesName == "Notropisstramineus" |
                                  SpeciesName =="Rhinichthysosculus" | SpeciesName == "Oreoleuciscuspotanini" |
                                  SpeciesName == "Semotilusatromaculatus" ~ "Leuciscidae",
                                SpeciesName == "Micropterusdolomieu" | SpeciesName == "Lepomiscyanellus" ~ "Centrarchidae",
                                SpeciesName == "Ameirusmelas" ~ "Ictaluridae",
                                SpeciesName == "Cottusbeldingii" ~ "Cottidae",
                                SpeciesName == "Gobiogobiocynocephalus" ~ "Gobionidae",
                                SpeciesName == "Lotalota" ~ "Lotidae",
                                SpeciesName == "Catostomusplatyrhynchus" | SpeciesName == "Catostomustahoensis" |
                                  SpeciesName == "Catostomuscommersonii" | SpeciesName == "Moxostomamacrolepidotum" ~ "Catostomidae",
                                SpeciesName == "Cyprinuscarpio" ~ "Cyprinidae",
                                SpeciesName == "Esoxlucius" ~ "Esocidae",
                                SpeciesName == "Barbatulaconilobus" | SpeciesName == "Barbatulatoni" ~ "Nemacheilidae",
                                SpeciesName == "Cobitismelanoleuca" ~ "Cobitidae",
                                TRUE ~ NA_character_))

# build models (one for each resource category)
habitat_model.aqua <- lmer(Aquatic*100 ~ CanCoverMid + CanCoverBank + WettedWidth + StreamOrder + (1 | Country.x/Ecoregion/Site) + (1 | FamilyName/SpeciesName), data = lmer.habitat_locs)
habitat_model.bact <- lmer(Bacterial*100 ~ CanCoverMid + CanCoverBank + WettedWidth + StreamOrder + (1 | Country.x/Ecoregion/Site) + (1 | FamilyName/SpeciesName), data = lmer.habitat_locs)
habitat_model.fung <- lmer(Fungal*100 ~ CanCoverMid + CanCoverBank + WettedWidth + StreamOrder + (1 | Country.x/Ecoregion/Site) + (1 | FamilyName/SpeciesName), data = lmer.habitat_locs)
habitat_model.terr <- lmer(Terrestrial*100 ~ CanCoverMid + CanCoverBank + WettedWidth + StreamOrder + (1 | Country.x/Ecoregion/Site) + (1 | FamilyName/SpeciesName), data = lmer.habitat_locs)

# view diagnostic plots for models
check_model(habitat_model.aqua)
check_model(habitat_model.bact)
check_model(habitat_model.fung)
check_model(habitat_model.terr)

# print model results (including p-values using Satterthwaite's method)
summary(habitat_model.aqua)
summary(habitat_model.bact)
summary(habitat_model.fung)
summary(habitat_model.terr)

# check 95% confidence intervals for each model term (overlap with zero indicates non-significant result)
confint(habitat_model.aqua, level = 0.95)
confint(habitat_model.bact, level = 0.95)
confint(habitat_model.fung, level = 0.95)
confint(habitat_model.terr, level = 0.95)