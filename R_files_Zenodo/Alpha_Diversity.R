########## packges load ###################################################################################
library(vegan)
library(car)
library (DESeq2); packageVersion("DESeq2")
library(Rfast) #to take row maximums
library (ggthemes)
library(data.table); packageVersion("data.table")
library(adaptiveGPCA)
library(ggrepel)
library(phyloseq)
library(plyr)
library(multcompView) # convert this table to a compact letter display
library(rcompanion)  # convert to table to use multcompLetters
############ load clean and rerified data for analysis ###############################################################
data_in <- readRDS(file="data_in_rarefied_8000.rds")

# Sample size
temp<-data.frame(sample_data(data_in))
tapply(temp$WeekSinceBreeding, temp$ControlGroupWeek, length)
# Pre_migration   Fall   Winter_fields   Spring_fields 
#     36            42            32            38 

tapply(temp$WeekSinceBreeding, temp$GroupWithControl, length)

Col2<-c("#ffa321","#48b823","#f7cdba","#d0eef5","#559bab","red4")
Col1<-c("#ffa321","#48b823","#559bab","red4")

###########################  Claculate Alpha-Diversity ##########################################################################
data_in.richness<- estimate_richness(data_in, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
data_in.richness.meta <- cbind(sample_data(data_in),data_in.richness)

############# Phylogenetic diversity ##########################
library(picante)
ps0.rar.asvtab <- as.data.frame(data_in@otu_table)
ps0.rar.tree <- data_in@phy_tree

# We first need to check if the tree is rooted or not 
data_in@phy_tree
## Rooted; includes branch lengths.
# it is a rooted tree
# run the diversity code
df.pd <- pd(ps0.rar.asvtab, ps0.rar.tree,include.root=T) 

############## Merge together ##########################
# Add the rownames to diversity table
df.pd$sam_name <- rownames(df.pd)
# Add the rownames as a new colum for easy integration later.
data_in.richness.meta$sam_name <- rownames(data_in.richness.meta)
# merge with the regular diversity data frame 
data_in.richness.meta <- merge(data_in.richness.meta,df.pd, by = "sam_name")

############# Add "Shape" variable for plotting ##########################
Shape<-0
Shape[data_in.richness.meta$GroupWithControl=="Control_Hula"]<-"ControlHula"
Shape[data_in.richness.meta$GroupWithControl=="Control_Other"]<-"ControlOther"
Shape[data_in.richness.meta$GroupWithControl!="Control_Hula" & data_in.richness.meta$GroupWithControl!="Control_Other"]<-"data"
Shape<-as.factor(Shape)
data_in.richness.meta$Shape<-Shape

# for plotting fall only
data_in.richness_Fall <- data_in.richness.meta[data_in.richness.meta$ControlGroupWeek=="Fall",]
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#=(1) Observed =======================================================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

############# Statitival analysis ##########################

# ++++++++++ (a) without  control 
#  we use Kruskal-Wallis (non-parametric equivalent of ANOVA)
kruskal.test(Observed ~ ControlGroupWeek, data=data_in.richness.meta)

# Wilcoxon-Mann-Whitney post hoc:==============================================
# Test pairwise within the groups with Wilcoxon Rank Sum Tests
pp<-pairwise.wilcox.test(data_in.richness.meta$Observed, data_in.richness.meta$ControlGroupWeek, 
                         p.adjust.method="bonferroni") # more conservative


PT = pp$p.value
PT1 = fullPTable(PT)
multcomp_observed <- multcompLetters(PT1,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)

# ++++++++++ (b) with  control
kruskal.test(Observed ~ GroupWithControl, data=data_in.richness.meta)
#Kruskal-Wallis chi-squared = 33.954, df = 3, p-value = 2.026e-07

# Wilcoxon-Mann-Whitney post hoc:==============================================
# Test pairwise within the groups with Wilcoxon Rank Sum Tests
pp<-pairwise.wilcox.test(data_in.richness.meta$Observed, data_in.richness.meta$GroupWithControl, 
                         p.adjust.method="bonferroni") # more conservative
PT = pp$p.value
PT1 = fullPTable(PT)
multcomp_observed_a <- multcompLetters(PT1,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)
 

############# Create plots ##########################

# ++++++++++ (a) without  control 
Observed_b<-ggplot(data_in.richness.meta, aes(x = factor(GroupWithControl), y = Observed,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape),alpha = 0.6,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  ylab("ASV Richness")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
        size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

# ++++++++++ (b) with  control 

Observed_a<-ggplot(data_in.richness.meta, aes(x = factor(ControlGroupWeek), y = Observed,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape), alpha=0.5,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  ylab("ASV Richness")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
        size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

# ++++++++++ (c) fall only
data_in.richness_Fall <- data_in.richness.meta[data_in.richness.meta$ControlGroupWeek=="Fall",]

Observed_fall<-ggplot(data_in.richness_Fall, aes(x = factor(GroupWithControl), y = PD,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape),alpha = 0.6,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values="#48b823")+
  scale_fill_manual(values="#48b823")+
  ylab("ASV Richness")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
        size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#=(2) Chao1 =======================================================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

############# Statitival analysis ##########################

# ++++++++++ (a) without  control 
#  we use Kruskal-Wallis (non-parametric equivalent of ANOVA)
kruskal.test(Chao1 ~ ControlGroupWeek, data=data_in.richness.meta)

# Wilcoxon-Mann-Whitney post hoc:==============================================
# Test pairwise within the groups with Wilcoxon Rank Sum Tests
pp<-pairwise.wilcox.test(data_in.richness.meta$Chao1, data_in.richness.meta$ControlGroupWeek, 
                         p.adjust.method="bonferroni") # more conservative


PT = pp$p.value
PT1 = fullPTable(PT)
multcomp_chao1 <- multcompLetters(PT1,
                                     compare="<",
                                     threshold=0.05,
                                     Letters=letters,
                                     reversed = FALSE)

# ++++++++++ (b) with  control
kruskal.test(Chao1 ~ GroupWithControl, data=data_in.richness.meta)
#Kruskal-Wallis chi-squared = 33.954, df = 3, p-value = 2.026e-07

# Wilcoxon-Mann-Whitney post hoc:==============================================
# Test pairwise within the groups with Wilcoxon Rank Sum Tests
pp<-pairwise.wilcox.test(data_in.richness.meta$Chao1, data_in.richness.meta$GroupWithControl, 
                         p.adjust.method="bonferroni") # more conservative
PT = pp$p.value
PT1 = fullPTable(PT)
multcomp_chao1_a <- multcompLetters(PT1,
                                       compare="<",
                                       threshold=0.05,
                                       Letters=letters,
                                       reversed = FALSE)


############# Create plots ##########################

# ++++++++++ (a) without  control 
Chao1_b<-ggplot(data_in.richness.meta, aes(x = factor(GroupWithControl), y = Chao1,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape),alpha = 0.6,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  ylab("Chao1")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

# ++++++++++ (b) with  control 

Chao1_a<-ggplot(data_in.richness.meta, aes(x = factor(ControlGroupWeek), y = Chao1,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape), alpha=0.5,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  ylab("Chao1")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

# ++++++++++ (c) fall only
kruskal.test(Chao1 ~ GroupWithControl, data=data_in.richness_Fall)
pairwise.wilcox.test(data_in.richness_Fall$Chao1, data_in.richness_Fall$GroupWithControl, 
                     p.adjust.method="bonferroni")

Chao1_fall<-ggplot(data_in.richness_Fall, aes(x = factor(GroupWithControl), y = Chao1,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape),alpha = 0.6,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values="#48b823")+
  scale_fill_manual(values="#48b823")+
  ylab("Chao1")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#=(3) Shannon =======================================================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

############# Statitival analysis ##########################

# ++++++++++ (a) without  control 
#  we use Kruskal-Wallis (non-parametric equivalent of ANOVA)
kruskal.test(Shannon ~ ControlGroupWeek, data=data_in.richness.meta)

# Wilcoxon-Mann-Whitney post hoc:==============================================
# Test pairwise within the groups with Wilcoxon Rank Sum Tests
pp<-pairwise.wilcox.test(data_in.richness.meta$Shannon, data_in.richness.meta$ControlGroupWeek, 
                         p.adjust.method="bonferroni") # more conservative


PT = pp$p.value
PT1 = fullPTable(PT)
multcomp_shannon <- multcompLetters(PT1,
                                  compare="<",
                                  threshold=0.05,
                                  Letters=letters,
                                  reversed = FALSE)

# ++++++++++ (b) with  control
kruskal.test(Shannon ~ GroupWithControl, data=data_in.richness.meta)
#Kruskal-Wallis chi-squared = 33.954, df = 3, p-value = 2.026e-07

# Wilcoxon-Mann-Whitney post hoc:==============================================
# Test pairwise within the groups with Wilcoxon Rank Sum Tests
pp<-pairwise.wilcox.test(data_in.richness.meta$Shannon, data_in.richness.meta$GroupWithControl, 
                         p.adjust.method="bonferroni") # more conservative
PT = pp$p.value
PT1 = fullPTable(PT)
multcomp_shannon_a <- multcompLetters(PT1,
                                    compare="<",
                                    threshold=0.05,
                                    Letters=letters,
                                    reversed = FALSE)


############# Create plots ##########################

# ++++++++++ (a) without  control 
Shannon_b<-ggplot(data_in.richness.meta, aes(x = factor(GroupWithControl), y = Shannon,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape),alpha = 0.6,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  ylab("Shannon")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
           size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

# ++++++++++ (b) with  control 

Shannon_a<-ggplot(data_in.richness.meta, aes(x = factor(ControlGroupWeek), y = Shannon,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape), alpha=0.5,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  ylab("Shannon")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
         size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

# ++++++++++ (c) fall only

Shannon_fall<-ggplot(data_in.richness_Fall, aes(x = factor(GroupWithControl), y = Shannon,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape),alpha = 0.6,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values="#48b823")+
  scale_fill_manual(values="#48b823")+
  ylab("Shannon")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
        size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

pp<-pairwise.wilcox.test(data_in.richness_Fall$Shannon, data_in.richness_Fall$GroupWithControl, 
                         p.adjust.method="bonferroni")
PT = pp$p.value
PT1 = fullPTable(PT)
multcompLetters(PT1,
                    compare="<",
                    threshold=0.05,
                    Letters=letters,
                    reversed = FALSE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#=(4) Simpson =======================================================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

############# Statitival analysis ##########################

# ++++++++++ (a) without  control 
#  we use Kruskal-Wallis (non-parametric equivalent of ANOVA)
kruskal.test(Simpson ~ ControlGroupWeek, data=data_in.richness.meta)

# Wilcoxon-Mann-Whitney post hoc:==============================================
# Test pairwise within the groups with Wilcoxon Rank Sum Tests
pp<-pairwise.wilcox.test(data_in.richness.meta$Simpson, data_in.richness.meta$ControlGroupWeek, 
                         p.adjust.method="bonferroni") # more conservative


PT = pp$p.value
PT1 = fullPTable(PT)
multcomp_simpson <- multcompLetters(PT1,
                                    compare="<",
                                    threshold=0.05,
                                    Letters=letters,
                                    reversed = FALSE)

# ++++++++++ (b) with  control
kruskal.test(Shannon ~ GroupWithControl, data=data_in.richness.meta)
#Kruskal-Wallis chi-squared = 33.954, df = 3, p-value = 2.026e-07

# Wilcoxon-Mann-Whitney post hoc:==============================================
# Test pairwise within the groups with Wilcoxon Rank Sum Tests
pp<-pairwise.wilcox.test(data_in.richness.meta$Simpson, data_in.richness.meta$GroupWithControl, 
                         p.adjust.method="bonferroni") # more conservative
PT = pp$p.value
PT1 = fullPTable(PT)
multcomp_simpson_a <- multcompLetters(PT1,
                                      compare="<",
                                      threshold=0.05,
                                      Letters=letters,
                                      reversed = FALSE)


############# Create plots ##########################

# ++++++++++ (a) without  control 
Simpson_b<-ggplot(data_in.richness.meta, aes(x = factor(GroupWithControl), y = Simpson,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape),alpha = 0.6,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  ylab("Simpson")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

# ++++++++++ (b) with  control 

Simpson_a<-ggplot(data_in.richness.meta, aes(x = factor(ControlGroupWeek), y = Simpson,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape), alpha=0.5,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  ylab("Simpson")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

# ++++++++++ (c) fall only
kruskal.test(Simpson ~ GroupWithControl, data=data_in.richness_Fall)

Simpson_fall<-ggplot(data_in.richness_Fall, aes(x = factor(GroupWithControl), y = Simpson,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape),alpha = 0.6,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values="#48b823")+
  scale_fill_manual(values="#48b823")+
  ylab("Shannon")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#=(5) Phylogenetic diversity =======================================================================================================
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

############# Statitival analysis ##########################

# ++++++++++ (a) without  control 
#  we use Kruskal-Wallis (non-parametric equivalent of ANOVA)
kruskal.test(PD ~ ControlGroupWeek, data=data_in.richness.meta)

# Wilcoxon-Mann-Whitney post hoc:==============================================
# Test pairwise within the groups with Wilcoxon Rank Sum Tests
pp<-pairwise.wilcox.test(data_in.richness.meta$PD, data_in.richness.meta$ControlGroupWeek, 
                         p.adjust.method="bonferroni") # more conservative


PT = pp$p.value
PT1 = fullPTable(PT)
multcomp_PD <- multcompLetters(PT1,
                                    compare="<",
                                    threshold=0.05,
                                    Letters=letters,
                                    reversed = FALSE)

# ++++++++++ (b) with  control
kruskal.test(PD ~ GroupWithControl, data=data_in.richness.meta)
#Kruskal-Wallis chi-squared = 33.954, df = 3, p-value = 2.026e-07

# Wilcoxon-Mann-Whitney post hoc:==============================================
# Test pairwise within the groups with Wilcoxon Rank Sum Tests
pp<-pairwise.wilcox.test(data_in.richness.meta$PD, data_in.richness.meta$GroupWithControl, 
                         p.adjust.method="bonferroni") # more conservative
PT = pp$p.value
PT1 = fullPTable(PT)
multcomp_PD_a <- multcompLetters(PT1,
                                      compare="<",
                                      threshold=0.05,
                                      Letters=letters,
                                      reversed = FALSE)


############# Create plots ##########################

# ++++++++++ (a) without  control 
PD_b<-ggplot(data_in.richness.meta, aes(x = factor(GroupWithControl), y = PD,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape),alpha = 0.6,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  ylab("PD")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

# ++++++++++ (b) with  control 

PD_a<-ggplot(data_in.richness.meta, aes(x = factor(ControlGroupWeek), y = PD,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape), alpha=0.5,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  ylab("PD")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

# ++++++++++ (c) fall only

kruskal.test(PD ~ GroupWithControl, data=data_in.richness_Fall)
pairwise.wilcox.test(data_in.richness_Fall$PD, data_in.richness_Fall$GroupWithControl, 
                     p.adjust.method="bonferroni")

PD_fall<-ggplot(data_in.richness_Fall, aes(x = factor(GroupWithControl), y = PD,color = ControlGroupWeek,fill = ControlGroupWeek)) +
  geom_boxplot(alpha = 0.4, outlier.shape=NA) +
  geom_point(aes(shape=Shape,size = Shape),alpha = 0.6,position=position_jitter(0.05))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values="#48b823")+
  scale_fill_manual(values="#48b823")+
  ylab("PD")+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"),legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())


## === Plot all

library(ggpubr)

# For figure 3
ggarrange(Chao1_a  + rremove("y.title"),
          Simpson_a + rremove("y.title"),
          PD_a  + rremove("y.title"),
          heights = c(5, 5),
          widths = c(2,2),
          labels = c("(a)", "(b)", "(c)"),
          ncol = 3, nrow = 1,
          align = "v")

ggarrange(Chao1_a  + rremove("y.title"),
          Simpson_a + rremove("y.title"),
          PD_a  + rremove("y.title"),
          heights = c(5,5,5),
          widths = c(3,3,3),
          ncol = 3, nrow = 1,
          align = "v")

ggarrange(Chao1_a  + rremove("y.title"),
          Shannon_a + rremove("y.title"),
          PD_a  + rremove("y.title"),
          heights = c(5,5,5),
          widths = c(3,3,3),
          ncol = 3, nrow = 1,
          align = "v")
# for supplementary 

ggarrange(Chao1_fall,Simpson_fall,PD_fall, 
          heights = c(3,3,3),
          widths = c(6,6,6),
          labels = c("(a)", "(b)", "(c)"),
          ncol = 1, nrow = 3,
          align = "v")


ggarrange(Observed_a,
          Chao1_a,
          Simpson_a, 
          PD_a,
          heights = c(3, 3),
          widths = c(2,2),
          labels = c("(a)", "(b)", "(c)","(d)"),
          ncol = 2, nrow = 2,
          align = "v")


ggarrange(Observed_b,
          Chao1_b,
          Simpson_b, 
          PD_b,
          heights = c(3, 3),
          widths = c(2,2),
          labels = c("(a)", "(b)", "(c)","(d)"),
          ncol = 2, nrow = 2,
          align = "v")



## make table of all statistics
library(dplyr)
WilcT_Simpson<-as.data.frame(multcomp_simpson_a$monospacedLetters)
WilcT_Simpson<-rename(WilcT_Simpson,Simpson = "multcomp_simpson_a$monospacedLetters")
WilcT_Shannon<-as.data.frame(multcomp_shannon_a$monospacedLetters)
WilcT_Shannon<-rename(WilcT_Shannon, Shannon="multcomp_shannon_a$monospacedLetters")
WilcT_Observed<-as.data.frame(multcomp_observed_a$monospacedLetters)
WilcT_Observed<-rename(WilcT_Observed, Observed="multcomp_observed_a$monospacedLetters")
WilcT_Caho1<-as.data.frame(multcomp_chao1_a$monospacedLetters)
WilcT_Caho1<-rename(WilcT_Caho1, Caho1="multcomp_chao1_a$monospacedLetters")
WilcT_PD<-as.data.frame(multcomp_PD_a$monospacedLetters)
WilcT_PD<-rename(WilcT_PD, PD="multcomp_PD_a$monospacedLetters")

wilcox_test_control <- cbind(WilcT_Simpson,WilcT_Shannon,WilcT_Observed,WilcT_Caho1,WilcT_PD)


WilcT_Simpson<-as.data.frame(multcomp_simpson$monospacedLetters)
WilcT_Simpson<-rename(WilcT_Simpson,Simpson = "multcomp_simpson$monospacedLetters")
WilcT_Shannon<-as.data.frame(multcomp_shannon$monospacedLetters)
WilcT_Shannon<-rename(WilcT_Shannon, Shannon="multcomp_shannon$monospacedLetters")
WilcT_Observed<-as.data.frame(multcomp_observed$monospacedLetters)
WilcT_Observed<-rename(WilcT_Observed, Observed="multcomp_observed$monospacedLetters")
WilcT_Caho1<-as.data.frame(multcomp_chao1$monospacedLetters)
WilcT_Caho1<-rename(WilcT_Caho1, Caho1="multcomp_chao1$monospacedLetters")
WilcT_PD<-as.data.frame(multcomp_PD$monospacedLetters)
WilcT_PD<-rename(WilcT_PD, PD="multcomp_PD$monospacedLetters")

wilcox_test <- cbind(WilcT_Simpson,WilcT_Shannon,WilcT_Observed,WilcT_Caho1,WilcT_PD)
