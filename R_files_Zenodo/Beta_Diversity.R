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
library(ggpubr)
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
############################ Beta - diversity ######################################################
#=(1)==============================================================================================================
#++++++++++++++++++ Wieghted UniFrac Beta-Diversity PLOT  Data +++++++++++++++++++++++++++++++++++++++++++++++++++++
#There is no diference if we choose method = "MDS" or method = "PCoA" - IT IS THE SAME
ord.bray.MDS <- phyloseq::ordinate(data_in, method = "PCoA", distance = "wunifrac")

#Scree plo
plot_scree(ord.bray.MDS, "Scree plot for Global Patterns, Wieghted UniFrac/PCoA")
#ordination plot
ord.plot.bray.MDS_12 <- plot_ordination(data_in, ord.bray.MDS, color = "ControlGroupWeek",axes = 1:2)
ord.plot.bray.MDS <- ord.plot.bray.MDS_12 + theme_bw() + ggtitle("Wieghted UniFrac Distance - PCoA")  +
  geom_point(size=4)+
  scale_color_manual(values=Col1)+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot(ord.plot.bray.MDS)

#------------CHECK HOMOGENEITY OF VARIANCE!!------------------------------
# see here: https://microbiome.github.io/tutorials/PERMANOVA.html
# and gere, for defult plotting options: https://www.fromthebottomoftheheap.net/2016/04/17/new-plot-default-for-betadisper/

metadata.BRAY.multi.anova.fall.data_in <- as(sample_data(data_in), "data.frame")
# Calculate distance matrices: Bray
data_in.bray <-phyloseq::distance(data_in, method = "wunifrac")
# Calculate the diferences in dispersion
ps.disper.bray <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$GroupWithControl)
permutest(ps.disper.bray, pairwise = TRUE)


#- make Plots---------------------------------------
# only run for the labels, we will re-run for the plot below 
ps.disper.bray <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$ControlGroupWeek) 
Centroids<-data.frame(ps.disper.bray[["centroids"]])
Vectors<-data.frame(ps.disper.bray[["vectors"]])
Distances<-data.frame(ps.disper.bray[["distances"]])
Lables<-data.frame(ps.disper.bray[["group"]])
LablesControl<-data.frame(ps.disper.brayWithControl[["group"]])
p.compHomogenity<-permutest(ps.disper.bray, pairwise = TRUE)
# make data-frame for boxplt
df <- data.frame(Distances,Lables)
names(df) <- c("Distances","Lables")

# make data-frame Ordination plot
ps.disper.brayWithControl <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$GroupWithControl)
Shape<-0
Shape[LablesControl=="Control_Hula"]<-"ControlHula"
Shape[LablesControl=="Control_Other"]<-"ControlOther"
Shape[LablesControl!="Control_Hula" & LablesControl!="Control_Other"]<-"data"
Shape<-as.factor(Shape)
dfOrd <-data.frame(Vectors$PCoA1,Vectors$PCoA2,Lables,LablesControl,Shape)
names(dfOrd) <- c("PCoA1","PCoA2","Lables","LablesControl","Shape")

# create the confidence interval ellipces
library(ellipse)
df_ell <- data.frame()
for(g in levels(dfOrd$Lables)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(dfOrd[dfOrd$Lables==g,], 
                      ellipse(cor(PCoA1, PCoA2), scale=c(sd(PCoA1),sd(PCoA2)),
                            centre=c(mean(PCoA1),mean(PCoA2))))),Lables=g))
}

# create centroid data-frame
dCent<-data.frame(Centroids$PCoA1,Centroids$PCoA2,rownames(Centroids))
names(dCent) <- c("PCoA1","PCoA2","Lables")


# Plot the ordination with centroids 
ggplot(dfOrd, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(colour = as.factor(LablesControl)),size = 3, shape=16)+
  scale_color_manual(values=Col2)+
  geom_path(data=df_ell, aes(x=x, y=y,colour=Lables), size=1, linetype=2)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))
OrdUnifrac <- ggplot(dfOrd, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill = Lables,color = Lables,shape=Shape,size = Shape))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  geom_path(data=df_ell, aes(x=x, y=y,colour=Lables), size=1, linetype=2)+
  geom_point(data=dCent,aes(x=PCoA1, y=PCoA2,colour = Lables),size = 5, shape=13,stroke = 2)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))+
  theme(legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())
# Boxplot of distances
dispersion_Unifrac <- ggplot(df, aes(x=Lables, y=Distances, color = Lables,fill = Lables))+
  geom_boxplot(alpha = 0.2, outlier.shape=NA)+
  geom_point(alpha = 0.6, size=3, position=position_jitter(0.05),aes(color=Lables))+
  scale_color_manual(values = Col1)+
  scale_fill_manual(values =Col1)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))+
theme(legend.position = "none",axis.text.x=element_blank(),
      axis.title.x=element_blank())

# ---Are the groups diferenet---------------------
# Perform a multivariate ANOVA, all the groups
adonis(data_in.bray ~ GroupWithControl, data = metadata.BRAY.multi.anova.fall.data_in, perm = 9999)

# use the pairwise comparison function dounloaded from here: https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
source("pairwise.adonis.R") # connect the function to this script
pairwiseControlGroup<- pairwise.adonis(data_in.bray,metadata.BRAY.multi.anova.fall.data_in$GroupWithControl
                                       , p.adjust("bonferroni"),perm=9999)

##  ---- for the main groups, without controls--------------------------------------
library(multcompView) # convert this table to a compact letter display
library(rcompanion)  # convert to table to use multcompLetters
# Make letters
pairwiseControlGroup<- pairwise.adonis(data_in.bray,metadata.BRAY.multi.anova.fall.data_in$ControlGroupWeek,
                                       p.adjust("bonferroni"),perm=9999)
AB<-pairwiseControlGroup[,6]
names(AB) <-c("PreMig-Fall", "PreMig-Winter", "PreMig-FeedSt", 
              "Fall-Winter","Fall-FeedSt","Winter-FeedSt")
multcompLetters(AB,compare="<",
                threshold=0.001,
                Letters=letters,
                reversed = FALSE)

#=(2)==============================================================================================================
#++++++++++++++++++ Un-Wieghted UniFrac Beta-Diversity PLOT  Data +++++++++++++++++++++++++++++++++++++++++++++++++++++
#There is no diference if we choose method = "MDS" or method = "PCoA" - IT IS THE SAME
ord.bray.MDS <- phyloseq::ordinate(data_in, method = "PCoA", distance = "uunifrac")

#Scree plo
plot_scree(ord.bray.MDS, "Scree plot for Global Patterns, Unweighted UniFrac/PCoA")
#ordination plot
ord.plot.bray.MDS_12 <- plot_ordination(data_in, ord.bray.MDS, color = "ControlGroupWeek",axes = 1:2)
ord.plot.bray.MDS <- ord.plot.bray.MDS_12 + theme_bw() + ggtitle("Unweighted UniFrac Distance - PCoA")  +
  geom_point(size=4)+
  scale_color_manual(values=Col1)+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot(ord.plot.bray.MDS)

#------------CHECK HOMOGENEITY OF VARIANCE!!------------------------------

metadata.BRAY.multi.anova.fall.data_in <- as(sample_data(data_in), "data.frame")
# Calculate distance matrices: Bray
data_in.bray <-phyloseq::distance(data_in, method = "uunifrac")
# Calculate the diferences in dispersion
ps.disper.bray <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$GroupWithControl)
permutest(ps.disper.bray, pairwise = TRUE)


#- make Plots---------------------------------------
# only run for the labels, we will re-run for the plot below 
ps.disper.bray <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$ControlGroupWeek) 
Centroids<-data.frame(ps.disper.bray[["centroids"]])
Vectors<-data.frame(ps.disper.bray[["vectors"]])
Distances<-data.frame(ps.disper.bray[["distances"]])
Lables<-data.frame(ps.disper.bray[["group"]])
LablesControl<-data.frame(ps.disper.brayWithControl[["group"]])
p.compHomogenity<-permutest(ps.disper.bray, pairwise = TRUE)
# make data-frame for boxplt
df <- data.frame(Distances,Lables)
names(df) <- c("Distances","Lables")

# make data-frame Ordination plot
ps.disper.brayWithControl <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$GroupWithControl)
Shape<-0
Shape[LablesControl=="Control_Hula"]<-"ControlHula"
Shape[LablesControl=="Control_Other"]<-"ControlOther"
Shape[LablesControl!="Control_Hula" & LablesControl!="Control_Other"]<-"data"
Shape<-as.factor(Shape)
dfOrd <-data.frame(Vectors$PCoA1,Vectors$PCoA2,Lables,LablesControl,Shape)
names(dfOrd) <- c("PCoA1","PCoA2","Lables","LablesControl","Shape")

# create the confidence interval ellipces
library(ellipse)
df_ell <- data.frame()
for(g in levels(dfOrd$Lables)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(dfOrd[dfOrd$Lables==g,], 
                                                   ellipse(cor(PCoA1, PCoA2), scale=c(sd(PCoA1),sd(PCoA2)),
                                                           centre=c(mean(PCoA1),mean(PCoA2))))),Lables=g))
}

# create centroid data-frame
dCent<-data.frame(Centroids$PCoA1,Centroids$PCoA2,rownames(Centroids))
names(dCent) <- c("PCoA1","PCoA2","Lables")


# Plot the ordination with centroids 
ggplot(dfOrd, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(colour = Lables),size = 3, shape=16)+
  scale_color_manual(values=c("tomato1","#559bab","#003f5c","red4"))+
  geom_path(data=df_ell, aes(x=x, y=y,colour=Lables), size=1, linetype=2)+
  geom_point(data=dCent,aes(x=PCoA1, y=PCoA2,colour = Lables),size = 5, shape=13,stroke = 2)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))
OrdUunifrac <- ggplot(dfOrd, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill = Lables,color = Lables,shape=Shape,size = Shape))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  geom_path(data=df_ell, aes(x=x, y=y,colour=Lables), size=1, linetype=2)+
  geom_point(data=dCent,aes(x=PCoA1, y=PCoA2,colour = Lables),size = 5, shape=13,stroke = 2)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))+
  theme(legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())
# Boxplot of distances
dispersion_Uunifrac <- ggplot(df, aes(x=Lables, y=Distances, color = Lables,fill = Lables))+
  geom_boxplot(alpha = 0.2, outlier.shape=NA)+
  geom_point(alpha = 0.6, size=3, position=position_jitter(0.05),aes(color=Lables))+
  scale_color_manual(values = Col1)+
  scale_fill_manual(values =Col1)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))+
  theme(legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())


# ---Are the groups diferenet---------------------
# Perform a multivariate ANOVA, all the groups
adonis(data_in.bray ~ GroupWithControl, data = metadata.BRAY.multi.anova.fall.data_in, perm = 9999)

# use the pairwise comparison function dounloaded from here: https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
source("pairwise.adonis.R") # connect the function to this script
pairwiseControlGroup<- pairwise.adonis(data_in.bray,metadata.BRAY.multi.anova.fall.data_in$GroupWithControl
                                       , p.adjust("bonferroni"),perm=9999)
##  ---- for the main groups, without controls--------------------------------------
library(multcompView) # convert this table to a compact letter display
library(rcompanion)  # convert to table to use multcompLetters
# Make letters
pairwiseControlGroup<- pairwise.adonis(data_in.bray,metadata.BRAY.multi.anova.fall.data_in$ControlGroupWeek,
                                       p.adjust("bonferroni"),perm=9999)
AB<-pairwiseControlGroup[,6]
names(AB) <-c("PreMig-Fall", "PreMig-Winter", "PreMig-FeedSt", 
              "Fall-Winter","Fall-FeedSt","Winter-FeedSt")
multcompLetters(AB,compare="<",
                threshold=0.001,
                Letters=letters,
                reversed = FALSE)

#=(3)==============================================================================================================
#++++++++++++++++++ Jaccard Beta-Diversity PLOT  Data +++++++++++++++++++++++++++++++++++++++++++++++++++++
#There is no diference if we choose method = "MDS" or method = "PCoA" - IT IS THE SAME
ord.bray.MDS <- phyloseq::ordinate(data_in, method = "PCoA",distance = "jaccard", binary = TRUE)

#Scree plo
plot_scree(ord.bray.MDS, "Scree plot for Global Patterns, jaccard/PCoA")
#ordination plot
ord.plot.bray.MDS_12 <- plot_ordination(data_in, ord.bray.MDS, color = "ControlGroupWeek",axes = 1:2)
ord.plot.bray.MDS <- ord.plot.bray.MDS_12 + theme_bw() + ggtitle("jaccard Distance - PCoA")  +
  geom_point(size=4)+
  scale_color_manual(values=Col1)+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot(ord.plot.bray.MDS)

#------------CHECK HOMOGENEITY OF VARIANCE!!------------------------------

metadata.BRAY.multi.anova.fall.data_in <- as(sample_data(data_in), "data.frame")
# Calculate distance matrices: Bray
data_in.bray <-phyloseq::distance(data_in,  method = "jaccard", binary = TRUE)
# Calculate the diferences in dispersion
ps.disper.bray <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$GroupWithControl)
permutest(ps.disper.bray, pairwise = TRUE)


#- make Plots---------------------------------------
# only run for the labels, we will re-run for the plot below 
ps.disper.bray <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$ControlGroupWeek) 
Centroids<-data.frame(ps.disper.bray[["centroids"]])
Vectors<-data.frame(ps.disper.bray[["vectors"]])
Distances<-data.frame(ps.disper.bray[["distances"]])
Lables<-data.frame(ps.disper.bray[["group"]])
LablesControl<-data.frame(ps.disper.brayWithControl[["group"]])
p.compHomogenity<-permutest(ps.disper.bray, pairwise = TRUE)
# make data-frame for boxplt
df <- data.frame(Distances,Lables)
names(df) <- c("Distances","Lables")

# make data-frame Ordination plot
ps.disper.brayWithControl <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$GroupWithControl)
Shape<-0
Shape[LablesControl=="Control_Hula"]<-"ControlHula"
Shape[LablesControl=="Control_Other"]<-"ControlOther"
Shape[LablesControl!="Control_Hula" & LablesControl!="Control_Other"]<-"data"
Shape<-as.factor(Shape)
dfOrd <-data.frame(Vectors$PCoA1,Vectors$PCoA2,Lables,LablesControl,Shape)
names(dfOrd) <- c("PCoA1","PCoA2","Lables","LablesControl","Shape")

# create the confidence interval ellipces
library(ellipse)
df_ell <- data.frame()
for(g in levels(dfOrd$Lables)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(dfOrd[dfOrd$Lables==g,], 
                                                   ellipse(cor(PCoA1, PCoA2), scale=c(sd(PCoA1),sd(PCoA2)),
                                                           centre=c(mean(PCoA1),mean(PCoA2))))),Lables=g))
}

# create centroid data-frame
dCent<-data.frame(Centroids$PCoA1,Centroids$PCoA2,rownames(Centroids))
names(dCent) <- c("PCoA1","PCoA2","Lables")


# Plot the ordination with centroids 
ggplot(dfOrd, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(colour = Lables),size = 3, shape=16)+
  scale_color_manual(values=c("tomato1","#559bab","#003f5c","red4"))+
  geom_path(data=df_ell, aes(x=x, y=y,colour=Lables), size=1, linetype=2)+
  geom_point(data=dCent,aes(x=PCoA1, y=PCoA2,colour = Lables),size = 5, shape=13,stroke = 2)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))
OrdJaccard <- ggplot(dfOrd, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill = Lables,color = Lables,shape=Shape,size = Shape))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  geom_path(data=df_ell, aes(x=x, y=y,colour=Lables), size=1, linetype=2)+
  geom_point(data=dCent,aes(x=PCoA1, y=PCoA2,colour = Lables),size = 5, shape=13,stroke = 2)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))+
  theme(legend.position = "none",
        axis.title.x=element_blank(),axis.title.y=element_blank())
# Boxplot of distances
dispersion_Jaccard <- ggplot(df, aes(x=Lables, y=Distances, color = Lables,fill = Lables))+
  geom_boxplot(alpha = 0.2, outlier.shape=NA)+
  geom_point(alpha = 0.6, size=3, position=position_jitter(0.05),aes(color=Lables))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))+
  theme(legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

lev <- levels(df$Lables) # get the variables
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
dispersion_Jaccard <- dispersion_Jaccard + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.01, 0.05, 1),
    symbols = c("**", "*", "n.s")
  )
)

# ---Are the groups diferenet---------------------
# Perform a multivariate ANOVA, all the groups
adonis(data_in.bray ~ GroupWithControl, data = metadata.BRAY.multi.anova.fall.data_in, perm = 9999)

# use the pairwise comparison function dounloaded from here: https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
source("pairwise.adonis.R") # connect the function to this script
pairwiseControlGroup<- pairwise.adonis(data_in.bray,metadata.BRAY.multi.anova.fall.data_in$GroupWithControl
                                       , p.adjust("bonferroni"),perm=9999)
##  ---- for the main groups, without controls--------------------------------------
library(multcompView) # convert this table to a compact letter display
library(rcompanion)  # convert to table to use multcompLetters
# Make letters
pairwiseControlGroup<- pairwise.adonis(data_in.bray,metadata.BRAY.multi.anova.fall.data_in$ControlGroupWeek,
                                       p.adjust("bonferroni"),perm=9999)
AB<-pairwiseControlGroup[,6]
names(AB) <-c("PreMig-Fall", "PreMig-Winter", "PreMig-FeedSt", 
              "Fall-Winter","Fall-FeedSt","Winter-FeedSt")
multcompLetters(AB,compare="<",
                threshold=0.001,
                Letters=letters,
                reversed = FALSE)
#=(4)==============================================================================================================
#++++++++++++++++++ Bray-Curtis Beta-Diversity PLOT  Data +++++++++++++++++++++++++++++++++++++++++++++++++++++
#There is no diference if we choose method = "MDS" or method = "PCoA" - IT IS THE SAME
ord.bray.MDS <- phyloseq::ordinate(data_in, method = "PCoA", distance = "bray")

#Scree plo
plot_scree(ord.bray.MDS, "Scree plot for Global Patterns, Bray-Curtis/PCoA")
#ordination plot
ord.plot.bray.MDS_12 <- plot_ordination(data_in, ord.bray.MDS, color = "ControlGroupWeek",axes = 1:2)
ord.plot.bray.MDS <- ord.plot.bray.MDS_12 + theme_bw() + ggtitle("Bray-Curtis Distance - PCoA")  +
  geom_point(size=4)+
  scale_color_manual(values=Col1)+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
plot(ord.plot.bray.MDS)

#------------CHECK HOMOGENEITY OF VARIANCE!!------------------------------

metadata.BRAY.multi.anova.fall.data_in <- as(sample_data(data_in), "data.frame")
# Calculate distance matrices: Bray
data_in.bray <-phyloseq::distance(data_in, method = "bray")
# Calculate the diferences in dispersion
ps.disper.bray <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$GroupWithControl)
permutest(ps.disper.bray, pairwise = TRUE)

#- make Plots---------------------------------------
# only run for the labels, we will re-run for the plot below 
ps.disper.bray <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$ControlGroupWeek) 
Centroids<-data.frame(ps.disper.bray[["centroids"]])
Vectors<-data.frame(ps.disper.bray[["vectors"]])
Distances<-data.frame(ps.disper.bray[["distances"]])
Lables<-data.frame(ps.disper.bray[["group"]])
LablesControl<-data.frame(ps.disper.brayWithControl[["group"]])
p.compHomogenity<-permutest(ps.disper.bray, pairwise = TRUE)
# make data-frame for boxplt
df <- data.frame(Distances,Lables)
names(df) <- c("Distances","Lables")

# make data-frame Ordination plot
ps.disper.brayWithControl <- betadisper(data_in.bray, metadata.BRAY.multi.anova.fall.data_in$GroupWithControl)
Shape<-0
Shape[LablesControl=="Control_Hula"]<-"ControlHula"
Shape[LablesControl=="Control_Other"]<-"ControlOther"
Shape[LablesControl!="Control_Hula" & LablesControl!="Control_Other"]<-"data"
Shape<-as.factor(Shape)
dfOrd <-data.frame(Vectors$PCoA1,Vectors$PCoA2,Lables,LablesControl,Shape)
names(dfOrd) <- c("PCoA1","PCoA2","Lables","LablesControl","Shape")

# create the confidence interval ellipces
library(ellipse)
df_ell <- data.frame()
for(g in levels(dfOrd$Lables)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(dfOrd[dfOrd$Lables==g,], 
                                                   ellipse(cor(PCoA1, PCoA2), scale=c(sd(PCoA1),sd(PCoA2)),
                                                           centre=c(mean(PCoA1),mean(PCoA2))))),Lables=g))
}

# create centroid data-frame
dCent<-data.frame(Centroids$PCoA1,Centroids$PCoA2,rownames(Centroids))
names(dCent) <- c("PCoA1","PCoA2","Lables")


# Plot the ordination with centroids 
ggplot(dfOrd, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(colour = Lables),size = 3, shape=16)+
  scale_color_manual(values=c("tomato1","#559bab","#003f5c","red4"))+
  geom_path(data=df_ell, aes(x=x, y=y,colour=Lables), size=1, linetype=2)+
  geom_point(data=dCent,aes(x=PCoA1, y=PCoA2,colour = Lables),size = 5, shape=13,stroke = 2)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))
OrdBray <- ggplot(dfOrd, aes(x=PCoA1, y=PCoA2))+
  geom_point(aes(fill = Lables,color = Lables,shape=Shape,size = Shape))+
  scale_shape_manual(values=c(1, 2, 16))+
  scale_size_manual(values=c(3, 2, 3))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  geom_path(data=df_ell, aes(x=x, y=y,colour=Lables), size=1, linetype=2)+
  geom_point(data=dCent,aes(x=PCoA1, y=PCoA2,colour = Lables),size = 5, shape=13,stroke = 2)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))+
  theme(legend.position = "none",
        axis.title.x=element_blank(),axis.title.y=element_blank())
# Boxplot of distances
dispersion_Bray <- ggplot(df, aes(x=Lables, y=Distances, color = Lables,fill = Lables))+
  geom_boxplot(alpha = 0.2, outlier.shape=NA)+
  geom_point(alpha = 0.6, size=3, position=position_jitter(0.05),aes(color=Lables))+
  scale_color_manual(values=Col1)+
  scale_fill_manual(values=Col1)+
  theme(text = element_text(size=18),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey", 
                                                                    size = 0.5, linetype = "solid"))+
  theme(legend.position = "none",axis.text.x=element_blank(),
        axis.title.x=element_blank())

lev <- levels(df$Lables) # get the variables
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
dispersion_Bray  <- dispersion_Bray  + stat_compare_means(
  comparisons = L.pairs,
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.01, 0.05, 1),
    symbols = c("**", "*", "n.s")
  )
)

# ---Are the groups diferenet---------------------
# Perform a multivariate ANOVA, all the groups
adonis(data_in.bray ~ GroupWithControl, data = metadata.BRAY.multi.anova.fall.data_in, perm = 9999)


# use the pairwise comparison function dounloaded from here: https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
source("pairwise.adonis.R") # connect the function to this script
pairwiseControlGroup<- pairwise.adonis(data_in.bray,metadata.BRAY.multi.anova.fall.data_in$GroupWithControl
                                       , p.adjust("bonferroni"),perm=9999)
##  ---- for the main groups, without controls--------------------------------------
library(multcompView) # convert this table to a compact letter display
library(rcompanion)  # convert to table to use multcompLetters
# Make letters
pairwiseControlGroup<- pairwise.adonis(data_in.bray,metadata.BRAY.multi.anova.fall.data_in$ControlGroupWeek,
                                       p.adjust("bonferroni"),perm=9999)
AB<-pairwiseControlGroup[,6]
names(AB) <-c("PreMig-Fall", "PreMig-Winter", "PreMig-FeedSt", 
              "Fall-Winter","Fall-FeedSt","Winter-FeedSt")
multcompLetters(AB,compare="<",
                threshold=0.001,
                Letters=letters,
                reversed = FALSE)

## === Plot all

library(ggpubr)

ggarrange(OrdBray,OrdJaccard,OrdUnifrac,OrdUunifrac, 
          heights = c(3, 3),
          widths = c(2,2),
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2,
          align = "v")

ggarrange(dispersion_Bray,dispersion_Jaccard,dispersion_Unifrac,dispersion_Uunifrac, 
          heights = c(3, 3,3,3),
          widths = c(2,2,2,2),
          labels = c("(a)", "(b)", "(c)","(d)"),
          ncol = 2, nrow = 2,
          align = "v")
