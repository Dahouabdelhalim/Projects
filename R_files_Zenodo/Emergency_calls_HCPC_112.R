library(devtools)
devtools::install_github("XanderHorn/autoEDA")
library(autoEDA)
library(ggplot2)
library(FactoMineR)
library(Factoshiny)
library(missMDA)
library(FactoInvestigate)
library(factoextra)
library(autoEDA)
library(cowplot) 
library(plotly)


setwd(' ') 
data.112 <- read.csv('dataset112_4Apr2021_PCA_1_3.csv', row.names=1, header=TRUE, sep=",")
#str(data.112)

###---------------------------------------------------------------------------
### o prima analiza e pentru clustere pe baza apelurilor grupate 
### pe speciiecele sase specii a doua analiza e pe baza categorilor de apeluri
###---------------------------------------------------------------------------


#----------------- HCPC pentru cele 6 specii de fauna
data.active <- data.112[1:318,5:10]
# Compute PCA with ncp = 3
res.pca <- PCA(data.active, ncp = 3, graph = TRUE)
# Compute hierarchical clustering on principal components
# HCPCshiny(res.pca)
res.hcpc <- HCPC(res.pca,nb.clust=4, graph = TRUE)

#------ grafice
fviz_dend(res.hcpc, main = "Cluster dendrogram PCA",
          cex = 0.2,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)
ggsave("PCA_dendrogram_map1_1_specii.tiff", width=40, height=20)
fviz_cluster(res.hcpc, geom = c("point", "text"), pointsize=0.8, labelsize= 16,
             font.x=c(20, "bold","black"), font.y=c(20, "bold","black"),
             legend=c("bottom"), legend.title="Cluster", font.legend=c(30, "plain","black"),
             font.tickslab=c(18,"plain","black"),
             repel = TRUE,                      # Avoid label overlapping
             show.clust.cent = TRUE,            # Show cluster centers
             palette = "lancet",                # Color palette see ?ggpubr::ggpar
             ggtheme = theme_gray(),
            )
ggsave("PCA_cluster_map1_1_specii.tiff", width=20, height=15)

#The function HCPC() returns a list containing:
head(res.hcpc$data.clust,318) #The original data with a supplementary column called class containing the partition.
res.hcpc$desc.var$quanti      #desc.var: The variables describing clusters
res.hcpc$desc.axes$quanti     #desc.ind: The more typical individuals of each cluster
res.hcpc$desc.ind$para        #desc.axes: The axes describing clusters


### cod realizare ploturi pentru clustere - preluat din https://nextjournal.com/pc-methods/hierarchical-clustering-pcs
### Vizualizare comparare clustere
names(res.hcpc$data.clust)[names(res.hcpc$data.clust) == 'clust'] <- 'Cluster'
autoEDA_results <- autoEDA(res.hcpc$data.clust, 
                           y = "Cluster", returnPlotList = TRUE,
                           outcomeType = "automatic", removeConstant = TRUE, 
                           removeZeroSpread = TRUE, removeMajorityMissing = TRUE, 
                           imputeMissing = TRUE, clipOutliers = FALSE, 
                           minLevelPercentage = 0.025, predictivePower = TRUE, 
                           outlierMethod = "tukey", lowPercentile = 0.05, 
                           upPercentile = 0.95, plotCategorical = "groupedBar", 
                           plotContinuous = "histogram", bins = 30, 
                           rotateLabels = TRUE, color = "#26A69A", 
                           verbose = FALSE) 
p <- plot_grid(plotlist = autoEDA_results$plots, ncol = 3)
tiff("grid PCA specii.tiff", width=2400, height=1200)
print(p)
dev.off()

### Predictive power cluster
res <- autoEDA_results$overview[autoEDA_results$overview$PredictivePower %in% c("High", "Medium", "Low"),]
res[, c('Feature', 'PredictivePowerPercentage', 'PredictivePower')]


#-------------------------------HCPC pentru cele 4 categori de apeluri HD,WD, RA si O
data.activec <- data.112[1:318,1:4]
# Compute PCA with ncp = 3
res.pca <- PCA(data.activec, ncp = 3, graph = TRUE)
# Compute hierarchical clustering on principal components
res.hcpc <- HCPC(res.pca,nb.clust=4, graph = FALSE)

#-------grafice
fviz_dend(res.hcpc, main = "Cluster dendrogram PCA",
          cex = 0.2, repel=TRUE,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)
ggsave("PCA_dendrogram_map1_1_categorii.tiff", width=30, height=15)
fviz_cluster(res.hcpc, geom = c("point", "text"), pointsize=0.8, labelsize= 16,
             font.x=c(20, "bold","black"), font.y=c(20, "bold","black"),
             legend=c("bottom"), legend.title="Cluster", font.legend=c(30, "plain","black"),
             font.tickslab=c(18,"plain","black"),
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE,            # Show cluster centers
             palette = "lancet",              # Color palette see ?ggpubr::ggpar
             ggtheme = theme_gray(),
)
ggsave("PCA_cluster_map1_1_categorii.tiff", width=20, height=15)

#The function HCPC() returns a list containing:
head(res.hcpc$data.clust,318) #The original data with a supplementary column called class containing the partition.
res.hcpc$desc.var$quanti      #desc.var: The variables describing clusters
res.hcpc$desc.axes$quanti     #desc.ind: The more typical individuals of each cluster
res.hcpc$desc.ind$para        #desc.axes: The axes describing clusters

### cod realizare ploturi pentru clustere - preluat din https://nextjournal.com/pc-methods/hierarchical-clustering-pcs
### Vizualizare comparare clustere
names(res.hcpc$data.clust)[names(res.hcpc$data.clust) == 'clust'] <- 'Cluster'
autoEDA_results <- autoEDA(res.hcpc$data.clust, 
                           y = "Cluster", returnPlotList = TRUE,
                           outcomeType = "automatic", removeConstant = TRUE, 
                           removeZeroSpread = TRUE, removeMajorityMissing = TRUE, 
                           imputeMissing = TRUE, clipOutliers = FALSE, 
                           minLevelPercentage = 0.025, predictivePower = TRUE, 
                           outlierMethod = "tukey", lowPercentile = 0.05, 
                           upPercentile = 0.95, plotCategorical = "groupedBar", 
                           plotContinuous = "histogram", bins = 30, 
                           rotateLabels = TRUE, color = "#26A69A", 
                           verbose = FALSE) 
p <- plot_grid(plotlist = autoEDA_results$plots, ncol = 3)
tiff("grid PCA categorii.tiff", width=2400, height=1200)
print(p)
dev.off()

### Predictive power cluster
res <- autoEDA_results$overview[autoEDA_results$overview$PredictivePower %in% c("High", "Medium", "Low"),]
res[, c('Feature', 'PredictivePowerPercentage', 'PredictivePower')]



# ------------------------------------- END ------------------------------






