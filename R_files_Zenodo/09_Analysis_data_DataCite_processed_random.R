
### Computational script-------------------------------------------------------------------------------------------------------------------------------------

# This script aims to analyze the data in the file "Data_DataCite_processed_random.txt".


## Load libraries--------------------------------------------------------------------------------------------------------------------------------------------

.libPaths("C:/R/library") # Modify the path to where the R libraries are available. If the packages are not yet installed, please install them first.

library('ggplot2')
library('vcd')
library('tidyr')
library('gmodels')
library('pheatmap')
library('cluster')
library('fpc')


## Load data------------------------------------------------------------------------------------------------------------------------------------------------

data_repo_final<-read.table(file=file.choose(), header=TRUE, sep="\\t", comment.char="", fill=TRUE, na.strings="character", quote="", row.names=NULL, stringsAsFactors = FALSE, strip.white = TRUE, encoding = "utf-8", blank.lines.skip = TRUE) # cf. the file "Data_DataCite_processed_random.txt"


## Bivariate analysis----------------------------------------------------------------------------------------------------------------------------------------

TABLE_1<-table(data_repo_final$Client_ID, data_repo_final$Identifiers_data_creators, dnn = c("Client_ID", "Identifiers_data_creators"))
TABLE_2<-table(data_repo_final$Client_ID, data_repo_final$Affiliations_data_creators, dnn = c("Client_ID", "Affiliations_data_creators"))
TABLE_3<-table(data_repo_final$Client_ID, data_repo_final$Identifiers_related_outputs, dnn = c("Client_ID", "Identifiers_related_outputs"))

mosaicplot(TABLE_1, shade = TRUE) # mosaic plot
assocstats(TABLE_1)
CrossTable(TABLE_1, chisq = TRUE, resid = TRUE, sresid = TRUE, digits = 2, format="SPSS")

mosaicplot(TABLE_2, shade = TRUE) # mosaic plot
assocstats(TABLE_2)
CrossTable(TABLE_2, chisq = TRUE, resid = TRUE, sresid = TRUE, digits = 2, format="SPSS")

mosaicplot(TABLE_3, shade = TRUE, main = "") # mosaic plot
assocstats(TABLE_3)
CrossTable(TABLE_3, chisq = TRUE, resid = TRUE, sresid = TRUE, digits = 2, format="SPSS")

jpeg(file = "mosaic_plot.jpeg", width = 748, height = 634, quality = 100, res = 100)
mosaicplot(TABLE_3, shade = TRUE, main = "", las = 2, xlab = "archiving organization", ylab = "identifiers related outputs") # mosaic plot
dev.off() 


## Detailed crosstabulations---------------------------------------------------------------------------------------------------------------------------------

MyTab2 <- xtabs(~ Identifiers_related_outputs + Client_ID + Affiliations_data_creators, data= data_repo_final)
MyTab3 <- xtabs(~ Identifiers_data_creators + Client_ID + Affiliations_data_creators, data = data_repo_final)
MyTab4 <- xtabs(~ Identifiers_data_creators + Client_ID + Affiliations_data_creators + Identifiers_related_outputs, data= data_repo_final)


## Detailed cotab plots--------------------------------------------------------------------------------------------------------------------------------------

cotabplot( ~ Identifiers_data_creators * Client_ID | Affiliations_data_creators, data = data_repo_final, split_vertical = FALSE, gp = shading_hcl, labeling_args = list(rot_labels = c(top = 80), offset_varnames = c(top = 3), offset_labels = c(top = 0.5)))

cotabplot( ~ Identifiers_related_outputs * Client_ID | Affiliations_data_creators, data = data_repo_final, split_vertical = FALSE, gp = shading_hcl, labeling_args = list(rot_labels = c(top = 80), offset_varnames = c(top = 3), offset_labels = c(top = 0.5)))


## Heatmap---------------------------------------------------------------------------------------------------------------------------------------------------

count_frame<-as.data.frame(MyTab4)
unite_DF <- count_frame %>% unite(New, Identifiers_data_creators, Affiliations_data_creators, Identifiers_related_outputs, sep = "__")
wide_DF <- unite_DF %>% spread(New, Freq)
wide_DF_2 <- data.frame(wide_DF[,-1], row.names=wide_DF[,1])
wide_DF_2 <- wide_DF_2 / 450 # sample size
wide_DF_3 <- wide_DF_2
colnames(wide_DF_3)<-c("yes_yes_no", "yes_yes_yes", "yes_no_no", "yes_no_yes", "no_yes_no", "no_yes_yes", "no_no_no", "no_no_yes")

heatmap_plot<-pheatmap(wide_DF_3, clustering_distance_rows = "euclidean", clustering_distance_cols = "canberra", clustering_method = "ward.D2", cutree_rows = 3, show_rownames = T, show_colnames = T, main = "", display_numbers = TRUE, angle_col = 45)

jpeg(file = "heatmap.jpeg", width = 748, height = 634, quality = 100, res = 100)
heatmap_plot
dev.off() 

wide_DF_3.dist<-dist(wide_DF_3, method = "euclidean")
wide_DF_3.hc <- hclust(wide_DF_3.dist, method = "ward.D2")
plot(wide_DF_3.hc, hang = -1)
rect.hclust(wide_DF_3.hc, k=3)
cluster.stats(wide_DF_3.dist, cutree(wide_DF_3.hc, k=3))
(asw<-sapply(2:8, function(x) summary(silhouette(cutree(wide_DF_3.hc, k=x), wide_DF_3.dist))$avg.width)) # Average Silhouette Width


## Plot visualizing percentages per repository----------------------------------------------------------------------------------------------------------------

unite_DF_2<-unite_DF
unite_DF_2$Freq <- unite_DF_2$Freq / 450
colnames(unite_DF_2)<-c("affiliation_info", "repository", "percentage")

(MS_plot <- ggplot(data = unite_DF_2, aes(x = affiliation_info, y = percentage, colour = repository, shape = repository)) + geom_point(size = 7) + theme_minimal() + scale_color_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#020202", "#a65628", "#f781bf", "#999999")) + scale_shape_manual(values = rep(15:17, len = 9)) + scale_y_continuous(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), limits = c(0.1,1)) + coord_flip()) # percentages below 0.1 are excluded


##-----------------------------------------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
