setwd("/Users/markh/Desktop/HFSP/")
#install.packages("BiocManager")
#BiocManager::install("devtools")
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library("ComplexHeatmap")
library("readxl")
library(circlize)

queen_genes <- read_excel("gene-list-aging.xlsx",  ##T0-T4, alternatives
                          sheet = "sheet2", col_names = TRUE)

queen_genes <- unique(queen_genes)
queen_genes$order <- 1:length(queen_genes$gene_id)
expression <- read.table("counts_new.stand.caste.DEG",header=T)
wgcna <- read.table("geneInfo_mnat_all.tab",header = T)

queen_data <- merge(queen_genes[,c(1:3,6,7)], expression,by = "gene_id")
queen_data <- merge(queen_data, wgcna[,1:2], by = "gene_id",all.x=T)
queen_data <- queen_data[,c(1:3,13,14,4:12,15:43)]
queen_data[queen_data$gene_id == "Mnat_00258","Dmel_Symbol"] <- "Ilp9"


queen_data <- queen_data[order(queen_data$order),]


################
ageing_data <- queen_data
subpathways <- factor(ageing_data$sub_pathway, levels = unique(ageing_data$sub_pathway))


row_ann <- rowAnnotation(
  WGCNA = ageing_data$moduleColor, 
  col = list(
    WGCNA = c("blue" = "blue", "darkorange2" = "orange", "lightcyan" = "lightcyan", "lightgreen" = "lightgreen", "lightyellow" = "lightyellow", "pink" = "pink", "red" = "red",  "yellow" = "yellow", "plum1" = "plum")), border=T,
  Gene = anno_text(gsub("[\\r\\n]", "",ageing_data$gene_name), gp = gpar(fontsize = 6)), #gsub("[\\r\\n]", "", x)
  Hsap = anno_text(ageing_data$Hsap_Symbol, gp = gpar(fontsize = 8)),
  Dmel = anno_text(ageing_data$Dmel_Symbol, gp = gpar(fontsize = 8))
)

pdf("Heatmap_ageing_expression.pdf",paper='A4', width = 21/2.54, height = 29.7/2.54)
par(mar=c(3,3,3,3))
hm <- Heatmap(as.matrix(ageing_data[,c(8,13:14)]),cluster_columns = F, cluster_rows = T,# rect_gp = gpar(col = "grey", lwd = 0.5),
              column_title = "ageing expression", column_title_gp = gpar(fontsize = 14, fontface = "bold"),
              row_split = subpathways, row_title_gp = gpar(fontsize = 6), cluster_row_slices = F,
              show_row_names = F,
              show_row_dend = F,right_annotation = row_ann,
              border = T,
              heatmap_legend_param = list(title = "expression"),
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              #cell_fun = function(j, i, x, y, width, height, fill) {
              #grid.text(deg_ann[i, j], x, y, gp = gpar(fontsize = 10))
              cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.1f", ageing_data[,c(8,13:14)][i, j]), x, y, gp = gpar(fontsize = 8))
              }
)
print(hm)
dev.off()
