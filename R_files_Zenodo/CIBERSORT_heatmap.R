#Figure 1D
cx_results <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/ISPY/CIBERSORTx_output/CIBERSORTx_Adjusted.txt", row.names =1, header = T)
cx_results <- cx_results[,1:10]

ispy <- read.table("/data1/ggulati/scrnaseq/BreastGrant/U01/LMO2_paper/ISPY/ISPY_pheno.txt", row.names = 1, header = T)
er_status <- as.character(unlist((ispy[,1])))

reorg <- c(2, 3, 1, 4, 10, 5, 6, 7, 8, 9)
df_pos <- cor(cx_results[which(er_status == "P"),])
df_pos <- df_pos[reorg, reorg]
df_neg <- cor(cx_results[which(er_status == "N"),])
df_neg <- df_neg[reorg, reorg]

#heatmap
library(ComplexHeatmap)
library(circlize)
color_range = c(-0.3, -0.2, -0.1,0, 0.1,  0.2, 0.3)
color_palette = c("black", "#352A86", "#343DAE", "#1422e5", "#2Db7a3", "#F8BA43","#F8FA0D")
col_fun <- colorRamp2(color_range, color_palette)
myheatmap = Heatmap(df_pos,
		    border = T,
		    rect_gp = gpar(col = "black", lwd = 0.5),
                    name = "z-score",
                    col = col_fun,
                    row_title = "",
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    show_row_names = TRUE,
                    show_column_names = TRUE,
                    show_column_dend = FALSE,
                    show_row_dend = FALSE,
                    heatmap_legend_param = list(legend_direction = "horizontal",
                                                legend_width = unit(5, "cm"),
                                                title_position = "lefttop"))

pdf("CIBERSORT_LMO2_ERpos.pdf", width = 5, height = 5, useDingbats = FALSE)
draw(myheatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
