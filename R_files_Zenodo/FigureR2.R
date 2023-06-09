#Heatmaps of luminal and basal markers
###NAVIN####
mat <- read.table("/data1/ggulati/mango/breastdata/Navin_chemoresistance/clonally_extinct_exprs.txt", row.names =1, header = T)
pheno <- read.table("/data1/ggulati/mango/breastdata/Navin_chemoresistance/clonally_extinct_pheno.txt", row.names =1, header = T)

#HUGO subset
mat <- data.frame(mat)
hugo <- read.table("/data1/ggulati/mberger/Cancer/SCE/HUGO.txt", header = T)
mat <- mat[as.character(unlist(hugo[,1])),]
mat <- mat[!is.na(rowSums(mat)),]
mat <- data.matrix(mat)

#expression normalization
mat <- t(t(mat)/apply(mat, 2, sum))*1000000
mat <- log(mat+1,2)

#Subset
genes <- c("KRT14", "KRT5", "KRT17", "KRT8", "KRT18", "KRT19")
type <- rep(c("Basal", "Luminal progenitor"), each = 3)
lmo2 <- mat["LMO2",which(mat["LMO2",]>0)]
mat <- mat[genes,which(mat["LMO2",]>0)]
basal <- as.numeric(apply(mat[1:3,], 2, mean)>0)
lp <- as.numeric(apply(mat[4:6,], 2, mean)>0)

ord <- apply(mat[1:3,], 2, mean) - apply(mat[4:6,], 2, mean)
mat <- mat[,order(ord)]


#heatmap
library(ComplexHeatmap)
library(circlize)
makeHaobject = function(legend, colors){
  annotlist = list()
  for (i in 1:length(legend)){
    classlist = legend
    colorlist = colors
    classi = classlist[i]
    colori = colorlist[i]
    annotlist[[i]] =  colori
    names(annotlist[[i]]) = classi
  }
  return(unlist(annotlist))
}

scale_apply <- function(x) {
        apply(x, 1, function(y) (y - mean(y))/sd(y))
    }

df <- t(scale_apply(mat))
row_split <- type
##create col bars and adjust color range and row bars??
color_range = c(-2, -1, -0.5,0, 0.5,  1, 2)
color_palette = c("black", "#352A86", "#343DAE", "#1422e5", "#2Db7a3", "#F8BA43","#F8FA0D")
col_fun <- colorRamp2(color_range, color_palette)

myheatmap = Heatmap(df,
                    name = "z-score",
                    col = col_fun,
                    row_title = "Genes",
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    show_row_names = TRUE,
                    show_column_names = FALSE,
                    show_column_dend = FALSE,
                    show_row_dend = FALSE,
                    row_split = row_split,
                    heatmap_legend_param = list(legend_direction = "horizontal",
                                                legend_width = unit(5, "cm"),
                                                title_position = "lefttop"))

