#Figure 1A
y <- readRDS("smallBRCADataset.rds")
mat <- y$data
labels <- y$pheno$SeuratCluster

#expression normalization
mat <- mat[,which(labels == "Basal")]
mat <- t(t(mat)/apply(mat, 2, sum))*1000000
mat <- log(mat+1,2)
mat <- data.frame(mat)

#HUGO subset
hugo <- read.table("/data1/ggulati/mberger/Cancer/SCE/HUGO.txt", header = T)
mat <- mat[as.character(unlist(hugo[,1])),]
mat <- mat[!is.na(rowSums(mat)),]
mat <- data.matrix(mat)

#overlap analysis
overlap <- sapply(1:nrow(mat), function(i) length(which(mat[i,]>0 & mat["VEGFA",]>0 & mat["THY1",]>0))/length(which(mat[i,]>0)))
names(overlap) <- rownames(mat)
counts <- rowSums(mat>0)
genes <- rev(sort(overlap[which(counts>5 & names(overlap) %in% hugo[,1] & names(overlap) %in% names(overlap)[grep("\\\\.",names(overlap),invert=T)])]))


write.table(genes, "Figure1_results.txt",sep = "\\t", quote = F)
