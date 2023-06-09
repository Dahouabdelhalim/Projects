source("http://bioconductor.org/biocLite.R")#Only run once to download
source("https://bioconductor.org/biocLite.R")
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade") 
biocLite("DESeq2")#Only run once to download 
library(DESeq2)
library(xlsx)
require("DESeq2")

read = read.csv2("read_count_egg_rem_filt4.csv",row.names = 1,sep=",",header=T)
head(read)
dds<-DESeq(read)
count<-counts(dds,normalize=TRUE)
write.csv(count, "count_Eggrem_filt4.csv")

#Egg removal
sample_table <- read.table("SampleInf_egg_rem.txt", sep=",", header=T)

#Create dataframe
Treat <- sample_table[,2]
Col <- factor(sample_table[,1]) #Important to convert this into factor otherwise Deseq2 thinks it's numerical
all <- list(Col,Treat)
Design = as.data.frame(all)
colnames(Design) <- c("Col","Treat")
row.names(Design)<-colnames(read)
dds <- DESeqDataSetFromMatrix(countData= read,colData = Design, design = ~ Col+Treat)
dds <- DESeq(dds)
plotPCA(dds,intgroup = c("Treat"))

results <- results(dds,contrast = c('Treat',"Egg","Ref"))
RSEM_filtered_results_sig_B <- subset(results, padj < 0.05)
install.packages("xlsx")
library(xlsx)
write.csv(RSEM_filtered_results_sig_B, "DGE_Eggrem_filt4.csv")

# Food restriction

read = read.csv2("read_count_food_rest_filt4.csv",row.names = 1,sep=",",header=T)
head(read)

Treat = c("Exp","Ctrl","Exp","Ctrl","Exp","Ctrl","Exp","Ctrl","Exp","Ctrl","Exp","Ctrl","Exp","Ctrl","Exp","Ctrl")
Design = as.data.frame(cbind(Treat))
row.names(Design)<-colnames(read)
dds <- DESeqDataSetFromMatrix(countData= read,colData = Design, design = ~ Treat)
dds <- DESeq(dds)
plotPCA(dds,intgroup = c("Treat"))

results <- results(dds,contrast = c('Treat',"Exp","Ctrl"))
RSEM_filtered_results_sig_B <- subset(results, padj < 0.05)
library(xlsx)
write.csv(RSEM_filtered_results_sig_B, "DGE_Food_Rest_fil4.csv")
