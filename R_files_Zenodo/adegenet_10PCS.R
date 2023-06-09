library(adegenet,quietly = T)

args<-commandArgs(trailingOnly = T)
if(length(args)==0){
  print("Input file missing")
  q()
}

fileName<-args[1]
data_all<-read.snp(fileName)
pca_all<- glPca(data_all,nf = 10)
pca_scores_all<-data.frame(pca_all$scores)
write.table(pca_scores_all, file = "cov.csv", sep = '\\t')
