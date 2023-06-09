#load the adegenet library 
library(adegenet,quietly = T)

#make sure an input is given
args<-commandArgs(trailingOnly = T)
if(length(args)==0){
  print("Input file missing")
  q()
}
fileName<-args[1]

#read in a SNP data file (i.e. with a '.snp' extension)
#the read.snp function converts the input into a genlight object
data<-read.snp(fileName)

#implement PCA for the genlight object, retaining the first 2 principal components
#results will be stored as a glPCA list
pca<- glPca(data,nf = 2)

#make new data frame for eigenvalues (which correspond to absolute variance) 
eig=data.frame(pca$eig)

#calculate percent variance explained by each principal component and save the data frame
eig$percent = (eig[, 1]/sum(eig$pca.eig))
outfile_name<-paste (fileName,".eig.csv",sep="")
write.table(eig, file = outfile_name, sep = '\\t', row.names = FALSE)

#make and save new data frame containing the coordinates of each individual on each PC axis
pca_scores<-data.frame(pca$scores)
pca_scores['samples'] <- NA
pca_scores$samples<-data$ind.names
outfile_name<-paste (fileName,".results.csv",sep="")
write.table(pca_scores, file = outfile_name, sep = '\\t', row.names = FALSE)