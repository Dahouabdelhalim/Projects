library(GenABEL,quietly = T)
args<-commandArgs(trailingOnly = T)
if(length(args)==0){
  print("Input file missing")
  q()
}
fileName<-args[1]
phenotype<-read.csv(file=fileName,header = TRUE, sep = ,)
phenotype$phenotype_residuals_RT<-rntransform(phenotype$phenotype_residuals)
write.csv(phenotype, file="phenotype_RT.csv")
