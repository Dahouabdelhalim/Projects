#adjust P values using the lambda genomic inflation factor
library(GenABEL)
args<-commandArgs(trailingOnly = T)

if(length(args)==0){
  
  print("Input file missing")
  q()
}

fileName<-args[1]
df<-read.table(fileName,as.is=TRUE,header=TRUE)
pvals<-df$P
genomic_inflation<-estlambda(pvals, plot=FALSE)
lambda<-genomic_inflation$estimate
chisq <- qchisq(pvals,1,lower.tail=FALSE)
newchisq <- chisq/lambda
df$Pc <- pchisq(newchisq, df=1,lower.tail=FALSE)
write.table(df,paste0(fileName,".Padjusted"),col=T,row=F,qu=F)