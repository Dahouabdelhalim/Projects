library(LEA)
###

###
for (k in c(1:200)){
vcf2lfmm(paste("Qa",k,"lfmm","vcf", sep="."))}


for (i in c(1:200)){

project2 = NULL
project2 = lfmm(paste("Qa",i,"lfmm","lfmm", sep="."),"pc1.load.env",K = 2,repetitions = 5,project = "new")
###
p = lfmm.pvalues(project2, K = 2)
pvalues = p$pvalues
print (pvalues)   ### total snp ,each snp correspond a pvalue
pdf(paste("K2","hist",i,"pdf", sep="."))
hist(pvalues, col = "lightblue")
dev.off()
write.table(x=pvalues,file=paste("K2pc1",i,"pvalues",sep="."),quote=F, sep="\\t", row.names=F)
###
library(qvalue)
for (alpha in c(0.01, 0.05)) {
# expected FDR
print(paste("Expected FDR:", alpha))
qval <- qvalue(pvalues)$qvalue
outliers <- which(qval < alpha)
write.table(x=outliers,file=paste("K2pc1",alpha,"outlier",i,"txt",sep="."),quote=F, sep="\\t", row.names=F)
}
}