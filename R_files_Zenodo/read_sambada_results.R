library(qvalue)
library(stringr)
library(data.table)

SAM_uni <- read.table('./mol-data-Phenig2-Out-1.txt', header=T)

SAM_uni$pvalueG <- 1-pchisq(SAM_uni$Gscore, df=1)
SAM_uni$qvalueG <- qvalue(SAM_uni$pvalueG)$qvalue

SAM_uni$scaffold <- as.numeric(str_split_fixed(SAM_uni$Marker, "_", 3)[,1])
SAM_uni$position <- as.numeric(str_split_fixed(SAM_uni$Marker, "_", 3)[,2])

sigUNI_all_g <- SAM_uni[SAM_uni$qvalueG<0.05,] 

sigUNI_all_g$position<-as.numeric(sigUNI_all_g$position)
sigUNI_all_g$scaffold<-as.numeric(sigUNI_all_g$scaffold)

