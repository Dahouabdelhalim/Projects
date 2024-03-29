load("Human_logitModel.RData")
test <- read.table(file="jurkat_chr22.ORF_info.tsv",sep="\\t",header=T)
test$Coding_prob <- predict(mylogit,newdata=test,type="response")
write.table(test, file="jurkat_chr22.ORF_prob.tsv", quote=F, sep="\\t",row.names=FALSE, col.names=TRUE)
