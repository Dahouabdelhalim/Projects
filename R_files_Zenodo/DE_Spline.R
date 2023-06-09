RawExpData<-as.matrix(read.table("Raw_Counts.txt",sep="\\t",header=T,row.names=1))
Nit<-c("noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit")
Time<-c(0,0,0,5,5,5,10,10,10,15,15,15,20,20,20,30,30,30,45,45,45,60,60,60,90,90,90,120,120,120,5,5,5,10,10,10,15,15,15,20,20,20,30,30,30,45,45,45,60,60,60,90,90,90,120,120,120)

library(splines)
library(limma)
X<-ns(Time,df=2)
design<-model.matrix(~X*Nit)


TotalCounts<-rowSums(RawExpData)
AveCounts<-(TotalCounts/60)
hist(TotalCounts,breaks=1000000,xlim=c(0,1000))
hist(AveCounts,breaks=500000,xlim=c(0,100))
FiltExp<- RawExpData[AveCounts>10,]

library(EDASeq)
fullmat <- betweenLaneNormalization(FiltExp, which="full",offset=FALSE)

write.table(fullmat,file="Normalized.FQ.txt",sep="\\t",quote=F)
fit<-lmFit(fullmat,design)
fit<-eBayes(fit)
Qle0.05<-topTable(fit,coef=5:6,adjust="BH",p.value=0.05,number=20000)

Quant.df2.0.05.ids<-rownames(Qle0.05)

write.table(format(Qle0.05,digits=3),file="DE_Genes.0.05.txt",sep="\\t",quote=F)
