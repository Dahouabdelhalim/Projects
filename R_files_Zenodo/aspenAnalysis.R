library("gbs2ploidy")

## getting our data ready in R
ids<-read.csv("DT_aspen_1Jul16.csv") ## NV17 was dropped = 129

het<-as.matrix(read.table("hetAlleleDepth.txt",header=F))
a<-seq(1,380,2)
b<-seq(2,380,2)
cov1<-het[,a]
cov2<-het[,b]

## drop = 129 to match ids
cov1<-cov1[,-129]
cov2<-cov2[,-129]

## drop conatminant, FSL 279 = 92 ids and cov
ids<-ids[-92,]
cov1<-cov1[,-92]
cov2<-cov2[,-92]

propOut<-estprops(cov1=cov1,cov2=cov2,props=c(0.25, 0.33, 0.5, 0.66, 0.75),mcmc.nchain=3,mcmc.steps=100,mcmc.burnin=10,mcmc.thin=5)

H<-apply(is.na(cov1)==FALSE,2,mean)
D<-apply(cov1+cov2,2,mean,na.rm=TRUE)

out<-estploidy(alphas=propOut,het=H,depth=D,train=FALSE,nclasses=2,ids=ids[,4])

