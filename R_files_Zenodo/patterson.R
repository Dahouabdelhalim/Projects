scale.mat<-function(dat,datID){
	#make vector of sqrt(p*(1-p))
	mv<-apply(dat,2,mean,na.rm=TRUE)
	vv<-mv/2
	vv<-sqrt(vv*(1-vv))
	#redo data and vectors to take out monomorphic sites
	dat<-dat[,vv>0]
	mv<-mv[vv>0]
	vv<-vv[vv>0]
	##M<-t((t(dat)-mv))
	M<-t((t(dat)-mv)/vv)
	M[is.na(M)]<-0
	rownames(M)<-as.character(datID[,1])
	return(M)
}

calc.TW<-function(lambda,m,idx){
	test<-read.table("twtable.txt")
	m<-m-idx+1
	lambda<-lambda[idx:length(lambda)]
	
	n<-(m+1)*sum(lambda)^2/((m-1)*sum(lambda^2)-sum(lambda)^2)  ##Equation 10

	sigma<-(sqrt(n-1)+sqrt(m))/n*(1/sqrt(n-1)+1/sqrt(m))^(1/3) ##Equation 6
	mu<-(sqrt(n-1)+sqrt(m))^2/n ##Equation 5

	L1<-(m-1)*lambda[1]/sum(lambda)		
	x<-(L1-mu)/sigma ##Equation 7
	
	dif<-abs(test[,1]-x)
	p<-test[dif==min(dif),2]
	
	return(c(x,p))
}

##Calculate significant PCA vectors using method from Patterson et al. PLoS Genetics (2006)
##clones<-NULL
##clones$id<-c(1:16,18:116)
##clones$pop<-trunc((clones$id-1)/10+1)
##clones<-as.data.frame(clones)

clones<-read.table("clones.txt",head=T)
datID <- data.frame(clones$pop)

dat <- read.table("defense_SNPs.txt",head=T)
M<-scale.mat(dat,datID)

m<-dim(M)[1] ##number of samples
n<-dim(M)[2] ##number of SNPs
X<-M%*%t(M)*1/n ##p. 2076

##Check for outlier samples
pc<-prcomp(M,scale=FALSE)
sdvec<-apply(pc$x,2,sd)
selmat<-abs(t(t(pc$x)/sdvec))

lambda<-svd(X)$d
##lambda<-eigen(X)$values

TW<-NULL
for(i in 1:10){
	TW<-rbind(TW,calc.TW(lambda,m,i))
}
TW

K<-min(which(TW[,2]>0.01))
##Plot cluster
d<-dist(svd(X)$u[,1:K])
x<-hclust(d,"ward")
plot(x)

x<-data.frame(cutree(x,(K)))
colnames(x)<-"classification"
un<-unique(as.vector(x$classification))
cols<-rainbow(length(un))
c<-match(x$classification,un)
cols<-cols[c]
##plot(svd(X)$u[,2]~clones$latitude,pch=19,col=cols)
plot(svd(X)$u[,1],svd(X)$u[,2],pch=19,col=cols)
