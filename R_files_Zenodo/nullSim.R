ests<-function(p0=NA,p1=NA){
	h<-0.5
	s<-(p1-p0)/(p0 * (1-p0) * (p0 + h * (1-2*p0)))
	s
}

wfsim<-function(N=NA,L=10000){
	p<-matrix(NA,nrow=L,ncol=2)
	pi<-runif(L,min=0.05,max=0.95)
	for(i in 1:L){
		g<-rbinom(n=N,size=2,prob=pi[i])
		survive<-sample(1:N,N*0.5,replace=FALSE)
		p[i,1]<-mean(g)/2
		p[i,2]<-mean(g[survive])/2
	}
	p
}

sim1<-wfsim(N=200)
sim2<-wfsim(N=500)
sim3<-wfsim(N=1000)
sim4<-wfsim(N=10000)

sims<-list(sim1,sim2,sim3,sim4)

sel<-matrix(NA,nrow=10000,ncol=4)
for(j in 1:4){
	sel[,j]<-ests(sims[[j]][,1],sims[[j]][,2])
}

library(RColorBrewer)

cs<-brewer.pal(n=9,"BuPu")[1+c(2,4,6,8)]

pdf("nulldist.pdf",width=5,height=11)
par(mfrow=c(2,1))
par(mar=c(4,5.5,3,0.5))
plot(density(sel[,4],from=0,to=1),main="",type='l',lwd=2,col=cs[4],xlim=c(0,1),ylab="density",xlab="selection coefficient",cex.lab=1.5,cex.axis=1.2)
for(j in 1:3){lines(density(sel[,j],from=0,to=1),lwd=2,col=cs[j])}
legend(0.69,11.7,c(200,500,1000,10000),lty=1,lwd=2,col=cs)
title(main="(a) s by sample size", adj=0,cex.main=1.2)

plot(sim4[,1],abs(sel[,4]),cex=0.8,col=cs[4],xlab="allele frequency",ylab="selection coefficient",cex.lab=1.5,cex.axis=1.2)
title(main="(b) effect of p on s",adj=0,cex.main=1.2)
dev.off()

