## Tests accounting for phylogenetic uncertainty

ff<-function(tree,x,y){
	pic.x<-pic.ortho(x,tree,intra=T)
	pic.y<-pic.ortho(y,tree,intra=T)
	fit<-lm(pic.y~pic.x-1)
	sum<-as.numeric(summary.lm(fit)$adj.r.squared)
	setNames(c(coef(fit),vcov(fit),sum),c("beta","var(beta)","R2"))
}
# now apply to all trees in your sample

BB<-t(sapply(trees,ff,x,y))


# total variance in beta estimated by
varBeta<-var(BB[,"beta"])+mean(BB[,"var(beta)"])
t.beta<-mean(BB[,"beta"])/sqrt(varBeta)
P.beta<-2*pt(abs(t.beta),df=length(trees[[1]]$tip)-2,lower.tail=FALSE)
