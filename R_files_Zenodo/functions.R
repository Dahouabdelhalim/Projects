######################################################################
###
### Local functions and necessary global code
###
######################################################################

suppressMessages(library(edgeR))
suppressMessages(library(IRanges))
suppressMessages(library(gplots))
suppressMessages(library(cluster))
library(MASS)
#source("legacy_functions.R")

### ------------------------------------------------------------------
### BabelAssess
###
### Wrapper to run the primary analysis
### ------------------------------------------------------------------

BabelAssess=function(si,dn,pval,verbose=FALSE) {
	rownames(dn)=dn$gene
	# Combine the invidual-condition p-values as two-sided
	
	pval.2=apply(pval,2,one.to.2); # Two-sided p-values
	sig=list()
	for(tt in unique(si$cell_type)) {; # For each cell type
		for(ex in unique(si$condition)) {; # And for each condition
			ii=grep(paste(".*",ex,sep=""),colnames(pval)); # Fix this in the long-term, fine for now
			if(length(ii)==0) next
			sig[[paste("pval_comb",ex,sep="_")]]=comb.p(pval[,ii[1]],pval[,ii[2]])
		}
	}
	sig=as.data.frame(sig)
	rownames(sig)=rownames(pval)
	q=apply(sig,2,function(p) p.adjust(one.to.2(p),method="BH"))
	colnames(q)=gsub("pval","qval",colnames(q))

	#
	#
	# Go back to Z-statistic, and see betweenPvalues function. Due to a bug discovered by 
	# AO in the Hotelling as of 20120639
	#
	#

	for(ct in unique(si$cell_type)) {
		for(ex1 in unique(si$condition)) {
			for(ex2 in unique(si$condition)) {
				if(ex1<ex2) {
					p1=sig[[paste("pval_comb",ex1,sep="_")]]
					p2=sig[[paste("pval_comb",ex2,sep="_")]]

					# Until a proper p- and q-value is calculated between conditions, calculate
					# a z-statistic based metric and rank
					p11=pval[[paste("pval_exp1",ex1,sep="_")]]; z11=qnorm(p11)
					p12=pval[[paste("pval_exp2",ex1,sep="_")]]; z12=qnorm(p12)
					p21=pval[[paste("pval_exp1",ex2,sep="_")]]; z21=qnorm(p21)
					p22=pval[[paste("pval_exp2",ex2,sep="_")]]; z22=qnorm(p22)
					zstat=((z11+z12)-(z21+z22))/sqrt(2)
					sig[[paste("zstat_te",ex1,ex2,sep="_")]]=zstat
					sig[[paste("zstat_call",ex1,ex2,sep="_")]]=ifelse(abs(zstat)>=3.5,TRUE,FALSE)
					#zpvals=betweenPvalues(p1,p2)	
					#sig[[paste("pval_te",ex1,ex2,sep="_")]]=zpvals
					#sig[[paste("qval_te",ex1,ex2,sep="_")]]=p.adjust(zpvals,method="BH")
					res=exprAnalysis(dn[,grep(paste("mrna_(",ex1,"|",ex2,")$",sep=""),colnames(dn))]); # Expression analysis
					midx=match(rownames(sig),rownames(res))
					sig[[paste("qval_rna",ex1,ex2,sep="_")]]=res$FDR[midx]
					sig[[paste("logf_rna",ex1,ex2,sep="_")]]=res$logFC[midx]
				}
			}
		}
	}

	#
	## Between-condition testing as of 06.27.2012 for the two-condition / two-replicate case
	#for(ct in unique(si$cell_type)) {
	#	for(ex1 in unique(si$condition)) {
	#		for(ex2 in unique(si$condition)) {
	#			if(ex1<ex2) {
	#				p11=sig[[paste("pval_exp1",ex1,sep="_")]]
	#				p12=sig[[paste("pval_exp2",ex1,sep="_")]]
	#				p21=sig[[paste("pval_exp1",ex2,sep="_")]]
	#				p22=sig[[paste("pval_exp2",ex2,sep="_")]]			
	#				bh=betweenHotelling(p11,p12,p21,p22,robust=TRUE)			
	#				res=exprAnalysis(dn[,grep(paste("mrna_(",ex1,"|",ex2,")$",sep=""),colnames(dn))]); # Expression analysis
	#				midx=match(rownames(sig),rownames(res))
	#				sig[[paste("pval_te",ex1,ex2,sep="_")]]=bh$p.value
	#				sig[[paste("qval_te",ex1,ex2,sep="_")]]=p.adjust(bh$p.value,method="BH")
	#				sig[[paste("qval_rna",ex1,ex2,sep="_")]]=res$FDR[midx]
	#				sig[[paste("logf_rna",ex1,ex2,sep="_")]]=res$logFC[midx]
	#			}
	#		}
	#	}
	#}
	#

	sig$symbol=dn$symbol[match(rownames(sig),dn$gene)]
	sig=sig[,c(ncol(sig),1:(ncol(sig)-1)),]
	if(!verbose) {
		return(sig[,grep("(symbol|qval_comb|qval_te|qval_rna|logf_rna)",colnames(sig))])
	}
	return(sig)
}

### ------------------------------------------------------------------
### Between-condition significance testing for datasets that aren't
### Stumpf's HeLa cell cycle. Takes to matrices of p-values within
### conditions and returns one-sided p-values
### ------------------------------------------------------------------

betweenPvalues.2=function(p1mat,p2mat) {
	p1z=qnorm(p1mat)
	p2z=qnorm(p2mat)
	p1sum=apply(p1z,1,sum)
	p2sum=apply(p2z,1,sum)
	p1diff=apply(p1z,1,diff)
	p2diff=apply(p2z,1,diff)
	pval=pnorm((p1sum-p2sum)/sqrt(var(p1diff)+var(p2diff)))
	pval
}

### ------------------------------------------------------------------
### betweenHotelling
###
### Revision of the between-condition testing where p11 is a vector of
### p-values for the first condition and replicate, p12 is the same for
### first condition second replicate, p21 is the same for second
### condition first replicate, and p22 is the same for second condition
### second replicate.
### ------------------------------------------------------------------

betweenHotelling=function(p11,p12,p21,p22,robust=FALSE,filter=FALSE,alpha=0.001) {
	n <- length(p11)
	z11 <- qnorm(p11)
	z12 <- qnorm(p12)
	z21 <- qnorm(p21)
	z22 <- qnorm(p22)
	mean.vector <- t(cbind(z11-z21,z12-z22))
	sigma <- matrix(NA,2,2)
	if(!robust) {
		sigma[1,1] <- sum((z11-z21)^2)/n - (sum(z11-z21)/n)^2
		sigma[2,2] <- sum((z12-z22)^2)/n - (sum(z12-z22)/n)^2
	} else {
		sigma[1,1] <- mad(z11-z21)
		sigma[2,2] <- mad(z12-z22)
	}
	sigma[1,2] <- sum((z11-z21)*(z21-z22))/n
	sigma[2,1] <- sigma[1,2]
	eigen.sigma <- eigen(sigma)
	sigma.negative.onehalf <- t(eigen.sigma$vectors)%*%diag(eigen.sigma$values^(-0.5))%*%eigen.sigma$vectors
	stat.original <- sigma.negative.onehalf%*%mean.vector
	stat <- apply(stat.original^2,2,sum)
	p.value <- 1-pchisq(stat,2)
	p.value[which((mean.vector[1,]>0 & mean.vector[2,]<0)|(mean.vector[1,]<0 & mean.vector[2,]>0))] <- 1
	if(filter) {
		keepers <- which(p.value>alpha)
		if(length(keepers)>0) {
			n2 <- length(keepers)
			z112 <- qnorm(p11[keepers])
			z122 <- qnorm(p12[keepers])
			z212 <- qnorm(p21[keepers])
			z222 <- qnorm(p22[keepers])
			sigma2 <- matrix(NA,2,2)
			if(!robust) {
				sigma2[1,1] <- sum((z112-z212)^2)/n2 - (sum(z112-z212)/n2)^2
				sigma2[2,2] <- sum((z122-z222)^2)/n2 - (sum(z122-z222)/n2)^2
			} else {
				sigma2[1,1] <- mad(z112-z212)
				sigma2[2,2] <- mad(z212-z222)
			}
			sigma2[1,2] <- sum((z112-z212)*(z212-z222))/n2
			sigma2[2,1] <- sigma2[1,2]
			eigen.sigma2 <- eigen(sigma2)
			sigma.negative.onehalf2 <- t(eigen.sigma2$vectors)%*%diag(eigen.sigma2$values^(-0.5))%*%eigen.sigma2$vectors
			stat.original2 <- sigma.negative.onehalf2%*%mean.vector
			stat2 <- apply(stat.original2^2,2,sum)
			p.value <- 1-pchisq(stat2,2)
			p.value[which((mean.vector[1,]>0 & mean.vector[2,]<0)|(mean.vector[1,]<0 & mean.vector[2,]>0))] <- 1
		}
	}
	if(!filter) {
		list(p.value=p.value,z11=z11,z12=z12,z21=z21,z22=z22,mean.vector=mean.vector,sigma=sigma)
	} else {
		list(p.value=p.value,z11=z11,z12=z12,z21=z21,z22=z22,mean.vector=mean.vector,sigma=sigma,sigma2=sigma2)
	}
}

### ------------------------------------------------------------------
### ------------------------------------------------------------------

betweenPvalues=function(p1,p2) {
	stat=(qnorm(p1)-qnorm(p2))/sqrt(2)
	pval=2*(1-pnorm(abs(stat)))
	pval
}

### ------------------------------------------------------------------
### ------------------------------------------------------------------

clusterPvalues <- function(rna,rp,nreps,rnadisp, method="ls", trim.x=NULL, trim.y=FALSE, seed=NULL) {
	if(is.null(ncol(rna))) rna <- matrix(rna)
	if(is.null(ncol(rp))) rp <- matrix(rp)
	n <- ncol(rna)
	smat <- matrix(NA,nrow(rp),n)
	for(i in 1:n) {
		if(!is.null(seed)) set.seed(seed)
		rna.vector <- rna[,i]
		rp.vector <- rp[,i]
		#rpdisp <- rpDisp(rna.vector,rp.vector)$disp
		rpdisp <- dpDisp(rna.vector,rp.vector)$disp
		unique.rna.vector <- unique(rna.vector)
		unique.positions <- match(unique.rna.vector,rna.vector)
		match.indices <- match(rna.vector,unique.rna.vector)
		rna.vector.reduced <- rna.vector[unique.positions]
		rp.vector.reduced <- rp.vector[unique.positions]
		len <- nreps*length(unique.positions)
		nbsample <- rnbinom(len,mu=rna.vector.reduced,size=1/rnadisp)
		fit <- fitter(rna.vector,rp.vector,method=method,trim.x=trim.x,trim.y=trim.y)
		fitted.fit <- nbsample*fit$coefficients
		ref <- (matrix(rnbinom(len,mu=fitted.fit,size=1/rpdisp),ncol=nreps))[match.indices,]
		smat[,i] <- apply(ref>rp[,i],1,sum)
	}
	smat
}

### ------------------------------------------------------------------
### colormap.rb
###
### Accessory function for producing a legible color range for heatmaps
### ------------------------------------------------------------------

colormap.rb <- function(levels) {
  levels = levels%/%2
  sup = (1:(levels))/levels
  sdn = 1-(0:(levels-1))/levels
  return(c(rgb(1,sup,sup),rgb(sdn,sdn,1)))
}

### ------------------------------------------------------------------
### ------------------------------------------------------------------

comb.p=function(p1,p2) ifelse((p1+p2)<1,0.5*(p1+p2)^2,0.5+(p1+p2-1)*(1-((p1+p2-1)/2)))

### ------------------------------------------------------------------
### estLibSize
###
### Accessory function for calculating the estimated library size
### per edgeR (Robinson and Oshlack)
### ------------------------------------------------------------------

estLibSize <- function(counts,group) {
	d <- DGEList(counts=counts,group=group)
	d <- calcNormFactors(d)
	libsize <- d$samples$lib.size*d$samples$norm.factors
	libsize
}

### ------------------------------------------------------------------
### ------------------------------------------------------------------

exprAnalysis=function(de) {
	require(edgeR)
	dge=suppressMessages(DGEList(counts=de,group=gsub("^.*_","",colnames(de))))
	dsp=estimateCommonDisp(dge)
	res=topTags(exactTest(dsp),n=nrow(de))$table
	return(res)
}

### ------------------------------------------------------------------
### ------------------------------------------------------------------

colormap.rb <- function(levels) {
  levels = levels%/%2
  sup = (1:(levels))/levels
  sdn = 1-(0:(levels-1))/levels
  return(c(rgb(1,sup,sup),rgb(sdn,sdn,1)))
}


### ------------------------------------------------------------------
### fitter
###
### Regression of RNA on RP
### ------------------------------------------------------------------

fitter=function(rna.vector, rp.vector, method, trim.x=NULL, trim.y=FALSE, disp=0.2, alpha=0.001) {
	if(method=="rlm") {
		fit <- rlm(rna.vector,rp.vector)
	} else {
		if(!is.null(trim.x)) {
			quantile.low <- quantile(rna.vector,trim.x)
			quantile.high <- quantile(rna.vector,1-trim.x)
			keepers.x <- which(rna.vector>=quantile.low & rna.vector<=quantile.high)
			rna.vector <- rna.vector[keepers.x]
			rp.vector <- rp.vector[keepers.x]
		}
		if(trim.y) {
			fit1 <- lsfit(rna.vector,rp.vector,intercept=FALSE)
			fitted.values <- rna.vector*fit1$coefficients
			pvalue.vector <- pnbinom(rp.vector,mu=fitted.values,size=1/disp)
			keepers.y <- which(pvalue.vector>alpha & pvalue.vector<(1-alpha))
			rna.vector <- rna.vector[keepers.y]
			rp.vector <- rp.vector[keepers.y]
		}
		if(method=="ls") {
			fit <- lsfit(rna.vector,rp.vector,intercept=FALSE)
		} else if(method=="blue") {
			fit1 <- lsfit(rna.vector,rp.vector,intercept=FALSE)
			fit1.values <- rna.vector*fit1$coefficients
			fit1.phi <- rpDisp(rna.vector,rp.vector)$disp
			fit1.variances <- fit1.values*(1+fit1.values*fit1.phi)
			fit <- lsfit(rna.vector,rp.vector,wt=1/fit1.variances,intercept=FALSE)
		}
	}
	fit
}

### ------------------------------------------------------------------
### glog2
###
### Accessory function for glog2 conversion.
### ------------------------------------------------------------------

glog2=function(x,p0=0,p1=1) return((asinh(p0+p1*x)-log(2*p1))/log(2))

### ------------------------------------------------------------------
### logC
###
### Accessory function for logging on the command line.
### ------------------------------------------------------------------

logC=function(status,cmd,msg) {
  ts=gsub(" [A-Z]+.*$","",Sys.time())
  cat(status,"\\t",ts,"\\t",cmd,"\\t",msg,"\\n")
  if(status=="ERROR") q(save="no")
}


### ------------------------------------------------------------------
### one.to.2
###
### Convert one- to two-sided p-values. This is an update to a version
### AO wrote, but this one can do it on vectors.
### ------------------------------------------------------------------

one.to.2=function(p) 2*ifelse(p<(1-p),p,1-p); # Changed to support passing a vector of p-values

### ------------------------------------------------------------------
### plotGene
###
### Make summary plots of a given gene in the experiments
### ------------------------------------------------------------------

plotGene=function(gene,si,dn,cond1,cond2) {
	gidx=which(dn$symbol==gene)
	xlims=c(1e-7,1e-2)
	ylims=c(1e-7,1e-2)
	
	if(length(grep(paste(cond1,"|",cond2,sep=""),unique(si$condition)))<2) {
		stop("Passing unsupported conditions for this analysis")
	}
		
	dev.new()
	plot(1,1,type="n",xlim=xlims,ylim=ylims,log="xy",xlab="mRNA",ylab="Ribosome",main=gene,sub="(fraction of library)",axes=F)
	axis(1,at=10^(-7:-2),las=2,cex.axis=0.75)
	axis(2,at=10^(-7:-2),las=2,cex.axis=0.75)
	box()
	
	cols=NULL
	cols=cbind(cols,c("orange3","red3"))
	cols=cbind(cols,c("dodgerblue","darkblue"))
	colnames(cols)=c(cond1,cond2)
	rownames(cols)=paste("exp",unique(si$replicate),sep="")
	
	tt=NULL
	for(cc in c(cond1,cond2)) {
		for(ex in unique(si$replicate)) {			
			x=dn[,paste("exp",ex,"_mrna_",cc,sep="")]
			y=dn[,paste("exp",ex,"_rp_",cc,sep="")]				
			fit=fitter(x,y,method="ls",trim.x=0.1,trim.y=FALSE)
					
			xs=sum(x)
			ys=sum(y)
			x=x/xs
			y=y/ys
			crep=cols[ex,cc]
	
			#points(x,y,pch=20,cex=0.35,col="darkgray")
			abline(fit,untf=TRUE,col=crep)
			lines(range(x)/xs,(fit$coef*range(x))/ys,col="darkblue",lwd=1.5)
			points(x[gidx],y[gidx],pch=20,cex=2,col=crep)
			abline(v=x[gidx],lty=3,col=crep)
			abline(h=y[gidx],lty=3,col=crep)	
			tt=c(tt,paste(ex,"-",cc))
		}	
	}
	legend("bottomright",1e-7,tt,text.col=as.vector(cols),cex=0.5,bty="n")
}

### ------------------------------------------------------------------
### rpDisp
###
### Estimate the dispersion for RP
### ------------------------------------------------------------------

rpDisp=function(x,y,nbins=20,min.bin=10,trim.x=0.1,fixdisp=0.2,fixalpha=0.01) {
	order.x <- order(x)
	x <- x[order.x]
	y <- y[order.x]
	x.lower <- quantile(x,trim.x)
	x.upper <- quantile(x,1-trim.x)
	keepers <- which(x>=x.lower & x<=x.upper)
	x <- x[keepers]
	y <- y[keepers]
	fit <- lsfit(x,y,intercept=FALSE); # Now fit by ls
	fitted.values <- x*fit$coefficients
	j=1; while(j<3) {; # Do two runs, one w/ a guess of alpha, the other with the estimated alpha
		if(j==1) pvalues <- pnbinom(y,mu=fitted.values,size=1/fixdisp)
		else if (j==2) pvalues <- pnbinom(y,mu=fitted.values,size=1/disp)
		pvalues <- sapply(pvalues,one.to.two)
		keepers <- rep(TRUE,length(pvalues))
		keepers[which(pvalues<fixalpha)] <- FALSE
		x.keepers <- x[keepers]
		y.keepers <- y[keepers]
		x.bins <- quantile(x.keepers,seq(0,1,length.out=nbins+1)); #Start/end of bins based on quantiles so same number in every bin
		v <- rep(NA,nbins)
		lambda <- rep(NA,nbins)
		count.bins <- rep(NA,nbins)
		for(i in 1:nbins) {
			which.bin <- which(x.keepers>=x.bins[i] & x.keepers<x.bins[i+1])
			count.bins[i] <- length(which.bin)
			if(length(which.bin)>=min.bin) {
				current.y <- y.keepers[which.bin]
				v[i] <- var(current.y)
				x.middle <- median(x.keepers[which.bin])
				lambda[i] <- x.middle*fit$coefficients
			}
		}
		lambda <- lambda[!is.na(lambda)]
		v <- v[!is.na(v)]
		term1 <- sum(v*lambda*(1+lambda))
		term2 <- sum(lambda^3)
		term3 <- sum(lambda^4)
		disp <- (term1-term2)/term3
		j <- j+1
	}
	list(disp=disp,lambda=lambda,v=v,x.bins=x.bins,count.bins=count.bins,pvalues=sort(pvalues),fitted.values=fitted.values,keepers=keepers,x=x,y=y)
}

### ------------------------------------------------------------------
### rpm
###
### Calculate the proportion of ribosomes engaged in translation of 
### each transcript
### ------------------------------------------------------------------

rpm=function(reads,scale=1e6,log=FALSE) {; # RPM = proportion of ribosomes enganged in translation of that transcript
	r=reads/(sum(reads)/scale)
	if(log) {
		return(log(r))
	} else {
		return(r)
	}
}

### ------------------------------------------------------------------
### rpkm
###
### Calculate rpkM per transcript in an experiment
### ------------------------------------------------------------------

rpkm=function(reads,size,log=FALSE) {; # Good ole-fashioned rpkM per transcript
	r=(1e9)*((reads+1)/(sum(reads+1)*as.numeric(size)))
	if(log) {
		return(log(r))
	} else {
		return(r)
	}
	return(ifelse(log,log(r),r))
}

### ------------------------------------------------------------------
### SummarizeReadData
###
### Given the sample information file and read counts on genes in the
### RNA and RP experiments, produce the Ingolia-like TE values and
### and related summary information.
### ------------------------------------------------------------------

SummarizeReadData=function(si,dn) {
	# Make a table of Ingolia-like TE value
	te=NULL
	te$gene=dn$gene
	te$symbol=dn$symbol
	te$size=dn$size
	for(ex in unique(si$condition)) {	
		e1=paste("exp1","mrna",ex,sep="_")
		e2=paste("exp2","mrna",ex,sep="_")
		te[[paste("mrna",ex,sep="_")]]=rowMeans(cbind(rpkm(dn[,e1],dn$size),rpkm(dn[,e2],dn$size)))
		e1=paste("exp1","rp",ex,sep="_")
		e2=paste("exp2","rp",ex,sep="_")
		te[[paste("rp",ex,sep="_")]]=rowMeans(cbind(rpkm(dn[,e1],dn$size),rpkm(dn[,e2],dn$size)))
	}
	te=as.data.frame(te)
	for(ex in unique(si$condition)) {
		id1=paste("mrna",ex,sep="_")
		id2=paste("rp",ex,sep="_")
		id3=gsub("mrna","te",id1)
		te[[id3]]=log2(te[[id2]]/te[[id1]])
	}
	return(te)
}
