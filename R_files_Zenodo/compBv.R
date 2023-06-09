## compare breeding values across genetic archs.
L<-17098
N<-c(329,344,347) # L1, L2, L14
gf<-list.files(pattern="tr")[c(2:3,1)]

G<-vector("list",3)
for(k in 1:3){
	G[[k]]<-matrix(scan(gf[[k]],n=L*N[k],sep=" "),nrow=L,ncol=N[k],byrow=TRUE)
}

## read genetic arch.
bbar<-matrix(scan("cmac_effects.txt",n=L*6,sep=" "),nrow=L,ncol=6,byrow=TRUE)

## weight
wmat<-vector("list",3)
## outer is population, inner is gwas
for(k in 1:3){
	wmat[[k]]<-matrix(NA,nrow=N[k],ncol=3)
	for(gwa in 1:3){
       		for(j in 1:N[k]){
			wmat[[k]][j,gwa]<-sum(G[[k]][,j] * bbar[,gwa],na.rm=TRUE)
		}
	}
}

## dev. time
dmat<-vector("list",3)
## outer is population, inner is gwas
for(k in 1:3){
	dmat[[k]]<-matrix(NA,nrow=N[k],ncol=3)
	for(gwa in 1:3){
       		for(j in 1:N[k]){
			dmat[[k]][j,gwa]<-sum(G[[k]][,j] * bbar[,gwa+3],na.rm=TRUE)
		}
	}
}

## summarize correlations, CIs and p-values
wcors<-matrix(NA,nrow=9,ncol=4)
dcors<-matrix(NA,nrow=9,ncol=4)
a<-c(1,1,2)
b<-c(2,3,3)

for(k in 1:3){for(cc in 1:3){
	i<-cc + (k-1) * 3
	o<-cor.test(wmat[[k]][,a[cc]],wmat[[k]][,b[cc]],na.rm=TRUE)
	wcors[i,]<-c(o$estimate,o$conf.int[1:2],o$p.value)
	o<-cor.test(dmat[[k]][,a[cc]],dmat[[k]][,b[cc]],na.rm=TRUE)
	dcors[i,]<-c(o$estimate,o$conf.int[1:2],o$p.value)
}}

## cors mixture of 0, pos, neg
## cors of bbar much higher for abs() then raw
## suggests that different LD, different SNPs tagging

save(list=ls(),file="bv.rdat")


