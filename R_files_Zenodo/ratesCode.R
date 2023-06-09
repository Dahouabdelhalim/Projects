Pf<-list.files(pattern="p_")
N<-length(Pf)
P<-vector("list",N)
for(i in 1:N){
    P[[i]]<-as.matrix(read.table(Pf[i],sep=",",header=FALSE))
    }
    
snps<-read.table("mcmc_M_L1F100.txt",header=FALSE,sep=",")  
L<-dim(snps)[1]
    
X<-read.table("X.contigs.txt",header=FALSE)
LGX<-(snps[,1] %in% X[,1])    
Y<-read.table("Y.contigs.txt",header=FALSE)
LGY<-(snps[,1] %in% Y[,1])     
LGA<-rep(1,L)
LGA[LGX==1 | LGY==1]<-0    
    
calcRate<-function(a=NA,b=NA){## a is ancestral, b is derived
    dp<-b-a
    v<-2*a*(1-a)
    sdp<-sqrt(v)
    dprat<-abs(dp)/sdp
    return(dprat)
}    
# mung = 22 vs. L1 F100 = 1

Comps<-5
erates<-matrix(NA,nrow=Comps,ncol=L)
# mung = 22 vs. L1 F100 = 1
erates[1,]<-calcRate(a=P[[22]][,1],b=P[[1]][,1])
# mung = 22 vs. L2 F87 = 19
erates[2,]<-calcRate(a=P[[22]][,1],b=P[[19]][,1])
## L1 F1-91 = 2 vs L1R F46 = 17
erates[3,]<-calcRate(a=P[[2]][,1],b=P[[17]][,1])
## L2 F78 = 18 vs L2R F45 = 21
erates[4,]<-calcRate(a=P[[18]][,1],b=P[[21]][,1])
## L 14 P = 5, L14A F16 = 6
erates[5,]<-calcRate(a=P[[5]][,1],b=P[[6]][,1])

## ordered, auto, X, Y
oerates<-erates[,c(which(LGA==1),which(LGX==1),which(LGY==1))]
lgs<-rep(c("A","X","Y"),c(sum(LGA),sum(LGX),sum(LGY)))

## fit HMM to rates of change
library(HiddenMarkov)

## Gaussian HMM with mean and sd fixed
Mstep.normc <- function(x, cond, pm, pn){
        nms <- sort(names(pm))
         n <- length(x)
         m <- ncol(cond$u)
         if (all(nms==c("mean", "sd"))){
             mean <- pm$mean
             sd <- pm$sd
             
             return(list(mean=mean, sd=sd))
         }
         if (all(nms == "mean")) {
             mean <- pm$mean
             return(list(mean = mean))
         }
         if (all(nms == "sd")) {
             sd <- pm$sd

             return(list(sd = sd))
         }
         stop("Invalid specification of parameters")

}

rnormc <- rnorm
dnormc <- dnorm
pnormc <- pnorm
qnormc <- qnorm

delta<-c(0.95,0.05)
PiA<-matrix(c(0.9,0.1,0.1,0.9),nrow=2,byrow=T)
PiB<-matrix(c(0.99,0.01,0.2,0.8),nrow=2,byrow=T)
Pi<-list(PiA,PiB)  

## variables for storing results
fitA<-vector("list",Comps)
fitB<-vector("list",Comps)
paramA<-vector("list",Comps)
paramB<-vector("list",Comps)

for(i in 1:Comps){
    mns<-quantile(erates[i,],probs=c(0.5,0.999))## median and 95th
    sds<-rep(sd(erates[i,]),2)
    
    ## rep A
    init<-dthmm(oerates[i,],PiA,delta,"normc",list(mean=mns,sd=sds),discrete=FALSE)
    paramA[[i]]<-BaumWelch(init,bwcontrol(maxiter = 500, tol = 1e-04,prt = TRUE, posdiff = FALSE,converge = expression(diff < tol)))
    fitA[[i]]<-Viterbi(paramA[[i]])

    ## rep B
    init<-dthmm(oerates[i,],PiB,delta,"normc",list(mean=mns,sd=sds),discrete=FALSE)
    paramB[[i]]<-BaumWelch(init,bwcontrol(maxiter = 500, tol = 1e-04,prt = TRUE, posdiff = FALSE,converge = expression(diff < tol)))
    fitB[[i]]<-Viterbi(paramB[[i]])
}

## plots
library(scales)
cs<-alpha(c("darkgray","red"),.65)
lb<-min(which(lgs=="X"));ub<-max(which(lgs=="X"))
tits<-c("(A) L1","(B) L2","(C) L1R","(D) L2R","(E) L14A") 
pdf("cmac_rates.pdf",width=7,height=11)
par(mfrow=c(5,1))
par(mar=c(4.5,4.5,2,1))
for(i in 1:Comps){
    plot(oerates[i,],col=cs[fitA[[i]]],pch=20,type='n',xlab="SNP",ylab="Change",cex.lab=1.5)
    yub<-max(oerates[i,],cex.axis=1.1)
    polygon(c(lb,lb,ub,ub),c(0,yub,yub,0),col=alpha("goldenrod",.4),border=NA)
    points(oerates[i,],col=cs[fitA[[i]]],pch=20)
    title(main=tits[i],cex.main=1.4)
    }
dev.off()

library(corrgram)
pdf("rats_corrgram.pdf",width=5,height=5)
corrgram(t(erates),upper.panel=panel.cor,labels=c("L1","L2","L1R","L2R","L14A"))
dev.off()

hmmstat<-as.matrix(cbind(fitA[[1]],fitA[[2]],fitA[[3]],fitA[[4]],fitA[[5]]))
pdf("hmm_corrgram.pdf",width=5,height=5)
corrgram(hmmstat,upper.panel=panel.cor,labels=c("L1","L2","L1R","L2R","L14A"))
dev.off()
## statistical summaries

## number and prop SNPs in high state (note two fits identical)
histat<-matrix(NA,nrow=Comps,ncol=2)
for(i in 1:Comps){
    histat[i,]<-c(sum(fitA[[i]]==2),mean(fitA[[i]]==2))
}
 histat
     [,1]        [,2]
#[1,]  113 0.005545740
#[2,]  109 0.005349431
#[3,]  174 0.008539458
#[4,]  150 0.007361602
#[5,]  122 0.005987436

hiwins<-histat
for(i in 1:Comps){
	o<-cbind(fitA[[i]][-1],fitA[[i]][-20376])
	hiwins[i,1]<-sum(o[,2]==1 & o[,1]==2)
	hiwins[i,2]<-histat[i,1]/hiwins[i,1]
}
#     [,1]     [,2]
#[1,]   87 1.298851
#[2,]   83 1.313253
#[3,]  118 1.474576
#[4,]   97 1.546392
#[5,]   95 1.284211

## transition probs
round(paramA[[1]]$Pi,3)
#      [,1]  [,2]
#[1,] 0.996 0.004
#[2,] 0.764 0.236
round(paramA[[2]]$Pi,3)
#      [,1]  [,2]
#[1,] 0.996 0.004
#[2,] 0.759 0.241
round(paramA[[3]]$Pi,3)
#      [,1]  [,2]
#[1,] 0.994 0.006
#[2,] 0.680 0.320
round(paramA[[4]]$Pi,3)
#      [,1]  [,2]
#[1,] 0.995 0.005
#[2,] 0.654 0.346
round(paramA[[5]]$Pi,3)
#      [,1]  [,2]
#[1,] 0.995 0.005

## cors in state and for snps in state
rcs<-matrix(NA,nrow=5,ncol=5)
pcs<-rcs
for(i in 1:5){for(j in 1:5){
	if(i > j){
		o<-cor.test(fitA[[i]],fitA[[j]])
		rcs[i,j]<-o$estimate
		pcs[i,j]<-o$p.value
	}
	if(j > i){
		kp<-which(fitA[[i]]==2 | fitA[[j]]==2)
		o<-cor.test(oerates[i,kp],oerates[j,kp])
		rcs[i,j]<-o$estimate
		pcs[i,j]<-o$p.value
	}
}}


## prop and enrichment on X
xenrich<-matrix(NA,nrow=Comps,ncol=3)
Nx<-sum(lgs=='X')
for(i in 1:Comps){
    null<-rep(NA,1000)
    obs<-sum(fitA[[i]][lgs=='X']==2)
    for(j in 1:1000){
        ix<-sample(1:L,Nx,replace=FALSE)
        null[j]<-sum(fitA[[i]][ix]==2)
    }
    xenrich[i,]<-c(obs,round(mean(null),1),mean(null>obs))
}   

xenrich
#      obs null    P
#[1,]   12  9.1 0.123
#[2,]   16  8.9 0.004
#[3,]   10 13.9 0.828
#[4,]    4 12.2 0.995
#[5,]   20  9.7 0.001

save(list=ls(),file="hmm.rdat")
     
## comparison with QTL results
oord<-c(which(LGA==1),which(LGX==1),which(LGY==1))
map<-read.table("mapping_loci.txt",sep=",",header=FALSE)
mapi<-rep(NA,dim(map)[1])
omapi<-mapi
for(i in 1:dim(map)[1]){
    mm<-which(snps[,1]==map[i,1] & snps[,2]==map[i,2])  
    if(length(mm)==1){
        mapi[i]<-mm
        omapi[i]<-which(oord==mm)
    }
}       
 

## L1 w, L2 w, L14 w, L1 dt, L2 dt, L14 dt
pips<-read.table("cmac_pips.txt",header=FALSE)
mav<-read.table("cmac_effects.txt",header=FALSE)

## comparison of QTL density
## L1 w, dL1
tapply(X=pips[,1],INDEX=fitA[[1]][omapi],mean,na.rm=TRUE)
#          1           2 
#0.004378731 0.005430733 
## L1 w, dL1R
tapply(X=pips[,1],INDEX=fitA[[3]][omapi],mean,na.rm=TRUE)
#          1           2 
#0.004385310 0.004316925 
## L1 dt, dL1
tapply(X=pips[,4],INDEX=fitA[[1]][omapi],mean,na.rm=TRUE)
#          1           2 
#0.003458553 0.003747900 
## L1 dt, dL1R
tapply(X=pips[,4],INDEX=fitA[[3]][omapi],mean,na.rm=TRUE)
#          1           2 
#0.003461266 0.003339871 

## L2 w, dL2
tapply(X=pips[,2],INDEX=fitA[[2]][omapi],mean,na.rm=TRUE)
#          1           2 
# 0.001859266 0.001808577 
## L2 w, dL2R
tapply(X=pips[,2],INDEX=fitA[[4]][omapi],mean,na.rm=TRUE)
#          1           2 
#0.001859349 0.001813688  
## L2 dt, dL2
tapply(X=pips[,5],INDEX=fitA[[2]][omapi],mean,na.rm=TRUE)
#          1           2 
#0.003933354 0.003771231  
## L2 dt, dL2R
tapply(X=pips[,5],INDEX=fitA[[4]][omapi],mean,na.rm=TRUE)
#          1           2 
#0.003931919 0.004019013

## L14 w, dL14
tapply(X=pips[,3],INDEX=fitA[[5]][omapi],mean,na.rm=TRUE)
#          1           2 
#0.003647183 0.003799913  
## L14 dt, dL14
tapply(X=pips[,6],INDEX=fitA[[5]][omapi],mean,na.rm=TRUE)
#          1           2 
#0.003543516 0.003560145

## put it all in a loop and run a null model
Xmap<-c(1,1,4,4,2,2,5,5,3,6)
Xdp<-c(1,3,1,3,2,4,2,4,5,5)
Xl<-10

obs<-rep(NA,Xl)
null<-matrix(NA,nrow=Xl,ncol=1000)
pv<-rep(NA,Xl)
for(i in 1:Xl){
    a<-pips[,Xmap[i]]
    b<-fitA[[Xdp[i]]][omapi]
    obs[i]<-mean(a[b==2],na.rm=TRUE)
    for(j in 1:1000){
        xshift<-sample(2:length(a),1)
        xnull<-c(xshift:length(a),1:(xshift-1))
        null[i,j]<-mean(a[xnull][b==2],na.rm=TRUE)
    }
    pv[i]<-mean(null[i,] >= obs[i])    
}

qs<-apply(null,1,quantile,probs=c(.5,.025,.975))

pdf("qtlXdp.pdf",width=5,height=5)
par(mar=c(5,5,.5,.5))
po<-c(1:2,5:6,9,3:4,7:8,10)
plot(qs[1,po],ylim=c(0,0.006),pch=20,axes=FALSE,xlab="Line",ylab="QTL density",cex.lab=1.4)
axis(2)
axis(1,at=1:10,rep(c("L1","L1R","L2","L2R","L14"),2),las=2)
box()
mtext("weight",3,adj=.2,line=-1.75,cex=1.3)
mtext("dev. time",3,adj=.8,line=-1.75,cex=1.3)
segments(1:10,qs[2,po],1:10,qs[3,po])
points(1:10,obs[po],col="red",pch=19)
abline(v=5.5,lty=2)
dev.off()

## in original order, not plot order
#       obs          pv
# [1,] 0.00543 0.002
# [2,] 0.00432 0.631
# [3,] 0.00375 0.006
# [4,] 0.00334 0.901
# [5,] 0.00181 0.600
# [6,] 0.00181 0.604
# [7,] 0.00377 0.752
# [8,] 0.00402 0.329
# [9,] 0.00380 0.107
#[10,] 0.00356 0.468



