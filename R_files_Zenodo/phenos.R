##### Read in phenotypes for each mapping line. First column is ID, second is weight (mg), third is development time (days).

L14<-read.table('L14_phenotypes.txt',header=FALSE,sep=' ')
L1<-read.table('L1_phenotypes.txt',header=FALSE,sep=' ')
L2<-read.table('L2_phenotypes.txt',header=FALSE,sep=' ')

t14<-t(as.data.frame(strsplit(as.character(L14[,1]), split ="-")))
L14<-data.frame(line=rep('L14',251),family=t14[,2],wght=L14[,2],dev=L14[,3])
rownames(L14) <- NULL

t1<-t(as.data.frame(strsplit(as.character(L1[,1]), split ="-")))
L1<-data.frame(line=rep('L1',241),family=t1[,2],wght=L1[,2],dev=L1[,3])
rownames(L1) <- NULL

t2<-t(as.data.frame(strsplit(as.character(L2[,1]), split ="-")))
L2<-data.frame(line=rep('L2',256),family=t2[,2],wght=L2[,2],dev=L2[,3])
rownames(L2) <- NULL

summary(L1)
summary(L2)
summary(L14)

l1f<-aggregate(L1$line,list(L1$family),length)
l2f<-aggregate(L2$line,list(L2$family),length)
l14f<-aggregate(L14$line,list(L14$family),length)

##### Mean number of individuals per family

mean(l1f[,2]) # 20.58333
sd(l1f[,2]) #16.273

mean(l2f[,2]) #36.09091
sd(l2f[,2]) #32.01392

mean(l14f[,2]) #27.33333
sd(l14f[,2]) #21.77711

##### Mean weight across families
aggregate(L1$wght,list(L1$family),mean,na.rm=T)
aggregate(L2$wght,list(L2$family),mean,na.rm=T)
aggregate(L14$wght,list(L14$family),mean,na.rm=T)

##### REML estimates of the proportion of variation in adult weight and development
##### time partitioned among families in each BC mapping population
library(lme4)
library(RLRsim)


test<-rbind(L1,L2,L14)
fm2 = lmer(wght ~ (1|line), data=test)



D<-L1

Nx<-dim(D)[2]
L<-length(unique(D[,2])) ## no. families

vars<-matrix(NA,nrow=Nx-2,ncol=3)
### first row weight, second row development time, proportion variation, statistic, p.value

for(i in 3:Nx){
    nas<-which(is.na(D[,i])==FALSE)
    dat<-data.frame(y=D[nas,i],family=as.numeric(as.factor(D[nas,2])))
    est<-lmer(y ~ 1 + (1 | family),data=dat)
    vc<-as.data.frame(VarCorr(est))
    vars[i-2,1]<-vc$vcov[1]/sum(vc$vcov)
    out<-exactRLRT(est)
    vars[i-2,3]<-out$p.value
    vars[i-2,2]<-out$statistic
}
vars
#            [,1]      [,2]   [,3]
# [1,] 0.07361105 10.418263 0.0004
# [2,] 0.06645068  2.551154 0.0374

D<-L2

Nx<-dim(D)[2]
L<-length(unique(D[,2])) ## no. families

vars<-matrix(NA,nrow=Nx-2,ncol=3)
### first row weight, second row development time, proportion variation, statistic, p.value

for(i in 3:Nx){
    nas<-which(is.na(D[,i])==FALSE)
    dat<-data.frame(y=D[nas,i],family=as.numeric(as.factor(D[nas,2])))
    est<-lmer(y ~ 1 + (1 | family),data=dat)
    vc<-as.data.frame(VarCorr(est))
    vars[i-2,1]<-vc$vcov[1]/sum(vc$vcov)
    out<-exactRLRT(est)
    vars[i-2,3]<-out$p.value
    vars[i-2,2]<-out$statistic
}

vars
#             [,1]        [,2]   [,3]
# [1,] 0.074317416 2.168706980 0.0494
# [2,] 0.001269539 0.004559697 0.3733

D<-L14

Nx<-dim(D)[2]
L<-length(unique(D[,2])) ## no. families

vars<-matrix(NA,nrow=Nx-2,ncol=3)
### first row weight, second row development time, proportion variation, statistic, p.value

for(i in 3:Nx){
    nas<-which(is.na(D[,i])==FALSE)
    dat<-data.frame(y=D[nas,i],family=as.numeric(as.factor(D[nas,2])))
    est<-lmer(y ~ 1 + (1 | family),data=dat)
    vc<-as.data.frame(VarCorr(est))
    vars[i-2,1]<-vc$vcov[1]/sum(vc$vcov)
    out<-exactRLRT(est)
    vars[i-2,3]<-out$p.value
    vars[i-2,2]<-out$statistic
}

vars
#             [,1]        [,2]   [,3]
# [1,] 0.097427554 11.11274731 0.0001
# [2,] 0.003854575  0.04428388 0.3429
