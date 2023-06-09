### R script for stochastic character state mapping

library(phytools)

## tr is the tree, from read.tree
tr<-read.tree("mmc3.txt")

#### resistance mutation data ##############

## val is the matrix of character prior probs. 0 vs. 1 when known, 0.5 0.5 when not
val<-read.table("charMatrixCor.txt",row.names=1,header=FALSE)

#### only taxa we have data for
missing<-row.names(val)[which(is.na(val[,1]))]
subtr<-drop.tip(phy=tr,tip=missing)
subval<-val[which(is.na(val[,1])==FALSE),]

write.tree(phy=subtr,file="subTree.tree",tree.names=TRUE)
write.table(subval,file="subVal.txt",row.names=T,col.names=F,quote=F)

sv<-subval[,1]
names(sv)<-rownames(subval)

substochout1<-make.simmap(subtr,x=sv,model="ARD",nsim=2000,Q="mcmc",burnin=10000,samplefreq=50)
submapsum1<-describe.simmap(substochout1,plot=FALSE)

substochout2<-make.simmap(subtr,x=sv,model="ARD",nsim=2000,Q="mcmc",burnin=10000,samplefreq=50)
submapsum2<-describe.simmap(substochout2,plot=FALSE)

####### same as above, but with different codings for the two different nucleotide states at codon 120, G vs. C ##################

val<-read.table("charMatrixAA.txt",row.names=1,header=FALSE)

#### only taxa we have data for
missing<-row.names(val)[which(is.na(val[,1]))]
subtr<-drop.tip(phy=tr,tip=missing)
subval<-val[which(is.na(val[,2])==FALSE),]

write.tree(phy=subtr,file="subTree.tree",tree.names=TRUE)
write.table(subval,file="subVal.txt",row.names=T,col.names=F,quote=F)

sv<-subval[,2]
names(sv)<-rownames(subval)

substochouta1<-make.simmap(subtr,x=sv,model="ARD",nsim=2000,Q="mcmc",burnin=10000,samplefreq=50)
submapsuma1<-describe.simmap(substochout1,plot=FALSE)

substochouta2<-make.simmap(subtr,x=sv,model="ARD",nsim=2000,Q="mcmc",burnin=10000,samplefreq=50)
submapsuma2<-describe.simmap(substochout2,plot=FALSE)


