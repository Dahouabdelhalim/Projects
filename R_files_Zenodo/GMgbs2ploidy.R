library("gbs2ploidy")


setwd("C:/Users/Mercedes/Desktop/PopGenProject/DiploidPaper/DraftFigures/Draft2/Draft3Figures/Submission/EcoEvolSubmission")
ids<-read.csv("GMploidyID2.csv",row.names=1) 

#rows with FACS ploidy assigned
mytrn<-c(2, 8, 10, 11, 12, 21, 22, 24, 25, 27, 29, 30, 31, 32, 33, 35, 38, 42, 43, 47, 48, 53, 59, 64, 70, 72, 76, 77, 82, 89, 91, 103, 104, 105, 112)

het<-as.matrix(read.table("HetAlleledepthGM2.txt",header=T))

###allele depth file has 116 samples * 2 columns allele data per sample = 232
a<-seq(1,232,2)
b<-seq(2,232,2)
cov1<-het[,a]
cov2<-het[,b]

propOut<-estprops(cov1=cov1,cov2=cov2,props=c(0.25, 0.5, 0.75),mcmc.nchain=3,mcmc.steps=5000,mcmc.burnin=500,mcmc.thin=5)
# ## calculate observed heterozygosity and depth of coverage from the allele count
# ## data
H<-apply(is.na(cov1)==FALSE,2,mean)
D<-apply(cov1+cov2,2,mean,na.rm=TRUE)

##Training set not used
out1<-estploidy(alphas=propOut,het=H,depth=D,train=FALSE,nclasses=2,ids=rownames(ids))
write.csv(out1$pp, "noTrainploidies.csv")

##Training set used
out2<-estploidy(alphas=propOut,het=H,depth=D,train=TRUE,pl=ids$ploidy_updated,set=mytrn,nclasses=2,ids=rownames(ids))
write.csv(out2$pp, "trainploidies.csv")

##Potential triploids in data
propOutT<-estprops(cov1=cov1,cov2=cov2,props=c(0.25, 0.33, 0.5, 0.66, 0.75),mcmc.nchain=3,mcmc.steps=5000,mcmc.burnin=500,mcmc.thin=5)
out3<-estploidy(alphas=propOutT,het=H,depth=D,train=FALSE,nclasses=3,ids=rownames(ids))
write.csv(out3$pp, "possibletriploid.csv")

##All globosum tetraploid
idsG<-read.csv("allG4MploidyID.csv",row.names=1) 
#rows with FACS ploidy assigned
mytrnG<-c(1,2,7,8,10,11,12,19,21,22,23,24,25,26,27,29,30,31,32,33,35,38,42,43,47,48,49,50,53,59,60,61,63,64,69,70,71,72,74,75,76,77,81,82,83,88,89,91,92,94,95,99,100,103,104,105,106,107,109,112,113,114,116)
out4<-estploidy(alphas=propOut,het=H,depth=D,train=TRUE,pl=idsG$ploidy_updated,set=mytrnG,nclasses=2,ids=rownames(ids))
write.csv(out4$pp, "allglobotetraploid.csv")



