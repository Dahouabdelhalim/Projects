#Sample analysis for Nolana using the genetic distance matrix insead of phylogeny

#load libraries that you will need
library(MCMCglmm)
library(MASS)


#Read genetic distance matrix and assign names
gd<-read.csv("nolanagd.csv",header=F)
mat<-as.matrix(gd[,2:12])


#Take generalized inverse and format for MCMCglmm
#generalized inverse of 1-gd matrix is required for the analysis
M1<-1-mat
M1a<- as(ginv(M1), "dgCMatrix")
dimnames(M1a)<-list(gd[,1],gd[,1])
#might get warning about supressing column names


#If using phylogeny instead
#Take generalized inverse and format for MCMCglmm
tree<-read.tree("nolana_ultrametric_tree.txt")
sp1inv<-inverseA(tree,nodes="TIPS")$Ainv
sp2inv<-inverseA(tree,nodes="TIPS")$Ainv
#use these matrices in place of M1a below


#Trait data
#Names in dataset need to match names for genetic distance matrix and/or phylogeny exactly
data<-read.csv("nolana.csv")

#Creating subsets of the data
#There should be no missing values for RI or the explanatory variables
total<-subset(data,!is.na(data$Total_Isolation))

#setting up uninformative priors based on y variable of interest
p.var<-var(total$Total_Isolation)

#Model potential correlation between fixed effects by making variance-covariance matrix
Bvar<-diag(4)*p.var
Bvar[Bvar==0]<-p.var/4

#prior for correlated fixed effects, residuals and two random effects "G"

prior1=list(B=list(mu=rep(0,4),V=Bvar),R=list(V=p.var/3, n=1), G=list(G1=list(V=p.var/3, n=1),G2=list(V=p.var/3,n=1)))


#Run two MCMC chains for each analysis so that we can look at diagnostics
#set seed so results can be reproduced
set.seed(1001)
m1a_1<-MCMCglmm(Total_Isolation~GeneticDistance+GeographicDistance+CorolDiam_Diff,random=~Species1+Species2,ginverse=list(Species1=M1a,Species2=M1a),prior=prior1,data=total,nitt=26000,thin=20,burnin=6000,verbose=F)

set.seed(3001)
m1a_2<-MCMCglmm(Total_Isolation~GeneticDistance+GeographicDistance+CorolDiam_Diff,random=~Species1+Species2,ginverse=list(Species1=M1a,Species2=M1a),prior=prior1,data=total,nitt=26000,thin=20,burnin=6000,verbose=F)

#Diagnostics based on gelman-rubin approach
#looking for values to be close to 1, values less than 1.1 acceptable
gelman.diag(mcmc.list(m1a_1$Sol,m1a_2$Sol))
gelman.diag(mcmc.list(m1a_1$VCV,m1a_2$VCV))

#Adding the two chains together to profile posterior distribution
m1asum<-as.mcmc(rbind(m1a_1$Sol,m1a_2$Sol))
HPDinterval(m1asum)


