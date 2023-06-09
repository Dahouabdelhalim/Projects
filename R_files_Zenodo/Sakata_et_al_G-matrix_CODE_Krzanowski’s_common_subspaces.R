#######################################################################################################################
This R code corresponds to the following article:		
#Title: 		Environmentally triggered variability in the genetic variance-covariance of herbivory resistance of an exotic plant Solidago altissima
#Authors: 	Yuzu Sakata, Shunsuke Utsumi, Timothy Craig, Joanne Itami, Mito Ikemoto, Takayuki Ohgushi	
#Journal: 	Ecology and Evolution

Below is the R codes for calculating G-matrix from plant resistance data from four groups and compare G-matrices using Krzanowski’s common subspaces following the code of Aguirre et al.(2014). 

For details please see supplementary tutorial in Aguirre et al. (2014) Comparing G: Multivariate Analysis of Genetic Variation in Multiple Populations. Heredity 112:21-29.

#######################################################################################################################

library(MCMCglmm)
library(QGglmm)
library(gdata);library(matrixcalc)

#read data#
data_jj<-read.csv("japan_j.csv")#Japanese plants in Japanese garden#
data_ju<-read.csv("japan_u.csv")#Japanese plants in USA garden#
data_uj<-read.csv("us_j.csv")#USA plants in Japanese garden#
data_uu<-read.csv("us_u.csv")#USA plants in USA garden#


#read pedigree data#
Ped_jj<-read.csv("Ped_japan_j.csv")
Ped_ju<-read.csv("Ped_japan_u.csv")
Ped_uj<-read.csv("Ped_us_j.csv")
Ped_uu<-read.csv("Ped_us_u.csv")


＃logarithm of the total number of leaves#
data_jj$logall7<-log(data_jj$all7)

#standardize
data_jj$s_damage7<-decostand(data_jj$damage7,method="max")
write.csv(data_jj,"japan_j_s.csv")

data_ju$r_damage7<-data_ju$damage7/data_ju$all7
data_ju$r_other7<-data_ju$other7/data_ju$all7

#animal model#data_jj for example#
phen.var<-matrix(c(var(data_jj$damage7,na.rm=T),0,0,var(data_jj$other7,na.rm=T)),2,2)

prioroffset<-list(G=list(G1=list(V=phen.var/4,n=2),G2=list(V=phen.var/4,n=2),G3=list(V=phen.var/4,n=2)),R=list(V=phen.var/4,n=2),B=list(mu=c(0,0,1,1),V=diag(4)*c(1e7,1e7,1e-7,1e-7)))

prior<-list(R=list(V=diag(2)/2,nu=2),G=list(G1=list(V=diag(2)/2,nu=2),G2=list(V=diag(2)/2,nu=2),G3=list(V=diag(2)/2,nu=2)))#offsetなし


#Bayesian bivariate model to obtain VCV including the G-matrix (Japanese plants in Japanese garden as an example)#
model<- MCMCglmm(cbind(damage7,other7)~trait-1+trait:logall7,random=~us(trait):animal+idh(trait):garden+idh(trait):population,rcov=~us(trait):units,data=data_jj,pedigree=Ped_jj,prior=prioroffset,family=c("poisson","poisson"),nitt=65000,thin=50000,burnin=15000)

#output VCV#
write.csv(model$VCV,"jj_vcv.csv")
posterior.mode(model$VCV)
HPDinterval(model$VCV,0.95)
plot(model$Sol)

#obtain the G-matrix in the observed scale (post)(Japanese plants in Japanese garden as an example)#
yhat <- predict(model,type="terms")
yhat2<-matrix(yhat,nrow=392,ncol=2)
df <- cbind(mu1=model[["Sol"]][,"traitdamage7"],mu2=model[["Sol"]][,"traitother7"],model[["VCV"]])
post <- apply(df,1,function(row){
  G <- matrix(c(row["traitdamage7:traitdamage7.animal"],row["traitother7:traitdamage7.animal"],row["traitdamage7:traitother7.animal"],row["traitother7:traitother7.animal"]),ncol=2)
  R <- matrix(c(row["traitdamage7:traitdamage7.units"],row["traitother7:traitdamage7.units"],row["traitdamage7:traitother7.units"],row["traitother7:traitother7.units"]),ncol=2)
  P <- G+R
  QGmvparams(predict=yhat2,vcv.G=G,vcv.P=P,
             model=c("Poisson.log","Poisson.log"),verbose=FALSE)
})






#read "post" data#
post_jj<-read.csv("post_jj.csv")#Japanese plants in Japanese gardens#
post_ju<-read.csv("post_ju.csv")#Japanese plants in USA gardens#
post_uj<-read.csv("post_uj.csv")#USA plants in Japanese gardens#
post_uu<-read.csv("post_uu.csv")#USA plants in USA gardens#

#Construct Garray (USA plants in USA gardens and USA plants in Japanese gardens for example)#
MCMCarray <- array(,c(4,10,2))
MCMCarray[,,1] <- as.matrix(post_uj)
MCMCarray[,,2] <- as.matrix(post_uu)
MCMCarray2 <- array(,c(10,4,2))
MCMCarray2[,,1] <- as.matrix(t(MCMCarray[,,1]))
MCMCarray2[,,2] <- as.matrix(t(MCMCarray[,,2]))
Garray <- array(,c(2,2,2,10))
dimnames(Garray) <- list(c("damage","other"),c("damage","other"),c("uj","uu"))
for (i in 1:2){
  for (j in 1:10){
    G <- matrix(MCMCarray2[j,1:(2^2),i],ncol= 2)
    Garray[,,i,j] <- G}}
    

    
#Construct randamized Garray#
rand.Garray <- array(,c(2,2,2,10))
dimnames(rand.Garray) <- list(c("damage","other"),c("damage","other"),c("uj","uu"))
for (i in 1:10){
  uj.bv<-rbv(Ped_uj,Garray[,,1,i])
  uu.bv<-rbv(Ped_uu,Garray[,,2,i])
  a.pop <- cumsum(c(dim(Ped_uj)[1],dim(Ped_uu)[1]))
  pop.bv <- rbind(uj.bv,uu.bv)
  rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1],replace=F),]
  rand.Garray[,,1,i] <- cov(rand.pop.bv[1:a.pop[1],])
  rand.Garray[,,2,i] <- cov(rand.pop.bv[(a.pop[1] + 1):a.pop[2],])
  
}
 

#Krzanowski’s common subspaces#
MCMCG.kr <- kr.subspace(Garray, vec = rep(2/2,2))
MCMCG.kr.rand <- kr.subspace(rand.Garray, vec = rep(2/2,2))

#obtain HPD interval#
HPD.H.val <- cbind(HPDinterval(as.mcmc(MCMCG.kr$MCMC.H.val)),HPDinterval(as.mcmc(MCMCG.kr.rand$MCMC.H.val)))



