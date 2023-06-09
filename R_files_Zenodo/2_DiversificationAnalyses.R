##############################################################################
#This script fits the diversification models to the Western Hemisphere clades#
##############################################################################

require(ape)
require(geiger)

final.div<-read.csv("WH_CladeData.csv",header=T)
final.div$Species<-as.numeric(final.div$Species)
final.div<-final.div[-46,] #Drop Loxia sinescriurus

###############################################
###Fit a yule (pure-birth) model over lambda###
sum.ape.yule.loglik<-function(lambda,ivals,tvals){
  sum(dyule(ivals,lambda,tvals,log = TRUE))
}

yule<-optim(par=c(0.1),sum.ape.yule.loglik, gr = NULL,
      ivals=final.div$Species,tvals=final.div$Stem,
      method = "L-BFGS-B",lower = c(0), upper = c(10),control =
        list(fnscale=-1), hessian = FALSE)
yule$par #Speciation rate
yule$value #Model lnLikelihood

#plot of likelihood surface
l.vals<-seq(0.1,0.3,length.out=200)
yule.loglik<-NULL
for(i in 1:200){yule.loglik[i]<-sum.ape.yule.loglik(l.vals[i],final.div$Species,final.div$Stem)}
plot(l.vals,yule.loglik,type="l")
delta.lnl<-yule$value-qchisq(0.95,1)
abline(h=delta.lnl,lty=3)

#approximate parameter CI
l.vals[c(min(which(yule.loglik>=delta.lnl,arr.ind=T)),max(which(yule.loglik>=delta.lnl,arr.ind=T)))]

aic<-function(k,ln){(2*k)-(2*ln)}
aic(1,yule$value) #AIC value for this model

###############################################


##########################################################
###Fit the constant birth-death model for stem age data###
calc.lm<-function(r,e){
  l<-r/(1-e)
  m<-e*l
  return(c(l,m))
}

r.vals<-seq(0.001,0.25,length.out=200)
e.vals<-seq(0,1,length.out=200)
ape.loglik<-matrix(NA,ncol=200,nrow=200)
for(i in 1:200){
  for(j in 1:200){
    pars<-calc.lm(r.vals[i],e.vals[j])
    ape.loglik[i,j]<-sum(dbd(final.div$Species,pars[1],pars[2],final.div$Stem,conditional=T,log = TRUE))
  }}

sum.ape.loglik<-function(par=c(r,e),ivals,tvals){
  pars<-calc.lm(par[1],par[2])
  sum(dbd(ivals,pars[1],pars[2],tvals,conditional=TRUE,log = TRUE))
}

bd<-optim(par=c(0.1,0.01),sum.ape.loglik, gr = NULL, ivals=final.div$Species,tvals=final.div$Stem,
      method = "L-BFGS-B",lower = c(0.000001,0.0000001), upper = c(10,0.99),control = list(fnscale=-1), hessian = FALSE)
bd$par #Net diversification (r) and extinction fraction (epsilon)
bd$value #lnLikelihood

#plot of likelihood surface
r.vals<-seq(0.001,0.25,length.out=200)
e.vals<-seq(0,0.999,length.out=200)
ape.loglik<-matrix(NA,ncol=200,nrow=200)
for(i in 1:200){
  for(j in 1:200){
    pars<-calc.lm(r.vals[i],e.vals[j])
    ape.loglik[i,j]<-sum(dbd(final.div$Species,pars[1],pars[2],final.div$Stem,conditional=T,log = TRUE))
  }}
contour(r.vals,e.vals,ape.loglik,xlab="Net Diversification Rate",ylab="Extinction Fraction",nlevels = 20, levels = seq(max(ape.loglik),max(ape.loglik)-10,length.out=20))

#approximate parameter CI
delta.lnl<-bd$value-qchisq(0.95,2)/2
r.vals[c(min(which(ape.loglik>=delta.lnl,arr.ind=T)[,1]),max(which(ape.loglik>=delta.lnl,arr.ind=T)[,1]))]
e.vals[c(min(which(ape.loglik>=delta.lnl,arr.ind=T)[,2]),max(which(ape.loglik>=delta.lnl,arr.ind=T)[,2]))]

aic(2,bd$value) #AIC value for this model



###Calculate birth and death rates from optimized r and epsilon###
calc.lm(bd$par[1],bd$par[2])

##########################################################


##############################################################
###Fit the constant birth-death model for crown data##########

#subsample to only include crown data
final.div.crownonly<-final.div[is.na(final.div$Crown)==F,]

#For crown data, we are using Raup equation A18: probability of clade size i at time t given l and m, starting with any number of clades
#scaled by probability of survival (1-alpha)--necessary because we are calculating for extant clades
prob_an_scaled<-function(a,n,t,l,m){
  #parameters are a starting lineages, n ending lineages, t time since start of process, l speciation rate m extinction rate
  p0.num<-m*(exp((l-m)*t)-1)
  p0.denom<-l*exp((l-m)*t)-m
  alpha<-p0.num/p0.denom #this is p0 or probability of clade extinction
  beta<-(alpha*l)/m
  pnt<-numeric(0)
  for(j in 0:min(a,n)){
    pnt[j+1]<-choose(a,j)*choose(a+n-j-1,a-1)*alpha^(a-j)*beta^(n-j)*(1-alpha-beta)^j
  }
  return(log(sum(pnt))-log(1-alpha))
}
#vectorized version
prob_an_scaled_vector<-Vectorize(prob_an_scaled)

#wrap this function for optimization
sum.raup.loglik<-function(par=c(r,e),nvals,tvals){
  pars<-calc.lm(par[1],par[2])
  val<-sum(prob_an_scaled_vector(2,nvals,tvals,pars[1],pars[2]))
  if(!is.finite(val)) return(-999999) else return(val)
}
bd_crown<-optim(par=c(0.1,0.01),sum.raup.loglik, gr = NULL, nvals=final.div.crownonly$Species,tvals=final.div.crownonly$Crown,
                method = "L-BFGS-B",lower = c(0.000001,0.0000001), upper = c(10,0.99),control = list(fnscale=-1,trace=6,REPORT=1), hessian = FALSE)
bd_crown$par #Net diversification (r) and extinction fraction (epsilon)
bd_crown$value #lnLikelihood
aic(2,bd_crown$value) #AIC value for this model

#plot of likelihood surface
r.vals<-seq(0.001,0.25,length.out=200)
e.vals<-seq(0,0.999,length.out=200)
raup.loglik<-matrix(NA,ncol=200,nrow=200)
for(i in 1:200){
  for(j in 1:200){
    raup.loglik[i,j]<-sum.raup.loglik(nvals=final.div.crownonly$Species,tvals=final.div.crownonly$Crown,par=c(r.vals[i],e.vals[j]))
  }}
contour(r.vals,e.vals,raup.loglik,xlab="Net Diversification Rate",ylab="Extinction Fraction",nlevels = 20, levels = seq(max(raup.loglik),max(raup.loglik)-10,length.out=20))

#approximate parameter CIs
delta.lnl<-bd_crown$value-qchisq(0.95,2)/2
r.vals[c(min(which(raup.loglik>=delta.lnl,arr.ind=T)[,1]),max(which(raup.loglik>=delta.lnl,arr.ind=T)[,1]))]
e.vals[c(min(which(raup.loglik>=delta.lnl,arr.ind=T)[,2]),max(which(raup.loglik>=delta.lnl,arr.ind=T)[,2]))]

##############################################################


#######################
###Generate Figure 1###
# Fig 1a. Stem ages
pdf(file='Figure1AB.pdf',height=6,width=11,useDingbats = F)
par(mfrow=c(1,2))
plot(final.div$Stem,log(final.div$Species),xlab="Stem Age (Ma)",ylab="ln(Species Richness)",
     xlim=c(0,38),ylim=c(0,7),pch=19,cex=1.25,cex.lab=1.5,bty='l')
text(3,6.5,labels="A",cex=3)
#need to plot expectation for clade size, CONDITIONAL on clade survivorship (eq. 5 of Magallon and Sanderson)
add.curve<-function(a,r,epsilon,maxage,...){ages<-seq(0,maxage,length.out=1000);divs<-NULL;for(i in 1:length(ages)){divs[i]<-a*exp(r*ages[i])/(1-(epsilon*(exp(r*ages[i])-1)/(exp(r*ages[i])-epsilon))^a)};points(ages,log(divs),type="l",...)}
add.curve(1,0.1253577,0.6365675,38,lwd=3,col='red')
xvals<-seq(1,1000)
tvals<-seq(0,38,length.out=1000)
pvals_stem<-matrix(0,1000,1000)
for(i in 1:length(xvals)){
  for(j in 1:length(tvals)){
    pvals_stem[i,j]<-stem.p(time=tvals[j],n=xvals[i],r=0.1253577,epsilon=0.6365675)}}
contour(tvals,log(xvals),t(pvals_stem),levels=c(0.025,0.975),drawlabels=F,add=T,col="gray80")

# Fig 1b. Crown ages
plot(final.div.crownonly$Crown,log(final.div.crownonly$Species),xlab="Crown Age (Ma)",ylab="ln(Species Richness)",
     xlim=c(0,38),ylim=c(0,7),pch=19,cex=1.25,cex.lab=1.5,bty='l')
text(3,6.5,labels="B",cex=3)
add.curve(2,0.1305868,0.8463151,38,lwd=3,col='red')
pvals_crown<-matrix(0,1000,1000)
for(i in 1:length(xvals)){
  for(j in 1:length(tvals)){
    pvals_crown[i,j]<-crown.p(time=tvals[j],n=xvals[i],r=0.1305868,epsilon=0.8463151)}}
contour(tvals,log(xvals),t(pvals_crown),levels=c(0.025,0.975),drawlabels=F,add=T,col="gray80")

dev.off()
par(mfrow=c(1,1))
#######################


#################################################################
###Fit the model with decelerating speciation and 0 extinction###
sum.dbdTime_noext.loglik<-function(par=c(l0,k),ivals,tvals){
  my.exp.lambda<-function(t,l0=par[1],k=par[2]){
    l0*exp(-k*t)
  }
  my.const.mu<-function(t,mu=0){
    mu
  }
  logliks<-NULL
  for(i in 1:length(ivals)){
    logliks[i]<-log(dbdTime(ivals[i],my.exp.lambda,my.const.mu,tvals[i],conditional=TRUE))
  }
  sum(logliks)
}

ds.noe<-optim(par=c(0.1,0.01),sum.dbdTime_noext.loglik, gr = NULL, ivals=final.div$Species,tvals=final.div$Stem,
      method = "L-BFGS-B",lower = c(0.000001,0.0000001), upper = c(10,0.99),control = list(fnscale=-1), hessian = FALSE)
ds.noe$par #Initial speciation rate (lamdba), and speciation decay constant (K)
ds.noe$value

aic(2,ds.noe$value) #AIC value for this model

#plot of likelihood surface
l0.vals<-seq(0.001,0.5,length.out=200)
k.vals<-seq(0,0.2,length.out=200)
ds.loglik<-matrix(NA,ncol=200,nrow=200)
for(i in 1:200){
  for(j in 1:200){
    ds.loglik[i,j]<-sum.dbdTime_noext.loglik(par=c(l0.vals[i],k.vals[j]),ivals=final.div$Species,tvals=final.div$Stem)
  }}
contour(l0.vals,k.vals,ds.loglik,xlab="Initial Diversification Rate",ylab="Rate of Diversification Decline",nlevels = 20, levels = seq(max(ds.loglik),max(ds.loglik)-10,length.out=20))

#approximate parameter CIs
delta.lnl<-ds.noe$value-qchisq(0.95,2)/2
l0.vals[c(min(which(ds.loglik>=delta.lnl,arr.ind=T)[,1]),max(which(ds.loglik>=delta.lnl,arr.ind=T)[,1]))]
k.vals[c(min(which(ds.loglik>=delta.lnl,arr.ind=T)[,2]),max(which(ds.loglik>=delta.lnl,arr.ind=T)[,2]))]

#################################################################


########################################################################
###Fit the model with decelerating speciation and constant extinction###
sum.dbdTime.loglik<-function(par=c(l0,k,e),ivals,tvals){
  par.trans<-calc.lm(par[1],par[3])
  my.exp.lambda<-function(t,l0=par.trans[1],k=par[2]){
    l0*exp(-k*t)
  }
  my.const.mu<-function(t,mu=par.trans[2]){
    mu
  }
  logliks<-NULL
  for(i in 1:length(ivals)){
    logliks[i]<-log(dbdTime(ivals[i],my.exp.lambda,my.const.mu,tvals[i],conditional=TRUE))
  }
  sum(logliks)
}

ds.e<-optim(par=c(0.1,0.01,0.01),sum.dbdTime.loglik, gr = NULL, ivals=final.div$Species,tvals=final.div$Stem,
      method = "L-BFGS-B",lower = c(0.000001,0.0000001,0), upper = c(10,0.99,0.99),control = list(fnscale=-1), hessian = FALSE)
ds.e$par #Net diversification rate (r), speciation decay constant (K), and extinction fraction (epsilon)

aic(3,ds.e$value) #AIC value for this model

#plot of likelihood surface
l0.vals<-seq(0,0.5,length.out=50)
k.vals<-seq(0,0.1,length.out=50)
e.vals<-seq(0,0.90,length.out=50)
ds.e.loglik<-array(NA,dim=c(50,50,50))
for(i in 1:50){
  for(j in 1:50){
    for(k in 1:50){
      ds.e.loglik[i,j,k]<-sum.dbdTime.loglik(par=c(l0.vals[i],k.vals[j],e.vals[k]),ivals=final.div$Species,tvals=final.div$Stem)
    }}}
contour(l0.vals,k.vals,ds.e.loglik[,,21],xlab="Initial Diversification Rate",ylab="Rate of Diversification Decline",levels = delta.lnl); points(0.18771166,0.01255653);
contourLines(l0.vals,k.vals,ds.e.loglik[,,21],levels = delta.lnl)

#approximate parameter CIs
delta.lnl<-ds.e$value-qchisq(0.95,3)/2
l0.vals[c(min(which(ds.e.loglik>=delta.lnl,arr.ind=T)[,1]),max(which(ds.e.loglik>=delta.lnl,arr.ind=T)[,1]))]
k.vals[c(min(which(ds.e.loglik>=delta.lnl,arr.ind=T)[,2]),max(which(ds.e.loglik>=delta.lnl,arr.ind=T)[,2]))]
e.vals[c(min(which(ds.e.loglik>=delta.lnl,arr.ind=T)[,3]),max(which(ds.e.loglik>=delta.lnl,arr.ind=T)[,3]))]

###Calculate initial speciation rate and constant extinction rates###
calc.lm(ds.e$par[1],ds.e$par[3])

########################################################################


#################################################
###Calculate expected diversity for 3 outliers###
require(TreeSim)

# Use estimates from uniform BD model
lambda<-0.3449271
mu<-0.2195694

# Emberizoidea
emb.sim<-vector()
emb.tree<-pbtree(lambda,mu,nsim=1000,t=final.div$Stem[1])
for(i in 1:1000){
  emb.sim[i]<-Ntip(emb.tree[[i]])
}
quantile(emb.sim,probs=c(0.025,0.975))
median(emb.sim)

# Turdus
turd.sim<-vector()
turd.tree<-pbtree(lambda,mu,nsim=1000,t=final.div$Stem[8])
for(i in 1:1000){
  turd.sim[i]<-Ntip(turd.tree[[i]])
}
quantile(turd.sim,probs=c(0.025,0.975))
median(turd.sim)

# Ptiliogonatidae
ptilio.sim<-vector()
ptilio.tree<-pbtree(lambda,mu,nsim=1000,t=final.div$Stem[17])
for(i in 1:1000){
  ptilio.sim[i]<-Ntip(ptilio.tree[[i]])
}
quantile(ptilio.sim,probs=c(0.025,0.975))

median(ptilio.sim)

#################################################


###############################################################
###Estimate diversification parameters across tree distribution

# Read in stem age distributions
stem.dist<-read.csv("StemAge_Dist.csv")
rownames(stem.dist)<-stem.dist$Clade
stem.dist<-stem.dist[complete.cases(stem.dist),-1]

dist.rich<-final.div[final.div$Lineage%in%rownames(stem.dist)==T,4]

## Pure-birth distribution
yule.dist<-vector(length=100)

for(i in 1:1000){
  tryCatch({
    yule<-optim(par=c(0.1),sum.ape.yule.loglik, gr = NULL,
            ivals=dist.rich,tvals=stem.dist[,i],
            method = "L-BFGS-B",lower = c(0), upper = c(10),control =
              list(fnscale=-1), hessian = FALSE)
    yule.dist[i]<-yule$par
  }, error=function(e){})
}
yule.dist<-yule.dist[yule.dist>0]
yule.dist<-na.om

#plot it
pdf("PB_Distribution.pdf",height=4,width=4,useDingbats = F)
plot(density(yule.dist),col='black',
     xlab=expression(paste("Pure-birth: ",lambda)),main=NA)
polygon(density(yule.dist),col=rgb(0,0,0,0.25))
abline(v=median(yule.dist),lwd=2,lty=2)
abline(v=0.1749,lwd=2,col='red')
dev.off()


## Birth-death distribution
bd.dist<-matrix(nrow=1000,ncol=2,dimnames=list(NULL,c("Lambda","Mu")))

for(i in 1:1000){
  bd<-optim(par=c(0.1,0.01),sum.ape.loglik, gr = NULL, ivals=dist.rich,tvals=stem.dist[,i],
            method = "L-BFGS-B",lower = c(0.000001,0.0000001), upper = c(10,0.99),control = list(fnscale=-1), hessian = FALSE)
  bd.dist[i,]<-calc.lm(bd$par[1],bd$par[2])
}

#plot it
pdf("BD_Distribution.pdf",height=4,width=4,useDingbats = F)
#lambda
plot(density(bd.dist[,1]),col='black',
     xlab=expression(paste("Birth-death: ",lambda)),main=NA)
polygon(density(bd.dist[,1]),col=rgb(0,0,0,0.25))
abline(v=median(bd.dist[,1]),lwd=2,lty=2)
abline(v=0.3449,lwd=2,col='red')
#mu
plot(density(bd.dist[,2]),col='black',
     xlab=expression(paste("Birth-death: ",mu)),main=NA)
polygon(density(bd.dist[,2]),col=rgb(0,0,0,0.25))
abline(v=median(bd.dist[,2]),lwd=2,lty=2)
abline(v=0.2196,lwd=2,col='red')
dev.off()

## Decelerating speciation distribution
ds.dist<-matrix(nrow=1000,ncol=2,dimnames=list(NULL,c("Lambda","K")))

for(i in 1:1000){
  ds.noe<-optim(par=c(0.1,0.01),sum.dbdTime_noext.loglik, gr = NULL, ivals=dist.rich,tvals=stem.dist[,i],
                method = "L-BFGS-B",lower = c(0.000001,0.0000001), upper = c(10,0.99),control = list(fnscale=-1), hessian = FALSE)
  ds.dist[i,]<-ds.noe$par
}

#plot it
#lambda
pdf("DS_Distribution.pdf",height=4,width=4,useDingbats = F)
plot(density(ds.dist[,1]),col='black',
     xlab=expression(paste("Decelerating Speciation: ",lambda)),main=NA)
polygon(density(ds.dist[,1]),col=rgb(0,0,0,0.25))
abline(v=median(ds.dist[,1]),lwd=2,lty=2)
abline(v=0.2491,lwd=2,col='red')
#k
plot(density(ds.dist[,2]),col='black',xlab="Decelerating Speciation: K",main=NA)
polygon(density(ds.dist[,2]),col=rgb(0,0,0,0.25))
abline(v=median(ds.dist[,2]),lwd=2,lty=2)
abline(v=0.0349,lwd=2,col='red')
dev.off()


## Decelerating speciation w/ extinction
dse.dist<-matrix(nrow=1000,ncol=3,dimnames=list(NULL,c("Lambda","Mu","K")))
  
for(i in 1:1000){
  tryCatch({
    ds.e<-optim(par=c(0.1,0.01,0.01),sum.dbdTime.loglik, gr = NULL, ivals=dist.rich,tvals=stem.dist[,i],
                method = "L-BFGS-B",lower = c(0.000001,0.0000001,0), upper = c(10,0.99,0.99),control = list(fnscale=-1), hessian = FALSE)
    dse.dist[i,1:2]<-calc.lm(ds.e$par[1],ds.e$par[3])
    dse.dist[i,3]<-ds.e$par[2]
  }, error=function(e){})
}
save.image("Diversification_Distributions.RData")

dse.dist<-dse.dist[complete.cases(dse.dist),]

#plot it
#lambda
pdf("DSE_Distribution.pdf",height=4,width=4,useDingbats = F)
plot(density(dse.dist[,1]),col='black',
     xlab=expression(paste("Decelerating Speciation with Extinction: ",lambda)),main=NA)
polygon(density(dse.dist[,1]),col=rgb(0,0,0,0.25))
abline(v=median(dse.dist[,1]),lwd=2,lty=2)
abline(v=0.2955,lwd=2,col='red')
#mu
plot(density(dse.dist[,2]),col='black',
     xlab=expression(paste("Decelerating Speciation with Extinction: ",mu)),main=NA)
polygon(density(dse.dist[,2]),col=rgb(0,0,0,0.25))
abline(v=median(dse.dist[,2]),lwd=2,lty=2)
abline(v=0.1078,lwd=2,col='red')
#k
plot(density(dse.dist[,3]),col='black',xlab="Decelerating Speciation with Extinction: K",main=NA)
polygon(density(dse.dist[,3]),col=rgb(0,0,0,0.25))
abline(v=median(dse.dist[,3]),lwd=2,lty=2)
abline(v=0.0126,lwd=2,col='red')
dev.off()

###############################################################