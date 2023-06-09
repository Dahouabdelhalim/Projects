library(plyr)
library(dplyr)
library(mvtnorm)
library(shapes)
library(geomorph)
library(doParallel)


registerDoParallel(cores = 4) #Registering parallel to increase speed of simulations

load("Ws.Rdata") #loading input data.
# Contains names of objects follows: 
# 1) W_/(representation)_(with or without size)- W matrices for all species on all representations.
# 2) W_bar- pooled within-group GM (with log centroid size) covariance for the full sample.
# 3) consensus-  the consensus shape for the full dataset.
# 4) msize- mean log(centroid size)
# 5) dists- vector of ild names
# 6) landmarks - vector of landmark names


# Function to calculate the relative eigenvalue variance
rEIGEN<-function(x) {
  eigenv<-eigen(x)$values
  eigenv<-eigenv[eigenv>1e-07] #excludes eigenvalues that are too small. Can be turned off
  m<-mean(eigenv)
  n<-length(eigenv)
  sqrd<-(eigenv-m)^2
  obs<-sum(sqrd)/n
  max<-(n-1)*sum(eigenv)^2/n^2
  obs/max
}
#

# Function to calculate the ild for each configuration.
# Bilaterally symmetric distances are calculated and averaged 
cdists<-function (X, dists, names=rownames(X)) {
  Dists<-matrix(NA,dim(X)[3],length(dists))
  
  dist.index<-t(matrix(unlist(strsplit(dists,"-")),2,length(dists)))
  
  l.index<- matrix(NA,length(names),2)
  l.names<-strsplit(names,"_")
  for (i in 1:length(names)) if (length(l.names[[i]])==1) l.index[i,1]<-l.names[[i]]   else  l.index[i,]<-l.names[[i]]
  
  l.inv<-apply(dist.index,1,function(x) sum(l.index[,1]==x[1]|l.index[,1]==x[2]))
  
  for(i in 1:length(dists)){
    if(l.inv[i]==2){
      l1<-which(l.index[,1]==dist.index[i,1])
      l2<-which(l.index[,1]==dist.index[i,2])
      Dists[,i]<- apply(X,3, function(x) sqrt(sum((x[l1,]-x[l2,])^2)))
    }
    
    if(l.inv[i]==3){
      l.1<-which(l.index[,1]==dist.index[i,1])
      l.2<-which(l.index[,1]==dist.index[i,2])
      if(length(l.1)==2) {
        l2<-l.1
        l1<-l.2
      }
      if(length(l.1)==1) {
        l1<-l.1
        l2<-l.2
      }
      m1<-apply(X,3, function(x) sqrt(sum((x[l1,]-x[l2[1],])^2)))
      m2<-apply(X,3, function(x) sqrt(sum((x[l1,]-x[l2[2],])^2)))
      Dists[,i]<-apply(cbind(m1,m2),1,function(x) mean(x,na.rm=T))
      
    }	
    if(l.inv[i]==4){	
      l1<-which(l.index[,1]==dist.index[i,1])
      l2<-which(l.index[,1]==dist.index[i,2])	
      
      m1<-apply(X,3, function(x) sqrt(sum((x[l1[l.index[l1,2]=="E"],]-x[l2[l.index[l2,2]=="E"],])^2)))
      m2<-apply(X,3, function(x) sqrt(sum((x[l1[l.index[l1,2]=="D"],]-x[l2[l.index[l2,2]=="D"],])^2)))
      Dists[,i]<-apply(cbind(m1,m2),1,function(x) mean(x,na.rm=T))
    }
  }
  return(Dists)
}


###################Empirical analysis#################
# Building dataframe for observed relative eigenvalue variance
rEIGEN_obs<-data.frame(gm=laply(W_gm_without_size,rEIGEN),
                       ild=laply(W_ild_without_size,rEIGEN),
                       GM=laply(W_gm_with_size,rEIGEN),
                       ILD=laply(W_ild_with_size,rEIGEN))

# Plotting empirical results
layout(matrix(c(1,2,3), 1, 3))
plot(rEIGEN_obs$gm,rEIGEN_obs$ild,
     asp=1,xlim=c(0.04,0.2),ylim=c(0.04,0.2),
     xlab = "Geometric morphometrics",
     ylab = "Interlandmark distances")
abline(a=0,b=1)
mtext("Form-With size")

plot(rEIGEN_obs$GM,rEIGEN_obs$ILD,
     asp=1,xlim=c(0.12,0.75),ylim=c(0.12,0.75),
     xlab = "Geometric morphometrics",
     ylab = "Interlandmark distances")
abline(a=0,b=1)
mtext("Shape-Without size")

plot(rEIGEN_obs$gm,rEIGEN_obs$ILD,
     asp=1,xlim=c(0.04,0.65),ylim=c(0.04,0.65),
     xlab = "Geometric morphometrics",
     ylab = "Interlandmark distances")
abline(a=0,b=1)
mtext("ShapexForm")


###################Simulation analysis#################

eigenW<-eigen(W_bar) # Eigenanalysis of pooled within-group covariance matrix

p<-seq(log(0.5),log(2.5),length.out = 20) %>% exp #stablishing scaling factors for the eigenvalues

Wsim<- alply(p, 1, function(i){                         #stablishing matrix for the simulations
  evs <- eigenW$values[1:35]^i                          #applying scaling to eigenvalues
  evs[evs<0]<-0                                         #ensuring there is no negative eingenvalues
  evs <- evs * (sum(eigenW$values[1:35])/sum(evs))      #scaling new eigenvalues to have the same sum as the original matrix
  eigenW$vectors[,1:35] %*% diag(evs[1:35]) %*% t(eigenW$vectors[,1:35])
                                                        # buiding covariance matrix
})


# Running simulations and obtained simulated values for the relative eigenvalue variance
rEIGEN_sim<-
  ldply(Wsim, function(m){
    rsamples<-rmvnorm(n = 10000, sigma = m,c(msize,t(consensus))) # obtaning all samples from a multivariate normal distribution
    sizes<-rsamples[,1]                                           # registring size
    rsamples<-arrayspecs(rsamples[,-1], dim(consensus)[1], 3)     # registring shape
    x<-adply(1:100, 1, function(i){                               # dividing the full sample into smaller samples
      samp_r<-rsamples[,,cut(1:10000,100)==levels(cut(1:10000,100))[i]]
                                                                  
      gpout<-gpagen(rsamples[,,cut(1:10000,100)==levels(cut(1:10000,100))[i]],
                    print.progress = F)                           # GPA
      samp<-gpout$coord
      sizs<-sizes[cut(1:10000,100)==levels(cut(1:10000,100))[i]]
      
      rownames(samp)<-landmarks
      r.gm<-var(two.d.array(samp))  # obtaining covariance matrix for GM shape
      # r_r.gm<-var(two.d.array(samp_r))  # obtaining covariance matrix for GM shape GM without the superposition step (can be turned on)
      r.ild<-var(cdists(samp, dists))     # obtaining covariance matrix for shape ILD 
      R.gm<-var(data.frame(sizs,two.d.array(samp))) # obtaining covariance matrix for form GM 
      # R_r.gmII<-var(data.frame(sizs,two.d.array(samp_r)))  # obtaining covariance matrix for form GM without the superposition step (can be turned on)
                       
      for(i in 1:dim(samp)[3]) {  #producing full sized configurations
        samp[,,i]<-samp[,,i]*(exp(sizs)[i]/centroid.size(samp[,,i]))
        }
      
      R.ild<-var(cdists(samp, dists)) #  obtaining covariance matrix for form ILD 
      
      data.frame(gm=rEIGEN(r.gm),
                 # gm_r=rEIGEN(r_r.gm),   #need to be turned on to calculate rEIGEN for unrotated data
                 ild=rEIGEN(r.ild),
                 GM=rEIGEN(R.gm),
                 # GM_r=rEIGEN(R_r.gmII), #need to be turned on to calculate rEIGEN for unrotated data
                 ILD=rEIGEN(R.ild))
    },.parallel = T)[,-1]
    return(x)
  })


#Plotting simulated results
layout(matrix(c(1,2,3), 1, 3))
plot(rEIGEN_sim$gm,rEIGEN_sim$ild,
     asp=1,xlim=c(0,1),ylim=c(0,1),
     xlab = "Geometric morphometrics",
     ylab = "Interlandmark distances")
abline(a=0,b=1)
mtext("Form-With size")

plot(rEIGEN_sim$GM,rEIGEN_sim$ILD,
     asp=1,xlim=c(0,1),ylim=c(0,1),
     xlab = "Geometric morphometrics",
     ylab = "Interlandmark distances")
abline(a=0,b=1)
mtext("Shape-Without size")

plot(rEIGEN_sim$gm,rEIGEN_sim$ILD,
     asp=1,xlim=c(0,1),ylim=c(0,1),
     xlab = "Geometric morphometrics",
     ylab = "Interlandmark distances")
abline(a=0,b=1)
mtext("ShapexForm")
