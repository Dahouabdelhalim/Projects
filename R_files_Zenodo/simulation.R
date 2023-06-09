# There are two different habitats A and B
# Male habitat preference locus M/female habitat preference locus F
# In the absence of sex limited gene expression for hab. pref., M determines hab. pref. for both sexes
# Imaginary number 0.5i corresponds to the "no-preference" allele

###library###
library(foreach)
library(doParallel)

##define functions #####
#######################

#x:vector (D+E), x=1 -> random sp. assign, x < 1 -> sp1:0, x > 1 -> sp2:1�@
#NOTE: this function is virtually doing nothing in the current version of the model
sp.identifier <- function(x) ifelse(x==1,rbinom(length(x),1,0.5),floor(2*x/3))

#x,y,z:vectors (M,F,sex) #other param: sex.limitation, habitat.pref
#When "sex.limitation <- 0", M determines the habitat choice of both sexes. Othewise M and F determine male and female
#habitat, respectively. Preferred habitat is preferred "habitat.pref.x" times more than the other one. i.e.habitat.pref.x=1 means no preference.
habitat.determiner <- function(x,y,z) {
  pot.hab  <- cbind(sample(c(0,1),length(x),replace=TRUE,prob=c(habitat.pref,1)), #habitats for M/F=0  +0i
                    sample(c(0,1),length(x),replace=TRUE,prob=c((habitat.pref+1)/2,1)), #             M/F=0  +0.5i
                    sample(c(0,1),length(x),replace=TRUE,prob=c(1,1)),                #             M/F=0.5+0i/0+1i
                    sample(c(0,1),length(x),replace=TRUE,prob=c(1,(habitat.pref+1)/2)), #             M/F=0.5+0.5i
                    sample(c(0,1),length(x),replace=TRUE,prob=c(1,habitat.pref)))     #             M/F=1  +0i
  ge <- if (sex.limitation==0) {cbind(ifelse(x==0,1,0),ifelse(x==0.5i,1,0),ifelse(x==0.5|x==1i,1,0),ifelse(x==0.5+0.5i,1,0),ifelse(x==1,1,0), #gene expression
                                      rep(0,length(x)),rep(0,length(x)),rep(0,length(x)),rep(0,length(x)),rep(0,length(x)))} else {
                                cbind(ifelse(z==1&x==0,1,0),ifelse(z==1&x==0.5i,1,0),ifelse(z==1&(x==0.5|x==1i),1,0),ifelse(z==1&x==0.5+0.5i,1,0),ifelse(z==1&x==1,1,0),
                                      ifelse(z==0&y==0,1,0),ifelse(z==0&y==0.5i,1,0),ifelse(z==0&(y==0.5|y==1i),1,0),ifelse(z==0&y==0.5+0.5i,1,0),ifelse(z==0&y==1,1,0))}
  habitat <- pot.hab[,1]*ge[,1]+pot.hab[,2]*ge[,2]+pot.hab[,3]*ge[,3]+pot.hab[,4]*ge[,4]+pot.hab[,5]*ge[,5]+
             pot.hab[,1]*ge[,6]+pot.hab[,2]*ge[,7]+pot.hab[,3]*ge[,8]+pot.hab[,4]*ge[,9]+pot.hab[,5]*ge[,10]
  habitat
}

#x,y:2data.frames (male,female) #other param:evolve.sp.recog, mate.pref.1, mate.pref.2
non.MAD.mating.matcher <- function(x,y) {
  z <- y[,"sec.hab"]
  fec <- floor((1-z)*2*R.inA + z*2*R.inB) + ifelse(z==0,rbinom(nrow(y),1,2*R.inA-floor(2*R.inA)),
                                                        rbinom(nrow(y),1,2*R.inB-floor(2*R.inB))) #determine female fecundity
  mom.genotype <- y[rep(1:nrow(y),fec),c(1:6,9,10)] #genotype(M,F,T,P) and mother's sp and sec.hab
  if (evolve.sp.recognition==1) {        #mating algorithm for models where mate choice evolves
    pref.matrix <- cbind(ifelse(x[,"T"]==0, mate.pref.1,     ifelse(x[,"T"]==0.5i,(mate.pref.1+1)/2,1)),  #male attractiveness for P=0        female of sp.1
                         ifelse(x[,"T"]==0,(mate.pref.1+1)/2,ifelse(x[,"T"]==0.5i,(mate.pref.1+3)/4,1)),                      #for P=0.5i     female of sp.1
                         rep(1,nrow(x)),                                                            #for all neutral female
                         ifelse(x[,"T"]==1,(mate.pref.1+1)/2,ifelse(x[,"T"]==0.5+0.5i,(mate.pref.1+3)/4,1)),                  #for P=0.5+0.5i female of sp.1
                         ifelse(x[,"T"]==1, mate.pref.1,     ifelse(x[,"T"]==0.5+0.5i,(mate.pref.1+1)/2,1)),                  #for P=1        female of sp.1
                         ifelse(x[,"T"]==0, mate.pref.2,     ifelse(x[,"T"]==0.5i,(mate.pref.2+1)/2,1)),                      #for P=0        female of sp.2
                         ifelse(x[,"T"]==0,(mate.pref.2+1)/2,ifelse(x[,"T"]==0.5i,(mate.pref.2+3)/4,1)),                      #for P=0.5i     female of sp.2
                         ifelse(x[,"T"]==1,(mate.pref.2+1)/2,ifelse(x[,"T"]==0.5+0.5i,(mate.pref.2+3)/4,1)),                  #for P=0.5+0.5i female of sp.2
                         ifelse(x[,"T"]==1 ,mate.pref.2,     ifelse(x[,"T"]==0.5+0.5i,(mate.pref.2+1)/2,1)))                  #for P=1        female of sp.2
    potential.dad.id <- cbind(sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,1]), #for P=0        female of sp.1
                              sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,2]), #for P=0.5i     female of sp.1
                              sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,3]), #for all neutral female
                              sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,4]), #for P=0.5+0.5i female of sp.1
                              sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,5]), #for P=1        female of sp.1
                              sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,6]), #for P=0        female of sp.2
                              sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,7]), #for P=0.5i     female of sp.2
                              sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,8]), #for P=0.5+0.5i female of sp.2
                              sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,9])) #for P=1        female of sp.2
    dad.id <- potential.dad.id[,1]*ifelse(mom.genotype[,"sp"]==0&mom.genotype[,"P"]==0,       1,0)+
              potential.dad.id[,2]*ifelse(mom.genotype[,"sp"]==0&mom.genotype[,"P"]==0.5i,    1,0)+
              potential.dad.id[,3]*ifelse(mom.genotype[,"P"]==0.5|mom.genotype[,"P"]==1i,     1,0)+
              potential.dad.id[,4]*ifelse(mom.genotype[,"sp"]==0&mom.genotype[,"P"]==0.5+0.5i,1,0)+
              potential.dad.id[,5]*ifelse(mom.genotype[,"sp"]==0&mom.genotype[,"P"]==1,       1,0)+
              potential.dad.id[,6]*ifelse(mom.genotype[,"sp"]==1&mom.genotype[,"P"]==0,       1,0)+
              potential.dad.id[,7]*ifelse(mom.genotype[,"sp"]==1&mom.genotype[,"P"]==0.5i,    1,0)+
              potential.dad.id[,8]*ifelse(mom.genotype[,"sp"]==1&mom.genotype[,"P"]==0.5+0.5i,1,0)+
              potential.dad.id[,9]*ifelse(mom.genotype[,"sp"]==1&mom.genotype[,"P"]==1,       1,0)
    dad.genotype <- x[dad.id,1:6]
    egg <- as.data.frame(cbind(0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,1])), #Re() only yet
                               0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,2])),
                               0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,3])),
                               0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,4])),
                               0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,5])),
                               0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,6]))))
    sperm <- as.data.frame(cbind(0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,1])), #Re() only yet
                                 0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,2])),
                                 0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,3])),
                                 0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,4])),
                                 0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,5])),
                                 0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,6]))))
    colnames(egg) <- c("M","F","T","P","D","E"); colnames(sperm) <- c("M","F","T","P","D","E")
    egg$M <- egg$M+ifelse(egg$M==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,1])/(1-Re(mom.genotype[,1])))) #add Im() component
    egg$F <- egg$F+ifelse(egg$F==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,2])/(1-Re(mom.genotype[,2]))))  
    egg$T <- egg$T+ifelse(egg$T==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,3])/(1-Re(mom.genotype[,3]))))
    egg$P <- egg$P+ifelse(egg$P==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,4])/(1-Re(mom.genotype[,4]))))
    sperm$M <- sperm$M+ifelse(sperm$M==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,1])/(1-Re(dad.genotype[,1])))) #add Im() component
    sperm$F <- sperm$F+ifelse(sperm$F==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,2])/(1-Re(dad.genotype[,2]))))  
    sperm$T <- sperm$T+ifelse(sperm$T==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,3])/(1-Re(dad.genotype[,3]))))
    sperm$P <- sperm$P+ifelse(sperm$P==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,4])/(1-Re(dad.genotype[,4]))))
    offspring <- egg+sperm
  } else {                           #mating algorithm for models where mate choice does not evolve
    pref.matrix <- cbind(ifelse(x[,"sp"]==0,mate.pref.1,1),ifelse(x[,"sp"]==1,mate.pref.2,1))     #male attractiveness for females of sp.1 and sp.2, respectively
    potential.dad.id <- cbind(sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,1]), #for sp.1 female
                              sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,2])) #for sp.2 female
    dad.id <- potential.dad.id[,1]*(1-mom.genotype[,"sp"])+potential.dad.id[,2]*mom.genotype[,"sp"]
    dad.genotype <- x[dad.id,1:6]
    egg <- as.data.frame(cbind(0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,1])), #Re() only yet
                               0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,2])),
                               0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,3])),
                               0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,4])),
                               0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,5])),
                               0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,6]))))
    sperm <- as.data.frame(cbind(0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,1])), #Re() only yet
                                 0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,2])),
                                 0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,3])),
                                 0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,4])),
                                 0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,5])),
                                 0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,6]))))
    colnames(egg) <- c("M","F","T","P","D","E"); colnames(sperm) <- c("M","F","T","P","D","E")
    egg$M <- egg$M+ifelse(egg$M==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,1])/(1-Re(mom.genotype[,1])))) #add Im() component
    egg$F <- egg$F+ifelse(egg$F==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,2])/(1-Re(mom.genotype[,2]))))  
    egg$T <- egg$T+ifelse(egg$T==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,3])/(1-Re(mom.genotype[,3]))))
    egg$P <- egg$P+ifelse(egg$P==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,4])/(1-Re(mom.genotype[,4]))))
    sperm$M <- sperm$M+ifelse(sperm$M==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,1])/(1-Re(dad.genotype[,1])))) #add Im() component
    sperm$F <- sperm$F+ifelse(sperm$F==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,2])/(1-Re(dad.genotype[,2]))))  
    sperm$T <- sperm$T+ifelse(sperm$T==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,3])/(1-Re(dad.genotype[,3]))))
    sperm$P <- sperm$P+ifelse(sperm$P==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,4])/(1-Re(dad.genotype[,4]))))
    offspring <- egg+sperm
  }
  offspring$nat.hab <- mom.genotype[,"sec.hab"]
  offspring
}

#x,y:vector (D,E), x+y=0 or x*y=1 -> 1, otherwise randomly 0 or 1 with prob=s.hybrid
fitness.assigner <- function(x,y) ifelse(x+y==0|x*y==1, rep(1,length(x)), rbinom(length(x),1,s.hybrid))

#x,y,z:2data.frames (male,female) and a [scaler] of habitat(0,1) #other param:evolve.sp.recog, mate.pref.1, mate.pref.2, neutral.f
MAD.mating.matcher <- function(x,y,z) {
   fec <- floor((1-z)*2*R.inA + z*2*R.inB) + rbinom(nrow(y),1,(1-z)*2*R.inA + z*2*R.inB-floor((1-z)*2*R.inA + z*2*R.inB)) #determine female fecundity
   mom.genotype <- y[rep(1:nrow(y),fec),c(1:6,9)] #mother's genotype + sp
      if (evolve.sp.recognition==1) {        #mating algorithm for models where mate choice evolves
        pref.matrix <- cbind(ifelse(x[,"T"]==0, mate.pref.1,     ifelse(x[,"T"]==0.5i,(mate.pref.1+1)/2,1)),  #male attractiveness for P=0        female of sp.1
                             ifelse(x[,"T"]==0,(mate.pref.1+1)/2,ifelse(x[,"T"]==0.5i,(mate.pref.1+3)/4,1)),                      #for P=0.5i     female of sp.1
                             rep(1,nrow(x)),                                                            #for all neutral female
                             ifelse(x[,"T"]==1,(mate.pref.1+1)/2,ifelse(x[,"T"]==0.5+0.5i,(mate.pref.1+3)/4,1)),                  #for P=0.5+0.5i female of sp.1
                             ifelse(x[,"T"]==1, mate.pref.1,     ifelse(x[,"T"]==0.5+0.5i,(mate.pref.1+1)/2,1)),                  #for P=1        female of sp.1
                             ifelse(x[,"T"]==0, mate.pref.2,     ifelse(x[,"T"]==0.5i,(mate.pref.2+1)/2,1)),                      #for P=0        female of sp.2
                             ifelse(x[,"T"]==0,(mate.pref.2+1)/2,ifelse(x[,"T"]==0.5i,(mate.pref.2+3)/4,1)),                      #for P=0.5i     female of sp.2
                             ifelse(x[,"T"]==1,(mate.pref.2+1)/2,ifelse(x[,"T"]==0.5+0.5i,(mate.pref.2+3)/4,1)),                  #for P=0.5+0.5i female of sp.2
                             ifelse(x[,"T"]==1 ,mate.pref.2,     ifelse(x[,"T"]==0.5+0.5i,(mate.pref.2+1)/2,1)))                  #for P=1        female of sp.2
        potential.dad.id <- cbind(sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,1]), #for P=0        female of sp.1
                                  sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,2]), #for P=0.5i     female of sp.1
                                  sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,3]), #for all neutral female
                                  sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,4]), #for P=0.5+0.5i female of sp.1
                                  sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,5]), #for P=1        female of sp.1
                                  sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,6]), #for P=0        female of sp.2
                                  sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,7]), #for P=0.5i     female of sp.2
                                  sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,8]), #for P=0.5+0.5i female of sp.2
                                  sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,9])) #for P=1        female of sp.2
        dad.id <- potential.dad.id[,1]*ifelse(mom.genotype[,"sp"]==0&mom.genotype[,"P"]==0,       1,0)+
                  potential.dad.id[,2]*ifelse(mom.genotype[,"sp"]==0&mom.genotype[,"P"]==0.5i,    1,0)+
                  potential.dad.id[,3]*ifelse(mom.genotype[,"P"]==0.5|mom.genotype[,"P"]==1i,     1,0)+
                  potential.dad.id[,4]*ifelse(mom.genotype[,"sp"]==0&mom.genotype[,"P"]==0.5+0.5i,1,0)+
                  potential.dad.id[,5]*ifelse(mom.genotype[,"sp"]==0&mom.genotype[,"P"]==1,       1,0)+
                  potential.dad.id[,6]*ifelse(mom.genotype[,"sp"]==1&mom.genotype[,"P"]==0,       1,0)+
                  potential.dad.id[,7]*ifelse(mom.genotype[,"sp"]==1&mom.genotype[,"P"]==0.5i,    1,0)+
                  potential.dad.id[,8]*ifelse(mom.genotype[,"sp"]==1&mom.genotype[,"P"]==0.5+0.5i,1,0)+
                  potential.dad.id[,9]*ifelse(mom.genotype[,"sp"]==1&mom.genotype[,"P"]==1,       1,0)
        dad.genotype <- x[dad.id,1:6]
        egg <- as.data.frame(cbind(0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,1])), #Re() only yet
                                   0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,2])),
                                   0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,3])),
                                   0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,4])),
                                   0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,5])),
                                   0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,6]))))
        sperm <- as.data.frame(cbind(0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,1])), #Re() only yet
                                     0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,2])),
                                     0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,3])),
                                     0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,4])),
                                     0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,5])),
                                     0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,6]))))
        colnames(egg) <- c("M","F","T","P","D","E"); colnames(sperm) <- c("M","F","T","P","D","E")
        egg$M <- egg$M+ifelse(egg$M==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,1])/(1-Re(mom.genotype[,1])))) #add Im() component
        egg$F <- egg$F+ifelse(egg$F==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,2])/(1-Re(mom.genotype[,2]))))  
        egg$T <- egg$T+ifelse(egg$T==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,3])/(1-Re(mom.genotype[,3]))))
        egg$P <- egg$P+ifelse(egg$P==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,4])/(1-Re(mom.genotype[,4]))))
        sperm$M <- sperm$M+ifelse(sperm$M==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,1])/(1-Re(dad.genotype[,1])))) #add Im() component
        sperm$F <- sperm$F+ifelse(sperm$F==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,2])/(1-Re(dad.genotype[,2]))))  
        sperm$T <- sperm$T+ifelse(sperm$T==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,3])/(1-Re(dad.genotype[,3]))))
        sperm$P <- sperm$P+ifelse(sperm$P==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,4])/(1-Re(dad.genotype[,4]))))
        offspring <- egg+sperm
      } else {                           #mating algorithm for models where mate choice does not evolve
         pref.matrix <- cbind(ifelse(x[,"sp"]==0,mate.pref.1,1),ifelse(x[,"sp"]==1,mate.pref.2,1))     #male attractiveness for females of sp.1 and sp.2, respectively
         potential.dad.id <- cbind(sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,1]), #for sp.1 female
                                   sample(1:nrow(x),nrow(mom.genotype),replace=TRUE,prob=pref.matrix[,2])) #for sp.2 female
         dad.id <- potential.dad.id[,1]*(1-mom.genotype[,"sp"])+potential.dad.id[,2]*mom.genotype[,"sp"]
         dad.genotype <- x[dad.id,1:6]
         egg <- as.data.frame(cbind(0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,1])), #Re() only yet
                                    0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,2])),
                                    0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,3])),
                                    0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,4])),
                                    0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,5])),
                                    0.5*rbinom(nrow(mom.genotype),1,Re(mom.genotype[,6]))))
         sperm <- as.data.frame(cbind(0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,1])), #Re() only yet
                                      0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,2])),
                                      0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,3])),
                                      0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,4])),
                                      0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,5])),
                                      0.5*rbinom(nrow(dad.genotype),1,Re(dad.genotype[,6]))))
         colnames(egg) <- c("M","F","T","P","D","E"); colnames(sperm) <- c("M","F","T","P","D","E")
         egg$M <- egg$M+ifelse(egg$M==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,1])/(1-Re(mom.genotype[,1])))) #add Im() component
         egg$F <- egg$F+ifelse(egg$F==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,2])/(1-Re(mom.genotype[,2]))))  
         egg$T <- egg$T+ifelse(egg$T==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,3])/(1-Re(mom.genotype[,3]))))
         egg$P <- egg$P+ifelse(egg$P==0.5,0,0.5i*rbinom(nrow(mom.genotype),1,Im(mom.genotype[,4])/(1-Re(mom.genotype[,4]))))
         sperm$M <- sperm$M+ifelse(sperm$M==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,1])/(1-Re(dad.genotype[,1])))) #add Im() component
         sperm$F <- sperm$F+ifelse(sperm$F==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,2])/(1-Re(dad.genotype[,2]))))  
         sperm$T <- sperm$T+ifelse(sperm$T==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,3])/(1-Re(dad.genotype[,3]))))
         sperm$P <- sperm$P+ifelse(sperm$P==0.5,0,0.5i*rbinom(nrow(dad.genotype),1,Im(dad.genotype[,4])/(1-Re(dad.genotype[,4]))))
         offspring <- egg+sperm
      }
  colnames(offspring) <- c("M","F","T","P","D","E")
  offspring
}

#x: data frame of adults before dispersal or mating; other param:sex.limitation, MAD, signal.surv.cost, mate.choice.surv.cost, habitat.choice.surv.cost
sin.of.preference <- function(x) {
  if (MAD==1) {
     if (sex.limitation==0) {
       rbinom(nrow(x),1,1-habitat.choice.surv.cost*ifelse(x[,"M"]==0.5,0,1-Im(x[,"M"]))-
                           signal.surv.cost*x[,"sex"]*ifelse(x[,"T"]==0.5,0,1-Im(x[,"T"]))-
                           mate.choice.surv.cost*(1-x[,"sex"])*ifelse(x[,"P"]==0.5,0,1-Im(x[,"P"])))
     } else {
       rbinom(nrow(x),1,1-habitat.choice.surv.cost*x[,"sex"]*ifelse(x[,"M"]==0.5,0,1-Im(x[,"M"]))-
                          habitat.choice.surv.cost*(1-x[,"sex"])*ifelse(x[,"F"]==0.5,0,1-Im(x[,"F"]))-
                          signal.surv.cost*x[,"sex"]*ifelse(x[,"T"]==0.5,0,1-Im(x[,"T"]))-
                          mate.choice.surv.cost*(1-x[,"sex"])*ifelse(x[,"P"]==0.5,0,1-Im(x[,"P"])))
     }
  } else {
    if (sex.limitation==0) {
      rbinom(nrow(x),1,1-habitat.choice.surv.cost*(1-x[,"sex"])*ifelse(x[,"M"]==0.5,0,1-Im(x[,"M"]))-
                          signal.surv.cost*x[,"sex"]*ifelse(x[,"T"]==0.5,0,1-Im(x[,"T"]))-
                          mate.choice.surv.cost*(1-x[,"sex"])*ifelse(x[,"P"]==0.5,0,1-Im(x[,"P"])))
    } else {
      rbinom(nrow(x),1,1-habitat.choice.surv.cost*(1-x[,"sex"])*ifelse(x[,"F"]==0.5,0,1-Im(x[,"F"]))-
               signal.surv.cost*x[,"sex"]*ifelse(x[,"T"]==0.5,0,1-Im(x[,"T"]))-
               mate.choice.surv.cost*(1-x[,"sex"])*ifelse(x[,"P"]==0.5,0,1-Im(x[,"P"])))
    }
  }
}


##Switches ############ the parameters below are for cost-free joint evolution model
#######################

#sex limited gene expression
sex.limitation <- 0 #0:M is expressed in both sex and F is not, 1:sex limited(M&F)

#mating before or after dispersal?
MAD <- 1 #0:mating before dispersal, 1:mating after dispersal

#allow habitat preference to evolve?
evolve.habitat.pref <- 1 #0:habitat pref do not evolve (habitat.pref=1), 1:habitat pref can evolve

#allow species recognition to evolve?
evolve.sp.recognition <- 1 #0:sp. recognition do not evolve (set "mate.pref.1","mate.pref.2"), 1:sp. recognition can evolve
#CAUTION!: even when sp.recognition does not evolve, assortative mating occurs unless mate.pref.x=1

#Run the program #########
##########################

t <- proc.time() #record time
iteration <- 100 #iteration per param. combination
pref.range <- c(1,2,4,8,16,32,64,128,256,512,1024)
s.hybrid.range <- c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
pref.range.length <- length(pref.range)
s.hybrid.range.length <- length(s.hybrid.range)

result <- NULL
param  <- NULL
cl <- makeCluster(detectCores())
registerDoParallel(cl)

for (j in 1:pref.range.length) {
  for (k in 1:s.hybrid.range.length) {
     param <- rbind(param, cbind(rep(pref.range[j],iteration), rep(s.hybrid.range[k],iteration)))
  }
}

result <- foreach(l = 1:nrow(param),.combine="rbind",.inorder=FALSE,.packages = "dplyr") %dopar% {
      
      ###initial conditions #########
      ###############################
      n01 <- 100 #initial population size of sp.1 (m + f)
      n02 <- 100 #initial population size of sp.2
      mA.1 <- 1/3 #initial allele freq. for male habitat pref. for B of sp.1
      mA.2 <- 1/3 #                                                     sp.2
      mB.1 <- 1/3 #initial allele freq. for male habitat pref. for B of sp.1
      mB.2 <- 1/3 #                                                     sp.2
      mN.1 <- 1-mA.1-mB.1
      mN.2 <- 1-mA.2-mB.2
      fA.1 <- 1/3 #initial allele freq. for female habitat pref. for B of sp.1
      fA.2 <- 1/3 #                                                       sp.2
      fB.1 <- 1/3
      fB.2 <- 1/3
      fN.1 <- 1-fA.1-fB.1
      fN.2 <- 1-fA.2-fB.2
      t0.1 <- 1/3 #initial allele freq. for male signaling trait for sp.1
      t0.2 <- 1/3 #                                                  sp.2
      t1.1 <- 1/3
      t1.2 <- 1/3
      tN.1 <- 1-t0.1-t1.1
      tN.2 <- 1-t0.2-t1.2
      p0.1 <- 1/3 #initial allele freq. for female mate preference for sp.1
      p0.2 <- 1/3 #                                                    sp.2
      p1.1 <- 1/3
      p1.2 <- 1/3
      pN.1 <- 1-p0.1-p1.1
      pN.2 <- 1-p0.2-p1.2
      d1 <- 0   #initial allele freq. for epistatic locus1 for sp.1
      d2 <- 1   #                                              sp.2
      e1 <- 0   #initial allele freq. for epistatic locus2 for sp.1
      e2 <- 1   #                                              sp.2
      
      #### define parameter ###############
      #####################################
      habitat.pref <- if(evolve.habitat.pref==0) {1} else {param[l,1]} #the extent of habitat preference(change only the latter)
      mate.pref.1 <- param[l,1] #the extent of mate preference for sp.1
      mate.pref.2 <- param[l,1]
      signal.surv.cost <- 0 #survival cost of male signal
      mate.choice.surv.cost <- 0 #survival cost of female mate choice
      habitat.choice.surv.cost <- 0 #survival cost of habitat pref (which is irrelevant to males of MBD model)
      s.hybrid <- param[l,2]
      R.inA <- 3.5
      R.inB <- 3.5
      K.inA <- 1000 #carryng capacity measured in female no.
      K.inB <- 1000
      
      G1 <- NULL; G2 <- NULL #the model begins with the invasions by adults to empty habitat(e.g. range expansion)
      G1 <- cbind(G1, 0.5*rbinom(n01,1,mB.1)+0.5*rbinom(n01,1,mB.1),  #M (only mB)
                      0.5*rbinom(n01,1,fB.1)+0.5*rbinom(n01,1,fB.1),  #F (only fB)
                      0.5*rbinom(n01,1,t1.1)+0.5*rbinom(n01,1,t1.1),  #T (only t1)
                      0.5*rbinom(n01,1,p1.1)+0.5*rbinom(n01,1,p1.1),  #P (only p1)
                      0.5*rbinom(n01,1,d1)+0.5*rbinom(n01,1,d1),      #D
                      0.5*rbinom(n01,1,e1)+0.5*rbinom(n01,1,e1),      #E
                      rep(0,n01),                                     #natal habitat
                      floor(runif(n01,0,2)))                          #sex 1:m, 0:f
      G2 <- cbind(G2, 0.5*rbinom(n02,1,mB.2)+0.5*rbinom(n02,1,mB.2),
                      0.5*rbinom(n02,1,fB.2)+0.5*rbinom(n02,1,fB.2),
                      0.5*rbinom(n02,1,t1.2)+0.5*rbinom(n02,1,t1.2),
                      0.5*rbinom(n02,1,p1.2)+0.5*rbinom(n02,1,p1.2),
                      0.5*rbinom(n02,1,d2)+0.5*rbinom(n02,1,d2),
                      0.5*rbinom(n02,1,e2)+0.5*rbinom(n02,1,e2),
                      rep(1,n02),
                      floor(runif(n02,0,2)))
      G1 <- as.data.frame(G1); G2 <- as.data.frame(G2)
      colnames(G1) <- c("M","F","T","P","D","E","nat.hab","sex")
      colnames(G2) <- c("M","F","T","P","D","E","nat.hab","sex")
      G1$sp <- sp.identifier(G1[,"D"]+G1[,"E"])
      G2$sp <- sp.identifier(G2[,"D"]+G2[,"E"])
      
      G1$M <- G1$M+ifelse(G1$M==1,0,0.5i*rbinom(n01,1,mN.1/(mN.1+mA.1)))+ifelse(G1$M!=0,0,0.5i*rbinom(n01,1,mN.1/(mN.1+mA.1))) #add mN and mA
      G1$F <- G1$F+ifelse(G1$F==1,0,0.5i*rbinom(n01,1,fN.1/(fN.1+fA.1)))+ifelse(G1$F!=0,0,0.5i*rbinom(n01,1,fN.1/(fN.1+fA.1))) #add fN and fB
      G1$T <- G1$T+ifelse(G1$T==1,0,0.5i*rbinom(n01,1,tN.1/(tN.1+t0.1)))+ifelse(G1$T!=0,0,0.5i*rbinom(n01,1,tN.1/(tN.1+t0.1))) #add tN and t0
      G1$P <- G1$P+ifelse(G1$P==1,0,0.5i*rbinom(n01,1,pN.1/(pN.1+p0.1)))+ifelse(G1$P!=0,0,0.5i*rbinom(n01,1,pN.1/(pN.1+p0.1))) #add pN and t0
      G2$M <- G2$M+ifelse(G2$M==1,0,0.5i*rbinom(n02,1,mN.2/(mN.2+mA.2)))+ifelse(G2$M!=0,0,0.5i*rbinom(n02,1,mN.2/(mN.2+mA.2)))
      G2$F <- G2$F+ifelse(G2$F==1,0,0.5i*rbinom(n02,1,fN.2/(fN.2+fA.2)))+ifelse(G2$F!=0,0,0.5i*rbinom(n02,1,fN.2/(fN.2+fA.2)))
      G2$T <- G2$T+ifelse(G2$T==1,0,0.5i*rbinom(n02,1,tN.2/(tN.2+t0.2)))+ifelse(G2$T!=0,0,0.5i*rbinom(n02,1,tN.2/(tN.2+t0.2)))
      G2$P <- G2$P+ifelse(G2$P==1,0,0.5i*rbinom(n02,1,pN.2/(pN.2+p0.2)))+ifelse(G2$P!=0,0,0.5i*rbinom(n02,1,pN.2/(pN.2+p0.2)))
      
      G1$sec.hab <- habitat.determiner(G1[,"M"],G1[,"F"],G1[,"sex"]) #secondary habitat, to which they disperse
      G2$sec.hab <- habitat.determiner(G2[,"M"],G2[,"F"],G2[,"sex"]) #secondary habitat, to which they disperse
      
      if (nrow(G1)==0) {stats1 <- c(NA,NA,NA,NA,0)} else {stats1 <- c(apply(G1[,1:4],2,mean), nrow(G1)) } #recording statistics
      if (nrow(G2)==0) {stats2 <- c(NA,NA,NA,NA,0)} else {stats2 <- c(apply(G2[,1:4],2,mean), nrow(G2)) } #freq.M,F,T,P & sp.abundance
      stats3 <- c(NA,NA,NA,NA,0)
      
      G.all <- rbind(G1,G2)
      G.all$fitness <- sin.of.preference(G.all) #fitness cost due to signal/mate pref/habitat pref
      G.all <- subset(G.all,fitness==1)[,1:10]
      
      generation <- 1
      extinct <- 0
      fix.M <- 0
      fix.P <- 0
      fix.P.or.M.prev <- 0
      fix.P.or.M.now <- 0
      partial.div <- NULL
      
      while (generation<=20000 && extinct==0 && (fix.M == 0 || fix.P == 0)) {
        generation <- generation + 1
        fix.P.or.M.prev <- fix.P.or.M.now
        
        if (MAD==0) {                                  #mating in mating-before-dispersal (non.MAD or MBD) model
          m.inA <- subset(G.all,nat.hab==0&sex==1)
          f.inA <- subset(G.all,nat.hab==0&sex==0)
          m.inB <- subset(G.all,nat.hab==1&sex==1)
          f.inB <- subset(G.all,nat.hab==1&sex==0)
          offspring.of.mom.from.A <- non.MAD.mating.matcher(m.inA,f.inA) #offspring production by females who mated in habitat A
          offspring.of.mom.from.B <- non.MAD.mating.matcher(m.inB,f.inB) #offspring production by females who mated in habitat B
          offspring <- rbind(offspring.of.mom.from.A,offspring.of.mom.from.B)
          offspring$fitness <- fitness.assigner(offspring[,"D"],offspring[,"E"]) #first round of mortality due to genetic incompatibility
          juvenile.inA <- subset(offspring,nat.hab==0&fitness==1)[,1:7] # this data.frame includes M, F, T, P, D, E, nat.hab(i.e. sec.hab for mothers)
          juvenile.inB <- subset(offspring,nat.hab==1&fitness==1)[,1:7]
        } else {                                       #mating in mating-after-dispersal (MAD) model
          m.inA <- subset(G.all,sec.hab==0&sex==1)
          f.inA <- subset(G.all,sec.hab==0&sex==0)
          m.inB <- subset(G.all,sec.hab==1&sex==1)
          f.inB <- subset(G.all,sec.hab==1&sex==0)
          offspring.inA <- MAD.mating.matcher(m.inA,f.inA,0) #offspring production by adults in A
          offspring.inB <- MAD.mating.matcher(m.inB,f.inB,1) #offspring production by adults in B
          offspring.inA$fitness <- fitness.assigner(offspring.inA[,"D"],offspring.inA[,"E"]) #first round of mortality due to genetic incompatibility
          offspring.inB$fitness <- fitness.assigner(offspring.inB[,"D"],offspring.inB[,"E"]) #first round of mortality due to genetic incompatibility
          juvenile.inA <- subset(offspring.inA,fitness==1)[,1:6]
          juvenile.inB <- subset(offspring.inB,fitness==1)[,1:6]
          juvenile.inA$nat.hab <- 0; juvenile.inB$nat.hab <- 1
        }
        
        if.inA <- nrow(juvenile.inA)/(2*R.inA) #imaginary reproducing female abundance in A for density regulation
        survivors.id.inA <- sample(1:nrow(juvenile.inA),2*R.inA*K.inA*if.inA/(K.inA+(R.inA-1)*if.inA),replace=FALSE) 
        G.inA <- juvenile.inA[survivors.id.inA,] #density regulation in A
        
        if.inB <- nrow(juvenile.inB)/(2*R.inB) #imaginary reproducing female abundance in B for density regulation
        survivors.id.inB <- sample(1:nrow(juvenile.inB),2*R.inB*K.inB*if.inB/(K.inB+(R.inB-1)*if.inB),replace=FALSE) 
        G.inB <- juvenile.inB[survivors.id.inB,] #density regulation in A
        
        G.all <- rbind(G.inA,G.inB)
        G.all$sex <- floor(runif(nrow(G.all),0,2))
        G.all$sp <- sp.identifier(G.all[,"D"]+G.all[,"E"])
        G.all$sec.hab <- habitat.determiner(G.all[,"M"],G.all[,"F"],G.all[,"sex"]) #newly matured adults
        
        G1.st <- subset(G.all,D+E==0) #taking stats
        G2.st <- subset(G.all,D*E==1)
        Ghyb.st <- subset(G.all,D+E!=0&D*E!=1)
        
        if (nrow(G1.st)==0) {stats1 <- rbind(stats1,c(NA,NA,NA,NA,0))} else {stats1 <- rbind(stats1,c(apply(G1.st[,1:4],2,mean), nrow(G1.st))) }
        if (nrow(G2.st)==0) {stats2 <- rbind(stats2,c(NA,NA,NA,NA,0))} else {stats2 <- rbind(stats2,c(apply(G2.st[,1:4],2,mean), nrow(G2.st))) }
        if (nrow(Ghyb.st)==0) {stats3 <- rbind(stats3,c(NA,NA,NA,NA,0))} else {stats3 <- rbind(stats3,c(apply(Ghyb.st[,1:4],2,mean), nrow(Ghyb.st))) }
        
        G.all$fitness <- sin.of.preference(G.all)
        G.all <- subset(G.all,fitness==1)[,1:10]
        
        extinct <- if (nrow(G1.st)==0||nrow(G2.st)==0) {1} else {0}
        fix.M <- if (fix.M == 1) {1} else {
                  if (generation >= 100 && 
                   (mean(Re(stats1[(nrow(stats1)-99):nrow(stats1),"M"]) + Im(stats1[(nrow(stats1)-99):nrow(stats1),"M"]), na.rm = T) <= 0.05 ||
                    mean(Re(stats1[(nrow(stats1)-99):nrow(stats1),"M"]), na.rm = T) >= 0.95 ||
                    mean(Im(stats1[(nrow(stats1)-99):nrow(stats1),"M"]), na.rm = T) >= 0.95 ) &&
                   (mean(Re(stats2[(nrow(stats2)-99):nrow(stats2),"M"]) + Im(stats2[(nrow(stats2)-99):nrow(stats2),"M"]), na.rm = T) <= 0.05 ||
                    mean(Re(stats2[(nrow(stats2)-99):nrow(stats2),"M"]), na.rm = T) >= 0.95 ||
                    mean(Im(stats2[(nrow(stats2)-99):nrow(stats2),"M"]), na.rm = T) >= 0.95)) {1} else {0} }
        fix.P <- if (fix.P == 1) {1} else {
                  if (generation >= 100 && 
                   (mean(Re(stats1[(nrow(stats1)-99):nrow(stats1),"P"]) + Im(stats1[(nrow(stats1)-99):nrow(stats1),"P"]), na.rm = T) <= 0.05 ||
                    mean(Re(stats1[(nrow(stats1)-99):nrow(stats1),"P"]), na.rm = T) >= 0.95 ||
                    mean(Im(stats1[(nrow(stats1)-99):nrow(stats1),"P"]), na.rm = T) >= 0.95 ) &&
                   (mean(Re(stats2[(nrow(stats2)-99):nrow(stats2),"P"]) + Im(stats2[(nrow(stats2)-99):nrow(stats2),"P"]), na.rm = T) <= 0.05 ||
                    mean(Re(stats2[(nrow(stats2)-99):nrow(stats2),"P"]), na.rm = T) >= 0.95 ||
                    mean(Im(stats2[(nrow(stats2)-99):nrow(stats2),"P"]), na.rm = T) >= 0.95)) {1} else {0} }
        if (fix.M == 1 || fix.P == 1) {fix.P.or.M.now <- 1} else {}
        if (fix.P.or.M.now - fix.P.or.M.prev == 1) {
           sp1.mean.0 <- stats1[nrow(stats1),]
           sp2.mean.0 <- stats2[nrow(stats2),]
           hyb.mean.0 <- stats3[nrow(stats3),]
           partial.div <- c(l,0,fix.M,fix.P,extinct,generation,habitat.pref,s.hybrid,sp1.mean.0,sp2.mean.0,hyb.mean.0) #taking the stats
        } else {}
      }
   
      #by removing #, you can plot genotypic changes and demographic dynamics in a single run
      #par(mfrow=c(3,1)) #habitat preference
      #plot(Re(stats1[,1])~c(1:nrow(stats1)),type="l",lwd=2,ylim=c(0,1)) #thiner lines for neutral alleles
      #lines(Im(stats1[,1])~c(1:nrow(stats1)))
      #lines(Re(stats1[,2])~c(1:nrow(stats1)),lty=3,lwd=2)
      #lines(Im(stats1[,2])~c(1:nrow(stats1)),lty=3)
      #lines(Re(stats2[,1])~c(1:nrow(stats2)),col=2,lwd=2)
      #lines(Im(stats2[,1])~c(1:nrow(stats2)),col=2)
      #lines(Re(stats2[,2])~c(1:nrow(stats2)),lty=3,col=2,lwd=2)
      #lines(Im(stats2[,2])~c(1:nrow(stats2)),lty=3,col=2)
      
      #plot(Re(stats1[,3])~c(1:nrow(stats1)),type="l",ylim=c(0,1),lwd=2) #mate preference
      #lines(Im(stats1[,3])~c(1:nrow(stats1)))
      #lines(Re(stats1[,4])~c(1:nrow(stats1)),lty=3,lwd=2)
      #lines(Im(stats1[,4])~c(1:nrow(stats1)),lty=3)
      #lines(Re(stats2[,3])~c(1:nrow(stats2)),col=2,lwd=2)
      #lines(Im(stats2[,3])~c(1:nrow(stats2)),col=2)
      #lines(Re(stats2[,4])~c(1:nrow(stats2)),lty=3,col=2,lwd=2)
      #lines(Im(stats2[,4])~c(1:nrow(stats2)),lty=3,col=2)
      
      #plot(Re(stats1[,5])~c(1:nrow(stats1)),type="l",ylim=c(0,4000)) #abundance
      #lines(Re(stats2[,5])~c(1:nrow(stats2)),col=2)
      #lines(Re(stats3[,5])~c(1:nrow(stats3)),col=4)
      
      sp1.mean <- stats1[nrow(stats1),]
      sp2.mean <- stats2[nrow(stats2),]
      hyb.mean <- stats3[nrow(stats3),]
      rbind(partial.div,c(l,1,fix.M,fix.P,extinct,generation,habitat.pref,s.hybrid,sp1.mean,sp2.mean,hyb.mean)) #taking the stats
}

stopCluster(cl)
result <- as.data.frame(result)

colnames(result) <- c("run.ID","end.of.run","fix.M","fix.P","extinct","generation","pref","s.hybrid","sp1.M","sp1.F","sp1.T","sp1.P","sp1.ab",
                                                                                                     "sp2.M","sp2.F","sp2.T","sp2.P","sp2.ab",
                                                                                                     "hyb.M","hyb.F","hyb.T","hyb.P","hyb.ab")
proc.time()-t #time count
rownames(result) <- 1:nrow(result)
joint.no.cost.two.termination <- result

write.table(joint.no.cost.two.termination,"joint.no.cost.two.termination.txt",quote=F) #export data

# end.of.run = 1 gives the end result of a simulation run
# end.of.run = 0 gives a middle state of a run, where one isolation mechanism diverged but the other did not.