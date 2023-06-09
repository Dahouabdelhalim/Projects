#######################################
#######################################
buildped<-function(founders=50,fert=2.0,ngen=2){
  nf<- founders
  sanc<-nf #sum of ancestors
  ni<-nf #number of inds each generation
  parents<-1:nf
  #ncoup<-nf %/% 2 #monogamy
  ped<-data.frame(ind=1:nf,mum=rep(NA,nf),dad=rep(NA,nf))
  
  draw.parents<-function(parents,fert,sanc){
    ncoup<-length(parents) %/% 2
    nbabs<-rpois(ncoup,fert)
    allbabs<-sum(nbabs)
    thisgen.ped<-data.frame(ind=integer(allbabs),
                            mum=integer(allbabs),
                            dad=integer(allbabs))
    for (ic in 1:ncoup){
      if(nbabs[ic]>0){  
        # no selfing, monogamy, but can reuse parents
        par.drawn<-sample(parents,size=2,replace=FALSE) 
        #  parents<-parents[-par.drawn]
        nbabsc<-c(0,cumsum(nbabs))
        x<-(nbabsc[ic]+1):nbabsc[ic+1]
        thisgen.ped[x,]<-cbind(sanc+x,rep(par.drawn[1],nbabs[ic]),
                               rep(par.drawn[2],nbabs[ic]))
      }
    }
    return(thisgen.ped)  
  }
  
  for (igen in 1:ngen){
    thisgen.ped<-draw.parents(parents,fert,sanc)
    ped<-rbind(ped,thisgen.ped)
    ni<-c(ni,dim(thisgen.ped)[1])
    sanc<-sum(ni)  
    parents<-(sum(ni[1:igen])+1):sanc
  }
  return(ped)
}
#######################################
#######################################
buildped.rm<-function(founders=50,fert=2.0,ngen=2){
  nf<- founders
  sanc<-nf #sum of ancestors
  ni<-nf #number of inds each generation
  parents<-1:nf
  ped<-data.frame(ind=1:nf,mum=rep(NA,nf),dad=rep(NA,nf))
  
  draw.parents<-function(parents,fert,sanc){
    ncoup<-length(parents) %/% 2
    nbabs<-rpois(ncoup,fert)
    allbabs<-sum(nbabs)
    thisgen.ped<-data.frame(ind=integer(allbabs),
                            mum=integer(allbabs),
                            dad=integer(allbabs))
    for (ic in 1:allbabs){
        # promiscuity, no selfing, 1 offs per mating
        par.drawn<-sample(parents,size=2,replace=FALSE) 
        #  parents<-parents[-par.drawn]
        #nbabsc<-c(0,cumsum(nbabs))
        thisgen.ped[ic,]<-cbind(sanc+ic,par.drawn[1],par.drawn[2])
      }
   
    return(thisgen.ped)  
  }
  
  for (igen in 1:ngen){
    thisgen.ped<-draw.parents(parents,fert,sanc)
    ped<-rbind(ped,thisgen.ped)
    ni<-c(ni,dim(thisgen.ped)[1])
    sanc<-sum(ni)  
    parents<-(sum(ni[1:igen])+1):sanc
  }
  return(ped)
}
##########################################################
draw.offsprings.ped.SNPs<-function(ped=ped,founders.genotypes=dat,ndigits=2){
  #assumes a pedigree (3 columns, ind, dam,sire)
  #and a data frame dat with founders genotypes
  #encoded as allelic dosage (0,1,2)
  #dam and sire of founders are entered as NA  
  
  nfounders<-dim(founders.genotypes)[1]
  nloc<-dim(founders.genotypes)[2]
  noffs<-dim(ped)[1]-nfounders
  genos<-data.frame(matrix(numeric(nloc*(nfounders+noffs)),ncol=nloc))
  names(genos)<-names(founders.genotypes)
  genos[1:nfounders,]<-founders.genotypes
  seqoffs<-(nfounders+1):(nfounders+noffs)
  mumid<-0
  dadid<-0
  for (io in seqoffs){
    if (ped[io,2]!=mumid | ped[io,3]!=dadid){
      mumid<-ped[io,2]
      dadid<-ped[io,3]           
      mum<-genos[mumid,]
      dad<-genos[dadid,]
      #      dat<-rbind(mum,dad)
      #only draw heterozygous positions. Not sure it speeds things up.   
      het.mum<-which(mum==1)
      lhm<-length(het.mum)
      het.dad<-which(dad==1)
      lhd<-length(het.dad)
    }
    tmp1<-rbinom(lhm,1,.5)
    tmp2<-rbinom(lhd,1,.5)
    gam.mum<-mum %/%2
    gam.mum[het.mum]<-tmp1
    gam.dad<-dad %/% 2
    gam.dad[het.dad]<-tmp2
    genos[io,]<-gam.mum+gam.dad 
  }
  genos
}
######################################
beta.coan.SNPs<-function(dat){
  #dat is a data frame with individuals in rows and allelic dosage for each locus in colums  
  #uses matching proba -same equation as for population i.e. Mij=[xiXj+(2-xi)(2-xj)]/4
  dat<-as.matrix(dat)
  nl<-dim(dat)[2]
  Mij<-(tcrossprod(dat)+tcrossprod(2-dat))/4
  diag(Mij)<-NA
  Mb<-mean(Mij,na.rm=T)
  (Mij-Mb)/(nl-Mb)
}
######################################
GCTA.SNPs<-function(dat){
  #weighted estimate of GCTA i.e sum_nl[(xil-2pl)(xjl-2pl)]/4 sum_nl[pl(1-pl)]
  dat<-as.matrix(dat)
  sfs<-colSums(dat)
  ni<-dim(dat)[1]
  nl<-dim(dat)[2]
  datc<-sweep(dat,2,sfs/ni,FUN="-")
  num<-tcrossprod(datc)
  diag(num)<-NA
  den<-sum(sfs/2/ni*(1-sfs/2/ni))*4
  num/den
}
GCTAu.SNPs<-function(dat){
  #unweighted estimate of GCTA i.e sum_nl[(xil-2pl)(xjl-2pl)]/4 sum_nl[pl(1-pl)]
  dat<-as.matrix(dat)
  sfs<-colMeans(dat,na.rm=TRUE)
  pol<-which(sfs>0 & sfs<2)
  dat<-dat[,pol]
  ni<-dim(dat)[1]
  nl<-dim(dat)[2]
  p<-sfs[pol]/2
  w<-(p*(1-p))^.5*2
 # ni<-dim(dat)[1]
 # nl<-dim(dat)[2]
  datc<-sweep(dat,2,2*p,FUN="-")
  datc[is.na(datc)]<-0.0
  datcr<-sweep(datc,2,w,FUN="/")
  res<-tcrossprod(datcr)/nl
  diag(res)<-NA
  res
  }
