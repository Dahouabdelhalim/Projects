
PasseyInverse<-function(Length1,dMeas1,depth1,finit,la,lm,maxlength,minlength,mindepth,df1,nsolxns){
  MEST<-list()
  DPE<-c()
  S<-list()
  for(n in 1:nsolxns){
    print(n)
    # SECTION 1: USER-ENTERED DATA
    finit<-finit
    la<-la1
    lm<-lm1
    openindx<-1
    avelength<-round(mean(Length1))
    maxlength<-maxlength1
    minlength <- minlength1
    mindepth <- mindepth1
    r1 <- .05
    r2 <- 1
    r3 <- 5
    df <- df1
    Length <-Length1
    depth<-depth1
    dMeas <- dMeas1
    maxratio <- max(dMeas)
    minratio <- min(dMeas)
    stdev <- sd(dMeas)
    numsam<-length(dMeas)
    # SECTION 2: ADDITION OF RANDOM ERROR
    rlength <- r2*rnorm(numsam)	
    Length <- round(rlength) + Length
    for(a in 1:numsam){
      if(Length[a]>maxlength){
        Length[a]<-maxlength
      }else{
        if(Length[a]<minlength){
          Length[a]<-minlength
        }else{ 
          Length[a]<-Length[a]
        }
      }
    }
    rdepth <- r3*rnorm(numsam)
    rdepth<-round(rdepth)
    depth2 <- round(rdepth) + depth
    for(z in 1:numsam){
      if(depth2[z]>la){
        depth2[z]<-la
      }
      if(depth2[a]<1.1){
        depth2[z]<-mindepth
      }
    }
    depth<-depth2
    # Copied from Matlab code associated with Passey et al: Determines the number of m's distal and proximal to those that directly correspond with d's
    #numbefore reflects the m's that are sampled into at the beginning of the profile because of sampling depth
    #numafter reflects the m's that contribute to the isotope values of the
    #final samples in open-ended cases.
    numbefore<-ceiling(la/avelength)
    numafter<-ceiling((lm-openindx)/avelength) +1
    lengthbefore<-avelength*rep(1,numbefore)
    lengthafter<-avelength*rep(1,numafter)
    fmat<- 1 - finit
    numcol<-numbefore + numsam + numafter
    Length2<-c(lengthbefore,Length,lengthafter,lm)
    depth<-c(depth,lengthbefore)
    # SECTION 3: CONSTRUCTION OF AVERAGING MATRIX 
    # SECTION 3.1 BINARY MATRIX 
    B<-matrix(nrow=Length2[1],ncol=numcol,0)
    B[,1]<-1
    for(m in 2:numcol){
      F<-matrix(nrow=Length2[m],ncol=numcol,0)
      F[,m]<-1
      B1<-rbind(B,F)
      B<-B1
    }
    F<-matrix(nrow=(((sum(Length2)-openindx)+1):sum(Length2))-nrow(B),ncol=numcol,0)# I think this is doing the right thing
    B<-rbind(B,F)
    B[((sum(Length2)-openindx)+1):sum(Length2),]<-c(0)
    B<-B[,1:(numcol-1)]
    # SECTION 3.2 MATURATION AVERAGE OF BINARY MATRIX 
    o<-1
    AB <- c(finit*(colMeans(B[o:(o+lm-1),])) + fmat*(colMeans(B[o:(o+lm-1),])))
    AB <- AB/(sum(B[o,]))
    for(o in 1:(sum(Length2)-lm)){
      p<-c(finit*(B[o,]) + fmat*(colMeans(B[o:(o+lm-1),])))
      if(sum(p)==0){
        p<-p
      }else{
        p<-(p/sum(p))
      }    
      AB<-rbind(AB,p)
    }
    AB<-AB[2:nrow(AB),] 
    # SECTION 3.3 CUMULATIVE LENGTH VECTOR 
    clength<-c(Length2[1])
    for(q in 2:numcol){
      cl<-Length2[q] + clength[q-1]
      clength<-c(clength, cl)
    }
    # SECTION 3.4 FINAL CALCULATION OF A  
    A<-AB[1,]
    for(k in (numbefore+1):(numsam+numbefore)){
      E<-AB[1,]
      for (j in 0:(depth[k-numbefore] - 1)){
        e<-colMeans(AB[(clength[k-1]-j+1):(clength[k] - j),])
        E<-rbind(E,e)
      }
      E<-E[2:nrow(E),] # this seems wrong but must be based on what I saw in matlab
      meanE<-colMeans(E)
      A <- rbind(A,meanE)
    }
    A<-A[2:nrow(A),]
    # SECTION 4: INVERSION
    I<-diag(numsam)
    dMeasr<-dMeas + r1*rnorm(numsam) #check other rnorm above
    NB<-numbefore
    Na<-numafter - 1
    mm <- matrix(nrow=numsam+numbefore+numafter-1,ncol=1,1)
    mm[1,]<-((maxratio-minratio)*runif(1))+minratio
    
    for(x in 2:(numsam+numbefore+numafter-1)){
      mm[x,]<-mm[x-1,]+stdev*rnorm(1)
      while (mm[x,]>maxratio){
        mm[x,]<-mm[x-1,]+stdev*rnorm(1)
      }
      while(mm[x,]<minratio){
        mm[x,]<-mm[x-1,]+stdev*rnorm(1)
      }
    }
    AA<-A%*%Conj(t(A))
    epsilon<-df*I
    AAep<-AA+epsilon
    library(matrixcalc)
    test<-matrix.power(AAep,-1)
    mEst<-mm + (Conj(t(A))%*%(test))%*%(dMeasr - A%*%mm)
    MEST[[n]]<-mEst 
    dPred<-A%*%mEst	
    dpe<-dPred - dMeas
    DPE[n] <- Conj(t(dpe))%*%dpe
    S[[n]]<- dpe           
  }
  
  dim <-nrow(MEST[[1]])
  
  # Plotting
  
  totallength<-matrix(1,nrow=numsam,ncol=1) 
  totallength[1,]<-Length2[1]  
  
  xx <- avelength*rep(1,numbefore)
  zz <- avelength*rep(1,numafter-1)
  totallength <- c(xx,totallength,zz)
  
  vec1 <- dMeasr[1]*rep(1,numbefore)
  vec2 <- dMeasr[numsam]*rep(1,numafter-1)
  dMeasd <- c(vec1,dMeasr,vec2)
  
  for(n in 2:(numsam+numbefore+numafter-1)){
    totallength[n] <- c(totallength[n-1]+Length2[n])
  }

  library(dplyr)
  temp<-MEST[[1]]
  for(i in 2:length(MEST)){
    temp<-cbind(temp,MEST[[i]])
  }
  mean<-rowMeans(temp)
  upper<-c()
  lower<-c()
  for(i in 1:nrow(temp)){
    lower[i]<-(quantile(temp[i,], c(0.025), type = 1,na.rm=TRUE))
    upper[i]<-(quantile(temp[i,], c(0.975), type = 1,na.rm=TRUE))
  }
  forPlot <- data.frame(totallength,mean,upper,lower)
  results_final<-list(DPE,forPlot)
  return(results_final)
}

