
PasseyEMeas1_1<-function(Length,dMeas,la,numtrials){
  allTrials<-list()
  Edist<-c()
  dMeasError<-c()
  # Equivalent to user entered data section from Passey et al. 2005
  numtrials<-numtrials1
  Length <-Length1
  dMeas <- dMeas1
  numsam<-length(dMeas)
  r1 <- .05
  r2 <- 1
  r3 <- 5
  la<-la1
  # This section of the code generates randomized error around the dMeas and Length values
  for(i in 1:numtrials){
    Isoslope<-rep(1,numsam)# check here
    for(n in 1:numsam){
      if(n==1){
        Isoslope[n]<-((abs(dMeas[n+1]-dMeas[n]))/(0.5*(Length[n+1]+Length[n])))
      }else{
        if(n==numsam){
          Isoslope[n]<-((abs(dMeas[n]-dMeas[n-1]))/(0.5*(Length[n]+Length[n-1])))
        }else{
          Isoslope[n]<-(((((abs(dMeas[n+1] - dMeas[n])) / ( 0.5*(Length[n+1]+Length[n])))   +   ((abs(dMeas[n] - dMeas[n-1])) / (0.5* (Length[n]+Length[n-1]))))/2))
        }
      }
    }
    RElength<-r2*rnorm(numsam)
    LengthError<-Isoslope*RElength
    LEdMeas<-LengthError + dMeas
    # This section makes sure the error values generated above are sensical
    # Checks that they are not above or below the highest or lowest adjacent values
    for(n in 1:numsam){
      if(n==1){
        if(LEdMeas[n]<min(dMeas[n:n+1])){
          LEdMeas[n]<-min(dMeas[n:n+1])
        }else{
          if(LEdMeas[n]>max(dMeas[n:n+1])){
            LEdMeas[n]<-max(dMeas[n:n+1])
          }else{
            LEdMeas[n]<-LEdMeas[n]
          }
          if(n==numsam){
            if(LEdMeas[n]<min(dMeas[n:(n-1)])){
              LEdMeas[n]<-min(dMeas[n:(n-1)])
            }else{
              if(LEdMeas[n]>max(dMeas[n:(n-1)])){
                LEdMeas[n]<-max(dMeas[n:(n-1)])
              }else{
                LEdMeas[n]<-LEdMeas[n]
              }             
              if(LEdMeas[n]<min(dMeas[n-1:(n+1)])){ 
                LEdMeas[n]<-min(dMeas[n-1:(n+1)])
              }else{
                if(LEdMeas[n]>max(dMeas[n-1:(n+1)])){
                  LEdMeas[n]<-max(dMeas[n-1:(n+1)])
                }else{
                  LEdMeas[n]<-LEdMeas[n]
                }
              }
            }
          }
        }
      }
    }
    LengthError<-LEdMeas - dMeas
    # This section calculates the depth-dependent isotope error
    # This is copied from the code associated with Passey et al. 2005
    # "The following section calculates the depth-dependent isotope error.  
    #It fits a cubic spline to the the measured data dMeas, and creates a vector of interpolated delta values for each unit length on the x-axis. 
    #A new vector is created that is shifted la %units, and the two vectors are subtracted to give DELTA-delta values reflecting the difference between the isotope ratio at the outside %enamel surface and the enamel-dentine junction.  These values are then multiplied by the sampling depth uncertainty to give a per mil uncertainty."
    totallength<-rep(1,numsam)
    totallength[1]<-Length[1]
    for(n in 2:numsam){
      totallength[n]<-totallength[n-1]+Length[n]
    }
    xx<-totallength[1]:totallength[numsam]
    library(pracma)
    dInterp <- try(interp1(totallength,dMeas, xx, method = "spline"),silent=T) # This is the closest equivalent of the spline function from Matlab in R
    if(inherits(dInterp, "try-error")){ 
      xx[1]<-xx[1]+0.0001
      xx[length(xx)]<-xx[length(xx)]-0.0001
      dInterp <- interp1(totallength,dMeas, xx, method = "spline")
    }
    addbefore<-dMeas[1]*rep(1,la)
    addafter<-dMeas[numsam]*rep(1,la)
    dInterpShift<-c(addbefore,t(dInterp))
    dInterp<-c(t(dInterp),addafter)
    Dd<-dInterpShift-dInterp
    Deltadelta<-rep(1,numsam)
    for(n in 1:numsam){ 	
      Deltadelta[n]<-Dd[totallength[n]]
    } 
    REdepth<-(r3*rnorm(numsam))/la
    DepthError<-Deltadelta*REdepth
    AnalysisError<-r1*rnorm(numsam)
    SqSumError<-(LengthError^2 + DepthError^2 + AnalysisError^2)^0.5
    E<-Conj(t(SqSumError))%*%SqSumError 
    Edist[i]<-E 
    allTrials[[i]]<-SqSumError
    dMeasError<-dMeas + SqSumError
  }
  totallength<-rep(1,numsam)
  totallength[1]<-Length[1]
  for(n in 2:numsam){
    totallength[n]<-totallength[n-1]+Length[n]
  }
  forPlot <- list(Edist,totallength,allTrials,dMeas,dMeasError)
  return(forPlot)
}

