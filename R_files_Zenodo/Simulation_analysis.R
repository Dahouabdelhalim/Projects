##  CODE FOR THE ANALYSIS OF RESULTS FROM SIMULATIONS CODED IN SIMULATIONS.R
setwd(getwd());

#parameters
tmax <- 40 # number of generations
# get results from a file, compute power, and place in array
nrep <- 200; # number of simulation repeats
alphalist <- c(0,0.1,0.25,0.5,0.75,0.9);
nlist <- c(50) #c(25,50)
sigmaz <- 1   # standard deviation of trait
Slist <- c(0,0.025,0.05,0.1,0.2); 
sigma_thetalist <-  c(0,0.5,1,2,4) #standard deviation of optimum, removing value 0.25, instead of c(0,0.25,0.5,1,2,4)
SDbetamax <- max(Slist)*max(sigma_thetalist)
sigma_betalist <- seq(0,SDbetamax,length.out=length(sigma_thetalist)); # Same range of betas than with fluctuations of optimum, and same number

param <- T; # whether or not to compute means of estimated parameters
sampling <- F; # boolean for whether or not we want to sample from the posterior distribution


pow <- array(NA,c(length(alphalist),length(nlist),length(sigma_thetalist),length(Slist),4))
if (param) 
  {vtheta <-array(NA,c(length(alphalist),length(nlist),length(sigma_thetalist),length(Slist)))
   rho <-array(NA,c(length(alphalist),length(nlist),length(sigma_thetalist),length(Slist)))
   CIvtheta <-array(NA,c(length(alphalist),length(nlist),length(sigma_thetalist),length(Slist),2))
   CIrho <-array(NA,c(length(alphalist),length(nlist),length(sigma_thetalist),length(Slist),2))
   
   # same, but conditional on AR being the best model
   vthetacond <-array(NA,c(length(alphalist),length(nlist),length(sigma_thetalist),length(Slist)))
   rhocond <-array(NA,c(length(alphalist),length(nlist),length(sigma_thetalist),length(Slist)))
   CIvthetacond <-array(NA,c(length(alphalist),length(nlist),length(sigma_thetalist),length(Slist),2))
   CIrhocond <-array(NA,c(length(alphalist),length(nlist),length(sigma_thetalist),length(Slist),2))
 
  }


for (i in 1:length(alphalist))
{alpha<-alphalist[i]
  for (j in 1:length(nlist) )
  {mean_n <- nlist[j]
    for (k in 1:length(sigma_thetalist))
    {sigma_theta <- sigma_thetalist[k]
      for(l in 1:length(Slist))
      {S <- Slist[l]
       
      ifelse((sigma_theta==0) & (S!=0),
      filename <- paste("Constopt S=", S,", n=",mean_n,".txt",sep=""),
      ifelse(S==0,
             # case with no quadratic term
             filename <- paste("noquad_V(beta)=",sigma_betalist[k],", alpha=",alpha,", n=",mean_n,".txt",sep="") 
             # general case, with quadratic term and fluctuations
             ,filename <- paste("S=",S,", V(theta)=",round(sigma_theta,2),", alpha=",alpha,", n=",mean_n,".txt",sep="")))
      res <-dget(filename)

      # Power: proportion of simulations with chosen criterion. 
      # Reminder: order of models: list(AR,whitenoise,ARnoquad,ARnoquadwhitenoise,constopt)
      # criterion: Better model means more than 2 DIC points
      
      # Stabilizing selection: best model with quadratic term is better than best with no quadratic terms
      pickstab <-  apply(res$DIC[,c(1:2,5)],1,min) < apply(res$DIC[,3:4],1,min)-2
      # variable selection: best model with changing selection (including non-stabilizing) better than model with constant optimum
      pickvarsel <- apply(res$DIC[,1:4],1,min) < res$DIC[,5] - 2
      # moving optimum: best model has a moving optimum
      pickmoveopt <- apply(res$DIC[,1:2],1,min) < apply(res$DIC[,3:5],1,min) - 2
      # autocorrelated optimum: best model has AR optimum
      pickAR <- res$DIC[,1] < apply(res$DIC[,2:5],1,min) - 2

      pow[i,j,k,l,] <- unlist(lapply(list(pickstab,pickvarsel,pickmoveopt,pickAR),sum))/nrep
      # Mean parameter values across replicates where the best model is AR
      if (param) 
        { vtheta[i,j,k,l] <- mean(unlist(lapply(res$hyperpar,function (x) x[1,1])))
          rho[i,j,k,l] <- mean(unlist(lapply(res$hyperpar,function (x) x[2,1])))
          # table with CI interval limits as rows
          allCIV <- do.call(rbind,lapply(res$hyperpar,function (x) x[1,c(3,5)]))
            CIvtheta[i,j,k,l,] <- colMeans(allCIV)
          allCIrho <- do.call(rbind,lapply(res$hyperpar,function (x) x[2,c(3,5)]))
            CIrho[i,j,k,l,] <- colMeans(allCIrho)
        
        
        # values conditional on AR being the best model
        vthetacond[i,j,k,l] <- mean(unlist(lapply(res$hyperpar,function (x) x[1,1]))[pickAR])
        rhocond[i,j,k,l] <- mean(unlist(lapply(res$hyperpar,function (x) x[2,1]))[pickAR])
        # table with CI interval limits as rows
        allCIVcond <- do.call(rbind,lapply(res$hyperpar,function (x) x[1,c(3,5)]))[pickAR,]
            CIvthetacond[i,j,k,l,] <- colMeans(allCIV)
        allCIrhocond <- do.call(rbind,lapply(res$hyperpar,function (x) x[2,c(3,5)]))[pickAR,]
            CIrhocond[i,j,k,l,] <- colMeans(allCIrho)
        }
      }
    }
  }
}



# Figure 1A
# power do detect stabilizing selection as a function of gamma
# gray level = magnitude of fluctuations. Darker means more autocorrelation
par(mar=c(5,5,0.5,1))
par(cex.lab=.8, cex.axis=.7)
plot(NULL,NULL,xlim=c(0,max(Slist)),ylim=c(0,1),xlab=expression(paste("S=",1/(omega^2+sigma[z]^2))),ylab="Pr(stabilizing selection)")   #expression(Pr(widehat(S) != 0))
for (i in 1:length(alphalist))
{for (j in 1:length(nlist) )
{for (k in 1:length(sigma_thetalist))
{lines(Slist,pow[i,j,k,,1],lwd=(k-1)/2,col=gray(1-k/(length(sigma_thetalist)+1)))}}}  #,col=gray(1-j/(length(nlist)+1)),lwd=0.7*(sigma_thetalist[k]+1)

# Figure 1B
# Power do detect fluctuating selection as a function of  SD of optimum
# gray level = omega. Darker means smaller omega, and stronger selection.
par(mar=c(5,5,0.5,1))
par(cex.lab=.8, cex.axis=.7)
plot(NULL,NULL,xlim=c(min(sigma_thetalist),max(sigma_thetalist)),ylim=c(0,1),xlab=expression(sigma[theta]),ylab="Pr(fluctuating optimum)") #Pr(widehat(sigma[theta]) != 0)
for (i in 1:length(alphalist))
{for (j in 1:length(nlist) )
{for(l in 1:length(Slist))
 {lines (sigma_thetalist,pow[i,j,,l,3],lwd=l/2,col=gray(1 - l/(length(Slist)+1)))}}}

# Figure 1C
# power to detect autocorrelated fluctuations as a function of autocorrelation (only in cases with actual fluctuations)
# gray level = variance in optimum
par(cex.lab=.8, cex.axis=.7)
plot(NULL,NULL,xlim=c(min(alphalist),max(alphalist)),ylim=c(0,1),xlab=expression(alpha),ylab="Pr(AR1 optimum)") #expression(Pr(widehat(alpha) != 0))
for (j in 1:length(nlist) )
{for (k in 2:length(sigma_thetalist))
{for(l in 2:length(Slist))
{lines (alphalist,pow[,j,k,l,4],lwd=(k-1)/2,col=gray(1-(k-1)/length(sigma_thetalist)))}}}

#Figure 1D
# power to detect autocorrelated fluctuations as a proportion of all fluctuations (only in cases with actual fluctuations)
plot(NULL,NULL,xlim=c(min(alphalist),max(alphalist)),ylim=c(0,1),xlab=expression(alpha),ylab="Pr(AR1 | fluct. optimum)" ) #expression(Pr(widehat(alpha) != 0)/Pr(widehat(sigma[theta]) != 0))
for (j in 1:length(nlist) )
{for (k in 2:length(sigma_thetalist))
{for(l in 2:length(Slist))
{lines (alphalist,pow[,j,k,l,4]/(pow[,j,k,l,3]),lwd=(k-1)/2,col=gray(1-(k-1)/length(sigma_thetalist)))}}}  
#lines(0:1,0:1,lty=2,lwd=2)

# Figure 2A
# Estimated autocorrelation in optimum versus actual autocorrelation, across all simuls
# Color indicates SD of fluctuations
par(mar=c(5,5,4,2))
par(cex.lab=.8, cex.axis=.7)
plot(NULL,NULL,xlim=c(0,1),ylim=c(-.05,1),xlab=expression(alpha),ylab=expression(widehat(alpha)))
for (j in 1:length(nlist) )
{for (k in 2:length(sigma_thetalist))
{for(l in 1:length(Slist))
{lines(alphalist,rho[,j,k,l],lwd=(k-1)/2,col=gray(1-(k-1)/(length(sigma_thetalist))))
 #points(alphalist,CIrho[,j,k,l,1],lty=4,lwd=1,col=gray(1-k/(length(sigma_thetalist)+1)))
 #points(alphalist,CIrho[,j,k,l,2],lty=4, lwd=1,col=gray(1-k/(length(sigma_thetalist)+1)))
}}}
lines(0:1,0:1,lty=2,lwd=2)

# Figure 2B
# width of CI of alpha as a function of SD of fluctuations. Color indicates alpha
par(mar=c(5,5,4,2))
par(cex.lab=.8, cex.axis=.7)
plot(NULL,NULL,xlim=c(min(sigma_thetalist),max(sigma_thetalist)),ylim=c(0,2),xlab=expression(sigma[theta]),ylab=expression(CI(widehat(alpha))))
for (i in 1:length(alphalist))
{for (j in 1:length(nlist) )
{for(l in 2:length(Slist))
{lines(sigma_thetalist,CIrho[i,j,,l,2]-CIrho[i,j,,l,1],lwd=i/2,col=gray(1-i/(length(alphalist)+1)))}}}
