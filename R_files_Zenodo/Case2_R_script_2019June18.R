# this script is partof the manuscript 'The Cyclostratigraphy Intercomparison Project (CIP): a consistency test for cyclostratigraphy and astrochronology' by M. Sinnesael et al. 

# ! Note that noise generation was not fixed, and that therefore this script deviates from the original CIP case. However, the Imbire modes is used as for the original CIP case.

# load the 'astrochron' package, which must be installed
library(astrochron)

#set time steps
Dt <- 1

# get orbital solution
La04 <- getLaskar()

# extract obliquity and precession
Obl <- etp(solution = La04, tmin=0, tmax=6000, dt=Dt, eWt=0, oWt=1, pWt=0)
Prec <- etp(solution = La04, tmin=0, tmax=6000, dt=Dt, eWt=0, oWt=0, pWt=1)

# compile p-0.5T mimicking NH insolation
PminhalfT <- Prec
PminhalfT[,2] <- (-Prec[,2]+.5*Obl[,2])

# scale to mimick insolation (needed for the Imbrie model)
PminhalfT <- s(PminhalfT)
PminhalfT[,2] <- PminhalfT[,2]*11.66+500

# apply the II model as in LRo4
# set parameters
b1=.6
b2=.3
Tm1=15
Tm2=5
Tc=.5
Dt=1


# initialize the I&I model by parameter setting for the 'imbrie' function
# set b and Tm parameters
# set b parameter
btpts <- matrix(NA,4,2)
btpts[1,1] <- 0
btpts[1,2] <- b1
btpts[2,1] <- 1500
btpts[2,2] <- b1
btpts[3,1] <- 3000
btpts[3,2] <- b2
btpts[4,1] <- 6000
btpts[4,2] <- b2

b <- linterp(btpts,Dt, verbose=F, genplot=F)

# set Tm parameter
Tmtpts <- matrix(NA,4,2)
Tmtpts[1,1] <- 0
Tmtpts[1,2] <- Tm1
Tmtpts[2,1] <- 1500
Tmtpts[2,2] <- Tm1
Tmtpts[3,1] <- 3000
Tmtpts[3,2] <- Tm2
Tmtpts[4,1] <- 6000
Tmtpts[4,2] <- Tm1

Tm <- linterp(Tmtpts,Dt, verbose=F, genplot=F)

# run Imbrie model
IIM<- imbrie(PminhalfT,b=b[,2],Tm=Tm[,2],output=T,genplot=1)

# plot result again (can be omitted)
dev.off()
plot(s(IIM[1:1500,]), type="l", col='red')

# standardize model
sIIM <- s(IIM)
# create AR1 noise
# ! random noise will differ fro run to run, therefore CIP case 2 can not be perfectly reproduced
noise <- ar1(npts=5500,dt=1,mean=0,sdev=1,rho=0.9,shuffle=F,nsim=1,genplot=F,verbose=T)
# create data series from 500-1000 and 1200-1500 ka
sIIM <- rbind(sIIM[500:1000,], sIIM[1200:1500,])
# place on new 'depth' scale
sIIM <- (cb(1:802,-sIIM[,2]))
# plot for check
plot(sIIM, type="l", xlab="Time in ka", ylab= "Imbrie & Imbrie ice model")
# plot line at time break
abline(v=501, lwd=2, col="red")
#apply ramping sedimentation rate
srII23 <- sedRamp(sIIM, srstart=1, srend=,1.5)
srII23 <- s(linterp(srII23, dt=1))
# plot for check
plot(srII23, type="l", xlab="Depth in cm", ylab= "Imbrie & Imbrie ice model", main="Sedimentation rate ramping from 1-1.5 cm/ka")

# combine model with noise
srIInoise <- cb(srII23[,1], srII23[,2]+noise[1:1002,2]*.6)

# plot for check
dev.off()
plot(srIInoise, type="l", xlab="Depth in cm", ylab="Artificial proxy data", main="Signal2, different noise representation")
abline(v=579, lwd=2, col="red")

# end script