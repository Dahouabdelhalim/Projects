####################################################################
##
## KittiwakeLifeTrajecPlot: Creates a plot of sample life trajectories
## for average and below-average quality kittwakes.
##
#######################################################################

rm(list=ls(all=TRUE)); graphics.off(); 
require(Matrix); require(statmod);

path ="/home/robin/ipm/drivers/Manuscripts/Luckpartition/code/"
out=try(setwd(path),silent=TRUE);

if(class(out)=="try-error") { 
    path=ifelse(.Platform$OS.type=="windows","c:/repos/drivers/Manuscripts/luckPartitioning","~/repos/drivers/Manuscripts/LuckPartitioning"); 
    setwd(path); 
  }
getwd(); 

source ("KittiwakeMatrices.R")
source("Standard Graphical Pars.R")
source("Utilities.R") 
source("MegamatrixFunctions.R");

## don't run all the code if we're just debugging
debugging = FALSE

## How many lives do we want to look at for each trait?
numBirds = 5

## Maximum lifespan
maxT = 25

mz = 25;   # size of the age/stage transition matrix 

## Prob. of breeding (is 1 for breeding classes and zero otherwise
## because the real probability of breeding is given by the transition
## rates into the breeding classes.)
pb = c(0, 0, 0, 1, 1) # prob. of breeding
pbMegavec = rep(pb, 5)

## Variance in per capita # of recruits. 
## Assumes those with 2 or 3 recruits have equal prob of 2 or 3. 
sigbsq = c(0, 0, 0, 0, 0.25)
sigbsqMegavec = rep(c(0, 0, 0, 0, 0.25), 5)

## Make cross-classified vectors
## expected per capita # of recruits
bMegavec = rep(b, 5)

## Set up birth size distribution - all are immature, age 1.  
c0 = c(1, rep(0, mz-1));

## function to calculate P.  THIS INCLUDES A STATE FOR THE DEAD. ######
makeP = function (newpars, breedingProbFactor) {
  
  ## Modify reproductive stage transition matrices
  newM3 = M3; newM4 = M4; newM5 = M5; 
  
  newM5[4:5,] = newM5[4:5,] * breedingProbFactor
  for (j in 1:5) {
    B = sum(newM5[4:5,j]); A = sum(newM5[1:3,j]); bfac = (1-B)/A; 
    newM5[1:3,j]=newM5[1:3,j]*bfac; 
  }
  
  s = c(1, exp(newpars)/(1 + exp(newpars)))
  P = matrix (0, 26, 26)
  P[5*1 + 1:5, 5*0 + 1:5] = M1 * matrix (rep(s1, 5), 5, 5, byrow=TRUE)
  P[5*2 + 1:5, 5*1 + 1:5] = M2 * matrix (rep(s2, 5), 5, 5, byrow=TRUE)
  P[5*3 + 1:5, 5*2 + 1:5] = newM3 * matrix (rep(s3, 5), 5, 5, byrow=TRUE)
  P[5*4 + 1:5, 5*3 + 1:5] = newM4 * matrix (rep(s, 5), 5, 5, byrow=TRUE)
  P[5*4 + 1:5, 5*4 + 1:5] = newM5 * matrix (rep(s, 5), 5, 5,
                                            byrow=TRUE)
  ## Probability of dying
  P[26,1:25] = 1 - apply (P[1:25, 1:25], 2, sum)
  ## Sometimes the above subtraction results in a terribly small
  ## negative number.
  foo = which ((P[26,] < 0) & (P[26,] > -1e-10))
  P[26, foo] = 0
  ## The dead stay dead.
  P[26, 26] = 1

  return (P)
}

## estimate CV of adult survivals from Cam et al. 2002, Fig. 1
q75 = log(0.9/0.1);   #75th percentile on inverse-logit scale ("quality")
q25 = log(0.7/0.3)    #25th percentile on inverse-logit scale ("quality")
IQR=q75-q25;
sigma=IQR/1.35;
mu = (q75+q25)/2;
CVSurv = sigma/mu;
sdSurvQuality = CVSurv * pars

## What is the estimated real CV for transitions to breeding classes?
## estimate CV of breeding prob. from Cam et al. 2002
q25 = 0.94
q75 = 0.98
IQR=q75-q25;
sigma=IQR/1.35;
mu = (q75+q25)/2;
CVBreedingProb = sigma/mu;
sdBreedingProb = sigma

## Record the trajectory of a life
trajecRecord = function (newpars, breedingProbFactor) {
  P = makeP (newpars, breedingProbFactor)
  ## Store reproductive stage (1--5) and clutch size vs time
  stage = clutchSize = rep (0, maxT+1)
  ## everyone starts as state 1 (pre-breeder) with no kids
  stage[1] = 1; clutchSize[1] = 0
  ## state and nextState use the full age x stage classification
  state = 1  ## age 1 pre-breeder
  ## You start off not dead.
  dead = FALSE
  
  t = 1
  while (!dead) {
    t = t + 1
    nextState = sample (1:26, 1, prob=P[,state])
    if (nextState == 26) {
      dead = TRUE
      lifespan = t - 1  ## minus 1 because age = t - 1
    } else {
      ## record reproductive stage
      stage[t] = (nextState-1) %% 5 + 1
      ## record clutch size
      if (stage[t] == 4) {
        clutchSize[t] = 1
      } else if (stage[t] == 5) {
        clutchSize[t] = rbinom (1, size=1, prob=0.5) + 2
      }
    }
    state = nextState
  }  ## end while loop
  out = list (stage=stage, clutchSize=clutchSize, lifespan=lifespan)
  return (out)
}

## place to record lives
averageBirdsStage = averageBirdsClutchSize = matrix (0, numBirds,
                                                     maxT+1)
wimpyBirdsStage = wimpyBirdsClutchSize = matrix (0, numBirds, maxT+1)
wimpyBirdsLifespan = averageBirdsLifespan = rep (0, numBirds)

## First, some trajectories for average birds
for (ii in 1:numBirds) {
  out = trajecRecord (pars, 1)
  averageBirdsStage[ii,] = out$stage
  averageBirdsClutchSize[ii,] = out$clutchSize
  averageBirdsLifespan[ii] = out$lifespan
}

## And now some trajectories for below-average quality birds
newpars = pars - sdSurvQuality
## Below-average quality in breeding probability too.
breedingProbFactor = 1 - 5*CVBreedingProb
for (ii in 1:numBirds) {
  out = trajecRecord (newpars, breedingProbFactor)
  wimpyBirdsStage[ii,] = out$stage
  wimpyBirdsClutchSize[ii,] = out$clutchSize
  wimpyBirdsLifespan[ii] = out$lifespan
}


############### Make plot ##############################

par (cex.lab=1.4)
maxLifespan = max(averageBirdsLifespan)
plot (0:(maxLifespan-1), rep(0, maxLifespan), ylim=c(0, 11),
      type="n", yaxt="n", xlab="Age", ylab="Reproductive stage")
axis (2, at=1:10, labels=rep(1:5,2))
colorVec = c("red", "orange", "green", "blue", "violet")

offset = -0.4
for (ii in 1:numBirds) {
  lines (0:(averageBirdsLifespan[ii]-1),
         averageBirdsStage[ii,1:averageBirdsLifespan[ii]] + 5 +
         offset, col=colorVec[ii])
  birthAges = which(averageBirdsClutchSize[ii,] > 0)
  points (birthAges-1, averageBirdsStage[ii,birthAges] + 5 + offset,
          pch=as.character(averageBirdsClutchSize[ii,birthAges]),
          col=colorVec[ii])
  offset = offset + 0.2
}

offset = -0.4
for (ii in 1:numBirds) {
  lines (0:(wimpyBirdsLifespan[ii]-1),
         wimpyBirdsStage[ii,1:wimpyBirdsLifespan[ii]] + offset,
         col=colorVec[ii], lty=2)
  birthAges = which(wimpyBirdsClutchSize[ii,] > 0)
  points (birthAges-1, wimpyBirdsStage[ii,birthAges] + offset,
          pch=as.character(wimpyBirdsClutchSize[ii,birthAges]),
          col=colorVec[ii])
  offset = offset + 0.2
}
