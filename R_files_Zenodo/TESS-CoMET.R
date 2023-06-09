
setwd('/Users/Emanuell/Documents/Carangaria/Revision/TESS/CoMET/')
library(TESS)
library(phytools)
library(geiger)
library(phangorn)


pruned.tree <- read.tree("myTree.tre")

#Run CoMET
numExpectedMassExtinctions <- 0
numExpectedRateChanges <- 2
samplingFraction = 0.9
# Specify the mean and standard deviation of the lognormal # prior on the speciation rate in real space 
speciationPriorMu <- 0.2
speciationPriorSigma <- 0.5
# Specify the mean and standard deviation of the lognormal # prior on the extinction rate in real space 
extinctionPriorMu <- 0.15
extinctionPriorSigma <- 0.5

# Transform the priors on the speciation rate into log space.
speciationRatePriorMean <- log((speciationPriorMu^2)
                               /sqrt(speciationPriorSigma^2+
                                       speciationPriorMu^2))
speciationRatePriorStDev <- sqrt( log(1+speciationPriorSigma^2
                                      /(speciationPriorMu^2)))
# Transform the priors on the extinction rate into log space.
extinctionRatePriorMean <- log((extinctionPriorMu^2)
                               /sqrt(extinctionPriorSigma^2+
                                       extinctionPriorMu^2))
extinctionRatePriorStDev <- sqrt( log(1+extinctionPriorSigma^2
                                      /(extinctionPriorMu^2)))

expectedSurvivalProbability <- 1

set.seed(123345)
tess.analysis(Tree,
              empiricalHyperPriors = FALSE,
              initialSpeciationRate = speciationPriorMu,
              speciationRatePriorMean = speciationRatePriorMean,
              speciationRatePriorStDev = speciationRatePriorStDev,
              initialExtinctionRate = extinctionPriorMu,
              extinctionRatePriorMean = extinctionRatePriorMean,
              extinctionRatePriorStDev = extinctionRatePriorStDev,
              samplingProbability = samplingFraction,
              samplingStrategy = "uniform",
              estimateNumberMassExtinctions = FALSE,
              MAX_ITERATIONS = 100000,
              BURNIN = 10000,
              MAX_TIME = 24*60*60,
              dir = "Uniform")



output <- tess.process.output("Uniform",
                              numExpectedRateChanges = numExpectedRateChanges,
                              numExpectedMassExtinctions = numExpectedMassExtinctions)

pdf('CoMET.pdf')
layout.mat <- matrix(1:4,nrow=2,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output,
                 fig.types = c("speciation rates","speciation shift times", "extinction rates", "extinction shift times"), las=2)
dev.off()
