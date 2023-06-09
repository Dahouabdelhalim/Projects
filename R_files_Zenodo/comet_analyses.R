
library(ape)
library(laser)
library(dplyr)
library(phytools)
library(TreePar)
library(TESS)

budI <- read.tree("bivalve_family_tree_budI_dates.tre")
budII <- read.tree("bivalve_family_tree_budII_dates.tre")
bifI <- read.tree("bivalve_family_tree_bifI_dates.tre")
bifII <- read.tree("bivalve_family_tree_bifII_dates.tre")

budI <- drop.tip(budI, c("Scaphopoda", "Gastropoda", "Cephalopoda", "Chitonidae", "Monoplacophora"))
budI <- multi2di(budI)
budI$edge.length[budI$edge.length==0] <- 0.001
budI <- force.ultrametric(budI)

budII <- drop.tip(budII, c("Scaphopoda", "Gastropoda", "Cephalopoda", "Chitonidae", "Monoplacophora"))
budII <- multi2di(budII)
budII$edge.length[budII$edge.length==0] <- 0.001
budII <- force.ultrametric(budII)

bifI <- drop.tip(bifI, c("Scaphopoda", "Gastropoda", "Cephalopoda", "Chitonidae", "Monoplacophora"))
bifI <- multi2di(bifI)
bifI$edge.length[bifI$edge.length==0] <- 0.001
bifI <- force.ultrametric(bifI)

bifII <- drop.tip(bifII, c("Scaphopoda", "Gastropoda", "Cephalopoda", "Chitonidae", "Monoplacophora"))
bifII <- multi2di(bifII)
bifII$edge.length[bifII$edge.length==0] <- 0.001
bifII <- force.ultrametric(bifII)

#############################
########### TESS ############
#############################

rho <- .85

ml.bd.estimate <- birthdeath(budI)
mu_mu = 0.05
mu_lambda = mu_mu + ml.bd.estimate$para[2]
std_lambda = 0.02
std_mu = 0.02

numExpectedMassExtinctions = 5
numExpectedRateChanges = 5

pMassExtinctionPriorShape2 <- 100
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 * expectedSurvivalProbability / (expectedSurvivalProbability - 1)


### Budding 1

tess.analysis( tree=budI,
               empiricalHyperPriors = TRUE,
               numExpectedRateChanges = numExpectedRateChanges,
               samplingProbability = rho,
               numExpectedMassExtinctions = numExpectedMassExtinctions,
               BURNIN = 2000,
               MAX_ITERATIONS = 10000,
               THINNING = 50,
               dir = "COMET/comet_budI_hp")


### Budding 2

tess.analysis( tree=budII,
               empiricalHyperPriors = TRUE,
               numExpectedRateChanges = numExpectedRateChanges,
               samplingProbability = rho,
               numExpectedMassExtinctions = numExpectedMassExtinctions,
               BURNIN = 2000,
               MAX_ITERATIONS = 10000,
               THINNING = 50,
               dir = "COMET/comet_budII_hp")


### Bifurcating I

tess.analysis( tree=bifI,
               empiricalHyperPriors = TRUE,
               numExpectedRateChanges = numExpectedRateChanges,
               samplingProbability = rho,
               numExpectedMassExtinctions = numExpectedMassExtinctions,
               BURNIN = 2000,
               MAX_ITERATIONS = 10000,
               THINNING = 50,
               dir = "COMET/comet_bifI_hp")


### Bifurcating II


tess.analysis( tree=bifII,
               empiricalHyperPriors = TRUE,
               numExpectedRateChanges = numExpectedRateChanges,
               samplingProbability = rho,
               numExpectedMassExtinctions = numExpectedMassExtinctions,
               BURNIN = 2000,
               MAX_ITERATIONS = 10000,
               THINNING = 50,
               dir = "COMET/comet_bifII_hp")



output.budI <- tess.process.output("COMET/comet_budII_hp",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)
output.budII <- tess.process.output("COMET/comet_budII_hp",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)
output.bifI <- tess.process.output("COMET/comet_bifI_hp",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)
output.bifII <- tess.process.output("COMET/comet_bifII_hp",
                numExpectedRateChanges = numExpectedRateChanges,
                numExpectedMassExtinctions = numExpectedMassExtinctions)

layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output,
fig.types = c("speciation rates",
	"speciation shift times",
	"extinction rates",
	"extinction shift times",
	"mass extinction Bayes factors",
	"mass extinction times"),
las=2)

png("COMETresults.png", width = 8, height=10, units="in", res=500)
par(mfrow=c(4,2))
par(mai=c(0.8,0.8,0.3,0.3))
tess.plot.output(output.budI, fig.types=c("speciation rates", "speciation shift times"))
tess.plot.output(output.budII, fig.types=c("speciation rates", "speciation shift times"))
tess.plot.output(output.bifI, fig.types=c("speciation rates", "speciation shift times"))
tess.plot.output(output.bifII, fig.types=c("speciation rates", "speciation shift times"))
dev.off()


