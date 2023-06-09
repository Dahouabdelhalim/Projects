##Creating the BAMMPriors
#Open libraries
library(ape)
library(BAMMtools)
#Read tree file and trait file, and setting the BAMM priors
phy <- read.tree("C_tree.nwk")
setBAMMpriors(phy = phy, traits = "Labellumlength.txt", outfile = "myPriors.txt",
              Nmax = 1000,
              suppressWarning = FALSE)

##Creating the control file
generateControlFile(file = 'divcontrol.txt', type = 'trait', params = list(
  treefile = 'C_tree.nwk',
  traitfile= 'Labellumlength.txt',
  numberOfGenerations = '1000000000',
  overwrite = '1',
  betaInitPrior = '0.25073301023299', 
  betaShiftPrior = '0.108641231382907',
  useObservedMinMaxAsTraitPriors = '1',
  expectedNumberOfShifts = '1'))

##POST-BAMM RUN ANALYSIS##
##Creating the bammdata object
edata <- getEventData(phy, eventdata = "event_data.txt", burnin=0.1, type = "trait")

#Assessing MCMC convergence
mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

##Discarding some run as burnin
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

##Checking the effective sample sizes of the log-likelihood and the number of shift events present in each sample
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

#Finding out the number of rate shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
post_probs['X'] / post_probs['Y']

##Summarizing the posterior distribution of the number of shifts using summary methods
library(BAMMtools)
phy <- read.tree("C_tree.nwk")
edata <- getEventData(phy, eventdata = "bammrun_eventdata.txt", burnin=0.1)
shift_probs <- summary(edata)

##Computing the Bayes Factor
postfile <- "post_mcmc_out.txt"
bfmat <- computeBayesFactors(mcmcout, expectedNumberOfShifts=1, burnin=0.1)
bfmat

##Mean phylorate plot
plot.bammdata(edata, lwd=2)
plot.bammdata(edata, lwd=2, legend=T)
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)
plot.credibleshiftset(css)

##Macroevolutionary cohort analysis
cmat <- getCohortMatrix(edata)
cohorts(cmat, edata)
