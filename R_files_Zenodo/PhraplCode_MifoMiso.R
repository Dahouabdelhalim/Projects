library(ape)
library(phytools)
library(phrapl)

setwd("WORKING/DIRECTORY")

# Phrapl ------------------------------------------------------------------
## In this example 3 populations and 3 free parameters are considered.
data(MigrationArray_3pops_maxMigrationK1)
migrationArray<-migrationArray[1:16]
migrationArrayMap<-GenerateMigrationArrayMap(migrationArray)

#Load my data
assignFile<-read.delim("cladeAssignments_newPipeline.txt", quote="", sep=" ")
trees<-read.tree(file="MifoMisoRooted_newPipeline.trees")

#Subsample trees, 10 iterations per locus
popAssignments<-list(c(2,3,3))

subsampledTrees<-PrepSubsampling(popAssignments=popAssignments, assignmentsGlobal=assignFile, observedTrees=trees, subsamplesPerGene=100, outgroup=T, outgroupPrune=T)

#Calculate weights for each subsampled tree
subsampleWeights<-GetPermutationWeightsAcrossSubsamples(popAssignments=popAssignments, observedTrees = subsampledTrees)

# Initial values for diverge time (tau)
collapseStarts=c(0.2, 1, 4, 8)

# Initial values for migration rate (m)
migrationStarts=c(0.01, 0.05, 0.1)

rm(result)
rm(gridList)
rm(totalData)
startTime<-as.numeric(Sys.time())
result<-GridSearch(migrationArray=migrationArray,
                   popAssignments=popAssignments,
                   subsampleWeights.df=subsampleWeights,
                   doWeights=F,
                   nTrees=10000,
                   observedTrees=subsampledTrees,
                   totalPopVector = c(6,12,14),
                   migrationStarts = migrationStarts,
                   collapseStarts = collapseStarts,
                   msPath="/Applications/msdir/ms",
                   print.results=TRUE,
                   debug=TRUE)
#Print summary results
print(result[[1]])

#Make dedicated grid list
gridList<-result[[1]]

#Get elapsed time
stopTime<-as.numeric(Sys.time()) #stop system timer
elapsedSecs<-stopTime - startTime #elapsed time in hours
elapsedHrs<-(elapsedSecs / 60) / 60 #convert to hours
elapsedDays<-elapsedHrs / 24 #convert to days

## Save results to an R object file
save(list=ls(), file="MifoMisoPhraplOutput_newPipeline.rda")
result$overall.results
############STEP5:Cull, print, and visualize results
#Generate an output table
rm(totalData)
totalData<-ConcatenateResults(rdaFiles="MifoMisoPhraplOutput_newPipeline.rda")

#Calculate model averaged parameter values
CalculateModelAverages(totalData,parmStartCol=9)
library(qpcR)
akaike.weights(totalData$AIC)

#Plot a 3-D movie of the best model
PlotModel(
  migrationIndividual=migrationArray[[13]],
  taxonNames=c("MifoN","MifoS","Miso"))
