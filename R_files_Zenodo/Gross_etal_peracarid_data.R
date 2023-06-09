library(picante)
library(vegan)
library(FD)
#Read in the species matrix
species_matrix <- read.csv("species_matrix.csv", header = T, row.names= 1)
species_matrix_sqrt <- sqrt(species_matrix)
Pacific_species_matrix <-  species_matrix_sqrt[c(1:4,13:18,22,23,26,27,31:34,41,42),]
Pacific_species_matrix <- Pacific_species_matrix[,-which(colSums(Pacific_species_matrix) == 0)]
Atlantic_species_matrix <- species_matrix_sqrt[-c(1:4,13:18,22,23,26,27,31:34,41,42),]
Atlantic_species_matrix <- Atlantic_species_matrix[,-which(colSums(Atlantic_species_matrix) == 0)]

#Read in the trait matrix
trait_matrix <- read.csv("trait_matrix.csv", header = T, row.names = 1)
#Code maximum fecundity as an ordered categorical trait.
trait_matrix$Maximum.fecundity..number.of.eggs. <- ordered(trait_matrix$Maximum.fecundity..number.of.eggs., levels = c("Very low","Low","Medium","High","Very high"))
#Code binary traits as numeric.
trait_matrix$Bioturbator <- as.numeric(trait_matrix$Bioturbator)
trait_matrix$eat.microalgae <- as.numeric(trait_matrix$eat.microalgae)
trait_matrix$eat.macroalgae <- as.numeric(trait_matrix$eat.macroalgae)
trait_matrix$eat.seagrass <- as.numeric(trait_matrix$eat.seagrass)
trait_matrix$eat.seagrass.detritus <- as.numeric(trait_matrix$eat.seagrass.detritus)
trait_matrix$suspension.feeder <- as.numeric(trait_matrix$suspension.feeder)
trait_matrix$detritivore.deposit.feeder <-  as.numeric(trait_matrix$detritivore.deposit.feeder)
trait_matrix$carnivore.parasite.scavenger <- as.numeric(trait_matrix$carnivore.parasite.scavenger)

#Read in the site data
site_data <- read.csv("site_data.csv", header = T)

#Calculate Gower distance matrix for all traits
trait_distances = as.matrix(gowdis(trait_matrix, ord =  "metric"))
#Calculate binary distance matrix for diet traits
feeding <- trait_matrix[,c(7:13)]
diet_distances <- as.matrix(dist(feeding, method = "binary", diag = T, upper = T))
rownames(diet_distances) <- rownames(trait_matrix)
colnames(diet_distances) <- rownames(trait_matrix)
diet_distances <- diet_distances[-c(40,103),-c(40,103)] #Remove taxa for which no data exists
#Calculate Gower distance matrix for microhabitat traits
microhab <- trait_matrix[,c(1,3:6)]
microhabitat_distances <- as.matrix(gowdis(microhab, ord = "metric"))
rownames(microhabitat_distances) <- rownames(trait_matrix)
colnames(microhabitat_distances) <- rownames(trait_matrix)

#Calculate standard effect sizes of mean pairwise distance based on the independent swap algorithm relative to the global species pool. 
set.seed(1993)
##Independent Swap model for all traits 
ses.traits.is.mpd <- ses.mpd(species_matrix_sqrt, trait_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for diet traits
ses.food.is.mpd <- ses.mpd(species_matrix_sqrt[,-c(40,103)], diet_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for microhabitat traits
ses.microhabitat.is.mpd <- ses.mpd(species_matrix_sqrt, microhabitat_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean pairwise distance based on the tip shuffle algorithm relative to the global species pool. 
##Tip Shuffle model for all traits 
ses.traits.ts.mpd <- ses.mpd(species_matrix_sqrt, trait_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for diet traits
ses.food.ts.mpd <- ses.mpd(species_matrix_sqrt[,-c(40,103)], diet_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for microhabitat traits
ses.microhabitat.ts.mpd <- ses.mpd(species_matrix_sqrt, microhabitat_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean nearest taxon distance based on the independent swap algorithm relative to the global species pool. 
set.seed(1993)
##Independent Swap model for all traits 
ses.traits.is.mntd <- ses.mntd(species_matrix_sqrt, trait_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for diet traits
ses.food.is.mntd <- ses.mntd(species_matrix_sqrt[,-c(40,103)], diet_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for microhabitat traits
ses.microhabitat.is.mntd <- ses.mntd(species_matrix_sqrt, microhabitat_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean nearest taxon distance based on the tip shuffle algorithm relative to the global species pool. 
##Tip Shuffle model for all traits 
ses.traits.ts.mntd <- ses.mntd(species_matrix_sqrt, trait_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for diet traits
ses.food.ts.mntd <- ses.mntd(species_matrix_sqrt[,-c(40,103)], diet_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for microhabitat traits
ses.microhabitat.ts.mntd <- ses.mntd(species_matrix_sqrt, microhabitat_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean pairwise distance based on the independent swap algorithm relative to the Pacific species pool. 
set.seed(1993)
##Independent Swap model for all traits 
ses.traits.is.mpd.pac <- ses.mpd(Pacific_species_matrix, trait_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for diet traits
ses.food.is.mpd.pac <- ses.mpd(Pacific_species_matrix[,-c(40,103)], diet_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for microhabitat traits
ses.microhabitat.is.mpd.pac <- ses.mpd(Pacific_species_matrix, microhabitat_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean pairwise distance based on the tip shuffle algorithm relative to the Pacific species pool. 
##Tip Shuffle model for all traits 
ses.traits.ts.mpd.pac <- ses.mpd(Pacific_species_matrix, trait_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for diet traits
ses.food.ts.mpd.pac <- ses.mpd(Pacific_species_matrix[,-c(40,103)], diet_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for microhabitat traits
ses.microhabitat.ts.mpd.pac <- ses.mpd(Pacific_species_matrix, microhabitat_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean nearest taxon distance based on the independent swap algorithm relative to the Pacific species pool. 
set.seed(1993)
##Independent Swap model for all traits 
ses.traits.is.mntd.pac <- ses.mntd(Pacific_species_matrix, trait_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for diet traits
ses.food.is.mntd.pac <- ses.mntd(Pacific_species_matrix[,-c(40,103)], diet_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for microhabitat traits
ses.microhabitat.is.mntd.pac <- ses.mntd(Pacific_species_matrix, microhabitat_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean nearest taxon distance based on the tip shuffle algorithm relative to the Pacific species pool. 
##Tip Shuffle model for all traits 
ses.traits.ts.mntd.pac <- ses.mntd(Pacific_species_matrix, trait_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for diet traits
ses.food.ts.mntd.pac <- ses.mntd(Pacific_species_matrix[,-c(40,103)], diet_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for microhabitat traits
ses.microhabitat.ts.mntd.pac <- ses.mntd(Pacific_species_matrix, microhabitat_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean pairwise distance based on the independent swap algorithm relative to the Atlantic species pool. 
set.seed(1993)
##Independent Swap model for all traits 
ses.traits.is.mpd.atl <- ses.mpd(Atlantic_species_matrix, trait_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for diet traits
ses.food.is.mpd.atl <- ses.mpd(Atlantic_species_matrix[,-c(40,103)], diet_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for microhabitat traits
ses.microhabitat.is.mpd.atl <- ses.mpd(Atlantic_species_matrix, microhabitat_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean pairwise distance based on the tip shuffle algorithm relative to the Atlantic species pool. 
##Tip Shuffle model for all traits 
ses.traits.ts.mpd.atl <- ses.mpd(Atlantic_species_matrix, trait_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for diet traits
ses.food.ts.mpd.atl <- ses.mpd(Atlantic_species_matrix[,-c(40,103)], diet_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for microhabitat traits
ses.microhabitat.ts.mpd.atl <- ses.mpd(Atlantic_species_matrix, microhabitat_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean nearest taxon distance based on the independent swap algorithm relative to the Atlantic species pool. 
set.seed(1993)
##Independent Swap model for all traits 
ses.traits.is.mntd.atl <- ses.mntd(Atlantic_species_matrix, trait_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for diet traits
ses.food.is.mntd.atl <- ses.mntd(Atlantic_species_matrix[,-c(40,103)], diet_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)
##Independent Swap model for microhabitat traits
ses.microhabitat.is.mntd.atl <- ses.mntd(Atlantic_species_matrix, microhabitat_distances, null.model = "independentswap", abundance.weighted = T, runs = 999, iterations = 1000)

#Calculate standard effect sizes of mean nearest taxon distance based on the tip shuffle algorithm relative to the Atlantic species pool. 
##Tip Shuffle model for all traits 
ses.traits.ts.mntd.atl <- ses.mntd(Atlantic_species_matrix, trait_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for diet traits
ses.food.ts.mntd.atl <- ses.mntd(Atlantic_species_matrix[,-c(40,103)], diet_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)
##Tip Shuffle model for microhabitat traits
ses.microhabitat.ts.mntd.atl <- ses.mntd(Atlantic_species_matrix, microhabitat_distances, null.model = "taxa.labels", abundance.weighted = T, runs = 999, iterations = 1000)