### Full R code for main-text results
###	From Evolution submission: "Testing adaptive radiation: a new test with Madagascar frogs"

### All analyses done on R version 4.0.2

### Load libraries used in this file
	library(tidyverse) # Version 1.3.0
	library(ape) # Version 5.4-1
	library(devtools) # Version 2.3.0
	# devtools::install_github('geomorphR/geomorph', ref = 'Develop') # Body-size analyses will only work with this version of 'geomorph' (Fall 2020)
	library(geomorph) # Version 3.3.1
	library(phytools) # Version 0.7-70
	library(geiger) # Version 2.0.7
	library(motmot) # Version 2.1.3

### Load data, functions
	setwd('') # You will need to set the directory to wherever you downloaded the data and R code appendices
	# Interspecific data for analysis (you will first need to save the data tabs of these as CSV files)
		mdat <- read.csv('SuppInfoS6.interspecific.residuals.foranalysis.csv',row.names = 1)
		pdat <- mdat[!is.na(mdat$cling.res),] # Making the performance data its own data frame
		mdat <- mdat[,1:13] # Stripping the performance variables from the morphology data frame
	# Phylogenies
		mtree <- read.tree('SuppInfoS4.217-species.phylogeny.tre')
		ptree <- read.tree('SuppInfoS5.80-species.phylogeny.tre')
	# Functions
		source('SuppInfoS7.AR.singlefunction.R') # Load power simulation function based on Adams 2014 Syst. Biol., Adams and Collyer 2018 Syst. Biol.
		source('SuppInfoS11.helper.functions.R') # Load functions used in analyses to assist with iterating analyses over datasets, as well as running and summarizing subsampling rate simulations. 
			# Also includes functions published by Adams 2014 Syst. Biol. 63:166–177 (phylogenetic transformation) and Moen et al. 2013 (two-block
			# partial least-squares analysis). 
	# Define analysis variables
		mvars <- grep('.res', colnames(mdat), value=T)[c(1,9,2:8)]
		pvars <- colnames(pdat)[14:ncol(pdat)]


###################################################################################################################################
### Part 1: Ecomorphology analyses (PGLS of microhabitat with morphology and performance; 2B-PLS of morphology and performance) ###
###################################################################################################################################

### A. Analyze fit of microhabitat, body shape, and performance

	# Modify data frame for analysis
			pp <- pdat
		# Mantellidae factor
			pp$mantellidae <- as.factor(ifelse(pp$family=='Mantellidae','yes','no'))
		# Make a merged microhabitat vector, given insufficient numbers of species per microhabitat within Mantellidae (explained in Supp. Info. S1)
			pp$ml <- pp$micro
			pp$ml[pp$ml == 'semi.aquatic'] <- 'aquatic'
			pp$ml[pp$ml == 'semi.burrowing'] <- 'burrowing'
			pp$ml[pp$ml == 'semi.arboreal'] <- 'arboreal'
			pp$ml <- as.factor(as.character(pp$ml))

	# 80-species dataset, fit of microhabitat and performance
		response <- pvars
		tree <- ptree
		dataset <- pp
		results <- procD.pgls(formula(paste(paste(response, collapse=" + "), paste(c('ml','mantellidae','ml:mantellidae'), 
					collapse=" + "), sep="~")), phy = tree, data = dataset)
		summary(results)

	# 80-species dataset, fit of microhabitat and body shape
		response <- mvars
		tree <- ptree
		dataset <- pp
		results <- procD.pgls(formula(paste(paste(response, collapse=" + "), paste(c('ml','mantellidae','ml:mantellidae'), 
					collapse=" + "), sep="~")), phy = tree, data = dataset)
		summary(results)
	
	# Clean up workspace
		rm(response, tree, dataset, results)

### B. Plot relationship among performance and body-shape variables (Fig. 4 a-d; all plots done with GGplot, then touched up in Adobe Illustrator)

	# Set up for GGplot
		pp <- pp # Just to remind ourselves that we are using the modified performance dataframe above.
		#	That data frame has the merged microhabitat factor, which we also use for easier visualization.
		# Set 'type' to indicate different symbols based on biogeography and clade
			pp$type <- 'all.other'
			pp$type[pp$family == 'Mantellidae'] <- 'Mantellidae'
			pp$type[pp$family != 'Mantellidae' & pp$location == 'Madagascar'] <- 'other.Malagasy'
		# Turn into tibble for somewhat easier manipulation
			pp$species <- rownames(pdat)
			pp <- tibble(pp) %>% select(species, type, ml, micro:swvel.res)			
		# Define colors and shapes
			colors <- c(aquatic = 'dodgerblue3', arboreal = 'forestgreen', burrowing = 'brown4', 
						terrestrial = 'darkgoldenrod1', torrential = 'black')
			shapes <- c(all.other = 21, Mantellidae = 23, other.Malagasy = 22)
		# Make plotting auxiliary function to cover the standard elements shared across plots
			frog.ecomorph.plot.extras <- function(){
				list(
					scale_fill_manual(
						name = 'Ecomorph',
						values = colors,
						labels = c('Aquatic / semiaquatic','Arboreal / semiarboreal','Burrowing / Semiburrowing','Terrestrial','Torrential')
					),
					scale_shape_manual(
						name = 'Origin',
						values = shapes,
						labels = c('All other frogs','Mantellidae','Other Malagasy frogs')
					),
					guides(fill = guide_legend(override.aes = list(color = colors))),
					coord_fixed(), # Ensure both axes have the same scaling, since plotting either residuals or PC axes
					theme_classic(),
					theme(
						axis.text = element_text(color = 'black'),
						axis.text.y = element_text(angle = 90, hjust = 0.5),
						axis.ticks = element_line(color = 'black'),
						panel.background = element_rect(color = 'black')
					)
				)
			}


	# Fig. 4a: Swimming velocity vs. jumping velocity
		ggplot(pp, aes(jpvel.res, swvel.res, fill = ml, shape = type)) +
			geom_point(size = 3) + 
			scale_x_continuous(
				name = 'Peak jumping velocity',
				limits = c(-1.0, 0.7),
				breaks = seq(-1,0.5,0.5)
			) +
			scale_y_continuous(
				name = 'Peak swimming velocity',
				limits = c(-1.2, 0.7),
				breaks = seq(-1,0.5,0.5)		
			) +
			frog.ecomorph.plot.extras()

		
	# Fig. 4b: Swimming velocity vs. maximum clinging angle
		ggplot(pp, aes(cling.res, swvel.res, fill = ml, shape = type)) +
			geom_point(size = 3) + 
			scale_x_continuous(
				name = 'Maximum clinging angle',
				limits = c(-0.9, 0.8),
				breaks = c(-0.5,0,0.5)
			) +
			scale_y_continuous(
				name = 'Peak swimming velocity',
				limits = c(-1.2, 0.7),
				breaks = seq(-1,0.5,0.5)		
			) +
			frog.ecomorph.plot.extras()

	# Figures 4c and 4d need phylogenetic PCA for data. PCA done with complete morphological dataset (217 species), 
		# but plotted for 80 species to correspond to performance and key analyses.
		mpca <- phyl.pca(mtree, as.matrix(mdat[,mvars]), mode='cov') # Covariance, since all traits are residuals
		tmp <- diag(mpca$Eval)
		eval.out <- rbind(tmp, tmp/sum(tmp), cumsum(tmp/sum(tmp)))
		rownames(eval.out) <- c('raw.eigenvals','prop.variance','cum.prop.var')
		eval.out # Top half of Table S7
		mpca$Evec # Bottom half of Table S7
		# We'll attach the PC scores in 'mpca$S' for morphological plotting
		unique(rownames(mpca$S) == rownames(mdat)) # Only 'TRUE', meaning the row orders match
		tmp <- data.frame(mpca$S[,str_c('PC',1:4)])
		tmp$species <- rownames(tmp)
		pp <- left_join(pp, tmp, by = 'species')
	
	# Fig. 4c: morphology pPC1 by pPC2
		ggplot(pp, aes(PC1, PC2, fill = ml, shape = type)) +
			geom_point(size = 3) + 
			scale_x_continuous(
				name = 'Shape pPC1 (-Webbing, -Toe\\n tip, -Finger tip)',
				limits = c(-3, 4),
				breaks = seq(-3,4,1)
			) +
			scale_y_continuous(
				name = 'Shape pPC2 (+Toe tip,\\n+Finger tip, -Webbing)',
				limits = c(-2,3),
				breaks = seq(-2,3,1)		
			) +
			frog.ecomorph.plot.extras()
	
	# Fig. 4d: morphology pPC3 by pPC4
		ggplot(pp, aes(PC3, PC4, fill = ml, shape = type)) +
			geom_point(size = 3) + 
			scale_x_continuous(
				name = 'Shape pPC3 (-Metatarsal tubercle)',
				limits = c(-1.3,1.3),
				breaks = seq(-1,1,0.5)
			) +
			scale_y_continuous(
				name = 'Shape pPC4 (-Leg muscle)',
				limits = c(-1.05,1),
				breaks = seq(-1,1,0.5),
				expand = c(0,0)		
			) +
			frog.ecomorph.plot.extras()

### C. Fig. 4e: 2-block partial least-squares analysis on all taxa
	
	# Set up variables for analysis
		dat <- pdat[,c(mvars, pvars)] # Remove all variables that won't be used
		nvars <- length(mvars) # User input for highest number of first variables, which in this case are morphology

	# First do phylogenetic transformation
		dat <- phylo.transform(ptree, dat)

	# Now two-block PLS
		alltax.2bpls <- tbpls(dat, nvars, correlation = T) # Basic analysis
		tmp <- data.frame(alltax.2bpls$scores)
		tmp$species <- rownames(tmp)
		pp <- full_join(pp, tmp, by = 'species')

	# Plot results for dimension 1 of 2B-PLS
		ggplot(pp, aes(set1.1, set2.1, fill = ml, shape = type)) +
			geom_point(size = 3) + 
			scale_x_continuous(
				name = 'Body shape dimension 1',
				limits = c(-0.17,0.16),
				breaks = round(seq(-0.15,0.15,0.05),2)
			) +
			scale_y_continuous(
				name = 'Performance dimension 1',
				limits = c(-0.15,0.12),
				breaks = round(seq(-0.15,0.10,0.05),2)
			) +
			frog.ecomorph.plot.extras()
	
	### D. 2-block partial least-squares analysis on Mantellidae alone and all other frogs alone (Tables S1, S2)
	
	# Set up data 	
		cont.vars <- c(mvars, pvars)
		dat <- pdat[,c('family',cont.vars)] # Remove all variables that won't be used
		nvars <- length(mvars)
		dat[,cont.vars] <- phylo.transform(ptree, dat[,cont.vars]) # Note: quantitative results
			# in 2B-PLS are nearly identical if taking taxon subsets prior to phylogenetic transform,
			# rather than doing the phylogenetic transform before taking subsets, as done here.

	# Make two datasets for two analyses (Mantellidae and non-Mantellidae)
		sp <- dat$family == 'Mantellidae' 
		data.list <- lapply(list(sp,!sp), function(x) dat[x,cont.vars]) 
		names(data.list) <- c('Mantellidae', 'nonMantellidae')

	# 2B-PLS analysis (most analysis)
		results <- lapply(data.list, tbpls, nvar = nvars, give.scores = F)
		results[['Mantellidae']]$two.block # Table S1, top rows
		results[['Mantellidae']]$axes.correlation # Table S1, penultimate row
		results[['nonMantellidae']]$two.block # Table S2, top rows
		results[['nonMantellidae']]$axes.correlation # Table S2, penultimate row
	
	# 2B-PLS permutation test of axis significance
		results.permute <- lapply(data.list, choose.axes.2BPLS, nvar = nvars)
		results.permute[['Mantellidae']] # Table S1, bottom row
		results.permute[['nonMantellidae']] # Table S2, bottom row

	# Vector correlation between first two (significant) axes in the Mantellidae and non-Mantellidae 2B-PLS
		# Merge matrices
			mm <- data.frame(rbind(results[['Mantellidae']]$two.block$u, results[['Mantellidae']]$two.block$v))
			nm <- data.frame(rbind(results[['nonMantellidae']]$two.block$u, results[['nonMantellidae']]$two.block$v))
			nm[,3] <- -1 * nm[,3] # Rigid rotation of this axis to make it comparable to mantellid 
			#	results (based on largest coefficients in non-Mantellid results: head length, 
			#	foot webbing, jump velocity, and swim velocity; see Table S2 legend for more justification)
	
		# Vector correlation function (following Collyer and Adams 2007 Ecology)
			vec.cor <- function(v1, v2){
				v1 <- as.numeric(matrix(v1,byrow=T)) # Force into column vector
				v2 <- as.numeric(matrix(v2,byrow=T)) # Force into column vector
				d1 <- sqrt(sum(v1^2)) # Length of vector
				d2 <- sqrt(sum(v2^2)) # Length of vector
				res <- t(v1/d1) %*% (v2/d2) # Vector correlation
				return(res)
			}

		# Calculate vector correlation between axes of each of the two datasets
			unlist(map2(mm, nm, vec.cor))
		
	### E. Clean workspace
		rm(mm, nm, vec.cor, results.permute, results, sp, data.list, colors, shapes, nvars, 
			tmp, pp, mpca, eval.out, dat, alltax.2bpls, frog.ecomorph.plot.extras)
		
###################################################
###### Part 2: Rates of phenotypic evolution ######
###################################################

### A. 80-species tree: analyses and power simulations to produce Figs. 5 and S1
	# Make Mantellidae factor for comparison with power simulation function ('candidate.adap.rad')
		mants <- rownames(pdat)[pdat$family == 'Mantellidae']
	# Set grouping factor for direct analysis of rates ('compare.evol.rates')
		gp <- rep('other', nrow(pdat))
		gp[pdat$family == 'Mantellidae'] <- 'focal'
		names(gp) <- rownames(pdat)
	# Analyze body size 
		# Simple analysis of just rates
			size.dat <- pdat[['SVL']] # Single traits must be input to 'compare.evol.rates' as a vector
			names(size.dat) <- rownames(pdat)
			size.cer <- compare.evol.rates(size.dat, ptree, gp, method = 'simulation')
			size.cer # Print most results presented in Table 2
			size.cer$sigma.d.all # Print overall rate for whole tree (first column of Table 2)			
		# Power analysis to produce Fig. 5a and part of Fig. S1. Manuscript results are 
			# based on 1000 simulations; 100 are done here to reduce computation time.
			size.out <- candidate.adap.rad(ptree, mants, size.dat, 
				rr = seq(1.1, 3.5, by = 0.1), sims = 100, power.verbose = T, snowfall.nodes = 2)
			size.out$results.summary # Overall results, including predictions based on power simulations	
		# Plotting function for power curves
			power.plot.fxn <- function(dat){
				ggplot(dat, aes(rate.ratio, power)) +
					geom_line()	+
					geom_point(fill = 'white', shape = 21, size = 2) +
					theme_classic() +
					scale_x_continuous(
						name = 'Rate ratio',
						breaks = seq(1.2,2.4,0.2)
					) +
					scale_y_continuous(
						name = 'Power',
						limits = c(0, 1),
						breaks = seq(0,1,0.2)		
					) +
					theme(
						axis.text.y = element_text(angle = 90, hjust = 0.5),
						axis.text = element_text(color = 'black')
					)			
			}
		# Plot power curve as in leftmost panel of Fig. S1 (80-taxon white curve)
			size.out$power.curve %>%
				filter(rate.ratio <= 2.5) %>% # Removing values outside range of Fig. S1
				power.plot.fxn()		
					
	# Analyze shape
		# Simple analysis of just rates
			shape.cer <- compare.evol.rates(pdat[,mvars], ptree, gp, method = 'simulation')
			shape.cer # Print most results presented in Table 2
			shape.cer$sigma.d.all # Print overall rate for whole tree (first column of Table 2)
		# Power analysis to produce Fig. 5b and part of Fig. S1. Manuscript results are 
			# based on 1000 simulations; 100 are done here to reduce computation time.
			shape.out <- candidate.adap.rad(ptree, mants, pdat[,mvars], 
				rr = seq(1.1, 2.5, by = 0.1), sims = 100, power.verbose = T, snowfall.nodes = 2)
			shape.out$results.summary # Overall results, including predictions based on power simulations
		# Plot power curve as in rightmost panel of Fig. S1 (80-taxon white curve)
			power.plot.fxn(shape.out$power.curve)

	# Analyze performance
		# Simple analysis of just rates
			perf.cer <- compare.evol.rates(pdat[,pvars], ptree, gp, method = 'simulation')
			perf.cer # Print most results presented in Table 2 (rate ratio in Table 2 is reciprocal: Mantellidae sigma / Other sigma)
			perf.cer$sigma.d.all # Print overall rate for whole tree (first column of Table 2)
		# Power analysis to produce Fig. 5b and part of Fig. S1. Manuscript results are 
			# based on 1000 simulations; 100 are done here to reduce computation time.
			perf.out <- candidate.adap.rad(ptree, mants, pdat[,pvars], 
				rr = seq(1.1, 2.5, by = 0.1), sims = 100, power.verbose = T, snowfall.nodes = 2)
			perf.out$results.summary # Overall results, including predictions based on power simulations
		# Plot power curve as in middle panel of Fig. S1 (80-taxon white curve)
			power.plot.fxn(perf.out$power.curve)

	# Plot bar graphs for Fig. 4a-c (4d not included here because it was based on a separate set of analyses; see section D below)
		# Extract information from previous analyses
			keep <- c('out.pred','AR.pred','out.obs','AR.obs')
			rates.plot <- tibble(
				trait = c(
					rep('Body size 80 sp', 4),
					rep('Body shape 80 sp',4),
					rep('Performance 80 sp', 4)
				),
				type = rep(c('Predicted\\nother','Predicted\\nMantellidae','Observed\\nother','Observed\\nMantellidae'), 3),
				rate = c(
					size.out$results.summary[keep],
					shape.out$results.summary[keep],
					perf.out$results.summary[keep]
				)
			) %>%
			mutate(
				type = factor(type, levels = c(
					'Predicted\\nother','Predicted\\nMantellidae','Observed\\nother','Observed\\nMantellidae')),
				trait = factor(trait, levels = c('Body size 80 sp','Body shape 80 sp','Performance 80 sp'))
			) 
		# Plot
			ggplot(rates.plot, aes(type, rate, fill = type)) +
				geom_bar(stat = 'identity', color = 'black') + 
				facet_wrap(~trait, scale = 'free_y', nrow = 2) + 
				scale_fill_grey(start = 1, end = 0) + 
				theme_classic() + 
				theme(
					legend.position = 'none',
					axis.title = element_text(size = 14),
					axis.text.x = element_text(size = 7, color = 'black'),
					axis.ticks.y = element_line(color = 'black')
					) + 
				xlab(NULL) + 
				ylab('Rate')
	
	# Clean workspace
		rm(rates.plot, keep, power.plot.fxn, perf.out, shape.out, size.out, perf.cer, 
			shape.cer, size.cer, size.dat, gp, mants)
	
### B. 217-species tree: analyses to produce Table 2 results

	# Set grouping factor for direct analysis of rates ('compare.evol.rates')
		gp <- rep('other', nrow(mdat))
		gp[mdat$family == 'Mantellidae'] <- 'focal'
		names(gp) <- rownames(mdat)
	# Analyze body size 
		size.dat <- mdat[['SVL']] # Single traits must be input to 'compare.evol.rates' as a vector
		names(size.dat) <- rownames(mdat)
		size.cer <- compare.evol.rates(size.dat, mtree, gp, method = 'simulation')
		size.cer # Print most results presented in Table 2
		size.cer$sigma.d.all # Print overall rate for whole tree (first column of Table 2)			
	# Analyze shape
		shape.cer <- compare.evol.rates(mdat[,mvars], mtree, gp, method = 'simulation')
		shape.cer # Print most results presented in Table 2
		shape.cer$sigma.d.all # Print overall rate for whole tree (first column of Table 2)
	# Clean workspace
		rm(gp, size.dat, size.cer, shape.cer)
		
### C. 69-species tree (no Pelodryadinae): analyses to produce Table 2 results

	# Remove Pelodryadinae
		pelo <- pdat$family == 'Hylidae' & pdat$location == 'Australia'
		keep.tax <- rownames(pdat)[!pelo]
		pdat.69sp <- pdat[keep.tax,]
		ptree.69sp <- keep.tip(ptree, keep.tax)
	# Set grouping factor for direct analysis of rates ('compare.evol.rates')
		gp <- rep('other', nrow(pdat.69sp))
		gp[pdat.69sp$family == 'Mantellidae'] <- 'focal'
		names(gp) <- rownames(pdat.69sp)
	# Analyze body size 
		size.dat <- pdat.69sp[['SVL']] # Single traits must be input to 'compare.evol.rates' as a vector
		names(size.dat) <- rownames(pdat.69sp)
		size.cer <- compare.evol.rates(size.dat, ptree.69sp, gp, method = 'simulation')
		size.cer # Print most results presented in Table 2
		size.cer$sigma.d.all # Print overall rate for whole tree (first column of Table 2)			
	# Analyze shape
		shape.cer <- compare.evol.rates(pdat.69sp[,mvars], ptree.69sp, gp, method = 'simulation')
		shape.cer # Print most results presented in Table 2
		shape.cer$sigma.d.all # Print overall rate for whole tree (first column of Table 2)
	# Analyze performance
		perf.cer <- compare.evol.rates(pdat.69sp[,pvars], ptree.69sp, gp, method = 'simulation')
		perf.cer # Print most results presented in Table 2
		perf.cer$sigma.d.all # Print overall rate for whole tree (first column of Table 2)
	# Clean workspace
		rm(pelo, keep.tax, pdat.69sp, ptree.69sp, gp, size.dat, size.cer, shape.cer, perf.cer)

### D. 206-species tree (no Pelodryadinae): analyses to produce Table 2 results
	# Remove Pelodryadinae
		pelo <- mdat$family == 'Hylidae' & mdat$location == 'Australia'
		keep.tax <- rownames(mdat)[!pelo]
		mdat.206sp <- mdat[keep.tax,]
		mtree.206sp <- keep.tip(mtree, keep.tax)
	# Set grouping factor for direct analysis of rates ('compare.evol.rates')
		gp <- rep('other', nrow(mdat.206sp))
		gp[mdat.206sp$family == 'Mantellidae'] <- 'focal'
		names(gp) <- rownames(mdat.206sp)
	# Analyze body size 
		size.dat <- mdat.206sp[['SVL']] # Single traits must be input to 'compare.evol.rates' as a vector
		names(size.dat) <- rownames(mdat.206sp)
		size.cer <- compare.evol.rates(size.dat, mtree.206sp, gp, method = 'simulation')
		size.cer # Print most results presented in Table 2
		size.cer$sigma.d.all # Print overall rate for whole tree (first column of Table 2)			
	# Analyze shape
		shape.cer <- compare.evol.rates(mdat.206sp[,mvars], mtree.206sp, gp, method = 'simulation')
		shape.cer # Print most results presented in Table 2
		shape.cer$sigma.d.all # Print overall rate for whole tree (first column of Table 2)
	# Note that Fig. 5d comes from a full power analysis of the shape data on the 206-species tree. 
	#	I did not replicate it here because it takes a long time to run, but here is how you would do it:
		# Make Mantellidae factor for comparison with power simulation function ('candidate.adap.rad')
			#	mants <- rownames(mdat.206sp)[mdat.206sp$family == 'Mantellidae']
		# Do power simulation
			#	shape.out <- candidate.adap.rad(mtree.206sp, mants, mdat.206sp[,mvars], 
			#		rr = seq(1.1, 2.5, by = 0.1), sims = 100, power.verbose = T, snowfall.nodes = 2)
			#	shape.out$results.summary # Overall results, including predictions based on power simulations
		# To add to bar graph (Fig. 5d), just add these results in the same way as other results 
		#	were added to the data frame 'rates.plot'

	# Clean workspace
		rm(pelo, keep.tax, mdat.206sp, mtree.206sp, gp, size.dat, size.cer, shape.cer,
			pdat, mdat, ptree, mtree, mvars, pvars, cont.vars)

###################################################
###### Part 3: Diversification-rate analyses ######
###################################################

# A. Load diversification-specific data
	dd <- read.csv("SuppInfoS9.diversification.data.csv", row.names = 1)
	dtree <- read.tree("SuppInfoS8.anuran.familylevel.tre")

# B. Set up comparison factor for phyANOVA 
	dd$gp <- rep('other', nrow(dd))
	dd$gp[rownames(dd) == 'Mantellidae'] <- 'focal'
	dd$gp <- factor(dd$gp)

# C. Do phylogenetic ANOVA across all sets of diversification rates using method-of-moments estimators, excluding rates that are NA or 0
	dvars <- unlist(lapply(c('AW.','Perl.'), function(x) grep(x, colnames(dd), value = T)))
	set.seed(134) # To ensure results here match those in the MS, if (and only if) using condition 'nsim = 1000'	
	results <- unlist(lapply(dvars, loop.ANOVA.pval, dat = dd, phy = dtree, sims = 100)) # Done at 'nsim = 1000' for paper; here reduced to speed computation
	names(results) <- dvars
	results # P-values for phylogenetic ANOVAs (probability that Mantellidae has a significantly higher diversification rate)
	
# D. Phylogenetic ANOVA for birth-death-based net diversification rates (done later, hence separate code and setting separate random-number seed)
	set.seed(1952) # To ensure reproducibility of results in MS
	loop.ANOVA.pval(dd, dtree, 'yule_rates')
	set.seed(1952) # To ensure reproducibility of results in MS
	loop.ANOVA.pval(dd, dtree, 'bd_rates')
	
# E. Graphically examine why P-values are generally large
	graph.dat <- tibble(dd) %>%
		mutate(family = rownames(dd)) %>%
		select(family, AW.JP.st.0.5) %>%
		rename(stem.rates = AW.JP.st.0.5) %>%
		filter(!is.na(stem.rates))
	# Number of decimal places to round to
		decimalplaces <- function(x) {
		   if ((x %% 1) != 0) {
				nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
			} else {
				return(0)
			}
		}
		interval <- 0.01 # Bar width on histogram below (this code was produced iteratively)
		ratebreaks <- seq(0,0.08,interval) # Set breaks for histogram, based on one unit above the largest rate
	# Line for Mantellidae
		mrate <- graph.dat %>% filter(family == 'Mantellidae') %>% pull(stem.rates)
		mant.line <- tibble(rate = rep(mrate, 2), height = c(9,11.25))
		n <- decimalplaces(interval)
		mant.line[['rate']] <- round(mant.line[['rate']], digits = n)
		ldat <- tibble(mant.line[1,1], mant.line[2,2] + 0.75, label = 'Mantellidae')
	ggplot(graph.dat, aes(stem.rates)) +
		geom_histogram(binwidth = interval) +
		geom_line(aes(rate, height), 
			data = mant.line, 
			arrow = arrow(length = unit(0.10, 'inches'), ends = 'first', type = 'closed'),
			size = 1.25
		) +
		geom_text(aes(rate, height, label = label), data = ldat, size = 5) +
		scale_y_continuous(expand = c(0,0), breaks = seq(0,14,2)) +
		scale_x_continuous(expand = c(0,0), breaks = ratebreaks) +
		xlab('Diversification rates (stem)') + 
		ylab('Number of families') +
		theme_classic() +
		theme(
			axis.text = element_text(size = 11, color = 'black'),
			axis.ticks.y = element_line(color = 'black'),
			axis.title = element_text(size = 12)
		)
	
# F. Clean workspace
	rm(dd, dvars, results, graph.dat, decimalplaces, interval, ratebreaks, mrate, mant.line, n, ldat)

# G. Diversification rates of Pelodryadinae (results presented in Discussion) 
	# Reload to start fresh
		dd <- read.csv("SuppInfoS9.diversification.data.csv", row.names = 1) 
	# Remove Hylidae, of which Pelodryadinae is a part. Also keep only variables analyzed in 
		# main MS (results qualitatively identical for other diversification rates).
		dd <- dd[rownames(dd) != 'Hylidae', c('n.AW', 'JP.stem.age', 'AW.JP.st.0', 'AW.JP.st.0.5', 'AW.JP.st.0.9')]
	# Add data for Pelodryadinae
		dd[nrow(dd) + 1,] <- c(213, 66.280032, 0.0808885, 0.07050131, 0.04677264)
		rownames(dd)[nrow(dd)] <- 'Pelodryadinae'
	# Add comparison factor
		dd$gp <- 'other'
		dd$gp[rownames(dd) == 'Pelodryadinae'] <- 'focal'
	# Adjust tree for these analyses
		dtree$tip.label <- gsub('Hylidae', 'Pelodryadinae', dtree$tip.label)
	# Do phylogenetic ANOVA across the three diversification rates
		dvars <- grep('AW.JP.st.', colnames(dd), value = T)
		results <- vector('double',length(dvars))
		names(results) <- dvars
		set.seed(134) # To ensure results here match those in the MS, if using condition 'nsim = 1000'	
		results <- unlist(lapply(dvars, loop.ANOVA.pval, dat = dd, phy = dtree, sims = 500)) # Done at 'nsim = 1000' for paper; here reduced to speed computation
		names(results) <- dvars
		results # P-values for phylogenetic ANOVAs (probability that Pelodryadinae has a significantly higher diversification rate)
	# Clean work space 
		rm(results, dvars, dtree, dd)
		
################################################################################
###### Part 4: Testing Brownian motion as a model of phenotypic evolution ######
################################################################################

### A. Remove rounding errors in branch lengths that cause error messages in subsequent tests
	mtree <- force.ultrametric(mtree, method = 'nnls')	
	ptree <- force.ultrametric(ptree, method = 'nnls')	
		
### B. Put all data into a list for looping across groups, trees, and traits
	# Looping across all six datasets (80- vs. 217-species trees, all vs. Mantellidae vs. others)
		dn <- paste0(rep(c('all', 'mant','other'), 2), rep(c('_80sp','_217sp'), each = 3))
		make_dataset <- function(x){
			if(grepl('80sp', x)){
			# Start with tree
				tree <- ptree 
				if(grepl('mant', x)) tree <- extract.clade(tree, 119)
				if(grepl('other', x)) tree <- drop.tip(tree, extract.clade(tree, 119)$tip.label)
			# Then the rest takes care of itself
				dat <- pdat[tree$tip.label,]
				traits <- c('SVL', grep('.res', colnames(dat), value = T))					
			} else {
				tree <- mtree
				if(grepl('mant', x)) tree <- extract.clade(tree, 360)
				if(grepl('other', x)) tree <- drop.tip(tree, extract.clade(tree, 360)$tip.label)
				dat <- mdat[tree$tip.label,]
				traits <- c('SVL', grep('.res', colnames(dat), value = T))					
			}
			full_dat <- list(tree = tree, dat = dat, traits = traits)
			return(full_dat)
		}
		all_data <- lapply(dn, make_dataset) # Make six datasets for looping Geiger
		names(all_data) <- dn		
	# Check to ensure data set up properly
		map(all_data, ~dim(.$dat)) # Looks good
		map(all_data, ~.$tree) # Also looks good

### C. Use 'phylosig' from Phytools to estimate lambda and confidence region (2 lnL units above and below MLE)
	# Function for looping analyses
		phytools_lambda <- function(tree, dat, trait){
			x <- setNames(dat[tree$tip.label, trait], tree$tip.label)
			out <- phylosig(tree, x, method = 'lambda', test = TRUE)
			p_val <- 1 - pchisq(2 * (out$logL - out$lik(1)), 1)
			vals <- seq(0, 1, 0.001)
			lik <- map_dbl(vals, out$lik)
			CI <- range(vals[lik >= (max(lik) - 2)]) # Extents of lambda values within 2 lnL units of MLE
			out_vec <- setNames(round(c(out$lambda, CI, p_val), 4), c('MLE','lowCI','highCI','p_from_1.0'))
			return(out_vec)
		}
	# Run for all traits (total of 64 calls to 'phytools_lambda()')
		lambda_out <- lapply(dn, function(z){
			x <- all_data[[z]]
			out <- sapply(x$traits, function(y) phytools_lambda(x$tree, x$dat, y))
			return(out)
		})
		names(lambda_out) <- dn
	# Save separate results
		save(lambda_out, file = paste0('all_lambda_results_',Sys.Date(),'.Rdata'))
		for (i in dn){
			x <- lambda_out[[i]]
			write.csv(x, file = paste0(i,'_phytools_lambda_',Sys.Date(),'.csv'))
		}
	# Now flip for easier compilation
		reformat_results <- function(d, dat = lambda_out){
			x <- t(dat[[d]]) # Transpose
			rownames(x) <- gsub('.res', '', rownames(x)) # Simplify
			colnames(x) <- paste0(colnames(x),'_',d)
			x <- as.data.frame(x) 
			x$trait <- rownames(x) 
			x <- as_tibble(x[,c(5,1:4)])
			return(x)
		}			
		new_out <- lapply(dn, reformat_results)			
		all_out <- full_join(new_out[[1]], new_out[[2]])
		for (i in 3:length(dn)) all_out <- full_join(all_out, new_out[[i]])		
		write_csv(all_out, path = paste0('compiled_lambda_MLE_95CI_',Sys.Date(),'.csv'))
	# Make Table S8 manually, using regex in BBedit

### D. Test four models of evolution, a la Harmon et al. 2010 Evolution
	# Function to loop Geiger's 'fit_continuous()' over multiple models per trait and dataset
		geiger_test <- function(tree, dat, trait, mod_names = c('BM','OU','EB','lambda')){
			x <- setNames(dat[tree$tip.label, trait], tree$tip.label)
			mods <- map(mod_names, ~fitContinuous(phy = tree, dat = x, model = .))
			lnL <- map_dbl(mods, ~.$opt$lnL)
			AICc <- map_dbl(mods, ~.$opt$aicc)
			aic.weight <- function(x, rd.digits = 3){
				l <- min(x) # Calculate the lowest AIC score
				dif <- x - l # Get differences
				d <- exp(-0.5*(dif)) # Exponentiate
				if(is.null(rd.digits)){
					wg <- d / sum(d)
				} else {
					wg <- round(d / sum(d), digits=rd.digits)
				}
				return(wg)
			}
			wgt <- aic.weight(AICc)
			out_mat <- tibble(mod_names, lnL, AICc, wgt)
			names(out_mat)[2:4] <- paste0(trait,'_',names(out_mat)[2:4])
			return(out_mat)
		}
	# Loop Geiger's 'fitContinuous()'. Three levels: first is one of 6 datasets via 'loop_fitCont()', 
	#	second is 10–13 traits within dataset via 'geiger_test()', and within 'geiger_test()' one loops over 4 models.
	#	Total number of 'fitContinuous()' runs = 276
		# Primary input is 'd', which corresponds to the different dataset names in 'dn'
		# Output is a list of model results (all, table of AICc, table of lnL, and table of AICc weights) for each dataset
		loop_fitCont <- function(d, dat = all_data){
			x <- dat[[d]]
			tmp <- lapply(x$traits, function(y) geiger_test(x$tree, x$dat, y))
			tmp_cat <- full_join(tmp[[1]], tmp[[2]])
			for(i in 3:length(tmp)) tmp_cat <- full_join(tmp_cat, tmp[[i]])
			kick_out <- list(complete = tmp_cat)
			for (i in c('lnL','AICc','wgt')){
				kick_out[[i]] <- kick_out[['complete']] %>%
					select(mod_names, ends_with(i))
			}
			return(kick_out)
		}
		# Run:
			mod_results <- lapply(dn, loop_fitCont)
			names(mod_results) <- dn
		# Save separate results
			save(mod_results, file = paste0('all_BM_results_',Sys.Date(),'.Rdata'))
			for (i in dn){
				x <- mod_results[[i]]
				for (j in c('lnL','AICc','wgt')){
					write_csv(x[[j]], path = paste0(i,'_BMmodtest_',j,'_',Sys.Date(),'.csv'))
				}
			}	
	# Now put all of the AICc weights together to make Table S9
		filenames <- dir(pattern = 'wgt')
		nm <- map_chr(filenames, ~paste(unlist(strsplit(., split = '_'))[1:2], collapse = '_'))
		wgt <- list()
		for (i in filenames) wgt[[i]] <- read_csv(i)
		names(wgt) <- nm
		for (i in nm){
			x <- wgt[[i]] %>% 
				mutate(dataset = i) %>%
				select(dataset, mod_names, everything())
			names(x) <- names(x) %>%
				str_replace('_wgt', '') %>%
				str_replace('.res', '')
			wgt[[i]] <- x
			if(i == nm[1]) all_wgt <- x else all_wgt <- full_join(all_wgt, x)
		}
		# Arrange by dataset, then finish table manually
		all_wgt <- all_wgt %>%
			separate(dataset, c('taxa','tree'), '_') %>%
			arrange(tree, taxa) %>%
			select(tree, everything())
		write_csv(all_wgt, path = paste0('AICc_weights_compiled_',Sys.Date(),'.csv'))						
