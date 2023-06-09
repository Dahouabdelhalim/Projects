#### R code to plot a sliding window and simulated plots.  Simulations are based on the birth-death model. 
#Code generated for Schenk and Steppan "The role of biogeography in adaptive radiation"
#Code generate by J. Schenk 2017
#Code dependent on the Geiger package (https://cran.r-project.org/web/packages/geiger/index.html; HARMON, L. J., J. T. WEIR, C. D. BROCK, R. E. GLOR, AND  W. CHALLENGER. 2008. GEIGER:  Investigating evolutionary radiations. Bioinformatics 24: 129â€“131.)
#User must include a chronogram, window size (2 million years was used in our publication), the number of simulations (default set to 100), speciation rate, extinction rate, the number of tips (n), and the age of the clade (t).  The last four settings are conditions under which the birth-death trees will be simulated.
 

Sliding.Window.Simulations <- function(window, phy, nsims = 100, b = 1, d = 0, n, t){
	if (class(phy) != "phylo") stop("\\'phy\\' must be an object of class \\'phylo\\'")
	if (!is.ultrametric(phy)) stop("This function is intended to be used on ultrametric trees only")
	if (!requireNamespace("geiger", quietly = TRUE)) stop("geiger package is needed for this function to work. Please install it.",
      call. = FALSE)
		#get branching times and sort them
		bt <- sort(branching.times(phy), decreasing = TRUE)
		#finds the age of the root, stores the number for the slider
		rootage <- max(bt)
		#creates the slider, the scale of which the window will progress across the entire length of the tree by increments of 0.02
		MinW <- sort(seq(0 + window, rootage, by = 0.02), decreasing = TRUE)
		n <- length(MinW)
		#make matrix
		Results <- matrix(data = NA, nrow = n, ncol = 5)
		#populate first column with min window values
		Results[ ,1] <- sort(seq(0 + window, rootage, by = 0.02), decreasing = TRUE)
		#populate second column with max window values
		Results[ ,2] <- Results[ ,1] - window
		#Count the number of nodes in the window, put in column 3
		for(i in 1:n)
			Results[i, 3] <- sum(Results[i, 1] > bt & Results[i, 2] < bt)
		#Estimate the number of lineages at the beginning of the window
		for(k in 1:n)
			Results[k, 4] <- sum(Results[k, 1] <= bt)	
		Results[ , 4] <- Results[ , 4] + 1
		#Estimate the number of nodes within a window divided by the number of lineages at the start of the window
		for(j in 1:n)
			Results[j, 5] <- Results[j, 3]/Results[j, 4]
		plot(Results[ , 1], Results[ , 5], type = "l", xlim = rev(range(Results[ , 1])), xlab = "Million of years before present", ylab = "New lineages / Lineages  at begining of  window", las = 1)	
	Sliding.Window.lines <- function(window, phy){
		bt <- sort(branching.times(phy), decreasing = TRUE)
		rootage <- max(bt)
		MinW <- sort(seq(0 + window, rootage, by = 0.02), decreasing = TRUE)
		n <- length(MinW)
		Results <- matrix(data = NA, nrow = n, ncol = 5)
		Results[ ,1] <- sort(seq(0 + window, rootage, by = 0.02), decreasing = TRUE)
		Results[ ,2] <- Results[ ,1] - window
		for(i in 1:n)
			Results[i,3] <- sum(Results[i,1] > bt & Results[i,2] < bt)
		for(k in 1:n)
			Results[k,4] <- sum(Results[k,1] <= bt)	
		Results[ ,4] <- Results[ ,4] + 1
		for(j in 1:n)
			Results[j,5] <- Results[j,3]/Results[j,4]
		lines(Results[ ,1], Results[ ,5], type = "l", col="gray", xlim = rev(range(Results[ ,1])))	
		}
	for(i in 1:nsims){
		sim <- sim.bdtree(b, d, stop = c("taxa", "time"), n = n, t = t, extinct = FALSE)
		dsim <- drop.extinct(sim)
		Sliding.Window.lines(window, dsim)	
		}
	Sliding.Window.Blacklines <- function(window, phy){
		bt <- sort(branching.times(phy), decreasing = TRUE)
		rootage <- max(bt)
		MinW <- sort(seq(0 + window, rootage, by = 0.02), decreasing = TRUE)
		n <- length(MinW)
		Results <- matrix(data = NA, nrow = n, ncol = 5)
		Results[ ,1] <- sort(seq(0 + window, rootage, by = 0.02), decreasing = TRUE)
		Results[ ,2] <- Results[ ,1] - window
		for(i in 1:n)
			Results[i,3] <- sum(Results[i,1] > bt & Results[i,2] < bt)
		for(k in 1:n)
			Results[k,4] <- sum(Results[k,1] <= bt)	
		Results[ ,4] <- Results[ ,4] + 1
		for(j in 1:n)
			Results[j,5] <- Results[j,3]/Results[j,4]
		lines(Results[ ,1], Results[ ,5], type = "l", col="black", xlim = rev(range(Results[ ,1])), lwd = 2)	
		}
	Sliding.Window.Blacklines(window = window, phy = phy)
}





#Example
library(geiger)
mytree <- sim.bdtree(b = 1, d = 0, stop = c("taxa", "time"), n = 100, t = 4, seed = 0, extinct = TRUE)
Sliding.Window.Simulations(2, mytree, nsims = 50, b = 1, d = 0, n, t)

