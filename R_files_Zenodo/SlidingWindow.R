# Sliding window analysis taken from MEREDITH, R. W., J. E. JANEčKA, J. GATESY, O. A. RYDER, C. A. FISHER, E. C. TEELING, A. GOODBLA, 
# et al. 2011. Impacts of the Cretaceous terrestrial revolution and KPg extinction on mammal diversification. Science 334: 521–524
# and coded in R for Steppan and Schenk "Muroid Rodent Phylogenetics: 900-species Tree Reveals Increasing Diversification Rates" 2017. 
# Coded by J. Schenk.


Sliding.Window <- function(window, phy){
	if (class(phy) != "phylo") stop("\\'phy\\' must be an object of class \\'phylo\\'")
	if (!is.ultrametric(phy)) stop("This function is intended to be used on ultrametric trees only")
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
	plot(Results[ , 1], Results[ , 5], type = "l", xlim = rev(range(Results[ , 1])), xlab = "Million of years before present", ylab = "New lineages / Lineages  at begining of  window")	
	}