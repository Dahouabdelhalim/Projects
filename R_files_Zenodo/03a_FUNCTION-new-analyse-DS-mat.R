analyse.DS.mat <- function(r_i = rel.impacts, studied.drivers = study.drivers, IR = indicator.reweights, summarise.evidence=FALSE, show.result=FALSE){
	#Although this could be done using more functions from EloRating, the DS.mat is instead constructed manually
	#The reason is that it makes weighting and different treatment of draws easier
	
	require(EloRating)

  names(r_i) <- make.names(names(r_i)) #Some of the original names are not syntactically valid

	# replace inferred relative impacts with NAs, using info from studied.drivers
	not.CC.studies <- studied.drivers$Study[studied.drivers$CC == FALSE]
	r_i$Direct.driver.importance.rescaled.Climate.change[r_i$UI %in% not.CC.studies] <- NA

	not.DE.studies <- studied.drivers$Study[studied.drivers$DE == FALSE]
	r_i$Direct.driver.importance.rescaled.Direct.exploitation[r_i$UI %in% not.DE.studies] <- NA

	not.IAS.studies <- studied.drivers$Study[studied.drivers$IAS == FALSE]
	r_i$Direct.driver.importance.rescaled.Invasive.alien.species[r_i$UI %in% not.IAS.studies] <- NA

	not.LU.studies <- studied.drivers$Study[studied.drivers$LU == FALSE]
	r_i$Direct.driver.importance.rescaled.Land.sea.use.change[r_i$UI %in% not.LU.studies] <- NA

	not.PO.studies <- studied.drivers$Study[studied.drivers$PO == FALSE]
	r_i$Direct.driver.importance.rescaled.Pollution[r_i$UI %in% not.PO.studies] <- NA

	not.OT.studies <- studied.drivers$Study[studied.drivers$OT == FALSE]
	r_i$Direct.driver.importance.rescaled.Other[r_i$UI %in% not.OT.studies] <- NA

	#Establish structure to hold the scores ready for conversion to pairgames, i.e., one row per driver within each study
	r_i$row.number <- c(1:nrow(r_i))
	r_i$unique.ref <- paste(r_i$UI, r_i$row.number, sep="_") #So 'study' will still be unique in bootstraps

	# How many of the 5 big direct drivers?
	r_i$Nd <- as.numeric(!is.na(r_i$Direct.driver.importance.rescaled.Climate.change)) +
	  as.numeric(!is.na(r_i$Direct.driver.importance.rescaled.Direct.exploitation)) +
	  as.numeric(!is.na(r_i$Direct.driver.importance.rescaled.Invasive.alien.species)) +
	  as.numeric(!is.na(r_i$Direct.driver.importance.rescaled.Land.sea.use.change)) +
	  as.numeric(!is.na(r_i$Direct.driver.importance.rescaled.Pollution))
	
	# Calculate case weights from geographic scale and number of drivers
  r_i$case.weight <- r_i$Scale.weight * (r_i$Nd-1) # case.weight is how much clout this row will have in the DS calculation
  r_i$old.cw <- r_i$Scale.weight * r_i$Ndrivers.weight # case weight in pre-July-2021 version of code; not the Ndriver.weight values are NOT the triangular numbers!
  
  # Set up data frame in which the indicators are rows and features of the evidence about them are the columns
	status <- as.data.frame(table(r_i$Indicator.name))
	names(status) <- c("Indicator.name", "Comparisons")
	status$Scale <- tapply(r_i$Scale.weight, r_i$Indicator.name, FUN="sum") #Scale holds the sum of the Scale.weight values for all comparisons for that indicator
	status$Ndrivers <- tapply(r_i$Nd, r_i$Indicator.name, FUN="sum") #Ndrivers holds the sum of the number of drivers across all the indicator's comparisons
	status$ndw <- tapply(r_i$Ndrivers.weight, r_i$Indicator.name, FUN="mean") #ndw holds the mean of the (old) Ndriver.weight values for the indicator's comparisons
	status$sum.case.wt <- tapply(r_i$case.weight, r_i$Indicator.name, FUN="sum") # sum.case.wt is the summed clout per indicator, without considering any indicator weighting

	status$sumscale.ir <- mean(status$Scale)/status$Scale
	status$comparison.ir <- mean(status$Comparisons)/status$Comparisons
	status$equalising.ir <- mean(status$sum.case.wt)/status$sum.case.wt
	
	if (summarise.evidence==TRUE){
	  print(status)
	  plot(status$Indicator.name, status$Comparisons, ylab="Number of comparisons", xlab="Indicator name", cex.axis=0.45,
	       main = paste("Range of values: ", min(status$Comparisons), "-", max(status$Comparisons),
	                    " (factor of ", round(max(status$Comparisons)/min(status$Comparisons),1), ")"), sep="")
	  plot(status$Indicator.name, status$Scale, ylab="Summed scale weight", xlab="Indicator name", cex.axis=0.45,
	       main = paste("Range of values: ", round(min(status$Scale), 3), "-", round(max(status$Scale), 3),
	                    " (factor of ", round(max(status$Scale)/min(status$Scale),1), ")"), sep="")
	  plot(status$Indicator.name, status$sum.case.wt, ylab="Sum of (scale wt * (Nd-1)", xlab="Indicator name", cex.axis=0.45,
	       main = paste("Range of values: ", round(min(status$sum.case.wt), 3), "-", round(max(status$sum.case.wt), 3),
	                    " (factor of ", round(max(status$sum.case.wt)/min(status$sum.case.wt),1), ")"), sep="")
	}

	# Set indicator weights based on value of IW passed to function, or the default which is "comparisons"
	if (IR == "comparisons") r_i$Indicator.reweight <- status$comparison.ir[match(r_i$Indicator.name, status$Indicator.name)] #Weights inversely proportional to no of comparisons
	if (IR == "none") r_i$Indicator.reweight <- 1 #Scale and Nd still weight comparisons but no mitigation of uneven evidence among indicators
	if (IR == "old_fix") r_i$Indicator.reweight <- r_i$Old.indicator.reweight #Passed from code block 09 - will be the same in all replicates!
	if (IR == "sumscale") r_i$Indicator.reweight <- status$sumscale.ir[match(r_i$Indicator.name, status$Indicator.name)]
	if (IR == "equalising") r_i$Indicator.reweight <- status$equalising.ir[match(r_i$Indicator.name, status$Indicator.name)] # Hopefully equalises contribution to DS calculation

  r_i$effective.weight <- r_i$Indicator.reweight * r_i$case.weight
  EW <- subset(status, select=c(Indicator.name, Comparisons, Scale, Ndrivers))
  EW$effective.weight <- tapply(r_i$effective.weight, r_i$Indicator.name, FUN="sum")
  EW$mean.reweight <- tapply(r_i$Indicator.reweight, r_i$Indicator.name, FUN="mean")
  
  if (show.result == TRUE){
    print(EW)
    plot(EW$Indicator.name, EW$effective.weight, ylab="Effective weight in analysis", xlab="Indicator name", cex.axis=0.45, 
       sub=paste("Analysed using indicator weights based on", IR),
       main = paste("Range of values: ", min(EW$effective.weight), "-", max(EW$effective.weight),
                    " (factor of ", round(max(EW$effective.weight)/min(EW$effective.weight),1), ")"), sep="")
  }
  
	# Set up matrix of zeroes to hold result of the head-to-heads
	DS.mat <- matrix(data=0, nrow=6, ncol=6, dimnames=list(c("CC", "DE", "IAS", "LU", "PO", "OT"), c("CC", "DE", "IAS", "LU", "PO", "OT")))
	
	# Populate the DS matrix using explicit coding
	DS.mat[1,2] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Climate.change > r_i$Direct.driver.importance.rescaled.Direct.exploitation), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Climate.change, r_i$Direct.driver.importance.rescaled.Direct.exploitation, tol=1e-5), na.rm=TRUE)
	DS.mat[2,1] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Climate.change < r_i$Direct.driver.importance.rescaled.Direct.exploitation), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Climate.change, r_i$Direct.driver.importance.rescaled.Direct.exploitation, tol=1e-5), na.rm=TRUE)
	DS.mat[1,3] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Climate.change > r_i$Direct.driver.importance.rescaled.Invasive.alien.species), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Climate.change, r_i$Direct.driver.importance.rescaled.Invasive.alien.species, tol=1e-5), na.rm=TRUE)
	DS.mat[3,1] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Climate.change < r_i$Direct.driver.importance.rescaled.Invasive.alien.species), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Climate.change, r_i$Direct.driver.importance.rescaled.Invasive.alien.species, tol=1e-5), na.rm=TRUE)
	DS.mat[1,4] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Climate.change > r_i$Direct.driver.importance.rescaled.Land.sea.use.change), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Climate.change, r_i$Direct.driver.importance.rescaled.Land.sea.use.change, tol=1e-5), na.rm=TRUE)
	DS.mat[4,1] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Climate.change < r_i$Direct.driver.importance.rescaled.Land.sea.use.change), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Climate.change, r_i$Direct.driver.importance.rescaled.Land.sea.use.change, tol=1e-5), na.rm=TRUE)
	DS.mat[1,5] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Climate.change > r_i$Direct.driver.importance.rescaled.Pollution), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Climate.change, r_i$Direct.driver.importance.rescaled.Pollution, tol=1e-5), na.rm=TRUE)
	DS.mat[5,1] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Climate.change < r_i$Direct.driver.importance.rescaled.Pollution), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Climate.change, r_i$Direct.driver.importance.rescaled.Pollution, tol=1e-5), na.rm=TRUE)
	DS.mat[1,6] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Climate.change > r_i$Direct.driver.importance.rescaled.Other), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Climate.change, r_i$Direct.driver.importance.rescaled.Other, tol=1e-5), na.rm=TRUE)
	DS.mat[6,1] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Climate.change < r_i$Direct.driver.importance.rescaled.Other), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Climate.change, r_i$Direct.driver.importance.rescaled.Other, tol=1e-5), na.rm=TRUE)

	DS.mat[2,3] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Direct.exploitation > r_i$Direct.driver.importance.rescaled.Invasive.alien.species), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Direct.exploitation, r_i$Direct.driver.importance.rescaled.Invasive.alien.species, tol=1e-5), na.rm=TRUE)
	DS.mat[3,2] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Direct.exploitation < r_i$Direct.driver.importance.rescaled.Invasive.alien.species), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Direct.exploitation, r_i$Direct.driver.importance.rescaled.Invasive.alien.species, tol=1e-5), na.rm=TRUE)
	DS.mat[2,4] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Direct.exploitation > r_i$Direct.driver.importance.rescaled.Land.sea.use.change), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Direct.exploitation, r_i$Direct.driver.importance.rescaled.Land.sea.use.change, tol=1e-5), na.rm=TRUE)
	DS.mat[4,2] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Direct.exploitation < r_i$Direct.driver.importance.rescaled.Land.sea.use.change), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Direct.exploitation, r_i$Direct.driver.importance.rescaled.Land.sea.use.change, tol=1e-5), na.rm=TRUE)
	DS.mat[2,5] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Direct.exploitation > r_i$Direct.driver.importance.rescaled.Pollution), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Direct.exploitation, r_i$Direct.driver.importance.rescaled.Pollution, tol=1e-5), na.rm=TRUE)
	DS.mat[5,2] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Direct.exploitation < r_i$Direct.driver.importance.rescaled.Pollution), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Direct.exploitation, r_i$Direct.driver.importance.rescaled.Pollution, tol=1e-5), na.rm=TRUE)
	DS.mat[2,6] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Direct.exploitation > r_i$Direct.driver.importance.rescaled.Other), na.rm=TRUE) +
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Direct.exploitation, r_i$Direct.driver.importance.rescaled.Other, tol=1e-5), na.rm=TRUE)
	DS.mat[6,2] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Direct.exploitation < r_i$Direct.driver.importance.rescaled.Other), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Direct.exploitation, r_i$Direct.driver.importance.rescaled.Other, tol=1e-5), na.rm=TRUE)

	DS.mat[3,4] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Invasive.alien.species > r_i$Direct.driver.importance.rescaled.Land.sea.use.change), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Invasive.alien.species, r_i$Direct.driver.importance.rescaled.Land.sea.use.change, tol=1e-5), na.rm=TRUE)
	DS.mat[4,3] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Invasive.alien.species < r_i$Direct.driver.importance.rescaled.Land.sea.use.change), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Invasive.alien.species, r_i$Direct.driver.importance.rescaled.Land.sea.use.change, tol=1e-5), na.rm=TRUE)
	DS.mat[3,5] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Invasive.alien.species > r_i$Direct.driver.importance.rescaled.Pollution), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Invasive.alien.species, r_i$Direct.driver.importance.rescaled.Pollution, tol=1e-5), na.rm=TRUE)
	DS.mat[5,3] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Invasive.alien.species < r_i$Direct.driver.importance.rescaled.Pollution), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Invasive.alien.species, r_i$Direct.driver.importance.rescaled.Pollution, tol=1e-5), na.rm=TRUE)
	DS.mat[3,6] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Invasive.alien.species > r_i$Direct.driver.importance.rescaled.Other), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Invasive.alien.species, r_i$Direct.driver.importance.rescaled.Other, tol=1e-5), na.rm=TRUE)
	DS.mat[6,3] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Invasive.alien.species < r_i$Direct.driver.importance.rescaled.Other), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Invasive.alien.species, r_i$Direct.driver.importance.rescaled.Other, tol=1e-5), na.rm=TRUE)

	DS.mat[4,5] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Land.sea.use.change > r_i$Direct.driver.importance.rescaled.Pollution), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Land.sea.use.change, r_i$Direct.driver.importance.rescaled.Pollution, tol=1e-5), na.rm=TRUE)
	DS.mat[5,4] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Land.sea.use.change < r_i$Direct.driver.importance.rescaled.Pollution), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Land.sea.use.change, r_i$Direct.driver.importance.rescaled.Pollution, tol=1e-5), na.rm=TRUE)
	DS.mat[4,6] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Land.sea.use.change > r_i$Direct.driver.importance.rescaled.Other), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Land.sea.use.change, r_i$Direct.driver.importance.rescaled.Other, tol=1e-5), na.rm=TRUE)
	DS.mat[6,4] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Land.sea.use.change < r_i$Direct.driver.importance.rescaled.Other), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Land.sea.use.change, r_i$Direct.driver.importance.rescaled.Other, tol=1e-5), na.rm=TRUE)


	DS.mat[5,6] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Pollution > r_i$Direct.driver.importance.rescaled.Other), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Pollution, r_i$Direct.driver.importance.rescaled.Other, tol=1e-5), na.rm=TRUE)
	DS.mat[6,5] <- sum(r_i$Scale.weight * r_i$Indicator.reweight * (r_i$Direct.driver.importance.rescaled.Pollution < r_i$Direct.driver.importance.rescaled.Other), na.rm=TRUE) + 
	  draw.score * sum(r_i$Scale.weight * r_i$Indicator.reweight * near(r_i$Direct.driver.importance.rescaled.Pollution, r_i$Direct.driver.importance.rescaled.Other, tol=1e-5), na.rm=TRUE)

	if (exclude.other.drivers == "Yes") DS.mat <- DS.mat[c(1:5), c(1:5)] 
	david <- DS(DS.mat)
	
	if (show.result == TRUE){
	  print(david)
	  print(DS.mat)
	}
	
	to.return <- list(david_scores = david, david_matrix = DS.mat, evidence=r_i, by.indicator=EW)
	
	#Add useful metadata
	attr(to.return, "analysis") <- "David's score (DS in EloRating package)"
	attr(to.return, "exclude_other_drivers") <- exclude.other.drivers
	attr(to.return, "number_of_rows") <- nrow(r_i)
	attr(to.return, "number_of_different_studies") <- length(unique(r_i$UI))
	attr(to.return, "combine_rank_magnitude") <- combine.rank.magnitude
	attr(to.return, "analysis_level") <- analysis.level
	attr(to.return, "selected_indicators") <- ifelse(indicator.selection == "Yes", toString(selected.indicator), "all")
	attr(to.return, "excluded_indicators") <- ifelse(indicator.exclusion == "Yes", toString(excluded.indicator), "none")
	attr(to.return, "selected_realms") <- ifelse(realm.selection == "Yes", selected.realm, "All realms")
	attr(to.return, "selected_regions") <- ifelse(region.selection == "Yes", selected.region, "All regions")
	attr(to.return, "selected_domains") <- ifelse(domain.selection == "Yes", selected.domain, "All domains")
	attr(to.return, "selected_scales") <- ifelse(scale.selection == "Yes", selected.scale, "All scales")
	attr(to.return, "selected_no_of_drivers") <- ifelse(n.drivers.selection == "Yes", toString(selected.n.drivers), toString(c(2:6)))
	attr(to.return, "manual_selection") <- manual.selection
	attr(to.return, "weighted_mean_procedure") <- weighted.mean.procedure
	attr(to.return, "scale.weights") <- toString(scale.weights)
	attr(to.return, "ndriver.weights") <- toString(ndriver.weights)
	attr(to.return, "indicator.reweights") <- IR
	attr(to.return, "draw.score") <- draw.score

	return(to.return)
}
