# Clark and Evans Index of Aggregation with Donnelly correction applied

# References
# Clark, P.J. & Evans, F.C. 1954. Distance to nearest neighbor as a measure of spatial relationships in populations. Ecology 35: 445-453.
# Donnelly, K.P. 1978. Simulations to determine the variance and edge effect of total nearest-neighbor distance. In: Hodder, I. (ed.) Simulation studies in archaeology, pp. 91-95. Cambridge University Press, Cambridge, London.

# For an example application, see: Fangliang, Legendre, LaFrankie (1997) J. Vegetation Science 8: 105 - 114.

# Depends on:
# library("vegan")
# library("spatstat")


ClarkEvans <- function(coords,quadratsize,boundary.strip=c(TRUE,FALSE),b.size=NULL){
	# coords are the observed x,y coordinates for each individual in the quadrat, including the boundary strip if present
	# quadratsize is the length of one side of the quadrat that includes the boundary strip (if present), in same units as coords. 
	# e.g., if have 150mm x 150mm quadrat, then quadratsize = 150. The focal quadrat will be within this.
	# If boundary.strip=FALSE, then the Donnelly (1978) correction is applied
	# b.size is the width of the boundary strip around the focal area.
	# b.size is only used if boundary.strip=TRUE
	# If b.size is not specified, then maximum of the minimum nearest neighbor distance is used
	
	# First, write function to get p-value from z score, following Clark Evans (1954)
	get.p.value<-function(z, one.sided=NULL) {
	    if(is.null(one.sided)) {
  	      pval = pnorm(-abs(z));
  	      pval = 2 * pval
  	  } else if(one.sided=="-") {
  	      pval = pnorm(z);
  	  } else {
  	      pval = pnorm(-z);                                                                                 
  	  }
  	  return(pval);
}

	coords <- as.matrix(coords)
	dis.mat <- as.matrix(vegdist(coords,method="euclidean"))	# needs to be a symetric distance matrix

	diag(dis.mat) <- 9999 			# an unreasonably high number so that the distance from an individual to itself doesn't get condused as a nearest neighbour distance in next step.
	
	NN <- apply(dis.mat, 2, min) 		# calculate the nearest neighbour (NN) distance for each individual (biased)

	if(boundary.strip == FALSE){
		ra <- mean(NN) 										# observed NN distance
		A <- quadratsize*quadratsize						# Surface Area of quadrat
		n <- length(NN) 									# number of NN measurements (number of individuals in study zone)
		p <- n/A 											# density per quadrat
		L <- quadratsize*4 									# length of the boundary of the whole study area
		re <- 1/(2 * sqrt(p)) 								# expected distance to NN NOT corrected (Clark Evans)
		rc <- re + (0.051 + (0.041/sqrt(n)))*(L/n)  		# expected distance to NN corrected for lack of a boundary strip (Clark Evans with Donnelly correction)
		sr <- (sqrt((0.07*A)+((0.037*L)*sqrt(A/n)))) / n	# standard error of the expected distance to nearest neighbor
		z  <- (ra - rc) / sr								# test statistic, standard normal deviate from random (R=1) expectation
		p.val.agg <- get.p.value(z=z,one.sided="-")			# test for significant aggregation.
		p.val.reg <- get.p.value(z=z,one.sided="")			# test for significant regularity

# Use the clarkevens() and clarkevens.test() function from the 'spatstat' package
		coords <- ppp(x=coords[,1],y=coords[,2], xrange=c(0, quadratsize),yrange=c(0, quadratsize))
		R.clarkevans <- clarkevans(coords,correction="Donnelly")
		p.val.agg.clarkevans <- clarkevans.test(coords,correction="Donnelly",alternative="clustered")$p.value
		p.val.reg.clarkevans <- clarkevans.test(coords,correction="Donnelly",alternative="regular")$p.value

		dat <- data.frame(	R					=	ra / rc,
							z.value				=	z,
							p.value.agg			=	p.val.agg,
							p.value.reg			=	p.val.reg,
							min.NN 				=	min(NN),
							max.NN				=	max(NN),
							med.NN				=	median(NN),
							focal.quadrat.Area 	= 	A,
							settler.number 		= 	n,
							settler.density 	= 	p,
							R.clarkevans 		=	R.clarkevans, # from the 'spatstat' package
							p.value.agg.clarkevans =	p.val.agg.clarkevans, # from the 'spatstat' package 	
							p.value.reg.clarkevans =	p.val.reg.clarkevans # from the 'spatstat' package
		)
	}

	if(boundary.strip == TRUE){

	    if(is.null(b.size)) {
			b.size = max(NN)
		} else {
  	      b.size = b.size                                                                                 
  	  }

	# index for individuals that fall within the focal quadrat
		focal.index <- which(coords[,1] >= b.size & coords[,1] <= quadratsize - b.size & coords[,2] >= b.size & coords[,2] <= quadratsize - b.size)

		NNmod <- NN[focal.index]			# only use NN for individuals within the focal area (unbiased)
		ra <- mean(NNmod) 					# unbiased observed NN distance
		A <- (quadratsize-(2*b.size))*(quadratsize-(2*b.size))		# Surface Area of quadrat
		n <- length(NNmod) 					# number of NN measurements (number of individuals in study zone)
		p <- n/A 							# density per quadrat
		re <- 1/(2 * sqrt(p)) 				# expected NN distrance 
		sr <- 0.26136 / sqrt(n*p) 			# standard error of the expected distance to nearest neighbor
		z  <- (ra - re) / sr				# test statistic, standard normal deviate from random (R=1) expectation
		p.val.agg <- get.p.value(z=z,one.sided="-")			# test for significant aggregation.
		p.val.reg <- get.p.value(z=z,one.sided="")			# test for significant regularity

# Use the clarkevens() and clarkevens.test() function from the 'spatstat' package
		coords <- ppp(x=coords[,1],y=coords[,2], xrange=c(0, quadratsize),yrange=c(0, quadratsize))
		R.clarkevans <- clarkevans(coords,correction="Donnelly")
		p.val.agg.clarkevans <- clarkevans.test(coords,correction="Donnelly",alternative="clustered")$p.value
		p.val.reg.clarkevans <- clarkevans.test(coords,correction="Donnelly",alternative="regular")$p.value

		dat <- data.frame(	R					=	ra / rc,
							z.value				=	z,
							p.value.agg			=	p.val.agg,
							p.value.reg			=	p.val.reg,
							min.NN 				=	min(NN),
							max.NN				=	max(NN),
							med.NN				=	median(NN),
							focal.quadrat.Area 	= 	A,
							settler.number 		= 	n,
							settler.density 	= 	p,
							R.clarkevans 		=	R.clarkevans, # from the 'spatstat' package
							p.value.agg.clarkevans =	p.val.agg.clarkevans, # from the 'spatstat' package 	
							p.value.reg.clarkevans =	p.val.reg.clarkevans # from the 'spatstat' package
		)
	}
	return(dat)
}

