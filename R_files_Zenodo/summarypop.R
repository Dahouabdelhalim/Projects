library(plyr) # to compute summary data

summarypop <- function(file) { 
	# A few helper functions
	ddtapply <- function(x, FUN=mean, ...) {
		tt <- tapply(x, list(ddind$Gener, ddind$Line), FUN, ...)
		tt["P",] <- tt["P","Parent"] # The Parent population needs to be duplicated for all selected lines
		c(tt[c("P",paste0("F", 1:4)), c("Up","Down","Ran")])
	}
	ddtapply.cov <- function(x1, x2) { # Specific code for the covariance between two columns
		tt <- by(cbind(x1, x2), list(ddind$Gener, ddind$Line), function(xx) var(xx, na.rm=TRUE)[1,2])
		# the rest is identical to the general function
		tt["P",] <- tt["P","Parent"] # The Parent population needs to be duplicated for all selected lines
		c(tt[c("P",paste0("F", 1:4)), c("Up","Down","Ran")])
	}
	ddtapply.sel <- function(x, FUN=mean, ...) {
		x.up <- x.do <- x
		x.up[!grepl(as.character(ddind$SEL), pattern="^Ut")] <- NA
		x.do[!grepl(as.character(ddind$SEL), pattern="^Dt")] <- NA
		tt.up <- ddtapply(x.up, FUN, ...)
		tt.do <- ddtapply(x.do, FUN, ...)
		tt <- ddtapply(x, FUN, ...)
		c(tt.up[1:5], tt.do[6:10], tt[11:15])
	}
	
	dd.raw <- read.table(file, header=TRUE)
	# from Elena, turnng the raw data.frame into summary data.frame
	dd <- ddply(dd.raw, .(ind, Line, Gener, GenLine, dam, sire, SEL), summarize, logGA=mean(logGA), logBA=mean(logUBA))
	
	
	# individual means (to be removed when the data file will be cleaned)
	ddind <-  do.call(rbind, by(dd, as.factor(dd$ind), function(xx) as.data.frame(lapply(xx, function(xxx) if (is.numeric(xxx)) mean(xxx, na.rm=TRUE) else xxx[1]))))
	
	ans <- data.frame(
		'Rep'=rep(c("Up","Down","Control"), each=5), 
		'Gen'=rep(1:5, 3), 
		'mu.x'= ddtapply(ddind$logGA, FUN=mean, na.rm=TRUE),
		'sig2.x'= ddtapply(ddind$logGA, FUN=var, na.rm=TRUE),
		'mus.x' = ddtapply.sel(ddind$logGA, FUN=mean, na.rm=TRUE),
		'sig2s.x'=ddtapply.sel(ddind$logGA, FUN=var, na.rm=TRUE),
		'mu.y'= ddtapply(ddind$logBA, FUN=mean, na.rm=TRUE),
		'sig2.y'= ddtapply(ddind$logBA, FUN=var, na.rm=TRUE),
		'sig.xy'= ddtapply.cov(ddind$logGA, ddind$logBA),
		'N'=ddtapply(ddind$logGA, FUN=function(x) sum(!is.na(x))))
	ans$se.x <- sqrt(ans$sig2.x/ans$N)
	ans$se.y <- sqrt(ans$sig2.y/ans$N)
	ans$S <- ans$mus.x-ans$mu.x
	ans$deltaV <- ans$sig2s.x-ans$sig2.x
	
	ans
}
