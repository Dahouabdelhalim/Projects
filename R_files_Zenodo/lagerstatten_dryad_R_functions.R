# function to average over a proxy curve
# x = x variable (time in this case)
# y = y variable
# times = vector of times to integrate over
# outputs average in y over time intervals specified by 'times'
integrate_proxy <- function(x, y, times) {
	f <- approxfun(x=x, y=y, rule=2)
	dt <- diff(times)[1]
	int <- sapply(1:length(times), function(x) {
		t1 <- times[x]
		t2 <- t1 + dt
		try(integrate(f, t1, t2))
	})
	as.numeric(t(int)[, "value"])/dt
}

# function to incorporate fossil uncertainty
# df is a dataset with columns minage_ma and maxage_ma
# function computes what bins these fall into
getcounts <- function(df, times=NULL) {
	bw <- diff(times)[1]
	df$bin1 <- times[findInterval(df$minage_ma, times)]
	df$bin2 <- times[findInterval(df$maxage_ma, times)]
	df$median_age <- apply(df[c("bin1","bin2")], 1, median)
	timeranges <- lapply(1:nrow(df), function(x) {seq(df$bin1[x], df$bin2[x], by=bw)})
	# timeranges <- lapply(1:nrow(df), function(x) {df$bin1[x]:(df$bin2[x]-1)})
	# probs of each lagerstatte being in a given time bin?
	res <- matrix(nrow=length(times), ncol=length(timeranges))
	# find fossils with possibility of lying within given time bin (specified in times)
	for (i in seq_along(times))
		for (j in seq_along(timeranges))
			res[i, j] <- ifelse(times[i] %in% timeranges[[j]], 1/length(timeranges[[j]]), 0)
	rowSums(res)
}
