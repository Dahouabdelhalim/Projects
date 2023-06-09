# sequence over a range
seq_range <- function(x, n=25) {
	x <- range(x, na.rm=TRUE)
	seq(x[1], x[2], length=n)
}