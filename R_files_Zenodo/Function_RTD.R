## Transform rownames to data
RTD <- function(dat){
	nr <- nrow(dat)
	tt <- length(strsplit(rownames(dat)[1], split=" ")[[1]])
	dat.names <- data.frame(matrix(, nrow=nr, ncol=tt))
	for(i in 1:nr){
		dat.names[i,] <- strsplit(rownames(dat)[i], split=" ")[[1]]
	}
	return(dat.names)
}

## Transform colnames to data
CTD <- function(dat){
	nc <- length(dat)
	tt <- length(strsplit(names(dat)[1], split=" ")[[1]]) 
	dat.names <- data.frame(matrix(, nrow=nc, ncol=tt))
	for(i in 1:nc){
		dat.names[i,] <- strsplit(names(dat)[i], split=" ")[[1]]
	}
	return(dat.names)
}