
##Function to perform Randomised Paired T test on two columns of data, where for each permutation the order of of the two datapoints within rows is assigned at random
	##Function accepted 4 arguments - the dataframe where the data are, the two column numbers where the data are held, and the number of permutations to run for

pair.randomise<-function(dframe,col1,col2,nsims=1000){
	
	dataframe<-dframe[,c(col1,col2)]
	
	#Empty vector for randomisations
	simresults<-numeric(nsims)
	
	for (k in 1:nsims){
		newdat<-t(apply(dataframe,1,function(x){(sample(x))}))
		simresults[k]<-t.test(newdat[,1],newdat[,2],paired=T)$statistic
	}
	realt<-t.test(dataframe[,1],dataframe[,2],paired=T)$statistic
	
	#List for data
	x<-list(t=realt,simvals=simresults,pval=ifelse(realt<0,2*mean(realt>simresults),2*mean(realt<simresults)))
	
	if(realt<0){cat("2 tailed p =",2*mean(realt>=simresults))} else {cat("2 tailed p =",2*mean(realt<=simresults))}
	
	cat("  // Number of permutations =", nsims)
	x
}

