#library(MASS)

numSimulation = 1000  # number of simulations to run with a set of parameter combinations

par(mfcol=c(1,2))

# R_corr=seq(from=0.0, to=0.8, by=0.2)	# slope relating Y to X
L2=seq(from=0, to=20, by=2)  # number of units to set lag, ranging from 0 to 20; higher values indicate more lag
numSpeciesPairs=seq(from=20, to=320, by=30)  # number of paired comparisons

output_model2 = matrix(-99, nrow=length(L2) * length(numSpeciesPairs), ncol=4)
colnames(output_model2)=c("L2", "numSpeciesPairs", "percentSignif","meanSlope")
output_model2=data.frame(output_model2)

overallCounter= 1

# for (R_corr_counter in 1:length(R_corr))
# {
	for (L2_counter in 1:length(L2))
	{
		for (numSpeciesPairs_counter in 1:length(numSpeciesPairs))
		{

		# corrMatrix=matrix(c(1, R_corr[R_corr_counter], R_corr[R_corr_counter], 1), nr=2, byrow=TRUE)  
		
		data_Xsp1=NULL  # vectors of simulated traits for unlagged trait X, lagged trait Y, for pairs of species sp1 and sp2		 
		data_Xsp2=NULL
		data_Ysp1=NULL
		data_Ysp2=NULL

		contrast_X=NULL # vectors of unstandardized contrasts for X and Y, between sp1 and sp2
		contrast_Y=NULL

		brLength=NULL  	# vector of branch lengths

		pvalList=NULL   # vector of statistics from test; length=numSimulation
		coefList=NULL
		resPvals=NULL
		lagTimes_list=NULL
		resPvals_list=NULL

		for (n in 1:numSimulation)
		{
			for (i in 1: numSpeciesPairs[numSpeciesPairs_counter])
				{
				brLength[i]=runif(1,0,1)   # branch length taken as a random number from a uniform distribution
	
				N=ceiling(brLength[i]*100)  # branch length segments - integer
	
		 		sampSp1 = rnorm(n=N, mean=0, sd=1)
		 		sampSp2 = rnorm(n=N, mean=0, sd=1)
	
				data_Xsp1[i]=sum(sampSp1)  # for X, just sum the random draws across all segments
				data_Xsp2[i]=sum(sampSp2)
						
				# data_Ysp1[i]=sum(sampSp1[1:L,2])
				# data_Ysp2[i]=sum(sampSp2[1:L,2])

					if (N> L2[L2_counter])   # for Y, have to only sum the ones excluding a lag
					{
					data_Ysp1[i]=sum(sampSp1[1:(N-L2[L2_counter])]) + rnorm(1,0,1)
					data_Ysp2[i]=sum(sampSp2[1:(N-L2[L2_counter])]) + rnorm(1,0,1)
					} 
					
					else  # in this case, species hasn't yet escaped lag period, so draw a random number.  Relevant evolution acquired prior to split, and so subtracted out.					
					{
					data_Ysp1[i]=rnorm(1,0,1)
					data_Ysp2[i]=rnorm(1,0,1)
					}
							
				if (data_Xsp1[i]>data_Xsp2[i])   # find direction of subtraction to make X contrasts positive
					{
					contrast_X[i]=data_Xsp1[i]-data_Xsp2[i]
					contrast_Y[i]=data_Ysp1[i]-data_Ysp2[i]
					}

				else
					{
					contrast_X[i]=data_Xsp2[i]-data_Xsp1[i]
					contrast_Y[i]=data_Ysp2[i]-data_Ysp1[i]
					}
				}
	        
		 	# plot(contrast_X, contrast_Y, cex=(brLength*2), main=paste("L2=",as.character(L2[L2_counter])))
		 	
		 	
		
			regressOut1=lm(contrast_Y ~ contrast_X - 1)  # regress contrasts through the origin
			regressOut2=lm(regressOut1$residuals ~ brLength)  # regress residuals from Out1 on branch length
					
			# plot(brLength, regressOut1$residuals, main=paste("p=", as.character(summary(regressOut2)$coefficients[8])))
					
			coefList[n]=summary(regressOut2)$coefficients[2]

			coefListSign[n]=sign(summary(regressOut2)$coefficients[2])   # store the coefficient's sign
		
			if(coefListSign[n]>0)  # positive slope, so record p-value
				pvalList[n]=summary(regressOut2)$coefficients[8]
			else
				pvalList[n]=99  # set the pvalue to 99 (nonsensical) if the coefficient was negative
			}
					
	#	output_model2$R_corr[overallCounter]=R_corr[R_corr_counter]		
		output_model2$L2[overallCounter]=L2[L2_counter]
		output_model2$numSpeciesPairs[overallCounter]=numSpeciesPairs[numSpeciesPairs_counter]
		output_model2$percentSignif[overallCounter]=sum(pvalList<0.1)/numSimulation # find the p-values less than 0.1, based on one tailed test
		output_model2$meanSlope[overallCounter]=mean(coefList)
		

		# print(sum(pvalList==99))

		overallCounter= overallCounter + 1
	
		} #end loop c
	} #end loop b
# } #end loop a
	


