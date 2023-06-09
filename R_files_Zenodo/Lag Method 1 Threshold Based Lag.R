# parameters that are fixed in all simulations

numSimulation = 1000  # number of simulations to run with each set of parameters; ideally >100
traitEvolVariance=0.03 # trait variance, held constant assuming that most change occurs immediately after speciation

# parameters that vary
B_slope=seq(from=0.5, to=3, by=0.5)	# slope relating Y to X
L1=seq(from=0, to=0.7, by=0.1)  # number of units to set lag, ranging from 0 to 1; higher values indicate more lag
numSpeciesPairs=seq(from=20, to=320, by=30)  # number of paired comparisons

output = matrix(-99, nrow=length(B_slope) * length(L1) * length(numSpeciesPairs), ncol=5)
colnames(output)=c("B_slope", "L1", "numSpeciesPairs", "percentSignif","meanSlope")
output=data.frame(output)

overallCounter=1


# start simulations, running through range of parameters

for (B_slope_counter in 1:length(B_slope))
	{
	for (L1_counter in 1:length(L1))
		{
		for (numSpeciesPairs_counter in 1:length(numSpeciesPairs))	
			{	
	
			data_Xsp1=NULL  # vectors of simulated traits for unlagged trait X, lagged trait Y, for pairs of species sp1 and sp2; 						length=numSpeciesPairs
			data_Xsp2=NULL
			data_Ysp1=NULL
			data_Ysp2=NULL

			contrast_X=NULL  # vectors of unstandardized contrasts for X and Y, between sp1 and sp2
			contrast_Y=NULL

			brLength=NULL  # vector of branch lengths

			pvalList=NULL  # vector of statistics from test; length=numSimulation
			coefList=NULL
			coefListSign=NULL


			for (n in 1:numSimulation)  
				{
			
				for (i in 1:numSpeciesPairs[numSpeciesPairs_counter])
					{
								
					brLength[i]=runif(1,0,1)  # draw branch length from a uniform distribution
			
				# traitEvolVariance = brLength[i]  # under the assumption of speciational change, trait evolution variance is constant 						across different branch lengths.  However, it is possible to also set the evolutionary variance as a function of 						branch length here, and it has a big effect on the results (give it a try).  I prefer to keep the variance fixed, as 					assumes a speciational model, and Model 2 effectively sets variance proportional to branch length.
	
					data_Xsp1[i]=rnorm(1,0, traitEvolVariance)  # draw trait value for X, sp1
					data_Xsp2[i]=rnorm(1,0, traitEvolVariance)  # draw trait value for X, sp2
		
					if (brLength[i]<L1[L1_counter])  # need to take into account lag
					{
					data_Ysp1[i]=rnorm(1,0, traitEvolVariance) + B_slope[B_slope_counter] * data_Xsp1[i] *brLength[i]/L1[L1_counter]  
    				data_Ysp2[i]=rnorm(1,0, traitEvolVariance) + B_slope[B_slope_counter] * data_Xsp2[i] *brLength[i]/L1[L1_counter]
					}
					else   # don't take into account lag; Y is a simpler function of X
					{	
					data_Ysp1[i]=rnorm(1,0, traitEvolVariance) + B_slope[B_slope_counter] * data_Xsp1[i]
    				data_Ysp2[i]=rnorm(1,0, traitEvolVariance) + B_slope[B_slope_counter] * data_Xsp2[i]
					}

					if (data_Xsp1[i]>data_Xsp2[i])  # find direction of subtraction to make X contrasts positive
					{
					contrast_X[i]=data_Xsp1[i]-data_Xsp2[i]  
					contrast_Y[i]=data_Ysp1[i]-data_Ysp2[i]
					}

					else
					{
					contrast_X[i]=data_Xsp2[i]-data_Xsp1[i]
					contrast_Y[i]=data_Ysp2[i]-data_Ysp1[i]
					}
				}  # end i loop
	
			# plot(contrast_X, contrast_Y, cex=(brLength*2), main=paste("L1=",as.character(L1[L1_counter]),", r=", R_corr[R_corr_counter]))
	
			regressOut1=lm(contrast_Y ~ contrast_X - 1)  # regress contrasts through the origin
			regressOut2=lm(regressOut1$residuals ~ brLength)  # regress residuals from Out1 on branch length
			coefList[n]=summary(regressOut2)$coefficients[2]
			coefListSign[n]=sign(summary(regressOut2)$coefficients[2])  # store the coefficient's sign
		
			if(coefListSign[n]>0)   # positive slope, so record p-value
				{
				pvalList[n]=summary(regressOut2)$coefficients[8]
				}
			else
				{
				pvalList[n]=99  # set the pvalue to 99 (nonsensical) if the coefficient was negative
				}
		
			} # end n loop

		output$B_slope[overallCounter]=B_slope[B_slope_counter]
		output$L1[overallCounter]=L1[L1_counter]
		output$numSpeciesPairs[overallCounter]=numSpeciesPairs[numSpeciesPairs_counter]
		output$percentSignif[overallCounter]=sum(pvalList<0.1)/numSimulation # find the p-values less than 0.1, based on a one tailed test
		output$meanSlope[overallCounter]=mean(coefList)

		# print(sum(pvalList==99))

		overallCounter= overallCounter + 1
	
		}  # end numSpeciesPairs_counter loop
	}  # end L1_counter loop
} # end B_slope_counter loop
	


	