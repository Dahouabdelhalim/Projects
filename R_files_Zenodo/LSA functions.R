#--------------------------------------------------------------------------------------
#
#  July 2022
#  Code associated with publication:
#  "Kinetic modulation of bacterial hydrolases by microbial community structure in coastal waters" 
#   Authors: N. Abad, A. Uranga, B. Ayo, J.M. Arrieta,Z. Baña, I. Artolozaga, I. Azúa, J. Iriberri, Santos J. González-Rojí and M. Unanue
#
#
#  Creative Commons Licence: Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
#--------------------------------------------------------------------------------------

### Adapted from Ruan et al (2006): doi:10.1093/bioinformatics/btl417

### The modified module from the original source is displayed below. The explanations of the original code 
### have been retained for the sake of better comprehension. Specific modifications related to a function are 
### indicated by the comment "#New: modified function".


################################
# LocalSimilarityAnalysis.R  #
################################
#
# The following functions are implemented:
#
# . normalTransform <- function(rawMatrix)
#
#	. tempChangeSeq(dataX)
#
#	. colSum2One(dataX)
#
#	. normalization(x, dims=1, na.rm=FALSE)
#
#	. LocalSimilarity3(lsTS1, lsTS2, maxDelay=1, numTimePoints=1, scale=F)  #New: modified function
#
#	. SigTestMaxDelay(data1, N=dim(data1)[[2]], delay=0, permu=1000)
#
#	. sigTestingCor1(lsTS1, lsTS2, method="pearson", permu=1000)
#
#	. sigTestingCor2(dataX, method="pearson", permu=1000)
#
#	. PlotOtuPair2(idx=15, ls.res=ls.res.all, ylab="rel. abund.(normalized)", data=dataXe)
#
#	. PlotOtuPair3(idx=15, ls.res=ls.res.all, ylab="rel. abund.", data=dataTe)
#
#	. CovaryOtuSif(tmp5) 
#
#
# Last updated: Apr 20, 2006, Quansong Ruan.  @All rights reserved.
#############################################################################

################################
# normalTransform(rawMatrix): 
# 
# For a rawMatrix[pxq] matrix, with each column as an OTU, and row as time points,
# this function perform a rank normal score transformation [Ker-chau Li, PNAS 2002]
# 
# Output: rankScoreMatrix: a normal score matrix of the same size as rawMatrix
############################### 

normalTransform <- function(rawMatrix)
{
	
dimRM <- dim(rawMatrix);
rowNum <- dimRM[[1]];
colNum <- dimRM[[2]];
rankScoreMatrix <- rawMatrix;

for(i in 1:colNum)
{
	rk <- rank(rawMatrix[,i]);
	rankScoreMatrix[,i] <- qnorm(rk/(rowNum+1));
}

rankScoreMatrix
}

###################################################
# temporalChangeSeq(dataX):
#
# This function calculate the temporal change sequence for each time series
#
# Arguments:
# =========
# 	dataX: a M by N non-negative data matrix. each row in dataX is a time series
#
# Return:
# ======
#    a M by (N-1) data matrix containing the temporal change sequences of 
# the time series in dataX.
#
# Note:
# =====
#
#
#  The temp change between a[i] and a[i+1] is defined as follows:
# 
#     Ta[i]  = 0, if a[i+1]=a[i]=0
#
#	     = (a[i+1]-a[i]) / (a[i+1]+a[i]), otherwise
#
#
# Mar 24, 2006
###################################################
tempChangeSeq <- function(dataX)
{
	data <- dataX
	data <- as.matrix(data)
	dim1 <- dim(data)
	
	dt1 <- data[,2:dim1[[2]]] - data[,1:(dim1[[2]]-1)]
	dt2 <- data[,2:dim1[[2]]] + data[,1:(dim1[[2]]-1)]


	# trick: if dt2[i,j]=0, then dt1[i,j] must be 0 as well
	# by definition, Ta[i,j] should be 0. so for efficient 
	# operation, if dt2[i,j]=0, set dt2[i,j] <- 108 (any 
	# non-zero number does not change the result)
	# 
	dt2[dt2<=0] <- 108

	return (dt1 / dt2) # entry-wise division
}





####################################################
# colSum2One(dataX): 
#
#    This function re-calculate the percentage for each OTU at the same time spot
#
#################################################### 
colSum2One <- function(dataX)
{

	tmp <- dataX
	dim1 <- dim(tmp)
	csum1 <- colSums(tmp)
	for(i in 1:dim1[[2]])
	{
		tmp[,i] <- tmp[,i]/csum1[i]*100	
	}
	
	return(tmp)
}


#########################################################################
# normalization <- function(x, dims=1, na.rm=FALSE)
###########################
#
#	This function performs normalization of given data (vector, array, 
# matrix, data frame)
# Arguments: 
# =========	
#	    x: an array of two or more dimensions, containing numeric, 
#	       complex, integer or logical values, or a numeric data frame.
#
#	 dims: which dimension is to normalize over. 1-row, 2-column
#
#	na.rm: remove NA's if TRUE, ignore NA's oterwise
#
# RETURN:
#    Row(column) normalized data of the same size as x
#   
##########################################################################
normalization <- function(x, dims=2, na.rm=FALSE)
{

 tmpX <- x; 
 dim1 <- dim(tmpX)

 if(is.null(dim1)) 	# integer or  1-dim array (vector)
 {
      tmpX <- (tmpX - mean(tmpX, na.rm=na.rm)) / sd(tmpX, na.rm=na.rm)
 }
 else   		# at least two dimensions
 {
   if(dims==1)   tmpX <- t(tmpX);
   
   dim1 <- dim(tmpX)

   Tmp.m <- colMeans(tmpX, na.rm=na.rm) 
   Tmp.sd <- apply(tmpX, 2, sd, na.rm=na.rm)    # standard deviation of each column
                  #sd(as.matrix(tmpX), na.rm=na.rm) 

   vect1 <- t(t(rep(1,dim1[[1]])));
   tmpX <- (tmpX - vect1 %*% Tmp.m) / (vect1 %*% Tmp.sd)
   
   if(dims==1) tmpX <- t(tmpX);

 }

   tmpX
} # normalization <- function(x, dims=2, na.rm=FALSE)

# a<- matrix(1:10, nrow=2,ncol=5)
# normalization(a, dim=2)
# normalization(a[1,])
# normalization(a[,1], dim=2)


##################################################################
# LocalSimilarity3<- function(lsTS1, lsTS2, maxDelay=1, numTimePoints=1, scale=F) 
#
#  This function computes the local similarity score for two sequences.
#
# INPUT:
# ======
#
#	lsTS1, lsTS2	: sequences to copute LS score
#	maxDelay	: maximum time shift allowed in computing LS score.
#	numTimePoints	: length of the sequences for computing LS score
#	scale		: If TRUE, perform normalization first; False, otherwise.
#
# RETURN:
# ======
# 
#  A five element vector contains: c(scoreMax, (startX-length), (startY-length), length, PosOrNeg)
#
LocalSimilarity3<- function(lsTS1, lsTS2, maxDelay=1, numTimePoints=1, scale=F) 
{

if(!is.array(lsTS1))
{
	ls.TS1x <- t(lsTS1)
	ls.TS2x <- t(lsTS2)
}

ls.TS1 <- ls.TS1x
ls.TS2 <- ls.TS2x


	
if(scale==T)
{
	lsTS1 <- (lsTS1-mean(lsTS1))/sd(lsTS1)  # Splus: stdev(lsTS1)
	lsTS2 <- (lsTS2-mean(lsTS2))/sd(lsTS2)  # Splus: stdev(lsTS2)	
}


scoreMatrixPos <- matrix(0, numTimePoints+1,numTimePoints+1);
scoreMatrixNeg <- matrix(0, numTimePoints+1, numTimePoints+1);

scoreMax <- 0.;
PosOrNeg <- 0.;
startX <- 0; # start of sub seq in dt1 from back to front
startY <- 0; # start of sub seq in dt2 from back to front

Thresh2 <- 0.000001

for(i in 2:(numTimePoints+1))
{
	for(j in 2:(numTimePoints+1))
	{
		if(abs(i-j) > maxDelay)
			next;
		
		scoreMatrixPos[i,j] <- scoreMatrixPos[i-1,j-1] + lsTS1[i-1] * lsTS2[j-1];
		if(scoreMatrixPos[i,j] < 0)
			scoreMatrixPos[i,j] <- 0;
			
		scoreMatrixNeg[i,j] <- scoreMatrixNeg[i-1,j-1] - lsTS1[i-1] * lsTS2[j-1];
		if(scoreMatrixNeg[i,j] < 0)
			scoreMatrixNeg[i,j] <- 0;
		
		if(scoreMatrixPos[i,j] > scoreMax)
		{
			scoreMax <- scoreMatrixPos[i,j];
			startX <- i;
			startY <- j;
			PosOrNeg <- 1;	
		}	
		
		if(scoreMatrixNeg[i,j] > scoreMax)
		{
			scoreMax <- scoreMatrixNeg[i,j];
			startX <- i;
			startY <- j;
			PosOrNeg <- 0;	
		}
		
	}
}

if(PosOrNeg == 1)
{
	for(i in 1:numTimePoints)
    	if(scoreMatrixPos[startX-i,startY-i]<=Thresh2)
		{    		break;
       }
}
else {
	for(i in 1:numTimePoints)
       if(scoreMatrixNeg[startX-i, startY-i]<=Thresh2)
		{     		break;
       }

}
length = i;
				
return(c(scoreMax, (startX-length), (startY-length), length, PosOrNeg));
}


######################
# compute the significance level of the LS score
#
# LocalSimilarity3<- function(lsTS1, lsTS2, maxDelay=1, numTimePoints=1, scale=F) 
######################
sigTesting3 <- function(stOtu1, stOtu2, numPermu=10, maxDelay=1, numTimePoints=1, scale=F, indexNai)
{
	scoreArray <- rep(0.0, numPermu+1);
	
	scoreMax1 <- LocalSimilarity3Nai(stOtu1, stOtu2, maxDelay, numTimePoints=length(stOtu1), scale)[indexNai]; #New: modified function
	scoreArray[1] <- scoreMax1;
	highScoreCnt <- 1;
	
	for(idx in 1:(numPermu-1)){
		dtt1 <- stOtu1[sample(numTimePoints)];
		dtt2 <- stOtu2[sample(numTimePoints)];
		scoreTmp <- LocalSimilarity3Nai(dtt1, dtt2, maxDelay, numTimePoints=length(dtt1), scale)[indexNai]; #New: modified function

		scoreArray[idx+1] <- scoreTmp;
		highScoreCnt <- highScoreCnt + (scoreTmp >= scoreMax1);
	}
	
	pValue <- 1.0 * highScoreCnt / numPermu;
	
	return(pValue);
	#scoreArray[numPermu+1] <- pValue;
	#return(scoreArray)
}

#ls01x15 <-sigTesting3(data1[,1], data1[,15], maxDelay=1, numTimePoints=length(data1[,1]), scale=T, 10)
#ls01x16 <- sigTesting3(data1[,1], data1[,16], maxDelay=1, numTimePoints=length(data1[,1]), scale=T, 100)
#ls01x16 <- sigTesting3(data1[,1], data1[,16], maxDelay=1, numTimePoints=length(data1[,1]), scale=T, 10)


####################################################################
#SigTestMaxDelay <- function(data1, N=dim(data1)[[2]], delay=0, permu=1000)
#
# 
# This function computes the p-value matrix for the corresponding LS score matrix
# computed by LocalSimilarity3() above.
#
# Arguments: 
# ==========
#	    data1: a M by N data matrix/data frame. NOTE: should be at least two columns
#
#	 delay:  maximal time shift allowed between two subsequences.
#
#	permu: = number of permutations for compuating the p-value of the LS score
#
# RETURN:
# =======
#    	lsMatrixN: A N x 10 matrix: c("rIdx", "cIdx", "LSscore", "startR", "startC", "length", 
#					"PorN", "pValue", "cor", "corpVal"), saved in file 
#	paste("lsMatrixN.R.MaxDelay", delay, ".ran", sample(10000:999999,1), ".txt", sep="")
#
# NOTE: 
# =====
#	dataX should have at least two columns	
#
#
#########################################################################
SigTestMaxDelay <- function(data1, N=dim(data1)[[2]], delay=0, permu=1000)
{

dim1 <- dim(data1)

if(is.null(dim1))
  return (0)

# lsMatrix to store pairwise results
lsMatrixN <- matrix(0, N*N/2, 2+5+1+2);

rIdx <- 0;
for(i in 1:(N-1))
{
	for(j in (i+1):N)
	{
		lsTmp <- LocalSimilarity3(data1[,i], data1[,j], maxDelay=delay, numTimePoints=length(data1[,1]), scale=F);
		
		rIdx <- rIdx + 1;
		lsMatrixN[rIdx, 1] <- i;				# 	first seq index
		lsMatrixN[rIdx, 2] <- j;				#  	second seq index
		lsMatrixN[rIdx, 3] <- lsTmp[1]; 			# 	LS score
		lsMatrixN[rIdx, 4] <- lsTmp[2]; 			# 	startX
		lsMatrixN[rIdx, 5] <- lsTmp[3];				# 	startY
		lsMatrixN[rIdx, 6] <- lsTmp[4]; 			#	length
		lsMatrixN[rIdx, 7] <- lsTmp[5];				#	PosOrNeg
		lsMatrixN[rIdx, 8] <- sigTesting3(data1[,i], data1[,j], numPermu=permu, maxDelay=delay, numTimePoints=length(data1[,1]), scale=F)

		corTmp <- cor(data1[,i], data1[,j]);
		lsMatrixN[rIdx,  9] <- corTmp;
		lsMatrixN[rIdx, 10] <- 0.5 + sign(corTmp) * (0.5 - pt(corTmp * sqrt((dim1[[1]]-1)/(1-corTmp^2)), df=(dim1[[1]]-1)) );

	}
}

# 'normalize' the LS scores
lsMatrixN[,3] <- lsMatrixN[,3]/dim1[[1]]

dimnames(lsMatrixN)[[2]] <- c("rIdx", "cIdx", "LSscore", "startR", "startC", "length", "PorN", "pValue", "cor", "corpVal")

# remove the empty rows
lsMatrixN <- lsMatrixN[lsMatrixN[,1]>=1,]

fileName <- paste("lsMatrixN.R.MaxDelay", delay, ".ran", sample(10000:999999,1), ".txt", sep="")
write.table(lsMatrixN, paste(fileName, sep=""))

dim(lsMatrixN)

############## To read data:
# lsMatrixN.rd <- read.table(paste(RWorkPath, fileName, sep=""), header=TRUE)
# lsMatrixN.rd[1:3,]
# lsMatrixN[1:3,]

} # end of SigTestMaxDelay() ...


#SigTestMaxDelay(data1=dataXe[,1:10], N=dim(dataXe[,1:10])[[2]],delay=1)

#SigTestMaxDelay(data1=DataSpots, N=dim(DataSpots)[[2]],delay=1)
#SigTestMaxDelay(data1=DataSpots, N=dim(DataSpots)[[2]],delay=0)



######################################
# PlotOtuPair2:
#
# . This function is used to plot the significant OTU pairs.
#
# Input:
# ======
# . idx: the index number of the pair to be plotted in ls.res
# . data: a 35 by 72 matrix
# . ls.res: a N by 11 result matrix, a sublist of ls.res.all
#
####################################
PlotOtuPair2<- function(idx=1,
			data = dataXe, 
			ls.res = ls.res.all, 
			xdisp=0,
			xcor= xdisp+0.25*dim(data)[[1]],
			ydisp=0,
			ycor=ydisp + 0.95*max(data[,c(ls.res[idx,1],ls.res[idx,2])]),
			xlab="time point",
			ylab= "rel. abund."
			)
{

i <- ls.res[idx,1]
j <- ls.res[idx,2]

dimnm <- dimnames(data)[[2]]
dimnm <- c(substr(dimnm[1:58], start=2, stop=5), dimnm[59:72])
otuNames.dpbin.all <- paste("otu ", 1:58, ": ", dimnm[1:58], sep="")
otuNames.dpbin.all <- c(otuNames.dpbin.all, paste("env ", 59:72, ": ", dimnm[59:72], sep=""))


ylim1 <- 1.05*min(data[,c(ls.res[idx,1],ls.res[idx,2])])
ylim2 <- 1.05*max(data[,c(ls.res[idx,1],ls.res[idx,2])])

plot(data[,j],  type="l", lty=2, lwd=0.5, ylim=c(ylim1, ylim2), xlab=xlab, ylab=ylab, cex.lab=1.8, cex.axis=1.5)
lines(data[,i], type="l", lty=1, lwd=0.5, col=1)

ix <- (ls.res[idx,5]):(ls.res[idx,5]+ls.res[idx,6]-1)
lines(ix, data[ix,j], type="l", lty=2, lwd=3., col=1)

ix <- (ls.res[idx,4]):(ls.res[idx,4]+ls.res[idx,6]-1)
lines(ix, data[ix,i], type="l", lwd=3., col=1)


# start of match points
startX <- c(ls.res[idx,4], ls.res[idx,5])
startY <- c(data[startX[1], i], data[startX[2],j])
points(startX, startY,pch='O', cex=1.2)

# end of match points
startX <- c(ls.res[idx,4], ls.res[idx,5])+ls.res[idx,6]-1
startY <- c(data[startX[1], i], data[startX[2],j])
points(startX, startY,pch=20, cex=2.5)

legend(xcor, ycor, legend=paste(otuNames.dpbin.all[c(i,j)], sep=""), 
	lty=c(1:2), lwd=2, cex=1.2)

#text(xi,yi,paste(i, sep=""), col=1, cex=1.6)
#text(xj,yj,paste(j,sep=""), col=1, cex=1.6)
grid()

}

#PlotOtuPair2(idx=15, ls.res=ls.res.all, ylab="rel. abund.(normalized)", data=dataXe)
#PlotOtuPair2(idx=15, ls.res=ls.res.all, ylab="rel. abund.", data=dataTe)


################################################
# PlotOtuPair3():
#
# This function is modified from PlotOtuPair2()
#
# This function plots the aligned view of the OTU pairs with
# significant local similarity scores with some time delay
# only the matching subsequences are ploted.
#
####################
PlotOtuPair3<- function(idx=1,
			data = dataTe, 
			ls.res = ls.res.all, 
			xlab="time point",
			ylab= "rel. abund."
			)
{
i <- ls.res[idx,1]
j <- ls.res[idx,2]

dimnm <- dimnames(data)[[2]]
dimnm <- c(substr(dimnm[1:58], start=2, stop=5), dimnm[59:72])
otuNames.dpbin.all <- paste("otu ", 1:58, ": ", dimnm[1:58], sep="")
otuNames.dpbin.all <- c(otuNames.dpbin.all, paste("env ", 59:72, ": ", dimnm[59:72], sep=""))


dt1 <- data[(ls.res[idx,4]):(ls.res[idx,4]+ls.res[idx,6]-1),ls.res[idx,1]]
dt2 <- data[(ls.res[idx,5]):(ls.res[idx,5]+ls.res[idx,6]-1),ls.res[idx,2]]

ylim1 <- 1.05*min(dt1,dt2)
ylim2 <- 1.05*max(dt1,dt2)

plot(dt2,  type="l", lty=2, lwd=2.5, ylim=c(ylim1, ylim2), xlab=xlab, ylab=ylab, cex.lab=1.8, cex.axis=1.5)
lines(dt1, type="l", lty=1, lwd=2.5, col=1)

xdisp=0
xcor= xdisp+0.15*length(dt1)
ycor=0.99*max(dt1,dt2)
legend(xcor, ycor, legend=paste(otuNames.dpbin.all[c(i,j)], sep=""), 
	lty=c(1:2), lwd=2, cex=1.2)

#text(xi,yi,paste(i, sep=""), col=1, cex=1.6)
#text(xj,yj,paste(j,sep=""), col=1, cex=1.6)
grid()

}

#PlotOtuPair3(idx=15, ls.res=ls.res.all, ylab="rel. abund.(normalized)", data=dataXe)
#PlotOtuPair3(idx=15, ls.res=ls.res.all, ylab="rel. abund.", data=dataTe)
#PlotOtuPair2(idx=15, ls.res=ls.res.all, ylab="rel. abund.(normalized)", data=dataXe)


############################################################
# CovaryOtuSif(data): 
#
# . This function create SIF file from LSA result. SIF file is needed 
# by Cytoscape (software) to generate Co-varying graphs
#
# Input:
#
# . data: a N by 11 matrix, a subset of ls.res.all
#
# Output:
#
# . lsla.covary.Notu.ranxxxxx.sif: a file needed by Cytoscape 
#               to draw covarying graph
#
##########################################  
CovaryOtuSif <- function(data, otu.names=otuNames.dpbin.all) {

tmp <- data
dim1 <- dim(tmp)
covary <- matrix(0, ncol=3, nrow=dim1[[1]])
#covary[1,2] <- "empty" 

for(idx in 1:dim1[[1]])
{
	covary[idx,1] <- otu.names[tmp[idx,1]]
	covary[idx,3] <- otu.names[tmp[idx,2]]

	if (tmp[idx,4]==tmp[idx,5]) { # no time delay, undirected
		if(tmp[idx,7] >= 1 )
		{
			covary[idx,2] <- "pu " #1 # pu: positive, undirected
		} else {
			covary[idx,2] <- "nu " #2 # nu: negative, undirected
		}
	} else if (tmp[idx,4] > tmp[idx,5])  { # otu2 is ahead of otu1, arrow points to otu1
		if(tmp[idx,7] >= 1 )
		{
			covary[idx,2] <- "pdl" #3 # pdl: pos, directed to the left
		} else {
			covary[idx,2] <- "ndl" #4 # ndl: neg, directed to the left
		}
	} else {  # otu1 is ahead of otu2, arrow points to otu2
		if(tmp[idx,7] >= 1 )
		{
			covary[idx,2] <- "pdr" #5 # pdr: pos, directed to the right
		} else {
			covary[idx,2] <- "ndr" #6 # ndr: neg, directed to the right
		}
	}
}

#covary

fileName <- paste("lsla.covary.", dim(covary)[[1]], "all.ran", sample(10000,1),".sif", sep="")
write.table(covary, file=fileName, quote=FALSE, col.names=FALSE, row.names=FALSE)

} # CovaryOtuSif <- function(data)

#CovaryOtuSif(tmp5)

#############################################################
# CovaryOtuSif2 <- function(data, otu.names=otuNames.dpbin.all)
#
# This function is modified from the one above. It creates SIF file and 
# also gives the LS score for each pair.
##########################################  
CovaryOtuSif2 <- function(data, otu.names=otuNames.dpbin.all) {

tmp <- data
dim1 <- dim(tmp)
covary <- matrix(0, ncol=5, nrow=dim1[[1]])
#covary[1,2] <- "empty" 

for(idx in 1:dim1[[1]])
{
	covary[idx,1] <- otu.names[tmp[idx,1]]
	covary[idx,3] <- otu.names[tmp[idx,2]]
	covary[idx,4] <- round(tmp[idx,3], digits=4)
	covary[idx,5] <- tmp[idx,4]-tmp[idx,5] # time shift

	if (tmp[idx,4]==tmp[idx,5]) { # no time delay, undirected
		if(tmp[idx,7] >= 1 )
		{
			covary[idx,2] <- "pu " #1 # pu: positive, undirected
		} else {
			covary[idx,2] <- "nu " #2 # nu: negative, undirected
		}
	} else if (tmp[idx,4] > tmp[idx,5])  { # otu2 is ahead of otu1, arrow points to otu1
		if(tmp[idx,7] >= 1 )
		{
			covary[idx,2] <- "pdl" #3 # pdl: pos, directed to the left
		} else {
			covary[idx,2] <- "ndl" #4 # ndl: neg, directed to the left
		}
	} else {  # otu1 is ahead of otu2, arrow points to otu2
		if(tmp[idx,7] >= 1 )
		{
			covary[idx,2] <- "pdr" #5 # pdr: pos, directed to the right
		} else {
			covary[idx,2] <- "ndr" #6 # ndr: neg, directed to the right
		}
	}
}

#covary

fileName <- paste("lsla.covary.", dim(covary)[[1]], "all.edge.ran", sample(10000,1),".sif", sep="")
write.table(covary, file=fileName, quote=FALSE, col.names=FALSE, row.names=FALSE)

} # CovaryOtuSif2 <- function(data)

#CovaryOtuSif2(tmp5)



