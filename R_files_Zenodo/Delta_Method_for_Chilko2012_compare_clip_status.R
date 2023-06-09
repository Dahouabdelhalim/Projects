#//******************************************************
#//	Programmer: Kintama Research Services (kintama.com)
#//	Project Name: Ken Jeffries et al. Immune response genes and pathogen presence predict migration survival in wild salmon smolts. Accepted by Molecular Ecology. 
#//	Program Date:  Aug 2013
#//	Version:1
#//	Comments:  Delta method for cumulative survival estimates for 2012 V7 tagged Chilko Lake sockeye salmon submitted to Dryad (http://datadryad.org/)
#//	
#//	
#//******************************************************/

# Compute the SE of the cumulative survival probabilites

# Used rMARK to generate survival rates for EACH interval
# Exported the estimates of phi and p and the VCV to Excel worksheets
#
# Original script assumed that the first entries were survival rates in order
# from first to last interval. 


# Read in the REAL estimates (phi and p) and the full (phi and p) vcv

AllEst<-read.csv("Estimates_V7_p67.csv", header=TRUE)
AllVCV<-read.csv("realVCmatrix_V7_p67.csv", header=TRUE)


V7GILL<-subset(AllEst, Group=="V7-2L_GILL")
V7<-subset(AllEst, Group=="V7-2L_NO_GILL")

#------------------V7GILL
# Extract the estimates of phi and the vcv of phi

phiEst<-V7GILL$Estimates

phiVCV<-AllVCV[1:6,1:6]

# Compute the cumulate survival rates
cumEst<-cumprod(phiEst)

# Compute the derivative matrix for the delta method
# See the course notes for the pattern

# The "outer(x,y,"*")" function (outer product) creates a matrix whose
#   (i,j)th element is x[i]*y[j]
# This gives the whole matrix, but the lower triangle must be zeroed
deriv <- outer(1/phiEst, cumEst, '*')

# Zero out the lower triangle part
deriv <- deriv * outer(1:length(phiEst),1:length(phiEst),"<=")

# Now to compute the delta method
#     vcv(cumproduct) = t(derive) vcv(phi) deriv
# where t() is the transpose operator

cumVCV <- t(deriv) %*% as.matrix(phiVCV) %*% deriv


# Create the final report
V7GILLoutput <- cbind(1:length(phiEst),
                phiEst,
                sqrt(diag(as.matrix(phiVCV))),   # the se of estPhi
                cumEst,
                sqrt(diag(cumVCV)))              # the se of the cum survivals
# label the columns
colnames(V7GILLoutput) <- c('i', 'Est_Phi', 'SE(phi)', 'Cum_Prod', 'SE(Cum_Prod)')
#----------------------

#------------------V7 Intact
# Extract the estimates of phi and the vcv of phi

phiEst<-V7$Estimates

phiVCV<-AllVCV[7:12,7:12]

# Compute the cumulate survival rates
cumEst<-cumprod(phiEst)

# Compute the derivative matrix for the delta method
# See the course notes for the pattern

# The "outer(x,y,"*")" function (outer product) creates a matrix whose
#   (i,j)th element is x[i]*y[j]
# This gives the whole matrix, but the lower triangle must be zeroed
deriv <- outer(1/phiEst, cumEst, '*')

# Zero out the lower triangle part
deriv <- deriv * outer(1:length(phiEst),1:length(phiEst),"<=")

# Now to compute the delta method
#     vcv(cumproduct) = t(derive) vcv(phi) deriv
# where t() is the transpose operator

cumVCV <- t(deriv) %*% as.matrix(phiVCV) %*% deriv


# Create the final report
V7output <- cbind(1:length(phiEst),
                phiEst,
                sqrt(diag(as.matrix(phiVCV))),   # the SE of estPhi
                cumEst,
                sqrt(diag(cumVCV)))              # the SE of the cum survivals
# label the columns
colnames(V7output) <- c('i', 'Est_Phi', 'SE(phi)', 'Cum_Prod', 'SE(Cum_Prod)')
#----------------------



V7GILLoutput
V7output

