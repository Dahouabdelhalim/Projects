# Script that performs a Monte-Carlo MFA Calculation
# =================================================================================================
# 
# Required input:     - TC : a matrix containing all the modal transfer coefficients from one 
#                       compartment to another (assuming that all these TC will have triangular 
#                       distributions). If one of these TC should not be defined this way, an NA
#                       value will enable to not create it.
# 
#                     - TC.uncertainty : a matrix of identical dimensions  than TC, containing the
#                       uncertainty to be used when creating a triangular distribution. Any NA
#                       value in this matrix will lead to a corresponding missing distribution
# 
#                     - ini.input : a vector containing the initial input for every compartment
#                       in the system
# 
#                     - ini.uncertainty : a vector of the same length than ini.input, containing the
#                       uncertainty to be used when creating a triangular distribution. Any NA
#                       value in this matrix will lead to a corresponding missing distribution
# 
#                     - possible additional distributions that need to be included in the code
#                       before the normalization step
# 
# Special instructions:   - little control of the input data is performed
#                         - the code was not tested for a use without row and column names
# 
# Output:   - Mass : a matrix containing the masses of the substances in each compartment
# 
# Version 2 with improved memory management
# 
# History of modifications: 
# 11.07.2017: Adapted for nano-silica and PEC calculation (selected code passage from the nano-SiO2 and IONP model of Yan et al. 2016a+b)
#                           
#                           
# Date of last modification: 11.07.2017
# =================================================================================================


##### INTRO #######################################################################################

# set number of simulation steps
SIM <- 100000
set.seed(150)

# source needed functions
source("functions.needed.R")

# package for importing Excel files (optional)
#library(xlsx)

# define arbitrary TC matrix (flows from column to rows)
TC <- read.table(file = "INSERT PATH TO FILE", header = TRUE, sep = ";", row.names = 1)
TC <- as.matrix(TC)


# take only part of the matrix (optional)
#TC <- TC[1:43,1:43]
# define arbitrary uncertainty matrix
TC.uncertainty <- matrix(0.5,dim(TC)[1], dim(TC)[2], dimnames = dimnames(TC))

# define the initial input vector, i.e. in which compartments there is an input (ProdVol in metric tonnes)
inp <- c(459100, rep(0,dim(TC)[1]-1))
# define uncertainty vector
inp.uncertainty <- rep(0.5, length(inp))
# give the vectors names (optional)
names(inp) <- names(inp.uncertainty) <- dimnames(TC)[[1]]

# first timer to check how long the simulation was running
timer <- proc.time()


# create a triangular distribution for every TC element using the uncertainty
# The output of this function will be a list with as many elements as there are compartments in
# the system. In each list element, the TC for the flows leaving the compartment will be stored.
# The do.trunc argument enables to truncate the distribution at 0 and 1 (in this case) or to
# any wished value. If no truncation is wished, write NA instead
TC.Distr <- triangulize.TC(modalvalues = TC,
                           uncertainty = TC.uncertainty,
                           N = SIM,
                           do.trunc = c(0,1))
						   
						   # create a triangular distribution for every TC element using the uncertainty
# The output of this function will be a 2-D matrix containing all the distributions of inputs.
inp.Distr <- triangulize.input(modalvalues = inp,
                               uncertainty = inp.uncertainty,
                               N = SIM)
							   
							   
							   # define custom distributions here (use str(TC.Distr) to see how the data is structured)

####WWTP#######################
###############################

## ##  ## ##  The Waste water entering WWTP but not treated due to overflows  ## ## ## ## ##

Over <- rnorm(SIM,0.032,0.004)
to_rem <- c(1-Over)


TC.Distr[["WWTP_con_rate_WWTP"]][["WWTP_over"]] <- Over

TC.Distr[["WWTP_con_rate_WWTP"]][["WWTP_rem"]] <- to_rem


##Waste water removal efficiency

STP_Jarvie <- rtriang.perc(0.94, perc = 0.5, 0.8*SIM, do.trunc=c(0,1))
STP_Robert <- rtriang.perc(0.97, perc = 0.5, 0.8*SIM, do.trunc = c(0,1))
STP_Liu <- rtriang.perc(0.99, perc = 0.5, 0.2*SIM, do.trunc = c(0,1))
STP_Den2006 <- rtriang.perc(0.90, perc = 0.5, 0.2*SIM, do.trunc = c(0,1))
STP_Den2005 <- rtriang.perc(0.988, perc = 0.5, 0.2*SIM, do.trunc = c(0,1))
STP_Pan2005 <- runif(0.2*SIM, 0.818, 0.927)
STP_Huang2004 <- runif(0.2*SIM, 0.945, 0.96)

STP_sum <- c(STP_Jarvie, STP_Robert, STP_Liu, STP_Den2006, STP_Den2005, STP_Pan2005, STP_Huang2004)
STP_rem <- sample(STP_sum, SIM, replace = FALSE)

#plot(density(STP_rem))
TC.Distr[["WWTP_rem"]][["Sewage_sludge_treat"]] <- STP_rem
TC.Distr[["WWTP_rem"]][["WWTP_effl"]] <- c(1-STP_rem)

### Waste incineration plant####
################################

#WIP Filtration efficiency (0.999) Walser et al.2012 and (0.995) Burtscher et al.2002
# Main fraction to filter ash, remains to filtered air (before acid washing)

filt_eff <- runif(SIM, 0.995, 0.999)

TC.Distr[["WIP_filter_eff"]][["WIP_filt_ash"]] <- filt_eff
TC.Distr[["WIP_filter_eff"]][["WIP_filt_air"]] <- c(1-filt_eff)

# WIP filtration to air (direct release before acid washing)

WIP_fil2_air <- rtriang.perc(0.5, perc = 0.5, N = SIM, do.trunc=c(0,1))

TC.Distr [["WIP_filt_air"]][["WIP_filt_air_air"]] <- WIP_fil2_air
TC.Distr[["WIP_filt_air"]][["WIP_acid_wash"]] <- c(1-WIP_fil2_air)

#WIP acid washing

WIP_acidw2air<- rtriang.perc(0.001, perc = 0.5, N = SIM, do.trunc=c(0,1))

TC.Distr[["WIP_acid_wash"]][["WIP_acid_wash_2air"]] <- WIP_acidw2air
TC.Distr[["WIP_acid_wash"]][["WIP_acid_wash_2LF"]] <- c(1-WIP_acidw2air) 


#Soil to surface water based on O'Brien et al 2010

#TC.Distr[["Soil"]][["Soil_2_SW"]] <-rtriang.perc(0.00549, perc = 0.5, N = SIM, do.trunc = c(0,1))

                                                            
# normalization step: make sure that all the flows going out of a compartment sum up to 1
TC.Distr.Norm <- normalize(TC.Distr)

##### SOLVE EQUATION ##############################################################################
# the solve.MC function will solve the equation system N times after having reconstructed the TC
# matrix for every iteration
Mass <- solve.MC(TC.Distr = TC.Distr.Norm,
                 inp.Distr = inp.Distr,
                 N = SIM)

# notify how long was needed for the whole calculation
message("Time needed for the simulation:")
print(proc.time() - timer) # second timer

##### GRAPHS ######################################################################################

# saves a pdf file of the histograms of the distributions to the working directory
pdf(file = "INSERTNAME.pdf",
    height = 7.5,
    width = 7.5,
    pointsize = 10)
par(mfrow = c(3,3), mar = c(3,3,3,1), mgp = c(1.5,0.5,0), xpd = F)
color <- adjustcolor(rainbow(dim(Mass)[1]), alpha.f = 0.6)
for(co in 1:dim(Mass)[1]){
  hist(Mass[co,], freq = F, xlab = paste(dimnames(Mass)[[1]][co], "(tonnes)"),
       ylab = "Probability density", main = dimnames(Mass)[[1]][co],
       col = color[co])
  legend("topright", c(paste("Mean:", round(mean(Mass[co,]), digits = 4)),
                       paste("SD:", round(sd(Mass[co,]), digits = 4))),
         bty = "n", cex = 0.85)
  box()
}
dev.off()

