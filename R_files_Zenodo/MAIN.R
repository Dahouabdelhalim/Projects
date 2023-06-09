#rm(list=ls())
library("boot", lib.loc="C:/Program Files/R/R-3.2.2/library")
#Combined.data <- read.csv("C:/Users/Derek/Dropbox/Fairoff/Cricket Experiment/DATA/Generation 3 data/Generation 3 Final Data/Combined.data2016", sep="")
setwd("C:/Users/Derek/Dropbox/Files for User Friendly Female Preference")
source("C:/Users/Derek/Dropbox/Files for User Friendly Female Preference/Functions.for.Female.preference.R")

###################################IMPORTANT###############################
# DATA MUST BE ENTERED AS WITH EACH MALE ON A SEPARATE LINE AND PAIRS FOLLOWING EACH OTHER
# Trait values for males are assumed to be positive
# First line is a header
######################### THINGS TO BE SET BY USER #########################  
Data <- Test.Data01           # Input data file Here assumed to have been read in
Zero.one.Data  <- "YES"     # The female preference metric is continuous. "YES" if 0,1 data
set.seed(1949)              # Value of seed to be set for bootstrap
Female.col      <- 3        # Column for measure of female preference for given male (here = number of approaches to male)
Cols.to.be.used <- 2        # Columns to be tested.  Here there are three columns.  They do not have to be sequential
N.BOOTS         <- 100       # Number of bootstraps.  Here set at 100 for demonstration.  Use at least 1000
##############################################################################
# These data are reorganised such that the first column holds the trait value of the "left" ("focal") male
# Column 2 holds the trait value of the other male
# The first lines here produce a file compliant with this
par(mfcol=c(3,2)) # Set up graphics page   
Left.males  <- Data[Data$Male=="Left",]
Right.males <- Data[Data$Male=="Right",]
if(Zero.one.Data=="NO")  Preference  <- Left.males[,Female.col]/(Left.males[,Female.col]+Right.males[,Female.col]) # preference metric
if(Zero.one.Data=="YES") Preference  <- Left.males[,Female.col] 
Names       <- colnames(Data, do.NULL = TRUE, prefix = "col")# Get plot labels
n           <- ncol(Data)  # Number of columns
Print.flag  <- "NO"  # YES = print anything else no
####################################################################################################
NCOLS        <- length(Cols.to.be.used)
Output.Table <- matrix(0,NCOLS,12)
Output.Table <- data.frame(Output.Table)
for (I in 1:NCOLS) #n)  # Iterate over traits
{    
  Trait.col         <- Cols.to.be.used[I]
  Label             <- Names[Trait.col]  # Get names of male traits from column
  Output.Table[I,1] <- Label
  Trait             <- cbind(Left.males[,Trait.col], Right.males[,Trait.col], Preference) # Col1 has choose =0 and col 2 has choose = 1
  Trait             <- na.omit(Trait)         # Omit pairs with missing values
  ################################## CALCULATE STATS##################################
  # Find the number of zeroes but do not eliminate
  # This is a check for bad data
  # call duration and vol can be zero but no others can be  
  n1.z            <- nrow(Trait)
  #Zero            <- Trait[Trait[,1]!=0,]
  #Zero            <- Zero [Zero [,2]!=0,]
  #n2.z            <- nrow(Zero )
  Male.Variance   <- var(c(Trait[,1], Trait[,2]) )#Data[,Trait.col], na.rm=TRUE)    # Variance of all males
  Male.SD         <- sqrt(Male.Variance)
  Female.Variance <- var(Trait[,2])                       # Variance of chosen males
  Male.Mean       <- mean(c(Trait[,1], Trait[,2]))        # Mean of all males
  Female.Mean     <- mean(Trait[,2])                      # Mean of chosen males
  Min.Trait       <- Male.Mean-3*Male.SD #min(Trait)      # Minimum value of trait
  Max.Trait       <- Male.Mean+3*Male.SD #max(Trait)      # Maximum value of trait
  Obs.min         <- min(c(Trait[,1], Trait[,2]))
  Obs.max         <- max(c(Trait[,1], Trait[,2]))
  N.X             <- 100                                  # Length of absolute preference values
  X.trial         <- seq(Min.Trait, Max.Trait, length=N.X)# Vector of absolute preference values
  x               <- rbind(Trait[,1], Trait[,2])
  hist(x,xlab=Label, main=Label)                          # Histogram of male trait values
  if(Zero.one.Data=="NO")
  {
    Direction      <- Directional.Preference(Trait, Label, Print.flag)
    Results        <- Stabilizing.Preference(X.trial,Trait, N.X, Min.Trait, Max.Trait, Label, N.pairs, Pred, Print.flag)
  }
  if(Zero.one.Data=="YES")
  {
    Direction      <- Directional.Preference01(Trait)
    Results        <- Stabilizing.Preference01(X.trial,Trait, N.X, Min.Trait, Max.Trait, Label, N.pairs, Pred)
  }
  
  Output         <- matrix(0,1, 12)
  Output[1,1:7]  <- c(Trait.col,Direction,Results) 
  OutputB        <- boot(data=Trait, statistic=BOOTSTRAP, R=N.BOOTS)
  boot.values    <- OutputB$t[,1]
  write.table(x=boot.values, file = "boot.values.txt", row.names=FALSE,append = FALSE ) # Output bootstrap vales     
  Boot.est        <- mean(boot.values)
  Boot.SE         <- sd(boot.values)
  hist( OutputB$t[,1], xlab="All bootstraps", main=Label) #Plot histogram of bootstrap values
  Output[1,8:9]   <- c(Boot.est, Boot.SE)
  #*****************Output Results*****************
  print(c("##########",Label,"#########"))
  #if(n1.z!=n2.z) print(c("Zeros deleted = ", n1.z,n2.z))
  print("Minimum and maximum values")
  print(c(Obs.min,Obs.max))
  print("Means(All males, Chosen males)")
  print(c(Male.Mean,Female.Mean))
  print("Ratio of means, Chosen by Female/Total ")
  AVG.RATIO <- Female.Mean/Male.Mean
  print(Female.Mean/Male.Mean)
  print("Variances(All males, Chosen males)")
  print(c(Male.Variance, Female.Variance))
  print("Ratio of variances, Chosen by Female/Total")
  VAR.RATIO <- Female.Variance/Male.Variance
  print(Female.Variance/Male.Variance)
  #print("Ratio of predicted value/Male value")
  #print(Results[1]/Male.Mean)
  #print(c("coeff of variation = ",sqrt(Male.Variance)/Male.Mean))
  # *************************************** Calculate best X***************************************
 
  print("*******************Analysis*******************")
  print(c("Sample size = ",n1.z))
  print("***Tests for Directional selection (P1, P2, slope)***")
  print(Output[1,2:4])
  print("***Test for Stabilizing selection (X, slope, P)***")
  print(Output[1,5:7])
  print("*** Bootstrap estimates (mu, SE)***")
  print(Output[1,8:9])
  print("*****************************************************************************")
  print("*****************************************************************************")
  #Output.Table[I,] <- c(Label,Male.Mean, AVG.RATIO,VAR.RATIO,Output[1,2:9] ) 
} # End of Traits.col loop 
#col.Names <- c("Trait", "Mean","Ratio of means", "Ratio of variances", "P.Pref.D.Pref", "P.Diff.D.Pref",
#               "Slope", "Est.X", "Slope","P", "MeanBoot", "SE")
#colnames(Output.Table)    <- col.Names
#write.table(x=Output.Table, file = "Output.Table.txt", row.names=FALSE,append = FALSE )