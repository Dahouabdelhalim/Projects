# Change Directory to where genepop files are stored.

# load diveRsity
# install.packages("diveRsity") in case it isn't installed.
library(diveRsity)

# load data
Pops <- readGenepop("MicrosCorrect_GENEPOP.txt")

# Calculation of Jost's D, Fst, Gst and G'st:
x <- diffCalc("MicrosCorrect_GENEPOP.txt",pairwise=TRUE,fst=TRUE,outfile = "Gen_Dist",bs_pairwise = TRUE,boots=1000)

# Testing for Correlation between within-population He and differentiation measures:
# Warning: fastDivPart is very heavy and wil take 10-15 minutes to run with 1000 boots.
y <- fastDivPart("Da,Db,S,J,A (Nuc)-GENEPOP.txt",gp=3,pairwise=TRUE,fst=TRUE,bs_pairwise = TRUE,boots=1000)

# Optional print to pdf
#pdf("Distance_Heterzygosity_Cor.pdf)

corPlot("MicrosCorrect_GENEPOP.txt.txt")

#dev.off()


#_______ Calculate Allelic Richness

# load hierfstat
#install.packages("hierfstat") in case it isn't installed
library(hierfstat)

# load data:
Pops2 <- read.fstat("MicrosCorrect_FSTAT-FSTAT.dat")

# Calculate Allelic Richness:
# Optional print to txt file:
#sink("Allelic_Richness.txt")

allelic.richness(Pops2,diploid=TRUE)

#sink()