##Load the packages
library(geiger)
##calculating and plotting disparity-through-time for a phylogenetic tree and phenotypic data
phy <- read.tree ("C_tree.nwk")
dat <- read.csv ("Laterallength.csv", header=T, row.names = 1)
##estimating the disparity through time
disparity(phy=NULL, dat, index = "avg.manhattan")

##Plotting the likelihood surface for the parameter ??
# 1. Open pdf file
pdf("Laterallength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
dtt_Laterallength <- dtt(phy, dat, index= "avg.manhattan",
    mdi.range=c(0,1), nsim=0, CI=0.95, plot=TRUE, calculateMDIp=F)
plot(dtt_Laterallength)
# 3. Close the file
dev.off()


disparity(phy=NULL, dat, index = c("avg.sq", "avg.manhattan", "num.states"))
dtt(phy, dat, index=c("avg.sq", "avg.manhattan", "num.states"),
    mdi.range=c(0,1), nsim=0, CI=0.95, plot=TRUE, calculateMDIp=F)