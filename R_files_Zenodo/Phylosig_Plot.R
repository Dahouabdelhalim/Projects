##Loading the libraries
library(phytools)
library(ape)

##Loading the tree and data files
phy <- read.tree ("C_tree.nwk")
dat <- read.csv ("Continuous.csv", header=T, row.names = 1)

##tip labels of the tree
phy$tip.label

##plot the tree
plotTree(phy)

#row names of the csv file
rownames(dat)

#checking if the names of tree terminals and row names are matching
name.check(phy,dat)

#checking the phylogenetic signal for Lablength
Lablength <- dat$Lablength
names(Lablength) <- rownames(dat)

## extract characters of interest
ln.Lablength <-log(Lablength)

## compute phylogenetic signal K
K.Lablength<-phylosig(phy,ln.Lablength,
                     test=TRUE)
print(K.Lablength)
plot(K.Lablength)

## compute phylogenetic signal lambda
lambda.Lablength<-phylosig(phy,ln.Lablength,
                          method="lambda",test=TRUE)
print(lambda.Lablength)

par(mar=c(4,4,4,4))
plot(lambda.Lablength)
