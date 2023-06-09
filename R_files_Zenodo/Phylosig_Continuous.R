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


##Lablength
#checking the phylogenetic signal for Lablength
Lablength <- dat$Lablength
names(Lablength) <- rownames(dat)

##Calculating the phylogenetic signal using Pagel's Lambda
L_Lablength <- phylosig(phy, Lablength, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
         control=list())
L_Lablength
# 1. Open pdf file
pdf("L_Lablength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(L_Lablength)
# 3. Close the file
dev.off()

##Calculating the phylogenetic signal using Blomberg's K
K_Lablength <- phylosig(phy, Lablength, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
         control=list())
K_Lablength
# 1. Open pdf file
pdf("K_Lablength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(K_Lablength)
# 3. Close the file
dev.off()


##Labwidth
#checking the phylogenetic signal for Labwidth
Labwidth <- dat$Labwidth
names(Labwidth) <- rownames(dat)

##Calculating the phylogenetic signal using Pagel's Lambda
L_Labwidth <- phylosig(phy, Labwidth, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
         control=list())
L_Labwidth
# 1. Open pdf file
pdf("L_Labwidth.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(L_Labwidth)
# 3. Close the file
dev.off()

##Calculating the phylogenetic signal using Blomberg's K
K_Labwidth <- phylosig(phy, Labwidth, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
         control=list())
K_Labwidth
# 1. Open pdf file
pdf("K_Labwidth.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(K_Labwidth)
# 3. Close the file
dev.off()


##Labnotch
#checking the phylogenetic signal for Labnotch
Labnotch <- dat$Labnotch
names(Labnotch) <- rownames(dat)

##Calculating the phylogenetic signal using Pagel's Lambda
L_Labnotch <- phylosig(phy, Labnotch, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                       control=list())
L_Labnotch
# 1. Open pdf file
pdf("L_Labnotch.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(L_Labnotch)
# 3. Close the file
dev.off()

##Calculating the phylogenetic signal using Blomberg's K
K_Labnotch <- phylosig(phy, Labnotch, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                       control=list())
K_Labnotch
# 1. Open pdf file
pdf("K_Labnotch.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(K_Labnotch)
# 3. Close the file
dev.off()


##Laterallength
#checking the phylogenetic signal for Laterallength
Laterallength <- dat$Laterallength
names(Laterallength) <- rownames(dat)

##Calculating the phylogenetic signal using Pagel's Lambda
L_Laterallength <- phylosig(phy, Laterallength, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                       control=list())
L_Laterallength
# 1. Open pdf file
pdf("L_Laterallength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(L_Laterallength)
# 3. Close the file
dev.off()

##Calculating the phylogenetic signal using Blomberg's K
K_Laterallength <- phylosig(phy, Laterallength, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                       control=list())
K_Laterallength
# 1. Open pdf file
pdf("K_Laterallength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(K_Laterallength)
# 3. Close the file
dev.off()


##Lateralwidth
#checking the phylogenetic signal for Lateralwidth
Lateralwidth <- dat$Lateralwidth
names(Lateralwidth) <- rownames(dat)

##Calculating the phylogenetic signal using Pagel's Lambda
L_Lateralwidth <- phylosig(phy, Lateralwidth, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                            control=list())
L_Lateralwidth
# 1. Open pdf file
pdf("L_Lateralwidth.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(L_Lateralwidth)
# 3. Close the file
dev.off()

##Calculating the phylogenetic signal using Blomberg's K
K_Lateralwidth <- phylosig(phy, Lateralwidth, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                            control=list())
K_Lateralwidth
# 1. Open pdf file
pdf("K_Lateralwidth.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(K_Lateralwidth)
# 3. Close the file
dev.off()


##FTLength
#checking the phylogenetic signal for FTLength
FTLength <- dat$FTLength
names(FTLength) <- rownames(dat)

##Calculating the phylogenetic signal using Pagel's Lambda
L_FTLength <- phylosig(phy, FTLength, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                           control=list())
L_FTLength
# 1. Open pdf file
pdf("L_FTLength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(L_FTLength)
# 3. Close the file
dev.off()

##Calculating the phylogenetic signal using Blomberg's K
K_FTLength <- phylosig(phy, FTLength, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                           control=list())
K_FTLength
# 1. Open pdf file
pdf("K_FTLength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(K_FTLength)
# 3. Close the file
dev.off()

##FTLength
#checking the phylogenetic signal for Fillength
Fillength <- dat$Fillength
names(Fillength) <- rownames(dat)

##Calculating the phylogenetic signal using Pagel's Lambda
L_Fillength <- phylosig(phy, Fillength, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                       control=list())
L_Fillength
# 1. Open pdf file
pdf("L_Fillength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(L_Fillength)
# 3. Close the file
dev.off()

##Calculating the phylogenetic signal using Blomberg's K
K_Fillength <- phylosig(phy, Fillength, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                       control=list())
K_Fillength
# 1. Open pdf file
pdf("K_Fillength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(K_Fillength)
# 3. Close the file
dev.off()

##Antherlength
#checking the phylogenetic signal for Antherlength
Antherlength <- dat$Antherlength
names(Antherlength) <- rownames(dat)

##Calculating the phylogenetic signal using Pagel's Lambda
L_Antherlength <- phylosig(phy, Antherlength, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL,
                        control=list())
L_Antherlength
# 1. Open pdf file
pdf("L_Antherlength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(L_Antherlength)
# 3. Close the file
dev.off()

##Calculating the phylogenetic signal using Blomberg's K
K_Antherlength <- phylosig(phy, Antherlength, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL,
                        control=list())
K_Antherlength
# 1. Open pdf file
pdf("K_Antherlength.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,4,4))
plot(K_Antherlength)
# 3. Close the file
dev.off()