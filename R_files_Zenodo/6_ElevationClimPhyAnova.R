# This R script tests whether there is an association between 
# pollination syndrome and species' median elevation or climate 
# niche using phylogenetic anova.  It uses the phylogenetic tree 
# file "chronogram_v8_2.tre" and the data file "species.elev.csv" 

library(geiger)
library(phytools) #for phylANOVA
library(car) #levene's test

setwd("~/Dropbox/Costus grant/Dena and Kathleen/FloralEvolution/Dryad_Deposit")

# chronogram
phy <- read.nexus("chronogram_v8_2.tre")
plot(phy)

#import species elev traits
phy.data <- read.csv("species.elev.csv")
rownames(phy.data) <- phy.data[,c("species.phyname")]

#occurrence record stats
occ.data <- read.csv("occurrences.elev.csv")
occ.data.summary<- as.data.frame(table(occ.data$acceptedScientificName))
names(occ.data.summary) <- c("acceptedScientificName", "N_occ")
tmp <- merge(phy.data, occ.data.summary, all.x=T, all.y=F, by="acceptedScientificName")
sum(tmp$N_occ) #total number of occurrences used
nrow(tmp) #number taxa
range(tmp$N_occ) #range per taxa
mean(tmp$N_occ) #mean per taxa

# use random forest predicted syndromes
phy.data["Costus_dirzoi_98079","syndrome"]<-"H"
phy.data["Costus_varzearum_19252","syndrome"]<-"H"

(name.check(phy, phy.data) -> phyOverlap)
drop.tip(phy, phyOverlap$tree_not_data) -> phyComparativeTree
plot(phyComparativeTree, cex=0.6)

attach(phy.data)

###
###PHY ANOVA on elevation and climate traits
###

#ELEV
dat = phy.data$elev_median
dat = as.data.frame(dat)
rownames(dat) = phy.data$species.phyname
dat = na.omit(dat)
grp = as.factor(phy.data$syndrome)
names(grp)=phy.data$species.phyname
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant
# Analysis of Variance Table
#
# Response: dat
# Df  Sum-Sq Mean-Sq F-value  Pr(>F) Pr(>F) given phy
# group      1  125105  125105 0.63968 0.42803           0.6863
# Residuals 45 8800820  195574 

hum <- subset(phy.data, syndrome=="H")
bee <- subset(phy.data, syndrome=="B")
range(hum$elev_median)
range(bee$elev_median)
table(phy.data$syndrome)

#climate PC1
dat = phy.data$PC1
dat = as.data.frame(dat)
rownames(dat) = phy.data$species.phyname
dat = na.omit(dat)
grp = as.factor(phy.data$syndrome)
names(grp)=phy.data$species.phyname
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant
# Analysis of Variance Table
#
# Response: dat
# Df  Sum-Sq Mean-Sq F-value  Pr(>F) Pr(>F) given phy
# group      1  1.1281 1.12806  2.1961 0.14533           0.4825
# Residuals 45 23.1147 0.51366       

#climate PC2
dat = phy.data$PC2
dat = as.data.frame(dat)
rownames(dat) = phy.data$species.phyname
dat = na.omit(dat)
grp = as.factor(phy.data$syndrome)
names(grp)=phy.data$species.phyname
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant
# Analysis of Variance Table
#
# Response: dat
# Df  Sum-Sq Mean-Sq F-value  Pr(>F) Pr(>F) given phy
# group      1  0.1287 0.12870 0.29961 0.58683           0.7912
# Residuals 45 19.3307 0.42957         

leveneTest(phy.data$elev_median ~ phy.data$syndrome)
leveneTest(phy.data$PC1 ~ phy.data$syndrome)
leveneTest(phy.data$PC2 ~ phy.data$syndrome)

####
## PLOT
####

pdf(file="Fig6.pdf", width=8,height=3.5)
par(mfrow=c(1,3))
boxplot(phy.data$elev_median ~ phy.data$syndrome, col=c("#1F639B","#ED553B"),ylab="Elevation (m)", horizontal=F, names=c("Orchid bee","Hummingbird"), varwidth=T, xlab="", cex.lab=1.2)
mtext("A", side=3, adj=0, cex=0.8)
boxplot(phy.data$PC1 ~ phy.data$syndrome, col=c("#1F639B","#ED553B"),ylab="Climate PC1", horizontal=F, names=c("Orchid bee","Hummingbird"), varwidth=T, xlab="", cex.lab=1.2)
mtext("B", side=3, adj=0, cex=0.8)
boxplot(phy.data$PC2 ~ phy.data$syndrome, col=c("#1F639B","#ED553B"),ylab="Climate PC2", horizontal=F, names=c("Orchid bee","Hummingbird"), varwidth=T, xlab="", cex.lab=1.2)
mtext("C", side=3, adj=0, cex=0.8)
dev.off()
