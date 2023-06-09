library(stats)
library(ggplot2)
library(vegan)



femaleraw <- read.csv("FemaleGenital_EMMs_PCA.csv", header = TRUE)
maleraw <- read.csv("MaleGenitalEMMs.csv", header = TRUE)
FST <- read.csv("FST.csv", header = TRUE)
row.names(FST) <- FST$X
FST <- FST[c(2:9)]
FST <- as.matrix(FST)
habitat <- read.csv("habitat.csv", header = TRUE)
row.names(habitat) <- habitat$X
habitat <- habitat[c(2:9)]
habitat <- as.matrix(habitat)


####All axes
femalerawdata <- data.frame(femaleraw$emmean.FemComp1, femaleraw$emmean.FemComp2, femaleraw$emmean.FemComp3)
row.names(femalerawdata) <- femaleraw$Site
femdist <- dist(femalerawdata, method = "euclidean", diag = TRUE, upper = FALSE)

malerawdata <- data.frame(maleraw[c(6:16)])
row.names(malerawdata) <- maleraw$Site
maldist <- dist(malerawdata, method = "euclidean", diag = TRUE, upper = FALSE)


#####Lock & Key
#Not corrected for habitat
mantel(femdist, FST, method="pearson", permutations=9999) #Significant positive relationship
mantel(maldist, FST, method="pearson", permutations=9999) #Barely NS positive relationship

#Corrected for habitat
mantel.partial(femdist, FST, habitat, method = "pearson", permutations = 9999) #Significant positive relationship
mantel.partial(maldist, FST, habitat, method = "pearson", permutations = 9999) #Barely NS positive relationship

######Sex stuff
mantel(maldist, femdist, method="pearson", permutations=9999) ##NS, but not by much
mantel.partial(femdist, maldist, FST, method = "pearson", permutations = 9999) ##NS
mantel.partial(femdist, maldist, habitat, method = "pearson", permutations = 9999) ##NS


####SUBSET of axes
femalerawdata <- data.frame(femaleraw$emmean.FemComp2, femaleraw$emmean.FemComp3)
row.names(femalerawdata) <- femaleraw$Site
femdist <- dist(femalerawdata, method = "euclidean", diag = TRUE, upper = FALSE)

malerawdata <- data.frame(maleraw[c(7, 14)])
row.names(malerawdata) <- maleraw$Site
maldist <- dist(malerawdata, method = "euclidean", diag = TRUE, upper = FALSE)


#####Lock & Key
#Not corrected for habitat
mantel(femdist, FST, method="pearson", permutations=9999) #Significant positive relationship
mantel(maldist, FST, method="pearson", permutations=9999) #Barely NS positive relationship

#Corrected for habitat
mantel.partial(femdist, FST, habitat, method = "pearson", permutations = 9999) #Significant positive relationship
mantel.partial(maldist, FST, habitat, method = "pearson", permutations = 9999) #Barely NS positive relationship

######Sex stuff
mantel(maldist, femdist, method="pearson", permutations=9999) ##NS, but not by much
mantel.partial(femdist, maldist, FST, method = "pearson", permutations = 9999) ##NS
mantel.partial(femdist, maldist, habitat, method = "pearson", permutations = 9999) ##NS


write.csv(femdist, file = "femdist.csv")

