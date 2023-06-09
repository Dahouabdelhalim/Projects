#call packages
library(phytools)
library(caper)
library(ape)
library(dplyr)

#call tree. 
mtree <-  read.newick("~/Desktop/Dataset_S5.tre")

###############################################################################################
######## Lines 14-60 is the code required to size correct data and gets PCA scores.       #####
######## If you want to quickly recreate results, jump to line 60 now and used Dataset.S4 #####
###############################################################################################

#call data. 
mdata <- read.csv("~/Desktop/Dataset_S3.csv") 
row.names(mdata) <- mdata$Species
#summary data file summarizes raw data (Datasets S1 and S2)

#sorted data into comparative dataset for caper package
weaponry <- comparative.data(phy = mtree, data = mdata, names.col = Species, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

#pgls using caper to correct traits for body size
model.1 <- pgls(Femur.width ~ PW, data = weaponry, lambda ="ML")
model.2 <- pgls(Femur.length ~ PW, data = weaponry, lambda ="ML")
model.3 <- pgls(Fulcrum ~ PW, data = weaponry, lambda ="ML")
model.4 <- pgls(Spine.length ~ PW, data = weaponry, lambda ="ML")
model.5 <- pgls(Tibial.length ~ PW, data = weaponry, lambda ="ML")
model.6 <- pgls(Location ~ PW, data = weaponry, lambda ="ML")
model.7 <- pgls(Mean.size.damage ~ PW, data = weaponry, lambda ="ML")
model.8 <- pgls(Median.size.damage ~ PW, data = weaponry, lambda ="ML")

#Table of size corrected data
weaponry.sc <- as.data.frame(cbind(model.1$residuals,model.2$residuals,model.3$residuals,model.4$residuals,model.5$residuals,model.6$residuals,model.7$residuals,model.8$residuals))
colnames(weaponry.sc) <- c("Femur.width", "Femur.length","Fulcrom", "Spine.length", "Tibial.length", "Location","Mean.size.damage", "Median.size.damage")
weaponry.sc.1 <- tibble::rownames_to_column(weaponry.sc, "Species")

#Add size corrected data to traits that are not size corrected, prominent spine number and proportion damaged
weaponry.sc.2 <- mdata[,c(1,9,10,13,14,15,16)]
weaponry.sc.3 <- merge(weaponry.sc.1, weaponry.sc.2, by="Species", all=TRUE)

#conduct phylogenetic PCA (pPCA)
PCA.data <- mdata[,c(3:9)]
pPCA <- phyl.pca(mtree, PCA.data, method="lambda",opt="ML")

#pPCA overview/summary
pPCA$Eval
summary(pPCA)

#reformat data
pPCA.b <- as.data.frame(pPCA$S)
pPCA.b <- tibble::rownames_to_column(pPCA.b, "Species")

#dataset used for analyses, includes size corrected traits and PC scores
dataset.4 <- merge(weaponry.sc.3, pPCA.b[,c(1:3)], by="Species", all=TRUE)
#write.csv(dataset.4, "~/Desktop/Dataset_S4.csv", row.names = TRUE)

################################################################################
######## To recreate main analyses reported in the manuscript pickup here  ######
################################################################################

#call data. 
dataset.4 <- read.csv("~/Desktop/Dataset_S4.csv") 
row.names(dataset.4) <- dataset.4$Species

#sorted data into comparative dataset for caper package
weaponry.2 <- comparative.data(phy = mtree, data = dataset.4, names.col = Species, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

### PGLS analyses using median damage severity
model.9 <- pgls(Femur.width ~ Median.size.damage, data = weaponry.2, lambda ="ML")
summary(model.9)
model.10 <- pgls(Femur.length ~ Median.size.damage, data = weaponry.2, lambda ="ML")
summary(model.10)
model.11 <- pgls(Fulcrom ~ Median.size.damage, data = weaponry.2, lambda ="ML")
summary(model.11)
model.12 <- pgls(Location ~ Median.size.damage, data = weaponry.2, lambda ="ML")
summary(model.12)
model.13 <- pgls(Spine.length ~ Median.size.damage, data = weaponry.2, lambda ="ML")
summary(model.13)
model.14 <- pgls(Tibial.length ~ Median.size.damage, data = weaponry.2, lambda ="ML")
summary(model.14)
model.15 <- pgls(Prominent.spine.number ~ Median.size.damage, data = weaponry.2, lambda ="ML")
summary(model.15)

model.16 <- pgls(PC2 ~ Median.size.damage, data = weaponry.2, lambda ="ML")
summary(model.16)

### PGLS analyses using median puncture number per damaged individual
model.17 <- pgls(Femur.width ~ Median.puncture.number, data = weaponry.2, lambda ="ML")
summary(model.17)
model.18 <- pgls(Femur.length ~ Median.puncture.number, data = weaponry.2, lambda ="ML")
summary(model.18)
model.19 <- pgls(Fulcrom ~ Median.puncture.number, data = weaponry.2, lambda ="ML")
summary(model.19)
model.20 <- pgls(Location ~ Median.puncture.number, data = weaponry.2, lambda ="ML")
summary(model.20)
model.21 <- pgls(Spine.length ~ Median.puncture.number, data = weaponry.2, lambda ="ML")
summary(model.21)
model.22 <- pgls(Tibial.length ~ Median.puncture.number, data = weaponry.2, lambda ="ML")
summary(model.22)
model.23 <- pgls(Prominent.spine.number ~ Median.puncture.number, data = weaponry.2, lambda ="ML")
summary(model.23)

model.24 <- pgls(PC2 ~ Median.puncture.number, data = weaponry.2, lambda ="ML")
summary(model.24)

### PGLS analyses using proportion damaged
model.25 <- pgls(Femur.width ~ Prop.male.damaged, data = weaponry.2, lambda ="ML")
summary(model.25)
model.26 <- pgls(Femur.length ~ Prop.male.damaged, data = weaponry.2, lambda ="ML")
summary(model.26)
model.27 <- pgls(Fulcrom ~ Prop.male.damaged, data = weaponry.2, lambda ="ML")
summary(model.27)
model.28 <- pgls(Location ~ Prop.male.damaged, data = weaponry.2, lambda ="ML")
summary(model.28)
model.29 <- pgls(Spine.length ~ Prop.male.damaged, data = weaponry.2, lambda ="ML")
summary(model.29)
model.30 <- pgls(Tibial.length ~ Prop.male.damaged, data = weaponry.2, lambda ="ML")
summary(model.30)
model.31 <- pgls(Prominent.spine.number ~ Prop.male.damaged, data = weaponry.2, lambda ="ML")
summary(model.31)

model.32 <- pgls(PC2 ~ Prop.male.damaged, data = weaponry.2, lambda ="ML")
summary(model.32)

### PGLS analyses investingating whether our 3 damage metrics are related
model.33 <- pgls(Median.puncture.number ~ Prop.male.damaged, data = weaponry.2, lambda ="ML")
summary(model.33)
model.34 <- pgls(Median.size.damage ~ Prop.male.damaged, data = weaponry.2, lambda ="ML")
summary(model.34)
model.35 <- pgls(Median.size.damage ~ Median.puncture.number, data = weaponry.2, lambda ="ML")
summary(model.35)

### PGLS analyses investigating whether variation in body size, variation in femur size, or absolute weapon size is associated with the frequency of damage observed among species
model.36 <- pgls(Prop.male.damaged~SD.PW, data = weaponry.2, lambda ="ML")
summary(model.36)
model.37 <- pgls(Prop.male.damaged~PC1, data = weaponry.2, lambda ="ML")
summary(model.37)
model.38 <- pgls(Prop.male.damaged~SD.Femur.width, data = weaponry.2, lambda ="ML")
summary(model.38)



