
# GCTA relatedness --------------------------------------------------------

# Canaries

relateGCTA_ID_CAN$number <- NA
relateGCTA_ID_CAN$number <- c(1:182)

relateGCTA_val_CAN$Indiv1 <- NA
relateGCTA_val_CAN$Indiv2 <- NA

# change the population names 
relateGCTA_ID_CAN$V1 <- gsub("C_01GRA", "GRA", relateGCTA_ID_CAN$V1)
relateGCTA_ID_CAN$V1 <- gsub("C_02LZ", "LZ", relateGCTA_ID_CAN$V1)
relateGCTA_ID_CAN$V1 <- gsub("C_03FV", "FV", relateGCTA_ID_CAN$V1)
relateGCTA_ID_CAN$V1 <- gsub("C_04GC", "GC", relateGCTA_ID_CAN$V1)
relateGCTA_ID_CAN$V1 <- gsub("C_05TF", "TF", relateGCTA_ID_CAN$V1)
relateGCTA_ID_CAN$V1 <- gsub("C_06TEID", "TEID", relateGCTA_ID_CAN$V1)
relateGCTA_ID_CAN$V1 <- gsub("C_07GOM", "GOM", relateGCTA_ID_CAN$V1)
relateGCTA_ID_CAN$V1 <- gsub("C_08LP", "LP", relateGCTA_ID_CAN$V1)
relateGCTA_ID_CAN$V1 <- gsub("C_09EH", "EH", relateGCTA_ID_CAN$V1)

relateGCTA_val_CAN$Indiv1 <- relateGCTA_ID_CAN$V2[match(relateGCTA_val_CAN$V1, relateGCTA_ID_CAN$number)]
relateGCTA_val_CAN$Indiv2 <- relateGCTA_ID_CAN$V2[match(relateGCTA_val_CAN$V2, relateGCTA_ID_CAN$number)]

relateGCTA_val_CAN$Pop1 <- relateGCTA_ID_CAN$V1[match(relateGCTA_val_CAN$V1, relateGCTA_ID_CAN$number)]
relateGCTA_val_CAN$Pop2 <- relateGCTA_ID_CAN$V1[match(relateGCTA_val_CAN$V2, relateGCTA_ID_CAN$number)]

# subset only within population comparisons but exclude same individual comparisons 
relateGCTA_val_CAN <- subset(relateGCTA_val_CAN, Pop1 == Pop2)
relateGCTA_val_CAN <- subset(relateGCTA_val_CAN, Indiv1 != Indiv2)

# Add population as factor and determine levels.
relateGCTA_val_CAN$Pop1 <- factor(relateGCTA_val_CAN$Pop1)
relateGCTA_val_CAN$Pop1 <- factor(relateGCTA_val_CAN$Pop1, levels = c("EH", "GOM", "LP", "TEID", "TF", "GC", "FV", "LZ", "GRA"))


# Maderian

relateGCTA_ID_MAD$number <- NA
relateGCTA_ID_MAD$number <- c(1:60)

relateGCTA_val_MAD$Indiv1 <- NA
relateGCTA_val_MAD$Indiv2 <- NA

relateGCTA_ID_MAD$V1 <- gsub("M_11M", "M", relateGCTA_ID_MAD$V1)
relateGCTA_ID_MAD$V1 <- gsub("M_12PS", "PS", relateGCTA_ID_MAD$V1)
relateGCTA_ID_MAD$V1 <- gsub("M_13DG", "DG", relateGCTA_ID_MAD$V1)

relateGCTA_val_MAD$Indiv1 <- relateGCTA_ID_MAD$V2[match(relateGCTA_val_MAD$V1, relateGCTA_ID_MAD$number)]
relateGCTA_val_MAD$Indiv2 <- relateGCTA_ID_MAD$V2[match(relateGCTA_val_MAD$V2, relateGCTA_ID_MAD$number)]

relateGCTA_val_MAD$Pop1 <- relateGCTA_ID_MAD$V1[match(relateGCTA_val_MAD$V1, relateGCTA_ID_MAD$number)]
relateGCTA_val_MAD$Pop2 <- relateGCTA_ID_MAD$V1[match(relateGCTA_val_MAD$V2, relateGCTA_ID_MAD$number)]

relateGCTA_val_MAD <- subset(relateGCTA_val_MAD, Pop1 == Pop2)
relateGCTA_val_MAD <- subset(relateGCTA_val_MAD, Indiv1 != Indiv2)

relateGCTA_val_MAD$Pop1 <- factor(relateGCTA_val_MAD$Pop1)
relateGCTA_val_MAD$Pop1 <- factor(relateGCTA_val_MAD$Pop1, levels = c("M", "PS", "DG"))

GCTA_MAD_0.2rel <- subset(relateGCTA_val_MAD, V4 > 0.2)


# Selvagens

relateGCTA_ID_SEL$number <- NA
relateGCTA_ID_SEL$number <- c(1:20)


relateGCTA_val_SEL$Indiv1 <- NA
relateGCTA_val_SEL$Indiv2 <- NA

relateGCTA_ID_SEL$V1 <- gsub("S_10SG", "SG", relateGCTA_ID_SEL$V1)


relateGCTA_val_SEL$Indiv1 <- relateGCTA_ID_SEL$V2[match(relateGCTA_val_SEL$V1, relateGCTA_ID_SEL$number)]
relateGCTA_val_SEL$Indiv2 <- relateGCTA_ID_SEL$V2[match(relateGCTA_val_SEL$V2, relateGCTA_ID_SEL$number)]

relateGCTA_val_SEL <- subset(relateGCTA_val_SEL, Indiv1 != Indiv2)


GCTA_SEL_0.2rel <- subset(relateGCTA_val_SEL, V4 > 0.2)


# Plink relatedness -------------------------------------------------------

# Canaries
colnames(relatePlink_val_CAN) <- relatePlink_ID_CAN$V2 # add the header names to the square pairwise matrix
rownames(relatePlink_val_CAN) <- relatePlink_ID_CAN$V2

# create table of values for each combination of individuals within one population. i.e. all La Palma birds with all other La Palma birds. 
x <- data.frame(t(combn(relatePlink_ID_CAN$V2,2)))
colnames(x) <- c("IND1","IND2")

x$rel <- NA
x$POP1 <- NA
x$POP2 <- NA


for(i in 1:nrow(x))
{
  x$rel[i] <- relatePlink_val_CAN[paste(x$IND1[i]),paste(x$IND2[i])]
  x$POP1[i] <- subset(relatePlink_ID_CAN,V2 == x$IND1[i])$V1
  x$POP2[i] <- subset(relatePlink_ID_CAN,V2 == x$IND2[i])$V1
}

x2 <- subset(x,POP1 == POP2)

# Change population names
x2$POP1 <- gsub("C_01GRA", "GRA", x2$POP1)
x2$POP1 <- gsub("C_02LZ", "LZ", x2$POP1)
x2$POP1 <- gsub("C_03FV", "FV", x2$POP1)
x2$POP1 <- gsub("C_04GC", "GC", x2$POP1)
x2$POP1 <- gsub("C_05TF", "TF", x2$POP1)
x2$POP1 <- gsub("C_06TEID", "TEID", x2$POP1)
x2$POP1 <- gsub("C_07GOM", "GOM", x2$POP1)
x2$POP1 <- gsub("C_08LP", "LP", x2$POP1)
x2$POP1 <- gsub("C_09EH", "EH", x2$POP1)

x2$POP1 <- factor(x2$POP1)
x2$POP1 <- factor(x2$POP1, levels = c("EH", "GOM", "LP", "TEID", "TF", "GC", "FV", "LZ", "GRA"))

Plink_relatedness_CAN <- x2

# Maderian

colnames(relatePlink_val_MAD) <- relatePlink_ID_MAD$V2 # add the deader names to the square pairwise matrix
rownames(relatePlink_val_MAD) <- relatePlink_ID_MAD$V2

x <- data.frame(t(combn(relatePlink_ID_MAD$V2,2)))
colnames(x) <- c("IND1","IND2")

x$rel <- NA
x$POP1 <- NA
x$POP2 <- NA


for(i in 1:nrow(x))
{
  x$rel[i] <- relatePlink_val_MAD[paste(x$IND1[i]),paste(x$IND2[i])]
  x$POP1[i] <- subset(relatePlink_ID_MAD,V2 == x$IND1[i])$V1
  x$POP2[i] <- subset(relatePlink_ID_MAD,V2 == x$IND2[i])$V1
}

x2 <- subset(x,POP1 == POP2)

x2$POP1 <- gsub("M_11M", "M", x2$POP1)
x2$POP1 <- gsub("M_12PS", "PS", x2$POP1)
x2$POP1 <- gsub("M_13DG", "DG", x2$POP1)

x2$POP1 <- factor(x2$POP1)
x2$POP1 <- factor(x2$POP1, levels = c("M", "PS", "DG"))

Plink_relatedness_MAD <- x2

Plink_MAD_0.2rel <- subset(Plink_relatedness_MAD, rel > 0.2)



# Selvagens
colnames(relatePlink_val_SEL) <- relatePlink_ID_SEL$V2 # add the header names to the square pairwise matrix
rownames(relatePlink_val_SEL) <- relatePlink_ID_SEL$V2

x <- data.frame(t(combn(relatePlink_ID_SEL$V2,2)))
colnames(x) <- c("IND1","IND2")

x$rel <- NA


for(i in 1:nrow(x))  
{
  x$rel[i] <- relatePlink_val_SEL[paste(x$IND1[i]),paste(x$IND2[i])]
}

Plink_relatedness_SEL <- x


# PCA ---------------------------------------------------------------------

colnames(pca_CAN) <- c("Population","Individual",paste0("PC",c(1:4)))

pca_CAN$Population <- gsub("C_01GRA", "GRA", pca_CAN$Population)
pca_CAN$Population <- gsub("C_02LZ", "LZ", pca_CAN$Population)
pca_CAN$Population <- gsub("C_03FV", "FV", pca_CAN$Population)
pca_CAN$Population <- gsub("C_04GC", "GC", pca_CAN$Population)
pca_CAN$Population <- gsub("C_05TF", "TF", pca_CAN$Population)
pca_CAN$Population <- gsub("C_06TEID", "TEID", pca_CAN$Population)
pca_CAN$Population <- gsub("C_07GOM", "GOM", pca_CAN$Population)
pca_CAN$Population <- gsub("C_08LP", "LP", pca_CAN$Population)
pca_CAN$Population <- gsub("C_09EH", "EH", pca_CAN$Population)

pca_CAN$Population <- factor(pca_CAN$Population)
pca_CAN$Population <- factor(pca_CAN$Population, levels = c("EH", "GOM", "LP", "TEID", "TF", "GC", "FV", "LZ", "GRA"))


colnames(pca_MAD) <- c("Population","Individual",paste0("PC",c(1:4)))

pca_MAD$Population <- gsub("M_11M", "M", pca_MAD$Population )
pca_MAD$Population <- gsub("M_12PS", "PS", pca_MAD$Population)
pca_MAD$Population <- gsub("M_13DG", "DG", pca_MAD$Population)

pca_MAD$Population <- factor(pca_MAD$Population)
pca_MAD$Population <- factor(pca_MAD$Population, levels = c("M", "PS", "DG"))



# LD by population --------------------------------------------------------

# calculate distance between two SNPs
LD_EH$BP_dis <- LD_EH$BP_B - LD_EH$BP_A
LD_EH <- subset(LD_EH, LD_EH$BP_dis > 0)
LD_LP$BP_dis <- LD_LP$BP_B - LD_LP$BP_A
LD_LP <- subset(LD_LP, LD_LP$BP_dis > 0)
LD_TF$BP_dis <- LD_TF$BP_B - LD_TF$BP_A
LD_TF <- subset(LD_TF, LD_TF$BP_dis > 0)
LD_TEID$BP_dis <- LD_TEID$BP_B - LD_TEID$BP_A
LD_TEID <- subset(LD_TEID, LD_TEID$BP_dis > 0)
LD_GOM$BP_dis <- LD_GOM$BP_B - LD_GOM$BP_A
LD_GOM <- subset(LD_GOM, LD_GOM$BP_dis > 0)
LD_GC$BP_dis <- LD_GC$BP_B - LD_GC$BP_A
LD_GC <- subset(LD_GC, LD_GC$BP_dis > 0)
LD_FV$BP_dis <- LD_FV$BP_B - LD_FV$BP_A
LD_FV <- subset(LD_FV, LD_FV$BP_dis > 0)
LD_LZ$BP_dis <- LD_LZ$BP_B - LD_LZ$BP_A
LD_LZ <- subset(LD_LZ, LD_LZ$BP_dis > 0)
LD_GRA$BP_dis <- LD_GRA$BP_B - LD_GRA$BP_A
LD_GRA <- subset(LD_GRA, LD_GRA$BP_dis > 0)

LD_DG$BP_dis <- LD_DG$BP_B - LD_DG$BP_A
LD_DG <- subset(LD_DG, LD_DG$BP_dis > 0)
LD_PS$BP_dis <- LD_PS$BP_B - LD_PS$BP_A
LD_PS <- subset(LD_PS, LD_PS$BP_dis > 0)
LD_M$BP_dis <- LD_M$BP_B - LD_M$BP_A
LD_M <- subset(LD_M, LD_M$BP_dis > 0)

LD_SG$BP_dis <- LD_SG$BP_B - LD_SG$BP_A
LD_SG <- subset(LD_SG, LD_SG$BP_dis > 0)

