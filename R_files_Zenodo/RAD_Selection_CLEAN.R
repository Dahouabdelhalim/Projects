
# Fst ---------------------------------------------------------------------

# Locus Fst ---------------------------------------------------------------

which(!(Locus_Fst_CAN$SNP %in% Locus_Fst_MAD$SNP))
which(!(Locus_Fst_MAD$SNP %in% Locus_Fst_CAN$SNP))

length(unique(Locus_Fst_CAN$SNP))
length(unique(Locus_Fst_MAD$SNP))

# Exclude SNPs that were unmapped
Locus_Fst_MAD_mapped <- subset(Locus_Fst_MAD, CHR >0)
Locus_Fst_CAN_mapped <- subset(Locus_Fst_CAN, CHR >0)

# just quickly check the Fst scores of the SNPs that were not able to be chrom mapped 
Locus_Fst_CAN_not_mapped <- subset(Locus_Fst_CAN, CHR == 0)
Locus_Fst_MAD_not_mapped <- subset(Locus_Fst_MAD, CHR == 0)


#length(unique(Locus_Fst_CAN_mapped$SNP))
#length(unique(Locus_Fst_MAD_mapped$SNP))


# Now must reorder the SNPs and remove SNPs that only occur in one data set 
inboth <- intersect(Locus_Fst_CAN_mapped$SNP,Locus_Fst_MAD_mapped$SNP)
Locus_Fst_CAN_reduced <- subset(Locus_Fst_CAN_mapped,SNP %in% inboth)
Locus_Fst_MAD_reduced <- subset(Locus_Fst_MAD_mapped,SNP %in% inboth)

#all(Locus_Fst_MAD_reduced$SNP == Locus_Fst_CAN_reduced$SNP)

##create a dataframe of both groups SNPs 
#str(Locus_Fst_CAN_reduced)
#str(Locus_Fst_MAD_reduced)
Fst <- data.frame(SNP_CAN = Locus_Fst_CAN_reduced$SNP, SNP_MAD = Locus_Fst_MAD_reduced$SNP,Fst_Canaries = Locus_Fst_CAN_reduced$FST, Fst_Maderia = Locus_Fst_MAD_reduced$FST)

#all(Fst$SNP_CAN == Fst$SNP_MAD)

# lower Fst bound of zero
Fst$Fst_Canaries[Fst$Fst_Canaries <0] <- 0
Fst$Fst_Maderia[Fst$Fst_Maderia <0] <- 0

#identify the loci that are at zero for one Arch
inoneonly <- setdiff(Locus_Fst_CAN_mapped$SNP,Locus_Fst_MAD_mapped$SNP)
#length(inoneonly)
Locus_Fst_CAN_mapped_not_in_MAD <- subset(Locus_Fst_CAN_mapped, SNP %in% inoneonly) # 1676 
Locus_Fst_MAD_mapped_not_in_CAN <- subset(Locus_Fst_MAD_mapped, SNP %in% inoneonly) # 0, this is probably not suprising considering it's likey that variation is lower in Madeira 

Locus_Fst_CAN_mapped_not_in_MAD$Fst_Madeira <- 0
Locus_Fst_CAN_mapped_not_in_MAD$Fst_Canaries <- Locus_Fst_CAN_mapped_not_in_MAD$FST
Locus_Fst_CAN_mapped_not_in_MAD$Fst_Canaries[Locus_Fst_CAN_mapped_not_in_MAD$Fst_Canaries <0] <- 0

additional_snps <- as.data.frame(Locus_Fst_CAN_mapped_not_in_MAD$SNP)
additional_snps$SNP_CAN <- additional_snps$`Locus_Fst_CAN_mapped_not_in_MAD$SNP`
additional_snps$Fst_Maderia <- Locus_Fst_CAN_mapped_not_in_MAD$Fst_Madeira
additional_snps$Fst_Canaries <- Locus_Fst_CAN_mapped_not_in_MAD$Fst_Canaries
colnames(additional_snps)[1] <- "SNP_MAD"

# add this data set to the bottom of Fst dataframe

All_SNPs <- rbind.data.frame(additional_snps, Fst)
#nrow(All_SNPs)


# Eigenvecotor PCA --------------------------------------------------------

colnames(egvec_CAN) <- c("Population","Individual",paste0("PC",c(1:4)))
colnames(egvec_MAD) <- c("Population","Individual",paste0("PC",c(1:4)))


# Eigenvalues Manhatton Plot ----------------------------------------------

colnames(eg1_CAN) <- c("SNP","CHR", "BP", "RefAllele", "AltAllele", "freq", "Beta", "SE", "Chi",  "P", "PGC", "n1", "freq1", "n2", "freq2", "Fst")
eg1_CAN[,c("Beta","SE","Chi","P","PGC","Fst")] <- apply(eg1_CAN[,c("Beta","SE","Chi","P","PGC","Fst")],2,as.numeric)
#head(eg1_CAN[is.na(eg1_CAN$Beta),])
eg1_CAN <- eg1_CAN[complete.cases(eg1_CAN),]

eg2_CAN[,c("Beta","SE","Chi","P","PGC","Fst")] <- apply(eg2_CAN[,c("Beta","SE","Chi","P","PGC","Fst")],2,as.numeric)
#head(eg2_CAN[is.na(eg2_CAN$Beta),])
eg2_CAN <- eg2_CAN[complete.cases(eg2_CAN),]


colnames(eg1_CAN) <- c("SNP","CHR", "BP", "RefAllele", "AltAllele", "freq", "Beta", "SE", "Chi",  "P", "PGC", "n1", "freq1", "n2", "freq2", "Fst")
eg1_MAD[,c("Beta","SE","Chi","P","PGC","Fst")] <- apply(eg1_MAD[,c("Beta","SE","Chi","P","PGC","Fst")],2,as.numeric)
#head(eg1_MAD[is.na(eg1_MAD$Beta),])
eg1_MAD <- eg1_MAD[complete.cases(eg1_MAD),]

eg2_MAD[,c("Beta","SE","Chi","P","PGC","Fst")] <- apply(eg2_MAD[,c("Beta","SE","Chi","P","PGC","Fst")],2,as.numeric)
#head(eg2_MAD[is.na(eg2_MAD$Beta),])
eg2_MAD <- eg2_MAD[complete.cases(eg2_MAD),]

