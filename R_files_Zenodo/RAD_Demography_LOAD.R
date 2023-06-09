library("ggplot2")
library("cowplot")
library("magrittr")
library("reshape")
library("Rmisc")
library("plyr")
library("lattice")
require(mgcv)

# GCTA relatedness --------------------------------------------------------

relateGCTA_val_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_noZ_GCTA_rel.grm.gz")
relateGCTA_ID_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_noZ_GCTA_rel.grm.id")

relateGCTA_val_MAD <- read.table("Berthelots_Madeira_subset_maf_trim_0.03_noZ_GCTA_rel.grm.gz")
relateGCTA_ID_MAD <- read.table("Berthelots_Madeira_subset_maf_trim_0.03_noZ_GCTA_rel.grm.id")

relateGCTA_val_SEL <- read.table("Berthelots_Selvagens_subset_maf_trim_0.03_noZ_GCTA_rel.grm.gz")
relateGCTA_ID_SEL <- read.table("Berthelots_Selvagens_subset_maf_trim_0.03_noZ_GCTA_rel.grm.id")


# Plink relatedness -------------------------------------------------------

relatePlink_val_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_noZ.rel",as.is = T)
relatePlink_ID_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_noZ.rel.id",as.is = T)

relatePlink_val_MAD <- read.table("Berthelots_Madeira_subset_maf_trim_0.03_noZ.rel",as.is = T)
relatePlink_ID_MAD <- read.table("Berthelots_Madeira_subset_maf_trim_0.03_noZ.rel.id",as.is = T)

relatePlink_val_SEL <- read.table("Berthelots_Selvagens_subset_maf_trim_0.03_noZ.rel",as.is = T)
relatePlink_ID_SEL <- read.table("Berthelots_Selvagens_subset_maf_trim_0.03_noZ.rel.id",as.is = T)


# PCA ---------------------------------------------------------------------

pca_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_noZ_Plinkrel0.2_PCA.eigenvec") 

pca_MAD <- read.table("Berthelots_Madeira_subset_maf_trim_0.03_noZ_Plinkrel0.2_PCA.eigenvec")



# LD by population --------------------------------------------------------

LD_EH <- read.table("LD_EH_chrs.ld", header = T)
LD_LP <- read.table("LD_LP_chrs.ld", header = T)
LD_TF <- read.table("LD_TF_chrs.ld", header = T)
LD_TEID <- read.table("LD_TEID_chrs.ld", header = T)
LD_GOM <- read.table("LD_GOM_chrs.ld", header = T)
LD_GC <- read.table("LD_GC_chrs.ld", header = T)
LD_FV <- read.table("LD_FV_chrs.ld", header = T)
LD_LZ <- read.table("LD_LZ_chrs.ld", header = T)
LD_GRA <- read.table("LD_GRA_chrs.ld", header = T)

LD_PS <- read.table("LD_PS_chrs.ld", header = T)
LD_DG <- read.table("LD_DG_chrs.ld", header = T)
LD_M <- read.table("LD_M_chrs.ld", header = T)

LD_SG <- read.table("LD_SG_chrs.ld", header = T)

