library("ggplot2")
library("Rmisc")
library("qqman")
library("devtools")
library("ggrepel")
library("ggman")
library("cowplot")
library("magrittr")
library("reshape")
library("plyr")
library("lattice")

# Fst ---------------------------------------------------------------------

# Genome wide Fst ---------------------------------------------------------

### Plot Locus Fst across the genome ###
Locus_Fst_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_Plinkrel0.2.fst", header = T,stringsAsFactors = F)
Locus_Fst_MAD <- read.table("Berthelots_Madeira_subset_maf_trim_0.03_Plinkrel0.2.fst", header = T,stringsAsFactors = F)

Locus_Fst_MAD_mapped <- subset(Locus_Fst_MAD, CHR >0)
Locus_Fst_CAN_mapped <- subset(Locus_Fst_CAN, CHR >0)


# Eigenvecotor PCA --------------------------------------------------------

egvec_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_Plinkrel0.2.eigenvec")
egvec_MAD <- read.table("Berthelots_Madeira_subset_maf_trim_0.03_Plinkrel0.2.eigenvec")


# Eigenvalues Manhatton Plot ----------------------------------------------

eg1_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_Plinkrel0.2.1.egwas", header = T,stringsAsFactors = F)
eg2_CAN <- read.table("Berthelots_Canaries_subset_maf_trim_0.03_Plinkrel0.2.2.egwas", header = T,stringsAsFactors = F)

eg1_MAD <- read.table("Berthelots_Madeira_subset_maf_trim_0.03_Plinkrel0.2.1.egwas", header = T, stringsAsFactors = F)
eg2_MAD <- read.table("Berthelots_Madeira_subset_maf_trim_0.03_Plinkrel0.2.2.egwas", header = T, stringsAsFactors = F)


