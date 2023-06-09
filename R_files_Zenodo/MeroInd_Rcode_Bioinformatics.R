suppressPackageStartupMessages(library(mjolnir))

# Define input fastq files (only names of R1 files are needed)
R1_filenames <-c("IL20_S7_L001_R1_001.fastq","IL21_S2_L001_R1_001.fastq","IL22_S3_L001_R1_001.fastq",
"IL23_S4_L001_R1_001.fastq","IL24_S5_L001_R1_001.fastq","IL25_S6_L001_R1_001.fastq",
"PLAA_S1_L001_R1_001.fastq","PLAB_S2_L001_R1_001.fastq","PLAC_S3_L001_R1_001.fastq",
"PLAD_S4_L001_R1_001.fastq","PLAE_S5_L001_R1_001.fastq","PLAF_S6_L001_R1_001.fastq",
"PLAG_S7_L001_R1_001.fastq","PLAH_S8_L001_R1_001.fastq",
"PLRA_S1_L001_R1_001.fastq","PLRB_S2_L001_R1_001.fastq","PLRC_S3_L001_R1_001.fastq",
"PLRD_S4_L001_R1_001.fastq","PLRE_S5_L001_R1_001.fastq")

# Define number of cores to be used in parallel. For best performance, the number of libraries to process x the number of cores should be less or equal than the total cores available.
cores <- 5

# Input names of the individual libraries to be used. It should be a 4-character name, matching the information of the ngsfilter files
lib_prefixes <- c("IL20","IL21","IL22","IL23","IL24","IL25","PLAA","PLAB","PLAC","PLAD","PLAE","PLAF","PLAG","PLAH",
"PLRA","PLRB","PLRC","PLRD","PLRE")

# Input name for the final combined library (should be a 4-character name)
lib <- "PLIN"

####################
# MJOLNIR pipeline #
####################

# RAN will distribute R1 & R2 fastq files into equal-sized pieces for parallel computing
mjolnir1_RAN(R1_filenames,cores,lib_prefixes,R1_motif="_R1_",R2_motif="_R2_") 
# FREYJA will do the paired-end alignment, demultiplexing & length filtering
mjolnir2_FREYJA(lib_prefixes,cores,Lmin=299,Lmax=320)
# HELA will remove chimaeric sequences in a sample-by-sample basis, will change identifiers of remaining unique sequences & will generate a table of their abundances in each sample & a fasta file with unique sequences and their total abundance for ODIN
mjolnir3_HELA(lib_prefixes,lib,cores)
# Here we change the number of cores to full computing power. All libraries will be treated as a combined one from here on
cores <- 60
# ODIN will do the clustering & will generate a table with the abundances of each MOTU in each sample
mjolnir4_ODIN(lib,cores,d=13)
# THOR will asign the taxonomy to the representative sequence of each MOTU
mjolnir5_THOR(lib,cores,tax_dir="~/taxo_NCBI",ref_db="DUFA_Leray_20200610.fasta",taxo_db="taxo_NCBI") 
# FRIGGA will integrate the information of MOTU abundances and taxonomy assignment from ODIN & THOR in a single table
mjolnir6_FRIGGA(lib)
# LOKI kill remove the pseudogenes and will keep track of the taxonomic information of the removed MOTUs
mjolnir7_LOKI(lib,min_id=.84)
# RAGNAROC will change the names of the samples to recover the original names and will remove unnecessary columns
mjolnir8_RAGNAROC(lib,"PLAY_metadata.csv","PLAY_final_file.csv")
mjolnir8_RAGNAROC(lib,"PLRX_metadata.csv","PLRX_final_file.csv")
mjolnir8_RAGNAROC(lib,"ILXX_metadata.csv","ILXX_final_file.csv")

