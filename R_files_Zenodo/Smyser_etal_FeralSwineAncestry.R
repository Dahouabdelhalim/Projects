# Note to users: You will need to first enter the path to your executable PLINK and ADMIXTURE
# files stored locally and the path to the working directory where data files are stored

# map to PLINK and ADMIXTURE program files
plinkPath <- "/home/tsmyser/bin/plink"
admixturePath <- "/home/tsmyser/bin/admixture/admixture"

# set working directory
setwd("/home/tsmyser/Documents/Smyser_FeralSwineAncestry")

# create "results" folder within working directory to compile results
call <- "mkdir results"
system(call)

# read in Feral Swine genotypes as *.bed, *.bim, *.fam file family to query against Reference Set
# Family ID (column 1 in *.fam file) must be '-' in order to conduct a supervised analysis in
# which feral swine of unknown ancestry are blindly queried against Reference Set
# 
FeralSwine <- "FeralSwine_ContinentalUS" 
# create "FeralSwine" as an R object for convenience in calling file later in code 
FeralSwine_fam <- read.table("FeralSwine_ContinentalUS.fam")

# specify column names for working in R as column names are not included in *.fam files
names(FeralSwine_fam) <- c("FamilyID", "PatientID", "FatherID", "MotherID", "Sex", "AffectionStatus")
head(FeralSwine_fam)

# read in Reference Set [domestic pigs and wild boar] with the defined reference cluster
# for each individual presented in the Family ID column (column 1 in *.fam file).
# Genotypes of Reference Set should be ordered based on the described reference clusters, as the 
# order of the reference clusters will correspond with the order of the columns for the ancestry
# vector for queried feral swine

RefSet <- "Sus_scrofa_ReferenceSet"
# create "RefSet" as an R object for convenience in calling file later in code
RefSet_fam <- read.table("Sus_scrofa_ReferenceSet.fam")
# specify column names as above
names(RefSet_fam) <- c("FamilyID", "PatientID", "FatherID", "MotherID", "Sex", "AffectionStatus")
head(RefSet_fam)
table(RefSet_fam$FamilyID)

###########################################################
# initialize empty data frames to compile results for iterative analyses 

# for this analysis, we are querying feral swine relative to k = 17 reference clusters

k <- length(unique(RefSet_fam$FamilyID))

SampleIDs <- as.data.frame(matrix(NA, nrow(FeralSwine_fam), 2))
AncestryResults_Q <- as.data.frame(matrix(0, nrow(FeralSwine_fam), k))
AncestryResults_SE <- as.data.frame(matrix(0, nrow(FeralSwine_fam), k))
AncestryResults_bias <- as.data.frame(matrix(0, nrow(FeralSwine_fam), k))

############################################################################################
# Beginning of first analysis to generate bootstrapped ancestry estimates of feral swine
# relative to 17 reference clusters for Sus scrofa
############################################################################################

for (i in 1:nrow(FeralSwine_fam)){
  
  # identify individual feral swine to be queried in iteration i
  # retain individual feral swine as QueryInd.ped, QueryInd.map
  Keep <- FeralSwine_fam[i,1:2]
  write.table(Keep, "KeepList", row.names = F, col.names = F, quote = F)
  call <- paste(plinkPath, "--bfile", FeralSwine, "--keep KeepList --recode --out QueryInd")
  system(call)
  
  # merge individual feral swine genotype (QueryInd.ped) with Reference Set genotypes for 
  # supervised analysis
  call <- paste(plinkPath, "--bfile", RefSet, "--merge QueryInd.ped QueryInd.map --make-bed --out Iteration")
  system(call)
  
  # ADMIXTURE requires a *.pop file to complement *.bed/*.bim/*.fam files in supervised analyses    
  pop <- read.table("Iteration.fam")[,1:2]
  write.table(pop, "Iteration.pop", quote = F, row.names = F, col.names = F)
  
  # conduct supervised analysis, querying the individual feral swine (QueryInd.ped) against 
  # Reference Set with desired number of bootstrapped iterations (option -B##) and number of 
  # threads to be allocated (-j#)
  # 
  call <- paste(admixturePath, "-j4 -B100 --supervised Iteration.bed", k) 
  system(call)
  
  # compile results files.  Given that feral swine, of unknown ancestry, have a Family ID of '-'
  # the QueryInd will be appended to top (row 1) of the Reference Set in the formation of the 
  # Iteration.bed/Iteration.bim/Iteration.fam.  Accordingly, the ancestry results for QueryInd
  # will also appear in row 1 of the ensuing Iteration.17.Q, Iteration.17.Q_se, and 
  # Iteration.17.Q_bias results files.
  SampleIDs[i,1] <- read.table("Iteration.fam", stringsAsFactors = F)[1,2]
  SampleIDs[i,2] <- read.table("Iteration.fam", stringsAsFactors = F)[1,1]
  
  AncestryResults_Q[i, ] <- read.table("Iteration.17.Q")[1,]
  AncestryResults_SE[i, ] <- read.table("Iteration.17.Q_se")[1,]
  AncestryResults_bias[i, ] <- read.table("Iteration.17.Q_bias")[1,]
  
  FeralSwine_Q <- cbind(SampleIDs, AncestryResults_Q)
  write.csv(FeralSwine_Q, "./results/FeralSwine_K17Ancestry_BootStrap_Q.csv", quote = FALSE)
  
  FeralSwine_SE <- cbind(SampleIDs, AncestryResults_SE)
  write.csv(FeralSwine_SE, "./results/FeralSwine_K17Ancestry_BootStrap_SE.csv", quote = FALSE)
  
  FeralSwine_bias <- cbind(SampleIDs, AncestryResults_bias)
  write.csv(FeralSwine_bias, "./results/FeralSwine_K17Ancestry_BootStrap_bias.csv", quote = FALSE)
  
}

#########################################
#########################################
#########################################

# identify significant ancestry groups based on results of first analysis of 17 reference 
# clusters 
Q <- FeralSwine_Q[, 3:19]
SE <- FeralSwine_SE[, 3:19]

# evaluate whether the confidence interval of the ancestry association for feral swine individual 'i' 
# overlaps with zero* - indicating a non-significant association
# *here we use 1e-05 in place of zero as this is the lowest association returned from ADMIXTURE

Zero <- Q-1.96*SE # assuming normal distribution of Q-values

Significance <- matrix(NA, nrow = dim(Q)[1], 17) # a table indicating significant Q-scores"

for (i in 1:dim(Q)[1]){
  for (j in 1:17){
    if(Zero[i,j] > 1e-05) {
      Significance[i,j] = 1} 
    else if (Zero[i,j] <= 1e-05) {
      Significance[i,j] <- 0 }
    else {Significance[i,j] <- NA}
  }}

Significance <- as.data.frame(Significance)
rownames(Significance) <- FeralSwine_Q$V1

write.csv(Significance, "FeralSwine_K17Ancestry_SignificanceMatrix.csv", quote = F, row.names = T)

# Identify and remove feral swine with significant association to a single reference cluster as subsequent 
# formation of 'an individually customized reference set composed of a subset of the 17 reference clusters'
# is uninformative for such individuals.
Significance$ClusterCount <-rowSums(Significance)
Purebred <- subset(Significance, Significance$ClusterCount == 1)
Purebred <- Purebred[, 1:17]
Admixed <- subset(Significance, Significance$ClusterCount >= 2)
Admixed <-Admixed[, 1:17]

if(nrow(Purebred) >0) {
  write.csv(Purebred, "Purebred_SigMatrix.csv", quote = F, row.names = T )
}
if(nrow(Admixed) >0) {
  write.csv(Admixed, "Admixed_SigMatrix.csv", quote = F, row.names = T )
}

rm(Purebred)

# populate matrix of reference clusters to be uniquely retained for the customized reference
# set for individual feral swine 'i'

m1 <- matrix(seq(1:k), k, nrow(Admixed))
m1 <- t(m1)
Combined_Ind_Ref <- Admixed*m1

###################################################################################
# Remove 'purebred' feral swine from genotype file, carrying forward feral swine 
# of admixed ancestry for analysis with individually customized reference set

Admixed_KeepList <- as.data.frame(matrix(NA, nrow(Admixed), 2))
Admixed_KeepList[,1] <- "-"
Admixed_KeepList[,2] <- rownames(Admixed)

write.table(Admixed_KeepList, "Admixed_KeepList", row.names = F, col.names = F, quote = F)

call <- paste(plinkPath,  " --bfile ", FeralSwine , " --keep Admixed_KeepList --make-bed --out Admixed_IndRefSet_", nrow(Admixed), "x29375", sep = "") 
system(call)

# create results files for individual reference set analysis
IndRefSet_SampleIDs <- as.data.frame(matrix(NA, nrow(Admixed_KeepList), 2))
IndRefSet_Results_Q <- as.data.frame(matrix(0, nrow(Admixed_KeepList), k))
IndRefSet_Results_SE <- as.data.frame(matrix(0, nrow(Admixed_KeepList), k))
IndRefSet_Results_bias <- as.data.frame(matrix(0, nrow(Admixed_KeepList), k))

############################################################################################
# Beginning of second analysis to generate bootstrapped ancestry estimates of feral swine
# relative to individually customized reference sets
############################################################################################

for (i in 1:nrow(Admixed_KeepList)){
  
  # define individual reference set
  Ind_ref <- Combined_Ind_Ref[i,]
  Ind_ref <- Ind_ref[ Ind_ref != 0 ]
  
  # Using 'keep' command in plink to only retain significant reference clusters for feral
  # swine 'i'
  
  Keep <- RefSet_fam[which(RefSet_fam$FamilyID == Ind_ref[1]),1:2]
  write.table(Keep, "Keep", col.names = F, row.names = F, quote = F)
  
  call <- paste(plinkPath, "--bfile", RefSet, "--keep Keep --make-bed --out Ind_RefSet") 
  system(call)
  
  
  for(r in 2:length(Ind_ref)){
    
    Keep <- RefSet_fam[which(RefSet_fam$FamilyID == Ind_ref[r]),1:2]
    write.table(Keep, "Keep", col.names = F, row.names = F, quote = F)
    
    call <- paste(plinkPath, "--bfile", RefSet, "--keep Keep --recode --out temp") 
    system(call)
    
    call <- paste(plinkPath, "--bfile Ind_RefSet --merge temp.ped temp.map --make-bed --out Ind_RefSet") 
    system(call)
  }
  
  # with individually customized reference set now assembled, pull out individual feral swine
  # genotype 'i' for analysis 
  
  Keep <- Admixed_KeepList[i,1:2]
  write.table(Keep, "Keep", col.names = F, row.names = F, quote = F)
  
  call <- paste(plinkPath, " --bfile Admixed_IndRefSet_", nrow(Admixed), "x29375 --keep Keep --recode --out FeSw_Ind", sep = "") 
  system(call)
  
  call <- paste(plinkPath, "--bfile Ind_RefSet --merge FeSw_Ind.ped FeSw_Ind.map --make-bed --out FeSw_Query")
  system(call)
  
  FeSw_Query <- read.table("FeSw_Query.fam")
  write.table(FeSw_Query[,1:2], "FeSw_Query.pop", col.names = F, row.names = F, quote = F)
  
  # initiate second supervised ADMIXTURE analysis with -BXXX bootstrap iterations 
  call <- paste(admixturePath, "-j4 -B100 --supervised FeSw_Query.bed ", length(unique(Ind_ref)))
  system(call)
  
  
  # compile results from individual reference set analysis and compile ancestry associates in
  # the appropriate columns
  
  Ind_Results <- as.numeric(read.table(paste("FeSw_Query.", length(unique(Ind_ref)), ".Q", sep = "" ))[1,])
  for (q in 1:length(unique(Ind_ref))){
    Column_Number <- Ind_ref[q]
    IndRefSet_Results_Q[i,Column_Number] <- Ind_Results[q]
  }
  
  Ind_SE <- as.numeric(read.table(paste("FeSw_Query.", length(unique(Ind_ref)), ".Q_se", sep = "" ))[1,])
  for (s in 1:length(unique(Ind_ref))){
    Column_Number <- Ind_ref[s]
    IndRefSet_Results_SE[i,Column_Number] <- Ind_SE[s]
  }
  
  Ind_BIAS <- as.numeric(read.table(paste("FeSw_Query.", length(unique(Ind_ref)), ".Q_bias", sep = "" ))[1,])
  for (b in 1:length(unique(Ind_ref))){
    Column_Number <- Ind_ref[b]
    IndRefSet_Results_bias[i,Column_Number] <- Ind_BIAS[b]
  }
  
  IndRefSet_SampleIDs[i,1] <- read.table("FeSw_Query.fam", stringsAsFactors = F)[1,2]
  IndRefSet_SampleIDs[i,2] <- read.table("FeSw_Query.fam", stringsAsFactors = F)[1,1]
  
  write.csv(cbind(IndRefSet_SampleIDs, IndRefSet_Results_Q), "./results/FeralSwine_Admixed_IndRefSet_Boostrap_Q.csv", quote = F)
  write.csv(cbind(IndRefSet_SampleIDs, IndRefSet_Results_SE), "./results/FeralSwine_Admixed_IndRefSet_Boostrap_SE.csv", quote = F)
  write.csv(cbind(IndRefSet_SampleIDs, IndRefSet_Results_bias), "./results/FeralSwine_Admixed_IndRefSet_Boostrap_bias.csv", quote = F)
}

