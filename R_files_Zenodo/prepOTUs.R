## prepOTUs.R
# Processing of OTU tables for analysis
# Author: Jun Ying Lim
# Reference: Lim et al. Climatic niche conservatism shapes the ecological assembly of Hawaiian arthropod communities
# NOTE: remember to set the working directory to the folder containing all the scripts and data

## PACKAGES ===============
library(reshape2)

# CLEAN UP OTU DATA ====================
# Import taxonomic reference table -------
taxonData <- read.csv("taxonData.csv")

# Define function to collapse matrices at different levels (zOTU vs OTU) -------
collapseOTU <- function(x, ref, collapse.by, to.collapse, subset.by){
  # Collapse zOTU tables into specified levels
  #
  # Arguments:
  #     x = matrix, e.g., zotu table to be collapsed; sites as columns, otus as rows
  #     ref = data.frame, reference table by which otu table
  #     collapse.by = string, specifies column in ref data.frame IDs to collapse by in input matrix, x (e.g., species)
  #     to.collapse = string, specifies column in ref data.frame IDs to be collapsed in input matrix, x (e.g., zOTU)
  #     subset.by = string, specifies column in ref data.frame IDs to be subset. If specified, output will be a list of matrices
  # Returns:
  #    matrix, or list of matrices, that have been collapsed by specified parameters
  
  if(missing(collapse.by)){
    stop("Missing argument: 'collapse.by' ")
  } else if(!collapse.by %in% names(ref)){
    stop("collapse.by not found in reference dataframe, ref")
  }
  if(missing(to.collapse)){
    stop("Missing argument: 'to.collapse' ")
  } else if(!to.collapse %in% names(ref)) {
    stop("to.collapse not found in reference dataframe, ref")
  }
  if(!missing(subset.by)){
    if(!subset.by %in% names(ref)){
      stop("subset.by not found in reference dataframe, ref")
    }
  }
  
  collapseByList <- unique(ref[[collapse.by]])
  resList <- list()
  for(i in 1:length(collapseByList)){
    tocollapseList <- ref[ref[collapse.by]==collapseByList[i],][[to.collapse]]
    if(length(tocollapseList) == 1){
      resList[[i]] <- x[rownames(x) %in% tocollapseList,]
    } else {
      resList[[i]] <- colSums(x[rownames(x) %in% tocollapseList,])
    }
  }
  res <- do.call("rbind", resList)
  rownames(res) <- collapseByList
  
  if(missing(subset.by)){
    return(res)
  } else{
    subsetres <- list()
    subsetList <- unique(as.vector(ref[[subset.by]]))
    for(j in 1:length(subsetList)){
      subsetCollapsed <- unique(ref[ref[subset.by]==subsetList[j],][[collapse.by]])
      subsetres[[subsetList[j]]] <- res[subsetCollapsed,]
    }
    return(subsetres)
  }
}

# Define some output folders -------
#masterOTU.dir <- paste0("data/", "masterOTU") 
OTU_native.dir <- paste0("data/", "OTU_native") #
zOTU_native.dir <- paste0("data/", "zOTU_native") # For 
OTU_native_order.dir <- paste0("data/", "OTU_native_order") # For each order separatel

# Create output folders
#if(!dir.exists(masterOTU.dir)){ dir.create(masterOTU.dir) }
if(!dir.exists(OTU_native.dir)){ dir.create(OTU_native.dir) }
if(!dir.exists(zOTU_native.dir)){ dir.create(zOTU_native.dir) }
if(!dir.exists(OTU_native_order.dir)){ dir.create(OTU_native_order.dir) }


otufiles <- list.files("data")[grep(list.files("data"), pattern = "combinedOTUdata_r")]
orderList <- c("Araneae", "Hemiptera", "Lepidoptera", "Psocoptera", "Coleoptera", "Orthoptera")

for(i in 1:100){
  print(i)
  
  # Create master zOTU table
  combData <- readRDS(paste0("data/", otufiles[i])) # import master files
  masterTable <- acast(zOTU_ID ~Site_ID, data = combData, fun.aggregate = sum, value.var = "rarefiedReadAbund")
  #saveRDS(masterTable, file.path(masterOTU.dir, paste0("master_zOTU_r", i, ".rds")))
  
  # Collapse ZOTU matrix into OTU and Species matrices
  taxonData <- taxonData[taxonData$zOTU_ID %in% combData$zOTU_ID,]
  OTUtab <- collapseOTU(x = masterTable, ref = taxonData, collapse.by = "OTU_ID", to.collapse = "zOTU_ID")
  
  # Exclude the following orders 
  toExclude <- subset(taxonData, TaxonomicOrder %in% c("Collembola", "Chilopoda", "Diplopoda", "Blattodea", "Isopoda", "Mantodea", "Phasmatodea"))
  masterTable_native <- masterTable[!rownames(masterTable) %in% toExclude$zOTU_ID, ]
  saveRDS(masterTable_native, file.path(zOTU_native.dir, paste0("zOTU_native_r", i, ".rds")))
  
  OTUtab_native <- OTUtab[!rownames(OTUtab) %in% toExclude$OTU_ID, ]
  saveRDS(OTUtab_native, file.path(OTU_native.dir, paste0("OTU_native_r", i, ".rds")))
  
  # Subset individual orders
  for(j in orderList){
    targetTaxon <- subset(taxonData, TaxonomicOrder == j)
    masterTable_subset <- masterTable_native[rownames(masterTable_native) %in% targetTaxon$zOTU_ID, ]
    OTUtab_subset <- OTUtab_native[rownames(OTUtab_native) %in% targetTaxon$OTU_ID, ]
    saveRDS(masterTable_subset, file.path(OTU_native_order.dir, paste0("zOTU_", j, "_r", i, ".rds")))
    saveRDS(OTUtab_subset, file.path(OTU_native_order.dir, paste0("OTU_", j, "_r", i,".rds")))
  }
}

# SUBSET RAREFIED DATA TABLES FOR RESISTANCE ANALYSIS ================
# select rarefied datasets to use
set.seed(12345)
otu_subsets <- sample(1:100, size = 5, replace = FALSE)
# 14, 51, 80, 90, 92

# create a new directory to store subsampled data
resistance.dir <- paste0("data/", "resistance_analysis")
if(!dir.exists(resistance.dir)){dir.create(resistance.dir)}

temp <- list.files(OTU_native_order.dir, pattern = paste(otu_subsets, collapse = "|"))
resistance.fname <- temp[grepl(x = temp, pattern = "^OTU")]

resistance_tables <- lapply(paste0(OTU_native_order.dir, "/", resistance.fname), FUN = readRDS)

siteData <- read.csv("siteData.csv")
laupahoehoe_siteIDs <- subset(siteData, site == "Laupahoehoe")$site.id
stainback_siteIDs <-subset(siteData, site == "Stainback")$site.id

for(i in 1:length(resistance_tables)){
  saveRDS(resistance_tables[[i]][,as.character(laupahoehoe_siteIDs)],
          file = file.path(resistance.dir, paste0("laupahoehoe_", resistance.fname[i])))
  saveRDS(resistance_tables[[i]][,as.character(stainback_siteIDs)],
          file = file.path(resistance.dir, paste0("stainback_", resistance.fname[i])))
}