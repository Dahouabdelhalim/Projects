## generateRarefiedMatrices.R
# Rarefies OTU tables for analysis
# Author: Jun Ying Lim
# Reference: Lim et al. Climatic niche conservatism shapes the ecological assembly of Hawaiian arthropod communities
# NOTE: remember to set the working directory to the folder containing all the scripts and data

## DIRECTORIES ===============
if(!dir.exists("data")){ dir.create("data") }

## PACKAGES ===============
library(stringr)
library(plyr)
library(reshape2)
library(vegan)

## IMPORT DATA ===============
# Import site and otu tables -------
siteData <- read.csv("siteData.csv")
dataARF <- read.delim("otuARF_clean.txt", header = TRUE, stringsAsFactors = FALSE)
dataMCO <- read.delim("otuMCO_clean.txt", header = TRUE, stringsAsFactors = FALSE)
names(dataARF)[1] <- "X.OTU.ID"
names(dataMCO)[1] <- "X.OTU.ID"

# Import specimen count table -------
specimenCounts <- read.delim("specimennumber.txt", header = TRUE)
names(specimenCounts) <- c("Site_ID", "Transect", "Island", "Year", "SizeCategory", "Count")

## DATA PREPARATION ===============
# Convert matrix into long-form data frame -------
dataARFmelt <- melt(dataARF, id.vars = "X.OTU.ID", value.name = "nReads", variable.name = "Site_SizeCategory")
dataMCOmelt <- melt(dataMCO, id.vars = "X.OTU.ID", value.name = "nReads", variable.name = "Site_SizeCategory")

# Clean up data -------
cleanOTUdata <- function(x){
  # Label OTUs by site and size categories
  temp <- str_split_fixed(x$Site_SizeCategory, pattern = "_", n = 3)
  x$Site_ID <- gsub(temp[,1], pattern = "X", replacement = "")
  x$SizeCategory <- temp[,2]
  temp <- str_split_fixed(x$X.OTU.ID, pattern = "_", n = 3)
  x$zOTU_ID <- temp[,1]
  x$OTU_ID <- temp[,2]
  return(x)
}

arfOTU <- cleanOTUdata(dataARFmelt)
mcoOTU <- cleanOTUdata(dataMCOmelt)

## RAREFY DATA BY SIZE CATEGORY AND SITE ===============

# Define rarefaction function -------
rarefyOTUbySpecimenCount <- function(df, counts, readsPerIndividual){
  # Rarefies category by specimen counts. This function is designed to be iterated using `ddply`
  #
  # Arguments:
  #   df: data.frame subset parsed from ddply
  #   counts: data.frame containing specimen counts
  #   readsPerIndividual: number of reads per individual you would like to rarefy with
  #   nReps: number of randomizations
  #
  # Returns:
  #   data,frame co
  
  # Extract the number of specimens for site and size category
  nSpecimen <- subset(counts, Site_ID == unique(df$Site_ID) &
                        SizeCategory == unique(df$SizeCategory))$Count
  
  rawReadAbund <- acast(df, X.OTU.ID~Site_SizeCategory, value.var = "nReads")
  
  rarefiedReadAbund <- rrarefy(rawReadAbund, sample = nSpecimen*readsPerIndividual) 
  temp <- melt(rarefiedReadAbund,
               value.name = "rarefiedReadAbund",
               varnames = c("Site_SizeCategory", "X.OTU.ID"))
  
  res <- merge(df, temp, by = c("X.OTU.ID", "Site_SizeCategory"))
  return(res)
}

# Rarefy datasets 
# Use the number of specimen counts in each size category to rarefy
# ARF and MCO are rarefied separately and then combined

nreps <- 100
arf_rarefied <- list()
mco_rarefied <- list()

arfOTU_summary <- ddply(.data =  arfOTU, .var = .(SizeCategory,Site_ID),.fun = summarize, arfTotalReads = sum(nReads))
mcoOTU_summary <- ddply(.data =  mcoOTU, .var = .(SizeCategory,Site_ID),.fun = summarize, mcoTotalReads = sum(nReads))

arfOTU_summary2 <- merge(arfOTU_summary, specimenCounts)
mcoOTU_summary2 <- merge(mcoOTU_summary, specimenCounts)

sum(mcoOTU_summary2$mcoTotalReads < (mcoOTU_summary2$Count * 7)) # all samples have enough reads for rarefaction
sum(mcoOTU_summary2$mcoTotalReads < (mcoOTU_summary2$Count * 8)) # all samples have enough reads for rarefaction

sum(arfOTU_summary2$arfTotalReads < (arfOTU_summary2$Count * 7)) # all samples have enough reads for rarefaction
sum(arfOTU_summary2$arfTotalReads < (arfOTU_summary2$Count * 8)) # 4 samples have too few reads for rarefaction at the same level as other samples

for(i in 1:nreps){
  pb = txtProgressBar(min = 0, max = nreps, style = 3, initial = 0) 
  setTxtProgressBar(pb, i)
  arf_rarefied[[i]] <- ddply(.data = arfOTU,
                             .variables = .(SizeCategory, Site_ID),
                             .fun = rarefyOTUbySpecimenCount,
                             counts = specimenCounts,
                             readsPerIndividual = 7)
                             
  mco_rarefied[[i]] <- ddply(.data = mcoOTU,
                             .variables = .(SizeCategory,Site_ID),
                             .fun = rarefyOTUbySpecimenCount,
                             counts = specimenCounts,
                             readsPerIndividual = 7)
}
close(pb)

saveRDS(arf_rarefied, file = "data/arf_rarefied.rds")
saveRDS(mco_rarefied, file = "data/mco_rarefied.rds")

## SUM OTU READ ABUNDANCES BY SITE ===============
for(i in 1:nreps){
  pb = txtProgressBar(min = 0, max = nreps, style = 3, initial = 0) 
  setTxtProgressBar(pb, i)
  temp <- rbind(arf_rarefied[[i]],
                mco_rarefied[[i]])
  
  saveRDS(cleanOTUdata(temp), paste0("data/combinedOTUdata_r", i, ".rds"))
}
close(pb)