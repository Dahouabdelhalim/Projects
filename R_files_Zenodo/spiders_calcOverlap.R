## SPIDER ANALYSIS: Hypervolume overlaps
# Written by: Jun Ying Lim and Susan Kennedy
# From publication:
# Kennedy, S. R., Lim, J. Y., Clavel, J., Krehenwinkel, H., and Gillespie, R. G. 2019. Spider webs, stable isotopes and molecular gut content analysis: Multiple lines of evidence support trophic niche differentiation in a community of Hawaiian spiders. Functional Ecology.

## DIRECTORIES ==============================
#Set "main.dir" to your own directory containing subfolders:
main.dir <- setwd() #Insert your file path in parentheses
data.dir <- file.path(main.dir, "raw_data")
res.dir <- file.path(main.dir, "results")

## PACKAGES ============================== 
library(vegan)
library(hypervolume)
source(file.path(main.dir, "hypervolumeTools.R"))

## INPUT FILE ==============================
webData <- read.csv(file.path(data.dir, "webs_2013_2014.csv"), header = TRUE)
isotopeData <- read.csv(file.path(data.dir, "isotopes_2014_2016.csv"), header = TRUE)
gutDataES <- read.csv(file.path(data.dir, "gutES.csv"), header = TRUE) # spider gut content data, exactly equal rep per species


## CLEAN AND SUBSET DATA ==============================
spListWeb <- c("acuta", "eurychasma", "filiciphilia", "stelarobusta", "trituberculata")
spListSpiny <- c("brevignatha", "kamakou", "quasimodo", "waikamoi")

webData <- subset(webData, species %in% spListWeb) # only include the following species

gutDataWebES <- subset(gutDataES, species %in% spListWeb)
gutDataSpinyES <- subset(gutDataES, species %in% spListSpiny)
write.csv(gutDataWebES, file.path(res.dir, "gutDataWebES.csv"))
write.csv(gutDataSpinyES, file.path(res.dir, "gutDataSpinyES.csv"))


isotopeData_Web <- subset(isotopeData, species %in% spListWeb)
isotopeData_Spiny <- subset(isotopeData, species %in% spListSpiny)
write.csv(isotopeData_Web, file.path(res.dir, "isotopeData_Web.csv"))
write.csv(isotopeData_Spiny, file.path(res.dir, "isotopeData_Spiny.csv"))

# Factorize or numericize some fields
webData$species<-as.factor(as.vector(webData$species))
webData$UID <- as.factor(webData$UID)

# Log transforming some variables
webData$log.CA <- log(webData$CA)
webData$log.CTL <- log(webData$CTL)

# Change row names
rownames(webData) <- webData$UID

# Define some column variables
site.col <- c("angle", "height")
web.col <- c("log.CA", "log.CTL", "MW", "SD1", "SD2","radii", "rows")
veg.col <- c("attveg1", "attveg1.group")
isotope.col <- c("percent.N", "delta.15.N", "percent.C", "delta.13.C")

gut.col <- c("Acari", "Araneae", "Coleoptera", "Collembola", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Psocoptera")

gut.colWeb <- c("Acari", "Araneae", "Collembola", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Psocoptera")

gut.colSpiny <- c("Acari", "Araneae", "Coleoptera", "Collembola", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera")

## PRINCIPAL COMPONENT ANALYSIS ==============================
# Principal component analysis with web.col

webPCA <- prcomp(~ log.CA + log.CTL + MW + SD1 + SD2 + radii + rows, data = webData, center = TRUE, scale. = TRUE)
webData[paste0("PCA", 1:7)] <- webPCA$x[,1:7]

webPCAaxes <- data.frame(webPCA$rotation)
webPCAaxes$label <- rownames(webPCAaxes)
summary(webPCA)
write.csv(webPCAaxes, file.path(res.dir, "webPCAaxes.csv"), row.names = FALSE)
write.csv(webData, file.path(res.dir, "webAnalysis.csv"), row.names = FALSE)

webPCAsummary <- rbind(webPCA$sdev^2, summary(webPCA)$importance[2,], webPCA$rotation)
webPCAsummary <- round(webPCAsummary, digits = 3)
write.csv(webPCAsummary, file.path(res.dir, "webPCAsummary.csv"), row.names = TRUE)

## CALCULATE HYPERVOLUMES ==============================

# Calculate web hypervolumes and pairwise hypervolume overlaps using ordination loads
webHypervol <- generateHypervolume(data = webData[, c("species", paste0("PCA", 1:2))], id.col = "species")
webHypervolOverlapJ <- pairwiseHypervolumeOverlap(webHypervol, method = "jaccard") #Jaccard

# Export hypervolume results
saveRDS(webHypervol, file.path(res.dir, "webHypervol.rds"))
saveRDS(webHypervolOverlapJ, file.path(res.dir, "webHypervolOverlapJ.rds"))

# Generate isotope hypervolumes and pairwise hypervolume overlaps
isotopeHypervolWeb <- generateHypervolume(data = isotopeData_Web[c("species", "delta.15.N", "delta.13.C")], id.col = "species")
isotopeHypervolOverlapWebJ <- pairwiseHypervolumeOverlap(isotopeHypervolWeb, method = "jaccard")

saveRDS(isotopeHypervolWeb, file.path(res.dir, "isotopeHypervolWeb.rds"))
saveRDS(isotopeHypervolOverlapWebJ, file.path(res.dir, "isotopeHypervolOverlapWebJ.rds"))

isotopeHypervolSpiny <- generateHypervolume(data = isotopeData_Spiny[c("species", "delta.15.N", "delta.13.C")], id.col = "species")
isotopeHypervolOverlapSpinyJ <- pairwiseHypervolumeOverlap(isotopeHypervolSpiny, method = "jaccard")

saveRDS(isotopeHypervolSpiny, file.path(res.dir, "isotopeHypervolSpiny.rds"))
saveRDS(isotopeHypervolOverlapSpinyJ, file.path(res.dir, "isotopeHypervolOverlapSpinyJ.rds"))

## CALCUALTE GUT BETA DIVERSITY

gutBetaWebES <- as.matrix(vegdist(x = gutDataWebES[gut.colWeb], method = "bray", binary = TRUE))

gutMeanBetaWebES <- matrix(nrow = length(spListWeb), ncol = length(spListWeb))
for(i in 1:length(spListWeb)){
  rowids <- rownames(subset(gutDataWebES, species == spListWeb[i]))
  for(j in 1:length(spListWeb)){
    colids <- rownames(subset(gutDataWebES, species == spListWeb[j]))
    gutMeanBetaWebES[i,j] <- mean(gutBetaWebES[rowids, colids])
  }
}
rownames(gutMeanBetaWebES) <- colnames(gutMeanBetaWebES) <- spListWeb
saveRDS(gutMeanBetaWebES, file.path(res.dir, "gutBetaWebES.rds"))

gutBetaSpinyES <- as.matrix(vegdist(x = gutDataSpinyES[gut.colSpiny], method = "bray", binary = TRUE))
gutMeanBetaSpinyES <- matrix(nrow = length(spListSpiny), ncol = length(spListSpiny))
for(i in 1:length(spListSpiny)){
  rowids <- rownames(subset(gutDataSpinyES, species == spListSpiny[i]))
  for(j in 1:length(spListSpiny)){
    colids <- rownames(subset(gutDataSpinyES, species == spListSpiny[j]))
    gutMeanBetaSpinyES[i,j] <- mean(gutBetaSpinyES[rowids, colids])
  }
}
rownames(gutMeanBetaSpinyES) <- colnames(gutMeanBetaSpinyES) <- spListSpiny
saveRDS(gutMeanBetaSpinyES, file.path(res.dir, "gutBetaSpinyES.rds"))

