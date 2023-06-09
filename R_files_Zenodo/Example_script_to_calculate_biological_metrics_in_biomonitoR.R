# install the devtools package and install the old version of biomonitoR
library(devtools)
install_github("alexology/biomonitoR", ref = "old_version")

# load the biomonitoR package
library(biomonitoR)


# load the taxa by samples file; we call it "taxaxsamples"
taxaxsamples <-  read.csv("insert_file_name.csv" , h = TRUE)

taxaxsamples[is.na(taxaxsamples)] <- 0

#convert to presence absence
taxaxsamples <- taxaxsamples
taxaxsamples[taxaxsamples > 0] <- 1

traits <- read.csv("traits.csv", h = TRUE)
#"traits" is where we look up the value each taxon gets for each trait

ref_def <-  read.csv("ref.csv", h = TRUE)
#"ref" is where we look up the taxon names

#remove sites with 0 taxa and species with 0 abundance
macro_tempG <- taxaxsamples[, -1]

macroG <- data.frame(taxaxsamples[rowSums(macro_tempG) > 0, 1,  drop = FALSE], macro_tempG[rowSums(macro_tempG) > 0,  colSums(macro_tempG) > 0])

macro_bioG <- asBiomonitor(macroG, dfref = ref_def)

macro_aggG <- aggregatoR(macro_bioG)

#richness
taxa_rich <- richness(macro_aggG, taxLev = "Taxa")
taxa_rich <- as.data.frame(taxa_rich)
colnames(taxa_rich) <- "Taxa"

fam_rich <- richness(macro_aggG, taxLev = "Family")
fam_rich <- as.data.frame(fam_rich)
colnames(fam_rich) <- "Families"
allind <- as.data.frame(merge(taxa_rich, fam_rich,by = "row.names"))

#BMWP
BMWP <- bmwp(macro_aggG, method = "uk")
BMWP <- as.data.frame(BMWP)
colnames(BMWP) <- "BMWP"
allind2 <- as.data.frame(merge(allind, BMWP, by.y = "row.names", by.x = "Row.names"))

#ASPT
ASPT <- aspt(macro_aggG, method = "uk")
ASPT <- as.data.frame(ASPT)
colnames(ASPT) <- "ASPT"
allind3 <- as.data.frame(merge(allind2, ASPT, by.y = "row.names", by.x = "Row.names"))

#IBMWP
IBMWP <- bmwp(macro_aggG, method = "spa")
IBMWP <- as.data.frame(IBMWP)
colnames(IBMWP) <- "IBMWP"
allind4 <- as.data.frame(merge(allind3, IBMWP, by.y = "row.names", by.x = "Row.names"))

#IASPT
IASPT <- aspt(macro_aggG, method = "spa")
IASPT <- as.data.frame(IASPT)
colnames(IASPT) <- "IASPT"
allind5 <- as.data.frame(merge(allind4, IASPT, by.y = "row.names", by.x = "Row.names"))

#WHPT
WHPT_Total <- whpt(macro_aggG, taxLev = "Family", type = "po", method = "bmwp", composite = FALSE)
WHPT_Total <- as.data.frame(WHPT_Total)
colnames(WHPT_Total) <- "WHPT"
allind6 <- as.data.frame(merge(allind5, WHPT_Total, by.y = "row.names", by.x = "Row.names"))

#WHPT ASPT
WHPT_ASPT <- whpt(macro_aggG, taxLev = "Family", type = "po", method = "aspt", composite = FALSE)
WHPT_ASPT <- as.data.frame(WHPT_ASPT)
colnames(WHPT_ASPT) <- "WHPT_ASPT"
allind7 <- as.data.frame(merge(allind6, WHPT_ASPT, by.y = "row.names", by.x = "Row.names"))

#need to create "traitDB = data_ts_avF" for functional (trait-based) metrics to work.
names(traits)[1] <- "Taxa"
data_tsF <- traitScaling(macro_aggG, traitDB = traits, dfref = ref_def)
data_ts_avF <- aggregate(. ~ Taxa, data = data_tsF[, -c(2:5)], FUN = mean)
CWMF<- as.data.frame(cwm(macro_aggG))
CMWF1<- as.data.frame(CWMF[, c(1:114)])
write.table(CMWF1, file="CMWF2.xls", sep = "\\t")
CWMFsel<-as.data.frame(CWMF[,c(1, 9:39, 54:71, 79:89, 112:114)])
colB = c(9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8)

#functional richness
ffrich <- ffrich(macro_aggG, traitDB = data_ts_avF, colB = colB, traceB = TRUE)
ffrich <- as.data.frame(ffrich$results)
colnames(ffrich) <- "ffrich"
allind8 <- as.data.frame(merge(allind7, ffrich, by.y = "row.names", by.x = "Row.names"))

#functional redundancy
ffred_int <- ffred(macro_aggG,  traitDB = data_ts_avF, colB = colB, traceB = TRUE)
ffred <- ffred_int$results[, c(2, 3)]
colnames(ffred) <- "ffred"
allind9 <- as.data.frame(merge(allind8, ffred, by.y = "row.names", by.x = "Row.names"))

write.table(allind9, file="metrics.xls", sep = "\\t")



