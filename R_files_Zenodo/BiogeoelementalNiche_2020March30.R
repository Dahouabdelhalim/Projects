# Set directory
setwd("~/Dropbox (HeckfordT)/PhD_Thesis/Chapter_StoichNiche/")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
###### Process elemental data ######
# Load foliar elemental data from the lab
# North Pen data
StoiNorthPen <- read.csv("Data/OriginalData/NorthPen/Original_SL/NorthPen_foliarCNP.csv")
StoiNorthPen <- read.csv(file.choose())
# Central forest data
StoiCenForest <- read.csv("Data/OriginalData/TNNP/2016/CentralForest_FoliarCNP.csv")
StoiCenForest  <- read.csv(file.choose())
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Throughout code: 
# Ecoregion --> the between ecoregion comparisons
# Region/Regional --> across ecoregion comparisons
#Local --> within and between ecoregion comparions
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# Find conspecific and heterospecific plots
library(dplyr)
library(BBmisc)
library(prob)
# North Pen data
StoiNorthPen$PlotID <- toupper(StoiNorthPen$PlotID) # To fix error when grouping below with lower case plot letter
StoiNorthPen = StoiNorthPen %>% 
  group_by(PlotID) %>% # Goup by PlotID
  mutate(Hetero = n() > 1) # If name occurs more than once = TRUE
StoiNorthPen$Hetero <- as.character(StoiNorthPen$Hetero) # Set field as character so we can get "TRUE"
StoiNorthPen$LocalComType <- ifelse(StoiNorthPen$Hetero == "TRUE", "NP: Heterospecific", "NP: Conspecific") # Set hetero and cons field
StoiNorthPen$RegComType <- ifelse(StoiNorthPen$Hetero == "TRUE", "Heterospecific", "Conspecific") # Set hetero and cons field
StoiNorthPen$Region <- "Northern Peninsula"
StoiNorthPen <- as.data.frame(subset(StoiNorthPen, select=c("Species", "Region", "LocalComType", "RegComType", "PercentC", "PercentN", "PercentP")))
colnames(StoiNorthPen) <- c("Species", "Region", "LocalComType", "RegComType", "C", "N", "P") # Convert column names to match central forest data
str(StoiNorthPen) # Check structure
# Central forest data
StoiCenForest = StoiCenForest %>% 
  group_by(PlotID) %>% # Goup by PlotID
  mutate(Hetero = n() > 1) # If name occurs more than once = TRUE
StoiCenForest$Hetero <- as.character(StoiCenForest$Hetero) # Set field as character so we can get "TRUE"
StoiCenForest$LocalComType <- ifelse(StoiCenForest$Hetero == "TRUE", "CF: Heterospecific", "CF: Conspecific") # Set hetero and cons field
StoiCenForest$RegComType <- ifelse(StoiCenForest$Hetero == "TRUE", "Heterospecific", "Conspecific") # Set hetero and cons field
StoiCenForest$Region <- "Central Forest"
StoiCenForest <- as.data.frame(subset(StoiCenForest, select=c("Species", "Region", "LocalComType", "RegComType", "C", "N", "P")))
str(StoiCenForest) # Check structure

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Combine data
FoliarCNP <- rbind(StoiNorthPen, StoiCenForest)
FoliarCNP$LocalComType <- factor(FoliarCNP$LocalComType, levels = c("NP: Conspecific","NP: Heterospecific", "CF: Conspecific","CF: Heterospecific"))
FoliarCNP$Region <- as.factor(FoliarCNP$Region)
FoliarCNP$LocalComType <- as.factor(FoliarCNP$LocalComType)
FoliarCNP$RegComType <- as.factor(FoliarCNP$RegComType)
str(FoliarCNP)

# Subset by species
FoliarCNP_bF <- subset(FoliarCNP, FoliarCNP$Species == "bF") # n=390
FoliarCNP_wB <- subset(FoliarCNP, FoliarCNP$Species == "wB") # n=229

# Subset by ecoregion
bF_NP <- subset(FoliarCNP, Species == "bF" & Region == "Northern Peninsula" )
bF_CF <- subset(FoliarCNP, Species == "bF" & Region == "Central Forest" )
wB_NP <- subset(FoliarCNP, Species == "wB" & Region == "Northern Peninsula" )
wB_CF <- subset(FoliarCNP, Species == "wB" & Region == "Central Forest" )

###### Get mean and SE #####
# Subset data, get number of samples, get mean and standard error and difference between north south con and heterospecific
# Balsam fir
ABBA_RegionalCon <- subset(FoliarCNP_bF, RegComType == "Conspecific")
str(ABBA_RegionalCon)
ABBA_RegionalHetero <- subset(FoliarCNP_bF, RegComType == "Heterospecific")
str(ABBA_RegionalHetero)

# Make table of means between C, N, P of NP and CF conspecific and heterospecific niches
ABBA_LocalNP_Con <- subset(FoliarCNP_bF, LocalComType == "NP: Conspecific")
str(ABBA_LocalNP_Con)
ABBA_Con_NP_C <- paste0(round(mean(ABBA_LocalNP_Con$C), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalNP_Con$C), digits=3))
ABBA_Con_NP_N <- paste0(round(mean(ABBA_LocalNP_Con$N), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalNP_Con$N), digits=3))
ABBA_Con_NP_P <- paste0(round(mean(ABBA_LocalNP_Con$P), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalNP_Con$P), digits=3))

ABBA_LocalNP_Hetero <- subset(FoliarCNP_bF, LocalComType == "NP: Heterospecific")
str(ABBA_LocalNP_Hetero)
ABBA_Hetero_NP_C <- paste0(round(mean(ABBA_LocalNP_Hetero$C), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalNP_Hetero$C), digits=3))
ABBA_Hetero_NP_N <- paste0(round(mean(ABBA_LocalNP_Hetero$N), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalNP_Hetero$N), digits=3))
ABBA_Hetero_NP_P <- paste0(round(mean(ABBA_LocalNP_Hetero$P), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalNP_Hetero$P), digits=3))

ABBA_LocalCF_Con <- subset(FoliarCNP_bF, LocalComType == "CF: Conspecific")
str(ABBA_LocalCF_Con)
ABBA_Con_CF_C <- paste0(round(mean(ABBA_LocalCF_Con$C), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalCF_Con$C), digits=3))
ABBA_Con_CF_N <- paste0(round(mean(ABBA_LocalCF_Con$N), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalCF_Con$N), digits=3))
ABBA_Con_CF_P <- paste0(round(mean(ABBA_LocalCF_Con$P), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalCF_Con$P), digits=3))

ABBA_LocalCF_Hetero <- subset(FoliarCNP_bF, LocalComType == "CF: Heterospecific")
str(ABBA_LocalCF_Hetero)
ABBA_Hetero_CF_C <- paste0(round(mean(ABBA_LocalCF_Hetero$C), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalCF_Hetero$C), digits=3))
ABBA_Hetero_CF_N <- paste0(round(mean(ABBA_LocalCF_Hetero$N), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalCF_Hetero$N), digits=3))
ABBA_Hetero_CF_P <- paste0(round(mean(ABBA_LocalCF_Hetero$P), digits=3), " ", "\\u00b1 ",round(se(ABBA_LocalCF_Hetero$P), digits=3))

# Get differences between means
ABBA_Con_C_Diff <- paste0(round(mean(ABBA_LocalNP_Con$C) - mean(ABBA_LocalCF_Con$C), digits=3))
ABBA_Con_N_Diff <- paste0(round(mean(ABBA_LocalNP_Con$N) - mean(ABBA_LocalCF_Con$N), digits=3))
ABBA_Con_P_Diff <- paste0(round(mean(ABBA_LocalNP_Con$P) - mean(ABBA_LocalCF_Con$P), digits=3))

ABBA_Hetero_C_Diff <- paste0(round(mean(ABBA_LocalNP_Hetero$C) - mean(ABBA_LocalCF_Hetero$C), digits=3))
ABBA_Hetero_N_Diff <- paste0(round(mean(ABBA_LocalNP_Hetero$N) - mean(ABBA_LocalCF_Hetero$N), digits=3))
ABBA_Hetero_P_Diff <- paste0(round(mean(ABBA_LocalNP_Hetero$P) - mean(ABBA_LocalCF_Hetero$P), digits=3))

# Balsam fir
ABBA_VolAvg_NP <- data.frame(cbind(ABBA_Con_NP_C, ABBA_Con_NP_N, ABBA_Con_NP_P), cbind(ABBA_Hetero_NP_C, ABBA_Hetero_NP_N, ABBA_Hetero_NP_P))
colnames(ABBA_VolAvg_NP) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
ABBA_VolAvg_NP$SiteName <- "Local NP"

ABBA_VolAvg_CF <- data.frame(cbind(ABBA_Con_CF_C, ABBA_Con_CF_N, ABBA_Con_CF_P), cbind(ABBA_Hetero_CF_C, ABBA_Hetero_CF_N, ABBA_Hetero_CF_P))
colnames(ABBA_VolAvg_CF) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
ABBA_VolAvg_CF$SiteName <- "Local CF"

ABBA_VolAvg_Diff <- data.frame(cbind(ABBA_Con_C_Diff, ABBA_Con_N_Diff, ABBA_Con_P_Diff), cbind(ABBA_Hetero_C_Diff, ABBA_Hetero_N_Diff, ABBA_Hetero_P_Diff))
colnames(ABBA_VolAvg_Diff) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
ABBA_VolAvg_Diff$SiteName <- "Difference"

ABBA_VolAvg <- data.frame(rbind(ABBA_VolAvg_NP, ABBA_VolAvg_CF, ABBA_VolAvg_Diff))


# White birch
BEPA_RegionalCon <- subset(FoliarCNP_wB, RegComType == "Conspecific")
str(BEPA_RegionalCon)
BEPA_RegionalHetero <- subset(FoliarCNP_wB, RegComType == "Heterospecific")
str(BEPA_RegionalHetero)

BEPA_LocalNP_Con <- subset(FoliarCNP_wB, LocalComType == "NP: Conspecific")
str(BEPA_LocalNP_Con)
BEPA_Con_NP_C <- paste0(round(mean(BEPA_LocalNP_Con$C), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalNP_Con$C), digits=3))
BEPA_Con_NP_N <- paste0(round(mean(BEPA_LocalNP_Con$N), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalNP_Con$N), digits=3))
BEPA_Con_NP_P <- paste0(round(mean(BEPA_LocalNP_Con$P), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalNP_Con$P), digits=3))

BEPA_LocalNP_Hetero <- subset(FoliarCNP_wB, LocalComType == "NP: Heterospecific")
str(BEPA_LocalNP_Hetero)
BEPA_Hetero_NP_C <- paste0(round(mean(BEPA_LocalNP_Hetero$C), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalNP_Hetero$C), digits=3))
BEPA_Hetero_NP_N <- paste0(round(mean(BEPA_LocalNP_Hetero$N), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalNP_Hetero$N), digits=3))
BEPA_Hetero_NP_P <- paste0(round(mean(BEPA_LocalNP_Hetero$P), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalNP_Hetero$P), digits=3))

BEPA_LocalCF_Con <- subset(FoliarCNP_wB, LocalComType == "CF: Conspecific")
str(BEPA_LocalCF_Con)
BEPA_Con_CF_C <- paste0(round(mean(BEPA_LocalCF_Con$C), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalCF_Con$C), digits=3))
BEPA_Con_CF_N <- paste0(round(mean(BEPA_LocalCF_Con$N), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalCF_Con$N), digits=3))
BEPA_Con_CF_P <- paste0(round(mean(BEPA_LocalCF_Con$P), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalCF_Con$P), digits=3))

BEPA_LocalCF_Hetero <- subset(FoliarCNP_wB, LocalComType == "CF: Heterospecific")
str(BEPA_LocalCF_Hetero)
BEPA_Hetero_CF_C <- paste0(round(mean(BEPA_LocalCF_Hetero$C), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalCF_Hetero$C), digits=3))
BEPA_Hetero_CF_N <- paste0(round(mean(BEPA_LocalCF_Hetero$N), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalCF_Hetero$N), digits=3))
BEPA_Hetero_CF_P <- paste0(round(mean(BEPA_LocalCF_Hetero$P), digits=3), " ", "\\u00b1 ",round(se(BEPA_LocalCF_Hetero$P), digits=3))

# Get differences between means
BEPA_Con_C_Diff <- paste0(round(mean(BEPA_LocalNP_Con$C) - mean(BEPA_LocalCF_Con$C), digits=3))
BEPA_Con_N_Diff <- paste0(round(mean(BEPA_LocalNP_Con$N) - mean(BEPA_LocalCF_Con$N), digits=3))
BEPA_Con_P_Diff <- paste0(round(mean(BEPA_LocalNP_Con$P) - mean(BEPA_LocalCF_Con$P), digits=3))

BEPA_Hetero_C_Diff <- paste0(round(mean(BEPA_LocalNP_Hetero$C) - mean(BEPA_LocalCF_Hetero$C), digits=3))
BEPA_Hetero_N_Diff <- paste0(round(mean(BEPA_LocalNP_Hetero$N) - mean(BEPA_LocalCF_Hetero$N), digits=3))
BEPA_Hetero_P_Diff <- paste0(round(mean(BEPA_LocalNP_Hetero$P) - mean(BEPA_LocalCF_Hetero$P), digits=3))

# White birch
BEPA_VolAvg_NP <- data.frame(cbind(BEPA_Con_NP_C, BEPA_Con_NP_N, BEPA_Con_NP_P), cbind(BEPA_Hetero_NP_C, BEPA_Hetero_NP_N, BEPA_Hetero_NP_P))
colnames(BEPA_VolAvg_NP) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
BEPA_VolAvg_NP$SiteName <- "Local NP"

BEPA_VolAvg_CF <- data.frame(cbind(BEPA_Con_CF_C, BEPA_Con_CF_N, BEPA_Con_CF_P), cbind(BEPA_Hetero_CF_C, BEPA_Hetero_CF_N, BEPA_Hetero_CF_P))
colnames(BEPA_VolAvg_CF) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
BEPA_VolAvg_CF$SiteName <- "Local CF"

BEPA_VolAvg_Diff <- data.frame(cbind(BEPA_Con_C_Diff, BEPA_Con_N_Diff, BEPA_Con_P_Diff), cbind(BEPA_Hetero_C_Diff, BEPA_Hetero_N_Diff, BEPA_Hetero_P_Diff))
colnames(BEPA_VolAvg_Diff) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
BEPA_VolAvg_Diff$SiteName <- "Difference"

BEPA_VolAvg <- data.frame(rbind(BEPA_VolAvg_NP, BEPA_VolAvg_CF, BEPA_VolAvg_Diff))

# Make table
Species <- data.frame(c("bF", "bF", "bF", "wB", "wB", "wB"))
colnames(Species) <- c("Species")

CNP_Diff <- data.frame(cbind(Species, rbind(ABBA_VolAvg, BEPA_VolAvg)))

write.csv(CNP_Diff, "Tables/CNP_Diff.csv")

# Get C, N, P mean and difference at a species level between the two ecoregions
# Balsam fir
bF_NP_meanC <- paste0(round(mean(bF_NP$C), digits=3), " ", "\\u00b1 ",round(se(bF_NP$C), digits=3))
bF_NP_meanN <- paste0(round(mean(bF_NP$N), digits=3), " ", "\\u00b1 ",round(se(bF_NP$N), digits=3))
bF_NP_meanP <- paste0(round(mean(bF_NP$P), digits=3), " ", "\\u00b1 ",round(se(bF_NP$P), digits=3))

bF_CF_meanC <- paste0(round(mean(bF_CF$C), digits=3), " ", "\\u00b1 ",round(se(bF_CF$C), digits=3))
bF_CF_meanN <- paste0(round(mean(bF_CF$N), digits=3), " ", "\\u00b1 ",round(se(bF_CF$N), digits=3))
bF_CF_meanP <- paste0(round(mean(bF_CF$P), digits=3), " ", "\\u00b1 ",round(se(bF_CF$P), digits=3))

bF_Ecoreg_meanC_Diff <- paste0(round(mean(bF_NP$C) - mean(bF_CF$C), digits=3))
bF_Ecoreg_meanN_Diff <- paste0(round(mean(bF_NP$N) - mean(bF_CF$N), digits=3))
bF_Ecoreg_meanP_Diff <- paste0(round(mean(bF_NP$P) - mean(bF_CF$P), digits=3))

# White birch
wB_NP_meanC <- paste0(round(mean(wB_NP$C), digits=3), " ", "\\u00b1 ",round(se(wB_NP$C), digits=3))
wB_NP_meanN <- paste0(round(mean(wB_NP$N), digits=3), " ", "\\u00b1 ",round(se(wB_NP$N), digits=3))
wB_NP_meanP <- paste0(round(mean(wB_NP$P), digits=3), " ", "\\u00b1 ",round(se(wB_NP$P), digits=3))

wB_CF_meanC <- paste0(round(mean(wB_CF$C), digits=3), " ", "\\u00b1 ",round(se(wB_CF$C), digits=3))
wB_CF_meanN <- paste0(round(mean(wB_CF$N), digits=3), " ", "\\u00b1 ",round(se(wB_CF$N), digits=3))
wB_CF_meanP <- paste0(round(mean(wB_CF$P), digits=3), " ", "\\u00b1 ",round(se(wB_CF$P), digits=3))

wB_Ecoreg_meanC_Diff <- paste0(round(mean(wB_NP$C) - mean(wB_CF$C), digits=3))
wB_Ecoreg_meanN_Diff <- paste0(round(mean(wB_NP$N) - mean(wB_CF$N), digits=3))
wB_Ecoreg_meanP_Diff <- paste0(round(mean(wB_NP$P) - mean(wB_CF$P), digits=3))

# Create table of means 
# Balsam fir
bF_NP_means <- data.frame(cbind(bF_NP_meanC, bF_NP_meanN, bF_NP_meanP))
colnames(bF_NP_means) <- c("C", "N", "P")
bF_NP_means$Ecoregion <- "Northern Peninsula"

bF_CF_means <- data.frame(cbind(bF_CF_meanC, bF_CF_meanN, bF_CF_meanP))
colnames(bF_CF_means) <- c("C", "N", "P")
bF_CF_means$Ecoregion <- "Central Forest"

bF_meansDiff <- data.frame(cbind(bF_Ecoreg_meanC_Diff, bF_Ecoreg_meanN_Diff, bF_Ecoreg_meanP_Diff))
colnames(bF_meansDiff) <- c("C", "N", "P")
bF_meansDiff$Ecoregion <- "Difference"

bF_NP_CF_means <- data.frame(rbind(bF_NP_means, bF_CF_means, bF_meansDiff))
bF_NP_CF_means$Species <- "bF"

# White birch
wB_NP_means <- data.frame(cbind(wB_NP_meanC, wB_NP_meanN, wB_NP_meanP))
colnames(wB_NP_means) <- c("C", "N", "P")
wB_NP_means$Ecoregion <- "Northern Peninsula"

wB_CF_means <- data.frame(cbind(wB_CF_meanC, wB_CF_meanN, wB_CF_meanP))
colnames(wB_CF_means) <- c("C", "N", "P")
wB_CF_means$Ecoregion <- "Central Forest"

wB_meansDiff <- data.frame(cbind(wB_Ecoreg_meanC_Diff, wB_Ecoreg_meanN_Diff, wB_Ecoreg_meanP_Diff))
colnames(wB_meansDiff) <- c("C", "N", "P")
wB_meansDiff$Ecoregion <- "Difference"

wB_NP_CF_means <- data.frame(rbind(wB_NP_means, wB_CF_means, wB_meansDiff))
wB_NP_CF_means$Species <- "wB"

bF_wB_Ecoregion_meanDiff <- data.frame(cbind(bF_NP_CF_means, wB_NP_CF_means))

write.csv(bF_wB_Ecoregion_meanDiff, "Tables/Ecoregion_CNP_Diff.csv")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
###### Perform PCA #####
library(FactoMineR)
library(factoextra)
library(ggsci) # for the colour palette
library(scales)

FoliarCNP_bF$Region <- factor(FoliarCNP_bF$Region, levels = c("Northern Peninsula", "Central Forest"))
FoliarCNP_wB$Region <- factor(FoliarCNP_wB$Region, levels = c("Northern Peninsula", "Central Forest"))

bF_PCA <- PCA(FoliarCNP_bF[5:7], scale.unit = TRUE, ncp = 2, graph = TRUE)
wB_PCA <- PCA(FoliarCNP_wB[5:7], scale.unit = TRUE, ncp = 2, graph = TRUE)

# Create PCA biplot
# Set PCA colours
RegCol <- c("Conspecific" = "#3E3E23FF", "Heterospecific" = "#CC8214FF") # Regional extent
LocalCol <- c("NP: Conspecific" = "#5B8FA8FF", "NP: Heterospecific" = "#616530FF", "CF: Conspecific" = "#0F425CFF", "CF: Heterospecific" = "#9A5324FF") # Local extent
show_col(pal_uchicago("dark")(9))

# PCA 
# Balsam fir
bF_Ecoreg_PCA_plot <- fviz_pca_biplot(bF_PCA, col.ind = FoliarCNP_bF$Region, addEllipses = TRUE, ellipse.level=0.95, label = "var", col.var = "black", repel = TRUE, legend.title = "Ecoregion", title = "Between ecoregions\\n(a) Balsam fir") + theme(text = element_text(size = 11)) + scale_color_manual(values=c("#9A5324FF", "#767676FF")) + scale_fill_manual(values=c("#9A5324FF", "#767676FF"))

bF_PCA_Reg <- fviz_pca_biplot(bF_PCA, col.ind = FoliarCNP_bF$RegComType, palette = RegCol, addEllipses = TRUE, ellipse.level=0.95, label = "var", col.var = "black", repel = TRUE, legend.title = "Community Type", title = "Community types across ecoregions\\n(c) Balsam fir", legend = "none") + theme(text = element_text(size = 11))
bF_PCA_Local <- fviz_pca_biplot(bF_PCA, col.ind = FoliarCNP_bF$LocalComType, palette = LocalCol, addEllipses = TRUE, ellipse.level=0.95, label = "var", col.var = "black", repel = TRUE, legend.title = "Community Type", title = "Community types within and between ecoregions\\n(e) Balsam fir", legend = "none") + theme(text = element_text(size = 11)) 


# White birch 
wB_Ecoreg_PCA_plot <- fviz_pca_biplot(wB_PCA, col.ind = FoliarCNP_wB$Region, addEllipses = TRUE, ellipse.level=0.95, label = "var", col.var = "black", repel = TRUE, legend.title = "Ecoregion", title = "\\n(b) White birch") + theme(text = element_text(size = 11)) + scale_color_manual(values=c("#9A5324FF", "#767676FF")) + scale_fill_manual(values=c("#9A5324FF", "#767676FF"))

wB_PCA_Reg <- fviz_pca_biplot(wB_PCA, col.ind = FoliarCNP_wB$RegComType, palette = RegCol, addEllipses = TRUE, ellipse.level=0.95, label = "var", col.var = "black", repel = TRUE, legend.title = "Community Type", title = "\\n(d) White birch", legend = "none") + theme(text = element_text(size = 11))
wB_PCA_Local <- fviz_pca_biplot(wB_PCA, col.ind = FoliarCNP_wB$LocalComType, palette = LocalCol, addEllipses = TRUE, ellipse.level=0.95, label = "var", col.var = "black", repel = TRUE, legend.title = "Community Type", title = "\\n(f) White birch", legend = "none") + theme(text = element_text(size = 11))


library(gridExtra)
library(ggpubr)
Ecoreg_ggarrange <- ggarrange(bF_Ecoreg_PCA_plot, wB_Ecoreg_PCA_plot, ncol=2, nrow=1, common.legend = T, legend = "bottom")
Reg_ggarrange <- ggarrange(bF_PCA_Reg, wB_PCA_Reg, ncol=2, nrow=1, common.legend = T, legend = "bottom")
Local_ggarrange <- ggarrange(bF_PCA_Local, wB_PCA_Local, ncol=2, nrow=1, common.legend = T, legend = "bottom")

tiff("Figures/PCA.tiff", width = 8, height = 10, units = "in", res = 600) 
ggarrange(Ecoreg_ggarrange, Reg_ggarrange, Local_ggarrange, ncol=1, nrow=3)
dev.off()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
###### Scatter plot ######
# Create ggplots of C, N, and P in 2 dimensional space

EcoregCol <- c("Northern Peninsula" = "#9A5324FF", "Central Forest" = "#767676FF") # Ecoregions
#test
library(ggpmisc)
my.formula <- y ~ x

FoliarCNP_bF$Region <- factor(FoliarCNP_bF$Region, levels = c("Northern Peninsula", "Central Forest"))
FoliarCNP_wB$Region <- factor(FoliarCNP_wB$Region, levels = c("Northern Peninsula", "Central Forest"))


#bF
bF_Ecoreg_CN <- ggplot(FoliarCNP_bF, aes(x=N, y=C, col=Region)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "N") + theme_linedraw() + ggtitle("Balsam fir\\n(a)") + theme(text = element_text(size=11)) + scale_color_manual(name = "Ecoregions", values=EcoregCol) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "right", label.y = "bottom", formula=formula, parse=T, size = 4) 

bF_Ecoreg_CP <- ggplot(FoliarCNP_bF, aes(x=P, y=C, col=Region)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "P") + theme_linedraw() + scale_color_manual(name = "Ecoregions", values=EcoregCol) + ggtitle("\\n(b)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "right", label.y = "top", formula=formula, parse=T, size = 4) 

bF_Ecoreg_NP <- ggplot(FoliarCNP_bF, aes(x=P, y=N, col=Region)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "N", x = "P") + theme_linedraw() + scale_color_manual(name = "Ecoregions", values=EcoregCol) + ggtitle("\\n(c)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", formula=formula, parse=T, size = 4) 

#wB
wB_Ecoreg_CN <- ggplot(FoliarCNP_wB, aes(x=N, y=C, col=Region)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "N") + theme_linedraw() + scale_color_manual(name = "Ecoregions", values=EcoregCol) +ggtitle("White birch\\n(d)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", formula=formula, parse=T, size = 4) 

wB_Ecoreg_CP <- ggplot(FoliarCNP_wB, aes(x=P, y=C, col=Region)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "P") + theme_linedraw() + scale_color_manual(name = "Ecoregions", values=EcoregCol) + ggtitle("\\n(e)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", formula=formula, parse=T, size = 4) 

wB_Ecoreg_NP <- ggplot(FoliarCNP_wB, aes(x=P, y=N, col=Region)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "N", x = "P") + theme_linedraw() + scale_color_manual(name = "Ecoregions", values=EcoregCol) + ggtitle("\\n(f)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", formula=formula, parse=T, size = 4) 


# Plot
library(cowplot)
Scatter_Ecoreg_legend <- get_legend(bF_Ecoreg_CN + guides(color=guide_legend(nrow=1)) + theme(legend.position = "bottom"))
Scatter_Ecoreg_title <- ggdraw() + draw_label("Between ecoregions", hjust=0, x=0.05) + theme(plot.margin = margin(4, 0, 7, 0))
Scatter_Ecoreg_plots <- plot_grid(bF_Ecoreg_CN, bF_Ecoreg_CP, bF_Ecoreg_NP, wB_Ecoreg_CN, wB_Ecoreg_CP, wB_Ecoreg_NP, ncol=3, nrow=2)

tiff("Figures/Scatter_BetweenEcoregions.tiff", width = 9, height = 8, units = "in", res = 600) 
plot_grid(Scatter_Ecoreg_title, Scatter_Ecoreg_plots, Scatter_Ecoreg_legend, ncol=1, rel_heights = c(0.09, 3, 0.15))
dev.off()



# Regional
# bF
bF_Reg_CN <- ggplot(FoliarCNP_bF, aes(x=N, y=C, col=RegComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "N") + theme_linedraw() + scale_color_manual(name = "Community Types", values=RegCol) + ggtitle("Balsam fir\\n(a)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "right", label.y = "bottom", formula=formula, parse=T, size = 4) 

bF_Reg_CP <- ggplot(FoliarCNP_bF, aes(x=P, y=C, col=RegComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=RegCol) + ggtitle("\\n(b)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "right", label.y = "bottom", formula=formula, parse=T, size = 4) 

bF_Reg_NP <- ggplot(FoliarCNP_bF, aes(x=P, y=N, col=RegComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "N", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=RegCol) + ggtitle("\\n(c)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", label.y = "top", formula=formula, parse=T, size = 4) 

# wB
wB_Reg_CN <- ggplot(FoliarCNP_wB, aes(x=N, y=C, col=RegComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "N") + theme_linedraw() + scale_color_manual(name = "Community Types", values=RegCol) + ggtitle("White birch\\n(d)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", label.y="top", formula=formula, parse=T, size = 4) 

wB_Reg_CP <- ggplot(FoliarCNP_wB, aes(x=P, y=C, col=RegComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=RegCol) + ggtitle("\\n(e)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", label.y="top", formula=formula, parse=T, size = 4) 

wB_Reg_NP <- ggplot(FoliarCNP_wB, aes(x=P, y=N, col=RegComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "N", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=RegCol) + ggtitle("\\n(f)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", label.y="top", formula=formula, parse=T, size = 4) 


Scatter_Reg_legend <- get_legend(wB_Reg_NP + guides(color=guide_legend(nrow=1)) + theme(legend.position = "bottom"))
Scatter_Reg_title <- ggdraw() + draw_label("Community types across ecoregions", hjust=0, x=0.05) + theme(plot.margin = margin(4, 0, 7, 0))
Scatter_Reg_plots <- plot_grid(bF_Reg_CN, bF_Reg_CP, bF_Reg_NP, wB_Reg_CN, wB_Reg_CP, wB_Reg_NP, ncol=3, nrow=2)

tiff("Figures/Scatter_AcrossEcoregions.tiff", width = 9, height = 8, units = "in", res = 600) 
plot_grid(Scatter_Reg_title, Scatter_Reg_plots, Scatter_Reg_legend, ncol=1, rel_heights = c(0.09, 3, 0.15))
dev.off()


# Local
bF_Local_CN <- ggplot(FoliarCNP_bF, aes(x=N, y=C, col=LocalComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "N") + theme_linedraw() + scale_color_manual(name = "Community Types", values=LocalCol) + ggtitle("Balsam fir\\n(a)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "right",  label.y = "bottom", formula=formula, parse=T, size = 4) 

bF_Local_CP <- ggplot(FoliarCNP_bF, aes(x=P, y=C, col=LocalComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=LocalCol) + ggtitle("\\n(b)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "right", label.y = "bottom", formula=formula, parse=T, size = 4) 

bF_Local_NP <- ggplot(FoliarCNP_bF, aes(x=P, y=N, col=LocalComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "N", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=LocalCol) + ggtitle("\\n(c)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", label.y = "top", formula=formula, parse=T, size = 4) 

# wB
wB_Local_CN <- ggplot(FoliarCNP_wB, aes(x=N, y=C, col=LocalComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "N") + theme_linedraw() + scale_color_manual(name = "Community Types", values=LocalCol) + ggtitle("White birch\\n(d)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", label.y = "top", formula=formula, parse=T, size = 4) 

wB_Local_CP <- ggplot(FoliarCNP_wB, aes(x=P, y=C, col=LocalComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=LocalCol) + ggtitle("\\n(e)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", label.y = "top", formula=formula, parse=T, size = 4) 

wB_Local_NP <- ggplot(FoliarCNP_wB, aes(x=P, y=N, col=LocalComType)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "N", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=LocalCol) + ggtitle("\\n(f)") + theme(text = element_text(size=11)) + theme(legend.position = "none") + geom_smooth(method="lm", formula = my.formula) + stat_poly_eq(aes(label=paste(..rr.label.., ..p.value.label.., sep="*','~")), rr.digits=3, p.digits = 3, label.x = "left", label.y = "top", formula=formula, parse=T, size = 4) 


Scatter_Local_legend <- get_legend(wB_Local_NP + guides(color=guide_legend(nrow=1)) + theme(legend.position = "bottom"))
Scatter_Local_title <- ggdraw() + draw_label("Community types within and between ecoregions", hjust=0, x=0.05) + theme(plot.margin = margin(4, 0, 7, 0))
Scatter_Local_plots <- plot_grid(bF_Local_CN, bF_Local_CP, bF_Local_NP, wB_Local_CN, wB_Local_CP, wB_Local_NP, ncol=3, nrow=2)

tiff("Figures/Scatter_WithinBetweenEcoregions.tiff", width = 9, height = 8, units = "in", res = 600) 
plot_grid(Scatter_Local_title, Scatter_Local_plots, Scatter_Local_legend, ncol=1, rel_heights = c(0.09, 3, 0.15))
dev.off()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
###### Hypervolume calculaution ######
# Hypervolume
library("hypervolume")
# Ecoregions
# Balsam fir
HV_bF_Ecoreg_NP <- hypervolume_gaussian(FoliarCNP_bF[which(FoliarCNP_bF$Region == "Northern Peninsula"),names(FoliarCNP_bF) %in% c("C","N","P")])
HV_bF_Ecoreg_CF <- hypervolume_gaussian(FoliarCNP_bF[which(FoliarCNP_bF$Region == "Central Forest"),names(FoliarCNP_bF) %in% c("C","N","P")])
# White birch
HV_wB_Ecoreg_NP <- hypervolume_gaussian(FoliarCNP_wB[which(FoliarCNP_wB$Region == "Northern Peninsula"),names(FoliarCNP_wB) %in% c("C","N","P")])
HV_wB_Ecoreg_CF <- hypervolume_gaussian(FoliarCNP_wB[which(FoliarCNP_wB$Region == "Central Forest"),names(FoliarCNP_wB) %in% c("C","N","P")])

# Regional balsam fir
HV_bF_Regional_Hetero <- hypervolume_gaussian(FoliarCNP_bF[which(FoliarCNP_bF$RegComType == "Heterospecific"),names(FoliarCNP_bF) %in% c("C","N","P")])
HV_bF_Regional_Cons <- hypervolume_gaussian(FoliarCNP_bF[which(FoliarCNP_bF$RegComType == "Conspecific"),names(FoliarCNP_bF) %in% c("C","N","P")])
# Local balsam fir
HV_bF_Local_NP_Hetero <- hypervolume_gaussian(FoliarCNP_bF[which(FoliarCNP_bF$LocalComType == "NP: Heterospecific"),names(FoliarCNP_bF) %in% c("C","N","P")])
HV_bF_Local_NP_Cons <- hypervolume_gaussian(FoliarCNP_bF[which(FoliarCNP_bF$LocalComType == "NP: Conspecific"),names(FoliarCNP_bF) %in% c("C","N","P")])
HV_bF_Local_CF_Hetero <- hypervolume_gaussian(FoliarCNP_bF[which(FoliarCNP_bF$LocalComType == "CF: Heterospecific"),names(FoliarCNP_bF) %in% c("C","N","P")])
HV_bF_Local_CF_Cons <- hypervolume_gaussian(FoliarCNP_bF[which(FoliarCNP_bF$LocalComType == "CF: Conspecific"),names(FoliarCNP_bF) %in% c("C","N","P")])

# Regional white birch
HV_wB_Regional_Hetero <- hypervolume_gaussian(FoliarCNP_wB[which(FoliarCNP_wB$RegComType == "Heterospecific"),names(FoliarCNP_wB) %in% c("C","N","P")])
HV_wB_Regional_Cons <- hypervolume_gaussian(FoliarCNP_wB[which(FoliarCNP_wB$RegComType == "Conspecific"),names(FoliarCNP_wB) %in% c("C","N","P")])
# Local white birch
HV_wB_Local_NP_Hetero <- hypervolume_gaussian(FoliarCNP_wB[which(FoliarCNP_wB$LocalComType == "NP: Heterospecific"),names(FoliarCNP_wB) %in% c("C","N","P")])
HV_wB_Local_NP_Cons <- hypervolume_gaussian(FoliarCNP_wB[which(FoliarCNP_wB$LocalComType == "NP: Conspecific"),names(FoliarCNP_wB) %in% c("C","N","P")])
HV_wB_Local_CF_Hetero <- hypervolume_gaussian(FoliarCNP_wB[which(FoliarCNP_wB$LocalComType == "CF: Heterospecific"),names(FoliarCNP_wB) %in% c("C","N","P")])
HV_wB_Local_CF_Cons <- hypervolume_gaussian(FoliarCNP_wB[which(FoliarCNP_wB$LocalComType == "CF: Conspecific"),names(FoliarCNP_wB) %in% c("C","N","P")])

# Group hypervolumes for between ecoregion level
HVS_bF_Ecoreg <- hypervolume_set(HV_bF_Ecoreg_NP, HV_bF_Ecoreg_CF, check.memory=FALSE)
HVS_wB_Ecoreg <- hypervolume_set(HV_wB_Ecoreg_NP, HV_wB_Ecoreg_CF, check.memory=FALSE)

# Group hypervolumes at the regional level
HVS_bF_Regional <- hypervolume_set(HV_bF_Regional_Hetero, HV_bF_Regional_Cons, check.memory=FALSE)
HVS_wB_Regional <- hypervolume_set(HV_wB_Regional_Hetero, HV_wB_Regional_Cons, check.memory=FALSE)

# Group hypervolumnes at the local level: balsam fir
HVS_bF_Local_NP <- hypervolume_set(HV_bF_Local_NP_Hetero, HV_bF_Local_NP_Cons, check.memory=FALSE)
HVS_bF_Local_CF <- hypervolume_set(HV_bF_Local_CF_Hetero, HV_bF_Local_CF_Cons, check.memory=FALSE)
HVS_bF_Local_Hetero <- hypervolume_set(HV_bF_Local_NP_Hetero, HV_bF_Local_CF_Hetero, check.memory=FALSE)
HVS_bF_Local_Cons <- hypervolume_set(HV_bF_Local_NP_Cons, HV_bF_Local_CF_Cons, check.memory=FALSE)

# Group hypervolumnes at the local level: white  birch
HVS_wB_Local_NP <- hypervolume_set(HV_wB_Local_NP_Hetero, HV_wB_Local_NP_Cons, check.memory=FALSE)
HVS_wB_Local_CF <- hypervolume_set(HV_wB_Local_CF_Hetero, HV_wB_Local_CF_Cons, check.memory=FALSE)
HVS_wB_Local_Hetero <- hypervolume_set(HV_wB_Local_NP_Hetero, HV_wB_Local_CF_Hetero, check.memory=FALSE)
HVS_wB_Local_Cons <- hypervolume_set(HV_wB_Local_NP_Cons, HV_wB_Local_CF_Cons, check.memory=FALSE)

# Comparisons
# Between ecoregions
HVS_bF_Ecoreg_Stats <- data.frame(hypervolume_overlap_statistics(HVS_bF_Ecoreg))
colnames(HVS_bF_Ecoreg_Stats) <- c("Index")
HVS_bF_Ecoreg_Stats$Comparison <- "bF: between ecoregions"

HVS_wB_Ecoreg_Stats <- data.frame(hypervolume_overlap_statistics(HVS_wB_Ecoreg))
colnames(HVS_wB_Ecoreg_Stats) <- c("Index")
HVS_wB_Ecoreg_Stats$Comparison <- "wB: between ecoregions"

# Regional
HVS_bF_Regional_Stats <- data.frame(hypervolume_overlap_statistics(HVS_bF_Regional))
colnames(HVS_bF_Regional_Stats) <- c("Index")
HVS_bF_Regional_Stats$Comparison <- "bF: regional con vs hetero"

HVS_wB_Regional_Stats <- data.frame(hypervolume_overlap_statistics(HVS_wB_Regional))
colnames(HVS_wB_Regional_Stats) <- c("Index")
HVS_wB_Regional_Stats$Comparison <- "wB: regional con vs hetero"

# Local
HVS_bF_Local_NP_Stats <- data.frame(hypervolume_overlap_statistics(HVS_bF_Local_NP))
colnames(HVS_bF_Local_NP_Stats) <- c("Index")
HVS_bF_Local_NP_Stats$Comparison <- "bF: local NP con vs hetero"

HVS_bF_Local_CF_Stats <- data.frame(hypervolume_overlap_statistics(HVS_bF_Local_CF))
colnames(HVS_bF_Local_CF_Stats) <- c("Index")
HVS_bF_Local_CF_Stats$Comparison <- "bF: local CF con vs hetero"

HVS_bF_Local_Hetero_Stats <- data.frame(hypervolume_overlap_statistics(HVS_bF_Local_Hetero))
colnames(HVS_bF_Local_Hetero_Stats) <- c("Index")
HVS_bF_Local_Hetero_Stats$Comparison <- "bF: local hetero"

HVS_bF_Local_Cons_Stats <- data.frame(hypervolume_overlap_statistics(HVS_bF_Local_Cons))
colnames(HVS_bF_Local_Cons_Stats) <- c("Index")
HVS_bF_Local_Cons_Stats$Comparison <- "bF: local con"

HVS_wB_Local_NP_Stats <- data.frame(hypervolume_overlap_statistics(HVS_wB_Local_NP))
colnames(HVS_wB_Local_NP_Stats) <- c("Index")
HVS_wB_Local_NP_Stats$Comparison <- "wB: local NP con vs hetero"

HVS_wB_Local_CF_Stats <- data.frame(hypervolume_overlap_statistics(HVS_wB_Local_CF))
colnames(HVS_wB_Local_CF_Stats) <- c("Index")
HVS_wB_Local_CF_Stats$Comparison <- "wB: local CF con vs hetero"

HVS_wB_Local_Hetero_Stats <- data.frame(hypervolume_overlap_statistics(HVS_wB_Local_Hetero))
colnames(HVS_wB_Local_Hetero_Stats) <- c("Index")
HVS_wB_Local_Hetero_Stats$Comparison <- "wB: local hetero"

HVS_wB_Local_Cons_Stats <- data.frame(hypervolume_overlap_statistics(HVS_wB_Local_Cons))
colnames(HVS_wB_Local_Cons_Stats) <- c("Index")
HVS_wB_Local_Cons_Stats$Comparison <- "wB: local con"

library("data.table")
# Create single table of hypervolumn results
HV_StatList <- list(HVS_bF_Ecoreg_Stats, HVS_wB_Ecoreg_Stats, HVS_bF_Regional_Stats, HVS_wB_Regional_Stats, HVS_bF_Local_NP_Stats, HVS_bF_Local_CF_Stats, HVS_bF_Local_Hetero_Stats, HVS_bF_Local_Cons_Stats, HVS_wB_Local_NP_Stats, HVS_wB_Local_CF_Stats, HVS_wB_Local_Hetero_Stats, HVS_wB_Local_Cons_Stats) # List of results in data frame
HV_StatsData <- rbindlist(lapply(HV_StatList, setDT, keep.rownames = TRUE)) # Bind lists together keeping row names
names(HV_StatsData)[names(HV_StatsData) == "rn"] <- "Statistic"

# Transpose for statistics as column names and comparison groups as row names
library("tidyr")
HV_StatsData_Spread <- spread(HV_StatsData, Statistic, Index)
HV_StatsData_Sub <- subset(HV_StatsData_Spread, select=c("Comparison", "jaccard", "sorensen", "frac_unique_1", "frac_unique_2"))
# Create vector with desired row order
HV_Order <- c("bF: between ecoregions", "wB: between ecoregions", "bF: regional con vs hetero", "bF: local NP con vs hetero", "bF: local CF con vs hetero", "bF: local con", "bF: local hetero", "wB: regional con vs hetero", "wB: local NP con vs hetero", "wB: local CF con vs hetero", "wB: local con", "wB: local hetero")
# re-order by HV_Order
HV_Stats_Data_Final <- HV_StatsData_Sub %>% mutate(Comparison =  factor(Comparison, levels = HV_Order)) %>% dplyr::arrange(Comparison)
# Export results
write.csv(HV_Stats_Data_Final, "Tables/Hypervolume_Stats.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
###### Multivariate normality test #####
library(mvShapiroTest)

# Ecoregion
bF_Ecoreg_NP <- mvShapiro.Test(data.matrix(FoliarCNP_bF[which(FoliarCNP_bF$Region == "Northern Peninsula"),names(FoliarCNP_bF) %in% c("C","N","P")]))
bF_Ecoreg_NP_Tab <- data.frame(c(bF_Ecoreg_NP$statistic, bF_Ecoreg_NP$p.value))
bF_Ecoreg_NP_Tab$Source <- "bF: ecoregion Northern Peninsula"
row.names(bF_Ecoreg_NP_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Ecoreg_NP_Tab$Metric <- row.names(bF_Ecoreg_NP_Tab) 
colnames(bF_Ecoreg_NP_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Ecoreg_NP_Tab) <- c()

bF_Ecoreg_CF <- mvShapiro.Test(data.matrix(FoliarCNP_bF[which(FoliarCNP_bF$Region == "Central Forest"),names(FoliarCNP_bF) %in% c("C","N","P")]))
bF_Ecoreg_CF_Tab <- data.frame(c(bF_Ecoreg_CF$statistic, bF_Ecoreg_CF$p.value))
bF_Ecoreg_CF_Tab$Source <- "bF: ecoregion Central Forest"
row.names(bF_Ecoreg_CF_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Ecoreg_CF_Tab$Metric <- row.names(bF_Ecoreg_CF_Tab) 
colnames(bF_Ecoreg_CF_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Ecoreg_CF_Tab) <- c()

wB_Ecoreg_NP <- mvShapiro.Test(data.matrix(FoliarCNP_wB[which(FoliarCNP_wB$Region == "Northern Peninsula"),names(FoliarCNP_wB) %in% c("C","N","P")]))
wB_Ecoreg_NP_Tab <- data.frame(c(wB_Ecoreg_NP$statistic, wB_Ecoreg_NP$p.value))
wB_Ecoreg_NP_Tab$Source <- "wB: ecoregion Northern Peninsula"
row.names(wB_Ecoreg_NP_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Ecoreg_NP_Tab$Metric <- row.names(wB_Ecoreg_NP_Tab) 
colnames(wB_Ecoreg_NP_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Ecoreg_NP_Tab) <- c()

wB_Ecoreg_CF <- mvShapiro.Test(data.matrix(FoliarCNP_wB[which(FoliarCNP_wB$Region == "Central Forest"),names(FoliarCNP_wB) %in% c("C","N","P")]))
wB_Ecoreg_CF_Tab <- data.frame(c(wB_Ecoreg_CF$statistic, wB_Ecoreg_CF$p.value))
wB_Ecoreg_CF_Tab$Source <- "wB: ecoregion Central Forest"
row.names(wB_Ecoreg_CF_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Ecoreg_CF_Tab$Metric <- row.names(wB_Ecoreg_CF_Tab) 
colnames(wB_Ecoreg_CF_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Ecoreg_CF_Tab) <- c()

# Regional
bF_Regional_Cons <- mvShapiro.Test(data.matrix(FoliarCNP_bF[which(FoliarCNP_bF$RegComType == "Conspecific"),names(FoliarCNP_bF) %in% c("C","N","P")]))
bF_Regional_Cons_Tab <- data.frame(c(bF_Regional_Cons$statistic, bF_Regional_Cons$p.value))
bF_Regional_Cons_Tab$Source <- "bF: regional con"
row.names(bF_Regional_Cons_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Regional_Cons_Tab$Metric <- row.names(bF_Regional_Cons_Tab) 
colnames(bF_Regional_Cons_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Regional_Cons_Tab) <- c()

bF_Regional_Hetero <- mvShapiro.Test(data.matrix(FoliarCNP_bF[which(FoliarCNP_bF$RegComType == "Heterospecific"),names(FoliarCNP_bF) %in% c("C","N","P")]))
bF_Regional_Hetero_Tab <- data.frame(c(bF_Regional_Hetero$statistic, bF_Regional_Hetero$p.value)) # Change values
bF_Regional_Hetero_Tab$Source <- "bF: regional hetero" # Change values
row.names(bF_Regional_Hetero_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Regional_Hetero_Tab$Metric <- row.names(bF_Regional_Hetero_Tab) 
colnames(bF_Regional_Hetero_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Regional_Hetero_Tab) <- c()

wB_Regional_Cons <- mvShapiro.Test(data.matrix(FoliarCNP_wB[which(FoliarCNP_wB$RegComType == "Conspecific"),names(FoliarCNP_wB) %in% c("C","N","P")]))
wB_Regional_Cons_Tab <- data.frame(c(wB_Regional_Cons$statistic, wB_Regional_Cons$p.value)) # Change values
wB_Regional_Cons_Tab$Source <- "wB: regional con" # Change values
row.names(wB_Regional_Cons_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Regional_Cons_Tab$Metric <- row.names(wB_Regional_Cons_Tab) 
colnames(wB_Regional_Cons_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Regional_Cons_Tab) <- c()

wB_Regional_Hetero <- mvShapiro.Test(data.matrix(FoliarCNP_wB[which(FoliarCNP_wB$RegComType == "Heterospecific"),names(FoliarCNP_wB) %in% c("C","N","P")]))
wB_Regional_Hetero_Tab <- data.frame(c(wB_Regional_Hetero$statistic, wB_Regional_Hetero$p.value)) # Change values
wB_Regional_Hetero_Tab$Source <- "wB: regional hetero" # Change values
row.names(wB_Regional_Hetero_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Regional_Hetero_Tab$Metric <- row.names(wB_Regional_Hetero_Tab) 
colnames(wB_Regional_Hetero_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Regional_Hetero_Tab) <- c()

# Local balsam fir
bF_Local_NP_Cons <- mvShapiro.Test(data.matrix(FoliarCNP_bF[which(FoliarCNP_bF$LocalComType == "NP: Conspecific"),names(FoliarCNP_bF) %in% c("C","N","P")]))
bF_Local_NP_Cons_Tab <- data.frame(c(bF_Local_NP_Cons$statistic, bF_Local_NP_Cons$p.value))
bF_Local_NP_Cons_Tab$Source <- "bF: local NP con"
row.names(bF_Local_NP_Cons_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Local_NP_Cons_Tab$Metric <- row.names(bF_Local_NP_Cons_Tab) 
colnames(bF_Local_NP_Cons_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Local_NP_Cons_Tab) <- c()

bF_Local_NP_Hetero <- mvShapiro.Test(data.matrix(FoliarCNP_bF[which(FoliarCNP_bF$LocalComType == "NP: Heterospecific"),names(FoliarCNP_bF) %in% c("C","N","P")]))
bF_Local_NP_Hetero_Tab <- data.frame(c(bF_Local_NP_Hetero$statistic, bF_Local_NP_Hetero$p.value))
bF_Local_NP_Hetero_Tab$Source <- "bF: local NP hetero"
row.names(bF_Local_NP_Hetero_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Local_NP_Hetero_Tab$Metric <- row.names(bF_Local_NP_Hetero_Tab) 
colnames(bF_Local_NP_Hetero_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Local_NP_Hetero_Tab) <- c()

bF_Local_CF_Cons <- mvShapiro.Test(data.matrix(FoliarCNP_bF[which(FoliarCNP_bF$LocalComType == "CF: Conspecific"),names(FoliarCNP_bF) %in% c("C","N","P")]))
bF_Local_CF_Cons_Tab <- data.frame(c(bF_Local_CF_Cons$statistic, bF_Local_CF_Cons$p.value))
bF_Local_CF_Cons_Tab$Source <- "bF: local CF con"
row.names(bF_Local_CF_Cons_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Local_CF_Cons_Tab$Metric <- row.names(bF_Local_CF_Cons_Tab) 
colnames(bF_Local_CF_Cons_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Local_CF_Cons_Tab) <- c()

bF_Local_CF_Hetero <- mvShapiro.Test(data.matrix(FoliarCNP_bF[which(FoliarCNP_bF$LocalComType == "CF: Heterospecific"),names(FoliarCNP_bF) %in% c("C","N","P")]))
bF_Local_CF_Hetero_Tab <- data.frame(c(bF_Local_CF_Hetero$statistic, bF_Local_CF_Hetero$p.value))
bF_Local_CF_Hetero_Tab$Source <- "bF: local CF hetero"
row.names(bF_Local_CF_Hetero_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Local_CF_Hetero_Tab$Metric <- row.names(bF_Local_CF_Hetero_Tab) 
colnames(bF_Local_CF_Hetero_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Local_CF_Hetero_Tab) <- c()

# Local white birch
# Can't compute wB NP cons because sample size is too low
wB_Local_NP_Cons <- mvShapiro.Test(data.matrix(FoliarCNP_wB[which(FoliarCNP_wB$LocalComType == "NP: Conspecific"),names(FoliarCNP_wB) %in% c("C","N","P")]))
wB_Local_NP_Cons_Tab <- data.frame(c(wB_Local_NP_Cons$statistic, wB_Local_NP_Cons$p.value))
wB_Local_NP_Cons_Tab$Source <- "wB: local NP con"
row.names(wB_Local_NP_Cons_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Local_NP_Cons_Tab$Metric <- row.names(wB_Local_NP_Cons_Tab) 
colnames(wB_Local_NP_Cons_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Local_NP_Cons_Tab) <- c()

wB_Local_NP_Hetero <- mvShapiro.Test(data.matrix(FoliarCNP_wB[which(FoliarCNP_wB$LocalComType == "NP: Heterospecific"),names(FoliarCNP_wB) %in% c("C","N","P")]))
wB_Local_NP_Hetero_Tab <- data.frame(c(wB_Local_NP_Hetero$statistic, wB_Local_NP_Hetero$p.value))
wB_Local_NP_Hetero_Tab$Source <- "wB: local NP hetero"
row.names(wB_Local_NP_Hetero_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Local_NP_Hetero_Tab$Metric <- row.names(wB_Local_NP_Hetero_Tab) 
colnames(wB_Local_NP_Hetero_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Local_NP_Hetero_Tab) <- c()

wB_Local_CF_Cons <- mvShapiro.Test(data.matrix(FoliarCNP_wB[which(FoliarCNP_wB$LocalComType == "CF: Conspecific"),names(FoliarCNP_wB) %in% c("C","N","P")]))
wB_Local_CF_Cons_Tab <- data.frame(c(wB_Local_CF_Cons$statistic, wB_Local_CF_Cons$p.value))
wB_Local_CF_Cons_Tab$Source <- "wB: local CF con"
row.names(wB_Local_CF_Cons_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Local_CF_Cons_Tab$Metric <- row.names(wB_Local_CF_Cons_Tab) 
colnames(wB_Local_CF_Cons_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Local_CF_Cons_Tab) <- c()

wB_Local_CF_Hetero <- mvShapiro.Test(data.matrix(FoliarCNP_wB[which(FoliarCNP_wB$LocalComType == "CF: Heterospecific"),names(FoliarCNP_wB) %in% c("C","N","P")]))
wB_Local_CF_Hetero_Tab <- data.frame(c(wB_Local_CF_Hetero$statistic, wB_Local_CF_Hetero$p.value))
wB_Local_CF_Hetero_Tab$Source <- "wB: local CF hetero"
row.names(wB_Local_CF_Hetero_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Local_CF_Hetero_Tab$Metric <- row.names(wB_Local_CF_Hetero_Tab) 
colnames(wB_Local_CF_Hetero_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Local_CF_Hetero_Tab) <- c()

MVN_Data <- rbind(bF_Ecoreg_NP_Tab, bF_Ecoreg_CF_Tab, wB_Ecoreg_NP_Tab, wB_Ecoreg_CF_Tab, bF_Regional_Cons_Tab, bF_Regional_Hetero_Tab, wB_Regional_Cons_Tab, wB_Regional_Hetero_Tab, bF_Local_NP_Cons_Tab, bF_Local_NP_Hetero_Tab, bF_Local_CF_Cons_Tab, bF_Local_CF_Hetero_Tab, wB_Local_NP_Hetero_Tab, wB_Local_CF_Cons_Tab, wB_Local_CF_Hetero_Tab)

library("tidyr")
MVN_Data_Spread <- spread(MVN_Data, Metric, Values)
MVN_Data_SpreadSub <- subset(MVN_Data_Spread, select=c("Source", "Shapiro_Wilk", "p_value"))
# Create vector with desired row order
MVN_Order <- c("bF: ecoregion Northern Peninsula", "bF: ecoregion Central Forest", "wB: ecoregion Northern Peninsula", "wB: ecoregion Central Forest", "bF: regional con", "bF: regional hetero", "bF: local NP con", "bF: local NP hetero", "bF: local CF con", "bF: local CF hetero", "wB: regional con", "wB: regional hetero", "wB: local NP hetero", "wB: local CF con", "wB: local CF hetero")
# re-order by HV_Order
MVN_Data_Final <- MVN_Data_SpreadSub %>% mutate(Source =  factor(Source, levels = MVN_Order)) %>% dplyr::arrange(Source)

write.csv(MVN_Data_Final, "Tables/MVN_Stats.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
###### MHD and PERMANOVA #######
# Subset by local extent
FoliarCNP_bF_NorthPen <- subset(FoliarCNP_bF, FoliarCNP_bF$Region == "Northern Peninsula") # n=295
FoliarCNP_bF_CenForest <- subset(FoliarCNP_bF, FoliarCNP_bF$Region == "Central Forest") # n=95
FoliarCNP_bF_Hetero <- subset(FoliarCNP_bF, FoliarCNP_bF$RegComType == "Heterospecific") # n=201
FoliarCNP_bF_Cons <- subset(FoliarCNP_bF, FoliarCNP_bF$RegComType == "Conspecific") # n=189

FoliarCNP_wB_NorthPen <- subset(FoliarCNP_wB, FoliarCNP_wB$Region == "Northern Peninsula") # n=158
FoliarCNP_wB_CenForest <- subset(FoliarCNP_wB, FoliarCNP_wB$Region == "Central Forest") # n=71
FoliarCNP_wB_Hetero <- subset(FoliarCNP_wB, FoliarCNP_wB$RegComType == "Heterospecific") #n=201
FoliarCNP_wB_Cons <- subset(FoliarCNP_wB, FoliarCNP_wB$RegComType == "Conspecific") # n=28

library(vegan)
# Compare ecoregions
PM_bF_Ecoreg_Dispr <- betadisper(d=vegdist(FoliarCNP_bF[5:7]), group=FoliarCNP_bF$Region)
permutest(PM_bF_Ecoreg_Dispr)
PM_bF_Ecoreg_Stats <- data.frame(adonis2(FoliarCNP_bF[5:7] ~ Region, data = FoliarCNP_bF))
PM_bF_Ecoreg_Stats$Source <- row.names(PM_bF_Ecoreg_Stats) 
PM_bF_Ecoreg_Stats$Source[PM_bF_Ecoreg_Stats$Source == "Region"] <- "bF: between ecoregions"

PM_wB_Ecoreg_Dispr <- betadisper(d=vegdist(FoliarCNP_wB[5:7]), group=FoliarCNP_wB$Region)
permutest(PM_wB_Ecoreg_Dispr)
PM_wB_Ecoreg_Stats <- data.frame(adonis2(FoliarCNP_wB[5:7] ~ Region, data = FoliarCNP_wB))
PM_wB_Ecoreg_Stats$Source <- row.names(PM_wB_Ecoreg_Stats) 
PM_wB_Ecoreg_Stats$Source[PM_wB_Ecoreg_Stats$Source == "Region"] <- "wB: between ecoregions"

# Compare regional conspecific and heterospecific groups
PM_bF_Regional_Dispr <- betadisper(d=vegdist(FoliarCNP_bF[5:7]), group=FoliarCNP_bF$RegComType)
permutest(PM_bF_Regional_Dispr)
PM_bF_Regional_Stats <- data.frame(adonis2(FoliarCNP_bF[5:7] ~ RegComType, data = FoliarCNP_bF))
PM_bF_Regional_Stats$Source <- row.names(PM_bF_Regional_Stats) 
PM_bF_Regional_Stats$Source[PM_bF_Regional_Stats$Source == "RegComType"] <- "bF: regional con vs hetero"

PM_wB_Regional_Dispr <- betadisper(d=vegdist(FoliarCNP_wB[5:7]), group=FoliarCNP_wB$RegComType)
permutest(PM_wB_Regional_Dispr)
PM_wB_Regional_Stats <- data.frame(adonis2(FoliarCNP_wB[5:7] ~ RegComType, data = FoliarCNP_wB))
PM_wB_Regional_Stats$Source <- row.names(PM_wB_Regional_Stats) 
PM_wB_Regional_Stats$Source[PM_wB_Regional_Stats$Source == "RegComType"] <- "wB: regional con vs hetero"


# Compare local conspecific and heterospecific: North Pen
# Balsam fir
PM_bF_Local_NP_Dispr <- betadisper(d=vegdist(FoliarCNP_bF_NorthPen[5:7]), group=FoliarCNP_bF_NorthPen$LocalComType)
permutest(PM_bF_Local_NP_Dispr)
PM_bF_Local_NP_Stats <- data.frame(adonis2(FoliarCNP_bF_NorthPen[5:7] ~ LocalComType, data = FoliarCNP_bF_NorthPen))
PM_bF_Local_NP_Stats$Source <- row.names(PM_bF_Local_NP_Stats) 
PM_bF_Local_NP_Stats$Source[PM_bF_Local_NP_Stats$Source == "LocalComType"] <- "bF: local NP con vs hetero"

PM_bF_Local_CF_Dispr <- betadisper(d=vegdist(FoliarCNP_bF_CenForest[5:7]), group=FoliarCNP_bF_CenForest$LocalComType)
permutest(PM_bF_Local_CF_Dispr)
PM_bF_Local_CF_Stats <- data.frame(adonis2(FoliarCNP_bF_CenForest[5:7] ~ LocalComType, data = FoliarCNP_bF_CenForest))
PM_bF_Local_CF_Stats$Source <- row.names(PM_bF_Local_CF_Stats) 
PM_bF_Local_CF_Stats$Source[PM_bF_Local_CF_Stats$Source == "LocalComType"] <- "bF: local CF con vs hetero"

PM_bF_Local_Hetero_Dispr <- betadisper(d=vegdist(FoliarCNP_bF_Hetero[5:7]), group=FoliarCNP_bF_Hetero$LocalComType)
permutest(PM_bF_Local_Hetero_Dispr)
PM_bF_Local_Hetero_Stats <- data.frame(adonis2(FoliarCNP_bF_Hetero[5:7] ~ LocalComType, data = FoliarCNP_bF_Hetero))
PM_bF_Local_Hetero_Stats$Source <- row.names(PM_bF_Local_Hetero_Stats) 
PM_bF_Local_Hetero_Stats$Source[PM_bF_Local_Hetero_Stats$Source == "LocalComType"] <- "bF: local hetero"

PM_bF_Local_Cons_Dispr <- betadisper(d=vegdist(FoliarCNP_bF_Cons[5:7]), group=FoliarCNP_bF_Cons$LocalComType)
permutest(PM_bF_Local_Cons_Dispr)
PM_bF_Local_Cons_Stats <- data.frame(adonis2(FoliarCNP_bF_Cons[5:7] ~ LocalComType, data = FoliarCNP_bF_Cons))
PM_bF_Local_Cons_Stats$Source <- row.names(PM_bF_Local_Cons_Stats) 
PM_bF_Local_Cons_Stats$Source[PM_bF_Local_Cons_Stats$Source == "LocalComType"] <- "bF: local con"

# White birch
PM_wB_Local_NP_Dispr <- betadisper(d=vegdist(FoliarCNP_wB_NorthPen[5:7]), group=FoliarCNP_wB_NorthPen$LocalComType)
permutest(PM_wB_Local_NP_Dispr)
PM_wB_Local_NP_Stats <- data.frame(adonis2(FoliarCNP_wB_NorthPen[5:7] ~ LocalComType, data = FoliarCNP_wB_NorthPen))
PM_wB_Local_NP_Stats$Source <- row.names(PM_wB_Local_NP_Stats) 
PM_wB_Local_NP_Stats$Source[PM_wB_Local_NP_Stats$Source == "LocalComType"] <- "wB: local NP con vs hetero"

PM_wB_Local_CF_Dispr <- betadisper(d=vegdist(FoliarCNP_wB_CenForest[5:7]), group=FoliarCNP_wB_CenForest$LocalComType)
permutest(PM_wB_Local_CF_Dispr)
PM_wB_Local_CF_Stats <- data.frame(adonis2(FoliarCNP_wB_CenForest[5:7] ~ LocalComType, data = FoliarCNP_wB_CenForest))
PM_wB_Local_CF_Stats$Source <- row.names(PM_wB_Local_CF_Stats) 
PM_wB_Local_CF_Stats$Source[PM_wB_Local_CF_Stats$Source == "LocalComType"] <- "wB: local CF con vs hetero"

PM_wB_Local_Hetero_Dispr <- betadisper(d=vegdist(FoliarCNP_wB_Hetero[5:7]), group=FoliarCNP_wB_Hetero$LocalComType)
permutest(PM_wB_Local_Hetero_Dispr)
PM_wB_Local_Hetero_Stats <- data.frame(adonis2(FoliarCNP_wB_Hetero[5:7] ~ LocalComType, data = FoliarCNP_wB_Hetero))
PM_wB_Local_Hetero_Stats$Source <- row.names(PM_wB_Local_Hetero_Stats) 
PM_wB_Local_Hetero_Stats$Source[PM_wB_Local_Hetero_Stats$Source == "LocalComType"] <- "wB: local hetero"

PM_wB_Local_Cons_Dispr <- betadisper(d=vegdist(FoliarCNP_wB_Cons[5:7]), group=FoliarCNP_wB_Cons$LocalComType)
permutest(PM_wB_Local_Cons_Dispr)
PM_wB_Local_Cons_Stats <- data.frame(adonis2(FoliarCNP_wB_Cons[5:7] ~ LocalComType, data = FoliarCNP_wB_Cons))
PM_wB_Local_Cons_Stats$Source <- row.names(PM_wB_Local_Cons_Stats) 
PM_wB_Local_Cons_Stats$Source[PM_wB_Local_Cons_Stats$Source == "LocalComType"] <- "wB: local con"

# List results by species
bF_List <- list(PM_bF_Ecoreg_Stats, PM_bF_Regional_Stats, PM_bF_Local_NP_Stats, PM_bF_Local_CF_Stats, PM_bF_Local_Cons_Stats, PM_bF_Local_Hetero_Stats)
wB_List <- list(PM_wB_Ecoreg_Stats, PM_wB_Regional_Stats, PM_wB_Local_NP_Stats, PM_wB_Local_CF_Stats, PM_wB_Local_Cons_Stats, PM_wB_Local_Hetero_Stats)
library(dplyr)
# Process table
bF_Row <- bind_rows(bF_List) # bind list together
bF_Rows <- bF_Row %>% dplyr::mutate_if(is.numeric, round, digits=4) # round numeric data to 3 decimal places
bF_RowsOrder <- subset(bF_Rows, select=c("Source", "Df", "SumOfSqs", "R2", "F", "Pr..F.")) # Subset data - change column order
colnames(bF_RowsOrder) <- c("Source", "Df", "SS", "R2", "F", "p-value") # Rename column names

wB_Row <- bind_rows(wB_List) # bind list together
wB_Rows <- wB_Row %>% mutate_if(is.numeric, round, digits=4) # round numeric data to 3 decimal places
wB_RowsOrder <- subset(wB_Rows, select=c("Source", "Df", "SumOfSqs", "R2", "F", "Pr..F.")) # Subset data - change column order
colnames(wB_RowsOrder) <- c("Source", "Df", "SS", "R2", "F", "p-value") # Rename column names

# Combine column
PM_Results <- cbind(bF_RowsOrder, wB_RowsOrder)
View(PM_Results)

# Export results
write.csv(PM_Results, "Tables/PERMANOVA_Stats.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
###### Niche volume metrics #######
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Gozalez metrics for calculating niche metrics
library(geometry)
library(rcdd)
# Subset data by species
# Subset by regional extent (i.e., species level)
FoliarCNP_bF <- subset(FoliarCNP, FoliarCNP$Species == "bF") # n=390
FoliarCNP_bF$RegComType <- as.factor(FoliarCNP_bF$RegComType)
FoliarCNP_wB <- subset(FoliarCNP, FoliarCNP$Species == "wB") # n=229
FoliarCNP_wB$RegComType <- as.factor(FoliarCNP_wB$RegComType)

# Calculate niche total available volumes
bF_NicheSpace <- as.matrix(FoliarCNP_bF[,c("C","N","P")])
bF_NicheSpace_TotVol <- CHVind(bF_NicheSpace)
bF_NicheSpace_TotVol <- as.numeric(bF_NicheSpace_TotVol$vol)

wB_NicheSpace <- as.matrix(FoliarCNP_wB[,c("C","N","P")])
wB_NicheSpace_TotVol <- CHVind(wB_NicheSpace)
wB_NicheSpace_TotVol <- as.numeric(wB_NicheSpace_TotVol$vol)

# Balsam fir subset for processing
bF_Spp <- subset(FoliarCNP_bF, select = c("C","N","P")) # Used to calculate total niche space
bF_Ecoreg <- subset(FoliarCNP_bF, select = c("Region","C","N","P" )) 
bF_RegCom <- subset(FoliarCNP_bF, select = c("RegComType","C","N","P" ))  # Regional extent
bF_LocalCom <- subset(FoliarCNP_bF, select = c("LocalComType","C","N","P" )) # Local extent

# White birch subset for processing
wB_Spp <- subset(FoliarCNP_wB, select = c("C","N","P"))
wB_Ecoreg <- subset(FoliarCNP_wB, select = c("Region","C","N","P" )) 
wB_RegCom <- subset(FoliarCNP_wB, select = c("RegComType","C","N","P" ))
wB_LocalCom <- subset(FoliarCNP_wB, select = c("LocalComType","C","N","P" ))



# Balsam fir ecoregion spatial extent
# This code generates a file that contains vol1, vol2, and volinter. Vol1= volume of the first population, vol2 = volume of the second population, and volinter = overlap between the two niches.
comparison=FoliarCNP_bF$Region[drop=T] # Change per variable
for(i in 1:length(levels(comparison[drop=T]))){
  loop=ifelse(comparison==levels(comparison[drop=T])[i],i,comparison)
}
cbindDB=cbind(bF_Ecoreg,loop) # Change per variable
vec_i=as.vector(1) # Change per level
vec_j=as.vector(2) # Change per level
comparison.names=paste(t(combn(levels(comparison[drop=T]),2))[,1],t(combn(levels(comparison[drop=T]),2))[,2],sep="-")
Ncomp<-matrix(NA,ncol=3,nrow=length(vec_i),dimnames=list(c(comparison.names),c("vol1","vol2","volinter")))
for(i in 1:length(vec_i)){
  for(j in 1:3) {
    Ncomp[i,j]<-round(CHVintersect(as.matrix(cbindDB[cbindDB$loop==vec_i[i],c("C","N","P")]),
                                   as.matrix(cbindDB[cbindDB$loop==vec_j[i],c("C","N","P")]))$vol[j],2)
  }
}
Ncomp=data.frame(ComType=rownames(Ncomp),Ncomp)
bF_Ecoreg_Vol=Ncomp%>%          # Change variable name                      
  mutate(voltot=vol1+vol2-volinter,                    
         vol.inter.rel=volinter/voltot*100,            
         minim=pmin(vol1,vol2),                          
         Turn=(2*minim-2*volinter)/(2*minim-volinter),
         Nest=(abs(vol1-vol2)/(vol1+vol2-volinter))*(volinter/(2*minim-volinter)))
mean(bF_Ecoreg_Vol$vol.inter.rel);se(bF_Ecoreg_Vol$vol.inter.rel) # mean and SE values of isotope volume # Change variable name   
mean(bF_Ecoreg_Vol$Turn);se(bF_Ecoreg_Vol$Turn) # Change variable name   
mean(bF_Ecoreg_Vol$Nest);se(bF_Ecoreg_Vol$Nest) # Change variable name   
bF_Ecoreg_Vol_Stats <- data.frame(bF_Ecoreg_Vol)
bF_Ecoreg_Vol_Stats$Species <- "bF"

# The following lines allow calculating the niche volumes of each population relative to the total niche volume
# occupied by all from all populations (i.e., Tot.vol)
bF_Ecoreg_VolRel=data.frame(Site.Species=c(levels(comparison[drop=TRUE])),Volume=c(bF_Ecoreg_Vol[1,2],bF_Ecoreg_Vol[1:length(levels(comparison[drop=TRUE]))-1,3])) # Change variable
bF_Ecoreg_VolRel=bF_Ecoreg_VolRel%>% # Change variable
  mutate(vol.sp.rel=(round((Volume/bF_NicheSpace_TotVol)*100,2))) # Change variable
mean(bF_Ecoreg_VolRel$vol.sp.rel);se(bF_Ecoreg_VolRel$vol.sp.rel) # Change variable
bF_EcoregVolRel <- data.frame(bF_Ecoreg_VolRel)
bF_EcoregVolRel$Species <- "bF"



# White birch ecoregion spatial extent
# This code generates a file that contains vol1, vol2, and volinter. Vol1= volume of the first population, vol2 = volume of the second population, and volinter = overlap between the two niches.
comparison=FoliarCNP_wB$Region[drop=T] # Change per variable
for(i in 1:length(levels(comparison[drop=T]))){
  loop=ifelse(comparison==levels(comparison[drop=T])[i],i,comparison)
}
cbindDB=cbind(wB_Ecoreg,loop) # Change per variable
vec_i=as.vector(1) # Change per level
vec_j=as.vector(2) # Change per level
comparison.names=paste(t(combn(levels(comparison[drop=T]),2))[,1],t(combn(levels(comparison[drop=T]),2))[,2],sep="-")
Ncomp<-matrix(NA,ncol=3,nrow=length(vec_i),dimnames=list(c(comparison.names),c("vol1","vol2","volinter")))
for(i in 1:length(vec_i)){
  for(j in 1:3) {
    Ncomp[i,j]<-round(CHVintersect(as.matrix(cbindDB[cbindDB$loop==vec_i[i],c("C","N","P")]),
                                   as.matrix(cbindDB[cbindDB$loop==vec_j[i],c("C","N","P")]))$vol[j],2)
  }
}
Ncomp=data.frame(ComType=rownames(Ncomp),Ncomp)
wB_Ecoreg_Vol=Ncomp%>%          # Change variable name                      
  mutate(voltot=vol1+vol2-volinter,                    
         vol.inter.rel=volinter/voltot*100,            
         minim=pmin(vol1,vol2),                          
         Turn=(2*minim-2*volinter)/(2*minim-volinter),
         Nest=(abs(vol1-vol2)/(vol1+vol2-volinter))*(volinter/(2*minim-volinter)))
mean(wB_Ecoreg_Vol$vol.inter.rel);se(wB_Ecoreg_Vol$vol.inter.rel) # mean and SE values of isotope volume # Change variable name   
mean(wB_Ecoreg_Vol$Turn);se(wB_Ecoreg_Vol$Turn) # Change variable name   
mean(wB_Ecoreg_Vol$Nest);se(wB_Ecoreg_Vol$Nest) # Change variable name   
wB_Ecoreg_Vol_Stats <- data.frame(wB_Ecoreg_Vol)
wB_Ecoreg_Vol_Stats$Species <- "wB"

# The following lines allow calculating the niche volumes of each population relative to the total niche volume
# occupied by all from all populations (i.e., Tot.vol)
wB_Ecoreg_VolRel=data.frame(Site.Species=c(levels(comparison[drop=TRUE])),Volume=c(wB_Ecoreg_Vol[1,2],wB_Ecoreg_Vol[1:length(levels(comparison[drop=TRUE]))-1,3])) # Change variable
wB_Ecoreg_VolRel=wB_Ecoreg_VolRel%>% # Change variable
  mutate(vol.sp.rel=(round((Volume/wB_NicheSpace_TotVol)*100,2))) # Change variable
mean(wB_Ecoreg_VolRel$vol.sp.rel);se(wB_Ecoreg_VolRel$vol.sp.rel) # Change variable
wB_EcoregVolRel <- data.frame(wB_Ecoreg_VolRel)
wB_EcoregVolRel$Species <- "wB"



# Balsam fir regional spatial extent
# This code generates a file that contains vol1, vol2, and volinter. Vol1= volume of the first population, vol2 = volume of the second population, and volinter = overlap between the two niches.
comparison=bF_RegCom$RegComType[drop=T] # Change per variable
for(i in 1:length(levels(comparison[drop=T]))){
  loop=ifelse(comparison==levels(comparison[drop=T])[i],i,comparison)
}
cbindDB=cbind(bF_RegCom,loop) # Change per variable
vec_i=as.vector(1) # Change per level
vec_j=as.vector(2) # Change per level
comparison.names=paste(t(combn(levels(comparison[drop=T]),2))[,1],t(combn(levels(comparison[drop=T]),2))[,2],sep="-")
Ncomp<-matrix(NA,ncol=3,nrow=length(vec_i),dimnames=list(c(comparison.names),c("vol1","vol2","volinter")))
for(i in 1:length(vec_i)){
  for(j in 1:3) {
    Ncomp[i,j]<-round(CHVintersect(as.matrix(cbindDB[cbindDB$loop==vec_i[i],c("C","N","P")]),
                                   as.matrix(cbindDB[cbindDB$loop==vec_j[i],c("C","N","P")]))$vol[j],2)
  }
}
Ncomp=data.frame(ComType=rownames(Ncomp),Ncomp)
bF_Reg_Vol=Ncomp%>%          # Change variable name                      
  mutate(voltot=vol1+vol2-volinter,                    
         vol.inter.rel=volinter/voltot*100,            
         minim=pmin(vol1,vol2),                          
         Turn=(2*minim-2*volinter)/(2*minim-volinter),
         Nest=(abs(vol1-vol2)/(vol1+vol2-volinter))*(volinter/(2*minim-volinter)))
mean(bF_Reg_Vol$vol.inter.rel);se(bF_Reg_Vol$vol.inter.rel) # mean and SE values of isotope volume # Change variable name   
mean(bF_Reg_Vol$Turn);se(bF_Reg_Vol$Turn) # Change variable name   
mean(bF_Reg_Vol$Nest);se(bF_Reg_Vol$Nest) # Change variable name   
bF_Reg_Vol_Stats <- data.frame(bF_Reg_Vol)
bF_Reg_Vol_Stats$Species <- "bF"

# The following lines allow calculating the niche volumes of each population relative to the total niche volume
# occupied by all   from all populations (i.e., Tot.vol)
bF_Reg_VolRel=data.frame(Site.Species=c(levels(comparison[drop=TRUE])),Volume=c(bF_Reg_Vol[1,2],bF_Reg_Vol[1:length(levels(comparison[drop=TRUE]))-1,3])) # Change variable
bF_Reg_VolRel=bF_Reg_VolRel%>% # Change variable
  mutate(vol.sp.rel=(round((Volume/bF_NicheSpace_TotVol)*100,2))) # Change variable
mean(bF_Reg_VolRel$vol.sp.rel);se(bF_Reg_VolRel$vol.sp.rel) # Change variable
bF_RegVolRel <- data.frame(bF_Reg_VolRel)
bF_RegVolRel$Species <- "bF"

# White birch regional
comparison=wB_RegCom$RegComType[drop=T] # Change per variable
for(i in 1:length(levels(comparison[drop=T]))){
  loop=ifelse(comparison==levels(comparison[drop=T])[i],i,comparison)
}
cbindDB=cbind(wB_RegCom,loop) # Change per variable
vec_i=as.vector(1) # Change per level
vec_j=as.vector(2) # Change per level
comparison.names=paste(t(combn(levels(comparison[drop=T]),2))[,1],t(combn(levels(comparison[drop=T]),2))[,2],sep="-")
Ncomp<-matrix(NA,ncol=3,nrow=length(vec_i),dimnames=list(c(comparison.names),c("vol1","vol2","volinter")))
for(i in 1:length(vec_i)){
  for(j in 1:3) {
    Ncomp[i,j]<-round(CHVintersect(as.matrix(cbindDB[cbindDB$loop==vec_i[i],c("C","N","P")]),
                                   as.matrix(cbindDB[cbindDB$loop==vec_j[i],c("C","N","P")]))$vol[j],2)
  }
}
Ncomp=data.frame(ComType=rownames(Ncomp),Ncomp)
wB_Reg_Vol=Ncomp%>%          # Change variable name                      
  mutate(voltot=vol1+vol2-volinter,                    
         vol.inter.rel=volinter/voltot*100,            
         minim=pmin(vol1,vol2),                          
         Turn=(2*minim-2*volinter)/(2*minim-volinter),
         Nest=(abs(vol1-vol2)/(vol1+vol2-volinter))*(volinter/(2*minim-volinter)))
mean(wB_Reg_Vol$vol.inter.rel);se(wB_Reg_Vol$vol.inter.rel) # mean and SE values of isotope volume # Change variable name   
mean(wB_Reg_Vol$Turn);se(wB_Reg_Vol$Turn) # Change variable name   
mean(wB_Reg_Vol$Nest);se(wB_Reg_Vol$Nest) # Change variable name   
wB_Reg_Vol_Stats <- data.frame(wB_Reg_Vol)
wB_Reg_Vol_Stats$Species <- "wB"
# The following lines allow calculating the niche volumes of each population relative to the total niche volume
# occupied by all   from all populations (i.e., Tot.vol)
wB_Reg_VolRel=data.frame(Site.Species=c(levels(comparison[drop=TRUE])),Volume=c(wB_Reg_Vol[1,2],wB_Reg_Vol[1:length(levels(comparison[drop=TRUE]))-1,3])) # Change variable
wB_Reg_VolRel=wB_Reg_VolRel%>% # Change variable
  mutate(vol.sp.rel=(round((Volume/wB_NicheSpace_TotVol)*100,2))) # Change variable
mean(wB_Reg_VolRel$vol.sp.rel);se(wB_Reg_VolRel$vol.sp.rel) # Change variable
wB_RegVolRel <- data.frame(wB_Reg_VolRel)
wB_RegVolRel$Species <- "wB"


# Balsam fir local spatial extent
comparison=bF_LocalCom$LocalComType[drop=T] # Change per variable
for(i in 1:length(levels(comparison[drop=T]))){
  loop=ifelse(comparison==levels(comparison[drop=T])[i],i,comparison)
}
cbindDB=cbind(bF_LocalCom,loop) # Change per variable
vec_i=as.vector(combn(1:length(levels(comparison[drop=T])),2)[1,])
vec_j=as.vector(combn(1:length(levels(comparison[drop=T])),2)[2,])
comparison.names=paste(t(combn(levels(comparison[drop=T]),2))[,1],t(combn(levels(comparison[drop=T]),2))[,2],sep="-")
Ncomp<-matrix(NA,ncol=3,nrow=length(vec_i),dimnames=list(c(comparison.names),c("vol1","vol2","volinter")))
for(i in 1:length(vec_i)){
  for(j in 1:3) {
    Ncomp[i,j]<-round(CHVintersect(as.matrix(cbindDB[cbindDB$loop==vec_i[i],c("C","N","P")]),
                                   as.matrix(cbindDB[cbindDB$loop==vec_j[i],c("C","N","P")]))$vol[j],2)
  }
}
Ncomp=data.frame(ComType=rownames(Ncomp),Ncomp)
bF_Local_Vol=Ncomp%>%          # Change variable name                      
  mutate(voltot=vol1+vol2-volinter,                    
         vol.inter.rel=volinter/voltot*100,            
         minim=pmin(vol1,vol2),                          
         Turn=(2*minim-2*volinter)/(2*minim-volinter),
         Nest=(abs(vol1-vol2)/(vol1+vol2-volinter))*(volinter/(2*minim-volinter)))
mean(bF_Local_Vol$vol.inter.rel);se(bF_Local_Vol$vol.inter.rel) # mean and SE values of isotope volume # Change variable name   
mean(bF_Local_Vol$Turn);se(bF_Local_Vol$Turn) # Change variable name   
mean(bF_Local_Vol$Nest);se(bF_Local_Vol$Nest) # Change variable name   
bF_Local_Vol_Stats <- data.frame(bF_Local_Vol) # Change variable
bF_Local_Vol_Stats$Species <- "bF"
# The following lines allow calculating the niche volumes of each population relative to the total niche volume
# occupied by all   from all populations (i.e., Tot.vol)
bF_Local_VolRel=data.frame(Site.Species=c(levels(comparison[drop=TRUE])),Volume=c(bF_Local_Vol[1,2],bF_Local_Vol[1:length(levels(comparison[drop=TRUE]))-1,3])) # Change variable
bF_Local_VolRel=bF_Local_VolRel%>% # Change variable
  mutate(vol.sp.rel=(round((Volume/bF_NicheSpace_TotVol)*100,2))) # Change variable
mean(bF_Local_VolRel$vol.sp.rel);se(bF_Local_VolRel$vol.sp.rel) # Change variable
bF_LocalVolRel <- data.frame(bF_Local_VolRel)
bF_LocalVolRel$Species <- "bF"

# White birch local spatial extent
comparison=wB_LocalCom$LocalComType[drop=T] # Change per variable
for(i in 1:length(levels(comparison[drop=T]))){
  loop=ifelse(comparison==levels(comparison[drop=T])[i],i,comparison)
}
cbindDB=cbind(wB_LocalCom,loop) # Change per variable
vec_i=as.vector(combn(1:length(levels(comparison[drop=T])),2)[1,])
vec_j=as.vector(combn(1:length(levels(comparison[drop=T])),2)[2,])
comparison.names=paste(t(combn(levels(comparison[drop=T]),2))[,1],t(combn(levels(comparison[drop=T]),2))[,2],sep="-")
Ncomp<-matrix(NA,ncol=3,nrow=length(vec_i),dimnames=list(c(comparison.names),c("vol1","vol2","volinter")))
for(i in 1:length(vec_i)){
  for(j in 1:3) {
    Ncomp[i,j]<-round(CHVintersect(as.matrix(cbindDB[cbindDB$loop==vec_i[i],c("C","N","P")]),
                                   as.matrix(cbindDB[cbindDB$loop==vec_j[i],c("C","N","P")]))$vol[j],2)
  }
}
Ncomp=data.frame(ComType=rownames(Ncomp),Ncomp)
wB_Local_Vol=Ncomp%>%          # Change variable name                      
  mutate(voltot=vol1+vol2-volinter,                    
         vol.inter.rel=volinter/voltot*100,            
         minim=pmin(vol1,vol2),                          
         Turn=(2*minim-2*volinter)/(2*minim-volinter),
         Nest=(abs(vol1-vol2)/(vol1+vol2-volinter))*(volinter/(2*minim-volinter)))
mean(wB_Local_Vol$vol.inter.rel);se(wB_Local_Vol$vol.inter.rel) # mean and SE values of isotope volume # Change variable name   
mean(wB_Local_Vol$Turn);se(wB_Local_Vol$Turn) # Change variable name   
mean(wB_Local_Vol$Nest);se(wB_Local_Vol$Nest) # Change variable name   
wB_Local_Vol_Stats <- data.frame(wB_Local_Vol) # Change variable
wB_Local_Vol_Stats$Species <- "wB"
# The following lines allow calculating the niche volumes of each population relative to the total niche volume
# occupied by all   from all populations (i.e., Tot.vol)
j <- wB_Local_Vol[1,2]
wB_Local_VolRel=data.frame(Site.Species=c(levels(comparison[drop=TRUE])),Volume=c(wB_Local_Vol[1,2],wB_Local_Vol[1:length(levels(comparison[drop=TRUE]))-1,3])) # Change variable
wB_Local_VolRel=wB_Local_VolRel%>% # Change variable
  mutate(vol.sp.rel=(round((Volume/wB_NicheSpace_TotVol)*100,2))) # Change variable
mean(wB_Local_VolRel$vol.sp.rel);se(wB_Local_VolRel$vol.sp.rel) # Change variable
wB_LocalVolRel <- data.frame(wB_Local_VolRel)
wB_LocalVolRel$Species <- "wB"

# Combine all results
Niche_Vol_Stats <- rbind(bF_Ecoreg_Vol_Stats, bF_Reg_Vol_Stats, bF_Local_Vol_Stats, wB_Ecoreg_Vol_Stats, wB_Reg_Vol_Stats, wB_Local_Vol_Stats)
write.csv(Niche_Vol_Stats, "Tables/StoicNiche_VolStats.csv")
Niche_VolRel_Stats <- rbind(bF_EcoregVolRel, bF_RegVolRel, bF_LocalVolRel, wB_EcoregVolRel, wB_RegVolRel, wB_LocalVolRel)
write.csv(Niche_VolRel_Stats, "Tables/StoicNiche_VolRelStats.csv")




###### Niche volume figures #######
# Create hypervolumne sphere figures using Gonzalez et. al., 2017 approach
library(rcdd)
library(geometry)
library(ade4)
library(manipulateWidget)
library(rgl)
library(vadr)
library(alphahull)
library(plyr)
library(dplyr)
library(tidyr)

# Use subsetted data from above - if starting from here; re-run load main dataset and data processing to get cons and hetero groups
FoliarCNP_bF <- subset(FoliarCNP, FoliarCNP$Species == "bF") # n=390
FoliarCNP_wB <- subset(FoliarCNP, FoliarCNP$Species == "wB") # n=229

# Run stoichiometric niche functions - R code before running the code below
# Create 3D Spheres

# Balsam fir ecoregion Northern Peninsula
bF_Ecoreg_NP_Niche=as.matrix(FoliarCNP_bF[FoliarCNP_bF[,"Region"]=="Northern Peninsula",c("C","N","P")])
Abies_Ecoreg_NP.Niche <- CHVind(bF_Ecoreg_NP_Niche)
# Balsam fir ecoregion Central Forest
bF_Ecoreg_CF_Niche=as.matrix(FoliarCNP_bF[FoliarCNP_bF[,"Region"]=="Central Forest",c("C","N","P")])
Abies_Ecoreg_CF.Niche <- CHVind(bF_Ecoreg_CF_Niche)
# Balsam fir - Regional extent conspecific
bF_O_Niche=as.matrix(FoliarCNP_bF[FoliarCNP_bF[,"RegComType"]=="Conspecific",c("C","N","P")])
Abies_O.Niche <- CHVind(bF_O_Niche)
# Balsam fir - Local extent conspecific Northern Peninsula population
bF_O_NP_Niche=as.matrix(FoliarCNP_bF[FoliarCNP_bF[,"LocalComType"]=="NP: Conspecific",c("C","N","P")])
Abies_O_NP.Niche <- CHVind(bF_O_NP_Niche)
# Balsam fir - Local extent conspecific Central Forest population
bF_O_CNLF_Niche=as.matrix(FoliarCNP_bF[FoliarCNP_bF[,"LocalComType"]=="CF: Conspecific",c("C","N","P")])
Abies_O_CNLF.Niche <- CHVind(bF_O_CNLF_Niche)
str(Abies_O_CNLF.Niche)
# Balsam fir - Regional extent heterospecific
bF_Co_Niche=as.matrix(FoliarCNP_bF[FoliarCNP_bF[,"RegComType"]=="Heterospecific",c("C","N","P")])
Abies_Co.Niche <- CHVind(bF_Co_Niche)
# Balsam fir - Local extent heterospecific Northern Peninsula population
bF_Co_NP_Niche=as.matrix(FoliarCNP_bF[FoliarCNP_bF[,"LocalComType"]=="NP: Heterospecific",c("C","N","P")])
Abies_Co_NP.Niche <- CHVind(bF_Co_NP_Niche)
# Balsam fir - Local extent heterospecific Central Forest population
bF_Co_CNLF_Niche=as.matrix(FoliarCNP_bF[FoliarCNP_bF[,"LocalComType"]=="CF: Heterospecific",c("C","N","P")])
Abies_Co_CNLF.Niche <- CHVind(bF_Co_CNLF_Niche)

# White birch
# White birch ecoregion Northern Peninsula
wB_Ecoreg_NP_Niche=as.matrix(FoliarCNP_wB[FoliarCNP_wB[,"Region"]=="Northern Peninsula",c("C","N","P")])
Betula_Ecoreg_NP.Niche <- CHVind(wB_Ecoreg_NP_Niche)
# White birch ecoregion Central Forest
wB_Ecoreg_CF_Niche=as.matrix(FoliarCNP_wB[FoliarCNP_wB[,"Region"]=="Central Forest",c("C","N","P")])
Betula_Ecoreg_CF.Niche <- CHVind(wB_Ecoreg_CF_Niche)
# White birch - Regional extent conspecific
wB_O_Niche=as.matrix(FoliarCNP_wB[FoliarCNP_wB[,"RegComType"]=="Conspecific",c("C","N","P")])
Betula_O.Niche <- CHVind(wB_O_Niche)
# White birch - Local extent conspecific Northern Peninsula populations
wB_O_NP_Niche=as.matrix(FoliarCNP_wB[FoliarCNP_wB[,"LocalComType"]=="NP: Conspecific",c("C","N","P")])
Betula_O_NP.Niche <- CHVind(wB_O_NP_Niche)
# White birch - Local extent conspecific Central Forest population
wB_O_CNLF_Niche=as.matrix(FoliarCNP_wB[FoliarCNP_wB[,"LocalComType"]=="CF: Conspecific",c("C","N","P")])
Betula_O_CNLF.Niche <- CHVind(wB_O_CNLF_Niche)

# White birch - Regional extent heterospecific
wB_Co_Niche=as.matrix(FoliarCNP_wB[FoliarCNP_wB[,"RegComType"]=="Heterospecific",c("C","N","P")])
Betula_Co.Niche <- CHVind(wB_Co_Niche)
# White birch - Local extent heterospecific Northern Peninsula population
wB_Co_NP_Niche=as.matrix(FoliarCNP_wB[FoliarCNP_wB[,"LocalComType"]=="NP: Heterospecific",c("C","N","P")])
Betula_Co_NP.Niche <- CHVind(wB_Co_NP_Niche)
# White birch - Local extent heterospecific Central Forest population
wB_Co_CNLF_Niche=as.matrix(FoliarCNP_wB[FoliarCNP_wB[,"LocalComType"]=="CF: Heterospecific",c("C","N","P")])
Betula_Co_CNLF.Niche <- CHVind(wB_Co_CNLF_Niche)

##### 3D Plots ######
# Balsam fir and White birch ecoregion 3D plot
Balsam_EcoregNP <- as.data.frame(Abies_Ecoreg_NP.Niche$vert_set1)
mean(Balsam_EcoregNP$V1)
mean(Balsam_EcoregNP$V2)
mean(Balsam_EcoregNP$V3)

Balsam_EcoregCF <- as.data.frame(Abies_Ecoreg_CF.Niche$vert_set1)
mean(Balsam_EcoregCF$V1)
mean(Balsam_EcoregCF$V2)
mean(Balsam_EcoregCF$V3)

Betula_EcoregNP <- as.data.frame(Betula_Ecoreg_NP.Niche$vert_set1)
mean(Betula_EcoregNP$V1)
mean(Betula_EcoregNP$V2)
mean(Betula_EcoregNP$V3)

Betula_EcoregCF <- as.data.frame(Betula_Ecoreg_CF.Niche$vert_set1)
mean(Betula_EcoregCF$V1)
mean(Betula_EcoregCF$V2)
mean(Betula_EcoregCF$V3)

# Set total range of data for 3d plot
bF_TotalNicheRange <- subset(FoliarCNP_bF, Species =="bF", select = c("Region", "C","N","P"))
wB_TotalNicheRange <- subset(FoliarCNP_wB, Species=="wB", select = c("Region", "C","N","P"))


# Create plot
mfrow3d(1,2)
rgl.clear()
par3d(windowRect = c(400, 100, 1600, 700), cex=1.2)

EcoregCol <- c("Northern Peninsula" = "#9A5324FF", "Central Forest" = "#767676FF") # Ecoregions

# Balsam fir plot ecoregion
plot3d(bF_TotalNicheRange[,2:4], size=2, type='p', axes=TRUE, edges="bbox", col=EcoregCol[bF_TotalNicheRange$Region], repel=F, xlab=" ", ylab=" ", zlab=" ")
axes3d(edges="bbox") # Remove box face
title3d(main = "Between ecoregions", line=5, level=2)
title3d( main="(a) Balsam fir", line=4) # Add title
# Add axis labels
mtext3d(text="Carbon", edge="x++", line=2.4)
mtext3d(text="Nitrogen", edge="y--", line=2)
mtext3d(text="Phosphorus", edge="z-+", line=-2) 
# Add spheres
spheres3d(apply(Abies_Ecoreg_NP.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Abies_Ecoreg_NP.Niche$vol[1]*3)/(4*pi))^(1/3),color="#9A5324FF")
spheres3d(apply(Abies_Ecoreg_CF.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Abies_Ecoreg_CF.Niche$vol[1]*3)/(4*pi))^(1/3),color="#767676FF")
# Add ablines
abclines3d(52.10333,0.9873333, 0.1120667, a=diag(3), col="#9A5324FF")
abclines3d(52.127,0.919, 0.0893, a=diag(3), col="#767676FF")


# White birch plot regional
plot3d(wB_TotalNicheRange[,2:4], size=2, type='p', axes=TRUE, edges="bbox", col=EcoregCol[wB_TotalNicheRange$Region], repel=F, xlab=" ", ylab=" ", zlab=" ")
axes3d(edges="bbox") # Remove box face
title3d( main="(b) White birch", line=4) # Add title
# Add axis labels
mtext3d(text="Carbon", edge="x++", line=2.4)
mtext3d(text="Nitrogen", edge="y--", line=2)
mtext3d(text="Phosphorus", edge="z-+", line=-2) 
# Add spheres
spheres3d(apply(Betula_Ecoreg_NP.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Betula_Ecoreg_NP.Niche$vol[1]*3)/(4*pi))^(1/3),color="#9A5324FF")
spheres3d(apply(Betula_Ecoreg_CF.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Betula_Ecoreg_CF.Niche$vol[1]*3)/(4*pi))^(1/3),color="#767676FF")
# Add drop lines
abclines3d(49.93636,2.82, 0.2937727, a=diag(3), col="#9A5324FF")
abclines3d(50.22778,1.747222, 0.1692222, a=diag(3), col="#767676FF")
# Add legend
legend3d("topright", legend = c("Northern Peninsula","Central Forest"), col = c("#9A5324FF","#767676FF"), pch=16, inset=0.07, cex=1.5)

# Save figure - need to find a way to make/set resolution (good as is)
snapshot3d("Figures/bFwB_BetweenEcoregions_Niche.png",fmt="png", top=TRUE) #To save the plot3D
rgl.postscript("Figures/bFwB_BetweenEcoregions_Niche.ps", fmt = 'ps')





# Balsam fir and white birch regional 3D plot
# Extract coordinates for drop lines - numbers added below in plot
Balsam_O <- as.data.frame(Abies_O.Niche$vert_set1)
mean(Balsam_O$V1)
mean(Balsam_O$V2)
mean(Balsam_O$V3)

Balsam_Co <- as.data.frame(Abies_Co.Niche$vert_set1)
mean(Balsam_Co$V1)
mean(Balsam_Co$V2)
mean(Balsam_Co$V3)

# Extract coordinates for drop lines
Birch_O <- as.data.frame(Betula_O.Niche$vert_set1)
mean(Birch_O$V1)
mean(Birch_O$V2)
mean(Birch_O$V3)

Birch_Co <- as.data.frame(Betula_Co.Niche$vert_set1)
mean(Birch_Co$V1)
mean(Birch_Co$V2)
mean(Birch_Co$V3)


# Regional spatial extent plots
# Set total range of data for 3d plot
bF_TotalNicheRange <- subset(FoliarCNP_bF, Species =="bF", select = c("RegComType", "C","N","P"))
wB_TotalNicheRange <- subset(FoliarCNP_wB, Species=="wB", select = c("RegComType", "C","N","P"))

# Create plot
mfrow3d(1,2)
rgl.clear()
par3d(windowRect = c(400, 100, 1600, 700), cex=1.2)

RegCol <- c("Conspecific" = "#3E3E23FF", "Heterospecific" = "#CC8214FF") # Regional extent

# Balsam fir plot regional
plot3d(bF_TotalNicheRange[,2:4], size=2, type='p', axes=TRUE, edges="bbox", col=RegCol[bF_TotalNicheRange$RegComType], repel=F, xlab=" ", ylab=" ", zlab=" ")
axes3d(edges="bbox") # Remove box face
title3d(main = "Community types across ecoregions", line=5, level=2)
title3d( main="(a) Balsam fir", line=4) # Add title
# Add axis labels
mtext3d(text="Carbon", edge="x++", line=2.4)
mtext3d(text="Nitrogen", edge="y--", line=2)
mtext3d(text="Phosphorus", edge="z-+", line=-2) 
# Add spheres
spheres3d(apply(Abies_O.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Abies_O.Niche$vol[1]*3)/(4*pi))^(1/3),color="#3E3E23FF")
spheres3d(apply(Abies_Co.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Abies_Co.Niche$vol[1]*3)/(4*pi))^(1/3),color="#CC8214FF")
# Add ablines
abclines3d(51.774,0.9384, 0.11112, a=diag(3), col="#3E3E23FF")
abclines3d(52.29,1.050833, 0.119875, a=diag(3), col="#CC8214FF")

# White birch plot regional
plot3d(wB_TotalNicheRange[,2:4], size=2, type='p', axes=TRUE, edges="bbox", col=RegCol[wB_TotalNicheRange$RegComType], repel=F, xlab=" ", ylab=" ", zlab=" ")
axes3d(edges="bbox") # Remove box face
title3d( main="(b) White birch", line=4) # Add title
# Add axis labels
mtext3d(text="Carbon", edge="x++", line=2.4)
mtext3d(text="Nitrogen", edge="y--", line=2)
mtext3d(text="Phosphorus", edge="z-+", line=-2) 
# Add spheres
spheres3d(apply(Betula_O.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Betula_O.Niche$vol[1]*3)/(4*pi))^(1/3),color="#3E3E23FF")
spheres3d(apply(Betula_Co.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Betula_Co.Niche$vol[1]*3)/(4*pi))^(1/3),color="#CC8214FF")
# Add drop lines
abclines3d(49.93273,2.231818, 0.2344545, a=diag(3), col="#3E3E23FF")
abclines3d(50.21409,2.541818, 0.2669091, a=diag(3), col="#CC8214FF")
# Add legend
legend3d("topright", legend = c("Conspecific","Heterospecific"), col = c("#3E3E23FF","#CC8214FF"), pch=16, inset=0.07, cex=1.5)

# Save figure - need to find a way to make/set resolution (good as is)
snapshot3d("Figures/bFwB_AcrossEcoregions_Niche.png",fmt="png", top=TRUE) #To save the plot3D
#rgl.postscript("Figures/bFwB_RegionalNiche.pdf",fmt="pdf", top=TRUE) #To save the plot3D
rgl.postscript("Figures/bFwB_AcrossEcoregions_Niche.ps", fmt = 'ps')



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Balsam fir and white birch local 3D plot
# Extract coordinates for drop lines - numbers added below in plot
# Put mean volume and standard errors plus differene between NP and CF into a table
# Balsam fir
Balsam_O_NP <- as.data.frame(Abies_O_NP.Niche$vert_set1)
bF_Con_NP_C <- paste0(round(mean(Balsam_O_NP$V1), digits=3), " ", "\\u00b1 ",round(se(Balsam_O_NP$V1), digits=3))
bF_Con_NP_N <- paste0(round(mean(Balsam_O_NP$V2), digits=3), " ", "\\u00b1 ",round(se(Balsam_O_NP$V2), digits=3))
bF_Con_NP_P <- paste0(round(mean(Balsam_O_NP$V3), digits=3), " ", "\\u00b1 ",round(se(Balsam_O_NP$V3), digits=3))

Balsam_Co_NP <- as.data.frame(Abies_Co_NP.Niche$vert_set1)
bF_Hetero_NP_C <- paste0(round(mean(Balsam_Co_NP$V1), digits=3), " ", "\\u00b1 ",round(se(Balsam_Co_NP$V1), digits=3))
bF_Hetero_NP_N <- paste0(round(mean(Balsam_Co_NP$V2), digits=3), " ", "\\u00b1 ",round(se(Balsam_Co_NP$V2), digits=3))
bF_Hetero_NP_P <- paste0(round(mean(Balsam_Co_NP$V3), digits=3), " ", "\\u00b1 ",round(se(Balsam_Co_NP$V3), digits=3))

Balsam_O_CF <- as.data.frame(Abies_O_CNLF.Niche$vert_set1)
bF_Con_CF_C <- paste0(round(mean(Balsam_O_CF$V1), digits=3), " ", "\\u00b1 ",round(se(Balsam_O_CF$V1), digits=3))
bF_Con_CF_N <- paste0(round(mean(Balsam_O_CF$V2), digits=3), " ", "\\u00b1 ",round(se(Balsam_O_CF$V2), digits=3))
bF_Con_CF_P <- paste0(round(mean(Balsam_O_CF$V3), digits=3), " ", "\\u00b1 ",round(se(Balsam_O_CF$V3), digits=3))

Balsam_Co_CF <- as.data.frame(Abies_Co_CNLF.Niche$vert_set1)
bF_Hetero_CF_C <- paste0(round(mean(Balsam_Co_CF$V1), digits=3), " ", "\\u00b1 ",round(se(Balsam_Co_CF$V1), digits=3))
bF_Hetero_CF_N <- paste0(round(mean(Balsam_Co_CF$V2), digits=3), " ", "\\u00b1 ",round(se(Balsam_Co_CF$V2), digits=3))
bF_Hetero_CF_P <- paste0(round(mean(Balsam_Co_CF$V3), digits=3), " ", "\\u00b1 ",round(se(Balsam_Co_CF$V3), digits=3))

# Get differences between means
bF_Con_C_Diff <- paste0(round(mean(Balsam_O_NP$V1) - mean(Balsam_O_CF$V1), digits=3))
bF_Con_N_Diff <- paste0(round(mean(Balsam_O_NP$V2) - mean(Balsam_O_CF$V2), digits=3))
bF_Con_P_Diff <- paste0(round(mean(Balsam_O_NP$V3) - mean(Balsam_O_CF$V3), digits=3))

bF_Hetero_C_Diff <- paste0(round(mean(Balsam_Co_NP$V1) - mean(Balsam_Co_CF$V1), digits=3))
bF_Hetero_N_Diff <- paste0(round(mean(Balsam_Co_NP$V2) - mean(Balsam_Co_CF$V2), digits=3))
bF_Hetero_P_Diff <- paste0(round(mean(Balsam_Co_NP$V3) - mean(Balsam_Co_CF$V3), digits=3))

# White birch
Birch_O_NP <- as.data.frame(Betula_O_NP.Niche$vert_set1)
wB_Con_NP_C <- paste0(round(mean(Birch_O_NP$V1), digits=3), " ", "\\u00b1 ",round(se(Birch_O_NP$V1), digits=3))
wB_Con_NP_N <- paste0(round(mean(Birch_O_NP$V2), digits=3), " ", "\\u00b1 ",round(se(Birch_O_NP$V2), digits=3))
wB_Con_NP_P <- paste0(round(mean(Birch_O_NP$V3), digits=3), " ", "\\u00b1 ",round(se(Birch_O_NP$V3), digits=3))


Birch_Co_NP <- as.data.frame(Betula_Co_NP.Niche$vert_set1)
wB_Hetero_NP_C <- paste0(round(mean(Birch_Co_NP$V1), digits=3), " ", "\\u00b1 ",round(se(Birch_Co_NP$V1), digits=3))
wB_Hetero_NP_N <- paste0(round(mean(Birch_Co_NP$V2), digits=3), " ", "\\u00b1 ",round(se(Birch_Co_NP$V2), digits=3))
wB_Hetero_NP_P <- paste0(round(mean(Birch_Co_NP$V3), digits=3), " ", "\\u00b1 ",round(se(Birch_Co_NP$V3), digits=3))


Birch_O_CF <- as.data.frame(Betula_O_CNLF.Niche$vert_set1)
wB_Con_CF_C <- paste0(round(mean(Birch_O_CF$V1), digits=3), " ", "\\u00b1 ",round(se(Birch_O_CF$V1), digits=3))
wB_Con_CF_N <- paste0(round(mean(Birch_O_CF$V2), digits=3), " ", "\\u00b1 ",round(se(Birch_O_CF$V2), digits=3))
wB_Con_CF_P <- paste0(round(mean(Birch_O_CF$V3), digits=3), " ", "\\u00b1 ",round(se(Birch_O_CF$V3), digits=3))


Birch_Co_CF <- as.data.frame(Betula_Co_CNLF.Niche$vert_set1)
wB_Hetero_CF_C <- paste0(round(mean(Birch_Co_CF$V1), digits=3), " ", "\\u00b1 ",round(se(Birch_Co_CF$V1), digits=3))
wB_Hetero_CF_N <- paste0(round(mean(Birch_Co_CF$V2), digits=3), " ", "\\u00b1 ",round(se(Birch_Co_CF$V2), digits=3))
wB_Hetero_CF_P <- paste0(round(mean(Birch_Co_CF$V3), digits=3), " ", "\\u00b1 ",round(se(Birch_Co_CF$V3), digits=3))

# Get differences between means
wB_Con_C_Diff <- paste0(round(mean(Birch_O_NP$V1) - mean(Birch_O_CF$V1), digits=3))
wB_Con_N_Diff <- paste0(round(mean(Birch_O_NP$V2) - mean(Birch_O_CF$V2), digits=3))
wB_Con_P_Diff <- paste0(round(mean(Birch_O_NP$V3) - mean(Birch_O_CF$V3), digits=3))

wB_Hetero_C_Diff <- paste0(round(mean(Birch_Co_NP$V1) - mean(Birch_Co_CF$V1), digits=3))
wB_Hetero_N_Diff <- paste0(round(mean(Birch_Co_NP$V2) - mean(Birch_Co_CF$V2), digits=3))
wB_Hetero_P_Diff <- paste0(round(mean(Birch_Co_NP$V3) - mean(Birch_Co_CF$V3), digits=3))

# Balsam fir
bF_VolAvg_NP <- data.frame(cbind(bF_Con_NP_C, bF_Con_NP_N, bF_Con_NP_P), cbind(bF_Hetero_NP_C, bF_Hetero_NP_N, bF_Hetero_NP_P))
colnames(bF_VolAvg_NP) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
bF_VolAvg_NP$SiteName <- "Local NP"

bF_VolAvg_CF <- data.frame(cbind(bF_Con_CF_C, bF_Con_CF_N, bF_Con_CF_P), cbind(bF_Hetero_CF_C, bF_Hetero_CF_N, bF_Hetero_CF_P))
colnames(bF_VolAvg_CF) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
bF_VolAvg_CF$SiteName <- "Local CF"

bF_VolAvg_Diff <- data.frame(cbind(bF_Con_C_Diff, bF_Con_N_Diff, bF_Con_P_Diff), cbind(bF_Hetero_C_Diff, bF_Hetero_N_Diff, bF_Hetero_P_Diff))
colnames(bF_VolAvg_Diff) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
bF_VolAvg_Diff$SiteName <- "Difference"

bF_VolAvg <- data.frame(rbind(bF_VolAvg_NP, bF_VolAvg_CF, bF_VolAvg_Diff))

# White birch
wB_VolAvg_NP <- data.frame(cbind(wB_Con_NP_C, wB_Con_NP_N, wB_Con_NP_P), cbind(wB_Hetero_NP_C, wB_Hetero_NP_N, wB_Hetero_NP_P))
colnames(wB_VolAvg_NP) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
wB_VolAvg_NP$SiteName <- "Local NP"

wB_VolAvg_CF <- data.frame(cbind(wB_Con_CF_C, wB_Con_CF_N, wB_Con_CF_P), cbind(wB_Hetero_CF_C, wB_Hetero_CF_N, wB_Hetero_CF_P))
colnames(wB_VolAvg_CF) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
wB_VolAvg_CF$SiteName <- "Local CF"

wB_VolAvg_Diff <- data.frame(cbind(wB_Con_C_Diff, wB_Con_N_Diff, wB_Con_P_Diff), cbind(wB_Hetero_C_Diff, wB_Hetero_N_Diff, wB_Hetero_P_Diff))
colnames(wB_VolAvg_Diff) <- c("Con_C", "Con_N", "Con_P", "Hetero_C", "Hetero_N", "Hetero_P")
wB_VolAvg_Diff$SiteName <- "Difference"

wB_VolAvg <- data.frame(rbind(wB_VolAvg_NP, wB_VolAvg_CF, wB_VolAvg_Diff))

Species <- data.frame(c("bF", "bF", "bF", "wB", "wB", "wB"))
colnames(Species) <- c("Species")

Niche_Vol_Diff <- data.frame(cbind(Species, rbind(bF_VolAvg, wB_VolAvg)))

write.csv(Niche_Vol_Diff, "Tables/NiceVol_Diff.csv")

# Local spatial extent plots
# Set total range of data for 3d plot
bF_TotalNicheRange_Local <- subset(FoliarCNP_bF, Species =="bF", select = c("LocalComType", "C","N","P"))
wB_TotalNicheRange_Local <- subset(FoliarCNP_wB, Species=="wB", select = c("LocalComType", "C","N","P"))

# Create plot
mfrow3d(1,2)
rgl.clear()
par3d(windowRect = c(400, 100, 1600, 700), cex=1.2)

LocalCol <- c("NP: Conspecific" = "#5B8FA8FF", "NP: Heterospecific" = "#616530FF", "CF: Conspecific" = "#0F425CFF", "CF: Heterospecific" = "#9A5324FF") # Local extent

# Balsam fir plot local
plot3d(bF_TotalNicheRange_Local[,2:4], size=2, type='p', axes=TRUE, edges="bbox", col=LocalCol[bF_TotalNicheRange_Local$LocalComType], repel=F, xlab=" ", ylab=" ", zlab=" ")
axes3d(edges="bbox") # Remove box face
title3d(main = "Community types within and between ecoregions", line=5, level=2)
title3d( main="(a) Balsam fir", line=4) # Add title
# Add axis labels
mtext3d(text="Carbon", edge="x++", line=2.4)
mtext3d(text="Nitrogen", edge="y--", line=3)
mtext3d(text="Phosphorus", edge="z-+", line=-2) 
# Add spheres
spheres3d(apply(Abies_O_NP.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Abies_O_NP.Niche$vol[1]*3)/(4*pi))^(1/3),color="#5B8FA8FF")
spheres3d(apply(Abies_Co_NP.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Abies_Co_NP.Niche$vol[1]*3)/(4*pi))^(1/3),color="#616530FF")
spheres3d(apply(Abies_O_CNLF.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Abies_O_CNLF.Niche$vol[1]*3)/(4*pi))^(1/3),color="#0F425CFF")
spheres3d(apply(Abies_Co_CNLF.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Abies_Co_CNLF.Niche$vol[1]*3)/(4*pi))^(1/3),color="#9A5324FF")
# Add drop lines
# Northern Peninsula population
abclines3d(51.6913,0.986087, 0.1126522, a=diag(3), col="#5B8FA8FF")
abclines3d(52.18182,1.060455, 0.1279091, a=diag(3), col="#616530FF")
# Central Forest population
abclines3d(52.23643,0.8764286, 0.08964286, a=diag(3), col="#0F425CFF")
abclines3d(52.18765,0.9441176, 0.08452941, a=diag(3), col="#9A5324FF")

# White birch loal plot
plot3d(wB_TotalNicheRange_Local[,2:4], size=2, type='p', axes=TRUE, edges="bbox", col=LocalCol[wB_TotalNicheRange_Local$LocalComType], repel=F, xlab=" ", ylab=" ", zlab=" ")
axes3d(edges="bbox") # Remove box face
title3d( main="(b) White birch", line=4) # Add title
# Add axis labels
mtext3d(text="Carbon", edge="x++", line=2.4)
mtext3d(text="Nitrogen", edge="y--", line=2)
mtext3d(text="Phosphorus", edge="z-+", line=-2) 
# Add spheres
spheres3d(apply(Betula_O_NP.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Betula_O_NP.Niche$vol[1]*3)/(4*pi))^(1/3),color="#5B8FA8FF")
spheres3d(apply(Betula_Co_NP.Niche$vert_set1,2,mean), alpha=0.7, radius = ((Betula_Co_NP.Niche$vol[1]*3)/(4*pi))^(1/3),color="#616530FF")
spheres3d(apply(Betula_O_CNLF.Niche$vert_set1,2,mean), radius = ((Betula_O_CNLF.Niche$vol[1]*3)/(4*pi))^(1/3),color="#0F425CFF")
spheres3d(apply(Betula_Co_CNLF.Niche$vert_set1,2,mean), radius = ((Betula_Co_CNLF.Niche$vol[1]*3)/(4*pi))^(1/3),color="#9A5324FF")
# Add drop lines
# Northern Peninsula population
abclines3d(49.86,3.07, 0.3342, a=diag(3), col="#5B8FA8FF")
abclines3d(49.93636,2.82, 0.2937727, a=diag(3), col="#616530FF")
# Central Forest population
abclines3d(49.84889,1.852222, 0.1875556, a=diag(3), col="#0F425CFF")
abclines3d(50.62267,1.555333, 0.1358667, a=diag(3), col="#9A5324FF")

# Add legend
legend3d("topright", legend = c("Northern Peninsula:", "Conspecific", "Heterospecific", "Central Forest:", "Conspecific", "Heterospecific"), col = c("white", "#5B8FA8FF","#616530FF","white", "#0F425CFF", "#9A5324FF"), pch=16, inset=0.025, cex=1.3)

# Save figure - need to find a way to make/set resolution (good as is)
snapshot3d("Figures/bFwB_WithinBetweenEcoregions_Niche.png",fmt="png", top=TRUE) #To save the plot3D
rgl.postscript("Figures/bFwB_WithinBetweenEcoregions_Niche.ps", fmt = 'ps')




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###### Process temporal data ######
# Temporal comparison
library(tidyr)
Elemental2017 <- read.csv("Data/OriginalData/TNNP/2017/CNP2017_Merge_FromLab.csv")
Elemental2017$Species <- as.character(Elemental2017$Species)
Elemental2017$PlotID <- as.character(Elemental2017$PlotID)
FoliarCNP_2017 <- subset(Elemental2017, Elemental2017$Species == "ABBA" | Elemental2017$Species == "BEPA")
str(FoliarCNP_2017)
FoliarCNP_2017$PlotID <- as.factor(FoliarCNP_2017$PlotID)
FoliarCNP_2017 = FoliarCNP_2017 %>% 
  group_by(PlotID) %>% # Goup by PlotID
  mutate(Hetero = n() > 1) # If name occurs more than once = TRUE
FoliarCNP_2017$Hetero <- as.character(FoliarCNP_2017$Hetero) # Set field as character so we can get "TRUE"
FoliarCNP_2017$LocalComType <- ifelse(FoliarCNP_2017$Hetero == "TRUE", "CF: Heterospecific", "CF: Conspecific") # Set hetero and cons field
FoliarCNP_2017$RegComType <- ifelse(FoliarCNP_2017$Hetero == "TRUE", "Heterospecific", "Conspecific") # Set hetero and cons field
FoliarCNP_2017$Region <- "Central Forest"
FoliarCNP_2017 <- as.data.frame(subset(FoliarCNP_2017, select=c("Species", "Region", "LocalComType", "RegComType", "TC", "N", "P")))
colnames(FoliarCNP_2017) <- c("Species", "Region", "LocalComType", "RegComType", "C", "N", "P") # Convert column names to match 
FoliarCNP_2017$Year <- "2017"
FoliarCNP_2017$Temporal <- paste0(FoliarCNP_2017$RegComType,"-", FoliarCNP_2017$Year)
str(FoliarCNP_2017)
# Pull data from above
FoliarCNP$Year <- "2016"
FoliarCNP_2016 <- subset(FoliarCNP, FoliarCNP$Region == "Central Forest") 
FoliarCNP_2016$Temporal <- paste0(FoliarCNP_2016$RegComType,"-", FoliarCNP_2016$Year)

str(FoliarCNP_2016)
# Combine 2016 and 2017 data
TemporalCNP <- rbind(FoliarCNP_2016, FoliarCNP_2017)

# Check 2017 data structure
ABBA_FoliarCNP_2017 <- subset(FoliarCNP_2017, Species == "ABBA")
str(ABBA_FoliarCNP_2017)
ABBA_FoliarCNP_2017_Con <- subset(ABBA_FoliarCNP_2017, LocalComType == "CF: Conspecific")
str(ABBA_FoliarCNP_2017_Con)
ABBA_FoliarCNP_2017_Hetero <- subset(ABBA_FoliarCNP_2017, LocalComType == "CF: Heterospecific")
str(ABBA_FoliarCNP_2017_Hetero)

BEPA_FoliarCNP_2017 <- subset(FoliarCNP_2017, Species == "BEPA")
str(BEPA_FoliarCNP_2017)
BEPA_FoliarCNP_2017_Con <- subset(BEPA_FoliarCNP_2017, LocalComType == "CF: Conspecific")
str(BEPA_FoliarCNP_2017_Con)
BEPA_FoliarCNP_2017_Hetero <- subset(BEPA_FoliarCNP_2017, LocalComType == "CF: Heterospecific")
str(BEPA_FoliarCNP_2017_Hetero)


# Set factor levels for plotting
TemporalCNP$Temporal <- factor(TemporalCNP$Temporal, levels = c("Conspecific-2016","Heterospecific-2016", "Conspecific-2017","Heterospecific-2017"))

# Subset into species
bF_TemporalCNP <- subset(TemporalCNP, Species=="bF" | Species == "ABBA")
ABBA_2016 <- subset(bF_TemporalCNP, Year == "2016")
str(ABBA_2016)
wB_TemporalCNP <- subset(TemporalCNP, Species=="wB" | Species == "BEPA")
BEPA_2016 <- subset(wB_TemporalCNP, Year == "2016")
str(BEPA_2016)

###### Temporal perform PCA #####
library(FactoMineR)
library(factoextra)
bF_Temp_PCA <- PCA(bF_TemporalCNP[5:7], scale.unit = TRUE, ncp = 2, graph = TRUE)
wB_Temp_PCA <- PCA(wB_TemporalCNP[5:7], scale.unit = TRUE, ncp = 2, graph = TRUE)

# Create PCA biplot
TempCol <- c("Conspecific-2016" = "orangered1", "Heterospecific-2016" = "darkgoldenrod1", "Conspecific-2017" = "darkolivegreen4", "Heterospecific-2017" = "burlywood4")
# PCA 

bF_PCA_Temp <- fviz_pca_biplot(bF_Temp_PCA, col.ind = bF_TemporalCNP$Temporal, palette = TempCol, addEllipses = TRUE, ellipse.level=0.95, label = "var", col.var = "black", repel = TRUE, legend.title = "Community Type", title = "(a) Balsam fir temporal", legend = "none") + theme(text = element_text(size = 11)) 

wB_PCA_Temp <- fviz_pca_biplot(wB_Temp_PCA, col.ind = wB_TemporalCNP$Temporal, palette = TempCol, addEllipses = TRUE, ellipse.level=0.95, label = "var", col.var = "black", repel = TRUE, legend.title = "Community Type", title = "(b) White birch temporal", legend = "none") + theme(text = element_text(size = 11)) 

library(gridExtra)
library(ggpubr)
# Export PCA plots
tiff("Figures/Temporal_PCA.tiff", width = 7, height = 6, units = "in", res = 600) 
ggarrange(bF_PCA_Temp, wB_PCA_Temp, ncol=2, nrow=1, common.legend = T, legend = "bottom")
dev.off()

###### Temporal scatter plot ######
# bF
bF_Temporal_CN <- ggplot(bF_TemporalCNP, aes(x=N, y=C, col=Temporal)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "N") + theme_linedraw() + scale_color_manual(name = "Community Types", values=TempCol) + ggtitle("Balsam fir temporal\\n(a)") + theme(text = element_text(size=11))

bF_Temporal_CP <- ggplot(bF_TemporalCNP, aes(x=P, y=C, col=Temporal)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=TempCol) + ggtitle("\\n(b)") + theme(text = element_text(size=11))

bF_Temporal_NP <- ggplot(bF_TemporalCNP, aes(x=P, y=N, col=Temporal)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "N", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=TempCol) + ggtitle("\\n(c)") + theme(text = element_text(size=11))

# wB
wB_Temporal_CN <- ggplot(wB_TemporalCNP, aes(x=N, y=C, col=Temporal)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "N") + theme_linedraw() + scale_color_manual(name = "Community Types", values=TempCol) + ggtitle("White birch temporal\\n(d)") + theme(text = element_text(size=11))

wB_Temporal_CP <- ggplot(wB_TemporalCNP, aes(x=P, y=C, col=Temporal)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "C", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=TempCol) + ggtitle("\\n(e)") + theme(text = element_text(size=11))

wB_Temporal_NP <- ggplot(wB_TemporalCNP, aes(x=P, y=N, col=Temporal)) +geom_point(size=0.5) + stat_density_2d(size= 0.5) + labs(y = "N", x = "P") + theme_linedraw() + scale_color_manual(name = "Community Types", values=TempCol) + ggtitle("\\n(f)") + theme(text = element_text(size=11))

tiff("Figures/Scatter_Temporal.tiff", width = 8, height = 8, units = "in", res = 600) 
ggarrange(bF_Temporal_CN, bF_Temporal_CP, bF_Temporal_NP, wB_Temporal_CN, wB_Temporal_CP, wB_Temporal_NP, ncol=3, nrow=2, common.legend = T, legend = "bottom")
dev.off()


###### Temporal hypervolume calculaution ######
library("hypervolume")
# balsam fir
HV_bF_Temporal_Cons2016 <- hypervolume_gaussian(bF_TemporalCNP[which(bF_TemporalCNP$Temporal == "Conspecific-2016"),names(bF_TemporalCNP) %in% c("C","N","P")])
HV_bF_Temporal_Hetero2016 <- hypervolume_gaussian(bF_TemporalCNP[which(bF_TemporalCNP$Temporal == "Heterospecific-2016"),names(bF_TemporalCNP) %in% c("C","N","P")])
HV_bF_Temporal_Cons2017 <- hypervolume_gaussian(bF_TemporalCNP[which(bF_TemporalCNP$Temporal == "Conspecific-2017"),names(bF_TemporalCNP) %in% c("C","N","P")])
HV_bF_Temporal_Hetero2017 <- hypervolume_gaussian(bF_TemporalCNP[which(bF_TemporalCNP$Temporal == "Heterospecific-2017"),names(bF_TemporalCNP) %in% c("C","N","P")])
HV_bF_Temporal_2016 <- hypervolume_gaussian(bF_TemporalCNP[which(bF_TemporalCNP$Year == "2016"),names(bF_TemporalCNP) %in% c("C","N","P")])
HV_bF_Temporal_2017 <- hypervolume_gaussian(bF_TemporalCNP[which(bF_TemporalCNP$Year == "2017"),names(bF_TemporalCNP) %in% c("C","N","P")])

# White  birch
HV_wB_Temporal_Cons2016 <- hypervolume_gaussian(wB_TemporalCNP[which(wB_TemporalCNP$Temporal == "Conspecific-2016"),names(wB_TemporalCNP) %in% c("C","N","P")])
HV_wB_Temporal_Hetero2016 <- hypervolume_gaussian(wB_TemporalCNP[which(wB_TemporalCNP$Temporal == "Heterospecific-2016"),names(wB_TemporalCNP) %in% c("C","N","P")])
HV_wB_Temporal_Cons2017 <- hypervolume_gaussian(wB_TemporalCNP[which(wB_TemporalCNP$Temporal == "Conspecific-2017"),names(wB_TemporalCNP) %in% c("C","N","P")])
HV_wB_Temporal_Hetero2017 <- hypervolume_gaussian(wB_TemporalCNP[which(wB_TemporalCNP$Temporal == "Heterospecific-2017"),names(wB_TemporalCNP) %in% c("C","N","P")])
HV_wB_Temporal_2016 <- hypervolume_gaussian(wB_TemporalCNP[which(wB_TemporalCNP$Year == "2016"),names(wB_TemporalCNP) %in% c("C","N","P")])
HV_wB_Temporal_2017 <- hypervolume_gaussian(wB_TemporalCNP[which(wB_TemporalCNP$Year == "2017"),names(wB_TemporalCNP) %in% c("C","N","P")])

# Hypervolume comparison
# Group hypervolumnes at the local level: white  birch
HVS_bF_Temporal_Cons <- hypervolume_set(HV_bF_Temporal_Cons2016, HV_bF_Temporal_Cons2017, check.memory=FALSE)
HVS_bF_Temporal_Hetero <- hypervolume_set(HV_bF_Temporal_Hetero2016, HV_bF_Temporal_Hetero2017, check.memory=FALSE)
HVS_bF_Temporal_2017 <- hypervolume_set(HV_bF_Temporal_Cons2017, HV_bF_Temporal_Hetero2017, check.memory=FALSE)
HVS_bF_Temporal <- hypervolume_set(HV_bF_Temporal_2016, HV_bF_Temporal_2017, check.memory=FALSE)

HVS_wB_Temporal_Cons <- hypervolume_set(HV_wB_Temporal_Cons2016, HV_wB_Temporal_Cons2017, check.memory=FALSE)
HVS_wB_Temporal_Hetero <- hypervolume_set(HV_wB_Temporal_Hetero2016, HV_wB_Temporal_Hetero2017, check.memory=FALSE)
HVS_wB_Temporal_2017 <- hypervolume_set(HV_wB_Temporal_Cons2017, HV_wB_Temporal_Hetero2017, check.memory=FALSE)
HVS_wB_Temporal <- hypervolume_set(HV_wB_Temporal_2016, HV_wB_Temporal_2017, check.memory=FALSE)


# Overlap statistics
# Balsam fir
HVS_bF_Temporal_Cons_Stats <- data.frame(hypervolume_overlap_statistics(HVS_bF_Temporal_Cons))
colnames(HVS_bF_Temporal_Cons_Stats) <- c("Index")
HVS_bF_Temporal_Cons_Stats$Comparison <- "bF: temporal con"

HVS_bF_Temporal_Hetero_Stats <- data.frame(hypervolume_overlap_statistics(HVS_bF_Temporal_Hetero))
colnames(HVS_bF_Temporal_Hetero_Stats) <- c("Index")
HVS_bF_Temporal_Hetero_Stats$Comparison <- "bF: temporal hetero"

HVS_bF_Temporal_2017_Stats <- data.frame(hypervolume_overlap_statistics(HVS_bF_Temporal_2017))
colnames(HVS_bF_Temporal_2017_Stats) <- c("Index")
HVS_bF_Temporal_2017_Stats$Comparison <- "bF: con vs hetero 2017"

HVS_bF_Temporal_Stats <- data.frame(hypervolume_overlap_statistics(HVS_bF_Temporal))
colnames(HVS_bF_Temporal_Stats) <- c("Index")
HVS_bF_Temporal_Stats$Comparison <- "bF: temporal by year"



# White birch
HVS_wB_Temporal_Cons_Stats <- data.frame(hypervolume_overlap_statistics(HVS_wB_Temporal_Cons))
colnames(HVS_wB_Temporal_Cons_Stats) <- c("Index")
HVS_wB_Temporal_Cons_Stats$Comparison <- "wB: temporal con"

HVS_wB_Temporal_Hetero_Stats <- data.frame(hypervolume_overlap_statistics(HVS_wB_Temporal_Hetero))
colnames(HVS_wB_Temporal_Hetero_Stats) <- c("Index")
HVS_wB_Temporal_Hetero_Stats$Comparison <- "wB: temporal hetero"

HVS_wB_Temporal_2017_Stats <- data.frame(hypervolume_overlap_statistics(HVS_wB_Temporal_2017))
colnames(HVS_wB_Temporal_2017_Stats) <- c("Index")
HVS_wB_Temporal_2017_Stats$Comparison <- "wB: con vs hetero 2017"

HVS_wB_Temporal_Stats <- data.frame(hypervolume_overlap_statistics(HVS_wB_Temporal))
colnames(HVS_wB_Temporal_Stats) <- c("Index")
HVS_wB_Temporal_Stats$Comparison <- "wB: temporal by year"


library("data.table")
# Create single table of hypervolumn results
HV_TemporalStatList <- list(HVS_bF_Temporal_Cons_Stats, HVS_bF_Temporal_Hetero_Stats, HVS_bF_Temporal_2017_Stats, HVS_bF_Temporal_Stats, HVS_wB_Temporal_Cons_Stats, HVS_wB_Temporal_Hetero_Stats, HVS_wB_Temporal_2017_Stats, HVS_wB_Temporal_Stats) # List of results in data frame
HV_TemporalStatsData <- rbindlist(lapply(HV_TemporalStatList, setDT, keep.rownames = TRUE)) # Bind lists together keeping row names
names(HV_TemporalStatsData)[names(HV_TemporalStatsData) == "rn"] <- "Statistic"

# Transpose for statistics as column names and comparison groups as row names
library("tidyr")
HV_TemporalStatsData_Spread <- spread(HV_TemporalStatsData, Statistic, Index)
HV_TemporalStatsData_Sub <- subset(HV_TemporalStatsData_Spread, select=c("Comparison", "jaccard", "sorensen", "frac_unique_1", "frac_unique_2"))
# Create vector with desired row order
HV_TemporalOrder <- c("bF: temporal con", "bF: temporal hetero", "bF: con vs hetero 2017", "bF: temporal by year", "wB: temporal con", "wB: temporal hetero", "wB: con vs hetero 2017", "wB: temporal by year" )
# re-order by HV_Order
HV_TemporalStats_Data_Final <- HV_TemporalStatsData_Sub %>% mutate(Comparison =  factor(Comparison, levels = HV_TemporalOrder)) %>% arrange(Comparison)
# Export results
write.csv(HV_TemporalStats_Data_Final, "Tables/Hypervolume_TemrporalStats.csv")


###### Temporal multivariate normality test #####
library(mvShapiroTest)
# Balsam fir
bF_Temporal_Cons2017 <- mvShapiro.Test(data.matrix(bF_TemporalCNP[which(bF_TemporalCNP$Temporal == "Conspecific-2017"),names(bF_TemporalCNP) %in% c("C","N","P")]))
bF_Temporal_Cons2017_Tab <- data.frame(c(bF_Temporal_Cons2017$statistic, bF_Temporal_Cons2017$p.value))
bF_Temporal_Cons2017_Tab$Source <- "bF: con 2017"
row.names(bF_Temporal_Cons2017_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Temporal_Cons2017_Tab$Metric <- row.names(bF_Temporal_Cons2017_Tab) 
colnames(bF_Temporal_Cons2017_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Temporal_Cons2017_Tab) <- c()

bF_Temporal_Hetero2017 <- mvShapiro.Test(data.matrix(bF_TemporalCNP[which(bF_TemporalCNP$Temporal == "Heterospecific-2017"),names(bF_TemporalCNP) %in% c("C","N","P")]))
bF_Temporal_Hetero2017_Tab <- data.frame(c(bF_Temporal_Hetero2017$statistic, bF_Temporal_Hetero2017$p.value))
bF_Temporal_Hetero2017_Tab$Source <- "bF: hetero 2017"
row.names(bF_Temporal_Hetero2017_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Temporal_Hetero2017_Tab$Metric <- row.names(bF_Temporal_Hetero2017_Tab) 
colnames(bF_Temporal_Hetero2017_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Temporal_Hetero2017_Tab) <- c()

bF_Temporal_2016 <- mvShapiro.Test(data.matrix(bF_TemporalCNP[which(bF_TemporalCNP$Year == "2016"),names(bF_TemporalCNP) %in% c("C","N","P")]))
bF_Temporal_2016_Tab <- data.frame(c(bF_Temporal_2016$statistic, bF_Temporal_2016$p.value))
bF_Temporal_2016_Tab$Source <- "bF: 2016"
row.names(bF_Temporal_2016_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Temporal_2016_Tab$Metric <- row.names(bF_Temporal_2016_Tab) 
colnames(bF_Temporal_2016_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Temporal_2016_Tab) <- c()

bF_Temporal_2017 <- mvShapiro.Test(data.matrix(bF_TemporalCNP[which(bF_TemporalCNP$Year == "2017"),names(bF_TemporalCNP) %in% c("C","N","P")]))
bF_Temporal_2017_Tab <- data.frame(c(bF_Temporal_2017$statistic, bF_Temporal_2017$p.value))
bF_Temporal_2017_Tab$Source <- "bF: 2017"
row.names(bF_Temporal_2017_Tab) <- c("Shapiro_Wilk", "p_value")
bF_Temporal_2017_Tab$Metric <- row.names(bF_Temporal_2017_Tab) 
colnames(bF_Temporal_2017_Tab) <- c("Values", "Source", "Metric") 
rownames(bF_Temporal_2017_Tab) <- c()

# White birch
wB_Temporal_Cons2017 <- mvShapiro.Test(data.matrix(wB_TemporalCNP[which(wB_TemporalCNP$Temporal == "Conspecific-2017"),names(wB_TemporalCNP) %in% c("C","N","P")]))
wB_Temporal_Cons2017_Tab <- data.frame(c(wB_Temporal_Cons2017$statistic, wB_Temporal_Cons2017$p.value))
wB_Temporal_Cons2017_Tab$Source <- "wB: con 2017"
row.names(wB_Temporal_Cons2017_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Temporal_Cons2017_Tab$Metric <- row.names(wB_Temporal_Cons2017_Tab) 
colnames(wB_Temporal_Cons2017_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Temporal_Cons2017_Tab) <- c()

wB_Temporal_Hetero2017 <- mvShapiro.Test(data.matrix(wB_TemporalCNP[which(wB_TemporalCNP$Temporal == "Heterospecific-2017"),names(wB_TemporalCNP) %in% c("C","N","P")]))
wB_Temporal_Hetero2017_Tab <- data.frame(c(wB_Temporal_Hetero2017$statistic, wB_Temporal_Hetero2017$p.value))
wB_Temporal_Hetero2017_Tab$Source <- "wB: hetero 2017"
row.names(wB_Temporal_Hetero2017_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Temporal_Hetero2017_Tab$Metric <- row.names(wB_Temporal_Hetero2017_Tab) 
colnames(wB_Temporal_Hetero2017_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Temporal_Hetero2017_Tab) <- c()

wB_Temporal_2016 <- mvShapiro.Test(data.matrix(wB_TemporalCNP[which(wB_TemporalCNP$Year == "2016"),names(wB_TemporalCNP) %in% c("C","N","P")]))
wB_Temporal_2016_Tab <- data.frame(c(wB_Temporal_2016$statistic, wB_Temporal_2016$p.value))
wB_Temporal_2016_Tab$Source <- "wB: 2016"
row.names(wB_Temporal_2016_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Temporal_2016_Tab$Metric <- row.names(wB_Temporal_2016_Tab) 
colnames(wB_Temporal_2016_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Temporal_2016_Tab) <- c()

wB_Temporal_2017 <- mvShapiro.Test(data.matrix(wB_TemporalCNP[which(wB_TemporalCNP$Year == "2017"),names(wB_TemporalCNP) %in% c("C","N","P")]))
wB_Temporal_2017_Tab <- data.frame(c(wB_Temporal_2017$statistic, wB_Temporal_2017$p.value))
wB_Temporal_2017_Tab$Source <- "wB: 2017"
row.names(wB_Temporal_2017_Tab) <- c("Shapiro_Wilk", "p_value")
wB_Temporal_2017_Tab$Metric <- row.names(wB_Temporal_2017_Tab) 
colnames(wB_Temporal_2017_Tab) <- c("Values", "Source", "Metric") 
rownames(wB_Temporal_2017_Tab) <- c()


# Collect results
MVN_TemporalData <- rbind(bF_Temporal_Hetero2017_Tab, bF_Temporal_2016_Tab, bF_Temporal_2017_Tab, wB_Temporal_Cons2017_Tab, wB_Temporal_Hetero2017_Tab, wB_Temporal_2016_Tab, wB_Temporal_2017_Tab)

library("tidyr")
MVN_TemporalData_Spread <- spread(MVN_TemporalData, Metric, Values)
MVN_TemporalData_SpreadSub <- subset(MVN_TemporalData_Spread, select=c("Source", "Shapiro_Wilk", "p_value"))
# Create vector with desired row order
MVN_TemporalOrder <- c("bF: hetero 2017", "bF: 2016", "bF: 2017", "wB: con 2017", "wB: hetero 2017", "wB: 2016", "wB: 2017")
# re-order by HV_Order
MVN_TemporalData_Final <- MVN_TemporalData_SpreadSub %>% mutate(Source =  factor(Source, levels = MVN_TemporalOrder)) %>% arrange(Source)

write.csv(MVN_TemporalData_Final, "Tables/MVN_TemporalStats.csv")


###### Temporal HMD and PERMANOVA #######
# PERMANOVA comparison 
# Subset data by community type
bF_TemporalCons <- subset(bF_TemporalCNP, bF_TemporalCNP$RegComType == "Conspecific")
bF_TemporalHetero <- subset(bF_TemporalCNP, bF_TemporalCNP$RegComType == "Heterospecific")
bF_TemporalYear2017 <- subset(bF_TemporalCNP, bF_TemporalCNP$Year == "2017")

wB_TemporalCons <- subset(wB_TemporalCNP, wB_TemporalCNP$RegComType == "Conspecific")
wB_TemporalHetero <- subset(wB_TemporalCNP, wB_TemporalCNP$RegComType == "Heterospecific")
wB_TemporalYear2017 <- subset(wB_TemporalCNP, wB_TemporalCNP$Year == "2017")

library(vegan)
# Compare regional conspecific and heterospecific groups between years
# Balsam fir
# 2016 vs 2017 conspecific comparison
PM_bF_TemporalCons_Dispr <- betadisper(d=vegdist(bF_TemporalCons[5:7]), group=bF_TemporalCons$Year)
permutest(PM_bF_TemporalCons_Dispr)
PM_bF_TemporalCons_Stats <- data.frame(adonis2(bF_TemporalCons[5:7] ~ Year, data = bF_TemporalCons))
PM_bF_TemporalCons_Stats$Source <- row.names(PM_bF_TemporalCons_Stats) 
PM_bF_TemporalCons_Stats$Source[PM_bF_TemporalCons_Stats$Source == "Year"] <- "bF: temporal con"

# 2016 vs 2017 heterospecific comparison
PM_bF_TemporalHetero_Dispr <- betadisper(d=vegdist(bF_TemporalHetero[5:7]), group=bF_TemporalHetero$Year)
permutest(PM_bF_TemporalHetero_Dispr)
PM_bF_TemporalHetero_Stats <- data.frame(adonis2(bF_TemporalHetero[5:7] ~ Year, data = bF_TemporalHetero))
PM_bF_TemporalHetero_Stats$Source <- row.names(PM_bF_TemporalHetero_Stats) 
PM_bF_TemporalHetero_Stats$Source[PM_bF_TemporalHetero_Stats$Source == "Year"] <- "bF: temporal hetero"

# 2017 conspecific vs heterospecific comparison
PM_bF_ConsHetero2017_Dispr <- betadisper(d=vegdist(bF_TemporalYear2017[5:7]), group=bF_TemporalYear2017$RegComType)
permutest(PM_bF_ConsHetero2017_Dispr)
PM_bF_ConsHetero2017_Stats <- data.frame(adonis2(bF_TemporalYear2017[5:7] ~ RegComType, data = bF_TemporalYear2017))
PM_bF_ConsHetero2017_Stats$Source <- row.names(PM_bF_ConsHetero2017_Stats) 
PM_bF_ConsHetero2017_Stats$Source[PM_bF_ConsHetero2017_Stats$Source == "RegComType"] <- "bF: 2017 con vs hetero"

# 2016 vs 2017 local comparison
PM_bF_2016_2017_Dispr <- betadisper(d=vegdist(bF_TemporalCNP[5:7]), group=bF_TemporalCNP$Year)
permutest(PM_bF_2016_2017_Dispr)
PM_bF_2016_2017_Stats <- data.frame(adonis2(bF_TemporalCNP[5:7] ~ Year, data = bF_TemporalCNP))
PM_bF_2016_2017_Stats$Source <- row.names(PM_bF_2016_2017_Stats) 
PM_bF_2016_2017_Stats$Source[PM_bF_2016_2017_Stats$Source == "Year"] <- "bF: 2016 vs 2017"


# White birch
# 2016 vs 2017 conspecific comparison
PM_wB_TemporalCons_Dispr <- betadisper(d=vegdist(wB_TemporalCons[5:7]), group=wB_TemporalCons$Year)
permutest(PM_wB_TemporalCons_Dispr)
PM_wB_TemporalCons_Stats <- data.frame(adonis2(wB_TemporalCons[5:7] ~ Year, data = wB_TemporalCons))
PM_wB_TemporalCons_Stats$Source <- row.names(PM_wB_TemporalCons_Stats) 
PM_wB_TemporalCons_Stats$Source[PM_wB_TemporalCons_Stats$Source == "Year"] <- "wB: temporal con"


# 2016 vs 2017 heterospecific comparison
PM_wB_TemporalHetero_Dispr <- betadisper(d=vegdist(wB_TemporalHetero[5:7]), group=wB_TemporalHetero$Year)
permutest(PM_wB_TemporalHetero_Dispr)
PM_wB_TemporalHetero_Stats <- data.frame(adonis2(wB_TemporalHetero[5:7] ~ Year, data = wB_TemporalHetero))
PM_wB_TemporalHetero_Stats$Source <- row.names(PM_wB_TemporalHetero_Stats) 
PM_wB_TemporalHetero_Stats$Source[PM_wB_TemporalHetero_Stats$Source == "Year"] <- "wB: temporal hetero"

# 2017 conspecific vs heterospecific comparison
PM_wB_ConsHetero2017_Dispr <- betadisper(d=vegdist(wB_TemporalYear2017[5:7]), group=wB_TemporalYear2017$RegComType)
permutest(PM_wB_ConsHetero2017_Dispr)
PM_wB_TemporalConsHetero2017_Stats <- data.frame(adonis2(wB_TemporalYear2017[5:7] ~ RegComType, data = wB_TemporalYear2017))
PM_wB_TemporalConsHetero2017_Stats$Source <- row.names(PM_wB_TemporalConsHetero2017_Stats) 
PM_wB_TemporalConsHetero2017_Stats$Source[PM_wB_TemporalConsHetero2017_Stats$Source == "RegComType"] <- "wB: 2017 con vs hetero"

# 2016 vs 2017 local comparison
PM_wB_2016_2017_Dispr <- betadisper(d=vegdist(wB_TemporalCNP[5:7]), group=wB_TemporalCNP$Year)
permutest(PM_wB_2016_2017_Dispr)
PM_wB_2016_2017_Stats1 <- data.frame(adonis2(wB_TemporalCNP[5:7] ~ Year, data = wB_TemporalCNP))
PM_wB_2016_2017_Stats$Source <- row.names(PM_wB_2016_2017_Stats) 
PM_wB_2016_2017_Stats$Source[PM_wB_2016_2017_Stats$Source == "Year"] <- "wB: 2016 vs 2017"

# List results by species
bF_TemporalList <- list(PM_bF_TemporalCons_Stats, PM_bF_TemporalHetero_Stats, PM_bF_ConsHetero2017_Stats, PM_bF_2016_2017_Stats)
wB_TemporalList <- list(PM_wB_TemporalCons_Stats, PM_wB_TemporalHetero_Stats, PM_wB_TemporalConsHetero2017_Stats, PM_wB_2016_2017_Stats)

# Process table
bF_TemporalRow <- bind_rows(bF_TemporalList) # bind list together
bF_TemporalRows <- bF_TemporalRow %>% mutate_if(is.numeric, round, digits=4) # round numeric data to 3 decimal places
bF_TemporalRowsOrder <- subset(bF_TemporalRows, select=c("Source", "Df", "SumOfSqs", "R2", "F", "Pr..F.")) # Subset data - change column order
colnames(bF_TemporalRowsOrder) <- c("Source", "Df", "SS", "R2", "F", "p-value") # Rename column names

wB_TemporalRow <- bind_rows(wB_TemporalList) # bind list together
wB_TemporalRows <- wB_TemporalRow %>% mutate_if(is.numeric, round, digits=4) # round numeric data to 3 decimal places
wB_TemporalRowsOrder <- subset(wB_TemporalRows, select=c("Source", "Df", "SumOfSqs", "R2", "F", "Pr..F.")) # Subset data - change column order
colnames(wB_TemporalRowsOrder) <- c("Source", "Df", "SS", "R2", "F", "p-value") # Rename column names

# Combine column
PM_TemporalResults <- cbind(bF_TemporalRowsOrder, wB_TemporalRowsOrder)
View(PM_TemporalResults)

# Export results
write.csv(PM_TemporalResults, "Tables/PERMANOVA_TemporalStats_v1.csv")



###### Temporal niche volume metrics #######
# Gozalez metrics for calculating niche metrics
library(rcdd)
library(geometry)
library(ade4)
library(rgl)
library(vadr)
library(alphahull)
library(plyr)
library(dplyr)
library(tidyr)


# Balsam fir subset for processing
bF_TemporalConHetero <- subset(bF_TemporalCNP, select = c("Temporal","C","N","P" ))  # Temporal compare of con vs hetero
# Subset further to remove conspecific 2017 group - not enough data to compute hypervolume metrics
bF_TemporalConHetero <- subset(bF_TemporalConHetero, Temporal == "Conspecific-2016" | Temporal == "Heterospecific-2016" | Temporal == "Heterospecific-2017")
bF_Year <- subset(bF_TemporalCNP, select = c("Year","C","N","P" )) # Temporal compare between years
bF_Year$Year <- as.factor(bF_Year$Year)

# White birch subset for processing
wB_TemporalConHetero <- subset(wB_TemporalCNP, select = c("Temporal","C","N","P" ))  # Regional extent
str(wB_TemporalConHetero)
wB_TemporalConHetero$Temporal <- as.factor(wB_TemporalConHetero$Temporal)
wB_Year <- subset(wB_TemporalCNP, select = c("Year","C","N","P" )) # Local extent
wB_Year$Year <- as.factor(wB_Year$Year)

# Calculate niche total available volumes
bF_NicheSpaceTemporal <- as.matrix(bF_TemporalCNP[,c("C","N","P")])
bF_NicheSpaceTemporal_TotVol <- CHVind(bF_NicheSpaceTemporal)
bF_NicheSpaceTemporal_TotVol <- as.numeric(bF_NicheSpaceTemporal_TotVol$vol)

wB_NicheSpaceTemporal <- as.matrix(wB_TemporalCNP[,c("C","N","P")])
wB_NicheSpaceTemporal_TotVol <- CHVind(wB_NicheSpaceTemporal)
wB_NicheSpaceTemporal_TotVol <- as.numeric(wB_NicheSpaceTemporal_TotVol$vol)


# Balsam fir temporal compare con vs hetero
# This code generates a file that contains vol1, vol2, and volinter. Vol1= volume of the first population, vol2 = volume of the second population, and volinter = overlap between the two niches.
comparison=bF_TemporalConHetero$Temporal[drop=T] # Change per variable
for(i in 1:length(levels(comparison[drop=T]))){
  loop=ifelse(comparison==levels(comparison[drop=T])[i],i,comparison)
}
cbindDB=cbind(bF_TemporalConHetero,loop) # Change per variable
vec_i=as.vector(combn(1:length(levels(comparison[drop=T])),2)[1,])
vec_j=as.vector(combn(1:length(levels(comparison[drop=T])),2)[2,])
comparison.names=paste(t(combn(levels(comparison[drop=T]),2))[,1],t(combn(levels(comparison[drop=T]),2))[,2],sep="-")
Ncomp<-matrix(NA,ncol=3,nrow=length(vec_i),dimnames=list(c(comparison.names),c("vol1","vol2","volinter")))
for(i in 1:length(vec_i)){
  for(j in 1:3) {
    Ncomp[i,j]<-round(CHVintersect(as.matrix(cbindDB[cbindDB$loop==vec_i[i],c("C","N","P")]),
                                   as.matrix(cbindDB[cbindDB$loop==vec_j[i],c("C","N","P")]))$vol[j],2)
  }
}
Ncomp=data.frame(ComType=rownames(Ncomp),Ncomp)
bF_TemporalConHetero_Vol=Ncomp%>%          # Change variable name                      
  mutate(voltot=vol1+vol2-volinter,                    
         vol.inter.rel=volinter/voltot*100,            
         minim=pmin(vol1,vol2),                          
         Turn=(2*minim-2*volinter)/(2*minim-volinter),
         Nest=(abs(vol1-vol2)/(vol1+vol2-volinter))*(volinter/(2*minim-volinter)))
mean(bF_TemporalConHetero_Vol$vol.inter.rel);se(bF_TemporalConHetero_Vol$vol.inter.rel) # mean and SE values of isotope volume # Change variable name   
mean(bF_TemporalConHetero_Vol$Turn);se(bF_TemporalConHetero_Vol$Turn) # Change variable name   
mean(bF_TemporalConHetero_Vol$Nest);se(bF_TemporalConHetero_Vol$Nest) # Change variable name   
bF_TemporalConHetero_Vol_Stats <- data.frame(bF_TemporalConHetero_Vol)
bF_TemporalConHetero_Vol_Stats$Species <- "bF"
#View(bF_TemporalConHetero_Vol_Stats)
# The following lines allow calculating the niche volumes of each population relative to the total niche volume
# occupied by all   from all populations (i.e., Tot.vol)
bF_TemporalConHetero_VolRel=data.frame(Site.Species=c(levels(comparison[drop=TRUE])),Volume=c(bF_TemporalConHetero_Vol[1,2],bF_TemporalConHetero_Vol[1:length(levels(comparison[drop=TRUE]))-1,3])) # Change variable
bF_TemporalConHetero_VolRel=bF_TemporalConHetero_VolRel%>% # Change variable
  mutate(vol.sp.rel=(round((Volume/bF_NicheSpaceTemporal_TotVol)*100,2))) # Change variable
mean(bF_TemporalConHetero_VolRel$vol.sp.rel);se(bF_TemporalConHetero_VolRel$vol.sp.rel) # Change variable
bF_TemporalConHetero_VolRel_data <- data.frame(bF_TemporalConHetero_VolRel)
bF_TemporalConHetero_VolRel_data$Species <- "bF"


# Balsam fir temporal between year
# This code generates a file that contains vol1, vol2, and volinter. Vol1= volume of the first population, vol2 = volume of the second population, and volinter = overlap between the two niches.
comparison=bF_Year$Year[drop=T] # Change per variable
for(i in 1:length(levels(comparison[drop=T]))){
  loop=ifelse(comparison==levels(comparison[drop=T])[i],i,comparison)
}
cbindDB=cbind(bF_Year,loop) # Change per variable
vec_i=as.vector(1) # Change per level
vec_j=as.vector(2) # Change per level
comparison.names=paste(t(combn(levels(comparison[drop=T]),2))[,1],t(combn(levels(comparison[drop=T]),2))[,2],sep="-")
Ncomp<-matrix(NA,ncol=3,nrow=length(vec_i),dimnames=list(c(comparison.names),c("vol1","vol2","volinter")))
for(i in 1:length(vec_i)){
  for(j in 1:3) {
    Ncomp[i,j]<-round(CHVintersect(as.matrix(cbindDB[cbindDB$loop==vec_i[i],c("C","N","P")]),
                                   as.matrix(cbindDB[cbindDB$loop==vec_j[i],c("C","N","P")]))$vol[j],2)
  }
}
Ncomp=data.frame(ComType=rownames(Ncomp),Ncomp)
bF_TemporalYear_Vol=Ncomp%>%          # Change variable name                      
  mutate(voltot=vol1+vol2-volinter,                    
         vol.inter.rel=volinter/voltot*100,            
         minim=pmin(vol1,vol2),                          
         Turn=(2*minim-2*volinter)/(2*minim-volinter),
         Nest=(abs(vol1-vol2)/(vol1+vol2-volinter))*(volinter/(2*minim-volinter)))
mean(bF_TemporalYear_Vol$vol.inter.rel);se(bF_TemporalYear_Vol$vol.inter.rel) # mean and SE values of isotope volume # Change variable name   
mean(bF_TemporalYear_Vol$Turn);se(bF_TemporalYear_Vol$Turn) # Change variable name   
mean(bF_TemporalYear_Vol$Nest);se(bF_TemporalYear_Vol$Nest) # Change variable name   
bF_TemporalYear_Vol_Stats <- data.frame(bF_TemporalYear_Vol)
bF_TemporalYear_Vol_Stats$Species <- "bF"
# The following lines allow calculating the niche volumes of each population relative to the total niche volume
# occupied by all   from all populations (i.e., Tot.vol)
bF_TemporalYear_VolRel=data.frame(Site.Species=c(levels(comparison[drop=TRUE])),Volume=c(bF_TemporalYear_Vol[1,2],bF_TemporalYear_Vol[1:length(levels(comparison[drop=TRUE]))-1,3])) # Change variable
bF_TemporalYear_VolRel=bF_TemporalYear_VolRel%>% # Change variable
  mutate(vol.sp.rel=(round((Volume/bF_NicheSpaceTemporal_TotVol)*100,2))) # Change variable
mean(bF_TemporalYear_VolRel$vol.sp.rel);se(bF_TemporalYear_VolRel$vol.sp.rel) # Change variable
bF_TemporalYear_VolRel_data <- data.frame(bF_TemporalYear_VolRel)
bF_TemporalYear_VolRel_data$Species <- "bF"



# White birch temporal compare con vs hetero
# This code generates a file that contains vol1, vol2, and volinter. Vol1= volume of the first population, vol2 = volume of the second population, and volinter = overlap between the two niches.
comparison=wB_TemporalConHetero$Temporal[drop=T] # Change per variable
for(i in 1:length(levels(comparison[drop=T]))){
  loop=ifelse(comparison==levels(comparison[drop=T])[i],i,comparison)
}
cbindDB=cbind(wB_TemporalConHetero,loop) # Change per variable
vec_i=as.vector(combn(1:length(levels(comparison[drop=T])),2)[1,])
vec_j=as.vector(combn(1:length(levels(comparison[drop=T])),2)[2,])
comparison.names=paste(t(combn(levels(comparison[drop=T]),2))[,1],t(combn(levels(comparison[drop=T]),2))[,2],sep="-")
Ncomp<-matrix(NA,ncol=3,nrow=length(vec_i),dimnames=list(c(comparison.names),c("vol1","vol2","volinter")))
for(i in 1:length(vec_i)){
  for(j in 1:3) {
    Ncomp[i,j]<-round(CHVintersect(as.matrix(cbindDB[cbindDB$loop==vec_i[i],c("C","N","P")]),
                                   as.matrix(cbindDB[cbindDB$loop==vec_j[i],c("C","N","P")]))$vol[j],2)
  }
}
Ncomp=data.frame(ComType=rownames(Ncomp),Ncomp)
wB_TemporalConHetero_Vol=Ncomp%>%          # Change variable name                      
  mutate(voltot=vol1+vol2-volinter,                    
         vol.inter.rel=volinter/voltot*100,            
         minim=pmin(vol1,vol2),                          
         Turn=(2*minim-2*volinter)/(2*minim-volinter),
         Nest=(abs(vol1-vol2)/(vol1+vol2-volinter))*(volinter/(2*minim-volinter)))
mean(wB_TemporalConHetero_Vol$vol.inter.rel);se(wB_TemporalConHetero_Vol$vol.inter.rel) # mean and SE values of isotope volume # Change variable name   
mean(wB_TemporalConHetero_Vol$Turn);se(wB_TemporalConHetero_Vol$Turn) # Change variable name   
mean(wB_TemporalConHetero_Vol$Nest);se(wB_TemporalConHetero_Vol$Nest) # Change variable name   
wB_TemporalConHetero_Vol_Stats <- data.frame(wB_TemporalConHetero_Vol)
wB_TemporalConHetero_Vol_Stats$Species <- "wB"
# The following lines allow calculating the niche volumes of each population relative to the total niche volume
# occupied by all   from all populations (i.e., Tot.vol)
wB_TemporalConHetero_VolRel=data.frame(Site.Species=c(levels(comparison[drop=TRUE])),Volume=c(wB_TemporalConHetero_Vol[1,2],wB_TemporalConHetero_Vol[1:length(levels(comparison[drop=TRUE]))-1,3])) # Change variable
wB_TemporalConHetero_VolRel=wB_TemporalConHetero_VolRel%>% # Change variable
  mutate(vol.sp.rel=(round((Volume/wB_NicheSpaceTemporal_TotVol)*100,2))) # Change variable
mean(wB_TemporalConHetero_VolRel$vol.sp.rel);se(wB_TemporalConHetero_VolRel$vol.sp.rel) # Change variable
wB_TemporalConHetero_VolRel_data <- data.frame(wB_TemporalConHetero_VolRel)
wB_TemporalConHetero_VolRel_data$Species <- "wB"


# White birch temporal between year
# This code generates a file that contains vol1, vol2, and volinter. Vol1= volume of the first population, vol2 = volume of the second population, and volinter = overlap between the two niches.
comparison=wB_Year$Year[drop=T] # Change per variable
for(i in 1:length(levels(comparison[drop=T]))){
  loop=ifelse(comparison==levels(comparison[drop=T])[i],i,comparison)
}
cbindDB=cbind(wB_Year,loop) # Change per variable
vec_i=as.vector(1) # Change per level
vec_j=as.vector(2) # Change per level
comparison.names=paste(t(combn(levels(comparison[drop=T]),2))[,1],t(combn(levels(comparison[drop=T]),2))[,2],sep="-")
Ncomp<-matrix(NA,ncol=3,nrow=length(vec_i),dimnames=list(c(comparison.names),c("vol1","vol2","volinter")))
for(i in 1:length(vec_i)){
  for(j in 1:3) {
    Ncomp[i,j]<-round(CHVintersect(as.matrix(cbindDB[cbindDB$loop==vec_i[i],c("C","N","P")]),
                                   as.matrix(cbindDB[cbindDB$loop==vec_j[i],c("C","N","P")]))$vol[j],2)
  }
}
Ncomp=data.frame(ComType=rownames(Ncomp),Ncomp)
wB_TemporalYear_Vol=Ncomp%>%          # Change variable name                      
  mutate(voltot=vol1+vol2-volinter,                    
         vol.inter.rel=volinter/voltot*100,            
         minim=pmin(vol1,vol2),                          
         Turn=(2*minim-2*volinter)/(2*minim-volinter),
         Nest=(abs(vol1-vol2)/(vol1+vol2-volinter))*(volinter/(2*minim-volinter)))
mean(wB_TemporalYear_Vol$vol.inter.rel);se(wB_TemporalYear_Vol$vol.inter.rel) # mean and SE values of isotope volume # Change variable name   
mean(wB_TemporalYear_Vol$Turn);se(wB_TemporalYear_Vol$Turn) # Change variable name   
mean(wB_TemporalYear_Vol$Nest);se(wB_TemporalYear_Vol$Nest) # Change variable name   
wB_TemporalYear_Vol_Stats <- data.frame(wB_TemporalYear_Vol)
wB_TemporalYear_Vol_Stats$Species <- "wB"
# The following lines allow calculating the niche volumes of each population relative to the total niche volume
# occupied by all   from all populations (i.e., Tot.vol)
wB_TemporalYear_VolRel=data.frame(Site.Species=c(levels(comparison[drop=TRUE])),Volume=c(wB_TemporalYear_Vol[1,2],wB_TemporalYear_Vol[1:length(levels(comparison[drop=TRUE]))-1,3])) # Change variable
wB_TemporalYear_VolRel=wB_TemporalYear_VolRel%>% # Change variable
  mutate(vol.sp.rel=(round((Volume/wB_NicheSpaceTemporal_TotVol)*100,2))) # Change variable
mean(wB_TemporalYear_VolRel$vol.sp.rel);se(wB_TemporalYear_VolRel$vol.sp.rel) # Change variable
wB_TemporalYear_VolRel_data <- data.frame(wB_TemporalYear_VolRel)
wB_TemporalYear_VolRel_data$Species <- "wB"


# Combine all results
TemporalNiche_Vol_Stats <- rbind(bF_TemporalConHetero_Vol_Stats, bF_TemporalYear_Vol_Stats, wB_TemporalConHetero_Vol_Stats, wB_TemporalYear_Vol_Stats)
write.csv(TemporalNiche_Vol_Stats, "Tables/StoicNiche_TemporalVolStats.csv")
TemporalNiche_VolRel_Stats <- rbind(bF_TemporalConHetero_VolRel_data, bF_TemporalYear_VolRel_data, wB_TemporalConHetero_VolRel_data, wB_TemporalYear_VolRel_data)
write.csv(TemporalNiche_VolRel_Stats, "Tables/StoicNiche_TemporalVolRelStats.csv")


###### Temporal niche volume figures #######
# Create hypervolumne sphere figures using Gonzalez et. al., 2017 approach

# Use subsetted data from above - if starting from here; re-run load main dataset and data processing to get cons and hetero groups
bF_TemporalCNP
wB_TemporalCNP 

# Run stoichiometric niche functions - R code before running the code below
# Create 3D Spheres
# Balsam fir
bF_con2016_Matrix=as.matrix(bF_TemporalCNP[bF_TemporalCNP[,"Temporal"]=="Conspecific-2016",c("C","N","P")])
bF_con2016_Niche <- CHVind(bF_con2016_Matrix)

#bF_con2017_Matrix=as.matrix(bF_TemporalCNP[bF_TemporalCNP[,"Temporal"]=="Conspecific-2017",c("C","N","P")]) # less than 4 pt no compute
#bF_con2017_Niche <- CHVind(bF_con2017_Matrix)

bF_hetero2016_Matrix=as.matrix(bF_TemporalCNP[bF_TemporalCNP[,"Temporal"]=="Heterospecific-2016",c("C","N","P")])
bF_hetero2016_Niche <- CHVind(bF_hetero2016_Matrix)

bF_hetero2017_Matrix=as.matrix(bF_TemporalCNP[bF_TemporalCNP[,"Temporal"]=="Heterospecific-2017",c("C","N","P")])
bF_hetero2017_Niche <- CHVind(bF_hetero2017_Matrix)

bF_2016_Matrix=as.matrix(bF_TemporalCNP[bF_TemporalCNP[,"Year"]=="2016",c("C","N","P")])
bF_2016_Niche <- CHVind(bF_2016_Matrix)

bF_2017_Matrix=as.matrix(bF_TemporalCNP[wB_TemporalCNP[,"Year"]=="2017",c("C","N","P")])
bF_2017_Niche <- CHVind(bF_2017_Matrix)

# White birch
wB_con2016_Matrix=as.matrix(wB_TemporalCNP[wB_TemporalCNP[,"Temporal"]=="Conspecific-2016",c("C","N","P")])
wB_con2016_Niche <- CHVind(wB_con2016_Matrix)

wB_con2017_Matrix=as.matrix(wB_TemporalCNP[wB_TemporalCNP[,"Temporal"]=="Conspecific-2017",c("C","N","P")])
wB_con2017_Niche <- CHVind(wB_con2017_Matrix)

wB_hetero2016_Matrix=as.matrix(wB_TemporalCNP[wB_TemporalCNP[,"Temporal"]=="Heterospecific-2016",c("C","N","P")])
wB_hetero2016_Niche <- CHVind(wB_hetero2016_Matrix)

wB_hetero2017_Matrix=as.matrix(wB_TemporalCNP[wB_TemporalCNP[,"Temporal"]=="Heterospecific-2017",c("C","N","P")])
wB_hetero2017_Niche <- CHVind(wB_hetero2017_Matrix)

wB_2016_Matrix=as.matrix(wB_TemporalCNP[wB_TemporalCNP[,"Year"]=="2016",c("C","N","P")])
wB_2016_Niche <- CHVind(wB_2016_Matrix)

wB_2017_Matrix=as.matrix(wB_TemporalCNP[wB_TemporalCNP[,"Year"]=="2017",c("C","N","P")])
wB_2017_Niche <- CHVind(wB_2017_Matrix)


# Balsam fir and white birch regional 3D plot
# Extract coordinates for drop lines - numbers added below in plot
bF_con2016_Data <- as.data.frame(bF_con2016_Niche$vert_set1)
mean(bF_con2016_Data$V1)
mean(bF_con2016_Data$V2)
mean(bF_con2016_Data$V3)

bF_hetero2016_Data <- as.data.frame(bF_hetero2016_Niche$vert_set1)
mean(bF_hetero2016_Data$V1)
mean(bF_hetero2016_Data$V2)
mean(bF_hetero2016_Data$V3)

bF_hetero2017_Data <- as.data.frame(bF_hetero2017_Niche$vert_set1)
mean(bF_hetero2017_Data$V1)
mean(bF_hetero2017_Data$V2)
mean(bF_hetero2017_Data$V3)

bF_2016_Data <- as.data.frame(bF_2016_Niche$vert_set1)
mean(bF_2016_Data$V1)
mean(bF_2016_Data$V2)
mean(bF_2016_Data$V3)

bF_2017_Data <- as.data.frame(bF_2017_Niche$vert_set1)
mean(bF_2017_Data$V1)
mean(bF_2017_Data$V2)
mean(bF_2017_Data$V3)

# White birch
wB_con2016_Data <- as.data.frame(wB_con2016_Niche$vert_set1)
mean(wB_con2016_Data$V1)
mean(wB_con2016_Data$V2)
mean(wB_con2016_Data$V3)

wB_hetero2016_Data <- as.data.frame(wB_hetero2016_Niche$vert_set1)
mean(wB_hetero2016_Data$V1)
mean(wB_hetero2016_Data$V2)
mean(wB_hetero2016_Data$V3)

wB_con2017_Data <- as.data.frame(wB_con2017_Niche$vert_set1)
mean(wB_con2017_Data$V1)
mean(wB_con2017_Data$V2)
mean(wB_con2017_Data$V3)

wB_hetero2017_Data <- as.data.frame(wB_hetero2017_Niche$vert_set1)
mean(wB_hetero2017_Data$V1)
mean(wB_hetero2017_Data$V2)
mean(wB_hetero2017_Data$V3)

wB_2016_Data <- as.data.frame(wB_2016_Niche$vert_set1)
mean(wB_2016_Data$V1)
mean(wB_2016_Data$V2)
mean(wB_2016_Data$V3)

wB_2017_Data <- as.data.frame(wB_2017_Niche$vert_set1)
mean(wB_2017_Data$V1)
mean(wB_2017_Data$V2)
mean(wB_2017_Data$V3)


TemporalYearCol <- c("2016" = "darkgreen", "2017" = "darkblue") # Regional extent
TemporalCol <- c("Conspecific-2016" = "orangered1", "Heterospecific-2016" = "darkgoldenrod1", "Conspecific-2017" = "darkolivegreen4", "Heterospecific-2017" = "burlywood4") # Local extent

# Create plot
mfrow3d(1,2)
rgl.clear()
par3d(windowRect = c(400, 100, 1600, 700), cex=1.2)

# Balsam fir plot local
plot3d(bF_TemporalCNP[,5:7], size=2, type='p', axes=TRUE, edges="bbox", col=TemporalCol[bF_TemporalCNP$Temporal], repel=F, xlab=" ", ylab=" ", zlab=" ")
axes3d(edges="bbox") # Remove box face 
# Add title/axis labels
mtext3d(text="(a) Balsam fir", edge="x-+", line=3.2)
mtext3d(text="Carbon", edge="x-+", line=2.3)
mtext3d(text="Nitrogen", edge="y+-", line=2)
mtext3d(text="Phosphorus", edge="z+-", line=-2) 
# Add spheres
spheres3d(apply(bF_con2016_Niche$vert_set1,2,mean), alpha=0.7, radius = ((bF_con2016_Niche$vol[1]*3)/(4*pi))^(1/3),color="orangered1")
spheres3d(apply(bF_hetero2016_Niche$vert_set1,2,mean), alpha=0.7, radius = ((bF_hetero2016_Niche$vol[1]*3)/(4*pi))^(1/3),color="darkgoldenrod1")
spheres3d(apply(bF_hetero2017_Niche$vert_set1,2,mean), alpha=0.7, radius = ((bF_hetero2017_Niche$vol[1]*3)/(4*pi))^(1/3),color="burlywood4")
spheres3d(apply(bF_2016_Niche$vert_set1,2,mean), alpha=0.7, radius = ((bF_2016_Niche$vol[1]*3)/(4*pi))^(1/3),color="darkgreen")
spheres3d(apply(bF_2017_Niche$vert_set1,2,mean), alpha=0.7, radius = ((bF_2017_Niche$vol[1]*3)/(4*pi))^(1/3),color="darkblue")
# Balsam fir con 2016
abclines3d(52.23643,0.8764286, 0.08964286, a=diag(3), col="orangered1")
# Balsam fir hetero 2016
abclines3d(52.18765,0.9441176, 0.08452941, a=diag(3), col="darkgoldenrod1")
# Balsam fir hetero 2017
abclines3d(50.71176,1.034706, 0.09870588, a=diag(3), col="burlywood4")
# Balsam fir 2016
abclines3d(52.127,0.919, 0.0893, a=diag(3), col="darkgreen")
# Balsam fir 2017
abclines3d(50.80053,0.9442105, 0.08889474, a=diag(3), col="darkblue")

# White birch plot local
plot3d(wB_TemporalCNP[,5:7], size=2, type='p', axes=TRUE, edges="bbox", col=TemporalCol[wB_TemporalCNP$Temporal], repel=F, xlab=" ", ylab=" ", zlab=" ")
axes3d(edges="bbox") # Remove box face 
# Add title/axis labels
mtext3d(text="(b) White birch", edge="x-+", line=3.2)
mtext3d(text="Carbon", edge="x-+", line=2.3)
mtext3d(text="Nitrogen", edge="y+-", line=3)
mtext3d(text="Phosphorus", edge="z+-", line=-2) 
# Add spheres
spheres3d(apply(wB_con2016_Niche$vert_set1,2,mean), alpha=0.7, radius = ((wB_con2016_Niche$vol[1]*3)/(4*pi))^(1/3),color="orangered1")
spheres3d(apply(wB_hetero2016_Niche$vert_set1,2,mean), alpha=0.7, radius = ((wB_hetero2016_Niche$vol[1]*3)/(4*pi))^(1/3),color="darkgoldenrod1")
spheres3d(apply(wB_con2017_Niche$vert_set1,2,mean), alpha=0.7, radius = ((wB_con2017_Niche$vol[1]*3)/(4*pi))^(1/3),color="darkolivegreen4")
spheres3d(apply(wB_hetero2017_Niche$vert_set1,2,mean), alpha=0.7, radius = ((wB_hetero2017_Niche$vol[1]*3)/(4*pi))^(1/3),color="burlywood4")
spheres3d(apply(wB_2016_Niche$vert_set1,2,mean), alpha=0.7, radius = ((wB_2016_Niche$vol[1]*3)/(4*pi))^(1/3),color="darkgreen")
spheres3d(apply(wB_2017_Niche$vert_set1,2,mean), alpha=0.7, radius = ((wB_2017_Niche$vol[1]*3)/(4*pi))^(1/3),color="darkblue")
# White birch con 2016
abclines3d(49.84889,1.852222, 0.1875556, a=diag(3), col="orangered1")
# White birch hetero 2016
abclines3d(50.62267,1.555333, 0.1358667, a=diag(3), col="darkgoldenrod1")
# White birch con 2017
abclines3d(49.22,1.63, 0.1824, a=diag(3), col="darkolivegreen4")
# White birch hetero 2017
abclines3d(49.16154,1.726923, 0.1416154, a=diag(3), col="burlywood4")
# White birch 2016
abclines3d(50.22778,1.747222, 0.1692222, a=diag(3), col="darkgreen")
# White birch 2017
abclines3d(48.45455,1.578182, 0.1745455, a=diag(3), col="darkblue")


# Add legend
legend3d("topright", legend = c("2016", "2017", "Conspecific-2016", "Heterospecific-2016", "Conspecifics-2017", "Heterospecifics-2017"), col = c("darkgreen", "darkblue","orangered1","darkgoldenrod1", "darkolivegreen4", "burlywood4"), pch=16, inset=0.015, cex=1.3)

# Save figure - need to find a way to make/set resolution (good as is)
snapshot3d("Figures/bFwB_TemporalNiche.png",fmt="png", top=TRUE) #To save the plot3D




#### Test for sample size effect on niche volume #####
getwd()
FoliarCNP <- read.csv("Foliar_CNP.csv")
View(FoliarCNP)
# Clear all previous data
rm(list=ls())

# Subset data by species
# Subset by regional extent (i.e., species level)
FoliarCNP_bF <- subset(FoliarCNP, FoliarCNP$Species == "bF") # n=390
FoliarCNP_bF$RegComType <- as.factor(FoliarCNP_bF$RegComType)
FoliarCNP_wB <- subset(FoliarCNP, FoliarCNP$Species == "wB") # n=229
FoliarCNP_wB$RegComType <- as.factor(FoliarCNP_wB$RegComType)

# Calculate niche total available volumes
bF_NicheSpace <- as.matrix(FoliarCNP_bF[,c("C","N","P")])
bF_NicheSpace_TotVol <- CHVind(bF_NicheSpace)
bF_NicheSpace_TotVol <- as.numeric(bF_NicheSpace_TotVol$vol)

wB_NicheSpace <- as.matrix(FoliarCNP_wB[,c("C","N","P")])
wB_NicheSpace_TotVol <- CHVind(wB_NicheSpace)
wB_NicheSpace_TotVol <- as.numeric(wB_NicheSpace_TotVol$vol)

# Balsam fir subsets for processing
bF_Ecoreg_NP <- subset(FoliarCNP_bF, Region == "Northern Peninsula", select=c("C", "N", "P")) # n = 295
bF_Ecoreg_CF <- subset(FoliarCNP_bF, Region == "Central Forest", select=c("C", "N", "P")) # n = 95
bF_Reg_Cons <- subset(FoliarCNP_bF, RegComType == "Conspecific", select=c("C", "N", "P")) # n = 189
bF_Reg_Hetero <- subset(FoliarCNP_bF, RegComType == "Heterospecific", select=c("C", "N", "P")) # n = 201
bF_LocalCF_Cons <- subset(FoliarCNP_bF, LocalComType == "CF: Conspecific", select=c("C", "N", "P")) # n= 47
bF_LocalCF_Hetero <- subset(FoliarCNP_bF, LocalComType == "CF: Heterospecific", select=c("C", "N", "P")) # n = 48
bF_LocalNP_Cons <- subset(FoliarCNP_bF, LocalComType == "NP: Conspecific", select=c("C", "N", "P")) # n = 142
bF_LocalNP_Hetero <- subset(FoliarCNP_bF, LocalComType == "NP: Heterospecific", select=c("C", "N", "P")) # n = 153

# White birch subset for processing
wB_Ecoreg_NP <- subset(FoliarCNP_wB, Region == "Northern Peninsula", select=c("C", "N", "P")) # n = 158
wB_Ecoreg_CF <- subset(FoliarCNP_wB, Region == "Central Forest", select=c("C", "N", "P")) # n = 71
wB_Reg_Cons <- subset(FoliarCNP_wB, RegComType == "Conspecific", select=c("C", "N", "P")) # n = 28
wB_Reg_Hetero <- subset(FoliarCNP_wB, RegComType == "Heterospecific", select=c("C", "N", "P")) # n = 201
wB_LocalCF_Cons <- subset(FoliarCNP_wB, LocalComType == "CF: Conspecific", select=c("C", "N", "P")) # n = 23
wB_LocalCF_Hetero <- subset(FoliarCNP_wB, LocalComType == "CF: Heterospecific", select=c("C", "N", "P")) # n = 48
wB_LocalNP_Cons <- subset(FoliarCNP_wB, LocalComType == "NP: Conspecific", select=c("C", "N", "P")) # n = 5
wB_LocalNP_Hetero <- subset(FoliarCNP_wB, LocalComType == "NP: Heterospecific", select=c("C", "N", "P")) # n = 153
library(rcdd)
library(geometry)
library(dplyr)
##### Generate random niche volumes #####
#####  bF ecoregion - North Peninsula ##### 
bF_Ecoreg_NP_Rand_pas01=rep(seq(5,295,5),each=999) # n = 295
bF_Ecoreg_NP_Rand<-matrix(NA,ncol=1,nrow=length(bF_Ecoreg_NP_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(bF_Ecoreg_NP_Rand_pas01)){
  bF_Ecoreg_NP_Rand[i,1]<-round(CHVind(as.matrix(bF_Ecoreg_NP[sample(seq(1,length(bF_Ecoreg_NP[,1]),1),bF_Ecoreg_NP_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(bF_Ecoreg_NP_Rand,"D:/RandomNiches/bF_Ecoreg_NP_Rand_v1.csv")
bF_Ecoreg_NP_Rand <- read.csv("RandomNiches/bF_Ecoreg_NP_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
bF_Ecoreg_NP_Rand_pas01=rep(seq(5,295,5),each=999)
bF_Ecoreg_NP_Rand_SP=cbind(bF_Ecoreg_NP_Rand,bF_Ecoreg_NP_Rand_pas01)
bF_Ecoreg_NP_Rand_SP$bF_Ecoreg_NP_Rand_pas01=as.factor(bF_Ecoreg_NP_Rand_SP$bF_Ecoreg_NP_Rand_pas01)
bF_Ecoreg_NP_Rand_SP1=bF_Ecoreg_NP_Rand_SP%>%
  group_by(bF_Ecoreg_NP_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
bF_Ecoreg_NP_Rand_SP1=bF_Ecoreg_NP_Rand_SP1[!duplicated(bF_Ecoreg_NP_Rand_SP1$sd.vol), ]
bF_Ecoreg_NP_Rand_pas.vol=seq(5,295,5)
bF_Ecoreg_NP_Rand_SP1=data.frame(bF_Ecoreg_NP_Rand_SP1,bF_Ecoreg_NP_Rand_pas.vol)

bF_Ecoreg_NP_Matrix <- as.matrix(bF_Ecoreg_NP[,c("C", "N", "P")])
bF_Ecoreg_NP_Obs<-CHVind(bF_Ecoreg_NP_Matrix)$vol


plot(bF_Ecoreg_NP_Rand_SP1$mean.vol~bF_Ecoreg_NP_Rand_SP1$bF_Ecoreg_NP_Rand_pas.vol,pch=16,ylim=c(0,0.6),
     main = "bF Northern Peninsula ecoregion\\n n = 295 (5,295,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_Ecoreg_NP_Rand_SP1$bF_Ecoreg_NP_Rand_pas.vol),c(bF_Ecoreg_NP_Rand_SP1$mean.vol-2*bF_Ecoreg_NP_Rand_SP1$sd.vol),c(bF_Ecoreg_NP_Rand_SP1$bF_Ecoreg_NP_Rand_pas.vol),c(bF_Ecoreg_NP_Rand_SP1$mean.vol+2*bF_Ecoreg_NP_Rand_SP1$sd.vol))
abline(h=bF_Ecoreg_NP_Obs)


###### bF Ecoregion - Central Forest #####
bF_Ecoreg_CF_Rand_pas01=rep(seq(5,95,5),each=999) # n = 95
bF_Ecoreg_CF_Rand<-matrix(NA,ncol=1,nrow=length(bF_Ecoreg_CF_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(bF_Ecoreg_CF_Rand_pas01)){
  bF_Ecoreg_CF_Rand[i,1]<-round(CHVind(as.matrix(bF_Ecoreg_CF[sample(seq(1,length(bF_Ecoreg_CF[,1]),1),bF_Ecoreg_CF_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(bF_Ecoreg_CF_Rand,"D:/RandomNiches/bF_Ecoreg_CF_Rand.csv")
bF_Ecoreg_CF_Rand <- read.csv("RandomNiches/bF_Ecoreg_CF_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
bF_Ecoreg_CF_Rand_pas01=rep(seq(5,95,5),each=999)
bF_Ecoreg_CF_Rand_SP=cbind(bF_Ecoreg_CF_Rand,bF_Ecoreg_CF_Rand_pas01)
bF_Ecoreg_CF_Rand_SP$bF_Ecoreg_CF_Rand_pas01=as.factor(bF_Ecoreg_CF_Rand_SP$bF_Ecoreg_CF_Rand_pas01)
bF_Ecoreg_CF_Rand_SP1=bF_Ecoreg_CF_Rand_SP%>%
  group_by(bF_Ecoreg_CF_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
bF_Ecoreg_CF_Rand_SP1=bF_Ecoreg_CF_Rand_SP1[!duplicated(bF_Ecoreg_CF_Rand_SP1$sd.vol), ]
bF_Ecoreg_CF_Rand_pas.vol=seq(5,95,5)
bF_Ecoreg_CF_Rand_SP1=data.frame(bF_Ecoreg_CF_Rand_SP1,bF_Ecoreg_CF_Rand_pas.vol)

bF_Ecoreg_CF_Matrix <- as.matrix(bF_Ecoreg_CF[,c("C", "N", "P")])
bF_Ecoreg_CF_Obs<-CHVind(bF_Ecoreg_CF_Matrix)$vol

plot(bF_Ecoreg_CF_Rand_SP1$mean.vol~bF_Ecoreg_CF_Rand_SP1$bF_Ecoreg_CF_Rand_pas.vol,pch=16,ylim=c(0,0.1),
     main = "bF Central Forest ecoregion\\n n = 95 (5,95,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_Ecoreg_CF_Rand_SP1$bF_Ecoreg_CF_Rand_pas.vol),c(bF_Ecoreg_CF_Rand_SP1$mean.vol-2*bF_Ecoreg_CF_Rand_SP1$sd.vol),c(bF_Ecoreg_CF_Rand_SP1$bF_Ecoreg_CF_Rand_pas.vol),c(bF_Ecoreg_CF_Rand_SP1$mean.vol+2*bF_Ecoreg_CF_Rand_SP1$sd.vol))
abline(h=bF_Ecoreg_CF_Obs)


##### bF Regional Conspecific ##### 
bF_Reg_Cons_Rand_pas01=rep(seq(7,189,7),each=999) # n = 189
bF_Reg_Cons_Rand<-matrix(NA,ncol=1,nrow=length(bF_Reg_Cons_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(bF_Reg_Cons_Rand_pas01)){
  bF_Reg_Cons_Rand[i,1]<-round(CHVind(as.matrix(bF_Reg_Cons[sample(seq(1,length(bF_Reg_Cons[,1]),1),bF_Reg_Cons_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(bF_Reg_Cons_Rand,"D:/RandomNiches/bF_Reg_Cons_Rand.csv")
bF_Reg_Cons_Rand <- read.csv("RandomNiches/bF_Reg_Cons_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
bF_Reg_Cons_Rand_pas01=rep(seq(7,189,7),each=999)
bF_Reg_Cons_Rand_SP=cbind(bF_Reg_Cons_Rand,bF_Reg_Cons_Rand_pas01)
bF_Reg_Cons_Rand_SP$bF_Reg_Cons_Rand_pas01=as.factor(bF_Reg_Cons_Rand_SP$bF_Reg_Cons_Rand_pas01)
bF_Reg_Cons_Rand_SP1=bF_Reg_Cons_Rand_SP%>%
  group_by(bF_Reg_Cons_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
bF_Reg_Cons_Rand_SP1=bF_Reg_Cons_Rand_SP1[!duplicated(bF_Reg_Cons_Rand_SP1$sd.vol), ]
bF_Reg_Cons_Rand_pas.vol=seq(7,189,7)
bF_Reg_Cons_Rand_SP1=data.frame(bF_Reg_Cons_Rand_SP1,bF_Reg_Cons_Rand_pas.vol)

bF_Reg_Cons_Matrix <- as.matrix(bF_Reg_Cons[,c("C", "N", "P")])
bF_Reg_Cons_Obs<-CHVind(bF_Reg_Cons_Matrix)$vol

plot(bF_Reg_Cons_Rand_SP1$mean.vol~bF_Reg_Cons_Rand_SP1$bF_Reg_Cons_Rand_pas.vol,pch=16,ylim=c(0,0.4),
     main = "bF Regional Conspecific\\n n = 189 (7,189,7)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_Reg_Cons_Rand_SP1$bF_Reg_Cons_Rand_pas.vol),c(bF_Reg_Cons_Rand_SP1$mean.vol-2*bF_Reg_Cons_Rand_SP1$sd.vol),c(bF_Reg_Cons_Rand_SP1$bF_Reg_Cons_Rand_pas.vol),c(bF_Reg_Cons_Rand_SP1$mean.vol+2*bF_Reg_Cons_Rand_SP1$sd.vol))
abline(h=bF_Reg_Cons_Obs)


##### bF Regional Heterospecific ##### 
bF_Reg_Hetero_Rand_pas01=rep(seq(5,200,5),each=999) # n = 201
bF_Reg_Hetero_Rand<-matrix(NA,ncol=1,nrow=length(bF_Reg_Hetero_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(bF_Reg_Hetero_Rand_pas01)){
  bF_Reg_Hetero_Rand[i,1]<-round(CHVind(as.matrix(bF_Reg_Hetero[sample(seq(1,length(bF_Reg_Hetero[,1]),1),bF_Reg_Hetero_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(bF_Reg_Hetero_Rand,"D:/RandomNiches/bF_Reg_Hetero_Rand_v1.csv")
bF_Reg_Hetero_Rand <- read.csv("RandomNiches/bF_Reg_Hetero_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
bF_Reg_Hetero_Rand_pas01=rep(seq(5,200,5),each=999)
bF_Reg_Hetero_Rand_SP=cbind(bF_Reg_Hetero_Rand,bF_Reg_Hetero_Rand_pas01)
bF_Reg_Hetero_Rand_SP$bF_Reg_Hetero_Rand_pas01=as.factor(bF_Reg_Hetero_Rand_SP$bF_Reg_Hetero_Rand_pas01)
bF_Reg_Hetero_Rand_SP1=bF_Reg_Hetero_Rand_SP%>%
  group_by(bF_Reg_Hetero_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
bF_Reg_Hetero_Rand_SP1=bF_Reg_Hetero_Rand_SP1[!duplicated(bF_Reg_Hetero_Rand_SP1$sd.vol), ]
bF_Reg_Hetero_Rand_pas.vol=seq(5,200,5)
bF_Reg_Hetero_Rand_SP1=data.frame(bF_Reg_Hetero_Rand_SP1,bF_Reg_Hetero_Rand_pas.vol)

bF_Reg_Hetero_Matrix <- as.matrix(bF_Reg_Hetero[,c("C", "N", "P")])
bF_Reg_Hetero_Obs<-CHVind(bF_Reg_Hetero_Matrix)$vol

plot(bF_Reg_Hetero_Rand_SP1$mean.vol~bF_Reg_Hetero_Rand_SP1$bF_Reg_Hetero_Rand_pas.vol,pch=16,ylim=c(0,0.6),
    main = "bF Regional Heterospecific\\n n = 201 (5,200,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_Reg_Hetero_Rand_SP1$bF_Reg_Hetero_Rand_pas.vol),c(bF_Reg_Hetero_Rand_SP1$mean.vol-2*bF_Reg_Hetero_Rand_SP1$sd.vol),c(bF_Reg_Hetero_Rand_SP1$bF_Reg_Hetero_Rand_pas.vol),c(bF_Reg_Hetero_Rand_SP1$mean.vol+2*bF_Reg_Hetero_Rand_SP1$sd.vol))
abline(h=bF_Reg_Hetero_Obs)


##### bF Local CF Conspecific #####
bF_LocalCF_Cons_Rand_pas01=rep(seq(5,45,5),each=999) # n = 47
bF_LocalCF_Cons_Rand<-matrix(NA,ncol=1,nrow=length(bF_LocalCF_Cons_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(bF_LocalCF_Cons_Rand_pas01)){
  bF_LocalCF_Cons_Rand[i,1]<-round(CHVind(as.matrix(bF_LocalCF_Cons[sample(seq(1,length(bF_LocalCF_Cons[,1]),1),bF_LocalCF_Cons_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(bF_LocalCF_Cons_Rand,"D:/RandomNiches/bF_LocalCF_Cons_Rand.csv")
bF_LocalCF_Cons_Rand <- read.csv("RandomNiches/bF_LocalCF_Cons_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
bF_LocalCF_Cons_Rand_pas01=rep(seq(5,45,5),each=999)
bF_LocalCF_Cons_Rand_SP=cbind(bF_LocalCF_Cons_Rand,bF_LocalCF_Cons_Rand_pas01)
bF_LocalCF_Cons_Rand_SP$bF_LocalCF_Cons_Rand_pas01=as.factor(bF_LocalCF_Cons_Rand_SP$bF_LocalCF_Cons_Rand_pas01)
bF_LocalCF_Cons_Rand_SP1=bF_LocalCF_Cons_Rand_SP%>%
  group_by(bF_LocalCF_Cons_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
bF_LocalCF_Cons_Rand_SP1=bF_LocalCF_Cons_Rand_SP1[!duplicated(bF_LocalCF_Cons_Rand_SP1$sd.vol), ]
bF_LocalCF_Cons_Rand_pas.vol=seq(5,45,5)
bF_LocalCF_Cons_Rand_SP1=data.frame(bF_LocalCF_Cons_Rand_SP1,bF_LocalCF_Cons_Rand_pas.vol)

bF_LocalCF_Cons_Matrix <- as.matrix(bF_LocalCF_Cons[,c("C", "N", "P")])
bF_LocalCF_Cons_Obs<-CHVind(bF_LocalCF_Cons_Matrix)$vol

plot(bF_LocalCF_Cons_Rand_SP1$mean.vol~bF_LocalCF_Cons_Rand_SP1$bF_LocalCF_Cons_Rand_pas.vol,pch=16,ylim=c(0,0.07),
     main = "bF Local CF Conspecific\\n n = 47 (5,45,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_LocalCF_Cons_Rand_SP1$bF_LocalCF_Cons_Rand_pas.vol),c(bF_LocalCF_Cons_Rand_SP1$mean.vol-2*bF_LocalCF_Cons_Rand_SP1$sd.vol),c(bF_LocalCF_Cons_Rand_SP1$bF_LocalCF_Cons_Rand_pas.vol),c(bF_LocalCF_Cons_Rand_SP1$mean.vol+2*bF_LocalCF_Cons_Rand_SP1$sd.vol))
abline(h=bF_LocalCF_Cons_Obs)


##### bF Local CF Heterospecific ##### 
bF_LocalCF_Hetero_Rand_pas01=rep(seq(4,48,4),each=999) # n = 48
bF_LocalCF_Hetero_Rand<-matrix(NA,ncol=1,nrow=length(bF_LocalCF_Hetero_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(bF_LocalCF_Hetero_Rand_pas01)){
  bF_LocalCF_Hetero_Rand[i,1]<-round(CHVind(as.matrix(bF_LocalCF_Hetero[sample(seq(1,length(bF_LocalCF_Hetero[,1]),1),bF_LocalCF_Hetero_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(bF_LocalCF_Hetero_Rand,"D:/RandomNiches/bF_LocalCF_Hetero_Rand.csv")
bF_LocalCF_Hetero_Rand <- read.csv("RandomNiches/bF_LocalCF_Hetero_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
bF_LocalCF_Hetero_Rand_pas01=rep(seq(4,48,4),each=999)
bF_LocalCF_Hetero_Rand_SP=cbind(bF_LocalCF_Hetero_Rand,bF_LocalCF_Hetero_Rand_pas01)
bF_LocalCF_Hetero_Rand_SP$bF_LocalCF_Hetero_Rand_pas01=as.factor(bF_LocalCF_Hetero_Rand_SP$bF_LocalCF_Hetero_Rand_pas01)
bF_LocalCF_Hetero_Rand_SP1=bF_LocalCF_Hetero_Rand_SP%>%
  group_by(bF_LocalCF_Hetero_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
bF_LocalCF_Hetero_Rand_SP1=bF_LocalCF_Hetero_Rand_SP1[!duplicated(bF_LocalCF_Hetero_Rand_SP1$sd.vol), ]
bF_LocalCF_Hetero_Rand_pas.vol=seq(4,48,4)
bF_LocalCF_Hetero_Rand_SP1=data.frame(bF_LocalCF_Hetero_Rand_SP1,bF_LocalCF_Hetero_Rand_pas.vol)

bF_LocalCF_Hetero_Matrix <- as.matrix(bF_LocalCF_Hetero[,c("C", "N", "P")])
bF_LocalCF_Hetero_Obs<-CHVind(bF_LocalCF_Hetero_Matrix)$vol

plot(bF_LocalCF_Hetero_Rand_SP1$mean.vol~bF_LocalCF_Hetero_Rand_SP1$bF_LocalCF_Hetero_Rand_pas.vol,pch=16,ylim=c(0,0.07),
     main = "bF Local CF Heterospecific\\n n = 48 (4,48,4)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_LocalCF_Hetero_Rand_SP1$bF_LocalCF_Hetero_Rand_pas.vol),c(bF_LocalCF_Hetero_Rand_SP1$mean.vol-2*bF_LocalCF_Hetero_Rand_SP1$sd.vol),c(bF_LocalCF_Hetero_Rand_SP1$bF_LocalCF_Hetero_Rand_pas.vol),c(bF_LocalCF_Hetero_Rand_SP1$mean.vol+2*bF_LocalCF_Hetero_Rand_SP1$sd.vol))
abline(h=bF_LocalCF_Hetero_Obs)



##### bF Local NP Conspecific ##### 
bF_LocalNP_Cons_Rand_pas01=rep(seq(4,140,4),each=999) # n = 142
bF_LocalNP_Cons_Rand<-matrix(NA,ncol=1,nrow=length(bF_LocalNP_Cons_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(bF_LocalNP_Cons_Rand_pas01)){
  bF_LocalNP_Cons_Rand[i,1]<-round(CHVind(as.matrix(bF_LocalNP_Cons[sample(seq(1,length(bF_LocalNP_Cons[,1]),1),bF_LocalNP_Cons_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(bF_LocalNP_Cons_Rand,"D:/RandomNiches/bF_LocalNP_Cons_Rand.csv")
bF_LocalNP_Cons_Rand <- read.csv("RandomNiches/bF_LocalNP_Cons_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
bF_LocalNP_Cons_Rand_pas01=rep(seq(4,140,4),each=999)
bF_LocalNP_Cons_Rand_SP=cbind(bF_LocalNP_Cons_Rand,bF_LocalNP_Cons_Rand_pas01)
bF_LocalNP_Cons_Rand_SP$bF_LocalNP_Cons_Rand_pas01=as.factor(bF_LocalNP_Cons_Rand_SP$bF_LocalNP_Cons_Rand_pas01)
bF_LocalNP_Cons_Rand_SP1=bF_LocalNP_Cons_Rand_SP%>%
  group_by(bF_LocalNP_Cons_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
bF_LocalNP_Cons_Rand_SP1=bF_LocalNP_Cons_Rand_SP1[!duplicated(bF_LocalNP_Cons_Rand_SP1$sd.vol), ]
bF_LocalNP_Cons_Rand_pas.vol=seq(4,140,4)
bF_LocalNP_Cons_Rand_SP1=data.frame(bF_LocalNP_Cons_Rand_SP1,bF_LocalNP_Cons_Rand_pas.vol)

bF_LocalNP_Cons_Matrix <- as.matrix(bF_LocalNP_Cons[,c("C", "N", "P")])
bF_LocalNP_Cons_Obs<-CHVind(bF_LocalNP_Cons_Matrix)$vol

plot(bF_LocalNP_Cons_Rand_SP1$mean.vol~bF_LocalNP_Cons_Rand_SP1$bF_LocalNP_Cons_Rand_pas.vol,pch=16,ylim=c(0,0.3),
     main = "bF Local NP Conspecific\\n n = 142 (4,140,4)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_LocalNP_Cons_Rand_SP1$bF_LocalNP_Cons_Rand_pas.vol),c(bF_LocalNP_Cons_Rand_SP1$mean.vol-2*bF_LocalNP_Cons_Rand_SP1$sd.vol),c(bF_LocalNP_Cons_Rand_SP1$bF_LocalNP_Cons_Rand_pas.vol),c(bF_LocalNP_Cons_Rand_SP1$mean.vol+2*bF_LocalNP_Cons_Rand_SP1$sd.vol))
abline(h=bF_LocalNP_Cons_Obs)


##### bF Local NP Heterospecific ##### 
bF_LocalNP_Hetero_Rand_pas01=rep(seq(9,153,9),each=999) # n = 153
bF_LocalNP_Hetero_Rand<-matrix(NA,ncol=1,nrow=length(bF_LocalNP_Hetero_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(bF_LocalNP_Hetero_Rand_pas01)){
  bF_LocalNP_Hetero_Rand[i,1]<-round(CHVind(as.matrix(bF_LocalNP_Hetero[sample(seq(1,length(bF_LocalNP_Hetero[,1]),1),bF_LocalNP_Hetero_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(bF_LocalNP_Hetero_Rand,"D:/RandomNiches/bF_LocalNP_Hetero_Rand.csv")
bF_LocalNP_Hetero_Rand <- read.csv("RandomNiches/bF_LocalNP_Hetero_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
bF_LocalNP_Hetero_Rand_pas01=rep(seq(9,153,9),each=999)
bF_LocalNP_Hetero_Rand_SP=cbind(bF_LocalNP_Hetero_Rand,bF_LocalNP_Hetero_Rand_pas01)
bF_LocalNP_Hetero_Rand_SP$bF_LocalNP_Hetero_Rand_pas01=as.factor(bF_LocalNP_Hetero_Rand_SP$bF_LocalNP_Hetero_Rand_pas01)
bF_LocalNP_Hetero_Rand_SP1=bF_LocalNP_Hetero_Rand_SP%>%
  group_by(bF_LocalNP_Hetero_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
bF_LocalNP_Hetero_Rand_SP1=bF_LocalNP_Hetero_Rand_SP1[!duplicated(bF_LocalNP_Hetero_Rand_SP1$sd.vol), ]
bF_LocalNP_Hetero_Rand_pas.vol=seq(9,153,9)
bF_LocalNP_Hetero_Rand_SP1=data.frame(bF_LocalNP_Hetero_Rand_SP1,bF_LocalNP_Hetero_Rand_pas.vol)

bF_LocalNP_Hetero_Matrix <- as.matrix(bF_LocalNP_Hetero[,c("C", "N", "P")])
bF_LocalNP_Hetero_Obs<-CHVind(bF_LocalNP_Hetero_Matrix)$vol

plot(bF_LocalNP_Hetero_Rand_SP1$mean.vol~bF_LocalNP_Hetero_Rand_SP1$bF_LocalNP_Hetero_Rand_pas.vol,pch=16,ylim=c(0,0.5),
     main = "bF Local NP Heterospecific\\n n = 153 (9,153,9)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_LocalNP_Hetero_Rand_SP1$bF_LocalNP_Hetero_Rand_pas.vol),c(bF_LocalNP_Hetero_Rand_SP1$mean.vol-2*bF_LocalNP_Hetero_Rand_SP1$sd.vol),c(bF_LocalNP_Hetero_Rand_SP1$bF_LocalNP_Hetero_Rand_pas.vol),c(bF_LocalNP_Hetero_Rand_SP1$mean.vol+2*bF_LocalNP_Hetero_Rand_SP1$sd.vol))
abline(h=bF_LocalNP_Hetero_Obs)

##### wB Ecoregion - North Peninsula ##### 
wB_Ecoreg_NP_Rand_pas01=rep(seq(4,156,4),each=999) # n = 158
wB_Ecoreg_NP_Rand<-matrix(NA,ncol=1,nrow=length(wB_Ecoreg_NP_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(wB_Ecoreg_NP_Rand_pas01)){
  wB_Ecoreg_NP_Rand[i,1]<-round(CHVind(as.matrix(wB_Ecoreg_NP[sample(seq(1,length(wB_Ecoreg_NP[,1]),1),wB_Ecoreg_NP_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(wB_Ecoreg_NP_Rand,"D:/RandomNiches/wB_Ecoreg_NP_Rand.csv")
wB_Ecoreg_NP_Rand <- read.csv("RandomNiches/wB_Ecoreg_NP_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
wB_Ecoreg_NP_Rand_pas01=rep(seq(4,156,4),each=999)
wB_Ecoreg_NP_Rand_SP=cbind(wB_Ecoreg_NP_Rand,wB_Ecoreg_NP_Rand_pas01)
wB_Ecoreg_NP_Rand_SP$wB_Ecoreg_NP_Rand_pas01=as.factor(wB_Ecoreg_NP_Rand_SP$wB_Ecoreg_NP_Rand_pas01)
wB_Ecoreg_NP_Rand_SP1=wB_Ecoreg_NP_Rand_SP%>%
  group_by(wB_Ecoreg_NP_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
wB_Ecoreg_NP_Rand_SP1=wB_Ecoreg_NP_Rand_SP1[!duplicated(wB_Ecoreg_NP_Rand_SP1$sd.vol), ]
wB_Ecoreg_NP_Rand_pas.vol=seq(4,156,4)
wB_Ecoreg_NP_Rand_SP1=data.frame(wB_Ecoreg_NP_Rand_SP1,wB_Ecoreg_NP_Rand_pas.vol)

wB_Ecoreg_NP_Matrix <- as.matrix(wB_Ecoreg_NP[,c("C", "N", "P")])
wB_Ecoreg_NP_Obs<-CHVind(wB_Ecoreg_NP_Matrix)$vol

plot(wB_Ecoreg_NP_Rand_SP1$mean.vol~wB_Ecoreg_NP_Rand_SP1$wB_Ecoreg_NP_Rand_pas.vol,pch=16,ylim=c(0,2.5),
     main = "wB Northern Peninsula ecoregion\\n n = 158 (4,156,4)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_Ecoreg_NP_Rand_SP1$wB_Ecoreg_NP_Rand_pas.vol),c(wB_Ecoreg_NP_Rand_SP1$mean.vol-2*wB_Ecoreg_NP_Rand_SP1$sd.vol),c(wB_Ecoreg_NP_Rand_SP1$wB_Ecoreg_NP_Rand_pas.vol),c(wB_Ecoreg_NP_Rand_SP1$mean.vol+2*wB_Ecoreg_NP_Rand_SP1$sd.vol))
abline(h=wB_Ecoreg_NP_Obs)



###### wB Ecoregion - Central Forest #####
wB_Ecoreg_CF_Rand_pas01=rep(seq(5,70,5),each=999) # n = 71
wB_Ecoreg_CF_Rand<-matrix(NA,ncol=1,nrow=length(wB_Ecoreg_CF_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(wB_Ecoreg_CF_Rand_pas01)){
  wB_Ecoreg_CF_Rand[i,1]<-round(CHVind(as.matrix(wB_Ecoreg_CF[sample(seq(1,length(wB_Ecoreg_CF[,1]),1),wB_Ecoreg_CF_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(wB_Ecoreg_CF_Rand,"D:/RandomNiches/wB_Ecoreg_CF_Rand.csv")
wB_Ecoreg_CF_Rand <- read.csv("RandomNiches/wB_Ecoreg_CF_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
wB_Ecoreg_CF_Rand_pas01=rep(seq(5,70,5),each=999)
wB_Ecoreg_CF_Rand_SP=cbind(wB_Ecoreg_CF_Rand,wB_Ecoreg_CF_Rand_pas01)
wB_Ecoreg_CF_Rand_SP$wB_Ecoreg_CF_Rand_pas01=as.factor(wB_Ecoreg_CF_Rand_SP$wB_Ecoreg_CF_Rand_pas01)
wB_Ecoreg_CF_Rand_SP1=wB_Ecoreg_CF_Rand_SP%>%
  group_by(wB_Ecoreg_CF_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
wB_Ecoreg_CF_Rand_SP1=wB_Ecoreg_CF_Rand_SP1[!duplicated(wB_Ecoreg_CF_Rand_SP1$sd.vol), ]
wB_Ecoreg_CF_Rand_pas.vol=seq(5,70,5)
wB_Ecoreg_CF_Rand_SP1=data.frame(wB_Ecoreg_CF_Rand_SP1,wB_Ecoreg_CF_Rand_pas.vol)

wB_Ecoreg_CF_Matrix <- as.matrix(wB_Ecoreg_CF[,c("C", "N", "P")])
wB_Ecoreg_CF_Obs<-CHVind(wB_Ecoreg_CF_Matrix)$vol

plot(wB_Ecoreg_CF_Rand_SP1$mean.vol~wB_Ecoreg_CF_Rand_SP1$wB_Ecoreg_CF_Rand_pas.vol,pch=16,ylim=c(0,0.8),
     main = "wB Central Forest ecoregion\\n n = 71 (5,70,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_Ecoreg_CF_Rand_SP1$wB_Ecoreg_CF_Rand_pas.vol),c(wB_Ecoreg_CF_Rand_SP1$mean.vol-2*wB_Ecoreg_CF_Rand_SP1$sd.vol),c(wB_Ecoreg_CF_Rand_SP1$wB_Ecoreg_CF_Rand_pas.vol),c(wB_Ecoreg_CF_Rand_SP1$mean.vol+2*wB_Ecoreg_CF_Rand_SP1$sd.vol))
abline(h=wB_Ecoreg_CF_Obs)


##### wB Regional Conspecific #####
wB_Reg_Cons_Rand_pas01=rep(seq(4,28,4),each=999) # n = 28
wB_Reg_Cons_Rand<-matrix(NA,ncol=1,nrow=length(wB_Reg_Cons_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(wB_Reg_Cons_Rand_pas01)){
  wB_Reg_Cons_Rand[i,1]<-round(CHVind(as.matrix(wB_Reg_Cons[sample(seq(1,length(wB_Reg_Cons[,1]),1),wB_Reg_Cons_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(wB_Reg_Cons_Rand,"D:/RandomNiches/wB_Reg_Cons_Rand.csv")
wB_Reg_Cons_Rand <- read.csv("RandomNiches/wB_Reg_Cons_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
wB_Reg_Cons_Rand_pas01=rep(seq(4,28,4),each=999)
wB_Reg_Cons_Rand_SP=cbind(wB_Reg_Cons_Rand,wB_Reg_Cons_Rand_pas01)
wB_Reg_Cons_Rand_SP$wB_Reg_Cons_Rand_pas01=as.factor(wB_Reg_Cons_Rand_SP$wB_Reg_Cons_Rand_pas01)
wB_Reg_Cons_Rand_SP1=wB_Reg_Cons_Rand_SP%>%
  group_by(wB_Reg_Cons_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
wB_Reg_Cons_Rand_SP1=wB_Reg_Cons_Rand_SP1[!duplicated(wB_Reg_Cons_Rand_SP1$sd.vol), ]
wB_Reg_Cons_Rand_pas.vol=seq(4,28,4)
wB_Reg_Cons_Rand_SP1=data.frame(wB_Reg_Cons_Rand_SP1,wB_Reg_Cons_Rand_pas.vol)

wB_Reg_Cons_Matrix <- as.matrix(wB_Reg_Cons[,c("C", "N", "P")])
wB_Reg_Cons_Obs<-CHVind(wB_Reg_Cons_Matrix)$vol

plot(wB_Reg_Cons_Rand_SP1$mean.vol~wB_Reg_Cons_Rand_SP1$wB_Reg_Cons_Rand_pas.vol,pch=16,ylim=c(0,1),
     main = "wB Regional Conspecific\\n n = 28 (4,28,1)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_Reg_Cons_Rand_SP1$wB_Reg_Cons_Rand_pas.vol),c(wB_Reg_Cons_Rand_SP1$mean.vol-2*wB_Reg_Cons_Rand_SP1$sd.vol),c(wB_Reg_Cons_Rand_SP1$wB_Reg_Cons_Rand_pas.vol),c(wB_Reg_Cons_Rand_SP1$mean.vol+2*wB_Reg_Cons_Rand_SP1$sd.vol))
abline(h=wB_Reg_Cons_Obs)


##### wB Regional Heterospecific ##### 
wB_Reg_Hetero_Rand_pas01=rep(seq(5,200,5),each=999) # n = 201
wB_Reg_Hetero_Rand<-matrix(NA,ncol=1,nrow=length(wB_Reg_Hetero_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(wB_Reg_Hetero_Rand_pas01)){
  wB_Reg_Hetero_Rand[i,1]<-round(CHVind(as.matrix(wB_Reg_Hetero[sample(seq(1,length(wB_Reg_Hetero[,1]),1),wB_Reg_Hetero_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(wB_Reg_Hetero_Rand,"D:/RandomNiches/wB_Reg_Hetero_Rand.csv")
wB_Reg_Hetero_Rand <- read.csv("RandomNiches/wB_Reg_Hetero_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
wB_Reg_Hetero_Rand_pas01=rep(seq(5,200,5),each=999)
wB_Reg_Hetero_Rand_SP=cbind(wB_Reg_Hetero_Rand,wB_Reg_Hetero_Rand_pas01)
wB_Reg_Hetero_Rand_SP$wB_Reg_Hetero_Rand_pas01=as.factor(wB_Reg_Hetero_Rand_SP$wB_Reg_Hetero_Rand_pas01)
wB_Reg_Hetero_Rand_SP1=wB_Reg_Hetero_Rand_SP%>%
  group_by(wB_Reg_Hetero_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
wB_Reg_Hetero_Rand_SP1=wB_Reg_Hetero_Rand_SP1[!duplicated(wB_Reg_Hetero_Rand_SP1$sd.vol), ]
wB_Reg_Hetero_Rand_pas.vol=seq(5,200,5)
wB_Reg_Hetero_Rand_SP1=data.frame(wB_Reg_Hetero_Rand_SP1,wB_Reg_Hetero_Rand_pas.vol)

wB_Reg_Hetero_Matrix <- as.matrix(wB_Reg_Hetero[,c("C", "N", "P")])
wB_Reg_Hetero_Obs<-CHVind(wB_Reg_Hetero_Matrix)$vol

plot(wB_Reg_Hetero_Rand_SP1$mean.vol~wB_Reg_Hetero_Rand_SP1$wB_Reg_Hetero_Rand_pas.vol,pch=16,ylim=c(0,3.5),
     main = "wB Regional Heterospecific\\n n = 201 (5,200,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_Reg_Hetero_Rand_SP1$wB_Reg_Hetero_Rand_pas.vol),c(wB_Reg_Hetero_Rand_SP1$mean.vol-2*wB_Reg_Hetero_Rand_SP1$sd.vol),c(wB_Reg_Hetero_Rand_SP1$wB_Reg_Hetero_Rand_pas.vol),c(wB_Reg_Hetero_Rand_SP1$mean.vol+2*wB_Reg_Hetero_Rand_SP1$sd.vol))
abline(h=wB_Reg_Hetero_Obs)


##### wB Local CF Conspecific ##### 
wB_LocalCF_Cons_Rand_pas01=rep(seq(7,21,7),each=999) # n = 23
wB_LocalCF_Cons_Rand<-matrix(NA,ncol=1,nrow=length(wB_LocalCF_Cons_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(wB_LocalCF_Cons_Rand_pas01)){
  wB_LocalCF_Cons_Rand[i,1]<-round(CHVind(as.matrix(wB_LocalCF_Cons[sample(seq(1,length(wB_LocalCF_Cons[,1]),1),wB_LocalCF_Cons_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(wB_LocalCF_Cons_Rand,"D:/RandomNiches/wB_LocalCF_Cons_Rand.csv")
wB_LocalCF_Cons_Rand <- read.csv("RandomNiches/wB_LocalCF_Cons_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
wB_LocalCF_Cons_Rand_pas01=rep(seq(7,21,7),each=999)
wB_LocalCF_Cons_Rand_SP=cbind(wB_LocalCF_Cons_Rand,wB_LocalCF_Cons_Rand_pas01)
wB_LocalCF_Cons_Rand_SP$wB_LocalCF_Cons_Rand_pas01=as.factor(wB_LocalCF_Cons_Rand_SP$wB_LocalCF_Cons_Rand_pas01)
wB_LocalCF_Cons_Rand_SP1=wB_LocalCF_Cons_Rand_SP%>% 
  group_by(wB_LocalCF_Cons_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
wB_LocalCF_Cons_Rand_SP1=wB_LocalCF_Cons_Rand_SP1[!duplicated(wB_LocalCF_Cons_Rand_SP1$sd.vol), ]
wB_LocalCF_Cons_Rand_pas.vol=seq(7,21,7)
wB_LocalCF_Cons_Rand_SP1=data.frame(wB_LocalCF_Cons_Rand_SP1,wB_LocalCF_Cons_Rand_pas.vol)

wB_LocalCF_Cons_Matrix <- as.matrix(wB_LocalCF_Cons[,c("C", "N", "P")])
wB_LocalCF_Cons_Obs<-CHVind(wB_LocalCF_Cons_Matrix)$vol

plot(wB_LocalCF_Cons_Rand_SP1$mean.vol~wB_LocalCF_Cons_Rand_SP1$wB_LocalCF_Cons_Rand_pas.vol,pch=16,ylim=c(0,0.5),
     main = "wB Local CF Conspecific\\n n = 23 (7,21,7)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_LocalCF_Cons_Rand_SP1$wB_LocalCF_Cons_Rand_pas.vol),c(wB_LocalCF_Cons_Rand_SP1$mean.vol-2*wB_LocalCF_Cons_Rand_SP1$sd.vol),c(wB_LocalCF_Cons_Rand_SP1$wB_LocalCF_Cons_Rand_pas.vol),c(wB_LocalCF_Cons_Rand_SP1$mean.vol+2*wB_LocalCF_Cons_Rand_SP1$sd.vol))
abline(h=wB_LocalCF_Cons_Obs)


##### wB Local CF Heterospecific #####
wB_LocalCF_Hetero_Rand_pas01=rep(seq(4,48,4),each=999) # n = 48
wB_LocalCF_Hetero_Rand<-matrix(NA,ncol=1,nrow=length(wB_LocalCF_Hetero_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(wB_LocalCF_Hetero_Rand_pas01)){
  wB_LocalCF_Hetero_Rand[i,1]<-round(CHVind(as.matrix(wB_LocalCF_Hetero[sample(seq(1,length(wB_LocalCF_Hetero[,1]),1),wB_LocalCF_Hetero_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(wB_LocalCF_Hetero_Rand,"D:/RandomNiches/wB_LocalCF_Hetero_Rand.csv")
wB_LocalCF_Hetero_Rand <- read.csv("RandomNiches/wB_LocalCF_Hetero_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
wB_LocalCF_Hetero_Rand_pas01=rep(seq(4,48,4),each=999)
wB_LocalCF_Hetero_Rand_SP=cbind(wB_LocalCF_Hetero_Rand,wB_LocalCF_Hetero_Rand_pas01)
wB_LocalCF_Hetero_Rand_SP$wB_LocalCF_Hetero_Rand_pas01=as.factor(wB_LocalCF_Hetero_Rand_SP$wB_LocalCF_Hetero_Rand_pas01)
wB_LocalCF_Hetero_Rand_SP1=wB_LocalCF_Hetero_Rand_SP%>%
  group_by(wB_LocalCF_Hetero_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
wB_LocalCF_Hetero_Rand_SP1=wB_LocalCF_Hetero_Rand_SP1[!duplicated(wB_LocalCF_Hetero_Rand_SP1$sd.vol), ]
wB_LocalCF_Hetero_Rand_pas.vol=seq(4,48,4)
wB_LocalCF_Hetero_Rand_SP1=data.frame(wB_LocalCF_Hetero_Rand_SP1,wB_LocalCF_Hetero_Rand_pas.vol)

wB_LocalCF_Hetero_Matrix <- as.matrix(wB_LocalCF_Hetero[,c("C", "N", "P")])
wB_LocalCF_Hetero_Obs<-CHVind(wB_LocalCF_Hetero_Matrix)$vol

plot(wB_LocalCF_Hetero_Rand_SP1$mean.vol~wB_LocalCF_Hetero_Rand_SP1$wB_LocalCF_Hetero_Rand_pas.vol,pch=16,ylim=c(0,0.35),
     main = "wB Local CF Heterospecific\\n n = 48 (4,48,4)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_LocalCF_Hetero_Rand_SP1$wB_LocalCF_Hetero_Rand_pas.vol),c(wB_LocalCF_Hetero_Rand_SP1$mean.vol-2*wB_LocalCF_Hetero_Rand_SP1$sd.vol),c(wB_LocalCF_Hetero_Rand_SP1$wB_LocalCF_Hetero_Rand_pas.vol),c(wB_LocalCF_Hetero_Rand_SP1$mean.vol+2*wB_LocalCF_Hetero_Rand_SP1$sd.vol))
abline(h=wB_LocalCF_Hetero_Obs)



##### wB Local NP Conspecific #####
wB_LocalNP_Cons_Rand_pas01=rep(seq(4,5,0.25),each=999) # n = 5
wB_LocalNP_Cons_Rand<-matrix(NA,ncol=1,nrow=length(wB_LocalNP_Cons_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(wB_LocalNP_Cons_Rand_pas01)){
  wB_LocalNP_Cons_Rand[i,1]<-round(CHVind(as.matrix(wB_LocalNP_Cons[sample(seq(1,length(wB_LocalNP_Cons[,1]),1),wB_LocalNP_Cons_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(wB_LocalNP_Cons_Rand,"D:/RandomNiches/wB_LocalNP_Cons_Rand.csv")
wB_LocalNP_Cons_Rand <- read.csv("RandomNiches/wB_LocalNP_Cons_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
wB_LocalNP_Cons_Rand_pas01=rep(seq(4,5,0.25),each=999)
wB_LocalNP_Cons_Rand_SP=cbind(wB_LocalNP_Cons_Rand,wB_LocalNP_Cons_Rand_pas01)
wB_LocalNP_Cons_Rand_SP$wB_LocalNP_Cons_Rand_pas01=as.factor(wB_LocalNP_Cons_Rand_SP$wB_LocalNP_Cons_Rand_pas01)
wB_LocalNP_Cons_Rand_SP1=wB_LocalNP_Cons_Rand_SP%>%
  group_by(wB_LocalNP_Cons_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
wB_LocalNP_Cons_Rand_SP1=wB_LocalNP_Cons_Rand_SP1[!duplicated(wB_LocalNP_Cons_Rand_SP1$sd.vol), ]
wB_LocalNP_Cons_Rand_pas.vol=seq(4,5,0.25)
wB_LocalNP_Cons_Rand_SP1=data.frame(wB_LocalNP_Cons_Rand_SP1,wB_LocalNP_Cons_Rand_pas.vol)

wB_LocalNP_Cons_Matrix <- as.matrix(wB_LocalNP_Cons[,c("C", "N", "P")])
wB_LocalNP_Cons_Obs<-CHVind(wB_LocalNP_Cons_Matrix)$vol

plot(wB_LocalNP_Cons_Rand_SP1$mean.vol~wB_LocalNP_Cons_Rand_SP1$wB_LocalNP_Cons_Rand_pas.vol,pch=16,ylim=c(0,0.02),
     main = "wB Local NP Conspecific\\n n = 5 (4,5,0.25)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_LocalNP_Cons_Rand_SP1$wB_LocalNP_Cons_Rand_pas.vol),c(wB_LocalNP_Cons_Rand_SP1$mean.vol-2*wB_LocalNP_Cons_Rand_SP1$sd.vol),c(wB_LocalNP_Cons_Rand_SP1$wB_LocalNP_Cons_Rand_pas.vol),c(wB_LocalNP_Cons_Rand_SP1$mean.vol+2*wB_LocalNP_Cons_Rand_SP1$sd.vol))
abline(h=wB_LocalNP_Cons_Obs)


##### wB Local NP Heterospecific #####
wB_LocalNP_Hetero_Rand_pas01=rep(seq(9,153,9),each=999) # n = 153
wB_LocalNP_Hetero_Rand<-matrix(NA,ncol=1,nrow=length(wB_LocalNP_Hetero_Rand_pas01),dimnames=list(c(),"vol1"))
for(i in 1:length(wB_LocalNP_Hetero_Rand_pas01)){
  wB_LocalNP_Hetero_Rand[i,1]<-round(CHVind(as.matrix(wB_LocalNP_Hetero[sample(seq(1,length(wB_LocalNP_Hetero[,1]),1),wB_LocalNP_Hetero_Rand_pas01[i]),]))$vol,2)
  #print(i)
}
write.csv(wB_LocalNP_Hetero_Rand,"D:/RandomNiches/wB_LocalNP_Hetero_Rand.csv")
wB_LocalNP_Hetero_Rand <- read.csv("RandomNiches/wB_LocalNP_Hetero_Rand.csv", header = T)

# Calculate mean, SD, and SE for the 999 random niches
wB_LocalNP_Hetero_Rand_pas01=rep(seq(9,153,9),each=999)
wB_LocalNP_Hetero_Rand_SP=cbind(wB_LocalNP_Hetero_Rand,wB_LocalNP_Hetero_Rand_pas01)
wB_LocalNP_Hetero_Rand_SP$wB_LocalNP_Hetero_Rand_pas01=as.factor(wB_LocalNP_Hetero_Rand_SP$wB_LocalNP_Hetero_Rand_pas01)
wB_LocalNP_Hetero_Rand_SP1=wB_LocalNP_Hetero_Rand_SP%>%
  group_by(wB_LocalNP_Hetero_Rand_pas01)%>%
  mutate(mean.vol=mean(vol1),median.vol=median(vol1),sd.vol=sd(vol1),se.vol=se(vol1))
wB_LocalNP_Hetero_Rand_SP1=wB_LocalNP_Hetero_Rand_SP1[!duplicated(wB_LocalNP_Hetero_Rand_SP1$sd.vol), ]
wB_LocalNP_Hetero_Rand_pas.vol=seq(9,153,9)
wB_LocalNP_Hetero_Rand_SP1=data.frame(wB_LocalNP_Hetero_Rand_SP1,wB_LocalNP_Hetero_Rand_pas.vol)

wB_LocalNP_Hetero_Matrix <- as.matrix(wB_LocalNP_Hetero[,c("C", "N", "P")])
wB_LocalNP_Hetero_Obs<-CHVind(wB_LocalNP_Hetero_Matrix)$vol

plot(wB_LocalNP_Hetero_Rand_SP1$mean.vol~wB_LocalNP_Hetero_Rand_SP1$wB_LocalNP_Hetero_Rand_pas.vol,pch=16,ylim=c(0,2.5),
     main = "wB Local NP Heterospecific\\n n = 153 (9,153,9)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_LocalNP_Hetero_Rand_SP1$wB_LocalNP_Hetero_Rand_pas.vol),c(wB_LocalNP_Hetero_Rand_SP1$mean.vol-2*wB_LocalNP_Hetero_Rand_SP1$sd.vol),c(wB_LocalNP_Hetero_Rand_SP1$wB_LocalNP_Hetero_Rand_pas.vol),c(wB_LocalNP_Hetero_Rand_SP1$mean.vol+2*wB_LocalNP_Hetero_Rand_SP1$sd.vol))
abline(h=wB_LocalNP_Hetero_Obs)

# Run all plots above to make sure they work
# The code below here will compile the plots into a single plot for the manuscript
setwd("~/Dropbox/PhD_Thesis/Chapter_StoichNiche/Figures/")
tiff("SamplingEffect.tiff", width = 12, height = 7, units = "in", res = 600) 
par(mfrow=c(4,4), mar = c(2, 2.1, 2.5, 0.6), oma=c(2,2,0,0), mgp = c(2,1,0))

####### Plots #######
####### Ecoregion ######
# bF North Pen ecoregion
plot(bF_Ecoreg_NP_Rand_SP1$mean.vol~bF_Ecoreg_NP_Rand_SP1$bF_Ecoreg_NP_Rand_pas.vol,pch=16,ylim=c(0,0.6),
     main = "bF Northern Peninsula ecoregion\\n n = 295 (5,295,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_Ecoreg_NP_Rand_SP1$bF_Ecoreg_NP_Rand_pas.vol),c(bF_Ecoreg_NP_Rand_SP1$mean.vol-2*bF_Ecoreg_NP_Rand_SP1$sd.vol),c(bF_Ecoreg_NP_Rand_SP1$bF_Ecoreg_NP_Rand_pas.vol),c(bF_Ecoreg_NP_Rand_SP1$mean.vol+2*bF_Ecoreg_NP_Rand_SP1$sd.vol))
abline(h=bF_Ecoreg_NP_Obs)


# bF Central Forest ecoregion
plot(bF_Ecoreg_CF_Rand_SP1$mean.vol~bF_Ecoreg_CF_Rand_SP1$bF_Ecoreg_CF_Rand_pas.vol,pch=16,ylim=c(0,0.1),
     main = "bF Central Forest ecoregion\\n n = 95 (5,95,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_Ecoreg_CF_Rand_SP1$bF_Ecoreg_CF_Rand_pas.vol),c(bF_Ecoreg_CF_Rand_SP1$mean.vol-2*bF_Ecoreg_CF_Rand_SP1$sd.vol),c(bF_Ecoreg_CF_Rand_SP1$bF_Ecoreg_CF_Rand_pas.vol),c(bF_Ecoreg_CF_Rand_SP1$mean.vol+2*bF_Ecoreg_CF_Rand_SP1$sd.vol))
abline(h=bF_Ecoreg_CF_Obs)


# wB North Pen ecoregion
plot(wB_Ecoreg_NP_Rand_SP1$mean.vol~wB_Ecoreg_NP_Rand_SP1$wB_Ecoreg_NP_Rand_pas.vol,pch=16,ylim=c(0,2.5),
     main = "wB Northern Peninsula ecoregion\\n n = 158 (4,156,4)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_Ecoreg_NP_Rand_SP1$wB_Ecoreg_NP_Rand_pas.vol),c(wB_Ecoreg_NP_Rand_SP1$mean.vol-2*wB_Ecoreg_NP_Rand_SP1$sd.vol),c(wB_Ecoreg_NP_Rand_SP1$wB_Ecoreg_NP_Rand_pas.vol),c(wB_Ecoreg_NP_Rand_SP1$mean.vol+2*wB_Ecoreg_NP_Rand_SP1$sd.vol))
abline(h=wB_Ecoreg_NP_Obs)


# wB Central Forest ecoregion
plot(wB_Ecoreg_CF_Rand_SP1$mean.vol~wB_Ecoreg_CF_Rand_SP1$wB_Ecoreg_CF_Rand_pas.vol,pch=16,ylim=c(0,0.8),
     main = "wB Central Forest ecoregion\\n n = 71 (5,70,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_Ecoreg_CF_Rand_SP1$wB_Ecoreg_CF_Rand_pas.vol),c(wB_Ecoreg_CF_Rand_SP1$mean.vol-2*wB_Ecoreg_CF_Rand_SP1$sd.vol),c(wB_Ecoreg_CF_Rand_SP1$wB_Ecoreg_CF_Rand_pas.vol),c(wB_Ecoreg_CF_Rand_SP1$mean.vol+2*wB_Ecoreg_CF_Rand_SP1$sd.vol))
abline(h=wB_Ecoreg_CF_Obs)

###### Across ecoregions #####
# bF Regional Cons
plot(bF_Reg_Cons_Rand_SP1$mean.vol~bF_Reg_Cons_Rand_SP1$bF_Reg_Cons_Rand_pas.vol,pch=16,ylim=c(0,0.4),
     main = "bF across ecoregion: conspecific\\n n = 189 (7,189,7)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_Reg_Cons_Rand_SP1$bF_Reg_Cons_Rand_pas.vol),c(bF_Reg_Cons_Rand_SP1$mean.vol-2*bF_Reg_Cons_Rand_SP1$sd.vol),c(bF_Reg_Cons_Rand_SP1$bF_Reg_Cons_Rand_pas.vol),c(bF_Reg_Cons_Rand_SP1$mean.vol+2*bF_Reg_Cons_Rand_SP1$sd.vol))
abline(h=bF_Reg_Cons_Obs)

# bF Regional Hetero
plot(bF_Reg_Hetero_Rand_SP1$mean.vol~bF_Reg_Hetero_Rand_SP1$bF_Reg_Hetero_Rand_pas.vol,pch=16,ylim=c(0,0.6),
     main = "bF across ecoregion: heterospecific\\n n = 201 (5,200,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_Reg_Hetero_Rand_SP1$bF_Reg_Hetero_Rand_pas.vol),c(bF_Reg_Hetero_Rand_SP1$mean.vol-2*bF_Reg_Hetero_Rand_SP1$sd.vol),c(bF_Reg_Hetero_Rand_SP1$bF_Reg_Hetero_Rand_pas.vol),c(bF_Reg_Hetero_Rand_SP1$mean.vol+2*bF_Reg_Hetero_Rand_SP1$sd.vol))
abline(h=bF_Reg_Hetero_Obs)

# wB Regional Cons
plot(wB_Reg_Cons_Rand_SP1$mean.vol~wB_Reg_Cons_Rand_SP1$wB_Reg_Cons_Rand_pas.vol,pch=16,ylim=c(0,1),
     main = "wB across ecoregion: conspecific\\n n = 28 (4,28,1)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_Reg_Cons_Rand_SP1$wB_Reg_Cons_Rand_pas.vol),c(wB_Reg_Cons_Rand_SP1$mean.vol-2*wB_Reg_Cons_Rand_SP1$sd.vol),c(wB_Reg_Cons_Rand_SP1$wB_Reg_Cons_Rand_pas.vol),c(wB_Reg_Cons_Rand_SP1$mean.vol+2*wB_Reg_Cons_Rand_SP1$sd.vol))
abline(h=wB_Reg_Cons_Obs)

# wB Regional Hetero
plot(wB_Reg_Hetero_Rand_SP1$mean.vol~wB_Reg_Hetero_Rand_SP1$wB_Reg_Hetero_Rand_pas.vol,pch=16,ylim=c(0,3.5),
     main = "wB across ecoregion: heterospecific\\n n = 201 (5,200,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_Reg_Hetero_Rand_SP1$wB_Reg_Hetero_Rand_pas.vol),c(wB_Reg_Hetero_Rand_SP1$mean.vol-2*wB_Reg_Hetero_Rand_SP1$sd.vol),c(wB_Reg_Hetero_Rand_SP1$wB_Reg_Hetero_Rand_pas.vol),c(wB_Reg_Hetero_Rand_SP1$mean.vol+2*wB_Reg_Hetero_Rand_SP1$sd.vol))
abline(h=wB_Reg_Hetero_Obs)


####### Within/between ecoregion: NP ######
# bF Local NP Cons
plot(bF_LocalNP_Cons_Rand_SP1$mean.vol~bF_LocalNP_Cons_Rand_SP1$bF_LocalNP_Cons_Rand_pas.vol,pch=16,ylim=c(0,0.3),
     main = "bF within/between NP: conspecific\\n n = 142 (4,140,4)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_LocalNP_Cons_Rand_SP1$bF_LocalNP_Cons_Rand_pas.vol),c(bF_LocalNP_Cons_Rand_SP1$mean.vol-2*bF_LocalNP_Cons_Rand_SP1$sd.vol),c(bF_LocalNP_Cons_Rand_SP1$bF_LocalNP_Cons_Rand_pas.vol),c(bF_LocalNP_Cons_Rand_SP1$mean.vol+2*bF_LocalNP_Cons_Rand_SP1$sd.vol))
abline(h=bF_LocalNP_Cons_Obs)

# bF Local NP Hetero
plot(bF_LocalNP_Hetero_Rand_SP1$mean.vol~bF_LocalNP_Hetero_Rand_SP1$bF_LocalNP_Hetero_Rand_pas.vol,pch=16,ylim=c(0,0.5),
     main = "bF within/between NP: heterospecific\\n n = 153 (9,153,9)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_LocalNP_Hetero_Rand_SP1$bF_LocalNP_Hetero_Rand_pas.vol),c(bF_LocalNP_Hetero_Rand_SP1$mean.vol-2*bF_LocalNP_Hetero_Rand_SP1$sd.vol),c(bF_LocalNP_Hetero_Rand_SP1$bF_LocalNP_Hetero_Rand_pas.vol),c(bF_LocalNP_Hetero_Rand_SP1$mean.vol+2*bF_LocalNP_Hetero_Rand_SP1$sd.vol))
abline(h=bF_LocalNP_Hetero_Obs)

# wB Local NP cons
plot(wB_LocalNP_Cons_Rand_SP1$mean.vol~wB_LocalNP_Cons_Rand_SP1$wB_LocalNP_Cons_Rand_pas.vol,pch=16,ylim=c(0,0.02),
     main = "wB within/between NP: conspecific\\n n = 5 (4,5,0.25)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_LocalNP_Cons_Rand_SP1$wB_LocalNP_Cons_Rand_pas.vol),c(wB_LocalNP_Cons_Rand_SP1$mean.vol-2*wB_LocalNP_Cons_Rand_SP1$sd.vol),c(wB_LocalNP_Cons_Rand_SP1$wB_LocalNP_Cons_Rand_pas.vol),c(wB_LocalNP_Cons_Rand_SP1$mean.vol+2*wB_LocalNP_Cons_Rand_SP1$sd.vol))
abline(h=wB_LocalNP_Cons_Obs)

# wB Local NP Hetero
plot(wB_LocalNP_Hetero_Rand_SP1$mean.vol~wB_LocalNP_Hetero_Rand_SP1$wB_LocalNP_Hetero_Rand_pas.vol,pch=16,ylim=c(0,2.5),
     main = "wB within/between NP: heterospecific\\n n = 153 (9,153,9)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_LocalNP_Hetero_Rand_SP1$wB_LocalNP_Hetero_Rand_pas.vol),c(wB_LocalNP_Hetero_Rand_SP1$mean.vol-2*wB_LocalNP_Hetero_Rand_SP1$sd.vol),c(wB_LocalNP_Hetero_Rand_SP1$wB_LocalNP_Hetero_Rand_pas.vol),c(wB_LocalNP_Hetero_Rand_SP1$mean.vol+2*wB_LocalNP_Hetero_Rand_SP1$sd.vol))
abline(h=wB_LocalNP_Hetero_Obs)


####### Within/between ecoregion: CF ######

# bF Local CF Cons
plot(bF_LocalCF_Cons_Rand_SP1$mean.vol~bF_LocalCF_Cons_Rand_SP1$bF_LocalCF_Cons_Rand_pas.vol,pch=16,ylim=c(0,0.07),
     main = "bF within/between CF: conspecific\\n n = 47 (5,45,5)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_LocalCF_Cons_Rand_SP1$bF_LocalCF_Cons_Rand_pas.vol),c(bF_LocalCF_Cons_Rand_SP1$mean.vol-2*bF_LocalCF_Cons_Rand_SP1$sd.vol),c(bF_LocalCF_Cons_Rand_SP1$bF_LocalCF_Cons_Rand_pas.vol),c(bF_LocalCF_Cons_Rand_SP1$mean.vol+2*bF_LocalCF_Cons_Rand_SP1$sd.vol))
abline(h=bF_LocalCF_Cons_Obs)


# bF Local CF Hetero
plot(bF_LocalCF_Hetero_Rand_SP1$mean.vol~bF_LocalCF_Hetero_Rand_SP1$bF_LocalCF_Hetero_Rand_pas.vol,pch=16,ylim=c(0,0.07),
     main = "bF within/between CF: heterospecific\\n n = 48 (4,48,4)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(bF_LocalCF_Hetero_Rand_SP1$bF_LocalCF_Hetero_Rand_pas.vol),c(bF_LocalCF_Hetero_Rand_SP1$mean.vol-2*bF_LocalCF_Hetero_Rand_SP1$sd.vol),c(bF_LocalCF_Hetero_Rand_SP1$bF_LocalCF_Hetero_Rand_pas.vol),c(bF_LocalCF_Hetero_Rand_SP1$mean.vol+2*bF_LocalCF_Hetero_Rand_SP1$sd.vol))
abline(h=bF_LocalCF_Hetero_Obs)

# wB Local CF Cons
plot(wB_LocalCF_Cons_Rand_SP1$mean.vol~wB_LocalCF_Cons_Rand_SP1$wB_LocalCF_Cons_Rand_pas.vol,pch=16,ylim=c(0,0.5),
     main = "wB within/between CF: conspecific\\n n = 23 (7,21,7)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_LocalCF_Cons_Rand_SP1$wB_LocalCF_Cons_Rand_pas.vol),c(wB_LocalCF_Cons_Rand_SP1$mean.vol-2*wB_LocalCF_Cons_Rand_SP1$sd.vol),c(wB_LocalCF_Cons_Rand_SP1$wB_LocalCF_Cons_Rand_pas.vol),c(wB_LocalCF_Cons_Rand_SP1$mean.vol+2*wB_LocalCF_Cons_Rand_SP1$sd.vol))
abline(h=wB_LocalCF_Cons_Obs)

# wB Local CF Hetero
plot(wB_LocalCF_Hetero_Rand_SP1$mean.vol~wB_LocalCF_Hetero_Rand_SP1$wB_LocalCF_Hetero_Rand_pas.vol,pch=16,ylim=c(0,0.35),
     main = "wB within/between CF: heterospecific\\n n = 48 (4,48,4)", xlab="Number of Individuals",ylab=expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")))
segments(c(wB_LocalCF_Hetero_Rand_SP1$wB_LocalCF_Hetero_Rand_pas.vol),c(wB_LocalCF_Hetero_Rand_SP1$mean.vol-2*wB_LocalCF_Hetero_Rand_SP1$sd.vol),c(wB_LocalCF_Hetero_Rand_SP1$wB_LocalCF_Hetero_Rand_pas.vol),c(wB_LocalCF_Hetero_Rand_SP1$mean.vol+2*wB_LocalCF_Hetero_Rand_SP1$sd.vol))
abline(h=wB_LocalCF_Hetero_Obs)

mtext("Number of Individuals", side = 1, line=1, outer=T, cex=1)

mtext(expression(paste("Mean niche volumes",phantom(...)%+-%phantom(..),"2*SD")), side=2, line=0.5, outer=T, cex=1)

dev.off()


