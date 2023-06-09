# Install Packages----
install.packages("ggpubr")
install.packages ("devtools")
install.packages ("tidyverse")
install.packages ("ggplot2")
install.packages ("dplyr")
install.packages ("ape")
install.packages ("devtools")
install.packages ("phyloseq")
install.packages ("gridExtra")
install.packages("microbiome")
install.packages ("reshape2")
install.packages("stringr")
install.packages ("DESeq2")
install.packages ("ade4")
install.packages ("data.table")
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
install.packages("igraph")
install.packages("yaml")
install.packages("reshape2")
install.packages("doParallel")
install.packages("DT")
install.packages("exactRankTests")
install.packages("foreach")
install.packages("ggplot2", dependencies = TRUE)
install.packages("Rcpp")
install.packages("shiny")
install.packages("coin")
install.packages("reshape2")
install.packages("Tax4Fun")
install.packages("qiimer")
install.packages("rlang")
install.packages("httpuv")
install.packages("XML")
install.packages("MuMIn")
install.packages("plotrix")
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
install.packages("ihs")
install.packages("scales")
install.packages("glue")
remove.packages("vegan")
uninstall("vegan")
install_version("vegan", version = "2.4-4", repos = "https://cran.r-project.org")
install.packages("installr")
install.packages("microbiome")
install.packages("wesanderson")
install.packages("plotrix")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq", version = "3.8")

library(devtools) # Load the devtools package
install_github("microbiome/microbiome") 

install.packages("ggpubr")
library("ggpubr")

install.packages("ggord")
install.packages("sf")
install.packages("Rcpp")
install.packages("phylobase")
install.packages("chronos")
install.packages("adegenet")
install.packages("phyloseq")

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
install.packages('MicEco')

install.packages("remotes")
library('remotes')
remotes::install_github("Russel88/MicEco")
install.packages("mefa")
install.packages("pracma")

remotes::install_github("twbattaglia/btools")

BiocManager::install("genefilter")
install.packages("phytools")
devtools::install_github(repo = "UVic-omics/selbal")
install.packages("phylosignal")

BiocManager::install("DESeq2")
remotes::install_github("twbattaglia/btools")
install.packages('car')
install.packages("standardize")
install.packages('MDMR')
install.packages('pedantics')
install.packages('Biostrings')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
install.packages("AICcmodavg")
install.packages('sp')
install.packages('spaa')
# Load Packages----
library("ggplot2")
library("rlang")
library ("phyloseq")
library("phytools")
library ("vegan")
library("glue")
library ("devtools")
library("installr")
library("reshape2")
library("yaml")
library("ape")
library("biomformat")
library("ggpubr")
library("wesanderson")
library("microbiome")
library("MuMIn")
library("genefilter")
library("dplyr")
library('picante')
library('MicEco')
library('selbal')
library('plotrix')
library("phylosignal")
library("phylobase")
library("tidyverse")
library("phangorn")
library("adegenet")
library("pracma")
library('mefa')
library("btools")
library('car')
library("standardize")
library("MDMR")
library("lmerTest")
library("installr")
library("pedantics")
library("AICcmodavg")
library('gridExtra')
library('spaa')
# Build and Filter Phyloseq Object----

Sable_Horses <- read_csv2phyloseq(otu.file = "ASV_Counts.csv", taxonomy.file = "ASV_Taxonomyy.csv", metadata.file = "Sample_Metadata.csv")
Sable_Horses

ASVTree <- read.newick('ASV_Phylogenetic_Tree.tre')
Sable_Horses <- merge_phyloseq(Sable_Horses, ASVTree)

colnames (tax_table (Sable_Horses))

Sable_Clean <- subset_taxa(Sable_Horses, 
                           Kingdom %in% 
                             c( "Bacteria"))

flist <- filterfun(kOverA(4, 1))
Sable_Filtered <- filter_taxa(Sable_Clean, flist, TRUE) 
Sable <- Sable_Filtered
Sable


# Scaling Spatio-temporal Variables-----

sample_data(Sable)$IslandDist <- -1*sample_data(Sable)$Longitude + max(sample_data(Sable)$Longitude) #Reverse to re-orient data west to east (dataset expressed as degrees westing)

sample_data(Sable)$IslandDistStandard <- scale(-1*sample_data(Sable)$Longitude) #Rescale Longitude to a mean of 0 and std dev of 1.
mean(sample_data(Sable)$IslandDistStandard)
sd(sample_data(Sable)$IslandDistStandard)

sample_data(Sable)$Distance.From.Centre <- sample_data(Sable)$Longitude - 59.9194622594016 #2014 Mean longitude

sample_data(Sable)$Centre.Dist <- abs(sample_data(Sable)$Distance.From.Centre)
sample_data(Sable)$Centre.Dist.Std <- scale(abs(sample_data(Sable)$Distance.From.Centre)) #rescale day of year.

sample_data(Sable)$Julian.Day.Std <- scale(sample_data(Sable)$julian.date) #rescale day of year

sample_data(Sable)$Sex <- "Female" 
# Habitat Classes----
# Relative Veg Area Calculations

#Total area of buffer in m2
Buffer_Size_150 <- 3.14159*150^2

#Calculate total occupiable terrestrial area
sample_data (Sable)$TerrArea <- Buffer_Size_150 - sample_data (Sable)$Ocean_2019 - sample_data (Sable)$Building_2009

sample_data (Sable)$GrasslandTerr <- (sample_data (Sable)$Dense_Grassland_2009 + sample_data(Sable)$Sparse_Grassland_2009)/sample_data (Sable)$TerrArea

sample_data (Sable)$HeathTerr <- (sample_data (Sable)$Dense_Heath_2009 + sample_data(Sable)$Sparse_Heath_2009)/sample_data (Sable)$TerrArea

sample_data (Sable)$WaterTerr <- (sample_data(Sable)$Water_area)/sample_data (Sable)$TerrArea

sample_data (Sable)$BPTerr <- (sample_data(Sable)$Beach_Pea_2009)/sample_data (Sable)$TerrArea

sample_data(Sable)$Sandwort <- sample_data(Sable)$Dense_Sandwort_2009 + sample_data(Sable)$Sparse_Sandwort_2009

sample_data(Sable)$SandwortAccess[sample_data(Sable)$Sandwort > 0] <- "1"
sample_data(Sable)$SandwortAccess[sample_data(Sable)$Sandwort == 0] <- "0"

## Alpha Diversity----

set.seed(123) #Set seed for reproducibility of rarifaction
Horses_Rared <- rarefy_even_depth(Sable, rngseed = TRUE) #Rarefy
sample_data(Sable)$ObsSV<-(estimate_richness(Horses_Rared, measures="Observed", split=TRUE))[,1]
sample_data(Sable)$ObsSV

mean(sample_data(Sable)$ObsSV)
sd(sample_data(Sable)$ObsSV)
shapiro.test(sample_data(Sable)$ObsSV)

ggplot() + geom_point(data = sample_data(Sable), aes(x = IslandDistStandard, y = ObsSV, colour = status, shape = SandwortAccess))


options(na.action = "na.fail")
fm1 <- glm(sample_data(Sable)$ObsSV  ~ sample_data(Sable)$SandwortAccess + sample_data(Sable)$status + poly(sample_data(Sable)$age, 2, raw = FALSE)  + sample_data(Sable)$IslandDistStandard + sample_data(Sable)$BPTerr + sample_data(Sable)$GrasslandTerr + sample_data(Sable)$HeathTerr + sample_data(Sable)$Centre.Dist.Std + sample_data(Sable)$Julian.Day.Std )

summary(fm1)
ms1 <- dredge(fm1)
dredgetab <- as.data.frame(ms1)
View(dredgetab)

par(mar = c(3,5,6,4))
plot(ms1, labAsExpr = TRUE)
avgmod.95p <- model.avg(ms1, subset = delta < 3)
summary(avgmod.95p)
confint(avgmod.95p)

## Beta-Diversity Sequence Variant Level CLR ----
    #CLR Transformation
Sable_CLR <- microbiome::transform(Sable, 'clr')
Sable_CLR

sv.euc.ord <- ordinate (Sable_CLR, 
                             method = "PCoA", distance = "euclidean") #stress 0.19


sv.ord.plot <- plot_ordination (Sable_CLR, 
                 sv.euc.ord, 
                 axes = c(1,2)) + 
  geom_point (aes(colour=SandwortAccess), size = 3) +
  theme(axis.text = element_text(size=15, 
                                 colour="black"), #y-axis number colour
        panel.border = element_rect(colour="black", size=2, fill=NA),  
        axis.line = element_line(colour = "black", size=0.1), #axis line colour
        axis.ticks = element_line(colour = "black", size=1),
        panel.background = element_rect(fill=NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size=15, 
                                  colour="black"),
        panel.grid.major = element_line(colour = "#EEEEEE", size = 0.1),
        panel.grid.minor = element_line(colour = "#EEEEEE", size = 0.1), 
        plot.tag = element_text(size = 15)) +
  scale_colour_manual(name="sandwort", 
                        labels=c("absent", "present"),
                        values = c("#37412B", "#ECCE62")) +
  stat_ellipse(aes(colour=SandwortAccess), linetype = 2, size = 1)  + 
  labs(tag = "(A)") +
  xlab("PCA 1 [5.3%]")+
  ylab("PCA 2 [4.3%]")

sv.ord.plot

grid.arrange(sv.ord.plot, sv.dist, ncol = 2)




sv.euc <- phyloseq::distance (Sable_CLR, 
                                 "euclidean") 

  #Backwards term selection
adonis2 (sv.euc ~ Centre.Dist.Std + IslandDistStandard + status + age + Julian.Day.Std + SandwortAccess + HeathTerr + GrasslandTerr + BPTerr, as(sample_data (Sable_CLR), "data.frame"), permutations=9999, by = 'margin') 

adonis2 (sv.euc ~ Centre.Dist.Std + IslandDistStandard + status + poly(age, 2, raw=FALSE) + Julian.Day.Std + SandwortAccess + HeathTerr + BPTerr, as(sample_data (Sable_CLR), "data.frame"), permutations=9999, by = 'margin') 

adonis2 (sv.euc ~ Centre.Dist.Std + IslandDistStandard + status + Julian.Day.Std + SandwortAccess + HeathTerr + BPTerr, as(sample_data (Sable_CLR), "data.frame"), permutations=9999, by = 'margin') 

adonis2 (sv.euc ~ Centre.Dist.Std + IslandDistStandard + Julian.Day.Std + SandwortAccess + HeathTerr + BPTerr, as(sample_data (Sable_CLR), "data.frame"), permutations=9999, by = 'margin') 

adonis2 (sv.euc ~ Centre.Dist.Std + IslandDistStandard + Julian.Day.Std + SandwortAccess + HeathTerr, as(sample_data (Sable_CLR), "data.frame"), permutations=9999, by = 'margin') 

adonis2 (sv.euc ~ Centre.Dist.Std + IslandDistStandard + Julian.Day.Std + SandwortAccess, as(sample_data (Sable_CLR), "data.frame"), permutations=9999, by = 'margin') 


#Testing for a band effect
adonis2 (sv.euc ~ band, as(sample_data (Sable_CLR), "data.frame"), permutations=999, by = 'margin')

#Mantel test for effect of longitudinal separation
Longitudes <- as.data.frame(sample_data(Sable_CLR)$IslandDistStandard)
Longitude.Matrix <- as.matrix(dist(Longitudes,method = "euclidean",diag = TRUE,upper = FALSE))


mantel.test.longitude <- mantel(Longitude.Matrix, sv.euc, method = 'pearson', permutations = 9999)
mantel.test.longitude

#Convert distance matrix to list
ID <- as.character(sample_names(Sable))
Band <- as.character(sample_data(Sable)$band)
SW <- as.character(sample_data(Sable)$SandwortAccess)
Dist <- as.numeric(sample_data(Sable)$IslandDistStandard)

ID2 <- as.character(sample_names(Sable))
Band2 <- as.character(sample_data(Sable)$band)
SW2 <- as.character(sample_data(Sable)$SandwortAccess)
Dist2 <- as.numeric(sample_data(Sable)$IslandDistStandard)

MetaTable <- as.data.frame(cbind(ID, Band, SW, Dist))

MetaTable2 <- as.data.frame(cbind(ID2, Band2, SW2, Dist2))

Euc.Dist <- phyloseq::distance (Sable_CLR, 
                              "euclidean") 
Euc.List <- dist2list(Euc.Dist)
Euc.List$EucDist <- Euc.List$value
Euc.List$PairID <- with(Euc.List, paste0(row,col))

Pairwise.List <- merge(Euc.List, MetaTable, by.x="col", by.y="ID", all=FALSE)
Pairwise.List <- merge(Pairwise.List, MetaTable2, by.x="row", by.y="ID2", all=FALSE)

Pairwise.List$dDist <- abs(as.numeric(as.character(Pairwise.List$Dist)) - as.numeric(as.character(Pairwise.List$Dist2)))

Pairwise.List$SandwortD <- as.character(as.numeric(as.character(Pairwise.List$SW)) + as.numeric(as.character(Pairwise.List$SW2)))

sv.dist <- ggplot() + geom_point(data = Pairwise.List, 
                      aes(x = dDist, y = EucDist, colour = SandwortD),
                      size = 1.5) + 
  ylim(min = 90, max = 150) +
  theme(axis.text = element_text(size=15, 
                                 colour="black"), #y-axis number colour
        panel.border = element_rect(colour="black", size=2, fill=NA),  
        axis.line = element_line(colour = "black", size=0.1), #axis line colour
        axis.ticks = element_line(colour = "black", size=1),
        panel.background = element_rect(fill=NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size=15, 
                                  colour="black"),
        panel.grid.major = element_line(colour = "#EEEEEE", size = 0.1),
        panel.grid.minor = element_line(colour = "#EEEEEE", size = 0.1), 
        plot.tag = element_text(size = 15)) +
  scale_colour_manual(name="sandwort", 
                      labels=c("absent", "both", "present"),
                      values = c("#37412B", "#6C8844", "#ECCE62")) +
  geom_smooth(data = Pairwise.List, aes(x = dDist, y = EucDist), 
              method = lm, colour = "black", size = 1) +
  ylab("Euclidean Distance") +
  xlab("Longitudinal Separation") + 
  labs(tag = "(B)")





WuniFracPartialMerge <- merge(Sable.wuni2, MetaTable, by.x="row", by.y="ID", all=FALSE)

Merged_Sable_Wunifrac <- merge(WuniFracPartialMerge, MetaTable2, by.x="col", by.y="ID2", all=FALSE)

Merged_Sable_Wunifrac$wunifracdist <- Merged_Sable_Wunifrac$dist
Merged_Sable_Wunifrac$PairID <- with(Merged_Sable_Wunifrac, paste0(row,col))


#WUnifrac Beta-Dist----
wuni.ord <- ordinate (Horses_Rared, 
                      method = "PCoA", distance = "wunifrac") 

plot_ordination (Horses_Rared, 
                 wuni.ord, 
                 axes = c(1,2)) + 
  geom_point (aes(colour=SandwortAccess), size = 4) +
  theme(axis.text = element_text(size=15, 
                                 colour="black"), #y-axis number colour
        panel.border = element_rect(colour="black", size=2, fill=NA),  
        axis.line = element_line(colour = "black", size=0.1), #axis line colour
        axis.ticks = element_line(colour = "black", size=1),
        panel.background = element_rect(fill=NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size=15, 
                                  colour="black"),
        panel.grid.major = element_line(colour = "#EEEEEE", size = 0.1),
        panel.grid.minor = element_line(colour = "#EEEEEE", size = 0.1), 
        plot.tag = element_text(size = 15)) +
  scale_colour_manual(name="sandwort", 
                      labels=c("absent", "present"),
                      values = c("#37412B", "#ECCE62")) +
  stat_ellipse(aes(colour=SandwortAccess), linetype = 2, size = 1) + 
  xlab("PCA 1 [16.9%]") +
  ylab("PCA 2 [10.7%]")

Sable.wuni <- phyloseq::distance (Horses_Rared, 
                                  "wunifrac") 

adonis2 (Sable.wuni ~ Centre.Dist.Std + IslandDistStandard + poly(age, 2, raw = TRUE) + Julian.Day.Std + SandwortAccess + HeathTerr + GrasslandTerr + BPTerr + status, as(sample_data (Sable), "data.frame"), permutations=9999, by = 'margin') 

adonis2 (Sable.wuni ~ Centre.Dist.Std + IslandDistStandard + poly(age, 2, raw = TRUE) + Julian.Day.Std + SandwortAccess + HeathTerr + GrasslandTerr + BPTerr, as(sample_data (Sable), "data.frame"), permutations=9999, by = 'margin') 

adonis2 (Sable.wuni ~ Centre.Dist.Std + IslandDistStandard + Julian.Day.Std + SandwortAccess + HeathTerr + GrasslandTerr + BPTerr, as(sample_data (Sable), "data.frame"), permutations=9999, by = 'margin') 

adonis2 (Sable.wuni ~ Centre.Dist.Std + IslandDistStandard + Julian.Day.Std + SandwortAccess + HeathTerr + BPTerr, as(sample_data (Sable), "data.frame"), permutations=9999, by = 'margin') 

adonis2 (Sable.wuni ~ Centre.Dist.Std + IslandDistStandard + Julian.Day.Std + SandwortAccess + BPTerr, as(sample_data (Sable), "data.frame"), permutations=9999, by = 'margin') 

adonis2 (Sable.wuni ~ IslandDistStandard + Julian.Day.Std + SandwortAccess + BPTerr, as(sample_data (Sable), "data.frame"), permutations=9999, by = 'margin') 

#Test for a band effect
adonis2 (Sable.wuni ~ band, as(sample_data (Sable), "data.frame"), permutations=999)


#Mantel test for effect of longitudinal separation
Longitudes <- as.data.frame(sample_data(Sable_CLR)$IslandDistStandard)
Longitude.Matrix <- as.matrix(dist(Longitudes,method = "euclidean",diag = TRUE,upper = FALSE))

mantel.test.longitude <- mantel(Longitude.Matrix, Sable.wuni, method = 'pearson', permutations = 9999)
mantel.test.longitude

#Plot
Sable.wuni <- phyloseq::distance (Horses_Rared, 
                                  "wunifrac") 

Wuni.List <- dist2list(Sable.wuni)
Wuni.List$WuniDist <- Wuni.List$value
Wuni.List$PairID <- with(Wuni.List, paste0(row,col))

Pairwise.List <- merge(Pairwise.List, Wuni.List, by.x="PairID", by.y="PairID", all=FALSE)

ggplot() + 
  geom_point(data = Pairwise.List, 
             aes(x = dDist, y = WuniDist, colour = SandwortD)) +
  theme(axis.text.x = element_text(size=12, 
                                   colour="black"), #x-axis number colour
        axis.text.y = element_text(size=12, 
                                   colour="black"), #y-axis number colour
        panel.border = element_rect(colour="black", size=2, fill=NA), #plot frame colour 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#4B6E79", size=0.1), #axis line colour
        axis.ticks = element_line(colour = "#4B6E79", size=1.5),
        panel.background = element_rect(fill=NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "top") +
  scale_shape_discrete(labels=c("Neither horse\\nwith sandwort", "One horse\\nwith sandwort", "Both horses\\nwith sandwort")) 


# Beta-Dispersion----

  #Ordination Plots
wuni <- phyloseq::distance (Horses_Rared, 
                               "wunifrac")
sable_sample_data <- data.frame(sample_data(Sable))
beta.wuni <- betadisper(wuni, sable_sample_data$Sex, type = "centroid", bias.adjust = TRUE)
beta.wuni$distances
dispers.wuni <- data.frame(beta.wuni$distances)
sample_data(Sable)$FracDisp <- dispers.wuni[,1]

CLR_SV <- microbiome::transform(Sable, 'clr')
sv.euc <- phyloseq::distance (Sable_CLR, 
                               "euclidean")
sable_sample_data <- data.frame(sample_data(Sable))
beta.sv <- betadisper(sv.euc, sable_sample_data$Sex, type = "centroid", bias.adjust = TRUE)
beta.sv$distances
dispers.sv <- data.frame(beta.sv$distances)
beta.sv$distances

sample_data(Sable)$SVDisp <- dispers.sv[,1]
sample_data(Sable)$FracDisp <- dispers.wuni[,1]


  # SV Beta-Dispersion
options(na.action = "na.fail")

fm1 <- glm(sample_data(Sable)$SVDisp ~ sample_data(Sable)$SandwortAccess + sample_data(Sable)$status + poly(sample_data(Sable)$age, 2, raw = FALSE) + sample_data(Sable)$IslandDistStandard + sample_data(Sable)$BPTerr + sample_data(Sable)$GrasslandTerr  + sample_data(Sable)$Centre.Dist.Std + sample_data(Sable)$HeathTerr + sample_data(Sable)$Julian.Day.Std)
ms1 <- dredge(fm1)
dredgetab <- as.data.frame(ms1)
View(dredgetab)
summary(fm1)

par(mar = c(3,5,6,4))
plot(ms1, labAsExpr = TRUE)
avgmod.95p <- model.avg(ms1, subset = delta < 3)
summary(avgmod.95p)
confint(avgmod.95p)

# Weighted UniFrac Beta-Dispersion
options(na.action = "na.fail")

fm1 <- glm(log(sample_data(Sable)$FracDisp) ~ sample_data(Sable)$SandwortAccess + sample_data(Sable)$status + poly(sample_data(Sable)$age, 2, raw = FALSE) + sample_data(Sable)$IslandDistStandard + sample_data(Sable)$BPTerr + sample_data(Sable)$GrasslandTerr  + sample_data(Sable)$Centre.Dist.Std + sample_data(Sable)$HeathTerr + sample_data(Sable)$Julian.Day.Std)
summary(fm1)
ms1 <- dredge(fm1)
dredgetab2 <- as.data.frame(ms1)
View(dredgetab2)

par(mar = c(3,5,6,4))
plot(ms1, labAsExpr = TRUE)
avgmod.95p <- model.avg(ms1, subset = delta < 3)
summary(avgmod.95p)
confint(avgmod.95p)





#Create a Phylogenetic Object with an Ultrametric Tree----
Sable_Ultra <- read_csv2phyloseq(otu.file = "ASVs_counts.csv", taxonomy.file = "ASVs_taxonomy.csv", metadata.file = "Sample_Metadata.csv")
Sable_Ultra

ASVTree <- read.newick("ASV_Phylogenetic_Tree.tre")
ChronoASVTree <- chronos(ASVTree)
Sable_Ultra <- merge_phyloseq(Sable_Ultra, ChronoASVTree)

colnames (tax_table (Sable_Ultra))
colnames (tax_table (Sable_Ultra)) <- c ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rank_names (Sable_Ultra) #Confirm
Sable_Ultra <- subset_taxa(Sable_Ultra, 
                           Kingdom %in% 
                             c( "Bacteria")) 

mean(sample_sums(Sable_Ultra))
flist <- filterfun(kOverA(4, 1))
Sable_Ultra <- filter_taxa(Sable_Ultra, flist, TRUE) 
Sable_Ultra <- Sable_Ultra
Sable_Ultra

sample_data(Sable_Ultra)$Sandwort <- sample_data(Sable_Ultra)$Dense_Sandwort_Prop + sample_data(Sable_Ultra)$Sparse_Sandwort_Prop

sample_data(Sable_Ultra)$SandwortAccess[sample_data(Sable_Ultra)$Sandwort > 0] <- "1"
sample_data(Sable_Ultra)$SandwortAccess[sample_data(Sable_Ultra)$Sandwort == 0] <- "0"

set.seed(123)
Sable_Ultra_R <- rarefy_even_depth(Sable_Ultra, rngseed = TRUE)

#Phylosignal----

set.seed(123) #Set seed for reproducibility of rarifaction
Sable_PhyloSignalR <- rarefy_even_depth(Sable_Ultra, rngseed = TRUE)


sample_data(Sable_PhyloSignalR)$SandwortNiche[sample_data(Sable_PhyloSignalR)$SandwortAccess == "1"] <- 1
sample_data(Sable_PhyloSignalR)$SandwortNiche[sample_data(Sable_PhyloSignalR)$SandwortAccess == "0"] <- -1

names(phy_tree(Sable_PhyloSignalR))
ASV_Ultra_Tree <- phy_tree(Sable_PhyloSignalR)
ASV_Proportional <- otu_table(Sable_PhyloSignalR)/taxa_sums(Sable_PhyloSignalR)


asv_sandwortniche <- ASV_Proportional * sample_data(Sable_PhyloSignal)$SandwortNiche
avg_ASV_sandwortniche <- rowSums(asv_sandwortniche)
hist(avg_ASV_sandwortniche)

InputSandwortPhylo4d <- phylo4d(ASV_Ultra_Tree, tip.data = avgASVsandwortniche)

test_pc <- phyloCorrelogram(InputSandwort)
plot(pc_Sable_Sandwort_Ultra)

#NTI----
Sable_Ultra_Merged <- merge_samples(Sable_Ultra, sample_data(Sable_Ultra)$Name, fun = mean)

set.seed(123)
Sable_Ultra_Rared <- rarefy_even_depth(Sable_Ultra_Merged, rngseed = TRUE)

write.csv(otu_table(Sable_Ultra_Rared), "Rared_Ultra_OTUs.csv")
Sable_Ultra_ASVs <- read.csv("Rared_Ultra_ASVs.csv", row.names = 1)

Sable_MNTD_Matrix_SES <- ses.mntd(Sable_Ultra_ASVs, cophenetic(phy_tree(Sable_Ultra_Rared)), abundance.weighted=TRUE, null.model = "taxa.labels")

Sable_MNTD_Matrix_SES$NTI <- Sable_MNTD_Ultra_SES$mntd.obs.z

Sable_MNTD_Matrix_SES$Name <- row.names(Sable_MNTD_Matrix_SES)

write.csv(sample_data(Sable), "Micro_DataTable.csv")
Micro_Samples <- read.csv("Micro_DataTable.csv")

Sable_NTI_Merged <- merge(Micro_Samples, Sable_MNTD_Matrix_SES, by.x="Name.x", by.y="Name", all=FALSE)
Sable_NTI_Merged$SandwortAccess <- as.character(Sable_NTI_Merged$SandwortAccess)

#Multi-Model Inference
options(na.action = "na.fail")

fm1 <- glm(NTI ~ SandwortAccess + status + poly(age, 2, raw = FALSE) + IslandDistStandard + BPTerr + Centre.Dist.Std + HeathTerr + GrasslandTerr + Julian.Day.Std, data = Sable_NTI_Merged)
summary(fm1)

ms1 <- dredge(fm1)
dredgetab <- as.data.frame(ms1)
View(dredgetab)

par(mar = c(3,5,6,4))
plot(ms1, labAsExpr = TRUE)
avgmod.95p <- model.avg(ms1, subset = delta < 3)
summary(avgmod.95p)
confint(avgmod.95p)

#BMNTD----
#Calculate BMNTD values

Sable_BMNTD_Ultra_SES <- ses.comdistnt(Sable_Ultra_ASVs, cophenetic(phy_tree(Sable_Ultra_Rared)), abundance.weighted=TRUE, null.model = "taxa.labels")

#All distances need to be positive for PERMANOVA analysis
BMNTD_Matrix_Pos <- Sable_BMNTD_Ultra_SES$comdistnt.obs.z + 3.5

#Ensure metadata is in the proper order
Sable_Data <- as.data.frame(sample_data(Sable))
Sable_Data_Reorder <- Sable_Data %>% arrange(Name.x)

#Backwards selection of terms 
adonis2(BMNTD_Matrix_Pos ~  GrasslandTerr + SandwortAccess + BPTerr + Centre.Dist.Std + IslandDistStandard + status + poly(age, 2, raw=TRUE) + HeathTerr + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=9999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matrix_Pos ~  GrasslandTerr + SandwortAccess + BPTerr + Centre.Dist.Std + poly(age, 2, raw=TRUE) + HeathTerr + status + Julian.Day.Std + WaterTerr, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matrix_Pos ~  GrasslandTerr + SandwortAccess + BPTerr + Centre.Dist.Std + poly(age, 2, raw=TRUE) + HeathTerr + status + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matrix_Pos ~  GrasslandTerr + SandwortAccess + BPTerr + Centre.Dist.Std + poly(age, 2, raw=TRUE) + HeathTerr + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matri_Pos ~  GrasslandTerr + SandwortAccess + BPTerr + poly(age, 2, raw=TRUE) + HeathTerr + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matrix_Pos ~  GrasslandTerr + SandwortAccess + BPTerr + HeathTerr + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

#Test for a band effect
adonis2(BMNTD_Matrix_Pos ~  band, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

#Consider absolute deviations
BMNTD_Matrix_Abs <- abs(Sable_BMNTD_Ultra_SES$comdistnt.obs.z)

adonis2(BMNTD_Matrix_Abs ~ GrasslandTerr + SandwortAccess + BPTerr + Centre.Dist.Std + IslandDistStandard + poly(age, 2, raw=TRUE) + HeathTerr + status + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=9999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matrix_Abs ~ GrasslandTerr + BPTerr + Centre.Dist.Std + IslandDistStandard + poly(age, 2, raw=TRUE) + HeathTerr + status + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matrix_Abs ~ GrasslandTerr + BPTerr + Centre.Dist.Std + IslandDistStandard + poly(age, 2, raw=TRUE) + HeathTerr + status, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matrix_Abs ~ GrasslandTerr + BPTerr + IslandDistStandard + poly(age, 2, raw=TRUE) + HeathTerr + status, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matrix_Abs ~ BPTerr + IslandDistStandard + poly(age, 2, raw=TRUE) + HeathTerr + status, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matrix_Abs ~ BPTerr + IslandDistStandard + poly(age, 2, raw=TRUE) + status, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(BMNTD_Matrix_Abs ~ BPTerr + IslandDistStandard + status, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

#Test for band effects
adonis2(BMNTD_Matrix_Abs ~ band, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

#Plotting BMNTD values
BMNTD_Matrix <- Sable_BMNTD_Ultra_SES$comdistnt.obs.z
Sable_BMNTD <- BMNTD_Matrix
Sable_BMNTD[lower.tri(Sable_BMNTD,diag=TRUE)]=NA 

Sable_BMNTD<-as.data.frame(as.table(Sable_BMNTD)) # as a dataframe
Sable_BMNTD<-na.omit(Sable_BMNTD)# remove NA


ID <- as.character(sample_names(Sable))
Name <- as.character(sample_data(Sable)$Name)
status <- as.character(sample_data(Sable)$status)
band <- as.character(sample_data(Sable)$band)
Dist <- as.numeric(sample_data(Sable)$IslandDistStandard)
SandwortCode <- as.numeric(as.character(sample_data(Sable)$SandwortAccess))
HeathCode <- as.numeric(sample_data(Sable)$HeathAccess)
Jul <- as.numeric(sample_data(Sable)$Julian.Day.Std)
Grass <- as.numeric(sample_data(Sable)$GrasslandTerr)
Heath <- as.numeric(sample_data(Sable)$HeathTerr)
BP <- as.numeric(sample_data(Sable)$BPTerr)

ID2 <- as.character(sample_names(Sable))
Name2 <- as.character(sample_data(Sable)$Name)
status2 <- as.character(sample_data(Sable)$status)
band2 <- as.character(sample_data(Sable)$band)
Dist2 <- as.numeric(sample_data(Sable)$IslandDistStandard)
SandwortCode2 <- as.numeric(as.character(sample_data(Sable)$SandwortAccess))
HeathCode2 <- as.numeric(sample_data(Sable)$HeathAccess)
Grass2 <- as.numeric(sample_data(Sable)$GrasslandTerr)
Heath2 <- as.numeric(sample_data(Sable)$HeathTerr)
Jul2 <- as.numeric(sample_data(Sable)$Julian.Day.Std)
BP2 <- as.numeric(sample_data(Sable)$BPTerr)

MetaTable <- as.data.frame(cbind(ID, Name, status, band, Dist, band, SandwortCode, age, Grass, Heath, Jul, BP))

MetaTable2 <- as.data.frame(cbind(ID2, Name2, status2, band2, Dist2,SandwortCode2, Grass2, Heath2, Jul2, BP2))


Merged_Sable_Partial <- merge(Sable_BMNTD, MetaTable, by.x="Var1", by.y="Name", all=FALSE)

Merged_Sable <- merge(Merged_Sable_Partial, MetaTable2, by.x="Var2", by.y="Name2", all=FALSE)

Merged_Sable$Dist <- as.numeric(as.character(Merged_Sable$Dist))
Merged_Sable$Dist2 <- as.numeric(as.character(Merged_Sable$Dist2))
Merged_Sable$aDist <- (as.numeric(as.character(Merged_Sable$Dist2)) + as.numeric(as.character(Merged_Sable$Dist)))/2
Merged_Sable$dDist <- abs((as.numeric(as.character(Merged_Sable$Dist2)) - as.numeric(as.character(Merged_Sable$Dist))))

Merged_Sable$BandMatch[Merged_Sable$band == Merged_Sable$band2] <- "Within"
Merged_Sable$BandMatch[Merged_Sable$band != Merged_Sable$band2] <- "Between"
Merged_Sable$SandwortD <- as.character(as.numeric(as.character(Merged_Sable$SandwortCode)) + as.numeric(as.character(Merged_Sable$SandwortCode2)))

Merged_Sable$aBP <- (as.numeric(as.character(Merged_Sable$BP)) + as.numeric(as.character(Merged_Sable$BP2)))/2
Merged_Sable$aHeath <- (as.numeric(as.character(Merged_Sable$Heath)) + as.numeric(as.character(Merged_Sable$Heath2)))/2
Merged_Sable$aGrass <- (as.numeric(as.character(Merged_Sable$Grass2)) + as.numeric(as.character(Merged_Sable$Grass)))/2
Merged_Sable$aJul <- (as.numeric(as.character(Merged_Sable$Jul2)) + as.numeric(as.character(Merged_Sable$Jul)))/2

Merged_Sable$StatusPairing <- with(Merged_Sable, paste0(status, status2))
Merged_Sable$StatusPairing [Merged_Sable$StatusPairing == "no_foalfoal"] <- "foalno_foal"


ggplot() + geom_point(data = Merged_Sable, 
                      aes(x = aGrass, 
                          y = Freq, 
                          colour = SandwortD),
                      alpha = 0.6,
                      size = 2) +
  xlab("Average Relative Area of Grassland") +
  ylab("Î²MNTD (effect size standardized)") +
  theme(axis.text = element_text(size=15, 
                                       colour="black"), #y-axis number colour
              panel.border = element_rect(colour="black", size=2, fill=NA),  
              axis.line = element_line(colour = "black", size=0.1), #axis line colour
              axis.ticks = element_line(colour = "black", size=1),
              panel.background = element_rect(fill=NA),
              plot.background = element_rect(fill="white"),
              legend.background = element_blank(),
              plot.title = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              legend.position = "none",
              axis.title = element_text(size=15, 
                                        colour="black"),
              panel.grid.major = element_line(colour = "#EEEEEE", size = 0.1),
              panel.grid.minor = element_line(colour = "#EEEEEE", size = 0.1)) + 
  scale_colour_manual(values = c("#37412B", "#6C8844", "#ECCE62")) + 
  facet_wrap(~SandwortD) +
  geom_smooth(data = Merged_Sable, aes(x = aGrass,y = Freq, group = ID), colour = "grey", se = F, method = lm, size = 0.5) +
  geom_smooth(data = Merged_Sable, aes(x = aGrass,y = Freq, colour = SandwortD), colour = "black", se = F, method = lm, size = 1.25) 



bpplot <- ggplot() + geom_point(data = Merged_Sable, aes(x = aBP, y = abs(Freq), colour = SandwortD), alpha = 0.8) + geom_smooth(data = Merged_Sable, aes(x = aBP, y = abs(Freq), group=ID), colour = "grey", se = F, method = lm, size = 0.5)+ xlab("average relative area of beach pea") +
  ylab("|BMNTD|\\n(effect size standardized)") + facet_wrap(~SandwortD) + geom_smooth(data = Merged_Sable, aes(x = aBP, y = abs(Freq)), se = T, colour = "black", method = lm) +
  theme(axis.text = element_text(size=15, 
                                 colour="black"), #y-axis number colour
        panel.border = element_rect(colour="black", size=2, fill=NA),  
        axis.line = element_line(colour = "black", size=0.1), #axis line colour
        axis.ticks = element_line(colour = "black", size=1),
        panel.background = element_rect(fill=NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "right",
        axis.title = element_text(size=12, 
                                  colour="black"),
        panel.grid.major = element_line(colour = "#EEEEEE", size = 0.1),
        panel.grid.minor = element_line(colour = "#EEEEEE", size = 0.1)) + 
  scale_colour_manual(values = c("#37412B", "#6C8844", "#ECCE62")) +
  scale_colour_discrete(name = "Sandwort Presence", labels=c("Neither horse\\nwith sandwort", "One horse\\nwith sandwort", "Both horses\\nwith sandwort")) 



heathplot <- ggplot() + geom_point(data = Merged_Sable, aes(x = aHeath, y = Freq, colour = SandwortD), alpha = 0.8) + geom_smooth(data = Merged_Sable, aes(x = aHeath, y = Freq, group=ID), colour = "grey", se = F, method = lm, size = 0.5)+ xlab("average relative area of heathland") +
  ylab("BMNTD\\n(effect size standardized)") + facet_wrap(~SandwortD) + geom_smooth(data = Merged_Sable, aes(x = aHeath, y = Freq), se = T, colour = "black", method = lm) +
  theme(axis.text = element_text(size=15, 
                                 colour="black"), #y-axis number colour
        panel.border = element_rect(colour="black", size=2, fill=NA),  
        axis.line = element_line(colour = "black", size=0.1), #axis line colour
        axis.ticks = element_line(colour = "black", size=1),
        panel.background = element_rect(fill=NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "right",
        axis.title = element_text(size=12, 
                                  colour="black"),
        panel.grid.major = element_line(colour = "#EEEEEE", size = 0.1),
        panel.grid.minor = element_line(colour = "#EEEEEE", size = 0.1)) + 
  scale_colour_manual(values = c("#37412B", "#6C8844", "#ECCE62")) +
  scale_colour_discrete(name = "Sandwort Presence", labels=c("Neither horse\\nwith sandwort", "One horse\\nwith sandwort", "Both horses\\nwith sandwort")) 

heathplot


statusplot <- ggplot() + geom_point(data = Merged_Sable, aes(x = aDist, y = abs(Freq), colour = StatusPairing), alpha = 0.8) + geom_smooth(data = Merged_Sable, aes(x = aDist, y = abs(Freq), group=ID), colour = "grey", se = F, method = lm, size = 0.5)+ xlab("average longitude") +
  ylab("|BMNTD|\\n(effect size standardized)") + facet_wrap(~StatusPairing) + geom_smooth(data = Merged_Sable, aes(x = aDist, y = abs(Freq)), se = T, colour = "black") +
  theme(axis.text = element_text(size=15, 
                                 colour="black"), #y-axis number colour
        panel.border = element_rect(colour="black", size=2, fill=NA),  
        axis.line = element_line(colour = "black", size=0.1), #axis line colour
        axis.ticks = element_line(colour = "black", size=1),
        panel.background = element_rect(fill=NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "right",
        axis.title = element_text(size=12, 
                                  colour="black"),
        panel.grid.major = element_line(colour = "#EEEEEE", size = 0.1),
        panel.grid.minor = element_line(colour = "#EEEEEE", size = 0.1)) + 
  scale_colour_manual(values = c("#37412B", "#6C8844", "#ECCE62")) +
  scale_colour_discrete(name = "Status Pairing", labels=c("both horses\\nwith foal", "one horse\\nwith foal", "neither horse\\nwith foal")) 

statusplot


grid.arrange(heathplot, bpplot, statusplot, ncol=1)

#BMNTD Mantel Test

Sandwort <- as.data.frame(Sable_Data_Reorder$SandwortAccess)
SandwortMatrix <- as.matrix(dist(Sandwort,method = "euclidean",diag = TRUE,upper = FALSE))

Longs <- as.data.frame(Sable_Data_Reorder$IslandDistStandard)
LongMatrix <- as.matrix(dist(Longs,method = "euclidean",diag = TRUE,upper = FALSE))


BMNTD_Matrix <- Sable_BMNTD_Ultra_SES$comdistnt.obs.z

mantel.test <- mantel.partial(LongMatrix, BMNTD_Matrix, SandwortMatrix, method = 'spearman')
mantel.test 


#RCBray----
Sable_RC_bray_Matrix <- raup_crick_abu_par(Sable_ASVs, reps=999, ncore=4, classic_metric=FALSE, split_ties=TRUE)

hist(Sable_RC_bray_Matrix)

#Distances Must be positive for PERMANOVA Analysis
Sable_RC_bray_Matrix_Pos <- Sable_RC_bray_Matrix + 1

adonis2(Sable_RC_bray_Matrix_Pos ~ GrasslandTerr + SandwortAccess + BPTerr + Centre.Dist.Std + IslandDistStandard + poly(age, 2, raw=TRUE) + HeathTerr + status + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(Sable_RC_bray_Matrix_Pos ~ SandwortAccess + BPTerr + Centre.Dist.Std + IslandDistStandard + poly(age, 2, raw=TRUE) + HeathTerr + status + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(Sable_RC_bray_Matrix_Pos ~ SandwortAccess + BPTerr + Centre.Dist.Std + IslandDistStandard + poly(age, 2, raw=TRUE) + HeathTerr + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=999, na.rm=TRUE, by = 'margin')

adonis2(Sable_RC_bray_Matrix_Pos ~ SandwortAccess + BPTerr + Centre.Dist.Std + IslandDistStandard + Julian.Day.Std, as(Sable_Data_Reorder, "data.frame"), permutations=9999, na.rm=TRUE, by = 'margin')


#Test for band effect
adonis(Sable_RC_bray_Matrix_Pos ~ band, as(Sable_Data_Reorder, "data.frame"), permutations=9999, na.rm=TRUE, by = 'margin')

# Mantel Test to test for an effect of longitudinal distance
Longs <- as.data.frame(Sable_Data_Reorder$IslandDistStandard)
LongMatrix <- as.matrix(dist(Longs,method = "euclidean",diag = TRUE,upper = FALSE))

BMNTD_Matrix <- Sable_BMNTD_Ultra_SES$comdistnt.obs.z

mantel.test <- mantel.partial(LongMatrix, Sable_RC_bray_Matrix, BMNTD_Matrix, method = 'spearman')
mantel.test 

Sable_RC_bray_Matrix_Manip <- Sable_RC_bray_Matrix

Sable_RC_bray_Matrix_Manip[lower.tri(Sable_RC_bray_Matrix_Manip,diag=TRUE)]=NA 
Sable_RC_bray_Matrix_Manip<-as.data.frame(as.table(Sable_RC_bray_Matrix_Manip)) # as a dataframe
Sable_RC_bray_Matrix_Manip<-na.omit(Sable_RC_bray_Matrix_Manip) # remove NA

Merged_Sable$PairedID <- with(Merged_Sable, paste0(Var1, Var2))
Sable_RC_bray_Matrix_Manip$PairedID <- with(Sable_RC_bray_Matrix_Manip, paste0(Var1, Var2))
Sable_Merged_BMNTD_RC <- merge(Merged_Sable, Sable_RC_bray_Matrix_Manip, by.x="PairedID", by.y="PairedID", all=FALSE)

Sable_Merged_BMNTD_RC$FreqType[abs(Sable_Merged_BMNTD_RC$Freq.x) >= 2] <- "Deterministic"
Sable_Merged_BMNTD_RC$FreqType[abs(Sable_Merged_BMNTD_RC$Freq.x) < 2] <- "Neutral"


Sable_Merged_BMNTD_RC$RCBi[Sable_Merged_BMNTD_RC$Freq.y >= 0.95] <- "1"
Sable_Merged_BMNTD_RC$RCBi[Sable_Merged_BMNTD_RC$Freq.y < 0.95] <- "0"
Sable_Merged_BMNTD_RC$RCBi <- as.numeric(Sable_Merged_BMNTD_RC$RCBi)

RC.Dist.Plot <- ggplot() + geom_jitter(data = Sable_Merged_BMNTD_RC, aes(x = dDist, y = Freq.y), size = 2, colour = "black", height = 0.05) +
  geom_smooth(data = Sable_Merged_BMNTD_RC, 
              method = "glm", 
              method.args = list(family = "binomial"),
              aes(x = dDist, 
                  y = RCBi,
                  linetype = FreqType), 
              size = 1.5,
              se = T, alpha = 0.3, 
              colour = "red") + 
  theme(axis.text = element_text(size=15, 
                                 colour="black"), #y-axis number colour
        panel.border = element_rect(colour="black", size=2, fill=NA),  
        axis.line = element_line(colour = "black", size=0.1), #axis line colour
        axis.ticks = element_line(colour = "black", size=1),
        panel.background = element_rect(fill=NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size=15, 
                                  colour="black"),
        panel.grid.major = element_line(colour = "#EEEEEE", size = 0.1),
        panel.grid.minor = element_line(colour = "#EEEEEE", size = 0.1),
        plot.tag = element_text(size = 15)) +
  ylab("Raup-Crick Bray Value") +
  xlab("Longitudinal Separation") +
  labs(tag = "(A)") 

RC.Dist.Plot

RC.Long.Plot <- ggplot() + geom_jitter(data = Sable_Merged_BMNTD_RC, aes(x = aDist, y = Freq.y), size = 2, colour = "black", height = 0.05) +
  geom_smooth(data = Sable_Merged_BMNTD_RC, 
              method = "glm", 
              method.args = list(family = "binomial"),
              aes(x = aDist, 
                  y = RCBi, 
                  colour = BandMatch,
                  fill = BandMatch), 
              size = 1.5,
              se = T, alpha = 0.3) + 
  theme(axis.text = element_text(size=15, 
                                 colour="black"), #y-axis number colour
        panel.border = element_rect(colour="black", size=2, fill=NA),  
        axis.line = element_line(colour = "black", size=0.1), #axis line colour
        axis.ticks = element_line(colour = "black", size=1),
        panel.background = element_rect(fill=NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size=15, 
                                  colour="black"),
        panel.grid.major = element_line(colour = "#EEEEEE", size = 0.1),
        panel.grid.minor = element_line(colour = "#EEEEEE", size = 0.1),
        plot.tag = element_text(size = 15)) +
  ylab("Raup-Crick Bray Value") +
  xlab("Average Longitude") + 
  scale_color_manual(name = "Band Membership", labels = c("Between", "Within"), values = c("#97A3D1", "#834686")) +
  scale_fill_manual(name = "Band Membership", labels = c("Between", "Within"), values = c("#97A3D1", "#834686")) +
  labs(tag = "(B)")
RC.Long.Plot

grid.arrange(RC.Dist.Plot, RC.Long.Plot, ncol = 2)




#Rarefaction Curve -----
library('plyr')

# Set seed so that rarefaction is reproducible
set.seed (123)

# Create function to subsample sequences at different depths
calculate_rarefaction_curves <- function (x, measures, depths) {
  # Create a function to subsample sequences at varying depths
  #
  # Args: 
  # x: A phyloseq, S4 format, object, with OTU sequence counts and sample metadata
  #
  # measures: An estimator of diversity
  #
  #depths: the size of sequence subsample
  #
  #Results:
  # A dataframe with variable sampling depth, sample id, diversity measure and diversity estimate
  #
  estimate_rarified_richness <- function (x, measures, depth) {
    if (max (sample_sums (x)) < depth) return ()
    x <- prune_samples (sample_sums (x) >= depth, x)
    
    rarified_Microbiome <- rarefy_even_depth (x, depth, verbose = FALSE)
    
    alpha_diversity_fecal <- estimate_richness (rarified_Microbiome, measures = c ("Observed"))
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity_fecal <- melt (as.matrix (alpha_diversity_fecal), varnames = c ('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity_fecal
  }
  
  names (depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data_fecal <- ldply (depths, estimate_rarified_richness, x = x, measures = measures, .id = 'Depth', .progress = ifelse (interactive (), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data_fecal$Depth <- as.numeric (levels (rarefaction_curve_data_fecal$Depth)) [rarefaction_curve_data_fecal$Depth]
  
  rarefaction_curve_data_fecal
}

Fecal_rarefaction <- calculate_rarefaction_curves (Horses_Rared, 
                                                   c ('Observed'), rep (c (5, 50, 500, 1000, 1500, 2000, 3000, 5000, 10000, 15000, 20000, 30000, 34280, 50000), each = 10))

# Average diversity across sampled replicates
Fecal_rarefaction <- ddply  (Fecal_rarefaction, 
                             c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

# Merge with sample data
Fecal_rarefaction <- merge (Fecal_rarefaction, data.frame (sample_data (Sable)), by.x = 'Sample', by.y = 'row.names')
head (Fecal_rarefaction)


# Get minimum sequence depth
fecal.sequence.depth <- min (sample_sums (Sable))
fecal.sequence.depth

# Plot rarefaction curve
fig <- ggplot() + 
  geom_line (data = Fecal_rarefaction, 
             aes (x = Depth, y = Alpha_diversity_mean, color = Sample), size=1) + 
  geom_vline (xintercept = fecal.sequence.depth, # Plots a vertical line at specified x - intercept
              linetype = 2) +# specifies plot aesthetics
  ylab ("# of ASVs") + # Assigns y - label
  xlab ("# of sequences per sample") +
  theme(axis.text = element_text(size=15, 
                                 colour="black"), #y-axis number colour
        panel.border = element_rect(colour="black", size=2, fill=NA),  
        axis.line = element_line(colour = "black", size=0.1), #axis line colour
        axis.ticks = element_line(colour = "black", size=1),
        panel.background = element_rect(fill=NA),
        plot.background = element_rect(fill="white"),
        legend.background = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size=15, 
                                  colour="black"),
        panel.grid.major = element_line(colour = "#EEEEEE", size = 0.1),
        panel.grid.minor = element_line(colour = "#EEEEEE", size = 0.1))

fig