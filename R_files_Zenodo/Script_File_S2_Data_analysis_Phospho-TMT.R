## # # # # # # # # # # # # # # # # # # # # # # # #
##
## Project: CDKL5-NLS Phosphoproteomics
##
## Script name: Script_File_S2_Data_analysis_Phospho-TMT.R
##
## Purpose of script: Process MaxQuant output and data analysis
##
## Author: Florian Weiland
##
## Date Created: 2020-04-03
##
## # # # # # # # # # # # # # # # # # # # # # # # #

## Libraries

suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(vsn, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(limma, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(seqinr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(plyr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(stringr, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggrepel, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggpointdensity, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(wesanderson, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(extrafont, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(scales, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(gtools, warn.conflicts = FALSE))

## Functions ----
source("./Functions/averageMaxQuant.R")

## Set fraction numbers and reporter intensity groups ----
## Match the way the data in the "Experiment" column of the MaxQuant "evidence.txt" file is labelled
fractions <- seq(1,20,1)

## WT and KD columns in MaxQuant "evidence.txt" file
wt.cols <- "^Reporter.intensity.corrected\\\\.[1,2,3,4,5]$"
kd.cols <- "^Reporter.intensity.corrected.[6,7,8,9]|Reporter.intensity.corrected.10$"

## Set folders ----

## Where is MaxQuant "evidence.txt" stored?
mq.dir <- "./Data/"

## Where are databases for metadata stored?
db.dir <- "./Databases/"

## Where do results go?
output.res.dir <- "./Output_R/"

### Set file names ----

## All TMT results
results.file <- "Phospho-TMT_results_complete.csv"

## Cleaned up TMT results
results.clean.file <- "Phospho-TMT_results_clean.csv"

### Plots and images

## VSN fit
vsn.fit.file <- "MeanSDplot_VSN_fit.png"

## QQ-plot
qqplot.file <- "QQ-plot.png"

## VSN calibration check
VSN.transform.TMT.file <- "VSN_calibration.png"

## Garden Sprinkler plot
gs.file <- "Garden_sprinkler_plot.png"

## Unlabelled Garden Sprinkler plot
gs.file.2 <- "Garden_sprinkler_plot_unlabelled.png"

## Descriptive statistics

stats.file <- "Descriptive_Statistics.csv"

## Set ggplot theme and palette ----

theme_set(theme_classic(base_family = "Arial"))
pal <- wes_palette("Darjeeling1")

# # # # # # # # # # # # # # # # # # # # # # # 

## Read in data ----
message(paste0("\\n\\nLoading ", mq.dir, "evidence.txt"))

data.mq <- read.csv(
  paste0(mq.dir, "evidence.txt"),
  sep = "\\t",
  stringsAsFactors = FALSE
)


## Remove Reverse hits, potential contaminants, PIF < 0.65 and PEP >= 0.05

rem.rev <- which(data.mq$Reverse == "+")

## MaxQuant kept a "REV__" hit as it has no Reverse flag?
rem.rev.2 <- grep("REV__", data.mq$Leading.razor.protein)
rem.rev <- union(rem.rev, rem.rev.2)

rem.cont <- which(data.mq$Potential.contaminant == "+")
rem.pep <- which(data.mq$PEP >= 0.05)
rem.pif <- which(data.mq$PIF < 0.60)
rem.pif.2 <- which(is.nan(data.mq$PIF) == TRUE)

data.preclean <- data.mq[-c(rem.rev, rem.cont, rem.pep, rem.pif, rem.pif.2), ]

## Count zero intensities for each quantified peptide, remove KD or WT rows with less than 2 datapoints

group.wt <- grep(
  wt.cols,
  colnames(data.preclean)
)

group.kd <- grep(
  kd.cols,
  colnames(data.preclean)
)

zero.kd <- rowSums(data.preclean[, group.kd] == 0)
zero.wt <- rowSums(data.preclean[, group.wt] == 0)

## Remove rows with >= 4 counts of 0 (= at least 2 datapoints in KD/WT group)

rem.zero.kd <- which(zero.kd >= 4)
rem.zero.wt <- which(zero.wt >= 4)

data.clean <- data.preclean[-c(rem.zero.kd, rem.zero.wt), ]

## Replace zeroes with NA in reporter intensities columns (MaxQuant puts zeroes instead NAs)

data.clean[, c(group.wt, group.kd)][data.clean[, c(group.wt, group.kd)] == 0] <- NA

### Data averaging ----

## Peptides with one observatiton are kept
## Peptides with two observations get averaged
## Peptide with three or more observation, peptides are averaged or peptide with median ratio WT/KD is taken
## This is being done for each fraction independently to avoid bias by fraction dependent co-isolation populations 

message("Averaging data")

data.av <- averageMaxQuant(data.clean, fractions)

## Data analysis ----
## VSN tansformation ----

message("\\nVSN Transform")

data.av.matrix <- as.matrix(data.av[ , c(group.wt, group.kd)])

## colMeans function used in averageMaxQuant returns NaN
## if all observations of a peptide in a TMT channel are NA.
## Replace with NA, VSN cannot handle otherwise
## Also is.nan() does not work in dataframes...

data.av.matrix[is.nan(data.av.matrix)] <- NA 

vsn.model <- vsn2(data.av.matrix)

## Check model fit

message(paste0("Writing ", output.res.dir, vsn.fit.file))

png(
  paste0(output.res.dir, vsn.fit.file),
  width = 6,
  height = 6,
  unit = "in",
  res = 600
)

vsn.plot <- meanSdPlot(vsn.model)

vsn.plot$gg +
  ggtitle("VSN fit of TMT data")

dev.off()

## VSN transform TMT data

data.trans <- predict(vsn.model, newdata = data.av.matrix)

## QQ plot to check normality

message(paste0("Writing ", output.res.dir, qqplot.file))

png(
  paste0(output.res.dir, qqplot.file),
  width = 6,
  height = 6,
  unit = "in",
  res = 600
)

par(mfrow = c(1,2))

qqnorm(
  data.av.matrix,
  main = "Normal Q-Q- Plot \\n TMT raw data"
)
qqline(data.av.matrix)

qqnorm(
  data.trans,
  main = "Normal Q-Q- Plot \\n TMT VSN transformed data"
)
qqline(data.trans)

dev.off()

## Boxplot VSN intensities to check intensity calibration

message(paste0("Writing ", output.res.dir, VSN.transform.TMT.file))

data.box <- melt(data.trans, id = colnames(data.trans))
colnames(data.box) <- c("Index", "TMT.label", "VSN.Intensity")

data.box$Group <- data.box$TMT.label
data.box$Group <- sub(wt.cols, "WT", data.box$Group)
data.box$Group <- sub(kd.cols, "KD", data.box$Group)
data.box$Group <- factor(data.box$Group, levels = c("WT", "KD"))

data.box$TMT.label <- sub("Reporter.intensity.corrected.", "", data.box$TMT.label)

tmt.labels <- c(
  "126",
  "127N",
  "127C",
  "128N",
  "128C",
  "129N",
  "129C",
  "130N",
  "130C",
  "131"
)

reporter <- seq(1,10,1)

for (cur.tmt in 1:length(tmt.labels)) {

data.box$TMT.label <- sub(paste0("^", reporter[cur.tmt], "$"), tmt.labels[cur.tmt], data.box$TMT.label)

}

data.box$TMT.label <- factor(data.box$TMT.label, levels = tmt.labels)

ggplot(data.box, aes(x = TMT.label, y = VSN.Intensity, fill = Group)) +
  geom_boxplot(
    lwd = 0.8,
    na.rm = TRUE
  ) +
  scale_fill_manual(values = c(pal[1], pal[5])) +
  theme(
    axis.text.x = element_text(family = "Arial", size = 14),
    axis.title = element_text(size = 16),
    axis.text.y = element_text(family = "Roboto Mono", size = 14),
    axis.line = element_line(
		  color = "black", 
      size = 1,
		  linetype = "solid"
	  ),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(.20, "cm")
  ) +
  labs(
    x = "TMT label",
    y = "VSN transformed Reporter.intensity.corrected"
  )

ggsave(
  paste0(output.res.dir, VSN.transform.TMT.file),
  dpi = 600
)

## Limma ----
## Build dataframe for limma
## KD as base, goes therefore first in data.frame

message("Limma calculation")

data.limma <- data.frame(
  Sequence = data.av$Sequence,
  Modified.sequence = data.av$Modified.sequence,
  Phospho.STY.probabilities = data.av$Phospho..STY..Probabilities,
  Phospho.STY = data.av$Phospho..STY.,
  Leading.razor.protein = data.av$Leading.razor.protein,
  Gene.names =  data.av$Gene.names,
  Fraction =  data.av$Experiment,
  data.trans[, c(6:10, 1:5)],
  stringsAsFactors = FALSE
)

## Remove non-phospho peptides

rem.non.phos <- which(data.limma[, "Phospho.STY"] == 0)

data.limma <- data.limma[-rem.non.phos, ]

## Implement counter to be sure limma output can be correctly merged with MaxQuant data

data.limma$Counter <- seq(1, nrow(data.limma), 1)
rownames(data.limma) <- data.limma$Counter

## Get Reporter Intensity columns

int.limma <- grep("Reporter", colnames(data.limma))

## Build model matrix

groups <- c(rep("KD",5), rep("WT",5)) 
groups <- factor(groups, levels = c("KD", "WT"))

design <- model.matrix(~ groups) 

## Limma
## eBayes with robust = TRUE to protect from hypervariable observations

model.limma <- lmFit(data.limma[, int.limma], design) 
fit <- eBayes(model.limma, robust = TRUE) 

TMT.results <- topTable(fit, number = nrow(data.limma))

## Combine limma dataframe with TMT data

TMT.results <- TMT.results[order(as.numeric(rownames(TMT.results))), ]
TMT.results[, "Counter"] <- rownames(TMT.results)

TMT.results <- merge(data.limma, TMT.results, by = "Counter")

## Add Metadata to each peptide ----
## Read in databases

message(paste0("Loading databases for metadata from ", db.dir))

uniprot.go <- read.csv(
  paste0(db.dir, "200405_Uniprot_human_GO.tab"),
  sep = "\\t",
  stringsAsFactors = FALSE
)

uniprot.phos <- read.fasta(
  file = paste0(db.dir, "180405_Uniprot_Homo+sapiens_canonical+isoforms.fasta"),
  seqtype = "AA",
  strip.desc = TRUE,
  as.string = TRUE
)

## Extract protein annotations , AA sequence and combine
annot <- lapply(uniprot.phos, attributes)
annot <- unlist(annot)
annot <- as.vector(annot[seq(2, length(annot), 3)])

uniprot.phos <- as.data.frame(ldply(uniprot.phos, rbind))
colnames(uniprot.phos) <- c("Protein", "AA.seq")
uniprot.phos$AA.seq <- as.character(uniprot.phos$AA.seq)

uniprot.phos <- cbind(uniprot.phos, annot)
uniprot.phos$annot <- as.character(uniprot.phos$annot)

## Extract Gene names
uniprot.phos$Gene.names <- str_extract(uniprot.phos$annot, "GN\\\\=.*? |GN\\\\=.*?$")
uniprot.phos$Gene.names <- gsub("GN=| ", "", uniprot.phos$Gene.names)

## Extract Uniprot Acc. Number

uniprot.phos$Uniprot.Acc <- str_extract(uniprot.phos$annot, "\\\\|.*?\\\\|")
uniprot.phos$Uniprot.Acc <- gsub("\\\\|", "", uniprot.phos$Uniprot.Acc)

## Extract Protein name

uniprot.phos$protein.name.long <- str_extract(uniprot.phos$annot, " .*? OS\\\\=")
uniprot.phos$protein.name.long <- gsub("^ | OS\\\\=", "", uniprot.phos$protein.name.long)

## Check if peptide is isoform specific, if not use canonical form for p-site number.
## Absolute position of phospho-peptide in AA seq
## MaxQuant decides for some reason that peptide is part of specific isoform of protein even when found in
## canonical form. We want to have position of p-site in canonical form.

TMT.results$Uniprot.Acc <- TMT.results$Leading.razor.protein
TMT.results$Uniprot.Acc <- gsub("-\\\\d{1,2}", "", TMT.results$Uniprot.Acc)

gene.cols <- which(colnames(uniprot.phos) == "Gene.names")
uni.merge <- uniprot.phos[, -gene.cols]

TMT.results <- merge(TMT.results, uni.merge, by = "Uniprot.Acc")

pep.pos <- str_locate(TMT.results$AA.seq, TMT.results$Sequence)
TMT.results[, "Pep.pos.in.AA.seq"] <- pep.pos[, "start"]

iso.pep <- which(is.na(TMT.results$Pep.pos.in.AA.seq) == TRUE)

TMT.results$Isoform.specific.peptide <- FALSE
TMT.results[iso.pep, "Isoform.specific.peptide"] <- TRUE

## Put isoform specific Uniprot Accession number, protein name and AA.seq where iso.pep indicates

TMT.results[iso.pep, "Uniprot.Acc"] <- TMT.results[iso.pep, "Leading.razor.protein"]

for (cur.iso.pep in iso.pep) {
  
  temp.row.uni <- which(TMT.results[cur.iso.pep, "Uniprot.Acc"] == uniprot.phos$Uniprot.Acc)
  TMT.results[cur.iso.pep, "AA.seq"] <- uniprot.phos[temp.row.uni, "AA.seq"]
  TMT.results[cur.iso.pep, "protein.name.long"] <- uniprot.phos[temp.row.uni, "protein.name.long"]
  
}

pep.pos <- str_locate(TMT.results$AA.seq, TMT.results$Sequence)
TMT.results[, "Pep.pos.in.AA.seq"] <- pep.pos[, "start"]

## Parent Entry for easier matching
TMT.results$Parent.Entry <- TMT.results$Leading.razor.protein
TMT.results$Parent.Entry <- gsub("-\\\\d{1,2}", "", TMT.results$Parent.Entry)

## Extract phospho-peptide sequence without indication of other modifications
TMT.results[, "Modified.sequence.2"] <- gsub("\\\\(ox\\\\)", "", TMT.results[, "Modified.sequence"])
TMT.results[, "Modified.sequence.2"] <- gsub("\\\\(ac\\\\)", "", TMT.results[, "Modified.sequence.2"])
TMT.results[, "Modified.sequence.2"] <- gsub("\\\\(de\\\\)", "", TMT.results[, "Modified.sequence.2"])
TMT.results[, "Modified.sequence.2"] <- gsub("\\\\(ca\\\\)", "", TMT.results[, "Modified.sequence.2"])
TMT.results[, "Modified.sequence.2"] <- gsub("_", "", TMT.results[, "Modified.sequence.2"])

## Add rest of metadata: 
## Position of p-site and is phosphorylation site part of RPX-S/T motif?

message("Adding metadata")

## Extract phosphorylation probabilities for each potential site within peptides of interest

TMT.results.meta <- data.frame()

progress.bar.2 <- txtProgressBar(min = 0, max = nrow(TMT.results), style = 3)

for(cur.row in 1:nrow(TMT.results)){
  
  setTxtProgressBar(progress.bar.2, cur.row)
  
  ## Several sequences might be averaged
  
  splits <- unlist(str_split(TMT.results$Phospho.STY.probabilities[cur.row], ";"))
  
  real.pos.all <- vector()
  
  ## Extract position of every S,T or Y in peptide which has a probability attached
  
  for (cur.seq in 1:length(splits)){
    
    temp.seq <- as.data.frame(str_locate_all(splits[cur.seq], "[S,T,Y]\\\\(.+?\\\\)"))
    temp.seq$Diff <- temp.seq$end - temp.seq$start
    temp.seq$Acc.Diff <- cumsum(temp.seq$Diff)
    
    real.pos.first <- temp.seq[1, "start"]
    
    real.pos.temp.2 <- vector()
    
    for (i in 1:nrow(temp.seq)){
      
      if(is.na(temp.seq$start[i+1]) == FALSE){
        
        ## Current position - additional non AA characters
        
        real.pos.temp <- temp.seq$start[i+1] - temp.seq$Acc.Diff[i]
        
        real.pos.temp.2 <- c(real.pos.temp.2, real.pos.temp)
        
      }
      
    }
    
    real.pos.temp.3 <- c(real.pos.first, real.pos.temp.2)
    
    real.pos.all <- c(real.pos.all, real.pos.temp.3)
    
  }
  
  ## Build dataframe with probability of p-site and localisation within peptide
  
  probs <- as.data.frame(str_extract_all(TMT.results$Phospho.STY.probabilities[cur.row], "[S,T,Y]\\\\(.+?\\\\)"))
  
  if (length(real.pos.all) != nrow(probs)) {
    
    message(paste0("In row", cur.row, " phosphosites do not match probabilities extracted"))
    
  }
  
  probs$Real.Pos <- real.pos.all
  
  ## Where are these sites within the AA seq?
  
  probs$Pep.pos.in.AA.seq <- TMT.results[cur.row, "Pep.pos.in.AA.seq"]
  probs$P.Site.in.AA.seq <- probs$Pep.pos.in.AA.seq + probs$Real.Pos - 1
  
  colnames(probs)[1] <- "Site.Probability"
  rem.1 <- grep("Pep.pos.in.AA.seq", colnames(probs))
  probs <- probs[, -rem.1]
  
  probs.final <- cbind(TMT.results[rep(cur.row, nrow(probs)), ], probs)
  
  TMT.results.meta <- rbind(TMT.results.meta, probs.final)
  
}

## Nicer p-site naming

TMT.results.meta$P.Site.harm <- paste0(
  substring(TMT.results.meta[, "Site.Probability"], 1, 1),
  TMT.results.meta$P.Site.in.AA.seq,
  "-p"
)

progress.bar.3 <- txtProgressBar(min = 0, max = nrow(TMT.results.meta), style = 3)

message("\\nSequence windows")

for (cur.row in 1:nrow(TMT.results.meta)) {
  
  setTxtProgressBar(progress.bar.3, cur.row)
  
  ## Lets get the sequence window
  
  seq.window <- vector()
  
  cur.phos.site <- TMT.results.meta[cur.row, "P.Site.in.AA.seq"]
  
  start.seq <- cur.phos.site - 7
  end.seq   <- cur.phos.site + 7
  
  prefix.counter <- 0
  suffix.counter <- 0
  
  ## Does the window fit? Add "_" to fill up (for RPX[S,T][A,G,P,S], works also for motif-x)
  
  if (start.seq < 1) {
    
    prefix.counter <- abs(start.seq) + 1
    start.seq <- 1
    
  }
  
  if (end.seq > nchar(TMT.results.meta[cur.row, "AA.seq"])) {
    
    suffix.counter <- end.seq - nchar(unique(TMT.results.meta[cur.row, "AA.seq"]))
    end.seq <- nchar(unique(TMT.results.meta[cur.row, "AA.seq"]))
    
  }
  
  seq.window.temp <- substring(TMT.results.meta[cur.row, "AA.seq"], start.seq, end.seq)
  
  prefix <- paste0(rep("_", prefix.counter), collapse = "")
  suffix <- paste0(rep("_", suffix.counter), collapse = "")
  
  seq.window.temp <- paste0(prefix, seq.window.temp, suffix)
  
  seq.window <- c(seq.window, seq.window.temp)
  
  ## Put into dataframe
  
  TMT.results.meta[cur.row, "Sequence.Window"] <- seq.window
  
}

## Add GO terms

col.uni.go <- c(
  "Entry",                              
  "Gene.names",                        
  "Gene.names...synonym..",             
  "Gene.ontology..biological.process.",
  "Gene.ontology.IDs"
)

uniprot.go.clean <- uniprot.go[, col.uni.go]
colnames(uniprot.go.clean)[1] <- "Parent.Entry"
colnames(uniprot.go.clean)[2] <- "Gene.Names.all"

TMT.results.final <- merge(TMT.results.meta, uniprot.go.clean, by = "Parent.Entry", all.x = TRUE)

## Adjust gene/protein names to synonyms to harmonize names in manuscript

TMT.results.final$Gene.names <- sub("YLPM1", "ZAP3", TMT.results.final$Gene.names)
TMT.results.final$Gene.names <- sub("MPLKIP", "TTDN1", TMT.results.final$Gene.names)
TMT.results.final$Gene.names <- sub("TCEB3", "ELOA", TMT.results.final$Gene.names)
TMT.results.final$Gene.names <- sub("MAPRE2", "EB2", TMT.results.final$Gene.names)


TMT.results.final$protein.name.long <- sub(
  "YLP motif-containing protein 1",
  "Nuclear protein ZAP 3",
  TMT.results.final$protein.name.long
)

TMT.results.final$protein.name.long <- sub(
  "M-phase-specific PLK1-interacting protein",
  "Tricothiodystrophy non-photosensitive 1 protein",
  TMT.results.final$protein.name.long
)

TMT.results.final$Gene.names...synonym.. <- gsub("^RP1$", "RP1 EB2", TMT.results.final$Gene.names...synonym..)
TMT.results.final$Gene.names...synonym.. <- gsub("^C14orf170 ZAP3$", "YLPM1 C14orf170", TMT.results.final$Gene.names...synonym..)
TMT.results.final$Gene.names...synonym.. <- gsub("^C7orf11 TTDN1$", "MPLKIP C7orf11", TMT.results.final$Gene.names...synonym..)

## Fold change 

TMT.results.final$Fold.change <- ifelse(
  TMT.results.final$logFC > 0,
  2^(TMT.results.final$logFC),
  as.numeric(paste0("-", 2^(abs(TMT.results.final$logFC))))
)

## Write TMT results file ----
## Remove AA sequence, increases file size unnecessarily

message(paste0("\\nWriting ", output.res.dir, results.file))

aa.seq.col <- grep("^AA.seq", colnames(TMT.results.final))

write.csv(
  TMT.results.final[, -aa.seq.col], 
  paste0(output.res.dir, results.file),
  row.names = FALSE
)

## Clean up TMT results ----
## Filter by adj.p.value <= 0.65 and sort by fold-change
## Flag 1% p-site localization FDR (>= 0.994)
## Flag >= 0.75 p-site probability

above.994 <- grep("[S,T,Y]\\\\(0\\\\.99[4-9]\\\\)|[S,T,Y]\\\\(1\\\\)", TMT.results.final$Site.Probability)
TMT.results.final$p.site.FDR.1.perc <- FALSE
TMT.results.final[above.994, "p.site.FDR.1.perc"] <- TRUE

above.075 <- grep("[S,T,Y]\\\\(0\\\\.7[5-9].*?\\\\)|[S,T,Y]\\\\(0\\\\.[8-9].*?\\\\)|[S,T,Y]\\\\(1\\\\)", TMT.results.final$Site.Probability)
TMT.results.final$Prob.75 <- FALSE
TMT.results.final[above.075, "Prob.75"] <- TRUE

## Phosphorylation site part of RPX[S,T][A,G,P,S] motif?

TMT.results.final$Part.of.RPX.ST.AGPS.motif.Prob.075 <- FALSE
TMT.results.final$Part.of.RPX.ST.AGPS.motif.1.FDR <- FALSE

rpxst <- grep("^.{4}RP[A-Z]{1}[S,T][A,G,P,S]", TMT.results.final[, "Sequence.Window"])
rpxst.75 <- which(TMT.results.final$Prob.75[rpxst] == TRUE )
rpxst.1fdr <- which(TMT.results.final$p.site.FDR.1.perc[rpxst] == TRUE) 

TMT.results.final[rpxst[rpxst.75], "Part.of.RPX.ST.AGPS.motif.Prob.075"] <- TRUE
TMT.results.final[rpxst[rpxst.1fdr], "Part.of.RPX.ST.AGPS.motif.1.FDR"] <- TRUE

TMT.results.clean <- TMT.results.final[which(TMT.results.final$adj.P.Val <= 0.65), ]
TMT.results.clean <- TMT.results.clean[order(TMT.results.clean$Fold.change, decreasing = TRUE), ]

## Write cleaned up TMT results ----
## Remove AA sequence, increases file size unnecessarily

message(paste0("Writing ", output.res.dir, results.clean.file))

aa.seq.col <- grep("^AA.seq", colnames(TMT.results.clean))

write.csv(
  TMT.results.clean[, -aa.seq.col],
  paste0(output.res.dir, results.clean.file),
  row.names = FALSE
)

## Garden Sprinkler Plot

message(paste0("Writing ", output.res.dir, gs.file))

## Keep peptides as detected by MS

TMT.results.gp <- TMT.results.final[!duplicated(TMT.results.final$Counter), ]

## Adjust y-axis scaling for better visualization
TMT.results.gp$gp.plot <- sqrt(1 - TMT.results.gp$adj.P.Val)

## Split away uninteresting peptides
TMT.results.base <- TMT.results.gp[which(TMT.results.gp$gp.plot <= 0.5), ]

## Split away peptides of interest
TMT.results.poi <- TMT.results.gp[which(TMT.results.gp$gp.plot >= 0.5), ]

## Split away peptides of interest with positive fold-change (i.e. more abundant in WT)
## and which were validated by XIC
TMT.results.poi.pos <- TMT.results.poi[which(TMT.results.poi$Fold.change > 0), ]

## Name of proteins validated by XIC

val.pep.names <- c(
  "CDKL5",
  "ZAP3",
  "EP400",
  "TTDN1",
  "ELOA",
  "ZNF592"
)

val.pep <- which(TMT.results.poi.pos$Gene.names %in% val.pep.names)

TMT.results.poi.pos <- TMT.results.poi.pos[val.pep, ]

##Sprinkler plot ----

pal.2 <- rev(wes_palette("Darjeeling1"))
pal.2 <- pal.2[c(1,4,3,2,5)]

ggplot(
  TMT.results.gp,
  aes(x = logFC, y = gp.plot)
  ) +
  geom_pointdensity() +
  scale_color_gradientn(
    colours = pal.2,
    name = "Number of \\nneighbouring \\npeptides",
    label = comma
  ) + 
  geom_point(
    data = TMT.results.poi.pos,
    colour = "black",
    size = 2
  ) +
  geom_text_repel(
    data = TMT.results.poi.pos,
    aes(label = Gene.names),
    min.segment.length = 0, 
    family = "Arial",
    fontface = "bold",
    size = 2
  ) +
  ylab(expression(sqrt(1-adj.~italic(p)-Value))) +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(log[2]~Fold~Change)) +
  theme(
    axis.title = element_text(size = 12),
    axis.text.x = element_text(family = "Arial", size = 10),
    axis.text.y = element_text(family = "Arial", size = 10),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = c(0.15, 0.77),
  )

ggsave(
  paste0(output.res.dir, gs.file),
  width = 109,
  height = 109,
  unit = "mm",
  dpi = 600
)

## Sprinkler plot without labels

ggplot(
  TMT.results.gp,
  aes(x = logFC, y = gp.plot)
) +
  geom_pointdensity() +
  scale_color_gradientn(
    colours = pal.2,
    name = "Number of \\nneighbouring \\npeptides",
    label = comma
  ) + 
  geom_point(
    data = TMT.results.poi.pos,
    colour = "black",
    size = 2
  ) +
  ylab("") +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(log[2]~Fold~Change)) +
  theme(
    axis.title = element_text(size = 12),
    axis.text.x = element_text(family = "Arial", size = 10),
    axis.text.y = element_text(family = "Arial", size = 10),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = c(0.15, 0.77),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

ggsave(
  paste0(output.res.dir, gs.file.2),
  width = 109,
  height = 109,
  unit = "mm",
  dpi = 1200,
  bg = "transparent"
)

## Table for Figure ----

cols.fig <- c(
  "Counter",
  "Gene.names",
  "protein.name.long",
  "Uniprot.Acc",
  "Fold.change",
  "P.Site.harm",
  "Phospho.STY",
  "Phospho.STY.probabilities",
  "Site.Probability",
  "Modified.sequence.2",
  "Sequence.Window",
  "Part.of.RPX.ST.AGPS.motif.Prob.075",
  "p.site.FDR.1.perc",
  "Prob.75"
)

table.fig <- TMT.results.clean[which(TMT.results.clean$Fold.change > 0), cols.fig]
table.fig$Fold.change <- round(table.fig$Fold.change, 2)
table.fig <- table.fig[order(-table.fig$Fold.change, table.fig$Counter, table.fig$Site.Probability), ]


## p-site >= 0.75 only and unique sites
table.fig <- table.fig[which(table.fig$Prob.75 == TRUE), ]
table.fig <- table.fig[!duplicated(table.fig$Counter,  fromLast = TRUE), ]
table.fig <- table.fig[!duplicated(table.fig$Sequence.Window), ]


## Extract phosphorylation motif

table.fig$motif <- substring(
  table.fig$Sequence.Window,
  5,
  9
)

## Clean-up and write table

table.fig.clean <- table.fig[table.fig$Phospho.STY == 1, ]

cols <- c(
  "Gene.names",
  "protein.name.long",
  "Uniprot.Acc",
  "Fold.change",
  "P.Site.harm",
  "Site.Probability",
  "Part.of.RPX.ST.AGPS.motif.Prob.075",
  "p.site.FDR.1.perc",
  "Prob.75",
  "motif"
)

table.fig.clean <- table.fig.clean[, cols]

colnames(table.fig.clean) <- c(
  "Protein",
  "",
  "Uniprot.Acc.No.",
  "Fold.Change",
  "Phosphorylation.Site(s)",
  "Best.phospho.STY.probability",
  "Part.of.RPX[S,T][A,G,P,S].motif",
  "Phospho.Site.within.1%.FDR",
  "Phospho.Site.probability.at.least.0.75",
  "Motif"
)

table.fig.clean$`Phosphorylation.Site(s)` <- sub("S", "Ser", table.fig.clean$`Phosphorylation.Site(s)`)
table.fig.clean$`Phosphorylation.Site(s)` <- sub("T", "Thr", table.fig.clean$`Phosphorylation.Site(s)`)
table.fig.clean$`Phosphorylation.Site(s)` <- sub("Y", "Tyr", table.fig.clean$`Phosphorylation.Site(s)`)

table.fig.clean$`Phosphorylation.Site(s)` <- sub("-p", "", table.fig.clean$`Phosphorylation.Site(s)`)


write.csv(
  table.fig.clean,
  paste0(output.res.dir, "Table_Figure_TMT.csv"),
  row.names = FALSE
)

## Supplementary Table

cols.2.keep <- c(
  "Counter",
  "Uniprot.Acc",
  "Gene.names",
  "Gene.names...synonym..",
  "protein.name.long",
  "Modified.sequence.2",
  "Phospho.STY.probabilities",
  "P.Site.in.AA.seq",
  "Site.Probability",
  "Prob.75",
  "p.site.FDR.1.perc",
  "Part.of.RPX.ST.AGPS.motif.Prob.075",
  "Part.of.RPX.ST.AGPS.motif.1.FDR",
  "adj.P.Val",
  "logFC",
  "Fold.change",
  "Phospho.STY",
  "Fraction",
  "Reporter.intensity.corrected.6",
  "Reporter.intensity.corrected.7",    
  "Reporter.intensity.corrected.8",
  "Reporter.intensity.corrected.9",    
  "Reporter.intensity.corrected.10",
  "Reporter.intensity.corrected.1",    
  "Reporter.intensity.corrected.2",
  "Reporter.intensity.corrected.3",    
  "Reporter.intensity.corrected.4",
  "Reporter.intensity.corrected.5",
  "t",
  "P.Value",
  "Isoform.specific.peptide"
)

sup.table <- TMT.results.final[, cols.2.keep]
sup.table <- sup.table[order(sup.table$adj.P.Val, sup.table$Uniprot.Acc, sup.table$P.Site.in.AA.seq), ]

write.csv(
  sup.table,
  paste0(output.res.dir, "Supplementary_Table_TMT.csv"),
  row.names = FALSE
)

message("Analysis done \\n")
