# The genetic basis of variation in sexual aggression: evolution versus social plasticity

# PLASTICITY - Forced copulation RNAseq DE Analysis using EdgeR
# Andrew M. Scott

#library(apeglm)
library(edgeR)
library(tximport)
library(tidyverse)
library(statmod)
library(grid)
library(gridExtra)

gene.symbols <- read_csv(file = 'fbgn_annotation_ID_fb_2020_04.csv')

# Functions for visualizing transcript counts by treatment

countPlot <- function(gene, count_file, group=group, group_color = Trt) {
  gene <- deparse(substitute(gene))
  gene.count <- count_file[gene, ]
  df <- data.frame(gene.count, group, group_color)
  gene.symbol <- gene.symbols[which(grepl(gene, gene.symbols$`primary_FBgn#`)), 1]
  plot <- ggplot(df, aes(y = gene.count, 
                         x = group, 
                         color = group_color)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(pch = 1, cex = 2.25, 
               position = position_jitter(width = 0.02, height = 0), 
               stroke = 1) +
    labs(y = "Counts or logCounts", title = paste(gene.symbol," ","(",gene,")", sep = "")) +
    theme_classic() +
    theme(legend.position = "none")
  if (length(unique(group_color))==2) {
    plot + scale_color_manual(values = c("blue", "red")) + scale_x_discrete(limits = c("E.C", "I.C", "E.T", "I.T"))
  } else {
    plot + scale_color_manual(values = c("blue", "black", "red")) + scale_x_discrete(limits = c("D.N", "C.N", "U.N", "D.T", "C.T", "U.T"))
  }
}

countPlot.eR <- function(gene, count_file, group=group, group_color = Trt) {
  gene <- deparse(substitute(gene))
  gene.count <- count_file[gene, ]
  df <- data.frame(gene.count, group, group_color)
  gene.symbol <- gene.symbols[which(grepl(gene, gene.symbols$`primary_FBgn#`)), 1]
  plot <- ggplot(df, aes(y = gene.count, 
                         x = group, 
                         color = group_color)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(pch = 1, cex = 2.25, 
               position = position_jitter(width = 0.02, height = 0), 
               stroke = 1) +
    labs(y = "Counts or logCounts", title = paste(gene.symbol," ","(",gene,")", sep = "")) +
    theme_classic() +
    theme(legend.position = "none")
  if (length(unique(group_color))==2) {
    plot + scale_color_manual(values = c("blue", "red")) + scale_x_discrete(limits = c("ConExp", "ConIso", "TenExp", "TenIso"))
  } else {
    plot + scale_color_manual(values = c("blue", "black", "red")) + scale_x_discrete(limits = c("ND", "NC", "NU", "TD", "TC", "TU"))
  }
}

# Make object of quant files with associated sample names

quant.files <- file.path("Plasticity_DE_Data/salmon_counts/quants", list.files("Plasticity_DE_Data/salmon_counts/quants"), "quant.sf")

samples <- str_match(quant.files, "-{1}(.*?)_")
samples <- samples[,2]

names(quant.files) <- samples

# Load in transcript-to-gene tsv

tx2gene <- read.table("txp_to_gene_alt.tsv", col.names = c("TXNAME", "GENEID"))

length(table(tx2gene$TXNAME))
length(table(tx2gene$GENEID))

# Read input using tximport (automatically sums results to gene level)

txi.eR <- tximport(quant.files,
                   type = "salmon",
                   tx2gene = tx2gene)

# Inspect txi

summary(txi.eR)
str(txi.eR)
head(txi.eR$counts)

##### Make DGEList object (taken from tximport vignette)

cts <- txi.eR$counts
normMat <- txi.eR$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
cts.dgelist <- DGEList(cts)
cts.dgelist <- scaleOffset(cts.dgelist, normMat)

# import design matrix
rna.design.init <- read.csv("Plasticity_DE_Data/plasticity_design.csv", header = T, stringsAsFactors = T)

# remove missing sample (ET1-1)
rna.design.init <- rna.design.init[-19,]

# make sample the rownames
rna.design.e <- rna.design.init[-1]
row.names(rna.design.e) <- rna.design.init$Sample

# add design to DGElist
samples <- merge(cts.dgelist$samples, rna.design.e, by = "row.names", all = T)
cts.dgelist$samples <- merge(cts.dgelist$samples, rna.design.e, by = "row.names", all = T)
cts.dgelist$samples <- samples[-1]
row.names(cts.dgelist$samples) <- samples$Row.names

# reorder column names in count and offest matrices to match sample df

cts.dgelist$counts <- cts.dgelist$counts[, order(colnames(cts.dgelist$counts))]
cts.dgelist$offset <- cts.dgelist$offset[, order(colnames(cts.dgelist$offset))]

# filtering out lowly expressed genes. Vignette says to do this by group you plan to look at DE of
keep <- filterByExpr(cts.dgelist, 
                     group = cts.dgelist$samples$Trt,
                     min.count = 5)

table(keep)

cts.dgelist <- cts.dgelist[keep, , keep.lib.sizes=FALSE]

# Create group variable for group contrasts

cts.dgelist$samples$group <- factor(paste0(cts.dgelist$samples$Teneral_presence,
                                           cts.dgelist$samples$Trt))

##### Plot mds

plotMDS(cts.dgelist, 
        labels = cts.dgelist$samples$Trt)

plotMDS(cts.dgelist, 
        labels = cts.dgelist$samples$Day)

plotMDS(cts.dgelist, 
        labels = cts.dgelist$samples$Teneral_presence)

##### Basic model ~Day - Estimate dispersion

design.1 <- model.matrix(~ Day, 
                         data = cts.dgelist$samples)

cts.dgelist.disp <- estimateDisp(cts.dgelist, 
                                 design.1)

plotBCV(cts.dgelist.disp)

# GLM

fit.1 <- glmQLFit(cts.dgelist.disp, 
                  design = design.1)
plotQLDisp(fit.1)

# test

day.test <- glmQLFTest(fit.1)

summary(decideTests(day.test))

topTags(day.test)

##### More complex model (trt main effect), with interaction

design.2 <- model.matrix(~ Vial_rep + Day + Teneral_presence + Trt + Teneral_presence:Trt, 
                         data = cts.dgelist$samples)

cts.dgelist.trt <- estimateDisp(cts.dgelist, 
                                design.2)
plotBCV(cts.dgelist.disp)

# GLM

fit.trt <- glmQLFit(cts.dgelist.trt, 
                    design = design.2)

plotQLDisp(fit.trt)

trt.test <- glmQLFTest(fit.trt, 
                       coef = "TrtIso")

summary(decideTests(trt.test)) 

plotMD(trt.test)

edger.trt <- topTags(trt.test, n=Inf)
edger.trt.222 <- topTags(trt.test, n = )

countPlot.eR(FBgn0031306, count_file = trt.test$fitted.values, group = trt.test$samples$group, group_color = trt.test$samples$Trt)

# interaction effect

int.test <- glmQLFTest(fit.trt, 
                       coef = "Teneral_presenceTen:TrtIso")

summary(decideTests(int.test)) # 0

plotMD(int.test)

countPlot.eR(FBgn0039714, count_file = int.test$fitted.values, group = int.test$samples$group, group_color = int.test$samples$Trt)

plast.eR.hits <- topTags(int.test, n = Inf)

##### Analysis of group contrasts

design.3 <- model.matrix(~ Vial_rep + Day + group,  
                         data = cts.dgelist$samples)

cts.dgelist.group <- estimateDisp(cts.dgelist, 
                                  design.3)
plotBCV(cts.dgelist.group)

# GLM

fit.group <- glmQLFit(cts.dgelist.group, 
                      design = design.3)
plotQLDisp(fit.group)

group.test <- glmQLFTest(fit.group)

summary(decideTests(group.test))

# Group contrasts

TIsovsTExp <- makeContrasts(groupTenIso-groupTenExp, levels = design.3)

TIsovsTExp.test <- glmQLFTest(fit.group, 
                              contrast = TIsovsTExp)

summary(decideTests(TIsovsTExp.test))  

hits.edger <- topTags(TIsovsTExp.test, n = Inf)

#### Effect of teneral exposure (main eff)

ten.test <- glmQLFTest(fit.trt, 
                       coef = "Teneral_presenceTen")

summary(decideTests(ten.test)) 

plotMD(ten.test)

edger.ten <- topTags(ten.test, n=Inf) # 0 


################################################
#### Plasticity DE analysis check with limma-voom

raw.counts <- txi.eR$counts

# Create DGEList

d0 <- DGEList(raw.counts)

# Calculate normalization factors

d0 <- calcNormFactors(d0)

# Filter low-expressed

cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

# Get design data from the column names

samples <- colnames(raw.counts)
Trt <- as.factor(substr(samples, 1, 1))
Teneral_presence <- as.factor(substr(samples, 2, 2))
Vial_rep <- as.factor(substr(samples, 3, 3))
Day <- as.factor(substr(samples, 5, 5))
group <- as.factor(interaction(Trt, Teneral_presence))

# Plot MDS

plotMDS(d, col = as.numeric(group))

# Voom transformation and calculation of variance weights

mm <- model.matrix(~0 + Day + Vial_rep + group)

y <- voom(d, mm, plot = T)

# Fitting linear models in limma

fit <- lmFit(y, mm)

# Make contrast of comp between Exp and Iso (tenerals present)

contr <- makeContrasts(groupI.T - groupE.T, levels = colnames(coef(fit)))
contr

# Estimate contrast for each gene

tmp <- contrasts.fit(fit, contr)

# Empirical Bayes smoothing of std errors

tmp <- eBayes(tmp)

# Top hits (adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value)

top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 30)

# How many hits

length(which(top.table$adj.P.Val < 0.05))

hits.limmav <- top.table[which(top.table$adj.P.Val < 0.05), ]
head(hits.limmav, 20)

# Plot counts for some genes as a test

countPlot(FBgn0031306, y$E, group)  # highly expressed in adult head

countPlot(FBgn0037835, y$E, group)  # highly expressed in adult head

countPlot(FBgn0031327, y$E, group)

countPlot(FBgn0000500, y$E, group) # Dsk - behavioural - highly expressed in adult head


# on voom estimates (log2 CPM scale) - logFC obtained from slope bw groups in contrast
countPlot(FBgn0000500, y$E, group)

### Main effect of Trt and Interaction

mm.2 <- model.matrix(~1 + Day + Vial_rep + Trt + Teneral_presence + Trt:Teneral_presence)

y.2 <- voom(d, mm.2, plot = T)

# Fitting linear models in limma

fit.2 <- lmFit(y.2, mm.2)


##### Main effect of Trt

# Empirical Bayes smoothing of std errors

tmp.2 <- eBayes(fit.2)

# Top hits (adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value)

top.table.trt <- topTable(tmp.2, coef = "TrtI", sort.by = "P", n = Inf)
head(top.table.trt, 30)

# How many hits

length(which(top.table.trt$adj.P.Val < 0.05))

hits.limmav.trt <- top.table.trt[which(top.table.trt$adj.P.Val < 0.05), ]
head(hits.limmav.trt, 20)


##### Interaction

# Top hits (adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value)

top.table.int <- topTable(tmp.2, coef = "TrtI:Teneral_presenceT", sort.by = "P", n = Inf)
head(top.table.int, 30)

# How many hits

length(which(top.table.int$adj.P.Val < 0.05))

###### Comparison of ExpC-ExpT and IsoC-IsoT contrasts

contr.IC.IT <- makeContrasts(groupI.T - groupI.C, levels = colnames(coef(fit)))
contr.IC.IT

# relevel for EC-ET contrast
group.rl <- relevel(group, ref = "I.C")
mm.EC.ET <- model.matrix(~0 + Day + Vial_rep + group.rl)
y.EC.ET <- voom(d, mm.EC.ET, plot = T)
fit.EC.ET <- lmFit(y.EC.ET, mm.EC.ET)

contr.EC.ET <- makeContrasts(group.rlE.T - group.rlE.C, levels = colnames(coef(fit.EC.ET)))
contr.EC.ET


##### IC-IT
tmp.IC.IT <- contrasts.fit(fit, contr.IC.IT)

# Empirical Bayes smoothing of std errors

tmp.IC.IT <- eBayes(tmp.IC.IT)

# Top hits (adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value)

top.table.IC.IT <- topTable(tmp.IC.IT, sort.by = "P", n = Inf)

head(top.table.IC.IT) # 0


##### EC-ET
tmp.EC.ET <- contrasts.fit(fit.EC.ET, contr.EC.ET)

# Empirical Bayes smoothing of std errors

tmp.EC.ET <- eBayes(tmp.EC.ET)

# Top hits (adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value)

top.table.EC.ET <- topTable(tmp.EC.ET, sort.by = "P", n = Inf)
head(top.table.EC.ET) # 0


#############################################################
#### Hits from limma-voom, edgeR

limma.hits <- top.table[which(top.table$adj.P.Val < 0.05), ]
limma.hits <- limma.hits %>% dplyr::rename(c(limma_logfc=logFC, 
                                             limma_AveExpr=AveExpr,
                                             limma_t=t,
                                             limma_P.value=P.Value,
                                             limma_adj.P.Val=adj.P.Val,
                                             limma_B=B))

# make df with just significant edgeR hits (according to decideTests)
edger.sig <- decideTests(TIsovsTExp.test)
hits.edger.int <- merge(hits.edger, edger.sig, by = "row.names", all.x = F)
hits.edger.final <- hits.edger.int[-1]
row.names(hits.edger.final) <- hits.edger.int$Row.names

edger.hits <- hits.edger.final[which(hits.edger.final$`-1*groupTenExp 1*groupTenIso` == 1
                                     | hits.edger.final$`-1*groupTenExp 1*groupTenIso` == -1), ]
edger.hits <- dplyr::rename(edger.hits, c("edger_logFC"="logFC",
                                          "edger_logCPM"="logCPM",
                                          "edger_F"="F",
                                          "edger_PValue"="PValue",
                                          "edger_FDR"="FDR",
                                          "edger_hit"="-1*groupTenExp 1*groupTenIso"))

both.hits <- merge(limma.hits, edger.hits, by = "row.names", all.x = F)
dim(both.hits) #149 both (edgeR had 184)
row.names(both.hits) <- both.hits$Row.names

plot(both.hits$limma_logfc, both.hits$edger_logFC)
cor.test(both.hits$limma_logfc, both.hits$edger_logFC)


#### Full models - Main effect of Trt

hits.limmav.trt <- dplyr::rename(hits.limmav.trt, c("limma_logfc"="logFC", 
                                                    "limma_AveExpr"="AveExpr",
                                                    "limma_t"="t",
                                                    "limma_P.value"="P.Value",
                                                    "limma_adj.P.Val"="adj.P.Val",
                                                    "limma_B"="B"))

edger.trt.sig <- decideTests(trt.test)
hits.edger.trt.int <- merge(edger.trt, edger.trt.sig, by = "row.names", all.x=F)
hits.edger.trt.final <- hits.edger.trt.int[-1]
row.names(hits.edger.trt.final) <- hits.edger.int$Row.names

edger.trt.hits <- hits.edger.trt.final[which(hits.edger.trt.final$TrtIso == 1
                                             | hits.edger.trt.final$TrtIso == -1), ]

edger.trt.hits <- dplyr::rename(edger.trt.hits, c("edger_logFC"="logFC",
                                                  "edger_logCPM"="logCPM",
                                                  "edger_F"="F",
                                                  "edger_PValue"="PValue",
                                                  "edger_FDR"="FDR",
                                                  "edger_hit"="TrtIso"))

both.trt.hits <- merge(hits.limmav.trt, edger.trt.hits, by = "row.names", all.x = F)
dim(both.trt.hits) # 177 both

plot(both.trt.hits$limma_logfc, both.trt.hits$edger_logFC)


#### Hits in both the overall IsoTvExpT and Trt main eff

both.both <- merge(both.hits, both.trt.hits, by = "Row.names", all.x = F)
dim(both.both) # 62 hits 



#####
#### Revision Jan 2022 Additions to analysis:

##### looking at the more complicated model without the interaction to see if it is much different

design.2.noint <- model.matrix(~ Vial_rep + Day + Teneral_presence + Trt, 
                               data = cts.dgelist$samples)

cts.dgelist.trt.noint <- estimateDisp(cts.dgelist, 
                                      design.2.noint)

# GLM

fit.trt.noint <- glmQLFit(cts.dgelist.trt.noint, 
                          design = design.2.noint)

plotQLDisp(fit.trt.noint)

trt.test.noint <- glmQLFTest(fit.trt.noint, 
                             coef = "TrtIso")

summary(decideTests(trt.test.noint)) 

plotMD(trt.test.noint)

edger.trt.noint <- topTags(trt.test.noint, n=Inf)

# correlation for ALL the genes
edger.trt.compare.int_vs_noint <- merge(edger.trt, edger.trt.noint, by = "row.names", all.x = F)
cor.test(edger.trt.compare.int_vs_noint$logFC.x, edger.trt.compare.int_vs_noint$logFC.y)

# correlation for the FOCAL genes of interest (i.e., the top ones)
edger.trt.compare.int_vs_noint_FOCAL <- merge(edger.trt.hits, edger.trt.noint, by = "row.names", all.x = F)
cor.test(edger.trt.compare.int_vs_noint_FOCAL$edger_logFC, edger.trt.compare.int_vs_noint_FOCAL$logFC) 
# Therefore: correlation between our focal (significant) trt effect estimates with and without interaction in the model: 0.979 95%CI[0.9728726, 0.9831687], p<2.2e-16

plot(edger.trt.compare.int_vs_noint_FOCAL$edger_logFC, edger.trt.compare.int_vs_noint_FOCAL$logFC) # note its just negative because its doing the opposite direction of the trt eff for one of them

#####
#### Test the top genes with glmmTMB as in the artificial selection analysis

library(glmmTMB)

tpm.plast <- txi.eR$abundance

trt.tmb <- as.factor(substr(colnames(tpm.plast),1,1))
teneral.tmb <- as.factor(substr(colnames(tpm.plast),2,2))
vial_rep.tmb <- as.factor(substr(colnames(tpm.plast), 3,3))
day.tmb <- as.factor(substr(colnames(tpm.plast), 5,5))

tmb.offset.plast <- y$targets$norm.factors

tmb.plast.genes.check <- list()

for (g in row.names(edger.trt.hits)) {
  tpm <- as.vector(tpm.plast[g, ])
  tpm[tpm == 0] <- 1e-06 # change zeroes to 1 transcript per million? (minimum nonzero in tpm.as is 3e-4)
  tmb.data <- data.frame(tpm, trt.tmb, teneral.tmb, day.tmb, vial_rep.tmb)
  tmb.mod <- glmmTMB(log2(tpm) ~ trt.tmb + 
                       teneral.tmb + 
                       trt.tmb:teneral.tmb + 
                       day.tmb + 
                       (1 | vial_rep.tmb), 
                     offset = tmb.offset.plast,
                     data = tmb.data)
  if (anyNA(as.vector(summary(tmb.mod)$coef$cond[,2:4])) == T) {
    tmb.mod.reduced <- glmmTMB(log2(tpm) ~ trt.tmb + 
                                 teneral.tmb + 
                                 trt.tmb:teneral.tmb + 
                                 day.tmb + vial_rep.tmb, 
                               offset = tmb.offset.plast,
                               data = tmb.data)
    tmb.plast.genes.check[[g]] <- list(summary(tmb.mod.reduced)$coef$cond, "Used Reduced Model")
  } else {
    tmb.plast.genes.check[[g]] <- list(summary(tmb.mod)$coef$cond, "Used Full Model")
  }
}

tmb.IvE.plast.estimates <- data.frame(FBgnID = character(), Estimate = double(), StdError = double(), P = double(), Model = character())

for (g in names(tmb.plast.genes.check)) {
  df <- data.frame(FBgnID = g, 
                   Estimate = tmb.plast.genes.check[g][[1]][[1]][2,1], 
                   StdError = tmb.plast.genes.check[g][[1]][[1]][2,2], 
                   P = tmb.plast.genes.check[g][[1]][[1]][2,4], 
                   Model = tmb.plast.genes.check[g][[1]][[2]])
  tmb.IvE.plast.estimates <- rbind(tmb.IvE.plast.estimates, df)
}

tmb.IvE.plast.estimates # All used the full model

# compare to edgeR estimates

row.names(tmb.IvE.plast.estimates) <- tmb.IvE.plast.estimates$FBgnID

plot(tmb.IvE.plast.estimates$Estimate, edger.trt.hits$edger_logFC)
cor.test(tmb.IvE.plast.estimates$Estimate, edger.trt.hits$edger_logFC)




