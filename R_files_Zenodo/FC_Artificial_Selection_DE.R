# The genetic basis of variation in sexual aggression: evolution versus social plasticity

# Artificial Selection - Forced copulation RNAseq DE Analysis using NEBULA
# Andrew M. Scott

#install.packages("devtools")
#library(devtools)
#install_github("lhe17/nebula")

library(apeglm)
library(nebula)
library(edgeR)
library(tximport)
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(grid)
library(gridExtra)

gene.symbols <- read_csv(file = 'fbgn_annotation_ID_fb_2020_04.csv')

# Functions for visualizing transcript counts by treatment and aggregating plots into a single figure

countPlot <- function(gene, count_file, group=group, group_color = Trt) {
  gene <- deparse(substitute(gene))
  gene.count <- count_file[gene, ]
  df <- data.frame(gene.count, group)
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
arrangePlots <- function(plotList, ncols = 4, LeftText = "Normalized Count", BottomText = "Group", TopText = "Gene Count Plots"){
  N <- length(plotList)
  for (p in seq_along(plotList)) {
    plotList[[p]] <- plotList[[p]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
  }
  grid.arrange(arrangeGrob(grobs = plotList, 
                           ncol = ncols, 
                           left = textGrob(LeftText, rot = 90, vjust = 1),
                           bottom = textGrob(BottomText),
                           top = textGrob(TopText, gp = gpar(fontsize = 18))))
}
multiCountPlot <- function(geneList, count_file, group, group_color) {
  listofplots <- list()
  for (g in geneList) {
    gene.count <- count_file[g, ]
    df <- data.frame(gene.count, group)
    gene.symbol <- gene.symbols[which(grepl(g, gene.symbols$`primary_FBgn#`)), 1]
    plot <- ggplot(df, aes(y = gene.count, 
                           x = group, 
                           color = group_color)) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(pch = 1, cex = 1.75, 
                 position = position_jitter(width = 0.02, height = 0), 
                 stroke = 0.75) +
      labs(y = "Counts or logCounts", title = paste(gene.symbol," ","(",g,")", sep = "")) +
      theme_classic() +
      theme(legend.position = "none")
    if (length(unique(group_color))==2) {
      plot <- plot + scale_color_manual(values = c("blue", "red")) + scale_x_discrete(limits = c("E.C", "I.C", "E.T", "I.T"))
    } else {
      plot <- plot + scale_color_manual(values = c("blue", "black", "red")) + scale_x_discrete(limits = c("D.N", "C.N", "U.N", "D.T", "C.T", "U.T"))
    }
    listofplots[[g]] <- plot
  }
  return(listofplots)
}

multiCountPlot_extended_in_prog <- function(geneList, count_file, group, group_color) {
  listofplots <- list()
  for (g in geneList) {
    gene.count <- count_file[g, ]
    if (anyNA(group) == T) {
      df.1 <- data.frame(gene.count, group_color)
    } else {
      df.2 <- data.frame(gene.count, group, group_color)
    }
    gene.symbol <- gene.symbols[which(grepl(g, gene.symbols$`primary_FBgn#`)), 1]
    if (anyNA(group) == T) {
      plot <- ggplot(df.1, aes(y = gene.count,
                               x = group_color,
                               color = group_color)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(pch = 1, cex = 1.75, 
                   position = position_jitter(width = 0.02, height = 0), 
                   stroke = 0.75) +
        labs(y = "Counts or logCounts", title = paste(gene.symbol," ","(",g,")", sep = "")) +
        theme_classic() +
        theme(legend.position = "none")
      if (length(unique(group_color))==2) {
        plot <- plot + scale_color_manual(values = c("blue", "red")) + scale_x_discrete(limits = c("Experienced", "Isolated"))
      } else {
        plot <- plot + scale_color_manual(values = c("blue", "black", "red")) + scale_x_discrete(limits = c("Down", "Control", "Up"))
      }
    } else {
      plot <- ggplot(df.2, aes(y = gene.count, 
                               x = group, 
                               color = group_color)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(pch = 1, cex = 1.75, 
                   position = position_jitter(width = 0.02, height = 0), 
                   stroke = 0.75) +
        labs(y = "Counts or logCounts", title = paste(gene.symbol," ","(",g,")", sep = "")) +
        theme_classic() +
        theme(legend.position = "none")
      if (length(unique(group_color))==2) {
        plot <- plot + scale_color_manual(values = c("blue", "red")) + scale_x_discrete(limits = c("E.C", "I.C", "E.T", "I.T"))
      } else {
        plot <- plot + scale_color_manual(values = c("blue", "black", "red")) + scale_x_discrete(limits = c("D.N", "C.N", "U.N", "D.T", "C.T", "U.T"))
      }
    }
    listofplots[[g]] <- plot
  }
  return(listofplots)
}

# Make object of quant files with associated sample names

quant.files.as <- file.path("AS_DE_Data/salmon_counts/AS", list.files("AS_DE_Data/salmon_counts/AS"), "quant.sf")

samples.as <- str_match(quant.files.as, "-{1}(.*?)_")
samples.as <- samples.as[,2]

names(quant.files.as) <- samples.as

# Load in transcript-to-gene tsv

tx2gene.as <- read.table("txp_to_gene_alt.tsv", col.names = c("TXNAME", "GENEID"))

length(table(tx2gene.as$TXNAME))
length(table(tx2gene.as$GENEID))

#### Read in as with edgeR, and make an offset object for use with Nebula

txi.eR.as <- tximport(quant.files.as,
                      type = "salmon",
                      tx2gene = tx2gene.as)

cts.as <- txi.eR.as$counts
normMat.as <- txi.eR.as$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat.as <- normMat.as/exp(rowMeans(log(normMat.as)))
normCts.as <- cts.as/normMat.as

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib.as <- calcNormFactors(normCts.as) * colSums(normCts.as)

# Combining effective library sizes with the length factors, and calculating
# offsets for use in Nebula
normMat.as <- sweep(normMat.as, 2, eff.lib.as, "*")
normMat.as <- log(normMat.as)

# Extracting variables for fixed, random effects

samples.neb <- colnames(cts.as)
lineage.neb <- as.factor(substr(samples.neb, 1, 2))
trt.neb <- as.factor(substr(samples.neb, 1, 1))
teneral.neb <- as.factor(substr(samples.neb, 3, 3))
day.neb <- as.factor(substr(samples.neb, 2, 2))
vial_rep.neb <- as.factor(substr(samples.neb, 4, 4))

fixedeff.neb <- data.frame(samples.neb, trt.neb, teneral.neb, day.neb)

# Generate offset for nebula object
offset.neb <- as.vector(colMeans(normMat.as))

# Make neb data object
neb.as.data <- group_cell(count = cts.as, 
                          id = lineage.neb, # Random effect 
                          pred = fixedeff.neb, # Fixed effects
                          offset = offset.neb)

# Make D baseline so we can see the U-D contrast
neb.as.data$pred <- as.data.frame(neb.as.data$pred)
neb.as.data$pred$trt.neb <- relevel(as.factor(neb.as.data$pred$trt.neb), ref = "D")

# Nebula model
neb.mat.1 <- model.matrix(~ trt.neb + teneral.neb + trt.neb:teneral.neb, # Don't have day because it is just a combination of lineage and Trt
                          data = as.data.frame(neb.as.data$pred))

# TAKES ~10min TO RUN (removes genes with low expression at this stage also)
#neb.mod.1 <- nebula(count = neb.as.data$count, 
#                    id = neb.as.data$id, 
#                    pred = neb.mat.1, 
#                    offset = neb.as.data$offset, 
#                    method = 'HL', 
#                    cpc = 5) # removes genes with counts per "cell"  less than 5

# save and load object so don't have to re-run
#save(neb.mod.1, file = "nebmod1.Rdata")
neb.mod.1 <- get(load(file = "AS_DE_Data/nebmod1.Rdata"))

rownames(neb.mod.1$summary) <- neb.mod.1$summary$gene

neb.mod.1$summary <- neb.mod.1$summary[order(neb.mod.1$summary$p_trt.nebU, decreasing = F),]
head(neb.mod.1$summary, n = 10)

length(which(neb.mod.1$summary$p_trt.nebU < 0.05))
neb.as.hits <- head(neb.mod.1$summary, n = 919)


# Quick check to see if it is different without the interaction
neb.mat.noint <- model.matrix(~ trt.neb + teneral.neb, 
                              data = as.data.frame(neb.as.data$pred))

#neb.mod.noint <- nebula(count = neb.as.data$count, 
#                    id = neb.as.data$id, 
#                    pred = neb.mat.noint, 
#                    offset = neb.as.data$offset, 
#                    method = 'HL', 
#                    cpc = 5) # removes genes with counts per "cell"  less than 5

neb.mod.noint <- get(load(file = "AS_DE_Data/nebmodnoint.Rdata"))
rownames(neb.mod.noint$summary) <- neb.mod.noint$summary$gene

neb.as.hits.comp_to_noint <- merge(neb.mod.noint$summary, neb.as.hits, by = "row.names", all.x = F)
plot(neb.as.hits.comp_to_noint$logFC_trt.neb1, neb.as.hits.comp_to_noint$logFC_trt.nebU)
cor.test(neb.as.hits.comp_to_noint$logFC_trt.neb1, neb.as.hits.comp_to_noint$logFC_trt.nebU)


# Also looking at Teneral presence main effect...
neb.mod.1.ten <- neb.mod.1$summary[order(neb.mod.1$summary$p_teneral.nebT, decreasing = F),]
head(neb.mod.1.ten, n = 10)

length(which(neb.mod.1.ten$p_teneral.nebT < 0.05))
neb.as.ten.hits <- head(neb.mod.1.ten, n = 260)
# write.csv(neb.as.ten.hits, file = "AS_ten_hits.csv")


########### Shrink estimates with apeglm
# Note - log(2) scaling to convert log 2 (fold change) to natural log scale, which apeglm uses
# https://www.bioconductor.org/packages/devel/bioc/vignettes/apeglm/inst/doc/apeglm.html

neb.mat.alpha <- model.matrix(~ trt.neb + teneral.neb + trt.neb:teneral.neb,
                              data = as.data.frame(neb.as.data$pred)[order(neb.as.data$pred$samples.neb), ])

# estimates we will be shrinking:

apeglm.input.est <- neb.mod.1$summary[, c("logFC_trt.nebU", "se_trt.nebU")]
apeglm.input.est <- apeglm.input.est[order(row.names(apeglm.input.est)), ]
apeglm.input.est <- apeglm.input.est[complete.cases(apeglm.input.est), ] # there was one NA val in the SEs

# also rm an extreme outlier
#countPlot(FBgn0053870, tpm.as, group = group.as, Trt.as)
apeglm.input.est <- apeglm.input.est[-which(apeglm.input.est == min(apeglm.input.est)), ]

# only looking at those genes that were not cut in the above NEBULA analysis:

apeglm.input.counts <- neb.as.data$count[, order(colnames(neb.as.data$count))]
apeglm.input.counts <- apeglm.input.counts[(row.names(apeglm.input.counts) %in% row.names(apeglm.input.est)), ]

# the counts have to be integers
apeglm.input.counts <- apply(apeglm.input.counts, c(1, 2), function(x) {
  (as.integer(round(x)))
})

# generate disp estimates for apeglm using edgeR's functions
# first, generate dgelist object
cts.dgelist.as <- DGEList(cts.as)
cts.dgelist.as <- scaleOffset(cts.dgelist.as, normMat.as)
rna.design.as.init <- read.csv("AS_DE_Data/artificial_selection_design.csv", header = T, stringsAsFactors = T)
rna.design.as <- rna.design.as.init[-1]
row.names(rna.design.as) <- rna.design.as.init$Sample
samples.as <- merge(cts.dgelist.as$samples, rna.design.as, by = "row.names", all = T)
cts.dgelist.as$samples <- merge(cts.dgelist.as$samples, rna.design.as, by = "row.names", all = T)
cts.dgelist.as$samples <- samples.as[-1]
row.names(cts.dgelist.as$samples) <- samples.as$Row.names

cts.dgelist.as$counts <- cts.dgelist.as$counts[, order(colnames(cts.dgelist.as$counts))]
cts.dgelist.as$offset <- cts.dgelist.as$offset[, order(colnames(cts.dgelist.as$offset))]

cts.dgelist.as.full <- cts.dgelist.as

neb.as.data.for.apeglm <- neb.as.data$count[(row.names(neb.as.data$count) %in% row.names(apeglm.input.est)), ]
cts.dgelist.as.full$counts <- cts.dgelist.as.full$counts[(row.names(cts.dgelist.as.full$counts) %in% row.names(apeglm.input.est)), ]
cts.dgelist.as.full$offset <- cts.dgelist.as.full$offset[(row.names(cts.dgelist.as.full$offset) %in% row.names(apeglm.input.est)), ]

# generate dispersion estimates using edgeR
disp.est.for.apeglm <- estimateDisp(cts.dgelist.as.full, 
                                    neb.mat.alpha)

# Takes ~5 mins to run. MAP estimates = maximum a posteriori estimates
apeglm.AS.fit <- apeglm(Y = apeglm.input.counts, 
                        x = neb.mat.alpha, 
                        log.lik = logLikNB,
                        param = disp.est.for.apeglm$tagwise.dispersion,
                        coef = 3, # corresponds to trt.nebU coeff
                        threshold = log(2)*1,
                        mle = log(2)*apeglm.input.est,
                        offset = disp.est.for.apeglm$offset)

# compare the shrunken estimates to the ones from NEBULA

plot(apeglm.input.est$logFC_trt.nebU, (log2(exp(1))*apeglm.AS.fit$map[,3]))
plot(apeglm.input.est$logFC_trt.nebU, (log2(exp(1))*apeglm.AS.fit$map[,3]), xlim = c(-10, 10), ylim = c(-10, 10))
plot(apeglm.input.est$logFC_trt.nebU, (log2(exp(1))*apeglm.AS.fit$map[,3]), xlim = c(-5, 5), ylim = c(-5, 5))

apeglm.AS.output <- data.frame(row.names = row.names(apeglm.AS.fit$map),
                               trt.UvD.logFC = (log2(exp(1))*apeglm.AS.fit$map[,3]),
                               trt.UvD.sval = apeglm.AS.fit$svalue[,1],
                               trt.UvD.FSR = apeglm.AS.fit$fsr[,1],
                               trt.UvD.2.5 = apeglm.AS.fit$interval[,1],
                               trt.UvD.97.5 = apeglm.AS.fit$interval[,2])

apeglm.AS.output <- apeglm.AS.output[complete.cases(apeglm.AS.output), ]

plot(apeglm.AS.output$trt.UvD.logFC ~ apeglm.AS.output$trt.UvD.sval)
plot(apeglm.AS.output$trt.UvD.logFC ~ apeglm.AS.output$trt.UvD.FSR)

# cor with p vals generated from nebula (pre shrinkage) with s values from apeglm

apeglm.AS.output <- merge(apeglm.AS.output, neb.mod.1$summary[, c(3, 15)], by = "row.names", all.x = F)
row.names(apeglm.AS.output) <- apeglm.AS.output$Row.names
apeglm.AS.output <- apeglm.AS.output[, -1]

plot(apeglm.AS.output$trt.UvD.sval ~ apeglm.AS.output$p_trt.nebU, xlim=c(0,.5), ylim=c(0,.2))
plot(-log10(apeglm.AS.output$trt.UvD.FSR), -log10(apeglm.AS.output$p_trt.nebU))

# determine "hits" from apeglm output

plot(apeglm.AS.output$trt.UvD.sval, apeglm.AS.output$trt.UvD.FSR)

# 0.01 s value cutoff

apeglm.AS.hits <- apeglm.AS.output[which(apeglm.AS.output$trt.UvD.sval < 0.01), ]
apeglm.AS.nonhits <- apeglm.AS.output[which(apeglm.AS.output$trt.UvD.sval >= 0.01), ]

plot(apeglm.AS.hits$trt.UvD.sval, apeglm.AS.hits$trt.UvD.logFC)


########### Shrink estimates for Ten presence effect

neb.mat.alpha <- model.matrix(~ trt.neb + teneral.neb + trt.neb:teneral.neb,
                              data = as.data.frame(neb.as.data$pred)[order(neb.as.data$pred$samples.neb), ])

# estimates we will be shrinking:

apeglm.input.est.ten <- neb.mod.1$summary[, c("logFC_teneral.nebT", "se_teneral.nebT")]
apeglm.input.est.ten <- apeglm.input.est.ten[order(row.names(apeglm.input.est.ten)), ]
apeglm.input.est.ten <- apeglm.input.est.ten[complete.cases(apeglm.input.est.ten), ]

# only looking at those genes that were not cut in the above NEBULA analysis:

apeglm.input.counts <- neb.as.data$count[, order(colnames(neb.as.data$count))]
apeglm.input.counts <- apeglm.input.counts[(row.names(apeglm.input.counts) %in% row.names(apeglm.input.est.ten)), ]

# counts have to be integers
apeglm.input.counts <- apply(apeglm.input.counts, c(1, 2), function(x) {
  (as.integer(round(x)))
})

neb.as.data.for.apeglm.ten <- neb.as.data$count[(row.names(neb.as.data$count) %in% row.names(apeglm.input.est.ten)), ]
cts.dgelist.as.full$counts <- cts.dgelist.as.full$counts[(row.names(cts.dgelist.as.full$counts) %in% row.names(apeglm.input.est.ten)), ]
cts.dgelist.as.full$offset <- cts.dgelist.as.full$offset[(row.names(cts.dgelist.as.full$offset) %in% row.names(apeglm.input.est.ten)), ]

disp.est.for.apeglm <- estimateDisp(cts.dgelist.as.full, 
                                    neb.mat.alpha)


# Takes ~5 mins to run. MAP estimates = maximum a posteriori estimates
apeglm.AS.fit.ten <- apeglm(Y = apeglm.input.counts, 
                            x = neb.mat.alpha, 
                            log.lik = logLikNB,
                            param = disp.est.for.apeglm$tagwise.dispersion,
                            coef = 4, # corresponds to teneral.nebT coeff
                            threshold = log(2)*1,
                            mle = log(2)*apeglm.input.est.ten,
                            offset = disp.est.for.apeglm$offset)

apeglm.AS.output.ten <- data.frame(row.names = row.names(apeglm.AS.fit.ten$map),
                                   trt.TvC.logFC = (log2(exp(1))*apeglm.AS.fit.ten$map[,3]),
                                   trt.TvC.sval = apeglm.AS.fit.ten$svalue[,1],
                                   trt.TvC.FSR = apeglm.AS.fit.ten$fsr[,1],
                                   trt.TvC.2.5 = apeglm.AS.fit.ten$interval[,1],
                                   trt.TvC.97.5 = apeglm.AS.fit.ten$interval[,2])

apeglm.AS.output.ten <- apeglm.AS.output.ten[complete.cases(apeglm.AS.output.ten), ]

# determine "hits" from apeglm output

plot(apeglm.AS.output.ten$trt.TvC.sval, apeglm.AS.output.ten$trt.TvC.FSR)

# 0.01 s value cutoff

apeglm.AS.hits.ten <- apeglm.AS.output.ten[which(apeglm.AS.output.ten$trt.TvC.sval < 0.01), ]
apeglm.AS.nonhits.ten <- apeglm.AS.output.ten[which(apeglm.AS.output.ten$trt.TvC.sval >= 0.01), ]

dim(apeglm.AS.hits.ten) # 12 hits


#####
# Back to analysis of nebula output
# Sorting hits by logFC brings up genes with most samples ~0, and a few high

# Use TPM from abundance:
tpm.as <- txi.eR.as$abundance

# Generate lots of countplots to inspect genes with largest effects
arrangePlots(multiCountPlot(neb.as.hits$gene[1:25], tpm.as, group.as, Trt.as), LeftText = "TPM", TopText = "1:25", ncols = 5)
arrangePlots(multiCountPlot(neb.as.hits$gene[26:50], tpm.as, group.as, Trt.as), LeftText = "TPM", TopText = "26:50", ncols = 5)
arrangePlots(multiCountPlot(neb.as.hits$gene[51:75], tpm.as, group.as, Trt.as), LeftText = "TPM", TopText = "51:75", ncols = 5)
arrangePlots(multiCountPlot(neb.as.hits$gene[76:100], tpm.as, group.as, Trt.as), LeftText = "TPM", TopText = "76:100", ncols = 5)
arrangePlots(multiCountPlot(neb.as.hits$gene[101:125], tpm.as, group.as, Trt.as), LeftText = "TPM", TopText = "101:125", ncols = 5)
arrangePlots(multiCountPlot(neb.as.hits$gene[126:150], tpm.as, group.as, Trt.as), LeftText = "TPM", TopText = "126:150", ncols = 5)
arrangePlots(multiCountPlot(neb.as.hits$gene[151:175], tpm.as, group.as, Trt.as), LeftText = "TPM", TopText = "151:175", ncols = 5)
arrangePlots(multiCountPlot(neb.as.hits$gene[176:200], tpm.as, group.as, Trt.as), LeftText = "TPM", TopText = "176:200", ncols = 5)

# Assemble list of hits
neb.as.hits.narrowed <- head(neb.as.hits, n = 200)
neb.as.rm <- c("FBgn0004431",	"FBgn0266442",	"FBgn0051668",	"FBgn0030775",	"FBgn0030334",	"FBgn0051792",
               "FBgn0053700",	"FBgn0262008",	"FBgn0250824",	"FBgn0034321",	"FBgn0038888",	"FBgn0067311",
               "FBgn0262894",	"FBgn0033774",	"FBgn0031412",	"FBgn0036951",	"FBgn0038715",	"FBgn0036563",
               "FBgn0284444",	"FBgn0262531",	"FBgn0052985",	"FBgn0040759")

neb.as.hits.narrowed <- neb.as.hits.narrowed[!row.names(neb.as.hits.narrowed) %in% neb.as.rm, ]


######
###### glmmTMB tests to verify results from NEBULA using a model that includes vial replicate and day

tpm.as <- txi.eR.as$abundance

trt.neb <- relevel(trt.neb, ref = "D")

# Use normalization factors from voom as offset
# generate limma-voom object for artificial selection
raw.counts.as <- txi.eR.as$counts
d0.as <- DGEList(raw.counts.as)
d0.as <- calcNormFactors(d0.as)
cutoff <- 5
drop.as <- which(apply(cpm(d0.as), 1, max) < cutoff)
d.as <- d0.as[-drop.as,] 
samples.as <- colnames(raw.counts.as)
Trt.as <- as.factor(substr(samples.as, 1, 1))
Trt.as <- relevel(Trt.as, ref = "D")
Teneral_presence.as <- as.factor(substr(samples.as, 3, 3))

mm.as.2a <- model.matrix(~ Trt.as + Teneral_presence.as + Trt.as:Teneral_presence.as) 

y.as.neb <- voom(d0.as, mm.as.2a, plot = T) #

tmb.offset <- y.as.neb$targets$norm.factors


# Takes about 5 mins to run. First runs the full model; if it gives NAs in the summary, runs a reduced model.
tmb.neb.genes.check <- list()

for (g in row.names(neb.as.hits.narrowed)) {
  tpm <- as.vector(tpm.as[g, ])
  tpm[tpm == 0] <- 1e-06 # change zeroes to 1 transcript per million? (minimum nonzero in tpm.as is 3e-4)
  tmb.data <- data.frame(tpm, samples.as, trt.neb, teneral.neb, day.neb, lineage.neb, vial_rep.neb)
  tmb.mod <- glmmTMB(log2(tpm) ~ trt.neb + 
                       teneral.neb + 
                       trt.neb:teneral.neb + 
                       (1 | day.neb) + 
                       (1 + teneral.neb | lineage.neb/vial_rep.neb), 
                     offset = tmb.offset,
                     data = tmb.data)
  if (anyNA(as.vector(summary(tmb.mod)$coef$cond[,2:4])) == T) {
    tmb.mod.reduced <- glmmTMB(log2(tpm) ~ trt.neb + 
                                 teneral.neb + 
                                 trt.neb:teneral.neb + 
                                 (1 | lineage.neb), 
                               offset = tmb.offset,
                               data = tmb.data)
    tmb.neb.genes.check[[g]] <- list(summary(tmb.mod.reduced)$coef$cond, "Used Reduced Model")
  } else {
    tmb.neb.genes.check[[g]] <- list(summary(tmb.mod)$coef$cond, "Used Full Model")
  }
}


# Extract estimates for inspection

tmb.UvD.estimates <- data.frame(FBgnID = character(), Estimate = double(), StdError = double(), P = double(), Model = character())

for (g in names(tmb.neb.genes.check)) {
  df <- data.frame(FBgnID = g, 
                   Estimate = tmb.neb.genes.check[g][[1]][[1]][3,1], 
                   StdError = tmb.neb.genes.check[g][[1]][[1]][3,2], 
                   P = tmb.neb.genes.check[g][[1]][[1]][3,4], 
                   Model = tmb.neb.genes.check[g][[1]][[2]])
  tmb.UvD.estimates <- rbind(tmb.UvD.estimates, df)
}
tmb.UvD.estimates

#write.csv(tmb.UvD.estimates, "tmb_UvD_estimates.csv")

# Compare estimates to logFC from Nebula

row.names(tmb.UvD.estimates) <- tmb.UvD.estimates$FBgnID
tmb.nebula.estimates <- merge(tmb.UvD.estimates, neb.as.hits.narrowed, by = "row.names")

plot(tmb.nebula.estimates$Estimate ~ tmb.nebula.estimates$logFC_trt.nebU)

cor.test(tmb.nebula.estimates$Estimate, tmb.nebula.estimates$logFC_trt.nebU) # r = 0.8


#### Top genes based on above, TMB, and GO analysis. Using voom logCPM values for the figures at this point!

neb.top.genes <- c("FBgn0038945", "FBgn0037562", "FBgn0039613", "FBgn0032416",	"FBgn0051148",	"FBgn0034866",	"FBgn0027106",	"FBgn0026399",	
                   "FBgn0031111",	"FBgn0033885",	"FBgn0036988",	"FBgn0035770",	"FBgn0034786",	
                   "FBgn0004620",	"FBgn0039747",	"FBgn0037976")

###### Prior EdgeR analysis for reference

##### Make DGEList object (taken from tximport vignette)

cts.as <- txi.eR.as$counts
normMat.as <- txi.eR.as$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat.as <- normMat.as/exp(rowMeans(log(normMat.as)))
normCts.as <- cts.as/normMat.as

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib.as <- calcNormFactors(normCts.as) * colSums(normCts.as)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat.as <- sweep(normMat.as, 2, eff.lib.as, "*")
normMat.as <- log(normMat.as)

# Creating a DGEList object for use in edgeR.
cts.dgelist.as <- DGEList(cts.as)
cts.dgelist.as <- scaleOffset(cts.dgelist.as, normMat.as)

# import design matrix
rna.design.as.init <- read.csv("AS_DE_Data/artificial_selection_design.csv", header = T, stringsAsFactors = T)

# make sample the rownames
rna.design.as <- rna.design.as.init[-1]
row.names(rna.design.as) <- rna.design.as.init$Sample

# add design to DGElist
samples.as <- merge(cts.dgelist.as$samples, rna.design.as, by = "row.names", all = T)
cts.dgelist.as$samples <- merge(cts.dgelist.as$samples, rna.design.as, by = "row.names", all = T)
cts.dgelist.as$samples <- samples.as[-1]
row.names(cts.dgelist.as$samples) <- samples.as$Row.names

# reorder column names in count and offest matrices to match sample df

cts.dgelist.as$counts <- cts.dgelist.as$counts[, order(colnames(cts.dgelist.as$counts))]
cts.dgelist.as$offset <- cts.dgelist.as$offset[, order(colnames(cts.dgelist.as$offset))]

cts.dgelist.as.full <- cts.dgelist.as
# Keep a copy of the whole list for use with APEGLM (see NEBULA script)

# filtering out lowly expressed genes. Vignette says to do this by group you plan to look at DE of
keep.as <- filterByExpr(cts.dgelist.as, 
                        group = cts.dgelist.as$samples$Trt,
                        min.count = 5)

table(keep.as)

cts.dgelist.as <- cts.dgelist.as[keep.as, , keep.lib.sizes=FALSE]

##### Plot mds

plotMDS(cts.dgelist.as, 
        labels = cts.dgelist.as$samples$Trt)

plotMDS(cts.dgelist.as, 
        labels = cts.dgelist.as$samples$Day)

plotMDS(cts.dgelist.as, 
        labels = cts.dgelist.as$samples$Teneral_presence)


##### More complex model (trt main effect)

cts.dgelist.as$samples$Trt <- relevel(cts.dgelist.as$samples$Trt, "D")

as.design.2 <- model.matrix(~ Day + Teneral_presence + Trt*Teneral_presence, 
                            data = cts.dgelist.as$samples)

cts.dgelist.trt.as <- estimateDisp(cts.dgelist.as, 
                                   as.design.2)
plotBCV(cts.dgelist.trt.as)

# GLM

as.fit.trt <- glmQLFit(cts.dgelist.trt.as, 
                       design = as.design.2)
plotQLDisp(as.fit.trt)

as.trt.test <- glmQLFTest(as.fit.trt, 
                          coef = 4:5)

summary(decideTests(as.trt.test)) 

edgeR.AS.trt <- topTags(as.trt.test, n = Inf)

##### Interaction

as.int.test <- glmQLFTest(as.fit.trt,
                          coef = 7)

summary(decideTests(as.int.test))  # 0

##### Treatment

as.trt.test.main.eff <- glmQLFTest(as.fit.trt,
                                   coef = 5)

summary(decideTests(as.trt.test.main.eff))

edgeR.AS.trt.fullList <- topTags(as.trt.test.main.eff, n = Inf)
edgeR.AS.trt.p.0.25.cutoff <- topTags(as.trt.test.main.eff, p.value = 0.25, n = Inf)

##### Make the "group" variable like we did with DESeq

cts.dgelist.as$samples$group <- factor(paste0(cts.dgelist.as$samples$Teneral_presence,
                                              cts.dgelist.as$samples$Trt))

as.design.3 <- model.matrix(~ 0 + Day + group,  
                            data = cts.dgelist.as$samples)

cts.dgelist.group.as <- estimateDisp(cts.dgelist.as, 
                                     as.design.3)
plotBCV(cts.dgelist.group.as)

# GLM

as.fit.group <- glmQLFit(cts.dgelist.group.as, 
                         design = as.design.3)
plotQLDisp(as.fit.group)

as.group.test <- glmQLFTest(as.fit.group)

summary(decideTests(as.group.test))

# Group contrasts

TUvsTD <- makeContrasts(groupTD-groupTU, levels = as.design.3)

TUvsTD.test <- glmQLFTest(as.fit.group, 
                          contrast = TUvsTD)

summary(decideTests(TUvsTD.test)) 

topTags(TUvsTD.test, n = 20)

