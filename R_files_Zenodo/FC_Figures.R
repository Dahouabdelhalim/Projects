# The genetic basis of variation in sexual aggression: evolution versus social plasticity

# Figures and Summary
# Andrew M. Scott

#### NOTE before running:
#### Script depends on objects generated from running FC_Vector_Corr, FC_Plasticity_DE and FC_Artificial_Selection_DE Scripts!

library(eulerr)
library(grid)
library(gridExtra)

#############################################################
#### MD plots (plasticity and artificial selection)

#####
# plasticity edgeR treatment (I-E) contrast MD
plast.edger.TRT.MD <- plotMD(trt.test) # Upregulated = higher expression in the Isolated (higher FC rate)

# Make plot by hand

sig.plast.trt.MD <- decideTests(trt.test)
plast.trt.MD <- edger.trt
plast.trt.MD <- merge(plast.trt.MD$table, sig.plast.trt.MD, by = "row.names")
rownames(plast.trt.MD) <- plast.trt.MD$Row.names
plast.trt.MD <- plast.trt.MD[, -1]
plast.trt.MD$Sig <- cut(plast.trt.MD$TrtIso, c(-Inf, -0.1, 0.1, Inf), c("Down", "NotSig", "Up"))
plast.trt.MD.hits <- plast.trt.MD[which(plast.trt.MD$Sig != "NotSig"), ]

# MD plot of edgeR logFC and significant hits. Note that edgeR does do a shrinkage step (see plotQLDisp - shows before and after shrinkage)

Plast.md.plot <- ggplot(data = plast.trt.MD %>% arrange(match(Sig, c("NotSig", "Down", "Up"))), 
                        aes(x = logCPM, 
                            y = logFC, 
                            color = Sig,
                            size = Sig)) +
  geom_point(shape = 16) +
  theme_classic() +
  labs(x = "logCPM", 
       y = "log Fold Change", 
       title = "Plasticity - edgeR") +
  scale_color_manual(limits = c("Up", "NotSig", "Down"), values = c("red", "black", "blue")) +
  scale_size_manual(limits = c("Up", "NotSig", "Down"), values = c(1.5, 0.5, 1.5)) +
  scale_x_continuous(limits = c(-3,15.5), breaks = c(0,5,10,15), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3.5,3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3), expand = c(0, 0)) +
  theme(plot.title = element_text(size = 12, color = "black", vjust = -5, hjust = 0.5),
        legend.position = "none", 
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 12, color = "black"))

# MD plot of apeglm shrunken logFC and sig hits

######
# AS Nebula treatment (U-D) contrast MD - have to make it from the data itself (nebula object not usable by plotMD())
# Again, upregulated (positive logFC) = higher expression in the Up (high FC) lineages
AS.trt.MD <- neb.mod.1$summary[, c(3, 15, 20)]

#obtain AS logCPM from edgeR, combine with estimates from NEBULA
AS.trt.MD <- merge(AS.trt.MD, edgeR.AS.trt$table[, -c(1,2,4:5)], by = "row.names", all.x = F)
row.names(AS.trt.MD) <- AS.trt.MD$Row.names
AS.trt.MD <- AS.trt.MD[, -1]
AS.trt.MD <- AS.trt.MD[order(AS.trt.MD$p_trt.nebU, decreasing = F), ]

# assign significant or non-significant variable indicator (1 = sig up, -1 = sig down, 0 = non sig)
AS.trt.MD.hits <- head(AS.trt.MD, n = length(which(AS.trt.MD$p_trt.nebU < 0.05)))
AS.trt.MD.hits <- AS.trt.MD.hits[order(AS.trt.MD.hits$logFC_trt.nebU, decreasing = F), ]
AS.trt.MD.down.hits <- head(AS.trt.MD.hits, n = length(which(AS.trt.MD.hits$logFC_trt.nebU < 0)))
AS.trt.MD.up.hits <- tail(AS.trt.MD.hits, n = length(which(AS.trt.MD.hits$logFC_trt.nebU > 0)))
AS.trt.MD.nonhits <- tail(AS.trt.MD, n = length(which(AS.trt.MD$p_trt.nebU >= 0.05)))
AS.trt.MD.down.hits$Sig <- "Down"
AS.trt.MD.up.hits$Sig <- "Up"
AS.trt.MD.nonhits$Sig <- "NotSig"
AS.trt.MD.hits <- rbind(AS.trt.MD.down.hits, AS.trt.MD.up.hits)
AS.trt.MD <- rbind(AS.trt.MD.hits, AS.trt.MD.nonhits)

# make df with shrunken estimates from apeglm
# also add in "sig" or "not sig" based on s-value from apeglm analysis

apeglm.AS.hits <- apeglm.AS.hits[order(apeglm.AS.hits$trt.UvD.logFC, decreasing = F), ]

Sig_shrunk.Up <- apeglm.AS.hits[which(apeglm.AS.hits$trt.UvD.logFC > 0), ]   
Sig_shrunk.Down <- apeglm.AS.hits[which(apeglm.AS.hits$trt.UvD.logFC < 0), ]
Sig_shrunk.Up$Sig <- "Up"
Sig_shrunk.Down$Sig <- "Down"
apeglm.AS.nonhits$Sig <- "NotSig"
apeglm.AS.hits <- rbind(Sig_shrunk.Up, Sig_shrunk.Down)

apeglm.AS.nonhits$gene <- rownames(apeglm.AS.nonhits)
apeglm.AS.hits$gene <- rownames(apeglm.AS.hits)

apeglm.AS.MD <- rbind(apeglm.AS.hits, apeglm.AS.nonhits)

# use edgeR to obtain logCPM values for AS

# filtering out lowly expressed genes
keep.as <- filterByExpr(cts.dgelist.as, 
                        group = cts.dgelist.as$samples$Trt,
                        min.count = 5)
cts.dgelist.as <- cts.dgelist.as[keep.as, , keep.lib.sizes=FALSE]
#trt main effect model
cts.dgelist.as$samples$Trt <- relevel(cts.dgelist.as$samples$Trt, "D")
as.design.2 <- model.matrix(~ Day + Teneral_presence + Trt*Teneral_presence, 
                            data = cts.dgelist.as$samples)
cts.dgelist.trt.as <- estimateDisp(cts.dgelist.as, 
                                   as.design.2)
# GLM
as.fit.trt <- glmQLFit(cts.dgelist.trt.as, 
                       design = as.design.2)
as.trt.test <- glmQLFTest(as.fit.trt, 
                          coef = 4:5)
edgeR.AS.trt <- topTags(as.trt.test, n = Inf)

# add logCPM from edgeR

apeglm.AS.MD <- merge(apeglm.AS.MD, edgeR.AS.trt$table[, -c(1,2,4:5)], by = "row.names", all.x = F)
row.names(apeglm.AS.MD) <- apeglm.AS.MD$Row.names
apeglm.AS.MD <- apeglm.AS.MD[, -1]


# create MD plots

AS.md.plot <- ggplot(data = AS.trt.MD %>% arrange(match(Sig, c("NotSig", "Down", "Up"))), 
                     aes(x = logCPM, 
                         y = logFC_trt.nebU, 
                         color = Sig,
                         size = Sig)) +
  geom_point(shape = 16) +
  theme_classic() +
  labs(x = "logCPM", 
       y = "log Fold Change", 
       title = "Artificial Selection - NEBULA") +
  scale_color_manual(limits = c("Up", "NotSig", "Down"), values = c("red", "black", "blue")) +
  scale_size_manual(limits = c("Up", "NotSig", "Down"), values = c(1.5, 0.5, 1.5)) +
  scale_x_continuous(limits = c(-3,15.5), breaks = c(0,5,10,15), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3.5,3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3), expand = c(0, 0)) +
  theme(plot.title = element_text(size = 12, color = "black", vjust = -5, hjust = 0.5),
        legend.position = "none", 
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 12, color = "black"))

AS.md.plot.shrunk <- ggplot(data = apeglm.AS.MD %>% arrange(match(Sig, c("NotSig", "Down", "Up"))), 
                            aes(x = logCPM, 
                                y = trt.UvD.logFC, 
                                color = Sig,
                                size = Sig)) +
  geom_point(shape = 16) +
  theme_classic() +
  labs(x = "logCPM", 
       y = "log Fold Change", 
       title = "Artificial Selection - APEGLM shrunken") +
  scale_color_manual(limits = c("Up", "NotSig", "Down"), values = c("red", "black", "blue")) +
  scale_size_manual(limits = c("Up", "NotSig", "Down"), values = c(1.5, 0.5, 1.5)) +
  scale_x_continuous(limits = c(-3,15.5), breaks = c(0,5,10,15), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-3.5,3.5), breaks = c(-3, -2, -1, 0, 1, 2, 3), expand = c(0, 0)) +
  theme(plot.title = element_text(size = 12, color = "black", vjust = -5, hjust = 0.5),
        legend.position = "none", 
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 12, color = "black"))

# make plot for both shrunken and non-shrunken AS

grid.arrange(arrangeGrob(AS.md.plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank()), 
                         AS.md.plot.shrunk + theme(axis.title.y = element_blank(), axis.title.x = element_blank()), 
                         ncol = 2, 
                         left = textGrob("log Fold Change", rot = 90, vjust = 0.5, gp = gpar(fontsize = 16)), 
                         bottom = textGrob("log Counts per Million", gp = gpar(fontsize = 16))))

# Final main MS plot

grid.arrange(arrangeGrob(AS.md.plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank()), 
                         AS.md.plot.shrunk + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                         Plast.md.plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                         ncol = 3, 
                         left = textGrob("log2 fold change", rot = 90, vjust = 0.5, gp = gpar(fontsize = 16)), 
                         bottom = textGrob("Average log2 counts per million", gp = gpar(fontsize = 16))))

###### Volcano plots
#### Plasticity

plast.volc.plot <- ggplot(data = plast.trt.MD %>% arrange(match(Sig, c("NotSig", "Down", "Up"))),
                          aes(x = logFC, 
                              y = -log10(PValue),
                              color = Sig,
                              size = Sig)) +
  geom_point(shape = 16) +
  theme_classic() +
  labs(x = "log Fold Change", 
       y = "-log10(P value)", 
       title = "Plasticity - edgeR") +
  scale_color_manual(limits = c("Up", "NotSig", "Down"), values = c("red", "black", "blue")) +
  scale_size_manual(limits = c("Up", "NotSig", "Down"), values = c(1.5, 0.5, 1.5)) +
  scale_x_continuous(limits = c(-3.1,3.1), breaks = c(-3,-2,-1,0,1,2,3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8), expand = c(0, 0)) +
  theme(plot.title = element_text(size = 12, color = "black", vjust = -5, hjust = 0.5),
        legend.position = "none", 
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 12, color = "black"))

#### Artificial selection

AS.volc.plot <- ggplot(data = AS.trt.MD %>% arrange(match(Sig, c("NotSig", "Down", "Up"))),
                       aes(x = logFC_trt.nebU, 
                           y = -log10(p_trt.nebU),
                           color = Sig,
                           size = Sig)) +
  geom_point(shape = 16) +
  theme_classic() +
  labs(x = "log Fold Change", 
       y = "-log10(P value)", 
       title = "Artificial Selection - NEBULA") +
  scale_color_manual(limits = c("Up", "NotSig", "Down"), values = c("red", "black", "blue")) +
  scale_size_manual(limits = c("Up", "NotSig", "Down"), values = c(1.5, 0.5, 1.5)) +
  scale_x_continuous(limits = c(-3.5,3.5), breaks = c(-3,-2,-1,0,1,2,3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0, 2, 4, 6, 8), expand = c(0, 0)) +
  theme(plot.title = element_text(size = 12, color = "black", vjust = -5, hjust = 0.5),
        legend.position = "none", 
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 12, color = "black"))

AS.volc.plot.shrunk <- ggplot(data = apeglm.AS.MD %>% arrange(match(Sig, c("NotSig", "Down", "Up"))),
                              aes(x = trt.UvD.logFC, 
                                  y = -log10(trt.UvD.sval),
                                  color = Sig,
                                  size = Sig)) +
  geom_point(shape = 16) +
  theme_classic() +
  labs(x = "log Fold Change", 
       y = "-log10(S value)", 
       title = "Artificial Selection - APEGLM shrunken") +
  scale_color_manual(limits = c("Up", "NotSig", "Down"), values = c("red", "black", "blue")) +
  scale_size_manual(limits = c("Up", "NotSig", "Down"), values = c(1.5, 0.5, 1.5)) +
  scale_x_continuous(limits = c(-3.5,3.5), breaks = c(-3,-2,-1,0,1,2,3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 10.5), breaks = c(0, 2, 4, 6, 8, 10), expand = c(0, 0)) +
  theme(plot.title = element_text(size = 12, color = "black", vjust = -5, hjust = 0.5),
        legend.position = "none", 
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 12, color = "black"))

grid.arrange(arrangeGrob(AS.volc.plot + theme(axis.title.x = element_blank()), 
                         AS.volc.plot.shrunk + theme(axis.title.x = element_blank()),
                         ncol = 2, 
                         bottom = textGrob("log Fold Change", gp = gpar(fontsize = 16))))


### Final MA + Volcano plots

grid.arrange(arrangeGrob(AS.md.plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank()), 
                         AS.md.plot.shrunk + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                         Plast.md.plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                         ncol = 3, 
                         left = textGrob("log2 fold change", rot = 90, vjust = 0.5, gp = gpar(fontsize = 16)), 
                         bottom = textGrob("Average log2 counts per million", gp = gpar(fontsize = 16))),
             arrangeGrob(AS.volc.plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank()), 
                         AS.volc.plot.shrunk + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                         plast.volc.plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                         ncol = 3,
                         left = textGrob("-log10(P value)", rot = 90, vjust = 0.5, gp = gpar(fontsize = 16)),
                         bottom = textGrob("log2 Fold Change", gp = gpar(fontsize = 16))))


#############################################################
#### Venn Diagram
#### Overlapping hits...

# Lists of hits: Based on the edgeR hits, and the either the shrunken or unshrunk AS estimates
overlap.AS.plast.figure <- merge(AS.trt.MD.hits, plast.trt.MD.hits, by = "row.names", all.x = F)
overlap.ASshrunk.plast.figure <- merge(apeglm.AS.hits, plast.trt.MD.hits, by = "row.names", all.x = F)

# Venn for each of the above sets
AS.trt.MD.hits$gene <- row.names(AS.trt.MD.hits)
apeglm.AS.hits$gene <- row.names(apeglm.AS.hits)
plast.trt.MD.hits$gene <- row.names(plast.trt.MD.hits)

venn.as.plast.data <- list(set1 = AS.trt.MD.hits$gene,
                           set2 = plast.trt.MD.hits$gene)

venn.as.plast <- plot(euler(venn.as.plast.data), quantities = T, 
                      labels = c("Artificial selection\\n(NEBULA)", "Plasticity (EdgeR -\\nIvE)"))

venn.as_shrunk.plast.data <- list(set1 = apeglm.AS.hits$gene,
                                  set2 = plast.trt.MD.hits$gene)

venn.as_shrunk.plast <- plot(euler(venn.as_shrunk.plast.data), quantities = T, 
                             labels = c("Artificial selection\\n(APEGLM shrunken)", "Plasticity (EdgeR -\\nIvE)"))

grid.arrange(arrangeGrob(venn.as.plast, venn.as_shrunk.plast,
                         ncol = 1))

#### second set of Venns that includes plasticity genes significant in EITHER IvE or ITvET contrasts
# (This mirrors our selection process for candidate genes)

#### For plasticity - combine hits from IvE main effect, and TIvTE contrast
# THIS IS FROM THE VECTOR CORRELATION SCRIPT

plast.hits.both.contrasts <- plast.hits.vc

# Overlapped hits

overlap.AS.plast_both.figure <- merge(AS.trt.MD.hits, plast.hits.vc, by = "row.names", all.x = F)
overlap.ASshrunk.plast_both.figure <- merge(apeglm.AS.hits, plast.hits.vc, by = "row.names", all.x = F)

# Make venn

plast.hits.both.contrasts$gene <- row.names(plast.hits.both.contrasts)

venn.as.plast.both.data <- list(set1 = AS.trt.MD.hits$gene,
                                set2 = plast.hits.both.contrasts$gene)

venn.as.plast_both <- plot(euler(venn.as.plast.both.data), quantities = T, 
                           labels = c("Artificial selection\\n(NEBULA)", "Plasticity (EdgeR)"))

venn.as_shrunk.plast.both.data <- list(set1 = apeglm.AS.hits$gene,
                                       set2 = plast.hits.both.contrasts$gene)

venn.as_shrunk.plast_both <- plot(euler(venn.as_shrunk.plast.both.data), quantities = T, 
                                  labels = c("Artificial selection\\n(APEGLM shrunken)", "Plasticity (EdgeR)"))

grid.arrange(arrangeGrob(venn.as.plast_both, venn.as_shrunk.plast_both,
                         ncol = 1))


#############################################################
#### Expression plots of overlapping hits. 

# Using the results from limma-voom for log2 Counts per million: y.as.neb$E, y$E

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

row.names(overlap.AS.plast_both.figure) <- overlap.AS.plast_both.figure$Row.names
overlap.AS.plast_both.figure <- overlap.AS.plast_both.figure[, -1]
overlap.AS.plast_both.figure <- overlap.AS.plast_both.figure[, c("gene", "logCPM", "logFC_trt.nebU", "edger_logFC")]
colnames(overlap.AS.plast_both.figure) <- c("gene", "logCPM", "logFC_AS", "logFC_Plast")

row.names(overlap.ASshrunk.plast_both.figure) <- overlap.ASshrunk.plast_both.figure$Row.names
overlap.ASshrunk.plast_both.figure <- overlap.ASshrunk.plast_both.figure[, -1]
overlap.ASshrunk.plast_both.figure <- overlap.ASshrunk.plast_both.figure[, c("gene", "edger_logCPM", "logFC_trt.nebU", "edger_logFC")]
colnames(overlap.ASshrunk.plast_both.figure) <- c("gene", "logCPM", "logFC_AS", "logFC_Plast")

overlapping.hits.expr.pre <- overlap.ASshrunk.plast_both.figure[!(row.names(overlap.ASshrunk.plast_both.figure) %in% row.names(overlap.AS.plast_both.figure)), ]
overlapping.hits.expr <- rbind(overlap.AS.plast_both.figure, overlapping.hits.expr.pre) # 27 total overlapping hits

# order by mean abs. log fold change between AS and plasticity

overlapping.hits.expr$mean_absLFC_across_exp <- rowMeans(abs(overlapping.hits.expr[ ,c("logFC_AS", "logFC_Plast")]))
overlapping.hits.expr <- overlapping.hits.expr[order(overlapping.hits.expr$mean_absLFC_across_exp, decreasing = T), ]

# Using the results from limma-voom for log2 Counts per million: y.as.neb$E, y$E

countPlot_finalFigs <- function(gene, count_file_AS, count_file_Plast) {
  set.seed(1234)
  sorted.count_file_AS <- count_file_AS[, order(colnames(count_file_AS))]
  sorted.count_file_Plast <- count_file_Plast[, order(colnames(count_file_Plast))]
  
  #gene <- deparse(substitute(gene))
  gene.count.as <- sorted.count_file_AS[gene, ]
  gene.count.plast <- sorted.count_file_Plast[gene, ]
  
  Trts.as <- as.factor(substr(names(gene.count.as), 1, 1))
  Trts.plast <- as.factor(substr(names(gene.count.plast), 1, 1))
  
  df.as <- data.frame("logCPM" = gene.count.as, "Trts" = Trts.as)
  df.as <- df.as[!df.as$Trts == "C", ] 
  df.as$Trts <- recode(df.as$Trts, U = "High_FC", D = "Low_FC")
  df.as$Exp <- "AS"
  
  df.plast <- data.frame("logCPM" = gene.count.plast, "Trts" = Trts.plast)
  df.plast$Trts <- recode(df.plast$Trts, I = "High_FC", E = "Low_FC")
  df.plast$Exp <- "Plast"
  
  df.fin <- rbind(df.as, df.plast)
  df.fin.means <- df.fin %>% group_by(Exp, Trts) %>% summarise(logCPM = mean(logCPM), .groups = "drop")
  
  gene.symbol <- gene.symbols[which(grepl(gene, gene.symbols$`primary_FBgn#`)), 1]
  plot <- ggplot(df.fin, aes(y = logCPM, 
                             x = Trts, 
                             color = Exp)) +
    geom_point(pch = 1, cex = 1, 
               position = position_jitterdodge(jitter.width = 0.02, dodge.width = 0.1), 
               stroke = 0.75) +
    geom_line(data = df.fin.means, aes(group = Exp), stat = "identity",
              size = 1) +
    labs(y = "log2 Counts per Million", title = paste(gene.symbol)) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_color_manual(limits = c("AS", "Plast"), values = c("black", "gray")) + 
    scale_x_discrete(limits = c("Low_FC", "High_FC"), labels = c("Low FC", "High FC"))
  plot
}

listofplots <- list()
for (g in overlapping.hits.expr$gene) {
  listofplots[[g]] <- countPlot_finalFigs(gene = g, y.as.neb$E, y$E)
}

arrangePlots(listofplots, ncols = 5, BottomText = "Treatment (Black = AS, Gray = Plasticity)", 
             LeftText = "log Counts per Million", TopText = NULL)  


# Top genes in plast, and AS (not overlapping)

# AS

countPlot_finalFigs_AS <- function(gene, count_file_AS) {
  set.seed(1234)
  sorted.count_file_AS <- count_file_AS[, order(colnames(count_file_AS))]
  
  #gene <- deparse(substitute(gene))
  gene.count.as <- sorted.count_file_AS[gene, ]
  
  Trts.as <- as.factor(substr(names(gene.count.as), 1, 1))
  
  df.as <- data.frame("logCPM" = gene.count.as, "Trts" = Trts.as)
  df.as <- df.as[!df.as$Trts == "C", ] 
  df.as$Trts <- recode(df.as$Trts, U = "High_FC", D = "Low_FC")
  df.as$Exp <- "AS"
  
  df.as.means <- df.as %>% group_by(Exp, Trts) %>% summarise(logCPM = mean(logCPM), .groups = "drop")
  
  gene.symbol <- gene.symbols[which(grepl(gene, gene.symbols$`primary_FBgn#`)), 1]
  plot <- ggplot(df.as, aes(y = logCPM, 
                            x = Trts)) +
    geom_point(pch = 1, cex = 1, 
               position = position_jitter(width = 0.02, height = 0), 
               stroke = 0.75) +
    geom_line(data = df.as.means, aes(group = Exp), stat = "identity",
              size = 1) +
    labs(y = "log2 Counts per Million", title = paste(gene.symbol)) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_discrete(limits = c("Low_FC", "High_FC"), labels = c("Low FC", "High FC"))
  plot
}

neb.top.genes

listofplots.as <- list()
for (g in neb.top.genes) {
  listofplots.as[[g]] <- countPlot_finalFigs_AS(gene = g, y.as.neb$E)
}

arrangePlots(listofplots.as, ncols = 4, BottomText = "Treatment", 
             LeftText = "log Counts per Million", TopText = "Artificial Selection (NEBULA) - Top 20 hits")  

# Plasticity

mm.2 <- model.matrix(~1 + Day + Vial_rep + Trt + Teneral_presence + Trt:Teneral_presence)
y.all.2 <- voom(d0, mm.2, plot = T) # d0 is the list without the cutoff

countPlot_finalFigs_Plast <- function(gene, count_file_Plast) {
  set.seed(1234)
  sorted.count_file_Plast <- count_file_Plast[, order(colnames(count_file_Plast))]
  
  #gene <- deparse(substitute(gene))
  gene.count.plast <- sorted.count_file_Plast[gene, ]
  
  Trts.plast <- as.factor(substr(names(gene.count.plast), 1, 1))
  
  df.plast <- data.frame("logCPM" = gene.count.plast, "Trts" = Trts.plast)
  df.plast$Trts <- recode(df.plast$Trts, I = "High_FC", E = "Low_FC")
  df.plast$Exp <- "Plast"
  
  df.fin <- df.plast
  df.fin.means <- df.fin %>% group_by(Exp, Trts) %>% summarise(logCPM = mean(logCPM), .groups = "drop")
  
  gene.symbol <- gene.symbols[which(grepl(gene, gene.symbols$`primary_FBgn#`)), 1]
  plot <- ggplot(df.fin, aes(y = logCPM, 
                             x = Trts)) +
    geom_point(pch = 1, cex = 1, 
               position = position_jitter(width = 0.02, height = 0), 
               stroke = 0.75) +
    geom_line(data = df.fin.means, aes(group = Exp), stat = "identity",
              size = 1) +
    labs(y = "log2 Counts per Million", title = paste(gene.symbol)) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_discrete(limits = c("Low_FC", "High_FC"), labels = c("Low FC", "High FC"))
  plot
}

plast.hits.both.contrasts.top16 <- head(plast.hits.both.contrasts[order(plast.hits.both.contrasts$edger_logFC, decreasing = T), ], n = 16)

listofplots.plast <- list()
for (g in row.names(plast.hits.both.contrasts.top16)) {
  listofplots.plast[[g]] <- countPlot_finalFigs_Plast(gene = g, y.all.2$E)
}

arrangePlots(listofplots.plast, ncols = 4, BottomText = "Treatment", 
             LeftText = "log Counts per Million", TopText = "Plasticity (EdgeR) - Top 20 hits")  


#########################################################################
######### Vector Correlation and Alpha Bar+line Plots
######### Non-shrunken AS

vec.cor.plot.data <- data.frame(Corr = c("Overlapping genes", "Sig Plast genes", "Sig AS genes"),
                                vector.cor = c(obs.vs.empir.dist.overlap[1,1], Plast.AS.VC.res.vs.empir.dist[1,1], Plast.AS.VC.res.vs.empir.dist.2[1,1]),
                                Lower.2.5 = c(obs.vs.empir.dist.overlap[2,1], Plast.AS.VC.res.vs.empir.dist[2,1], Plast.AS.VC.res.vs.empir.dist.2[2,1]),
                                Upper.97.5 = c(obs.vs.empir.dist.overlap[4,1], Plast.AS.VC.res.vs.empir.dist[4,1], Plast.AS.VC.res.vs.empir.dist.2[4,1]))

random_vectors_overlap.plot <- as.data.frame(random_vectors_overlap[,1]) # overlap list of random vectors
random_vectors_overlap.plot$Corr <- "Overlapping genes"
base::colnames(random_vectors_overlap.plot) <- c("vector.cor", "Corr")

random_vectors_Plast_AS.plot <- as.data.frame(random_vectors_Plast_AS[,1]) # sig plast genes
random_vectors_Plast_AS.plot$Corr <- "Sig Plast genes"
base::colnames(random_vectors_Plast_AS.plot) <- c("vector.cor", "Corr")

random_vectors_AS_hits.plot <- as.data.frame(random_vectors_AS_hits[,1]) # sig as genes
random_vectors_AS_hits.plot$Corr <- "Sig AS genes"
base::colnames(random_vectors_AS_hits.plot) <- c("vector.cor", "Corr")

random_vectors_plot <- rbind(random_vectors_overlap.plot, random_vectors_Plast_AS.plot, random_vectors_AS_hits.plot)

vec.cor.plot <- ggplot(data = vec.cor.plot.data, 
                       aes(x = Corr, 
                           y = vector.cor, 
                           color = Corr)) +
  # geom_point(data = random_vectors_plot, aes(x = Corr, y = vector.cor, color = Corr), alpha = 0.25, 
  #            stroke = 0, size = 0.08, position= position_jitter(width = 0.1)) +
  geom_crossbar(aes(ymin = Lower.2.5, ymax = Upper.97.5), width = 0.5, size = 0.85) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black")) +
  scale_colour_manual(limits = c("Overlapping genes", "Sig AS genes", "Sig Plast genes"), values = c("red", "black", "gray45")) +
  scale_x_discrete(limits = c("Overlapping genes", "Sig AS genes", "Sig Plast genes")) +
  scale_y_continuous(limits = c(0,0.9), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, .6, .7, .8), expand = c(0, 0)) +
  labs(x = NULL,
       y = "Vector Correlation")

alpha.plot.data <- data.frame(Corr = c("Overlapping genes", "Sig Plast genes", "Sig AS genes"),
                              alpha = c(obs.vs.empir.dist.overlap[1,3], Plast.AS.VC.res.vs.empir.dist[1,3], Plast.AS.VC.res.vs.empir.dist.2[1,3]),
                              Lower.2.5 = c(obs.vs.empir.dist.overlap[2,3], Plast.AS.VC.res.vs.empir.dist[2,3], Plast.AS.VC.res.vs.empir.dist.2[2,3]),
                              Upper.97.5 = c(obs.vs.empir.dist.overlap[4,3], Plast.AS.VC.res.vs.empir.dist[4,3], Plast.AS.VC.res.vs.empir.dist.2[4,3]))

random_alphas_overlap.plot <- as.data.frame(random_vectors_overlap[,3]) # overlap list of random vectors
random_alphas_overlap.plot$Corr <- "Overlapping genes"
base::colnames(random_alphas_overlap.plot) <- c("alpha", "Corr")

random_alphas_Plast_AS.plot <- as.data.frame(random_vectors_Plast_AS[,3]) # sig plast genes
random_alphas_Plast_AS.plot$Corr <- "Sig Plast genes"
base::colnames(random_alphas_Plast_AS.plot) <- c("alpha", "Corr")

random_alphas_AS_hits.plot <- as.data.frame(random_vectors_AS_hits[,3]) # sig as genes
random_alphas_AS_hits.plot$Corr <- "Sig AS genes"
base::colnames(random_alphas_AS_hits.plot) <- c("alpha", "Corr")

random_alphas_plot <- rbind(random_alphas_overlap.plot, random_alphas_Plast_AS.plot, random_alphas_AS_hits.plot)

alpha.plot <- ggplot(data = alpha.plot.data, ####### Note for alpha graph - the alpha values are always Plast/AS, regardless of the gene set (so vals >1 mean plast vector has a larger magnitude)
                     aes(x = Corr, 
                         y = alpha, 
                         color = Corr)) +
  #  geom_point(data = random_alphas_plot, aes(x = Corr, y = alpha, color = Corr), alpha = 0.25, 
  #             stroke = 0, size = 0.08, position= position_jitter(width = 0.1)) +
  geom_crossbar(aes(ymin = Lower.2.5, ymax = Upper.97.5), width = 0.5, size = 0.85) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black")) +
  scale_colour_manual(limits = c("Overlapping genes", "Sig AS genes", "Sig Plast genes"), values = c("red", "black", "gray45")) +
  scale_x_discrete(limits = c("Overlapping genes", "Sig AS genes", "Sig Plast genes")) +
  scale_y_continuous(limits = c(0,3.25), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3), expand = c(0, 0)) +
  labs(x = NULL,
       y = "Alpha (Plasticity/AS)")


grid.arrange(arrangeGrob(vec.cor.plot + theme(axis.title.x = element_blank()), 
                         alpha.plot + theme(axis.title.x = element_blank()),
                         top = "Overlapping genes and Plasticity (edgeR) set vs Artificial Selection (NEBULA) set",
                         ncol = 2))


######### Bar+line Plots
######### Shrunken AS

vec.cor.plot.data_AS_shrunk <- data.frame(Corr = c("Sig Plast genes", "Sig AS genes"),
                                          vector.cor = c(Plast.AS_shrunk.VC.res.vs.empir.dist[1,1], AS_shrunk_plast.VC.res.vs.empir.dist[1,1]),
                                          Lower.2.5 = c(Plast.AS_shrunk.VC.res.vs.empir.dist[2,1], AS_shrunk_plast.VC.res.vs.empir.dist[2,1]),
                                          Upper.97.5 = c(Plast.AS_shrunk.VC.res.vs.empir.dist[4,1], AS_shrunk_plast.VC.res.vs.empir.dist[4,1]))

vec.cor.plot_shrunk <- ggplot(data = vec.cor.plot.data_AS_shrunk, 
                              aes(x = Corr, 
                                  y = vector.cor, 
                                  color = Corr)) +
  geom_crossbar(aes(ymin = Lower.2.5, ymax = Upper.97.5), width = 0.5, size = 0.85) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black")) +
  scale_colour_manual(limits = c("Sig AS genes", "Sig Plast genes"), values = c("black", "gray45")) +
  scale_x_discrete(limits = c("Sig AS genes", "Sig Plast genes")) +
  scale_y_continuous(limits = c(0,0.55), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), expand = c(0, 0)) +
  labs(x = NULL,
       y = "Vector Correlation")

alpha.plot.data_shrunk <- data.frame(Corr = c("Sig Plast genes", "Sig AS genes"),
                                     alpha = c(Plast.AS_shrunk.VC.res.vs.empir.dist[1,3], AS_shrunk_plast.VC.res.vs.empir.dist[1,3]),
                                     Lower.2.5 = c(Plast.AS_shrunk.VC.res.vs.empir.dist[2,3], AS_shrunk_plast.VC.res.vs.empir.dist[2,3]),
                                     Upper.97.5 = c(Plast.AS_shrunk.VC.res.vs.empir.dist[4,3], AS_shrunk_plast.VC.res.vs.empir.dist[4,3]))

alpha.plot_shrunk <- ggplot(data = alpha.plot.data_shrunk, ####### Note for alpha graph - the alpha values are always Plast/AS, regardless of the gene set (so vals >1 mean plast vector has a larger magnitude)
                            aes(x = Corr, 
                                y = alpha, 
                                color = Corr)) +
  geom_crossbar(aes(ymin = Lower.2.5, ymax = Upper.97.5), width = 0.5, size = 0.85) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black")) +
  scale_colour_manual(limits = c("Sig AS genes", "Sig Plast genes"), values = c("black", "gray45")) +
  scale_x_discrete(limits = c("Sig AS genes", "Sig Plast genes")) +
  scale_y_continuous(limits = c(0,5.25), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), expand = c(0, 0)) +
  labs(x = NULL,
       y = "Alpha (Plasticity/AS)")


grid.arrange(arrangeGrob(vec.cor.plot_shrunk + theme(axis.title.x = element_blank()), 
                         alpha.plot_shrunk + theme(axis.title.x = element_blank()),
                         top = "Plasticity (edgeR) vs Artificial Selection shrunken (APEGLM)",
                         ncol = 2))

######### Bar+line Plots
######### edgeR AS vs edgeR Plasticity

vec.cor.plot.data_AS_edgeR <- data.frame(Corr = c("Sig Plast genes", "Sig AS genes"),
                                         vector.cor = c(Plast.AS_edgeR.VC.res.vs.empir.dist[1,1], AS_edgeR_plast.VC.res.vs.empir.dist[1,1]),
                                         Lower.2.5 = c(Plast.AS_edgeR.VC.res.vs.empir.dist[2,1], AS_edgeR_plast.VC.res.vs.empir.dist[2,1]),
                                         Upper.97.5 = c(Plast.AS_edgeR.VC.res.vs.empir.dist[4,1], AS_edgeR_plast.VC.res.vs.empir.dist[4,1]))

vec.cor.plot_ASedgeR <- ggplot(data = vec.cor.plot.data_AS_edgeR, 
                               aes(x = Corr, 
                                   y = vector.cor, 
                                   color = Corr)) +
  geom_crossbar(aes(ymin = Lower.2.5, ymax = Upper.97.5), width = 0.5, size = 0.85) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black")) +
  scale_colour_manual(limits = c("Sig AS genes", "Sig Plast genes"), values = c("black", "gray45")) +
  scale_x_discrete(limits = c("Sig AS genes", "Sig Plast genes"), labels = c("AS genes (p < 0.25)", "Sig Plast genes")) +
  scale_y_continuous(limits = c(0,0.55), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), expand = c(0, 0)) +
  labs(x = NULL,
       y = "Vector Correlation")


alpha.plot.data_ASedgeR <- data.frame(Corr = c("Sig Plast genes", "Sig AS genes"),
                                      alpha = c(Plast.AS_edgeR.VC.res.vs.empir.dist[1,3], AS_edgeR_plast.VC.res.vs.empir.dist[1,3]),
                                      Lower.2.5 = c(Plast.AS_edgeR.VC.res.vs.empir.dist[2,3], AS_edgeR_plast.VC.res.vs.empir.dist[2,3]),
                                      Upper.97.5 = c(Plast.AS_edgeR.VC.res.vs.empir.dist[4,3], AS_edgeR_plast.VC.res.vs.empir.dist[4,3]))

alpha.plot_ASedgeR <- ggplot(data = alpha.plot.data_ASedgeR, ####### Note for alpha graph - the alpha values are always Plast/AS, regardless of the gene set (so vals >1 mean plast vector has a larger magnitude)
                             aes(x = Corr, 
                                 y = alpha, 
                                 color = Corr)) +
  geom_crossbar(aes(ymin = Lower.2.5, ymax = Upper.97.5), width = 0.5, size = 0.85) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 16, color = "black")) +
  scale_colour_manual(limits = c("Sig AS genes", "Sig Plast genes"), values = c("black", "gray45")) +
  scale_x_discrete(limits = c("Sig AS genes", "Sig Plast genes"), labels = c("AS genes (p < 0.25)", "Sig Plast genes")) +
  scale_y_continuous(limits = c(0,5.25), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), expand = c(0, 0)) +
  labs(x = NULL,
       y = "Alpha (Plasticity/AS)")


grid.arrange(arrangeGrob(vec.cor.plot_ASedgeR + theme(axis.title.x = element_blank()), 
                         alpha.plot_ASedgeR + theme(axis.title.x = element_blank()),
                         top = "Plasticity (edgeR) vs Artificial Selection (edgeR)",
                         ncol = 2))


################################################################################
#### FOR SUPPLEMENT:
#### LogFC reaction norm plot

#### Significant plasticity genes and corresponding AS genes, by logFC

react.norm.sigplast.data <- vc.sig.plast.genes[, c(1, 8, 25)]
colnames(react.norm.sigplast.data) <- c("Plasticity_logFC", "ArtificialSel_logFC", "Gene")
react.norm.sigplast.data <- gather(react.norm.sigplast.data, Experiment, logFC, Plasticity_logFC:ArtificialSel_logFC)

cor(react.norm.sigplast.data$edger_logFC, react.norm.sigplast.data$logFC_trt.nebU)

react.norm.sigplast <- ggplot(data = react.norm.sigplast.data, 
                              aes(x = Experiment, 
                                  y = logFC)) + 
  geom_point(alpha = 0.3, size = 0.25) + 
  geom_line(aes(group = Gene), size = 0.25, alpha = 0.3) +
  theme_classic() +
  labs(title = "Significant plasticity genes (edgeR) and\\ncorresponding set in AS (NEBULA)") +
  scale_x_discrete(expand = c(0, 0.2), 
                   limits = c("ArtificialSel_logFC", "Plasticity_logFC"),
                   labels = c("Artificial Selection", "Plasticity")) +
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"))

react.norm.sigplast.abs <- ggplot(data = react.norm.sigplast.data, 
                                  aes(x = Experiment, 
                                      y = abs(logFC))) + 
  geom_point(alpha = 0.3, size = 0.25) + 
  geom_line(aes(group = Gene), size = 0.25, alpha = 0.3) +
  theme_classic() +
  labs(title = "Significant plasticity genes (edgeR) and\\ncorresponding set in AS (NEBULA)") +
  scale_x_discrete(expand = c(0, 0.2), 
                   limits = c("ArtificialSel_logFC", "Plasticity_logFC"),
                   labels = c("Artificial Selection", "Plasticity")) +
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"))


#### Significant AS genes and corresponding plast genes, by logFC

react.norm.sigAS.data <- vc.sig.AS.genes[, c(3, 20, 21)]
colnames(react.norm.sigAS.data) <- c("ArtificialSel_logFC", "Gene", "Plasticity_logFC")
react.norm.sigAS.data <- gather(react.norm.sigAS.data, Experiment, logFC, c(ArtificialSel_logFC, Plasticity_logFC))

react.norm.sigAS <- ggplot(data = react.norm.sigAS.data, 
                           aes(x = Experiment, 
                               y = logFC)) + 
  geom_point(alpha = 0.3, size = 0.25) + 
  geom_line(aes(group = Gene), size = 0.25, alpha = 0.3) +
  theme_classic() +
  labs(title = "Significant AS genes (NEBULA) and\\ncorresponding set in plasticity (edgeR)") +
  scale_x_discrete(expand = c(0, 0.2), 
                   limits = c("ArtificialSel_logFC", "Plasticity_logFC"),
                   labels = c("Artificial Selection", "Plasticity")) +
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"))

react.norm.sigAS.abs <- ggplot(data = react.norm.sigAS.data, 
                               aes(x = Experiment, 
                                   y = abs(logFC))) + 
  geom_point(alpha = 0.3, size = 0.25) + 
  geom_line(aes(group = Gene), size = 0.25, alpha = 0.3) +
  theme_classic() +
  labs(title = "Significant AS genes (NEBULA) and\\ncorresponding set in plasticity (edgeR)") +
  scale_x_discrete(expand = c(0, 0.2), 
                   limits = c("ArtificialSel_logFC", "Plasticity_logFC"),
                   labels = c("Artificial Selection", "Plasticity")) +
  theme(axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"))

grid.arrange(arrangeGrob(react.norm.sigplast + theme(axis.title.x = element_blank()), 
                         react.norm.sigAS + theme(axis.title.x = element_blank()),
                         react.norm.sigplast.abs + theme(axis.title.x = element_blank(), plot.margin = unit(c(1,0,0.1,0), "cm")),
                         react.norm.sigAS.abs + theme(axis.title.x = element_blank(), plot.margin = unit(c(1,0,0.1,0), "cm")),
                         ncol = 2))


################################################################################
#### Expression reaction norm plot
#### 2 plots for each of Selection, Plasticity, comparing Down-Up, and Exp-Iso

# y.as.neb$E, y$E

# looking at all genes, so re-run voom w/o cutoff to get logCPM for all genes (for plast - this was already done for AS)

mm.2 <- model.matrix(~1 + Day + Vial_rep + Trt + Teneral_presence + Trt:Teneral_presence)
y.all.2 <- voom(d0, mm.2, plot = T) # d0 is the list without the cutoff

gene.counts.as <- as.data.frame(y.as.neb$E[, order(colnames(y.as.neb$E))])
gene.counts.plas <- as.data.frame(y.all.2$E[, order(colnames(y.all.2$E))])

gene.counts.as$Gene <- row.names(gene.counts.as)
gene.counts.as <- gene.counts.as %>% gather(Sample, Log2CPM, 1:54)
gene.counts.as$Treatment <- as.factor(substr(gene.counts.as$Sample, 1, 1))
gene.counts.as <- gene.counts.as[!gene.counts.as$Treatment == "C", ]
row.names(gene.counts.as) <- 1:length(row.names(gene.counts.as))

gene.counts.plas$Gene <- row.names(gene.counts.plas)
gene.counts.plas <- gene.counts.plas %>% gather(Sample, Log2CPM, 1:23)
gene.counts.plas$Treatment <- as.factor(substr(gene.counts.plas$Sample, 1, 1))
row.names(gene.counts.plas) <- 1:length(row.names(gene.counts.plas))

# summarize counts into average Log2CPM

summarized.gene.counts.as <- gene.counts.as %>% group_by(Gene, Treatment) %>% summarise(meanlog2CPM = mean(Log2CPM), .groups = "drop")
summarized.gene.counts.plast <- gene.counts.plas %>% group_by(Gene, Treatment) %>% summarise(meanlog2CPM = mean(Log2CPM), .groups = "drop")

# plot of ALL genes (not readable - do not run)
# expr.react.norm.as <- ggplot(data = summarized.gene.counts.as,
#                              aes(x = Treatment,
#                                  y = meanlog2CPM)) +
#   geom_point() +
#   geom_line(aes(group = Gene), stat = "identity",
#             size = 0.5, alpha = 0.5) + 
#   theme_classic()

summarized.sig.genes.as <- summarized.gene.counts.as[summarized.gene.counts.as$Gene %in% row.names(plast.hits.vc), ]

expr.react.norm.as <- ggplot(data = summarized.sig.genes.as,
                             aes(x = Treatment,
                                 y = meanlog2CPM)) +
  geom_line(aes(group = Gene), stat = "identity",
            size = 0.2, alpha = 0.2) + 
  theme_classic()

summarized.sig.genes.plast <- summarized.gene.counts.plast[summarized.gene.counts.plast$Gene %in% row.names(AS.trt.MD.hits), ]

expr.react.norm.plast <- ggplot(data = summarized.sig.genes.plast,
                                aes(x = Treatment,
                                    y = meanlog2CPM)) +
  geom_line(aes(group = Gene), stat = "identity",
            size = 0.2, alpha = 0.2) + 
  theme_classic()

grid.arrange(arrangeGrob(expr.react.norm.as + theme(axis.title.x = element_blank(), axis.title.y = element_blank()), 
                         expr.react.norm.plast + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                         left = "Mean log2 CPM",
                         bottom = "Treatment",
                         ncol = 2))

#################################################################
#################################################################
#### FIGURES FOR REVISION (Jan 2022)

### 1. Editor comment - figure showing fold change of all genes in both experiments

# Data frames with full sets of genes with logFC values for plasticity, AS (unshrunk and shrunk)

plast.trt.logFCcomp <- plast.trt.MD
AS.trt.logFCcomp <- AS.trt.MD
AS_apeglm.trt.logFCcomp <- apeglm.AS.MD

base::colnames(plast.trt.logFCcomp) <- base::paste(base::colnames(plast.trt.logFCcomp), "plast", sep = "_")
base::colnames(AS.trt.logFCcomp) <- base::paste(base::colnames(AS.trt.logFCcomp), "AS_not_shrunk", sep = "_")
base::colnames(AS_apeglm.trt.logFCcomp) <- base::paste(base::colnames(AS_apeglm.trt.logFCcomp), "AS_shrunk", sep = "_")

# as(unsrunk) vs. plasticity logFC

AS_unshrunk_vs_plast_logFC <- merge(plast.trt.logFCcomp, AS.trt.logFCcomp, by = "row.names", all.x = F)
row.names(AS_unshrunk_vs_plast_logFC) <- AS_unshrunk_vs_plast_logFC$Row.names
AS_unshrunk_vs_plast_logFC <- AS_unshrunk_vs_plast_logFC[, -1]

AS_unshrunk_vs_plast_logFC$Sig_plast <- recode_factor(AS_unshrunk_vs_plast_logFC$Sig_plast, Down = "SigPlast", Up = "SigPlast")
AS_unshrunk_vs_plast_logFC$Sig_AS_not_shrunk <- recode_factor(as.factor(AS_unshrunk_vs_plast_logFC$Sig_AS_not_shrunk), Down = "SigAS", Up = "SigAS")

AS_unshrunk_vs_plast_logFC$Sig_overall <- with(AS_unshrunk_vs_plast_logFC, interaction(Sig_plast, Sig_AS_not_shrunk))


AS_unshrunk_vs_plast_logFC_plot <- ggplot(data = AS_unshrunk_vs_plast_logFC %>% 
                                            arrange(match(Sig_overall, c("NotSig.NotSig", "SigPlast.NotSig", "NotSig.SigAS", "SigPlast.SigAS"))),
                                          aes(x = logFC_trt.nebU_AS_not_shrunk,
                                              y = logFC_plast,
                                              colour = Sig_overall)) +
  geom_point(stat = "identity", shape = 16, alpha = 0.8, stroke = 0, size = 1.5) +
  scale_y_continuous(limits = c(-3.2, 3.2), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  scale_x_continuous(limits = c(-3.2, 3.2), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  scale_color_manual(limits = c("NotSig.NotSig", "SigPlast.NotSig", "NotSig.SigAS", "SigPlast.SigAS"), 
                     values = c("black", "blue", "red", "green2")) +
  theme_classic() +
  labs(y = "log2 fold change (Plasticity, Isolated vs. Experienced)", 
       x = "log2 fold change (Artificial Selection, High vs. Low FC)") +
  theme(legend.position = "none", 
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 12, color = "black"))

# as(apeglm shrunk) vs. plasticity logFC

AS_shrunk_vs_plast_logFC <- merge(plast.trt.logFCcomp, AS_apeglm.trt.logFCcomp, by = "row.names", all.x = F)
row.names(AS_shrunk_vs_plast_logFC) <- AS_shrunk_vs_plast_logFC$Row.names
AS_shrunk_vs_plast_logFC <- AS_shrunk_vs_plast_logFC[, -1]

AS_shrunk_vs_plast_logFC_plot <- ggplot(data = AS_shrunk_vs_plast_logFC,
                                        aes(x = trt.UvD.logFC_AS_shrunk,
                                            y = logFC_plast)) +
  geom_point(stat = "identity", shape = 16, alpha = 0.2, stroke = 0, size = 1.5) +
  scale_y_continuous(limits = c(-3.2, 3.2), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  scale_x_continuous(limits = c(-3.2, 3.2), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  theme_classic() +
  labs(y = "log2 fold change (Plasticity)", x = "log2 fold change (AS APEGLM Shrunken)") +
  theme(legend.position = "none", 
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 12, color = "black"))


#### Teneral presence/absence effect - evolution vs plasticity

plast.ten.logFCcomp <- edger.ten$table
AS.ten.logFCcomp <- neb.mod.1$summary

base::colnames(plast.ten.logFCcomp) <- base::paste(base::colnames(plast.ten.logFCcomp), "plast", sep = "_")
base::colnames(AS.ten.logFCcomp) <- base::paste(base::colnames(AS.ten.logFCcomp), "AS_not_shrunk", sep = "_")

AS_unshrunk_vs_plast_ten_logFC <- merge(plast.ten.logFCcomp, AS.ten.logFCcomp, by = "row.names", all.x = F)
row.names(AS_unshrunk_vs_plast_ten_logFC) <- AS_unshrunk_vs_plast_ten_logFC$Row.names
AS_unshrunk_vs_plast_ten_logFC <- AS_unshrunk_vs_plast_ten_logFC[, -1]

tenpres.hits <- merge(neb.as.ten.hits, AS_unshrunk_vs_plast_ten_logFC, by = "row.names", all.x = F)
row.names(tenpres.hits) <- tenpres.hits$Row.names
tenpres.hits <- tenpres.hits[, -1]

AS_shrunk_vs_plast_logFC_ten_plot <- ggplot(data = AS_unshrunk_vs_plast_ten_logFC,
                                            aes(x = logFC_teneral.nebT_AS_not_shrunk,
                                                y = logFC_plast)) +
  geom_point(stat = "identity", shape = 16, alpha = 0.8, stroke = 0, size = 1.5) +
  geom_point(data = tenpres.hits, aes(x = logFC_teneral.nebT, y = logFC_plast), stat = "identity", shape = 16, alpha = 0.8, stroke = 0, size = 1.5, color = "red") +
  scale_y_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  scale_x_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  theme_classic() +
  labs(y = "log2 fold change (Plasticity, Teneral exposed vs. not exposed)", 
       x = "log2 fold change (Artificial Selection, Teneral exposed vs. not exposed)") +
  theme(legend.position = "none", 
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 12, color = "black"))

grid.arrange(arrangeGrob(AS_unshrunk_vs_plast_logFC_plot, 
                         AS_shrunk_vs_plast_logFC_ten_plot,
                         ncol = 2))



################################################################
##### Reviewer comment - testing how many genes we would expect to get overlapping just by chance

plast.rand <- row.names(as.data.frame(cts.dgelist$counts))
as.rand <- row.names(as.data.frame(neb.mod.1$summary))

overlap.samples <- c()

for (i in 1:10000) {
  plast.sample <- sample(plast.rand, size = 375, replace = FALSE)
  as.sample <- sample(as.rand, size = 1030, replace = FALSE)
  overlap.num <- length(base::intersect(plast.sample, as.sample))
  overlap.samples <- rbind(overlap.samples, overlap.num)
}

quantile(overlap.samples[,1], probs = c(0.025, 0.5, 0.975))

