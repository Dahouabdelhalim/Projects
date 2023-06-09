library(tidyverse)
library(edgeR)
library(limma)
library(splines)
library(ggExtra)
library(nlme)
library(patchwork)

#load gene counts file
counts <- read.delim("submission_files/RNAseq_all_counts.txt", row.names=1)

#create object for limma and edgeR
d <- DGEList(counts)

#load metadata
metadat <- read.table("resubmission_files/sample_metadata_v2.txt", header=T, sep='\\t') %>%
  drop_na(novogene_id_RNAseq) %>%
  arrange(novogene_id_RNAseq,rownames(d$samples)) #sort by DGEList object

#confirm that sample order is the same btwn counts file and metadata
identical(rownames(d$samples),metadat$novogene_id_RNAseq)
#add age & colony metadata to the DGEList object
d$samples$colony <- metadat$Colony
d$samples$age <- metadat$Age_days

# transform counts
cpm <- cpm(d)
lcpm <- cpm(d, log=TRUE) #default offset added prior to taking log is 2/L where L = ave library size in M

#filter genes 
keep.exprs <- filterByExpr(d)
d.filt <- d[keep.exprs, keep.lib.sizes=FALSE]
dim(d); dim(d.filt)
par(mfrow=c(1,2))
voom(d, plot=T)
voom(d.filt, plot=T)

#ensure that expression distribution (histogram of log-cpm) is similar across samples
L <- mean(d$samples$lib.size) * 1e-6
M <- median(d$samples$lib.size) * 1e-6
L; M
lcpm.cutoff <- log2(10/M + 2/L)

nsamples <- ncol(d)
par(mfrow=c(1,2))
#plot raw data
plot(density(lcpm[,1]), lwd=2, ylim=c(0,0.26), las=2, main="", xlab=""); title(main="Raw data", xlab="Log-cpm"); abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, lwd=2, col="dodgerblue")
}
#plot filtered data
lcpm.filt <- cpm(d.filt, log=TRUE)
plot(density(lcpm.filt[,1]), lwd=2, ylim=c(0,0.26), las=2, main="", xlab=""); title(main="Filtered data", xlab="Log-cpm"); abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm.filt[,i])
  lines(den$x, den$y, lwd=2, col="maroon")
}

# calculate normalization factors by TMM
d.norm <- calcNormFactors(d.filt, method = "TMM")
d.norm$samples$norm.factors # closer to 1.0 means less effect of normalization

par(mfrow=c(1,2))
#plot unnormalized
boxplot(lcpm.filt, las=2, main="")
title(main="Unnormalized data", ylab="Log-cpm")
#plot normalized
lcpm.norm <- cpm(d.norm, log=TRUE)
boxplot(lcpm.norm, las=2, main="")
title(main="Normalized data", ylab="Log-cpm")



##### visualize overall gene expression profile similarity
#####

#PCoA: manually extract coordinates to plot in ggplot
#by default uses  gene.selection = "pairwise"
expr.plot <- plotMDS(d.norm, col = as.numeric(d.norm$samples$colony))

dim1 <- expr.plot$x
dim2 <- expr.plot$y
df_MDS <- data.frame(dim1,dim2) %>%
  mutate(novogene_id_RNAseq = rownames(d.norm$samples)) %>%
  full_join(metadat, by='novogene_id_RNAseq')

#add a discrete age scale
#groups of appx equal number of bees each
discrete_ages <- cut(df_MDS$Age_days, breaks=c(-Inf,1,19,43, Inf),
                     labels=c("new","young","middle","old"))
table(discrete_ages)
df_MDS$discrete_ages <- discrete_ages

#age distribution by category
dplyr::select(df_MDS, c(discrete_ages,Age_days)) %>% arrange(Age_days) %>% distinct()

metadat$discrete_ages <- discrete_ages
d.norm$samples$discrete_ages <- discrete_ages

#invert along x axis to match the young -> old direction of the 16S ordination
df_MDS$dim1_inv <- -1*df_MDS$dim1

plot_by_age_discrete <- ggplot(df_MDS, aes(x=dim1_inv, y=dim2, color=discrete_ages, fill=discrete_ages)) +
  geom_point(size=3, pch=21, alpha=0.8) +
  scale_color_manual(values=c("new" = "#762A83",  #matches the strain network colors
                             "young" = "#E7D4E8",
                             "middle"="#D9F0D3",
                             "old"="#1B7837"),
                     name="Bee age (days)") +
  scale_fill_manual(values=c("new" = "#762A83", 
                              "young" = "#E7D4E8",
                              "middle"="#D9F0D3",
                              "old"="#1B7837"),
                     name="Bee age (days)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  xlab(paste(expr.plot$axislabel," 1 (",100*round(expr.plot$var.explained[1],2),"%)",sep="")) +
  ylab(paste(expr.plot$axislabel," 2 (",100*round(expr.plot$var.explained[2],2),"%)",sep=""))
plot_by_age_discrete



## continue with analyses using discrete age groups using edgeR

colony <- as.factor(d.norm$samples$colony)

#creating the design matrix
design1 <- model.matrix(~0+discrete_ages+ colony) #removes intercept from ages
colnames(design1) <- gsub("discrete_ages", "", colnames(design1))

#set up contrasts using limma function
contr.matrix <- makeContrasts(
  new_v_young = new-young,
  young_v_middle = young-middle,
  middle_v_old = middle-old,
  levels = colnames(design1))
contr.matrix

#remove heteroscedascity and fit linear models
v <- voom(d.norm, design1, plot=T) # plot variance ~ mean before applying voom precision weighting
vfit <- lmFit(v, design1)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)

#check number of DEGs for contrasts X_v_Y
#"down" or "up"  = down/upregulated in X (left) relative to Y (right)
dt <- decideTests(efit)
summary(dt) #by default, adjusted p < 0.05

plotMD(efit, column=1, status=dt[,1], main=colnames(efit)[1], xlim=c(-8,13)) #new vs young
plotMD(efit, column=2, status=dt[,2], main=colnames(efit)[2], xlim=c(-8,13)) #young vs middle
plotMD(efit, column=1, status=dt[,3], main=colnames(efit)[3], xlim=c(-8,13)) #middle vs old

#ranked DEGs
rankedDEGs <- topTable(efit, n = Inf) #all genes sorted by P
rankedDEGs$Gene <- rownames(rankedDEGs)
rankedDEGs$rank <- 1:length(rankedDEGs[,1])


## Gene annotation data

#downloaded all B. imp gene IDs (entrez ids), to get gene names in EntrezDirect
genenames <- read.table("resubmission_files/AllGenesEntrezIDs.GeneNames.txt",sep='\\t',
                        quote = "", 
                        row.names = NULL, 
                        stringsAsFactors = FALSE)
names(genenames) <- c("entrezID","Gene","genename")

#merge w/ DEGs
rankedDEGs_w_genename <- full_join(rankedDEGs, genenames, by="Gene")
#merge with contrasts DE stats
contrast_DE <- as.data.frame(dt)
colnames(contrast_DE) <- gsub("^","sig_",colnames(contrast_DE))
contrast_DE$Gene <- rownames(contrast_DE)
rankedDEGs_w_genename <- full_join(rankedDEGs_w_genename, contrast_DE, by="Gene")

#select candidate immunity effector genes
imm_genes_list <- c("abaecin","apidaecin","defensin","hymenoptaecin","dual oxidase","catalase")
imm_genes <- genenames %>% filter(str_detect(genename, imm_genes_list[1]) |
                                    str_detect(genename, imm_genes_list[2]) |
                                    str_detect(genename, imm_genes_list[3]) |
                                    str_detect(genename, imm_genes_list[4]) |
                                    str_detect(genename, imm_genes_list[5]) |
                                    str_detect(genename, imm_genes_list[6]))


### plot expression (logCPM) of imm genes ~ discrete age categories (Fig. 6B)

#prep metadata and gene name data
d.norm$samples$Sample <- rownames(d.norm$samples)
imm_genes_subset <- imm_genes %>%
  mutate(genename = gsub("s type 73","",genename)) %>%
  mutate(genename = gsub("-1","",genename)) %>%
  mutate(genename = gsub("-like","",genename)) %>%
  filter(genename != "dual oxidase maturation factor 1")

lcpm.df.imm.genes <- as.data.frame(lcpm.norm) %>%
  rownames_to_column(var="Gene") %>% 
  pivot_longer(-Gene,names_to = "Sample", values_to = "logCPM") %>%
  inner_join(d.norm$samples, by="Sample") %>%
  inner_join(imm_genes_subset, by="Gene")

lcpm.df.imm.genes$genename <- factor(lcpm.df.imm.genes$genename, c("apidaecin",
                                                                   "abaecin",
                                                                   "dual oxidase",
                                                                   "hymenoptaecin",
                                                                   "defensin",
                                                                   "catalase"))

ggplot(lcpm.df.imm.genes, aes(x=discrete_ages,y=logCPM)) +
  facet_wrap(~genename, scales = "free_y") +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("")

# to show significant pairwise contrasts:
rankedDEGs_w_genename %>% filter(Gene %in% imm_genes_subset$Gene)
# manually add as dashed lines to plot





###### IMD and Toll pathway analyses

#select genes
IMD_Toll_list <- c("cactus",
                   "embryonic polarity protein dorsal","immune deficiency",
                   "nuclear factor NF-kappa-B p100 subunit")
IMD_Toll_genes <- genenames %>% filter(str_detect(genename, IMD_Toll_list[1]) |
                                      str_detect(genename, IMD_Toll_list[2]) |
                                      str_detect(genename, IMD_Toll_list[3]) |
                                      str_detect(genename, IMD_Toll_list[4]))

### plot expression (logCPM) of IMD/Toll genes ~ discrete age categories (Fig. S10)

#prep metadata and gene name data
lcpm.df.IMDToll.genes <- as.data.frame(lcpm.norm) %>%
  rownames_to_column(var="Gene") %>% 
  pivot_longer(-Gene,names_to = "Sample", values_to = "logCPM") %>%
  inner_join(d.norm$samples, by="Sample") %>%
  inner_join(IMD_Toll_genes, by="Gene")

ggplot(lcpm.df.IMDToll.genes, aes(x=discrete_ages,y=logCPM)) +
  facet_wrap(~genename, scales = "free_y") +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("")

# to indicate significant contrasts:
rankedDEGs_w_genename %>% filter(Gene %in% IMD_Toll_genes$Gene)
# manually add to plot as dashed lines