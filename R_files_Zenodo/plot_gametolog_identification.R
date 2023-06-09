#load the required libraries
library(dplyr)
library(ggplot2)
library(see)

#manhattan plot genomewide using ggplot2 (code adapted from https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/)
data_GWAS<-read.csv(file="FL_Cuba.F.M.gwas_for_plotting.txt", header=TRUE, sep=' ')
data_chr7<-data_GWAS[data_GWAS$CHR == '7',]

#manhattan plot genomewide
nCHR <- length(unique(data_GWAS$CHR))
data_GWAS$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(data_GWAS$CHR)){
  nbp[i] <- max(data_GWAS[data_GWAS$CHR == i,]$BP)
  data_GWAS[data_GWAS$CHR == i,"BPcum"] <- data_GWAS[data_GWAS$CHR == i,"BP"] + s
  s <- s + nbp[i]
}
axis.set <- data_GWAS %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- abs(floor(log10(min(data_GWAS$P)))) + 3
sig <- 0.05 / nrow(data_GWAS)

ggplot(data_GWAS, aes(x=BPcum, y=-log10(P))) +
        geom_point(aes(color=as.factor(CHR)), alpha = 0.75, size = 1.25) +
  scale_color_manual(values = rep(c("#a32126ff","#e09b45ff"), nCHR)) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0,ylim)) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  labs(x = NULL, y = "-log10(p)") + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, vjust = 0.5))

#manhattan plot for CHR7 only using ggplot2
ggplot(data_chr7, aes(x=BP, y=-log10(P))) +
  geom_point(colour="#a32126ff", alpha = 0.75, size = 1.25) +
  scale_y_continuous(expand = c(0,0), limits = c(0,ylim)) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  labs(x = "\\nChromosome 7 (Mbp)", y = "-log10(p)") + 
  geom_vline(xintercept = 21000000)+
  geom_vline(xintercept = 92000000)+
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())


##dotplot and violin plot for gametolog SNP heterozygosity (separate for females and males)
heterozygosity<-read.csv(file="FM.het_rate.txt", header=TRUE, sep=' ')

ggplot(heterozygosity, aes(x=Sex, y=Heterozygosity)) +
  geom_violindot(size_dots = 0.3,binwidth = 0.04)+
  labs(x = "\\nSex", y = "Heterozygosity (%)\\n") + 
  theme_bw()+
  theme(axis.text =element_text(size=12),
        axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank())

#test for heterozygosity differences between females and males at gametolog SNPs, using a Wilcoxon test
females <- heterozygosity[heterozygosity$Sex == 'F',]
het_females<-females$Heterozygosity

males <- heterozygosity[heterozygosity$Sex == 'M',]
het_males<-males$Heterozygosity

wilcox.test(het_females, het_males)

#plot results of PCA based on gametolog SNPs (used to confirm the sex of genotyped samples)
#the input file was obtained using the "adegenet_PCA.R" script; before plotting, we added the field-based sex assignment for each sample

pca_all<-read.csv(file="chr7_gametolog_snps.PCA_results.csv",header = TRUE, sep = ',')
library(ggplot2)
ggplot(pca_all, aes(PC1, PC2))+
  geom_point(aes(colour=Field_sex_assignment), size=2)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_colour_manual(values = c("black","lightgrey","firebrick"))
