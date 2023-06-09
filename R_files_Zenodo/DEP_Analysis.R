# Load required packages 

library(DEP)
library(tidyverse)

# Load dataset - MaxQuant output

mQ.proteins <- read.table("combined/txt/proteinGroups.txt", 
           sep = "\\t", header = TRUE, stringsAsFactors = FALSE)

# Establish experimental design

design <- data.frame(label = c("DMSO_R1", "DMSO_R2", "DMSO_R3",
                               "UNC6934_R1", "UNC6934_R2", "UNC6934_R3",
                               "UNC7145_R1", "UNC7145_R2", "UNC7145_R3"),
           condition = c("DMSO", "DMSO", "DMSO", 
                         "UNC6934", "UNC6934", "UNC6934",  
                         "UNC7145", "UNC7145", "UNC7145"),
           replicate = c(1,2,3,1,2,3,1,2,3), stringsAsFactors = FALSE)

# Create input dataset for DEP differential analysis using import_MaxQuant function 

data_se <- import_MaxQuant(proteins = mQ.proteins,
                expdesign = design,
                filter = c("Reverse", "Potential.contaminant"),
                names = "Gene.names",
                ids = "Protein.IDs") 

# Protein frequency plot - most protein identifications are found in only one sample, indicating most signal is likely background. 

plot_frequency(data_se)


# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 1)

# Plot number of proteins identified per sample & independent experiment - Generally, fewer proteins were capture for biological replicate 1. 

plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data & plot
data_norm <- normalize_vsn(data_filt)

plot_normalization(data_se, data_norm)

plot_missval(data_se)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)


# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_norm, type = "control", control = "DMSO")

# Denote significant proteins based on user defined cutoffs

dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(2))


# Plot a volcano plot for the contrast "UNC6934 vs DMSO""

plot_volcano(dep, 
             contrast = "UNC6934_vs_DMSO", 
             label_size = 2, 
             add_names = TRUE, 
             plot = FALSE, adjusted = FALSE) -> UNC6934

UNC6934 %>% 
  mutate(labels = ifelse(significant == "TRUE", as.character(protein), "")) %>% 
  ggplot(aes(log2_fold_change, `p_value_-log10`, 
             fill = significant, label = labels, size = significant, 
             alpha = significant)) +
  geom_point(pch =21) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = -1.5, linetype = 2) +
  geom_vline(xintercept = 0) +
  scale_alpha_discrete(range = c("0.5", "1")) +
  scale_size_discrete(range = c(2, 3)) +
  scale_fill_manual(values = c("grey", "red")) +
  ggrepel::geom_text_repel(point.padding = 0.5, size =3) +
  theme_classic() +
  ylab("-log10(pvalue)") +
  xlab("log2(Fold Change)") +
  theme(legend.position = "none")
  

plot_single(dep, proteins = c("WHSC1"), type = 'centered')

# Plot a volcano plot for bot UNC7145 & UNC6934

plot_volcano(dep, 
             contrast = "UNC7145_vs_DMSO", 
             label_size = 2, 
             add_names = TRUE,
             plot = FALSE, adjusted = FALSE) -> UNC7145

UNC6934 %>% mutate(Sample = "UNC6934") -> UNC6934

UNC7145 %>%  mutate(Sample = "UNC7145") %>% 
  bind_rows(UNC6934) %>% 
  mutate(labels = ifelse(protein =="WHSC1" & Sample == "UNC6934", 
                         "NSD2", "")) %>% 
  mutate(hit = ifelse(significant == TRUE & log2_fold_change < -1.5,
                      "Yes", "No")) %>% 
  ggplot(aes(log2_fold_change, `p_value_-log10`, 
             fill = Sample, label = labels, size = hit, 
             alpha = hit)) + 
  geom_rect(aes(xmin = -Inf, xmax = -1.5, 
                ymin =  1.3010, ymax = Inf), 
            fill = "#f2f2f2") +
  geom_point(pch =21) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = -1.5, linetype = 2) +
  geom_vline(xintercept = 0) +
  scale_alpha_discrete(range = c("0.5", "1")) +
  scale_size_discrete(range = c(2, 3)) +
  scale_fill_manual(values = c("blue", "red")) +
  ggrepel::geom_text_repel(point.padding = 0.5, size =3) +
  theme_classic() +
  ylab("-log10(pvalue)") +
  xlab("log2(Fold Change)") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 
  
