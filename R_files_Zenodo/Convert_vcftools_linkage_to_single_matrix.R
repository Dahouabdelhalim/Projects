#!/usr/bin/env Rscript

#####
# C. Tepolt, WHOI, 21 October 2020
# Convert pairwise LD values from vcftools to matrix
#####

library(data.table)
library(reshape2)
library(stringr)

options(scipen = 999)

within = fread('[FILENAME].out.geno.ld')
between = fread('[FILENAME].out.interchrom.geno.ld')

# If you want to omit some SNPs from the table, specify the file here
# File format should be:
# Contig1_SNP1
# Contig1_SNP2
# Contig2_SNP1
#to_remove = scan("[LIST_TO_REMOVE].txt", what = list(SNP = character()))

w_new = cbind(paste(within$CHR, within$POS1, sep='_'), paste(within$CHR, within$POS2, sep='_'), within$'R^2', str_extract(within$CHR, "[0-9]+"), within$POS1, str_extract(within$CHR, "[0-9]+"), within$POS2)
b_new = cbind(paste(between$CHR1, between$POS1, sep='_'), paste(between$CHR2, between$POS2, sep='_'), between$'R^2', str_extract(between$CHR1, "[0-9]+"), between$POS1, str_extract(between$CHR2, "[0-9]+"), between$POS2)

zaphod = rbind(w_new, b_new)
colnames(zaphod) = c("SNP1","SNP2","R2", "CHR1", "POS1", "CHR2", "POS2")
zaphod = as.data.table(zaphod)
zaphod = zaphod[order(as.numeric(CHR1), as.numeric(POS1), as.numeric(CHR2), as.numeric(POS2))]

#zaphod = zaphod[!(SNP1 %in% to_remove$SNP) & !(SNP2 %in% to_remove$SNP)]

zaphod_df = data.frame(zaphod[,0:3])

SNP_order = c(unique(zaphod$SNP1), tail(zaphod$SNP2, n = 1))

zaphod_df$SNP1 = factor(zaphod_df$SNP1, SNP_order)
zaphod_df$SNP2 = factor(zaphod_df$SNP2, SNP_order)

#if (sum(is.na(zaphod_df$R2 == "NA")) > 0)
#{
# zaphod_df[zaphod_df$R2 == "NA",]$R2 = 0
#}
#if (sum(is.nan(zaphod_df$R2 == "NaN")) > 0)
#{
# zaphod_df[zaphod_df$R2 == "NaN",]$R2 = 0
#}

mat = dcast(zaphod_df, SNP1~SNP2, drop = FALSE)
rownames(mat) = mat[,1]
mat$SNP1 = NULL

write.csv(mat, file = "LD_matrix.csv")
