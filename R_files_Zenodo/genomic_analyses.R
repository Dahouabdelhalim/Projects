###Analyze FST data from sympatric benthic-limnetic species pairs and allopatric populations from small and large lakes, GWAS data for 28 allopatric populations 

library(data.table)
library(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(qqman)
library(reshape2)
library(qqman)
library(purrr)
library(utils)

#Load small-large data with all FST values
fst_small_large <- read.csv("fst_small_large_lakes.csv", sep = ";")

#Define 95th percentile
norm_allo_95 <- quantile(fst_small_large$mean_norm, na.rm = TRUE, probs = 0.95)
lowcamp_lilmud_95 <- quantile(fst_small_large$lowcamp_lilmud, na.rm = TRUE, probs = 0.95)
uppcamp_ormund_95 <- quantile(fst_small_large$uppcamp_ormund, na.rm = TRUE, probs = 0.95)
stella_lilgoose_95 <- quantile(fst_small_large$stella_lilgoose, na.rm = TRUE, probs = 0.95)

#FST outliers across all small-large lake population pairs
small_large_outliers <- fst_small_large$lg_bp_group[which(fst_small_large$mean_norm >= norm_allo_95)]
small_large_outliers <- as.character(small_large_outliers)

#Calculate mean genome-wide FST for each population pair
mean(fst_small_large$lowcamp_lilmud, na.rm = TRUE)
mean(fst_small_large$uppcamp_ormund, na.rm = TRUE)
mean(fst_small_large$stella_lilgoose, na.rm = TRUE)

#Calculate mean genome-wide FST for each populations pair, excluding FST outlier (95th percentile)
lowcamp_lilmud_no_outlier <- fst_small_large$lowcamp_lilmud[which(fst_small_large$lowcamp_lilmud < lowcamp_lilmud_95)]
mean(lowcamp_lilmud_no_outlier, na.rm = TRUE)

uppcamp_ormund_no_outlier <- fst_small_large$uppcamp_ormund[which(fst_small_large$uppcamp_ormund < uppcamp_ormund_95)]
mean(uppcamp_ormund_no_outlier, na.rm = TRUE)

stella_lilgoose_no_outlier <- fst_small_large$stella_lilgoose[which(fst_small_large$stella_lilgoose < stella_lilgoose_95)]
mean(stella_lilgoose_no_outlier, na.rm = TRUE)

#Load sympatric species pair data with mean FST values for each window
fst_species_pairs <- read.csv("fst_species_pairs.csv", sep = ";")

#Define 95th percentile
norm_sym_95 <- quantile(fst_species_pairs$mean_norm, na.rm = TRUE, probs = 0.95)
pax_95 <- quantile(fst_species_pairs$pax, na.rm = TRUE, probs = 0.95)
pri_95 <- quantile(fst_species_pairs$pri, na.rm = TRUE, probs = 0.95)
lq_95 <- quantile(fst_species_pairs$lq, na.rm = TRUE, probs = 0.95)

#FST outliers across all benthic-limnetic species pairs
ben_lim_outliers <- fst_species_pairs$lg_bp_group[which(fst_species_pairs$mean_norm >= norm_sym_95)]
ben_lim_outliers <- as.character(ben_lim_outliers)

#Calculate mean genome-wide FST for each species pair
mean(fst_species_pairs$pax, na.rm = TRUE)
mean(fst_species_pairs$pri, na.rm = TRUE)
mean(fst_species_pairs$lq, na.rm = TRUE)

#Calculate mean genome-wide FST for each species pair, excluding FST outlier (95th percentile)
pax_no_outlier <- fst_species_pairs$pax[which(fst_species_pairs$pax < pax_95)]
mean(pax_no_outlier, na.rm = TRUE)

pri_no_outlier <- fst_species_pairs$pri[which(fst_species_pairs$pri < pri_95)]
mean(pri_no_outlier, na.rm = TRUE)

lq_no_outlier <- fst_species_pairs$lq[which(fst_species_pairs$lq < lq_95)]
mean(lq_no_outlier, na.rm = TRUE)

##Plot FST correlation across lakes within each data set
#Small-large lakes
plot(fst_small_large$stella_lilgoose~fst_small_large$lowcamp_lilmud, pch = 20, col = "grey58", ylim = c(0,1), xlim = c(0,1))
abline(lm(fst_small_large$stella_lilgoose~fst_small_large$lowcamp_lilmud))

plot(fst_small_large$stella_lilgoose~fst_small_large$uppcamp_ormund, pch = 20, col = "grey58", ylim = c(0,1), xlim = c(0,1))
abline(lm(fst_small_large$stella_lilgoose~fst_small_large$uppcamp_ormund))

plot(fst_small_large$uppcamp_ormund~fst_small_large$lowcamp_lilmud, pch = 20, col = "grey58", ylim = c(0,1), xlim = c(0,1))
abline(lm(fst_small_large$uppcamp_ormund~fst_small_large$lowcamp_lilmud))

#Benthic-limnetic species pairs
plot(fst_species_pairs$pax ~ fst_species_pairs$pri, main = "Paxton-Priest", pch = 20, col = "grey58", ylim = c(0,1))
abline(lm(fst_species_pairs$pax ~ fst_species_pairs$pri))

plot(fst_species_pairs$pax ~ fst_species_pairs$lq, main = "Paxton-Little Quarry", pch = 20, col = "grey58", ylim = c(0,1))
abline(lm(fst_species_pairs$pax ~ fst_species_pairs$lq))

plot(fst_species_pairs$pri ~ fst_species_pairs$lq, main = "Priest-Little Quarry", pch = 20, col = "grey58", ylim = c(0,1))
abline(lm(fst_species_pairs$pri ~ fst_species_pairs$lq))

#Calculate Spearman correlations of FST values (for 50-kb windows)

#Benthic-limnetic species pairs
spearman_rho1 <- NA
spearman_p1 <- NA
comp_names1 <- NA

ben_lim_names <- colnames(fst_species_pairs)[4:6]

for (i in 1:2){
  count <- 3-i
  for (g in 1:count){
    cor_test <- cor.test(fst_species_pairs[,i+3], fst_species_pairs[,i+3+g], method = "spearman")
    spearman_rho1 <- c(spearman_rho1, cor_test$estimate)
    spearman_p1 <- c(spearman_p1, cor_test$p.value)
    name <- sprintf("%s_%s", ben_lim_names[i], ben_lim_names[i+g])
    comp_names1 <- c(comp_names1, name)
  }
}

spearman_rho1 <- spearman_rho1[-1]
spearman_p1 <- spearman_p1[-1]
comp_names1 <- comp_names1[-1]

#Small-large lakes
spearman_rho2 <- NA
spearman_p2 <- NA
comp_names2 <- NA

small_large_names <- colnames(fst_small_large)[4:6]

for (i in 1:2){
  count <- 3-i
  for (g in 1:count){
    cor_test <- cor.test(fst_small_large[,i+3], fst_small_large[,i+3+g], method = "spearman")
    spearman_rho2 <- c(spearman_rho2, cor_test$estimate)
    spearman_p2 <- c(spearman_p2, cor_test$p.value)
    name <- sprintf("%s_%s", small_large_names[i], small_large_names[i+g])
    comp_names2 <- c(comp_names2, name)
  }
}

spearman_rho2 <- spearman_rho2[-1]
spearman_p2 <- spearman_p2[-1]
comp_names2 <- comp_names2[-1]

#Create dataframes with correlation coefficients and p-values
benlim_spearman <- data.frame(data.frame(matrix(NA, nrow = 3),
                                              comp = comp_names1,
                                              spearman_rho = spearman_rho1,
                                              spearman_p = spearman_p1,
                                              cat = "sym",
                                              stringsAsFactors=FALSE))

benlim_spearman <- benlim_spearman[,-1]

allo_spearman <- data.frame(data.frame(matrix(NA, nrow = 3),
                                         comp = comp_names2,
                                         spearman_rho = spearman_rho2,
                                         spearman_p = spearman_p2,
                                         cat = "allo",
                                         stringsAsFactors=FALSE))

allo_spearman <- allo_spearman[,-1]

#Combine dataframes
fst_spearman <- rbind(benlim_spearman, allo_spearman)

#Calculate correlations across data sets (allopatric vs. sympatric)

fst_sym_allo <- fst_small_large[which(fst_small_large$lg_bp_group %in% fst_species_pairs$lg_bp_group),]
fst_sym_allo$pax <- fst_species_pairs$pax[match(fst_sym_allo$lg_bp_group, fst_species_pairs$lg_bp_group)]
fst_sym_allo$pri <- fst_species_pairs$pri[match(fst_sym_allo$lg_bp_group, fst_species_pairs$lg_bp_group)]
fst_sym_allo$lq <- fst_species_pairs$lq[match(fst_sym_allo$lg_bp_group, fst_species_pairs$lg_bp_group)]

spearman_rho <- NA
spearman_p <- NA
comp_names <- NA

for (i in 1:length(small_large_names)){
  cor_test1 <- cor.test(fst_sym_allo[,i+3], fst_sym_allo$pax, method = "spearman")
  cor_test2 <- cor.test(fst_sym_allo[,i+3], fst_sym_allo$pri, method = "spearman")
  cor_test3 <- cor.test(fst_sym_allo[,i+3], fst_sym_allo$lq, method = "spearman")
  spearman_rho <- c(spearman_rho, cor_test1$estimate, cor_test2$estimate, cor_test3$estimate)
  spearman_p <- c(spearman_p, cor_test1$p.value, cor_test2$p.value, cor_test3$p.value)
  name1 <- sprintf("%s_%s", small_large_names[i], "pax")
  name2 <- sprintf("%s_%s", small_large_names[i], "pri")
  name3 <- sprintf("%s_%s", small_large_names[i], "lq")
  comp_names <- c(comp_names, name1, name2, name3)
}

spearman_rho <- spearman_rho[-1]
spearman_p <- spearman_p[-1]
comp_names <- comp_names[-1]

spearman_rho <- as.numeric(spearman_rho)

sym_allo_spearman <- data.frame(data.frame(matrix(NA, nrow = 9),
                                           comp = comp_names,
                                           spearman_rho = spearman_rho,
                                           spearman_p = spearman_p,
                                           cat = "sym_allo",
                                           stringsAsFactors=FALSE))

sym_allo_spearman <- sym_allo_spearman[,-1]

spearman_combined <- rbind(fst_spearman, sym_allo_spearman)
spearman_combined$spearman_rho <- as.numeric(spearman_combined$spearman_rho)
spearman_combined$cat <- as.factor(spearman_combined$cat)

plot(spearman_combined$spearman_rho ~ spearman_combined$cat, pch = 16, ylim = c(-0.1,0.6))
points(spearman_combined$spearman_rho ~ spearman_combined$cat, cex = 0.7, pch = 16)

#Calculate mean correlation coefficient for each category
mean(spearman_combined$spearman_rho[which(spearman_combined$cat == "allo")], na.rm = TRUE)
mean(spearman_combined$spearman_rho[which(spearman_combined$cat == "sym")], na.rm = TRUE)
mean(spearman_combined$spearman_rho[which(spearman_combined$cat == "sym_allo")], na.rm = TRUE)

#Parallelism estimate
#Determine repeated outlier windows

#Allopatric population pairs
fst_small_large$lowcamp_lilmud_sig <- fst_small_large$lowcamp_lilmud >= lowcamp_lilmud_95
fst_small_large$uppcamp_ormund_sig <- fst_small_large$uppcamp_ormund >= uppcamp_ormund_95
fst_small_large$stella_lilgoose_sig <- fst_small_large$stella_lilgoose >= stella_lilgoose_95

fst_small_large$overalltrue <- rowSums(fst_small_large =="TRUE", na.rm=TRUE)
fst_small_large$overallfalse <- rowSums(fst_small_large =="FALSE", na.rm=TRUE)
fst_small_large$overall_total <- fst_small_large$overallfalse + fst_small_large$overalltrue

fst_small_large$overallparallel <- (((fst_small_large$overalltrue) * (fst_small_large$overalltrue -1))/((fst_small_large$overall_total)*(fst_small_large$overall_total - 1)))

fst_parallel_small_large <- fst_small_large$lg_bp_group[which(fst_small_large$overallparallel == 1)]
fst_parallel_small_large <- as.character(fst_parallel_small_large)

#Sympatric species pairs
fst_species_pairs$pax_sig <- fst_species_pairs$pax >= pax_95
fst_species_pairs$pri_sig <- fst_species_pairs$pri >= pri_95
fst_species_pairs$lq_sig <- fst_species_pairs$lq >= lq_95

fst_species_pairs$overalltrue <- rowSums(fst_species_pairs =="TRUE",na.rm=TRUE)
fst_species_pairs$overallfalse <- rowSums(fst_species_pairs =="FALSE",na.rm=TRUE)
fst_species_pairs$overall_total <- fst_species_pairs$overallfalse + fst_species_pairs$overalltrue

fst_species_pairs$overallparallel <- ((fst_species_pairs$overalltrue) * (fst_species_pairs$overalltrue -1))/((fst_species_pairs$overall_total)*(fst_species_pairs$overall_total - 1))

fst_parallel_benlim <- fst_species_pairs$lg_bp_group[which(fst_species_pairs$overallparallel == 1)]
fst_parallel_benlim <- as.character(fst_parallel_benlim)

#Significant repeated windows shared across both datasets
fst_parallel_combined <- fst_parallel_benlim[which(fst_parallel_benlim %in% fst_parallel_small_large)]

#Significant outlier windows based on mean_norm FST across both datasets
fst_outlier_combined <- ben_lim_outliers[which(ben_lim_outliers %in% small_large_outliers)]

##Manhattan plots

#Small-large lakes
manhattan(fst_small_large, chr="lg", bp="bp_group", snp="lg_bp_group", p= "lowcamp_lilmud", logp = FALSE, suggestiveline = FALSE, genomewideline = lowcamp_lilmud_95, highlight = fst_parallel_small_large, main= "Lower Campbell-Lil Mud" , ylim = c(0,1))
manhattan(fst_small_large, chr="lg", bp="bp_group", snp="lg_bp_group", p= "uppcamp_ormund", logp = FALSE, suggestiveline = FALSE, genomewideline = uppcamp_ormund_95, highlight = fst_parallel_small_large, main= "Upper Campbell-Ormund" , ylim = c(0,1))
manhattan(fst_small_large, chr="lg", bp="bp_group", snp="lg_bp_group", p= "stella_lilgoose", logp = FALSE, suggestiveline = FALSE, genomewideline = stella_lilgoose_95, highlight = fst_parallel_small_large, main= "Stella-Lil Goose" , ylim = c(0,1))

#Benthic-limnetic species pairs
manhattan(fst_species_pairs, chr="lg", bp="bp_group", snp="lg_bp_group", p="pax", logp = FALSE, suggestiveline = FALSE, genomewideline = pax_95 , highlight = fst_parallel_benlim, main= "Paxton_ben-lim" , ylim = c(0,1))
manhattan(fst_species_pairs, chr="lg", bp="bp_group", snp="lg_bp_group", p="pri", logp = FALSE, suggestiveline = FALSE, genomewideline = pri_95 , highlight = fst_parallel_benlim, main= "Priest_ben-lim" , ylim = c(0,1))
manhattan(fst_species_pairs, chr="lg", bp="bp_group", snp="lg_bp_group", p="lq", logp = FALSE, suggestiveline = FALSE, genomewideline = lq_95 , highlight = fst_parallel_benlim, main= "Little Quarry_ben-lim" , ylim = c(0,1))

#Load lake size GWAS data
gwas <- read.csv("gwas_lakesize.csv", sep = ";")
gwas$lg_bp_group <- as.character(gwas$lg_bp_group)

gwas_sig_vec <- gwas$lg_bp_group[which(gwas$pvalue < 0.05)]

#Combined dataset for windows present in FST ben-lim and GWAS lake size datasets
fst_gwas <- fst_species_pairs[which(fst_species_pairs$lg_bp_group %in% gwas$lg_bp_group),]
fst_gwas$gwas_pvalue <- gwas$pvalue[match(fst_gwas$lg_bp_group, gwas$lg_bp_group)]

#Multiple testing correction
fst_gwas$gwas_pvalue_bonf <- p.adjust(fst_gwas$gwas_pvalue, method = "bonferroni", n = 4985)
fst_gwas$gwas_pvalue_fdr <- p.adjust(fst_gwas$gwas_pvalue, method = "fdr", n = 4985)
#No windows remained significant after correction, uncorrected p-values were used for all further analyses

fstgwas_95 <- quantile(fst_gwas$mean_norm, na.rm = TRUE, probs = 0.95)

fst_gwas_doublesig <- fst_gwas$lg_bp_group[which(fst_gwas$mean_norm >= fstgwas_95 & fst_gwas$gwas_pvalue <= 0.05)]
fst_gwas_doublesig <- as.character(fst_gwas_doublesig)

#Combined dataset for windows present in FST ben-lim and FST small-large lakes: fst_sym_allo
fst_sym_allo$mean_norm_sym <- fst_species_pairs$mean_norm_sym[match(fst_sym_allo$lg_bp_group, fst_species_pairs$lg_bp_group)]

fst_allo_sub_95 <- quantile(fst_sym_allo$mean_norm_allo, na.rm = TRUE, probs = 0.95) 
fst_sym_sub_95 <- quantile(fst_sym_allo$mean_norm_sym, na.rm = TRUE, probs = 0.95) 

##Enrichment analysis

#Outlier/significant windows for different datasets
fst_species_pairs$sig <- fst_species_pairs$mean_norm >= norm_sym_95
fst_small_large$sig <- fst_small_large$mean_norm >= norm_allo_95
gwas$sig <- gwas$pvalue < 0.05
fst_gwas$sig <- fst_gwas$mean_norm >= fstgwas_95 & fst_gwas$gwas_pvalue <= 0.05
fst_sym_allo$sig <- fst_sym_allo$mean_norm_allo >= fst_allo_sub_95 & fst_sym_allo$mean_norm_sym >= fst_sym_sub_95

#Pick datasets accordingly

#FST benthic-limnetic outlier, exclude FST-GWAS double outlier
data <- fst_species_pairs

#FST small-large lakes
data <- fst_small_large

#Lake size GWAS
data <- gwas

#FST-GWAS Double outlier
data <- fst_gwas

#FST-FST Double outlier
data <- fst_sym_allo

#Enrichment calculation
emp <- data%>%
  group_by(lg)%>%
  summarise(freq=sum(sig == TRUE))%>%
  ungroup()%>%
  data.frame()

temp <- emp

perm <- 10000

for(i in 1:perm){
  a <- sample(data$sig, replace=FALSE)
  b <- data$lg
  c <- as.data.frame(cbind(a,b))
  colnames(c)[1] <-c("pval")
  colnames(c)[2] <-c("lg")
  d <- c%>%
    group_by(lg)%>%
    summarise(freq=sum(pval))%>%
    ungroup()%>%
    data.frame()
  temp[,i+2] <- d[,2]
}

temp[1:10,1:10]

pval_fun<-function(perm,h){(length(perm[perm>=h])+1)/(10001)} 
pval<-data.frame()
for (r in 1:nrow(temp)) {
  perm<-temp[r,3:10002]
  h<-temp[r,2] #empirical value 
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  

colnames(pval)<-c("pval")
perm_repeat <-data.frame(emp,pval)
head(perm_repeat)

#P value detection .05 threshold 
perm_repeat$parallel_0.05<-perm_repeat$pval<=0.05 
table(perm_repeat$parallel_0.05)

#final output
perm_repeat_final <- perm_repeat

final <- subset(perm_repeat_final, perm_repeat_final$parallel_0.05 == TRUE)

final

#Enrichment plots, make sure which data is used above
chr_vec <- NA
for (i in 1:21){
  chr <- data[which(data$lg == i),]
  sig <- length(chr$lg[which(chr$sig == TRUE)])
  prop <- sig/length(chr$lg)
  chr_vec <- c(chr_vec, prop)
}

chr_vec <- chr_vec[-1]

chr_perm_df <- data.frame(data.frame(matrix(NA, nrow = 21),
                                     chrom=c(1:21),
                                     prop=as.numeric(0),
                                     stringsAsFactors=FALSE))

chr_perm_df$prop <- chr_vec

plot(chr_perm_df$prop ~ chr_perm_df$chrom, type = "lines", xlab = "Chromosome", ylab = "Proportion")
mapply(axis, side = 1, at = c(1:21), labels = c(1:21), cex.axis = 1)

#Visualization of different categories (non-significant, FST outlier, GWAS significant, double outlier)
#FST species pairs & lake size GWAS

#Define which windows are significant for either FST, GWAS or both
fst_gwas$non_sig <- fst_gwas$mean_norm_sym < fstgwas_95 & fst_gwas$gwas_pvalue > 0.05
fst_gwas$gwas_sig <- fst_gwas$mean_norm_sym < fstgwas_95 & fst_gwas$gwas_pvalue <= 0.05
fst_gwas$fst_sig <- fst_gwas$mean_norm_sym >= fstgwas_95 & fst_gwas$gwas_pvalue > 0.05
fst_gwas$double_sig <- fst_gwas$mean_norm_sym >= fstgwas_95 & fst_gwas$gwas_pvalue <= 0.05

fst_gwas$gwas_log <- -log10(fst_gwas$gwas_pvalue)

#Create separate dataframes for each of the four categories
non_sig <- fst_gwas[which(fst_gwas$mean_norm_sym < fstgwas_95 & fst_gwas$gwas_pvalue > 0.05),]
gwas_sig <- fst_gwas[which(fst_gwas$mean_norm_sym < fstgwas_95 & fst_gwas$gwas_pvalue <= 0.05),]
fst_sig <- fst_gwas[which(fst_gwas$mean_norm_sym >= fstgwas_95 & fst_gwas$gwas_pvalue > 0.05),]
double_sig <- fst_gwas[which(fst_gwas$mean_norm_sym >= fstgwas_95 & fst_gwas$gwas_pvalue <= 0.05),]

#All repeated windows in fst_gwas dataframe
fst_parallel_benlim2 <- fst_gwas[which(fst_gwas$lg_bp_group %in% fst_parallel_benlim),]

#Check for overlap between FSTs in 95th percentile in all three lakes with categories above
parallel_benlim_fst <- fst_parallel_benlim2[which(fst_parallel_benlim2$lg_bp_group %in% fst_sig$lg_bp_group),]
parallel_benlim_gwas <- fst_parallel_benlim2[which(fst_parallel_benlim2$lg_bp_group %in% gwas_sig$lg_bp_group),]
parallel_benlim_double <- fst_parallel_benlim2[which(fst_parallel_benlim2$lg_bp_group %in% double_sig$lg_bp_group),]

plot(non_sig$gwas_log ~ non_sig$mean_norm, ylim = c(0,4.2), xlim = c(-1,3.5), col = "grey", main = "Shared genomic loci", ylab = "-log10(p)", xlab = "normalized FST")
points(gwas_sig$gwas_log ~ gwas_sig$mean_norm, ylim = c(0,4.2), xlim = c(-1,3.5), col = "blue", pch = 16)
points(fst_sig$gwas_log ~ fst_sig$mean_norm, ylim = c(0,4.2), xlim = c(-1,3.5), col = "red", pch = 16)
points(double_sig$gwas_log ~ double_sig$mean_norm, ylim = c(0,4.2), xlim = c(-1,3.5), col = "magenta2", pch = 17)

abline(h = 1.30103, lty = 2)
abline(v = 1.93, lty = 2)

#Chi square test
length(non_sig$lg_bp_group)
length(gwas_sig$lg_bp_group)
length(fst_sig$lg_bp_group)
length(double_sig$lg_bp_group)

dat <- matrix(c(55,195,960,3775), ncol=2)

chisq.test(dat)$expected

chisq.test(dat)

#Test for correlation between GWAS p value & FST for all windows
plot(fst_gwas$gwas_log ~ fst_gwas$mean_norm)
abline(lm(fst_gwas$gwas_log ~ fst_gwas$mean_norm))
cor.test(fst_gwas$gwas_log, fst_gwas$mean_norm, method = "spearman")

#Test whether correlation between GWAS p value & FST is different among different categories
plot(non_sig$gwas_log ~ non_sig$mean_norm)
abline(lm(non_sig$gwas_log ~ non_sig$mean_norm))
cor.test(non_sig$gwas_log, non_sig$mean_norm, method = "spearman")

plot(fst_sig$gwas_log ~ fst_sig$mean_norm)
abline(lm(fst_sig$gwas_log ~ fst_sig$mean_norm))
cor.test(fst_sig$gwas_log, fst_sig$mean_norm, method = "spearman")

plot(gwas_sig$gwas_log ~ gwas_sig$mean_norm)
abline(lm(gwas_sig$gwas_log ~ gwas_sig$mean_norm))
cor.test(gwas_sig$gwas_log, gwas_sig$mean_norm, method = "spearman")

plot(double_sig$gwas_log ~ double_sig$mean_norm)
abline(lm(double_sig$gwas_log ~ double_sig$mean_norm))
cor.test(double_sig$gwas_log, double_sig$mean_norm, method = "spearman")

#Plot Spearman's rho across categories
group <- as.factor(c("a_non_sig", "b_fst_sig", "c_gwas_sig", "d_double_sig"))
rho <- c(-0.07919248, -0.046361, 0.05961251, 0.06998557)
plot(rho ~ group)
points(rho ~ group, type = "points")

#FST species pairs & FST small-large lakes

#Define which windows are significant for either FST, GWAS or both
fst_sym_allo$non_sig <- fst_sym_allo$mean_norm_sym < fst_sym_sub_95 & fst_sym_allo$mean_norm_allo < fst_allo_sub_95
fst_sym_allo$sym_sig <- fst_sym_allo$mean_norm_sym >= fst_sym_sub_95 & fst_sym_allo$mean_norm_allo < fst_allo_sub_95
fst_sym_allo$allo_sig <- fst_sym_allo$mean_norm_sym < fst_sym_sub_95 & fst_sym_allo$mean_norm_allo >= fst_allo_sub_95
fst_sym_allo$double_fst_sig <- fst_sym_allo$mean_norm_sym >= fst_sym_sub_95 & fst_sym_allo$mean_norm_allo >= fst_allo_sub_95

#Create separate dataframes for each of the four categories
non_sig <- fst_sym_allo[which(fst_sym_allo$mean_norm_sym < fst_sym_sub_95 & fst_sym_allo$mean_norm_allo < fst_allo_sub_95),]
sym_sig <- fst_sym_allo[which(fst_sym_allo$mean_norm_sym >= fst_sym_sub_95 & fst_sym_allo$mean_norm_allo < fst_allo_sub_95),]
allo_sig <- fst_sym_allo[which(fst_sym_allo$mean_norm_sym < fst_sym_sub_95 & fst_sym_allo$mean_norm_allo >= fst_allo_sub_95),]
double_fst_sig <- fst_sym_allo[which(fst_sym_allo$mean_norm_sym >= fst_sym_sub_95 & fst_sym_allo$mean_norm_allo >= fst_allo_sub_95),]

#All repeated windows
fst_parallel_benlim3 <- fst_sym_allo[which(fst_sym_allo$lg_bp_group %in% fst_parallel_benlim),]

#Check for overlap between FSTs in 95th percentile in all three lakes with categories above
parallel_benlim_sym <- fst_parallel_benlim3[which(fst_parallel_benlim3$lg_bp_group %in% sym_sig$lg_bp_group),]
parallel_benlim_allo <- fst_parallel_benlim3[which(fst_parallel_benlim3$lg_bp_group %in% allo_sig$lg_bp_group),]
parallel_benlim_double <- fst_parallel_benlim3[which(fst_parallel_benlim3$lg_bp_group %in% double_fst_sig$lg_bp_group),]

plot(non_sig$mean_norm_allo ~ non_sig$mean_norm_sym, ylim = c(-0.87,-0.45), xlim = c(-1,3.5), col = "grey", main = "Shared genomic loci", ylab = "normalized FST allo", xlab = "normalized FST sym")
points(sym_sig$mean_norm_allo ~ sym_sig$mean_norm_sym, ylim = c(-0.87,-0.45), xlim = c(-1,3.5), col = "red", pch = 16)
points(allo_sig$mean_norm_allo ~ allo_sig$mean_norm_sym, ylim = c(-0.87,-0.45), xlim = c(-1,3.5), col = "blue", pch = 16)
points(double_fst_sig$mean_norm_allo ~ double_fst_sig$mean_norm_sym, ylim = c(-0.87,-0.45), xlim = c(-1,3.5), col = "magenta2", pch = 17)

abline(h = -0.704, lty = 2)
abline(v = 1.91, lty = 2)

#Chi square test
length(non_sig$lg_bp_group)
length(sym_sig$lg_bp_group)
length(allo_sig$lg_bp_group)
length(double_fst_sig$lg_bp_group)

dat <- matrix(c(15,202,202,3919), ncol=2)

chisq.test(dat)$expected

chisq.test(dat)

#Test for overlap of windows from different categories (fst-sig, gwas-sig, double-sig) with QTL data

#Load published QTL data
qtl <- read.csv("qtl_ben_lim.csv", sep = ";")

# set window size to use for QTL data
win=50000
#assign window ids 
qtl$bp_group <- as.integer(qtl$position_to_use/win)+1
qtl$lg_bp_group <- sprintf("%s_%s", qtl$lg, qtl$bp_group)

qtl_win <- as.character(qtl$lg_bp_group)

#Identify overlap between significant windows & QTL data

#FST species pairs
fst_species_pairs_qtl_sig <- qtl[qtl$lg_bp_group %in% ben_lim_outliers,]

#FST small-large lakes
fst_small_large_qtl_sig <- qtl[qtl$lg_bp_group %in% small_large_outliers,]

#Lake size GWAS
gwas_qtl_sig <- qtl[qtl$lg_bp_group %in% gwas_sig_vec,]

#FST-GWAS double outlier
double_qtl_sig <- qtl[qtl$lg_bp_group %in% double_sig$lg_bp_group,]

#FST-FST double outlier
double_fst_qtl_sig <- qtl[qtl$lg_bp_group %in% fst_outlier_combined,]

#Manhattan Plots
#Plot p values for GWAS and FST values

#GWAS
manhattan(gwas, chr="lg", bp="bp_group", snp="lg_bp_group", p="pvalue", logp = TRUE, suggestiveline = FALSE, genomewideline = 1.30103, main="GWAS windows", ylim = c(0,5))

#FST species pairs
manhattan(fst_species_pairs, chr="lg", bp="bp_group", snp="lg_bp_group", p="mean_norm_sym", logp = FALSE, suggestiveline = FALSE, genomewideline = norm_sym_95, highlight = fst_parallel_benlim, ylim = c(-1,4), main="FST windows", ylab = "normalized FST-species pairs")

#FST small-large lakes
manhattan(fst_small_large, chr="lg", bp="bp_group", snp="lg_bp_group", p="mean_norm_allo", logp = FALSE, suggestiveline = FALSE, genomewideline = norm_allo_95, main="FST windows", ylab = "normalized FST-small&large lakes", highlight = fst_parallel_small_large, ylim = c(-0.9,-0.4))

##Test whether overlap between repeatedly differentiated windows and FST outliers is significant

#Permutation test

#Species pairs
data <- fst_species_pairs

perm <- 10000

vec <- NA

for(i in 1:perm){
  a <- sample(data$lg_bp_group, 287, replace=FALSE)
  b <- sum(a %in% fst_parallel_benlim)
  vec <- c(vec, b)
}

vec <- vec[-1]

vec2 <- outer(67, vec, ">")
vec2 <- as.character(vec2)
pvalue <- 1-sum(vec2 == "TRUE")/10000

#Small-large lakes
data <- fst_small_large

perm <- 10000

vec <- NA

for(i in 1:perm){
  a <- sample(data$lg_bp_group, 292, replace=FALSE)
  b <- sum(a %in% fst_parallel_small_large)
  vec <- c(vec, b)
}

vec <- vec[-1]

vec2 <- outer(19, vec, ">")
vec2 <- as.character(vec2)
pvalue <- 1-sum(vec2 == "TRUE")/10000