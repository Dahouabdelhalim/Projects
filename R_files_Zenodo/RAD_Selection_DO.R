# Run the FUNCTIONS, LOAD and CLEAN R scripts prior to producing the outputs
source("RAD_Selection_FUNCTIONS.R")
source("RAD_Selection_LOAD.R")
source("RAD_Selection_CLEAN.R")

# Fst ---------------------------------------------------------------------

# Genome wide Fst ---------------------------------------------------------

### Plot Locus Fst across the genome * FIGURE S10 * ###

manhattan(na.omit(Locus_Fst_MAD_mapped), chr = "CHR", bp = "POS", p = "FST", snp = "SNP", logp = FALSE, ylab = expression("F"[ST]),
          chrlabs = c(1:15,"","","","","","","","","","", "", "","","","", "", ""))
unique(Locus_Fst_MAD_mapped$CHR)

manhattan(na.omit(Locus_Fst_CAN_mapped), chr = "CHR", bp = "POS", p = "FST", snp = "SNP", logp = FALSE, ylab = expression("F"[ST]),
          chrlabs = c(1:15,"","","","","","","","","","", "", "","","","", "", ""))
unique(Locus_Fst_CAN_mapped$CHR)



# Locus Fst ---------------------------------------------------------------

# plot without SNPs fixed in one arch
ggplot(Fst, aes(x= Fst_Maderia, y = Fst_Canaries))+
  ylab(expression(paste("Canary Islands", " F"[ST])))+
  xlab(expression(paste("Madeiran archipelago", " F"[ST])))+
  xlim(0, 0.6)+
  ylim(0, 0.25)+
  geom_point(size = 1.2)+
  geom_point(col = ifelse(Fst$Fst_Canaries > 0.23 | Fst$Fst_Maderia > 0.47, "red", "black"))+
  geom_text_repel(aes(label=ifelse(Fst$Fst_Maderia>0.47 | Fst$Fst_Canaries > 0.23, as.character(Fst$SNP_CAN),'')), nudge_x =-0.00, nudge_y=0.03, segment.size = 0.3, size = 3.5)+
  theme_bw()

# plot with all SNPs, including those with MAF = 0 in one of the archs. This is the plot presented in the paper.Figure 5.
# * FIGURE 5 *
ggplot(All_SNPs, aes(x= Fst_Maderia, y = Fst_Canaries))+
  ylab(expression(paste("Canary Islands", " F"[ST])))+
  xlab(expression(paste("Madeiran archipelago", " F"[ST])))+
  xlim(0, 0.6)+
  ylim(0, 0.27)+
  geom_point(size = 1.2)+
  geom_text_repel(aes(label=ifelse(All_SNPs$Fst_Maderia > 0.47 | All_SNPs$Fst_Canaries > 0.23, as.character(All_SNPs$SNP_CAN),'')), nudge_x =-0.00, nudge_y=0.03, segment.size = 0.3, size = 3.5)+
  geom_point(col = ifelse(All_SNPs$Fst_Canaries > 0.23 | All_SNPs$Fst_Maderia > 0.47, "red", "black"))+
  #geom_point(col = ifelse(All_SNPs$Fst_Canaries <0.0001 | All_SNPs$Fst_Maderia <0.0001, "grey", "black"))+ # just checking to see which SNPs are fixed in one of the archs 
  theme_claudia()



cor.test(Locus_Fst_CAN_reduced$FST,Locus_Fst_MAD_reduced$FST, use = "complete.obs")
cor.test(Fst$Fst_Canaries, Fst$Fst_Maderia, use = "complete.obs") # when the lower Fst bound is set to zero
cor.test(All_SNPs$Fst_Maderia, All_SNPs$Fst_Canaries, use = "complete.obs") # just check that this is the still not significant when we include all of the SNPs which are fixed in one arch

## spearmans rank
cor.test(Fst$Fst_Canaries, Fst$Fst_Maderia, method = "spearman") # lower bound zero 
cor.test(All_SNPs$Fst_Canaries, All_SNPs$Fst_Maderia, method = "spearman")


# Extracting values for Fst table -----------------------------------------

subset(Fst, SNP_CAN == "2549s24")
subset(Locus_Fst_CAN, SNP == "1458s107")
subset(Locus_Fst_MAD, SNP == "855s106")

subset(Locus_Fst_CAN, FST > 0.2) # identify which loci have high Fst in Canaries 

#       CHR   SNP     POS    NMISS   FST
# 2037   5  219s24 58990367   174 0.234473
# 3264  17 2294s59  3527112   175 0.217134
# 4003  28  933s86  2504629   172 0.204064

# Get a value for the Fst of these SNPs then in the Maderian arch 
subset(Locus_Fst_MAD, SNP == "219s24") # missing, is this because it is fixed for one of the alleles in Maderia?
# check the allele frequencies 
subset(Locus_Fst_MAD, SNP == "2294s59")
subset(Locus_Fst_MAD, SNP == "933s86") # missing 

# do the same for Maderia 

subset(Locus_Fst_MAD, FST > 0.4)

#       CHR    SNP     POS  NMISS  FST
# 1446   6  1585s94 33504910    49 0.561289
# 1447   6 1585s112 33504928    49 0.561289
# 1604   9 1458s107  5312107    46 0.422414
# 2476  24   790s54   718487    48 0.473106
# 2731  29  351s116 65539694    49 0.448589

subset(Locus_Fst_CAN, SNP == "1585s94")
subset(Locus_Fst_CAN, SNP == "1585s112")
subset(Locus_Fst_CAN, SNP == "1458s107")
subset(Locus_Fst_CAN, SNP == "790s54")
subset(Locus_Fst_CAN, SNP == "351s116")

# just check these values with those that have been plotted

subset(Fst, SNP_CAN == "219s24") # not plotted as expected
subset(Fst, SNP_CAN == "2294s59") # correct
subset(Fst, SNP_CAN == "933s86") # Not plotted
subset(Fst, SNP_CAN == "1585s94") # correct
subset(Fst, SNP_CAN == "1585s112") # correct 
subset(Fst, SNP_CAN == "1458s107") # not plotted
subset(Fst, SNP_CAN == "790s54") # correct
subset(Fst, SNP_CAN == "351s116") # correct 


# Eigenvecotor PCA --------------------------------------------------------

F1 <- ggplot(egvec_CAN,aes(x = PC1,y = PC2,col = Population))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.spacing = unit(0.8, "lines"),
        axis.title.x = element_text(size = 10, colour = "black", vjust = 1),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.line = element_blank(),
        legend.title = element_blank())+
  scale_shape_manual(values=c(2,17,6,19,1,10,15,0,9))+
  scale_color_manual(values=c(dark_blue, dark_blue,dark_blue, fire_orange,fire_orange,fire_orange,purple,purple,purple))+
  geom_point(aes(shape=Population), size = 1.5)

F1


F2 <- ggplot(egvec_MAD,aes(x = PC1,y = PC2,col = Population))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.spacing = unit(0.8, "lines"),
        axis.title.x = element_text(size = 10, colour = "black", vjust = 1),
        axis.title.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.line = element_blank(),
        legend.title = element_blank())+
  scale_shape_manual(values=c(2,17,6))+
  scale_color_manual(values=c(dark_blue, fire_orange,purple))+
  geom_point(aes(shape=Population), size = 1.5)

F2

# Eigenvalues Manhatton Plot * FIGURE 4 A&B, S8 and S9 *----------------------------------------------

# Canaries 
# Bonferroni line calculation 
# wc -l PIP_BER_parseV-90-3flag.Allsnp.pop.nl.trimmed_Canaries_subset_maf_trim_0.03_rel_trim_Plink_0.2.bim 
# p/number of SNPs. 0.05/4470= 0.00001119

pdf("eg1_CAN.pdf", width = 9, height = 5)
manhattan(eg1_CAN, chr = "CHR",p = "PGC",bp = "BP",
          suggestiveline = F, genomewideline = F,
          chrlabs = c("Un", 1:13,"", "", "","","","","","","", "", "","","","", "1A", "", "4A", "LGE22", "Z", ""), cex.axis = 1, cex = 1)
text(x=570000000, y=8.5, "219s24", cex = 1.2)
suggestiveline <- -log10(1.119e-5)
abline(h= suggestiveline, col = "red")
dev.off()

unique(eg1_CAN$CHR)

qq(eg1_CAN$PGC)# the qq function has been hidden by another package so run "qqman::qq(eg1$PGC)" so it knows where to source qq from 


## now for the Eigenvector 2 plots ###

manhattan(eg2_CAN, chr = "CHR",p = "PGC",bp = "BP",
          suggestiveline = F,
          chrlabs = c("Un", 1:13,"", "", "","","","","","","", "", "","","","", "", "", "", "", "", ""), cex.axis = 0.9, cex = 0.7)
suggestiveline <- -log10(1.119e-5)
abline(h= suggestiveline, col = "red")

unique(eg1_MAD$CHR)

qq(eg2_CAN$PGC)



# Maderia 
# Bonferroni line calculation 
# wc -l PIP_BER_parseV-90-3flag.Allsnp.pop.nl.trimmed_Maderia_subset_maf_trim_0.03_rel_trim_Plink_0.2.bim 
# p/number of SNPs. 0.05/2938 = 0.00001702

pdf("eg1_MAD.pdf", width = 9, height = 5)
manhattan(eg1_MAD, chr = "CHR",p = "PGC",bp = "BP",
          suggestiveline = F, genomewideline = F,
          chrlabs = c("Un", 1:13,"", "", "","","","","","","", "", "","","","", "1A", "", "", "", "Z"), cex.axis = 1, cex = 1)
text(x=420000000, y=5.55, "  1585s94 & 
     1585s112", cex = 1.2)
text(x=925000000, y=5.2, "790s54", cex = 1.2)
suggestiveline <- -log10(1.702e-5)
abline(h= suggestiveline, col = "red")
segments(x0 = 825000000, y0 = 5, x1 = 870000000, y1 = 5.2, lwd=1) 
segments(x0 = 540000000, y0 = 5.8, x1 = 500000000, y1 = 5.8, lwd=1) 
dev.off()

qq(eg1_MAD$PGC)


# now for EV2 

manhattan(eg2_MAD, chr = "CHR",p = "PGC",bp = "BP",
          suggestiveline = F, genomewideline = F,
          chrlabs = c("Un", 1:13, "", "","","","","","","", "", "","","","", "", "", "", "", "", ""), cex.axis = 0.9, cex = 0.7)
suggestiveline <- -log10(1.702e-5)
abline(h= suggestiveline, col = "red")

qq(eg2_MAD$PGC)
