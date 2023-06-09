#load the required library 
library(ggplot2)

#read in data files for heterozygosity rate and deleterious mutation proportions per population
#heterozygosity rates and deleterious mutation proportions were obtained using the "per_sample_heterozygosity.sh" and "estimate_deleterious_snp_proportion.sh" scripts; additional information for each population was added for plotting ("Category" and "Range")
#for the heterozygosity rate input file, DNA and temporal replicates were excluded; heterozygosity was averaged per population
het<-read.csv(file="heterozygosity_rate.csv",header = TRUE, sep = ',')
del<-read.csv(file="percent_HIGH_SNPs.csv",header = TRUE, sep = ',')

#populations grouped by range and ordered within each range from low to high genome-wide heterozygosity
het$Population = factor(het$Population, levels=c("ALA","ORA","SAR","DUV","POL","STJ","BRE","MARi","VOL","FLA","LEV","LAK","LOW","LEE","HIG","MAN","TAM","CMB","HAM","COL","CIT","TIF","GLA","HEN","MIA","BRO","STL","MON","PAL","SOR","STP","LHA","LHB","MAR","GUA","JIC","CAB","POR","ESM","BAH"))
het$Category = factor(het$Category, levels=c("Invasive","SOR","MAR","LHA","LHB","GUA","JIC","CAB","BAH","ESM","POR"))
het$Range = factor(het$Range, levels=c("Invasive","Native"))

del$Population = factor(del$Population, levels=c("ALA","ORA","SAR","DUV","POL","STJ","BRE","MARi","VOL","FLA","LEV","LAK","LOW","LEE","HIG","MAN","TAM","CMB","HAM","COL","CIT","TIF","GLA","HEN","MIA","BRO","STL","MON","PAL","SOR","STP","LHA","LHB","MAR","GUA","JIC","CAB","POR","ESM","BAH"))
del$Category = factor(del$Category, levels=c("Invasive","SOR","MAR","LHA","LHB","GUA","JIC","CAB","BAH","ESM","POR"))
del$Range = factor(del$Range, levels=c("Invasive","Native"))

###heterozygosity per population, ordered from low to high (Fig. S8)
ggplot(het, aes(Population,Heterozygosity, colour=Category, fill=Category))+
  geom_jitter(width = 0.15)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  xlab("\\nPopulation")+
  ylab("Heterozygosity\\n")+
  ylim(0,0.12)+
  coord_flip()+
  scale_colour_manual(values = c("grey70","#e09148ff","#e09148ff","#e09148ff","#e09148ff","#e09148ff","#b53e35ff","#9dc8e9ff","#3f919bff","#415b9eff","#853789ff"))

###deleterious SNP proportion per population, ordered from low to high (Fig. S8)
ggplot(del, aes(Population,Percent_high_effect_SNPs, group=Category, fill=Category))+
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  xlab("\\nPopulation")+
  ylab("% Deleterious SNPs\\n")+
  ylim(0,0.05)+
  coord_flip()+
  scale_fill_manual(values = c("grey70","#e09148ff","#e09148ff","#e09148ff","#e09148ff","#e09148ff","#b53e35ff","#9dc8e9ff","#3f919bff","#415b9eff","#853789ff"))
  
#heterozygosity per range histogram (Fig. S8)
ggplot(het, aes(x=Mean_heterozygosity))+
  geom_histogram(aes(fill=Range), binwidth = 0.0025, position="identity")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  xlab("\\nHeterozygosity")+
  xlim(0,0.1)+
  scale_fill_manual(values = c("grey60","black"))

#deleterious SNP proportion per range histogram (Fig. S8)
ggplot(del, aes(x=Percent_high_effect_SNPs))+
  geom_histogram(aes(fill=Range), binwidth = 0.0015, position="identity")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  xlab("\\n% Deleterious SNPs\\n")+
  xlim(0,0.05)+
  scale_fill_manual(values = c("grey60","black"))

#heterozygosity per range (Fig. 1)
ggplot(het, aes(x=Range, y=Mean_heterozygosity))+
  geom_jitter(aes(colour=Category), size=6, width=0.2, alpha=0.8)+
  theme_bw()+
  theme(axis.title=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  ylim(0,0.1)+
  scale_colour_manual(values = c("grey70","#e09148ff","#e09148ff","#e09148ff","#e09148ff","#e09148ff","#b53e35ff","#9dc8e9ff","#3f919bff","#415b9eff","#853789ff"))

#deleterious SNP proportion per range (Fig. 1)
ggplot(del, aes(x=Range, y=Percent_high_effect_SNPs))+
  geom_jitter(aes(colour=Category), size=6, width=0.2, alpha=0.8)+
  theme_bw()+
  theme(axis.title=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  ylim(0,0.05)+
  scale_colour_manual(values = c("grey70","#e09148ff","#e09148ff","#e09148ff","#e09148ff","#e09148ff","#b53e35ff","#9dc8e9ff","#3f919bff","#415b9eff","#853789ff"))

#test for differences between ranges in the frequency of deleterious mutations using a Wilcoxon test
del_native <- del[del$Range == 'Native',]
del_native_proportion<-del_native$Percent_high_effect_SNPs

del_invasive <- del[del$Range == 'Invasive',]
del_invasive_proportion<-del_invasive$Percent_high_effect_SNPs

wilcox.test(del_native_proportion, del_invasive_proportion)


#test for differences between ranges using a Wilcoxon test
het_native <- na.omit(het[het$Range == 'Native',])
mean_heterozygosity_native<-het_native$Mean_heterozygosity

het_invasive <- na.omit(het[het$Range == 'Invasive',])
mean_heterozygosity_invasive<-het_invasive$Mean_heterozygosity

wilcox.test(mean_heterozygosity_native, mean_heterozygosity_invasive)

