#load relevant libraries 
library(SNPRelate)
library(reshape2)
library(ggplot2)

#read in the filtered VCF file
my_vcf <-"FL_Cuba_IRN.with_replicates.70d.maf1.vcf"
snpgdsVCF2GDS(my_vcf, "all_with_replicates.gds", method="biallelic.only")
genofile <- snpgdsOpen('all_with_replicates.gds')

#calculate IBS
ibs <- snpgdsIBS(genofile, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, num.thread = 1, verbose = TRUE)

#output sample names (used to add to the IBS matrix) and the IBS matrix (IBS.txt)
write.table(ibs$sample.id,file="sample_names.txt", col.names=FALSE)
write.table(ibs$ibs,file="IBS.txt",col.names = FALSE)

#get single-column IBS values; the output ("IBS_one_column.temp") is used as input for the trim_IBS_one_column.sh script
dist<-ibs$ibs
distances<-as.dist(dist)
df <- melt(as.matrix(distances), varnames = c("row", "col"))
write.table(df$value, file="IBS_one_column.temp")

#read in the single-column IBS data (output from trim_IBS_one_column.sh)
IBS_one_column <- read.csv(file="IBS_one_column.txt",header = FALSE, sep = ,)

#plot IBS histogram
ggplot(data = IBS_one_column, aes(x=V1)) +
  geom_histogram(binwidth=.0002, colour="black", fill="black")+
  xlab("\\nPairwise IBS values")+
  #0.99296152 is the lowest IBS value observed among DNA replicates (this is used 
  #as a threshold to identify erroneously duplicated samples)
  geom_vline(xintercept = 0.99296152)+
  ylab("Count\\n")+
  xlim(0.8,1)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
