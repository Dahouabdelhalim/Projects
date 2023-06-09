#load the SNPRelate library 
library(SNPRelate)

#calculate identity-by-state separately for each population and also for each of the native and invasive ranges
population_list<-c("ALA","BAH","BRE","BRO","CAB","CCH","CIT","CMB","COL","DUV","ESM","FLA","GLA","GRC","GUA","HAM","HEN","HIG","HOU","JIC","LAK","LEE","LEV","LHA","LHB","LOW","MAN","MAR","MARi","MIA","MON","NOR","ORA","PAL","POL","POR","SAR","SOR","STJ","STL","STP","TAM","TIF","VOL","native","invasive")
for (population in population_list) {
  infile<-paste (population,".vcf",sep="")
  outfile<-paste (population,"_one_col.temp",sep="")
  vcf<-infile
  snpgdsVCF2GDS(vcf, population, method="biallelic.only")
  genofile <- snpgdsOpen(population)
  #ibs
  ibs <- snpgdsIBS(genofile, autosome.only = FALSE, remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, num.thread = 1, verbose = TRUE)
  dist<-ibs$ibs
  distances<-as.dist(dist)
  library(reshape2)
  df <- melt(as.matrix(distances), varnames = c("row", "col"))
  write.table(df$value, file=outfile)
}