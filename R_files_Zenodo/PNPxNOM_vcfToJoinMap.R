# vcf to JoinMap for PNPxNOM 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(here)
library(stringr)

# read in vcf file 
PNPxNOM_vcf<-read.table("PNPxNOM_7more_f1_r1_f2_f3_homfixF1het.vcf", comment.char="", skip=6929, header=1, check.names=F)

# check column (and with it individual) names
names(PNPxNOM_vcf)

# make a vector containing all the positions of all snps
PNPxNOM_positions<-with(PNPxNOM_vcf,paste(PNPxNOM_vcf$`#CHROM`,POS,sep="_"))

# specify the "core" of the vcf-file (only containing each individual's data at each snp)
PNPxNOM_vcfcore<-PNPxNOM_vcf[,10:length(PNPxNOM_vcf)]
PNPxNOM_vcfstrip<-PNPxNOM_vcfcore

# in the vcfstrip extract only info about snp (1/1, 1/0, 0/1 etc.)
for (i in 1:length(PNPxNOM_vcfcore)){
  PNPxNOM_vcfstrip[,i]<-as.vector(str_extract(PNPxNOM_vcfcore[,i],pattern="./."))
}

# specify parental genotypes
PNPxNOM_grandmother<-PNPxNOM_vcfstrip$`63120.F0fP1`
PNPxNOM_grandfather<-PNPxNOM_vcfstrip$`42957.F0mP1`

# define genotypes of F2s (u=undetermined, h=heterozygous, a=like grandmother, b=like grandfather)
PNPxNOM_vcfstrip2<-PNPxNOM_vcfstrip
for (i in 1:length(PNPxNOM_vcfstrip2)){
  for(j in 1:nrow(PNPxNOM_vcfstrip2)){
    if (PNPxNOM_vcfstrip2[j,i]=="0/1") PNPxNOM_vcfstrip2[j,i]<-"h"
    if (PNPxNOM_vcfstrip2[j,i]=="./.") PNPxNOM_vcfstrip2[j,i]<-"u"
    if (PNPxNOM_vcfstrip2[j,i]=="./1") PNPxNOM_vcfstrip2[j,i]<-"u"
    if (PNPxNOM_vcfstrip2[j,i]=="./0") PNPxNOM_vcfstrip2[j,i]<-"u"
    if (PNPxNOM_vcfstrip2[j,i]==PNPxNOM_grandmother[j]) PNPxNOM_vcfstrip2[j,i]<-"a"
    if (PNPxNOM_vcfstrip2[j,i]==PNPxNOM_grandfather[j]) PNPxNOM_vcfstrip2[j,i]<-"b"
  }}

# remove grandparents
PNPxNOM_joinmapcore<-PNPxNOM_vcfstrip2[,!names(PNPxNOM_vcfstrip2) %in% c("63120.F0fP1","42957.F0mP1")]
names(PNPxNOM_joinmapcore)

# make a dataframe out of it and add pos and chr info
PNPxNOM_joinmap<-data.frame(Pos=as.character(PNPxNOM_positions),Class=as.character(rep("(a,h,b)",nrow(PNPxNOM_joinmapcore))),PNPxNOM_joinmapcore,stringsAsFactors=FALSE)

# write output file 
sink("PNPxNOM_7more_JM_input.txt")
cat("name = PNPxNOM","\\n")
cat("popt = F2","\\n")
cat(paste("nloc =",nrow(PNPxNOM_joinmapcore)),"\\n")
cat(paste("nind =",length(PNPxNOM_joinmapcore)),"\\n","\\n")
for(i in 1:nrow(PNPxNOM_joinmap)){
  cat(paste(PNPxNOM_joinmap[i,]),"\\n")}
cat("\\n")
cat("individual names:","\\n","\\n")
for(i in 1:length(PNPxNOM_joinmapcore)){
  cat(paste(names(PNPxNOM_joinmapcore)[i]),"\\n")}
sink()
 
 
 
# Linkage map generated with JoinMap
# individuals excluded from linkage map construction due to too much missing data:
# 194143.F2m2,194148.F2m2,194160.F2f1,194204.F2m2,194180.F2f1

