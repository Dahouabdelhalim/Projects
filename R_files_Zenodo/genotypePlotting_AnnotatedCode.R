### VCF was produced using pipeline for output from GATK
### i then applied some additional filtering steps to the vcf:

# only keep sites where the philodice grandmothers are fixed for the non-reference allele
awk '($10 ~ /^1/)' ColiasAll.beagle.zchr.vcf > test1.vcf
awk '($11 ~ /^1/)' test1.vcf > test2.vcf

# only keep sites where the hybrid uncle is a heterozygote
awk '($13 ~ /^0\\/1/ || $13 ~ /^0\\|1/ )' test2.vcf > test1.vcf

# only keep sites where the eurytheme grandma is homozygous reference
awk '($14 ~ /^0\\/0/ || $14 ~ /^0\\|0/ )' test1.vcf > test2.vcf

# only keep sites where the eurytheme grandma is homozygous reference
awk '($14 ~ /^0\\/0/ || $14 ~ /^0\\|0/ )' test1.vcf > test2.vcf

# after filtering, re-add the vcf header
head -154 ColiasAll.beagle.zchr.vcf > header.txt
cat header.txt test2.vcf > test3.vcf

# drop individuals 00-11-041 and 00-79-062
vcftools --remove-indv 00-11-041 --remove-indv 00-79-062 --vcf test3.vcf --recode --out ColiasAll.beagle.zchr.filtered.vcf

# and then bgzip and tabix the file to prepare it for vidualisation by genotype_plot
bgzip -c ColiasAll.beagle.zchr.filtered.vcf > ColiasAll.beagle.zchr.filtered.vcf.gz
tabix -f -p vcf ColiasAll.beagle.zchr.filtered.vcf.gz 

### read genotype_plot.R, from https://github.com/JimWhiting91
source("./genotype_plot/R/genotype_plot.R")

# read a file with a list of individuals in one column and a list of populations in the ohter (in this case column 1 == column 2)
popmap <- read.table('popmap_ColiasAll.txt',h=T)

new_plot <- genotype_plot(vcf="ColiasAll.zchrFUSED_filtered.vcf.gz",
                          chr="ChrZ",
                          start=1,
                          end=15299928,
                          popmap=popmap,
                          cluster=FALSE,
                          snp_label_size=100000,
                          colour_scheme=c('#b10079', '#ff9d2b', '#00ade8',"grey","grey","grey","grey","grey","grey"))
new_plot

### make version for main figure
popmap2 <- head(popmap,24)
popmap2 <- popmap2[popmap2$ind != '00-79-055',]
## filtering 
# only keep sites where the philodice grandmas are homozygous alt
awk '($10 ~ /^1\\/1/ || $10 ~ /^1\\|1/ )' ColiasAll.zchrFUSED_massonly.vcf.gz.recode.vcf > test1.vcf
awk '($11 ~ /^1\\/1/ || $11 ~ /^1\\|1/ )' test1.vcf > test2.vcf
# only keep sites where the 00-79-055 is hom ref
awk '($27 ~ /^0\\/0/ || $27 ~ /^0\\|0/ )' test2.vcf > test1.vcf

# sort test1
sort -k2 -n test1.vcf > test2.vcf

# after filtering, re-add the vcf header
head -11 ColiasAll.zchrFUSED_massonly.vcf.gz.recode.vcf > header.txt
cat header.txt test2.vcf > test3.vcf

vcftools --remove-indv 00-79-055 --vcf test3.vcf --recode --out test4.vcf

bgzip -c test4.vcf.recode.vcf > ColiasAll.zchrFUSED_massonly.vcf.gz
tabix -f -p vcf ColiasAll.zchrFUSED_massonly.vcf.gz 

new_plot <- genotype_plot(vcf="ColiasAll.zchrFUSED_massonly.vcf.gz",
                          chr="ChrZ",
                          start=1,
                          end=15299928,
                          popmap=popmap2,
                          cluster=FALSE,
                          snp_label_size=100000,
                          colour_scheme=c('#b10079', '#ff9d2b', '#00ade8',"grey","grey","grey","grey","grey","grey"))
new_plot

