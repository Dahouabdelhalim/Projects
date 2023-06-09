#R-script for estimating the inbreeding coefficient FROH from Runs of Homozygosity (ROHs) produced by PLINK using the script "Estimating_inbreeding_coefficients.txt"
#Alina Niskanen, original script October 2017, modified May 2020
#alina.niskanen@gmail.com

#Bring ROHs generated before in PLINK to find inbreeding that happened at most approximately 10 generations ago
ROH_10_generations_ind <- read.table("Helgeland_roh10gen.hom", header = T, stringsAsFactors = F)
ROH_10_total <- read.table("Helgeland_roh10gen.hom.indiv", header = T, stringsAsFactors = F)

#The genotyped length of genome in bp for LD-pruned SNPs
SNP_coverage <- 914814900

#Estimate the length and number of ROHs for the whole data
min(ROH_10_generations_ind$KB)
max(ROH_10_generations_ind$KB)
mean(ROH_10_generations_ind$KB)
max(sort(table(ROH_10_generations_ind$KB), decreasing = T))

#Histogram for all ROH lengths
hist(ROH_10_generations_ind$KB, breaks = 200, xlim=c(0,80000), xlab = "ROH kb, 10 generations estimate", ylab = "Frequency", main = "Histogram of individual ROH lengths for Helgeland sparrows")

#Inbreeding estimate FROH (ROH is in KB and SNP-coverage in basepairs, thus divide by 1000)
ROH_10_total$FROH <- ROH_10_total$KB/(SNP_coverage/1000)

#Output a table with FROH, length of total genome in ROHs and the number of ROH segments
write.table(ROH_10_total[,c(2,4:7)], file="FROH_10_generations.txt", col.names = T, row.names = F, quote = F)


