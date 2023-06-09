# sequence analysis
#EE25 and EE24 amplicon deep seq data

#file names:

#...R1... fwd reads
#...R2... rv reads (see Geneious for overview)

#S1...PCR amplicon 1 ... pooled samples EE24 E (ctr no gal)
#S2...PCR amplicon 2 ... pooled samples EE24.12CIS (IS-)(first 4columns of 96-well plate together)
#S3...PCR amplicon 3 ... pooled samples EE24.12CIS (IS-)(second 4cols)
#S4...PCR amplicon 4 ... pooled samples EE24.12CIS (IS-)(third 4cols)
#S5...PCR amplicon 2 ... pooled samples EE24.12C30 (IS+) (first 4cols)
#S6...PCR amplicon 3 ... pooled samples EE24.12C30 (IS+)(second 4cols)
#S7...PCR amplicon 4 ... pooled samples EE24.12C30 (IS+)(third 4cols)
#S8,9,10: pooled (p0+p-2) EE25B (IS+) [Wells in column 1-4, 5-8,9-12] 	
#S11,12,13: pool (p0+p-2) EE25B (IS-) [Wells in column 1-4,5-8,9-12] 
#S14: pool EE24A30 (IS+) (pooled all populations) 
#S15: pool EE24AIS (IS-) (pool all populations)	


#DO A COMPARISON OF SNPS IN END OF P0 VERSUS (average # of) single SNPS IN GALK read

#package ShortRead 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ShortRead")
library(ShortRead)

#Biostrings: application downsteram of ShortRead
#https://web.stanford.edu/class/bios221/labs/biostrings/lab_1_biostrings.html   #tutorial
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")



setwd("/Users/itomanek/amplicon_p0_sequencing")

### read in files in a loop
allfiles = list.files(pattern = "*.fastq")
for (i in 1:length(allfiles)){ 
  x <- readFastq(allfiles[i]) 
  assign(paste0(allfiles[i]), x)
}#end loop


#### Analyze sequence patterns (mutations spotted previously from Sanger sequencing (Geneious))

### analyze p0-2 gaps in EE25B - IS+ and IS-
#### read in fastq files
#IS+ pops
S8=readFastq("S8_S8_R1_001.fastq.gz")
S9=readFastq("S9_S9_R1_001.fastq.gz")
S10=readFastq("S10_S10_R1_001.fastq.gz")
#IS- pops
S11=readFastq("S11_S11_R1_001.fastq.gz")
S12=readFastq("S12_S12_R1_001.fastq.gz")
S13=readFastq("S13_S13_R1_001.fastq.gz")
####

#p0-2 sequence patterns (39 bases) - Fw (R1) reads
anc=DNAString("GGTGAGATAACTCCGTAGTTGACTACGCGTCCCTCTAGG")
gap=DNAString("GGTGAGATAACTCCGTAGTTGACTACGCGTCCCTTGACG") #12bp deletion ("gap")
anc==gap
compareStrings(anc,gap) #3 SNPs -gap is at the end of the read 

#p0 seq pattern Fw (R1) reads
p0=DNAString("TCAAGATCCTCTTAATAAGCCCCCGTCACTGTTGGTTGT") 
#p0 seq pattern Rv (R2) reads - reverse complemented to match the raw reads
H5_CT=reverseComplement(DNAString("GCGTTACCTTGCAGGAATTGAGGCCGTCTGTTAATTTCC")) 
H5_TA=reverseComplement(DNAString("GCGTTACCTTGCAGGAATTGAGGCCGTCCGTTAATATCC")) 
H5=reverseComplement(DNAString("GCGTTACCTTGCAGGAATTGAGGCCGTCTGTTAATATCC")) 
p0R2=reverseComplement(DNAString("GCGTTACCTTGCAGGAATTGAGGCCGTCCGTTAATTTCC"))
galK=reverseComplement(DNAString("atgagtctgaaagaaaaaacacaatctctgtttgccaac"))

compareStrings(H5,p0R2) #2SNPs

#count gap pattern in multiple reads
g_S8=vcountPattern(gap,sread(S8)) #vcountPattern (works on DNAstring set)
g_S9=vcountPattern(gap,sread(S9))
g_S10=vcountPattern(gap,sread(S10))
g_S11=vcountPattern(gap,sread(S11))
g_S12=vcountPattern(gap,sread(S12))
g_S13=vcountPattern(gap,sread(S13))

#length of read sets compared to gap and anc
length(g_S8)
length(g_S9)
length(g_S10)
length(g_S11)
length(g_S12)
length(g_S13)

gaplist=list(g_S8,g_S9,g_S10,g_S11,g_S12,g_S13)

#count p0-2 anc pattern in multiple reads
a_S8=vcountPattern(anc,sread(S8)) #vcountPattern (works on DNAstring set)
a_S9=vcountPattern(anc,sread(S9))
a_S10=vcountPattern(anc,sread(S10))
a_S11=vcountPattern(anc,sread(S11))
a_S12=vcountPattern(anc,sread(S12))
a_S13=vcountPattern(anc,sread(S13))

#length of read sets compared to gap and anc
length(a_S8)
length(a_S9)
length(a_S10)
length(a_S11)
length(a_S12)
length(a_S13)

anclist=list(a_S8,a_S9,a_S10,a_S11,a_S12,a_S13)

#count reads with p0 pattern
p_S8=vcountPattern(p0,sread(S8)) #vcountPattern (works on DNAstring set)
p_S9=vcountPattern(p0,sread(S9))
p_S10=vcountPattern(p0,sread(S10))
p_S11=vcountPattern(p0,sread(S11))
p_S12=vcountPattern(p0,sread(S12))
p_S13=vcountPattern(p0,sread(S13))

p0list=list(p_S8,p_S9,p_S10,p_S11,p_S12,p_S13)



## for both gap and anc pattern counting the SD is comparable to the mean --> random data. 

### gap/anc+gap, otherwise IÂ´m also counting p0

#count gaps and ancs and p0s - i.e. 1s in  boolean matrices (loop)
GA_results=NULL
absolute_GA_results=matrix(0,nrow=3,ncol=6)
for (i in 1: 6) {
  g=length(which(unlist(gaplist[i]) ==1)) #all reads with p0-2 gaps
  a=length(which(unlist(anclist[i]) ==1)) #all reads with p0-2 ancs
  p=length(which(unlist(p0list[i])==1))#all reads with p0 (ancestral) pattern
   GA_results[i]<- g/(a+g)
absolute_GA_results[1,i]=g
absolute_GA_results[2,i]=a+g
absolute_GA_results[3,i]=p
 }

plot(GA_results)
mean(GA_results[1:3]) #0.37 IS+
mean(GA_results[4:6]) #0.60 IS-   % of p0-2 with anc sequence pattern at the gap region

sd(GA_results[1:3]) # IS+
sd(GA_results[4:6]) # IS-

dev.off()
barplot(absolute_GA_results,beside=TRUE,col=c(" green","dark green","grey"),legend=c("gap", "p02","p0"),args.legend = list(x = "topleft"),names.arg=c("IS+","IS+","IS+","IS-","IS-","IS-"))


library(matrixStats)

df <- data.frame(bar =c(mean(GA_results[1:3]),mean(GA_results[4:6])) , error = c(sd(GA_results[1:3])),sd(GA_results[4:6]))
    foo <- barplot(df$bar,border=NA, ylim=c(0,1), names.arg=c("IS+","IS-"),ylab="gap reads/ancestral reads") #arrows become errorbars stackoverflow
arrows(x0=foo,y0=df$bar+df$error,y1=df$bar-df$error,angle=90,code=3,length=0.1)
                         
#############################################
#############################################

###test other amplicons for the prevalence of gap, anc(p02) and p0
#test also the p0+p02 pooled samples from above!
### fw reads R1
#ctr E
S1=S1_S1_R1_001.fastq.gz
# EE24CIS (IS-)
S2=S2_S2_R1_001.fastq.gz
S3=S3_S3_R1_001.fastq.gz
S4=S4_S4_R1_001.fastq.gz
#EE24C30 (IS+)
S5=S5_S5_R1_001.fastq.gz
S6=S6_S6_R1_001.fastq.gz
S7=S7_S7_R1_001.fastq.gz
#EE24A30 (IS+)
S14=S14_S14_R1_001.fastq.gz
#EE24AIS (IS-)
S15=S15_S15_R1_001.fastq.gz

##reverse reads R2
#ctr E
S1r=S1_S1_R2_001.fastq.gz
# EE24CIS (IS-)
S2r=S2_S2_R2_001.fastq.gz
S3r=S3_S3_R2_001.fastq.gz
S4r=S4_S4_R2_001.fastq.gz
#EE24C30 (IS+)
S5r=S5_S5_R2_001.fastq.gz
S6r=S6_S6_R2_001.fastq.gz
S7r=S7_S7_R2_001.fastq.gz
#EE25B (IS+)
S8r=S8_S8_R2_001.fastq.gz
S9r=S9_S9_R2_001.fastq.gz
S10r=S10_S10_R2_001.fastq.gz
#EE25B (IS-)
S11r=S11_S11_R2_001.fastq.gz
S12r=S12_S12_R2_001.fastq.gz
S13r=S13_S13_R2_001.fastq.gz
#EE24A30 (IS+)
S14r=S14_S14_R2_001.fastq.gz
#EE24AIS (IS-)
S15r=S15_S15_R2_001.fastq.gz



#count pattern p0
p_S1=vcountPattern(p0,sread(S1)) #vcountPattern (works on DNAstring set)
p_S2=vcountPattern(p0,sread(S2))
p_S3=vcountPattern(p0,sread(S3))
p_S4=vcountPattern(p0,sread(S4))
p_S5=vcountPattern(p0,sread(S5))
p_S6=vcountPattern(p0,sread(S6))
p_S7=vcountPattern(p0,sread(S7)) #vcountPattern (works on DNAstring set)
p_S8=vcountPattern(p0,sread(S8))
p_S9=vcountPattern(p0,sread(S9))
p_S10=vcountPattern(p0,sread(S10))
p_S11=vcountPattern(p0,sread(S11))
p_S12=vcountPattern(p0,sread(S12))
p_S13=vcountPattern(p0,sread(S13))
p_S14=vcountPattern(p0,sread(S14))
p_S15=vcountPattern(p0,sread(S15))
p0list=list(p_S1,p_S2,p_S3,p_S4,p_S5,p_S6,p_S7,p_S14,p_S15)
p0list_all=list(p_S1,p_S2,p_S3,p_S4,p_S5,p_S6,p_S7,p_S8,p_S9,p_S10,p_S11,p_S12,p_S13,p_S14,p_S15)

#count pattern anc (=ancestral p02)
a_S1=vcountPattern(anc,sread(S1)) #vcountPattern (works on DNAstring set)
a_S2=vcountPattern(anc,sread(S2))
a_S3=vcountPattern(anc,sread(S3))
a_S4=vcountPattern(anc,sread(S4))
a_S5=vcountPattern(anc,sread(S5))
a_S6=vcountPattern(anc,sread(S6))
a_S7=vcountPattern(anc,sread(S7)) #vcountPattern (works on DNAstring set)
a_S8=vcountPattern(anc,sread(S8))
a_S9=vcountPattern(anc,sread(S9))
a_S10=vcountPattern(anc,sread(S10))
a_S11=vcountPattern(anc,sread(S11))
a_S12=vcountPattern(anc,sread(S12))
a_S13=vcountPattern(anc,sread(S13))
a_S14=vcountPattern(anc,sread(S14))
a_S15=vcountPattern(anc,sread(S15))
anclist=list(a_S1,a_S2,a_S3,a_S4,a_S5,a_S6,a_S7,a_S14,a_S15)
anclist_all=list(a_S1,a_S2,a_S3,a_S4,a_S5,a_S6,a_S7,a_S8,a_S9,a_S10,a_S11,a_S12,a_S13,a_S14,a_S15)

#count pattern gap
g_S1=vcountPattern(gap,sread(S1)) #vcountPattern (works on DNAstring set)
g_S2=vcountPattern(gap,sread(S2))
g_S3=vcountPattern(gap,sread(S3))
g_S4=vcountPattern(gap,sread(S4))
g_S5=vcountPattern(gap,sread(S5))
g_S6=vcountPattern(gap,sread(S6))
g_S7=vcountPattern(gap,sread(S7)) #vcountPattern (works on DNAstring set)
g_S8=vcountPattern(gap,sread(S8))
g_S9=vcountPattern(gap,sread(S9))
g_S10=vcountPattern(gap,sread(S10))
g_S11=vcountPattern(gap,sread(S11))
g_S12=vcountPattern(gap,sread(S12))
g_S13=vcountPattern(gap,sread(S13))
g_S14=vcountPattern(gap,sread(S14))
g_S15=vcountPattern(gap,sread(S15))
gaplist=list(g_S1,g_S2,g_S3,g_S4,g_S5,g_S6,g_S7,g_S14,g_S15)
gaplist_all=list(g_S1,g_S2,g_S3,g_S4,g_S5,g_S6,g_S7,g_S8,g_S9,g_S10,g_S11,g_S12,g_S13,g_S14,g_S15)


##count gaps and ancs and p0s - i.e. 1s in  boolean matrices (loop)
GA_results=NULL
absolute_GA_results=matrix(0,nrow=3,ncol=length(gaplist_all))
for (i in 1: length(gaplist_all)) {
  g=length(which(unlist(gaplist_all[i]) ==1)) #all reads with p0-2 gaps
  a=length(which(unlist(anclist_all[i]) ==1)) #all reads with p0-2 ancs
  p=length(which(unlist(p0list_all[i])==1))#all reads with p0 (ancestral) pattern
  GA_results[i]<- g/(a+g)
  absolute_GA_results[1,i]=g #nr of gaps in p02 reads
  absolute_GA_results[2,i]=a+g #nr of gaps+ancestral p02 reads
  absolute_GA_results[3,i]=p #nr of ancestral p0 reads
}

barplot(absolute_GA_results,beside=TRUE,col=c(" green","dark green","grey"),legend=c("gap", "p02","p0"),args.legend = list(x = "topleft"),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"))
#This shows that except for S5,S6,S7 (EE24C30), all samples are highly contaminated with p02 reads
#moreover, all sequences have high amounts of gap reads.

#Figure 5 Supplement/Supplementary Figure 4
#number of reads and contamination issue
barplot(absolute_GA_results[2:3,],beside=TRUE,col=c("dark green","grey"),legend=c( "p02","p0"),args.legend = list(x = "topleft"),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"))

#add number of total reads (lengths)
lengths=c(length(S1),length(S2),length(S3),length(S4),length(S5),length(S6),length(S7),length(S8),length(S9),length(S10),length(S11),length(S12),length(S13),length(S14),length(S15))

#plot number of total reads, nr of p02 and p0 
absolute_GA_results[1,]=lengths
barplot(log(absolute_GA_results),beside=TRUE,col=c("grey","dark green","green"),legend=c( "all reads","p02","p0"),args.legend = list(x = "bottomright"),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"))

##reverse reads (R2) matching p0 or H5
#count pattern p0
p_S1r=vcountPattern(p0R2,sread(S1r)) #vcountPattern (works on DNAstring set)
p_S2r=vcountPattern(p0R2,sread(S2r))
p_S3r=vcountPattern(p0R2,sread(S3r))
p_S4r=vcountPattern(p0R2,sread(S4r))
p_S5r=vcountPattern(p0R2,sread(S5r))
p_S6r=vcountPattern(p0R2,sread(S6r))
p_S7r=vcountPattern(p0R2,sread(S7r)) #vcountPattern (works on DNAstring set)
p_S8r=vcountPattern(p0R2,sread(S8r))
p_S9r=vcountPattern(p0R2,sread(S9r))
p_S10r=vcountPattern(p0R2,sread(S10r))
p_S11r=vcountPattern(p0R2,sread(S11r))
p_S12r=vcountPattern(p0R2,sread(S12r))
p_S13r=vcountPattern(p0R2,sread(S13r))
p_S14r=vcountPattern(p0R2,sread(S14r))
p_S15r=vcountPattern(p0R2,sread(S15r))
p0R2list_all=list(p_S1r,p_S2r,p_S3r,p_S4r,p_S5r,p_S6r,p_S7r,p_S8r,p_S9r,p_S10r,p_S11r,p_S12r,p_S13r,p_S14r,p_S15r)

h5_S1r=vcountPattern(H5,sread(S1r)) #vcountPattern (works on DNAstring set)
h5_S2r=vcountPattern(H5,sread(S2r))
h5_S3r=vcountPattern(H5,sread(S3r))
h5_S4r=vcountPattern(H5,sread(S4r))
h5_S5r=vcountPattern(H5,sread(S5r))
h5_S6r=vcountPattern(H5,sread(S6r))
h5_S7r=vcountPattern(H5,sread(S7r)) #vcountPattern (works on DNAstring set)
h5_S8r=vcountPattern(H5,sread(S8r))
h5_S9r=vcountPattern(H5,sread(S9r))
h5_S10r=vcountPattern(H5,sread(S10r))
h5_S11r=vcountPattern(H5,sread(S11r))
h5_S12r=vcountPattern(H5,sread(S12r))
h5_S13r=vcountPattern(H5,sread(S13r))
h5_S14r=vcountPattern(H5,sread(S14r))
h5_S15r=vcountPattern(H5,sread(S15r))
H5list_all=list(h5_S1r,h5_S2r,h5_S3r,h5_S4r,h5_S5r,h5_S6r,h5_S7r,h5_S8r,h5_S9r,h5_S10r,h5_S11r,h5_S12r,h5_S13r,h5_S14r,h5_S15r)

h5_TA_S1r=vcountPattern(H5_TA,sread(S1r)) #vcountPattern (works on DNAstring set)
h5_TA_S2r=vcountPattern(H5_TA,sread(S2r))
h5_TA_S3r=vcountPattern(H5_TA,sread(S3r))
h5_TA_S4r=vcountPattern(H5_TA,sread(S4r))
h5_TA_S5r=vcountPattern(H5_TA,sread(S5r))
h5_TA_S6r=vcountPattern(H5_TA,sread(S6r))
h5_TA_S7r=vcountPattern(H5_TA,sread(S7r)) #vcountPattern (works on DNAstring set)
h5_TA_S8r=vcountPattern(H5_TA,sread(S8r))
h5_TA_S9r=vcountPattern(H5_TA,sread(S9r))
h5_TA_S10r=vcountPattern(H5_TA,sread(S10r))
h5_TA_S11r=vcountPattern(H5_TA,sread(S11r))
h5_TA_S12r=vcountPattern(H5_TA,sread(S12r))
h5_TA_S13r=vcountPattern(H5_TA,sread(S13r))
h5_TA_S14r=vcountPattern(H5_TA,sread(S14r))
h5_TA_S15r=vcountPattern(H5_TA,sread(S15r))
H5TAlist_all=list(h5_TA_S1r,h5_TA_S2r,h5_TA_S3r,h5_TA_S4r,h5_TA_S5r,h5_TA_S6r,h5_TA_S7r,h5_TA_S8r,h5_TA_S9r,h5_TA_S10r,h5_TA_S11r,h5_TA_S12r,h5_TA_S13r,h5_TA_S14r,h5_TA_S15r)

h5_CT_S1r=vcountPattern(H5_CT,sread(S1r)) #vcountPattern (works on DNAstring set)
h5_CT_S2r=vcountPattern(H5_CT,sread(S2r))
h5_CT_S3r=vcountPattern(H5_CT,sread(S3r))
h5_CT_S4r=vcountPattern(H5_CT,sread(S4r))
h5_CT_S5r=vcountPattern(H5_CT,sread(S5r))
h5_CT_S6r=vcountPattern(H5_CT,sread(S6r))
h5_CT_S7r=vcountPattern(H5_CT,sread(S7r))
h5_CT_S8r=vcountPattern(H5_CT,sread(S8r))
h5_CT_S9r=vcountPattern(H5_CT,sread(S9r))
h5_CT_S10r=vcountPattern(H5_CT,sread(S10r))
h5_CT_S11r=vcountPattern(H5_CT,sread(S11r))
h5_CT_S12r=vcountPattern(H5_CT,sread(S12r))
h5_CT_S13r=vcountPattern(H5_CT,sread(S13r))
h5_CT_S14r=vcountPattern(H5_CT,sread(S14r))
h5_CT_S15r=vcountPattern(H5_CT,sread(S15r))
H5CTlist_all=list(h5_CT_S1r,h5_CT_S2r,h5_CT_S3r,h5_CT_S4r,h5_CT_S5r,h5_CT_S6r,h5_CT_S7r,h5_CT_S8r,h5_CT_S9r,h5_CT_S10r,h5_CT_S11r,h5_CT_S12r,h5_CT_S13r,h5_CT_S14r,h5_CT_S15r)

#galK wt patterns (39 first nucleotides ATG...)
galK_wt_S1r=vcountPattern(galK,sread(S1r),max.mismatch=0) #vcountPattern (works on DNAstring set)
galK_wt_S2r=vcountPattern(galK,sread(S2r),max.mismatch=0)
galK_wt_S3r=vcountPattern(galK,sread(S3r),max.mismatch=0)
galK_wt_S4r=vcountPattern(galK,sread(S4r),max.mismatch=0)
galK_wt_S5r=vcountPattern(galK,sread(S5r),max.mismatch=0)
galK_wt_S6r=vcountPattern(galK,sread(S6r),max.mismatch=0)
galK_wt_S7r=vcountPattern(galK,sread(S7r),max.mismatch=0)
galK_wt_S8r=vcountPattern(galK,sread(S8r),max.mismatch=0)
galK_wt_S9r=vcountPattern(galK,sread(S9r),max.mismatch=0)
galK_wt_S10r=vcountPattern(galK,sread(S10r),max.mismatch=0)
galK_wt_S11r=vcountPattern(galK,sread(S11r),max.mismatch=0)
galK_wt_S12r=vcountPattern(galK,sread(S12r),max.mismatch=0)
galK_wt_S13r=vcountPattern(galK,sread(S13r),max.mismatch=0)
galK_wt_S7r=vcountPattern(galK,sread(S7r),max.mismatch=0) #vcountPattern (works on DNAstring set)
galK_wt_S14r=vcountPattern(galK,sread(S14r),max.mismatch=0)
galK_wt_S15r=vcountPattern(galK,sread(S15r),max.mismatch=0)
galKwtlist_all=list(galK_wt_S1r,galK_wt_S2r,galK_wt_S3r,galK_wt_S4r,galK_wt_S5r,galK_wt_S6r,galK_wt_S7r,galK_wt_S8r,galK_wt_S9r,galK_wt_S10r,galK_wt_S11r,galK_wt_S12r,galK_wt_S13r,galK_wt_S14r,galK_wt_S15r)

#galK 1SNP anywhere patterns (39 first nucleotides ATG...)
galK_1SNP_S1r=vcountPattern(galK,sread(S1r),max.mismatch=1) #vcountPattern (works on DNAstring set)
galK_1SNP_S2r=vcountPattern(galK,sread(S2r),max.mismatch=1)
galK_1SNP_S3r=vcountPattern(galK,sread(S3r),max.mismatch=1)
galK_1SNP_S4r=vcountPattern(galK,sread(S4r),max.mismatch=1)
galK_1SNP_S5r=vcountPattern(galK,sread(S5r),max.mismatch=1)
galK_1SNP_S6r=vcountPattern(galK,sread(S6r),max.mismatch=1)
galK_1SNP_S7r=vcountPattern(galK,sread(S7r),max.mismatch=1)
galK_1SNP_S8r=vcountPattern(galK,sread(S8r),max.mismatch=1)
galK_1SNP_S9r=vcountPattern(galK,sread(S9r),max.mismatch=1)
galK_1SNP_S10r=vcountPattern(galK,sread(S10r),max.mismatch=1)
galK_1SNP_S11r=vcountPattern(galK,sread(S11r),max.mismatch=1)
galK_1SNP_S12r=vcountPattern(galK,sread(S12r),max.mismatch=1)
galK_1SNP_S13r=vcountPattern(galK,sread(S13r),max.mismatch=1)
galK_1SNP_S7r=vcountPattern(galK,sread(S7r),max.mismatch=1) #vcountPattern (works on DNAstring set)
galK_1SNP_S14r=vcountPattern(galK,sread(S14r),max.mismatch=1)
galK_1SNP_S15r=vcountPattern(galK,sread(S15r),max.mismatch=1)
galK1SNPlist_all=list(galK_1SNP_S1r,galK_1SNP_S2r,galK_1SNP_S3r,galK_1SNP_S4r,galK_1SNP_S5r,galK_1SNP_S6r,galK_1SNP_S7r,galK_1SNP_S8r,galK_1SNP_S9r,galK_1SNP_S10r,galK_1SNP_S11r,galK_1SNP_S12r,galK_1SNP_S13r,galK_1SNP_S14r,galK_1SNP_S15r)

##count H5 and p0 - i.e. 1s in  boolean matrices (loop)
GA_results=matrix(0,nrow=4,ncol=length(p0R2list_all))
absolute_GA_results=matrix(0,nrow=7,ncol=length(p0R2list_all))
for (i in 1: length(p0R2list_all)) {
  p=length(which(unlist(p0R2list_all[i]) ==1)) #all reads with p0 
  h=length(which(unlist(H5list_all[i]) ==1)) #all reads with H5
  hta=length(which(unlist(H5TAlist_all[i]) ==1)) #all reads with H5-TA mutation
  hct=length(which(unlist(H5CTlist_all[i]) ==1)) #all reads with H5-CT mutation
  galKwt=length(which(unlist(galKwtlist_all[i])==1))
  galK1SNP=length(which(unlist(galK1SNPlist_all[i])==1))
  GA_results[1,i]=h/p
  GA_results[2,i]=hta/p
  GA_results[3,i]=hct/p
  GA_results[4,i]=galK1SNP/galKwt/39
  absolute_GA_results[1,i]=p #number of ancestral p0
  absolute_GA_results[2,i]=h #number of h5r SNP
  absolute_GA_results[3,i]=hta #nr of h5r-T>A SNP
  absolute_GA_results[4,i]=hct #nr of h5r-C>T SNP
  absolute_GA_results[5,i]=galK1SNP/39 #mean nr of reads with 1SNPs in galK (39bp pattern)
  absolute_GA_results[6,i]=galKwt #nr of galKwt reads (no SNPs for 39bp pattern)
  absolute_GA_results[7,i]=length(unlist(p0R2list_all[i])) #nr all reads
}

#add +1 to allow logplot (some have 0 H5r instances -> now 1)
barplot(absolute_GA_results+1,beside=TRUE,col=c("grey","green","dark green", "blue","orange","red","black"),legend=c("p0", "H5","TA","CT","1SNPgalK","galK","all reads"),args.legend = list(x = "right"),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"),log="y")
barplot(absolute_GA_results+1,beside=TRUE,col=c("grey","green","dark green", "blue","orange","red","black"),args.legend = list(x = "right"),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"),log="y")
#legend=c("p0", "H5","TA","CT","1SNPgalK","galK","all reads"),
#non-log
barplot(absolute_GA_results,beside=TRUE,col=c("grey","green","dark green", "blue","orange","red","black"),args.legend = list(x = "right"),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"))
#legend=c("p0", "H5","TA","CT","1SNPgalK","galK","all reads"),

barplot(GA_results,names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"), main="SNPs/p0")
barplot((absolute_GA_results[3,])/(absolute_GA_results[1,]+1),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"), main="T>A SNPs/p0")
barplot((absolute_GA_results[2,])/(absolute_GA_results[1,]+1),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"), main="H5r SNPs/p0")
barplot((absolute_GA_results[4,])/(absolute_GA_results[1,]+1),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"), main="C>T SNPs/p0")
barplot(GA_results[3,],names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"), main="C>T SNPs/p0")#same plot as above (check ok)
#all SNPs relative to p0 -nice
barplot(GA_results[1:4,],beside=TRUE,col=c("green","dark green","blue","red"),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"))
#legend=c("H5r",TA","CT","mean1SNP")

# mean frequency of a single SNP/all galK reads
barplot(GA_results[4,],beside=TRUE,col=c("grey"),args.legend = list(x = "right"),names.arg=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15"))
#> by chance we expect 0.028 of all reads to have a certain SNP 
mean(GA_results[4,])

#this does not normalize for multiple p0 copies with one SNP though

### SUMMARY PLOTS WITH ERRORBARS
library(matrixStats)

## plot for EE24 *C*IS and *C*30
#with ctr #of galK 1SNPs over all samples   ...mean(GA_results[4,1:15])
##H5r SNPs (very low overall nr. = 1/100 of T>A)
df <- data.frame(bar =c(mean(GA_results[4,1:15]),mean(GA_results[1,2:4]),mean(GA_results[1,5:7])) , error = c(sd(GA_results[4,1:15]),sd(GA_results[1,2:4]),sd(GA_results[1,5:7])))
foo <- barplot(df$bar,border=NA,ylim=c(0,0.05), names.arg=c("ctr1SNP","IS-","IS+"),ylab="H5r mutations/ancestral p0") #arrows become errorbars stackoverflow
arrows(x0=foo,y0=df$bar+df$error,y1=df$bar-df$error,angle=90,code=3,length=0.1)
##T>A SNPs 
df <- data.frame(bar =c(mean(GA_results[4,1:15]),mean(GA_results[2,2:4]),mean(GA_results[2,5:7])) , error = c(sd(GA_results[4,1:15]),sd(GA_results[2,2:4]),sd(GA_results[2,5:7])))
foo <- barplot(df$bar,border=NA, names.arg=c("ctr1SNP","IS-","IS+"),ylim=c(0,0.15),ylab="T>A mutations/ancestral p0") #arrows become errorbars stackoverflow
arrows(x0=foo,y0=df$bar+df$error,y1=df$bar-df$error,angle=90,code=3,length=0.1)
##C>T SNPs 
df <- data.frame(bar =c(mean(GA_results[4,1:15]),mean(GA_results[3,2:4]),mean(GA_results[3,5:7])) , error = c(sd(GA_results[4,1:15]),sd(GA_results[3,2:4]),sd(GA_results[3,5:7])))
foo <- barplot(df$bar,border=NA,ylim=c(0,0.15), names.arg=c("ctr1SNP","IS-","IS+"),ylab="C>T mutations/ancestral p0") #arrows become errorbars stackoverflow
arrows(x0=foo,y0=df$bar+df$error,y1=df$bar-df$error,angle=90,code=3,length=0.1)


## plot for EE25*B* IS+ and IS-  ... no significant difference
##H5r SNPs (very low overall, none in IS-)
df <- data.frame(bar =c(mean(GA_results[1,8:10]),mean(GA_results[1,11:13])) , error = c(sd(GA_results[1,7:9])),sd(GA_results[1,10:13]))
foo <- barplot(df$bar,border=NA,ylim=c(0,0.0001), names.arg=c("IS+","IS-"),ylab="H5r mutations/ancestral p0") #arrows become errorbars stackoverflow
arrows(x0=foo,y0=df$bar+df$error,y1=df$bar-df$error,angle=90,code=3,length=0.1)
##T>A SNPs 
df <- data.frame(bar =c(mean(GA_results[2,8:10]),mean(GA_results[2,11:13])) , error = c(sd(GA_results[2,7:9])),sd(GA_results[2,10:13]))
foo <- barplot(df$bar,border=NA, names.arg=c("IS+","IS-"),ylim=c(0,0.08),ylab="T>A mutations/ancestral p0") #arrows become errorbars stackoverflow
arrows(x0=foo,y0=df$bar+df$error,y1=df$bar-df$error,angle=90,code=3,length=0.1)
##C>T SNPs 
df <- data.frame(bar =c(mean(GA_results[3,8:10]),mean(GA_results[3,11:13])) , error = c(sd(GA_results[3,7:9])),sd(GA_results[3,10:13]))
foo <- barplot(df$bar,border=NA,ylim=c(0,0.008), names.arg=c("IS+","IS-"),ylab="C>T mutations/ancestral p0") #arrows become errorbars stackoverflow
arrows(x0=foo,y0=df$bar+df$error,y1=df$bar-df$error,angle=90,code=3,length=0.1)

## plot for EE24*A* IS+ and IS-  ... just 1 sample each 
#with ctr no gal (E)
##H5r SNPs  
df <- data.frame(bar =c((GA_results[1,1]),(GA_results[1,14]),(GA_results[1,15])) )
foo <- barplot(df$bar,border=NA,ylim=c(0,0.0004), names.arg=c("ctr","IS+","IS-"),ylab="H5r mutations/ancestral p0") #arrows become errorbars stackoverflow
##T>A SNPs 
df <- data.frame(bar =c(GA_results[2,1],(GA_results[2,14]),(GA_results[2,15])) )
foo <- barplot(df$bar,border=NA, names.arg=c("ctr","IS+","IS-"),ylim=c(0,0.2),ylab="T>A mutations/ancestral p0") #arrows become errorbars stackoverflow
##C>T SNPs 
df <- data.frame(bar =c(GA_results[3,1],(GA_results[3,14]),(GA_results[3,15])) )
foo <- barplot(df$bar,border=NA,ylim=c(0,0.005), names.arg=c("ctr","IS+","IS-"),ylab="C>T mutations/ancestral p0") #arrows become errorbars stackoverflow


##Plot all of these together - paper plots
library(ggplot2) 
library(reshape2)#for melt
library(plyr)#for summary stats
?ddply()

# (to see the full errorbars plot from -0.01-0.15)

plot_df=(GA_results)
colnames(plot_df)=c("E","CIS-","CIS-","CIS-","CIS+","CIS+","CIS+","BIS+","BIS+","BIS+","BIS-","BIS-","BIS-","AIS+","AIS-")
rownames(plot_df)=c("H5r","TA","CT","1SNP")
plot_df=melt(plot_df)
colnames(plot_df)=c("mutation","condition","frequency")
sum <- ddply(plot_df,c("mutation","condition"),summarise,mfrequency=mean(frequency),sdf=sd(frequency,na.rm=TRUE))

#all points
ggplot()+
  geom_point(data=plot_df,aes(x=condition,y=frequency,col=mutation))+
  theme_bw()+ ylim(0,0.15)+ ggtitle("read frequency")+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))

#all bars
ggplot(data=plot_df,aes(x=(condition),y=frequency,col=mutation),group=condition)+
  geom_bar(stat="identity", position=position_dodge2())+
  #geom_text(aes(label=round(frequency,digits=2)), vjust=-0.3, size=3.5)+
  theme_bw()+ ggtitle("read frequency")+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))

#bars sum stats  
ggplot(data=sum,aes(x=(condition),y=mfrequency,col=mutation,fill=mutation))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(data=na.omit(sum),aes(x=condition,y=mfrequency,ymin=mfrequency-sdf,ymax=mfrequency+sdf),position=position_dodge(),col="black")+
 # geom_text(aes(label=round(mfrequency,digits=2)), vjust=-0.3, size=3.5)+
  theme_bw()+ ylim(-0.01,0.15)+ ggtitle("read frequency")+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))

#all points sum stats
ggplot(data=sum,aes(x=(condition),y=mfrequency,col=mutation,fill=mutation))+
  geom_point()+ 
  geom_errorbar(data=(sum),aes(x=condition,ymin=mfrequency-sdf,ymax=mfrequency+sdf),col="grey",width=0.2)+
  theme_bw()+ggtitle("read frequency")+ylim(-0.01,0.15)+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))

#bars sum stats - mutations 
ggplot(data=subset(sum,mutation=="TA"),aes(x=(condition),y=mfrequency,col=mutation,fill=mutation))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(x=condition,y=mfrequency,ymin=mfrequency-sdf,ymax=mfrequency+sdf),position=position_dodge(0.9),width=0.2,col="black")+
  # geom_text(aes(label=round(mfrequency,digits=2)), vjust=-0.3, size=3.5)+
  theme_bw()+ ggtitle("read frequency")+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))


#Figure 5 A plots
#bars sum stats - mutations 
ggplot(data=subset(sum,condition=="CIS-"|condition=="CIS+"),aes(x=(condition),y=mfrequency,col=mutation,fill=mutation))+
  geom_bar(stat="identity", position=position_dodge(),col="black")+
  geom_errorbar(aes(x=condition,y=mfrequency,ymin=mfrequency-sdf,ymax=mfrequency+sdf),position=position_dodge(0.9),width=0.2,col="black")+
  # geom_text(aes(label=round(mfrequency,digits=2)), vjust=-0.3, size=3.5)+
  ylim(-0.010,0.15)+theme_bw()+ ggtitle("read frequency")+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))


#Figure 5 B plots
#bars sum stats - mutations 
ggplot(data=subset(sum,condition=="BIS-"|condition=="BIS+"),aes(x=(condition),y=mfrequency,col=mutation,fill=mutation))+
  geom_bar(stat="identity", position=position_dodge(),col="black")+
  geom_errorbar(aes(x=condition,y=mfrequency,ymin=mfrequency-sdf,ymax=mfrequency+sdf),position=position_dodge(0.9),width=0.2,col="black")+
  # geom_text(aes(label=round(mfrequency,digits=2)), vjust=-0.3, size=3.5)+
  ylim(-0.01,0.15)+theme_bw()+ ggtitle("read frequency")+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))


#controls
ggplot(data=subset(sum,condition=="E"),aes(x=(condition),y=mfrequency,col=mutation,fill=mutation))+
  geom_bar(stat="identity", position=position_dodge(),col="black")+
  geom_errorbar(aes(x=condition,y=mfrequency,ymin=mfrequency-sdf,ymax=mfrequency+sdf),position=position_dodge(),col="black")+
  # geom_text(aes(label=round(mfrequency,digits=2)), vjust=-0.3, size=3.5)+
  ylim(-0.01,0.15)+theme_bw()+ ggtitle("read frequency")+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))


