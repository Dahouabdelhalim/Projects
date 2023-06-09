# Version 3.2 
# 10/23/2011
rm(list = ls())
gc()

## Set parameters for run ##

params <- list(
  run.name = "Bphy-snps",
  folder.name = "~/FASTQ_by_species/Bphy",
  ref.fname = "Bphy-genbank-SNP-sequences.fa",
  pipeline = "snp",
  aligner = "bwa",
  num.cores = 2
)
# run.name - Unique name for the run used to name folder that results will be stored in.
# folder.name - Full path location of fastq NGS reads.
# ref.fname - Full path of fasta file containing reference sequence(s).
# pipeline - Which pipeline to run: "mtdna" or "snp".
# aligner - Which alignment engine to use: "bwa", "soap", "bowtie".
# num.cores - Number of CPUs to distribute runs for each sample over (also sets default for 'bwa.num.threads').


qa.params <- list(qa.run = F, qa.q1.trim.qual = 5, qa.min.length = 10, qa.filter.qual = 5, qa.filter.pct = 95, qa.mask.qual = 2)
# qa.run : Do QA filtering and trimming (default = F)
# qa.q1.trim.qual : Trim all reads to the last cycle that has its 1st quartile >= q1.trim.qual (default = 10)
# qa.min.length (-l) : Remove all trimmed reads that are shorter than min.length (default = 20)
# qa.filter.qual (-q), qa.filter.pct (-p) : Remove any read that does not have at least filter.pct 
#     of its nucleotides >= filter.qual (default = 15, 50)
# qa.mask.qual (-q) : Change any remaining nucleotide <= 'mask.qual' to 'N' (default = 10)


align.params <- list(bwa.max.edit.dist = 0.04, bwa.max.gap.opens = 3, 
  bwa.max.gap.extensions = 1, bwa.deletion.disallow = 16, bwa.indel.disallow = 5, 
  bwa.seed.length = 35, bwa.max.seed.edit.distance = 2, bwa.num.threads = params$num.cores, 
  bwa.mismatch.penalty = 3, bwa.gap.open.penalty = 11, bwa.gap.extension.penalty = 4, 
  bwa.read.trimming = 0, remove.duplicate.reads = T
)   
# bwa.max.edit.dist (-n) : Maximum edit distance if the value is INT, or the fraction of missing 
#     alignments given 2% uniform base error rate if FLOAT. In the latter case, the 
#     maximum edit distance is automatically chosen for different read lengths. [0.04]
# bwa.max.gap.opens (-o) : Maximum number of gap opens [1]
# bwa.max.gap.exentsions (-e) : Maximum number of gap extensions, -1 for 
#     k-difference mode (disallowing long gaps) [-1]
# bwa.deletion.disallow (-d) : Disallow a long deletion within INT bp towards the 3’-end [16]
# bwa.indel.disallow (-i) : Disallow an indel within INT bp towards the ends [5]
# bwa.seed.length (-l) : Take the first INT subsequence as seed. If INT is larger than the query 
#     sequence, seeding will be disabled. For long reads, this option is typically 
#     ranged from 25 to 35 for ‘-k 2’. [inf]
# bwa.max.seed.edit.distance (-k) : Maximum edit distance in the seed [2]
# bwa.num.threads (-t) : Number of threads (multi-threading mode) [1]
# bwa.mismatch.penalty (-M) : Mismatch penalty. BWA will not search for suboptimal hits with
#     a score lower than (bestScore-misMsc). [3]
# bwa.gap.open.penalty (-O) : Gap open penalty [11]
# bwa.gap.extension.penalty (-E) : Gap extension penalty [4] 
# bwa.read.trimming (-q) : Parameter for read trimming. BWA trims a read down to 
#     argmax_x{\\sum_{i=x+1}^l(INT-q_i)} if q_l<INT where l is the original read length. [0]
# remove.duplicate.reads : If (T)rue, BAM file will only contain unique reads mapped to reference.


realign.params <- list(ra.min.reads = 4, ra.window.size = 10, ra.mismatch.frac = 0.15, 
  ra.max.interval.size = 500, ra.LOD = 5
)
#  ra.min.reads (--minReadsAtLocus) : The minimum coverage at a locus for the entropy calculation to be enabled; default=4
#  ra.window.size (--windowSize) : Any two SNP calls and/or high entropy positions are considered clustered when they occur 
#      no more than N basepairs apart; default=10
#  ra.mismatch.frac (--mismatchFraction) : Fraction of total sum of base qualities at a position that need to mismatch for the 
#      position to be considered to have high entropy; default=0.15; to disable, set to <= 0 or > 1
#      Note that this fraction should be adjusted based on your particular data set. For deep coverage 
#      and/or when looking for indels with low allele frequency, this number should be smaller.
#  ra.max.interval.size (--maxIntervalSize) : Max size in bp of intervals that we'll pass to the realigner; default=500
#      Because the realignment algorithm is N^2, allowing too large an interval might take too long 
#      to completely realign.
#  ra.LOD (-LOD) : LOD threshold above which the realigner will proceed to realign; default=5.0
#      This term is equivalent to "significance" - i.e. is the improvement significant enough to 
#      merit realignment? Note that this number should be adjusted based on your particular data set.
#      For low coverage and/or when looking for indels with low allele frequency, this number should be smaller.


summary.params <- list(smry.min.map.qual = 0, smry.min.base.qual = 0)
# smry.min.map.qual - Ignores reads with mapping quality below x (not counted in coverage or base counts)
# smry.min.base.qual - Ignores bases with base quality below x (not counted in coverage or base counts)


mpileup.params <- list(mp.disable.prob.realign = T, mp.mismatch.qual.downgrade = 50, 
  mp.min.map.qual = 0, mp.min.base.qual = 13, mp.gap.open.error.prob = 40, 
  mp.gap.extension.error.prob = 20, mp.homopolymer.error.coef = 100, mp.no.indels = F, 
  mp.indel.skip.depth = 5000
)
# mp.disable.prob.realign (-B) : Disable probabilistic realignment for the computation of base alignment quality (BAQ). 
#     BAQ is the Phred-scaled probability of a read base being misaligned. 
#     Applying this option greatly helps to reduce false SNPs caused by misalignments (default = T)
# mp.mismatch.qual.downgrade (-C) : Coefficient for downgrading mapping quality for reads containing excessive mismatches. 
#     Given a read with a phred-scaled probability q of being generated from the mapped position, 
#     the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables 
#     this functionality; if enabled, the recommended value for BWA is 50 (default = 0)
# mp.min.map.qual (-q) : Minimum mapping quality for an alignment to be used in the consensus (default = 0)
# mp.min.base.qual (-Q) : Minimum base quality for a base to be considered in the consensus (default = 13)
# mp.gap.open.error.prob (-o) : Phred-scaled gap open sequencing error probability.
#     Reducing INT leads to more indel calls. (default = 40)
# mp.gap.extension.error.prob (-e) : Phred-scaled gap extension sequencing error probability. 
#     Reducing leads to longer indels. (default = 20)
# mp.homopolymer.error.coef (-h) : Coefficient for modeling homopolymer errors. Given an l-long homopolymer run, 
#     the sequencing error of an indel of size s is modeled as INT * s/l. (default = 100)
# mp.no.indels (-I) : Do not perform INDEL calling (default = F)
# mp.indel.skip.depth (-L) : Skip INDEL calling if the average per-sample depth is above INT. (default = 250)


bcftools.params <- list(bcf.indel.snp.ratio = 0.15, bcf.variant.prob.thresh = 0.5, bcf.scaled.mut.rate = 0.001)
# bcf.indel.snp.ratio (-i) : Ratio of INDEL-to-SNP mutation rate [0.15]
# bcf.variant.prob.thresh (-p) :  A site is considered to be a variant if P(ref|D)<FLOAT [0.5]
# bcf.scaled.mut.rate (-t) : Scaled mutation rate for variant calling [0.001]


vcf2fq.params <- list(vfq.min.depth = 3, vfq.max.depth = 100000, vfq.min.RMS.map.qual = 10, vfq.indel.filtering.window = 5)
# vfq.minimum.depth (-d) : Minimum depth for fastq file
# vfq.maximum.depth (-D) : Maximum depth for fastq file
# vfq.min.RMS.map.qual (-Q) : Minimum root mean square mapping quality for fastq file (default = 10)
# vfq.indel.filtering.window (-l)


fasta.params <- list(fasta.min.cov = 2, fasta.min.freq = 0.6, fasta.min.freq.cov = 5)
# fasta.min.cov : Minimum coverage for base calling (sites with coverage below this are assigned N's).
# fasta.min.freq : Minimum frequency of either the reference or alternate base for calling. 
#      If both bases are below this frequency, an N is assigned.
# fasta.min.freq.cov : Minimum coverage above which fasta.min.freq is applied. Sites below this and >= than 
#      fasta.min.cov will only be called if all reads agree.


snp.params <- list(snp.min.site.call.qual = 4, snp.min.site.emit.qual = 30,
  snp.min.base.qual = 17, snp.min.map.qual = 20, snp.max.coverage = 250
)
# snp.min.site.call.qual - the minimum phred-scaled Qscore threshold to separate high
#     confidence from low confidence calls.  Only genotypes with confidence >= this threshold are 
#     emitted as called sites. A reasonable threshold is 30 for high-pass calling (this is the default). 
#     Note that the confidence (QUAL) values for multi-sample low-pass (e.g. 4x per sample) calling 
#     might be significantly smaller with the new EXACT model than with our older GRID_SEARCH model, 
#     as the latter tended to over-estimate the confidence; for low-pass calling we tend to use much 
#     smaller thresholds (e.g. 4).
# snp.min.site.emit.qual - the minimum phred-scaled Qscore threshold to emit low confidence calls. 
#     Genotypes with confidence >= this but less than the calling threshold are emitted but marked as filtered.
#     The default value is 30.
# snp.min.base.qual - specifies the minimum base quality required to consider a base for calling. Default is 17.
# snp.min.map.qual - specifies the minimum read mapping quality required to consider a read for calling. 
#     Only used in indel calling. Default is 20.
# snp.max.coverage - specifies the maximum coverage at a locus for any given sample. Downsampled reads are 
#     randomly selected from all possible reads at a locus. The default is 250 reads.


## Run pipeline ##
load("~/RCode/ngs.funcs.rdata")
params <- c(params, qa.params, align.params, realign.params, summary.params, 
  mpileup.params, bcftools.params, vcf2fq.params, fasta.params, snp.params
)
params$znum.list <- extract.ngs.znums(params$folder.name)
run.pipeline(params)