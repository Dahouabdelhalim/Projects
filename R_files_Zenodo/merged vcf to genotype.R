#converts .merged.vcf output file from SNP pipeline into a genotypes .csv file that can be put into a database and queried for genotypes by locus for validation and/or genotyping. 

merged.vcf.to.genotype <- function(merged.vcf.fname, params) {
  options(stringsAsFactors = F)	 
  vcf <- scan(merged.vcf.fname, what = "character", sep = "\\n", quiet = T)
  first.line <- grep("#CHROM\\t", vcf)
  if(is.null(first.line)) return(NULL)
  vcf <- scan(merged.vcf.fname, what = "character", sep = "\\n", skip = first.line - 1, quiet = T)
  headers <- gsub("#CHROM", "CHROM", vcf[1])
  headers <- sapply(headers, function(x) strsplit(x, "\\t")[[1]])
  vcf <- vcf[-1]
  vcf <- data.frame(t(sapply(vcf, function(x) strsplit(x, "\\t")[[1]])))
  if(ncol(vcf) == 0) return(NULL)
  colnames(vcf) <- headers  
  
  vcf$POS <- as.numeric(vcf$POS)
  vcf$QUAL <- as.numeric(vcf$QUAL)  
  
  parse.format <- function(x, fmt) {
    x.list <- as.list(strsplit(x, ":")[[1]])
    if(x.list[[1]] == "./.") return(c(allele1.count = 0, allele2.count = 0, total.reads = 0, gtype.qual = 0))
    names(x.list) <- strsplit(fmt, ":")[[1]]
    allele.counts <- as.numeric(strsplit(x.list$AD, ",")[[1]])
    c(allele1.count = allele.counts[1], allele2.count = allele.counts[2], 
      total.reads = as.numeric(x.list$DP), gtype.qual = as.numeric(x.list$GQ))
  }
  
  extract.znum <- function(x) {
    znum <- strsplit(x, "_")[[1]][1]   
    as.numeric(gsub("z", "", znum))
  }
  
  result <- do.call(rbind, lapply(10:ncol(vcf), function(col) {
    labid <- extract.znum(colnames(vcf)[col])
    gtype.pos <- do.call(rbind, lapply(1:nrow(vcf), function(row) {
      ref <- vcf[row, "REF"]
      alt <- vcf[row, "ALT"]
      x <- parse.format(vcf[row, col], vcf[row, "FORMAT"])
      sum.ref.alt <- sum(x[c("allele1.count", "allele2.count")])
      pct.ref <- x["allele1.count"] / sum.ref.alt
      alleles <- if(sum.ref.alt >= params$min.reads.het) {
        if(pct.ref > params$pct.ref.hom) {
          c(ref, ref)
        } else if(pct.ref < (1 - params$pct.ref.hom)) {
          c(alt, alt)
        } else if(pct.ref >= params$pct.ref.het.min & pct.ref <= params$pct.ref.het.max) {
          c(ref, alt)
        } else c(NA, NA)
      } else if(sum.ref.alt >= params$min.reads.hom) {
        if(pct.ref == 1) {
          c(ref, ref)
        } else if(pct.ref == 0) {
          c(alt, alt)
        } else c(NA, NA)
      } else c(NA, NA)
      data.frame(locus = vcf[row, "CHROM"], position = vcf[row, "POS"], ref = ref, alt = alt, 
                 gtype.qual = x["gtype.qual"], labid = labid, allele1 = alleles[1], allele2 = alleles[2], 
                 allele1.count = x["allele1.count"], allele2.count = x["allele2.count"],
                 total.reads = x["total.reads"])
    }))
  }))
  rownames(result) <- NULL
  result
}


params <- list(min.reads.het = 7, min.reads.hom = 5, pct.ref.hom = 0.8, pct.ref.het.min = 0.3, pct.ref.het.max = 0.7)

# min.reads.het = minimum # of reads to call a heterozygote
# min.reads.hom = minimum # of reads to call a homozygotes
# pct.ref.het.min = minimum % ref allele to call a heterozygote
# pct.ref.het.max = maximum % ref allele to call a heterozygote


merged.vcf.fname <- "Oorc_Transient-snps_20120911_1028.merged.vcf"
result <- merged.vcf.to.genotype(merged.vcf.fname, params)
write.csv(result, paste(merged.vcf.fname, ".genotypes.csv", sep = ""), row.names = F)
