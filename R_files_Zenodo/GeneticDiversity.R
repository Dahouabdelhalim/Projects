library(adegenet)
library(pegas)
library(vcfR)
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ggdendro)
library(hierfstat)
library(geiger)
library(caper)
library(poppr)

snpgdsVCF2GDS("laredoensis_diversity.vcf", "test.gds", method="copy.num.of.ref", ignore.chr.prefix="RAD_")
genofile <- snpgdsOpen("test.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- scan("laredoensis_pops.txt", what=character())

flag <- pop_code %in% c("gularis", "laredoensisB", "laredoensisA", "sexlineatus")
samp.sel <- sample.id[flag]
pop.sel <- pop_code[flag]
Fst_results <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel), method="W&C84")
Fst_results$FstSNP
Fst_results$MeanFst

flag <- pop_code %in% c("gularis", "sexlineatus")
samp.sel <- sample.id[flag]
pop.sel <- pop_code[flag]
Fst_results <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel), method="W&C84")

snpgdsVCF2GDS("laredoensis_only.vcf", "test.gds", method="copy.num.of.ref", ignore.chr.prefix="RAD_")
genofile <- snpgdsOpen("test.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- scan("laredoensis_only_pops.txt", what=character())
flag <- pop_code %in% c("laredoensisB", "laredoensisA")
samp.sel <- sample.id[flag]
pop.sel <- pop_code[flag]
Fst_results <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel), method="W&C84")
write.table(Fst_results$FstSNP, file="laredoensis_Fst_allsnps.txt",row.names=FALSE)


data1 <- read.vcfR("laredoensis_diversity.vcf",nrows=22218)
my_genind <- vcfR2genind(data1)
pops <- c("gularis","gularis","gularis","gularis","gularis","gularis","gularis","laredoensisB","laredoensisB","laredoensisA","laredoensisA","laredoensisA","laredoensisA","laredoensisA","laredoensisA","laredoensisB","laredoensisB","laredoensisB","laredoensisB","laredoensisB","laredoensisA","sexlineatus","sexlineatus","sexlineatus","sexlineatus","sexlineatus","sexlineatus","sexlineatus")
my_genind@pop <- as.factor(pops)

basicstat <- basic.stats(my_genind, diploid = TRUE, digits = 2)
taxa <- colnames(basicstat$Ho)
df_total = data.frame()
for (i in 1:length(taxa)) {
  df <- c(mean(basicstat$Ho[,i],na.rm=TRUE),mean(basicstat$Hs[,i],na.rm=TRUE))
  df_total <- rbind(df_total,df)
}
colnames(df_total) <- c("Ho", "Hs")
row.names(df_total) <- taxa
df_total
write.table(df_total,file = "mean_summary.txt",sep = "\\t",quote = FALSE,row.names = TRUE, col.names = TRUE)