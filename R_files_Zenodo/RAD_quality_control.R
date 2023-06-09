##A recent special issue in Molecular Ecology Resources provides a nice overview of the arsenal of tools available in R (Paradis et al., 2017). New tools have become available in R for analyzing HTS data including adegenet (Jombart, 2008), ape (Paradis, Claude & Strimmer, 2004), vcfR (Knaus & Gr??nwald, 2017), and poppr (Kamvar, Tabima & Gr??nwald, 2014; Kamvar, Brooks & Gr??nwald, 2015) among others. 
##pipe in R: ">", see https://rforcats.net/#pipes

##install.packages(c("poppr", "mmod", "magrittr", "treemap"), repos = "http://cran.rstudio.com", dependencies = TRUE)
library(vcfR)
library(pinfsc50)
library(ape)
library(reshape2)
#setwd("/home/ljcao/eb_rad/3_filter/miss0.2")
getwd()

########################Part 1: import and read vcf file into R, check data#################
eb.vcf.raw <- vcfR::read.vcfR("populations.snps.vcf", convertNA = TRUE,verbose = FALSE) #time-costing

#save.image("/home/ljcao/eb_rad/2_stack/3_callsnp/image.filter")

eb.vcf <- eb.vcf.raw
#strwrap(eb.vcf@meta[1:10])
#queryMETA(eb.vcf)
#queryMETA(eb.vcf, element = 'DP')
#queryMETA(eb.vcf, element = 'FORMAT=<ID=DP')
#head(getFIX(eb.vcf))
eb.vcf@gt[1:10, 1:10]

#########################################Part 2: Quality Control############################
##check summary data, including missing data
eb.vcf

##extract depth of genotype#
dp <- extract.gt(eb.vcf, element = "DP", as.numeric=TRUE) #time-costing

##barplot of missingness for every sample#
myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) }) ###MARGIN = 2, missingness for every sample.data from bwa and bowtie2, miss data will represent with NA, while data from bowtie, miss data will be indecated by 0.
#myMiss <- apply(dp, MARGIN = 2, function(x){ sum(x==0) })   ####data from bowtie?
library(RColorBrewer)
myMiss <- myMiss/nrow(eb.vcf)
pdf(file="eb_ddRAD_missing_barplot11.pdf", width=30, height=8)  ##histogram of missingness for every sample#
palette(brewer.pal(n=12, name = 'Set3'))
par(mar = c(16,4,4,2))
barplot(myMiss, las = 2, col = 1:16)
title(ylab = "Missingness (%)")
dev.off()

myMiss <- apply(dp, MARGIN = 1, function(x){ sum(is.na(x)) }) ###MARGIN = 1, missingness for snps. data from bwa, miss data will represent with NA, while data from bowtie, miss data will be indecated by 0.
#myMiss <- apply(dp, MARGIN = 1, function(x){ sum(x==0) })   ####data from bowtie
myMiss <- myMiss/ncol(eb.vcf@gt[,-1])
pdf(file="eb_ddRAD_missing_histogram.pdf", width=10, height=8)  ####histogram of missingness for snps#
hist(myMiss, col = "#8DD3C7", xlab = "Missingness (%)", main = "")
dev.off()

##draw barplot of depth for every sample#
pdf(file="eb_ddRAD_depth_boxplot.pdf", width=30, height=8)
par(mar=c(12,4,4,2))
boxplot(dp, col=2:8, las=3)
title(ylab = "Depth (DP)")
dev.off()

##draw violin plot of depth#
library(reshape2)
library(ggplot2) 
library(cowplot)

##Melt our matrix into a long form data.frame.
dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
#dpf <- dpf[ dpf$Depth > 5,]
#str(dpf)
##Create a row designator, You may want to adjust this
samps_per_row <- 25                               ####it is better to be divided evenly by numbel of individuals
myRows <- ceiling(length(levels(dpf$Sample))/samps_per_row)
myList <- vector(mode = "list", length = myRows)

for(i in 1:myRows){
  myIndex <- c(i*samps_per_row - samps_per_row + 1):c(i*samps_per_row)
  myIndex <- myIndex[myIndex <= length(levels(dpf$Sample))]
  myLevels <- levels(dpf$Sample)[myIndex]
  myRegex <- paste(myLevels, collapse = "$|^")
  myRegex <- paste("^", myRegex, "$", sep = "")
  myList[[i]] <- dpf[grep(myRegex, dpf$Sample),]
  myList[[i]]$Sample <- factor(myList[[i]]$Sample)
}

# Create the plot.
	pdf("eb_ddRAD_depth_violin.pdf", width = 10, height=8)
	myPlots <- vector(mode = "list", length = myRows)
	for(i in 1:myRows){
	  myPlots[[i]] <- ggplot(myList[[i]], aes(x=Sample, y=Depth)) + 
	    geom_violin(fill="#8dd3c7", adjust=1.0, scale = "count", trim=TRUE)
	  
	  myPlots[[i]] <- myPlots[[i]] + theme_bw()
	  myPlots[[i]] <- myPlots[[i]] + theme(axis.title.x = element_blank(), 
	                                       axis.text.x = element_text(angle = 60, hjust = 1))
	  myPlots[[i]] <- myPlots[[i]] + scale_y_continuous(trans=scales::log2_trans(), 
	                                                    breaks=c(1, 10, 100, 800),
	                                                    minor_breaks=c(1:10, 2:10*10, 2:8*100))
	  myPlots[[i]] <- myPlots[[i]] + theme( panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6) )
	  myPlots[[i]] <- myPlots[[i]] + theme( panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2) )
	}
	myPlots
	dev.off()

##heatmap of depth#
pdf("eb_ddRAD_depth_raw_variant_heatmap.pdf", width = 40, height=40)
heatmap.bp((dp[1:1000,])^(1/4),rlabels = F)           #########some sequenced samples are so deep that it was hard to read by heatmap, so the depth was divided by the square root---Lijun
dev.off()

##ommit data, replace low and high coverage SNPs (DP value) with NA#. MARGIN=2, missingness for every samples.
dp.raw <- extract.gt(eb.vcf.raw, element = "DP", as.numeric=TRUE)
dpf.raw <- melt(dp.raw, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
dpf.raw <- dpf.raw[which(dpf.raw$Depth<1000),]

pdf("Distribution of deepth.19pops.pdf",width = 16,height = 6)
ggplot(dpf.raw,aes(x=dpf.raw$Depth))+geom_histogram(binwidth =10,fill="lightblue",colour="lightblue")
dev.off()

#sum(rel.dep.snp>=quants[6]) ###number of SNPs higher than 4th quantile.

#snp.vcf  <- snp.vcf[rel.dep.snp<quants[6], ] #change value to keep final number of SNPS.
#snp.vcf.high.DP  <- snp.vcf[rel.dep.snp>=quants[4], ] 
dp <- extract.gt(eb.vcf, element = "DP", as.numeric=TRUE) #time-costing
#dp.high <- extract.gt(eb.vcf.high.DP, element = "DP", as.numeric=TRUE) #time-costing

dp[dp < 4] <- NA  #select SNPsat least 5X coverage, canbe adjusted.
dp[dp>500] <- NA

####ommited snps with higher missing data

eb.vcf@gt[,-1][ is.na(dp) == TRUE ] <- NA    #use -1 to omit the first column of the gt matrix,which is the ???FORMAT??? column
#eb.vcf@gt
#eb.vcf.raw@gt
myMiss <- apply(dp, MARGIN = 1, function(x){ sum( is.na(x) ) } ) 
myMiss <- myMiss / ncol(dp)
eb.vcf  <- eb.vcf[myMiss < 0.2, ] #change value to keep final number of SNPS.
eb.vcf
dp <- extract.gt(eb.vcf, element = "DP", as.numeric=TRUE)

####ommited individuals
myMiss <- apply(dp, MARGIN = 2, function(x){ sum( is.na(x) ) } ) 
myMiss <- myMiss / nrow(dp)
eb.vcf@gt<- eb.vcf@gt[, c(TRUE, myMiss < 0.2)] #myMiss is key to remove samples, determined according ps_ddRAD_depth_ommit_sample_barplot.tiff
eb.vcf
dp <- extract.gt(eb.vcf, element = "DP", as.numeric=TRUE)

pdf("eb_ddRAD_depth_ommit_sample0.1_heatmap1.pdf", width = 40, height=40)
heatmap.bp((dp[1:1000,])^(1/4), rlabels = FALSE)
dev.off()

write.vcf(eb.vcf, file = "eb.ddRAD.258.4-500X.28762snp.vcf.gz")
#write.vcf(eb.vcf.high.DP, file = "cs.ddRAD.218.high.DP.2706snp.vcf.gz")

system(command="vcftools --gzvcf eb.ddRAD.258.4-500X.28762snp.vcf.gz --mac 2 --recode --out eb.ddRAD.258.4-500X.87058snp.mac2")

eb.vcf.filter <- vcfR::read.vcfR("eb.ddRAD.258.4-500X.28577snp.mac2.recode.vcf", convertNA = TRUE,verbose = FALSE) #replace ./.:.:.:.:. with ./. before input
##position <- getFIX(cs.vcf.filter)
##write.table(position, file = "eb.ddRAD.215ind.106640snp.ID.txt", quote = F, col.names = T, row.names = F)
eb.vcf.filter
dp <- extract.gt(eb.vcf.filter, element = "DP", as.numeric=TRUE)
myMiss <- apply(dp, MARGIN = 1, function(x){ sum( is.na(x) ) } )
myMiss <- myMiss / ncol(dp)
myMiss[1:30]
eb.vcf.filter  <- eb.vcf.filter[myMiss < 0.2, ] #change value to keep final number of SNPS.
eb.vcf.filter

dp <- extract.gt(eb.vcf.filter, element = "DP", as.numeric=TRUE)
myMiss <- apply(dp, MARGIN = 2, function(x){ sum( is.na(x) ) } ) 
myMiss <- myMiss / nrow(dp)
eb.vcf.filter@gt<- eb.vcf.filter@gt[, c(TRUE, myMiss < 0.2)] #myMiss is key to remove samples, determined according ps_ddRAD_depth_ommit_sample_barplot.tiff
eb.vcf.filter

dp <- extract.gt(eb.vcf.filter, element = "DP", as.numeric=TRUE)
#
myMiss <- apply(dp, MARGIN = 2, function(x){ sum( is.na(x) ) } ) 
myMiss <- myMiss / nrow(dp)
write.table(myMiss, file = "eb.ddRAD.89738snp.miss.each.ind.more.strict.txt", quote = F, col.names = T, row.names = F)
sum(myMiss<0.05)          ###samples had < 5 % of missing data
max(myMiss)       ###the sample with the greatest number of missing loci 
myMiss[which(myMiss>0.05)]  ### the samples had > 5 % of missing data

sample.dp <- apply(dp,2, mean, na.rm=TRUE)  ###mean depth of every sample
eb.dp <- apply(dp,1, mean, na.rm=TRUE)  ###mean depth of every SNP
mean(sample.dp)
mean(eb.dp)
sum(eb.dp>30)
sum(sample.dp>20)

####we don't need following two steps if none SNP or sample was deleted in the second filter
write.vcf(eb.vcf.filter, file = "eb.ddRAD.258.4-500X.28571snp.mac2.recode.vcf.gz")
system(command="vcftools --gzvcf eb.ddRAD.258.4-500X.28571snp.mac2.recode.vcf.gz --maf 0.05 --recode --out eb.ddRAD.258.4-500X.28571snp.mac2")

dp.raw <- extract.gt(eb.vcf.raw, element = "DP", as.numeric=TRUE)
dpf.raw <- melt(dp.raw, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
dpf.raw <- dpf.raw[which(dpf.raw$Depth<1000),]

dp <- extract.gt(eb.vcf.filter, element = "DP", as.numeric=TRUE) 
dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)

pdf("Distribution.of.depth.befor.and.after.filter.pdf",width = 16,height = 6)
ggplot(dpf.raw,aes(x=dpf.raw$Depth))+geom_histogram(binwidth =10,fill="lightblue",colour="lightblue")
ggplot(dpf,aes(x=dpf$Depth))+geom_histogram(binwidth =10,fill="lightblue",colour="lightblue")
dev.off()

position <- getFIX(eb.vcf.filter)
write.table(position, file = "eb.snp.ID.txt", quote = F, col.names = T, row.names = F)

