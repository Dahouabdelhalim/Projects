#install.packages("LEA-master", repos=NULL, type="source")
library(raster)

library("LEA")
vcf2lfmm("gm.north.6pops.maf0.05.recode.vcf", force=T) ###convert vcf to lfmm
cs_pca = pca("gm.north.6pops.maf0.05.recode.lfmm", scale=T) ###Create a pcaProject object: cs_pca.

tw = tracy.widom(cs_pca) # Perfom Tracy-Widom tests on all eigenvalues.

tw$pvalues[1:10]

pdf(file="gm.ddRAD.pca.explained.by.each.component.pdf", width=8, height=8)
plot(tw$percentage,xlim = c(1,15))  # plot the percentage of variance explained by each component
dev.off()

# main options, K: (the number of ancestral populations), entropy: calculate the cross-entropy criterion, CPU: the number of CPUs.
# Runs with K between 1 and 10 with cross-entropy and 10 repetitions.
project = NULL
project = snmf("gm.north.6pops.maf0.05.recode.lfmm", K=1:10, entropy = TRUE, repetitions = 10, project = "new", CPU = 50)

#pdf("SNMF.pdf")
#plot(project, col = "blue", pch = 19, cex = 1.2)
#dev.off()
#
Kn = 10 ### K 
ce <- matrix(nrow = 1, ncol = Kn)
for (i in 1:Kn){
  ce[,i] = mean(cross.entropy(project, K = i))
  }

## select the K with the lowest mean cross-entropy
best = which.min(ce)
best
#
#####if we look for loci that are extremely differentiated among populations, 
#####rather than loci that are extremely correlated with environmental gradients.
####you can run the following code:
#p = snmf.pvalues(project, entropy = TRUE, ploidy = 2, K = 6)
#pvalues = p$pvalues
#pdf("Differentiation_Outliers_K7.pdf",width = 16,height = 4)
#plot(-log10(pvalues), pch = 19, col = "blue", cex = .7, xlab = "SNP (ordered by contig arbitrarily)")
#dev.off()
#
##compute the 5% threshold
#P=-log10(pvalues)
#pod.thresh5=quantile(P,probs=0.95)
#pod.thresh5
##compute the 1% threshold
#pod.thresh1=quantile(P,probs=0.99)
#pod.thresh1
##compute the 0.1% threshold
#pod.thresh11=quantile(P,probs=0.999)
#pod.thresh11
##add the thresh to the actual XtX plot
#
#pdf(file="cs.snmf.outliers.pdf", width=16, height=4)
#plot(P,pch = 19, col = "blue", cex = .7, xlab = "SNP (ordered by contig arbitrarily)")
#abline(h=pod.thresh5,lty=2,col="red")
#abline(h=pod.thresh1,lty=2,col="red")
#abline(h=pod.thresh11,lty=2,col="red")
#dev.off()
##########get loci ID of outliers####
#library("vcfR")
#vcf <- vcfR::read.vcfR("gm.ddRAD.325.5-500X.0.999.33627snp.maf0.05.recode.vcf", verbose = FALSE) #prepare loci ID
#head(getID(vcf))
#head(getCHROM(vcf))
#head(getPOS(vcf))
##loci_ID <- getID(vcf)
#ID <- getID(vcf)
#CHROM <- getCHROM(vcf)
#POS <- getPOS(vcf)
#loci_ID <- cbind.data.frame(CHROM, POS, ID)
#cs.core <- cbind.data.frame(loci_ID, P)
#outlier.thresh <- pod.thresh11
##cs.bp.core.outliers <- cs.core[, "loci_ID"][which(cs.core[, "P"]> outlier.thresh)]  ###get loci ID of outliers
##write.table(cs.bp.core.outliers, file = "cs_ddrad_cs.bp.core_outliers.txt", quote = F, col.names = F, row.names = F)
#outliers.CHROM <- cs.core[, "CHROM"][which(cs.core[, "P"]> outlier.thresh)]
#outliers.POS <- cs.core[, "POS"][which(cs.core[, "P"]> outlier.thresh)]
#outliers.ID <- cs.core[, "ID"][which(cs.core[, "P"]> outlier.thresh)]
#write.table(cbind.data.frame(outliers.CHROM, outliers.POS, outliers.ID), file = "cs.snmf.outliers.0.001.txt", quote = F, col.names = T, row.names = F)
#
## plot cross-entropy criterion of all runs of the project
#
#pdf(file="cs.ddRAD.cross-entropy criterion of all runs of the project.pdf", width=8, height=8)
#plot(project, lwd = 5, col = "black", pch=1)
#dev.off()

project = NULL
project = lfmm("gm.north.6pops.maf0.05.recode.lfmm", "gm.north.env", K = best, ####using the best K
               repetitions = 4, project = "new", iterations = 10000, burnin = 5000, CPU = 60)

######if we have manys env. variables, we will use following commands to calculate P values
p.bio.adj <- matrix(nrow =47202, ncol =4 ) ####number of SNP, number of bios
colnames(p.bio.adj) <- c("bio04","bio05","bio08","bio18")
dim(p.bio.adj)
for (i in 1:4){
  z.bio.10rep = z.scores(project, K = 1, d = i)
  z.bio.mean <- apply(z.bio.10rep, 1, median)
  lambda.bios = median(z.bio.mean^2)/qchisq(0.5, df = 1)
  p.bio.adj[,i] = pchisq(z.bio.mean ^2/lambda.bios, df = 1, lower = FALSE)
}


q.bios <- matrix(nrow =47202, ncol =4) ####number of SNP, number of bios
dim(q.bios)
for (i in 1:4){
  q.bios[,i] <- qvalue(p.bio.adj[,i])$qvalues
}

##visually summarize large numbers of association tests is using a Manhattan plots,

pdf("LFMM_Manhattan.pdf",width = 16,height = 60)
par(mfrow = c(6,1))
for (i in 1:6){
  plot(-log10(q.bios[,i]), pch = 19, col = "blue", cex = .7, xlab = "SNP (ordered by contig arbitrarily)")
}
dev.off()


#######write the loci####
library("vcfR")
vcf <- vcfR::read.vcfR("gm.north.6pops.maf0.05.recode.vcf", verbose = FALSE) #prepare loci ID
head(getID(vcf))
head(getCHROM(vcf))
head(getPOS(vcf))
head(getFIX(vcf))
position <- getFIX(vcf)
write.table(position, file = "gm.snp.ID.txt", quote = F, col.names = T, row.names = F)

position <- read.table("gm.snp.ID.txt",header = T)
cs.bios <- cbind.data.frame(position, q.bios)
str(cs.bios)
head(cs.bios)
alpha <- 0.01

outliers.bio03 <- cs.bios[,c(1:5,8)][which(cs.bios[,8]< alpha),]
outliers.bio07 <- cs.bios[,c(1:5,9)][which(cs.bios[,9]< alpha),]
outliers.bio11 <- cs.bios[,c(1:5,10)][which(cs.bios[,10]< alpha),]
outliers.bio13 <- cs.bios[,c(1:5,11)][which(cs.bios[,11]< alpha),]
outliers.bio15 <- cs.bios[,c(1:5,12)][which(cs.bios[,12]< alpha),]

write.table(outliers.bio03, file = "outlies.bio04.0.01.txt", quote = F, col.names = T, row.names = F)
write.table(outliers.bio07, file = "outlies.bio05.0.01.txt", quote = F, col.names = T, row.names = F)
write.table(outliers.bio11, file = "outlies.bio08.0.01.txt", quote = F, col.names = T, row.names = F)
write.table(outliers.bio13, file = "outlies.bio18.0.01.txt", quote = F, col.names = T, row.names = F)
write.table(outliers.bio15, file = "outlies.bio13.0.01.txt", quote = F, col.names = T, row.names = F)

