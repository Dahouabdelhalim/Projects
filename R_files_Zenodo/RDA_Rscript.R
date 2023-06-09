setwd('')

library(vegan)
library(BBmisc)
library(caret)
library(ggplot2)
library(reshape2)
library(factoextra)
library(psych)
library(plyr)
library(qqman)
library(HDMD)

#-----------------------------------------------------------------------------------------------
## Select variables
#-----------------------------------------------------------------------------------------------
geometrics <- read.table('ClimCoords.txt', header = T, sep = '\\t')
geometrics <- geometrics[order(geometrics$Pop, geometrics$ID),]

## ORJU: Keeping species and climatic variables:
clim <- geometrics[ ,c("Pop", "bio_1", "bio_12", "bio_10", "bio_18", "bio_3", "bio_4", "bio_15", "tree_r", "ndvi_r3", "ndvistd_r", "srtm")]
colnames(clim)[1] <- 'Group'
clim$Group <- as.factor(clim$Group)
clim <- normalize(clim, method ='scale')
names(clim)
levels(clim$Group)

## Explore variability
# Corr panels
pairs.panels(clim[,2:12], scale=F)

# Boxplots
clim.melt <- melt(clim)

for (i in levels(clim.melt$variable)){
print(ggplot(clim.melt[clim.melt$variable==i,], aes(x = Group, y = value, fill = variable)) +
  geom_boxplot()
  )
}

# Distances
var.maha <- data.frame(row.names = colnames(clim[,2:12]))
h=0

for (i in colnames(clim[,2:12])){
	h <- h+1
	maha.envtest <- pairwise.mahalanobis(as.data.frame(clim[,i]), grouping = clim$Group, digits = 4 )
	D.mahalanobis.env <- sqrt(maha.envtest$distance)
	colnames(D.mahalanobis.env) <- row.names(maha.envtest$means)
	row.names(D.mahalanobis.env) <- row.names(maha.envtest$means)
	print(D.mahalanobis.env)
	var.maha$Value[h] <- sum(D.mahalanobis.env[3,])
}
var.maha$Variable <- row.names(var.maha)
summary(var.maha$Value)
row.names(var.maha[var.maha$Value>=mean(var.maha$Value),])

clim <- normalize(geometrics[ ,c("ID", "Pop", "bio_1", "bio_10", "bio_4", "ndvi_r3", "srtm")])

#-----------------------------------------------------------------------------------------------
## Redundacy Analysis
#-----------------------------------------------------------------------------------------------
# Load 012 dataset
snps.dataset <- read.table('ORJU06_biall_dp450_q40_hwe000001_KNOWNmac3_ppmiss05.012', header = F, row.names = 1, sep = '\\t', na.strings = '-1')
sum(is.na(snps.dataset))
snps.indv <- read.table('ORJU06_biall_dp450_q40_hwe000001_KNOWNmac3_ppmiss05.012.indv')
snps.dataset <- apply(snps.dataset, 2, function(x){ 
  x[is.na(x)] <- names(which.max(table(x)))
  return(x) })
sum(is.na(snps.dataset))

write.table(snps.dataset, 'ORJU06_biall_dp450_q40_hwe000001_KNOWNmac3_ppmiss05_Imp.012', sep = '\\t', col.names = F, quote = F, row.names = F)
snps.dataset <- read.table('ORJU06_biall_dp450_q40_hwe000001_KNOWNmac3_ppmiss05_Imp.012', header = F, sep = '\\t', stringsAsFactors = FALSE)

str(snps.dataset)
row.names(snps.dataset) <- snps.indv$V1
snps.dataset <- snps.dataset[order(row.names(snps.dataset)), ]
snps.dataset.patt <- snps.dataset

#-----------------------------------------------------------------------------------------------
## Environmental Data
#-----------------------------------------------------------------------------------------------
## Validating datset for RDA
clim.model <- clim
row.names(clim.model) <- clim.model$ID
clim.model <- na.omit(clim.model[row.names(snps.dataset),])
clim.model <- clim.model[order(row.names(clim.model)),]
identical(row.names(clim.model), row.names(snps.dataset.patt))
clim.set <- clim.model[,-(1:2)]

## RDA 
rda.clim.stepf <- rda(snps.dataset.patt ~., clim.set)
vif.cca(rda.clim.stepf)
rda.clim.stepf <- rda(snps.dataset.patt ~ bio_10 + bio_4 + ndvi_r3 + srtm, clim.set)

RsquareAdj(rda.clim.stepf)
summary(rda.clim.stepf)
climstepf.summary <- summary(rda.clim.stepf)
capture.output(climstepf.summary, file ='climstepf_RDA_summaryreport.txt')
permutation_sig.simpleRDA <- anova.cca(rda.clim.stepf, step = 1000)
permutation_axissig.simpleRDA <- anova.cca(rda.clim.stepf, by = 'axis', step = 1000)

# RDA1 vs RDA2
color.figure <- read.table('colors_figure.txt', header = T, row.names = 1, sep = '\\t')
color.figure <- color.figure[order(row.names(color.figure)),]
colvec <- rgb(red = color.figure, names = row.names(color.figure), maxColorValue = 255, alpha = 255)

plot(rda.clim.stepf, type="n", scaling=3)
points(rda.clim.stepf, display="sites",
       col = colvec[clim.model$Species], pch = 21, bg = colvec[clim.model$Pop],
       scaling = 3, cex = 1.3)
text(rda.clim.stepf, scaling=3, display="bp", col="black", cex=1)
legend("bottomright", legend = levels(clim.model$Pop), bty="n", col = colvec, pch = 21, pt.bg = colvec)

#-----------------------------------------------------------------------------------------------
## Identify candidate SNPs involved in local adaptation
#-----------------------------------------------------------------------------------------------
load.rda <- scores(rda.clim.stepf, choices=c(1:2), display="species")  # Species scores for the first three constrained axes

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
summary(load.rda[,1])
summary(load.rda[,2])
sd.RD1 <- sd(load.rda[,1])
sd.RD2 <- sd(load.rda[,2])

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

## Build dataset
cand1 <- outliers(load.rda[,1],6)
cand2 <- outliers(load.rda[,2],6)
ncand <- length(cand1) + length(cand2)

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))

colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2)
cand$snp <- as.character(cand$snp)

## Add correlations
foo <- matrix(nrow=(ncand), ncol=4)  # 4 columns for 4 predictors
colnames(foo) <- c("bio_10", "bio_4", "ndvi_r3", "srtm")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- snps.dataset.patt[,nam]
  foo[i,] <- apply(clim.set[,c("bio_10", "bio_4", "ndvi_r3", "srtm")],2,function(x) cor(x,snp.gen))
}

## Complete the table
cand <- cbind.data.frame(cand,foo)
cand <- cand[!duplicated(cand$snp),]

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,8] <- names(which.max(abs(bar[4:7]))) # gives the variable
  cand[i,9] <- max(abs(bar[4:7]))              # gives the correlation
}

colnames(cand)[8] <- "predictor"
colnames(cand)[9] <- "correlation"

table(cand$predictor)
write.table(cand, "RDA_CandidateSNPs_table.txt", col.names = T, row.names = F, sep = '\\t', quote = F)

#-----------------------------------------------------------------------------------------------
## vcf file treatment and regions to blast
#-----------------------------------------------------------------------------------------------
# Extract vcf header
lines <- readLines("ORJU06_biall_dp450_q40_hwe000001_KNOWNmac3_ppmiss05.recode.vcf") 
id <- grep("^#.*", lines) 
lines.head<-lines[id]

# Writing subsets
vcf.table <- read.table("ORJU06_biall_dp450_q40_hwe000001_KNOWNmac3_ppmiss05.recode.vcf", header = F)

snpssel.vector <- as.numeric(gsub( "V", "", as.character(cand$snp)))
vcf.table.selected <- vcf.table[snpssel.vector,]
write.table(vcf.table.selected,
            file = 'ORJU06_biall_dp450_q40_hwe000001_KNOWNmac3_ppmiss05_selENM.vcf',
            quote = F, sep = '\\t', row.names = F, col.names = F)

# To extract regions
snps.toblast <- vcf.table.selected[,1:2]
snps.toblast$V3 <- snps.toblast$V2-1000
snps.toblast$V4 <- snps.toblast$V2+1000
write.table(snps.toblast, 'snps_toblast_2000pb.txt', col.names = F, row.names = F, sep = '\\t', quote = F)

# Readding header
fConn <- file('ORJU06_biall_dp450_q40_hwe000001_KNOWNmac3_ppmiss05_selENM.vcf', 'r+') 
Lines <- readLines(fConn) 
writeLines(c(lines.head, Lines), con = fConn) 
close(fConn)
