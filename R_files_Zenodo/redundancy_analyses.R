##data preparation
#plink --vcf ps.all.203.obs0.8.maf0.05.19634.recode.vcf --make-bed --allow-extra-chr --out ps.203ind  ###for stacks 2.0
#plink --bfile ps.203ind --recode A --allow-extra-chr --out ps.203ind   ###for stacks 2.0
#install.packages(c("psych","vegan"), dependencies=TRUE)
#install.packages("adegenet")

library(psych)    
library(vegan) 
library(raster)
library(rgdal)

sample.coord <-read.table("cs.sample.location.txt", header=T, stringsAsFactors=F)
sample.coord

cs.snp<- read.table("cs.218ind.raw",header = T,row.names = 2)[,-(1:5)] #read into R

dim(cs.snp)    #check
sum(is.na(cs.snp)) #NAs in the matrix 
#write.table(cs.snp,"cs.snp.maf0.05.7235.txt",quote = F)
###replace NAs with the most common genotype at each SNP across all individuals
cs.snp.imp <- apply(cs.snp,2,function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(cs.snp.imp))

#rownames(cs.snp.imp)
clim.points <- read.table("cs.env", header = T)
cs.env <- cbind(sample.coord,data.frame(clim.points))
cs.env$IID <- as.character(cs.env$IID) # Make individual names characters (not factors)

coord <- cs.env[,3:4]

pca2 <- prcomp(cs.snp.imp, scale=T)

cs.env <- cbind(cs.env, pca2$x[,1:2])
str(cs.env)
# Confirm that genotypes and environmental data are in the same order, it must be TURE
identical(rownames(cs.snp.imp), cs.env[,2]) 

#remove strong correlated varables

pred <- cs.env[,5:ncol(cs.env)]
pred <- subset(pred, select = -c(bio01,bio03,bio04,bio05,bio07,bio08,bio09,bio11,bio12,bio16,bio17,bio18,bio19))
str(pred)

##full
##partial
cs.rda.p <- rda(cs.snp.imp ~ bio02+bio06+bio10+bio13+bio14+bio15+Condition(PC1+PC2), data = pred, scale = T)
cs.rda.p  ###show results
###The proportion of the variance explained by the environmental predictors is 
###given under the ???Proportion??? column for ???Constrained???; this is equivalent to 
###the R2 of a multiple regression. Just like in multiple regression, 
###this R2 will be biased and should be adjusted based on the number of predictors. 
###We can calculate the adjusted R2 using:

RsquareAdj(cs.rda.p)
####Our constrained ordination explains about 13.5% of the variation (OFM data)
###SNP loadings from the three significant constrained axes:
load.rda.p <- scores(cs.rda.p, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

### a function to identify SNPs that load in the tails of these distributions.
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

##apply it to each significant constrained axis:
cand1 <- outliers(load.rda.p[,1],4) # 32
cand2 <- outliers(load.rda.p[,2],4) # 5
cand3 <- outliers(load.rda.p[,3],4) # 75
####standard deviation cutoff can be 2.5, 3, or 3.5. higher sd identify loci under stronger selection.

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand
###For OFM data, we have 32 candidates on axis 1, 5 on axis 2, and 75 on axis 3, for a total of 112 candidate SNPs

###organize our results by making one data frame with the axis, SNP name, loading, & correlation with each predictor:
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

###add in the correlations of each candidate SNP with the eight environmental predictors:
ofm <- matrix(nrow=(ncand), ncol=8)  # 10 columns for 10 predictors
colnames(ofm) <- c("bio02","bio06","bio10","bio13","bio14","bio15","PC1","PC2")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- cs.snp.imp[,nam]
  ofm[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,ofm)  
head(cand)


###Now we have a data frame of 141 candidate SNPs and their correlation with our 8 environmental predictors.

###looking for duplicate detections. These are SNPs that are identified as candidates on more than one RDA axis.
length(cand$snp[duplicated(cand$snp)])

ofm <- cbind(cand$axis, duplicated(cand$snp)) 
table(ofm[ofm[,1]==1,2]) # duplicates on axis 1
table(ofm[ofm[,1]==2,2]) # duplicates on axis 2
table(ofm[ofm[,1]==3,2]) # duplicates on axis 3

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections, if any
head(cand)

###which of the predictors each candidate SNP is most strongly correlated with:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,12] <- names(which.max(abs(bar[4:11]))) # gives the variable, j in cand[i,j] is 3+biovars+PCNMvars+1, m and n in bar[m:n] is m=4, n= 4+vasrs(biovars + PCNMvars)
  cand[i,13] <- max(abs(bar[4:11]))              # gives the correlation,j in cand[i,j] is 3+biovars+PCNMvars+2, m and n in bar[m:n] is m=4, n= 4+vasrs(biovars + PCNMvars)
}

colnames(cand)[12] <- "predictor"
colnames(cand)[13] <- "correlation"

table(cand$predictor) 
write.csv(cand,"snp.candiates.sd4.full.csv")
