# Parentage analysis of POND 6 fish using clean snplist

# setwd("/Users/schluter/zoologyCloud/stk/g.r.e.a.t/dryad/")

x <- read.csv("genotypes all crosses 454 markers.csv", stringsAsFactors = FALSE)
x1 <- split(x, x$SNP.Name)

# -----------------------------------
# Parentage of F2's using MasterBayes

# extract the genotypes (use genotype)
# decide on a threshold probability for genotypes (use 0.99)
prob.min <- 0.99
z1 <- lapply(x1,function(x){
		x1 <- x[,c("id","genotype","cross","prob")]
		x1$genotype[x1$prob < prob.min] <- NA
		x2 <- x1[x$cross=="F1" | x$cross=="F2",]
		return(x2)
		})
		
names(z1) <- as.integer(as.factor(names(z1))) # simplify marker names
id <- z1[[1]]$id

# Split genotypes to get two alleles for each marker for each individual
z2 <- lapply(z1,function(x){
		z <- strsplit(x$genotype,split="")
		a <- sapply(z,function(z){return(z[1])})
		b <- sapply(z,function(z){return(z[2])})
		a[a=="-"] <- NA
		b[b=="-"] <- NA
		return(data.frame(a=a,b=b,stringsAsFactors=FALSE))
		})

# Put allele data into a genotype data frame for MasterBayes
z3 <- do.call("cbind.data.frame",z2)
fishG <- cbind.data.frame(id,z3)
names(fishG) # check allele names
fishG[1:15,1:10] # peek at data file
	         # id 1.a 1.b 2.a 2.b 3.a 3.b 4.a 4.b 5.a
	# 1  GEF1fem1   A   A   C   C   A   C   A   G   A
	# 2  GEF1fem3   A   G   C   C   A   C   A   G   A
	# 3  GEF1fem4   A   G   C   C   A   C   G   G   A
	# 4  GEF1fem6   A   A   C   C   A   C   A   G   A
	# 5  GEF1fem7   A   G   C   C   A   C   A   A   G
	# 6  GEF1fem8   A   G   C   C   A   C   A   G   A
	# 7  GEF1mal3   A   A   C   G   A   C   A   A   A
	# 8  GEF1mal4   A   A   C   G   A   C   A   G   G
	# 9  GEF1mal5   A   A   C   G   A   C   A   A   A
	# 10 GEF1mal8   A   A   C   G   A   C   A   A   A
	# 11   F2.011   A   A   C   C   A   C   A   A   A
	# 12   F2.013   A   G   C   C   C   C   A   G   A
	# 13   F2.014   A   A   C   C   A   C   A   G   A
	# 14   F2.015   A   G   C   C   A   A   A   A   A
	# 15   F2.016   A   G   C   C   A   C   A   G   A

# Create the "phenotype" data frame for MasterBayes
fishP <- data.frame(id=fishG$id)
fishP$offspring <- rep(1,nrow(fishP))
fishP$offspring[grep("F1",fishP$id)] <- 0
fishP$sex <- rep(NA,nrow(fishP))
fishP$sex[grep("fem",fishP$id)] <- "Female"
fishP$sex[grep("mal",fishP$id)] <- "Male"
fishP$sex <- as.factor(fishP$sex)
fishP[1:15,]
	         # id offspring    sex
	# 1  GEF1fem1         0 Female
	# 2  GEF1fem3         0 Female
	# 3  GEF1fem4         0 Female
	# 4  GEF1fem6         0 Female
	# 5  GEF1fem7         0 Female
	# 6  GEF1fem8         0 Female
	# 7  GEF1mal3         0   Male
	# 8  GEF1mal4         0   Male
	# 9  GEF1mal5         0   Male
	# 10 GEF1mal8         0   Male
	# 11   F2.011         1   <NA>
	# 12   F2.013         1   <NA>
	# 13   F2.014         1   <NA>
	# 14   F2.015         1   <NA>
	# 15   F2.016         1   <NA>

library(MasterBayes) # uses the obsolete 'genetics' package
res1fish <- expression(varPed(x="offspring",restrict=0))
PdPfish <- PdataPed(formula=list(res1fish), data=fishP) # offspring, sex and id are implicit
# PdPfish <- PdataPed(formula=list(res1fish), data=fishP, USsire=TRUE, USdam=TRUE) # try!


# Create the genotype data object
# perlocus=FALSE means that different error rates for different loci are NOT modeled
#   (changing to TRUE didn't change anything)
# marker.type="MSW" refers to Wang's model of genotyping error for codominant markers, whereas
# marker.type="MSC" refers to CERVUS's model of genotyping error (this gave nonsensical results)
# 	(="AFLP" is for dominant markers)
GdPfish <- GdataPed(fishG, perlocus=FALSE, marker.type="MSW") 

# Set up the design matrix for model
# mm.tol is the number of genotype mismatches tolerated for potential parents
# 	Program gave me an error, no possible dams, if mm.tol set too low
gX <- getXlist(PdPfish, GdPfish, E1=0.02, E2=0.02, mm.tol=40)

# Estimate the pedigree
F2pedigree.ml <- MLE.ped(gX)
# F2pedigree.ml <- MLE.ped(gX, USdam=TRUE, USsire=TRUE, nUSdam=2, nUSsire=2) # try to add missing parents!

table(F2pedigree.ml$P[,2],F2pedigree.ml$P[,3]) # the most likely parent combinations
           GEF1mal3 GEF1mal4 GEF1mal5 GEF1mal8


x <- cbind.data.frame(F2pedigree.ml$P[grep("F2",fishP$id),],F2pedigree.ml$prob[grep("F2",fishP$id)],
		stringsAsFactors=FALSE)
names(x) <- c("F2id","F1fem","F1mal","prob")
hist(x$prob) # histogram of the probabilities of the most likely parent combination

# The most uncertain parentage calls:
# (the following two fish 2 were indicated in my homemade analysis as matching more than 1 mother,
#  and having a sample size of just 3 and 1 informative markers for assessing maternity, respectively)
# (we could try using a less stringent genotyping call and tolerate a higher error rate to see if
#  we can increase the amount of information for these individuals)
x[x$prob<0.99,]
	#      F2id    F1fem    F1mal      prob
	# 71 F2.162 GEF1fem8 GEF1mal5 0.5129392

# Set to NA the parents of the two fish that were not assigned to the correct F1's
x$F1fem[is.element(x$F2id,c("F2.162","F2.336"))] <- NA
x$F1mal[is.element(x$F2id,c("F2.162","F2.336"))] <- NA
x$prob[is.element(x$F2id,c("F2.162","F2.336"))] <- 0
nrow(x)
[1] 232

f2pedigree <- x

# Check on marker compatibility assuming this pedigree is correct:
# Determine if the F2 genotype is one of the possible genotypes from the F1xF1 parents
# Generates an array of TRUE or FALSEs, indicating compatibility of F2 with parents for each marker.
# Remember that by this stage F2 genotypes might be NA but none of the F1 parents will be
z1 <- lapply(x1,function(x){
	x$genotype[x$prob < prob.min] <- NA
	f2geno <- x$genotype[x$cross=="F2"]
	f2names <- x$id[x$cross=="F2"]
	femgeno <- x$genotype[match(f2pedigree$F1fem,x$id)]
	malgeno <- x$genotype[match(f2pedigree$F1mal,x$id)]
	z1 <- strsplit(femgeno,split="")
	z2 <- strsplit(malgeno,split="")
	result <- vector()
	for(i in 1:length(f2geno)){
	  if(!is.na(f2geno[i])){
		possible <- unique(c(paste(sort(c(z1[[i]][1],z2[[i]][1])),collapse=""),
					paste(sort(c(z1[[i]][1],z2[[i]][2])),collapse=""),
					paste(sort(c(z1[[i]][2],z2[[i]][1])),collapse=""),
					paste(sort(c(z1[[i]][2],z2[[i]][2])),collapse="")))
		if(all(nchar(possible) > 1))
				result[i] <- is.element(f2geno[i],possible)
			else result[i] <- NA
		}  else result[i] <- NA
	  }
	return(result)	
	})
z1 <- do.call("cbind.data.frame",z1)
rownames(z1) <- f2pedigree$F2id

# Here's what the array looks like:
z1[1:20,1:3]


# Which markers were most often in contradiction of parentage?
z2 <- apply(z1,2,function(x){
	sum(as.numeric(x),na.rm=TRUE)/length(x[!is.na(x)])
	})
mean(z2)

hist(z2,nclass=20)
cbind((sort(z2)[1:20]))


# Which markers had the most missing data?
z2 <- apply(z1,2,function(x){
	length(x[is.na(x)])
	})
median(z2)

hist(z2,nclass=20)
cbind((sort(z2, decreasing=TRUE)[1:20]))


# Which F2 genotypes conflicted most with their assigned F1 parents??
# Fraction of non-missing genotypes that were FALSE
z3 <- apply(z1,1,function(x){
	sum(as.numeric(x),na.rm=TRUE)/length(x[!is.na(x)])
	})
hist(z3,nclass=20)
z4 <- cbind.data.frame(f2pedigree[,2:3],fracTRUE=z3,stringsAsFactors=FALSE)
z4[order(z4$fracTRUE),][1:20,]


# Total number of FALSE calls:
z3 <- apply(z1,1,function(x){
	sum(as.numeric(!x),na.rm=TRUE)
	})
hist(z3,nclass=20)
z4 <- cbind.data.frame(f2pedigree[,2:3],nFALSE=z3,stringsAsFactors=FALSE)
z4[order(z4$nFALSE,decreasing=TRUE),][1:20,]


# -------------------------------------------------------------
# Parentage of F3's using MasterBayes
# This analysis includes the optional cleanup steps at the beginning
# No F2's have been dropped from this analysis
length(x1)
	# [1] 454

# Set up the data files.
# grab the genotypes (use genotype)
z1 <- lapply(x1,function(x){
		x1 <- x[,c("id","genotype","cross")]
		x1 <- x[,c("id","genotype","cross","prob")]
		x1$genotype[x1$prob<0.99] <- NA
		x2 <- x1[x$cross=="F2" | x$cross=="F3",]
		return(x2)
		})
names(z1) <- as.integer(as.factor(names(z1))) # simplify marker names
id <- z1[[1]]$id

# Split genotypes to get two alleles for each marker for each individual
z2 <- lapply(z1,function(x){
		z <- strsplit(x$genotype,split="")
		a <- sapply(z,function(z){return(z[1])})
		b <- sapply(z,function(z){return(z[2])})
		a[a=="-"] <- NA
		b[b=="-"] <- NA
		return(data.frame(a=a,b=b,stringsAsFactors=FALSE))
		})
		
z3 <- do.call("cbind.data.frame",z2)
fishG <- cbind.data.frame(id,z3)
names(fishG) # check allele names
fishG[1:15,1:10] # peek at data file
	       # id 1.a 1.b 2.a 2.b 3.a 3.b 4.a 4.b 5.a
	# 1  F2.011   A   A   C   C   A   C   A   A   A
	# 2  F2.013   A   G   C   C   C   C   A   G   A
	# 3  F2.014   A   A   C   C   A   C   A   G   A
	# 4  F2.015   A   G   C   C   A   A   A   A   A
	# 5  F2.016   A   G   C   C   A   C   A   G   A
	# 6  F2.017   A   G   C   C   A   C   A   A   A
	# 7  F2.024   A   G   C   C   A   C   A   G   A
	# 8  F2.025   A   A   C   C   A   A   A   A   A
	# 9  F2.027   A   A   C   C   A   C   A   G   A
	# 10 F2.029   A   A   C   C   A   A   A   G   A
	# 11 F2.032   A   A   C   C   A   C   A   G   A
	# 12 F2.035   A   A   C   C   C   C   A   G   G
	# 13 F2.036   A   A   C   C   A   A   A   G   A
	# 14 F2.041   A   A   C   C   A   C   A   A   A
	# 15 F2.051   A   A   C   C   A   A   A   G   A

fishP <- data.frame(id=fishG$id)
fishP$offspring <- rep(1,nrow(fishP))
fishP$offspring[grep("F2",fishP$id)] <- 0
fishP$sex <- rep(NA,nrow(fishP))
fishP$sex[grep("F2",fishP$id)] <- "Female"
fishP$sex <- as.factor(fishP$sex)
fishP[1:15,]
	       # id offspring    sex
	# 1  F2.011         0 Female
	# 2  F2.013         0 Female
	# 3  F2.014         0 Female
	# 4  F2.015         0 Female
	# 5  F2.016         0 Female
	# 6  F2.017         0 Female
	# 7  F2.024         0 Female
	# 8  F2.025         0 Female
	# 9  F2.027         0 Female
	# 10 F2.029         0 Female
	# 11 F2.032         0 Female
	# 12 F2.035         0 Female
	# 13 F2.036         0 Female
	# 14 F2.041         0 Female
	# 15 F2.051         0 Female

library(MasterBayes)
res1fish <- expression(varPed(x="offspring",restrict=0))
# Unsampled females and males noted
PdPfish <- PdataPed(formula=list(res1fish), data=fishP, USdam=TRUE, USsire=TRUE) 

# Create the genotype data object
# perlocus=FALSE means that different error rates for different loci are NOT modeled
# marker.type="MSW" refers to Wang's model of genotyping error for codominant markers, whereas
# marker.type="MSC" refers to CERVUS's model of genotyping error (this gave nonsensical results)
# 	(="AFLP" is for dominant markers)
GdPfish <- GdataPed(fishG, perlocus=FALSE, marker.type="MSW") 

# Set up the design matrix for model
# mm.tol is the number of genotype mismatches tolerated for potential parents
# 	Program gave me an error, no possible dams, if mm.tol set too low (????)
gX <- getXlist(PdPfish, GdPfish, E1=0.01, E2=0.01, mm.tol=20)

# Find the MLE pedigree
F3pedigree.ml <- MLE.ped(gX, USdam=TRUE, USsire=TRUE, nUSdam=10, nUSsire=200)

x <- cbind.data.frame(F3pedigree.ml$P[grep("F3",fishP$id),1:2],F3pedigree.ml$prob[grep("F3",fishP$id)],
		stringsAsFactors=FALSE)
names(x) <- c("F3id","F2fem","prob")
hist(x$prob) # histogram of the probabilities of the most likely parent combination

# Distribution of numbers of F2 offspring per F2 female (not including 0's)
barplot(sort(table(x$F2fem),decreasing=TRUE),col="red",cex.names=0.5,ylab="No. F3 offspring")

f3pedigree <- x

# ----
# Calculate the number of F3 offspring per F2 female
x <- f3pedigree[f3pedigree$prob >= 0.6,]  # At least a 60% chance of maternity (can change this)
z1 <- tapply(x$prob,x$F2fem,length)

z <- as.character(fishP$id[grep("F2.",fishP$id, fixed=TRUE)])
z2 <- z1[match(z,names(z1))]
z2[is.na(z2)] <- 0
hist(z2)
z3 <- cbind.data.frame(z,z2,stringsAsFactors=FALSE)
names(z3) <- c("id","noffspring.p0.6")
f2fitness <- z3
