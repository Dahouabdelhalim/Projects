# Genotype changes between F2's and F3's

# Run the "rqtl analysis.r" script file in R before this point because we use the rqtl objects again here

# ------
# Genotypes of F2 females at markers near Eda

# fitness peak qtl
z <- pull.geno(rqtlAll.gp, 4)[,"chrIV:12815024"] 
z1 <- table(z) # 1=AA, 2=AB, 3=BB
z1
	 # 1  2  3 
	# 55 99 54

# allele frequency (low allele)
(2 * z1[3] + z1[2]) / (2 * sum(z1))
# 0.4975962 

# ------
# Genotypes based on genotype probabilities at peak marker for fitness

x <- pull.genoprob(rqtlAll.gp, chr = 4)
z <- x[, grep("chrIV:12815024", colnames(x))]
head(z)
	            # chrIV:12815024:AA chrIV:12815024:AB      chrIV:12815024:BB
	# F2.027 0.00000000378746991632 0.999999992424943 0.00000000378758731936
	# F2.029 0.00000000378746993790 0.999999992425060 0.00000000378746993790
	# F2.032 0.00000000378746992902 0.999999992424989 0.00000000378754037196
	# F2.051 0.00000000000007174365 0.000000003788052 0.99999999621187551924
	# F2.063 0.99999999621187640741 0.000000003788052 0.00000000000007174365
	# F2.071 0.00000000378746993898 0.999999992425061 0.00000000378746993790

z1 <- apply(z, 1, max)
range(z1)
	# [1] 0.8028931 1.0000000
z1 <- table( apply(z, 1, function(x){which(x == max(x))}) )
z1
	  # AA  AB  BB # B is the freshwater (low) allele
	  # 1   2   3 
	 # 60 102  62 
expect <- c(.25,.5,.25)*sum(z1)
expect
	# [1]  56 112  56 # expected frequencies

chisq.test(z1, p = c(.25,.5,.25))
	# data:  z1
	# X-squared = 1.8214, df = 2, p-value = 0.4022

# allele frequency (low allele)
(2 * z1[3] + z1[2]) / (2 * sum(z1))
	# 0.5044643 

# ==========================================================================
# F3 genotype frequencies at fitness peak marker

f3genotypes <- read.csv("genotypes of f3s.csv", row.names = 1, stringsAsFactors = FALSE)

z1 <- table(unlist(f3genotypes[rownames(f3genotypes) == "chrIV:12815024", ]), useNA = "always")
z1
  # AA   AC   CC <NA> 
  # 66  241  139   28 
  
  
# Grab F3 markers that are on the linkage map

linkmap <- pull.map(rqtlAll.gp, chr = 4, as.table = TRUE)
linkmap[rownames(linkmap) == "chrIV:12815024", ]
	               # chr    pos
	# chrIV:12815024   4 49.685

x <- f3genotypes[match(rownames(linkmap), rownames(f3genotypes)), ]
nrow(x)
	# [1] 39

# -----
# Generate a consensus marker state for missing genotypes using genotypes in the vicinity of the fitness marker peak

# grab F3 genotypes for markers within 3 cM of the peak marker for fitness
x1 <- x[which(linkmap$pos > linkmap$pos[rownames(linkmap) == "chrIV:12815024"] - 3 & 
				linkmap$pos < linkmap$pos[rownames(linkmap) == "chrIV:12815024"] + 3), ]

# set F3 genotypes to 1, 2, or 3 depending on whether they are AA, AB, or BB
for(i in 1:nrow(x1)){
	# i <- 4
	z <- unlist(f0genotypes[grep(rownames(x1)[i], rownames(f0genotypes)), ])
	z1 <- unlist(strsplit(z[1], split = ""))
	z2 <- unlist(strsplit(z[2], split = ""))
	if(length(unique(z1)) > 1 | length(unique(z2)) > 1 | length(intersect(z1,z2)) > 0) x1[i, ] <- NA else{
		x1[i, x1[i, ] == z[1] & !is.na(x1[i, ])] <- "1"
		x1[i, x1[i, ] == z[2] & !is.na(x1[i, ])] <- "3"
		x1[i, (x1[i, ] == paste0(z1[1], z2[1]) | x1[i, ] == paste0(z2[1], z1[1])) & !is.na(x1[i, ])] <- "2"
		}
	}
x2 <- x1[!apply(x1, 1, function(x){all(is.na(x))}), ]
x2
	               # F3.001 F3.002 F3.003 F3.004 F3.005 F3.006 F3.007 F3.008 F3.009 F3.010 F3.011 F3.012
	# chrIV:10960835      2      2   <NA>      2      2      3      2      3      2   <NA>      3      2
	# chrIV:12811933      2      2      2      2      2      3      2      3      2   <NA>   <NA>      2
	# chrIV:12815024      2      2      2      2      2      3      2   <NA>      2      2      2      2
	# chrIV:15052901   <NA>   <NA>   <NA>   <NA>   <NA>      3   <NA>   <NA>   <NA>      2   <NA>   <NA>
	               # F3.013 F3.014 F3.015 F3.016 F3.017 F3.018 F3.019 F3.020 F3.022 F3.023 F3.024 F3.025
	# chrIV:10960835      1      3      2      1      2      1      1      3      2      2      2      2
	# chrIV:12811933      1      3      2      1      2      1   <NA>      3      2      2      2      2
	# chrIV:12815024      1      3      2      1      2      1   <NA>   <NA>      2      2      2      2
	# chrIV:15052901      2      3   <NA>   <NA>      2      2      2      3      2   <NA>   <NA>   <NA>
	# ...

# Generate the missing consensus marker (return NA if there is no consensus)
peakmarker <- rownames(x2) == "chrIV:12815024"
for(j in 1:ncol(x2)){
	# j <- 19
	if( is.na(x2[peakmarker, j]) ){
		z <- table(x2[!peakmarker, j])
		z1 <- names(z)[z==max(z)]
		if(length(z1) == 1) x2[peakmarker, j] <- z1
		}
	}

z2 <- table(unlist(x2[peakmarker, ]))
z2
  # 1   2   3 
 # 67 260 146 

# allele frequency (low allele)
pB <- unname( (2 * z2[3] + z2[2]) / (2 * sum(z2)) )
pB
# 0.5835095 

# HW expectation:
expect <- c(pB^2, 2*pB*(1-pB), (1-pB)^2)
chisq.test(z2, p = expect)	# Chi-squared test for given probabilities
	# X-squared = 108.71, df = 2, p-value < 2.2e-16

# -------
# SUMMARY

# Actual F2 genotypes at peak marker for fitness
	# 1=AA, 2=AB, 3=BB
	 # 1  2  3 
	# 55 99 54
# allele frequency (B: low allele)
# 0.4975962 

# F2 genotypes based on genotype probabilities
	  # 1   2   3 
	 # 60 102  62 
# allele frequency (B: low allele)
# 0.5044643 
sqrt(0.5044643*(1-0.5044643)/(60+102+62))
# [1] 0.001115982 SE

# Actual F3 genotypes at peak marker for fitness
  # AA   AC   CC <NA> 
  # 66  241  139   28 
# allele frequency (low allele is CC)
# 0.5474684

# F3 genotypes using consensus genotypes for NA values
  # 1   2   3 
 # 67 260 146 
# allele frequency (low allele)
# 0.5818386 
sqrt(0.5818386*(1-0.5818386)/(67+260+146))
# [1] 0.0005143815 SE

# ---
# Barplot with error bars
# graph the frequencies of the F2's and F3's
# SE's for a finite population of F2's are calculated as
#	sqrt( p(1-p)/n * (N-n)/(N-1) )
# Where n is sample size and N is known population size
# 636 F2's originally introduced, so 636/2 = 318 females (assume)
N1 <- 318
p1 <- c(60,102,62) # F2
n1 <- sum(p1)
	# [1] 224
p2 <- c(67,260,146) # F3
n2 <- sum(p2)
	# [1] 473
p1 <- p1/n1
p2 <- p2/n2
se1 <- sqrt( (p1*(1-p1)/n1) * ((N1-n1)/(N1-1)) )
se2 <- sqrt( p2*(1-p2)/n2)
z <- barplot(matrix(c(p1,p2), nrow = 2, byrow = TRUE), beside = TRUE, 
	col = c("white","firebrick"), space = c(0,0.6), cex.names = 1.2, 
	names.arg=c("MM","MF","FF"), xlab = "Genotype", ylab = "Relative frequency",
	las = 1, bty = "l", ylim = c(0,0.6))
# add error bars
segments(x0 = z[1,], y0 = p1-se1, y1 = p1+se1)
segments(x0 = z[2,], y0 = p2-se2, y1 = p2+se2)



# -------
# Compare allele frequency change with that predicted from differences in reproductive success of F2 females

library(emmeans)
library(visreg)

# From Hoekstra & Linnen (2010, Cold Spring Harbor):
	# If p and q are the frequencies of alleles A and B; 
	# wAA, wAB, and wBB are the relative fitnesses of genotypes AA, AB, and BB; 
	# andw is the mean fitness of the population, 
	# then change in the frequency of the A allele (Δp = p′ – p), is given by
	# Δp = pq[p(wAA -  wAB) + q(wAB - wBB)] /w 

deltaP <- function(p, freq, fitness){
	# p is the frequency of the A allele
	# freq is a vector of frequencies of the 3 genotypes
	# fitness is a vector if average fitnesses of the 3 genotypes 
	q <- 1 - p
	wAA <- fitness[1]
	wAB <- fitness[2]
	wBB <- fitness[3]
	wbar <- sum(freq * fitness) / sum(freq)
	dP <- p*q*( p*(wAA - wAB) + q*(wAB - wBB) ) / wbar
	names(dP) <- "change in allele A"
	return(dP)
	}

# Mean number of offspring for **KNOWN** genotypes
x <- factor(pull.geno(rqtlAll.gp, chr=4)[, "chrIV:12815024"])
freq <- table(x)
pA = (2 * freq[1] + freq[2]) / (2 * sum(freq))
fitness <- tapply(rqtlAll$pheno$noffspring.p0.6, x, mean, na.rm = TRUE)
fitness
#       AA       AB       BB
       # 1        2        3 
# 1.618182 1.777778 3.185185

# or take family into account
y <- rqtlAll$pheno$noffspring.p0.6
z <- lm(y ~ family + x)
emmeans(z, "x")
 # x   emmean        SE  df  lower.CL upper.CL
 # 1 1.566820 0.3080872 200 0.9593036 2.174336
 # 2 1.687694 0.2302505 200 1.2336637 2.141724
 # 3 3.083757 0.3121458 200 2.4682378 3.699276
visreg(z, xvar = "x")
visreg(z, xvar = "x", scale = "response", ylim = range(y), rug=FALSE)

deltaP(pA, freq, fitness)
	# change in allele A 
	       # -0.09287225

# Predicted frequency of allele B in the next generation:
unname( (1 - pA) - deltaP(pA, freq, fitness) )
# [1] 0.5904684


# Mean number of offspring for **INFERRED** genotypes from genotype probabilities
x <- pull.genoprob(rqtlAll.gp, chr = 4)
x1 <- x[, grep("chrIV:12815024", colnames(x))]
z1 <- apply(x1, 1, function(x){which(x == max(x))})
z1 <- recode(z1, c("1", "2", "3"), c("MM", "MF", "FF"))
z1 <- factor(z1, levels = c("MM", "MF", "FF"))
freq <- table(z1)
pA = (2 * freq[1] + freq[2]) / (2 * sum(freq))
fitness <- tapply(rqtlAll$pheno$noffspring.p0.6, z1, mean)
fitness
	      # MM       MF       FF 
	# 1.600000 1.754902 3.209677 
deltaP(pA, freq, fitness)
	# change in allele A 
	       # -0.09576438 
# Predicted frequency of allele B (low allele) in the next generation:
unname( (1 - pA) - deltaP(pA, freq, fitness) )
# [1] 0.6002287

# or take family into account
y <- rqtlAll$pheno$noffspring.p0.6
z <- lm(y ~ family + z1)
emmeans(z, "z1")
	 # z1   emmean        SE  df  lower.CL upper.CL
	 # MM 1.499391 0.2891321 216 0.9295098 2.069273
	 # MF 1.668428 0.2228812 216 1.2291271 2.107728
	 # FF 3.099231 0.2852085 216 2.5370829 3.661379
visreg(z, xvar = "z1")
visreg(z, xvar = "z1", scale = "response", ylim = range(y), rug=FALSE)
# Trick to superimpose the points
visreg(z, xvar = "z1", scale = "response", rug = FALSE, ylim = range(y), bty = "l", 
	xlab = "F2 female genotype at peak marker", ylab = "Estimated number of offspring")
	points(y ~ I(0.145 + (as.integer(z1) - 1)*0.715/2 + runif(length(z1), 
		min = -0.13, max = .13)), cex = 1.2)

fitness <- summary(lsmeans(z, "z1"))$lsmean
deltaP(pA, freq, fitness)
	# change in allele A 
	       # -0.09972984 

# Predicted frequency of allele B (low allele) in the next generation:
unname( (1 - pA) - deltaP(pA, freq, fitness) )
# [1] 0.6041941


# ------------
# Bootstrap SE's for "s" and for predicted delta p.
# Based on predicted genotypes, not on actual genotypes
# Resampling conditions on both family and genotype, so "freq" and "pA" stay the same

# setup
B <- 1000 # number of replicates
x <- pull.genoprob(rqtlAll.gp, chr = 4)
y <- rqtlAll$pheno$noffspring.p0.6
x1 <- x[, grep("chrIV:12815024", colnames(x))]
z1 <- apply(x1, 1, function(x){which(x == max(x))})
z1 <- recode(z1, c("1", "2", "3"), c("MM", "MF", "FF"))
z1 <- factor(z1, levels = c("MM", "MF", "FF"))
freq <- table(z1)
pA = (2 * freq[1] + freq[2]) / (2 * sum(freq)) # using observed p rather than fixing at 0.5
deltaP <- function(p, freq, fitness){
	# p is the frequency of the A allele
	# freq is a vector of frequencies of the 3 genotypes
	# fitness is a vector if average fitnesses of the 3 genotypes 
	q <- 1 - p
	wAA <- fitness[1]
	wAB <- fitness[2]
	wBB <- fitness[3]
	wbar <- sum(freq * fitness) / sum(freq)
	dP <- p*q*( p*(wAA - wAB) + q*(wAB - wBB) ) / wbar
	names(dP) <- "change in allele A"
	return(dP)
	}

s <- vector()
h <- vector()
s.em <- vector()
h.em <- vector()
predictLowAllele <- vector()

for(i in 1:B){
	# bootstrap replicate
	z <- tapply(y, INDEX = list(z1, family), FUN = sample, replace = TRUE)
	   # fem1.mal4  fem3.mal8  fem4.mal5  fem6.mal3  fem7.mal5  fem8.mal8 
	# MM Numeric,16 Numeric,12 Numeric,9  Numeric,5  Numeric,8  Numeric,10
	# MF Numeric,14 Numeric,27 Numeric,22 Numeric,13 Numeric,14 Numeric,12
	# FF Numeric,15 Numeric,15 Numeric,11 Numeric,9  Numeric,5  Numeric,7 
	z <- unlist(data.frame(z)) 
	
	# head(z)
		# fem1.mal4.MM1 fem1.mal4.MM2 fem1.mal4.MM3 fem1.mal4.MM4 fem1.mal4.MM5 fem1.mal4.MM6 
		            # 6             1             2             6             6             2 
	bootfam <- sub("([A-z0-9]+[.][A-z0-9]+)[.].*", "\\\\1", names(z))
	bootgeno <- factor(sub(".*([MF][MF])[0-9]+", "\\\\1", names(z)), levels = c("MM", "MF", "FF"))
	bootmean <- tapply(z, bootgeno, mean) # not accounting for family
	relFit <- bootmean/max(bootmean)
	s[i] <- 1 - relFit[1]
	h[i] <- (1 - relFit[2])/s[i]
	bootmean.em <- summary( emmeans(lm(z ~ bootfam + bootgeno), "bootgeno") )$emmean
	relFit.em <- bootmean.em/max(bootmean.em)
	s.em[i] <- 1 - relFit.em[1]
	h.em[i] <- (1 - relFit.em[2])/s.em[i]
	
	# Predicted frequency of allele B (low allele) in the next generation (no accounting for family)
	predictLowAllele[i] <- unname( (1 - pA) - deltaP(pA, freq, bootmean) )
	}

# results
mean(s)
	# [1] 0.4993793
sd(s)
	# [1] 0.08606801
mean(h)
	# [1] 0.9038182
sd(h)
	# [1] 0.181026

# bootstrap replicates
mean(s.em)
	# [1] 0.5139701
sd(s.em)
	# [1] 0.08825952
# 95% CI
quantile(s.em, probs = c(0.025,0.975)) # calculated separately in a second run
	     # 2.5%     97.5% 
	# 0.3149128 0.6609153

mean(h.em)
	# [1] 0.8939136
sd(h.em)
	# [1] 0.1887849
# 95% CI
quantile(h.em, probs = c(0.025,0.975)) # calculated separately in a second run
	     # 2.5%     97.5% 
	# 0.5589373 1.3021864 


mean(predictLowAllele)
	# [1] 0.6006851
sd(predictLowAllele)
	# [1] 0.02216601
quantile(predictLowAllele, probs = c(0.025,0.975)) # calculated separately in a second run
	     # 2.5%     97.5% 
	# 0.5541283 0.6427198

# ----
# body size vs INFERRED genotypes from genotype probabilities
x <- pull.genoprob(rqtlAll.gp, chr = 4)
x1 <- x[, grep("chrIV:12815024", colnames(x))]
z1 <- apply(x1, 1, function(x){which(x == max(x))})
z1 <- recode(z1, c("1", "2", "3"), c("MM", "MF", "FF"))
z1 <- factor(z1, levels = c("MM", "MF", "FF"))

# or take family into account
y <- rqtlAll$pheno$stdl
y[ y<2.6 ] <- NA # optional removal of outlier
z <- lm(y ~ family + z1)
emmeans(z, "z1")
	 # z1   emmean         SE  df lower.CL upper.CL
	 # MM 4.709447 0.04039910 215 4.629818 4.789076
	 # MF 4.860560 0.03114682 215 4.799167 4.921952
	 # FF 4.854922 0.04026212 215 4.775563 4.934281
	# Results are averaged over the levels of: family 

visreg(z, xvar = "z1", points = list(pch = 1, cex = 1.3),
	ylab = "F2 female standard length (cm)", xlab = "F2 female genotype at peak marker")

# include day
day <- rqtlAll$pheno$day
z <- lm(y ~ family + z1 + day)
emmeans(z, "z1")
	 # z1   emmean         SE  df lower.CL upper.CL
	 # MM 4.727459 0.03775790 214 4.653034 4.801884
	 # MF 4.857569 0.02901692 214 4.800374 4.914765
	 # FF 4.855510 0.03750314 214 4.781587 4.929433
	# Results are averaged over the levels of: family 

visreg(z, xvar = "z1", points = list(pch = 1, cex = 1.3),
	ylab = "F2 female standard length (cm)", xlab = "F2 female genotype at peak marker")

