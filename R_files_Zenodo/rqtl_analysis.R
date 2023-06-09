# Rqtl mapping of reproductive success and other traits of F2 females
# -------------------------------------------------------------------

# Notes:

# From Stanford Genome Browser:
# EDA is located at chrIV:12800313-12810434 (Human Proteins Mapped by Chained tBLASTn)
# WNT7B is located at chrIV:19900894-19905379 (Human Proteins Mapped by Chained tBLASTn)

# From Ensembl March 23 2018
# Gene: eda ENSGACG00000018311
# Descriptionectodysplasin A [Source:ZFIN;Acc:ZDB-GENE-050107-6]
# Synonyms si:ch73-223d24.5, nackt, nkt
# Location groupIV: 12,800,220-12,810,446 forward strand.
# This gene has 2 transcripts (splice variants), 92 orthologues and is a member of 1 Ensembl protein family.
# Nearest markers on map:
	# chrIV:10960835   49.036  4         4      4 10960835
	# chrIV:12811933   49.370  4         4      4 12811933
	# chrIV:12815024   49.685  4         4      4 12815024

# Gene: wnt7ba ENSGACG00000018959
# Description: wingless-type MMTV integration site family, member 7Ba [Source:ZFIN;Acc:ZDB-GENE-041210-178]
# Synonyms: Zwnt[c], si:ch211-239e6.2, wnt[c], wnt7, wnt7b
# Location: groupIV: 19,899,773-19,905,382 forward strand.
# This gene has 1 transcript (splice variant), 80 orthologues, 16 paralogues and is a member of 1 Ensembl protein family.
# Nearest marker on map:
	# chrIV:19271805          23.124

# --------------------------------------------------------------
# Read files into Rqtl

library(qtl)

# Read data from individual families into Rqtl
great1 <- read.cross(format="csvs", genfile="fem6.mal3.geno.csv", 
	phefile="fem6.mal3.pheno.csv", genotypes=NULL)
great2 <- read.cross(format="csvs", genfile="fem1.mal4.geno.csv", 
	phefile="fem1.mal4.pheno.csv", genotypes=NULL)
great3 <- read.cross(format="csvs", genfile="fem4.mal5.geno.csv", 
	phefile="fem4.mal5.pheno.csv", genotypes=NULL)
great4 <- read.cross(format="csvs", genfile="fem7.mal5.geno.csv", 
	phefile="fem7.mal5.pheno.csv", genotypes=NULL)
great5 <- read.cross(format="csvs", genfile="fem3.mal8.geno.csv", 
	phefile="fem3.mal8.pheno.csv", genotypes=NULL)
great6 <- read.cross(format="csvs", genfile="fem8.mal8.geno.csv", 
	phefile="fem8.mal8.pheno.csv", genotypes=NULL)

# Combine data frames of 6 families to a list
great <- list(fem6.mal3=great1, fem1.mal4=great2, fem4.mal5=great3, 
			fem7.mal5=great4, fem3.mal8=great5, fem8.mal8=great6)

# Read the data files for all families combined into Rqtl
rqtlAll <- read.cross(format="csvs", genfile="rqtlAll.geno.csv", 
	phefile="rqtlAll.pheno.csv", genotypes=NULL)

# Plot  linkagemap
plotMap(great[[1]], alternate.chrid=TRUE)

# Description of rqtl cross object
rqtlAll
	  # This is an object of class "cross".
	  # It is too complex to print, so we provide just this summary.
	    # F2 intercross
	
	    # No. individuals:    224 
	
	    # No. phenotypes:     5 
	    # Percent phenotyped: 100 100 100 100 100 
	
	    # No. chromosomes:    21 
	        # Autosomes:      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 
	
	    # Total markers:      388 
	    # No. markers:        22 16 22 39 17 16 26 16 21 21 20 12 14 17 13 19 24 16 2 16 19 
	    # Percent genotyped:  79.6 
	    # Genotypes (%):      AA:13.3  AB:28.6  BB:14.1  not BB:21.6  not AA:22.3 


# Missing genotype info
# "plot.missing" creates a grid of missing cases. Black indicates missing
plotMissing(rqtlAll, alternate.chrid=TRUE) 
	
# "plot.info" gives a measure of information missing
plotInfo(rqtlAll, map.function="kosambi", step=0, alternate.chrid=TRUE)

# Calculate and plot the fraction of genotypes missing for each marker
z <- pull.geno(rqtlAll)
z1 <- apply(z,2,function(x){length(x[is.na(x)])})

# Some trickery to get Rqtl to plot the fraction missing for each marker
z2 <- scanone(rqtlAll, pheno.col=3, method="mr")[,-3]
z2$fraction.missing <- z1/nrow(z)
for(i in names(rqtlAll$geno)){
	plot(z2, chr=i, ylim=c(0,1),show.marker.names=TRUE, col="blue",
			main = paste("chr",i))		
	}

# --------------------------------------------------------------
# Create family indicator variable

# Family sizes
lapply(great, nind)
	# $fem6.mal3
	# [1] 27
	
	# $fem1.mal4
	# [1] 45
	
	# $fem4.mal5
	# [1] 42
	
	# $fem7.mal5
	# [1] 27
	
	# $fem3.mal8
	# [1] 54
	
	# $fem8.mal8
	# [1] 29
	
# Create family variable
famsize <- sapply(great,function(x){nrow(x$pheno)})
family <- rep(names(great),famsize)

# ===============================================================================================
# Phenotypic associations

phenodat <- data.frame(rqtlAll$pheno, family, stringsAsFactors = FALSE)
phenodat$plate.morph <- recode(phenodat$platemorph, c("0","1","2"), c("low","partial","complete"))
phenodat$plate.morph <- factor(phenodat$plate.morph, levels = c("complete", "partial", "low"))

phenodat1 <- phenodat[phenodat$stdl > 3, ] # removes the stdl outlier

library(car)
library(visreg)
library(emmeans)
library(lattice)

# -----
# Body size differences between families
# ** crosses made with with male8 were July hatched rather than May hatched in 2006 **
# The male8 families are still the smallest
boxplot(stdl ~ family, data = phenodat, cex.axis = 0.8, ylab = "F2 female standard length (mm)", xlab = "F1 parent crosses") 

# -----
# Stats on number of offspring
mean(rqtlAll$pheno$noffspring.p0.6)
# [1] 2.116071

sd(rqtlAll$pheno$noffspring.p0.6)
# [1] 2.39673

se(rqtlAll$pheno$noffspring.p0.6)
# [1] 0.1601383


# --------------------------------------------------------------
# Plots and model fits of reproductive success vs size
# Using visreg() to plot model fits
# (further below we repeat this for all 3 genotypes at the peak fitness marker)

z1 <- lm (noffspring.p0.6 ~ stdl, data= phenodat1) # excludes outlier
z  <- lm (noffspring.p0.6 ~ stdl, data= phenodat)

visreg(z1, "stdl")
visreg(z, "stdl")

plot(jitter(noffspring.p0.6, amount = 0.1) ~ stdl, data=phenodat, 
		main="reproductive success and size (dashed removes outlier)",
		las = 1, ylab = "Estimated number of offspring", xlab = "Female standard length (mm)", bty = "l")
abline(z$coef)
abline(z1$coef, lty=2)

# Models with "family" and "day"
z <- lm(noffspring.p0.6 ~ family + stdl + day, data=phenodat)
visreg(z, "stdl")                   
visreg(z, "stdl", by = "family")                   

anova(z)
	#            Df  Sum Sq Mean Sq F value    Pr(>F)    
	# family      5  132.45  26.490  5.1981 0.0001585 ***
	# stdl        1   47.39  47.385  9.2983 0.0025803 ** 
	# day         1    0.38   0.384  0.0753 0.7840454    
	# Residuals 216 1100.76   5.096                      

Anova(z, type = 3)
	             # Sum Sq  Df F value   Pr(>F)   
	# (Intercept)   14.37   1  2.8194 0.094578 . 
	# family        69.76   5  2.7379 0.020191 * 
	# stdl          38.19   1  7.4933 0.006709 **
	# day            0.38   1  0.0753 0.784045   
	# Residuals   1100.76 216 

# -------
# Fiddling to get publication quality plot of fitness vs body size in visreg

col2rgb("firebrick")
	      # [,1]
	# red    178
	# green   34
	# blue    34

z1 <- lm(noffspring.p0.6 ~ family + stdl + day, data=phenodat1) # outlier removed
visreg(z1, "stdl")
visreg(z1, "stdl", by = "family")
visreg(z1, "stdl", scale = "response", rug = FALSE, xlim = range(phenodat1$stdl), 
	ylim = range(phenodat$noffspring.p0.6), bty = "l", 
	xlab = "F2 female standard length (cm)", ylab = "Estimated number of offspring")
	# points(jitter(noffspring.p0.6, amount = 0.1) ~ stdl, data=phenodat, cex = 1.2)
	points(jitter(noffspring.p0.6, amount = 0.1) ~ stdl, data=phenodat, cex = 1.3, pch = 16,
		col = rgb(178,  34,  34, maxColorValue = 255, alpha = 100))
	z2 <- visreg(glm(noffspring.p0.6 ~ family + stdl + day, data=phenodat1, family="quasipoisson"), plot = FALSE)[[2]]$fit 
	lines( exp(visregFit) ~ stdl, data = z2, lty = 2 )

anova(z1)
	#                                 Df  Sum Sq Mean Sq F value    Pr(>F)    
	#family[rqtlAll$pheno$stdl > 3]   5  129.29  25.859  5.0561 0.0002108 ***
	#stdl                              1   47.35  47.352  9.2588 0.0026356 ** 
	#day                               1    0.26   0.259  0.0506 0.8222534    
	#Residuals                       215 1099.58   5.114

Anova(z1, type = 3)
	                                 # Sum Sq  Df F value   Pr(>F)   
	# (Intercept)                       15.49   1  3.0287 0.083232 . 
	# family[rqtlAll$pheno$stdl > 3]   69.28   5  2.7091 0.021344 * 
	# stdl                              38.25   1  7.4783 0.006765 **
	# day                                0.26   1  0.0506 0.822253   
	# Residuals                       1099.58 215                    

# ----
# Plots and model fits of reproductive success vs size using glm
z1 <- glm(noffspring.p0.6 ~ family + stdl + day, data=phenodat1, family="quasipoisson")

visreg(z1, "stdl", scale = "response", rug = FALSE, xlim = range(phenodat$stdl), 
	ylim = range(phenodat$noffspring.p0.6), bty = "l", 
	xlab = "F2 female standard length (mm)", ylab = "Estimated number of offspring")
	points(jitter(noffspring.p0.6, amount = 0.1) ~ stdl, data=phenodat, cex = 1.2)

visreg(z1, "stdl", by = "family", scale = "response", ylim = range(phenodat1$noffspring.p0.6), rug = FALSE)

summary(z1)
# (Dispersion parameter for quasipoisson family taken to be 2.221627)

anova(z1, test="F") # F test is recommended for quasi
	       # Df Deviance Resid. Df Resid. Dev       F    Pr(>F)    
	# NULL                     222     566.00                      
	# family  5   59.771       217     506.23  5.3921 0.0001079 ***
	# stdl    1   24.314       216     481.92 10.9671 0.0010878 ** 
	# day     1    0.061       215     481.86  0.0274 0.8686998    

Anova(z1, type = 1, test = "F") # doesn't exist                 

Anova(z1, type = 3, test = "F")
	
# # 	              SS  Df      F  Pr(>F)   
	# family     30.01   5 2.7070 0.02143 * 
	# stdl       20.12   1 9.0760 0.00290 **
	# day         0.06   1 0.0274 0.86870   
	# Residuals 476.63 215                  


# -----
# Reproductive success by plate morph

par(bty = "l")
stripchart(jitter(noffspring.p0.6, amount = 0.1) ~ plate.morph, data=phenodat, 
 	vertical=TRUE, method="jitter", pch=1, main="reproductive success and plate morph", las = 1,
 	ylab = "Estimated number of offspring", xlab="Female plate morph")
z <- with(phenodat, tapply(noffspring.p0.6, plate.morph, mean))
lines(z ~ c(1:3), lwd=2)
z
# complete  partial      low 
# 1.761905 1.500000 3.140625 
 
z <- lm(noffspring.p0.6 ~ plate.morph + family + day, data = phenodat)
visreg(z, "plate.morph")                   
visreg(z, "plate.morph", by = "family")                   
z1 <- emmeans(z, "plate.morph")
	 # plate.morph   emmean        SE  df  lower.CL upper.CL
	 # complete    1.604527 0.2008094 215 1.2087194 2.000334
	 # partial     1.717994 0.4099828 215 0.9098936 2.526094
	 # low         3.018481 0.2812457 215 2.4641296 3.572833
visreg(z, "plate.morph", scale = "response", rug = FALSE, ylim = range(phenodat$noffspring.p0.6), bty = "l", 
	xlab = "F2 female lateral plate morph", ylab = "Estimated number of offspring")
	points(noffspring.p0.6 ~ I(0.145 + (as.integer(plate.morph) - 1)*0.715/2 + runif(nrow(phenodat), 
		min = -0.13, max = .13)), data=phenodat, cex = 1.2)

anova(z)
	             # Df  Sum Sq Mean Sq F value    Pr(>F)    
	# plate.morph   2   95.89  47.945  9.8074 8.393e-05 ***
	# family        5  124.64  24.928  5.0991 0.0001935 ***
	# day           1    9.39   9.387  1.9202 0.1672757    
	# Residuals   215 1051.07   4.889   

# same using glm
z <- glm(noffspring.p0.6 ~ plate.morph + family + day, family="quasipoisson", data = phenodat)
visreg(z, "plate.morph")                   
visreg(z, "plate.morph", scale = "response", rug = FALSE)                   
visreg(z, "plate.morph", by = "family")                   
z1 <- emmeans(z, "plate.morph")
	 # plate.morph    emmean        SE  df  asymp.LCL asymp.UCL
	 # complete    0.4210739 0.1062826 Inf 0.21276387 0.6293839
	 # partial     0.4350445 0.2165112 Inf 0.01069043 0.8593986
	 # low         1.0132265 0.1115714 Inf 0.79455054 1.2319024
visreg(z, "plate.morph", scale = "response", rug = FALSE, ylim = range(phenodat$noffspring.p0.6), bty = "l", 
	xlab = "F2 female lateral plate morph", ylab = "Estimated number of offspring")
	points(noffspring.p0.6 ~ I(0.145 + (as.integer(plate.morph) - 1)*0.715/2 + runif(nrow(phenodat), 
		min = -0.13, max = .13)), data=phenodat, cex = 1.2)


summary(z)
	# (Dispersion parameter for quasipoisson family taken to be 2.109607)


anova(z, test = "F")
	            # Df Deviance Resid. Df Resid. Dev       F    Pr(>F)    
	# NULL                          223     570.24                      
	# plate.morph  2   42.311       221     527.93 10.0281 6.858e-05 ***
	# family       5   57.648       216     470.28  5.4653 9.324e-05 ***
	# day          1    4.306       215     465.98  2.0412    0.1545    

# -----
# Reproductive success vs day (number of days past March 1, 2007) - numeric variable

plot(noffspring.p0.6 ~ jitter(day, amount=1), data= phenodat, 
		xlab="Day collected after March 1, 2007", main="reproductive success and collection day")
z <- lm(noffspring.p0.6 ~ day, data= phenodat)
abline(z)

boxplot(noffspring.p0.6 ~ day, data= phenodat, 
		xlab="Day collected after March 1, 2007", main="reproductive success and collection day")

anova(z)
	# Response: noffspring.p0.6
	           # Df  Sum Sq Mean Sq F value Pr(>F)
	# day         1   12.22 12.2223  2.1386  0.145
	# Residuals 222 1268.76  5.7151               


# -----
# Size vs day (number of days past March 1, 2007) - no outlier

plot(stdl ~ jitter(day,amount=1), data= phenodat1, 
		xlab="Day collected after March 1, 2007", main="size and collection date")
z <- lm(stdl ~ day, data= phenodat)
abline(z)
summary(z)
	            # Estimate Std. Error t value Pr(>|t|)    
	# (Intercept) 4.711015   0.044228  106.52  < 2e-16 ***
	# day         0.014166   0.004036    3.51 0.000543 ***

boxplot(stdl ~ day, data= phenodat1, 
		xlab="Day collected after March 1, 2007", main="size and collection date")

z <- lm(stdl ~ day, data= phenodat1)
visreg(z)
z <- lm(stdl ~ day + family, data= phenodat1)
visreg(z)
visreg(z, "day", by = "family")
anova(z)
	           # Df  Sum Sq Mean Sq F value    Pr(>F)    
	# day         1  2.3727  2.3727  27.992 2.983e-07 ***
	# family      5 23.3219  4.6644  55.028 < 2.2e-16 ***
	# Residuals 216 18.3090  0.0848                      

# Log size instead
plot(log(stdl) ~ jitter(day,amount=1), data= phenodat1, 
		xlab="Day collected after March 1, 2007", main="log size and collection date")
z <- lm(log(stdl) ~ day, data= phenodat)
abline(z)
summary(z)
	             # Estimate Std. Error t value Pr(>|t|)    
	# (Intercept) 1.5445129  0.0098491  156.82  < 2e-16 ***
	# day         0.0030019  0.0008988    3.34 0.000983 ***

z <- lm(log(stdl) ~ day + family, data= phenodat)
visreg(z)
visreg(z, "day", by = "family")
anova(z)
           # Df  Sum Sq  Mean Sq F value    Pr(>F)    
# day         1 0.11614 0.116135  23.901 1.971e-06 ***
# family      5 1.25686 0.251373  51.734 < 2.2e-16 ***
# Residuals 217 1.05438 0.004859                      

# ------
# Size vs plate morph

z <- lm(stdl ~ plate.morph + family + day, data = phenodat1)
visreg(z, "plate.morph", points = list(cex=1.3, pch=1), bty = "l", 
	xlab = "Lateral plate morph", ylab = "F2 female standard length (mm)")                   
visreg(z, "plate.morph", by = "family")                   
emmeans(z, "plate.morph")
	 # plate.morph   emmean         SE  df lower.CL upper.CL
	 # complete    4.793787 0.02639783 214 4.741754 4.845820
	 # partial     4.864036 0.05390967 214 4.757774 4.970298
	 # low         4.860178 0.03734871 214 4.786560 4.933797
	# Results are averaged over the levels of: family 

anova(z)
	             # Df  Sum Sq Mean Sq F value    Pr(>F)    
	# plate.morph   2  0.5012  0.2506  2.9663   0.05361 .  
	# family        5 22.3798  4.4760 52.9839 < 2.2e-16 ***
	# day           1  3.0444  3.0444 36.0377 8.191e-09 ***
	# Residuals   214 18.0782  0.0845                      

Anova(z)
	            # Sum Sq  Df   F value    Pr(>F)    
	# (Intercept) 798.79   1 9455.5907 < 2.2e-16 ***
	# plate.morph   0.23   2    1.3655    0.2575    
	# family       23.07   5   54.6152 < 2.2e-16 ***
	# day           3.04   1   36.0377 8.191e-09 ***
	# Residuals    18.08 214                        


# -----
# plate morph vs date
plot(jitter(platemorph, amount=.1) ~ jitter(day,amount=1), data= phenodat,
	xlab="Day collected after March 1, 2007", ylab="Plate morph", 
	main="plate morph and collection date")
z <- with(phenodat, tapply(platemorph, day, mean))
lines(as.numeric(names(z)),z, lwd=1)
z
#     4.5       23       37 
#1.286486 1.166667 2.000000


dev.off()

# ======================================================================
# QTL mapping - rqtlAll

# Analyze all families simultaneously, using family as a covariate

# Genotype probabilities
rqtlAll.gp <- calc.genoprob(rqtlAll, step=2, map.function="kosambi")

family.matrix <- model.matrix(~family)[,-1] # [,-1] removes the unwanted intercept column

# Same but with body size outlier deleted
rqtlAll2.gp <- subset(rqtlAll.gp, ind = rqtlAll.gp$pheno$stdl > 3.0)

family.matrix2 <- family.matrix[rqtlAll.gp$pheno$stdl > 3.0, ]

family2 <- family[rqtlAll.gp$pheno$stdl > 3.0]

# ----------------------------------------------------------------------
# Test of "omnigenic" model
# Make a variable indicating each indicating how limnetic and how freshwater each individual is on each chromosome
# AA is marine, AB is het, BB is freshwater (Cranby)

# Each chromosome one at a time (except 19, the sex chromosome)
pfresh <- list()
for(i in c(1:18,20,21)){
	# i <- 1
	# Grab a chromosome
	z <- as.data.frame(pull.genoprob(rqtlAll.gp, chr = i), stringsAsFactors = FALSE)

	# Keep only the real markers
	z <- z[, substr(colnames(z), 1, 3) == "chr"]

	# Number of markers on chromosome
	nmarkers <- ncol(z)/3
	# nmarkers
		# [1] 22
	# sumBB is the sum of probabilities that genotype is BB across markers
	sumAA <- apply(z[, grep("AA", names(z))], 1, sum)
	sumAB <- apply(z[, grep("AB", names(z))], 1, sum)
	sumBB <- apply(z[, grep("BB", names(z))], 1, sum)
	# head(sumAA + sumAB + sumBB)
		# F2.027 F2.029 F2.032 F2.051 F2.063 F2.071 
		    # 22     22     22     22     22     22 
	# hist(sumAA, breaks = 100)
	# hist(sumAB, breaks = 100)
	# hist(sumBB, breaks = 100)
	
	# pfresh is the proportion of freshwater alleles possessed by an individual
	pfresh[[i]] <- ( 0.5*sumAB + sumBB ) / (ncol(z)/3)
	# hist(pfresh[[i]])
	plot(rqtlAll.gp$pheno$noffspring.p0.6 ~ pfresh[[i]], main = paste("chr", i, "nmarkers =",nmarkers),
		ylab = "Number of surviving offspring", xlab = "Proportion feshwater alleles")
	z1 <- lm(rqtlAll.gp$pheno$noffspring.p0.6 ~ pfresh[[i]])
	z2 <- loess(rqtlAll.gp$pheno$noffspring.p0.6[order(pfresh[[i]])] ~ pfresh[[i]][order(pfresh[[i]])], span = 0.75)
	abline(z1)
	lines(pfresh[[i]][order(pfresh[[i]])], predict(z2))
	z3 <- lm(rqtlAll.gp$pheno$noffspring.p0.6 ~ family + pfresh[[i]])
	cat(paste("\\nchr",i,"nmarkers =",nmarkers, "\\n"))
	print(cbind(anova(z3)))
	}

	# chr 1 nmarkers = 22 
	             # Df      Sum Sq   Mean Sq  F value       Pr(>F)
	# family        5  132.450269 26.490054 5.031644 0.0002203852
	# pfresh[[i]]   1    6.093879  6.093879 1.157500 0.2831793577
	# Residuals   217 1142.437994  5.264691       NA           NA
	
	# chr 2 nmarkers = 16 
	             # Df      Sum Sq   Mean Sq  F value       Pr(>F)
	# family        5  132.450269 26.490054 5.009815 0.0002301978
	# pfresh[[i]]   1    1.115865  1.115865 0.211033 0.6464188118
	# Residuals   217 1147.416009  5.287631       NA           NA
	
	# chr 3 nmarkers = 22 
	             # Df      Sum Sq   Mean Sq   F value      Pr(>F)
	# family        5  132.450269 26.490054 5.0113210 0.000229507
	# pfresh[[i]]   1    1.460733  1.460733 0.2763377 0.599648472
	# Residuals   217 1147.071141  5.286042        NA          NA
	
	# chr 4 nmarkers = 39 
	             # Df     Sum Sq   Mean Sq   F value       Pr(>F)
	# family        5  132.45027 26.490054  5.330474 0.0001213818
	# pfresh[[i]]   1   70.13969 70.139688 14.113893 0.0002211138
	# Residuals   217 1078.39219  4.969549        NA           NA
	
	# chr 5 nmarkers = 17 
	             # Df       Sum Sq    Mean Sq   F value       Pr(>F)
	# family        5  132.4502691 26.4900538 5.0059078 0.0002319996
	# pfresh[[i]]   1    0.2203404  0.2203404 0.0416384 0.8385018359
	# Residuals   217 1148.3115334  5.2917582        NA           NA
	
	# chr 6 nmarkers = 16 
	             # Df     Sum Sq   Mean Sq  F value       Pr(>F)
	# family        5  132.45027 26.490054 5.051124 0.0002119821
	# pfresh[[i]]   1   10.49977 10.499766 2.002096 0.1585152834
	# Residuals   217 1138.03211  5.244388       NA           NA
	
	# chr 7 nmarkers = 26 
	             # Df       Sum Sq    Mean Sq    F value       Pr(>F)
	# family        5  132.4502691 26.4900538 5.00714457 0.0002314277
	# pfresh[[i]]   1    0.5039702  0.5039702 0.09526034 0.7578894046
	# Residuals   217 1148.0279036  5.2904512         NA           NA
	
	# chr 8 nmarkers = 16 
	             # Df        Sum Sq     Mean Sq       F value       Pr(>F)
	# family        5  132.45026911 26.49005382 5.00494858172 0.0002324441
	# pfresh[[i]]   1    0.00025837  0.00025837 0.00004881564 0.9944317895
	# Residuals   217 1148.53161538  5.29277242            NA           NA
	
	# chr 9 nmarkers = 21 
	             # Df      Sum Sq   Mean Sq   F value       Pr(>F)
	# family        5  132.450269 26.490054 5.0130077 0.0002287358
	# pfresh[[i]]   1    1.846696  1.846696 0.3494708 0.5550281188
	# Residuals   217 1146.685178  5.284263        NA           NA
	
	# chr 10 nmarkers = 21 
	             # Df      Sum Sq   Mean Sq   F value       Pr(>F)
	# family        5  132.450269 26.490054 5.0123692 0.0002290274
	# pfresh[[i]]   1    1.700607  1.700607 0.3217839 0.5711239108
	# Residuals   217 1146.831266  5.284937        NA           NA
	
	# chr 11 nmarkers = 20 
	             # Df      Sum Sq   Mean Sq  F value       Pr(>F)
	# family        5  132.450269 26.490054 5.034328 0.0002192083
	# pfresh[[i]]   1    6.702794  6.702794 1.273839 0.2602939470
	# Residuals   217 1141.829080  5.261885       NA           NA
	
	# chr 12 nmarkers = 12 
	             # Df       Sum Sq    Mean Sq    F value      Pr(>F)
	# family        5  132.4502691 26.4900538 5.00554186 0.000232169
	# pfresh[[i]]   1    0.1363862  0.1363862 0.02577145 0.872609092
	# Residuals   217 1148.3954875  5.2921451         NA          NA
	
	# chr 13 nmarkers = 14 
	             # Df       Sum Sq    Mean Sq   F value       Pr(>F)
	# family        5  132.4502691 26.4900538 5.0081401 0.0002309684
	# pfresh[[i]]   1    0.7321838  0.7321838 0.1384247 0.7102151452
	# Residuals   217 1147.7996899  5.2893995        NA           NA
	
	# chr 14 nmarkers = 17 
	             # Df         Sum Sq      Mean Sq      F value       Pr(>F)
	# family        5  132.450269111 26.490053822 5.0049576861 0.0002324398
	# pfresh[[i]]   1    0.002347632  0.002347632 0.0004435552 0.9832165480
	# Residuals   217 1148.529526113  5.292762793           NA           NA
	
	# chr 15 nmarkers = 13 
	             # Df      Sum Sq   Mean Sq  F value       Pr(>F)
	# family        5  132.450269 26.490054 5.032905 0.0002198313
	# pfresh[[i]]   1    6.380149  6.380149 1.212179 0.2721207573
	# Residuals   217 1142.151725  5.263372       NA           NA
	
	# chr 16 nmarkers = 19 
	             # Df      Sum Sq   Mean Sq   F value       Pr(>F)
	# family        5  132.450269 26.490054 5.0135305 0.0002284973
	# pfresh[[i]]   1    1.966258  1.966258 0.3721356 0.5424807604
	# Residuals   217 1146.565616  5.283713        NA           NA
	
	# chr 17 nmarkers = 24 
	             # Df       Sum Sq    Mean Sq    F value       Pr(>F)
	# family        5  132.4502691 26.4900538 5.00592872 0.0002319899
	# pfresh[[i]]   1    0.2251365  0.2251365 0.04254492 0.8367779352
	# Residuals   217 1148.3067373  5.2917361         NA           NA
	
	# chr 18 nmarkers = 16 
	             # Df         Sum Sq      Mean Sq     F value       Pr(>F)
	# family        5  132.450269111 26.490053822 5.004990987 0.0002324244
	# pfresh[[i]]   1    0.009989452  0.009989452 0.001887392 0.9653874526
	# Residuals   217 1148.521884294  5.292727577          NA           NA
	
	# chr 20 nmarkers = 16 
	             # Df     Sum Sq   Mean Sq  F value       Pr(>F)
	# family        5  132.45027 26.490054 5.039211 0.0002170823
	# pfresh[[i]]   1    7.80937  7.809370 1.485579 0.2242270077
	# Residuals   217 1140.72250  5.256786       NA           NA
	
	# chr 21 nmarkers = 19 
	             # Df     Sum Sq   Mean Sq  F value       Pr(>F)
	# family        5  132.45027 26.490054 5.084709 0.0001982397
	# pfresh[[i]]   1   18.01648 18.016484 3.458225 0.0642912446
	# Residuals   217 1130.51539  5.209748       NA           NA
	
# All chromosomes combined (except 4 and 19)
# First, no weighting of chromosomes
pfreshsum <- rep(0, length(pfresh[[1]]))
for(i in c(1:3,5:18,20,21)){
	pfreshsum <- pfreshsum + pfresh[[i]]
	}
pfreshsum <- pfreshsum/19
plot(rqtlAll.gp$pheno$noffspring.p0.6 ~ pfreshsum,
	ylab = "Number of surviving offspring", xlab = "Proportion freshwater alleles")
z1 <- lm(rqtlAll.gp$pheno$noffspring.p0.6 ~ pfreshsum)
z2 <- loess(rqtlAll.gp$pheno$noffspring.p0.6[order(pfreshsum)] ~ pfreshsum[order(pfreshsum)], span = 0.75)
abline(z1)
lines(pfreshsum[order(pfreshsum)], predict(z2))
z3 <- lm(rqtlAll.gp$pheno$noffspring.p0.6 ~ family + pfreshsum)
anova(z3)
	# Response: rqtlAll.gp$pheno$noffspring.p0.6
	           # Df  Sum Sq Mean Sq F value    Pr(>F)    
	# family      5  132.45 26.4901  5.0323 0.0002201 ***
	# pfreshsum   1    6.25  6.2501  1.1873 0.2770777    
	# Residuals 217 1142.28  5.2640                      

# Repeat, but weight by length of chromosomes
# Chromosome length by Glazer et al 2015 Table 2

chr <- 1:21
Mb <- c(29.63,23.7,17.8,34.14,15.56,18.85,30.84,20.53,20.58,18.03,17.64,20.76,20.74,16.17,17.32,19.52,20.25,15.99,20.61,20.45,17.35)
x <- data.frame(chr, Mb)
x

pfreshsum <- rep(0, length(pfresh[[1]]))
genomesum <- 0
for(i in c(1:3,5:18,20,21)){
	pfreshsum <- pfreshsum + pfresh[[i]] * x$Mb[i]
	genomesum <- genomesum + x$Mb[i]
	}
pfreshsum <- pfreshsum/genomesum
plot(rqtlAll.gp$pheno$noffspring.p0.6 ~ pfreshsum,
	ylab = "Number of surviving offspring", xlab = "Proportion freshwater alleles")
z1 <- lm(rqtlAll.gp$pheno$noffspring.p0.6 ~ pfreshsum)
z2 <- loess(rqtlAll.gp$pheno$noffspring.p0.6[order(pfreshsum)] ~ pfreshsum[order(pfreshsum)], span = 0.75)
abline(z1)
lines(pfreshsum[order(pfreshsum)], predict(z2))
z3 <- lm(rqtlAll.gp$pheno$noffspring.p0.6 ~ family + pfreshsum)
anova(z3)
	# Response: rqtlAll.gp$pheno$noffspring.p0.6
	           # Df  Sum Sq Mean Sq F value    Pr(>F)    
	# family      5  132.45 26.4901  5.0365 0.0002183 ***
	# pfreshsum   1    7.20  7.1982  1.3686 0.2433382    
	# Residuals 217 1141.33  5.2596                      

# Repeat this for each of the Eda genotypes separately
# Obtain the best-guess genotypes for Eda marker
tempcross <- argmax.geno(rqtlAll.gp, map.function="kosambi")
Eda <- pull.argmaxgeno(tempcross, 4)[,"chrIV:12811933"]
cbind(Eda,pull.geno(rqtlAll.gp, 4)[,"chrIV:12811933"]) # compare real to best guess - same
table(Eda, pull.pheno(rqtlAll.gp, "platemorph"))
       # platemorph
	# Eda  0  1  2
	  # 1  0  1 60 = MM
	  # 2  2 33 66 = MF
	  # 3 62  0  0 = FF
	  
plot(rqtlAll.gp$pheno$noffspring.p0.6[Eda == 1] ~ pfreshsum[Eda == 1], main = "Eda = MM",
	ylab = "Number of surviving offspring", xlab = "Proportion freshwater alleles")
anova(lm(rqtlAll.gp$pheno$noffspring.p0.6[Eda == 1] ~ family[Eda == 1] + pfreshsum[Eda == 1]))
	# Response: rqtlAll.gp$pheno$noffspring.p0.6[Eda == 1]
	                    # Df  Sum Sq Mean Sq F value Pr(>F)
	# family[Eda == 1]     5  24.637  4.9275  1.7543 0.1380
	# pfreshsum[Eda == 1]  1   2.604  2.6038  0.9270 0.3399
	# Residuals           54 151.677  2.8088               

plot(rqtlAll.gp$pheno$noffspring.p0.6[Eda == 2] ~ pfreshsum[Eda == 2], main = "Eda = MF",
	ylab = "Number of surviving offspring", xlab = "Proportion freshwater alleles")
anova(lm(rqtlAll.gp$pheno$noffspring.p0.6[Eda == 2] ~ family[Eda == 2] + pfreshsum[Eda == 2]))
	# Response: rqtlAll.gp$pheno$noffspring.p0.6[Eda == 2]
	                    # Df Sum Sq Mean Sq F value Pr(>F)
	# family[Eda == 2]     5  17.60  3.5202  1.0177 0.4118
	# pfreshsum[Eda == 2]  1   3.02  3.0239  0.8742 0.3522
	# Residuals           94 325.14  3.4589               

plot(rqtlAll.gp$pheno$noffspring.p0.6[Eda == 3] ~ pfreshsum[Eda == 3], main = "Eda = FF",
	ylab = "Number of surviving offspring", xlab = "Proportion freshwater alleles")
anova(lm(rqtlAll.gp$pheno$noffspring.p0.6[Eda == 3] ~ family[Eda == 3] + pfreshsum[Eda == 3]))
	# Response: rqtlAll.gp$pheno$noffspring.p0.6[Eda == 3]
	                    # Df Sum Sq Mean Sq F value   Pr(>F)   
	# family[Eda == 3]     5 192.29  38.458  4.6153 0.001378 **
	# pfreshsum[Eda == 3]  1   1.69   1.686  0.2023 0.654633   
	# Residuals           55 458.30   8.333                    



# --------------------------------------------------------------------
# Main QTL analysis

# Create an Eda genotype matrix
z <- pull.genoprob(rqtlAll.gp, chr = 4, omit.first.prob = TRUE)
eda.matrix <- z[,grep("chrIV:12815024", colnames(z))]

# Create Wnt variable here too:
z <- pull.genoprob(rqtlAll.gp, chr = 4, omit.first.prob = TRUE)
wnt.matrix <- z[,grep("chrIV:19271805", colnames(z))]

# Genome scan for qtl holding family identity as a covariate
rqtlAll.hk <- scanone(rqtlAll.gp, pheno.col=c(2:4), method="hk", addcovar= family.matrix)
head(rqtlAll.hk)
	               # chr   pos noffspring.p0.6        stdl platemorph
	# chrI:913033      1 0.000       0.3105307 0.044092405 0.04164234
	# c1.loc2          1 2.000       0.2499719 0.015321436 0.05167422
	# chrUn:37631434   1 3.954       0.2161730 0.001432632 0.06192761
	# c1.loc4          1 4.000       0.2162101 0.001246312 0.06213832
	# c1.loc6          1 6.000       0.2172834 0.009961621 0.06473771
	# c1.loc8          1 8.000       0.2139186 0.054251535 0.04667767

summary(rqtlAll.hk, pheno.col = "noffspring.p0.6")
plot(rqtlAll.hk, lodcolumn=1, alternate.chrid=TRUE, bandcol="gray90", 
 		main="no.offspring, family is covariate", ylab = "lod") # bty = "l"

# Where the peaks are
lodpeaks <- summary(rqtlAll.hk, format = "tabByCol", threshold = 3.5)
lodpeaks
	# noffspring.p0.6:
	               # chr  pos ci.low ci.high  lod
	# chrIV:12815024   4 49.7     42      72 4.33
	
	# stdl:
	               # chr  pos ci.low ci.high  lod
	# chrIV:31350187   4 64.1   51.1   71.97 3.53
	# c8.loc4          8  4.0    0.0    7.66 5.72
	
	# platemorph:
	               # chr  pos ci.low ci.high  lod
	# chrIV:12811933   4 49.4     49    49.7 95.7

# Lod score plots
plot(rqtlAll.hk, lodcolumn=1, chr = 4, main="", ylab = "lod", bty = "l", col = "firebrick")
plot(rqtlAll.hk, lodcolumn=2, chr = 4, main="", ylab = "lod", bty = "l", add = TRUE, col = "goldenrod1")
plot(rqtlAll.hk, lodcolumn=3, chr = 4, main="", ylab = "lod", bty = "l")

# markers on chrIV within 1.5 lod of lod peak
rqtlAll.hk[rqtlAll.hk$noffspring.p0.6 >= 4.33 - 1.5, ]
	               # chr    pos noffspring.p0.6     stdl platemorph
	# c4.loc44         4 44.000        2.934677 1.471995   42.76698
	# chrIV:4034002    4 44.039        2.939860 1.475557   42.81763
	# c4.loc46         4 46.000        3.228475 1.397277   54.33834
	# chrIV:6128193    4 47.766        3.204576 1.232787   57.61190
	# c4.loc48         4 48.000        3.622364 1.387663   63.21568
	# chrIV:8545605    4 48.284        4.101923 1.568555   68.82018
	# chrIV:10997988   4 48.375        4.246754 1.623763   70.15406
	# chrIV:10960835   4 49.036        3.932244 1.550100   81.95437
	# chrIV:12811933   4 49.370        4.161201 1.698145   95.70913
	# chrIV:12815024   4 49.685        4.327829 1.737419   91.88453
	# c4.loc50         4 50.000        4.206871 1.729891   91.67866
	# chrIV:11367975   4 50.875        3.796518 1.626819   86.28825
	# chrIV:9309735    4 51.072        3.729692 1.902212   84.18304
	# chrIV:15052901   4 51.355        3.517276 2.196336   78.94690
	# chrIV:13931030   4 51.886        3.456400 2.169753   78.32334
	# c4.loc52         4 52.000        3.424387 2.192867   76.69714
	# chrIV:15721538   4 52.560        3.165987 2.258592   65.14136
	# chrIV:15530121   4 52.807        3.035711 1.942964   60.85282
	# c4.loc62         4 62.000        3.062502 3.235179   31.23610
	# chrIV:29763654   4 63.612        3.207395 3.410527   25.54031
	# c4.loc64         4 64.000        3.199177 3.500850   24.55479
	# chrIV:31350187   4 64.123        3.193291 3.526232   24.22169
	# c4.loc66         4 66.000        3.319402 3.358659   21.35956
	# chrIV:32033500   4 67.123        3.158106 3.051179   18.22116
	# chrIV:30568387   4 67.656        3.116864 2.957303   17.30649
	# c4.loc68         4 68.000        3.136162 2.867496   16.92758
	# c4.loc70         4 70.000        3.046072 2.300470   13.62905

# Association between genotype and plate morph
table( pull.geno(rqtlAll.gp, 4)[,"chrIV:12811933"], pull.pheno(rqtlAll.gp, "platemorph") , useNA = "ifany")
	        # 0  1  2
	  # 1     0  1 58
	  # 2     1 31 63
	  # 3    62  0  0
	  # <NA>  1  2  5

# -----
# Are there any QTL when Eda is controlled statistically? No!
rqtlAll5.hk <- scanone(rqtlAll.gp, pheno.col=c(2), method="hk", addcovar=cbind(family.matrix, eda.matrix))
lodpeaks5 <- summary(rqtlAll5.hk, format = "tabByCol", threshold = 2.5) # use a low threshold
lodpeaks5
    # There were no LOD peaks above the threshold.
plot(rqtlAll5.hk, lodcolumn=1, alternate.chrid=TRUE, bandcol="gray90", 
 		main="no.offspring, family is covariate", ylab = "lod") # bty = "l"

# ---
# Redo after removing body size outlier
rqtlAll2.gp <- subset(rqtlAll.gp, ind = rqtlAll.gp$pheno$stdl > 3.0)
family.matrix2 <- family.matrix[rqtlAll.gp$pheno$stdl > 3.0, ]
family2 <- family[rqtlAll.gp$pheno$stdl > 3.0]
wnt.matrix2 <- wnt.matrix[rqtlAll.gp$pheno$stdl > 3.0, ]
eda.matrix2 <- eda.matrix[rqtlAll.gp$pheno$stdl > 3.0, ]

rqtlAll2.hk <- scanone(rqtlAll2.gp, pheno.col=c(2:4), method="hk", addcovar=family.matrix2)

lodpeaks2 <- summary(rqtlAll2.hk, format = "tabByCol", threshold = 3.5)
lodpeaks2
	# noffspring.p0.6:
	               # chr  pos ci.low ci.high lod
	# chrIV:12815024   4 49.7     42      72 4.5
	
	# stdl:
	                # chr   pos ci.low ci.high  lod
	# chrIV:31350187    4 64.12     58   71.97 4.40
	# chrVIII:1929053   8  3.52      0    7.66 7.12
	
	# platemorph:
	               # chr  pos ci.low ci.high  lod
	# chrIV:12811933   4 49.4     49    49.7 94.7

# Including more covariates - family, stdl
rqtlAll3.hk <- scanone(rqtlAll2.gp, pheno.col=c(2:4), method="hk", 
		addcovar=cbind(family.matrix2, rqtlAll2.gp$pheno$stdl))

lodpeaks3 <- summary(rqtlAll3.hk, format = "tabByCol", threshold = 3.5)
lodpeaks3[-2]
	# $noffspring.p0.6
	               # chr    pos ci.low ci.high      lod
	# chrIV:12815024   4 49.685     42      70 4.196753
	
	# $platemorph
	               # chr   pos ci.low ci.high      lod
	# chrIV:12811933   4 49.37 49.036  49.685 94.23854

# Including more covariates - family, day, stdl
rqtlAll3.hk <- scanone(rqtlAll2.gp, pheno.col=c(2:4), method="hk", 
		addcovar=cbind(family.matrix2, rqtlAll2.gp$pheno$day, rqtlAll2.gp$pheno$stdl))

lodpeaks3 <- summary(rqtlAll3.hk, format = "tabByCol", threshold = 3.5)
lodpeaks3[-2]
	# $noffspring.p0.6
	               # chr    pos ci.low ci.high      lod
	# chrIV:12815024   4 49.685     42      70 4.203018
	
	# $platemorph
	               # chr   pos ci.low ci.high      lod
	# chrIV:12811933   4 49.37 49.036  49.685 94.43911

# ---
# Test of Eda effects after taking WNT into account
# Redo after including Wnt7b (chrIV near 19,900,000) marker as covariate
rqtlAll4.hk <- scanone(rqtlAll2.gp, pheno.col=c(2:4), method="hk", 
					addcovar=cbind(family.matrix2, wnt.matrix2))
lodpeaks4 <- summary(rqtlAll4.hk, format = "tabByCol", threshold = 2.5) # use a low threshold
lodpeaks4
	# noffspring.p0.6:
	               # chr pos ci.low ci.high  lod
	# chrIV:27614532   4  54   47.8    55.6 2.71
	
	# stdl:
	                 # chr   pos ci.low ci.high  lod
	# c4.loc66           4 66.00   58.0   73.95 3.20
	# chrVIII:1929053    8  3.52    0.0    6.00 7.01
	# c16.loc32         16 32.00   20.3   53.60 2.53
	# c17.loc10         17 10.00    0.0   34.00 2.77
	# chrXVIII:4836241  18  5.39    0.0    7.26 3.48
	
	# platemorph:
	                # chr   pos ci.low ci.high   lod
	# chrIV:12811933    4 49.37  49.04    49.7 40.72
	# chrVII:26448674   7 27.97  10.00    28.0  3.03
	# chrXXI:8268451   21  7.96   6.34    30.0  4.52


# Test using a linear model whether eda is still significant after including wnt.
z <- lm(rqtlAll.gp$pheno$noffspring.p0.6 ~ family.matrix + eda.matrix)
anova(z)
	# Response: rqtlAll.gp$pheno$noffspring.p0.6
	               # Df  Sum Sq Mean Sq F value     Pr(>F)    
	# family.matrix   5  132.45  26.490  5.4455 0.00009675 ***
	# eda.matrix      2   97.78  48.888 10.0498 0.00006710 ***
	# Residuals     216 1050.76   4.865                       

z1 <- lm(rqtlAll.gp$pheno$noffspring.p0.6 ~ family.matrix + wnt.matrix)
z2 <- lm(rqtlAll.gp$pheno$noffspring.p0.6 ~ family.matrix + wnt.matrix + eda.matrix)
anova(z1,z2)
	  # Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
	# 1    216 1085.5                                
	# 2    214 1036.6  2    48.941 5.0518 0.007182 **


# ---------------------------------------------------------------------
# Lod thresholds using built-in permutation procedure, within-family permuting

z <- scanone(rqtlAll.gp, pheno.col = c(2:4), method = "hk", addcovar = family.matrix,
	n.perm = 10000, perm.strata = family)
summary(z, alpha=c(0.05, 0.04, 0.03, 0.02, 0.01, 0.009, 0.008)) 
	# LOD thresholds (10000 permutations)
	     # noffspring.p0.6 stdl platemorph
	# 5%              3.69 3.73       3.70
	# 4%              3.79 3.84       3.81
	# 3%              3.94 4.00       3.97
	# 2%              4.14 4.20       4.14
	# 1%              4.40 4.57       4.48
	# 0.9%            4.44 4.63       4.51
	# 0.8%            4.52 4.68       4.57

# Again, after removing body size outlier
hist(rqtlAll.gp$pheno$stdl, col = "red", breaks = 20)
z <- scanone(rqtlAll2.gp, pheno.col = c(2:4), method = "hk", addcovar = family.matrix2,
	n.perm = 1000, perm.strata = family2)
summary(z, alpha=c(0.05, 0.04, 0.03, 0.02, 0.01, 0.009, 0.008)) # new, using redos
# redo!

# noffspring and plates but stdl is a covariate and stdl outlier is removed
z <- scanone(rqtlAll2.gp, pheno.col = c(2,4), method = "hk", 
	addcovar = cbind(family.matrix2, rqtlAll2.gp$pheno$stdl), n.perm = 10000, perm.strata = family2)
summary(z, alpha=c(0.05, 0.04, 0.03, 0.02, 0.01, 0.009, 0.008)) # new, using redos
	# LOD thresholds (10000 permutations)
	     # noffspring.p0.6 platemorph
	# 5%              3.73       3.70
	# 4%              3.87       3.81
	# 3%              4.00       3.97
	# 2%              4.21       4.16
	# 1%              4.55       4.47
	# 0.9%            4.59       4.56
	# 0.8%            4.65       4.62


# ---
# More plots of lod profiles

#
# NOFFSPRING
# noffspring whole-genome qtl scan
plot(rqtlAll.hk, lodcolumn=1, alternate.chrid=TRUE, bandcol="gray90", 
 		main="no.offspring, family is covariate", ylab = "lod") # bty = "l"
lines(c(3.65, 3.65) ~ c(0,10000), lty = 2) # genome-wide lod threshold (see permutation test)

# noffspring, chromosome 4 only
plot(rqtlAll.hk, chr=4, lodcolumn=1, main="chr4", bty = "l")  # hint of two peaks on chromosome 4
lines(c(3.65, 3.65) ~ c(0,10000), lty = 2) # genome-wide lod threshold (see permutation test)

rqtlAll.hk[rqtlAll.hk$noffspring.p0.6 == max(rqtlAll.hk$noffspring.p0.6),1:3]
	               # chr    pos noffspring.p0.6
	# chrIV:12815024   4 49.685        4.327829
	

# noffspring with stdl as a covariate and stdl outlier is removed
x1 <- scanone(rqtlAll2.gp, pheno.col=2, method="hk", addcovar= cbind(family.matrix2, rqtlAll2.gp$pheno$stdl) )
plot(x1, chr = 4, main="no.offspring, stdl covariate (also family) no outlier", alternate.chrid=TRUE, bandcol="gray90")
lines(c(3.65, 3.65) ~ c(0,10000), lty = 2) # genome-wide lod threshold (see permutation test)

# noffspring with plate morph as a covariate (eliminates lod peak on chrIV)
x2 <- scanone(rqtlAll.gp, pheno.col=2, method="hk", 
	addcovar= cbind(family.matrix, model.matrix(~rqtlAll.gp$pheno$platemorph)[,-1]))
plot(x2, main="no.offspring, plate morph covariate (also family)", ylim = c(0,4), 
 			alternate.chrid=TRUE, bandcol="gray90")

#
# STDL
# stdl whole-genome qtl scan, after remove outlier at about 2.5 cm
plot(rqtlAll2.hk, lodcolumn=2, alternate.chrid=TRUE, bandcol="gray90", ylab = "lod",
				main="stdl no outlier, family is covariate")
lines(c(3.65, 3.65) ~ c(0,10000), lty = 2) # genome-wide lod threshold (see permutation test below)


#
# PLATES
plot(rqtlAll.hk, lodcolumn=3, alternate.chrid=TRUE, bandcol="gray90", 
 		main="plates, family is covariate", ylab = "lod") # bty = "l"
lines(c(3.65, 3.65) ~ c(0,10000), lty = 2) # genome-wide lod threshold (see permutation test below)
 
# plates with stdl as a covariate, no outlier
x1 <- scanone(rqtlAll2.gp, pheno.col=4, method="hk", addcovar= cbind(family.matrix2, rqtlAll2.gp$pheno$stdl) )

plot(x1, main="plates, stdl covariate (also family) no outlier", alternate.chrid=TRUE, bandcol="gray90")
lines(c(3.65, 3.65) ~ c(0,10000), lty = 2) # genome-wide lod threshold (see permutation test below)

dev.off()


# Plot phenotype against real + imputed genotypes at marker nearest the qtl
plotPXG(rqtlAll.gp, marker = find.marker(rqtlAll.gp, chr = lodpeaks[["noffspring.p0.6"]]$chr, 
		pos = lodpeaks[["noffspring.p0.6"]]$pos), pheno.col= "noffspring.p0.6", jitter=1, infer=TRUE)


# ---------
# Plots of noffspring for peak marker genotypes

# Mean number of offspring for each genotype at peak fitness marker
# Taking account of family and body size stdl.

noffspring <- pull.pheno(rqtlAll.gp, "noffspring.p0.6")
stdl <- pull.pheno(rqtlAll.gp, "stdl")

mean(noffspring)
	# [1] 2.116071

se(noffspring)

nearestmarker <- find.marker(rqtlAll.gp, chr=lodpeaks[["noffspring.p0.6"]]$chr, 
					pos=lodpeaks[["noffspring.p0.6"]]$pos)
nearestmarker
	# [1] "chrIV:12815024"

# Using known genotypes at best marker
gtype <- pull.geno(rqtlAll.gp, chr = lodpeaks[["noffspring.p0.6"]]$chr)[,nearestmarker]
tapply(noffspring, gtype, mean, na.rm=TRUE)
       # 1        2        3 
# 1.618182 1.777778 3.185185

tapply(noffspring, gtype, se, na.rm=TRUE) # standard error
        # 1         2         3 
# 0.2371135 0.1886031 0.4571614 

# Eda genotype frequencies in the F2 using the "chrIV:12815024" marker
table(gtype)
	 # 1  2  3 
	# 55 99 54

# ----
# Generate genotypes at the best marker - fill in missing genotypes using genoprob
# This uses all the genotype probabilities
x <- pull.genoprob(rqtlAll.gp, chr = 4)
gprob <- x[, grep(nearestmarker, colnames(x))] # genotype probabilities at peak marker
x <- round(gprob, 0) # assumes that the genotypes are known, ie rounds the probabilities
head(x)
	       # chrIV:12815024:AA chrIV:12815024:AB chrIV:12815024:BB
	# F2.027                 0                 1                 0
	# F2.029                 0                 1                 0
	# F2.032                 0                 1                 0
	# F2.051                 0                 0                 1
	# F2.063                 1                 0                 0
	# F2.071                 0                 1                 0

apply(x, 2, sum)
	# chrIV:12815024:AA chrIV:12815024:AB chrIV:12815024:BB 
	               # 60               102                62

summary( lm(noffspring ~ x - 1) )
                   # Estimate Std. Error     
# xchrIV:12815024:AA   1.6000     0.2980
# xchrIV:12815024:AB   1.7549     0.2286
# xchrIV:12815024:BB   3.2097     0.2932

# Mean number of offspring is 
(1.6000*60+ 1.7549*102 + 3.2097*62)/(60 + 102 + 62)
	# [1] 2.116077

gtype.imp <- apply(x, 1, function(x){which(x==1)}) # imputed genotype
# gt <- recode(gtype.imp, c(1,2,3), c("AA","Aa","aa"))
gt <- recode(gtype.imp, c(1,2,3), c("MM", "MF", "FF"))
gt <- factor(gt, levels = c("MM", "MF", "FF"))

# noffspring at marker chr4 imputed
par(bty = "l")
stripchart(jitter(noffspring, amount = 0.1) ~ gt, 
 	vertical=TRUE, method="jitter", pch=1, main="includes imputed genotypes", las = 1,
 	ylab = "Estimated number of offspring", xlab="Genotype at peak marker chr IV")
z <- tapply(noffspring, gt, mean)
lines(c(1:3), z, lwd=1)

z
	      # MM       MF       FF 
	# 1.600000 1.754902 3.209677 

# or 

# pub quality plot
col2rgb("firebrick")
	      # [,1]
	# red    178
	# green   34
	# blue    34
z <- lm(noffspring ~ family + gt) # take family into account
emmeans(z, "gt")
 # gt emmean    SE  df lower.CL upper.CL
 # MM   1.50 0.289 216     0.93     2.07
 # MF   1.67 0.223 216     1.23     2.11
 # FF   3.10 0.285 216     2.54     3.66
visreg(z, xvar = "gt") # conditional plot
visreg(z, xvar = "gt", scale = "response", ylim = range(y), rug=FALSE)
# Trick to superimpose the points - not conditioning on family
# visreg(z, xvar = "gt", scale = "response", rug = FALSE, ylim = range(noffspring), bty = "l", 
	# xlab = "F2 female genotype at peak marker", ylab = "Estimated number of offspring")
	# points(noffspring ~ I(0.145 + (as.integer(gt) - 1)*0.715/2 + runif(length(gt), 
		# min = -0.13, max = .13)), cex = 1.2)
visreg(z, xvar = "gt", scale = "response", rug = FALSE, ylim = range(noffspring), bty = "l", 
	xlab = "F2 female genotype", ylab = "Number of surviving offspring")
	points(jitter(noffspring, amount = 0.1) ~ I(0.145 + (as.integer(gt) - 1)*0.715/2 + runif(length(gt), 
		min = -0.13, max = .13)), cex = 1.3, pch = 16, col = rgb(178,  34,  34, maxColorValue = 255, alpha = 100))

summary( lm(noffspring ~ gt - 1) )
	     # Estimate Std. Error t value Pr(>|t|)    
	# gtMM   1.6000     0.2980   5.369 2.00e-07 ***
	# gtMF   1.7549     0.2286   7.678 5.16e-13 ***
	# gtFF   3.2097     0.2932  10.949  < 2e-16 ***

# Same using "emmeans"
library(emmeans)
z <- lm(noffspring ~ gt)
emmeans(z, ~gt)
	 # gt emmean    SE  df lower.CL upper.CL
	 # MM   1.60 0.298 221     1.01     2.19
	 # MF   1.75 0.229 221     1.30     2.21
	 # FF   3.21 0.293 221     2.63     3.79

z <- glm(noffspring ~ gt, family = quasipoisson(link = "log"))
lsmean <- emmeans(z, ~gt)
summary(lsmean)
	 # gt emmean    SE  df asymp.LCL asymp.UCL
	 # MM  0.470 0.155 Inf     0.165     0.775
	 # MF  0.562 0.114 Inf     0.339     0.785
	 # FF  1.166 0.108 Inf     0.955     1.378
	# Confidence level used: 0.95 

exp( summary(lsmean)[, c(2,5,6)])
	    # emmean asymp.LCL asymp.UCL
	# 1 1.600000  1.179923  2.169633
	# 2 1.754902  1.404072  2.193392
	# 3 3.209677  2.597733  3.965777

# Same but also controlling for stdl
stdl[ stdl < 2.6 ] <- NA # optional removal of outlier
z <- lm(noffspring ~ family + stdl + gt) # take family into account

# conditional plot
visreg(z, xvar = "gt", bty = "l", line = list(col = "black"),
	xlab = "Female F2 genotype", ylab = "Number of surviving offspring",
	points = list(cex = 1.2, pch = 1, col = "black", lwd = 1.5))

emmeans(z, "gt")
	 # gt emmean    SE  df lower.CL upper.CL
	 # MM   1.67 0.291 214     1.09     2.24
	 # MF   1.64 0.220 214     1.20     2.07
	 # FF   3.11 0.284 214     2.55     3.67

anova(z)
	           # Df  Sum Sq Mean Sq F value     Pr(>F)    
	# family      5  129.29  25.859  5.4998 0.00008729 ***
	# stdl        1   47.35  47.352 10.0713   0.001728 ** 
	# gt          2   93.67  46.835  9.9611 0.00007305 ***
	# Residuals 214 1006.17   4.702                       

# Many plots
# Make it a conditional plot
visreg(z, xvar = "gt", bty = "l", line = list(col = "black"),
		xlab = "Female F2 genotype", ylab = "Number of surviving offspring",
		points = list(cex = 1.2, pch = 1, col = "black", lwd = 1.5))

# Show the effect of stdl on number of offspring for the three genotypes
z <- lm(noffspring ~ family + stdl + gt) # take family into account
visreg(z, xvar = "stdl", by = "gt", layout = c(3,1), line = list(col = "black"),
	xlab = "F2 female standard length (cm)", ylab = "Number of surviving offspring (with family)",
	points = list(cex = 1.2, pch = 1, lwd = 1.5), band = FALSE)

# without family
z <- lm(noffspring ~ stdl + gt) # take family into account
visreg(z, xvar = "stdl", by = "gt", layout = c(3,1), line = list(col = "black"),
	xlab = "F2 female standard length (cm)", ylab = "Number of surviving offspring (without family)",
	points = list(cex = 1.3, pch = 1, lwd = 1.5), band = FALSE)

# Jitter=TRUE didn't work so this extra plot tried another strategy (don't use the lines)
z <- lm(jitter(noffspring, amount=0.1) ~ stdl + gt) # take family into account
visreg(z, xvar = "stdl", by = "gt", layout = c(3,1), line = list(col = "black"),
	xlab = "F2 female standard length (cm)", ylab = "Number of surviving offspring (without family)",
	points = list(cex = 1.3, pch = 16, lwd = 1.5), band = FALSE)

# include interaction
z <- lm(noffspring ~ stdl * gt) # take family into account
anova(z)
	           # Df  Sum Sq Mean Sq F value      Pr(>F)    
	# stdl        1  106.49 106.491 21.7744 0.000005364 ***
	# gt          2   89.49  44.747  9.1494   0.0001532 ***
	# stdl:gt     2   19.22   9.612  1.9653   0.1425964    
	# Residuals 217 1061.28   4.891                        
visreg(z, xvar = "stdl", by = "gt", layout = c(3,1), line = list(col = "black"),
	xlab = "F2 female standard length (cm)", ylab = "Number of surviving offspring (without family)",
	points = list(cex = 1.2, pch = 1, lwd = 1.5), band = FALSE)

z <- lm(noffspring ~ family2 + stdl * gt) # take family into account
	# Error in model.frame.default(formula = noffspring ~ family2 + stdl * gt,  : 
	  # variable lengths differ (found for 'family2')
anova(z)

# glm
z <- glm(noffspring ~ stdl + gt, family="quasipoisson")
visreg(z, xvar = "stdl", by = "gt", layout = c(3,1), line = list(col = "black"), scale = "response", rug = "FALSE",
	xlab = "F2 female standard length (cm)", ylab = "Number of surviving offspring (without family)",
	ylim=c(0,11), points = list(cex = 1.2, pch = 1, lwd = 1.5), band = FALSE)

z <- glm(noffspring ~ stdl * gt, family="quasipoisson")
anova(z, test = "Chisq")
	        # Df Deviance Resid. Df Resid. Dev Pr(>Chi)    
	# NULL                      222      566.0             
	# stdl     1    53.89       221      512.1 3.77e-07 ***
	# gt       2    37.11       219      475.0 0.000138 ***
	# stdl:gt  2     0.96       217      474.0 0.795058    

z <- glm(noffspring ~ family2 + stdl * gt, family="quasipoisson")
	# Error in model.frame.default(formula = noffspring ~ family2 + stdl * gt,  : 
	  # variable lengths differ (found for 'family2')
anova(z, test = "Chisq")


# Plot of stdl against peak marker genotype
# USes above trick to plot points on visreg plot

z <- glm(stdl ~ gt)
visreg(z, xvar = "gt", scale = "response", rug = FALSE, ylim = range(stdl, na.rm = TRUE), bty = "l", 
	xlab = "F2 female genotype", ylab = "F2 female standard length (cm)")
points(stdl ~ I(0.145 + (as.integer(gt) - 1)*0.715/2 + runif(length(gt), 
	min = -0.13, max = .13)), cex = 1.3, pch = 16, col = rgb(178,  34,  34, maxColorValue = 255, alpha = 100))

# ---------
# Models to estimate effects of best marker and body size markers

# Repeating the above (and with body size outlier removed)
z <- pull.genoprob(rqtlAll.gp, chr = 4, omit.first.prob = TRUE)
eda.matrix <- z[,grep("chrIV:12815024", colnames(z))]
eda.matrix2 <- eda.matrix[rqtlAll.gp$pheno$stdl > 3.0, ]
family2 <- family[rqtlAll.gp$pheno$stdl > 3.0]
family.matrix2 <- family.matrix[rqtlAll.gp$pheno$stdl > 3.0, ]
# rqtlAll2.hk <- scanone(rqtlAll2.gp, pheno.col=c(2:4), method="hk", addcovar=family.matrix2)
lodpeaks2 <- summary(rqtlAll2.hk, format = "tabByCol", threshold = 3.5)

noffspring <- pull.pheno(rqtlAll2.gp, "noffspring.p0.6")
stdl <- pull.pheno(rqtlAll2.gp, "stdl")
peakmarker <- find.marker(rqtlAll2.gp, chr=lodpeaks2[["noffspring.p0.6"]]$chr, # 4
					pos=lodpeaks2[["noffspring.p0.6"]]$pos) # 52.247
peakmarker
	# [1] "chrIV:12815024"
stdlmarker1 <- find.marker(rqtlAll2.gp, chr=lodpeaks2[["stdl"]]$chr[1], # 4
					pos=lodpeaks2[["stdl"]]$pos[1]) # 64.123 (modified in redo result)
z <- pull.genoprob(rqtlAll2.gp, chr = 4, omit.first.prob = TRUE)
stdlmarker1.matrix <- z[,grep(stdlmarker1, colnames(z))]
stdlmarker1
	# [1] "chrIV:32033500" - old result from better markers
	# [1] "chrIV:31350187" - new result from redo

stdlmarker2 <- find.marker(rqtlAll2.gp, chr=lodpeaks2[["stdl"]]$chr[2], # 8
					pos=lodpeaks2[["stdl"]]$pos[2]) # 3.524 (modified in redo result)
z <- pull.genoprob(rqtlAll2.gp, chr = 8, omit.first.prob = TRUE)
stdlmarker2.matrix <- z[,grep(stdlmarker2, colnames(z))]
stdlmarker2
	# [1] "chrVIII:1929053" (unchanged in redo results)

# Generate genotypes at the best marker - again but dropping body size outlier
x <- pull.genoprob(rqtlAll2.gp, chr = 4)
gprob <- x[, grep(nearestmarker, colnames(x))] # genotype probabilities at peak marker
x <- round(gprob, 0) # assumes that the genotypes are known, ie rounds the probabilities
apply(x, 2, sum)
	# chrIV:12815024:AA chrIV:12815024:AB chrIV:12815024:BB 
	               # 60               102                61 # not 62 because outlier dropped
gtype.imp <- apply(x, 1, function(x){which(x==1)}) # imputed genotype
gt <- recode(gtype.imp, c(1,2,3), c("MM", "MF", "FF"))
gt <- factor(gt, levels = c("MM", "MF", "FF"))
table(gt)
	 # MM  MF  FF 
	 # 60 102  61 # not 62 because outlier dropped

# Linear models with peak marker (size outlier not included)
z <- lm(noffspring ~ family2 + gprob)
anova(z)
	           # Df  Sum Sq Mean Sq F value    Pr(>F)    
	# family2     5  129.29  25.859  5.3183  0.000125 ***
	# gprob       2  101.81  50.906 10.4697 4.582e-05 ***
	# Residuals 215 1045.38   4.862                      

z <- lm(noffspring ~ family2 + stdl + gprob)
anova(z)
	           # Df  Sum Sq Mean Sq F value    Pr(>F)    
	# family2     5  129.29  25.859  5.4869 8.956e-05 ***
	# stdl        1   47.35  47.352 10.0477  0.001749 ** 
	# gprob       2   91.31  45.653  9.6871 9.389e-05 ***



# Test of contribution of size markers to noffspring
z <- lm(noffspring ~ family2 + eda.matrix2)
anova(z)
# old result from better markers
	             # Df  Sum Sq Mean Sq F value    Pr(>F)    
	# family2       5  129.29  25.859  5.3181 0.0001251 ***
	# eda.matrix2   2  101.77  50.883 10.4645 4.604e-05 ***
	# Residuals   215 1045.42   4.862  
# New result from redo	
	             # Df  Sum Sq Mean Sq F value     Pr(>F)    
	# family2       5  129.29  25.859  5.3183   0.000125 ***
	# eda.matrix2   2  101.81  50.906 10.4697 0.00004582 ***
	# Residuals   215 1045.38   4.862                       
                    

z <- lm(noffspring ~ family2 + eda.matrix2 + stdlmarker1.matrix + stdlmarker2.matrix)
anova(z)
# old result from better markers
	                    # Df  Sum Sq Mean Sq F value    Pr(>F)    
	# family2              5  129.29  25.859  5.2837 0.0001353 ***
	# eda.matrix2          2  101.77  50.883 10.3968 4.939e-05 ***
	# stdlmarker1.matrix   2    9.87   4.933  1.0079 0.3667519    
	# stdlmarker2.matrix   2    2.90   1.452  0.2968 0.7435196    
	# Residuals          211 1032.65   4.894                      
# new result from redos
	                    # Df  Sum Sq Mean Sq F value     Pr(>F)    
	# family2              5  129.29  25.859  5.2881  0.0001342 ***
	# eda.matrix2          2  101.81  50.906 10.4102 0.00004879 ***
	# stdlmarker1.matrix   2    8.31   4.154  0.8494  0.4291224    
	# stdlmarker2.matrix   2    5.28   2.639  0.5396  0.5837672    
	# Residuals          211 1031.79   4.890                       


# Effects of without including size to fit mean noffspring of genotypes at peak marker
options(digits=4)
z <- lm(noffspring ~ gt) # just gt
emmeans(z, "gt")
	 # gt emmean    SE  df lower.CL upper.CL
	 # MM   1.60 0.297 220     1.01     2.19
	 # MF   1.75 0.228 220     1.31     2.20
	 # FF   3.26 0.295 220     2.68     3.84
3.26 - 1.60
	# [1] 1.66
(1.60*60 + 1.75*102 + 3.26*61)/(60+102+61)
	# [1] 2.12
# s
1.6623/(3.26)
	# [1] 0.509908
summary(z)
	            # Estimate Std. Error t value Pr(>|t|)    
	# (Intercept)   1.6000     0.2973   5.381 1.89e-07 ***
	# gtMF          0.1549     0.3747   0.413     0.68    
	# gtFF          1.6623     0.4188   3.969 9.76e-05 ***
                    ******

# Adding family
z <- lm(noffspring ~ family2 + gt) # family included
emmeans(z, "gt") 
	 # gt emmean    SE  df lower.CL upper.CL
	 # MM   1.50 0.289 215     0.93     2.07
	 # MF   1.67 0.223 215     1.23     2.11
	 # FF   3.14 0.288 215     2.58     3.71
3.14 - 1.5
	# [1] 1.64
(1.50*60 + 1.67*102 + 3.14*61)/(60+102+61)
	# [1] 2.026368
# so add 0.094 to fitness of each genotype before calculation s so that mean fitness is the same
# s
1.6433/(3.14 + 0.094)
	# [1] 0.5081323

summary(z)
	                 # Estimate Std. Error t value Pr(>|t|)    
	# (Intercept)        2.1540     0.3948   5.456 1.33e-07 ***
	# family2fem3.mal8  -0.9747     0.4485  -2.173  0.03084 *  
	# family2fem4.mal5   0.6345     0.4771   1.330  0.18490    
	# family2fem6.mal3  -1.4147     0.5397  -2.621  0.00939 ** 
	# family2fem7.mal5  -0.7702     0.5401  -1.426  0.15527    
	# family2fem8.mal8  -1.4017     0.5319  -2.635  0.00902 ** 
	# gtMF               0.1730     0.3633   0.476  0.63437    
	# gtFF               1.6433     0.4035   4.073 6.53e-05 ***
                         ******

# Including size
z <- lm(noffspring ~ family2 + stdl + gt) # stdl too
anova(z)
	           # Df  Sum Sq Mean Sq F value    Pr(>F)    
	# family2     5  129.29  25.859  5.4998 8.729e-05 ***
	# stdl        1   47.35  47.352 10.0713  0.001728 ** 
	# gt          2   93.67  46.835  9.9611 7.305e-05 ***
emmeans(z, "gt")
	 # gt emmean    SE  df lower.CL upper.CL
	 # MM   1.67 0.291 214     1.09     2.24
	 # MF   1.64 0.220 214     1.20     2.07
	 # FF   3.11 0.284 214     2.55     3.67
3.11 - 1.67
	# [1] 1.44 # So adding size hardly diminishes the difference.
(1.67*60 + 1.64*102 + 3.11*61)/(60+102+61)
	# [1] 2.050179
# so add 0.07 to fitness of each genotype before calculating s so that mean fitness is the same as above, 2.12
# s
1.44760/(3.11 + 0.07)
	# [1] 0.4552201

summary(z)
	                 # Estimate Std. Error t value Pr(>|t|)    
	# (Intercept)      -4.52532    2.41632  -1.873 0.062458 .  
	# family2fem3.mal8 -0.38241    0.48954  -0.781 0.435573    
	# family2fem4.mal5  0.58693    0.46997   1.249 0.213076    
	# family2fem6.mal3 -1.48277    0.53191  -2.788 0.005787 ** 
	# family2fem7.mal5 -0.38842    0.54885  -0.708 0.479901    
	# family2fem8.mal8 -0.19987    0.67697  -0.295 0.768095    
	# stdl              1.34536    0.48037   2.801 0.005566 ** 
	# gtMF             -0.03029    0.36491  -0.083 0.933931    
	# gtFF              1.44760    0.40330   3.589 0.000411 ***
	                    *******
# Same, leaving out family
z <- lm(noffspring ~ stdl + gt) # stdl too
anova(z)
emmeans(z, "gt")
	 # gt emmean    SE  df lower.CL upper.CL
	 # MM   1.79 0.290 219     1.22     2.36
	 # MF   1.70 0.220 219     1.27     2.14
	 # FF   3.16 0.285 219     2.60     3.72
(1.79*60 + 1.70*102 + 3.16*61)/(60+102+61)
	# [1] 2.12359
# s
1.3694/(3.16)
	# [1] 0.433354
summary(z)
	            # Estimate Std. Error t value Pr(>|t|)    
	# (Intercept)  -5.1084     1.6269   -3.14  0.00192 ** 
	# stdl          1.4275     0.3408    4.19  4.1e-05 ***
	# gtMF         -0.0876     0.3660   -0.24  0.81112    
	# gtFF          1.3694     0.4099    3.34  0.00098 ***
                    ******

# ----------------
# Estimate effects and PVE at lod peaks - outlier is not removed

z <- summary(rqtlAll.hk,format="tabByCol", threshold=3.4)
z
	# noffspring.p0.6:
	               # chr  pos ci.low ci.high  lod
	# chrIV:12815024   4 49.7     42      72 4.33
	
	# stdl:
	               # chr  pos ci.low ci.high  lod
	# chrIV:31350187   4 64.1   51.1   71.97 3.53
	# c8.loc4          8  4.0    0.0    7.66 5.72
	
	# platemorph:
	               # chr  pos ci.low ci.high  lod
	# chrIV:12811933   4 49.4     49    49.7 95.7


lodmin <- 3.5
for(i in 1:length(z)){
		z[[i]]$trait <- names(z)[i]
		z[[i]]$name <- rownames(z[[i]])
		}
z1 <- do.call("rbind",z)
z1 <- z1[z1$lod >= lodmin,]
z1
	                    # chr    pos ci.low ci.high       lod           trait           name
	# noffspring.p0.6       4 49.685 42.000  71.970  4.327829 noffspring.p0.6 chrIV:12815024
	# stdl.chrIV:31350187   4 64.123 51.072  71.970  3.526232            stdl chrIV:31350187
	# stdl.c8.loc4          8  4.000  0.000   7.657  5.715246            stdl        c8.loc4
	# platemorph            4 49.370 49.036  49.685 95.709127      platemorph chrIV:12811933


# Broman and Sen say that estimated percent variance explained is 100 * ( 1 - 10^(-2*LOD/n) )

100 * ( 1 - 10^(-2 * z1$lod/nrow(rqtlAll.gp$pheno)) )
[1]  8.513151  6.992969 11.085822 86.021704

# Redo after removing stdl outlier
z <- summary(rqtlAll2.hk, format="tabByCol")
for(i in 1:length(z)){
		z[[i]]$trait <- names(z)[i]
		z[[i]]$name <- rownames(z[[i]])
		}
z1 <- do.call("rbind",z)
z1 <- z1[z1$lod >= lodmin,]
z1
	                               # chr    pos ci.low ci.high       lod           trait
	# noffspring.p0.6.chrIV:12815024   4 49.685 42.000  71.970  4.500375 noffspring.p0.6
	# stdl.chrIV:31350187              4 64.123 58.000  71.970  4.403281            stdl
	# stdl.chrVIII:1929053             8  3.524  0.000   7.657  7.123274            stdl
	# platemorph.chrIV:12811933        4 49.370 49.036  49.685 94.664221      platemorph
	                                          # name
	# noffspring.p0.6.chrIV:12815024  chrIV:12815024
	# stdl.chrIV:31350187             chrIV:31350187
	# stdl.chrVIII:1929053           chrVIII:1929053
	# platemorph.chrIV:12811933       chrIV:12811933


# PVE %var
100 * ( 1 - 10^(-2 * z1$lod/nrow(rqtlAll.gp$pheno)) )
	# [1]  8.837110  8.654956 13.622756 85.718174
