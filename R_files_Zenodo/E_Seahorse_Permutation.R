
setwd("C:/Users/agjon/OneDrive/Documents/Work/Rexamples/Lotta_Data")
sh_data <- read.csv("D_Seahorse_Overlap.csv")

library(lattice)

#Take a look at the data. 

dev.new(width=4.5,height=4.5)
trellis.par.set(strip.background=list(col="lightgrey"))

histogram(~ overlap_m2 | Sexes, data=sh_data,
	ylab = list("Percentage of pairs", cex=1), 
	xlab = list("Home range overlap", cex=1),
	scales=list(cex=c(1,1)),
	strip=strip.custom(var.name="Comparison",
		factor.levels=c("Female-Female","Female-Male","Mated Pairs","Male-Male"),
		strip.levels=rep(TRUE,4)),
	par.settings=simpleTheme(col="gray")
	)

#The pairs (FMP) have fewer overlaps of zero and more small range
#overlaps. They do not have a lot of mid-range overlaps.
#It also looks like the female-female overlaps tend to have
#more pairs with extreme range overlaps.

#Do the overlaps need to be corrected for mean home-range size?
#Females may have larger overlaps because they have larger home-ranges.


# Compare means using a permutation test
# First compare male-female pairs to random male-female overlap

sh_data_mf <- subset(sh_data, Sexes == "FM" | Sexes == "FMP")

obs_mean <- tapply(sh_data_mf$overlap_m2, sh_data_mf$pair, FUN=mean) 
obs_diff <- obs_mean[2] - obs_mean[1]
obs_diff
#The difference between known pairs and random pairs is 10.14014


nPerm <- 10000
permResult <- vector()

for (i in 1:nPerm) {
	permSample <- sample(sh_data_mf$overlap_m2, replace=FALSE)
	permMeans <- tapply(permSample, sh_data_mf$pair, FUN=mean)
	permResult[i] <- permMeans[2] - permMeans[1]
}

dev.new(width=6,height=4)
par(lwd=3)
hist(permResult, right=FALSE, main="Mated Pair Permutation Test",
	xlab=list("Difference in Home-Range Overlap", cex=1.2),
	ylab=list("No. Permutations", cex=1.2), lwd=3, cex.axis=1.2,
	xlim=c(-5,15)
)

proportion_below <- sum(as.numeric(permResult <= obs_diff))/nPerm
perm_p_val <- 2*(1 - proportion_below)
perm_p_val
#The mean overlap for mated pairs is significantly greater than the 
#mean overlap for randomly chosen male-female pairs, p=0.0004


#Use a permutation test to compare mean male-male overlap to
#mean female-female overlap

sh_data_ffmm <- subset(sh_data, Sexes == "FF" | Sexes == "MM")
sh_data_ffmm$Sexes <- factor(sh_data_ffmm$Sexes)
obs_mean <- tapply(sh_data_ffmm$overlap_m2, sh_data_ffmm$Sexes, FUN=mean)
obs_diff <- obs_mean[1] - obs_mean[2]
obs_diff

nPerm <- 10000
permResult <- vector()

for (i in 1:nPerm) {
	permSample <- sample(sh_data_ffmm$overlap_m2, replace=FALSE)
	permMeans <- tapply(permSample, sh_data_ffmm$Sexes, FUN=mean)
	permResult[i] <- permMeans[1] - permMeans[2]
}

hist(permResult, right=FALSE)

dev.new(width=6,height=4)
par(lwd=3)
hist(permResult, right=FALSE, main="Female-Male Permutation Test",
	xlab=list("Difference in Home-Range Overlap", cex=1.2),
	ylab=list("No. Permutations", cex=1.2), lwd=3, cex.axis=1.2,
	xlim=c(-3,3)
)


proportion_below <- sum(as.numeric(permResult <= obs_diff))/nPerm
perm_p_val <- 2*(1 - proportion_below)
perm_p_val

#The mean overlap among females is significantly greater than the mean
#overlap among males, p = 0, so p < 0.0001.


