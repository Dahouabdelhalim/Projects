################################################################################
################## PROTEIN ANALYSIS R SCRIPT ###################################
################################################################################


########## INSTALL REQUIRED PACKAGES ###########################################

if(!("NADA" %in% list.files(.libPaths()))){install.packages("NADA")}
if(!("vegan" %in% list.files(.libPaths()))){install.packages("vegan")}
if(!("mvtnorm" %in% list.files(.libPaths()))){install.packages("mvtnorm")}
if(!("ellipse" %in% list.files(.libPaths()))){install.packages("ellipse")}

########## DEFINE FUNCTIONS ####################################################

# FDR takes a  vector of (sorted from smallest to largest) p-values
# and returns the number of significant p-values after controlling
# for the false discovery rate
#(see Benjamini & Hochberg 1995)
FDR <- function(v, alpha = 0.05){
  m <- length(v)
  counter <- 1
  while(v[counter] <= (counter * alpha/m)) counter <- counter + 1
  return(counter - 1)
}

#glmmer fits lognormal distn to data and compares against null model using Chi
# squared analysis of deviance test
glmmer <- function(v, X){
  fit0 <- glm(v ~ 1, family = gaussian(link="identity"))
  fit1 <- glm(v ~ X, family = gaussian(link = "identity"))
  fit.anova <- anova(fit0, fit1, test = "Chisq")
  overall.pval <- fit.anova[[5]][2]
  names(overall.pval) <- "overall.pval"
  fit.aov <- aov(v ~ X)
  pw <- TukeyHSD(fit.aov)
  pw.pvals <- pw$X[, 4]
  names(pw.pvals) <- paste0(names(pw.pvals), ".pval")
  pw.mags <- exp(pw$X[, 1] + fit.aov$coefficients[1])/exp(fit.aov$coefficients[1])
  names(pw.mags) <- paste0(names(pw.mags), ".foldchange")
  res <- c(overall.pval, pw.pvals, pw.mags)
  return(res)
}

########## SET WORKING DIRECTORY AND LOAD DATA FILE ############################

# point to data file to set the working directory to its parent folder
setwd(dirname(file.choose()))

# Note: changes to files include removing columns and rows not representing measurements
# and changing cell types to Number with 0 decimal places (read.csv struggled with excel's
# scientific notation)

Y <- read.csv("totalTIC.csv", header = F, colClasses = 'numeric')
Y <- t(as.matrix(Y))
treatment <- paste(rep(c("B", "N", "C", "S"), each = 12), rep(1:6, each = 2, times = 4), sep = "")

# order by factor column
Y <- Y[order(treatment),]
treatment <- sort(treatment)
n <- nrow(Y)
p <- ncol(Y)
rownames(Y) <- 1:n
colnames(Y) <- 1:p

############ BODY ##############################################################

# calculate and apply scale factor to ensure all samples have equal sampling depth
tmp <- Y[, apply(Y, 2, function(v) all(v > 0))] #only retain columns without nondetects
fun <- Vectorize(function(y, z) median(tmp[y, ]/tmp[z, ]))
scale.mat <- outer(1:n, 1:n, fun)
scale.totals <- apply(scale.mat, 1, sum)
most.undersampled <- which.min(scale.totals)
print(paste("most undersampled replicate:", most.undersampled))
# should be replicate 13
Y <- Y * scale.mat[most.undersampled, ]

# remove any columns with more than 1/3 non-detects
include <- apply(Y, 2, function(v) sum(v == 0) < 0.33 * n)
Y <- Y[, include]
print(paste((p - sum(include)), "columns removed"))
p <- ncol(Y)
no.imputed <- apply(Y, 2, function(v) sum(v == 0))

# impute missing values using robust regression on order statistics (ROS)
# note these are ordered so will need to be resampled
Y.obs <- as.vector(Y)
Y.cens <- Y.obs == 0
no.censored <- sum(Y.cens)
Y.ros <- NADA::ros(obs = Y.obs, censored = Y.cens)
Y.mod <- Y.ros$modeled

# resample ROS modeled values and replace nondetects in Y matrix
set.seed(123456789)
Y.replacements <- sample(Y.mod[1:no.censored])
Y.i <- Y.obs
Y.i[Y.cens] <- Y.replacements
Y.i <- matrix(Y.i, nrow = n)
dimnames(Y.i) <- dimnames(Y)

# merge technical replicates and calculate p-values/model parameters
odds <- seq(1, 48, by = 2)
evens <- seq(2, 48, by = 2)
Y.i <- (Y.i[odds, ] + Y.i[evens, ])/2
rownames(Y.i) <- 1:(n/2)
Y.log <- log(Y.i)

# calculate omnibus p-value using the vegan package
X <- factor(rep(c("B", "C", "N", "S"), each = 6))
X <- relevel(X, ref = "N")
veg.dist <- vegan::vegdist(Y.log)
veg.adonis <- vegan::adonis(veg.dist ~ X)
veg.p <- veg.adonis$aov.tab[1,6]

# calculate univariate p-values and fold changes, and output to table
anovas <- apply(Y.log, 2, glmmer, X)
out <- t(anovas)
Protein = data.frame(Protein = as.integer(colnames(Y.log)))
out <- cbind(Protein, out)
out <- out[order(out$overall.pval),]
number.signif <- FDR(out$overall.pval, alpha = 0.1)
out <- out[1:number.signif,]
out <- out[order(out$Protein),]
out$no.imputed.techreps <- no.imputed[names(no.imputed) %in% out$Protein]

### output imputed totals
out2 <- Y.i[, colnames(Y.i) %in% out$Protein]
rownames(out2) <- treatment[odds]
write.csv(out2, "totals.csv")

## visualize with ordination diagram
view.protein <- colnames(Y.log) %in% out$Protein
prot.MDS <- vegan::metaMDS(Y.log, autotransform = FALSE)

plot(prot.MDS, type="n", ylim = c(-0.035, 0.035))
points(prot.MDS, display = "species", select = view.protein, cex = 0.6, pch = 3)
points(prot.MDS$points[1:6,], col = 'orange', pch = 20, cex = 1.5)
points(prot.MDS$points[7:12,], col = 'red', pch = 20, cex = 1.5)
points(prot.MDS$points[13:18,], col = 'blue', pch = 20, cex = 1.5)
points(prot.MDS$points[19:24,], col = 'green', pch = 20, cex = 1.5)

library(ellipse)
lines(ellipse(cov(prot.MDS$points[1:6,]),
              centre = c(mean(prot.MDS$points[1:6,][,1]), mean(prot.MDS$points[1:6,][,2])),
              level = 0.95), col='orange')
lines(ellipse(cov(prot.MDS$points[7:12,]),
              centre = c(mean(prot.MDS$points[7:12,][,1]), mean(prot.MDS$points[7:12,][,2])),
              level = 0.95), col='red')
lines(ellipse(cov(prot.MDS$points[13:18,]),
              centre = c(mean(prot.MDS$points[13:18,][,1]), mean(prot.MDS$points[13:18,][,2])),
              level = 0.95), col='blue')
lines(ellipse(cov(prot.MDS$points[19:24,]),
              centre = c(mean(prot.MDS$points[19:24,][,1]), mean(prot.MDS$points[19:24,][,2])),
              level = 0.95), col='green')
legend("topright", legend= paste("Stress =", round(prot.MDS$stress, 3)), bty="n")



## create NMDS ordination biplot as PDF document in working directory
pdf(file = "proteins.pdf", width = 8, height = 6)
plot(prot.MDS, type="n", ylim = c(-0.035, 0.035))
points(prot.MDS, display = "species", select = view.protein, cex = 0.6, pch = 3)
points(prot.MDS$points[1:6,], col = 'orange', pch = 20, cex = 1.5)
points(prot.MDS$points[7:12,], col = 'red', pch = 20, cex = 1.5)
points(prot.MDS$points[13:18,], col = 'blue', pch = 20, cex = 1.5)
points(prot.MDS$points[19:24,], col = 'green', pch = 20, cex = 1.5)

lines(ellipse(cov(prot.MDS$points[1:6,]),
              centre = c(mean(prot.MDS$points[1:6,][,1]), mean(prot.MDS$points[1:6,][,2])),
              level = 0.95), col='orange')
lines(ellipse(cov(prot.MDS$points[7:12,]),
              centre = c(mean(prot.MDS$points[7:12,][,1]), mean(prot.MDS$points[7:12,][,2])),
              level = 0.95), col='red')
lines(ellipse(cov(prot.MDS$points[13:18,]),
              centre = c(mean(prot.MDS$points[13:18,][,1]), mean(prot.MDS$points[13:18,][,2])),
              level = 0.95), col='blue')
lines(ellipse(cov(prot.MDS$points[19:24,]),
              centre = c(mean(prot.MDS$points[19:24,][,1]), mean(prot.MDS$points[19:24,][,2])),
              level = 0.95), col='green')
legend("topright", legend= paste("Stress =", round(prot.MDS$stress, 3)), bty="n")
dev.off()

## iterate 100 times to ensure imputation didn't affect the outcome
## this step takes a while (around 15 mins)
nrep <- 100
univars <- matrix(NA, ncol = ncol(Y.log), nrow = nrep)
multivars <- rep(NA, nrep)
set.seed(123456789)
for(i in 1:nrep){
  ## resample ROS modeled values and replace nondetects in Y matrix
  Y.replacements <- sample(Y.mod[1:no.censored])
  Y.i <- Y.obs
  Y.i[Y.cens] <- Y.replacements
  Y.i <- matrix(Y.i, nrow = n)
  dimnames(Y.i) <- dimnames(Y)
  ## average over technical replicates
  Y.i <- (Y.i[odds, ] + Y.i[evens, ])/2
  rownames(Y.i) <- 1:(n/2)
  Y.log <- log(Y.i)
  glms <- apply(Y.log, 2, glmmer, X)
  #univars[i,] = sapply(glm.list, function(x) x$overall.pval)
  univars[i, ] <- glms[1,]
  veg.dist.i <- vegan::vegdist(Y.log)
  veg.adonis.i <- vegan::adonis(veg.dist.i ~ X)
  veg.p.i <- veg.adonis.i$aov.tab[1,6]
  multivars[i] <- veg.p.i
}
colnames(univars) <- colnames(Y.log)
prop.sig <- apply(univars, 2, function(v) sum(v < 0.05))
names(prop.sig) <- colnames(univars)
prop.sig <- prop.sig[view.protein]
out$resamp.percent.signif <- prop.sig
write.csv(out, file = "pvals2.csv", row.names = FALSE)
# what percentage of multivariate pvals are significant?
sum(multivars < 0.05)/nrep
################################# END ##########################################
