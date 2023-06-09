# Activate Rqtl
library(qtl)
library(nortest)


rm(list=ls())
dev.off(dev.list()["RStudioGD"])

# upload .csv files using read.cross function (p195) 

phd <- read.cross(format=c("csv"), 
                  dir="D:/PhD/Crossing/Analyses/Rqtl/P.lanceolata", 
                  file="ALL_DATA.qtl.csv", na.strings=c("NA"), 
                  genotypes=NULL)

phd_N <- read.cross(format=c("csv"), 
                    dir="D:/PhD/Crossing/Analyses/Rqtl/P.lanceolata", 
                    file="ALL_DATA_N.qtl.csv", na.strings=c("NA"), 
                    genotypes=NULL)

phd_S <- read.cross(format=c("csv"), 
                    dir="D:/PhD/Crossing/Analyses/Rqtl/P.lanceolata", 
                    file="ALL_DATA_S.qtl.csv", na.strings=c("NA"), 
                    genotypes=NULL)

# Summarize data set
summary(phd)
summary(phd$pheno)


# Calculate genotype probabilities and add to file:
phd <- calc.genoprob(phd, step=2, off.end=0, error.prob=0.001,
                     map.function="kosambi", stepwidth="fixed")

# Calculate genotype probabilities and add to file:
phd_N <- calc.genoprob(phd_N, step=2, off.end=0, error.prob=0.001,
                       map.function="kosambi", stepwidth="fixed")

# Calculate genotype probabilities and add to file:
phd_S <- calc.genoprob(phd_S, step=2, off.end=0, error.prob=0.001,
                       map.function="kosambi", stepwidth="fixed")



#  Simulate genotypes (allows use of pseudomarkers in effects plot):
phd <- sim.geno(phd, n.draws=128, step=2, off.end=0, error.prob=0.001,
                map.function="kosambi", stepwidth="fixed")
geno.table(phd)

#  Simulate genotypes (allows use of pseudomarkers in effects plot):
phd_N <- sim.geno(phd_N, n.draws=128, step=2, off.end=0, error.prob=0.001,
                  map.function="kosambi", stepwidth="fixed")
geno.table(phd_N)

#  Simulate genotypes (allows use of pseudomarkers in effects plot):
phd_S <- sim.geno(phd_S, n.draws=128, step=2, off.end=0, error.prob=0.001,
                  map.function="kosambi", stepwidth="fixed")
geno.table(phd_S)


# First plot missing genotypes, then genetic maps, 
#       then plot.cross gives a summary of these data
#       reorder=TRUE ==> reorder indiviudals according to 
#       the sum of their phenotypes
plotMissing(phd)
plotMissing(phd, reorder=TRUE)
plotMap(phd, show.marker.names=TRUE, horizontal=FALSE, shift=TRUE, alternate.chrid=TRUE)

phd <- drop.nullmarkers(phd)

# Pull out genotype data as single big matrix
pull.geno(phd)

# Pull out results of calc.genoprob from a cross as a matrix
pull.genoprob(phd, omit.first.prob=FALSE, include.pos.info=FALSE, rotate=FALSE)

# Return the genetic map as a table with chromosome assignments & marker names
pull.map(phd, as.table=TRUE)

#This eliminates the margins so the plots will fit--still too many traits
#How do you plot only one column of interest in R??
op<-par(mar=rep(0,4))
plot.new()

#Reset to default plotting parameters
par(op)

# Plot genetic map of marker locations for all chromosomes
plotMap(phd, horozizontal=FALSE, shift=FALSE, show.marker.names=TRUE, alternate.chrid=TRUE)
plot(phd, auto.layout=TRUE, pheno.col="FT.Plasticity", alternate.chrid=TRUE)

# Generate vector of phenotypes:
Reflectance.Plast <- pull.pheno(phd, "Reflectance.Plasticity")
Reflectance.Cold <- pull.pheno(phd, "Reflectance.Cold")
Reflectance.Warm <- pull.pheno(phd, "Reflectance.Warm")

# Phenotypic Histograms
plot(phd, auto.layout=FALSE, pheno.col="Reflectance.Plasticity", alternate.chrid=TRUE)
plot(phd, auto.layout=FALSE, pheno.col="Reflectance.Cold", alternate.chrid=TRUE)
plot(phd, auto.layout=FALSE, pheno.col="Reflectance.Warm", alternate.chrid=TRUE)

## Reflectance ---------------------------------------------------------------
Reflectance.Plast.PermAC <- scanone(phd, pheno.col="Reflectance.Plasticity", model="normal", 
                                    method="hk", n.perm=1000)
summary(Reflectance.Plast.PermAC, thr=2, df=TRUE, format="allpeaks")
plot(Reflectance.Plast.PermAC)

Reflectance.Cold.PermAC <- scanone(phd, pheno.col="Reflectance.Cold", model="normal", 
                                   method="hk", n.perm=1000)
summary(Reflectance.Cold.PermAC, thr=2, df=TRUE, format="allpeaks")
plot(Reflectance.Cold.PermAC)

Reflectance.Warm.PermAC <- scanone(phd, pheno.col="Reflectance.Warm", model="normal", 
                                   method="hk", n.perm=1000)
summary(Reflectance.Warm.PermAC, thr=2, df=TRUE, format="allpeaks")
plot(Reflectance.Warm.PermAC)

## Reflectance

# REFLECTANCE PLASTICITY

Reflectance.Plasticity.one <- scanone(phd, pheno.col="Reflectance.Plasticity", model="normal", method="hk")
summary(Reflectance.Plasticity.one, thr=3, df=TRUE)
plot(Reflectance.Plasticity.one, ylab="Reflectance Plasticity", ylim=c(0,20))
add.threshold(Reflectance.Plasticity.one, perms=Reflectance.Plast.PermAC, alpha=0.05)

# Make some QTLs
REF.PLAST.1 <-  makeqtl(phd, 6, 18, what= "prob")
summary (REF.PLAST.1)

Recip <- as.numeric(pull.pheno(phd, "Cytoplasm"))-1

chr <- 6
pos <- "loc18"
trait <- "Reflectance.Plasticity"
covariate <- Recip

traitval <- trait
traitval <- pull.pheno(phd,trait)
covval <- covariate
if(is.character(covariate)) {covval <- pull.pheno(cross,covariate)}
geno <- data.frame(phd$geno[[chr]]$prob[,pos,])
geno$a <- geno$AC - geno$BD
geno$d <- (geno$BC + geno$AD - geno$AC - geno$BD)/2
geno$i <- geno$BC - geno$AD

testmod <- glm(traitval ~ a + d + i, data=geno, family=gaussian)
summary(testmod)

testmod2 <- glm(traitval  ~ a + d + i + covval,
                data=geno, family=gaussian)
summary(testmod2)

testmod3 <- glm(traitval  ~ (a + d + i) * covval,
                data=geno, family=gaussian)
summary(testmod3)

# Test significance of the added terms.
anova(testmod,testmod2,testmod3, test=c("Chisq"))
# ---------------------------------------------
cross <- phd
trait <- "Reflectance.Plasticity"
traitlab <- "Reflectance.Plasticity"
scaneffect <- function(cross,chr,trait,traitlab,covariate) 
{
  traitval <- trait
  if(is.character(trait)) {traitval <- pull.pheno(cross,trait)}
  covval <- covariate
  if(is.character(covariate)) {covval <- pull.pheno(cross,covariate)}
  s1.chrtrait <- scanone(cross,pheno.col=trait,chr=chr,model="normal",
                         method="hk", addcovar=covval)
  s1.chrtrait$a <- 0
  s1.chrtrait$d <- 0
  s1.chrtrait$i <- 0
  
  for(j in 1:length(cross$geno[[chr]]$prob[1,,"AC"])) {
    geno <- data.frame(cross$geno[[chr]]$prob[,j,])
    geno$a <- geno$AC - geno$BD
    geno$d <- (geno$BC + geno$AD - geno$AC - geno$BD)/2
    geno$i <- geno$BC - geno$AD
    mod <- lm(traitval ~ a + d + i,
              data=geno)
    s1.chrtrait$a[[j]] <- mod$coefficients[[2]]
    s1.chrtrait$d[[j]] <- mod$coefficients[[3]]
    s1.chrtrait$i[[j]] <- mod$coefficients[[4]]
  }
  maxtrait <- max(c(max(s1.chrtrait$lod),max(s1.chrtrait$a),max(s1.chrtrait$d),
                    max(s1.chrtrait$i)))
  mintrait <- min(c(min(s1.chrtrait$lod),min(s1.chrtrait$a),min(s1.chrtrait$d),
                    min(s1.chrtrait$i)))
  
  plot(s1.chrtrait, lodcolumn=1,col="black",
       sub=paste("Chromosome",chr,", ",traitlab),
       ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))
  
  rect(10,mintrait,20,maxtrait, lty=5, col="lightgrey")
  add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
  plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
  plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
  plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
  plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
  
}
#---------------------//\\\\-------------------------------------
#-----------------//////\\\\\\\\\\\\----------------------------------
### /////////////  EFFECTS PLOT \\\\\\\\\\\\\\\\\\\\\\\\ ### ----------------

# SUMMARY PLOT
par(mfrow=c(2,2))
# Genome Scan
plot(Reflectance.Plasticity.one, ylab="Reflectance Plasticity", ylim=c(0,20))
add.threshold(Reflectance.Plasticity.one, perms=Reflectance.Plast.PermAC, alpha=0.05)

# A,D,I Scan
plot(s1.chrtrait, lodcolumn=1,col="black",
     sub=paste("Chromosome",chr,", ",traitlab),
     ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))

rect(10,mintrait,20,maxtrait, lty=5, col="lightgrey")
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)

# PxG
mar2 <- find.marker(phd, chr=6, pos=18)
plotPXG(phd, main="TOGETHER", pheno.col="Reflectance.Plasticity", marker=mar2)

# EFFECT
effectplot(phd, main="TOGETHER", pheno.col="Reflectance.Plasticity", mname1="6@18")

#
#
#
#


REF.P_all_cyto <- fitqtl(phd, qtl=REF.PLAST.1, pheno.col="Reflectance.Plasticity",
                         covar=data.frame(covval), formula=y~Q1+covval, method="hk", get.ests=TRUE)
summary(REF.P_all_cyto)

# DROP CYTOPLASM TERM
REF.P_all_ <- fitqtl(phd, qtl=REF.PLAST.1, pheno.col="Reflectance.Plasticity",
                     covar=data.frame(covval), formula=y~Q1, method="hk", get.ests=TRUE)
summary(REF.P_all_)
#
#
##  REF.PLAST COMPLETE
###---------------------
##
#

##  REF.COLD  
###---------------------

Reflectance.Cold.one <- scanone(phd, pheno.col="Reflectance.Cold", model="normal", method="hk")
summary(Reflectance.Cold.one, thr=3, df=TRUE)
plot(Reflectance.Cold.one, ylab="Reflectance Cold", ylim=c(0,20))
add.threshold(Reflectance.Cold.one, perms=Reflectance.Cold.PermAC, alpha=0.05)

# Make some QTLs
REF.COLD.1 <-  makeqtl(phd, 6, 18, what= "prob")
summary (REF.COLD.1)

Recip <- as.numeric(pull.pheno(phd, "Cytoplasm"))-1

chr <- 6
pos <- "loc18"
trait <- "Reflectance.Cold"
covariate <- Recip

traitval <- trait
traitval <- pull.pheno(phd,trait)
covval <- covariate
if(is.character(covariate)) {covval <- pull.pheno(cross,covariate)}
geno <- data.frame(phd$geno[[chr]]$prob[,pos,])
geno$a <- geno$AC - geno$BD
geno$d <- (geno$BC + geno$AD - geno$AC - geno$BD)/2
geno$i <- geno$BC - geno$AD

testmod <- glm(traitval ~ a + d + i, data=geno, family=gaussian)
summary(testmod)

testmod2 <- glm(traitval  ~ a + d + i + covval,
                data=geno, family=gaussian)
summary(testmod2)

testmod3 <- glm(traitval  ~ (a + d + i) * covval,
                data=geno, family=gaussian)
summary(testmod3)

# Test significance of the added terms.
anova(testmod,testmod2,testmod3, test=c("Chisq"))
# ---------------------------------------------
cross <- phd
trait <- "Reflectance.Cold"
traitlab <- "Reflectance.Cold"
scaneffect <- function(cross,chr,trait,traitlab,covariate) 
{
  traitval <- trait
  if(is.character(trait)) {traitval <- pull.pheno(cross,trait)}
  covval <- covariate
  if(is.character(covariate)) {covval <- pull.pheno(cross,covariate)}
  s1.chrtrait <- scanone(cross,pheno.col=trait,chr=chr,model="normal",
                         method="hk", addcovar=covval)
  s1.chrtrait$a <- 0
  s1.chrtrait$d <- 0
  s1.chrtrait$i <- 0
  
  for(j in 1:length(cross$geno[[chr]]$prob[1,,"AC"])) {
    geno <- data.frame(cross$geno[[chr]]$prob[,j,])
    geno$a <- geno$AC - geno$BD
    geno$d <- (geno$BC + geno$AD - geno$AC - geno$BD)/2
    geno$i <- geno$BC - geno$AD
    mod <- lm(traitval ~ a + d + i,
              data=geno)
    s1.chrtrait$a[[j]] <- mod$coefficients[[2]]
    s1.chrtrait$d[[j]] <- mod$coefficients[[3]]
    s1.chrtrait$i[[j]] <- mod$coefficients[[4]]
  }
  maxtrait <- max(c(max(s1.chrtrait$lod),max(s1.chrtrait$a),max(s1.chrtrait$d),
                    max(s1.chrtrait$i)))
  mintrait <- min(c(min(s1.chrtrait$lod),min(s1.chrtrait$a),min(s1.chrtrait$d),
                    min(s1.chrtrait$i)))
  
  plot(s1.chrtrait, lodcolumn=1,col="black",
       sub=paste("Chromosome",chr,", ",traitlab),
       ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))
  
  rect(10,mintrait,20,maxtrait, lty=5, col="lightgrey")
  add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
  plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
  plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
  plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
  plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
  
}
#---------------------//\\\\-------------------------------------
#-----------------//////\\\\\\\\\\\\----------------------------------
### /////////////  EFFECTS PLOT \\\\\\\\\\\\\\\\\\\\\\\\ ### ----------------

# SUMMARY PLOT
par(mfrow=c(2,2))
# Genome Scan
plot(Reflectance.Cold.one, ylab="Reflectance Cold", ylim=c(0,20))
add.threshold(Reflectance.Cold.one, perms=Reflectance.Cold.PermAC, alpha=0.05)

# A,D,I Scan
plot(s1.chrtrait, lodcolumn=1,col="black",
     sub=paste("Chromosome",chr,", ",traitlab),
     ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))

rect(10,mintrait,20,maxtrait, lty=5, col="lightgrey")
add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05, col="black", lty=2 )
plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)

# PxG
mar2 <- find.marker(phd, chr=6, pos=18)
plotPXG(phd, main="TOGETHER", pheno.col="Reflectance.Cold", marker=mar2)

# EFFECT
effectplot(phd, main="TOGETHER", pheno.col="Reflectance.Cold", mname1="6@18")

##
#
#

REF.C_all_cyto <- fitqtl(phd, qtl=REF.COLD.1, pheno.col="Reflectance.Cold",
                         covar=data.frame(covval), formula=y~Q1+covval, method="hk", get.ests=TRUE)
summary(REF.C_all_cyto)

# DROP CYTOPLASM TERM
REF.C_all_ <- fitqtl(phd, qtl=REF.COLD.1, pheno.col="Reflectance.Cold",
                     covar=data.frame(covval), formula=y~Q1, method="hk", get.ests=TRUE)
summary(REF.C_all_)

#95% Bayesian interval
bayesint(Reflectance.Cold.one, chr=6)
# Markers flanking interval
bayesint(Reflectance.Cold.one, chr=6, expandtomarkers=TRUE)

##  REF.COLD COMPLETE
###---------------------
##
#


##  REF.WARM  
###---------------------

Reflectance.Warm.one <- scanone(phd, pheno.col="Reflectance.Warm", model="normal", method="hk")
summary(Reflectance.Warm.one, thr=3, df=TRUE)
plot(Reflectance.Warm.one, ylab="Reflectance Warm", ylim=c(0,7))
add.threshold(Reflectance.Warm.one, perms=Reflectance.Warm.PermAC, alpha=0.05)

# Make some QTLs
REF.WARM.1 <-  makeqtl(phd, 6, 46.9, what= "prob")
summary (REF.WARM.1)

Recip <- as.numeric(pull.pheno(phd, "Cytoplasm"))-1

chr <- 6
pos <- "329"
trait <- "Reflectance.Warm"
covariate <- Recip

traitval <- trait
traitval <- pull.pheno(phd,trait)
covval <- covariate
if(is.character(covariate)) {covval <- pull.pheno(cross,covariate)}
geno <- data.frame(phd$geno[[chr]]$prob[,pos,])
geno$a <- geno$AC - geno$BD
geno$d <- (geno$BC + geno$AD - geno$AC - geno$BD)/2
geno$i <- geno$BC - geno$AD

testmod <- glm(traitval ~ a + d + i, data=geno, family=gaussian)
summary(testmod)

testmod2 <- glm(traitval  ~ a + d + i + covval,
                data=geno, family=gaussian)
summary(testmod2)

testmod3 <- glm(traitval  ~ (a + d + i) * covval,
                data=geno, family=gaussian)
summary(testmod3)

# Test significance of the added terms.
anova(testmod,testmod2,testmod3, test=c("Chisq"))
# ---------------------------------------------
cross <- phd
trait <- "Reflectance.Warm"
traitlab <- "Reflectance.Warm"
scaneffect <- function(cross,chr,trait,traitlab,covariate) 
{
  traitval <- trait
  if(is.character(trait)) {traitval <- pull.pheno(cross,trait)}
  covval <- covariate
  if(is.character(covariate)) {covval <- pull.pheno(cross,covariate)}
  s1.chrtrait <- scanone(cross,pheno.col=trait,chr=chr,model="normal",
                         method="hk", addcovar=covval)
  s1.chrtrait$a <- 0
  s1.chrtrait$d <- 0
  s1.chrtrait$i <- 0
  
  for(j in 1:length(cross$geno[[chr]]$prob[1,,"AC"])) {
    geno <- data.frame(cross$geno[[chr]]$prob[,j,])
    geno$a <- geno$AC - geno$BD
    geno$d <- (geno$BC + geno$AD - geno$AC - geno$BD)/2
    geno$i <- geno$BC - geno$AD
    mod <- lm(traitval ~ a + d + i,
              data=geno)
    s1.chrtrait$a[[j]] <- mod$coefficients[[2]]
    s1.chrtrait$d[[j]] <- mod$coefficients[[3]]
    s1.chrtrait$i[[j]] <- mod$coefficients[[4]]
  }
  maxtrait <- max(c(max(s1.chrtrait$lod),max(s1.chrtrait$a),max(s1.chrtrait$d),
                    max(s1.chrtrait$i)))
  mintrait <- min(c(min(s1.chrtrait$lod),min(s1.chrtrait$a),min(s1.chrtrait$d),
                    min(s1.chrtrait$i)))
  
  plot(s1.chrtrait, lodcolumn=1,col="black",
       sub=paste("Chromosome",chr,", ",traitlab),
       ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))
  
  rect(28,mintrait,46.88,maxtrait, lty=5, col="lightgrey")
  add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
  plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
  plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
  plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
  plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
  
}
#---------------------//\\\\-------------------------------------
#-----------------//////\\\\\\\\\\\\----------------------------------
### /////////////  EFFECTS PLOT \\\\\\\\\\\\\\\\\\\\\\\\ ### ----------------

# SUMMARY PLOT
par(mfrow=c(2,2))
# Genome Scan
plot(Reflectance.Warm.one, ylab="Reflectance Warm", ylim=c(0,20))
add.threshold(Reflectance.Warm.one, perms=Reflectance.Warm.PermAC, alpha=0.05)

# A,D,I Scan
plot(s1.chrtrait, lodcolumn=1,col="black",
     sub=paste("Chromosome",chr,", ",traitlab),
     ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))

rect(28,mintrait,46.88,maxtrait, lty=5, col="lightgrey")
add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05, col="black", lty=2 )
plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)

# PxG
mar2 <- find.marker(phd, chr=6, pos=46.88)
plotPXG(phd, main="TOGETHER", pheno.col="Reflectance.Warm", marker=mar2)

# EFFECT
effectplot(phd, main="TOGETHER", pheno.col="Reflectance.Warm", mname1="6@46.88", ylim=c(85,95))

##
#
#


# FULL MODEL

REF.W_all_cyto <- fitqtl(phd, qtl=REF.WARM.1, pheno.col="Reflectance.Warm",
                         covar=data.frame(covval), formula=y~Q1+covval+Q1:covval, method="hk", get.ests=TRUE)
summary(REF.W_all_cyto)

# DROP CYTOPLASM TERM
REF.W_all_ <- fitqtl(phd, qtl=REF.WARM.1, pheno.col="Reflectance.Warm",
                     covar=data.frame(covval), formula=y~Q1, method="hk", get.ests=TRUE)
summary(REF.W_all_)

#95% Bayesian interval
bayesint(Reflectance.Warm.one, chr=6)
# Markers flanking interval
bayesint(Reflectance.Warm.one, chr=6, expandtomarkers=TRUE)
##  REF.WARM COMPLETE
###---------------------
##
#
##  REF.SCANS  
###---------------------
dev.off(dev.list()["RStudioGD"])

Reflectance.Plasticity.one <- scanone(phd, pheno.col="Reflectance.Plasticity", model="normal", method="hk")
summary(Reflectance.Plasticity.one, thr=3, df=TRUE)
plot(Reflectance.Plasticity.one, ylab="Reflectance Plasticity", ylim=c(0,20))
add.threshold(Reflectance.Plasticity.one, perms=Reflectance.Plast.PermAC, alpha=0.05)

Reflectance.Cold.one <- scanone(phd, pheno.col="Reflectance.Cold", model="normal", method="hk")
summary(Reflectance.Cold.one, thr=3, df=TRUE)
plot(Reflectance.Cold.one, ylab="Reflectance Cold", ylim=c(0,20))
add.threshold(Reflectance.Cold.one, perms=Reflectance.Cold.PermAC, alpha=0.05)

Reflectance.Warm.one <- scanone(phd, pheno.col="Reflectance.Warm", model="normal", method="hk")
summary(Reflectance.Warm.one, thr=3, df=TRUE)
plot(Reflectance.Warm.one, ylab="Reflectance Warm", ylim=c(0,7))
add.threshold(Reflectance.Warm.one, perms=Reflectance.Warm.PermAC, alpha=0.05)

Reflectance. <- plot(Reflectance.Plasticity.one, Reflectance.Cold.one, Reflectance.Warm.one,
                     ylab="lod", main="Reflectance")
add.threshold(Reflectance.Plasticity.one, perms=Reflectance.Plast.PermAC, alpha=0.05, col="black", lty=2)
add.threshold(Reflectance.Cold.one, perms=Reflectance.Cold.PermAC, alpha=0.05, col="blue", lty=3)
add.threshold(Reflectance.Warm.one, perms=Reflectance.Warm.PermAC, alpha=0.05, col="red", lty=4)

maxtrait <- max(c(max(Reflectance.Plasticity.one$lod), max(Reflectance.Cold.one$lod), max(Reflectance.Warm.one$lod)))
mintrait <- min(c(max(Reflectance.Plasticity.one$lod), min(Reflectance.Cold.one$lod), min(Reflectance.Warm.one$lod)))


Reflectance.Plast.6 <- rect(503.2, mintrait, 513.2, maxtrait, lty=1, col="black", density =30, angle=90)

Reflectance.Cold.6 <- rect(503.2, mintrait, 513.2, maxtrait, lty=2, col="blue", density =10, angle=45)

Reflectance.Warm.6 <- rect(521.2, mintrait, 540.08, maxtrait, lty=4, col="red", density =10, angle=135)

Reflectance.Plast.6 <- rect(503.2, mintrait, 513.2, maxtrait, lty=1, col="black", density =0,lwd=2, angle=90)

Reflectance.Cold.6 <- rect(503.2, mintrait, 513.2, maxtrait, lty=2, col="blue", density =0,lwd=2, angle=45)

Reflectance.Warm.6 <- rect(521.2, mintrait, 540.08, maxtrait, lty=4, col="red", density =0, lwd=2, angle=135)

plot(Reflectance.Plasticity.one, lodcolumn=1,col="black",add=TRUE)
plot(Reflectance.Cold.one, lodcolumn=1,col="blue",add=TRUE)
plot(Reflectance.Warm.one, lodcolumn=1,col="red",add=TRUE)
###################################################---------------------------
















