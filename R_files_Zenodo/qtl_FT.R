# Activate Rqtl
library(qtl)
library(nortest)


rm(list=ls())
dev.off(dev.list()["RStudioGD"])

setwd("D:/PhD/Crossing/Analyses/Rqtl/P.lanceolata")
getwd()

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

# CHROMOSOMES 4 & 6 SPLIT FOR 95% CI
# phd <- read.cross(format=c("csv"), 
#                  dir="D:/PhD/Crossing/Analyses/Rqtl/P.lanceolata", 
#                  file="INFORMED_DATA_SPLIT.csv", na.strings=c("NA"), 
#                  genotypes=NULL)

getwd()
ls()
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
FT.Plast <- pull.pheno(phd, "FT.Plasticity")
FT.Cold <- pull.pheno(phd, "FT.Cold")
FT.Warm <- pull.pheno(phd, "FT.Warm")

# Phenotypic Histograms
plot(phd, auto.layout=FALSE, pheno.col="FT.Plasticity", alternate.chrid=TRUE)
plot(phd, auto.layout=FALSE, pheno.col="FT.Cold", alternate.chrid=TRUE)
plot(phd, auto.layout=FALSE, pheno.col="FT.Warm", alternate.chrid=TRUE)

# Permutation tests  ( _ df):
## FLOWERING TIME ---------------------------------------------------------------
FT.Plast.PermAC <- scanone(phd, pheno.col="FT.Plasticity", model="normal", 
                           method="hk", n.perm=1000)
summary(FT.Plast.PermAC, thr=2, df=TRUE, format="allpeaks")
plot(FT.Plast.PermAC)

FT.Cold.PermAC <- scanone(phd, pheno.col="FT.Cold", model="normal", 
                          method="hk", n.perm=1000)
summary(FT.Cold.PermAC, thr=2, df=TRUE, format="allpeaks")
plot(FT.Cold.PermAC)

FT.Warm.PermAC <- scanone(phd, pheno.col="FT.Warm", model="normal", 
                          method="hk", n.perm=1000)
summary(FT.Warm.PermAC, thr=2, df=TRUE, format="allpeaks")
plot(FT.Warm.PermAC)



##--------------QTL MAPPING --------------------------------------------------------------
## Flowering Time
dev.off(dev.list()["RStudioGD"])

FT.Plasticity.one <- scanone(phd, pheno.col="FT.Plasticity", model="normal", method="hk")
summary(FT.Plasticity.one, thr=3, df=TRUE)
FT.Plast.Scan <- plot(FT.Plasticity.one, ylab="Flowering Time Plasticity", ylim=c(0,7))
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05) 

# Make some QTLs
FT.PLAST.1 <-  makeqtl(phd, 2, 94, what= "prob")
summary (FT.PLAST.1)

# Obtain estimates of additive and dominance deviation for each locus.
# ---------------------------------------------------------------
# USE LM TO TEST a,d,i:  [NOTE:  Current version below uses glm, not lm!!]
#
# Fill in chromosome # (chr), position (pos, quoted), trait, and covariate:
# use position name, not chromosome location
#
#
##  FT.PLAST.2 QTL1
###---------------------
Recip <- as.numeric(pull.pheno(phd, "Cytoplasm"))-1

chr <- 2
pos <- "loc94"
trait <- "FT.Plasticity"
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

# Default version using fitqtl:
# FULL MODEL
all_cyto <- fitqtl(phd, qtl=FT.PLAST.1, pheno.col="FT.Plasticity",
                   covar=data.frame(covval), formula=y~Q1+covval+Q1:covval, method="hk", get.ests=TRUE)
summary(all_cyto)

all_cyto <- fitqtl(phd, qtl=FT.PLAST.1, pheno.col="FT.Plasticity",
                   covar=data.frame(covval), formula=y~Q1+covval, method="hk", get.ests=TRUE)
summary(all_cyto)

all <- fitqtl(phd, qtl=FT.PLAST.1, pheno.col="FT.Plasticity",
              formula=y~Q1, method="hk", get.ests=TRUE)
summary(all)

### \\\\\\\\\\\\\\\\\\\\\\\\  EFFECTS PLOT //////////// ### ----------------
# -------------------------------------------------------------
#-----------------\\\\\\\\\\\\//////----------------------------------
#---------------------\\\\//-------------------------------------

# ---------------------------------------------------------------

# ---------------------------------------------------------------

FUNCTION TO RUN EFFECTS SCAN:
  
  # This does something similar to effectscan in R/qtl for individual linkage groups,
  #  but handles the 4-way data in outcross F2 parameterization.  a, d, and i are plotted.
  #      I.e:  a = difference of homozygotes from mid-homozygote value
  #            d = difference of heterozygotes from mid-homozygote value
  #            i = difference of each heterozygote class from mid-heterozygote value
  # This was not used directly in the manuscript, but was used by the authors
  #  to visualize QTL effects across individual chromosomes.
  
  # NOTE:  If no covariate, remove addcovar statements from s1.chrtrait
  #  and mod lines in function!!!
cross <- phd
traitlab <- "FT.Plasticity"
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
    mod <- lm(traitval ~ a + d + i ,
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
  rect(74,mintrait,97.7,maxtrait, lty=5, col="lightgrey")
  add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
  plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
  plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
  plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
  plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
  
}
#---------------------//\\\\-------------------------------------
#-----------------//////\\\\\\\\\\\\----------------------------------
### /////////////  EFFECTS PLOT \\\\\\\\\\\\\\\\\\\\\\\\ ### ----------------
# -------------------------------------------------------------

plot(FT.Plasticity.one)
FT. <- plot(FT.Plasticity.one,
            ylab="Flowering Time")
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05, col="blue", lty=3 )
add.threshold(FT.Warm.one, perms=FT.Warm.PermAC, alpha=0.05, col="red", lty=4 )

plot(FT.Plasticity.one)

# ---------------------------------------------------------------

## EFFECTS OF QTL
max(FT.Plasticity.one)
mar <- find.marker(phd, chr=2, pos=94)
plotPXG(phd, marker=mar, pheno.col="FT.Plasticity", infer=TRUE)
plotPXG(phd, marker=mar, pheno.col="FT.Plasticity", infer=TRUE)
plotPXG(phd, marker=mar, pheno.col="FT.Plasticity", infer=TRUE)


plotPXG(phd, marker=mar, pheno.col="FT.Plasticity")
plotPXG(phd, marker=mar, pheno.col="FT.Plasticity", infer=FALSE)

#
#
par(mfrow=c(2,2))

# Genome Scan
FT.Plast.Scan <- plot(FT.Plasticity.one, ylab="Flowering Time Plasticity", ylim=c(0,7))
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05) 

# A,D,I Scan
plot(s1.chrtrait, lodcolumn=1,col="black",
     sub=paste("Chromosome",chr,", ",traitlab),
     ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))
rect(74,mintrait,97.7,maxtrait, lty=5, col="lightgrey")
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)

# PxG
plotPXG(phd, marker=mar, pheno.col="FT.Plasticity", infer=TRUE)

# Effect plot
effectplot(phd, main="TOGETHER", pheno.col="FT.Plasticity", mname1="2@94")

##  FT.PLAST COMPLETE
###----------------------------------

##
#
##  FT.COLD 
###---------------------
FT.Cold.one <- scanone(phd, pheno.col="FT.Cold", model="normal", method="hk")
summary(FT.Cold.one, thr=3, df=TRUE)
plot(FT.Cold.one, ylab="Flowering Time Cold", ylim=c(0,7))
FT.Cold.Scan <- add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05)

# Make some QTLs
FT.COLD.1 <-  makeqtl(phd, 2, 90, what= "prob")
summary (FT.COLD.1)

Recip <- as.numeric(pull.pheno(phd, "Cytoplasm"))-1

chr <- 2
pos <- "loc90"
trait <- "FT.Cold"
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
}
### \\\\\\\\\\\\\\\\\\\\\\\\  EFFECTS PLOT //////////// ### ----------------
# -------------------------------------------------------------
#-----------------\\\\\\\\\\\\//////----------------------------------
#---------------------\\\\//-------------------------------------

# ---------------------------------------------------------------

# ---------------------------------------------------------------

FUNCTION TO RUN EFFECTS SCAN:
  
  # This does something similar to effectscan in R/qtl for individual linkage groups,
  #  but handles the 4-way data in outcross F2 parameterization.  a, d, and i are plotted.
  #      I.e:  a = difference of homozygotes from mid-homozygote value
  #            d = difference of heterozygotes from mid-homozygote value
  #            i = difference of each heterozygote class from mid-heterozygote value
  # This was not used directly in the manuscript, but was used by the authors
  #  to visualize QTL effects across individual chromosomes.
  
  # NOTE:  If no covariate, remove addcovar statements from s1.chrtrait
  #  and mod lines in function!!!
cross <- phd
trait <- "FT.Cold"
traitlab <- "FT.Cold"
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
  
  rect(74,mintrait,96,maxtrait, lty=5, col="lightgrey")
  add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
  plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
  plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
  plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
  plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
  
}
#---------------------//\\\\-------------------------------------
#-----------------//////\\\\\\\\\\\\----------------------------------
### /////////////  EFFECTS PLOT \\\\\\\\\\\\\\\\\\\\\\\\ ### ----------------

## EFFECTS OF QTL
max(FT.Cold.one)
mar <- find.marker(phd, chr=2, pos=90)
plotPXG(phd, marker=mar, pheno.col="FT.Cold", infer=TRUE)
plotPXG(phd, marker=mar, pheno.col="FT.Cold", infer=TRUE)
plotPXG(phd, marker=mar, pheno.col="FT.Cold", infer=TRUE)


plotPXG(phd, marker=mar, pheno.col="FT.Cold")
plotPXG(phd, marker=mar, pheno.col="FT.Cold", infer=FALSE)

effectplot(phd, main="TOGETHER", pheno.col="FT.Cold", mname1="2@90")

# -------------------------------------------------------------
par(mfrow=c(2,2))
# Genome Scan
plot(FT.Cold.one, ylab="Flowering Time Cold", ylim=c(0,7))
FT.Cold.Scan <- add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05)
# A,D,I Scan
plot(s1.chrtrait, lodcolumn=1,col="black",
     sub=paste("Chromosome",chr,", ",traitlab),
     ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))

rect(74,mintrait,96,maxtrait, lty=5, col="lightgrey")
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)

# PxG
plotPXG(phd, marker=mar, pheno.col="FT.Cold", infer=TRUE)

# EFFECT PLOT
effectplot(phd, main="TOGETHER", pheno.col="FT.Cold", mname1="2@90")

# -end QTL 1------------------------------------------------------------------

FT.COLD.2 <-  makeqtl(phd, 4, 6, what= "prob")
summary (FT.COLD.2)

Recip <- as.numeric(pull.pheno(phd, "Cytoplasm"))-1

chr <- 4
pos <- "loc6"
trait <- "FT.Cold"
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

### \\\\\\\\\\\\\\\\\\\\\\  EFFECTS PLOT ////////// ### ----------------
cross <- phd
trait <- "FT.Cold"
traitlab <- "FT.Cold"
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
    mod <- lm(traitval ~ (a + d + i )* covval,
              data=geno)
    s1.chrtrait$a[[j]] <- mod$coefficients[[2]]
    s1.chrtrait$d[[j]] <- mod$coefficients[[3]]
    s1.chrtrait$i[[j]] <- mod$coefficients[[4]]
    s1.chrtrait$c[[j]] <- mod$coefficients[[5]]
    s1.chrtrait$ci[[j]] <- mod$coefficients[[8]]
    
  }
  maxtrait <- max(c(max(s1.chrtrait$lod),max(s1.chrtrait$a),max(s1.chrtrait$d),
                    max(s1.chrtrait$i),max(s1.chrtrait$ci)))
  mintrait <- min(c(min(s1.chrtrait$lod),min(s1.chrtrait$a),min(s1.chrtrait$d),
                    min(s1.chrtrait$i),min(s1.chrtrait$ci)))
  
  plot(s1.chrtrait, lodcolumn=1,col="black",
       sub=paste("Chromosome",chr,", ",traitlab),
       ylab="lod (black), a (blue), d (red), i (green), c:i (yellow)", ylim=c(mintrait,maxtrait))
  
  rect(0,mintrait,12,(maxtrait+.2), lty=5, col="lightgrey")
  add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05, col="black", lty=2 )
  plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
  plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
  plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
  plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
  plot(s1.chrtrait, lodcolumn=6,col="yellow", add=TRUE)
  
}
#---------------------//\\\\-------------------------------------
#-----------------//////\\\\\\\\\\\\----------------------------------
### /////////////  EFFECTS PLOT \\\\\\\\\\\\\\\\\\\\\\\\ ### ----------------
# INTERACTION PLOTS:
par(mfrow=c(2,2))

max(FT.Cold.one, chr=4)
mar2 <- find.marker(phd, chr=4, pos=6)
# BOTH N & S TOGETHER
plotPXG(phd, pheno.col="FT.Cold", marker=mar2)
plotPXG(phd, main="TOGETHER", pheno.col="FT.Cold", marker=mar2)
effectplot(phd, main="TOGETHER", pheno.col="FT.Cold", mname1="4@6", ylim=c(60,90))

# N CYTOPLASM ONLY
plotPXG(phd_N, main="N Cytoplasm", pheno.col="FT.Cold", marker=mar2)
effectplot(phd_N, main="N Cytoplasm", pheno.col="FT.Cold", mname1="251", ylim=c(60,90))

# S CYTOPLASM ONLY
plotPXG(phd_S, main="S Cytoplasm", pheno.col="FT.Cold", marker=mar2)
effectplot(phd_S, main="S Cytoplasm", pheno.col="FT.Cold", mname1="4@6", ylim=c(60,90))
par()

par(mfrow=c(2,2))
# Genome Scan
plot(FT.Cold.one, ylab="Flowering Time Cold", ylim=c(0,7))
FT.Cold.Scan <- add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05)

# A,D,I Scan
plot(s1.chrtrait, lodcolumn=1,col="black",
     sub=paste("Chromosome",chr,", ",traitlab),
     ylab="lod (black), a (blue), d (red), i (green), c:i (yellow)", ylim=c(mintrait,maxtrait))

rect(0,mintrait,12,(maxtrait+.2), lty=5, col="lightgrey")
add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05, col="black", lty=2 )
plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
plot(s1.chrtrait, lodcolumn=6,col="yellow", add=TRUE)

# PxG
plotPXG(phd, main="TOGETHER", pheno.col="FT.Cold", marker=mar2)

# EFFECT
effectplot(phd, main="TOGETHER", pheno.col="FT.Cold", mname1="4@6", ylim=c(60,90))


## INTERACTION PLOT
par(mfrow=c(1,3))

effectplot(phd, main="TOGETHER", pheno.col="FT.Cold", mname1="4@6", ylim=c(60,90))
effectplot(phd_N, main="N Cytoplasm", pheno.col="FT.Cold", mname1="4@6", ylim=c(60,90))
effectplot(phd_S, main="S Cytoplasm", pheno.col="FT.Cold", mname1="4@6", ylim=c(60,90))

print(effectplot(phd, main="TOGETHER", pheno.col="FT.Cold", mname1="4@6", ylim=c(60,90)))
print(effectplot(phd_N, main="N Cytoplasm", pheno.col="FT.Cold", mname1="4@6", ylim=c(60,90)))
print(effectplot(phd_S, main="S Cytoplasm", pheno.col="FT.Cold", mname1="4@6", ylim=c(60,90)))

plotPXG(phd, main="TOGETHER", pheno.col="FT.Cold", marker=mar2, ylim=c(40,125))
plotPXG(phd_N, main="N Cytoplasm", pheno.col="FT.Cold", marker=mar2, ylim=c(40,125))
plotPXG(phd_S, main="S Cytoplasm", pheno.col="FT.Cold", marker=mar2, ylim=c(40,125))

# - END QTL 2 ------------------------------------------------------------
# ---------------------------------------------
# ---------------------------------------------
FT.COLD.3 <-  makeqtl(phd, 6, 38.6, what= "prob")
summary (FT.COLD.3)

Recip <- as.numeric(pull.pheno(phd, "Cytoplasm"))-1

chr <- 6
pos <- "32"
trait <- "FT.Cold"
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
#--------------------------------------
cross <- phd
trait <- "FT.Cold"
traitlab <- "FT.Cold"
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
  
  rect(0,mintrait,46.88,maxtrait, lty=5, col="lightgrey")
  add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
  plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
  plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
  plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
  plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
  
}
#---------------------//\\\\-------------------------------------
#-----------------//////\\\\\\\\\\\\----------------------------------
### /////////////  EFFECTS PLOT \\\\\\\\\\\\\\\\\\\\\\\\ ### ----------------
  
max(FT.Cold.one, chr=6)
mar2 <- find.marker(phd, chr=6, pos=38.61)

# SUMMARY PLOT
par(mfrow=c(2,2))
# Genome Scan
plot(FT.Cold.one, ylab="Flowering Time Cold", ylim=c(0,7))
FT.Cold.Scan <- add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05)

# A,D,I Scan
plot(s1.chrtrait, lodcolumn=1,col="black",
     sub=paste("Chromosome",chr,", ",traitlab),
     ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))

rect(0,mintrait,46.88,maxtrait, lty=5, col="lightgrey")
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)


# PxG
plotPXG(phd, main="TOGETHER", pheno.col="FT.Cold", marker=mar2)

# EFFECT
effectplot(phd, main="TOGETHER", pheno.col="FT.Cold", mname1="6@38.61", ylim=c(60,90))



# MAKE MULTI-QTL PARAMETER
FT.COLD.1.2.3 <- makeqtl(phd, c(2,4,6), c(90,6,38.6), what= "prob")
summary (FT.COLD.1.2.3)

# Default version using fitqtl:
#FULL MODEL
all_it_cyto <- fitqtl(phd, qtl=FT.COLD.1.2.3, pheno.col="FT.Cold",
                      covar=data.frame(covval), formula=y~Q1+Q2+Q3+Q1:Q2+Q1:Q3+Q2:Q3+covval, method="hk", get.ests=TRUE)
summary(all_it_cyto)


all_int_cyto <- fitqtl(phd, qtl=FT.COLD.1.2.3, pheno.col="FT.Cold",
                       covar=data.frame(covval), formula=y~Q1+Q2+Q3+Q1:Q2+Q1:Q3+Q2:Q3+covval, method="hk", get.ests=TRUE)
summary(all_int_cyto)

# DROP INTERACTION TERMS
all_cyto <- fitqtl(phd, qtl=FT.COLD.1.2.3, pheno.col="FT.Cold",
                   covar=data.frame(covval), formula=y~Q1+Q2+Q3+covval, method="hk", get.ests=TRUE)
summary(all_cyto)

# DROP CYTOPLASM TERM
all_ <- fitqtl(phd, qtl=FT.COLD.1.2.3, pheno.col="FT.Cold",
               formula=y~Q1+Q2+Q3, method="hk", get.ests=TRUE)
summary(all_)




# BEST MODEL flowering time cold -->  
all_cyto <- fitqtl(phd, qtl=FT.COLD.1.2.3, pheno.col="FT.Cold",
                   covar=data.frame(covval), formula=y~Q1+Q2+Q3+covval+Q2:covval, method="hk", get.ests=TRUE)
summary(all_cyto)


#
#
##  FT.COLD COMPLETE
###---------------------
##
#
##  FT.WARM 
###---------------------
#
FT.Warm.one <- scanone(phd, pheno.col="FT.Warm", model="normal", method="hk")
summary(FT.Warm.one, thr=3, df=TRUE)
plot(FT.Warm.one, ylab="Flowering Time Warm", ylim=c(0,9))
FT.Warm.Scan <- add.threshold(FT.Warm.one, perms=FT.Warm.PermAC, alpha=0.05)

# Make some QTLs
FT.WARM.1 <-  makeqtl(phd, 4, 2, what= "prob")
summary (FT.WARM.1)

Recip <- as.numeric(pull.pheno(phd, "Cytoplasm"))-1

chr <- 4
pos <- "loc2"
trait <- "FT.Warm"
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
trait <- "FT.Warm"
traitlab <- "FT.Warm"
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
  
  rect(0,mintrait,12,maxtrait, lty=5, col="lightgrey")
  rect(39,mintrait,54.13,maxtrait, lty=5, col="lightgrey")
  add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
  plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
  plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
  plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
  plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
  
}
#---------------------//\\\\-------------------------------------
#-----------------//////\\\\\\\\\\\\----------------------------------
### /////////////  EFFECTS PLOT \\\\\\\\\\\\\\\\\\\\\\\\ ### ----------------



mar2 <- find.marker(phd, chr=4, pos=2)

# SUMMARY PLOT
par(mfrow=c(2,2))
# Genome Scan
plot(FT.Warm.one, ylab="Flowering Time Warm", ylim=c(0,9))
FT.Warm.Scan <- add.threshold(FT.Warm.one, perms=FT.Warm.PermAC, alpha=0.05)

# A,D,I Scan
plot(s1.chrtrait, lodcolumn=1,col="black",
     sub=paste("Chromosome",chr,", ",traitlab),
     ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))

rect(0,mintrait,12,maxtrait, lty=5, col="lightgrey")
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)

# PxG
plotPXG(phd, main="TOGETHER", pheno.col="FT.Warm", marker=mar2)

# EFFECT
par(mfrow=c(1,3))

effectplot(phd, main="TOGETHER", pheno.col="FT.Warm", mname1="4@2", ylim=c(0,60))
effectplot(phd_N, main="North", pheno.col="FT.Warm", mname1="4@2", ylim=c(0,60))
effectplot(phd_S, main="South", pheno.col="FT.Warm", mname1="4@2", ylim=c(0,60))

print(effectplot(phd, main="TOGETHER", pheno.col="FT.Warm", mname1="4@2", ylim=c(0,60)))
print(effectplot(phd_N, main="North", pheno.col="FT.Warm", mname1="4@2", ylim=c(0,60)))
print(effectplot(phd_S, main="South", pheno.col="FT.Warm", mname1="4@2", ylim=c(0,60)))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FT.WARM.2 <-  makeqtl(phd, 4, 473, what= "prob")
summary (FT.WARM.2)

Recip <- as.numeric(pull.pheno(phd, "Cytoplasm"))-1

chr <- 4
pos <- "473"
trait <- "FT.Warm"
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
trait <- "FT.Warm"
traitlab <- "FT.Warm"
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
  
  rect(42.83,mintrait,54.13,maxtrait, lty=5, col="lightgrey")
  add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
  plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
  plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
  plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
  plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
  
}
#---------------------//\\\\-------------------------------------
#-----------------//////\\\\\\\\\\\\----------------------------------
### /////////////  EFFECTS PLOT \\\\\\\\\\\\\\\\\\\\\\\\ ### ----------------

mar2 <- find.marker(phd, chr=4, pos=54.13)

# SUMMARY PLOT
par(mfrow=c(2,2))
# Genome Scan
plot(FT.Warm.one, ylab="Flowering Time Warm", ylim=c(0,9))
FT.Warm.Scan <- add.threshold(FT.Warm.one, perms=FT.Warm.PermAC, alpha=0.05)

# A,D,I Scan
plot(s1.chrtrait, lodcolumn=1,col="black",
     sub=paste("Chromosome",chr,", ",traitlab),
     ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))

rect(42.83,mintrait,54.13,maxtrait, lty=5, col="lightgrey")
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)

# PxG
plotPXG(phd, main="TOGETHER", pheno.col="FT.Warm", marker=mar2)

# EFFECT
effectplot(phd, main="TOGETHER", pheno.col="FT.Warm", mname1="4@54.13", ylim=c(0,60))


# LOC 6 CHROMOSOME 4 IN WARM \\/\\/\\/\\/\\/
FT.WARM.2 <-  makeqtl(phd, 4, 473, what= "prob")
summary (FT.WARM.2)

Recip <- as.numeric(pull.pheno(phd, "Cytoplasm"))-1

chr <- 4
pos <- "loc6"
trait <- "FT.Warm"
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
trait <- "FT.Warm"
traitlab <- "FT.Warm"
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
  
  rect(42.83,mintrait,54.13,maxtrait, lty=5, col="lightgrey")
  add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
  plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
  plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
  plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
  plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)
  
}
#---------------------//\\\\-------------------------------------
#-----------------//////\\\\\\\\\\\\----------------------------------
### /////////////  EFFECTS PLOT \\\\\\\\\\\\\\\\\\\\\\\\ ### ----------------

mar2 <- find.marker(phd, chr=4, pos=54.13)

# SUMMARY PLOT
par(mfrow=c(2,2))
# Genome Scan
plot(FT.Warm.one, ylab="Flowering Time Warm", ylim=c(0,9))
FT.Warm.Scan <- add.threshold(FT.Warm.one, perms=FT.Warm.PermAC, alpha=0.05)

# A,D,I Scan
plot(s1.chrtrait, lodcolumn=1,col="black",
     sub=paste("Chromosome",chr,", ",traitlab),
     ylab="lod (black), a (blue), d (red), i (green)", ylim=c(mintrait,maxtrait))

rect(42.83,mintrait,54.13,maxtrait, lty=5, col="lightgrey")
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
plot(s1.chrtrait, lodcolumn=2,col="blue",add=TRUE)
plot(s1.chrtrait, lodcolumn=3,col="red",add=TRUE)
plot(s1.chrtrait, lodcolumn=4,col="green",add=TRUE)
plot(s1.chrtrait, lodcolumn=1,col="black", add=TRUE)

# PxG
plotPXG(phd, main="TOGETHER", pheno.col="FT.Warm", marker=mar2)

# EFFECT
par(mfrow=c(1,3))

effectplot(phd, main="TOGETHER", pheno.col="FT.Warm", mname1="4@6", ylim=c(17,42))
effectplot(phd_N, main="North", pheno.col="FT.Warm", mname1="4@6", ylim=c(17,42))
effectplot(phd_S, main="South", pheno.col="FT.Warm", mname1="4@6", ylim=c(17,42))

print(effectplot(phd, main="TOGETHER", pheno.col="FT.Warm", mname1="4@6", ylim=c(17,42)))
print(effectplot(phd_N, main="North", pheno.col="FT.Warm", mname1="4@6", ylim=c(17,42)))
print(effectplot(phd_S, main="South", pheno.col="FT.Warm", mname1="4@6", ylim=c(17,42)))

# MAKE MULTI-QTL PARAMETER
FT.WARM.1.2 <- makeqtl(phd, c(4,4), c(2,54.13), what= "prob")
summary (FT.WARM.1.2)

# Default version using fitqtl:
FTW_all_cyto <- fitqtl(phd, qtl=FT.WARM.1.2, pheno.col="FT.Warm",
                       covar=data.frame(covval), formula=y~Q1+Q2+Q1:Q2+covval, method="hk", get.ests=TRUE)
summary(FTW_all_cyto)

# DROP INTERACTION TERM
FTW_all_cyto <- fitqtl(phd, qtl=FT.WARM.1.2, pheno.col="FT.Warm",
                       covar=data.frame(covval), formula=y~Q1+Q2+covval, method="hk", get.ests=TRUE)
summary(FTW_all_cyto)

# DROP CYTOPLASM TERM
FTW_all_ <- fitqtl(phd, qtl=FT.WARM.1.2, pheno.col="FT.Warm",
                   formula=y~Q1+Q2, method="hk", get.ests=TRUE)
summary(FTW_all_)

# Default version using fitqtl, ONE QTL ONLY
FTW_ONE_cyto <- fitqtl(phd, qtl=FT.WARM.1.2, pheno.col="FT.Warm",
                       covar=data.frame(covval), formula=y~Q1+covval, method="hk", get.ests=TRUE)
summary(FTW_ONE_cyto)

# DROP INTERACTION TERM
FTW_ONE_ <- fitqtl(phd, qtl=FT.WARM.1.2, pheno.col="FT.Warm",
                   covar=data.frame(covval), formula=y~Q1, method="hk", get.ests=TRUE)
summary(FTW_ONE_)

# DROP CYTOPLASM TERM
FTW_all_ <- fitqtl(phd, qtl=FT.WARM.1.2, pheno.col="FT.Warm",
                   formula=y~Q1+Q2, method="hk", get.ests=TRUE)
summary(FTW_all_)
#
#

##  FT.WARM COMPLETE
###---------------------
##
#
##  FT.QTL SCANS 
###---------------------
dev.off(dev.list()["RStudioGD"])

FT.Plasticity.one <- scanone(phd, pheno.col="FT.Plasticity", model="normal", method="hk")
summary(FT.Plasticity.one, thr=3, df=TRUE)
FT.Plast.Scan <- plot(FT.Plasticity.one, ylab="Flowering Time Plasticity", ylim=c(0,7))
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05) 

FT.Cold.one <- scanone(phd, pheno.col="FT.Cold", model="normal", method="hk")
summary(FT.Cold.one, thr=3, df=TRUE)
plot(FT.Cold.one, ylab="Flowering Time Cold", ylim=c(0,7))
FT.Cold.Scan <- add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05)

FT.Warm.one <- scanone(phd, pheno.col="FT.Warm", model="normal", method="hk")
summary(FT.Warm.one, thr=3, df=TRUE)
plot(FT.Warm.one, ylab="Flowering Time Warm", ylim=c(0,9))
FT.Warm.Scan <- add.threshold(FT.Warm.one, perms=FT.Warm.PermAC, alpha=0.05)

FT. <- plot(FT.Plasticity.one, FT.Cold.one, FT.Warm.one,
            ylab="lod", main="Flowering Time")
add.threshold(FT.Plasticity.one, perms=FT.Plast.PermAC, alpha=0.05, col="black", lty=2 )
add.threshold(FT.Cold.one, perms=FT.Cold.PermAC, alpha=0.05, col="blue", lty=3)
add.threshold(FT.Warm.one, perms=FT.Warm.PermAC, alpha=0.05, col="red", lty=4)

maxtrait <- max(c(max(FT.Plasticity.one$lod), max(FT.Cold.one$lod), max(FT.Warm.one$lod)))
mintrait <- min(c(max(FT.Plasticity.one$lod), min(FT.Cold.one$lod), min(FT.Warm.one$lod)))
                                

FT.Plast.2 <- rect(204.62, mintrait, 228.32, maxtrait, lty=1,  col="black", density =30, angle=90)

FT.Cold.2 <- rect(204.62, mintrait, 226.62, maxtrait, lty=2,  col="blue", density =10, angle=45)
FT.Cold.4 <- rect(341.11, mintrait, 353.11, maxtrait, lty=2,  col="blue", density =10, angle=45)
FT.Cold.6 <- rect(493.2, mintrait, 540.08, maxtrait, lty=2,  col="blue", density =10, angle=45)
  
FT.Warm.4 <- rect(341.11, mintrait, 353.11, maxtrait, lty=4, col="red", density =10, angle=135)
  
FT.Plast.2 <- rect(204.62, mintrait, 228.32, maxtrait, lty=1,  col="black", density =0, lwd=2, angle=90)

FT.Cold.2 <- rect(204.62, mintrait, 226.62, maxtrait, lty=2,  col="blue", density =0, lwd=2, angle=45)
FT.Cold.4 <- rect(341.11, mintrait, 353.11, maxtrait, lty=2,  col="blue", density =0, lwd=2, angle=45)
FT.Cold.6 <- rect(493.2, mintrait, 540.08, maxtrait, lty=2,  col="blue", density =0, lwd=2, angle=45)

FT.Warm.4 <- rect(341.11, mintrait, 353.11, maxtrait, lty=4, col="red", density =0, lwd=2, angle=135)

plot(FT.Plasticity.one, lodcolumn=1,col="black",add=TRUE)
plot(FT.Cold.one, lodcolumn=1,col="blue",add=TRUE)
plot(FT.Warm.one, lodcolumn=1,col="red",add=TRUE)



## Interval estimates of QTL positions
# 1.5-LOD support interval
# lodint(FT.Cold.one, chr=6)
#95% Bayesian interval
bayesint(FT.Plasticity.one, chr=2)
bayesint(FT.Cold.one, chr=2)
bayesint(FT.Cold.one, chr=4)
bayesint(FT.Cold.one, chr=6)
bayesint(FT.Warm.one, chr=4)
# Markers flanking interval
# lodint(FT.Cold.one, chr=6, expandtomarkers=TRUE)
bayesint(FT.Plasticity.one, chr=2, expandtomarkers=TRUE)
bayesint(FT.Cold.one, chr=2, expandtomarkers=TRUE)
bayesint(FT.Cold.one, chr=4, expandtomarkers=TRUE)
bayesint(FT.Cold.one, chr=6, expandtomarkers=TRUE)
bayesint(FT.Warm.one, chr=4, expandtomarkers=TRUE)

## - - - - CONFIDENCE INTERVALS--SPLIT CHROMOSOME 4 INTO {4 AND 5}
#95% Bayesian interval
bayesint(FT.Warm.one, chr=4)
bayesint(FT.Warm.one, chr=5)
# lodint(FT.Cold.one, chr=6, expandtomarkers=TRUE)
bayesint(FT.Warm.one, chr=4, expandtomarkers=TRUE)
bayesint(FT.Warm.one, chr=5, expandtomarkers=TRUE)


# 2.0-LOD support interval
lodint(out.hk, chr=7, drop=2)
# 99% Bayes credible interval
bayesint(FT.Cold.one, chr=6, prob=0.99)
###################################################---------------------------
##----------------------------------------------------------------------------