
#### h2 and GWAS in RepeatABEL ####  

#### PART 1: h2 analysis

library(MCMCglmm)
library(MASS) 
library(matrixcalc) 
library(Matrix)
library(magrittr)
library(corpcor)
library(lqmm)
library(data.table)

### genomic approach

# note: this code describes the final step of the analysis
# certain additional input files need to be generated first (GRM, file with phenotypic data,...)

# load the GRM ("G1a")
# GRM can be constructed in R, PLINK, GCTA,...
# most likely loaded as filename.N.bin and filename.bin

# load the phenotypic data "data"

# loading the IDs
names <- read.table("anom_523_hihi_family_data.tfam", header = FALSE)

G1a_inv <- ginv(G1a)   #remember G1a is the GRM
rownames(G1a_inv) <- names[,2]
colnames(G1a_inv) <- names[,2]

G1a_inv_step2 <- make.positive.definite(G1a_inv)
G1a_inv_final <- as(G1a_inv_step2, "dgCMatrix")

# setting up MCMCglmm prior
prior <- list(R=list(V=1,nu=1),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),
                                      G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G3=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),
                                      G4=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G5=list(V=1,nu=1,alpha.mu=0,alpha.V=1000)))

prior2 <- list(R=list(V=1,nu=1),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),
                                      G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G3=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),
                                      G4=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G5=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),
						              G6=list(V=1,nu=1,alpha.mu=0,alpha.V=1000)))

# example run MCMCglmm
Tarsus_genomic <- MCMCglmm(FldgTarsus ~ Sex + ClutchSize, 
					random=~animal+Dam.ID+Social.Sire.ID+MonthYear+Year,
					ginverse=list(animal=G1a_inv_final), data=data, prior= prior,       ## the data file contains IDs and trait columns (Sex, ...)
					nitt=503000, thin=10, burnin=3000)

## example run MCMCglmm: chromosome partitioning
Tarsus_chr_part_chr_1 <- MCMCglmm(FldgTarsus ~ Sex + ClutchSize, 
		 random=~animal+animal1+Dam.ID+Social.Sire.ID+MonthYear+Year,
		 ginverse=list(animal=G1a_inv_final, animal1=G2a_inv_final), data=data1, prior= prior2,
	     nitt=50, thin=1, burnin=3)
#		 nitt=503000, thin=10, burnin=3000)

#### run MCMCglmm: pedigree approach
Tarsus_ped <- MCMCglmm(FldgTarsus ~ Sex + ClutchSize,
                            random=~animal+Dam.ID+Social.Sire.ID+MonthYear+Year,
                            pedigree=pedigree, data=data1, prior= prior,
                            nitt=503000, thin=10, burnin=3000)
							
### PART 2: GWAS

## 
library(RepeatABEL)
# make genotype file based on transposed plink format .tped and .tmap
convert.snp.tped('anom_523_hihi_genotypes.tped', 'anom_523_hihi_family_data.tfam', 'Hihi_Date_Example_Name.out', strand = "u", bcast = 10000)

# MAKE GWAA FILE (GENABEL DATA FILE)
# combine the genabel genotype file with the phenotypic data using the "load.gwaa.data" function
# you need to specify : the name of the phenotype file and the name of the genotype file
phen_gwaa <-load.gwaa.data(phenofile = 'anom_523_hihi_phenotypes_revised.phen', 
                                     genofile = 'Hihi_Date_Example_Name.out', 
                                     force = TRUE, makemap = FALSE, sort = TRUE)

# load phenotypic data
PD <- read.table("anom_523_hihi_phenotypes_revised.phen", header = T)
PD

# SOME USEFUL FUNCTIONS
descriptives.trait(phen_gwaa) # Summary of the phenotypic data
summary <- summary(gtdata(phen_gwaa)) # Summary of the genotypic data (chromosomeID, call rates, HWE)
summary[1:5, ]
#
mean(summary$CallRate)
#

## make GRM
GRM <- compute.GRM(phen_gwaa, method = "GenABEL")

#####
Mod_1 <- preFitModel(Tarsus ~ sex + ClutchSize, 
                     random=~1|id + 1|Year + 1|Dam.ID + 1|Social.Sire.ID + 1|MonthYear,
                     genabel.data = phen_gwaa, phenotype.data = PD,
                     corStruc=list( id=list("GRM") , Year=list("Ind") , 
                                    Dam.ID=list("Ind") , Social.Sire.ID=list("Ind") , MonthYear=list("Ind")) )
GWAS_1 <- rGLS(Tarsus ~ sex + ClutchSize, genabel.data = phen_gwaa,
               phenotype.data = PD, V = Mod_1$V)                         #If V is not specified, then a model including random polygenic effects and permanent environmental effects is fitted (using the hglm package) to compute V

### following the online documentation of RepeatABEL

## extract estimated variance components for the prefitted model
est.hglm <- Mod_1$fitted.hglm
cat("Genotypic and permanent env. variance components:","\\n",
    est.hglm$varRanef,", resp.","\\n",
    "The residual variance is", est.hglm$varFix,".","\\n")

##  computations for the genotypic variance
logVCE <- est.hglm$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\\n", logVCE[1],"and",logVCE[2],", resp.", "\\n \\n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\\n")

## heritability:
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + sum(est.hglm$varRanef)),
    "\\n")

## test for inflation
lambda(GWAS_1)
estlambda(GWAS_1[, "P1df"], plot=TRUE)

