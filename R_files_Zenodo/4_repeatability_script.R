# --------------------------------------------------------------------------------
# Install and load packages

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("geomorph", "RRPP")
ipak(packages)

# global options
options(scipen = 20)


# --------------------------------------------------------------------------------
# Load data

setwd("path-to-directory") # edit file path for working directory

CVdat <- read.csv("linear_cv_data.csv", stringsAsFactors=TRUE)
LEdat <- read.csv("fox_linear_volume_data.csv", stringsAsFactors=TRUE)
str(LEdat)
str(CVdat)

#---------------------------------------------------------------------------------

##Linear Measurement repeatability
##Find replicate variation between 3 fox skulls in all linear and volumetric metrics
zsd1 <- sd(CVdat[1,2:16])
zsd2 <- sd(CVdat[2,2:16])
zsd3 <- sd(CVdat[3,2:16])
ZSD <- mean(c(zsd1,zsd2,zsd3))
cwsd1 <- sd(CVdat[1,17:31])
cwsd2 <- sd(CVdat[2,17:31])
cwsd3 <- sd(CVdat[3,17:31])
CVW <- mean(c(cwsd1,cwsd2,cwsd3))
tssd1 <- sd(CVdat[1,32:46])
tssd2 <- sd(CVdat[2,32:46])
tssd3 <- sd(CVdat[3,32:46])
TSL <- mean(c(tssd1,tssd2,tssd3))
slsd1 <- sd(CVdat[1,47:61])
slsd2 <- sd(CVdat[2,47:61])
slsd3 <- sd(CVdat[3,47:61])
SL <- mean(c(slsd1,slsd2,slsd3))
chsd1 <- sd(CVdat[1,62:76])
chsd2 <- sd(CVdat[2,62:76])
chsd3 <- sd(CVdat[3,62:76])
CVH <- mean(c(chsd1,chsd2,chsd3))
ujsd1 <- sd(CVdat[1,77:91])
ujsd2 <- sd(CVdat[2,77:91])
ujsd3 <- sd(CVdat[3,77:91])
UJW <- mean(c(ujsd1,ujsd2,ujsd3))
esd1 <- sd(CVdat[1,92:106])
esd2 <- sd(CVdat[2,92:106])
esd3 <- sd(CVdat[3,92:106])
ECV <- mean(c(esd1,esd2,esd3))

ReplicateSD <- data.frame(ZSD, CVW, TSL, SL, CVH, UJW, ECV)

##SD comparison for linear measurements
DF.sd <- subset(LEdat, population == 'Domesticated' & sex == 'F', select = c(zygomatic_width, cranial_vault_width, total_skull_length, snout_length, cranial_vault_height, upper_jaw_width, endocranial_volume))
DF.sd <- apply(DF.sd, 2, sd)

DF.sd<- t(as.data.frame(DF.sd))
Rat <- (DF.sd-ReplicateSD)/(ReplicateSD)

#-------------------------------------------------------------------------------------
##Load Data and read morphologika data into R

precision.data.23 <- read.morphologika("DM23_morphologika.txt")
dm.precision.data <- read.morphologika("DM_morphologika.txt")
precision.data.476 <- read.morphologika("DF476_morphologika.txt")
df.precision.data <- read.morphologika("DF_morphologika.txt")
precision.data.1058 <- read.morphologika("UF1058_morphologika.txt")
uf.precision.data <- read.morphologika("UF_morphologika.txt")

## Microscribe repeatability
## read morphologika data into R

## DM23
precision.data <- geomorph.data.frame(coords = precision.data.23)
PrecisionGPA.23 <- gpagen(A = precision.data$coords) ## create a GPA

## Procrustes distance among replicates
ProcrustesDistance <- PrecisionGPA.23$procD
rm1 <- mean(ProcrustesDistance)

##DM Foxes
dm.precision.data <- geomorph.data.frame(coords = dm.precision.data)
dm.PrecisionGPA <- gpagen(A = dm.precision.data$coords) ## create a GPA

## Procrustes distance among replicates
dm.ProcrustesDistance <- dm.PrecisionGPA$procD
pm1 <- mean(dm.ProcrustesDistance)

#---------------------------
## DF476
precision.data.476 <- geomorph.data.frame(coords = precision.data.476)
PrecisionGPA.476 <- gpagen(A = precision.data.476$coords) ## create a GPA

## Procrustes distance among replicates
ProcrustesDistance2 <- PrecisionGPA.476$procD
rm2 <- mean(ProcrustesDistance2)

##DF Foxes
df.precision.data <- geomorph.data.frame(coords = df.precision.data)
df.PrecisionGPA <- gpagen(A = df.precision.data$coords) ## create a GPA

## Procrustes distance among replicates
df.ProcrustesDistance <- df.PrecisionGPA$procD
pm2 <- mean(df.ProcrustesDistance)

#----------------------
## UF1058
precision.data.1058 <- geomorph.data.frame(coords = precision.data.1058)
PrecisionGPA.1058 <- gpagen(A = precision.data.1058$coords) ## create a GPA

## Procrustes distance among replicates
ProcrustesDistance3 <- PrecisionGPA.1058$procD
rm3 <- mean(ProcrustesDistance3)

##UF Foxes
uf.precision.data <- geomorph.data.frame(coords = uf.precision.data)
uf.PrecisionGPA <- gpagen(A = uf.precision.data$coords) ## create a GPA

## Procrustes distance among replicates
uf.ProcrustesDistance <- uf.PrecisionGPA$procD
pm3 <- mean(uf.ProcrustesDistance)

##Average sensitivty ratio
rat.gm <- mean(c((pm1-rm1)/rm1,(pm2-rm2)/rm2,(pm3-rm3)/rm3))




