######################################################################
###
### Analysis of CH4 flux data
###
### Authors: Elvira Schnyder (1), Paul L.E. Bodelier (2), Martin
### Hartmann (3), Ruth Henneberger (4), Pascal A. Niklaus (1*)
###
### (1) Department of Evolutionary Biology and Environmental Studies,
### University of Zürich, Winterthurerstrasse 190, CH-8057 Zürich
###
### (2) Department of Microbial Ecology, Netherlands Institute of
### Ecology (NIOO-KNAW), Droevendaalsesteeg 10, 6708PB Wageningen, the
### Netherlands.
###
### (3) Institute of Agricultural Sciences, Department of
### Environmental Systems Science, ETH Zürich, Universitätstrasse 2,
### 8092 Zürich, Switzerland
###
### (4) Institute of Molecular Health Science, ETH Zürich,
### Otto-Stern-Weg 7, CH-8093 Zürich
###
### * corresponding author
###
### see Schnyder et al. (20xx), Ecology, XXX, XXX for details.

library(pascal) # https://github.com/pascal-niklaus/pascal
library(asreml) # https://vsni.co.uk/software/asreml-r

d <- read.csv("data/CH4_fluxes.csv", stringsAsFactors = TRUE)

######################################################################
###
### Notes:
###
### 1) DNA extraction was poor for the least-diluted diversity level
### (10^-1). Therefore, all analyses including DNA-analysis-based
### diversity metrics were performed excluding this diversity
### level. See manuscript for details.
###
### 2) The 10^-7 dilution has substantially (ca. 8x) larger
### within-group variance than the other dilution levels. Therefore, a
### separate variance component is estimated for this group to
### fullfill the assumption of variance heterogeneity.
###
### 3) The output of the statistical tests is copied as comment into
### this source file.

######################################################################
###
### Independent variable: Dilution gradient dil.exp

d.asr <- asreml(CH4 ~ dil.exp * temperature,
                random =~ id(series):idh(fdil7):id(fdilx),
                data = d)
test.asreml(d.asr)
## ---- Wald tests:
##                     Df denDF F.inc      Pr
## (Intercept)          1  10.2  1190 7.2e-12 ***
## dil.exp              1  11.4    74 2.6e-06 ***
## temperature          1  13.0     0    0.89
## dil.exp:temperature  1  13.0     0    0.51

######################################################################
###
### Independent variable: OTU richness S

d.asr <- asreml(CH4 ~ S * temperature,
                random =~ id(series):idh(fdil7):id(fdilx),
                data = d)
test.asreml(d.asr)
## ---- Wald tests:
##               Df denDF F.inc      Pr
## (Intercept)    1  13.0 184.3 4.8e-09 ***
## S              1  17.9   1.4    0.26
## temperature    1  12.2   0.0    0.97
## S:temperature  1  12.6   0.0    0.84

## ----- excluding 10^-1 dilution because of poor DNA extraction -----

d.asr <- asreml.nvc(CH4 ~ S * temperature,
                random =~ id(series):idh(fdil7):id(fdilx),
                data = d,
                subset = dil.exp > 1)
test.asreml(d.asr)
## ---- Wald tests:
##               Df denDF F.inc      Pr
## (Intercept)    1   7.1  1702 9.2e-10 ***
## S              1   9.3    16  0.0028 **
## temperature    1   9.9     0  0.8527
## S:temperature  1  10.8     0  0.9675

######################################################################
###
### Independent variable: Shannon Index H

d.asr <- asreml(CH4 ~ H * temperature,
                random =~ id(series):idh(fdil7):id(fdilx),
                data = d)
test.asreml(d.asr)
## ---- Wald tests:
##               Df denDF F.inc      Pr
## (Intercept)    1  13.1   161 9.6e-09 ***
## H              1  16.8     0    0.97
## temperature    1  12.7     0    0.89
## H:temperature  1  12.6     0    0.88

## ----- excluding 10^-1 dilution because of poor DNA extraction -----

d.asr <- asreml.nvc(CH4 ~ H * temperature,
                random =~ id(series):idh(fdil7):id(fdilx),
                data = d,
                subset = dil.exp > 1)
test.asreml(d.asr)
## ---- Wald tests:
##               Df denDF F.inc      Pr
## (Intercept)    1   7.4  1743 5.2e-10 ***
## H              1   9.5    18   0.002 **
## temperature    1  10.0     0   0.755
## H:temperature  1  10.2     0   0.930

######################################################################
###
### Independent varibale: Phylogenetic diversity PD

d.asr <- asreml(CH4 ~ PD * temperature,
                random =~ id(series):idh(fdil7):id(fdilx),
                data = d)
test.asreml(d.asr)
## ---- Wald tests:
##                Df denDF F.inc      Pr
## (Intercept)     1  13.2 196.9 2.6e-09 ***
## PD              1  13.5   3.0    0.11
## temperature     1  12.4   0.0    0.96
## PD:temperature  1  12.8   0.5    0.50

## ----- excluding 10^-1 dilution because of poor DNA extraction -----

d.asr <- asreml.nvc(CH4 ~ PD * temperature,
                random =~ id(series):idh(fdil7):id(fdilx),
                data = d,
                subset = dil.exp > 1)
test.asreml(d.asr)
## ---- Wald tests:
##                Df denDF F.inc      Pr
## (Intercept)     1   7.1  1493 1.8e-09 ***
## PD              1  11.9    14  0.0032 **
## temperature     1   9.8     0  0.7741
## PD:temperature  1  10.6     0  0.7406
