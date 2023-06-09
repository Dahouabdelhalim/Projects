######################################################################
###
### Analysis of diversity data
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

######################################################################
###
### Notes:
###
### 1) DNA extraction was poor for the least-diluted diversity level
### (10^-1). Therefore, this level was excluded before analysis. See
### manuscript for details.
###
### 2) The output of the statistical tests is copied as comment into
### this source file.

for (gene in c("16S", "pmoA")) {
    heading(sprintf("Gene = %s", gene), center = TRUE)

    d <- read.csv(sprintf("data/%s_diversity.csv", gene),
                  stringsAsFactors = TRUE)
    d <- transform(d,
                   fday = factor(sprintf("day-%02d", day)),
                   dil.exp = -log(dilution) / log(10))
    d <- transform(d,
                   fdil = factor(sprintf("%s:%.0f", series, dil.exp)))

    for (v in c("Srare", "Hrare", "PDrare")) {
        simple.heading(sprintf("%s", v), pre = 1)
        y <- d[[v]]
        d_asr <- asreml(y ~ fday + dil.exp + temperature,
                        random = ~ fdil,
                        subset = (dil.exp > 1),
                        trace = FALSE,
                        data = d)
        test.asreml(d_asr)
    }
}

## ################################################################################
## #
## #                                  Gene = 16S
## #

## ##################################  Srare  #####################################

## ---- Wald tests:
##             Df denDF   F.inc        Pr
## (Intercept)  1  10.0 1094.00 1.446e-11 ***
## fday         3  79.1    9.40 2.205e-05 ***
## dil.exp      1  10.1  273.50 1.225e-08 ***
## temperature  1  79.1    0.09    0.7698

## ##################################  Hrare  #####################################

## ---- Wald tests:
##             Df denDF   F.inc        Pr
## (Intercept)  1    10 1782.00 1.322e-12 ***
## fday         3    79   48.17 < 2.2e-16 ***
## dil.exp      1    10  128.80 4.816e-07 ***
## temperature  1    79    0.00    0.9846

## ##################################  PDrare  ####################################

## ---- Wald tests:
##             Df denDF   F.inc        Pr
## (Intercept)  1  10.0 1308.00 5.895e-12 ***
## fday         3  79.1    9.68 1.624e-05 ***
## dil.exp      1  10.1  273.30 1.192e-08 ***
## temperature  1  79.1    0.12    0.7256

## ################################################################################
## #
## #                                 Gene = pmoA
## #

## ##################################  Srare  #####################################

## ---- Wald tests:
##             Df denDF  F.inc        Pr
## (Intercept)  1  10.1 764.50 8.020e-11 ***
## fday         3  71.4  14.31 2.133e-07 ***
## dil.exp      1  10.3 120.10 5.111e-07 ***
## temperature  1  71.1   0.84    0.3638

## ##################################  Hrare  #####################################

## ---- Wald tests:
##             Df denDF   F.inc        Pr
## (Intercept)  1  10.0 311.500 7.138e-09 ***
## fday         3  71.6   0.395    0.7568
## dil.exp      1  10.4  55.650 1.710e-05 ***
## temperature  1  71.1   0.097    0.7565

## ##################################  PDrare  ####################################

## ---- Wald tests:
##             Df denDF  F.inc        Pr
## (Intercept)  1  10.1 3254.0 5.245e-14 ***
## fday         3  71.6   30.6 7.752e-13 ***
## dil.exp      1  10.5  154.0 1.379e-07 ***
## temperature  1  71.2    2.3    0.1338
