######################################################################
###
### Analysis of qPCR data
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

######################################################################
###
### Notes:
###
### 1) DNA extraction was poor for the least-diluted diversity level
### (10^-1). Therefore, all analyses including DNA-analysis-based
### diversity metrics were performed excluding this diversity
### level. See manuscript for details.
###
### 2) The output of the statistical tests is copied as comment into
### this source file.

library(pascal) # https://github.com/pascal-niklaus/pascal
library(asreml) # https://vsni.co.uk/software/asreml-r

d <- read.csv("data/qPCR.csv", stringsAsFactors = TRUE)

d <- transform(d,
               fday=sprintf("day-%.0f", day),
               dil.exp = -log(dilution) / log(10))
d <- transform(d,
               fdil=sprintf("dil-%.0f", dil.exp))

d <- subset(d, dil.exp > 1)

### aggregate data by sample across days
d2 <- aggr(d,
           c("sample", "temperature", "series", "dil.exp", "fdil"),
           c("pmoA=mean(log(pmoA))", "S16=mean(log(S16))"),
           keep.numerics = TRUE)

######################################################################
###
### pmoA copy numbers

d2.asr <- asreml(pmoA ~ dil.exp + temperature,
    random = ~ series:fdil,
    data = d2
)
test.asreml(d2.asr)
## ---- Wald tests:
##             Df denDF   F.inc       Pr
## (Intercept)  1    10 2625.00 1.94e-13 ***
## dil.exp      1    10    0.64  0.44238
## temperature  1    11    4.05  0.06922 .


######################################################################
###
### 16S copy numbers

d2.asr <- asreml(S16 ~ dil.exp + temperature,
    random = ~ series:fdil,
    data = d2
)
test.asreml(d2.asr)
## ---- Wald tests:
##             Df denDF   F.inc      Pr
## (Intercept)  1    10 19910.0 < 2e-16 ***
## dil.exp      1    10     0.1 0.74694
## temperature  1    11     8.1 0.01593 *
