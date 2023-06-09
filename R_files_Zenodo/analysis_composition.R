######################################################################
###
### Analysis of community composition data
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

library(pascal)
library(vegan)

toWide <- function(d) {
    d.wide <- splt(d, by = "otu", factors = "sample", to.split = "number")
    names(d.wide) <- gsub("number.", "", names(d.wide))
    n <- d.wide$sample
    d.wide <- d.wide[, grepl("otu", names(d.wide))]
    d.wide <- as.data.frame(NAtozero(d.wide))
    rownames(d.wide) <- n
    d.wide
}

for (gene in c("pmoA", "16S")) {
    for (dy in c(0, 31, 58, 86)) {
        heading(sprintf("Gene = %s, day = %d", gene, dy), pre = 1)

        div <- read.csv(sprintf("data/%s_diversity.csv", gene),
                        stringsAsFactors = TRUE)
        div <- subset(div, day == 86)
        div$dil.exp <- -log(div$dilution, base = 10)
        div <- transform(div, dil.fac = factor(sprintf("dil-%d", dil.exp)))
        div <- subset(div, dil.fac != "dil-1")
        div <- droplevels(div)

        d <- read.csv(sprintf("data/%s_otu_long.csv", gene))

        dd <- toWide(subset(d, day == dy & d$dilution != 0.1))
        stopifnot(match(rownames(dd), div$sample) == 1:nrow(dd))

        dd.sqrtbray <- vegdist(sqrt(dd), method = "bray")

        ## Diution factor: permute plots within series

        simple.heading("Dilution", pre = 1)

        ctrl <- how(
            blocks = div$series,
            plots = Plots(strata = div$dil.fac, type = "free"),
            within = Within(type = "none"),
            nperm = 1000
        )

        dd.ado <- adonis(
            dd.sqrtbray ~ dil.fac,
            data = div,
            permutations = ctrl,
            parallel = 4
        )
        print(dd.ado)

        ## Temperature:
        ## permute temperature treatments within keeping dilution:series

        simple.heading("Temperature", pre = 1)

        div$strata <- factor(paste(div$dil.fac, div$series, sep = ":"))
        dd.ado <- adonis(
            dd.sqrtbray ~ temperature,
            data = div,
            strata = div$strata,
            permutations = 10000,
            parallel = 4
        )
        print(dd.ado)
    }
}

stop()

################################################################################
#
#  Gene = pmoA, day = 0
#

#################################  Dilution  ###################################

Call:
adonis(formula = dd.sqrtbray ~ dil.fac, data = div, permutations = ctrl)

Blocks:  div$series
Plots: div$dil.fac, plot permutation: free
Permutation: none
Number of permutations: 1000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model     R2  Pr(>F)
dil.fac    3    3.4634 1.15447  10.284 0.6067 0.01099 *
Residuals 20    2.2452 0.11226         0.3933
Total     23    5.7086                 1.0000

###############################  Temperature  ##################################


Call:
adonis(formula = dd.sqrtbray ~ temperature, data = div, permutations = 10000, strata = div$strata)

Blocks:  strata
Permutation: free
Number of permutations: 4095

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
temperature  1    0.0809 0.080851 0.31606 0.01416 0.1606
Residuals   22    5.6278 0.255807         0.98584
Total       23    5.7086                  1.00000

################################################################################
#
#  Gene = pmoA, day = 31
#


#################################  Dilution  ###################################


Call:
adonis(formula = dd.sqrtbray ~ dil.fac, data = div, permutations = ctrl)

Blocks:  div$series 
Plots: div$dil.fac, plot permutation: free
Permutation: none
Number of permutations: 1000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
dil.fac    3    2.0395 0.67984  7.6625 0.53475 0.002997 **
Residuals 20    1.7745 0.08872         0.46525
Total     23    3.8140                 1.00000

###############################  Temperature  ##################################


Call:
adonis(formula = dd.sqrtbray ~ temperature, data = div, permutations = 10000, strata = div$strata)

Blocks:  strata
Permutation: free
Number of permutations: 4095

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
temperature  1    0.0344 0.034384 0.20014 0.00902 0.9902
Residuals   22    3.7796 0.171799         0.99098
Total       23    3.8140                  1.00000

################################################################################
#
#  Gene = pmoA, day = 58
#


#################################  Dilution  ###################################


Call:
adonis(formula = dd.sqrtbray ~ dil.fac, data = div, permutations = ctrl)

Blocks:  div$series 
Plots: div$dil.fac, plot permutation: free
Permutation: none
Number of permutations: 1000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
dil.fac    3    2.1431 0.71435  9.7064 0.59283 0.002997 **
Residuals 20    1.4719 0.07360         0.40717
Total     23    3.6150                 1.00000

###############################  Temperature  ##################################


Call:
adonis(formula = dd.sqrtbray ~ temperature, data = div, permutations = 10000, strata = div$strata)

Blocks:  strata
Permutation: free
Number of permutations: 4095

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
temperature  1    0.0458 0.045785 0.28221 0.01267 0.4565
Residuals   22    3.5692 0.162235         0.98733
Total       23    3.6150                  1.00000

################################################################################
#
#  Gene = pmoA, day = 86
#


#################################  Dilution  ###################################


Call:
adonis(formula = dd.sqrtbray ~ dil.fac, data = div, permutations = ctrl, parallel = 4)

Blocks:  div$series 
Plots: div$dil.fac, plot permutation: free
Permutation: none
Number of permutations: 1000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
dil.fac    3   1.91213 0.63738  15.665 0.70147 0.001998 **
Residuals 20   0.81378 0.04069         0.29853
Total     23   2.72591                 1.00000

###############################  Temperature  ##################################


Call:
adonis(formula = dd.sqrtbray ~ temperature, data = div, permutations = 10000, strata = div$strata)

Blocks:  strata
Permutation: free
Number of permutations: 4095

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs  F.Model      R2 Pr(>F)
temperature  1   0.00765 0.007647 0.061888 0.00281  0.437
Residuals   22   2.71826 0.123557          0.99719
Total       23   2.72591                   1.00000

################################################################################
#
#  Gene = 16S, day = 0
#


#################################  Dilution  ###################################


Call:
adonis(formula = dd.sqrtbray ~ dil.fac, data = div, permutations = ctrl)

Blocks:  div$series 
Plots: div$dil.fac, plot permutation: free
Permutation: none
Number of permutations: 1000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model     R2   Pr(>F)
dil.fac    3    2.4330 0.81102  9.5342 0.5885 0.001998 **
Residuals 20    1.7013 0.08506         0.4115
Total     23    4.1343                 1.0000

###############################  Temperature  ##################################


Call:
adonis(formula = dd.sqrtbray ~ temperature, data = div, permutations = 10000, strata = div$strata)

Blocks:  strata
Permutation: free
Number of permutations: 4095

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
temperature  1    0.0316 0.031568 0.16928 0.00764 0.0957 .
Residuals   22    4.1028 0.186489         0.99236
Total       23    4.1343                  1.00000

################################################################################
#
#  Gene = 16S, day = 31
#


#################################  Dilution  ###################################


Call:
adonis(formula = dd.sqrtbray ~ dil.fac, data = div, permutations = ctrl)

Blocks:  div$series 
Plots: div$dil.fac, plot permutation: free
Permutation: none
Number of permutations: 1000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
dil.fac    3    2.3337 0.77790  6.1045 0.47799 0.001998 **
Residuals 20    2.5486 0.12743         0.52201
Total     23    4.8823                 1.00000

###############################  Temperature  ##################################


Call:
adonis(formula = dd.sqrtbray ~ temperature, data = div, permutations = 10000, strata = div$strata)

Blocks:  strata
Permutation: free
Number of permutations: 4095

Terms added sequentially (first to last)

            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
temperature  1    0.0431 0.04308 0.19585 0.00882 0.8872
Residuals   22    4.8393 0.21997         0.99118
Total       23    4.8823                 1.00000

################################################################################
#
#  Gene = 16S, day = 58
#


#################################  Dilution  ###################################


Call:
adonis(formula = dd.sqrtbray ~ dil.fac, data = div, permutations = ctrl)

Blocks:  div$series 
Plots: div$dil.fac, plot permutation: free
Permutation: none
Number of permutations: 1000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
dil.fac    3    2.4622 0.82073  7.5265 0.53029 0.007992 **
Residuals 20    2.1809 0.10905         0.46971
Total     23    4.6431                 1.00000

###############################  Temperature  ##################################


Call:
adonis(formula = dd.sqrtbray ~ temperature, data = div, permutations = 10000, strata = div$strata)

Blocks:  strata
Permutation: free
Number of permutations: 4095

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
temperature  1    0.0401 0.040132 0.19181 0.00864 0.1431
Residuals   22    4.6030 0.209227         0.99136
Total       23    4.6431                  1.00000

################################################################################
#
#  Gene = 16S, day = 86
#


#################################  Dilution  ###################################


Call:
adonis(formula = dd.sqrtbray ~ dil.fac, data = div, permutations = ctrl)

Blocks:  div$series 
Plots: div$dil.fac, plot permutation: free
Permutation: none
Number of permutations: 1000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)
dil.fac    3    2.3932 0.79772  6.8226 0.50578 0.004995 **
Residuals 20    2.3385 0.11692         0.49422
Total     23    4.7316                 1.00000

###############################  Temperature  ##################################


Call:
adonis(formula = dd.sqrtbray ~ temperature, data = div, permutations = 10000, strata = div$strata)

Blocks:  strata
Permutation: free
Number of permutations: 4095

Terms added sequentially (first to last)

            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
temperature  1    0.0232 0.023234 0.10856 0.00491 0.6089
Residuals   22    4.7084 0.214017         0.99509
Total       23    4.7316                  1.00000
