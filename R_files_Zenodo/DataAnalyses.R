# This script file contains all statistical analyses for the paper "Effects of
# social experience, aggressiveness and comb size on contest success in male
# domestic fowl". The results are presented in the same order as they appear in
# the paper.

# Some functions to process data

# function to z-standardize a variable
zstd <- function(x) {
    avx <- mean(x[!is.na(x)])
    sdx <- sd(x[!is.na(x)])
    (x - avx)/sdx
}

# function to get proportion wins
pWin <- function(x) if (length(x) > 0) sum(x)/length(x) else NA

# function to extract a data frame with individuals as rows
IndData <- function(dx) {
    # split data frame dx into a list of data frames according to levels of Id
    ldx <- split(dx, dx$Id)
    # take first row of data frames in list
    dxi <- ldx[[1]][1, ] # first row of first data frame in list
    # append first row of remaining data frames in list
    for (i in 2:length(levels(dxi$Id))) {
        dxi <- rbind(dxi, ldx[[i]][1, ])
    }
    dxi
}

######################################################################
# read in data, construct data frames and perform different analyses

# all data used for the analyses come from this data file (see Readme.txt file
# for description of variables)
dat0 <- read.delim("DuelData.tsv", na.strings = "", stringsAsFactors = TRUE)


##### First analysis ##################################################
# in the first analysis, we only use data for group and single (not dom or sub)
datGS <- subset(dat0, Treatment != "dom" & Treatment != "sub")
# exclude draws from the data set
datGS <- subset(datGS, !is.na(WonLost))
# see to it that the levels for Treatment and Id are right
datGS$Treatment <- droplevels(datGS$Treatment)
datGS$Id <- droplevels(datGS$Id)

# z-standardised focal - opponent differences
datGS$Aggd <- datGS$AggPreDuel - datGS$OppAgg
datGS$Aggdstd <- zstd(datGS$Aggd)
datGS$Comd <- datGS$Com - datGS$OppCom
datGS$Comdstd <- zstd(datGS$Comd)
datGS$Weid <- datGS$Wei - datGS$OppWei
datGS$Weidstd <- zstd(datGS$Weid)

# Treatment (group vs single) and focal - opponent differences in aggression,
# comb length and weight as explanatory variables
# levels of Treatment (in order): group, single

library(lme4)     # for glmer and lmer

fmTab1 <- glmer(WonLost ~ Treatment + Aggdstd + Comdstd + Weidstd +
                Batch + (1|Id), family = binomial, data = datGS)
# this analysis is reported in Table 1 in the main text
summary(fmTab1)
# Random effects:
#  Groups Name        Variance Std.Dev.
#  Id     (Intercept) 0.3015   0.5491
# Number of obs: 186, groups:  Id, 16

# Fixed effects:
#                 Estimate Std. Error z value Pr(>|z|)
# (Intercept)       1.8727     0.4612   4.060  4.9e-05 ***
# Treatmentsingle  -1.2631     0.4910  -2.573 0.010093 *
# Aggdstd           0.8871     0.2698   3.288 0.001008 **
# Comdstd           0.7279     0.2250   3.235 0.001215 **
# Weidstd          -0.2720     0.2151  -1.265 0.206003
# Batchy2014       -1.8403     0.5014  -3.670 0.000242 ***


##### Second analysis #################################################
# data for dominant - subordinate treatment comparisons; note that there are
# some "special individuals": m455 in pair p4 only has data from the beginning
# of the first round (autumn) and m461 and m479 in pair p7 changed status from
# autumn to spring.

datDS <- subset(dat0, Treatment == "dom" | Treatment == "sub")
# exclude draws from the data set
datDS <- subset(datDS, !is.na(WonLost))
datDS$Treatment <- droplevels(datDS$Treatment)
datDS$Id <- droplevels(datDS$Id)

# z-standardised focal - opponent differences
datDS$Aggd <- datDS$AggPreDuel - datDS$OppAgg
datDS$Aggdstd <- zstd(datDS$Aggd)
datDS$Comd <- datDS$Com - datDS$OppCom
datDS$Comdstd <- zstd(datDS$Comd)
datDS$Weid <- datDS$Wei - datDS$OppWei
datDS$Weidstd <- zstd(datDS$Weid)

# Treatment (group vs single) and focal - opponent differences in aggression,
# comb length and weight as explanatory variables
# levels of Treatment (in order): dom, sub

fmTab2 <- glmer(WonLost ~ Treatment + Aggdstd + Comdstd + Weidstd +
                Batch + (1|Id), family = binomial, data = datDS)
# this analysis is reported in Table 2 in the main text
summary(fmTab2)
# Random effects:
#  Groups Name        Variance Std.Dev.
#  Id     (Intercept) 1.178    1.085
# Number of obs: 116, groups:  Id, 13

# Fixed effects:
#              Estimate Std. Error z value Pr(>|z|)
# (Intercept)    3.7449     1.2125   3.088 0.002012 **
# Treatmentsub  -1.7981     1.1012  -1.633 0.102498
# Aggdstd        1.1586     0.6377   1.817 0.069245 .
# Comdstd        0.6143     0.4164   1.475 0.140098
# Weidstd       -0.2764     0.4728  -0.585 0.558833
# Batchy2014    -4.5055     1.2010  -3.751 0.000176 ***


##### Third analysis ##################################################
# Check if aggressiveness (AggMean) can be explained by personality variables
# measured in the novel arena; use the average NA measurements

# individual level data
dat0i <- IndData(dat0)
dat0i$Expstd <- zstd(dat0i$Exp1 + dat0i$Exp2)
dat0i$Vigstd <- zstd(dat0i$Vig1 + dat0i$Vig2)
dat0i$Crostd <- zstd(dat0i$Cro1 + dat0i$Cro2)

fmTab3 <- lm(AggMean ~ Expstd + Vigstd + Crostd, data = dat0i)
# this analysis is reported in Table 3 in the main text
summary(fmTab3)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  4.94483    0.19251  25.686   <2e-16 ***
# Expstd       0.28474    0.22653   1.257    0.220
# Vigstd       0.02778    0.20330   0.137    0.892
# Crostd      -0.27424    0.21997  -1.247    0.224


##### Additional analyses ##############################################

# check if social treatments (group vs. single) differ in personality traits

datGSi <- IndData(datGS)

fmGSiAgg <- lm(AggMean ~ Treatment, data = datGSi)
summary(fmGSiAgg)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
# (Intercept)       5.4375     0.3743  14.529 7.76e-10 ***
# Treatmentsingle  -0.6250     0.5293  -1.181    0.257
# Residual standard error: 1.059 on 14 degrees of freedom
# Multiple R-squared:  0.09058,   Adjusted R-squared:  0.02562
# F-statistic: 1.394 on 1 and 14 DF,  p-value: 0.2573

datGSi$ExpMean <- with(datGSi, (Exp1 + Exp2)/2)
fmGSiExp <- lm(ExpMean ~ Treatment, data = datGSi)
summary(fmGSiExp)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
# (Intercept)       0.3178     0.2821   1.127    0.279
# Treatmentsingle  -0.2309     0.3989  -0.579    0.572
# Residual standard error: 0.7979 on 14 degrees of freedom
# Multiple R-squared:  0.02337,   Adjusted R-squared:  -0.04639
# F-statistic: 0.335 on 1 and 14 DF,  p-value: 0.5719

datGSi$VigMean <- with(datGSi, (Vig1 + Vig2)/2)
fmGSiVig <- lm(VigMean ~ Treatment, data = datGSi)
summary(fmGSiVig)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
# (Intercept)      0.85250    0.02576  33.099 1.07e-14 ***
# Treatmentsingle  0.08563    0.03642   2.351   0.0339 *
# Residual standard error: 0.07285 on 14 degrees of freedom
# Multiple R-squared:  0.283, Adjusted R-squared:  0.2318
# F-statistic: 5.526 on 1 and 14 DF,  p-value: 0.03391

datGSi$CroMean <- with(datGSi, (Cro1 + Cro2)/2)
fmGSiCro <- lm(CroMean ~ Treatment, data = datGSi)
summary(fmGSiCro)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
# (Intercept)       11.562      3.128   3.697  0.00239 **
# Treatmentsingle    1.250      4.423   0.283  0.78162
# Residual standard error: 8.847 on 14 degrees of freedom
# Multiple R-squared:  0.005672,  Adjusted R-squared:  -0.06535
# F-statistic: 0.07986 on 1 and 14 DF,  p-value: 0.7816


# correlation aggression and comb size
cor.test(datGSi$AggMean, (datGSi$ComOct + datGSi$ComApr)/2, method = "spearman")
#     Spearman's rank correlation rho
# data:  datGSi$AggMean and (datGSi$ComOct + datGSi$ComApr)/2
# S = 800.6, p-value = 0.5111
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho
# -0.1773524


# investigate consistency of focal winning between autumn and spring;
# only use data for group, single and dom (not sub: no contests for autumn)
datWL <- subset(dat0, Treatment != "sub" & !is.na(WonLost))
# see to it that the levels for Treatment and Id are right
datWL$Treatment <- droplevels(datWL$Treatment)
datWL$Id <- droplevels(datWL$Id)

# add data for proportions of wins for autumn and spring
datWL$pW1 <- NA
datWL$pW2 <- NA
for (i in 1:length(levels(datWL$Id))) {
    iid <- datWL$Id == levels(datWL$Id)[i]
    iaut <- iid & datWL$Season == "aut"
    ispr <- iid & datWL$Season == "spr"
    datWL$pW1[iid] <- pWin(datWL$WonLost[iaut])
    datWL$pW2[iid] <- pWin(datWL$WonLost[ispr])
}

# won lost consistency at individual level
datWLi <- IndData(datWL)
# "long form" data set with pW1 and pW2 in the same column
datWLiL <- reshape(datWLi, dir = "long", varying = c("pW1", "pW2"),
                    idvar = "Id", sep = "")
datWLiL$time = as.factor(datWLiL$time)

library(rptR)     # for rptGaussian
reppW <- rptGaussian(pW ~ (1 | Id), grname = "Id", data = datWLiL)
summary(reppW)
# 3 rows containing missing values were removed
# Repeatability estimation using the lmm method
# Data: 45 observations
# Id (24 groups)
# Repeatability estimation overview:
#       R     SE   2.5%  97.5% P_permut  LRT_P
#   0.615   0.14  0.273  0.815       NA  0.001
# Bootstrapping and Permutation test:
#             N   Mean Median   2.5%  97.5%
# boot     1000  0.601  0.619  0.273  0.815
# permut      1     NA     NA     NA     NA
# Likelihood ratio test:
# logLik full model = -11.51168
# logLik red. model = -16.38918
# D  = 9.75, df = 1, P = 0.000894

cor.test(datWLi$pW1, datWLi$pW2, method = "spearman")
#     Spearman's rank correlation rho
# S = 390.22, p-value = 0.0001013
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho
# 0.7466109
