# Repeated measures ANOVA across the last four preshift sessions for SPCCon tunnel mice to assess
# whether they had plateaued. Using data file TunnelSPCConPreshiftRM.csv.

install.packages('rstatix')
install.packages('ggpubr')
install.packages('tidyverse')
install.packages('psych')

TunnelSPCConPreshiftRM$Treatment <-as.factor(TunnelSPCConPreshiftRM$Treatment)
TunnelSPCConPreshiftRM$Trial <-as.factor(TunnelSPCConPreshiftRM$Trial)

# Quick visualisation of the data 
bxp <- ggboxplot(TunnelSPCConPreshiftRM, x = "Trial", y = "LCS500", add = "point")
bxp

# descriptive stats 
library(psych)
describeBy(TunnelSPCConPreshiftRM$LCS500, 
           group = TunnelSPCConPreshiftRM$Trial,
           digits= 4)

# Descriptive statistics by group 
# group: 1
# vars n  mean    sd median trimmed  mad   min max range skew kurtosis   se
# X1    1 8 22.46 11.64  19.49   22.46 6.52 11.29  48 36.71 1.16     0.06 4.12
# --------------------------------------------------------------------------- 
#   group: 2
# vars n  mean  sd median trimmed  mad   min   max range skew kurtosis   se
# X1    1 8 22.48 3.7  22.27   22.48 3.14 16.62 27.44 10.82 0.03    -1.38 1.31
# --------------------------------------------------------------------------- 
#   group: 3
# vars n  mean    sd median trimmed  mad   min max range skew kurtosis   se
# X1    1 8 26.68 12.11  25.86   26.68 14.7 13.67  48 34.33 0.48    -1.33 4.28
# --------------------------------------------------------------------------- 
#   group: 4
# vars n  mean   sd median trimmed  mad   min   max range  skew kurtosis   se
# X1    1 8 31.92 9.47  33.05   31.92 6.29 15.77 47.33 31.56 -0.13    -0.98 3.35

# normality assumption 
library(rstatix)
TunnelSPCConPreshiftRM %>%
  group_by(Trial) %>%
  shapiro_test(LCS500)

# # A tibble: 4 x 4
# Trial variable statistic      p
# <fct> <chr>        <dbl>  <dbl>
#   1 1     LCS500       0.830 0.0592
# 2 2     LCS500       0.946 0.667 
# 3 3     LCS500       0.922 0.449 
# 4 4     LCS500       0.981 0.969 

# q q plot 
library(ggpubr)
ggqqplot(TunnelSPCConPreshiftRM, "LCS500", facet.by = "Trial")
# looks good and well distributed. 

# now to look at the stats anova 

res.aov <- anova_test(data = TunnelSPCConPreshiftRM, dv = LCS500, wid = Mouse, within = Trial)
get_anova_table(res.aov)

# ANOVA Table (type III tests)
# 
# Effect DFn DFd     F     p p<.05   ges
# 1  Trial   3  21 2.121 0.128       0.152