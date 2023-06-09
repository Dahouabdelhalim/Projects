# Repeated measures ANOVA across the last four preshift sessions for SPCCon tail mice to assess
# whether they had plateaued. Using data file TailSPCConPreshiftRM.csv.

install.packages('rstatix')
install.packages('ggpubr')
install.packages('tidyverse')
install.packages('psych')

TailSPCConPreshiftRM$Treatment <-as.factor(TailSPCConPreshiftRM$Treatment)
TailSPCConPreshiftRM$Trial <-as.factor(TailSPCConPreshiftRM$Trial)

# Quick visualisation of the data 
bxp <- ggboxplot(TailSPCConPreshiftRM, x = "Trial", y = "LCS500", add = "point")
bxp

# descriptive stats 
library(psych)
describeBy(TailSPCConPreshiftRM$LCS500, 
           group = TailSPCConPreshiftRM$Trial,
           digits= 4)

# Descriptive statistics by group 
# group: 1
# vars n  mean  sd median trimmed  mad min   max range skew kurtosis   se
# X1    1 8 18.39 6.7  16.95   18.39 5.09  11 32.85 21.85 1.02    -0.03 2.37
# --------------------------------------------------------------------------- 
#   group: 2
# vars n  mean   sd median trimmed  mad min   max range skew kurtosis   se
# X1    1 8 17.87 8.28  17.41   17.87 4.03   7 35.25 28.25 0.82    -0.16 2.93
# --------------------------------------------------------------------------- 
#   group: 3
# vars n mean   sd median trimmed  mad min   max range  skew kurtosis   se
# X1    1 8 15.5 2.96  15.77    15.5 3.49  11 19.14  8.14 -0.16    -1.56 1.05
# --------------------------------------------------------------------------- 
#   group: 4
# vars n  mean    sd median trimmed mad   min   max range skew kurtosis   se
# X1    1 8 28.17 14.31  22.17   28.17 9.5 15.34 54.29 38.95 0.67    -1.28 5.06

# normality assumption 
library(rstatix)
TailSPCConPreshiftRM %>%
  group_by(Trial) %>%
  shapiro_test(LCS500)

# # A tibble: 4 x 4
# Trial variable statistic     p
# <fct> <chr>        <dbl> <dbl>
#   1 1     LCS500       0.870 0.150
# 2 2     LCS500       0.893 0.248
# 3 3     LCS500       0.942 0.632
# 4 4     LCS500       0.858 0.115

# q q plot 
library(ggpubr)
ggqqplot(TailSPCConPreshiftRM, "LCS500", facet.by = "Trial")
# looks good and well distributed. 

# now to look at the stats anova 

res.aov <- anova_test(data = TailSPCConPreshiftRM, dv = LCS500, wid = Mouse, within = Trial)
get_anova_table(res.aov)

# ANOVA Table (type III tests)
# 
# Effect  DFn   DFd     F     p p<.05   ges
# 1  Trial 1.72 12.07 2.957 0.095       0.248