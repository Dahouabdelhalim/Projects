# Repeated measures ANOVA across the last four preshift sessions for SNCCon tail mice to assess
# whether they had plateaued. Using data file TailSNCConPreshiftRM.csv.

install.packages('rstatix')
install.packages('ggpubr')
install.packages('tidyverse')
install.packages('psych')

TailSNCConPreshiftRM$Treatment <-as.factor(TailSNCConPreshiftRM$Treatment)
TailSNCConPreshiftRM$Trial <-as.factor(TailSNCConPreshiftRM$Trial)

# Quick visualisation of the data 
bxp <- ggboxplot(TailSNCConPreshiftRM, x = "Trial", y = "LCS500", add = "point")
bxp

# descriptive stats 
library(psych)
describeBy(TailSNCConPreshiftRM$LCS500, 
           group = TailSNCConPreshiftRM$Trial,
           digits= 4)

# Descriptive statistics by group 
# group: 1
# vars n mean   sd median trimmed  mad   min   max range skew kurtosis   se
# X1    1 8 17.1 4.12   16.8    17.1 4.05 11.94 24.36 12.42 0.31    -1.23 1.46
# --------------------------------------------------------------------------- 
#   group: 2
# vars n  mean   sd median trimmed  mad   min   max range skew kurtosis   se
# X1    1 8 17.38 4.48  16.61   17.38 3.29 11.38 25.25 13.87 0.43    -1.18 1.58
# --------------------------------------------------------------------------- 
#   group: 3
# vars n  mean   sd median trimmed  mad  min  max range  skew kurtosis  se
# X1    1 8 16.28 5.36  17.15   16.28 5.43 7.38 22.9 15.52 -0.35     -1.5 1.9
# --------------------------------------------------------------------------- 
#   group: 4
# vars n  mean   sd median trimmed  mad   min   max range skew kurtosis   se
# X1    1 8 18.25 3.83  17.73   18.25 5.09 12.81 23.79 10.98 0.03    -1.61 1.35

# normality assumption 
library(rstatix)
TailSNCConPreshiftRM %>%
  group_by(Trial) %>%
  shapiro_test(LCS500)

# A tibble: 4 x 4
# Trial variable statistic     p
# <int> <chr>        <dbl> <dbl>
# 1     1 LCS500       0.962 0.825
# 2     2 LCS500       0.953 0.746
# 3     3 LCS500       0.954 0.749
# 4     4 LCS500       0.964 0.843

# q q plot 
library(ggpubr)
ggqqplot(TailSNCConPreshiftRM, "LCS500", facet.by = "Trial")
# looks good and well distributed. 

# now to look at the stats anova 

res.aov <- anova_test(data = TailSNCConPreshiftRM, dv = LCS500, wid = Mouse, within = Trial)
get_anova_table(res.aov)

# ANOVA Table (type III tests)
# 
# Effect DFn DFd     F     p p<.05   ges
# 1  Trial   3  21 0.396 0.757       0.027