# Repeated measures ANOVA across the last four preshift sessions for SNCCon tunnel mice to assess
# whether they had plateaued. Using data file TunnelSNCConPreshift.csv.

install.packages('rstatix')
install.packages('ggpubr')
install.packages('tidyverse')
install.packages('psych')

TunnelSNCConPreshiftRM$Treatment <-as.factor(TunnelSNCConPreshiftRM$Treatment)
TunnelSNCConPreshiftRM$Trial <-as.factor(TunnelSNCConPreshiftRM$Trial)

# Quick visualisation of the data 
bxp <- ggboxplot(TunnelSNCConPreshiftRM, x = "Trial", y = "LCS500", add = "point")
bxp

# descriptive stats 
library(psych)
describeBy(TunnelSNCConPreshiftRM$LCS500, 
           group = TunnelSNCConPreshiftRM$Trial,
           digits= 4)

# Descriptive statistics by group 
# group: 1
# vars n  mean   sd median trimmed mad  min  max range skew kurtosis   se
# X1    1 8 22.86 7.95  20.88   22.86 5.1 11.9 36.8  24.9 0.45    -1.16 2.81
# --------------------------------------------------------------------------- 
#   group: 2
# vars n mean   sd median trimmed  mad   min  max range  skew kurtosis   se
# X1    1 8   19 4.18  20.14      19 2.25 11.33 23.5 12.17 -0.74     -1.1 1.48
# --------------------------------------------------------------------------- 
#   group: 3
# vars n  mean   sd median trimmed  mad   min  max range skew kurtosis   se
# X1    1 8 21.84 3.87  21.02   21.84 2.42 16.38 27.5 11.12 0.32    -1.41 1.37
# --------------------------------------------------------------------------- 
#   group: 4
# vars n  mean sd median trimmed  mad   min   max range skew kurtosis  se
# X1    1 8 26.56 15  21.23   26.56 5.76 16.49 62.42 45.93 1.62     1.13 5.3
# 
# # normality assumption 
library(rstatix)
TunnelSNCConPreshiftRM %>%
  group_by(Trial) %>%
  shapiro_test(LCS500)

# # A tibble: 4 x 4
# Trial variable statistic        p
# <fct> <chr>        <dbl>    <dbl>
#   1 1     LCS500       0.948 0.688   
# 2 2     LCS500       0.862 0.127   
# 3 3     LCS500       0.916 0.400   
# 4 4     LCS500       0.659 0.000758

# q q plot 
library(ggpubr)
ggqqplot(TunnelSNCConPreshiftRM, "LCS500", facet.by = "Trial")
# looks good and well distributed. 

# now to look at the stats anova 

res.aov <- anova_test(data = TunnelSNCConPreshiftRM, dv = LCS500, wid = Mouse, within = Trial)
get_anova_table(res.aov)

# ANOVA Table (type III tests)
# 
# Effect  DFn  DFd     F     p p<.05   ges
# 1  Trial 1.27 8.91 1.224 0.314       0.095