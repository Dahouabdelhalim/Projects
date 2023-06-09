library(readxl);
library("dplyr")        

#repeated measures anova for compartment and age-dependent differences in proportion ROIs participating
library(rstatix)

base <- "C:/Users/Josh Tworig/Desktop/Tworig_eLife_2021_figures+analysis/figure 2/"  #change to appropriate file location
fname <- paste(base, "Figure 2-source data 1.csv",sep="") 

fig2C_allWaves <- read.csv(fname)[1:624,1:5]

colnames(fig2C_allWaves) <- c("prop_resp","group","age","compartment","id")
fig2C_allWaves$age[which(fig2C_allWaves$age==1)] <- "p8_9"
fig2C_allWaves$age[which(fig2C_allWaves$age==2)] <- "p11_12"
fig2C_allWaves$compartment[which(fig2C_allWaves$compartment==1)] <- "stalks"
fig2C_allWaves$compartment[which(fig2C_allWaves$compartment==2)] <- "processes"

#overall anova
fig2C_allWaves %>%
  group_by(compartment, age) %>%
  get_summary_stats(prop_resp, type = "mean_se")
res.aov <- anova_test(data = fig2C_allWaves, dv = prop_resp, wid = id, within = compartment, between = age)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig2C_allWaves %>%
  group_by(age) %>%
  t_test(
    prop_resp ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig2C_allWaves %>%
  group_by(compartment) %>%
  t_test(
    prop_resp ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

## now do same but by FOV

fig2C_allWaves <- read.csv(fname)[1:72,6:10]

colnames(fig2C_allWaves) <- c("prop_resp","group","age","compartment","id")
fig2C_allWaves$age[which(fig2C_allWaves$age==1)] <- "p8_9"
fig2C_allWaves$age[which(fig2C_allWaves$age==2)] <- "p11_12"
fig2C_allWaves$compartment[which(fig2C_allWaves$compartment==1)] <- "stalks"
fig2C_allWaves$compartment[which(fig2C_allWaves$compartment==2)] <- "processes"

#overall anova
fig2C_allWaves %>%
  group_by(compartment, age) %>%
  get_summary_stats(prop_resp, type = "mean_se")
res.aov <- anova_test(data = fig2C_allWaves, dv = prop_resp, wid = id, within = compartment, between = age)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig2C_allWaves %>%
  group_by(age) %>%
  t_test(
    prop_resp ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig2C_allWaves %>%
  group_by(compartment) %>%
  t_test(
    prop_resp ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc


#figure 2D - latency comparisons
fig2D_allWaves <- read.csv(fname)[1:1582,11:12]

colnames(fig2D_allWaves) <- c("latency","age")

fig2D_allWaves %>%
  group_by(age) %>%
  get_summary_stats(latency, show = c("median","q1","q3"))


#test differences by age (unpaired)
pwc <- fig2D_allWaves %>%
  wilcox_test(
    latency ~ age) %>%  
  adjust_pvalue(method = "none") 
pwc

#test for difference from zero
osc <- fig2D_allWaves %>%
  group_by(age) %>%
  wilcox_test(
    latency ~ 1, mu = 0)  %>%  
  adjust_pvalue(method = "BH")
osc

#figure 2D - latency comparisons by FOV
fig2D_allWaves <- read.csv(fname)[1:21,13:14]

colnames(fig2D_allWaves) <- c("latency","age")

#test differences by age (unpaired)
pwc <- fig2D_allWaves %>%
  wilcox_test(
    latency ~ age) %>%  
  adjust_pvalue(method = "none") 
pwc

#test for difference from zero
osc <- fig2D_allWaves %>%
  group_by(age) %>%
  wilcox_test(
    latency ~ 1, mu = 0)  %>%  
  adjust_pvalue(method = "BH")
osc
