library(rstatix)
library(readxl);
library("dplyr") 

#figure 4 - proportion participating ACSF, gabazine, gabazine+pirenzepine; pooled ages
base <- "C:/Users/Josh Tworig/Desktop/Tworig_eLife_2021_figures+analysis/figure 4/"  #change to appropriate file location
fname <- paste(base, "Figure 4-source data 1.csv",sep="") 

fig4B_allWaves <- read.csv(fname)[1:912,1:4]
colnames(fig4B_allWaves) <- c("proportion","drug","compartment","id")

#overall anova
fig4B_allWaves %>%
  group_by(drug,compartment) %>%
  get_summary_stats(proportion, type = "mean_se")
res.aov <- anova_test(data = fig4B_allWaves, dv = proportion, between = c(drug), within = compartment, wid = id)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig4B_allWaves %>%
  group_by(drug) %>%
  t_test(
    proportion ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by drug (unpaired)
pwc <- fig4B_allWaves %>%
  group_by(compartment) %>%
  t_test(
    proportion ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc


#figure 4 - proportion participating; separate ages

fig4B_allWaves <- read.csv(fname)[1:912,5:9]
colnames(fig4B_allWaves) <- c("proportion","drug","compartment","age","id")

#overall anova
fig4B_allWaves %>%
  group_by(drug,age,compartment) %>%
  get_summary_stats(proportion, type = "mean_se")
res.aov <- anova_test(data = fig4B_allWaves, dv = proportion, between = c(drug,age), within = compartment, wid = id)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig4B_allWaves %>%
  group_by(drug,age) %>%
  t_test(
    proportion ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by drug (unpaired)
pwc <- fig4B_allWaves %>%
  group_by(compartment,age) %>%
  t_test(
    proportion ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig4B_allWaves %>%
  group_by(compartment,drug) %>%
  t_test(
    proportion ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc


#figure 4 - FC proportion participating; separate ages; acsf to gbz

fig4B_allWaves <- read.csv(fname)[1:30,10:13]

colnames(fig4B_allWaves) <- c("fc","compartment","age","fov")

#overall anova
fig4B_allWaves %>%
  group_by(compartment,age) %>%
  get_summary_stats(fc, type = "mean_se")
res.aov <- anova_test(data = fig4B_allWaves, dv = fc, between = age, within = compartment, wid = fov)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig4B_allWaves %>%
  group_by(age) %>%
  t_test(
    fc ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig4B_allWaves %>%
  group_by(compartment) %>%
  t_test(
    fc ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 4 - FC proportion participating; separate ages; gbz to gbz+pir

fig4B_allWaves <- read.csv(fname)[1:26,14:17]
colnames(fig4B_allWaves) <- c("fc","compartment","age","fov")

#overall anova
fig4B_allWaves %>%
  group_by(compartment,age) %>%
  get_summary_stats(fc, type = "mean_se")
res.aov <- anova_test(data = fig4B_allWaves, dv = fc, between = age, within = compartment, wid = fov)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig4B_allWaves %>%
  group_by(age) %>%
  t_test(
    fc ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig4B_allWaves %>%
  group_by(compartment) %>%
  t_test(
    fc ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 4 - wave-associated transient; pooled ages
fig4C_allWaves <- read.csv(fname)[1:43052,18:20]
colnames(fig4C_allWaves) <- c("transient","drug","compartment")

#overall anova
fig4C_allWaves %>%
  group_by(drug,compartment) %>%
  get_summary_stats(transient, type = "mean_se")
res.aov <- anova_test(data = fig4C_allWaves, dv = transient, between = c(drug,compartment))
get_anova_table(res.aov)

#test differences by compartment (unpaired)
pwc <- fig4C_allWaves %>%
  group_by(drug) %>%
  t_test(
    transient ~ compartment, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by drug (unpaired)
pwc <- fig4C_allWaves %>%
  group_by(compartment) %>%
  t_test(
    transient ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 4 - wave-associated transient ; separate ages
fig4C_allWaves <- read.csv(fname)[1:43052,21:24]
colnames(fig4C_allWaves) <- c("transient","drug","compartment","age")

#overall anova
fig4C_allWaves %>%
  group_by(drug,age,compartment) %>%
  get_summary_stats(transient, type = "mean_se")
res.aov <- anova_test(data = fig4C_allWaves, dv = transient, between = c(drug,compartment,age))
get_anova_table(res.aov)

#test differences by compartment (unpaired)
pwc <- fig4C_allWaves %>%
  group_by(drug,age) %>%
  t_test(
    transient ~ compartment, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by drug (unpaired)
pwc <- fig4C_allWaves %>%
  group_by(compartment,age) %>%
  t_test(
    transient ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc


#figure 4 - FC amplitude; separate ages; acsf to gbz
fig4C_allWaves <- read.csv(fname)[1:30,25:28]
colnames(fig4C_allWaves) <- c("fc","compartment","age","fov")

#overall anova
fig4C_allWaves %>%
  group_by(compartment,age) %>%
  get_summary_stats(fc, type = "mean_se")
res.aov <- anova_test(data = fig4C_allWaves, dv = fc, between = age, within = compartment, wid = fov)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig4C_allWaves %>%
  group_by(age) %>%
  t_test(
    fc ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig4C_allWaves %>%
  group_by(compartment) %>%
  t_test(
    fc ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 4 - FC proportion participating; separate ages; gbz to gbz+pir

fig4C_allWaves <- read.csv(fname)[1:26,29:32]
colnames(fig4C_allWaves) <- c("fc","compartment","age","fov")

#overall anova
fig4C_allWaves %>%
  group_by(compartment,age) %>%
  get_summary_stats(fc, type = "mean_se")
res.aov <- anova_test(data = fig4C_allWaves, dv = fc, between = age, within = compartment, wid = fov)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig4C_allWaves %>%
  group_by(age) %>%
  t_test(
    fc ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig4C_allWaves %>%
  group_by(compartment) %>%
  t_test(
    fc ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc


#figure 4 - latency comparisons
fig4D_allWaves <- read.csv(fname)[1:3354,33:34]
colnames(fig4D_allWaves) <- c("latency","drug")
#overall kruskal-wallis
fig4D_allWaves %>%
  group_by(drug) %>%
  get_summary_stats(latency, show = c("median","q1","q3"))
res <- kruskal_test(data = fig4D_allWaves, latency ~ drug)

#test differences by drug
pwc <- fig4D_allWaves %>%
  wilcox_test(
    latency ~ drug) %>%  #    latency ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#one-sample tests for difference from 0
osc <- fig4D_allWaves %>%
  group_by(drug) %>%
  wilcox_test(
    latency ~ 1, mu = 0) %>% #latency ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
osc


#figure 4 - latency comparisons
fig4D_allWaves <- read.csv(fname)[1:3662,41:43]
colnames(fig4D_allWaves) <- c("latency","drug","indicator")

#test differences by calcium indicator (unpaired)
pwc <- fig4D_allWaves %>%
  group_by(drug) %>%
  wilcox_test(
    latency ~ indicator) %>% #latency ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by calcium drug (unpaired)
pwc <- fig4D_allWaves %>%
  group_by(indicator) %>%
  wilcox_test(
    latency ~ drug) %>% #latency ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by indicator (unpaired)
pwc <- fig4D_allWaves %>%
  group_by(drug,indicator) %>%
  get_summary_stats(latency, show = c("median","q1","q3"))
pwc

osc <- fig4D_allWaves %>%
  group_by(drug,indicator) %>%
  wilcox_test(
    latency ~ 1, mu = 0) %>% #latency ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
osc


#figure 4E - ephys data
fig4E_allWaves <- read.csv(fname)[1:30,35:37]
colnames(fig4E_allWaves) <- c("iwi","drug","id")

#overall anova
fig4E_allWaves %>%
  group_by(drug) %>%
  get_summary_stats(iwi, type = "mean_se")
res.aov <- anova_test(data = fig4E_allWaves, dv = iwi, within = drug, wid = id)
get_anova_table(res.aov)

#test differences by drug (paired)
pwc <- fig4E_allWaves %>%
  t_test(
    iwi ~ drug, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 4 - ephys data - epsc
fig4E_allWaves <- read.csv(fname)[1:21,38:40]
colnames(fig4E_allWaves) <- c("epsc","drug","id")

#overall anova
fig4E_allWaves %>%
  group_by(drug) %>%
  get_summary_stats(epsc, type = "mean_se")
res.aov <- anova_test(data = fig4E_allWaves, dv = epsc, within = drug, wid = id)
get_anova_table(res.aov)

#test differences by drug (paired)
pwc <- fig4E_allWaves %>%
  t_test(
    epsc ~ drug, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc