#statistics for Figure 3

#figure 3 - proportion participating; pooled ages

base <- "C:/Users/Josh Tworig/Desktop/Tworig_eLife_2021_figures+analysis/figure 3/"  #change to appropriate file location
fname <- paste(base, "Figure 3-source data 1.csv",sep="") 

fig3B_allWaves <- read.csv(fname)[1:936,1:4]

colnames(fig3B_allWaves) <- c("proportion","drug","compartment","id")

#overall anova
fig3B_allWaves %>%
  group_by(drug,compartment) %>%
  get_summary_stats(proportion, type = "mean_se")
res.aov <- anova_test(data = fig3B_allWaves, dv = proportion, between = c(drug), within = compartment, wid = id)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig3B_allWaves %>%
  group_by(drug) %>%
  t_test(
    proportion ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by drug (unpaired)
pwc <- fig3B_allWaves %>%
  group_by(compartment) %>%
  t_test(
    proportion ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 3 - proportion participating; separate ages
fig3B_allWaves <- read.csv(fname)[1:936,5:9]
colnames(fig3B_allWaves) <- c("proportion","drug","compartment","age","id")

#overall anova
fig3B_allWaves %>%
  group_by(drug,compartment,age) %>%
  get_summary_stats(proportion, type = "mean_se")
res.aov <- anova_test(data = fig3B_allWaves, dv = proportion, between = c(drug,age), within = compartment, wid = id)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig3B_allWaves %>%
  group_by(drug,age) %>%
  t_test(
    proportion ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by drug (unpaired)
pwc <- fig3B_allWaves %>%
  group_by(compartment,age) %>%
  t_test(
    proportion ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig3B_allWaves %>%
  group_by(compartment,drug) %>%
  t_test(
    proportion ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc


#figure 3 - FC proportion participating; separate ages

fig3B_allWaves <- read.csv(fname)[1:24,10:13]
colnames(fig3B_allWaves) <- c("fc","compartment","age","fov")

#overall anova
fig3B_allWaves %>%
  group_by(compartment,age) %>%
  get_summary_stats(fc, type = "mean_se")
res.aov <- anova_test(data = fig3B_allWaves, dv = fc, between = age, within = compartment, wid = fov)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig3B_allWaves %>%
  group_by(age) %>%
  t_test(
    fc ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig3B_allWaves %>%
  group_by(compartment) %>%
  t_test(
    fc ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 3C - latency comparisons

fig3C_allWaves <- read.csv(fname)[1:1399,14:17]
colnames(fig3C_allWaves) <- c("latency","drug","age","fov")

#overall anova
fig3C_allWaves %>%
  group_by(drug,age) %>%
  get_summary_stats(latency, show = c("median","q1","q3"))
res.aov <- anova_test(data = fig3C_allWaves, dv = latency, between = c(age,drug))
get_anova_table(res.aov)

#test differences by drug
pwc <- fig3C_allWaves %>%
  group_by(age) %>%
  wilcox_test(
    latency ~ drug) %>%  #    latency ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig3C_allWaves %>%
  group_by(drug) %>%
  wilcox_test(
    latency ~ age) %>% #latency ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

osc <- fig3C_allWaves %>%
  group_by(drug,age) %>%
  wilcox_test(
    latency ~ 1, mu = 0) %>% #latency ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
osc

#figure 3 - wave-associated transient; pooled ages
fig3D_allWaves <- read.csv(fname)[1:31377,18:20]
colnames(fig3D_allWaves) <- c("transient","drug","compartment")

#overall anova
fig3D_allWaves %>%
  group_by(drug,compartment) %>%
  get_summary_stats(transient, type = "mean_se")
res.aov <- anova_test(data = fig3D_allWaves, dv = transient, between = c(drug,compartment))
get_anova_table(res.aov)

#test differences by compartment (unpaired)
pwc <- fig3D_allWaves %>%
  group_by(drug) %>%
  t_test(
    transient ~ compartment, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by drug (unpaired)
pwc <- fig3D_allWaves %>%
  group_by(compartment) %>%
  t_test(
    transient ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 3 - wave-associated transient ; separate ages
fig3D_allWaves <- read.csv(fname)[1:31377,21:24]
colnames(fig3D_allWaves) <- c("transient","drug","compartment","age")

#overall anova
fig3D_allWaves %>%
  group_by(drug,compartment,age) %>%
  get_summary_stats(transient, type = "mean_se")
res.aov <- anova_test(data = fig3D_allWaves, dv = transient, between = c(drug,compartment,age))
get_anova_table(res.aov)

#test differences by compartment (unpaired)
pwc <- fig3D_allWaves %>%
  group_by(drug,age) %>%
  t_test(
    transient ~ compartment, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by drug (unpaired)
pwc <- fig3D_allWaves %>%
  group_by(compartment,age) %>%
  t_test(
    transient ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc


#figure 3 - FC amplitude; separate ages
fig3D_allWaves <- read.csv(fname)[1:24,25:28]

colnames(fig3D_allWaves) <- c("fc","compartment","age","fov")

#overall anova
fig3D_allWaves %>%
  group_by(compartment,age) %>%
  get_summary_stats(fc, type = "mean_se")
res.aov <- anova_test(data = fig3D_allWaves, dv = fc, between = age, within = compartment, wid = fov)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig3D_allWaves %>%
  group_by(age) %>%
  t_test(
    fc ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by age (unpaired)
pwc <- fig3D_allWaves %>%
  group_by(compartment) %>%
  t_test(
    fc ~ age, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc


