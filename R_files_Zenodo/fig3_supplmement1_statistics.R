#statistics for Figure 3-supplement 1

base <- "C:/..."  #change to appropriate file location
fname <- paste(base, "Figure 3-source data 2.csv",sep="") 

#A: proportion participating in DNQX; pooled ages
fig3s7_allWaves <- read.csv(fname)[1:88,1:4]
colnames(fig3s7_allWaves) <- c("proportion","drug","compartment","id")

#overall anova
fig3s7_allWaves %>%
  group_by(drug,compartment) %>%
  get_summary_stats(proportion, type = "mean_se")
res.aov <- anova_test(data = fig3s7_allWaves, dv = proportion, between = c(drug), within = compartment, wid = id)
get_anova_table(res.aov)

#test differences by compartment (paired)
pwc <- fig3s7_allWaves %>%
  group_by(drug) %>%
  t_test(
    proportion ~ compartment, paired = TRUE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by drug (unpaired)
pwc <- fig3s7_allWaves %>%
  group_by(compartment) %>%
  t_test(
    proportion ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure B - wave-associated transient; pooled ages
fig3s7_allWaves <- read.csv(fname)[1:1657,5:7]
colnames(fig3s7_allWaves) <- c("transient","drug","compartment")

#overall anova
fig3s7_allWaves %>%
  group_by(drug,compartment) %>%
  get_summary_stats(transient, type = "mean_se")
res.aov <- anova_test(data = fig3s7_allWaves, dv = transient, between = c(drug,compartment))
get_anova_table(res.aov)

#test differences by compartment (unpaired)
pwc <- fig3s7_allWaves %>%
  group_by(drug) %>%
  t_test(
    transient ~ compartment, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#test differences by drug (unpaired)
pwc <- fig3s7_allWaves %>%
  group_by(compartment) %>%
  t_test(
    transient ~ drug, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc