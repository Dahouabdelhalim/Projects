#the following runs statistical analyses on figure 6 source data for sholl analysis
#sholl analysis - P12 - XY

base <- "C:/..."  #change to appropriate file location
fname <- paste(base, "Figure 6-source data 1.csv",sep="") 

sholl_xy <- read.csv(fname)[1:3000,1:5]

colnames(sholl_xy) <- c("intersections","radius","geno","age","id")

#overall anova
sholl_xy %>%
  group_by(radius,geno) %>%
  get_summary_stats(intersections, type = "mean_se")
res.aov <- anova_test(data = sholl_xy, dv = intersections, between = c(geno), within = radius, wid = id)
get_anova_table(res.aov)

#test differences by genotype
pwc <- sholl_xy %>%
  group_by(radius) %>%
  t_test(
    intersections ~ geno, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 6 - sholl analysis - P30 - XY
sholl_xy <- read.csv(fname)[3001:4860,1:5]
colnames(sholl_xy) <- c("intersections","radius","geno","age","id")

#overall anova
sholl_xy %>%
  group_by(radius,geno) %>%
  get_summary_stats(intersections, type = "mean_se")
res.aov <- anova_test(data = sholl_xy, dv = intersections, between = c(geno), within = radius, wid = id)
get_anova_table(res.aov)

#test differences by genotype
pwc <- sholl_xy %>%
  group_by(radius) %>%
  t_test(
    intersections ~ geno, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 6 - sholl analysis - P12 - XZ
sholl_xz <- read.csv(fname)[1:2750,6:10]
colnames(sholl_xz) <- c("intersections","radius","geno","age","id")

#overall anova
sholl_xz %>%
  group_by(radius,geno) %>%
  get_summary_stats(intersections, type = "mean_se")
res.aov <- anova_test(data = sholl_xz, dv = intersections, between = c(geno), within = radius, wid = id)
get_anova_table(res.aov)

#test differences by genotype
pwc <- sholl_xz %>%
  group_by(radius) %>%
  t_test(
    intersections ~ geno, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 6 - sholl analysis - P30 - XZ
sholl_xz <- read.csv(fname)[2751:4455,6:10]
colnames(sholl_xz) <- c("intersections","radius","geno","age","id")

#overall anova
sholl_xz %>%
  group_by(radius,geno) %>%
  get_summary_stats(intersections, type = "mean_se")
res.aov <- anova_test(data = sholl_xz, dv = intersections, between = c(geno), within = radius, wid = id)
get_anova_table(res.aov)

#test differences by genotype
pwc <- sholl_xz %>%
  group_by(radius) %>%
  t_test(
    intersections ~ geno, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 6 - sholl analysis - P12 - YZ
sholl_yz <- read.csv(fname)[1:2750,11:15]
colnames(sholl_yz) <- c("intersections","radius","geno","age","id")

#overall anova
sholl_yz %>%
  group_by(radius,geno) %>%
  get_summary_stats(intersections, type = "mean_se")
res.aov <- anova_test(data = sholl_yz, dv = intersections, between = c(geno), within = radius, wid = id)
get_anova_table(res.aov)

#test differences by genotype
pwc <- sholl_yz %>%
  group_by(radius) %>%
  t_test(
    intersections ~ geno, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc

#figure 6 - sholl analysis - P30 - YZ
sholl_yz <- read.csv(fname)[2751:4455,11:15]
colnames(sholl_yz) <- c("intersections","radius","geno","age","id")

#overall anova
sholl_yz %>%
  group_by(radius,geno) %>%
  get_summary_stats(intersections, type = "mean_se")
res.aov <- anova_test(data = sholl_yz, dv = intersections, between = c(geno), within = radius, wid = id)
get_anova_table(res.aov)

#test differences by genotype
pwc <- sholl_yz %>%
  group_by(radius) %>%
  t_test(
    intersections ~ geno, paired = FALSE) %>%
  adjust_pvalue(method = "BH")
pwc