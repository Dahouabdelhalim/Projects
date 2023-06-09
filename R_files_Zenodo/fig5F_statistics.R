#Figure 5F: comparing proportion stable across conditions
library(rstatix)
library(readxl);
library("dplyr") 

base <- "C:/..."  #change to appropriate file location
fname <- paste(base, "Figure 5-source data 1.csv",sep="") 

# test for differences in proportion stable in glutamate puff, P9
glu <- read.csv(fname)[1:3,2:3]
#compare overall proportions of motility categories, p9
acsf_glu_p9 <- chisq.test(glu[1:2,])
get_summary_stats(glu, show = c("median","q1","q3"))

# test for differences in proportion stable in glutamate puff, P24
glu <- read.csv(fname)[1:3,5:6]
#compare overall proportions of motility categories, p24
acsf_glu_p24 <- chisq.test(glu[1:2,])

#do pairwise tests for difference in stable

glu_props <- read.csv(fname)[1:2,8:9]
colnames(glu_props) <- c("ctrl","drug")
get_summary_stats(glu_props, show = c("median","q1","q3"))

t <- wilcox.test(glu_props$ctrl, glu_props$drug,paired = TRUE) #P9

glu_props <- read.csv(fname)[3:4,8:9]
colnames(glu_props) <- c("ctrl","drug")
get_summary_stats(glu_props, show = c("median","q1","q3"))

t <- wilcox.test(glu_props$ctrl, glu_props$drug,paired = TRUE) #P24

glu_props <- read.csv(fname)[1:4,8:9]
colnames(glu_props) <- c("ctrl","drug")

t <- wilcox.test(glu_props$ctrl, glu_props$drug,paired = TRUE) #note this combines p9 and p24 data for glutamate puff


# test for differences in proportion stable in bapta
bapta <- read.csv(fname)[1:3,11:12]
#compare overall proportions of motility categories
acsf_bapta <- chisq.test(bapta[1:2,])

#do pairwise tests for difference in stable
bapta_props <- read.csv(fname)[1:6,14:15]
colnames(bapta_props) <- c("ctrl","drug")
get_summary_stats(bapta_props, show = c("median","q1","q3"))

t <- wilcox.test(bapta_props$ctrl, bapta_props$drug,paired = FALSE)

# test for differences in proportion stable in carbachol
carb <- read.csv(fname)[1:3,17:18]
#compare overall proportions of motility categories
acsf_carb <- chisq.test(carb[1:2,])

#do pairwise test for carbachol
carb_props <- read.csv(fname)[1:16,20:21]
colnames(carb_props) <- c("ctrl","drug")
t <- wilcox.test(carb_props$ctrl, carb_props$drug,paired = TRUE)
get_summary_stats(carb_props, show = c("median","q1","q3"))

# test for differences in proportion stable in gbz,pir

gbz <- read.csv(fname)[1:3,23:24]
get_summary_stats(gbz,'median_iqr')
#compare overall proportions of motility categories
acsf_gbz <- chisq.test(gbz[1:2,1:2])

pir <- read.csv(fname)[1:3,35:36]
ctrl_pir <- chisq.test(pir[1:2,1:2])
#p_corr <- p.adjust(c(acsf_gbz$p.value, ctrl_pir$p.value), "BH")

#do pairwise test for proportion stable in gabazine, gbz+pirenzepine
gbz_pir_props <- read.csv(fname)[1:14,c(26:27,38:39)]
colnames(gbz_pir_props) <- c("ctrl_gbz","gbz","ctrl_pir","pir")
get_summary_stats(gbz_pir_props, show = c("median","q1","q3"))

t1 <- wilcox.test(gbz_pir_props$ctrl_gbz, gbz_pir_props$gbz,paired = TRUE)
t2 <- wilcox.test(gbz_pir_props$ctrl_pir, gbz_pir_props$pir,paired = TRUE)
p_corr <- p.adjust(c(t1$p.value,t2$p.value), "BH")

# test for differences in proportion stable in tboa, dnqx/ap5
acsf_tboa <- read.csv(fname)[1:3,29:30]
#compare overall proportions of motility categories
tboa <- chisq.test(acsf_tboa[1:2,1:2])

ctrl_da <- read.csv(fname)[1:3,41:42]
dnqx <- chisq.test(ctrl_da[1:2,1:2])
#p_corr <- p.adjust(c(tboa$p.value, dnqx$p.value), "BH")

#do pairwise test for tboa, dnqx/ap5
tboa_da_props <- read.csv(fname)[1:13,c(32:33,44:45)]
colnames(tboa_da_props) <- c("ctrl_tboa","tboa","ctrl_dnqx","dnqx")
get_summary_stats(tboa_da_props, show = c("median","q1","q3"))

t1 <- wilcox.test(tboa_da_props$ctrl_tboa, tboa_da_props$tboa,paired = TRUE)
t2 <- wilcox.test(tboa_da_props$ctrl_dnqx, tboa_da_props$tboa_dnqx,paired = TRUE)
p_corr <- p.adjust(c(t1$p.value,t2$p.value), "BH")
