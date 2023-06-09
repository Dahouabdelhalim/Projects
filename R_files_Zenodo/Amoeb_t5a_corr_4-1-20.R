### Amoeb and gene expression laboratory values: Correlations

# 3-2-20
# Data: 2015 Puerto Rico co-infection experiment
# note - nfkb is used to refer to IKB because it is involved in nkfb pathways 

data<-read.csv("amoeb_t5a_corr_no9.csv") # colony 9 excluded like in other lab analyses
head(data)

# None are normal
hist(data$perc_amoeb)
hist(data$t5a_ab)
hist(data$mmp_ab)
hist(data$nfkb_ab)

shapiro.test(data$perc_amoeb)
shapiro.test(data$t5a_ab)
shapiro.test(data$mmp_ab)
shapiro.test(data$nfkb_ab)

### Exploratory plots ###
plot(perc_amoeb~t5a_ab,data=data)
plot(perc_amoeb~mmp_ab,data=data)
plot(perc_amoeb~nfkb_ab,data=data)
 
# genes vs. e/o
plot(mmp_ab~t5a_ab,data=data)
plot(mmp_ab~nfkb_ab,data=data)
plot(nfkb_ab~t5a_ab,data=data)


### Spearman's rho correlation for non-parametric ###

# Amoeb vs. genes
cor.test(data$perc_amoeb,data$t5a_ab,method="spearman",cotinuity=FALSE,conf.level=0.95) # positive corr, significant 

cor.test(data$perc_amoeb,data$mmp_ab,method="spearman",cotinuity=FALSE,conf.level=0.95) # positive corr, significant 

cor.test(data$perc_amoeb,data$nfkb_ab,method="spearman",cotinuity=FALSE,conf.level=0.95) # NS

# Genes vs. each other
cor.test(data$nfkb_ab,data$t5a_ab,method="spearman",cotinuity=FALSE,conf.level=0.95) # NS

cor.test(data$nfkb_ab,data$mmp_ab,method="spearman",cotinuity=FALSE,conf.level=0.95) # NS

cor.test(data$mmp_ab,data$t5a_ab,method="spearman",cotinuity=FALSE,conf.level=0.95) # T5A and mmp correlate positively 

# Bonferroni
shapiro_p<- c(0.00038,0.007219,0.7886,0.6502,0.6261,0.0006691)
# Order: amoeb with t5a, mmp, nkfb; nkfb and t5a; nfkb and mmp; mmp and t5a
p.adjust(shapiro_p, method = "holm", n = 6) 
# Still significant: same 3 that were originally

