# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

####################################################################
## Microbiome similarity between individual sifaka
####################################################################
rm(list=ls())
graphics.off()

## set your working directory to where files are located
## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
## set your working directory to where files are located
setwd(dir)

library(tidyr)
library(dplyr)
library(cowplot)
library(glmmTMB)
library(MuMIn)
library(broom)
library(gtsummary)
library(broom.mixed)
library(performance)

df_both_beh = readRDS("Rdata/pairwise_social_behavior_and_microbiome_df_lim.rds")

length(union(unique(df_both_beh$Name1),unique(df_both_beh$Name2)))#47
length(union(unique(df_both_beh$iso1),unique(df_both_beh$iso2)))#253
head(df_both_beh)
names(df_both_beh)
####################################################################
## Beta GLMMs: predictors of between individual pairwise microbial dissimilarity
####################################################################

####################################################################
## Model Set 1: Social group, social network distance, relatedness
####################################################################
social_related_glmm = glmmTMB(bray_curtis ~ Group + scale(Groom_PL_12mo) + scale(Prox_PL_12mo)+ Related +  (1|Year1) +(1|Name1)+ (1|Name2)+
                                (1|Extraction),data=df_both_beh,family=beta_family(link="logit"))
summary(social_related_glmm)
check_collinearity(social_related_glmm,component = "all")
# Low Correlation
#                     Term  VIF Increased SE Tolerance
# scale(Prox_PL_12mo)  3.63         1.91      0.28
# Related              1.36         1.17      0.73
# 
# Moderate Correlation
# 
#                     Term  VIF Increased SE Tolerance
# Group                7.14         2.67      0.14
# scale(Groom_PL_12mo) 8.02         2.83      0.12

## VIFs are >3, so test social network covariates in separate models

groom_only_glmm = glmmTMB(bray_curtis ~ scale(Groom_PL_12mo) + Related + (1|Year1)+(1|Name1)+ (1|Name2) + (1|Extraction),
                          data=df_both_beh,family=beta_family(link="logit")) 
summary(groom_only_glmm)
prox_only_glmm = glmmTMB(bray_curtis ~ scale(Prox_PL_12mo) + Related + (1|Year1)+(1|Name1)+ (1|Name2) + (1|Extraction),
                         data=df_both_beh,family=beta_family(link="logit")) 
summary(prox_only_glmm)
group_only_glmm = glmmTMB(bray_curtis ~ Group + Related + (1|Year1)+(1|Name1)+ (1|Name2) + (1|Extraction),
                          data=df_both_beh,family=beta_family(link="logit")) 
summary(group_only_glmm)

check_collinearity(group_only_glmm,component = "all")
check_collinearity(prox_only_glmm,component = "all")
check_collinearity(groom_only_glmm,component = "all")
# all VIFs are now <2


model.sel(groom_only_glmm,prox_only_glmm,group_only_glmm,rank = AICc)
# Model selection table 
#                 cnd((Int)) dsp((Int)) cnd(Rlt) cnd(scl(Grm_PL_12m)) cnd(scl(Prx_PL_12m)) cnd(Grp)      family df   logLik    AICc  delta weight
# group_only_glmm    0.30290          +        +                                                  + beta(logit)  8 3976.913 -7937.8   0.00      1
# groom_only_glmm    0.03657          +        +                0.252                               beta(logit)  8 3898.570 -7781.1 156.69      0
# prox_only_glmm    -0.09530          +        +                                    0.3338          beta(logit)  8 3851.789 -7687.5 250.25      0
performance::compare_performance(groom_only_glmm,prox_only_glmm,group_only_glmm,metrics=c("AIC","BIC","R2"),rank=T)
# Name            |   Model |       AIC |       BIC | R2 (cond.) | R2 (marg.) | Performance-Score
# -----------------------------------------------------------------------------------------------
# group_only_glmm | glmmTMB | -7937.827 | -7887.165 |      0.978 |      0.293 |            74.31%
# groom_only_glmm | glmmTMB | -7781.141 | -7730.478 |      0.979 |      0.295 |            45.17%
# prox_only_glmm  | glmmTMB | -7687.578 | -7636.916 |      0.991 |      0.203 |            25.00%

glance(group_only_glmm)
tbl1 = tbl_regression(group_only_glmm,
                      estimate_fun = function(x) style_ratio(x, digits = 2),
                      pvalue_fun = function(x) style_pvalue(x, digits = 3),
                      label = list(Group ~ "Group Membership", Related ~ "Relatedness"))
tbl1

broom::glance(group_only_glmm)

ft = tbl1 %>%  
  add_glance_source_note(
    label = list(df.residual  ~ "Residual df", sigma ~ "\\U03C3"),
    fmt_fun = df.residual ~ style_number,
    include = c(AIC, BIC,sigma, df.residual)
  )%>%
  bold_labels() %>%
  bold_p() %>%
  as_flex_table()
ft

tbl1 = tbl_regression(groom_only_glmm,include=c("scale(Groom_PL_12mo)","Related"),
                      estimate_fun = function(x) style_ratio(x, digits = 2),
                      pvalue_fun = function(x) style_pvalue(x, digits = 3),
                      label = list("scale(Groom_PL_12mo)" ~ "Grooming Path Length", Related ~ "Relatedness"))

ft = tbl1 %>%  
  add_glance_source_note(
    label = list(df.residual  ~ "Residual df", sigma ~ "\\U03C3"),
    fmt_fun = df.residual ~ style_number,
    include = c(AIC, BIC,sigma, df.residual)
  )%>%
  bold_labels() %>%
  bold_p() %>%
  as_flex_table()
ft

glance(prox_only_glmm)
tbl1 = tbl_regression(prox_only_glmm,include=c("scale(Prox_PL_12mo)","Related"),
                      estimate_fun = function(x) style_ratio(x, digits = 2),
                      pvalue_fun = function(x) style_pvalue(x, digits = 3),
                      label = list("scale(Prox_PL_12mo)" ~ "Proximity Path Length", Related ~ "Relatedness"))

tbl1

broom::glance(prox_only_glmm)

ft = tbl1 %>%  
  add_glance_source_note(
    label = list(df.residual  ~ "Residual df", sigma ~ "\\U03C3"),
    fmt_fun = df.residual ~ style_number,
    include = c(AIC, BIC,sigma, df.residual)
  )%>%
  bold_labels() %>%
  bold_p() %>%
  as_flex_table()
ft

####################################################################
## Model Set 2: Social group, social network distance, relatedness, and dietary distance
####################################################################
diet_df = df_both_beh %>% filter(!is.na(food_distance)&Groom_PL_12mo!=Inf&Prox_PL_12mo!=Inf)%>%droplevels()%>%
  mutate(Age_binary = ifelse(Age1==Age2,"Same Age","Different Age"))
nrow(diet_df)#1056
names(diet_df)
length(union(unique(diet_df$Name1),unique(diet_df$Name2)))#25
length(union(unique(diet_df$iso1),unique(diet_df$iso2)))#141

unique(diet_df$food_distance)
diet_glmm = glmmTMB(bray_curtis ~ food_distance +scale(Groom_PL_12mo) +scale(Prox_PL_12mo) + Group +
                      Related +
                      (1|Year1) +(1|Name1)+ (1|Name2)+ (1|Extraction),
                    data=diet_df,family=beta_family(link="logit"))
check_collinearity(diet_glmm,component = "all")
# Low Correlation
# 
#                     Term  VIF Increased SE Tolerance
# food_distance       1.11         1.06      0.90
# Related             1.40         1.18      0.71
# 
# Moderate Correlation
# 
#                     Term  VIF Increased SE Tolerance
# scale(Groom_PL_12mo) 7.96         2.82      0.13
# 
# High Correlation
# 
#                     Term   VIF Increased SE Tolerance
# scale(Prox_PL_12mo) 11.18         3.34      0.09
# Group               16.67         4.08      0.06

diet_glmm = glmmTMB(bray_curtis ~ food_distance + Group + 
                       (1|Year1) +(1|Name1)+ (1|Name2)+ (1|Extraction),
                     data=diet_df,family=beta_family(link="logit"))

diet_glmm1 = glmmTMB(bray_curtis ~ food_distance + Group + Related + 
                       (1|Year1) +(1|Name1)+ (1|Name2)+ (1|Extraction),
                     data=diet_df,family=beta_family(link="logit"))
diet_glmm2 = glmmTMB(bray_curtis ~ food_distance +
                       scale(Groom_PL_12mo) +
                       Related + 
                       (1|Year1) +(1|Name1)+ (1|Name2)+ (1|Extraction),
                     data=diet_df,family=beta_family(link="logit"))
diet_glmm3 = glmmTMB(bray_curtis ~ food_distance + 
                       scale(Prox_PL_12mo)+ 
                       Related + 
                       (1|Year1) +(1|Name1)+ (1|Name2)+ (1|Extraction),
                     data=diet_df,family=beta_family(link="logit"))
diet_glmm4 = glmmTMB(bray_curtis ~ food_distance +
                       scale(Groom_PL_12mo) +
                       (1|Year1) +(1|Name1)+ (1|Name2)+ (1|Extraction),
                     data=diet_df,family=beta_family(link="logit"))
diet_glmm5 = glmmTMB(bray_curtis ~ food_distance + 
                       scale(Prox_PL_12mo)+ 
                       (1|Year1) +(1|Name1)+ (1|Name2)+ (1|Extraction),
                     data=diet_df,family=beta_family(link="logit"))

check_collinearity(diet_glmm)
check_collinearity(diet_glmm1)
check_collinearity(diet_glmm2)
check_collinearity(diet_glmm3)
model.sel(diet_glmm, diet_glmm1,diet_glmm2,diet_glmm3,diet_glmm4,diet_glmm5,rank=AICc)
#             cnd((Int)) dsp((Int)) cnd(fod_dst) cnd(Grp) cnd(Rlt) cnd(scl(Grm_PL_12m)) cnd(scl(Prx_PL_12m))      family df  logLik    AICc delta weight
# diet_glmm      0.1403          +       0.4047        +                                                    beta(logit)  8 912.722 -1809.3  0.00  0.692
# diet_glmm1     0.1392          +       0.4079        +        +                                           beta(logit)  9 912.723 -1807.3  2.03  0.251
# diet_glmm5    -0.2185          +       0.4850                                                      0.1640 beta(logit)  8 909.771 -1803.4  5.90  0.036
# diet_glmm3    -0.1917          +       0.4186                 +                                    0.1583 beta(logit)  9 910.265 -1802.4  6.95  0.021
# diet_glmm4    -0.1870          +       0.4602                                 0.1992                      beta(logit)  8 901.592 -1787.0 22.26  0.000
# diet_glmm2    -0.1681          +       0.4221                 +               0.1938                      beta(logit)  9 901.769 -1785.4 23.94  0.000

performance::compare_performance(diet_glmm, diet_glmm1,diet_glmm2,diet_glmm3,diet_glmm4,diet_glmm5, rank=T)

summary(diet_glmm)
tbl1 = tbl_regression(diet_glmm,include=c("Group","food_distance"),
                      estimate_fun = function(x) style_ratio(x, digits = 2),
                      pvalue_fun = function(x) style_pvalue(x, digits = 3))
tbl1

ft = tbl1 %>%  
  add_glance_source_note(
    label = list(df.residual  ~ "Residual df", sigma ~ "\\U03C3"),
    fmt_fun = df.residual ~ style_number,
    include = c(AIC, BIC,sigma, df.residual)
  )%>%
  bold_labels() %>%
  bold_p() %>%
  as_flex_table()
ft
