# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

####################################################################
## Predictors of within-host gut microbiome dynamics
####################################################################

rm(list=ls())
graphics.off()

## set your working directory to where files are located
## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
setwd(dir)

library(ggplot2)
library(cowplot)
library(dplyr)
library(glmmTMB)
library(MuMIn)
library(FSA)
library(broom)
library(broom.mixed)
library(ggplot2)
library(performance)
library(gtsummary)
library(ggpubr)

load("Rdata/pairwise_predictors_same_ind_lim.RData")#df_dist_same_ind 
head(df_dist_same_ind_lim)
names(df_dist_same_ind_lim)

## limit to Groups I to V, samples more than 5 days apart
df_dist_same_ind_red = df_dist_same_ind_lim %>%
  filter(!(Group1%in%c("VI","Bachelor_1","Bachelor_2")) & !(Group2%in%c("VI","Bachelor_1","Bachelor_2")) ) %>%
  filter(Date_Diff > 5) %>% #samples collected more than 5 days apart
  droplevels()

nrow(df_dist_same_ind_red)

cols <- c("Age","Sex","Name1","Group","Extraction","Recent_Dispersal","season_diff", "iso1","iso2")
df_dist_same_ind_red[cols] <- lapply(df_dist_same_ind_red[cols], factor)
df_dist_same_ind_red$Recent_Dispersal = factor(df_dist_same_ind_red$Recent_Dispersal,levels=c("Both_Resident","Res_Imm","Both_Immigrant"))
df_dist_same_ind_red$season_diff = factor(df_dist_same_ind_red$season_diff,levels=c("within_season","one_year","three_years","four_years"))
df_dist_same_ind_red$season_binary = if_else(df_dist_same_ind_red$season_diff =="within_season","same_season","diff_season")

nrow(df_dist_same_ind_lim)#879
nrow(df_dist_same_ind_red) #719
unique(df_dist_same_ind_red$season_diff)

unique(df_dist_same_ind_red$Name1)

df_dist_same_ind_red %>%
  group_by(season_diff,Name1)%>%
  tally() %>% 
  filter(n>1) %>%
  group_by(season_diff)%>%
  tally()

## animals with more than two samples
df_dist_same_ind_red %>%
  group_by(Name1)%>%
  tally() %>% 
  filter(n>2) 

more_samples = df_dist_same_ind_red %>%
  group_by(Name1)%>%
  tally() %>% 
  filter(n>2) 
more_samples

keep = as.vector(more_samples$Name1)
length(keep)#43
df_dist_same_ind_red = as.data.frame(df_dist_same_ind_red)

names_df = table(df_dist_same_ind_red$Name1)[table(df_dist_same_ind_red$Name1)>2]
vec = names(names_df)
length(vec)#43
red = df_dist_same_ind_red %>% filter(Name1 %in% vec) %>% droplevels()
nrow(red)#709
length(union(unique(red$iso1),unique(red$iso2)))#254

levels(red$Age) <- c("No change","Adult_Juvenile", "Adult_Subadult","No change","Juvenile_Subadult","No change")
red$Age = factor(red$Age,levels=c("No change","Juvenile_Subadult","Adult_Subadult","Adult_Juvenile"))
levels(red$Recent_Dispersal)<- c("Both_Resident","Recent Immigrant","Recent Immigrant")

red$Age_binary = if_else(red$Age %in% c("Adult_Juvenile","Adult_Subadult","Juvenile_Subadult"),"Change","No change")

within = red %>% filter(season_diff=="within_season") %>% distinct(Name1)
nrow(within) #43

one = red %>% filter(season_diff=="one_year") %>% distinct(Name1)
nrow(one) #19

three = red %>% filter(season_diff=="three_years") %>% distinct(Name1)
nrow(three) #12

four = red %>% filter(season_diff=="four_years") %>% distinct(Name1)
nrow(four) #12

#all years
length(Reduce(intersect, list(one$Name1,three$Name1,four$Name1)))#10

#at least two years
length(unique(c(as.character(one$Name1),as.character(three$Name1),as.character(four$Name1))))#23

red %>%
  group_by(season_diff,Name1)%>%
  tally() %>% 
  group_by(season_diff)%>%
  tally()
# season_diff       n
# <fct>         <int>
# 1 within_season    43
# 2 one_year         19
# 3 three_years      12
# 4 four_years       12

red2 =red %>% filter(Name1 %in% keep) %>% droplevels()
nrow(red2)
length(union(unique(red2$iso1),unique(red2$iso2)))#252
length(unique(red2$Name1))

####################################################################
## Beta GLMMs: predictors of within-host gut microbiome dynamics
####################################################################

fit_beta_ind_turnover1 <- glmmTMB(bray_curtis ~ Date_Diff + Age + Sex + Recent_Dispersal + (1|Extraction) + (1|Name1), data=red2,family=beta_family(link = "logit"))
summary(fit_beta_ind_turnover1)## warning message
length(unique(intersect(red$iso1,red$iso2)))
length(unique(intersect(red$Name1,red$Name2)))

fit_beta_ind_turnover2 <- glmmTMB(bray_curtis ~ season_diff + Age + Sex + Recent_Dispersal +  (1|Extraction) + (1|Name1), data=red2, family=beta_family(link = "logit"))
summary(fit_beta_ind_turnover2)

## seasonal difference performs better than number of days
model.sel(fit_beta_ind_turnover1,fit_beta_ind_turnover2)
# Model selection table 
#                         cnd((Int)) dsp((Int)) cnd(Age) cnd(Dat_Dff) cnd(Rcn_Dsp) cnd(Sex) cnd(ssn_dff)      family df  logLik    AICc delta weight
# fit_beta_ind_turnover2    -0.6334          +        +                         +        +            + beta(logit) 12 637.682 -1250.9  0.00      1
# fit_beta_ind_turnover1    -0.6249          +        +    0.0002392            +        +              beta(logit) 10 627.701 -1235.1 15.83      0
# Models ranked by AICc(x) 
# Random terms (all models): 
#   ‘cond(1 | Extraction)’, ‘cond(1 | Name1)’

check_collinearity(fit_beta_ind_turnover2) ## all VIFs <2

check_singularity(fit_beta_ind_turnover2)
check_model(fit_beta_ind_turnover2)
summary(fit_beta_ind_turnover2)

tbl1 = tbl_regression(fit_beta_ind_turnover2,include=c("season_diff","Age" ,"Sex","Recent_Dispersal"),estimate_fun = function(x) style_ratio(x, digits = 2))
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
# print(ft, preview = "docx")

####################################################################
## Figure S20: individual predictors of within host gut microbiome turnover
####################################################################

diff_yr_df = red %>% filter(Year=="Different") 
table(diff_yr_df$Age)

med = diff_yr_df %>% filter(Age=="No change") %>% summarise(med=median(bray_curtis))
age_plot = ggplot(diff_yr_df,aes(x=Age,y=bray_curtis))+
  geom_jitter(alpha=0.5,fill="grey",size=2,pch=21,width = 0.25)+
  geom_boxplot(outlier.shape = NA,alpha=0)+
  stat_compare_means(label = "..p.signif..", method = "wilcox.test",
                     ref.group = "No change",size=7,method.args = list(p.adjust.method="BH"),label.y.npc = 0.8,vjust=-1)  + # Add pairwise comparisons p-value
  geom_hline(yintercept = med$med,lty="dashed",alpha=0.7)+
  ylab("Microbiome Dissimilarity\\n(Bray-Curtis)")+
  xlab("Shift in Age Class between Microbiome Samples\\n(Between-Year Comparisons)")+
  theme_bw(base_size = 18)+
  background_grid("xy")+
  scale_x_discrete(breaks=c("No change","Juvenile_Subadult", "Adult_Subadult","Adult_Juvenile"),
                   labels=c("No change", "Juvenile to\\n Subadult", 
                            "Subadult to\\n Adult","Juvenile to\\n Adult")) + 
  theme(legend.position = "none")
age_plot


unique(red$season_diff)
table(red$season_diff)
red$season_diff = factor(red$season_diff,levels=c("within_season","one_year","three_years","four_years"))

table(red$season_diff)

med = red %>% filter(season_diff=="within_season") %>% summarise(med=median(bray_curtis))
time_plot = ggplot(red,aes(x=season_diff,y=bray_curtis))+
  geom_jitter(alpha=0.5,fill="grey",size=3,pch=21,width = 0.25)+
  geom_boxplot(outlier.shape = NA,alpha=0)+
  stat_compare_means(label = "..p.signif..", method = "wilcox.test",
                     ref.group = "within_season",size=7,method.args = list(p.adjust.method="BH"),label.y.npc = 0.8,vjust=-1)  + # Add pairwise comparisons p-value
  geom_hline(yintercept = med$med,lty="dashed",alpha=0.7)+
  ylab("Microbiome Dissimilarity\\n(Bray-Curtis)")+  
  xlab("Time between Microbiome Samples")+
  theme_bw(base_size = 18)+
  background_grid("xy")+
  ylim(0.1,0.8)+
  scale_x_discrete(breaks=c("within_season","one_year", "three_years", "four_years"),
                   labels=c("same\\nseason\\n", "one\\n year\\n", 
                            "three\\n years\\n","four\\n years\\n")) + 
  theme(legend.position = "none")
time_plot

names(diff_yr_df)
table(red$Recent_Dispersal)

dispersal_plot = ggplot(red,aes(x=Recent_Dispersal,y=bray_curtis))+
  geom_jitter(alpha=0.5,fill="grey",size=3,pch=21,width = 0.25)+
  geom_boxplot(outlier.shape = NA,alpha=0)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     size=7,label.x = 1.5, 
                     label.y.npc = 0.8,vjust=-1)  + # Add pairwise comparisons p-value
  ylab("Microbiome Dissimilarity\\n(Bray-Curtis)")+
  xlab("Dispersal Status between Microbiome Samples")+
  theme_bw(base_size = 18)+
  background_grid("xy")+
  scale_x_discrete(breaks=c("Both_Resident","Recent Immigrant"),
                   labels=c("Group Resident", "Recent Immigrant")) + 
  theme(legend.position = "none")
dispersal_plot

table(red$Sex)

sex_plot = ggplot(red,aes(x=Sex,y=bray_curtis))+
  geom_jitter(alpha=0.5,fill="grey",size=3,pch=21,width = 0.25)+
  geom_boxplot(outlier.shape = NA,alpha=0)+
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     size=7,label.x = 1.5, label.y.npc = 0.8,vjust=-1)  + # Add pairwise comparisons p-value
  ylab("Microbiome Dissimilarity\\n(Bray-Curtis)")+
  xlab("Sex")+
  theme_bw(base_size = 18)+
  background_grid("xy")+
  scale_x_discrete(breaks=c("Female_Female","Male_Male"),
                   labels=c("Female", "Male")) + 
  theme(legend.position = "none")
sex_plot

all_plots = plot_grid(time_plot,age_plot,dispersal_plot,sex_plot,align = "hv",labels=c("a","b","c","d"))
all_plots
save_plot(plot=all_plots,filename = "figures/SFig20_bray_curtis_within_host_plots.pdf", base_width = 16,base_height = 12)
