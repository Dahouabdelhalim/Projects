# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

####################################################################
## microbial similarity between long-term residents and immigrants
####################################################################

rm(list=ls())
graphics.off()

## set your working directory to where files are located
## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
setwd(dir)

library(phyloseq)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggpubr)
library(ggformula)
library(glmmTMB)
library(performance)
library(gtsummary)

## limit analysis to individuals within the same group
df_dist_bw = readRDS("Rdata/imm_res_diss_df_lim.rds")
head(df_dist_bw)
names(df_dist_bw)
############################################################################################
## Figures 5 and S21: Resident-Recent Immigrant Microbial Dissimilarity
############################################################################################

p0 <- ggplot(df_dist_bw %>% filter(Recent_Dispersal!="Both_Immigrant") %>% droplevels(),aes(factor(Recent_Dispersal),bray_curtis))
p2 <- p0 + geom_violin(fill="grey",alpha=0.3, trim=F) +
  geom_boxplot(width=0.2, color="black")+
  ylab("Microbiome Dissimilarity\\n(Bray-Curtis)")+
  xlab("") +
  scale_x_discrete(labels=c("Both Residents", "Immigrant-\\nResident"))
imm_violin_bc=p2 +stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test", label.x.npc = 0.4,label.y.npc = 1,
                          size=10)+
              theme_cowplot()+ background_grid(major="xy")
imm_violin_bc

p0 <- ggplot(df_dist_bw %>% filter(Recent_Dispersal!="Both_Immigrant") %>% droplevels(),aes(factor(Recent_Dispersal),w_unifrac))
p2 <- p0 + geom_violin(fill="grey",alpha=0.3, trim=F) +
  geom_boxplot(width=0.2, color="black")+
  ylab("Microbiome Dissimilarity\\n(Weighted Unifrac)")+
  xlab("") +
  scale_x_discrete(labels=c("Both Residents", "Immigrant-\\nResident", "Both Immigrants"))
imm_violin_wu=p2 +stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test", 
                                     label.x.npc = 0.4,label.y.npc = 1,size=6)+
              theme_cowplot()+ 
              background_grid(major="xy")

group_tenure = readRDS("Rdata/imm_group_tenure_df_lim.rds")
head(group_tenure)
names(group_tenure)

## limit to samples collected within two weeks of each other
ind_df = group_tenure %>%
  filter(Date_Diff<=14)

rec_imm_bc <- ggplot(ind_df%>% filter(tenure_mo<=12), aes(x=tenure_mo, y=bray_curtis)) +
  geom_jitter(pch=21,fill="grey",alpha=0.8,size=4) +
  # stat_cor(method="spearman",aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),label.x.npc ="right",hjust=1)+
  xlab("Length of Group Tenure (Months)")+
  ylab("Microbiome Dissimilarity to Residents\\n(Bray-Curtis)") +
  geom_smooth(method="lm",se=F,color="black",lwd=0.75,lty="dashed")+
  theme(legend.position = "none")+
  theme_bw(base_size = 12)

rec_imm_wu <- ggplot(ind_df%>% filter(tenure_mo<=12), aes(x=tenure_mo, y=w_unifrac)) +
  geom_jitter(pch=21,fill="grey",alpha=0.8,size=4) +
  # stat_cor(method="spearman",aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),label.x.npc ="right",hjust=1)+
  xlab("Length of Group Tenure (Months)")+
  ylab("Microbiome Dissimilarity to Residents\\n(Weighted Unifrac)") +
  geom_smooth(method="lm",se=F,color="black",lwd=0.75,lty="dashed")+
  theme(legend.position = "none")+
  theme_bw(base_size = 12)
rec_imm_wu

combined_group_tenure_bc = plot_grid(imm_violin_bc,rec_imm_bc,labels = c("a","b"))
combined_group_tenure_bc
save_plot(combined_group_tenure_bc,filename = "figures/Fig5_bray_curtis_group_tenure_vs_microbiome_diss.pdf",base_width = 10,base_height = 5)
save_plot(combined_group_tenure_bc,filename = "figures/Fig5_bray_curtis_group_tenure_vs_microbiome_diss.png",base_width = 10,base_height = 5,dpi=600)

combined_group_tenure_wu = plot_grid(imm_violin_wu,rec_imm_wu,labels = c("a","b"))
save_plot(combined_group_tenure_wu,filename = "figures/SFig21_weighted_unifrac_group_tenure_vs_microbiome_diss.pdf",base_width = 10,base_height = 5)

############################################################################################
## Figure S22: Resident - Past Immigrant Microbial Dissimilarity
############################################################################################
past_imm_wu <- ggplot(group_tenure%>% filter(tenure_year>=1 & tenure_year<=6), aes(x=tenure_year, y=w_unifrac)) + 
  geom_jitter(pch=21,fill="grey",alpha=0.8,size=4) +
  # geom_jitter(pch=21,fill="black",alpha=0.5,size=4) +
  xlab("Length of Group Tenure (Years)")+
  ylab("Microbiome Dissimilarity to Residents\\n(Weighted Unifrac)") +
  geom_smooth(method="lm",se=F,color="black",lwd=0.75,lty="dashed")+
  theme_bw(base_size = 14)

past_imm_bc <- ggplot(group_tenure %>% filter(tenure_year>=1 & tenure_year<6), aes(x=tenure_year, y=bray_curtis)) + 
  geom_jitter(pch=21,fill="grey",alpha=0.8,size=4) +
  xlab("Length of Group Tenure (Years)")+
  geom_smooth(method="lm",se=F,color="black",lwd=0.75,lty="dashed")+
  ylab("Microbiome Dissimilarity to Residents\\n(Bray-Curtis)") +
  theme_bw(base_size = 14)
past_imm_bc

combined_immigrant = plot_grid(past_imm_bc,past_imm_wu,labels = c("a","b"))
combined_immigrant
save_plot(combined_immigrant,filename = "figures/SFig22_immigrant_one_year_plus.pdf",base_width = 12,base_height = 5)

ind_df = ind_df%>% 
  rowwise() %>%      # for each row
  mutate(name_pairs = paste(sort(c(as.character(Name1), as.character(Name2))), collapse = "_"))
  
recent = ind_df%>% filter(tenure_mo<=12)
recent %>%distinct(name_pairs)%>%nrow()#23
length(union(unique(recent$iso1),unique(recent$iso2)))#65
past = ind_df %>% filter(tenure_year>=1 & tenure_year<6) 
past %>%distinct(name_pairs)%>%nrow()#45
length(union(unique(past$iso1),unique(past$iso2)))#172
############################################################################################
## Beta GLMMs: Resident - Recent Immigrant Microbial Dissimilarity vs Group Tenure
############################################################################################
### recent immigrants
fit_beta_ind_turnover1 <- glmmTMB(bray_curtis ~ scale(tenure_mo) + (1|Name1) + (1|Name2) + (1|Extraction), data=ind_df%>% filter(tenure_mo<=12), family= beta_family(link="logit"))
summary(fit_beta_ind_turnover1)
# Dispersion parameter for beta family (): 17.2 
# 
# Conditional model:
#                   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)      -0.23499    0.07370  -3.189  0.00143 **
# scale(tenure_mo) -0.11700    0.04081  -2.867  0.00414 **

fit_beta_ind_turnover2 <- glmmTMB(w_unifrac ~ scale(tenure_mo) + (1|Name1) + (1|Name2) + (1|Extraction), data=ind_df%>% filter(tenure_mo<=12), family= beta_family(link="logit"))
summary(fit_beta_ind_turnover2)
# Dispersion parameter for beta family (): 53.2 
# 
# Conditional model:
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -1.66458    0.05902 -28.203   <2e-16 ***
# scale(tenure_mo) -0.08905    0.03836  -2.322   0.0203 *  

############################################################################################
## Beta GLMMs: Resident - Past Immigrant Microbial Dissimilarity vs Group Tenure
############################################################################################
## past immigrants
fit_beta_ind_turnover1 <- glmmTMB(bray_curtis ~ scale(tenure_mo) + (1|Name1) + (1|Name2) + (1|Extraction), data=ind_df%>% filter(tenure_year>=1 & tenure_year<=6), family= beta_family(link="logit"))
summary(fit_beta_ind_turnover1)
fit_beta_ind_turnover2 <- glmmTMB(w_unifrac ~ scale(tenure_mo) + (1|Name1) + (1|Name2) + (1|Extraction), data=ind_df%>% filter(tenure_year>=1 & tenure_year<=6), family= beta_family(link="logit"))
summary(fit_beta_ind_turnover2)
##########################################
### social network distance between male recent immigrants and female residents
rec_imm_df_MF = group_tenure  %>% filter(tenure_mo<=12 & Sex=="Female_Male") %>% droplevels()
prox_MF <- rec_imm_df_MF %>% filter(Prox_PL_12mo!=Inf)
groom_MF <- rec_imm_df_MF %>% filter(Groom_PL_12mo!=Inf)

## male/female proximity
bc_MF = glmmTMB(bray_curtis~ scale(Prox_PL_12mo)+ (1|Name1) + (1|Name2) + (1|Extraction),data=prox_MF,family= beta_family(link="logit"))
summary(bc_MF)
wu_MF = glmmTMB(w_unifrac~ scale(Prox_PL_12mo) + (1|Name1) + (1|Name2) + (1|Extraction),data=prox_MF,family= beta_family(link="logit"))
summary(wu_MF)

## male/female grooming
bc_MF = glmmTMB(bray_curtis~ scale(Groom_PL_12mo)+(1|Name1) + (1|Name2) + (1|Extraction),data=groom_MF,family= beta_family(link="logit"))
summary(bc_MF)
wu_MF = glmmTMB(w_unifrac~ scale(Groom_PL_12mo)+(1|Name1) + (1|Name2) + (1|Extraction),data=groom_MF,family= beta_family(link="logit"))
summary(wu_MF)