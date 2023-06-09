# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

####################################################################
## Social Partner Stability and Within Host Gut Microbiome Dynamics
####################################################################

rm(list=ls())
graphics.off()

## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
setwd(dir)

library(phyloseq)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(performance)
library(MuMIn)
library(jcolors)
library(gtsummary)
library(ggplot2)
library(cowplot)

both_beh = readRDS("Rdata/PSI_vs_microbiome_df_lim.rds")
nrow(both_beh)#199
length(unique(both_beh$Name))#12
head(both_beh)
####################################################################
## Beta GLMMs: grooming PSI and proximity PSI vs within host turnover
####################################################################
fit_beta_ind_turnover1 <- glmmTMB(bray_curtis ~  groom_PSI + prox_PSI + Sex + (1|season_diff)+ (1|Extraction), data=both_beh, family= beta_family(link="logit"))
summary(fit_beta_ind_turnover1)
check_collinearity(fit_beta_ind_turnover1)
psi_dredge = dredge(fit_beta_ind_turnover1)
psi_dredge # sex does not improve model fit

fit_beta_ind_turnover1 <- glmmTMB(bray_curtis ~  groom_PSI + prox_PSI + (1|season_diff)+ (1|Extraction), data=both_beh, family= beta_family(link="logit"))
summary(fit_beta_ind_turnover1)
check_collinearity(fit_beta_ind_turnover1)
psi_dredge = dredge(fit_beta_ind_turnover1)
psi_dredge # only proximity retained in best fitting model

tbl1 = tbl_regression(fit_beta_ind_turnover1,include=c("groom_PSI","prox_PSI"),
                      estimate_fun = function(x) style_ratio(x, digits = 2),
                      pvalue_fun = function(x) style_pvalue(x, digits = 3))
tbl1
ft = tbl1 %>%  
  add_glance_source_note(
    label = list(df.residual  ~ "Residual df", sigma ~ "\\U03C3"),
    fmt_fun = df.residual ~ style_number,
    include = c(AIC, BIC,sigma, df.residual)
  )%>%
  # bold_labels() %>%
  bold_p() %>%
  as_flex_table()
ft

both_beh$season_diff = factor(both_beh$season_diff,levels=c("one year","three years","four years"))

####################################################################
## Figure 6
####################################################################

full_prox = ggplot(both_beh,aes(x=prox_PSI,y=bray_curtis))+
  geom_jitter(size=4,pch=21,fill="grey",alpha=0.8)+
  geom_smooth(method="glm",se=F,lwd=0.7,lty="dashed",color="black")+
  ylab("Microbiome Dissimilarity to Self\\n(Bray-Curtis)")+
  xlab("Proximity Partner Stability Index (PSI)")+
  scale_color_manual(values="grey")+
  theme_bw(base_size=14)+
  scale_fill_jcolors(palette = "pal8")+
  theme(legend.position = "none",legend.title.align = 0.5,plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=10),legend.title = element_text(size=12))
full_prox

full_groom = ggplot(both_beh,aes(x=groom_PSI,y=bray_curtis))+
  geom_jitter(size=4,pch=21,fill="grey",alpha=0.8)+
  geom_smooth(method="glm",se=F,lwd=0.7,lty="dashed",color="black")+
  ylab("Microbiome Dissimilarity to Self\\n(Bray-Curtis)")+
  xlab("Grooming Partner Stability Index (PSI)")+
  scale_color_manual(values="grey")+
  theme_bw(base_size=14)+
  theme(legend.position = "none",legend.title.align = 0.5,plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=10),legend.title = element_text(size=12))
full_groom

combined_psi_plots = plot_grid(full_prox,full_groom,nrow=1,labels=c("a","b"))
save_plot(combined_psi_plots,filename = "figures/Fig6_combined_years_psi_plots.pdf",base_width = 12,base_height = 6)
save_plot(combined_psi_plots,filename = "figures/Fig6_combined_years_psi_plots.png",dpi=600,base_width = 12,base_height = 6)

# ## individual seasonal comparisons
# q = ggplot(both_beh,aes(x=prox_PSI,y=bray_curtis))+
#   geom_point(size=4,pch=21,aes(fill=season_diff),alpha=0.8)+
#   geom_smooth(method="glm",se=F,lwd=0.7,lty="dashed",color="black")+
#   ylab("Microbiome Dissimilarity to Self\\n(Bray-Curtis)")+
#   xlab("Proximity Partner Stability Index (PSI)")+
#   scale_color_manual(values="grey")+
#   theme_bw(base_size=14)+
#   scale_fill_jcolors(palette = "pal8")+
#   theme(legend.position = "none",legend.title.align = 0.5,plot.title = element_text(hjust = 0.5),
#         legend.text=element_text(size=10),legend.title = element_text(size=12))+
#   facet_wrap(~season_diff,scales="free_x")
# q
# r = ggplot(both_beh,aes(x=groom_PSI,y=bray_curtis))+
#   geom_point(size=4,pch=21,aes(fill=season_diff),alpha=0.8)+
#   geom_smooth(method="glm",se=F,lwd=0.75,lty="dashed",color="black")+
#   ylab("Microbiome Dissimilarity to Self\\n(Bray-Curtis)")+
#   xlab("Grooming Partner Stability Index (PSI)")+
#   scale_color_manual(values="grey")+
#   theme_bw(base_size=14)+
#   scale_fill_jcolors("pal8")+
#   theme(legend.position = "none",legend.title.align = 0.5,plot.title = element_text(hjust = 0.5),
#         legend.text=element_text(size=10),legend.title = element_text(size=12))+
#   facet_wrap(~season_diff,scales="free_x")
# 
# r