### absolute values and ranges

# This script creates Figure S3 and Figure 2.
# Calculation of means +- sd and ranges as reported in the Results section.

library(tidyverse)
library(RColorBrewer)
library(ggpmisc)
library(ggpubr)
library(cowplot)
#--------------------------
# load script, which loads needed data
source("load_data_script.R")

# for cld letters load 
source("absolute_linear_mixed_effects_model_script.R")


#### incubation over time figures ####

#### Figure S3 - cumulative resp over time ####
ggplot(data=incub_resp, aes(x=d_of_incub, y=(cumul_resp_total_time/1000000/(Corg_per_mb_start/100)), group = ID, colour = treatment)) +
  geom_point(aes(shape = as.factor(replicate)), size = 2, alpha = 0.4) +
  stat_smooth(aes(group=treatment, fill=treatment),method="loess", se=F) +
  theme_bw(base_size = 13) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text()) +
  guides(shape = "none") +
  theme(legend.position = "bottom") +
  xlab(expression(days~of~incubation)) +
  ylab(expression(cumulative~respiration~"["~mg~CO[2]-C~g^-1~SOC~"]")) +
  scale_colour_manual(values = c("#999999", "#EFC000FF", "#0073C2FF"),
                      labels = c("10°C", "20°C", "20°C + N")) +
  scale_fill_manual(values = c("#999999", "#EFC000FF", "#0073C2FF"),
                    labels = c("10°C", "20°C", "20°C + N")) +
  facet_wrap(site~landuse, scales="free")


# additional plot
#----------------
# respiration rate over time
ggplot(data=na.omit(incub_resp), aes(x=d_since_start_num, y=respiration_ngC_g_h_n_at_20.C, group = ID, colour = treatment)) +
  geom_point(aes(shape = as.factor(replicate)), size = 2, alpha = 0.4) +
  stat_smooth(aes(group=treatment, fill=treatment),method="loess", se=F) +
  theme_bw(base_size = 13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(days~of~incubation)) +
  ylab(expression(respiration~rate~"["~ng~C~g^-1~h^-1~"]")) +
  scale_colour_manual(values = c("#999999", "#EFC000FF", "#0073C2FF")) +
  facet_grid(site~landuse, scales="free")


#### incubation experiment ####

# cumulative respiration per g SOC
cumulresp.Corg.absolute.fig <- ggplot(data=sample_data, aes(x=landuse, y=cumul_resp_Corg, color = site))+#, color = treatment)) +
  geom_boxplot(aes(fill=site), alpha = 0.3) +
  stat_summary(aes(group=site), geom = "point", fun = "mean", alpha=1, size = 2.5, color = "black")+  
  theme_bw(base_size = 13) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.y = element_text(angle = 45)) +
  theme(strip.background = element_blank()) +
  guides(shape="none", fill="none", color="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(NULL) +
  ylab(expression(cumulative~respiration~"["~mg~CO[2]-C~g^-1~SOC~"]")) +
  scale_color_manual(values = c("#8F7700FF", "#999999", "#EFC000FF")) +
  scale_fill_manual(values = c("#8F7700FF", "#999999", "#EFC000FF")) +
  facet_grid(.~site) +
  geom_text(data = letters.cumulrespCorg.interactive, aes(group = site, x = landuse, y = 0, label = group),
            position = position_dodge(width = 0.6), size = 4) 

cumulresp.Corg.absolute.fig




# additional plot
#----------------
# cumulative respiration per g soil per site, land use and treatment
ggplot(data=sample_data, aes(x=treatment, y=cumul_resp_incub_ngC_g_total_time/1000000, colour = treatment)) +
  geom_point() +
  theme_bw(base_size = 13) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(treatment)) +
  ylab(expression(cumulative~respiration~"["~mg~CO[2]-C~g^-1~"]")) +
  guides(colour = "none") +
  scale_colour_manual(values = c("#999999", "#EFC000FF", "#0073C2FF")) +
  facet_grid(site~landuse, scales="free")




#### CUE experiment ####

# create ggplot containing boxplots of all three treatments (n=9) and mean values per land use
# results differ largely depending whether CUE_NPOC or CUE (CD and NO as DIFF and SI as NPOC)
CUE.absolute.fig <- ggplot(data=sample_data, aes(x=landuse, y=CUE, color = site))+
  geom_boxplot(aes(fill=site), alpha = 0.3) +
  stat_summary(aes(group=site), geom = "point", fun = "mean", alpha=1, size = 2.5, color = "black")+
  theme_bw(base_size = 13) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.y = element_text(angle = 45)) +
  theme(strip.background = element_blank()) +
  guides(shape="none", fill="none", color="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(NULL) +
  ylab(expression(CUE)) +
  scale_color_manual(values = c("#8F7700FF", "#999999", "#EFC000FF")) +
  scale_fill_manual(values = c("#8F7700FF", "#999999", "#EFC000FF")) + 
  facet_grid(.~site) +
  geom_text(data = letters.CUE.interactive, aes(group = site, x = landuse, y = 0, label = group),
            position = position_dodge(width = 0.6), size = 4) 
CUE.absolute.fig


#### Figure 2 - absolute values per site and land use ####

plot_grid(cumulresp.Corg.absolute.fig +
            theme(legend.position = "none",
                  axis.title.y = element_text(size=13),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  strip.text = element_text(size=13)),
          CUE.absolute.fig +
            theme(legend.position = "none",
                  axis.title.y = element_text(size=13),
                  axis.text.x = element_text(size=13),
                  strip.text = element_blank()),    
          labels = c("a","b"),
          vjust = 1,
          ncol = 1,
          nrow = 2,
          align = "hv", 
          axis = "rl", 
          rel_heights = c(1,1))
# export as svg at 800x600 and adjust in inkscape



#### means +- sd, ranges ####

#### meansdallincub ####
# create data frame of mean and sd values for selected params of incub all data per landuse x site (including all treatments)
meansdallincub <- sample_data %>%
  group_by(site, landuse) %>%
  summarise_if(is.numeric, list(mean=mean, sd=sd)) %>%
  ungroup()

# calculate mean +- sd of turnover for each treatment (for results section)
sample_data %>%
  group_by(treatment) %>%
  summarise(mean_turnover = mean(turnover_d),
            sd_turnover = sd(turnover_d))  %>%
  ungroup()


# sample_data ranges
#----------------------
sample_data %>%
  group_by(site) %>% 
  summarise(min = min(cumul_resp_incub_ngC_g_total_time/1000000),
            max = max(cumul_resp_incub_ngC_g_total_time/1000000),
            mean = mean(cumul_resp_incub_ngC_g_total_time/1000000),
            sd = sd(cumul_resp_incub_ngC_g_total_time/1000000))

sample_data %>%
  group_by(site) %>% 
  summarise(min = min(CUE),
            max = max(CUE),
            mean = mean(CUE),
            sd = sd(CUE))

sample_data %>%
  group_by(site, landuse) %>% 
  summarise(mean = mean(Cmic_Corg_perc),
            sd = sd(Cmic_Corg_perc))

sample_data %>%
  group_by(site, landuse) %>% 
  summarise(mean = mean(Cmic_V1),
            sd = sd(Cmic_V1))

# summarise general microbial parameters for table of results per plot
#---------------------------------------------------------------------
# B_cn_gsoil
# B.A
# F.B
# Cmic_V1_ugC_gsoil
# Nmic_V1_ugN_gsoil
# Cmic_V1_ugC_gsoil/Nmic_V1_ugN_gsoil
# CUE_V1
# Cmic_Corg_perc
gensoil %>% group_by(site, landuse) %>%
  summarise(min = min(Nmic_V1_ugN_gsoil),
            max = max(Nmic_V1_ugN_gsoil),
            mean = mean(Nmic_V1_ugN_gsoil),
            sd = sd(Nmic_V1_ugN_gsoil))

gensoil %>%
  group_by(site, landuse) %>% 
  summarise(min = min(F_cn_gsoil),
            max = max(F_cn_gsoil),
            mean = mean(F_cn_gsoil),
            sd = sd(F_cn_gsoil))

gensoil %>%
  group_by(site, landuse) %>% 
  summarise(min = min(F.B),
            max = max(F.B),
            mean = mean(F.B),
            sd = sd(F.B))

compsoil %>%
  group_by(site, landuse) %>% 
  summarise(CP_ratio = CP_ratio)

# calculate median value of F:B ratio (for results section)
median(gensoil$F.B) 
