### plot correlation between absolute values and soil parameters

# This script creates Figure 3. 
# Relationships were tested based on indicated Spearman's coefficients
# of correlation according to the correlogram of Figure S1. 

library(tidyverse)
library(RColorBrewer)
library(ggpmisc)
library(ggpubr)
library(cowplot)
#--------------------------

# load script, which loads needed data
source("correlograms_absolute_soil_params_script.R")  


my.formula <- y~x
my.formula.log <- y ~ log(x)
my.formula.exp <- y ~ exp(x)
my.formula.pol <- y ~ poly(x, 2, raw = TRUE)

# create ordered data.frame of mean +- sd values per plot for plotting
meansd.slu <- sample_data %>%
  group_by(site, landuse) %>%
  summarise_if(is.numeric, list(mean=mean, sd=sd)) %>%
  ungroup()
meansd_compsoil <- merge(meansd.slu, compsoil) 
meansd_compsoil_ordered <- meansd_compsoil
meansd_compsoil_ordered$landuse <- ordered(meansd_compsoil_ordered$landuse, 
                                                   levels = c("Forest", "Grassland", "Cropland"))

# create  data.frame of mean +- sd values per sample for plotting
meansd <- sample_data %>%
  group_by(site, landuse, replicate) %>%
  summarise_if(is.numeric, list(mean=mean,sd=sd)) %>%
  ungroup()
meansd_gensoil <- merge(meansd, gensoil)



#### Figure 3 - CUE ~ pH ####
ggplot(data=meansd_compsoil_ordered, aes(x=pH_1.5wvH2O, y=CUE_mean))+ 
  geom_errorbar(aes(ymin=CUE_mean-CUE_sd, 
                    ymax=CUE_mean+CUE_sd), color = "grey35", width = 0) +
  stat_smooth(method="lm", formula = my.formula, se=TRUE, color = "grey50", linetype = "dotted") +
  geom_point(aes(fill=site, shape = landuse), size = 3, color = "grey35") +
  stat_poly_eq(formula = my.formula, 
               aes(label =  paste(stat(adj.rr.label), 
                                  stat(p.value.label), 
                                  sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(pH[H2O])) +
  ylab(expression(CUE)) +
  scale_shape_manual(values = c(24,21,22)) +
  scale_fill_manual(values = c("#8F7700FF", "#999999", "#EFC000FF"))
# export svg at 600x460 and adjust legend in inkscape


#### CUE ~ F:B ####
ggplot(data=meansd_gensoil, aes(x=F.B, y=CUE_mean)) +
  geom_point(aes(colour=site, shape = landuse), size = 2) +
  stat_smooth(method = "lm", formula = my.formula.pol, se=T, fill=NA, color = "grey50", linetype = "dotted") +
  stat_poly_eq(formula = my.formula.pol, 
               aes(label =  paste(stat(eq.label),
                                  stat(rr.label), 
                                  stat(p.value.label), 
                                  sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  geom_errorbar(aes(ymin=CUE_mean-CUE_sd, 
                    ymax=CUE_mean+CUE_sd)) +
  theme_bw(base_size = 13) + 
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression("F:B"~ratio)) +
  ylab(expression(CUE)) +
  scale_shape_manual(values = c(17,16,15)) +
  scale_colour_manual(values = c("#8F7700FF", "#999999", "#EFC000FF"))
# 730x530


####  CUE ~ gene copy numbers of Fungi ####
ggplot(data=meansd_gensoil, aes(x=F_cn_gsoil, y=CUE_mean))+ 
  geom_point(aes(colour=site, shape = landuse)) +
  stat_smooth(method = "lm", formula = my.formula.pol, se=T, fill=NA, color = "grey50", linetype = "dotted") +
  stat_poly_eq(formula = my.formula.pol, 
               aes(label =  paste(stat(eq.label),
                                  stat(rr.label), 
                                  stat(p.value.label), 
                                  sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  geom_errorbar(aes(ymin=CUE_mean-CUE_sd, 
                    ymax=CUE_mean+CUE_sd)) +
  theme_bw(base_size = 13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab(expression(mean~CUE)) +
  xlab(expression(Fungi~"["~gene~copies~g^-1~soil~DW~"]")) +
  scale_shape_manual(values = c(17,16,15)) +
  scale_colour_manual(values = c("#8F7700FF", "#999999", "#EFC000FF"))


#### CUE ~ clay content ####
ggplot(data=meansd_compsoil, aes(x=clay_perc, y=CUE_mean)) + 
  geom_point(aes(colour=site, shape = landuse)) +
  stat_smooth(method = "lm", formula = my.formula, se=T, fill=NA, color = "grey50", linetype = "dotted") +
  stat_poly_eq(formula = my.formula, 
               aes(label =  paste(stat(eq.label),
                                  stat(rr.label), 
                                  stat(p.value.label), 
                                  sep = "*\\", \\"*")),
               rr.digits = 3, coef.digits = 4, 
               parse = TRUE, label.y = 0.98) +
  geom_errorbar(aes(ymin=CUE_mean-CUE_sd, 
                    ymax=CUE_mean+CUE_sd), width = 0.01) +
  theme_bw(base_size = 13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab(expression(clay~"%")) +
  ylab(expression(CUE)) + 
  scale_colour_manual(values = c("#8F7700FF", "#999999", "#EFC000FF"))
