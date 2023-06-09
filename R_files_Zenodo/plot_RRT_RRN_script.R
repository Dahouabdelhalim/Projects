### plot RRT and RRN

# This script creates Figure 4.

library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(ggpmisc) # needed to add fromula to plot (stat_poly_eq_())
library(ggpubr) # can be used to combine and arrange multiple plots
library(cowplot)
#--------------------------
# Run the following script to create cld letters and save them in environment 
source("RRT_RRN_linear_mixed_effects_models_script.R")


#### cumul resp ####
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.cumulresp <- multipleRRTN %>%
  group_by(RR, landuse) %>%
  summarise(n = length(cumul_resp_incub_ngC_g_total_time),
            mean = mean(cumul_resp_incub_ngC_g_total_time),
            sd = sd(cumul_resp_incub_ngC_g_total_time),
            se = sd(cumul_resp_incub_ngC_g_total_time)/sqrt(length(cumul_resp_incub_ngC_g_total_time)),
            P = t.test(cumul_resp_incub_ngC_g_total_time, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(cumul_resp_incub_ngC_g_total_time)) %>%
  ungroup()

sum.multipleRRTN.cumulresp$landuse <- ordered(sum.multipleRRTN.cumulresp$landuse, 
                                              levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.cumulresp <- merge(sum.multipleRRTN.cumulresp, letters.cld.cumulresp, by = c("RR","landuse"))


#### Cresp ####
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.Cresp <- multipleRRTN %>%
  group_by(RR, landuse) %>%
  summarise(n = length(Cresp_CUE_ngC_gDW_h),
            mean = mean(Cresp_CUE_ngC_gDW_h),
            sd = sd(Cresp_CUE_ngC_gDW_h),
            se = sd(Cresp_CUE_ngC_gDW_h)/sqrt(length(Cresp_CUE_ngC_gDW_h)),
            P = t.test(Cresp_CUE_ngC_gDW_h, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(Cresp_CUE_ngC_gDW_h)) %>%
  ungroup()

sum.multipleRRTN.Cresp$landuse <- ordered(sum.multipleRRTN.Cresp$landuse,  
                                          levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.Cresp <- merge(sum.multipleRRTN.Cresp, letters.cld.Cresp, by = c("RR","landuse"))


#### CUE ####
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.CUE <- multipleRRTN %>%
  group_by(RR, landuse) %>%
  summarise(n = length(CUE),
            mean = mean(CUE),
            sd = sd(CUE),
            se = sd(CUE)/sqrt(length(CUE)),
            P = t.test(CUE, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(CUE)) %>%
  ungroup()

sum.multipleRRTN.CUE$landuse <- ordered(sum.multipleRRTN.CUE$landuse,  
                                        levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.CUE <- merge(sum.multipleRRTN.CUE, letters.cld.CUE, by = c("RR","landuse"))

t.test(multipleRRT$CUE, mu=1)
t.test(multipleRRN$CUE, mu=1)


#### figure for Eurosoil poster ####
# plot RRT and RRN per landuse
plot.sum.multipleRRTN.CUE <- sum.multipleRRTN.CUE
plot.sum.multipleRRTN.CUE$landuse <- ordered(plot.sum.multipleRRTN.CUE$landuse,  
                                             levels = c("Forest", "Grassland", "Cropland"))
plot.sum.multipleRRTN.CUE$RR <- as.factor(plot.sum.multipleRRTN.CUE$RR)

ggplot(plot.sum.multipleRRTN.CUE, aes(x= landuse, y=mean, fill=RR)) +
  geom_pointrange(aes(ymin = mean-sd, ymax = mean+sd, group=RR, shape = landuse), position = position_dodge(0.3)) +
  theme_bw(base_size = 13) +
  geom_hline(aes(yintercept = 1.0), colour = "#999999", linetype = "dashed")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(.~RR, scales = "free") +
  coord_flip()+
  theme(strip.background = element_rect(fill = "#999999")) +
  theme(strip.background = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(RR~CUE)) +
  scale_shape_manual(values = c(24,21,22)) +
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF")) +
  geom_text(data = letters.cld.CUE, aes(group = RR, x = landuse, y = 0.5, label = group),
            position = position_dodge(width = 0.8), size = 4) +
  geom_text(aes(label = Sig, y = mean+sd + 0.2), size = 4,
            data = plot.sum.multipleRRTN.CUE)



#### CGrowth ####
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.Cgrowth <- multipleRRTN %>%
  group_by(RR, landuse) %>%
  summarise(n = length(Cgrowth_ngC_gDW_h),
            mean = mean(Cgrowth_ngC_gDW_h),
            sd = sd(Cgrowth_ngC_gDW_h),
            se = sd(Cgrowth_ngC_gDW_h)/sqrt(length(Cgrowth_ngC_gDW_h)),
            P = t.test(Cgrowth_ngC_gDW_h, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(Cgrowth_ngC_gDW_h)) %>%
  ungroup()

sum.multipleRRTN.Cgrowth$landuse <- ordered(sum.multipleRRTN.Cgrowth$landuse,  
                                            levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.Cgrowth <- merge(sum.multipleRRTN.Cgrowth, letters.cld.Cgrowth, by = c("RR","landuse"))


#### Cmic at end of incubation ####
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.Cmicend <- na.omit(multipleRRTN) %>%  # one missing value has to be removed (CD_CM_mb_2)
  group_by(RR, landuse) %>%
  summarise(n = length(Cmic_ugC_gDW),
            mean = mean(Cmic_ugC_gDW),
            sd = sd(Cmic_ugC_gDW),
            se = sd(Cmic_ugC_gDW)/sqrt(length(Cmic_ugC_gDW)),
            P = t.test(Cmic_ugC_gDW, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(Cmic_ugC_gDW)) %>%
  ungroup()

sum.multipleRRTN.Cmicend$landuse <- ordered(sum.multipleRRTN.Cmicend$landuse,  
                                            levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.Cmicend <- merge(sum.multipleRRTN.Cmicend, letters.cld.Cmicend, by = c("RR","landuse"))



#### Nmic at end of incubation ####                   
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.Nmicend <- multipleRRTN %>%  
  group_by(RR, landuse) %>%
  summarise(n = length(Nmic_ugN_gDW),
            mean = mean(Nmic_ugN_gDW),
            sd = sd(Nmic_ugN_gDW),
            se = sd(Nmic_ugN_gDW)/sqrt(length(Nmic_ugN_gDW)),
            P = t.test(Nmic_ugN_gDW, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(Nmic_ugN_gDW)) %>%
  ungroup()

sum.multipleRRTN.Nmicend$landuse <- ordered(sum.multipleRRTN.Nmicend$landuse,  
                                            levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.Nmicend <- merge(sum.multipleRRTN.Nmicend, letters.cld.Nmicend, by = c("RR","landuse"))


#### mass specific growth rate ####
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.massspec <- multipleRRTN %>%
  group_by(RR, landuse) %>%
  summarise(n = length(mass_specific_growth_rate_1perd),
            mean = mean(mass_specific_growth_rate_1perd),
            sd = sd(mass_specific_growth_rate_1perd),
            se = sd(mass_specific_growth_rate_1perd)/sqrt(length(mass_specific_growth_rate_1perd)),
            P = t.test(mass_specific_growth_rate_1perd, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(mass_specific_growth_rate_1perd)) %>%
  ungroup()

sum.multipleRRTN.massspec$landuse <- ordered(sum.multipleRRTN.massspec$landuse,  
                                             levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.massspec <- merge(sum.multipleRRTN.massspec, letters.cld.massspec, by = c("RR","landuse"))


#### turnover ####
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.turnover <- multipleRRTN %>%
  group_by(RR, landuse) %>%
  summarise(n = length(turnover_d),
            mean = mean(turnover_d),
            sd = sd(turnover_d),
            se = sd(turnover_d)/sqrt(length(turnover_d)),
            P = t.test(turnover_d, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(turnover_d)) %>%
  ungroup()

sum.multipleRRTN.turnover$landuse <- ordered(sum.multipleRRTN.turnover$landuse,  
                                             levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.turnover <- merge(sum.multipleRRTN.turnover, letters.cld.turnover, by = c("RR","landuse"))

t.test(multipleRRT$turnover_d, mu=1)
t.test(multipleRRN$turnover_d, mu=1)


#### extractable C at the end of incubation nF ####
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.extractableCnF <- multipleRRTN %>%
  group_by(RR, landuse) %>%
  summarise(n = length(nF_ugC_gDW),
            mean = mean(nF_ugC_gDW),
            sd = sd(nF_ugC_gDW),
            se = sd(nF_ugC_gDW)/sqrt(length(nF_ugC_gDW)),
            P = t.test(nF_ugC_gDW, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(nF_ugC_gDW)) %>%
  ungroup()

sum.multipleRRTN.extractableCnF$landuse <- ordered(sum.multipleRRTN.extractableCnF$landuse,  
                                                   levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.extractableCnF <- merge(sum.multipleRRTN.extractableCnF, letters.cld.extractableCnF, by = c("RR","landuse"))


#### extractable N at the end of incubation nF ####          
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.extractableNnF <- multipleRRTN %>%
  group_by(RR, landuse) %>%
  summarise(n = length(nF_ugN_gDW),
            mean = mean(nF_ugN_gDW),
            sd = sd(nF_ugN_gDW),
            se = sd(nF_ugN_gDW)/sqrt(length(nF_ugN_gDW)),
            P = t.test(nF_ugN_gDW, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(nF_ugN_gDW)) %>%
  ungroup()

sum.multipleRRTN.extractableNnF$landuse <- ordered(sum.multipleRRTN.extractableNnF$landuse,  
                                                   levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.extractableNnF <- merge(sum.multipleRRTN.extractableNnF, letters.cld.extractableNnF, by = c("RR","landuse"))


#### Cmic:Nmic at end of incubation ####
# summarise RRN and RRT to prepare plotting
sum.multipleRRTN.CNmic <- multipleRRTN %>%
  group_by(RR, landuse) %>%
  summarise(n = length(CN_mic),
            mean = mean(CN_mic),
            sd = sd(CN_mic),
            se = sd(CN_mic)/sqrt(length(CN_mic)),
            P = t.test(CN_mic, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", "ns"),
            MaxWidth = max(CN_mic)) %>%
  ungroup()

sum.multipleRRTN.CNmic$landuse <- ordered(sum.multipleRRTN.CNmic$landuse,  
                                          levels = c("Forest", "Grassland", "Cropland"))

# add cld letters from emmeans analysis of lmem
sum.multipleRRTN.CNmic <- merge(sum.multipleRRTN.CNmic, letters.cld.CNmic, by = c("RR","landuse"))



#### plotting several RRs in same plot ####
t1 <- sum.multipleRRTN.cumulresp
t1$cat <- c("cumulative respiration")
t2 <- sum.multipleRRTN.Cresp
t2$cat <- c("C respiration")
t3 <- sum.multipleRRTN.CUE
t3$cat <- c("CUE")
t4 <- sum.multipleRRTN.Cgrowth
t4$cat <- c("C growth")
t5 <- sum.multipleRRTN.Cmicend
t5$cat <- c("microbial biomass C")
t6 <- sum.multipleRRTN.massspec
t6$cat <- c("mass specific growth rate")
t7 <- sum.multipleRRTN.turnover
t7$cat <- c("turnover")
t8 <- sum.multipleRRTN.extractableCnF
t8$cat <- c("extractable C nF")
t9 <- sum.multipleRRTN.Nmicend
t9$cat <- c("microbial biomass N")
t10 <- sum.multipleRRTN.extractableNnF
t10$cat <- c("extractable N nF")

# combine parameters, which shall be plotted together and define order of parameters
CUE.fig <- rbind(t2, t3, t4)
CUE.fig$cat <- ordered(CUE.fig$cat, levels=c("C respiration","C growth","CUE"))

extract.fig <- rbind(t8, t10)
turnover.fig <- t7 

incub.fig <- rbind(t1, t5)
incub.fig$cat <- ordered(incub.fig$cat, levels=c("cumulative respiration","microbial biomass C"))


#### Figure 4 - RRT and RRN #### 
# combine factors, which shall be plotted together and define order of parameters
RR.fig.all <- rbind(t1,t5,t2,t3,t4)
RR.fig.all$cat <- ordered(RR.fig.all$cat, levels=c("cumulative respiration","microbial biomass C", "C respiration","C growth","CUE"))
RRT.fig.all <- RR.fig.all[RR.fig.all=="T",]
RRN.fig.all <- RR.fig.all[RR.fig.all=="N",]


# plot RRT for all params (cumulative respiration, Cmic, Cresp, Cgrowth, CUE)
RRT.fig.plot <- ggplot(RRT.fig.all, aes(x= cat, y=mean, fill=RR)) +
  geom_pointrange(aes(ymin = mean-sd, ymax = mean+sd, group=landuse, shape = landuse), 
                  position = position_dodge(0.6), color = "grey35") +
  theme_bw(base_size = 13) +
  geom_hline(aes(yintercept = 1.0), colour = "#999999", linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(RR)) +
  scale_fill_manual(values = c("#EFC000FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  ylab(expression(RR[T])) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,5.5)) +
  geom_text(data = RRT.fig.all, aes(group = landuse, x = cat, y = 0, label = group),
            position = position_dodge(width = 0.6), size = 4) +
  geom_text(aes(group = landuse, label = Sig, y = mean+sd + 0.4), size = 4, position = position_dodge(width = 0.6),
            data = RRT.fig.all)
RRT.fig.plot

# plot RRT for all params (cumulative respiration, Cmic, Cresp, Cgrowth, CUE)
RRN.fig.plot <- ggplot(RRN.fig.all, aes(x= cat, y=mean, fill=RR)) +
  geom_pointrange(aes(ymin = mean-sd, ymax = mean+sd, group=landuse, shape = landuse), 
                  position = position_dodge(0.6), color = "grey35") +
  theme_bw(base_size = 13) +
  geom_hline(aes(yintercept = 1.0), colour = "#999999", linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(RR)) +
  scale_fill_manual(values = c("#0073C2FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  ylab(expression(RR[N])) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,2.5)) +
  geom_text(data = RRN.fig.all, aes(group = landuse, x = cat, y = 0, label = group),
            position = position_dodge(width = 0.6), size = 4) +
  geom_text(aes(group = landuse, label = Sig, y = mean+sd + 0.4), size = 4, position = position_dodge(width = 0.6),
            data = RRN.fig.all)
RRN.fig.plot

legend <- get_legend(RRT.fig.plot + theme(legend.position = "bottom") + guides(fill = "none"))
plotA <- plot_grid(RRT.fig.plot +
                     theme(legend.position = "none",
                           axis.title.y = element_text(size=13),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank()),
                   RRN.fig.plot +
                     theme(legend.position = "none",
                           axis.title.y = element_text(size=13),
                           axis.text.x = element_text(size=13)),    
                   labels = c("a","b"),
                   vjust = 1,
                   ncol = 1,
                   nrow = 2,
                   align = "hv", 
                   axis = "rl",
                   rel_widths = c(1,1.3))
plot_grid(plotA,
          legend,
          nrow = 2,
          rel_heights = c(1, 0.2))


#### Figure S2 - RRT and RRN of microbial turnover ####
# combine factors, which shall be plotted together and define order of parameters
RR.fig.turnover <- t7
RRT.fig.turnover <- RR.fig.turnover[RR.fig.turnover$RR=="T",]
RRN.fig.turnover <- RR.fig.turnover[RR.fig.turnover$RR=="N",]


# plot RRT for all params (cumulative respiration, Cmic, Cresp, Cgrowth, CUE)
RRT.fig.turnover <- ggplot(RRT.fig.turnover, aes(x= cat, y=mean, fill=RR)) +
  geom_pointrange(aes(ymin = mean-sd, ymax = mean+sd, group=landuse, shape = landuse), 
                  position = position_dodge(0.6), color = "grey35") +
  theme_bw(base_size = 13) +
  geom_hline(aes(yintercept = 1.0), colour = "#999999", linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(RR)) +
  scale_fill_manual(values = c("#EFC000FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  ylab(expression(RR[T])) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,2.5)) +
  geom_text(data = RRT.fig.turnover, aes(group = landuse, x = cat, y = 0, label = group),
            position = position_dodge(width = 0.6), size = 4) +
  geom_text(aes(group = landuse, label = Sig, y = mean+sd + 0.4), size = 4, position = position_dodge(width = 0.6),
            data = RRT.fig.turnover)
RRT.fig.turnover

# plot RRT for all params (cumulative respiration, Cmic, Cresp, Cgrowth, CUE)
RRN.fig.turnover <- ggplot(RRN.fig.turnover, aes(x= cat, y=mean, fill=RR)) +
  geom_pointrange(aes(ymin = mean-sd, ymax = mean+sd, group=landuse, shape = landuse), 
                  position = position_dodge(0.6), color = "grey35") +
  theme_bw(base_size = 13) +
  geom_hline(aes(yintercept = 1.0), colour = "#999999", linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(RR)) +
  scale_fill_manual(values = c("#0073C2FF")) +
  scale_shape_manual(values = c(24,21,22)) +
  ylab(expression(RR[N])) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,2.5)) +
  geom_text(data = RRN.fig.turnover, aes(group = landuse, x = cat, y = 0, label = group),
            position = position_dodge(width = 0.6), size = 4) +
  geom_text(aes(group = landuse, label = Sig, y = mean+sd + 0.4), size = 4, position = position_dodge(width = 0.6),
            data = RRN.fig.turnover)
RRN.fig.turnover

legend <- get_legend(RRN.fig.turnover + theme(legend.position = "bottom") + guides(fill = "none"))
plotA <- plot_grid(RRT.fig.turnover +
                     theme(legend.position = "none",
                           axis.title.y = element_text(size=13),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank()),
                   RRN.fig.turnover +
                     theme(legend.position = "none",
                           axis.title.y = element_text(size=13),
                           axis.text.x = element_text(size=13)),    
                   labels = c("a","b"),
                   vjust = 1,
                   ncol = 1,
                   nrow = 2,
                   align = "hv", 
                   axis = "rl",
                   rel_widths = c(1,1.3))
plot_grid(plotA,
          legend,
          nrow = 2,
          rel_heights = c(1, 0.2))
