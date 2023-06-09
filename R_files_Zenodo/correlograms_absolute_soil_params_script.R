### dependence of absolute values of microbial metabolism on soil parameters 

# Script to create correlogram of absolute microbial parameters and general soil parameters.

library(tidyverse)
library(RColorBrewer)
library(ggpmisc) # needed to add fromula to plot (stat_poly_eq_())
library(ggpubr) # can be used to combine and arrange multiple plots
                # needed to add stat_cor to plot
library(corrplot) # to plot correlogram
library(cowplot)
#--------------------------

# load script, which loads needed data
source("load_data_script.R")


#### create data.frames ####

# gensoil       data available for each sample_ID (site, landuse, replicate) = 27 values
# compsoil      data available for each plot 

#### meansdallincub_compsoil ####
# create data frame of mean and sd values for selected params of incub all data per landuse x site (including all treatments)
meansdallincub.slu <- sample_data %>%
  group_by(site, landuse) %>%
  summarise_if(is.numeric, list(mean=mean)) %>%
  ungroup()

meansdallincub_compsoil <- merge(meansdallincub.slu, compsoil) 

#### meansdallincub_gensoil ####
# create data frame of mean and sd values for selected params of incub all data per landuse x site (including all treatments)
meansdallincub <- sample_data %>%
  group_by(site, landuse, replicate) %>%
  summarise_if(is.numeric, list(mean=mean)) %>%
  ungroup()

meansdallincub_gensoil <- merge(meansdallincub, gensoil) # n = 27

#### gensoil_compsoil ####
# create mean of meansdallincub_gensoil to have n = 9 also for gensoil data
meansdallincub_gensoil_mean <- meansdallincub_gensoil %>% group_by(site, landuse) %>%
  summarise_if(is.numeric, list(mean=mean)) %>% 
  ungroup()

gensoil_compsoil <- merge(meansdallincub_gensoil_mean, compsoil, by=c("site","landuse"))

## overview new data.frames
# meansdallincub_compsoil     # n = 9 (incl. fractionation data for n = 5)
# meansdallincub_gensoil      # n = 27
# gensoil_compsoil            # n = 9 (incl. mean of gensoil data)



#### correlograms ####

# create color palette
col.jco <- colorRampPalette(c("#EFC000FF", "#0073C2FF"))
my.col.all <- col.jco(50)

#-------------------------
## correlogram - all_incub 
#-------------------------
# define parameters to be included in correlogram
params.allincub <- c("CUE", "Cgrowth_ngC_gDW_h", "Cresp_CUE_ngC_gDW_h", "mass_specific_growth_rate_1perd",
            "turnover_d", "cumul_resp_incub_ngC_g_total_time", "Cmic_ugC_gDW", "nF_ugC_gDW", 
            "Nmic_ugN_gDW","nF_ugN_gDW", "CN_mic")

# calculate correlation coefficient
corallincub <- cor(sample_data[params.allincub], method = "spearman") 

# change row and column names
colnames(corallincub) <- rownames(corallincub) <- c("CUE",
                                                    "Cgrowth",
                                                    "Cresp", 
                                                    "mass.spec",
                                                    "turnover",
                                                    "cumulresp",
                                                    "Cmic.end", 
                                                    "extr.C.nF", 
                                                    "Nmic.end", 
                                                    "extr.N.nF", 
                                                    "Cmic.Nmic.end") 

# calculate p-values of given matrix
corallincub.p <- cor.mtest(sample_data[params.allincub], method = "spearman", exact = FALSE)

# plot correlogram including information if p-values are significant
corrplot(corallincub, 
         method = "circle",
         col = my.col.all,
         bg = "white",
         tl.col = "black",
         tl.pos = "lt",
         type = "upper",
         diag = TRUE, 
         p.mat= corallincub.p$p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black") 
corrplot(corallincub, 
         add = TRUE,
         method = "number",
         number.cex = .7,
         col = "black",
         cl.pos = "n",
         bg = "white",
         tl.pos = "n",
         type = "lower",
         diag = TRUE,
         p.mat= corallincub.p$p,
         insig = "blank")
# export as 600x600

#------------------------------
## correlogram - compsoil (n=9)
#------------------------------
# define parameters to be included in correlogram
params.compsoil <- c("CUE_mean", "Cgrowth_ngC_gDW_h_mean", "Cresp_CUE_ngC_gDW_h_mean", "mass_specific_growth_rate_1perd_mean",
                     "turnover_d_mean", "cumul_resp_incub_ngC_g_total_time_mean", "Cmic_ugC_gDW_mean", "nF_ugC_gDW_mean", 
                     "Nmic_ugN_gDW_mean","nF_ugN_gDW_mean", "CN_mic_mean", "P_perc", "CP_ratio", "NP_ratio", "WHC", 
                     "sand_perc", "silt_perc", "clay_perc", "pH_1.5wvH2O", "EC_H2O_uS_cm")


# calculate correlation coefficient
cormeansincubcompsoil <- cor(meansdallincub_compsoil[params.compsoil], method = "spearman")

# change row and column names
colnames(cormeansincubcompsoil) <- rownames(cormeansincubcompsoil) <- c("CUE", 
                                                                        "Cgrowth", 
                                                                        "Cresp",
                                                                        "mass.spec", 
                                                                        "turnover",
                                                                        "cumulresp",
                                                                        "Cmic.end",
                                                                        "extr.C.nF",
                                                                        "Nmic.end", 
                                                                        "extr.N.nF", 
                                                                        "Cmic.Nmic.end", 
                                                                        "OlsenP", 
                                                                        "C:P", 
                                                                        "N:P",
                                                                        "WHC",
                                                                        "sand",
                                                                        "silt",
                                                                        "clay",
                                                                        "pH",
                                                                        "EC")

# calculate p-values of given matrix
cormeansincubcompsoil.p <- cor.mtest(meansdallincub_compsoil[params.compsoil],
                                     method = "spearman", exact = FALSE)

# plot correlogram including information if p-values are significant
corrplot(cormeansincubcompsoil, 
         method = "circle",
         col = my.col.all,
         bg = "white",
         tl.col = "black",
         tl.pos = "lt",
         type = "upper",
         diag = TRUE, 
         p.mat= cormeansincubcompsoil.p$p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black") 
corrplot(cormeansincubcompsoil, 
         add = TRUE,
         method = "number",
         number.cex = .7,
         col = "black",
         cl.pos = "n",
         bg = "white",
         tl.pos = "n",
         type = "lower",
         diag = TRUE,
         p.mat= cormeansincubcompsoil.p$p,
         insig = "blank")



#-------------------------------
## correlogram - gensoil  (n=27)
#-------------------------------
# define parameters to be included in correlogram
params.gensoil <- c("CUE_mean", "Cgrowth_ngC_gDW_h_mean", "Cresp_CUE_ngC_gDW_h_mean", "mass_specific_growth_rate_1perd_mean",
                    "turnover_d_mean", "cumul_resp_incub_ngC_g_total_time_mean", "Cmic_ugC_gDW_mean", "nF_ugC_gDW_mean", 
                    "Nmic_ugN_gDW_mean","nF_ugN_gDW_mean", "CN_mic_mean", "B_cn_gsoil", "A_cn_gsoil", "F_cn_gsoil", "B.A", 
                    "F.A", "F.B", "Cmic_V1_ugC_gsoil", "Nmic_V1_ugN_gsoil", "Cmic.Nmic", "Cmic_Corg_perc",
                    "Corg_perc", "Ntotal_perc", "CN_ratio")


# calculate correlation coefficient
cormeansincubgensoil <- cor(meansdallincub_gensoil[params.gensoil], method = "spearman")

# change row and column names
colnames(cormeansincubgensoil) <- rownames(cormeansincubgensoil) <- c("CUE",
                                                                      "Cgrowth", 
                                                                      "Cresp",
                                                                      "mass.spec", 
                                                                      "turnover",
                                                                      "cumulresp",
                                                                      "Cmic.end", 
                                                                      "extr.C.nF",  
                                                                      "Nmic.end", 
                                                                      "extr.N.nF", 
                                                                      "Cmic.Nmic.end", 
                                                                      "Bacteria", 
                                                                      "Archaea", 
                                                                      "Fungi",
                                                                      "B:A",
                                                                      "F:A",
                                                                      "F:B",
                                                                      "Cmic.start", 
                                                                      "Nmic.start", 
                                                                      "Cmic:Nmic", 
                                                                      "Cmic/Corg", 
                                                                      "Corg", 
                                                                      "Ntotal", 
                                                                      "C:N")



# calculate p-values of given matrix
cormeansincubgensoil.p <- cor.mtest(meansdallincub_gensoil[params.gensoil], 
                                    method = "spearman", exact = FALSE)

# plot correlogram including information if p-values are significant
corrplot(cormeansincubgensoil, 
         method = "circle",
         col = my.col.all,
         bg = "white",
         tl.col = "black",
         tl.pos = "lt",
         type = "upper",
         diag = TRUE, 
         p.mat= cormeansincubgensoil.p$p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black") 
corrplot(cormeansincubgensoil, 
         add = TRUE,
         method = "number",
         number.cex = .7,
         col = "black",
         cl.pos = "n",
         bg = "white",
         tl.pos = "n",
         type = "lower",
         diag = TRUE,
         p.mat= cormeansincubgensoil.p$p,
         insig = "blank")
# export as 800x800

#----------------------------------
## correlogram - gensoil - compsoil
#----------------------------------
# define parameters to be included in correlogram
params.gensoil.compsoil <- c("CUE_mean_mean", "Cgrowth_ngC_gDW_h_mean_mean", "Cresp_CUE_ngC_gDW_h_mean_mean", "mass_specific_growth_rate_1perd_mean_mean",
                             "turnover_d_mean_mean", "cumul_resp_incub_ngC_g_total_time_mean_mean", "Cmic_ugC_gDW_mean_mean", "nF_ugC_gDW_mean_mean", 
                             "Nmic_ugN_gDW_mean_mean","nF_ugN_gDW_mean_mean", "CN_mic_mean_mean", "B_cn_gsoil_mean", "A_cn_gsoil_mean", "F_cn_gsoil_mean", "B.A_mean", 
                             "F.A_mean", "F.B_mean", "Cmic_V1_ugC_gsoil_mean", "Nmic_V1_ugN_gsoil_mean", "Cmic.Nmic_mean", "Cmic_Corg_perc_mean",
                             "TOC_perc", "N_perc", "CN_ratio",
                             "P_perc", "CP_ratio", "NP_ratio", "WHC", 
                             "sand_perc", "silt_perc", "clay_perc", "pH_1.5wvH2O", "EC_H2O_uS_cm")

# calculate correlation coefficient
corgensoilcompsoil <- cor(gensoil_compsoil[params.gensoil.compsoil], 
                          method = "spearman")

colnames(corgensoilcompsoil) <- rownames(corgensoilcompsoil) <- c("CUE",
                                                                  "Cgrowth", 
                                                                  "Cresp",
                                                                  "mass.spec", 
                                                                  "turnover",
                                                                  "cumulresp",
                                                                  "Cmic.end", 
                                                                  "extr.C.nF",  
                                                                  "Nmic.end", 
                                                                  "extr.N.nF",
                                                                  "Cmic.Nmic.end", 
                                                                  "Bacteria", 
                                                                  "Archaea", 
                                                                  "Fungi",
                                                                  "B:A",
                                                                  "F:A",
                                                                  "F:B",
                                                                  "Cmic.start", 
                                                                  "Nmic.start", 
                                                                  "Cmic:Nmic", 
                                                                  "Cmic/Corg", 
                                                                  "Corg", 
                                                                  "Ntotal", 
                                                                  "C:N", 
                                                                  "OlsenP", 
                                                                  "C:P", 
                                                                  "N:P",
                                                                  "WHC",
                                                                  "sand",
                                                                  "silt",
                                                                  "clay",
                                                                  "pH",
                                                                  "EC")

# calculate p-values of given matrix
corgensoilcompsoil.p <- cor.mtest(gensoil_compsoil[params.gensoil.compsoil], 
                                  method = "spearman", exact = FALSE)

# plot correlogram including information if p-values are significant
corrplot(corgensoilcompsoil, 
         method = "circle",
         col = my.col.all,
         bg = "white",
         tl.col = "black",
         tl.pos = "lt",
         type = "upper",
         diag = TRUE, 
         p.mat= corgensoilcompsoil.p$p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black") 
corrplot(corgensoilcompsoil, 
         add = TRUE,
         method = "number",
         number.cex = .7,
         col = "black",
         cl.pos = "n",
         bg = "white",
         tl.pos = "n",
         type = "lower",
         diag = TRUE,
         p.mat= corgensoilcompsoil.p$p,
         insig = "blank")
# export as 800x800



#### Figure S1 - correlogram absolute values ####

## combine matrices to get one correlogram for all data pair combinations
# only works if same parameters are at same position. Check!

# corgensoilcompsoil  n=9 - 1:30
# cormeansincubgensoil n=27 - 1:21
# corallincub n=81 - 1:8

mat <- corgensoilcompsoil
sub_mat_1 <- cormeansincubgensoil
sub_mat_2 <- corallincub

mat[1:24,1:24] <- sub_mat_1 #1:21 without Nmic
mat[1:11,1:11] <- sub_mat_2 #1:8 without Nmic


mat.p <- corgensoilcompsoil.p$p
sub_mat_1.p <- cormeansincubgensoil.p$p
sub_mat_2.p <- corallincub.p$p

mat.p[1:24,1:24] <- sub_mat_1.p
mat.p[1:11,1:11] <- sub_mat_2.p

mat.without.Nmic.end <- mat[c(1:7,12:33),c(1:7,12:33)]
mat.p.without.Nmic.end <- mat.p[c(1:7,12:33),c(1:7,12:33)]

# plot correlogram including information if p-values are significant
corrplot(mat.without.Nmic.end, 
         method = "circle",
         col = my.col.all,
         bg = "white",
         tl.col = "black",
         tl.pos = "lt",
         type = "upper",
         diag = TRUE, 
         p.mat= mat.p.without.Nmic.end,
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black") 
corrplot(mat.without.Nmic.end, 
         add = TRUE,
         method = "number",
         number.cex = .7,
         col = "black",
         cl.pos = "n",
         bg = "white",
         tl.pos = "n",
         type = "lower",
         diag = TRUE,
         p.mat= mat.p.without.Nmic.end,
         insig = "blank")
#export 1500x1500