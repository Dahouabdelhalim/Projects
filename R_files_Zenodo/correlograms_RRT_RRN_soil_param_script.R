### dependence of temperature sensitivity and nitrogen response on soil parameters 

# Script to create correlogram of response ratios and general soil parameters.

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

# load necessary scripts creating multiple RRT0 and RRN data.frames
source("RRT_RRN_calculation_script.R")


#### create data.frame ####


# create data.frame which contains RRT, RRN and general soil parameters as well as qPCR data
gensoilRRTN <- merge(gensoil, multipleRRTN, by = c("site", "landuse", "replicate"))

# create two sub.data.frames for correlation analysis
gensoilRRTN.sub.N <- gensoilRRTN %>% filter(RR=="N")
gensoilRRTN.sub.T <- gensoilRRTN %>% filter(RR=="T")


# create combined data.frame with general soil parameters on composite samples 
# and mean RRT and RRN values

# summarise three field replicates of multipleRRTN
multipleRRTN.meansd <- multipleRRTN %>% 
  group_by(site, landuse, RR) %>%
  summarise_if(is.numeric, list(mean = mean, sd=sd)) %>%
  ungroup()

multipleRRTN.meansd.sub.T <- multipleRRTN.meansd %>% filter(RR=="T")
multipleRRTN.meansd.sub.N <- multipleRRTN.meansd %>% filter(RR=="N")

compsoilmeanRRT <- merge(compsoil, multipleRRTN.meansd.sub.T, by = c("site", "landuse"))  
compsoilmeanRRN <- merge(compsoil, multipleRRTN.meansd.sub.N, by = c("site", "landuse"))  

compsoilmeanRRTN <- merge(compsoil, multipleRRTN.meansd, by = c("site", "landuse"))  

# summarise three field replicates of gensoil
meangensoil <- gensoil %>% 
  group_by(site, landuse) %>%
  summarise_if(is.numeric, list(mean=mean, sd=sd)) %>%
  ungroup()

compsoilmeanRRTNmeangensoil <- merge(compsoilmeanRRTN, meangensoil,  by=c("site","landuse"))

compsoilmeanRRTNmeangensoil.subT <- compsoilmeanRRTNmeangensoil[compsoilmeanRRTNmeangensoil$RR=="T",]
compsoilmeanRRTNmeangensoil.subN <- compsoilmeanRRTNmeangensoil[compsoilmeanRRTNmeangensoil$RR=="N",]


#### correlograms for RRT and RRN ####
# create color palettes
col.jco.yellow <- colorRampPalette(c("#999999","white", "#EFC000FF")) # grey to yellow (= T treatment)
my.col.RRT <- col.jco.yellow(50)

col.jco.blue <- colorRampPalette(c("#999999","white", "#0073C2FF"))  # grey to blue (= N treatment)
my.col.RRN <- col.jco.blue(50)


### RRT and RRN - gensoil - per sample ####

# define parameters to be included in correlogram
params.RRTN.gensoil <- c("site", "landuse","CUE", "Cgrowth_ngC_gDW_h", "Cresp_CUE_ngC_gDW_h", "mass_specific_growth_rate_1perd",
                         "turnover_d", "cumul_resp_incub_ngC_g_total_time", "Cmic_ugC_gDW", "nF_ugC_gDW", 
                         "Nmic_ugN_gDW","nF_ugN_gDW", "CN_mic", "B_cn_gsoil", "A_cn_gsoil", "F_cn_gsoil", "B.A", 
                         "F.A", "F.B", "Cmic_V1_ugC_gsoil", "Nmic_V1_ugN_gsoil", "Cmic.Nmic", "Cmic_Corg_perc",
                         "Corg_perc", "Ntotal_perc", "CN_ratio")


#### correlogram - RRT - mean per sample ####

sub.T <- gensoilRRTN.sub.T[params.RRTN.gensoil]
colnames(sub.T) <- c("site", "landuse",
                     "RRT_CUE",
                     "RRT_Cgrowth", 
                     "RRT_Cresp",
                     "RRT_mass.spec", 
                     "RRT_turnover",
                     "RRT_cumulresp",
                     "RRT_Cmic.end", 
                     "RRT_extr.C.nF", 
                     "RRT_Nmic.end", 
                     "RRT_extr.N.nF", 
                     "RRT_Cmic.Nmic.end", 
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

gensoilRRT.matrix.cor <- cor(sub.T[c(3:10,
                                     11:13, # Nmic.end, extr.N.nF, Cmic.Nmic.end
                                     14:20,
                                     21:22, # Nmic.start, Cmic:Nmic
                                     23:26)], method = "spearman")
gensoilRRT.matrix.cor.p <- cor.mtest(sub.T[c(3:10,11:13,14:20,21:22,23:26)], method = "spearman",exact = FALSE)

# plot correlogram including information if p-values are significant
corrplot(gensoilRRT.matrix.cor[1:11,], 
         method = "circle",
         col = my.col.RRT,
         bg = "white",
         tl.col = "black",
         tl.pos = "td",
         type = "upper",
         diag = TRUE, 
         p.mat= gensoilRRT.matrix.cor.p$p[1:11,],
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black") 


#### correlogram - RRN - mean per sample ####

sub.N <- gensoilRRTN.sub.N[params.RRTN.gensoil]
colnames(sub.N) <- c("site", "landuse",
                     "RRN_CUE",
                     "RRN_Cgrowth", 
                     "RRN_Cresp",
                     "RRN_mass.spec", 
                     "RRN_turnover",
                     "RRN_cumulresp",
                     "RRN_Cmic.end", 
                     "RRN_extr.C.nF", 
                     "RRN_Nmic.end", 
                     "RRN_extr.N.nF", 
                     "RRN_Cmic.Nmic.end",
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

gensoilRRN.matrix.cor <- cor(sub.N[c(3:10,
                                     11:13, # Nmic.end, extr.N.nF, Cmic.Nmic.end
                                     14:20,
                                     21:22, # Nmic.start, Cmic:Nmic
                                     23:26)], method = "spearman")
gensoilRRN.matrix.cor.p <- cor.mtest(sub.N[c(3:10,11:13,14:20,21:22,23:26)], method = "spearman",exact = FALSE)


# plot correlogram including information if p-values are significant
corrplot(gensoilRRN.matrix.cor[1:11,], 
         method = "circle",
         col = my.col.RRN,
         bg = "white",
         tl.col = "black",
         tl.pos = "td",
         type = "upper",
         diag = TRUE, 
         p.mat= gensoilRRN.matrix.cor.p$p[1:11,],
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black") 



#### RRT and RRN - gensoil - compsoil - per plot ####

# correlation analysis for mean values (n=9) of RRT and RRN and gensoil data with compsoil data
# including gensoil data is needed so that smaller correlogram of RRT and RRN with gensoil (n=27) can fit into bigger matrix

# define parameters to be included in correlogram
params.RRTN.gensoil.compsoil <- c("CUE_mean", "Cgrowth_ngC_gDW_h_mean", "Cresp_CUE_ngC_gDW_h_mean", "mass_specific_growth_rate_1perd_mean",
                                  "turnover_d_mean", "cumul_resp_incub_ngC_g_total_time_mean", "Cmic_ugC_gDW_mean", "nF_ugC_gDW_mean", 
                                  "Nmic_ugN_gDW_mean","nF_ugN_gDW_mean", "CN_mic_mean", "B_cn_gsoil_mean", "A_cn_gsoil_mean", "F_cn_gsoil_mean", "B.A_mean", 
                                  "F.A_mean", "F.B_mean", "Cmic_V1_ugC_gsoil_mean", "Nmic_V1_ugN_gsoil_mean", "Cmic.Nmic_mean", "Cmic_Corg_perc_mean",
                                  "Corg_perc_mean", "Ntotal_perc_mean", "CN_ratio_mean",
                                  "P_perc", "CP_ratio", "NP_ratio", "WHC", 
                                  "sand_perc", "silt_perc", "clay_perc", "pH_1.5wvH2O")

#### correlogram - RRT - mean per plot ####
# calculate correlation coefficient
corRRTgensoilcompsoil <- cor(compsoilmeanRRTNmeangensoil.subT[params.RRTN.gensoil.compsoil], method = "spearman")

colnames(corRRTgensoilcompsoil) <- rownames(corRRTgensoilcompsoil) <- c("RRT_CUE",
                                                                        "RRT_Cgrowth", 
                                                                        "RRT_Cresp",
                                                                        "RRT_mass.spec", 
                                                                        "RRT_turnover",
                                                                        "RRT_cumulresp",
                                                                        "RRT_Cmic.end", 
                                                                        "RRT_extr.C.nF", 
                                                                        "RRT_Nmic.end", 
                                                                        "RRT_extr.N.nF", 
                                                                        "RRT_Cmic.Nmic.end", 
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
                                                                        "pH")
                                                                        
                                                                   



# calculate p-values of given matrix
corRRTgensoilcompsoil.p <- cor.mtest(compsoilmeanRRTNmeangensoil.subT[params.RRTN.gensoil.compsoil], method = "spearman",exact = FALSE)

# plot correlogram including information if p-values are significant
corrplot(corRRTgensoilcompsoil[1:11,], 
         method = "circle",
         col = my.col.RRT,
         bg = "white",
         tl.col = "black",
         tl.pos = "td",
         type = "upper",
         diag = TRUE, 
         p.mat= corRRTgensoilcompsoil.p$p[1:11,],
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black") 


#### correlogram - RRN - mean per plot ####
# calculate correlation coefficient
corRRNgensoilcompsoil <- cor(compsoilmeanRRTNmeangensoil.subN[params.RRTN.gensoil.compsoil], method = "spearman")

colnames(corRRNgensoilcompsoil) <- rownames(corRRNgensoilcompsoil) <- c("RRN_CUE",
                                                                        "RRN_Cgrowth", 
                                                                        "RRN_Cresp",
                                                                        "RRN_mass.spec", 
                                                                        "RRN_turnover",
                                                                        "RRN_cumulresp",
                                                                        "RRN_Cmic.end", 
                                                                        "RRN_extr.C.nF", 
                                                                        "RRN_Nmic.end", 
                                                                        "RRN_extr.N.nF", 
                                                                        "RRN_Cmic.Nmic.end",
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
                                                                        "pH")


# calculate p-values of given matrix
corRRNgensoilcompsoil.p <- cor.mtest(compsoilmeanRRTNmeangensoil.subN[params.RRTN.gensoil.compsoil], 
                                     method = "spearman",exact = FALSE)


# plot correlogram including information if p-values are significant
corrplot(corRRNgensoilcompsoil[1:11,], 
         method = "circle",
         col = my.col.RRN,
         bg = "white",
         tl.col = "black",
         tl.pos = "td",
         type = "upper",
         diag = TRUE, 
         p.mat= corRRNgensoilcompsoil.p$p[1:11,],
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black") 



# Figure 5 - RRT and RRN correlograms ####

#### RRT ####
# combine correlation table for n=27 and n=9
# without Nmic data
c <- gensoilRRT.matrix.cor[c(1:7),c(1:7,12:24)]
c.p <- gensoilRRT.matrix.cor.p$p[c(1:7),c(1:7,12:24)]

d <- corRRTgensoilcompsoil[c(1:7),25:32]
d.p <- corRRTgensoilcompsoil.p$p[c(1:7),25:32]

correloRRTall <- cbind(c,d)
correloRRTall.p <- cbind(c.p, d.p)

# plot correlogram including information if p-values are significant
corrplot(correloRRTall, 
         method = "circle",
         col = my.col.RRT,
         bg = "white",
         tl.col = "black",
         tl.pos = "td",
         type = "upper",
         diag = TRUE, 
         p.mat= correloRRTall.p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black")



#### RRN ####
# combine correlation table for n=27 and n=9
# without Nmic data
a <- gensoilRRN.matrix.cor[1:7,c(1:7,12:24)]
a.p <- gensoilRRN.matrix.cor.p$p[1:7,c(1:7,12:24)]

b <- corRRNgensoilcompsoil[1:7,25:32]
b.p <- corRRNgensoilcompsoil.p$p[1:7,25:32]


correloRRNall <- cbind(a,b)
correloRRNall.p <- cbind(a.p, b.p)


# plot correlogram including information if p-values are significant
corrplot(correloRRNall, 
         method = "circle",
         col = my.col.RRN,
         bg = "white",
         tl.col = "black",
         tl.pos = "td",
         type = "upper",
         diag = TRUE, 
         p.mat= correloRRNall.p,
         insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black")
# 1000x500