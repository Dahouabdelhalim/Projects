### 0. import data -------------------------------------------------------------
setwd("~/Documents/Post_Doc/01_Reproductive_Isolation_Review/4_Manuscript/files_for_Dryad_REVISED/")
list.files(pattern = ".csv")
dat <- read.csv(file = "RI_data_FINAL.csv", header = T, stringsAsFactors = F)


### 01. find average barrier strength for each taxa pair -----------------------

# find average barrier strength for both taxa
dat$Ecogeo_mean <- apply(dat[, c("Ecogeo1", "Ecogeo2")], MARGIN = 1, mean, na.rm = T)
dat$ImmigrantInviability_mean <- apply(dat[, c("ImmigrantInviability1", "ImmigrantInviability2")], MARGIN = 1, mean, na.rm = T)
dat$Pheno_mean <- apply(dat[, c("Pheno1", "Pheno2")], MARGIN = 1, mean, na.rm = T)

# combine "Mating System" and "Differential Pollen (Production)" columns 
dat$MatingSystem_mean <- apply(dat[, c("MatingSystem1", "MatingSystem2", "DifferentialPollen1", "DifferentialPollen2")], MARGIN = 1, mean, na.rm = T)

dat$FloralIsolation_mean <- apply(dat[, c("FloralIsolation1", "FloralIsolation2")], MARGIN = 1, mean, na.rm = T)
dat$PollenPistil_mean <- apply(dat[, c("PollenPistil1", "PollenPistil2")], MARGIN = 1, mean, na.rm = T)
dat$FruitProduction_mean <- apply(dat[, c("FruitProduction1", "FruitProduction2")], MARGIN = 1, mean, na.rm = T)
dat$SeedProduction_mean <- apply(dat[, c("SeedProduction1", "SeedProduction2")], MARGIN = 1, mean, na.rm = T)
dat$F1Germination_mean <- apply(dat[, c("F1Germination1", "F1Germination2")], MARGIN = 1, mean, na.rm = T)
dat$F1Viability_mean <- apply(dat[, c("F1Viability1", "F1Viability2")], MARGIN = 1, mean, na.rm = T)

# combine "F1 Pollen Sterility" and "F1 Ovule Fertility" columns 
dat$F1Sterility_mean <- apply(dat[, c("F1PollenSterility1", "F1PollenSterility2", "F1OvuleFertility1", "F1OvuleFertility2")], MARGIN = 1, mean, na.rm = T)

dat$ExtrinsicPost_mean <- apply(dat[, c("ExtrinsicPost1", "ExtrinsicPost2")], MARGIN = 1, mean, na.rm = T)


# add individual sterility columns (for later asymmetry analyses)
dat$F1Sterility1 <- apply(dat[, c("F1PollenSterility1", "F1OvuleFertility1")], MARGIN = 1, mean, na.rm = T)
dat$F1Sterility2 <- apply(dat[, c("F1PollenSterility2", "F1OvuleFertility2")], MARGIN = 1, mean, na.rm = T)

# add individual mating system columns (for later asymmetry analyses)
dat$MatingSystem1 <- apply(dat[, c("MatingSystem1", "DifferentialPollen1")], MARGIN = 1, mean, na.rm = T)
dat$MatingSystem2 <- apply(dat[, c("MatingSystem2", "DifferentialPollen2")], MARGIN = 1, mean, na.rm = T)


# find number of studies for each barrier
sum(!is.na(dat$Ecogeo_mean)) # 26
sum(!is.na(dat$ImmigrantInviability_mean)) # 22
sum(!is.na(dat$Pheno_mean)) # 52
sum(!is.na(dat$MatingSystem_mean)) # 10
sum(!is.na(dat$FloralIsolation_mean)) # 59
sum(!is.na(dat$PollenPistil_mean)) # 30
sum(!is.na(dat$FruitProduction_mean)) # 40
sum(!is.na(dat$SeedProduction_mean)) # 62
sum(!is.na(dat$F1Germination_mean)) # 33
sum(!is.na(dat$F1Viability_mean)) # 25
sum(!is.na(dat$F1Sterility_mean)) # 21
sum(!is.na(dat$ExtrinsicPost_mean)) # 20


# find average barrier strengths for Ecogeography and Floral Isolation subcategories
dat$Ecogeo_ecogeo_mean <- apply(dat[, c("Ecogeo1_ecogeo", "Ecogeo2_ecogeo")], MARGIN = 1, mean, na.rm = T)
dat$Ecogeo_micro_mean <- apply(dat[, c("Ecogeo1_micro", "Ecogeo2_micro")], MARGIN = 1, mean, na.rm = T)

dat$FloralIsolation_pollinators_mean <- apply(dat[, c("FloralIsolation1_pollinators", "FloralIsolation2_pollinators")], MARGIN = 1, mean, na.rm = T)
dat$FloralIsolation_transitions_mean <- apply(dat[, c("FloralIsolation1_transitions", "FloralIsolation2_transitions")], MARGIN = 1, mean, na.rm = T)
dat$FloralIsolation_deposition_mean <- apply(dat[, c("FloralIsolation1_deposition", "FloralIsolation2_deposition")], MARGIN = 1, mean, na.rm = T)

# find number of studies for each barrier
sum(!is.na(dat$Ecogeo_ecogeo_mean)) # 22
sum(!is.na(dat$Ecogeo_micro_mean)) # 6

sum(!is.na(dat$FloralIsolation_pollinators_mean)) # 25
sum(!is.na(dat$FloralIsolation_transitions_mean)) # 26
sum(!is.na(dat$FloralIsolation_deposition_mean)) # 18


### 04. format data for ggplot -------------------------------------------------
RI_df <- data.frame(
    RI_value =
    c(dat$Ecogeo_mean,
    dat$ImmigrantInviability_mean,
    dat$Pheno_mean,
    dat$MatingSystem_mean,
    dat$FloralIsolation_mean,
    dat$PollenPistil_mean,
    dat$FruitProduction_mean,
    dat$SeedProduction_mean,
    dat$F1Germination_mean,
    dat$F1Viability_mean,
    dat$F1Sterility_mean,
    dat$ExtrinsicPost_mean), 
    
    Barrier = rep(c("Ecogeographic Isolation", "Immigrant Inviability", 
                    "Phenology","Mating System", "Floral Isolation",
                    "Pollen Pistil Interactions", "Fruit Production", "Seed Production", 
                    "F1 Germination","F1 Viability", "F1 Sterility",
                    "Extrinsic Postzygotic Isolation"), each = nrow(dat)),
    
    in_Lowry = rep(c(dat$in_Lowry_2008), times = 12),
    
    geography = rep(dat$geography, times = 12),
    
    taxa_type = rep(dat$Taxa_type2, times = 12))

# 12 barriers; including 5 premating, 1 postmating prezygotic, 6 postzygotic (5 intrinsic, 1 extrinsic)
RI_df$Barrier <- factor(RI_df$Barrier, 
    levels = c("Ecogeographic Isolation", "Immigrant Inviability", 
               "Phenology","Mating System", "Floral Isolation",
               "Pollen Pistil Interactions", "Fruit Production", "Seed Production",
               "F1 Germination","F1 Viability", "F1 Sterility",
               "Extrinsic Postzygotic Isolation"))

Barrier <- factor(RI_df$Barrier, 
                  levels = c("Ecogeographic Isolation", "Immigrant Inviability", 
                             "Phenology","Mating System", "Floral Isolation",
                             "Pollen Pistil Interactions", "Fruit Production", "Seed Production",
                             "F1 Germination","F1 Viability", "F1 Sterility",
                             "Extrinsic Postzygotic Isolation"))


### 05. summarize barrier strengths (TABLE 1) ----------------------------------
tapply(RI_df$RI_value, RI_df$Barrier, mean, na.rm = T)
tapply(RI_df$RI_value, RI_df$Barrier, median, na.rm = T)
tapply(RI_df$RI_value, RI_df$Barrier, min, na.rm = T)
tapply(RI_df$RI_value, RI_df$Barrier, max, na.rm = T)

tapply(RI_df$RI_value, RI_df$Barrier, function(x)  quantile(x, probs = c(0.25, 0.75), na.rm = T))

# summarize barrier strengths for sub-categories 
# of Ecogeographic Isolation and Floral Isolation

summary(dat$Ecogeo_ecogeo_mean, na.rm = T) 
summary(dat$Ecogeo_micro_mean, na.rm = T) 

summary(dat$FloralIsolation_pollinators_mean, na.rm = T) 
summary(dat$FloralIsolation_transitions_mean, na.rm = T) 
summary(dat$FloralIsolation_deposition_mean, na.rm = T) 


###  06. test for differences in barrier strength among barriers ---------------
library(car)

m1 <- lm(RI_value ~ Barrier, data = RI_df)
Anova(m1, type = "3")

# test for pairwise differences among barriers
library(emmeans)
library(ggplot2)

m1_emm <- emmeans(m1, "Barrier")
pairs(m1_emm, adjust = "tukey")


dev.new()
# https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html#graphical
p1 <- plot(m1_emm, comparisons = T, xlab = "Barrier strength",
     colors = c("black", "black", "red", "blue"), horizontal = F)
# gray bars are confidence intervals for EMMS
# blue arrows are comparisons among them
# if arrow overlap, difference is not significant based on Tukey's adjustment

p1 + 
    theme(axis.text.x = element_text(angle = 50, vjust = 0.95, hjust = 1)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(face = "bold", size = 14)) 
    
# quartz.save(file = "FigureS1_Barrier_Comparisons_05_22_2022_REVISED.jpg", type = "jpg", dpi = 300)


# pairwise matrix (upper triangle) of p-values
pwpm(m1_emm) # matrix of pairwise p-values

# plotting customization
my.aes <- list(point = list(shape = "circle"), 
               segment = list(linetype = "solid", color = "blue"),
               label = list(family = "serif", fontface = "bold"))

my.pal <- rep("black", times = 12)

p2 <- pwpp(m1_emm, aes = my.aes, values = T) 

p2 + 
    scale_color_manual(values = my.pal) +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_text(face = "bold", size = 14)) 

# quartz.save(file = "FigureS2_Barrier_Significance_05_22_2022_REVISED.jpg", type = "jpg", dpi = 300)


# 19 pairs of barriers are significantly different after correcting for multiple comparisons
# 66 pairwise comparisons
19/66 # 28.8% of pairwise comparisons have significantly different barrier strengths


### 07. plot data (Figure 2) ---------------------------------------------------

# plot data using ggplot
library(ggplot2)

dev.new()

# add sample size to x-axis labels (custom labels)
my_labels <- c(
    "Ecogeographic Isolation (n = 26)", 
    "Immigrant Inviability (n = 22)", 
    "Phenology (n = 52)",
    "Mating System (n = 10)", 
    "Floral Isolation (n = 59)",
    "Pollen Pistil Interactions (n = 30)", 
    "Fruit Production (n = 40)", 
    "Seed Production (n = 62)",
    "F1 Germination (n = 33)",
    "F1 Viability (n = 25)", 
    "F1 Sterility (n = 21)",
    "Extrinsic Postzygotic Isolation (n = 20)"
    )


# violin plot
ggplot(data = RI_df, aes(x = Barrier, y = RI_value)) +
    geom_hline(yintercept = 0, linetype = "solid") + 
    geom_violin(trim = TRUE) +
    scale_y_continuous(limits = c(-0.25,1), breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
    scale_x_discrete(labels = my_labels) +
    geom_jitter(width = 0.075, aes(shape = geography, fill = taxa_type), size = 2, alpha = 0.7) + 
    scale_shape_manual(values = c(24,22,21), na.translate = F) +
    stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 55, vjust = 0.8, hjust = 0.75),
        axis.text = element_text(size = 11, face = "bold")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(face = "bold", size = 16)) +
    labs(y = "Barrier Strength") +
    geom_segment(mapping = aes(x = 6.5, y = -0.25, xend = 6.5, yend = 1), size = 0.5, color = "black", linetype = "dashed") +
    guides(fill = guide_legend(override.aes=list(shape=23))) +
    labs(fill = "Taxa type", shape = "Geography") +
    theme(legend.title = element_text(size = 14, face = "bold")) +
    theme(legend.text = element_text(size = 12))

# quartz.save(file = "Figure2_FINAL.jpg", type = "jpg", dpi = 600)


### 09. plot barrier strength by taxon type ------------------------------------
dev.new()

ggplot(data = RI_df, aes(x = taxa_type, y = RI_value, fill = taxa_type)) +
    geom_violin() +
    scale_fill_discrete(name = "Taxa type", 
                        labels = c(
                            "cytotypes", 
                            "ecotypes",
                            "species")) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.35, size = 0.5, color = "black", show.legend = F) +
    scale_shape_manual(values = c(24,22,21), na.translate = F) +
    geom_jitter(width = 0.075, aes(shape = geography), color = "black", fill = "black", size = 2, alpha = 0.35) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ylim(-0.5,1) +
    facet_wrap(~Barrier, nrow = 3) +
    theme_light() +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(face = "bold", size = 14)) +
    labs(y = "Barrier Strength") +
    theme(legend.title = element_text(size = 12, face = "bold")) +
    guides(fill = guide_legend(override.aes=list(shape=23))) +
    labs(fill = "Taxa type", shape = "Geography") +
    theme(legend.title = element_text(size = 10, face = "bold")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.75, vjust = 0.75))
    
# quartz.save(file = "FigureS3_barrier_strengths_by_taxon_type_05_22_2022.jpg", type = "jpg", dpi = 300)



### 10. plot RI for sub-categories of ecogeographic and floral isolation -------
# data plotted, but not shown or discussed in the manuscript 
# beyond summary statistics provided in Table 1

### Ecogeographic isolation
RI_df_ecogeo <- data.frame(
    RI_value =
        c(dat$Ecogeo_ecogeo_mean,
          dat$Ecogeo_micro_mean), 
    Barrier = rep(c("Geographic Isolation", "Microhabitat Isolation"), each = nrow(dat)),
    geography = rep(dat$geography, times = 2),
    taxa_type = rep(dat$Taxa_type2, times = 2))

# violin plot
ggplot(data = RI_df_ecogeo, aes(x = Barrier, y = RI_value)) +
    geom_hline(yintercept = 0, linetype = "solid") + 
    geom_violin(trim = TRUE) +
    scale_y_continuous(limits = c(-0.25,1), breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
    geom_jitter(width = 0.075, aes(shape = geography, fill = taxa_type), size = 2, alpha = 0.7) +
    scale_shape_manual(values = c(24,22,21), na.translate = F) +
    stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 50, vjust = 0.8, hjust = 0.75),
        axis.text = element_text(size = 9, face = "bold")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(face = "bold", size = 14)) +
    labs(y = "Barrier Strength") +
    guides(fill = guide_legend(override.aes=list(shape=23))) +
    labs(fill = "Taxa type", shape = "Geography") +
    theme(legend.title = element_text(size = 10, face = "bold")) 


### Floral isolation
RI_df_floral <- data.frame(
    RI_value =
        c(dat$FloralIsolation_pollinators_mean,
          dat$FloralIsolation_transitions_mean,
          dat$FloralIsolation_deposition_mean), 
    Barrier = rep(c("Pollinator assemblages", "Pollinator transitions", "Pollen deposition"), each = nrow(dat)),
    geography = rep(dat$geography, times = 3),
    taxa_type = rep(dat$Taxa_type2, times = 3))

# define levels of factor
RI_df_floral$Barrier <- factor(RI_df_floral$Barrier, 
                        levels = c("Pollinator assemblages", "Pollinator transitions", "Pollen deposition"))


# violin plot
ggplot(data = RI_df_floral, aes(x = Barrier, y = RI_value)) +
    geom_hline(yintercept = 0, linetype = "solid") + 
    geom_violin(trim = TRUE) +
    scale_y_continuous(limits = c(-0.25,1), breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
    geom_jitter(width = 0.075, aes(shape = geography, fill = taxa_type), size = 2, alpha = 0.7) +
    scale_shape_manual(values = c(24,22,21), na.translate = F) +
    stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 50, vjust = 0.8, hjust = 0.75),
        axis.text = element_text(size = 9, face = "bold")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(face = "bold", size = 14)) +
    labs(y = "Barrier Strength") +
    guides(fill = guide_legend(override.aes=list(shape=23))) +
    labs(fill = "Taxa type", shape = "Geography") +
    theme(legend.title = element_text(size = 10, face = "bold")) 



### 11. (for Discussion) -------------------------------------------------------
### assess correlation between immigrant inviability and extrinsic post --------

II_df <- data.frame(
    II = dat$ImmigrantInviability_mean,
    EP = dat$ExtrinsicPost_mean)


II_df <- II_df[complete.cases(II_df),]
nrow(II_df)

cor.test(II_df$II, II_df$EP) # cor = 0.58, p = 0.02















