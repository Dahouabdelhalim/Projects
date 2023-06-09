### 0. import data -------------------------------------------------------------
setwd("~/Documents/Post_Doc/01_Reproductive_Isolation_Review/4_Manuscript/files_for_Dryad_REVISED/")
list.files(pattern = ".csv")
dat <- read.csv(file = "RI_data_FINAL.csv", header = T, stringsAsFactors = F)


### 01. Calculate asymmetry ----------------------------------------------------

# add individual sterility columns
dat$F1Sterility1 <- apply(dat[, c("F1PollenSterility1", "F1OvuleFertility1")], MARGIN = 1, mean, na.rm = T)
dat$F1Sterility2 <- apply(dat[, c("F1PollenSterility2", "F1OvuleFertility2")], MARGIN = 1, mean, na.rm = T)

# add individual mating system columns
dat$MatingSystem1 <- apply(dat[, c("MatingSystem1", "DifferentialPollen1")], MARGIN = 1, mean, na.rm = T)
dat$MatingSystem2 <- apply(dat[, c("MatingSystem2", "DifferentialPollen2")], MARGIN = 1, mean, na.rm = T)


# calculate asymmetry for each barrier
dat$Ecogeo_asym <- abs(dat$Ecogeo1 - dat$Ecogeo2)
dat$ImmigrantInviability_asym <- abs(dat$ImmigrantInviability1 - dat$ImmigrantInviability2)
dat$Pheno_asym <- abs(dat$Pheno1 - dat$Pheno2)
dat$MatingSystem_asym <- abs(dat$MatingSystem1 - dat$MatingSystem2)
dat$FloralIsolation_asym <- abs(dat$FloralIsolation1 - dat$FloralIsolation2)
dat$PollenPistil_asym <- abs(dat$PollenPistil1 - dat$PollenPistil2)
dat$FruitProduction_asym <- abs(dat$FruitProduction1 - dat$FruitProduction2)
dat$SeedProduction_asym <- abs(dat$SeedProduction1 - dat$SeedProduction2)
dat$F1Germination_asym <- abs(dat$F1Germination1 - dat$F1Germination2)
dat$F1Viability_asym <- abs(dat$F1Viability1 - dat$F1Viability2)
dat$F1Sterility_asym <- abs(dat$F1Sterility1 - dat$F1Sterility2)
dat$ExtrinsicPost_asym <- abs(dat$ExtrinsicPost1 - dat$ExtrinsicPost2)


# compile all data into single data.frame
asym_df <- data.frame(
    asymmetry = 
        c(dat$Ecogeo_asym,
        dat$ImmigrantInviability_asym,
        dat$Pheno_asym,
        dat$MatingSystem_asym,
        dat$FloralIsolation_asym,
        dat$PollenPistil_asym,
        dat$FruitProduction_asym,
        dat$SeedProduction_asym,
        dat$F1Germination_asym,
        dat$F1Viability_asym,
        dat$F1Sterility_asym,
        dat$ExtrinsicPost_asym),
    barrier = 
        rep(c("ecogeo", "immigrant_inviability", "pheno", "mating_system", "floral_isolation",
      "pollen_pistil", "fruit_production", "seed_production", "F1_germination",
      "F1_viability", "F1_sterility", "extrinsic_post"), each = nrow(dat)),
    pre_post = c(rep("prezygotic", times = nrow(dat)*6), rep("postzygotic", times = nrow(dat)*6)),
    taxa_type = rep(dat$Taxa_type2, times = 12),
    geography = rep(dat$geography, times = 12)
    )


### 02. calculate summary statistics -------------------------------------------

# convert barrier names to factors (for plotting order)
asym_df$barrier2 = factor(asym_df$barrier,
                          levels = c("ecogeo", "immigrant_inviability", "pheno", "mating_system", "floral_isolation",
                                     "pollen_pistil", "fruit_production", "seed_production", "F1_germination",
                                     "F1_viability", "F1_sterility", "extrinsic_post"))

asym_df$pre_post = factor(asym_df$pre_post,
                          levels = c("prezygotic", "postzygotic"))

# summary statistics
tapply(asym_df$asymmetry, asym_df$barrier2, range, na.rm = T)

tapply(asym_df$asymmetry, asym_df$barrier2, mean, na.rm = T)
sort(tapply(asym_df$asymmetry, asym_df$barrier2, mean, na.rm = T))

tapply(asym_df$asymmetry, asym_df$pre_post, mean, na.rm = T)

# test for differencs in asymmetry between pre- and post-
pre_asym <- asym_df[asym_df$pre_post == "prezygotic", "asymmetry"]
post_asym <- asym_df[asym_df$pre_post == "postzygotic", "asymmetry"]

mean(pre_asym, na.rm = T) # 0.2480052
mean(post_asym, na.rm = T) # 0.2224243

wilcox.test(pre_asym, post_asym)



### 03. Figure 4 ---------------------------------------------------------------
library(ggplot2)

y_labels <- c(
      "Ecogeographic Isolation", 
      "Immigrant Inviability", 
      "Phenology",
      "Mating System", 
      "Floral Isolation",
      "Pollen Pistil Interactions", 
      "Fruit Production", 
      "Seed Production",
      "F1 Germination",
      "F1 Viability", 
      "F1 Sterility",
      "Extrinsic Postzygotic Isolation")


### plot data
dev.new()

ggplot(data = asym_df, aes(x = barrier2, y = asymmetry)) +
      geom_violin(trim = TRUE) +
      geom_jitter(width = 0.075, aes(shape = geography, fill = taxa_type), size = 2, alpha = 0.7) + 
      scale_shape_manual(values = c(24,22,21), na.translate = F) +
      stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
      theme_minimal() +
      theme(
            axis.text.x = element_text(angle = 55, vjust = 0.8, hjust = 0.75),
            axis.text = element_text(size = 11, face = "bold")) +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_text(face = "bold", size = 16)) +
      labs(y = "Barrier asymmetry") +
      scale_x_discrete(labels = y_labels) +
      geom_segment(mapping = aes(x = 6.5, y = 0, xend = 6.5, yend = max(asym_df$asymmetry, na.rm = T)), size = 0.5, color = "black", linetype = "dashed") +
      guides(fill = guide_legend(override.aes=list(shape=23))) +
      labs(fill = "Taxa type", shape = "Geography") +
      theme(legend.title = element_text(size = 14, face = "bold")) +
      theme(legend.text = element_text(size = 12))

# quartz.save(file = "Figure4_FINAL.jpg", type = "jpg", dpi = 600)





