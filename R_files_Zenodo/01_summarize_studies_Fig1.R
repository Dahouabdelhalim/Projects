### 0. import data -------------------------------------------------------------
setwd("~/Documents/Post_Doc/01_Reproductive_Isolation_Review/4_Manuscript/files_for_Dryad_REVISED/")
list.files(pattern = ".csv")
dat <- read.csv(file = "RI_data_FINAL.csv", header = T, stringsAsFactors = F)


### 01. summarize data ---------------------------------------------------------
nrow(dat) # 89 taxa pairs
length(unique(dat$Title)) # 70 studies
range(dat$Year) # 1996-2021

# remove studies that have more than one species pair for descriptive stats
dat2 <- dat[!duplicated(dat$Title),] 

# find most common journals
sort(table(dat2$Journal), decreasing = T)[1:5]
    # Evolution, AJB, Ecology and Evolution, Plant Biology, Annals of Botany

# how many studies were included in Lowry et al. 2008? 
table(dat2$in_Lowry_2008) # 60 NO, 10 YES

# how many taxa pairs were included in Lowry et al. 2008?
table(dat$in_Lowry_2008) # 79 NO, 10 YES

# how many instances of polyploidy? 
table(dat$parents_different_ploidy) # 6 taxa pairs represent cytotypes

# taxonomic summary
length(unique(dat$Family)) # 32 flowering plant families
length(unique(dat$Genus)) # 50 genera

sort(table(dat$Family)) / sum(table(dat$Family))
sort(table(dat$Genus)) / sum(table(dat$Genus))

# taxa type summary
table(dat$Taxa_type) / sum(table(dat$Taxa_type))

# life form 
table(dat$Lifeform) / sum(table(dat$Life_History))

# life history
table(dat$Life_History2) / sum(table(dat$Life_History2))

# geography 
table(dat$geography) / sum(table(dat$geography))

# how many temperate vs. tropical (<23.5 degrees) taxa pairs?
table(dat$temperate_tropical) # 73 temperate pairs, 16 tropical pairs
16/89 # 18% tropical


# how many barriers did studies typically quantify?

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

 # combine "Pollen Sterility" and "Ovule Fertility" columns 
dat$F1Sterility_mean <- apply(dat[, c("F1PollenSterility1", "F1PollenSterility2", "F1OvuleFertility1", "F1OvuleFertility2")], MARGIN = 1, mean, na.rm = T)

dat$ExtrinsicPost_mean <- apply(dat[, c("ExtrinsicPost1", "ExtrinsicPost2")], MARGIN = 1, mean, na.rm = T)

# create new data.frame to count number of barriers each study measured
n_barrier_df <- rbind(
    dat$Ecogeo_mean,
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
    dat$ExtrinsicPost_mean)

barriers_measured <- apply(X = n_barrier_df, MARGIN = 2, FUN = function(x) table(is.na(x))[1])
summary(barriers_measured)
table(barriers_measured)


# how many studies published since 2015 have used Sobel and Chen RI metrics?
dat3 <- subset(dat2, Year > 2014)
sort(table(dat3$SC_metrics), decreasing = T) # 20 YES, 18 NO
sort(table(dat3$SC_metrics), decreasing = T) / sum(table(dat3$SC_metrics))
    # 53% yes, 47% no



### 02. plot all summary statistics (FIGURE 1) ---------------------------------

### A) plot publications over time ---------------------------------------------
dev.new()

par(mfrow = c(3,3))
par(mar = c(4, 5, 4, 4))
par(mgp = c(2.5, 1, 0))

table(dat2$Year)

bp1 <- barplot(table(dat2$Year), xaxt= "n", ylim = c(-3, 12),
    xlab = "Year", ylab = "n publications", font.lab = 2,
    cex.lab = 1.75, xaxt = "n", yaxt = "n")

axis(side = 2, at = c(0,2,4,6,8,10,12), labels = T)

text(x = bp1 - 0.3, y = -1.75, labels = as.numeric(names(table(dat2$Year))), 
     cex = 1, srt = 70)
legend("topright", "A", cex = 2, bty = "n", text.font = 2)
box(which = "figure", lwd = 1, col = "gray")


### B) families ----------------------------------------------------------------
bp2 <- barplot(sort(table(dat$Family), decreasing = T)[1:10], 
               xaxt= "n", yaxt = "n", ylim = c(-22.5,15), xlim = c(-0.25, 11.5),
               xlab = "Family", ylab = "n taxa pairs", 
               font.lab = 2, cex.lab = 1.75)

my_families <- names(sort(table(dat$Family), decreasing = T))[1:10]

text(x = bp2 - 0.6, y = -12, 
     labels = my_families, 
     cex = 1.15, srt = 70)

axis(side = 2, at = c(0,5,10,15), labels = T)

legend("topright", "B", cex = 2, bty = "n", text.font = 2)
box(which = "figure", lwd = 1, col = "gray")


### C) Genus -------------------------------------------------------------------
bp3 <- barplot(sort(table(dat$Genus), decreasing = T)[1:10], 
               xaxt= "n", yaxt = "n", ylim = c(-8,10),
               xlab = "Genus", ylab = "n taxa pairs", 
               font.lab = 2, cex.lab = 1.75)

my_genera <- names(sort(table(dat$Genus), decreasing = T))[1:10]

text(x = bp3 - 0.4, y = -4, 
     labels = my_genera, 
     cex = 1.15, srt = 70)

axis(side = 2, at = c(0,2,4,6,8), labels = T)

legend("topright", "C", cex = 2, bty = "n", text.font = 2)
box(which = "figure", lwd = 1, col = "gray")


### D) taxa type ---------------------------------------------------------------
bp4 <- barplot(sort(table(dat$Taxa_type), decreasing = T), 
               xaxt= "n", yaxt = "n", ylim = c(-57.5, 75),
               xlab = "Taxa type", ylab = "n taxa pairs", 
               font.lab = 2, cex.lab = 1.75)

text(x = bp4 - 0.15, y = -30, 
     labels = names(sort(table(dat$Taxa_type), decreasing = T)), 
     cex = 1.25, srt = 70)

axis(side = 2, at = c(0,25,50,75), labels = T)

legend("topright", "D", cex = 2, bty = "n", text.font = 2)
box(which = "figure", lwd = 1, col = "gray")


### E) Lifeform ----------------------------------------------------------------
bp5 <- barplot(sort(table(dat$Lifeform), decreasing = T), 
               xaxt= "n", yaxt = "n", ylim = c(-25,75),
               xlab = "Lifeform", ylab = "n taxa pairs", 
               font.lab = 2, cex.lab = 1.75)

axis(side = 2, at = c(0,25,50,75), labels = T)

my_lifeforms <- names(sort(table(dat$Lifeform), decreasing = T))

text(x = bp5 - 0.1, y = -15, 
     labels = my_lifeforms, 
     cex = 1.25, srt = 70)

legend("topright", "E", cex = 2, bty = "n", text.font = 2)
box(which = "figure", lwd = 1, col = "gray")


### F) Life History ------------------------------------------------------------
bp6 <- barplot(sort(table(dat$Life_History2), decreasing = T), 
               xaxt= "n", yaxt = "n", ylim = c(-70,75),
               xlab = "Life History", ylab = "n taxa pairs", 
               font.lab = 2, cex.lab = 1.75)

axis(side = 2, at = c(0,25,50,75), labels = T)

my_life_histories <- names(sort(table(dat$Life_History2), decreasing = T))

text(x = bp6 - 0.1, y = -35, 
     labels = my_life_histories, 
     cex = 1.25, srt = 70)

legend("topright", "F", cex = 2, bty = "n", text.font = 2)
box(which = "figure", lwd = 1, col = "gray")


### G) Geography ---------------------------------------------------------------

# geography
sort(table(dat$geography), decreasing = T)

bp7 <- barplot(sort(table(dat$geography), decreasing = T), 
               xaxt= "n", yaxt = "n", ylim = c(-55,75),
               xlab = "Geography", ylab = "n taxa pairs", 
               font.lab = 2, cex.lab = 1.75)

axis(side = 2, at = c(0,25,50,75), labels = T)

my_geogs <- names(sort(table(dat$ge), decreasing = T))

text(x = bp6 - 0.125, y = -30, 
     labels = my_geogs, 
     cex = 1.25, srt = 70)

legend("topright", "G", cex = 2, bty = "n", text.font = 2)
box(which = "figure", lwd = 1, col = "gray")


### H) number of barriers assessed ---------------------------------------------
bp8 <- barplot(table(barriers_measured), 
               xaxt= "n", yaxt = "n", ylim = c(-5,25),
               xlab = "n isolating barriers quantified", ylab = "n studies", 
               font.lab = 2, cex.lab = 1.75)

axis(side = 2, at = c(0,5,10,15,20,25), labels = T)

text(x = bp8, y = -2, 
     labels = 2:8, 
     cex = 1.25, srt = 0)

legend("topright", "H", cex = 2, bty = "n", text.font = 2)
box(which = "figure", lwd = 1, col = "gray")


### I) Sobel and Chen metrics --------------------------------------------------
bp9 <- barplot(sort(table(dat3$SC_metrics), decreasing = T), 
               xaxt= "n", yaxt = "n", ylim = c(-8,30),
               xlab = "Sobel and Chen RI metrics", ylab = "n studies (since 2015)", 
               font.lab = 2, cex.lab = 1.75)

axis(side = 2, at = c(0,5,10,15,20, 25), labels = T)

text(x = bp9, y = -5, 
     labels = c("yes", "no"), 
     cex = 1.25, srt = 70)

legend("topright", "I", cex = 2, bty = "n", text.font = 2)
box(which = "figure", lwd = 1, col = "gray")


# save figure
# quartz.save(file = "Figure1_FINAL.jpg", type = "jpg", dpi = 600)

