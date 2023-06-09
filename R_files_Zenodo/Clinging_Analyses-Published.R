# Clinging analyses and Figures
library(geomorph)
library(ape)

footshape <- readland.tps("footshape.tps", specID = "ID")
full_dataset <- read.csv("full_dataset.csv")
phylo <- read.tree("phylo.tre")
phylo$tip.label <- stringr::str_replace_all(phylo$tip.label, "_", " ") 

GDF <- geomorph.data.frame(foot = footshape, cling_0 = full_dataset$Max_angle_0µm, 
                           cling_2000 = full_dataset$Max_angle_2000µm,
                           microhab_strict = as.factor(full_dataset$Strict), 
                           centroid = full_dataset$centroid_pruned....1., 
                           cent_corrected = full_dataset$centroid_pruned....1./full_dataset$fsa_mass_mean,
                           fsa_corrected = full_dataset$fsa_mean/full_dataset$fsa_mass_mean,
                           mass = full_dataset$fsa_mass_mean,
                           phy = phylo)

## Cling ~ microhab

anov <- procD.pgls(cling_0~microhab_strict, data = GDF, phy = phy)
summary(anov) # z = 0.15041, p = 0.454

anov <- procD.pgls(cling_2000~microhab_strict, data = GDF, phy = phy)
summary(anov) # z = 0.1561, p = 0.428

## Landmark-based Shape analyses

### smooth
cling <- full_dataset$Max_angle_0µm
cling <- as.numeric(cling)
names(cling) <- full_dataset$Species
identical(GDF$phy$tip.label, dimnames(GDF$foot)[[3]])

results_full <- phylo.integration(cling, GDF$foot, GDF$phy) # same math as two.b.pls but can add phylo
results_full

cling_a <- cling[which(GDF$microhab_strict == "A")]
foot_a <- GDF$foot[,,which(GDF$microhab_strict == "A")]
phy_a <- keep.tip(GDF$phy, GDF$phy$tip.label[which(GDF$microhab_strict == "A")])

results_a <- phylo.integration(cling_a, foot_a, phy_a) # same math as two.b.pls but can add phylo
results_a

cling_t <- cling[which(GDF$microhab_strict == "T")]
foot_t <- GDF$foot[,,which(GDF$microhab_strict == "T")]
phy_t <- keep.tip(GDF$phy, GDF$phy$tip.label[which(GDF$microhab_strict == "T")])

results_t <- phylo.integration(cling_t, foot_t, phy_t) # same math as two.b.pls but can add phylo
results_t

cling_w <- cling[which(GDF$microhab_strict == "W")]
foot_w <- GDF$foot[,,which(GDF$microhab_strict == "W")]
phy_w <- keep.tip(GDF$phy, GDF$phy$tip.label[which(GDF$microhab_strict == "W")])

results_w <- phylo.integration(cling_w, foot_w, phy_w) # same math as two.b.pls but can add phylo
results_w

smooth_results <- rbind(c(results_full$r.pls, results_full$Z, results_full$P.value),
                         c(results_a$r.pls, results_a$Z, results_a$P.value),
                         c(results_t$r.pls, results_t$Z, results_t$P.value),
                         c(results_w$r.pls, results_w$Z, results_w$P.value))
### rough
cling <- full_dataset$Max_angle_2000µm
cling <- as.numeric(cling)
names(cling) <- full_dataset$Species

results <- phylo.integration(cling, GDF$foot, GDF$phy) # same math as two.b.pls but can add phylo
results_full

cling_a <- cling[which(GDF$microhab_strict == "A")]
foot_a <- GDF$foot[,,which(GDF$microhab_strict == "A")]
phy_a <- keep.tip(GDF$phy, GDF$phy$tip.label[which(GDF$microhab_strict == "A")])

results_a <- phylo.integration(cling_a, foot_a, phy_a) # same math as two.b.pls but can add phylo
results_a

cling_t <- cling[which(GDF$microhab_strict == "T")]
foot_t <- GDF$foot[,,which(GDF$microhab_strict == "T")]
phy_t <- keep.tip(GDF$phy, GDF$phy$tip.label[which(GDF$microhab_strict == "T")])

results_t <- phylo.integration(cling_t, foot_t, phy_t) # same math as two.b.pls but can add phylo
results_t

cling_w <- cling[which(GDF$microhab_strict == "W")]
foot_w <- GDF$foot[,,which(GDF$microhab_strict == "W")]
phy_w <- keep.tip(GDF$phy, GDF$phy$tip.label[which(GDF$microhab_strict == "W")])

results_w <- phylo.integration(cling_w, foot_w, phy_w) # same math as two.b.pls but can add phylo
results_w

rough_results <- rbind(c(results_full$r.pls, results_full$Z, results_full$P.value),
                        c(results_a$r.pls, results_a$Z, results_a$P.value),
                        c(results_t$r.pls, results_t$Z, results_t$P.value),
                        c(results_w$r.pls, results_w$Z, results_w$P.value))


### Table 2 
all_results <- rbind(smooth_results, rough_results)
all_results <- cbind(rep(c("Full", "A", "T", "W"), 2), all_results)
all_results <- cbind(rep(c("Smooth", "Rough"), each = 4), all_results)
colnames(all_results) <- c("Substrate", "Data Subset", "r-PLS", "Effect Size (Z)", "P-value")
all_results

## Centroid Size
            
### smooth
results_full <- procD.pgls(cling_0 ~ cent_corrected*microhab_strict, data = GDF, phy = phy) # same math as two.b.pls but can add phylo
summary(results_full)
results_full$pgls.coefficients
results_full$aov.table # Table S1 Csize Smooth

results_full <- procD.pgls(cling_0 ~ cent_corrected*microhab_strict + mass, data = GDF, phy = phy) # same math as two.b.pls but can add phylo
summary(results_full)
results_full$pgls.coefficients
results_full$aov.table # Table 3 Csize Smooth

microhab_col <- as.character(GDF$microhab_strict)
microhab_col[which(microhab_col=="T")] <- "black"
microhab_col[which(microhab_col=="A")] <- "chartreuse4"
microhab_col[which(microhab_col=="W")] <- "cornflowerblue"

plot(GDF$cling_0 ~ GDF$cent_corrected, col = microhab_col, pch = 19) # Figure 1B

### rough
results_full <- procD.pgls(cling_2000 ~ cent_corrected*microhab_strict, data = GDF, phy = phy) # same math as two.b.pls but can add phylo
summary(results_full)
results_full$pgls.coefficients
results_full$aov.table # Table S1 Csize Rough

results_full <- procD.pgls(cling_2000 ~ cent_corrected*microhab_strict + mass, data = GDF, phy = phy) # same math as two.b.pls but can add phylo
summary(results_full)
results_full$pgls.coefficients
results_full$aov.table # Table 3 Csize Rough

microhab_col <- as.character(GDF$microhab_strict)
microhab_col[which(microhab_col=="T")] <- "black"
microhab_col[which(microhab_col=="A")] <- "chartreuse4"
microhab_col[which(microhab_col=="W")] <- "cornflowerblue"

plot(GDF$cling_2000 ~ GDF$cent_corrected, col = microhab_col, pch = 19) # not published - version of figure 1B with ROUGH data
            
## FSA
            
### smooth
results_full <- procD.pgls(cling_0 ~ fsa_corrected*microhab_strict, data = GDF, phy = phy) # same math as two.b.pls but can add phylo
summary(results_full)
results_full$pgls.coefficients
results_full$aov.table # Table S1 FSA Smooth

results_full <- procD.pgls(cling_0 ~ fsa_corrected*microhab_strict + mass, data = GDF, phy = phy) # same math as two.b.pls but can add phylo
summary(results_full)
results_full$pgls.coefficients
results_full$aov.table # Table 3 FSA Smooth

microhab_col <- as.character(GDF$microhab_strict)
microhab_col[which(microhab_col=="T")] <- "black"
microhab_col[which(microhab_col=="A")] <- "chartreuse4"
microhab_col[which(microhab_col=="W")] <- "cornflowerblue"

plot(GDF$cling_0 ~ GDF$fsa_corrected, col = microhab_col, pch = 19) # Fig 1C

### rough
results_full <- procD.pgls(cling_2000 ~ fsa_corrected*microhab_strict, data = GDF, phy = phy) # same math as two.b.pls but can add phylo
summary(results_full)
results_full$pgls.coefficients
results_full$aov.table # Table S1 FSA Rough

results_full <- procD.pgls(cling_2000 ~ fsa_corrected*microhab_strict + mass, data = GDF, phy = phy) # same math as two.b.pls but can add phylo
summary(results_full)
results_full$pgls.coefficients
results_full$aov.table # Table 3 FSA Rough

microhab_col <- as.character(GDF$microhab_strict)
microhab_col[which(microhab_col=="T")] <- "black"
microhab_col[which(microhab_col=="A")] <- "chartreuse4"
microhab_col[which(microhab_col=="W")] <- "cornflowerblue"

plot(GDF$cling_2000 ~ GDF$fsa_corrected, col = microhab_col, pch = 19) # not published - version of figure 1C with ROUGH data


## Fig 1
            
plot_this <- gm.prcomp(phy = GDF$phy, A = GDF$foot)
lks <- cbind((1:20), (2:21))

## assigning colors to tips
ArbCol <- as.character(GDF$microhab_strict)
ArbCol[which(ArbCol == "A")] <- "chartreuse4"
ArbCol[which(ArbCol == "W")] <- "cornflowerblue"
ArbCol[which(ArbCol == "T")] <- "black"
            
png("Fig1A.png", res = 100, height = 790, width = 1000)
            
par(mfrow = c(1,1), mar = c(5,4,1,1))
plot(plot_this, phylo = T, cex = GDF$cling_0^2/10000, pch = 19, col = ArbCol,
     phylo.par = list(tip.labels = F, node.labels = F, node.cex = 0.8)) # pca of foot shape to identify species names
legend(.085,-.08, legend = c("90°","135°","180°"), pch = 19, col = "gray",
       pt.cex = c(0.81,1.8225,3.24),x.intersp = c(1,1,1.3), border = "black",
       y.intersp = 1.2, horiz = F, title = "Max Cling Angle")
legend(.14,-.08, legend = c("A","T","W"), fill = c("chartreuse4", "black", "cornflowerblue"),
       pt.cex = c(0.81,1.8225,3.24),x.intersp = c(1,1,1.3), border = "black",
       y.intersp = 1.2, horiz = F, title = "Microhabitat")
par(fig = c(0,0.35, 0.07, 0.42), new = T)  
plotRefToTarget(M1 = GDF$foot[,,1], M2 = mshape(GDF$foot), links = lks, add = T, mag = 1.5) # a maculatum
par(fig = c(0.65,1, 0.45, 0.8), new = T)  
plotRefToTarget(M1 = GDF$foot[,,11], M2 = mshape(GDF$foot), links = lks, add = T, mag = 1.5) # b franklini
par(fig = c(0.03,0.38, 0.6, 0.95), new = T, mar = c(0,0,6,5))  
plotRefToTarget(M1 = GDF$foot[,,6], M2 = mshape(GDF$foot), links = lks, add = T, mag = 1.5) # d ocoee

dev.off()


## assigning colors 
microhab_col <- as.character(GDF$microhab_strict)
microhab_col[which(microhab_col=="T")] <- "black"
microhab_col[which(microhab_col=="A")] <- "chartreuse4"
microhab_col[which(microhab_col=="W")] <- "cornflowerblue"

png("Fig1BC.png", res = 100, height = 790, width = 500)
par(mfrow = c(2,1), mar = c(4,4,1,1))
plot(GDF$cling_0 ~ GDF$cent_corrected, col = microhab_col,cex = 1.5,  pch = 19, ylab = "Max Cling Angle", xlab = "Centroid Size/Mass")
plot(GDF$cling_0 ~ GDF$fsa_corrected, col = microhab_col, cex = 1.5, pch = 19, ylab = "Max Cling Angle", xlab = "FSA/Mass")
dev.off()


## Figure S1
png("FigS1.png", res = 100, height = 790, width = 600)
par(mfrow = c(1,1), mar = c(2,2,2,2))
plot(GDF$phy)
dev.off()
