library(geomorph)

load("data/procrustes_INT.Rdata")
load("data/procrustes_EXT.Rdata")
md <- read.csv("data/int_mean.csv")
tree <- read.tree("data/McGee2020_tree.tre")

# intersection of internal and external
name_overlap <- intersect(names(int_mean), names(ext_mean))
int_mean <- int_mean[match(name_overlap, names(int_mean))]
ext_mean <- ext_mean[match(name_overlap, names(ext_mean))]
tree <- keep.tip(tree, name_overlap)

# convert to an array
int_mean <- array(as.numeric(unlist(int_mean)),
                  dim = c(dim(int_mean[[1]]),
                          length(int_mean)))
dimnames(int_mean)[[3]] <- name_overlap
ext_mean <- array(as.numeric(unlist(ext_mean)),
                  dim = c(dim(ext_mean[[1]]),
                          length(ext_mean)))
dimnames(ext_mean)[[3]] <- name_overlap

# match names
md <- md[match(name_overlap, md$mcgee_label), ]
inxn <- paste(md$reproduction, md$feeding, sep = " & ")
mb_idx <- which(md$reproduction == "mouthbrooding" &
                  md$feeding == "other")

# run partial-block least squares
pbls <- phylo.integration(int_mean, ext_mean, tree, iter = 999)

png("../draft/images/S01_PLS.png", width = 6, height = 6,
    units = "in",
    res = 150)

# plot results; one point = one species
# external and internal head shape are certainly correlated
plot(pbls, type = "n", 
     cex.axis = 1.3,
     cex.lab = 1.3)

# highlight the mouthbrooding, non-winnowing species
cols <- c("#002CA5", grey(0.7), "goldenrod1", "#39D69A")
for (i in 1:length(unique(inxn))) {
  idx <- which(inxn == unique(inxn)[i])
  points(pbls$XScores[idx, 1],
         pbls$YScores[idx, 1],
         pch = 19, cex = 1.2,
         col = cols[i])
}

legend(-.05, .02, c("Neither behavior",
                    "Mouthbrooding only",
                    "Winnowing only",
                    "Both behaviors"),
       fill = cols[c(2, 3, 4, 1)])
dev.off()


# all have positive residuals - larger buccal cavities than predicted
# given their head shape
