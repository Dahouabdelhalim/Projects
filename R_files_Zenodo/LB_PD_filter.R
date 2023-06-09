#!/usr/bin/env Rscript
require(ape)

options(digits=10)

nzmean <- function(x){
    if (all(x==0)) 0 else mean(x[x!=0])
}

cat("Gene", "\\t", "LB_score_Heterogeneity", "\\t", "Average_PD", "\\t", "PD_range", "\\n", file="LB_scores_per_tree.txt", append=FALSE)

file_list <- Sys.glob("*.tre")
    for (file_name in file_list){
    tree <- read.tree(file_name)
    tip_to_tip <- cophenetic(tree)
    PD_range <- range(tip_to_tip[tip_to_tip!=0])[2] - range(tip_to_tip[tip_to_tip!=0])[1]
    tip_means <- apply(tip_to_tip, 1, nzmean)
    LB <- c()
    for (x in tip_means){
        score <- (x/nzmean(tip_to_tip)-1) * 100
        LB <- c(LB, score)
    }
    cat(file_name, "\\t", sd(LB), "\\t", nzmean(tip_to_tip), "\\t", PD_range, "\\n", file="LB_scores_per_tree.txt", append=TRUE)
}

scores <- read.table("LB_scores_per_tree.txt", sep="\\t", header=TRUE)
LBH = scores$LB_score_Heterogeneity
APD = scores$Average_PD

out_LBH = median(LBH) + (1.5 * (IQR(LBH)))
out_APD = median(APD) + (1.5 * (IQR(APD)))

rejected <- subset(scores, LBH >= out_LBH | APD >= out_APD, select=c(Gene))
write.table(rejected, file="LBH_APD_outliers.txt", append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE)