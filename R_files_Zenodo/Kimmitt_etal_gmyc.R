library(readxl)
library(dplyr)
library(ape)
library(splits)
library(ggplot2)
library(ggtree)

rm(list=ls())

setwd("")

##### FUNCTIONS #####
# Define gmyc.summary to scriptable
# Modifies output text to distinguish between "cluster" and "entity" CIs for unique grep search
gmyc.summary <- function (object, second.peak = FALSE, ...) 
{
  if (second.peak == TRUE) {
    tmp <- table(cummax(object$likelihood))
    lik.peaks <- names(tmp[tmp > 20])
    peak <- which(object$likelihood == lik.peaks[(length(lik.peaks) - 
                                                    1)])
  }
  cat("Result of GMYC species delimitation\\n")
  cat("\\n\\tmethod:\\t", object[["method"]], sep = "")
  cat("\\n\\tlikelihood of null model:\\t", object$likelihood[1], 
      sep = "")
  if (second.peak == FALSE) {
    cat("\\n\\tmaximum likelihood of GMYC model:\\t", max(object$likelihood), 
        sep = "")
  }
  else {
    cat("\\n\\tmaximum likelihood of GMYC model:\\t", object$likelihood[peak], 
        sep = "")
  }
  if (second.peak == FALSE) {
    LR <- 2 * (max(object$likelihood) - object$likelihood[1])
  }
  else {
    LR <- 2 * (object$likelihood[peak] - object$likelihood[1])
  }
  cat("\\n\\tlikelihood ratio:\\t", LR, sep = "")
  pvalue <- 1 - pchisq(LR, 2)
  cat("\\n\\tresult of LR test:\\t", pvalue, ifelse(pvalue < 0.001, 
                                                 "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, 
                                                                                           "*", "n.s."))), sep = "")
  if (second.peak == FALSE) {
    cat("\\n\\n\\tnumber of ML clusters:\\t", object$cluster[which.max(object$likelihood)], 
        sep = "")
    tmp <- object$cluster[object$likelihood > (max(object$likelihood) - 
                                                 2)]
    cat("\\n\\tcluster confidence interval:\\t", paste(min(tmp), max(tmp), 
                                                    sep = "-"), sep = "")
    cat("\\n\\n\\tnumber of ML entities:\\t", object$entity[which.max(object$likelihood)], 
        sep = "")
    tmp <- object$entity[object$likelihood > (max(object$likelihood) - 
                                                2)]
    cat("\\n\\tentity confidence interval:\\t", paste(min(tmp), max(tmp), 
                                                   sep = "-"), sep = "")
    if (object[["method"]] == "single") {
      cat("\\n\\n\\tthreshold time:\\t", object$threshold.time[which.max(object$likelihood)], 
          "\\n", sep = "")
    }
    else if (object[["method"]] == "multiple" || object[["method"]] == 
             "exhaustive") {
      cat("\\n\\n\\tthreshold time:\\t", object$threshold.time[[which.max(object$likelihood)]], 
          "\\n", sep = " ")
    }
    cat("\\n")
  }
  else {
    cat("\\n\\n\\tnumber of ML clusters:\\t", object$cluster[peak], 
        sep = "")
    cat("\\n\\tnumber of ML entities:\\t", object$entity[peak], 
        sep = "")
    if (object[["method"]] == "single") {
      cat("\\n\\tthreshold time:\\t", object$threshold.time[peak], 
          "\\n", sep = "")
    }
    else if (object[["method"]] == "multiple" || object[["method"]] == 
             "exhaustive") {
      cat("\\n\\tthreshold time:\\t", object$threshold.time[[peak]], 
          "\\n", sep = " ")
    }
    cat("\\n")
  }
}

##### GMYC Runs #####
# Read in data and set up dataframe to store run results
metadata <- #load sample list with populations assigned 
sp_summary <- #load outgroup assignments for each species
newcols <- c("LR_Pvalue","N_Clusters","CI_Cluster_Low", "CI_Cluster_High","Threshold")
sp_summary[newcols] <- NA

sp_list <- sort(unique(sp_summary$species_code))
# Run GMYC on loop - no outgroup
for(i in sp_list){
  # read in mcc.ca.tre
  current_spp <- i 
  spp_data <- dplyr::filter(metadata, species_code==paste(current_spp))
  file_name <- paste0(current_spp,"/",current_spp,"_concatenated_outgroup_combined.mcc.ca.tre")
  mcc_tr <-read.nexus(file_name)
  
  # drop outgroup
  outgroup <- sp_summary$outgroup_individual[sp_summary$species_code==current_spp]
  mcc_tr_pruned <- ape::drop.tip(mcc_tr, outgroup)
  
  # check that tree is fully bifurcating and ultrametric
  is.binary(mcc_tr_pruned)
  is.ultrametric(mcc_tr_pruned)
  
  # run GMYC
  gmyc_result <- gmyc(mcc_tr_pruned) # gmyc algorithm
  gmyc_output <- capture.output(gmyc.summary(gmyc_result)) # store p-value, confidence intervals, threshold number, etc
  clusters <- spec.list(gmyc_result) # store assignment of individual samples to clusters
  pops <- dplyr::select(spp_data, fasta_label, population)
  colnames(clusters) <- c("GMYC_spec","fasta_label")
  clusters_pops <- left_join(pops, clusters)
  
  # write out species assignment df to csv
  csv_pathname <- paste0("gmyc_results/",current_spp,"_gmyc.csv")
  write.csv(clusters_pops, file=csv_pathname)
  
  # put p-value, total number of clusters, cluster confidence interval, and threshold time
  LR_test <- grep("result of LR test", gmyc_output, value=TRUE)
  LR_test <- gsub("\\tresult of LR test:\\t", "", LR_test)
  LR_test <- gsub("([*]+)", "", LR_test)
  LR_test <- gsub("n.s.", "", LR_test)
  sp_summary$LR_Pvalue[sp_summary$species_code==i] <- LR_test
  N_MLclust <- grep("number of ML clusters", gmyc_output, value=TRUE)
  sp_summary$N_Clusters[sp_summary$species_code==i] <- gsub("\\tnumber of ML clusters:\\t", "", N_MLclust)
  if(sp_summary$LR_Pvalue[sp_summary$species_code==i] > 0.05) {
    sp_summary$N_Clusters[sp_summary$species_code==i] = 1}
  Cluster_CI <- grep("cluster confidence interval", gmyc_output, value=TRUE)
  Cluster_CI <- gsub("\\tcluster confidence interval:\\t", "", Cluster_CI)
  sp_summary$CI_Cluster_Low[sp_summary$species_code==i] <- gsub("(-.+)", "", Cluster_CI) # lower bound of CI
  sp_summary$CI_Cluster_High[sp_summary$species_code==i] <- gsub("(.+-)", "", Cluster_CI) # upper bound of CI
  Thresh_time <- grep("threshold time", gmyc_output, value=TRUE)
  sp_summary$Threshold[sp_summary$species_code==i] <- gsub("\\tthreshold time:\\t", "", Thresh_time)
}

write.csv(sp_summary, "gmyc_results/gmyc_all_spp.csv")

# Visualize trees and cluster assignments
pdf("gmyc_trees_allspp.pdf")
for(i in sp_list){
  # read in mcc.ca.tre
  current_spp <- i 
  csv_pathname <- paste0("gmyc_results/",current_spp,"_gmyc.csv")
  clusters_pops <- read.csv(csv_pathname)
  tree_name <- paste0(current_spp,"/",current_spp,"_concatenated_outgroup_combined.mcc.ca.tre")
  mcc_tr <-read.nexus(tree_name)
  
  # cluster, p-value, threshold variable
  nclusters <- sp_summary$N_Clusters[sp_summary$species_code==i]
  pvalue <- sp_summary$LR_Pvalue[sp_summary$species_code==i]
  threshold <- as.numeric(sp_summary$Threshold[sp_summary$species_code==i])
  
  # prepare clusters_pops to add geographic data to tree
  clusters_pops$X <- NULL
  clusters_pops$GMYC_spec <- as.character(clusters_pops$GMYC_spec)
  
  plasma_colors <- c("#fcce25","#ea7457","#b12a90","#5601a4")
  p <- ggtree(mcc_tr) %<+% clusters_pops # attach data about tips to the tree
  p2 <- p + theme_tree2() + # theme_tree2() adds timescale
    geom_tiplab(as_ylab=TRUE, size=4) +
    geom_tippoint(aes(color=population), size=2, position=position_nudge(x=3.0)) + # adds tip pops
    scale_color_manual(values = plasma_colors, name="geographic\\npopulation") +
    ggtitle(paste0(current_spp,"; Clusters = ",nclusters,"; p = ",pvalue)) 
  p3 <- revts(p2) + scale_x_continuous(labels=abs) +
    geom_vline(xintercept=(threshold), color="darkgray", linetype="solid", lwd=0.5) # make timescale project into past Ma
  print(p3)
}
dev.off()

# for subset of species with significant p-values, annotate trees with geographic pop and genetic cluster
sig_sp <- sort(sp_summary$species_code[sp_summary$N_Clusters>1])
pdf("gmyc_trees_clusterpops.pdf")
for(i in sig_sp){
  # read in mcc.ca.tre
  current_spp <- i 
  csv_pathname <- paste0("gmyc_results/",current_spp,"_gmyc.csv")
  clusters_pops <- read.csv(csv_pathname)
  tree_name <- paste0(current_spp,"/",current_spp,"_concatenated_outgroup_combined.mcc.ca.tre")
  mcc_tr <-read.nexus(tree_name)
  
  # cluster, p-value, threshold variable
  nclusters <- sp_summary$N_Clusters[sp_summary$species_code==i]
  pvalue <- sp_summary$LR_Pvalue[sp_summary$species_code==i]
  threshold <- as.numeric(sp_summary$Threshold[sp_summary$species_code==i])
  
  # prepare two dfs with geographic population and genetic cluster to annotate tips
  clusters_pops$X <- NULL
  geog_pops <- clusters_pops[,c(1,2)]
  gmyc_clust <- clusters_pops[,c(1,3)]
  rownames(geog_pops) <- geog_pops$fasta_label
  geog_pops$fasta_label <- NULL
  rownames(gmyc_clust) <- gmyc_clust$fasta_label
  gmyc_clust$fasta_label <- NULL
  gmyc_clust$GMYC_spec <- as.character(gmyc_clust$GMYC_spec)
  
  # set up custom color palette
  plasma_colors <- c("#fcce25","#ea7457","#b12a90","#5601a4")
  # get outgroup label
  outgroup <- mcc_tr[["tip.label"]][!(mcc_tr[["tip.label"]] %in% clusters_pops$fasta_label)]
  mcc_tr[["tip.label"]][mcc_tr[["tip.label"]]==outgroup] <- "Outgroup"
  # make trees
  p <- ggtree(mcc_tr) %<+% clusters_pops
  p2 <- p +
    geom_tippoint(aes(color=population), size=2, position=position_nudge(x=3.0)) +
    geom_tiplab(aes(subset=(mcc_tr[["tip.label"]] == "Outgroup")), size=3) +
    scale_color_manual(values = plasma_colors, name="geographic\\npopulation") +
    ggtitle(paste(sp_summary$species[sp_summary$species_code==i]))
  p3 <- gheatmap(p2, gmyc_clust, offset=0, width=0.05, colnames=FALSE) +
    scale_fill_viridis_d(option="G", name="genetic\\ncluster")
  print(p3)
}
dev.off()

