########################################################################################
## FUNCTIONS FOR DATA ANALYSIS AND NETWORK CONSTRUCTION
##
## These script contains the functions used to reproduce results generated in paper
## "A systems level analysis of epileptogenesis-associated proteome alterations" by
## M. Keck and G. Androsova et al. (2016).
## 
## Script author: Ganna Androsova, ganna.androsova@uni.lu
## 
########################################################################################

# ipak function was developed by Steven Worthington and deposited at 
# https://gist.github.com/stevenworthington/3178163
ipak = function(pkg){
  new.pkg = pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(c("openxlsx", "FactoMineR", "qvalue", "WGCNA", "flashClust", "made4", "limma"))

if (!require("factoextra")) {
  install_github("kassambara/factoextra")
  library("factoextra")
}

merge_original_tables = function(original_directory, excel_files, brain_region){
  
  original_data = lapply(excel_files, function(file){
    input_file = read.xlsx(file.path(original_directory, file), 
                           sheet = 1, colNames = TRUE, startRow = 3)
    
    #Transform excel output into data frame and keep only expression values
    df = data.frame(input_file[which(input_file[,7]!="N/A"),18:27])
    rownames(df) = input_file[which(input_file[,7]!="N/A"),7]
    colnames(df) = c(rep("Ctrl",5), rep("SE", 5))
    
    return(df)
  })
  
  merged_full_data = original_data[[1]]
  for(i in c(2:length(original_data))){
    merge2 = merge(merged_full_data, original_data[[i]], by="row.names", all = TRUE)
    merged_full_data = merge2[,2:ncol(merge2)]
    rownames(merged_full_data) = merge2[,1]
  }
  colnames(merged_full_data) <- c(rep("Ctrl_2d",5), rep("SE_2d", 5), 
                                  rep("Ctrl_10d",5), rep("SE_10d", 5), 
                                  rep("Ctrl_8w",5), rep("SE_8w", 5))
  
  png(file = paste(output_directory, "/", brain_region, "_expression_before_norm.png", sep = ""), 
      height = 5, width = 7, units="in", res=300)
  boxplot(merged_full_data, xlab="Samples", ylab="Protein abundances", 
          main="Protein expression values before normalization")
  dev.off()
  
  return(list(merged_full_data, original_data))
}


plot_PCA = function(filtered_measurements, output_directory, brain_region, type){
  #The proteins with highest contribution to the variability of data tend to be on horizontal axis (of principle component 1) thus it is difficult to have a clear separation
  res.pca = PCA(filtered_measurements, quali.sup=ncol(filtered_measurements), graph=FALSE)
  
  png(file = paste(output_directory, "/", brain_region, "_", type, "_PCA.png", sep = ""), height = 5, width = 5, units="in", res=300)
  print(fviz_pca_ind(res.pca, pointsize=4, invisible="quali", label="none", habillage=ncol(filtered_measurements)) + labs(title = brain_region) + scale_color_brewer(palette="Set1") + theme_minimal())
  dev.off()
  print(fviz_pca_ind(res.pca, pointsize=4, invisible="quali", label="none", habillage=ncol(filtered_measurements)) + labs(title = brain_region) + scale_color_brewer(palette="Set1") + theme_minimal())
  
  top_10_contributors = sapply(1:ncol(res.pca$var$contrib), function(x){names(sort(res.pca$var$contrib[,x], decreasing = T)[1:10])})
  colnames(top_10_contributors) = colnames(res.pca$var$contrib)
  write.table(top_10_contributors, file=paste(output_directory, "/", brain_region, "_", type, "_top10_PCA_proteins.txt", sep = ""), sep="\\t", row.names=F, quote=FALSE)
}


plot_separate_PCA = function(separate_data, output_directory, brain_region, time_label){
  sapply(c(1:length(separate_data)), function(i){
    time_point = t(asinh(separate_data[[i]]))
    rownames(time_point) = NULL
    time_point = as.data.frame(time_point)
    filtered_measurements = time_point[, colSums(is.na(time_point)) != nrow(time_point)]
    filtered_measurements$group = c(rep("Ctrl", 5), rep("SE", 5))
    plot_PCA(filtered_measurements, output_directory, brain_region, time_label[i])
  })
}


normalization = function(merged_data, input_directory, output_directory, brain_region, data.type){
  
  #Arcsin transformation and median centering
  asin_df = as.data.frame(asinh(merged_data))
  preprocessed_data = t(apply(asin_df,2,function(x){x-median(x[!is.na(x)])}))
  
  png(file = paste(output_directory, "/", brain_region, "_expression_after_norm.png", sep = ""), height = 6, width = 8, units="in", res=300)
  boxplot(t(preprocessed_data), xlab="Samples", ylab="Protein abundances", main="Protein expression values after normalization")
  dev.off()
  outliers = NULL
  for (i in 1:nrow(preprocessed_data)){
    outliers = c(outliers, names(which(preprocessed_data[i,]<(-10))))
  }
  outliers = unique(outliers)
  
  write.csv(preprocessed_data, file = file.path(input_directory, paste0(brain_region, "_preprocessed_data.csv")), quote=F)
  write.table(outliers, file = file.path(output_directory, paste0(brain_region, "_outliers.txt", sep="")), quote=F, row.names=F, col.names=F)
  return(preprocessed_data)
}


DE_prots = function(data, output_directory, time_label, brain_region){
  data = t(data)
  design = model.matrix(~factor(c(rep(1,5), rep(2, 5))))
  stat_results = lapply(seq(0, 20, by=10), function(i){
    subset = data[,(i+1):(i+10)]
    subset = subset[!!rowSums(!is.na(subset)),]
    fit = lmFit(subset, design)
    fit_ebayes = eBayes(fit)
    logFC = fit_ebayes$coefficients[, 2]
    FC = exp(logFC)
    t_tset = fit_ebayes$t[, 2]
    p_values = fit_ebayes$p.value[, 2]
    q_values = qvalue(p_values)$q
    results = data.frame(FC, logFC, t_tset, p_values, q_values)
  })
  names(stat_results) = time_label
  
  DE_table = stat_results[[1]][,c("p_values", "FC")]
  for(i in c(2:length(stat_results))){
    m1=merge(DE_table, stat_results[[i]][,c("p_values", "FC")], by="row.names", all = T)
    DE_table=m1[,2:ncol(m1)]
    rownames(DE_table) = m1[,1]
  }
  colnames(DE_table) <- c(paste(rep(time_label, each=2), c("p-values", "Fold Change")))
  DE_table = round(DE_table,3)
  write.csv2(DE_table, file = file.path(output_directory, paste0(brain_region, "_DEG.csv")))
  
  return(stat_results)
}


hist_FC = function(stat_results, output_directory, brain_region){
  png(file = paste(output_directory, "/", brain_region, "_FC_distribution.png", sep = ""), height = 3, width = 6, units="in", res=300)
  par(mfrow=c(1,3))
  for(i in 1:length(stat_results)){
    hist(stat_results[[i]]$logFC, main = paste(brain_region, names(stat_results)[i]), breaks=20, xlab="log(fold change)")
  }
  dev.off()
  
  #Plotting into Markdown
  par(mfrow=c(1,3))
  for(i in 1:length(stat_results)){
    hist(stat_results[[i]]$logFC, main = paste(brain_region, names(stat_results)[i]), breaks=20, xlab="log(fold change)")
  }
}


volcano_DE = function(stat_results, output_directory, brain_region){
  png(file = paste(output_directory, "/", brain_region, "_DE_volcano.png", sep = ""), height = 3, width = 7, units="in", res=300)
  par(mfrow=c(1,3))
  for(i in 1:length(stat_results)){
    ry <- c(0, ceiling(max(-log10(stat_results[[i]]$p_values))))
    plot(stat_results[[i]]$logFC, -log10(stat_results[[i]]$p_values), pch=21, bg="lightgrey", cex=0.9, 
         xlab="log(fold change)", ylab="-log10(p-value)", xlim=c(-2,2), main = paste(brain_region, names(stat_results)[i]))
    abline(h=(-log10(0.05)), col="grey", lty="dotted")
    
    # Add colored points: red if p-value < 0.05, orange of log2FC > 1, green if both)
    sub1 = rownames(subset(stat_results[[i]], p_values<.05))
    points(stat_results[[i]][sub1, "logFC"], -log10(stat_results[[i]][sub1, "p_values"]), pch=20, col="red")
    sub2 = rownames(subset(stat_results[[i]], stat_results[[i]]$FC<=0.67|stat_results[[i]]$FC>=1.5))
    points(stat_results[[i]][sub2, "logFC"], -log10(stat_results[[i]][sub2, "p_values"]), pch=20, col="orange")
    sub3 = rownames(subset(stat_results[[i]], p_values<.05 & (stat_results[[i]]$FC<=0.67|stat_results[[i]]$FC>1.5)))
    points(stat_results[[i]][sub3, "logFC"], -log10(stat_results[[i]][sub3, "p_values"]), pch=20, col="green")
  }
  dev.off()
  
  #Plotting into Markdown
  par(mfrow=c(1,3))
  for(i in 1:length(stat_results)){
    ry <- c(0, ceiling(max(-log10(stat_results[[i]]$p_values))))
    plot(stat_results[[i]]$logFC, -log10(stat_results[[i]]$p_values), pch=21, bg="lightgrey", cex=0.9, 
         xlab="log(fold change)", ylab="-log10(p-value)", xlim=c(-2,2), main = paste(brain_region, names(stat_results)[i]))
    abline(h=(-log10(0.05)), col="grey", lty="dotted")
    
    # Add colored points: red if p-value < 0.05, orange of log2FC > 1, green if both)
    sub1 = rownames(subset(stat_results[[i]], p_values<.05))
    points(stat_results[[i]][sub1, "logFC"], -log10(stat_results[[i]][sub1, "p_values"]), pch=20, col="red")
    sub2 = rownames(subset(stat_results[[i]], stat_results[[i]]$FC<=0.67|stat_results[[i]]$FC>=1.5))
    points(stat_results[[i]][sub2, "logFC"], -log10(stat_results[[i]][sub2, "p_values"]), pch=20, col="orange")
    sub3 = rownames(subset(stat_results[[i]], p_values<.05 & (stat_results[[i]]$FC<=0.67|stat_results[[i]]$FC>1.5)))
    points(stat_results[[i]][sub3, "logFC"], -log10(stat_results[[i]][sub3, "p_values"]), pch=20, col="green")
  }
}


barplot_DE = function(stat_results, output_directory, brain_region){
  list_of_DE_prots = lapply(stat_results, function(x){
    prot=subset(x, p_values<.05 & (x$FC<=0.67|x$FC>=1.5))
    write.table(prot,
                file=file.path(output_directory, 
                               paste0(brain_region, " DE at ", names(stat_results[parent.frame()$i[]]), ".txt")), 
                quote = F, sep = "\\t")
    nrow(prot)
  })
  
  up_down_regulated = lapply(stat_results, function(x){
    sub1 = nrow(subset(x, p_values<.05 & x$FC>=1.5))
    sub2 = nrow(subset(x, p_values<.05 & x$FC<=0.67))
    up_down_regulated = paste(paste("(+)", sub1), paste("(-)",sub2), sep="\\n")
  })
  
  DEs =  unlist(list_of_DE_prots)
  colours = c("black", "dimgrey", "gray")
  
  png(file = paste(output_directory, "/", brain_region, "_DE_proteins.png", sep = ""), height = 5, width = 4, units="in", res=300)
  par(mar=c(4,5,2,2))
  bp1 = barplot(DEs[1:3], col=colours, names=names(DEs[1:3]), ylim=c(0,max(DEs[1:3])+100), ylab ="Differentially expressed proteins", main=brain_region)
  text(bp1, DEs[1:3], up_down_regulated[1:3], pos=3)
  dev.off()
  
  par(mar=c(4,5,2,2))
  bp1 = barplot(DEs[1:3], col=colours, names=names(DEs[1:3]), ylim=c(0,max(DEs[1:3])+100), ylab ="Differentially expressed proteins", main=brain_region)
  text(bp1, DEs[1:3], up_down_regulated[1:3], pos=3)
}


get_common_prots = function(original_data, input_directory, 
                            output_directory, brain_region){
  
  merged_common_data = original_data[[1]]
  for(i in c(2:length(original_data))){
    merge1 = merge(merged_common_data, original_data[[i]], by="row.names", all = FALSE)
    merged_common_data = merge1[,2:ncol(merge1)]
    rownames(merged_common_data) = merge1[,1]
  }
  colnames(merged_common_data) <- c(rep("Ctrl_2d",5), rep("SE_2d", 5), 
                                    rep("Ctrl_10d",5), rep("SE_10d", 5), 
                                    rep("Ctrl_8w",5), rep("SE_8w", 5))
  normalized_data = normalization(merged_common_data, input_directory, 
                output_directory, brain_region)
  return(normalized_data)
}

visualize_correlation = function(correlation_matrix, output_directory, brain_region){
  hmcols = colorRampPalette(c("white", "red"))(256)
  png(file = paste(output_directory, "/", brain_region, "_correlation_heatmap.png", sep = ""), height = 9, width = 9, units="in", res=300)
  heatmap(correlation_matrix, col=hmcols, main = brain_region)
  dev.off()
  heatmap(correlation_matrix, col=hmcols, main = brain_region)
}


estimate_correlation = function(preprocessed_data){
  # Indicate the method for estimation of the protein expression similarity
  options = "use = 'p', method = 'spearman'"
  
  # Create adjacency matrix required for network construction
  adjacency_matrix = adjacency(preprocessed_data, corFnc = "cor", corOptions = options, type = "unsigned", power = 6)
  adjacency_matrix[is.na(adjacency_matrix)] = 0
  return(adjacency_matrix)
}


estimate_network_topology = function(adjacency_matrix, output_directory, brain_region){
  topological_table = matrix(rep(0,5), nrow = 1, ncol = 5)
  colnames(topological_table) = c("Density", "Centralization", "Heterogeneity", 
                                  "Mean clustering coefficient", "Mean scaled connectivity")
  
  # Estimation of the network topological parameters
  statistics = fundamentalNetworkConcepts(adjacency_matrix, GS = NULL)
  topological_table[,1] = statistics$Density
  topological_table[,2] = statistics$Centralization
  topological_table[,3] = statistics$Heterogeneity
  topological_table[,4] = mean(statistics$ClusterCoef)
  topological_table[,5] = mean(statistics$ScaledConnectivity)
  
  # Write-down the table
  write.table(topological_table, file=file.path(output_directory, paste(brain_region, "_network_topological_parameters.txt", sep = "")), sep="\\t", row.names=F, quote=FALSE)
  return(topological_table)
}


detect_modules = function(adjacency_matrix, output_directory, brain_region){
  # Calculation of the topological overlap matrix
  TOM = TOMsimilarity(adjacency_matrix)
  # Calculation of topological overlap dissimilarity
  dissTOM = 1-TOM
  
  # Call the hierarchical clustering function
  geneTree = flashClust(as.dist(dissTOM), method = "average");
  
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, 
                              pamRespectsDendro = FALSE, minClusterSize = 20);
  print(table(dynamicMods))
  
  # Convert numeric lables into colors
  module_asignment = labels2colors(dynamicMods)
  print(table(module_asignment))
  
  # Save the plotted dendrogram into a JPEG file
  title = paste(brain_region,"_network_dendrogram_and_module_assignment", sep="")
  png(file = file.path(output_directory, paste(title, ".png",sep="")), height = 3.7, width = 5.5, units="in", res=300)
  plotDendroAndColors(geneTree, module_asignment, "Module colors", dendroLabels = FALSE, 
                      hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = title)
  dev.off()
  plotDendroAndColors(geneTree, module_asignment, "Module colors", dendroLabels = FALSE, 
                      hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = title)
  
  # Save the network, detected modules and modules colors into a RData file
  return(module_asignment)
}


#Modified WGCNA plotMat to correct colors for color-blind people
plotMat = function (x, nrgcols = 50, rlabels = FALSE, clabels = FALSE, rcols = 1, ccols = 1, title = "", ...){
  n <- nrow(x)
  p <- ncol(x)
  image(1:p, 1:n, t(x[n:1, ]), col = blueWhiteRed(nrgcols), 
        axes = FALSE, xlab = "", ylab = "", ...)
  if (length(ccols) == 1) {
    axis(3, at = 1:p, labels = clabels, las = 2, cex.axis = 0.6, 
         col.axis = ccols)
  }
  if (length(ccols) == p) {
    cols <- unique(ccols)
    for (i in 1:length(cols)) {
      which <- (1:p)[ccols == cols[i]]
      axis(3, at = which, labels = clabels[which], las = 2, 
           cex.axis = 0.6, col.axis = cols[i])
    }
  }
  if (length(rcols) == 1) {
    axis(2, at = n:1, labels = rlabels, las = 2, cex.axis = 0.6, 
         col.axis = rcols)
  }
  if (length(rcols) == n) {
    cols <- unique(rcols)
    for (i in 1:length(cols)) {
      which <- (1:n)[rcols == cols[i]]
      axis(2, at = (n:1)[which], labels = rlabels[which], 
           las = 2, cex.axis = 0.6, col.axis = cols[i])
    }
  }
  mtext(title, side = 3, line = 3)
  box()
}


plot_eigengenes = function(preprocessed_data, module_asignment, output_directory, brain_region){
  module_eigengenes = moduleEigengenes(preprocessed_data, module_asignment)$eigengenes
  no = 1
  for (color in standardColors(length(table(module_asignment))-1)){
    ME=module_eigengenes[, paste("ME",color, sep="")]
    png(file = paste(output_directory, "/", brain_region, " module ", no, " expression.png", sep = ""), width=8, height=10, units="in", res=300)
    layout(matrix(c(1,2)), heights=c(2,1))
    par(mar=c(0.3, 5.5, 6, 2))
    plotMat(t(scale(preprocessed_data[,module_asignment==color]) ),
            nrgcols=50,rlabels=colnames(preprocessed_data[,module_asignment==color]),clabels = rownames(preprocessed_data), cex.main=2)
    title(main = paste(brain_region, " module ", no, sep = ""), line = 4.5)
    par(mar=c(5, 4.2, 1, 0.7))
    barplot(ME, col=c(rep("lightpink",5), rep("red",5), rep("lightskyblue",5), rep("blue",5), rep("moccasin",5), rep("orange",5)), main="", cex.main=2,
            ylab="eigengene expression",xlab="array sample")
    dev.off()
    no = no+1
  }
}


plot_overlap_between_modules = function(HC_preprocessed, HC_modules, PHC_preprocessed, PHC_modules, output_directory){
  no = 1
  HC_present_prots = lapply(standardColors(length(table(HC_modules))-1), function(color)
    colnames(HC_preprocessed)[which(HC_modules == color)]
  )
  PHC_present_prots = lapply(standardColors(length(table(PHC_modules))-1), function(color)
    colnames(PHC_preprocessed)[which(PHC_modules == color)]
  )
  
  Jaccard_matrix = matrix(0, ncol = (length(table(PHC_modules))-1), nrow = (length(table(HC_modules))-1), 
                          dimnames=list(paste("Module", 1:(length(table(HC_modules))-1)), 
                                        paste("Module", 1:(length(table(PHC_modules))-1))))
  
  for(i in 1:(length(table(HC_modules))-1)){
    list1 = HC_present_prots[[i]]
    for(j in 1:(length(table(PHC_modules))-1)){
      list2 = PHC_present_prots[[j]]
      x = comparelists(list1,list2)
      intersection = x$intersec
      union = union(list1,list2)
      
      Jaccard_index = length(intersection)/length(union)
      Jaccard_matrix[i,j] = round(Jaccard_index, digits = 2)
    }
  }
  png(file = paste(output_directory, "/", "Jaccard_module_overlap.png", sep = ""), height = 6, width = 8, units="in", res=300)
  labeledHeatmap(Matrix = Jaccard_matrix, 
                 xLabels = colnames(Jaccard_matrix),
                 yLabels = rownames(Jaccard_matrix), 
                 colorLabels = F, colors = blueWhiteRed(100)[50:100],
                 cex.text = 0.8, textMatrix = Jaccard_matrix,
                 cex.lab = 0.8, main = "Jaccard similarity between modules", cex.main = 1, xLabelsAngle=25)
  dev.off()
  labeledHeatmap(Matrix = Jaccard_matrix, 
                 xLabels = colnames(Jaccard_matrix),
                 yLabels = rownames(Jaccard_matrix), 
                 colorLabels = F, colors = blueWhiteRed(100)[50:100],
                 cex.text = 0.8, textMatrix = Jaccard_matrix,
                 cex.lab = 0.8, main = "Jaccard similarity between modules", cex.main = 1, xLabelsAngle=25)
}


#Modified function from WGCNA
plotModuleSignificance = function(GeneSignificance, module_asignment, brain_region, boxplot = FALSE, ylab = "Gene Significance", xlab = "Modules", ...){
  no.colors = length(names(table(module_asignment)))
  pp = try(kruskal.test(GeneSignificance, factor(module_asignment))$p.value)
  title = paste(brain_region, " gene significance across modules,", " p-value=", signif(pp, 2), sep = "")
  
  if (boxplot != TRUE) {
    means1 = as.vector(tapply(GeneSignificance, module_asignment, mean, 
                              na.rm = TRUE))
    se1 = as.vector(tapply(GeneSignificance, module_asignment, stdErr))
    merged = matrix(c(means1, se1), ncol = 2)
    rownames(merged) = names(table(module_asignment))
    merged = merged[rownames(merged) != "grey",]
    merged = merged[standardColors(length(table(module_asignment))-1),]
    
    barplot(merged[,1], names.arg = c(1:nrow(merged)), col = rownames(merged), 
            ylab = ylab, xlab = xlab, main = title, ylim = c(0, max(merged[,1])+0.2), ...)
    addErrorBars(merged[,1], as.vector(1.96 * merged[,2]), 
                 two.side = TRUE)
  }
  else {
    boxplot(split(GeneSignificance, module_asignment), notch = T, varwidth = T, 
            col = names(table(module_asignment)), ylab = ylab, main = title, 
            ...)
  }
}


calculate_module_significance = function(trait, preprocessed_data, module_asignment, output_directory, brain_region){
  
  GeneSignificance = abs(as.numeric(cor(trait, preprocessed_data, use="p", method="spearman")))
  
  # Module significance is defined as average gene significance.
  ModuleSignificance=tapply(GeneSignificance, module_asignment, mean, na.rm=T)
  
  #Plot module significance
  png(file = paste(output_directory, "/", brain_region, " module significance.png", sep = ""), width=6, height=5, units="in", res=300)
  plotModuleSignificance(GeneSignificance, module_asignment, brain_region)
  dev.off()
  
  return(GeneSignificance)
}



get_module_trait_relationship = function(preprocessed_data, module_asignment, output_directory, time_label, brain_region){
  nGenes = ncol(preprocessed_data)
  nSamples = nrow(preprocessed_data)
  trait = matrix(c(rep(0,5), rep(1,5), rep(NA, 30), rep(0,5), rep(1,5), rep(NA, 30), rep(0,5), rep(1,5)), ncol=3, dimnames=list(rownames(preprocessed_data), time_label))
  
  #Calculate module Eigengenes
  MEs0 = moduleEigengenes(preprocessed_data, module_asignment, excludeGrey = TRUE)$eigengenes
  
  for(i in 1:ncol(MEs0)){
    names(MEs0)[i] = paste("Module", which(standardColors(length(table(module_asignment))-1) == strsplit(names(MEs0)[i], "ME")[[1]][2]))
  }
  MEs = MEs0[,order(as.numeric(gsub(".* ", "", names(MEs0))))]
  
  plot_module_trait_relationship(trait, MEs, nSamples, "time points separated", output_directory, brain_region)
}


plot_module_trait_relationship = function(trait, MEs, nSamples, type, output_directory, brain_region){
  moduleTraitCor = cor(MEs, trait, use = "p", method="spearman")
  moduleTraitCor = moduleTraitCor[which(rownames(moduleTraitCor) != "MEgrey"),]
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  textMatrix = paste(signif(moduleTraitCor, 2), "\\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  png(file = paste(output_directory, "/", brain_region, " module-trait relationship ", type, ".png", sep = ""), width=7, height=9, units="in", res=300)
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(trait),
                 yLabels = paste0("ME", standardColors(nrow(moduleTraitCor))),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = TRUE,
                 cex.text = 1,
                 zlim = c(-1,1),
                 main = paste(brain_region, "module-trait relationships"))
  dev.off()
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(trait),
                 yLabels = paste0("ME", standardColors(nrow(moduleTraitCor))),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = TRUE,
                 cex.text = 0.9,
                 zlim = c(-1,1),
                 main = paste(brain_region, "module-trait relationships"))
}


get_intramodular_connectivity = function(adjacency_matrix, module_asignment, output_directory, brain_region){
  protein_intramodular_connectivity = intramodularConnectivity(adjacency_matrix, module_asignment, scaleByMax = TRUE)[,-1]
  protein_intramodular_connectivity[,2] = module_asignment
  hubs = lapply(1:(length(table(module_asignment))-1), function(module_no){
    color = standardColors(length(table(module_asignment))-1)[module_no]
    filtered_module = protein_intramodular_connectivity[protein_intramodular_connectivity[,2] == color,]
    filtered_module = filtered_module[order(filtered_module$kWithin, decreasing = TRUE),]
    filtered_module = filtered_module["kWithin"]
    colnames(filtered_module) = paste(brain_region, "module", module_no)
    print(head(filtered_module))
    write.table(filtered_module, file = paste(output_directory, "/", brain_region, " module ", module_no, " intramodular hubs.txt", sep = ""), row.names = T, col.names = F, sep = "\\t", quote = FALSE)
    return(rownames(filtered_module)[1:15])
  })
  return(hubs)
}


compare_hubs_with_PCA = function(hubs, time_label, brain_regions){
  hub_table = NULL
  for (i in 1:length(hubs)){
    brain_region = brain_regions[i]
    sub = data.frame(paste(brain_region, rep(1:length(hubs[[i]]), each=15)), row.names=unlist(hubs[i]))
    hub_table = rbind(hub_table, sub)
  }
  for (brain_region in brain_regions){
    for (time in time_label){
      prots = rep(0, nrow(hub_table))
      file = paste(brain_region, time, "top10_PCA_proteins.txt", sep="_")
      PC1 = read.table(file.path(output_directory, file), sep="\\t", header=TRUE)[,1]
      hub_table = cbind(hub_table, i=ifelse(rownames(hub_table) %in% PC1==TRUE, "+", "-"))
    }
  }
  colnames(hub_table) = c("Module", paste(rep(brain_regions, each=length(time_label)), time_label))
  write.table(hub_table, file = file.path(output_directory, "Hubs_and_PC1.txt"), row.names = T, col.names = T, sep = "\\t", quote = FALSE)
  return(hub_table)
}