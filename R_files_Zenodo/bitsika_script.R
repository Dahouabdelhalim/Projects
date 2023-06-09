########################################################################################
## WORKFLOW TO ANALYSE BITSIKA DATA
##
## These script contains a workflow used to reproduce results generated in paper
## "A systems level analysis of epileptogenesis-associated proteome alterations" by
## M. Keck and G. Androsova et al. (2016).
## 
## Script author: Ganna Androsova, ganna.androsova@uni.lu
## 
########################################################################################
source("Data_analysis_and_network_construction.R")

#General settings
options(stringsAsFactors = FALSE)

#Indicate the directory with original data
original_directory = "/Original_data"

#Indicate the directory for input files
input_directory = "/Preprocessed_data"

#Indicate the output directory
output_directory = "/Generated_results"

#Indicate time points of the data
time_label = c("1 day", "3 days", "30 days")

#Indicate the names of sheets to read from the original Bitsika's excel files
sheets_to_read = c("1dpi-all proteins", "3dpi-all proteins", "30dpi-all proteins")

################################### Data preprocessing ################################
#Read-in original data
Bitsika_original_data = lapply(sheets_to_read, function(sheet_name){
  input_file = read.xlsx(file.path(original_directory, "pr6b00003_si_005.xlsx"), 
                         sheet = sheet_name, colNames = TRUE, startRow = 2)
  #Transform excel output into data frame and keep only expression values
  df = data.frame(input_file[which(unique(input_file[,2]) == input_file[,2]),c(12:16,7:11)])
  rownames(df) = input_file[which(unique(input_file[,2]) == input_file[,2]),2]
  return(df)
})

#Merge the time points
Bitsika_merged_data = Bitsika_original_data[[1]]
for(i in c(2:length(Bitsika_original_data))){
  merge1 = merge(Bitsika_merged_data, Bitsika_original_data[[i]], 
                 by="row.names", all = FALSE)
  Bitsika_merged_data = merge1[,2:ncol(merge1)]
  rownames(Bitsika_merged_data) = merge1[,1]
}
colnames(Bitsika_merged_data) = c(rep("NaCl_1d", 5), rep("KA_1d",5), 
                                  rep("NaCl_3d", 5), rep("KA_3d",5), 
                                  rep("NaCl_30d", 5), rep("KA_30d",5))

png(file = paste(output_directory, "/Bitsika_expression_before_norm.png", sep = ""), 
    height = 5, width = 7, units="in", res=300)
boxplot(Bitsika_merged_data, xlab="Samples", ylab="Protein abundances", 
        main="Protein expression values before normalization")
dev.off()

#PCA for each seperate timepoint with all detected proteins
plot_separate_PCA(Bitsika_original_data, output_directory, "Bitsika", time_label)

#Normalize the data
Bitsika_preprocessed = normalization(Bitsika_merged_data, input_directory, 
                                   output_directory, "Bitsika")

#Get protein expression values
Bitsika_post_stat = DE_prots(Bitsika_preprocessed, output_directory, time_label, "Bitsika")

#Gen names of differentially expressed proteins
Bitsika_list_of_DE_prots = lapply(Bitsika_post_stat, function(x){
  prot=subset(x, q_values<.05 & (x$FC<=0.67|x$FC>=1.5))
  return(rownames(prot))
})

#Calculate the overlap of DE proteins between Bitsika data and HC/PHC
DE_overlap = matrix(0, ncol = 6, nrow = 3, dimnames = list(
  c(paste("Bitsika", time_label)), 
  c("HC 2 days", "HC 10 days", "HC 8 weeks", "PHC 2 days", "PHC 10 days", "PHC 8 weeks")))
for(row in c(1:3)){
  for(col in c(1:3)){
    DE_overlap[row,col] = length(Reduce(intersect, list(Bitsika_list_of_DE_prots[[row]], HC_list_of_DE_prots[[col]])))
    DE_overlap[row,col+3] = length(Reduce(intersect, list(Bitsika_list_of_DE_prots[[row]], PHC_list_of_DE_prots[[col]])))
  }
}

#Plot the overlap of DE proteins between Bitsika data and HC/PHC
png(file = paste(output_directory, "/", "DE_protein_overlap.png", sep = ""), height = 3.5, width = 7, units="in", res=300)
labeledHeatmap(Matrix = DE_overlap, 
               xLabels = colnames(DE_overlap),
               yLabels = rownames(DE_overlap), 
               colorLabels = F, colors = blueWhiteRed(100)[50:100],
               cex.text = 0.8, textMatrix = DE_overlap,
               cex.lab = 0.8, xLabelsAngle=25)
dev.off()

#Plot the differentially expressed proteins
barplot_DE(Bitsika_post_stat, output_directory, "Bitsika")

#Plot the correlation matrix
Bitsika_correlation_matrix = cor(t(Bitsika_preprocessed), 
                                 method="spearman", use="complete.obs")
visualize_correlation(Bitsika_correlation_matrix, output_directory, "Bitsika")

#Create adjacency matrix
Bitsika_adjacency_matrix = estimate_correlation(Bitsika_preprocessed)

#Calculate network topology
Bitsika_topological_table = estimate_network_topology(Bitsika_adjacency_matrix, 
                                                      output_directory, "Bitsika")

#Module detection
Bitsika_modules = detect_modules(Bitsika_adjacency_matrix, output_directory, "Bitsika")

#Export the network to Cytoscape
exportNetworkToCytoscape(Bitsika_adjacency_matrix, weighted = TRUE, threshold = 0.26, 
                         nodeNames = colnames(Bitsika_preprocessed), 
                         nodeAttr = Bitsika_modules,
                         edgeFile = file.path(output_directory, "Bitsika_network_edges.txt"),
                         nodeFile = file.path(output_directory, "Bitsika_network_nodes.txt"))

#Plot module eigengenes
plot_eigengenes(Bitsika_preprocessed, Bitsika_modules, output_directory, "Bitsika")

#Module-trait relationship
#Create the trait matrix
trait = matrix(rep(c(rep(0,5), rep(1,5)),3))
rownames(trait) = rownames(Bitsika_preprocessed)

#Get protein significance
Bitsika_protein_significance = calculate_module_significance(trait, Bitsika_preprocessed, Bitsika_modules, output_directory, "Bitsika")

#Plot figures of module significance
plotModuleSignificance(Bitsika_protein_significance, Bitsika_modules, "Bitsika")

#Calculate module-trail relationship
get_module_trait_relationship(Bitsika_preprocessed, Bitsika_modules, output_directory, time_label, "Bitsika")

#Intramodular hub proteins
#Get tables with intramodular connectivity of the proteins
Bitsika_hubs = get_intramodular_connectivity(Bitsika_adjacency_matrix, Bitsika_modules, output_directory, "Bitsika")

#Compare intramodular hubs to PCA contributing proteins
compare_hubs_with_PCA(list(Bitsika_hubs), time_label, c("Bitsika"))
