########################################################################################
## WORKFLOW TO REPRODUCE MANUSCRIPT RESULTS
##
## These script contains a workflow used to reproduce results generated in paper
## "A systems level analysis of epileptogenesis-associated proteome alterations" by
## M. Keck and G. Androsova et al. (2016).
## 
## Script author: Ganna Androsova, ganna.androsova@uni.lu
## 
########################################################################################

### Call scripts with functions and configurations
source("Data_analysis_and_network_construction.R")
source("Configurations.R")

################################### Data preprocessing #################################
#Read-in original data
HC_data = merge_original_tables(original_directory, HC_excel_files, "HC")
PHC_data = merge_original_tables(original_directory, PHC_excel_files, "PHC")

#PCA for each seperate timepoint with all detected proteins
plot_separate_PCA(HC_data[[2]], output_directory, "HC", time_label)
plot_separate_PCA(PHC_data[[2]], output_directory, "PHC", time_label)

#Arcsinh transformation of median-centering
HC_preprocessed_all = normalization(HC_data[[1]], input_directory, 
                                    output_directory, "HC")
PHC_preprocessed_all = normalization(PHC_data[[1]], input_directory, 
                                     output_directory, "PHC")

#Identify differentially expressed proteins
HC_DEGs = DE_prots(HC_preprocessed_all, output_directory, time_label, "HC")
PHC_DEGs = DE_prots(PHC_preprocessed_all, output_directory, time_label, "PHC")

HC_list_of_DE_prots = lapply(HC_DEGs, function(x){
  prot=subset(x, p_values<.05 & (x$FC<=0.67|x$FC>=1.5))
  return(rownames(prot))
})
PHC_list_of_DE_prots = lapply(PHC_DEGs, function(x){
  prot=subset(x, p_values<.05 & (x$FC<=0.67|x$FC>=1.5))
  return(rownames(prot))
})

hist_FC(HC_DEGs, output_directory, "HC")
hist_FC(PHC_DEGs, output_directory, "PHC")

volcano_DE(HC_DEGs, output_directory, "HC")
volcano_DE(PHC_DEGs, output_directory, "PHC")

barplot_DE(HC_DEGs, output_directory, "HC")
barplot_DE(PHC_DEGs, output_directory, "PHC")

#Get the merged datasets with proteins present at all timepoints
HC_merged_data = get_common_prots(HC_data[[2]], input_directory, output_directory, "HC_common")
PHC_merged_data = get_common_prots(PHC_data[[2]], input_directory, output_directory, "PHC_common")

#Identify differential expression among common proteins
DE_prots(HC_merged_data, output_directory, time_label, "HC_common")
DE_prots(PHC_merged_data, output_directory, time_label, "PHC_common")

HC_correlation_matrix = cor(t(HC_merged_data), method="spearman", use="complete.obs")
PHC_correlation_matrix = cor(t(HC_merged_data), method="spearman", use="complete.obs")

visualize_correlation(HC_correlation_matrix, output_directory, "HC")
visualize_correlation(PHC_correlation_matrix, output_directory, "PHC")

############################## Network construction ####################################

#Create adjacency matrix
HC_adjacency_matrix = estimate_correlation(HC_merged_data)
PHC_adjacency_matrix = estimate_correlation(PHC_merged_data)

#Calculate network topology
HC_topological_table = estimate_network_topology(HC_adjacency_matrix, output_directory, "HC")
PHC_topological_table = estimate_network_topology(PHC_adjacency_matrix, output_directory, "PHC")
topology = rbind(HC_topological_table, PHC_topological_table)
rownames(topology) = c("Hippocampus", "Parahippocampus")

#Module detection
HC_modules = detect_modules(HC_adjacency_matrix, output_directory, "HC")
PHC_modules = detect_modules(PHC_adjacency_matrix, output_directory, "PHC")

#Export network to Cytoscape
exportNetworkToCytoscape(HC_adjacency_matrix, weighted = TRUE, threshold = 0.26, 
                         nodeNames = colnames(HC_adjacency_matrix), nodeAttr = HC_modules,
                         edgeFile = file.path(output_directory, "HC_network_edges.txt"),
                         nodeFile = file.path(output_directory, "HC_node_color_attribute.txt"))
exportNetworkToCytoscape(PHC_adjacency_matrix, weighted = TRUE, threshold = 0.26, 
                         nodeNames = colnames(PHC_adjacency_matrix), nodeAttr = PHC_modules,
                         edgeFile = file.path(output_directory, "PHC_network_edges.txt"),
                         nodeFile = file.path(output_directory, "PHC_node_color_attribute.txt"))

#Plot module eigengenes
plot_eigengenes(HC_merged_data, HC_modules, output_directory, "HC")
plot_eigengenes(PHC_merged_data, PHC_modules, output_directory, "PHC")

#Module overlap
plot_overlap_between_modules(HC_merged_data, HC_modules, PHC_merged_data, PHC_modules, output_directory)

############################## Module-trait relationship ##############################
#Create the trait matrix
trait = matrix(rep(c(rep(0,5), rep(1,5)),3))
rownames(trait) = rownames(HC_merged_data)
source("Data_analysis_and_network_construction.R")
#Get protein significance
HC_protein_significance = calculate_module_significance(trait, HC_merged_data, 
                                                        HC_modules, output_directory, "HC")
rownames(trait) = rownames(PHC_merged_data)
PHC_protein_significance = calculate_module_significance(trait, PHC_merged_data, 
                                                         PHC_modules, output_directory, "PHC")

#Calculate module-trail relationship
get_module_trait_relationship(HC_merged_data, HC_modules, 
                              output_directory, time_label, "HC")
get_module_trait_relationship(PHC_merged_data, PHC_modules, 
                              output_directory, time_label, "PHC")

#Intramodular hub proteins
#Get tables with intramodular connectivity of the proteins
HC_hubs = get_intramodular_connectivity(HC_adjacency_matrix, HC_modules, output_directory, "HC")
PHC_hubs = get_intramodular_connectivity(PHC_adjacency_matrix, PHC_modules, output_directory, "PHC")

#Compare intramodular hubs to PCA contributing proteins
hubs_vs_PCA = compare_hubs_with_PCA(list(HC_hubs, PHC_hubs), time_label, c("HC", "PHC"))

######################### Neurological disease association ########################
diseases_filtered = read.xlsx(file.path(output_directory, "DisGeNET_results_filtered.xlsx"), 
                     sheet = 1, colNames = TRUE)

ggplot(diseases_filtered, aes(Symbol, Disease.Name)) + 
  geom_tile(aes(fill = Score), colour = "white") + 
  scale_fill_gradient(low = "pink", high = "red") + 
  theme_grey(base_size = 9) + 
  labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(filename = file.path(output_directory, "Disease-associations_neuro_filtered.png"), 
       scale = 1, width = par("din")[1], height = par("din")[2], units = "in")