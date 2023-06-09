################################################################################ 

# SCRIPT 5: NODE-SPECIFIC RANDOM FOREST ANALYSIS PREDICTORS OF ANCESTRAL STATE
# ESTIMATION ACCURACY, PRECISION, AND SENSITIVITY TO BRANCH LENGTH CHOICE.

# AUTHOR: Nicolás Mongiardino Koch

# ARTICLE: Wilson et al. 2022 - Chronogram or phylogram for ancestral state
# estimation? Model-fit statistics indicate the branch lengths underlying a
# binary character’s evolution

################################################################################


library(phangorn)
library(castor)
library(randomForest)
library(tidyverse)
library(phytools)
load('trees&characters&metrics.Rdata')


# GENERATING PREDICTOR AND RESPONSE VARIABLES-----------------------------------

for(i in 1:nrow(data)) {
  this_phylo = data$phylogram[[i]]
  this_chrono = data$chronogram[[i]]
  
  root_tip_chrono = max(node.depth.edgelength(this_chrono))
  root_tip_phylo = max(node.depth.edgelength(this_phylo))
  
  desc = node.depth(this_phylo, method = 1)
  desc = desc[(length(this_phylo$tip.label)+1):length(desc)]
  desc_prop = desc/length(this_phylo$tip.label)
  
  node_depth_chrono = node.depth.edgelength(this_chrono)/root_tip_chrono
  node_depth_chrono = node_depth_chrono[(length(this_phylo$tip.label)+1):length(node_depth_chrono)]
  
  this_tree_data = data.frame(tree = i, node = 1:this_phylo$Nnode, 
                              desc, desc_prop, node_depth_chrono)
  
  change_branches_rel = branches_phylo = branches_chrono = 
    mean_dist_phylo = min_dist_phylo = diff_dist_mean = 
    diff_dist_min = vector(length = this_phylo$Nnode)
  
  for(j in 1:this_phylo$Nnode) {
    
    desc_branches = which(this_phylo$edge[,1] == (length(this_phylo$tip.label)+j))
    this_branches_phylo = this_phylo$edge.length[desc_branches]
    this_branches_chrono = this_chrono$edge.length[desc_branches]
    
    if(j != 1) {
      anc_branch = which(this_phylo$edge[,2] == (length(this_phylo$tip.label)+j))
      this_branches_phylo = c(this_branches_phylo, this_phylo$edge.length[anc_branch])
      this_branches_chrono = c(this_branches_chrono, this_chrono$edge.length[anc_branch])
    }
    
    branches_phylo[j] = sum(this_branches_phylo)/root_tip_phylo
    branches_chrono[j] = sum(this_branches_chrono)/root_tip_chrono
    change_branches_rel[j] = sum(abs(this_branches_phylo - this_branches_chrono)) / 
      sum(c(sum(this_branches_phylo),sum(this_branches_chrono)))
    
    desc_tips = unlist(Descendants(this_phylo, (length(this_phylo$tip.label)+j), 
                                   type = 'tips'))
    dist_phylo = get_pairwise_distances(this_phylo, 
                                        rep((length(this_phylo$tip.label)+j), 
                                            length(desc_tips)), 
                                        desc_tips)
    
    dist_chrono = get_pairwise_distances(this_chrono, 
                                         (length(this_phylo$tip.label)+j), 
                                         desc_tips[1])
    
    mean_dist_phylo[j] = mean(dist_phylo/root_tip_phylo)
    min_dist_phylo[j] = min(dist_phylo/root_tip_phylo)
    
    diff_dist_mean[j] = abs((mean(dist_phylo) - dist_chrono)/sum(c(mean(dist_phylo), dist_chrono)))
    diff_dist_min[j] = abs((min(dist_phylo) - dist_chrono)/sum(c(min(dist_phylo), dist_chrono)))
  }
  
  chrono = F
  if(data$character_branch_lengths[i] == 'chronogram') chrono = T
  
  truth = data$node_states[i][[1]]
  truth = sapply(1:2, "==", truth) * 1
  estimated_chrono = data$chrono_ARD_marginals[i][[1]]
  estimated_phylo = data$phylo_ARD_marginals[i][[1]]
  
  if(chrono) {
    accuracy = 1 - abs(truth - estimated_chrono)[,1]
    precision = abs(matrix(0.5, ncol = 2, nrow = nrow(estimated_chrono)) - estimated_chrono)[,1]
  } else {
    accuracy = 1 - abs(truth - estimated_phylo)[,1]
    precision = abs(matrix(0.5, ncol = 2, nrow = nrow(estimated_phylo)) - estimated_phylo)[,1]
  }
  
  precision = precision/0.5
  sensitivity = abs(estimated_chrono - estimated_phylo)[,1]
  
  this_tree_data = cbind(this_tree_data, branches_phylo, branches_chrono, 
                         change_branches_rel, mean_dist_phylo, 
                         min_dist_phylo, diff_dist_mean, diff_dist_min, 
                         accuracy, precision, sensitivity)
  
  if(i == 1) {
    all_node_prop = this_tree_data
  } else {
    all_node_prop = rbind(all_node_prop, this_tree_data)
  }
}

save(all_node_prop, file = "node_property_dataset.Rdata")


# RUN RANDOM FOREST ANALYSES----------------------------------------------------


ntrees = 1000

accuracy_forest = randomForest(x = all_node_prop[,3:14], y = all_node_prop$accuracy, 
                               nodesize = 5, sampsize = floor(nrow(all_node_prop) * 0.01), 
                               data = all_node_prop[, 3:15], 
                               method = anova, proximity = F, 
                               ntree = ntrees, imtest_BItance = TRUE,
                               type = 'regression', replace = TRUE, importance = T)

# Print the model output                             
imp_accuracy = as.data.frame(round(importance(accuracy_forest, conditional = T), 2), optional = T)
imp_accuracy$variable = rownames(imp_accuracy)
imp_accuracy$type = 'accuracy'

precision_forest = randomForest(x = all_node_prop[,3:14], y = all_node_prop$precision, 
                                nodesize = 5, sampsize = floor(nrow(all_node_prop) * 0.01), 
                                data = all_node_prop[, c(3:14, 16)], 
                                method = anova, proximity = F, 
                                ntree = ntrees, imtest_BItance = TRUE,
                                type = 'regression', replace = TRUE, importance = T)

# Print the model output                             
imp_precision = as.data.frame(round(importance(precision_forest, conditional = T), 2), optional = T)
imp_precision$variable = rownames(imp_precision)
imp_precision$type = 'precision'

sensitivity_forest = randomForest(x = all_node_prop[,3:14], y = all_node_prop$sensitivity, 
                                  nodesize = 5, sampsize = floor(nrow(all_node_prop) * 0.01), 
                                  data = all_node_prop[, c(3:14, 17)], 
                                  method = anova, proximity = F, 
                                  ntree = ntrees, imtest_BItance = TRUE,
                                  type = 'regression', replace = TRUE, importance = T)

# Print the model output                             
imp_sensitivity = as.data.frame(round(importance(sensitivity_forest, conditional = T), 2), optional = T)
imp_sensitivity$variable = rownames(imp_sensitivity)
imp_sensitivity$type = 'sensitivity'

rownames(imp_accuracy) = rownames(imp_precision) = 
  rownames(imp_sensitivity) = NULL

results = tibble(rbind(imp_precision, imp_accuracy, imp_sensitivity))[,2:4]

results$color_Purity = c((imp_precision$IncNodePurity - min(imp_precision$IncNodePurity))/max(imp_precision$IncNodePurity - min(imp_precision$IncNodePurity)), 
                         (imp_accuracy$IncNodePurity - min(imp_accuracy$IncNodePurity))/max(imp_accuracy$IncNodePurity - min(imp_accuracy$IncNodePurity)), 
                         (imp_sensitivity$IncNodePurity - min(imp_sensitivity$IncNodePurity))/max(imp_sensitivity$IncNodePurity - min(imp_sensitivity$IncNodePurity)))

results = results %>% group_by(type) %>% arrange(desc(IncNodePurity), .by_group = T) %>% ungroup()
results$type = factor(results$type, levels = unique(results$type)[c(2, 1, 3)])

#Plot figure 5
ggplot(results, aes(x = factor(variable, levels = variable[10:1]), y = IncNodePurity)) + 
  geom_bar(stat='identity', aes(fill = color_Purity), color='black') + 
  scale_fill_gradient(high = met.brewer(name="Hokusai1",type="discrete")[3], low = 'white') +
  coord_flip() + theme(legend.position="none") +
  facet_wrap(~ type, scales = 'free_x')