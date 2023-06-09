
library(statnet)
library(GGally)
library(ggplot2)
library(igraph)
library(tnet)
library(corrplot)

# read in data

law_oc_mat <- as.matrix(read.csv("law_occurence.csv", row.names = 1)) # load occurence of issues in laws
proj_oc_mat <- as.matrix(read.csv("issue_occurence.csv", row.names = 1)) # load occurence of issues in actor activity

# Issue association strength based on law overlap

# ochiai similiarity

issue_issue_law_mat <- t(law_oc_mat) %*% law_oc_mat

colsums_x_rowsums_mat_docs <- outer(colSums(issue_issue_law_mat), rowSums(issue_issue_law_mat))
issue_issue_law_mat_norm <- (issue_issue_law_mat / sqrt(colsums_x_rowsums_mat_docs))
diag(issue_issue_law_mat_norm) <- 0 # remove self-ties

# visualize

issue_issue_law_net <- network(issue_issue_law_mat_norm,
                               ignore.eval = F,
                               names.eval = "weights")
coords <- gplot.layout.fruchtermanreingold(issue_issue_law_net, NULL)

issue_issue_law_viz <-
  ggnet2(issue_issue_law_net,
         mode = coords,
         node.label = T,
         edge.size = network::get.edge.attribute(issue_issue_law_net,"weights")*10,
         label.size = 3, node.alpha = 0.7, node.color = "#6698FF",
         size = "degree") +
  ggtitle("Issue association based on law co-occurence") +
  theme(legend.position="none")
issue_issue_law_viz

# Issue association strength based on actor overlap

# actor_projection based on newman actor_projection

el_bip_actor_issue <- as_edgelist(graph.incidence(proj_oc_mat,directed = F))
tnet_rep <- as.tnet(el_bip_actor_issue)

actor_projection <- projecting_tm(net = cbind(tnet_rep$p, tnet_rep$i),method = "Newman")

el_bip_actor_issue_newman <- cbind(as.character(factor(actor_projection$i, labels = levels(tnet_rep$p))),
                                   as.character(factor(actor_projection$j, labels = levels(tnet_rep$p))),
                                   as.numeric(actor_projection$w))


unique_issues_1 <- levels(tnet_rep$p)

issue_issue_actors_mat <- matrix(0,nrow = length(unique_issues_1), ncol = length(unique_issues_1))
colnames(issue_issue_actors_mat) <- unique_issues_1
rownames(issue_issue_actors_mat) <- unique_issues_1

issue_issue_actors_mat[el_bip_actor_issue_newman[,1:2]] <- as.numeric(as.character(el_bip_actor_issue_newman[,3]))

issue_issue_mat_actors_norm <- issue_issue_actors_mat

# visualize

#reorder based on ordering of law association mat:

issue_issue_actor_viz_net <- network(issue_issue_mat_actors_norm,
                                     ignore.eval = F,
                                     names.eval = "weights")

issue_issue_actor_viz <-
  ggnet2(issue_issue_actor_viz_net,
         mode = coords,
         node.label = T,
         edge.size = as.numeric(network::get.edge.attribute(issue_issue_actor_viz_net,"weights"))/10,
         label.size = 3, node.alpha = 0.7, node.color = "#6698FF",
         node.size = sna::degree(issue_issue_actor_viz_net)) +
  ggtitle("Issue association based on actor activity") +
  theme(legend.position="none")
issue_issue_actor_viz

issue_issue_mat_actors_norm <- # same ordering for both measures
  issue_issue_mat_actors_norm[rownames(issue_issue_law_mat_norm),colnames(issue_issue_law_mat_norm)]

# compare measures to identify misfit

# make a graph out of both matrices, get edgelists and create combined edgelist

# for laws
issue_issue_law_graph <- graph.adjacency(issue_issue_law_mat_norm,
                                         weighted = T,
                                         mode = "undirected")
# for actors
issue_issue_actors_graph <- graph.adjacency(issue_issue_mat_actors_norm,
                                            weighted = T,
                                            mode = "undirected")
# combine the two
combined_graph <- union(issue_issue_actors_graph, issue_issue_law_graph, byname = T)

issue_issue_el <- igraph::get.edgelist(combined_graph)
actor_overlap <- igraph::get.edge.attribute(combined_graph,"weight_1")
actor_overlap[is.na(actor_overlap)] <- 0
docs_overlap <- igraph::get.edge.attribute(combined_graph,"weight_2")
docs_overlap[is.na(docs_overlap)] <- 0
issue_issue_el <- igraph::get.edgelist(combined_graph)
issue_issue_el <- data.frame(cbind(issue_issue_el,
                                   actor_overlap,
                                   docs_overlap), stringsAsFactors = F)
colnames(issue_issue_el) <- c("sender", "receiver","actor_overlap","docs_overlap")
issue_issue_el$actor_overlap <- as.numeric(as.character(issue_issue_el$actor_overlap))
issue_issue_el$docs_overlap <- as.numeric(as.character(issue_issue_el$docs_overlap))
head(issue_issue_el)

# scale both between 0 and 1

scale_01 <- function(x){(x-min(x))/(max(x)-min(x))}

issue_issue_el$actor_overlap_scaled <- scale_01(issue_issue_el$actor_overlap)
issue_issue_el$docs_overlap_scaled <- scale_01(issue_issue_el$docs_overlap)

# create fit measure

issue_issue_el$fit_measure <- issue_issue_el$actor_overlap_scaled - issue_issue_el$docs_overlap_scaled

# make heatmaps

unique_issues <- unique(c(as.character(issue_issue_el$sender), as.character(issue_issue_el$receiver)))
fit_heatmap <- matrix(0,nrow = length(unique_issues), ncol = length(unique_issues))
colnames(fit_heatmap) <- unique_issues
rownames(fit_heatmap) <- unique_issues
fit_heatmap[cbind(issue_issue_el$sender,issue_issue_el$receiver)] <- issue_issue_el$fit_measure
#symmetrize matrix
fit_heatmap[lower.tri(fit_heatmap)] <- fit_heatmap[lower.tri(fit_heatmap)] + t(fit_heatmap)[lower.tri(t(fit_heatmap))]
fit_heatmap[upper.tri(fit_heatmap)] <- t(fit_heatmap)[upper.tri(t(fit_heatmap))]
isSymmetric(fit_heatmap)

heatmap <-
  corrplot(fit_heatmap, type = "lower", method = "color",
           is.corr = F, order = "original",
           tl.srt = 90,tl.cex = 1,
           mar = c(0, 0, 3, 3),
           tl.col = "black", diag = T)
heatmap
