#install.packages('igraph')
#install.packages('ggplot2')
#install.packages('reshape2')
#install.packages('ggplot')
library('igraph')
# create network visualisation with igraph
#install.packages("data.table")

#install.packages("RColorBrewer")
library(RColorBrewer)

# read in my own data

#edges <- read.csv(file="P:/thesis/thesis_final_data/produced_data/clusters_network/clusters_edges.csv", header=TRUE, sep=",")

#edges <- read.csv(file="P:/thesis/thesis_final_data/produced_data/clusters_network_nodes_edges_ch7/all_SC_cat_norm_cosine_similarity/clusters_edges_all_SC_cat_norm_cosine_similarity.csv", header=TRUE, sep=",")
edges <- read.csv(file="P:/thesis/thesis_final_data/produced_data/clusters_network_nodes_edges_ch7/all_SC_cat_norm_waltmann_etal/clusters_edges_all_SC_cat_norm_waltmann_etal.csv", header=TRUE, sep=",")
#edges <- read.csv(file="P:/thesis/thesis_final_data/produced_data/clusters_network_nodes_edges_ch7/comp_plus_soc_norm_cosine_similarity/clusters_edges_comp_plus_soc_norm_cosine_similarity.csv", header=TRUE, sep=",")
#edges <- read.csv(file="P:/thesis/thesis_final_data/produced_data/clusters_network_nodes_edges_ch7/comp_plus_soc_norm_waltmann_etal/clusters_edges_comp_plus_soc_norm_waltmann_etal.csv", header=TRUE, sep=",")

edges  <- subset(edges, select = -c(X))
edges  <- subset(edges, select = -c(index))
str(edges)

#nodes <- read.csv(file="P:/thesis/thesis_final_data/produced_data/clusters_network/clusters_nodes.csv", header=TRUE, sep=",")

#nodes <- read.csv(file="P:/thesis/thesis_final_data/produced_data/clusters_network_nodes_edges_ch7/all_SC_cat_norm_cosine_similarity/clusters_nodes_all_SC_cat_norm_cosine_similarity.csv", header=TRUE, sep=",")
nodes <- read.csv(file="P:/thesis/thesis_final_data/produced_data/clusters_network_nodes_edges_ch7/all_SC_cat_norm_waltmann_etal/clusters_nodes_all_SC_cat_norm_waltmann_etal.csv", header=TRUE, sep=",")
#nodes <- read.csv(file="P:/thesis/thesis_final_data/produced_data/clusters_network_nodes_edges_ch7/comp_plus_soc_norm_cosine_similarity/clusters_nodes_comp_plus_soc_norm_cosine_similarity.csv", header=TRUE, sep=",")
#nodes <- read.csv(file="P:/thesis/thesis_final_data/produced_data/clusters_network_nodes_edges_ch7/comp_plus_soc_norm_waltmann_etal/clusters_nodes_comp_plus_soc_norm_waltmann_etal.csv", header=TRUE, sep=",")


nodes  <- subset(nodes, select = -c(X))
str(nodes)
nodes

# colors
cols <- colorRampPalette(brewer.pal(8,"Dark2"))(8)
pal2 <- rainbow(8, alpha=.5)  
plot(x=1:10, y=1:10, pch=19, cex=5, col=pal2)

# create igraph object

net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
vertex_attr(net, 'size', index = V(net)) <- nodes$size / 10
vertex_attr(net, 'size', index = V(net)) <- nodes$size / 15
edge_attr(net, 'weight', index = E(net)) <- edges$weight

#vertex_attr(net, 'color', index = V(net)) <- c("gray50", "tomato", "gold", "blue", "green", "orange", "blue", "orange")

pal1 <- heat.colors(8, alpha=.5) 
pal2 <- rainbow(8, alpha=.5)  
vertex_attr(net, 'color', index = V(net)) <- pal2

# plot

plot(net)


# UNUSED CODE


my_layout = layout_with_kk
my_layout = layout_with_mds
my_layout = layout_with_lgl
my_layout = layout_with_dri
my_layout = layout_on_sphere
my_layout = layout.lgl
my_layout = layout.reingold.tilford(g, circular=T)
my_layout = component_wise
my_layout = layout_components
my_layout = layout_as_tree
my_layout = layout_with_dh

# 
#plot(net, edge.width = E(net)$weight, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)

# vertex.label=NA,
E(net)$weight

# plot only connnected component
dg <- decompose.graph(net) # returns a list of three graphs
dg

# set color
display.brewer.all(5)

library('colorspace')
pal1 <- sequential_hcl(12)
pal2 <- rainbow(12, alpha=.5)
pal3 <- heat_hcl(12, h = c(0, -100), l = c(75, 40), c = c(40, 80), power = 1)
pal3
pal4 <- diverge_hcl(12, c = 100, l = c(50, 90), power = 1)

vertex_attr(dg[[1]], 'color', index = V(dg[[1]])) <- pal4

# plot
my_layout = layout_with_dh
plot(dg[[1]], edge.width = E(net)$weight, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)

nodes_in_connected_comp = V(dg[[1]])
nodes_in_connected_comp

# plot subject categories in each connected community
# https://stackoverflow.com/questions/23971974/single-barplot-for-each-row-of-dataframe

# nodes_in_connected_comp
# https://igraph.org/r/doc/as_ids.html

str(nodes)
nodes_in_connected_comp
a <- as.numeric(as_ids(nodes_in_connected_comp))
print(a)

connected_nodes <- nodes[is.element(nodes$Id, a),]
str(connected_nodes)

# transpose data frame
library(reshape2)
str(nodes)
m <- melt(connected_nodes, id = 'Id')
m

# or all plots in one window
ggplot(m, aes(variable, value,  fill = variable)) + 
  facet_wrap(~ Id) +
  geom_bar(stat="identity", show.legend=FALSE) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
