#install.packages('igraph')
#install.packages('ggplot2')
#install.packages('reshape2')
#install.packages('ggplot')
#install.packages('qgraph')
library('igraph')
library('qgraph')
# create network visualisation with igraph
#install.packages("data.table")


# VIS COORDINATES INFO
# https://rstudio-pubs-static.s3.amazonaws.com/341807_7b9d1a4e787146d492115f90ad75cd2d.html


#install.packages("RColorBrewer")
library(RColorBrewer)

# EXPERIMENT WITH COLORS

# colors
# cols <- colorRampPalette(brewer.pal(8,"Dark2"))(8)
#pal2 <- rainbow(8, alpha=.5)  
#plot(x=1:10, y=1:10, pch=19, cex=5, col=pal2)





# 1. read in nodes and edges data

# GET EDGES

#edges <- read.csv(file="P:/thesis/thesis_final_data/produced_data/heterogeneous_networks/machine_learning/cos_sim_0_1/edges_np.csv", header=TRUE, sep=",")

#edges <- read.csv(file="P:/thesis/thesis_final_data/produced_data/heterogeneous_networks/machine_learning/cos_sim_0_0/machine_learning_edges.csv", header=TRUE, sep=",")
edges <- read.csv(file="P:/thesis/thesis_final_data/produced_data/heterogeneous_networks/social_network_analysis/cos_sim_0_0/social_network_analysis_edges.csv", header=TRUE, sep=",")
str(edges)

# optional: filter edges
#filt_cos_value = 0.01
#edges_filt <- edges[edges$Weight > filt_cos_value, ] 
#str(edges_filt)
#str(edges)



# GET NODES

#nodes <- read.csv(file="P:/thesis/thesis_final_data/produced_data/heterogeneous_networks/machine_learning/cos_sim_0_1/nodes_np.csv", header=TRUE, sep=",")

# machine learning
#nodes <- read.csv(file="P:/thesis/thesis_final_data/produced_data/heterogeneous_networks/machine_learning/cos_sim_0_0/machine_learning_nodes.csv", header=TRUE, sep=",")

# social network analysis
nodes <- read.csv(file="P:/thesis/thesis_final_data/produced_data/heterogeneous_networks/social_network_analysis/cos_sim_0_0/social_network_analysis_nodes.csv", header=TRUE, sep=",")

str(nodes)

nodes$actant <- as.character(nodes$actant)
nodes$diff_binary <- as.character(nodes$diff_binary)
nodes$actant_type <- as.character(nodes$actant_type)
nodes$shape <- as.character(nodes$shape)

str(nodes)

nodes_id <- subset(nodes, select = c(Id))
str(nodes_id)



# CREATE IGRAPH OBJECT

# filtered edges - I'm not sure if this works!! I think edges would
# need to be filtered in Python and given appropritate ids and node ids

#net <- graph_from_data_frame(d=edges_filt, vertices=nodes_id, directed=F) 
#edge_attr(net, 'weight', index = E(net)) <- edges_filt$weight

net <- graph_from_data_frame(d=edges, vertices=nodes_id, directed=F) 
edge_attr(net, 'weight', index = E(net)) <- edges$Weight

# add edge attributes
#edge_attr(net, 'weight', index = E(net)) <- edges$Weight
#edge_attr(net, 'weight', index = E(net)) <- edges_filt$weight
is_weighted(net)


# add node attributes
vertex_attr(net, 'actant', index = V(net)) <- nodes$actant
vertex_attr(net, 'type', index = V(net)) <- nodes$actant_type
vertex_attr(net, 'diff', index = V(net)) <- nodes$diff_binary
vertex_attr(net, 'shape', index = V(net)) <- nodes$shape
vertex_attr(net, 'betweenness', index = V(net)) <- betweenness(net)
vertex_attr(net, 'betweenness_norm', index = V(net)) <- betweenness(net, normalized = TRUE)

# add color to nodes: MANUAL CHOICE
# red and yellow
pal1 <- heat.colors(2, alpha=.5) 
pal1
vertex_attr(net, 'color', index = V(net)) <- pal1[as.numeric(as.factor(vertex_attr(net, "diff")))]

# red and blue
pal2 <- rainbow(2, alpha=.5)  
vertex_attr(net, 'color', index = V(net)) <- pal2[as.numeric(as.factor(vertex_attr(net, "diff")))]


#plot(net, vertex.color = pal2[as.numeric(as.factor(vertex_attr(net, "diff")))])

# experiment with layouts
#my_layout = layout_with_kk
#my_layout = layout_with_mds
#my_layout = layout_with_lgl
#my_layout = layout_with_dri
#my_layout = layout_on_sphere
#my_layout = layout.lgl
#my_layout = layout.reingold.tilford(g, circular=T)
#my_layout = layout_components
#my_layout = layout_as_tree
#my_layout = layout_with_dh

my_layout <- layout_with_fr(net, dim = 3, niter = 5) # Fruchterman-Reingold layout algorithm, good but need to come up with a way for nodes not to overlap

# PLOT WITH DIFFEREN LABELS

# plot with vertex label
plot(net, edge.width = E(net)$weight, vertex.shape=V(net)$shape, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)

# plot without vertex labesl, same node size
plot(net, edge.width = E(net)$weight, vertex.size=3, vertex.shape=V(net)$shape, vertex.label=NA, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)

# plot with vertex labesl, same node size
plot(net, edge.width = E(net)$weight, vertex.size=3, vertex.shape=V(net)$shape, vertex.label=V(net)$diff, vertex.label.cex=0.25, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)

# plot with vertex labesl, same node size
plot(net, edge.width = E(net)$weight, vertex.size=3, vertex.shape=V(net)$shape, vertex.label=V(net)$actant, vertex.label.cex=0.25, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)



# plot without vertex label and suitable vertex size
vertex_attr(net, 'size', index = V(net)) <- degree(net, normalized = TRUE)
vertex_attr(net, 'size_2', index = V(net)) <- V(net)$size * 30
plot(net, edge.width = E(net)$weight, vertex.size=6, vertex.shape=V(g)$shape, vertex.label=NA, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)
#plot(net, edge.width = E(net)$weight, vertex.size=V(net)$size_2, vertex.shape=V(net)$shape, vertex.label=NA, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)

# plot with vertex label and suitable vertex size
plot(net, edge.width = E(net)$weight, vertex.size=V(net)$size_2, vertex.shape=V(net)$shape, vertex.label=V(net)$diff, vertex.label.cex = 0.25, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)


# this created the funny.pdf
# plot without vertex label
#plot(net, edge.width = E(net)$weight, vertex.size=betweenness(net), vertex.shape=V(g)$shape, vertex.label=NA, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)

# plot without vertex label, normalised betweenness centrality
vertex_attr(net, 'size', index = V(net)) <- betweenness(net, normalized = TRUE)
vertex_attr(net, 'size_2', index = V(net)) <- V(net)$size * 30
plot(net, edge.width = E(net)$weight, vertex.size= V(net)$size_2, vertex.shape=V(net)$shape, vertex.label=NA, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)



# plot without vertex labesl, same node size
plot(net, edge.width = E(net)$weight, vertex.size=3, vertex.shape=V(net)$shape, vertex.label=NA, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)

my_layout <- layout_with_fr(net, dim = 3, niter = 1) # Fruchterman-Reingold layout algorithm, good but need to come up with a way for nodes not to overlap
my_layout <- layout_with_fr(net, dim = 3)
plot(net, edge.width = E(net)$weight, vertex.size=4, vertex.shape=V(net)$shape, vertex.label=V(net)$Id, vertex.label.cex=0.3, vertex.frame.color="gray", vertex.label.color=c("black"), layout = my_layout)




# add legend for shape

legend_cat_3 <- data.frame(type_n = c('author', 'keyword', 'new' , 'not new'), 
                           shape_n = c(21, 22, 23, 23),
                           color = c('#696969', '#696969', "#FF000080", "#FFFF0080"))

legend_cat_3$type_n <- as.character(legend_cat_3$type_n )
legend_cat_3$shape_n <- as.numeric(legend_cat_3$shape_n )
legend_cat_3$color <- as.character(legend_cat_3$color )

str(legend_cat_3)

legend(x = "bottomleft",      ## position, also takes x,y coordinates
       legend = legend_cat_3$type_n,
       pch = legend_cat_3$shape_n, ## legend symbols see ?points,
       col = legend_cat_3$color,
       cex=0.7,
       bty = "n",
       title = "Node type")

legend(1, 95,       ## position, also takes x,y coordinates
       legend = legend_cat_3$type_n,
       pch = legend_cat_3$shape_n, ## legend symbols see ?points,
       col = legend_cat_3$color,
       cex=0.7,
       bty = "n",
       title = "Node type")


# save vertices as csv

vertex_attr(net, 'name', index = V(net)) <- nodes$word_or_AU
vertex_attr(net, 'id_num', index = V(net)) <- nodes$Id

vertices_df <- as_data_frame(net, what="vertices")

str(vertices_df)

vertices_df$diff[vertices_df$diff == 'previously_used'] <- 'old'
vertices_df$diff[vertices_df$diff == 'first_time_with_keyword'] <- 'new'

str(vertices_df)

vertices_df_2 <-subset(vertices_df, select = -c(shape))
vertices_df_2 <-subset(vertices_df_2, select = -c(diff))
vertices_df_2 <-subset(vertices_df_2, select = -c(betweenness))
str(vertices_df_2)

# reorder columns
col_order <- c("id_num", "name", "type", "betweenness_norm")
vertices_df_3 <- vertices_df_2[, col_order]
str(vertices_df_3)

write.csv(vertices_df_3,"P:/thesis_final_data/produced_data/heterogeneous_networks/machine_learning/cos_sim_0_01_nodes_from_r.csv", row.names = FALSE)










# add legend for novelty

legend_cats <- data.frame(attr = unique(V(net)$diff),
                          color = unique(V(net)$color))

legend_cats$attr <- as.character(legend_cats$attr)
legend_cats$color <- as.character(legend_cats$color)

str(legend_cats)

legend_cats <- legend_cats[order(legend_cats$attr), c(1, 2)]

legend(x = "bottomleft",      ## position, also takes x,y coordinates
       legend = legend_cats$attr,
       pch = 23,              ## legend symbols see ?points
       col = legend_cats$color,
       cex=0.7,
       bty = "n",
       title = "Node novelty")


# add legend for shape

legend_cat_2 <- data.frame(type_n = c('author', 'keyword'), shape_n = c(21, 22))

legend_cat_2$type_n <- as.character(legend_cat_2$type_n )
str(legend_cat_2)

legend(x = "topleft",      ## position, also takes x,y coordinates
       legend = legend_cat_2$type_n,
       pch = legend_cat_2$shape_n,              ## legend symbols see ?points,
       cex=0.7,
       bty = "n",
       title = "Node type")



