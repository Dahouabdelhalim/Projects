#Karoo Assemblage Zone test

#Each section is repeated to create a node and edge list for each AZ community. 
#At the end are optional steps to plt different route, directed, undirected or weighted networks.
#set your own working directory and file location before starting.
# the working directory and file location for Viglietti is as follows: ("C:/Users/pviglietti/Dropbox/Field Museum/2019 postdoc/Publications/Bipartite Network Analysis/Assemblage zone test")

#********************************
#Network 1: Rubidge
#********************************

AZ_Rubidge95<-read.csv("C:/Users/pviglietti/Dropbox/Field Museum/2019 postdoc/Publications/Bipartite Network Analysis/Assemblage zone test/AZ_Rubidge95.csv")

library(tidyverse)


#Optional step: creating our node objects. Our nodes are stratigraphic interval (bin) and species (taxa).

Bin95 <- distinct(AZ_Rubidge95, Bin_Interval)

Taxa95 <- distinct(AZ_Rubidge95, Taxa)

#Step 2: Combine node objects to create node list

# Need to give both columns in the object were given the same "label" so they can be joined in a node list
# Note the 2nd step above can be skipped as you can name the columns of interest "label" directly from main dataset

Taxa95 <- AZ_Rubidge95 %>%
  distinct(Taxa) %>%
  rename(label = Taxa)

Bin95 <- AZ_Rubidge95 %>%
  distinct(Bin_Interval) %>%
  rename(label = Bin_Interval)

#Step 3 Create Node list
# Now that both columns of interest were given the same "label" they could be joined in a node list

nodes.AZ_Rubidge95 <- full_join(Taxa95,Bin95, by = "label")

#We now have an object called nodes with is our node list for the network analysis.
# Now give give each node a unique id by adding another column that will give it a unique number.
nodes.Rubidge95 <- nodes.AZ_Rubidge95 %>% rowid_to_column("id")

# Step 4 create Edge list
#Now we need a "weight" column which is basically a taxon abundance count for each bin.

per_route.Rubidge95 <- AZ_Rubidge95 %>%  
  group_by(Taxa, Bin_Interval) %>%
  summarise(weight = n()) %>% 
  ungroup()

#Repeat the above steps for making edge list of the other communities by substituting other community dataset.
#The column names are the same for each dataset.
#How to export edge list CSV file
write.csv(per_route.Rubidge95, "edgelist_Rubidge95.csv")


#*********************************
#Network 2: Member
#*********************************

AZ_Member<-read.csv("C:/Users/pviglietti/Dropbox/Field Museum/2019 postdoc/Publications/Bipartite Network Analysis/Assemblage zone test/AZ_Member.csv")

#Step 1: Combine node objects to create node list

Taxa_M <- AZ_Member %>%
  distinct(Taxa) %>%
  rename(label = Taxa)

Bin_M <- AZ_Member %>%
  distinct(Bin_Interval) %>%
  rename(label = Bin_Interval)

#Step 3 Create Node list

nodes.AZ_Member <- full_join(Taxa_M,Bin_M, by = "label")

nodes.Member <- nodes.AZ_Member %>% rowid_to_column("id")

# Step 4 create Edge list

per_route.Member <- AZ_Member %>%  
  group_by(Taxa, Bin_Interval) %>%
  summarise(weight = n()) %>% 
  ungroup()

write.csv(per_route.Member, "edgelist_Member.csv")

#*********************************
#Network 3: Viglietti et al 2021
#*********************************

AZ_Viglietti21<-read.csv("C:/Users/pviglietti/Dropbox/Field Museum/2019 postdoc/Publications/Bipartite Network Analysis/Assemblage zone test/AZ_Viglietti21.csv")

#Step 1: Combine node objects to create node list

Taxa21 <- AZ_Viglietti21 %>%
  distinct(Taxa) %>%
  rename(label = Taxa)

Bin21 <- AZ_Viglietti21 %>%
  distinct(Bin_Interval) %>%
  rename(label = Bin_Interval)

#Step 2 Create Node list

nodes.AZ_Viglietti21 <- full_join(Taxa21,Bin21, by = "label")

nodes.Viglietti21 <- nodes.AZ_Viglietti21 %>% rowid_to_column("id")

# Step 3 create Edge list

per_route.Viglietti21 <- AZ_Viglietti21 %>%  
  group_by(Taxa, Bin_Interval) %>%
  summarise(weight = n()) %>% 
  ungroup()

write.csv(per_route.Viglietti21, "edgelist_Viglietti21.csv")

#****************************
#Network 4: Formation
#****************************

AZ_Formation<-read.csv("C:/Users/pviglietti/Dropbox/Field Museum/2019 postdoc/Publications/Bipartite Network Analysis/Assemblage zone test/AZ_Formation.csv")

#Step 1: Combine node objects to create node list

TaxaFm <- AZ_Formation %>%
  distinct(Taxa) %>%
  rename(label = Taxa)

BinFm <- AZ_Formation %>%
  distinct(Bin_Interval) %>%
  rename(label = Bin_Interval)

#Step 2 Create Node list

nodes.AZ_Formation <- full_join(TaxaFm,BinFm, by = "label")

nodes.Formation <- nodes.AZ_Formation %>% rowid_to_column("id")

# Step 3 create Edge list

per_route.Formation <- AZ_Formation %>%  
  group_by(Taxa, Bin_Interval) %>%
  summarise(weight = n()) %>% 
  ungroup()

write.csv(per_route.Formation, "edgelist_Formation.csv")

#****************************
#Network 5: Broom (1906)
#****************************

AZ_Broom06<-read.csv("C:/Users/pviglietti/Dropbox/Field Museum/2019 postdoc/Publications/Bipartite Network Analysis/Assemblage zone test/AZ_Broom06.csv")

#Step 1: Combine node objects to create node list

TaxaBroom06 <- AZ_Broom06 %>%
  distinct(Taxa) %>%
  rename(label = Taxa)

BinBroom06 <- AZ_Broom06 %>%
  distinct(Bin_Interval) %>%
  rename(label = Bin_Interval)

#Step 2 Create Node list

nodes.AZ_Broom06 <- full_join(TaxaBroom06,BinBroom06, by = "label")

nodes.Broom06 <- nodes.AZ_Broom06 %>% rowid_to_column("id")

# Step 3 create Edge list

per_route.Broom06 <- AZ_Broom06 %>%  
  group_by(Taxa, Bin_Interval) %>%
  summarise(weight = n()) %>% 
  ungroup()

write.csv(per_route.Broom06, "edgelist_Broom06.csv")

#****************************
#Network 6: Gastaldo (2021)
#****************************

AZ_gastaldo<-read.csv("C:/Users/pviglietti/Dropbox/Field Museum/2019 postdoc/Publications/Bipartite Network Analysis/Assemblage zone test/AZ_gastaldo.csv")

#Step 1: Combine node objects to create node list

Taxagastaldo <- AZ_gastaldo %>%
  distinct(Taxa) %>%
  rename(label = Taxa)

Bingastaldo <- AZ_gastaldo %>%
  distinct(Bin_Interval) %>%
  rename(label = Bin_Interval)

#Step 2 Create Node list

nodes.AZ_gastaldo <- full_join(Taxagastaldo,Bingastaldo, by = "label")

nodes.gastaldo <- nodes.AZ_gastaldo %>% rowid_to_column("id")

# Step 3 create Edge list

per_route.gastaldo <- AZ_gastaldo %>%  
  group_by(Taxa, Bin_Interval) %>%
  summarise(weight = n()) %>% 
  ungroup()

write.csv(per_route.gastaldo, "edgelist_gastaldo.csv")

# *****************************************
#Optional steps for making network plots***
#******************************************

#Currently "source" and "destination" columns contain labels rather than ids.
#Now we link ids assigned in the nodes to each location in the source and destination.
#This is accomplished by another join function that links our S ("source") and F ("destination") object to node ids.

#*******************************************
#Network 1: Rubidge
#*******************************************

edges.AZ_Rubidge95 <- per_route.Rubidge95 %>% 
  left_join(nodes.Rubidge95, by = c("Taxa" = "label")) %>% 
  rename(from = id)

edges.Rubidge95 <- edges.AZ_Rubidge95 %>% 
  left_join(nodes.Rubidge95, by = c("Bin_Interval" = "label")) %>% 
  rename(to = id)

#when we did the previous analysis the per route object data frame was on the left of the columns
#Now we redorder so that the "to" and "from" lists are on the left of the data frame
#our "source" and "destination" columns are also removed as the "to" and "from" serve this purpose (numerically).
edges.Rubidge95_final <- select(edges.Rubidge95, from, to, weight)
edges.Rubidge95_final


#Step 6
#Repeat steps for different bin intervals

#Step 7 plotting routes network

#Making a network object

library(network)

routes_Rubidge95 <- network(edges.Rubidge95, vertex.attr = nodes.Rubidge95, matrix.type = "edgelist", ignore.eval = FALSE)
#can view the type of class the routes network object is by using the class() function.

class(routes_Rubidge95)
summary(routes_Rubidge95)

#We can now plot a rudimentary network graph
plot(routes_Rubidge95, vertex.cex = 3)

#Graph cleanup and igraph.Ths step removes routes_Rubidge95 object so don't do step if you want to do more plots.

detach(package:network)
rm(routes_Rubidge95)

#Install package igraph

library(igraph)

routes_Rubidge95 <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

plot(routes_Rubidge95, edge.arrow.size = 0.2)

plot(routes_Rubidge95, layout = layout_with_graphopt, edge.arrow.size = 0.2)

# Load Tidygraph and ggraph
#Always start by loading necessary packages.
library(tidygraph)
library(ggraph)
# Going to create a tbl_graph using Tidygraph which uses an edge and node tibble.
# tbl_graphs are essentially igraph objects.
routes_tidy <- tbl_graph(nodes = nodes.Rubidge95, edges = edges.Rubidge95, directed = TRUE)

routes_tidy

#Plot new graph with ggraph

ggraph(routes_tidy) + geom_edge_link() + geom_node_point() + theme_graph()
# We can also show the weight of the edges and make graph more informative.

ggraph(routes_tidy, layout = "graphopt") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Letters") +
  theme_graph()

#Arc graphs
#Indicate directionality of the edges

ggraph(routes_igraph, layout = "linear") + 
  geom_edge_arc(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "Letters") +
  theme_graph()

#Interactive network graphs using visNetwork and networkD3

library(visNetwork)
library(networkD3)

#visNetwork function uses a nodes list and edges list to create an interactive graph. 
#The nodes list must include an "id" column, and the edge list must have "from" and "to" columns. 

visNetwork(nodes.Rubidge95, edges.Rubidge95)

#The mutate() function allows us to create a graph with variable edge widths.

edges <- mutate(edges.Rubidge95, width = weight/5 + 1)

visNetwork(nodes.Rubidge95, edges) %>% 
  visIgraphLayout(layout = "layout_with_fr") %>% 
  visEdges(arrows = "middle")

#NetworkD3

#Some preparation required to present data in networkD3 graph.
#The edge and node list requires that the IDs be a series of numeric integers that begin with 0. 
#Currently, the node IDs for our data begin with 1, and so we have to do a bit of data manipulation.
#Once again, this can be done with the mutate() function. 
#The goal is to recreate the current columns, while subtracting 1 from each ID. 
#making new node and edge objects that are manipulated to start ids at 0 rather than 1.

nodes_d3 <- mutate(nodes.Rubidge95, id = id - 1)
edges_d3 <- mutate(edges, from = from - 1, to = to - 1)

#It is now possible to plot a networkD3 graph.

forceNetwork(Links=edges_d3, Nodes=nodes_d3, Source="from", Target="to", 
             NodeID="label", Group="id", Value="weight", 
             opacity= 1, fontSize= 16, zoom= TRUE)

#Making a Sankey diagram
#Good for when there are not too many nodes...
#uses the sankeyNetwork function which takes the same arguments as forceNetwork.
sankeyNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
              NodeID = "label", Value = "weight", fontSize = 16, unit = "Letter(s)")



#*********************************
#Network 2: Member
#*********************************

edges.AZ_Member <- per_route.Member %>% 
  left_join(nodes.Member, by = c("Taxa" = "label")) %>% 
  rename(from = id)

edges.Member <- edges.AZ_Member %>% 
  left_join(nodes.Member, by = c("Bin_Interval" = "label")) %>% 
  rename(to = id)

#when we did the previous analysis the per route object data frame was on the left of the columns
#Now we redorder so that the "to" and "from" lists are on the left of the data frame
#our "source" and "destination" columns are also removed as the "to" and "from" serve this purpose (numerically).
edges.Member_final <- select(edges.Member, from, to, weight)
edges.Member_final


#Step 6
#Repeat steps for different bin intervals

#Step 7 plotting routes network

#Making a network object

library(network)

routes_Member <- network(edges.Member, vertex.attr = nodes.Member, matrix.type = "edgelist", ignore.eval = FALSE)
#can view the type of class the routes network object is by using the class() function.

class(routes_Member)
summary(routes_Member)

#We can now plot a rudimentary network graph
plot(routes_Member, vertex.cex = 3)

#Graph cleanup and igraph.Ths step removes routes_Rubidge95 object so don't do step if you want to do more plots.

detach(package:network)
rm(routes_Member)

#Install package igraph

library(igraph)

routes_Member <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

plot(routes_Member, edge.arrow.size = 0.2)

plot(routes_Member, layout = layout_with_graphopt, edge.arrow.size = 0.2)

# Load Tidygraph and ggraph
#Always start by loading necessary packages.
library(tidygraph)
library(ggraph)
# Going to create a tbl_graph using Tidygraph which uses an edge and node tibble.
# tbl_graphs are essentially igraph objects.
routes_tidy <- tbl_graph(nodes = nodes.Member, edges = edges.Member, directed = TRUE)

routes_tidy

#Plot new graph with ggraph

ggraph(routes_tidy) + geom_edge_link() + geom_node_point() + theme_graph()
# We can also show the weight of the edges and make graph more informative.

ggraph(routes_tidy, layout = "graphopt") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Letters") +
  theme_graph()

#Arc graphs
#Indicate directionality of the edges

ggraph(routes_igraph, layout = "linear") + 
  geom_edge_arc(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "Letters") +
  theme_graph()

#Interactive network graphs using visNetwork and networkD3

library(visNetwork)
library(networkD3)

#visNetwork function uses a nodes list and edges list to create an interactive graph. 
#The nodes list must include an "id" column, and the edge list must have "from" and "to" columns. 

visNetwork(nodes.Member, edges.Member)

#The mutate() function allows us to create a graph with variable edge widths.

edges <- mutate(edges.Member, width = weight/5 + 1)

visNetwork(nodes.Member, edges) %>% 
  visIgraphLayout(layout = "layout_with_fr") %>% 
  visEdges(arrows = "middle")

#NetworkD3

#Some preparation required to present data in networkD3 graph.
#The edge and node list requires that the IDs be a series of numeric integers that begin with 0. 
#Currently, the node IDs for our data begin with 1, and so we have to do a bit of data manipulation.
#Once again, this can be done with the mutate() function. 
#The goal is to recreate the current columns, while subtracting 1 from each ID. 
#making new node and edge objects that are manipulated to start ids at 0 rather than 1.

nodes_d3 <- mutate(nodes.Member, id = id - 1)
edges_d3 <- mutate(edges, from = from - 1, to = to - 1)

#It is now possible to plot a networkD3 graph.

forceNetwork(Links=edges_d3, Nodes=nodes_d3, Source="from", Target="to", 
             NodeID="label", Group="id", Value="weight", 
             opacity= 1, fontSize= 16, zoom= TRUE)

#Making a Sankey diagram
#Good for when there are not too many nodes...
#uses the sankeyNetwork function which takes the same arguments as forceNetwork.
sankeyNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
              NodeID = "label", Value = "weight", fontSize = 16, unit = "Letter(s)")


#*********************************
#Network 3: Viglietti
#*********************************

edges.AZ_Viglietti21 <- per_route.Viglietti21 %>% 
  left_join(nodes.Viglietti20, by = c("Taxa" = "label")) %>% 
  rename(from = id)

edges.Viglietti21 <- edges.AZ_Viglietti20 %>% 
  left_join(nodes.Viglietti20, by = c("Bin_Interval" = "label")) %>% 
  rename(to = id)

#when we did the previous analysis the per route object data frame was on the left of the columns
#Now we redorder so that the "to" and "from" lists are on the left of the data frame
#our "source" and "destination" columns are also removed as the "to" and "from" serve this purpose (numerically).
edges.Viglietti21_final <- select(edges.Viglietti21, from, to, weight)
edges.Viglietti21_final


#Step 6
#Repeat steps for different bin intervals

#Step 7 plotting routes network

#Making a network object

library(network)

routes_Viglietti21 <- network(edges.Viglietti21, vertex.attr = nodes.Viglietti21, matrix.type = "edgelist", ignore.eval = FALSE)
#can view the type of class the routes network object is by using the class() function.

class(routes_Viglietti21)
summary(routes_Viglietti21)

#We can now plot a rudimentary network graph
plot(routes_Viglietti21, vertex.cex = 3)

#Graph cleanup and igraph.Ths step removes routes_Rubidge95 object so don't do step if you want to do more plots.

detach(package:network)
rm(routes_Viglietti21)

#Install package igraph

library(igraph)

routes_Viglietti21 <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

plot(routes_Viglietti21, edge.arrow.size = 0.2)

plot(routes_Viglietti21, layout = layout_with_graphopt, edge.arrow.size = 0.2)

# Load Tidygraph and ggraph
#Always start by loading necessary packages.
library(tidygraph)
library(ggraph)
# Going to create a tbl_graph using Tidygraph which uses an edge and node tibble.
# tbl_graphs are essentially igraph objects.
routes_tidy <- tbl_graph(nodes = nodes.Viglietti21, edges = edges.Viglietti21, directed = TRUE)

routes_tidy

#Plot new graph with ggraph

ggraph(routes_tidy) + geom_edge_link() + geom_node_point() + theme_graph()
# We can also show the weight of the edges and make graph more informative.

ggraph(routes_tidy, layout = "graphopt") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Letters") +
  theme_graph()

#Arc graphs
#Indicate directionality of the edges

ggraph(routes_igraph, layout = "linear") + 
  geom_edge_arc(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "Letters") +
  theme_graph()

#Interactive network graphs using visNetwork and networkD3

library(visNetwork)
library(networkD3)

#visNetwork function uses a nodes list and edges list to create an interactive graph. 
#The nodes list must include an "id" column, and the edge list must have "from" and "to" columns. 

visNetwork(nodes.Viglietti21, edges.Viglietti21)

#The mutate() function allows us to create a graph with variable edge widths.

edges <- mutate(edges.Viglietti21, width = weight/5 + 1)

visNetwork(nodes.Viglietti20, edges) %>% 
  visIgraphLayout(layout = "layout_with_fr") %>% 
  visEdges(arrows = "middle")

#NetworkD3

#Some preparation required to present data in networkD3 graph.
#The edge and node list requires that the IDs be a series of numeric integers that begin with 0. 
#Currently, the node IDs for our data begin with 1, and so we have to do a bit of data manipulation.
#Once again, this can be done with the mutate() function. 
#The goal is to recreate the current columns, while subtracting 1 from each ID. 
#making new node and edge objects that are manipulated to start ids at 0 rather than 1.

nodes_d3 <- mutate(nodes.Viglietti21, id = id - 1)
edges_d3 <- mutate(edges, from = from - 1, to = to - 1)

#It is now possible to plot a networkD3 graph.

forceNetwork(Links=edges_d3, Nodes=nodes_d3, Source="from", Target="to", 
             NodeID="label", Group="id", Value="weight", 
             opacity= 1, fontSize= 16, zoom= TRUE)

#Making a Sankey diagram
#Good for when there are not too many nodes...
#uses the sankeyNetwork function which takes the same arguments as forceNetwork.
sankeyNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
              NodeID = "label", Value = "weight", fontSize = 16, unit = "Letter(s)")

#*********************************
#Network 4: Formation
#*********************************

edges.AZ_Formation <- per_route.Formation %>% 
  left_join(nodes.Formation, by = c("Taxa" = "label")) %>% 
  rename(from = id)

edges.Formation <- edges.AZ_Formation %>% 
  left_join(nodes.Formation, by = c("Bin_Interval" = "label")) %>% 
  rename(to = id)

#when we did the previous analysis the per route object data frame was on the left of the columns
#Now we redorder so that the "to" and "from" lists are on the left of the data frame
#our "source" and "destination" columns are also removed as the "to" and "from" serve this purpose (numerically).
edges.Formation_final <- select(edges.Formation, from, to, weight)
edges.Formation_final


#Step 6
#Repeat steps for different bin intervals

#Step 7 plotting routes network

#Making a network object

library(network)

routes_Formation <- network(edges.Formation, vertex.attr = nodes.Formation, matrix.type = "edgelist", ignore.eval = FALSE)
#can view the type of class the routes network object is by using the class() function.

class(routes_Formation)
summary(routes_Formation)

#We can now plot a rudimentary network graph
plot(routes_Formation, vertex.cex = 3)

#Graph cleanup and igraph.Ths step removes routes_Formation object so don't do step if you want to do more plots.

detach(package:network)
rm(routes_Formation)

#Install package igraph

library(igraph)

routes_Formation <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

plot(routes_Formation, edge.arrow.size = 0.2)

plot(routes_Formation, layout = layout_with_graphopt, edge.arrow.size = 0.2)

# Load Tidygraph and ggraph
#Always start by loading necessary packages.
library(tidygraph)
library(ggraph)
# Going to create a tbl_graph using Tidygraph which uses an edge and node tibble.
# tbl_graphs are essentially igraph objects.
routes_tidy <- tbl_graph(nodes = nodes.Formation, edges = edges.Viglietti20, directed = TRUE)

routes_tidy

#Plot new graph with ggraph

ggraph(routes_tidy) + geom_edge_link() + geom_node_point() + theme_graph()
# We can also show the weight of the edges and make graph more informative.

ggraph(routes_tidy, layout = "graphopt") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Letters") +
  theme_graph()

#Arc graphs
#Indicate directionality of the edges

ggraph(routes_igraph, layout = "linear") + 
  geom_edge_arc(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "Letters") +
  theme_graph()

#Interactive network graphs using visNetwork and networkD3

library(visNetwork)
library(networkD3)

#visNetwork function uses a nodes list and edges list to create an interactive graph. 
#The nodes list must include an "id" column, and the edge list must have "from" and "to" columns. 

visNetwork(nodes.Formation, edges.Formation)

#The mutate() function allows us to create a graph with variable edge widths.

edges <- mutate(edges.Formation, width = weight/5 + 1)

visNetwork(nodes.Formation, edges) %>% 
  visIgraphLayout(layout = "layout_with_fr") %>% 
  visEdges(arrows = "middle")

#NetworkD3

#Some preparation required to present data in networkD3 graph.
#The edge and node list requires that the IDs be a series of numeric integers that begin with 0. 
#Currently, the node IDs for our data begin with 1, and so we have to do a bit of data manipulation.
#Once again, this can be done with the mutate() function. 
#The goal is to recreate the current columns, while subtracting 1 from each ID. 
#making new node and edge objects that are manipulated to start ids at 0 rather than 1.

nodes_d3 <- mutate(nodes.Formation, id = id - 1)
edges_d3 <- mutate(edges, from = from - 1, to = to - 1)

#It is now possible to plot a networkD3 graph.

forceNetwork(Links=edges_d3, Nodes=nodes_d3, Source="from", Target="to", 
             NodeID="label", Group="id", Value="weight", 
             opacity= 1, fontSize= 16, zoom= TRUE)

#Making a Sankey diagram
#Good for when there are not too many nodes...
#uses the sankeyNetwork function which takes the same arguments as forceNetwork.
sankeyNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
              NodeID = "label", Value = "weight", fontSize = 16, unit = "Letter(s)")

