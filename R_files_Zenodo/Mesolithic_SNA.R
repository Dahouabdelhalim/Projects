##Script used in the article “Reconstructing social networks on the Iberian Peninsula using ornaments”.


## INSTALL THE PACKAGE

install.packages("igraph") 

# LOAD IGRAPH
library(igraph)



# CONSTRUCT THE NETWORK

# Convert the similarity matrix unto an Igraph object
my_data <-as.matrix(read.csv(file.choose(), header=T, row.names=1)) # The command 'file.choose()' allows to select an input file. In this case, one should select either 'Early_Meso.csv' or 'Late_Meso.csv'. 

my_data
g <- graph_from_adjacency_matrix (my_data, mode = "undirected", weighted=TRUE, diag=FALSE) # create an 'igraph object'
g 
V(g)
E(g)
E(g)$weight 


#PLOT THE NETWORK 
# The netowk is plotted according to a driving-force algorithm that defines the appearance of nodes and links (colour, labels and size)

plot.igraph(g,vertex.size=10,vertex.label.color="black", arrow.size= 0.05, vertex.label.font=1,vertex.color="lightblue",edge.color="black", main="Early Mesolithic")


#GLOBAL NETWORK METRICS 

# Number of nodes

vcount (g)

# Number of edges

ecount(g)

# Density 

Density3 <- ecount(g)/(vcount(g)*(vcount(g)-1)/2)
Density3

# Average degree 

mean(degree(g))

#Average weighted degree

mean(E(g)$weight)

#Average path lenght 
mean_distance(g)

#Assortativity

#It needs a table (Geo_unit)that includes the ID of the assemblages and the geographical unit they belong to (GU)
Geo_unit <-read.csv('Geo_Units.csv', header=T, row.names=1) 
assortativity_nominal (g, as.integer(Geo_unit$GU),directed=TRUE)


# NODE METRICS

# Degree  

degree<-degree(g, mode="total", loops= FALSE)
degree

#Betweenness

Betweenness<-betweenness(g, directed = FALSE, weights = NULL)
Betweenness

#Weighted degree 

All_Strength<-strength(g, mode="all") ## weighted by default
All_Strength

