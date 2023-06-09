#load packages
require(gstudio)
require(popgraph)
require(ggmap)
require(vegan)
require(igraph)
require(ggplot2)
require(fields)
require(GGally)
library("adegenet")
library("hierfstat")
library("pegas")

#load mtDNA data
spiderAll1 <- read_population("AllSites1Phxb.csv",type="snp",header=TRUE,sep=",",phased=FALSE,locus.columns=seq(6,131))

#make a popgraph
data <- to_mv(spiderAll1)
pops <- spiderAll1$Pop
cluster <- spiderAll1$Cluster
spgraph <- popgraph(x=data, groups=pops)
plot(spgraph)

#make popgraph look better
layout <- layout.fruchterman.reingold( spgraph )
plot( spgraph, layout=layout)

t <- table(spiderAll1$Cluster, spiderAll1$Pop)
nodes <- V(spgraph)$name
color <- rep("#f4a582", length(nodes))
# the A clade
color[t[1, ] != 0] <- "yellow"
# the B clade
color[t[2, ] != 0] <- "blue"
V(spgraph)$color <- color
plot(spgraph, vertex.label = "",layout=layout)

#add data to graph
spgraph<-decorate_graph(spgraph,spiderAll1, stratum="Pop")

#get geographic distances
pDist <- rdist.earth( cbind( V(spgraph)$Longitude, V(spgraph)$Latitude ) )
#get cGD as genetic distance
cGD <- to_matrix( spgraph, mode="shortest path")

df <- data.frame( cGD=cGD[upper.tri(cGD)], Phys=pDist[upper.tri(pDist)])
#test for IBD
cor.test( df$Phys, df$cGD, method="spearman")
qplot( Phys, cGD, geom="point", data=df) + stat_smooth(method="lm") + xlab("Physical Distance") + ylab("Conditional Genetic Distance")

#get social network params
df.nodes <- data.frame(Pop=V(spgraph)$name)
df.nodes$closeness <- closeness(spgraph)
df.nodes$betweenness <- betweenness(spgraph)
df.nodes$degree <- degree( spgraph)
df.nodes$eigenCent <- evcent( spgraph )$vector
df.nodes$Type<-factor(V(spgraph)$color)

#basic summary of SN params
summary(df.nodes, color="Type")
ggpairs(df.nodes, columns=2:6,color="Type")

#write these to csv, can play with them in excel, user friendly
write.csv(df.nodes, file="Mtnodespopgraph.csv")

#run a PCA
x.2 <- to_mv(spiderAll1, drop.allele = TRUE)
fit.pca.2 <- princomp(x.2)
pred.2 <- predict(fit.pca.2)
df.2 <- data.frame(PC1 = pred.2[, 1], PC2 = pred.2[, 2],  PC3 = pred.2[, 3], PC4 = pred.2[, 4], 
                   Type = spiderAll1$Cluster, Pop = spiderAll1$Pop)
ggplot(df.2) + geom_point(aes(x = PC1, y = PC2, color = Pop),  size = 3, alpha = 0.75)
#get PC eigenvalues
summary(fit.pca.2)
#color code PCA by type
ggplot(df.2) + geom_point(aes(x = PC3, y = PC4, color = Type), 
                          size = 3, alpha = 0.75)

mtdiv<-genetic_diversity(spiderAll1, stratum = "Pop", mode=c("He","Ho","Fis"))
###########################################################
###load nuclear data
spnu <- read_population("A1strD.csv",type="separated",header=T,sep=",",locus.columns=seq(6,5006))

#making a popgraph
data1 <- to_mv(spnu)
pop1 <- spnu$Pop
clust1 <- spnu$Type
graph1 <- popgraph(x=data1, groups=pop1)
plot(graph1)

layout <- layout.fruchterman.reingold( graph1 )
plot( graph1, layout=layout)

t1 <- table(spnu$Type, spnu$Pop)
nodes <- V(graph1)$name
color <- rep("#f4a582", length(nodes))
# the A clade
color[t1[1, ] != 0] <- "yellow"
# the B clade
color[t1[2, ] != 0] <- "blue"
V(graph1)$color <- color
plot(graph1,layout=layout)

#PCA
x.2b <- to_mv(spnu, drop.allele = TRUE)
fit.pca.2b <- prcomp(x.2b)
pred.2b <- predict(fit.pca.2b)
df.2b <- data.frame(PC1 = pred.2b[, 1], PC2 = pred.2b[, 2], PC3 = pred.2b[, 3], 
                    PC4 = pred.2b[, 4],PC5 = pred.2b[, 5],PC6 = pred.2b[, 6],
                    PC7 = pred.2b[, 7],PC8 = pred.2b[, 8],PC9 = pred.2b[, 9],PC10 = pred.2b[, 10],
                    Type = spnu$Type, Pop = spnu$Pop)
ggplot(df.2b) + geom_point(aes(x = PC1, y = PC2, color = Type ), 
                           size = 3, alpha = 0.75)

ggplot(df.2b) + geom_point(aes(x = PC1, y = PC2, color = Pop, shape=Type ), 
                           size = 3, alpha = 0.75)
summary(fit.pca.2b)
#PC % variance
#PC1=16.06, PC2=10.83, PC3= 4.94, PC4=2.21, PC5=1.41,
#PC6=1.494, PC7=1.378, PC8=1.037, PC9=1.032, PC10=.998


ggplot(df.2b) + geom_point(aes(x = PC3, y = PC4, color = Type ), 
                           size = 3, alpha = 0.75)

ggplot(df.2b) + geom_point(aes(x = PC3, y = PC4, color= Pop, shape= Type ), 
                           size = 3, alpha = 0.75)
ggplot(df.2b) + geom_point(aes(x = PC5, y = PC6, color = Pop, shape=Type ), 
                           size = 3, alpha = 0.75)


graph1<-decorate_graph(graph1,spnu, stratum="Pop")

pDist1 <- rdist.earth( cbind( V(graph1)$Longitude, V(graph1)$Latitude ) )
cGD1 <- to_matrix( graph1, mode="shortest path")

df1 <- data.frame( cGD=cGD1[upper.tri(cGD1)], Phys=pDist1[upper.tri(pDist1)])
cor.test( df1$Phys, df1$cGD, method="spearman")
#Spearman's rank correlation rho
#data:  df1$Phys and df1$cGD
#S = 1392000, p-value = 0.1563
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#rho 
#0.09814217
#So there is no significant IBD

qplot( Phys, cGD, geom="point", data=df1) + stat_smooth(method="loess") + xlab("Physical Distance") + ylab("Conditional Genetic Distance")



cong <- congruence_topology(spgraph,graph1)
plot(cong)
test_congruence(spgraph,graph1,method="combinatorial")

#.0097


###########################################
#popgraph params
df.edge <- data.frame( Weight=E(spgraph)$weight )
df.edge$betweenness <- edge.betweenness(spgraph)

ggpairs(df.edge)

df.nodes <- data.frame(Pop=V(spgraph)$name)
df.nodes$closeness <- closeness(spgraph)
df.nodes$betweenness <- betweenness(spgraph)
df.nodes$degree <- degree( spgraph)
df.nodes$eigenCent <- evcent( spgraph )$vector
df.nodes$Type<-factor(V(spgraph)$color)
summary(df.nodes, color="Type")
ggpairs(df.nodes, columns=2:6,color="Type")


df.edge1 <- data.frame( Weight=E(graph1)$weight )
df.edge1$betweenness <- edge.betweenness(graph1)
ggpairs(df.edge1)

df.nodes1 <- data.frame(Pop=V(graph1)$name)
df.nodes1$closeness <- closeness(graph1)
df.nodes1$betweenness <- betweenness(graph1)
df.nodes1$degree <- degree( graph1)
df.nodes1$eigenCent <- evcent( graph1 )$vector
df.nodes1$Type<-factor(V(graph1)$color)
summary(df.nodes1, color="Type")
ggpairs(df.nodes1, columns=2:6,color="Type")



adj1 <-to_matrix(spgraph, mode="adjacency")
EW1 <-to_matrix(spgraph, mode="edge weight")
cgd1 <-to_matrix(spgraph, mode="shortest path")

adj2 <-to_matrix(graph1, mode="adjacency")
EW2 <-to_matrix(graph1, mode="edge weight")
cgd2 <-to_matrix(graph1, mode="shortest path")

dfcgd <- data.frame( cGDmt=cgd1[upper.tri(cgd1)], cGDnu=cgd2[upper.tri(cgd2)])
cor.test( dfcgd$cGDnu, dfcgd$cGDmt, method="spearman")
#Spearman's rank correlation rho
#data:  dfcgd$cGDnu and dfcgd$cGDmt
#S = 1187200, p-value = 0.0007493
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.2308417 
qplot( cGDnu, cGDmt, geom="point", data=dfcgd) + stat_smooth(method="loess") + xlab("cGD nuclear") + ylab("cGD mt")



dfadj <- data.frame( adjmt=adj1[upper.tri(adj1)], adjnu=adj2[upper.tri(adj2)])
cor.test( dfadj$adjnu, dfadj$adjmt, method="spearman")
#Spearman's rank correlation rho
#data:  dfadj$adjnu and dfadj$adjmt
#S = 1270700, p-value = 0.0103
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.1766948 

dfEW <- data.frame( EWmt=EW1[upper.tri(EW1)], EWnu=EW2[upper.tri(EW2)])
cor.test( dfEW$EWnu, dfEW$EWmt, method="spearman")
#Spearman's rank correlation rho
#data:  dfEW$EWnu and dfEW$EWmt
#S = 1321400, p-value = 0.03719
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#rho 
#0.1438936 

write.csv(df.edge, file="Mtedgepopgraph.csv")
write.csv(df.edge1, file="Nuedgepopgraph.csv")
write.csv(df.nodes, file="Mtnodespopgraph.csv")
write.csv(df.nodes1, file="Nunodespopgraph.csv")
write.csv(dfcgd, file="cGDmatrix.csv")
write.csv(cgd1, file="Mtcgd.csv")
write.csv(cgd2, file="Nucgd.csv")

dfnodes2 <-merge(df.nodes,df.nodes1, by="Pop")

wilcox.test(dfnodes2$betweenness.x~ dfnodes2$Type.x, data=dfnodes2)
#Wilcoxon rank sum test with continuity correction
#data:  dfnodes2$betweenness.x by dfnodes2$Type.x
#W = 56.5, p-value = 0.9415
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(dfnodes2$betweenness.y~ dfnodes2$Type.x, data=dfnodes2)
#Wilcoxon rank sum test with continuity correction
#data:  dfnodes2$betweenness.y by dfnodes2$Type.x
#W = 44, p-value = 0.4542
#alternative hypothesis: true location shift is not equal to 0


