setwd("~/EMBER/Graph_Experiments")
dd<-read.table(file="Discovery_TDA_Features_Rows.csv",header=TRUE,sep=",")

library(igraph)
library(lsa)
library(WGCNA)
library(dunn.test)
library(netcom)
library(bluster)
d_cut<-0.5

#All data Graph
dd_ALL<-dd
dd_c<-t(dd_ALL)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_ALL$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-apply(dd_c,2,as.numeric)
d<-1-abs(cor(dd_c))
d[d>d_cut]<-0
ALL_g<-graph.adjacency(d,mode="undirected",weighted=TRUE,diag=FALSE)


dd_ALL<-dd
dd_c<-t(dd_ALL)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_ALL$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-apply(dd_c,2,as.numeric)
pickSoftThreshold(dd_c,dataIsExpr = TRUE)
d<-abs(cor(dd_c))^4
TOM_sim<-TOMsimilarity(d)
colnames(TOM_sim)<-colnames(d)
rownames(TOM_sim)<-rownames(d)	
#ALL_g<-graph.adjacency(d,mode="undirected",weighted=TRUE,diag=FALSE)
ALL_g<-graph.adjacency(TOM_sim,mode="undirected",weighted=TRUE,diag=FALSE)
plot(ALL_g, edge.label=round(E(ALL_g)$weight, 3))

ALL_g_clust<-cluster_louvain(ALL_g,weights = edge_attr(ALL_g)$weight)
ALL_g_clust<-cluster_leading_eigen(ALL_g,weights = edge_attr(ALL_g)$weight)
ALL_g_clust<-cluster_optimal(ALL_g,weights = edge_attr(ALL_g)$weight)
ALL_g_clust<-cluster_fast_greedy(ALL_g,weights = edge_attr(ALL_g)$weight)
plot(ALL_g_clust,ALL_g,mark.groups=NULL,main="ALL")
communities(ALL_g_clust)

dd_ALL_t<-t(dd_ALL)
dd_ALL_t<-data.frame(dd_ALL_t)
names(dd_ALL_t)<-dd_ALL_t[nrow(dd_ALL_t),]
dd_ALL_t<-dd_ALL_t[-nrow(dd_ALL_t),]

Grp<-c()
for(i in 1:nrow(dd_ALL_t)){
Grp[i]<-paste(gsub(substr(rownames(dd_ALL_t)[i],1,8),"",rownames(dd_ALL_t)[i]))
}

dd_ALL_t<-data.frame(apply(dd_ALL_t,2,as.numeric))

eigen_features<-list()
for(i in 1:length(communities(ALL_g_clust))){
if(length(communities(ALL_g_clust)[[i]]) >= 2)	{eigen_features[[i]]<-list(prcomp(dd_ALL_t[,names(dd_ALL_t) %in% communities(ALL_g_clust)[[i]]],center=TRUE,scale.=FALSE)$x[,1],length(communities(ALL_g_clust)[[i]]))}
if(length(communities(ALL_g_clust)[[i]]) < 2)	{eigen_features[[i]]<-list(dd_ALL_t[,names(dd_ALL_t) %in% communities(ALL_g_clust)[[i]]],length(communities(ALL_g_clust)[[i]]))}   
}

dd_eigen_features<-NULL
for(i in 1:length(communities(ALL_g_clust))){
dd_eigen_features<-cbind(dd_eigen_features,eigen_features[[i]][[1]])
}	
dd_eigen_features<-data.frame(dd_eigen_features)
dd_eigen_features$Grp<-Grp

p_vals<-apply(dd_eigen_features[,-ncol(dd_eigen_features)], 2, function(x) kruskal.test(x,dd_eigen_features[,ncol(dd_eigen_features)])$p.value)
which(p_vals<0.05)
			  
for(i in which(p_vals<0.05)){
print(communities(ALL_g_clust)[[i]])
}
for(i in which(p_vals<0.05)){			  
dunn.test(x=dd_eigen_features[,i], g= dd_eigen_features[,ncol(dd_eigen_features)], method = 'bonferroni', altp = TRUE)
}
			  
dd_feature_red<-data.frame(dd_eigen_features[,3],dd_eigen_features$Grp)
names(dd_feature_red)<-c("Eigen_Feature","Group")	
library(ggpubr)			  
ggviolin(dd_feature_red, x = "Group",
          y = c("Eigen_Feature"),
          combine = TRUE, 
          color = "Group", palette = "jco",
          ylab = "Eigen Feature Score", 
          add = "median_iqr")
library(pROC)
sc1<-roc(as.numeric(dd_eigen_features$Grp=="Healthy"),dd_eigen_features[,3])
plot(sc1,print.auc=TRUE)
			  
library(TDAstats)
eigen_phom<-calculate_homology(dd_eigen_features[,-ncol(dd_eigen_features)])			  
plot_barcode(eigen_phom)			  
plot_persist(eigen_phom)	

		  
#By Group Graphs
dd_H<-dd[,c(grep("Healthy",names(dd)),140)]
dd_c<-t(dd_H)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_H$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-apply(dd_c,2,as.numeric)
d<-1-abs(cor(dd_c))
d[d>d_cut]<-0
TOM_sim<-TOMsimilarity(d)
colnames(TOM_sim)<-colnames(d)
rownames(TOM_sim)<-rownames(d)			  
#H_g<-graph.adjacency(d,mode="undirected",weighted=TRUE,diag=FALSE)
H_g<-graph.adjacency(TOM_sim,mode="undirected",weighted=TRUE,diag=FALSE)
coords_H_g<-layout.fruchterman.reingold(H_g)		  

dd_C<-dd[,c(grep("COPD",names(dd)),140)]
dd_c<-t(dd_C)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_C$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-apply(dd_c,2,as.numeric)
d<-1-abs(cor(dd_c))
d[d>d_cut]<-0
TOM_sim<-TOMsimilarity(d)
colnames(TOM_sim)<-colnames(d)
rownames(TOM_sim)<-rownames(d)			  
#C_g<-graph.adjacency(d,mode="undirected",weighted=TRUE,diag=FALSE)
C_g<-graph.adjacency(TOM_sim,mode="undirected",weighted=TRUE,diag=FALSE)			  
coords_C_g<-layout.fruchterman.reingold(C_g)
			  
dd_A<-dd[,c(grep("Asthma",names(dd)),140)]	  
dd_c<-t(dd_A)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_A$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-apply(dd_c,2,as.numeric)
d<-1-abs(cor(dd_c))
d[d>d_cut]<-0
TOM_sim<-TOMsimilarity(d)
colnames(TOM_sim)<-colnames(d)
rownames(TOM_sim)<-rownames(d)			  
#A_g<-graph.adjacency(d,mode="undirected",weighted=TRUE,diag=FALSE)
A_g<-graph.adjacency(TOM_sim,mode="undirected",weighted=TRUE,diag=FALSE)
coords_A_g<-layout.fruchterman.reingold(A_g)			  

dd_P<-dd[,c(grep("Pneumonia",names(dd)),140)]
dd_c<-t(dd_P)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_P$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-apply(dd_c,2,as.numeric)
d<-1-abs(cor(dd_c))
d[d>d_cut]<-0
TOM_sim<-TOMsimilarity(d)
colnames(TOM_sim)<-colnames(d)
rownames(TOM_sim)<-rownames(d)			  
#P_g<-graph.adjacency(d,mode="undirected",weighted=TRUE,diag=FALSE)
P_g<-graph.adjacency(TOM_sim,mode="undirected",weighted=TRUE,diag=FALSE)			  
coords_P_g<-layout.fruchterman.reingold(P_g)
			  
dd_HF<-dd[,c(grep("HF",names(dd)),140)]
dd_c<-t(dd_HF)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_HF$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-apply(dd_c,2,as.numeric)
d<-1-abs(cor(dd_c))
d[d>d_cut]<-0
TOM_sim<-TOMsimilarity(d)
colnames(TOM_sim)<-colnames(d)
rownames(TOM_sim)<-rownames(d)			  
#HF_g<-graph.adjacency(d,mode="undirected",weighted=TRUE,diag=FALSE)
HF_g<-graph.adjacency(TOM_sim,mode="undirected",weighted=TRUE,diag=FALSE)			  
coords_HF_g<-layout.fruchterman.reingold(HF_g)

#Louvain
setwd("~/EMBER/Graph_Experiments/Louvain_clusters")
H_g_clust<-cluster_louvain(H_g,weights = edge_attr(H_g)$weight)
P_g_clust<-cluster_louvain(P_g,weights = edge_attr(P_g)$weight)
A_g_clust<-cluster_louvain(A_g,weights = edge_attr(A_g)$weight)
C_g_clust<-cluster_louvain(C_g,weights = edge_attr(C_g)$weight)
HF_g_clust<-cluster_louvain(HF_g,weights = edge_attr(HF_g)$weight)
tiff(filename = "Louvain_FR_cut_05.tiff" ,units="in", width=17, height=11, res=300,compression = "lzw")
par(mfrow=c(2,3))
plot(H_g_clust,H_g,mark.groups=NULL,main="Healthy",layout=coords_H_g)
plot(P_g_clust,P_g,mark.groups=NULL,main="Pneumonia",,layout=coords_P_g)
plot(A_g_clust,A_g,mark.groups=NULL,main="Asthma",layout=coords_A_g)
plot(C_g_clust,C_g,mark.groups=NULL,main="COPD",layout=coords_C_g)
plot(HF_g_clust,HF_g,mark.groups=NULL,main="HF",layout=coords_HF_g)
dev.off()
			  
Healthy_g_comm<-communities(H_g_clust)
Pneumonia_g_comm<-communities(P_g_clust)
Asthma_g_comm<-communities(A_g_clust)
COPD_g_comm<-communities(C_g_clust)
HF_g_comm<-communities(HF_g_clust)
capture.output(Healthy_g_comm, file = "Healthy_Clusters.txt")
capture.output(Pneumonia_g_comm, file = "Pneumonia_Clusters.txt")
capture.output(Asthma_g_comm, file = "Asthma_Clusters.txt")
capture.output(COPD_g_comm, file = "COPD_Clusters.txt")
capture.output(HF_g_comm, file = "HF_Clusters.txt")			  

library(alphahull)	
library(geometry)
library(shotGroups)
library(dplyr)
library(concaveman)			  
library(sf)	
library(lwgeom)
#Asthma cluster 6 in Asthma			  
#dd_convex<-data.frame(cbind(get.vertex.attribute(A_g)$name,layout_(A_g,with_fr())))		  
dd_convex<-data.frame(cbind(get.vertex.attribute(A_g)$name,coords_A_g))	
names(dd_convex)<-c("Feature","X","Y")
communities(A_g_clust)$'6'
dd_convex<-dd_convex[dd_convex$Feature %in% communities(A_g_clust)$'6',][,c(2,3)]			  
pnts <- dd_convex %>%
  st_as_sf(coords = c("X", "Y"))
polygon_A_g <- concaveman(pnts)

plot(polygon_A_g, reset = FALSE,main="Asthma Graph, Asthma Orange Cluster")
plot(pnts,add = TRUE)
plot(st_minimum_bounding_circle(polygon_A_g),add = TRUE)			  
plot(st_make_grid(st_bbox(polygon_A_g),n=1),add=TRUE)		  
A_g_circle_compact<-st_area(polygon_A_g)/st_area(st_minimum_bounding_circle(polygon_A_g))			  			  
A_g_box_compact<-(st_bbox(polygon_A_g)$xmax-st_bbox(polygon_A_g)$xmin)/(st_bbox(polygon_A_g)$ymax-st_bbox(polygon_A_g)$ymin)
			  
#Asthma cluster 6 in Healthy				  			  
#dd_convex<-data.frame(cbind(get.vertex.attribute(H_g)$name,layout_(H_g,with_fr())))	
dd_convex<-data.frame(cbind(get.vertex.attribute(H_g)$name,coords_H_g))				  
names(dd_convex)<-c("Feature","X","Y")
communities(A_g_clust)$'6'
dd_convex<-dd_convex[dd_convex$Feature %in% communities(A_g_clust)$'6',][,c(2,3)]
pnts <- dd_convex %>%
  st_as_sf(coords = c("X", "Y"))
polygon_H_g <- concaveman(pnts)

plot(polygon_H_g, reset = FALSE,main="Healthy Graph, Asthma Orange Cluster")
plot(pnts, add = TRUE)
plot(st_minimum_bounding_circle(polygon_H_g),add = TRUE)
plot(st_make_grid(st_bbox(polygon_H_g),n=1),add=TRUE)			  
H_g_circle_compact<-st_area(polygon_H_g)/st_area(st_minimum_bounding_circle(polygon_H_g))			  			  
H_g_box_compact<-(st_bbox(polygon_H_g)$xmax-st_bbox(polygon_H_g)$xmin)/(st_bbox(polygon_H_g)$ymax-st_bbox(polygon_H_g)$ymin)	
			  
#Asthma cluster 6 in COPD				  			  
#dd_convex<-data.frame(cbind(get.vertex.attribute(C_g)$name,layout_(C_g,with_fr())))		  
dd_convex<-data.frame(cbind(get.vertex.attribute(C_g)$name,coords_C_g))		  
names(dd_convex)<-c("Feature","X","Y")
communities(A_g_clust)$'6'
dd_convex<-dd_convex[dd_convex$Feature %in% communities(A_g_clust)$'6',][,c(2,3)]
pnts <- dd_convex %>%
  st_as_sf(coords = c("X", "Y"))
polygon_C_g <- concaveman(pnts)

plot(polygon_C_g, reset = FALSE,main="COPD Graph, Asthma Orange Cluster")
plot(pnts, add = TRUE)
plot(st_minimum_bounding_circle(polygon_C_g),add = TRUE)
plot(st_make_grid(st_bbox(polygon_C_g),n=1),add=TRUE)			  
C_g_circle_compact<-st_area(polygon_C_g)/st_area(st_minimum_bounding_circle(polygon_C_g))			  			  
C_g_box_compact<-(st_bbox(polygon_C_g)$xmax-st_bbox(polygon_C_g)$xmin)/(st_bbox(polygon_C_g)$ymax-st_bbox(polygon_C_g)$ymin)
			  
#Asthma cluster 6 in Pneumonia				  			  
#dd_convex<-data.frame(cbind(get.vertex.attribute(P_g)$name,layout_(P_g,with_fr())))		  
dd_convex<-data.frame(cbind(get.vertex.attribute(P_g)$name,coords_P_g))	
names(dd_convex)<-c("Feature","X","Y")
communities(A_g_clust)$'6'
dd_convex<-dd_convex[dd_convex$Feature %in% communities(A_g_clust)$'6',][,c(2,3)]
pnts <- dd_convex %>%
  st_as_sf(coords = c("X", "Y"))
polygon_P_g <- concaveman(pnts)

plot(polygon_P_g, reset = FALSE,main="Pneumonia Graph, Asthma Orange Cluster")
plot(pnts, add = TRUE)
plot(st_minimum_bounding_circle(polygon_P_g),add = TRUE)
plot(st_make_grid(st_bbox(polygon_P_g),n=1),add=TRUE)			  
P_g_circle_compact<-st_area(polygon_P_g)/st_area(st_minimum_bounding_circle(polygon_P_g))			  			  
P_g_box_compact<-(st_bbox(polygon_P_g)$xmax-st_bbox(polygon_P_g)$xmin)/(st_bbox(polygon_P_g)$ymax-st_bbox(polygon_P_g)$ymin)
			  
#Asthma cluster 6 in HF				  			  
#dd_convex<-data.frame(cbind(get.vertex.attribute(HF_g)$name,layout_(HF_g,with_fr())))
dd_convex<-data.frame(cbind(get.vertex.attribute(HF_g)$name,coords_HF_g))			  
names(dd_convex)<-c("Feature","X","Y")
communities(A_g_clust)$'6'
dd_convex<-dd_convex[dd_convex$Feature %in% communities(A_g_clust)$'6',][,c(2,3)]
pnts <- dd_convex %>%
  st_as_sf(coords = c("X", "Y"))
polygon_HF_g <- concaveman(pnts)

plot(polygon_HF_g, reset = FALSE,main="HF Graph, Asthma Orange Cluster")
plot(pnts, add = TRUE)
plot(st_minimum_bounding_circle(polygon_HF_g),add = TRUE)
plot(st_make_grid(st_bbox(polygon_HF_g),n=1),add=TRUE)			  
HF_g_circle_compact<-st_area(polygon_HF_g)/st_area(st_minimum_bounding_circle(polygon_HF_g))			  			  
HF_g_box_compact<-(st_bbox(polygon_HF_g)$xmax-st_bbox(polygon_HF_g)$xmin)/(st_bbox(polygon_HF_g)$ymax-st_bbox(polygon_HF_g)$ymin)
			  
dd_circle_compact<-data.frame(A_g_circle_compact,H_g_circle_compact,P_g_circle_compact,C_g_circle_compact,HF_g_circle_compact)
dd_box_compact<-data.frame(A_g_box_compact,H_g_box_compact,P_g_box_compact,C_g_box_compact,HF_g_box_compact)			  

library(diceR)			  
dd_H_clust<-dd_H[dd_H$Feature %in% communities(A_g_clust)$'6',]			  
rownames(dd_H_clust)<-dd_H_clust$Feature
H_g_dice_compact<-compactness(dd_H_clust[,-ncol(dd_H_clust)],rep("Orange",nrow(dd_H_clust)))
			  
dd_A_clust<-dd_A[dd_A$Feature %in% communities(A_g_clust)$'6',]			  
rownames(dd_A_clust)<-dd_A_clust$Feature
A_g_dice_compact<-compactness(dd_A_clust[,-ncol(dd_A_clust)],rep("Orange",nrow(dd_A_clust)))	
			  
dd_C_clust<-dd_C[dd_C$Feature %in% communities(A_g_clust)$'6',]			  
rownames(dd_C_clust)<-dd_C_clust$Feature
C_g_dice_compact<-compactness(dd_C_clust[,-ncol(dd_C_clust)],rep("Orange",nrow(dd_C_clust)))	
			  
dd_P_clust<-dd_P[dd_P$Feature %in% communities(A_g_clust)$'6',]			  
rownames(dd_P_clust)<-dd_P_clust$Feature
P_g_dice_compact<-compactness(dd_P_clust[,-ncol(dd_P_clust)],rep("Orange",nrow(dd_P_clust)))			  

dd_HF_clust<-dd_HF[dd_HF$Feature %in% communities(A_g_clust)$'6',]			  
rownames(dd_HF_clust)<-dd_HF_clust$Feature
HF_g_dice_compact<-compactness(dd_HF_clust[,-ncol(dd_HF_clust)],rep("Orange",nrow(dd_HF_clust)))
			  
dd_dice_compact<-data.frame(A_g_dice_compact,H_g_dice_compact,P_g_dice_compact,C_g_dice_compact,HF_g_dice_compact)			  
which.min(dd_dice_compact[1,])
			  
tiff(filename = "Louvain_FR_cut_05_Asthma_Clusters_Imposed.tiff" ,units="in", width=17, height=11, res=300,compression = "lzw")
par(mfrow=c(2,3))
plot(A_g_clust,H_g,mark.groups=NULL,vertex.color=membership(A_g_clust),main="Healthy")
plot(A_g_clust,P_g,mark.groups=NULL,vertex.color=membership(A_g_clust),main="Pneumonia")
plot(A_g_clust,A_g,mark.groups=NULL,vertex.color=membership(A_g_clust),main="Asthma")
plot(A_g_clust,C_g,mark.groups=NULL,vertex.color=membership(A_g_clust),main="COPD")
plot(A_g_clust,HF_g,mark.groups=NULL,vertex.color=membership(A_g_clust),main="HF")
dev.off()
			  
tiff(filename = "Louvain_FR_cut_05_Healthy_Clusters_Imposed.tiff" ,units="in", width=17, height=11, res=300,compression = "lzw")
par(mfrow=c(2,3))
plot(H_g_clust,H_g,mark.groups=NULL,vertex.color=membership(H_g_clust),main="Healthy")
plot(H_g_clust,P_g,mark.groups=NULL,vertex.color=membership(H_g_clust),main="Pneumonia")
plot(H_g_clust,A_g,mark.groups=NULL,vertex.color=membership(H_g_clust),main="Asthma")
plot(H_g_clust,C_g,mark.groups=NULL,vertex.color=membership(H_g_clust),main="COPD")
plot(H_g_clust,HF_g,mark.groups=NULL,vertex.color=membership(H_g_clust),main="HF")
dev.off()	  

tiff(filename = "Louvain_FR_cut_05_Pneumonia_Clusters_Imposed.tiff" ,units="in", width=17, height=11, res=300,compression = "lzw")
par(mfrow=c(2,3))
plot(P_g_clust,H_g,mark.groups=NULL,vertex.color=membership(P_g_clust),main="Healthy")
plot(P_g_clust,P_g,mark.groups=NULL,vertex.color=membership(P_g_clust),main="Pneumonia")
plot(P_g_clust,A_g,mark.groups=NULL,vertex.color=membership(P_g_clust),main="Asthma")
plot(P_g_clust,C_g,mark.groups=NULL,vertex.color=membership(P_g_clust),main="COPD")
plot(P_g_clust,HF_g,mark.groups=NULL,vertex.color=membership(P_g_clust),main="HF")
dev.off()	
			  
tiff(filename = "Louvain_FR_cut_05_COPD_Clusters_Imposed.tiff" ,units="in", width=17, height=11, res=300,compression = "lzw")
par(mfrow=c(2,3))
plot(C_g_clust,H_g,mark.groups=NULL,vertex.color=membership(C_g_clust),main="Healthy")
plot(C_g_clust,P_g,mark.groups=NULL,vertex.color=membership(C_g_clust),main="Pneumonia")
plot(C_g_clust,A_g,mark.groups=NULL,vertex.color=membership(C_g_clust),main="Asthma")
plot(C_g_clust,C_g,mark.groups=NULL,vertex.color=membership(C_g_clust),main="COPD")
plot(C_g_clust,HF_g,mark.groups=NULL,vertex.color=membership(C_g_clust),main="HF")
dev.off()
			  
tiff(filename = "Louvain_FR_cut_05_HF_Clusters_Imposed.tiff" ,units="in", width=17, height=11, res=300,compression = "lzw")
par(mfrow=c(2,3))
plot(HF_g_clust,H_g,mark.groups=NULL,vertex.color=membership(HF_g_clust),main="Healthy")
plot(HF_g_clust,P_g,mark.groups=NULL,vertex.color=membership(HF_g_clust),main="Pneumonia")
plot(HF_g_clust,A_g,mark.groups=NULL,vertex.color=membership(HF_g_clust),main="Asthma")
plot(HF_g_clust,C_g,mark.groups=NULL,vertex.color=membership(HF_g_clust),main="COPD")
plot(HF_g_clust,HF_g,mark.groups=NULL,vertex.color=membership(HF_g_clust),main="HF")
dev.off()				  
			  
x_mod<-c(modularity(H_g_clust),modularity(P_g_clust),modularity(A_g_clust),modularity(C_g_clust),modularity(HF_g_clust))
x_graph<-c("Healthy","Pneumonia","Asthma","COPD","HF")
tiff(filename = "Modularity_Louvain_FR_cut_05.tiff" ,units="in", width=17, height=11, res=300,compression = "lzw")			  
barplot(x_mod~x_graph,xlab="Clustering on Graph",ylab="Modularity")
dev.off()



		  			
u<-c("A_g","H_g","P_g","C_g","HF_g")
v<-u
graphs_list<-list(A_g,H_g,C_g,P_g,HF_g)
names(graphs_list)<-u
dd_graphs<-expand.grid(u,v,stringsAsFactors=FALSE)	
			  
for(i in 1:nrow(dd_graphs)){
cat(grep(dd_graphs$Var1[i],names(graphs_list)),grep(dd_graphs$Var2[i],names(graphs_list)),fill=TRUE)
dd_graphs$x_g[i]<-grep(dd_graphs$Var1[i],names(graphs_list))
dd_graphs$y_g[i]<-grep(dd_graphs$Var2[i],names(graphs_list))	
}			  

cons_struc<-c()
cons_struc_flip<-c()			  
for(i in 1:nrow(dd_graphs)){
x <- graphs_list[[dd_graphs$x_g[i]]]
y <- graphs_list[[dd_graphs$y_g[i]]]			  
cons_struc[i]<-eval(substitute(ics(as.matrix(get.adjacency(xx)),as.matrix(get.adjacency(yy)), align(as.matrix(get.adjacency(xx)),as.matrix(get.adjacency(yy)),characterization = "gini")$alignment),list(xx=x,yy=y)))
cons_struc_flip[i]<-eval(substitute(ics(as.matrix(get.adjacency(xx)),as.matrix(get.adjacency(yy)), align(as.matrix(get.adjacency(xx)),as.matrix(get.adjacency(yy)),characterization = "gini")$alignment,flip=TRUE),list(xx=x,yy=y)))
}	
dd_graphs$Conserved_Structure<-cons_struc
write.table(dd_graphs,file="Graphs_Conserved_Structure.csv",row.names=FALSE,sep=",")			  
dd_graphs$Conserved_Structure_Flip<-cons_struc_flip	
			  
#All Graphs with own clustering	
par(mfrow=c(2,3))
ratio_A_g <- pairwiseModularity(A_g, A_g_clust$membership, as.ratio=TRUE)			  
library(pheatmap)
p1<-pheatmap(log2(ratio_A_g+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100),main="Asthma")			  
			  
ratio_H_g <- pairwiseModularity(H_g, H_g_clust$membership, as.ratio=TRUE)			  
library(pheatmap)
p2<-pheatmap(log2(ratio_H_g+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100),main="Healthy")				  

ratio_C_g <- pairwiseModularity(C_g, C_g_clust$membership, as.ratio=TRUE)			  
library(pheatmap)
p3<-pheatmap(log2(ratio_C_g+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100),main="COPD")
			  
ratio_P_g <- pairwiseModularity(P_g, P_g_clust$membership, as.ratio=TRUE)			  
library(pheatmap)
p4<-pheatmap(log2(ratio_P_g+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100),main="Pneumonia")	
			  
ratio_HF_g <- pairwiseModularity(HF_g, HF_g_clust$membership, as.ratio=TRUE)			  
library(pheatmap)
p5<-pheatmap(log2(ratio_HF_g+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100),main="HF")

sum(log2(ratio_A_g+1)[upper.tri(log2(ratio_A_g+1))],na.rm=TRUE)
sum(log2(ratio_H_g+1)[upper.tri(log2(ratio_H_g+1))],na.rm=TRUE)
sum(log2(ratio_C_g+1)[upper.tri(log2(ratio_C_g+1))],na.rm=TRUE)
sum(log2(ratio_P_g+1)[upper.tri(log2(ratio_P_g+1))],na.rm=TRUE)
sum(log2(ratio_HF_g+1)[upper.tri(log2(ratio_HF_g+1))],na.rm=TRUE)			  

#Asthma Imposed
ratio_A_g <- pairwiseModularity(A_g, A_g_clust$membership, as.ratio=TRUE)			  
library(pheatmap)
pheatmap(log2(ratio_A_g+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100))			  
			  
ratio_H_g <- pairwiseModularity(H_g, A_g_clust$membership, as.ratio=TRUE)			  
library(pheatmap)
pheatmap(log2(ratio_H_g+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100))				  

ratio_C_g <- pairwiseModularity(C_g, A_g_clust$membership, as.ratio=TRUE)			  
library(pheatmap)
pheatmap(log2(ratio_C_g+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100))
			  
ratio_P_g <- pairwiseModularity(P_g, A_g_clust$membership, as.ratio=TRUE)			  
library(pheatmap)
pheatmap(log2(ratio_P_g+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100))	
			  
ratio_HF_g <- pairwiseModularity(HF_g, A_g_clust$membership, as.ratio=TRUE)			  
library(pheatmap)
pheatmap(log2(ratio_HF_g+1), cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("white", "blue"))(100))			  
			  
sum(log2(ratio_H_g+1)[upper.tri(log2(ratio_H_g+1))],na.rm=TRUE)-sum(log2(ratio_A_g+1)[upper.tri(log2(ratio_A_g+1))],na.rm=TRUE)
sum(log2(ratio_C_g+1)[upper.tri(log2(ratio_C_g+1))],na.rm=TRUE)-sum(log2(ratio_A_g+1)[upper.tri(log2(ratio_A_g+1))],na.rm=TRUE)
sum(log2(ratio_P_g+1)[upper.tri(log2(ratio_P_g+1))],na.rm=TRUE)-sum(log2(ratio_A_g+1)[upper.tri(log2(ratio_A_g+1))],na.rm=TRUE)
sum(log2(ratio_HF_g+1)[upper.tri(log2(ratio_HF_g+1))],na.rm=TRUE)-sum(log2(ratio_A_g+1)[upper.tri(log2(ratio_A_g+1))],na.rm=TRUE)			  
			  
#By Group Graphs
dd_H<-dd[,c(grep("Healthy",names(dd)),140)]
dd_c<-t(dd_H)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_H$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-dd_c[,names(dd_c) %in% Asthma_g_comm$'8']			  
dd_c<-apply(dd_c,2,as.numeric)
dd_pca_H<-prcomp(dd_c)			  

dd_C<-dd[,c(grep("COPD",names(dd)),140)]
dd_c<-t(dd_C)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_C$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-dd_c[,names(dd_c) %in% Asthma_g_comm$'8']			  
dd_c<-apply(dd_c,2,as.numeric)
dd_pca_C<-prcomp(dd_c)			  

dd_A<-dd[,c(grep("Asthma",names(dd)),140)]	  
dd_c<-t(dd_A)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_A$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-dd_c[,names(dd_c) %in% Asthma_g_comm$'8']			  
dd_c<-apply(dd_c,2,as.numeric)
dd_pca_A<-prcomp(dd_c)	

dd_P<-dd[,c(grep("Pneumonia",names(dd)),140)]
dd_c<-t(dd_P)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_P$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-dd_c[,names(dd_c) %in% Asthma_g_comm$'8']			  
dd_c<-apply(dd_c,2,as.numeric)
dd_pca_P<-prcomp(dd_c)			  

dd_HF<-dd[,c(grep("HF",names(dd)),140)]
dd_c<-t(dd_HF)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_HF$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-dd_c[,names(dd_c) %in% Asthma_g_comm$'8']			  
dd_c<-apply(dd_c,2,as.numeric)
dd_pca_HF<-prcomp(dd_c)

			  			 
par(mfrow=c(2,3))			  
hist(dd_pca_H$x[,1],xlim=c(-10,10),xlab="PC Score Grey Cluster",main="Healthy")
abline(v=median(dd_pca_H$x[,1]),col="red")			  
hist(dd_pca_C$x[,1],xlim=c(-10,10),xlab="PC Score Grey Asthma Cluster",main="COPD")	
abline(v=median(dd_pca_C$x[,1]),col="red")			  
hist(dd_pca_A$x[,1],xlim=c(-10,10),xlab="PC Score Grey Asthma Cluster",main="Asthma")
abline(v=median(dd_pca_A$x[,1]),col="red")			  
hist(dd_pca_P$x[,1],xlim=c(-10,10),xlab="PC Score Grey Asthma Cluster",main="Pneumonia")
abline(v=median(dd_pca_P$x[,1]),col="red")			  
hist(dd_pca_HF$x[,1],xlim=c(-10,10),xlab="PC Score Grey Asthma Cluster",main="HF")	
abline(v=median(dd_pca_HF$x[,1]),col="red")			  
plot(A_g_clust,A_g,mark.groups=NULL,main="Asthma")

dd_pc_score<-c(median(dd_pca_H$x[,1]),median(dd_pca_C$x[,1]),median(dd_pca_A$x[,1]),median(dd_pca_P$x[,1]),median(dd_pca_HF$x[,1]))
x_graph<-c("Healthy","COPD","Asthma","Pneumonia","HF")
par(mfrow=c(2,2))			  
barplot(dd_pc_score~x_graph,xlab="Clustering on Graph",ylab="PC Score Grey Asthma Cluster")			  
plot(A_g_clust,A_g,mark.groups=NULL,main="Asthma")
			  
pc_score<-c(dd_pca_H$x[,1],dd_pca_C$x[,1],dd_pca_A$x[,1],dd_pca_P$x[,1],dd_pca_HF$x[,1])
Grp<-c(rep("Healthy",27),rep("COPD",29),rep("Asthma",33),rep("Pneumonia",28),rep("HF",22))	

library(pROC)
sc1<-roc(as.numeric(Grp=="Asthma"),pc_score)
plot(sc1,print.auc=TRUE)


			  
#Spectral
setwd("~/EMBER/Graph_Experiments/Spectral_clusters")
H_g_clust<-cluster_leading_eigen(H_g,weights = edge_attr(H_g)$weight)
P_g_clust<-cluster_leading_eigen(P_g,weights = edge_attr(P_g)$weight)
A_g_clust<-cluster_leading_eigen(A_g,weights = edge_attr(A_g)$weight)
C_g_clust<-(cluster_leading_eigenC_g,weights = edge_attr(C_g)$weight)
HF_g_clust<-cluster_leading_eigen(HF_g,weights = edge_attr(HF_g)$weight)
tiff(filename = "Leading_Eigen_FR_cut_05.tiff" ,units="in", width=17, height=11, res=300,compression = "lzw")
par(mfrow=c(2,3))
plot(H_g_clust,H_g,mark.groups=NULL,main="Healthy")
plot(P_g_clust,P_g,mark.groups=NULL,main="Pneumonia")
plot(A_g_clust,A_g,mark.groups=NULL,main="Asthma")
plot(C_g_clust,C_g,mark.groups=NULL,main="COPD")
plot(HF_g_clust,HF_g,mark.groups=NULL,main="HF")
dev.off()

x_mod<-c(modularity(H_g_clust),modularity(P_g_clust),modularity(A_g_clust),modularity(C_g_clust),modularity(HF_g_clust))
x_graph<-c("Healthy","Pneumonia","Asthma","COPD","HF")
barplot(x_mod~x_graph,xlab="Clustering on Graph",ylab="Modularity")

Healthy_g_comm<-communities(H_g_clust)
Pneumonia_g_comm<-communities(P_g_clust)
Asthma_g_comm<-communities(A_g_clust)
COPD_g_comm<-communities(C_g_clust)
HF_g_comm<-communities(HF_g_clust)
capture.output(Healthy_g_comm, file = "Healthy_Clusters.txt")
capture.output(Pneumonia_g_comm, file = "Pneumonia_Clusters.txt")
capture.output(Asthma_g_comm, file = "Asthma_Clusters.txt")
capture.output(COPD_g_comm, file = "COPD_Clusters.txt")
capture.output(HF_g_comm, file = "HF_Clusters.txt")


#Optimal
setwd("~/EMBER/Graph_Experiments/Optimal_clusters")
H_g_clust<-cluster_optimal(H_g,weights = edge_attr(H_g)$weight)
P_g_clust<-cluster_optimal(P_g,weights = edge_attr(P_g)$weight)
A_g_clust<-cluster_optimal(A_g,weights = edge_attr(A_g)$weight)
C_g_clust<-cluster_optimal(C_g,weights = edge_attr(C_g)$weight)
HF_g_clust<-cluster_optimal(HF_g,weights = edge_attr(HF_g)$weight)
tiff(filename = "Optimal_FR_cut_05.tiff" ,units="in", width=17, height=11, res=300,compression = "lzw")
par(mfrow=c(2,3))
plot(H_g_clust,H_g,mark.groups=NULL,main="Healthy")
plot(P_g_clust,P_g,mark.groups=NULL,main="Pneumonia")
plot(A_g_clust,A_g,mark.groups=NULL,main="Asthma")
plot(C_g_clust,C_g,mark.groups=NULL,main="COPD")
plot(HF_g_clust,HF_g,mark.groups=NULL,main="HF")
dev.off()

x_mod<-c(modularity(H_g_clust),modularity(P_g_clust),modularity(A_g_clust),modularity(C_g_clust),modularity(HF_g_clust))
x_graph<-c("Healthy","Pneumonia","Asthma","COPD","HF")
barplot(x_mod~x_graph,xlab="Clustering on Graph",ylab="Modularity")

Healthy_g_comm<-communities(H_g_clust)
Pneumonia_g_comm<-communities(P_g_clust)
Asthma_g_comm<-communities(A_g_clust)
COPD_g_comm<-communities(C_g_clust)
HF_g_comm<-communities(HF_g_clust)
capture.output(Healthy_g_comm, file = "Healthy_Clusters.txt")
capture.output(Pneumonia_g_comm, file = "Pneumonia_Clusters.txt")
capture.output(Asthma_g_comm, file = "Asthma_Clusters.txt")
capture.output(COPD_g_comm, file = "COPD_Clusters.txt")
capture.output(HF_g_comm, file = "HF_Clusters.txt")

#Greedy
setwd("~/EMBER/Graph_Experiments/Greedy_clusters")
H_g_clust<-cluster_fast_greedy(H_g,weights = edge_attr(H_g)$weight)
P_g_clust<-cluster_fast_greedy(P_g,weights = edge_attr(P_g)$weight)
A_g_clust<-cluster_fast_greedy(A_g,weights = edge_attr(A_g)$weight)
C_g_clust<-cluster_fast_greedy(C_g,weights = edge_attr(C_g)$weight)
HF_g_clust<-cluster_fast_greedy(HF_g,weights = edge_attr(HF_g)$weight)
tiff(filename = "Greedy_FR_cut_05.tiff" ,units="in", width=17, height=11, res=300,compression = "lzw")
par(mfrow=c(2,3))
plot(H_g_clust,H_g,mark.groups=NULL,main="Healthy")
plot(P_g_clust,P_g,mark.groups=NULL,main="Pneumonia")
plot(A_g_clust,A_g,mark.groups=NULL,main="Asthma")
plot(C_g_clust,C_g,mark.groups=NULL,main="COPD")
plot(HF_g_clust,HF_g,mark.groups=NULL,main="HF")
dev.off()

x_mod<-c(modularity(H_g_clust),modularity(P_g_clust),modularity(A_g_clust),modularity(C_g_clust),modularity(HF_g_clust))
x_graph<-c("Healthy","Pneumonia","Asthma","COPD","HF")
barplot(x_mod~x_graph,xlab="Clustering on Graph",ylab="Modularity")

Healthy_g_comm<-communities(H_g_clust)
Pneumonia_g_comm<-communities(P_g_clust)
Asthma_g_comm<-communities(A_g_clust)
COPD_g_comm<-communities(C_g_clust)
HF_g_comm<-communities(HF_g_clust)
capture.output(Healthy_g_comm, file = "Healthy_Clusters.txt")
capture.output(Pneumonia_g_comm, file = "Pneumonia_Clusters.txt")
capture.output(Asthma_g_comm, file = "Asthma_Clusters.txt")
capture.output(COPD_g_comm, file = "COPD_Clusters.txt")
capture.output(HF_g_comm, file = "HF_Clusters.txt")
