#display an annotated tree on map
#this script imagines that geographic range position evolves according to minimum squared change parsimony
library(combinat)
library(ICSNP)
library(plyr)
library(dplyr)
library(phylotate)
library(ape)
library(phytools)
library(geosphere)
library(maps)
library(usmap)
library(ggplot2)
library(igraph)
library(grid)

#set working directory where annotated nexus file resides
setwd("~/Documents")

#phylotate functions to annote nodes
tree <- read_annotated("May3tree", format="nexus")
tree <- makeNodeLabel(tree, method = "number", prefix=NULL, nodeList = list())
tree$tip.label <- gsub("_.*","",tree$tip.label)
plot(tree, order=c("increasing"), show.node.label = TRUE)

collapse_identical_tips <- function(phy,tip_label){
  matching_tips <- which(phy$tip.label==tip_label)
  nt <- length(phy$tip.label) # number of tips in tree
  nm <- length(matching_tips) # Number of tips matching the label
  keep <- numeric(nm)
  
  cur_tip <- 1
  while(cur_tip<=nm){
    if(cur_tip == nm){
      keep[cur_tip] <- 1
      break
    }
    next_tip <- cur_tip + 1
    mrca_ <- getMRCA(phy,c(matching_tips[cur_tip],matching_tips[next_tip]))
    descendants <- getDescendants(phy, mrca_)
    descendant_tips <- descendants[descendants<=nt]
    if(all(descendant_tips %in% matching_tips)){
      keep[cur_tip] <- 1
      cur_tip <- cur_tip + length(descendant_tips)
    }else{
      keep[cur_tip] <- 1
      cur_tip <- cur_tip + 1
    }
  }
  to_drop <- matching_tips[!keep]
  new_phy <- drop.tip(phy,to_drop)
  return(new_phy)
}


pops <- tree$tip.label

for(p in pops) {
tree <- collapse_identical_tips(tree, p)
}

plot(tree)

#tree$edge is dataframe of all branches connecting nodes
#numbers 1-25 are the terminal taxa, in order of tree$tip.label
#numbers 26 and above are nodes in the tree (==node.label+25)

#file with lat_long
sites <- read.csv("Rad_Pop_LATLON_CP.csv")
pops <- c("AL2012", "ALBG", "AL79", "AR125", "AR56", "IA10", "IN46", "IN77", "KS60", "KY51", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OH119", "OH64", "OK61", "PA27", "TN19", "WI128", "WV72")
sites <- sites[sites$pop %in% pops,]
rownames(sites)<-sites$pop
sites <- sites[,2:3]
sites$pop <- rownames(sites)
sites[nrow(sites)+1,]<-c(39.61374,-79.11579,"MD5")
sites[nrow(sites)+1,]<-c(37.35337,-80.55216, "VA73")

#tree tip labels
tree$tip.label <- make.unique(tree$tip.label)
sites_big <- data.frame(ID=tree$tip.label, pop=gsub("\\\\..*","",tree$tip.label))
sites_big <- merge(sites_big, sites, by = c("pop", "pop"), all = TRUE)
rownames(sites_big) <- sites_big$ID
sites_big <- sites_big[,3:4]
sites_big$lat <- as.numeric(sites_big$lat)
sites_big$long <- as.numeric(sites_big$long)

#Reconstructed lat, long of each ancestral node  (ape package, 'ace' function, ML option)
old  <- phylomorphospace(tree, sites_big)

#root to tip distance begins with df$node==26
#tree on the map 
obj<- phylo.to.map(tree, sites_big, type="direct", database="state", xlim=c(-96,-78), ylim=c(32,46))
plot(obj, rotate=T, show.node.label = TRUE)

#pretty ggplot map
#how to place tree on it?
#still in progress
nodes <- df[26:nrow(df),]
leaves <- df[1:25,]
states <- map_data("state")
apple <- read.csv("apple_poly.csv")
ggplot(data = subset(states,long > -100 & long < -75 & lat > 25 & lat < 50)) + theme_bw() + ylim(c(25,50)) + xlim(c(-100,-75))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_polygon(data=apple, aes(long,lat), fill="dimgrey", alpha=1) +
  geom_polygon(data=extent, aes(long,lat), fill="deepskyblue2", alpha=0.2) +
  geom_point(data=nodes, aes(long, lat), color="black") +
  geom_point(data=sites, aes(long, lat), color="red") +
  coord_fixed(1.3) +
  xlab("longitude") +
  ylab("latitude") +
  theme(legend.position="none") +
  guides(fill=FALSE)  # do this to leave off the color legend


