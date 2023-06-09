#Written by: Dylan DeBaun

#load libraries
library('MSCquartets')
library('ape')
library('gtools')
library('stringr')

#load species tree
snake_tree <-read.tree("Tree_Dated_Point_PL_Astral_topology_130")

#load gene trees
genedata_ind = "allindiv.treefile"
gtrees_ind <- read.tree(genedata_ind)

#INPUTS
species = read.delim("",header  = F) #list of species, with the last species being the outgroup species
outgroup = species[length(species)]
indivs = read.delim("",header = F) # list of individuals
setwd("") #where we want the outfiles to go
group = "" #prefix for the clade

#Create the quartet file for NANUQ
b <- quartetTable(gtrees_ind,taxonnames=indivs,random = 0) #make quartet table at individual level
b <- b[rowSums(b) !=4, colSums(b) !=0] #edit it to look the way we want for snaq
# quartet file at the individual level
write.csv(b, paste0(group,"_fullindivCFs.csv"), row.names = F) 

#Create the quartet file for SnaQ
x<-b
y<- matrix(rep(0,dim(x)[1]*8),ncol = 8,nrow=dim(x)[1])
colnames(y)<-c("t1","t2","t3","t4","CF12_34","CF13_24","CF14_23","ngenes")
for(i in 1:dim(x)[1]){
  count=0
  for(j in 1:dim(x)[2]){
    if(x[i,j] == 1){
      count = count+1
      y[i,count] = colnames(x)[j]
    }
  }
  y[i,"ngenes"] = x[i,dim(x)[2]-2] + x[i,dim(x)[2]-1]+x[i,dim(x)[2]]
  y[i,"CF12_34"] = as.numeric(x[i,dim(x)[2]-2])/as.numeric(y[i,"ngenes"])
  y[i,"CF13_24"] = as.numeric(x[i,dim(x)[2]-1])/as.numeric(y[i,"ngenes"])
  y[i,"CF14_23"] = as.numeric(x[i,dim(x)[2]])/as.numeric(y[i,"ngenes"])
}
write.csv(y, paste0(group,"_individualsCFs.csv"), row.names = F)

#create the clade's tree from the species tree for SnaQ
tree_new = keep.tip(snake_tree, species)
tree_new$edge.length<-NULL #we dont need branch length information/etc
tree_new$node.label<-NULL
write.tree(tree_new,paste0(group,".tre"))


