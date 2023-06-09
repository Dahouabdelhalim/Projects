library(phangorn)
library(caper)
library(adephylo)
library(phangorn)
library(stringr)
library(apTreeshape)
​
clist <- list()
cmatrix <- matrix(, nrow = 1000, ncol = 1)
colnames(cmatrix) <- c( "Treeshape")
Pro<-"Processing tree:"
store_i = 0L
​
# Begin forloop for 999 replicates
for(i in 0:999)
{
  store_i<-i;
  pad<-str_pad(i, 3, pad = "0")
  cat(Pro, (sprintf("%d", i)))
  cat("\\n")  # read nexus tree files
  readSim <- paste("TREvoSim_tree_", pad, ".nex", sep = "")
  Mtree <- read.nexus(readSim) #  iterate over every sim tree, read nexus and calculate treeshape (colless index)
  Tshape <- as.treeshape(Mtree)
  Ttips <- Ntip(Mtree)
  shapevalue <- colless(Tshape)
  maxvalue <- ((Ttips-1)*(Ttips-2))/2
  N.shapevalue <- shapevalue/maxvalue
  cmatrix[i,]<- N.shapevalue
  cat("Shape is, ",sprintf("%g", N.shapevalue),"\\n")
}
​
write.table(clist, file = "clist,n.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
write.csv(clist, file = "clist.n.csv", row.names=FALSE)
write.csv(cmatrix, file = "cmatrix.n.csv", row.names = FALSE)
Collapse
