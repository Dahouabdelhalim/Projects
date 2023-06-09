library(phangorn)
library("caper")
library(adephylo)
library(phangorn)
library(stringr)

dlist <- list()
dmatrix <- matrix(, nrow = 1000, ncol = 3)
colnames(dmatrix) <- c("n_Extinct", "n_Extant", "D_est")
Pro<-"Processing treee:"

int store_i;

# Begin forloop for 999 replicates
for(i in 0:999)
{
store_i<-i;
pad<-str_pad(i, 3, pad = "0")
cat(Pro, (sprintf("%d", i)))
cat("\\n")  # read nexus tree files
readSim <- paste("TREvoSim_tree_", pad, ".nex", sep = "")
Mtree <- read.nexus(readSim) #  iterate over every sim tree, read nexus  #identify extinct and extant taxa (0/1)
dist_taxa <- distRoot(Mtree)
maxt<-max(distRoot(Mtree))
extant<-vector(length=(length(dist_taxa)))

for (j in 1:length(dist_taxa))
                         {
                             if(dist_taxa[j]<maxt)
                               {
                               extant[j]<-0
                               }
                             else
                               {
                                extant[j]<-1
                               }
                         }
extant<-as.data.frame(extant)
species<-as.data.frame(Mtree$tip.label)
colnames(species)<-c("taxon")
chardat_this<-cbind(species,extant)  #Estimate "D"
tab <- as.data.frame(table(unlist(extant[])))

if(tab[1,2]>2)
    {
    Mcomparative<-comparative.data(phy=Mtree, data=chardat_this, names.col=taxon, vcv=TRUE, force.root=TRUE, na.omit=TRUE, warn.dropped=FALSE)
    phyloD<-phylo.d(data=Mcomparative,binvar=extant,permut=1000)
    D<-phyloD$DEstimate

    #Numbering of list goes from 1 - need to add one to get out of C numbering
    dmatrix[i,3]<-D
    dmatrix[i,1]<-tab[1,2]
    dmatrix[i,2]<-tab[2,2]
    dlist[i+1] <- D
    names(dlist)[i+1]<-(sprintf("Run_%d", i+1))
    #dlist[i+1].colnames(sprintf("%d", i))
    }
else
   {
            dlist[i] <- NULL
   }
}

mean(unlist(dlist))
dlist[store_i+2] <-mean(unlist(dlist))
names(dlist)[store_i+2]<-(sprintf("Mean"))
Dmean<-mean(dmatrix[,3])
#write.table(dlist, file = "dlist.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
write.csv(dlist, file = "dlist.csv", row.names=FALSE)
write.csv(dmatrix, file = "dmatrix.csv", row.names = FALSE)

rm(dlist)
rm(Mtree)
rm(Mcomparative)
rm(phyloD)
rm(D)
rm(i)
rm(j)
rm(pad)
rm(readSim)
rm(species)
rm(tab)
rm(dist_taxa)
rm(chardat_this)
rm(extant)
rm(maxt)
