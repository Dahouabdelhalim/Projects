library(phangorn)
library("caper")
library(adephylo)
library(phangorn)
library(stringr)


dlist <- list()

# read nexus tree file
treeFile <- file.choose()
tree <- read.nexus(treeFile)

#Sort out case if multiphylo
if(class(tree)=="multiPhylo")
{
  cat("Which tree? Enter a number from 1 to ",length(tree))

  treeNumber <- readline()
  treeNumber <- as.integer(treeNumber)

  tree<-tree[[treeNumber]]
}

#identify extinct and extant taxa (0/1)
dist_taxa <- distRoot(tree)
maxt<-max(distRoot(tree))
extant<-vector(length=(length(dist_taxa)))

#Count number of fossils - note L ensures this is an integer.
extantCount=0L;

for (j in 1:length(dist_taxa))
{
 if(dist_taxa[j]<maxt)
 {
   extant[j]<-0
 }
 else
 {
   extant[j]<-1
   extantCount<-extantCount+1
 }
}

#If only one extant taxon it's probably not a time tree
if(extantCount<2)
{
 cat("It looks like you have only one extant taxon - so this is probably not a time calibrated tree. For each terminal please enter whether the terminal is a fossil or not.\\n")
 plot(tree)
  for (i in 1:length(tree$tip.label))
   {
     cat("Is",tree$tip.label[i],"a fossil? (y/n)\\n")
     repeat
         {
           answer <- readline()
           if(answer=="y"){extant[i]<-0; break;}
           else if (answer=="n"){extant[i]<-1; break;}
           else cat ("Error try again\\n")
        }
   }

}

extant<-as.data.frame(extant)
species<-as.data.frame(tree$tip.label)
colnames(species)<-c("taxon")
chardat_this<-cbind(species,extant)

 #Estimate "D"
 tab <- as.data.frame(table(unlist(extant[])))

 if(tab[1,2]>2)
 {
   Mcomparative<-comparative.data(phy=tree, data=chardat_this, names.col=taxon, vcv=TRUE, force.root=TRUE, na.omit=TRUE, warn.dropped=FALSE)
   phyloD<-phylo.d(data=Mcomparative,binvar=extant,permut=1000)
   D<-phyloD$DEstimate
   cat("D is ",D,"\\n")

   sink("D_output.txt", append=TRUE)
   cat("Results for tree file",treeFile,"\\n\\n")
   cat("Taxon","Extant\\n")

   for (i in 1:length(tree$tip.label))
   {
     cat(tree$tip.label[i],"\\t",extant$extant[i],"\\n")
   }
   cat("D is ",D,"\\n\\n")
   sink()
 }

rm(dlist)
rm(treeFile)
rm(tree)
rm(extantCount)
rm(Mcomparative)
rm(phyloD)
rm(D)
rm(i)
rm(j)
rm(species)
rm(tab)
rm(dist_taxa)
rm(chardat_this)
rm(extant)
rm(maxt)
rm(answer)
rm(treeNumber)
