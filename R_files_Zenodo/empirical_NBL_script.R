#A script for calculating the D-value (Fritz & Purvis) for empirical phylogenies. The script checks for multiphylo, asks user to select which tree if so, and then:
#1. Identifies fossils based on branch lengths
#2. If only one extant taxon is present it assumes a non-time calibrated tree and then loops through the tips asking whether they are fossils or not
#3. Calculates the D value of the fossils on the tree.
#Written by Rob Sansom & Russell Garwood

library(phangorn)
library("caper")
library(adephylo)
library(phangorn)
library(stringr)


remove_data <- function()
{
  rm(dlist, pos=1)
  rm(treeFile, pos=1)
  rm(tree, pos=1)
  rm(Mcomparative, pos=1)
  rm(phyloD, pos=1)
  rm(D, pos=1)
  rm(i, pos=1)
  rm(species, pos=1)
  rm(tab, pos=1)
  rm(dist_taxa, pos=1)
  rm(chardat_this, pos=1)
  rm(extant, pos=1)
  rm(answer, pos=1)
}

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


#Set all branch lengths to 1
tree <- compute.brlen(tree, 1)

#Create extant/extinct vector
dist_taxa <- distRoot(tree)
extant<-vector(length=(length(dist_taxa)))

#Enter whether fossil or not
cat("All branch lengths have been converted to a length of one. Now, for each terminal please enter whether the terminal is a fossil or not.\\n")
plot(tree)
for (i in 1:length(tree$tip.label))
 {
   cat("Is",tree$tip.label[i],"a fossil? (y/n)\\n")
   repeat
       {
         answer <- readline()
         if(answer=="y"){extant[i]<-0; break;}
         else if (answer=="n"){extant[i]<-1; break;}
         else if (answer=="q")
         {
           remove_data()
           stop(quitting)
         }
         else cat ("Error try again\\n")
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
   remove_data()
  }
