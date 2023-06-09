################################################################################
##########Using the functions in BarcodingR as a base, we wrote our own loops
##########to generate pairwise inter-and-intraspecific genetic distance 
##########calculations.  If using this script, please cite(BarcodingR) as well!
################################################################################

library(BarcodingR)
#This R script assumes that a fasta alignment of barcodes has been imported to R
#and named as the object ref . We used read.dna for this purpose.

#extract the species names from the imported alignment
sampleSpeNames <- attr(ref, "dimnames")[[1]]
#We'll use gsub to find the unique species name references.
mpattern <- ".+,"
Spp <- gsub(mpattern, "", sampleSpeNames)
#we are now going to convert this list of species names into a factor.
f <- factor(Spp)
#Let's extract the different species names
levels(f)->tax
#Now we generate every pairwise combination of species comparisons.
combn(tax,m=2,simplify=FALSE) -> combos

#We now set up some empty dataframes to loop over.
inters <- data.frame(NULL)
intras <- data.frame(NULL)
total <- data.frame(NULL)



#calculate interspecific values - this is based entirely on the barcoding.gap
#function in R.

for (i in 1:length(combos)){
  combos[[i]][[1]] -> temp1
  combos[[i]][[2]] -> temp2
  name_to_subset1 <-paste0('\\\\b',temp1,'\\\\b')
  name_to_subset2 <-paste0('\\\\b',temp2,'\\\\b')
  ref_adults[grep(name_to_subset1,rownames(ref_adults)),] -> sp1
  ref_adults[grep(name_to_subset2,rownames(ref_adults)),] -> sp2
  n1 <- dim(sp1)[1]
  n2 <- dim(sp2)[1]
  sp12 <- rbind(sp1, sp2)
  dist <- dist.dna(sp12, model = "K80")
  dist <- as.matrix(dist)
  diag(dist) <- NA
  inter <- dist[(n1 + 1):(n1 + n2), 1:n1]
  inter <- list(inter)
  inter <- unlist(inter)
  mean(inter) -> mean_temp
  data.frame(temp1,temp2,mean_temp) -> inter_temp
  colnames(inter_temp) <- c("sp1","sp2","mean")
  rbind(inters,inter_temp)->inters
}


#Now let's calculate the intraspecific values as well. Again, based on barcoding.gap
total_intra <-data.frame(NULL)
summary_intra <- data.frame(NULL)
for (i in 1:length(tax)){
  tax[i] -> temp1
  name_to_subset1 <-paste0('\\\\b',temp1,'\\\\b')
  ref_adults[grep(name_to_subset1,rownames(ref_adults)),] -> sp1
  n1 <- dim(sp1)[1]
  intra1 <- dist.dna(sp1, model = "k80", as.matrix = TRUE)
  diag(intra1) <- NA
  #cbind(melt(intra1),temp1) -> test2
  intra1[upper.tri(intra1)] <- NA
  cbind(melt(intra1),temp1)->temp_intra
  colnames(temp_intra) <- c("ind1","ind2","dist","sp")
  length(levels(temp_intra$ind1)) -> N_temp
  mean(temp_intra$dist,na.rm=TRUE) -> mean_temp
  min(temp_intra$dist,na.rm=TRUE) -> min_temp
  max(temp_intra$dist,na.rm=TRUE) -> max_temp
  if (N_temp==1) {
    mean_temp <- NA
    min_temp <- NA
    max_temp <- NA
  }
  data.frame(temp1,N_temp,mean_temp,min_temp,max_temp) -> summary_temp
  colnames(summary_temp) <- c("sp","N","intra_mean","intra_min","intra_max")
  rbind(summary_intra,summary_temp) -> summary_intra
  na.exclude(temp_intra) -> temp_intra
  if (N_temp==1) {
    cbind(melt(intra1),temp1)->temp_intra
    colnames(temp_intra) <- c("ind1","ind2","dist","sp")
  }
  rbind(total_intra,temp_intra) -> total_intra
}


