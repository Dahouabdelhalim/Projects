#R-script for estimating the pedigree based inbreeding coefficient using the SNP pedigree
#Alina Niskanen, original script made in 2018, modified May 2020
#alina.niskanen@gmail.com


library(dplyr)
library(pedigree)

#Bring in the pedigree
pedigree_temp <- read.table("pedigree.txt", header = T, stringsAsFactors = F) 


#Remove dummy individuals from the pedigree
#Remove dummies from the id column
individuals <- pedigree_temp$id[! pedigree_temp$id %in% grep("^F", pedigree_temp$id, value = TRUE)]
individuals <- individuals[! individuals %in% grep("^M", individuals, value = TRUE)]
pedigree_temp <- pedigree_temp[pedigree_temp$id %in% individuals,]
#Remove dummy individuals from the parental columns
pedigree_temp$dam[which(pedigree_temp$dam %in% grep("^F", pedigree_temp$dam, value = TRUE))] <- NA
pedigree_temp$sire[which(pedigree_temp$sire %in% grep("^M", pedigree_temp$sire, value = TRUE))] <- NA


#Calculate F for the whole pedigree
ped <- pedigree_temp

# First order the pedigree and then calculate F
ped <- ped[order(orderPed(ped)),]
F_ped <- calcInbreeding(ped)
pedWF <- cbind(ped,F_ped)


#Include only the individuals with at least 2 full generations in the pedigree
##Count grandparents
x <- pedigree_temp
x$dam[is.na(x$dam)] <- "0"
x$sire[is.na(x$sire)] <- "0"

x["grandparent_count"] <- NA

i<-1
for (i in 1:nrow(x)){
  if(x[i,2]==0 & x[i,3]==0){
    x[i,4] <- 0
  } else if(x[i,2]==0 & x[i,3]!=0){
    x[i,4] <- 2 - sum(x[x$id==x[i,3],2:3]==0)
  } else if(x[i,2]!=0 & x[i,3]==0){
    x[i,4] <- 2 - sum(x[x$id==x[i,2],2:3]==0)
  } else{
    x[i,4] <- 4 - (sum(x[x$id==x[i,2],2:3]==0, na.rm=TRUE)+sum(x[x$id==x[i,3],2:3]==0, na.rm=TRUE))
  }
}


#Make a separate column of pedigree F for those individuals that have 4 grandparents in the pedigree (Fped_2_gen)
x <- full_join(x, pedWF[,c(1,4)], by="id")
x <- x[,c(1:5,5)]
colnames(x)[6]<- "Fped_2_gen"
x$Fped_2_gen[which(x$grandparent_count != 4)] <- NA

write.table(x[,c(1,6)], file = "pedigree_F_coefficients.txt", sep=" ", quote=F, col.names = T, row.names = F)

