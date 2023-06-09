# getting R/qtl input file for PNPxNOM QTL mapping 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(here)

# LGs have been renamed to match reference genome numbering
Map<-read.table(here("PNPxNOM_7more_map_LGadj.txt"),header=T) 
map<-Map[,2:4] # remove the unnecessary first column

# read in the individual IDs
IDs<-read.table(here("PNPxNOM_7more_ids.txt"),header=T)
ids<-matrix(as.character(IDs[,1]), ncol=1, byrow=TRUE)[,1]

# read in the F2 genotypes (the JM input loc file, without the header and list of individuals)
genos<-read.table(here("PNPxNOM_7more_genos.txt"))
# minus the (a,h,b)-column
genos<-genos[,-2] 

# add the colnames, individuals' ids 
colnames(genos)<-c("Locus",ids) 

# merge the map and the genotypes by locus
Geno<-merge(map,genos,by="Locus",all.x=T,sort=F) 
GenoT<-as.data.frame(t(Geno)) # switches rows nd columns
GenoT2<-cbind(id=rownames(GenoT),GenoT) # adjust rownames

# write out
write.table(GenoT2,file="allGenotypes.csv",sep=",",col.names=F,row.names=F)

# read in pheno data (already prepped and clean (see 'PNPxNOM_pheno_data_prep.R'))
Phenodata<-read.table(here("allPhenos.csv"),header=T, sep=",")

# combine them with info file 
phenos<-merge(IDs[1], Phenodata, by="id", all=T)
# (now still includes 7 individuals to be excluded, but do later)

#check and adjust variables
str(phenos)
phenos$id<-as.character(phenos$id)
phenos$sex<-as.factor(phenos$sex)
phenos$family<-as.factor(phenos$family)

# extract the ids from GenoT (not the first three rows, they are not ids); 
# as I need them in the right order, like the genotypes file
IDgeno<-as.data.frame(GenoT2$id[4:length(GenoT2$id)])
colnames(IDgeno)<-"id" # adjust the colnames

# merge phenotypic dataframe and ids (only those that have genotypes)
phenoIDs<-merge(IDgeno,phenos,by="id",all.x=F,sort=F) 

# write out
write.table(phenoIDs,file="allPhenotypes.csv",sep=",",col.names=T,row.names=F)

# join geno and pheno and adjust input 
GenoTypes<-read.table(here("allGenotypes.csv"), sep=",", header=T)
colnames(GenoTypes)[1]<-"id"
PhenoTypes<-read.table(here("allPhenotypes.csv"), sep=",", header=T)

GenoPheno<-merge(PhenoTypes, GenoTypes, by="id",all.y=T,sort=F)
GenoPheno<-GenoPheno[c(169,170,1:168),] # re-arrange rows
GenoPheno<-GenoPheno[,c(4:43,2,3,1,44:1338)] # re-arrange columns

# now remove the 7 individuals
# "194145", "194154" , "194155", "194156", "194157", "194158", "194205"
GenoPhenoOUT<-GenoPheno[!GenoPheno$id %in% c("194145", "194154" , "194155", "194156", "194157", "194158", "194205") ,]

# write out
write.table(GenoPhenoOUT,file="GenoPheno.csv",sep=",",col.names=T,row.names=F)
# and in excel make final adjustments
# >remove contents of rows 2&3 in phenotype columns
# >replace NA and u with -;  make a,b,h capital letters (!match case, whole cells!)
# > also make tricuspid yes=1, no=0 

