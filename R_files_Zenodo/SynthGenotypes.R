## Function for generating null distribution of genotypes using observed allele frequencies and their possible linkage.
##   The function compares the distribution of the number of supertypes in the individuals between the real population 
##   and the null distribution in the generated set of genotypes.
## Variables to be introduced by user:
##   popul - matrix of MHC genotypes with individuals in rows and alleles in columns (0/1 - allele absent/present)
##   iter - number of genotypes to be generated (positive integer, preferably large)
##   linked - a vector of linked allele pairs (using allele order of popul matrix columns)
##           e.g. linked <- c(1,2,3,6,4,10,5,14,11,12) where 1 and 2 are linked pair, 3 and 6 are linked pair, etc.
##   supertype - a vector of supertype assignment of alleles provided in the popul matrix (the order must be the same)
##  
## The function prints out:
##   Pop, mean number of supertypes in the analysed population
##   Synth (nonlinked), mean number of supertypes in the synthetic genotypes with alleles drawn at random
##   Synth (linked), mean number of supertypes in the synthetic genotypes created assuming association between some allele pairs
##   All p-values refer to the results of Kolmogorov-Smirnov test (from package dgof) of difference in distribution 
##   of the number of supertypes between the real and synthetic populations
##     

SynthGenotypes <- function(popul,iter = 99999,linked = NULL,supertype) {
  SynthNotLinked <- matrix(NA,nrow = iter,ncol = ncol(popul))
  #creating synthetic genotypes
  for (j in 1:iter){
    for (i in 1:ncol(popul)) {
      allel.i <- sample(popul[,i],1)
      SynthNotLinked[j,i]<-allel.i
    }
  }
  if (!is.null(linked)){
    SynthLinked <- matrix(NA,nrow = iter,ncol = ncol(popul))
    firstAllele <- linked[seq(1,length(linked),2)]
    secondAllele <- linked[seq(2,length(linked),2)] 
    for (j in 1:iter){
      randomLociOrder <- sample(1:ncol(popul),ncol(popul))
      ToDoList <- randomLociOrder
      for (i in randomLociOrder) {
        if (i %in% ToDoList){
          allel.i <- sample(popul[,i],1)
          SynthLinked[j,i]<-allel.i
          if (i %in% linked){
            if (i %in% firstAllele) xlink <- secondAllele[which(firstAllele==i)] else
              xlink <- firstAllele[which(secondAllele==i)]
            if (xlink %in% ToDoList){
              SynthLinked[j,xlink]<-allel.i
              ToDoList<- ToDoList[-(which(ToDoList ==xlink))]
            }
          }
          ToDoList <- ToDoList[-1]
        }
        }
    }
    }
  nSuperSynthNotLinked <- c()
  for (j in 1:iter){
    nSuperSynthNotLinked <- append(nSuperSynthNotLinked,length(unique(as.numeric(supertype[which(SynthNotLinked[j,]==1)]))))}
  if (sum(nSuperSynthNotLinked==0)!=0){
    nSuperSynthNotLinked <- nSuperSynthNotLinked[-which(nSuperSynthNotLinked==0)]}  
  
  nSuperPopul <- c()
  for (j in 1:nrow(popul)){
    nSuperPopul <- append(nSuperPopul,length(unique(as.numeric(supertype[which(popul[j,]==1)]))))}
  
  nSuperSynth <- c()
  for (j in 1:iter){nSuperSynth <- append(nSuperSynth,
                                          length(unique(as.numeric(supertype[which(SynthLinked[j,]==1)]))))}
  if (sum(nSuperSynth==0)!=0){
    nSuperSynth <- nSuperSynth[-which(nSuperSynth==0)]}  
  
  cat("Pop, ",mean(nSuperPopul),"\\n")
  cat("Synt(nonlinked), ",mean(nSuperSynthNotLinked), "p = ",dgof::ks.test(nSuperPopul,ecdf(nSuperSynthNotLinked))$p.value,"\\n")
  cat("Synt(linked), ",mean(nSuperSynth),", p = ",dgof::ks.test(nSuperPopul,ecdf(nSuperSynth))$p.value,"\\n")
  
}
