####
######## Working Directory and Libraries ###########
####

#set WD to a folder containing the allele frequencies (a csv called allele_freqs.csv), and include the datfile_part1.txt and datfile_part2.txt files if you want to specify the COLONY parameters and then have the script generate a formatted .dat file to be run directly in COLONY.  
setwd("/Users/kloope/Dropbox/Projects/Madison_Vpenn/simulating_data_for_colony/Supp Info/")

#Libraries
library(tidyr)
library(ggplot2)

####
######## Functions ################################
####

#***************************************************
#SAMPLE_ALLELE: make a sampler for any number of loci with any number of alleles, input is a data frame with columns titled "marker", "allele","freq"
sample_allele<-function(df){
  markers<-df[,c("marker")]
  alleles<-df[,c("allele")]
  freq<-df[,c("freq")]
  umarkers<-unique(markers) #make a list of the markers
  outvector<-vector(length=length(umarkers)) #make an empty output vector of the right length
  for(i in c(1:length(umarkers))){ #for each of the markers...
    sub<-markers==umarkers[i] #these are the subset of the rows that are for the current marker
    outvector[i]<-sample(x = alleles[sub], 1, replace = T, prob = freq[sub])  #sample the alleles from the current marker with probability according to the freqs of current marker (with replacement)
  }
  outvector #return the vector
}


#___________________________________________________



#***************************************************
#NEWQUEENS: function that generates n random queens in a list from a specified allele frequency df
newqueens <- function(n,allele_freq){
  outlist<- vector(mode = "list", length = n)
  for(i in c(1:n)){
    outlist[[i]]<-as.matrix(rbind(sample_allele(allele_freq),sample_allele(allele_freq)))
  }
  outlist
}
#___________________________________________________



#***************************************************
#NEWMALES: function that generates n random males in a list from a specified alllele frequency df
newmales <- function(n,allele_freq){
  outlist<- vector(mode = "list", length = n)
  for(i in c(1:n)){
    outlist[[i]]<-as.matrix(sample_allele(allele_freq))
  }
  outlist
}
#___________________________________________________


#***************************************************
#SAMPLE_QUEEN: function to sample randomly from a queen's alleles to generate offspring types
sample_queen<-function(queen){
  apply(queen, 2, function(x) x[sample(c(1,2),1)]) #randomly choose element 1 or 2 from a column, then apply across columns for the matrix queen
}
#___________________________________________________



#***************************************************
#MAKE_DAUGHTER: make a daughter from a queen and a male
make_daughter<-function(queen,male){
  cbind(as.matrix(sample_queen(queen)),male) #makes a daughter (2 column format!)
}
#___________________________________________________


#***************************************************
#MAKE_HALF_SIBS: This function generates n half-siblings from a random mother given an allele_freq df.  use this to generate sib queens (output of females is in 2 column format)
make_half_sibs<-function(n,allele_freq){
  mother<-newqueens(1,allele_freq)[[1]] #generate a random mother from the allele freqs
  siblist<- vector(mode = "list", length = n)
  for (i in c(1:n)){
    siblist[[i]]<-t(make_daughter(mother,sample_allele(allele_freq)))
  }
  siblist
}
#___________________________________________________



#***************************************************
#MAKE_FATHERS: make fathers for a given list of queens.  this will make a list of lists, each "outer" list corresponds to a queen, and the inner list is a list of each of her male mates.  requires an allele freq df as well
make_fathers<-function(queenlist,nmates=5,allele_freq){
  l<-length(queenlist)
  outlist<-vector(mode="list",length=l) #make a list of lists...one for each female
  for (i in c(1:l)){
    outlist[[i]]<-newmales(nmates,allele_freq)
  }
  outlist
}
#___________________________________________________


#***************************************************
#MAKE_COLONY: here's a function to make a colony of n daughters from a randomly sampled queen from a pregenerated list of queens (using newqueens() or make_half_sibs()) and a randomly sampled male mate from a pre-generated list of male mates (function make_fathers())
make_colony<-function(ndaughter,nson,queens,mates,colnum){
  outlist<-vector(mode="list",length=ndaughter+nson)
  if(ndaughter>0){ #only make daughters if ndaughter > 0
    for(i in c(1:ndaughter)){
      q<-sample(c(1:length(queens)),1) #sample a queen
      m<-sample(c(1:length(mates[[q]])),1) #sample a mate of that queen
      outlist[[i]]<-c(paste0("C",colnum),paste0("daughter",sprintf("%03d", i)),paste0("Q",q),paste0("M",m),t(make_daughter(queens[[q]],mates[[q]][[m]])))#make a daughter from those parents and add the colony, queen and male info at the start.  transpose the makedaughter output so it makes one long list of alleles with pairs from each marker next to each other in order of markers.
    }
  }
  if(nson>0){ #only make sons if nson > 0
    for(i in c(1:nson)){
      q<-sample(c(1:length(queens)),1) #sample a queen
      m<-"NA" #no mate..it's a son
      qgeno<-sample_queen(queens[[q]])
      outlist[[ndaughter+i]]<-c(paste0("C",colnum),paste0("son",sprintf("%03d", i)),paste0("Q",q),paste0("M",m),c(rbind(qgeno,rep("-99",length(qgeno))))) #make a son from one of the queensa nd add the colony, queen and male info at the start.  make second allele "-99" for each marker, so it's identified as haploid by colony program)
    }
  }
  outlist
  
}
#___________________________________________________



#***************************************************
#MAKE_COLONIES: this makes multiple colonies, each with: noff= # of offpsring, queens = # of queens, mates = # mates per queen, ncol = number of colonies to make
make_colonies<-function(ndaughter,nson,queens,mates,ncol){
  outlist<-vector(mode="list",length=ncol)
  for(i in c(1:ncol)) { 
    outlist[[i]]<-make_colony(ndaughter,nson,queens,mates,i) #makes the colonies
    outlist[[i]]<-as.data.frame(do.call(rbind,outlist[[i]],)) #combines the lists into data frames, one for each colony
    outlist[[i]][]<-lapply(outlist[[i]], as.character) # converts factors to characters
  }
  outlist<-do.call("rbind", outlist) #rbinds all of the colony data frames together into one data frame
  nmarker<-(ncol(outlist)-4)/2
  rnamesmarkers<-c()  
  for(i in 1:nmarker){
    a1<-paste0("mk",i,"-1")
    a2<-paste0("mk",i,"-2")
    rnamesmarkers<-c(rnamesmarkers,a1,a2)
  }
  names(outlist)<-c(c("colony","offspringID","queen","father"),rnamesmarkers)
  outlist<-within(outlist, id <- paste(colony,offspringID,queen,father,sep='-'))
  x<-ncol(outlist)
  outlist<-outlist[,c(x,1:(x-1))]
  outlist
}
#___________________________________________________



#***************************************************
#FORMAT_QUEENS: takes a queen list and turns it into a data frame ready for export
format_queens<-function(queens){
  out<-lapply(queens,"as.vector")#make them into single, interleaved genotype vectors
  out<-do.call("rbind",out)# bind them into one matrix
  out<-cbind(paste0("Q",c(1:length(queens))),out) #add the numbers to the first column
  out<-data.frame(out,stringsAsFactors = F) #turn them into a data frame
  nmarker<-(ncol(out)-1)/2
  rnamesmarkers<-c()  
  for(i in 1:nmarker){
    a1<-paste0("mk",i,"-1")
    a2<-paste0("mk",i,"-2")
    rnamesmarkers<-c(rnamesmarkers,a1,a2)
  }
  names(out)<-c("qID",rnamesmarkers)
  out
}
#___________________________________________________



#***************************************************
#FORMAT_FATHERS: 
format_fathers<-function(fathers){
  out<-vector(mode="list")
  for(i in c(1:length(fathers))){
    out[[i]]<-do.call("rbind",fathers[[i]]) #collapse each of the lists into a matrix
  }
  nfathers<-length(fathers[[1]])
  nqueens<-length(fathers)
  #make two lists of names, one of queens, one of fathers
  qID<-vector()
  mID<-vector()
  for(i in c(1:length(fathers))){
    for(j in c(1:length(fathers[[i]]))){
      qID<-c(qID,paste0("Q",i)) 
      mID<-c(mID,paste0("M",j))
    }
  }
  out<-do.call("rbind",out)
  #######PROBLEM BELOW HERE!!!####
  #out<-cbind(qID,mID,) # collapse the list into a bigger matrix and bind the queen names and father names
  
  out<-cbind(qID,mID,out)
  out<-data.frame(out,stringsAsFactors = F) #turn them into a data frame
  names(out)<-c("qID","mID",paste0("mk",c(1:(ncol(out)-2))))
  out
  #qID
  #mID
}
#___________________________________________________



#***************************************************
#EXCLUSIONSLISTER: Here's a function that ouptuts appropriate exlusions list from two vectors, one of individuals, the other of their colonies.  exclusions are maternal sibship exclusions listing all individuals in the analysis that are from other colonies (and thus cannot be siblings) of each individual
exclusionslister<-function(individuals,colonies){
  out<-list()
  temp1<-data.frame(Individual=individuals,colony=colonies)
  for (i in c(1:length(individuals))){
    #loop through the individuals...for each individual, make a vector of the names of all the individuals not in the same colony
    temp<-as.vector(temp1[temp1$colony!=temp1$colony[i],]$Individual)
    #to the ith element of the list, add the individuals name, the number of excluded individuals, and the vector of excluded individual names
    out[[i]]<-c(toString(temp1$Individual[i]),length(temp),temp)
  }
  out
}
#___________________________________________________



#***************************************************
#ERROR_MISSING: Here's a function to apply an error rate and drop samples to create missing data to match what normal data would look like.  input is a df from make_colonies function
error_missing<-function(error,missing,df){
  df_front<-df[,c(1:5)] #save the first 5 columns for appending later
  df<-df[,-c(1:5)] #just get the genotype data
  for(i in c(1:nrow(df))){ #first loop through every allele to see if it should be changed, and change it.
    for(j in c(1:ncol(df))){
      change<-sample(c(T,F),1,replace=F,prob=c(error,1-error))#should the genotype change?
      if(change==T){ #if it should be changed...
        if(df[i,j]==1) { #..and the genotype is 1...
          df[i,j]<-2 #set it to 2.
        } else { #othrwise...
          if(df[i,j]==2) {df[i,j]<-1}} #..if it's equal to 2, set it to 1.  Slow and ugly, but easy to write.
      }
    }
  }
  
  #now loop through each locus to see if it should be missing.  
  for(i in c(1:nrow(df))){
    for(j in c(1:(ncol(df)/2))){
      ismissing<-sample(c(T,F),1,replace=F,prob=c(missing,1-missing))#should the (biallilic) genotype be missing?
      if(ismissing==T) {
        df[i,2*j-1]<-0 #set the first genotype to zero
        if(df[i,2*j]!=-99) { df[i,2*j] <-0} #also set the other genotype from that pair to zero if it's not a male
      }
    }
  }
  
  cbind(df_front,df)
}



######### IMPORT ALLELE FREQS ##########
#import allele frequency data
allele_freq<-read.csv("allele_freqs.csv",stringsAsFactors=F)

names(allele_freq)<-c("marker","allele","freq") #change names so the columns are recongized by functions later
allele_freq$allele<-as.character(allele_freq$allele) #make the alleles characters instead of integers



########## SIMULATE COLONIES ###########

#first, import the standard .dat settings for the colony.dat file.  the lines that will be specified by the function are:
#dat1[1]<-name1
#dat1[2]<-name2
#dat1[3]<-number of offspring
#dat1[4]<-number of loci
#dat1[23]<-marker prefix (default to "mk")
#dat1[25]<-error rate1
#dat1[26]<-error rate2

#and for the second part:
#dat2[19] <-number of offspring with excluded maternal sibships

fileConnection<-file("/Users/kloope/Dropbox/Projects/Madison_Vpenn/simulating_data_for_colony/datfile_part1.txt")
dat1<-readLines(fileConnection)
close(fileConnection)

fileConnection<-file("/Users/kloope/Dropbox/Projects/Madison_Vpenn/simulating_data_for_colony/datfile_part2.txt")
dat2<-readLines(fileConnection)
close(fileConnection)


#I'm going to write a function here that generates all of the output files produced line by line below.  It'll have a lot of variables to specify.

export_all<-function(allele_freq,nq,related,nmates,ndaughter,nson,ncol,error,missing,mkname="mk",likelihood_mode="1",run_length="3",dat1file=dat1,dat2file=dat2,seed="1234") {
  if(related==F){queens<-newqueens(nq,allele_freq)}
  if(related==T){queens<-make_half_sibs(nq,allele_freq)}
  fathers<-make_fathers(queens,nmates,allele_freq)
  colonies<-make_colonies(ndaughter,nson,queens,fathers,8)
  colonies_em<-error_missing(error,missing,colonies)
  exclusions<-exclusionslister(colonies$id,colonies$colony)
  if(related==F)(rel<-"u")
  if(related==T)(rel<-"r")
  
  run_name<-paste0(rel,"_q",nq,"_m",nmates,"_d",ndaughter,"_s",nson,"_ncol",ncol,"_e",error,"_m",missing)
  
  q_name<-paste0("qtypes-",run_name,".csv")
  m_name<-paste0("mtypes-",run_name,".csv")
  c_name<-paste0("colonies-",run_name,".txt")
  em_name<-paste0("em_colonies-",run_name,".txt")
  ex_name<-paste0("excl-",run_name,".txt")
  datfile_name<-paste0("datfile-",run_name,".dat")
  
  write.table(queens,file=q_name,sep=",",quote=F,row.names=F,col.names=T,append=F)
  write.table(fathers,file=m_name,sep=",",quote=F,row.names=F,col.names=T,append=F)
  write.table(colonies,file=c_name,sep="\\t",quote=F,row.names=F,col.names=T,append=F)
  write.table(colonies_em,file=em_name,sep="\\t",quote=F,row.names=F,col.names=T,append=F)
  

  fileConn<-file(ex_name)
  writeLines(unlist(lapply(exclusions, paste, collapse="\\t")),fileConn)
  close(fileConn)
  
  #now modify the dat files and export the final datfile
  
  dat1file[1]<-paste0("'",run_name,"'") #name
  dat1file[2]<-paste0("'",run_name,"'") #name
  dat1file[3]<-paste0(nrow(colonies),dat1file[3]) #number of offspring
  dat1file[4]<-paste0((ncol(colonies)-5)/2,dat1file[4]) #number of loci
  dat1file[5]<-paste0(seed,dat1[5]) #seed
  dat1file[16]<-paste0(run_length,dat1[16])
  dat1file[20]<-paste0(likelihood_mode,dat1file[20])
  dat1file[23]<-paste0(mkname,dat1file[23]) # (default to "mk")
  dat1file[25]<-paste0(error,dat1file[25]) #error rate 1
  dat1file[26]<-paste0(error,dat1file[26]) #error rate 2
  
  dat2file[19]<-length(exclusions) #name
  
  #ok, now create the .dat file
  fileConn<-file(datfile_name,"w")
  writeLines(dat1file,fileConn) #write the first part of the dat file
  write.table(colonies_em[,-(2:5)],file=fileConn,sep="\\t",quote=F,row.names=F,col.names=F,append=T) #add the genotypes
  writeLines(dat2file,fileConn) #write the second part
  writeLines(unlist(lapply(exclusions, paste, collapse="\\t")),fileConn) #write the exclusions
  writeLines("",fileConn)
  close(fileConn)
}


#now we can generate a lot of datasets!
#here are the unrelated queen datasets
export_all(allele_freq,nq=1,related=F,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=2,related=F,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=5,related=F,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=10,related=F,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=15,related=F,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=20,related=F,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=25,related=F,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=30,related=F,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=40,related=F,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=50,related=F,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)

#here are the related queen datasets
export_all(allele_freq,nq=2,related=T,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=5,related=T,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=10,related=T,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=15,related=T,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=20,related=T,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=25,related=T,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=30,related=T,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=40,related=T,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)
export_all(allele_freq,nq=50,related=T,nmates=5,ndaughter=60,nson=40,ncol=8,error=0.03,missing=0.065)

#now, run those .dat files in COLONY.  


#Here's example code to create them manually, rather than using export_all
##########case 1: 5 queens, unrelated #########
u5q<-newqueens(5,allele_freq) # make 5 unrelated queens
u5m<-make_fathers(u5q,nmates=5,allele_freq) # make 5 unrelated male mates per queen
u5col<-make_colony(ndaughter=100,nson=100,u5q,u5m,1) # make one colony as a test

#make 8 colonies
u5cols<-make_colonies(60,40,u5q,u5m,8) #make 8 colonies

#add errors and missing values.  I calculated 6.5% missing data rate for the SNPs dataset.  We're assuming 3%
u5cols_em<-error_missing(0.03,.065,u5cols)

fileConn<-file("testout.txt")
writeLines(u5cols_em[,-(2:5)],sep="\\t")
close(fileConn)


#add genotype data to outfile
write.csv(u5cols,"testout.csv")



#generate and export exclusion list (will work for any of the datasets since they have the same individuals)
u5exclusions<-exclusionslister(u5cols$id,u5cols$colony)
length(u5exclusions)
#export the list to a space delimited file, which you then have to find-replace in textwrangler to get tabs in there
sink("exclusions_test.txt")
writeLines(unlist(lapply(exclusions, paste, collapse=" ")))
sink()



#now we need to output the queen and male genotypes so we know them later
u5qout<-format_queens(u5q)
u5qout #looks good!

u5mout<-format_fathers(u5m)
u5mout

#export these, as well as the exclusions list




########## Analyzing COLONY results ################
#now we want to import the results from the COLONY run and process them to see how well they did

#import an example file (pick a best config file from COLONY, convert to csv, then import)
data<-read.csv("~/Supp_Info/r_q2_m5_d60_s40_ncol8_e0.03_m0.065-1.BestConfig.csv",header=T,stringsAsFactors = F)


#first we want to exract offspring data from the ID name

splitparents<-function(df){
  deets<-as.data.frame(do.call(rbind, strsplit(df[,c("OffspringID")],"-")))
  
  names(deets)<-c("colony","individual","truemother","truefather")
  deets$truemf<-paste0(deets$truemother,deets$truefather)
  data<-cbind(deets,df)
  data$estmf<-paste0(data$MotherID,data$FatherID)
  data
}

#now we want to make a set of matrices, one for each colony, describing the pairwise relationships between all individuals (full, half, unrelated)
truematrix<-function(data){
listresults<-list()
for(i in c(1:length(unique(data$colony)))){
  data_colony<-data[data$colony==unique(data$colony)[i],]
  
  result<-matrix(nrow=nrow(data_colony),ncol=nrow(data_colony))
  result[,]<-"unrelated"
  for(j in c(1:nrow(result))){
    result[j,data_colony$truemother[j]==data_colony$truemother]<-"half"
    result[j,data_colony$truemf[j]==data_colony$truemf]<-"full"
   
  }
  listresults[[i]]<-result
}
listresults
}

#now make a function to do the same for the colony estimates of parentage
estmatrix<-function(data){
  listresults<-list()
  for(i in c(1:length(unique(data$colony)))){
    data_colony<-data[data$colony==unique(data$colony)[i],]
    
    result<-matrix(nrow=nrow(data_colony),ncol=nrow(data_colony))
    result[,]<-"unrelated"
    for(j in c(1:nrow(result))){
      result[j,data_colony$MotherID[j]==data_colony$MotherID]<-"half"
      result[j,data_colony$estmf[j]==data_colony$estmf]<-"full"
      
    }
    listresults[[i]]<-result
  }
  listresults
}

#test those functions out on the imported COLONY output file
truematrix(splitparents(data))
estmatrix(splitparents(data))

#now we want to compare the two data sets (the true data and the estimated colony result).  We can paste them together and then summarize the result to get the transitions

true1<-truematrix(splitparents(data))
est1<-estmatrix(splitparents(data))

table(paste(true1[[2]],est1[[2]],sep=">"))

#write a function to get a table of the counts for each of the "transitions" from the true relationships to the estimated relationshpis
tabletrans<-function(matrix1,matrix2){
 x<-vector(length=9)
 matrixsum<-paste(matrix1,matrix2,sep=">")
x[1]<-sum(matrixsum=="full>full")-ncol(matrix1)
x[2]<-sum(matrixsum=="full>half")
x[3]<-sum(matrixsum=="full>unrelated")
x[4]<-sum(matrixsum=="half>full")
x[5]<-sum(matrixsum=="half>half")
x[6]<-sum(matrixsum=="half>unrelated")
x[7]<-sum(matrixsum=="unrelated>full")
x[8]<-sum(matrixsum=="unrelated>half")
x[9]<-sum(matrixsum=="unrelated>unrelated")
x/2 #divide by 2 to get the number of non-self relationships (since there are two entries per relationship in the matrix, except for the self relationship which is already subtracted out in the x[1] sum)
}

#now do tabletrans() on lists of colonies from input files, and tabluate number of queens and mates in each colony
tablescolonies<-function(df){
  splitdata<-splitparents(df)
  matrix1list<-truematrix(splitdata)
  matrix2list<-estmatrix(splitdata)
  outlist<-list()
  qtrue<-vector(length=length(matrix1list))
  mtrue<-vector(length=length(matrix1list))
  qest<-vector(length=length(matrix1list))
  mest<-vector(length=length(matrix1list))
  colonies<-unique(splitdata$colony)
  
  for(i in c(1:length(matrix1list))){
    outlist[[i]]<-tabletrans(matrix1list[[i]],matrix2list[[i]])
    
    qtrue[i]<-length(unique(splitdata[splitdata$colony==colonies[i],c("truemother")]))
    mtrue[i]<-length(unique(splitdata[splitdata$colony==colonies[i],c("truemf")]))-sum(grepl("MNA",unique(splitdata[splitdata$colony==colonies[i],c("truemf")]))) #here we have to subtract the combined mother-father names that correspond to qeuens who produced males and thus have a m-f name including MNA.  we can't just count the fathers since they have redundant names across queens
    qest[i]<-length(unique(splitdata[splitdata$colony==colonies[i],c("MotherID")]))
    mest[i]<-length(unique(splitdata[splitdata$colony==colonies[i],c("FatherID")]))-sum(grepl("*-",unique(splitdata[splitdata$colony==colonies[i],c("FatherID")]))) #here we just have to subtract the symbol generated by colony to represent hte father of a male (if there is a male in the colony)
  }
  outdf<-as.data.frame(do.call(rbind,outlist))
  colnames(outdf)<-c("full-full","full-half","full-unrelated","half-full","half-half","half-unrelated","unrelated-full","unrelated-half","unrelated-unrelated")
  cbind(outdf,qtrue,mtrue,qest,mest)

}



#### ok now import all of the data from the BestConfig.csv files (use the tocsv script i made to convert them first, or find/replace to get in csv format).  replace the "Runs_finished_..." with the folder containing the BestConfig.csv files from Colony runs.
best_config_files<-list.files("Runs_finished_3_31_2020", pattern="*BestConfig.csv", full.names=F)

bcs<-lapply(paste0("Runs_finished_3_31_2020/",best_config_files), read.csv,stringsAsFactors=F)

tableresults<-lapply(bcs,tablescolonies)

tableresults

#generate the average percentage of each error type for each colony in each run


#this will get the percentage of each type of error out of the total number of the true relationships of each type.  so it will give you the percent of true full sib relationships that were estimated to be true sibs,  half sibs, and unrelated, and same for half and for unrelated.  It will then also get whether the queens were related or not, and how many queens were in the whole colony (not just the sample)
get_percents<-function(df,name){
  sum_matrix<-matrix(nrow=nrow(df),ncol=ncol(df)) #make a matrix that totals the number of each type of relationship
  sum_matrix[,c(1:3)]<-rowSums(df[,c(1:3)])  
  sum_matrix[,c(4:6)]<-rowSums(df[,c(4:6)])
  sum_matrix[,c(7:9)]<-rowSums(df[,c(7:9)])
  output<-as.data.frame(as.matrix(df)/sum_matrix)
  output[,c("qtrue","mtrue","qest","mest")]<-df[,c("qtrue","mtrue","qest","mest")]
  output$related<-strsplit(name,"_")[[1]][1] #get the relatedness (r or u)
  output$queen<-as.numeric(substr(name,start=4,stop=5)) #get the number of queens
  output$colony<-paste0("c",c(1:nrow(output)))
  output
}

#now make a big data frame with all of the summary percentages from each run in one big df

get_all_percents<-function(tableresults,best_config_file){
summary<-get_percents(tableresults[[1]],best_config_files[1]) #do the first one
for(i in c(2:length(best_config_files))){ #loop across the others...
  newrows<-get_percents(tableresults[[i]],best_config_files[i])
  summary<-rbind(summary,newrows)#...adding runs 2-19
}
summary
}

#now get them all!
summary<-get_all_percents(tableresults = tableresults,best_config_file=best_config_file)

#subset just the related queens
summary_r<-summary[summary$related=="r",]

#and unrealted qeuens
summary_u<-summary[summary$related=="u",]




sum_long<-pivot_longer(summary,cols=c(1:9),names_to="change")
sum_long$original<-as.vector(sapply(sum_long$change,function(x) strsplit(x,"-")[[1]][1]))
sum_long$estimated<-as.vector(sapply(sum_long$change,function(x) strsplit(x,"-")[[1]][2]))


#Plot a stacked chart with faceting:


ggplot(subset(sum_long,related=="u"), aes(x = colony, y = value, fill = change)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(original ~ queen) + scale_fill_manual(values=c("black", "gray","red","gray","black","red","red","gray","black"))


#OK, now let's visualize it from the perspective of the estimated individuals.  Instead of the percentage of true full and half sibs that are correct, let's look a the percentage of the estimated relationships that are correct, and where the errors come from, since this will make it easier to inerpret the output of real colony runs.  This is the basis for the supplementary figure S1.

#first we need a different percentage:
get_percents_est<-function(df,name){
  sum_matrix<-matrix(nrow=nrow(df),ncol=ncol(df)) #make a matrix that totals the number of each type of relationship
  sum_matrix[,c(1,4,7)]<-rowSums(df[,c(1,4,7)])  
  sum_matrix[,c(2,5,8)]<-rowSums(df[,c(2,5,8)])
  sum_matrix[,c(3,6,9)]<-rowSums(df[,c(3,6,9)])
  output<-as.data.frame(as.matrix(df)/sum_matrix)
  output[,c("qtrue","mtrue","qest","mest")]<-df[,c("qtrue","mtrue","qest","mest")]
  output$related<-strsplit(name,"_")[[1]][1] #get the relatedness (r or u)
  output$queen<-as.numeric(substr(name,start=4,stop=5)) #get the number of queens
  output$colony<-paste0("c",c(1:nrow(output)))
  output
}



get_percents_est(test1,name1)
#now make a big data frame with all of the summary percentages from each run in one big df

get_all_percents_est<-function(tableresults,best_config_files){
  summary<-get_percents_est(tableresults[[1]],best_config_files[1]) #do the first one
  for(i in c(2:length(best_config_files))){ #loop across the others...
    newrows<-get_percents_est(tableresults[[i]],best_config_files[i])
    summary<-rbind(summary,newrows)#...adding runs 2-19
  }
  summary
}

#now calculate them all
summary_est<-get_all_percents_est(tableresults,best_config_files)


#now split the names so we ahve a original and an estimated relationship
sum_long_est<-pivot_longer(summary_est,cols=c(1:9),names_to="change")

sum_long_est$original<-as.vector(sapply(sum_long_est$change,function(x) strsplit(x,"-")[[1]][1]))

sum_long_est$estimated<-as.vector(sapply(sum_long_est$change,function(x) strsplit(x,"-")[[1]][2]))

qlabels<- c("1 Q", "2 Qs","5 Qs","10 Qs","15 Qs","20 Qs","25 Qs","30 Qs","40 Qs", "50 Qs")
names(qlabels)<-c(1,2,5,10,15,20,25,30,40,50)
#now plot them (this is Figure S1)
#unrelated queens
ggplot(subset(sum_long_est,related=="u"), aes(x = colony, y = value, fill = original)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(estimated ~ queen, labeller=labeller(queen = qlabels)) + scale_fill_manual(values=c("black", "gray","red"),name="true relationship",labels=c("full sib", "half sib", "unrelated")) +labs(x= "colony", y = "fraction of estimated pairiwse relationships")+ theme( axis.text.x = element_blank(), axis.ticks=element_blank())+ ggtitle("unrelated queens")

#relatd queens
ggplot(subset(sum_long_est,related=="r"), aes(x = colony, y = value, fill = original)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(estimated ~ queen, labeller=labeller(queen = qlabels)) + scale_fill_manual(values=c("black", "gray","red"),name="true relationship",labels=c("full sib", "half sib", "unrelated")) +labs(x= "colony", y = "fraction of estimated pairiwse relationships") + theme( axis.text.x = element_blank(), axis.ticks=element_blank())+ ggtitle("half-sibling queens")


#ok now make a version with just sibs (half or full) vs non-sibs for figure S2
get_percents_sibs_est<-function(df,name){
  sib_matrix<-matrix(data=0,nrow=nrow(df),ncol=4) #this will be the sibship matrix with 4 columns
  
  #names(sib_matrix)<-c("sib-sib","sib-non","non-sib","non-non")
  sib_matrix[,1]<-df[,1]+df[,2]+df[,4]+df[,5]
  sib_matrix[,2]<-df[,3]+df[,6]
  sib_matrix[,3]<-df[,7]+df[,8]
  sib_matrix[,4]<-df[,9]
  sum_matrix<-matrix(nrow=nrow(df),ncol=4) #make a matrix that totals the number of each type of relationship
  sum_matrix[,c(1,3)]<-rowSums(sib_matrix[,c(1,3)])  
  sum_matrix[,c(2:4)]<-rowSums(sib_matrix[,c(2,4)])

  output<-as.data.frame(sib_matrix/sum_matrix)
  names(output)<-c("sib-sib","sib-non","non-sib","non-non")
  output<-cbind(output,df[,c("qtrue","mtrue","qest","mest")])
  output$related<-strsplit(name,"_")[[1]][1] #get the relatedness (r or u)
  output$queen<-as.numeric(substr(name,start=4,stop=5)) #get the number of queens
  output$colony<-paste0("c",c(1:nrow(output)))
  output
}


#do it for all of the results
get_all_percents_sibs_est<-function(tableresults,best_config_files){
  summary<-get_percents_sibs_est(tableresults[[1]],best_config_files[1]) #do the first one
  for(i in c(2:length(best_config_files))){ #loop across the others...
    newrows<-get_percents_sibs_est(tableresults[[i]],best_config_files[i])
    summary<-rbind(summary,newrows)#...adding runs 2-19
  }
  summary
}


summary_sibs_est<-get_all_percents_sibs_est(tableresults,best_config_files)
summary_sibs_est

#now just plot the percentage of sibling relationshps that are correct (Figure S2)
par(mar=c(5.1,5.1,1.1,2.1))
plot(summary_sibs_est[summary_sibs_est$related=="u",c("qest","sib-sib")],ylim=c(0,1), xlab="estimated queen number",ylab="accuracy of estimated sibling relationships \\n (% estimated correct)",col="red")
points(summary_sibs_est[summary_sibs_est$related=="r",c("qest","sib-sib")],ylim=c(0,1),col="blue")
legend(2, .2, legend=c("unrelated Qs", "half sibling Qs"),
       col=c("red", "blue"), pch=1,cex=1,box.lty=0)



#now look at true colony queen number as a function of estimated queen number (Figure S3)
plot(jitter(summary_sibs_est[summary_sibs_est$related=="u",c("queen")],2)~summary_sibs_est[summary_sibs_est$related=="u",c("qest")],ylim=c(0,50),xlab="estimated queen number observed",ylab="true colony queen number",col="red")
points(jitter(summary_sibs_est[summary_sibs_est$related=="r",c("queen")],2)~summary_sibs_est[summary_sibs_est$related=="r",c("qest")],ylim=c(0,50),col="blue")
abline(c(0,0),1)
legend(3, 50, legend=c("unrelated Qs", "half sibling Qs"),
       col=c("red", "blue"), pch=1,cex=1,box.lty=0)
