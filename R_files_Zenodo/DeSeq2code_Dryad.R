# This code provides the DEseq analysis code for the Burghardt et al Molecular Ecology transcriptome paper

# if you don't have DEseq package....
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
# Note that DESeq2 does not work in the most recent versions of R

setwd("~/PATH")
library("DESeq2")

#### filtering and DEseq analysis ############ 

# read in the raw count data from Htseq
M_counts = read.table("RNAseqcounts_allconditions_head.txt", header=T, row.names =1)

# Generate seperate nodule (rhizo) and root (uninoc) datasets most analysis is conducted on the nodule data set but the root data set is used to identify genes that are expressed only in nodules
counts.rhizo<-cbind(M_counts[,grep("Med",colnames(M_counts))],M_counts[,grep("Mel",colnames(M_counts))])
counts.uninoc<-M_counts[,grep("Uni",colnames(M_counts))]

# Remove the HM101 reps that were mixed innoculums
counts.rhizo<-counts.rhizo[,-which(colnames(counts.rhizo)=="X101_Med_3")]
counts.rhizo<-counts.rhizo[,-which(colnames(counts.rhizo)=="X101_Mel_3")]

#########Calculate normalization scalars for samples from UNINOCULATED ROOTS ############

#### Generate normalized count data for downstream filtering steps for expressed genes
subset<-counts.uninoc

# Explain design of experiment to DEseq in order to extract sizeFactors from the experiment
colData = data.frame(row.names = colnames(subset), 
                     host = factor(substr(colnames(subset),start = 2, stop=4))) 


dds<-DESeqDataSetFromMatrix(countData=subset,
                            colData=colData,
                            design=~host)

dds<-estimateSizeFactors(dds)
UniSF<-sizeFactors(dds) # store size factors for Roots
counts.Uni<-counts(dds,normalize=T) ###### store normalized count data for root expression

#### Require that two or more reps of a genotype have exp levels above 5 to be considered an expressed gene in uninoculated roots. Then compile all hosts together ######
counts.Uni<-data.frame(counts.Uni)

Exp<-c(row.names(counts.Uni[rowSums(counts.Uni[,grep("101",colnames(counts.Uni))] > 5) > 1,]),
       row.names(counts.Uni[rowSums(counts.Uni[,grep("340",colnames(counts.Uni))] > 5) > 1,]),
       row.names(counts.Uni[rowSums(counts.Uni[,grep("034",colnames(counts.Uni))] > 5) > 1,]),
       row.names(counts.Uni[rowSums(counts.Uni[,grep("056",colnames(counts.Uni))] > 5) > 1,]))

ExpUni<-names(table(Exp))# List of all the genes that are expressed in at least one host

#### If desired this is code to run the model to determine differential expression across hosts for genes expressed in roots using size factors calculated based on all genes####
# This data was not used in the paper but it was used to confirm that the root samples have similar amounts of host specific expression as nodules ###

#subset down to genes that meet expression criteria
subset<-counts.uninoc[rownames(counts.uninoc) %in% ExpUni,]

# Explain design of experiment to DEseq
colData = data.frame(row.names = colnames(subset), 
                     host = factor(substr(colnames(subset),start = 2, stop=4))) 

dds<-DESeqDataSetFromMatrix(countData=subset,
                            colData=colData,
                            design=~host)
sizeFactors(dds)<-UniSF

dds<-DESeq(dds,test="LRT",reduced=~1)

### Extract the model results using Benjamin Hoochberg p value adjustment and do not allow DeSeq to perform independent filtering
res.Uni<-data.frame(results(dds ,pAdjustMethod = "BH",independentFiltering = FALSE))
write.table(res.Uni,"FullModelRootResults.txt",sep="\\t", row.names=TRUE)

# Number of genes influenced by host in the roots is comparable to the nodule samples
dim(res.Uni[res.Uni$padj<.05,])

#### RUN FULL Host x Sym  MODEL FOR NODULE SAMPLES########

#### Define the nodule samples as the subset of interest and calculated size factors based on all genes ####
subset<-counts.rhizo

# Explain design of experiment to DEseq
colData = data.frame(row.names = colnames(subset), 
                     host = factor(substr(colnames(subset),start = 2, stop=4)), 
                     sym = factor(substr(colnames(subset),start = 6, stop=8))) 

dds<-DESeqDataSetFromMatrix(countData=subset,
                            colData=colData,
                            design=~host+sym+host*sym)

dds<-estimateSizeFactors(dds)
NodSF<-sizeFactors(dds)  # Save the sizeFactors for use throughout subsequent analysis
counts.normal<-counts(dds,normalize=T) ###### generate normalized counts data for filtering steps....

#### Require that all reps of at least 1 genotype have exp levels above 5 to be considered an expressed gene for the treatment. Then compile all treaments together ######
counts.normal<-data.frame(counts.normal,median=apply(counts.normal,1,median),mean=apply(counts.normal,1,mean,na.rm=TRUE))
Exp<-c(row.names(counts.normal[rowSums(counts.normal[,grep("101_Mel",colnames(counts.normal))] > 5) > 1,]),
       row.names(counts.normal[rowSums(counts.normal[,grep("101_Med",colnames(counts.normal))] > 5) > 1,]),
       row.names(counts.normal[rowSums(counts.normal[,grep("340_Mel",colnames(counts.normal))] > 5) > 1,]),
       row.names(counts.normal[rowSums(counts.normal[,grep("340_Med",colnames(counts.normal))] > 5) > 1,]),
       row.names(counts.normal[rowSums(counts.normal[,grep("034_Mel",colnames(counts.normal))] > 5) > 1,]),
       row.names(counts.normal[rowSums(counts.normal[,grep("034_Med",colnames(counts.normal))] > 5) > 1,]),
       row.names(counts.normal[rowSums(counts.normal[,grep("056_Mel",colnames(counts.normal))] > 5) > 1,]),
       row.names(counts.normal[rowSums(counts.normal[,grep("056_Med",colnames(counts.normal))] > 5) > 1,]))

ExpNew<-names(table(Exp)) ##### Number of treatments that genes are expressed in
counts.sub<-counts.normal[rownames(counts.normal) %in% ExpNew,]

# Pull out those genes only expressed in nodules if interested
NodOnly<-ExpNew[!as.character(ExpNew) %in% ExpUni]

###### Run the full model on all the expressed genes ###### 

# Read in workhorse function for summarizing DEseq result output #### 
### This takes it aguments the output of DeSeq (DESeqob,), the list of genes tested (annotation), the different results from DEseq (contrasts), and prefix to identify each test

DESeqResults<-function(DESeqobj, annotation, contrasts, prefix){
  out<-annotation
  for(i in 1:length(prefix)){
    res <- data.frame(results(DESeqobj, contrast=contrasts[i], pAdjustMethod = "BH",independentFiltering = FALSE))
    colnames(res)<-paste(prefix[i],colnames(res),sep=".")
    dat<-res[,grep("padj", colnames(res))]
    res[,gsub("padj","sig",colnames(res)[grep("padj", colnames(res))])]<-ifelse(dat>0.05 | is.na(dat),"","sig")
    out<-cbind(out,res)
    hist(dat, main=prefix[i], breaks=100,xlim=c(0,1))
    cat(prefix[i], "n alpha<0.05 = ",length(dat[dat<0.05]),"\\n")
  }
  return(out)
}

######### Conduct likelihood ratio testing on all genes Expressed in nodules #####

# Only test genes that are on Expressed Gene list
subset<-counts.rhizo[rownames(counts.rhizo) %in% ExpNew,]

# Explain design of experiment to DEseq
colData = data.frame(row.names = colnames(subset), 
                     host = factor(substr(colnames(subset),start = 2, stop=4)), 
                     sym = factor(substr(colnames(subset),start = 6, stop=8))) 

dds<-DESeqDataSetFromMatrix(countData=subset,
                            colData=colData,
                            design=~host+sym+host*sym)

sizeFactors(dds)<-as.vector(NodSF)


####### Actually Run the full Model ############ 

# Calculate GxE significance via LRT
dds<-DESeq(dds,test="LRT",reduced=~host+sym)
res.HxS<-DESeqResults(DESeqobj=dds,
                      annotation=data.frame(rownames(dds)), 
                      contrasts=sapply(resultsNames(dds)[4], list),
                      prefix=c("hostsym"))

#Calculate treatment level significance via LRT
dds<-DESeqDataSetFromMatrix(countData=subset,
                            colData=colData,
                            design=~host+sym)
sizeFactors(dds)<-NodSF
dds<-DESeq(dds,test="LRT",reduced=~host)

res.S<-DESeqResults(DESeqobj=dds,
                    annotation=data.frame(rownames(dds)), 
                    contrasts=sapply(resultsNames(dds)[3], list),
                    prefix=c("sym"))

#Calculate genotype level significance via LRT
dds<-DESeqDataSetFromMatrix(countData=subset,
                            colData=colData,
                            design=~sym+host)

sizeFactors(dds)<-as.vector(NodSF)
dds<-DESeq(dds,test="LRT",reduced=~sym)

res.H<-DESeqResults(DESeqobj=dds,
                    annotation=data.frame(rownames(dds)), 
                    contrasts=sapply(resultsNames(dds)[3], list),
                    prefix=c("host"))

# Bind together all of the model results 
res.all<-cbind(res.H,res.S,res.HxS)
res.all<-data.frame(res.all)

# label all the significant genes
ms<-res.all[,grep(".sig",colnames(res.all))]
ms$sym.sig[ms$sym.sig=="sig"]<-"S"
ms$host.sig[ms$host.sig=="sig"]<-"H"
ms$hostsym.sig[ms$hostsym.sig=="sig"]<-"HxS"

# Combine the two sets of critera and find those genes that are H, S, H+S, or HxS

sumlist<-data.frame(host=ms$host.sig =="H",
                    sym=ms$sym.sig =="S" ,
                    hostsym=ms$hostsym.sig == "HxS")

rownames(sumlist)<-rownames(ms)

my.sig<-paste(sumlist[,1],sumlist[,2],sumlist[,3],sep="")
my.sig[my.sig=="FALSEFALSEFALSE"]<-"none"
my.sig[my.sig=="TRUEFALSEFALSE"]<-"H"
my.sig[my.sig=="FALSETRUEFALSE"]<-"S"
my.sig[my.sig=="TRUETRUEFALSE"]<-"H+S"
my.sig[my.sig=="TRUEFALSETRUE"]<-"HxS"
my.sig[my.sig=="TRUETRUETRUE"]<-"HxS"
my.sig[my.sig=="FALSETRUETRUE"]<-"HxS"
my.sig[my.sig=="FALSEFALSETRUE"]<-"HxS"

res.all$sig.cat<-my.sig
table(res.all$my.sig) # number of genes in each expression category
res.all$rownames.dds.<-NULL
res.all$rownames.dds..1<-NULL
res.all$rownames.dds..2<-NULL

write.table(res.all,"FullModelHostSymResults.txt",sep="\\t", row.names=TRUE)


############# Conduct pairwise analysis on each host contrast ###############

# Create sub datasets to analyze in all pairwise comparisons of genotypes for HxS
counts.034.056<-cbind(counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("034",colnames(counts.rhizo))],counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("056",colnames(counts.rhizo))])
counts.101.034<-cbind(counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("034",colnames(counts.rhizo))],counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("101",colnames(counts.rhizo))])
counts.034.340<-cbind(counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("034",colnames(counts.rhizo))],counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("340",colnames(counts.rhizo))])
counts.101.056<-cbind(counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("101",colnames(counts.rhizo))],counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("056",colnames(counts.rhizo))])
counts.056.340<-cbind(counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("340",colnames(counts.rhizo))],counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("056",colnames(counts.rhizo))])
counts.101.340<-cbind(counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("101",colnames(counts.rhizo))],counts.rhizo[rownames(counts.rhizo) %in% ExpNew,grep("340",colnames(counts.rhizo))])

#Generate results for all sub datasets (first load in function GetDESeqHxS.LRT below)
res.all.034.340<-GetDESeqHxS.LRT(subset=counts.034.340,NodSF,counts.sub)
res.all.034.056<-GetDESeqHxS.LRT(counts.034.056,NodSF,counts.sub)
res.all.101.034<-GetDESeqHxS.LRT(counts.101.034,NodSF,counts.sub)
res.all.101.056<-GetDESeqHxS.LRT(counts.101.056,NodSF,counts.sub)
res.all.056.340<-GetDESeqHxS.LRT(counts.056.340,NodSF,counts.sub)
res.all.101.340<-GetDESeqHxS.LRT(counts.101.340,NodSF,counts.sub)

# Summarize the total # of genes in each catagory for each comparison
Summary<-data.frame(rbind(table(res.all.101.056$sig.cat),
                          table(res.all.101.034$sig.cat),
                          table(res.all.101.340$sig.cat),
                          table(res.all.034.056$sig.cat),
                          table(res.all.056.340$sig.cat),
                          table(res.all.034.340$sig.cat)))

rownames(Summary)<-c("101v056","101v034","101v340","056v034","056v340","034v340")
print(Summary) # number of genes in each category in each comparison

write.table(res.all.034.340,file="PairwiseResults.034.340.txt",sep="\\t",row.names=FALSE)
write.table(res.all.034.056,file="PairwiseResults.034.056.txt",sep="\\t",row.names=FALSE)
write.table(res.all.101.034,file="PairwiseResults.101.034.txt",sep="\\t",row.names=FALSE)
write.table(res.all.101.056,file="PairwiseResults.101.056.txt",sep="\\t",row.names=FALSE)
write.table(res.all.034.340,file="PairwiseResults.034.340.txt",sep="\\t",row.names=FALSE)
write.table(res.all.101.340,file="PairwiseResults.101.340.txt",sep="\\t",row.names=FALSE)

# List of all Genes with HxS in any comparison
HxSlist<-table(c(rownames(res.all.034.340[res.all.034.340$sig.cat=="HxS",]),
                 rownames(res.all.034.056[res.all.034.056$sig.cat=="HxS",]),
                 rownames(res.all.101.034[res.all.101.034$sig.cat=="HxS",]),
                 rownames(res.all.101.056[res.all.101.056$sig.cat=="HxS",]),
                 rownames(res.all.056.340[res.all.056.340$sig.cat=="HxS",]),
                 rownames(res.all.101.340[res.all.101.340$sig.cat=="HxS",])))

# Narrow down to those found in more than one comparison
FocalHxSGenes<-rownames(HxSlist[HxSlist>=2])

#### FUNCTION THAT CALCUTES SIGNIFICANCE AND LOG FOLD CHANGE FOR EACH sym CONTRAST #########
GetDESeqHxS.LRT<-function(subset,NodSF,counts.sub){

  # Explain design of experiment to DEseq
  colData = data.frame(row.names = colnames(subset), 
                       geno = factor(substr(colnames(subset),start = 2, stop=4)), 
                       rhizo = factor(substr(colnames(subset),start = 6, stop=8))) 
  
  dds<-DESeqDataSetFromMatrix(countData=subset,
                              colData=colData,
                              design=~geno+rhizo+geno*rhizo)
  
  SFsub<-NodSF[names(NodSF) %in% colnames(dds)] ### Use size factors calculated from all sym contrasts for consistancy
  
  # Calculate log fold changes for each genotype and rhizobia contrast for later filtering
  Mean.Med_Geno1<-rowMeans(counts.sub[,grep(paste(levels(colData[,1])[1],levels(colData[,2])[1],sep="_"),colnames(counts.sub))])
  Mean.Med_Geno2<-rowMeans(counts.sub[,grep(paste(levels(colData[,1])[2],levels(colData[,2])[1],sep="_"),colnames(counts.sub))])
  Mean.Mel_Geno1<-rowMeans(counts.sub[,grep(paste(levels(colData[,1])[1],levels(colData[,2])[2],sep="_"),colnames(counts.sub))])
  Mean.Mel_Geno2<-rowMeans(counts.sub[,grep(paste(levels(colData[,1])[2],levels(colData[,2])[2],sep="_"),colnames(counts.sub))])
  
  sym.lfc<-log2(rowMeans(cbind(counts.sub[,grep(paste(levels(colData[,1])[1],levels(colData[,2])[1],sep="_"),colnames(counts.sub))],counts.sub[,grep(paste(levels(colData[,1])[2],levels(colData[,2])[1],sep="_"),colnames(counts.sub))]))+1)-
    log2(rowMeans(cbind(counts.sub[,grep(paste(levels(colData[,1])[1],levels(colData[,2])[2],sep="_"),colnames(counts.sub))],counts.sub[,grep(paste(levels(colData[,1])[2],levels(colData[,2])[2],sep="_"),colnames(counts.sub))]))+1)
  
  host.lfc<-log2(rowMeans(counts.sub[,grep(levels(colData[,1])[1],colnames(counts.sub))])+1)-log2(rowMeans(counts.sub[,grep(levels(colData[,1])[2],colnames(counts.sub))])+1)
  
  Med.lfc<-log2(rowMeans(counts.sub[,grep(paste(levels(colData[,1])[1],levels(colData[,2])[1],sep="_"),colnames(counts.sub))])+1)-log2(rowMeans(counts.sub[,grep(paste(levels(colData[,1])[2],levels(colData[,2])[1],sep="_"),colnames(counts.sub))])+1)
  Mel.lfc<-log2(rowMeans(counts.sub[,grep(paste(levels(colData[,1])[1],levels(colData[,2])[2],sep="_"),colnames(counts.sub))])+1)-log2(rowMeans(counts.sub[,grep(paste(levels(colData[,1])[2],levels(colData[,2])[2],sep="_"),colnames(counts.sub))])+1)
  hostsym.lfc<-abs(Med.lfc-Mel.lfc)
  
  # Calculate HxS significance via LRT
  sizeFactors(dds)<-SFsub
  dds<-DESeq(dds,test="LRT",reduced=~geno+rhizo)
  res.HxS<-DESeqResults(DESeqobj=dds,
                        annotation=data.frame(rownames(dds)), 
                        contrasts=sapply(resultsNames(dds)[4], list),
                        prefix=c("hostsym"))
  
  #Calculate treatment level significance via LRT
  
  dds<-DESeqDataSetFromMatrix(countData=subset,
                              colData=colData,
                              design=~geno+rhizo)
  sizeFactors(dds)<-SFsub
  dds<-DESeq(dds,test="LRT",reduced=~geno)
  
  res.S<-DESeqResults(DESeqobj=dds,
                      annotation=data.frame(rownames(dds)), 
                      contrasts=sapply(resultsNames(dds)[3], list),
                      prefix=c("sym"))
  
  #Calculate genotype level significance via LRT
  
  dds<-DESeqDataSetFromMatrix(countData=subset,
                              colData=colData,
                              design=~rhizo+geno)
  sizeFactors(dds)<-SFsub
  dds<-DESeq(dds,test="LRT",reduced=~rhizo)
  
  res.H<-DESeqResults(DESeqobj=dds,
                      annotation=data.frame(rownames(dds)), 
                      contrasts=sapply(resultsNames(dds)[3], list),
                      prefix=c("host"))
  
  # Bind together all of the model results 
  res.all<-cbind(res.H,res.S,res.HxS)
  res.all<-data.frame(res.all)
  res.all<-cbind(res.all,host.lfc,sym.lfc,hostsym.lfc)
  
  # label all the significant genes
  ms<-res.all[,grep(".sig",colnames(res.all))]
  ms$sym.sig[ms$sym.sig=="sig"]<-"S"
  ms$host.sig[ms$host.sig=="sig"]<-"H"
  ms$hostsym.sig[ms$hostsym.sig=="sig"]<-"HxS"
  
  # Label all the log fold changes that have an absolute value >1
  mlfc<-data.frame(host.lfc,sym.lfc,hostsym.lfc)
  mlfc$sym.lfc[abs(mlfc$sym.lfc)>=1&!is.na(mlfc$sym.lfc)]<-"S"
  mlfc$sym.lfc[abs(as.numeric(mlfc$sym.lfc))<1|is.na(mlfc$sym.lfc)]<-""
  mlfc$host.lfc[abs(mlfc$host.lfc)>=1&!is.na(mlfc$host.lfc)]<-"H"
  mlfc$host.lfc[abs(as.numeric(mlfc$host.lfc))<1|is.na(mlfc$host.lfc)]<-""
  mlfc$hostsym.lfc[abs(mlfc$hostsym.lfc)>=1&!is.na(mlfc$hostsym.lfc)]<-"HxS"
  mlfc$hostsym.lfc[abs(as.numeric(mlfc$hostsym.lfc))<1|is.na(mlfc$hostsym.lfc)]<-""
  
  # Combine the two sets of critera and find those genes that are both sig and meet the log fold change requirements
  mylist<-cbind(ms,mlfc)
  sumlist<-data.frame(host=mylist$host.sig =="H"& mylist$host.lfc=="H",
                      sym=mylist$sym.sig =="S" & mylist$sym.lfc=="S" ,
                      hostsym=mylist$hostsym.sig == "HxS" & mylist$hostsym.lfc== "HxS")
  rownames(sumlist)<-rownames(mylist)
  
  my.sig<-paste(sumlist[,1],sumlist[,2],sumlist[,3],sep="")
  my.sig[my.sig=="FALSEFALSEFALSE"]<-"none"
  my.sig[my.sig=="TRUEFALSEFALSE"]<-"H"
  my.sig[my.sig=="FALSETRUEFALSE"]<-"S"
  my.sig[my.sig=="TRUETRUEFALSE"]<-"H+S"
  my.sig[my.sig=="TRUEFALSETRUE"]<-"HxS"
  my.sig[my.sig=="TRUETRUETRUE"]<-"HxS"
  my.sig[my.sig=="FALSETRUETRUE"]<-"HxS"
  my.sig[my.sig=="FALSEFALSETRUE"]<-"HxS"
  
  res.all$sig.cat<-my.sig
  res.all<-cbind(res.all,counts.sub,Mean.Med_Geno1,Mean.Med_Geno2,Mean.Mel_Geno1,Mean.Mel_Geno2)
  
  ### Indicate lowly expressed genes
  res.all$LowExp<-"No"
  res.all[res.all$Mean.Med_Geno1<5&res.all$Mean.Mel_Geno1<5,]$LowExp<-"Yes"
  res.all[res.all$Mean.Med_Geno2<5&res.all$Mean.Mel_Geno2<5,]$LowExp<-"Yes"
  return(res.all)
}
