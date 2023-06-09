#########################   EXAMPLE FOR FORMATING DATA FOR THEIR USE IN THE
#########################   RECONSTRUCTION METHOD - CASE OF INFLUENZA DATA

## load package used for reading fasta files
library(seqinr)


## function for reading sequence data and formatting them in numeric format
##     (read.fasta() from package seqinr is embedded in this function)
read.spos=function(name){
	## address of the file (i.e. "directory/filename" or only "filename" if the file is in the working directory)
    data.spos = read.fasta(name,forceDNAtolower=TRUE,strip.desc=TRUE)
    n.spos=length(data.spos)
    ID=NULL
    gen=NULL
    freq=NULL
    for(i in 1:n.spos){
        temp0=as.character(data.spos[[i]])
        temp0[tolower(temp0)=="a"]="1"
        temp0[tolower(temp0)=="c"]="2"
        temp0[tolower(temp0)=="g"]="3"
        temp0[tolower(temp0)=="t"]="4"
        temp0[! (tolower(temp0)=="1" | tolower(temp0)=="2" | tolower(temp0)=="3" |
        tolower(temp0)=="4")]="5"
        gen=rbind(gen,as.integer(temp0))
        info=strsplit(attr(data.spos[[i]],"name"),"_")[[1]]
        ID=c(ID,info[1])
        freq=c(freq,info[3])
    }
    gen=cbind(gen[,apply(gen,2,function(u) sum(u==5)==0)])
    return(list(gen=rbind(gen),ID=ID,freq=as.numeric(freq)))
}

## function aggregating similar sequences in a set of sequences and providing the weight 
##    of each unique sequence
aggregate.spos=function(spos){
	## set of sequences formatted as the output of read.spos()
	if(nrow(spos$gen)>1){
		tobeaggregated=duplicated(spos$gen)
		M=ncol(spos$gen)
		gen1=rbind(spos$gen[!tobeaggregated,])
		weight1=spos$weight[!tobeaggregated]
		if(!is.null(spos$IDseq)){
			IDseq1=spos$IDseq[!tobeaggregated]
		}
		for(i in (1:nrow(spos$gen))[tobeaggregated]){
			j=1
			while(sum(gen1[j,]==spos$gen[i,])<M){
				j=j+1
			}
			weight1[j]=weight1[j]+spos$weight[i]
		}
		if(is.null(spos$IDseq)){
			return(list(gen=rbind(gen1),weight=weight1))
		} else {
			return(list(gen=rbind(gen1),IDseq=IDseq1,weight=weight1))
		}
	} else {
		return(spos)
	}
}

## basic host information
HOSTS=rbind(cbind(113,2:4,0,1),cbind(115,2:4,0,2),
			cbind(104,4:6,1,1),cbind(116,5:6,1,2),
			cbind(109,7:8,2,1),cbind(111,7:8,2,2),
			cbind(105,9,3,1),cbind(108,9:10,3,2),
			cbind(106,15,4,1),cbind(112,c(12,15),4,2))
colnames(HOSTS)=c("ID","time","generation","num")


## build host.table
HOSTS.G=NULL
seq.spos=NULL
for(i in 1:nrow(HOSTS)){
	pig0=HOSTS[i,1]
	day0=HOSTS[i,2]
	spos0=read.spos(paste("influenza-raw-genomic-data/",pig0,"D",day0,".fasta",sep=""))
	spos0$freq=rep(1,length(spos0$freq))
	names(spos0)=c("gen","IDseq","weight")
	seq.spos[[i]]=aggregate.spos(spos0)
	HOSTS.G=rbind(HOSTS.G,c(HOSTS[i,1:2],dim(spos0$gen)))
}
colnames(HOSTS.G)[3:4]=c("N0","M")
names(seq.spos)=paste(HOSTS[,1],HOSTS[,2],sep="-")
HOSTS.G=as.data.frame(HOSTS.G)
HOSTS.G$ID=paste(HOSTS.G$ID,HOSTS.G$time,sep="_")

## build the R list containing data required to reconstruct epidemiological links 
swine=list(readme=paste("This R list is made of the following items: readme contains the present description; host.table is a dataframe containing metadata about pigs in the transmission chain: ID (pig identifier), time (day of data collection), N0 (sequencing depth), M (length of sequence fragments); set.of.sequences is a list with 1 item per pig x time, each item is a list with three items: gen is a matrix providing the haplotypes that have been collected from the pig at the sampling time (1->A, 2->C, 3->G, 4->T, 5->other), IDseq gives sequence identifiers, weight is a vector providing the number of repetitions of each haplotype.",sep=""),
host.table=HOSTS.G,set.of.sequences=seq.spos)

## save swine as an rds file
#saveRDS(swine,file="filename.rds")

## load swine contained in the rds file
#swine=readRDS(file="filename.rds")

