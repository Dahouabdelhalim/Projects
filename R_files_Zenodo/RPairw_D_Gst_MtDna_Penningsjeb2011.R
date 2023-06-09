#This R script basically calculates the statistics for fig 4 of Pennings et al 2011 (JEB)
#There is one difference: in my paper I calculated the statistics using all populations, but in this script I calculate all pairwise JostD values (and pvalues) 
#the first 2 characters of the name of the individual (in the fastafile) are its location, that is why I can use substring(getName.SeqFastadna(filefas),1,2) to get the population of an individual. 
NR = 20 
# I ran the calculations 20 time, i.e., I choose a random stretch of 100 bp to calculate D and calculated a pvalue for it, and repeated this procedure 20 times, do fewer to try out the program. 
NUMRAN = 1000
# do 1000 permutations to calculate p-values, or fewer to try out the program. 

#load necessary libraries to read fasta files
library(seqinr)    
library(ape)    
#read fasta file
filefas_all<-read.fasta(file="TlongSequences.fas",set.attributes=FALSE,as.string=TRUE)  
#store the important part of the data in a dataframe

#prepare to choose two populations to compare
poplist_all<-unique(substr(names(filefas_all),1,2))
#the first 2 characters of the name of the individual are its location, that is why I can use substring(getName.SeqFastadna(filefas),1,2) to get the population of an individual. 
#use two populations at a time to get pairwise values
for (popnum1 in 1:(length(poplist_all)-1)){
#for (popnum1 in 1){ #use this line instead of the previous to try it out for one population 
	for (popnum2 in (popnum1+1):length(poplist_all)){
#	for (popnum2 in (popnum1+1)){ #use this line instead of the previous to try it out for one population 

pop1 = poplist_all[popnum1]
		pop2 = poplist_all[popnum2]
		print(paste("pop1",pop1))
		print(paste("pop2",pop2))
		poplist=poplist_all[c(popnum1,popnum2)]
# if you want to calculate D values for all populations together, uncomment the next line
#set poplist=poplist_all		

#make filefas with only two populations (using first two letters of names of inds)
filefas = filefas_all[		
		which(substr(names(filefas_all),1,2)==poplist[1]|substr(names(filefas_all),1,2)==poplist[2])]
		
#create dataframe with important information per sequence
DataInds<-data.frame("IND"=getName.SeqFastadna(filefas),"POP"=substring(getName.SeqFastadna(filefas),1,2),"POPRAN"=numeric(length(filefas[])),"SEQ"=numeric(length(filefas[])),"FRAG"=numeric(length(filefas[])),"Mt"=numeric(length(filefas[])))   
for (ind in 1:length(filefas[])){
	DataInds$SEQ[ind]=getSequence.SeqFastadna(filefas)[ind][[1]]}

#REMOVE gaps and missing data (whenever there is an "n" or a "-", make all seqs "n" at that nucleotide)
aa<-array(0,c(length(filefas[])))
for (nuc in 1:nchar(DataInds$SEQ[1])){
	aa=substring(DataInds$SEQ,nuc,nuc)
	if(length(which(aa=="n"))>0){substring(DataInds$SEQ,nuc,nuc)="n"}
	if(length(which(aa=="-"))>0){substring(DataInds$SEQ,nuc,nuc)="n"}
	}

#make POPS dataframe
numpops=length(poplist)
Pops=data.frame("POP"=poplist,"NUMIND"=numeric(length(poplist)),"NUMALL"=numeric(length(poplist)),"ENA"=numeric(length(poplist)),"HET"=numeric(length(poplist)),"HOM"=numeric(length(poplist)))
#COUNT INDIVIDUALS PER POP (or rather, number of sequences per pop)
for (pop in 1: numpops){
	Pops$NUMIND[which(Pops$POP==poplist[pop])]=sum(DataInds$POP==poplist[pop])
}

#calculate summary stats for varying length of seqs (including the whole length of the sequence)
numreps=NR 
#I ran each length twenty times for the paper
lengths=c(100,300,500,700,900,1100,nchar(DataInds$SEQ[1])) #to calculate D for different sequence lengths. 
#if you want to test only one length (e.g., the total length you have), use this line instead
#lengths=nchar(DataInds$SEQ[1])
numlengths=length(lengths)

#make data_frame to put the summary statistics for each length and each repeat
ll=numreps*numlengths
SumStats=data.frame("SPEC"=numeric(ll),"LENGTH"=numeric(ll),"HETS"=numeric(ll),"HOMS"=numeric(ll),"ENAS"=numeric(ll),"HETT"=numeric(ll),"HOMT"=numeric(ll),"ENAT"=numeric(ll),"NTILDE"=numeric(ll),"HS_est"=numeric(ll),"HT_est"=numeric(ll),"D_est"=numeric(ll),"D_est_pval"=numeric(ll),"G_est"=numeric(ll))

#here a loop starts which calculates all summary startistics for all lengths of DNA and a given number of times. 
		run=0
for (rep in 1:numreps) {
	print(paste("rep",rep,"of",numreps))
	for (len in lengths) {
		(run=run+1)
		print(paste("length",len))
		SumStats$LENGTH[run]=len

#fill DataInds$FRAG with fragment of genotype and $Mt with haplotype. Choose a random starting point in the sequence. If length is long, this starting point has to be at the beginning of the total sequence, otherwise it doesnt fit. 
		begin=sample(1:(nchar(DataInds$SEQ[1])-len+1),1)
		end=begin+(len)-1
		DataInds$FRAG=substring(DataInds$SEQ,begin,end)

#find all unique haplotypes and give individuals that have the same haplotype the same number.	
		haplolist<-unique(DataInds$FRAG)
		for (unihaplo in 1:length(haplolist)){
			for (ind in 1:length(DataInds[,1])){
				if(DataInds$FRAG[ind]==haplolist[unihaplo])
				{DataInds$Mt[ind]=unihaplo} 
			}}

#COUNT ALLELES PER POP AND ENA, HET, HOM
#counts the number of alleles per population and 
#calculates the effective number of alleles (ENA), heterozygosity (HET), homozygosity (HOM), and an unbiased estimator of within pop heterozygosity (HET_est) and an unbiased estimator of the effective number of alleles (ENA_est). 
for (pop in as.vector(poplist)){
	Pops$NUMALL[which(Pops$POP==pop)]=length(unique(DataInds$Mt[which(DataInds$POP==pop)]))
	alleles=DataInds$Mt[which(DataInds$POP==pop)]
	unique(alleles)->uni_alleles
	hom=0
	for (all in uni_alleles){
		freq=length(which(alleles==all))/length(alleles)
		hom=hom+(freq)^2
	}
	Pops$HOM[which(Pops$POP==pop)]=hom
	Pops$ENA[which(Pops$POP==pop)]=1/hom
	Pops$HET[which(Pops$POP==pop)]=1-hom
	}
Pops$HS_est=Pops$HET*(Pops$NUMIND/(Pops$NUMIND-1))
Pops$ENA_est=(1-Pops$HS_est)^(-1)

#ENA, HET, HOM total and ave over pops 
#I average freqs over pops, then square them to get homozygosity
alleles=DataInds$Mt
unique(alleles)->uni_alleles
hom=0
for (all in uni_alleles){
meanfreq=0
for (po in  as.vector(poplist)){
freq=length(which((DataInds$POP==po)&(DataInds$Mt==all)))/length(which(DataInds$POP==po))
meanfreq=meanfreq+freq}
meanfreq=meanfreq/numpops
hom=hom+meanfreq^2
}
SumStats$HOMT[run]=hom
SumStats$ENAT[run]=1/hom
SumStats$HETT[run]=1-hom
SumStats$HOMS[run]=mean(Pops$HOM)
SumStats$ENAS[run]=mean(Pops$ENA)
SumStats$HETS[run]=mean(Pops$HET)

#With the calculated homozygosities and heterozygosities etc, we can now get D_est and G_est
samplesizes<-Pops$NUMIND
samplesizes[samplesizes>0]->samplesizes #remove the zeros
SumOneOverN=0
for (sam in samplesizes){
SumOneOverN=SumOneOverN+(1/sam)}
SumStats$NTILDE[run]=length(samplesizes)/(SumOneOverN)
SumStats$HS_est[run]=(SumStats$NTILDE[run]/(SumStats$NTILDE[run]-1))*SumStats$HETS[run]
SumStats$HT_est[run]=SumStats$HETT[run]+(SumStats$HS_est[run]/(SumStats$NTILDE[run]*numpops))
SumStats$D_est[run]=((SumStats$HT_est[run]-SumStats$HS_est[run])/(1-SumStats$HS_est[run]))*((numpops)/(numpops-1))		
#if all sequences are different (no two inds carry same MtDNA, then both HS_est and HT_est are 1, and this leads to D_est = 2. This is clearly not correct. D_est should be 1, and pval should also be 1 in this case. 		
		if (!is.na(SumStats$D_est[run])){if (SumStats$D_est[run]==2){SumStats$D_est[run]=1} }		
		
SumStats$G_est[run]=((SumStats$HT_est[run]-SumStats$HS_est[run])/(SumStats$HT_est[run]))

#get p-value for D_est
#I have to permute individuals over populations and calculate HS_est (all else stays the same)
numran=NUMRAN  		
#create an array to store "numran" D_est values
D_est_Random<-array(0,c(numran))
#loop "numran" times
for (r in 1:numran){
#shuffle the ind over the pops
	DataInds$POPRAN=sample(DataInds$POP,replace=FALSE) #permute POPlabels
	for (pop in as.vector(poplist)){  #look at each of the population to get its heterozygosity. 
		indlist = which(DataInds$POPRAN==pop) #make a list of the individuals in that pop
		alleles=DataInds$Mt[indlist] # and of the alleles in that population
		unique(alleles)->uni_alleles # list of unique alleles in that population
		hom=0 #homozygosity = freq1^2 + freq2^2 etc for all alleles in a population
		for (all in uni_alleles){ # loop over the unique alleles in the population
			freq=length(which(alleles==all))/length(alleles) #calculate their frequency in the population
			hom=hom+(freq)^2 # add  the square of the freq to hom
		}
		Pops$HET[which(Pops$POP==pop)]=1-hom #heterozygosity for the pop is 1-hom 
	}
	HS_est_Random=(SumStats$NTILDE[run]/(SumStats$NTILDE[run]-1))*mean(Pops$HET) 
# get unbiased heterozygosity averaged over all pops 
	D_est_Random[r]=((SumStats$HT_est[run]-HS_est_Random)/(1-HS_est_Random))*((numpops)/(numpops-1)) 
# get D_est (at the end you ll have 1000 random D values
	}
#with 1000 random D values calculate a p value	(how often is the random D larger than the real D?
SumStats$D_est_pval[run]=sum(SumStats$D_est[run]<=D_est_Random)/length(D_est_Random)
}}

#D and G values below zero donot make sense, I replace them with zero, this has no effect on the p values, but it can change slightly the mean D or G value for shorter MtDNA sequences.
for (l in 1:ll){
	if (!is.na(SumStats$D_est[l])){	
		if (SumStats$D_est[l]<0){SumStats$D_est[l]=0}}
	if (!is.na(SumStats$G_est[l])){	
		if (SumStats$G_est[l]<0){SumStats$G_est[l]=0}}
	}

#Calculate means for each sequence length 
MeanSumStats=data.frame("SPEC"=numeric(numlengths),"LENGTH"=lengths,"HETS"=numeric(numlengths),"HOMS"=numeric(numlengths),"ENAS"=numeric(numlengths),"HETT"=numeric(numlengths),"HOMT"=numeric(numlengths),"ENAT"=numeric(numlengths),"NTILDE"=numeric(numlengths),"HS_est"=numeric(numlengths),"HT_est"=numeric(numlengths),"D_est"=numeric(numlengths),"D_est_pval"=numeric(numlengths),"G_est"=numeric(numlengths))

#calculate all the means (for each length of MtDNA you used)
for (i in 1:numlengths){
	MeanSumStats[i,]=
	colMeans(SumStats[which(SumStats$LENGTH==MeanSumStats$LENGTH[i]),])} 

#add column for the power of the D startistic
MeanSumStats$PowerD<-0
for (i in 1:numlengths){
MeanSumStats$PowerD[i]=length(which((SumStats$D_est_pval<0.05)&(SumStats$LENGTH==lengths[i])))/length(which((SumStats$LENGTH==lengths[i])))
}

#make a figures (like fig 4 in Pennings et al JEB 2011)
png(paste("Figure4",pop1,"_",pop2,".png",sep=""))
plot(MeanSumStats$LENGTH,MeanSumStats$D_est,ylim=c(0,1),xlim=c(100,1400),cex=1.5,yaxt="n", xaxt="n",xlab="length of MtDNA fragment",ylab="D, Gst, Power",pch="D",main=expression(paste("Effect of sequence length ", italic("P. americanus"), " MtDNA")),col=2)
# pch = 16 and co=1 in the paper
mtext(side=3,paste(pop1,"vs",pop2))
axis(1,labels=seq(100, 1300, by=200) , at=seq(100, 1300, by=200),cex.axis=0.8)
axis(2,labels=c(0,0.2,0.4,0.6,0.8,1.0), at=c(0,0.2,0.4,0.6,0.8,1.0),cex.axis=0.8)
points(lengths[1:8]+30,MeanSumStats$PowerD[1:8],cex=1.5,pch="P")
# pch = "*" in the paper
points(MeanSumStats$LENGTH,MeanSumStats$G_est,cex=1.5,pch="G",col="darkgrey") 
# pch = 17 in the paper
x=0.2
dev.off()

#keep your data
write.csv(MeanSumStats,paste("EffectOfLengthFigure_Stats",pop1,"_",pop2,".csv",sep=""))
	}
}		
		

