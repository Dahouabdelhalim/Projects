#Data<-read.table("P.am_tandem.csv", sep="\\t", dec=".",header=T,as.is = T) 
#Data<-read.table("T.lon_tandem.csv", sep="\\t", dec=".",header=T,as.is = T)
Data<-read.table("T.cur_tandem.csv", sep="\\t", dec=".",header=T,as.is = T)

species<-names(Data)[1]
numloci<-(length(Data[1,])-1)/2
for (i in 1:numloci){
names(Data)[i*2]<-Data[3,i*2]
names(Data)[i*2+1]<-paste(Data[3,i*2],"_2",sep="")
}
Data<-Data[-(1:4),]
#Data[,1]<-substring(Data[,1],1)
rowsOH2=which(substr(Data[,1],1,3)=="OH2")
if(length(rowsOH2)>0)substr(Data[rowsOH2,1],1,2)<-"O2"

loci<-names(Data)[(1:numloci)*2]
rowspopnames<-which(substr(Data[,1],1,3)=="pop")
Data<-Data[-rowspopnames,]
Data$POP<-substr(Data[,1],1,2)
poplist<-unique(Data$POP)
names(Data)[1]<-"Colony"

for (i in 1:numloci){
Data[,i*2]<-as.numeric(Data[,i*2])}

Pops=data.frame("POP"=rep(poplist,length(loci)),"LOCUS"=rep(loci,each=length(poplist)),"NUMIND"=numeric(length(poplist)*length(loci)),"NUMALL"=numeric(length(poplist)*length(loci)),"ENA"=numeric(length(poplist)*length(loci)),"HET"=numeric(length(poplist)*length(loci)),"HOM"=numeric(length(poplist)*length(loci)),"HETRAN"=numeric(length(poplist)*length(loci)),"ENA_est"=numeric(length(poplist)*length(loci)),"SPEC"=numeric(length(poplist)*length(loci)))

Pops$SPEC=species

#COUNT INDIVIDUALS
for (locus in 1:length(loci)){
ind=which(!is.na(Data[,locus*2]))
for (pop in poplist){
Pops$NUMIND[which((Pops$POP==pop)&(Pops$LOCUS==loci[locus]))]=sum(Data$POP[ind]==pop)
}}

#COUNT ALLELES & FREQUENCIES
cols=2:(numloci*2+1)
alleles=0
for (a in cols){alleles=c(alleles,Data[,a])}
alleles<-alleles[which(alleles>0)]
uni_alleles<-sort(unique(alleles[!is.na(alleles)]))
FreqsPop<-array(0,c(length(loci),length(poplist),length(uni_alleles)))

for (locus in 1:length(loci)){
for (pop in 1:length(poplist)){
popind = which(Data$POP==poplist[pop])
popalleles=c(Data[popind,locus*2],Data[popind,locus*2+1])
popalleles[!is.na(popalleles)]->popalleles
Pops$NUMALL[which((Pops$POP==poplist[pop])&(Pops$LOCUS==loci[locus]))]=length(unique(popalleles))
for (allele in 1:length(uni_alleles)){
FreqsPop[locus,pop,allele]=length(which(popalleles==uni_alleles[allele]))/length(popalleles)
}}}

#ENA, HET, HOM
for (locus in 1:length(loci)){
for (pop in 1:length(poplist)){
hom=sum(FreqsPop[locus,pop,]^2)
Pops$HOM[which((Pops$POP==poplist[pop])&(Pops$LOCUS==loci[locus]))]=hom
Pops$ENA[which((Pops$POP==poplist[pop])&(Pops$LOCUS==loci[locus]))]=1/hom
Pops$HET[which((Pops$POP==poplist[pop])&(Pops$LOCUS==loci[locus]))]=1-hom
}}
Pops$HET_est=(Pops$HET*2*Pops$NUMIND)/(2*Pops$NUMIND-1)
Pops$ENA_est=(1-Pops$HET_est)^(-1)

DataLoci=data.frame("SPEC"=numeric(length(loci)),"LOCUS"=loci,"HETS"=numeric(length(loci)),"HOMS"=numeric(length(loci)),"ENAS"=numeric(length(loci)),"HETT"=numeric(length(loci)),"HOMT"=numeric(length(loci)),"ENAT"=numeric(length(loci)),"NTILDE"=numeric(length(loci)),"HS_est"=numeric(length(loci)),"HT_est"=numeric(length(loci)),"D"=numeric(length(loci)),"D_pval"=numeric(length(loci)),"Gst"=numeric(length(loci)),"RANGE"=numeric(length(loci)),"NUMALLELES"=numeric(length(loci)),"M"=numeric(length(loci)), "MeanDRan"=numeric(length(loci)), "MedianDRan"=numeric(length(loci)))

DataLoci$SPEC=species

for (locus in 1:length(loci)){
localleles=unique(c(Data[,locus*2],Data[,locus*2+1]))
localleles[!is.na(localleles)]->localleles
DataLoci$RANGE[locus]=range(localleles)[2]-range(localleles)[1]
DataLoci$NUMALLELES[locus]=length(localleles)
DataLoci$M[locus]=DataLoci$NUMALLELES[locus]/DataLoci$RANGE[locus]
}

#ENA, HET, HOM total and ave over pops #average freqs over pops, then square them to get hom
FreqsAve<-array(0,c(length(loci),length(uni_alleles)))
for (locus in 1:length(loci)){
for (allele in 1:length(uni_alleles)){
FreqsAve[locus,allele]=mean(FreqsPop[locus,,allele])
}}
for (locus in 1:length(loci)){
hom=sum(FreqsAve[locus,]^2)
DataLoci$HOMT[locus]=hom
DataLoci$ENAT[locus]=1/hom
DataLoci$HETT[locus]=1-hom
DataLoci$HOMS[locus]=mean(Pops$HOM[which(Pops$LOCUS==loci[locus])],na.rm=T)
DataLoci$ENAS[locus]=mean(Pops$ENA[which(Pops$LOCUS==loci[locus])],na.rm=T)
DataLoci$HETS[locus]=mean(Pops$HET[which(Pops$LOCUS==loci[locus])],na.rm=T)
}

#D and G_st
for (locus in 1:length(loci)){
samplesizes<-Pops$NUMIND[which(Pops$LOCUS==loci[locus])]
samplesizes[samplesizes>0]->samplesizes #remove the zeros
n=length(samplesizes)
SumOneOverN=0
for (sam in samplesizes){
SumOneOverN=SumOneOverN+(1/sam)}
DataLoci$NTILDE[locus]=length(samplesizes)/(SumOneOverN)
DataLoci$HS_est[locus]=(2*DataLoci$NTILDE[locus]/(2*DataLoci$NTILDE[locus]-1))*DataLoci$HETS[locus]
DataLoci$HT_est[locus]=DataLoci$HETT[locus]+(DataLoci$HS_est[locus]/(2*DataLoci$NTILDE[locus]*n))
DataLoci$D[locus]=((DataLoci$HT_est[locus]-DataLoci$HS_est[locus])/(1-DataLoci$HS_est[locus]))*((n)/(n-1))
if(DataLoci$D[locus]<0){DataLoci$D[locus]=0}
DataLoci$Gst[locus]=((DataLoci$HT_est[locus]-DataLoci$HS_est[locus])/(DataLoci$HT_est[locus]))
if(DataLoci$Gst[locus]<0){DataLoci$Gst[locus]=0}
}

#Permutation test. Permute - Recalculate HS_est and Gst and D 
repeats=1000
RandomD=data.frame("L5"=numeric(repeats),"GT1"=numeric(repeats),"L18"=numeric(repeats),"Myrt3"=numeric(repeats),"GT218"=numeric(repeats),"L4"=numeric(repeats),"GT223"=numeric(repeats),"MS86"=numeric(repeats))
RandomGst=data.frame("L5"=numeric(repeats),"GT1"=numeric(repeats),"L18"=numeric(repeats),"Myrt3"=numeric(repeats),"GT218"=numeric(repeats),"L4"=numeric(repeats),"GT223"=numeric(repeats),"MS86"=numeric(repeats))
RandomHETS<-array(0,c(length(loci)))
RandomHS_est<-array(0,c(length(loci)))

#Permute - Recalculate HS_est and Gest and Dest
for (rep in 1:repeats){
Data$PopRAN<-sample(Data$POP,replace=FALSE)
for (locus in 1:length(loci)){
for (pop in 1:length(poplist)){
popind = which(Data$PopRAN==poplist[pop])
col=cols[locus*2-1]
#Data[popind,col]<-sample(Data[popind,col],replace=FALSE)#to randomize alleles!!
popalleles=c(Data[popind,col],Data[popind,col+1])
popalleles[!is.na(popalleles)]->popalleles
for (allele in 1:length(uni_alleles)){
FreqsPop[locus,pop,allele]=length(which(popalleles==uni_alleles[allele]))/length(popalleles)}
hom=sum(FreqsPop[locus,pop,]^2)
Pops$HETRAN[which((Pops$POP==poplist[pop])&(Pops$LOCUS==loci[locus]))]=1-hom}

RandomHETS[locus]=mean(Pops$HETRAN[which(Pops$LOCUS==loci[locus])],na.rm=T)
RandomHS_est[locus]=(2*DataLoci$NTILDE[locus]/(2*DataLoci$NTILDE[locus]-1))*RandomHETS[locus]
RandomD[rep,locus]=((DataLoci$HT_est[locus]-RandomHS_est[locus])/(1-RandomHS_est[locus]))*((n)/(n-1))
if (RandomD[rep,locus]<0){RandomD[rep,locus]=0}
RandomGst[rep,locus]=((DataLoci$HT_est[locus]-RandomHS_est[locus])/(DataLoci$HT_est[locus]))
if (RandomGst[rep,locus]<0){RandomGst[rep,locus]=0}
}}

for (locus in 1:length(loci)){
DataLoci$MeanDRan[locus]=mean(RandomD[,locus])
DataLoci$MedianDRan[locus]=median(RandomD[,locus])
}

# Determine p-values
Dpvalue<-array(0,c(length(loci)))
Gstpvalue<-array(0,c(length(loci)))
for (locus in 1:length(loci)){
Dpvalue[locus]=1-(length(which(DataLoci$D[locus]>RandomD[,locus]))/length(RandomD[,locus]))
DataLoci$D_pval[locus]=Dpvalue[locus]
Gstpvalue[locus]=1-(length(which(DataLoci$Gst[locus]>RandomGst[,locus]))/length(RandomGst[,locus]))
}

#write.csv(DataLoci,paste("DataLoci_",species,".csv",sep=""))
#write.csv(Pops,paste("Pops_",species,".csv",sep=""))
