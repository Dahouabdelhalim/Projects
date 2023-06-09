

#Load Paired T Test Randomisations
source('Paired T Randomisation .R', chdir = TRUE)	

#Genotype Data
markers<-read.table('Weaver Final Genotype Tables/Assembled Genotypes v2.txt',header=T)
mark<-markers[,-c(2,3)]

#Mean Heterozygosity Data
meanhetdata<-read.csv('Mean Heterozygosity by Locus.csv',header=T)

#################
# Hetsim function - requires following arguments 
	# 'dframe' - dataframe where data are contained
	# 'noff' - number of offspring to simulate
	# 'mumcolumn' - column of maternal ID
	# 'egpdadcolumn' - column of EG sire ID  
	# 'igpdadcolumn' - column of dominant male ID 
	
	#The function will then use the 'marker data' dataframe to drag out the genotype of the 3 breeders and simulate offspring genotypes following Mendelian inheritance of alleles 
	#Note that to make this code custom to your dataset, you will have to adjust the 'coords of markers' section as this is currently set up for 10 typed loci
	#Contact xav.harrison@gmail.com for assistance with this. I'd be happy to help
	
	#Toy example is provided at the bottom of this code
#################
hetsim2<-function(dframe,noff,mumcolumn,egpdadcolumn,igpdadcolumn){


	
	#Coords of markers
	cmat<-matrix(0,nrow=10,ncol=2)
	cmat[,1]<-seq(1,19,2)
	cmat[,2]<-seq(2,20,2)
	dadlocs<-seq(1,19,2)
	mumlocs<-seq(2,20,2)

#Blank Matrix for mean Het Values
hetoutput<-matrix(0,nrow=nrow(dframe),ncol=2)
colnames(hetoutput)<-c("igpmeanhet","egpmeanhet")


for (k in 1:nrow(dframe)){
	
	#Get Parental genotypes
	mum<-subset(mark,mark[,1] == as.character(dframe[k,mumcolumn]))
	mum<-mum[,-1]
	egpdad<-subset(mark,mark[,1] == as.character(dframe[k,egpdadcolumn]))
	egpdad<-egpdad[,-1]
	igpdad<-subset(mark,mark[,1] == as.character(dframe[k,igpdadcolumn]))
	igpdad<-igpdad[,-1]
		
		egphetvals<-numeric(noff)
		igphetvals<-numeric(noff)

for (j in 1:noff){

########IGP
	#Randomly Draw Haplotypes
	igpdadalleles<-as.numeric(igpdad[1,apply(cmat,1,function(x)sample(x,1))])
	mumalleles<-as.numeric(mum[1,apply(cmat,1,function(x)sample(x,1))])
	
		#Find Out which are zeroes and set male equivalent to zero
		igpdadalleles[which(mumalleles==0)]<-0
		mumalleles[which(igpdadalleles==0)]<-0
		
		#Change 1 of the haplotypes with missing data to 999 - will flag up after substraction
		mumalleles[which(mumalleles==0)]<-999
		
		#Subtract to Work out Which are Heterozygous
		gendif<-igpdadalleles-mumalleles
		
		#Flag 1 for heterozygous, 0 if not - keep -999
		gendif2<-ifelse(gendif==-999,-999,ifelse(gendif==0,0,1))
				
				#Add mean hets to missing data
				gendif2[which(gendif2==-999)]<-meanhetdata[which(gendif2==-999),3]
				igphetvals[j]<-sum(gendif2)/10
			
########EGP
	#Randomly Draw Haplotypes
	egpdadalleles<-as.numeric(egpdad[1,apply(cmat,1,function(x)sample(x,1))])
	mumalleles<-as.numeric(mum[1,apply(cmat,1,function(x)sample(x,1))])
	
		#Find Out which are zeroes and set male equivalent to zero
		egpdadalleles[which(mumalleles==0)]<-0
		mumalleles[which(egpdadalleles==0)]<-0
	
		#Change 1 of the haplotypes with missing data to 999 - will flag up after substraction
		mumalleles[which(mumalleles==0)]<-999
		
		#Subtract to Work out Which are Heterozygous
		gendif<-egpdadalleles-mumalleles
		
		#Flag 1 for heterozygous, 0 if not - keep -999
		gendif2<-ifelse(gendif==-999,-999,ifelse(gendif==0,0,1))
				
				#Add mean hets to missing data
				gendif2[which(gendif2==-999)]<-meanhetdata[which(gendif2==-999),3]
				egphetvals[j]<-sum(gendif2)/10
				
					} #offspring loop end

hetoutput[k,1]<-mean(igphetvals)
hetoutput[k,2]<-mean(egphetvals)

}
hetoutput
}

									######################################################
									#
									# 			Toy Example and Analysis
									#
									######################################################

#Simulate Nonsense genotypes (most individuals will be heterozygous)
mark2<-as.data.frame(matrix(rpois(200,100),ncol=20))
mark<-cbind(paste("bird",seq(10),sep=""),mark2)

#Assemble Dataframe of Candidate Parents
toyframe<-data.frame(mum=c("bird1","bird2","bird3"),dom=c("bird5","bird4","bird6"),egp=c("bird7","bird8","bird9"))

#Load a 'meanhetdata' file (wont be used for this example but required for function to run)
meanhetdata<-data.frame(X=seq(10),markers=paste("marker",seq(10),sep=""),meanhets=runif(10))

#Simulate 20 offspring from dataframe
x<-hetsim2(toyframe,20,"mum","dom","egp")
x #almost everyone is heterozygous

#Test for significant differences (will throw a p value of 0 for this small dataset)
x2<-pair.randomise(x,1,2,10000)
