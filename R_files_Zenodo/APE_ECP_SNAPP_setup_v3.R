### Set up SNAPP files for testing sensitivity to geographic sampling in APE_ECP drainage ###
### Libraries ###
require("adegenet")

### remove scientific notation ###
options(scipen=999)

### Read in bigger GenePop flie
og<-read.genepop("~/Dropbox/InsectTaxa/genepop/Ape.gen")
v1<-which(sapply(strsplit(rownames(og@tab),"_"),function(x) x[3]) =="ECA4248")
v2<-grep("ECP",sapply(strsplit(rownames(og@tab),"_"),function(x) x[3]))
gp<-og[c(v1,v2),]

### Read in CLUMPP K 2 file to get Q scores for taxon assignment ###
ecp<-read.genepop("~/Dropbox/InsectTaxa/genepop/ApeECP.gen")
ecp_k2<-readLines(file("~/Dropbox/InsectTaxa/dapc_and_structure/Ape_ECP_Struc/StructureHarvesterOutput 2/K2_CLUMPP.OUT"))

q<-sapply(strsplit(ecp_k2,":"),function(x) x[length(x)])
q_df<-data.frame(lapply(strsplit(q," "),function(x) x[3:4]),stringsAsFactors=F)
colnames(q_df)<-NULL
ecp_k2_q_df<-as.matrix(q_df)
Q_vec<-as.numeric(ecp_k2_q_df[1,])
names(Q_vec)<-row.names(ecp@tab)

### SET UP THIS DATA FRAME TO INCLUDE THE OUTGRUP !!!!
ass<-data.frame(row.names=rownames(ecp@tab),taxon=paste("taxon",as.numeric(Q_vec>0.5)+1,sep=""))
ass<-rbind(ass,data.frame(row.names=rownames(og@tab)[v1],taxon=rep("taxon3",length(v1))) )
ass<-cbind(ass,site=sapply(strsplit(row.names(ass),"_"),function(x) x[3]))

## Remove outermost 4 populations to facilitate running SNAPP ###
ass<-ass[!ass[,2] %in% levels(ass[,2])[c(2,3,13,14)],]
ass[,2]<-factor(ass[,2])

### Set up basic parameters
n_ingroup<-36
n_per_pop<-4
n_loci<-100

### Create list for subsampling scenarios ###
scen_list<-list()
scen_list[[1]]<-levels(ass[,2])[-1][rep(T,9)]
scen_list[[2]]<-levels(ass[,2])[-1][c(rep(T,4),F,rep(T,4))]
scen_list[[3]]<-levels(ass[,2])[-1][c(rep(T,3),rep(F,3),rep(T,3))]
scen_list[[4]]<-levels(ass[,2])[-1][c(rep(T,2),rep(F,5),rep(T,2))]

#GenePop2SNAPP<-function(inFile,outFile,assignmentFile,n_per_pop,n_loci,writeFile=T){
#gp<-read.genepop(inFile)

genotypes<-gp@tab
genotypes[is.na(genotypes)]<-"?"
	
npop<-length(unique(ass[,2]))

### Remove individuals with > 50% missing data, reexamine sample sizes per site 
missing_data<-vector()
for(i in 1:nrow(genotypes)){
	missing_data[i]<-table(genotypes[i,])["?"]/sum(table(genotypes[i,]))
}

missing_data

gp_red<-gp[which(missing_data<0.5),]

table(ass[which(missing_data<0.5),2]) ### Still enough to get random samples 

### Create n_rep Replicates of Random samples with n_per_pop for each population and all of outgroup ###
og_samp<-sample(rownames(ass)[ass$site=="ECA4248"],8) #Keep same OGs for all sampling scenarios

### SET UP PARENTAL POPS TO BE PURE PARENTALS ###
parental_list<-list()
for(i in 1:length(scen_list[[4]])){
	foo<-Q_vec[sapply(strsplit(names(Q_vec),"_"),function(x) x[3]) %in% scen_list[[4]][i]]
	foo<-foo[names(foo) %in% row.names(gp_red@tab)]
	parental_list[[i]]<-names(foo[foo < 0.025 | foo > 0.975])
}
names(parental_list)<-scen_list[[4]]

### Set up admixed indivuals to have at least some level of admixture ###
admixed_list<-list()
admixed_pops<-scen_list[[1]][!scen_list[[1]]%in%scen_list[[4]]]
for(i in 1:length(admixed_pops)){
	foo<-Q_vec[sapply(strsplit(names(Q_vec),"_"),function(x) x[3]) %in% admixed_pops[i]]
	admixed_list[[i]]<-names(sort(abs(foo-0.5)))[names(sort(abs(foo-0.5))) %in% row.names(gp_red@tab)]
}
names(admixed_list)<-admixed_pops

scen_list

samp_list <-list()
samp_list[[1]]<-list()
samp_list[[1]][[1]]<-og_samp
for(j in 2:length(levels(ass[,2]))){
	if(levels(ass[,2])[j] %in% names(parental_list)){
		samp_list[[1]][[j]]<-sample(parental_list[[levels(ass[,2])[j]]], 9) 
		}else{
			if(levels(ass[,2])[j] %in% names(admixed_list)){
				samp_list[[1]][[j]]<-admixed_list[[levels(ass[,2])[j]]][1:9]
			}
		}	
	}
				
names(samp_list[[1]])<-levels(ass[,2])

for(i in 2:4){
	samp_list[[i]]<-list()
	samp_list[[i]]<-samp_list[[1]]
}

for(i in 1:length(scen_list)){
	nperpop<-floor(n_ingroup/length(which(names(samp_list[[i]])[-1] %in% c(scen_list[[i]]))))
	for(j in 2:length(levels(ass[,2]))){
		if(names(samp_list[[i]])[j] %in% scen_list[[i]]){
			samp_list[[i]][[names(samp_list[[i]])[j]]]<-samp_list[[i]][[names(samp_list[[i]])[j]]][1:nperpop]
			samp_list[[i]][[names(samp_list[[i]])[j]]]<-sample(samp_list[[i]][[names(samp_list[[i]])[j]]],length(samp_list[[i]][[names(samp_list[[i]])[j]]]))
		}else{
			samp_list[[i]][[names(samp_list[[i]])[j]]]<-NA
		}	
	}
	samp_list[[i]]<-samp_list[[i]][!sapply(samp_list[[i]],function(x) is.na(x[1]))]
}	
sapply(samp_list,function(x) length(unlist(x)))

gp_list<-list()
for(i in 1:length(samp_list)){
	gp_list[[i]]<-gp_red[unlist(samp_list[[i]]),]
}

### Subsample loci to maximize data availability and minor allele frequencies ###
gp_red2<-gp_red[unlist(samp_list)[!duplicated(unlist(samp_list))],]

missing_data<-list()
maf<-vector()
for(j in 1:length(levels(gp_red2@loc.fac))){
	missing_data[[j]]<-vector()
	for(k in 1:length(levels(gp_red2@pop))){
		foo<-gp_red2@tab[which(gp_red2@pop ==levels(gp_red2@pop)[k]),which(gp_red2@loc.fac == levels(gp_red@loc.fac)[j])]
		missing_data[[j]][k]<-length(which(is.na(foo[,1])))/nrow(foo)
	}
	foo<-gp_red2@tab[,which(gp_red2@loc.fac == levels(gp_red@loc.fac)[j])]
	bar<-apply(foo,2,function(x) sum(x,na.rm=T))
	maf[j]<-min(bar)/sum(bar)
}


loci_filt<-sample(levels(gp_list[[i]]@loc.fac)[which(!sapply(missing_data,function(x) any(x==1)) & maf > 0.05)],100)

### Filter by missing data and MAF for whole group ###
gp_locisel_list<-list()
for(i in 1:length(gp_list)){
	#Constrain to biallelic loci#
	gp_locisel_list[[i]]<-gp_list[[i]][,gp_list[[i]]@loc.fac %in% names(gp_list[[i]]@loc.n.all)[which(gp_list[[i]]@loc.n.all == 2)]]
	
	#Filter to those loci included above#
	gp_locisel_list[[i]]<-gp_list[[i]][,gp_list[[i]]@loc.fac %in% loci_filt]
	}
			
rownames(gp_locisel_list[[1]]@tab) == 	unlist(samp_list[[1]]) #ALL TRUE ALL GOOD
	
### We've subsampled loci for each replicate ### Now to create xml output files for each replicate, each geosampling scenario, and each speciesdelim scenario
	
# ### Create subsample data sets ###
# ### Create subsampling figure ###
# pdf(file="~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/Figures/sampscheme_v2.pdf",width=3.25,height=2)
# #quartz()
# par(mfrow=c(4,1))
# par(mar=c(0,0,0,0))
# plot(1:length(levels(ass[,2])[-1]), rep(9,length(levels(ass[,2])[-1])),cex=3,axes=F,xlab="",ylab="",pch=21,bg="black",col="black")
# mtext("Sampling Scenario 1",side=1,cex=0.7,line=-1)
# plot(1:length(levels(ass[,2])[-1]), rep(9,length(levels(ass[,2])[-1])),cex=3,axes=F,xlab="",ylab="",pch=21,bg=c(rep("black",4),"white",rep("black",4)), col="black")
# mtext("Sampling Scenario 2",side=1,cex=0.7,line=-1)
# plot(1:length(levels(ass[,2])[-1]), rep(9,length(levels(ass[,2])[-1])),cex=3,axes=F,xlab="",ylab="",pch=21,bg=c(rep("black",3),rep("white",3),rep("black",3)), col="black")
# mtext("Sampling Scenario 3",side=1,cex=0.7,line=-1)
# plot(1:length(levels(ass[,2])[-1]), rep(9,length(levels(ass[,2])[-1])),cex=3,axes=F,xlab="",ylab="",pch=21,bg=c(rep("black",2),rep("white",5),rep("black",2)), col="black")
# mtext("Sampling Scenario 4",side=1,cex=0.7,line=-1)
# dev.off()
### Create nexus file to generate appropriate priors ###

for(m in 1:length(gp_locisel_list)){
	nexus_loci_list<-list()
	for(i in 1:length(levels(gp_locisel_list[[m]]@loc.fac))){
		foo<-gp_locisel_list[[m]]@tab[,gp_locisel_list[[m]]@loc.fac==levels(gp_locisel_list[[m]]@loc.fac)[i]][,1]
		foo[is.na(foo)]<-"?"
		nexus_loci_list[[i]]<-foo
	}
	
	nexus_df<-data.frame(nexus_loci_list,stringsAsFactors=F)
	colnames(nexus_df)<-names(gp_locisel_list[[m]]@all.names)
	
	sink(paste0("~/Desktop/SNAPP_SamplingScenario_",m,".nex"))
	cat("#NEXUS\\n")
	cat("\\n")
	cat("BEGIN DATA;")
	cat("\\n")
	cat(paste0("\\tDIMENSIONS NTAX=",nrow(nexus_df)," NCHAR=",ncol(nexus_df),";\\n"))
	cat("\\tFORMAT DATATYPE=INTEGER SYMBOLS=\\"012\\" MISSING=? GAP=- INTERLEAVE=NO;\\n")
	cat("\\tMATRIX\\n")
	
	for(i in 1:nrow(nexus_df)){
			cat(paste0("\\t\\t",rownames(nexus_df)[i],"\\t\\t\\t\\t",paste(nexus_df[i,],collapse=""),"\\n"))	
	}
	cat("\\t\\t;\\n")
	cat("END;")
	sink()
}

gp_locisel_list[[1]]

### For each replicate, each sampling scenario, and each species delimitation scenario, generate an xml file for path sampling analysis ###
### Read template xml ###
template<-readLines(file("~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/SNAPP/Empirical/Input/Rep_1_sampscen_4_spDscen_2_100loci_v13.xml"))

#### Set up variables ###
id_name<-"test"
n_steps_bfd<-24
alpha<-0.3
chainLength<-500000
burnInPercentage<-20
preBurnin<-50000

template<-gsub("test",id_name,template)


#for(k in 1:length(gp_locisel_list)){
	k<-1
	
for(m in 1:length(gp_locisel_list)){

### Create separate species assignment runs
for(n in 1:2){

sink(file=paste0("~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/SNAPP/Rep_",k,"_sampscen_",m,"_spDscen_",n,"_100loci_v15.xml"))

cat(template[1:4],sep="\\n")

### print sequence data for all individuals included ###
for(i in 1:nrow(gp_locisel_list[[m]]@tab)){
	seqfoo<-vector()
	for(j in 1:n_loci){
		seqfoo <-c(seqfoo,gp_locisel_list[[m]]@tab[i,gp_locisel_list[[m]]@loc.fac==names(gp_locisel_list[[m]]@all.names)[j]][1])
	}
	
	seqfoo[is.na(seqfoo)]<-"?"	
	seqfoo<-paste(seqfoo,collapse="")
	grep("[?]",strsplit(seqfoo,"")[[1]])
	cat("\\t\\t\\t")
	cat("<sequence id=\\"")
	cat(paste("seq", rownames(gp_locisel_list[[m]]@tab)[i], ass[rownames(gp_locisel_list[[m]]@tab)[i],1],sep="_"))
	cat("\\" taxon=\\"")
	cat(paste(rownames(gp_locisel_list[[m]]@tab)[i], ass[rownames(gp_locisel_list[[m]]@tab)[i],1],sep="_"))
	cat("\\" totalcount=\\"3\\" value=\\"")
	cat(seqfoo)
	cat("\\"/>")
	cat("\\n")
	
}

cat("\\n")
cat("\\t</data>")
cat("\\n\\n")

cat(template[grep("<map name",template)[1]:grep("<map name",template)[length(grep("<map name",template))]],sep="\\n")

### Set up BFD Parameters ###
cat("\\n")

cat("<run spec=\\"beast.inference.PathSampler\\"\\n")
cat(paste0("\\tchainLength=\\"",as.character(chainLength),"\\"\\n"))
cat(paste0("\\talpha=\\"",alpha,"\\"\\n"))
cat(paste0("\\trootdir=\\""))
cat(paste0("/workdir/nam232/Rep_",k,"_sampscen_",m,"_spDscen_",n,"\\""))
cat("\\n")
cat(paste0("\\tburnInPercentage=\\"", burnInPercentage,"\\"\\n"))
cat(paste0("\\tpreBurnin=\\"", preBurnin,"\\"\\n"))
cat("\\tdeleteOldLogs=\\"true\\"\\n")
cat(paste0("\\tnrOfSteps=\\"", n_steps_bfd,"\\">\\n"))
cat("\\n")

cat(template[grep("cd ",template):grep("<rawdata",template)],sep="\\n")

if(n==1){
	cat("\\t\\t\\t\\t<taxonset id=\\"taxon_three\\" spec=\\"TaxonSet\\">\\n")
	og<-grep("ECA",rownames(gp_locisel_list[[m]]@tab),value=T)
		for(i in 1:length(og)){
			cat("\\t\\t\\t\\t\\t<taxon id=\\"")
			cat(paste(og[i],ass[og[i],1],sep="_"))
			cat("\\" spec=\\"Taxon\\"/>")
			cat("\\n")
		}
	cat("\\t\\t\\t\\t</taxonset>\\n")
	cat("\\t\\t\\t\\t<taxonset id=\\"taxon_onetwocomb\\" spec=\\"TaxonSet\\">\\n")
	
	ig<-rownames(gp_locisel_list[[m]]@tab)[!rownames(gp_locisel_list[[m]]@tab) %in% og]
	
	for(i in 1:length(ig)){
			cat("\\t\\t\\t\\t\\t<taxon id=\\"")
			cat(paste(ig[i],ass[ig[i],1],sep="_"))
			cat("\\" spec=\\"Taxon\\"/>")
			cat("\\n")
		}
			
	cat("\\t\\t\\t\\t</taxonset>\\n")
	}else{
		cat("\\t\\t\\t\\t<taxonset id=\\"taxon_three\\" spec=\\"TaxonSet\\">\\n")
		og<-grep("ECA",rownames(gp_locisel_list[[m]]@tab),value=T)
			for(i in 1:length(og)){
				cat("\\t\\t\\t\\t\\t<taxon id=\\"")
				cat(paste(og[i],ass[og[i],1],sep="_"))
				cat("\\" spec=\\"Taxon\\"/>")
				cat("\\n")
			}
		cat("\\t\\t\\t\\t</taxonset>\\n")
		cat("\\t\\t\\t\\t<taxonset id=\\"taxon_one\\" spec=\\"TaxonSet\\">\\n")
	
		t1<-rownames(ass)[ass[,1]=="taxon1"] [rownames(ass)[ass[,1]=="taxon1"] %in% rownames(gp_locisel_list[[m]]@tab)]
			for(i in 1:length(t1)){
				cat("\\t\\t\\t\\t\\t<taxon id=\\"")
				cat(paste(t1[i],ass[t1[i],1],sep="_"))
				cat("\\" spec=\\"Taxon\\"/>")
				cat("\\n")
			}
		
		cat("\\t\\t\\t\\t</taxonset>\\n")
		cat("\\t\\t\\t\\t<taxonset id=\\"taxon_two\\" spec=\\"TaxonSet\\">\\n")
		
		t2<-rownames(ass)[ass[,1]=="taxon2"] [rownames(ass)[ass[,1]=="taxon2"] %in% rownames(gp_locisel_list[[m]]@tab)]
			for(i in 1:length(t2)){
				cat("\\t\\t\\t\\t\\t<taxon id=\\"")
				cat(paste(t2[i],ass[t2[i],1],sep="_"))
				cat("\\" spec=\\"Taxon\\"/>")
				cat("\\n")
			}
			cat("\\t\\t\\t\\t</taxonset>\\n")

						
	}
	
	cat(template[grep("</taxa>",template):length(template)],sep="\\n")
	sink()
}	

}

### Create batch running script ###
### Medium memory has 40 cores, use 32 / 8 = 4 per run ###

sink("~/Desktop/bfd_v6.sh")
for(m in 1:length(gp_locisel_list)){
	for(n in 1:2){
		cat(paste0("nohup java -jar /programs/beast/lib/beast.jar -threads 4 /workdir/nam232/Rep_",k,"_sampscen_",m,"_spDscen_",n,"_100loci_v15.xml > nohup_",m,"_",n,".out &\\n"))
	}
}
sink()

### Examine Q Scores for individuals in this run ###
foo<-na.omit(Q_vec[row.names(gp_locisel_list[[1]]@tab)])
bg_vec<-c(rep("red",8),rep("purple",20),rep("blue",8))

seq(1,length(foo),4)

### Heuristic figure ###
pdf(file="~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/Figures/sampscheme_v4.pdf",width=6.5,height=4)
quartz(width=6.5,height=4)
layout(matrix(c(1,1,1,1,2,3,4,5),ncol=2))
par(mar=c(4,4,0.5,0.5))
plot(1:length(foo),foo,pch=21,bg=bg_vec,cex=2,xlab="Sampling Site",ylab="Q score")

par(mar=c(0,2,0,2))
par(xpd=NA)
plot(1:length(levels(ass[,2])[-1]), rep(9,length(levels(ass[,2])[-1])),cex=3,axes=F,xlab="",ylab="",pch=21,bg=c(rep("red",2),rep("purple",5),rep("blue",2)), col="black")
mtext("Sampling Scenario 1",side=1,cex=0.7,line=-1)
plot(1:length(levels(ass[,2])[-1]), rep(9,length(levels(ass[,2])[-1])),cex=3,axes=F,xlab="",ylab="",pch=21,bg=c(rep("red",2),"purple","purple","white","purple","purple",rep("blue",2)), col="black")
mtext("Sampling Scenario 2",side=1,cex=0.7,line=-1)

plot(1:length(levels(ass[,2])[-1]), rep(9,length(levels(ass[,2])[-1])),cex=3,axes=F,xlab="",ylab="",pch=21,bg=c(rep("red",2),"purple","white","white","white","purple",rep("blue",2)), col="black")
mtext("Sampling Scenario 3",side=1,cex=0.7,line=-1)
plot(1:length(levels(ass[,2])[-1]), rep(9,length(levels(ass[,2])[-1])),cex=3,axes=F,xlab="",ylab="",pch=21,bg=c(rep("red",2),"white","white","white","white","white",rep("blue",2)), col="black")
mtext("Sampling Scenario 4",side=1,cex=0.7,line=-1)


dev.off()


admixed_list[1]
rep_list[[1]]
