require(ape)
require(xlsx)
require(phytools)
require(BAMMtools)

### Reading in song data ###
t.song<-read.xlsx("~/Desktop/Manuscripts/Oscinevssuboscinerates/SongData/Thraupidae_AppendixA.xlsx",sheetName="Thraupidae_AppendixA.csv",row.names=1)

### Reading in VP indices ###
t.vp<-read.xlsx("~/Desktop/Manuscripts/Oscinevssuboscinerates/SongData/Thraupidae_AppendixA.xlsx",sheetName="Vocal Performance - upper bound",row.names=1,stringsAsFactors=F)

t.song<-cbind(t.song,t.vp[8]) # ADD VP to song df
rownames(t.song)<-gsub(" ","_",rownames(t.song))

f.song<-read.xlsx("~/Desktop/Manuscripts/Oscinevssuboscinerates/SongData/Furn_song.xlsx",sheetName="Main dataset.xls",row.names=1)

### Matching character names
colnames(f.song)[c(4,5,6,7,8,9,10,11)]<-colnames(t.song)[c(1,11,19,13,2,18,10,9)]

song<-rbind(t.song[c(1,11,19,13,2,18,10,9,22)],f.song[c(4,5,6,7,8,9,10,11,12)])

### Fixing rownames of song to match tree
#song<-t.song
rownames(song)<-sapply(strsplit(row.names(song),"_"),function(x) paste(x[1:2],collapse="_"))
rownames(song)<-gsub("bouvreuil","pileata",row.names(song))
rownames(song)<-gsub("Buthraupis","Cnemathraupis",row.names(song))
rownames(song)<-gsub("Cnemathraupis_montana","Buthraupis_montana",row.names(song))
rownames(song)<-gsub("Cnemathraupis_wetmorei","Buthraupis_wetmorei",row.names(song))
rownames(song)<-gsub("bairdii","bairdi",row.names(song))
rownames(song)<-gsub("Thraupis_bonariensis","Pipraeidea_bonariensis",row.names(song))
rownames(song)<-gsub("Piezorhina_cinerea","Piezorina_cinerea",row.names(song))
rownames(song)<-gsub("Delothraupis_castaneoventris","Dubusia_castaneoventris",row.names(song))
## Tanagers are taken care of at this point
### Writing out csv file for Jonathan Drury of tanager song with updated taxonomy ###
#write.csv(song,file="~/Desktop/Thraupidae_song_03Jun2016.csv",quote=F)

rownames(song)<-gsub("'","",row.names(song))
rownames(song)<-gsub("Syndactyla_dimidiatum","Syndactyla_dimidiata",row.names(song))
rownames(song)<-gsub("Certhiasomus_stictolaemus","Deconychura_stictolaema",row.names(song)

# tree.trim<-list()
# for(i in 1:length(tree)){
	# tree.trim[[i]]<-tree[[i]]
	# tree.trim[[i]]<-drop.tip(tree.trim[[i]],tree.trim[[i]]$tip.label[!tree.trim[[i]]$tip.label%in% rownames(song)])
# }

# ### Song data and tree match! ###
# sapply(tree.trim,function(x) rownames(song) %in% x$tip.label)
# sapply(tree.trim, function(x) x$tip.label %in% rownames(song))
# class(tree.trim)<-"multiPhylo"

### Write out trimmed phylogeny ###
# write.nexus(tree.trim,file="~/Desktop/Manuscripts/Oscinevssuboscinerates/Manuscript/SupplementaryFiles/ThraupidFurnariid500TrimmedTrees.nex")

### March 26, 2015 Adding Mass to the Mix
ob.morph<-read.xlsx("~/Desktop/Manuscripts/Oscinevssuboscinerates/SongData/ovenbird_morphdata.xlsx",sheetName="Main dataset.xls",stringsAsFactors=F) #Good to Go
rownames(ob.morph)<-ob.morph$Species
colnames(ob.morph)[6]<-"Mean"

## Reading in masses ##
mass <- read.csv("~/Desktop/Manuscripts/TanagerBillsandSong/Data/AppendixC.csv", row.names = 1, stringsAsFactors = F)
rownames(mass) <- gsub(" ", "_", rownames(mass))
mass <- mass[!is.na(mass$Mean), ] #Remove NAs
mass.n<-mass[2]
mass<-mass[3]

## Combine thraupid furnariid masses
mass<-rbind(mass,ob.morph[6])

### Match Mass and Song dataset
song.trim<-song[rownames(song) %in% rownames(mass),]
mass.trim<-mass[rownames(mass) %in% rownames(song),]

### combine mass and song dataset ###
song_mass<-merge(mass,song,by=0)
rownames(song_mass)<-song_mass$Row.names
song_mass$Row.names<-NULL

### Create supplementary appendix with song and mass data ###
charlab<-gsub("[.]"," ",colnames(song_mass))
charlab<-gsub("VP","Vocal\\nperformance",charlab)
charlab<-gsub("Song frequency range","Song\\nfrequency range",charlab)
charlab<-gsub("Low ","Minimum\\n",charlab)
charlab<-gsub("High ","Maximum\\n",charlab)
charlab[1]<-"Mass"

song_mass2<-song_mass
colnames(song_mass2)<-charlab
charlab[grep("frequency",charlab)]<-paste("log(",charlab[grep("frequency",charlab)],")",sep="")

write.xlsx(song_mass2,file="/Users/NickMason/Desktop/Manuscripts/Oscinevssuboscinerates/AdditionalSupplementaryFiles/ThraupidFurnariid_SongMass_v2.xlsx")

song_mass2<-read.xlsx("/Users/NickMason/Desktop/Manuscripts/Oscinevssuboscinerates/AdditionalSupplementaryFiles/ThraupidFurnariid_SongMass_v2.xlsx",sheetName="Sheet1",stringsAsFactors=F,row.names=1)

### Read in tree ###
mcc<-read.nexus("~/Desktop/Manuscripts/Oscinevssuboscinerates/Trees/ThraupidFurnariidSupertreeTrimmed.nex")
mcc.trim<-drop.tip(mcc,mcc$tip.label[!mcc$tip.label %in% rownames(song_mass2)])
colnames(song_mass2)

## Ladderize phylogeny, revisions for 22 Nov 2016
mcc.trim<-ladderize(mcc.trim)

### Check out phylogeny and node numbers corresponding to tanagers / furnariids
plot(mcc.trim,show.tip.label=F)
nodelabels(cex=0.5)
### Number of songs, taxa, sd, etc
Ntip(mcc.trim)

### Tanagers
Ntip(extract.clade(mcc.trim,858)) #Number of tanagers included
sum(song_mass2[rownames(song_mass2) %in% extract.clade(mcc.trim,858)$tip.label,]$n)
sd(song_mass2[rownames(song_mass2) %in% extract.clade(mcc.trim,858)$tip.label,]$n)
mean(song_mass2[rownames(song_mass2) %in% extract.clade(mcc.trim,858)$tip.label,]$n)/sqrt(Ntip(extract.clade(mcc.trim,858)))

### Ovenbirds
Ntip(extract.clade(mcc.trim,583)) #Number of tanagers included
sum(song_mass2[rownames(song_mass2) %in% extract.clade(mcc.trim, 583)$tip.label,]$n)
sd(song_mass2[rownames(song_mass2) %in% extract.clade(mcc.trim, 583)$tip.label,]$n)
mean(song_mass2[rownames(song_mass2) %in% extract.clade(mcc.trim, 583)$tip.label,]$n)/sqrt(Ntip(extract.clade(mcc.trim, 583)))

sum(song_mass2[rownames(song_mass2) %in% extract.clade(mcc.trim,858)$tip.label,]$n)+sum(song_mass2[rownames(song_mass2) %in% extract.clade(mcc.trim, 583)$tip.label,]$n)

### NEW AS OF JAN 2016 ###
### Log transform characters with frequency units ###
freq_char<-grep("frequency",colnames(song_mass))

for(i in 1:length(freq_char)){
	song_mass[,freq_char[i]]<-log(song_mass[,freq_char[i]])
}

#Remove n from data frame
song_mass$n<-NULL

### genereating phylogenetic PCAs ###
### Examine correlations with mass and take residuals as needed ###
require(nlme)
output.list<-list()
head(song_mass)
for(i in 1:(ncol(song_mass)-1)){
	output.list[[i]]<-gls(as.formula(paste(colnames(song_mass)[i+1]," ~ Mean",sep="")),data=song_mass,correlation=corPagel(1,mcc.trim))
}

output.summary<-lapply(output.list,summary)
char_correlations<-sapply(output.summary,function(x) x$tTable[2,4])<0.05
corr_vec<-which(char_correlations)

### Create supplementary table for mass correlations ###
beta_se<-paste(round(sapply(output.summary,function(x) x$tTable[2,1]),3),round(sapply(output.summary,function(x) x$tTable[2,2]),3),sep=" Â± ")
pagel_lambda<-round(sapply(output.summary,function(x) x$modelStruct$corStruct[1]),2)
p_values<-round(sapply(output.summary,function(x) x$tTable[2,4]),3)

pgls.out<-data.frame(character=colnames(song_mass)[-1],lambda=pagel_lambda,beta=beta_se,p=p_values)
rownames(pgls.out)<-pgls.out$character
pgls.out$character<-NULL

pgls.out$p[as.character(pgls.out$p)=="0"]<-"<0.001"
pgls.out
write.xlsx(pgls.out,"/Users/NickMason/Desktop/Manuscripts/Oscinevssuboscinerates/Tables/Mass_PGLS_out_v2.xlsx")

#### Set up mass vector ####
mass.vec<-song_mass[,1]
names(mass.vec)<-row.names(song_mass)
as.matrix(song_mass[corr_vec+1])

### This was reversed but is now fixed as of 08 Nov 2015 and correspondence with Liam Revell
song_corr_resid<-phyl.resid(mcc.trim, mass.vec,as.matrix(song_mass[corr_vec+1]),method="lambda")
save(song_corr_resid,file="song_corr_resid_v2.Rdata")
load("song_corr_resid_v2.Rdata")

### Combine residual data set with raw data from characters that are uncorrelated with mass
song_mass_resid<-song_mass
for(i in 1:length(corr_vec)){
	song_mass_resid[,corr_vec[i]+1]<-song_corr_resid$resid[,i]
}

#song.massresid.ppca<-phyl.pca(mcc.trim,scale(song_mass_resid[-1]),method="lambda") #This takes a long time
#song.ppca<-phyl.pca(mcc.trim,scale(song_mass[-1]),method="lambda") #This takes a long time

save(song.massresid.ppca,file="song.massresid.ppca.v2.Rdata")
save(song.ppca,file="song.ppca.v4.Rdata")

load("song.massresid.ppca.v2.Rdata")
load("song.ppca.v4.Rdata")

### Loadings of PCA axes ###
cumsum(diag(song.massresid.ppca$Eval))/sum(diag(song.massresid.ppca$Eval))#PC1-5 are close to 95% of variation
cumsum(diag(song.ppca$Eval))/sum(diag(song.ppca$Eval))#PC1-5 are close to 95% of variation

round(diag(song.massresid.ppca$Eval)/sum(diag(song.massresid.ppca$Eval))*100,2)
round(diag(song.ppca$Eval)/sum(diag(song.ppca$Eval))*100,2)

song.massresid.ppca$L
song.ppca$L

song.massresid.ppca$lambda
song.ppca$lambda

### PCA Loading Table ###
write.xlsx(song.massresid.ppca$L[,1:6],file="/Users/NickMason/Desktop/Manuscripts/Oscinevssuboscinerates/Tables/MassresidPCLoadings_v5.xlsx")
write.xlsx(song.ppca$L[,1:6],file="/Users/NickMason/Desktop/Manuscripts/Oscinevssuboscinerates/Tables/PCLoadings_v5.xlsx")

### Setting up BAMM input files
is.ultrametric(mcc.trim)
is.binary.tree(mcc.trim)
min(mcc.trim$edge.length) ##all branches must be positive

### Writing out data files for each song character after extracting residuals
for(i in c(1,3:ncol(song_mass_resid))){
	sink(file=paste("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version3/",gsub("[.]","",colnames(song_mass_resid))[i],".txt",sep=""))
	for(j in 1:nrow(song_mass_resid)){
		cat(rownames(song_mass_resid)[j])
		cat("\\t")
		cat(scale(song_mass_resid)[j,i])
		cat("\\n")
	}
	sink()
	dir.create(paste("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version3/Output/",gsub("[.]","",colnames(song_mass_resid))[i],sep=""))### Create folders for output ###
}
	
### Writing out data files for PC1--PC6 for both massresid PCA and normal PCA
for(i in 1:6){
	sink(file=paste("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version3/massresidPC",i,".txt",sep=""))
	for(j in 1:nrow(song.massresid.ppca$S)){
		cat(rownames(song.massresid.ppca$S)[j])
		cat("\\t")
		cat(song.massresid.ppca$S[j,i])
		cat("\\n")
	}
	sink()
	dir.create(paste("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version3/Output/massresidPC",i,sep=""))### Create folders for output ###
}

### Writing out data files for each non-residual PCA axis 1 - 6
for(i in 1:6){
	sink(file=paste("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version3/normalPC",i,".txt",sep=""))
	for(j in 1:nrow(song.ppca$S)){
		cat(rownames(song.ppca $S)[j])
		cat("\\t")
		cat(song.ppca $S[j,i])
		cat("\\n")
	}
	sink()
	dir.create(paste("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version3/Output/normalPC",i,sep=""))### Create folders for output ###
}

### NEW AS OF JAN 2016 ###
### Generate appropriate BAMM priors for each trait ###
581/677#Setting for incomplete taxonomic variable in control files

setBAMMpriors(mcc.trim,total.taxa=677) #375 tanagers + 302 furnariids speciation extinction priors

files2<-list.files("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version3",full.names=T)
files2<-files2[grep("[.]txt",files2)]
files2_names<-gsub("[.]txt","",sapply(strsplit(files2,"/"),function(x) x[length(x)]))

for(i in 1:length(files2)){
	foo<-paste("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version3/setBAMMpriors_output/",files2_names[i],"BAMMpriors.txt",sep="")
	setBAMMpriors(mcc.trim,traits=files2[i],total.taxa=677,outfile=foo)
	bar<-readLines(file(foo))
	bar2<-bar[7:9]
	unlink(bar)
	foo2<-readLines(file(paste("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version3/ControlFiles/",files2_names[i],".txt",sep="")))
	foo2[71]<-bar2[1]
	foo2[75]<-bar2[3]
	write(foo2,file=paste("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version3/ControlFiles/",files2_names[i],".txt",sep=""))
}

### Write Ladderized TREE FOR V4 ### !!!!
write.tree(mcc.trim,file="~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version4/mcctrim.tre")