### Dean Adams Multivariate Method ###
require(ape)
require(xlsx)
require(phytools)
require(BAMMtools)
source("AdamsMultivariateScript.R")
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
rownames(song)<-gsub("'","",row.names(song))
rownames(song)<-gsub("Syndactyla_dimidiatum","Syndactyla_dimidiata",row.names(song))
rownames(song)<-gsub("Certhiasomus_stictolaemus","Deconychura_stictolaema",row.names(song))

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

song.trim<-song[rownames(song) %in% rownames(mass),]

### combine mass and song dataset ###
song_mass<-merge(mass,song,by=0)
rownames(song_mass)<-song_mass$Row.names
song_mass$Row.names<-NULL

### Read in tree ###
mcc<-read.nexus("~/Desktop/Manuscripts/Oscinevssuboscinerates/Trees/ThraupidFurnariidSupertreeTrimmed.nex")
mcc.trim<-drop.tip(mcc,mcc$tip.label[!mcc$tip.label %in% rownames(song_mass)])
colnames(song_mass)

### Number of songs, taxa, sd, etc
Ntip(mcc.trim)

### Try ADAMS SCRIPT 2014 MULTIVARIATE ###
song_data<-song_mass[-c(1,2)]

### Log transform frequency units ###
colnames(song_data)
log.vec<-c(1,3,6,7)
for(i in 1:length(log.vec)){
	song_data[,log.vec[i]]<-log(song_data[,log.vec[i]])
}

group<-rep(NA,Ntip(mcc.trim))
names(group)<-mcc.trim$tip.label

group[extract.clade(mcc.trim,583)$tip.label]<-"furn"
group[extract.clade(mcc.trim,858)$tip.label]<-"thraup"

output<-CompareRates.sigma.d(mcc.trim,song_data,group)

output

### Jacknife to remove groups that are part of the different evolutionary regimes ###

### 100 random jackknifing replicates ###
### Create list of 100 randomizations of furanriid and thraupid taxa ###
furn_group<-extract.clade(mcc.trim,583)$tip.label
thraup_group<-extract.clade(mcc.trim,858)$tip.label

samp_list<-list()
for(i in 1:500){
	samp_list[[i]]<-furn_group[sample(1:length(furn_group),floor(length(furn_group)/10))]
	samp_list[[i]]<-c(samp_list[[i]],thraup_group[sample(1:length(thraup_group),floor(length(thraup_group)/10))])
	}
samp_list[1:2]

length(sample(1:length(furn_group),floor(length(furn_group)/10)))
length(sample(1:length(thraup_group),floor(length(thraup_group)/10)))

#rando_jk_output_10perc<-list()
for(i in (length(rando_jk_output_10perc)+1):length(samp_list)){
	print(paste("JackKNIFE Replicate!! #",i,sep="  :  "))
	mcc.foo<-drop.tip(mcc.trim,samp_list[[i]])
	group_vec<-rep("furn",Ntip(mcc.foo))
	group_vec[mcc.foo$tip.label %in% thraup_group]<-"thraup"
	names(group_vec)<-mcc.foo$tip.label
	song_data_foo<-song_data[mcc.foo$tip.label,]
	rando_jk_output_10perc[[i]]<-CompareRates.sigma.d(mcc.foo,song_data_foo,group_vec)
}


save(rando_jk_output_10perc,file="rando_jk_output_10perc.Rdata")
sigma.ratio<-sapply(rando_jk_output_10perc,function(x) x$sigma.d.all[2]/x$sigma.d.all[1])

quartz(width=6.5,height=6.5)
par(mar=c(5,5,0,0))
hist(log10(sigma.ratio[1:100]),breaks=20,main="",xlab="",ylab="")
abline(v=mean(log10(sigma.ratio[1:100])),col="red",lty=2,lwd=3)

sigma.jk.mean<-mean(sigma.ratio[1:100]) #Statistic to report in manuscript
c(sigma.jk.mean - var(sigma.ratio[1:100]) * 1.96,sigma.jk.mean + var(sigma.ratio[1:100]) * 1.96) #95% confidence intervals

sig_vec<-sapply(rando_jk_output_10perc,function(x) x$significance)
rando_jk_output_10perc[sig_vec>0.05]

### Remove the clades that have the highest speciation rates ###
load("taxa_list_bursts.Rdata")

taxa_list_bursts
mcc.bursts<-drop.tip(mcc.trim,unlist(taxa_list_bursts))
group_vec_bursts<-rep("furn",Ntip(mcc.bursts))
group_vec_bursts[mcc.bursts$tip.label %in% thraup_group]<-"thraup"
names(group_vec_bursts)<-mcc.bursts$tip.label
song_data_bursts<-song_data[mcc.bursts$tip.label,]

burst_jk_output<-CompareRates.sigma.d(mcc.bursts,song_data_bursts,group_vec_bursts)

plot(mcc.bursts,cex=0.3)
axisPhylo()