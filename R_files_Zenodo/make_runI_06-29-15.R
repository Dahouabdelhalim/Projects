#make_runI_10-22-14.R

#how many taxa

library(paleotree)

nChar<-54	#66 chars - 12 chars that aren't pars informative with only share taxa = 54 chars
nrun<-10

############################

#what tree to use? -let's use the majRule morph-only tree from PAUP (unnamed run)
setwd("C:/dave/workspace/mrbayes")
morphTree<-read.nexus("runYmajRule_forRunI_06-29-15.nex")

#This tree is already rooted correctly
plot(morphTree)

#####

#GET MOLECULAR MATRIX FROM RUN D

files<-list.files()    
nexFile<-grep(x=grep(x=files,pattern="^runD",value=TRUE),
	pattern=".nex$",value=TRUE)
origNex<-scan(nexFile,what="character",sep="\\n",blank.lines.skip=FALSE)
molMat<-origNex[grep(origNex,pattern="^Abyss"):(grep(origNex,pattern="^Tethy"))]
splitMolMat<-lapply(strsplit(molMat,"\\\\s+",fixed=FALSE),function(x)
	 c(x[1],paste0(x[2:length(x)],collapse="")))
#make sure they are all length 2
if(any(sapply(splitMolMat,length)!=2)){stop("WARNING WILL ROBINSON")}
#now seperate them
molTaxa<-sapply(splitMolMat,function(x) x[[1]])
molMat<-unlist(lapply(splitMolMat,function(x) x[[2]]))

###########

#get file for MrBayes run F (maximize information settings)
	#Settings identical to run F (maximize information, shared taxa)
nexFile<-grep(x=grep(x=files,pattern="^runF",value=TRUE),
	pattern=".nex$",value=TRUE)
origNex<-scan(nexFile,what="character",sep="\\n",blank.lines.skip=FALSE)

#####################

#identify droppers as taxon names

droppers<-strsplit(strsplit(strsplit(grep(origNex,
	pattern="taxset unsharedDelete =",value=TRUE),
	"=")[[1]][2],";")[[1]][1]," ")[[1]]
droppers<-as.numeric(droppers[nchar(droppers)>0])

#get taxon names
mat<-origNex[grep(origNex,pattern="^Abyss"):(grep(origNex,pattern="^Tethy"))]
taxaF<-sapply(strsplit(mat[-grep(mat,pattern="^\\t")]," "),function(x) x[1])
dropTaxa<-taxaF[droppers]

#drop from morphTree
morphTreeDrop<-drop.tip(morphTree,dropTaxa)

plot(morphTreeDrop)

#dummy matrix for adding to simMat later
dropMat<-matrix("?",length(dropTaxa),nChar)
rownames(dropMat)<-dropTaxa

####################

#modify run F file

#replace line defining number of characters
origNex[grep(origNex,pattern="^dimensions")]<-paste0(
	"dimensions ntax=",Ntip(morphTreeDrop)+length(dropTaxa),
	" nchar=",3435+nChar,";")

#replace format line
origNex[grep(origNex,pattern="^format")]<-paste0(
	"format datatype=mixed( dna: 1-3435, standard: 3436-",3435+nChar,
	") gap=- missing=?;")
#remove charlabels
#origNex<-origNex[-grep(origNex,pattern="charstatelabels")]

#edit 'DATA' section
	#remove all sections but taxon-dropping section
origNex<-origNex[-((grep(origNex,pattern="^.DATA")+1):
	(grep(origNex,pattern="^.delete all taxa")-1))]
origNex<-origNex[-((grep(origNex,pattern="means there are now")):
	(grep(origNex,pattern="\\t\\texclude")))]

#break into parts
nexStart<-origNex[2:grep(origNex,pattern="^matrix")]
nexGut<-origNex[(grep(origNex,pattern="^;")[1]):grep(origNex,pattern="mcmcp")]
nexEnd<-origNex[(1+grep(origNex,pattern="mcmcp")):length(origNex)]

#replace log start line in nexGut
nexGut[grep(nexGut,pattern="log start")]<-'log start filename="I_log.out" replace;'


bayesNex<-c("#NEXUS")
for(i in 1:nrun){
	simMorph<-perfectParsCharTree(morphTreeDrop,nchar=nChar)	
	#add dropped taxa back to simMat with ??? codings
		#rbind(rep(paste0(rep("?",nChar),collapse=""),length(dropTaxa))
	simMorph1<-rbind(simMorph,dropMat)
	#
	simMat1<-simMat<-simMorph1[order(rownames(simMorph1)),]
	colnames(simMat1)<-rownames(simMat1)<-NULL
	#need to make simMat look like char block from nexus file
	#
	#check order with molmat taxa
	if(!identical(row.names(simMat),molTaxa)){stop("molTaxa doesn't match simMat")}
	#
	#okay, now combine with molMat
	combMat<-paste0(molMat,apply(simMat1,1,paste0,collapse=""))
	#
	#first, appropriate # of spaces
	spaceBlock<-max(nchar(row.names(simMat))+10)-nchar(row.names(simMat))
	spaceBlock<-sapply(spaceBlock,function(x)
		paste0(rep(" ",x),collapse=""))
	#
	#next, make the matrix as char strings MWAHAHAH
	simMid<-paste0(row.names(simMat),spaceBlock,combMat)
	#
	#need new mcmcp line for filename
	newName<-paste0("runI_sim-",letters[i],"_rynch_10-24-14.nex")
	fileLine<-paste0("\\tmcmcp filename = ",newName," ; [new file name]")
	#
	#use log lines or not?
	nexGut1<-if(i==1){nexGut}else{nexGut[-grep(nexGut,pattern="log start")]}
	nexEnd1<-if(i==nrun){nexEnd}else{nexEnd[-grep(nexEnd,pattern="log stop")]}
	#
	#now use regexp to put it into .nex files for MrBayes
	newBayes<-c(nexStart,"","",simMid,"",nexGut1,fileLine,nexEnd1)
	bayesNex<-c(bayesNex,newBayes)
	}

write(bayesNex,file= "runI_rynchComb_06-29-14.nex")
