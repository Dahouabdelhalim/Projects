require(ape)
require(coda)
require(BAMMtools)
require(fBasics)
require(RColorBrewer)
require(plotrix)
require(phytools)
require(png)
require(raster)
require(xlsx)
require(shape)
require(aspace)
source("exponentialRate.R")
packageVersion("BAMMtools")

tree<-read.tree("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version4/mcctrim.tre")

outputdir_fix<-list.files("~/Desktop/Manuscripts/Oscinevssuboscinerates/BAMM/Version4/Output",full.names=T)
outputdir_fix<-outputdir_fix[-grep("Mass",outputdir_fix)] #Remove body mass
outputdir_fix<-outputdir_fix[-grep("normal",outputdir_fix)] #Remove normal PCAs

outputdir<-outputdir_fix[1]

bamm.output<-function(outputdir,tree){
	if(length(grep("Speciation",outputdir))==1){
		edata <- getEventData(tree, eventdata = paste(outputdir,"/event_data.txt",sep=""), burnin=0.1,type="diversification")
		}else{
			edata <- getEventData(tree, eventdata = paste(outputdir,"/event_data.txt",sep=""), burnin=0.1,type="trait")
		}
	cmat <- getCohortMatrix(edata)
	#foo.cohorts<-cohorts(cmat,edata,lwd=0.5)
	mcmcout <- read.csv(paste(outputdir,"/mcmc_out.txt",sep=""),header=T)
	#priorfile <- read.csv(paste(outputdir,"/prior_probs.txt",sep=""),header=T)
	burnstart <- floor(0.1 * nrow(mcmcout))
	postburn <- mcmcout[burnstart:nrow(mcmcout), ]
	post_probs <- table(postburn$N_shifts) / nrow(postburn)
	# bf.matrix<-computeBayesFactors(mcmcout,priorfile)
	# if("Inf" %in% bf.matrix[,1]){
		# bf.matrix<-bf.matrix[-which(bf.matrix[,1]=="Inf"),]
	# }
	ess.shifts<-effectiveSize(postburn$N_shifts)
	ess.ll<-effectiveSize(postburn$logLik)
	shift_probs <- summary(edata)
	#priorshifts<-getBranchShiftPriors(tree,2)
	
	sc <- distinctShiftConfigurations(edata, 2, threshold = 20)
	bestSamples <- sc$samplesets[[1]]
	bested <- subsetEventData(edata, index=bestSamples)
	best2<-getBestShiftConfiguration(edata,expectedNumberOfShifts=2,threshold=20)
	
	css <- credibleShiftSet(edata,2, threshold = 20, set.limit = 0.95)
	
	?credibleShiftSet
	class(css)<-"bammdata"
	msc.set <- maximumShiftCredibility(edata, maximize='product')
	msc.config <- subsetEventData(edata, index = msc.set$sampleindex)
	
	### Using Pascal's script ###
	subb <- subsetEventData(css, index = css$indices[[1]])
	coreshifts <- c((length(css$tip.label) + 1), css$coreshifts)
	coreshifts <- intersect(subb$eventData[[1]]$node, coreshifts)
	for (i in 1:length(subb$eventData)) {
		if (i == 1) {
			ff <- subb$eventData[[i]]
		}
		ff <- rbind(ff, subb$eventData[[i]])
	}
	xn <- numeric(length(coreshifts))
	xc <- character(length(coreshifts))
	dff <- data.frame(generation = xn, leftchild = xc, rightchild = xc, abstime = xn, betainit = xn, betashift = xn, stringsAsFactors = F)
	for (i in 1:length(coreshifts)) {
		if (coreshifts[i] <= length(css$tip.label)) {
			dset <- c(css$tip.label[coreshifts[i]], NA)
		} else {
		tmp <- extract.clade(as.phylo(css), node = coreshifts[i])
		dset <- tmp$tip.label[c(1, length(tmp$tip.label))]
		}
		tmp2 <- ff[ff$node == coreshifts[i], ]
		dff$leftchild[i] <- dset[1]
		dff$rightchild[i] <- dset[2]
		dff$abstime[i] <- mean(tmp2$time)
		dff$betainit[i] <- mean(tmp2$lam1)
		dff$betashift[i] <- mean(tmp2$lam2)
		}
	best_ed3 <- getEventData(as.phylo(css), eventdata = dff, type = "trait")
	
	bested3 <- best_ed3
	
	shiftnodes <- getShiftNodesFromIndex(bested3, 1)
	shiftnode_parents <- bested3 $edge[match(shiftnodes, bested3$edge[, 2], nomatch = 0), 1]
	root <- (shiftnode_parents == (bested3$Nnode + 2))
	if (sum(root) > 0) {
		isShiftNodeParent <- integer(length(shiftnodes))
		isShiftNodeParent[root] <- 1
		isShiftNodeParent[!root] <- bested3$eventVectors[[1]][match(shiftnode_parents[!root], bested3$edge[, 2])]
		} else {
		isShiftNodeParent <- bested3$eventVectors[[1]][match(shiftnode_parents, bested3$edge[, 2])]
		}
	isShiftNode <- match(shiftnodes, bested3$eventData[[1]]$node)
	time <- bested3$eventData[[1]][isShiftNode, 2] - bested3$eventData[[1]][isShiftNodeParent, 2]
	
	lam1 <- bested3$eventData[[1]][isShiftNodeParent, 3]
	lam2 <- bested3$eventData[[1]][isShiftNodeParent, 4]
	
	rateBeforeShift <- exponentialRate(time, lam1, lam2)
	rateAfterShift <- bested3$eventData[[1]][isShiftNode, 3]
	
	res <- cbind(bested3$eventData[[1]], beforeShift=c(NA, rateBeforeShift))
	
	colvec <- rep('#0000FF80', nrow(res) - 1)
	restmp <- res[complete.cases(res),]
	colvec[restmp$lam1 > restmp$beforeShift] <- '#FF000080'
	
	return(list(edata=edata,ess.shifts=ess.shifts,ess.ll=ess.ll,bested=best2,msc.config=msc.config,css=css,colvec=colvec))
}

output.list<-list()
for(i in length(output.list)+1:length(outputdir_fix)){
	if(i ==1){
		output.list[[i]]<-bamm.output(outputdir_fix[i],tree)
	}else{
		output.list[[i]]<-bamm.output(outputdir_fix[i],tree)
	}
}	

names(output.list)<-sapply(strsplit(outputdir_fix,"/"),function(x) x[length(x)])

save(output.list,file="output.list_v2.Rdata")
load("output.list_v2.Rdata")

### Identify taxa in speciation rate shifts ###
se_eventdata<-output.list[[14]]$bested$eventData[[1]]

taxa1<-extract.clade(tree,se_eventdata[4,1])$tip.label
taxa2<-extract.clade(tree,se_eventdata[5,1])$tip.label
taxa3<-extract.clade(tree,se_eventdata[6,1])$tip.label

taxa_list_bursts<-list(taxa1,taxa2,taxa3)

save(taxa_list_bursts,file="taxa_list_bursts.Rdata")

### Split output.list into two: 1. resid-corrected PCA; 2. resid-corrected characters
output.list.mPC<-c(14,grep("massresid",names(output.list)))
output.list.c<-c(14,(1:length(output.list))[!1:length(output.list) %in% c(output.list.mPC)])

names(output.list)[output.list.c] #Mass is not included here in v5

output.list.mPC<-output.list[output.list.mPC]
output.list.c<-output.list[output.list.c]

### Reorder items in output.list
names(output.list.c)
output.list.c<-output.list.c[c(1,3,2,6,7,8,4,5,9)]

## Rename characters to match rest of the manuscript
names(output.list.c)<-c("Speciation rate","Minimum frequency","Maximum frequency","Peak frequency","Song frequency range","Song length","Note count","Note rate","Vocal performance")

#names(output.list.mPC)<-c("Speciation Rate",paste("PCA",1:6))

#################
### BAMM Plots ###
#################
colpal<-rgb(colorRamp(c("black","purple4","blue","limegreen","yellow","orange","red"),space="rgb")(0:100/100),maxColorValue=255)
colpal2<-c("#0F0002","#3A013A","#100063","#005D8C","#00B75B","#4AE000","#FFD400","#FF0000")

#Read in PNGs for figure
sbill<-readPNG(source="~/Desktop/Manuscripts/Oscinevssuboscinerates/Plates/Scythebill_rot1.png")
tanchil<-readPNG(source="~/Desktop/Manuscripts/Oscinevssuboscinerates/Plates/TangaraChilensis.png")

## Figure out which nodes are part of diversification rate shifts to address reviewer comments

se_nodes<-output.list.mPC[[1]]$bested$eventData[[1]]$node[-1]
desc_nodes<-lapply(se_nodes, function(x) c(x,getDescendants(tree,x)))

for(i in 3:5){
	desc_nodes[[1]]<-desc_nodes[[1]][!desc_nodes[[1]] %in% desc_nodes[[i]]]
	desc_nodes[[2]]<-desc_nodes[[2]][!desc_nodes[[2]] %in% desc_nodes[[i]]]
}

plot(tree,show.tip.label=F,type="fan")
nodelabels(node=se_nodes)
layout.mat<-matrix(c(rep(1:4,each=2),8,rep(5:7,each=2),9,rep(10,8)),nrow=3,byrow=T)

pch_vec<-c(23,21,21,23,21) #bg furn, bg tan, spor, cran, DarFin
bor_vec<-c("#FF0000","#0000FF","#0000FF","#FF0000","#0000FF")
bg_vec<-c("#FF0000","#0000FF","#FFFF00","#FFFF00","#FF0000")
bg_vec<-paste(bg_vec,"90",sep="")
cex_vec<-c(0.6,0.6,1.2,1.2,1.2)
cex_text_vec<-c(0.3,0.3,0.6,0.6,0.6)

### PCA plot in the supplementary materials ###
names(output.list.mPC)[1] <- "Speciation rate"
names(output.list.mPC)<-gsub("massresid","",names(output.list.mPC))

pdf(width=6.5,height=6.5*(2/4)+(6.5*(2/4)*0.05),file="~/Desktop/Manuscripts/Oscinevssuboscinerates/Figures/PhyloRatePlots_massresid_v5.pdf")

#quartz(width=6.5,height=6.5*(2/4)+(6.5*(2/4)*0.05))

layout(layout.mat,heights=c(rep(0.95/2,2),0.05))
par(mar=c(0.5,0.5,0.5,0.5))

for(i in 1:length(output.list.mPC)){
	par(xpd=T)
	par(mar=c(0.5,0.5,0.5,0.5))
	plot.bammdata(output.list.mPC[[i]]$bested,pal=colpal,tau=0.01,labels=F,cex=0.1,method="polar",lwd=0.5)
	rasterImage(sbill,xleft=0.15,ybottom=-0.35,xright=0.55,ytop=0.15)
	rasterImage(tanchil,xleft=-0.55,ybottom=0.15,xright=-0.2,ytop=0.5)
	title(main=names(output.list.mPC)[i],cex.main=0.8,line=-0.1,font.main=1)
	
	if(length(output.list.mPC[[i]]$best$eventData[[1]]$node[-1])!=0){
		tf_vec1<-vector()
		for(j in 1:length(output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)])){
			tf_vec1<-c(tf_vec1,which(sapply(desc_nodes, function(x) output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j] %in% x)))
			if(tf_vec1[j] %in% c(1,2) & output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j] > Ntip(tree) ){
				if(any(getDescendants(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j]) %in% desc_nodes[[3]]) & Ntip(extract.clade(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j])) < 50){
					tf_vec1[j]<-3
				}
				if(any(getDescendants(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j]) %in% desc_nodes[[4]]) & Ntip(extract.clade(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j])) < 50){
					tf_vec1[j]<-4
				}
				if(any(getDescendants(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j]) %in% desc_nodes[[5]]) & Ntip(extract.clade(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j])) < 50){
					tf_vec1[j]<-5
					}
				}
			}
				
		nodelabels(node=output.list.mPC[[i]]$best$eventData[[1]]$node[-1],pch=pch_vec[tf_vec1],lwd=0.5,bg=bg_vec[tf_vec1],col=bor_vec[tf_vec1],cex=cex_vec[tf_vec1])
		nodelabels(text=1:length(output.list.mPC[[i]]$best$eventData[[1]]$node[-1]),node= output.list.mPC[[i]]$best$eventData[[1]]$node[-1],frame="none",cex=cex_text_vec[tf_vec1])
	}
	
	corner.label(label=LETTERS[i],figcorner=T,font=2,cex=0.8)
	if(i == 1){
		par(xpd=NA)		
		text(x=-0.42,y=0.14,"Thraupidae",cex=0.6,srt=326)
		text(x=0.30,y=-0.22,"Furnariidae",cex=0.6,srt=338)
		arctext("Cranioleuca",radius=1.15,clockwise=50,middle=2.45*pi/2,col='black',cex=0.35,font=3)
		arctext("Darwin's Finches",radius=1.25,clockwise=50,middle=0.25*pi/2,col='black',cex=0.35,font=1)
		arctext("Sporophilinae",radius=1.15,clockwise=50,middle=0.1*pi/2,col='black',cex=0.35,font=1)
	}
}

## Fill in empty plot slots ###
par(mar=c(0,0,0,0))
plot(0:10,0:10, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
plot(0:10,0:10, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
plot(0:10,0:10, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")

foo.ind<-seq(3,7,length.out=length(colpal)+1)
par(xpd=NA)

for(i in 1:length(colpal)){
	rect(foo.ind[i],6,foo.ind[i+1],11,border=NA,col=colpal[i])
}
text(x=5,y=13,label="Evolutionary Rates",cex=0.8,adj=c(0.5,0))
text(x=3,y=1,label="Low",cex=0.8,adj=c(0.5,0))
text(x=7,y=1,label="High",cex=0.8,adj=c(0.5,0))

dev.off()

### Song character (not PCA) phylo rate figures (i.e. Figure 1 in manuscript) ###
pdf(width=6.5,height=6.5+(6.5*0.05),file="~/Desktop/Manuscripts/Oscinevssuboscinerates/Figures/PhyloRatePlot_characters_v5.pdf")

#quartz(width=6.5,height=6.5+(6.5*0.05))
layout(matrix(c(1:9,rep(10,3)),nrow=4,byrow=T),heights=c(rep(0.95/3,3),0.05))
par(mar=c(0.5,0.5,0.5,0.5))

for(i in 1:length(output.list.c)){
	par(xpd=T)		
	par(mar=c(0.5,0.5,0.5,0.5))
	plot.bammdata(output.list.c[[i]]$best,pal=colpal,tau=0.01,labels=F,cex=0.1,method="polar",lwd=0.5)
	rasterImage(sbill,xleft=0.15,ybottom=-0.35,xright=0.55,ytop=0.15)
	rasterImage(tanchil,xleft=-0.55,ybottom=0.15,xright=-0.2,ytop=0.5)
	title(main=names(output.list.c)[i],cex.main=0.8,line=-0.1,font.main=1)

	if(length(output.list.c[[i]]$best$eventData[[1]]$node[-1])!=0){
		tf_vec1<-vector()
		for(j in 1:length(output.list.c[[i]]$best$eventData[[1]]$node[-(1)])){
			tf_vec1<-c(tf_vec1,which(sapply(desc_nodes, function(x) output.list.c[[i]]$best$eventData[[1]]$node[-(1)][j] %in% x)))
			if(tf_vec1[j] %in% c(1,2) & output.list.c[[i]]$best$eventData[[1]]$node[-(1)][j] > Ntip(tree) ){
				if(any(getDescendants(tree, output.list.c[[i]]$best$eventData[[1]]$node[-(1)][j]) %in% desc_nodes[[3]]) & Ntip(extract.clade(tree, output.list.c[[i]]$best$eventData[[1]]$node[-(1)][j])) < 50){
					tf_vec1[j]<-3
				}
				if(any(getDescendants(tree, output.list.c[[i]]$best$eventData[[1]]$node[-(1)][j]) %in% desc_nodes[[4]]) & Ntip(extract.clade(tree, output.list.c[[i]]$best$eventData[[1]]$node[-(1)][j])) < 50){
					tf_vec1[j]<-4
				}
				if(any(getDescendants(tree, output.list.c[[i]]$best$eventData[[1]]$node[-(1)][j]) %in% desc_nodes[[5]]) & Ntip(extract.clade(tree, output.list.c[[i]]$best$eventData[[1]]$node[-(1)][j])) < 50){
					tf_vec1[j]<-5
					}
				}
			}
		nodelabels(node= output.list.c[[i]]$best$eventData[[1]]$node[-1],pch=pch_vec[tf_vec1],lwd=0.5,bg=bg_vec[tf_vec1],col=bor_vec[tf_vec1],cex=cex_vec[tf_vec1])
		nodelabels(text=1:length(output.list.c[[i]]$best$eventData[[1]]$node[-1]),node= output.list.c[[i]]$best$eventData[[1]]$node[-1],frame="none",cex=cex_text_vec[tf_vec1])
	}	
	corner.label(label=LETTERS[i],figcorner=T,font=2,cex=0.8)
	
	#if(i==1){
	par(xpd=NA)		
	text(x=-0.42,y=0.14,"Thraupidae",cex=0.6,srt=326)
	text(x=0.30,y=-0.22,"Furnariidae",cex=0.6,srt=338)
	arctext("Cranioleuca",radius=1.125,clockwise=50,middle=2.45*pi/2,col='black',cex=0.35,font=3)
	arctext("Darwin's Finches",radius=1.225,clockwise=50,middle=0.28*pi/2,col='black',cex=0.35,font=1)
	arctext("Sporophilinae",radius=1.15,clockwise=50,middle=0.1*pi/2,col='black',cex=0.35,font=1)
	#}
}

par(mar=c(0,0,0,0))
plot(0:10,0:10, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")

foo.ind<-seq(3,7,length.out=length(colpal)+1)
par(xpd=NA)

for(i in 1:length(colpal)){
	rect(foo.ind[i],4,foo.ind[i+1],9,border=NA,col=colpal[i])
}
text(x=5,y=10,label="Evolutionary Rates",cex=0.8,adj=c(0.5,0))
text(x=3,y=1,label="Low",cex=0.8,adj=c(0.5,0))
text(x=7,y=1,label="High",cex=0.8,adj=c(0.5,0))

dev.off()

### Create extra supplementary figure that has speciation rate plot alone with tip labels ###
pdf(width=6.5,height=6.5,file="~/Desktop/Manuscripts/Oscinevssuboscinerates/Figures/BigTreeSupplementary_v1.pdf")

layout(matrix(c(1,2),nrow=2),heights=c(0.9,0.1))
par(mar=c(0.25,0.25,0.25,0.25))

for(i in 1:1){
	par(xpd=T)		
	plot.bammdata(output.list.c[[i]]$best,pal=colpal,tau=0.01,labels=T,cex=0.2,method="polar",lwd=1)
	rasterImage(sbill,xleft=0.15,ybottom=-0.35,xright=0.55,ytop=0.15)
	rasterImage(tanchil,xleft=-0.55,ybottom=0.15,xright=-0.2,ytop=0.5)
	#title(main=names(output.list.c)[i],cex.main=0.8,line=-0.1,font.main=1)

	if(length(output.list.mPC[[i]]$best$eventData[[1]]$node[-1])!=0){
		tf_vec1<-vector()
		for(j in 1:length(output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)])){
			tf_vec1<-c(tf_vec1,which(sapply(desc_nodes, function(x) output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j] %in% x)))
			if(tf_vec1[j] %in% c(1,2) & output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j] > Ntip(tree) ){
				if(any(getDescendants(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j]) %in% desc_nodes[[3]]) & Ntip(extract.clade(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j])) < 50){
					tf_vec1[j]<-3
				}
				if(any(getDescendants(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j]) %in% desc_nodes[[4]]) & Ntip(extract.clade(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j])) < 50){
					tf_vec1[j]<-4
				}
				if(any(getDescendants(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j]) %in% desc_nodes[[5]]) & Ntip(extract.clade(tree,output.list.mPC[[i]]$best$eventData[[1]]$node[-(1)][j])) < 50){
					tf_vec1[j]<-5
					}
				}
			}
		nodelabels(node= output.list.c[[i]]$best$eventData[[1]]$node[-1],pch=pch_vec[tf_vec1],lwd=1.25,cex=1.75,bg=bg_vec[tf_vec1], col=bor_vec[tf_vec1])
		nodelabels(text=1:length(output.list.c[[i]]$best$eventData[[1]]$node[-1]),node= output.list.c[[i]]$best$eventData[[1]]$node[-1],frame="none",cex=0.8)
	}	
	
	#if(i==1){
	par(xpd=NA)		
	text(x=-0.42,y=0.14,"Thraupidae",cex=1.2,srt=326)
	text(x=0.30,y=-0.22,"Furnariidae",cex=1.2,srt=338)
	arctext("Cranioleuca",radius=1.3,clockwise=50,middle=2.45*pi/2,col='black',cex=1,font=3)
	arctext("Darwin's Finches",radius=1.4,clockwise=50,middle=0.28*pi/2,col='black',cex=1,font=1)
	arctext("Sporophilinae",radius=1.3,clockwise=50,middle=0.1*pi/2,col='black',cex=1,font=1)
	#}
}

par(mar=c(0,0,0,0))
plot(0:10,0:10, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")

foo.ind<-seq(3,7,length.out=length(colpal)+1)
par(xpd=NA)

for(i in 1:length(colpal)){
	rect(foo.ind[i],4,foo.ind[i+1],9,border=NA,col=colpal[i])
}
text(x=5,y=10,label="Evolutionary Rates",cex=0.8,adj=c(0.5,0))
text(x=3,y=1,label="Low",cex=0.8,adj=c(0.5,0))
text(x=7,y=1,label="High",cex=0.8,adj=c(0.5,0))

dev.off()
	
### Barplots of shift probability scores for mass resid###
pdf(width=6.5,height=6.5*(2/4)+6.5*0.05,file="~/Desktop/Manuscripts/Oscinevssuboscinerates/Figures/NumberOfRegimes_Massresid_v3.pdf")
quartz(width=6.5,height=6.5*(2/4)+6.5*0.05)

layout(matrix(c(rep(1:4,each=2),5,rep(6:8,each=2),9,rep(10,8)),nrow=3,byrow=T),heights=c(rep(0.90/2,2),0.1))
for(i in 1:length(output.list.mPC)){
	par(mar=c(1.5,2,0.5,0))
	foo<-rep(0,100)
	names(foo)<-1:100

	post_prob<-summary(output.list.mPC[[i]]$edata)
	foo[post_prob[,1]+1]<-post_prob[,2]
	
	bp<-barplot(foo,names.arg=F,axes=F,ylim=c(0,(range(post_prob[,2])[2]-range(post_prob[,2])[1])*1.15),col=colpal[round(seq(from=1,to=length(colpal),length.out=100))],border=colpal[round(seq(from=1,to=length(colpal),length.out=100))])
	
	axis(1,at=c(0,bp[seq(0,100,5),1]),labels=seq(0,100,5),cex.axis=0.5,padj=-3.5,tck=-0.025)
	axis(2,cex.axis=0.5,padj=2.5, tck=-0.025)
	mtext("Number of Evolutionary Regimes",side=1,cex=0.5,line=0.65)
	mtext("Frequency",side=2,cex=0.5,line=1)
	box()
	
	corner.label(names(output.list.mPC)[i],x=1,y=1,cex=0.8)
	corner.label(LETTERS[i],figcorner=F,font=2,cex=0.8)
	par(xpd=NA)
	
	if(i %in% c(4,7)){
		par(mar=c(0,0,0,0))
		plot(0:10,0:10, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
	}
}
par(mar=c(0,0,0,0))
plot(0:10,0:10, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")

foo.ind<-seq(3,7,length.out=length(colpal)+1)
par(xpd=NA)

for(i in 1:length(colpal)){
	rect(foo.ind[i],3,foo.ind[i+1],5,border=NA,col=colpal[i])
}
text(x=5,y=6,label="Evolutionary Heterogeneity",cex=0.8,adj=c(0.5,0))
text(x=3,y=0.5,label="Low",cex=0.8,adj=c(0.5,0))
text(x=7,y=0.5,label="High",cex=0.8,adj=c(0.5,0))

dev.off()

### Barplots of shift probability scores for non-PCA characters###
pdf(width=6.5,height=6.5+6.5*0.05,file="~/Desktop/Manuscripts/Oscinevssuboscinerates/Figures/NumberOfRegimes_characters_v3.pdf")
quartz(width=6.5,height=6.5+6.5*0.05)

layout(matrix(c(1:9,rep(10,3)),nrow=4,byrow=T),heights=c(rep(0.95/3,3),0.05))
for(i in 1:length(output.list.c)){
	par(mar=c(1.5,2,0.5,0))
	foo<-rep(0,100)
	names(foo)<-1:100

	post_prob<-summary(output.list.c[[i]]$edata)
	foo[post_prob[,1]+1]<-post_prob[,2]
	bp<-barplot(foo,names.arg=F,axes=F,ylim=c(0,(range(post_prob[,2])[2]-range(post_prob[,2])[1])*1.15),col=colpal[round(seq(from=1,to=length(colpal),length.out=100))],border=colpal[round(seq(from=1,to=length(colpal),length.out= 100))])
	
	axis(1,at=c(0,bp[seq(0,100,5),1]),labels=seq(0,100,5),cex.axis=0.5,padj=-3.5,tck=-0.025)
	axis(2,cex.axis=0.5,padj=2.5, tck=-0.025)
	mtext("Number of Evolutionary Regimes",side=1,cex=0.5,line=0.65)
	mtext("Frequency",side=2,cex=0.5,line=1)
	box()
	
	corner.label(names(output.list.c)[i],x=1,y=1,cex=0.8)
	corner.label(LETTERS[i],figcorner=F,font=2,cex=0.8)
	par(xpd=NA)
	
	# if(i %in% c(4,7)){
		# par(mar=c(0,0,0,0))
		# plot(0:10,0:10, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
	# }
}
par(mar=c(0,0,0,0))
plot(0:10,0:10, type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")

foo.ind<-seq(3,7,length.out=length(colpal)+1)
par(xpd=NA)

for(i in 1:length(colpal)){
	rect(foo.ind[i],3,foo.ind[i+1],5,border=NA,col=colpal[i])
}
text(x=5,y=6,label="Evolutionary Heterogeneity",cex=0.8,adj=c(0.5,0))
text(x=3,y=0.5,label="Low",cex=0.8,adj=c(0.5,0))
text(x=7,y=0.5,label="High",cex=0.8,adj=c(0.5,0))

dev.off()

#############################
##### BAMM PLOT SET UP #####
#############################
###########################################
### Compare diversification rates with other rates mPC rates !###
### STRAPP 7/14/15 ###
tip.rates<-list()
for(i in 1:length(output.list.mPC)){
	foo<-credibleShiftSet(output.list.mPC[[i]]$edata,2,threshold=20)
	class(foo)<-"bammdata"
	if(i==1){
		tip.rates[[i]]<-getTipRates(foo)$lambda.avg
	}else{
		tip.rates[[i]]<-getTipRates(foo)$beta.avg
	}
}
i<-1
names(output.list.mPC[[1]])

### Create data frame of tip rates###
tr.df.all<-as.data.frame(tip.rates)

### RUNN STRAPP for mass PCs
#?traitDependentBAMM

speciationSet<-credibleShiftSet(output.list.mPC[[1]]$edata,2,threshold=20)
class(speciationSet)<-"bammdata"

strapp.list<-list()
for(i in 2:ncol(tr.df.all)){
	print(i-1)
	foo<-tr.df.all[,i]
	names(foo)<-row.names(tr.df.all)
	strapp.list[[i-1]]<-traitDependentBAMM(output.list.mPC[[1]]$bested, foo,reps=1000,method="s",two.tailed=F,traitorder="p",rate='speciation')
}
strapp.df<-matrix(unlist(strapp.list),nrow=ncol(tr.df.all)-1,byrow=T)
strapp.df

# Create color vector for points
furn.bg<-extract.clade(tree, 583)$tip.label
tan.bg<-extract.clade(tree, 858)$tip.label

point.labs<-rbind(cbind(furn.bg,"furn.bg"),cbind(tan.bg,"tan.bg"))
rownames(point.labs)<-point.labs[,1]
point.labs<-point.labs[,-1]

extract.clade(tree, output.list.mPC[[1]]$best$eventData[[1]][2,1])$tip.label

#point.labs[!names(point.labs) %in% c(extract.clade(tree, output.list.mPC[[1]]$best$eventData[[1]][2,1])$tip.label,tan.bg)]<-"geo/scleur"

point.labs[extract.clade(tree, output.list.mPC[[1]]$best$eventData[[1]][4,1])$tip.label]<-"sporo"
point.labs[extract.clade(tree, output.list.mPC[[1]]$best$eventData[[1]][5,1])$tip.label]<-"cranio"
point.labs[extract.clade(tree, output.list.mPC[[1]]$best$eventData[[1]][6,1])$tip.label]<-"darwins"

point.labs<-as.factor(point.labs[rownames(tr.df.all)])
point.labs<-factor(point.labs,levels=levels(point.labs)[c(5,4,2,3,1)])

point.pch<-c(21,21,21,23,23)
border.col<-c("#0000FF90","#0000FF90","#0000FF90","#FF000090","#FF000090")
inside.col<-c("#0000FF90","#FFFF0090","#FF000090","#FF000090","#FFFF0090")

### DO IT FOR JUST THRAUPIDAE ###
tan.output.list.mPC<-list()
tip.rates<-list()
for(i in 1:length(output.list.mPC)){
	tan.output.list.mPC[[i]]<-subtreeBAMM(output.list.mPC[[i]]$edata,tips=tan.bg)
	class(tan.output.list.mPC[[i]])<-"bammdata"
	foo<-credibleShiftSet(tan.output.list.mPC[[i]],2,threshold=20)
	class(foo)<-"bammdata"
	if(i==1){
		tip.rates[[i]]<-getTipRates(tan.output.list.mPC[[i]])$lambda.avg
	}else{
		tip.rates[[i]]<-getTipRates(tan.output.list.mPC[[i]])$beta.avg
	}
}
sapply(tip.rates, length)

### Create data frame of tip rates###
tr.df.t<-as.data.frame(tip.rates)

speciationSet<-credibleShiftSet(tan.output.list.mPC[[1]],2,threshold=20)
class(speciationSet)<-"bammdata"

strapp.list<-list()
for(i in 2:ncol(tr.df.t)){
	print(i-1)
	foo<-tr.df.t[,i]
	names(foo)<-row.names(tr.df.t)
	strapp.list[[i-1]]<-traitDependentBAMM(speciationSet, foo,reps=1000,method="s",two.tailed=F,traitorder="p")
}
strapp.df.tan<-matrix(unlist(strapp.list),nrow=ncol(tr.df.t)-1,byrow=T)
strapp.df.tan

### DO IT FOR JUST FURNARIIDAE ###
furn.output.list.mPC<-list()
for(i in 1:length(output.list.mPC)){
	furn.output.list.mPC[[i]]<-subtreeBAMM(output.list.mPC[[i]]$css,tips=furn.bg)
	class(furn.output.list.mPC[[i]])<-"bammdata"
	if(i==1){
		tip.rates[[i]]<-getTipRates(furn.output.list.mPC[[i]])$lambda.avg
	}else{
		tip.rates[[i]]<-getTipRates(furn.output.list.mPC[[i]])$beta.avg
	}
}

### Create data frame of tip rates###
tr.df.f<-as.data.frame(tip.rates)

### STRAPP 7/14/15 ###
strapp.list<-list()
for(i in 2:ncol(tr.df.f)){
	print(i-1)
	foo<-tr.df.f[,i]
	names(foo)<-row.names(tr.df.f)
	strapp.list[[i-1]]<-traitDependentBAMM(furn.output.list.mPC[[1]], foo,reps=1000,method="s",two.tailed=F,traitorder="p")
}
strapp.df.furn<-matrix(unlist(strapp.list),nrow=ncol(tr.df.f)-1,byrow=T)
strapp.df.furn

############################################################
### Scatter plots of lineage diversification vs mPC rates  ###
############################################################
pdf(file="~/Desktop/Manuscripts/Oscinevssuboscinerates/Figures/DivRateScattPlotRates_massresid_v4.pdf",width=6.5,height=6.5*(2/4))
quartz(width=6.5,height=6.5*(2/4))
par(xpd=F)
layout(matrix(c(rep(1:4,each=2),5,rep(6:8,each=2),9),nrow=2,byrow=T))
par(mar=c(2,2.5,1,0.5))

for(i in 1:(length(output.list.mPC)-1)){
	x.label<-names(output.list.mPC)[i+1]
	if(as.numeric(strapp.df[,2][i])<=0.05){
	plot(log(tr.df.all[,i+1]),log(tr.df.all[,1]),xlab=x.label,ylab="log(Speciation Rate)",axes=F,pch=point.pch[as.numeric(point.labs)],col=border.col[as.numeric(point.labs)],bg=inside.col[as.numeric(point.labs)],ylim=c(range(log(tr.df.all[,1]))[1],diff(range(log(tr.df.all[,1])))*0.3+range(log(tr.df.all[,1]))[2]))
		l_line<-loess(log(tr.df.all[,1])~log(tr.df.all[,i+1]),span=diff(range(log(tr.df.all[,i+1])))/2)
		Q<-order(tr.df.all[,i+1])
		lines(log(tr.df.all[,i+1])[Q],l_line$fitted[Q],col="black",lty=1,lwd=2)	
	}else{
	plot(tr.df.all[,i+1],log(tr.df.all[,1]),xlab=x.label,ylab="log(Speciation Rate)",axes=F,pch=point.pch[as.numeric(point.labs)],col=border.col[as.numeric(point.labs)],bg=inside.col[as.numeric(point.labs)],ylim=c(range(log(tr.df.all[,1]))[1],diff(range(log(tr.df.all[,1])))*0.3+range(log(tr.df.all[,1]))[2]))
		}
	axis(1,padj=-3,cex.axis=0.6,tck=-0.025)
	mtext(paste("log(",names(output.list.mPC)[i+1]," rate)",sep=""),1,line=0.75,cex=0.6)
	axis(2,padj=1.5,cex.axis=0.6,tck=-0.025)
	mtext("log(Speciation rate)",2,line=1.2,cex=0.6)
	box()
	corner.label(LETTERS[[i]],figcorner=F,font=2)
	corner.label(bquote(paste(rho," = ",.(round(as.numeric(strapp.df[i,1]),2)),"; P = ",.(round(as.numeric(strapp.df[i,2]),3)))),x=1,cex=0.8)
	corner.label(bquote(paste(rho," = ",.(round(as.numeric(strapp.df.tan[i,1]),2)),"; P = ",.(round(as.numeric(strapp.df.tan[i,2]),3)))),x=1,cex=0.8,col='blue',yoff=0.3)
	corner.label(bquote(paste(rho," = ",.(round(as.numeric(strapp.df.furn[i,1]),2)),"; P = ",.(round(as.numeric(strapp.df.furn[i,2]),3)))),x=1,cex=0.8,col='red',yoff=0.5)
		if(i %in% c(4,7)){
		plot(1:10,1:10,axes=F,col="white")
	}
	}

plot(1:10,1:10,axes=F,col="white")
points(rep(1,6),c(9:7,4:2),pch=point.pch[1:6],col=border.col[1:6],bg=inside.col[1:6])
text(rep(1.5,6),c(9:7,4:2),labels=c("Thraupidae background","Sporophilinae","Darwin's Finches","Furnariidae background","Cranioleuca",""),cex=0.8,adj=c(0,0.5))

dev.off()

##############################################################
### Jacknifing mPC dataset to see if results hold up with resampling taxa ###
##############################################################
jack.mPC.tip.rates<-list()
for(i in 1:length(output.list.mPC)){
	foo<-credibleShiftSet(output.list.mPC[[i]]$edata,2,threshold=20)
	class(foo)<-"bammdata"
	if(i==1){
		jack.mPC.tip.rates[[i]]<-getTipRates(foo)$lambda.avg
	}else{
		jack.mPC.tip.rates[[i]]<-getTipRates(foo)$beta.avg
	}
}

### Create data frame of tip rates###
jack.mPC.tr.df.all<-as.data.frame(jack.mPC.tip.rates)
jack.mPC.array<-replicate(100,sample(1:nrow(jack.mPC.tr.df.all),(nrow(jack.mPC.tr.df.all)/2)),)

test<-output.list.mPC[[1]]$bested
class(test)<-"phylo"
test.trim<-drop.tip(test,output.list.mPC[[1]]$edata$tip.label[jack.mPC.array[,i]])

str(test.trim)

output.list.mPC[[1]]$edata

### RUNN STRAPP for mass PCs
#?traitDependentBAMM

speciationSet<-credibleShiftSet(output.list.mPC[[1]]$edata,2,threshold=20)
class(speciationSet)<-"bammdata"

strapp.list<-list()
for(i in 2:ncol(tr.df.all)){
	print(i-1)
	foo<-tr.df.all[,i]
	names(foo)<-row.names(tr.df.all)
	strapp.list[[i-1]]<-traitDependentBAMM(output.list.mPC[[1]]$bested, foo,reps=1000,method="s",two.tailed=F,traitorder="p",rate='speciation')
}
strapp.df<-matrix(unlist(strapp.list),nrow=ncol(tr.df.all)-1,byrow=T)
strapp.df

#######################################################
### Compare diversification rates with other CHARACTER rates ###
#######################################################
tip.rates<-list()
for(i in 1:length(output.list.c)){
	foo<-credibleShiftSet(output.list.c[[i]]$edata,20)
	class(foo)<-"bammdata"
	if(i==1){
		tip.rates[[i]]<-getTipRates(foo)$lambda.avg
	}else{
		tip.rates[[i]]<-getTipRates(foo)$beta.avg
	}
}

### Create data frame of tip rates###
tr.df.all<-as.data.frame(tip.rates)

strapp.list<-list()
for(i in 2:ncol(tr.df.all)){
	print(i-1)
	foo<-tr.df.all[,i]
	names(foo)<-row.names(tr.df.all)
	strapp.list[[i-1]]<-traitDependentBAMM(output.list.c[[1]]$bested, foo,reps=1000,method="s",two.tailed=F,traitorder="p")
}
strapp.df<-matrix(unlist(strapp.list),nrow=ncol(tr.df.all)-1,byrow=T)
strapp.df

### DO IT FOR JUST THRAUPIDAE ###
tan.output.list.c<-list()
for(i in 1:length(output.list.c)){
	tan.output.list.c[[i]]<-subtreeBAMM(output.list.c[[i]]$edata,tips=tan.bg)
	class(tan.output.list.c[[i]])<-"bammdata"
	if(i==1){
		tip.rates[[i]]<-getTipRates(tan.output.list.c[[i]])$lambda.avg
	}else{
		tip.rates[[i]]<-getTipRates(tan.output.list.c[[i]])$beta.avg
	}
}

### Create data frame of tip rates###
tr.df.t<-as.data.frame(tip.rates)

### STRAPP 7/14/15 ###
strapp.list<-list()
for(i in 2:ncol(tr.df.t)){
	print(i-1)
	foo<-tr.df.t[,i]
	names(foo)<-row.names(tr.df.t)
	strapp.list[[i-1]]<-traitDependentBAMM(tan.output.list.c[[1]], foo,reps=1000,method="s",two.tailed=F,traitorder="p")
}
strapp.df.tan<-matrix(unlist(strapp.list),nrow=ncol(tr.df.t)-1,byrow=T)
strapp.df.tan

### DO IT FOR JUST FURNARIIDAE ###
furn.output.list.c<-list()
for(i in 1:length(output.list.c)){
	furn.output.list.c[[i]]<-subtreeBAMM(output.list.c[[i]]$edata,tips=furn.bg)
	class(furn.output.list.c[[i]])<-"bammdata"
	if(i==1){
		tip.rates[[i]]<-getTipRates(furn.output.list.c[[i]])$lambda.avg
	}else{
		tip.rates[[i]]<-getTipRates(furn.output.list.c[[i]])$beta.avg
	}
}

### Create data frame of tip rates###
tr.df.f<-as.data.frame(tip.rates)

### STRAPP 7/14/15 ###
strapp.list<-list()
for(i in 2:ncol(tr.df.f)){
	print(i-1)
	foo<-tr.df.f[,i]
	names(foo)<-row.names(tr.df.f)
	strapp.list[[i-1]]<-traitDependentBAMM(furn.output.list.c[[1]], foo,reps=1000,method="s",two.tailed=F,traitorder="p")
}
strapp.df.furn<-matrix(unlist(strapp.list),nrow=ncol(tr.df.f)-1,byrow=T)
strapp.df.furn

############################################################
### Scatter plots of lineage diversification vs individual character rates ###
############################################################
# Generate scatter plots of PGLS Rates #
pdf(file="~/Desktop/Manuscripts/Oscinevssuboscinerates/Figures/DivRateScattPlotRates_characters_v3.pdf",width=6.5,height=6.5*(2/4)+(6.5*(2/4)*0.1))
#quartz(width=6.5,height=6.5*(2/4)+(6.5*(2/4)*0.15))
par(xpd=F)
layout(matrix(c(1:8,rep(9,4)),nrow=3,byrow=T),heights=c(rep(0.90/2,2),0.1))
par(mar=c(2,2.5,1,0.5))

for(i in 1:(length(output.list.c)-1)){
	x.label<-names(output.list.c)[i+1]
	if(as.numeric(strapp.df[,2][i])<=0.05){
	plot(log(tr.df.all[,i+1]),log(tr.df.all[,1]),xlab=x.label,ylab="log(Speciation Rate)",axes=F,pch=point.pch[as.numeric(point.labs)],col=border.col[as.numeric(point.labs)],bg=inside.col[as.numeric(point.labs)],ylim=c(range(log(tr.df.all[,1]))[1],diff(range(log(tr.df.all[,1])))*0.3+range(log(tr.df.all[,1]))[2]))
		l_line<-loess(log(tr.df.all[,1])~log(tr.df.all[,i+1]),span=diff(range(log(tr.df.all[,i+1])))/2)
		Q<-order(tr.df.all[,i+1])
		lines(log(tr.df.all[,i+1])[Q],l_line$fitted[Q],col="black",lty=1,lwd=2)	
	}else{
	plot(tr.df.all[,i+1],log(tr.df.all[,1]),xlab=x.label,ylab="log(Speciation Rate)",axes=F,pch=point.pch[as.numeric(point.labs)],col=border.col[as.numeric(point.labs)],bg=inside.col[as.numeric(point.labs)],ylim=c(range(log(tr.df.all[,1]))[1],diff(range(log(tr.df.all[,1])))*0.3+range(log(tr.df.all[,1]))[2]))
		}
	axis(1,padj=-3,cex.axis=0.6,tck=-0.025)
	mtext(paste("log(",names(output.list.c)[i+1]," rate)",sep=""),1,line=0.75,cex=0.6)
	axis(2,padj=1.5,cex.axis=0.6,tck=-0.025)
	mtext("log(Speciation rate)",2,line=1.2,cex=0.6)
	box()
	corner.label(LETTERS[[i]],figcorner=F,font=2)
	corner.label(bquote(paste(rho," = ",.(round(as.numeric(strapp.df[i,1]),2)),"; P = ",.(round(as.numeric(strapp.df[i,2]),3)))),x=1,cex=0.8)
	corner.label(bquote(paste(rho," = ",.(round(as.numeric(strapp.df.tan[i,1]),2)),"; P = ",.(round(as.numeric(strapp.df.tan[i,2]),3)))),x=1,cex=0.8,col='blue',yoff=0.3)
	corner.label(bquote(paste(rho," = ",.(round(as.numeric(strapp.df.furn[i,1]),2)),"; P = ",.(round(as.numeric(strapp.df.furn[i,2]),3)))),x=1,cex=0.8,col='red',yoff=0.5)
	}

par(mar=c(0,0,0,0))
plot(1:10,1:10,axes=F,col="white")
points(rep(seq(1,9,length.out=3),2),c(8,8,8,4,4,4),pch=point.pch[1:6],col=border.col[1:6],bg=inside.col[1:6])
text(rep(seq(1,9,length.out=3),2)+0.1,c(8,8,8,4,4,4),labels=c("Thraupidae background","Sporophilinae","Darwin's Finches","Furnariidae background","Cranioleuca",""),cex=0.8,adj=c(0,0.5))

dev.off()

###
names(output.list)

