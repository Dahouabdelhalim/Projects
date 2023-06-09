## script for counting the number of morphological characters in a combined matrix that are parsimony informative

morphMat<-read.fwf(file.choose(),widths=rep(1,66))

## FULL MATRIX

#any invariant?
inVar<-which(apply(morphMat,2,function(x){
	xx<-x[x!="?"]	#remove question marks
	length(unique(xx))==1
	}))

3435+inVar


notParsInf<-which(apply(morphMat,2,function(x){
	xx<-x[x!="?"]	#remove question marks
	(sum(sapply(unique(xx),function(z) sum(z==xx)>1))<2 &
		length(unique(xx))<3)
	}))

3435+notParsInf

# 3440 3449 3458 3466 3475 3479 3488

## SHARED TAXA ONLY

morphMatSh<-morphMat[-c(1,3,7,11,13,21,22,23),]

#any invariant?
inVar<-which(apply(morphMatSh,2,function(x){
	xx<-x[x!="?"]	#remove question marks
	length(unique(xx))==1
	}))

inVar
3435+inVar


notParsInf<-which(apply(morphMatSh,2,function(x){

	#xx<-x[x!="?"]	#remove question marks
	#sum(sapply(unique(xx),function(z) sum(z==xx)>1))<2
	#}))

	xx<-x[x!="?"]	#remove question marks
	(sum(sapply(unique(xx),function(z) sum(z==xx)>1))<2 &
		length(unique(xx))<3)
	}))

notParsInf
3435+notParsInf












