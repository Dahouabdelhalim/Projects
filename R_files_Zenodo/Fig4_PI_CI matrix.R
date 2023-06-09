library(RColorBrewer)
require("gplots")
require("spatstat")
require("ggplot2")
library(reshape2)      
library(grid)  

######################################################################
#import data sumarizing mutated genes found in founder clones
folder=c("~/Dropbox/R/Wunsche/") 
data.file<-c("Fig4_parallelism.txt") 
open.file<-paste(folder,data.file,sep="")
data=read.table(open.file, header=T)

#Genes mutated relative to ancestor in founder clones
data.file<-c("Fig4_mutations_progenitor_also_in_replay.txt") 
open.file<-paste(folder,data.file,sep="")
data=read.table(open.file, header=T)
mut=as.data.frame(data)


#For evolvability data
data.file<-c("Fig4_replay_clone_mutations.txt") 
open.file<-paste(folder,data.file,sep="")
data=read.table(open.file, header=T, skip=1)
#remove genes with preexisting mutations in some marker variants
data=data[,-which(colnames(data)=="ara")]
data=data[,-which(colnames(data)=="eutB")]
data=data[,-which(colnames(data)=="ompR")]


#List of genes with mutations in early founders that have fixed mutations inlater founders
exclude_mut=c("pykF", "malT", "spoT", "kup.InsJ.5", "rbsD.yieO")

#Set global variable for exclusion of above genes:
remove.mutations = TRUE


#######################################################################
#Plot mutations per founder line and averages

dataP=data
#use this to control order genotypes plotted.
genotypes=c("anc","r","rt","rts","rtsg","rtsgp","K3","K10")
geno.lab=c("anc",expression(Ev^1),expression(Ev^2),expression(Ev^3),expression(Ev^4),expression(Ev^5),expression(Ev^9),expression(Ev^25))

#There are differences in mutation number by founder genotype, but does this correlate with initial fitness 
fitness=c(1,1.012,1.103,1.197,1.228,1.328,1.397,1.551) #in order: "Anc","r","rt","rts","rtsg","rtsgp","3K","10K"


#####################################################################
#convergence and parallelism indices of replicate populations
#After Kryazhimskiy et al. Science 2014

#mutation matrix -- evolved lines by mutations
m = data[,-which(names(data) %in% c("gene","strain"))]

#Observed CI
CIindex = function (dataM){
	#number of lines
	N = nrow(dataM)
	#vector to store ij sums
	CIij = vector(length = (N-1))
	#transpose mut matrix - makes much faster!
	dataT=t(dataM)
	for (i in 1:(N-1)){
		ith = dataT[,i]
		#vector to store ij values
		temp2 = vector(length = (N-i))
		for(j in (i+1):N){
			jth = dataT[,j] 
			temp = ith + jth
			ij = length(na.omit(temp))
						#add ij to a vector
			temp2[j-i] = ij
		}
	CIij[i] = sum(temp2)	
	} 
	return(sum(CIij))
}
#1/denominator number of potential interactions
rows = nrow(m)
df = rows*(rows-1) 

#alter df to account for some mutations occurring in early founders already being present in later founders
#these should be removed as a potential convergence event. 
#for each later occurring mutation, number of lines it occurs in/mutation rows * number of lines that already have it
p = 20/75 * 33
t = 1/75 * 57
y = 1/75 * 11
ma = 31/75 * 11
s = 3/75 * 52
k = 3/75 * 11
r = 8/75 * 66
h = 1/75 * 23
a = p+t+y+ma+s+k+r+h
CI = CIindex(dataM = m) * (2/(df-a))


####
#PI per founder
temp3=NULL
for(i in 1:length(genotypes)){
	temp=dataP[dataP$gene==genotypes[i],]
	temp2 = CIindex(dataM = temp[,-c(1:2)])
	temp3=c(temp3,temp2)
}
PI_individual=temp3*(2/c(72,72,20,72,90,90,132,110))



#loop through founder genotypes 
PIindex=function(d, founders){
	PI = NULL
		for(i in 1:length(founders)){
			temp=CIindex(dataM=d[d$gene == founders[i], 3:ncol(d)])
			PI = c(PI, temp)	
		}
		
	return(PI)
}

PIraw = PIindex(d = dataP, founders = genotypes)

#calculate df to use for PI calculations
count=table(dataP$gene)
count2 = 2/sum(count*(count-1))
PI = sum(PIraw * count2)


#log ratio of PI:CI for each founder pair
#Note: mean of CI's calculated here will be different from the global estimate because founder mutations are treated differently. 
	n = length(genotypes)
	CIij = matrix(nrow=n, ncol=n)
	PIij = matrix(nrow=n, ncol=n)
	PIii = matrix(nrow=n, ncol=n) #matrix to collect individual PI genotypes in each comparison

	#i - loop through focal founder genotypes
	for(i in 1:(n-1)){
		ith = dataP[dataP$gene == genotypes[i], 1:ncol(dataP)]
		#j - loop through founders to pair with focal
		for(j in (i+1):n){
			jth = data=dataP[dataP$gene == genotypes[j], 1:ncol(dataP)]
			#make new data set with founder pair
			ij_temp = rbind(ith, jth)
			
			#remove columns from genotypes i & j matrix that correspond to mutations fixed in j (i.e., the later genotype)
			j.founder.mut = which(mut[,genotypes[j]]==1)
			j.founder.mut.genes = as.matrix(mut)[j.founder.mut,1]
			ij = ij_temp[,-which(colnames(ij_temp) %in% j.founder.mut.genes)]
						
			#CI for founder pair
			#Should be corrected for pre-exisintg mutations if remove.mutations=FALSE
			N = nrow(ij)
			CIij[i,j] = CIindex(ij[,3:ncol(ij)]) * (2/(N*(N-1)))
			#PI for founder pair
			PIraw = PIindex(ij, founders = c(genotypes[i], genotypes[j]))
			count=c(length(ij[ij$gene== genotypes[i],1]),length(ij[ij$gene== genotypes[j],1]))
			individual.PI = PIraw * (2/(count*(count-1)))
			PIij[i,j] = mean(individual.PI)
			#PI for i only -- except when i = 7, then j as well
			PIii[i,j] = individual.PI[1]
		}
	}
	
	PI.CI = log(PIij/CIij,2)
	CI.PI = log(CIij/PIij,2)
#add individual PI estimates -- with no gene omisssions - to Pii
for (i in 1:8){
	PIii[i,i] = PI_individual[i]
}
#Rows of PIii show how parallelism for a founder declines in comparisons with other founders as genes already in 
#those other founders are omitted frm the analysis. Later founders have more of those genes, so parallelism is the same or
#lower in a comparison with late vs. early founders. 


x = 1:length(genotypes)
y = 1:length(genotypes)
z = as.matrix(PI.CI)
z[8,8]=0 #provides baseline for colorscale

#z2 for random epistasis
z2=matrix(NA, ncol=8, nrow=8)
for(i in 1:7){
	for(j in (i+1):8){
		z2[i,j] = runif(1,min=0,max=0.2)
	}
}
#baseline for color scale
z2[8,8]=1 #provides baseline for colorscale


my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 149)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,0.4,length=50),  # for red
  seq(0.4,0.7,length=50),              # for yellow
  seq(0.7,1,length=50))              # for green

fig<-paste("Fig4_PI.CI.pdf",sep="")
fig.file<-paste(folder,fig,sep="")
pdf(fig.file,onefile=F,width=3.5,height=4)

 
par(fig=c(0,1,0,0.8)) 
par(mar=c(4,4,1,1))          
#PI:CI          
image(x,y,z,
	col=my_palette,#"white", #-use to plot PI only
	xlab="", ylab="",
	xaxt='n', yaxt='n')          

axis(1, at = 1:8, labels = FALSE)
text(1:8, par("usr")[3]-0.4 , labels = geno.lab, srt = 0, pos = 1, xpd = TRUE,cex=0.8,font=1)
mtext("Founder genotype", side = 1, line = 2)
axis(2, at = 1:8, labels = FALSE)
text(y=1:8, par("usr")[1]-.4 , labels = geno.lab, srt = 0, pos = 2, xpd = TRUE,cex=0.8,font=1)
mtext("Founder genotype", side = 2, line = 3)
        
par(new=TRUE)
#PI along diagonal
PIdiag=matrix(NA, nrow=8, ncol=8)
for(i in 1:8){
	PIdiag[i,i]=PI_individual[i]
}
PIdiag[7,5]=0 #provides baseline for blue color

image(x,y,PIdiag,
	col=c("white",(colorRampPalette(brewer.pal(9,"Blues"))(500))),
	xlab="", ylab="",
	xaxt='n', yaxt='n')       

#gridlines
grid(nx=8, col="white", lty=1, lwd=2)
box()

mtext(expression('log'[2]*' parallelism:convergence'),side=3, line=3.8, adj=1.05, cex=0.7)
mtext("Parallelism",side=3, line=1.1, adj=0.57, cex=0.7)


#color bars
par(fig=c(0,1,0.75,0.9), new=T) 
par(mar=c(1,4,0.7,8))

#PI color bar
image(c(1:50),1,(as.matrix(seq(0,1.5,length.out=50))), 
	col = c("white",(colorRampPalette(brewer.pal(9,"Blues"))(50))),
	xlab="",ylab="",xaxt='n',yaxt='n', useRaster=T)
axis(1, at = c(1,50), labels = FALSE, line=NA, tck=-0.15)
text(c(1,50),par("usr")[1]+0.2, labels=c(0,2), srt = 0, pos = 1, xpd = TRUE,cex=0.65,font=1)
	
par(fig=c(0,1,0.85,1), new=T) 
par(mar=c(1.5,4,0.2,8))
	
#PI:CI color bar
image(c(1:200),1,(as.matrix(seq(2,4,length.out=200))), col=my_palette,
	xlab="",ylab="",xaxt='n',yaxt='n', useRaster=T)
axis(1, at = c(5,195), labels = FALSE, line=NA, tck=-0.2)
text(c(5,195),par("usr")[1]+0.15, labels=c(0,1), srt = 0, pos = 1, xpd = TRUE,cex=0.65,font=1)

     
dev.off()

###############
#Mantel test to evaluate correlation between genetic distance and PI:CI
#mutation diff matrix
mut1 = mut[,2:9] #drop gene name row
mut1[is.na(mut1)] = 0#convert NA to 0
n = length(genotypes)
ma = matrix(nrow = n, ncol = n)
for(i in 1:(length(genotypes)-1)){
	for(j in (i+1):length(genotypes)){
		ith = mut1[,i] 
		jth = mut1[,j]
		#ij = length(na.omit(jth)-length(na.omit(ith)))
		#Hamming distance
		ij = sum(ith != jth)
		ma[i,j] = ij
	}
}

#fitness matrix
library(ecodist) #for Mantel test
fit = distance(fitness)

#make mutation matrix symmetric
ma2 = ma
ma2[is.na(ma2)] = 0
ma2 = ma2 + t(ma2)

#make symmetric PI:CI matrix
PI.CI2 = PI.CI
PI.CI2[is.na(PI.CI2)] = 0
PI.CI2 = PI.CI2 + (t(PI.CI2))

#ecodist mantel test
#Note: pval1 = 1-tailed null r>0; pval2 = 1-tailed r<0; pval3 = 2-tailed r=0 
mantelMA = mantel(lower(PI.CI2) ~ lower(ma2), mrank=F, nperm=100000) # " mrank = TRUE' for Spearman
mantelF = mantel(lower(PI.CI2) ~ lower(fit), mrank=F, nperm=10000)


