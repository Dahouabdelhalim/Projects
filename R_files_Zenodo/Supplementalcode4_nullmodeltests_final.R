#### Supplemental Code for Madsen, Vander Meiden & Shizuka: Social partners and temperature jointly affect morning foraging activity of small birds in winter#####

### PART 3: Permutations of networks#####

##required libraries
library(asnipe)
library(igraph)
library(iterators) #for parallel processing
library(foreach) #for parallel processing
library(parallel) #for parallel processing
library(doParallel) #for parallel processing


#import files from Gaussian Mixture Models
gmmDOWO=readRDS("conspecificDOWOflocks_final.rds") #import gmm results file. Will need this for permuting group-by-individual matrices.
gmmWBNU=readRDS("conspecificWBNUflocks_final.rds") #import gmm results file. Will need this for permuting group-by-individual matrices.

#extract group-by-individual matrices
gbi_dowo=gmmDOWO$gbi
gbi_wbnu=gmmWBNU$gbi

#make adjacency matrices
dowoadj=get_network(gbi_dowo)
wbnuadj=get_network(gbi_wbnu)

diag(dowoadj)=NA #make diagonal of adjacency matrices NA so we don't count these when normalizing values later.
diag(wbnuadj)=NA #make diagonal of adjacency matrices NA so we don't count these when normalizing values later.

###import activity correlation matrix (simmat) and spatial overlap matrix (logmat) from Supplemental Code 3. Note that this matrix has both species together, and they will be separated in the next set of codes

simmat=as.matrix(read.csv("simmat.csv", row.names = 1, check.names = F)) #import correlation matrix of daily activity (Z-scores)
logmat=as.matrix(read.csv("logmat.csv", row.names = 1, check.names = F)) #import correlation matrix of proportional feeder use

#extract & sort the activity correlation matrix & spatial overlap matrix for DOWO so rows/columns match adjacency matrix
dowosim=simmat[match(rownames(dowoadj), rownames(simmat)), match(rownames(dowoadj), rownames(simmat))] 
dowospat=logmat[match(rownames(dowoadj), rownames(logmat)), match(rownames(dowoadj), rownames(logmat))] 


#same for WBNU
wbnusim=simmat[match(rownames(wbnuadj), rownames(simmat)), match(rownames(wbnuadj), rownames(simmat))]
wbnuspat=logmat[match(rownames(wbnuadj), rownames(logmat)), match(rownames(wbnuadj), rownames(logmat))]

# 
#function to normalize values
normalize_matrix=function(m){
  (m-min(m, na.rm=T))/(max(m, na.rm=T)-min(m, na.rm=T))
}

dowosim.norm=normalize_matrix(dowosim)
dowoadj.norm=normalize_matrix(dowoadj)
dowospat.norm=normalize_matrix(dowospat)

###DOWO permutations

#setup
dowo.ids=rownames(dowoadj) #get DOWO ids
dowogbi=gbi_dowo[,which(colnames(gbi_dowo)%in%dowo.ids)] #get gbi with only DOWOs
dowogbi.filt=dowogbi[which(rowSums(dowogbi)>0),] #remove groups that no DOWOs belong to.
dowometadata.filt=gmmDOWO$metadata[which(rowSums(dowogbi)>0),] #remove the same groups in the group metadata

#using the filtered metadata, we can extract the feeder & date of the group. This will be useful when we constrain the permutation by day
dowo.locations=as.numeric(as.factor(substr(dowometadata.filt$Location, start=1, stop=8)))
dowo.days=as.numeric(as.factor(substr(dowometadata.filt$Location, start=10, stop=13)))


#store the results of MRQAP with empirical network
emp.mod.dowo=mrqap.dsp(dowosim.norm~dowoadj.norm+dowospat.norm)
emp.coef.dowo=emp.mod.dowo$coefficients[2]



### DOWO permutations with parallel processing
times=1000 #run 10 times to check. In full analysis, change to 1000
n.cores=detectCores()
system.time({
  registerDoParallel(n.cores)
  dowo.results.parallel=foreach(i = 1:times) %dopar% as.matrix(network_swap(dowogbi.filt, data_format="GBI", swaps=10000, association_matrix = dowoadj, locations=dowo.locations, days=dowo.days, within_day=TRUE, within_location=FALSE)$Association_index)
})
stopImplicitCluster()

dowoperm.adjs=dowo.results.parallel

dowoperm.mrqap.coefs=vector(length=times)
for(j in 1:times){
  adj=dowoperm.adjs[[j]]
  diag(adj)=NA
  adj.norm=normalize_matrix(adj)
  dowoperm.mrqap.coefs[j]=mrqap.dsp(dowosim.norm~adj.norm+dowospat.norm, randomisations=1)$coefficients[2]
}

dowoperm.mrqap.coefs

p.dowo=(length(which(dowoperm.mrqap.coefs>=emp.coef.dowo))+1)/(times+1)
p.dowo #p-value

ci.dowo=quantile(dowoperm.mrqap.coefs, probs=c(0.025, 0.975))
ci.dowo #confidence interval of MRQAP coefficients from null model 

#run this when you want to save the results
save(dowoperm.adjs, dowoperm.mrqap.coefs, p.dowo, ci.dowo, file="dowo_results_201005.rdata") 

### WBNU with parallel processing!
wbnusim.norm=normalize_matrix(wbnusim)
wbnuadj.norm=normalize_matrix(wbnuadj)
wbnuspat.norm=normalize_matrix(wbnuspat)

###wbnu permutations

#setup
wbnu.ids=rownames(wbnuadj) #get wbnu ids
wbnugbi=gbi_wbnu[,which(colnames(gbi_wbnu)%in%wbnu.ids)] #get gbi with only wbnus
wbnugbi.filt=wbnugbi[which(rowSums(wbnugbi)>0),] #remove groups that no wbnus belong to.
wbnumetadata.filt=gmmWBNU$metadata[which(rowSums(wbnugbi)>0),] #remove the same groups in the group metadata

#using the filtered metadata, we can extract the feeder & date of the group. This will be useful when we constrain the permutation by day
wbnu.locations=as.numeric(as.factor(substr(wbnumetadata.filt$Location, start=1, stop=8)))
wbnu.days=as.numeric(as.factor(substr(wbnumetadata.filt$Location, start=10, stop=13)))


#store the results of MRQAP with empirical network
emp.mod.wbnu=mrqap.dsp(wbnusim.norm~wbnuadj.norm+wbnuspat.norm)
emp.coef.wbnu=emp.mod.wbnu$coefficients[2]


library(foreach)
library(parallel)
library(doParallel)
times=1000
n.cores=detectCores()
system.time({
  registerDoParallel(n.cores)
wbnu.results.parallel=foreach(i = 1:times) %dopar% as.matrix(network_swap(wbnugbi.filt, data_format="GBI", swaps=10000, association_matrix = wbnuadj, locations=wbnu.locations, days=wbnu.days, within_day=TRUE, within_location=FALSE)$Association_index)
})
stopImplicitCluster()

wbnuperm.adjs=wbnu.results.parallel

wbnuperm.mrqap.coefs=vector(length=times)
for(j in 1:times){
  adj=wbnuperm.adjs[[j]]
  diag(adj)=NA
  adj.norm=normalize_matrix(adj)
  wbnuperm.mrqap.coefs[j]=mrqap.dsp(wbnusim.norm~adj.norm+wbnuspat.norm, randomisations=1)$coefficients[2]
}

wbnuperm.mrqap.coefs

p.wbnu=(length(which(wbnuperm.mrqap.coefs>=emp.coef.wbnu))+1)/(times+1)
p.wbnu #p-value

ci.wbnu=quantile(wbnuperm.mrqap.coefs, probs=c(0.025, 0.975))
ci.wbnu #confidence interval of MRQAP coefficients from null model 

#run this when you want to save the results
save(wbnuperm.adjs, wbnuperm.mrqap.coefs, p.wbnu, ci.wbnu, file="wbnu_results_201005.rdata")
