#This script conducts comparative analyses of reproductive isolation from male mate recognition ~ female wing coloration
Â£
#For phylogenetic analyses, see Script2.R

all.ri<-read.csv("Dataset1.csv")

################################
######## with photo data, raw data analyses ########
################################

c1<-cor.test(x=all.ri$prop.hs.misidentified.as.focal.by.dfa.PHOTO1,y=all.ri$isolation.index,method="spearman")

c2<-cor.test(x=all.ri$prop.hs.misidentified.as.focal.by.dfa.PHOTO2,y=all.ri$isolation.index,method="spearman")




################################
######## with spec data, raw data analyses ########
################################
all.ri.spec<-all.ri[which(!is.na(all.ri$diff.weighted.spec)),]
all.ri.spec.dfa<-all.ri[which(!is.na(all.ri$prop.hs.misidentified.as.focal.by.dfa.SPEC)),]

c3<-cor.test(x=all.ri.spec.dfa$prop.hs.misidentified.as.focal.by.dfa.SPEC,y=all.ri.spec.dfa$isolation.index,method="spearman")

##############################
###### permutation approach for remaining variables (to deal with non-independence arising from the similarity indices being identical for both species in a comparison)
##############################

#Permutation approach (see below), where each run is 1e4 permutations, col 1 is spearman's rho, and col 2 is p value

######
##difference in photographic lightness index
######

res.mat<-matrix(nrow=1E4,ncol=2)
colnames(res.mat)<-c("rho","p")

f.data<-all.ri
comps<-unique(all.ri$comparison)

for(i in 1:1E4){
	
	#create permuted dataset
	for(j in 1:length(comps)){
		cn<-comps[j]
		if(j==1){
		perm.mat<-f.data[which(f.data$comparison==cn)[sample(1:dim(f.data[which(f.data$comparison==cn),])[1],1)],]
		} else {
		hold<-f.data[which(f.data$comparison==cn)[sample(1:dim(f.data[which(f.data$comparison==cn),])[1],1)],]
		perm.mat<-rbind(perm.mat,hold)
		}
	}

	#run spearman correlation
	
	c1<-cor.test(x=perm.mat$diff.photographic.lightness,y=perm.mat$isolation.index,method="spearman")
	
	#store summary statistics
	res.mat[i,]<-c(c1$estimate,c1$p.value)
	
}
mean(res.mat[,1])
sd(res.mat[,1])
quantile(res.mat[,2],c(0.025,0.975))

######
##overlap index
#######

res.mat<-matrix(nrow=1E4,ncol=2)
colnames(res.mat)<-c("rho","p")

f.data<-all.ri
comps<-unique(all.ri$comparison)

for(i in 1:1E4){
	
	#create permuted dataset
	for(j in 1:length(comps)){
		cn<-comps[j]
		if(j==1){
		perm.mat<-f.data[which(f.data$comparison==cn)[sample(1:dim(f.data[which(f.data$comparison==cn),])[1],1)],]
		} else {
		hold<-f.data[which(f.data$comparison==cn)[sample(1:dim(f.data[which(f.data$comparison==cn),])[1],1)],]
		perm.mat<-rbind(perm.mat,hold)
		}
	}

	#run spearman correlation
	
	c1<-cor.test(x=perm.mat$overlap.photo.lightness,y=perm.mat$isolation.index,method="spearman")
	
	#store summary statistics
	res.mat[i,]<-c(c1$estimate,c1$p.value)
	
}
mean(res.mat[,1])
sd(res.mat[,1])
quantile(res.mat[,2],c(0.025,0.975))



######
## difference in spectral intensity
########

res.mat<-matrix(nrow=1E4,ncol=2)
colnames(res.mat)<-c("rho","p")

f.data<-all.ri.spec
comps<-unique(all.ri.spec$comparison)

for(i in 1:1E4){
	
	#create permuted dataset
	for(j in 1:length(comps)){
		cn<-comps[j]
		if(j==1){
		perm.mat.spec<-f.data[which(f.data$comparison==cn)[sample(1:dim(f.data[which(f.data$comparison==cn),])[1],1)],]
		} else {
		hold<-f.data[which(f.data$comparison==cn)[sample(1:dim(f.data[which(f.data$comparison==cn),])[1],1)],]
		perm.mat.spec<-rbind(perm.mat.spec,hold)
		}
	}

	#run spearman correlation
	
	c1<-cor.test(x=perm.mat.spec$diff.weighted.spec.lightness,y=perm.mat.spec$isolation.index,method="spearman")
	
	#store summary statistics
	res.mat[i,]<-c(c1$estimate,c1$p.value)
	
}
mean(res.mat[,1])
sd(res.mat[,1])
quantile(res.mat[,2],c(0.025,0.975))



