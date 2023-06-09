###phylogenetic simulation analyses

require(RPANDA)
require(mvMORPH)
require(MASS)
require(MCMCglmm)
require(overlapping)
nsim=1000

########################################
##### load in phylogeny and trait data
########################################

phylo<-read.tree("Hetaerina_phylogeny.tree")
trait<-read.csv("Dataset2A.csv")

#for each trait
#create data matrices with means and

##photoDFA1

photoDFA1_mean<-data.frame(trait$hw.basal.Red_mean,trait$hw.basal.Green_mean,trait$hw.basal.Blue_mean,trait$hw.tip.Red_mean,trait$hw.tip.Green_mean,trait$hw.tip.Blue_mean,trait$fw.basal.Red_mean,trait$fw.basal.Green_mean,trait$fw.basal.Blue_mean,trait$fw.tip.Red_mean,trait$fw.tip.Green_mean,trait$fw.tip.Blue_mean)
rownames(photoDFA1_mean)<-trait$treename
photoDFA1_mean<-photoDFA1_mean[phylo$tip.label,]

##photoDFA2

photoDFA2_mean<-data.frame(trait$basal.hw.weighted_mean,trait$distal.hw.weighted_mean,trait$basal.fw.weighted_mean,trait$distal.fw.weighted_mean)
rownames(photoDFA2_mean)<-trait$treename
photoDFA2_mean<-photoDFA2_mean[phylo$tip.label,]

##specDFA

specDFA_mean<-data.frame(trait$tot.int.basal.hw_mean,trait$tot.int.mid.hw_mean,trait$tot.int.tip.hw_mean,trait$tot.int.basal.fw_mean,trait$tot.int.mid.fw_mean,trait$tot.int.tip.fw_mean)
rownames(specDFA_mean)<-trait$treename
specDFA_mean<-specDFA_mean[phylo$tip.label,]
specDFA_mean<-specDFA_mean[complete.cases(specDFA_mean),]

##photographic lightness (both.weighted)

photo.lightness_mean<-trait$both.weighted_mean
names(photo.lightness_mean)<-trait$treename
photo.lightness_mean<-photo.lightness_mean[phylo$tip.label]

##weighted total intensity (spec measure of lightess)

spec.intensity_mean<-trait$int.weighted_mean[complete.cases(trait$int.weighted_mean)]
names(spec.intensity_mean)<-trait$treename[complete.cases(trait$int.weighted_mean)]
spec.intensity_mean<-spec.intensity_mean[drop.tip(phylo,c("H_capitalis","H_majuscula"))$tip.label]

########################################
##### fit evolutionary trait models
########################################

#fit multivariate trait models (for DFA analyses of different wing patches)


BM.photoDFA1<-fit_t_pl(Y=as.matrix(photoDFA1_mean),tree=phylo)
OU.photoDFA1<-fit_t_pl(Y=as.matrix(photoDFA1_mean),model="OU",tree=phylo)

BM.photoDFA2<-fit_t_pl(Y=as.matrix(photoDFA2_mean),tree=phylo)
OU.photoDFA2<-fit_t_pl(Y=as.matrix(photoDFA2_mean),tree=phylo,model="OU")

BM.specDFA<-fit_t_pl(tree=drop.tip(phylo,c("H_capitalis","H_majuscula")),Y=as.matrix(specDFA_mean))
OU.specDFA<-fit_t_pl(tree=drop.tip(phylo,c("H_capitalis","H_majuscula")),Y=as.matrix(specDFA_mean),model="OU")

BM.photoDFA1.sim <- mvSIM(BM.photoDFA1$corrstruct$phy, nsim=nsim, model="BM1", param=list(sigma=BM.photoDFA1$R$R,theta=ancestral(BM.photoDFA1)$root))
OU.photoDFA1.sim <- mvSIM(OU.photoDFA1$corrstruct$phy, nsim=nsim, model="BM1", param=list(sigma=OU.photoDFA1$R$R,theta=ancestral(OU.photoDFA1)$root))

BM.photoDFA2.sim <- mvSIM(BM.photoDFA2$corrstruct$phy, nsim=nsim, model="BM1", param=list(sigma=BM.photoDFA2$R$R,theta=ancestral(BM.photoDFA2)$root))
OU.photoDFA2.sim <- mvSIM(OU.photoDFA2$corrstruct$phy, nsim=nsim, model="BM1", param=list(sigma=OU.photoDFA2$R$R,theta=ancestral(OU.photoDFA2)$root))

BM.specDFA.sim <- mvSIM(BM.specDFA$corrstruct$phy, nsim=nsim, model="BM1", param=list(sigma=BM.specDFA$R$R,theta=ancestral(BM.specDFA)$root))
OU.specDFA.sim <- mvSIM(OU.specDFA$corrstruct$phy, nsim=nsim, model="BM1", param=list(sigma=OU.specDFA$R$R,theta=ancestral(OU.specDFA)$root))

#fit univariate trait models (for datasets using univariate indices of lightness)

BM.photoLightness<-mvBM(tree=phylo,data=photo.lightness_mean,model="BM1",method='pic')
OU.photoLightness<-mvOU(tree=phylo,data=photo.lightness_mean,model="OU1",method='inverse')

BM.specIntensity<-mvBM(tree=drop.tip(phylo,c("H_capitalis","H_majuscula")),data=spec.intensity_mean,model="BM1",method='pic')
OU.specIntensity<-mvOU(tree=drop.tip(phylo,c("H_capitalis","H_majuscula")),data=spec.intensity_mean,model="OU1",method='inverse')

BM.photoLightness.sim<-simulate(BM.photoLightness,nsim=nsim,tree=phylo)
OU.photoLightness.sim<-simulate(OU.photoLightness,nsim=nsim,tree=phylo)

BM.specIntensity.sim<-simulate(BM.specIntensity,nsim=nsim,tree=drop.tip(phylo,c("H_capitalis","H_majuscula")))
OU.specIntensity.sim<-simulate(OU.specIntensity,nsim=nsim,tree=drop.tip(phylo,c("H_capitalis","H_majuscula")))


########################################
##### simulate data for each site
########################################

#for DFA analyses, overall goal is to have a list for each DFA, each with N entries
#each list element is a data.frame with the DFA misrecognition probability and isolation index for each species X site combo

sitewise<-read.csv("Dataset2B.csv")
photos<-data.frame(sitewise$site,sitewise$treename,sitewise$n.photos)
photos<-subset(photos,sitewise.treename!="H_caja")
specs<-data.frame(sitewise$site,sitewise$treename,sitewise$n.scans)
specs<-subset(specs,sitewise.n.scans!="NA" & sitewise.treename!="H_caja")
photo.sites<-unique(photos$sitewise.site)
spec.sites<-unique(specs$sitewise.site)

all.ri<-read.csv("Dataset1.csv")
comps<-unique(all.ri$comparison)

BM.sim.photoDFA1.list<-list()
OU.sim.photoDFA1.list<-list()
BM.sim.photoDFA2.list<-list()
OU.sim.photoDFA2.list<-list()
BM.sim.specDFA.list<-list()
OU.sim.specDFA.list<-list()
BM.sim.photoDIFF.list<-list()
OU.sim.photoDIFF.list<-list()
BM.sim.photoOVERLAP.list<-list()
OU.sim.photoOVERLAP.list<-list()
BM.sim.specDIFF.list<-list()
OU.sim.specDIFF.list<-list()

#function takes a vector of mean values and a vector of sd values (of equal length)
#and returns a matrix simulated values from a truncated normal distribution using rtnorm in MCMCglmm

mrtnorm<-function(n.sim,mean.vec,sd.vec,lower=-Inf,upper=Inf){
	if(length(mean.vec)!=length(sd.vec)){stop("mean and standard deviation vectors of unequal length")}
	out=matrix(nrow=n.sim,ncol=length(mean.vec))
	for(i in 1:length(mean.vec)){
		out[,i]=rtnorm(n=n.sim,mean=mean.vec[i],sd=sd.vec[i],lower=lower,upper=upper)
	}
	return(out)
}

for(j in 1:nsim){
	
	BM.sim.photoDFA1.j=BM.photoDFA1.sim[[j]]
	OU.sim.photoDFA1.j=OU.photoDFA1.sim[[j]]
	BM.sim.photoDFA2.j=BM.photoDFA2.sim[[j]]
	OU.sim.photoDFA2.j=OU.photoDFA2.sim[[j]]
	BM.sim.specDFA.j=BM.specDFA.sim[[j]]
	OU.sim.specDFA.j=OU.specDFA.sim[[j]]
	BM.sim.photoLightness.j=BM.photoLightness.sim[,j]
	OU.sim.photoLightness.j=OU.photoLightness.sim[,j]
	BM.sim.specIntensity.j=BM.specIntensity.sim[,j]
	OU.sim.specIntensity.j=OU.specIntensity.sim[,j]

	hold1=matrix(ncol=5)
	hold2=matrix(ncol=5)
	hold3=matrix(ncol=5)
	hold4=matrix(ncol=5)
	hold5=matrix(ncol=5)
	hold6=matrix(ncol=5)
	hold7=matrix(ncol=5)
	hold8=matrix(ncol=5)
	hold9=matrix(ncol=5)
	hold10=matrix(ncol=5)
	hold11=matrix(ncol=5)
	hold12=matrix(ncol=5)

	colnames(hold1)<-c("site","treename","treename.HS","prop.hs.misidentified.as.focal.by.dfa.PHOTO1","isolation.index")
	colnames(hold2)<-c("site","treename","treename.HS","prop.hs.misidentified.as.focal.by.dfa.PHOTO1","isolation.index")
	colnames(hold3)<-c("site","treename","treename.HS","prop.hs.misidentified.as.focal.by.dfa.PHOTO2","isolation.index")
	colnames(hold4)<-c("site","treename","treename.HS","prop.hs.misidentified.as.focal.by.dfa.PHOTO2","isolation.index")
	colnames(hold5)<-c("site","treename","treename.HS","prop.hs.misidentified.as.focal.by.dfa.SPEC","isolation.index")
	colnames(hold6)<-c("site","treename","treename.HS","prop.hs.misidentified.as.focal.by.dfa.SPEC","isolation.index")
	colnames(hold7)<-c("site","treename","treename.HS","diff.photographic.lightness","isolation.index")
	colnames(hold8)<-c("site","treename","treename.HS","diff.photographic.lightness","isolation.index")
	colnames(hold9)<-c("site","treename","treename.HS","overlap.photo.lightness","isolation.index")
	colnames(hold10)<-c("site","treename","treename.HS","overlap.photo.lightness","isolation.index")
	colnames(hold11)<-c("site","treename","treename.HS","diff.weighted.spec.lightness","isolation.index")
	colnames(hold12)<-c("site","treename","treename.HS","diff.weighted.spec.lightness","isolation.index")

	#create permuted dataset
	for(p in 1:length(comps)){
		cn<-comps[p]
		if(p==1){
		perm.mat<-all.ri[which(all.ri$comparison==cn)[sample(1:dim(all.ri[which(all.ri$comparison==cn),])[1],1)],]
		} else {
		hold<-all.ri[which(all.ri$comparison==cn)[sample(1:dim(all.ri[which(all.ri$comparison==cn),])[1],1)],]
		perm.mat<-rbind(perm.mat,hold)
		}
	}

	for(i in 1:length(photo.sites)){
		site=subset(photos,sitewise.site==photo.sites[i])
		site.data<-matrix(ncol=35)
		site.listBM<-list()
		site.listOU<-list()
		for(n in 1:dim(site)[1]){
			n.samp=site[n,]$sitewise.n.photos	
							
			sp=as.character(rep(site[n,2],n.samp))
			
			simdat=cbind(sp,mrtnorm(n.sim=n.samp,mean.vec=BM.sim.photoDFA1.j[as.character(site[n,2]),],sd.vec=trait[which(trait$treename==as.character(site[n,2])),c(15,17,19,21,23,25,27,29,31,33,35,37)],lower=0),mrtnorm(n.sim=n.samp,mean.vec=OU.sim.photoDFA1.j[as.character(site[n,2]),],sd.vec=trait[which(trait$treename==as.character(site[n,2])),c(15,17,19,21,23,25,27,29,31,33,35,37)],lower=0),mrtnorm(n.sim=n.samp,mean.vec=BM.sim.photoDFA2.j[as.character(site[n,2]),],sd.vec=trait[which(trait$treename==as.character(site[n,2])),c(5,7,9,11)],lower=0),mrtnorm(n.sim=n.samp,mean.vec=OU.sim.photoDFA2.j[as.character(site[n,2]),],sd.vec=trait[which(trait$treename==as.character(site[n,2])),c(5,7,9,11)],lower=0),rtnorm(n=n.samp,mean=BM.sim.photoLightness.j[as.character(site[n,2])],sd=trait[which(trait$treename==as.character(site[n,2])),13],lower=0),rtnorm(n=n.samp,mean=OU.sim.photoLightness.j[as.character(site[n,2])],sd=trait[which(trait$treename==as.character(site[n,2])),13],lower=0))
			site.data=rbind(site.data,simdat)
			site.listBM[[n]]<-as.numeric(as.character(simdat[,34]))
			site.listOU[[n]]<-as.numeric(as.character(simdat[,35]))
			#this loop should produce a matrix, where each row is an individual
			#and columns are sp.ID and simulated trait values that will be passed to DFA analyses
		}
		site.data<-as.data.frame(site.data[-1,])
		site.data[,2:35]<-sapply(site.data[,2:35],function(x)as.numeric(as.character(x)))	

		BM.dfa1<-lda(site.data[,1]~site.data[,2]+site.data[,3]+site.data[,4]+site.data[,5]+site.data[,6]+site.data[,7]+site.data[,8]+site.data[,9]+site.data[,10]+site.data[,11]+site.data[,12]+site.data[,13])
		OU.dfa1<-lda(site.data[,1]~site.data[,14]+site.data[,15]+site.data[,16]+site.data[,17]+site.data[,18]+site.data[,19]+site.data[,20]+site.data[,21]+site.data[,22]+site.data[,23]+site.data[,24]+site.data[,25])
		BM.dfa2<-lda(site.data[,1]~site.data[,26]+site.data[,27]+site.data[,28]+site.data[,29])
		OU.dfa2<-lda(site.data[,1]~site.data[,30]+site.data[,31]+site.data[,32]+site.data[,33])

		names(site.listBM)<-site$sitewise.treename
		names(site.listOU)<-site$sitewise.treename
		BM.ov<-overlap(site.listBM)$OV
		OU.ov<-overlap(site.listOU)$OV

		freqtable.BM.dfa1<-table(site.data[,1],predict(BM.dfa1)$class)
		freqtable.OU.dfa1<-table(site.data[,1],predict(OU.dfa1)$class)
		freqtable.BM.dfa2<-table(site.data[,1],predict(BM.dfa2)$class)
		freqtable.OU.dfa2<-table(site.data[,1],predict(OU.dfa2)$class)

		#calculate misrecognition probability
		target.comps=subset(all.ri,site==photo.sites[i])
		for(p in 1:dim(target.comps)[1]){
			int1=c(as.character(site[1,1]),as.character(target.comps[p,3]),as.character(target.comps[p,5]),freqtable.BM.dfa1[as.character(target.comps[p,5]),as.character(target.comps[p,3])]/rowSums(freqtable.BM.dfa1)[as.character(target.comps[p,5])],target.comps[p,22])
			hold1=rbind(hold1,int1)
			int2=c(as.character(site[1,1]),as.character(target.comps[p,3]),as.character(target.comps[p,5]),freqtable.OU.dfa1[as.character(target.comps[p,5]),as.character(target.comps[p,3])]/rowSums(freqtable.OU.dfa1)[as.character(target.comps[p,5])],target.comps[p,22])
			hold2=rbind(hold2,int2)
			int3=c(as.character(site[1,1]),as.character(target.comps[p,3]),as.character(target.comps[p,5]),freqtable.BM.dfa2[as.character(target.comps[p,5]),as.character(target.comps[p,3])]/rowSums(freqtable.BM.dfa2)[as.character(target.comps[p,5])],target.comps[p,22])
			hold3=rbind(hold3,int3)
			int4=c(as.character(site[1,1]),as.character(target.comps[p,3]),as.character(target.comps[p,5]),freqtable.OU.dfa2[as.character(target.comps[p,5]),as.character(target.comps[p,3])]/rowSums(freqtable.OU.dfa2)[as.character(target.comps[p,5])],target.comps[p,22])
			hold4=rbind(hold4,int4)
		}
		
		target.comps2=subset(perm.mat,site==photo.sites[i])
		for(p in 1:dim(target.comps2)[1]){
			int7=c(as.character(site[1,1]),as.character(target.comps2[p,3]),as.character(target.comps2[p,5]),abs(BM.sim.photoLightness.j[as.character(target.comps2[p,5])]-BM.sim.photoLightness.j[as.character(target.comps2[p,3])]),target.comps2[p,22])
			hold7=rbind(hold7,int7)
			int8=c(as.character(site[1,1]),as.character(target.comps2[p,3]),as.character(target.comps2[p,5]),abs(OU.sim.photoLightness.j[as.character(target.comps2[p,5])]-OU.sim.photoLightness.j[as.character(target.comps2[p,3])]),target.comps2[p,22])
			hold8=rbind(hold8,int8)
			int9=c(as.character(site[1,1]),as.character(target.comps2[p,3]),as.character(target.comps2[p,5]),BM.ov[intersect(grep(as.character(target.comps2[p,5]),names(BM.ov)),grep(as.character(target.comps2[p,3]),names(BM.ov)))],target.comps2[p,22])
			hold9=rbind(hold9,int9)
			int10=c(as.character(site[1,1]),as.character(target.comps2[p,3]),as.character(target.comps2[p,5]),OU.ov[intersect(grep(as.character(target.comps2[p,5]),names(OU.ov)),grep(as.character(target.comps2[p,3]),names(OU.ov)))],target.comps2[p,22])
			hold10=rbind(hold10,int10)
		}
		
	}

	
	BM.sim.photoDFA1.list[[j]]<- as.data.frame(hold1[-1,],row.names=FALSE)
	OU.sim.photoDFA1.list[[j]]<- as.data.frame(hold2[-1,],row.names=FALSE)
	BM.sim.photoDFA2.list[[j]]<- as.data.frame(hold3[-1,],row.names=FALSE)
	OU.sim.photoDFA2.list[[j]]<- as.data.frame(hold4[-1,],row.names=FALSE)
	BM.sim.photoDIFF.list[[j]]<- as.data.frame(hold7[-1,],row.names=FALSE)
	OU.sim.photoDIFF.list[[j]]<- as.data.frame(hold8[-1,],row.names=FALSE)
	BM.sim.photoOVERLAP.list[[j]]<- as.data.frame(hold9[-1,],row.names=FALSE)
	OU.sim.photoOVERLAP.list[[j]]<- as.data.frame(hold10[-1,],row.names=FALSE)

	for(i in 1:length(spec.sites)){
		site=subset(specs,sitewise.site==spec.sites[i])
		site.data<-matrix(ncol=13)
		for(n in 1:dim(site)[1]){
			n.samp=site[n,]$sitewise.n.scans
			sp=as.character(rep(site[n,2],n.samp))
			#simdat=cbind(sp,mvrnorm(n=n.samp,mu=BM.sim.specDFA.j[as.character(site[n,2]),],Sigma=Sigma.spec),mvrnorm(n=n.samp,mu=OU.sim.specDFA.j[as.character(site[n,2]),],Sigma=Sigma.spec))
			simdat=cbind(sp,mrtnorm(n.sim=n.samp,mean.vec=BM.sim.specDFA.j[as.character(site[n,2]),],sd.vec=trait[which(trait$treename==as.character(site[n,2])),c(40,42,44,46,48,50)],lower=0),mrtnorm(n.sim=n.samp,mean.vec=OU.sim.specDFA.j[as.character(site[n,2]),],sd.vec=trait[which(trait$treename==as.character(site[n,2])),c(40,42,44,46,48,50)],lower=0))			
			site.data=rbind(site.data,simdat)
		}
		site.data<-as.data.frame(site.data[-1,])
		site.data[,2:13]<-sapply(site.data[,2:13],function(x)as.numeric(as.character(x)))	
		BM.dfas<-lda(site.data[,1]~site.data[,2]+site.data[,3]+site.data[,4]+site.data[,5]+site.data[,6]+site.data[,7])
		OU.dfas<-lda(site.data[,1]~site.data[,8]+site.data[,9]+site.data[,10]+site.data[,11]+site.data[,12]+site.data[,13])
		freqtable.BM.dfas<-table(site.data[,1],predict(BM.dfas)$class)
		freqtable.OU.dfas<-table(site.data[,1],predict(OU.dfas)$class)
		
		target.comps=subset(all.ri,site==spec.sites[i] & prop.hs.misidentified.as.focal.by.dfa.SPEC!="NA")
		for(p in 1:dim(target.comps)[1]){
			int5=c(as.character(site[1,1]),as.character(target.comps[p,3]),as.character(target.comps[p,5]),freqtable.BM.dfas[as.character(target.comps[p,5]),as.character(target.comps[p,3])]/rowSums(freqtable.BM.dfas)[as.character(target.comps[p,5])],target.comps[p,22])
			hold5=rbind(hold5,int5)
			int6=c(as.character(site[1,1]),as.character(target.comps[p,3]),as.character(target.comps[p,5]),freqtable.OU.dfas[as.character(target.comps[p,5]),as.character(target.comps[p,3])]/rowSums(freqtable.OU.dfas)[as.character(target.comps[p,5])],target.comps[p,22])
			hold6=rbind(hold6,int6)
		}
		
		target.comps2=subset(perm.mat,site==spec.sites[i])
		for(p in 1:dim(target.comps2)[1]){
			int11=c(as.character(site[1,1]),as.character(target.comps2[p,3]),as.character(target.comps2[p,5]),abs(BM.sim.specIntensity.j[as.character(target.comps2[p,5])]-BM.sim.specIntensity.j[as.character(target.comps2[p,3])]),target.comps2[p,22])
			hold11=rbind(hold11,int11)
			int12=c(as.character(site[1,1]),as.character(target.comps2[p,3]),as.character(target.comps2[p,5]),abs(OU.sim.specIntensity.j[as.character(target.comps2[p,5])]-OU.sim.specIntensity.j[as.character(target.comps2[p,3])]),target.comps2[p,22])
			hold12=rbind(hold12,int12)
		}
	}
		
	BM.sim.specDFA.list[[j]]<-as.data.frame(hold5[-1,],row.names=FALSE)
	OU.sim.specDFA.list[[j]]<-as.data.frame(hold6[-1,],row.names=FALSE)
	BM.sim.specDIFF.list[[j]]<-as.data.frame(hold11[-1,],row.names=FALSE)
	OU.sim.specDIFF.list[[j]]<-as.data.frame(hold12[-1,],row.names=FALSE)

	if(j%%20==0){print(j)}
}
	

########################################
##### run spearman correlations on each list element
########################################

res.mat<-matrix(ncol=25,nrow=nsim)
colnames(res.mat)=c("no.","BM.photoDFA1.rho","BM.photoDFA1.p","OU.photoDFA1.rho","OU.photoDFA1.p","BM.photoDFA2.rho","BM.photoDFA2.p","OU.photoDFA2.rho","OU.photoDFA2.p","BM.specDFA.rho","BM.specDFA.p","OU.specDFA.rho","OU.specDFA.p","BM.photoDIFF.rho","BM.photoDIFF.p","OU.photoDIFF.rho","OU.photoDIFF.p","BM.photoOVERLAP.rho","BM.photoOVERLAP.p","OU.photoOVERLAP.rho","OU.photoOVERLAP.p","BM.specDIFF.rho","BM.specDIFF.p","OU.specDIFF.rho","OU.specDIFF.p")

#load in each of the simulated datasets
#run correlation

for(i in 1:nsim){
	
	BM.ph1.cor=cor.test(x=as.numeric(as.character(BM.sim.photoDFA1.list[[i]]$prop.hs.misidentified.as.focal.by.dfa.PHOTO1)),y=as.numeric(as.character(BM.sim.photoDFA1.list[[i]]$isolation.index)),method="spearman")
	OU.ph1.cor=cor.test(x=as.numeric(as.character(OU.sim.photoDFA1.list[[i]]$prop.hs.misidentified.as.focal.by.dfa.PHOTO1)),y=as.numeric(as.character(OU.sim.photoDFA1.list[[i]]$isolation.index)),method="spearman")
	BM.ph2.cor=cor.test(x=as.numeric(as.character(BM.sim.photoDFA2.list[[i]]$prop.hs.misidentified.as.focal.by.dfa.PHOTO2)),y=as.numeric(as.character(BM.sim.photoDFA2.list[[i]]$isolation.index)),method="spearman")
	OU.ph2.cor=cor.test(x=as.numeric(as.character(OU.sim.photoDFA2.list[[i]]$prop.hs.misidentified.as.focal.by.dfa.PHOTO2)),y=as.numeric(as.character(OU.sim.photoDFA2.list[[i]]$isolation.index)),method="spearman")
	BM.spe.cor=cor.test(x=as.numeric(as.character(BM.sim.specDFA.list[[i]]$prop.hs.misidentified.as.focal.by.dfa.SPEC)),y=as.numeric(as.character(BM.sim.specDFA.list[[i]]$isolation.index)),method="spearman")
	OU.spe.cor=cor.test(x=as.numeric(as.character(OU.sim.specDFA.list[[i]]$prop.hs.misidentified.as.focal.by.dfa.SPEC)),y=as.numeric(as.character(OU.sim.specDFA.list[[i]]$isolation.index)),method="spearman")
	BM.phd.cor=cor.test(x=as.numeric(as.character(BM.sim.photoDIFF.list[[i]]$diff.photographic.lightness)),y=as.numeric(as.character(BM.sim.photoDIFF.list[[i]]$isolation.index)),method="spearman")
	OU.phd.cor=cor.test(x=as.numeric(as.character(OU.sim.photoDIFF.list[[i]]$diff.photographic.lightness)),y=as.numeric(as.character(OU.sim.photoDIFF.list[[i]]$isolation.index)),method="spearman")
	BM.pho.cor=cor.test(x=as.numeric(as.character(BM.sim.photoOVERLAP.list[[i]]$overlap.photo.lightness)),y=as.numeric(as.character(BM.sim.photoOVERLAP.list[[i]]$isolation.index)),method="spearman")
	OU.pho.cor=cor.test(x=as.numeric(as.character(OU.sim.photoOVERLAP.list[[i]]$overlap.photo.lightness)),y=as.numeric(as.character(OU.sim.photoOVERLAP.list[[i]]$isolation.index)),method="spearman")
	BM.spd.cor=cor.test(x=as.numeric(as.character(BM.sim.specDIFF.list[[i]]$diff.weighted.spec.lightness)),y=as.numeric(as.character(BM.sim.specDIFF.list[[i]]$isolation.index)),method="spearman")
	OU.spd.cor=cor.test(x=as.numeric(as.character(OU.sim.specDIFF.list[[i]]$diff.weighted.spec.lightness)),y=as.numeric(as.character(OU.sim.specDIFF.list[[i]]$isolation.index)),method="spearman")

	res.mat[i,]<-c(i,BM.ph1.cor$estimate,BM.ph1.cor$p.value,OU.ph1.cor$estimate,OU.ph1.cor$p.value,BM.ph2.cor$estimate,BM.ph2.cor$p.value,OU.ph2.cor$estimate,OU.ph2.cor$p.value,BM.spe.cor$estimate,BM.spe.cor$p.value,OU.spe.cor$estimate,OU.spe.cor$p.value,BM.phd.cor$estimate,BM.phd.cor$p.value,OU.phd.cor$estimate,OU.phd.cor$p.value,BM.pho.cor$estimate,BM.pho.cor$p.value,OU.pho.cor$estimate,OU.pho.cor$p.value,BM.spd.cor$estimate,BM.spd.cor$p.value,OU.spd.cor$estimate,OU.spd.cor$p.value)
	if(i%%100==0){print(i)}
}

resm<-as.data.frame(res.mat)


########################################
##### run one-sample t-tests against phylo simulations
########################################

#mu values are mean empirical rho values from Script1.R
t.test(resm$BM.photoDFA1.rho,mu=-0.52)
t.test(resm$OU.photoDFA1.rho,mu=-0.52)
t.test(resm$BM.photoDFA2.rho,mu=-0.59)
t.test(resm$OU.photoDFA2.rho,mu=-0.59)
t.test(resm$BM.specDFA.rho,mu=-0.75)
t.test(resm$OU.specDFA.rho,mu=-0.75)
t.test(resm$BM.photoDIFF.rho,mu=0.54)
t.test(resm$OU.photoDIFF.rho,mu=0.54)
t.test(resm$BM.photoOVERLAP.rho,mu=-0.64)
t.test(resm$OU.photoOVERLAP.rho,mu=-0.64)
t.test(resm$BM.specDIFF.rho,mu=0.72)
t.test(resm$OU.specDIFF.rho,mu=0.72)
