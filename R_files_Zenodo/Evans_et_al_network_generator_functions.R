#SYMMAT UTILITY FUNC
symmat = function (bool,dyads){
	return(rbind(dyads[bool,],dyads[bool,c(2,1)]))
}

#MAIN WRAPPER FUNCTION
do_networksim = function(groups,mean.group.size,max.group.size,d.eff,i.dens,o.dens,m.i.eff,m.o.eff,sex.eff,
	obs.eff.v,timesteps,floaterprob,probnorm,
	nreps,exportdir){
	#generate a network and sample at a variety of observation efforts
	#returns both observation based network and grouping event based network
	#groups,mean.group.size,max.group.size,d.eff,i.dens,o.dens,m.i.eff,m.o.eff,sex.eff - see network.generator
	#nreps = number of networks to simulate with these parameters
	#obs.eff.v = a vector of observation efforts
	#timesteps, floaterprob, probnorm - see networkobs
	
	#result folder names
	resultfolders=c("popdat","obsnet","obsgbimat","obsgbigroups","interactions")
	
	#if export directory does not exist, create it
	if(!dir.exists(as.character(exportdir))){
		dir.create(as.character(exportdir))
	}
	#always export the "true" network
	obs.eff.v=c(1,obs.eff.v)

	for (rep in 1:nreps){
		#generate a network
		simulated.networks=network.generator(groups,mean.group.size,max.group.size,d.eff,o.dens,i.dens,sex.eff,m.i.eff,m.o.eff)
		#for each observation effort
		for(obs.eff.c in obs.eff.v){
			parvec=paste(d.eff,i.dens,o.dens,m.i.eff,sex.eff,obs.eff.c,timesteps,sep="_")#parameters of interest
			#observe the network
			obs.sim.networks=networkobs(simulated.networks,timesteps, obs.eff.c, floaterprob,probnorm)
			
			#if parameter set folder does not exist, create
			if(!dir.exists(file.path(exportdir,parvec))){
				dir.create(file.path(exportdir,parvec))
			}


			#if export subfolders do not exist, create
			for(folder in resultfolders){
				fulldir=file.path(exportdir,parvec,folder,sep="/")
				if(!dir.exists(fulldir)){
					dir.create(fulldir)
				}
			}

			#export - pop info, observed networks, observed gbimats,observed gbigroups, 
			write.csv(simulated.networks$ind_data,paste(exportdir,"/",parvec,"/popdat/",rep,".csv",sep=""),row.names=F)
			write.csv(obs.sim.networks$obsnetwork,paste(exportdir,"/",parvec,"/obsnet/",rep,".csv",sep=""),row.names=F)
			write.csv(obs.sim.networks$obsgbigroups,paste(exportdir,"/",parvec,"/obsgbigroups/",rep,".csv",sep=""),row.names=F)
			write.csv(obs.sim.networks$obsgbi,paste(exportdir,"/",parvec,"/obsgbimat/",rep,".csv",sep=""),row.names=F)
		}
	}

}

###NETWORK GENERATION FUNCTION####
network.generator<-function(groups,mean.group.size,max.group.size,d.eff,o.dens,i.dens,sex.eff=NA,m.i.eff=NA,m.o.eff=NA){
	#####
	#groups = number of groups
	#mean.group.size = average groups
	#max.group.size = maximum group size
	#Network generator: a function to generate simulated animal networks
	#i.dens = density of within group associations
	#o.dens = density of outside group associations
	#d.eff = effect of distance between group centroids on association strength
	#m.i.eff = effect of being a male on within group assoc
	#m.o.eff = effect of being a male on outside group assoc
	#sex.eff = effect of individuals being same sex on their assoc strength
	#
	#####

	require(igraph)

	if(is.na(m.i.eff)){
		m.i.eff=0
	}

	if(is.na(m.o.eff)){
		m.o.eff=m.i.eff # if no outside sex effect included, same as inside
	}

	if(is.na(sex.eff)){
		sex.eff=1
	}

	#Scale m effects by i.eff and o.eff

	m.i.eff=m.i.eff*i.dens
	m.o.eff=m.o.eff*o.dens

	###SETUP BASE POPULATION####
	#and use this to calculate the overall size of the population
	n.indivs<-mean.group.size*groups

	#create the individuals
	indivs<-seq(1,n.indivs,1)

##and sample individuals into groups
	#size of possible groups is equal to max group size, then sample the overall population size from this
	poss.groups<-rep(seq(1,groups,1),each=max.group.size)
	indiv.groups<-sample(poss.groups,n.indivs,replace=F)

	#arrange in data frame and order by group
	inds<-data.frame(indivs,groups=indiv.groups)
	inds<-inds[order(indiv.groups),]
	
	#Generate sex vector. This sounds rude.
	inds$sex=sample(c("M","F"),nrow(inds),replace=T)

	###SPATIAL LOCATION OF GROUPS/CLUSTERS####

	#create a dataframe of group locations
	#points end up on a grid in space
	group.id<-seq(1,groups,1)

	group.x=rep(seq(1,groups,1),each=groups)[sample(1:groups^2,groups,replace=F)]
	group.y=rep(seq(1,groups,1),each=groups)[sample(1:groups^2,groups,replace=F)]

	#create dataframe of spatial information
	spat<-data.frame(group.id,group.x,group.y)

	#create distance matrix for groups
	dists<-as.matrix(dist(spat[,2:3]))
	rownames(dists)<-colnames(dists)<-group.id

	#standardise and invert (so 1 is closest and 0.001 is furthest away)
	dists2<-dists/max(dists)
	dists3<-1.001-dists2
	diag(dists3)<-1

	##could try this instead? Neater
	#dists3= 1/(1+dists) # neater, but does not give exactly the same results
	

	inds$x=group.x[inds[,2]]
	inds$y=group.y[inds[,2]]
	#-----------------------------------------------------------------------------------------------------------------

	#####NETWORK STUFF#####
	#create empty network in association matrix form
	net.d<-array(NA,dim=rep(nrow(inds),2))
	colnames(net.d)<-rownames(net.d)<-inds[,1]

	#create network info
	dyads=which(upper.tri(net.d),arr.ind=T)
	dsex=sapply(1:2,function (x) inds[dyads[,x],"sex"])
	dsites=sapply(1:2,function (x) inds[dyads[,x],2])

	distsv=sapply(1:nrow(dsites),function (x) dists3[dsites[x,1],dsites[x,2]])

	samesex=dsex[,1]==dsex[,2]
	samesite=dsites[,1]==dsites[,2]
	ismale=dsex[,1]=="M"|dsex[,2]=="M"
	

	
	#WITHIN GROUP EDGES####

	#FF
	net.d[symmat(samesite&!ismale,dyads)]=sapply(which(samesite&!ismale,dyads),function (x) round(rnbinom(1,size=i.dens*(distsv[x])^d.eff,prob=0.3))+
		round(rnbinom(1,size=i.dens*(distsv[x])^d.eff,prob=0.3)))
	#MF
	net.d[symmat(samesite&ismale&!samesex,dyads)]=sapply(which(samesite&ismale&!samesex),function (x) round(rnbinom(1,size=i.dens*(distsv[x])^d.eff,prob=0.3))+
		round(rnbinom(1,size=(m.i.eff+i.dens)*((distsv[x])^d.eff),prob=0.3)))
	#MM
	net.d[symmat(samesite&ismale&samesex,dyads)]=sapply(which(samesite&ismale&samesex),function (x) round(rnbinom(1,size=(m.i.eff+i.dens)*((distsv[x])^d.eff),prob=0.3))+
		round(rnbinom(1,size=(m.i.eff+i.dens)*((distsv[x])^d.eff),prob=0.3)))
	
	#OUTSIDE GROUP EDGES####

	#FF
	net.d[symmat(!samesite&!ismale,dyads)]=sapply(which(!samesite&!ismale,dyads),function (x) round(rnbinom(1,size=o.dens*(distsv[x])^d.eff,prob=0.3))+
		round(rnbinom(1,size=o.dens*(distsv[x])^d.eff,prob=0.3)))
	#MF
	net.d[symmat(!samesite&ismale&!samesex,dyads)]=sapply(which(!samesite&ismale&!samesex),function (x) round(rnbinom(1,size=o.dens*(distsv[x])^d.eff,prob=0.3))+
		round(rnbinom(1,size=(m.o.eff+o.dens)*((distsv[x])^d.eff),prob=0.3)))
	#MM
	net.d[symmat(!samesite&ismale&samesex,dyads)]=sapply(which(!samesite&ismale&samesex),function (x) round(rnbinom(1,size=(m.o.eff+o.dens)*((distsv[x])^d.eff),prob=0.3))+
		round(rnbinom(1,size=(m.o.eff+o.dens)*((distsv[x])^d.eff),prob=0.3)))	
	
	#SEX HOMOPHILY EFFECT####
	net.d[symmat(samesex,dyads)]= round(net.d[symmat(samesex,dyads)]*sex.eff)

	diag(net.d)=0
	net.d[net.d<0]=0

	pop.dat<-list(ind_data=inds,network=net.d,distmat=dists3)

	return(pop.dat)
}

###NETWORK SAMPLING FUNCTION####		      
		      
networkobs<-function(pop.dat,timesteps,obseff,floaterprob=0.01,probnorm=NA){
	####
	#pop.dat - pop.dat object produced by network_generator
	#timesteps - number of grouping events
	#obseff - number from 0 to 1 representing how much is seen
	#probnorm - a value to normalise associations
	# 		if not given is obtained from the network being considered
	####

	#get objects
	inds=pop.dat$ind_data
	network=pop.dat$network
	distmat=pop.dat$distmat
	
	#data frame of dyads
	dyads=which(upper.tri(network),arr.ind=T)
	assocs=network[dyads]
	names1=as.numeric(row.names(network)[dyads[,1]])
	names2=as.numeric(colnames(network)[dyads[,2]])
	dyads=data.frame(names1,names2,dyads,assocs=assocs)

	#observed contact network - adds noise based on observation effort. 
	#network
	obsnetwork=network
	obsnetwork[symmat(T,as.matrix(dyads[,3:4]))]=sapply(dyads$assoc,function (x) assocnoise(x,obseff))
	
	#generate GBI mat. Select a random individual - check chance of being seen with other
	#indiv based on their associations. Repeat this for any individual added.
	gbid=t(sapply(1:timesteps,function (x) makeevent(inds,dyads,floaterprob,probnorm)))
	gbimat=do.call(rbind,gbid[,1])
	colnames(gbimat)=inds$indivs
	gbigroups=do.call(c,gbid[,2])


	#randomly see only certain groups based on observation effort (could make it so larger groups are more likely to be seen - but, simple for now)
	obsevents=sort(sample(1:timesteps,round(obseff*timesteps)))
	obsgbimat=gbimat[obsevents,]
	obsgbigroups=gbigroups[obsevents]
	#similarly we could throw away individuals from the GBI here based on their sociality - but keeping it simple for now
	
	return(list(truegbimat=gbimat,truegbigroups=gbigroups,obsgbi=obsgbimat,obsgbigroups=obsgbigroups,obsnetwork=obsnetwork))

}

assocnoise<-function(x,obseff){
	if(x==0){
		return(0)
	} else {
		xvec=seq(0,1,length.out=x+1)
		probs=(-(obseff-1)*(obseff*(xvec-1)+1))^(2-(2*xvec))
		obsassoc=sample(0:x,1,1,prob=probs)
		return(obsassoc)
	}
}

makeevent<-function(inds,dyads,floaterprob=0.01,probnorm=NA,verbose=F){
	if(is.na(probnorm)){
		probnorm=max(dyads$assoc)+1
	}
	#randomly pick an individual to start
	gbirow=matrix(0,1,length(inds$indivs))
	colnames(gbirow)=inds$indivs
	seed=sample(inds$indivs,1)
	seedsite=inds$groups[inds$indivs==seed]
	seeddegree=abs(round(rnorm(1,sum(focalpotentials(seed,dyads)$assocs>0),2)))
	gbirow[colnames(gbirow)==seed]=1
	
	if(verbose==T){
		print(paste("seed is",seed))
		print(paste("seed degree is",seeddegree))
	}	
	fassoc=focalassoc(seed,dyads,gbirow,probnorm,floaterprob,T)
	todo=fassoc
	checked=seed
	if(verbose==T){
		print(paste("todo:",paste(todo,collapse=" ")))
		print(paste("checked:",paste(checked,collapse=" ")))
	}	
	gbirow[colnames(gbirow)%in%fassoc]=1
	
	#check the associations of the individuals we are adding until we are no longer adding members
	while(length(todo)>0){

		#randomise order of todo so we're not biased by the order they were added to matrix
		todo=todo[sample(1:length(todo))]
		
		if(verbose==T){
			print(paste("groupsize:",sum(gbirow)))
			print(paste("focal individual:",todo[1]))
		}	
		
		fassoc=focalassoc(todo[1],dyads,gbirow,probnorm,floaterprob)
		gbirow[colnames(gbirow)%in%fassoc]=1		
		checked=c(checked,todo[1])
		todo=c(todo[!todo%in%checked],unique(fassoc[!fassoc%in%checked&!fassoc%in%todo]))
		if(verbose==T){
			print(paste("todo:",paste(todo,collapse=" ")))
			print(paste("checked:",paste(checked,collapse=" ")))
		}	

	}
	return(list(gbirow=gbirow,group=seedsite))	
}

focalpotentials=function(focal,dyads){
	potentials=dyads[dyads$names1==focal|dyads$names2==focal,]
	#swap the names around for readability
	potentials$names2[potentials$names1!=focal]=potentials$names1[potentials$names1!=focal]
	potentials$names1[potentials$names1!=focal]=focal
	return(potentials)
}


focalassoc<-function(focal,dyads,currevent,probnorm,floaterprob=0.01,forcemulti=F){
	####
	# focal = ID of individual whose associations we are using
	# dyads = all possible dyads
	# forcemulti = Force function to return at least one interaction?
	# probnorm = value used for the normalisation of associations probabilities - if not given is set to 
	# 	
	####
	

	
	if(is.na(probnorm)){
		probnorm=max(dyads$assoc)
	}
	
	potentials=focalpotentials(focal,dyads)
	
	eventmembers=colnames(currevent)[currevent==1&(colnames(currevent)!=focal)]
	probs=potentials$assoc/(probnorm+1)
	if(sum(probs)==0){
		return()
	}
	potids=potentials$names2[probs>0]
	#consider these potentials vs existing group members - if they have a zero association with 
	#an individual already in the group reduce the probability of them being in the grouping event
	probs[probs>0][sapply(potids, function (i) {potassocs=focalpotentials(i,dyads);sum(potassocs[potassocs$assocs==0,"names2"]%in%eventmembers)>0})]=floaterprob

	probs=probs^2
	if(forcemulti==T){
		gbirow=rep(0,nrow(potentials))
		while(sum(gbirow)==0){

			gbirow=rbinom(potentials$assoc, 1, probs)
		}
	} else {
		gbirow=rbinom(potentials$assoc, 1,probs )

	}
	return(potentials$names2[gbirow==1])
}




