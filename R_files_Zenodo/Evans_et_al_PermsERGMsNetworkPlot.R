library(igraph)
library(viridis)
source("Evans_et_al_network_generator_functions.R")
library(vegan)
library(ggthemes)

require(parallel)
require(doParallel)
corn=detectCores()-1
cl=registerDoParallel(cores=corn)

donetplot=function(cnet,ewmulti=5,l1=NULL,bline=F,bpad=0.5,tpad=0.1,fw=1){
	cnet=graph_from_adjacency_matrix(cnet,"undirected",weighted=T)
	cnetdf=as_data_frame(cnet)
	cnetdf=cnetdf[order(cnetdf$weight),]
	cnetdf=cnetdf[cnetdf$weight>0.01,]
	cnet=graph_from_data_frame(cnetdf,directed=F,vertices=ctrue$ind_data)

	E(cnet)$samegroup=F
	E(cnet)$group=0
	for(i in unique(V(cnet)$groups)){
		E(cnet)[(V(cnet)$groups==i)%--%(V(cnet)$groups==i)]$samegroup=T
		E(cnet)[(V(cnet)$groups==i)%--%(V(cnet)$groups==i)]$group=i
	}
	cnetdf=as_data_frame(cnet)
	cnetdf=rbind(cnetdf[cnetdf$samegroup==0,],cnetdf[cnetdf$samegroup==1,])
	cnet=graph_from_data_frame(cnetdf,directed=F,vertices=ctrue$ind_data)

	cols=gdocs_pal()(10)
	groups=c(1:max(ctrue$ind_data$groups))
	colvec2=c()
	grvec=c()
	colvec2=c(colvec2,cols[1])
	grvec=c(grvec,groups[1])
	colend=F
	cols=cols[!cols%in%colvec2]
	groups=groups[!groups%in%grvec]
	while(length(groups)>0){
		cdist=ctrue$distmat[,grvec[length(grvec)]]
		cdist=cdist[!names(cdist)%in%grvec]
		closest=as.numeric(names(cdist)[cdist==max(cdist)])[1]
		grvec=c(grvec,closest)
		if(colend){
			colvec2=c(colvec2,cols[length(cols)])
			colend=F
		}else{
			colvec2=c(colvec2,cols[1])
		}
		cols=cols[!cols%in%colvec2]
		groups=groups[!groups%in%grvec]
	}
	lastcol=colvec2[length(colvec2)]
	firstcol=colvec2[1]
	colvec2[1]=lastcol
	colvec2[length(colvec2)]=firstcol
	colvec3=colvec2[match(ctrue$ind_data$groups,grvec)]

	estr=E(cnet)$weight
	colvec1=rgb(colorRamp(c("white","dimgrey"))(estr),maxColorValue=256)
	
	for(i in which(E(cnet)$group==0)){
		colvec1[i]=lighten("darkgrey",((1-E(cnet)$weight[i])*0.25)+0.35)
	}

	lighten <- function(color, factor = 0.5) {
  		col <- col2rgb(color)
  		col <- col + (255 - col)*factor
  		col <- rgb(t(col), maxColorValue=255)
  		col
	}
	for(i in which(E(cnet)$group>0)){
		colvec1[i]=lighten(colvec2[match(E(cnet)$group[i],grvec)],((1-E(cnet)$weight[i])*0.25)+0.3)
	}

	vshapes=c("fsquare","fcircle")[(V(cnet)$sex=="F")+1]
	
	if(is.null(l1)){
		starts=ctrue$ind_data[,c("x","y")]
		starts$x=jitter(starts$x,amount=0.1)
		starts$y=jitter(starts$y,amount=0.1)
		l1=layout_with_graphopt(cnet,start=as.matrix(starts),mass=30,charge=0.01,niter = 100,max.sa.movement=0.01,spring.constant=5)
	}

	par(mar=c(0,0,0,0))
	plot(cnet,layout=l1,vertex.size=max(as.vector(l1))*2,edge.width=(estr)*ewmulti,edge.color=colvec1,vertex.color=colvec3,vertex.label="",vertex.shape=vshapes,rescale=F,asp=0,
	xlim=c(min(l1[,1])-0.5,max(l1[,1])+0.5),ylim=c(min(l1[,2])-bpad,max(l1[,2])+tpad),frame=F,vertex.frame.width=fw,vertex.frame.color="black")
	if(bline){
		abline(h=round(min(l1[,2])-bpad),lwd=2)
	}
	return(l1)

}

parameters=read.csv("parametersets.csv")
parameters$intfreq=NULL
#generate networks for each parameter set
result1=foreach(job=1:nrow(parameters))%dopar%{
	source("network_generator_functions.R")
	do_networksim = function(groups,mean.group.size,max.group.size,d.eff,i.dens,o.dens,m.i.eff,m.o.eff,sex.eff,
		
		obs.eff.v,timesteps,floaterprob,probnorm,
		nreps,exportdir){
		obs.eff.v=c(1,obs.eff.v)
	
		allresults=list()

		for (rep in 1:nreps){

			#generate a network
			simulated.networks=network.generator(groups,mean.group.size,max.group.size,d.eff,o.dens,i.dens,sex.eff,m.i.eff,m.o.eff)
			#for each observation effort
			obs.sim.networks=list()
			for(obs.eff.c in obs.eff.v){

				parvec=paste(d.eff,i.dens,o.dens,m.i.eff,sex.eff,obs.eff.c,timesteps,sep="_")#parameters of interest
				#observe the network
				obs.sim.networks=c(obs.sim.networks,list(networkobs(simulated.networks,timesteps, obs.eff.c ,floaterprob,probnorm)))
			}
		
			allresults=c(allresults,list(list(simulated.networks=simulated.networks,obs.sim.networks=obs.sim.networks)))
		}
		return(allresults)
	}

	
	library(vegan)
	library(igraph)
	library(asnipe)
	currpar=parameters[job,]

	obs.effvec=as.numeric(strsplit(as.character(currpar$obs.eff)," ")[[1]])

	#version of currpar with all obs.effs for checking for completion
	currpar2=do.call(rbind,(lapply(obs.effvec,function (i) {
		cpar1=currpar
		cpar1$obs.eff=i
		return(cpar1)
	})))

	currpar3=currpar
	currpar3$nreps=10
	currpar3=as.list(currpar3)
	currpar3$obs.eff=rev(obs.effvec)

	results=do.call(do_networksim,currpar3)
	
	for(j in 1:length(results)){
		cresults=results[[j]]
		ctrue=cresults$simulated.networks
		callobs=cresults$obs.sim.networks

		rescale=function(x){(x-min(x))/(max(x)-min(x))}

		ctruenet=ctrue$network
		ctruenet2=rescale(ctruenet)


		for(i in 2:4){
			currobs=callobs[[i]]

			cobsnet=currobs$obsnetwork
			cobsnet2=rescale(cobsnet)

			cgbinet=asnipe::get_network(currobs$obsgbi)
			cgbinet2=rescale(cgbinet)

			m1=mantel(ctruenet2,cobsnet2)
			m2=mantel(ctruenet2,cgbinet2)
			
			results[[j]]$obs.sim.networks[[i]]$mtests=list(m1=m1,m2=m2)
		}
		
	}
	results
	
}

#chosen examples
chosen=c(
which(parameters$i.dens==0.8&parameters$o.dens==0.4&parameters$d.eff==4&parameters$m.i.eff==0),
which(parameters$i.dens==1.2&parameters$o.dens==0.1&parameters$d.eff==4&parameters$m.i.eff==0&parameters$sex.eff==1),
which(parameters$i.dens==0.4&parameters$o.dens==0.4&parameters$d.eff==4&parameters$m.i.eff==0&parameters$sex.eff==1)
)

for(j in chosen){
cresults=result1[[j]][[1]]
currpars=parameters[j,c(4,5,6,7,9)]
ctrue=cresults$simulated.networks
callobs=cresults$obs.sim.networks

rescale=function(x){(x-min(x))/(max(x)-min(x))}

ctruenet=ctrue$network
ctruenet2=rescale(ctruenet)

currpars2=unlist(lapply(1:length(currpars),function(x){c(names(currpars)[x],currpars[x])}))

tiff(paste(c("netplot",currpars2,".tiff"),collapse="_"),units="mm",width=177,height=120,res=600)
layout(matrix (c(9,8,9,9,1,8,2,3,1,8,4,5,1,8,6,7),nrow=4,byrow=T),widths=c(1,0.05,0.5,0.5),heights=c(0.1,1,1,1),respect=F)
par(mar=c(0,0,0,0))
l=donetplot(ctruenet2,bline=F,tpad=1,bpad=1)
for(i in 2:4){
	currobs=callobs[[i]]

	cobsnet=currobs$obsnetwork
	cobsnet2=rescale(cobsnet)

	cgbinet=asnipe::get_network(currobs$obsgbi)
	cgbinet2=rescale(cgbinet)

	m1=mantel(ctruenet2,cobsnet2)
	m2=mantel(ctruenet2,cgbinet2)

	nl=donetplot(cobsnet2,2.5,l,bline=(i!=4),bpad=1.5,fw=0.5)
	text(round(median(l[,1]))-0.5,round(min(l[,2])-1),paste0("r = ",round(m1$statistic,2)),offset=0,adj=c(0.5,0.5))
	nl=donetplot(cgbinet2,2.5,l,bline=(i!=4),bpad=1.5,fw=0.5)
	text(round(median(l[,1]))-0.5,round(min(l[,2])-1),paste0("r = ",round(m2$statistic,2)),offset=0,adj=c(0.5,0.5))
}
plot(0,0,type="n",axes=F)
text(0.25,c(0.4,-0.3,-1)-0.025,c("a)","b)","c)"),cex=1.25)
plot(0,0,type="n",axes=F)
text(-0.6,-0.7,"'True' network",xpd=T,cex=1,offset=0,adj=c(0.5,0))
text(0.3,-0.7,"Dyad-based networks",xpd=T,cex=1,offset=0,adj=c(0.5,0))
text(0.8,-0.7,"Grouping-based networks",xpd=T,cex=1,offset=0,adj=c(0.5,0))
dev.off()
}




