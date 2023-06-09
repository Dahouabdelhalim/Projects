# Code for inferring displacement events from temporal data, from:
# Evans JC, Devost I, Jones TB, Morand-Ferron J. 2018 Inferring dominance interactions from automatically recorded temporal data. Ethology 124, 188â€“195. 
# (doi:10.1111/eth.12720)


detectdisplacement = function(arrivetime,departtime=NA,ids,displacetime,remaintime,totaltime=NA,difftime=NA,sites=NA,divisions=NA){
	require(EloRating)

	##
	#arrivetime = time of arrival in seconds 
	#departtime = time of departure in seconds -not required if difftime/totaltime is provided
	#
	#ids = vector of corresponding ids
	#displace time = difference between one individual leaving and another arriving to be considered a displacement
	#
	#totaltime = total time individual spends at antennae. Will be calculated automatically if not provided
	#difftime = difference in time between an individual departing and an individual arriving. Will be calculated automatically if not provided
	#		NOTE: It is reccomended that this is calculated prior, and visits by the same individual within a short period are merged.
	#
	#sites = vector of sites - if not provided it will be assumed that all individuals  are present at the same site
	#divisons = vector of by which to divide the data. This can be less computationally intense than calculating a whole site at once. 
	#		For example seperating data into days,sub-sites, or using gmmevents to divide into bursts of activity.
	#
	#Returns: A dominance matrix based, list of detected interactions.
	##

	if(((anyNA(arrivetime)&length(arrivetime)==1)|(anyNA(departtime)&length(departtime)==1))&(anyNA(difftime)&length(difftime)==1)){
		stop("Unable to obtain time differences")
	}

	#if totaltime is not available, calculate it here

	if((length(totaltime==1)&anyNA(totaltime))){
		totaltime=departtime-arrivetime
	}

	#if sites are not provided, make one here
	if((length(sites==1)&anyNA(sites))){
		sites=rep(1,length(arrivetime))
	}

	#if divisions are not provided, make one here
	if((length(divisions==1)&anyNA(divisions))){
		divisions=rep(1,length(arrivetime))
	}

	#if difftime is not available, calculate it here
	if((length(difftime)==1&anyNA(difftime))){
		difftime=rep(NA,length(arrivetime))
		for(site in unique(sites)){
			for(d in unique(divisions[sites==site])){
				currind=which(sites==site&divisions==d)
				if(length(currind)>1){
					currind=currind[2:length(currind)]
					difftime[currind]=sapply(currind,function (x) arrivetime[x]-departtime[x-1])
				}
			}
		}
	}

	winlossdate<-data.frame(matrix(NA,nrow=0,ncol=4))
	for(sn in unique(sites)){
		for(burst in unique(divisions[sites==sn])){
			tempselect=divisions==burst&sites==sn
			difftime2=difftime[tempselect]
			totaltime2=totaltime[tempselect]
			ids2=ids[tempselect]
			arrivetime2=arrivetime[tempselect]
			pwin<-which(difftime2<displacetime&!is.na(difftime2)&totaltime2>=remaintime)
			if(length(pwin)>0){
				dyads<-lapply(pwin,function (i) (c(as.character(ids2[i]),as.character(ids2[i-1]))))
				diffid=which(unlist(lapply(dyads,function (x) x[1] != x[2])))
				pwin=pwin[diffid]
				dyads=dyads[diffid]
				if(is.list(dyads)&length(dyads)>0){		
					wld<-data.frame(do.call(rbind,dyads),arrivetime2[pwin],sn)
					winlossdate<-rbind(winlossdate,wld)
				} else if (length(dyads)>0) {
					wld<-data.frame(dyads,arrivetime[pwin],sn)
					winlossdate<-rbind(winlossdate,wld)
				}
			}
		}
		
	}

	names(winlossdate)<-c("win","loss","date","Site")
	winlossdate$dt<-winlossdate$date
	winlossdate$date<-as.Date(as.POSIXct(winlossdate$date,origin="1970-01-01"))
	dommat<-data.frame(matrix(nrow=0,ncol=5))
	if(nrow(winlossdate)==0){
		#insufficient interactions!
		return(list("dominance matrix"=dommat,"interactions"=winlossdate))
	}

	for (i in unique(winlossdate$Site)){
		wld<-winlossdate[winlossdate$Site==i,]
		if(nrow(wld)>5){
			SEQ <- elo.seq(winner=wld$win, loser=wld$loss, Date=wld$date,progressbar=F)
			mat <- creatematrix(SEQ)
			dommat2<-DS(mat)
			dommat2$Site<-rep(i,nrow(dommat2))
			dommat2$ranked<-rank(dommat2[,2])
			dommat<-rbind(dommat,dommat2)
		}else{
			dommat2<-data.frame(ID=unique(ids[sites==i]))
			dommat2$DS<-rep(NA,nrow(dommat2))
			dommat2$normDS<-dommat2$DS
			dommat2$Site<-rep(i,nrow(dommat2))
			dommat2$ranked<-rank(dommat2[,2])						
			dommat<-rbind(dommat,dommat2)
			}
		}
	dommat<-dommat[order(dommat$Site),]
	return(list("dominance_matrix"=dommat,"interactions"=winlossdate))
}

comparedom<-function(dommat1,dommat2){
	#function to calculate the correlation between two dominance hierarchies
	#Where dommat1 and dommat2 are both dataframes with the following columns
	#ID - individual ID
	#Site - site individual ID is at.
	#DS - David's score, used for ranking individuals.
	#
	#Outputs - the comparable dominance matrix, phi (The correlation between the ranks of the hierarchies, controlling for comparable individuals),
	#		Rho value, p value, a per site phi, a per site rho, a per site p value.
	op <- options(warn = (-1)) # suppress warnings 
	alldoom<-dommat1[dommat1$ID%in%dommat2$ID,]
	alldoom<-alldoom[order(alldoom[,1],decreasing=T),]
	alldoom2<-data.frame(matrix(nrow=0,ncol=3))

	for (i in unique(alldoom$Site)){
		alldom<-data.frame(matrix(NA,ncol=0,nrow=length(alldoom$ID[alldoom$Site==i])))
		alldom$id<-alldoom$ID[alldoom$Site==i]
		alldom$DS<-alldoom$DS[alldoom$Site==i]
		alldom$ranked<-rank(alldoom$DS[alldoom$Site==i])
		alldom$Site<-rep(i,nrow(alldom))
		alldoom2<-rbind(alldoom2,alldom)
	}

	names(alldoom2)[1]<-"ID"

	dommat3<-dommat2[dommat2$ID%in%alldoom2$ID,]

	dommat3$ranked<-unlist(sapply(unique(dommat3$Site),function (x) rank(dommat3$DS[dommat3$Site==x])))

	alldoom3<-alldoom2[order(alldoom2$Site),]

	doom4<-merge(alldoom3,dommat3,by="ID")
	doom4<-doom4[order(doom4$Site.x),]
	persitephi=sapply(unique(doom4$Site.x),function (x) as.numeric(cor.test(doom4$ranked.x[doom4$Site.x==x],doom4$ranked.y[doom4$Site.x==x],method="spearman")[[4]])*((nrow(doom4[doom4$Site.x==x,])/nrow(dommat1[dommat1$Site==x,]))))
	persiterho=sapply(unique(doom4$Site.x),function (x) as.numeric(cor.test(doom4$ranked.x[doom4$Site.x==x],doom4$ranked.y[doom4$Site.x==x],method="spearman")[[4]]))
	persitep=sapply(unique(doom4$Site.x),function (x) as.numeric(cor.test(doom4$ranked.x[doom4$Site.x==x],doom4$ranked.y[doom4$Site.x==x],method="spearman")[[3]]))

	cortest=cor.test(doom4$ranked.x,doom4$ranked.y,method="spearman")
	rho=cortest[[4]]
	pvalue=cortest[[3]]
	
	r2z=as.numeric(fisher.r2z(persiterho))
	r2z[r2z==Inf]=2.6467
	r2z[r2z==-Inf]=-2.6467
	tpval=t.test(r2z)$p.value
	meanz=mean(r2z)
	phi2=meanz*pmin((nrow(doom4)/nrow(dommat1)),1)

	options(op)
	


	return(list(dominance_matrix=doom4,phi=phi2,phi2pval=as.numeric(tpval),overall_rho=rho,scaled_rho=as.numeric(rho*pmin((nrow(doom4)/nrow(dommat1)),1)),rho_pval=pvalue,per_site_rho=persiterho,per_site_rho_scaled=persitephi,per_site_p=persitep))


}

fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }

