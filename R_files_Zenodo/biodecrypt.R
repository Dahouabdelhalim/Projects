biodecrypt.view<-function(mat,id,alpha=NULL, fraction=0.95, partCount=10, map=NULL,main=NULL, xlim=NULL,ylim=NULL,clipToCoast="terrestrial",cex0=0.2,cex1=0.25){
	res<-NULL
	taxa<-max(id)
	if(is.null(alpha)){
		alpha=rep(8,taxa)
	}
	if(is.null(xlim)){
		xlim<-c(min(mat[,1]),max(mat[,1]))
	}
	if(is.null(ylim)){
		ylim<-c(min(mat[,2]),max(mat[,2]))
	}
	plot(cbind(xlim,ylim),type="n")
	if(!is.null(map)){
		plot(map,add=T)
	}
	points(mat,col="grey",cex=cex0)
	points(mat,col=id,cex=cex1)
	hulls<-list()
	areas<-NULL
	oldw <- getOption("warn")
		options(warn = -1)
		for (sp in 1:taxa){
			taxsp<-which(id==sp)
			hulla<-mat[taxsp,]
			hullas<-hulla[!duplicated(hulla), ]
			hull<-rangeBuilder::getDynamicAlphaHull(hullas[,c(1,2)],fraction = fraction, partCount = partCount, buff = 0, clipToCoast=clipToCoast,initialAlpha=alpha[sp])[[1]]
			plot(hull,border=sp,add=T)
			hulls[sp]<-hull 
			areas[sp]<-area(hulls[[sp]])/1000000
		}
		intersect<-matrix(NA,taxa,taxa)
		sympatry<-intersect
		for(k in 1:(taxa-1)){
			for(c in 2:taxa){
				inter<-raster::intersect(hulls[[k]],hulls[[c]])
				if(!is.null(inter)){
					intersect[k,c]<-(raster::area(inter)/1000000)
					sympatry[k,c]<-(raster::area(inter)/1000000)/(areas[c]+areas[k]-area(inter)/1000000)
					intersect[c,k]<-intersect[k,c]
					sympatry[c,k]<-sympatry[k,c]
				}
				if(is.null(inter)){
					intersect[k,c]<-0
					sympatry[k,c]<-0
					intersect[c,k]<-intersect[k,c]
					sympatry[c,k]<-sympatry[k,c]
				}		
			}
		}	
		options(warn = oldw)
		res$areas<-areas
		res$intersections<-intersect
		res$sympatry<-sympatry
		return(res)
}





biodecrypt<-function(mat,id,alpha=NULL,ratio=2.5,buffer=90,fraction=0.95, partCount = 10, clipToCoast="terrestrial", checkdist=T, proj = "+proj=longlat +datum=WGS84", minimum=7, map=NULL,xlim=NULL,ylim=NULL,main=NULL){
	res<-NULL
	res$type<-"sep"
	borders<-NULL
	taxa<-max(id)
	colnames(mat)<-c("Long","Lat")
	distances<-matrix(0,nrow(mat), taxa)
	distances2<-matrix(0,nrow(mat), taxa)
	if(is.null(alpha)){
		alpha=rep(8,taxa)
	}
	if(is.null(xlim)){
		xlim<-c(min(mat[,1]),max(mat[,1]))
	}
	if(is.null(ylim)){
		ylim<-c(min(mat[,2]),max(mat[,2]))
	}
	#plot(cbind(xlim,ylim),type="n",main=main)
	#if(!is.null(map)){
	#	plot(map,add=T)
	#}
	vectab<-matrix(0,nrow(mat),taxa+1)
	#points(mat,col="grey",cex=0.2)
	#points(mat,col=id,cex=0.25)
	hulls<-list()
	areas<-NULL
	oldw <- getOption("warn")
	options(warn = -1)
	for (spec in 1:taxa){
		#spec<-1
		taxsp<-which(id==spec)
		hulla<-mat[taxsp,]
		hullas<-hulla[!duplicated(hulla), ]
		if(nrow(hullas)>=minimum){
			hull<-rangeBuilder::getDynamicAlphaHull(hullas[,c(1,2)],fraction=fraction, partCount = partCount, buff = 0, clipToCoast=clipToCoast, proj=proj, initialAlpha=alpha[spec])[[1]]
	  	}
		if(nrow(hullas)<minimum){
			coordinates(hullas) <- ~Long+Lat
			hull <- rgeos::gConvexHull(hullas)
			proj4string(hull) <- proj 
			data(gshhs, envir = environment())
			if (clipToCoast != "no") {
				data(gshhs, envir = environment())
				gshhs <- sp::spTransform(gshhs, CRS(proj4string(hull)))
				if (clipToCoast == "terrestrial") {
					hull <- rgeos::gIntersection(hull, gshhs)
				}
			  	else {
			    	hull <- rgeos::gDifference(hull, gshhs)
				}
			}
		}
		hulls[spec]<-hull 
		areas[spec]<-raster::area(hulls[[spec]])/1000000
		#plot(hull,border=spec,add=T)
		vectab[prevR::point.in.SpatialPolygons(mat[,1], mat[,2],hull),spec]<-1	
		fuo<-which(id==spec & vectab[,spec]==0)
		fuori<-mat[fuo,]
		distneed<-which(vectab[,spec]==0 & id==0)
		if(length(distneed>0)){
			needdist<-mat[which(vectab[,spec]==0 & id==0),] 		  
			disthull<-geosphere::dist2Line(needdist,hull)[,1]/1000
		}
		if(length(fuo)>1 & length(which(vectab[,spec]==0 & id==0))==1){
			distpoints<-(geosphere::distm(rbind(fuori, needdist), fun=distGeo)/1000)
			distpointsneed<-distpoints[(length(fuo)+1):nrow(distpoints),(1:length(fuo))]
			distot<-c(disthull, distpointsneed)
			distances[which(vectab[,spec]==0 & id==0),spec]<-min(distot)
		}
		if(length(fuo)>1 & length(which(vectab[,spec]==0 & id==0))>1){
			distpoints<-(geosphere::distm(rbind(fuori, needdist), fun=distGeo)/1000)
			distpointsneed<-distpoints[(length(fuo)+1):nrow(distpoints),(1:length(fuo))]
			distot<-cbind(disthull, distpointsneed)
			distmin<-apply(distot, 1, FUN=min)
			distances[which(vectab[,spec]==0 & id==0),spec]<-distmin
		}
		if(length(fuo)==1 & length(which(vectab[,spec]==0 & id==0))>1 ){
			distpoints<-(geosphere::distm(rbind(fuori, needdist), fun=distGeo)/1000)
			distpointsneed<-distpoints[2:nrow(distpoints),1]
			distot<-cbind(disthull, distpointsneed)
			distmin<-apply(distot, 1, FUN=min)
			distances[which(vectab[,spec]==0 & id==0),spec]<-distmin	
		}
		if(length(fuo)==1 & length(which(vectab[,spec]==0 & id==0))==1 ){
			distpoints<-(geosphere::distm(rbind(fuori, needdist), fun=distGeo)/1000)
			distpointsneed<-distpoints[2:nrow(distpoints),1]
			distot<-c(disthull, distpointsneed)			
			distances[which(vectab[,spec]==0 & id==0),spec]<-min(distot)	
		}
		if(length(fuo)==0){
			distances[which(vectab[,spec]==0 & id==0),spec]<-disthull
		}	 
	}	
	vectab[,ncol(vectab)]<-rowSums(vectab[,1:taxa])
	id2<-id
	uncertain1<-which(vectab[,ncol(vectab)]>1 & id==0)
	uncertain2<-which(vectab[,ncol(vectab)]==0 & id==0)
	inside<-which(vectab[,ncol(vectab)]==1 & id==0)
	if (length(uncertain2)>1){
		distancesunc<-distances[uncertain2,]
		order<-matrix(NA,length(uncertain2),taxa)
		for (unc2 in 1:length(uncertain2)){
			wh<-uncertain2[unc2]
			order[unc2,1:taxa]<-c(1:taxa)[order(distancesunc[unc2,1:taxa])]
			}
		attribution<-order[,1]	
		for (k in 1:length(uncertain2)){
			dist<-distancesunc[k,]
			ordereddist<-dist[order(dist)]
			if(ordereddist[2]>buffer & (ordereddist[2]/ordereddist[1])>ratio){
				id2[uncertain2[k]]<-attribution[k]
			}
		}
	}
	if (length(uncertain2)==1){
		distancesunc<-distances[uncertain2,]
		order<-matrix(NA,1,taxa)
		order[1,1:taxa]<-c(1:taxa)[order(distancesunc[1:taxa])]
		attribution<-order[,1]
		ordereddist<-distancesunc[order(distancesunc)]
		if(ordereddist[2]>buffer & (ordereddist[2]/ordereddist[1])>ratio){
				id2[uncertain2]<-attribution

		}
	
	}
	if (length(inside)>1){
		distancesunc<-distances[inside,]
		for (k in 1: length(inside)){
			#k<-50
			attr<-which(vectab[inside[k],1:taxa]==1)
			diste<-distancesunc[k,]
			diste<-diste[-attr]
			if(min(diste)>buffer){
				id2[inside[k]]<-attr
			}
		}
	}	
	if (length(inside)==1){
		distancesunc<-distances[inside,]
		attr<-which(vectab[inside,1:taxa]==1)
		diste<-distancesunc[-attr]
		if(min(diste)>buffer){
			id2[inside]<-attr	
		}	
	}
	if(checkdist){
	  check<-which(id2>0 & id==0)
	  if(length(check>0)){
	    for(ch in 1:length(check)){
	 	attrp<-id2[check[ch]]
	      dist1<-geosphere::distGeo(mat[check[ch],], mat[which(id>0),])
	      minimum<-aggregate(dist1~id[which(id>0)],FUN="min")
	      mini<-which(minimum[,2]==min(minimum[,2]))
	      if(attrp!=mini){
	        id2[check[ch]]<-0
	      }
	    }
	  }
	}
	intersect<-matrix(NA,taxa,taxa)
	sympatry<-intersect
	for(k in 1:(taxa-1)){
		for(c in 2:taxa){
			inter<-raster::intersect(hulls[[k]],hulls[[c]])
			if(!is.null(inter)){
				intersect[k,c]<-(raster::area(inter)/1000000)
				sympatry[k,c]<-(raster::area(inter)/1000000)/(areas[c]+areas[k]-area(inter)/1000000)
				intersect[c,k]<-intersect[k,c]
				sympatry[c,k]<-sympatry[k,c]
			}
			if(is.null(inter)){
				intersect[k,c]<-0
				sympatry[k,c]<-0
				intersect[c,k]<-intersect[k,c]
				sympatry[c,k]<-sympatry[k,c]
			}		
		}
	}	
	options(warn = oldw)
	res$areas<-areas
	res$intersections<-intersect
	res$sympatry<-sympatry
	res$NUR<-(length(which(id2==0))/length(which(id==0)))*100
	res$table<-cbind(mat,id2,id)
	return(res)
	points(mat,col=id2,cex=0.5)
}


biodecrypt.cross<-function(mat,id,alpha=NULL,ratio=2.5,buffer=90,fraction=0.95, partCount = 10, checkdist=T, clipToCoast="terrestrial", proj = "+proj=longlat +datum=WGS84", minimum=7,map=NULL,xlim=NULL,ylim=NULL,main=NULL,runs=10,test=T){
	res<-NULL
	res$type<-"cross"
	taxa<-max(id)
	idr<-id[which(id>0)]
	matrix<-mat[which(id>0),]
	if(is.null(alpha)){
		alpha=rep(8,taxa)
	}
	q<-aggregate(rep(1,length(idr))~idr,FUN="sum")
	min<-which(q[,2]==min(q[,2]))
	if(min(q[,2])<4){
		stop("The cross validation procedure requires a minimum of 4 distinct points")
	}
	if(test){
	test_run<-biodecrypt(mat,id,alpha=alpha,ratio=ratio,buffer=buffer,map=map,fraction=fraction, partCount=partCount, checkdist=checkdist, clipToCoast=clipToCoast, proj = proj, xlim=xlim,ylim=ylim,main=main)
	res$NUR<-test_run$NUR
	res$areas<-test_run$areas
	res$intersections<-test_run$intersections
	res$sympatry<-test_run$sympatry
	res$table<-test_run$table
	}
	species<-unique(idr)
	which_sp<-NULL
	how_may<-NULL
	order<-NULL
	for (sp in 1:length(species)){
		which_sp[[sp]]<-which(idr==species[sp])[sample(1:length(which(idr==species[sp])))]
		how_may[[sp]]<-length(which_sp[[sp]])/runs
		order<-c(order,which_sp[[sp]])
		}
	matrixnewbs<-cbind(matrix[order,],idr[order])
	attr<-rep(NA, nrow(matrix))
	for (giro in 1:runs){
		first<-NULL
		last<-NULL
		memberbs<-matrixnewbs[,3]
		for (spe in 1:length(species)){
			#spe<-2
			first[spe]<-round(1+((giro-1)*how_may[spe]))
			last[spe]<-round(giro*how_may[spe],0)
			qualitot<-which(matrixnewbs[,3]==species[spe])
			if(first[spe]<=last[spe]){
			memberbs[qualitot[first[spe]:last[spe]]]<-0
			if(giro==runs){
			  memberbs[qualitot[first[spe]:qualitot[length(qualitot)]]]<-0
			}
			}
		}
		finalebs<-biodecrypt(matrixnewbs[,c(1,2)], memberbs, alpha=alpha, minimum=minimum, ratio=ratio,buffer=buffer,fraction=fraction, partCount=partCount, checkdist=checkdist,clipToCoast=clipToCoast, proj = proj,map=map,xlim=xlim,ylim=ylim,main=main) 
		attri<-which(memberbs==0)
		attr[attri]<-finalebs$table[attri,3]
	}
	res$cross<-cbind(matrixnewbs[,3],attr,rep(0,nrow(matrix)),rep(0,nrow(matrix)), matrixnewbs[,1:2])
	colnames(res$cross)<-c("original","predicted","MIR","NIR","Long","Lat")
	res$cross[which(res$cross[,2]>0 & res$cross[,1]!=res$cross[,2]),3]<-1
	res$cross[which(res$cross[,2]==0),4]<-1
	res$MIR<-(sum(res$cross[,3])/nrow(res$cross))*100
	res$NIR<-(sum(res$cross[,4])/nrow(res$cross))*100
	return(res)
}




biodecrypt.wrap<-function(mat,id,alpha=c(1,5,10,15),alphamat=NULL,ratio=c(2,3,4,5),buffer=c(0,40,80,120,160),fraction=0.95, partCount=10, checkdist=T, clipToCoast="terrestrial", proj="+proj=longlat +datum=WGS84", minimum=7, map=NULL,xlim=NULL,ylim=NULL,main=NULL,save=T,name="res_cross.txt",runs=10){
	res<-NULL
	taxa<-max(id)
	al<-length(alpha)
	if(!is.null(alphamat)){
		al<-ncol(alphamat)
	}
	rat<-length(ratio)
	buf<-length(buffer)
	res_cross<-matrix(NA,al*rat*buf,6)
	if(!is.null(alphamat)){
		res_cross<-matrix(NA,al*rat*buf,(5+nrow(alphamat)))
	}

	riga<-1
	for(alphav in 1:al){
		for(ratiov in 1:rat){
			for (bufferv in 1:buf){
				alphause<-rep(alpha[alphav],taxa)
				if(!is.null(alphamat)){
					alphause<-alphamat[,alphav]
					addcol<-ncol(alphamat)
				}
				print(c(alphav,ratiov,bufferv))
				cross<-biodecrypt.cross(mat, id, ratio=ratio[ratiov],buffer=buffer[bufferv],alpha=alphause, checkdist=checkdist,map=map, fraction=fraction, partCount=partCount, clipToCoast=clipToCoast, proj = proj,runs=runs, test=T) 
				res_cross[riga,4]<-ratio[ratiov]
				res_cross[riga,5]<-buffer[bufferv]
				if(is.null(alphamat)){
					res_cross[riga,6]<-alphause[1]
					addcol<-1
				}
				if(!is.null(alphamat)){
					res_cross[riga,6:ncol(res_cross)]<-alphause
				}
				res_cross[riga,1]<-cross$MIR
				res_cross[riga,2]<-cross$NIR
				res_cross[riga,3]<-cross$NUR
				if(save){
					write.table(res_cross,name)
				}
				riga<-riga+1
			}
		}
	}
colnames(res_cross)<-c("MIR","NIR","NUR","ratio","buffer",rep("alpha",addcol))
res$table<-res_cross
return(res)
}


biodecrypt.optimise<-function(tab,coef=c(2,1,1), penalty=10){
  	res<-NULL
	val<-tab[,1]^coef[1]+tab[,2]^coef[2]+tab[,3]^coef[3]
	mini<-min(val)
	quali<-which((val-mini)<penalty)
	if(length(quali)>1){
		tot<-sum(1/val[quali])
		weight<-tab[quali,(4:ncol(tab))]/val[quali]
		res$ratio<-(sum(weight[,1]))/tot
		res$buffer<-(sum(weight[,2]))/tot
		res$MIR<-mean(tab[quali,1])
		res$NIR<-mean(tab[quali,2])
		res$NUR<-mean(tab[quali,3])
		quantialpha<-ncol(weight)-2
		res$alpha<-NULL
			for(q in 1:quantialpha){
			res$alpha[q]<-(sum(weight[,(2+q)]))/tot
		}
	}
	if(length(quali)==1){
		res$ratio<-tab[quali,4]
		res$buffer<-tab[quali,5]
		res$MIR<-tab[quali,1]
		res$NIR<-tab[quali,2]
		res$NUR<-tab[quali,3]
		quantialpha<-ncol(tab)-5
		res$alpha<-NULL
			for(q in 1:quantialpha){
			res$alpha[q]<-tab[quali,(5+q)]	
		}
	}
	return(res)
}



plot.biodecrypt<-function(x,minsize=0.3,pchid=1,cexid=0.1,square=0.001,col=c("red","darkgreen","blue","purple"), attributed=c("fade","points"), NUR="black", fading=50){
	if(x$type=="sep"){
		data<-as.data.frame(x$table)
		newcol<-as.character(paste(data[,1],data[,2],sep="_"))
		data<-cbind(data,newcol)
		quanti<-aggregate(rep(1,nrow(data))~ data[,5],FUN="sum")
		quanti<-quanti[which(quanti[,2]>1),]
		zeri<-which(data[,5]%in%quanti[,1] & data[,4]==0)
		if(length(zeri)>0){
			data2<-data[-zeri,]}else{
			data2<-data
		}
		matcol<-cbind(matrix(1,nrow(data2),2),matrix(NA, nrow(data2),3))
		for(sp in 1:max(as.numeric(data2[,3]))){
			#sp<-2
			which<-which(as.numeric(data2[,3])==sp)
			matcol[which,3]<-as.vector(col2rgb(col[sp]))[1]
			matcol[which,4]<-as.vector(col2rgb(col[sp]))[2]
			matcol[which,5]<-as.vector(col2rgb(col[sp]))[3]
		}
		if(attributed=="fade"){
			lower<-which(data2[,4]==0 & data2[,3]>0)
			for (c in 1:length(lower)){
				for(colo in 3:5){
					matcol[lower[c],colo]<-matcol[lower[c],colo]+((255-matcol[lower[c],colo])/(100/fading))
				}
			}
		}
		matcol[which(is.na(matcol[,3])),3:5]<-col2rgb(NUR)
		recluster.plot.pie(data2[,1],data2[,2],mat=as.matrix(matcol),square=square,minsize=minsize,proportional=F,add=T)
		if(attributed=="points"){
		lat<-data2[,2]
		long<-data2[,1]
		latsqo<-floor(lat/square)*square
		longsqo<-floor(long/square)*square
		newcoord0<-cbind(aggregate(long~longsqo+latsqo, 	FUN="mean"),aggregate(lat~longsqo+latsqo, FUN="mean")[,3]) 
		newcoord0[,5]<-paste(newcoord0[,1],newcoord0[,2])
	lat<-data2[which(data2[,4]>0),2]
	long<-data2[which(data2[,4]>0),1]
	latsqo<-floor(lat/square)*square
	longsqo<-floor(long/square)*square
	newcoord1<-cbind(aggregate(long~longsqo+latsqo, FUN="mean"),aggregate(lat~longsqo+latsqo, FUN="mean")[,3])
	newcoord1[,5]<-paste(newcoord1[,1],newcoord1[,2])
	newcoord3<-newcoord0[which(newcoord0[,5] %in% newcoord1[,5]),]
	points(newcoord3[,3:4],cex=cexid,pch=pchid)
	}
	}
	if(x$type=="cross"){
		data2<-x$cross
		matcol<-cbind(data2[,c(5,6)],matrix(NA,nrow(data2),3))
		for(sp in 1:max(data2[,1])){
			which<-which(data2[,2]==sp)
			matcol[which,3]<-as.vector(col2rgb(col[sp]))[1]
			matcol[which,4]<-as.vector(col2rgb(col[sp]))[2]
			matcol[which,5]<-as.vector(col2rgb(col[sp]))[3]
		}
		matcol[which(data2[,3]==1),3:5]<-col2rgb("black")
		matcol[which(data2[,4]==1),3:5]<-col2rgb("white")
		matcol<-matcol[complete.cases(matcol),]
		recluster.plot.pie(matcol[,1],matcol[,2],mat=as.matrix(matcol),square=square,minsize=minsize,proportional=F,add=T)	
	}
}


recluster.plot.pie<-function(long, lat, mat=NULL, distance=NULL, loc=NULL, areas=NULL, square=2,map=NULL,add=F,minsize=NULL,proportional=T,xlim=NULL,ylim=NULL,main=NULL,xlab=NULL,ylab=NULL,...){
	if(is.null(mat) & is.null(distance)){
 		stop("A distance matrix or a colour matrix from recluster.col must be provided")
	}
	if(is.null(loc)){
		if(is.null(areas)){
			areas<-rep(1,length(lat))
		}
		latsq<-floor(lat/square)*square
		longsq<-floor(long/square)*square
		newcoord<-cbind(aggregate(long~longsq+latsq+areas, FUN="mean"),aggregate(lat~longsq+latsq+areas, FUN="mean")[,4])
		for (k in 1:nrow(newcoord)){
			quali1<-c(1:length(lat))[which(latsq==newcoord[k,2])]
			quali2<-quali1[which(longsq[quali1]==newcoord[k,1])]
			quali3<-quali2[which(areas[quali2]==newcoord[k,3])]
				if(length(quali3)>0){
				loc[quali3]<-k
			}
		}
	}
	if(!is.null(distance)){
		pcoall<-cmdscale(distance)
		mat<-recluster.col(pcoall,st=F,rot=F)
	}
	if(is.null(xlim)){
		xlim<-range(long)
		}
	if(is.null(ylim)){
		ylim<-range(lat)
		}
	xylim<-cbind(xlim,ylim)
	if(!(add)){
		plot(cbind(xylim[1:2],xylim[3:4]),type="n",main=main,xlab=xlab,ylab=ylab)
	}
	if(!is.null(map)){
		plot(map, asp = 2,add=T,xlab=xlab,ylab=ylab)
	}
	if(is.null(minsize)){
		minsize<-min(abs(range(long)[1]-range(long)[2]),abs(range(lat)[1]-range(lat)[2]))/30
	}
	for (i in 1:max(loc)){
		specim<-which(loc==i)
		if(length(specim)==1){
			specimens<-c(long[specim], lat[specim],mat[specim,])
			plotrix::floating.pie(specimens[1],specimens[2],1,radius=minsize,border=NA,col=rgb(specimens[5], specimens[6], specimens[7], maxColorValue = 255))
			plotrix::draw.circle(specimens[1],specimens[2],radius=minsize)
		}
	if(length(specim)>1){
		specimens<-cbind(long[specim], lat[specim],mat[specim,])
		if(length(specim)>3){
			dista<-dist(mat[specim,1:2])
			if(sum(dista)>0){
				mds<-cbind(c(1:length(specim)),cmdscale(dist(mat[specim,1:2]),k=1))
				specimens<-specimens [order(mds[,2]),]
				}
			}
		if(proportional){
			rad<-minsize*(length(specim))^0.25
		}else{
			rad<-minsize
		}
		plotrix::floating.pie(mean(specimens[,1]),mean(specimens[,2]),rep(1,length(specim)),border=NA,radius=rad,col=rgb(specimens[, 5], specimens[, 6], specimens[, 7], maxColorValue = 255))
		plotrix::draw.circle(mean(specimens[,1]),mean(specimens[,2]),radius=rad)
		}
	}
}


