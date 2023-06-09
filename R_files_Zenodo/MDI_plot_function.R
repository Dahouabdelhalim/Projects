#DTT accounting for phylogenetic uncertainty



mdi_phy<-function(tree,x){
 		disp<-dtt(tree,x,index="avg.sq",nsim=1000,plot=FALSE)
 		mdi<-disp$MDI
 		val<-disp$dtt
 		time<-disp$times
 		sim<-t(disp$sim)
 		return(list(MDI=mdi,time=time,DTT=val,SIM=sim))
 }

plots2<-function(x,y) {
	x.x<-apply(x,1,function (x) data.frame(x))
	y.y<-apply(y,1,function (y) data.frame(y))
	comb<-mapply(data.frame, y.y, x.x, SIMPLIFY=FALSE)
	lapply(comb,function (x) lines(as.data.frame(x),col=rgb(1,0,0,0.1)))
	rm(x.x)
	rm(y.y)
	rm(comb)
}

elems<-c("MDI","time","DTT","SIM") 

mdi<-lapply(posterior_trees,mdi_phy,trait_data)

res<-sapply(elems,function(E){do.call(rbind,lapply(mdi,function(X){X[[E]]}))},simplify=FALSE)
 
plot(NA,NA,xlim=c(0.0,1.0),ylim=c(0.0,1.25),xlab="Time (my bp)",ylab="Disparity (MDI)")
axis(1,at=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
axis(2,at=c(0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0))


s<-data.frame(res$SIM)
t<-data.frame(res$time)
v<-data.frame(res$DTT)
x.MDI<-mean(res$MDI)
x.MDI


sim.median<-apply(s,2,median)
time.median<-apply(t,2,median)
val.median<-apply(v,2,median)
sd<-apply(s,2,sd)
mean<-apply(s,2,mean)
upper<-mean+(2*sd)
lower<-mean-(2*sd)
s.u<-data.frame(time.median,upper)
s.l<-data.frame(rev(time.median),rev(lower))
colnames(s.l)<-colnames(s.u)
poly<-rbind(s.u,s.l)
polygon(poly,col=grey(0.1,0.3),border=NA)

plots2(v,t)

obs<-dtt(median_tree,trait_data,index="avg.sq",nsim=1000,plot=FALSE)
obs.dtt<-obs$dtt
obs.time<-obs$times
lines(time.median,sim.median,col=rgb(1,1,0,1),lty=2,lwd=4)
lines(obs.time,obs.dtt,col=rgb(0,0,0,1),lwd=4)
lines(time.median,val.median,col=rgb(0,0.9,1,1),lty=2,lwd=4)


dev.copy2pdf(file="trait_dtt.pdf")
