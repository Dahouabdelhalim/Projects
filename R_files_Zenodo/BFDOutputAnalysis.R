require(RColorBrewer)
require(plotrix)

### Read in Empirical Data ###
empirical_out<-list.files("~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/SNAPP/Empirical/Output_v2",pattern=".out",full.names=T)

empirical_mL<-vector()
empirical_bf<-vector()
for(i in 1:length(empirical_out)){
	foo<-readLines(file(empirical_out[i]))
	bar<-strsplit(grep("marginal L",foo,value=T)," ")[[1]]
	empirical_mL[i]<-as.numeric(bar[length(bar)])
}
ngss<-length(empirical_mL)/2 #ngs = Number of Geographic Sampling Schemes

for(i in 1: ngss){
	foo<-empirical_mL[(2*((1:ngss)-1)+1)[i]:(2*((1:ngss)-1)+2)[i]]
	empirical_bf[i]<-2*(foo[1]-foo[2])
}
empirical_mL
empirical_bf

### Analyze SNAPP simulated output for BFD and generate figure ###
m_vec<-c(0.5,1,5,10,50,100)

### Read in output files ###
output_list<-list()
for(i in 1:length(m_vec)){
output_list[[i]]<-list.files("~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/SNAPP/Simulated/Output",pattern="SpDelim[12].out",full.names=T)
output_list[[i]]<-output_list[[i]][grep(paste("_",m_vec[i],"_",sep=""),output_list[[i]])]
}
output_list

mL_list<-list()
bf_list<-list()

for(j in 1:length(output_list)){	
	mL_list[[j]]<-vector()
	bf_list[[j]]<-vector()
	for(i in 1:length(output_list[[j]])){
		foo<-readLines(file(output_list[[j]][i]))
		bar<-strsplit(grep("marginal L",foo,value=T)," ")[[1]]
		mL_list[[j]][i]<-as.numeric(bar[length(bar)])
	}
	ngss<-length(mL_list[[j]])/2 #ngs = Number of Geographic Sampling Schemes

	for(i in 1: ngss){
		foo<-mL_list[[j]][(2*((1:ngss)-1)+1)[i]:(2*((1:ngss)-1)+2)[i]]
		bf_list[[j]][i]<-2*(foo[1]-foo[2])
	}
	
}

mL_list
bf_list

### Make plot with all levels of m investigated and included here ###
png(file="~/Desktop/Manuscripts/SpeciesDelimitationSensitivity/Figures/BFD_results_v5.png",width=3.25,height=8,units="in",res=500)
#quartz(width=3.25,height=8)

layout(matrix(c(1,2,3),nrow=3),heights=c(0.35,0.35,0.3))

### Simulated data scatterplot ###
par(mar=c(1,4,1,1))
plot(0,0,type="n",ylim=c(-810,810),xlim=c(1,4),xlab=" ",ylab=" ",axes=F)
axis(1,at=1:4,labels=F,cex=0.8, tck=-0.025)
axis(2,at=seq(-800,800,100),cex.axis=0.8,mgp=c(3,0.5,0),tck=-0.025)
mtext("Bayes Factor",2,line=2.5,cex=1)
par(xpd=NA)
text("One Species",srt=90,cex=1,x=0.55,y=500,col="forest green")
text("Two Species",srt=90,cex=1,x=0.55,y=-500,col="orange")

box()

par(xpd=F)
abline(h=0,lty=2)

for(i in 1:length(bf_list)){
	points(bf_list[[i]],type="l",col=brewer.pal(9,"Blues")[-c(1,2,3)][i],lwd=0.75)
	points(bf_list[[i]],bg=paste(brewer.pal(9,"Blues")[-c(1,2,3)][i],"90",sep=""),col="black",cex=2,lwd=0.5,pch=c(21:25,21:25,21:25)[i])
}
legend("topright",legend=rev(paste("m = ",m_vec,sep="")),pch=rev(c(21:25,21:25,21:25)[1:length(bf_list)]),pt.bg=paste(rev(brewer.pal(9,"Blues")[-c(1,2,3)][1:6]),"90",sep=""),pt.cex=1.5,pt.lwd=0.5,cex=0.8,y.intersp=1.2)

par(xpd=NA)
text(x=par()$usr[1]*0.45,y=par()$usr[4]*1.0,label="A",cex=2,font=2)
text(x=par()$usr[1],y=par()$usr[4]-(par()$usr[4]-par()$usr[3])*0.025,label=" Simulated",cex=1,font=1,adj=c(0,0.5))


### Empirical Scatterplot ###
par(mar=c(1,4,1,1))
plot(0,0,type="n",ylim=c(-2500,2500),xlim=c(1,4),xlab=" ",ylab=" ",axes=F)
axis(1,at=1:4,labels=F,cex=0.8, tck=-0.025)
axis(2,at=seq(-2500,2500,250),cex.axis=0.5,mgp=c(3,0.5,0),tck=-0.025,labels=F)
axis(2,at=seq(-2500,2500,500),cex.axis=0.6,mgp=c(3,0.5,0),tck=-0.025,labels=seq(-2500,2500,500))

mtext("Bayes Factor",2,line=2.5,cex=1)
par(xpd=NA)
text("One Species",srt=90,cex=1,x=0.55,y=1250,col="forest green")
text("Two Species",srt=90,cex=1,x=0.55,y=-1250,col="orange")

box()

par(xpd=F)
abline(h=0,lty=2)

points(empirical_bf, type="l",lwd=0.75)
points(empirical_bf, pch=21,lwd=0.75,cex=2)

par(xpd=NA)
text(x=par()$usr[1]*0.45,y=par()$usr[4]*1.0,label="B",cex=2,font=2)
text(x=par()$usr[1],y=par()$usr[4]-(par()$usr[4]-par()$usr[3])*0.025,label=" Empirical",cex=1,font=1,adj=c(0,0.5))

### Geographic sampling scenario panels ###
plot(0,0,type="n",ylim=c(-10,10),xlim=c(1,4),xlab=" ",ylab=" ",axes=F)

col1<-colorRampPalette(colors=c("red","blue"))(9)
col4<-col3<-col2<-col1
col2[5]<-"white"
col3[4:6]<-"white"
col4[3:7]<-"white"

points(rep(1,9),seq(10,-5,length.out=9),cex=3,pch=21,bg=col1, col="black")
points(rep(2,9), seq(10,-5,length.out=9),cex=3,pch=21,bg=col2,col="black")
points(rep(3,9), seq(10,-5,length.out=9),cex=3,pch=21,bg=col3, col="black")
points(rep(4,9), seq(10,-5,length.out=9),cex=3,pch=21,bg=col4, col="black")


par(xpd=NA)
text(x=par()$usr[1]*0.45,y=par()$usr[4]*1.0,label="C",cex=2,font=2)

arrows(x0=0.8,y0=-9,x1=4.2,y1=-9,code=3,angle=15)
text(x=2.5,y= -8,label="Admixed Populations",font=1,cex=1)
text(x=1,y= -11,label="More",font=1,cex=0.8)
text(x=4,y= -11,label="Fewer",font=1,cex=0.8)

dev.off()
