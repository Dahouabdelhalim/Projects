library(ggplot2)
library(ggthemes)

ResponsiveGenes<-data.frame(t(read.table("Root_DEGenes_NormalizedCounts.txt",header=T,sep="\\t",row.names=1)))
ResponsiveGenes + 1
genes<-colnames(ResponsiveGenes)

genesLocus<-read.table(file=paste("Genenames", ".txt", sep=""))
genesLocusup<-toupper(genesLocus[,1])
genesDescrip<-readLines("GeneDescrip.txt")

#Shoot
#ResponsiveGenes$Condition <- c("noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit")
#ResponsiveGenes$TimePoint <-c(0,0,0,5,5,5,10,10,10,15,15,15,20,20,20,30,30,30,45,45,45,60,60,60,90,90,120,120,120,0,0,0,5,5,5,10,10,10,15,15,15,20,20,20,30,30,30,45,45,45,60,60,60,90,90,90,120,120,120)
#Root
ResponsiveGenes$Condition <- c("noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNot","noNit","noNit","noNit","noNit","noNit","noNit","noNit","noNit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit","Nit")
ResponsiveGenes$TimePoint <-c(0,0,0,5,5,5,10,10,10,15,15,15,20,20,20,30,30,30,45,45,45,60,60,60,90,90,90,120,120,120,0,0,0,5,5,5,10,10,10,15,15,15,20,20,20,30,30,30,45,45,45,60,60,60,90,90,90,120,120,120)

pdf("Root_N_Responsive.pdf")

for( gene in genes ){ 
  cat(gene)
  if (gene %in% genesLocusup){
  indx_gene<-which(genesLocusup%in%toupper(gene))
  genedesc<-genesDescrip[indx_gene]
  } else{
	genedesc<-"Uknown protein"
  }
	d<-ggplot(data=ResponsiveGenes,aes(x=TimePoint,y=ResponsiveGenes[,c(gene)],colour=Condition))
	d <- d + labs(x="Time in Minutes",y="log10( Normalized Read Count )") 
	d <- d + ggtitle(gene,genedesc)
	d <- d + geom_point()
	d <- d + scale_x_continuous(breaks=c(0,5,10,15,20,30,45,60,90,120)) 
	
	maxval <- max(ResponsiveGenes[[gene]])
	minval <- min(ResponsiveGenes[[gene]])
	if (maxval > 10000) {
	   newmax <- round(maxval,digits = -4)
	   d <- d + scale_y_continuous(trans="log2",limits=c(1,newmax),breaks=c(100,1000,10000))
	} else {
	   d <- d + scale_y_continuous(trans="log2",limits=c(1,10000),breaks=c(100,1000,10000))
	}
	d <- d + stat_smooth(se=FALSE,method="loess",span=0.5)
	d <- d + theme(plot.title = element_text(size=18, face="bold", hjust=0.5,vjust=1, lineheight=0.6), axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),legend.position = "top", legend.direction = "horizontal", legend.box = "vertical",legend.title=element_blank(),legend.text=element_text(size=18))
	d <- d + theme(plot.subtitle = element_text(size=12, hjust=0.5,face="italic"))
	print(d)
}
dev.off()
