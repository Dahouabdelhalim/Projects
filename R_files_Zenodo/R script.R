Figure 2

rfc<-read.csv(file.choose())
str(rfc)
rfc.melt<-melt(rfc, id.vars=c("numsampled"), na.rm=TRUE)
colnames(rfc.melt)[2] = "colonies"
colnames(rfc.melt)[3] = "otus"

str(rfc.melt)

theme_format <- theme_bw()+
  theme(axis.text.x  = element_text(size=10, colour = "black"))+
  theme(axis.ticks = element_line(colour="black"))+
  theme(panel.grid.minor=element_blank())+
  theme(panel.grid.major=element_blank())

ggplot(rfc.melt,aes(numsampled, otus))+   
 #geom_line(data = rfc.melt, aes(numsampled, otus, fill=colonies))+
 geom_smooth(data = rfc.melt, aes(numsampled, otus, colour=colonies), se=FALSE, size=0.5)+
  #guides(colour = guide_legend("Region"), fill="white")+
  guides(colour= FALSE)+
  theme_format+
  coord_cartesian(ylim = c(0, 220))+
  scale_y_continuous(breaks=c(0,10,20,30,40),name= "Estimated number of unique OTUs")+
  coord_cartesian(xlim = c(0, 50000))+
  scale_x_continuous(breaks=c(0, 20000, 40000),name = "Sequences sampled")



Figure 3

colonies<-read.csv(file.choose())

str(colonies)

colonies.matrix<-data.matrix(colonies) #for data matrix

colonies.matrix.prop<-prop.table(colonies.matrix,2)
colonies.matrix.prop.sqrt<-sqrt(colonies.matrix.prop)

colonies.matrix.sqrt[colonies.matrix.sqrt==0]

heatmap.2(scott.matrix.sqrt, Rowv=NA, na.col="white",dendrogram="none",Colv=NA,key=TRUE, keysize=1.2, main=NULL, sepwidth=c(0.001,0.001),sepcolor="dark grey",colsep=0:14,rowsep=0:34, col=colors,trace="none", tracecol="black",margins=c(10,10),symm=F,symkey=F,scale="none",xlab="",ylab="", cexCol=1.0)

