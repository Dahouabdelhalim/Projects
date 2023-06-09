#setwd

#load package
library(ggplot2)
library(ggrepel)
library(scales)
#load data
AP<-read.csv("Anterior-Posterior_3.13.21VolcanoPlot.csv", header=TRUE)
is.numeric(AP$pvalue)
APclean<-AP[with(AP, ave(pvalue, GeneID, FUN = min)==pvalue),]

#initial graph
ggplot(data=APclean, aes(x=Log2FoldChange, y=pvalue)) + geom_point()

p<-ggplot(data=APclean, aes(x=Log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()
p2<-p + geom_vline(xintercept = c(-0.6, 0.6), col="red") + geom_hline(yintercept = -log10(0.05), col="red")    

APclean$diffexpressed<-"Not Significant"
APclean$diffexpressed[APclean$Log2FoldChange > 0.6 & APclean$pvalue < 0.05] <- "Posterior"
APclean$diffexpressed[APclean$Log2FoldChange < -0.6 & APclean$pvalue < 0.05] <- "Anterior"

p<-ggplot(data=APclean, aes(x=Log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
p2<-p + geom_vline(xintercept = c(-0.6, 0.6), col="red") + geom_hline(yintercept = -log10(0.05), col="red")    
p2

p3<-p2 +scale_colour_manual(values=c("blue", "black", "magenta"))
p3
mycolors<-c("#00BFC4","#F8766D", "black")
mycolors2<-c("darkslateblue", "gold2", "black")
names(mycolors2)<-c("Anterior", "Posterior", "Not Significant")
p3<-p2 + scale_colour_manual(values=mycolors2)
p3

APclean$delabel<-NA
APclean$delabel[APclean$diffexpressed != "Not Significant"] <- APclean$GeneID[APclean$diffexpressed != "Not Significant"]

ggplot(data=APclean, aes(x=Log2FoldChange, y=-log10(p-value), col=diffexpressed, label=delabel))+ geom_point() + theme_minimal() + geom_text()

#install.packages("ggrepel")



tiff("Volcano Plot Anterior-Posterior.tiff", units="in", width=10, height=7, res=200)
ggplot(data=APclean, aes(x=Log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel))+ 
  geom_point() + 
  theme_classic() + 
  labs(y=expression('-log'[10]*'(p-value)'), x=expression('Log'[2]*'Fold Change'), col="Differentially Expressed") + 
  geom_text_repel(max.overlaps = 30) +
  scale_colour_manual(values=c("navy", "grey60", "lightskyblue")) 
dev.off()
  




 # geom_vline(xintercept = c(-0.6, 0.6), col="red") +
  #geom_hline(yintercept = -log10(0.05), col="red")












