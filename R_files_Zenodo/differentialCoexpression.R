# Import necessary libraries
library(WGCNA)
library(preprocessCore)
library(RColorBrewer)

allowWGCNAThreads()

# Import necessary datasets
load("/nethome/sxp720/Raw_Data/exprNoneTSubset.RData")
load("/nethome/sxp720/Raw_Data/exprBothTSubset.RData")
load("/nethome/sxp720/Raw_Data/exprMTSubset.RData")
load("/nethome/sxp720/Raw_Data/exprRTSubset.RData")


# Supporting functions
plotHeatmap<-function(colorh1C1C2,AdjMat1C1,AdjMat1C2, datC1, datC2,ordering=NULL,file="DifferentialPlot.png")
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(datC1[,union(union(which(colorh1C1C2!="grey"), which(colorh1C1C2!="blue")), which(colorh1C1C2!="royalblue"))],colorh1C1C2[union(union(which(colorh1C1C2!="grey"), which(colorh1C1C2!="blue")), which(colorh1C1C2!="royalblue"))],rbind(datC1,datC2)[,union(union(which(colorh1C1C2!="grey"), which(colorh1C1C2!="blue")), which(colorh1C1C2!="royalblue"))])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1C1C2 ==c))
    }
  }
  mat_tmp<-(AdjMat1C1[ordering,ordering])
  mat_tmp[which(row(mat_tmp)>col(mat_tmp))]<-(AdjMat1C2[ordering,ordering][which(row(mat_tmp)>col(mat_tmp))])
  diag(mat_tmp)<-0
  mat_tmp<-sign(mat_tmp)*abs(mat_tmp)^(1/2)
  png(file=file,height=1000,width=1000)
  image(mat_tmp,col=rev(brewer.pal(11,"RdYlBu")),axes=F,asp=1,breaks=seq(-1,1,length.out=12))
  dev.off()
  unique(colorh1C1C2[ordering])
}

plotExprChange<-function(datC1,datC2, colorh1C1C2,ordering=NULL)
{
  if (is.null(ordering))
  {
    h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(datC1[,which(colorh1C1C2!="grey")],colorh1C1C2[which(colorh1C1C2!="grey")],rbind(datC1,datC2)[,which(colorh1C1C2!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1C1C2 ==c))
    }
  }
  mycolors<-colorh1C1C2[ordering]
  plot(x=0:length(which(mycolors!="grey")),y=rep(1,length(which(mycolors!="grey"))+1),col="white",axes=F,xlab="",ylab="",ylim=c(0,1))
  rr=c(244,239,225,215,209,193,181,166,151,130,110)
  gg=c(228,204,174,160,146,117,94,58,44,45,45)
  bb=c(176,140,109,105,102,91,84,74,70,68,66)
  MyColours<-NULL
  for ( i in 1:11)
  {
    MyColours=c(MyColours,rgb(rr[i],gg[i],bb[i],maxColorValue=255)  )
  }
  exprDiff<-NULL
  l<-0
  for (c in setdiff(unique(mycolors),"grey"))
  {
    meanC1<-mean(t(datC1)[colnames(datC1)[which(colorh1C1C2 == c)],])
    meanC2<-mean(t(datC2)[colnames(datC2)[which(colorh1C1C2 == c)],])
    exprDiff<-rbind(exprDiff,c(meanC1,meanC2))
    r<-l+length(which(mycolors==c))
    rect(l,0.85,r,1,col=c,border=F)
    rect(l,0,r,.4,col=MyColours[floor(meanC2*2)-10],border="white",lwd=2)
    rect(l,0.4,r,.8,col=MyColours[floor(meanC1*2)-10],border="white",lwd=2)
    l<-r
  }
  exprDiff
}

getEigenGeneValues<-function(datRef,colorh1,datAll)
{
  eigenGenesCoef<-list()
  i<-0
  for (c in unique(colorh1))
  {
    i<-i+1
    eigenGenesCoef[[i]]<-prcomp(scale(datRef[,which(colorh1 == c)]))$rotation[,1]
  }
  names(eigenGenesCoef)<-unique(colorh1)
  values<-NULL
  for( c in unique(colorh1))
  {
    v<-rbind(datAll)[,which(colorh1 == c)] %*%  eigenGenesCoef[[c]]
    values<-cbind(values,sign(mean(v))*v)
  }
  colnames(values)<-unique(colorh1)
  values
}

dispersionModule2Module<-function(c1,c2,datC1,datC2,colorh1C1C2)
{
    if (c1==c2)
    {
       difCor<-(cor(datC1[,which(colorh1C1C2 == c1)],method="spearman")-
       cor(datC2[,which(colorh1C1C2 == c1)],method="spearman"))^2
       n<-length(which(colorh1C1C2  ==c1))
      (1/((n^2 -n)/2)*(sum(difCor)/2))^(.5)
    }
    else if (c1!=c2)
    {
      difCor<-(cor(datC1[,which(colorh1C1C2 == c1)],datC1[,which(colorh1C1C2==c2)],method="spearman")-
              cor(datC2[,which(colorh1C1C2 == c1)],datC2[,which(colorh1C1C2==c2)],method="spearman"))^2
     n1<-length(which(colorh1C1C2  ==c1))
     n2<-length(which(colorh1C1C2  ==c2))
     (1/((n1*n2))*(sum(difCor)))^(.5)
    }
}

permutationProcedureModule2Module<-function(permutation,d,c1,c2,colorh1C1C2)
{
  d1<-d[permutation,]
  d2<-d[-permutation,]
  dispersionModule2Module(c1,c2,d1,d2,colorh1C1C2)
}


beta1=4 #user defined parameter for soft thresholding

# Assign datasets to more general names
datC1<-exprBothTSubset
datC2<-exprMTSubset
datC3<-exprRTSubset
datC4<-exprNoneTSubset

# Calculate adjacency to be used in later calculations
AdjMatC1<-sign(cor(datC1,method="spearman"))*(cor(datC1,method="spearman"))^2
AdjMatC2<-sign(cor(datC2,method="spearman"))*(cor(datC2,method="spearman"))^2
AdjMatC3<-sign(cor(datC3,method="spearman"))*(cor(datC3,method="spearman"))^2
AdjMatC4<-sign(cor(datC4,method="spearman"))*(cor(datC4,method="spearman"))^2

diag(AdjMatC1)<-0
diag(AdjMatC2)<-0
diag(AdjMatC3)<-0
diag(AdjMatC4)<-0

# Calculate average of the matrices for use in distance matrix calculations
AdjMatC0<-1/4*(AdjMatC1+AdjMatC2+AdjMatC3+AdjMatC4)

# Calculate distance matrix
D <- (1/6 * (abs(AdjMatC1 - AdjMatC0) + abs(AdjMatC2 - AdjMatC0) + abs(AdjMatC3 - AdjMatC0) + abs(AdjMatC4 - AdjMatC0)))^(beta1/2)

# Convert distance matrix to topoligical overlap matrix
dissTOMC1C2C3C4 <- TOMdist(D)

# Save the distance matrix as it is very computationally expensive to generate
save(dissTOMC1C2C3C4, file = "dissTOMC1C2C3C4beta4.RData")

# Clear garbage
collectGarbage()

# load("/nethome/sxp720/Scripts/dissTOMC1C2C3C4beta4.RData")
#Hierarchical clustering is performed using the Topological Overlap of the adjacency difference as input distance matrix
geneTreeC1C2C3C4= hclust(as.dist(dissTOMC1C2C3C4), method = "average");

png(file="hierarchicalTree.png",height=2000,width=2000)
plot(geneTreeC1C2C3C4, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
dev.off()

dynamicModsHybridC1C2C3C4 <- cutreeDynamic(dendro = geneTreeC1C2C3C4, distM = dissTOMC1C2C3C4, method="hybrid", deepSplit = T, pamRespectsDendro = FALSE, minClusterSize = 100, cutHeight = .980)

#Every module is assigned a color. Note that GREY is reserved for genes which do not belong to any differentially coexpressed module
dynamicColorsHybridC1C2C3C4 = labels2colors(dynamicModsHybridC1C2C3C4)

# Merge close modules according to WGCNA's instructions
mergedColorC1C2C3C4<-mergeCloseModules(rbind(datC1,datC2,datC3,datC4),dynamicColorsHybridC1C2C3C4,cutHeight=.2)$color
colorh1C1C2C3C4<-mergedColorC1C2C3C4

# Export the modules and the genes in the modules
modules <- data.frame(colnames(datC1), colorh1C1C2C3C4)
colnames(modules) <- c("GeneName", "Module")
write.table(modules, file = "ModuleExport.txt")

ordering <- c()
h<-hclust(as.dist(1-abs(cor(getEigenGeneValues(datC1[,which(colorh1C1C2C3C4!="grey")],colorh1C1C2C3C4[which(colorh1C1C2C3C4!="grey")],rbind(datC1,datC2)[,which(colorh1C1C2C3C4!="grey")])))))
    for (c in h$label[h$order])
    {
      ordering<-c(ordering,which(colorh1C1C2C3C4 ==c))
    }

# Export plots
# plotHeatmap(colorh1C1C2C3C4,AdjMatC1,AdjMatC2, datC1, datC2, file="BothvsMHeatmap.png", ordering = ordering)
# plotHeatmap(colorh1C1C2C3C4,AdjMatC1,AdjMatC3, datC1, datC3, file="BothvsRHeatmap.png", ordering = ordering)
# plotHeatmap(colorh1C1C2C3C4,AdjMatC1,AdjMatC4, datC1, datC4, file="BothvsNoneHeatmap.png", ordering = ordering)

plotHeatmap(colorh1C1C2C3C4,AdjMatC1,AdjMatC1, datC1, datC1, file="BothHeatmap1.png", ordering = ordering)
plotHeatmap(colorh1C1C2C3C4,AdjMatC2,AdjMatC2, datC2, datC2, file="MHeatmap1.png", ordering = ordering)
plotHeatmap(colorh1C1C2C3C4,AdjMatC3,AdjMatC3, datC3, datC3, file="RHeatmap1.png", ordering = ordering)
plotHeatmap(colorh1C1C2C3C4,AdjMatC4,AdjMatC4, datC4, datC4, file="NoneHeatmap1.png", ordering = ordering)


# png(file="exprChangeBothvsM.png",height=2000,width=2000)
# plotExprChange(datC1,datC2,colorh1C1C2C3C4, ordering = ordering)
# dev.off()

# png(file="exprChangeBothvsR.png",height=2000,width=2000)
# plotExprChange(datC1,datC3,colorh1C1C2C3C4, ordering = ordering)
# dev.off()

# png(file="exprChangeBothvsNone.png",height=2000,width=2000)
# plotExprChange(datC1,datC4,colorh1C1C2C3C4, ordering = ordering)
# dev.off()

# #save.image("differentialExpressionbeta4.RData");

# #Functional Enrichment analysis 

# # Permutation Analysis

# permutations<-NULL

# colorh1C1C2 <- diffCoExModules$Module

# for (i in 1:1000)
# {
   # permutations<-rbind(permutations,sample(1:(nrow(datC1)+nrow(datC2)),nrow(datC1)))
# }
# d<-rbind(scale(datC1),scale(datC2))

# dispersionMatrix<-matrix(nrow=length(unique(colorh1C1C2))-1,ncol=length(unique(colorh1C1C2))-1)
# nullDistrib<-list()
# i<-j<-0
# for (c1 in setdiff(unique(colorh1C1C2),"grey"))
# {
  # i<-i+1
  # j<-0
  # nullDistrib[[c1]]<-list()
  # for (c2 in setdiff(unique(colorh1C1C2),"grey"))
  # {
    # j<-j+1
	# if(c1 == c2){
    # dispersionMatrix[i,j]<-dispersionModule2Module(c1,c2,datC1,datC2,colorh1C1C2)
    # nullDistrib[[c1]][[c2]]<-apply(permutations,1,permutationProcedureModule2Module,d,c2,c1,colorh1C1C2)
	# }
  # }
# }


# permutations<-NULL
# for (i in 1:1000)
# {
   # permutations<-rbind(permutations,sample(1:(nrow(datC1)+nrow(datC3)),nrow(datC1)))
# }
# d<-rbind(scale(datC1),scale(datC3))
# dispersionMatrix2<-matrix(nrow=length(unique(colorh1C1C2))-1,ncol=length(unique(colorh1C1C2))-1)
# nullDistrib2<-list()
# i<-j<-0

# for (c1 in setdiff(unique(colorh1C1C2),"grey"))
# {
  # i<-i+1
  # j<-0
  # nullDistrib2[[c1]]<-list()
  # for (c2 in setdiff(unique(colorh1C1C2),"grey"))
  # {
    # j<-j+1
    # if(c1 == c2){
	# dispersionMatrix2[i,j]<-dispersionModule2Module(c1,c2,datC1,datC3,colorh1C1C2)
    # nullDistrib2[[c1]][[c2]]<-apply(permutations,1,permutationProcedureModule2Module,d,c2,c1,colorh1C1C2)
	# }
  # }
# }

# permutations<-NULL
# for (i in 1:1000)
# {
   # permutations<-rbind(permutations,sample(1:(nrow(datC1)+nrow(datC4)),nrow(datC1)))
# }
# d<-rbind(scale(datC1),scale(datC4))
# dispersionMatrix3<-matrix(nrow=length(unique(colorh1C1C2))-1,ncol=length(unique(colorh1C1C2))-1)
# nullDistrib3<-list()
# i<-j<-0

# for (c1 in setdiff(unique(colorh1C1C2),"grey"))
# {
  # i<-i+1
  # j<-0
  # nullDistrib3[[c1]]<-list()
  # for (c2 in setdiff(unique(colorh1C1C2),"grey"))
  # {
    # j<-j+1
    # if(c1 == c2){
	# dispersionMatrix3[i,j]<-dispersionModule2Module(c1,c2,datC1,datC4,colorh1C1C2)
    # nullDistrib3[[c1]][[c2]]<-apply(permutations,1,permutationProcedureModule2Module,d,c2,c1,colorh1C1C2)
	# }
  # }
# }

# permutations<-NULL
# for (i in 1:1000)
# {
   # permutations<-rbind(permutations,sample(1:(nrow(datC2)+nrow(datC3)),nrow(datC2)))
# }
# d<-rbind(scale(datC2),scale(datC3))
# dispersionMatrix4<-matrix(nrow=length(unique(colorh1C1C2))-1,ncol=length(unique(colorh1C1C2))-1)
# nullDistrib4<-list()
# i<-j<-0

# for (c1 in setdiff(unique(colorh1C1C2),"grey"))
# {
  # i<-i+1
  # j<-0
  # nullDistrib4[[c1]]<-list()
  # for (c2 in setdiff(unique(colorh1C1C2),"grey"))
  # {
    # j<-j+1
    # if(c1 == c2){
	# dispersionMatrix4[i,j]<-dispersionModule2Module(c1,c2,datC2,datC3,colorh1C1C2)
    # nullDistrib4[[c1]][[c2]]<-apply(permutations,1,permutationProcedureModule2Module,d,c2,c1,colorh1C1C2)
	# }
  # }
# }

# permutations<-NULL
# for (i in 1:1000)
# {
   # permutations<-rbind(permutations,sample(1:(nrow(datC2)+nrow(datC4)),nrow(datC2)))
# }
# d<-rbind(scale(datC2),scale(datC4))
# dispersionMatrix5<-matrix(nrow=length(unique(colorh1C1C2))-1,ncol=length(unique(colorh1C1C2))-1)
# nullDistrib5<-list()
# i<-j<-0

# for (c1 in setdiff(unique(colorh1C1C2),"grey"))
# {
  # i<-i+1
  # j<-0
  # nullDistrib5[[c1]]<-list()
  # for (c2 in setdiff(unique(colorh1C1C2),"grey"))
  # {
    # j<-j+1
    # if(c1 == c2){
	# dispersionMatrix5[i,j]<-dispersionModule2Module(c1,c2,datC2,datC4,colorh1C1C2)
    # nullDistrib5[[c1]][[c2]]<-apply(permutations,1,permutationProcedureModule2Module,d,c2,c1,colorh1C1C2)
	# }
  # }
# }



# permutations<-NULL
# for (i in 1:1000)
# {
   # permutations<-rbind(permutations,sample(1:(nrow(datC3)+nrow(datC4)),nrow(datC2)))
# }
# d<-rbind(scale(datC3),scale(datC4))
# dispersionMatrix6<-matrix(nrow=length(unique(colorh1C1C2))-1,ncol=length(unique(colorh1C1C2))-1)
# nullDistrib6<-list()
# i<-j<-0

# for (c1 in setdiff(unique(colorh1C1C2),"grey"))
# {
  # i<-i+1
  # j<-0
  # nullDistrib6[[c1]]<-list()
  # for (c2 in setdiff(unique(colorh1C1C2),"grey"))
  # {
    # j<-j+1
    # if(c1 == c2){
	# dispersionMatrix6[i,j]<-dispersionModule2Module(c1,c2,datC3,datC4,colorh1C1C2)
    # nullDistrib6[[c1]][[c2]]<-apply(permutations,1,permutationProcedureModule2Module,d,c2,c1,colorh1C1C2)
	# }
  # }
# }

