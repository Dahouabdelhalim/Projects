setwd("~/EMBER/Reintegrated_Data_Nov2020")
dd_clinical<-read.table(file="FINAL_DATASET.csv",header=TRUE,sep=",")
for(i in 1:nrow(dd_clinical)){
if(dd_clinical$Diagnosis[i]==1) {dd_clinical$Grp[i]<-"Healthy"}
else if(dd_clinical$Diagnosis[i]==2) {dd_clinical$Grp[i]<-"Acute Asthma"}
else if(dd_clinical$Diagnosis[i]==3) {dd_clinical$Grp[i]<-"Acute COPD"}
else if(dd_clinical$Diagnosis[i]==4) {dd_clinical$Grp[i]<-"HF"}
else if(dd_clinical$Diagnosis[i]==5) {dd_clinical$Grp[i]<-"Pneumonia"}
}
for(i in 1:nrow(dd_clinical)){	
if(is.na(dd_clinical$V1.BNP[i])) {dd_clinical$V1.BNP[i]<-median(dd_clinical$V1.BNP,na.rm=TRUE)}	
if(!is.na(dd_clinical$V1.BNP[i])) {dd_clinical$V1.BNP[i]<-dd_clinical$V1.BNP[i]}
if(is.na(dd_clinical$V1.CRP[i])) {dd_clinical$V1.CRP[i]<-median(dd_clinical$V1.CRP,na.rm=TRUE)}	
if(!is.na(dd_clinical$V1.CRP[i])) {dd_clinical$V1.CRP[i]<-dd_clinical$V1.CRP[i]}	
if(is.na(dd_clinical$V1.Eos.count[i])) {dd_clinical$V1.Eos.count[i]<-median(dd_clinical$V1.Eos.count,na.rm=TRUE)}	
if(!is.na(dd_clinical$V1.Eos.count[i])) {dd_clinical$V1.Eos.count[i]<-dd_clinical$V1.Eos.count[i]}	
}
	
for(i in 1:nrow(dd_clinical)){	
if(is.na(dd_clinical$V1.BNP[i])) {dd_clinical$BNP_Grp[i]<-NA}	
if(dd_clinical$V1.BNP[i]>=500 & !is.na(dd_clinical$V1.BNP[i])) {dd_clinical$BNP_Grp[i]<-"High"}
if(dd_clinical$V1.BNP[i]<500 & !is.na(dd_clinical$V1.BNP[i])) {dd_clinical$BNP_Grp[i]<-"Low"}	
if(is.na(dd_clinical$V1.CRP[i])) {dd_clinical$CRP_Grp[i]<-NA}	
if(dd_clinical$V1.CRP[i]>=50 & !is.na(dd_clinical$V1.CRP[i])) {dd_clinical$CRP_Grp[i]<-"High"}
if(dd_clinical$V1.CRP[i]<50 & !is.na(dd_clinical$V1.CRP[i])) {dd_clinical$CRP_Grp[i]<-"Low"}
if(is.na(dd_clinical$V1.Eos.count[i])) {dd_clinical$EOS_Grp[i]<-NA}	
if(dd_clinical$V1.Eos.count[i]>=0.3 & !is.na(dd_clinical$V1.Eos.count[i])) {dd_clinical$EOS_Grp[i]<-"High"}
if(dd_clinical$V1.Eos.count[i]<0.3 & !is.na(dd_clinical$V1.Eos.count[i])) {dd_clinical$EOS_Grp[i]<-"Low"}		
}

dd_clinical_grps<-data.frame(dd_clinical$Sample_ID,dd_clinical$BNP_Grp,dd_clinical$CRP_Grp,dd_clinical$EOS_Grp,dd_clinical$Grp)
names(dd_clinical_grps)<-c("Sample_ID","BNP_Grp","CRP_Grp","EOS_Grp","Grp")
	

setwd("~/EMBER/Graph_Experiments")
dd_disc<-read.table(file="Discovery_TDA_Features_Rows.csv",header=TRUE,sep=",")
dd_rep<-read.table(file="Replication_TDA_Features_Rows.csv",header=TRUE,sep=",")
dd_tmp<-rbind(data.frame(t(dd_disc)),data.frame(t(dd_rep)))
dd<-t(dd_tmp)
dd<-data.frame(dd)
dd<-dd[,-grep("Healthy",names(dd))]


library(igraph)
library(lsa)
library(WGCNA)
library(dunn.test)
library(netcom)
library(bluster)


#All data Graph
dd_ALL<-dd
dd_c<-t(dd_ALL)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_ALL$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-apply(dd_c,2,as.numeric)
d<-abs(cor(dd_c))
d_cut<-0.5
d[d>d_cut]<-0
ALL_g<-graph.adjacency(d,mode="undirected",weighted=TRUE,diag=FALSE)
plot(ALL_g, edge.label=round(E(ALL_g)$weight, 3))
ALL_g_clust<-cluster_louvain(ALL_g,weights = edge_attr(ALL_g)$weight)
plot(ALL_g_clust,ALL_g,mark.groups=NULL,main="ALL")
communities(ALL_g_clust)

#Unsigned
dd_ALL<-dd
dd_c<-t(dd_ALL)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_ALL$Feature
dd_c<-dd_c[-113,]
dd_c<-dd_c[-nrow(dd_c),]	
dd_c<-apply(dd_c,2,as.numeric)
pickSoftThreshold(dd_c,dataIsExpr = TRUE)
#Beta = 5 for Disc
beta_pow<-5
#Beta = 7 for Rep
#beta_pow<-7
d<-abs(cor(dd_c))^beta_pow
k=as.vector(apply(d,2,sum, na.rm=T))
scaleFreePlot(k, main="Check scale free topology\\n")
#d<-abs((1+cor(dd_c))/2)^3
#TOM_sim<-TOMsimilarity(d)
#colnames(TOM_sim)<-colnames(d)
#rownames(TOM_sim)<-rownames(d)	
ALL_g<-graph.adjacency(d,mode="undirected",weighted=TRUE,diag=FALSE)
#ALL_g<-graph.adjacency(TOM_sim,mode="undirected",weighted=TRUE,diag=FALSE)
plot(ALL_g, edge.label=round(E(ALL_g)$weight, 3))


ALL_g_clust<-cluster_louvain(ALL_g,weights = edge_attr(ALL_g)$weight)
plot(ALL_g_clust,ALL_g,mark.groups=NULL,main="ALL")

tiff(filename = "Disc_Rep_Graph_Louvain.tiff" ,units="in", width=13, height=9, res=400,compression = "lzw")	
plot(ALL_g_clust,ALL_g,vertex.label="",mark.groups=NULL,main="Discovery and Replication, 101 Features, Graph with Louvain Clusters")
dev.off()
communities(ALL_g_clust)

	
	
	
#Signed
dd_ALL<-dd
dd_c<-t(dd_ALL)
dd_c<-data.frame(dd_c)
names(dd_c)<-dd_ALL$Feature
dd_c<-dd_c[-nrow(dd_c),]
dd_c<-apply(dd_c,2,as.numeric)
pickSoftThreshold(dd_c,dataIsExpr = TRUE, networkType = "signed")
#d<-abs(cor(dd_c))^3
d<-abs((1+cor(dd_c))/2)^3
#TOM_sim<-TOMsimilarity(d)
#colnames(TOM_sim)<-colnames(d)
#rownames(TOM_sim)<-rownames(d)	
ALL_g<-graph.adjacency(d,mode="undirected",weighted=TRUE,diag=FALSE)
#ALL_g<-graph.adjacency(TOM_sim,mode="undirected",weighted=TRUE,diag=FALSE)
plot(ALL_g, edge.label=round(E(ALL_g)$weight, 3))

ALL_g_clust<-cluster_louvain(ALL_g,weights = edge_attr(ALL_g)$weight)
plot(ALL_g_clust,ALL_g,mark.groups=NULL,main="ALL")
communities(ALL_g_clust)





library(rlist)
featureSets<-list()
for(i in 1:length(communities(ALL_g_clust))){
featureSets<-list.append(featureSets,communities(ALL_g_clust)[[i]])	
}	
#names(featureSets)<-c("Set1","Set2","Set3","Set4","Set5")
names(featureSets)<-paste(rep("Set",length(featureSets)),seq(1:length(featureSets)),sep="")
capture.output(featureSets, file = "All_Disease_Louvain_Clusters.txt")

dd_ALL_t<-t(dd_ALL)
dd_ALL_t<-data.frame(dd_ALL_t)
names(dd_ALL_t)<-dd_ALL_t[nrow(dd_ALL_t),]
dd_ALL_t<-dd_ALL_t[-113,]	
dd_ALL_t<-dd_ALL_t[-nrow(dd_ALL_t),]

library(data.table)
Sample_ID_tmp<-substr(rownames(dd_ALL_t),1,7)
dd_Sample_ID_tmp<-data.frame(Sample_ID_tmp)
tbl_Sample_ID<-setDT(dd_Sample_ID_tmp)
tbl_clinical<-setDT(dd_clinical_grps)
setkey(tbl_Sample_ID,Sample_ID_tmp)
setkey(tbl_clinical,Sample_ID)
tbl_samp_clinical <-tbl_Sample_ID[tbl_clinical, nomatch=0]
tbl_samp_clinical_sorted<-tbl_samp_clinical[match(substr(rownames(dd_ALL_t),1,7),tbl_samp_clinical$Sample_ID_tmp),]

#Use BNP,CRP, EOS groups
dd_ALL<-dd
Grp<-tbl_samp_clinical_sorted$BNP_Grp	
Grp<-tbl_samp_clinical_sorted$CRP_Grp	
Grp<-tbl_samp_clinical_sorted$EOS_Grp	
Grp<-as.factor(Grp)
rownames(dd_ALL)<-dd_ALL$Feature
dd_ALL<-dd_ALL[,-113]
dd_ALL<-dd_ALL[,-ncol(dd_ALL)]
dd_ALL<-as.matrix(dd_ALL)
dd_ALL_tmp<-apply(dd_ALL,2,as.numeric)
dimnames(dd_ALL_tmp)<-dimnames(dd_ALL)		
mod_mat_tmp<-data.frame(Group=Grp,row.names=dimnames(dd_ALL_tmp)[[2]])
mod_mat<-model.matrix(~Group, data = mod_mat_tmp)
mod_mat<-model.matrix(~0+Group, data = mod_mat_tmp)
colnames(mod_mat) <- levels(Grp) 	
#
	
#Use Disease groups, Asthma,COPD,Pneumonia, HF	
Grp<-c()
for(i in 1:nrow(dd_ALL_t)){
Grp[i]<-paste(gsub(substr(rownames(dd_ALL_t)[i],1,8),"",rownames(dd_ALL_t)[i]))
}
Grp<-as.factor(Grp)
rownames(dd_ALL)<-dd_ALL$Feature
dd_ALL<-dd_ALL[,-113]
dd_ALL<-dd_ALL[,-ncol(dd_ALL)]
dd_ALL<-as.matrix(dd_ALL)
dd_ALL_tmp<-apply(dd_ALL,2,as.numeric)
dimnames(dd_ALL_tmp)<-dimnames(dd_ALL)	
mod_mat_tmp<-data.frame(Group=Grp,row.names=dimnames(dd_ALL_tmp)[[2]])
mod_mat<-model.matrix(~Group, data = mod_mat_tmp)
mod_mat<-model.matrix(~0+Group, data = mod_mat_tmp)
colnames(mod_mat) <- levels(Grp) 
#
	
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)

	

fit <- lmFit(dd_ALL_tmp, mod_mat)
fit <- eBayes(fit)
gsva_es <- gsva(dd_ALL_tmp, featureSets, mx.diff=TRUE)
#gsva_es <- gsva(dd_ALL, featureSets, method=c("ssgsea"),mx.diff=TRUE)

fit <- lmFit(gsva_es, mod_mat)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.05)
summary(res)
topTable(fit)
dd_top_Asthma<-topTable(fit,coef="Acute.Asthma")
dd_top_COPD<-topTable(fit,coef="Acute.COPD")
dd_top_HF<-topTable(fit,coef="HF")
dd_top_Pneumonia<-topTable(fit,coef="Pneumonia")
dd_res_gsva<-data.frame(t(gsva_es),Grp)

write.table(dd_top_Asthma,file="Asthma_Enrichment_Sets.csv",sep=",")
write.table(dd_top_COPD,file="COPD_Enrichment_Sets.csv",sep=",")
write.table(dd_top_HF,file="HF_Enrichment_Sets.csv",sep=",")
write.table(dd_top_Pneumonia,file="Pneumonia_Enrichment_Sets.csv",sep=",")

	
aggregate(dd_res_gsva[,-ncol(dd_res_gsva)],list(dd_res_gsva$Grp),mean)
library(pROC)
plot(roc(dd_res_gsva$Grp=="High",dd_res_gsva$Set5),print.auc=TRUE,main="BNP High versus Low")
plot(roc(dd_res_gsva$Grp=="High",dd_res_gsva$Set2),print.auc=TRUE,main="CRP High versus Low")	
plot(roc(dd_res_gsva$Grp=="High",dd_res_gsva$Set1),print.auc=TRUE,main="EOS High versus Low")
	
sc1<-roc(dd_res_gsva$Grp=="Acute.Asthma",dd_res_gsva$Set1)
sc2<-roc(dd_res_gsva$Grp=="Acute.Asthma",dd_res_gsva$Set2)
sc3<-roc(dd_res_gsva$Grp=="Acute.Asthma",dd_res_gsva$Set3)
sc4<-roc(dd_res_gsva$Grp=="Acute.Asthma",dd_res_gsva$Set4)
sc5<-roc(dd_res_gsva$Grp=="Acute.Asthma",dd_res_gsva$Set5)
plot(sc1,print.auc=TRUE)
plot(sc2, add=TRUE,print.auc=TRUE,col=c("blue"))
plot(sc3, add=TRUE,print.auc=TRUE,col=c("green"))
plot(sc4, add=TRUE,print.auc=TRUE,col=c("red"))
plot(sc5, add=TRUE,print.auc=TRUE,col=c("orange"))
legend("bottomright", legend=c("Set1","Set2","Set3","Set4","Set5"),col=c("black", "blue","green","red","orange"),lwd=0.5)

plot(roc(dd_res_gsva$Grp=="HF",dd_res_gsva$Set2),print.auc=TRUE)
dd_set2<-data.frame(dd_res_gsva$Set2,dd_res_gsva$Grp)
names(dd_set2)<-c("Enricment_Set2","Grp")
for(i in 1:length(dd_set2$Grp)){
if(dd_set2$Grp[i]=="HF") {dd_set2$Grp2[i]<-"HF"}
if(dd_set2$Grp[i]!="HF") {dd_set2$Grp2[i]<-"Non HF"}	
}
library(dabestr)
unpaired_mean_diff <- dabest(dd_set2, Grp2, Enricment_Set2,
                             idx = c("HF", "Non HF"),
                             paired = FALSE) %>% 
                      mean_diff()
plot(unpaired_mean_diff)

plot(roc(dd_res_gsva$Grp=="Acute.Asthma",dd_res_gsva$Set2),print.auc=TRUE)
plot(roc(dd_res_gsva$Grp=="Acute.COPD",dd_res_gsva$Set6),print.auc=TRUE)
plot(roc(dd_res_gsva$Grp=="Pneumonia",dd_res_gsva$Set3),print.auc=TRUE)

plot(roc(dd_res_gsva$Grp=="Acute.Asthma",dd_res_gsva$Set3),print.auc=TRUE)
plot(roc(dd_res_gsva$Grp=="Acute.Asthma",dd_res_gsva$Set1),print.auc=TRUE)

plot(roc(dd_res_gsva$Grp=="HF",dd_res_gsva$Set3),print.auc=TRUE)	
plot(roc(dd_res_gsva$Grp=="HF",dd_res_gsva$Set5),print.auc=TRUE)	
plot(roc(dd_res_gsva$Grp=="HF",dd_res_gsva$Set2),print.auc=TRUE)		

plot(roc(dd_res_gsva$Grp=="Acute.COPD",dd_res_gsva$Set3),print.auc=TRUE)
plot(roc(dd_res_gsva$Grp=="Pneumonia",dd_res_gsva$Set2),print.auc=TRUE)	
#
