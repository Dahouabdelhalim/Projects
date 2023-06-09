library(plyr)

data1#bacteria data
data2#fungi data
# bactname<-"C:\\\\Users\\\\33281\\\\Desktop\\\\cp.xlsx"
# bactname<-read.xlsx(bactname,1)
# funginame<-"C:\\\\Users\\\\33281\\\\Desktop\\\\cp.xlsx"
# funginame<-read.xlsx(funginame,1)
# yin#设定分组变量
# funginame.<-as.vector(funginame[,1])
# funginame.phy<-substr(funginame.,28,38)

# bactname.phy
# funginame.phy
# 
# fix(fungidata)
# fix(bactdata)

#restoring data; removing rare species
names(data1)<-c(1:length(data1))
names(data2)<-c(1:length(data2))

c<-c()
a<-1
while(a<=length(data1)){
  alpha<-length(which(data1[,a]==0))
  if(alpha<=20){
    c<-c(c,a)}
  a<-a+1
}

data1.<-data1[,c]
length(data1.)

c<-c()
a<-1
while(a<=length(data2)){
  alpha<-length(which(data2[,a]==0))
  if(alpha<=20){
    c<-c(c,a)}
  a<-a+1
}

data2.<-data2[,c]
length(data2.)

#fungi RII
fa<-data2/rowSums(data2)

fungi.Qilan.loose.in<-(fa[1,]+fa[2,]+fa[3,])/3
fungi.Qilan.loose.out<-(fa[4,]+fa[5,]+fa[6,])/3
fungi.Qilan.tight.in<-(fa[7,]+fa[8,]+fa[9,])/3
fungi.Qilan.tight.out<-(fa[10,]+fa[11,]+fa[12,])/3
fungi.Tian.loose.in<-(fa[13,]+fa[14,]+fa[15,])/3
fungi.Tian.loose.out<-(fa[16,]+fa[17,]+fa[18,])/3
fungi.Tian.tight.in<-(fa[19,]+fa[20,]+fa[21,])/3
fungi.Tian.tight.out<-(fa[22,]+fa[23,]+fa[24,])/3

RII.qilian.loose<-(fungi.Qilan.loose.in-fungi.Qilan.loose.out)/(fungi.Qilan.loose.in+fungi.Qilan.loose.out)
RII.qilian.tight<-(fungi.Qilan.tight.in-fungi.Qilan.tight.out)/(fungi.Qilan.tight.in+fungi.Qilan.tight.out)
RII.tian.loose<-(fungi.Tian.loose.in-fungi.Tian.loose.out)/(fungi.Tian.loose.in+fungi.Tian.loose.out)
RII.tian.tight<-(fungi.Tian.tight.in-fungi.Tian.tight.out)/(fungi.Tian.tight.in+fungi.Tian.tight.out)

fungi.RII<-rbind(RII.qilian.loose,RII.qilian.tight,RII.tian.loose,RII.tian.tight)
na.<-is.na(fungi.RII)
fungi.RII.re<-t(na.omit(t(fungi.RII)))#去除稀有种


#bacteria RII
ba<-data1/rowSums(data1)

bact.Qilan.loose.in<-(ba[1,]+ba[2,]+ba[3,])/3
bact.Qilan.loose.out<-(ba[4,]+ba[5,]+ba[6,])/3
bact.Qilan.tight.in<-(ba[7,]+ba[8,]+ba[9,])/3
bact.Qilan.tight.out<-(ba[10,]+ba[11,]+ba[12,])/3
bact.Tian.loose.in<-(ba[13,]+ba[14,]+ba[15,])/3
bact.Tian.loose.out<-(ba[16,]+ba[17,]+ba[18,])/3
bact.Tian.tight.in<-(ba[19,]+ba[20,]+ba[21,])/3
bact.Tian.tight.out<-(ba[22,]+ba[23,]+ba[24,])/3

RII.qilian.loose<-(bact.Qilan.loose.in-bact.Qilan.loose.out)/(bact.Qilan.loose.in+bact.Qilan.loose.out)
RII.qilian.tight<-(bact.Qilan.tight.in-bact.Qilan.tight.out)/(bact.Qilan.tight.in+bact.Qilan.tight.out)
RII.tian.loose<-(bact.Tian.loose.in-bact.Tian.loose.out)/(bact.Tian.loose.in+bact.Tian.loose.out)
RII.tian.tight<-(bact.Tian.tight.in-bact.Tian.tight.out)/(bact.Tian.tight.in+bact.Tian.tight.out)

bact.RII<-rbind(RII.qilian.loose,RII.qilian.tight,RII.tian.loose,RII.tian.tight)
na.<-is.na(bact.RII)
bact.RII.re<-t(na.omit(t(bact.RII)))#去除稀有种


#clust
fungi.clust<-t(fungi.RII.re)
dis<-dist(fungi.clust, method = "euclidean") 
fit.clust<-hclust(dis, method="ward.D2")
plot(fit.clust)
rect.hclust(fit.clust, k=6, border="red")
groups.fungi <- cutree(fit.clust, k=6)
groups.fungi<-as.data.frame(groups.fungi)
order(groups.fungi)
result.fungi<-data.frame(groups.fungi,fungi.clust)


bact.clust<-t(bact.RII.re)
dis<-dist(bact.clust, method = "euclidean") 
fit.clust<-hclust(dis, method="ward.D2")
plot(fit.clust)
rect.hclust(fit.clust, k=6, border="red")
groups.bact <- cutree(fit.clust, k=6)
groups.bact<-as.data.frame(groups.bact)
order(groups.bact)
results.bact<-data.frame(groups.bact,bact.clust)


# #PCA
# library("FactoMineR")
# library("factoextra")
# 
# fungi.RII.pca<-PCA(fungi.clust,scale.unit = TRUE,ncp=2, graph = TRUE)
# eigenvalues<-fungi.RII.pca$eig
# fungi.RII.pca1<-fungi.RII.pca$ind$coord
# fungi.RII.pca2<-fungi.RII.pca$var$coord
# 
# barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
#         main = "Variances",
#         xlab = "Principal Components",
#         ylab = "Percentage of variances",
#         col ="steelblue")
# lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
#       type="b", pch=19, col = "red")
# 
# 
# fungi.RII.pca<-PCA(bact.clust,scale.unit = TRUE,ncp=2, graph = TRUE)
# eigenvalues<-fungi.RII.pca$eig
# fungi.RII.pca1<-fungi.RII.pca$ind$coord
# fungi.RII.pca2<-fungi.RII.pca$var$coord
# 
# barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
#         main = "Variances",
#         xlab = "Principal Components",
#         ylab = "Percentage of variances",
#         col ="steelblue")
# lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
#       type="b", pch=19, col = "red")



#################################################################################################3
#fungi group data
names.fungi.g<-row.names(result.fungi)
fun.1<-names.fungi.g[which(result.fungi$groups.fungi==1)]
fun.2<-names.fungi.g[which(result.fungi$groups.fungi==2)]
fun.3<-names.fungi.g[which(result.fungi$groups.fungi==3)]
fun.4<-names.fungi.g[which(result.fungi$groups.fungi==4)]
fun.5<-names.fungi.g[which(result.fungi$groups.fungi==5)]
fun.6<-names.fungi.g[which(result.fungi$groups.fungi==6)]

fun.1<-as.numeric(fun.1)
fun.2<-as.numeric(fun.2)
fun.3<-as.numeric(fun.3)
fun.4<-as.numeric(fun.4)
fun.5<-as.numeric(fun.5)
fun.6<-as.numeric(fun.6)

fungi1<-data2[,fun.1]
fungi1.<-data.frame(rowSums(fungi1))
fungi2<-data2[,fun.2]
fungi2.<-data.frame(rowSums(fungi2))
fungi3<-data2[,fun.3]
fungi3.<-data.frame(rowSums(fungi3))
fungi4<-data2[,fun.4]
fungi4.<-data.frame(rowSums(fungi4))
fungi5<-data2[,fun.5]
fungi5.<-data.frame(rowSums(fungi5))
fungi6<-data2[,fun.6]
fungi6.<-data.frame(rowSums(fungi6))

fr1<-rich(fungi1)
fr2<-rich(fungi2)
fr3<-rich(fungi3)
fr4<-rich(fungi4)
fr5<-rich(fungi5)
fr6<-rich(fungi6)

#bacteria group data
names.bact.g<-row.names(results.bact)
bac.1<-names.bact.g[which(results.bact$groups.bact==1)]
bac.2<-names.bact.g[which(results.bact$groups.bact==2)]
bac.3<-names.bact.g[which(results.bact$groups.bact==3)]
bac.4<-names.bact.g[which(results.bact$groups.bact==4)]
bac.5<-names.bact.g[which(results.bact$groups.bact==5)]
bac.6<-names.bact.g[which(results.bact$groups.bact==6)]

bac.1<-as.numeric(bac.1)
bac.2<-as.numeric(bac.2)
bac.3<-as.numeric(bac.3)
bac.4<-as.numeric(bac.4)
bac.5<-as.numeric(bac.5)
bac.6<-as.numeric(bac.6)

bact1<-data1[,bac.1]
bact1.<-data.frame(rowSums(bact1))
bact2<-data1[,bac.2]
bact2.<-data.frame(rowSums(bact2))
bact3<-data1[,bac.3]
bact3.<-data.frame(rowSums(bact3))
bact4<-data1[,bac.4]
bact4.<-data.frame(rowSums(bact4))
bact5<-data1[,bac.5]
bact5.<-data.frame(rowSums(bact5))
bact6<-data1[,bac.6]
bact6.<-data.frame(rowSums(bact6))

br1<-rich(bact1)
br2<-rich(bact2)
br3<-rich(bact3)
br4<-rich(bact4)
br5<-rich(bact5)
br6<-rich(bact6)

# ##############
# bact1.phy<-bactname.phy[as.numeric(names(bact1))]
# bact1.ab<-colSums(bact1)
# bact1.phy<-data.frame(bact1.phy,bact1.ab)
# 
# bact2.phy<-bactname.phy[as.numeric(names(bact2))]
# bact2.ab<-colSums(bact2)
# bact2.phy<-data.frame(bact2.phy,bact2.ab)
# 
# bact3.phy<-bactname.phy[as.numeric(names(bact3))]
# bact3.ab<-colSums(bact3)
# bact3.phy<-data.frame(bact3.phy,bact3.ab)
# 
# bact4.phy<-bactname.phy[as.numeric(names(bact4))]
# bact4.ab<-colSums(bact4)
# bact4.phy<-data.frame(bact4.phy,bact4.ab)
# 
# bact5.phy<-bactname.phy[as.numeric(names(bact5))]
# bact5.ab<-colSums(bact5)
# bact5.phy<-data.frame(bact5.phy,bact5.ab)
# 
# bact6.phy<-bactname.phy[as.numeric(names(bact6))]
# bact6.ab<-colSums(bact6)
# bact6.phy<-data.frame(bact6.phy,bact6.ab)
# #
# bact1.phy.<-bact1.phy[order(bact1.phy$bact1.phy),]
# bact1.phy.data<-data.frame()
# for(i in 1:length(levels(bact1.phy$bact1.phy))){
#   a<-bact1.phy[which(bact1.phy$bact1.phy==levels(bact1.phy$bact1.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(bact1.phy$bact1.phy)[i]
#   bact1.phy..<-data.frame(n,b)
#   bact1.phy.data<-rbind(bact1.phy.data,bact1.phy..)
# }
# 
# bact2.phy.<-bact2.phy[order(bact2.phy$bact2.phy),]
# bact2.phy.data<-data.frame()
# for(i in 1:length(levels(bact2.phy$bact2.phy))){
#   a<-bact2.phy[which(bact2.phy$bact2.phy==levels(bact2.phy$bact2.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(bact2.phy$bact2.phy)[i]
#   bact2.phy..<-data.frame(n,b)
#   bact2.phy.data<-rbind(bact2.phy.data,bact2.phy..)
# }
# 
# bact3.phy.<-bact3.phy[order(bact3.phy$bact3.phy),]
# bact3.phy.data<-data.frame()
# for(i in 1:length(levels(bact3.phy$bact3.phy))){
#   a<-bact3.phy[which(bact3.phy$bact3.phy==levels(bact3.phy$bact3.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(bact3.phy$bact3.phy)[i]
#   bact3.phy..<-data.frame(n,b)
#   bact3.phy.data<-rbind(bact3.phy.data,bact3.phy..)
# }
# 
# bact4.phy.<-bact4.phy[order(bact4.phy$bact4.phy),]
# bact4.phy.data<-data.frame()
# for(i in 1:length(levels(bact4.phy$bact4.phy))){
#   a<-bact4.phy[which(bact4.phy$bact4.phy==levels(bact4.phy$bact4.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(bact4.phy$bact4.phy)[i]
#   bact4.phy..<-data.frame(n,b)
#   bact4.phy.data<-rbind(bact4.phy.data,bact4.phy..)
# }
# 
# bact5.phy.<-bact5.phy[order(bact5.phy$bact5.phy),]
# bact5.phy.data<-data.frame()
# for(i in 1:length(levels(bact5.phy$bact5.phy))){
#   a<-bact5.phy[which(bact5.phy$bact5.phy==levels(bact5.phy$bact5.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(bact5.phy$bact5.phy)[i]
#   bact5.phy..<-data.frame(n,b)
#   bact5.phy.data<-rbind(bact5.phy.data,bact5.phy..)
# }
# 
# bact6.phy.<-bact6.phy[order(bact6.phy$bact6.phy),]
# bact6.phy.data<-data.frame()
# for(i in 1:length(levels(bact6.phy$bact6.phy))){
#   a<-bact6.phy[which(bact6.phy$bact6.phy==levels(bact6.phy$bact6.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(bact6.phy$bact6.phy)[i]
#   bact6.phy..<-data.frame(n,b)
#   bact6.phy.data<-rbind(bact6.phy.data,bact6.phy..)
# }
# 
# bact.phy.data<-rbind(bact1.phy.data,bact2.phy.data,bact3.phy.data,bact4.phy.data,bact5.phy.data,bact6.phy.data)
# group.phy<-c(rep(1,nrow(bact1.phy.data)),rep(2,nrow(bact2.phy.data)),rep(3,nrow(bact3.phy.data)),
#   rep(4,nrow(bact4.phy.data)),rep(5,nrow(bact5.phy.data)),rep(6,nrow(bact6.phy.data)))
# bact.phy.data<-data.frame(group.phy,bact.phy.data)
# 
# bact1.phy.data$n
# duplicated(levels(bact1.phy.data$n),levels(bact2.phy.data$n))
# ####################
# #bact.phy
# names(bact.phy.data)<-c("group.phy","name.phy","b" )
# 
# i<-1
# for(i in 1:6){
#   a<-which(bact.phy.data$group.phy==i)
#   b<-bact.phy.data$b[a]
#   bact.phy.data$b[a]<-b/sum(bact.phy.data$b[a])
# }
# ("JL???ETNP???Z","SHA-109; ","SM2F11; c","WCHB1-60;")
# fix(bact.phy.data)
# bact.phy.data$name.phy[c(14,15,11,17,33,35,52,53,55,64,79,92)]<-as.factor(rep("Others",12))
# library(ggplot2)
# bact.phy<-ggplot()+
#   geom_bar(data=bact.phy.data,aes(x=group.phy,y=b,fill=factor(name.phy)),stat="identity",position="stack",colour="black",width = 0.75)+ 
#   geom_hline(yintercept=0)+ 
#   geom_vline(xintercept=0)+ 
#   #scale_fill_continuous(low="red",high = "white")+
#   scale_x_continuous(limits = c(0,7),
#                      breaks=c(1,2,3,4,5,6),
#                      labels=c("Group 1","Group 2","Group 3","Group 4","Group 5","Group 6"),
#                      expand = c(0,0))+ 
#   scale_y_continuous(expand = c(0,0))+
#   labs(title=NULL, x=NULL, y="Bacteria abundancce")+
#   theme(panel.grid=element_blank(),
#         panel.background=element_blank(),
#         panel.border=element_blank(),
#         axis.line=element_line(color = "black"),
#         legend.position="right",
#         legend.key = element_blank(),
#         legend.title= element_blank(),
#         legend.text=element_text(size=14),
#         axis.title=element_text(size=20,face="bold"),
#         axis.text=element_text(size=15,face = "bold"),
#         axis.text.x = element_text(angle=45, hjust=1))
# 
# ggsave("W:\\\\论文发表\\\\bact_phy.pdf",bact.phy,width = 10,height = 6,dpi=300)
# 
# ####fungi phy group
# fungi1.phy<-funginame.phy[as.numeric(names(fungi1))]
# fungi1.ab<-colSums(fungi1)
# fungi1.phy<-data.frame(fungi1.phy,fungi1.ab)
# 
# fungi2.phy<-funginame.phy[as.numeric(names(fungi2))]
# fungi2.ab<-colSums(fungi2)
# fungi2.phy<-data.frame(fungi2.phy,fungi2.ab)
# 
# fungi3.phy<-funginame.phy[as.numeric(names(fungi3))]
# fungi3.ab<-colSums(fungi3)
# fungi3.phy<-data.frame(fungi3.phy,fungi3.ab)
# 
# fungi4.phy<-funginame.phy[as.numeric(names(fungi4))]
# fungi4.ab<-colSums(fungi4)
# fungi4.phy<-data.frame(fungi4.phy,fungi4.ab)
# 
# fungi5.phy<-funginame.phy[as.numeric(names(fungi5))]
# fungi5.ab<-colSums(fungi5)
# fungi5.phy<-data.frame(fungi5.phy,fungi5.ab)
# 
# fungi6.phy<-funginame.phy[as.numeric(names(fungi6))]
# fungi6.ab<-colSums(fungi6)
# fungi6.phy<-data.frame(fungi6.phy,fungi6.ab)
# #
# fungi1.phy.<-fungi1.phy[order(fungi1.phy$fungi1.phy),]
# fungi1.phy.data<-data.frame()
# for(i in 1:length(levels(fungi1.phy$fungi1.phy))){
#   a<-fungi1.phy[which(fungi1.phy$fungi1.phy==levels(fungi1.phy$fungi1.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(fungi1.phy$fungi1.phy)[i]
#   fungi1.phy..<-data.frame(n,b)
#   fungi1.phy.data<-rbind(fungi1.phy.data,fungi1.phy..)
# }
# 
# fungi2.phy.<-fungi2.phy[order(fungi2.phy$fungi2.phy),]
# fungi2.phy.data<-data.frame()
# for(i in 1:length(levels(fungi2.phy$fungi2.phy))){
#   a<-fungi2.phy[which(fungi2.phy$fungi2.phy==levels(fungi2.phy$fungi2.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(fungi2.phy$fungi2.phy)[i]
#   fungi2.phy..<-data.frame(n,b)
#   fungi2.phy.data<-rbind(fungi2.phy.data,fungi2.phy..)
# }
# 
# fungi3.phy.<-fungi3.phy[order(fungi3.phy$fungi3.phy),]
# fungi3.phy.data<-data.frame()
# for(i in 1:length(levels(fungi3.phy$fungi3.phy))){
#   a<-fungi3.phy[which(fungi3.phy$fungi3.phy==levels(fungi3.phy$fungi3.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(fungi3.phy$fungi3.phy)[i]
#   fungi3.phy..<-data.frame(n,b)
#   fungi3.phy.data<-rbind(fungi3.phy.data,fungi3.phy..)
# }
# 
# fungi4.phy.<-fungi4.phy[order(fungi4.phy$fungi4.phy),]
# fungi4.phy.data<-data.frame()
# for(i in 1:length(levels(fungi4.phy$fungi4.phy))){
#   a<-fungi4.phy[which(fungi4.phy$fungi4.phy==levels(fungi4.phy$fungi4.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(fungi4.phy$fungi4.phy)[i]
#   fungi4.phy..<-data.frame(n,b)
#   fungi4.phy.data<-rbind(fungi4.phy.data,fungi4.phy..)
# }
# 
# fungi5.phy.<-fungi5.phy[order(fungi5.phy$fungi5.phy),]
# fungi5.phy.data<-data.frame()
# for(i in 1:length(levels(fungi5.phy$fungi5.phy))){
#   a<-fungi5.phy[which(fungi5.phy$fungi5.phy==levels(fungi5.phy$fungi5.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(fungi5.phy$fungi5.phy)[i]
#   fungi5.phy..<-data.frame(n,b)
#   fungi5.phy.data<-rbind(fungi5.phy.data,fungi5.phy..)
# }
# 
# fungi6.phy.<-fungi6.phy[order(fungi6.phy$fungi6.phy),]
# fungi6.phy.data<-data.frame()
# for(i in 1:length(levels(fungi6.phy$fungi6.phy))){
#   a<-fungi6.phy[which(fungi6.phy$fungi6.phy==levels(fungi6.phy$fungi6.phy)[i]),2]
#   b<-sum(a)
#   n<-levels(fungi6.phy$fungi6.phy)[i]
#   fungi6.phy..<-data.frame(n,b)
#   fungi6.phy.data<-rbind(fungi6.phy.data,fungi6.phy..)
# }
# 
# fungi.phy.data<-rbind(fungi1.phy.data,fungi2.phy.data,fungi3.phy.data,fungi4.phy.data,fungi5.phy.data,fungi6.phy.data)
# group.phy<-c(rep(1,nrow(fungi1.phy.data)),rep(2,nrow(fungi2.phy.data)),rep(3,nrow(fungi3.phy.data)),
#              rep(4,nrow(fungi4.phy.data)),rep(5,nrow(fungi5.phy.data)),rep(6,nrow(fungi6.phy.data)))
# fungi.phy.data<-data.frame(group.phy,fungi.phy.data)
# 

####################
# #bact.phy
# names(fungi.phy.data)<-c("group.phy","name.phy","b" )
# display.brewer.all(n=5, exact.n=FALSE)
# brewer.pal.info
# 
# i<-1
# for(i in 1:6){
#   a<-which(fungi.phy.data$group.phy==i)
#   b<-fungi.phy.data$b[a]
#   fungi.phy.data$b[a]<-b/sum(fungi.phy.data$b[a])
# }
# 
# levels(fungi.phy.data$name.phy)
# 
# fungi.phy<-ggplot() +
#   geom_bar(data=fungi.phy.data,aes(x=group.phy,y=b,fill=factor(name.phy)),stat="identity",position="stack",colour="black",width = 0.75)+ 
#   geom_hline(yintercept=0)+ 
#   geom_vline(xintercept=0)+ 
#   scale_fill_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00')，
#                     labels = c())+
#   scale_x_continuous(limits = c(0,7),
#                      breaks=c(1,2,3,4,5,6),
#                      labels=c("Group 1","Group 2","Group 3","Group 4","Group 5","Group 6"),
#                      expand = c(0,0))+ 
#   scale_y_continuous(expand = c(0,0))+
#   labs(title=NULL, x=NULL, y="Bacteria abundancce")+
#   theme(panel.grid=element_blank(),
#         panel.background=element_blank(),
#         panel.border=element_blank(),
#         axis.line=element_line(color = "black"),
#         legend.position="right",
#         legend.key = element_blank(),
#         legend.title= element_blank(),
#         legend.text=element_text(size=14),
#         axis.title=element_text(size=20,face="bold"),
#         axis.text=element_text(size=15,face = "bold"),
#         axis.text.x = element_text(angle=45, hjust=1))
# 
# ggsave("W:\\\\论文发表\\\\fungi_phy.pdf",fungi.phy,width = 10,height = 6,dpi=300)



#fungi RII
RII.qilian.loose.1<-(fungi1.[c(1:3),]-fungi1.[c(4:6),])/(fungi1.[c(1:3),]+fungi1.[c(4:6),])
RII.qilian.tight.1<-(fungi1.[c(7:9),]-fungi1.[c(10:12),])/(fungi1.[c(7:9),]+fungi1.[c(10:12),])
RII.tian.loose.1<-(fungi1.[c(13:15),]-fungi1.[c(16:18),])/(fungi1.[c(13:15),]+fungi1.[c(16:18),])
RII.tian.tight.1<-(fungi1.[c(19:21),]-fungi1.[c(22:24),])/(fungi1.[c(19:21),]+fungi1.[c(22:24),])
f.1<-c(RII.qilian.loose.1,RII.qilian.tight.1,RII.tian.loose.1,RII.tian.tight.1)

RII.qilian.loose.2<-(fungi2.[c(1:3),]-fungi2.[c(4:6),])/(fungi2.[c(1:3),]+fungi2.[c(4:6),])
RII.qilian.tight.2<-(fungi2.[c(7:9),]-fungi2.[c(10:12),])/(fungi2.[c(7:9),]+fungi2.[c(10:12),])
RII.tian.loose.2<-(fungi2.[c(13:15),]-fungi2.[c(16:18),])/(fungi2.[c(13:15),]+fungi2.[c(16:18),])
RII.tian.tight.2<-(fungi2.[c(19:21),]-fungi2.[c(22:24),])/(fungi2.[c(19:21),]+fungi2.[c(22:24),])
f.2<-c(RII.qilian.loose.2,RII.qilian.tight.2,RII.tian.loose.2,RII.tian.tight.2)

RII.qilian.loose.3<-(fungi3.[c(1:3),]-fungi3.[c(4:6),])/(fungi3.[c(1:3),]+fungi3.[c(4:6),])
RII.qilian.tight.3<-(fungi3.[c(7:9),]-fungi3.[c(10:12),])/(fungi3.[c(7:9),]+fungi3.[c(10:12),])
RII.tian.loose.3<-(fungi3.[c(13:15),]-fungi3.[c(16:18),])/(fungi3.[c(13:15),]+fungi3.[c(16:18),])
RII.tian.tight.3<-(fungi3.[c(19:21),]-fungi3.[c(22:24),])/(fungi3.[c(19:21),]+fungi3.[c(22:24),])
f.3<-c(RII.qilian.loose.3,RII.qilian.tight.3,RII.tian.loose.3,RII.tian.tight.3)

RII.qilian.loose.4<-(fungi4.[c(1:3),]-fungi4.[c(4:6),])/(fungi4.[c(1:3),]+fungi4.[c(4:6),])
RII.qilian.tight.4<-(fungi4.[c(7:9),]-fungi4.[c(10:12),])/(fungi4.[c(7:9),]+fungi4.[c(10:12),])
RII.tian.loose.4<-(fungi4.[c(13:15),]-fungi4.[c(16:18),])/(fungi4.[c(13:15),]+fungi4.[c(16:18),])
RII.tian.tight.4<-(fungi4.[c(19:21),]-fungi4.[c(22:24),])/(fungi4.[c(19:21),]+fungi4.[c(22:24),])
f.4<-c(RII.qilian.loose.4,RII.qilian.tight.4,RII.tian.loose.4,RII.tian.tight.4)

RII.qilian.loose.5<-(fungi5.[c(1:3),]-fungi5.[c(4:6),])/(fungi5.[c(1:3),]+fungi5.[c(4:6),])
RII.qilian.tight.5<-(fungi5.[c(7:9),]-fungi5.[c(10:12),])/(fungi5.[c(7:9),]+fungi5.[c(10:12),])
RII.tian.loose.5<-(fungi5.[c(13:15),]-fungi5.[c(16:18),])/(fungi5.[c(13:15),]+fungi5.[c(16:18),])
RII.tian.tight.5<-(fungi5.[c(19:21),]-fungi5.[c(22:24),])/(fungi5.[c(19:21),]+fungi5.[c(22:24),])
f.5<-c(RII.qilian.loose.5,RII.qilian.tight.5,RII.tian.loose.5,RII.tian.tight.5)

RII.qilian.loose.6<-(fungi6.[c(1:3),]-fungi6.[c(4:6),])/(fungi6.[c(1:3),]+fungi6.[c(4:6),])
RII.qilian.tight.6<-(fungi6.[c(7:9),]-fungi6.[c(10:12),])/(fungi6.[c(7:9),]+fungi6.[c(10:12),])
RII.tian.loose.6<-(fungi6.[c(13:15),]-fungi6.[c(16:18),])/(fungi6.[c(13:15),]+fungi6.[c(16:18),])
RII.tian.tight.6<-(fungi6.[c(19:21),]-fungi6.[c(22:24),])/(fungi6.[c(19:21),]+fungi6.[c(22:24),])
f.6<-c(RII.qilian.loose.6,RII.qilian.tight.6,RII.tian.loose.6,RII.tian.tight.6)

#fungi richness RII
RII.qilian.loose.1<-(fr1[c(1:3)]-fr1[c(4:6)])/(fr1[c(1:3)]+fr1[c(4:6)])
RII.qilian.tight.1<-(fr1[c(7:9)]-fr1[c(10:12)])/(fr1[c(7:9)]+fr1[c(10:12)])
RII.tian.loose.1<-(fr1[c(13:15)]-fr1[c(16:18)])/(fr1[c(13:15)]+fr1[c(16:18)])
RII.tian.tight.1<-(fr1[c(19:21)]-fr1[c(22:24)])/(fr1[c(19:21)]+fr1[c(22:24)])
f.1.<-c(RII.qilian.loose.1,RII.qilian.tight.1,RII.tian.loose.1,RII.tian.tight.1)

RII.qilian.loose.2<-(fr2[c(1:3)]-fr2[c(4:6)])/(fr2[c(1:3)]+fr2[c(4:6)])
RII.qilian.tight.2<-(fr2[c(7:9)]-fr2[c(10:12)])/(fr2[c(7:9)]+fr2[c(10:12)])
RII.tian.loose.2<-(fr2[c(13:15)]-fr2[c(16:18)])/(fr2[c(13:15)]+fr2[c(16:18)])
RII.tian.tight.2<-(fr2[c(19:21)]-fr2[c(22:24)])/(fr2[c(19:21)]+fr2[c(22:24)])
f.2.<-c(RII.qilian.loose.2,RII.qilian.tight.2,RII.tian.loose.2,RII.tian.tight.2)

RII.qilian.loose.3<-(fr3[c(1:3)]-fr3[c(4:6)])/(fr3[c(1:3)]+fr3[c(4:6)])
RII.qilian.tight.3<-(fr3[c(7:9)]-fr3[c(10:12)])/(fr3[c(7:9)]+fr3[c(10:12)])
RII.tian.loose.3<-(fr3[c(13:15)]-fr3[c(16:18)])/(fr3[c(13:15)]+fr3[c(16:18)])
RII.tian.tight.3<-(fr3[c(19:21)]-fr3[c(22:24)])/(fr3[c(19:21)]+fr3[c(22:24)])
f.3.<-c(RII.qilian.loose.3,RII.qilian.tight.3,RII.tian.loose.3,RII.tian.tight.3)

RII.qilian.loose.4<-(fr4[c(1:3)]-fr4[c(4:6)])/(fr4[c(1:3)]+fr4[c(4:6)])
RII.qilian.tight.4<-(fr4[c(7:9)]-fr4[c(10:12)])/(fr4[c(7:9)]+fr4[c(10:12)])
RII.tian.loose.4<-(fr4[c(13:15)]-fr4[c(16:18)])/(fr4[c(13:15)]+fr4[c(16:18)])
RII.tian.tight.4<-(fr4[c(19:21)]-fr4[c(22:24)])/(fr4[c(19:21)]+fr4[c(22:24)])
f.4.<-c(RII.qilian.loose.4,RII.qilian.tight.4,RII.tian.loose.4,RII.tian.tight.4)

RII.qilian.loose.5<-(fr5[c(1:3)]-fr5[c(4:6)])/(fr5[c(1:3)]+fr5[c(4:6)])
RII.qilian.tight.5<-(fr5[c(7:9)]-fr5[c(10:12)])/(fr5[c(7:9)]+fr5[c(10:12)])
RII.tian.loose.5<-(fr5[c(13:15)]-fr5[c(16:18)])/(fr5[c(13:15)]+fr5[c(16:18)])
RII.tian.tight.5<-(fr5[c(19:21)]-fr5[c(22:24)])/(fr5[c(19:21)]+fr5[c(22:24)])
f.5.<-c(RII.qilian.loose.5,RII.qilian.tight.5,RII.tian.loose.5,RII.tian.tight.5)

RII.qilian.loose.6<-(fr6[c(1:3)]-fr6[c(4:6)])/(fr6[c(1:3)]+fr6[c(4:6)])
RII.qilian.tight.6<-(fr6[c(7:9)]-fr6[c(10:12)])/(fr6[c(7:9)]+fr6[c(10:12)])
RII.tian.loose.6<-(fr6[c(13:15)]-fr6[c(16:18)])/(fr6[c(13:15)]+fr6[c(16:18)])
RII.tian.tight.6<-(fr6[c(19:21)]-fr6[c(22:24)])/(fr6[c(19:21)]+fr6[c(22:24)])
f.6.<-c(RII.qilian.loose.6,RII.qilian.tight.6,RII.tian.loose.6,RII.tian.tight.6)

#
pheno<-rep(rep(c("A","B","A","B"),each=3),6)
site<-rep(rep(c("A","B"),each=6),6)
group<-rep(c("A","B","C","D","E","F"),each=12)
pheno.<-rep(c("A","B","A","B"),each=3)
site.<-rep(c("A","B"),each=6)

fungi.ab<-c(f.1,f.2,f.3,f.4,f.5,f.6)
fungi.ab.<-data.frame(site,pheno,group,fungi.ab)
fungi.ab.<-replace(fungi.ab.,is.na(fungi.ab.),0)
fungi.abs<-data.frame(site.,pheno.,f.1,f.2,f.3,f.4)
fungi.abs<-replace(fungi.abs,is.na(fungi.abs),0)


fungi.r<-c(f.1.,f.2.,f.3.,f.4.,f.5.,f.6.)
fungi.r.<-data.frame(site,pheno,group,fungi.r)
fungi.r.<-replace(fungi.r.,is.na(fungi.r.),0)
fungi.rs<-data.frame(site.,pheno.,f.1.,f.2.,f.3.,f.4.,f.5.,f.6.)
fungi.rs<-replace(fungi.rs,is.na(fungi.rs),0)

#bacteria RII
RII.qilian.loose.1<-(bact1.[c(1:3),]-bact1.[c(4:6),])/(bact1.[c(1:3),]+bact1.[c(4:6),])
RII.qilian.tight.1<-(bact1.[c(7:9),]-bact1.[c(10:12),])/(bact1.[c(7:9),]+bact1.[c(10:12),])
RII.tian.loose.1<-(bact1.[c(13:15),]-bact1.[c(16:18),])/(bact1.[c(13:15),]+bact1.[c(16:18),])
RII.tian.tight.1<-(bact1.[c(19:21),]-bact1.[c(22:24),])/(bact1.[c(19:21),]+bact1.[c(22:24),])
f.1<-c(RII.qilian.loose.1,RII.qilian.tight.1,RII.tian.loose.1,RII.tian.tight.1)

RII.qilian.loose.2<-(bact2.[c(1:3),]-bact2.[c(4:6),])/(bact2.[c(1:3),]+bact2.[c(4:6),])
RII.qilian.tight.2<-(bact2.[c(7:9),]-bact2.[c(10:12),])/(bact2.[c(7:9),]+bact2.[c(10:12),])
RII.tian.loose.2<-(bact2.[c(13:15),]-bact2.[c(16:18),])/(bact2.[c(13:15),]+bact2.[c(16:18),])
RII.tian.tight.2<-(bact2.[c(19:21),]-bact2.[c(22:24),])/(bact2.[c(19:21),]+bact2.[c(22:24),])
f.2<-c(RII.qilian.loose.2,RII.qilian.tight.2,RII.tian.loose.2,RII.tian.tight.2)

RII.qilian.loose.3<-(bact3.[c(1:3),]-bact3.[c(4:6),])/(bact3.[c(1:3),]+bact3.[c(4:6),])
RII.qilian.tight.3<-(bact3.[c(7:9),]-bact3.[c(10:12),])/(bact3.[c(7:9),]+bact3.[c(10:12),])
RII.tian.loose.3<-(bact3.[c(13:15),]-bact3.[c(16:18),])/(bact3.[c(13:15),]+bact3.[c(16:18),])
RII.tian.tight.3<-(bact3.[c(19:21),]-bact3.[c(22:24),])/(bact3.[c(19:21),]+bact3.[c(22:24),])
f.3<-c(RII.qilian.loose.3,RII.qilian.tight.3,RII.tian.loose.3,RII.tian.tight.3)

RII.qilian.loose.4<-(bact4.[c(1:3),]-bact4.[c(4:6),])/(bact4.[c(1:3),]+bact4.[c(4:6),])
RII.qilian.tight.4<-(bact4.[c(7:9),]-bact4.[c(10:12),])/(bact4.[c(7:9),]+bact4.[c(10:12),])
RII.tian.loose.4<-(bact4.[c(13:15),]-bact4.[c(16:18),])/(bact4.[c(13:15),]+bact4.[c(16:18),])
RII.tian.tight.4<-(bact4.[c(19:21),]-bact4.[c(22:24),])/(bact4.[c(19:21),]+bact4.[c(22:24),])
f.4<-c(RII.qilian.loose.4,RII.qilian.tight.4,RII.tian.loose.4,RII.tian.tight.4)

RII.qilian.loose.5<-(bact5.[c(1:3),]-bact5.[c(4:6),])/(bact5.[c(1:3),]+bact5.[c(4:6),])
RII.qilian.tight.5<-(bact5.[c(7:9),]-bact5.[c(10:12),])/(bact5.[c(7:9),]+bact5.[c(10:12),])
RII.tian.loose.5<-(bact5.[c(13:15),]-bact5.[c(16:18),])/(bact5.[c(13:15),]+bact5.[c(16:18),])
RII.tian.tight.5<-(bact5.[c(19:21),]-bact5.[c(22:24),])/(bact5.[c(19:21),]+bact5.[c(22:24),])
f.5<-c(RII.qilian.loose.5,RII.qilian.tight.5,RII.tian.loose.5,RII.tian.tight.5)

RII.qilian.loose.6<-(bact6.[c(1:3),]-bact6.[c(4:6),])/(bact6.[c(1:3),]+bact6.[c(4:6),])
RII.qilian.tight.6<-(bact6.[c(7:9),]-bact6.[c(10:12),])/(bact6.[c(7:9),]+bact6.[c(10:12),])
RII.tian.loose.6<-(bact6.[c(13:15),]-bact6.[c(16:18),])/(bact6.[c(13:15),]+bact6.[c(16:18),])
RII.tian.tight.6<-(bact6.[c(19:21),]-bact6.[c(22:24),])/(bact6.[c(19:21),]+bact6.[c(22:24),])
f.6<-c(RII.qilian.loose.6,RII.qilian.tight.6,RII.tian.loose.6,RII.tian.tight.6)

#bacterria richness RII
RII.qilian.loose.1<-(br1[c(1:3)]-br1[c(4:6)])/(br1[c(1:3)]+br1[c(4:6)])
RII.qilian.tight.1<-(br1[c(7:9)]-br1[c(10:12)])/(br1[c(7:9)]+br1[c(10:12)])
RII.tian.loose.1<-(br1[c(13:15)]-br1[c(16:18)])/(br1[c(13:15)]+br1[c(16:18)])
RII.tian.tight.1<-(br1[c(19:21)]-br1[c(22:24)])/(br1[c(19:21)]+br1[c(22:24)])
f.1.<-c(RII.qilian.loose.1,RII.qilian.tight.1,RII.tian.loose.1,RII.tian.tight.1)

RII.qilian.loose.2<-(br2[c(1:3)]-br2[c(4:6)])/(br2[c(1:3)]+br2[c(4:6)])
RII.qilian.tight.2<-(br2[c(7:9)]-br2[c(10:12)])/(br2[c(7:9)]+br2[c(10:12)])
RII.tian.loose.2<-(br2[c(13:15)]-br2[c(16:18)])/(br2[c(13:15)]+br2[c(16:18)])
RII.tian.tight.2<-(br2[c(19:21)]-br2[c(22:24)])/(br2[c(19:21)]+br2[c(22:24)])
f.2.<-c(RII.qilian.loose.2,RII.qilian.tight.2,RII.tian.loose.2,RII.tian.tight.2)

RII.qilian.loose.3<-(br3[c(1:3)]-br3[c(4:6)])/(br3[c(1:3)]+br3[c(4:6)])
RII.qilian.tight.3<-(br3[c(7:9)]-br3[c(10:12)])/(br3[c(7:9)]+br3[c(10:12)])
RII.tian.loose.3<-(br3[c(13:15)]-br3[c(16:18)])/(br3[c(13:15)]+br3[c(16:18)])
RII.tian.tight.3<-(br3[c(19:21)]-br3[c(22:24)])/(br3[c(19:21)]+br3[c(22:24)])
f.3.<-c(RII.qilian.loose.3,RII.qilian.tight.3,RII.tian.loose.3,RII.tian.tight.3)

RII.qilian.loose.4<-(br4[c(1:3)]-br4[c(4:6)])/(br4[c(1:3)]+br4[c(4:6)])
RII.qilian.tight.4<-(br4[c(7:9)]-br4[c(10:12)])/(br4[c(7:9)]+br4[c(10:12)])
RII.tian.loose.4<-(br4[c(13:15)]-br4[c(16:18)])/(br4[c(13:15)]+br4[c(16:18)])
RII.tian.tight.4<-(br4[c(19:21)]-br4[c(22:24)])/(br4[c(19:21)]+br4[c(22:24)])
f.4.<-c(RII.qilian.loose.4,RII.qilian.tight.4,RII.tian.loose.4,RII.tian.tight.4)

RII.qilian.loose.5<-(br5[c(1:3)]-br5[c(4:6)])/(br5[c(1:3)]+br5[c(4:6)])
RII.qilian.tight.5<-(br5[c(7:9)]-br5[c(10:12)])/(br5[c(7:9)]+br5[c(10:12)])
RII.tian.loose.5<-(br5[c(13:15)]-br5[c(16:18)])/(br5[c(13:15)]+br5[c(16:18)])
RII.tian.tight.5<-(br5[c(19:21)]-br5[c(22:24)])/(br5[c(19:21)]+br5[c(22:24)])
f.5.<-c(RII.qilian.loose.5,RII.qilian.tight.5,RII.tian.loose.5,RII.tian.tight.5)

RII.qilian.loose.6<-(br6[c(1:3)]-br6[c(4:6)])/(br6[c(1:3)]+br6[c(4:6)])
RII.qilian.tight.6<-(br6[c(7:9)]-br6[c(10:12)])/(br6[c(7:9)]+br6[c(10:12)])
RII.tian.loose.6<-(br6[c(13:15)]-br6[c(16:18)])/(br6[c(13:15)]+br6[c(16:18)])
RII.tian.tight.6<-(br6[c(19:21)]-br6[c(22:24)])/(br6[c(19:21)]+br6[c(22:24)])
f.6.<-c(RII.qilian.loose.6,RII.qilian.tight.6,RII.tian.loose.6,RII.tian.tight.6)

#
pheno<-rep(rep(c("A","B","A","B"),each=3),6)
site<-rep(rep(c("A","B"),each=6),6)
group<-rep(c("A","B","C","D","E","F"),each=12)
pheno.<-rep(c("A","B","A","B"),each=3)
site.<-rep(c("A","B"),each=6)

bact.ab<-c(f.1,f.2,f.3,f.4,f.5,f.6)
bact.ab.<-data.frame(site,pheno,group,bact.ab)
bact.ab.<-replace(bact.ab.,is.na(bact.ab.),0)
bact.abs<-data.frame(site.,pheno.,f.1,f.2,f.3,f.4,f.5,f.6)
bact.abs<-replace(bact.abs,is.na(bact.abs),0)

bact.r<-c(f.1.,f.2.,f.3.,f.4.,f.5.,f.6.)
bact.r.<-data.frame(site,pheno,group,bact.r)
bact.r.<-replace(bact.r.,is.na(bact.r.),0)
bact.rs<-data.frame(site.,pheno.,f.1.,f.2.,f.3.,f.4.,f.5.,f.6.)
bact.rs<-replace(bact.rs,is.na(bact.rs),0)

#ANOVA fungi
fit.1<-aovp(f.1~site.*pheno.,data=fungi.abs,perm=999)
summary.aov(fit.1)
fit.2<-aovp(f.2~site.*pheno.,data=fungi.abs,perm=999)
summary.aov(fit.2)
fit.3<-aovp(f.3~site.*pheno.,data=fungi.abs,perm=999)
summary.aov(fit.3)
fit.4<-aov(f.4~site.*pheno.,data=fungi.abs,perm=999)
summary.aov(fit.4)
fit.5<-aov(f.5~site.*pheno.,data=fungi.abs,perm=999)
summary.aov(fit.5)
fit.6<-aov(f.6~site.*pheno.,data=fungi.abs,perm=999)
summary.aov(fit.6)
fit.t<-aov(fungi.ab~site*group*pheno,data=fungi.ab.)
summary(fit.t)
cld(emmeans(fit.t,~site*group*pheno,type="response"),Letters=letters)
cld(emmeans(fit.t,~site*group,type="response"),Letters=letters)
cld(emmeans(fit.t,~group,type="response"),Letters=letters)

fit.1<-aov(f.1.~site.*pheno.,data=fungi.rs)
summary.aov(fit.1)
cld(emmeans(fit.1,~site.*pheno.,type="response"),Letters=letters)
fit.2<-aov(f.2.~site.*pheno.,data=fungi.rs)
summary.aov(fit.2)
cld(emmeans(fit.2,~site.*pheno.,type="response"),Letters=letters)
fit.3<-aov(f.3.~site.*pheno.,data=fungi.rs)
summary.aov(fit.3)
cld(emmeans(fit.3,~site.*pheno.,type="response"),Letters=letters)
fit.4<-aov(f.4.~site.*pheno.,data=fungi.rs)
summary.aov(fit.4)
cld(emmeans(fit.4,~site.*pheno.,type="response"),Letters=letters)
fit.5<-aov(f.5.~site.*pheno.,data=fungi.rs)
summary.aov(fit.5)
cld(emmeans(fit.5,~site.*pheno.,type="response"),Letters=letters)
fit.6<-aov(f.6.~site.*pheno.,data=fungi.rs)
summary.aov(fit.6)
cld(emmeans(fit.6,~site.*pheno.,type="response"),Letters=letters)
summary(aovp(fungi.r~site*group*pheno,data=fungi.r.,perm=999))
fanova<-aov(fungi.r~site*group*pheno,data=fungi.r.)
cld(emmeans(fanova,~site*group*pheno,type="response"),Letters=letters)

#ANOVA bacteria
fit.1<-aov(f.1~site.*pheno.,data=bact.abs)
summary.aov(fit.1)
cld(emmeans(fit.1,~site.*pheno.,type="response"),Letters=letters)
fit.2<-aov(f.2~site.*pheno.,data=bact.abs)
summary.aov(fit.2)
cld(emmeans(fit.2,~site.*pheno.,type="response"),Letters=letters)
fit.3<-aov(f.3~site.*pheno.,data=bact.abs)
summary.aov(fit.3)
cld(emmeans(fit.3,~site.*pheno.,type="response"),Letters=letters)
fit.4<-aov(f.4~site.*pheno.,data=bact.abs)
summary.aov(fit.4)
cld(emmeans(fit.4,~site.*pheno.,type="response"),Letters=letters)
fit.5<-aov(f.5~site.*pheno.,data=bact.abs)
summary.aov(fit.5)
cld(emmeans(fit.5,~site.*pheno.,type="response"),Letters=letters)
fit.6<-aov(f.6~site.*pheno.,data=bact.abs)
summary.aov(fit.6)
fit.t<-aov(bact.ab~site*group*pheno,data=bact.ab.)
summary(fit.t)
cld(emmeans(fit.t,~site*group*pheno,type="response"),Letters=letters)
cld(emmeans(fit.t,~site*group,type="response"),Letters=letters)
cld(emmeans(fit.t,~group,type="response"),Letters=letters)

fit.1<-aov(f.1.~site.*pheno.,data=bact.rs)
summary.aov(fit.1)
cld(emmeans(fit.1,~site.*pheno.,type="response"),Letters=letters)
fit.2<-aov(f.2.~site.*pheno.,data=bact.rs)
summary.aov(fit.2)
cld(emmeans(fit.2,~site.*pheno.,type="response"),Letters=letters)
fit.3<-aov(f.3.~site.*pheno.,data=bact.rs)
summary.aov(fit.3)
cld(emmeans(fit.3,~site.*pheno.,type="response"),Letters=letters)
fit.4<-aov(f.4.~site.*pheno.,data=bact.rs)
summary.aov(fit.4)
cld(emmeans(fit.4,~site.*pheno.,type="response"),Letters=letters)
fit.5<-aov(f.5.~site.*pheno.,data=bact.rs)
summary.aov(fit.5)
cld(emmeans(fit.5,~site.*pheno.,type="response"),Letters=letters)
fit.6<-aov(f.6.~site.*pheno.,data=bact.rs)
summary.aov(fit.6)
cld(emmeans(fit.6,~site.*pheno.,type="response"),Letters=letters)
summary(aovp(bact.r~site*group*pheno,data=bact.r.,perm=999))
banova<-aov(bact.r~site*group*pheno,data=bact.r.)
cld(emmeans(banova,~site*group*pheno,type="response"),Letters=letters)

#plot
final.fungi
final.bact

final.plot.fungi<-ddply(fungi.ab.,.(group,site,pheno),summarize,mean1=round(mean(fungi.ab),2),sd1=round(sd(fungi.ab)/sqrt(3),2))
final.plot.fungi2<-ddply(fungi.r.,.(group,site,pheno),summarize,mean1=round(mean(fungi.r),2),sd1=round(sd(fungi.r)/sqrt(3),2))

final.plot.bact<-ddply(bact.ab.,.(group,site,pheno),summarize,mean1=round(mean(bact.ab),2),sd1=round(sd(bact.ab)/sqrt(3),2))
final.plot.bact2<-ddply(bact.r.,.(group,site,pheno),summarize,mean1=round(mean(bact.r),2),sd1=round(sd(bact.r)/sqrt(3),2))

final.plot.fungi$group<-c(1.1,1.9,3.1,3.9,6.1,6.9,8.1,8.9,11.1,11.9,13.1,13.9,16.1,16.9,18.1,18.9,21.1,21.9,23.1,23.9,26.1,26.9,28.1,28.9)
final.plot.fungi2$group<-c(1.1,1.9,3.1,3.9,6.1,6.9,8.1,8.9,11.1,11.9,13.1,13.9,16.1,16.9,18.1,18.9,21.1,21.9,23.1,23.9,26.1,26.9,28.1,28.9)
final.plot.bact$group<-c(1.1,1.9,3.1,3.9,6.1,6.9,8.1,8.9,11.1,11.9,13.1,13.9,16.1,16.9,18.1,18.9,21.1,21.9,23.1,23.9,26.1,26.9,28.1,28.9)
final.plot.bact2$group<-c(1.1,1.9,3.1,3.9,6.1,6.9,8.1,8.9,11.1,11.9,13.1,13.9,16.1,16.9,18.1,18.9,21.1,21.9,23.1,23.9,26.1,26.9,28.1,28.9)



a<-3
while(a<=8){
  t.1<-t.test(fungi.rs[c(1:3),a])
  t.2<-t.test(fungi.rs[c(4:6),a])
  t.3<-t.test(fungi.rs[c(7:9),a])
  t.4<-t.test(fungi.rs[c(10:12),a])
  show(t.1)
  show(t.2)
  show(t.3)
  show(t.4)
  a<-a+1
}

a<-3
while(a<=8){
  t.1<-t.test(fungi.rs[c(1:3),a])
  t.2<-t.test(fungi.rs[c(4:6),a])
  t.3<-t.test(fungi.rs[c(7:9),a])
  t.4<-t.test(fungi.rs[c(10:12),a])
  show(t.1)
  show(t.2)
  show(t.3)
  show(t.4)
  a<-a+1
}

a<-3
while(a<=8){
  t.1<-t.test(bact.rs[c(1:3),a])
  t.2<-t.test(bact.rs[c(4:6),a])
  t.3<-t.test(bact.rs[c(7:9),a])
  t.4<-t.test(bact.rs[c(10:12),a])
  show(t.1)
  show(t.2)
  show(t.3)
  show(t.4)
  a<-a+1
}

text.fungi<-c("***","","","","","","","","***","***","","","","","","","","","","***","","**","**","")
text.fungi.<-c("","***","","","","***","","","","","","","","","(*)","*","***","*","***","","***","","","")
text.fungi2<-c("","","","","","","","","a","b","","a","b","","","","","","","b","","b","b","b")
text.fungi2.<-c("","","","","","","","","","","a","","","b","a","ab","a","a","a","","a","","","")

final.plot.fungi.<-data.frame(text.fungi,text.fungi.,final.plot.fungi)
final.plot.fungi2.<-data.frame(text.fungi,text.fungi.,final.plot.fungi2)


text.bact<-c("*","(*)","","","","","","","","","","","","","","(*)","","","","","","","","")
text.bact.<-c("","","","","","*","","","","","","","","","","","","","","","(*)","","","")
text.bact2<-c("ab","b","","","","","","","ab","ab","","b","","","","","","","","","","b","","")
text.bact2.<-c("","","a","ab","","","","","","","a","","","","","","","","","","a","","ab","ab")

final.plot.bact.<-data.frame(text.bact,text.bact.,final.plot.bact)
final.plot.bact2.<-data.frame(text.bact,text.bact.,final.plot.bact2)

bact.group<-ggplot() +
  geom_bar(data=final.plot.bact.,aes(x=group,y=mean1,fill=factor(pheno),color=pheno),stat="identity",position="dodge",colour="black",width = 0.7)+ 
  geom_errorbar(data=final.plot.bact.,aes(x=group,ymin=mean1-sd1,ymax=mean1+sd1,fill=pheno),stat="identity",position=position_dodge(width = 0.75),width=0.1)+
  geom_text(data=final.plot.bact.,aes(x=group,y=0.1 + mean1 + sd1,label=text.bact),color="black",fontface="bold",size=7.5)+
  geom_text(data=final.plot.bact.,aes(x=group,y=mean1-0.15 - sd1,label=text.bact.),color="black",fontface="bold",size=7.5)+
  geom_text(data=final.plot.bact.,aes(x=group,y=0.05 + mean1 + sd1,label=text.bact2),color="black",fontface="bold",size=5)+
  geom_text(data=final.plot.bact.,aes(x=group,y=mean1-0.05 - sd1,label=text.bact2.),color="black",fontface="bold",size=5)+
  geom_hline(yintercept=0)+ 
  geom_hline(yintercept=1.25)+ 
  geom_segment(aes(x=5,y=-1.25,xend=5,yend=1.5),linetype=1)+
  geom_segment(aes(x=10,y=-1.25,xend=10,yend=1.5),linetype=1)+
  geom_segment(aes(x=15,y=-1.25,xend=15,yend=1.5),linetype=1)+
  geom_segment(aes(x=20,y=-1.25,xend=20,yend=1.5),linetype=1)+
  geom_segment(aes(x=25,y=-1.25,xend=25,yend=1.5),linetype=1)+
  geom_segment(aes(x=30,y=-1.25,xend=30,yend=1.5),linetype=1)+
  geom_segment(aes(x=2.5,y=-1.25,xend=2.5,yend=1.25),linetype=2)+
  geom_segment(aes(x=7.5,y=-1.25,xend=7.5,yend=1.25),linetype=2)+
  geom_segment(aes(x=12.5,y=-1.25,xend=12.5,yend=1.25),linetype=2)+
  geom_segment(aes(x=17.5,y=-1.25,xend=17.5,yend=1.25),linetype=2)+
  geom_segment(aes(x=22.5,y=-1.25,xend=22.5,yend=1.25),linetype=2)+
  geom_segment(aes(x=27.5,y=-1.25,xend=27.5,yend=1.25),linetype=2)+
  geom_vline(xintercept=0)+
  scale_fill_manual(values = c("white","black"),
                    labels=c("loose phenotype","tight phenotype"))+
  scale_x_continuous(limits = c(0,30),
                     breaks=c(1.25,3.75,6.25,8.75,11.25,13.75,16.25,18.75,21.25,23.75,26.25,28.75),
                     labels=c("QL","TS","QL","TS","QL","TS","QL","TS","QL","TS","QL","TS"),
                     expand = c(0,0))+ 
  scale_y_continuous(limits=c(-1.25,1.5),  
                     breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),
                     labels=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),
                     expand = c(0,0))+
  labs(title=NULL, x=NULL, y="RIIabundance")+
  annotate("text",label = "Group 1", x = 2.5, y = 1.375,size=5)+
  annotate("text",label = "Group 2", x = 7.5, y = 1.375,size=5)+
  annotate("text",label = "Group 3", x = 12.5, y = 1.375,size=5)+
  annotate("text",label = "Group 4", x = 17.5, y = 1.375,size=5)+
  annotate("text",label = "Group 5", x = 22.5, y = 1.375,size=5)+
  annotate("text",label = "Group 6", x = 27.5, y = 1.375,size=5)+
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill="white", color="black"),
        panel.border=element_blank(),
        axis.line=element_line(color = "black"),
        legend.position="right",
        legend.key = element_blank(),
        legend.title= element_blank(),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=15,face = "bold"))
ggsave("E:/paper/bactgroupab.pdf",bact.group,width = 10,height = 6,dpi=300)


fungi.group<-ggplot() +
  geom_bar(data=final.plot.fungi.,aes(x=group,y=mean1,fill=factor(pheno),color=pheno),stat="identity",position="dodge",colour="black",width = 0.7)+ 
  geom_errorbar(data=final.plot.fungi.,aes(x=group,ymin=mean1-sd1,ymax=mean1+sd1,fill=pheno),stat="identity",position=position_dodge(width = 0.75),width=0.1)+
  geom_text(data=final.plot.fungi.,aes(x=group,y=0.1 + mean1 + sd1,label=text.fungi),color="black",fontface="bold",size=7.5)+
  geom_text(data=final.plot.fungi.,aes(x=group,y=mean1- sd1- 0.15,label=text.fungi.),color="black",fontface="bold",size=7.5)+
  geom_text(data=final.plot.fungi.,aes(x=group,y=0.05 + mean1 + sd1,label=text.fungi2),color="black",fontface="bold",size=5)+
  geom_text(data=final.plot.fungi.,aes(x=group,y=mean1-0.05 - sd1,label=text.fungi2.),color="black",fontface="bold",size=5)+
  geom_hline(yintercept=0)+ 
  geom_hline(yintercept=1.25)+ 
  geom_segment(aes(x=5,y=-1.25,xend=5,yend=1.5),linetype=1)+
  geom_segment(aes(x=10,y=-1.25,xend=10,yend=1.5),linetype=1)+
  geom_segment(aes(x=15,y=-1.25,xend=15,yend=1.5),linetype=1)+
  geom_segment(aes(x=20,y=-1.25,xend=20,yend=1.5),linetype=1)+
  geom_segment(aes(x=25,y=-1.25,xend=25,yend=1.5),linetype=1)+
  geom_segment(aes(x=30,y=-1.25,xend=30,yend=1.5),linetype=1)+
  geom_segment(aes(x=2.5,y=-1.25,xend=2.5,yend=1.25),linetype=2)+
  geom_segment(aes(x=7.5,y=-1.25,xend=7.5,yend=1.25),linetype=2)+
  geom_segment(aes(x=12.5,y=-1.25,xend=12.5,yend=1.25),linetype=2)+
  geom_segment(aes(x=17.5,y=-1.25,xend=17.5,yend=1.25),linetype=2)+
  geom_segment(aes(x=22.5,y=-1.25,xend=22.5,yend=1.25),linetype=2)+
  geom_segment(aes(x=27.5,y=-1.25,xend=27.5,yend=1.25),linetype=2)+
  geom_vline(xintercept=0)+
  scale_fill_manual(values = c("white","black"),
                    labels=c("loose phenotype","tight phenotype"))+
  scale_x_continuous(limits = c(0,30),
                     breaks=c(1.25,3.75,6.25,8.75,11.25,13.75,16.25,18.75,21.25,23.75,26.25,28.75),
                     labels=c("QL","TS","QL","TS","QL","TS","QL","TS","QL","TS","QL","TS"),
                     expand = c(0,0))+ 
  scale_y_continuous(limits=c(-1.25,1.5),  
                     breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),
                     labels=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),
                     expand = c(0,0))+
  labs(title=NULL, x=NULL, y="RIIabundance")+
  annotate("text",label = "Group 1", x = 2.5, y = 1.375,size=5)+
  annotate("text",label = "Group 2", x = 7.5, y = 1.375,size=5)+
  annotate("text",label = "Group 3", x = 12.5, y = 1.375,size=5)+
  annotate("text",label = "Group 4", x = 17.5, y = 1.375,size=5)+
  annotate("text",label = "Group 5", x = 22.5, y = 1.375,size=5)+
  annotate("text",label = "Group 6", x = 27.5, y = 1.375,size=5)+
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill="white", color="black"),
        panel.border=element_blank(),
        axis.line=element_line(color = "black"),
        legend.position="right",
        legend.key = element_blank(),
        legend.title= element_blank(),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=15,face = "bold"))
ggsave("E:/paper/fungigrouprich.pdf",fungi.group,width = 10,height = 6,dpi=300)
ggsave("E:/paper/fungigroupab.pdf",fungi.group,width = 10,height = 6,dpi=300)


###################################################
