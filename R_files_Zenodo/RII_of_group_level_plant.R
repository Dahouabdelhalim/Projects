library(plyr)
library(ggplot2)
library(lmPerm)
library(multcomp)
library(emmeans)
data3
names(data3)<-c(1:length(data3))

#richness
rich<-function(data){
  n<-rownames(data)
  rich<-c()
  for (i in 1:length(n)) {
    rich<-c(rich,length(which(data[i,]!=0)))
  }
  
  return(rich)
}


##removing rare species
c<-c()
a<-1
while(a<=length(data3)){
  alpha<-length(which(data3[,a]==0))
  if(alpha<=300){
  c<-c(c,a)}
  a<-a+1
}

data3.<-data3[,c]
pn<-names(data3.1)[c+1]

#diversity
pd.<-diversity(data3.,"shannon")

#relative abundance
pa<-data3./rowSums(data3.)
na.<-is.na(pa)
pa<-replace(pa,na.,0)
pr<-rich(pa)


d1<-pa[c(1:40),]
d2<-pa[c(41:80),]
d3<-pa[c(81:120),]
d4<-pa[c(121:160),]
d5<-pa[c(161:200),]
d6<-pa[c(201:240),]
d7<-pa[c(241:280),]
d8<-pa[c(281:320),]

a<-1
c1<-c()
c2<-c()
c3<-c()
c4<-c()
c5<-c()
c6<-c()
c7<-c()
c8<-c()
while(a<=length(pa)){
  b1<-sum(d1[,a])/40
  b2<-sum(d2[,a])/40
  b3<-sum(d3[,a])/40
  b4<-sum(d4[,a])/40
  b5<-sum(d5[,a])/40
  b6<-sum(d6[,a])/40
  b7<-sum(d7[,a])/40
  b8<-sum(d8[,a])/40
  c1<-c(c1,b1)
  c2<-c(c2,b2)
  c3<-c(c3,b3)
  c4<-c(c4,b4)
  c5<-c(c5,b5)
  c6<-c(c6,b6)
  c7<-c(c7,b7)
  c8<-c(c8,b8)
  a<-a+1}


dt<-rbind(c1,c2,c3,c4,c5,c6,c7,c8)
dt<-replace(dt,is.na(dt),0)


RII.QL<-(dt[1,]-dt[2,])/(dt[1,]+dt[2,])
RII.QT<-(dt[3,]-dt[4,])/(dt[3,]+dt[4,])
RII.TL<-(dt[5,]-dt[6,])/(dt[5,]+dt[6,])
RII.TT<-(dt[7,]-dt[8,])/(dt[7,]+dt[8,])

RII.total<-rbind(RII.QL,RII.QT,RII.TL,RII.TT)
na.<-is.na(RII.total)
plant.RII<-replace(RII.total,na.,0)

RII.total.re<-data.frame(plant.RII)


# #PCA
# library("FactoMineR")
# library("factoextra")
# 
# 
# 
# fungi.RII.pca<-PCA(bact.clust,scale.unit = TRUE,ncp=2, graph = TRUE)
# eigenvalues<-fungi.RII.pca$eig
# fungi.RII.pca1<-fungi.RII.pca$ind$coord
# fungi.RII.pca2<-fungi.RII.pca$var$coord
# 
# color<-c(2,2,2,4,1,3,4,4,2,3,3,4,1,4)
# group.pca<-data.frame(color,fungi.RII.pca1)
# 
# ggplot()+
#   #geom_point(data=plot.rda.,aes(pca1m ,pca2m,color=factor(site)),size=3)+
#   geom_point(data= group.pca,aes(Dim.1, Dim.2,color=factor(color)),size=6)+
#   #scale_color_manual(values = c("red","blue"))+
#   scale_color_manual(values = c("red","blue","yellow","green"),
#                      labels=c("group1",
#                               "group2",
#                               "group3",
#                               "group4"))+
#   #guides(shape=guide_legend(ncol=2))+
#   #geom_errorbar(aes(ymin=pca2m-pca2se,ymax=pca2m+pca2se),
#                 #width=0.01)+
#   #geom_errorbarh(aes(y=pca2m,xmin=pca1m-pca1se,xmax=pca1m+pca1se),
#                  #height=0.01)+
#   geom_hline(yintercept=0,size=.2)+
#   geom_vline(xintercept=0,size=.2)+
#   #coord_cartesian(xlim = c(-0.3,0.3),ylim = c(-0.2,0.3))+
#   labs( x="PC1" , y="PC2"  )+
#   theme_bw()+
#   theme(panel.grid=element_blank(),
#         plot.title = element_blank(),
#         legend.position="none",
#         legend.key = element_blank(),
#         legend.title = element_blank(),
#         legend.text=element_text(size=12),
#         axis.title=element_text(size=20,face="bold"),
#         axis.text=element_text(size=15,face = "bold"))
# 
# 


#clust
clust<-t(RII.total.re)
row.names(clust)<-names(data3.)
dis<-dist(clust, method = "euclidean") 
fit.clust<-hclust(dis, method="ward.D2")
plot(fit.clust)
rect.hclust(fit.clust, k=4, border="red")
groups <- cutree(fit.clust, k=4)
groups<-as.data.frame(groups)
order(groups)
cresult<-data.frame(groups,clust)




#
names.g<-row.names(cresult)
g1<-names.g[which(cresult$groups==1)]
g2<-names.g[which(cresult$groups==2)]
g3<-names.g[which(cresult$groups==3)]
g4<-names.g[which(cresult$groups==4)]

g1<-as.numeric(g1)
g2<-as.numeric(g2)
g3<-as.numeric(g3)
g4<-as.numeric(g4)

plant1<-data3[,g1]
plant1.<-data.frame(rowSums(plant1))
plant2<-data3[,g2]
plant2.<-data.frame(rowSums(plant2))
plant3<-data3[,g3]
plant3.<-data.frame(rowSums(plant3))
plant4<-data3[,g4]
plant4.<-data.frame(rowSums(plant4))

pr1<-rich(plant1)
pr2<-rich(plant2)
pr3<-rich(plant3)
pr4<-rich(plant4)

#RII abundance
RII.QL<-(plant1.[c(1:40),]-plant1.[c(41:80),])/(plant1.[c(1:40),]+plant1.[c(41:80),])
RII.QT<-(plant1.[c(81:120),]-plant1.[c(121:160),])/(plant1.[c(81:120),]+plant1.[c(121:160),])
RII.TL<-(plant1.[c(161:200),]-plant1.[c(201:240),])/(plant1.[c(161:200),]+plant1.[c(201:240),])
RII.TT<-(plant1.[c(241:280),]-plant1.[c(281:320),])/(plant1.[c(241:280),]+plant1.[c(281:320),])

RII.QL<-replace(RII.QL,is.na(RII.QL),0)
RII.QT<-replace(RII.QT,is.na(RII.QT),0)
RII.TL<-replace(RII.TL,is.na(RII.TL),0)
RII.TT<-replace(RII.TT,is.na(RII.TT),0)

f.1<-c(RII.QL,RII.QT,RII.TL,RII.TT)

RII.QL<-(plant2.[c(1:40),]-plant2.[c(41:80),])/(plant2.[c(1:40),]+plant2.[c(41:80),])
RII.QT<-(plant2.[c(81:120),]-plant2.[c(121:160),])/(plant2.[c(81:120),]+plant2.[c(121:160),])
RII.TL<-(plant2.[c(161:200),]-plant2.[c(201:240),])/(plant2.[c(161:200),]+plant2.[c(201:240),])
RII.TT<-(plant2.[c(241:280),]-plant2.[c(281:320),])/(plant2.[c(241:280),]+plant2.[c(281:320),])

RII.QL<-replace(RII.QL,is.na(RII.QL),0)
RII.QT<-replace(RII.QT,is.na(RII.QT),0)
RII.TL<-replace(RII.TL,is.na(RII.TL),0)
RII.TT<-replace(RII.TT,is.na(RII.TT),0)

f.2<-c(RII.QL,RII.QT,RII.TL,RII.TT)

RII.QL<-(plant3.[c(1:40),]-plant3.[c(41:80),])/(plant3.[c(1:40),]+plant3.[c(41:80),])
RII.QT<-(plant3.[c(81:120),]-plant3.[c(121:160),])/(plant3.[c(81:120),]+plant3.[c(121:160),])
RII.TL<-(plant3.[c(161:200),]-plant3.[c(201:240),])/(plant3.[c(161:200),]+plant3.[c(201:240),])
RII.TT<-(plant3.[c(241:280),]-plant3.[c(281:320),])/(plant3.[c(241:280),]+plant3.[c(281:320),])

RII.QL<-replace(RII.QL,is.na(RII.QL),0)
RII.QT<-replace(RII.QT,is.na(RII.QT),0)
RII.TL<-replace(RII.TL,is.na(RII.TL),0)
RII.TT<-replace(RII.TT,is.na(RII.TT),0)

f.3<-c(RII.QL,RII.QT,RII.TL,RII.TT)


RII.QL<-(plant4.[c(1:40),]-plant4.[c(41:80),])/(plant4.[c(1:40),]+plant4.[c(41:80),])
RII.QT<-(plant4.[c(81:120),]-plant4.[c(121:160),])/(plant4.[c(81:120),]+plant4.[c(121:160),])
RII.TL<-(plant4.[c(161:200),]-plant4.[c(201:240),])/(plant4.[c(161:200),]+plant4.[c(201:240),])
RII.TT<-(plant4.[c(241:280),]-plant4.[c(281:320),])/(plant4.[c(241:280),]+plant4.[c(281:320),])

RII.QL<-replace(RII.QL,is.na(RII.QL),0)
RII.QT<-replace(RII.QT,is.na(RII.QT),0)
RII.TL<-replace(RII.TL,is.na(RII.TL),0)
RII.TT<-replace(RII.TT,is.na(RII.TT),0)

f.4<-c(RII.QL,RII.QT,RII.TL,RII.TT)

#RII richness
RII.QL<-(pr1[c(1:40)]-pr1[c(41:80)])/(pr1[c(1:40)]+pr1[c(41:80)])
RII.QT<-(pr1[c(81:120)]-pr1[c(121:160)])/(pr1[c(81:120)]+pr1[c(121:160)])
RII.TL<-(pr1[c(161:200)]-pr1[c(201:240)])/(pr1[c(161:200)]+pr1[c(201:240)])
RII.TT<-(pr1[c(241:280)]-pr1[c(281:320)])/(pr1[c(241:280)]+pr1[c(281:320)])

RII.QL<-replace(RII.QL,is.na(RII.QL),0)
RII.QT<-replace(RII.QT,is.na(RII.QT),0)
RII.TL<-replace(RII.TL,is.na(RII.TL),0)
RII.TT<-replace(RII.TT,is.na(RII.TT),0)

f.1.<-c(RII.QL,RII.QT,RII.TL,RII.TT)

RII.QL<-(pr2[c(1:40)]-pr2[c(41:80)])/(pr2[c(1:40)]+pr2[c(41:80)])
RII.QT<-(pr2[c(81:120)]-pr2[c(121:160)])/(pr2[c(81:120)]+pr2[c(121:160)])
RII.TL<-(pr2[c(161:200)]-pr2[c(201:240)])/(pr2[c(161:200)]+pr2[c(201:240)])
RII.TT<-(pr2[c(241:280)]-pr2[c(281:320)])/(pr2[c(241:280)]+pr2[c(281:320)])

RII.QL<-replace(RII.QL,is.na(RII.QL),0)
RII.QT<-replace(RII.QT,is.na(RII.QT),0)
RII.TL<-replace(RII.TL,is.na(RII.TL),0)
RII.TT<-replace(RII.TT,is.na(RII.TT),0)

f.2.<-c(RII.QL,RII.QT,RII.TL,RII.TT)

RII.QL<-(pr3[c(1:40)]-pr3[c(41:80)])/(pr3[c(1:40)]+pr3[c(41:80)])
RII.QT<-(pr3[c(81:120)]-pr3[c(121:160)])/(pr3[c(81:120)]+pr3[c(121:160)])
RII.TL<-(pr3[c(161:200)]-pr3[c(201:240)])/(pr3[c(161:200)]+pr3[c(201:240)])
RII.TT<-(pr3[c(241:280)]-pr3[c(281:320)])/(pr3[c(241:280)]+pr3[c(281:320)])

RII.QL<-replace(RII.QL,is.na(RII.QL),0)
RII.QT<-replace(RII.QT,is.na(RII.QT),0)
RII.TL<-replace(RII.TL,is.na(RII.TL),0)
RII.TT<-replace(RII.TT,is.na(RII.TT),0)

f.3.<-c(RII.QL,RII.QT,RII.TL,RII.TT)


RII.QL<-(pr4[c(1:40)]-pr4[c(41:80)])/(pr4[c(1:40)]+pr4[c(41:80)])
RII.QT<-(pr4[c(81:120)]-pr4[c(121:160)])/(pr4[c(81:120)]+pr4[c(121:160)])
RII.TL<-(pr4[c(161:200)]-pr4[c(201:240)])/(pr4[c(161:200)]+pr4[c(201:240)])
RII.TT<-(pr4[c(241:280)]-pr4[c(281:320)])/(pr4[c(241:280)]+pr4[c(281:320)])

RII.QL<-replace(RII.QL,is.na(RII.QL),0)
RII.QT<-replace(RII.QT,is.na(RII.QT),0)
RII.TL<-replace(RII.TL,is.na(RII.TL),0)
RII.TT<-replace(RII.TT,is.na(RII.TT),0)

f.4.<-c(RII.QL,RII.QT,RII.TL,RII.TT)

#ANOVA
pheno<-rep(rep(c("A","B","A","B"),each=40),4)
site<-rep(rep(c("A","B"),each=80),4)
group<-rep(c("A","B","C","D"),each=160)
pheno.<-rep(c("A","B","A","B"),each=40)
site.<-rep(c("A","B"),each=80)

plant.ab<-c(f.1,f.2,f.3,f.4)
plant.ab.<-data.frame(site,pheno,group,plant.ab)
plant.ab.<-replace(plant.ab.,is.na(plant.ab.),0)
plant.abs<-data.frame(site.,pheno.,f.1,f.2,f.3,f.4)
plant.abs<-replace(plant.abs,is.na(plant.abs),0)

plant.r<-c(f.1.,f.2.,f.3.,f.4.)
plant.r.<-data.frame(site,pheno,group,plant.r)
plant.r.<-replace(plant.r.,is.na(plant.r.),0)
plant.rs<-data.frame(site.,pheno.,f.1.,f.2.,f.3.,f.4.)
plant.rs<-replace(plant.rs,is.na(plant.rs),0)

#ANOVA test
fit.1<-aovp(f.1~site.*pheno.,data=plant.abs,perm=999)
summary.aov(fit.1)
fit.2<-aovp(f.2~site.*pheno.,data=plant.abs,perm=999)
summary.aov(fit.2)
fit.3<-aovp(f.3~site.*pheno.,data=plant.abs,perm=999)
summary.aov(fit.3)
fit.4<-aov(f.4~site.*pheno.,data=plant.abs,perm=999)
summary.aov(fit.4)
fit.t<-aov(plant.ab~site*group*pheno,data=plant.ab.)
summary(fit.t)
cld(emmeans(fit.t,~site*group*pheno,type="response"),Letters=letters)
cld(emmeans(fit.t,~site*group,type="response"),Letters=letters)
cld(emmeans(fit.t,~group,type="response"),Letters=letters)

fit.1<-aov(f.1.~site.*pheno.,data=plant.rs)
summary.aov(fit.1)
cld(emmeans(fit.1,~site.*pheno.,type="response"),Letters=letters)
fit.2<-aov(f.2.~site.*pheno.,data=plant.rs)
summary.aov(fit.2)
cld(emmeans(fit.2,~site.*pheno.,type="response"),Letters=letters)
fit.3<-aov(f.3.~site.*pheno.,data=plant.rs)
summary.aov(fit.3)
cld(emmeans(fit.3,~site.*pheno.,type="response"),Letters=letters)
fit.4<-aov(f.4.~site.*pheno.,data=plant.rs)
summary.aov(fit.4)
cld(emmeans(fit.4,~site.*pheno.,type="response"),Letters=letters)
summary(aovp(plant.r~site*group*pheno,data=plant.r.,perm=999))
summary(aov(plant.r~site*group*pheno,data=plant.r.))
panova<-aov(plant.r~site*group*pheno,data=plant.r.)
cld(emmeans(panova,~site*group*pheno,type="response"),Letters=letters)

#plot
final.plot.plant<-ddply(plant.ab.,.(group,site,pheno),summarize,mean1=round(mean(plant.ab),2),sd1=round(sd(plant.ab)/sqrt(40),2))
final.plot.plant2<-ddply(plant.r.,.(group,site,pheno),summarize,mean1=round(mean(plant.r),2),sd1=round(sd(plant.r)/sqrt(40),2))

final.plot.plant$group<-c(1.1,1.9,3.1,3.9,6.1,6.9,8.1,8.9,11.1,11.9,13.1,13.9,16.1,16.9,18.1,18.9)
final.plot.plant2$group<-c(1.1,1.9,3.1,3.9,6.1,6.9,8.1,8.9,11.1,11.9,13.1,13.9,16.1,16.9,18.1,18.9)


a<-3
while(a<=6){
  t.1<-t.test(plant.abs[c(1:40),a])
  t.2<-t.test(plant.abs[c(41:80),a])
  t.3<-t.test(plant.abs[c(81:120),a])
  t.4<-t.test(plant.abs[c(121:160),a])
  show(t.1)
  show(t.2)
  show(t.3)
  show(t.4)
  a<-a+1
}

a<-3
while(a<=6){
  t.1<-t.test(plant.rs[c(1:40),a])
  t.2<-t.test(plant.rs[c(41:80),a])
  t.3<-t.test(plant.rs[c(81:120),a])
  t.4<-t.test(plant.rs[c(121:160),a])
  show(t.1)
  show(t.2)
  show(t.3)
  show(t.4)
  a<-a+1
}


text.plant<-c("","","***","","","","","","","","*","","***","","","")
text.plant.<-c("","*","","","","","","","","","","","","","","")
text.plant2<-c("","","","","","","","","","","","","","","","")
text.plant2.<-c("","","","","","","","","","","","","","","","")
final.plot.plant.<-data.frame(text.plant,text.plant.,final.plot.plant)
final.plot.plant2.<-data.frame(text.plant,text.plant.,final.plot.plant2)

plant.group<-ggplot() +
  geom_bar(data=final.plot.plant.,aes(x=group,y=mean1,fill=factor(pheno)),stat="identity",position="dodge",colour="black",width = 0.7)+ 
  geom_errorbar(data=final.plot.plant.,aes(x=group,ymin=mean1-sd1,ymax=mean1+sd1,fill=factor(pheno)),stat="identity",position=position_dodge(width = 0.75),width=0.1)+
  geom_text(data=final.plot.plant.,aes(x=group,y=0.1 + mean1 + sd1,label=text.plant),color="black",fontface="bold",size=7.5)+
  geom_text(data=final.plot.plant.,aes(x=group,y=mean1-0.15 - sd1,label=text.plant.),color="black",fontface="bold",size=7.5)+
  geom_text(data=final.plot.plant.,aes(x=group,y=0.05 + mean1 + sd1,label=text.plant2),color="black",fontface="bold",size=5)+
  geom_text(data=final.plot.plant.,aes(x=group,y=mean1-0.05 - sd1,label=text.plant2.),color="black",fontface="bold",size=5)+
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
  geom_segment(aes(x=17.5,y=-0.75,xend=17.5,yend=1),linetype=2)+
  scale_x_continuous(limits = c(0,20),
                     breaks=c(1.25,3.75,6.25,8.75,11.25,13.75,16.25,18.75),
                     labels=c("QL","TS","QL","TS","QL","TS","QL","TS"),
                     expand = c(0,0))+ 
  scale_y_continuous(limits=c(-1.25,1.5),  
                     breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),
                     labels=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),
                     expand = c(0,0))+
  scale_fill_manual(values=c("white","black"),
                    labels=c("loose phenotype","tight phenotype"))+
  labs(title=NULL, x=NULL, y="RIIrichness")+
  annotate("text",label = "Group 1", x = 2.5, y = 1.375,size=5)+
  annotate("text",label = "Group 2", x = 7.5, y = 1.375,size=5)+
  annotate("text",label = "Group 3", x = 12.5, y = 1.375,size=5)+
  annotate("text",label = "Group 4", x = 17.5, y = 1.375,size=5)+
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

ggsave("E:/paper/plantgrouprich.pdf",plant.group,width = 10,height = 6,dpi=300)
ggsave("E:/paper/plantgroupab.pdf",plant.group,width = 10,height = 6,dpi=300)


####################################################################
####################################################################

