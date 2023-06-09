library(xlsx)
library(vegan)
library(ggplot2)
library(data.table)

#richness
rich<-function(data){
  n<-rownames(data)
  rich<-c()
  for (i in 1:length(n)) {
    rich<-c(rich,length(which(data[i,]!=0)))
  }
  
  return(rich)
}

#data
data1.10<-read.xlsx("bacteriadata.xlsx",1)
data2.10<-read.xlsx("fungidata.xlsx",1)
data3.1<-read.xlsx("plantdata.xlsx",1)
data1<-data1.10[,-1]#bacteriadata
data2<-data2.10[,-1]#fungidata
data3<-data3.1[,-1]#plantdata

bd<-diversity(data1,index = "shannon")#bacteria diversity
fd<-diversity(data2,index = "shannon")#fungi diversity
pd<-diversity(data3,index = "shannon")#plant diversity

br<-rich(data1)#bacteria richness
fr<-rich(data2)#fungi richness
pr<-rich(data3)#plant richness

ba<-data1/rowSums(data1)
fa<-data2/rowSums(data2)
pa<-data3/rowSums(data3)

#plant RII
d1<-pd[c(1:40)]
d2<-pd[c(41:80)]
d3<-pd[c(81:120)]
d4<-pd[c(121:160)]
d5<-pd[c(161:200)]
d6<-pd[c(201:240)]
d7<-pd[c(241:280)]
d8<-pd[c(281:320)]

RII.QL<-(d1-d2)/(d1+d2)
RII.QT<-(d3-d4)/(d3+d4)
RII.TL<-(d5-d6)/(d5+d6)
RII.TT<-(d7-d8)/(d7+d8)

plant.community.RII<-rbind(RII.QL,RII.QT,RII.TL,RII.TT)
na.<-is.na(plant.community.RII)
plant.community.RII<-replace(plant.community.RII,na.,0)

d1<-pr[c(1:40)]
d2<-pr[c(41:80)]
d3<-pr[c(81:120)]
d4<-pr[c(121:160)]
d5<-pr[c(161:200)]
d6<-pr[c(201:240)]
d7<-pr[c(241:280)]
d8<-pr[c(281:320)]

RII.QL<-(d1-d2)/(d1+d2)
RII.QT<-(d3-d4)/(d3+d4)
RII.TL<-(d5-d6)/(d5+d6)
RII.TT<-(d7-d8)/(d7+d8)

plant.community.RII.R<-rbind(RII.QL,RII.QT,RII.TL,RII.TT)
na.<-is.na(plant.community.RII.R)
plant.community.RII.R<-replace(plant.community.RII.R,na.,0)

#bacteria RII
d1<-bd[c(1:3)]
d2<-bd[c(4:6)]
d3<-bd[c(7:9)]
d4<-bd[c(10:12)]
d5<-bd[c(13:15)]
d6<-bd[c(16:18)]
d7<-bd[c(19:21)]
d8<-bd[c(22:24)]

RII.QL<-(d1-d2)/(d1+d2)
RII.QT<-(d3-d4)/(d3+d4)
RII.TL<-(d5-d6)/(d5+d6)
RII.TT<-(d7-d8)/(d7+d8)

bact.community.RII<-rbind(RII.QL,RII.QT,RII.TL,RII.TT)
na.<-is.na(bact.community.RII)
bact.community.RII<-replace(bact.community.RII,na.,0)

d1<-br[c(1:3)]
d2<-br[c(4:6)]
d3<-br[c(7:9)]
d4<-br[c(10:12)]
d5<-br[c(13:15)]
d6<-br[c(16:18)]
d7<-br[c(19:21)]
d8<-br[c(22:24)]

RII.QL<-(d1-d2)/(d1+d2)
RII.QT<-(d3-d4)/(d3+d4)
RII.TL<-(d5-d6)/(d5+d6)
RII.TT<-(d7-d8)/(d7+d8)

bact.community.RII.R<-rbind(RII.QL,RII.QT,RII.TL,RII.TT)
na.<-is.na(bact.community.RII.R)
bact.community.RII.R<-replace(bact.community.RII.R,na.,0)

#fungi RII
d1<-fd[c(1:3)]
d2<-fd[c(4:6)]
d3<-fd[c(7:9)]
d4<-fd[c(10:12)]
d5<-fd[c(13:15)]
d6<-fd[c(16:18)]
d7<-fd[c(19:21)]
d8<-fd[c(22:24)]

RII.QL<-(d1-d2)/(d1+d2)
RII.QT<-(d3-d4)/(d3+d4)
RII.TL<-(d5-d6)/(d5+d6)
RII.TT<-(d7-d8)/(d7+d8)

fungi.community.RII<-rbind(RII.QL,RII.QT,RII.TL,RII.TT)
na.<-is.na(fungi.community.RII)
fungi.community.RII<-replace(fungi.community.RII,na.,0)

d1<-fr[c(1:3)]
d2<-fr[c(4:6)]
d3<-fr[c(7:9)]
d4<-fr[c(10:12)]
d5<-fr[c(13:15)]
d6<-fr[c(16:18)]
d7<-fr[c(19:21)]
d8<-fr[c(22:24)]

RII.QL<-(d1-d2)/(d1+d2)
RII.QT<-(d3-d4)/(d3+d4)
RII.TL<-(d5-d6)/(d5+d6)
RII.TT<-(d7-d8)/(d7+d8)

fungi.community.RII.R<-rbind(RII.QL,RII.QT,RII.TL,RII.TT)
na.<-is.na(fungi.community.RII.R)
fungi.community.RII.R<-replace(fungi.community.RII.R,na.,0)

#t-test
t.test(plant.community.RII[1,])
t.test(plant.community.RII[2,])
t.test(plant.community.RII[3,])
t.test(plant.community.RII[4,])

t.test(bact.community.RII[1,])
t.test(bact.community.RII[2,])
t.test(bact.community.RII[3,])
t.test(bact.community.RII[4,])

t.test(fungi.community.RII[1,])
t.test(fungi.community.RII[2,])
t.test(fungi.community.RII[3,])
t.test(fungi.community.RII[4,])


t.test(plant.community.RII.R[1,])
t.test(plant.community.RII.R[2,])
t.test(plant.community.RII.R[3,])
t.test(plant.community.RII.R[4,])

t.test(bact.community.RII.R[1,])
t.test(bact.community.RII.R[2,])
t.test(bact.community.RII.R[3,])
t.test(bact.community.RII.R[4,])

t.test(fungi.community.RII.R[1,])
t.test(fungi.community.RII.R[2,])
t.test(fungi.community.RII.R[3,])
t.test(fungi.community.RII.R[4,])

#ANOVA
tr<-c("A","A","B","B")
tr2<-c("A","B","A","B")

plant.community.RII.<-data.table(tr,tr2,plant.community.RII)
plant.community.RII.aov<-melt(plant.community.RII., id=1:2)
summary(aov(value~tr*tr2,data=plant.community.RII.aov))
plant.community.RII.R.<-data.table(tr,tr2,plant.community.RII.R)
plant.community.RII.R.aov<-melt(plant.community.RII.R., id=1:2)
summary(aov(value~tr*tr2,data=plant.community.RII.R.aov))

bact.community.RII.<-data.table(tr,tr2,bact.community.RII)
bact.community.RII.aov<-melt(bact.community.RII., id=1:2)
summary(aov(value~tr*tr2,data=bact.community.RII.aov))
bact.community.RII.R.<-data.table(tr,tr2,bact.community.RII.R)
bact.community.RII.R.aov<-melt(bact.community.RII.R., id=1:2)
summary(aov(value~tr*tr2,data=bact.community.RII.R.aov))

fungi.community.RII.<-data.table(tr,tr2,fungi.community.RII)
fungi.community.RII.aov<-melt(fungi.community.RII., id=1:2)
summary(aov(value~tr*tr2,data=fungi.community.RII.aov))
fungi.community.RII.R.<-data.table(tr,tr2,fungi.community.RII.R)
fungi.community.RII.R.aov<-melt(fungi.community.RII.R., id=1:2)
summary(aov(value~tr*tr2,data=fungi.community.RII.R.aov))

#plot
plant.plot.mean<-c()
plant.plot.sd<-c()
a<-1
while(a<=4){
  b<-mean(t(plant.community.RII)[,a])
  c<-sd(t(plant.community.RII)[,a])/sqrt(length(t(plant.community.RII)[,a])-1)
  plant.plot.mean<-c(plant.plot.mean,b)
  plant.plot.sd<-c(plant.plot.sd,c)
  a<-a+1
}

bact.plot.mean<-c()
bact.plot.sd<-c()
a<-1
while(a<=4){
  b<-mean(t(bact.community.RII)[,a])
  c<-sd(t(bact.community.RII)[,a])/sqrt(length(t(bact.community.RII)[,a])-1)
  bact.plot.mean<-c(bact.plot.mean,b)
  bact.plot.sd<-c(bact.plot.sd,c)
  a<-a+1
}

fungi.plot.mean<-c()
fungi.plot.sd<-c()
a<-1
while(a<=4){
  b<-mean(t(fungi.community.RII)[,a])
  c<-sd(t(fungi.community.RII)[,a])/sqrt(length(t(fungi.community.RII)[,a])-1)
  fungi.plot.mean<-c(fungi.plot.mean,b)
  fungi.plot.sd<-c(fungi.plot.sd,c)
  a<-a+1
}

mean<-c(plant.plot.mean,fungi.plot.mean,bact.plot.mean)
sd<-c(plant.plot.sd,fungi.plot.sd,bact.plot.sd)
group<-c(1.1,1.9,3.1,3.9,5.1,5.9,7.1,7.9,9.1,9.9,11.1,11.9)
type<-c(1,2,1,2,1,2,1,2,1,2,1,2)
site<-c(1,1,2,2,1,1,2,2,1,1,2,2)
text<-c("***","",".","","","","",".","","","","")
text2<-c("","","","","","","","","","","","")
plot.data<-data.frame(group,type,site,mean,sd,text,text2)

community<-ggplot()+ 
  geom_bar(data=plot.data,aes(x=group,y=mean,fill=factor(type)),stat="identity",position="stack",colour="black",width = 0.75)+ 
  geom_errorbar(data=plot.data,aes(x=group,ymin=mean-sd,ymax=mean+sd),position=position_dodge(width = 0.2),width=0.1)+
  geom_text(data=plot.data,aes(x=group,y=1.1*(mean + sd),label=text),color="black",fontface="bold",size=7.5)+
  geom_text(data=plot.data,aes(x=group,y=1.1 * (mean - sd),label=text2),color="black",fontface="bold",size=7.5)+
  geom_hline(yintercept=0)+ 
  geom_vline(xintercept=0)+
  geom_hline(yintercept = 0.4)+
  geom_vline(xintercept = 4.5,linetype=1)+
  geom_vline(xintercept = 8.5,linetype=1)+
  geom_segment(aes(x = 2.5,xend=2.5,y=-0.25,yend=0.4),linetype=2)+
  geom_segment(aes(x = 6.5,xend=6.5,y=-0.25,yend=0.4),linetype=2)+
  geom_segment(aes(x = 10.5,xend=10.5,y=-0.25,yend=0.4),linetype=2)+
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks=c(1.5,3.5,5.5,7.5,9.5,11.5),
                     labels=c("QL","TS","QL","TS","QL","TS"),
                     expand = c(0,0))+ 
  scale_y_continuous(limits=c(-0.25,0.55),  
                     breaks=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4),
                     labels=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4),
                     expand = c(0,0))+
  scale_fill_manual(values=c("white","black"),
                    labels=c("loose","tight"))+
  labs(title=NULL, x=NULL, y="RII of species richness")+
  annotate("text",label = "plant", x = 2.5, y = 0.475,size=5,color="gray30")+
  annotate("text",label = "fungi", x = 6.5, y = 0.475,size=5,color="gray30")+
  annotate("text",label = "bactreia", x = 10.5, y = 0.475,size=5,color="gray30")+
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill="white", color="black"),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        panel.ontop = FALSE,
        legend.position="right",
        legend.title=element_blank(),
        axis.ticks=element_blank(),
        legend.text=element_text(size=12,face = "bold"),
        #axis.text.x = element_text(angle=45, hjust=1),
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=15,face = "bold"))

ggsave("W:/论文发表/communityrich.pdf",community,width = 10,height = 6,dpi=300)


#plot richness
plant.plot.mean<-c()
plant.plot.sd<-c()
a<-1
while(a<=4){
  b<-mean(t(plant.community.RII.R)[,a])
  c<-sd(t(plant.community.RII.R)[,a])/sqrt(length(t(plant.community.RII.R)[,a]))
  plant.plot.mean<-c(plant.plot.mean,b)
  plant.plot.sd<-c(plant.plot.sd,c)
  a<-a+1
}

bact.plot.mean<-c()
bact.plot.sd<-c()
a<-1
while(a<=4){
  b<-mean(t(bact.community.RII.R)[,a])
  c<-sd(t(bact.community.RII.R)[,a])/sqrt(length(t(bact.community.RII.R)[,a]))
  bact.plot.mean<-c(bact.plot.mean,b)
  bact.plot.sd<-c(bact.plot.sd,c)
  a<-a+1
}

fungi.plot.mean<-c()
fungi.plot.sd<-c()
a<-1
while(a<=4){
  b<-mean(t(fungi.community.RII.R)[,a])
  c<-sd(t(fungi.community.RII.R)[,a])/sqrt(length(t(fungi.community.RII.R)[,a]))
  fungi.plot.mean<-c(fungi.plot.mean,b)
  fungi.plot.sd<-c(fungi.plot.sd,c)
  a<-a+1
}

mean<-c(plant.plot.mean,fungi.plot.mean,bact.plot.mean)
sd<-c(plant.plot.sd,fungi.plot.sd,bact.plot.sd)
group<-c(1.1,1.9,3.1,3.9,5.1,5.9,7.1,7.9,9.1,9.9,11.1,11.9)
type<-c(1,2,1,2,1,2,1,2,1,2,1,2)
site<-c(1,1,2,2,1,1,2,2,1,1,2,2)
text<-c("***","","(*)","","","","","(*)","","","","")
text2<-c("","","","","","","","","","","","")
plot.data<-data.frame(group,type,site,mean,sd,text,text2)

community<-ggplot()+ 
  geom_bar(data=plot.data,aes(x=group,y=mean,fill=factor(type)),stat="identity",position="stack",colour="black",width = 0.7)+ 
  geom_errorbar(data=plot.data,aes(x=group,ymin=mean-sd,ymax=mean+sd),position=position_dodge(width = 0.2),width=0.1)+
  geom_text(data=plot.data,aes(x=group,y=1.1*(mean+sd),label=text),color="black",fontface="bold",size=7.5)+
  geom_text(data=plot.data,aes(x=group,y=1.1*(mean-sd),label=text2),color="black",fontface="bold",size=7.5)+
  geom_hline(yintercept=0)+ 
  geom_vline(xintercept=0)+
  geom_hline(yintercept = 0.4)+
  geom_vline(xintercept = 4.5,linetype=1)+
  geom_vline(xintercept = 8.5,linetype=1)+
  geom_segment(aes(x = 2.5,xend=2.5,y=-0.25,yend=0.4),linetype=2)+
  geom_segment(aes(x = 6.5,xend=6.5,y=-0.25,yend=0.4),linetype=2)+
  geom_segment(aes(x = 10.5,xend=10.5,y=-0.25,yend=0.4),linetype=2)+
  scale_x_continuous(limits = c(0.5,12.5),
                     breaks=c(1.5,3.5,5.5,7.5,9.5,11.5),
                     labels=c("QL","TS","QL","TS","QL","TS"),
                     expand = c(0,0))+ 
  scale_y_continuous(limits=c(-0.25,0.55),  
                     breaks=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4),
                     labels=c(-0.2,-0.1,0,0.1,0.2,0.3,0.4),
                     expand = c(0,0))+
  scale_fill_manual(values=c("white","black"),
                    labels=c("loose phenotype","tight phenotype"))+
  labs(title=NULL, x=NULL, y="RIIrichness")+
  annotate("text",label = "plant", x = 2.5, y = 0.475,size=5,color="gray30")+
  annotate("text",label = "fungi", x = 6.5, y = 0.475,size=5,color="gray30")+
  annotate("text",label = "bactreia", x = 10.5, y = 0.475,size=5,color="gray30")+
  theme(panel.grid=element_blank(),
        panel.background=element_rect(fill="white", color="black"),
        panel.border = element_rect(linetype = "dashed", fill = NA),
        panel.ontop = FALSE,
        legend.position="right",
        legend.title=element_blank(),
        axis.ticks=element_blank(),
        legend.text=element_text(size=12,face = "bold"),
        #axis.text.x = element_text(angle=45, hjust=1),
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=15,face = "bold"))

ggsave("E:/paper/communityrich.pdf",community,width = 10,height = 6,dpi=300)
