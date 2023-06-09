library(NST)
library(vegan)
library(lmPerm)
library(data.table)
library(emmeans)
library(multcomp)


group1<-data.frame(substr(data1.10[,1],1,4))
group2<-data.frame(substr(data2.10[,1],1,4))
group3<-data.frame(substr(data3.1[,1],1,4))
g<-paste(t1,t2,t3)
g2<-substr(data3.1[,1],1,4)
set.seed(1234)

rownames(data1.)<-rownames(group1)
rownames(data2.)<-rownames(group2)
rownames(data3.)<-rownames(group3)

set.seed(1234)
#bacteria
tnst <- tNST(comm = data1.,group=group1, dist.method = 'bray', null.model = 'PF', between.group=F,
             rand = 999, nworker = 1)

nst_bact_t <- tnst$index.pair.grp
nst_bact<-nst_bact_t$MST.ij.bray
t1<-substr(nst_bact_t$group,1,2)
t2<-substr(nst_bact_t$group,3,3)
t3<-substr(nst_bact_t$group,4,4)
summary(aovp(nst_bact_t$MST.ij.bray~t1*t2*t3,data=nst_bact_t))

d1<-nst_bact[c(1:3)]
d2<-nst_bact[c(4:6)]
d3<-nst_bact[c(7:9)]
d4<-nst_bact[c(10:12)]
d5<-nst_bact[c(13:15)]
d6<-nst_bact[c(16:18)]
d7<-nst_bact[c(19:21)]
d8<-nst_bact[c(22:24)]

tr<-c(1,2,4,5,7,8,10,11)
t<-c(1,2,1,2,1,2,1,2)
mps<-c(mean(d1),mean(d2),mean(d3),mean(d4),mean(d5),mean(d6),mean(d7),mean(d8))
sdps<-c(sd(d1)/sqrt(length(d1)-1),sd(d2)/sqrt(length(d2)-1),sd(d3)/sqrt(length(d3)-1),sd(d4)/sqrt(length(d4)-1),
        sd(d5)/sqrt(length(d5)-1),sd(d6)/sqrt(length(d6)-1),sd(d7)/sqrt(length(d7)-1),sd(d8)/sqrt(length(d8)-1))
plot.data.ps<-data.frame(tr,mps,sdps)

bs<-ggplot()+ 
  geom_bar(data=plot.data.ps,aes(x=tr,y=mps,fill=factor(t)),stat="identity",position="stack",colour="black",width = 0.75)+ 
  geom_errorbar(data=plot.data.ps,aes(x=tr,ymin=mps-sdps,ymax=mps+sdps),position=position_dodge(width = 0.2),width=0.1)+
  # geom_text(data=plot.data,aes(x=group,y=1.1*(mean + sd),label=text),color="black",fontface="bold",size=7.5)+
  # geom_text(data=plot.data,aes(x=group,y=1.1 * (mean - sd),label=text2),color="black",fontface="bold",size=7.5)+
  geom_hline(yintercept=0)+ 
  geom_vline(xintercept=0)+
  geom_hline(yintercept = 0.5,linetype=2,color="red")+
  geom_vline(xintercept = 6,linetype=1)+
  geom_vline(xintercept = 12,linetype=1)+
  geom_hline(yintercept = 1,linetype=1)+
  geom_hline(yintercept = 1.2,linetype=1)+
  geom_segment(aes(x = 3,xend=3,y=0,yend=1),linetype=2)+
  geom_segment(aes(x = 9,xend=9,y=0,yend=1),linetype=2)+
  scale_x_continuous(limits = c(0,12),
                     breaks=c(1.5,4.5,7.5,11.5),
                     labels=c("L","T","L","T"),
                     expand = c(0,0))+ 
  scale_y_continuous(limits=c(0,1.2),  
                     breaks=c(0.25,0.5,0.75,1),
                     labels=c(0.25,0.5,0.75,1),
                     expand = c(0,0))+
  scale_fill_manual(values=c("white","black"),
                    labels=c("In","Out"))+
  labs(title=NULL, x=NULL, y="Stochasticity %")+
  annotate("text",label = "Qilian", x = 3, y = 1.1,size=5,color="gray30")+
  annotate("text",label = "Tianshan", x = 9, y = 1.1,size=5,color="gray30")+
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

ggsave("bs.pdf",bs,width = 10,height = 6,dpi=300)

b.RII.QL<-(d1-d2)/(d1+d2)
b.RII.QT<-(d3-d4)/(d3+d4)
b.RII.TL<-(d5-d6)/(d5+d6)
b.RII.TT<-(d7-d8)/(d7+d8)

b.RII.s<-c(b.RII.QL,b.RII.QT,b.RII.TL,b.RII.TT)
ts<-c(rep(1,length(b.RII.QL)),rep(1,length(b.RII.QT)),rep(2,length(b.RII.TL)),rep(2,length(b.RII.TT)))
tp<-c(rep(1,length(b.RII.QL)),rep(2,length(b.RII.QT)),rep(1,length(b.RII.TL)),rep(2,length(b.RII.TT)))
summary(aovp(b.RII.s~ts*tp,perm=999))
cld(emmeans(aovp(b.RII.s~ts*tp),~ts*tp),Letters=letters)

#fungi
tnst <- tNST(comm = data2.,group=group2, dist.method = 'bray', null.model = 'PF', between.group=F,
             rand = 999, nworker = 1)

nst_fungi_t <- tnst$index.pair.grp
nst_fungi<-nst_fungi_t$MST.ij.bray
t1<-substr(nst_fungi_t$group,1,2)
t2<-substr(nst_fungi_t$group,3,3)
t3<-substr(nst_fungi_t$group,4,4)
summary(aovp(nst_fungi_t$MST.ij.bray~t1*t2*t3,perm=999))

d1<-nst_fungi[c(1:3)]
d2<-nst_fungi[c(4:6)]
d3<-nst_fungi[c(7:9)]
d4<-nst_fungi[c(10:12)]
d5<-nst_fungi[c(13:15)]
d6<-nst_fungi[c(16:18)]
d7<-nst_fungi[c(19:21)]
d8<-nst_fungi[c(22:24)]

tr<-c(1,2,4,5,7,8,10,11)
t<-c(1,2,1,2,1,2,1,2)
text<-c("b","a","a","ab","a","ab","a","ab")
mps<-c(mean(d1),mean(d2),mean(d3),mean(d4),mean(d5),mean(d6),mean(d7),mean(d8))
sdps<-c(sd(d1)/sqrt(length(d1)-1),sd(d2)/sqrt(length(d2)-1),sd(d3)/sqrt(length(d3)-1),sd(d4)/sqrt(length(d4)-1),
        sd(d5)/sqrt(length(d5)-1),sd(d6)/sqrt(length(d6)-1),sd(d7)/sqrt(length(d7)-1),sd(d8)/sqrt(length(d8)-1))
plot.data.fs<-data.frame(tr,mps,sdps,text)

fs<-ggplot()+ 
  geom_bar(data=plot.data.fs,aes(x=tr,y=mps,fill=factor(t)),stat="identity",position="stack",colour="black",width = 0.75)+ 
  geom_errorbar(data=plot.data.fs,aes(x=tr,ymin=mps-sdps,ymax=mps+sdps),position=position_dodge(width = 0.2),width=0.1)+
  geom_text(data=plot.data.fs,aes(x=tr,y=0.025+(mps+sdps),label=text),color="black",fontface="bold",size=7.5)+
  geom_hline(yintercept=0)+ 
  geom_vline(xintercept=0)+
  geom_hline(yintercept = 0.5,linetype=2,color="red")+
  geom_vline(xintercept = 6,linetype=1)+
  geom_vline(xintercept = 12,linetype=1)+
  geom_hline(yintercept = 1,linetype=1)+
  geom_hline(yintercept = 1.2,linetype=1)+
  geom_segment(aes(x = 3,xend=3,y=0,yend=1),linetype=2)+
  geom_segment(aes(x = 9,xend=9,y=0,yend=1),linetype=2)+
  scale_x_continuous(limits = c(0,12),
                     breaks=c(1.5,4.5,7.5,11.5),
                     labels=c("L","T","L","T"),
                     expand = c(0,0))+ 
  scale_y_continuous(limits=c(0,1.2),  
                     breaks=c(0.25,0.5,0.75,1),
                     labels=c(0.25,0.5,0.75,1),
                     expand = c(0,0))+
  scale_fill_manual(values=c("white","black"),
                    labels=c("In","Out"))+
  labs(title=NULL, x=NULL, y="Stochasticity %")+
  annotate("text",label = "Qilian", x = 3, y = 1.1,size=5,color="gray30")+
  annotate("text",label = "Tianshan", x = 9, y = 1.1,size=5,color="gray30")+
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

ggsave("fs.pdf",fs,width = 10,height = 6,dpi=300)

f.RII.QL<-(d1-d2)/(d1+d2)
f.RII.QT<-(d3-d4)/(d3+d4)
f.RII.TL<-(d5-d6)/(d5+d6)
f.RII.TT<-(d7-d8)/(d7+d8)

f.RII.s<-c(f.RII.QL,f.RII.QT,f.RII.TL,f.RII.TT)
ts<-rep(c(1,2),each=6)
tp<-rep(c(1,2,1,2),each=3)
summary(aovp(f.RII.s~ts*tp,perm=999))
cld(emmeans(aovp(f.RII.s~ts*tp),~ts*tp),Letters=letters)

#plant
data3r<-data3[-which(rowSums(data3)==0),]
group3r<-group3[-which(rowSums(data3)==0)]

tnst <- tNST(comm = data3r,group=group3r, dist.method = 'bray', null.model = 'PF', between.group=F,
             rand = 999, nworker = 1)

nst_plant_t <- tnst$index.pair.grp
nst_plant<-nst_plant_t$MST.ij.bray
t1<-substr(nst_plant_t$group,1,2)
t2<-substr(nst_plant_t$group,3,3)
t3<-substr(nst_plant_t$group,4,4)
summary(aovp(nst_plant_t$MST.ij.bray~t1*t2*t3,perm=999))

d1<-nst_plant[which(nst_plant_t$group=="QLLI")]
d2<-nst_plant[which(nst_plant_t$group=="QLLO")]
d3<-nst_plant[which(nst_plant_t$group=="QLTI")]
d4<-nst_plant[which(nst_plant_t$group=="QLTO")]
d5<-nst_plant[which(nst_plant_t$group=="TSLI")]
d6<-nst_plant[which(nst_plant_t$group=="TSLO")]
d7<-nst_plant[which(nst_plant_t$group=="TSLI")]
d8<-nst_plant[which(nst_plant_t$group=="TSLO")]

tr<-c(1,2,4,5,7,8,10,11)
t<-c(1,2,1,2,1,2,1,2)
text<-c("f","d","e","c","b","a","ab","ab")
mps<-c(mean(d1),mean(d2),mean(d3),mean(d4),mean(d5),mean(d6),mean(d7),mean(d8))
sdps<-c(sd(d1)/sqrt(length(d1)-1),sd(d2)/sqrt(length(d2)-1),sd(d3)/sqrt(length(d3)-1),sd(d4)/sqrt(length(d4)-1),
        sd(d5)/sqrt(length(d5)-1),sd(d6)/sqrt(length(d6)-1),sd(d7)/sqrt(length(d7)-1),sd(d8)/sqrt(length(d8)-1))
plot.data.ps<-data.frame(tr,mps,sdps)

ps<-ggplot()+ 
  geom_bar(data=plot.data.ps,aes(x=tr,y=mps,fill=factor(t)),stat="identity",position="stack",colour="black",width = 0.75)+ 
  geom_errorbar(data=plot.data.ps,aes(x=tr,ymin=mps-sdps,ymax=mps+sdps),position=position_dodge(width = 0.2),width=0.1)+
  geom_text(data=plot.data.ps,aes(x=tr,y=0.025+(mps+sdps),label=text),color="black",fontface="bold",size=7.5)+
  geom_hline(yintercept=0)+ 
  geom_vline(xintercept=0)+
  geom_hline(yintercept = 0.5,linetype=2,color="red")+
  geom_vline(xintercept = 6,linetype=1)+
  geom_vline(xintercept = 12,linetype=1)+
  geom_hline(yintercept = 1,linetype=1)+
  geom_hline(yintercept = 1.2,linetype=1)+
  geom_segment(aes(x = 3,xend=3,y=0,yend=1),linetype=2)+
  geom_segment(aes(x = 9,xend=9,y=0,yend=1),linetype=2)+
  scale_x_continuous(limits = c(0,12),
                     breaks=c(1.5,4.5,7.5,11.5),
                     labels=c("L","T","L","T"),
                     expand = c(0,0))+ 
  scale_y_continuous(limits=c(0,1.2),  
                     breaks=c(0.25,0.5,0.75,1),
                     labels=c(0.25,0.5,0.75,1),
                     expand = c(0,0))+
  scale_fill_manual(values=c("white","black"),
                    labels=c("In","Out"))+
  labs(title=NULL, x=NULL, y="Stochasticity %")+
  annotate("text",label = "Qilian", x = 3, y = 1.1,size=5,color="gray30")+
  annotate("text",label = "Tianshan", x = 9, y = 1.1,size=5,color="gray30")+
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

ggsave("ps.pdf",ps,width = 10,height = 6,dpi=300)


p.RII.QL<-(d1-d2)/(d1+d2)
p.RII.QT<-(d3-d4)/(d3+d4)
p.RII.TL<-(d5-d6)/(d5+d6)
p.RII.TT<-(d7-d8)/(d7+d8)
p.RII.QL<-replace(p.RII.QL,is.na(p.RII.QL),0)
p.RII.QT<-replace(p.RII.QT,is.na(p.RII.QT),0)
p.RII.TL<-replace(p.RII.TL,is.na(p.RII.TL),0)
p.RII.TT<-replace(p.RII.TT,is.na(p.RII.TT),0)

p.RII.s<-c(p.RII.QL,p.RII.QT,p.RII.TL,p.RII.TT)
ts<-c(rep(1,length(p.RII.QL)),rep(1,length(p.RII.QT)),rep(2,length(p.RII.TL)),rep(2,length(p.RII.TT)))
tp<-c(rep(1,length(p.RII.QL)),rep(2,length(p.RII.QT)),rep(1,length(p.RII.TL)),rep(2,length(p.RII.TT)))
summary(aovp(p.RII.s~ts*tp,perm=999))
cld(emmeans(aovp(p.RII.s~ts*tp),~ts*tp),Letters=letters)


#

tr<-c(1,2,4,5,7,8,10,11,13,14,16,17)
t<-c(1,2,1,2,1,2,1,2,1,2,1,2)
mps<-c(mean(p.RII.QL),mean(p.RII.QT),mean(p.RII.TL),mean(p.RII.TT),
       mean(f.RII.QL),mean(f.RII.QT),mean(f.RII.TL),mean(f.RII.TT),
       mean(b.RII.QL),mean(b.RII.QT),mean(b.RII.TL),mean(b.RII.TT))
sdps<-c(sd(p.RII.QL)/sqrt(length(p.RII.QL)-1),sd(p.RII.QT)/sqrt(length(p.RII.QT)-1),
        sd(p.RII.TL)/sqrt(length(p.RII.TL)-1),sd(p.RII.TT)/sqrt(length(p.RII.TT)-1),
        sd(f.RII.QL)/sqrt(length(f.RII.QL)-1),sd(f.RII.QT)/sqrt(length(f.RII.QT)-1),
        sd(f.RII.TL)/sqrt(length(f.RII.TL)-1),sd(f.RII.TT)/sqrt(length(f.RII.TT)-1),
        sd(b.RII.QL)/sqrt(length(b.RII.QL)-1),sd(b.RII.QT)/sqrt(length(b.RII.QT)-1),
        sd(b.RII.TL)/sqrt(length(b.RII.TL)-1),sd(b.RII.TT)/sqrt(length(b.RII.TT)-1))
text<-c("b","a","a","a","b","","","","","","","")
text2<-c("","","","","","a","a","ab","","","","")
text3<-c("***","***","**","**","**","","","","","","","")
text4<-c("","","","","","***","(*)","","","","","")
plot.data.st<-data.frame(tr,mps,sdps,text,text2,text3,text4)

stoch<-ggplot()+ 
  geom_bar(data=plot.data.st,aes(x=tr,y=mps,fill=factor(t)),stat="identity",position="stack",colour="black",width = 0.75)+ 
  geom_errorbar(data=plot.data.st,aes(x=tr,ymin=mps-sdps,ymax=mps+sdps),position=position_dodge(width = 0.2),width=0.1)+
  geom_text(data=plot.data.st,aes(x=tr,y=0.1+(mps+sdps),label=text),color="black",fontface="bold",size=7.5)+
  geom_text(data=plot.data.st,aes(x=tr,y=(mps-sdps)-0.085,label=text2),color="black",fontface="bold",size=7.5)+
  geom_text(data=plot.data.st,aes(x=tr,y=0.025+(mps+sdps),label=text3),color="black",fontface="bold",size=7.5)+
  geom_text(data=plot.data.st,aes(x=tr,y=(mps-sdps)-0.05,label=text4),color="black",fontface="bold",size=7.5)+
  geom_hline(yintercept=0)+ 
  geom_hline(yintercept=0.75)+ 
  geom_vline(xintercept=0)+
  geom_vline(xintercept = 6,linetype=1)+
  geom_vline(xintercept = 12,linetype=1)+
  geom_vline(xintercept = 18,linetype=1)+
  geom_segment(aes(x = 3,xend=3,y=-0.6,yend=0.75),linetype=2)+
  geom_segment(aes(x = 9,xend=9,y=-0.6,yend=0.75),linetype=2)+
  geom_segment(aes(x = 15,xend=15,y=-0.6,yend=0.75),linetype=2)+
  scale_x_continuous(limits = c(0,18),
                     breaks=c(1.5,4.5,7.5,11.5,14.5,17.5),
                     labels=c("QL","TS","QL","TS","QL","TS"),
                     expand = c(0,0))+ 
  scale_y_continuous(limits = c(-0.6,1),
                     breaks=c(-0.5,-0.25,0,0.25,0.5,0.75),
                     expand = c(0,0))+
  scale_fill_manual(values=c("white","black"),
                    labels=c("loose","tight"))+
  labs(title=NULL, x=NULL, y="RIIstochasticity")+
  annotate("text",label = "Plant", x = 3, y = 0.875,size=5,color="gray30")+
  annotate("text",label = "Fungi", x = 9, y = 0.875,size=5,color="gray30")+
  annotate("text",label = "Bacteria", x = 15, y = 0.875,size=5,color="gray30")+
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
ggsave("stoch.pdf",stoch,width = 10,height = 6,dpi=300)


##relationship
bstoch<-data.frame(rbind(b.RII.QL,b.RII.QT,b.RII.TL,b.RII.TT))
fstoch<-data.frame(rbind(f.RII.QL,f.RII.QT,f.RII.TL,f.RII.TT))
pstoch<-data.frame(rbind(p.RII.QL,p.RII.QT,p.RII.TL,p.RII.TT))

bsd<-vegdist(scale(bstoch), "euclid")
fsd<-vegdist(scale(fstoch), "euclid")
psd<-vegdist(scale(pstoch), "euclid",na.rm = T)
mantel(bsd,psd,method = "spear")
mantel(fsd,psd,method = "spear")

# cld(emmeans(aovp(nst_bact_t$MST.ij.bray~g,perm=999),~g),Letters=letters)
# cld(emmeans(aov(nst_fungi_t$MST.ij.bray~group,data=nst_fungi_t),~group),Letters=letters)
# cld(emmeans(aov(nst_plant_t$MST.ij.bray~group,data=nst_plant_t),~group),Letters=letters)
