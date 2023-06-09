####fungi phy group
fungi1.phy<-funginame.phy[as.numeric(names(fungi1))]
fungi1.ab<-colSums(fungi1)
fungi1.phy<-data.frame(fungi1.phy,fungi1.ab)

fungi2.phy<-funginame.phy[as.numeric(names(fungi2))]
fungi2.ab<-colSums(fungi2)
fungi2.phy<-data.frame(fungi2.phy,fungi2.ab)

fungi3.phy<-funginame.phy[as.numeric(names(fungi3))]
fungi3.ab<-colSums(fungi3)
fungi3.phy<-data.frame(fungi3.phy,fungi3.ab)

fungi4.phy<-funginame.phy[as.numeric(names(fungi4))]
fungi4.ab<-colSums(fungi4)
fungi4.phy<-data.frame(fungi4.phy,fungi4.ab)

fungi5.phy<-funginame.phy[as.numeric(names(fungi5))]
fungi5.ab<-colSums(fungi5)
fungi5.phy<-data.frame(fungi5.phy,fungi5.ab)

fungi6.phy<-funginame.phy[as.numeric(names(fungi6))]
fungi6.ab<-colSums(fungi6)
fungi6.phy<-data.frame(fungi6.phy,fungi6.ab)
#
fungi1.phy.<-fungi1.phy[order(fungi1.phy$fungi1.phy),]
fungi1.phy.data<-data.frame()
for(i in 1:length(levels(fungi1.phy$fungi1.phy))){
  a<-fungi1.phy[which(fungi1.phy$fungi1.phy==levels(fungi1.phy$fungi1.phy)[i]),2]
  b<-sum(a)
  n<-levels(fungi1.phy$fungi1.phy)[i]
  fungi1.phy..<-data.frame(n,b)
  fungi1.phy.data<-rbind(fungi1.phy.data,fungi1.phy..)
}

fungi2.phy.<-fungi2.phy[order(fungi2.phy$fungi2.phy),]
fungi2.phy.data<-data.frame()
for(i in 1:length(levels(fungi2.phy$fungi2.phy))){
  a<-fungi2.phy[which(fungi2.phy$fungi2.phy==levels(fungi2.phy$fungi2.phy)[i]),2]
  b<-sum(a)
  n<-levels(fungi2.phy$fungi2.phy)[i]
  fungi2.phy..<-data.frame(n,b)
  fungi2.phy.data<-rbind(fungi2.phy.data,fungi2.phy..)
}

fungi3.phy.<-fungi3.phy[order(fungi3.phy$fungi3.phy),]
fungi3.phy.data<-data.frame()
for(i in 1:length(levels(fungi3.phy$fungi3.phy))){
  a<-fungi3.phy[which(fungi3.phy$fungi3.phy==levels(fungi3.phy$fungi3.phy)[i]),2]
  b<-sum(a)
  n<-levels(fungi3.phy$fungi3.phy)[i]
  fungi3.phy..<-data.frame(n,b)
  fungi3.phy.data<-rbind(fungi3.phy.data,fungi3.phy..)
}

fungi4.phy.<-fungi4.phy[order(fungi4.phy$fungi4.phy),]
fungi4.phy.data<-data.frame()
for(i in 1:length(levels(fungi4.phy$fungi4.phy))){
  a<-fungi4.phy[which(fungi4.phy$fungi4.phy==levels(fungi4.phy$fungi4.phy)[i]),2]
  b<-sum(a)
  n<-levels(fungi4.phy$fungi4.phy)[i]
  fungi4.phy..<-data.frame(n,b)
  fungi4.phy.data<-rbind(fungi4.phy.data,fungi4.phy..)
}

fungi5.phy.<-fungi5.phy[order(fungi5.phy$fungi5.phy),]
fungi5.phy.data<-data.frame()
for(i in 1:length(levels(fungi5.phy$fungi5.phy))){
  a<-fungi5.phy[which(fungi5.phy$fungi5.phy==levels(fungi5.phy$fungi5.phy)[i]),2]
  b<-sum(a)
  n<-levels(fungi5.phy$fungi5.phy)[i]
  fungi5.phy..<-data.frame(n,b)
  fungi5.phy.data<-rbind(fungi5.phy.data,fungi5.phy..)
}

fungi6.phy.<-fungi6.phy[order(fungi6.phy$fungi6.phy),]
fungi6.phy.data<-data.frame()
for(i in 1:length(levels(fungi6.phy$fungi6.phy))){
  a<-fungi6.phy[which(fungi6.phy$fungi6.phy==levels(fungi6.phy$fungi6.phy)[i]),2]
  b<-sum(a)
  n<-levels(fungi6.phy$fungi6.phy)[i]
  fungi6.phy..<-data.frame(n,b)
  fungi6.phy.data<-rbind(fungi6.phy.data,fungi6.phy..)
}

fungi.phy.data<-rbind(fungi1.phy.data,fungi2.phy.data,fungi3.phy.data,fungi4.phy.data,fungi5.phy.data,fungi6.phy.data)
group.phy<-c(rep(1,nrow(fungi1.phy.data)),rep(2,nrow(fungi2.phy.data)),rep(3,nrow(fungi3.phy.data)),
             rep(4,nrow(fungi4.phy.data)),rep(5,nrow(fungi5.phy.data)),rep(6,nrow(fungi6.phy.data)))
fungi.phy.data<-data.frame(group.phy,fungi.phy.data)


####################
names(fungi.phy.data)<-c("group.phy","name.phy","b" )

i<-1
for(i in 1:6){
  a<-which(fungi.phy.data$group.phy==i)
  b<-fungi.phy.data$b[a]
  fungi.phy.data$b[a]<-b/sum(fungi.phy.data$b[a])
}


#dismiss error value
fungi.phy.data$b[17]<-fungi.phy.data$b[17]+fungi.phy.data$b[18]
fungi.phy.data<-fungi.phy.data[-18,]
fungi.phy.data$name.phy[17]<-factor("Ascomycota;",levels = "Ascomycota;")

name.fungi<-c("Others **","Ascomycota **","Basidiomycota ***","Zygomycota ***")

fungi.phy<-ggplot() +
  geom_bar(data=fungi.phy.data,aes(x=group.phy,y=b,fill=factor(name.phy)),stat="identity",position="stack",width = 0.75)+ 
  geom_hline(yintercept=0)+ 
  geom_vline(xintercept=0)+ 
  scale_fill_manual(values= brewer.pal(9,"Set1")[c(1,5,3,2)],
                    labels = name.fungi)+
  scale_x_continuous(limits = c(0,7),
                     breaks=c(1,2,3,4,5,6),
                     labels=c("Group 1","Group 2","Group 3","Group 4","Group 5","Group 6"),
                     expand = c(0,0))+ 
  scale_y_continuous(expand = c(0,0))+
  labs(title=NULL, x=NULL, y="Fungi relative abundancce")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(color = "black"),
        legend.position="right",
        legend.key = element_blank(),
        legend.title= element_blank(),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=15,face = "bold"),
        axis.text.x = element_text(angle=45, hjust=1))

ggsave("W:\\\\论文发表\\\\fungi_phy.pdf",fungi.phy,width = 10,height = 6,dpi=300)





#bact.phy
names(bact.phy.data)<-c("group.phy","name.phy","b" )

i<-1
for(i in 1:6){
  a<-which(bact.phy.data$group.phy==i)
  b<-bact.phy.data$b[a]
  bact.phy.data$b[a]<-b/sum(bact.phy.data$b[a])
}

fix(bact.phy.data)
bact.phy.data$name.phy<-as.character(bact.phy.data$name.phy)
a<-nchar(bact.phy.data$name.phy)

bact.phy.data$name.phy[c(which(a<9),which(is.na(a)==TRUE))]<-c("Others")
bact.phy.data$name.phy<-as.factor(bact.phy.data$name.phy)

name.bact<-c("Others","Acidobacteria","Actinobacteria","Bacteroidetes","Candidate","Chlorobi",
  "Chloroflexi","Cyanobacteria","Firmicutes","Gemmatimonadetes","Planctomycetes",
  "Proteobacteria","Verrucomicrobia","Armatimonadetes","Deinococcaceae","Elusimicr",
  "Nitrospirae","Thermotogae")
length(name.bact)

display.brewer.all()

c(brewer.pal(9,"Reds")[7],)
color<-c(brewer.pal(9,"YlOrRd")[7],brewer.pal(9,"YlOrBr")[7],brewer.pal(9,"YlGnBu")[7],
         brewer.pal(9,"YlGn")[7],brewer.pal(9,"Reds")[7],brewer.pal(9,"RdPu")[7],
         brewer.pal(9,"Purples")[7],brewer.pal(9,"PuRd")[7],brewer.pal(9,"PuBuGn")[7],
         brewer.pal(9,"PuBu")[7],brewer.pal(9,"OrRd")[7],brewer.pal(9,"Oranges")[7],
         brewer.pal(9,"Greys")[7],brewer.pal(9,"Greens")[7],brewer.pal(9,"GnBu")[7],
         brewer.pal(9,"BuPu")[7],brewer.pal(9,"BuGn")[7],brewer.pal(9,"Blues")[7])

brewer.pal(12, "Accent")
bact.phy<-ggplot()+
  geom_bar(data=bact.phy.data,aes(x=group.phy,y=b,fill=factor(name.phy)),stat="identity",position="stack",width = 0.75)+ 
  geom_hline(yintercept=0)+ 
  geom_vline(xintercept=0)+ 
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(colourCount), 
                    labels = name.bact)+
  scale_x_continuous(limits = c(0,7),
                     breaks=c(1,2,3,4,5,6),
                     labels=c("Group 1","Group 2","Group 3","Group 4","Group 5","Group 6"),
                     expand = c(0,0))+ 
  scale_y_continuous(expand = c(0,0))+
  labs(title=NULL, x=NULL, y="Bacteria abundancce")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        legend.position="right",
        legend.key = element_blank(),
        legend.title= element_blank(),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=15,face = "bold"),
        axis.text.x = element_text(angle=45, hjust=1))

ggsave("W:\\\\论文发表\\\\bact_phy_Plus.pdf",bact.phy,width = 10,height = 6,dpi=300)

library(RColorBrewer)
library(ggplot2)
display.brewer.all()
RColorBrewer
colourCount<-18
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
 ggplot()+
  geom_bar(data=bact.phy.data,aes(x=group.phy,y=b,fill=factor(name.phy)),stat="identity",position="stack",colour="black",width = 0.75)+ 
  geom_hline(yintercept=0)+ 
  geom_vline(xintercept=0)+ 
  scale_fill_manual(values = getPalette(colourCount), 
                    labels = name.bact)
 
 
 color<-getPalette(18)
 ggplot()+
   geom_bar(data=bact.phy.data,aes(x=group.phy,y=b,color=color),stat="identity",position="stack",colour="black",width = 0.75)+ 
   geom_hline(yintercept=0)+ 
   geom_vline(xintercept=0)
   

###############################################
###plant 
gr<-rep(c(1:4),each=320)
grass<-c(grass1,grass2,grass3,grass4)
grass<-data.frame(gr,grass)
library(lmPerm)
grass$grass<-log(grass$grass+1)
hist(grass$grass)
grass$grass<-log(grass$grass+1)
summary(aov(grass~gr,data=grass))

legume<-c(legume1,legume2,legume3,legume4)
legume<-data.frame(gr,legume)
legume$legume<-log(legume$legume+1)
hist(legume$legume)
legume$legume<-log(legume$legume+1)
summary(aov(legume~gr,data=legume))

forb<-c(forb1,forb2,forb3,forb4)
forb<-data.frame(gr,forb)
forb$forb<-log(forb$forb+1)
hist(forb$forb)
forb$forb<-log(forb$forb+1)
summary(aov(forb~gr,data=forb))
 
#
Gramineae<-c("1","3","4","5")
Cyperaceae<-c("7","18")
Compositae<-c("8","11","15")
Cruciferae<-c("13")
Leguminosae<-c("23")
Rosaceae<-c("24","25")
Gentianaceae<-c("26")
pt<-list(Gramineae,Cyperaceae,Leguminosae,Compositae,Cruciferae,Rosaceae,Gentianaceae)

g<-c()
g.<-c()
for (i in 1:length(pt)) {
  if(length(intersect(names(plant1),pt[[i]]))!=0){
    a<-sum(plant1[,which(names(plant1)==intersect(names(plant1),pt[[i]]))])
    if(length(intersect(names(plant1),pt[[i]]))== 1){
      a.<-plant1[,which(names(plant1)==intersect(names(plant1),pt[[i]]))]
    }else{
      a.<-rowSums(plant1[,which(names(plant1)==intersect(names(plant1),pt[[i]]))])
    }
  }else{
    a<-0
    a.<-rep(0,320)
  }
  g<-c(g,a)
  g.<-c(g.,a.)
}
g1<-g
g1.<-g.

g<-c()
g.<-c()
for (i in 1:length(pt)) {
  if(length(intersect(names(plant2),pt[[i]]))!=0){
    a<-sum(plant2[,which(names(plant2)==intersect(names(plant2),pt[[i]]))])
    if(length(intersect(names(plant2),pt[[i]]))== 1){
      a.<-plant2[,which(names(plant2)==intersect(names(plant2),pt[[i]]))]
    }else{
      a.<-rowSums(plant2[,which(names(plant2)==intersect(names(plant2),pt[[i]]))])
    }
  }else{
    a<-0
    a.<-rep(0,320)
  }
  g<-c(g,a)
  g.<-c(g.,a.)
}
g2<-g
g2.<-g.

g<-c()
g.<-c()
for (i in 1:length(pt)) {
  if(length(intersect(names(plant3),pt[[i]]))!=0){
    a<-sum(plant3[,which(names(plant3)==intersect(names(plant3),pt[[i]]))])
    if(length(intersect(names(plant3),pt[[i]]))== 1){
      a.<-plant3[,which(names(plant3)==intersect(names(plant3),pt[[i]]))]
    }else{
      a.<-rowSums(plant3[,which(names(plant3)==intersect(names(plant3),pt[[i]]))])
    }
  }else{
    a<-0
    a.<-rep(0,320)
  }
  g<-c(g,a)
  g.<-c(g.,a.)
}
g3<-g
g3.<-g.

g<-c()
g.<-c()
for (i in 1:length(pt)) {
  if(length(intersect(names(plant4),pt[[i]]))!=0){
    a<-sum(plant4[,which(names(plant4)==intersect(names(plant4),pt[[i]]))])
    if(length(intersect(names(plant4),pt[[i]]))== 1){
      a.<-plant4[,which(names(plant4)==intersect(names(plant4),pt[[i]]))]
    }else{
      a.<-rowSums(plant4[,which(names(plant4)==intersect(names(plant4),pt[[i]]))])
    }
  }else{
    a<-0
    a.<-rep(0,320)
  }
  g<-c(g,a)
  g.<-c(g.,a.)
}
g4<-g
g4.<-g.

tr<-rep(c("A","B","C","D"),each=320)
for(i in 1:7){
  i1<-320*(i-1)+1
  i2<-i1+319
  pg<-c(g1.[i1:i2],g2.[i1:i2],g3.[i1:i2],g4.[i1:i2])
  show(summary(aovp(pg~tr,perm=999)))
}


l<-c(g1/sum(g1),g2/sum(g2),g3/sum(g3),g4/sum(g4))
taxa<-rep(c("Gramineae","Cyperaceae","Leguminosae","Compositae","Cruciferae","Rosaceae","Gentianaceae"),4)
group<-rep(c(1:4),each=7)
plant.phy.data<-data.frame(l,taxa,group)

##
library(RColorBrewer)

color<-c(brewer.pal(3,"OrRd")[c(3,2)],brewer.pal(3,"Greens")[2],brewer.pal(8,"Blues")[c(7,5,3,2)])

plant.phy<-ggplot()+
  geom_bar(data=plant.phy.data,aes(x=group,y=l,fill=taxa),stat="identity",position="stack",width = 0.75)+ 
  geom_hline(yintercept=0)+ 
  geom_vline(xintercept=0)+ 
  scale_fill_manual(values = color, 
                    labels = taxa)+
  scale_x_continuous(limits = c(0,5),
                     breaks=c(1,2,3,4),
                     labels=c("Group 1","Group 2","Group 3","Group 4"),
                     expand = c(0,0))+ 
  scale_y_continuous(expand = c(0,0))+
  labs(title=NULL, x=NULL, y="Plant taxa relative abundancce")+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        axis.line=element_blank(),
        legend.position="right",
        legend.key = element_blank(),
        legend.title= element_blank(),
        legend.text=element_text(size=14),
        axis.title=element_text(size=20,face="bold"),
        # axis.text.x = element_text(angle=45, hjust=1),
        axis.text=element_text(size=15,face = "bold"))

ggsave("E:/paper/planttaxa.pdf",plant.phy,width = 10,height = 6,dpi=300)

