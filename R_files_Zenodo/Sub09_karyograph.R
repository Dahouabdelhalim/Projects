###### Library ######
library(diversitree)
library(ggplot2)
library(ggpubr)
library(extrafont)

###### Font preparation ######
subset(fonttable(), FamilyName == "Arial") ##Check not more than one "Arial" in your system
loadfonts(device = "postscript")

###### Function for change between karyotype state vs. (arm no., chr no.) ######


#(chr number, arm number) -> state number
#this is (y,x) but not (x,y).
scal<-function(d,j){d*(d+1)/2+j-d}

#state number -> c(chr number, arm number)
chr_arm_cal<-function(s){
  d<- 1
  while(s-d*(d+1)/2 > d){
    d<- d+1
  }
  return(c(d,s-d*(d-1)/2))
}

#state vector -> matrix with two columns of chr number and arm number of each state
chr_arm_vec<-function(svec){
  chr_arm_mat<-c()  
  for(i in 1:length(svec)){
    chr_arm_mat<-rbind(chr_arm_mat,chr_arm_cal(svec[i]))
  }
  return(chr_arm_mat)
}

#state number -> c(arm number,chr number)
arm_chr_cal<-function(s){
  d<- 1
  while(s-d*(d+1)/2 > d){
    d<- d+1
  }
  return(c(s-d*(d-1)/2,d))
}

#state number -> karyotype "(x,y)"
karyotype_cal<-function(s){
  d<- 1
  while(s-d*(d+1)/2 > d){
    d<- d+1
  }
  return(paste("(",s-d*(d-1)/2,",",d,")",sep=""))
}

## Function, chr_arm_list
### vector of states -> list with [[1]] chr no. vector and [[2]] arm no. vector

chr_arm_list<-function(state){
  maxstate<-max(state)
  chr<-rep(0,length=length(state))
  arm<-rep(0,length=length(state))
  s<-1
  d<-1
  while(s <= maxstate){
    for(j in d:(2*d)){
      idents<-which(state==s)
      if(length(idents)>0){
        chr[idents]<-rep(d,length=length(idents))
        arm[idents]<-rep(j,length=length(idents))
      }
      s<-s+1
    }
    d<-d+1
  }
  return(list(chr,arm))
}



###### Fish data preparation for kartyograph ######

chrinfo<-read.table("Teleostei_karyotype_data_YK2021.csv",header=T,sep=",",na.strings="")
load("neot_state_YK2021.Robj")#817
load("otop_state_YK2021.Robj")#519

## Take teleostei except used one in PCM to do random sampling bellow
### This random sampling is for deletion of overlapped reports in one species
### Eurypterygii and Otophysi species not found in the phylogeny are also contained in this samples randomly sampled.

delindex<-which(is.na(match(chrinfo$Species_name,names(neot.state)))==F)
delindex<-c(delindex,which(is.na(match(chrinfo$Species_name,names(otop.state)))==F))
chrinfo2<-chrinfo[-delindex,]
val_species<-chrinfo2$Species_name
val_species<-unique(val_species) #1268 species

## Random sampling of one karyotype from one species
index_choice<-function(x,y){
  c_ind<-c()
  match_list<-charmatch(x,y)
  for(i in 1:length(x)){
    if(match_list[i]==0){
      c_ind[i]<-sample(which(y==x[i]),1)
    }else{
      c_ind[i]<-match_list[i]
    }
  }
  c_ind
}

val_index<-index_choice(val_species,chrinfo2$Species_name)

## Index of Used one in PCM
usedindex<-match(names(neot.state),chrinfo$Species_name)
usedindex<-c(usedindex,match(names(otop.state),chrinfo$Species_name))

## Whole valid karyotypes of teleostei

val_chrinfo<-rbind(chrinfo2[val_index,],chrinfo[usedindex,]) #37 orders 2604 species in total
val_chrinfo$Chr_num<-round(val_chrinfo$Chr_num/2)
val_chrinfo$Arm_num<-round(val_chrinfo$Arm_num/2)
val_chrinfo$State<-scal(val_chrinfo$Chr_num,val_chrinfo$Arm_num)
val_chrinfo$State[1269:2604]<-c(neot.state,otop.state) 

## Save val_chrinfo to fix the random sampling for another use.
save(val_index,file="val_chrinfo_XXXXXX.Robj")

teleostoutlist<-c("POLYPTERIFORMES","ACIPENSERIFORMES","LEPISOSTEIFORMES","AMIIFORMES")
outlist<-c("POLYPTERIFORMES","HIODONTIFORMES","ACIPENSERIFORMES","LEPISOSTEIFORMES","AMIIFORMES","OSTEOGLOSSIFORMES",
           "ELOPIFORMES","ANGUILLIFORMES","CLUPEIFORMES","GONORYNCHIFORMES","OSMERIFORMES",
           "SALMONIFORMES","ESOCIFORMES","ARGENTINIFORMES")
otoplist<-c("CYPRINIFORMES","CHARACIFORMES","SILURIFORMES","GYMNOTIFORMES")
neotlist<-c("AULOPIFORMES","MYCTOPHIFORMES","GADIFORMES","OPHIDIIFORMES",
            "MUGILIFORMES","ATHERINIFORMES","BELONIFORMES","CYPRINODONTIFORMES","STEPHANOBERYCIFORMES",
            "BERYCIFORMES","GASTEROSTEIFORMES","BATRACHOIDIFORMES","SYNBRANCHIFORMES","SCORPAENIFORMES",
            "PERCIFORMES","PLEURONECTIFORMES","LOPHIIFORMES","TETRAODONTIFORMES","ZEIFORMES","OPHIDIIFORMES")
sum(val_chrinfo$State>860)

#For whole teleostei
count_obj1<-table(val_chrinfo$State[is.na(match(val_chrinfo$Order,teleostoutlist))])

#For Eurypterygii
count_obj2<-table(val_chrinfo$State[is.na(match(val_chrinfo$Order,neotlist))==F])

#For Otophysi
count_obj3<-table(val_chrinfo$State[is.na(match(val_chrinfo$Order,otoplist))==F])

## (sub) different order
count_obj4<-table(val_chrinfo$State[val_chrinfo$Order=="CYPRINIFORMES"])
count_obj5<-table(val_chrinfo$State[val_chrinfo$Order=="CHARACIFORMES"])
count_obj6<-table(val_chrinfo$State[val_chrinfo$Order=="SILURIFORMES"])
count_obj7<-table(val_chrinfo$State[val_chrinfo$Order=="GYMNOTIFORMES"])
count_obj8<-table(val_chrinfo$State[val_chrinfo$Order=="BELONIFORMES"])
count_obj9<-table(val_chrinfo$State[val_chrinfo$Order=="GASTEROSTEIFORMES"])
count_obj10<-table(val_chrinfo$State[val_chrinfo$Order=="CYPRINODONTIFORMES"])
count_obj11<-table(val_chrinfo$State[val_chrinfo$Order=="CLUPEIFORMES"])


###### Drawing karyograph for whole Teleostei (S5 and S6 Figs) ######

pdata<-counting.karyotype(count_obj1,40)
pdata$pt<-ifelse(pdata$pz>=1,pdata$pz,NA)
pdata$ps<-ifelse(pdata$pz>=1,1,0)
pdata$pl<-log10(pdata$pz+1)
xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
labtitle<-expression(paste("log"[10],"(",italic("count")," + 1)"))

postscript("Fig_karyograph_plot_teleostei_log_XXXXXX.eps", height = 4, width = 7.4,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
ggplot(pdata)+geom_point(aes(x=px,y=py,size=pl),colour="black") +
  labs(x=xtitle,y=ytitle,size=labtitle) +
  annotate("text", x=6.5,y=39,size=6,label="Teleostei")+
  scale_size(range=c(0.02,2.6),limits=c(0,3)) +
  theme(text=element_text(size=13),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        legend.position=c(0.85,0.3))
dev.off()

postscript("Fig_karyograph_plot_teleostei_number_XXXXXX.eps", height = 4, width = 7.4,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
ggplot()+geom_point(aes(x=pdata$px,y=pdata$py,size=pdata$ps),colour="black")+
  scale_size(range=c(0.02,2)) +
  annotate("text", x=6.5,y=39,size=6,label="Teleostei")+
  geom_text(aes(x=pdata$px,y=pdata$py,label=pdata$pt),size=0.88,colour="white")+
  labs(x=xtitle,y=ytitle)+
  theme(text=element_text(size=13),legend.position="none",
        panel.border=element_rect(colour="black",size=0.6,fill=NA))

dev.off()


###### Karyograph for two groups with log-scale points (Fig2) ######

pdata<-counting.karyotype(count_obj2,40)
pdata$pl<-log10(pdata$pz+1)
xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
labtitle<-expression(paste("log"[10],"(",italic("count")," + 1)"))
g1<-ggplot(pdata)+geom_point(aes(x=px,y=py,size=pl),colour="black") +
  labs(x=xtitle,y=ytitle,size=labtitle) +
  annotate("text", x=10,y=39,size=5,label="Eurypterygii")+
  scale_size(range=c(0.02,2.4),limits=c(0,3)) +
  theme(text=element_text(size=10),legend.position="none",
        panel.border=element_rect(colour="black",size=0.6,fill=NA))

pdata<-counting.karyotype(count_obj3,40)
pdata$pl<-log10(pdata$pz+1)
g2<-ggplot(pdata)+geom_point(aes(x=px,y=py,size=pl),colour="black") +
  labs(x=xtitle,y=ytitle,size=labtitle) +
  annotate("text", x=7.5,y=39,size=5,label="Otophysi")+
  scale_size(range=c(0.02,2.4),limits=c(0,3)) +
  theme(text=element_text(size=10),legend.position=c(0.85,0.3),
        panel.border=element_rect(colour="black",size=0.6,fill=NA))

postscript("Fig_karyograph_plot_two_groups_XXXXXX.eps", height = 5.8, width = 5.2,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
ggarrange(g1,g2,ncol=1,nrow=2,labels=c("A","B"),widths=c(1,2),font.label=list(size=22,familty="Arial",face="plain"),hjust=-0.1,vjust=1.2)
dev.off()

###### Karyograph for two groups with the count number (S8 Fig) ######

pdata<-counting.karyotype(count_obj2,40)
pdata$pt<-ifelse(pdata$pz>=1,pdata$pz,NA)
pdata$ps<-ifelse(pdata$pz>=1,1,0)
xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
g1<-ggplot(pdata)+geom_point(aes(x=px,y=py,size=ps),colour="black")+
  scale_size(range=c(0.02,2)) +  
  geom_text(aes(x=px,y=py,label=pt),size=0.9,colour="white")+
  annotate("text", x=8.5,y=39,size=6,label="Eurypterygii")+
  labs(x=xtitle,y=ytitle)+
  theme(text=element_text(size=13),legend.position="none",
        panel.border=element_rect(colour="black",size=0.6,fill=NA))

pdata<-counting.karyotype(count_obj3,40)
pdata$pt<-ifelse(pdata$pz>=1,pdata$pz,NA)
pdata$ps<-ifelse(pdata$pz>=1,1,0)
g2<-ggplot(pdata)+geom_point(aes(x=px,y=py,size=ps),colour="black")+
  scale_size(range=c(0.02,2)) +  
  geom_text(aes(x=px,y=py,label=pt),size=0.9,colour="white")+
  annotate("text", x=6.5,y=39,size=6,label="Otophysi")+
  labs(x=xtitle,y=ytitle)+
  theme(text=element_text(size=13),legend.position="none",
        panel.border=element_rect(colour="black",size=0.6,fill=NA))

postscript("Fig_karyograph_plot_two_groups_number_XXXXXX.eps", height = 8.4, width = 7.4,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
ggarrange(g1,g2,ncol=1,nrow=2,labels=c("A","B"),widths=c(1,2),font.label=list(size=25,familty="Helvetica Neue",face="plain"),hjust=-0.1,vjust=1.2)
dev.off()

###### Histogram for the % acrocentric (S7 Fig) ######

max(val_chrinfo$Chr_num) #61
pdata<-counting.karyotype(count_obj1,61)
pdata$pa<-(2*pdata$py-pdata$px)/pdata$py
pvec1<-c()
for(i in 1:nrow(pdata)){
  pvec1<-c(pvec1,rep(pdata$pa[i],pdata$pz[i]))
}
length(pvec1)
plottitle<-expression(paste("Teleostei (",italic("N"),"= 2,587)"))
g1<-ggplot()+geom_histogram(aes(x=pvec1*100,y=(..count../sum(..count..))*100),binwidth=10)+
  coord_cartesian(ylim=c(0,50))+
 labs(x="% acrocentric chromosomes",y="% species",title=plottitle)+
  theme(text=element_text(size=13),plot.title=element_text(size=12),
        panel.border=element_rect(colour="black",size=0.6,fill=NA))

pdata<-counting.karyotype(count_obj2,61)
pdata$pa<-(2*pdata$py-pdata$px)/pdata$py
pvec2<-c()
for(i in 1:nrow(pdata)){
  pvec2<-c(pvec2,rep(pdata$pa[i],pdata$pz[i]))
}
length(pvec2)
plottitle<-expression(paste("Eurypterygii (",italic("N"),"= 1,368)"))
g2<-ggplot()+geom_histogram(aes(x=pvec2*100,y=(..count../sum(..count..))*100),binwidth=10)+
  coord_cartesian(ylim=c(0,50))+
  labs(x="% acrocentric chromosomes",y="% species",title=plottitle)+
  theme(text=element_text(size=13),plot.title=element_text(size=12),
        panel.border=element_rect(colour="black",size=0.6,fill=NA))

pdata<-counting.karyotype(count_obj3,61)
pdata$pa<-(2*pdata$py-pdata$px)/pdata$py
pvec3<-c()
for(i in 1:nrow(pdata)){
  pvec3<-c(pvec3,rep(pdata$pa[i],pdata$pz[i]))
}
length(pvec3)
plottitle<-expression(paste("Otophysi (",italic("N"),"= 1,030)"))
g3<-ggplot()+geom_histogram(aes(x=pvec3*100,y=(..count../sum(..count..))*100),binwidth=10)+
  coord_cartesian(ylim=c(0,50))+
  labs(x="% acrocentric chromosomes",y="% species",title=plottitle)+
  theme(text=element_text(size=13),plot.title=element_text(size=12),
        panel.border=element_rect(colour="black",size=0.6,fill=NA))

postscript("Fig_histogram_pa_XXXXXX.eps", height = 8.7, width = 4.35,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
ggarrange(g1,g2,g3, ncol=1,nrow=3,labels=c("A","B","C"),
          font.label=list(size=25,familty="Helvetica Neue",face="plain"),
          hjust=-0.1,vjust=1.2)
dev.off()
#ggsave(filename="FigS3_histogram_pa_for_XXXXXX.pdf",units="cm",width=10,height=20,device="pdf")

## Function, counting.karyotype
### make plot data, pdata, in different groups

counting.karyotype<-function(count_obj,chrmax){
  snum<-(chrmax+1)*(chrmax+2)/2-1
  pdata<-data.frame(px=rep(0,length=snum),py=rep(0),pz=rep(0))
  pvec<-c()
  avec<-c()
  mvec<-c()
  s<-1
  for(d in 1:chrmax){
    for(j in d:(2*d)){
      pdata$px[s]<-j
      pdata$py[s]<-d
      sind<-which(names(count_obj)==s)
      pdata$pz[s]<-ifelse(length(sind)>0,count_obj[sind],0)
      s<-s+1
    }
  }
  return(pdata)
}

###### Brassicacae karyograph (S15 Fig) ######

kdata<-read.table(file="Brassicaceae_karyotype_data_2021.txt",header=T,sep="\\t")
kdata$state<-scal(kdata$Chromosome_haploid,kdata$Arm_haploid)
ktable<-table(kdata$state)

chrmax<-25
snum<-(chrmax+1)*(chrmax+2)/2-1

scount<-rep(0,snum)
scount[as.numeric(names(ktable))]<-as.vector(ktable)

pdata<-data.frame(cbind(chr_arm_vec(1:snum),scount))
colnames(pdata)<-c("py","px","ps")
subset(fonttable(), FamilyName == "Arial") ##Check not more than one "Arial" in your system
loadfonts(device = "postscript")

postscript("Fig_bras_kryograph_XXXXXX.eps", height = 2.5, width = 5.2,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
ggplot(pdata)+
  geom_point(aes(x=px,y=py,size=ps))+
  scale_size_continuous(range=c(0,2),limits=c(0,4))+
  labs(x=xtitle,y=ytitle,size="Count")+
  theme(text=element_text(size=10),
        panel.border=element_rect(colour="black",size=0.6,fill=NA))
dev.off()

table(pdata$ps)
pdata$px[pdata$ps==4]
pdata$py[pdata$ps==4]
ktable #33 & 34