###### Parameters ######
#Determine chromosome number limit
chrmax<-35

###### Library ######
library(dplyr)
library(diversitree)
library(ggplot2)
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



###### Additional Functions ######
snum<-(chrmax+1)*(chrmax+2)/2-1

#This function colors areas with given cumultive probability.
#This coloring needs notice when one state has high probability.
#Different color just overlap and show only "50%".
coloring_prob<-function(ps){
  ordered_asr<-sort(ps,decreasing=T)
  range50<-sum(cumsum(ordered_asr)<0.50)+1
  range75<-sum(cumsum(ordered_asr)<0.75)+1
  range90<-sum(cumsum(ordered_asr)<0.90)+1
  pc<-rep("100%",length(ps))
  pc[ps>=ordered_asr[range90]]<-rep("90%")
  pc[ps>=ordered_asr[range75]]<-rep("75%")
  pc[ps>=ordered_asr[range50]]<-rep("50%")
  return(pc)
}

#Calculation of 95% range
range95cal<-function(probs){
  cumprob<-0
  range95<-c()
  chrindex<-1:length(probs)
  while(cumprob<0.95){
    index<-which(probs==max(probs))
    range95<-c(range95,chrindex[index])
    cumprob<-cumprob+probs[index]
    chrindex<-chrindex[-index]
    probs<-probs[-index]
  }
  lower<-min(range95)
  upper<-max(range95)
  return(c(lower,upper))
}

###### M0 Data preparation ######
## Eurypterygii (M0) pdata1

load("st_nodes_neot_M0_YK2021.Robj")

pdata1<-data.frame(chr_arm_vec(1:snum))
#px - arm number, py - chromosome number
colnames(pdata1)<-c("py","px")
#ps - size
pdata1$ps<-st.nodes.neot[,1]
#pc - color
pdata1$pc<-coloring_prob(st.nodes.neot[,1])
pdata1$pc<-factor(pdata1$pc,level=c("50%","75%","90%","100%"))
pdata1$Taxon<-rep("Eurypterygii (M0)")

## Otophysi (M0) pdata2

load("st_nodes_otop_M0_YK2021.Robj")

pdata2<-data.frame(chr_arm_vec(1:snum))
colnames(pdata2)<-c("py","px")
pdata2$ps<-st.nodes.otop[,1]
pdata2$pc<-coloring_prob(st.nodes.otop[,1])
pdata2$pc<-factor(pdata2$pc,level=c("50%","75%","90%","100%"))
pdata2$Taxon<-rep("Otophysi (M0)")


## Cyprinodontiformes (M0) pdata3

node<-1260-815

pdata3<-data.frame(chr_arm_vec(1:snum))
colnames(pdata3)<-c("py","px")
pdata3$ps<-st.nodes.neot[,node]
pdata3$pc<-coloring_prob(st.nodes.neot[,node])
pdata3$pc<-factor(pdata3$pc,level=c("50%","75%","90%","100%"))
pdata3$Taxon<-rep("Cyprinodontiformes (M0)")

pdata3 %>%
  group_by(py) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal
pdata3 %>%
  group_by(px) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal

## Goodeidae (M0) pdata4

node<-1277-815

pdata4<-data.frame(chr_arm_vec(1:snum))
colnames(pdata4)<-c("py","px")
pdata4$ps<-st.nodes.neot[,node]
pdata4$pc<-coloring_prob(st.nodes.neot[,node])
pdata4$pc<-factor(pdata4$pc,level=c("50%","75%","90%","100%"))
pdata4$Taxon<-rep("Goodeidae (M0)")

pdata4 %>%
  group_by(py) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal
pdata4 %>%
  group_by(px) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal


###### M0 Drawing MRCA of 3 taxons in Eurypterygii (S9 Fig) ######

pdata<-rbind(pdata1,pdata3,pdata4)
pdata$Taxon<-factor(pdata$Taxon,level=c("Eurypterygii (M0)","Cyprinodontiformes (M0)","Goodeidae (M0)"))

postscript("Fig_MRCA_neot_3taxon_M0_XXXXXX.eps", height = 6, width = 5,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
ggplot(pdata)+geom_point(aes(x=px,y=py,size=ps,col=pc)) +
  scale_color_manual(values=c("firebrick3","darkorange2","steelblue1","grey60"))+
  labs(x=xtitle,y=ytitle,size="Probability",color="Top") +
  scale_size(range=c(0,2),limits=c(0,1)) +
  facet_grid(Taxon~.,scales="free",space="free") +
  theme(text=element_text(size=10),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        strip.background = element_rect(color="black",size=0.5),
        legend.position="right",
        legend.spacing=unit(0,"in"),
        legend.key.size=unit(0.15,"in"),
        legend.box.spacing=unit(0,"in"),
        plot.margin=margin(5,5,5,5))
dev.off()

###### M0 Drawing MRCA only Otophysi (S10 Fig) ######

postscript("Fig_MRCA_otop_M0_XXXXXX.eps", height = 2, width = 5,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
ggplot(pdata2)+geom_point(aes(x=px,y=py,size=ps,col=pc)) +
  scale_color_manual(values=c("firebrick3","darkorange2","steelblue1","grey60"))+
  labs(x=xtitle,y=ytitle,size="Probability",color="Top") +
  scale_size(range=c(0,1.5),limits=c(0,0.1),breaks=c(0,0.02,0.04,0.06)) +
  facet_grid(Taxon~.,scales="free",space="free") +
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2))+
  theme(text=element_text(size=10),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        strip.background = element_rect(color="black",size=0.5),
        legend.margin=margin(3,0,0,10),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing=unit(0,"in"),
        legend.key.size=unit(0.15,"in"),
        legend.box.spacing=unit(0,"in"),
        plot.margin=margin(5,5,5,5))
dev.off()

###### M0 Calculation of 95% range of x or y axis ######

pdata1 %>%
  group_by(py) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal %>% print
pdata1 %>%
  group_by(px) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal %>% print
pdata2 %>%
  group_by(py) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal %>% print
pdata2 %>%
  group_by(px) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal %>% print




###### M2 Data preparation ######
## Eurypterygii M2 pdata1

load("st_nodes_neot_M2_YK2021.Robj")

pdata1<-data.frame(chr_arm_vec(1:snum))
colnames(pdata1)<-c("py","px")
pdata1$ps<-st.nodes.neot[,1]
pdata1$pc<-coloring_prob(st.nodes.neot[,1])
pdata1$pc<-factor(pdata1$pc,level=c("50%","75%","90%","100%"))
pdata1$Taxon<-rep("Eurypterygii (M2)")

## Otophysi M2 pdata2

load("st_nodes_otop_M2_YK2021.Robj")

pdata2<-data.frame(chr_arm_vec(1:snum))
colnames(pdata2)<-c("py","px")
pdata2$ps<-st.nodes.otop[,1]
pdata2$pc<-coloring_prob(st.nodes.otop[,1])
pdata2$pc<-factor(pdata2$pc,level=c("50%","75%","90%","100%"))
pdata2$Taxon<-rep("Otophysi (M2)")

## Cyprinodontiformes (M2) pdata3
node<-1260-815

pdata3<-data.frame(chr_arm_vec(1:snum))
colnames(pdata3)<-c("py","px")
pdata3$ps<-st.nodes.neot[,node]
pdata3$pc<-coloring_prob(st.nodes.neot[,node])
pdata3$pc<-factor(pdata3$pc,level=c("50%","75%","90%","100%"))
pdata3$Taxon<-rep("Cyprinodontiformes (M2)")

## Goodeidae (M2) pdata4
node<-1277-815

pdata4<-data.frame(chr_arm_vec(1:snum))
colnames(pdata4)<-c("py","px")
pdata4$ps<-st.nodes.neot[,node]
pdata4$pc<-coloring_prob(st.nodes.neot[,node])
pdata4$pc<-factor(pdata4$pc,level=c("50%","75%","90%","100%"))
pdata4$Taxon<-rep("Goodeidae (M2)")

###### M2 Drawing MRCA of 3 taxons in Eurypterygii (S12 Fig) ######
pdata<-rbind(pdata1,pdata3,pdata4)
pdata$Taxon<-factor(pdata$Taxon,level=c("Eurypterygii (M2)","Cyprinodontiformes (M2)","Goodeidae (M2)"))
postscript("Fig_MRCA_neot_3taxon_M2_XXXXXX.eps", height = 6, width = 5,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
ggplot(pdata)+geom_point(aes(x=px,y=py,size=ps,col=pc)) +
  scale_color_manual(values=c("firebrick3","darkorange2","steelblue1","grey60"),limits=c("50%","75%","90%","100%"))+
  labs(x=xtitle,y=ytitle,size="Probability",color="Top") +
  scale_size(range=c(0,2),limits=c(0,1)) +
  facet_grid(Taxon~.,scales="free",space="free") +
  theme(text=element_text(size=10),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        strip.background = element_rect(color="black",size=0.5),
        legend.position="right",
        legend.spacing=unit(0,"in"),
        legend.key.size=unit(0.15,"in"),
        legend.box.spacing=unit(0,"in"),
        plot.margin=margin(5,5,5,5))
dev.off()

###### M2 Drawing MRCA only otop (S13 Fig) ######

postscript("Fig_MRCA_otop_M2_XXXXXX.eps", height = 2, width = 5,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
xtitle<-expression(paste("Arm number, ",italic("x")))
ytitle<-expression(paste("Chromosome number, ",italic("y")))
ggplot(pdata2)+geom_point(aes(x=px,y=py,size=ps,col=pc)) +
  scale_color_manual(values=c("firebrick3","darkorange2","steelblue1","grey60"))+
  labs(x=xtitle,y=ytitle,size="Probability",color="Top") +
  scale_size(range=c(0,1.5),limits=c(0,0.1),breaks=c(0,0.02,0.04,0.06)) +
  facet_grid(Taxon~.,scales="free",space="free") +
  guides(size = guide_legend(order = 1),
         color = guide_legend(order = 2))+
  theme(text=element_text(size=10),
        panel.border=element_rect(colour="black",size=0.6,fill=NA),
        strip.background = element_rect(color="black",size=0.5),
        legend.margin=margin(3,0,0,10),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.spacing=unit(0,"in"),
        legend.key.size=unit(0.15,"in"),
        legend.box.spacing=unit(0,"in"),
        plot.margin=margin(5,5,5,5))
dev.off()


###### M2 Calculation of 95% range of x or y axis ######

pdata1 %>%
  group_by(py) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal
pdata1 %>%
  group_by(px) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal
pdata2 %>%
  group_by(py) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal
pdata2 %>%
  group_by(px) %>% summarise(probs=sum(ps)) %>% pull(probs) %>% range95cal

